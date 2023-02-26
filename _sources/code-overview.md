# Code overview

## Time and spatial domains

The KO code uses a 2nd order time centered finite difference method. The central calculation involves 4 half time steps and nj spatial nodes. The spatial domain is comprised of mesh edges (even j nodes), where position and velocity are defined, and mesh centers (odd j nodes), where thermodynamic parameters are defined. The main hydro arrays have size (nt,nj) with nt(=4) time indices and nj(=number of nodes) spatial indices.

In Wilkins notation, delta_t<sup>n</sup> = t<sup>n+1/2</sup> - t<sup>n-1/2</sup> is used with the pressure field to advance the Lagrangian nodes from time t<sup>n-1/2</sup> to t<sup>n+1/2</sup>.
The position of mesh edges are advanced delta_t<sup>n+1/2</sup> = t<sup>n+1</sup>-t<sup>n</sup> == one true time step, which is determined from the stability conditions with delta_t<sup>n</sup> = (1/2)(delta_t<sup>n+1/2</sup> + delta_t<sup>n-1/2</sup>). Time centering is upset when stability conditions permit deviation from these relationships.

Mapping the Wilkins notation to the python time array indices:
* t(0) = n-1/2
* t(1) = n
* t(2) = n+1/2
* t(3) = n+1
* True time step = t(3)-t(1)

## Principal classes

* <b>RunClass</b>: This class holds all of the simulation parameters. The configuration file parameters are read into this class. 
    * <i>checkinput</i> function: prints the simulation parameters for user inspection
* <b>DomainClass</b>: This class holds the primary (time,length) arrays for the hydrocode calculation.
    * <i>makegrid</i> function: initializes the hydro arrays in the DomainClass object based on RunClass information, filling the t(2) and t(3) positions
    * <i>advancetime</i> function: advances the hydro arrays one time step
    * <i>shiftdata</i> function: at the beginning of a time step, the time indices t(2) and t(3) are shifted to positions t(0) and t(1)
    * <i>createinteriorboundary</i> function: inserts a pair of free surfaces and void cell in the spatial domain
    * <i>closeinteriorboundary</i> function: removes a pair of free surfaces and void cell in the spatial domain
    * <i>checkforcontact</i> function: if interior free surface(s) is(are) present, check for contact during the current time step
    * <i>writefirstoutput</i> function: creates ascii or binary output file and writes the problem configuration at time zero
    * <i>appendoutput</i> function: appends ascii snapshot data to output file
    * <i>binaryoutputpint</i> function: appends binary snapshot data (OutputClass object) to output file using pickle and pint
    
## Material properties and boundary conditions

Equations of state
* <b>IDGClass</b>: ideal gas material parameters
* <b>MGRClass</b>: Mie-Grueneisen material parameters and functions
* <b>SESClass</b>: Class to hold tabular EOS and associated functions; requires the eos_table.py module.
    * <i>readtable</i> function: Reads in Stewart group style NEW-SESAME-STD.TXT and NEW-SESAME-EXT.TXT tables.
    * <i>sesp</i>: Main run-time function for bilinear interpolation on spatial domain vector to extract table P and T from rho and U.
    * <i>onept(rho,u)</i>: Bilinear interpolation to extract table one point P and T from rho and U. Used to help define initial conditions.
    * <i>onepu(rho,t)</i>: Bilinear interpolation to extract table one point P and U from rho and T. Used to help define initial conditions.

Strength
* <b>HydroClass</b>: Class to indicate hydrodynamic material
* <b>VonMisesClass</b>: Class to hold strength parameters
* <b>FractureClass</b>: Class to hold fracture parameters

Boundaries
* <b>BCClass</b>: Class to hold boundary condition parameters


## Input and output functions and classes
* <b>readinput_borg</b> function: reads formatted ascii problem initialization file that is compatible with modified fortran KO. Populates the RunClass object.
* <b>readinput_yaml</b> function: reads yaml formatted configuration file into a dictionary object (named config) and then populates the main RunClass with simulation parameters (object named run). The input parameter units are converted to code units using pint. Tabular EOS are converted to code units. 
* <b>run</b> function: primary wrapper function for I/O and calling main hydro time steps
* <b>OutputClass</b>: This class holds the data arrays for binary snapshot output using pickle. Create a new OutputClass object to load data from the binary pickle file. The OutputClass data arrays contains only mesh-centered values (odd j), with position and velocity interpolated from the mesh edges (even j). Interior void spaces are included as an empty mesh center.

### Customizing the binary snapshots
The default values in OutputClass can be customized by editing the arrays in the class and DomainClass.binaryoutputpint function. Several DomainClass variables are not included in the default binary snapshot. For pandas compatibility, all arrays in the OutputClass must be the same length (all oddj).

```
class OutputClass:
    """ Problem domain data structure for output dumps """
    def __init__(self):
        """ Initialize the main arrays for the problem domain """
        # unknown number of cells initially. Initialize data types with zero length.
        # variables in space at a particular snapshot in time
        self.stepn    = np.zeros(0,dtype='int') # counter for the number of time steps; int
        self.time     = np.zeros(0) # same time for each cell; used to make pandas conversion easier
        self.mat      = np.zeros(0,dtype='int') # material id number
        self.pos      = np.zeros(0) # position in 1D geometry (x,r_cyl,r_sph)
        self.rho0     = np.zeros(0) # initial density rho0 g/cm3
        self.rho      = np.zeros(0) # local density g/cm3
        self.up       = np.zeros(0) # particle velocity
        #self.iev0     = np.zeros(0) # internal energy per initial volume 1e12 erg/cm3
        self.ie       = np.zeros(0) # specific internal energy 
        self.pres     = np.zeros(0) # pressure Mbar
        self.mass     = np.zeros(0) # mass in node g
        self.temp     = np.zeros(0) # temperature K
        self.sigmar   = np.zeros(0) # total radial stress sigma_r, only evaluated at n=1 current time
        self.sigmao   = np.zeros(0) # total tangential stress sigma_theta, only evaluated at n=1 current time
        self.etot     = np.zeros(0) # total IE + total KE
        self.j        = np.zeros(0) # node index for comparisons to fKO
        #self.vr       = np.zeros(0) # relative volume to initial state=rho0/rho=1/eta in Wilkins [dimless]
        #self.phi      = np.zeros(0) # rho0*(pos_j+1-pos_j)/vr_j+1/2 = (pos_j+1-pos_j)/rho_j+1/2=phi in Wilkins
        #self.q        = np.zeros(0) # artificial viscosity Mbar
        #self.eps1     = np.zeros(0) # velocity strain dup/dpos, per microsec
        #self.eps2     = np.zeros(0) # velocity strain up/pos, per microsec
        # these variables are evaluated at n=3; full time step ahead
        #self.beta     = np.zeros(0) # (sigmar-sigmao)/(0.5 dpos)*rho = beta in Wilkins
        #self.s1       = np.zeros(0) # s1 stress deviator deriv = 2G(deps1-dvr/vr/3)
        #self.s2       = np.zeros(0) # s2 stress deviator deriv
        #self.s3       = np.zeros(0) # s3 stress deviator deriv; not needed in 1D; kept for equation clarity
        #self.entropy  = np.zeros(0) # entropy [DEBUG check eu/cm3 units]
        # variables in time only - n
        #self.ibc      = np.zeros(0,dtype='int') # boundary condition id number, integer
        #self.yld      = np.zeros(0) # yield stress
        #self.pfrac    = np.zeros(0) # fracture stress
        #self.deltaz   = np.zeros(0) # distortion energy Mbar
        #self.alocal   = np.zeros(0) # intermediate variable for AV
```

## Units

pyKO uses the <a href="https://pint.readthedocs.io/en/stable/">pint package</a> to convert user-specified input units to a self-consistent set of internal code units. The data in the binary output files are converted from code units to the same user-specified input units. The example input files use mks units.

The pint package is also used to convert any tabular EOS files into code units. Refer to the configuration file section on [](config:units).

The original KO code units are suitable for calculations of lab-scale gas gun impact experiments: microseconds, cm, g, and 1 eu = 10<sup>12</sup> ergs. The code normalizes volume to the initial volume (v0) so that changes in density are internalized as relative volume. To maintain self-consistent units, all specific quantities are multiplied by the initial density, so that specific internal energy (eu/g) becomes internal energy per initial volume (eu/cm<sup>3</sup>).

````
    time            : 'microseconds'
    length          : 'cm'
    mass            : 'g'
    density         : 'g/cm^3'
    relative_volume : 'dimensionless' 
    velocity        : 'cm/microsecond'
    pressure        : 'megabar'
    temperature     : 'K'
    energy          : 'eu'
    sp_energy       : 'eu/g'
    sp_heat_cap     : 'eu/K/g'
    sp_entropy      : 'eu/K/g'
    ie_perv0        : 'eu/cm^3'
    cv_perv0        : 'eu/cm^3/K'
````
