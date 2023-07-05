---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
   display_name: Python 3
   language: python
   name: python3
---

(page:eos)=
# Equations of state

The following pressure-volume-energy-temperature equations of state (EOS) models are implemented in pyKO:
* [](eos:idg)
* [](eos:mgr)
* [](eos:ses)
* [](eos:til)

(eos:idg)=
## Ideal gas law (IDG)

The ideal gas is defined by the initial conditions, the ratio of specific heat capacities ($\gamma$), and the specific heat capacity at constant volume ($c_v$).  The initial conditions must specify the initial density, pressure, specific internal energy and particle velocity. The speed of sound is

$$
c_s = \left( \frac{P}{\rho} \right)^{1/2}.
$$

The ideal gas form for the sound speed is used in the Wilkins artificial viscosity formulation. The material EOS sound speeds are used to determine the time steps.

The specific heat capacity is assumed to be constant and used to determine temperature by

$$
T = \frac{E}{c_v}.
$$ 

The initial temperature is determined from $E_0$ and $c_v$.

Example configuration file entry (mks):
```
mat1:
    init:
        up0    : 0.0
        rho0   : 10.0E3
        p0     : 100.0E11
        e0     : 2.5E9
    eos: 
        name   : 'IDG 1'
        type   : 'IDG'
        gamma  : 1.4
        cv     : 2.5E8
```

(eos:mgr)=
## Mie-Gr&uuml;neisen EOS (MGR)

The Mie-Gr&uuml;neisen model uses a reference curve, usually the principal Hugoniot, and a thermal parameter to describe a pressure-volume-energy EOS surface (e.g., see Chapter 5 in {cite}`Forbes2012` and <a href="https://impactswiki.net/impact-tools-book/im/manual.html">notes for impedance matching</a>). 

The MGR material model requires a reference density ($\rho_0$), linear $U_S-u_p$ Hugoniot parameters ($c,s$), the reference state Gr&uuml;neisen gamma ($\gamma_0$), initial temperature ($T_0$), and specific heat capacity ($c_v$). The Wilkins implementation of the MGR EOS only uses a linear Hugoniot and constant $\gamma/V$. $e_0$ is calculated from the initial temperature and specific heat capacity. The pressure and energy equations are 4th order polynomials in strain and are only valid for a small range of shock pressures (typically 10's GPa).

Required initial conditions: density, pressure, energy, temperature and particle velocity. To do: allow either energy or temperature for input.

Example configuration file entry (mks). Note input keys $c0 = c$ and $s1=s$ in the equations below.
```
mat1:
    init:
        up0    : 0.0
        rho0   : 8930.0
        p0     : 0.0
        e0     : 0.0
    eos:
        name   : 'Cu test'
        type   : 'MGR'
        rhoref : 8930.0
        c0     : 3900.0
        s1     : 1.49
        gamma0 : 1.99
        cv     : 1.0E6
```

The thermal Gr&uuml;neisen parameter $\gamma$ is given by

$$ 
\gamma(V) = V \left[\frac{dP}{dE}\right]_V, 
$$

where $\gamma$ is often assumed to only depend on volume. Note that parentheses are used to denote dependent variables.

In pyKO, it is assumed that $\gamma$ is solely dependent on V and

$$
\frac{\gamma(V)}{V} = \frac{\gamma_0}{V_0} = {\rm constant.}
$$

A word of caution: these simple equations for the shock Hugoniot and the Gr&uuml;neisen parameter are unlikely to be thermodynamically self-consistent over wide ranges of volume. Beware using the Mie-Gr&uuml;neisen model over wide ranges of pressure and density. 

The Mie-Gr&uuml;neisen equation of state defines a thermal pressure term, $P(V,E)$, at specific volume $V$ and specific internal energy $E$ referenced to a known state at $V$, which in this case is the principal shock Hugoniot with pressure $P_H(V)$ and specific internal energy $E_H(V)$. Then the thermal pressure term is given by 

$$ 
P(V,E) = P_H(V)+\frac{\gamma(V)}{V} \left[E(V)-E_H(V)\right]. 
$$

A single shock must satisfy the Rankine-Hugoniot conservation equations between an initial state, $s_0$, and the shock state, $s_H$. The initial state $s_0$ is defined by the state variables: specific volume $V_0$, presure $P_0$, and specific energy $E_0$ (which then specifies the initial specific entropy $S_0$ and temperature $T_0$). Note the tabulated initial value for the Gr&uuml;neisen parameter $\gamma_0$ must correspond to state $s_0$. The energy conservation equation is 

$$ 
E_H(V)=E_0+\left[\frac{P_0+P_H(V)}{2}\right]\left[V_0-V\right]. 
$$

Combining these two equations defines the Mie-Gr&uuml;neisen equation of state surface, which is limited to the volume range where $P_H(V)$ is defined,

$$
P(V,E) & =  P_H(V)+\frac{\gamma(V)}{V}\left[E(V)-E_0\right] - \left[P_0+P_H(V)\right]\frac{\gamma(V)}{2V}[V_0-V]  \\
       & =  \frac{\gamma(V)}{V}\left[E(V)-E_0\right] + P_H(V)\left\{ 1- \frac{\gamma(V)}{2V}\left[V_0-V\right] \right\}
        - \frac{P_0\gamma(V)}{2V}\left[V_0-V\right],
$$

and the last term is dropped when $P_0=0$. 

The MGR EOS model sound speed was not calculated in {cite}`Wilkins1999` (he used the ideal gas equation for sound speed with a minimum of the reference sound speed). The adiabatic sound speed $c_s$ is calculated here from the adiabatic bulk modulus $K_s$:

$$
c_s^2 = \frac{K_s}{\rho} = - V^2 \frac{dP}{dV},
$$

where specific volume $V=1/\rho$. $P$ is a function of specific volume $V$ and specific internal energy $E$; expand the total derivative

$$
c_s^2 = - V^2 \left[ \left. \frac{\partial P}{\partial V} \right|_E \frac{\partial V}{\partial V} + \left. \frac{\partial P}{\partial E} \right|_V \frac{\partial E}{\partial V} \right].
$$

Along an isentrope (adiabat), $dS=0$. From the reversible form of the first law of thermodynamics,

$$
dE = TdS - PdV = -P dV, \\
\left. \frac{\partial E}{\partial V} \right|_S = -P.
$$

Substitute density for specific volume,

$$
\partial V = \partial \left( \frac{1}{\rho} \right) = - \frac{1}{\rho^2} \partial \rho.
$$

Then,

$$
c_s^2 = - \frac{1}{\rho^2} \left[ - \rho^2 \left. \frac{\partial P}{\partial \rho} \right|_E - P \left. \frac{\partial P}{\partial E} \right|_V \right] \\
c_s^2 = \left. \frac{\partial P}{\partial \rho} \right|_E + \frac{P}{\rho^2} \left. \frac{\partial P}{\partial E} \right|_V
$$

Wilkins uses the following expansion for $P(V,E)$:

$$
P(V,E) = P_0 + \frac{\gamma_0}{V_0} [E - E_0]\\
P(V,E) = \rho_0 c^2 \left\{ x + \left[ 2s - \frac{\gamma_0}{2} \right] x^2 + s [3s - \gamma_0] x^3 \right\} + \frac{\gamma_0}{V_0} E \\
P(V,E) = k_1 x + k_2 x^2 + k_3 x^3 + \frac{\gamma_0}{V_0} E \\
$$

where $U_s = c + s u_p$, $E$ is specific internal energy, $x=1-V/V_0$ is strain.

Temperature is determined from 

$$
T(V,E) = \frac{1}{c_v} [ E(V)-E_{T_0}(V) ]
$$

where $E_{T_0}(V)$ is the specific internal energy along the initial temperature $T_0$ isotherm:

$$
E_{T_0}(V) = \epsilon_{00} + \epsilon_{01} x + \epsilon_{02} x^2 + \epsilon_{03} x^3 + \epsilon_{04} x^4,\\
\epsilon_{00} = -3R T_0 = E_0,\\
\epsilon_{01} = \gamma_0 \epsilon_{00}\\
\epsilon_{02} = [1/2][c^2 + \gamma_0^2 \epsilon_{00}]\\
\epsilon_{03} = [1/6][4sc^2 + \gamma_0^3 \epsilon_{00}]\\
\epsilon_{04} = [1/24][18s^2c^2 - 2\gamma_0 s c^2 + \gamma_0^4 \epsilon_{00}]
$$

Finally, the pressure along the initial temperature isotherm is

$$
P_{T_0}(V) = - \left. \frac{\partial E}{\partial V} \right|_{T=T_0}
$$

(eos:ses)=
## Tabular EOS (SES)

A tabular EOS contains the thermodynamic variables for the material over a rectangular grid of two independent variables, typically density and temperature. The pressure-volume-energy relationships are interpolated during problem execution to solve the energy conservation equation. pyKO uses bilinear interpolation on the density-temperature grid, but other tabular EOS schemes are easy to implement.

The SESClass.readses function assumes STD and EXT formats compatible with the eos_table.py module from the Stewart group.

* SESAME standard 201+301 format. This table format specficiation is described in {cite:t}`Lyon1992` (<a href="https://github.com/ststewart/aneos-forsterite-2019/blob/master/EOS-docs/Lyon-Johnson-1992-SESAME-database.pdf">PDF file</a>).

The units in the tabular EOS must be specified in the pyKO configuration file (see [](config:units)).

Stewart Group ANEOS tables (https://github.com/ststewart?tab=repositories) have the following format:

* NEW-SESAME-STD.TXT: Standard length SESAME file with 201 table and 301 table (density, temperature, pressure, sp. internal energy, Helmholtz free energy). 301 table units: g/cm<sup>3</sup>, K, GPa, MJ/kg, MJ/kg.

* NEW-SESAME-EXT.TXT: SESAME-style table with extra variables from ANEOS. Contains the standard 201 table and non-standard 301-extra-variables EOS table. The 301 table has: density grid values, temperature grid values, sp. entropy(T,rho), sound speed(T,rho), sp. heat capacity(T,rho), KPA flag(T,rho). 2-D arrays list all densities, looping over each temperature. 301 table units: g/cm<sup>3</sup>, K, MJ/K/kg, cm/s, MJ/K/kg, integer flag. The KPA flag is an ANEOS output with phase information.

* See ANEOS developments and documentation in this directory: https://github.com/ststewart/aneos-forsterite-2019/tree/master/EOS-docs

To use a tabular EOS, specify the path, std and ext file names and the units for the file (in the units section). To use a different tabular format, customize the readses function to interface with the eos_tables.py module.

Initialization of materials with a tabular EOS requires some care by the user. There are some functions provided to aid in the initial density and internal energy settings. See the tutorial notebooks for examples. Initialization near zero pressure is tricky.

EOS table interpolation is bilinear with simple search criteria for finding positions in the table. Future development is needed for EOS tables with discontinuities/degeneracies (e.g., with multiple solid phases).
```
mat1:
    init:
        up0    : 0.0
        rho0   : 3220.0
        p0     : 1.e5
        e0     : 73277.974035
    eos: 
        name   : 'Forsterite'
        type   : 'SES'
        path    : '../aneos-forsterite-2019/'
        filestd : 'NEW-SESAME-STD.TXT'
        fileext : 'NEW-SESAME-EXT.TXT'
```


(eos:til)=
## Tillotson EOS (TIL)

The equation of state model by {cite}`Tillotson1962` can span a larger pressure-volume space compared to the MGR model. The EOS is divided into compressed and expanded forms with an interpolated region in between. See {cite}`Brundage2013` and {cite}`Melosh1987` (Appendix A) for more details about the Tillotson EOS.

The Tillotson EOS requires 11 parameters: $\rho_0$, $E_0$, $E_{IV}$, $E_{CV}$, $A$, $B$, $a$, $b$, $\alpha$, $\beta$, $c_v$
The mks units are kg/m$^3$, J/kg, J/kg, J/kg, Pa, Pa, and $a$, $b$, $\alpha$, and $\beta$ are dimensionless. $c_v$ is J/K/kg.

There are 4 regions in the Tillotson EOS: 1-compressed, 2-interpolated, 3-expanded, 4-low energy expansion.

For compressed states (region 1), where $\rho \ge \rho_0$ and $0<E<E_{IV}$:

$$
P_1(\rho,E) = \left[ a + \frac{b}{ \frac{E}{E_0 \eta^2} + 1 } \right] \rho E + A \mu + B \mu^2
$$

where $\eta=\rho/\rho_0$ and $\mu=\eta - 1$.

For expanded states (region 3), where $\rho \le \rho_0$ and $E>E_{CV}$:

$$
P_3(\rho,E) = a \rho E + \left\{ \frac{ b \rho E }{ \frac{E}{E_0 \eta^2} + 1 } + A \mu e^{- \beta [ [ \rho_0 / \rho]-1 ]} \right\}  e^{- \alpha [[ \rho_0 / \rho ] - 1 ]^2}
$$

For intermediate states, these equations are interpolated (region 2). When $\rho < \rho_0$ and $E_{IV}<E<E_{CV}$:

$$
P_2(\rho,E) = \frac{ [E-E_{IV}] P_3 + [E_{CV}-E] P_1 } {E_{CV}-E_{IV}}
$$

{cite}`Brundage2013` identifies a cold expanded state (region 4) where $\rho_0 < \rho < \rho_{IV}$ and $E \le E_{IV}$. In this case, the equation for $P_1$ is used with $B=0$. Here we have not calculated the $\rho_{IV}$ and the fourth region is entered when $\rho > \rho_0$ and $E \le E_{IV}$.

The initial temperature isotherm is given by

$$
E(\rho,T_0) = E_c(\rho) + E_{T0}
$$

where $E_{T0}=c_v T_0$ and the cold curve is calculated with a fourth-order Runge-Jutta integration of

$$
\frac{dE_c}{d\rho} = \frac{P_1(\rho,E_c)}{\rho^2}
$$

with initial conditions: $\rho=\rho_0$ and $E_c=0$. The temperature is determined relative to this cold curve in the same manner as for the MGR EOS:

$$
T(\rho) = T_0 + \frac{ E(\rho) - E_c(\rho) } {c_v}
$$

The Grueneisen gamma is $\gamma_0 = a+b$, where $a$ is the high compression limit (usually 0.5; or 2/3 for a Fermi electron gas). The linear shock velocity coefficients are related to the Tillotson parameters by

$$
c = [A/\rho_0]^{1/2}\\
s = \{ 1 + B/A + [a+b]/2 \} / 2
$$

Thus $A$ is the bulk modulus at zero pressure. The sound speed is calculated by taking the derivative of pressure with respect to volume as described in the MGR EOS section above.

$E_{IV}$ and $E_{CV}$ are the incipient and complete vaporization energy. $E_0$ is not the initial specific internal energy; it is a value often close to $E_{IV}$.


```
# parameters: [rho0, E0, EIV, ECV, AA, BB, a, b, alpha, beta]
# units:    [kg/m3, J/kg, J/kg, J/kg, Pa, Pa, [-]x4]
# these olivine parameters are from Marinova et al. 2011 Icarus 
# dunitetill = [3500.0, 550.0e+6, 4.500e+6, 14.50e+6, 131.00e+9,  49.00e+9, 0.5, 1.4, 5.0, 5.0]
# Basalt parameters from iSALE -- from where? Benz?
# basalttill = [2650.0, 4.87E8, 4.72E6, 18.2E6, 5.3E10, 5.3E10, 0.6, 0.6, 5., 5.]
```

Example configuration file entry (mks) from {cite}`Marinova2011`:
```
mat1:
    init:
        up0    : 0.0
        rho0   : 3500.0
        p0     : 0.0
        e0     : 0.0
        t0     : 298.0
    eos:
        name   : 'Olivine'
        type   : 'TIL'
        rhoref : 3500.0
        E0     : 550.0E6
        EIV    : 4.5E6
        ECV    : 14.5E6
        AA     : 131.0E9
        BB     : 49.0E9
        a      : 0.5
        b      : 1.4
        alpha  : 5.0
        beta   : 5.0
        cv     : 830.0
```

