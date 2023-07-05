;=========================================================================
; make a CTH SESAME EOS table for H2O
; liquid, vapor, ice Ih, ice VI, ice VII
;
; created STSM 11/20/2006
; v2 done 11/29/2006
; v3 12/6/2006 -- add tension region: extrapolate each density to low
;                         temp using 3rd order Birch Murnaghan eq.,
; v3.1 12/7/2006 -- fixed zero density column to extrapolate energy, entropy
;                   from adjacent values; extrapolate energy and entropy in
;                   tension region (this is a fudge); add more points to ice
;                   Ih region
; v4.0 12/8/2006 -- make dP/dT a small positive number for Ih-VI and VI-VII
;                   phase boundaries
; v4.3 12/20/2006 -- latest meshing and tension region
; v4.3_clean  3/4/2007  -- clean version for sharing
; v4.3  5/9/2007 -- bug around 2.8 g/cc needs to be fixed
;                -- kink in ice VII melting curve - from NIST to Frank et al.?
;                -- some overlapping phase space with the non-vertical P-T ice
;                   VI and VII phase boundaries
;                -- melting curve e and s are bogus
;
; v5.0  5/9/2007 -- working on it
;                -- ice VI-VII boundary fixed dP/dT
;                -- 2 cases of '7' should be ' 7'
;                -- NIST TMELTP program was giving 0's for E and S! - fixed now
;                -- NIST EOSTP program was giving 0's for E and S! - fixed now
;                -- NIST PSATT program was giving 0's for E and S! - fixed now
;                -- NIST EOSTD check added
;                -- check sublimation curve - vapor side ; little gap
;                -- 'ES' extrapolation now >3000K - FIXED extrapolation
;                   code - fixed kink at 2.7 g/cc
;                -- 'L6' region below 273K fixed interpolation problems
;                -- ' 1' entropy overlaps with 'L6' in T-S space -
;                   right? yes
;                -- ' 7' low temp region overlapping temp in T-S space
;                   - real? no -- hard to calc ds with integral of cp
;                -- checked interpolated ' S' region (line 1019) - much better
;                   now
;
; v6.0 4/21/2008 -- correcting some errors with v5 table
;                   -- ice Ih-vapor equilibrium field has mass weighted
;                      averaging errors for densities <0.5(rho_ice) for
;                      entropy and energy (found by Robert Marcus) - FIXED
;                   -- need to add tension region for liquid water
;                   -- try new definition for some phase boundaries:
;                      ice Ih to ice VI: remains ice Ih at increasing
;                      densities until abrupt transformation to ice VI
;                      ice VI to ice VII: remains ice VI at increasing
;                      densities until abrupt transformation to ice
;                      VII -- THIS IS BAD; UNREALISTIC HIGH PRESSURES
;                -- changed ice sublimation curve low T extrapolation
;                   in FIRST setup run
; v6.1 6/11/2008 -- keep normal phase boundaries as mixed phases
;                   weighted by boundary densities of each phase
;                -- option A for ice Ih in tension
; v6.2 6/13/2008 -- add liquid tension region (reduces field of LV
;                   mixture) for densities of 600 to the L phase
;                   boundary
;                -- this uses OPTION A for extrapolating E and S --
;                   this is probably bogus, but will run stably for
;                   now -- FIX LATER
;                -- 'LT' and '1T' flags used for tension regions
; v6.3 6/13/2008 -- option added for flag for no tension regions - for
;                   robert marcus
; v6.4 6/15/2008 -- fixing some low temp errors:
;                   * need to fix low T ice VII entropy - fudged to be
;                                                         the same as
;                                                         ice VI
;                -- added flag to phasearr for zero density and zero
;                   temperature columns 'ZZ' for zero boundary
;                -- added temps to zero density column - was incorrect
;                   in tijarr - fine for the actual table
;                -- low T extrapolation of the ice sublimation curve
;                   done
;                -- AT THIS POINT EQUIL TABLE (no tension) LOOKS
;                   REASONABLE IN ALL PLACES (even though approximated
;                   in places). no known major problems.
;                -- add plots for tension regions
; v6.5 7/4/2008  -- Robert Marcus found that my fudge for the ice VII
;                   entropy at low pressures was a problem. I
;                   decreased the discontinuity in entropy at 251 K in
;                   the ice VII region, but there is still a small kink.
; v7.0 11/29/2008-- try to improve low temp region for vapor
;                   v7.0test -- change sublimation curve to truncate
;                               and linearly extrapolate to zero below
;                               P=1.e-20 MPa and T=100 K
;                   v7.0   -- improved ice Ih mesh resolution
;                          -- adjusted sub curve pressures to prevent
;                             very low pressure values for T<150K
;                          -- adjusted sub curve energies not to
;                             overlap with vapor
; v7.1 12/01/2008-- add zero temperature row
; v7.2 12/09/2008-- fixed mass fraction weighting for energy and
;                   entropy on phase boundaries --  was accidentally
;                                                  using volume fraction
;                   weighting
; v8.0 04/07/2010-- limit tensile strength of ice region to -100 MPa
;                   (allow for overshoot!)
;                   (estimate is about -10 MPa from Ahrens et al. DPS 2009)
;                -- limit tensile strength of liquid region to -50 MPa
;                   (estimate is about -10 MPa from Boteler &
;                   Sutherland JAP 2004)
;                -- note about triple points: in the density-temp
;                   grid, there will be a gap in entropy and energy
;                   across the phase boundary (dS and dE at constant T)
; v8.2 11/26/2013-- better grid near STP for Lock                             
; v8.3 4/9/2014  -- higher temps for Lock
; v8.3c 2/16/2023 -- added phase boundary to output file to port to
;                    ANEOS style SESAME tables and pyKO
;                   
;=========================================================================
;
; CAREFUL OF UNITS!! all the programs are different!
;
; H2O
; liquid and vapor
; NIST Steam code: WAGNER AND PRUSS, J. PHYS. CHEM. REF. DATA, 22, 783 (1993)
;
; ice Ih
; FEISTEL AND WAGNER, J. PHYS. CHEM. REF. DATA, (2005)
;
; ice VI and VII
; Stewart and Ahrens JGR 2004, Frank et al GCA 2004
;
; liquid-solid melting curve
; 0 to 3 GPa: NIST
; 3 to 60 GPa: Frank et al GCA 2004
;
; ice Ih-vapor sublimation curve
; FEISTEL AND WAGNER, J. PHYS. CHEM. REF. DATA, (2005)
;
;===============================================================
;===============================================================
pro make_h2o_ses
;
common setup, first
common boundaries, ttp6,ptp6,ttp,ptp,tcrit,pcrit,ttp3,ptp3
;
; USER SETTINGS
first=0 ; set this to 1 if starting from a cleaned directory (no *.sav files)
tension=0 ; default 1, yes tension region otherwise equilibrium boundaries used for ice Ih and liquid
;
liquid_tension_limit = -100.d0 ; MPa
iceIh_tension_limit = -100.d0 ; MPa
;
; tablename is set in the mesh specfications below
;
;
; if you already have a table *.sav file,
; then you can just plot or save a CTH mesh:
; tablename='h2o_table_v5.0'
; tablename='h2o_table_v6.4NT'
; tablename='h2o_table_v6.4'
; tablename='h2o_table_v6.2'
; tablename = 'h2o_table_test'
; tablename = 'h2o_table_test_LT'
;tablename = 'h2o_table_test'
;tablename='h2o_table_v6.5'
;tablename='h2o_table_v6.6testNT'
;tablename='h2o_table_v6.5NT'
;tablename='h2o_table_v7.0testNT'
;tablename='h2o_table_v7.0NT'
;tablename='h2o_table_v8.3bNT'
tablename='h2o_table_v8.3cNT'
;
; restore,tablename+'.sav'
; goto, plotmesh ; just plot then write txt file from sav file
; goto, writemesh ; just write the txt file from sav file
;
;===============================================================
; DEFINE REFERENCE POINTS
; Triple point: 273.16 K, 611.7 Pa (Petrenko and Whitworth 1999)
  ptp = 0.000611657d0 ; MPa Feistel and Wagner 2005
  ttp = 273.16d0 ; K Feistel and Wagner 2005
;
; Critical point: 647.096 K, 22.064 MPa, 332 kg/m3 IAPWS
; www.lsbu.ac.uk/water/data.html
  pcrit = 22.064d0 ; MPa
  tcrit = 647.096d0 ; K
  dcrit = 332.0d0 ; kg/m3
;
; Ice Ih-III-L triple point - transition from negative slope melting
;                             curve to positive slope melting curve
; (Petrenko and Whitworth 1999)
  ptp3 = 207.5d0 ; MPa
  ttp3 = ttp-22.d0 ; K
;
; Ice VI-VII-L triple point
; (Petrenko and Whitworth 1999)
  ptp6 = 2200.d0 ; MPa
  ttp6 = 81.6d0+ttp ; K
;
;===============================================================
;===============================================================
; NIST's Water program assigns entropy and energy=0 at the triple
; point.
; apply the following to shift the values to absolute values
; referenced to 0K
;
; enthalpy of fusion at tp
LmeltIh = 333.43d0  ; kJ/kg at IVL Feistel and Wagner 2005
LmeltVI = 352.701d0 ; kJ/kg at VI-VII-L 81.5 and 2150 MPa ref 535 www.lsbu.ac.uk/water/data.html
;
; before I found Feistel and Wagner 2005:
; Sshift = 3.510d0 ; kJ/kg/K Petrenko and Whitworth 1999
; Sshift = 3.522010d0 ; kJ/kg/K at 0C 1 bar 63.45 J/mol/K ref 869 www.lsbu.ac.uk/water
; Sicetp = Sshift-1.22063 ; kJ/kg/K
; Hicetp = 307.571 ; kJ/kg - integrate GAigue and Stout cp over T

; ice Ih properties at the TP (Feistel and Wagner 2005)
tmp = IhEOS(ttp,ptp) ; in K and MPa
Dicetp = tmp(0) ; kg/m3
Sicetp = tmp(1) ; kJ/kg/K
Eicetp = tmp(2) ; kJ/kg
Hicetp = tmp(3) ; kJ/kg

; liquid properties at the TP
cmd='./eostp '+strtrim(string(ttp),2)+' '+strtrim(string(ptp)) ; K MPa
spawn,cmd,res
tmp = strsplit(res(0),/extract)
if tmp(3) ne 0 then stop ; error in density calc
Dliqtp = tmp(0) ; kg/m3
Hliqtp = Hicetp+LmeltIh ; kJ/kg

; dP/dT = dS/dV
; dSiceliq = 1.22063 ; kJ/kg web page
; get second P,T,V point for calc slope dP/dT
ptp2=1.d0 ; MPa
cmd='./tmeltp '+strtrim(string(ptp2),2)
spawn,cmd,res
tmp = strsplit(res(0),/extract)
if tmp(4) ne 0 or tmp(5) ne 0 then stop ; errors in EOS
ttp2 = tmp(0)
; calc dS
dt_0 = ttp2-ttp ; K
dp_0 = (ptp2-ptp)*1.d6 ; Pa
dv_0 = 1.d0/Dliqtp - 1.d0/Dicetp ; m3/kg
ds_0 = dv_0*dp_0/dt_0 ; J/kg
dSiceliq = ds_0/1.d3 ; kJ/kg ; 1.2275766 kJ/kg

; These are the values to shift NIST to 0K reference point
;Eshift = (Hliqtp*1.e3-ttp*(ptp*1.d6))/1.e3 ; kJ/kg  465.03426 ; E = H-PV
Eshift = Hliqtp ; kJ/kg ; this works
Sshift = Sicetp+dSiceliq ; kJ/kg 3.5233498
print,'Eshift,Sshift=',Eshift,Sshift
;stop
;===============================================================
;===============================================================
; SESAME TABLE MESH SPECIFICATION
;
; mesh: temperature and density
;
;  test low T ice Ih, VI, VII
if 0 then begin
  tmarkers =    [10.,300.]
  tnum     = 1.*[40.]
  ttype    =    [1]
  dmarkers =    [1200.,5000.] ; kg/m3, boundaries for meshing
  dnum     = 1.*[19] ; number of mesh points between markers
  dtype    =    [1,1] ; type of meshing: 1 for linear, 2 for log

  tablename = 'h2o_table_test'
end
;  test liquid tension v6.2, v6.4
if 0 then begin
  tmarkers =    [274.,700.]
  tnum     = 1.*[40.]
  ttype    =    [1]
  dmarkers =    [300.,1500.] ; kg/m3, boundaries for meshing
  dnum     = 1.*[40] ; number of mesh points between markers
  dtype    =    [1] ; type of meshing: 1 for linear, 2 for log

  tablename = 'h2o_table_test_LT'
end

; test 1V mixing v6.0
if 0 then begin
  tmarkers =    [10.,150.d0,250.d0,270.d0,276.]
  tnum     = 0.5*[6,  16,    6,     6]
  ttype    =     [1,  1,     1,     1]
  dmarkers =    [1.,500.,1000.] ; kg/m3, boundaries for meshing
  dnum     = 1.*[10,   4] ; number of mesh points between markers
  dtype    =    [1,   1] ; type of meshing: 1 for linear, 2 for log

  tablename = 'h2o_table_test_1V'
end
; ; test coverage 
if 0 then begin
   tmarkers =    [10.,150.d0,250.d0,270.d0,276.,355.,647.,1000.d0,5000.d0,2.d4,2.d5] ; K, boundaries for meshing
  tnum     = 0.5*[6,  16,    6,     6,     10,   20,   20,  10,      10,   4] ; number of mesh points between markers
  ttype    =    [1,  1,     1,     1,     1,     1,   1,   1,      1,    2] ; type of meshing: 1 for linear, 2 for log
;  dmarkers =    [2000.,3000.,5000.] ; kg/m3, boundaries for meshing
  dmarkers =    [1.,500.,5000.] ; kg/m3, boundaries for meshing
  dnum     = 1.*[20,   20] ; number of mesh points between markers
  dtype    =    [1,   1] ; type of meshing: 1 for linear, 2 for log

  tablename = 'h2o_table_test'
end
;
; VERSION 5.0 and 6.0, v6.2, v6.3, v6.4, v6.5
if 0 then begin
   tmarkers =    [10.,150.d0,250.d0,270.d0,276.,355.,647.,1000.d0,5000.d0,2.d4,2.d5] ; K, boundaries for meshing
  tnum     = 2.*[6,  16,    6,     6,     10,   20,   20,  10,      10,   4] ; number of mesh points between markers
  ttype    =    [1,  1,     1,     1,     1,     1,   1,   1,      1,    2] ; type of meshing: 1 for linear, 2 for log
  dmarkers =    [0.0, 1.e-9,100.,917,930.,1200.,1400.,2000.,3000.,5000.] ; kg/m3, boundaries for meshing
;  dnum     = .5*[2,   22,   8,   16,  10,  16,   20,   20,   6] ; number of mesh points between markers
  dnum     = 2.*[1,   22,   8,   16,  10,  16,   20,   20,   6] ; number of mesh points between markers
  dtype    =    [1,    2,    1,   1,   1,   1,   1,    1,   1] ; type of meshing: 1 for linear, 2 for log
  tablename = 'h2o_table_v6.5'
end

; VERSION 6.6, 7.0 test
if 0 then begin
  tmarkers =    [250.d0,270.d0,276.,300]
  tnum     = 1.*[6,     6,     10]
  ttype    =    [1,     1,     1]
;  dmarkers =    [0.0, 1.e-9,100.,917,930.,1200.]
;  dnum     = 1.*[1,   22,   8,   16,  10]
;  dtype    =    [1,    2,    1,   1,   1]
  dmarkers =    [1.e-10, 100.,910,950.,1200.] ; 1.e-10 lower limit for density
  dnum     = 1.*[5, 5,  5,  5]
  dtype    =    [2,  1,   1,   1]
  tablename = 'h2o_table_v7.2testNT'
end

; VERSION 7.0, 7.1, 7.2, 8.0
if 0 then begin
   tmarkers =    [10.,150.d0,250.d0,270.d0,276.,355.,647.,1000.d0,5000.d0,2.d4,2.d5] ; K, boundaries for meshing
  tnum     = 2.*[6,  16,    6,     6,     10,   20,   20,  10,      10,   4] ; number of mesh points between markers
  ttype    =    [1,  1,     1,     1,     1,     1,   1,   1,      1,    2] ; type of meshing: 1 for linear, 2 for log
  dmarkers =    [0.0, 1.e-10,100.,910,950.,1200.,1400.,2000.,3000.,5000.] ; kg/m3, boundaries for meshing
  dnum     = 2.*[1,   22,    22,  20, 10,  16,   20,   20,   6] ; number of mesh points between markers
  dtype    =    [1,    2,    1,   1,   1,   1,   1,    1,   1] ; type of meshing: 1 for linear, 2 for log
;  tablename = 'h2o_table_v7.0NT'
;  tablename = 'h2o_table_v7.1'
;  tablename = 'h2o_table_v7.2'
  tablename = 'h2o_table_v8.0NT'
end

; VERSION 8.2
if 0 then begin
   tmarkers =    [10.,150.d0,250.d0,270.d0,276.,300.,355.,647.,1000.d0,5000.d0,2.d4,2.d5] ; K, boundaries for meshing
  tnum     = 2.*[6,  16,    6,     6,     10,   20,   20,  10,      10,   4] ; number of mesh points between markers
  ttype    =    [1,  1,     1,     1,     1,     1,   1,   1,      1,    2] ; type of meshing: 1 for linear, 2 for log
  dmarkers =    [0.0, 1.e-10,100.,910,950.,995,1023.,1200.,1400.,2000.,3000.,5000.] ; kg/m3, boundaries for meshing
  dnum     = 2.*[1,   22,    22,  20, 15,  14, 16,   16,   20,   20,   6] ; number of mesh points between markers
  dtype    =    [1,    2,    1,   1,   1,   1, 1,    1,    1,    1,   1] ; type of meshing: 1 for linear, 2 for log
  tablename = 'h2o_table_v8.2NT'
end

; VERSION 8.3
if 0 then begin
;   tmarkers =    [10.,150.d0,250.d0,270.d0,276.,300.,355.,647.,1000.d0,5000.d0,2.d4,6.d5] ; K, boundaries for meshing
;  tnum     = 2.*[6,  16,    6,     6,     10,   20,   20,  10,      10,   10, 10] ; number of mesh points between markers
;  ttype    =    [1,  1,     1,     1,     1,     1,   1,   1,      1,    2, 2] ; type of meshing: 1 for linear, 2 for log
   tmarkers =    [10.,150.d0,250.d0,270.d0,276.,300.,355.,647.,1000.d0,5000.d0,1.d5] ; K, boundaries for meshing
  tnum     = 2.*[6,  16,    6,     6,     10,   20,   20,  10,      10,   10] ; number of mesh points between markers
  ttype    =    [1,  1,     1,     1,     1,     1,   1,   1,      1,    2] ; type of meshing: 1 for linear, 2 for log
  dmarkers =    [0.0, 1.e-10,100.,910,950.,995,1023.,1200.,1400.,2000.,3000.,5000.] ; kg/m3, boundaries for meshing
  dnum     = 2.*[1,   22,    22,  20, 15,  14, 16,   16,   20,   20,   6] ; number of mesh points between markers
  dtype    =    [1,    2,    1,   1,   1,   1, 1,    1,    1,    1,   1] ; type of meshing: 1 for linear, 2 for log
  tablename = 'h2o_table_v8.3bNT'
end

; VERSION 8.3c
if 1 then begin
;   tmarkers =    [10.,150.d0,250.d0,270.d0,276.,300.,355.,647.,1000.d0,5000.d0,2.d4,6.d5] ; K, boundaries for meshing
;  tnum     = 2.*[6,  16,    6,     6,     10,   20,   20,  10,      10,   10, 10] ; number of mesh points between markers
;  ttype    =    [1,  1,     1,     1,     1,     1,   1,   1,      1,    2, 2] ; type of meshing: 1 for linear, 2 for log
   tmarkers =    [10.,150.d0,250.d0,270.d0,276.,300.,355.,647.,1000.d0,5000.d0,1.d5] ; K, boundaries for meshing
  tnum     = 4.*[6,  16,    6,     6,     10,   20,   20,  10,      10,   10] ; number of mesh points between markers
  ttype    =    [1,  1,     1,     1,     1,     1,   1,   1,      1,    2] ; type of meshing: 1 for linear, 2 for log
  dmarkers =    [0.0, 1.e-10,100.,910,950.,995,1023.,1200.,1400.,2000.,3000.,10000.] ; kg/m3, boundaries for meshing
  dnum     = 4.*[1,   22,    22,  20, 15,  14, 16,   16,   20,   20,   20] ; number of mesh points between markers
  dtype    =    [1,    2,    1,   1,   1,   1, 1,    1,    1,    1,   1] ; type of meshing: 1 for linear, 2 for log
  tablename = 'h2o_table_v8.3cNT'
end

; VERSION 4.3 MESH
if 0 then begin
  tmarkers =    [10.,150.d0,250.d0,270.d0,276.,355.,647.,1000.d0,5000.d0,2.d4,2.d5] ; K, boundaries for meshing
  tnum     = 2.*[6,  16,    6,     6,     10,   20,   20,  10,      10,   4] ; number of mesh points between markers
  ttype    =    [1,  1,     1,     1,     1,     1,   1,   1,      1,    2] ; type of meshing: 1 for linear, 2 for log
  dmarkers =    [0.0, 1.e-9,100.,917,930.,1200.,1400.,2000.,3000.,5000.] ; kg/m3, boundaries for meshing
  dnum     = 2.*[2,   22,   8,   16,  10,  16,   20,   20,   6] ; number of mesh points between markers
  dtype    =    [1,    2,    1,   1,   1,   1,   1,    1,   1] ; type of meshing: 1 for linear, 2 for log
  tablename = 'h2o_table_v5.0_v4.3mesh'
end

; FIRST RUN - don't make a huge mesh
;if 1 then begin
if first then begin
  tmarkers = [200.,300.]
  tnum     = [2]
  ttype    = [1]
  dmarkers = [1200,1500.]
  dnum     = [4]
  dtype    = [1]
  tablename = 'FIRSTRUN'
end

zzz = where(tnum lt 1.)
if zzz(0) ne -1 then tnum(zzz) = 1.
zzz = where(dnum lt 1.)
if zzz(0) ne -1 then dnum(zzz) = 1.

  tsize = total(tnum) +1
  dsize = total(dnum) +1

  tarr = dblarr(tsize)
  darr = dblarr(dsize)

  phasearr = strarr(tsize,dsize)
  testarr = intarr(tsize,dsize)
  parr = dblarr(tsize,dsize) ; pressure in MPa
  earr = dblarr(tsize,dsize) ; internal energy in kJ/kg
  sarr = dblarr(tsize,dsize) ; specific entropy in kJ/kg/K
  tijarr = dblarr(tsize,dsize) ; specific entropy in kJ/kg/K
  dijarr = dblarr(tsize,dsize) ; specific entropy in kJ/kg/K

  ind = 0
  for i=0,n_elements(tnum)-1 do begin
      if i gt 0 then ind = ind+tnum(i-1)
      if ttype(i) eq 1 then  begin
          step = (tmarkers(i+1)-tmarkers(i))/tnum(i)
          tarr(ind:ind+tnum(i)-1) = findgen(tnum(i))*step+tmarkers(i)
      end
      if ttype(i) eq 2 then begin
          step = (alog10(tmarkers(i+1))-alog10(tmarkers(i)))/tnum(i)
          tarr(ind:ind+tnum(i)-1) = 10.d0^(findgen(tnum(i))*step+alog10(tmarkers(i)))
      end
  end
  tarr(tsize-1)=tmarkers(i); add the last point

  ind = 0
  for i=0,n_elements(dnum)-1 do begin
      if i gt 0 then ind = ind+dnum(i-1)
      if dtype(i) eq 1 then  begin
          step = (dmarkers(i+1)-dmarkers(i))/dnum(i)
          darr(ind:ind+dnum(i)-1) = findgen(dnum(i))*step+dmarkers(i)
      end
      if dtype(i) eq 2 then begin
          step = (alog10(dmarkers(i+1))-alog10(dmarkers(i)))/dnum(i)
          darr(ind:ind+dnum(i)-1) = 10.d0^(findgen(dnum(i))*step+alog10(dmarkers(i)))
      end
  end
  darr(dsize-1)=dmarkers(i); add the last point

  ; check the mesh
  if 1 and not(first) then begin
      print,tarr
      print,darr
      print,n_elements(tarr),n_elements(darr)

      window,0
      !p.multi=[0,1,3]
      !p.charsize=2.
      plot,tarr,psym=-1,title='MESH TEMP',xtit='index',ytit='Temp (K)'
      plot,darr,psym=-1,title='MESH DENS LIN',xtit='index',ytit='Density (kg/m3)'
      plot,darr,psym=-1,/ylog,title='MESH DENS LOG',xtit='index',ytit='Density (kg/m3)'
      print,'check the mesh for table ',tablename
      print,'then continue'
      print,'tension = ',tension
      stop
  end

if first then print,'FIRST PASS. MAKING INITIAL CALCS (BE PATIENT):'

;===============================================================
;===============================================================
; CALC THE MELTING CURVE AS A FN OF PRESSURE: LIQUID-SOLID BOUNDARY
; energy and entropy need to be fixed v5
nneg = 21
if first then begin
;if 1 then begin
    print,'calculating liquid-solid melting curve'
;stop
;    pmeltall =(findgen(1000))*100.+ptp ; MPa ; up to 100 GPa
    pmeltall_tmp =[findgen(nneg)/(nneg-1)*(ptp3-ptp)+ptp, (findgen(1000)+1)*100.+ptp3] ; MPa ; up to 100 GPa
    zzz = where(pmeltall_tmp le 1.e5) ; 100 GPa
    pmeltall = pmeltall_tmp(zzz)
    zzz = where(pmeltall lt ptp)
    if zzz(0) gt -1 then pmeltall(zzz) = ptp
    dmeltall = dblarr(n_elements(pmeltall)) ; kg/m3
    tmeltall = dblarr(n_elements(pmeltall)) ; K
    smeltall = dblarr(n_elements(pmeltall)) ; kJ/kg/K
    emeltall = dblarr(n_elements(pmeltall)) ; kJ/kg
    dice1meltall = dblarr(nneg); kg/m3
    sice1meltall = dblarr(nneg) ; kJ/kg/K
    eice1meltall = dblarr(nneg) ; kJ/kg
    for i=0,n_elements(pmeltall)-1 do begin
        if pmeltall(i) lt 3000. then begin ; ORIG
;        if pmeltall(i) lt 5000. then begin ; v5 - not a big diff
            ; use NIST melting curve
            cmd='./tmeltp '+strtrim(String(pmeltall(i)))
            spawn,cmd,res
            tmp = strsplit(res(0),/extract)
            if tmp(4) ne 0 or tmp(5) ne 0 then stop ; errors in EOS
            tmeltall(i) = tmp(0)
            dmeltall(i) = tmp(1) ; liquid
            emeltall(i) = tmp(2)+Eshift ; liquid
            smeltall(i) = tmp(3)+Sshift ; liquid
;            stop
        end else begin
                                ; if P>3GPa use ice VII melting curve
                                ; from Frank et al GCA 2004
                                ; good to 60 GPa
                                ; P in GPa
            tmeltall(i) = $
              (((pmeltall(i)/1.e3 - 2.17)/0.764 + 1.d0)^(1./4.32))*355.d0

            cmd='./eostp '+strtrim(String(tmeltall(i)))+' '+strtrim(String(pmeltall(i)))
            spawn,cmd,res
            tmp = strsplit(res(0),/extract)
            if tmp(3) ne 0 then stop ; error in density calc
            dmeltall(i) = tmp(0) ; liquid
            emeltall(i) = tmp(1)+Eshift ; liquid
            smeltall(i) = tmp(2)+Sshift ; liquid
        end
        if i lt nneg then begin
            ; also calc ice Ih values
            tmpi = IhEOS(tmeltall(i),pmeltall(i)) ; wants K and MPa
            dice1meltall(i) = tmpi(0)
            sice1meltall(i) = tmpi(1)
            eice1meltall(i) = tmpi(2)
        end
    end
    save,filename='meltcurve.sav',tmeltall,dmeltall,pmeltall,smeltall,emeltall,dice1meltall,sice1meltall,eice1meltall,nneg
end else restore,'meltcurve.sav'
;stop
;===============================================================
;===============================================================
; CALC THE SUBLIMATION CURVE AS A FN OF TEMPERATURE: VAPOR-SOLID BOUNDARY
; Wagner and Pruss sublimation curve is only good >190 K (psubt)
; Use Feistel and Wagner 2005 sublimation curve for ice Ih since it is newer
;if 1 then begin
if first then begin
    print,'calculating vapor-solid sublimation curve'

;    tsuball =(findgen(273)+1.1) ; K
    tsuball =(findgen(273)+1.1) ; K
    psuball = dblarr(n_elements(tsuball)) ; MPa vapor (not used - use ice instead)
    dsuball = dblarr(n_elements(tsuball)) ; kg/m3 vapor
    ssuball = dblarr(n_elements(tsuball)) ; kJ/kg/K vapor
    esuball = dblarr(n_elements(tsuball)) ; kJ/kg vapor
    dicesuball = dblarr(n_elements(tsuball)) ; kg/m3
    picesuball = dblarr(n_elements(tsuball)) ; MPa
    sicesuball = dblarr(n_elements(tsuball)) ; kJ/kg/K
    eicesuball = dblarr(n_elements(tsuball)) ; kJ/kg

    for i=n_elements(tsuball)-1,0,-1 do begin
        ; calc ice Ih sublimation pressure, then get EOS
        picesuball(i) = icesub(tsuball(i)) ; returns MPa, limited to 1.e-20 MPa
        tmpi = IhEOS(tsuball(i),picesuball(i)) ; wants K and MPa
        dicesuball(i) = tmpi(0)
        sicesuball(i) = tmpi(1) ; already in absolute scale
        eicesuball(i) = tmpi(2) ; already in absolute scale

        ; get vapor side from liquid/steam EOS
        cmd='./eostp '+strtrim(String(tsuball(i)))+' '+strtrim(String(picesuball(i)))
        spawn,cmd,res
        tmp = strsplit(res(0),/extract)
        if tmp(3) ne 0 then begin
            ; the EOS is out of range for T<=100K
            print,'line 345 - VAPOR-ICE SUB eostp out of range: ',tsuball(i),picesuball(i),tmp(0),tmp(1),tmp(2)
            ;stop                ; error in density calc
        end
        dsuball(i) = tmp(0)    ; vapor
        esuball(i) = tmp(1)+Eshift ; vapor
        ssuball(i) = tmp(2)+Sshift ; vapor
    end

;stop
    ; there are bogus vapor values at 100 K and lower

   if 0 then begin
    ; v5.0 - THIS DOES NOT WORK - ENTROPIES ARE TOO HIGH
    ; interpolate over these values
    zzzbad = findgen(101)
    zzzgood = findgen(100)+100.
    esuball(zzzbad) = interpol(esuball(zzzgood),tsuball(zzzgood),tsuball(zzzbad))
    ssuball(zzzbad) = interpol(ssuball(zzzgood),tsuball(zzzgood),tsuball(zzzbad))
    dsuball(zzzbad) = interpol(dsuball(zzzgood),tsuball(zzzgood),tsuball(zzzbad))
                                ; THIS IS OBVIOUSLY TOO LARGE OF DENSITIES _
                                ; SHOULD BE EXPONENTIAL BUT SHOULD BE GOOD ENOUGH
    zzzgood = where(dsuball gt 1.e-5,complement=zzzbad)
;    dsuball(zzzbad) = 10.d0^interpol(alog10(dsuball(zzzgood)),tsuball(zzzgood),tsuball(zzzbad))
;    dsuball(zzzbad) =
;    10.d0^interpol(alog10(dsuball(zzzgood)),tsuball(zzzgood),tsuball(zzzbad))
     end ; v5.0 attempt

    ; v6.0 attempt
    ; fix low temperature values to the last reasonable value from EOSTP
    zzzbad = findgen(201)
    ; v6.0 -- energies are very slightly too high if help constant here    
;    esuball(zzzbad) = esuball(201) ; 2908.5291
    ; v6.4 fix lowest value to minimum in the vapor field - linear to 200K
    ; esuball(zzzbad) = esuball(201) - (esuball(201)-2862.d0)*(200.-tsuball(zzzbad))/200.; 
    ; v7.0
                                ; v7.0 -- there was a slight overlap
                                ;         with esuball and vapor field
                                ;         so change the slope slightly
    esuball(zzzbad) = esuball(201) - (esuball(201)-2862.d0)*(200.-tsuball(zzzbad))/35.; 
    ; v6.0
    ; leave pressure (picesuball) alone
    ; v7.0
    ; limit low values for pressure
    ; linearly extrapolate to zero below 150 K
      zzz = where(tsuball le 150.)
      picesuball(zzz) = (picesuball(max(zzz)))/150.*tsuball(zzz)
    ; v6.0
    ; density points are bad below 200 K
    ; use log-linear extrapolation
    ; this looks ok
    zzzgood = findgen(50)+210
    zzzbad = findgen(210)
    dsuball(zzzbad) = 10.d0^(interpol(alog10(dsuball(zzzgood)),tsuball(zzzgood),tsuball(zzzbad)))
    ;
    ; v6.0 -- entropies are too low if held constant
;DONT USE    ssuball(zzzbad) = ssuball(201)
    ; v6.4
    ; try log-linear extrapolation for entropy
    ; want to get close to the liquid side of the curve
    zzzgood = findgen(50)+210
    zzzbad = findgen(210)
    ssuball(zzzbad) = 10.d0^(interpol(alog10(ssuball(zzzgood)),tsuball(zzzgood),tsuball(zzzbad)))

    !p.multi=[0,1,4]
    plot,tsuball,picesuball,/ylog,tit='p'
    plot,tsuball,dsuball,/ylog,tit='d'
    plot,tsuball,esuball,/yno,tit='e'
    plot,tsuball,ssuball,tit='s'
;    stop

        ; vapor side
;        stop

        ; old code to run if only using the Wagner and Pruss equations
;        ; psubt is not valid below 190 K
;        if tsuball(i) ge 100 then begin
;            cmd='./psubt '+strtrim(string(tsuball(i)))
;            spawn,cmd,res
;            tmp = strsplit(res(0),/extract)
;            print,tsuball(i),tmp
;            err1 = tmp(4)
;            err2 = tmp(5)
;            psuball(i) = tmp(0) ; vapor
;            dsuball(i) = tmp(1) ; vapor
;            esuball(i) = tmp(2)+Eshift ; vapor
;            ssuball(i) = tmp(3)+Sshift ; vapor
;        end else begin
;            psuball(i) = tmp(0) ; vapor
;            dsuball(i) = tmp(1) ; vapor
;            esuball(i) = tmp(2)+Eshift ; vapor
;            ssuball(i) = tmp(3)+Sshift ; vapor
;        end

    save,filename='subcurve.sav',tsuball,dsuball,picesuball,ssuball,esuball,dicesuball,picesuball,sicesuball,eicesuball
end else restore,'subcurve.sav'
;stop
;stop
;===============================================================
;===============================================================
; CALC THE ICE Ih-ICE VI BOUNDARY AS A FN OF TEMPERATURE
; Feistel and Wagner 2005
; assume boundary is at constant pressure at the ice Ih-III-L triple point
if first then begin
;if 1 then begin
    print,'calculating ice Ih-ice VI boundary'

    t16all =(findgen(251)+1.16) ; K
;    p16all = dblarr(n_elements(t16all))+ptp3 ; MPa
    p16all = dblarr(n_elements(t16all))+b16(t16all) ; MPa
    d16all = dblarr(n_elements(t16all)) ; kg/m3 Ih
    s16all = dblarr(n_elements(t16all)) ; kJ/kg/K Ih
    e16all = dblarr(n_elements(t16all)) ; kJ/kg Ih

    for i=n_elements(t16all)-1,0,-1 do begin
        tmpi = IhEOS(t16all(i),p16all(i)) ; wants K and MPa
        d16all(i) = tmpi(0)
        s16all(i) = tmpi(1) ; already in absolute scale
        e16all(i) = tmpi(2) ; already in absolute scale
    end
    save,filename='ice16curve.sav',t16all,d16all,p16all,s16all,e16all
end else restore,'ice16curve.sav'
;stop
;===============================================================
;===============================================================
; CALC THE SATURATION VAPOR CURVE AS A FN OF TEMPERATURE: VAPOR-LIQUID BOUNDARY
if first then begin
    print,'calculating saturation vapor curve (vapor-liquid boundary)'

    tsatall = [findgen(tcrit-ttp)+ttp,tcrit-.001] ; K ; lower highest value slightly for psatt **
    dlsatall = dblarr(n_elements(tsatall)) ; kg/m3
    dvsatall = dblarr(n_elements(tsatall)) ; kg/m3
    psatall = dblarr(n_elements(tsatall)) ; MPa
    elsatall = dblarr(n_elements(tsatall))
    evsatall = dblarr(n_elements(tsatall))
    slsatall = dblarr(n_elements(tsatall))
    svsatall = dblarr(n_elements(tsatall))
    for i=n_elements(tsatall)-1,0,-1 do begin
        cmd='./psatt '+strtrim(string(tsatall(i)))
        spawn,cmd,res
        tmp = strsplit(res,/extract)
        if tmp(5) ne 0 then stop ; ERROR IN SAT CURVE
        psatall(i) = tmp(0)     ; L-V
        dlsatall(i) = tmp(1)     ; liquid
        dvsatall(i) = tmp(4)     ; vapor
        ;print,tsatall(i),tmp
        ; get e and s for liquid side
        cmd= './eostd '+strtrim(String(tsatall(i)),2)+' '+strtrim(string(dlsatall(i)),2)
        spawn,cmd,res
        tmp = strsplit(res,/extract)
        if tmp(1) eq 0.0 and tmp(2) eq 0.0 then stop ; ERROR IN EOS
        elsatall(i)= tmp(1)+Eshift
        slsatall(i)= tmp(2)+Sshift
        ; get e and s for vapor side
        cmd= './eostd '+strtrim(String(tsatall(i)),2)+' '+strtrim(string(dvsatall(i)),2)
        spawn,cmd,res
        tmp = strsplit(res,/extract)
        if tmp(1) eq 0.0 and tmp(2) eq 0.0 then stop ; ERROR IN EOS
        evsatall(i)= tmp(1)+Eshift
        svsatall(i)= tmp(2)+Sshift

    end
;    stop
    save,filename='satcurve.sav',tsatall,dlsatall,dvsatall,psatall,$
      evsatall,elsatall,svsatall,slsatall
end else restore,'satcurve.sav'

;===================================================================
;===================================================================
; CALC TABLE POINTS
; check for degeneracy between 250 and 273 K
if darr(0) eq 0.0 then begin ; v6.4
   jstart = 1 
   tijarr(*,0) = tarr(*) ; fill in tijarr zero density column with temps
end else jstart = 0
  for i = 0, tsize-1 do begin
      print,i,'temp=',tarr(i)
      for j = jstart, dsize-1 do begin
          ; work our way up in temperature to check for different phases
          notdone=1

          tijarr(i,j) = tarr(i)
          dijarr(i,j) = darr(j)

          ;-------------------------------------------------------------
          ; ------------- T<=TP -------------------------------------
          ;-------------------------------------------------------------
          ; below the triple point
          if (tarr(i) le ttp) and notdone then begin

              ; -------------- VAPOR ---------------------------
                                ; A. check if density less than vapor
                                ; on sublimation curve
              dvapsub = interpol(dsuball,tsuball,tarr(i))
              dicesub = interpol(dicesuball,tsuball,tarr(i))
              if darr(j) le dvapsub and notdone then begin
                  ; T<TP, D<Dvapmax --> Vapor
                  cmd= './eostd '+strtrim(String(tarr(i)),2)+' '+strtrim(string(darr(j)),2)
                  spawn,cmd,res
                  tmp = strsplit(res,/extract)
                  if tmp(1) eq 0.0 and tmp(2) eq 0.0 then print,'line 494 EOSTD: ',tarr(i),darr(j) ; ERROR IN EOS
                  phasearr(i,j)= ' V'
                  if darr(j) eq 0.d0 then parr(i,j)=0.d0 else parr(i,j)= tmp(0)
                  earr(i,j)= tmp(1)+Eshift
                  sarr(i,j)= tmp(2)+Sshift
                  ;print,cmd, res

                  notdone=0
              end


              ; -------------- V-Ih ---------------------------
              ; V3: This is the tension region below the triple point
              if darr(j) gt dvapsub and darr(j) le dicesub and notdone then begin
                                ; T<TP, Dvapmax<D<Dice1min --> 
                                ; Ice Ih in tension 
                                ; OR 
                                ; V-Ih equilibrium
                                ; calc E and S by mixing properties
                                ; weighted by mass fraction
                  if tension then begin
                     if darr(j) gt dicesub/2.d0 then begin
                                ; format table as ice Ih in tension
                                ; use Birch-Murnaghan equation to
                                ; determine tension region
                      kt0 = 8.4d9 ; Pa stewart and ahrens
                      ktp = 3.0d0
                      v0 = 1./dicesub ; mks
                      vneg = 1./darr(j) ; mks
                      xarr  = (vneg/v0)^(1./3.d0)
                      pneg = 3.d0*kt0*(1.-xarr)*xarr^(-2) * exp(3./2.*(ktp-1.)*(1.-xarr)) ; 3rd order B-M eq.

                      ; v8.0
                      ; only fill the array if pneg < limit
                      if pneg/1.e6 gt iceIh_tension_limit then begin
                        parr(i,j) = pneg/1.e6
                                  ;tmpi = IhEOS(tarr(i),pneg)  ; does not work
                        earr(i,j) = -1.d0
                        sarr(i,j) = -1.d0
                        phasearr(i,j)= '1T' ; v6.2, previously V1 flag
                        notdone=0
                      end         ; else must be ice-vapor region

                   end ; right density for tension region
                  end ; check for tension flag
                  if notdone then begin
                      ; must be in the ice-vapor mix region
                      ; format table as ice-vapor mix
                      phasearr(i,j)= 'V1'
                      notdone=0
;volume               fraction = (darr(j)-dvapsub) / (dicesub-dvapsub) ; of ice
;v7.2 fixed to mass fraction of ice
                      vfraction = (darr(j)-dvapsub) / (dicesub-dvapsub) ; of ice
                      fraction = vfraction*dicesub / $
                                 (vfraction*dicesub + (1.d0-vfraction)*dvapsub)
;print,fraction,vfraction,darr(j),dicesub,dvapsub
;stop
                     ;print,fraction,darr(j), (1.-fraction)*dvapsub+(fraction)*dicesub

                      earr(i,j) = (1.-fraction)*interpol(esuball,tsuball,tarr(i))+$
                        fraction *interpol(eicesuball,tsuball,tarr(i))
                      sarr(i,j) = (1.-fraction)*interpol(ssuball,tsuball,tarr(i))+$
                        fraction *interpol(sicesuball,tsuball,tarr(i))
                      parr(i,j) = interpol(picesuball,tsuball,tarr(i))
;                                stop

                   end          ; end V1 mix region

               end ; end V1 mix density range

              ; -------------- Ih ---------------------------------
              ; check for pure ice Ih
                                ; if temp is < ice III TP, then max
                                ; density on iceI-VI boundary
; v5.0 and older, v6.1
              dice16max = interpol(d16all,t16all,tarr(i))
; v6.0 extend ice Ih to the ice VI density
;              p6min = b16(tarr(i)) ; now small slope
;              tmp6 = icevi(tarr(i),p6min) ; in K and MPa
;              dice16max = tmp6(0) ; density of ice VI on Ih-VI boundary at tarr(i)

              if tarr(i) le ttp3 then d1max = dice16max
                                ; if TTP3<T<TTP, then max density on
                                ; ice I-Liquid boundary
              dice1melt = interpol(dice1meltall,tmeltall(0:nneg-1),tarr(i))
              if tarr(i) gt ttp3 then d1max = dice1melt
              ; Check if density is greater than
              ; ice Ih on sublimation curve and lt d1max
              if darr(j) gt dicesub and darr(j) le d1max and notdone then begin
                  ; T<TP, dice1min<D<Dice1max --> pure ice Ih
                  pfit = [1.e-5,1.e-4,1.e-3,1.e-2,0.1,1.,10.,100.,150.,ptp3]
                  tmpi = IhEOS(tarr(i),pfit)
                  dfit = tmpi(0:n_elements(pfit)-1)
                  pint = interpol(pfit,dfit,darr(j))
                  phasearr(i,j) = ' 1'
                  tmpi = IhEOS(tarr(i),pint)
                  parr(i,j) = pint
                  earr(i,j) = tmpi(2)
                  sarr(i,j) = tmpi(1)
                  notdone=0
;                  print,'Ih:',tarr(i),darr(j),pint,tmpi
;                    stop
              end

              ; -------------- Ih-VI ---------------------------------
              ;p6min = ptp3      ; using ice III boundary pressure with no temp dependence
              p6min = b16(tarr(i)) ; now small slope
              tmp6 = icevi(tarr(i),p6min) ; in K and MPa
              d61min = tmp6(0) ; density of ice VI on Ih-VI boundary at tarr(i)
              if tarr(i) le ttp3 and darr(j) gt d1max and darr(j) le d61min and notdone then begin
                                ; T<TP3, Dice1max<D<Dice6min
                                ; Ih-VI boundary

                  s61min  = tmp6(1)
                  e61min  = tmp6(2)
                  eice16max = interpol(e16all,t16all,tarr(i))
                  sice16max = interpol(s16all,t16all,tarr(i))
; v7.2 volume and mass fraction
                  vfraction = (darr(j)-d1max) / (d61min-d1max) ; of ice 6
                  fraction = vfraction*d61min / $
                             (vfraction*d61min + (1.d0-vfraction)*d1max)

                  parr(i,j) = p6min
                  earr(i,j) = (1.-fraction)*eice16max+fraction*e61min
                  sarr(i,j) = (1.-fraction)*sice16max+fraction*s61min
                  phasearr(i,j) = '16'
                  ;print,tarr(i),parr(i,j),darr(j),fraction,earr(i,j),sarr(i,j)
                  notdone=0

              end

              ; -------------- Ih-Liquid ---------------------------------
              dliqmelt = interpol(dmeltall(0:nneg-1),tmeltall(0:nneg-1),tarr(i))
              if tarr(i) gt ttp3 and darr(j) gt dice1melt and darr(j) le dliqmelt and notdone then begin
                  ; mix Ih-liquid below the TP and above the ice III TP
                  pliqmelt = interpol(pmeltall(0:nneg-1),tmeltall(0:nneg-1),tarr(i))
                  eliqmelt = interpol(emeltall(0:nneg-1),tmeltall(0:nneg-1),tarr(i));+Eshift
                  sliqmelt = interpol(smeltall(0:nneg-1),tmeltall(0:nneg-1),tarr(i));+Sshift
                  eice1melt = interpol(eice1meltall(0:nneg-1),tmeltall(0:nneg-1),tarr(i))
                  sice1melt = interpol(sice1meltall(0:nneg-1),tmeltall(0:nneg-1),tarr(i))

;v7.2 volume and mass fraction
                  vfraction = (darr(j)-dice1melt) / (dliqmelt-dice1melt) ; of liquid
                  fraction = vfraction*dliqmelt / $
                             (vfraction*dliqmelt + (1.d0-vfraction)*dice1melt)

                  parr(i,j) = pliqmelt
                  earr(i,j) = (1.-fraction)*eice1melt+fraction*eliqmelt
                  sarr(i,j) = (1.-fraction)*sice1melt+fraction*sliqmelt
                  phasearr(i,j) = '1L'
             ;         print,tarr(i),pliqmelt,darr(j),fraction,earr(i,j),sarr(i,j)
                  notdone=0
             ;         stop
              end

             ; ---------------- Liquid -------------------------------
             ; next pure liquid below TP
              ; density of liquid on the L-VI boundary
              dliq6melt = interpol(dmeltall(nneg:*),tmeltall(nneg:*),tarr(i))
              if tarr(i) gt ttp3 and darr(j) gt dliqmelt and darr(j) le dliq6melt and notdone then begin
                                ; liquid below the TP and above the ice III TP
                  cmd= './eostd '+strtrim(String(tarr(i)),2)+' '+strtrim(string(darr(j)),2)
                  spawn,cmd,res
                  tmp = strsplit(res,/extract)
                  if tmp(1) eq 0.0 and tmp(2) eq 0.0 then print,'line 628 EOSTD: ',tarr(i),darr(j) ; ERROR IN EOS
;                  if tmp(1) eq 0.0 and tmp(2) eq 0.0 then stop ; ERROR IN EOS
                  phasearr(i,j)= ' L'
                  parr(i,j)= tmp(0)
                  earr(i,j)= tmp(1)+Eshift
                  sarr(i,j)= tmp(2)+Sshift
                  ;print,cmd, res
                  notdone=0
                  ;stop
              end

             ; ---------------- Liquid - VI -------------------------------
              p6melt = interpol(pmeltall(nneg:*),tmeltall(nneg:*),tarr(i))
              tmp6 = icevi(tarr(i),p6melt) ; in K and MPa
              dice6melt = tmp6(0) ; density of ice VI on the L-VI boundary at tarr(i)
              if tarr(i) gt ttp3 and darr(j) gt dliq6melt and darr(j) le dice6melt and notdone then begin
                                ; mix liquid-ice 6 below the TP and above the ice III TP
               ;   stop

                  ; mix liquid-ice VI below the TP and above the ice III TP
                  dliqmelt = interpol(dmeltall(nneg:*),tmeltall(nneg:*),tarr(i)) ; new v5.0
                  pliqmelt = interpol(pmeltall(nneg:*),tmeltall(nneg:*),tarr(i))
                  eliqmelt = interpol(emeltall(nneg:*),tmeltall(nneg:*),tarr(i));+Eshift
                  sliqmelt = interpol(smeltall(nneg:*),tmeltall(nneg:*),tarr(i));+Sshift

                  sice6melt = tmp6(1)
                  eice6melt = tmp6(2)

;v7.2 volume and mass fraction
                  vfraction = (darr(j)-dliqmelt) / (dice6melt-dliqmelt) ; of ice6 
                  dcheck = (1.-vfraction)*dliqmelt+vfraction*dice6melt
                  fraction = vfraction*dice6melt / $
                             (vfraction*dice6melt + (1.d0-vfraction)*dliqmelt)

                  parr(i,j) = pliqmelt
                  earr(i,j) = (1.-fraction)*eliqmelt+fraction*eice6melt
                  sarr(i,j) = (1.-fraction)*sliqmelt+fraction*sice6melt

                  phasearr(i,j) = 'L6'
;                  print,'L6: ',tarr(i),parr(i,j),darr(j),fraction,earr(i,j),sarr(i,j)
                  notdone=0
;                  print,darr(j),dcheck,dliqmelt,dice6melt
                 ; stop
              end

             ; ---------------- Ih - VI -------------------------------
;              if tarr(i) le ttp3 and darr(j) gt d1max and darr(j) le d61min and notdone then begin
;                  ; mix ice Ih-ice 6 below the ice III TP
;                  phasearr(i,j) = '16'
;                  notdone=0
;                                ;  stop
;              end

             ; ---------------- VI -------------------------------
              if tarr(i) le ttp3 then d6min = d61min ; ice VI on the I-VI boundary
              if tarr(i) gt ttp3 then d6min = dice6melt ; ice VI on the L=VI boundary
; old -- vertical boundary
              ;tmp6 = icevi(tarr(i),ptp6) ; K MPa
; v5.0: sloped boundary, v6.1
              tmp6 = icevi(tarr(i),b67(tarr(i))) ; K MPa
              d6max = tmp6(0)

; v6.0 -- try extending ice VI all the way to ice VII density
;              tmp6 = icevii(tarr(i),b67(tarr(i)))
;              d6max = tmp6(0)

              if darr(j) gt d6min and darr(j) le d6max and notdone then begin
                  ; ice VI
                  ; fit density/pressure
                  pfit = [ptp3,300.,400.,500.,600.,ptp6]
                  tmp6 = icevi(tarr(i),pfit)
                  dfit = tmp6(0:n_elements(pfit)-1)
                  pint = interpol(pfit,dfit,darr(j))

                  tmp6 = icevi(tarr(i),pint)
                  parr(i,j) = pint
                  sarr(i,j) = tmp6(1)
                  earr(i,j) = tmp6(2)
                  phasearr(i,j) = ' 6'
                  notdone=0
                  ;print,'6: ',tarr(i),pint,tmp6
                  ;stop
              end

             ; ---------------- VI-VII -------------------------------
              ;tmp7 = icevii(tarr(i),ptp6)
              tmp7 = icevii(tarr(i),b67(tarr(i)))
              d7min = tmp7(0)
              if darr(j) gt d6max and darr(j) le d7min and notdone then begin
                  ; ice VI-VII mix

                  s6 = tmp6(1)
                  e6 = tmp6(2)
                  s7 = tmp7(1)
                  e7 = tmp7(2)

; v7.2 volume and mass fraction
                  vfraction = (darr(j)-d6max) / (d7min-d6max) ; of ice 7
                  fraction = vfraction*d7min / $
                             (vfraction*d7min + (1.d0-vfraction)*d6max)

;                  parr(i,j) = ptp6 ; incorrect from v4.3
                  parr(i,j) = b67(tarr(i)); use sloped boundary
                  earr(i,j) = (1.-fraction)*e6+fraction*e7
                  sarr(i,j) = (1.-fraction)*s6+fraction*s7

                  ;print,'67: ',tarr(i),darr(j),parr(i,j),fraction,earr(i,j),sarr(i,j)
                  phasearr(i,j) = '67'
                  notdone=0
                  ;stop
              end

             ; ---------------- VII -------------------------------
              if darr(j) gt d7min and notdone then begin
                  ; ice VII
                  ; fit density/pressure
                  pfit = [ptp6,700.,800.,900.,1000.,2500.,5000.,7500.,10000.,25000.,50000.,75000.,100000.]
                  tmp7 = icevii(tarr(i),pfit)
                  dfit = tmp7(0:n_elements(pfit)-1)
                  pint = interpol(pfit,dfit,darr(j))

                  tmp7 = icevii(tarr(i),pint)
                  parr(i,j) = pint
                  sarr(i,j) = tmp7(1)
                  earr(i,j) = tmp7(2)
                  phasearr(i,j) = ' 7'
                  notdone=0
;                  print,'7: ',tarr(i),darr(j),tmp7
;                  stop
              end

          end ; T<TP

          ;-------------------------------------------------------------
          ; ------------- TP < T < TCRIT -------------------------------
          ;-------------------------------------------------------------
          if tarr(i) gt ttp and tarr(i) le tcrit and notdone then begin
              ; above the triple point and below the critical point
              dvapsat = interpol(dvsatall,tsatall,tarr(i))
              dliqsat = interpol(dlsatall,tsatall,tarr(i))
              dliqmelt = interpol(dmeltall,tmeltall,tarr(i))
              dice1melt = interpol(dice1meltall,tmeltall(0:nneg-1),tarr(i))
              ;dice6melt defined below
              ;dice7melt defined below

              ; ----------------------- VAPOR -----------------------------
              if darr(j) le dvapsat and notdone then begin
                                ; pure vapor - call EOSTD
                  cmd= './eostd '+strtrim(String(tarr(i)),2)+' '+strtrim(string(darr(j)),2)
                  spawn,cmd,res
                  tmp = strsplit(res,/extract)
                  if tmp(1) eq 0.0 and tmp(2) eq 0.0 then print,'line 762 EOSTD: ',tarr(i),darr(j) ; ERROR IN EOS
;                  if tmp(1) eq 0.0 and tmp(2) eq 0.0 then stop ; ERROR IN EOS
                  phasearr(i,j)= ' V'
                  parr(i,j)= tmp(0)
                  earr(i,j)= tmp(1)+Eshift
                  sarr(i,j)= tmp(2)+Sshift
                  ;print,'V: ',cmd, res
                  ;print,sarr(i,j),earr(i,j)
                  notdone=0
                  ;stop
               end

              ; ------------- VAPOR-LIQUID -------------------------------
              if darr(j) gt dvapsat and darr(j) le dliqsat and notdone then begin

; ADD LIQUID TENSION REGION HERE v6.2
; allow liquid to reach 600 kg/m3 and see what happens
; careful near the critical point!
                if tension then begin
                  if darr(j) gt 600.d0 then begin 
                  ; LIQUID IN TENSION v6.2
                  ; use Birch-Murnaghan equation to determine tension region
                      kt0 = 2.19d9 ; Pa at 20C based on sound speed and density
                      ktp = 4.76d0 ; from hugoniot ktp=4s-1 (probably not terribly accurate for the pressure range here)
                      v0 = 1./dliqsat ; mks
                      vneg = 1./darr(j) ; mks
                      xarr  = (vneg/v0)^(1./3.d0)
                      pneg = 3.d0*kt0*(1.-xarr)*xarr^(-2) * exp(3./2.*(ktp-1.)*(1.-xarr)) ; 3rd order B-M eq.

                  ; v8.0
                  ; only fill array if pneg < limit
                    if pneg/1.e6 gt liquid_tension_limit then begin
                        parr(i,j) = pneg/1.e6
                        earr(i,j) = -1.d0
                        sarr(i,j) = -1.d0
                        phasearr(i,j)= 'LT'
                        notdone=0
                    end ; else must be vapor-liquid region
                end ; end density for tension region
               end ; end flag for tension region

               if notdone then begin
                                ; must be mixed vapor-liquid
                  phasearr(i,j)= 'VL'
                  notdone=0
                  cmd= './eostd '+strtrim(String(tarr(i)),2)+' '+strtrim(string(darr(j)),2)
                  spawn,cmd,res
                  tmp = strsplit(res,/extract)
                  if tmp(1) eq 0.0 and tmp(2) eq 0.0 then print,'line 780 EOSTD: ',tarr(i),darr(j) ; ERROR IN EOS
;                  if tmp(1) eq 0.0 and tmp(2) eq 0.0 then stop ; ERROR IN EOS
                  parr(i,j)= tmp(0)
                  earr(i,j)= tmp(1)+Eshift ; these don't seem to be correct
                  sarr(i,j)= tmp(2)+Sshift ; these don't seem to be correct
                  ;print,cmd, res

                                ; weight by mass fraction using the
                                ; saturation vapor curves
                  dvv = interpol(dvsatall,tsatall,tarr(i))
                  dll = interpol(dlsatall,tsatall,tarr(i))
                  evv = interpol(evsatall,tsatall,tarr(i))
                  ell = interpol(elsatall,tsatall,tarr(i))
                  svv = interpol(svsatall,tsatall,tarr(i))
                  sll = interpol(slsatall,tsatall,tarr(i))

; v7.2 volume and mass fraction 
                  vfraction = (darr(j)-dvv) / (dll-dvv) ; of liquid
                  fraction = vfraction*dll / $
                             (vfraction*dll + (1.d0-vfraction)*dvv)

                  earr(i,j) = (1.-fraction)*evv+fraction*ell
                  sarr(i,j) = (1.-fraction)*svv+fraction*sll

                  ;print,dvv,dll,fraction,earr(i,j),sarr(i,j)
;                  stop
               end ; mixed vapor-liquid region
            end ; density range for mixed vapor-liquid region
              ; -------------- LIQUID ------------------------------------
              if darr(j) gt dliqsat and darr(j) le dliqmelt and notdone then begin
                  ; pure liquid
                  cmd= './eostd '+strtrim(String(tarr(i)),2)+' '+strtrim(string(darr(j)),2)
                  spawn,cmd,res
                  tmp = strsplit(res,/extract)
                  if tmp(1) eq 0.0 and tmp(2) eq 0.0 then print,'line 812 EOSTD: ',tarr(i),darr(j) ; ERROR IN EOS
;                  if tmp(1) eq 0.0 and tmp(2) eq 0.0 then stop ; ERROR IN EOS
                  phasearr(i,j)= ' L'
                  parr(i,j)= tmp(0)
                  earr(i,j)= tmp(1)+Eshift
                  sarr(i,j)= tmp(2)+Sshift
                  ;print,cmd, res
                  notdone=0
              end

              ; ------------- LIQUID - ICE VI -------------------------------
              p6melt = interpol(pmeltall(nneg:*),tmeltall(nneg:*),tarr(i))
              if tarr(i) le ttp6 then begin
                  tmp6 = icevi(tarr(i),p6melt) ; in K and MPa
                  dice6melt = tmp6(0) ; density of ice VI on the L-VI boundary at tarr(i)
              end else dice6melt = 1.e20 ; temp too high for ice 6
              if tarr(i) le ttp6 and darr(j) gt dliqmelt and darr(j) le dice6melt and notdone then begin
                  ; temp above tpt
                  pliqmelt = interpol(pmeltall(nneg:*),tmeltall(nneg:*),tarr(i))
                  eliqmelt = interpol(emeltall(nneg:*),tmeltall(nneg:*),tarr(i));+Eshift
                  sliqmelt = interpol(smeltall(nneg:*),tmeltall(nneg:*),tarr(i));+Sshift

                  sice6melt = tmp6(1)
                  eice6melt = tmp6(2)

; v7.2 volume and mass fraction
                  vfraction = (darr(j)-dliqmelt) / (dice6melt-dliqmelt) ; of ice6
                  fraction = vfraction*dice6melt / $
                             (vfraction*dice6melt + (1.d0-vfraction)*dliqmelt)

                  parr(i,j) = pliqmelt
                  earr(i,j) = (1.-fraction)*eliqmelt+fraction*eice6melt
                  sarr(i,j) = (1.-fraction)*sliqmelt+fraction*sice6melt

                  ; mixed liquid-solid
                  phasearr(i,j)= 'L6'
                  ;print,'L6 t>tpt: ',tarr(i),parr(i,j),darr(j),fraction,earr(i,j),sarr(i,j)
                  notdone=0
                  ;stop
              end

              ; ------------- ICE VI -------------------------------
              if tarr(i) le ttp6 then begin
                  ;tmp6 = icevi(tarr(i),ptp6) ; K MPa
; v5.0 - sloped boundary, v6.1
                  tmp6 = icevi(tarr(i),b67(tarr(i))) ; K MPa
                  d6max = tmp6(0)
; v6.0 -- try extending ice VI all the way to ice VII density
;                  tmp6 = icevii(tarr(i),b67(tarr(i)))
;                  d6max = tmp6(0)
              end else d6max = 1.d20 ; temp too high for ice 6

              if tarr(i) le ttp6 and darr(j) gt dice6melt and darr(j) le d6max and notdone then begin
                  ; pure ice VI
                  ; fit density/pressure
                  pfit = [ptp3,300.,400.,500.,600.,ptp6]
                  tmp6 = icevi(tarr(i),pfit)
                  dfit = tmp6(0:n_elements(pfit)-1)
                  pint = interpol(pfit,dfit,darr(j))

                  tmp6 = icevi(tarr(i),pint)
                  parr(i,j) = pint
                  sarr(i,j) = tmp6(1)
                  earr(i,j) = tmp6(2)
                  ;print,'6 t>tpt: ',tarr(i),pint,tmp6
                  phasearr(i,j)= ' 6'
                  notdone=0
              end

             ; ---------------- VI-VII -------------------------------
              ;tmp7 = icevii(tarr(i),ptp6)
              tmp7 = icevii(tarr(i),b67(tarr(i)))
              d7min = tmp7(0)
              if tarr(i) le ttp6 and darr(j) gt d6max and darr(j) le d7min and notdone then begin
                  ; ice VI-VII mix

                  s6 = tmp6(1)
                  e6 = tmp6(2)
                  s7 = tmp7(1)
                  e7 = tmp7(2)

; v7.2 volume and mass fraction
                  vfraction = (darr(j)-d6max) / (d7min-d6max) ; of ice 7
                  fraction = vfraction*d7min / $
                             (vfraction*d7min + (1.d0-vfraction)*d6max)

;                  parr(i,j) = ptp6 ; this is bug from v4.3
                  parr(i,j) = b67(tarr(i)) ; use the sloped boundary
                  earr(i,j) = (1.-fraction)*e6+fraction*e7
                  sarr(i,j) = (1.-fraction)*s6+fraction*s7

                  ;print,'67 t>tpt: ',tarr(i),darr(j),parr(i,j),fraction,earr(i,j),sarr(i,j)
                  phasearr(i,j) = '67'
                  notdone=0
                  ;stop
              end

              ; ------------- LIQUID - ICE VII -------------------------------
              p7melt = interpol(pmeltall(nneg:*),tmeltall(nneg:*),tarr(i))
              tmp7 = icevii(tarr(i),p7melt) ; in K and MPa
              dice7melt = tmp7(0) ; density of ice VII on the L-VII boundary at tarr(i)
              if tarr(i) gt ttp6 and darr(j) gt dliqmelt and darr(j) le dice7melt and notdone then begin
                  ; temp above tpt
                  ; temp above tpt
                  pliqmelt = interpol(pmeltall(nneg:*),tmeltall(nneg:*),tarr(i))
                  eliqmelt = interpol(emeltall(nneg:*),tmeltall(nneg:*),tarr(i));+Eshift
                  sliqmelt = interpol(smeltall(nneg:*),tmeltall(nneg:*),tarr(i));+Sshift

                  sice7melt = tmp7(1)
                  eice7melt = tmp7(2)

; v7.2 mass and volume fraction
                  vfraction = (darr(j)-dliqmelt) / (dice7melt-dliqmelt) ; of ice7
                  fraction = vfraction*dice7melt / $
                             (vfraction*dice7melt + (1.d0-vfraction)*dliqmelt)

                  parr(i,j) = pliqmelt
                  earr(i,j) = (1.-fraction)*eliqmelt+fraction*eice7melt
                  sarr(i,j) = (1.-fraction)*sliqmelt+fraction*sice7melt

                  ; mixed liquid-solid
                  phasearr(i,j)= 'L7'
                  ;print,'L7: ',tarr(i),parr(i,j),darr(j),fraction,earr(i,j),sarr(i,j)
                  notdone=0
                  ;stop
              end

              ; ------------- ICE VII -------------------------------
              ; first if above the L-VI-VII triple point
              if tarr(i) gt ttp6 and darr(j) gt dice7melt and notdone then begin
                  ; pure ice VII
                  ; fit density/pressure
                  pfit = [ptp6,700.,800.,900.,1000.,2500.,5000.,7500.,10000.,25000.,50000.,75000.,100000.]
                  tmp7 = icevii(tarr(i),pfit)
                  dfit = tmp7(0:n_elements(pfit)-1)
                  pint = interpol(pfit,dfit,darr(j))

                  tmp7 = icevii(tarr(i),pint)
                  parr(i,j) = pint
                  sarr(i,j) = tmp7(1)
                  earr(i,j) = tmp7(2)
                  phasearr(i,j) = ' 7'
                  ;print,'7: ',tarr(i),pint,tmp6
                  notdone=0
              end
              ; check if denser than ice VI
              if tarr(i) le ttp6 and darr(j) gt d7min and notdone then begin
                  ; pure ice VII
                  ; fit density/pressure
                  pfit = [ptp6,700.,800.,900.,1000.,2500.,5000.,7500.,10000.,25000.,50000.,75000.,100000.]
                  tmp7 = icevii(tarr(i),pfit)
                  dfit = tmp7(0:n_elements(pfit)-1)
                  pint = interpol(pfit,dfit,darr(j))

                  tmp7 = icevii(tarr(i),pint)
                  parr(i,j) = pint
                  sarr(i,j) = tmp7(1)
                  earr(i,j) = tmp7(2)
                  phasearr(i,j) = ' 7'
                  ;print,'7: ',tarr(i),pint,tmp6
                  notdone=0
              end

          end ; TP<T<TCRIT

          ;-------------------------------------------------------------
          ; ------------- T > TCRIT ------------------------------------
          ;-----CHECK EXTRAPOLATION----------------
          ;-------------------------------------------------------------
          if tarr(i) gt tcrit and notdone then begin
              ; above the critical point

              ; -------------- VAPOR/SUPERCRITICAL LIQUID --------------
              dliqmelt = interpol(dmeltall,tmeltall,tarr(i))
              if darr(j) le dliqmelt and notdone then begin
                                ; supercritical fluid
                  if tarr(i) lt 3000.d0 then begin ; v5.0
;                  if tarr(i) lt 5000.d0 then begin ; v5.0
;                  if tarr(i) lt 4000.d0 then begin ; v5.0
;v4.3                  if tarr(i) lt 1300.d0 then begin ; or (tarr(i) lt 5000. and darr(j) lt 950.d0) then begin
                                ; code only good to 5000 K for liquid
                                ; densities
                      ; lower temp limit for high densities
                      cmd= './eostd '+strtrim(String(tarr(i)),2)+' '+strtrim(string(darr(j)),2)
                      spawn,cmd,res
                      tmp = strsplit(res,/extract)
                      isflg=tmp(3) ; check for solid
                      i2ph=tmp(4) ; if 1 then superheated vapor

                      phasearr(i,j)= ' S'
                      parr(i,j)= tmp(0)
                      earr(i,j)= tmp(1)+Eshift
                      sarr(i,j)= tmp(2)+Sshift
  ;                    print,cmd, res

;                      if tmp(1) eq 0.0 and tmp(2) eq 0.0 then print,'line 1019 EOSTD: ',tarr(i),darr(j) ; ERROR IN EOS
;                      if tmp(1) eq 0.0 and tmp(2) eq 0.0 then stop ; ERROR IN EOS
                      if  tmp(0) eq 0.d0 or (not finite(tmp(2))) then begin
                          ; extrapolate - errors here in energy v4.3
;                          tmp = where(phasearr(0:i-1,j) eq ' 7') ; v4.3
                           ;tmp = where(phasearr(i,*) eq ' S') ; v5
;                          testarr(i,j)= 1
;                          tmp = where(phasearr(0:i-1,j) eq ' S') ; v5
; CHECK CHECK STSM
;                          if n_elements(tmp) eq 1 and tmp(0) ne -1 then tmp = findgen(i)
;                          if tmp(0) eq -1 then ind = findgen(i) else ind = tmp


                          if j ne 0 then begin
                            if i2ph eq 1 then zzz = where(phasearr(i,0:j-1) eq ' S') else stop ; v5 - if i2ph=1 then supercritical region
                            parr(i,j)= interpol(reform(parr(i,zzz)),darr(zzz),darr(j))
                            earr(i,j)= interpol(reform(earr(i,zzz)),darr(zzz),darr(j))
                            sarr(i,j)= interpol(reform(sarr(i,zzz)),darr(zzz),darr(j))
                            print,'LINE1036 CHECK S:',i,j,parr(i,j),earr(i,j),sarr(i,j)
                        end else begin
                            ; interpolate in temp
                            if i2ph eq 1 then zzz = where(phasearr(0:i,j) eq ' S') else stop ; v5 - if i2ph=1 then supercritical region
                            parr(i,j)= interpol(reform(parr(zzz,j)),tarr(zzz),tarr(i))
                            earr(i,j)= interpol(reform(earr(zzz,j)),tarr(zzz),tarr(i))
                            sarr(i,j)= interpol(reform(sarr(zzz,j)),tarr(zzz),tarr(i))
                            print,'LINE1045 TEXP CHECK S:',i,j,parr(i,j),earr(i,j),sarr(i,j)
;                          stop
                        end
                    end
;                      print,parr(i,j),earr(i,j),sarr(i,j)
                  end else begin
                      ; extrapolate to higher temps
                      phasearr(i,j)= 'ES'
;                      stop ; use 0:i-1
;                      tmp = where(phasearr(0:i-1,j) eq ' 7')
;                      if tmp(0) eq -1 then ind = findgen(i) else ind = tmp
;                      tmp = where(phasearr(0:i-1,j) eq ' S')
;                      if tmp(0) eq -1 then stop else ind = tmp
;                      parr(i,j)= interpol(reform(parr(ind,j)),tarr(ind),tarr(i))
;                      earr(i,j)= interpol(reform(earr(ind,j)),tarr(ind),tarr(i))
;                      sarr(i,j)= interpol(reform(sarr(ind,j)),tarr(ind),tarr(i))
;                      if parr(i,j) lt 0.d0 or earr(i,j) lt 0.d0 or sarr(i,j) lt 0.d0 then stop

;                      zzz = where(phasearr(i,0:j-1) eq ' S')
;                      if zzz(0) eq -1 then stop
;                      parr(i,j)= interpol(reform(parr(i,zzz)),darr(zzz),darr(i))
;                      earr(i,j)= interpol(reform(earr(i,zzz)),darr(zzz),darr(j))
;                      sarr(i,j)= interpol(reform(sarr(i,zzz)),darr(zzz),darr(j))

                      zzz = where(phasearr(0:i-1,j) eq ' S' or phasearr(0:i-1,j) eq 'ES')
                      if zzz(0) eq -1 then stop else ind = zzz 
                      if n_elements(zzz) eq 1 then begin
                           PRINT,'LINE1074: one element'
                           ind = [zzz(0)-1,zzz(0)]
                      end
                      parr(i,j)= interpol(reform(parr(ind,j)),tarr(ind),tarr(i))
                      earr(i,j)= interpol(reform(earr(ind,j)),tarr(ind),tarr(i))
                      sarr(i,j)= interpol(reform(sarr(ind,j)),tarr(ind),tarr(i))

                      ;stop
                  end


                  notdone=0
              end

              ; ------------ Supercritical Liquid - ICE VII -------------
              p7melt = interpol(pmeltall(nneg:*),tmeltall(nneg:*),tarr(i))
              tmp7 = icevii(tarr(i),p7melt) ; in K and MPa
              dice7melt = tmp7(0) ; density of ice VII on the L-VII boundary at tarr(i)
              if darr(j) gt dliqmelt and darr(j) le dice7melt and notdone then begin
                  ; temp above tpt
                  pliqmelt = interpol(pmeltall(nneg:*),tmeltall(nneg:*),tarr(i))
                  eliqmelt = interpol(emeltall(nneg:*),tmeltall(nneg:*),tarr(i));+Eshift
                  sliqmelt = interpol(smeltall(nneg:*),tmeltall(nneg:*),tarr(i));+Sshift

                  sice7melt = tmp7(1)
                  eice7melt = tmp7(2)

; v7.2 mass and volume fraction
                  vfraction = (darr(j)-dliqmelt) / (dice7melt-dliqmelt) ; of ice7
                  fraction  = vfraction*dice7melt / $
                              (vfraction*dice7melt + (1.d0-vfraction)*dliqmelt)

                  parr(i,j) = pliqmelt
                  earr(i,j) = (1.-fraction)*eliqmelt+fraction*eice7melt
                  sarr(i,j) = (1.-fraction)*sliqmelt+fraction*sice7melt

                  ; mixed supercritical fluid-solid
                  phasearr(i,j)= 'S7'
                  ;print,'S7: ',tarr(i),parr(i,j),darr(j),fraction,earr(i,j),sarr(i,j)
                  notdone=0
;                  stop
              end

              ; ------------ ICE VII -------------
              if darr(j) gt dice7melt and notdone then begin
                  ; pure ice VII
                  ; fit density/pressure
                  pfit = [ptp6,700.,800.,900.,1000.,2500.,5000.,7500.,10000.,25000.,50000.,75000.,100000.]
                  tmp7 = icevii(tarr(i),pfit)
                  dfit = tmp7(0:n_elements(pfit)-1)
                  pint = interpol(pfit,dfit,darr(j))

                  tmp7 = icevii(tarr(i),pint)
                  parr(i,j) = pint
                  sarr(i,j) = tmp7(1)
                  earr(i,j) = tmp7(2)
                  phasearr(i,j) = ' 7'
                  ;print,'7: ',tarr(i),pint,tmp6
                  notdone=0
              end
          end ; T>TCRIT

          ;-------------------------------------------------------------
          ; ------- CHECK IF FORGOT SOME T-D SPACE ---------------------
          ;-------------------------------------------------------------
          if notdone then begin
              print,'PHASE SPACE NOT CLASSIFIED!', tarr(i),darr(j)
              stop
          end

      end ; density loop
  end ; temperature loop

;-----------------------------------------------------------
; SAVE ORIGINAL PRESSURE ARRAY
porig = parr

;-----------------------------------------------------------
; GET RID OF NaN and Inf
zzz = where(finite(parr),complement=nnn)
if nnn(0) ne -1 then parr(nnn) =0.d0
zzz = where(finite(sarr),complement=nnn)
if nnn(0) ne -1 then sarr(nnn) =0.d0
zzz = where(finite(earr),complement=nnn)
if nnn(0) ne -1 then earr(nnn) =0.d0


;-----------------------------------------------------------
; TENSION REGION - OPTION B: fudged real calc
;                  not doing the full integral (P) dV
; pressure set above, calc energy and entropy along isotherm
; THIS IS NOT WORKING YET
if 0 then begin
for i=0,tsize-1 do begin
    zzz = where(parr(i,1:*) le +0.d0)
    if zzz(0) ne -1 and tarr(i) le ttp then begin
        zzz = zzz+1
        if tarr(i) le ttp3 then begin
            sref = interpol(sicesuball,tsuball,tarr(i))
            eref = interpol(eicesuball,tsuball,tarr(i))
            dref = interpol(dicesuball,tsuball,tarr(i))
        end else begin
            sref = interpol(sice1meltall,tmeltall(0:nneg-1),tarr(i))
            eref = interpol(eice1meltall,tmeltall(0:nneg-1),tarr(i))
            dref = interpol(dice1meltall,tmeltall(0:nneg-1),tarr(i))
        end
        tmpi = iheos(tarr(i),0)
        ds = tmpi(4)*tmpi(5)*(1./dijarr(0,zzz)-1./dref); J/kg/K
        sarr(i,zzz) = (sref + ds/1.e3) ; kJ/kg/K
        earr(i,zzz) = (eref*1.e3 + tarr(i)*ds - 1.e6*parr(i,zzz)*(1./dijarr(0,zzz)-1./dref))/1.e3 ; kJ/kg

        ind = findgen(max(zzz+10))
        plot,dijarr(i,ind),sarr(i,ind),psym=-7
        plot,dijarr(i,ind),earr(i,ind),psym=-7
        plot,dijarr(i,ind),parr(i,ind),psym=-7
        stop
    end
end
end

;-----------------------------------------------------------
; TENSION REGION - OPTION A: extrapolate
; pressure set above, extrapolate energy and entropy along isotherm

;if 0 then begin ; TEMP FOR TESTING v6.0
if tension then begin ; v6.1, v6.2, v6.3
 for i=0,tsize-1 do begin
    zzz = where(parr(i,1:*) le +0.d0)

    ; this is for ice Ih tension
    if zzz(0) ne -1 and tarr(i) le ttp then begin
        zzz = zzz+1
        if tarr(i) le ttp3 then $
                                ; note that the extension into 16 and
                                ; 1L is because near the temp minimum,
                                ; there may be very few points in the
                                ; ice Ih field
          yyy = where((phasearr(i,*) eq ' 1' or phasearr(i,*) eq '16') and parr(i,*) gt 0.d0) else $
          yyy = where((phasearr(i,*) eq ' 1' or phasearr(i,*) eq '1L') and parr(i,*) gt 0.d0)
        sarr(i,zzz) = interpol(sarr(i,yyy),dijarr(i,yyy),dijarr(i,zzz))
        earr(i,zzz) = interpol(earr(i,yyy),dijarr(i,yyy),dijarr(i,zzz))

        nnn = where(sarr(i,zzz) lt 0.d0)
        if nnn(0) ne -1 then sarr(i,zzz(nnn)) = min(sarr(i,yyy))
        nnn = where(earr(i,zzz) lt 0.d0)
        if nnn(0) ne -1 then earr(i,zzz(nnn)) = min(earr(i,yyy))

        if tarr(i) gt ttp3 then begin
          ; problem at 270K
            sarr(i,zzz) = min(sarr(i,yyy))
            earr(i,zzz) = min(sarr(i,yyy))
        end
        if 0 then begin
            ind = findgen(max(yyy))
            plot,dijarr(i,ind),sarr(i,ind),psym=-7
            plot,dijarr(i,ind),earr(i,ind),psym=-7
            plot,dijarr(i,ind),parr(i,ind),psym=-7
            if tarr(i) ge ttp3 then stop
        end
     end ; end of ice Ih tension

    ; this is for liquid water in tension v6.2
    if zzz(0) ne -1 and tarr(i) gt ttp then begin
        zzz = zzz+1
        yyy = where((phasearr(i,*) eq ' L') and parr(i,*) gt 0.d0) 
        sarr(i,zzz) = interpol(sarr(i,yyy),dijarr(i,yyy),dijarr(i,zzz))
        earr(i,zzz) = interpol(earr(i,yyy),dijarr(i,yyy),dijarr(i,zzz))

        nnn = where(sarr(i,zzz) lt 0.d0)
        if nnn(0) ne -1 then sarr(i,zzz(nnn)) = min(sarr(i,yyy))
        nnn = where(earr(i,zzz) lt 0.d0)
        if nnn(0) ne -1 then earr(i,zzz(nnn)) = min(earr(i,yyy))

        if 0 then begin
            ind = findgen(max(yyy))
            plot,dijarr(i,ind),sarr(i,ind),psym=-7
            plot,dijarr(i,ind),earr(i,ind),psym=-7
            plot,dijarr(i,ind),parr(i,ind),psym=-7
            if tarr(i) ge ttp then stop
         end
     end
                       ; end of LIQUID IN tension

 end ; loop over temps
end  ; use option A for tension


if 0 then begin
ind = findgen(20)
plot,dijarr(10,ind),sarr(10,ind),psym=-7
oplot,dijarr(5,ind),sarr(5,ind),psym=-7
plot,dijarr(10,ind),earr(10,ind),psym=-7
oplot,dijarr(5,ind),earr(5,ind),psym=-7
plot,dijarr(10,ind),parr(10,ind),psym=-7
oplot,dijarr(5,ind),parr(5,ind),psym=-7
end

;stop
;-----------------------------------------------------------
; ZERO DENSITY COLUMN: different if vapor or ice Ih in tension
; pressure set to zero, extrapolate energy and entropy along isotherm

if darr(0) eq 0.d0 then begin
    parr(*,0) = 0.d0
    phasearr(*,0) = 'ZZ' ; zero density column v6.4
    for i=0,tsize-1 do begin
        zzz = where(phasearr(i,*) eq phasearr(i,0))
        if n_elements(zzz) lt 3 then zzz = [zzz,zzz+1,zzz+2]
        zzz = zzz(1:*)
        sarr(i,0) = interpol(sarr(i,zzz),dijarr(i,zzz),0.d0)
        earr(i,0) = interpol(earr(i,zzz),dijarr(i,zzz),0.d0)
        if sarr(i,0) eq -1.d0 then stop ;sarr(i,0) = 0.d0
        if earr(i,0) eq -1.d0 then stop ;earr(i,0) = 0.d0
    end
end

;-----------------------------------------------------------
; v7.1
; ZERO TEMPERATURE ROW: 
; pressure -- set to zero to try to prevent extrapolation to negative
; pressures
; energy -- set to zero to try to prevent extrapolation to negative values
; entropy --- set to zero for now (probably should extrapolate along isotherms)
if 1 then begin ; add zero temp row
   ; update these variables
                                ; tsize, tarr, tijarr,
                                ; dijarr, parr, sarr, earr, phasearr
   ; (i,j) is (temp,density)
   tarr = [0.,tarr]
   tsize = tsize+1
   tmpt = tijarr
   tijarr = dblarr(tsize,dsize)
   tijarr(1:tsize-1,*) = tmpt 
   tijarr(0,*) = 0.d0
   tmpd = dijarr
   dijarr = dblarr(tsize,dsize)
   dijarr(1:tsize-1,*) = tmpd
   dijarr(0,*) = darr
   tmpph = phasearr
   phasearr = strarr(tsize,dsize)
   phasearr(1:tsize-1,*) = tmpph
   phasearr(0,*) = 'ZZ'
   tmpp = parr
   parr = dblarr(tsize,dsize)
   parr(1:tsize-1,*) = tmpp
   parr(0,*) = 0.d0
   tmpe = earr
   earr = dblarr(tsize,dsize)
   earr(1:tsize-1,*) = tmpe
   earr(0,*) = 0.d0
   tmps = sarr
   sarr = dblarr(tsize,dsize)
   sarr(1:tsize-1,*) = tmps
   sarr(0,*) = 0.d0
end                             ; add zero temp row
;-----------------------------------------------------------
; CHECK FOR very small pos/neg values - set to zero to prevent floating point errors
vmin = 1.e-20
fix=0
print,'SMALL VAL CHECK:'
zzz = where(sarr le vmin and sarr ge -vmin)
if zzz(0) ne -1 then print,'small val s:', dijarr(zzz),tijarr(zzz),parr(zzz),sarr(zzz)
if zzz(0) ne -1 then if fix then sarr(zzz) = 0.d0
zzz = where(earr le vmin  and earr ge -vmin)
if zzz(0) ne -1 then print,'small val e:', dijarr(zzz),tijarr(zzz),parr(zzz),earr(zzz)
if zzz(0) ne -1 then if fix then earr(zzz) = 0.d0
zzz = where(parr le vmin  and parr ge -vmin)
if zzz(0) ne -1 then print,'small val p:', dijarr(zzz),tijarr(zzz),parr(zzz)
if zzz(0) ne -1 then if fix then parr(zzz) = 0.d0

;print,tarr
;print,darr
;print,phasearr


if 0 then begin
;v4.3
; fix solid extrapolation
sorig = sarr
eorig = earr

for i=0,tsize-1 do begin
    stop; v4.3
    zzz = where(dijarr(i,*) gt 1000.d0 and dijarr(i,*) lt 6000. and phasearr(i,*) eq 'ES',complement=czzz)
    if zzz(0) ne -1 and czzz(0) ne -1 then begin ; v5
        sarr(i,zzz) = interpol(reform(sorig(i,czzz)),reform(dijarr(i,czzz)),reform(dijarr(i,zzz)))
        earr(i,zzz) = interpol(reform(eorig(i,czzz)),reform(dijarr(i,czzz)),reform(dijarr(i,zzz)))
    end

end
end

;-----------------------------------------------------------
; GET RID OF NaN and Inf ONE LAST TIME
zzz = where(finite(parr),complement=nnn)
if nnn(0) ne -1 then parr(nnn) =0.d0
zzz = where(finite(sarr),complement=nnn)
if nnn(0) ne -1 then sarr(nnn) =0.d0
zzz = where(finite(earr),complement=nnn)
if nnn(0) ne -1 then earr(nnn) =0.d0

if first then begin
  print, 'FIRST run: not saving sav file'
  print,'DONE MAKING FIRST FILES'
  print,'re-run to make table'
  stop
  stop
end else begin
  ; v7.1 get rid of sorig, eorig, porig
  save,filename=tablename+'.sav',tarr,darr,tijarr,dijarr,sarr,earr,parr,tsize,dsize,phasearr
;  save,filename=tablename+'.sav',tarr,darr,tijarr,dijarr,sarr,earr,parr,tsize,dsize,phasearr,sorig,eorig,porig
  print,'SAVING TABLE FILE: ',tablename,+'.sav'
end

if 0 then begin
  ; extra checks only
; check 0 density column
!p.multi=0
;plot,sarr,tijarr,psym=7,yr=[0,1000],xr=[0,10]
plot,sarr,tijarr,psym=7,yr=[0,10000],xr=[0,45],xtit='Entropy',ytit='Temp'
;oplot,sarr(*,1),tijarr(*,1),psym=4,col=200
zzz = where(parr lt 0.d0)
oplot,sarr(zzz),tijarr(zzz),psym=4,col=150
oplot,sarr(*,0),tijarr(*,0),psym=4,col=100

kkk = where(phasearr eq 'ES')
oplot,sarr(kkk),tijarr(kkk),psym=2,col=200
zzz2 = where(dijarr gt 1000.d0 and dijarr lt 6000.d0 and phasearr eq 'ES',complement=czzz)
oplot,sarr(zzz2),tijarr(zzz2),psym=2,col=150
stop

!p.multi=0
plot,earr,tijarr,psym=7,xtit='Energy',ytit='Temp';,yr=[0,1000],xr=[0,10]
;zzz = where(parr lt 0.d0)
;oplot,earr(zzz),tijarr(zzz),psym=4,col=150
;zzz = where(phasearr eq 'ES')
zzz = where(dijarr gt 1000.d0 and dijarr lt 6000. and phasearr eq 'ES',complement=czzz)
oplot,earr(zzz),tijarr(zzz),psym=4,col=200
stop
end

PLOTMESH:
restore,'meltcurve.sav'
restore,'subcurve.sav'
restore,'ice16curve.sav'
restore,'satcurve.sav'

; make colors to denote different phases
colarr = intarr(tsize,dsize)
zzz = where(phasearr eq ' 1')
if zzz(0) ne -1 then colarr(zzz) = 250
if zzz(0) ne -1 then print,' 1: ',n_elements(zzz) else print,' 1: 0'
zzz = where(phasearr eq '1T')
if zzz(0) ne -1 then colarr(zzz) = 240
if zzz(0) ne -1 then print,'1T: ',n_elements(zzz) else print,'1T: 0'
zzz = where(phasearr eq 'V1')
if zzz(0) ne -1 then colarr(zzz) = 10
if zzz(0) ne -1 then print,'VI: ',n_elements(zzz) else print,'VI: 0'
zzz = where(phasearr eq ' V')
if zzz(0) ne -1 then colarr(zzz) = 220
if zzz(0) ne -1 then print,'V: ',n_elements(zzz) else print,'V: 0'
zzz = where(phasearr eq 'VL')
if zzz(0) ne -1 then colarr(zzz) = 20
if zzz(0) ne -1 then print,'VL: ',n_elements(zzz) else print,'VL: 0'
zzz = where(phasearr eq ' L')
if zzz(0) ne -1 then colarr(zzz) = 200
if zzz(0) ne -1 then print,'L: ',n_elements(zzz) else print,'L: 0'
zzz = where(phasearr eq 'LT')
if zzz(0) ne -1 then colarr(zzz) = 190
if zzz(0) ne -1 then print,'LT: ',n_elements(zzz) else print,'LT: 0'
zzz = where(phasearr eq ' S')
if zzz(0) ne -1 then colarr(zzz) = 40
if zzz(0) ne -1 then print,'S: ',n_elements(zzz) else print,'S: 0'
zzz = where(phasearr eq 'ES')
if zzz(0) ne -1 then colarr(zzz) = 240
if zzz(0) ne -1 then print,'ES: ',n_elements(zzz) else print,'ES: 0'
zzz = where(phasearr eq '1L')
if zzz(0) ne -1 then colarr(zzz) = 60
if zzz(0) ne -1 then print,'1L: ',n_elements(zzz) else print,'1L: 0'
zzz = where(phasearr eq 'L6')
if zzz(0) ne -1 then colarr(zzz) = 5
if zzz(0) ne -1 then print,'L6: ',n_elements(zzz) else print,'L6: 0'
zzz = where(phasearr eq '16')
if zzz(0) ne -1 then colarr(zzz) = 245
if zzz(0) ne -1 then print,'16: ',n_elements(zzz) else print,'16: 0'
zzz = where(phasearr eq ' 6')
if zzz(0) ne -1 then colarr(zzz) = 150
if zzz(0) ne -1 then print,' 6: ',n_elements(zzz) else print,' 6: 0'
zzz = where(phasearr eq '67')
if zzz(0) ne -1 then colarr(zzz) = 100
if zzz(0) ne -1 then print,'67: ',n_elements(zzz) else print,'67: 0'
zzz = where(phasearr eq 'L7')
if zzz(0) ne -1 then colarr(zzz) = 140
if zzz(0) ne -1 then print,'L7: ',n_elements(zzz) else print,'L7: 0'
zzz = where(phasearr eq ' 7')
if zzz(0) ne -1 then colarr(zzz) = 120
if zzz(0) ne -1 then print,' 7: ',n_elements(zzz) else print,' 7: 0'
zzz = where(phasearr eq 'S7')
if zzz(0) ne -1 then colarr(zzz) = 100
if zzz(0) ne -1 then print,'S7: ',n_elements(zzz) else print,'S7: 0'
zzz = where(phasearr eq 'ZZ') ; v6.4
if zzz(0) ne -1 then colarr(zzz) = 1 ; near black
if zzz(0) ne -1 then print,'ZZ: ',n_elements(zzz) else print,'ZZ: 0'

zzz = where(colarr eq 0)
print,zzz

;stop
colorit=1
phaselab=0
printit=1
if colorit then loadct,37 else loadct, 0

if printit then begin
    lzrpc,tablename
    !p.multi=0
    !p.charsize=1.5
    nt=2
    !p.thick=nt
    !p.charthick=nt
    !x.thick=nt
    !y.thick=nt
end else begin
    window,0, retain=1
    device, decomposed=0
    !p.multi=[0,1,1]
    !p.charsize=1.5
    nt=1
    !p.thick=nt
    !p.charthick=nt
    !x.thick=nt
    !y.thick=nt
end

if not(printit) then window,/free
plot_io,darr,earr,/nodata,xtit='Temperature (K)',ytit='Energy (kJ/kg)',yr=[1.e-2,1.e6],xr=[min(tarr),300.]
oplot,tmeltall,emeltall
oplot,tsuball,esuball
oplot,tsuball,esuball
;oplot,dlsatall,esatall,thick=3
;oplot,dvsatall,esatall,thick=3
;if phaselab then for i=0,tsize-1 do for j=0,dsize-1 do $
;  xyouts,darr(j),earr(i,j),phasearr(i,j), alignment=0.5,charsize=1.
if colorit then  for i=0,tsize-1 do for j=0,dsize-1 do $
  plots,tarr(i),earr(i,j),psym=6,col=colarr(i,j),symsize=.2
if not(printit) then stop

; v7.1 for Laurel density-energy log plot for vapor boundary
if not(printit) then window,/free
plot_oo,darr,earr,/nodata,xtit='Density (kg/m3)',ytit='Energy (kJ/kg)',yr=[1.e-2,1.e6],xr=[min(darr),max(darr)]
oplot,dmeltall,emeltall
oplot,dsuball,esuball
oplot,dicesuball,esuball
;oplot,dlsatall,esatall,thick=3
;oplot,dvsatall,esatall,thick=3
;if phaselab then for i=0,tsize-1 do for j=0,dsize-1 do $
;  xyouts,darr(j),earr(i,j),phasearr(i,j), alignment=0.5,charsize=1.
if colorit then  for i=0,tsize-1 do for j=0,dsize-1 do $
  plots,darr(j),earr(i,j),psym=6,col=colarr(i,j),symsize=.2
if not(printit) then stop

if not(printit) then window,/free
plot_oo,darr,earr,/nodata,xtit='Density (kg/m3)',ytit='Energy (kJ/kg)',yr=[1000,1.e4],xr=[1.,1000.]
oplot,dmeltall,emeltall
oplot,dsuball,esuball
oplot,dicesuball,esuball
;oplot,dlsatall,esatall,thick=3
;oplot,dvsatall,esatall,thick=3
;if phaselab then for i=0,tsize-1 do for j=0,dsize-1 do $
;  xyouts,darr(j),earr(i,j),phasearr(i,j), alignment=0.5,charsize=1.
if colorit then  for i=0,tsize-1 do for j=0,dsize-1 do $
  plots,darr(j),earr(i,j),psym=6,col=colarr(i,j),symsize=.2
if not(printit) then stop

if not(printit) then window,/free
plot_io,darr,earr,/nodata,xtit='Density (kg/m3)',ytit='Energy (kJ/kg)',yr=[1,1.e4],xr=[100.,1000.]
oplot,dmeltall,emeltall
oplot,dsuball,esuball
oplot,dicesuball,esuball
;oplot,dlsatall,esatall,thick=3
;oplot,dvsatall,esatall,thick=3
;if phaselab then for i=0,tsize-1 do for j=0,dsize-1 do $
;  xyouts,darr(j),earr(i,j),phasearr(i,j), alignment=0.5,charsize=1.
if colorit then  for i=0,tsize-1 do for j=0,dsize-1 do $
  plots,darr(j),earr(i,j),psym=6,col=colarr(i,j),symsize=.2
if not(printit) then stop

; for robert density-entropy
if not(printit) then window,/free
plot,darr,sarr,/nodata,ytit='Entropy (kJ/kg/K)',xtit='Density (kg/m3)',/xsty,/ysty, xr=[0,5000.],yr=[0,20.],tit='Entropy-density'
oplot,dmeltall,smeltall
oplot,dsuball,ssuball
oplot,dlsatall,slsatall,thick=3
oplot,dvsatall,svsatall,thick=3
if phaselab then for i=0,tsize-1 do for j=0,dsize-1 do $
  xyouts,darr(j),sarr(i,j),phasearr(i,j), alignment=0.5,charsize=1.
if colorit then for i=0,tsize-1 do for j=0,dsize-1 do $
  plots,darr(j),sarr(i,j),psym=6,color=colarr(i,j),symsize=0.2
if not(printit) then stop

; for robert density-entropy zoom
if not(printit) then window,/free
plot,darr,sarr,/nodata,ytit='Entropy (kJ/kg/K)',xtit='Density (kg/m3)',/xsty,/ysty, xr=[0,2000.],yr=[0,15.],tit='Entropy-density zoom'
oplot,dmeltall,smeltall
oplot,dsuball,ssuball
oplot,dlsatall,slsatall,thick=3
oplot,dvsatall,svsatall,thick=3
if phaselab then for i=0,tsize-1 do for j=0,dsize-1 do $
  xyouts,darr(j),sarr(i,j),phasearr(i,j), alignment=0.5,charsize=1.
if colorit then for i=0,tsize-1 do for j=0,dsize-1 do $
  plots,darr(j),sarr(i,j),psym=6,color=colarr(i,j),symsize=0.2
if not(printit) then stop

; for robert density-entropy zoom no V1
if not(printit) then window,/free
plot,darr,sarr,/nodata,ytit='Entropy (kJ/kg/K)',xtit='Density (kg/m3)',/xsty,/ysty, xr=[0,2000.],yr=[0,15.],tit='With no VL and V1 mixed phases'
oplot,dmeltall,smeltall
oplot,dsuball,ssuball
oplot,dlsatall,slsatall,thick=3
oplot,dvsatall,svsatall,thick=3
if phaselab then for i=0,tsize-1 do for j=0,dsize-1 do $
  xyouts,darr(j),sarr(i,j),phasearr(i,j), alignment=0.5,charsize=1.
if colorit then for i=0,tsize-1 do for j=0,dsize-1 do $
  if phasearr(i,j) ne 'V1' then $
   if phasearr(i,j) ne 'VL' then $
    plots,darr(j),sarr(i,j),psym=6,color=colarr(i,j),symsize=0.2
if not(printit) then stop
;if printit then clzr else stop

if tension then begin
if not(printit) then window,/free
zzz = where(phasearr eq 'LT')
plot,darr,parr,/nodata,xtit='Density (kg/m3)',ytit='Pressure (MPa)',/xsty,/ysty,tit='Liquid Tension',xr=[min(dijarr(zzz)),max(dijarr(zzz))],yr=[min(parr(zzz)),100.]
;oplot,dmeltall,pmeltall
;oplot,dsuball,psuball
;oplot,dicesuball,psuball
;oplot,dlsatall,psatall,thick=3
;oplot,dvsatall,psatall,thick=3
oplot,[0,10000],[0,0]
if phaselab then for i=0,tsize-1 do for j=0,dsize-1 do $
  xyouts,darr(j),parr(i,j),phasearr(i,j), alignment=0.5,charsize=1.
if colorit then  for i=0,tsize-1 do for j=0,dsize-1 do $
  if phasearr(i,j) eq 'LT' then plots,darr(j),parr(i,j),psym=6,col=colarr(i,j),symsize=.2
if colorit then  for i=0,tsize-1 do for j=0,dsize-1 do $
  if phasearr(i,j) eq ' L' then plots,darr(j),parr(i,j),psym=6,col=colarr(i,j),symsize=.2
if not(printit) then stop

if not(printit) then window,/free
zzz = where(phasearr eq 'LT')
plot,earr,parr,/nodata,xtit='Energy (kJ/kg)',ytit='Pressure (MPa)',/xsty,/ysty,tit='Liquid Tension',xr=[min(earr(zzz)),max(earr(zzz))],yr=[min(parr(zzz)),100.]
;oplot,dmeltall,pmeltall
;oplot,dsuball,psuball
;oplot,dicesuball,psuball
;oplot,dlsatall,psatall,thick=3
;oplot,dvsatall,psatall,thick=3
oplot,[0,10000],[0,0]
if phaselab then for i=0,tsize-1 do for j=0,dsize-1 do $
  xyouts,earr(j),parr(i,j),phasearr(i,j), alignment=0.5,charsize=1.
if colorit then  for i=0,tsize-1 do for j=0,dsize-1 do $
  if phasearr(i,j) eq 'LT' then plots,earr(i,j),parr(i,j),psym=6,col=colarr(i,j),symsize=.2
if colorit then  for i=0,tsize-1 do for j=0,dsize-1 do $
  if phasearr(i,j) eq ' L' then plots,earr(i,j),parr(i,j),psym=6,col=colarr(i,j),symsize=.2
if not(printit) then stop

if not(printit) then window,/free
zzz = where(phasearr eq 'LT')
plot,sarr,parr,/nodata,xtit='Entropy (kJ/kg/K)',ytit='Pressure (MPa)',/xsty,/ysty,tit='Liquid Tension',xr=[min(sarr(zzz)),max(sarr(zzz))],yr=[min(parr(zzz)),100.]
;oplot,dmeltall,pmeltall
;oplot,dsuball,psuball
;oplot,dicesuball,psuball
;oplot,dlsatall,psatall,thick=3
;oplot,dvsatall,psatall,thick=3
oplot,[0,10000],[0,0]
if phaselab then for i=0,tsize-1 do for j=0,dsize-1 do $
  xyouts,earr(j),parr(i,j),phasearr(i,j), alignment=0.5,charsize=1.
if colorit then  for i=0,tsize-1 do for j=0,dsize-1 do $
  if phasearr(i,j) eq 'LT' then plots,sarr(i,j),parr(i,j),psym=6,col=colarr(i,j),symsize=.2
if colorit then  for i=0,tsize-1 do for j=0,dsize-1 do $
  if phasearr(i,j) eq ' L' then plots,sarr(i,j),parr(i,j),psym=6,col=colarr(i,j),symsize=.2
if not(printit) then stop


if not(printit) then window,/free
zzz = where(phasearr eq '1T')
plot,darr,parr,/nodata,xtit='Density (kg/m3)',ytit='Pressure (MPa)',/xsty,/ysty,tit='Ice Ih Tension',xr=[min(dijarr(zzz)),max(dijarr(zzz))],yr=[min(parr(zzz)),100.]
;oplot,dmeltall,pmeltall
;oplot,dsuball,psuball
;oplot,dicesuball,psuball
;oplot,dlsatall,psatall,thick=3
;oplot,dvsatall,psatall,thick=3
oplot,[0,10000],[0,0]
if phaselab then for i=0,tsize-1 do for j=0,dsize-1 do $
  xyouts,darr(j),parr(i,j),phasearr(i,j), alignment=0.5,charsize=1.
if colorit then  for i=0,tsize-1 do for j=0,dsize-1 do $
  if phasearr(i,j) eq '1T' then plots,darr(j),parr(i,j),psym=6,col=colarr(i,j),symsize=.2
if colorit then  for i=0,tsize-1 do for j=0,dsize-1 do $
  if phasearr(i,j) eq ' 1' then plots,darr(j),parr(i,j),psym=6,col=colarr(i,j),symsize=.2
if not(printit) then stop

if not(printit) then window,/free
zzz = where(phasearr eq '1T')
plot,earr,parr,/nodata,xtit='Energy (kJ/kg)',ytit='Pressure (MPa)',/xsty,/ysty,tit='Ice Ih Tension',xr=[min(earr(zzz)),max(earr(zzz))],yr=[min(parr(zzz)),100.]
;oplot,dmeltall,pmeltall
;oplot,dsuball,psuball
;oplot,dicesuball,psuball
;oplot,dlsatall,psatall,thick=3
;oplot,dvsatall,psatall,thick=3
oplot,[0,10000],[0,0]
if phaselab then for i=0,tsize-1 do for j=0,dsize-1 do $
  xyouts,earr(j),parr(i,j),phasearr(i,j), alignment=0.5,charsize=1.
if colorit then  for i=0,tsize-1 do for j=0,dsize-1 do $
  if phasearr(i,j) eq '1T' then plots,earr(i,j),parr(i,j),psym=6,col=colarr(i,j),symsize=.2
if colorit then  for i=0,tsize-1 do for j=0,dsize-1 do $
  if phasearr(i,j) eq ' 1' then plots,earr(i,j),parr(i,j),psym=6,col=colarr(i,j),symsize=.2
if not(printit) then stop

if not(printit) then window,/free
zzz = where(phasearr eq '1T')
plot,sarr,parr,/nodata,xtit='Entropy (kJ/kg/K',ytit='Pressure (MPa)',/xsty,/ysty,tit='Ice Ih Tension',xr=[min(sarr(zzz)),max(sarr(zzz))],yr=[min(parr(zzz)),100.]
;oplot,dmeltall,pmeltall
;oplot,dsuball,psuball
;oplot,dicesuball,psuball
;oplot,dlsatall,psatall,thick=3
;oplot,dvsatall,psatall,thick=3
oplot,[0,10000],[0,0]
if phaselab then for i=0,tsize-1 do for j=0,dsize-1 do $
  xyouts,sarr(j),parr(i,j),phasearr(i,j), alignment=0.5,charsize=1.
if colorit then  for i=0,tsize-1 do for j=0,dsize-1 do $
  if phasearr(i,j) eq '1T' then plots,sarr(i,j),parr(i,j),psym=6,col=colarr(i,j),symsize=.2
if colorit then  for i=0,tsize-1 do for j=0,dsize-1 do $
  if phasearr(i,j) eq ' 1' then plots,sarr(i,j),parr(i,j),psym=6,col=colarr(i,j),symsize=.2
if not(printit) then stop

end ; ;END OF TENSION PLOTS

if not(printit) then window,/free
plot,darr,tarr,/nodata,xtit='Entropy (kJ/kg/K)',ytit='Temperature (K)',xr=[0,4],yr=[0,500.]
oplot,smeltall,tmeltall
oplot,ssuball,tsuball
oplot,slsatall,tsatall,thick=3
oplot,svsatall,tsatall,thick=3
if phaselab then for i=0,tsize-1 do for j=0,dsize-1 do $
  xyouts,sarr(i,j),tarr(i),phasearr(i,j), alignment=0.5,charsize=1.
if colorit then for i=0,tsize-1 do for j=0,dsize-1 do $
  plots,sarr(i,j),tarr(i),psym=6,color=colarr(i,j),symsize=0.2
if not(printit) then stop


if not(printit) then window,/free
plot,darr,tarr,/nodata,xtit='Density (kg/m3)',ytit='Temperature (K)',yr=[0.,700],xr=[0,2000]
oplot,dmeltall,tmeltall
oplot,dsuball,tsuball
oplot,dicesuball,tsuball
oplot,dlsatall,tsatall,thick=3
oplot,dvsatall,tsatall,thick=3
oplot,d16all,t16all,thick=2
;oplot,b67(findgen(ttp6)),findgen(ttp6),thick=2
if phaselab then for i=0,tsize-1 do for j=0,dsize-1 do $
  xyouts,darr(j),tarr(i),phasearr(i,j), alignment=0.5,charsize=1.
if colorit then  for i=0,tsize-1 do for j=0,dsize-1 do $
  plots,darr(j),tarr(i),psym=6,col=colarr(i,j),symsize=.2

if not(printit) then stop
if not(printit) then window,/free

;plot,darr,tarr,/nodata,xtit='Density (kg/m3)',ytit='Temperature (K)',yr=[0,4000.],xr=[0,max(darr)]
plot,darr,tarr,/nodata,xtit='Density (kg/m3)',ytit='Temperature (K)',yr=[0,300.],xr=[0,max(darr)]
oplot,dmeltall,tmeltall
oplot,dsuball,tsuball
oplot,dicesuball,tsuball
oplot,dlsatall,tsatall,thick=3
oplot,dvsatall,tsatall,thick=3
if phaselab then for i=0,tsize-1 do for j=0,dsize-1 do $
  xyouts,darr(j),tarr(i),phasearr(i,j), alignment=0.5,charsize=1.
if colorit then  for i=0,tsize-1 do for j=0,dsize-1 do $
  plots,darr(j),tarr(i),psym=6,col=colarr(i,j),symsize=.2
if not(printit) then stop

if not(printit) then window,/free
plot_oo,darr,tarr,/nodata,xtit='Density (kg/m3)',ytit='Temperature (K)',yr=[min(tarr),max(tarr)],xr=[min(darr),max(darr)]
oplot,dmeltall,tmeltall
oplot,dsuball,tsuball
oplot,dicesuball,tsuball
oplot,dlsatall,tsatall,thick=3
oplot,dvsatall,tsatall,thick=3
if phaselab then for i=0,tsize-1 do for j=0,dsize-1 do $
  xyouts,darr(j),tarr(i),phasearr(i,j), alignment=0.5,charsize=1.
if colorit then  for i=0,tsize-1 do for j=0,dsize-1 do $
  plots,darr(j),tarr(i),psym=6,col=colarr(i,j),symsize=.2

if not(printit) then stop
if not(printit) then window,/free
plot_oi,darr,tarr,/nodata,xtit='Density (kg/m3)',ytit='Temperature (K)',yr=[0,1000.],xr=[min(darr),max(darr)]
oplot,dmeltall,tmeltall
oplot,dsuball,tsuball
oplot,dicesuball,tsuball
oplot,dlsatall,tsatall,thick=3
oplot,dvsatall,tsatall,thick=3
if phaselab then for i=0,tsize-1 do for j=0,dsize-1 do $
  xyouts,darr(j),tarr(i),phasearr(i,j), alignment=0.5,charsize=1.
if colorit then  for i=0,tsize-1 do for j=0,dsize-1 do $
  plots,darr(j),tarr(i),psym=6,col=colarr(i,j),symsize=.2
if not(printit) then stop

if not(printit) then window,/free
;plot_io,darr,parr,/nodata,xtit='Density (kg/m3)',ytit='Pressure (MPa)',yr=[10,20000.],xr=[0,2000]
plot_io,darr,parr,/nodata,xtit='Density (kg/m3)',ytit='Pressure (MPa)',yr=[1.e-6,1.e3],xr=[0,2000]
oplot,dmeltall,pmeltall
oplot,dsuball,psuball
oplot,dicesuball,psuball
oplot,dlsatall,psatall,thick=3
oplot,dvsatall,psatall,thick=3
if phaselab then for i=0,tsize-1 do for j=0,dsize-1 do $
  xyouts,darr(j),parr(i,j),phasearr(i,j), alignment=0.5,charsize=1.
if colorit then  for i=0,tsize-1 do for j=0,dsize-1 do $
  plots,darr(j),parr(i,j),psym=6,col=colarr(i,j),symsize=.2
if not(printit) then stop

if not(printit) then window,/free
; v6.4 - check ice Ih line near bottom of plot
plot_io,darr,earr,/nodata,xtit='Density (kg/m3)',ytit='Energy (kJ/kg)',yr=[1.e-2,1.e6],xr=[1,max(darr)]
oplot,dmeltall,emeltall
oplot,dsuball,esuball
oplot,dicesuball,esuball
;oplot,dlsatall,esatall,thick=3
;oplot,dvsatall,esatall,thick=3
if phaselab then for i=0,tsize-1 do for j=0,dsize-1 do $
  xyouts,darr(j),earr(i,j),phasearr(i,j), alignment=0.5,charsize=1.
if colorit then  for i=0,tsize-1 do for j=0,dsize-1 do $
  plots,darr(j),earr(i,j),psym=6,col=colarr(i,j),symsize=.2
if not(printit) then stop

if not(printit) then window,/free
plot_io,darr,earr,/nodata,xtit='Density (kg/m3)',ytit='Energy (kJ/kg)',yr=[1.e1,1.e4],xr=[1,2500]
oplot,dmeltall,emeltall
oplot,dsuball,esuball
oplot,dicesuball,esuball
;oplot,dlsatall,esatall,thick=3
;oplot,dvsatall,esatall,thick=3
if phaselab then for i=0,tsize-1 do for j=0,dsize-1 do $
  xyouts,darr(j),earr(i,j),phasearr(i,j), alignment=0.5,charsize=1.
if colorit then  for i=0,tsize-1 do for j=0,dsize-1 do $
  plots,darr(j),earr(i,j),psym=6,col=colarr(i,j),symsize=.2
if not(printit) then stop

if not(printit) then window,/free
;plot_oo,parr,tarr,/nodata,xtit='Pressure (MPa)',ytit='Temperature (K)',xr=[1.e-10,max(parr)],yr=[min(tarr),max(tarr)],/xsty,/ysty
plot_oo,parr,tarr,/nodata,xtit='Pressure (MPa)',ytit='Temperature (K)',xr=[min(parr),max(parr)],yr=[min(tarr),max(tarr)],/xsty,/ysty
oplot,pmeltall,tmeltall
oplot,picesuball,tsuball
oplot,psatall,tsatall,thick=3
if phaselab then for i=0,tsize-1 do for j=0,dsize-1 do $
  xyouts,parr(i,j),tarr(i),phasearr(i,j), alignment=0.5,charsize=1.
if colorit then  for i=0,tsize-1 do for j=0,dsize-1 do $
  plots,parr(i,j),tarr(i),psym=6,col=colarr(i,j),symsize=.2
if not(printit) then stop

if not(printit) then window,/free
plot,parr,tarr,/nodata,xtit='Pressure (MPa)',ytit='Temperature (K)',xr=[0,15.e3],yr=[0,800.]
oplot,pmeltall,tmeltall
oplot,picesuball,tsuball
oplot,psatall,tsatall,thick=3
if phaselab then for i=0,tsize-1 do for j=0,dsize-1 do $
  xyouts,parr(i,j),tarr(i),phasearr(i,j), alignment=0.5,charsize=1.
if colorit then  for i=0,tsize-1 do for j=0,dsize-1 do $
  plots,parr(i,j),tarr(i),psym=6,col=colarr(i,j),symsize=.2

if not(printit) then stop

if not(printit) then window,/free
plot_oo,parr,tarr,/nodata,xtit='Pressure (MPa)',ytit='Temperature (K)',xr=[1.e-13,max(parr)],yr=[min(tarr),max(tarr)],/xsty,/ysty
oplot,pmeltall,tmeltall
oplot,picesuball,tsuball
oplot,psatall,tsatall,thick=3
if phaselab then for i=0,tsize-1 do for j=0,dsize-1 do $
  xyouts,parr(i,j),tarr(i),phasearr(i,j), alignment=0.5,charsize=1.
if colorit then  for i=0,tsize-1 do for j=0,dsize-1 do $
  plots,parr(i,j),tarr(i),psym=6,col=colarr(i,j),symsize=.2
if not(printit) then stop

if not(printit) then window,/free
plot_oi,parr,tarr,/nodata,xtit='Pressure (MPa)',ytit='Temperature (K)',xr=[1.e-13,20.e3],yr=[0,1000.],/xsty,/ysty
oplot,pmeltall,tmeltall
oplot,picesuball,tsuball
oplot,psatall,tsatall,thick=3
if phaselab then for i=0,tsize-1 do for j=0,dsize-1 do $
  xyouts,parr(i,j),tarr(i),phasearr(i,j), alignment=0.5,charsize=1.
if colorit then  for i=0,tsize-1 do for j=0,dsize-1 do $
  plots,parr(i,j),tarr(i),psym=6,col=colarr(i,j),symsize=.2
if not(printit) then stop

if not(printit) then window,/free
plot,sarr,parr,/nodata,xtit='Entropy (kJ/kg/K)',ytit='Pressure (MPa)',xr=[0,max(sarr)],yr=[0,max(parr)],/xsty,/ysty
oplot,smeltall,pmeltall
oplot,ssuball,picesuball
oplot,slsatall,psatall,thick=3
oplot,svsatall,psatall,thick=3
if phaselab then for i=0,tsize-1 do for j=0,dsize-1 do $
  xyouts,sarr(i,j),parr(i,j),phasearr(i,j), alignment=0.5,charsize=1.
if colorit then for i=0,tsize-1 do for j=0,dsize-1 do $
  plots,sarr(i,j),parr(i,j),psym=6,color=colarr(i,j),symsize=0.2
if not(printit) then stop

if not(printit) then window,/free
plot,sarr,parr,/nodata,xtit='Entropy (kJ/kg/K)',ytit='Pressure (MPa)',xr=[0,10.],yr=[0,15.e3]
oplot,smeltall,pmeltall
oplot,ssuball,picesuball
oplot,slsatall,psatall,thick=3
oplot,svsatall,psatall,thick=3
if phaselab then for i=0,tsize-1 do for j=0,dsize-1 do $
  xyouts,sarr(i,j),parr(i,j),phasearr(i,j), alignment=0.5,charsize=1.
if colorit then for i=0,tsize-1 do for j=0,dsize-1 do $
  plots,sarr(i,j),parr(i,j),psym=6,color=colarr(i,j),symsize=0.2
if not(printit) then stop

if not(printit) then window,/free
plot_io,sarr,parr,/nodata,xtit='Entropy (kJ/kg/K)',ytit='Pressure (MPa)',xr=[0,max(sarr)],yr=[1.e-13,max(parr)],/xsty,/ysty
oplot,smeltall,pmeltall
oplot,ssuball,picesuball
oplot,slsatall,psatall,thick=3
oplot,svsatall,psatall,thick=3
if phaselab then for i=0,tsize-1 do for j=0,dsize-1 do $
  xyouts,sarr(i,j),parr(i,j),phasearr(i,j), alignment=0.5,charsize=1.
if colorit then for i=0,tsize-1 do for j=0,dsize-1 do $
  plots,sarr(i,j),parr(i,j),psym=6,color=colarr(i,j),symsize=0.2
if not(printit) then stop

if not(printit) then window,/free
plot_io,sarr,parr,/nodata,xtit='Entropy (kJ/kg/K)',ytit='Pressure (MPa)',xr=[0,20],yr=[1.e-13,25.e5],/xsty,/ysty
oplot,smeltall,pmeltall
oplot,ssuball,picesuball
oplot,slsatall,psatall,thick=3
oplot,svsatall,psatall,thick=3
if phaselab then for i=0,tsize-1 do for j=0,dsize-1 do $
  xyouts,sarr(i,j),parr(i,j),phasearr(i,j), alignment=0.5,charsize=1.
if colorit then for i=0,tsize-1 do for j=0,dsize-1 do $
  plots,sarr(i,j),parr(i,j),psym=6,color=colarr(i,j),symsize=0.2
if not(printit) then stop

if not(printit) then window,/free
plot,parr,earr,/nodata,xtit='Pressure (MPa)',ytit='Energy (kJ/kg)',xr=[0,max(parr)],yr=[0,max(earr)],/xsty,/ysty
oplot,pmeltall,emeltall
oplot,picesuball,esuball
oplot,psatall,elsatall,thick=3
oplot,psatall,evsatall,thick=3
if phaselab then for i=0,tsize-1 do for j=0,dsize-1 do $
  xyouts,parr(i,j),earr(i,j),phasearr(i,j), alignment=0.5,charsize=1.
if colorit then for i=0,tsize-1 do for j=0,dsize-1 do $
  plots,parr(i,j),earr(i,j),psym=6,color=colarr(i,j),symsize=0.2
if not(printit) then stop

if not(printit) then window,/free
plot,parr,earr,/nodata,xtit='Pressure (MPa)',ytit='Energy (kJ/kg)',xr=[0,15.e3],yr=[0,3000.]
oplot,pmeltall,emeltall
oplot,picesuball,esuball
oplot,psatall,elsatall,thick=3
oplot,psatall,evsatall,thick=3
if phaselab then for i=0,tsize-1 do for j=0,dsize-1 do $
  xyouts,parr(i,j),earr(i,j),phasearr(i,j), alignment=0.5,charsize=1.
if colorit then for i=0,tsize-1 do for j=0,dsize-1 do $
  plots,parr(i,j),earr(i,j),psym=6,color=colarr(i,j),symsize=0.2
if not(printit) then stop

if not(printit) then window,/free

plot_oi,parr,earr,/nodata,xtit='Pressure (MPa)',ytit='Energy (kJ/kg)',xr=[1.e-13,max(parr)],yr=[0,max(earr)],/xsty,/ysty
oplot,pmeltall,emeltall
oplot,picesuball,esuball
oplot,psatall,elsatall,thick=3
oplot,psatall,evsatall,thick=3
if phaselab then for i=0,tsize-1 do for j=0,dsize-1 do $
  xyouts,parr(i,j),earr(i,j),phasearr(i,j), alignment=0.5,charsize=1.
if colorit then for i=0,tsize-1 do for j=0,dsize-1 do $
  plots,parr(i,j),earr(i,j),psym=6,color=colarr(i,j),symsize=0.2
if not(printit) then stop

if not(printit) then window,/free

plot_oi,parr,earr,/nodata,xtit='Pressure (MPa)',ytit='Energy (kJ/kg)',xr=[1.e-13,10.e3],yr=[0,4000.]
oplot,pmeltall,emeltall
oplot,picesuball,esuball
oplot,psatall,elsatall,thick=3
oplot,psatall,evsatall,thick=3
if phaselab then for i=0,tsize-1 do for j=0,dsize-1 do $
  xyouts,parr(i,j),earr(i,j),phasearr(i,j), alignment=0.5,charsize=1.
if colorit then for i=0,tsize-1 do for j=0,dsize-1 do $
  plots,parr(i,j),earr(i,j),psym=6,color=colarr(i,j),symsize=0.2
if not(printit) then stop

plot_io,earr,parr,/nodata,ytit='Pressure (MPa)',xtit='Energy (kJ/kg)',yr=[1.e-13,10.e3],xr=[0,4000.]
oplot,emeltall,pmeltall
oplot,esuball,picesuball
oplot,elsatall,psatall,thick=3
oplot,evsatall,psatall,thick=3
if phaselab then for i=0,tsize-1 do for j=0,dsize-1 do $
  xyouts,earr(i,j),parr(i,j),phasearr(i,j), alignment=0.5,charsize=1.
if colorit then for i=0,tsize-1 do for j=0,dsize-1 do $
  plots,earr(i,j),parr(i,j),psym=6,color=colarr(i,j),symsize=0.2
if not(printit) then stop

if not(printit) then window,/free
plot_io,sarr,tarr,/nodata,xtit='Entropy (kJ/kg/K)',ytit='Temperature (K)',xr=[0,max(sarr)],yr=[min(tarr),max(tarr)],/xsty,/ysty
oplot,smeltall,tmeltall
oplot,ssuball,tsuball
oplot,slsatall,tsatall,thick=3
oplot,svsatall,tsatall,thick=3
if phaselab then for i=0,tsize-1 do for j=0,dsize-1 do $
  xyouts,sarr(i,j),tarr(i),phasearr(i,j), alignment=0.5,charsize=1.
if colorit then for i=0,tsize-1 do for j=0,dsize-1 do $
  plots,sarr(i,j),tarr(i),psym=6,color=colarr(i,j),symsize=0.2
if not(printit) then stop

if not(printit) then window,/free
plot,darr,tarr,/nodata,xtit='Entropy (kJ/kg/K)',ytit='Temperature (K)',xr=[0,20],yr=[0,800.]
oplot,smeltall,tmeltall
oplot,ssuball,tsuball
oplot,slsatall,tsatall,thick=3
oplot,svsatall,tsatall,thick=3
if phaselab then for i=0,tsize-1 do for j=0,dsize-1 do $
  xyouts,sarr(i,j),tarr(i),phasearr(i,j), alignment=0.5,charsize=1.
if colorit then for i=0,tsize-1 do for j=0,dsize-1 do $
  plots,sarr(i,j),tarr(i),psym=6,color=colarr(i,j),symsize=0.2
if not(printit) then stop

if not(printit) then window,/free
plot,darr,tarr,/nodata,xtit='Entropy (kJ/kg/K)',ytit='Temperature (K)',xr=[0,4],yr=[0,500.]
oplot,smeltall,tmeltall
oplot,ssuball,tsuball
oplot,slsatall,tsatall,thick=3
oplot,svsatall,tsatall,thick=3
if phaselab then for i=0,tsize-1 do for j=0,dsize-1 do $
  xyouts,sarr(i,j),tarr(i),phasearr(i,j), alignment=0.5,charsize=1.
if colorit then for i=0,tsize-1 do for j=0,dsize-1 do $
  plots,sarr(i,j),tarr(i),psym=6,color=colarr(i,j),symsize=0.2
if not(printit) then stop

if not(printit) then window,/free
plot,darr,tarr,/nodata,xtit='Entropy (kJ/kg/K)',ytit='Temperature (K)',xr=[0,4],yr=[0,500.]
oplot,smeltall,tmeltall
oplot,ssuball,tsuball
oplot,slsatall,tsatall,thick=3
oplot,svsatall,tsatall,thick=3
if phaselab then for i=0,tsize-1 do for j=0,dsize-1 do $
  xyouts,sarr(i,j),tarr(i),phasearr(i,j), alignment=0.5,charsize=1.
if colorit then for i=0,tsize-1 do for j=0,dsize-1 do $
  if phasearr(i,j) ne ' 7' then $
  if phasearr(i,j) ne '67' then $
  plots,sarr(i,j),tarr(i),psym=6,color=colarr(i,j),symsize=0.2
if not(printit) then stop

if printit then clzr
print,'done plotting mesh'
print,'next: save ASCII SESAME FILE'
;stop

WRITEMESH:
;restore save file above
; UNITS: t K, rho kg/m3, P MPa, E kJ/kg, S kJ/kg/K, phase arr string

if 0 then begin
    ; if didn't get rid of small values originally, do it now
    vmin = 1.d-20
    zzz = where(parr gt -vmin and parr lt vmin)
    if zzz(0) ne -1 then parr(zzz) = 0.d0
    zzz = where(earr gt -vmin and earr lt vmin)
    if zzz(0) ne -1 then earr(zzz) = 0.d0
    zzz = where(sarr gt -vmin and sarr lt vmin)
    if zzz(0) ne -1 then sarr(zzz) = 0.d0

    zzz = where(earr lt 0.d0)
    if zzz(0) ne -1 then earr(zzz) = 0.d0
end

if 1 then begin
;-----------------------------------------------------------
; GET RID OF NaN and Inf ONE LAST TIME
    zzz = where(finite(parr),complement=nnn)
    if nnn(0) ne -1 then parr(nnn) =0.d0
    zzz = where(finite(sarr),complement=nnn)
    if nnn(0) ne -1 then sarr(nnn) =0.d0
    zzz = where(finite(earr),complement=nnn)
    if nnn(0) ne -1 then earr(nnn) =0.d0
end

MATID = 3 ; CHANGE IF NECESSARY FOR CTH FILE
nwds = 2+tsize+dsize+tsize*dsize*3
date = 11282006
version = 1
ntables = 2
t1 = 201
t2 = 301
nwds1 = 5
nwds0 = 9

; water parameters, reference point STP
; might want to change for Mars....
fmn = 10.d0
fmw = 18.d0
T0_ref = 298.d0 ; K
; pick P or rho to define ref state
P0_ref = 1.e5/1.e6 ; 1 bar in MPa
tmpw = nisteos(t0_ref,P0_ref) ; K, MPa
rho0_ref = tmpw(0)
; cs^2= dP/drho
; take isotherm through T0 and P0
tt1 = nisteos(t0_ref,P0_ref*0.99) ; K, MPa
tt2 = nisteos(t0_ref,P0_ref*1.01) ; K, MPa
dd1 = tt1(0)
dd2 = tt2(0)
cs = sqrt((p0_ref*1.01-p0_ref*0.99)*1.e6/(dd2-dd1)) ; m/s
C0_ref = cs
K0_ref = (rho0_ref*cs^2.)/1.d6 ; MPa
print,'reference point rho, K0, T0: ',rho0_ref,k0_ref,t0_ref

; UNITS
; up until this point everything in NIST units
; K, MPa, kg/m3, kJ/kg

; now convert to BCAT/SESAME units:
; K, GPa, g/cm3, MJ/kg, MJ/kg/K
; tarr is the same
darr = darr/1.d3 ; g/cm3
parr = parr/1.d3 ; GPa
earr = earr/1.d3 ; MJ/kg
sarr = sarr/1.d3 ; MJ/kg/K
rho0_ref = rho0_ref/1.d3 ; g/cm3
k0_ref = k0_ref / 1.d3 ; GPa
; t0_ref is the same

openw,unit,tablename+'.txt',/get_lun
;openw,unit2,'check '+tablename+'.txt',/get_lun
; WRITE OUT INDEX SECTION
printf,unit,matid,nwds0,format="(' INDEX',5X,'MATID =',I7,5X,'NWDS = ',I5)"
printf,unit,matid,date,date,version,ntables,format="(5('  ',E14.8E2))"
printf,unit,t1,t2,nwds1,nwds,format="(4('  ',E14.8E2))"

; WRITE OUT 201 TABLE SECTION
printf,unit,t1,nwds1,format="(' RECORD',5X,'TYPE =',I5,5X,'NWDS = ',I5)"
printf,unit,fmn,fmw,rho0_ref,k0_ref,t0_ref,format="(5('  ',E14.8E2))"

; WRITE OUT 301 TABLE SECTION
printf,unit,t2,nwds,format="(' RECORD',5X,'TYPE =',I5,5X,'NWDS = ',I6)"

; first: # density points, # temperature points in mesh
printf,unit,dsize,format="('  ',E14.8E2,$)"
printf,unit,tsize,format="('  ',E14.8E2,$)"

; second: density array, temperature array
nw = 3.d0
for j=0,dsize-1 do begin
    if nw mod 5. eq 0. then printf,unit,darr(j),format="('  ',E14.8E2)" else $
      printf,unit,darr(j),format="('  ',E14.8E2,$)"
    nw=nw+1.
end ; d loop
for j=0,tsize-1 do begin
    if nw mod 5. eq 0. then printf,unit,tarr(j),format="('  ',E14.8E2)" else $
      printf,unit,tarr(j),format="('  ',E14.8E2,$)"
    nw=nw+1.
end ; d loop
; third: pressure, loop over temperture, loop over density
for i=0,tsize-1 do begin
    for j=0,dsize-1 do begin
        if nw mod 5. eq 0. then begin
            if parr(i,j) ge 0.d0 then printf,unit,parr(i,j),format="('  ',E14.8E2)"
            if parr(i,j) lt 0.d0 then printf,unit,parr(i,j),format="(' ',E15.8E2)"
        end else begin
            if parr(i,j) ge 0.d0 then printf,unit,parr(i,j),format="('  ',E14.8E2,$)"
            if parr(i,j) lt 0.d0 then printf,unit,parr(i,j),format="(' ',E15.8E2,$)"
        end
        nw=nw+1.
    end ; d loop
end ; t loop
; third: energy, loop over temperture, loop over density
for i=0,tsize-1 do begin
    for j=0,dsize-1 do begin
        if nw mod 5. eq 0. then begin
            print,nw,nw mod 5.,earr(i,j)
            if earr(i,j) ge 0.d0 then printf,unit,earr(i,j),format="('  ',E14.8E2)"
            if earr(i,j) lt 0.d0 then printf,unit,earr(i,j),format="(' ',E15.8E2)"
;            printf,unit2,'e ',i,j,earr(i,j)
        end else begin
            if earr(i,j) ge 0.d0 then printf,unit,earr(i,j),format="('  ',E14.8E2,$)"
            if earr(i,j) lt 0.d0 then printf,unit,earr(i,j),format="(' ',E15.8E2,$)"
        end
        nw=nw+1.
    end ; d loop
end ; t loop
; fourth: entropy, loop over temperture, loop over density
for i=0,tsize-1 do begin
    for j=0,dsize-1 do begin
        if nw mod 5. eq 0. then begin
            if sarr(i,j) ge 0.d0 then printf,unit,sarr(i,j),format="('  ',E14.8E2)"
            if sarr(i,j) lt 0.d0 then printf,unit,sarr(i,j),format="(' ',E15.8E2)"
        end else begin
            if sarr(i,j) ge 0.d0 then printf,unit,sarr(i,j),format="('  ',E14.8E2,$)"
            if sarr(i,j) lt 0.d0 then printf,unit,sarr(i,j),format="(' ',E15.8E2,$)"
        end
        nw=nw+1.
    end                         ; d loop
end ; t loop
printf,unit,' END'
close,unit
free_lun,unit

print,nw,nwds
print,'finished writing CTH SESAME mesh'
print,'CHECK END OF FILE FOR CORRECT LAST LINE'
print,'END OF MAKE_H2O_SES'
stop
stop
end



;=========================================================================================
;=========================================================================================
; ICE VI EOS
;
; based on Stewart and Ahrens 2005 with slight modifications to
; thermal expansion and using NIST water EOS as reference boundary
; points for S and E
;
; inputs in K and MPa
; t is a scalar
; pmpa can be an array
function icevi, t, pmpa
common setup, first
  ;-------------------------------------------------------------------------------------
  ; FIRST CALC V(T,P)
  ; reference isotherm T=230. K
  tref = 230.d0                   ; K
  v0=1./1300.d0                  ; m3/kg
  Kt0=13.08d9                     ; Pa
  Ktp=6.8d0

  ; estimates - not very sensitive
  a0 = 3.e-5 ; 1/K -- based on limited measurements Mishima et al JCP 1979
  a1 = 1.e-7 ; 1/K^2  - guess
  eta = 1.d  ; not sensitive to exponent
  if first then begin
      ta0arr = [findgen(400)+2.]
      a0max = 4.e-5
      Tdeb = 250.d0             ; K
      xarr = Tdeb/ta0arr
      dfint = (xarr^4.d0)*exp(xarr)/(exp(xarr)-1.d0)^2.d0
      darr = fltarr(n_elements(xarr))
      for i=0,n_elements(xarr)-2 do begin
                                ; normalized 3rd order Debye function
          darr(i) = 3.*(xarr(i)^(-3.))*int_tabulated(xarr(i:*),dfint(i:*),/sort,/double)
      end
      a0arr = darr*a0max        ; assume cp=cv for ice VII
      ; make cv stay constant to very high pressrues
      zzz = where(darr eq max(darr))
      a0arr(zzz:*) = max(darr)*a0max

;      plot,ta0arr,a0+a1*ta0arr
;      oplot,ta0arr,a0arr,lines=2
      save,filename='alpha6.sav',ta0arr,a0arr
;      stop
  end else restore,'alpha6.sav'

  ; range of real densities about 1.3 to 1.5 g/cm3, 0.66 to 0.77 cm3/g
  vref = v0-(findgen(101)/100.)*0.11*1.d-3
  pref = (kt0/ktp)*( (v0/vref)^ktp - 1.d0)
  ;plot,vref*1.e3,pref/1.e9

  ; calc V at a given P, T
  vi = dblarr(n_elements(pmpa))
  for i =0, n_elements(pmpa)-1 do begin
      ; loop over possible array of P values
      p = pmpa(i)*1.d6             ; Pa
      vpt0 = interpol(vref,pref,p)
      if t ne tref then begin
          tarrtmp = (findgen(101)/100.d0)*(t-tref)+tref
          apt0 = interpol(a0arr,ta0arr,tarrtmp)
          apt = apt0 * (1.d0 + Ktp*p/Kt0)^(-eta)
;          apt = (a0+a1*tarrtmp) * (1.d0 + Ktp*p/Kt0)^(-eta)
          tmp = int_tabulated(tarrtmp,apt,/double,/sort)
          if t lt tref then vpt = vpt0*exp(-tmp) else vpt = vpt0*exp(tmp)
      end else begin
          vpt = vpt0
      end
      vi(i) = vpt
  end

  ;-------------------------------------------------------------------------------------
                                ; next calc S - S is almost constant
                                ;               along isotherms
                                ;               because of the small
                                ;               pressure region

  ; latent heat across VI-liquid - from Dorsey book
  tlt = [52.5,57.2,66.0,73.8,80.8,81.6]+273.15 ; K
  llt = [333.5,336.5,339.1,345.3,352.0,354.5]*1000. ; J/kg
  ;pdvlt =[74.7,74.9,74.7,73.8,72.3,72.4]*1000. ; J/kg
  plt = [14518.,15485.,17421.,19357.,21292.,21680.]*1.0143*1.e5 ; Pa

  ; extend to I-III-L triple point at 251.16 K
  tlt = [251.16,273.5,tlt] ; K
  llt = [333.43d3,333.43d3,llt] ; J/kg use Ice Ih latent heat in extrapolation
  plt = [207.5d6,625.d6,plt] ; Pa

  v6lt = dblarr(n_elements(plt))
  vllt = dblarr(n_elements(plt))
  ;for i=0,n_elements(plt)-1 do v6lt(i) = iceVI(tlt(i),plt(i)/1.e6)
  ;restore,'meltcurve.sav'
  ;for i=0,n_elements(plt)-1 do vllt(i) = 1./interpol(dmeltall,pmeltall,plt(i)/1.e6)
  ;IDL> print,v6lt
  v6lt = [  0.00075820870 ,  0.00073885498,   0.00070880721 ,  0.00070587115,   0.00070025577, $
            0.00069496090 ,   0.00068995640,    0.00068897919]
  ;IDL> print,vllt
  vllt = [0.00091676372,   0.00084205879,   0.00077717935,   0.00077172491,   0.00076156368, $
          0.00075226851,   0.00074370507,   0.00074205618]
  dvlt = vllt-v6lt

  slt1 = llt/tlt
  slt2 = deriv(tlt,plt)*dvlt
  ;print,slt1,slt2
  ;stop
                                ; use slt1 --> makes llt go down in
                                ; the extrapolated region (ices II and V

  ; if 251.16 < t < 354.76 K then use the S along the phase boundary
  ; if t<251.16 then calc dS using cp at the same pressure
                                ; use cp of ice Ih for ice VI b/c
                                ; little ice VI available and seems to
                                ; work OK
  dsi = dblarr(n_elements(pmpa))
  si = dblarr(n_elements(pmpa))
  tcp = findgen(355)+1.d0
  ; cp/t for integration
  cpp = (-47.9d0+9.59*tcp+0.0125*tcp*tcp-2.78d-4*(tcp)^3.+1.35d-6*(tcp)^4.-2.65d-9*(tcp)^5.+$
         1.86d-12*(tcp)^6.)
  cpp(274:*) = cpp(273)
  cpt = cpp/tcp

  restore,'meltcurve.sav'

  ; don't need to loop over pressure here b/c S approx const on isotherms
  if t ge 251.16 and t le 354.76 then begin
      dsi(*) = interpol(slt1,tlt,t) ; use S on the phase boundary b/c isotherm approx. eq. isentrope
      sliq = interpol(smeltall(nneg:*),tmeltall(nneg:*),t) ; in MPa
      si(*) = sliq*1.d3-dsi(*)        ; J/kg/K
  end else begin
      ; use V-VI-L tp reference
      if t lt 251.16 then begin
          zzz = where(tcp ge t and tcp le tlt(1))
          dsi(*) = slt1(1) + int_tabulated(tcp(zzz),cpt(zzz),/double,/sort)
          sliq = interpol(smeltall(nneg:*),pmeltall(nneg:*),plt(1)/1.d6) ; in MPa
          si(*) = sliq*1.e3 - dsi(*) ; J/kg/K
      end else begin
          print,'TEMP TOO HIGH FOR ICE VI'
          stop
      end
  end

;  plot,tcp,cpt*tcp
;stop

  ;-------------------------------------------------------------------------------------
                                ; get internal energy from the latent
                                ; heat and the heat capacity at the
                                ; specified pressure
  ei = dblarr(n_Elements(pmpa))
  dei = dblarr(n_Elements(pmpa))
  eliq = interpol(emeltall,pmeltall,pmpa)*1.e3 ;J/kg
  tliq = interpol(tmeltall,pmeltall,pmpa) ; K
  latentheat =  interpol(llt,plt,pmpa*1.d6)
  for i=0,n_elements(pmpa)-1 do begin
      zzz = where(tcp ge t and tcp le tliq(i))
      if zzz(0) ne -1 and n_elements(zzz) gt 1 then dcp = int_tabulated(tcp(Zzz),cpt(Zzz)*tcp(zzz),/double) else dcp = 0.d0
      dei(i) = latentheat(i) + dcp
      ei(i) = eliq(i)-dei(i)             ; J/kg
  end

  ;-------------------------------------------------------------------------------------
  ; CHECK FOR NEGATIVE VALUES!

  ;print,t,pmpa
  ;print,vi,1./vi,dsi,si,eliq,latentheat,dei,ei ; mks
  ;stop

  zzz = where(si le 0.d0)
  if zzz(0) ne -1 then si(zzz) = 1.d-17
  zzz = where(ei le 0.d0)
  if zzz(0) ne -1 then ei(zzz) = 1.d-17

  ; each variable is an array of length of the input pmpa array
  return, [1./vi,si/1.e3,ei/1.e3] ; kg/m3, kJ/kg/K, kJ/kg
  ;return, vi
end

;=========================================================================================
;=========================================================================================
; ICE VII EOS
;
; based on Stewart and Ahrens 2005 using NIST water EOS as reference boundary
; points for S and E
;
; inputs in K and MPa
; t is a scalar
; pmpa can be an array
function icevii, t, pmpa
common setup, first
  ;-------------------------------------------------------------------------------------
  ; FIRST CALC V(T,P)
  ; reference isotherm T=300 K
  tref = 300.d0                   ; K
  v0=1./1460.d0                  ; m3/kg
  Kt0=21.1d9                     ; Pa
  Ktp=4.4d0

  ; thermal expansion from Frank et al 2005
  a0 = 4.2e-4 ; 1/K CHECK SIGN - Frank et al have this as negative - must be an error
  a1 = 1.56e-6 ; 1/K^2
  eta = 1.1d0
                                ; extrapolate alpha(T) at 0 P to zero
                                ; K using the debye function using
                                ; same theta=850k as cp

  if first then begin
      tcparr = [findgen(1000)+2.,(findgen(1500-1000)+2.)*20.+1000.]
      alpha00 = a0+a1*tcparr
  end else begin
      restore,'debye.sav'
      alpha00 = darr/darr(298)*(a0+a1*300.)
  end

  ; range of real densities about 1.4 to 3. g/cm3, 0.33 to 0.7 cm3/g
  vref = v0-(findgen(101)/100.)*0.4*1.d-3
  ; pref = (kt0/ktp)*( (v0/vref)^ktp - 1.d0) ; Murnahan eq.
  xarr  = (vref/v0)^(1./3.d0)
  pref = 3.d0*kt0*(1.-xarr)*xarr^(-2) * exp(3./2.*(ktp-1.)*(1.-xarr)) ; 3rd order B-M eq.
  ;plot,vref*1.e3,pref/1.e9
  ;oplot,vref*1.e3,pref2/1.d9,lines=2
  ; need B-M eq b/c going to high pressures

  ; calc V at a given P, T
  vi = dblarr(n_elements(pmpa))
  for i =0, n_elements(pmpa)-1 do begin
      ; loop over possible array of P values
      p = pmpa(i)*1.d6             ; Pa
      vpt0 = interpol(vref,pref,p)
      if t ne tref then begin
          tarrtmp = (findgen(101)/100.d0)*(t-tref)+tref
          apt0 = interpol(alpha00,tcparr,tarrtmp) ; use alpha0 to extrapolate to high T
          apt = apt0 * (1.d0 + Ktp*p/Kt0)^(-eta)
          tmp = int_tabulated(tarrtmp,apt,/double,/sort)
          if t lt tref then vpt = vpt0*exp(-tmp) else vpt = vpt0*exp(tmp)
      end else begin
          vpt = vpt0
      end
      vi(i) = vpt
  end

  ;-------------------------------------------------------------------------------------
  ; heat capacity of ice VII - Debye function fit by Stewart and Ahrens 2005
; if 1 then begin
 if first then begin
      nnn = 1500 ;v4.3
      tcparr = [findgen(1000)+2.,(findgen(nnn-1000)+2.)*20.+1000.]
                                ; GOOD FOR ZERO INTEGRAL ENTROPY GOOD FOR ENTROPY FIGURE
                                ; 2700 J/kg/K
                                ; need higher cV (Effective) to match different entropy segments on hugoniot
      cvmax = 4600.d0
      Tdeb = 850.d0             ;  K v4.3
      xarr = Tdeb/tcparr
      dfint = (xarr^4.d0)*exp(xarr)/(exp(xarr)-1.d0)^2.d0
      darr = fltarr(n_elements(xarr))
      for i=0,n_elements(xarr)-2 do $;begin
                                ; normalized 3rd order Debye function
          darr(i) = 3.*(xarr(i)^(-3.))*int_tabulated(xarr(i:nnn-1),dfint(i:nnn-1),/sort,/double)
;      end
      cparr = darr*cvmax        ; assume cp=cv for ice VII
      ; make cv stay constant to very high pressrues
      zzz = where(darr eq max(darr))
      if zzz(0) ne -1 then cparr(zzz:*) = cvmax
      save,filename='debye.sav',tcparr,darr
      save,filename='icevii_cp.sav',tcparr,cparr

      ;!p.multi=[0,1,3]
      ;plot,tcparr,cparr,xr=[0,1000],xtit='Temp',ytit='Cv Debye function J/kg/K'
      ;oplot,tarrtmp,(a0+a1*tarrtmp)*1.e6,lines=2
      alpha00 = darr/darr(298)*(a0+a1*300.)
      ;plot,tcparr,darr/darr(298)*(a0+a1*300.),xr=[0,1000]
      ;oplot,tarrtmp,(a0+a1*tarrtmp),lines=2
      save,filename='alpha7_00.sav',tcparr,alpha00
      ;stop
  end else restore,'icevii_cp.sav'
  ; from Fei et al 1993 - cp fit between 300 and 600 K
  ;cpa = 72.49/0.018
  ;cpb = (3.024e-3)/0.018
  ;cpc = (-1.442e6)/0.018
  ;cparr = cpa+cpb*tarr+cpc*tarr^(-2.) ; EXPSTS
  ;plot,tarr,cparr,xr=[0,1000],yr=[0,5000],xtit='Temp',ytit='cV',tit='Ice VII'
  ;oplot,tarr,darr*cvmax,lines=2
  ;oplot,[300.,300],[0,1.e4]
  ;oplot,[600.,600],[0,1.e4]
  ;print,'should be about 2700 J/kg/K'
  ;zzz = where(tarr le 354.d0)
  ;print,int_Tabulated(tarr(zzz),cvmax*darr(zzz)/tarr(zzz),/double,/sort)

  ;-------------------------------------------------------------------------------------
  ; calc S: need to calc S along isotherms

  ; latent heat across VII-liquid - from Dorsey book
  tlt = [81.6,95.3,110.3,124.1,137.1,149.5,161.1,172.1,182.5,192.3]+273.15 ; K
  llt = [354.5,398.0,444.4,474.6,500.1,526.9,554.1,582.6,610.2,642.4]*1000. ; J/kg
  ;pdvlt =[200.0,206.9,215.7,224.6,231.4,239.3,246.1,253.0,258.9,264.8]*1000. ; J/kg
  plt = [21680.,23228.,25164.,27100.,29035.,30971, 32906., 34842.,36778.,38714.]*1.0143*1.e5 ; Pa

  ; extend to very high pressures, temperatures
  tlt = [tlt] ; K
  llt = [llt] ; J/kg
  plt = [plt] ; Pa

  v7lt = dblarr(n_elements(plt))
  vllt = dblarr(n_elements(plt))
  ;for i=0,n_elements(plt)-1 do v7lt(i) = iceVII(tlt(i),plt(i)/1.e6)
  restore,'meltcurve.sav'
  for i=0,n_elements(plt)-1 do vllt(i) = 1./interpol(dmeltall,pmeltall,plt(i)/1.e6)
  ;IDL> print,v7lt
   v7lt = [0.00064976614,   0.00065153325,   0.00065279120,   0.00065348520,   0.00065379853, $
           0.00065382055,   0.00065348720,   0.00065288582,   0.00065202932,   0.00065092763]
  ;IDL> print,vllt
  ;vllt = [0.00074205618,   0.00073720537,   0.00073162643,   0.00072628558,   0.00072175819, $
  ;        0.00071759778,   0.00071245401,   0.00070744022,   0.00070256679,   0.00069783603]

  dvlt = vllt-v7lt

  slt1 = llt/tlt
  slt2 = deriv(tlt,plt)*dvlt
  ;print,slt1,slt2
                                ; these values are similar. Looks like
                                ; dS = 1000 J/kg/K and approximately
                                ; constant

  ; dS across the liquid phase boundary
  ; USE dS=1000 J/kg/K to extrapolate to high pressures
  dsi = dblarr(n_elements(pmpa))
  si = dblarr(n_elements(pmpa))
  dsi0 = 1000.d0 ; APPROX CONSTANT TO EXTRAPOLATE TO HIGH PRESSURES

  ; Calc dS referenced to the melt boundary at the requested P
  restore,'meltcurve.sav'
  sliq = interpol(smeltall(nneg:*),pmeltall(nneg:*),pmpa) ; in MPa
  tliq = interpol(tmeltall(nneg:*),pmeltall(nneg:*),pmpa) ; in MPa

  for i=0,n_elements(pmpa)-1 do begin
      if t gt tliq(i) then begin
          dsi(i) = dsi0         ; outside of phase space
          si(i) = sliq(i)*1.e3 - dsi(i)
      end else begin
          ; need to integrate from the boundary to lower temp
          zzz = where(tcparr ge t and tcparr le tliq(i))
          if zzz(0) ne -1 and n_elements(zzz) gt 1 then begin
              cparr_p = cparr * (1.d0 + Ktp*(pmpa(i)*1.d6)/Kt0)^(-0.3)
              ds_int =  int_Tabulated(tcparr(zzz),cparr_p(zzz)/tcparr(zzz),/double)
          end else ds_int = 0.d0
          dsi(i) = dsi0+ds_int
; v6.4 check if this is the low t problem - not the problem
; this must be a high temp thing....
          if dsi(i) gt sliq(i)*1.e3 then begin
              print,'VII A: ',t,pmpa(0)
              ds_int = (1.-t/tliq(i))*(sliq(i)*1.e3-dsi0)
          end

;          dsi(i) = dsi0+ds_int
          ;v6.4 add check for negative dsi -- not the problem
          if dsi(i) lt 0.d0 then begin
             print,'VII B:', t, pmpa(0)
             dsi(i)=0.d0
          end
          si(i) = sliq(i)*1.e3 - dsi(i)
;          stop
          
; v6.4 THIS IS A FUDGE THAT PREVENTS CROSSING LINES OF ENTROPY in S-T space
         if t lt 251.16 then begin ;v6.4
;            ; use ice VI as a reference???
            ttmp = t
            ptmp = pmpa(0)
            test = icevi(ttmp,ptmp)
; v6.4: this causes a discontinuity in entropy for low temp ice VII
;            si(i) = test(1)*1.e3
; v6.5: try makeing it an interpolation between the ice VI boundary
; and the low end of entropy -- use density as a guide?
            tmpss = (test(1)-.4*(1.d0/vi(i))/5000.)*1.e3
            if tmpss gt 0.d0 then si(i) = tmpss else si(i) = 1.e-20
            
;            print,ttmp,ptmp,test
;            stop
         end
       end
  end

  ;-------------------------------------------------------------------------------------
                                ; get internal energy from the latent
                                ; heat and the heat capacity at the
                                ; specified pressure
  ei = dblarr(n_Elements(pmpa))
  dei = dblarr(n_Elements(pmpa))
  eliq = interpol(emeltall,pmeltall,pmpa)*1.e3 ;J/kg
  tliq = interpol(tmeltall,pmeltall,pmpa) ; K
  latentheat =  1000.d0*tliq ; fixed by the dS=1000 kj/kg on the melting curve
  for i=0,n_elements(pmpa)-1 do begin
      zzz = where(tcparr ge t and tcparr le tliq(i))
      if zzz(0) ne -1 and n_elements(zzz) gt 1 then begin
          cparr_p = cparr * (1.d0 + Ktp*(pmpa(i)*1.d6)/Kt0)^(-0.3)
          dcp = int_tabulated(tcparr(zzz),cparr_p(Zzz),/double)
      end else dcp = 0.d0
      ; check for negative energy v5
      dei(i) = latentheat(i) + dcp
      if dei(i) gt eliq(i) then dcp = (1.-t/tliq(i))*(eliq(i)-latentheat(i)) ; assume proportional to temp
; stop
      dei(i) = latentheat(i) + dcp
      ei(i) = eliq(i)-dei(i)             ; J/kg
  end

  ;-------------------------------------------------------------------------------------
  ; CHECK FOR NEGATIVE VALUES!

  zzz = where(ei le 0.d0)
if zzz(0) ne -1 then begin
  print,t,pmpa
  print,vi,1./vi,sliq*1.e3,dsi,si,eliq,dei,ei ; mks
  stop
end

  zzz = where(si le 0.d0)
  if zzz(0) ne -1 then si(zzz) = 1.d-17
  zzz = where(ei le 0.d0)
  if zzz(0) ne -1 then ei(zzz) = 1.d-17

  ; each variable is an array of length of the input pmpa array
  return, [1./vi,si/1.e3,ei/1.e3] ; kg/m3, kJ/kg/K, kJ/kg
  ;return,vi
end

;=========================================================================================
;=========================================================================================

function IhEOS, t, p ; in K and MPa
ppa = p*1.e6 ; convert to Pa for calc

; From Feistel and Wagner 2005
g00 = -632020.233449497d0 ; J kg1
g01 = 0.655022213658955d0 ; J kg1
g02 = -1.89369929326131D-08 ; J kg1
g03 = 3.39746123271053D-15 ; J kg1
g04 = -5.56464869058991D-22 ; J kg1
s0  = 189.13d0 ; absolute J kg1 K1
;s0  = -3327.33756492168d0 ; IAPWS-95 J kg1 K1
t1  = complex(3.68017112855051E-02, 5.10878114959572E-02,/double)
r1  = complex(44.7050716285388, 65.6876847463481,/double); J kg1 K1
t2  = complex(0.337315741065416, 0.335449415919309, /double)
r20 = complex(-72.597457432922, -78.100842711287,/double) ; J kg1 K1
r21 = complex(-5.57107698030123E-05, 4.64578634580806E-05,/double); J kg1 K1
r22 = complex(2.34801409215913E-11, -2.85651142904972E-11,/double) ; J kg1 K1

;density(T,P) = 1.d0/gp
pi0=101325.d0/611.657d0
pt = 611.657d0 ; Pa
pi = ppa/pt
tt = 273.16d0 ; K
tau = t/tt

g0p = g01*(1.d0/pt)*(pi-pi0)^(1.d0-1.d0) + g02*(2.d0/pt)*(pi-pi0)^(2.d0-1.d0) + g03*(3.d0/pt)*(pi-pi0)^(3.d0-1.d0) + g04*(4.d0/pt)*(pi-pi0)^(4.d0-1.d0)
r2p = r21*(1.d0/pt)*(pi-pi0)^(1.d0-1.d0) + r22*(2.d0/pt)*(pi-pi0)^(2.d0-1.d0)
gp  = g0p + tt*real_part(r2p*( (t2-tau)*alog(t2-tau) + (t2+tau)*alog(t2+tau) - 2.d0*t2*alog(t2) - (tau^2.d0)/t2))
density = 1.d0/gp

; specific entropy = -gt
r2 = r20*(pi-pi0)^0.d0 + r21*(pi-pi0)^1.d0 + r22*(pi-pi0)^2.d0
icegt = -s0 + real_part(r1*(-alog(t1-tau)+alog(t1+tau)-2.d0*tau/t1) + r2*(-alog(t2-tau)+alog(t2+tau)-2.d0*tau/t2))
entropy = -icegt

; specific internal energy = g-t*gt-p*gp
g0 = g00*(pi-pi0)^0.d0 + g01*(pi-pi0)^1.d0 + g02*(pi-pi0)^2.d0 + g03*(pi-pi0)^3.d0 + g04*(pi-pi0)^4.d0
g  = g0 - s0*tt*tau + tt*real_part(r1*( (t1-tau)*alog(t1-tau) + (t1+tau)*alog(t1+tau) - 2.d0*t1*alog(t1) - (tau^2.)/t1) + $
      r2*( (t2-tau)*alog(t2-tau) + (t2+tau)*alog(t2+tau) - 2.d0*t2*alog(t2) - (tau^2.)/t2) )
energy = g - t*icegt - ppa*gp
enthalpy = g - t*icegt

Hshift=632128.7427016902d0 ; value at TP
Hmelt = 333430.d0 ; J/kg from Feistel and Wagener 2005
;print,t,p,density,entropy,energy+hshift,enthalpy+hshift
;;stop

; volume expansion coef
; alpha = gtp/gp
gtp = real_part(r2p*( -alog(t2-tau)+alog(t2+tau)-2.d0*tau/t2))
alpha = gtp/gp

; isothermal modulus
; -gpp/gp
g0pp = g02*2.*(2.-1.d0)/pt*(pi-pi0)^(2.d0-2.d0) + $
       g03*3.*(3.-1.d0)/pt*(pi-pi0)^(3.d0-2.d0) + $
       g04*2.*(4.-1.d0)/pt*(pi-pi0)^(4.d0-2.d0)
gpp = g0pp + tt*real_part( $
                           r1*(-alog(t1-tau)+alog(t1+tau)-2.d0*tau/t1) + $
                           r2*(-alog(t2-tau)+alog(t2+tau)-2.d0*tau/t2) )
kt = -gpp/gp

; units: kg/m3, kJ/kg/K, kJ/kg, kJ/kg, 1./K, 1./Pa
return,[density,entropy/1.d3,(energy+hshift)/1.d3,(enthalpy+hshift)/1.d3,alpha,kt]
end

;=========================================================================================
;=========================================================================================

; Feistel and Wagner 2005
function icesub, t ; in K
  pt = 611.657d0 ; Pa, triple point pressure
  tt = 273.16d0 ; K, triple point temp
  r  = 461.52364d0 ; J/kg/K
  dht = 2834.4D3 ; J/kg

  psub = (pt * exp( (dht/r) * (1.d0/tt - 1.d0/(t*1.d0)) ))/1.D6 ; MPa
  zzz = where(psub lt 1.D-30)
  if zzz(0) ne -1 then psub(zzz)=1.D-30
  return, psub ; in MPa
end

;=========================================================================================
;=========================================================================================

function nisteos, t, p ; K, MPa
  ;Eshift,Sshift=       632.11448       3.5233498
  cmd='./eostp '+strtrim(string(t),2)+' '+strtrim(string(p)) ; K MPa
  spawn,cmd,res
  tmp = strsplit(res(0),/extract)
  if tmp(3) ne 0 then stop ; error in density calc
  d = 1.d0*tmp(0) ; kg/m3
  e = 1.d0*tmp(1)+632.11448d0 ; kJ/kg
  s = 1.d0*tmp(2)+3.5233498d0 ; kJ/kg/K
  return, [d,e,s]
end

;=========================================================================================
;=========================================================================================
function b67, t
  ; calc the pressure of the ice VI-VII phase boundary at T (K)
  common boundaries, ttp6,ptp6,ttp,ptp,tcrit,pcrit,ttp3,ptp3
  p = 1.*(t-ttp6)+ptp6 ; MPa
  return, p ; in MPa
end

function b16, t
  ; calc the pressure of the ice Ih-VI phase boundary at T (K)
  common boundaries, ttp6,ptp6,ttp,ptp,tcrit,pcrit,ttp3,ptp3
  p = 0.2*(t-ttp3)+ptp3 ; MPa
  return, p ; in MPa
end

;=========================================================================================
;=========================================================================================
; END OF FILE
;=========================================================================================
;=========================================================================================
