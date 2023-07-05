function readses, tablename

openr,unit,tablename,/get_lun
; READ INDEX SECTION
line=''
;readf,unit,matid,nwds0,format="(' INDEX     MATID =',I7,'     NWDS = ',I5)"
readf,unit,line
readf,unit,matid,date,date,version,ntables,format="(5(E16.8E2))"
readf,unit,t1,t2,nwds1,nwds,format="(4(E16.8E2))"

; WRITE OUT 201 TABLE SECTION
;readf,unit,t1,nwds1,format="(' RECORD',5X,'TYPE =',I5,5X,'NWDS = ',I5)"
readf,unit,line
readf,unit,fmn,fmw,rho0_ref,k0_ref,t0_ref,format="(5(E16.8E2))"

; WRITE OUT 301 TABLE SECTION
;readf,unit,t2,nwds,format="(' RECORD',5X,'TYPE =',I5,5X,'NWDS = ',I6)"
readf,unit,line

; read in everything and split up after
linarr = dblarr(5)
allarr = [-1.d0]
readf,unit,line
while not EOF(unit) do begin
    linarr = double(strsplit(line,/extract))
    if n_elements(linarr) ne 5 then begin
        tmp =n_elements(allarr)
        print,tmp,tmp-allarr(0)-allarr(1)-2
        stop
    end
    allarr = [allarr,linarr]
    readf,unit,line
end
allarr = allarr[1:*]
;plot,allarr
;stop

; first: # density points, # temperature points in mesh
;readf,unit,dsize,format="(E16.8E2,$)"
;readf,unit,tsize,format="(E16.8E2,$)"
dsize = allarr(0)
tsize = allarr(1)

darr = dblarr(dsize) ; density g/cm3
tarr = dblarr(tsize) ; K
parr = dblarr(tsize,dsize) ; pressure in GPa
earr = dblarr(tsize,dsize) ; internal energy in MJ/kg
sarr = dblarr(tsize,dsize) ; specific entropy in MJ/kg/K

; second: density array, temperature array
inext = 2
darr(*) = allarr(inext:inext+dsize-1)
inext = inext+dsize
tarr(*) = allarr(inext:inext+tsize-1)
inext = inext+tsize

;print,darr,tarr
;stop
; third: pressure, loop over temperture, loop over density
for i=0,tsize-1 do begin
;    for j=0,dsize-1 do begin
        parr(i,*) = allarr(inext:inext+dsize-1)
        inext = inext+dsize
;    end ; d loop
end ; t loop
; third: energy, loop over temperture, loop over density
for i=0,tsize-1 do begin
;    for j=0,dsize-1 do begin
        earr(i,*) = allarr(inext:inext+dsize-1)
        inext = inext+dsize
;    end ; d loop
end ; t loop
if inext lt n_elements(allarr) then begin
; fourth: entropy, loop over temperture, loop over density
    for i=0,tsize-1 do begin
;    for j=0,dsize-1 do begin
        sarr(i,*) = allarr(inext:inext+dsize-1)
        inext = inext+dsize
;    end                         ; d loop
    end                         ; t loop
end
close,unit
free_lun,unit

tijarr = dblarr(tsize,dsize) 
dijarr = dblarr(tsize,dsize) 
for i=0,tsize-1 do begin
    tijarr(i,*) = tarr(i)
end
for j=0,dsize-1 do begin
    dijarr(*,j) = darr(j)
end


print,inext, n_elements(allarr)


!p.multi=[0,1,2]
plot,darr,tarr,/nodata,xtit='Density (g/cm3)',ytit='Temperature (K)',yr=[0,2000],xr=[0,4]
for i=0,tsize-1 do for j=0,dsize-1 do $
  plots,darr(j),tarr(i),psym=6,symsize=.2

plot,darr,parr,/nodata,xtit='Density (g/cm3)',ytit='Pressure (GPa)',yr=[-10.,200.],xr=[0,4]
for i=0,tsize-1 do for j=0,dsize-1 do $
  plots,darr(j),parr(i,j),psym=6,symsize=.2
zzz = where(parr lt 0.d0)
;print,dijarr(zzz)
;print,tijarr(zzz)

plot,dijarr(zzz),parr(zzz),xtit='Density (g/cm3)',ytit='Pressure (GPa)',psym=6,symsize=.2;,yr=[-10.,200.],xr=[0,4]


plot,tijarr(zzz),earr(zzz),xtit='Temperature (K)',ytit='Energy',psym=6,symsize=.2;,yr=[-10.,200.],xr=[0,4]

!p.multi=0
device,decomposed=0
loadct,37
zzz = where(parr lt 0.d0)
plot,sarr(zzz),tijarr(zzz),ytit='Temperature (K)',xtit='Entropy',psym=6,symsize=.2,col=255,yr=[0,500];,xr=[0,0.01]
yyy = where(parr gt 0.d0 and parr le 0.5d0)
oplot,sarr(yyy),tijarr(yyy),psym=4,symsize=.2,col=50


!p.multi=0
device,decomposed=0
loadct,37
zzz = where(parr lt 0.d0)
plot,dijarr(zzz),parr(zzz),ytit='pressure',xtit='density',psym=6,symsize=.2,col=255,yr=[-5,5.];,xr=[0,0.01]
yyy = where(parr gt 0.d0 and parr le 5.d0)
oplot,dijarr(yyy),parr(yyy),psym=4,symsize=.2,col=50


!p.multi=0
device,decomposed=0
loadct,37
zzz = where(parr eq 0.d0)
plot_io,sarr(zzz),tijarr(zzz),ytit='Temperature (K)',xtit='Entropy',psym=6,symsize=2,col=255;,yr=[0,2000],xr=[0,0.01]
yyy = where(parr gt 0.d0 and parr le 0.5d0)
oplot,sarr(yyy),tijarr(yyy),psym=4,symsize=.2,col=150

kkk = where(tijarr(zzz) eq 900.)
xxx = zzz(kkk)
!p.multi=[0,1,3]
plot,dijarr(xxx),parr(xxx),psym=-7
plot,dijarr(xxx),earr(xxx),psym=-7,/yno
plot,dijarr(xxx),sarr(xxx),psym=-7,/yno
stop



d=2.d0
;d=0.0010000000d0
;d=1.d-6
;d=3.d0
zzz = where(parr le 0.d0)
yyy = where(parr gt 0.d0 and parr le 5.d0)
kkk = where(dijarr(zzz) eq d)
xxx = zzz(kkk)
kkk = where(dijarr(yyy) eq d)
www = yyy(kkk)
xxx = [xxx,www]
!p.multi=[0,1,3]
!p.charsize=2.
plot_oi,tijarr(xxx),parr(xxx),psym=-7,tit='Density = '+string(d),xr=[100.,1.2*(max(tijarr(xxx)))]
oplot,tijarr(www),parr(www),psym=-4,col=100
plot_oi,tijarr(xxx),earr(xxx),psym=-7,/yno,xr=[100,1.2*(max(tijarr(xxx)))]
oplot,tijarr(www),earr(www),psym=-4,col=100
plot_oi,tijarr(xxx),sarr(xxx),psym=-7,/yno,xr=[100,1.2*(max(tijarr(xxx)))]
oplot,tijarr(www),sarr(www),psym=-4,col=100


;------------
d=2.d0
;d=0.0010000000d0
;d=1.d-6
;d=3.d0
zzz = where(parr le 0.d0)
yyy = where(parr gt 0.d0 and parr le 5.d0)
kkk = where(dijarr(zzz) eq d)
xxx = zzz(kkk)
kkk = where(dijarr(yyy) eq d)
www = yyy(kkk)
xxx = [xxx,www]
device,decomposed=0
!p.multi=[0,1,3]
!p.charsize=2.
loadct,37
plot_oi,tijarr(xxx),parr(xxx),psym=-7,tit='Density = '+string(d),xr=[100.,1.2*(max(tijarr(xxx)))]
oplot,tijarr(www),parr(www),psym=-4,col=100

d=0.53182960d0
zzz = where(parr le 0.d0)
yyy = where(parr gt 0.d0 and parr le 5.d0)
kkk = where(dijarr(zzz) eq d)
xxx = zzz(kkk)
kkk = where(dijarr(yyy) eq d)
www = yyy(kkk)
xxx = [xxx,www]
oplot,tijarr(xxx),parr(xxx),psym=7,col=200
oplot,tijarr(www),parr(www),psym=-4,col=200

d=1.0313386d0
zzz = where(parr le 0.d0)
yyy = where(parr gt 0.d0 and parr le 5.d0)
kkk = where(dijarr(zzz) eq d)
xxx = zzz(kkk)
kkk = where(dijarr(yyy) eq d)
www = yyy(kkk)
xxx = [xxx,www]
oplot,tijarr(xxx),parr(xxx),psym=7,col=50
oplot,tijarr(www),parr(www),psym=-4,col=50

;d=2.5999999d0
;d=2.6480000d0
d=1.0313386d0
zzz = where(parr le 0.d0)
yyy = where(parr gt 0.d0 and parr le 5.d0)
kkk = where(dijarr(zzz) eq d)
xxx = zzz(kkk)
kkk = where(dijarr(yyy) eq d)
www = yyy(kkk)
xxx = [xxx,www]
oplot,tijarr(xxx),parr(xxx),psym=7,col=150
oplot,tijarr(www),parr(www),psym=-4,col=150


plot_oi,tijarr(xxx),earr(xxx),psym=-7,/yno,xr=[100,1.2*(max(tijarr(xxx)))]
oplot,tijarr(www),earr(www),psym=-4,col=100
plot_oi,tijarr(xxx),sarr(xxx),psym=-7,/yno,xr=[100,1.2*(max(tijarr(xxx)))]
oplot,tijarr(www),sarr(www),psym=-4,col=100


stop
stop

end
