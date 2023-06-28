c     Fortran code by John Borg
c     v11c - updated by STS for I/O comparisons to pyKOv0.1; initial minimalist code match
c     v11d - checking d=2, d=3 and ibc +- 4 with pyKOv0.2
c     v11d - I/O tweaks for clean test cases with pyKOv0.3; added etot to main output
c     January 30, 2023, updated February 2, 2023; updated 2/10/23
c     KO see Springer book by Wilkins, M.
c     "Computer simulations of Dynamic Phenomena"
      PROGRAM fKO
       parameter(jj = 502)       !total number of nodes allowed
                                 !this number should be greater than the number
                                 !of nodes that you want so that nodes can be added 
                                 !if fracture occurs
       parameter(jmat=7)        !number of materials
       parameter(nmax=4)        !number of placeholders for time: t(0) and t(1) current, t(2) and t(3) new; keeps memory very small 
!! fortran is crazy and lets you declare the first index 0 instead of 1
!     ! all these deckared variables have first index 0
       integer narg
       integer n, j,imat,ibc(0:jj), icontact,idebug,jdebug,iskip,ipoint
       integer j1(jj),iEOS(2,0:jj),jjj,jzz,ndebug,i,amr
       integer icompact(0:nmax,jj), mat(0:jj), stepn

       real*8  m(0:jj), r(0:nmax,0:jj), U(0:nmax,0:jj)
       real*8  t(0:nmax),pfrac(0:jj),pvoid
       real*8  phi(0:nmax,0:jj), sigmar(0:nmax,0:jj)
       real*8  sigmao(0:nmax,0:jj),Temp(0:nmax,0:jj)
       real*8  beta(0:nmax,0:jj), P(0:nmax,0:jj), q(0:nmax,0:jj)
       real*8  s1(0:nmax,0:jj),s2(0:nmax,0:jj), s3(0:nmax,0:jj)
       real*8  rho(0:jj), V(0:nmax,0:jj), entropy(0:nmax,0:jj)
       real*8  epsi1(0:nmax,0:jj),epsi2(0:nmax,0:jj),epdt1(0:nmax,0:jj)
       real*8  E(0:nmax,0:jj), K(0:nmax,0:jj), Y(0:jj),dtminj(0:jj)
       real*8  deltaZ, Vdot ,dt_min, dr_min,delt_temp,deltaZarr(0:jj)
       real*8  deltar,  qbar
       real*8  a, b, Co, CL, d ,a_min,b_min,rho_min
       real*8  alocal(0:jj),rholocal(0:jj)
       real*8  EOS(0:jmat,13),init(jmat,16),bc(-9:9,5),LL(jj),xstart(jj)
       real*8  xx,xa,xb,deltat_0,deltat,delt,tstop
c       real*8  U_0,P_0,V_0,E_0,rho_0,a_0,up
c       real*8  U_1,P_1,V_1,E_1,a,
       real *8 zero
       real*8  aa,ww,rr,bb,dti0,dti1,Usave
       real*8  phiv,phij,betav,betaj,dthalf
       real*8  qtotal,mvtotal,ketotal,ietotal,etotal,ke3total
       real*8  bs1,bs2,bs3,bs4,bs5,bs6,bs7,bs8,bs9,bs10
       real*8  tskip,dtskip
       real*8  rho_local,stemp,ctemp,gtemp,v0,v00,vv
c       real*8  alpha
c      real*8  pe,ps,ap,ae
       real*8  k1,k2,k3,gamma0
       real*8  Us,up,PH,EH,TH,strain,P0,E0,T0
       real*8  En2j1,diffE
       real*8  startTime,endTime
       character*20 name
       character*172 text
       character*40 fin, fout,arg,foutlong

       fin = 'none'
       fout = 'none'
       foutlong = 'none'
       DO i = 1, iargc()
          CALL getarg(i, arg)
C          WRITE (*,*) arg
          if (i .eq. 1) then
            fin = trim(arg)
          endif
          if (i .eq. 2) then
             fout = trim(arg)
          endif
          if (i .eq. 3) then
             read(arg,'(f5.3)')dtskip
          endif
       END DO
       if (fin .eq. 'none') then
          fin = 'ko.in'
       endif
       if (fout .eq. 'none') then
          fout = 'ko-fort.dat'
       endif
       write (*,*) fin, fout, dtskip
       !read(*,*) foo
c
c     Initalize some key varables
       Co        = 2.0d0  !artifical viscoity constants, default=2 Von Neumann formulation
       CL        = 1.0d0  !artifical viscoity constants, default=1
       pvoid     = 0.d0   !pressure in the void
       icontact  = 0      !contact flag this should not be changed but i think it could
                          ! be elimated if the code were cleaned up.
       idebug    = 0      ! turn on debug
       ndebug    = 1      ! what time steup to you want to print to screen
       jdebug    = 2      ! what node do you want to look at when debugging
       d         = 0.d0   !defines geometry (1-1D planar, 2-1D cylinderical-never tested), see willkins
       deltat_0  = 1.0d-3 ! initial time step to get going microsec
c
       iskip     = 1000    !number of iterations to skip to echo to screen
       dtskip    = 0.01  ! This is the amount of time skipped between writes to file this is 3rd argument
       tskip     = dtskip-deltat_0     ! this is the value at which data starts getting written to file
       amr       = 0      ! this is an adaptive mesh refinement (1 is on - 0 is off)
c---------------------------------------------------------------------
c     zero all variables
c
c      print *,'Initializing variables ...'
      zero = 0.d0
      stepn=0
       do n = 0,nmax
          do j = 0,jj
             ibc(j)      = 9  ! these are non-used nodes
             U(n,j)      = zero
             r(n,j)      = zero
             t(n)        = zero ! sts changed deltat_0*dfloat(n)
             deltat      = deltat_0
             phi(n,j)    = zero
             sigmar(n,j) = zero
             sigmao(n,j) = zero
             beta(n,j)   = zero
             q(n,j)      = zero
             s1(n,j)     = zero
             s2(n,j)     = zero
             s3(n,j)     = zero
             epsi1(n,j)  = zero
             epsi2(n,j)  = zero
             K(n,j)      = zero
             Y(j)        = zero
             deltaZ      = zero
             pfrac(j)    = zero
             Temp(n,j)   = zero
             entropy(n,j)= zero
             mat(j)      = 0 ! stsm int
             dtminj(j)   = zero ! stsm
             alocal(j)   = zero ! stsm
             rholocal(j) = zero ! stsm
            enddo
       enddo
c
      do j=1,jmat
       do i=1,7
        eos(j,i)= 0.d0
       enddo
       do i=1,4
        init(j,i)= 0.d0
       enddo
      enddo
      do j=1,jj
       do jzz=1,2
        ieos(jzz,j)= 0
       enddo
      enddo
c---------------------------------------------------------------------
c     Read input Data file
c
c      print *,'running'
c      read(*,*) foo
      imat = 0
      jsum = 0

c      write(*,*) "test"
c      open (unit=31,file='ko.in',form='formatted',status='unknown')
      open (unit=31,file=fin,form='formatted',status='unknown')
c      open (unit=35, file='fko-energy.dat', 
c     &      form='formatted', status='unknown')
      print *,'input file: ',fin
c      read(*,*) foo
      Read(31,'(A149)') text
c      print *,text
      Read(31,'(A149)') text
c      print *,text
c     read(*,*) foo
      print *,'EOS input'      
 5    imat=imat+1
c      print *,'reading',imat
c      read(*,*) foo
      Read(31,'(A149)') text
c      print *,Text
c             print *,'testing text'
c             read(*,*) foo
      if(text(1:7) .eq. '       ') goto 6
c     print *,'text passed.'
      iEOS(1,jj-(imat-1)) = imat
      Read(text,'(2I7,7e7.3,3e9.3,10e7.3)') !g77 ! reads until there is a blank line and then goto exits the loop
c      Read(text,*)  !g90
     & iEOS(2,jj-(imat-1)),j1(imat),LL(imat),xstart(imat),
     & init(imat,1),init(imat,2),init(imat,3),
     & init(imat,4),init(imat,5),
     & eos(imat,1),eos(imat,2),eos(imat,3),eos(imat,4),
     & eos(imat,5),eos(imat,6),eos(imat,7),eos(imat,8),eos(imat,9),
     & eos(imat,10),eos(imat,11),eos(imat,12),eos(imat,13)
      if ( int(j1(imat)/2) .ne. float(j1(imat))/2.) then
       print *,'Error: nodes specificed for a material must be even'
        read(*,*) foo
      endif
      mat(jsum:jsum+j1(imat)) = imat
      jsum = jsum + j1(imat)
      if (jsum .ge. jj) then
       print *,'Error: Increase jj or reduce the # nodes in kov3.in'
       read(*,*) foo
      endif
C      print *,'data at j',imat,eos(imat,7)
C      print *,'gg',eos(imat,4)
      write (*,'(3I4,13f11.4,2e12.4e2)')
     & iEOS(1,jj-(imat-1)),
     & iEOS(2,jj-(imat-1)),j1(imat),LL(imat),xstart(imat),
     & init(imat,1),init(imat,2),init(imat,3),
     & init(imat,4),init(imat,5),
     & eos(imat,1),eos(imat,2),eos(imat,3),eos(imat,4),
     &     eos(imat,5),eos(imat,6),eos(imat,7),eos(imat,8)
c     & eos(imat,9),
c     & eos(imat,10),eos(imat,11),eos(imat,12)
c      read(*,*) foo
      goto 5
 6    Print *,' '
      imat = imat-1
c      read(*,*) foo
c
c---------------------------------------------------------------------
c     Input boundary conditions
c
      Read(31,'(A113)') text
c      print *,text
c      read(*,*) foo
      Read(31,'(A113)') text
c      print *,text
c      read(*,*) foo
      Read(31,'(A113)') text
c      print *,text
      Read(text,'(I7,5e7.3)') ! g77
c      Read(text,*) ! g90
     &     ibc(0),bc(-1,1),bc(-1,2),bc(-1,3),bc(-1,4),bc(-1,5)
      print *,'Boundary Conditions'
      write(*,'(I4,5f7.3)')  ! g77
     & ibc(0),bc(-1,1),bc(-1,2),bc(-1,3),bc(-1,4),bc(-1,5) !note that ibc is stored at jjj just temporarily
c      print *,'read Boundary Conditions...',ibc(0)
c      read(*,*) foo
      Read(31,'(A150)') text
c      print *,text
      Read(text,'(I7,5e7.3)')  ! g77
     & ibc(jj),bc( 1,1),bc( 1,2),bc( 1,3),bc( 1,4),bc( 1,5) !note that ibc is stored at jjj just temporarily
c     print *,'read Boundary Conditions...',ibc(jj)
      write(*,'(I4,5f7.3)')  ! g77
     & ibc(jj),bc( 1,1),bc( 1,2),bc( 1,3),bc( 1,4),bc( 1,5) !note that ibc is stored at jjj just temporarily
c
      Read(31,'(A113)') text
c      print *,text
      Read(31,'(A113)') text
c      print *,text
      Read(31,'(e20.3)') tstop ! g77
c      print *,'tstop = ',tstop
      Read(31,'(A113)') text
c      print *,text
      Read(31,'(e20.3)') d    ! g77 geometry
c      print *,'d = ',d
      if (d .ne. 1.0) then
         if (d .ne. 2.0) then
            if (d .ne. 3.0) then
               print *,'INVALID GEOMETRY d=',d
               print *,'STOP'
               read(*,*) foo
            endif
         endif
      endif
      !read(*,*) foo ! a debugging stop
      close(31)
c
c     STS output file format for comparison to pyKO 
      print *,'output file: ',fout
      open (unit=33, file=fout, form='formatted', status='unknown')
      write(33,*) "step j ibc mat time r(0) r(1) r(2) pos 
     & up vr rho rho0 iev0 pres sigmar s1 s2 q dtminj phi beta
     & eps1 eps2 deltaz aj rhoj sigmao temp etot"
c      open (unit=34, file=foutlong, form='formatted', status='unknown')
c      write(34,*) "step j ibc mat time r(0) r(1) r(2) pos 
c     & up vr rho rho0 iev0 pres sigmar s1 s2 q dtminj phi beta
c     & eps1 eps2 deltaz aj rhoj sigmao temp etot"
cc---------------------------------------------------------------------
c      Distritazation LOOP -- assign initial mass
c
C      print *,'distrize ..'
c      read(*,*) foo
        r(0,1) = xstart(1) ! stsm these 2 statement are overwritten below
        r(1,1) = xstart(1)
      do jjj=1,imat
C       print *,'distrize ..',jjj
c       read(*,*) foo
       deltar = Ll(jjj)/dfloat(j1(jjj)/2)
C       write(*,'(a6,1x,3e24.16)') 'deltar',deltar
c       read(*,*) foo
       if (jjj .eq. 1) then                     ! this if assigns the ipoint to the first node of the material
                                                !  and checks to see it materials are initally in contact.
        ipoint = 0
       elseif (abs(r(0,ipoint+j1(jjj-1))-xstart(jjj)) .lt. 1.e-5) then
c       no gap between materials
        r(0,ipoint+j1(jjj-1)) = xstart(jjj)
        r(1,ipoint+j1(jjj-1)) = xstart(jjj)
        r(2,ipoint+j1(jjj-1)) = xstart(jjj)
        r(3,ipoint+j1(jjj-1)) = xstart(jjj)
        ibc(ipoint+j1(jjj-1)) = 0
        ipoint                = ipoint+j1(jjj-1)
       elseif ( r(0,ipoint+j1(jjj-1)) .lt. xstart(jjj)) then
c       inital gap between materials
        ipoint      = ipoint+j1(jjj-1)+2
        ibc(ipoint  ) = -2
        ieos(2,ipoint-1)= 0
        ibc(ipoint-2) =  2
       else
       print *,'Input Geometry Error!'
       read(*,*) foo
       endif
c         ! here if ipoint=0 for first material
          r(0,ipoint)  = xstart(jjj) ! jjj=1
          r(1,ipoint)  = xstart(jjj)
          r(2,ipoint)  = xstart(jjj) ! jjj=1
          r(3,ipoint)  = xstart(jjj)
          ibc(1)       = 0
          U(0,ipoint)  = init(jjj,2)
          U(1,ipoint)  = init(jjj,2)
          U(2,ipoint)  = init(jjj,2)
          U(3,ipoint)  = init(jjj,2)
       do j=ipoint+2,ipoint+j1(jjj),2      !node definitions
          ibc(j)       = 0      !this defines the node as a centeral difference
          r(0,j)       = r(0,j-2) + deltar
          r(1,j)       = r(1,j-2) + deltar
          r(2,j)       = r(2,j-2) + deltar
          r(3,j)       = r(3,j-2) + deltar
          U(0,j)       = init(jjj,2)
          U(1,j)       = init(jjj,2)
          U(2,j)       = init(jjj,2)
          U(3,j)       = init(jjj,2)
       enddo
c       read(*,*) foo
       do j=ipoint+1,ipoint+j1(jjj)-1,2    !cell definitions initial conditions
          ibc(j)       = 0      !this defines the cell as a centeral difference
          r(0,j)       = 0.5d0*(r(0,j+1)+r(0,j-1))
          r(1,j)       = 0.5d0*(r(1,j+1)+r(1,j-1))
          P(0,j)       = init(jjj,1)
          P(1,j)       = init(jjj,1)
          P(2,j)       = init(jjj,1) ! added by sts
          P(3,j)       = init(jjj,1) ! added by sts
          rho(j)       = init(jjj,3)  !rho is not updated in time therefore it is always rho0
          E(0,j)       = init(jjj,4)  ! this gets overwritten below
          E(1,j)       = init(jjj,4)
          ieos(1,j)    = ieos(1,jj-(jjj-1))  !info was stored there (RHS) just temp.
          ieos(2,j)    = ieos(2,jj-(jjj-1))  !info was stored there (RHS) just temp.
          Y(j)         = eos(jjj,5)
          pfrac(j-1)   = eos(jjj,7)
          pfrac(j)     = eos(jjj,7)
          pfrac(j+1)   = eos(jjj,7)
       enddo
          ieos(1,jj-(jjj-1))     = 0   !this assigns values to the last cell center
          ieos(2,jj-(jjj-1))     = 0
          ibc(j1(jjj)+ipoint)    = 2
          ibc(j1(jjj)+1+ipoint)  = 9
      enddo  ! jjj loop through imat
      ibc(j1(imat)+ipoint)    = ibc(jj)  !this assigns values to the last node
c      rho(j1(imat)+ipoint+1)    = init(imat,3) !rho is not updated in time therefore it is always rho0

c      print *,'ipoint',j1(imat)+ipoint,ibc(jj),
c     c   init(imat-1,3),init(imat,3)
c
C      print *,'Assign mass to nodes'
c
        do j=0,jj-2,2
             if (ibc(j+1) .eq. 0) then
              m(j+1) = rho(j+1)*((r(0,j+2)**d- r(0,j)**d)/d)
             endif
        enddo
c        read(*,*) foo
c---------------------------------------------------------------------
c set initial specific volume
      do j=0,jj-2,2
         if (ibc(j+1) .eq. 0) then
        V(0,j+1)=rho(j+1)*((r(0,j+2)**d-r(0,j)**d)/d)/m(j+1)
        V(1,j+1)=rho(j+1)*((r(0,j+2)**d-r(0,j)**d)/d)/m(j+1)
c
        E(0,j+1)=P(0,j+1)/(v(0,j+1)*eos(ieos(1,j+1),4)-1)
        E(1,j+1)=P(1,j+1)/(v(1,j+1)*eos(ieos(1,j+1),4)-1)
        Temp(0,j+1) =E(0,j+1)/(rho(j+1)*eos(ieos(1,j+1),8))
        Temp(1,j+1) =E(1,j+1)/(rho(j+1)*eos(ieos(1,j+1),8))
        entropy(0,j+1) = 6.8d-5 ! i got this from a VT website EES thing
        entropy(1,j+1) = 6.8d-5 ! units mbar-cc/K/g

        ! sts added for first time step output
        V(2,j+1)=rho(j+1)*((r(0,j+2)**d-r(0,j)**d)/d)/m(j+1)
        V(3,j+1)=rho(j+1)*((r(0,j+2)**d-r(0,j)**d)/d)/m(j+1)
        E(2,j+1)=P(0,j+1)/(v(0,j+1)*eos(ieos(1,j+1),4)-1)
        E(3,j+1)=P(1,j+1)/(v(1,j+1)*eos(ieos(1,j+1),4)-1)
        Temp(2,j+1) =E(0,j+1)/(rho(j+1)*eos(ieos(1,j+1),8)) ! this needs a loop for each material or an array with cvs STS DEBUG; currently wrong
        Temp(3,j+1) =E(1,j+1)/(rho(j+1)*eos(ieos(1,j+1),8))
C        print *,'Pres [Mbar]',P(0,j+1)
C        print *,'Engr [MBar-cc/cc]',E(0,j+1)
C        print *,'Volm [cc]',((r(0,j+2)**d-r(0,j)**d)/d)
C        print *,'temp*[K]',Temp(0,j+1),rho(j+1),m(j+1),v(0,j+1)
C     print *,'entropy',entropy(0,j+1),eos(ieos(1,j+1),8)
c        print *,'rho(j+1)',rho(j+1),Temp(3,j+1),E(3,j+1)
c        read(*,*) foo
C        print *,'here1',j
c        icompact(0,j-2)  = 0 ! compaction flag for use with snow plow type models initilized to zero
c      print *,'here2',j
c        icompact(1,j-2)  = 0 ! compaction flag for use with snow plow type models initilized to zero
c      print *,'here3',j
c        icompact(2,j-2)  = 0 ! compaction flag for use with snow plow type models initilized to zero
c        print *,'i am here'
c	read(*,*) foo
       endif
      enddo
c---------------------------------------------------------------------
c     Output Initial Conditions to screen
c
      n=1
C      write(*,98)
C     &'j','ibc',' rad',' rad','Vel','Vel'
C     &,'pres','Vol'
c      do j=0,jj-2,2
C      do j=0,jj-2,1
c       write(*,'(I4,0I2,10f15.7)')
C       write(*,*)
C     &  j,  !ieos(1,j+1),ieos(2,j+1),
C     &  !ibc(j),ibc(j+1),ibc(j+2),
C     &  ibc(j),r(n,j),r(n,j+1),r(n,j+2),U(n,j),P(n,j+1)!,
C     &  !v(n,j+1),m(j+1),rho(j+1),Y(j+1)
C      enddo
c
c     STS write initial conditions as 0th time step
         qtotal = 0.
         mvtotal = 0.
         ketotal = 0.
         ietotal = 0.
         etotal = 0.
         do j = 0,jj-2,2
            if (ibc(j+1) .eq. 0) then
              qtotal  = qtotal  + q(n+2,j)
              mvtotal = mvtotal + m(j+1)*(  u(n+1,j) + u(n+1,j+2) )/2.d0
              ketotal = ketotal + 0.5*m(j+1)*
     &                  (0.5*((u(n+1,j) + u(n+1,j+2))))**2
              ietotal = ietotal + E(n+2,j+1)/rho(j+1)*m(j+1)
              etotal  = ietotal + ketotal
           endif
         enddo
      
         do j=0,jj-1,2
             rhoout= 0.
             if (ibc(j) .eq. 0) then
                rhoout = rho(j+1)/v(n+2,j+1)
             endif  
           write(33,'(4I6,26e17.8e3)')
     &           stepn, j, ibc(j), mat(j),
     &           t(n+2), r(n-1,j), r(n,j),r(n+1,j),r(n+2,j),
     &           U(n+1,j), V(n+2,j+1), rhoout, rho(j+1),    
     &           E(n+2,j+1), P(n+2,j+1), sigmar(n,j+1), s1(n+2,j+1),
     &           s2(n+2,j+1),q(n+1,j+1),dtminj(j),phi(n,j),beta(n,j),
     &           epsi1(n+1,j+1),epsi2(n+1,j+1),deltaZarr(j),
     &           alocal(j+1),rholocal(j+1),sigmao(n,j+1),Temp(n+1,j+1),
     &           etotal
         enddo ! this is j loop
c
c        do j=0,jj-1,1
c           write(34,'(4I6,26e17.8e3)')
c     &           stepn, j, ibc(j), mat(j),
c     &           t(n+2), r(n-1,j), r(n,j),r(n+1,j),r(n+2,j),
c     &           U(n+1,j), V(n+2,j+1), rhoout, rho(j+1),    
c     &           E(n+2,j+1), P(n+2,j+1), sigmar(n,j+1), s1(n+2,j+1),
c     &           s2(n+2,j+1),q(n+1,j+1),dtminj(j),phi(n,j),beta(n,j),
c     &           epsi1(n+1,j+1),epsi2(n+1,j+1),deltaZarr(j),
c     &           alocal(j+1),rholocal(j+1),sigmao(n,j+1),Temp(n+2,j+1),
c     &           etotal
c         enddo ! this is j loop

       stepn=1 ! advance for the first time step
       print *,'Start Main Loop,  goto tstop = ',tstop
c      read(*,*) foo
c
 98   format(2A3,11x,A9,A9,A9,A10,a8)
c99   format(2I3,e9.2,9e9.2,2f7.4,6e9.2)
c
c*****************************************************************************
c******************************** MAIN LOOP **********************************
c*****************************************************************************
c
       call cpu_time(startTime)
       ncount = -1
c       if (t(n) .le. tstop) then
 1      ncount = ncount + 2 ! when this use to count to an integer, ncount was the max integer
        n     =  1    ! this use to be the time step counter when all n data was stored
!                         n must remain n=1 because the time data is no longer stored
      !if( m(1) .ne. m(3) ) then
      ! print *,'Masses do not match',n,m(1),m(3)
       !read(*,*) foo ! keep going to see if cylindrical and spherical geometries will work
      !endif
c
cccccccccccccccc   Contact Check   ccccccccccccccccccccccc
c
c     calculate exact time of contact between 2 nodes to meet at one point
c     removes the extra node and then joins
c     every node is assigned a value; if zero then a central different node
c     ibc is -1 or 1; second order left stencil or second order right stencil -- could be a rigid wall or inflow or outflow
c     every node has a boundary condition
c     -2 and -3; +2 +3 are internal boundaries (planes or fractures)
c     CONTACT IS NOT IN WILKINS BOOK - this is Borg code
c     
      if (1 .eq. 0) then ! if you have a gap -- when contact need to eliminate a surface node and join the 2 surfaces together
      do j=0,jj,2
      if (ibc(j) .ne. 9) then
        delt     =  (t(n+1)-t(n-1))/2.d0
        r(n+2,j) =  r(n,j)+U(n-1,j)*delt*2.0d0
        r(n+1,j) = (r(n,j)+r(n+2,j))/2.d0
      endif
      enddo
c
      do j=2,jj-2,2
        if ( ibc(j) .ne. 9 .and. r(n+2,j) .le. r(n+2,j-2) .and.
     &       ibc(j-1) .eq. 9) then
          icontact = 2
          ibc(j) = -3  !this is a flag to let the contact stuff after the
                       ! momentum step know that
                       ! the delta t was recaclulated so that the two nodes
                       ! that are about to collide will touch perfectely after
                       ! this time step.
         print *,'Contact Eminent!',n,j
         print *,'Previous half Time step',t(n)-t(n-1)
         print *,'For time step 2 x Deltat',2.d0*delt
         Print *,' Estimated collision:'
         print *,'Void Node at ',r(n,j-2),' will move to ',
     &              r(n,j-2)+U(n-1,j-2)*2.d0*deltat
         print *,'Void node velocity',U(n-1,j-2)
         print *,'J Node at ',r(n,j),' will move to ',
     &              r(n,j)  +U(n-1,j  )*2.d0*delt
         print *,'J node velocity',U(n-1,j)

        jjv       = j-2
        jjj       = j
        ww        = U(n-1,jjj) - U(n-1,jjv)
        rr        = r(n  ,jjj) - r(n  ,jjv)
c
c this is the ibc(j)= 1 for the inside (v) side of the boundary
        sigmar(n,jjv-1) = (-(P(n,jjv-1)+q(n-1,jjv-1))+s1(n,jjv-1))
        sigmao(n,jjv-1) = (-(P(n,jjv-1)+q(n-1,jjv-1))+s2(n,jjv-1))
      phiv      = (rho(jjv-1)*((r(n,jjv)-r(n,jjv-2))/V(n,jjv-1)))/2.d0
        betav     = (sigmar(n,jjv-1)-sigmao(n,jjv-1))*V(n,jjj-1)/
     &              (r(n,jjv-1)*rho(jjv-1))
c this is the ibc(j)=-1 for the outside (j) side of the boundary
        sigmar(n,jjj+1) = (-(P(n,jjj+1)+q(n-1,jjj+1))+s1(n,jjj+1))
        sigmao(n,jjj+1) = (-(P(n,jjj+1)+q(n-1,jjj+1))+s2(n,jjj+1))
      phij      = (rho(jjj+1)*((r(n,jjj+2)-r(n,jjj))/V(n,jjj+1)))/2.d0
        betaj     = (sigmar(n,jjj+1)-sigmao(n,jjj+1))*V(n,jjj+1)/
     &              (r(n,jjj+1)*rho(jjj+1))
c
      aa = sigmar(n,jjj+1)/phij + sigmar(n,jjv-1)/phiv
     &   + (betav + betaj)*(d-1.d0)
      dthalf = t(n) - t(n-1)
      bb   = (2.d0*ww + aa * dthalf)
      dti0 = 0.d0
      dti1 = dti0-(aa*dti0*dti0+bb*dti0+2.d0*rr)/(2.d0*aa*dti0+bb)
      print *,'sig',sigmar(n,jjj+1),phij,sigmar(n,jjv-1),phiv
      print *,'bet',betav,betaj,d
      print *,aa,bb,rr,ww
      print *,rho(jjv-1),r(n,jjv),r(n,jjv-2),V(n,jjv-1)
      print *,rho(jjj+1),r(n,jjj+2),r(n,jjj),V(n,jjj+1)
      print *,'betav',sigmar(n,jjv-1),sigmao(n,jjv-1),
     &  r(n,jjj),r(n,jjj-2),V(n,jjj-1),rho(jjj-1)
c
      jj1=0
      do while (abs(dti1-dti0) .gt. 1.d-16)
      jj1 = jj1 + 1
      dti0 = dti1
      dti1 = dti0-(aa*dti0*dti0+bb*dti0+2.d0*rr)/(2.d0*aa*dti0+bb)
      print *,aa,bb,rr,ww
      print *,'iteration',dti0,dti1
      if (jj1 .gt.1) read(*,*) foo
      enddo
      print *,'Whole Time step adjusted from ',2.d0*delt
      deltat = dti1/2.d0                            !this might need to be adjusted because i changed the way time steping works
      print *,'To ',(deltat+dti1)
      print *,'where ',dti1,' is the time to impact'
      print *,'and ',deltat,' is an arbatrary small number.'
       read(*,*) foo
       endif !contact if
      enddo
       read(*,*) foo
      endif   !contact check if ! need this to do a 1d mesosecale
c
ccccccccccccccccccccccc  Advance Time step  cccccccccccccccccccccccccccccc
c
c     advance time step, for first calculation timestep is set by deltat_0
c     after that the time step is calculated after the end of the current time step, deltat
c
      if (icontact .eq. 2) then      ! if there was contact time step adjusted to perfectly
         t(n+1) = t(n) + dti1/2.d0   ! have nodes touch at the end of the time step
         t(n+2) = t(n+1) + deltat
         icontact = 0
         if  (dti1 - deltat .lt. 0. .or. t(n+2)-t(n+1) .lt. 0.) then
          print *,'Time step Error from contact'
          print *, dti1,deltat,dti1 - deltat , t(n+2)-t(n+1)
 10       read(*,*) foo
          go to 10
         endif
      else
         t(n+1) = t(n) + deltat/2.d0
         t(n+2) = t(n) + deltat
      endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     2.    Conservation of Momumtum  - start j
c
c  Boundary Conditions
c   -/+ 1 Free End, stress free
c   -/+ 2 and 3 Part of contact check
c   -/+ 4 Fixed end
c
c System is failing in this do loop line 450	
      do j=0,jj-1,2
      if (ibc(j+1) .eq. 0) then !this is true for central difference cells
       sigmar(n,j+1) = (-(P(n,j+1)+q(n-1,j+1))+s1(n,j+1))
       sigmao(n,j+1) = (-(P(n,j+1)+q(n-1,j+1))+s2(n,j+1))
      endif
      enddo
c	
      do j=0,jj,2
c
       delt = t(n+1)-t(n-1)    !this is different than deltat, deltat is the true timestep
      if (ibc(j) .ne. 9) then

       if (ibc(j) .eq. -1) then
        phi(n,j) = (rho(j+1)*((r(n,j+2)-r(n,j))/V(n,j+1)))/2.d0
        beta(n,j) = (sigmar(n,j+1)-sigmao(n,j+1))*V(n,j+1)
     &             /(r(n,j+1)*rho(j+1)) ! STSM this is just the value betewen j and j+2 as written in wilkins
        U(n+1,j) = U(n-1,j)+(delt/phi(n,j))
     &  *( sigmar(n,j+1)+bc(-1,1))+delt*beta(n,j)*(d-1.d0)
c       print *,'BC thing',bc(-1,1)
       elseif (ibc(j) .eq. 1) then
        phi(n,j) =  rho(j-1)*((r(n,j)-r(n,j-2))/V(n,j-1))/2.d0
        beta(n,j) = (sigmar(n,j-1)-sigmao(n,j-1))*V(n,j-1)
     &             /(r(n,j-1)*rho(j-1))
        U(n+1,j) = U(n-1,j)+(delt/phi(n,j))
     &  *(-bc( 1,1)-sigmar(n,j-1))+delt*beta(n,j)*(d-1.d0)
c       U(n+1,j) = 0.d0
c
c this is the new way which treats the voids like inner or outer bc's
       elseif (ibc(j) .eq. -2 .or. ibc(j) .eq. -3) then
        phi(n,j) = (rho(j+1)*((r(n,j+2)-r(n,j))/V(n,j+1)))/2.d0
        beta(n,j) = (sigmar(n,j+1)-sigmao(n,j+1))*V(n,j+1)
     &             /(r(n,j+1)*rho(j+1))
        U(n+1,j) = U(n-1,j)+(delt/phi(n,j))
     &  *(sigmar(n,j+1)+pvoid)+delt*beta(n,j)*(d-1.d0)
       elseif (ibc(j) .eq. 2) then
        phi(n,j) =  rho(j-1)*((r(n,j)-r(n,j-2))/V(n,j-1))/2.d0
        beta(n,j) = (sigmar(n,j-1)-sigmao(n,j-1))*V(n,j-1)
     &             /(r(n,j-1)*rho(j-1))
        U(n+1,j) = U(n-1,j)+(delt/phi(n,j))
     &  *(pvoid-sigmar(n,j-1))+delt*beta(n,j)*(d-1.d0)
       elseif (ibc(j) .eq. 0) then ! this is a central node and will follow equations B-2 for central nodes !! MAIN CENTRAL CELLS
        phi(n,j) = (0.5d0)*(rho(j+1)*((r(n,j+2)-r(n,j))/V(n,j+1))+
     &   rho(j-1)*((r(n,j)-r(n,j-2))/V(n,j-1)))
        beta(n,j)=( (sigmar(n,j+1)-sigmao(n,j+1))*(V(n,j+1)/rho(j+1))/
     &  (0.5d0*(r(n,j+2)+r(n,j  ))) +
     &            (sigmar(n,j-1)-sigmao(n,j-1))*(V(n,j-1)/rho(j-1))/
     &  (0.5d0*(r(n,j  )+r(n,j-2))) )/2.d0
        U(n+1,j) = U(n-1,j)+(delt/phi(n,j))
     &  *(sigmar(n,j+1)-sigmar(n,j-1))+delt*beta(n,j)*(d-1.d0)
       elseif (ibc(j) .eq. -4) then ! fixed boundary condition -4 or 4
        U(n+1,j) = 0.d0
c       if ( U(n+1,j) .lt. 1.d-5) U(n+1,j) = 0.d0
       elseif (ibc(j) .eq. 4) then
        U(n+1,j) = 0.d0
c       if ( U(n+1,j) .lt. 1.d-5) U(n+1,j) = 0.d0
       else
        print *,'Momentum ERROR!!!'
        print *,j,ibc(j)
        read(*,*) foo
       endif
c
      endif
      enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     3.  Position  - within j
c
      do j=0,jj,2
      if (ibc(j) .ne. 9) then ! 9 is reserved for extra nodes usually add 10 or 20 extra nodes; nodes at r=0 with 0 stress and pressure -- extra nodes used for spalling
        r(n+2,j) =  r(n,j)+U(n+1,j)*deltat
        r(n+1,j) = (r(n,j)+r(n+2,j))/2.d0
      endif
      enddo
c debug block
        if (idebug .eq. 1 .and. ncount .ge. ndebug) then
        j=jdebug
        print *,'j=',j,'n=',n,'ncount=',ncount
        print *,'time',t(n-1),t(n),t(n+1),t(n+2)
        print *,'deltat',deltat,delt
        print *,'vel',U(n,j),U(n+1,j)
        print *,'r',r(n,j-2),r(n,j),r(n,j+2)
        print *,'V',V(n,j+1),V(n,j),V(n,j-1)
        print *,'rho',rho(j-1),rho(j),rho(j+1)
        print *,'sigmar',sigmar(n,j+1),sigmar(n,j-1)
        print *,'sigmao',sigmao(n,j-1),sigmao(n,j+1)
        print *,'4',phi(n,j),beta(n,j),U(n-1,j)
        print *,'5',(deltat/phi(n,j))*(sigmar(n,j+1)-0.d0)
        read(*,*) foo
        endif
c	print *,'Line 537'

c
cccccc     Begin the connection of nodes due to contact check
c
      do j=2,jj-2,2 ! print info for the node contact
        if ( ibc(j) .eq. -3 ) then
          print *,'Joining Nodes',j-2,' and ',j,' at n= ',n
          print *,'Before Joining:'
          print *,'The two nodes (VOID (left)and J(right)) at:'
          print *,'left: r(void,n-2)=',r(n,j-4),' U= ',U(n-1,j-4)
          print *,'VOID: r(void,n-2)=',r(n,j-2),' U= ',U(n-1,j-2)
          print *,'J:    r(j   ,n-2)=',r(n,j)  ,' U= ',U(n-1,j)
          print *,'right:r(j   ,n-2)=',r(n,j+2),' U= ',U(n-1,j+2)
          print *,'then stepped to:'
          print *,'left: r(-2  )=',r(n+2,j-4),' U= ',U(n+1,j-4)
          print *,'VOID: r(void)=',r(n+2,j-2),' U= ',U(n+1,j-2)
          print *,'J:    r(j   )=',r(n+2,j)  ,' U= ',U(n+1,j)
          print *,'right:r(+2  )=',r(n+2,j+2),' U= ',U(n+1,j+2)
          print *,'where they were joined'
                   ibc(j)   = 0
       Usave = (m(j+1)*U(n+1,j)+m(j-3)*U(n+1,j-2))
     &                      /(m(j-3)+m(j+1))
       print *,'Usave',usave,m(j+1),U(n+1,j),m(j-3),U(n+1,j-2)
     &                      ,m(j-3),m(j+1)
          do jjj = j-2,jj-2 ! shifts all the nodes together
           do nz = n,n+2
              ibc(    jjj  )  =    ibc(    jjj+2)
              ibc(    jjj+1)  =    ibc(    jjj+3)
             ieos(1,  jjj+1)  =   ieos(1,  jjj+3)
             ieos(2,  jjj+1)  =   ieos(2,  jjj+3)
                m(    jjj  )  =      m(    jjj+2)
              rho(    jjj+1)  =    rho(    jjj+3)

                r(nz ,jjj  )  =     r(nz ,jjj+2)
                U(nz ,jjj  )  =     U(nz ,jjj+2)
              phi(nz ,jjj  ) =    phi(nz ,jjj+2)
             beta(nz ,jjj  ) =   beta(nz ,jjj+2)
           sigmar(nz ,jjj+1) = sigmar(nz ,jjj+3)
           sigmao(nz ,jjj+1) = sigmao(nz ,jjj+3)
                V(nz ,jjj+1) =      V(nz ,jjj+3)
               s1(nz ,jjj+1) =     s1(nz ,jjj+3)
               s2(nz ,jjj+1) =     s2(nz ,jjj+3)
               s3(nz ,jjj+1) =     s3(nz ,jjj+3)
                E(nz ,jjj+1) =      E(nz ,jjj+3)
                Y(    jjj+1) =      Y(    jjj+3)
            pfrac(    jjj  ) =  pfrac(    jjj+2)  !node value
            pfrac(    jjj+1) =  pfrac(    jjj+3)  !Cell value
           enddo
          enddo
                U(n+1 ,j  -2) = Usave
             pfrac(    j  -2) = 1.d-2
c                Y(     j  -1) = 0.d0   !ie the materials are not welded after impact
c               r(n+2 ,j  -2) = r(n+1,j-2)+Usave*(t(n+2)-t(n+1))
c          if (abs(r(0,ipoint+j1(jjj-1))-xstart(jjj)) .lt. 1.e-5) then
c           Print *,'Warning:  Joined moved to same spot'
c          endif
              print *,'After Joining:'
              print *,'left: r(j-2)=',r(n+2,j-4),' U= ',u(n+1,j-4)
              print *,'New:  r(j) = ',r(n+2,j-2),' U= ',usave
              print *,'right:r(j+2)=',r(n+2,j  ),' U= ',u(n+1,j)
               read(*,*) foo
      do jzz=1,jj-2,2
       write(*,'(I4,5I2,3f7.3,2e9.3,5f7.3)')
     &  jzz,ieos(1,jzz),ieos(2,jzz),
     &  ibc(jzz-1),ibc(jzz),ibc(jzz+1),
     &  r(n+2,jzz-1),r(n+2,jzz),r(n+2,jzz+2),
     &  U(n+1,jzz-1),U(n+1,jzz+2),U(n+1,jzz+1)
c     &  ,v(n+2,jzz+1),m(jzz+1),rho(jzz+1),y(jzz+1),pfrac(jzz+2)
      enddo
      read(*,*) foo
        endif
      enddo
c	print *,'line 608'
cccccccccccccccccccccc End of contact ccccccccccccccccccccccccccccc
c     Node passes it's neighbor check
c     check that no nodes are overlapping
c
      do j=2,jj-2,2
        if ( ibc(j) .ne. 9 .and. r(n+2,j) .lt. r(n+2,j-2) ) then
          print *,'Contact Check Failed',t(n),j
          print *,'Try lowering time step'
          print *,'Void',r(n+2,j-2),U(n+1,j-2)
          print *,'j   ',r(n+2,j)  ,U(n+1,j)
          print *,'deltat',deltat
          print *,bc(-1,1),bc(1,1)
          read(*,*) foo
       endif
      enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     3.  Relative VOLUME   - within j
c      B-3 conservation of mass
c      deltat     =  t(n+2)-t(n)
c      print *,'line 626'
c	read(*,*) foo
      do j=0,jj-2,2
      if (ibc(j+1) .eq. 0) then
c
        r(n+1,j+1) = ( r(n+1,j)+r(n+1,j+2) ) /2.d0 ! average data in the intermediate nodes - central; staggered mesh in time and space
        r(n+2,j+1) = ( r(n+2,j)+r(n+2,j+2) ) /2.d0
c
        V(n+2,j+1)=rho(j+1)*((r(n+2,j+2)**d-r(n+2,j)**d)/d)/m(j+1)
        V(n+1,j+1)=rho(j+1)*((r(n+1,j+2)**d-r(n+1,j)**d)/d)/m(j+1)
c
        if ( V(n+2,j+1) .eq. 0.) then
        print *,'------------- Zero volume error! #1-----------------'
        print *,'volume',n,j+1,r(n+2,j+2),r(n+2,j),V(n+2,j+1)
        read(*,*) foo
        endif
        if ( V(n+1,j+1) .eq. 0.) then
        print *,'------------- Zero volume error! #2-----------------'
        print *,'volume',n,j+1,r(n+1,j+2),r(n+1,j),V(n+1,j+1)
        read(*,*) foo
        endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     4. VELOCITY STRAINS  - finish j
c     half steps in time and space
         epsi1(n+1,j+1) = (U(n+1,j+2)-U(n+1,j))/(r(n+1,j+2)-r(n+1,j))
         epdt1(n+1,j+1) = (epsi1(n+1,j+1)-epsi1(n-1,j+1))/deltat !this is 1st order accurate ! stsm not used here
        if (d .eq. 1.) then
         epsi2(n+1,j+1) = 0.d0
        else
         epsi2(n+1,j+1) = (U(n+1,j+2)+U(n+1,j))/(r(n+1,j+2)+r(n+1,j))
        endif
      endif

      enddo  ! end j loop for volume
c--------------------  AMR Feature  -------------------------------

        if (idebug .eq. 1 .and. n .ge. ndebug) then
        j=jdebug
        print *,'r',j,r(n+2,j+1),m(j+1)
        print *,'epsil',j,V(n+1,j+1),V(n+2,j+1),
     &                    epsi1(n+1,j+1),epsi2(n+1,j+1)
        read(*,*) foo
        endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     5a. STRESSES  - start j
c     B-5 stress deviators and then pressure done below 5b
c      print *,'line 679'
      do j=0,jj-2,2
      if (ibc(j+1) .eq. 0) then
c       deltat = t(n+2)-t(n)
       if (iEOS(2,j+1) .ne. 3) then  
        s1(n+2,j+1) = s1(n,j+1) + 2.d0*eos(ieos(1,j+1),6)
     &             * ( epsi1(n+1,j+1)*deltat
     &                 - (V(n+2,j+1)-V(n,j+1))/(3.d0*V(n+1,j+1)) )
        s2(n+2,j+1) = s2(n,j+1) + 2.d0*eos(ieos(1,j+1),6)
     &             * ( epsi2(n+1,j+1)*deltat
     &                 - (V(n+2,j+1)-V(n,j+1))/(3.d0*V(n+1,j+1)) )
        s3(n+2,j+1) = -(s1(n+2,j+1)+s2(n+2,j+1))
c
        s1(n+1,j+1) = ( s1(n+2,j+1) + s1(n,j+1) )/2.d0 ! fill in the staggered mesh averages
        s2(n+1,j+1) = ( s2(n+2,j+1) + s2(n,j+1) )/2.d0
        s3(n+1,j+1) = ( s3(n+2,j+1) + s3(n,j+1) )/2.d0
       else     !Gamma law ideal gas Newtonian Stress
        s1(n+2,j+1) = 4.d0*1.8d-10*epsi1(n+1,j+1)/3.d0
     &             -1.387e-6* (V(n+2,j+1)-V(n,j+1))/V(n+1,j+1) 
c        s1(n+2,j+1) = s1(n,j+1)+4.d0*1.8d-10*epdt1(n+1,j+1)*deltat/3.d0
c     &             -1.387e-6* (V(n+2,j+1)-V(n,j+1))/(3.d0*V(n+1,j+1)) 
c
c        print *,epdt1(n+1,j+1),epsi1(n+1,j+1),epsi1(n-1,j+1)
c
        s1(n+1,j+1) = ( s1(n+2,j+1) + s1(n,j+1) )/2.d0
      endif ! iEOS check
c
      endif
      enddo
        if (idebug .eq. 1 .and. n .ge. ndebug) then
      j=jdebug
      j=199
      print *,'1',j,eos(ieos(1,j+1),6),V(n+2,j+1),V(n,j+1),V(n+1,j+1)
      print *,'s',s1(n+2,j+1),s2(n+2,j+1),s3(n+2,j+1),
     & epsi1(n+1,j+1),deltat
        endif
c      if (idebug .eq. 1 .and. n .ge. ndebug) read(*,*) foo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     6. VON MISES YIELD CONDITION - start j
c
      do j=0,jj-2,2
      if (ibc(j+1) .eq. 0) then
       if (iEOS(2,j+1) .ne. 3) then  
c calculate the deviatoric strain at n+1 and compare it to the yield strength
c 
       k(n+1,j+1) =(s1(n+1,j+1)**2 + s2(n+1,j+1)**2 + s3(n+1,j+1)**2)
     &           -(2.d0/3.d0)*Y(j+1)**2
c       print *,'Yield check',Y(j+1),k(n+2,j+1)
       if (k(n+1,j+1) .gt. 0.) then
c       Print *,'material yielded',j,k(n+2,j+1)
c        read(*,*) foo
c        print *,s1(n+2,j+1),s2(n+2,j+1),s3(n+2,j+1)
         xx = sqrt(s1(n+1,j+1)**2 + s2(n+1,j+1)**2 + s3(n+1,j+1)**2)
         xx = sqrt(2.d0/3.d0)*Y(j+1)/xx
         s1(n+1,j+1) =   xx*s1(n+1,j+1)
         s2(n+1,j+1) =   xx*s2(n+1,j+1)
         s3(n+1,j+1) =   xx*s3(n+1,j+1)
       endif
c calculate the deviatoric strain at n+2 and compare it to the yield strength

         k(n+2,j+1) =(s1(n+2,j+1)**2 + s2(n+2,j+1)**2 + s3(n+2,j+1)**2)
     &           -(2.d0/3.d0)*Y(j+1)**2
c       print *,'Yield check',Y(j+1),k(n+2,j+1)
       if (k(n+2,j+1) .gt. 0.) then
c       Print *,'material yielded',j,k(n+2,j+1)
c        read(*,*) foo
c        print *,s1(n+2,j+1),s2(n+2,j+1),s3(n+2,j+1)
         xx = sqrt(s1(n+2,j+1)**2 + s2(n+2,j+1)**2 + s3(n+2,j+1)**2)
         xx = sqrt(2.d0/3.d0)*Y(j+1)/xx
         s1(n+2,j+1) =   xx*s1(n+2,j+1)
         s2(n+2,j+1) =   xx*s2(n+2,j+1)
         s3(n+2,j+1) =   xx*s3(n+2,j+1)
c         s1(n+1,j+1) = ( s1(n+2,j+1) + s1(n,j+1) )/2.d0
c         s2(n+1,j+1) = ( s2(n+2,j+1) + s2(n,j+1) )/2.d0
c         s3(n+1,j+1) = ( s3(n+2,j+1) + s3(n,j+1) )/2.d0
c        Print *,'material yielded',j,k(n+2,j+1),Y(j+1)
        else !iEOS = 3, i.e. gamma law gas
        endif
c
c        print *,s1(n+2,j+1),s2(n+2,j+1),s3(n+2,j+1)
         if(idebug .eq. 1 ) then
          Print *,'Yield Check',Y(j+1),k(n+2,j+1)
          read(*,*) foo
         endif ! debug if
       endif  !ibc(j+1) .eq. 0
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     7. ARTIFICIAL VISCOSITY - within j
c
         if (1 .eq. 1) then
c           print *,j,u(n+1,j+2),u(n+1,j)
            if( u(n+1,j+2) .lt. u(n+1,j) .and.
     &          V(n+2,j+1)-V(n,j+1) .lt. 0. ) then
              xx = rho(j+1)/V(n+1,j+1)
              if (p(n,j+1) .lt. 0.) then
c                print *,'Negative Pressure!'
                 a  = dsqrt(-P(n,j+1)/xx)
              else
                 a  = dsqrt( P(n,j+1)/xx)
              endif
              if (dabs(P(n,j+1)) .lt. 0.002d0) a = 0.d0
              rholocal(j+1) = xx
              alocal(j+1) = a
c             print *,'Sound Speed',j,P(n,j+1),xx,a
              q(n+1,j+1) = Co*Co*xx*    (U(n+1,j+2)-U(n+1,j))**2 
     &               + CL*a *xx*dabs(U(n+1,j+2)-U(n+1,j))
            else
              q(n+1,j+1)=0.d0
            endif
         endif ! if vel strain and volume check is on
         endif ! if 1=1
      enddo !j counter
        if (idebug .eq. 1 .and. n .ge. ndebug) then
         j=jdebug
       print *,'Q',j,q(n+1,j+1),xx,(U(n+1,j+2)-U(n+1,j)),
     &    a,rho(j+1),v(n,j+1),P(n,j+1)
        endif
             if(idebug .eq.1 .and. n .ge. ndebug) read(*,*) foo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     8.  ENERGY  - start j
c      EOS provides the internal energy at V and P
c     iEOS(2,imat)=1 - Mie Grunisen
c     iEOS(2,imat)=2 - Gamma law ideal gas
c     iEOS(2,imat)=3 - Gamma law ideal gas with Newtonian Stress and heat transfer
c     iEOS(2,imat)=4 - snow plow model with KO inputs for Hugoniot
c     iEOS(2,imat)=5 - snow plow model with anamolous hugoniot
c     iEOS(2,imat)=6 - P-alpha Model
c
      do j=0,jj-2,2
      if (ibc(j+1) .eq. 0) then
c
       if (iEOS(2,j+1) .eq. 1) then                    !Mie Grunisen
        gamma0 = eos(ieos(1,j+1),4)
        k1  = rho(j+1)*eos(ieos(1,j+1),1)**2
        k2  =  k1*(2.d0*eos(ieos(1,j+1),2)-gamma0/2.d0)
        k3  =  k1*(3.d0*eos(ieos(1,j+1),2)-gamma0)*eos(ieos(1,j+1),2)
        qbar = (q(n+1,j+1)+q(n-1,j+1))/2.d0
        deltaZ = V(n+1,j+1)*(s1(n+1,j+1)*epsi1(n+1,j+1)
     &       +(d-1.d0)*s2(n+1,j+1)*epsi2(n+1,j+1))*deltat
        deltaZarr(j) = deltaZ
        strain = 1.d0 - V(n+2,j+1)
        xb = gamma0 !gamma
        if(strain .lt. 0.d0) then
         xa = k1*strain + k3*strain**3
        else
         xa = k1*strain + k2*strain**2 + k3*strain**3
        endif
        E(n+2,j+1)=
     &  ( E(n,j+1)-((xa+P(n,j+1))/2.d0+ qbar)*(V(n+2,j+1)-V(n,j+1))
     &            + deltaZ )
     &   / (1.d0 + xb*(V(n+2,j+1)-V(n,j+1))/2.d0)
        E(n+1,j+1)= ( E(n+2,j+1) + E(n,j+1) )/2.d0
c      if (n .eq. 33) then
c        write (*,'(a2,I3,9e9.2)')
c     &   'E ',j,E(n+2,j+1),E(n+2,j+1),xa,P(n,j+1),qbar,deltaZ
c     &   ,V(n+2,j+1),V(n,j+1)
cc      read(*,*) foo
c      endif
c
       elseif (iEOS(2,j+1) .eq. 2) then                !Gamma law ideal gas
        qbar   = (q(n+1,j+1)+q(n-1,j+1))/2.d0
        delt   = t(n+1)-t(n)
        deltaZ = V(n+1,j+1)*(s1(n+1,j+1)*epsi1(n+1,j+1)
     &              +(d-1.d0)*s2(n+1,j+1)*epsi2(n+1,j+1))*delt
c
        xa   = 0.d0
        xb   = (eos(ieos(1,j+1),4)-1.d0)/V(n+2,j+1)  !gamma
        E(n+2,j+1)= ( E(n,j+1)
     &   -((xa+P(n,j+1))/2.d0+ qbar)*(V(n+2,j+1)-V(n,j+1)) + deltaZ )/
     &  (1.d0+xb*(V(n+2,j+1)-V(n,j+1))/2.d0)
        E(n+1,j+1)= ( E(n+2,j+1) + E(n,j+1) )/2.d0
c      if (j .eq. 0 ) then
c        if (n .lt.2. .or. n .eq. 201)    E(n+2,j+1)=0.d0
c        xx=V(n+2,j+1)-V(n,j+1)
c        print *,E(n,j+1)   /(1.*(1.+xb*(V(n+2,j+1)-V(n,j+1))/2.))
c        print *,xa*xx      /(2.*(1.+xb*(V(n+2,j+1)-V(n,j+1))/2.))
c        print *,P(n,j+1)*xx/(2.*(1.+xb*(V(n+2,j+1)-V(n,j+1))/2.))
c        print *,qbar*xx    /(1.*(1.+xb*(V(n+2,j+1)-V(n,j+1))/2.))
c        print *,deltaZ  /(1.*(1.+xb*(V(n+2,j+1)-V(n,j+1))/2.))
c      endif
c
       elseif (iEOS(2,j+1) .eq. 3) then                !Gamma law ideal gas     ???
        qbar = (q(n+1,j+1)+q(n-1,j+1))/2.d0
        delt   = t(n+1)-t(n)
        deltaZ = V(n+1,j+1)*(s1(n+1,j+1)*epsi1(n+1,j+1)
     &              +(d-1.0d0)*s2(n+1,j+1)*epsi2(n+1,j+1))*delt
c
        xx   = rho(j+1)/V(n+2,j+1)
        xa   = 0.d0
        xb   = (eos(ieos(1,j+1),4)-1.d0)/V(n+2,j+1)
        E(n+2,j+1)= ( E(n,j+1)
     &   -((xa+P(n,j+1))/2.d0+ qbar)*(V(n+2,j+1)-V(n,j+1))+ deltaZ )/
     &  (1.d0+xb*(V(n+2,j+1)-V(n,j+1))/2.d0)
c
       elseif (iEOS(2,j+1) .eq. 4) then                !Snow Plow
        qbar = (q(n+1,j+1)+q(n-1,j+1))/2.d0
        deltaZ = V(n+1,j+1)*(s1(n+1,j+1)*epsi1(n+1,j+1)
     &             +(d-1.d0)*s2(n+1,j+1)*epsi2(n+1,j+1))*deltat
        xb = eos(ieos(1,j+1),4) !gamma
        v0 = eos(ieos(1,j+1),12)
        xx = 1.d0 - (V(n+2,j+1)/rho(j+1))/V0   !this makes the compression relative to the compacted material hugoniot
        if ( v(n+2,j+1)/rho(j+1) .le. v0) icompact(n+2,j+1) = 1
        if(     icompact(n+2,j+1) .eq. 1 .and. xx .lt. 0.d0) then
         xa = eos(ieos(1,j+1),1)*xx+eos(ieos(1,j+1),3)*xx**3
        elseif (icompact(n+2,j+1) .eq. 1 .and. xx .ge. 0.d0) then
         xa = eos(ieos(1,j+1),1)*xx+eos(ieos(1,j+1),2)*xx**2
     &       +eos(ieos(1,j+1),3)*xx**3
        elseif( icompact(n+2,j+1) .eq. 0) then
         xa = 0.d0   !snow plow model - this zeros out pressure until threshold is reached
        else
         print *,'compaction error in energy'
         print *,icompact(n+2,j+1),xx
         read(*,*) foo
        endif
        E(n+2,j+1)=
     &  ( E(n,j+1)-((xa+P(n,j+1))/2.d0+ qbar)*(V(n+2,j+1)-V(n,j+1))
     &            + deltaZ )
     &   / (1.d0 + xb*(V(n+2,j+1)-V(n,j+1))/2.d0)
        E(n+1,j+1)= ( E(n+2,j+1) + E(n,j+1) )/2.d0
c
       elseif (iEOS(2,j+1) .eq. 5) then                !Snow Plow with anaomolus hugoniot
        qbar = (q(n+1,j+1)+q(n-1,j+1))/2.d0
        deltaZ = V(n+1,j+1)*(s1(n+1,j+1)*epsi1(n+1,j+1)
     &             +(d-1.d0)*s2(n+1,j+1)*epsi2(n+1,j+1))*deltat
        stemp = 0.1048d0      !slope
        ctemp = 0.5124d0      !bulk sound speed
        gtemp = 0.9d0         !gamma
        v0 = eos(ieos(1,j+1),8)
        V00= 1.d0/rho(j+1)
        vv = V(n+2,j+1)/rho(j+1)    !v_local
        if ( v(n+2,j+1)/rho(j+1) .le. v0) icompact(n+2,j+1) = 1
        if (icompact(n+2,j+1) .eq. 1) then  !icompact is the flag, once compact always compact
         xa = ((2.d0*vv-gtemp*(v0 -vv))*ctemp*ctemp*(v0-vv))     !porous hugoniot, see meyers pg 141
     &       /((2.d0*vv-gtemp*(v00-vv))*((v0-stemp*(v0-vv))**2) )
         xa = dabs(xa) !not sure if i need this but it was add so that re-denstended materials don't have negative pressure
        else
         xa = 0.d0
        endif
        E(n+2,j+1)=
     &  ( E(n,j+1)-((xa+P(n,j+1))/2.d0+ qbar)*(V(n+2,j+1)-V(n,j+1))
     &            + deltaZ )
     &   / (1.d0 + xb*(V(n+2,j+1)-V(n,j+1))/2.d0)
        E(n+1,j+1)= ( E(n+2,j+1) + E(n,j+1) )/2.d0
c
       elseif (iEOS(2,j+1) .eq. 6) then                 !p-alpha
        qbar = (q(n+1,j+1)+q(n-1,j+1))/2.d0
        deltaZ = V(n+1,j+1)*(s1(n+1,j+1)*epsi1(n+1,j+1)
     &             +(d-1.d0)*s2(n+1,j+1)*epsi2(n+1,j+1))*deltat
        xx = 1.d0 - V(n+2,j+1)
        xb = eos(ieos(1,j+1),4) !gamma
        if(xx .lt. 0.) then
        xa = eos(ieos(1,j+1),1)*xx+0.d0*xx**2+eos(ieos(1,j+1),3)*xx**3
        else
        xa = eos(ieos(1,j+1),1)*xx+eos(ieos(1,j+1),2)*xx**2
     &      +eos(ieos(1,j+1),3)*xx**3
        endif
        E(n+2,j+1)=
     &  ( E(n,j+1)-((xa+P(n,j+1))/2.d0+ qbar)*(V(n+2,j+1)-V(n,j+1))
     &            + deltaZ )
     &   / (1.d0 + xb*(V(n+2,j+1)-V(n,j+1))/2.d0)
        E(n+1,j+1)= ( E(n+2,j+1) + E(n,j+1) )/2.d0
c
c
       else
        print *,'EOS error!'
        print *,n,j,ieos(2,j+1)
        read(*,*) foo
       endif
      endif
c
      enddo
c
      if(idebug .eq.1  .and. n .ge. ndebug) then
        j=jdebug
        print *,'Energy',j,E(n+2,j+1),E(n,j+1),xx,
     & (2.d0*(1.d0+xb*(V(n+2,j+1)-V(n,j+1))/2.d0))
      endif
      if(idebug .eq.1  .and. n .gt. ndebug) read(*,*) foo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     5b PRESSURE - start j MOVED DOWN BECAUSE IT NEEDS ENERGY
c     DEPENDS ON EOS
c
c      print *,'line 935'
      do j = 0,jj-2,2
      if (ibc(j+1) .eq. 0) then
       if (iEOS(2,j+1) .eq. 1) then   !     Mie Gruneisen EOS as formulated in Wilkins
        gamma0 = eos(ieos(1,j+1),4)
        k1  = rho(j+1)*eos(ieos(1,j+1),1)**2
        k2  =  k1*(2.d0*eos(ieos(1,j+1),2)-gamma0/2.d0)
        k3  =  k1*(3.d0*eos(ieos(1,j+1),2)-gamma0)*eos(ieos(1,j+1),2)
c        if (j .eq. 2840) then
c        print *,eos(ieos(1,j+1),1),eos(ieos(1,j+1),2),eos(ieos(1,j+1),3)
c        read(*,*) foo
c        endif
        strain      = 1.d0 - V(n+1,j+1)
        if (strain .lt. 0.d0) then
         P(n+1,j+1) = k1*strain + k3*strain**3 + gamma0*E(n+1,j+1)
        else
         P(n+1,j+1) = k1*strain + k2*strain**2 + k3*strain**3
     &                                         + gamma0*E(n+1,j+1)
        endif
        strain      = 1.d0 - V(n+2,j+1)
        if (strain .lt. 0.d0) then
         P(n+2,j+1) = k1*strain + k3*strain**3 + gamma0*E(n+2,j+1)
        else
         P(n+2,j+1) = k1*strain + k2*strain**2 + k3*strain**3
     &                                         + gamma0*E(n+2,j+1)
        endif
c        print *,'Mie Grun EOS 1',P(n+1,j+1),P(n+2,j+1)
c        read(*,*) foo
c       if (icontact .eq. 1) then
c       if (P(n+2,j+1) .lt. 0.) then
c       print *,'Neg Pres',n+1,j+1,xx,v(n+1,j+1),E(n+1,j+1),p(n+1,j+1)
c       print *,'Neg Pres',n+2,j+1,xx,v(n+2,j+1),E(n+2,j+1),p(n+2,j+1)
cc       read(*,*) foo
c       endif
c       endif
       elseif (iEOS(2,j+1) .eq. 2) then            ! Gamma Law (perfect gas)
         P(n+2,j+1) = (eos(ieos(1,j+1),4)-1.d0)*E(n+2,j+1)/V(n+2,j+1)
         P(n+1,j+1) = (P(n,j+1)+P(n+2,j+1))/2.d0
       elseif (ieos(2,j+1) .eq. 3) then            ! Gamma Law
         P(n+2,j+1) = (eos(ieos(1,j+1),4)-1.d0)*E(n+2,j+1)/V(n+2,j+1)
         P(n+1,j+1) = (P(n,j+1)+P(n+2,j+1))/2.d0
c        xa          = (u(n+2,j)+u(n+2,j+2))/2.d0
c        xx          = rho(j+1)/V(n+2,j+1)
c         P(n+2,j+1) = P(n,j+1)
c     &     +(eos(ieos(1,j+1),4)+1.d0)*xx*xa*xa/2.d0
c       print *,'pres',rho_0,gamma+1,xa,xx
       elseif (ieos(2,j+1) .eq. 4) then             !Snow Plow
        v0  = eos(ieos(1,j+1),12)
        V00 = 1.d0/rho(j+1)
        vv  = v(n+1,j+1)/rho(j+1)
        if (v0 .eq. 0.) then
         print *,'v0 can not equal zero'
         print *,j,ieos(1,j+1),eos(ieos(1,j+1),12)
         read(*,*) foo
        endif
c
        if (vv .le. V0) icompact(n+1,j+1) = 1
         xx         = 1.d0 - (V(n+1,j+1)/rho(j+1))/V0   !this makes the compression relative to the compacted material hugoniot
        if (icompact(n+1,j+1) .eq. 1 .and. xx .lt. 0.) then
          P(n+1,j+1) = eos(ieos(1,j+1),1)*xx+
     &                 eos(ieos(1,j+1),3)*xx**3+
     &                 eos(ieos(1,j+1),4)*E(n+1,j+1)
        elseif(icompact(n+1,j+1) .eq. 1 .and. xx .ge. 0.d0) then
          P(n+1,j+1) = eos(ieos(1,j+1),1)*xx+
     &                 eos(ieos(1,j+1),2)*xx**2+
     &                 eos(ieos(1,j+1),3)*xx**3+
     &                 eos(ieos(1,j+1),4)*E(n+1,j+1)
        elseif (icompact(n+1,j+1) .eq. 0) then
          P(n+1,j+1) = 0.d0*eos(ieos(1,j+1),4)*E(n+1,j+1)
        else
         Print *,'pressure compaction error'
         read(*,*) foo
        endif
c        if ( icompact(n+1,j+1) .eq. 1) then
c        print *,t(n),j+1,icompact(n+1,j+1),vv,p(n+1,j+1)
c        print *,eos(ieos(1,j+1),1)*xx,eos(ieos(1,j+1),3)*xx**3
c     &         ,eos(ieos(1,j+1),4)*E(n+1,j+1)
c        print *,xx,eos(ieos(1,j+1),1),eos(ieos(1,j+1),3)
c     &         ,eos(ieos(1,j+1),4)
c        read(*,*) foo
c        else
c        print *,t(n),j+1,icompact(n+1,j+1),vv,p(n+1,j+1)
c        endif
c
        if (v(n+2,j+1)/rho(j+1) .le. V0) icompact(n+2,j+1) = 1
         xx         = 1.d0 - (V(n+2,j+1)/rho(j+1))/V0     !this makes the compression relative to the compacted material hugoniot
         if (icompact(n+2,j+1) .eq. 1 .and. xx .lt. 0.) then
          P(n+2,j+1) = eos(ieos(1,j+1),1)*xx+
     &                 eos(ieos(1,j+1),3)*xx**3+
     &                 eos(ieos(1,j+1),4)*E(n+2,j+1)
         elseif (icompact(n+2,j+1) .eq. 1 .and. xx .ge. 0.) then
          P(n+2,j+1) = eos(ieos(1,j+1),1)*xx+
     &                 eos(ieos(1,j+1),2)*xx**2+
     &                 eos(ieos(1,j+1),3)*xx**3+
     &                 eos(ieos(1,j+1),4)*E(n+2,j+1)
        elseif (icompact(n+2,j+1) .eq. 0) then  !v .gt. v0
          P(n+2,j+1) = 0.d0*eos(ieos(1,j+1),4)*E(n+2,j+1)
        else
         Print *,'pressure compaction error'
         read(*,*) foo
        endif
c
       elseif (ieos(2,j+1) .eq. 5) then   ! Snow Plow with anaomolous Hugoniot
         V0 = eos(ieos(1,j+1),8)
         V00= 1.d0/rho(j+1)

         vv = V(n+1,j+1)/rho(j+1)    !v_local
         if (VV .le. V0) icompact(n+1,j+1) = 1
         if (icompact(n+1,j+1) .eq. 1) then
          xa = ((2.d0*vv-gtemp*(v0 -vv))*ctemp*ctemp*(v0-vv))    !porous hugoniot, see meyers pg 141
     &        /((2.d0*vv-gtemp*(v00-vv))*((v0-stemp*(v0-vv))**2) )
          if (xa .lt. 0.d0) then
           print *,'Compression less than v0 on anamolous hugoniot, n+1'
           xa=dabs(xa)
          endif
         else
          xa = 0.d0
         endif

         P(n+1,j+1) = xa + eos(ieos(1,j+1),4)*E(n+1,j+1)

         vv = V(n+2,j+1)/rho(j+1)    !v_local
         if (VV .le. V0) icompact(n+2,j+1) = 1
         if (icompact(n+2,j+1) .eq. 1) then
         xa = ((2.d0*vv-gtemp*(v0 -vv))*ctemp*ctemp*(v0-vv))     !porous hugoniot, see meyers pg 141
     &       /((2.d0*vv-gtemp*(v00-vv))*((v0-stemp*(v0-vv))**2) )
          if (xa .lt. 0.d0) then
           print *,'Compression less than v0 on anamolous hugoniot, n+2'
           xa=dabs(xa)
          endif
         else
          xa = 0.d0
         endif
c
       elseif (ieos(2,j+1) .eq. 6) then  !p-alpha model
c
c
       elseif  (ieos(2,j+1) .eq. 7) then  ! Mie Gruneisen CTH like formulation
        P0 = 0.d0
        E0 = 0.d0
        T0 = 293.d0
        strain        = 1.d0 - V(n+1,j+1)
        rho_local = rho(j+1)/V(n+1,j+1)
        up = (U(n+1,j)+U(n+1,j+2))/2.d0
        Us = eos(ieos(1,j+1),1)           ! local shock speed
     &     + eos(ieos(1,j+1),2)*up
     &     + eos(ieos(1,j+1),3)/eos(ieos(1,j+1),1)*up**2
        PH = P0 + rho_local*Us*up
        EH = E0 + up*up/2.d0
        TH = dexp( eos(ieos(1,j+1),4)*strain)*T0
c       E  = EH + eos(ieos(1,j+1),8)*(T-TH)
c       P(n+2,j+1) = PH + eos(ieos(1,j+1),4)*rho_local*(
c         if (Us*up .ne. Us*Us*strain) then
c          print *,'EOS does not jive'
c          read(*,*) foo
c         endif
        xx         = 1.d0 - V(n+2,j+1)
        P(n+2,j+1) = eos(ieos(1,j+1),1)*xx+
     &               eos(ieos(1,j+1),2)*xx**2+
     &               eos(ieos(1,j+1),3)*xx**3+
     &               eos(ieos(1,j+1),4)*E(n+2,j+1)
       else    !this is ieos(2,j+1) else
      Print *,'eos error in setup file'
        read(*,*) foo
       endif  ! IEOS
      endif !ibc
      enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c This is the Temperature and Entropy Calculation
c      print *,'line 1104'
      do j=0,jj-2,2
c	print *,'line 1106',j,nc
c	print *,'temp',temp(n+2c,j+1)
c	print *,'E   ',e(n+2,j+1)
c	print *,'rho ',rho(j+1)
c	print *,'ieos',ieos(1,j+1)
c	print *,'eos ',eos(ieos(1,j+1),8)
c
c       Temp(n+2,j+1)=E(n+2,j+1)/(rho(j+1)*eos(ieos(1,j+1),8)) ! heat capacity energy change
       Temp(n+2,j+1)=E(n+2,j+1)/eos(ieos(1,j+1),8) ! heat capacity energy change ! STS fixed 2/5/23
       Temp(n+1,j+1)=0.5*(Temp(n,j+1)+Temp(n+2,j+1))
c       print *,Temp(n+2,j+1),E(n+2,j+1),rho(j+1),v(n+2,j+1)
c     &  ,eos(ieos(1,j+1),8)
c      print *,'line 1109',j,n
       entropy(n+2,j+1)=entropy(n,j+1)+
     &  s1(n+2,j+1)/Temp(n+2,j+1)  ! entropy from deviatoric strain
c       entropy(n+2,j+1)=entropy(n,j+1)+
c       print *,entropy(n+2,j+1),entropy(n,j+1),
c     &  s1(n+2,j+1),epsi1(n+2,j+1),Temp(n+2,j+1)
      enddo
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c This is the Pressure and Energy convergence check
c     print *,'line 1118'
      do j=0,jj-2,2
       if (iEOS(2,j+1) .eq. 1) then                    !Mie Grunisen
        gamma0 = eos(ieos(1,j+1),4)
        k1  = rho(j+1)*eos(ieos(1,j+1),1)**2
        k2  =  k1*(2.d0*eos(ieos(1,j+1),2)-gamma0/2.d0)
        k3  =  k1*(3.d0*eos(ieos(1,j+1),2)-gamma0)*eos(ieos(1,j+1),2)
        qbar = (q(n+1,j+1)+q(n-1,j+1))/2.d0
        deltaZ = V(n+1,j+1)*(s1(n+1,j+1)*epsi1(n+1,j+1)
     &             +(d-1.d0)*s2(n+1,j+1)*epsi2(n+1,j+1))*deltat
        strain = 1.d0 - V(n+2,j+1)
        xb = gamma0 !gamma
        if(strain .lt. 0.d0) then
        xa = k1*strain + k3*strain**3
        else
        xa = k1*strain + k2*strain**2 + k3*strain**3
        endif
        En2j1=
     &  ( E(n,j+1)-((xa+P(n,j+1))/2.d0+ qbar)*(V(n+2,j+1)-V(n,j+1))
     &            + deltaZ )
     &   / (1.d0 + xb*(V(n+2,j+1)-V(n,j+1))/2.d0)
         diffE=    E(n+2,j+1)-En2j1
        if (dabs(diffE) .gt. 0.d-6 ) then
         print *,'Energy not converged',En2j1,E(n+2,j+1)
c         if (
        endif
       elseif (  iEOS(2,j+1) .eq. 0) then
       !this is for not used nodes
       else
c        print *,'Non-MG EOS',j+1,iEOS(2,j+1)
c        read(*,*) foo
       endif
      enddo
cccccccccccccccccccccccccc   Spall   ccccccccccccccccccccccccccccccccccccccccccc
c     Check for tensile fracture, ie does the pressure/stress exceed pfrac at j?
c     check for physical separation
c     central node split - at interface node and set boundary conditions and shift variables
c
      if (1 .eq. 0) then
      do j=0,jj-2,2
       if (-p(n+2,j+1)+s1(n+2,j+1) .gt. pfrac(j)) then
        print *,j,n,-p(n+2,j+1)+s1(n+2,j+1),pfrac(j)
        print *,"Tensile Fracture!!"
        print *,'Separating Node',j,' at n= ',n
        print *,'Before Separating:'
        print *,'The node at:'
        print *,'left: r(j-2 ,n+2)=',r(n+2,j-2),' U(n-1)=',U(n+1,j-2)
        print *,'J:    r(j   ,n+2)=',r(n+2,j)  ,' U(n-1)=',U(n+1,j+0)
        print *,'right:r(j+2 ,n+2)=',r(n+2,j+2),' U(n-1)=',U(n+1,j+2)
        print *,'      r(j+4 ,n+2)=',r(n+2,j+4),' U(n-1)=',U(n+1,j+4)
       print *,'Frac',-p(n+2,j+1)+s1(n+2,j+1),pfrac(j)
        do jjj = jj-4,j+2,-2
         do nz = n,n+2
              ibc(    jjj-0)  =    ibc(    jjj-2)
              ibc(    jjj-1)  =    ibc(    jjj-3)
             ieos(1,  jjj-1)  =   ieos(1,  jjj-3)
             ieos(2,  jjj-1)  =   ieos(2,  jjj-3)
             ieos(1,  jjj-0)  =   ieos(1,  jjj-2)
             ieos(2,  jjj-0)  =   ieos(2,  jjj-2)
                m(    jjj-1)  =      m(    jjj-3)
              rho(    jjj-1)  =    rho(    jjj-3)

                r(nz ,jjj-0)  =     r(nz ,jjj-2)
                r(nz ,jjj-1)  =     r(nz ,jjj-3)
                U(nz ,jjj-0)  =     U(nz ,jjj-2)
              phi(nz ,jjj-0) =    phi(nz ,jjj-2)
             beta(nz ,jjj-0) =   beta(nz ,jjj-2)
           sigmar(nz ,jjj-1) = sigmar(nz ,jjj-3)
           sigmao(nz ,jjj-1) = sigmao(nz ,jjj-3)
                V(nz ,jjj-1) =      V(nz ,jjj-3)
               s1(nz ,jjj-1) =     s1(nz ,jjj-3)
               s2(nz ,jjj-1) =     s2(nz ,jjj-3)
               s3(nz ,jjj-1) =     s3(nz ,jjj-3)
                q(nz ,jjj-1) =      q(nz ,jjj-3)
            epsi1(nz ,jjj-1) =  epsi1(nz ,jjj-3)
            epsi2(nz ,jjj-1) =  epsi2(nz ,jjj-3)
                E(nz ,jjj-1) =      E(nz ,jjj-3)
                Y(    jjj-1) =      Y(    jjj-3)
                k(nz ,jjj-1) =      k(nz ,jjj-3)
            pfrac(    jjj-0) =  pfrac(    jjj-2)  !node value
            pfrac(    jjj-1) =  pfrac(    jjj-3)  !Cell value
           enddo
          enddo
c
          do nz = n,n+2  !void cell center
                 U(nz,j+1) = zero
                 U(nz,j+1) = zero
              phi(nz ,j+1) = zero
             beta(nz ,j+1) = zero
           sigmar(nz ,j+1) = zero
           sigmao(nz ,j+1) = zero
                V(nz ,j+1) = zero
               s1(nz ,j+1) = zero
               s2(nz ,j+1) = zero
               s3(nz ,j+1) = zero
                E(nz ,j+1) = zero
                Y(    j+1) = zero
                q(nz ,j+1) = zero
            pfrac(    j+1) = zero
            pfrac(    j+1) = zero
          enddo
                ibc(j)   = 2 !outer
                ibc(j+1) = 9 !void
                ibc(j+2) =-2 !inner
              P(n+1,j+1) = pvoid
              P(n+2,j+1) = pvoid
              P(n+1,j-1) = pvoid
              P(n+2,j-1) = pvoid
              P(n+1,j+3) = pvoid
              P(n+2,j+3) = pvoid
              u(n+1,j+2) = u(n+1,j+4)
              r(n+1,j+1) = r(n+1,j)
              r(n+2,j+1) = r(n+2,j)
              pfrac(j+2) =  pfrac(j+4)  !node value
              pfrac(j+3) =  pfrac(j+4)  !node value
        print *,'After Separating:',n,j
        print *,'r(j-2)=',r(n+2,j-2),'U= ',u(n+1,j-2),'P=',P(n+2,j-2)
        print *,'r(j-1)=',r(n+2,j-1),'U= ',u(n+1,j-1),'P=',P(n+2,j-1)
        print *,'r(j)  =',r(n+2,j  ),'U= ',u(n+1,j-0),'P=',P(n+2,j+0)
        print *,'r(j+1)=',r(n+2,j+1),'U= ',u(n+1,j+1),'P=',P(n+2,j+1)
        print *,'r(j+2)=',r(n+2,j+2),'U= ',u(n+1,j+2),'P=',P(n+2,j+2)
        print *,'r(j+3)=',r(n+2,j+3),'U= ',u(n+1,j+3),'P=',P(n+2,j+3)
        print *,'r(j+4)=',r(n+2,j+4),'U= ',u(n+1,j+4),'P=',P(n+2,j+4)
               read(*,*) foo
c
      do jzz=1,jj-2,1
c       write(*,'(I4,5I2,10f7.3)')
       write(*,'(I4,5I2,3f7.3,2e11.3,5f7.3)')
     &  jzz,ieos(1,jzz),ieos(2,jzz),
     &  ibc(jzz-1),ibc(jzz),ibc(jzz+1),
     &  r(n+2,jzz-1),r(n+2,jzz),r(n+2,jzz+1),
     &  U(n+1,jzz-1),U(n+1,jzz),U(n+1,jzz+1),
     &  p(n+2,jzz)
c     &  ,v(n+2,jzz+1),m(jzz+1),rho(jzz+1),y(jzz+1),pfrac(jzz+2)
      enddo
      read(*,*) foo
        endif
      enddo ! this is the fracture check
      endif
c
c******************  Spall End  ***************************
c       TIME-STEP
c     Wilkins section 9
c     Needs local sound speed
c      print *,'line 1258'
       dt_min =  1.d9
       dr_min   =  r(n+2,0+2)-r(n+2,0  )  !this is just to get us started
c
       jzz    = -1
      do j = 0, jj-2, 2
      if (ibc(j+1) .eq. 0) then
c       deltat =     t(n+2)    -t(n)
       Vdot   =   ( V(n+2,j+1)-V(n  ,j+1) )/deltat
       deltar = abs(r(n+2,j+2)-r(n+2,j  ) )
c
       b      = 8.d0*(Co**2+CL)*deltar*(Vdot/V(n+1,j+1))
       if (Vdot/V(n+1,j+1) .ge. 0.) b = 0.d0
c
       rho_local = rho(j+1)/V(n+2,j+1)
c       print *,'rho_local',j,rho_local,rho(j+1),1.d0/v(n+2,j+1)
       if (rho_local .gt. 1.e10) then
       print *,'Error: density too high'
       print *,j,ncount,t(n)
       print *,rho(j+1),v(n+2,j+1)
       print *,j,r(n+2,j+2),r(n+2,j)
       print *,V(n+2,j+1),V(n,j+1)
       print *,P(n+1,j),p(n+1,j+1),P(n+2,j),p(n+2,j+1)
       read(*,*) foo
       endif
c       print *,'rho local',rho_local
       a  = sqrt(abs(P(n+2,j+1))/ rho_local ) !I checked units
c       print *,'ts',j,a,b,deltar,Vdot,V(n+1,j+1)
       if ( (a**2+b**2) .ne. 0.d0 ) then
          delt_temp = (2.d0/3.d0)*(deltar/sqrt(a**2+b**2))
          dtminj(j) = delt_temp
        if (delt_temp .lt. dt_min) then
         jzz    = j
         dt_min = delt_temp
         a_min = a
         b_min = b
         rho_min=rho_local
c        print *,'rr',rho_min,deltat,jzz,a_min,b_min,p(n+2,j+1)
c        read(*,*) foo
        endif
       endif
       endif !ibc(j+1)
      enddo !j loop

      dt_min=dt_min/6.d0    !Here is where you can arbitrarly lower time step
                            !     the demomenator is suppose to be 1.d0
                            ! STS finds that a factor of 6 keeps free ends from spurious oscillations in planar impacts
      delt = t(n+1)-t(n-1)
      if (dt_min .gt. 1.1d0*delt ) then
        dt_min = 1.1d0 * delt
      endif
      ! comment this out for constant time step
      deltat = dt_min ! STS DEBUG this was commented out in online v11 -- so the time step is constant in this code!
c
c      t(n+1) = t(n) + deltat/2           ! this is now advanced along with all the
c      t(n+2) = t(n) + deltat !+ deltat  ! other state variables
c
      if (idebug .eq. 1) then
      if (icontact .eq. 1) then
      Print *,'End of first contact'
             if(idebug .eq.1 )read(*,*) foo
      endif
      endif
c        if (deltat .lt. 0.1d-5) then  ! this is just larry swabies suggestion for stability
c        print *,rho_min,deltat,jzz,a_min,b_min
c        endif
c
c
c ******************* Write solution to File  *************************
c       write(*,*) 'Enter file name'
c       read(*,*) foo name
c       name='al.txt'
c       open (unit=33, file=name, form='formatted', status='unknown')
c     do n=1,nn,2
c     output is the newly calculated time step n+2
c      if (t(n) .eq. 0. .or. t(n) .ge. tskip ) then
      if (t(n) .ge. tskip .and. stepn .gt. 1) then
         qtotal = 0.
         mvtotal = 0.
         ketotal = 0.
         ietotal = 0.
         etotal = 0.
         do j = 0,jj-2,2
            if (ibc(j+1) .eq. 0) then
              qtotal  = qtotal  + q(n+2,j)
              mvtotal = mvtotal + m(j+1)*(  u(n+1,j) + u(n+1,j+2) )/2.d0
              ketotal = ketotal + 0.5*m(j+1)*
     &                  (0.5*((u(n+1,j) + u(n+1,j+2))))**2
              ietotal = ietotal + E(n+2,j+1)/rho(j+1)*m(j+1)
              etotal  = ietotal + ketotal
           endif
         enddo
c        STS OUTPUT TO MATCH pyKO
         do j=0,jj-1,2
             rhoout= 0.
             if (ibc(j) .eq. 0) then
                rhoout = rho(j+1)/v(n+2,j+1)
             endif  
           write(33,'(4I6,26e17.8e3)')
     &           stepn, j, ibc(j), mat(j),
     &           t(n+2), r(n-1,j), r(n,j),r(n+1,j),r(n+2,j),
     &           U(n+1,j), V(n+2,j+1), rhoout, rho(j+1),    
     &           E(n+2,j+1), P(n+2,j+1), sigmar(n,j+1), s1(n+2,j+1),
     &           s2(n+2,j+1),q(n+1,j+1),dtminj(j),phi(n,j),beta(n,j),
     &           epsi1(n+1,j+1),epsi2(n+1,j+1),deltaZarr(j),
     &           alocal(j+1),rholocal(j+1),sigmao(n,j+1),Temp(n+2,j+1),
     &           etotal
         enddo ! this is j loop

c         do j=0,jj-1,1
c             rhoout= 0.
c             if (ibc(j) .eq. 0) then
c                rhoout = rho(j+1)/v(n+2,j+1)
c             endif  
c           write(34,'(4I6,26e17.8e3)')
c     &           stepn, j, ibc(j), mat(j),
c     &           t(n+2), r(n-1,j), r(n,j),r(n+1,j),r(n+2,j),
c     &           U(n+1,j), V(n+2,j+1), rhoout, rho(j+1),    
c     &           E(n+2,j+1), P(n+2,j+1), sigmar(n,j+1), s1(n+2,j+1),
c     &           s2(n+2,j+1),q(n+1,j+1),dtminj(j),phi(n,j),beta(n,j),
c     &           epsi1(n+1,j+1),epsi2(n+1,j+1),deltaZarr(j),
c     &           alocal(j+1),rholocal(j+1),sigmao(n,j+1),Temp(n+2,j+1),
c     &           etotal
c         enddo ! this is j loop
c
         tskip = tskip + dtskip
      endif
c
cccccccccccccccccccccc Update solution ccccccccccccccccccccccccccccc
c
c     shift full timestep into the future
c     current becomes back time step
      do j=0,jj
      if (ibc(j) .ne. 9) then
c
        U(n-1,j)       = U(n+1,j)
        U(n  ,j)       = U(n+2,j)
        phi(n-1,j)     = phi(n+1,j)
        phi(n  ,j)     = phi(n+2,j)
        sigmar(n-1,j)  = sigmar(n+1,j)
        sigmar(n  ,j)  = sigmar(n+2,j)
        sigmao(n-1,j)  = sigmao(n+1,j)
        sigmao(n  ,j)  = sigmao(n+2,j)
        beta(n-1,j)    = beta(n+1,j)
        beta(n  ,j)    = beta(n+2,j)
        V(n-1,j)       = V(n+1,j)
        V(n  ,j)       = V(n+2,j)
        r(n-1,j)       = r(n+1,j)
        r(n  ,j)       = r(n+2,j)
        epsi1(n-1,j)   = epsi1(n+1,j)
        epsi1(n  ,j)   = epsi1(n+2,j)
        epsi2(n-1,j)   = epsi2(n+1,j)
        epsi2(n  ,j)   = epsi2(n+2,j)
        s1(n-1,j)      = s1(n+1,j)
        s1(n  ,j)      = s1(n+2,j)
        s2(n-1,j)      = s2(n+1,j)
        s2(n  ,j)      = s2(n+2,j)
        s3(n-1,j)      = s3(n+1,j)
        s3(n  ,j)      = s3(n+2,j)
        P(n-1,j)       = P(n+1,j)
        P(n  ,j)       = P(n+2,j)
        q(n-1,j)       = q(n+1,j)
        q(n  ,j)       = q(n+2,j)
        t(n-1)         = t(n+1)
        t(n  )         = t(n+2)
        E(n-1,j)       = E(n+1,j)
        E(n  ,j)       = E(n+2,j)
        K(n-1,j)       = K(n+1,j)
        K(n  ,j)       = K(n+2,j)
        Temp(n-1,j)    = Temp(n+1,j)
        Temp(n  ,j)    = Temp(n+2,j)
        entropy(n-1,j) = entropy(n+1,j)
        entropy(n  ,j) = entropy(n+2,j)
c        icompact(n-1,j)= icompact(n+1,j)
c        icompact(n  ,j)= icompact(n+2,j)

      endif
      enddo
cccccccccccccccccccccc  Main loop closed ccccccccccccccccccccccccccccccccccccc
c     write(*,'(e20.5)') 'time=', t(n)
        stepn = stepn+1
c        if (stepn .lt. 10) goto 1
        if (t(n) .le. tstop) goto 1
c      enddo  ! n or time loop

      call cpu_time(endTime)
c     write LAST time step to file
      do j=0,jj-1,2
          rhoout= 0.
          if (ibc(j) .eq. 0) then
             rhoout = rho(j+1)/v(n+2,j+1)
          endif
  
          write(33,'(4I6,26e17.8e3)')
     &           stepn, j, ibc(j), mat(j),
     &           t(n+2), r(n-1,j), r(n,j),r(n+1,j),r(n+2,j),
     &           U(n+1,j), V(n+2,j+1), rhoout, rho(j+1),    
     &           E(n+2,j+1), P(n+2,j+1), sigmar(n,j+1), s1(n+2,j+1),
     &           s2(n+2,j+1),q(n+1,j+1),dtminj(j),phi(n,j),beta(n,j),
     &           epsi1(n+1,j+1),epsi2(n+1,j+1),deltaZarr(j),
     &           alocal(j+1),rholocal(j+1),sigmao(n,j+1),Temp(n+1,j+1),
     &           etotal
      enddo ! this is j loop

c      do j=0,jj-1,1
c          rhoout= 0.
c          if (ibc(j) .eq. 0) then
c             rhoout = rho(j+1)/v(n+2,j+1)
c          endif
c  
c          write(34,'(4I6,26e17.8e3)')
c     &           stepn, j, ibc(j), mat(j),
c     &           t(n+2), r(n-1,j), r(n,j),r(n+1,j),r(n+2,j),
c     &           U(n+1,j), V(n+2,j+1), rhoout, rho(j+1),    
c     &           E(n+2,j+1), P(n+2,j+1), sigmar(n,j+1), s1(n+2,j+1),
c     &           s2(n+2,j+1),q(n+1,j+1),dtminj(j),phi(n,j),beta(n,j),
c     &           epsi1(n+1,j+1),epsi2(n+1,j+1),deltaZarr(j),
c     &           alocal(j+1),rholocal(j+1),sigmao(n,j+1),Temp(n+1,j+1),
c     &           etotal
c      enddo ! this is j loop
c     
      close(33)
c      close(34)
c      close(35)
      print *,'Fortran KO time in main loop = ',endTime-startTime
      print *,'*******************  Finished!!  ******************'
c      read(*,*) foo
      end
