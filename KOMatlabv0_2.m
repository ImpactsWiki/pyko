%% Slower KO program for ease of coding
% Adapted from Mark L. Wilkins KO Program
% By Nathaniel Helminiak
% Fortran Version Maintained by John Borg
%https://www.eng.mu.edu/shockphysics/KO/
% C Version Maintained by David Helminiak
%https://github.com/Yatagarasu50469/KO-Hydrocode
%-------------------Version History-------------------
% V0_1 Coded Version
% V0_2 Updated Timesetting and added output
%-------------------SUGGESTED TO DO-------------------
% Add better comments and clean
% Add input deck
% Add more strength models such as JC
% Add material library
% Add local timestep based on local info speed
% Add parallelization

%% Clean Workspace
clc
clear
close all

% c      KO see Springer book by Wilkins, M.
% c              "Computer simulations of Dynamic Phenomena"
%        parameter(jj = 1510)       !total number of nodes allowed
%                                  !this number should be greater than the number
%                                  !of nodes that you want so that nodes can be added 
%                                  !if fracture occurs
%        parameter(jmat=7)        !number of materials
%        parameter(nmax=4)        !number of placeholders for time: t(0) and t(1) current, t(2) and t(3) new
%        integer n, j,imat,ibc(0:jj), icontact,idebug,jdebug,iskip,ipoint
%        integer j1(jj),iEOS(2,0:jj),jjj,jzz,ndebug,i,amr
%        integer icompact(0:nmax,jj)
% 
%        real*8  m(0:jj), r(0:nmax,0:jj), U(0:nmax,0:jj)
%        real*8  t(0:nmax),pfrac(0:jj),pvoid
%        real*8  phi(0:nmax,0:jj), sigmar(0:nmax,0:jj)
%        real*8  sigmao(0:nmax,0:jj),Temp(0:nmax,0:jj)
%        real*8  beta(0:nmax,0:jj), P(0:nmax,0:jj), q(0:nmax,0:jj)
%        real*8  s1(0:nmax,0:jj),s2(0:nmax,0:jj), s3(0:nmax,0:jj)
%        real*8  rho(0:jj), V(0:nmax,0:jj), entropy(0:nmax,0:jj)
%        real*8  epsi1(0:nmax,0:jj),epsi2(0:nmax,0:jj),epdt1(0:nmax,0:jj)
%        real*8  E(0:nmax,0:jj), K(0:nmax,0:jj), Y(0:jj)
%        real*8  deltaZ, Vdot ,dt_min, dr_min,delt_temp
%        real*8  deltar,  qbar
%        real*8  a, b, Co, CL, d ,a_min,b_min,rho_min
%        real*8  EOS(0:jmat,13),init(jmat,16),bc(-9:9,5),LL(jj),xstart(jj)
%        real*8  xx,xa,xb,deltat_0,deltat,delt,tstop
% c       real*8  U_0,P_0,V_0,E_0,rho_0,a_0,up
% c       real*8  U_1,P_1,V_1,E_1,a,
%        real *8 zero
%        real*8  aa,ww,rr,bb,dti0,dti1,Usave
%        real*8  phiv,phij,betav,betaj,dthalf
%        real*8  qtotal,mvtotal,ketotal,ietotal,etotal,ke3total
%        real*8  bs1,bs2,bs3,bs4,bs5,bs6,bs7,bs8,bs9,bs10
%        real*8  tskip,dtskip
%        real*8  rho_local,stemp,ctemp,gtemp,v0,v00,vv
% c       real*8  alpha
% c      real*8  pe,ps,ap,ae
%        real*8  k1,k2,k3,gamma0
%        real*8  Us,up,PH,EH,TH,strain,P0,E0,T0
%        real*8  En2j1,diffE
%        character*20 name
%        character*172 text

%% Initial Variables
tstop=2        %Stop Time

jj = 20*6;     %total number of nodes allowed this number should be greater than the number
               %of nodes that you want so that nodes can be added if fracture occurs
jmat=7;        %number of materials
nmax=4;        %NSH number of placeholders for time: t(1) and t(2) current, t(3) and t(4) new

Co        = 2;  %artifical viscoity constants, default=2
CL        = 1;  %artifical viscoity constants, default=1
pvoid     = 0;  %pressure in the void
icontact  = 0;  %contact flag this should not be changed but i think it could
                    % be elimated if the code were cleaned up.
idebug    = 0;      % turn on debug
ndebug    = 1;      % what time steup to you want to print to screen
jdebug    = 2;      % what node do you want to look at when debugging
d         = 1;      %defines geometry (1-1D planar, 2-1D cylinderical-never tested), see willkins
deltat_0  = 3.0d-3;
delt      = 3.0d-3;
iskip     = 1;      %number of iterations to skip to echo to screen
tskip     = 0;      % this is the value at which data starts getting written to file
dtskip    = 0.01;   % This is the amount of time skipped between writes to file
amr       = 0;      % this is an adaptive mesh refinement (1 is on - 0 is off)
flag      = 0;
init=zeros(jmat,16);
P=zeros(nmax+1,jj+1);
m=zeros(1,jj+1);
E=zeros(nmax+1,jj+1);
V=zeros(nmax+1,jj+1);
rho=zeros(1,jj+1);

%% Zero All Variables
fprintf('Initializing variables ...\n')
zero = 0;
for n = 0+1:nmax+1
    for j = 0+1:jj+1
        ibc(j)      = 9;  % these are non-used nodes
        U(n,j)      = zero;
        U2(n,j)      = zero;
        r(n,j)      = zero;
        t(n)        = deltat_0*single(n-1);
        deltat      = deltat_0;
        phi(n,j)    = zero;
        sigmar(n,j) = zero;
        sigmao(n,j) = zero;
        beta(n,j)   = zero;
        q(n,j)      = zero;
        s1(n,j)     = zero;
        s2(n,j)     = zero;
        s3(n,j)     = zero;
        epsi1(n,j)  = zero;
        epsi2(n,j)  = zero;
        K(n,j)      = zero;
        Y(j)        = zero;
        deltaZ      = zero;
        pfrac(j)    = zero;
        Temp(n,j)   = zero;
        entropy(n,j)= zero;
    end
end

for j=1:jmat
    for i=1+1:7+1
        eos(j,i)= 0.d0;
    end
    for i=1:4
        init(j,i)= 0.d0;
    end
end
for j=1+1:jj+1
    for jzz=1:2
        ieos(jzz,j)= 0;
    end
end

%% Read In Materials
fprintf('running\n')
      imat = 0; %number of materials
      jsum = 0;
%          jGeometry:            cm     cm     Mbar  cm/us   g/cc          g/cc    cm/us unitless                Megabar Megbar Megbar Mbar-cc/K-g
%                   1      2      3      4       5      6      7      8      9       10      11       12     13      14     15     16          17     18     19     20     21     22
%                iEOS  Nodes Length xstart     P_0    U_0  Rho_0    E_0 Rho_00       Co       s       s2  gamma  Yeld,Y   mu,G  pfrac          Cv    a_p    a_e    p_e    p_s  rho_0               
%MAT1
MATPROP(1,:)=[      1  16    0.2       -0.21  0.0E-6 0.0200   8.0  0.000    8.0    0.460    1.49     0.00   2.17   2.5e-3     .79  1.0e9      1.0e-9   0.00   0.00   0.00   0.00   0.00];
%MAT2
MATPROP(2,:)=[      1  32   0.4        0.0  0.0E-6  0.0000   6.0  0.000    6.0    0.460    1.49     0.00   2.17   2.5e-3     .79  0.25e-1      1.0e-9   0.00   0.00   0.00   0.00   0.00];


temp=size(MATPROP);
imat=temp(1);
      for i=1:imat
        ieos(1,jj-i+2) = i
      end

      for imat=1:size(MATPROP)
          ieos(2,jj-(imat-1)+1) = MATPROP(imat,1)
          j1(imat)            = MATPROP(imat,2);
          LL(imat)            = MATPROP(imat,3);
          xstart(imat)        = MATPROP(imat,4);
          init(imat,1)        = MATPROP(imat,5);
          init(imat,2)        = MATPROP(imat,6);
          init(imat,3)        = MATPROP(imat,7);
          init(imat,4)        = MATPROP(imat,8);
          init(imat,5)        = MATPROP(imat,9);
          eos(imat,1)         = MATPROP(imat,10);
          eos(imat,2)         = MATPROP(imat,11);
          eos(imat,3)         = MATPROP(imat,12);
          eos(imat,4)         = MATPROP(imat,13);
          eos(imat,5)         = MATPROP(imat,14);
          eos(imat,6)         = MATPROP(imat,15);
          eos(imat,7)         = MATPROP(imat,16);
          eos(imat,8)         = MATPROP(imat,17);
          eos(imat,9)         = MATPROP(imat,18);
          eos(imat,10)        = MATPROP(imat,19);
          eos(imat,11)        = MATPROP(imat,20);
          eos(imat,12)        = MATPROP(imat,21);
          eos(imat,13)        = MATPROP(imat,22);
      end
            
       if ( int32(j1(imat)/2) ~= single(j1(imat))/2.)
        fprintf('Error: nodes specificed for a material must be even\n')
       else
           fprintf('OK\n');
       end
       jsum = jsum + j1(imat);
       if (jsum >= jj)
        fprintf('Error: Increase jj or reduce the # nodes in kov3.in');
       end
       fprintf('data at j %2.0f  CV %2.0f\n',imat,eos(imat,7));
       fprintf('gg Yield %2.0f\n',eos(imat,4));
       fprintf('End of material data input file...\n');


%% Boundary Conditions  % This ONE is A bit strange
                                                                                                                                                                        
% Boundary Conditions:
%       1      2      3      4      5      6
%    iend    P_0    U_0  Rho_0    E_0    V_0                                                                                                                        
%      -1 0.0e-1     .0   2.70  0.000    1.0                                                                                                                                  
%       1 0.0e-1     .0   2.70  0.000    1.0  
BOUNDARY(1,:)=[      -1 0.0e-1     .0   8  0.000    1.0];
BOUNDARY(2,:)=[       1 0.0e-1     .0   6  0.000    1.0  ];

 ibc(1) =BOUNDARY(1,1);
 bc(1,1)=BOUNDARY(1,2);
 bc(1,2)=BOUNDARY(1,3);
 bc(1,3)=BOUNDARY(1,4);
 bc(1,4)=BOUNDARY(1,5);
 bc(1,5)=BOUNDARY(1,6);
 ibc(jj+1)=BOUNDARY(2,1);
 bc(2,1)=BOUNDARY(2,2);
 bc(2,2)=BOUNDARY(2,3);
 bc(2,3)=BOUNDARY(2,4);
 bc(2,4)=BOUNDARY(2,5);
 bc(2,5)=BOUNDARY(2,6);
 
 
 %% Discritization LOOP
 fprintf('distrize ..')
 r(0+1,1) = xstart(1);
 r(1+1,1) = xstart(1);

 for jjj=1:imat
     fprintf('distrize ..%2.0f\n',jjj)
     deltar = LL(jjj)/(j1(jjj)/2);
     %write(*,'(a6,1x,3e24.16)') 'deltar',deltar
     if (jjj == 1)                     % this if assigns the ipoint to the first node of the material
         %  and checks to see it materials are initally in contact.
         ipoint = 1;
     elseif (abs(r(0+1,ipoint+j1(jjj-1))-xstart(jjj)) < 1.e-5)
         %       no gap between materials
         r(0+1,ipoint+j1(jjj-1)) = xstart(jjj);
         r(1+1,ipoint+j1(jjj-1)) = xstart(jjj);
         ibc(ipoint+j1(jjj-1)) = 0;
         ipoint                = ipoint+j1(jjj-1);
     elseif ( r(0+1,ipoint+j1(jjj-1)) < xstart(jjj))
         %       inital gap between materials
         ipoint      = ipoint+j1(jjj-1)+2;
         ibc(ipoint  ) = -2;
         ieos(2,ipoint-1)= 0;
         ibc(ipoint-2) =  2;
     else
         fprintf('Input Geometry Error!')
     end
     
     r(1,ipoint)  = xstart(jjj);
     r(1+1,ipoint)  = xstart(jjj);
     ibc(1+1)       = 0;
     U(0+1,ipoint)  = init(jjj,2);
     U(1+1,ipoint)  = init(jjj,2);
     for j=ipoint+2:2:ipoint+j1(jjj)      %node definitions
         ibc(j)       = 0;      %this defines the node as a centeral difference
         r(1,j)       = r(0+1,j-2) + deltar;
         r(2,j)       = r(1+1,j-2) + deltar;
         U(1,j)       = init(jjj,2);
         U(2,j)       = init(jjj,2);
     end
     for j=ipoint+1:2:ipoint+j1(jjj)-1   %cell definitions initial conditions
         ibc(j)       = 0;      %this defines the cell as a centeral difference
         r(1,j)       = 0.5*(r(0+1,j+1)+r(0+1,j-1));
         r(2,j)       = 0.5d0*(r(1+1,j+1)+r(1+1,j-1));
         P(1,j)       = init(jjj,1);
         P(2,j)       = init(jjj,1);
         rho(j)       = init(jjj,3);  %rho is not updated in time therefore it is always rho0
         E(1,j)       = init(jjj,4);  % this gets overwritten below
         E(2,j)       = init(jjj,4);
         ieos(1,j)    = ieos(1,jj-(jjj-1)+1);  %info was stored there (RHS) just temp.
         ieos(2,j)    = ieos(2,jj-(jjj-1)+1);  %info was stored there (RHS) just temp.
         Y(j)         = eos(jjj,5);
         pfrac(j-1)   = eos(jjj,7);
         pfrac(j)     = eos(jjj,7);
         pfrac(j+1)   = eos(jjj,7);
     end
     ieos(1,jj-(jjj-1)+1)     = 0;   %this assigns values to the last cell center
     ieos(2,jj-(jjj-1)+1)     = 0;
     ibc(j1(jjj)+ipoint)    = 2;
     ibc(j1(jjj)+1+ipoint)  = 9;
 end  % jjj loop through imat
 %NATHAIEL ibc(jj) should be 1
 ibc(j1(imat)+ipoint)    = ibc(jj+1);  %this assigns values to the last node
  ibc(jj+1)=9;
fprintf('ipoint %2.0f %2.0f\n',j1(imat)+ipoint,ibc(jj+1));
 fprintf('Assign mass to nodes\n');
 for j=0+1:2:jj-2+1
     if (ibc(j+1) == 0)
         m(j+1) = rho(j+1)*((r(0+1,j+2)^d - r(0+1,j)^d)/d);
     end
 end
 
%% Set the Initial Volume

%ieos screwed up j=302
for j=0+1:2:jj-2+1
    if (ibc(j+1) == 0)
        V(0+1,j+1)=rho(j+1)*((r(0+1,j+2)^d-r(0+1,j)^d)/d)/m(j+1);
        V(1+1,j+1)=rho(j+1)*((r(0+1,j+2)^d-r(0+1,j)^d)/d)/m(j+1);
        E(0+1,j+1)=P(0+1,j+1)/(V(0+1,j+1)*eos(ieos(1,j+1),4)-1);
        E(1+1,j+1)=P(1+1,j+1)/(V(1+1,j+1)*eos(ieos(1,j+1),4)-1);
        Temp(0+1,j+1) =E(0+1,j+1)/(rho(j+1)*eos(ieos(1,j+1),8));
        Temp(1+1,j+1) =E(1+1,j+1)/(rho(j+1)*eos(ieos(1,j+1),8));
        entropy(0+1,j+1) = 6.8d-5; %! i got this from a VT website EES thing
        entropy(1+1,j+1) = 6.8d-5; %! units mbar-cc/K/g
        fprintf('Pres [Mbar] %2.0f \n',P(0+1,j+1));
        fprintf('Engr [MBar-cc/cc] %2.0f \n',E(0+1,j+1));
        fprintf('Volm [cc] %2.0f \n',((r(0+1,j+2)^d-r(0+1,j)^d)/d));
        fprintf('temp*[K] %2.0f %2.0f %2.0f %2.0f \n',Temp(0+1,j+1),rho(j+1),m(j+1),V(0+1,j+1));
        fprintf('entropy %2.0f %2.0f \n',entropy(0+1,j+1),eos(ieos(1+1,j+1),8));
    end
end

%% Output Inital Conditions to screen
n=1+1;
         

 
%% MAIN LOOP
%%%%%%%%%%%%
ncount = 0
timesteps=1

%%
while (t(n) <= tstop)
%for timesteps=1:2000%147%147
    
%      clc
%      clear
%      load('SaveTimeMassCrash')
%      timesteps=234
    %fprintf('%2.0f \n',timesteps);
%     ncount = ncount + 2 % when this use to count to an integer, ncount was the max integer
%     n     =  1    % this use to be the time step counter when all n data was stored
%     % n must remain n=1 because the time data is no longer stored
    if( m(1) ~= m(3) )
        fprintf('Masses do not match %2.0f %2.0f %2.0f',n,m(1),m(3))
    else
        %disp('ok')
    end
     
    %% Contact Check
    
    if (1 == 1)
        for j=0+1:2:jj+1
            if (ibc(j) ~= 9)
                %delt     =  (t(n+2)-t(n+1))/2;
                r(n+2,j) =  r(n,j)+U(n-1,j)*deltat;
                r(n+1,j) = (r(n,j)+r(n+2,j))/2;
            end
        end
        % % c
        
        % Loop below needs checking!
        % begin 144 -1     0     0     0     0     0     0     0     0     0     0     0     1     9     9     9     9     9     1     9     1
        for j=2+1:2:jj-2+1
            if ( ibc(j) ~= 9 && r(n+2,j) <= r(n+2,j-2) &&  ibc(j-1) == 9)
                j
                icontact = 2;
                ibc(j) = -3  %this is a flag to let the contact stuff after the
                % momentum step know that
                % the delta t was recaclulated so that the two nodes
                % that are about to collide will touch perfectely after
                % this time step.
                fprintf('Contact Eminent! %2.4f %2.4f\n',n,j);
                fprintf('Previous half Time step %2.4f\n',t(n)-t(n-1));
                fprintf('For time step 2 x Deltat %2.4f\n',2.d0*delt);
                fprintf(' Estimated collision:\n');
                fprintf('Void Node at %2.4f will move to %2.4f\n',r(n,j-2), r(n,j-2)+U(n-1,j-2)*2.d0*deltat);
                fprintf('Void node velocity %2.4f\n',U(n-1,j-2));
                fprintf('J Node at %2.4f will move to %2.4f\n',r(n,j),r(n,j)+U(n-1,j)*2.d0*delt);
                fprintf('J node velocity %2.4f\n',U(n-1,j));
                
                jjv       = j-2;
                jjj       = j;
                ww        = U(n-1,jjj) - U(n-1,jjv);
                rr        = r(n  ,jjj) - r(n  ,jjv);
                %
                % this is the ibc(j)= 1 for the inside (v) side of the boundary
                sigmar(n,jjv-1) = (-(P(n,jjv-1)+q(n-1,jjv-1))+s1(n,jjv-1));
                sigmao(n,jjv-1) = (-(P(n,jjv-1)+q(n-1,jjv-1))+s2(n,jjv-1));
                phiv      = (rho(jjv-1)*((r(n,jjv)-r(n,jjv-2))/V(n,jjv-1)))/2.d0;
                betav     = (sigmar(n,jjv-1)-sigmao(n,jjv-1))*V(n,jjj-1)/(r(n,jjv-1)*rho(jjv-1));
                % this is the ibc(j)=-1 for the outside (j) side of the boundary
                sigmar(n,jjj+1) = (-(P(n,jjj+1)+q(n-1,jjj+1))+s1(n,jjj+1));
                sigmao(n,jjj+1) = (-(P(n,jjj+1)+q(n-1,jjj+1))+s2(n,jjj+1));
                phij      = (rho(jjj+1)*((r(n,jjj+2)-r(n,jjj))/V(n,jjj+1)))/2.d0;
                betaj     = (sigmar(n,jjj+1)-sigmao(n,jjj+1))*V(n,jjj+1)/(r(n,jjj+1)*rho(jjj+1));
                aa = sigmar(n,jjj+1)/phij + sigmar(n,jjv-1)/phiv + (betav + betaj)*(d-1.d0);
                dthalf = t(n) - t(n-1);
                bb   = (2.d0*ww + aa * dthalf);
                dti0 = 0.d0;
                dti1 = dti0-(aa*dti0*dti0+bb*dti0+2.d0*rr)/(2.d0*aa*dti0+bb);
                fprintf('sig %2.4f %2.4f %2.4f %2.4f \n',sigmar(n,jjj+1),phij,sigmar(n,jjv-1),phiv);
                fprintf('bet %2.4f %2.4f %2.4f \n',betav,betaj,d);
                fprintf('%2.4f %2.4f %2.4f %2.4f \n',aa,bb,rr,ww);
                fprintf('%2.4f %2.4f %2.4f %2.4f \n',rho(jjv-1),r(n,jjv),r(n,jjv-2),V(n,jjv-1));
                fprintf('%2.4f %2.4f %2.4f %2.4f \n',rho(jjj+1),r(n,jjj+2),r(n,jjj),V(n,jjj+1));
                fprintf('betav %2.4f %2.4f %2.4f %2.4f %2.4f %2.4f',sigmar(n,jjv-1),sigmao(n,jjv-1),r(n,jjj),r(n,jjj-2),V(n,jjj-1),rho(jjj-1));
                % % c
                jj1=0;
                while (abs(dti1-dti0) > 1.d-16)
                    jj1 = jj1 + 1;
                    dti0 = dti1;
                    dti1 = dti0-(aa*dti0*dti0+bb*dti0+2.d0*rr)/(2.d0*aa*dti0+bb);
                    fprintf('%2.4f %2.4f %2.4f %2.4f\n',aa,bb,rr,ww);
                    fprintf('iteration\n',dti0,dti1);
                    if (jj1 > 1)
                        %read(*,*) foo
                    end
                end
                fprintf('Whole Time step adjusted from %2.4f \n',2.d0*delt);
                deltat = dti1/2.d0;   %this might need to be adjusted because i changed the way time steping works
                fprintf('To %2.4f\n',(deltat+dti1));
                fprintf('where %2.4f is the time to impact \n',dti1);
                fprintf('and %2.4f is an arbatrary small number. \n',deltat);
                %read(*,*) foo
            end %contact if
        end
        %read(*,*) foo
    end   %!contact check if
    %%
    
% % c
% % ccccccccccccccccccccccc  Advance Time step  cccccccccccccccccccccccccccccc
% % c
% % c     advance time step, for first calculation timestep is set by deltat_0
% % c     after that the time step is calculated after the end of the current time step, deltat
% % c
if (icontact == 2)     % if there was contact time step adjusted to perfectly
    t(n+1) = t(n) + dti1/2.d0;   % have nodes touch at the end of the time step
    t(n+2) = t(n+1) + deltat;
    icontact = 0;
    flag=1;
    if  (dti1 - deltat < 0. || t(n+2)-t(n+1) < 0.)
        fprintf('Time step Error from contact'\n);
        fprintf('%2.0f %2.0f %2.0f',dti1,deltat,dti1 - deltat , t(n+2)-t(n+1));
    end
else
    %deltat=deltat_0;
    %t(n+1) = t(n) + (deltat/2.d0)/lower
    %t(n+2) = t(n) + (deltat)/lower
    t(n+1) = t(n) + deltat/2.d0; %Higher
    t(n+2) = t(n) + deltat;  %Higher
    
end

%%
%ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
%c     2.    Conservation of Momumtum  - start j
%c
%c  Boundary Conditions
%c   -/+ 1 Free End, stress free
%c   -/+ 2 and 3 Part of contact check
%c   -/+ 4 Fixed end
%c
%c System is failing in this do loop line 450	
      for j=0+1:2:jj-1+1
      if (ibc(j+1) == 0)  %this is true for central difference cells
       sigmar(n,j+1) = (-(P(n,j+1)+q(n-1,j+1))+s1(n,j+1));
       sigmao(n,j+1) = (-(P(n,j+1)+q(n-1,j+1))+s2(n,j+1));
      end
      end
      
      %%
             delt = t(n+1)-t(n-1);   %this is different than deltat, deltat is the true timestep
%c	
      for j=0+1:2:jj+1
%c

      if (ibc(j) ~= 9) %13
       if (ibc(j) == -1)
                              %%%%%%%%<Problem HERE NSH
        phi(n,j) = (rho(j+1)*((r(n,j+2)-r(n,j))/V(n,j+1)))/2.d0;
        beta(n,j) = (sigmar(n,j+1)-sigmao(n,j+1))*V(n,j+1)/(r(n,j+1)*rho(j+1));
        if (delt/phi(n,j))*(sigmar(n,j+1)+0.000000000000)+delt*beta(n,j)*(d-1.d0) == 0
        U(n+1,j) = U(n-1,j);
        else
        U(n+1,j) = U(n-1,j)+(delt/phi(n,j))*(sigmar(n,j+1)+0.000000000000)+delt*beta(n,j)*(d-1.d0);
        end
%c       print *,'BC thing',bc(-1,1)
% What the heck is bc(-1,1)
       elseif (ibc(j) == 1)
        phi(n,j) =  rho(j-1)*((r(n,j)-r(n,j-2))/V(n,j-1))/2.d0;
        beta(n,j) = (sigmar(n,j-1)-sigmao(n,j-1))*V(n,j-1)/(r(n,j-1)*rho(j-1));
        if (delt/phi(n,j))*(-bc( 1,1)-sigmar(n,j-1))+delt*beta(n,j)*(d-1.d0) == 0
        U(n+1,j) = U(n-1,j);
        else
        U(n+1,j) = U(n-1,j)+(delt/phi(n,j))*(-bc( 1,1)-sigmar(n,j-1))+delt*beta(n,j)*(d-1.d0);            
        end
%c       U(n+1,j) = 0.d0
%c
%c this is the new way which treats the voids like inner or outer bc's
       elseif (ibc(j) == -2 || ibc(j) == -3) 
        phi(n,j) = (rho(j+1)*((r(n,j+2)-r(n,j))/V(n,j+1)))/2.d0;
        beta(n,j) = (sigmar(n,j+1)-sigmao(n,j+1))*V(n,j+1)/(r(n,j+1)*rho(j+1));
        U(n+1,j) = U(n-1,j)+(delt/phi(n,j))*(sigmar(n,j+1)+pvoid)+delt*beta(n,j)*(d-1.d0);
       elseif (ibc(j) == 2) 
        phi(n,j) =  rho(j-1)*((r(n,j)-r(n,j-2))/V(n,j-1))/2.d0;
        beta(n,j) = (sigmar(n,j-1)-sigmao(n,j-1))*V(n,j-1)/(r(n,j-1)*rho(j-1));
        U(n+1,j) = U(n-1,j)+(delt/phi(n,j))*(pvoid-sigmar(n,j-1))+delt*beta(n,j)*(d-1.d0);
       elseif (ibc(j) == 0)
        phi(n,j) = (0.5d0)*(rho(j+1)*((r(n,j+2)-r(n,j))/V(n,j+1))+rho(j-1)*((r(n,j)-r(n,j-2))/V(n,j-1)));
        beta(n,j)=( (sigmar(n,j+1)-sigmao(n,j+1))*(V(n,j+1)/rho(j+1))/(0.5d0*(r(n,j+2)+r(n,j  ))) +(sigmar(n,j-1)-sigmao(n,j-1))*(V(n,j-1)/rho(j-1))/(0.5d0*(r(n,j  )+r(n,j-2))) )/2.d0;
        if (delt/phi(n,j))*(sigmar(n,j+1)-sigmar(n,j-1))+delt*beta(n,j)*(d-1.d0) == 0
        U(n+1,j) = U(n-1,j);
        else
        U(n+1,j) = U(n-1,j)+(delt/phi(n,j))*(sigmar(n,j+1)-sigmar(n,j-1))+delt*beta(n,j)*(d-1.d0);            
        end
       elseif (ibc(j) == -4) 
        U(n+1,j) = 0.d0;
%c       if ( U(n+1,j) .lt. 1.d-5) U(n+1,j) = 0.d0
       elseif (ibc(j) == 4)
        U(n+1,j) = 0.d0;
%c       if ( U(n+1,j) .lt. 1.d-5) U(n+1,j) = 0.d0
       else
        fprintf('Momentum ERROR!!! \n');
        fprintf('%2.0f %2.0f',j,ibc(j));
       end
%c
      end
      end    
      
      
%% Increase or decrease timestep based on convergence?
if flag == 0
             delt2 = (t(n+1)-t(n-1))/4;   %NSH for self adjusting timestep
%c	
      for j=0+1:2:jj+1
%c

      if (ibc(j) ~= 9) %13
       if (ibc(j) == -1)
                              %%%%%%%%<Problem HERE NSH
        phi2(n,j) = (rho(j+1)*((r(n,j+2)-r(n,j))/V(n,j+1)))/2.d0;
        beta2(n,j) = (sigmar(n,j+1)-sigmao(n,j+1))*V(n,j+1)/(r(n,j+1)*rho(j+1));
        if (delt2/phi2(n,j))*(sigmar(n,j+1)+0.000000000000)+delt2*beta2(n,j)*(d-1.d0) == 0
        U2(n+2,j) = U(n-1,j);
        else
        U2(n+2,j) = U(n-1,j)+(delt2/phi2(n,j))*(sigmar(n,j+1)+0.000000000000)+delt2*beta2(n,j)*(d-1.d0);
        end
%c       print *,'BC thing',bc(-1,1)
% What the heck is bc(-1,1)
       elseif (ibc(j) == 1)
        phi2(n,j) =  rho(j-1)*((r(n,j)-r(n,j-2))/V(n,j-1))/2.d0;
        beta2(n,j) = (sigmar(n,j-1)-sigmao(n,j-1))*V(n,j-1)/(r(n,j-1)*rho(j-1));
        if (delt2/phi2(n,j))*(-bc( 1,1)-sigmar(n,j-1))+delt2*beta2(n,j)*(d-1.d0) == 0
        U2(n+2,j) = U(n-1,j);
        else
        U2(n+2,j) = U(n-1,j)+(delt2/phi2(n,j))*(-bc( 1,1)-sigmar(n,j-1))+delt2*beta2(n,j)*(d-1.d0);            
        end
%c       U(n+2,j) = 0.d0
%c
%c this is the new way which treats the voids like inner or outer bc's
       elseif (ibc(j) == -2 || ibc(j) == -3) 
        phi2(n,j) = (rho(j+1)*((r(n,j+2)-r(n,j))/V(n,j+1)))/2.d0;
        beta2(n,j) = (sigmar(n,j+1)-sigmao(n,j+1))*V(n,j+1)/(r(n,j+1)*rho(j+1));
        U2(n+2,j) = U(n-1,j)+(delt2/phi2(n,j))*(sigmar(n,j+1)+pvoid)+delt2*beta2(n,j)*(d-1.d0);
       elseif (ibc(j) == 2) 
        phi2(n,j) =  rho(j-1)*((r(n,j)-r(n,j-2))/V(n,j-1))/2.d0;
        beta2(n,j) = (sigmar(n,j-1)-sigmao(n,j-1))*V(n,j-1)/(r(n,j-1)*rho(j-1));
        U2(n+2,j) = U(n-1,j)+(delt2/phi2(n,j))*(pvoid-sigmar(n,j-1))+delt2*beta2(n,j)*(d-1.d0);
       elseif (ibc(j) == 0)
        phi2(n,j) = (0.5d0)*(rho(j+1)*((r(n,j+2)-r(n,j))/V(n,j+1))+rho(j-1)*((r(n,j)-r(n,j-2))/V(n,j-1)));
        beta2(n,j)=( (sigmar(n,j+1)-sigmao(n,j+1))*(V(n,j+1)/rho(j+1))/(0.5d0*(r(n,j+2)+r(n,j  ))) +(sigmar(n,j-1)-sigmao(n,j-1))*(V(n,j-1)/rho(j-1))/(0.5d0*(r(n,j  )+r(n,j-2))) )/2.d0;
        if (delt2/phi2(n,j))*(sigmar(n,j+1)-sigmar(n,j-1))+delt2*beta2(n,j)*(d-1.d0) == 0
        U2(n+2,j) = U(n-1,j);
        else
        U2(n+2,j) = U(n-1,j)+(delt2/phi2(n,j))*(sigmar(n,j+1)-sigmar(n,j-1))+delt2*beta2(n,j)*(d-1.d0);            
        end
       elseif (ibc(j) == -4) 
        U2(n+2,j) = 0.d0;
%c       if ( U(n+2,j) .lt. 1.d-5) U(n+2,j) = 0.d0
       elseif (ibc(j) == 4)
        U2(n+2,j) = 0.d0;
%c       if ( U(n+2,j) .lt. 1.d-5) U(n+2,j) = 0.d0
       else
        fprintf('Momentum ERROR!!! \n');
        fprintf('%2.0f %2.0f',j,ibc(j));
       end
%c
      end
      end    
      
      
    if  sum((U(n+1,:)-U2(n+2,:))/nnz(U(n+1,:))) <  1e-4 % This is current Recomendation 1e-7
        fprintf('increasetime\n')
        deltat=deltat*1.01;
        t(n+1) = t(n) + deltat/2.d0; %Higher
        t(n+2) = t(n) + deltat;  %Higher
                  delt = t(n+1)-t(n-1);   %this is different than deltat, deltat is the true timestep
%c	
      for j=0+1:2:jj+1
%c

      if (ibc(j) ~= 9) %13
       if (ibc(j) == -1)
                              %%%%%%%%<Problem HERE NSH
        phi(n,j) = (rho(j+1)*((r(n,j+2)-r(n,j))/V(n,j+1)))/2.d0;
        beta(n,j) = (sigmar(n,j+1)-sigmao(n,j+1))*V(n,j+1)/(r(n,j+1)*rho(j+1));
        if (delt/phi(n,j))*(sigmar(n,j+1)+0.000000000000)+delt*beta(n,j)*(d-1.d0) == 0
        U(n+1,j) = U(n-1,j);
        else
        U(n+1,j) = U(n-1,j)+(delt/phi(n,j))*(sigmar(n,j+1)+0.000000000000)+delt*beta(n,j)*(d-1.d0);
        end
%c       print *,'BC thing',bc(-1,1)
% What the heck is bc(-1,1)
       elseif (ibc(j) == 1)
        phi(n,j) =  rho(j-1)*((r(n,j)-r(n,j-2))/V(n,j-1))/2.d0;
        beta(n,j) = (sigmar(n,j-1)-sigmao(n,j-1))*V(n,j-1)/(r(n,j-1)*rho(j-1));
        if (delt/phi(n,j))*(-bc( 1,1)-sigmar(n,j-1))+delt*beta(n,j)*(d-1.d0) == 0
        U(n+1,j) = U(n-1,j);
        else
        U(n+1,j) = U(n-1,j)+(delt/phi(n,j))*(-bc( 1,1)-sigmar(n,j-1))+delt*beta(n,j)*(d-1.d0);            
        end
%c       U(n+1,j) = 0.d0
%c
%c this is the new way which treats the voids like inner or outer bc's
       elseif (ibc(j) == -2 || ibc(j) == -3) 
        phi(n,j) = (rho(j+1)*((r(n,j+2)-r(n,j))/V(n,j+1)))/2.d0;
        beta(n,j) = (sigmar(n,j+1)-sigmao(n,j+1))*V(n,j+1)/(r(n,j+1)*rho(j+1));
        U(n+1,j) = U(n-1,j)+(delt/phi(n,j))*(sigmar(n,j+1)+pvoid)+delt*beta(n,j)*(d-1.d0);
       elseif (ibc(j) == 2) 
        phi(n,j) =  rho(j-1)*((r(n,j)-r(n,j-2))/V(n,j-1))/2.d0;
        beta(n,j) = (sigmar(n,j-1)-sigmao(n,j-1))*V(n,j-1)/(r(n,j-1)*rho(j-1));
        U(n+1,j) = U(n-1,j)+(delt/phi(n,j))*(pvoid-sigmar(n,j-1))+delt*beta(n,j)*(d-1.d0);
       elseif (ibc(j) == 0)
        phi(n,j) = (0.5d0)*(rho(j+1)*((r(n,j+2)-r(n,j))/V(n,j+1))+rho(j-1)*((r(n,j)-r(n,j-2))/V(n,j-1)));
        beta(n,j)=( (sigmar(n,j+1)-sigmao(n,j+1))*(V(n,j+1)/rho(j+1))/(0.5d0*(r(n,j+2)+r(n,j  ))) +(sigmar(n,j-1)-sigmao(n,j-1))*(V(n,j-1)/rho(j-1))/(0.5d0*(r(n,j  )+r(n,j-2))) )/2.d0;
        if (delt/phi(n,j))*(sigmar(n,j+1)-sigmar(n,j-1))+delt*beta(n,j)*(d-1.d0) == 0
        U(n+1,j) = U(n-1,j);
        else
        U(n+1,j) = U(n-1,j)+(delt/phi(n,j))*(sigmar(n,j+1)-sigmar(n,j-1))+delt*beta(n,j)*(d-1.d0);            
        end
       elseif (ibc(j) == -4) 
        U(n+1,j) = 0.d0;
%c       if ( U(n+1,j) .lt. 1.d-5) U(n+1,j) = 0.d0
       elseif (ibc(j) == 4)
        U(n+1,j) = 0.d0;
%c       if ( U(n+1,j) .lt. 1.d-5) U(n+1,j) = 0.d0
       else
        fprintf('Momentum ERROR!!! \n');
        fprintf('%2.0f %2.0f',j,ibc(j));
       end
%c
      end
      end    
      
        
    else
        deltat=deltat/4;
        t(n+1) = t(n) + deltat/2.d0; %Higher
        t(n+2) = t(n) + deltat;  %Higher
                  delt = t(n+1)-t(n-1);   %this is different than deltat, deltat is the true timestep
%c	
      for j=0+1:2:jj+1
%c

      if (ibc(j) ~= 9) %13
       if (ibc(j) == -1)
                              %%%%%%%%<Problem HERE NSH
        phi(n,j) = (rho(j+1)*((r(n,j+2)-r(n,j))/V(n,j+1)))/2.d0;
        beta(n,j) = (sigmar(n,j+1)-sigmao(n,j+1))*V(n,j+1)/(r(n,j+1)*rho(j+1));
        if (delt/phi(n,j))*(sigmar(n,j+1)+0.000000000000)+delt*beta(n,j)*(d-1.d0) == 0
        U(n+1,j) = U(n-1,j);
        else
        U(n+1,j) = U(n-1,j)+(delt/phi(n,j))*(sigmar(n,j+1)+0.000000000000)+delt*beta(n,j)*(d-1.d0);
        end
%c       print *,'BC thing',bc(-1,1)
% What the heck is bc(-1,1)
       elseif (ibc(j) == 1)
        phi(n,j) =  rho(j-1)*((r(n,j)-r(n,j-2))/V(n,j-1))/2.d0;
        beta(n,j) = (sigmar(n,j-1)-sigmao(n,j-1))*V(n,j-1)/(r(n,j-1)*rho(j-1));
        if (delt/phi(n,j))*(-bc( 1,1)-sigmar(n,j-1))+delt*beta(n,j)*(d-1.d0) == 0
        U(n+1,j) = U(n-1,j);
        else
        U(n+1,j) = U(n-1,j)+(delt/phi(n,j))*(-bc( 1,1)-sigmar(n,j-1))+delt*beta(n,j)*(d-1.d0);            
        end
%c       U(n+1,j) = 0.d0
%c
%c this is the new way which treats the voids like inner or outer bc's
       elseif (ibc(j) == -2 || ibc(j) == -3) 
        phi(n,j) = (rho(j+1)*((r(n,j+2)-r(n,j))/V(n,j+1)))/2.d0;
        beta(n,j) = (sigmar(n,j+1)-sigmao(n,j+1))*V(n,j+1)/(r(n,j+1)*rho(j+1));
        U(n+1,j) = U(n-1,j)+(delt/phi(n,j))*(sigmar(n,j+1)+pvoid)+delt*beta(n,j)*(d-1.d0);
       elseif (ibc(j) == 2) 
        phi(n,j) =  rho(j-1)*((r(n,j)-r(n,j-2))/V(n,j-1))/2.d0;
        beta(n,j) = (sigmar(n,j-1)-sigmao(n,j-1))*V(n,j-1)/(r(n,j-1)*rho(j-1));
        U(n+1,j) = U(n-1,j)+(delt/phi(n,j))*(pvoid-sigmar(n,j-1))+delt*beta(n,j)*(d-1.d0);
       elseif (ibc(j) == 0)
        phi(n,j) = (0.5d0)*(rho(j+1)*((r(n,j+2)-r(n,j))/V(n,j+1))+rho(j-1)*((r(n,j)-r(n,j-2))/V(n,j-1)));
        beta(n,j)=( (sigmar(n,j+1)-sigmao(n,j+1))*(V(n,j+1)/rho(j+1))/(0.5d0*(r(n,j+2)+r(n,j  ))) +(sigmar(n,j-1)-sigmao(n,j-1))*(V(n,j-1)/rho(j-1))/(0.5d0*(r(n,j  )+r(n,j-2))) )/2.d0;
        if (delt/phi(n,j))*(sigmar(n,j+1)-sigmar(n,j-1))+delt*beta(n,j)*(d-1.d0) == 0
        U(n+1,j) = U(n-1,j);
        else
        U(n+1,j) = U(n-1,j)+(delt/phi(n,j))*(sigmar(n,j+1)-sigmar(n,j-1))+delt*beta(n,j)*(d-1.d0);            
        end
       elseif (ibc(j) == -4) 
        U(n+1,j) = 0.d0;
%c       if ( U(n+1,j) .lt. 1.d-5) U(n+1,j) = 0.d0
       elseif (ibc(j) == 4)
        U(n+1,j) = 0.d0;
%c       if ( U(n+1,j) .lt. 1.d-5) U(n+1,j) = 0.d0
       else
        fprintf('Momentum ERROR!!! \n');
        fprintf('%2.0f %2.0f',j,ibc(j));
       end
%c
      end
      end    
    end

SaveTimeStep(timesteps)=deltat;
fprintf('TimeStep = %4.4d\n',deltat)
end
flag=0;

%%
      
%ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
%c     3.  Position  - within j
%c
for j=0+1:2:jj+1
    if (ibc(j) ~= 9)
        r(n+2,j) =  r(n,j)+U(n+1,j)*deltat;
        r(n+1,j) = (r(n,j)+r(n+2,j))/2.d0;
    end
end

%%
% c debug block
%          if 1==1; %(idebug == 1 && ncount >= ndebug)
%          j=jdebug
%          fprintf('j= %2.4f n= %2.0f ncount= %2.4f \n',j,n,ncount);
%          fprintf('time %2.4f %2.4f %2.4f %2.4f \n',t(n-1),t(n),t(n+1),t(n+2));
%          fprintf('deltat %2.4f %2.4f \n',deltat,delt);
%          fprintf('vel %2.4f %2.4f \n',U(n,j),U(n+1,j));
%          fprintf('r %2.4f %2.4f %2.4f \n',r(n,j-2),r(n,j),r(n,j+2));
%          fprintf('V %2.4f %2.4f %2.4f \n',V(n,j+1),V(n,j),V(n,j-1));
%          fprintf('rho %2.4f %2.4f %2.4f \n',rho(j-1),rho(j),rho(j+1));
%          fprintf('sigmar %2.4f %2.4f \n',sigmar(n,j+1),sigmar(n,j-1));
%          fprintf('sigmao %2.4f %2.4f \n',sigmao(n,j-1),sigmao(n,j+1));
%          fprintf('4 %2.4f %2.4f %2.4f \n',phi(n,j),beta(n,j),U(n-1,j));
%          fprintf('5 %2.4f \n',(deltat/phi(n,j))*(sigmar(n,j+1)-0.d0));
% %         read(*,*) foo
%          end
% c	print *,'Line 537'
%% 
% c
% cccccc     Begin the connection of nodes due to contact check
% c
%% ISSUE HERE WITH IBC
       for j=2+1:2:jj-2+1
         if ( ibc(j) == -3 )
           fprintf('Joining Nodes %2.0f and %2.4f at n= %2.4f \n',j-2,j,n);
           fprintf('Before Joining: \n');
           fprintf('The two nodes (VOID (left)and J(right)) at: \n');
           fprintf('left: r(void,n-2)= %2.4f U= %2.4f \n',r(n,j-4),U(n-1,j-4));
           fprintf('VOID: r(void,n-2)=%2.4f U= %2.4f \n',r(n,j-2),U(n-1,j-2));
           fprintf('J:    r(j   ,n-2)=%2.4f U= %2.4f \n',r(n,j)  ,U(n-1,j));
           fprintf('right:r(j   ,n-2)=%2.4f U= %2.4f \n',r(n,j+2),U(n-1,j+2));
           fprintf('then stepped to: \n');
           fprintf('left: r(-2  )=%2.4f U= %2.4f  \n',r(n+2,j-4),U(n+1,j-4));
           fprintf('VOID: r(void)=%2.4f U= %2.4f  \n',r(n+2,j-2),U(n+1,j-2));
           fprintf('J:    r(j   )=%2.4f U= %2.4f  \n',r(n+2,j)  ,U(n+1,j));
           fprintf('right:r(+2  )=%2.4f U= %2.4f  \n',r(n+2,j+2),U(n+1,j+2));
           fprintf('where they were joined \n');
                    ibc(j)   = 0;
        Usave = (m(j+1)*U(n+1,j)+m(j-3)*U(n+1,j-2))/(m(j-3)+m(j+1));
        fprintf('Usave %2.4f %2.4f %2.4f %2.4f %2.4f %2.4f %2.4f \n',Usave,m(j+1),U(n+1,j),m(j-3),U(n+1,j-2),m(j-3),m(j+1));
           for jjj = j-2:jj-2
            for nz = n:n+2
               ibc(    jjj  )  =    ibc(    jjj+2);
               ibc(    jjj+1)  =    ibc(    jjj+3);
              ieos(1,  jjj+1)  =   ieos(1,  jjj+3);
              ieos(2,  jjj+1)  =   ieos(2,  jjj+3);
                 m(    jjj  )  =      m(    jjj+2);
               rho(    jjj+1)  =    rho(    jjj+3);
 
                 r(nz ,jjj  )  =     r(nz ,jjj+2);
                 U(nz ,jjj  )  =     U(nz ,jjj+2);
               phi(nz ,jjj  ) =    phi(nz ,jjj+2);
              beta(nz ,jjj  ) =   beta(nz ,jjj+2);
            sigmar(nz ,jjj+1) = sigmar(nz ,jjj+3);
            sigmao(nz ,jjj+1) = sigmao(nz ,jjj+3);
                 V(nz ,jjj+1) =      V(nz ,jjj+3);
                s1(nz ,jjj+1) =     s1(nz ,jjj+3);
                s2(nz ,jjj+1) =     s2(nz ,jjj+3);
                s3(nz ,jjj+1) =     s3(nz ,jjj+3);
                 E(nz ,jjj+1) =      E(nz ,jjj+3);
                 Y(    jjj+1) =      Y(    jjj+3);
             pfrac(    jjj  ) =  pfrac(    jjj+2);  %!node value
             pfrac(    jjj+1) =  pfrac(    jjj+3);  %!Cell value
            end
            end
                 U(n+1 ,j  -2) = Usave;
              pfrac(    j  -2) = 1.d-2;
% c                Y(     j  -1) = 0.d0   !ie the materials are not welded after impact
% c               r(n+2 ,j  -2) = r(n+1,j-2)+Usave*(t(n+2)-t(n+1))
% c          if (abs(r(0,ipoint+j1(jjj-1))-xstart(jjj)) .lt. 1.e-5) then
% c           Print *,'Warning:  Joined moved to same spot'
% c          endif
               fprintf('After Joining:');
               fprintf('left: r(j-2)=%2.4f U= %2.4f\n',r(n+2,j-4),U(n+1,j-4));
               fprintf('New:  r(j) = %2.4f U= %2.4f\n',r(n+2,j-2),Usave);
               fprintf('right:r(j+2)=%2.4f U= %2.4f\n',r(n+2,j  ),U(n+1,j));
%                read(*,*) foo
%       for jzz=1,jj-2,2
%        write(*,'(I4,5I2,3f7.3,2e9.3,5f7.3)')
%      &  jzz,ieos(1,jzz),ieos(2,jzz),
%      &  ibc(jzz-1),ibc(jzz),ibc(jzz+1),
%      &  r(n+2,jzz-1),r(n+2,jzz),r(n+2,jzz+2),
%      &  U(n+1,jzz-1),U(n+1,jzz+2),U(n+1,jzz+1)
% c     &  ,v(n+2,jzz+1),m(jzz+1),rho(jzz+1),y(jzz+1),pfrac(jzz+2)
%       end
       %read(*,*) foo
         end
       end

       
% c	print *,'line 608'
% cccccccccccccccccccccc End of contact ccccccccccccccccccccccccccccc
% c     Node passes it's neighbor check
% c
       for j=2+1:2:jj-2+1
         if ( ibc(j) ~= 9 && r(n+2,j) < r(n+2,j-2) )
           fprintf('Contact Check Failed %2.4f %2.4f\n',t(n),j);
           fprintf('Try lowering time step\n');
           fprintf('Void %2.4f %2.4f\n',r(n+2,j-2),U(n+1,j-2));
           fprintf('j  %2.4f %2.4f\n',r(n+2,j)  ,U(n+1,j));
           fprintf('deltat %2.4f\n',deltat);
           fprintf('%2.4f %2.4f\n',bc(-1,1),bc(1,1));
%           read(*,*) foo
         end
         end

 %%     
% ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
% c     3.  Relative VOLUME   - within j
% c
% c      deltat     =  t(n+2)-t(n)
% c      print *,'line 626'
% c	read(*,*) foo
       for j=0+1:2:jj-2+1
       if (ibc(j+1) == 0)
% c
         r(n+1,j+1) = ( r(n+1,j)+r(n+1,j+2) ) /2.d0;
         r(n+2,j+1) = ( r(n+2,j)+r(n+2,j+2) ) /2.d0;
% c
         V(n+2,j+1)=rho(j+1)*((r(n+2,j+2)^d-r(n+2,j)^d)/d)/m(j+1);
         V(n+1,j+1)=rho(j+1)*((r(n+1,j+2)^d-r(n+1,j)^d)/d)/m(j+1);
% c
         if ( V(n+2,j+1) == 0.)
         fprintf('------------- Zero volume error! #1-----------------');
         fprintf('volume %2.0f %2.0f %2.0f %2.0f %2.0f',n,j+1,r(n+2,j+2),r(n+2,j),V(n+2,j+1));
         %read(*,*) foo
         end
         if ( V(n+1,j+1) == 0.)
         fprintf('------------- Zero volume error! #2-----------------');
         fprintf('volume %2.0f %2.0f %2.0f %2.0f %2.0f',n,j+1,r(n+1,j+2),r(n+1,j),V(n+1,j+1));
         %read(*,*) foo
         end
%ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
%c     4. VELOCITY STRAINS  - finish j
%c
         epsi1(n+1,j+1) = (U(n+1,j+2)-U(n+1,j))/(r(n+1,j+2)-r(n+1,j));
         epdt1(n+1,j+1) = (epsi1(n+1,j+1)-epsi1(n-1,j+1))/deltat; %!this is 1st order accurate
        if (d == 1.)
         epsi2(n+1,j+1) = 0.d0;
        else
         epsi2(n+1,j+1) = (U(n+1,j+2)+U(n+1,j))/(r(n+1,j+2)+r(n+1,j));
        end
       end

       end  %! end j loop for volume

      %%
% c--------------------  AMR Feature  -------------------------------
% 
%         if (idebug .eq. 1 .and. n .ge. ndebug) then
%         j=jdebug
%         print *,'r',j,r(n+2,j+1),m(j+1)
%         print *,'epsil',j,V(n+1,j+1),V(n+2,j+1),
%      &                    epsi1(n+1,j+1),epsi2(n+1,j+1)
%         read(*,*) foo
%         endif

     %% 
%ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
%c     5a. STRESSES  - start j
%c ibc(241) means j =239
%c      print *,'line 679'
      for j=0+1:2:jj-2+1
      if (ibc(j+1) == 0)
%c       deltat = t(n+2)-t(n)
       if (ieos(2,j+1) ~= 3)
        s1(n+2,j+1) = s1(n,j+1) + 2.d0*eos(ieos(1,j+1),6)* ( epsi1(n+1,j+1)*deltat- (V(n+2,j+1)-V(n,j+1))/(3.d0*V(n+1,j+1)) );
        s2(n+2,j+1) = s2(n,j+1) + 2.d0*eos(ieos(1,j+1),6)* ( epsi2(n+1,j+1)*deltat- (V(n+2,j+1)-V(n,j+1))/(3.d0*V(n+1,j+1)) );
        s3(n+2,j+1) = -(s1(n+2,j+1)+s2(n+2,j+1));
%c
        s1(n+1,j+1) = ( s1(n+2,j+1) + s1(n,j+1) )/2.d0;
        s2(n+1,j+1) = ( s2(n+2,j+1) + s2(n,j+1) )/2.d0;
        s3(n+1,j+1) = ( s3(n+2,j+1) + s3(n,j+1) )/2.d0;
       else     %!Gamma law ideal gas Newtonian Stress
        s1(n+2,j+1) = 4.d0*1.8d-10*epsi1(n+1,j+1)/3.d0-1.387e-6* (V(n+2,j+1)-V(n,j+1))/V(n+1,j+1);
%c        s1(n+2,j+1) = s1(n,j+1)+4.d0*1.8d-10*epdt1(n+1,j+1)*deltat/3.d0
%c     &             -1.387e-6* (V(n+2,j+1)-V(n,j+1))/(3.d0*V(n+1,j+1)) 
%c
%c        print *,epdt1(n+1,j+1),epsi1(n+1,j+1),epsi1(n-1,j+1)
%c
        s1(n+1,j+1) = ( s1(n+2,j+1) + s1(n,j+1) )/2.d0;
       end %! iEOS check
%c


      end
      end
        if (idebug == 1 && n >= ndebug)
      j=jdebug;
      fprintf('1 %2.0f %2.0f %2.0f %2.0f %2.0f',j,eos(ieos(1,j+1),6),V(n+2,j+1),V(n,j+1),V(n+1,j+1));
      fprintf('s %2.0f %2.0f %2.0f %2.0f',s1(n+2,j+1),s2(n+2,j+1),s3(n+2,j+1),epsi1(n+1,j+1));
        end

      %%
% ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
% c     6. VON MISES YIELD CONDITION - start j
% c
for j=0+1:2:jj-2+1
    if (ibc(j+1) == 0)
        if (ieos(2,j+1) ~= 3)
            % c calculate the deviatoric strain at n+1 and compare it to the yield strength
            % c
            k(n+1,j+1) =(s1(n+1,j+1)^2 + s2(n+1,j+1)^2 + s3(n+1,j+1)^2)-(2.d0/3.d0)*Y(j+1)^2;
            % c       print *,'Yield check',Y(j+1),k(n+2,j+1)
            if (k(n+1,j+1) > 0.)
                %fprintf('material yeilded 1\n')
                % c       Print *,'material yielded',j,k(n+2,j+1)
                % c        read(*,*) foo
                % c        print *,s1(n+2,j+1),s2(n+2,j+1),s3(n+2,j+1)
                xx = sqrt(s1(n+1,j+1)^2 + s2(n+1,j+1)^2 + s3(n+1,j+1)^2);
                xx = sqrt(2.d0/3.d0)*Y(j+1)/xx;
                s1(n+1,j+1) =   xx*s1(n+1,j+1);
                s2(n+1,j+1) =   xx*s2(n+1,j+1);
                s3(n+1,j+1) =   xx*s3(n+1,j+1);
            end
            % c calculate the deviatoric strain at n+2 and compare it to the yield strength
            
            k(n+2,j+1) =(s1(n+2,j+1)^2 + s2(n+2,j+1)^2 + s3(n+2,j+1)^2)-(2.d0/3.d0)*Y(j+1)^2;
            % c       print *,'Yield check',Y(j+1),k(n+2,j+1)
            if (k(n+2,j+1) > 0.)
               % fprintf('material yeilded 2\n')
                % c       Print *,'material yielded',j,k(n+2,j+1)
                % c        read(*,*) foo
                % c        print *,s1(n+2,j+1),s2(n+2,j+1),s3(n+2,j+1)
                xx = sqrt(s1(n+2,j+1)^2 + s2(n+2,j+1)^2 + s3(n+2,j+1)^2);
                xx = sqrt(2.d0/3.d0)*Y(j+1)/xx;
                s1(n+2,j+1) =   xx*s1(n+2,j+1);
                s2(n+2,j+1) =   xx*s2(n+2,j+1);
                s3(n+2,j+1) =   xx*s3(n+2,j+1);
                % c         s1(n+1,j+1) = ( s1(n+2,j+1) + s1(n,j+1) )/2.d0
                % c         s2(n+1,j+1) = ( s2(n+2,j+1) + s2(n,j+1) )/2.d0
                % c         s3(n+1,j+1) = ( s3(n+2,j+1) + s3(n,j+1) )/2.d0
                % c        Print *,'material yielded',j,k(n+2,j+1),Y(j+1)
            else %!iEOS = 3, i.e. gamma law gas
            end
            % c
            % c        print *,s1(n+2,j+1),s2(n+2,j+1),s3(n+2,j+1)
            if(idebug == 1 )
                fprintf('Yield Check %2.0f %2.0f',Y(j+1),k(n+2,j+1))
                %           read(*,*) foo
            end %! debug if           
        end  %!ibc(j+1) .eq. 0
%         
%         if s1(n+1,j+1) < 1e-8 &&  s2(n+1,j+1) < 1e-8 &&   s3(n+1,j+1) < 1e-8
%             s1(n+1,j+1) =   0;
%             s2(n+1,j+1) =   0;
%             s3(n+1,j+1) =   0;
%         end
%         if s1(n+2,j+1) < 1e-8 &&  s2(n+2,j+1) < 1e-8 &&   s3(n+2,j+1) < 1e-8
%             s1(n+2,j+1)=0;
%             s2(n+2,j+1)=0;
%             s3(n+2,j+1)=0;
%         end


        %ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        %c     7. ARTIFICIAL VISCOSITY - within j
        %c
        if (1 == 1)
            %c       print *,j,u(n+1,j+2),u(n+1,j)
            if( U(n+1,j+2) < U(n+1,j) && V(n+2,j+1)-V(n,j+1) < 0. )
                xx = rho(j+1)/V(n+1,j+1);
                if (P(n,j+1) < 0.)
                    %c         print *,'Negative Pressure!'
                    a  = sqrt(-P(n,j+1)/xx);
                else
                    a  = sqrt( P(n,j+1)/xx);
                end
                if (abs(P(n,j+1)) < 0.002d0)
                    a = 0.d0;
                end
                %c       print *,'Sound Speed',j,P(n,j+1),xx,a
                q(n+1,j+1) = Co*Co*xx*    (U(n+1,j+2)-U(n+1,j))^2+ CL*a *xx*abs(U(n+1,j+2)-U(n+1,j));
            else
                q(n+1,j+1)=0.d0;
            end
        end
    end
end %!j counter

%if (k(n+1,:) > 0.)
%fprintf('material yeilded\n')
%end
        
        
        %%
        
        
        
        
        
        
        
        
        
        
        if (idebug == 1 && n >= ndebug)
%         j=jdebug;
%       print *,'Q',j,q(n+1,j+1),xx,(U(n+1,j+2)-U(n+1,j)),a,rho(j+1),v(n,j+1),P(n,j+1)
        end
%             if(idebug .eq.1 .and. n .ge. ndebug) read(*,*) foo

%cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
%c     8.  ENERGY  - start j
%c     iEOS(2,imat)=1 - Mie Grunisen
%c     iEOS(2,imat)=2 - Gamma law ideal gas
%c     iEOS(2,imat)=3 - Gamma law ideal gas with Newtonian Stress and heat transfer
%c     iEOS(2,imat)=4 - snow plow model with KO inputs for Hugoniot
%c     iEOS(2,imat)=5 - snow plow model with anamolous hugoniot
%c     iEOS(2,imat)=6 - P-alpha Model
%c
      for j=0+1:2:jj-2+1
      if (ibc(j+1) == 0)
%c
       if (ieos(2,j+1) == 1)                  %!Mie Grunisen
        gamma0 = eos(ieos(1,j+1),4);
        k1  = rho(j+1)*eos(ieos(1,j+1),1)^2;
        k2  =  k1*(2.d0*eos(ieos(1,j+1),2)-gamma0/2.d0);
        k3  =  k1*(3.d0*eos(ieos(1,j+1),2)-gamma0)*eos(ieos(1,j+1),2);
        qbar = (q(n+1,j+1)+q(n-1,j+1))/2.d0;
        deltaZ = V(n+1,j+1)*(s1(n+1,j+1)*epsi1(n+1,j+1)+(d-1.d0)*s2(n+1,j+1)*epsi2(n+1,j+1))*deltat;
        strain = 1.d0 - V(n+2,j+1);
        xb = gamma0; %!gamma
        if(strain < 0.d0)
         xa = k1*strain + k3*strain^3;
        else
         xa = k1*strain + k2*strain^2 + k3*strain^3;
        end
        E(n+2,j+1)=( E(n,j+1)-((xa+P(n,j+1))/2.d0+ qbar)*(V(n+2,j+1)-V(n,j+1))+ deltaZ )/ (1.d0 + xb*(V(n+2,j+1)-V(n,j+1))/2.d0);
        E(n+1,j+1)= ( E(n+2,j+1) + E(n,j+1) )/2.d0;
%c      if (n .eq. 33) then
%c        write (*,'(a2,I3,9e9.2)')
%c     &   'E ',j,E(n+2,j+1),E(n+2,j+1),xa,P(n,j+1),qbar,deltaZ
%c     &   ,V(n+2,j+1),V(n,j+1)
%cc      read(*,*) foo
%c      endif
%c
       elseif (ieos(2,j+1) == 2)              %!Gamma law ideal gas
        qbar   = (q(n+1,j+1)+q(n-1,j+1))/2.d0;
        delt   = t(n+1)-t(n);
        deltaZ = V(n+1,j+1)*(s1(n+1,j+1)*epsi1(n+1,j+1)+(d-1.d0)*s2(n+1,j+1)*epsi2(n+1,j+1))*delt;
%c
        xa   = 0.d0;
        xb   = (eos(ieos(1,j+1),4)-1.d0)/V(n+2,j+1);  %!gamma
        E(n+2,j+1)= ( E(n,j+1)-((xa+P(n,j+1))/2.d0+ qbar)*(V(n+2,j+1)-V(n,j+1)) + deltaZ )/(1.d0+xb*(V(n+2,j+1)-V(n,j+1))/2.d0);
        E(n+1,j+1)= ( E(n+2,j+1) + E(n,j+1) )/2.d0;
%c      if (j .eq. 0 ) then
%c        if (n .lt.2. .or. n .eq. 201)    E(n+2,j+1)=0.d0
%c        xx=V(n+2,j+1)-V(n,j+1)
%c        print *,E(n,j+1)   /(1.*(1.+xb*(V(n+2,j+1)-V(n,j+1))/2.))
%c        print *,xa*xx      /(2.*(1.+xb*(V(n+2,j+1)-V(n,j+1))/2.))
%c        print *,P(n,j+1)*xx/(2.*(1.+xb*(V(n+2,j+1)-V(n,j+1))/2.))
%c        print *,qbar*xx    /(1.*(1.+xb*(V(n+2,j+1)-V(n,j+1))/2.))
%c        print *,deltaZ  /(1.*(1.+xb*(V(n+2,j+1)-V(n,j+1))/2.))
%c      endif
%c
       elseif (iEOS(2,j+1) == 3)                %!Gamma law ideal gas     ???
        qbar = (q(n+1,j+1)+q(n-1,j+1))/2.d0;
        delt   = t(n+1)-t(n);
        deltaZ = V(n+1,j+1)*(s1(n+1,j+1)*epsi1(n+1,j+1)+(d-1.0d0)*s2(n+1,j+1)*epsi2(n+1,j+1))*delt;
%c
        xx   = rho(j+1)/V(n+2,j+1);
        xa   = 0.d0;
        xb   = (eos(ieos(1,j+1),4)-1.d0)/V(n+2,j+1);
        E(n+2,j+1)= ( E(n,j+1)-((xa+P(n,j+1))/2.d0+ qbar)*(V(n+2,j+1)-V(n,j+1))+ deltaZ )/(1.d0+xb*(V(n+2,j+1)-V(n,j+1))/2.d0);
%c
       elseif (iEOS(2,j+1) == 4)               %!Snow Plow
        qbar = (q(n+1,j+1)+q(n-1,j+1))/2.d0;
        deltaZ = V(n+1,j+1)*(s1(n+1,j+1)*epsi1(n+1,j+1)+(d-1.d0)*s2(n+1,j+1)*epsi2(n+1,j+1))*deltat;
        xb = eos(ieos(1,j+1),4); %!gamma
        v0 = eos(ieos(1,j+1),12);
        xx = 1.d0 - (V(n+2,j+1)/rho(j+1))/V0;  %!this makes the compression relative to the compacted material hugoniot
        if ( v(n+2,j+1)/rho(j+1) <= v0) 
            icompact(n+2,j+1) = 1
        end
        if(     icompact(n+2,j+1) == 1 && xx < 0.d0)
         xa = eos(ieos(1,j+1),1)*xx+eos(ieos(1,j+1),3)*xx^3;
        elseif (icompact(n+2,j+1) == 1 && xx >= 0.d0)
         xa = eos(ieos(1,j+1),1)*xx+eos(ieos(1,j+1),2)*xx^2+eos(ieos(1,j+1),3)*xx^3;
        elseif( icompact(n+2,j+1) == 0)
         xa = 0.d0;   %!snow plow model - this zeros out pressure until threshold is reached
        else
         fprintf('compaction error in energy')
         fprintf('%2.0f %2.0f',icompact(n+2,j+1),xx)
        end
        E(n+2,j+1)=( E(n,j+1)-((xa+P(n,j+1))/2.d0+ qbar)*(V(n+2,j+1)-V(n,j+1))+ deltaZ )/ (1.d0 + xb*(V(n+2,j+1)-V(n,j+1))/2.d0);
        E(n+1,j+1)= ( E(n+2,j+1) + E(n,j+1) )/2.d0;
%c
       elseif (iEOS(2,j+1) == 5)               %!Snow Plow with anaomolus hugoniot
        qbar = (q(n+1,j+1)+q(n-1,j+1))/2.d0;
        deltaZ = V(n+1,j+1)*(s1(n+1,j+1)*epsi1(n+1,j+1)+(d-1.d0)*s2(n+1,j+1)*epsi2(n+1,j+1))*deltat;
        stemp = 0.1048d0;      %!slope
        ctemp = 0.5124d0;      %!bulk sound speed
        gtemp = 0.9d0;         %!gamma
        v0 = eos(ieos(1,j+1),8);
        V00= 1.d0/rho(j+1);
        vv = V(n+2,j+1)/rho(j+1);    %!v_local
        if ( v(n+2,j+1)/rho(j+1) <= v0) 
            icompact(n+2,j+1) = 1
        end
        if (icompact(n+2,j+1) == 1)  %!icompact is the flag, once compact always compact
         xa = ((2.d0*vv-gtemp*(v0 -vv))*ctemp*ctemp*(v0-vv))/((2.d0*vv-gtemp*(v00-vv))*((v0-stemp*(v0-vv))^2) );     %!porous hugoniot, see meyers pg 141
         xa = dabs(xa); %!not sure if i need this but it was add so that re-denstended materials don't have negative pressure
        else
         xa = 0.d0;
        end
        E(n+2,j+1)=( E(n,j+1)-((xa+P(n,j+1))/2.d0+ qbar)*(V(n+2,j+1)-V(n,j+1))+ deltaZ )/ (1.d0 + xb*(V(n+2,j+1)-V(n,j+1))/2.d0);
        E(n+1,j+1)= ( E(n+2,j+1) + E(n,j+1) )/2.d0;
%c
       elseif (iEOS(2,j+1) == 6)                %!p-alpha
        qbar = (q(n+1,j+1)+q(n-1,j+1))/2.d0;
        deltaZ = V(n+1,j+1)*(s1(n+1,j+1)*epsi1(n+1,j+1)+(d-1.d0)*s2(n+1,j+1)*epsi2(n+1,j+1))*deltat;
        xx = 1.d0 - V(n+2,j+1);
        xb = eos(ieos(1,j+1),4); %!gamma
        if(xx < 0.)
        xa = eos(ieos(1,j+1),1)*xx+0.d0*xx^2+eos(ieos(1,j+1),3)*xx^3;
        else
        xa = eos(ieos(1,j+1),1)*xx+eos(ieos(1,j+1),2)*xx^2+eos(ieos(1,j+1),3)*xx^3;
        end
        E(n+2,j+1)=( E(n,j+1)-((xa+P(n,j+1))/2.d0+ qbar)*(V(n+2,j+1)-V(n,j+1))+ deltaZ )/ (1.d0 + xb*(V(n+2,j+1)-V(n,j+1))/2.d0);
        E(n+1,j+1)= ( E(n+2,j+1) + E(n,j+1) )/2.d0;
%c
%c
       else
        fprintf('EOS error!');
        fprintf('%2.0f %2.0f %2.0f',n,j,ieos(2,j+1));
        %read(*,*) foo
        end
        end
%c
        end
%c
      if(idebug == 1  && n >= ndebug)
        j=jdebug;
        fprintf('Energy %2.0f %2.0f %2.0f %2.0f %2.0f',j,E(n+2,j+1),E(n,j+1),xx,(2.d0*(1.d0+xb*(V(n+2,j+1)-V(n,j+1))/2.d0)));
      end
%      if(idebug .eq.1  .and. n .gt. ndebug) read(*,*) foo

        
        %%
        
% cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
% c     5b PRESSURE - start j
% c
% c      print *,'line 935'
       for j = 0+1:2:jj-2+1
       if (ibc(j+1) == 0)
        if (ieos(2,j+1) == 1)   %!     Mie Gruneisen EOS as formulated in Wilkins
         gamma0 = eos(ieos(1,j+1),4);
         k1  = rho(j+1)*eos(ieos(1,j+1),1)^2;
         k2  =  k1*(2.d0*eos(ieos(1,j+1),2)-gamma0/2.d0);
         k3  =  k1*(3.d0*eos(ieos(1,j+1),2)-gamma0)*eos(ieos(1,j+1),2);
% c        if (j .eq. 2840) then
% c        print *,eos(ieos(1,j+1),1),eos(ieos(1,j+1),2),eos(ieos(1,j+1),3)
% c        read(*,*) foo
% c        endif
         strain      = 1.d0 - V(n+1,j+1);
         if (strain < 0.d0)
          P(n+1,j+1) = k1*strain + k3*strain^3 + gamma0*E(n+1,j+1);
         else
          P(n+1,j+1) = k1*strain + k2*strain^2 + k3*strain^3 + gamma0*E(n+1,j+1);
         end
         strain      = 1.d0 - V(n+2,j+1);
         if (strain < 0.d0)
          P(n+2,j+1) = k1*strain + k3*strain^3 + gamma0*E(n+2,j+1);
         else
          P(n+2,j+1) = k1*strain + k2*strain^2 + k3*strain^3+ gamma0*E(n+2,j+1);
         end
% c        print *,'Mie Grun EOS 1',P(n+1,j+1),P(n+2,j+1)
% c        read(*,*) foo
% c       if (icontact .eq. 1) then
% c       if (P(n+2,j+1) .lt. 0.) then
% c       print *,'Neg Pres',n+1,j+1,xx,v(n+1,j+1),E(n+1,j+1),p(n+1,j+1)
% c       print *,'Neg Pres',n+2,j+1,xx,v(n+2,j+1),E(n+2,j+1),p(n+2,j+1)
% cc       read(*,*) foo
% c       endif
% c       endif
        elseif (iEOS(2,j+1) == 2)           %! Gamma Law (perfect gas)
          P(n+2,j+1) = (eos(ieos(1,j+1),4)-1.d0)*E(n+2,j+1)/V(n+2,j+1);
          P(n+1,j+1) = (P(n,j+1)+P(n+2,j+1))/2.d0;
        elseif (ieos(2,j+1) == 3)            %! Gamma Law
          P(n+2,j+1) = (eos(ieos(1,j+1),4)-1.d0)*E(n+2,j+1)/V(n+2,j+1);
          P(n+1,j+1) = (P(n,j+1)+P(n+2,j+1))/2.d0;
% c        xa          = (u(n+2,j)+u(n+2,j+2))/2.d0
% c        xx          = rho(j+1)/V(n+2,j+1)
% c         P(n+2,j+1) = P(n,j+1)
% c     &     +(eos(ieos(1,j+1),4)+1.d0)*xx*xa*xa/2.d0
% c       print *,'pres',rho_0,gamma+1,xa,xx
        elseif (ieos(2,j+1) == 4)             %!Snow Plow
         v0  = eos(ieos(1,j+1),12);
         V00 = 1.d0/rho(j+1);
         vv  = v(n+1,j+1)/rho(j+1);
         if (v0 == 0.)
          fprintf('v0 can not equal zero');
          fprintf('%2.0f %2.0f %2.0f',j,ieos(1,j+1),eos(ieos(1,j+1),12));
          %read(*,*) foo
         end
% c
         if (vv <= V0) 
             icompact(n+1,j+1) = 1
         end
          xx         = 1.d0 - (V(n+1,j+1)/rho(j+1))/V0;   %!this makes the compression relative to the compacted material hugoniot
         if (icompact(n+1,j+1) == 1 && xx < 0.)
           P(n+1,j+1) = eos(ieos(1,j+1),1)*xx+eos(ieos(1,j+1),3)*xx^3+eos(ieos(1,j+1),4)*E(n+1,j+1);
         elseif(icompact(n+1,j+1) == 1 && xx >= 0.d0)
           P(n+1,j+1) = eos(ieos(1,j+1),1)*xx+eos(ieos(1,j+1),2)*xx^2+eos(ieos(1,j+1),3)*xx^3+eos(ieos(1,j+1),4)*E(n+1,j+1);
         elseif (icompact(n+1,j+1) == 0)
           P(n+1,j+1) = 0.d0*eos(ieos(1,j+1),4)*E(n+1,j+1);
         else
          fprintf('pressure compaction error');
%          read(*,*) foo
         end
% c        if ( icompact(n+1,j+1) .eq. 1) then
% c        print *,t(n),j+1,icompact(n+1,j+1),vv,p(n+1,j+1)
% c        print *,eos(ieos(1,j+1),1)*xx,eos(ieos(1,j+1),3)*xx**3
% c     &         ,eos(ieos(1,j+1),4)*E(n+1,j+1)
% c        print *,xx,eos(ieos(1,j+1),1),eos(ieos(1,j+1),3)
% c     &         ,eos(ieos(1,j+1),4)
% c        read(*,*) foo
% c        else
% c        print *,t(n),j+1,icompact(n+1,j+1),vv,p(n+1,j+1)
% c        endif
% c
         if (v(n+2,j+1)/rho(j+1) <= V0) 
             icompact(n+2,j+1) = 1;
         end
          xx         = 1.d0 - (V(n+2,j+1)/rho(j+1))/V0     %!this makes the compression relative to the compacted material hugoniot
          if (icompact(n+2,j+1) == 1 && xx < 0.)
           P(n+2,j+1) = eos(ieos(1,j+1),1)*xx+eos(ieos(1,j+1),3)*xx^3+eos(ieos(1,j+1),4)*E(n+2,j+1)
          elseif (icompact(n+2,j+1) == 1 && xx >= 0.)
           P(n+2,j+1) = eos(ieos(1,j+1),1)*xx+eos(ieos(1,j+1),2)*xx^2+eos(ieos(1,j+1),3)*xx^3+eos(ieos(1,j+1),4)*E(n+2,j+1);
         elseif (icompact(n+2,j+1) == 0) %!v .gt. v0
           P(n+2,j+1) = 0.d0*eos(ieos(1,j+1),4)*E(n+2,j+1);
         else
          fprintf('pressure compaction error');
          %read(*,*) foo
          end
% c
        elseif (ieos(2,j+1) == 5)   %! Snow Plow with anaomolous Hugoniot
          V0 = eos(ieos(1,j+1),8);
          V00= 1.d0/rho(j+1); 
          vv = V(n+1,j+1)/rho(j+1);    %!v_local
          if (VV <= V0)
              icompact(n+1,j+1) = 1
          end
          if (icompact(n+1,j+1) == 1)
           xa = ((2.d0*vv-gtemp*(v0 -vv))*ctemp*ctemp*(v0-vv))/((2.d0*vv-gtemp*(v00-vv))*((v0-stemp*(v0-vv))^2) );    %!porous hugoniot, see meyers pg 141
           if (xa < 0.d0)
            fprintf('Compression less than v0 on anamolous hugoniot, n+1');
            xa=abs(xa);
           end
          else
           xa = 0.d0;
          end
 
          P(n+1,j+1) = xa + eos(ieos(1,j+1),4)*E(n+1,j+1);
 
          vv = V(n+2,j+1)/rho(j+1);%    !v_local
          if (VV <= V0)
              icompact(n+2,j+1) = 1
          end
          if (icompact(n+2,j+1) == 1)
          xa = ((2.d0*vv-gtemp*(v0 -vv))*ctemp*ctemp*(v0-vv))/((2.d0*vv-gtemp*(v00-vv))*((v0-stemp*(v0-vv))^2) );     %!porous hugoniot, see meyers pg 141
           if (xa < 0.d0)
            fprintf('Compression less than v0 on anamolous hugoniot, n+2');
            xa=abs(xa);
           end
          else
           xa = 0.d0;
          end
% c
        elseif (ieos(2,j+1) == 6)  %!p-alpha model
% c
% c
        elseif  (ieos(2,j+1) == 7)  %! Mie Gruneisen CTH like formulation
         P0 = 0.d0;
         E0 = 0.d0;
         T0 = 293.d0;
         strain        = 1.d0 - V(n+1,j+1);
         rho_local = rho(j+1)/V(n+1,j+1);
         up = (U(n+1,j)+U(n+1,j+2))/2.d0;
         Us = eos(ieos(1,j+1),1)+ eos(ieos(1,j+1),2)*up+ eos(ieos(1,j+1),3)/eos(ieos(1,j+1),1)*up^2;          %! local shock speed
         PH = P0 + rho_local*Us*up;
         EH = E0 + up*up/2.d0;
         TH = dexp( eos(ieos(1,j+1),4)*strain)*T0;
% c       E  = EH + eos(ieos(1,j+1),8)*(T-TH)
% c       P(n+2,j+1) = PH + eos(ieos(1,j+1),4)*rho_local*(
% c         if (Us*up .ne. Us*Us*strain) then
% c          print *,'EOS does not jive'
% c          read(*,*) foo
% c         endif
         xx         = 1.d0 - V(n+2,j+1);
         P(n+2,j+1) = eos(ieos(1,j+1),1)*xx+eos(ieos(1,j+1),2)*xx^2+eos(ieos(1,j+1),3)*xx^3+eos(ieos(1,j+1),4)*E(n+2,j+1)
        else    %!this is ieos(2,j+1) else
       fprintf('eos error in setup file')
         %read(*,*) foo
        end  %! IEOS
       end %!ibc
       end

       
       %%
% for j = 0+1:2:jj-2+1
%          if (P(n+1,j+1) < 1e-9)
%           P(n+1,j+1) = 0;
%          end
%          if (beta(n+1,j+1) < 1e-9)
%           beta(n+1,j+1) = 0;
%          end
%          
% end
       %%
       
% cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
% c This is the Temperature and Entropy Calculation
% c      print *,'line 1104'
       for j=0+1:2:jj-2+1
% c	print *,'line 1106',j,nc
% c	print *,'temp',temp(n+2c,j+1)
% c	print *,'E   ',e(n+2,j+1)
% c	print *,'rho ',rho(j+1)
% c	print *,'ieos',ieos(1,j+1)
% c	print *,'eos ',eos(ieos(1,j+1),8)
% c
        Temp(n+2,j+1)=E(n+2,j+1)/(rho(j+1)*eos(ieos(1,j+1)+1,8));
        % Double Check the line above!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
% c       print *,Temp(n+2,j+1),E(n+2,j+1),rho(j+1),v(n+2,j+1)
% c     &  ,eos(ieos(1,j+1),8)
% c      print *,'line 1109',j,n
        entropy(n+2,j+1)=entropy(n,j+1)+s1(n+2,j+1)/Temp(n+2,j+1);
% c       entropy(n+2,j+1)=entropy(n,j+1)+
% c       print *,entropy(n+2,j+1),entropy(n,j+1),
% c     &  s1(n+2,j+1),epsi1(n+2,j+1),Temp(n+2,j+1)
       end
% c


%%

% cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
% c This is the Pressure and Energy convergence check
% c     print *,'line 1118'
       for j=0+1:2:jj-2+1
        if (ieos(2,j+1) == 1)                    %!Mie Grunisen
         gamma0 = eos(ieos(1,j+1),4);
         k1  = rho(j+1)*eos(ieos(1,j+1),1)^2;
         k2  =  k1*(2.d0*eos(ieos(1,j+1),2)-gamma0/2.d0);
         k3  =  k1*(3.d0*eos(ieos(1,j+1),2)-gamma0)*eos(ieos(1,j+1),2);
         qbar = (q(n+1,j+1)+q(n-1,j+1))/2.d0;
         deltaZ = V(n+1,j+1)*(s1(n+1,j+1)*epsi1(n+1,j+1)+(d-1.d0)*s2(n+1,j+1)*epsi2(n+1,j+1))*deltat;
         strain = 1.d0 - V(n+2,j+1);
         xb = gamma0; %!gamma
         if(strain < 0.d0)
         xa = k1*strain + k3*strain^3;
         else
         xa = k1*strain + k2*strain^2 + k3*strain^3;
         end
         En2j1=( E(n,j+1)-((xa+P(n,j+1))/2.d0+ qbar)*(V(n+2,j+1)-V(n,j+1))+ deltaZ )/ (1.d0 + xb*(V(n+2,j+1)-V(n,j+1))/2.d0);
          diffE=    E(n+2,j+1)-En2j1;
         if (abs(diffE) > 0.d-6 )
          fprintf('Energy not converged %2.0f %2.0f',En2j1,E(n+2,j+1));
% c         if (
         end
        elseif (  ieos(2,j+1) == 0)
        %!this is for not used nodes
        else
% c        print *,'Non-MG EOS',j+1,iEOS(2,j+1)
% c        read(*,*) foo
         end
        end
    
      
%%      
% cccccccccccccccccccccccccc   Spall   ccccccccccccccccccccccccccccccccccccccccccc
% c     Check for tensile fracture, ie does the pressure/stress exceed pfrac at j?
% c     check for physical separation
% c j=5
       if (1 == 1)
       for j=0+1:2:jj-2+1
%            if ibc(j+2) == 2 || ibc(j+2)==1
%                disp('Single Node No Seperation Possible')
%            end
        if (-P(n+2,j+1)+s1(n+2,j+1) > pfrac(j)) && ibc(j+2)~=2 && ibc(j+2)~=1 && ibc(j-1)~=9
         fprintf('j=%2.4f n=%2.4f P=%2.4f Pfrac=%2.4f\n',j,n,-P(n+2,j+1)+s1(n+2,j+1),pfrac(j));
         fprintf('Tensile Fracture!!\n');
         fprintf('Separating Node %2.4f at n= %2.4f\n',j,n);
         fprintf('Before Separating:\n');
         fprintf('The node at:\n');
         fprintf('left: r(j-2 ,n+2)= %2.4f U(n-1)=%2.4f\n',r(n+2,j-2),U(n+1,j-2));
         fprintf('J:    r(j   ,n+2)= %2.4f U(n-1)=%2.4f\n',r(n+2,j)  ,U(n+1,j+0));
         fprintf('right:r(j+2 ,n+2)= %2.4f U(n-1)=%2.4f\n',r(n+2,j+2),U(n+1,j+2));
         fprintf('      r(j+4 ,n+2)= %2.4f U(n-1)=%2.4f\n',r(n+2,j+4),U(n+1,j+4));
         fprintf('Frac %2.4f\n',-P(n+2,j+1)+s1(n+2,j+1),pfrac(j));
         %-1     0     0     0     0     0     0     0     0     0     0     0     1     9     9     9     9     9     9     9     9
         %-1     0     0     0     0     0     0     0     0     0     0     0     0     0     1     9     9     9     9     9     9
         %-1     0     0     0     2     9    -2     0     0     0     0     0     0     0     1     9     9     9     9     9     9
         % 1     1     1     1     1
         for jjj = jj-4:-2:j+2
          for nz = n:n+2
               ibc(    jjj-0)  =    ibc(    jjj-2);
               ibc(    jjj-1)  =    ibc(    jjj-3);
               
              ieos(1,  jjj-0)  =   ieos(1,  jjj-2);
              ieos(1,  jjj-1)  =   ieos(1,  jjj-3);
              
              ieos(2,  jjj-0)  =   ieos(2,  jjj-2);
              ieos(2,  jjj-1)  =   ieos(2,  jjj-3);
              
                 m(    jjj-0)  =      m(    jjj-2);
                 m(    jjj-1)  =      m(    jjj-3);
                 
               rho(    jjj-0)  =    rho(    jjj-2);  
               rho(    jjj-1)  =    rho(    jjj-3);
 
                 r(nz ,jjj-0)  =     r(nz ,jjj-2);
                 r(nz ,jjj-1)  =     r(nz ,jjj-3);
                 
                 U(nz ,jjj-0)  =     U(nz ,jjj-2);
                 U(nz ,jjj-1)  =     U(nz ,jjj-3);
                 
               phi(nz ,jjj-0) =    phi(nz ,jjj-2);
               phi(nz ,jjj-1) =    phi(nz ,jjj-3);
               
              beta(nz ,jjj-0) =   beta(nz ,jjj-2);
              beta(nz ,jjj-1) =   beta(nz ,jjj-3);
              
            sigmar(nz ,jjj-0) = sigmar(nz ,jjj-2);
            sigmar(nz ,jjj-1) = sigmar(nz ,jjj-3);
            
            sigmao(nz ,jjj-0) = sigmao(nz ,jjj-2);
            sigmao(nz ,jjj-1) = sigmao(nz ,jjj-3);
            
                 V(nz ,jjj-0) =      V(nz ,jjj-2);            
                 V(nz ,jjj-1) =      V(nz ,jjj-3);
                 
                s1(nz ,jjj-0) =     s1(nz ,jjj-2);                 
                s1(nz ,jjj-1) =     s1(nz ,jjj-3);
                
                s2(nz ,jjj-0) =     s2(nz ,jjj-2);                
                s2(nz ,jjj-1) =     s2(nz ,jjj-3);
                
                s3(nz ,jjj-0) =     s3(nz ,jjj-2);                
                s3(nz ,jjj-1) =     s3(nz ,jjj-3);
                
                 q(nz ,jjj-0) =      q(nz ,jjj-2);                
                 q(nz ,jjj-1) =      q(nz ,jjj-3);
                 
             epsi1(nz ,jjj-0) =  epsi1(nz ,jjj-2);                 
             epsi1(nz ,jjj-1) =  epsi1(nz ,jjj-3);
             
             epsi2(nz ,jjj-0) =  epsi2(nz ,jjj-2);             
             epsi2(nz ,jjj-1) =  epsi2(nz ,jjj-3);
             
                 E(nz ,jjj-0) =      E(nz ,jjj-2);             
                 E(nz ,jjj-1) =      E(nz ,jjj-3);
                 
                 Y(    jjj-0) =      Y(    jjj-2);                 
                 Y(    jjj-1) =      Y(    jjj-3);
                 
                 K(nz ,jjj-0) =      K(nz ,jjj-2);                
                 K(nz ,jjj-1) =      K(nz ,jjj-3);
               
               
%               ieos(1,  jjj-1)  =   ieos(1,  jjj-3);
%               ieos(2,  jjj-1)  =   ieos(2,  jjj-3);
%               ieos(1,  jjj-0)  =   ieos(1,  jjj-2);
%               ieos(2,  jjj-0)  =   ieos(2,  jjj-2);
%                  m(    jjj-1)  =      m(    jjj-3);
%                rho(    jjj-1)  =    rho(    jjj-3);
%  
%                  r(nz ,jjj-0)  =     r(nz ,jjj-2);
%                  r(nz ,jjj-1)  =     r(nz ,jjj-3);
%                  U(nz ,jjj-0)  =     U(nz ,jjj-2);
%                phi(nz ,jjj-0) =    phi(nz ,jjj-2);
%               beta(nz ,jjj-0) =   beta(nz ,jjj-2);
%             sigmar(nz ,jjj-1) = sigmar(nz ,jjj-3);
%             sigmao(nz ,jjj-1) = sigmao(nz ,jjj-3);
%                  V(nz ,jjj-1) =      V(nz ,jjj-3);
%                 s1(nz ,jjj-1) =     s1(nz ,jjj-3);
%                 s2(nz ,jjj-1) =     s2(nz ,jjj-3);
%                 s3(nz ,jjj-1) =     s3(nz ,jjj-3);
%                  q(nz ,jjj-1) =      q(nz ,jjj-3);
%              epsi1(nz ,jjj-1) =  epsi1(nz ,jjj-3);
%              epsi2(nz ,jjj-1) =  epsi2(nz ,jjj-3);
%                  E(nz ,jjj-1) =      E(nz ,jjj-3);
%                  Y(    jjj-1) =      Y(    jjj-3);
%                  K(nz ,jjj-1) =      K(nz ,jjj-3);
               
               

             pfrac(    jjj-0) =  pfrac(    jjj-2);  %!node value
             pfrac(    jjj-1) =  pfrac(    jjj-3);  %!Cell value
          end
         end
% c
           for nz = n:n+2  %!void cell center
                  U(nz,j+1) = 0;
                  U(nz,j+1) = 0;
               phi(nz ,j+1) = 0;
              beta(nz ,j+1) = 0;
            sigmar(nz ,j+1) = 0;
            sigmao(nz ,j+1) = 0;
                 V(nz ,j+1) = 0;
                s1(nz ,j+1) = 0;
                s2(nz ,j+1) = 0;
                s3(nz ,j+1) = 0;
                 E(nz ,j+1) = 0;
                 Y(    j+1) = 0;
                 q(nz ,j+1) = 0;
             pfrac(    j+1) = 0;
             pfrac(    j+1) = 0;
           end
                 ibc(j)   = 2; %!outer
                 ibc(j+1) = 9; %!void
                 ibc(j+2) =-2; %!inner
                 ieos(1,j+1) =0
                 ieos(2,j+1) =0
               P(n+1,j+1) = pvoid;
               P(n+2,j+1) = pvoid;
               P(n+1,j-1) = pvoid;
               P(n+2,j-1) = pvoid;
               P(n+1,j+3) = pvoid;
               P(n+2,j+3) = pvoid;
               U(n+1,j+2) = U(n+1,j+4);
               r(n+1,j+1) = r(n+1,j);
               r(n+2,j+1) = r(n+2,j);
               %pfrac(j+2) =  pfrac(j+4);  %!node value
               %pfrac(j+3) =  pfrac(j+4);  %!node value
               pfrac(j) =  pfrac(j-1);  %!node value NSH CHANGE
               pfrac(j+2) =  pfrac(j+3);  %!node value NSH CHANGE
         fprintf('After Separating: %2.4f %2.4f\n',n,j);
         fprintf('r(j-2)=%2.4f  U= %2.4f  P=%2.4f\n',r(n+2,j-2),U(n+1,j-2),P(n+2,j-2));
         fprintf('r(j-1)=%2.4f  U= %2.4f  P=%2.4f\n',r(n+2,j-1),U(n+1,j-1),P(n+2,j-1));
         fprintf('r(j)  =%2.4f  U= %2.4f  P=%2.4f\n',r(n+2,j  ),U(n+1,j-0),P(n+2,j+0));
         fprintf('r(j+1)=%2.4f  U= %2.4f  P=%2.4f\n',r(n+2,j+1),U(n+1,j+1),P(n+2,j+1));
         fprintf('r(j+2)=%2.4f  U= %2.4f  P=%2.4f\n',r(n+2,j+2),U(n+1,j+2),P(n+2,j+2));
         fprintf('r(j+3)=%2.4f  U= %2.4f  P=%2.4f\n',r(n+2,j+3),U(n+1,j+3),P(n+2,j+3));
         fprintf('r(j+4)=%2.4f  U= %2.4f  P=%2.4f\n',r(n+2,j+4),U(n+1,j+4),P(n+2,j+4));
%                read(*,*) foo
% c
       for jzz=1:1:jj-2
% c       write(*,'(I4,5I2,10f7.3)')
%        write(*,'(I4,5I2,3f7.3,2e11.3,5f7.3)')
%      &  jzz,ieos(1,jzz),ieos(2,jzz),
%      &  ibc(jzz-1),ibc(jzz),ibc(jzz+1),
%      &  r(n+2,jzz-1),r(n+2,jzz),r(n+2,jzz+1),
%      &  U(n+1,jzz-1),U(n+1,jzz),U(n+1,jzz+1),
%      &  p(n+2,jzz)
% c     &  ,v(n+2,jzz+1),m(jzz+1),rho(jzz+1),y(jzz+1),pfrac(jzz+2)
        end
%       read(*,*) foo
        end
       end %! this is the fracture check
       end
       
       
       
% c
% c******************  Spall End  ***************************
% c       TIME-STEP
% c
% c      print *,'line 1258'
        dt_min =  1.d9;
        dr_min   =  r(n+2,0+2+1)-r(n+2,0+1); % !this is just to get us started
% c
        jzz    = -1;
       for j = 0+1:2:jj-2+1
       if (ibc(j+1) == 0)
% c       deltat =     t(n+2)    -t(n)
        Vdot   =   ( V(n+2,j+1)-V(n  ,j+1) )/deltat;
        deltar = abs(r(n+2,j+2)-r(n+2,j  ) );
% c
        b      = 8.d0*(Co^2+CL)*deltar*(Vdot/V(n+1,j+1));
        if (Vdot/V(n+1,j+1) >= 0.) 
            b = 0.d0;
        end
% c
        rho_local = rho(j+1)/V(n+2,j+1);
% c       print *,'rho_local',j,rho_local,rho(j+1),1.d0/v(n+2,j+1)
        if (rho_local > 1.e10)
        fprintf('Error: density too high\n');
        fprintf('%2.4f %2.4f %2.4f\n',j,ncount,t(n));
        fprintf('%2.4f %2.4f\n',rho(j+1),v(n+2,j+1));
        fprintf('%2.4f %2.4f %2.4f\n',j,r(n+2,j+2),r(n+2,j));
        fprintf('%2.4f %2.4f\n',V(n+2,j+1),V(n,j+1));
        fprintf('%2.4f %2.4f %2.4f %2.4f\n',P(n+1,j),p(n+1,j+1),P(n+2,j),p(n+2,j+1));
 %       read(*,*) foo
        end
% c       print *,'rho local',rho_local
        a  = sqrt(abs(P(n+2,j+1))/ rho_local ); %!I checked units
% c       print *,'ts',j,a,b,deltar,Vdot,V(n+1,j+1)
        if ( (a^2+b^2) ~= 0.d0 )
         delt_temp = (2.d0/3.d0)*(deltar/sqrt(a^2+b^2));
         if (delt_temp < dt_min)
          jzz    = j;
          dt_min = delt_temp;
          a_min = a;
          b_min = b;
          rho_min=rho_local;
% c        print *,'rr',rho_min,deltat,jzz,a_min,b_min,p(n+2,j+1)
% c        read(*,*) foo
         end
         end
         end %!ibc(j+1)
         end %!j loop


%         deltat = (dt_min)/20.d0
if deltat>(dt_min)/6
    deltat=(dt_min)/6
end

DeltatSave(timesteps)=dt_min;
 
%        dt_min=dt_min/5.d0;    %!Here is where you can arbitrarly lower time step
%                               %!the demomenator is suppose to be 1.d0
%        delt = t(n+1)-t(n-1);
%        if (dt_min > 1.1d0*delt )
%          dt_min = 1.1d0 * delt;
%        end
% ccccc      deltat = (dt_min+delt)/2.d0
% c
% c      t(n+3) = t(n+2) + deltat           ! this is now advanced along with all the
% c      t(n+4) = t(n+2) + deltat + deltat  ! other state variables
% c
       if (idebug == 1)
       if (icontact == 1)
       fprintf('End of first contact');
              if(idebug == 1 )
              end
              end
       end
% c        if (deltat .lt. 0.1d-5) then  ! this is just larry swabies suggestion for stability
% c        print *,rho_min,deltat,jzz,a_min,b_min
% c        endif


%%

% c*******************OUTPUT SOLUTION to Screen**********************
% c
       if (idebug == 1 && ncount > ndebug)
       for j = 1+1:2:jj-2+1
         %write(*,'(2I4,11f7.4)')
      %& n+1,j,t(n),r(n+1,j-1),r(n+1,j+1),U(n+1,j-1),U(n+1,j+1)
      %&,P(n+2,j),q(n+2,j),e(n+2,j),V(n+2,j)
       end
       %read(*,*) foo
       end
% c
       jzz = ncount;
       if (int32((jzz-1)/iskip) == single(jzz-1)/single(iskip) || t(n+2) >= tstop)
% c      if (n .eq. 1 .or. t(n) .ge. tskip ) then
       qtotal  = 0.d0;  %! total artifical viscosity
       mvtotal = 0.d0;  %! total momentum
       etotal  = 0.d0;  %! total energy
       ketotal = 0.d0;  %! total kinetic energy
       ietotal = 0.d0;  %! total internal energy
       ke3total= 0.d0;  %! total kinetic energy
% c
       for j = 0+1:2:jj-2+1
        if (ibc(j+1) == 0)
        qtotal =qtotal  + q(n+2,j);
        mvtotal=mvtotal + m(j+1)*(  U(n+1,j) + U(n+1,j+2) )/2.d0;
        ketotal=ketotal + m(j+1)*(((U(n+1,j) + U(n+1,j+2))/2.d0)^2)/2.d0;
        etotal =etotal  +(E(n+2,j+1)-E(n,j+1))/deltat+ P(n+1,j+1)*(V(n+2,j+1)-V(n,j+1))/deltat;
        if (ieos(1,j-1+1) == 2)
        ke3total=ke3total + m(j+1)*(((U(n+1,j) + U(n+1,j+2))/2.d0)^2)/2.d0;
        end
        end
       end
       fprintf('%1.0f %2.4f %2.4f %2.6f %2.4f %2.4f \n',timesteps,t(n+1),mvtotal,ketotal,ietotal,etotal);
% 
%         write(*,'(1I6,5e15.6)')
%      & n,t(n+1),mvtotal,KEtotal,ietotal,etotal
% 
% c       ioffset = 10
% c       do j=98-ioffset,98,2
% c         write(*,99)
% c     & n+2,j+1,t(n),r(n+1,j),r(n+1,j+2),U(n+1,j),U(n+1,j+2)
% c     &  ,P(n+2,j+1),e(n+2,j+1),q(n+1,j+1)
% c       enddo
% c       do j=98,98+ioffset,2
% c         write(*,99)
% c     & n+2,j+1,t(n),r(n+1,j),r(n+1,j+2),U(n+1,j),U(n+1,j+2)
% c     &  ,P(n+2,j+1),e(n+2,j+1),q(n+1,j+1)
% c       enddo
% c      print *, ' '
% c
% c      skip = (jj/2)/256
% c      skip = nn/skip
% c      do j=ioffset,nn,skip
% c      write(*,98)
% c     &'n','Cell','dt','rad','rad','Vel','Vel'
% c     &,'pres','Vol'
% c      print *,' '
% c      enddo
% c      write(*,'(2I4,6e9.2,6f7.4)')
% c     & n,j,deltat,r(n,j-1),r(n,j+1),U(n,j-1),U(n,j+1),P(n,j)
% c     & ,V(n,j),rho(j)
       end


       
       %%
       
       
% c
% c ******************* Write solution to File  *************************
% c       write(*,*) 'Enter file name'
% c       read(*,*) foo name
% c       name='al.txt'
% c       open (unit=33, file=name, form='formatted', status='unknown')
% c       do n=1,nn,2
         if (t(n) == 0. || t(n) >= tskip )
           for j=0+1:2:jj-4
          bs1 =     V(n  ,j+1); %!rho(j+1)
          bs2 =     U(n+1,j  );
          bs3 =     P(n  ,j+1);
          bs4 =     E(n  ,j+1);
          bs5 =     q(n+1,j+1);
          bs6 =    s1(n  ,j+1);
          bs7 =    s2(n  ,j+1);
          bs8 = epsi1(n-1,j+1);
          bs10=  Temp(n  ,j+1);
 
          if (abs(bs1)   < 1.d-99)
              bs1 = 0.d0;
          end%!this is because if the expoent is
          if (abs(U(n+1,j  ))    < 1.d-99)
              bs2 = 0.d0;
          end%! larger than 99 it will not write
          if (abs(bs3)    < 1.d-99)
              bs3 = 0.d0;
          end%! to the file correctly, so I just
          if (abs(bs4)    < 1.d-99)
              bs4 = 0.d0;
          end%! zero it out
          if (abs(bs5)    < 1.d-99)
              bs5 = 0.d0;
          end
          if (abs(bs6)    < 1.d-99)
              bs6 = 0.d0;
          end
          if (abs(bs7)    < 1.d-99)
              bs7 = 0.d0;
          end
          if (abs(bs8)    < 1.d-99)
              bs8 = 0.d0;
          end
          bs9 = -bs3 + bs6;
 
% c      epsi1(n+1,j+1) = (U(n+1,j+2)-U(n+1,j))/(r(n+1,j+2)-r(n+1,j))
% c      epsi2(n+1,j+1)
% c       if (j .eq. 3140) print *,'u',U(n+1,j)
%             write(33,'(2I8,3I3,18e26.18)')
%      &     ncount, j+1   , ibc(j)  ,ibc(j+1),ibc(j+2),
%      &     t(n)  , r(n,j), r(n,j+1),r(n,j+2),     bs1,
%      &     bs2   ,    bs3, rho(j+1),     bs4,     bs5,
%      &     bs6   ,    bs7, k(n,j+1), ketotal,     bs9,
%      &     bs10  ,   sigmar(n,j+1), p(n,j+1)-s1(n,j+1)
% c             write(33,'(4d26.18)')
% c     &   t(n),bs3,r(n,j),bs2
% c
           end  %! this is j loop
            % write(35,'(d26.18)') ke3total
       tskip = tskip + dtskip;
         end

       
       
       
       %%
       
% c
% cccccccccccccccccccccc Update solution ccccccccccccccccccccccccccccc
% c
       for j=0+1:jj+1
       if (ibc(j) ~= 9)
% c
         U(n-1,j)       = U(n+1,j);
         U(n  ,j)       = U(n+2,j);
         phi(n-1,j)     = phi(n+1,j);
         phi(n  ,j)     = phi(n+2,j);
         sigmar(n-1,j)  = sigmar(n+1,j);
         sigmar(n  ,j)  = sigmar(n+2,j);
         sigmao(n-1,j)  = sigmao(n+1,j);
         sigmao(n  ,j)  = sigmao(n+2,j);
         beta(n-1,j)    = beta(n+1,j);
         beta(n  ,j)    = beta(n+2,j);
         V(n-1,j)       = V(n+1,j);
         V(n  ,j)       = V(n+2,j);
         r(n-1,j)       = r(n+1,j);
         r(n  ,j)       = r(n+2,j);
         epsi1(n-1,j)   = epsi1(n+1,j);
         epsi1(n  ,j)   = epsi1(n+2,j);
         epsi2(n-1,j)   = epsi2(n+1,j);
         epsi2(n  ,j)   = epsi2(n+2,j);
         s1(n-1,j)      = s1(n+1,j);
         s1(n  ,j)      = s1(n+2,j);
         s2(n-1,j)      = s2(n+1,j);
         s2(n  ,j)      = s2(n+2,j);
         s3(n-1,j)      = s3(n+1,j);
         s3(n  ,j)      = s3(n+2,j);
         P(n-1,j)       = P(n+1,j);
         P(n  ,j)       = P(n+2,j);
         q(n-1,j)       = q(n+1,j);
         q(n  ,j)       = q(n+2,j);
         t(n-1)         = t(n+1);
         t(n  )         = t(n+2);
         E(n-1,j)       = E(n+1,j);
         E(n  ,j)       = E(n+2,j);
         K(n-1,j)       = K(n+1,j);
         K(n  ,j)       = K(n+2,j);
         Temp(n-1,j)    = Temp(n+1,j);
         Temp(n  ,j)    = Temp(n+2,j);
         entropy(n-1,j) = entropy(n+1,j);
         entropy(n  ,j) = entropy(n+2,j);
% c        icompact(n-1,j)= icompact(n+1,j)
% c        icompact(n  ,j)= icompact(n+2,j)
% 
       end
       end
% U = round(U,6,'significant');
% phi = round(phi,6,'significant');
% sigmar = round(sigmar,6,'significant');
% sigmao = round(sigmao,6,'significant');
% % beta = round(beta,6,'significant');
% V = round(V,8,'significant');
% % r = round(r,6,'significant');
% % epsi1 = round(epsi1,6,'significant');
% % epsi2 = round(epsi2,6,'significant');
% % s1 = round(s1,6,'significant');
% % s2 = round(s2,6,'significant');
% % s3 = round(s3,6,'significant');
% P = round(P,6,'significant');
% q = round(q,6,'significant');
% E = round(E,6,'significant');
% K = round(K,6,'significant');
% Temp = round(Temp,6,'significant');
% entropy = round(entropy,6,'significant');
       
       
% cccccccccccccccccccccc  Main loop closed ccccccccccccccccccccccccccccccccccccc
%        if (t(n) .le. tstop) goto 1
% c      enddo  ! n or time loop
%        close(33)
%       print *,'*******************  Finished!!  ******************'
% c      read(*,*) foo
    clf

    subplot(4,2,1)
    plot(r(1,1:2:jj+1),U(1,1:2:jj+1)*10000,'k.');
    ylabel('velocity (m/s)')
    
    
    subplot(4,2,3)
    plot(r(1,1:2:jj),P(1,2:2:jj),'k.');    
    ylabel('pressure')
    
    
    subplot(4,2,5)
    plot(r(1,1:2:jj),rho(1,2:2:jj)./V(1,2:2:jj),'k.'); 
    ylabel('density')
    
    
    subplot(4,2,[2,4,6])
    plot(r(1,1:2:jj),s1(1,2:2:jj)); 
    hold on
    plot(r(1,1:2:jj),s2(1,2:2:jj),'linewidth',3); 
    %plot(r(1,1:2:jj),s2(1,2:2:jj),'linewidth',3); 
    plot(r(1,1:2:jj),s3(1,2:2:jj)); 
    ylabel('s')

    subplot(4,2,[7,8])
    for i=1:2:jj
        if i==1
            x=[r(4,1);r(4,3);r(4,3);r(4,1)];
            y=[0;0;1/V(4,2);1/V(4,2)];
            c=[U(1,1).*10000];
        else
            x=[[x],[r(4,i);r(4,2+i);r(4,2+i);r(4,i)]];
            y=[[y],[0;0;1/V(4,1+i);1/V(4,1+i)]]   ;
            c=[[c],[U(1,i).*10000]]  ;
        end
    end

    
colormap(jet)
patch(x,y,c)
colorbar
title(sprintf('t=%4.4f',t(1)))
%('Ticks',[-1,1],'Ticklabels',{'-1','1'})
caxis([-100 100])
      
  pause(0.00001)    
       
 Tracery(timesteps)=U(1,49);
 Tracerx(timesteps)=t(1);
end
