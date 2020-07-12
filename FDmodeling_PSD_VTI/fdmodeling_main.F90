program fdmodeling
!this is the main program of finite-difference elastic wave equation modeling using standard staggered grid
!10th order in space, 2nd order in time
!this program will read 3 input files (VP, VS, DEN), and then calculate the elastic tensors
!one can define an anisotropic zone and input Thomsen anisotropic parameters, and the aisotropic zone is either VTI or HTI.
  use paramod
  implicit none
  integer(kind=INP) :: sctf,sctp,fspx,fspy,fspz,fgvx,fgvy,fgvz,fdcof
  integer(kind=INP) :: NX,NY,NZ,NT,boudp,SNX,SNY,SNZ,snapx,snapy,snapz,snapdt,snhmdt,NREC,ANIX(2),ANIY(2),ANIZ(2),ANITP
  real(kind=RLP)    :: dx,dy,dz,dt,freq,M0,Mxx,Myy,Mzz,Mxy,Mxz,Myz,vpmax,vsmin,vapsln,gamma,delta

  integer(kind=INP),allocatable,dimension(:)  :: RECX,RECY,RECZ
  real(kind=RLP),allocatable,dimension(:,:,:) :: VP,VS,DEN
  real(kind=RLP),allocatable,dimension(:,:,:) :: TXX,TYY,TZZ,TXY,TXZ,TYZ,VX,VY,VZ,ABR
  real(kind=RLP),allocatable :: Soutf(:)

  !local parameters
  integer(kind=INP) :: dhmax,dhmin,ii,it,ix,iy,iz
  real(kind=RLP)    :: fcmax,nsmp,nlim,dtlim,hfac,tstart,tend,b1,b2,b3,b4,b5,dtdx,dtdy,dtdz,ddv,ssignal,tempvl
  real(kind=RLP)    :: C11,C12,C13,C22,C23,C33,C44,C55,C66
  real(kind=RLP)    :: C44YP1,C44ZP1,C44ZYP1,C55XP1,C55ZP1,C55ZXP1,C66XP1,C66YP1,C66YXP1
  real(kind=RLP)    :: TXX_x,TXY_y,TXZ_z,TXY_x,TYY_y,TYZ_z,TXZ_x,TYZ_y,TZZ_z,denx,deny,denz
  real(kind=RLP)    :: VX_x,VX_y,VX_z,VY_x,VY_y,VY_z,VZ_x,VZ_y,VZ_z,c44yz,c55xz,c66xy
  character(20)     :: namchar
  real(kind=RLP),allocatable    :: tempo(:,:)

  !read input file for modeling parameters
  open(unit=10,file='input.dat')
  read(10,*)                    !Head line
  read(10,*) NX, NY, NZ
  read(10,*) ANIX
  read(10,*) ANIY
  read(10,*) ANIZ
  read(10,*) ANITP
  read(10,*) vapsln, gamma, delta
  read(10,*) dx, dy, dz
  read(10,*) NT
  read(10,*) dt
  read(10,*) boudp
  read(10,*) fdcof
  read(10,*)                    !Head line
  read(10,*) SNX, SNY, SNZ
  read(10,*) sctf
  read(10,*) sctp
  read(10,*) freq
  read(10,*) M0
  read(10,*) Mxx
  read(10,*) Myy
  read(10,*) Mzz
  read(10,*) Mxy
  read(10,*) Mxz
  read(10,*) Myz
  read(10,*)                    !Head line
  read(10,*) fspx, fspy, fspz
  read(10,*) snapx, snapy, snapz
  read(10,*) snapdt
  read(10,*) fgvx, fgvy, fgvz
  read(10,*) snhmdt
  read(10,*)                    !Head line
  read(10,*) vpmax, vsmin
  close(unit=10)

  !output some important modeling parameters
  write(*,*) "Modeling size: NX=", NX, ", NY=", NY, ", NZ=", NZ
  write(*,*) "Anisotropic zone X:", ANIX(1), "-", ANIX(2)
  write(*,*) "Anisotropic zone Y:", ANIY(1), "-", ANIY(2)
  write(*,*) "Anisotropic zone Z:", ANIZ(1), "-", ANIZ(2)
  if (ANITP==1) then
    write(*,*) "In this anisotropic zone, media have VTI anisotropy."
  else
    write(*,*) "In this anisotropic zone, media have HTI anisotropy."
  end if
  write(*,*) "Thomsen anisotropy parameters: ", vapsln, gamma, delta 
  write(*,*) "grid interval: dx=", dx, "m, dy=", dy, "m, dz=", dz, "m"
  write(*,*) "time steps:", NT
  write(*,*) "time sampling inteval:", dt
  write(*,*) "boundary points:", boudp
  if (fdcof==1) then
    write(*,*) "Adopt Taylor FD coefficients"
  else
    write(*,*) "Adopt Holberg FD coefficients"
  end if
  write(*,*) "source location: SNX=", SNX, ", SNY=", SNY, ", SNZ=", SNZ
  if(sctf==1) then
    write(*,*) "Using Ricker wavelet and main frequency of wavelet is", freq, "Hz."
  else
    stop "Error: incorrect input for source time function!"
  end if
  if(sctp==1) then
    write(*,*) "Using compressional source"
  else if(sctp==2) then
    write(*,*) "Using momtent tensor source"
    write(*,*) "Mxx=", Mxx
    write(*,*) "Myy=", Myy
    write(*,*) "Mzz=", Mzz
    write(*,*) "Mxy=", Mxy
    write(*,*) "Mxz=", Mxz
    write(*,*) "Myz=", Myz
  else
    stop "Error: incorrect input for source type!"
  end if

  !decide the limit for Grid dispersion and Courant Stability in specific spatial derivative order of FD. Refer to Sofi3D manual Table 2 & 3.
  if (fdcof==1) then
    !limit for Taylor coefficients
    nlim=5.0_8
    hfac=53089_8/40320_8
  else
    !limit for Holberg coefficients
    nlim=3.19_8
    hfac=1.38766_8
  end if

  !check Grid Dispersion (numerical artefacts)
  !for Taylor finite difference coefficient and 10th order in space: nsmp should larger than 5. Refer to Sofi3D manual Table 2.
  dhmax=max(dx,dy,dz)
  fcmax=2.0_8*freq
  nsmp=vsmin/(dhmax*fcmax)
  if (nsmp<nlim) then
    write(*,*) "Modeling has grid dispersion! For minimal S-wave wavelength, spatial sampling points per wavelength is", nsmp, " It should >", nlim
    stop 
  else 
    write(*,*) "Spatial sampling points for minimal S-wave wavelength is", nsmp
  end if

  !check Courant Stability
  !for Taylor finite difference coefficient and 10th order in space: hfac is 53089/40320. Refer to Sofi3D manual Table 3.
  dhmin=min(dx,dy,dz)
  dtlim=dhmin/(sqrt(3.0_8)*hfac*vpmax)
  if (dt>dtlim) then
    write(*,*) "Modeling unstable! dt should less than", dtlim, "s under current parameters."
    stop
  end if

  !allocate memory for arrays
  allocate(VP(NZ,NY,NX),VS(NZ,NY,NX),DEN(NZ,NY,NX))
  allocate(TXX(NZ,NY,NX),TYY(NZ,NY,NX),TZZ(NZ,NY,NX),TXY(NZ,NY,NX),TXZ(NZ,NY,NX),TYZ(NZ,NY,NX),VX(NZ,NY,NX),VY(NZ,NY,NX),VZ(NZ,NY,NX),ABR(NZ,NY,NX))
  allocate(Soutf(NT))
  allocate(tempo(NZ,NY))

  !Initialize some arrays
  TXX=0D0
  TYY=0D0
  TZZ=0D0
  TXY=0D0
  TXZ=0D0
  TYZ=0D0
  VX =0D0
  VY =0D0
  VZ =0D0
  ABR=1D0

  !normalize the input moment tensor, let sum(m_ij^2)=1
  !and also apply seismic moment (M0) on the moment tensor
  tempvl=sqrt(Mxx*Mxx+Myy*Myy+Mzz*Mzz+2.0_8*(Mxy*Mxy+Mxz*Mxz+Myz*Myz))
  Mxx=M0*Mxx/tempvl
  Myy=M0*Myy/tempvl
  Mzz=M0*Mzz/tempvl
  Mxy=M0*Mxy/tempvl
  Mxz=M0*Mxz/tempvl
  Myz=M0*Myz/tempvl

  !generate source time function
  if (sctf==1) then
    call Ricker()
  end if
  !scaled source time function with respect to spatical and temporal discretization
  !note the temporal term 'dt' disappers as we need to calculate the time derivative of source time function later
  ddv=2.0_8*dx*dy*dz
  Soutf=Soutf/ddv

  !input P-wave and S-wave velocity and density of the model
  open(unit=11,file='VP',form='unformatted',status='old',access='stream',action='read')
  read(11) VP
  close(unit=11)
  open(unit=12,file='VS',form='unformatted',status='old',access='stream',action='read')
  read(12) VS
  close(unit=12)
  open(unit=13,file='DEN',form='unformatted',status='old',access='stream',action='read')
  read(13) DEN
  close(unit=13)

  write(*,*) minval(VP), maxval(VP)
  write(*,*) minval(VS), maxval(VS)
  write(*,*) minval(DEN), maxval(DEN)

  !input receiver parameters
  open(unit=90,file='receiver.dat',form='formatted',status='old',action='read')
  read(90,*) NREC
  write(*,*) "Number of receivers:", NREC
  allocate(RECX(NREC),RECY(NREC),RECZ(NREC))
  do ii=1,NREC
    read(90,*) RECX(ii),RECY(ii),RECZ(ii)
  end do
  close(unit=90)

  !calculate parameters for boundary conditions, here we use absorbing boundary condition
  call absorbcdt()
  open(unit=99,file='abr.rec',form='unformatted',status='replace',access='stream',action='write')
  tempo=ABR(:,:,NX/2)
  write(99) tempo
  close(unit=99)

  !determine the finite-difference coefficient, here we use 10th order in space
  if(fdcof==1) then
    !Taylor coefficients
    b1=19845.0_8/16384.0_8
    b2=-735.0_8/8192.0_8
    b3=567.0_8/40960.0_8
    b4=-405.0_8/229376.0_8
    b5=35.0_8/294912.0_8
  else
    !Holberg coefficients
    b1=1.2415_8
    b2=-0.11231_8
    b3=0.026191_8
    b4=-0.0064682_8
    b5=0.001191_8
  end if

  !precalculate some parameters which will be used in modeling
  dtdx=dt/dx
  dtdy=dt/dy
  dtdz=dt/dz

  !open file for storing received trace data
  open(unit=101,file='vx.rec',form='unformatted',status='replace',access='stream',action='write')
  open(unit=102,file='vy.rec',form='unformatted',status='replace',access='stream',action='write')
  open(unit=103,file='vz.rec',form='unformatted',status='replace',access='stream',action='write')

  !open file for storing wavefield snapshot profile
  if (fspx==1) then
    open(unit=111,file='vxpx.snap',form='unformatted',status='replace',access='stream',action='write')
    open(unit=112,file='vypx.snap',form='unformatted',status='replace',access='stream',action='write')
    open(unit=113,file='vzpx.snap',form='unformatted',status='replace',access='stream',action='write')
  end if
  if (fspy==1) then
    open(unit=114,file='vxpy.snap',form='unformatted',status='replace',access='stream',action='write')
    open(unit=115,file='vypy.snap',form='unformatted',status='replace',access='stream',action='write')
    open(unit=116,file='vzpy.snap',form='unformatted',status='replace',access='stream',action='write')
  end if
  if (fspz==1) then
    open(unit=117,file='vxpz.snap',form='unformatted',status='replace',access='stream',action='write')
    open(unit=118,file='vypz.snap',form='unformatted',status='replace',access='stream',action='write')
    open(unit=119,file='vzpz.snap',form='unformatted',status='replace',access='stream',action='write')
  end if

  !start finite-difference modeling and timing
  call cpu_time(tstart)
  do it=1,NT
  !loop for time steps

    !$OMP PARALLEL DEFAULT(SHARED)
    !$OMP DO PRIVATE(TXX_x,TXY_y,TXZ_z,TXY_x,TYY_y,TYZ_z,TXZ_x,TYZ_y,TZZ_z,denx,deny,denz,ix,iy,iz) SCHEDULE(DYNAMIC)
    !update particle velocity
    do ix=6,NX-5
      do iy=6,NY-5
        do iz=6,NZ-5
          !calculate spatial derivatives of stress
          TXX_x=b1*(TXX(iz,iy,ix+1)-TXX(iz,iy,ix))+b2*(TXX(iz,iy,ix+2)-TXX(iz,iy,ix-1))+b3*(TXX(iz,iy,ix+3)-TXX(iz,iy,ix-2))+b4*(TXX(iz,iy,ix+4)-TXX(iz,iy,ix-3))+b5*(TXX(iz,iy,ix+5)-TXX(iz,iy,ix-4))
          TXY_y=b1*(TXY(iz,iy,ix)-TXY(iz,iy-1,ix))+b2*(TXY(iz,iy+1,ix)-TXY(iz,iy-2,ix))+b3*(TXY(iz,iy+2,ix)-TXY(iz,iy-3,ix))+b4*(TXY(iz,iy+3,ix)-TXY(iz,iy-4,ix))+b5*(TXY(iz,iy+4,ix)-TXY(iz,iy-5,ix))
          TXZ_z=b1*(TXZ(iz,iy,ix)-TXZ(iz-1,iy,ix))+b2*(TXZ(iz+1,iy,ix)-TXZ(iz-2,iy,ix))+b3*(TXZ(iz+2,iy,ix)-TXZ(iz-3,iy,ix))+b4*(TXZ(iz+3,iy,ix)-TXZ(iz-4,iy,ix))+b5*(TXZ(iz+4,iy,ix)-TXZ(iz-5,iy,ix))

          TXY_x=b1*(TXY(iz,iy,ix)-TXY(iz,iy,ix-1))+b2*(TXY(iz,iy,ix+1)-TXY(iz,iy,ix-2))+b3*(TXY(iz,iy,ix+2)-TXY(iz,iy,ix-3))+b4*(TXY(iz,iy,ix+3)-TXY(iz,iy,ix-4))+b5*(TXY(iz,iy,ix+4)-TXY(iz,iy,ix-5))
          TYY_y=b1*(TYY(iz,iy+1,ix)-TYY(iz,iy,ix))+b2*(TYY(iz,iy+2,ix)-TYY(iz,iy-1,ix))+b3*(TYY(iz,iy+3,ix)-TYY(iz,iy-2,ix))+b4*(TYY(iz,iy+4,ix)-TYY(iz,iy-3,ix))+b5*(TYY(iz,iy+5,ix)-TYY(iz,iy-4,ix))
          TYZ_z=b1*(TYZ(iz,iy,ix)-TYZ(iz-1,iy,ix))+b2*(TYZ(iz+1,iy,ix)-TYZ(iz-2,iy,ix))+b3*(TYZ(iz+2,iy,ix)-TYZ(iz-3,iy,ix))+b4*(TYZ(iz+3,iy,ix)-TYZ(iz-4,iy,ix))+b5*(TYZ(iz+4,iy,ix)-TYZ(iz-5,iy,ix))

          TXZ_x=b1*(TXZ(iz,iy,ix)-TXZ(iz,iy,ix-1))+b2*(TXZ(iz,iy,ix+1)-TXZ(iz,iy,ix-2))+b3*(TXZ(iz,iy,ix+2)-TXZ(iz,iy,ix-3))+b4*(TXZ(iz,iy,ix+3)-TXZ(iz,iy,ix-4))+b5*(TXZ(iz,iy,ix+4)-TXZ(iz,iy,ix-5))
          TYZ_y=b1*(TYZ(iz,iy,ix)-TYZ(iz,iy-1,ix))+b2*(TYZ(iz,iy+1,ix)-TYZ(iz,iy-2,ix))+b3*(TYZ(iz,iy+2,ix)-TYZ(iz,iy-3,ix))+b4*(TYZ(iz,iy+3,ix)-TYZ(iz,iy-4,ix))+b5*(TYZ(iz,iy+4,ix)-TYZ(iz,iy-5,ix))
          TZZ_z=b1*(TZZ(iz+1,iy,ix)-TZZ(iz,iy,ix))+b2*(TZZ(iz+2,iy,ix)-TZZ(iz-1,iy,ix))+b3*(TZZ(iz+3,iy,ix)-TZZ(iz-2,iy,ix))+b4*(TZZ(iz+4,iy,ix)-TZZ(iz-3,iy,ix))+b5*(TZZ(iz+5,iy,ix)-TZZ(iz-4,iy,ix))

          !interpolate density
          denx=(DEN(iz,iy,ix)+DEN(iz,iy,ix+1))/2.0_8
          deny=(DEN(iz,iy,ix)+DEN(iz,iy+1,ix))/2.0_8
          denz=(DEN(iz,iy,ix)+DEN(iz+1,iy,ix))/2.0_8

          !calculate particle velocity of the next time step
          VX(iz,iy,ix)=(VX(iz,iy,ix)+(TXX_x*dtdx+TXY_y*dtdy+TXZ_z*dtdz)/denx)*ABR(iz,iy,ix)
          VY(iz,iy,ix)=(VY(iz,iy,ix)+(TXY_x*dtdx+TYY_y*dtdy+TYZ_z*dtdz)/deny)*ABR(iz,iy,ix)
          VZ(iz,iy,ix)=(VZ(iz,iy,ix)+(TXZ_x*dtdx+TYZ_y*dtdy+TZZ_z*dtdz)/denz)*ABR(iz,iy,ix)
        end do
      end do
    end do
    !$OMP END DO

    !$OMP DO PRIVATE(VX_x,VX_y,VX_z,VY_x,VY_y,VY_z,VZ_x,VZ_y,VZ_z,c44yz,c55xz,c66xy,ix,iy,iz) &
    !$OMP &  PRIVATE(C11,C12,C13,C22,C23,C33,C44,C55,C66,C44YP1,C44ZP1,C44ZYP1,C55XP1,C55ZP1,C55ZXP1,C66XP1,C66YP1,C66YXP1) &
    !$OMP &  SCHEDULE(DYNAMIC)
    !update stress
    do ix=6,NX-5
      do iy=6,NY-5
        do iz=6,NZ-5
          !calculate spatial derivatives of particle velocity
          VX_x=b1*(VX(iz,iy,ix)-VX(iz,iy,ix-1))+b2*(VX(iz,iy,ix+1)-VX(iz,iy,ix-2))+b3*(VX(iz,iy,ix+2)-VX(iz,iy,ix-3))+b4*(VX(iz,iy,ix+3)-VX(iz,iy,ix-4))+b5*(VX(iz,iy,ix+4)-VX(iz,iy,ix-5))
          VX_y=b1*(VX(iz,iy+1,ix)-VX(iz,iy,ix))+b2*(VX(iz,iy+2,ix)-VX(iz,iy-1,ix))+b3*(VX(iz,iy+3,ix)-VX(iz,iy-2,ix))+b4*(VX(iz,iy+4,ix)-VX(iz,iy-3,ix))+b5*(VX(iz,iy+5,ix)-VX(iz,iy-4,ix))
          VX_z=b1*(VX(iz+1,iy,ix)-VX(iz,iy,ix))+b2*(VX(iz+2,iy,ix)-VX(iz-1,iy,ix))+b3*(VX(iz+3,iy,ix)-VX(iz-2,iy,ix))+b4*(VX(iz+4,iy,ix)-VX(iz-3,iy,ix))+b5*(VX(iz+5,iy,ix)-VX(iz-4,iy,ix))
          VY_x=b1*(VY(iz,iy,ix+1)-VY(iz,iy,ix))+b2*(VY(iz,iy,ix+2)-VY(iz,iy,ix-1))+b3*(VY(iz,iy,ix+3)-VY(iz,iy,ix-2))+b4*(VY(iz,iy,ix+4)-VY(iz,iy,ix-3))+b5*(VY(iz,iy,ix+5)-VY(iz,iy,ix-4))
          VY_y=b1*(VY(iz,iy,ix)-VY(iz,iy-1,ix))+b2*(VY(iz,iy+1,ix)-VY(iz,iy-2,ix))+b3*(VY(iz,iy+2,ix)-VY(iz,iy-3,ix))+b4*(VY(iz,iy+3,ix)-VY(iz,iy-4,ix))+b5*(VY(iz,iy+4,ix)-VY(iz,iy-5,ix))
          VY_z=b1*(VY(iz+1,iy,ix)-VY(iz,iy,ix))+b2*(VY(iz+2,iy,ix)-VY(iz-1,iy,ix))+b3*(VY(iz+3,iy,ix)-VY(iz-2,iy,ix))+b4*(VY(iz+4,iy,ix)-VY(iz-3,iy,ix))+b5*(VY(iz+5,iy,ix)-VY(iz-4,iy,ix))
          VZ_x=b1*(VZ(iz,iy,ix+1)-VZ(iz,iy,ix))+b2*(VZ(iz,iy,ix+2)-VZ(iz,iy,ix-1))+b3*(VZ(iz,iy,ix+3)-VZ(iz,iy,ix-2))+b4*(VZ(iz,iy,ix+4)-VZ(iz,iy,ix-3))+b5*(VZ(iz,iy,ix+5)-VZ(iz,iy,ix-4))
          VZ_y=b1*(VZ(iz,iy+1,ix)-VZ(iz,iy,ix))+b2*(VZ(iz,iy+2,ix)-VZ(iz,iy-1,ix))+b3*(VZ(iz,iy+3,ix)-VZ(iz,iy-2,ix))+b4*(VZ(iz,iy+4,ix)-VZ(iz,iy-3,ix))+b5*(VZ(iz,iy+5,ix)-VZ(iz,iy-4,ix))
          VZ_z=b1*(VZ(iz,iy,ix)-VZ(iz-1,iy,ix))+b2*(VZ(iz+1,iy,ix)-VZ(iz-2,iy,ix))+b3*(VZ(iz+2,iy,ix)-VZ(iz-3,iy,ix))+b4*(VZ(iz+3,iy,ix)-VZ(iz-4,iy,ix))+b5*(VZ(iz+4,iy,ix)-VZ(iz-5,iy,ix))

          !calculate the elastic tensors
          if ((iz<ANIZ(1)).or.(iz>ANIZ(2)).or.(iy<ANIY(1)).or.(iy>ANIY(2)).or.(ix<ANIX(1)).or.(ix>ANIX(2)))  then
              !isotropic media
              C11=VP(iz,iy,ix)*VP(iz,iy,ix)*DEN(iz,iy,ix)
              C22=C11
              C33=C11
              C44=VS(iz,iy,ix)*VS(iz,iy,ix)*DEN(iz,iy,ix)
              C55=C44
              C66=C44
              C12=C11-2.0_8*C44
              C13=C12
              C23=C12
              C44YP1=VS(iz,iy+1,ix)*VS(iz,iy+1,ix)*DEN(iz,iy+1,ix)
              C44ZP1=VS(iz+1,iy,ix)*VS(iz+1,iy,ix)*DEN(iz+1,iy,ix)
              C44ZYP1=VS(iz+1,iy+1,ix)*VS(iz+1,iy+1,ix)*DEN(iz+1,iy+1,ix)
              C55XP1=VS(iz,iy,ix+1)*VS(iz,iy,ix+1)*DEN(iz,iy,ix+1)
              C55ZP1=C44ZP1
              C55ZXP1=VS(iz+1,iy,ix+1)*VS(iz+1,iy,ix+1)*DEN(iz+1,iy,ix+1)
              C66XP1=C55XP1
              C66YP1=C44YP1
              C66YXP1=VS(iz,iy+1,ix+1)*VS(iz,iy+1,ix+1)*DEN(iz,iy+1,ix+1)
          else if (ANITP==1) then
              !VTI
              C33=VP(iz,iy,ix)*VP(iz,iy,ix)*DEN(iz,iy,ix)
              C44=VS(iz,iy,ix)*VS(iz,iy,ix)*DEN(iz,iy,ix)
              C11=2.0_8*vapsln*C33+C33
              C66=2.0_8*gamma*C44+C44
              C13=sqrt(2.0_8*delta*C33*(C33-C44)+(C33-C44)*(C33-C44))-C44
              C22=C11
              C55=C44
              C12=C11-2.0_8*C66
              C23=C13
              C44YP1=VS(iz,iy+1,ix)*VS(iz,iy+1,ix)*DEN(iz,iy+1,ix)
              C44ZP1=VS(iz+1,iy,ix)*VS(iz+1,iy,ix)*DEN(iz+1,iy,ix)
              C44ZYP1=VS(iz+1,iy+1,ix)*VS(iz+1,iy+1,ix)*DEN(iz+1,iy+1,ix)
              C55XP1=VS(iz,iy,ix+1)*VS(iz,iy,ix+1)*DEN(iz,iy,ix+1)
              C55ZP1=C44ZP1
              C55ZXP1=VS(iz+1,iy,ix+1)*VS(iz+1,iy,ix+1)*DEN(iz+1,iy,ix+1)
              C66XP1=2.0_8*gamma*(VS(iz,iy,ix+1)*VS(iz,iy,ix+1)*DEN(iz,iy,ix+1))+(VS(iz,iy,ix+1)*VS(iz,iy,ix+1)*DEN(iz,iy,ix+1))
              C66YP1=2.0_8*gamma*(VS(iz,iy+1,ix)*VS(iz,iy+1,ix)*DEN(iz,iy+1,ix))+(VS(iz,iy+1,ix)*VS(iz,iy+1,ix)*DEN(iz,iy+1,ix))
              C66YXP1=2.0_8*gamma*(VS(iz,iy+1,ix+1)*VS(iz,iy+1,ix+1)*DEN(iz,iy+1,ix+1))+(VS(iz,iy+1,ix+1)*VS(iz,iy+1,ix+1)*DEN(iz,iy+1,ix+1))
          else
              !HTI
              call htic6r(VP(iz,iy,ix),VS(iz,iy,ix),DEN(iz,iy,ix),vapsln,gamma,delta,C11,C12,C13,C22,C23,C33,C44,C55,C66)
              call htic6rc44(VP(iz,iy+1,ix),VS(iz,iy+1,ix),DEN(iz,iy+1,ix),vapsln,gamma,delta,C44YP1)
              call htic6rc44(VP(iz+1,iy,ix),VS(iz+1,iy,ix),DEN(iz+1,iy,ix),vapsln,gamma,delta,C44ZP1)
              call htic6rc44(VP(iz+1,iy+1,ix),VS(iz+1,iy+1,ix),DEN(iz+1,iy+1,ix),vapsln,gamma,delta,C44ZYP1)
              call htic6rc55(VP(iz,iy,ix+1),VS(iz,iy,ix+1),DEN(iz,iy,ix+1),vapsln,gamma,delta,C55XP1)
              call htic6rc55(VP(iz+1,iy,ix),VS(iz+1,iy,ix),DEN(iz+1,iy,ix),vapsln,gamma,delta,C55ZP1)
              call htic6rc55(VP(iz+1,iy,ix+1),VS(iz+1,iy,ix+1),DEN(iz+1,iy,ix+1),vapsln,gamma,delta,C55ZXP1)
              C66XP1=C55XP1
              call htic6rc66(VP(iz,iy+1,ix),VS(iz,iy+1,ix),DEN(iz,iy+1,ix),vapsln,gamma,delta,C66YP1)
              call htic6rc66(VP(iz,iy+1,ix+1),VS(iz,iy+1,ix+1),DEN(iz,iy+1,ix+1),vapsln,gamma,delta,C66YXP1)
          end if

          !interpolate elastic parameter
          c44yz=4.0_8/(1.0_8/C44+1.0_8/C44YP1+1.0_8/C44ZP1+1.0_8/C44ZYP1)
          c55xz=4.0_8/(1.0_8/C55+1.0_8/C55XP1+1.0_8/C55ZP1+1.0_8/C55ZXP1)
          c66xy=4.0_8/(1.0_8/C66+1.0_8/C66XP1+1.0_8/C66YP1+1.0_8/C66YXP1)

          !calculate stress of the next time step
          TXX(iz,iy,ix)=(TXX(iz,iy,ix)+C11*VX_x*dtdx+C12*VY_y*dtdy+C13*VZ_z*dtdz)*ABR(iz,iy,ix)
          TYY(iz,iy,ix)=(TYY(iz,iy,ix)+C12*VX_x*dtdx+C22*VY_y*dtdy+C23*VZ_z*dtdz)*ABR(iz,iy,ix)
          TZZ(iz,iy,ix)=(TZZ(iz,iy,ix)+C13*VX_x*dtdx+C23*VY_y*dtdy+C33*VZ_z*dtdz)*ABR(iz,iy,ix)
          TYZ(iz,iy,ix)=(TYZ(iz,iy,ix)+c44yz*(VY_z*dtdz+VZ_y*dtdy))*ABR(iz,iy,ix)
          TXZ(iz,iy,ix)=(TXZ(iz,iy,ix)+c55xz*(VZ_x*dtdx+VX_z*dtdz))*ABR(iz,iy,ix)
          TXY(iz,iy,ix)=(TXY(iz,iy,ix)+c66xy*(VY_x*dtdx+VX_y*dtdy))*ABR(iz,iy,ix)
        end do
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    !apply temporal derivative to source time function. Here we use incremental stress as source excitation, and apply on velocity-stress scheme.
    !so the temporal derivative of source time function is used. Reference ("The finite-difference modelling of earthquake motions" Moczo p245) &
    !("3D elastic finite-difference modelling of seismic motion using staggered grids with nonuniform spacing" Arben Pitarka 1999).
    if (it==1) then
      ssignal=Soutf(2)
    else if (it<NT) then
      ssignal=Soutf(it+1)-Soutf(it-1)
    else
      ssignal=-Soutf(it-1)
    end if

    !add source
    if (sctp==1) then
      !compressional source
      TXX(SNZ,SNY,SNX)=TXX(SNZ,SNY,SNX)-ssignal
      TYY(SNZ,SNY,SNX)=TYY(SNZ,SNY,SNX)-ssignal
      TZZ(SNZ,SNY,SNX)=TZZ(SNZ,SNY,SNX)-ssignal
    else
      !moment tensor source
      TXX(SNZ,SNY,SNX)=TXX(SNZ,SNY,SNX)-Mxx*ssignal
      TYY(SNZ,SNY,SNX)=TYY(SNZ,SNY,SNX)-Myy*ssignal
      TZZ(SNZ,SNY,SNX)=TZZ(SNZ,SNY,SNX)-Mzz*ssignal
      TXY(SNZ,SNY,SNX)=TXY(SNZ,SNY,SNX)-0.25_8*Mxy*ssignal
      TXY(SNZ,SNY,SNX-1)=TXY(SNZ,SNY,SNX-1)-0.25_8*Mxy*ssignal
      TXY(SNZ,SNY-1,SNX)=TXY(SNZ,SNY-1,SNX)-0.25_8*Mxy*ssignal
      TXY(SNZ,SNY-1,SNX-1)=TXY(SNZ,SNY-1,SNX-1)-0.25_8*Mxy*ssignal
      TXZ(SNZ,SNY,SNX)=TXZ(SNZ,SNY,SNX)-0.25_8*Mxz*ssignal
      TXZ(SNZ,SNY,SNX-1)=TXZ(SNZ,SNY,SNX-1)-0.25_8*Mxz*ssignal
      TXZ(SNZ-1,SNY,SNX)=TXZ(SNZ-1,SNY,SNX)-0.25_8*Mxz*ssignal
      TXZ(SNZ-1,SNY,SNX-1)=TXZ(SNZ-1,SNY,SNX-1)-0.25_8*Mxz*ssignal
      TYZ(SNZ,SNY,SNX)=TYZ(SNZ,SNY,SNX)-0.25_8*Myz*ssignal
      TYZ(SNZ-1,SNY,SNX)=TYZ(SNZ-1,SNY,SNX)-0.25_8*Myz*ssignal
      TYZ(SNZ,SNY-1,SNX)=TYZ(SNZ,SNY-1,SNX)-0.25_8*Myz*ssignal
      TYZ(SNZ-1,SNY-1,SNX)=TYZ(SNZ-1,SNY-1,SNX)-0.25_8*Myz*ssignal
    end if

    !write received trace data
    do ii=1,NREC
      write(101) VX(RECZ(ii),RECY(ii),RECX(ii))
      write(102) VY(RECZ(ii),RECY(ii),RECX(ii))
      write(103) VZ(RECZ(ii),RECY(ii),RECX(ii))
    end do

    !write wavefield snapshot profile
    if (mod(it-1,snapdt)==0) then
      if (fspx==1) then
        write(111) VX(:,:,snapx)
        write(112) VY(:,:,snapx)
        write(113) VZ(:,:,snapx)
      end if
      if (fspy==1) then
        write(114) VX(:,snapy,:)
        write(115) VY(:,snapy,:)
        write(116) VZ(:,snapy,:)
      end if
      if (fspz==1) then
        write(117) VX(snapz,:,:)
        write(118) VY(snapz,:,:)
        write(119) VZ(snapz,:,:)
      end if
    end if

    !write wavefield snapshot of the whole wavefield for vx, vy and vz compoment
    if (mod(it-1,snhmdt)==0) then
      write(namchar,*) it-1
      if (fgvx==1) then
        open(121,file='vx_'//trim(adjustl(namchar))//'.snap',form='unformatted',status='replace',access='stream',action='write')
        write(121) VX
        close(121)
      end if
      if (fgvy==1) then
        open(122,file='vy_'//trim(adjustl(namchar))//'.snap',form='unformatted',status='replace',access='stream',action='write')
        write(122) VY
        close(122)
      end if
      if (fgvz==1) then
        open(123,file='vz_'//trim(adjustl(namchar))//'.snap',form='unformatted',status='replace',access='stream',action='write')
        write(123) VZ
        close(123)
      end if
    end if

    write(*,*) "Time step:", it
  end do!end loop for time steps
  call cpu_time(tend)
  write(*,*) "Computing time for finite-difference wave equation modeling:", (tend-tstart)/60.0, "minute."

  !close opened received trace data file
  close(unit=101)
  close(unit=102)
  close(unit=103)

  !close opened wavefield snapshot profile file
  if (fspx==1) then
    close(unit=111)
    close(unit=112)
    close(unit=113)
  end if
  if (fspy==1) then
    close(unit=114)
    close(unit=115)
    close(unit=116)
  end if
  if (fspz==1) then
    close(unit=117)
    close(unit=118)
    close(unit=119)
  end if

  !free memory
  deallocate(VP,VS,DEN)
  deallocate(TXX,TYY,TZZ,TXY,TXZ,TYZ,VX,VY,VZ,ABR,Soutf)
  deallocate(RECX,RECY,RECZ)

contains
!internal subroutines
  subroutine Ricker()
  !generate Ricker wavelet
  use paramod
  implicit none
    real(kind=RLP) :: factor,t0,ft
    integer(kind=INP) :: i
    t0 = 1.1_8/freq
    write(*,*) 'Using Ricker wavelet as souce time function, time delay: ', t0, 's'

    do i=0,NT-1
       ft = dt*i-t0
       factor = (PI*freq*ft)*(PI*freq*ft)
       Soutf(i+1) = (1.0_8-2.0_8*factor)*exp(-factor)         
    end do
    open(unit=222,file='ricker.o',form='unformatted',status='replace',access='stream',action='write')
    write(222) Soutf
    close(unit=222)
  end subroutine Ricker

  subroutine  absorbcdt()
  !calculate damp factor of absorbing boundary condition
  use paramod
  implicit none
    integer(kind=INP) :: icx,icy,icz
    real(kind=RLP) :: afc,cdamp(boudp)

    afc=log(0.92)
    do icx=1,boudp
      cdamp(icx)=exp(afc*(boudp-icx+1)*(boudp-icx+1)/(boudp*boudp))
    end do

    do icx=1,boudp
      do icy=1,NY
        do icz=1,NZ
          ABR(icz,icy,icx)=ABR(icz,icy,icx)*cdamp(icx)
          ABR(icz,icy,NX-icx+1)=ABR(icz,icy,NX-icx+1)*cdamp(icx)
        end do
      end do
    end do

    do icx=1,NX
      do icy=1,boudp
        do icz=1,NZ
          ABR(icz,icy,icx)=ABR(icz,icy,icx)*cdamp(icy)
          ABR(icz,NY-icy+1,icx)=ABR(icz,NY-icy+1,icx)*cdamp(icy)
        end do
      end do
    end do

    do icx=1,NX
      do icy=1,NY
        do icz=1,boudp
          ABR(icz,icy,icx)=ABR(icz,icy,icx)*cdamp(icz)
          ABR(NZ-icz+1,icy,icx)=ABR(NZ-icz+1,icy,icx)*cdamp(icz)
        end do
      end do
    end do
  end subroutine absorbcdt

end program fdmodeling 

!external subroutines
subroutine htic6r(vp,vs,den,vapsln,gamma,delta,c11,c12,c13,c22,c23,c33,c44,c55,c66)
!obtain elastic tensor of HTI medium by rotate the VTI medium.
!the HTI medium are constructed by rotating VTI medium anticlockwise(??) along Y-axis by pi/2.
!note here: the 'vapsln, gamma, delta' means the anisotropic parameters in the corresponding VTI medium.
use paramod
implicit none
real(kind=RLP),intent(in)  :: vp,vs,den,vapsln,gamma,delta
real(kind=RLP),intent(out) :: c11,c12,c13,c22,c23,c33,c44,c55,c66
real(kind=RLP) :: rotang, c6(6,6), rc3(3,3,3,3), c3(3,3,3,3)
real(kind=RLP) :: rotmax(3,3)=0.0
integer        :: m,n,r,s,i,j,k,l

  rc3=0.0
  c3=0.0
  c6=0.0

  ! set the rot angle to be pi/2
  rotang=atan(1.0)*2.0

  !first use the input parameters to construct elastic tensor of VTI medium
  c6(3,3)=vp*vp*den
  c6(4,4)=vs*vs*den
  c6(1,1)=2.0*vapsln*c6(3,3)+c6(3,3)
  c6(6,6)=2.0*gamma*c6(4,4)+c6(4,4)
  c6(1,3)=sqrt(2.0*delta*c6(3,3)*(c6(3,3)-c6(4,4))+(c6(3,3)-c6(4,4))*(c6(3,3)-c6(4,4)))-c6(4,4)
  c6(2,2)=c6(1,1)
  c6(5,5)=c6(4,4)
  c6(1,2)=c6(1,1)-2.0*c6(6,6)
  c6(2,3)=c6(1,3)
  c6(2,1)=c6(1,2)
  c6(3,1)=c6(1,3)
  c6(3,2)=c6(2,3)
  !then calculate the corresponding fourth-order elastic tensor of this VTI medium
  rc3(1,1,1,1)=c6(1,1)
  rc3(2,2,2,2)=c6(2,2)
  rc3(3,3,3,3)=c6(3,3)

  rc3(2,3,2,3)=c6(4,4)
  rc3(3,2,2,3)=c6(4,4)
  rc3(2,3,3,2)=c6(4,4)
  rc3(3,2,3,2)=c6(4,4)

  rc3(1,3,1,3)=c6(5,5)
  rc3(3,1,1,3)=c6(5,5)
  rc3(1,3,3,1)=c6(5,5)
  rc3(3,1,3,1)=c6(5,5)

  rc3(1,2,1,2)=c6(6,6)
  rc3(2,1,1,2)=c6(6,6)
  rc3(1,2,2,1)=c6(6,6)
  rc3(2,1,2,1)=c6(6,6)

  rc3(1,1,2,2)=c6(1,2)
  rc3(1,1,3,3)=c6(1,3)
  rc3(2,2,3,3)=c6(2,3)

  rc3(2,2,1,1)=c6(2,1)
  rc3(3,3,1,1)=c6(3,1)
  rc3(3,3,2,2)=c6(3,2)

  !construct the rotation matrix
  rotmax(1,1)=cos(rotang)
  rotmax(1,3)=-sin(rotang)
  rotmax(2,2)=1.0
  rotmax(3,1)=sin(rotang)
  rotmax(3,3)=cos(rotang)

  !then rotate VTI elastic tensor to TTI elastic tensor
  do s=1,3,1
    do r=1,3,1
      do n=1,3,1
        do m=1,3,1
          do l=1,3,1
            do k=1,3,1
              do j=1,3,1
                do i=1,3,1
                  c3(m,n,r,s)=c3(m,n,r,s)+rotmax(m,i)*rotmax(n,j)*rotmax(r,k)*rotmax(s,l)*rc3(i,j,k,l)
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end do

  c11=c3(1,1,1,1)
  c12=c3(1,1,2,2)
  c13=c3(1,1,3,3)
  c22=c3(2,2,2,2)
  c23=c3(2,2,3,3)
  c33=c3(3,3,3,3)
  c44=c3(2,3,2,3)
  c55=c3(1,3,1,3)
  c66=c3(1,2,1,2)
end subroutine htic6r

subroutine htic6rc44(vp,vs,den,vapsln,gamma,delta,c44)
!obtain elastic tensor of HTI medium by rotate the VTI medium.
!the HTI medium are constructed by rotating VTI medium anticlockwise(??) along Y-axis by pi/2.
!note here: the 'vapsln, gamma, delta' means the anisotropic parameters in the corresponding VTI medium.
!only calculate c44
use paramod
implicit none
real(kind=RLP),intent(in)  :: vp,vs,den,vapsln,gamma,delta
real(kind=RLP),intent(out) :: c44
real(kind=RLP) :: c11,c12,c13,c22,c23,c33,c55,c66
real(kind=RLP) :: rotang, c6(6,6), rc3(3,3,3,3), c3(3,3,3,3)
real(kind=RLP) :: rotmax(3,3)=0.0
integer        :: m,n,r,s,i,j,k,l

  rc3=0.0
  c3=0.0
  c6=0.0

  ! set the rot angle to be pi/2
  rotang=atan(1.0)*2.0

  !first use the input parameters to construct elastic tensor of VTI medium
  c6(3,3)=vp*vp*den
  c6(4,4)=vs*vs*den
  c6(1,1)=2.0*vapsln*c6(3,3)+c6(3,3)
  c6(6,6)=2.0*gamma*c6(4,4)+c6(4,4)
  c6(1,3)=sqrt(2.0*delta*c6(3,3)*(c6(3,3)-c6(4,4))+(c6(3,3)-c6(4,4))*(c6(3,3)-c6(4,4)))-c6(4,4)
  c6(2,2)=c6(1,1)
  c6(5,5)=c6(4,4)
  c6(1,2)=c6(1,1)-2.0*c6(6,6)
  c6(2,3)=c6(1,3)
  c6(2,1)=c6(1,2)
  c6(3,1)=c6(1,3)
  c6(3,2)=c6(2,3)
  !then calculate the corresponding fourth-order elastic tensor of this VTI medium
  rc3(1,1,1,1)=c6(1,1)
  rc3(2,2,2,2)=c6(2,2)
  rc3(3,3,3,3)=c6(3,3)

  rc3(2,3,2,3)=c6(4,4)
  rc3(3,2,2,3)=c6(4,4)
  rc3(2,3,3,2)=c6(4,4)
  rc3(3,2,3,2)=c6(4,4)

  rc3(1,3,1,3)=c6(5,5)
  rc3(3,1,1,3)=c6(5,5)
  rc3(1,3,3,1)=c6(5,5)
  rc3(3,1,3,1)=c6(5,5)

  rc3(1,2,1,2)=c6(6,6)
  rc3(2,1,1,2)=c6(6,6)
  rc3(1,2,2,1)=c6(6,6)
  rc3(2,1,2,1)=c6(6,6)

  rc3(1,1,2,2)=c6(1,2)
  rc3(1,1,3,3)=c6(1,3)
  rc3(2,2,3,3)=c6(2,3)

  rc3(2,2,1,1)=c6(2,1)
  rc3(3,3,1,1)=c6(3,1)
  rc3(3,3,2,2)=c6(3,2)

  !construct the rotation matrix
  rotmax(1,1)=cos(rotang)
  rotmax(1,3)=-sin(rotang)
  rotmax(2,2)=1.0
  rotmax(3,1)=sin(rotang)
  rotmax(3,3)=cos(rotang)

  !then rotate VTI elastic tensor to TTI elastic tensor
  do s=1,3,1
    do r=1,3,1
      do n=1,3,1
        do m=1,3,1
          do l=1,3,1
            do k=1,3,1
              do j=1,3,1
                do i=1,3,1
                  c3(m,n,r,s)=c3(m,n,r,s)+rotmax(m,i)*rotmax(n,j)*rotmax(r,k)*rotmax(s,l)*rc3(i,j,k,l)
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end do

  c11=c3(1,1,1,1)
  c12=c3(1,1,2,2)
  c13=c3(1,1,3,3)
  c22=c3(2,2,2,2)
  c23=c3(2,2,3,3)
  c33=c3(3,3,3,3)
  c44=c3(2,3,2,3)
  c55=c3(1,3,1,3)
  c66=c3(1,2,1,2)
end subroutine htic6rc44

subroutine htic6rc55(vp,vs,den,vapsln,gamma,delta,c55)
!obtain elastic tensor of HTI medium by rotate the VTI medium.
!the HTI medium are constructed by rotating VTI medium anticlockwise(??) along Y-axis by pi/2.
!note here: the 'vapsln, gamma, delta' means the anisotropic parameters in the corresponding VTI medium.
!only calculate c55
use paramod
implicit none
real(kind=RLP),intent(in)  :: vp,vs,den,vapsln,gamma,delta
real(kind=RLP),intent(out) :: c55
real(kind=RLP) :: c11,c12,c13,c22,c23,c33,c44,c66
real(kind=RLP) :: rotang, c6(6,6), rc3(3,3,3,3), c3(3,3,3,3)
real(kind=RLP) :: rotmax(3,3)=0.0
integer        :: m,n,r,s,i,j,k,l

  rc3=0.0
  c3=0.0
  c6=0.0

  ! set the rot angle to be pi/2
  rotang=atan(1.0)*2.0

  !first use the input parameters to construct elastic tensor of VTI medium
  c6(3,3)=vp*vp*den
  c6(4,4)=vs*vs*den
  c6(1,1)=2.0*vapsln*c6(3,3)+c6(3,3)
  c6(6,6)=2.0*gamma*c6(4,4)+c6(4,4)
  c6(1,3)=sqrt(2.0*delta*c6(3,3)*(c6(3,3)-c6(4,4))+(c6(3,3)-c6(4,4))*(c6(3,3)-c6(4,4)))-c6(4,4)
  c6(2,2)=c6(1,1)
  c6(5,5)=c6(4,4)
  c6(1,2)=c6(1,1)-2.0*c6(6,6)
  c6(2,3)=c6(1,3)
  c6(2,1)=c6(1,2)
  c6(3,1)=c6(1,3)
  c6(3,2)=c6(2,3)
  !then calculate the corresponding fourth-order elastic tensor of this VTI medium
  rc3(1,1,1,1)=c6(1,1)
  rc3(2,2,2,2)=c6(2,2)
  rc3(3,3,3,3)=c6(3,3)

  rc3(2,3,2,3)=c6(4,4)
  rc3(3,2,2,3)=c6(4,4)
  rc3(2,3,3,2)=c6(4,4)
  rc3(3,2,3,2)=c6(4,4)

  rc3(1,3,1,3)=c6(5,5)
  rc3(3,1,1,3)=c6(5,5)
  rc3(1,3,3,1)=c6(5,5)
  rc3(3,1,3,1)=c6(5,5)

  rc3(1,2,1,2)=c6(6,6)
  rc3(2,1,1,2)=c6(6,6)
  rc3(1,2,2,1)=c6(6,6)
  rc3(2,1,2,1)=c6(6,6)

  rc3(1,1,2,2)=c6(1,2)
  rc3(1,1,3,3)=c6(1,3)
  rc3(2,2,3,3)=c6(2,3)

  rc3(2,2,1,1)=c6(2,1)
  rc3(3,3,1,1)=c6(3,1)
  rc3(3,3,2,2)=c6(3,2)

  !construct the rotation matrix
  rotmax(1,1)=cos(rotang)
  rotmax(1,3)=-sin(rotang)
  rotmax(2,2)=1.0
  rotmax(3,1)=sin(rotang)
  rotmax(3,3)=cos(rotang)

  !then rotate VTI elastic tensor to TTI elastic tensor
  do s=1,3,1
    do r=1,3,1
      do n=1,3,1
        do m=1,3,1
          do l=1,3,1
            do k=1,3,1
              do j=1,3,1
                do i=1,3,1
                  c3(m,n,r,s)=c3(m,n,r,s)+rotmax(m,i)*rotmax(n,j)*rotmax(r,k)*rotmax(s,l)*rc3(i,j,k,l)
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end do

  c11=c3(1,1,1,1)
  c12=c3(1,1,2,2)
  c13=c3(1,1,3,3)
  c22=c3(2,2,2,2)
  c23=c3(2,2,3,3)
  c33=c3(3,3,3,3)
  c44=c3(2,3,2,3)
  c55=c3(1,3,1,3)
  c66=c3(1,2,1,2)
end subroutine htic6rc55

subroutine htic6rc66(vp,vs,den,vapsln,gamma,delta,c66)
!obtain elastic tensor of HTI medium by rotate the VTI medium.
!the HTI medium are constructed by rotating VTI medium anticlockwise(??) along Y-axis by pi/2.
!note here: the 'vapsln, gamma, delta' means the anisotropic parameters in the corresponding VTI medium.
!only calculate c66
use paramod
implicit none
real(kind=RLP),intent(in)  :: vp,vs,den,vapsln,gamma,delta
real(kind=RLP),intent(out) :: c66
real(kind=RLP) :: c11,c12,c13,c22,c23,c33,c44,c55
real(kind=RLP) :: rotang, c6(6,6), rc3(3,3,3,3), c3(3,3,3,3)
real(kind=RLP) :: rotmax(3,3)=0.0
integer        :: m,n,r,s,i,j,k,l

  rc3=0.0
  c3=0.0
  c6=0.0

  ! set the rot angle to be pi/2
  rotang=atan(1.0)*2.0

  !first use the input parameters to construct elastic tensor of VTI medium
  c6(3,3)=vp*vp*den
  c6(4,4)=vs*vs*den
  c6(1,1)=2.0*vapsln*c6(3,3)+c6(3,3)
  c6(6,6)=2.0*gamma*c6(4,4)+c6(4,4)
  c6(1,3)=sqrt(2.0*delta*c6(3,3)*(c6(3,3)-c6(4,4))+(c6(3,3)-c6(4,4))*(c6(3,3)-c6(4,4)))-c6(4,4)
  c6(2,2)=c6(1,1)
  c6(5,5)=c6(4,4)
  c6(1,2)=c6(1,1)-2.0*c6(6,6)
  c6(2,3)=c6(1,3)
  c6(2,1)=c6(1,2)
  c6(3,1)=c6(1,3)
  c6(3,2)=c6(2,3)
  !then calculate the corresponding fourth-order elastic tensor of this VTI medium
  rc3(1,1,1,1)=c6(1,1)
  rc3(2,2,2,2)=c6(2,2)
  rc3(3,3,3,3)=c6(3,3)

  rc3(2,3,2,3)=c6(4,4)
  rc3(3,2,2,3)=c6(4,4)
  rc3(2,3,3,2)=c6(4,4)
  rc3(3,2,3,2)=c6(4,4)

  rc3(1,3,1,3)=c6(5,5)
  rc3(3,1,1,3)=c6(5,5)
  rc3(1,3,3,1)=c6(5,5)
  rc3(3,1,3,1)=c6(5,5)

  rc3(1,2,1,2)=c6(6,6)
  rc3(2,1,1,2)=c6(6,6)
  rc3(1,2,2,1)=c6(6,6)
  rc3(2,1,2,1)=c6(6,6)

  rc3(1,1,2,2)=c6(1,2)
  rc3(1,1,3,3)=c6(1,3)
  rc3(2,2,3,3)=c6(2,3)

  rc3(2,2,1,1)=c6(2,1)
  rc3(3,3,1,1)=c6(3,1)
  rc3(3,3,2,2)=c6(3,2)

  !construct the rotation matrix
  rotmax(1,1)=cos(rotang)
  rotmax(1,3)=-sin(rotang)
  rotmax(2,2)=1.0
  rotmax(3,1)=sin(rotang)
  rotmax(3,3)=cos(rotang)

  !then rotate VTI elastic tensor to TTI elastic tensor
  do s=1,3,1
    do r=1,3,1
      do n=1,3,1
        do m=1,3,1
          do l=1,3,1
            do k=1,3,1
              do j=1,3,1
                do i=1,3,1
                  c3(m,n,r,s)=c3(m,n,r,s)+rotmax(m,i)*rotmax(n,j)*rotmax(r,k)*rotmax(s,l)*rc3(i,j,k,l)
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end do

  c11=c3(1,1,1,1)
  c12=c3(1,1,2,2)
  c13=c3(1,1,3,3)
  c22=c3(2,2,2,2)
  c23=c3(2,2,3,3)
  c33=c3(3,3,3,3)
  c44=c3(2,3,2,3)
  c55=c3(1,3,1,3)
  c66=c3(1,2,1,2)
end subroutine htic6rc66