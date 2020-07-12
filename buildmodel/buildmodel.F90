program buildmodel
!this program is used to build a layered model for wave equation modeling
!each layer could be an isotropic, VTI or HTI medium
  use paramod
  implicit none
  integer(kind=INP) :: NX,NY,NZ,NLY,nzp
  integer(kind=INP),allocatable,dimension(:) :: DPTH,ANITP
  real(kind=RLP),allocatable,dimension(:) :: VP,VS,RHO,vapsln,gamma,delta
  real(kind=RLP),allocatable,dimension(:,:,:) :: C11,C12,C13,C22,C23,C33,C44,C55,C66,DEN
  integer(kind=INP) :: ii,igz

  open(unit=100,file="model.dat")
  read(100,*) NX, NY, NZ
  read(100,*) NLY

  allocate(DPTH(1:NLY+1),VP(NLY),VS(NLY),RHO(NLY),vapsln(NLY),gamma(NLY),delta(NLY),ANITP(NLY))
  allocate(C11(NZ,NY,NX),C12(NZ,NY,NX),C13(NZ,NY,NX),C22(NZ,NY,NX),C23(NZ,NY,NX),C33(NZ,NY,NX),C44(NZ,NY,NX),C55(NZ,NY,NX),C66(NZ,NY,NX),DEN(NZ,NY,NX))

  read(100,*) DPTH(2:NLY)
  read(100,*) VP
  read(100,*) VS
  read(100,*) RHO
  read(100,*) vapsln
  read(100,*) gamma
  read(100,*) delta
  read(100,*) ANITP
  close(unit=100)

  DPTH(1)=1
  DPTH(NLY+1)=NZ
  write(*,*) "Model size: NX=", NX, ", NY=", NY, ", NZ=", NZ
  write(*,*) "Number of layers:", NLY
  write(*,*) "Bounday of layers:", DPTH
  DPTH(1)=0
  do ii=1,NLY
    write(*,*) "--------------------------------------------------------------"
    write(*,*) "Layer", ii, "depth:", DPTH(ii)+1, "--", DPTH(ii+1)
    write(*,*) "VP=", VP(ii), "VS=", VS(ii), "DEN=", RHO(ii)
    write(*,*) "vapsln=", vapsln(ii), "gamma=", gamma(ii), "delta=", delta(ii)
    igz=DPTH(ii)+1
    if (ANITP(ii)==0) then
      call isoc6(VP(ii),VS(ii),RHO(ii),C11(igz,1,1),C12(igz,1,1),C13(igz,1,1),C22(igz,1,1),C23(igz,1,1),C33(igz,1,1),C44(igz,1,1),C55(igz,1,1),C66(igz,1,1))
      write(*,*) "This layer is an isotropic medium."
      write(*,*) "C11=", C11(igz,1,1)
      write(*,*) "C12=", C12(igz,1,1)
      write(*,*) "C13=", C13(igz,1,1)
      write(*,*) "C22=", C22(igz,1,1)
      write(*,*) "C23=", C23(igz,1,1)
      write(*,*) "C33=", C33(igz,1,1)
      write(*,*) "C44=", C44(igz,1,1)
      write(*,*) "C55=", C55(igz,1,1)
      write(*,*) "C66=", C66(igz,1,1)
    else if(ANITP(ii)==1) then
      call vtic6(VP(ii),VS(ii),RHO(ii),vapsln(ii),gamma(ii),delta(ii),C11(igz,1,1),C12(igz,1,1),C13(igz,1,1),C22(igz,1,1),C23(igz,1,1),C33(igz,1,1),C44(igz,1,1),C55(igz,1,1),C66(igz,1,1))
      write(*,*) "This layer is a VTI medium."
      write(*,*) "C11=", C11(igz,1,1)
      write(*,*) "C12=", C12(igz,1,1)
      write(*,*) "C13=", C13(igz,1,1)
      write(*,*) "C22=", C22(igz,1,1)
      write(*,*) "C23=", C23(igz,1,1)
      write(*,*) "C33=", C33(igz,1,1)
      write(*,*) "C44=", C44(igz,1,1)
      write(*,*) "C55=", C55(igz,1,1)
      write(*,*) "C66=", C66(igz,1,1)
    else if(ANITP(ii)==2) then
      call htic6r(VP(ii),VS(ii),RHO(ii),vapsln(ii),gamma(ii),delta(ii),C11(igz,1,1),C12(igz,1,1),C13(igz,1,1),C22(igz,1,1),C23(igz,1,1),C33(igz,1,1),C44(igz,1,1),C55(igz,1,1),C66(igz,1,1))
      write(*,*) "This layer is a HTI medium."
      write(*,*) "C11=", C11(igz,1,1)
      write(*,*) "C12=", C12(igz,1,1)
      write(*,*) "C13=", C13(igz,1,1)
      write(*,*) "C22=", C22(igz,1,1)
      write(*,*) "C23=", C23(igz,1,1)
      write(*,*) "C33=", C33(igz,1,1)
      write(*,*) "C44=", C44(igz,1,1)
      write(*,*) "C55=", C55(igz,1,1)
      write(*,*) "C66=", C66(igz,1,1)
    else
      stop "Error! Incorrect input for anisotropy type!"
    end if
    C11(igz:DPTH(ii+1),:,:)=C11(igz,1,1)
    C12(igz:DPTH(ii+1),:,:)=C12(igz,1,1)
    C13(igz:DPTH(ii+1),:,:)=C13(igz,1,1)
    C22(igz:DPTH(ii+1),:,:)=C22(igz,1,1)
    C23(igz:DPTH(ii+1),:,:)=C23(igz,1,1)
    C33(igz:DPTH(ii+1),:,:)=C33(igz,1,1)
    C44(igz:DPTH(ii+1),:,:)=C44(igz,1,1)
    C55(igz:DPTH(ii+1),:,:)=C55(igz,1,1)
    C66(igz:DPTH(ii+1),:,:)=C66(igz,1,1)
    DEN(igz:DPTH(ii+1),:,:)=RHO(ii)
    write(*,*)
  end do

  open(unit=11,file='C11',form='unformatted',status='replace',access='stream')
  open(unit=12,file='C12',form='unformatted',status='replace',access='stream')
  open(unit=13,file='C13',form='unformatted',status='replace',access='stream')
  open(unit=22,file='C22',form='unformatted',status='replace',access='stream')
  open(unit=23,file='C23',form='unformatted',status='replace',access='stream')
  open(unit=33,file='C33',form='unformatted',status='replace',access='stream')
  open(unit=44,file='C44',form='unformatted',status='replace',access='stream')
  open(unit=55,file='C55',form='unformatted',status='replace',access='stream')
  open(unit=66,file='C66',form='unformatted',status='replace',access='stream')
  open(unit=88,file='DEN',form='unformatted',status='replace',access='stream')

  write(unit=11) C11
  write(unit=12) C12
  write(unit=13) C13
  write(unit=22) C22
  write(unit=23) C23
  write(unit=33) C33
  write(unit=44) C44
  write(unit=55) C55
  write(unit=66) C66
  write(unit=88) DEN

  close(unit=11)
  close(unit=12)
  close(unit=13)
  close(unit=22)
  close(unit=23)
  close(unit=33)
  close(unit=44)
  close(unit=55)
  close(unit=66)
  close(unit=88)
  deallocate(DPTH,VP,VS,RHO,vapsln,gamma,delta)
  deallocate(C11,C12,C13,C22,C23,C33,C44,C55,C66,DEN)
end program buildmodel

subroutine isoc6(vp,vs,den,c11,c12,c13,c22,c23,c33,c44,c55,c66)
!this subroutine is used to calculate the elastic tensor in an isotropic medium
  use paramod
  implicit none
  real(kind=RLP),intent(in)  :: vp,vs,den
  real(kind=RLP),intent(out) :: c11,c12,c13,c22,c23,c33,c44,c55,c66

  c11=vp*vp*den
  c22=c11
  c33=c11
  c44=vs*vs*den
  c55=c44
  c66=c44
  c12=c11-2.0_8*c44
  c13=c12
  c23=c12
end subroutine isoc6

subroutine vtic6(vp,vs,den,vapsln,gamma,delta,c11,c12,c13,c22,c23,c33,c44,c55,c66)
!this subroutine is used to calculate the elastic tensor in a VTI medium
  use paramod
  implicit none
  real(kind=RLP),intent(in)  :: vp,vs,den,vapsln,gamma,delta
  real(kind=RLP),intent(out) :: c11,c12,c13,c22,c23,c33,c44,c55,c66

  c33=vp*vp*den
  c44=vs*vs*den
  c11=2.0_8*vapsln*c33+c33
  c66=2.0_8*gamma*c44+c44
  c13=sqrt(2.0_8*delta*c33*(c33-c44)+(c33-c44)*(c33-c44))-c44
  c22=c11
  c55=c44
  c12=c11-2.0_8*c66
  c23=c13
end subroutine vtic6

subroutine htic6(vp,vs,den,vapsln,gamma,delta,c11,c12,c13,c22,c23,c33,c44,c55,c66)
!this subroutine is used to calculate the elastic tensor in a HTI medium
  use paramod
  implicit none
  real(kind=RLP),intent(in)  :: vp,vs,den,vapsln,gamma,delta
  real(kind=RLP),intent(out) :: c11,c12,c13,c22,c23,c33,c44,c55,c66

  c33=vp*vp*den
  c44=vs*vs*den
  c11=2.0_8*vapsln*c33+c33
  c66=2.0_8*gamma*c44+c44
  c13=sqrt(2.0_8*delta*c33*(c33-c44)+(c33-c44)*(c33-c44))-c44
  c22=c33
  c55=c66
  c12=c13
  c23=c33-2.0_8*c44
end subroutine htic6

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

subroutine c3tovoigt(c6,c3)
!transform the standard fourth-order elastic tensor c_ijkl to voigt form c_mn
use paramod
implicit none
real(kind=RLP),intent(in)  :: c3(3,3,3,3)
real(kind=RLP),intent(out) :: c6(6,6)

  !first calculate the upper triangular part of c6
  c6(1,1)=c3(1,1,1,1)
  c6(1,2)=c3(1,1,2,2)
  c6(1,3)=c3(1,1,3,3)
  c6(1,4)=c3(1,1,2,3)
  c6(1,5)=c3(1,1,1,3)
  c6(1,6)=c3(1,1,1,2)
  c6(2,2)=c3(2,2,2,2)
  c6(2,3)=c3(2,2,3,3)
  c6(2,4)=c3(2,2,2,3)
  c6(2,5)=c3(2,2,1,3)
  c6(2,6)=c3(2,2,1,2)
  c6(3,3)=c3(3,3,3,3)
  c6(3,4)=c3(3,3,2,3)
  c6(3,5)=c3(3,3,1,3)
  c6(3,6)=c3(3,3,1,2)
  c6(4,4)=c3(2,3,2,3)
  c6(4,5)=c3(2,3,1,3)
  c6(4,6)=c3(2,3,1,2)
  c6(5,5)=c3(1,3,1,3)
  c6(5,6)=c3(1,3,1,2)
  c6(6,6)=c3(1,2,1,2)
  !then apply the symmetry of c6
  c6(2,1)=c6(1,2)
  c6(3,1)=c6(1,3)
  c6(3,2)=c6(2,3)
  c6(4,1)=c6(1,4)
  c6(4,2)=c6(2,4)
  c6(4,3)=c6(3,4)
  c6(5,1)=c6(1,5)
  c6(5,2)=c6(2,5)
  c6(5,3)=c6(3,5)
  c6(5,4)=c6(4,5)
  c6(6,1)=c6(1,6)
  c6(6,2)=c6(2,6)
  c6(6,3)=c6(3,6)
  c6(6,4)=c6(4,6)
  c6(6,5)=c6(5,6)
end subroutine c3tovoigt
