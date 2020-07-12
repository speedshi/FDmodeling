program buildreceiver
!this program is used to build the receiver file for wave equation modeling
!receiver can be placed at different zones.
  use paramod
  implicit none
  integer(kind=INP) :: NC,ii,NRC,ix,iy,iz
  integer(kind=INP),allocatable,dimension(:) :: RXS,RXE,RXD,RYS,RYE,RYD,RZS,RZE,RZD

  open(unit=11,file="recein.dat")
  read(11,*) NC

  allocate(RXS(NC),RXE(NC),RXD(NC),RYS(NC),RYE(NC),RYD(NC),RZS(NC),RZE(NC),RZD(NC))

  do ii=1,NC
    read(11,*) RXS(ii),RXE(ii),RXD(ii)
    read(11,*) RYS(ii),RYE(ii),RYD(ii)
    read(11,*) RZS(ii),RZE(ii),RZD(ii)
  end do
  close(unit=11)

  NRC=0
  do ii=1,NC
    NRC=NRC+((RXE(ii)-RXS(ii))/RXD(ii)+1)*((RYE(ii)-RYS(ii))/RYD(ii)+1)*((RZE(ii)-RZS(ii))/RZD(ii)+1)
  end do

  open(unit=12,file='receiver.dat',form='formatted',status='replace')
  write(12,*) NRC
  do ii=1,NC
    do ix=RXS(ii),RXE(ii),RXD(ii)
      do iy=RYS(ii),RYE(ii),RYD(ii)
        do iz=RZS(ii),RZE(ii),RZD(ii)
          write(12,*) ix,iy,iz
        end do
      end do
    end do
  end do
  close(unit=12)

  deallocate(RXS,RXE,RXD,RYS,RYE,RYD,RZS,RZE,RZD)  
end program buildreceiver