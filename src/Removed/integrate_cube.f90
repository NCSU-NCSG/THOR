program integrate_cube
!***********************************************************************
!
! Integrates AHOT-C fluxes into 27 'canonical' 1 mfp x 1 mfp x 1 mfp 
! cubes
!
! Author: R.M. Ferrer
!
!***********************************************************************

! Use user derived types

  use types

  implicit none

! Declare variables
  integer(kind=li) :: n, i, j, k, ii, jj, kk, indexii, indexi, tet
  integer(kind=li) :: num_elements, alloc_stat
  real(kind=d_t) :: reference_flux, subcube_volume
  real(kind=d_t), dimension(3,3,3) :: ref_flux
  real(kind=d_t), dimension(:), allocatable :: volume, fine_flux

! Open problem AHOT-C output and read

  open(unit=10,file='flux_file',status='unknown',action='read')

! Read flux file

  read(10,*) num_elements

  n=INT(EXP(LOG(num_elements/135.0)/3.0))

  allocate(volume(num_elements),fine_flux(num_elements),&
  stat=alloc_stat)
  if(alloc_stat /= 0) stop "*** not enough memory ***" 

  do i=1, num_elements
     read(10,*) volume(i), fine_flux(i)
  end do

  indexii=0
  indexi=0

  do kk=1, 3
     do jj=1, 3
        do ii=1, 3
           indexii=5*n*(ii-1)+5*3*n*n*(jj-1)+5*9*n*n*n*(kk-1)
           reference_flux=0.0
           subcube_volume=0.0
           do k=1, n
              do j=1, n
                 do i=1, n
                    do tet=1, 5
                       indexi=tet+indexii+5*(i-1)+5*3*n*(j-1)+&
                            5*9*n*n*(k-1)
                       reference_flux=reference_flux+fine_flux(indexi)*&
                            volume(indexi)
                       subcube_volume=subcube_volume+volume(indexi)
                    end do
                 end do
              end do
           end do
        ref_flux(ii,jj,kk)=reference_flux/subcube_volume
        end do
     end do
  end do

  close(10)

! Write reference solution into output

  open(unit=20,file='ahot_flux',status='unknown',action='write')

  do kk=1, 3
     do jj=1, 3
        do ii=1, 3
           write(20,*) ii, jj, kk, ref_flux(ii,jj,kk)
        end do
     end do
  end do

  close(20)

end program integrate_cube
