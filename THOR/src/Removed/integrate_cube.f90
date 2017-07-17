PROGRAM integrate_cube
  !***********************************************************************
  !
  ! Integrates AHOT-C fluxes into 27 'canonical' 1 mfp x 1 mfp x 1 mfp
  ! cubes
  !
  ! Author: R.M. Ferrer
  !
  !***********************************************************************

  ! Use user derived types

  USE types

  IMPLICIT NONE

  ! Declare variables
  INTEGER(kind=li) :: n, i, j, k, ii, jj, kk, indexii, indexi, tet
  INTEGER(kind=li) :: num_elements, alloc_stat
  REAL(kind=d_t) :: reference_flux, subcube_volume
  REAL(kind=d_t), DIMENSION(3,3,3) :: ref_flux
  REAL(kind=d_t), DIMENSION(:), ALLOCATABLE :: volume, fine_flux

  ! Open problem AHOT-C output and read

  OPEN(unit=10,file='flux_file',status='unknown',action='read')

  ! Read flux file

  READ(10,*) num_elements

  n=INT(EXP(LOG(num_elements/135.0)/3.0))

  ALLOCATE(volume(num_elements),fine_flux(num_elements),&
        stat=alloc_stat)
  IF(alloc_stat /= 0) STOP "*** not enough memory ***"

  DO i=1, num_elements
    READ(10,*) volume(i), fine_flux(i)
  END DO

  indexii=0
  indexi=0

  DO kk=1, 3
    DO jj=1, 3
      DO ii=1, 3
        indexii=5*n*(ii-1)+5*3*n*n*(jj-1)+5*9*n*n*n*(kk-1)
        reference_flux=0.0
        subcube_volume=0.0
        DO k=1, n
          DO j=1, n
            DO i=1, n
              DO tet=1, 5
                indexi=tet+indexii+5*(i-1)+5*3*n*(j-1)+&
                      5*9*n*n*(k-1)
                reference_flux=reference_flux+fine_flux(indexi)*&
                      volume(indexi)
                subcube_volume=subcube_volume+volume(indexi)
              END DO
            END DO
          END DO
        END DO
        ref_flux(ii,jj,kk)=reference_flux/subcube_volume
      END DO
    END DO
  END DO

  CLOSE(10)

  ! Write reference solution into output

  OPEN(unit=20,file='ahot_flux',status='unknown',action='write')

  DO kk=1, 3
    DO jj=1, 3
      DO ii=1, 3
        WRITE(20,*) ii, jj, kk, ref_flux(ii,jj,kk)
      END DO
    END DO
  END DO

  CLOSE(20)

END PROGRAM integrate_cube
