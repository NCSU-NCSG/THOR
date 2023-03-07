!THOR is licensed under the MIT License.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!> @brief Adjoint module.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
MODULE adjoint_module
  !***********************************************************************
  ! This module contains subroutines for the adjoint solver.
  !***********************************************************************
  USE globals
  USE mpi
  IMPLICIT NONE
  PRIVATE
  !
  PUBLIC :: transpose_xs, reverse_odd_mom

CONTAINS

  !------------------------------------------------------------------------------------------------------------!
  !------------------------------------------------------------------------------------------------------------!

  !**********************************************************************
  !
  ! Transposes XS production matrix by energy for the adjoint solver
  !
  !**********************************************************************
  SUBROUTINE transpose_xs
    INTEGER :: m,g,gp,j
    REAL(8),ALLOCATABLE :: temp_mat(:,:)

    ALLOCATE(temp_mat(egmax,egmax))
    !loop over all materials
    DO m=1,num_mat
      !loop over angular moments for scattering
      DO j=1,xs_ord+1
        temp_mat=xs_mat(m)%sigma_scat(j,:,:)
        !transpose the matrix
        DO g=1,egmax
          DO gp=1,egmax
            xs_mat(m)%sigma_scat(j,g,gp)=temp_mat(gp,g)
          ENDDO
        ENDDO
      ENDDO
      !swap nusigmaf and chi
      temp_mat=0
      temp_mat(1,:)=xs_mat(m)%chi(:)
      temp_mat(2,:)=xs_mat(m)%nusig_f(:)
      xs_mat(m)%chi(:)=temp_mat(2,:)
      xs_mat(m)%nusig_f(:)=temp_mat(1,:)
    ENDDO
    DEALLOCATE(temp_mat)
  ENDSUBROUTINE transpose_xs

  !**********************************************************************
  !
  ! Reverses odd flux moments
  !
  !**********************************************************************
  SUBROUTINE reverse_odd_mom(flux)
    REAL(8),INTENT(INOUT) :: flux(num_moments_v,namom,num_cells,egmax,niter)
    INTEGER :: i

    !reverse odd flux moments (they are indexed with the even indeces since it's indexed from 1)
    DO i=2,num_moments_v,2
      flux(i,:,:,:,:)=-flux(i,:,:,:,:)
    ENDDO
  ENDSUBROUTINE reverse_odd_mom
END MODULE adjoint_module
