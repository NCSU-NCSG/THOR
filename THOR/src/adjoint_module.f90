MODULE adjoint_module
  !***********************************************************************
  ! This module contains subroutines for the adjoint solver.
  !***********************************************************************
  USE globals
  USE mpi
  IMPLICIT NONE
  PRIVATE
  !
  PUBLIC :: transpose_xs

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
    !stop 'transpose_xs not complete'
  ENDSUBROUTINE transpose_xs
END MODULE adjoint_module
