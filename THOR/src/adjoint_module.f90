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
    INTEGER :: m,g
    REAL(8),ALLOCATABLE :: temp_mat(:,:)

    ALLOCATE(temp_mat(egmax,egmax))
    !loop over all materials
    DO m=1,num_mat
      DO g=1,egmax
      ENDDO
    ENDDO
    DEALLOCATE(temp_mat)
    stop 'transpose_xs not complete'
  ENDSUBROUTINE transpose_xs
END MODULE adjoint_module
