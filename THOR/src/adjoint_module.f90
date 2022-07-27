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

  SUBROUTINE transpose_xs
    !**********************************************************************
    !
    ! Transposes XS matrix by energy for the adjoint solver
    !
    !**********************************************************************
    stop 'transpose_xs not complete'
  ENDSUBROUTINE transpose_xs
END MODULE adjoint_module
