!THOR is licensed under the MIT License.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!> @brief Execution module.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
MODULE execution_module
  !***********************************************************************
  !
  ! Execution module contains all subroutines to run problem
  !
  !***********************************************************************

  ! Use derived-type modules

  USE types
  USE parameter_types
  USE filename_types
  USE vector_types
  USE cross_section_types
  USE geometry_types
  USE angle_types
  USE multindex_types
  USE globals

  ! Use modules that pertain setting up problem

  USE solver_module

  IMPLICIT NONE

CONTAINS

  !> Is a simple calling subroutine that starts the timer and calls the fixed source solver
  SUBROUTINE execute_ext(flux)
    !**********************************************************************
    !
    ! Subroutine execution calls routines that performs computation
    !
    !**********************************************************************

    ! Declare scalar flux

    REAL(kind=d_t) :: flux(num_moments_v,namom,num_cells,egmax,niter)

    ! Call solver to initiate problem execution

    CALL CPU_TIME(start)

    CALL solver_ext(flux)

    CALL CPU_TIME(finish)

  END SUBROUTINE execute_ext

  !-----------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------

  !> Is a simple calling subroutine that starts the timer and calls the eigenvalue solver
  SUBROUTINE execute_eig(flux,keff)
    !**********************************************************************
    !
    ! Subroutine execution calls routines that performs computation
    !
    !**********************************************************************

    ! Declare k-effective and scalar flux

    REAL(kind=d_t), INTENT(inout) :: keff
    REAL(kind=d_t) :: flux(num_moments_v,namom,num_cells,egmax,niter)

    ! Call solver to initiate problem execution

    CALL CPU_TIME(start)

    CALL solver_eig(flux,keff)

    CALL CPU_TIME(finish)

  END SUBROUTINE execute_eig

  !-----------------------------------------------------------------------------------------
  ! End
  !-----------------------------------------------------------------------------------------

END MODULE execution_module
