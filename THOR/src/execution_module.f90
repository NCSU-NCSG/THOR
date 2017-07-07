module execution_module
!***********************************************************************
!
! Execution module contains all subroutines to run problem
!
!***********************************************************************

! Use derived-type modules

  use types
  use parameter_types
  use filename_types
  use vector_types
  use cross_section_types
  use geometry_types
  use angle_types
  use multindex_types
  use global_variables

! Use modules that pertain setting up problem

  use solver_module  

  implicit none

contains

  !> Is a simple calling subroutine that starts the timer and calls the fixed source solver
  subroutine execute_ext(flux)
  !**********************************************************************
  !
  ! Subroutine execution calls routines that performs computation
  !
  !**********************************************************************

  ! Declare scalar flux 

    real(kind=d_t) :: flux(num_moments_v,namom,num_cells,egmax,niter) 

  ! Call solver to initiate problem execution

    call cpu_time(start)    

    call solver_ext(flux)

    call cpu_time(finish)

  end subroutine execute_ext

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

  !> Is a simple calling subroutine that starts the timer and calls the eigenvalue solver
  subroutine execute_eig(flux,keff)
  !**********************************************************************
  !
  ! Subroutine execution calls routines that performs computation
  !
  !**********************************************************************

  ! Declare k-effective and scalar flux

    real(kind=d_t), intent(inout) :: keff
    real(kind=d_t) :: flux(num_moments_v,namom,num_cells,egmax,niter) 
    
  ! Call solver to initiate problem execution

    call cpu_time(start)    

    call solver_eig(flux,keff)

    call cpu_time(finish)

  end subroutine execute_eig

!-----------------------------------------------------------------------------------------
! End
!-----------------------------------------------------------------------------------------  
  
end module execution_module

