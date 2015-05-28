program ahot_c_ug
!***********************************************************************
!
! Arbitrarily High Order Transport-
! Characteristic Type
! in Unstructered Tetrahedral Grids
!
! Author: R.M. Ferrer
!
!
!***********************************************************************

! Use derived-type modules

!FIX ME - Remove all non-communication MPI calls from all files and abstract rank & process count to global
!Fix ME - Remove timing and assign work (in setup) subroutine and move to utils.mod or such

!#########
  use mpi
!#########  
  use types
  use parameter_types
  use filename_types
  use vector_types
  use cross_section_types
  use geometry_types
  use angle_types
  use multindex_types
  use global_variables

! Use modules that contain necessary subroutines and functions to 
! execute transport code

  use termination_module
  use setup_module
  use execution_module
  use wrapup_module

  implicit none

! Declare scalar flux

  real(kind=d_t),allocatable :: flux(:,:,:,:,:)

! Declare the eigenvalue

  real(kind=d_t) :: keffective

! Declare temporary variables

  integer(kind=li) :: alloc_stat, eg, i, ii, l, it, n, m
  
  !######### MPI variables 
  integer :: mpi_err , mpi_sig_size(MPI_STATUS_SIZE), num_p
  !#########  
  
  ! timing variables
  integer:: do_timing = 0
  
  logical          :: existence
  character(100):: temp

  
  call GET_COMMAND_ARGUMENT(2,temp)
  if (trim(temp) .eq. '-t') do_timing = 1
  
  !#########  
  call MPI_INIT(mpi_err)
  
  if (do_timing .eq. 1) parallel_timing(1,1) = MPI_WTIME()
  
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, num_p, mpi_err)
  !#########  
  
  if (rank .EQ. 0) then
    
  ! Print banner for THOR

    write(6,*) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    write(6,*) 
    write(6,*) "   TTTTTTT  HH     HH  OOOOO  RRRRRR "
    write(6,*) "     TTT    HH     HH OOOOOOO RRRRRRR"
    write(6,*) "     TTT    HH     HH OO   OO RR   RR"
    write(6,*) "     TTT    HHHHHHHHH OO   OO RRRRRR "
    write(6,*) "     TTT    HHHHHHHHH OO   OO RRRR   "
    write(6,*) "     TTT    HH     HH OO   OO RR RR  "
    write(6,*) "     TTT    HH     HH 0000000 RR  RR "
    write(6,*) "     TTT    HH     HH  OOOOO  RR   RR"
    write(6,*) 
    write(6,*) "   Tetrahedral High Order Radiation Transport Code"
    write(6,*)
    write(6,*) "   By R. M. Ferrer"
    write(6,*) 
    write(6,*) "   Version 1.0 BETA - Update 05/10/2012"
    write(6,*) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    write(6,*)
  end if

!***********************************************************************
! Setup module reads input and mesh file, prepares input for execution
!***********************************************************************
  if (do_timing .eq. 1) parallel_timing(2,1) = MPI_WTIME()
  call setup
  if (do_timing .eq. 1) parallel_timing(2,2) = MPI_WTIME()
!***********************************************************************
! Allocate scalar flux 
!***********************************************************************

  allocate(flux(num_moments_v,namom,num_cells,egmax,niter),stat=alloc_stat)
  if(alloc_stat /= 0) call stop_thor(2_li)
  flux = zero

  
!***********************************************************************
! Do some preliminary stuff
!***********************************************************************

  if(print_conv.eq.1 .and. rank .eq. 0) then
    inquire(file = "thor.convergence", exist = existence)  
    if(existence) then
      open (unit = 21, file = "thor.convergence", status = "OLD", action = "WRITE")
    else
      open (unit = 21, file = "thor.convergence", status = "NEW", action = "WRITE")
    end if
    write(21,*) '========================================================'
    write(21,*) '   Begin outer iterations.'
    write(21,*) '========================================================'
  end if


  
  
!***********************************************************************
! Call execution to perform the actual computation
!***********************************************************************
  if (do_timing .eq. 1) parallel_timing(3,1) = MPI_WTIME()
  if (problem==0) then
    call execute_ext(flux)
  else
    call execute_eig(flux,keffective)
  end if
  if (do_timing .eq. 1) parallel_timing(3,2) = MPI_WTIME()
!***********************************************************************
! Call wrapup to finish up post-processing and output results
!***********************************************************************
  
  if(print_conv.eq.1) close(unit=21)
  
  if (do_timing .eq. 1) parallel_timing(4,1) = MPI_WTIME()
  call wrapup(flux = flux,keff = keffective, unit_number = 6, suffix = "")
  if (do_timing .eq. 1) parallel_timing(4,2) = MPI_WTIME()
!***********************************************************************
! Cleanup
!***********************************************************************

  ! Flux

  deallocate(flux)

!***********************************************************************
! Terminate execution
!***********************************************************************
  
  if (do_timing .eq. 1) parallel_timing(1,2) = MPI_WTIME()
  call MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
  
  call write_timing
  call MPI_FINALIZE(mpi_err)
  call stop_thor(1_li)
  
end program ahot_c_ug

subroutine write_timing

  use mpi
  use global_variables
  implicit none

  integer :: i, num_p, mpi_err, err_size(MPI_STATUS_SIZE), do_timing=0
  real*8:: print_timing(4,2)
	character(100):: temp  
  call GET_COMMAND_ARGUMENT(2,temp)
  if (trim(temp) .eq. '-t') do_timing = 1
  flush(6)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, num_p, mpi_err)
  call MPI_SEND(parallel_timing, 8, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD, mpi_err)
  
  if (rank .eq. 0 .and. do_timing .eq. 1) then 
    do i =0, num_p-1
      call MPI_RECV(print_timing, 8, MPI_DOUBLE, i, i, MPI_COMM_WORLD, err_size, mpi_err )
      write(*,*)
      write(*,*)'/====================================================='
      write(*,*)'| Timing Data for process: ', i
      write(*,*)'|-----------------------------------------------------'
      write(*,*)'| Total Measured Time:      ', print_timing(1,2)- print_timing(1,1)
      write(*,*)'| Total Setup Time:         ', print_timing(2,2)- print_timing(2,1)
      write(*,*)'| Total Execute Time:       ', print_timing(3,2)- print_timing(3,1)
      write(*,*)'| Total Wrapup Time:        ', print_timing(4,2)- print_timing(4,1)
      write(*,*)'| Total Non-Accounted Time: ',(print_timing(1,2)- print_timing(1,1))- &
                                               (print_timing(2,2)- print_timing(2,1)) - &
                                               (print_timing(3,2)- print_timing(3,1)) - &
                                               (print_timing(4,2)- print_timing(4,1))
      write(*,*)'\====================================================='
    end do
  end if

end subroutine


