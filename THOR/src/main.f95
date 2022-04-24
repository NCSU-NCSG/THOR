PROGRAM ahot_c_ug
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
  USE mpi
  USE ISO_FORTRAN_ENV
  !#########
  USE types
  USE parameter_types
  USE filename_types
  USE vector_types
  USE cross_section_types
  USE geometry_types
  USE angle_types
  USE multindex_types
  USE global_variables

  ! Use modules that contain necessary subroutines and functions to
  ! execute transport code

  USE error_module
  USE setup_module
  USE execution_module
  USE wrapup_module

  IMPLICIT NONE

  ! Declare scalar flux

  REAL(kind=d_t),ALLOCATABLE :: flux(:,:,:,:,:)

  ! Declare the eigenvalue

  REAL(kind=d_t) :: keffective

  ! Declare temporary variables

  INTEGER(kind=li) :: alloc_stat,converge_unit

  !######### MPI variables
  INTEGER :: mpi_err, num_p
  !#########

  ! timing variables
  INTEGER:: do_timing = 0
  LOGICAL        :: existence
  CHARACTER(100) :: temp

  stdout_unit=OUTPUT_UNIT

  CALL GET_COMMAND_ARGUMENT(2,temp)
  IF (TRIM(temp) .EQ. '-t') do_timing = 1

  !#########
  CALL MPI_INIT(mpi_err)

  IF (do_timing .EQ. 1) parallel_timing(1,1) = MPI_WTIME()

  CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, num_p, mpi_err)
  !#########

  IF (rank .EQ. 0) THEN

    ! Print banner for THOR

    WRITE(stdout_unit,*) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    WRITE(stdout_unit,*)
    WRITE(stdout_unit,*) "   TTTTTTT  HH     HH  OOOOO  RRRRRR "
    WRITE(stdout_unit,*) "     TTT    HH     HH OOOOOOO RRRRRRR"
    WRITE(stdout_unit,*) "     TTT    HH     HH OO   OO RR   RR"
    WRITE(stdout_unit,*) "     TTT    HHHHHHHHH OO   OO RRRRRR "
    WRITE(stdout_unit,*) "     TTT    HHHHHHHHH OO   OO RRRR   "
    WRITE(stdout_unit,*) "     TTT    HH     HH OO   OO RR RR  "
    WRITE(stdout_unit,*) "     TTT    HH     HH 0000000 RR  RR "
    WRITE(stdout_unit,*) "     TTT    HH     HH  OOOOO  RR   RR"
    WRITE(stdout_unit,*)
    WRITE(stdout_unit,*) "   Tetrahedral High Order Radiation Transport Code"
    WRITE(stdout_unit,*)
    WRITE(stdout_unit,*) "   By R. M. Ferrer"
    WRITE(stdout_unit,*)
    WRITE(stdout_unit,*) "   Version 1.0 BETA - Update 05/10/2012"
    WRITE(stdout_unit,*) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    WRITE(stdout_unit,*)
  END IF

  !***********************************************************************
  ! Setup module reads input and mesh file, prepares input for execution
  !***********************************************************************
  IF (do_timing .EQ. 1) parallel_timing(2,1) = MPI_WTIME()
  CALL setup
  IF (do_timing .EQ. 1) parallel_timing(2,2) = MPI_WTIME()
  !***********************************************************************
  ! Allocate scalar flux
  !***********************************************************************

  ALLOCATE(flux(num_moments_v,namom,num_cells,egmax,niter),stat=alloc_stat)
  IF(alloc_stat /= 0) CALL raise_fatal_error("*** Not enough memory ***")
  flux = zero


  !***********************************************************************
  ! Do some preliminary stuff
  !***********************************************************************

  IF(print_conv.EQ.1 .AND. rank .EQ. 0) THEN
    INQUIRE(file = "thor.convergence", exist = existence)
    converge_unit=21
    IF(existence) THEN
      OPEN(unit = converge_unit, file = "thor.convergence", status = "OLD", action = "WRITE")
    ELSE
      OPEN(unit = converge_unit, file = "thor.convergence", status = "NEW", action = "WRITE")
    END IF
    WRITE(converge_unit,*) '========================================================'
    WRITE(converge_unit,*) '   Begin outer iterations.'
    WRITE(converge_unit,*) '========================================================'
  END IF




  !***********************************************************************
  ! Call execution to perform the actual computation
  !***********************************************************************
  IF (do_timing .EQ. 1) parallel_timing(3,1) = MPI_WTIME()
  IF (problem==0) THEN
    CALL execute_ext(flux)
  ELSE
    CALL execute_eig(flux,keffective)
  END IF
  IF (do_timing .EQ. 1) parallel_timing(3,2) = MPI_WTIME()
  !***********************************************************************
  ! Call wrapup to finish up post-processing and output results
  !***********************************************************************

  IF(print_conv.EQ.1) CLOSE(unit=21)

  IF (do_timing .EQ. 1) parallel_timing(4,1) = MPI_WTIME()
  CALL wrapup(flux = flux,keff = keffective, unit_number = 6, suffix = "", is_final = .TRUE.)
  IF (do_timing .EQ. 1) parallel_timing(4,2) = MPI_WTIME()
  !***********************************************************************
  ! Cleanup
  !***********************************************************************

  ! Flux

  DEALLOCATE(flux)

  !***********************************************************************
  ! Terminate execution
  !***********************************************************************

  IF (do_timing .EQ. 1) parallel_timing(1,2) = MPI_WTIME()
  CALL MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
  IF (do_timing .EQ. 1) CALL write_timing
  CALL MPI_FINALIZE(mpi_err)
  CALL thor_success

END PROGRAM ahot_c_ug

SUBROUTINE write_timing

  USE mpi
  USE global_variables
  IMPLICIT NONE

  INTEGER :: i, num_p, mpi_err, err_size(MPI_STATUS_SIZE), do_timing=0
  REAL*8:: print_timing(4,2)
  CHARACTER(100) :: temp
  CALL GET_COMMAND_ARGUMENT(2,temp)
  IF (TRIM(temp) .EQ. '-t') do_timing = 1
  flush(6)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, num_p, mpi_err)
  CALL MPI_SEND(parallel_timing, 8, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD, mpi_err)

  IF (rank .EQ. 0 .AND. do_timing .EQ. 1) THEN
    DO i =0, num_p-1
      CALL MPI_RECV(print_timing, 8, MPI_DOUBLE, i, i, MPI_COMM_WORLD, err_size, mpi_err )
      WRITE(stdout_unit,*)
      WRITE(stdout_unit,*)'/====================================================='
      WRITE(stdout_unit,*)'| Timing Data for process: ', i
      WRITE(stdout_unit,*)'|-----------------------------------------------------'
      WRITE(stdout_unit,*)'| Total Measured Time:      ', print_timing(1,2)- print_timing(1,1)
      WRITE(stdout_unit,*)'| Total Setup Time:         ', print_timing(2,2)- print_timing(2,1)
      WRITE(stdout_unit,*)'| Total Execute Time:       ', print_timing(3,2)- print_timing(3,1)
      WRITE(stdout_unit,*)'| Total Wrapup Time:        ', print_timing(4,2)- print_timing(4,1)
      WRITE(stdout_unit,*)'| Total Non-Accounted Time: ',(print_timing(1,2)- print_timing(1,1))- &
                                               (print_timing(2,2)- print_timing(2,1)) - &
                                               (print_timing(3,2)- print_timing(3,1)) - &
                                               (print_timing(4,2)- print_timing(4,1))
      WRITE(stdout_unit,*)'\====================================================='
    END DO
  END IF

END SUBROUTINE
