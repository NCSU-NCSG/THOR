!THOR is licensed under the MIT License.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!> @brief Program driver.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
PROGRAM ahot_c_ug
  !***********************************************************************
  !
  ! Arbitrarily High Order Transport-
  ! Characteristic Type
  ! in Unstructered Tetrahedral Grids
  !
  ! Author: R.M. Ferrer, Y.Y. Azmy, S. Schunert, R.A. Yessayan, and N.F. Herring
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
  USE globals

  ! Use modules that contain necessary subroutines and functions to
  ! execute transport code

  USE error_module
  USE setup_module
  USE execution_module
  USE wrapup_module
  USE adjoint_module

  IMPLICIT NONE

  ! Declare scalar flux

  REAL(kind=d_t),ALLOCATABLE :: flux(:,:,:,:,:)

  ! Declare the eigenvalue

  REAL(kind=d_t) :: keffective

  ! Declare temporary variables

  INTEGER(kind=li) :: alloc_stat

  !######### MPI variables
  INTEGER :: mpi_err, num_p
  !#########

  ! timing variables
  INTEGER:: do_timing = 0,imain
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

  !set log filename
  CALL GET_COMMAND_ARGUMENT(1,jobname)
  !find extension start
  DO imain=LEN_TRIM(jobname),1,-1
    IF(jobname(imain:imain) .EQ. '.')EXIT
  ENDDO
  !if it has an extension, check if it's an input extension and cut it from the logname
  temp=TRIM(jobname)
  IF(imain .GE. 2)THEN
    temp=jobname(imain:LEN_TRIM(jobname))
    SELECTCASE(TRIM(temp))
      CASE('.i','in','.inp')
        temp=jobname(1:imain-1)
      CASE DEFAULT
        temp=TRIM(jobname)
    ENDSELECT
  ENDIF
  jobname=TRIM(temp)
  OPEN(unit = log_unit, file = TRIM(ADJUSTL(jobname))//'.log', status = "REPLACE", action = "WRITE")

  IF (rank .EQ. 0) THEN

    ! Print banner for THOR

    CALL printlog("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    CALL printlog('')
    CALL printlog("   TTTTTTT  HH     HH  OOOOO  RRRRRR ")
    CALL printlog("     TTT    HH     HH OOOOOOO RRRRRRR")
    CALL printlog("     TTT    HH     HH OO   OO RR   RR")
    CALL printlog("     TTT    HHHHHHHHH OO   OO RRRRRR ")
    CALL printlog("     TTT    HHHHHHHHH OO   OO RRRR   ")
    CALL printlog("     TTT    HH     HH OO   OO RR RR  ")
    CALL printlog("     TTT    HH     HH 0000000 RR  RR ")
    CALL printlog("     TTT    HH     HH  OOOOO  RR   RR")
    CALL printlog('')
    CALL printlog("   Tetrahedral High Order Radiation Transport Code")
    CALL printlog('')
    CALL printlog("   By R.M. Ferrer, Y.Y. Azmy, S. Schunert, R.A. Yessayan, and N.F. Herring")
    CALL printlog('')
    CALL printlog("   Version 1.0.0 - Update 09/07/2022")
    CALL printlog("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    CALL printlog('')
  END IF

  !***********************************************************************
  ! Setup module reads input and mesh file, prepares input for execution
  !***********************************************************************
  IF (do_timing .EQ. 1) parallel_timing(2,1) = MPI_WTIME()
  CALL setup
  IF (do_timing .EQ. 1) parallel_timing(2,2) = MPI_WTIME()
  IF(num_p .GT. quad_ord*4)THEN
    CALL raise_warning(TRIM(str(num_p))//" MPI processes requested. But the angular quadrature &
        &order is "//TRIM(str(quad_ord))//".")
    CALL raise_warning("THOR performance is poor for more than 4*n (where n is the quadrature &
        &order) processes.")
    CALL raise_warning(TRIM(str(quad_ord*4))//" or less MPI processes are recommended for this &
        &input.")
  ENDIF
  ! If execution is not desired then stop here
  IF(execution .EQ. 0)THEN
    call sleep(2)
    CALL printlog("*****************************************************************************")
    CALL printlog("*****************************************************************************")
    CALL printlog("*****************************************************************************")
    CALL printlog("User specified no execution. Finalizing THOR and stopping.")
    CALL printlog("*****************************************************************************")
    CALL MPI_FINALIZE(mpi_err)
    CALL thor_success
  ENDIF

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
    INQUIRE(file = converge_filename, exist = existence)
    OPEN(unit = 21, file = converge_filename, status = "REPLACE", action = "WRITE")
  END IF




  !***********************************************************************
  ! Call execution to perform the actual computation
  !***********************************************************************
  IF(adjoint_opt)THEN
    IF(eig_switch .EQ. 1)STOP 'jfnk not supported for adjoint solutions'
    CALL transpose_xs()
    most_thermal=1
  ENDIF
  IF (do_timing .EQ. 1) parallel_timing(3,1) = MPI_WTIME()
  IF (problem==0) THEN
    CALL execute_ext(flux)
  ELSE
    CALL execute_eig(flux,keffective)
  END IF
  IF (do_timing .EQ. 1) parallel_timing(3,2) = MPI_WTIME()
  IF(adjoint_opt)THEN
    CALL transpose_xs()
    CALL reverse_odd_mom(flux)
  ENDIF
  !***********************************************************************
  ! Call wrapup to finish up post-processing and output results
  !***********************************************************************


  IF(print_conv.EQ.1 .AND. rank .EQ. 0) THEN
    CLOSE(unit=21)
  END IF

  IF (do_timing .EQ. 1) parallel_timing(4,1) = MPI_WTIME()
  CALL wrapup(flux = flux,keff = keffective, unit_number = stdout_unit, suffix = "", is_final = .TRUE.)
  CALL wrapup(flux = flux,keff = keffective, unit_number = log_unit, suffix = "", is_final = .TRUE.)
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
  USE globals
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
      CALL printlog('')
      CALL printlog('/=====================================================')
      WRITE(amsg,'(A,I0)')'| Timing Data for process: ', i
      CALL printlog(amsg)
      CALL printlog('|-----------------------------------------------------')
      WRITE(amsg,'(A,ES12.4)')'| Total Measured Time:      ', print_timing(1,2)- print_timing(1,1)
      CALL printlog(amsg)
      WRITE(amsg,'(A,ES12.4)')'| Total Setup Time:         ', print_timing(2,2)- print_timing(2,1)
      CALL printlog(amsg)
      WRITE(amsg,'(A,ES12.4)')'| Total Execute Time:       ', print_timing(3,2)- print_timing(3,1)
      CALL printlog(amsg)
      WRITE(amsg,'(A,ES12.4)')'| Total Wrapup Time:        ', print_timing(4,2)- print_timing(4,1)
      CALL printlog(amsg)
      WRITE(amsg,'(A,ES12.4)')'| Total Non-Accounted Time: ',(print_timing(1,2)- print_timing(1,1))- &
                                               (print_timing(2,2)- print_timing(2,1)) - &
                                               (print_timing(3,2)- print_timing(3,1)) - &
                                               (print_timing(4,2)- print_timing(4,1))
      CALL printlog(amsg)
      CALL printlog('\=====================================================')
    END DO
  END IF

END SUBROUTINE
