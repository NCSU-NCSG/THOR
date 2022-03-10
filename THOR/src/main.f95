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

  USE termination_module
  USE setup_module
  USE execution_module
  USE wrapup_module
  USE ITMM_module

  IMPLICIT NONE

  ! Declare scalar flux

  REAL(kind=d_t),ALLOCATABLE :: flux(:,:,:,:,:)

  ! Declare the eigenvalue

  REAL(kind=d_t) :: keffective, IOtime_start

  ! Declare temporary variables

  INTEGER(kind=li) :: alloc_stat, eg, i, ii, l, it, n, m,j

  !######### MPI variables
  INTEGER :: mpi_err , mpi_sig_size(MPI_STATUS_SIZE)
  !#########

  INTEGER :: global_maxreg, global_minreg

  ! timing variables
  INTEGER:: do_timing = 0
  LOGICAL        :: existence
  CHARACTER(250) :: temp

  CALL GET_COMMAND_ARGUMENT(2,temp)
  IF (TRIM(temp) .EQ. '-t') do_timing = 1

  !#########
  CALL MPI_INIT(mpi_err)

  !! ADD REMOVED - OCT 2019
  !IF (do_timing .EQ. 1) parallel_timing(1,1) = MPI_WTIME()
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, PBJrank, mpi_err)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, num_proc, mpi_err)
  !! ADD REMOVED - OCT 2019
  !#########

  IF (PBJrank .EQ. 0) THEN

    ! Print banner for THOR

    WRITE(6,*) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    WRITE(6,*)
    WRITE(6,*) "   TTTTTTT HH     HH  OOOOO  RRRRRR "
    WRITE(6,*) "     TTT   HH     HH OOOOOOO RRRRRRR"
    WRITE(6,*) "     TTT   HH     HH OO   OO RR   RR"
    WRITE(6,*) "     TTT   HHHHHHHHH OO   OO RRRRRR "
    WRITE(6,*) "     TTT   HHHHHHHHH OO   OO RRRR   "
    WRITE(6,*) "     TTT   HH     HH OO   OO RR RR  "
    WRITE(6,*) "     TTT   HH     HH OOOOOOO RR  RR "
    WRITE(6,*) "     TTT   HH     HH  OOOOO  RR   RR"
    WRITE(6,*)
    WRITE(6,*) "   Tetrahedral High Order Radiation Transport Code"
    WRITE(6,*)
    WRITE(6,*) "   By R. M. Ferrer"
    WRITE(6,*)
    WRITE(6,*) "   Version 1.0 BETA - Update 05/10/2012"
    WRITE(6,*) "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    WRITE(6,*)
    WRITE(6,*)'System Clock Precision (s) : ',MPI_WTICK()
  END IF


  !***********************************************************************
  ! Setup module reads input and mesh file, prepares input for execution
  !***********************************************************************
  !! ADD REMOVED - OCT 2019
  !IF (do_timing .EQ. 1) parallel_timing(2,1) = MPI_WTIME()
  !! ADD REMOVED - OCT 2019
  IOtime_start=MPI_WTIME()
  CALL setup
  IOtime=IOtime+MPI_WTIME()-IOtime_start
  !! ADD REMOVED - OCT 2019
  !IF (do_timing .EQ. 1) parallel_timing(2,2) = MPI_WTIME()
  !! ADD REMOVED - OCT 2019
  !***********************************************************************
  ! Allocate scalar flux
  !***********************************************************************

  ALLOCATE(flux(num_moments_v,namom,num_cells,egmax,niter),stat=alloc_stat)
  IF(alloc_stat /= 0) CALL stop_thor(2_li)
  flux = zero

  !***********************************************************************
  ! Allocate ITMM matriindinces it we're performing construction
  !***********************************************************************

  IF (ITMM.EQ.2) THEN
  	ALLOCATE(Jphi(num_cells,num_cells,egmax), &
  	Jpsi(4*nangle*N_side_SDbound,num_cells,egmax), &
  	Kphi(num_cells,4*nangle*N_side_SDbound,egmax))
  	Jphi=zero
  	Jpsi=zero
  	Kphi=zero
  	Kpsi_reallocate=4*nangle*N_side_SDbound
  END IF
  ALLOCATE(ITMMKindex(4*nangle*N_side_SDbound,4,2))
  ITMMKindex=-99
  ALLOCATE(psiin(4*nangle*N_side_SDbound,egmax),psiout(4*nangle*N_side_SDbound,egmax))
  psiin=0.0d0



  !***********************************************************************
  ! Do some preliminary stuff
  !***********************************************************************

  IF(print_conv.EQ.1 .AND. PBJrank .EQ. 0) THEN
    INQUIRE(file = "thor.convergence", exist = existence)
    IF(existence) THEN
      OPEN (unit = 21, file = "thor.convergence", status = "OLD", action = "WRITE")
    ELSE
      OPEN (unit = 21, file = "thor.convergence", status = "NEW", action = "WRITE")
    END IF
    WRITE(21,*) '========================================================'
    WRITE(21,*) '   Begin outer iterations.'
    WRITE(21,*) '========================================================'
  END IF

89 CONTINUE



   !If we just finished GFIC (or IPBJ setup), we need to redo the readmesh routine,
   !to reset BCs that were artificially set to vacuum

   IF ((ITMM.EQ.1).OR.(ITMM.EQ.0)) THEN

    !Store global minreg and maxreg, since they will get messed up in readmesh
    global_maxreg=maxreg
    global_minreg=minreg

    IOtime_start=MPI_WTIME()
    !Reread mesh
    CALL read_tetmesh

    maxreg=global_maxreg
    minreg=global_minreg

    !Redo refl BC
    IF (rside_cells.GT.0) THEN
        CALL check_reflective
        CALL classify_reflective
    END IF
    IOtime=IOtime+MPI_WTIME()-IOtime_start

    !Initial communcation of ITMM unpack instructions
    comm_instructions_start_time = MPI_WTIME()
    CALL ITMM_communicate_instructions
    comm_instructions_end_time = MPI_WTIME()

   END IF

  !***********************************************************************
  ! Call execution to perform the actual computation
  !***********************************************************************
  !! ADD REMOVED - OCT 2019
  !IF (do_timing .EQ. 1) parallel_timing(3,1) = MPI_WTIME()
  !! ADD REMOVED - OCT 2019
  IF (problem==0.OR.ITMM.EQ.2.OR.ITMM.EQ.3) THEN
    CALL execute_ext(flux)
  ELSE
    CALL execute_eig(flux,keffective)
  END IF
  !! ADD REMOVED - OCT 2019
  !IF (do_timing .EQ. 1) parallel_timing(3,2) = MPI_WTIME()
  !! ADD REMOVED - OCT 2019
  !***********************************************************************
  ! If we just completed GFIC, set ITMM to 1 or 0 and return to beginning
  !***********************************************************************
  IF ((ITMM.EQ.2) .OR. (ITMM .EQ. 3)) THEN
  	IF (ITMM .EQ. 2) THEN
  	    ITMM=1
  	    !DEALLOCATE(Kpsi)
    END IF
  	IF (ITMM .EQ. 3) ITMM=0
  	max_outer=max_outer_temp
  	max_inner=max_inner_temp
  	flux=zero
    construction_inners_G = tot_nInners
    tot_nInners = 0
  	GOTO 89
  END IF


  !***********************************************************************
  ! Call wrapup to finish up post-processing and output results
  !***********************************************************************

  IF(print_conv.EQ.1) CLOSE(unit=21)
  !! ADD REMOVED - OCT 2019
  !IF (do_timing .EQ. 1) parallel_timing(4,1) = MPI_WTIME()
  !! ADD REMOVED - OCT 2019
  CALL wrapup(flux = flux,keff = keffective, unit_number = 6, suffix = "", is_final = .TRUE.)
  !! ADD REMOVED - OCT 2019
  !IF (do_timing .EQ. 1) parallel_timing(4,2) = MPI_WTIME()
  !! ADD REMOVED - OCT 2019
  !***********************************************************************
  ! Cleanup
  !***********************************************************************

  ! Flux

  DEALLOCATE(flux)

  !***********************************************************************
  ! Terminate execution
  !***********************************************************************
  !! ADD REMOVED - OCT 2019
  !IF (do_timing .EQ. 1) parallel_timing(1,2) = MPI_WTIME()
  !CALL MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
  !IF (do_timing .EQ. 1) CALL write_timing
  CALL write_timing()
  !! ADD REMOVED - OCT 2019
  !Temporary flux writing
  do i=1,num_cells
    !write(*,*)SDD_cells_l2g_G(i),flux(1,1,i,1,niter)
  end do

  CALL MPI_FINALIZE(mpi_err)
  IF (PBJrank.EQ.0) CALL stop_thor(1_li)

END PROGRAM ahot_c_ug

SUBROUTINE write_timing

  USE mpi
  USE global_variables
  USE SDD_global_variables
  IMPLICIT NONE

  INTEGER :: i, num_p, mpi_err, err_size(MPI_STATUS_SIZE), do_timing=0
  REAL*8:: print_timing(4,2)
	CHARACTER(100):: temp, temp2, temp3
  INTEGER:: values(8)
  REAL(8):: max_iter_time, min_comm_time, min_iter_time, total_comm_time, total_constr_time
  REAL(8):: total_factor_time, total_instr_time, total_IO_time, total_solve_time
  REAL(8):: avg_comm_time, avg_iter_time, comm_std_dev, iter_std_dev, max_comm_time
  INTEGER:: min_comm_loc, min_iter_loc, max_comm_loc, max_iter_loc
  CHARACTER(4) :: mode='NONE'

  CALL GET_COMMAND_ARGUMENT(2,temp)
  IF (TRIM(temp) .EQ. '-t') do_timing = 1
  flush(6)
  !! ADD REMOVED - OCT 2019
  !CALL MPI_COMM_SIZE(MPI_COMM_WORLD, num_p, mpi_err)
  num_p = 1
  rank=0
  !CALL MPI_SEND(parallel_timing, 8, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD, mpi_err)

  IF (PBJrank .EQ. 0 .AND. do_timing .EQ. 1) THEN
    DO i =0, num_p-1
      !CALL MPI_RECV(print_timing, 8, MPI_DOUBLE, i, i, MPI_COMM_WORLD, err_size, mpi_err )
      WRITE(*,*)
      WRITE(*,*)'/====================================================='
      WRITE(*,*)'| Timing Data for process: ', i
      WRITE(*,*)'|-----------------------------------------------------'
      WRITE(*,*)'| Total Measured Time:      ', print_timing(1,2)- print_timing(1,1)
      WRITE(*,*)'| Total Setup Time:         ', print_timing(2,2)- print_timing(2,1)
      WRITE(*,*)'| Total Execute Time:       ', print_timing(3,2)- print_timing(3,1)
      WRITE(*,*)'| Total Wrapup Time:        ', print_timing(4,2)- print_timing(4,1)
      WRITE(*,*)'| Total Non-Accounted Time: ',(print_timing(1,2)- print_timing(1,1))- &
                                               (print_timing(2,2)- print_timing(2,1)) - &
                                               (print_timing(3,2)- print_timing(3,1)) - &
                                               (print_timing(4,2)- print_timing(4,1))
      WRITE(*,*)'\====================================================='
    END DO
  END IF

  !SDD timings
  IF (PBJrank .EQ. 0) THEN
      !CALL MPI_RECV(print_timing, 8, MPI_DOUBLE, i, i, MPI_COMM_WORLD, err_size, mpi_err )
    ! WRITE(*,*)
    ! WRITE(*,*)'/====================================================='
    ! WRITE(*,*)'| Timing data given for root processor'
    ! WRITE(*,*)'|-----------------------------------------------------'
    ! WRITE(*,*)'| Total Construction Time:      ', Construction_end_time - Construction_start_time
    ! IF (ITMM.EQ.1) WRITE(*,*)'| Total Factorization Time:         ', Factorization_end_time - Construction_end_time
    ! WRITE(*,*)'| Total Instruction Transmit Time:       ', comm_instructions_end_time - comm_instructions_start_time
    ! WRITE(*,*)'| Total Solver Time:        ', total_solver_end_time - total_solver_start_time
    ! WRITE(*,*)'| Total Inner Time: ', solver_inner_time
    ! WRITE(*,*)'\====================================================='


    ! !TEMPORARY WRITING. REMOVE LATER!!!!!!!!!!!!!
    ! WRITE(*,*)' This should not be here after Raffi and Dylan do Godiva scaling!'
    ! WRITE(*,*)' REMOVE THIS!!!!!!!'
    ! WRITE(*,*)outer, tot_nInners, Construction_end_time - Construction_start_time &
    ! & , Factorization_end_time - Construction_end_time, &
    ! & comm_instructions_end_time - comm_instructions_start_time, &
    ! & total_solver_end_time - total_solver_start_time, &
    ! & solver_inner_time
    !!==========================================================================

    IF (ITMM.EQ.1) mode='ITMM'
    IF (ITMM.EQ.0) mode='IPBJ'
    CALL DATE_AND_TIME(VALUES=values)
    WRITE(temp3,'(I0,4(A,I0),I0)') num_proc,"-",values(2),"_",values(3),"-",values(5),"_",values(6)
    temp2 = TRIM(mesh_root)//'_'//mode//'_'//TRIM(ADJUSTL(temp3))//".time"
    OPEN (UNIT = 55, FILE = temp2, ACTION = "WRITE", STATUS = "REPLACE")
    !!85 col max width
    !!Line 1 - run info
    !! temp2
    !!Line 2 - Solution info
    !! outers, tot_inners, k
    !!Line 3 - Dylan's Timers (above)
    !!Line 4 - min/max/avg iter times + iter num
    !!Line 5 - min/max/avg solve + nums
    !!Line 6 - min/max/avg comm + nums
!   1	            2	           3	            4      	    5	          6
!1	mesh_root    	num_proc	   Month	        Day	        Hour	      Min
!2	outers	      inners	     k	            reg1	      reg2	      reg3
!3	constr. Time	Factor time	 trasnmit time	solve time	inner time	I/O Time
!4	min iter	    avg iter	   max iter	      Min #	      Max #	      std
!6	min comm	    avg comm	   max comm	      Min #	      Max #	      std

    total_constr_time = Construction_end_time - Construction_start_time
    total_factor_time = Factorization_end_time - Construction_end_time
    total_instr_time  = comm_instructions_end_time - comm_instructions_start_time
    total_solve_time  = total_solver_end_time - total_solver_start_time
    total_comm_time   = comm_time
    total_IO_time     = IOtime

    min_iter_time     = MINVAL(per_iter_time(1:MIN(tot_nInners-1, num_dbg_records)))
    min_comm_time     = MINVAL(per_comm_time(1:MIN(tot_nInners-1, num_dbg_records)))

    max_iter_time     = MAXVAL(per_iter_time(1:MIN(tot_nInners-1, num_dbg_records)))
    max_comm_time     = MAXVAL(per_comm_time(1:MIN(tot_nInners-1, num_dbg_records)))

    avg_iter_time     = SUM(per_iter_time(1:MIN(tot_nInners-1, num_dbg_records))) / MIN(tot_nInners-1, num_dbg_records)
    avg_comm_time     = SUM(per_comm_time(1:MIN(tot_nInners-1, num_dbg_records))) / MIN(tot_nInners-1, num_dbg_records)

    min_iter_loc      = MINLOC(per_iter_time(1:MIN(tot_nInners-1, num_dbg_records)), DIM=1)
    min_comm_loc      = MINLOC(per_comm_time(1:MIN(tot_nInners-1, num_dbg_records)), DIM=1)

    max_iter_loc      = MAXLOC(per_iter_time(1:MIN(tot_nInners-1, num_dbg_records)), DIM=1)
    max_comm_loc      = MAXLOC(per_comm_time(1:MIN(tot_nInners-1, num_dbg_records)), DIM=1)

    iter_std_dev      = 0
    comm_std_dev      = 0
    DO i = 1, MIN(tot_nInners-1, num_dbg_records)
      iter_std_dev = iter_std_dev + (per_iter_time(i) - avg_iter_time)**2
      comm_std_dev = comm_std_dev + (per_comm_time(i) - avg_comm_time)**2
    END DO
    iter_std_dev = SQRT(iter_std_dev / (MIN(tot_nInners-1, num_dbg_records)-1))
    comm_std_dev = SQRT(comm_std_dev / (MIN(tot_nInners-1, num_dbg_records)-1))



    WRITE(55,'(2A,I0,3X,A,3X,3(I0,A), I0)') TRIM(mesh_root), " - ", num_proc,mode, values(2),'/', &
                           & values(3),'  ', values(5),':', values(6)
    WRITE(55,'(2(I12,2X),4(ES13.6,2X))') outer, tot_nInners, k_print, reg_dbg_info(1), reg_dbg_info(2), reg_dbg_info(3)
    WRITE(55,'(6(ES13.6,2X))') total_constr_time, total_factor_time, total_instr_time, total_solve_time, &
                           &  total_comm_time, total_IO_time
    WRITE(55,'(4(ES13.6,2X),I12, 2X, I12)') avg_iter_time, min_iter_time, max_iter_time, iter_std_dev, min_iter_loc, max_iter_loc
    WRITE(55,'(4(ES13.6,2X),I12, 2X, I12)') avg_comm_time, min_comm_time, max_comm_time, comm_std_dev, min_comm_loc, max_comm_loc



  END IF
  !! ADD REMOVED - OCT 2019
END SUBROUTINE
