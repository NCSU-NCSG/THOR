!>Outer iteration module contains the subroutines necessary to perform
!>multi-group iteration operations for both eigenvalue and fixed
!>source problems.
!>It calls inner iteration
MODULE outer_iteration_module
  !***********************************************************************
  !
  ! Outer iteration module contains multi-group procedure and calls inner
  ! iteration
  !
  !***********************************************************************

  ! User derived-type modules

  USE types
  USE parameter_types
  USE filename_types
  USE vector_types
  USE cross_section_types
  USE geometry_types
  USE angle_types
  USE multindex_types
  USE global_variables
  USE wrapup_module

  ! Use modules that pertain setting up problem

  USE termination_module
  USE inner_iteration_module
  USE dump_inguess_module

  IMPLICIT NONE

CONTAINS
  !=============================================================================
  !Subroutine > outer_iteration_ext
  !=============================================================================
  !> This subroutine performs an outer iteration for an external source
  !> problem (compare Alg. 4 in the primer). Note: It does not perform
  !> thermal iterations. The group sweep is explicitly handled by a loop
  !> in the subroutine.
  SUBROUTINE outer_iteration_ext(flux,LL,U,Lf,Uf)
    !**********************************************************************
    !
    ! Subroutine outer iteration calls inner iteration and loops over all
    ! energy groups
    !
    !**********************************************************************

    ! Declare scalar flux types used globally

    REAL(kind=d_t) :: flux(num_moments_v,namom,num_cells,egmax,niter)

    ! Pass pre-computed matrices

    REAL(kind=d_t), DIMENSION(num_moments_v,num_moments_v) :: LL, U
    REAL(kind=d_t), DIMENSION(num_moments_f,num_moments_f) :: Lf, Uf


    ! Define temporary variables

    INTEGER(kind=li)               :: alloc_stat, eg, egg, q, octant, i, l, &
          order, n, m, ii, k, indx,face,f,j
    REAL(kind=d_t)                 :: t_error
    LOGICAL                        :: existence

    ! Define reflected flux

    INTEGER(kind=li)                                 :: rs,rg,mpi_err
    REAL(kind=d_t),DIMENSION(:,:,:,:,:), ALLOCATABLE :: reflected_flux

    ! distributed source

    REAL(kind=d_t) :: src(num_moments_v,namom,num_cells,egmax),SDD_errors_send(2),SDD_errors_recv(2)

    ! Prepare array reflected_flux
    rs=MAX(1_li,rside_cells)
    IF(page_refl.EQ.0_li) THEN
      rg=egmax
    ELSE
      rg=1_li
    END IF
    ALLOCATE( reflected_flux(num_moments_f,rs,8,nangle,rg),stat=alloc_stat )
    IF(alloc_stat /=0 .AND. ITMM.NE.1 .AND. ITMM.NE.0) CALL stop_thor(2_li)
    reflected_flux=0.0_d_t

    !  Allocate group-dependent maximum spatial error array

    ALLOCATE(max_error(egmax),stat=alloc_stat)
    IF(alloc_stat /= 0 .AND. ITMM.NE.1 .AND. ITMM.NE.0) CALL stop_thor(2_li)

    ! Open file to page out reflective BC if desired

    IF(page_refl.EQ.1_li) THEN

      INQUIRE(file="reflected_flux.pg",exist=existence)
      IF( existence .EQV. .TRUE.) THEN
        OPEN(unit=98,file='reflected_flux.pg',status='replace',form='unformatted',access='direct',recl=64*nangle*rs*num_moments_f)
      ELSE
        OPEN(unit=98,file='reflected_flux.pg',status='new'    ,form='unformatted',access='direct',recl=64*nangle*rs*num_moments_f)
      END IF

      DO eg = 1,egmax
        WRITE(98,rec=eg) reflected_flux(:,:,:,:,1)
      END DO

    END IF

    ! Keep track of iterations
    if (PBJrank .EQ. 0 .AND. ITMM .NE. 2 .AND. ITMM.NE.3) THEN
      WRITE(6,*) '========================================================'
      WRITE(6,*) '   Begin outer iterations.'
      WRITE(6,*) '========================================================'
    END IF

    if (PBJrank .EQ. 0 .AND. ITMM .EQ. 2) THEN
      WRITE(6,*) '========================================================'
      WRITE(6,*) '   Begin ITMM Construction Phase via GFIC'
      WRITE(6,*) '========================================================'
    END IF

    if (PBJrank .EQ. 0 .AND. ITMM .EQ. 3) THEN
      WRITE(6,*) '========================================================'
      WRITE(6,*) '   Begin IPBJ Communication Pre-Process'
      WRITE(6,*) '========================================================'
    END IF

    IF ((ITMM.EQ.2).OR.(ITMM.EQ.3)) Construction_start_time = MPI_WTIME()

    ! Begin outer iteration

    DO outer=1, max_outer

      ! Initialize the src with external source ...

      src = zero
      !If we are constructing ITMM, we leave the source as zero
      IF (ITMM.NE.2 .AND. ITMM.NE.3) THEN
      	DO eg=1,egmax
        	DO i=1,num_cells
              	! imposed internal source contribution
              	DO l=1,num_moments_v
                	src(l,1,i,eg) = src_str(cells(i)%src,eg)*src_m(l,cells(i)%src,eg)
              	END DO
        	END DO
      	END DO


      	! ... and upscattering

      	CALL compute_upscattering(flux, src)

      END IF

      DO eg=1, egmax

        ! Read binflow for group eg if necessary

        IF (page_iflw.EQ.1_li) THEN
          DO q=1,nangle
            DO octant=1,8
              DO f=1,fside_cells
                READ(97,*) face
                face = b_cells(face)%ptr
                READ(97,*) (binflx(m,face,octant,q,1),m=1,num_moments_f)
              END DO
            END DO
          END DO
        END IF

        ! Call inner iteration
        IF ((ITMM.EQ.0).OR.(ITMM.EQ.1)) solver_inner_start_time=MPI_WTIME()
        IF      (page_refl.EQ.0_li) THEN
          CALL inner_iteration(eg,flux(:,:,:,eg,niter),src(:,:,:,eg),LL,U,Lf,Uf,rs,reflected_flux(:,:,:,:,eg),.TRUE.)
        ELSE IF (page_refl.EQ.1_li) THEN
          READ (98,rec=eg) reflected_flux(:,:,:,:,1)
          CALL inner_iteration(eg,flux(:,:,:,eg,niter),src(:,:,:,eg),LL,U,Lf,Uf,rs,reflected_flux(:,:,:,:,1) ,.TRUE.)
          WRITE(98,rec=eg) reflected_flux(:,:,:,:,1)
        ELSE
          reflected_flux(:,:,:,:,1)=0.0_d_t
          CALL inner_iteration(eg,flux(:,:,:,eg,niter),src(:,:,:,eg),LL,U,Lf,Uf,rs,reflected_flux(:,:,:,:,1),.TRUE.)
        END IF
        IF ((ITMM.EQ.0).OR.(ITMM.EQ.1)) THEN
          solver_inner_end_time=MPI_WTIME()
          solver_inner_time=solver_inner_time+solver_inner_end_time-solver_inner_start_time
        END IF

        ! Compute downscattering from eg to egg

		!Don't compute downscattering if we are constucting ITMM matrices
		IF (ITMM.NE.2 .AND. ITMM.NE.3) THEN
    	    DO egg=eg+1, egmax
    	      DO i=1, num_cells
    	        DO l=0,scatt_ord
    	          DO m=0,l
    	            indx=1_li+m+(l+1_li)*l/2_li
    	            DO k=1, num_moments_v
    	              src(k,indx,i,egg) = src(k,indx,i,egg)                                 +&
    	                    scat_mult(l,m)*sigma_scat(reg2mat(cells(i)%reg),l+1,egg,eg)%xs    *&
    	                    dens_fact(cells(i)%reg)*flux(k,indx,i,eg,niter)
    	            END DO
    	          END DO
    	        END DO
    	        ! odd contributions
    	        DO l=1,scatt_ord
    	          DO m=1,l
    	            indx=neven+m+(l-1_li)*l/2_li
    	            DO k=1, num_moments_v
    	              src(k,indx,i,egg) = src(k,indx,i,egg)                                 +&
    	                    scat_mult(l,m)*sigma_scat(reg2mat(cells(i)%reg),l+1,egg,eg)%xs    *&
    	                    dens_fact(cells(i)%reg)*flux(k,indx,i,eg,niter)
    	            END DO
    	          END DO
    	        END DO
    	      END DO
	        END DO
        END IF

        ! Rewind unit=97 (finflow) if page_iflw == 1

        IF(page_iflw.EQ.1_li) REWIND(unit=97)
      END DO

      !If this was ITMM construction, exit
      IF ((ITMM.EQ.2).OR.(ITMM.EQ.3)) EXIT

      ! Compute error ...

      max_outer_error = zero

      DO eg =1,egmax
        DO i=1,num_cells
          IF( ABS(flux(1,1,i,eg,niter)) >  1.0e-12_d_t) THEN
            t_error = ABS( flux(1,1,i,eg,niter)-flux(1,1,i,eg,niter-1)) / &
                  flux(1,1,i,eg,niter)
          ELSE
            t_error = ABS( flux(1,1,i,eg,niter)-flux(1,1,i,eg,niter-1))
          END IF
          IF( t_error > max_outer_error ) THEN
            max_outer_error=t_error
          END IF
        END DO
      END DO

      ! ... and copy over iterates

      DO ii=2,niter
        DO eg=1, egmax
          DO i=1,num_cells
            DO n=1,namom
              DO l=1,num_moments_v
                flux(l,n,i,eg,ii-1)=flux(l,n,i,eg,ii)
              END DO
            END DO
          END DO
        END DO
      END DO
      if (PBJrank .EQ. 0) THEN
        WRITE(6,*)   '---------------------------------------'
        WRITE(6,103) '---itn i-itn   max error   max error---'
        WRITE(6,102) outer,tot_nInners, max_outer_error, MAXVAL(max_error),' %% '
        WRITE(6,*)   '---------------------------------------'
        flush(6)
      END IF
      IF(print_conv.EQ.1 .AND. PBJrank .EQ. 0) THEN
        WRITE(21,*)   '---------------------------------------'
        WRITE(21,103) '---itn i-itn   max error   max error---'
        WRITE(21,102) outer,tot_nInners, max_outer_error, MAXVAL(max_error),' %% '
        WRITE(21,*)   '---------------------------------------'
        flush(21)
      END IF
103   FORMAT(1X,A)
102   FORMAT(1X,2I6,2ES12.4,A)

      ! Convergence check
      SDD_errors_send(1)=MAXVAL(max_error)
      SDD_errors_send(2)=max_outer_error
      IF ((ITMM.NE.2).AND.(ITMM.NE.3)) THEN
        comm_start_time=MPI_WTIME()
        CALL MPI_ALLREDUCE(SDD_errors_send,SDD_errors_recv,2 ,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,mpi_err)
        comm_time=comm_time+MPI_WTIME()-comm_start_time
      END IF
                        
      max_error=SDD_errors_recv(1)
      max_outer_error=SDD_errors_recv(2)

      IF(most_thermal==0) THEN ! no upscattering
        IF(MAXVAL(max_error)<inner_conv) THEN
          conv_flag=1
          go to 10
        END IF
      ELSE
        IF(MAXVAL(max_error)<inner_conv .AND. max_outer_error<outer_conv) THEN
          conv_flag=1
          go to 10
        END IF
      END IF

    END DO

    outer=outer-one

10  CONTINUE

    if (PBJrank .EQ. 0 .AND. ITMM .NE. 2 .AND. ITMM.NE.3) THEN
      WRITE(6,*) '========================================================'
      WRITE(6,*) '   End outer iterations.'
      WRITE(6,*) '========================================================'
    END IF

    if (PBJrank .EQ. 0 .AND. ITMM .EQ. 2) THEN
      WRITE(6,*) '========================================================'
      WRITE(6,*) '   End GFIC ITMM Pre-Process'
      WRITE(6,*) '========================================================'
    END IF

    if (PBJrank .EQ. 0 .AND. ITMM .EQ. 3) THEN
      WRITE(6,*) '========================================================'
      WRITE(6,*) '   End IPBJ Communication Pre-Process'
      WRITE(6,*) '========================================================'
    END IF

    ! Close reflected flux file if page_ref .eq. 1

    IF(page_refl .EQ. 1_li) CLOSE(unit=98)

    IF( ALLOCATED(reflected_flux) ) DEALLOCATE(reflected_flux)

  END SUBROUTINE outer_iteration_ext

  !=============================================================================
  !Subroutine > outer_iteration_eig
  !=============================================================================
  !> This subroutine performs a power iteration for an eigenvalue
  !> problem (compare Alg. 5 in the primer). The subroutine name is a misnomer.
  !> Note: It does not perform
  !> thermal iterations. The group sweep is explicitly handled by a loop
  !> in the subroutine.
  SUBROUTINE outer_iteration_eig(flux,keff,LL,U,Lf,Uf)

    ! Pass eigenvalue

    REAL(kind=d_t), INTENT(inout) :: keff

    ! Declare angular and scalar flux types used globally

    REAL(kind=d_t) :: flux(num_moments_v,namom,num_cells,egmax,niter)

    ! Pass pre-computed matrices

    REAL(kind=d_t), DIMENSION(num_moments_v,num_moments_v) :: LL, U
    REAL(kind=d_t), DIMENSION(num_moments_f,num_moments_f) :: Lf, Uf

    ! Define temporary variables

    INTEGER(kind=li)               :: alloc_stat, eg, egg, q, octant, i, l, &
          order, ii, n, m, indx, k
    REAL(kind=d_t)                 :: fiss_den_old, fiss_den_new,SDD_fiss_den_new, &
          keff_error, keff_old, keff_new, fiss_error, fiss_dist_error(2),&
          flux_error,ts,te
    REAL(kind=d_t)                 :: t_error
    LOGICAL                        :: existence

    ! Define reflected flux

    INTEGER(kind=li)                                 :: rs,rg
    REAL(kind=d_t),DIMENSION(:,:,:,:,:), ALLOCATABLE :: reflected_flux

    ! Define fission source

    REAL(kind=d_t) :: fiss_src(num_moments_v,num_cells,2)

    ! Error mode extrapolation variables

    REAL(kind=d_t)   :: theta(3),thet
    INTEGER(kind=li) :: extra_flag
    REAL(kind=d_t)   :: a,b,c

    ! Define source that is passed into inner iteration

    REAL(kind=d_t) :: src(num_moments_v,namom,num_cells,egmax)

    ! Variables for calling wrapup

    CHARACTER(100)   :: suffix

    !Chebychev Acceleration Variables///////////////////////////////////////////
    !---------------------------------------------------------------------------
    !integer:
    !---p: Index of current chebychev iteration
    !---power_iter_count: Count of total # of corrective PIs performed
    !---cheby_pi: Number of PIs to perform per corrective cycle
    !---cheby_pi_rem: Number of PIs remaining to perform this corrective cycle
    !real:
    !---alpha, beta, gamma: Coefficients of Chebychev Acceleration term
    !---chebychev_error: Iteration error ratio / entry error ratio
    !---entry_error: Ratio of error at iteration n_entry
    !---entry_theta: Spectral radius at p=0
    !---theor_err: Maximum theoretical error reduction, f(theta, p)
    !
    !TEMPORARY // TO BE REMOVED
    !temp:
    !---temp_arg: number of command line args (to override accel mode)
    !---temp_accel: value of override argument
    !TEMPORARY // TO BE REMOVED
    !
    !---------------------------------------------------------------------------
    INTEGER (kind=li) :: p = 0_li
    INTEGER (kind=li) :: power_iter_count = 0_li, cheby_pi=5_li, cheby_pi_rem
    REAL(kind=d_t)    :: alpha = zero, beta = zero, gamma = zero
    REAL(kind=d_t)    :: chebychev_error = zero, entry_error = zero
    REAL(kind=d_t)    :: entry_theta = zero, theor_err = zero, SDD_errors_send(3)=0.0,SDD_errors_recv(3)=0.0

    INTEGER:: temp_arg, mpi_err
    CHARACTER:: temp_accel
    LOGICAL :: intermediate_flux_print

    ! Prepare array reflected_flux
    rs=MAX(1_li,rside_cells)
    IF(page_refl.EQ.0_li) THEN
      rg=egmax
    ELSE
      rg=1_li
    END IF
    ALLOCATE( reflected_flux(num_moments_f,rs,8,nangle,rg),stat=alloc_stat )
    IF(alloc_stat /=0 .AND. ITMM.NE.1 .AND. ITMM.NE.0) CALL stop_thor(2_li)
    reflected_flux=0.0_d_t

    ! Initialize fiss_src

    fiss_den_old=0.0_d_t
    fiss_den_new=0.0_d_t
    SDD_fiss_den_new=0.0_d_t
    keff=1.0_d_t
    keff_old=1.0_d_t
    keff_new=1.0_d_t

    ! Assume flat source density (for eigenvalue calculation only)

    flux=zero

    DO ii=1, niter
      DO eg=1, egmax
        DO i=1, num_cells
          flux(1,1,i,eg,ii)=one
        END DO
      END DO
    END DO

    fiss_src = zero

    DO eg=1, egmax
      DO i=1, num_cells
        DO l=1, num_moments_v
          fiss_src(l,i,1)=fiss_src(l,i,1)                                         +&
                nu(reg2mat(cells(i)%reg),eg)%xs*fiss(reg2mat(cells(i)%reg),eg)%xs  *&
                dens_fact(cells(i)%reg)*flux(l,1,i,eg,niter)
        END DO
      END DO
    END DO

    !  Allocate group-dependent maximum spatial error array

    ALLOCATE(max_error(egmax),stat=alloc_stat)
    IF(alloc_stat /= 0 .AND. ITMM.NE.1 .AND. ITMM.NE.0) CALL stop_thor(2_li)

    ! Open file to page out reflective BC if desired

    IF(page_refl.EQ.1_li) THEN

      INQUIRE(file="reflected_flux.pg",exist=existence)
      IF( existence .EQV. .TRUE.) THEN
        OPEN(unit=98,file='reflected_flux.pg',status='replace',form='unformatted',access='direct',recl=64*nangle*rs*num_moments_f)
      ELSE
        OPEN(unit=98,file='reflected_flux.pg',status='new'    ,form='unformatted',access='direct',recl=64*nangle*rs*num_moments_f)
      END IF

      DO eg = 1,egmax
        WRITE(98,rec=eg) reflected_flux(:,:,:,:,1)
      END DO

    END IF

    ! if inguess_flag == 1

    IF (inguess_flag==1) THEN

      CALL inguess_eig(flux,keff,niter)
      keff_new=keff
      keff_old=keff
      if (PBJrank .EQ. 0) THEN
        WRITE(6,*)
      END IF

      DO ii=1,niter-1
        DO eg =1,egmax
          DO i=1,num_cells
            DO n=1,namom
              DO l=1,num_moments_v
                flux(l,n,i,eg,ii) = flux(l,n,i,eg,niter)
              END DO
            END DO
          END DO
        END DO
      END DO

      fiss_src = zero

      DO eg=1, egmax
        DO i=1, num_cells
          DO l=1, num_moments_v
            fiss_src(l,i,1)=fiss_src(l,i,1)                                       +&
                  nu(reg2mat(cells(i)%reg),eg)%xs*fiss(reg2mat(cells(i)%reg),eg)%xs *&
                  dens_fact(cells(i)%reg)*flux(l,1,i,eg,niter)
          END DO
        END DO
      END DO

    END IF
    !! ADD REMOVED - OCT 2019
    !CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)
    !! ADD REMOVED - OCT 2019
    ! Keep track of iterations
    if (PBJrank .EQ. 0 .AND. ITMM .NE. 2 .AND. ITMM.NE.3) THEN
      WRITE(6,*) '========================================================'
      WRITE(6,*) '   Begin outer iterations.'
      WRITE(6,*) '========================================================'
    END IF

    if (PBJrank .EQ. 0 .AND. ITMM .EQ. 2) THEN
      WRITE(6,*) '========================================================'
      WRITE(6,*) '   Begin ITMM Construction Phase via GFIC'
      WRITE(6,*) '========================================================'
    END IF

    if (PBJrank .EQ. 0 .AND. ITMM .EQ. 3) THEN
      WRITE(6,*) '========================================================'
      WRITE(6,*) '   Begin IPBJ Communication Pre-Process'
      WRITE(6,*) '========================================================'
    END IF

    ! Set error mode extrapolation parameters
    theta=0.0_d_t
    extra_flag=0_li

    ! Compute fission density
    fiss_den_new=zero
    IF ((ITMM.NE.1).AND.(ITMM.NE.0)) THEN
      DO i=1, num_cells
        fiss_den_new=fiss_den_new+cells(i)%volume*fiss_src(1,i,1)
      END DO
    ELSE
      fiss_den_new=1.0d0
    END IF

    !===========================================================================
    ! Begin outer iteration
    !===========================================================================

    temp_arg = command_argument_count()                                         !FIXME - REMOVE AND IMPLEMENT SWITCH INTO INPUT FILE
    IF (temp_arg .EQ. 2) THEN
      CALL get_command_argument(2, temp_accel)
      IF ( temp_accel .EQ. '1') THEN
        outer_acc = 1_li
        WRITE(*,*) "NO ACCELERATION"
      ELSE IF ( temp_accel .EQ. '2') THEN
        outer_acc = 2_li
        WRITE(*,*) "FISSION SOURCE ACCELERATION"
      ELSE IF ( temp_accel .EQ. '3') THEN
        outer_acc = 3_li
        WRITE(*,*) "CHEBYCHEV ACCELERATION"
      ELSE
        WRITE(*,*) "Invalid command line acceleration parameter"
        WRITE(*,*) "Using value from input file"
      END IF
    END IF                                                                      !FIXME - REMOVE AND IMPLEMENT SWITCH INTO INPUT FILE

    DO outer=1, max_outer

      !=========================================================================
      ! Start timer
      !=========================================================================
      CALL CPU_TIME(ts)

      !=========================================================================
      ! Compute fission and set all angular moments >0 to zero...
      !=========================================================================
      src = zero
      DO eg=1,egmax
        DO i=1,num_cells
          DO l=1,num_moments_v
            src(l,1,i,eg)  =one/keff_new     *&
                  chi(reg2mat(cells(i)%reg),eg)%xs*fiss_src(l,i,1)
          END DO
        END DO
      END DO

      !=========================================================================
      ! Compute upscattering
      !=========================================================================
      CALL compute_upscattering(flux, src)

      !=========================================================================
      !Begin Group sweep
      !=========================================================================
      DO eg=1, egmax

        !=======================================================================
        ! Call inner iteration
        !=======================================================================
        IF ((ITMM.EQ.0).OR.(ITMM.EQ.1)) solver_inner_start_time=MPI_WTIME()
        IF      (page_refl.EQ.0_li) THEN
          CALL inner_iteration(eg,flux(:,:,:,eg,niter),src(:,:,:,eg),LL,U,Lf,Uf,rs,reflected_flux(:,:,:,:,eg),.TRUE.)
        ELSE IF  (page_refl.EQ.1_li) THEN
          READ (98,rec=eg) reflected_flux(:,:,:,:,1)
          CALL inner_iteration(eg,flux(:,:,:,eg,niter),src(:,:,:,eg),LL,U,Lf,Uf,rs,reflected_flux(:,:,:,:,1) ,.TRUE.)
          WRITE(98,rec=eg) reflected_flux(:,:,:,:,1)
        ELSE
          reflected_flux(:,:,:,:,1)=0.0_d_t
          CALL inner_iteration(eg,flux(:,:,:,eg,niter),src(:,:,:,eg),LL,U,Lf,Uf,rs,reflected_flux(:,:,:,:,1) ,.TRUE.)
        END IF
        IF ((ITMM.EQ.0).OR.(ITMM.EQ.1)) THEN
          solver_inner_end_time=MPI_WTIME()
          solver_inner_time=solver_inner_time+solver_inner_end_time-solver_inner_start_time
          
        END IF
        
        !=====================================================================
        ! Compute downscattering from eg to egg
        !=====================================================================
        CALL update_downscattering(eg, flux, src)
      END DO

      !========================================================================
      ! Check to see if current cycle meets acceleration criteria
      !========================================================================
      IF(outer_acc.EQ.2 .AND. outer.GE.3) THEN
        a=ABS( ( theta(1) - theta(2) )/theta(1) )
        b=ABS( ( theta(2) - theta(3) )/theta(2) )
        IF( MAX(a,b) < extol .AND. extra_flag.EQ.0_li) THEN
          extra_flag=1_li
        ELSE
          extra_flag=0_li
        END IF
      ELSE IF (outer_acc.EQ.3 .AND. outer.GT.5) THEN                                       !If more than 5 iterations, accelerate
        IF (extra_flag .EQ. 0_li .AND. theta(3) .GT. .4_d_t .AND. theta(3) .LT. 1_li) THEN !And if spectral radius is large enough
          p=1
          entry_theta = theta(3)                                              !Set sigma hat to previous guess of sp. rad
          extra_flag=1_li
        END IF
      END IF

      !========================================================================
      ! If yes, perform acceleration
      !========================================================================
      IF(extra_flag .EQ. 1_li .AND. outer_acc.EQ.2) THEN
        ! make sure that the fractional extrapolation is not larger than exmax
        thet=MIN(exmax/fiss_error,theta(3)/(1.0_d_t-theta(3)))
        ! extrapolation
        DO eg=1, egmax
          DO i=1,num_cells
            DO n=1,namom
              DO l=1,num_moments_v
                flux(l,n,i,eg,niter)=flux(l,n,i,eg,niter)+thet*(flux(l,n,i,eg,niter)-flux(l,n,i,eg,niter-1))
              END DO
            END DO
          END DO
        END DO
      ELSE IF (extra_flag .EQ. 1_li .AND. outer_acc.EQ.3) THEN
        IF (p .EQ. 1) THEN                                                      !First iteration only
          entry_error = zero                                                    !Entry error is the 2 norm of the (n - (n-1)) flux error
          DO eg=1, egmax
            DO i=1,num_cells
              entry_error=entry_error+ cells(i)%volume*(flux(1,1,i,eg,niter)- flux(1,1,i,eg,niter-1))**2
            END DO
          END DO
          entry_error = SQRT(entry_error)

          alpha = two/(two-entry_theta)                                         !Hebert Eq. B.5
          beta  = zero                                                          !^
          gamma = zero
          chebychev_error = zero
          theor_err = zero
        ELSE                                                                    !For all iteration p!=1

          gamma = dacosh((two/entry_theta) - one)                               !Hebert Eq. B.5
          alpha = (four/entry_theta) * dcosh((p-one)*gamma)/dcosh(p*gamma)      !Hebert Eq. B.5
          beta  = (one-(entry_theta/two))-(one/alpha)                           !Hebert Eq. B.5

        END IF

        chebychev_error = zero                                                  !Chebychev_err is the two norm of (n - (n-1)) flux error
        DO eg=1, egmax                                                          !divided by the enty error (error === 1 for p=1)
          DO i=1,num_cells
            chebychev_error=chebychev_error+cells(i)%volume  *&
                  (flux(1,1,i,eg,niter)- flux(1,1,i,eg,niter-1))**2
          END DO
        END DO
        chebychev_error = SQRT(chebychev_error)
        chebychev_error = chebychev_error/entry_error


        DO eg=1, egmax                                                          !Accelerate flux
          DO i=1,num_cells
            DO n=1,namom
              DO l=1,num_moments_v
                flux(l,n,i,eg,niter)=flux(l,n,i,eg,niter-1)             + &
                      alpha*(flux(l,n,i,eg,niter)-flux(l,n,i,eg,niter-1))+ &
                      beta*(flux(l,n,i,eg,niter-1)-flux(l,n,i,eg,niter-2))
              END DO
            END DO
          END DO
        END DO


        theor_err = (dcosh((REAL(p,d_t)-one)*dacosh((two - entry_theta) /&      !Theretical error
              entry_theta)))**(-one)

        IF(chebychev_error > theor_err .AND. &
              p .GE. 3*(power_iter_count+3)) THEN                                !If insufficient decrease in error,
          !write(*,*) "RESTART FLAG THROWN"                                      !reset chebychev
          extra_flag = 0_li
          !write(*,*) "BEGIN SPECTRAL RADIUS BACKTRACK"
          extra_flag = -1_li                                                    !Set flag to acceleration interrupt
          cheby_pi_rem = cheby_pi
        END IF

        p=p+1                                                                   !Increment p
      ELSE IF(extra_flag .EQ. -1_li) THEN                                       !If acceleration is interrupted
        !write(*,*) "Performing Between-Chebychev-Cycle Power Iteration"
        power_iter_count = power_iter_count + 1_li
        cheby_pi_rem = cheby_pi_rem -1_li                                       !Track remaining interrupts
        p=0
        IF (cheby_pi_rem .EQ. 0) extra_flag = 0_li                              !Reactivate acceleration
      END IF

      !========================================================================
      ! Zero fission source. Then, recompute fission source using new update
      !========================================================================
      DO i=1, num_cells
        DO l=1, num_moments_v
          fiss_src(l,i,2)=zero
        END DO
      END DO
      DO eg=1, egmax
        DO i=1, num_cells
          DO l=1, num_moments_v
            fiss_src(l,i,2)=fiss_src(l,i,2)                                        +&
                  nu(reg2mat(cells(i)%reg),eg)%xs*fiss(reg2mat(cells(i)%reg),eg)%xs *&
                  dens_fact(cells(i)%reg)*flux(l,1,i,eg,niter)
          END DO
        END DO
      END DO

      !========================================================================
      ! Compute new fission distribution and density
      !========================================================================
      fiss_den_old=fiss_den_new
      fiss_den_new=zero
      DO i=1, num_cells
        fiss_den_new=fiss_den_new+cells(i)%volume*fiss_src(1,i,2)
      END DO
!write(*,*)'A', outer,pbjrank,fiss_den_new, keff_new
      comm_start_time=MPI_WTIME()
      CALL MPI_ALLREDUCE(fiss_den_new,SDD_fiss_den_new,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpi_err)
      comm_time=comm_time+MPI_WTIME()-comm_start_time
!write(*,*)'B', outer,pbjrank,SDD_fiss_den_new, keff_new
      fiss_den_new=SDD_fiss_den_new

      !========================================================================
      ! Compute new eigenvalue based on fission densities
      !========================================================================
      keff_new=keff_old*(fiss_den_new/fiss_den_old)
      !write(*,*)pbjrank,keff_old,fiss_den_new,fiss_den_old
      !call mpi_barrier(mpi_comm_world,mpi_err)
      !stop
!write(*,*)'C', outer,pbjrank, keff_new, keff_old, fiss_den_new, fiss_den_old
!read(*,*)

      !========================================================================
      ! compute error in keff and copy new to old iterate
      !========================================================================
      keff_error=ABS(keff_new-keff_old)/keff_new
      keff_old=keff_new

      !========================================================================
      ! compute convergence based on fission source
      !========================================================================
      fiss_dist_error(1)=fiss_dist_error(2)
      IF (outer .LE. 1) fiss_dist_error(1) = 1
      fiss_dist_error(2)=0.0_d_t
      fiss_error        =0.0_d_t

      DO i=1, num_cells
        fiss_dist_error(2)=fiss_dist_error(2)+cells(i)%volume*    &
              ABS( fiss_src(1,i,2) - fiss_src(1,i,1) )
        IF(ABS(fiss_src(1,i,2)) > 1e-12)THEN
          t_error=ABS(( fiss_src(1,i,2) - fiss_src(1,i,1) )/ fiss_src(1,i,2))
        ELSE
          t_error=ABS(  fiss_src(1,i,2) - fiss_src(1,i,1) )
        END IF
        IF(t_error > fiss_error)THEN
          fiss_error=t_error
        END IF
      END DO

      !========================================================================
      ! compute the convergence based on error in group fluxes
      !========================================================================
      flux_error=0.0_d_t
      DO eg =1,egmax
        DO i=1,num_cells
          IF( ABS(flux(1,1,i,eg,niter)) > 1.0e-12_d_t) THEN
            t_error = ABS( flux(1,1,i,eg,niter)-flux(1,1,i,eg,niter-1))/&
                  flux(1,1,i,eg,niter)
          ELSE
            t_error = ABS( flux(1,1,i,eg,niter)-flux(1,1,i,eg,niter-1))
          END IF
          IF( t_error > flux_error ) THEN
            flux_error=t_error
          END IF
        END DO
      END DO

      !========================================================================
      ! compute theta (Spectral Radius)
      !========================================================================
      IF (p.EQ.0) theta(1)=theta(2)
      IF (p.EQ.0) theta(2)=theta(3)
      IF (p.EQ.0) theta(3)=fiss_dist_error(2)/fiss_dist_error(1)

      !========================================================================
      ! stop timer
      !========================================================================
      CALL CPU_TIME(te)

      !========================================================================
      ! write information for outer iteration
      !========================================================================
      if (PBJrank .EQ. 0) THEN
        WRITE(6,*)   '---------------------------------------------------------------------------------------------------'
        WRITE(6,101) '---itn i-itn        keff    err-keff    err-fiss     err-flx      Sp Rad      extrap        time---'
        WRITE(6,102) outer,tot_nInners ,keff_new, keff_error,fiss_error,flux_error,theta(3),extra_flag,te-ts,' %% '
        WRITE(6,*)   '---------------------------------------------------------------------------------------------------'
        flush(6)
        IF(print_conv.EQ.1) THEN
          WRITE(21,*)   '---------------------------------------------------------------------------------------------------'
          WRITE(21,101) '---itn i-itn        keff    err-keff    err-fiss     err-flx      Sp Rad      extrap        time---'
          WRITE(21,102) outer,tot_nInners ,keff_new, keff_error,fiss_error,flux_error,theta(3),extra_flag,te-ts,' %% '
          WRITE(21,*)   '---------------------------------------------------------------------------------------------------'
          flush(21)
        END IF
      END IF
101   FORMAT(1X,A)
102   FORMAT(1X,2I6,5ES12.4,I12,ES12.4,A)
103   FORMAT(1X,A,ES12.4)
104   FORMAT(1X,A,I5,3ES12.4)

      !========================================================================
      ! write iteration results to file if desired
      !========================================================================
      IF(dump_flag==1) THEN
        CALL dump_PI(flux,keff_old)
      END IF

      !========================================================================
      ! call wrapup even/odd
      !========================================================================
      intermediate_flux_print=.FALSE.
      IF (intermediate_flux_print.EQV..TRUE.) THEN
          IF (MOD(outer,2) .EQ. 0) THEN
            suffix = "even"
            OPEN(unit=49,file='intermediate_output_even.dat',status='unknown',action='write')
          ELSE
            suffix = "odd"
            OPEN(unit=49,file='intermediate_output_odd.dat',status='unknown',action='write')
          END IF
          CALL wrapup(flux = flux ,keff = keff_new, unit_number = 49, suffix = suffix, is_final = .FALSE.)
          CLOSE(unit=49)
      END IF

      !========================================================================
      ! check convergence and quit if criteria are satisfied
      !========================================================================

      !All reduce to make errors global
      SDD_errors_send(1)=fiss_error
      SDD_errors_send(2)=flux_error
      SDD_errors_send(3)=keff_error
      comm_start_time=MPI_WTIME()
      CALL MPI_ALLREDUCE(SDD_errors_send,SDD_errors_recv,3,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,mpi_err)
      comm_time=comm_time+MPI_WTIME()-comm_start_time
      fiss_error=SDD_errors_recv(1)
      flux_error=SDD_errors_recv(2)
      keff_error=SDD_errors_recv(3)

      IF(outer > 2 .AND. (fiss_error < outer_conv .OR. flux_error < outer_conv) .AND. &
            keff_error < k_conv) THEN
        conv_flag=1
        go to 10
      END IF

      !========================================================================
      ! Copy over iterates of flux
      !========================================================================
      DO ii=2,niter
        DO eg=1, egmax
          DO i=1,num_cells
            DO n=1,namom
              DO l=1,num_moments_v
                flux(l,n,i,eg,ii-1)=flux(l,n,i,eg,ii)
              END DO
            END DO
          END DO
        END DO
      END DO

      !========================================================================
      ! copy over the fission source: Note, that if the flux is accelerated the
      ! newly computed fission source will also be accelerated
      !========================================================================
      DO i=1, num_cells
        DO l=1, num_moments_v
          fiss_src(l,i,1)=fiss_src(l,i,2)
        END DO
      END DO

    END DO

    outer=outer-1

10  CONTINUE

    k_error=keff_error
    f_error=fiss_error
    max_outer_error = flux_error
    keff=keff_old

    if (PBJrank .EQ. 0 .AND. ITMM .NE. 2 .AND. ITMM.NE.3) THEN
      WRITE(6,*) '========================================================'
      WRITE(6,*) '   End outer iterations.'
      WRITE(6,*) '========================================================'
    END IF

    if (PBJrank .EQ. 0 .AND. ITMM .EQ. 2) THEN
      WRITE(6,*) '========================================================'
      WRITE(6,*) '   End GFIC ITMM Pre-Process'
      WRITE(6,*) '========================================================'
    END IF

    if (PBJrank .EQ. 0 .AND. ITMM .EQ. 3) THEN
      WRITE(6,*) '========================================================'
      WRITE(6,*) '   End IPBJ Communication Pre-Process'
      WRITE(6,*) '========================================================'
    END IF

    ! Close reflected flux file if page_ref .eq. 1

    IF(page_refl .EQ. 1_li) CLOSE(unit=98)

    ! Deallocate temporary arrays

    IF( ALLOCATED(reflected_flux) ) DEALLOCATE(reflected_flux)
  END SUBROUTINE outer_iteration_eig

  !=============================================================================
  !Subroutine > compute_upscattering
  !=============================================================================
  !> Computes the iteration source upscattering component based on the passed
  !> flux and  various global scattering & materials parameters.
  !> This subroutine is for both the eigenvalue and fixed source solvers
  SUBROUTINE compute_upscattering(flux, src)

    INTEGER:: eg, egg, i, l, m, k, indx
    REAL(kind=d_t) :: flux(num_moments_v,namom,num_cells,egmax,niter)
    REAL(kind=d_t) :: src(num_moments_v,namom,num_cells,egmax)

    DO eg=most_thermal,egmax     ! eg is the group that it is scattered to
      DO i=1,num_cells
        ! even contributions
        DO l=0,scatt_ord
          DO m=0,l
            indx=1_li+m+(l+1_li)*l/2_li
            DO egg=eg+1,egmax  ! egg is the group that is scattered from
              DO k=1, num_moments_v
                src(k,indx,i,eg) =  src(k,indx,i,eg) +&
                      scat_mult(l,m)*sigma_scat(reg2mat(cells(i)%reg),l+1,eg,egg)%xs  *  &
                      dens_fact(cells(i)%reg)*flux(k,indx,i,egg,niter)
              END DO
            END DO
          END DO
        END DO
        ! odd contributions
        DO l=1,scatt_ord
          DO m=1,l
            indx=neven+m+(l-1_li)*l/2_li
            DO egg=eg+1,egmax  ! egg is the group that is scattered from
              DO k=1, num_moments_v
                src(k,indx,i,eg) = src(k,indx,i,eg)                                   +&
                      scat_mult(l,m)*sigma_scat(reg2mat(cells(i)%reg),l+1,eg,egg)%xs    *&
                      dens_fact(cells(i)%reg)*flux(k,indx,i,egg,niter)
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO
  END SUBROUTINE compute_upscattering

  !=============================================================================
  !Subroutine > update_downscattering
  !=============================================================================
  !> Updates all downscattering sources for groups slower than group eg, using
  !> the flux that was just computed for group g.
  !> This subroutine is for both the eigenvalue and fixed source solvers
  SUBROUTINE update_downscattering(eg, flux, src)

    INTEGER:: eg, eg_from, i, l, m, k, index
    REAL(kind=d_t) :: flux(num_moments_v,namom,num_cells,egmax,niter)
    REAL(kind=d_t) :: src(num_moments_v,namom,num_cells,egmax)

    DO eg_from = eg + 1, egmax
      DO i = 1, num_cells
        ! Even contributions
        DO l = 0, scatt_ord
          DO m = 0, l
            index = 1_li + m + (l + 1_li) * l / 2_li
            DO k = 1, num_moments_v
              src(k,index,i,eg_from) = src(k,index,i,eg_from)                        +&
                scat_mult(l,m) * sigma_scat(reg2mat(cells(i)%reg),l+1,eg_from,eg)%xs *&
                dens_fact(cells(i)%reg) * flux(k,index,i,eg,niter)
            END DO
          END DO
        END DO
        ! odd contributions
        DO l = 1, scatt_ord
          DO m = 1, l
            index = neven + m + (l - 1_li) * l / 2_li
            DO k = 1, num_moments_v
              src(k,index,i,eg_from) = src(k,index,i,eg_from)                      +&
                scat_mult(l,m)*sigma_scat(reg2mat(cells(i)%reg),l+1,eg_from,eg)%xs *&
                dens_fact(cells(i)%reg)*flux(k,index,i,eg,niter)
            END DO
          END DO
        END DO
      END DO
    END DO
  END SUBROUTINE update_downscattering

END MODULE outer_iteration_module
!-----------------------------------------------------------------------------------------
! End
!-----------------------------------------------------------------------------------------
