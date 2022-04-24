MODULE jfnk_module
  !***********************************************************************
  !
  ! module contains all functions and subroutines required for the
  ! jacobian-free-newton-krylov(jfnk) algorithm to solve the k
  ! eigenvalue problem
  !
  ! Created by Sebastian Schunert at INL 05/26/2010
  ! *** test version (first attempt) ***
  !
  ! Revisions and Milestones
  ! --------------------------------------
  !
  ! 1. | 05/26/2010 | a. Creation of the File        | Sebastian Schunert
  ! 2. | 06/14/2010 | a. Frozen Version created      | Sebastian Schunert
  !    |            |    w/o using variables unknown |
  ! 3. | 06/22/2010 | a. Added optional PI before    | Sebastian Schunert
  !    |            |    JFNK                        |
  ! 4. | 07/13/2011 | a. Added restart and dump optio| Sebastian Schunert
  !    |            |    ns. Minor changes in variabl|
  !    |            |    e names                     |
  ! 5. | 07/21/2011 | a. Major changes of source and | Sebastian Schunert
  !    |            |    flux data structures        |
  ! 6. | 03/22/2012 | a. Make most variables global  | Sebastian Schunert
  ! 7. | 04/09/2012 | a. Change call to sweep to     | Sebastian Schunert
  !    |            |    accomodate reflected_flux   |
  !    |            | Note: all reflective BC does not work !!!
  ! 8. | 01/28/2013 | a. Change flux+src data struct.| Sebastian Schunert
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
  USE global_variables

  ! Use modules that pertain setting up problem

  USE outer_iteration_module
  USE inner_iteration_module
  USE dump_inguess_module

  IMPLICIT NONE

  ! Variables specific to jfnk: integer

  INTEGER(kind=li) :: method   ! flag to switch between 3 formulations
  ! of F(u): 1 - lagged upscattering
  !          2 - flat
  !          3 - flat with grp sweep
  INTEGER(kind=li) :: max_rest ! max # restarts of gmres/krylov
  INTEGER(kind=li) :: restart  ! max # iterations between restarts
  INTEGER(kind=li) :: max_kit  ! max # of krylov iterations - # matrix*vector
  INTEGER(kind=li) :: num_var  ! number of variables
  INTEGER(kind=li) :: iit_count! inner iteration count
  INTEGER(kind=li) :: dflag    ! determines if dump file is written
  INTEGER(kind=li) :: igflag   ! determines if initial guess is read
  INTEGER(kind=li) :: max_nit  ! maximum Newton iterations
  INTEGER(kind=li) :: grs      ! =min(#reflective faces, 1)

  ! Variables specific to jfnk: real

  REAL(kind=d_t),DIMENSION(:),ALLOCATABLE :: residual ! residual of most recent newton step
  REAL(kind=d_t),DIMENSION(:),ALLOCATABLE :: du       ! newton step, solution increment
  REAL(kind=d_t) :: jac_eps                           ! finite difference perturbation
  REAL(kind=d_t) :: nres                              ! newton residual nres=||residual||_p
  REAL(kind=d_t) :: nit_conv                          ! newton iteration convergence


CONTAINS

  !-----------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------

  SUBROUTINE set_jfnk
    !**********************************************************************
    !
    ! Subroutine jfnk_setup allocates and sets the jfnk variables
    !
    !**********************************************************************

    ! local variables

    INTEGER(kind=li) :: alloc_stat

    ! set method, krest and kit

    method = rd_method
    restart = rd_restart
    max_kit = rd_max_kit
    max_rest = INT(max_kit/restart)

    ! calculate length of unknown vector

    num_var = egmax*num_cells*num_moments_v*namom

    ! allocate residual and du

    IF(.NOT.ALLOCATED(du)) ALLOCATE(du(num_var+1),stat=alloc_stat)
    IF(alloc_stat /= 0) CALL raise_fatal_error("*** Not enough memory ***")
    IF(.NOT.ALLOCATED(residual)) ALLOCATE(residual(num_var+1),stat=alloc_stat)
    IF(alloc_stat /= 0) CALL raise_fatal_error("*** Not enough memory ***")

    ! set dflag, igflag

    dflag =dump_flag
    igflag=inguess_flag

    ! set jac_eps

    jac_eps = 1.4901d-08

    ! set max_nit and nit_conv

    max_nit=max_outer
    nit_conv=outer_conv

    ! make sure that niter .eq. 2

    IF (niter .NE. 2) CALL raise_fatal_error("JFNK module requires niter to be equal to 2. This is a coding mistake!")

    ! Allocate max_error

    ALLOCATE(max_error(egmax),stat=alloc_stat)
    IF(alloc_stat /= 0) CALL raise_fatal_error("*** Not enough memory ***")

    ! Set grs

    grs=MAX(1_li,rside_cells)

  END SUBROUTINE set_jfnk

  !-----------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------

  SUBROUTINE clean_jfnk
    !**********************************************************************
    !
    ! subroutine clean_jfnk deallocates variables used in jfnk algorithm
    !
    !**********************************************************************

    ! deallocate residual and du

    IF(ALLOCATED(du)) DEALLOCATE(du)
    IF(ALLOCATED(residual)) DEALLOCATE(residual)

  END SUBROUTINE clean_jfnk

  SUBROUTINE do_jfnk(flux,keff,LL,U,Lf,Uf)
    !**********************************************************************
    !
    ! Subroutine do_jfnk is the driver routine for the jfnk method
    !
    !**********************************************************************

    ! Declare scalar flux types used globally+keff

    REAL(kind=d_t) :: flux(num_moments_v,namom,num_cells,egmax,niter)
    REAL(kind=d_t) :: keff

    ! Define temporary variables

    REAL(kind=d_t), DIMENSION(num_moments_v,num_moments_v) :: LL, U
    REAL(kind=d_t), DIMENSION(num_moments_f,num_moments_f) :: Lf, Uf

    ! Local iteration counter

    INTEGER(kind=li) :: kit

    ! Local variables

    REAL(kind=d_t) :: kerr,tol,tfiss1,tfiss2
    INTEGER(kind=li) :: i,l,eg,m,n
    INTEGER(kind=li) :: one,two,three
    REAL(kind=d_t) :: tot_vol
    REAL(kind=d_t) :: ts,te,keff_error,max_outer_error,tmp

    ! Reflective flux in case ipiter > 0

    REAL(kind=d_t), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: reflected_flux

    ! set one two three

    one = 1_li
    two = 2_li
    three = 3_li

    ! Set cell volumes
    tot_vol = 0.0_d_t
    DO i=1, num_cells
      tot_vol=tot_vol+cells(i)%volume
    END DO

    ! start newton iteration

    nit       = 0
    tot_kit   = 0
    iit_count = 0

    ! set up initial guess: perform initial power iteration ??

    IF(igflag==1) THEN
      CALL inguess_eig(flux,keff,one)
    ELSE
      flux=0.0_d_t
      DO eg =1,egmax
        DO i=1,num_cells
          flux(1,1,i,eg,1)=1.0_d_t
          flux(1,1,i,eg,2)=1.0_d_t
        END DO
      END DO
      keff = 1.0_d_t
    END IF
    ! Perform ipow initial power iterations
    IF (ipow .GT. 0 .AND. rank .EQ. 0) THEN
      WRITE(stdout_unit,*)
      WRITE(stdout_unit,*) '========================================================'
      WRITE(stdout_unit,*) '   Begin Initial Power Iterations.'
      WRITE(stdout_unit,*) '========================================================'
    END IF

    ! note: we do not need to duplicate reflected_flux on all processors
    ! but for simplicity we do here.
    IF (ipow .GT. 0) THEN
      ALLOCATE(reflected_flux(num_moments_f,grs,8,nangle,egmax))
      reflected_flux = 0.0_d_t
    END IF

    DO m = 1,ipow

      ! start timer

      CALL CPU_TIME(ts)

      ! call single outer iteration - it does not update keff

      CALL single_outer_iteration(flux,keff,reflected_flux,LL,U,Lf,Uf,.TRUE.)

      ! update keff

      tfiss1 = tfissions (flux,two)
      tfiss2 = tfissions (flux,one)
      keff_error = ABS(keff * tfiss1 / tfiss2-keff)/keff
      keff = keff * tfiss1 / tfiss2

      ! copy iterate #2 into working iterate #1
      max_outer_error=0.0_d_t
      DO eg =1,egmax
        DO i=1,num_cells
          tmp=ABS(flux(1,1,i,eg,2)-flux(1,1,i,eg,1)) / &
                flux(1,1,i,eg,1)
          IF(tmp>max_outer_error) max_outer_error=tmp
          DO n =1,namom
            DO l=1,num_moments_v
              flux(l,n,i,eg,1)=flux(l,n,i,eg,2)
            END DO
          END DO
        END DO
      END DO

      ! stop timer

      CALL CPU_TIME(te)

      ! write convergence monitor
      IF (rank .EQ. 0) THEN
        WRITE(stdout_unit,*)   '---------------------------------------------------------------------------'
        WRITE(stdout_unit,401) '---itn        keff    err-keff     err-flx        time---'
        WRITE(stdout_unit,402) m,keff, keff_error,max_outer_error,te-ts,' %i'
        WRITE(stdout_unit,*)   '---------------------------------------------------------------------------'
        IF(print_conv.EQ.1) THEN
          WRITE(21,*)   '---------------------------------------------------------------------------'
          WRITE(21,401) '---itn        keff    err-keff     err-flx        time---'
          WRITE(21,402) m,keff, keff_error,max_outer_error,te-ts,' %i'
          WRITE(21,*)   '---------------------------------------------------------------------------'
        END IF
      END IF
401   FORMAT(1X,A)
402   FORMAT(1X,I6,4ES12.4,A)

    END DO

    IF( ALLOCATED(reflected_flux) ) DEALLOCATE(reflected_flux)

    IF (ipow .GT. 0 .AND. rank .EQ. 0) THEN
      WRITE(stdout_unit,*)
      WRITE(stdout_unit,*) '========================================================'
      WRITE(stdout_unit,*) '   End Initial Power Iterations.'
      WRITE(stdout_unit,*) '========================================================'
    END IF

    IF (rank .EQ. 0) THEN
      WRITE(stdout_unit,*)
      WRITE(stdout_unit,*) '========================================================'
      WRITE(stdout_unit,*) '   Begin JFNK iterations.'
      WRITE(stdout_unit,*) '========================================================'
    END IF

    ! evaluate residual

    residual=evaluate_residual(flux,keff,LL,U,Lf,Uf)

    ! compute current p-norm of newton residual

    nres = pnorm(residual,12_li)

    ! compute forcing factor

    tol = forcing_factor()

    ! increment newton iteration counter

    nit = nit + 1

    DO i = 1, max_nit

      ! Start Timer

      CALL CPU_TIME(ts)

      ! call krylov solver

      CALL krylovM(-residual,du,tol,kit,kerr,&
                                ! overhead variables
            flux,keff,LL,U,Lf,Uf)

      ! update unknown

      CALL newton_update(flux,keff)

      ! dump result into file if dump_flag==1

      IF(dflag==1) THEN
        CALL dump_jfnk(flux,keff)
      END IF

      ! evaluate residual

      residual=evaluate_residual(flux,keff,LL,U,Lf,Uf)

      ! compute current p-norm of newton residual

      nres = pnorm(residual,12_li)

      ! compute forcing factor

      tol = forcing_factor()

      ! stop timer

      CALL CPU_TIME(te)

      ! output convergence progress
      IF (rank .EQ. 0) THEN
        WRITE(stdout_unit,*)   '---------------------------------------------------------------------------'
        WRITE(stdout_unit,102) '---nitn  kitn        keff     max-res        time---'
        WRITE(stdout_unit,101) nit,kit,keff,nres,te-ts," %n"
        WRITE(stdout_unit,*)   '---------------------------------------------------------------------------'
        IF(print_conv.EQ.1) THEN
          WRITE(21,*)   '---------------------------------------------------------------------------'
          WRITE(21,102) '---nitn  kitn        keff     max-res        time---'
          WRITE(21,101) nit,kit,keff,nres,te-ts," %n"
          WRITE(21,*)   '---------------------------------------------------------------------------'
        END IF
      END IF
101   FORMAT (1X,I7,I6,3ES12.4,A)
102   FORMAT (1X,A)

      ! Check for convergence

      IF(nres<nit_conv .and. rank .EQ. 0) THEN
        WRITE(stdout_unit,*)   '---------------------------------------------------------------------------'
        WRITE(stdout_unit,104) 'Newton Iteration Convergence achieved after',nit,' iterations.'
        WRITE(stdout_unit,103) 'with final residual ',nres
        WRITE(stdout_unit,*)   '---------------------------------------------------------------------------'
        conv_flag=1
      END IF
      IF(nres<nit_conv)EXIT
103   FORMAT (1X,A,ES12.4)
104   FORMAT (1X,A44,I4,A)

      ! increment iteration counters

      nit = nit + 1
      tot_kit = tot_kit + kit


    END DO

    nit=nit-1

    ! normalize eigenmode

    CALL norm_eigmode(flux,tot_vol)

    IF (rank .EQ. 0) THEN
      WRITE(stdout_unit,*) '========================================================'
      WRITE(stdout_unit,*) '   End JFNK iterations.'
      WRITE(stdout_unit,*) '========================================================'
    END IF

    ! newton iteration end, copy unknown to flux

    DO eg = 1, egmax
      DO i = 1, num_cells
        DO n=1,namom
          DO l = 1, num_moments_v
            flux(l,n,i,eg,2)=flux(l,n,i,eg,1)
          END DO
        END DO
      END DO
    END DO

    ! assign total # inner iterations

    tot_nInners = iit_count

    ! calculate max_error

    CALL get_max_error


  END SUBROUTINE do_jfnk

  !-----------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------

  FUNCTION evaluate_residual(flux,keff,LL,U,Lf,Uf)
    !**********************************************************************
    !
    ! Function evaluate_residual evaluates the non-linear function
    ! r = F(u)
    !
    !**********************************************************************

    ! Pass flux

    REAL(kind=d_t) :: flux(num_moments_v,namom,num_cells,egmax,niter)

    ! pass multiplication factor

    REAL(kind=d_t) :: keff

    ! Precomputed matrices

    REAL(kind=d_t), DIMENSION(num_moments_v,num_moments_v) :: LL, U
    REAL(kind=d_t), DIMENSION(num_moments_f,num_moments_f) :: Lf, Uf

    ! define return value

    REAL(kind=d_t), DIMENSION(num_var+1) :: evaluate_residual

    ! Local variables

    INTEGER(kind=li) :: eg,i,l,alloc_stat,n,indx,zero,one,two,three
    REAL(kind=d_t)   :: f1,f2

    ! Source: src is a general source used for all methods

    REAL(kind=d_t) :: src(num_moments_v,namom,num_cells,egmax)

    ! Reflected flux array

    REAL(kind=d_t),ALLOCATABLE ::  reflected_flux(:,:,:,:,:)

    ! if method = 2 or 3 allocate and set outward_normal. allocate src

    IF (method==1) THEN

      ! initialize reflected_flux

      ALLOCATE(reflected_flux(num_moments_f,grs,8,nangle,egmax),stat=alloc_stat)
      IF(alloc_stat /=0) CALL raise_fatal_error("*** Not enough memory ***")
      reflected_flux = 0.0_d_t

    ELSE IF(method==2 .OR. method==3) THEN

      ! initialize reflected_flux

      ALLOCATE(reflected_flux(num_moments_f,grs,8,nangle,1),stat=alloc_stat)
      IF(alloc_stat /=0) CALL raise_fatal_error("*** Not enough memory ***")
      reflected_flux = 0.0_d_t

      ! initialize source

      src = 0.0_d_t

    END IF

    ! set 1,2,3

    zero = 0_li
    one = 1_li
    two = 2_li
    three = 3_li

    ! initialize new iterates

    IF(method .EQ. 1) THEN
      DO eg = 1, egmax
        DO i = 1, num_cells
          DO n = 1, namom
            DO l = 1, num_moments_v
              flux(l,n,i,eg,2)=flux(l,n,i,eg,1)
            END DO
          END DO
        END DO
      END DO
    ELSE IF (method .EQ. 2 .OR. method .EQ. 3) THEN
      DO eg = 1, egmax
        DO i = 1, num_cells
          DO n = 1, namom
            DO l = 1, num_moments_v
              flux(l,n,i,eg,2) = 0.0_d_t
            END DO
          END DO
        END DO
      END DO
    END IF


    ! evaluate residual based on the chosen method

    IF (method==1) THEN

      ! Outer iteration with lagged upscattering

      CALL single_outer_iteration(flux,keff,reflected_flux,LL,U,Lf,Uf,.FALSE.)

    ELSE IF (method==2) THEN

      ! Flat - Sweep on fixed source
      ! 1. Build fixed source

      CALL bfsrc (flux,src,keff,one)

      DO eg = 1, egmax
        CALL bscsrc (flux,src,one,one,eg)
      END DO

      ! 2. Sweep on source
      DO eg = 1, egmax
        reflected_flux=0.0_d_t
        CALL sweep(eg,flux(:,:,:,eg,2),src(:,:,:,eg),grs,reflected_flux,LL,U,Lf,Uf)
      END DO

    ELSE IF (method==3) THEN

      ! Sweep through groups with updated downscattering
      ! 1. Build fission source

      CALL bfsrc (flux,src,keff,one)

      ! Calculate self-scattering source into group 1
      eg = 1
      CALL bscsrc (flux,src,one,one,eg)

      ! Sweep on fixed source in group 1
      reflected_flux=0.0_d_t
      CALL sweep (eg,flux(:,:,:,eg,2),src(:,:,:,eg),grs,reflected_flux,LL,U,Lf,Uf)

      ! 2. Sweep through grps 2..egmax and re-evaluate downscatter source
      DO eg = 2, egmax
        CALL bscsrc (flux,src,two,one,eg)
        !
        reflected_flux=0.0_d_t
        CALL sweep(eg,flux(:,:,:,eg,2),src(:,:,:,eg),grs,reflected_flux,LL,U,Lf,Uf)
      END DO
    ELSE
      CALL raise_fatal_error("Method needs to be 1, 2 or 3.")
    END IF

    ! Compute evaluate_residual

    DO eg = 1, egmax
      DO i = 1, num_cells
        DO n=1,namom
          DO l = 1, num_moments_v
            indx = (eg-1)*num_cells*namom*num_moments_v+&
                  (i-1)           *namom*num_moments_v+&
                  (n-1)                 *num_moments_v+&
                  l
            evaluate_residual(indx) =  flux(l,n,i,eg,1)-flux(l,n,i,eg,2)
          END DO
        END DO
      END DO
    END DO

    ! assign [num_var+1] using FR approach

    f1 = tfissions (flux,one)
    f2 = tfissions (flux,two)

    evaluate_residual(num_var+1) = keff*(1.0_d_t-f2/f1)

    ! Deallocate

    IF(method==2 .OR. method==3) THEN

      ! deallocate reflected flux
      DEALLOCATE(reflected_flux)

    END IF

  END FUNCTION evaluate_residual

  !-----------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------

  SUBROUTINE bfsrc (flx,src,keff,flag)
    !**********************************************************************
    !
    ! Subroutine bfsrc adds the fission source to src.
    !
    ! upon outout src := src + fsrc
    ! this subroutine can be easily modified to house upscattering
    ! as well
    !
    !**********************************************************************

    ! Pass input parameters

    INTEGER(kind=li), INTENT(in) :: flag

    ! Pass flux and source and keff

    REAL(kind=d_t) :: flx (num_moments_v,namom,num_cells,egmax,niter)
    REAL(kind=d_t) :: src (num_moments_v,namom,num_cells,egmax)
    REAL(kind=d_t) :: keff

    ! Define fission source

    REAL(kind=d_t) :: tfsrc(num_moments_v,num_cells)

    ! Local variables

    INTEGER(kind=li) :: eg,m,i,l,egg,k,indx,mat_indx

    ! Initialize tfsrc

    tfsrc = 0.0_d_t

    ! Compute total fission source

    DO i = 1, num_cells
      mat_indx=material_ids(reg2mat(cells(i)%reg))
      DO eg = 1, egmax
        DO l = 1, num_moments_v
          tfsrc(l,i) = tfsrc(l,i) + xs_mat(mat_indx)%nu(eg) &
            *xs_mat(mat_indx)%sigma_f(eg)* &
                dens_fact(cells(i)%reg)*flx(l,1,i,eg,flag)
        END DO
      END DO
    END DO

    ! Add to src

    DO i = 1, num_cells
      mat_indx=material_ids(reg2mat(cells(i)%reg))
      DO eg = 1, egmax
        DO l = 1, num_moments_v
          src(l,1,i,eg) = src(l,1,i,eg) + xs_mat(mat_indx)%chi(eg)*tfsrc(l,i)/keff
        END DO
      END DO
    END DO

    ! Add upscattering

    DO i=1,num_cells
      mat_indx=material_ids(reg2mat(cells(i)%reg))
      DO eg=most_thermal,egmax    ! eg is the group that it is scattered to
        ! even contribution
        DO l=0,scatt_ord
          DO m=0,l
            indx=1_li+m+(l+1_li)*l/2_li
            DO egg=eg+1,egmax  ! egg is the group that is scattered from
              DO k=1, num_moments_v
                src(k,indx,i,eg) = src(k,indx,i,eg)                                  +&
                      scat_mult(l,m)*xs_mat(mat_indx)%sigma_scat(l+1,eg,egg)   *&
                      dens_fact(cells(i)%reg)*flx(k,indx,i,egg,flag)
              END DO
            END DO
          END DO
        END DO
        ! odd contribution
        DO l=1,scatt_ord
          DO m=1,l
            indx=neven+m+(l-1_li)*l/2_li
            DO egg=eg+1,egmax  ! egg is the group that is scattered from
              DO k=1, num_moments_v
                src(k,indx,i,eg)=src(k,indx,i,eg)                                    +&
                      scat_mult(l,m)*xs_mat(mat_indx)%sigma_scat(l+1,eg,egg)*&
                      dens_fact(cells(i)%reg)*flx(k,indx,i,egg,flag)
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE bfsrc

  !-----------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------

  SUBROUTINE bscsrc (flx,src,flag1,flag2,egg)
    !**********************************************************************
    !
    ! Subroutine bscsrc adds the down and and selfscattering source
    ! into group egg to src so upon output
    !
    ! src = src + bscsrc
    !
    ! the parameter flag switches between using old or new
    ! iterate value of the scalar flux
    !
    ! flag1 - flag for downscattering source
    ! flag2 - flag for self-scattering source
    !
    !**********************************************************************

    ! Pass input parameters

    INTEGER(kind=li), INTENT(in) :: flag1,flag2

    ! Pass flux and source and keff

    REAL(kind=d_t) :: flx (num_moments_v,namom,num_cells,egmax,niter)
    REAL(kind=d_t) :: src (num_moments_v,namom,num_cells,egmax)

    ! Pass group

    INTEGER(kind=li) :: egg

    ! Local variables

    INTEGER(kind=li) :: eg,i,l,m,k,indx,mat_indx

    ! Start adding downscattering source

    DO i = 1, num_cells
      mat_indx=material_ids(reg2mat(cells(i)%reg))
      DO eg = 1, egg-1
        ! even contribution
        DO l=0, scatt_ord
          DO m = 0,l
            indx=1_li+m+(l+1_li)*l/2_li
            DO k = 1, num_moments_v
              src(k,indx,i,egg)=src(k,indx,i,egg)                               +&
                    scat_mult(l,m)*xs_mat(mat_indx)%sigma_scat(l+1,egg,eg)*&
                    dens_fact(cells(i)%reg)*flx(k,indx,i,eg,flag1)
            END DO
          END DO
        END DO
        ! odd contribution
        DO l=1, scatt_ord
          DO m = 1,l
            indx=neven+m+(l-1_li)*l/2_li
            DO k = 1, num_moments_v
              src(k,indx,i,egg)=src(k,indx,i,egg)                               +&
                    scat_mult(l,m)*xs_mat(mat_indx)%sigma_scat(l+1,egg,eg)*&
                    dens_fact(cells(i)%reg)*flx(k,indx,i,eg,flag1)
            END DO
          END DO
        END DO
      END DO
    END DO

    ! Add selfscattering source

    DO i = 1, num_cells
      mat_indx=material_ids(reg2mat(cells(i)%reg))
      ! Even contributions
      DO l = 0, scatt_ord
        DO m = 0,l
          indx=1_li+m+(l+1_li)*l/2_li
          DO k = 1, num_moments_v
            src(k,indx,i,egg)=src(k,indx,i,egg)                                +&
                  scat_mult(l,m)*xs_mat(mat_indx)%sigma_scat(l+1,egg,egg)*&
                  dens_fact(cells(i)%reg)*flx(k,indx,i,eg,flag2)
          END DO
        END DO
      END DO
      ! odd contributions
      DO l = 1, scatt_ord
        DO m = 1,l
          indx=neven+m+(l-1_li)*l/2_li
          DO k = 1, num_moments_v
            src(k,indx,i,egg)=src(k,indx,i,egg)                                +&
                  scat_mult(l,m)*xs_mat(mat_indx)%sigma_scat(l+1,egg,egg)*&
                  dens_fact(cells(i)%reg)*flx(k,indx,i,eg,flag2)
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE bscsrc

  !-----------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------

  FUNCTION tfissions (flux,flag)
    !**********************************************************************
    !
    ! function tfission computes the total fissions in the domain
    !
    ! Sum_{#groups,#cells} nsigf(grp,cell#)*phi_000
    !
    !**********************************************************************

    ! Pass input parameters

    INTEGER(kind=li), INTENT(in) :: flag

    ! pass flux

    REAL(kind=d_t) :: flux(num_moments_v,namom,num_cells,egmax,niter)

    ! define return value

    REAL(kind=d_t) :: tfissions

    ! local variables

    INTEGER(kind=li) :: eg,i,mat_indx

    ! initialize tfissions

    tfissions = 0.0_d_t

    ! compute tfission

    DO i = 1, num_cells
      mat_indx=material_ids(reg2mat(cells(i)%reg))
      DO eg = 1, egmax
        tfissions = tfissions+ &
              xs_mat(mat_indx)%nu(eg)*xs_mat(mat_indx)%sigma_f(eg)* &
              cells(i)%volume*dens_fact(cells(i)%reg)*flux(1,1,i,eg,flag)
      END DO
    END DO

  END FUNCTION tfissions

  !-----------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------

  SUBROUTINE krylovM(rhs,solution,tol,kit,error,          &
                                ! overhead variables
        flux,keff,LL,U,Lf,Uf)
    !**********************************************************************
    !
    ! Subroutine krylovM is the interface to sparskit
    !
    !**********************************************************************

    ! Pass input variables

    REAL(kind=d_t),DIMENSION(num_var+1) :: rhs,solution
    REAL(kind=d_t) :: tol
    INTEGER(kind=li) :: kit
    REAL(kind=d_t) :: error

    ! Pass flux

    REAL(kind=d_t) :: flux(num_moments_v,namom,num_cells,egmax,niter)

    ! pass multiplication factor

    REAL(kind=d_t) :: keff

    ! Pre-computed matrices

    REAL(kind=d_t), DIMENSION(num_moments_v,num_moments_v) :: LL, U
    REAL(kind=d_t), DIMENSION(num_moments_f,num_moments_f) :: Lf, Uf

    ! Local variables
    ! ws    - work-size, dimension of work array
    ! ipar  - integer parameter krylov method
    ! work  - work array
    ! fpar  - real parameters krylov method

    INTEGER(kind=li) :: ws
    INTEGER,DIMENSION(16) :: ipar
    REAL(kind=d_t),DIMENSION(:),ALLOCATABLE :: work
    REAL(kind=d_t),DIMENSION(16) :: fpar
    INTEGER :: alloc_stat
    REAL(kind=d_t) :: ts,te,nrhs
    INTEGER :: k

    ! Calculate ws(for gmres)

    ws = (num_var+4)*(restart+2)+(restart+1)*restart/2

    ! Allocate and initialize  work

    ALLOCATE(work(ws),stat=alloc_stat)
    IF(alloc_stat /= 0) CALL raise_fatal_error("*** Not enough memory ***")
    work = 0.0_d_t

    ! Initialize solution, ipar and fpar

    solution = 0.0_d_t
    ipar = 0_li
    fpar = 0.0_d_t

    ! Set ipar array

    ipar(1)=0          ! indicates start of iterative solver
    ipar(2)=0          ! 0/1/2/3 --> no/left/right/both preconditioning
    ipar(3)=2          ! see iters.f for definitions
    ipar(4)=ws         ! size of workspace, minimum defined in gmres subroutine
    ipar(5)=restart    ! number of GMRES iterations between restarts
    ipar(6)=max_kit    ! total number of matvecs allowed

    fpar(1)=tol        ! relative tolerance
    fpar(2)=0.0_d_t    ! absolute tolerance
    fpar(11)=0.0_d_t   ! Clearing FLOPS counter

    ! Compute the norm of the rhs

    nrhs =pnorm(rhs,2_li)

    ! Print header for Krylov convergence monitor

    IF(rank .EQ. 0) WRITE(stdout_unit,102) '---itn         err    trgt-res        time ---'

    ! Start gmres iteration

    k=0
    DO

      ! increment counter

      k=k+1

      ! start timer

      CALL CPU_TIME(ts)

      CALL gmres(num_var+1,rhs,solution,ipar,fpar,work)

      ! Go through output cases

      IF(ipar(1).EQ.1) THEN           ! this means code wants Ax
        work(ipar(9):ipar(9)+num_var)=jacvec(work(ipar(8):ipar(8)+num_var),&
                                ! overhead variables
              flux,keff,LL,U,Lf,Uf)
      ELSEIF(ipar(1).EQ.2) THEN    ! this means code wants A^T x
        CALL raise_fatal_error("GMRES error")
      ELSEIF(ipar(1).EQ.3) THEN    ! this means code wants P(l)^{-1} z
        CALL raise_fatal_error("GMRES error")
      ELSEIF(ipar(1).EQ.3) THEN    ! this means code wants P(l)^{-T} z
        CALL raise_fatal_error("GMRES error")
      ELSEIF(ipar(1).EQ.5) THEN    ! this means code wants P(r)^{-1} z
        CALL raise_fatal_error("GMRES error")
      ELSEIF(ipar(1).EQ.6) THEN    ! this means code wants P(r)^{-T} z
        CALL raise_fatal_error("GMRES error")
      ELSEIF(ipar(1).EQ.10) THEN   ! call self-supplied stopping test
        CALL raise_fatal_error("GMRES error")
      ELSEIF(ipar(1).GT.0) THEN    ! shouldn't happen
        CALL raise_fatal_error("GMRES error")
      ELSEIF(ipar(1).EQ.0) THEN    ! successful solve
        error=fpar(6)/nres         ! pass back error
        kit=ipar(7)                ! pass back number of iters
        EXIT
      ELSEIF(ipar(1).EQ.-1) THEN   ! maximum number of its reached
        error=fpar(6)/nres         ! pass back error
        kit=ipar(7)                ! pass back number of iters
        EXIT
      ELSEIF(ipar(1).EQ.-2) THEN   ! insufficient workspace
        CALL raise_fatal_error("GMRES error")
      ELSEIF(ipar(1).LT.-2) THEN   ! some other error, investigate further
        WRITE(stdout_unit,*) 'SPARSKIT error is ', ipar(1)
      ENDIF

      ! stop timer

      CALL CPU_TIME(te)

      ! write convergence monitor
      IF (rank .EQ. 0) WRITE(stdout_unit,101) k,fpar(5),tol*nrhs,te-ts,' %%k'

    END DO

    ! End gmres iteration

    ! Deallocate work

    DEALLOCATE(work)

101 FORMAT(1X,I6,3ES12.4,A)
102 FORMAT(1X,A)
  END SUBROUTINE krylovM

  !-----------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------

  FUNCTION jacvec(u_vec,&
                                ! overhead variables
        flux,keff,LL,U,Lf,Uf)

    ! Pass input variables

    REAL(kind=d_t), DIMENSION(num_var+1) :: u_vec

    ! Pass flux

    REAL(kind=d_t) :: flux(num_moments_v,namom,num_cells,egmax,niter)

    ! pass multiplication factor

    REAL(kind=d_t) :: keff

    ! Pre-computed matrices

    REAL(kind=d_t), DIMENSION(num_moments_v,num_moments_v) :: LL, U
    REAL(kind=d_t), DIMENSION(num_moments_f,num_moments_f) :: Lf, Uf

    ! Define return

    REAL(kind=d_t), DIMENSION(num_var+1) :: jacvec

    ! define local

    REAL(kind=d_t)   :: dflux(num_moments_v,namom,num_cells,egmax,niter)
    REAL(kind=d_t)   :: dkeff
    INTEGER(kind=li) :: eg,i,l,one,indx,n

    ! allocate dflx and initialize

    one = 1_li

    ! compute dflx and dkeff from

    DO eg=1,egmax
      DO i=1,num_cells
        DO n=1,namom
          DO l=1,num_moments_v
            indx = (eg-1)*num_cells*namom*num_moments_v+&
                  (i-1)           *namom*num_moments_v+&
                  (n-1)                 *num_moments_v+&
                  l
            dflux(l,n,i,eg,1)=flux(l,n,i,eg,1)+jac_eps*u_vec(indx)
          END DO
        END DO
      END DO
    END DO
    dkeff = keff+jac_eps*u_vec(num_var+1)


    ! Compute jacvec

    jacvec = ( evaluate_residual(dflux,dkeff,LL,U,Lf,Uf) - residual ) &
          / jac_eps

  END FUNCTION jacvec

  !-----------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------

  FUNCTION converged(conv)
    !**********************************************************************
    !
    ! Function converged determines whether newton iteration is converged
    ! or not
    !
    !**********************************************************************

    ! pass input variables

    REAL(kind=d_t) :: conv
    LOGICAL :: converged

    ! convergence test

    converged = .FALSE.
    IF(MAXVAL(ABS(residual))<conv) converged = .TRUE.

  END FUNCTION converged

  !-----------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------

  SUBROUTINE newton_update(flux,keff)
    !**********************************************************************
    !
    ! Subroutine computes the new vector of unknowns using the du from last
    ! newton step
    !
    !**********************************************************************

    ! pass arguments

    REAL(kind=d_t) :: flux(num_moments_v,namom,num_cells,egmax,niter)
    REAL(kind=d_t) :: keff

    ! local variables

    INTEGER(kind=li) :: eg,i,l,one,indx,n

    ! standard update

    one = 1_li

    ! update

    DO eg = 1, egmax
      DO i = 1, num_cells
        DO n=1,namom
          DO l = 1, num_moments_v
            indx = (eg-1)*num_cells*namom*num_moments_v+&
                  (i-1)           *namom*num_moments_v+&
                  (n-1)                 *num_moments_v+&
                  l
            flux(l,n,i,eg,one) = flux(l,n,i,eg,one) + du(indx)
          END DO
        END DO
      END DO
    END DO
    keff = keff + du(num_var+1)

  END SUBROUTINE newton_update

  !-----------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------

  FUNCTION forcing_factor()
    !**********************************************************************
    !
    ! Function forcing_factor computes the forcing factor
    !
    ! cnres - current norm of residual vector
    ! cnit  - current newton iteration
    !
    !**********************************************************************

    ! pass input arguments

    REAL(kind=d_t) :: forcing_factor

    ! easiest implementation

    !      dummy = 1/(real(cnit,d_t)+2.0_d_t)
    !      forcing_factor = min(dummy,cnres)
    forcing_factor = 0.01_d_t
  END FUNCTION forcing_factor

  !-----------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------

  FUNCTION pnorm(cresidual,pp)
    !**********************************************************************
    !
    ! Function forcing_factor computes the forcing factor
    !
    ! cnres - current residual vector
    ! pp    - type of norm
    !
    ! for pp >= 10 inf norm is computed
    !
    !**********************************************************************

    ! pass input arguments

    REAL(kind=d_t),DIMENSION(:),INTENT(in)   :: cresidual
    INTEGER(kind=li),INTENT(in) :: pp
    REAL(kind=d_t) :: pnorm

    ! length of vector

    INTEGER(kind=li) :: length

    ! local variables

    INTEGER :: i

    ! determine length of cresidual

    length = SIZE(cresidual)

    ! compute pnorm

    pnorm = 0.0_d_t
    IF (pp<10)  THEN
      DO i = 1, length
        pnorm = pnorm + cresidual(i)**(REAL(pp,d_t))
      END DO
      pnorm = pnorm**(1.0_d_t/REAL(pp,d_t))
    ELSE
      pnorm = MAXVAL(ABS(cresidual))
    END IF

  END FUNCTION pnorm

  !-----------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------

  SUBROUTINE norm_eigmode(flux,tot_vol)
    !**********************************************************************
    !
    ! Subroutine norm_eigmode normalizes the eigenmode to 1
    !
    !**********************************************************************

    ! Pass flux

    REAL(kind=d_t) :: flux(num_moments_v,namom,num_cells,egmax,niter)
    REAL(kind=d_t),INTENT(in) :: tot_vol

    ! Local variable

    REAL(kind=d_t) :: s,value
    INTEGER(kind=li) :: i,eg,l,n

    value = REAL(egmax,d_t) * tot_vol

    ! Add up average scalar fluxes

    s = 0.0_d_t
    l = 1_li
    DO eg=1,egmax
      DO i=1,num_cells
        s = s + cells(i)%volume*flux(l,1,i,eg,1)
      END DO
    END DO

    ! Normalize

    DO eg=1,egmax
      DO i=1,num_cells
        DO n=1,namom
          DO l = 1,num_moments_v
            flux(l,n,i,eg,1)=flux(l,n,i,eg,1)*value/s
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE norm_eigmode

  !-----------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------

  SUBROUTINE get_max_error
    !**********************************************************************
    !
    ! Subroutine get_max_error calculates the max residual
    ! per group, it only considers average flux moments
    !
    !**********************************************************************

    ! local arguments

    INTEGER(kind=li) :: i,eg,indx
    REAL(kind=d_t),DIMENSION(num_cells) :: temp

    ! loop over group

    DO eg=1,egmax
      DO i=1,num_cells
        indx=(eg-1)*num_cells*namom*num_moments_v+&
              (i-1)           *namom*num_moments_v+1
        temp(i)=ABS(residual(indx))
      END DO
      max_error(eg)=MAXVAL(temp)
    END DO

  END SUBROUTINE get_max_error

  !-----------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------

  SUBROUTINE single_outer_iteration(flx,keff,reflected_flux,LL,U,Lf,Uf,prnt)
    !**********************************************************************
    !
    ! Subroutine single_outer_iteration performs a single outer
    ! iteration with lagged upscattering. This subroutine is used
    ! for JFNK method 1
    !
    !**********************************************************************

    ! Pass eigenvalue

    REAL(kind=d_t), INTENT(inout) :: keff

    ! Declare scalar flux

    REAL(kind=d_t) :: flx(num_moments_v,namom,num_cells,egmax,niter)

    !  Pre-computed matrices

    REAL(kind=d_t), DIMENSION(num_moments_v,num_moments_v) :: LL, U
    REAL(kind=d_t), DIMENSION(num_moments_f,num_moments_f) :: Lf, Uf

    ! Print flag for inner convergence monitor

    LOGICAL :: prnt

    ! Define fission source

    REAL(kind=d_t) :: fiss_src(num_moments_v,num_cells)

    ! Define source that is passed into inner iteration

    REAL(kind=d_t) :: src(num_moments_v,namom,num_cells,egmax)

    ! local variabels

    INTEGER(kind=li) :: eg,i,l
    INTEGER(kind=li) :: zero,one,two,three
    INTEGER(kind=li) :: mat_indx

    ! dummy for max_errors

    REAL(kind=d_t),DIMENSION(:),ALLOCATABLE :: dmax_error

    ! Define reflected flux

    REAL(kind=d_t),DIMENSION(num_moments_f,grs,8,nangle,egmax) :: reflected_flux

    ! set one,two, three

    one   = 1_li
    two   = 2_li
    three = 3_li
    zero  = 0_li

    ! Allocate dmax_error

    ALLOCATE(dmax_error(egmax))

    ! Initialize fiss_src and src

    src      = 0.0_d_t
    fiss_src = 0.0_d_t

    ! build fission source and initialize src with chi*fiss_src

    DO i = 1, num_cells
      mat_indx=material_ids(reg2mat(cells(i)%reg))
      DO eg = 1, egmax
        DO l = 1, num_moments_v
          fiss_src(l,i) =fiss_src(l,i)+xs_mat(mat_indx)%nu(eg)*xs_mat(mat_indx)%sigma_f(eg) *&
                dens_fact(cells(i)%reg)*flx(l,1,i,eg,1)
        END DO
      END DO
    END DO
    !
    DO i = 1, num_cells
      mat_indx=material_ids(reg2mat(cells(i)%reg))
      DO eg = 1, egmax
        DO l = 1, num_moments_v
          src(l,1,i,eg)=xs_mat(mat_indx)%chi(eg)*fiss_src(l,i)/keff
        END DO
      END DO
    END DO

    !=========================================================================
    ! Compute upscattering
    !=========================================================================
    CALL compute_upscattering(flx, src)

    ! sweep through all groups

    DO eg=1,egmax
      IF  (page_refl.NE.0_li) THEN
        WRITE(stdout_unit, *) "Reflective Boundary values must currently be saved in JFNK."
      END IF

      CALL inner_iteration(eg,flx(:,:,:,eg,2),src(:,:,:,eg),LL,U,Lf,Uf,grs,reflected_flux(:,:,:,:,eg),prnt)

      !=====================================================================
      ! Compute downscattering from eg to egg
      !=====================================================================
      CALL update_downscattering(eg, flx, src)

      ! increase iteration count
      iit_count = iit_count + inner
    END DO

    ! deallocate dmax_error

    DEALLOCATE(dmax_error)

  END SUBROUTINE single_outer_iteration

  !-----------------------------------------------------------------------------------------
  ! THE END
  !-----------------------------------------------------------------------------------------

END MODULE jfnk_module
