!>Outer iteration module contains the subroutines necessary to perform
!>multi-group iteration operations for both eigenvalue and fixed
!>source problems.
!>It calls inner iteration
module outer_iteration_module
  !***********************************************************************
  !
  ! Outer iteration module contains multi-group procedure and calls inner
  ! iteration
  !
  !***********************************************************************

  ! User derived-type modules

  use types
  use parameter_types
  use filename_types
  use vector_types
  use cross_section_types
  use geometry_types
  use angle_types
  use multindex_types
  use global_variables
  use wrapup_module

  ! Use modules that pertain setting up problem

  use termination_module
  use inner_iteration_module
  use dump_inguess_module

  implicit none

contains
  !=============================================================================
  !Subroutine > outer_iteration_ext
  !=============================================================================
  !> This subroutine performs an outer iteration for an external source
  !> problem (compare Alg. 4 in the primer). Note: It does not perform
  !> thermal iterations. The group sweep is explicitly handled by a loop
  !> in the subroutine.
  subroutine outer_iteration_ext(flux,LL,U,Lf,Uf)
    !**********************************************************************
    !
    ! Subroutine outer iteration calls inner iteration and loops over all
    ! energy groups
    !
    !**********************************************************************

    ! Declare scalar flux types used globally

    real(kind=d_t) :: flux(num_moments_v,namom,num_cells,egmax,niter)

    ! Pass pre-computed matrices

    real(kind=d_t), dimension(num_moments_v,num_moments_v) :: LL, U
    real(kind=d_t), dimension(num_moments_f,num_moments_f) :: Lf, Uf


    ! Define temporary variables

    integer(kind=li)               :: alloc_stat, eg, egg, q, octant, i, l, &
         order, n, m, ii, k, indx,face,f
    real(kind=d_t)                 :: t_error
    logical                        :: existence

    ! Define reflected flux

    integer(kind=li)                                 :: rs,rg
    real(kind=d_t),dimension(:,:,:,:,:), allocatable :: reflected_flux

    ! distributed source

    real(kind=d_t) :: src(num_moments_v,namom,num_cells,egmax)

    ! Prepare array reflected_flux
    rs=max(1_li,rside_cells)
    if(page_refl.eq.0_li) then
      rg=egmax
    else
      rg=1_li
    end if
    allocate( reflected_flux(num_moments_f,rs,8,nangle,rg),stat=alloc_stat )
    if(alloc_stat /=0) call stop_thor(2_li)
    reflected_flux=0.0_d_t

    !  Allocate group-dependent maximum spatial error array

    allocate(max_error(egmax),stat=alloc_stat)
    if(alloc_stat /= 0) call stop_thor(2_li)

    ! Open file to page out reflective BC if desired

    if(page_refl.eq.1_li) then

      inquire(file="reflected_flux.pg",exist=existence)
      if( existence .eqv. .true.) then
        open(unit=98,file='reflected_flux.pg',status='replace',form='unformatted',access='direct',recl=64*nangle*rs*num_moments_f)
      else
        open(unit=98,file='reflected_flux.pg',status='new'    ,form='unformatted',access='direct',recl=64*nangle*rs*num_moments_f)
      end if

      do eg = 1,egmax
        write(98,rec=eg) reflected_flux(:,:,:,:,1)
      end do

    end if

    ! Keep track of iterations
    if (rank .eq. 0) then
      write(6,*) '========================================================'
      write(6,*) '   Begin outer iterations.'
      write(6,*) '========================================================'
    end if
    ! Begin outer iteration

    do outer=1, max_outer

      ! Initialize the src with external source ...

      src = zero
      do eg=1,egmax
        do i=1,num_cells
          ! imposed internal source contribution
          do l=1,num_moments_v
            src(l,1,i,eg) = src_str(cells(i)%src,eg)*src_m(l,cells(i)%src,eg)
          end do
        end do
      end do

      ! ... and upscattering

      call compute_upscattering(flux, src)

      do eg=1, egmax

        ! Read binflow for group eg if necessary

        if (page_iflw.eq.1_li) then
          do q=1,nangle
            do octant=1,8
              do f=1,fside_cells
                read(97,*) face
                face = b_cells(face)%ptr
                read(97,*) (binflx(m,face,octant,q,1),m=1,num_moments_f)
              end do
            end do
          end do
        end if

        ! Call inner iteration

        if      (page_refl.eq.0_li) then
          call inner_iteration(eg,flux(:,:,:,eg,niter),src(:,:,:,eg),LL,U,Lf,Uf,rs,reflected_flux(:,:,:,:,eg),.true.)
        else if (page_refl.eq.1_li) then
          read (98,rec=eg) reflected_flux(:,:,:,:,1)
          call inner_iteration(eg,flux(:,:,:,eg,niter),src(:,:,:,eg),LL,U,Lf,Uf,rs,reflected_flux(:,:,:,:,1) ,.true.)
          write(98,rec=eg) reflected_flux(:,:,:,:,1)
        else
          reflected_flux(:,:,:,:,1)=0.0_d_t
          call inner_iteration(eg,flux(:,:,:,eg,niter),src(:,:,:,eg),LL,U,Lf,Uf,rs,reflected_flux(:,:,:,:,1),.true.)
        end if

        ! Compute downscattering from eg to egg

        do egg=eg+1, egmax
          do i=1, num_cells
            do l=0,scatt_ord
              do m=0,l
                indx=1_li+m+(l+1_li)*l/2_li
                do k=1, num_moments_v
                  src(k,indx,i,egg) = src(k,indx,i,egg)                                 +&
                       scat_mult(l,m)*sigma_scat(reg2mat(cells(i)%reg),l+1,egg,eg)%xs    *&
                       dens_fact(cells(i)%reg)*flux(k,indx,i,eg,niter)
                end do
              end do
            end do
            ! odd contributions
            do l=1,scatt_ord
              do m=1,l
                indx=neven+m+(l-1_li)*l/2_li
                do k=1, num_moments_v
                  src(k,indx,i,egg) = src(k,indx,i,egg)                                 +&
                       scat_mult(l,m)*sigma_scat(reg2mat(cells(i)%reg),l+1,egg,eg)%xs    *&
                       dens_fact(cells(i)%reg)*flux(k,indx,i,eg,niter)
                end do
              end do
            end do
          end do
        end do

        ! Rewind unit=97 (finflow) if page_iflw == 1

        if(page_iflw.eq.1_li) rewind(unit=97)
      end do

      ! Compute error ...

      max_outer_error = zero

      do eg =1,egmax
        do i=1,num_cells
          if( abs(flux(1,1,i,eg,niter)) >  1.0e-12_d_t) then
            t_error = abs( flux(1,1,i,eg,niter)-flux(1,1,i,eg,niter-1)) / &
                 flux(1,1,i,eg,niter)
          else
            t_error = abs( flux(1,1,i,eg,niter)-flux(1,1,i,eg,niter-1))
          end if
          if( t_error > max_outer_error ) then
            max_outer_error=t_error
          end if
        end do
      end do

      ! ... and copy over iterates

      do ii=2,niter
        do eg=1, egmax
          do i=1,num_cells
            do n=1,namom
              do l=1,num_moments_v
                flux(l,n,i,eg,ii-1)=flux(l,n,i,eg,ii)
              end do
            end do
          end do
        end do
      end do
      if (rank .eq. 0) then
        write(6,*)   '---------------------------------------'
        write(6,103) '---itn i-itn   max error   max error---'
        write(6,102) outer,tot_nInners, max_outer_error, maxval(max_error),' %% '
        write(6,*)   '---------------------------------------'
        flush(6)
      end if
      if(print_conv.eq.1 .and. rank .eq. 0) then
        write(21,*)   '---------------------------------------'
        write(21,103) '---itn i-itn   max error   max error---'
        write(21,102) outer,tot_nInners, max_outer_error, maxval(max_error),' %% '
        write(21,*)   '---------------------------------------'
        flush(21)
      end if
      103 format(1X,A)
      102 format(1X,2I6,2ES12.4,A)

      ! Convergence check

      if(most_thermal==0) then ! no upscattering
        if(maxval(max_error)<inner_conv) then
          conv_flag=1
          go to 10
        end if
      else
        if(maxval(max_error)<inner_conv .and. max_outer_error<outer_conv) then
          conv_flag=1
          go to 10
        end if
      end if

    end do

    outer=outer-one

    10 continue

    if (rank .eq. 0) then
      write(6,*) '========================================================'
      write(6,*) '   End outer iterations.'
      write(6,*) '========================================================'
    end if
    ! Close reflected flux file if page_ref .eq. 1

    if(page_refl .eq. 1_li) close(unit=98)


    if( allocated(reflected_flux) ) deallocate(reflected_flux)
  end subroutine outer_iteration_ext

  !=============================================================================
  !Subroutine > outer_iteration_eig
  !=============================================================================
  !> This subroutine performs a power iteration for an eigenvalue
  !> problem (compare Alg. 5 in the primer). The subroutine name is a misnomer.
  !> Note: It does not perform
  !> thermal iterations. The group sweep is explicitly handled by a loop
  !> in the subroutine.
  subroutine outer_iteration_eig(flux,keff,LL,U,Lf,Uf)

    ! Pass eigenvalue

    real(kind=d_t), intent(inout) :: keff

    ! Declare angular and scalar flux types used globally

    real(kind=d_t) :: flux(num_moments_v,namom,num_cells,egmax,niter)

    ! Pass pre-computed matrices

    real(kind=d_t), dimension(num_moments_v,num_moments_v) :: LL, U
    real(kind=d_t), dimension(num_moments_f,num_moments_f) :: Lf, Uf

    ! Define temporary variables

    integer(kind=li)               :: alloc_stat, eg, egg, q, octant, i, l, &
         order, ii, n, m, indx, k
    real(kind=d_t)                 :: fiss_den_old, fiss_den_new,       &
         keff_error, keff_old, keff_new, fiss_error, fiss_dist_error(2),&
         flux_error,ts,te
    real(kind=d_t)                 :: t_error
    logical                        :: existence

    ! Define reflected flux

    integer(kind=li)                                 :: rs,rg
    real(kind=d_t),dimension(:,:,:,:,:), allocatable :: reflected_flux

    ! Define fission source

    real(kind=d_t) :: fiss_src(num_moments_v,num_cells,2)

    ! Error mode extrapolation variables

    real(kind=d_t)   :: theta(3),thet
    integer(kind=li) :: extra_flag
    real(kind=d_t)   :: a,b,c

    ! Define source that is passed into inner iteration

    real(kind=d_t) :: src(num_moments_v,namom,num_cells,egmax)

    ! Variables for calling wrapup

    character(100)   :: suffix

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
    integer (kind=li) :: p = 0_li
    integer (kind=li) :: power_iter_count = 0_li, cheby_pi=5_li, cheby_pi_rem
    real(kind=d_t)    :: alpha = zero, beta = zero, gamma = zero
    real(kind=d_t)    :: chebychev_error = zero, entry_error = zero
    real(kind=d_t)    :: entry_theta = zero, theor_err = zero

    integer:: temp_arg
    character:: temp_accel

    ! Prepare array reflected_flux
    rs=max(1_li,rside_cells)
    if(page_refl.eq.0_li) then
      rg=egmax
    else
      rg=1_li
    end if
    allocate( reflected_flux(num_moments_f,rs,8,nangle,rg),stat=alloc_stat )
    if(alloc_stat /=0) call stop_thor(2_li)
    reflected_flux=0.0_d_t

    ! Initialize fiss_src

    fiss_den_old=0.0_d_t
    fiss_den_new=0.0_d_t
    keff=0.0_d_t
    keff_old=1.0_d_t
    keff_new=1.0_d_t

    ! Assume flat source density (for eigenvalue calculation only)

    flux=zero

    do ii=1, niter
      do eg=1, egmax
        do i=1, num_cells
          flux(1,1,i,eg,ii)=one
        end do
      end do
    end do

    fiss_src = zero

    do eg=1, egmax
      do i=1, num_cells
        do l=1, num_moments_v
          fiss_src(l,i,1)=fiss_src(l,i,1)                                         +&
               nu(reg2mat(cells(i)%reg),eg)%xs*fiss(reg2mat(cells(i)%reg),eg)%xs  *&
               dens_fact(cells(i)%reg)*flux(l,1,i,eg,niter)
        end do
      end do
    end do

    !  Allocate group-dependent maximum spatial error array

    allocate(max_error(egmax),stat=alloc_stat)
    if(alloc_stat /= 0) call stop_thor(2_li)

    ! Open file to page out reflective BC if desired

    if(page_refl.eq.1_li) then

      inquire(file="reflected_flux.pg",exist=existence)
      if( existence .eqv. .true.) then
        open(unit=98,file='reflected_flux.pg',status='replace',form='unformatted',access='direct',recl=64*nangle*rs*num_moments_f)
      else
        open(unit=98,file='reflected_flux.pg',status='new'    ,form='unformatted',access='direct',recl=64*nangle*rs*num_moments_f)
      end if

      do eg = 1,egmax
        write(98,rec=eg) reflected_flux(:,:,:,:,1)
      end do

    end if

    ! if inguess_flag == 1

    if (inguess_flag==1) then

      call inguess_eig(flux,keff,niter)
      keff_new=keff
      keff_old=keff
      if (rank .eq. 0) then
        write(6,*)
      end if

      do ii=1,niter-1
        do eg =1,egmax
          do i=1,num_cells
            do n=1,namom
              do l=1,num_moments_v
                flux(l,n,i,eg,ii) = flux(l,n,i,eg,niter)
              end do
            end do
          end do
        end do
      end do

      fiss_src = zero

      do eg=1, egmax
        do i=1, num_cells
          do l=1, num_moments_v
            fiss_src(l,i,1)=fiss_src(l,i,1)                                       +&
                 nu(reg2mat(cells(i)%reg),eg)%xs*fiss(reg2mat(cells(i)%reg),eg)%xs *&
                 dens_fact(cells(i)%reg)*flux(l,1,i,eg,niter)
          end do
        end do
      end do

    end if

    ! Keep track of iterations
    if (rank .eq. 0) then
      write(6,*) '========================================================'
      write(6,*) '   Begin outer iterations.'
      write(6,*) '========================================================'
    end if

    ! Set error mode extrapolation parameters
    theta=0.0_d_t
    extra_flag=0_li

    ! Compute fission density
    fiss_den_new=zero
    do i=1, num_cells
      fiss_den_new=fiss_den_new+cells(i)%volume*fiss_src(1,i,1)
    end do

    !===========================================================================
    ! Begin outer iteration
    !===========================================================================

    temp_arg = command_argument_count()                                         !FIXME - REMOVE AND IMPLEMENT SWITCH INTO INPUT FILE
    if (temp_arg .eq. 2) then
      call get_command_argument(2, temp_accel)
      if ( temp_accel .eq. '1') then
        outer_acc = 1_li
        write(*,*) "NO ACCELERATION"
      else if ( temp_accel .eq. '2') then
        outer_acc = 2_li
        write(*,*) "FISSION SOURCE ACCELERATION"
      else if ( temp_accel .eq. '3') then
        outer_acc = 3_li
        write(*,*) "CHEBYCHEV ACCELERATION"
      else
        write(*,*) "Invalid command line acceleration parameter"
        write(*,*) "Using value from input file"
      end if
    end if                                                                      !FIXME - REMOVE AND IMPLEMENT SWITCH INTO INPUT FILE

    do outer=1, max_outer

      !=========================================================================
      ! Start timer
      !=========================================================================
      call cpu_time(ts)

      !=========================================================================
      ! Compute fission and set all angular moments >0 to zero...
      !=========================================================================
      src = zero
      do eg=1,egmax
        do i=1,num_cells
          do l=1,num_moments_v
            src(l,1,i,eg)  =one/keff_new     *&
                 chi(reg2mat(cells(i)%reg),eg)%xs*fiss_src(l,i,1)
          end do
        end do
      end do

      !=========================================================================
      ! Compute upscattering
      !=========================================================================
      call compute_upscattering(flux, src)

      !=========================================================================
      !Begin Group sweep
      !=========================================================================
      do eg=1, egmax

        !=======================================================================
        ! Call inner iteration
        !=======================================================================
        if      (page_refl.eq.0_li) then
          call inner_iteration(eg,flux(:,:,:,eg,niter),src(:,:,:,eg),LL,U,Lf,Uf,rs,reflected_flux(:,:,:,:,eg),.true.)
        else if  (page_refl.eq.1_li) then
          read (98,rec=eg) reflected_flux(:,:,:,:,1)
          call inner_iteration(eg,flux(:,:,:,eg,niter),src(:,:,:,eg),LL,U,Lf,Uf,rs,reflected_flux(:,:,:,:,1) ,.true.)
          write(98,rec=eg) reflected_flux(:,:,:,:,1)
        else
          reflected_flux(:,:,:,:,1)=0.0_d_t
          call inner_iteration(eg,flux(:,:,:,eg,niter),src(:,:,:,eg),LL,U,Lf,Uf,rs,reflected_flux(:,:,:,:,1) ,.true.)
        end if

        !=====================================================================
        ! Compute downscattering from eg to egg
        !=====================================================================
        call compute_downscattering(eg, flux, src)
      end do

      !========================================================================
      ! Check to see if current cycle meets acceleration criteria
      !========================================================================
      if(outer_acc.eq.2 .and. outer.ge.3) then
        a=abs( ( theta(1) - theta(2) )/theta(1) )
        b=abs( ( theta(2) - theta(3) )/theta(2) )
        if( max(a,b) < extol .and. extra_flag.eq.0_li) then
          extra_flag=1_li
        else
          extra_flag=0_li
        end if
      else if (outer_acc.eq.3 .and. outer.gt.5) then                                       !If more than 5 iterations, accelerate
        if (extra_flag .eq. 0_li .and. theta(3) .gt. .4_d_t .and. theta(3) .lt. 1_li) then !And if spectral radius is large enough
          p=1
          entry_theta = theta(3)                                              !Set sigma hat to previous guess of sp. rad
          extra_flag=1_li
        end if
      end if

      !========================================================================
      ! If yes, perform acceleration
      !========================================================================
      if(extra_flag .eq. 1_li .and. outer_acc.eq.2) then
        ! make sure that the fractional extrapolation is not larger than exmax
        thet=min(exmax/fiss_error,theta(3)/(1.0_d_t-theta(3)))
        ! extrapolation
        do eg=1, egmax
          do i=1,num_cells
            do n=1,namom
              do l=1,num_moments_v
                flux(l,n,i,eg,niter)=flux(l,n,i,eg,niter)+thet*(flux(l,n,i,eg,niter)-flux(l,n,i,eg,niter-1))
              end do
            end do
          end do
        end do
      else if (extra_flag .eq. 1_li .and. outer_acc.eq.3) then
        if (p .eq. 1) then                                                      !First iteration only
          entry_error = zero                                                    !Entry error is the 2 norm of the (n - (n-1)) flux error
          do eg=1, egmax
            do i=1,num_cells
              entry_error=entry_error+ cells(i)%volume*(flux(1,1,i,eg,niter)- flux(1,1,i,eg,niter-1))**2
            end do
          end do
          entry_error = sqrt(entry_error)

          alpha = two/(two-entry_theta)                                         !Hebert Eq. B.5
          beta  = zero                                                          !^
          gamma = zero
          chebychev_error = zero
          theor_err = zero
        else                                                                    !For all iteration p!=1

          gamma = dacosh((two/entry_theta) - one)                               !Hebert Eq. B.5
          alpha = (four/entry_theta) * dcosh((p-one)*gamma)/dcosh(p*gamma)      !Hebert Eq. B.5
          beta  = (one-(entry_theta/two))-(one/alpha)                           !Hebert Eq. B.5

        end if

        chebychev_error = zero                                                  !Chebychev_err is the two norm of (n - (n-1)) flux error
        do eg=1, egmax                                                          !divided by the enty error (error === 1 for p=1)
          do i=1,num_cells
            chebychev_error=chebychev_error+cells(i)%volume  *&
                 (flux(1,1,i,eg,niter)- flux(1,1,i,eg,niter-1))**2
          end do
        end do
        chebychev_error = sqrt(chebychev_error)
        chebychev_error = chebychev_error/entry_error


        do eg=1, egmax                                                          !Accelerate flux
          do i=1,num_cells
            do n=1,namom
              do l=1,num_moments_v
                flux(l,n,i,eg,niter)=flux(l,n,i,eg,niter-1)             + &
                     alpha*(flux(l,n,i,eg,niter)-flux(l,n,i,eg,niter-1))+ &
                     beta*(flux(l,n,i,eg,niter-1)-flux(l,n,i,eg,niter-2))
              end do
            end do
          end do
        end do


        theor_err = (dcosh((real(p,d_t)-one)*dacosh((two - entry_theta) /&      !Theretical error
             entry_theta)))**(-one)

        if(chebychev_error > theor_err .and. &
             p .ge. 3*(power_iter_count+3)) then                                !If insufficient decrease in error,
          !write(*,*) "RESTART FLAG THROWN"                                      !reset chebychev
          extra_flag = 0_li
          !write(*,*) "BEGIN SPECTRAL RADIUS BACKTRACK"
          extra_flag = -1_li                                                    !Set flag to acceleration interrupt
          cheby_pi_rem = cheby_pi
        end if

        p=p+1                                                                   !Increment p
      else if(extra_flag .eq. -1_li) then                                       !If acceleration is interrupted
        !write(*,*) "Performing Between-Chebychev-Cycle Power Iteration"
        power_iter_count = power_iter_count + 1_li
        cheby_pi_rem = cheby_pi_rem -1_li                                       !Track remaining interrupts
        p=0
        if (cheby_pi_rem .eq. 0) extra_flag = 0_li                              !Reactivate acceleration
      end if

      !========================================================================
      ! Zero fission source. Then, recompute fission source using new update
      !========================================================================
      do i=1, num_cells
        do l=1, num_moments_v
          fiss_src(l,i,2)=zero
        end do
      end do
      do eg=1, egmax
        do i=1, num_cells
          do l=1, num_moments_v
            fiss_src(l,i,2)=fiss_src(l,i,2)                                        +&
                 nu(reg2mat(cells(i)%reg),eg)%xs*fiss(reg2mat(cells(i)%reg),eg)%xs *&
                 dens_fact(cells(i)%reg)*flux(l,1,i,eg,niter)
          end do
        end do
      end do

      !========================================================================
      ! Compute new fission distribution and density
      !========================================================================
      fiss_den_old=fiss_den_new
      fiss_den_new=zero
      do i=1, num_cells
        fiss_den_new=fiss_den_new+cells(i)%volume*fiss_src(1,i,2)
      end do

      !========================================================================
      ! Compute new eigenvalue based on fission densities
      !========================================================================
      keff_new=keff_old*(fiss_den_new/fiss_den_old)

      !========================================================================
      ! compute error in keff and copy new to old iterate
      !========================================================================
      keff_error=abs(keff_new-keff_old)/keff_new
      keff_old=keff_new

      !========================================================================
      ! compute convergence based on fission source
      !========================================================================
      fiss_dist_error(1)=fiss_dist_error(2)
      fiss_dist_error(2)=0.0_d_t
      fiss_error        =0.0_d_t

      do i=1, num_cells
        fiss_dist_error(2)=fiss_dist_error(2)+cells(i)%volume*    &
             abs( fiss_src(1,i,2) - fiss_src(1,i,1) )
        if(abs(fiss_src(1,i,2)) > 1e-12)then
          t_error=abs(( fiss_src(1,i,2) - fiss_src(1,i,1) )/ fiss_src(1,i,2))
        else
          t_error=abs(  fiss_src(1,i,2) - fiss_src(1,i,1) )
        end if
        if(t_error > fiss_error)then
          fiss_error=t_error
        end if
      end do

      !========================================================================
      ! compute the convergence based on error in group fluxes
      !========================================================================
      flux_error=0.0_d_t
      do eg =1,egmax
        do i=1,num_cells
          if( abs(flux(1,1,i,eg,niter)) > 1.0e-12_d_t) then
            t_error = abs( flux(1,1,i,eg,niter)-flux(1,1,i,eg,niter-1))/&
                 flux(1,1,i,eg,niter)
          else
            t_error = abs( flux(1,1,i,eg,niter)-flux(1,1,i,eg,niter-1))
          end if
          if( t_error > flux_error ) then
            flux_error=t_error
          end if
        end do
      end do

      !========================================================================
      ! compute theta (Spectral Radius)
      !========================================================================
      if (p.eq.0) theta(1)=theta(2)
      if (p.eq.0) theta(2)=theta(3)
      if (p.eq.0) theta(3)=fiss_dist_error(2)/fiss_dist_error(1)

      !========================================================================
      ! stop timer
      !========================================================================
      call cpu_time(te)

      !========================================================================
      ! write information for outer iteration
      !========================================================================
      if (rank .eq. 0) then
        write(6,*)   '---------------------------------------------------------------------------------------------------'
        write(6,101) '---itn i-itn        keff    err-keff    err-fiss     err-flx      Sp Rad      extrap        time---'
        write(6,102) outer,tot_nInners ,keff_new, keff_error,fiss_error,flux_error,theta(3),extra_flag,te-ts,' %% '
        write(6,*)   '---------------------------------------------------------------------------------------------------'
        flush(6)
        if(print_conv.eq.1) then
          write(21,*)   '---------------------------------------------------------------------------------------------------'
          write(21,101) '---itn i-itn        keff    err-keff    err-fiss     err-flx      Sp Rad      extrap        time---'
          write(21,102) outer,tot_nInners ,keff_new, keff_error,fiss_error,flux_error,theta(3),extra_flag,te-ts,' %% '
          write(21,*)   '---------------------------------------------------------------------------------------------------'
          flush(21)
        end if
      end if
      101 format(1X,A)
      102 format(1X,2I6,5ES12.4,I12,ES12.4,A)
      103 format(1X,A,ES12.4)
      104 format(1X,A,I5,3ES12.4)

      !========================================================================
      ! write iteration results to file if desired
      !========================================================================
      if(dump_flag==1) then
        call dump_PI(flux,keff_old)
      end if

      !========================================================================
      ! call wrapup even/odd
      !========================================================================
      if (mod(outer,2) .eq. 0) then
        suffix = "even"
        open(unit=49,file='intermediate_output_even.dat',status='unknown',action='write')
      else
        suffix = "odd"
        open(unit=49,file='intermediate_output_odd.dat',status='unknown',action='write')
      end if
      call wrapup(flux = flux ,keff = keff_new, unit_number = 49, suffix = suffix, is_final = .false.)
      close(unit=49)

      !========================================================================
      ! check convergence and quit if criteria are satisfied
      !========================================================================
      if(outer > 2 .and. (fiss_error < outer_conv .or. flux_error < outer_conv) .and. &
           keff_error < k_conv) then
        conv_flag=1
        go to 10
      end if

      !========================================================================
      ! Copy over iterates of flux
      !========================================================================
      do ii=2,niter
        do eg=1, egmax
          do i=1,num_cells
            do n=1,namom
              do l=1,num_moments_v
                flux(l,n,i,eg,ii-1)=flux(l,n,i,eg,ii)
              end do
            end do
          end do
        end do
      end do

      !========================================================================
      ! copy over the fission source: Note, that if the flux is accelerated the
      ! newly computed fission source will also be accelerated
      !========================================================================
      do i=1, num_cells
        do l=1, num_moments_v
          fiss_src(l,i,1)=fiss_src(l,i,2)
        end do
      end do

    end do

    outer=outer-1

    10 continue

    k_error=keff_error
    f_error=fiss_error
    max_outer_error = flux_error
    keff=keff_old

    if (rank .eq. 0) then
      write(6,*) '========================================================'
      write(6,*) '   End outer iterations.'
      write(6,*) '========================================================'
    end if

    ! Close reflected flux file if page_ref .eq. 1

    if(page_refl .eq. 1_li) close(unit=98)

    ! Deallocate temporary arrays

    if( allocated(reflected_flux) ) deallocate(reflected_flux)
  end subroutine outer_iteration_eig

  !=============================================================================
  !Subroutine > compute_upscattering
  !=============================================================================
  !> Computes the iteration source upscattering component based on the passed
  !> flux and  various global scattering & materials parameters.
  !> This subroutine is for both the eigenvalue and fixed source solvers
  subroutine compute_upscattering(flux, src)

    integer:: eg, egg, i, l, m, k, indx
    real(kind=d_t) :: flux(num_moments_v,namom,num_cells,egmax,niter)
    real(kind=d_t) :: src(num_moments_v,namom,num_cells,egmax)

    do eg=most_thermal,egmax     ! eg is the group that it is scattered to
      do i=1,num_cells
        ! even contributions
        do l=0,scatt_ord
          do m=0,l
            indx=1_li+m+(l+1_li)*l/2_li
            do egg=eg+1,egmax  ! egg is the group that is scattered from
              do k=1, num_moments_v
                src(k,indx,i,eg) =  src(k,indx,i,eg) +&
                     scat_mult(l,m)*sigma_scat(reg2mat(cells(i)%reg),l+1,eg,egg)%xs  *  &
                     dens_fact(cells(i)%reg)*flux(k,indx,i,egg,niter)
              end do
            end do
          end do
        end do
        ! odd contributions
        do l=1,scatt_ord
          do m=1,l
            indx=neven+m+(l-1_li)*l/2_li
            do egg=eg+1,egmax  ! egg is the group that is scattered from
              do k=1, num_moments_v
                src(k,indx,i,eg) = src(k,indx,i,eg)                                   +&
                     scat_mult(l,m)*sigma_scat(reg2mat(cells(i)%reg),l+1,eg,egg)%xs    *&
                     dens_fact(cells(i)%reg)*flux(k,indx,i,egg,niter)
              end do
            end do
          end do
        end do
      end do
    end do
  end subroutine compute_upscattering

  !=============================================================================
  !Subroutine > compute_downscattering
  !=============================================================================
  !> Computes the iteration source downscattering component based on the passed
  !> flux and  various global scattering & materials parameters.
  !> This subroutine is for both the eigenvalue and fixed source solvers
  subroutine compute_downscattering(eg, flux, src)

    integer:: eg, egg, i, l, m, k, indx
    real(kind=d_t) :: flux(num_moments_v,namom,num_cells,egmax,niter)
    real(kind=d_t) :: src(num_moments_v,namom,num_cells,egmax)

    do egg=eg+1, egmax
      do i=1, num_cells
        ! Even contributions
        do l=0,scatt_ord
          do m=0,l
            indx=1_li+m+(l+1_li)*l/2_li
            do k=1, num_moments_v
              src(k,indx,i,egg) = src(k,indx,i,egg)                                +&
                   scat_mult(l,m)*sigma_scat(reg2mat(cells(i)%reg),l+1,egg,eg)%xs   *&
                   dens_fact(cells(i)%reg)*flux(k,indx,i,eg,niter)
            end do
          end do
        end do
        ! odd contributions
        do l=1,scatt_ord
          do m=1,l
            indx=neven+m+(l-1_li)*l/2_li
            do k=1, num_moments_v
              src(k,indx,i,egg) = src(k,indx,i,egg)                                 +&
                   scat_mult(l,m)*sigma_scat(reg2mat(cells(i)%reg),l+1,egg,eg)%xs    *&
                   dens_fact(cells(i)%reg)*flux(k,indx,i,eg,niter)
            end do
          end do
        end do
      end do
    end do
  end subroutine compute_downscattering
end module outer_iteration_module
!-----------------------------------------------------------------------------------------
! End
!-----------------------------------------------------------------------------------------
