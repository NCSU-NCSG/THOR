module jfnk_module
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

  use outer_iteration_module
  use inner_iteration_module
  use dump_inguess_module

  implicit none

  ! Variables specific to jfnk: integer

  integer(kind=li) :: method   ! flag to switch between 3 formulations
                               ! of F(u): 1 - lagged upscattering
                               !          2 - flat
                               !          3 - flat with grp sweep
  integer(kind=li) :: max_rest ! max # restarts of gmres/krylov
  integer(kind=li) :: restart  ! max # iterations between restarts
  integer(kind=li) :: max_kit  ! max # of krylov iterations - # matrix*vector
  integer(kind=li) :: num_var  ! number of variables
  integer(kind=li) :: iit_count! inner iteration count
  integer(kind=li) :: dflag    ! determines if dump file is written
  integer(kind=li) :: igflag   ! determines if initial guess is read
  integer(kind=li) :: max_nit  ! maximum Newton iterations
  integer(kind=li) :: grs      ! =min(#reflective faces, 1)

  ! Variables specific to jfnk: real

  real(kind=d_t),dimension(:),allocatable :: residual ! residual of most recent newton step
  real(kind=d_t),dimension(:),allocatable :: du       ! newton step, solution increment
  real(kind=d_t) :: jac_eps                           ! finite difference perturbation
  real(kind=d_t) :: nres                              ! newton residual nres=||residual||_p
  real(kind=d_t) :: nit_conv                          ! newton iteration convergence


contains

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

  subroutine set_jfnk
  !**********************************************************************
  !
  ! Subroutine jfnk_setup allocates and sets the jfnk variables
  !
  !**********************************************************************

  ! local variables

    integer(kind=li) :: alloc_stat

  ! set method, krest and kit

    method = rd_method
    restart = rd_restart
    max_kit = rd_max_kit
    max_rest = int(max_kit/restart)

  ! calculate length of unknown vector

    num_var = egmax*num_cells*num_moments_v*namom

  ! allocate residual and du

    if(.not.allocated(du)) allocate(du(num_var+1),stat=alloc_stat)
    if(alloc_stat /= 0) call stop_thor(2_li)
    if(.not.allocated(residual)) allocate(residual(num_var+1),stat=alloc_stat)
    if(alloc_stat /= 0) call stop_thor(2_li)

  ! set dflag, igflag

    dflag =dump_flag
    igflag=inguess_flag

  ! set jac_eps

    jac_eps = 1.4901d-08

  ! set max_nit and nit_conv

    max_nit=max_outer
    nit_conv=outer_conv

  ! make sure that niter .eq. 2

    if (niter .ne. 2) call stop_thor(25_li)

  ! Allocate max_error

    allocate(max_error(egmax),stat=alloc_stat)
    if(alloc_stat /= 0) call stop_thor(2_li)

  ! Set grs

    grs=max(1_li,rside_cells)

  end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

  subroutine clean_jfnk
  !**********************************************************************
  !
  ! subroutine clean_jfnk deallocates variables used in jfnk algorithm
  !
  !**********************************************************************

  ! deallocate residual and du

    if(allocated(du)) deallocate(du)
    if(allocated(residual)) deallocate(residual)

  end subroutine

  subroutine do_jfnk(flux,keff,LL,U,Lf,Uf)
  !**********************************************************************
  !
  ! Subroutine do_jfnk is the driver routine for the jfnk method
  !
  !**********************************************************************

  ! Declare scalar flux types used globally+keff

    real(kind=d_t) :: flux(num_moments_v,namom,num_cells,egmax,niter)
    real(kind=d_t) :: keff

  ! Define temporary variables

    real(kind=d_t), dimension(num_moments_v,num_moments_v) :: LL, U
    real(kind=d_t), dimension(num_moments_f,num_moments_f) :: Lf, Uf

  ! Local iteration counter

    integer(kind=li) :: kit

  ! Local variables

    real(kind=d_t) :: kerr,tol,tfiss1,tfiss2
    integer(kind=li) :: i,j,l,indx,eg,alloc_stat,m,n,k
    integer(kind=li) :: one,two,three
    real(kind=d_t) :: tot_vol
    real(kind=d_t) :: ts,te,keff_error,max_outer_error,tmp

  ! Reflective flux in case ipiter > 0

    real(kind=d_t), allocatable, dimension(:,:,:,:,:) :: reflected_flux

  ! set one two three

    one = 1_li
    two = 2_li
    three = 3_li

  ! Set cell volumes
    tot_vol = 0.0_d_t
    do i=1, num_cells
       tot_vol=tot_vol+cells(i)%volume
    end do

  ! start newton iteration

    nit       = 0
    tot_kit   = 0
    iit_count = 0

  ! set up initial guess: perform initial power iteration ??

    if(igflag==1) then
      call inguess_eig(flux,keff,one)
    else
      flux=0.0_d_t
      do eg =1,egmax
        do i=1,num_cells
          flux(1,1,i,eg,1)=1.0_d_t
          flux(1,1,i,eg,2)=1.0_d_t
        end do
      end do
      keff = 1.0_d_t
    end if
  ! Perform ipow initial power iterations
    if (ipow .gt. 0 .and. rank .eq. 0) then
      write(6,*)
      write(6,*) '========================================================'
      write(6,*) '   Begin Initial Power Iterations.'
      write(6,*) '========================================================'
    end if

    ! note: we do not need to duplicate reflected_flux on all processors
    ! but for simplicity we do here.
    if (ipow .gt. 0) then
      allocate(reflected_flux(num_moments_f,grs,8,nangle,egmax))
      reflected_flux = 0.0_d_t
    end if

    do m = 1,ipow

    ! start timer

      call cpu_time(ts)

    ! call single outer iteration - it does not update keff

      call single_outer_iteration(flux,keff,reflected_flux,LL,U,Lf,Uf,.true.)

    ! update keff

      tfiss1 = tfissions (flux,two)
      tfiss2 = tfissions (flux,one)
      keff_error = abs(keff * tfiss1 / tfiss2-keff)/keff
      keff = keff * tfiss1 / tfiss2

    ! copy iterate #2 into working iterate #1
      max_outer_error=0.0_d_t
      do eg =1,egmax
        do i=1,num_cells
          tmp=abs(flux(1,1,i,eg,2)-flux(1,1,i,eg,1)) / &
                  flux(1,1,i,eg,1)
          if(tmp>max_outer_error) max_outer_error=tmp
          do n =1,namom
            do l=1,num_moments_v
              flux(l,n,i,eg,1)=flux(l,n,i,eg,2)
            end do
          end do
        end do
      end do

    ! stop timer

      call cpu_time(te)

    ! write convergence monitor
      if (rank .eq. 0) then
        write(6,*)   '---------------------------------------------------------------------------'
        write(6,401) '---itn        keff    err-keff     err-flx        time---'
        write(6,402) m,keff, keff_error,max_outer_error,te-ts,' %i'
        write(6,*)   '---------------------------------------------------------------------------'
        if(print_conv.eq.1) then
          write(21,*)   '---------------------------------------------------------------------------'
          write(21,401) '---itn        keff    err-keff     err-flx        time---'
          write(21,402) m,keff, keff_error,max_outer_error,te-ts,' %i'
          write(21,*)   '---------------------------------------------------------------------------'
        end if
      end if
      401 FORMAT(1X,A)
      402 FORMAT(1X,I6,4ES12.4,A)

    end do

    if( allocated(reflected_flux) ) deallocate(reflected_flux)

    if (ipow .gt. 0 .and. rank .eq. 0) then
      write(6,*)
      write(6,*) '========================================================'
      write(6,*) '   End Initial Power Iterations.'
      write(6,*) '========================================================'
    end if

    if (rank .eq. 0) then
      write(6,*)
      write(6,*) '========================================================'
      write(6,*) '   Begin JFNK iterations.'
      write(6,*) '========================================================'
    end if

  ! evaluate residual

    residual=evaluate_residual(flux,keff,LL,U,Lf,Uf)

  ! compute current p-norm of newton residual

    nres = pnorm(residual,12_li)

  ! compute forcing factor

    tol = forcing_factor(nres,nit)

  ! increment newton iteration counter

    nit = nit + 1

    do i = 1, max_nit

      ! Start Timer

        call cpu_time(ts)

      ! call krylov solver

        call krylovM(-residual,du,tol,kit,kerr,&
               ! overhead variables
               flux,keff,LL,U,Lf,Uf)

      ! update unknown

        call newton_update(flux,keff)

      ! dump result into file if dump_flag==1

        if(dflag==1) then
          call dump_jfnk(flux,keff)
        end if

      ! evaluate residual

        residual=evaluate_residual(flux,keff,LL,U,Lf,Uf)

      ! compute current p-norm of newton residual

        nres = pnorm(residual,12_li)

      ! compute forcing factor

        tol = forcing_factor(nres,nit)

      ! stop timer

        call cpu_time(te)

      ! output convergence progress
        if (rank .eq. 0) then
          write(6,*)   '---------------------------------------------------------------------------'
          write(6,102) '---nitn  kitn        keff     max-res        time---'
          write(6,101) nit,kit,keff,nres,te-ts," %n"
          write(6,*)   '---------------------------------------------------------------------------'
          if(print_conv.eq.1) then
            write(21,*)   '---------------------------------------------------------------------------'
            write(21,102) '---nitn  kitn        keff     max-res        time---'
            write(21,101) nit,kit,keff,nres,te-ts," %n"
            write(21,*)   '---------------------------------------------------------------------------'
          end if
        end if
        101 format (1X,I7,I6,3ES12.4,A)
        102 format (1X,A)

      ! Check for convergence

      if(nres<nit_conv) then
          write(6,*)   '---------------------------------------------------------------------------'
          write(6,104) 'Newton Iteration Convergence achieved after',nit,' iterations.'
          write(6,103) 'with final residual ',nres
          write(6,*)   '---------------------------------------------------------------------------'
          conv_flag=1
          go to 10
        end if
        103 format (1X,A,ES12.4)
        104 format (1X,A44,I4,A)

      ! increment iteration counters

        nit = nit + 1
        tot_kit = tot_kit + kit


    end do

   nit=nit-1

10 continue

    ! normalize eigenmode

      call norm_eigmode(flux,tot_vol)

      write(6,*) '========================================================'
      write(6,*) '   End JFNK iterations.'
      write(6,*) '========================================================'

    ! newton iteration end, copy unknown to flux

      do eg = 1, egmax
        do i = 1, num_cells
          do n=1,namom
            do l = 1, num_moments_v
               flux(l,n,i,eg,2)=flux(l,n,i,eg,1)
            end do
          end do
        end do
      end do

    ! assign total # inner iterations

      tot_nInners = iit_count

    ! calculate max_error

      call get_max_error


  end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

  function evaluate_residual(flux,keff,LL,U,Lf,Uf)
  !**********************************************************************
  !
  ! Function evaluate_residual evaluates the non-linear function
  ! r = F(u)
  !
  !**********************************************************************

    ! Pass flux

      real(kind=d_t) :: flux(num_moments_v,namom,num_cells,egmax,niter)

    ! pass multiplication factor

      real(kind=d_t) :: keff

    ! Precomputed matrices

      real(kind=d_t), dimension(num_moments_v,num_moments_v) :: LL, U
      real(kind=d_t), dimension(num_moments_f,num_moments_f) :: Lf, Uf

    ! define return value

      real(kind=d_t), dimension(num_var+1) :: evaluate_residual

    ! Local variables

      integer(kind=li) :: eg,i,j,l,alloc_stat,m,n,indx,zero,one,two,three,&
                          q,octant
      real(kind=d_t)   :: f1,f2

    ! Source: src is a general source used for all methods

      real(kind=d_t) :: src(num_moments_v,namom,num_cells,egmax)

    ! Reflected flux array

      real(kind=d_t),allocatable ::  reflected_flux(:,:,:,:,:)

    ! if method = 2 or 3 allocate and set outward_normal. allocate src

      if (method==1) then

      ! initialize reflected_flux

        allocate(reflected_flux(num_moments_f,grs,8,nangle,egmax),stat=alloc_stat)
        if(alloc_stat /=0) call stop_thor(2_li)
        reflected_flux = 0.0_d_t

      else if(method==2 .or. method==3) then

      ! initialize reflected_flux

        allocate(reflected_flux(num_moments_f,grs,8,nangle,1),stat=alloc_stat)
        if(alloc_stat /=0) call stop_thor(2_li)
        reflected_flux = 0.0_d_t

      ! initialize source

        src = 0.0_d_t

      end if

    ! set 1,2,3

      zero = 0_li
      one = 1_li
      two = 2_li
      three = 3_li

    ! initialize new iterates

      if(method .eq. 1) then
         do eg = 1, egmax
           do i = 1, num_cells
             do n = 1, namom
               do l = 1, num_moments_v
                 flux(l,n,i,eg,2)=flux(l,n,i,eg,1)
               end do
             end do
           end do
         end do
      else if (method .eq. 2 .or. method .eq. 3) then
         do eg = 1, egmax
           do i = 1, num_cells
             do n = 1, namom
               do l = 1, num_moments_v
                 flux(l,n,i,eg,2) = 0.0_d_t
               end do
             end do
           end do
         end do
      end if


    ! evaluate residual based on the chosen method

      if (method==1) then

        ! Outer iteration with lagged upscattering

          call single_outer_iteration(flux,keff,reflected_flux,LL,U,Lf,Uf,.false.)

      else if (method==2) then

        ! Flat - Sweep on fixed source
        ! 1. Build fixed source

          call bfsrc (flux,src,keff,one)

          do eg = 1, egmax
            call bscsrc (flux,src,one,one,eg)
          end do

        ! 2. Sweep on source
          do eg = 1, egmax
            reflected_flux=0.0_d_t
            call sweep(eg,flux(:,:,:,eg,2),src(:,:,:,eg),grs,reflected_flux,LL,U,Lf,Uf)
          end do

      else if (method==3) then

        ! Sweep through groups with updated downscattering
        ! 1. Build fission source

          call bfsrc (flux,src,keff,one)

        ! Calculate self-scattering source into group 1
          eg = 1
          call bscsrc (flux,src,one,one,eg)

        ! Sweep on fixed source in group 1
          reflected_flux=0.0_d_t
          call sweep (eg,flux(:,:,:,eg,2),src(:,:,:,eg),grs,reflected_flux,LL,U,Lf,Uf)

        ! 2. Sweep through grps 2..egmax and re-evaluate downscatter source
          do eg = 2, egmax
            call bscsrc (flux,src,two,one,eg)
        !
            reflected_flux=0.0_d_t
            call sweep(eg,flux(:,:,:,eg,2),src(:,:,:,eg),grs,reflected_flux,LL,U,Lf,Uf)
          end do
      else
        call stop_thor(26_li)
      end if

    ! Compute evaluate_residual

        do eg = 1, egmax
          do i = 1, num_cells
             do n=1,namom
               do l = 1, num_moments_v
                 indx = (eg-1)*num_cells*namom*num_moments_v+&
                        (i-1)           *namom*num_moments_v+&
                        (n-1)                 *num_moments_v+&
                        l
                 evaluate_residual(indx) =  flux(l,n,i,eg,1)-flux(l,n,i,eg,2)
               end do
             end do
          end do
        end do

      ! assign [num_var+1] using FR approach

        f1 = tfissions (flux,one)
        f2 = tfissions (flux,two)

        evaluate_residual(num_var+1) = keff*(1.0_d_t-f2/f1)

    ! Deallocate

      if(method==2 .or. method==3) then

        ! deallocate reflected flux
        deallocate(reflected_flux)

      end if

  end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

  subroutine bfsrc (flx,src,keff,flag)
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

      integer(kind=li), intent(in) :: flag

    ! Pass flux and source and keff

      real(kind=d_t) :: flx (num_moments_v,namom,num_cells,egmax,niter)
      real(kind=d_t) :: src (num_moments_v,namom,num_cells,egmax)
      real(kind=d_t) :: keff

    ! Define fission source

     real(kind=d_t) :: tfsrc(num_moments_v,num_cells)

    ! Local variables

      integer(kind=li) :: alloc_stat,eg,m,n,i,l,egg,k,indx

    ! Initialize tfsrc

      tfsrc = 0.0_d_t

    ! Compute total fission source

       do eg = 1, egmax
         do i = 1, num_cells
           do l = 1, num_moments_v
             tfsrc(l,i) = tfsrc(l,i) + nu(reg2mat(cells(i)%reg),eg)%xs*fiss(reg2mat(cells(i)%reg),eg)%xs * &
                                       dens_fact(cells(i)%reg)*flx(l,1,i,eg,flag)
           end do
         end do
       end do

    ! Add to src

      do eg = 1, egmax
        do i = 1, num_cells
          do l = 1, num_moments_v
            src(l,1,i,eg) = src(l,1,i,eg) + chi(reg2mat(cells(i)%reg),eg)%xs*tfsrc(l,i)/keff
          end do
        end do
      end do

    ! Add upscattering

      do eg=most_thermal,egmax    ! eg is the group that it is scattered to
        do i=1,num_cells
          ! even contribution
          do l=0,scatt_ord
            do m=0,l
              indx=1_li+m+(l+1_li)*l/2_li
              do egg=eg+1,egmax  ! egg is the group that is scattered from
                do k=1, num_moments_v
                  src(k,indx,i,eg) = src(k,indx,i,eg)                                  +&
                      scat_mult(l,m)*sigma_scat(reg2mat(cells(i)%reg),l+1,eg,egg)%xs   *&
                      dens_fact(cells(i)%reg)*flx(k,indx,i,egg,flag)
                end do
              end do
            end do
          end do
          ! odd contribution
          do l=1,scatt_ord
            do m=1,l
              indx=neven+m+(l-1_li)*l/2_li
              do egg=eg+1,egmax  ! egg is the group that is scattered from
                do k=1, num_moments_v
                  src(k,indx,i,eg)=src(k,indx,i,eg)                                    +&
                      scat_mult(l,m)*sigma_scat(reg2mat(cells(i)%reg),l+1,eg,egg)%xs   *&
                      dens_fact(cells(i)%reg)*flx(k,indx,i,egg,flag)
                end do
              end do
            end do
          end do
        end do
     end do

  end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

  subroutine bscsrc (flx,src,flag1,flag2,egg)
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

      integer(kind=li), intent(in) :: flag1,flag2

    ! Pass flux and source and keff

      real(kind=d_t) :: flx (num_moments_v,namom,num_cells,egmax,niter)
      real(kind=d_t) :: src (num_moments_v,namom,num_cells,egmax)

    ! Pass group

      integer(kind=li) :: egg

    ! Local variables

      integer(kind=li) :: eg,i,l,m,n,k,indx

    ! Start adding downscattering source

      do eg = 1, egg-1
        do i = 1, num_cells
          ! even contribution
          do l=0, scatt_ord
            do m = 0,l
              indx=1_li+m+(l+1_li)*l/2_li
              do k = 1, num_moments_v
                src(k,indx,i,egg)=src(k,indx,i,egg)                               +&
                  scat_mult(l,m)*sigma_scat(reg2mat(cells(i)%reg),l+1,egg,eg)%xs  *&
                  dens_fact(cells(i)%reg)*flx(k,indx,i,eg,flag1)
              end do
            end do
          end do
          ! odd contribution
          do l=1, scatt_ord
            do m = 1,l
              indx=neven+m+(l-1_li)*l/2_li
              do k = 1, num_moments_v
                src(k,indx,i,egg)=src(k,indx,i,egg)                               +&
                  scat_mult(l,m)*sigma_scat(reg2mat(cells(i)%reg),l+1,egg,eg)%xs  *&
                  dens_fact(cells(i)%reg)*flx(k,indx,i,eg,flag1)
              end do
            end do
          end do
        end do
      end do

    ! Add selfscattering source

      do i = 1, num_cells
        ! Even contributions
        do l = 0, scatt_ord
          do m = 0,l
            indx=1_li+m+(l+1_li)*l/2_li
            do k = 1, num_moments_v
              src(k,indx,i,egg)=src(k,indx,i,egg)                                +&
                scat_mult(l,m)*sigma_scat(reg2mat(cells(i)%reg),l+1,egg,egg)%xs  *&
                dens_fact(cells(i)%reg)*flx(k,indx,i,eg,flag2)
            end do
          end do
        end do
        ! odd contributions
        do l = 1, scatt_ord
          do m = 1,l
            indx=neven+m+(l-1_li)*l/2_li
            do k = 1, num_moments_v
              src(k,indx,i,egg)=src(k,indx,i,egg)                                +&
                scat_mult(l,m)*sigma_scat(reg2mat(cells(i)%reg),l+1,egg,egg)%xs  *&
                dens_fact(cells(i)%reg)*flx(k,indx,i,eg,flag2)
            end do
          end do
        end do
      end do

  end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

  function tfissions (flux,flag)
   !**********************************************************************
   !
   ! function tfission computes the total fissions in the domain
   !
   ! Sum_{#groups,#cells} nsigf(grp,cell#)*phi_000
   !
   !**********************************************************************

     ! Pass input parameters

       integer(kind=li), intent(in) :: flag

     ! pass flux

       real(kind=d_t) :: flux(num_moments_v,namom,num_cells,egmax,niter)

     ! define return value

       real(kind=d_t) :: tfissions

     ! local variables

       integer(kind=li) :: eg,i,indx

     ! initialize tfissions

       tfissions = 0.0_d_t

     ! compute tfission

       do eg = 1, egmax
         do i = 1, num_cells
           tfissions = tfissions+ &
                     nu(reg2mat(cells(i)%reg),eg)%xs*fiss(reg2mat(cells(i)%reg),eg)%xs* &
                     cells(i)%volume*dens_fact(cells(i)%reg)*flux(1,1,i,eg,flag)
         end do
       end do

  end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

  subroutine krylovM(rhs,solution,tol,kit,error,          &
       ! overhead variables
       flux,keff,LL,U,Lf,Uf)
  !**********************************************************************
  !
  ! Subroutine krylovM is the interface to sparskit
  !
  !**********************************************************************

    ! Pass input variables

      real(kind=d_t),dimension(num_var+1) :: rhs,solution
      real(kind=d_t) :: tol
      integer(kind=li) :: kit
      real(kind=d_t) :: error

    ! Pass flux

      real(kind=d_t) :: flux(num_moments_v,namom,num_cells,egmax,niter)

    ! pass multiplication factor

      real(kind=d_t) :: keff

    ! Pre-computed matrices

      real(kind=d_t), dimension(num_moments_v,num_moments_v) :: LL, U
      real(kind=d_t), dimension(num_moments_f,num_moments_f) :: Lf, Uf

    ! Local variables
    ! ws    - work-size, dimension of work array
    ! ipar  - integer parameter krylov method
    ! work  - work array
    ! fpar  - real parameters krylov method

      integer(kind=li) :: ws
      integer,dimension(16) :: ipar
      real(kind=d_t),dimension(:),allocatable :: work
      real(kind=d_t),dimension(16) :: fpar
      integer :: alloc_stat
      real(kind=d_t) :: ts,te,nrhs
      integer :: k

    ! Calculate ws(for gmres)

      ws = (num_var+4)*(restart+2)+(restart+1)*restart/2

    ! Allocate and initialize  work

      allocate(work(ws),stat=alloc_stat)
      if(alloc_stat /= 0) call stop_thor(2_li)
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

      write(6,102) '---itn         err    trgt-res        time ---'

    ! Start gmres iteration

      k=0
      do

        ! increment counter

        k=k+1

        ! start timer

        call cpu_time(ts)

        call gmres(num_var+1,rhs,solution,ipar,fpar,work)

        ! Go through output cases

          if(ipar(1).EQ.1) then           ! this means code wants Ax
            work(ipar(9):ipar(9)+num_var)=jacvec(work(ipar(8):ipar(8)+num_var),&
                                          ! overhead variables
                                          flux,keff,LL,U,Lf,Uf)
          elseif(ipar(1).EQ.2) then    ! this means code wants A^T x
            call stop_thor(27_li)
          elseif(ipar(1).EQ.3) then    ! this means code wants P(l)^{-1} z
            call stop_thor(27_li)
          elseif(ipar(1).EQ.3) then    ! this means code wants P(l)^{-T} z
            call stop_thor(27_li)
          elseif(ipar(1).EQ.5) then    ! this means code wants P(r)^{-1} z
            call stop_thor(27_li)
          elseif(ipar(1).EQ.6) then    ! this means code wants P(r)^{-T} z
            call stop_thor(27_li)
          elseif(ipar(1).EQ.10) then   ! call self-supplied stopping test
            call stop_thor(27_li)
          elseif(ipar(1).GT.0) then    ! shouldn't happen
            call stop_thor(27_li)
          elseif(ipar(1).EQ.0) then    ! successful solve
            error=fpar(6)/nres         ! pass back error
            kit=ipar(7)                ! pass back number of iters
            exit
          elseif(ipar(1).EQ.-1) then   ! maximum number of its reached
            error=fpar(6)/nres         ! pass back error
            kit=ipar(7)                ! pass back number of iters
            exit
          elseif(ipar(1).EQ.-2) then   ! insufficient workspace
            call stop_thor(27_li)
          elseif(ipar(1).LT.-2) then   ! some other error, investigate further
            write(6,*) 'SPARSKIT error is ', ipar(1)
          endif

        ! stop timer

        call cpu_time(te)

        ! write convergence monitor
        write(6,101) k,fpar(5),tol*nrhs,te-ts,' %%k'

      end do

    ! End gmres iteration

    ! Deallocate work

      deallocate(work)

101 FORMAT(1X,I6,3ES12.4,A)
102 FORMAT(1X,A)
  end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

  function jacvec(u_vec,&
       ! overhead variables
       flux,keff,LL,U,Lf,Uf)

    ! Pass input variables

      real(kind=d_t), dimension(num_var+1) :: u_vec

    ! Pass flux

      real(kind=d_t) :: flux(num_moments_v,namom,num_cells,egmax,niter)

    ! pass multiplication factor

      real(kind=d_t) :: keff

    ! Pre-computed matrices

      real(kind=d_t), dimension(num_moments_v,num_moments_v) :: LL, U
      real(kind=d_t), dimension(num_moments_f,num_moments_f) :: Lf, Uf

    ! Define return

      real(kind=d_t), dimension(num_var+1) :: jacvec

    ! define local

      real(kind=d_t)   :: dflux(num_moments_v,namom,num_cells,egmax,niter)
      real(kind=d_t)   :: dkeff
      integer(kind=li) :: eg,i,l,j,one,indx,ii,n,m

    ! allocate dflx and initialize

      one = 1_li

    ! compute dflx and dkeff from

      do eg=1,egmax
        do i=1,num_cells
          do n=1,namom
            do l=1,num_moments_v
              indx = (eg-1)*num_cells*namom*num_moments_v+&
                     (i-1)           *namom*num_moments_v+&
                     (n-1)                 *num_moments_v+&
                     l
              dflux(l,n,i,eg,1)=flux(l,n,i,eg,1)+jac_eps*u_vec(indx)
            end do
          end do
        end do
      end do
      dkeff = keff+jac_eps*u_vec(num_var+1)


    ! Compute jacvec

      jacvec = ( evaluate_residual(dflux,dkeff,LL,U,Lf,Uf) - residual ) &
       / jac_eps

  end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

  function converged(conv)
  !**********************************************************************
  !
  ! Function converged determines whether newton iteration is converged
  ! or not
  !
  !**********************************************************************

  ! pass input variables

    real(kind=d_t) :: conv
    logical :: converged

  ! convergence test

    converged = .false.
    if(maxval(abs(residual))<conv) converged = .true.

  end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

  subroutine newton_update(flux,keff)
  !**********************************************************************
  !
  ! Subroutine computes the new vector of unknowns using the du from last
  ! newton step
  !
  !**********************************************************************

    ! pass arguments

      real(kind=d_t) :: flux(num_moments_v,namom,num_cells,egmax,niter)
      real(kind=d_t) :: keff

    ! local variables

      integer(kind=li) :: eg,i,l,one,indx,n,m

    ! standard update

      one = 1_li

    ! update

      do eg = 1, egmax
        do i = 1, num_cells
          do n=1,namom
            do l = 1, num_moments_v
              indx = (eg-1)*num_cells*namom*num_moments_v+&
                     (i-1)           *namom*num_moments_v+&
                     (n-1)                 *num_moments_v+&
                     l
              flux(l,n,i,eg,one) = flux(l,n,i,eg,one) + du(indx)
            end do
          end do
        end do
      end do
      keff = keff + du(num_var+1)

  end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

  function forcing_factor(cnres,cnit)
  !**********************************************************************
  !
  ! Function forcing_factor computes the forcing factor
  !
  ! cnres - current norm of residual vector
  ! cnit  - current newton iteration
  !
  !**********************************************************************

    ! pass input arguments

      real(kind=d_t)   :: cnres
      integer(kind=li) ::cnit
      real(kind=d_t) :: forcing_factor
      real(kind=d_t) :: dummy

    ! easiest implementation

!      dummy = 1/(real(cnit,d_t)+2.0_d_t)
!      forcing_factor = min(dummy,cnres)
      forcing_factor = 0.01_d_t
  end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

  function pnorm(cresidual,pp)
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

      real(kind=d_t),dimension(:),intent(in)   :: cresidual
      integer(kind=li),intent(in) :: pp
      real(kind=d_t) :: pnorm

    ! length of vector

      integer(kind=li) :: length

    ! local variables

      integer :: i

    ! determine length of cresidual

      length = size(cresidual)

    ! compute pnorm

      pnorm = 0.0_d_t
      if (pp<10)  then
        do i = 1, length
          pnorm = pnorm + cresidual(i)**(real(pp,d_t))
        end do
        pnorm = pnorm**(1.0_d_t/real(pp,d_t))
      else
        pnorm = maxval(abs(cresidual))
      end if

  end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

  subroutine norm_eigmode(flux,tot_vol)
  !**********************************************************************
  !
  ! Subroutine norm_eigmode normalizes the eigenmode to 1
  !
  !**********************************************************************

    ! Pass flux

      real(kind=d_t) :: flux(num_moments_v,namom,num_cells,egmax,niter)
      real(kind=d_t),intent(in) :: tot_vol

    ! Local variable

      real(kind=d_t) :: s,value
      integer(kind=li) :: i,eg,l,n,m

      value = real(egmax,d_t) * tot_vol

    ! Add up average scalar fluxes

      s = 0.0_d_t
      l = 1_li
      do eg=1,egmax
        do i=1,num_cells
          s = s + cells(i)%volume*flux(l,1,i,eg,1)
        end do
      end do

    ! Normalize

      do eg=1,egmax
        do i=1,num_cells
          do n=1,namom
            do l = 1,num_moments_v
              flux(l,n,i,eg,1)=flux(l,n,i,eg,1)*value/s
            end do
          end do
        end do
      end do

  end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

  subroutine get_max_error
  !**********************************************************************
  !
  ! Subroutine get_max_error calculates the max residual
  ! per group, it only considers average flux moments
  !
  !**********************************************************************

    ! local arguments

      integer(kind=li) :: i,eg,indx
      real(kind=d_t),dimension(num_cells) :: temp

    ! loop over group

      do eg=1,egmax
        do i=1,num_cells
          indx=(eg-1)*num_cells*namom*num_moments_v+&
               (i-1)           *namom*num_moments_v+1
          temp(i)=abs(residual(indx))
        end do
        max_error(eg)=maxval(temp)
      end do

  end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------

  subroutine single_outer_iteration(flx,keff,reflected_flux,LL,U,Lf,Uf,prnt)
  !**********************************************************************
  !
  ! Subroutine single_outer_iteration performs a single outer
  ! iteration with lagged upscattering. This subroutine is used
  ! for JFNK method 1
  !
  !**********************************************************************

  ! Pass eigenvalue

    real(kind=d_t), intent(inout) :: keff

  ! Declare scalar flux

    real(kind=d_t) :: flx(num_moments_v,namom,num_cells,egmax,niter)

  !  Pre-computed matrices

    real(kind=d_t), dimension(num_moments_v,num_moments_v) :: LL, U
    real(kind=d_t), dimension(num_moments_f,num_moments_f) :: Lf, Uf

  ! Print flag for inner convergence monitor

    logical :: prnt

  ! Define fission source

    real(kind=d_t) :: fiss_src(num_moments_v,num_cells)

  ! Define source that is passed into inner iteration

    real(kind=d_t) :: src(num_moments_v,namom,num_cells,egmax)

  ! local variabels

    integer(kind=li) :: eg,eeg,i,j,l,m,n,alloc_stat,indx,k
    integer(kind=li) :: zero,one,two,three
    integer(kind=li) :: q,octant,egg

  ! dummy for max_errors

    real(kind=d_t),dimension(:),allocatable :: dmax_error

  ! Define reflected flux

    real(kind=d_t),dimension(num_moments_f,grs,8,nangle,egmax) :: reflected_flux

  ! set one,two, three

    one   = 1_li
    two   = 2_li
    three = 3_li
    zero  = 0_li

  ! Allocate dmax_error

    allocate(dmax_error(egmax))

  ! Initialize fiss_src and src

     src      = 0.0_d_t
     fiss_src = 0.0_d_t

  ! build fission source and initialize src with chi*fiss_src

    do eg = 1, egmax
      do i = 1, num_cells
        do l = 1, num_moments_v
          fiss_src(l,i) =fiss_src(l,i)                                          +&
              nu(reg2mat(cells(i)%reg),eg)%xs*fiss(reg2mat(cells(i)%reg),eg)%xs *&
              dens_fact(cells(i)%reg)*flx(l,1,i,eg,1)
        end do
      end do
    end do
    !
    do eg = 1, egmax
      do i = 1, num_cells
        do l = 1, num_moments_v
          src(l,1,i,eg)=chi(reg2mat(cells(i)%reg),eg)%xs*fiss_src(l,i)/keff
        end do
      end do
    end do

    !=========================================================================
    ! Compute upscattering
    !=========================================================================
    call compute_upscattering(flx, src)

  ! sweep through all groups

    do eg=1,egmax
      if  (page_refl.ne.0_li) then
        write(6, *) "Reflective Boundary values must currently be saved in JFNK."
      end if

      call inner_iteration(eg,flx(:,:,:,eg,2),src(:,:,:,eg),LL,U,Lf,Uf,grs,reflected_flux(:,:,:,:,eg),prnt)

      !=====================================================================
      ! Compute downscattering from eg to egg
      !=====================================================================
      call compute_downscattering(eg, flx, src)

      ! increase iteration count
      iit_count = iit_count + inner
    end do

  ! deallocate dmax_error

    deallocate(dmax_error)

  end subroutine

!-----------------------------------------------------------------------------------------
! THE END
!-----------------------------------------------------------------------------------------

  end module jfnk_module
