module solver_module
!***********************************************************************
!
! Solver module initializes fluxes and sources, then determines type of
! problem and calls outer iteration
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

! Use modules that pertain setting up problem

  use ahotc_matrix_module
  use sph_harmonics_module
  use outer_iteration_module
  use jfnk_module
  use termination_module
  
  implicit none

contains

  !> Prepares some data and allocates memory, then hands the solve over to
  !> the outer_iteration_ext subroutine
  subroutine solver_ext(flux)
  !*********************************************************************
  !
  ! Subroutine solver calls outer iteration subroutine
  !
  !*********************************************************************

  ! Pass scalar flux 

    real(kind=d_t) :: flux(num_moments_v,namom,num_cells,egmax,niter) 

  ! Define temporary variables

    integer(kind=li) :: alloc_stat, k
    real(kind=d_t), dimension(:,:), allocatable :: M, LL, U, &
         Mf, Lf, Uf

  ! Allocate and initialize all necessary derived types

    allocate(M(num_moments_v,num_moments_v),&
         LL(num_moments_v,num_moments_v),&
         U(num_moments_v,num_moments_v),&
         Mf(num_moments_f,num_moments_f),&
         Lf(num_moments_f,num_moments_f),&
         Uf(num_moments_f,num_moments_f),&
         stat=alloc_stat);M=zero;&
         LL=zero;U=zero;Mf=zero;Lf=zero;Uf=zero;
    if(alloc_stat /= 0) call stop_thor(2_li)

    
    allocate(Ysh(nangle,8,namom),stat=alloc_stat)
    if(alloc_stat /=0) call stop_thor(2_li) 
    
  ! Pre-compute and apply LU decomposition to 'mass matrices'

       call ahotc_matrix(num_moments_v,num_moments_f,index_v,&
            index_f,M,LL,U,Mf,Lf,Uf)

  ! Pre-compute the spherical harmonics basis functions

       call spherical_harmonics(scatt_ord,nangle,quadrature,namom,Ysh)

  ! Initate outer iteration
      if (rank .eq. 0) then
        write(6,*) '-- Commencing fixed source computation.'
      end if
      call outer_iteration_ext(flux,LL,U,Lf,Uf)

  ! Deallocate temporary arrays

      deallocate(M,Mf,LL,U,Lf,Uf)

  end subroutine solver_ext

!-----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------- 

  !> Prepares some data and allocates memory, then hands the solve over to
  !> the outer_iteration_eig subroutine or the JFNK solver. JFNK is a different
  !> solution method developed by Dan Gill. 
  subroutine solver_eig(flux,keff)
  !*********************************************************************
  !
  ! Subroutine solver calls outer iteration subroutine
  !
  !*********************************************************************



  ! Pass eigenvalue and scalar flux

    real(kind=d_t), intent(inout) :: keff
    real(kind=d_t) :: flux(num_moments_v,namom,num_cells,egmax,niter) 

  ! Define temporary variables

    integer(kind=li) :: alloc_stat, k
    real(kind=d_t), dimension(:,:), allocatable :: M, LL, U, &
         Mf, Lf, Uf

  ! Allocate and initialize all necessary derived types

    allocate(M(num_moments_v,num_moments_v),&
         LL(num_moments_v,num_moments_v),&
         U(num_moments_v,num_moments_v),&
         Mf(num_moments_f,num_moments_f),&
         Lf(num_moments_f,num_moments_f),&
         Uf(num_moments_f,num_moments_f),&
         stat=alloc_stat);M=zero;&
         LL=zero;U=zero;Mf=zero;Lf=zero;Uf=zero;
    if(alloc_stat /= 0) call stop_thor(2_li)

    allocate(Ysh(nangle,8,namom),stat=alloc_stat)
    if(alloc_stat /=0) call stop_thor(2_li) 
    
  ! Pre-compute and apply LU decomposition to 'mass matrices'

       call ahotc_matrix(num_moments_v,num_moments_f,index_v,&
            index_f,M,LL,U,Mf,Lf,Uf)

  ! Pre-compute the spherical harmonics basis functions

       call spherical_harmonics(scatt_ord,nangle,quadrature,namom,Ysh)

  ! Initate outer iteration
       if (eig_switch .eq. 0) then
        if (rank .eq. 0) then
          write(6,*) '-- Commencing eigenvalue computation.'
          write(6,*) '-- A power iteration computation is executed.'
        end if
         call outer_iteration_eig(flux,keff,LL,U,Lf,Uf) 
       else
        if (rank .eq. 0) then
          write(6,*) '-- Commencing eigenvalue computation.'
          write(6,*) '-- A JFNK computation is executed.'
        end if
         call set_jfnk
    
         call do_jfnk(flux,keff,LL,U,Lf,Uf)    
    
         call clean_jfnk         
             
       end if

  ! Deallocate temporary arrays

       deallocate(M,Mf,LL,U,Lf,Uf)

  end subroutine solver_eig

!-----------------------------------------------------------------------------------------
! End
!-----------------------------------------------------------------------------------------  
  
end module solver_module
