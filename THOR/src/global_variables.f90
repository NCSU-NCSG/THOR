module global_variables
!***********************************************************************
! This module stores global variables that are used throughout the
! THOR code. Variables are only defined as global variables if they
! do not change throughout the iteration procedure, e.g. fluxes are
! not defined to be global.
!***********************************************************************
use types
use parameter_types
use filename_types
use vector_types
use cross_section_types
use geometry_types
use angle_types
use multindex_types
implicit none

  ! from parameters

    integer(kind=li) :: conv_flag
    integer(kind=li) :: tot_nInners
    real(kind=d_t)   :: k_error
    real(kind=d_t)   :: f_error
    integer(kind=li) :: namom
    integer(kind=li) :: niter
    integer(kind=li) :: most_thermal
    real(kind=d_t)   :: max_outer_error

  ! General input variables

    integer(kind=li)  :: problem, space_ord, scatt_ord,quad_ord, xs_ord,       &
                         egmax, num_mat, max_outer, max_inner, num_src_mat,    &
                         vtk_flux_output, finflow ,num_moments_v,num_moments_f,&
                         vtk_mat_output,vtk_reg_output,execution,print_xs_flag,&
                         upscattering,vtk_src_output,multiplying,print_conv

    real(kind=d_t)    :: inner_conv, outer_conv

  ! Postprocessing variables

    logical :: glob_do_cartesian_mesh
    integer(kind=li) :: glob_cmap_nx, glob_cmap_ny, glob_cmap_nz
    real(kind=d_t) :: glob_cmap_min_x, glob_cmap_min_y, glob_cmap_min_z
    real(kind=d_t) :: glob_cmap_max_x, glob_cmap_max_y, glob_cmap_max_z

  ! Source types

    integer, dimension(:), allocatable :: src_mat
    real(kind=d_t),allocatable         :: src_m  (:,:,:)
    real(kind=d_t),allocatable         :: src_str(:,:)

  ! Inflow flux derived type

    real(kind=d_t),dimension(:,:,:,:,:),allocatable :: binflx

!    type(bsource_mat), dimension(:), allocatable :: bsrc_mat
!    type(bsource_str), dimension(:,:), allocatable :: bsrc_str
!    type(bsource_moments), dimension(:), allocatable :: bsrc_m


  ! Cross Sections

    type(cross_section_mat), dimension(:)      , allocatable :: xs_mat
    type(cross_section)    , dimension(:,:)    , allocatable :: chi,eg_bounds, fiss, nu, sigma_t
    type(cross_section)    , dimension(:,:)    , allocatable :: tsigs
    type(cross_section)    , dimension(:,:,:,:), allocatable :: sigma_scat

  ! Geometry types

    integer(kind=li) :: num_vert, num_cells, num_cell_blk,num_side_sets, &
                        side_cells, adjacent_cells
    integer(kind=li) :: vside_cells    ! # vacuum       boundary faces
    integer(kind=li) :: rside_cells    ! # reflective   boundary faces
    integer(kind=li) :: fside_cells    ! # fixed inflow boundary faces
    type(vertex)       , dimension(:)  , allocatable :: vertices
    type(cell)         , dimension(:)  , allocatable :: cells
    type(boundary_cell), dimension(:)  , allocatable :: b_cells       !  all boundary faces
    type(boundary_cell), dimension(:)  , allocatable :: rb_cells      !  reflective boundary faces
    type(boundary_cell), dimension(:)  , allocatable :: vb_cells      !  vacuum boundary faces
    type(boundary_cell), dimension(:)  , allocatable :: fb_cells      !  fixed inflow boundary faces
    integer(kind=li)   , dimension(:)  , allocatable :: refl_face_tpe !  supplementary array indicating  type of
                                                                      !  reflective bc
    type(list)         , dimension(:,:), allocatable :: adjacency_list

  ! Quadrature

    type(ordinate), dimension(:), allocatable :: quadrature

  ! Spatial moment indices

    type(indices_v), dimension(:), allocatable :: index_v
    type(indices_f), dimension(:), allocatable :: index_f

  ! New input variables

    integer(kind=li) :: dump_flag,quad_tpe,nangle,inguess_flag
    real(kind=d_t )  :: k_conv
    integer(kind=li) :: eig_switch,rd_restart,rd_max_kit,rd_method,ipow
    integer(kind=li) :: tot_kit, inner, outer, nit

  ! groupwise maximum flux error

    real(kind=d_t), dimension(:), allocatable :: max_error

  ! start and finish times

    real(kind=d_t)   :: start, finish

  ! Spherical harmonics variables

    real(kind=d_t), dimension(:,:,:), allocatable :: Ysh

  ! Outward normal vectors

    type(vector), dimension(:,:), allocatable :: outward_normal

  ! mu, eta and xi-mates

    integer(kind=li),dimension(8) :: mu_mate  = (/2,1,4,3,6,5,8,7/), &
                                     eta_mate = (/4,3,2,1,8,7,6,5/), &
                                     xi_mate  = (/5,6,7,8,1,2,3,4/)

  ! ordering of the octants for sweep

    integer(kind=li), dimension(8) :: ordering=(/7,6,8,5,3,2,4,1/)


  ! Additional input switches

    integer(kind=li) :: sweep_tpe         ! 1. Pre-computed
    integer(kind=li) :: outer_acc         ! 1. no acceleration, 2. Error mode extrapolation

  ! Variables for error mode extrapolation

    real(kind=d_t) :: extol,exmax

  ! Sweep path stores the order of the sweep for all directions

    integer(kind=li),allocatable :: sweep_path(:,:,:)
    integer(kind=li) :: page_sweep

  ! Region to material ID map

    integer(kind=li),allocatable :: reg2mat(:)
    integer(kind=li) :: minreg,maxreg

  ! Scattering XS multiplier

    real(kind=d_t),allocatable :: scat_mult(:,:)
    integer(kind=li) :: scat_mult_flag,neven

  ! Memorize the connections severed for elimination of cycles

    type(cycle_breaker),allocatable :: eldep(:,:,:)
    integer(kind=li),allocatable    :: neldep(:,:)

  ! Memory consumption estimate

    real(kind=d_t) :: mem_req

  ! Flags to page out specific angular information

    integer(kind=li) :: page_refl
    integer(kind=li) :: page_iflw
    integer(kind=li) :: eg_iflw

  ! Flag for determining if face is in cycle
    integer(kind=1),allocatable :: is_cycle(:,:)

  ! density factors

    real(kind=d_t),allocatable  :: dens_fact(:)
    real(kind=d_t),allocatable  :: reg_vol(:)
    integer(kind=li)            :: dfact_opt

  ! Parallel Data

    integer, allocatable:: parallel_map_g2l(:,:)
    integer, allocatable:: parallel_map_l2g(:,:)
    integer:: rank
    real*8:: parallel_timing(4,2)
end module
