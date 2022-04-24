MODULE global_variables
  !***********************************************************************
  ! This module stores global variables that are used throughout the
  ! THOR code. Variables are only defined as global variables if they
  ! do not change throughout the iteration procedure, e.g. fluxes are
  ! not defined to be global.
  !***********************************************************************
  USE types
  USE parameter_types
  USE filename_types
  USE vector_types
  USE cross_section_types
  USE geometry_types
  USE angle_types
  USE multindex_types
  USE integer_array_tools
  USE source_types
  IMPLICIT NONE

  ! from parameters

  INTEGER(kind=li) :: conv_flag
  INTEGER(kind=li) :: tot_nInners
  REAL(kind=d_t)   :: k_error
  REAL(kind=d_t)   :: f_error
  INTEGER(kind=li) :: namom
  INTEGER(kind=li) :: niter
  INTEGER(kind=li) :: most_thermal
  REAL(kind=d_t)   :: max_outer_error

  ! General input variables

  INTEGER(kind=li)  :: problem, space_ord, scatt_ord,quad_ord, xs_ord,&
        egmax, num_mat, max_outer, max_inner, num_src_mat,&
        vtk_flux_output, finflow ,num_moments_v,num_moments_f,&
        vtk_mat_output,vtk_reg_output,execution,print_xs_flag,&
        upscattering,vtk_src_output,multiplying,print_conv

  REAL(kind=d_t) :: inner_conv, outer_conv

  ! Postprocessing variables

  LOGICAL :: glob_do_cartesian_mesh
  INTEGER(kind=li) :: glob_cmap_nx, glob_cmap_ny, glob_cmap_nz
  REAL(kind=d_t) :: glob_cmap_min_x, glob_cmap_min_y, glob_cmap_min_z
  REAL(kind=d_t) :: glob_cmap_max_x, glob_cmap_max_y, glob_cmap_max_z

  INTEGER(kind=li) :: number_point_flux_locations
  REAL(kind=d_t), ALLOCATABLE :: point_flux_locations(:,:)

  !External Sources
  TYPE(source_type), ALLOCATABLE :: ext_src(:)
  INTEGER(kind=li)      , DIMENSION(:)      , ALLOCATABLE :: source_ids
  INTEGER(kind=li) :: src_id_min,src_id_max

  ! Inflow flux derived type

  REAL(kind=d_t), DIMENSION(:,:,:,:,:), ALLOCATABLE :: binflx

  !    type(bsource_mat), dimension(:), allocatable :: bsrc_mat
  !    type(bsource_str), dimension(:,:), allocatable :: bsrc_str
  !    type(bsource_moments), dimension(:), allocatable :: bsrc_m


  ! Cross Sections
  TYPE(xs_material_type) , DIMENSION(:)      , ALLOCATABLE :: xs_mat
  REAL(kind=d_t)         , DIMENSION(:)      , ALLOCATABLE :: eg_bounds
  INTEGER(kind=li)      , DIMENSION(:)      , ALLOCATABLE :: material_ids
  INTEGER(kind=li) :: mat_id_min,mat_id_max

  ! Geometry types

  INTEGER(kind=li) :: num_vert, num_cells, num_cell_blk,num_side_sets, &
        side_cells, adjacent_cells
  INTEGER(kind=li) :: vside_cells    ! # vacuum       boundary faces
  INTEGER(kind=li) :: rside_cells    ! # reflective   boundary faces
  INTEGER(kind=li) :: fside_cells    ! # fixed inflow boundary faces
  TYPE(vertex)       , DIMENSION(:)  , ALLOCATABLE :: vertices
  TYPE(cell)         , DIMENSION(:)  , ALLOCATABLE :: cells
  TYPE(boundary_cell), DIMENSION(:)  , ALLOCATABLE :: b_cells       !  all boundary faces
  TYPE(boundary_cell), DIMENSION(:)  , ALLOCATABLE :: rb_cells      !  reflective boundary faces
  TYPE(boundary_cell), DIMENSION(:)  , ALLOCATABLE :: vb_cells      !  vacuum boundary faces
  TYPE(boundary_cell), DIMENSION(:)  , ALLOCATABLE :: fb_cells      !  fixed inflow boundary faces
  INTEGER(kind=li)   , DIMENSION(:)  , ALLOCATABLE :: refl_face_tpe !  supplementary array indicating  type of
  !  reflective bc
  TYPE(list)         , DIMENSION(:,:), ALLOCATABLE :: adjacency_list

  ! Quadrature

  TYPE(ordinate), DIMENSION(:), ALLOCATABLE :: quadrature

  ! Spatial moment indices

  TYPE(indices_v), DIMENSION(:), ALLOCATABLE :: index_v
  TYPE(indices_f), DIMENSION(:), ALLOCATABLE :: index_f

  ! New input variables

  INTEGER(kind=li) :: dump_flag,quad_tpe,nangle,inguess_flag
  REAL(kind=d_t )  :: k_conv
  INTEGER(kind=li) :: eig_switch,rd_restart,rd_max_kit,rd_method,ipow
  INTEGER(kind=li) :: tot_kit, inner, outer, nit

  ! groupwise maximum flux error

  REAL(kind=d_t), DIMENSION(:), ALLOCATABLE :: max_error

  ! start and finish times

  REAL(kind=d_t)   :: start, finish

  ! Spherical harmonics variables

  REAL(kind=d_t), DIMENSION(:,:,:), ALLOCATABLE :: Ysh

  ! Outward normal vectors

  TYPE(vector), DIMENSION(:,:), ALLOCATABLE :: outward_normal

  ! mu, eta and xi-mates

  INTEGER(kind=li),DIMENSION(8) :: mu_mate  = (/2,1,4,3,6,5,8,7/), &
        eta_mate = (/4,3,2,1,8,7,6,5/), &
        xi_mate  = (/5,6,7,8,1,2,3,4/)

  ! ordering of the octants for sweep

  INTEGER(kind=li), DIMENSION(8) :: ordering=(/7,6,8,5,3,2,4,1/)
  INTEGER(kind=li), DIMENSION(8) :: octants_to_sweep, ordered_octants_to_sweep

  ! Additional input switches

  INTEGER(kind=li) :: sweep_tpe         ! 1. Pre-computed
  INTEGER(kind=li) :: outer_acc         ! 1. no acceleration, 2. Error mode extrapolation

  ! Variables for error mode extrapolation

  REAL(kind=d_t) :: extol,exmax

  ! Sweep path stores the order of the sweep for all directions

  INTEGER(kind=li),ALLOCATABLE :: sweep_path(:,:,:)
  INTEGER(kind=li) :: page_sweep

  ! Region to material ID map

  INTEGER(kind=li),ALLOCATABLE :: reg2mat(:)
  INTEGER(kind=li) :: minreg,maxreg

  ! Scattering XS multiplier

  REAL(kind=d_t),ALLOCATABLE :: scat_mult(:,:)
  INTEGER(kind=li) :: scat_mult_flag,neven

  ! Memorize the connections severed for elimination of cycles

  TYPE(cycle_breaker),ALLOCATABLE :: eldep(:,:,:)
  INTEGER(kind=li),ALLOCATABLE    :: neldep(:,:)

  ! Memory consumption estimate

  REAL(kind=d_t) :: mem_req

  ! Flags to page out specific angular information

  INTEGER(kind=li) :: page_refl
  INTEGER(kind=li) :: page_iflw
  INTEGER(kind=li) :: eg_iflw

  ! Flag for determining if face is in cycle
  INTEGER(kind=1),ALLOCATABLE :: is_cycle(:,:)

  ! density factors

  REAL(kind=d_t),ALLOCATABLE  :: dens_fact(:)
  REAL(kind=d_t),ALLOCATABLE  :: reg_vol(:)
  INTEGER(kind=li)            :: dfact_opt

  ! Parallel Data

  INTEGER, ALLOCATABLE:: parallel_map_g2l(:,:)
  INTEGER, ALLOCATABLE:: parallel_map_l2g(:,:)
  INTEGER:: rank
  REAL*8:: parallel_timing(4,2)

  !Input/output units
  INTEGER :: stdout_unit

  !log variables
  INTEGER, PARAMETER :: log_unit=9999
END MODULE global_variables
