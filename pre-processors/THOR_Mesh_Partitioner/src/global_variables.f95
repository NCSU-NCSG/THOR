MODULE global_variables
  IMPLICIT NONE

  CHARACTER(60):: root_name, thrm_file, metis_in_file, out_mesh_file
  CHARACTER(60):: metis_elem_out_file, metis_node_out_file, vtk_file
  CHARACTER(60):: num_procs_str
  CHARACTER(60):: metis_in_path = "./metis_in/"
  CHARACTER(60):: thrm_path = "./meshes/"
  CHARACTER(60):: stats=''
  CHARACTER:: com_check
  INTEGER::num_procs, num_elem, num_node, num_BC_faces, ierr
  INTEGER:: i,j,k,q, temp_int, num_adj, temp_proc_self, temp_proc_neighbor
  INTEGER, PARAMETER :: k11 = selected_int_kind(12)
  INTEGER(kind = k11)::mem_usage = 0
  REAL:: t_start, t_end1, t_end2
  LOGICAL:: temp_log
  INTEGER:: adj_array_size, elem_array_size, node_array_size, BC_array_size
  INTEGER:: num_side_sets = 0, num_cell_blocks = 0

  INTEGER:: temp_adj_line(4), temp_BC_line(3)
  DOUBLE PRECISION,ALLOCATABLE:: temp_node_line(:,:)

  INTEGER:: load_MIN, load_MAX
  DOUBLE PRECISION:: ext_int_ratio_MIN, ext_int_ratio_MAX
  INTEGER:: num_spikes
  INTEGER:: neigh_MIN, neigh_MAX
  INTEGER, ALLOCATABLE:: neigh_count(:)
  DOUBLE PRECISION:: load_AVG, ext_int_ratio_AVG, neigh_AVG
  LOGICAL:: contig, refl_faces_present

  DOUBLE PRECISION, ALLOCATABLE:: node_proc_map(:,:,:)
  INTEGER, ALLOCATABLE:: elem_list(:,:)
  INTEGER, ALLOCATABLE:: elem_proc_indexes(:)
  INTEGER, ALLOCATABLE:: node_proc_indexes(:)
  INTEGER, ALLOCATABLE:: elem_node_proc_map(:)
  INTEGER, ALLOCATABLE:: element_proc_map(:,:,:)
  INTEGER, ALLOCATABLE:: src_proc_map(:,:,:)
  INTEGER, ALLOCATABLE:: src_proc_indexes(:)
  INTEGER, ALLOCATABLE:: adj_proc_map(:,:,:)
  INTEGER, ALLOCATABLE:: adj_proc_indexes(:)
  INTEGER, ALLOCATABLE:: BC_proc_map(:,:,:,:)
  INTEGER, ALLOCATABLE:: BC_proc_indexes(:,:)
  INTEGER, ALLOCATABLE:: proc_map(:)
  
  CHARACTER(60):: seed

END MODULE
