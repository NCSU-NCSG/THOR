MODULE SDD_global_variables
USE geometry_types
IMPLICIT NONE


INTEGER,ALLOCATABLE :: SDD_cells_g2l_G(:)
INTEGER,ALLOCATABLE :: SDD_cells_l2g_G(:)
INTEGER,ALLOCATABLE :: SDD_vert_g2l_G(:)
INTEGER,ALLOCATABLE :: SDD_vert_l2g_G(:)
INTEGER,ALLOCATABLE :: unique_neighs(:)
INTEGER,ALLOCATABLE :: num_ang_tosend(:)
INTEGER,ALLOCATABLE :: comm_instructions(:,:,:)
INTEGER,ALLOCATABLE :: pack_instructions(:,:)
INTEGER,ALLOCATABLE :: unpack_instructions(:,:)
INTEGER,ALLOCATABLE :: SDD_refl_BC_instructions(:,:)
INTEGER,ALLOCATABLE :: SDD_send_handles_G(:,:), SDD_recv_handles_G(:,:)
REAL(8),ALLOCATABLE :: psi_comm(:,:),psi_comm_send(:,:)

INTEGER :: tot_num_vert_G
INTEGER :: tot_num_cells_G
INTEGER :: SDD_side_cells
INTEGER :: max_ang_comm
INTEGER :: num_proc
INTEGER :: switch=1
INTEGER :: ind_k

REAL(8) :: reg_dbg_info(3)=0.0d0

INTEGER, PARAMETER :: num_dbg_records=25000
REAL(8) :: per_iter_time(num_dbg_records)=0.0d0, per_comm_time(num_dbg_records)=0.0d0
REAL(8) :: IOtime=0.0d0, comm_time=0.0d0,comm_start_time

!Multi processor adjacency list mapping. Maps local cell to global neighbor
!Includes boundary and off processor neighbor data
TYPE(SDD_list_entry), DIMENSION(:,:), ALLOCATABLE :: SDD_adjacency_list

INTEGER(kind=li),ALLOCATABLE :: IPVTtemp(:)
REAL(8),ALLOCATABLE :: Jphitemp(:,:),Kphitemp(:,:),Jpsitemp(:,:),IJtemp(:,:)

END MODULE SDD_global_variables
