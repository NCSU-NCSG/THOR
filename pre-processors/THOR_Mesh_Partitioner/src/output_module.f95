MODULE output_module

  USE global_variables

  IMPLICIT NONE

CONTAINS

  SUBROUTINE write_files()
    !!Write output files
    DO i = 1, num_procs
      WRITE(out_mesh_file, '(A,I0)') "./part_meshes/"//TRIM(root_name)//".thrm.",i
      OPEN(UNIT = 24, FILE = out_mesh_file, ACTION = "WRITE", STATUS="REPLACE")
      WRITE(24,'(2(I0,4x),6x, A)') node_proc_indexes(i)-1, num_node, "!! Local Node Count, Global Node Count"
      WRITE(24,'(2(I0,4x),6x, A)') elem_proc_indexes(i)-1, num_elem, "!! Local Element Count, Global Node Count"
      WRITE(24,'(I0,10x, A)') 1, "!! # of Cell Blocks"
      WRITE(24,'(I0,10x, A)') num_side_sets, "!! # of Side Sets"
      WRITE(24,'(A)') "!! NODE LIST"
      DO j = 1, node_array_size
        IF (INT(node_proc_map(i,j,1)) .EQ. 0) EXIT
        WRITE(24, "(I8,X,3(D24.12))") INT(node_proc_map(i,j,1)), node_proc_map(i,j,2:)
      END DO
      WRITE(24,'(A)') "!! REGION & SRC LIST"
      DO j = 1, elem_array_size
        IF (INT(src_proc_map(i,j,1)) .EQ. 0) EXIT
        WRITE(24, "(5I8)") src_proc_map(i,j,:)
      END DO
      WRITE(24,'(A)') "!! ELEMENT LIST"
      DO j = 1, elem_array_size
        IF (INT(element_proc_map(i,j,1)) .EQ. 0) EXIT
        WRITE(24, "(5I8)") element_proc_map(i,j,:)
      END DO
      WRITE(24,'(A)') "!! BC FACE LIST"
      DO k = 1, num_side_sets
        WRITE(24,'(I0)') BC_proc_indexes(i,k) - 1
        DO j = 1, BC_proc_indexes(i,k) - 1
          WRITE(24,'(3I8)') BC_proc_map(i,k,j,:)
        END DO
      END DO
      WRITE(24,'(A)') "!! ADJACENCY LIST"
      WRITE(24,'(I0)') adj_proc_indexes(i) - 1
      DO j = 1, adj_proc_indexes(i)-1
        WRITE(24, "(5I8)") adj_proc_map(i,j,:)
      END DO
      CLOSE(24)
    END DO
    WRITE(*,*)
    WRITE(*,*) TRIM(out_mesh_file(1:INDEX(out_mesh_file,".",BACK=.TRUE.)))//"*"
  END SUBROUTINE

END MODULE
