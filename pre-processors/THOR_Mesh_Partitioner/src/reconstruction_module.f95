MODULE reconstruction_module

  USE global_variables

  IMPLICIT NONE

CONTAINS

  SUBROUTINE reconstruct_meshes()
    IMPLICIT NONE

    !!Generate Temp Build Files
    ! Open a node, element, BC, Mat, Adj File for Each Proc
    ! Stripe out relevant lines
    elem_array_size = CEILING(num_elem/(num_procs*1.0) * 2.0)
    node_array_size = num_elem/num_procs * 4 + 10


    ALLOCATE(element_proc_map(num_procs, elem_array_size, 5))
    ALLOCATE(node_proc_map(num_procs, node_array_size, 4))
    ALLOCATE(elem_node_proc_map(num_node))
    ALLOCATE(elem_proc_indexes(num_procs), node_proc_indexes(num_procs))
    ALLOCATE(src_proc_map(num_procs,elem_array_size, 3))
    ALLOCATE(src_proc_indexes(num_procs))
    ALLOCATE(proc_map(num_elem))
    ALLOCATE(temp_node_line(num_node,4))

    element_proc_map = 0
    elem_node_proc_map = 0
    src_proc_map = 0
    node_proc_map = 0
    elem_proc_indexes = 1
    node_proc_indexes = 1
    src_proc_indexes = 1

    OPEN(UNIT = 23, FILE = metis_elem_out_file, ACTION = "READ", STATUS= "OLD")
    OPEN(UNIT = 21, FILE = thrm_file, ACTION = "READ", STATUS= "OLD")
    DO i = 1, num_elem
      READ(23,*) proc_map(i)
    END DO
    REWIND(23)
    REWIND(21)
    DO i = 1,4+num_elem+num_node
      READ(21,*)com_check
      IF (com_check.NE.'!') BACKSPACE(21)
      READ(21,*)
    END DO

    !!Map original elements to new proc locations
    READ(21,*)com_check
    IF (com_check.NE.'!') BACKSPACE(21)
    DO i = 1, num_elem
      READ(23,*) temp_int
      READ(21,*) element_proc_map(temp_int+1, elem_proc_indexes(temp_int+1), :)
      elem_proc_indexes(temp_int+1) = elem_proc_indexes(temp_int+1) + 1
    END DO

    !! Read in Node Data
    REWIND(21)
    DO i = 1, 4
      READ(21,*)
    END DO
    READ(21,*)com_check
    IF (com_check.NE.'!') BACKSPACE(21)
    DO i =1, num_node
      READ(21,*) temp_node_line(i,:)
    END DO

    !!Get Node Lists Based on Elements to remove duplication
    temp_int = SIZE(element_proc_map(1,:,1))
    DO i = 1, num_procs
      elem_node_proc_map=0
      DO j = 1, temp_int
        DO k = 2, 5
          IF (element_proc_map(i,j,k) .NE. 0) elem_node_proc_map(element_proc_map(i,j,k)) = 1
        END DO
      END DO
      !Before progressing to the next processor, we get the node list for this proc,
      !so that we can overwrite element_node_proc_map
      DO q =1, num_node
        IF (elem_node_proc_map(INT(temp_node_line(q,1))) .EQ. 1) THEN
          node_proc_map(i, node_proc_indexes(i),:) = temp_node_line(q,:)
          node_proc_indexes(i) = node_proc_indexes(i) + 1
        END IF
      END DO
    END DO

    !!Read in Region / SRC
    REWIND(23)
    READ(21,*)com_check
    IF (com_check.NE.'!') BACKSPACE(21)
    DO i = 1, num_elem
      READ(23,*) temp_int
      READ(21,*) src_proc_map(temp_int+1, src_proc_indexes(temp_int+1), :)
      src_proc_indexes(temp_int+1) = src_proc_indexes(temp_int+1) + 1
    END DO

    !! Read in Adjacency and Compute Proc Adjacency Lookup
    adj_array_size = CEILING(8 * (num_elem / num_procs) * 1.25)
    ALLOCATE(adj_proc_map(num_procs,adj_array_size, 5))
    ALLOCATE(adj_proc_indexes(num_procs))
    adj_proc_map = -1
    adj_proc_indexes = 1
    READ(21,*)com_check
    IF (com_check.NE.'!') BACKSPACE(21)
    DO i = 1, num_elem
      READ(21,*)
    END DO
    READ(21,*)com_check
    IF (com_check.NE.'!') BACKSPACE(21)
    DO j = 1, num_side_sets
      READ(21,*) num_BC_faces
      DO i = 1, num_BC_faces
        READ(21,*)
      END DO
    END DO
    READ(21,*)com_check
    IF (com_check.NE.'!') BACKSPACE(21)
    READ(21,*) num_adj
    DO i=1, num_adj
      READ(21,*) temp_adj_line
      temp_proc_self = proc_map(temp_adj_line(1))+1
      IF (temp_adj_line(3) .NE. 0) THEN
        temp_proc_neighbor = proc_map(temp_adj_line(3))+1
      ELSE
        temp_proc_neighbor = 0
      END IF
      adj_proc_map(temp_proc_self, adj_proc_indexes(temp_proc_self),1:4) = temp_adj_line
      adj_proc_map(temp_proc_self, adj_proc_indexes(temp_proc_self),5) = temp_proc_neighbor
      adj_proc_indexes(temp_proc_self) = adj_proc_indexes(temp_proc_self) + 1
    END DO

    !! Read in BC Faces & Compute New Ones
    BC_array_size = CEILING(4 * (num_elem / num_procs) * 1.5)
    ALLOCATE(BC_proc_map(num_procs,num_side_sets, BC_array_size, 3 ))
    ALLOCATE(BC_proc_indexes(num_procs, num_side_sets))
    BC_proc_map = -1
    BC_proc_indexes = 1
    REWIND(21)
    temp_int = num_node+ 2*num_elem + 4
    DO i =1, temp_int
      READ(21,*)com_check
      IF (com_check.NE.'!') BACKSPACE(21)
      READ(21,*)
    END DO
    !!Get BCs from orginal file
    READ(21,*)com_check
    IF (com_check.NE.'!') BACKSPACE(21)
    DO j = 1, num_side_sets
      READ(21,*) num_BC_faces
      refl_faces_present = .FALSE.
      DO i = 1, num_BC_faces
        READ(21,*) temp_BC_line
        IF(temp_BC_line(3) .EQ. 1 .AND. .NOT. refl_faces_present) THEN
          refl_faces_present = .TRUE.
        END IF
        temp_proc_self = proc_map(temp_BC_line(1)) + 1
        BC_proc_map(temp_proc_self,j, BC_proc_indexes(temp_proc_self,j),:) = temp_BC_line
        !BC_proc_map(temp_proc_self, BC_proc_indexes(temp_proc_self),3) = 0
        BC_proc_indexes(temp_proc_self,j) = BC_proc_indexes(temp_proc_self,j) + 1
      END DO
    END DO
    !!Build new BCs from adjacency list
    DO i = 1, num_procs
      DO j = 1, adj_proc_indexes(i)-1
        IF (adj_proc_map(i,j,5) .EQ. 0 .OR. adj_proc_map(i,j,5) .EQ. i ) THEN
        ELSE
          BC_proc_map(i, 1, BC_proc_indexes(i,1),1) = adj_proc_map(i,j,1)
          BC_proc_map(i, 1, BC_proc_indexes(i,1),2) = adj_proc_map(i,j,2)
          BC_proc_map(i, 1, BC_proc_indexes(i,1),3) = -1
          BC_proc_indexes(i,1) = BC_proc_indexes(i,1) + 1
        END IF
      END DO
    END DO

  END SUBROUTINE

END MODULE
