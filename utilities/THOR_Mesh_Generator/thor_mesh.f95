!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! THOR MESH GENERATION UTILITY
!   THOR_mesh module:
!     Handles output of *.thrm (thor_mesh) files
!     This is intended to enable an *.e (exodus II) to *.thm toolchain for
!     probelm development purposes
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
MODULE thor_mesh

USE globals
IMPLICIT NONE
CONTAINS

  SUBROUTINE outputThorMesh()

    INTEGER :: i
    INTEGER :: j
    INTEGER :: position
    INTEGER :: element
    INTEGER :: neighbor
    INTEGER :: current_face
    INTEGER :: current_neighbor_face
    INTEGER :: local_node_list_gmesh(3), order_gmesh(3)
    INTEGER :: local_node_list_adjacency(4), order_adjacency(4)

    !Open File
    OPEN(UNIT = out_unit, file = out_file, ACTION = "write", STATUS = "replace")

    !Header
    WRITE(out_unit, '(I0)') node_count
    WRITE(out_unit, '(I0)') element_count
    !TODO: Find out what these 1 entries are for:
    WRITE(out_unit, '(I0)') 1
    WRITE(out_unit, '(I0)') 1

    !Block I: Node Information
    DO i = 1, node_count
      WRITE(out_unit, '(I11, 3ES26.16)') i, node_list(i,:)
    END DO

    !Block II: Element int int (Source Def?)
    DO i = 1, element_count
      position = mapIndexOf(block_id(i), block_id_map(:, 1))
      WRITE(out_unit, '(I0,2(X,I0))') i, block_id_map(position, 2),&
                                      source_id_map(position, 2)
    END DO

    !Block III: Element Composition
    DO i = 1, element_count
      WRITE(out_unit, '(I0, 4(X,I0))') i, element_list(i,:)
    END DO

    !# Of BC Definitions
    WRITE(out_unit, '(I0)') bc_count

    ! TODO: this search can be costly for many bcs; in this case we need to
    ! come up with a better lookup [consider string-hashes]
    IF (SIZE(boundary_element_list) .NE. bc_count) &
      CALL generateErrorMessage(err_code_default, err_fatal, &
        'bc count from gmesh file and adjacency list are inconsistent')

    ! TODO: side set reassignment
    ! TODO: better algorithm
    ! TODO: if no better algorithm at least make it a pretty function
    DO i = 1, bc_count
      local_node_list_gmesh = bc_list(i, :)
      CALL quickSortInteger(local_node_list_gmesh, order_gmesh)
      DO j = 1, bc_count
        element = boundary_element_list(j)
        current_face = boundary_face_list(j)
        local_node_list_adjacency = element_list(element, :)
        local_node_list_adjacency(current_face + 1) = -1
        CALL quickSortInteger(local_node_list_adjacency, order_adjacency)
        IF (local_node_list_gmesh(1) .EQ. local_node_list_adjacency(2) .AND. &
            local_node_list_gmesh(2) .EQ. local_node_list_adjacency(3) .AND. &
            local_node_list_gmesh(3) .EQ. local_node_list_adjacency(4)) THEN
          WRITE(out_unit, '(I0, 2(X,I0))') element, current_face, side_set(i)
        END IF
      END DO
    END DO

    ! write adjacency list
    WRITE(out_unit, '(I0)') 4 * element_count
    DO i = 1, element_count
      DO j = 1, 4
        CALL adjacencyListEntry(i, j, neighbor, current_face, current_neighbor_face)
        WRITE(out_unit, '(I0, 3(X, I0))') i, current_face, neighbor, current_neighbor_face
      END DO
    END DO

    CLOSE(out_unit)

  END SUBROUTINE outputThorMesh

  SUBROUTINE adjacencyListEntry(element, adjacency_entry, neighbor, current_face,&
                                current_neighbor_face)

    INTEGER, INTENT(IN) :: element
    INTEGER, INTENT(IN) :: adjacency_entry
    INTEGER, INTENT(OUT) :: neighbor
    INTEGER, INTENT(OUT) :: current_face
    INTEGER, INTENT(OUT) :: current_neighbor_face

    INTEGER :: inv_adjacency_entry

    neighbor = adjacency_map(element, adjacency_entry)
    IF (neighbor .EQ. -1) THEN
      neighbor = 0
      current_face = adjacency_entry - 1
      current_neighbor_face = 0
    ELSE
      ! from the adjaceny_map we know that neighbor is the neighbor of element
      ! across face adjacency_entry; now we need to find inv_adjacency_entry
      ! which is the face of neighbor adajcent to element
      inv_adjacency_entry = indexOf(element, adjacency_map(neighbor, :))
      current_face = adjacency_entry - 1
      current_neighbor_face = inv_adjacency_entry - 1
    END IF

  END SUBROUTINE adjacencyListEntry

END MODULE thor_mesh
