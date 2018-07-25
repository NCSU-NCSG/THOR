!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! THOR MESH GENERATION UTILITY
!   THOR Mesh Module:
!
!> This module contains the functionality necessary to output a file in the
!! Thor_mesh format (.thrm)
!
!> @author Raffi Yessayan
!> @author Sebastian Schunert
!> @version 1.0
!> @date July, 2017
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
MODULE thor_mesh

  USE globals
  USE tet_mesh_tet_neighbors
  IMPLICIT NONE
CONTAINS

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  !> Updates old to new THOR format
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  SUBROUTINE oldToNewTHORFormat()
    INTEGER :: i

    DO i = 1, bc_count
      !Modify boundary IDs
      side_set(i) = side_set(i) - 1_li

      !Modify face id: this is necessary because the original
      !format shuffled these in a weird way
      IF (boundary_face_list(i) .eq. 1_li) THEN
        boundary_face_list(i) = 2_li
      ELSE IF (boundary_face_list(i) .eq. 2_li) THEN
        boundary_face_list(i) = 0_li
      ELSE IF (boundary_face_list(i) .eq. 3_li) THEN
        boundary_face_list(i) = 1_li
      ELSE IF (boundary_face_list(i) .eq. 4_li) THEN
        boundary_face_list(i) = 3_li
      END IF

    END DO

  END SUBROUTINE

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  !> Extracts the total number of sides from THOR mesh file
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  SUBROUTINE getTHORNumberSides(n_instructions, side_count)
    INTEGER :: n_instructions, side_count
    INTEGER :: n_node, n_elem, local_nside
    INTEGER :: i, j

    OPEN(UNIT = in_unit, file = in_file, ACTION = "read")

    !Header
    READ(in_unit, *) n_node
    READ(in_unit, *) n_elem
    READ(in_unit, *)
    READ(in_unit, *) n_instructions

    DO i = 1, n_node + 2 * n_elem
      READ (in_unit, *)
    END DO

    DO i = 1, n_instructions
      READ(in_unit, *) local_nside
      side_count = side_count + local_nside
      DO j = 1, local_nside
        READ(in_unit, *)
      END DO
    END DO

    CLOSE(in_unit)
  END SUBROUTINE

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  !> Reads a THOR mesh file
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  SUBROUTINE ingestThorMesh()

    INTEGER :: i, j, l, k, s
    INTEGER :: number_side_instructions

    !old THOR format allows more than one side instructions
    !need to get the actual number of side faces
    CALL getTHORNumberSides(number_side_instructions, bc_count)

    !Open File
    OPEN(UNIT = in_unit, file = in_file, ACTION = "read")

    !Header
    READ(in_unit, *) node_count
    READ(in_unit, *) element_count
    !TODO: Find out what these 1 entries are for:
    READ(in_unit, *)
    READ(in_unit, *)

    !Allocate node, element, and block_id array
    ALLOCATE(node_list(node_count,3),element_list(element_count,4),block_id(element_count),&
             source_id(element_count))

    !READ node_list
    DO i = 1, node_count
      READ (in_unit, *) j, node_list(i, :)
    END DO

    !Read block_id and source_id
    DO i = 1, element_count
      READ (in_unit, *) j, block_id(i), source_id(i)
    END DO

    !Read element_list
    DO i = 1, element_count
      READ (in_unit, *) j, element_list(i, :)
    END DO

    !Read bc
    ALLOCATE(boundary_element_list(bc_count), boundary_face_list(bc_count),&
             side_set(bc_count))
    s = 1
    DO i = 1, number_side_instructions
      READ (in_unit, *) j
      DO l = 1, j
        READ (in_unit, *) boundary_element_list(s), boundary_face_list(s), side_set(s)
        s = s + 1
      END DO
    END DO

    !Read adjacency list information
    READ (in_unit, *) j
    ALLOCATE(adjacency_map(element_count,4))
    DO i = 1, 4 * element_count
      READ (in_unit, *) l, s, k, j
      IF (k .eq. 0) THEN
        adjacency_map(l, s + 1) = -1
      ELSE
        adjacency_map(l, s + 1) = k
      END IF
    END DO

    CLOSE(in_unit)

  END SUBROUTINE

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  !> Writes a THOR mesh file
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  SUBROUTINE outputThorMesh()

    INTEGER :: i
    INTEGER :: j
    INTEGER :: position
    INTEGER :: element
    INTEGER :: neighbor
    INTEGER :: current_face
    INTEGER :: current_neighbor_face
    INTEGER :: local_node_list_gmesh(3)
    INTEGER :: order_gmesh(3)
    INTEGER :: local_node_list_adjacency(4)
    INTEGER :: order_adjacency(4)

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

    !if we simply reformat the mesh, then just do a primitive printing of data
    !no additional work is necessary
    IF (execution_mode .eq. 4_li) THEN
      WRITE(out_unit, '(I0)') bc_count
      DO j = 1, bc_count
        WRITE(out_unit, '(I0, 2(X,I0))')  boundary_element_list(j),  boundary_face_list(j),&
                                          side_set(j)
      END DO
    ELSE
      ! If boundary conditions are determined from adjacency go into first branch
      IF (boundary_condition_from_adjacency) THEN
        WRITE(out_unit, '(I0)') SIZE(boundary_element_list)
        DO j = 1, SIZE(boundary_element_list)
          element = boundary_element_list(j)
          current_face = boundary_face_list(j)
          WRITE(out_unit, '(I0, 2(X,I0))') element, current_face, single_boundary_condition_type
        END DO
      ELSE
        !# Of BC Definitions
        WRITE(out_unit, '(I0)') bc_count

        ! TODO: this search can be costly for many bcs; in this case we need to
        ! come up with a better lookup [consider string-hashes]
        IF (SIZE(boundary_element_list) .NE. bc_count) &
              CALL generateErrorMessage(err_code_default, err_fatal, &
              'bc count from gmesh file and adjacency list are inconsistent')

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
              position = mapIndexOf(side_set(i), boundary_id_map(:, 1))
              WRITE(out_unit, '(I0, 2(X,I0))') element, current_face, boundary_id_map(position, 2)
            END IF
          END DO
        END DO
      END IF
    END IF

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

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  !> Computes entries in the Thor adjacency list format from the adjacency_map
  !! Thor format is element - face : neighbor - neighbor_face
  !
  !! @param element element operated upon
  !! @param adjacency_entry index of adacent element
  !! @param neighbor neighboring element
  !! @param current_face face of element being operated upon
  !! @current_neighbor_face face of neighbor being compared against
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
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
