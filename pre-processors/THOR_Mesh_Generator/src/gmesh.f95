!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! THOR MESH GENERATION UTILITY
!   Gmesh Module:
!
!> This module contains the functionality necessary to ingest a file in the
!! Gmesh 3.0 format (.msh)
!
!> @author Raffi Yessayan
!> @author Sebastian Schunert
!> @version 1.0
!> @date July, 2017
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
MODULE gmesh
CONTAINS

  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------
  !> Extracts elements, node, bc, and block_id data from the gmesh file
  !-------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------
  SUBROUTINE ingestGmesh()
    USE globals
    IMPLICIT NONE

    INTEGER :: temp_int_1
    INTEGER :: temp_int_2
    INTEGER :: i, j
    REAL :: temp_real_1
    CHARACTER(200):: line = ""

    INTEGER :: element_instruction_count
    INTEGER :: boundary_instruction_count
    INTEGER :: total_instruction_count
    INTEGER :: elem_read_buffer(8)
    INTEGER :: boundary_read_buffer(4)
    INTEGER :: hex_to_tet_split(4, 6)
    INTEGER :: quad_to_tri_split(3, 2)
    INTEGER :: prism_to_tet_split(4, 3)
    INTEGER :: counter, counterElem, counterBC
    INTEGER :: num_msh_data
    INTEGER :: temp_id

    ! when splitting hex we do not know if we get 5 or 6 tets
    ! hence we need to buffer tets first and then discard dummies because
    ! at read time we don't know how many hexes split into 5 tets
    INTEGER :: temporary_element_count
    INTEGER, ALLOCATABLE :: temporary_element_list(:,:)
    INTEGER, ALLOCATABLE :: temporary_block_id(:)
    INTEGER, ALLOCATABLE :: msh_data_entries(:)

    !Open the gmesh file
    OPEN(UNIT = in_unit, FILE = TRIM(in_file), ACTION = 'READ', STATUS = 'OLD', IOSTAT = io_status)

    !Move to $MeshFormat block
    DO WHILE (TRIM(ADJUSTL(line)) .NE. "$MeshFormat")
      READ(in_unit, *, IOSTAT=io_status) line
      IF (io_status .NE. 0) &
            CALL generateErrorMessage(io_status, err_fatal, 'Issue finding $MeshFormat block of input')
    END DO

    !Ensure float size is 8
    !NOTE:: Should check if float size is sizeof(system float)
    READ(in_unit, *, IOSTAT=io_status) temp_real_1, temp_int_1, temp_int_2
    IF (temp_int_2 .NE. 8) &
          CALL generateErrorMessage(io_status, err_fatal, 'System float width does not match input file')

    !Move to $Nodes block
    !NOTE:: Node numbering must be sequential in this naive implementation
    DO WHILE (TRIM(ADJUSTL(line)) .NE. "$Nodes")
      READ(in_unit, *, IOSTAT=io_status) line
      IF (io_status .NE. 0) &
            CALL generateErrorMessage(io_status, err_fatal, 'Issue finding $Nodes block of input')
    END DO

    !Get node count & allocate
    READ(in_unit, *, IOSTAT=io_status) node_count
    ALLOCATE(node_list(node_count, 3))

    !Read in node <x,y,z> positions
    DO i = 1, node_count
      READ(in_unit, *, IOSTAT=io_status) temp_int_1, node_list(i, :)
      IF (io_status .NE. 0) &
            CALL generateErrorMessage(io_status, err_fatal, 'Issue reading $Nodes block of input')
    END DO

    !Move to $Elements block
    !NOTE:: Element numbering must be sequential in this naive implementation
    !NOTE:: BC data must follow element data
    DO WHILE (TRIM(ADJUSTL(line)) .NE. "$Elements")
      READ(in_unit, *, IOSTAT=io_status) line
      IF (io_status .NE. 0) &
            CALL generateErrorMessage(io_status, err_fatal, 'Issue finding $Elements block of input')
    END DO
    READ(in_unit, *, IOSTAT=io_status) total_instruction_count
    IF (io_status .NE. 0) &
          CALL generateErrorMessage(io_status, err_fatal, 'Issue reading element + BC count')

    !Get the true number of elements by looking at each element instruction
    !accounting for the need to split elements & boundaries
    !This is necessary because GMSH considers boundary faces as d-1 dimensional
    !elements
    element_instruction_count = 0
    boundary_instruction_count = 0
    temporary_element_count = 0
    bc_count = 0
    DO i = 1, total_instruction_count
      READ(in_unit, *, IOSTAT=io_status) temp_int_1, temp_int_2, num_msh_data
      IF (temp_int_2 .EQ. 4) THEN
        ! tetrahedron
        element_instruction_count = element_instruction_count + 1
        temporary_element_count = temporary_element_count + 1
      ELSE IF (temp_int_2 .EQ. 5) THEN
        ! hexahedron
        element_instruction_count = element_instruction_count + 1
        temporary_element_count = temporary_element_count + 6
      ELSE IF (temp_int_2 .EQ. 6) THEN
        ! triangular prism
        element_instruction_count = element_instruction_count + 1
        temporary_element_count = temporary_element_count + 3
      ELSE IF (temp_int_2 .EQ. 2) THEN
        ! triangle [boundary]
        bc_count = bc_count + 1
        boundary_instruction_count = boundary_instruction_count + 1
      ELSE IF (temp_int_2 .EQ. 3) THEN
        ! quadrilaterial [boundary]
        bc_count = bc_count + 2
        boundary_instruction_count = boundary_instruction_count + 1
      ELSE
        CALL generateErrorMessage(io_status, err_fatal,&
            'Invalid shape: mesh converter allows tet, tri-prism, hex [volume] and triangle, quadrilaterial [boundary]')
      END IF
    END DO
    ALLOCATE(temporary_element_list(temporary_element_count, 4), temporary_block_id(temporary_element_count))
    ALLOCATE(bc_list(bc_count,3), side_set(bc_count))
    ALLOCATE(msh_data_entries(num_msh_data))
    line = ""

    !Return to start of $Elements block
    REWIND(in_unit)
    DO WHILE (TRIM(ADJUSTL(line)) .NE. "$Elements")
      READ(in_unit, *, IOSTAT=io_status) line
      IF (io_status .NE. 0) &
            CALL generateErrorMessage(io_status, err_fatal, 'Issue finding $Elements block of input')
    END DO
    READ(in_unit, *, IOSTAT=io_status)
    IF (io_status .NE. 0) &
          CALL generateErrorMessage(io_status, err_fatal, 'Issue parsing input file')

    !Read in elements and BCs
    counterElem = 1
    counterBC = 1
    DO i = 1, element_instruction_count+boundary_instruction_count
      READ(in_unit, '(A200)') line
      READ(line, *, IOSTAT=io_status) temp_int_1, temp_int_2, temp_int_1, temp_id

      !Make sure this is a valid element shape and the read was successful
      IF (io_status .NE. 0) &
            CALL generateErrorMessage(io_status, err_fatal, 'Issue reading $Elements block of input')
      IF (temp_int_2 .NE. 2 .AND. temp_int_2 .NE. 3 .AND. temp_int_2 .NE. 4 .AND. temp_int_2 .NE. 5 .AND. temp_int_2 .NE. 6) &
            CALL generateErrorMessage(err_code_default, err_fatal, 'Invalid Element Shape')

      !For 2-D faces
      IF (temp_int_2 .EQ. 2) THEN
        READ(line, *, IOSTAT=io_status) temp_int_1, temp_int_1, temp_int_1, msh_data_entries, bc_list(counterBC,:)
        IF (io_status .NE. 0) &
            CALL generateErrorMessage(io_status, err_fatal, &
            'Issue reading $Elements block of input during BC ingestion processing triangle')
        side_set(counterBC) = temp_id
        counterBC = counterBC + 1
      ELSE IF (temp_int_2 .EQ. 3) THEN
        READ(line, *, IOSTAT=io_status) temp_int_1, temp_int_1, temp_int_1, msh_data_entries, boundary_read_buffer
        IF (io_status .NE. 0) &
            CALL generateErrorMessage(io_status, err_fatal, &
            'Issue reading $Elements block of input during BC ingestion processing quadrilateral')
        CALL splitQuadIntoTri(boundary_read_buffer, quad_to_tri_split)
        DO j = 1, 2
          bc_list(counterBC,:) = quad_to_tri_split(:, j)
          side_set(counterBC) = temp_id
          counterBC = counterBC + 1
        END DO
      END IF

      !For volumetric elements
      IF (temp_int_2 .EQ. 4) THEN
        ! tetrahdron -> just save it
        READ(line, *, IOSTAT=io_status) temp_int_1, temp_int_1, temp_int_1, msh_data_entries, temporary_element_list(counterElem, :)
        IF (io_status .NE. 0) &
            CALL generateErrorMessage(io_status, err_fatal, 'Issue reading $Elements when processing tet info')
        temporary_block_id(counterElem) = temp_id
        counterElem = counterElem + 1
      ELSE IF (temp_int_2 .EQ. 5) THEN
        ! hexahedron splits into 6 tets
        READ(line, *, IOSTAT=io_status) temp_int_1, temp_int_1, temp_int_1, msh_data_entries, elem_read_buffer
        IF (io_status .NE. 0) &
            CALL generateErrorMessage(io_status, err_fatal, 'Issue reading $Elements when processing hex')
        CALL splitHexIntoTet(elem_read_buffer, hex_to_tet_split)
        DO j = 1, 6
          temporary_element_list(counterElem, :) = hex_to_tet_split(:, j)
          temporary_block_id(counterElem) = temp_id
          counterElem = counterElem + 1
        END DO
      ELSE IF (temp_int_2 .EQ. 6) THEN
        ! triangular prism splits into 3 tets
        READ(line, *, IOSTAT=io_status) temp_int_1, temp_int_1, temp_int_1, msh_data_entries, elem_read_buffer(1:6)
        IF (io_status .NE. 0) &
            CALL generateErrorMessage(io_status, err_fatal, 'Issue reading $Elements when processing prism')
        CALL splitPrismToTet(elem_read_buffer(1:6), prism_to_tet_split)
        DO j = 1, 3
          temporary_element_list(counterElem, :) = prism_to_tet_split(:, j)
          temporary_block_id(counterElem) = temp_id
          counterElem = counterElem + 1
        END DO
      END IF
    END DO
    CLOSE(in_unit)

    ! Finally we need to set the permanent block_id and element_list
    element_count = 0
    DO i = 1, temporary_element_count
      IF (temporary_element_list(i, 1) .NE. -1) element_count = element_count + 1
    END DO
    ALLOCATE(element_list(element_count, 4), block_id(element_count))
    counter = 1
    DO i = 1, temporary_element_count
      IF (temporary_element_list(i, 1) .NE. -1) THEN
        element_list(counter, :) = temporary_element_list(i, :)
        block_id(counter) = temporary_block_id(i)
        counter = counter + 1
      END IF
    END DO
    DEALLOCATE(temporary_element_list, temporary_block_id)

  END SUBROUTINE ingestGmesh
END MODULE gmesh
