!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! THOR MESH GENERATION UTILITY
!   Gmesh module:
!     Handles ingestion of *.gmsh (Gmesh) files.
!     This is intended to enable an *.e (exodus II) to *.thm toolchain for
!     probelm development purposes
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

MODULE gmesh
CONTAINS

  SUBROUTINE ingestGmesh()
    USE globals
    IMPLICIT NONE

    INTEGER:: temp_int_1, temp_int_2
    INTEGER:: i
    REAL:: temp_real_1
    CHARACTER(200):: line = ""


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
    READ(in_unit, *, IOSTAT=io_status) element_count
    IF (io_status .NE. 0) &
         CALL generateErrorMessage(io_status, err_fatal, 'Issue reading element + BC count')

    !Get the true number of elements (Gmesh element count - BC count)
    DO i = 1, element_count
       READ(in_unit, *, IOSTAT=io_status) temp_int_1, temp_int_2
       IF (temp_int_2 .NE. 4) THEN
          IF (temp_int_2 .EQ. 2) EXIT
          IF (temp_int_2 .NE. 2) &
               CALL generateErrorMessage(io_status, err_fatal, 'Invalid shape (non-tet element, non-tri face) found in input')
       END IF
    END DO
    bc_count =  element_count - (temp_int_1 - 1)
    element_count = (temp_int_1 - 1)
    ALLOCATE(element_list(element_count, 4), block_id(element_count))
    ALLOCATE(bc_list(bc_count,3), side_set(bc_count))
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

    !Read in element block id and component nodes
    DO i = 1, element_count
       READ(in_unit, *, IOSTAT=io_status) temp_int_1, temp_int_2, temp_int_2, &
            block_id(i), temp_int_2, temp_int_2, element_list(i,:)
       IF (io_status .NE. 0) &
            CALL generateErrorMessage(io_status, err_fatal, 'Issue reading $Elements block of input')
    END DO

    !Read in boundary condition side set id and component nodes
    DO i = 1, bc_count
       READ(in_unit, *, IOSTAT=io_status) temp_int_1, temp_int_2, temp_int_2, &
          side_set(i), temp_int_2, temp_int_2, bc_list(i,:)
       IF (io_status .NE. 0) &
            CALL generateErrorMessage(io_status, err_fatal, 'Issue reading $Elements block of input during BC ingestion')
    END DO

    CLOSE(in_unit)

  END SUBROUTINE ingestGmesh
END MODULE gmesh
