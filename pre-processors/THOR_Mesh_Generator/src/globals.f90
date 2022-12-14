!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! THOR MESH GENERATION UTILITY
!   Globals Module:
!
!>    This module contains variables and functions used to store and manipulate
!!    data common to both the input mesh and output mesh.
!
!> @author Raffi Yessayan
!> @author Sebastian Schunert
!> @version 1.0
!> @date July, 2017
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
MODULE globals
  USE integer_array_tools
  USE tet_mesh_tet_neighbors
  IMPLICIT NONE

  !> Unit number opened for primary input file
  INTEGER :: std_in_unit = 19
  !> Unit number supplied for mesh input file
  INTEGER :: in_unit = 20
  !> Unit number supplied for mesh output file
  INTEGER :: out_unit = 21
  !> Unit number supplied for region id input file
  INTEGER :: region_id_unit = 22
  !> Unit number supplied for source id input file
  INTEGER :: source_id_unit = 23
  !> Unit number supplied for boundary id input file
  INTEGER :: boundary_id_unit = 24
  !> Mesh refinement level
  INTEGER :: uniform_refine = 0
  !> Used to track I/O operation status codes
  INTEGER :: io_status = 0
  !> Default error code when generating an error message without a more specific
  !! error type
  INTEGER :: err_code_default = -1
  !> Error severity level requiring output to user, but no other action
  INTEGER :: err_warning = 0
  !> Error severity level requiring immediate termination of program
  INTEGER :: err_fatal = 1
  !> execution_mode: 1 - gmsh, 2 - exodus, 3 - unv, 4 - old THOR mesh
  INTEGER :: execution_mode
  !> Number of <x,y,z> nodes in the mesh
  INTEGER :: node_count
  !> Number of elements in the mesh
  INTEGER :: element_count
  !> Number of boundary conditions faces defined in the mesh
  INTEGER :: bc_count
  !> Number of block ids defined in the mesh
  INTEGER :: nblocks
  !> Number of side set ids defined in the mesh
  INTEGER :: n_side_sets
  !> Paramter defined int, width 8
  INTEGER, PARAMETER :: li = SELECTED_INT_KIND(8)
  !> Parameter defined double precision real
  INTEGER, PARAMETER :: d_t = SELECTED_REAL_KIND(15,307)
  !> Mapping from input mesh block ids to output mesh block ids
  !! Defaults to a -> a
  INTEGER, ALLOCATABLE :: block_id_map(:,:)
  !> Map giving the 4 elements adjacent to each element
  !! If no neighbor exists on a face, cell face is boundary
  !! and has adjacency = -1
  INTEGER, ALLOCATABLE :: adjacency_map(:,:)
  !> Mapping from input source ids to output source ids
  !! Defaults to a -> a
  INTEGER, ALLOCATABLE :: source_id_map(:,:)
  !> Mapping from input sideset IDS to output boundary condition types
  !! Defaults to a -> 0
  INTEGER, ALLOCATABLE :: boundary_id_map(:,:)
  !> List of elements with a face on a boundary
  !! Element is repeated for each boundary face
  INTEGER, ALLOCATABLE :: boundary_element_list(:)
  !> List containing the face number corresponding to the element with
  !! the same index in boundary_element_list
  INTEGER, ALLOCATABLE :: boundary_face_list(:)
  !> A list of the 4 nodes composing each tetrahedral element
  INTEGER, ALLOCATABLE :: element_list(:,:)
  !>The mapping of element to block id from the input file
  INTEGER, ALLOCATABLE :: block_id(:)
  !>The mapping of element to source id from the input file
  INTEGER, ALLOCATABLE :: source_id(:)
  !> The nodes that make up each boundary face
  INTEGER, ALLOCATABLE :: bc_list(:,:)
  !> The boundary condition assigned to the boundary face
  INTEGER, ALLOCATABLE :: side_set(:)
  !> The <x,y,z> position of each node
  REAL, ALLOCATABLE :: node_list(:,:)
  !> Used to skip user driven region mapping if input file not found or
  !! otherwise disabled
  LOGICAL :: skip_region_map
  !> Used to skip user driven source mapping if input file not found or
  !! otherwise disabled
  LOGICAL :: skip_source_map
  !> Used to skip user driven boundary conditions mapping if input file not
  !! found or otherwise disabled
  LOGICAL :: skip_boundary_map
  !> If true boundary conditions are all of a single type and determined by
  !! by adjacency list
  LOGICAL :: boundary_condition_from_adjacency
  !> If boundary_condition_from_adjacency is true this boundary condition type
  !! is used for all exterior faces determined by adjacency list
  INTEGER :: single_boundary_condition_type
  !> Holds the name of the standard input config file
  CHARACTER(200) :: std_in_file = ""
  !> Holds the name of the standard input mesh file
  CHARACTER(200) :: in_file = ""
  !> Holds the original name of the standard input mesh file; can be changed during execution
  CHARACTER(200) :: original_in_file = ""
  !> Holds the name of the standard output mesh file
  CHARACTER(200) :: out_file = ""
  !> Holds the name of the standard input region id file
  CHARACTER(200) :: region_id_file = ""
  !> Holds the name of the standard input source id file
  CHARACTER(200) :: source_id_file = ""
  !> Holds the name of the standard input boundary id file
  CHARACTER(200) :: boundary_id_file = ""

CONTAINS

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  !> Used to check if an I/O operation returned a specific status.
  !! If the status is unexpected, generates a fatal error
  !!
  !! @param status Status returned by I/O operation
  !! @param expected Expected status from I/O operation
  !! @param message Message for generateErrorMessage if actual != expected
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  SUBROUTINE checkStatus(status, expected, message)
    INTEGER, INTENT(IN):: status, expected
    CHARACTER(*), INTENT(IN):: message
    IF (status .NE. expected) &
          CALL generateErrorMessage(err_code_default, err_fatal, message)
  END SUBROUTINE checkStatus

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  !> Used to check if a given filename represents a file on the system.
  !! If it does not, generates a fatal error message via generateErrorMessage.
  !!
  !! @param filename filename to check for existance
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  SUBROUTINE checkFileExists(filename)
    CHARACTER(len=*), INTENT(in):: filename
    LOGICAL:: file_exists
    INQUIRE(FILE = TRIM(filename), EXIST = file_exists)
    IF (.NOT. file_exists) CALL generateErrorMessage(err_code_default, err_fatal, 'File does not exist')
  END SUBROUTINE checkFileExists

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  !> Used to generate error messages. Allows for differing actions based on
  !! severity of error.
  !!
  !! @param error_code The error code generated by the action triggering the
  !!        error
  !! @param severity Defines the severity of the error. See err_fatal, etc.
  !! @param message Message for user explaining error cause or location
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  SUBROUTINE generateErrorMessage(error_code, severity, message)
    INTEGER, INTENT(IN):: error_code, severity
    CHARACTER(*), INTENT(IN):: message
    WRITE(*,'(A)') '============================================================'
    WRITE(*,'(A)')'AN UNEXPECTED ERROR OCCURRED'
    WRITE(*,'(A,I0)') 'Error Code: ', error_code
    WRITE(*,'(A)')'Error Message: ', TRIM(message)
    WRITE(*,'(A)') '============================================================'
    IF (severity .EQ. err_fatal) STOP
  END SUBROUTINE generateErrorMessage

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  !> Opens a filename and assigns it to unit. Then reads the first entry to
  !! extract the length of the file
  !!
  !! @param filename name of file to be opened
  !! @param unit unit number to be assigned to file
  !! @return The number of entries contained in the file exclusive of the
  !!                 entry count
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  INTEGER FUNCTION lengthTextFile(filename, unit)
    CHARACTER(len=*), INTENT(IN):: filename
    INTEGER, INTENT(IN):: unit

    INTEGER :: status

    CALL checkFileExists(filename)
    OPEN(UNIT = unit, FILE = TRIM(filename), ACTION = 'READ', STATUS = 'OLD',&
          IOSTAT = status)
    CALL checkStatus(status, 0, 'File cannot be opened')
    READ(unit, *, IOSTAT = status) lengthTextFile
    CALL checkStatus(status, 0, 'Reading the number of entries in file failed')
    CLOSE(unit)
    RETURN
  END FUNCTION lengthTextFile

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  !> Extends the functionality of lengthTextFile by allowing N entries to be
  !! read in into an INT array in place of simply returning the number of
  !! entries.
  !! Still expects the file to have num_entries on the first line
  !!
  !! @param filename name of file to be opened
  !! @param unit unit number to be assigned to file
  !! @param array array into which to load data
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  SUBROUTINE loadTextFile(filename, unit, array)
    CHARACTER(len=*), INTENT(IN):: filename
    INTEGER, INTENT(IN):: unit
    INTEGER, DIMENSION (:, :), INTENT(OUT) :: array

    INTEGER :: j, s1, nentries
    INTEGER :: status

    s1 = SIZE(array(:, 1))

    CALL checkFileExists(filename)
    OPEN(UNIT = unit, FILE = TRIM(filename), ACTION = 'READ', STATUS = 'OLD', IOSTAT = status)
    CALL checkStatus(status, 0, 'File cannot be opened')
    READ(unit, *, IOSTAT = status) nentries
    CALL checkStatus(status, 0, 'Reading the number of entries in file failed')
    CALL checkStatus(nentries, s1, 'Number of file entries is different than allocate array size')
    DO j = 1, s1
      READ(unit, *, IOSTAT = status) array(j, :)
      CALL checkStatus(status, 0, 'Reading instructions in file failed')
    END DO
    CLOSE(unit)
  END SUBROUTINE loadTextFile

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  !> Ingests the standard user input to identify mesh in, mesh out, region
  !! mapping, and source mapping files.
  !! If the source or region mapping files are missing, sets skip_region_map
  !! and skip_source_map appropriately
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  SUBROUTINE ingestInput()
    ! open the standard input file
    CHARACTER(len = 200):: line
    LOGICAL :: temp_logical

    CALL checkFileExists(std_in_file)
    OPEN(UNIT = std_in_unit, FILE = TRIM(std_in_file), ACTION = 'READ', STATUS = 'OLD', IOSTAT = io_status)
    CALL checkStatus(io_status, 0, 'Standard input file cannot be opened')

    ! mesh infile
    READ(std_in_unit, '(A200)', IOSTAT = io_status) line
    CALL checkStatus(io_status, 0, 'Reading input mesh file name failed')
    in_file = TRIM(line)

    ! mesh outfile
    READ(std_in_unit, '(A200)', IOSTAT = io_status) line
    CALL checkStatus(io_status, 0, 'Reading output mesh file name failed')
    out_file = TRIM(line)

    ! region id mapping file
    READ(std_in_unit, '(A200)', IOSTAT = io_status) line
    temp_logical = .FALSE.
    IF (io_status .EQ. 0) &
          INQUIRE(FILE = TRIM(ADJUSTL(line)), EXIST = temp_logical)
    ! region id mapping can be skipped
    IF (.NOT. temp_logical) THEN
      skip_region_map = .TRUE.
    ELSE
      skip_region_map = .FALSE.
      region_id_file = TRIM(ADJUSTL(line))
    END IF

    ! src id mapping file
    READ(std_in_unit, '(A200)', IOSTAT = io_status) line
    temp_logical = .FALSE.
    IF (io_status .EQ. 0) &
          INQUIRE(FILE = TRIM(ADJUSTL(line)), EXIST = temp_logical)
    ! region id mapping can be skipped
    IF (.NOT. temp_logical) THEN
      skip_source_map = .TRUE.
    ELSE
      skip_source_map = .FALSE.
      source_id_file = TRIM(ADJUSTL(line))
    END IF

    READ(std_in_unit, '(A200)', IOSTAT = io_status) line
    temp_logical = .FALSE.
    IF (io_status .EQ. 0) &
          INQUIRE(FILE = TRIM(ADJUSTL(line)), EXIST = temp_logical)
    ! region id mapping can be skipped
    IF (.NOT. temp_logical) THEN
      skip_boundary_map = .TRUE.
    ELSE
      skip_boundary_map = .FALSE.
      boundary_id_file = TRIM(ADJUSTL(line))
    END IF
    CLOSE(std_in_unit)
  END SUBROUTINE ingestInput

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  !> Establishes default source and region mapping from mesh input and then,
  !! if enabled, overrides these with user specified instructions from the
  !! appropriate input files
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  SUBROUTINE setupOptionalReMapping()

    INTEGER :: j, position
    INTEGER, DIMENSION(:), ALLOCATABLE :: unique_ids
    INTEGER, DIMENSION(:), ALLOCATABLE :: order
    INTEGER :: nentries_region, nentries_source, nentries_bc
    INTEGER, DIMENSION(:, :), ALLOCATABLE :: region_instructions
    INTEGER, DIMENSION(:, :), ALLOCATABLE :: source_instructions
    INTEGER, DIMENSION(:, :), ALLOCATABLE :: boundary_instructions

    ! Set the default mapping
    nblocks = numUniqueEntries(block_id)
    n_side_sets = numUniqueEntries(side_set)

    ALLOCATE(block_id_map(nblocks, 2), source_id_map(nblocks, 2), boundary_id_map(n_side_sets, 2))

    ! sort unique ids to make block_id_map and source_id_map sorted
    ALLOCATE(unique_ids(nblocks), order(nblocks))
    CALL uniqueEntries(block_id, unique_ids)
    CALL quickSortInteger(unique_ids, order)
    ! assign defaults to block_id_map and source_id_map
    DO j = 1, nblocks
      block_id_map(j, :) = unique_ids(j)
      source_id_map(j, :) = unique_ids(j)
    END DO

    ! assign defaults to boundary_id_map
    DEALLOCATE(unique_ids, order)
    ALLOCATE(unique_ids(n_side_sets), order(n_side_sets))
    CALL uniqueEntries(side_set, unique_ids)
    CALL quickSortInteger(unique_ids, order)
    DO j = 1, n_side_sets
      boundary_id_map(j, 1) = unique_ids(j)
      boundary_id_map(j, 2) = 0
    END DO

    ! Read region id and source id from files and override block_id_map default
    IF (.NOT. skip_region_map) THEN
      nentries_region = lengthTextFile(region_id_file, region_id_unit)
      ALLOCATE(region_instructions(nentries_region, 2))
      CALL loadTextFile(region_id_file, region_id_unit, region_instructions)

      ! work instructions into block_id_map
      IF (.NOT. hasUniqueEntries(region_instructions(:, 1))) &
            CALL generateErrorMessage(err_code_default, err_fatal, &
            'regions instruction keys are not unique')

      DO j = 1, nentries_region
        position = mapIndexOf(region_instructions(j, 1), block_id_map(:, 1))
        block_id_map(position, 2) = region_instructions(j, 2)
      END DO
    END IF

    ! Read region id and source id from files and override block_id_map default
    IF (.NOT. skip_source_map) THEN
      nentries_source = lengthTextFile(source_id_file, source_id_unit)
      ALLOCATE(source_instructions(nentries_source, 2))
      CALL loadTextFile(source_id_file, source_id_unit, source_instructions)
      ! work instructions into source_id_map
      IF (.NOT. hasUniqueEntries(source_instructions(:, 1))) &
            CALL generateErrorMessage(err_code_default, err_fatal, &
            'source instruction keys are not unique')

      DO j = 1, nentries_source
        position = mapIndexOf(source_instructions(j, 1), source_id_map(:, 1))
        source_id_map(position, 2) = source_instructions(j, 2)
      END DO
    END IF

    ! Read boundary id and override block_id_map default
    boundary_condition_from_adjacency = .FALSE.
    IF (.NOT. skip_boundary_map) THEN
      nentries_bc = lengthTextFile(boundary_id_file, boundary_id_unit)
      ALLOCATE(boundary_instructions(nentries_bc, 2))
      CALL loadTextFile(boundary_id_file, boundary_id_unit, boundary_instructions)
      ! work instructions into source_id_map
      IF (.NOT. hasUniqueEntries(boundary_instructions(:, 1))) &
            CALL generateErrorMessage(err_code_default, err_fatal, &
            'boundary instruction keys are not unique')

      DO j = 1, nentries_bc
        position = mapIndexOf(boundary_instructions(j, 1), boundary_id_map(:, 1))
        boundary_id_map(position, 2) = boundary_instructions(j, 2)
        ! FIXME: currently if any of the keys is < 0 we set the boundary_condition_from_adjacency
        ! to true and get single_boundary_condition_type. That is not really logical but the easist
        ! until better input handling is in place.
        IF (boundary_instructions(j, 1) < 0) THEN
          boundary_condition_from_adjacency = .TRUE.
          single_boundary_condition_type = boundary_instructions(j, 2)
        END IF
      END DO
    END IF

    DEALLOCATE(unique_ids, order)
    IF (ALLOCATED(region_instructions)) DEALLOCATE(region_instructions)
    IF (ALLOCATED(source_instructions)) DEALLOCATE(source_instructions)
    IF (ALLOCATED(boundary_instructions)) DEALLOCATE(boundary_instructions)
  END SUBROUTINE setupOptionalReMapping

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  !> Populates the tet_mesh adjacency_map by calling out to the contrib
  !! subroutine tet_mesh_neighbor_tets in tet_mesh_tet_neighbors
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  SUBROUTINE computeAdjacencyList()
    ALLOCATE(adjacency_map(element_count, 4))
    CALL tet_mesh_neighbor_tets(4, element_count, TRANSPOSE(element_list),&
          TRANSPOSE(adjacency_map))
  END SUBROUTINE computeAdjacencyList

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  !> Given an adjacency map has been established, returns a pair of arrays
  !! detailing boundary faces and elements.
  !!
  !! See boundary_face_list and boundary_element_list
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  SUBROUTINE getBoundaryElements()
    INTEGER :: i
    INTEGER :: j
    INTEGER :: counter
    ALLOCATE(boundary_element_list(bc_count), boundary_face_list(bc_count))
    counter = 1
    DO i = 1, element_count
      DO j = 1, 4
        IF (adjacency_map(i, j) .EQ. -1) THEN
          ! check if we, by error, exceeded bc_count
          IF (counter .GT. bc_count .AND. .NOT. boundary_condition_from_adjacency)  &
                CALL generateErrorMessage(err_code_default, err_fatal, &
                'Boundary elements and adjacency list are inconsistent: too many -1 in adjacency list')
          boundary_element_list(counter) = i
          boundary_face_list(counter) = j - 1
          counter = counter + 1
        END IF
      END DO
    END DO
    IF (counter - 1 .NE. bc_count  .AND. .NOT. boundary_condition_from_adjacency)  &
          CALL generateErrorMessage(err_code_default, err_fatal, &
          'Boundary elements and adjacency list are inconsistent: not enough -1 in adjacency list')
  END SUBROUTINE getBoundaryElements

  SUBROUTINE printProgramHeader()
    ! TODO: A unified program name header?
    WRITE(6, '(A)') "Program Name: THOR_MESH_GENERATOR"
  END SUBROUTINE printProgramHeader

  SUBROUTINE echoIngestedInput()
    INTEGER :: j
    INTEGER :: s

    WRITE(6, *)
    WRITE(6, '(A)') "------------------------------------ Input ------------------------------------"
    WRITE(6, '(A)') "Infile mesh file name: " // TRIM(original_in_file)
    WRITE(6, '(A)') "Output mesh file name: " // TRIM(out_file)
    WRITE(6, '(A23, I4)') "Mesh refinement level: ", uniform_refine

    WRITE(6, *)
    IF (skip_region_map) THEN
      WRITE(6, '(A)') "No region ID edits are provided"
    ELSE
      WRITE(6, '(A, A)') "Reporting Region ID reassignments from file: ", TRIM(region_id_file)
      IF (.NOT. ALLOCATED(block_id_map)) &
            CALL generateErrorMessage(err_code_default, err_fatal, &
            'block id map not allocated,  echoIngestedInput called too early')
      s = SIZE(block_id_map(:, 1))
      DO j = 1, s
        WRITE(6, '(A,I0,A,I0)') "Old Region ID ", block_id_map(j, 1), " New Region ID ",&
              block_id_map(j, 2)
      END DO
    END IF

    IF (skip_source_map) THEN
      WRITE(6, '(A)') "No source ID edits are provided"
    ELSE
      IF (.NOT. ALLOCATED(source_id_map)) &
            CALL generateErrorMessage(err_code_default, err_fatal, &
            'source id map not allocated,  echoIngestedInput called too early')
      WRITE(6, '(A,A)') "Reporting Source ID reassignments from file: ", TRIM(source_id_file)
      s = SIZE(source_id_map(:, 1))
      DO j = 1, s
        WRITE(6, '(A,I0,A,I0)') "Old Region ID ", source_id_map(j, 1), " Source ID ",&
              source_id_map(j, 2)
      END DO
    END IF
    IF (boundary_condition_from_adjacency) THEN
      WRITE(6, '(A)') "Boundary ID assignment is overriden. Boundaries determined by adjacency list."
      WRITE(6, '(A,I4)') "All exterior boundaries are assigned boundary type ", single_boundary_condition_type
    ELSE
      IF (skip_boundary_map) THEN
        WRITE(6, '(A)') "No boundary ID edits are provided"
      ELSE
        IF (.NOT. ALLOCATED(boundary_id_map)) &
              CALL generateErrorMessage(err_code_default, err_fatal, &
              'boundary id map not allocated,  echoIngestedInput called too early')
        WRITE(6, '(A,A)') "Reporting Source Boundary reassignments from file: ", TRIM(boundary_id_file)
        s = SIZE(boundary_id_map(:, 1))
        DO j = 1, s
          WRITE(6, '(A,I0,A,I0)') "Sideset ID ", boundary_id_map(j, 1), " Boundary Type ",&
                boundary_id_map(j, 2)
        END DO
      END IF
    END IF
    WRITE(6, *)
  END SUBROUTINE echoIngestedInput

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  !> Splits a hexahedron into six tetrahedra
  !! Algorithm described in: HOW TO SUBDIVIDE PYRAMIDS, PRISMS
  !! AND HEXAHEDRA INTO TETRAHEDRA. Dompierre, J. et al.
  !! NOTE: local hex ordering required in paper identical to GMESH ordering
  !!
  !! @param vertex_id a vector of length 8 containing the global ids of the vertices
  !! @param tet_split 4x6 array containing the 4 global ids of the 6 tetrahedra
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  SUBROUTINE splitHexIntoTet(vertex_id, tet_split)
    INTEGER, INTENT(in) :: vertex_id(8)
    INTEGER, INTENT(inout) :: tet_split(4, 6)

    ! Define local variables
    INTEGER :: smallest_pos, smallest_value, i, n_diags, bits(3), temp
    INTEGER :: table(8)

    ! Find smallest vertex id
    smallest_pos = 1
    smallest_value = vertex_id(1)
    DO i = 2, 8
      IF (smallest_value .GT. vertex_id(i)) THEN
        smallest_value = vertex_id(i)
        smallest_pos = i
      END IF
    END DO

    IF (smallest_pos .EQ. 1) THEN
      table = (/1,2,3,4,5,6,7,8/)
    ELSE IF (smallest_pos .EQ. 2) THEN
      table = (/2,1,5,6,3,4,8,7/)
    ELSE IF (smallest_pos .EQ. 3) THEN
      table = (/3,2,6,7,4,1,5,8/)
    ELSE IF (smallest_pos .EQ. 4) THEN
      table = (/4,1,2,3,8,5,6,7/)
    ELSE IF (smallest_pos .EQ. 5) THEN
      table = (/5,1,4,8,6,2,3,7/)
    ELSE IF (smallest_pos .EQ. 6) THEN
      table = (/6,2,1,5,7,3,4,8/)
    ELSE IF (smallest_pos .EQ. 7) THEN
      table = (/7,3,2,6,8,4,1,5/)
    ELSE IF (smallest_pos .EQ. 8) THEN
      table = (/8,4,3,7,5,1,2,6/)
    END IF

    ! Compute number of diagonal through vertex_id(table(7))
    n_diags = 0
    bits = 0
    IF (min(vertex_id(table(2)), vertex_id(table(7))) < min(vertex_id(table(3)), vertex_id(table(6)))) THEN
      n_diags = n_diags + 1
      bits(1) = 1
    END IF
    IF (min(vertex_id(table(4)), vertex_id(table(7))) < min(vertex_id(table(3)), vertex_id(table(8)))) THEN
      n_diags = n_diags + 1
      bits(1) = 2
    END IF
    IF (min(vertex_id(table(5)), vertex_id(table(7))) < min(vertex_id(table(6)), vertex_id(table(8)))) THEN
      n_diags = n_diags + 1
      bits(3) = 1
    END IF

    ! Rotation
    IF ((bits(1) .EQ. 0 .AND. bits(2) .EQ. 0 .AND. bits(3) .EQ. 1) .OR. &
        (bits(1) .EQ. 1 .AND. bits(2) .EQ. 1 .AND. bits(3) .EQ. 0)) THEN
      ! Rotation by 120 degrees
      temp = table(2)
      table(2) = table(5)
      table(5) = table(4)
      table(4) = temp

      temp = table(6)
      table(6) = table(8)
      table(8) = table(3)
      table(3) = temp
    ELSE IF ((bits(1) .EQ. 0 .AND. bits(2) .EQ. 1 .AND. bits(3) .EQ. 0) .OR. &
             (bits(1) .EQ. 1 .AND. bits(2) .EQ. 0 .AND. bits(3) .EQ. 1)) THEN
      ! Rotation by 240 degrees
      temp = table(2)
      table(2) = table(4)
      table(4) = table(5)
      table(5) = temp

      temp = table(6)
      table(6) = table(3)
      table(3) = table(8)
      table(8) = temp
    END IF

    ! now split by n_diags
    IF (n_diags .EQ. 0) THEN
      ! tet 1
      tet_split(1, 1) = vertex_id(table(1))
      tet_split(2, 1) = vertex_id(table(2))
      tet_split(3, 1) = vertex_id(table(3))
      tet_split(4, 1) = vertex_id(table(6))

      ! tet 2
      tet_split(1, 2) = vertex_id(table(1))
      tet_split(2, 2) = vertex_id(table(3))
      tet_split(3, 2) = vertex_id(table(8))
      tet_split(4, 2) = vertex_id(table(6))

      ! tet 3
      tet_split(1, 3) = vertex_id(table(1))
      tet_split(2, 3) = vertex_id(table(3))
      tet_split(3, 3) = vertex_id(table(4))
      tet_split(4, 3) = vertex_id(table(8))

      ! tet 4
      tet_split(1, 4) = vertex_id(table(1))
      tet_split(2, 4) = vertex_id(table(6))
      tet_split(3, 4) = vertex_id(table(8))
      tet_split(4, 4) = vertex_id(table(5))

      ! tet 5
      tet_split(1, 5) = vertex_id(table(3))
      tet_split(2, 5) = vertex_id(table(8))
      tet_split(3, 5) = vertex_id(table(6))
      tet_split(4, 5) = vertex_id(table(7))

      ! tet 6
      tet_split(1, 6) = -1
      tet_split(2, 6) = -1
      tet_split(3, 6) = -1
      tet_split(4, 6) = -1

    ELSE IF (n_diags .EQ. 1) THEN
      ! tet 1
      tet_split(1, 1) = vertex_id(table(1))
      tet_split(2, 1) = vertex_id(table(6))
      tet_split(3, 1) = vertex_id(table(8))
      tet_split(4, 1) = vertex_id(table(5))

      ! tet 2
      tet_split(1, 2) = vertex_id(table(1))
      tet_split(2, 2) = vertex_id(table(2))
      tet_split(3, 2) = vertex_id(table(8))
      tet_split(4, 2) = vertex_id(table(6))

      ! tet 3
      tet_split(1, 3) = vertex_id(table(2))
      tet_split(2, 3) = vertex_id(table(7))
      tet_split(3, 3) = vertex_id(table(8))
      tet_split(4, 3) = vertex_id(table(6))

      ! tet 4
      tet_split(1, 4) = vertex_id(table(1))
      tet_split(2, 4) = vertex_id(table(8))
      tet_split(3, 4) = vertex_id(table(3))
      tet_split(4, 4) = vertex_id(table(4))

      ! tet 5
      tet_split(1, 5) = vertex_id(table(1))
      tet_split(2, 5) = vertex_id(table(8))
      tet_split(3, 5) = vertex_id(table(2))
      tet_split(4, 5) = vertex_id(table(3))

      ! tet 6
      tet_split(1, 6) = vertex_id(table(2))
      tet_split(2, 6) = vertex_id(table(8))
      tet_split(3, 6) = vertex_id(table(7))
      tet_split(4, 6) = vertex_id(table(3))

    ELSE IF (n_diags .EQ. 2) THEN
      ! tet 1
      tet_split(1, 1) = vertex_id(table(1))
      tet_split(2, 1) = vertex_id(table(5))
      tet_split(3, 1) = vertex_id(table(6))
      tet_split(4, 1) = vertex_id(table(7))

      ! tet 2
      tet_split(1, 2) = vertex_id(table(1))
      tet_split(2, 2) = vertex_id(table(4))
      tet_split(3, 2) = vertex_id(table(8))
      tet_split(4, 2) = vertex_id(table(7))

      ! tet 3
      tet_split(1, 3) = vertex_id(table(1))
      tet_split(2, 3) = vertex_id(table(8))
      tet_split(3, 3) = vertex_id(table(5))
      tet_split(4, 3) = vertex_id(table(7))

      ! tet 4
      tet_split(1, 4) = vertex_id(table(1))
      tet_split(2, 4) = vertex_id(table(2))
      tet_split(3, 4) = vertex_id(table(3))
      tet_split(4, 4) = vertex_id(table(6))

      ! tet 5
      tet_split(1, 5) = vertex_id(table(1))
      tet_split(2, 5) = vertex_id(table(4))
      tet_split(3, 5) = vertex_id(table(7))
      tet_split(4, 5) = vertex_id(table(3))

      ! tet 6
      tet_split(1, 6) = vertex_id(table(1))
      tet_split(2, 6) = vertex_id(table(7))
      tet_split(3, 6) = vertex_id(table(6))
      tet_split(4, 6) = vertex_id(table(3))

    ELSE IF (n_diags .EQ. 3) THEN
      ! tet 1
      tet_split(1, 1) = vertex_id(table(1))
      tet_split(2, 1) = vertex_id(table(3))
      tet_split(3, 1) = vertex_id(table(4))
      tet_split(4, 1) = vertex_id(table(7))

      ! tet 2
      tet_split(1, 2) = vertex_id(table(1))
      tet_split(2, 2) = vertex_id(table(4))
      tet_split(3, 2) = vertex_id(table(8))
      tet_split(4, 2) = vertex_id(table(7))

      ! tet 3
      tet_split(1, 3) = vertex_id(table(1))
      tet_split(2, 3) = vertex_id(table(8))
      tet_split(3, 3) = vertex_id(table(5))
      tet_split(4, 3) = vertex_id(table(7))

      ! tet 4
      tet_split(1, 4) = vertex_id(table(1))
      tet_split(2, 4) = vertex_id(table(6))
      tet_split(3, 4) = vertex_id(table(7))
      tet_split(4, 4) = vertex_id(table(5))

      ! tet 5
      tet_split(1, 5) = vertex_id(table(2))
      tet_split(2, 5) = vertex_id(table(6))
      tet_split(3, 5) = vertex_id(table(7))
      tet_split(4, 5) = vertex_id(table(1))

      ! tet 6
      tet_split(1, 6) = vertex_id(table(2))
      tet_split(2, 6) = vertex_id(table(7))
      tet_split(3, 6) = vertex_id(table(3))
      tet_split(4, 6) = vertex_id(table(1))
    END IF

  END SUBROUTINE splitHexIntoTet

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  !> Splits a triangular prism into three tets using the rule in
  !! Algorithm described in: HOW TO SUBDIVIDE PYRAMIDS, PRISMS
  !! AND HEXAHEDRA INTO TETRAHEDRA. Dompierre, J. et al.
  !! A quadrilateral face is subdivided into two triangular
  !! faces by the diagonal issuing from the smallest vertex
  !! of the face
  !!
  !! @param vertex_id a vector of length 6 containing the global ids of the vertices
  !! @param tri_split 4x3 array containing the 4 global ids of the 3 tets
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  SUBROUTINE splitPrismToTet(vertex_id, tet_split)
    INTEGER, INTENT(in) :: vertex_id(6)
    INTEGER, INTENT(inout) :: tet_split(4, 3)

    ! Define local variables
    INTEGER :: smallest_pos, smallest_value, i
    INTEGER :: table(6)

    ! Find smallest vertex id
    smallest_pos = 1
    smallest_value = vertex_id(1)
    DO i = 2, 6
      IF (smallest_value .GT. vertex_id(i)) THEN
        smallest_value = vertex_id(i)
        smallest_pos = i
      END IF
    END DO

    IF (smallest_pos .EQ. 1) THEN
      table = (/1,2,3,4,5,6/)
    ELSE IF (smallest_pos .EQ. 2) THEN
      table = (/2,3,1,5,6,4/)
    ELSE IF (smallest_pos .EQ. 3) THEN
      table = (/3,1,2,6,4,5/)
    ELSE IF (smallest_pos .EQ. 4) THEN
      table = (/4,6,5,1,3,2/)
    ELSE IF (smallest_pos .EQ. 5) THEN
      table = (/5,4,6,2,1,3/)
    ELSE IF (smallest_pos .EQ. 6) THEN
      table = (/6,5,4,3,2,1/)
    END IF

    IF (min(vertex_id(table(2)), vertex_id(table(6))) < min(vertex_id(table(3)), vertex_id(table(5)))) THEN
      ! tet 1
      tet_split(1, 1) = vertex_id(table(1))
      tet_split(2, 1) = vertex_id(table(2))
      tet_split(3, 1) = vertex_id(table(3))
      tet_split(4, 1) = vertex_id(table(6))

      ! tet 2
      tet_split(1, 2) = vertex_id(table(1))
      tet_split(2, 2) = vertex_id(table(2))
      tet_split(3, 2) = vertex_id(table(6))
      tet_split(4, 2) = vertex_id(table(5))

      ! tet 3
      tet_split(1, 3) = vertex_id(table(1))
      tet_split(2, 3) = vertex_id(table(5))
      tet_split(3, 3) = vertex_id(table(6))
      tet_split(4, 3) = vertex_id(table(4))
    ELSE
      ! tet 1
      tet_split(1, 1) = vertex_id(table(1))
      tet_split(2, 1) = vertex_id(table(2))
      tet_split(3, 1) = vertex_id(table(3))
      tet_split(4, 1) = vertex_id(table(5))

      ! tet 2
      tet_split(1, 2) = vertex_id(table(1))
      tet_split(2, 2) = vertex_id(table(5))
      tet_split(3, 2) = vertex_id(table(3))
      tet_split(4, 2) = vertex_id(table(6))

      ! tet 3
      tet_split(1, 3) = vertex_id(table(1))
      tet_split(2, 3) = vertex_id(table(5))
      tet_split(3, 3) = vertex_id(table(6))
      tet_split(4, 3) = vertex_id(table(4))
    END IF

  END SUBROUTINE splitPrismToTet

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  !> Splits a quadrilaterial into two triangles using the rule in
  !! Algorithm described in: HOW TO SUBDIVIDE PYRAMIDS, PRISMS
  !! AND HEXAHEDRA INTO TETRAHEDRA. Dompierre, J. et al.
  !! A quadrilateral face is subdivided into two triangular
  !! faces by the diagonal issuing from the smallest vertex
  !! of the face
  !!
  !! @param vertex_id a vector of length 4 containing the global ids of the vertices
  !! @param tri_split 3x2 array containing the 3 global ids of the 2 triangles
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  SUBROUTINE splitQuadIntoTri(vertex_id, tri_split)
    INTEGER, INTENT(in) :: vertex_id(4)
    INTEGER, INTENT(inout) :: tri_split(3, 2)

    ! Define local variables
    INTEGER :: smallest_pos, smallest_value, i

    ! Find smallest vertex id
    smallest_pos = 1
    smallest_value = vertex_id(1)
    DO i = 2, 4
      IF (smallest_value .GT. vertex_id(i)) THEN
        smallest_value = vertex_id(i)
        smallest_pos = i
      END IF
    END DO

    IF (smallest_pos .EQ. 1 .OR. smallest_pos .EQ. 3) THEN
      tri_split(1, 1) = vertex_id(1)
      tri_split(2, 1) = vertex_id(3)
      tri_split(3, 1) = vertex_id(4)

      tri_split(1, 2) = vertex_id(1)
      tri_split(2, 2) = vertex_id(2)
      tri_split(3, 2) = vertex_id(3)
    ELSE IF (smallest_pos .EQ. 2 .OR. smallest_pos .EQ. 4) THEN
      tri_split(1, 1) = vertex_id(1)
      tri_split(2, 1) = vertex_id(2)
      tri_split(3, 1) = vertex_id(4)

      tri_split(1, 2) = vertex_id(2)
      tri_split(2, 2) = vertex_id(3)
      tri_split(3, 2) = vertex_id(4)
    END IF

  END SUBROUTINE splitQuadIntoTri

END MODULE globals
