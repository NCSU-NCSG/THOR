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
  !> Used to track I/O operation status codes
  INTEGER :: io_status = 0
  !> Default error code when generating an error message without a more specific
  !! error type
  INTEGER :: err_code_default = -1
  !> Error severity level requiring output to user, but no other action
  INTEGER :: err_warning = 0
  !> Error severity level requiring immediate termination of program
  INTEGER :: err_fatal = 1
  !> Number of <x,y,z> nodes in the mesh
  INTEGER :: node_count
  !> Number of elements in the mesh
  INTEGER :: element_count
  !> Number of boundary conditions faces defined in the mesh
  INTEGER :: bc_count
  !> Number of block ids defined in the mesh
  INTEGER:: nblocks
  !> Paramter defined int, width 8
  INTEGER, PARAMETER :: li = selected_int_KIND(8)
  !> Parameter defined double precision real
  INTEGER, PARAMETER :: d_t = selected_real_KIND(15,307)
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
  !> Holds the name of the standard input config file
  CHARACTER(200) :: std_in_file = ""
  !> Holds the name of the standard input mesh file
  CHARACTER(200) :: in_file = ""
  !> Holds the name of the standard output mesh file
  CHARACTER(200) :: out_file = ""
  !> Holds the name of the standard input region id file
  CHARACTER(200) :: region_id_file = ""
  !> Holds the name of the standard input source id file
  CHARACTER(200) :: source_id_file = ""

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
    ! region id mapping can be skipped
    IF (io_status .NE. 0) THEN
      skip_region_map = .TRUE.
    ELSE
      skip_region_map = .FALSE.
      region_id_file = TRIM(line)
    END IF

    ! src id mapping file
    READ(std_in_unit, '(A200)', IOSTAT = io_status) line

    ! region id mapping can be skipped
    IF (io_status .NE. 0) THEN
      skip_source_map = .TRUE.
    ELSE
      skip_source_map = .FALSE.
      source_id_file = TRIM(line)
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
  SUBROUTINE setupRegionAndSourceMapping()

    INTEGER :: j, position
    INTEGER, DIMENSION(:), ALLOCATABLE :: unique_ids
    INTEGER, DIMENSION(:), ALLOCATABLE :: order
    INTEGER :: nentries_region, nentries_source
    INTEGER, DIMENSION(:, :), ALLOCATABLE :: region_instructions
    INTEGER, DIMENSION(:, :), ALLOCATABLE :: source_instructions

    ! Set the default mapping
    nblocks = numUniqueEntries(block_id)
    ALLOCATE(block_id_map(nblocks, 2), source_id_map(nblocks, 2), unique_ids(nblocks),&
          order(nblocks))
    CALL uniqueEntries(block_id, unique_ids)
    ! sort unique ids to make block_id_map and source_id_map sorted
    CALL quickSortInteger(unique_ids, order)

    ! assign defaults to block_id_map and source_id_map
    DO j = 1, nblocks
      block_id_map(j, :) = unique_ids(j)
      source_id_map(j, :) = unique_ids(j)
    END DO

    ! Read region id and source id from files and override block_id_map default
    IF (.NOT. skip_region_map) THEN
      nentries_region = lengthTextFile(region_id_file, region_id_unit)
      ALLOCATE(region_instructions(nentries_region, 2))
      CALL loadTextFile(region_id_file, region_id_unit, region_instructions)
    END IF

    ! work instructions into block_id_map
    IF (.NOT. hasUniqueEntries(region_instructions(:, 1))) &
          CALL generateErrorMessage(err_code_default, err_fatal, &
          'regions instruction keys are not unique')

    DO j = 1, nentries_region
      position = mapIndexOf(region_instructions(j, 1), block_id_map(:, 1))
      block_id_map(position, 2) = region_instructions(j, 2)
    END DO

    ! Read region id and source id from files and override block_id_map default
    IF (.NOT. skip_source_map) THEN
      nentries_source = lengthTextFile(source_id_file, source_id_unit)
      ALLOCATE(source_instructions(nentries_source, 2))
      CALL loadTextFile(source_id_file, source_id_unit, source_instructions)
    END IF

    ! work instructions into source_id_map
    IF (.NOT. hasUniqueEntries(source_instructions(:, 1))) &
          CALL generateErrorMessage(err_code_default, err_fatal, &
          'source instruction keys are not unique')

    DO j = 1, nentries_source
      position = mapIndexOf(source_instructions(j, 1), source_id_map(:, 1))
      source_id_map(position, 2) = source_instructions(j, 2)
    END DO

    DO j = 1, nentries_region
      WRITE(*,*) 'Reg instructions ', region_instructions(j, :)
    END DO
    DO j = 1, nblocks
      WRITE(*,*) 'Reg map ', block_id_map(j, :)
    END DO

    DO j = 1, nentries_source
      WRITE(*,*) 'Src instructions ', source_instructions(j, :)
    END DO
    DO j = 1, nblocks
      WRITE(*,*) 'Src map ', source_id_map(j, :)
    END DO

    DEALLOCATE(unique_ids, order)
    IF (ALLOCATED(region_instructions)) DEALLOCATE(region_instructions)
    IF (ALLOCATED(source_instructions)) DEALLOCATE(source_instructions)
  END SUBROUTINE setupRegionAndSourceMapping

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
  !> Given an adjacency map has been establsihed, returns a pair of arrays
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
          boundary_element_list(counter) = i
          boundary_face_list(counter) = j - 1
          counter = counter + 1
        END IF
      END DO
    END DO
  END SUBROUTINE getBoundaryElements

  SUBROUTINE printProgramHeader()
    ! TODO: A unified program name header?
    WRITE(6, '(A)') "Program Name: THOR_MESH_GENERATOR"
  END

  SUBROUTINE echoIngestedInput()
    INTEGER :: j
    INTEGER :: s

    WRITE(6, *)
    WRITE(6, '(A)') "------------------------------------ Input ------------------------------------"
    WRITE(6, '(A)') "Infile mesh file name: " // TRIM(in_file)
    WRITE(6, '(A)') "Output mesh file name: " // TRIM(out_file)

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

    WRITE(6, *)
    IF (skip_source_map) THEN
      WRITE(6, '(A)') "No source ID edits are provided"
    ELSE
      IF (.NOT. ALLOCATED(source_id_map)) &
        CALL generateErrorMessage(err_code_default, err_fatal, &
          'source id map not allocated,  echoIngestedInput called too early')
      WRITE(6, '(A,A)') "Reporting Source ID reassignments from file: ", TRIM(source_id_file)
      IF (.NOT. ALLOCATED(source_id_map)) &
        CALL generateErrorMessage(err_code_default, err_fatal, &
          'source id map not allocated,  echoIngestedInput called too early')
      s = SIZE(source_id_map(:, 1))
      DO j = 1, s
        WRITE(6, '(A,I0,A,I0)') "Old Region ID ", source_id_map(j, 1), " Source ID ",&
                                source_id_map(j, 2)
      END DO
    END IF
  END SUBROUTINE echoIngestedInput

END MODULE globals
