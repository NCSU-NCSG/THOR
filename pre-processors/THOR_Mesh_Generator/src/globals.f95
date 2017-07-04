!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! THOR MESH GENERATION UTILITY
!   globals module:
!     Conatins variables common to thor_mesh and input files
!     This is intended to enable an *.e (exodus II) to *.thm toolchain for
!     probelm development purposes
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
MODULE globals
  USE integer_array_tools
  IMPLICIT NONE

  INTEGER:: std_in_unit = 19, in_unit = 20, out_unit = 21, region_id_unit = 22,&
            source_id_unit = 23, io_status = 0

  LOGICAL:: skip_region_map, skip_source_map

  INTEGER:: err_code_default = -1, err_warning = 0, err_fatal = 1

  CHARACTER(200):: std_in_file = "", in_file = "", out_file = "",&
                   region_id_file = "", source_id_file = ""

  INTEGER:: node_count, element_count, bc_count

  INTEGER, ALLOCATABLE:: element_list(:,:), block_id(:)

  INTEGER, ALLOCATABLE:: bc_list(:,:), side_set(:)

  REAL, ALLOCATABLE:: node_list(:,:)

  ! a list of three-tuples, each entry is (old_block_id, new_block_id)
  ! defaults to (a, a) where a is from the block_id array
  INTEGER:: nblocks
  INTEGER, ALLOCATABLE, DIMENSION(:, :) :: block_id_map, source_id_map, adjacency_map
  INTEGER, ALLOCATABLE :: boundary_element_list(:), boundary_face_list(:)

  INTEGER, PARAMETER :: li = selected_int_kind(8)
  INTEGER, PARAMETER :: d_t = selected_real_kind(15,307)

CONTAINS

  SUBROUTINE checkStatus(status, expected, message)
    INTEGER, INTENT(IN):: status, expected
    CHARACTER(*), INTENT(IN):: message
    IF (status .NE. expected) &
      CALL generateErrorMessage(err_code_default, err_fatal, message)
  END SUBROUTINE checkStatus

  SUBROUTINE checkFileExists(filename)
    CHARACTER(len=*), intent(in):: filename
    LOGICAL:: file_exists
    INQUIRE(FILE = TRIM(filename), EXIST = file_exists)
    IF (.NOT. file_exists) CALL generateErrorMessage(err_code_default, err_fatal, 'File does not exist')
  END SUBROUTINE checkFileExists

  SUBROUTINE generateErrorMessage(error_code, fatal, message)
    INTEGER, INTENT(IN):: error_code, fatal
    CHARACTER(*), INTENT(IN):: message
    WRITE(*,'(A)') '============================================================'
    WRITE(*,'(A)')'AN UNEXPECTED ERROR OCCURRED'
    WRITE(*,'(A,I0)') 'Error Code: ', error_code
    WRITE(*,'(A)')'Error Message: ', TRIM(message)
    WRITE(*,'(A)') '============================================================'
    IF (fatal .EQ. 1) STOP
  END SUBROUTINE generateErrorMessage

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
      write(*,*) 'Reg instructions ', region_instructions(j, :)
    END DO
    DO j = 1, nblocks
      write(*,*) 'Reg map ', block_id_map(j, :)
    END DO

    DO j = 1, nentries_source
      write(*,*) 'Src instructions ', source_instructions(j, :)
    END DO
    DO j = 1, nblocks
      write(*,*) 'Src map ', source_id_map(j, :)
    END DO

    DEALLOCATE(unique_ids, order)
    IF (ALLOCATED(region_instructions)) DEALLOCATE(region_instructions)
    IF (ALLOCATED(source_instructions)) DEALLOCATE(source_instructions)
  END SUBROUTINE setupRegionAndSourceMapping

  SUBROUTINE computeAdjacencyList()
    INTEGER :: j
    ALLOCATE(adjacency_map(element_count, 4))
    CALL tet_mesh_neighbor_tets(4, element_count, TRANSPOSE(element_list),&
                                TRANSPOSE(adjacency_map))
  END SUBROUTINE computeAdjacencyList

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

END MODULE globals
