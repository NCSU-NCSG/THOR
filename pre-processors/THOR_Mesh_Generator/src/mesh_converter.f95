!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! THOR MESH GENERATION UTILITY
!   Generate Mesh Module:
!
!> @mainpage THOR Mesh Generation Pre-Processor
!!    This pre-processor provides mesh format conversion capabilities. It is
!!    designed to allow for a simple workflow from Exodus II to Thor Mesh
!!    using Gmesh 3.0 as an intermediary format.
!! @todo document how to run this code
!
!>   This driver provides the control structure for all entities within the
!!   THOR_MESH_GENERATOR utility.It parses command line input, opens the input
!!   file, and calls the appropriate conversion functions based on file
!!   extension.
!
!> @author Raffi Yessayan
!> @author Nicholas Herring
!> @author Sebastian Schunert
!> @version 1.0
!> @date Mar, 2018
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
PROGRAM generateMesh

  USE globals
  USE gmesh
  USE unv
  USE thor_mesh
  USE integer_array_tools

  !> The number of provided command line arguments
  INTEGER :: arg_count = 0
  !> An integer for status feedback
  INTEGER :: stat
  !> Loop variable
  INTEGER :: j
  !> String for temporarily holding input data for validation prior to storing
  !! in globals
  CHARACTER(200) :: temp_string = ""
  !> Holds the extension of the input mesh name specified in the input file
  CHARACTER(200) :: in_extension = ""
  !> Holds the extension of the output mesh name specified in the input file
  CHARACTER(200) :: out_extension = ""
  !> Holds the libmesh build directory env var
  CHARACTER(500) :: libmesh_dir
  !> Holds command line commands
  CHARACTER(500) :: command

  ! Print some information
  WRITE(6, *) "TODO: Add a general header printer in lib common to all modules"
  CALL printProgramHeader()

  !Get number of command line args and check for proper invocation
  arg_count = COMMAND_ARGUMENT_COUNT()

  IF (arg_count .LT. 2 .OR. MOD(arg_count, 2) .NE. 0) &
        CALL generateErrorMessage(err_code_default, err_fatal, 'Must have an even, nonzero number of command line args.')

  DO j = 1, arg_count, 2
    !Loop through all odd command line args
    CALL GET_COMMAND_ARGUMENT(j, temp_string)

    IF (TRIM(ADJUSTL(temp_string)) .EQ. '-i' ) THEN
      CALL GET_COMMAND_ARGUMENT(j + 1, std_in_file)
    ELSE IF (TRIM(ADJUSTL(temp_string)) .EQ. '-r' ) THEN
      CALL GET_COMMAND_ARGUMENT(j + 1, temp_string)
      READ(temp_string, *, iostat = stat)  uniform_refine
      IF (stat .NE. 0) &
          CALL generateErrorMessage(err_code_default, err_fatal, 'Argument following -r is not an integer.')
    ELSE
      CALL generateErrorMessage(err_code_default, err_fatal, 'Unknown argument. Invocation must be: ./prog -i input.in -r <n>')
    END IF

  END DO

  !At least the input file must have been set
  IF (TRIM(std_in_file) .EQ. "") &
      CALL generateErrorMessage(err_code_default, err_fatal, 'Standard input &
      file must be provided following -i command line argument.')

  !Get the input parameters first
  CALL ingestInput()

  !Get input file extension
  DO i = 200, 1, -1
    in_extension(i:i) = in_file(i:i)
    IF (in_file(i:i) .EQ. '.' ) EXIT
    IF(i .EQ. 1) &
          CALL generateErrorMessage(err_code_default, err_fatal, 'Input file must have an extension to be processed')
  END DO
  in_extension = ADJUSTL(in_extension)

  !Get output file extension
  DO i = 200, 1, -1
    out_extension(i:i) = out_file(i:i)
    IF (out_file(i:i) .EQ. '.' ) EXIT
    IF(i .EQ. 1) &
          CALL generateErrorMessage(err_code_default, err_fatal, 'Output file must have an extension to be processed')
  END DO
  out_extension = ADJUSTL(out_extension)

  !Ingest File
  original_in_file = TRIM(in_file)
  SELECT CASE (in_extension)
  CASE ('.msh') !Gmesh 3.0

    IF (uniform_refine .GT. 0) THEN
      CALL get_environment_variable("THOR_LIBMESH_DIRECTORY", libmesh_dir)
      WRITE(temp_string, *) "-r ", uniform_refine
      command = TRIM(libmesh_dir) // "/meshtool-opt -i " // TRIM(in_file) &
                // " -o internally_generated_gmesh_file.msh " // TRIM(temp_string)
      in_file = "internally_generated_gmesh_file.msh"
      CALL execute_command_line(TRIM(command), exitstat = stat)
    END IF

    CALL ingestGmesh()

    !Clean up the internally generated file if it exists
    IF (uniform_refine .GT. 0) THEN
      command = "rm internally_generated_gmesh_file.msh"
      CALL execute_command_line(TRIM(command), exitstat = stat)
    END IF

    execution_mode = 1
  CASE ('.e') !Exodus II
    !Call libtool & then call ingestGmesh on the new file
    CALL get_environment_variable("THOR_LIBMESH_DIRECTORY", libmesh_dir)

    temp_string = TRIM(libmesh_dir) // "/meshtool-opt -i " // TRIM(in_file) // " -o internally_generated_gmesh_file.msh"
    IF (uniform_refine .GT. 0) THEN
      WRITE(command, *) temp_string, " ", "-r ", uniform_refine
    ELSE
      command = TRIM(temp_string)
    END IF
    in_file = "internally_generated_gmesh_file.msh"
    CALL execute_command_line(TRIM(command), exitstat = stat)

    CALL ingestGmesh()

    !Clean up the internally generated file
    command = "rm internally_generated_gmesh_file.msh"
    CALL execute_command_line(TRIM(command), exitstat = stat)

    execution_mode = 2
  CASE('.UNV','.unv')
    CALL ingestUNV()
    execution_mode = 3
  CASE('.LTHRM','.lthrm')
    CALL ingestThorMesh()
    execution_mode = 4
  CASE DEFAULT
    CALL generateErrorMessage(err_code_default, err_fatal, "Input file extension not supported")
  END SELECT

  ! Read in remapping for region, source, and boundary ids
  CALL setupOptionalReMapping()

  ! Print information about the input to screen
  CALL echoIngestedInput()

  IF (execution_mode .eq. 4_li) THEN
    CALL oldToNewTHORFormat()
  ELSE
    ! Actual reformatting work
    CALL computeAdjacencyList()
    CALL getBoundaryElements()
  END IF

  !Output File
  SELECT CASE (out_extension)
  CASE ('.thrm') !THOR_Mesh
    CALL outputThorMesh()
  CASE DEFAULT
    CALL generateErrorMessage(err_code_default, err_fatal, "Input file extension not supported")
  END SELECT

  ! Safe deallocation
  IF (ALLOCATED(block_id_map)) DEALLOCATE(block_id_map)
  IF (ALLOCATED(source_id_map)) DEALLOCATE(source_id_map)
  IF (ALLOCATED(adjacency_map)) DEALLOCATE(adjacency_map)
  IF (ALLOCATED(boundary_element_list)) DEALLOCATE(boundary_element_list)
  IF (ALLOCATED(boundary_face_list)) DEALLOCATE(boundary_face_list)
  IF (ALLOCATED(block_id)) DEALLOCATE(block_id)
  IF (ALLOCATED(source_id)) DEALLOCATE(source_id)
  IF (ALLOCATED(side_set)) DEALLOCATE(side_set)

  ! TODO: replace with a termination subroutine residing in the lib folder, printing
  ! exec time and a general message that all went well
  WRITE(6, *) "TODO We need a successful termination message here. Maybe have that in lib too?"
END PROGRAM generateMesh
