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
!> @author Sebastian Schunert
!> @version 1.0
!> @date July, 2017
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
PROGRAM generateMesh

  USE globals
  USE gmesh
  USE thor_mesh
  USE integer_array_tools

  !> The number of provided command line arguments
  INTEGER :: arg_count = 0
  !> String for temporarily holding input data for validation prior to storing
  !! in globals
  CHARACTER(200) :: temp_string = ""
  !> Holds the extension of the input mesh name specified in the input file
  CHARACTER(200) :: in_extension = ""
  !> Holds the extension of the output mesh name specified in the input file
  CHARACTER(200) :: out_extension = ""

  ! Print some information
  WRITE(6, *) "TODO: Add a general header printer in lib common to all modules"
  CALL printProgramHeader()

  !Get number of command line args and check for proper invocation
  arg_count = COMMAND_ARGUMENT_COUNT()

  IF (arg_count .NE. 2) &
        CALL generateErrorMessage(err_code_default, err_fatal, 'Invocation must be: ./prog -i input.in')

  !Get the input command line arg & ensure it specifies input file
  CALL GET_COMMAND_ARGUMENT(1, temp_string)

  IF(TRIM(ADJUSTL(temp_string)) .NE. '-i' ) THEN
    CALL generateErrorMessage(err_code_default, err_fatal, 'Invocation must be: ./prog -i input.in')
  ELSE
    CALL GET_COMMAND_ARGUMENT(2, std_in_file)
  END IF

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
  SELECT CASE (in_extension)
  CASE ('.msh') !Gmesh 3.0
    CALL ingestGmesh()
  CASE ('.e') !Exodus II
    !Call libtool & then call ingestGmesh on the new file
    CALL generateErrorMessage(err_code_default, err_fatal, "Exodus II file support is not yet complete")
  CASE DEFAULT
    CALL generateErrorMessage(err_code_default, err_fatal, "Input file extension not supported")
  END SELECT

  ! Read in remapping for region, source, and boundary ids
  CALL setupOptionalReMapping()

  ! Print information about the input to screen
  CALL echoIngestedInput()

  ! Actual reformatting work
  CALL computeAdjacencyList()
  CALL getBoundaryElements()

  !Output File
  SELECT CASE (out_extension)
  CASE ('.thrm') !THOR_Mesh
    CALL outputThorMesh()
  CASE DEFAULT
    CALL generateErrorMessage(err_code_default, err_fatal, "Input file extension not supported")
  END SELECT

  ! Raffi: Shoudn't we deallocate?
  DEALLOCATE(block_id_map, source_id_map, adjacency_map, boundary_element_list,&
        boundary_face_list)

  ! TODO: replace with a termination subroutine residing in the lib folder, printing
  ! exec time and a general message that all went well
  WRITE(6, *) "TODO We need a successful termination message here. Maybe have that in lib too?"
END PROGRAM generateMesh
