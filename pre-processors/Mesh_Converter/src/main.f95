!-------------------------------------------------------------------------------
! THOR MESH converter to go from gmsh to thrm
!> @author Nicholas F. Herring
!-------------------------------------------------------------------------------
PROGRAM thor_mesh_converter
  USE globals
  USE read_gmsh
  USE output_thrm
  USE boundary_conditions
  IMPLICIT NONE

  !> String for temporarily holding input data for validation prior to storing
  !! in globals
  CHARACTER(200) :: temp_string = ""

  !> The number of provided command line arguments
  INTEGER :: arg_count,i,j

  !Get number of command line args and check for proper invocation
  arg_count = COMMAND_ARGUMENT_COUNT()

  IF (arg_count .LT. 1)STOP 'must at least give mesh file'

  !get mesh in file name
  CALL GET_COMMAND_ARGUMENT(1, mesh_infile)
  !get boundary conditions
  i=2
  DO
    CALL GET_COMMAND_ARGUMENT(i, temp_string)
    IF(temp_string .EQ. '-bc')THEN
      DO j=1,6
        i=i+1
        CALL GET_COMMAND_ARGUMENT(i, temp_string)
        READ(temp_string,*)side_bc(j)
      ENDDO
    ENDIF
    IF(i .GE. arg_count)EXIT
    i=i+1
  ENDDO

  CALL read_gmsh_file()

  CALL adjacency_calc()

  CALL output_gmsh_file()
  stop 'testing converter'
END PROGRAM thor_mesh_converter
