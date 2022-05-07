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
      !boundary conditions are ordered -x, +x, -y, +y, -z, +z
      !0 is vacuum, 1 is reflective, and 2 is fixed inflow
      DO j=1,6
        i=i+1
        CALL GET_COMMAND_ARGUMENT(i, temp_string)
        READ(temp_string,*)side_bc(j)
        IF(side_bc(j) .GT. 2 .OR. side_bc(j) .LT. 0)STOP 'Invalid BC'
      ENDDO
    ENDIF
    IF(i .GE. arg_count)EXIT
    i=i+1
  ENDDO

  WRITE(*,'(A)')'---------------------Reading in gmsh:'
  WRITE(*,*)
  CALL read_gmsh_file()


  WRITE(*,'(A)')'---------------------Calculating Adjacencies:'
  WRITE(*,*)
  CALL adjacency_calc()

  WRITE(*,'(A)')'---------------------Outputting thrm file:'
  WRITE(*,*)
  CALL output_thrm_file()

  WRITE(*,'(A)')'---------------------Calculating volumes:'
  WRITE(*,*)
  CALL calcvols()

  WRITE(*,*)
  WRITE(*,'(A)')'**********************************************************************************'
  WRITE(*,'(A)')'**********************************************************************************'
  WRITE(*,'(A)')'**********************************************************************************'
  WRITE(*,'(A)')'**************************THOR mesh converter sucessful.**************************'
  WRITE(*,'(2A)')'--------------- Output written to ',TRIM(ADJUSTL(mesh_infile))//'_out.thrm'
END PROGRAM thor_mesh_converter
