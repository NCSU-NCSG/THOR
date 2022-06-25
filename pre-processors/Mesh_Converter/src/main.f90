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

  !set log filename
  mesh_outfile=TRIM(ADJUSTL(mesh_infile))
  !find extension start
  DO i=LEN_TRIM(mesh_outfile),1,-1
    IF(mesh_outfile(i:i) .EQ. '.')EXIT
  ENDDO
  !if it has an extension, check if it's an input extension and cut it from the logname
  temp_string=TRIM(mesh_outfile)
  IF(i .GE. 2)THEN
    temp_string=mesh_outfile(i:LEN_TRIM(mesh_outfile))
    SELECTCASE(TRIM(temp_string))
      CASE('.msh')
        temp_string=mesh_outfile(1:i-1)
      CASE DEFAULT
        temp_string=TRIM(mesh_outfile)
    ENDSELECT
  ENDIF
  mesh_outfile=TRIM(temp_string)//'.thrm'

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

  WRITE(*,'(A)')'----------------------- Reading in gmsh:'
  CALL read_gmsh_file()


  WRITE(*,'(A)')'----------------------- Calculating Adjacencies:'
  CALL adjacency_calc()

  WRITE(*,'(A)')'----------------------- Outputting thrm file:'
  CALL output_thrm_file()

  WRITE(*,'(A)')'----------------------- Calculating volumes:'
  CALL calcvols()

  WRITE(*,'(A)')'--------------------------------------------------------------------------------'
  WRITE(*,'(A)')'--------------------------------------------------------------------------------'
  WRITE(*,'(A)')'--------------------------------------------------------------------------------'
  WRITE(*,'(A)')'------------------------ THOR mesh converter successful ------------------------'
  WRITE(*,'(2A)')'----------------------- Output written to ',TRIM(ADJUSTL(mesh_outfile))
END PROGRAM thor_mesh_converter
