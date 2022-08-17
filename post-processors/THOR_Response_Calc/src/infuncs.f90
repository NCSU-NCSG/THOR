!input functions
MODULE infuncs
  USE globals
  USE string_module
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: readinp
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Gets command line arguments and read input file
  SUBROUTINE readinp()
    CHARACTER(64) :: response_inp
    INTEGER(ki4) :: arg_count

    arg_count = COMMAND_ARGUMENT_COUNT()

    IF(arg_count .NE. 1) STOP 'only give input file argument'

    CALL GET_COMMAND_ARGUMENT(1,response_inp)

    CALL read_input_file(response_inp)
  ENDSUBROUTINE readinp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !read in the input file
  SUBROUTINE read_input_file(response_inp)
    CHARACTER(*),INTENT(IN) :: response_inp
    CHARACTER(64) :: t_char,t_char2
    INTEGER(ki4) :: in_unit=22,t_int

    OPEN(UNIT=in_unit, FILE=response_inp, STATUS='OLD', ACTION = "READ", IOSTAT=t_int, IOMSG=t_char)
    IF(t_int .NE. 0)THEN
      WRITE(*,'(A)')t_char
      STOP 'file error'
    ENDIF

    !get the flux data
    REWIND(in_unit)
    DO
      READ(in_unit,*,IOSTAT=t_int)t_char,t_char2
      IF(t_int .NE. 0)STOP 'not found'
      IF(t_char .EQ. 'flux_file')THEN
        CALL read_in_flux(t_char2)
        EXIT
      ENDIF
    ENDDO

    !get the response type
    REWIND(in_unit)
    DO
      READ(in_unit,*,IOSTAT=t_int)t_char,response_type
      IF(t_int .NE. 0)STOP 'not found'
      IF(t_char .EQ. 'response_type')THEN
        IF(response_type .NE. 'cell_wise')THEN
          STOP 'only cell_wise supported for now'
        ENDIF
        EXIT
      ENDIF
    ENDDO

    !get the response function
    REWIND(in_unit)
    DO
      READ(in_unit,*,IOSTAT=t_int)t_char,t_char2
      IF(t_int .NE. 0)STOP 'not found'
      IF(t_char .EQ. 'response_func')THEN
        CALL read_in_map(t_char2)
        EXIT
      ENDIF
    ENDDO

    CLOSE(in_unit)
  ENDSUBROUTINE read_input_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !read in the flux data
  SUBROUTINE read_in_flux(filename)
    CHARACTER(*),INTENT(IN) :: filename
    CHARACTER(100000) :: t_char
    INTEGER(ki4) :: in_unit=23,t_int,nwords,i
    CHARACTER(40) :: words(10000)

    OPEN(UNIT=in_unit, FILE=filename, STATUS='OLD', ACTION = "READ", IOSTAT=t_int, IOMSG=t_char)
    IF(t_int .NE. 0)THEN
      WRITE(*,'(A)')t_char
      STOP 'file error'
    ENDIF

    !get the number of cells
    READ(in_unit,*)num_cells

    !get the number of energy groups
    READ(in_unit,'(A)')t_char
    t_char=TRIM(ADJUSTL(t_char))
    CALL parse(t_char,' ',words,nwords)
    num_groups=nwords-1

    ALLOCATE(flux(num_cells,num_groups),volume(num_cells),resp_func(num_cells,num_groups))

    REWIND(in_unit)

    READ(in_unit,*)nwords
    DO i=1,num_cells
      READ(in_unit,*)volume(i),flux(i,:)
    ENDDO

    CLOSE(in_unit)
  ENDSUBROUTINE read_in_flux

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !read in the cell based response map function
  SUBROUTINE read_in_map(filename)
    CHARACTER(*),INTENT(IN) :: filename
    INTEGER(ki4) :: num_resp_cells,in_unit=23,t_int,i
    CHARACTER(64) :: t_char

    OPEN(UNIT=in_unit, FILE=filename, STATUS='OLD', ACTION = "READ", IOSTAT=t_int, IOMSG=t_char)
    IF(t_int .NE. 0)THEN
      WRITE(*,'(A)')t_char
      STOP 'file error'
    ENDIF

    !get the number of cells
    READ(in_unit,*)t_char,num_resp_cells

    !read in the response function
    resp_func=0.0D0
    DO i=1,num_resp_cells
      READ(in_unit,*)t_int
      BACKSPACE(in_unit)
      READ(in_unit,*)t_char,resp_func(t_int,:)
    ENDDO

    CLOSE(in_unit)
  ENDSUBROUTINE read_in_map
END MODULE infuncs























