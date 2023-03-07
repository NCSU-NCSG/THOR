!THOR is licensed under the MIT License.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!> @brief Input functions.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
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
    INTEGER(ki4) :: arg_count

    arg_count = COMMAND_ARGUMENT_COUNT()

    IF(arg_count .NE. 1) STOP 'only give input file argument'

    CALL GET_COMMAND_ARGUMENT(1,response_inp)

    CALL read_input_file()
  ENDSUBROUTINE readinp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !read in the input file
  SUBROUTINE read_input_file()
    CHARACTER(200) :: t_char,t_char2
    INTEGER(ki4) :: in_unit=22,t_int,num_flux_files,i

    OPEN(UNIT=in_unit, FILE=response_inp, STATUS='OLD', ACTION = "READ", IOSTAT=t_int, IOMSG=t_char)
    IF(t_int .NE. 0)THEN
      WRITE(*,'(A)')TRIM(t_char)
      STOP 'file error'
    ENDIF

    !get the flux data
    REWIND(in_unit)
    DO
      READ(in_unit,*,IOSTAT=t_int)t_char
      IF(t_int .NE. 0)STOP 'flux_files not found'
      IF(t_char .EQ. 'flux_files')THEN
        BACKSPACE(in_unit)
        READ(in_unit,*)t_char,t_char2
        READ(t_char2,*)num_flux_files
        IF(num_flux_files .LT. 1)STOP 'need at least one flux file'
        !get the problem size data, groups and number of cells as well as the volumes
        READ(in_unit,'(A)')t_char2
        CALL get_vols(t_char2)
        BACKSPACE(in_unit)
        !read in all the flux data and sum it up
        DO i=1,num_flux_files
          READ(in_unit,'(A)')t_char2
          CALL read_in_flux(t_char2)
        ENDDO
        EXIT
      ENDIF
    ENDDO

    !get the response type
    REWIND(in_unit)
    DO
      READ(in_unit,*,IOSTAT=t_int)t_char
      IF(t_int .NE. 0)STOP 'response_type not found'
      IF(t_char .EQ. 'response_type')THEN
        BACKSPACE(in_unit)
        READ(in_unit,*)t_char,response_type
        SELECT CASE(response_type)
          CASE('cell_wise','region_wise')
          CASE DEFAULT
            STOP 'only cell_wise, region_wise, or mat_wise supported for now'
        ENDSELECT
        EXIT
      ENDIF
    ENDDO

    !get the response function
    REWIND(in_unit)
    DO
      READ(in_unit,*,IOSTAT=t_int)t_char
      IF(t_int .NE. 0)STOP 'response_func not found'
      IF(t_char .EQ. 'response_func')THEN
        BACKSPACE(in_unit)
        READ(in_unit,'(A)')t_char
        t_char=TRIM(ADJUSTL(t_char))
        t_char2=t_char(14:200)
        t_char2=TRIM(ADJUSTL(t_char2))
        SELECT CASE(response_type)
          CASE('cell_wise')
            CALL read_in_map_cell(t_char2)
          CASE('region_wise')
            CALL read_in_map_reg(t_char2)
        ENDSELECT
        EXIT
      ENDIF
    ENDDO

    CLOSE(in_unit)
  ENDSUBROUTINE read_input_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !read in the flux data
  SUBROUTINE get_vols(filename)
    CHARACTER(*),INTENT(IN) :: filename
    CHARACTER(100000) :: t_char
    INTEGER(ki4) :: in_unit=23,t_int,nwords,i
    CHARACTER(40) :: words(10000)

    OPEN(UNIT=in_unit, FILE=filename, STATUS='OLD', ACTION = "READ", IOSTAT=t_int, IOMSG=t_char)
    IF(t_int .NE. 0)THEN
      WRITE(*,'(A)')TRIM(t_char)
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
    flux=0.0D0

    REWIND(in_unit)

    READ(in_unit,*)nwords
    DO i=1,num_cells
      READ(in_unit,*)volume(i)
    ENDDO

    CLOSE(in_unit)
  ENDSUBROUTINE get_vols

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !read in the flux data
  SUBROUTINE read_in_flux(filename)
    CHARACTER(*),INTENT(IN) :: filename
    CHARACTER(100000) :: t_char
    INTEGER(ki4) :: in_unit=23,t_int,i
    REAL(kr8),ALLOCATABLE :: tline_arr(:)

    OPEN(UNIT=in_unit, FILE=filename, STATUS='OLD', ACTION = "READ", IOSTAT=t_int, IOMSG=t_char)
    IF(t_int .NE. 0)THEN
      WRITE(*,'(A)')TRIM(t_char)
      STOP 'file error'
    ENDIF

    ALLOCATE(tline_arr(num_groups))

    READ(in_unit,*)t_char
    DO i=1,num_cells
      READ(in_unit,*)t_char,tline_arr(:)
      flux(i,:)=flux(i,:)+tline_arr
    ENDDO

    DEALLOCATE(tline_arr)

    CLOSE(in_unit)
  ENDSUBROUTINE read_in_flux

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !read in the cell based response map function
  SUBROUTINE read_in_map_cell(filename)
    CHARACTER(*),INTENT(IN) :: filename
    INTEGER(ki4) :: num_resp_cells,in_unit=23,t_int,i
    CHARACTER(64) :: t_char

    OPEN(UNIT=in_unit, FILE=filename, STATUS='OLD', ACTION = "READ", IOSTAT=t_int, IOMSG=t_char)
    IF(t_int .NE. 0)THEN
      WRITE(*,'(A)')TRIM(t_char)
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
  ENDSUBROUTINE read_in_map_cell

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !read in the region based response map function
  SUBROUTINE read_in_map_reg(filename)
    CHARACTER(*),INTENT(IN) :: filename
    INTEGER(ki4) :: in_unit=23,in_unit2=24,t_int,i,regs(num_cells),num_regions
    CHARACTER(200) :: t_char,mesh_file
    REAL(kr8),ALLOCATABLE :: reg_resp(:,:)

    OPEN(UNIT=in_unit, FILE=filename, STATUS='OLD', ACTION = "READ", IOSTAT=t_int, IOMSG=t_char)
    IF(t_int .NE. 0)THEN
      WRITE(*,'(A)')TRIM(t_char)
      STOP 'file error'
    ENDIF

    READ(in_unit,'(A)')mesh_file

    OPEN(UNIT=in_unit2, FILE=mesh_file, STATUS='OLD', ACTION = "READ", IOSTAT=t_int, IOMSG=t_char)
    IF(t_int .NE. 0)THEN
      WRITE(*,'(A)')TRIM(t_char)
      STOP 'file error'
    ENDIF
    regs=0
    READ(in_unit2,*)num_regions
    READ(in_unit2,*)t_char
    READ(in_unit2,*)t_char
    READ(in_unit2,*)t_char
    !get past the vertices
    DO i=1,num_regions
      READ(in_unit2,*)t_char
    ENDDO
    !read in the regions
    DO i=1,num_cells
      READ(in_unit2,*)t_char,regs(i),t_char
    ENDDO

    CLOSE(in_unit2)

    READ(in_unit,*)t_char,num_regions
    ALLOCATE(reg_resp(num_regions,num_groups))

    !read in the response function based on regions
    reg_resp=0.0D0
    DO i=1,num_regions
      READ(in_unit,*)t_int
      BACKSPACE(in_unit)
      READ(in_unit,*)t_char,reg_resp(t_int,:)
    ENDDO

    CLOSE(in_unit)

    !now actually assign the cell based response function
    resp_func=0
    DO i=1,num_cells
      resp_func(i,:)=reg_resp(regs(i),:)
    ENDDO

    DEALLOCATE(reg_resp)
  ENDSUBROUTINE read_in_map_reg
END MODULE infuncs























