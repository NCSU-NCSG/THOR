MODULE read_source_module
  !***********************************************************************
  !
  ! Read source module contains all subroutines needed to read source file
  !
  !***********************************************************************
  USE mpi
  USE types
  USE parameter_types
  USE filename_types
  USE multindex_types
  USE global_variables
  USE termination_module

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: read_src

!local opening unit
  INTEGER :: local_unit

CONTAINS

  SUBROUTINE read_src
    !*********************************************************************
    !
    ! Subroutine reads cross-section in 'unique' ahot format
    ! (could be adapted for other formats, of course)
    !
    !*********************************************************************

    INTEGER(kind=li) :: rank,ios,mpi_err,i
    CHARACTER(200) :: src_format

    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)
    local_unit=200+rank

    ! Open and read source file
    OPEN(unit=local_unit,file=TRIM(source_filename),status='old',action='read',IOSTAT=ios)
    IF(ios .NE. 0)THEN
      WRITE(*,*)'error opening ',TRIM(source_filename)
      STOP 'fatal error'
    ENDIF

    DO
      READ(local_unit,*,IOSTAT=ios)src_format
      IF(ios .NE. 0)THEN
        !legacy original THOR source format with no format indicater
        REWIND(local_unit)
        CALL src_read_legacyv0()
        EXIT
      ENDIF
      src_format=TRIM(ADJUSTL(src_format))
      IF(src_format .EQ. 'THOR_SRC_V1')THEN
        !current THOR source format
        REWIND(local_unit)
        CALL src_read_current()
        EXIT
      ENDIF
    ENDDO

    !close src file
    CLOSE(local_unit)

    src_id_min=9999999
    src_id_max=-9999999
    !find max and min id values
    DO i=1,num_src_mat
      IF(ext_src(i)%src_id .GE. src_id_max)src_id_max=ext_src(i)%src_id
      IF(ext_src(i)%src_id .LE. src_id_min)src_id_min=ext_src(i)%src_id
    ENDDO
    !allocate and assign pointer values
    !everywhere else will be 0
    ALLOCATE(src_pointer(src_id_min:src_id_max))
    src_pointer=0
    DO i=1,num_src_mat
      src_pointer(ext_src(i)%src_id)=i
    ENDDO

  END SUBROUTINE read_src

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE src_read_legacyv0
    !*********************************************************************
    !
    ! Subroutine reads source in 'unique' ahot format
    ! (could be adapted for other formats, of course)
    !
    !*********************************************************************

    ! Declare temporary variables

    INTEGER(kind=li) :: alloc_stat, eg, m, l

    REAL(kind=d_t) :: src_str

    ! Read source strength and moments from file

    READ(local_unit,*) num_src_mat

    !allocate the external source. For legacy input, only spatial moments are allowed
    ALLOCATE(ext_src(num_src_mat),stat=alloc_stat)
    IF(alloc_stat /= 0) CALL stop_thor(2_li)
    ext_src(:)%src_id=0.0
    DO m=1,num_src_mat
      ALLOCATE(ext_src(m)%mom(num_moments_v,namom,egmax),stat=alloc_stat)
      IF(alloc_stat /= 0) CALL stop_thor(2_li)
      ext_src(m)%mom(:,:,:)=0.0
    ENDDO

    DO m=1, num_src_mat
      READ(local_unit,*) ext_src(m)%src_id
      IF (ext_src(m)%src_id .ge. num_src_mat) CALL stop_thor(-1_li,&
                                       "source region ID is 0-indexed and must be < # source regions")
      DO eg=1, egmax
        src_str=0.0
        READ(local_unit,*) src_str
        READ(local_unit,*) (ext_src(m)%mom(l,1,eg), l = 1, num_moments_v)
        DO l=1,num_moments_v
          ext_src(m)%mom(l,1,eg)=src_str*ext_src(m)%mom(l,1,eg)
        ENDDO
      END DO
    END DO

  END SUBROUTINE src_read_legacyv0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE src_read_current

  STOP 'src_read_current not finished'

  END SUBROUTINE src_read_current

END MODULE read_source_module
