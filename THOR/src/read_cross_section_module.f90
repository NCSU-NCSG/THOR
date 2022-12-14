MODULE read_cross_section_module
  !***********************************************************************
  !
  ! Read cross-section module  contains all subroutines needed to read
  ! cross-section file
  !
  !***********************************************************************

  USE mpi
  USE types
  USE parameter_types
  USE filename_types
  USE cross_section_types
  USE globals
  USE error_module
  USE stringmod

  IMPLICIT NONE
  PRIVATE
  !
  ! List of public members
  PUBLIC :: read_xs

  INTEGER :: local_unit

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE read_xs
    !*********************************************************************
    !
    ! Subroutine reads cross-section in 'unique' ahot format
    ! (could be adapted for other formats, of course)
    !
    !*********************************************************************

    INTEGER(kind=li) :: rank,ios,mpi_err,i,g,gp,j
    CHARACTER(200) :: xs_format

    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)
    local_unit=100+rank
    ! Open and read mesh file

    OPEN(unit=local_unit,file=TRIM(cross_section_filename),status='old',action='read',IOSTAT=ios)
    IF(ios .NE. 0)THEN
      CALL raise_fatal_error('error opening '//TRIM(cross_section_filename))
    ENDIF

    DO
      READ(local_unit,*,IOSTAT=ios)xs_format
      IF(ios .NE. 0)THEN
        !legacy original THOR xs format with no format indicater
        REWIND(local_unit)
        CALL xs_read_legacy_v0()
        EXIT
      ENDIF
      xs_format=TRIM(ADJUSTL(xs_format))
      IF(xs_format .EQ. 'THOR_XS_V1')THEN
        !current THOR xs format
        REWIND(local_unit)
        CALL xs_read_current()
        EXIT
      ENDIF
    ENDDO

    !close xs file
    CLOSE(local_unit)

    ! Set most_thermal to no upscattering
    most_thermal=egmax+1
    ! Determine most_thermal
    DO i=1,num_mat
      IF(upscattering .NE. 0) THEN
        DO g=1,egmax
          DO gp=g+1,egmax
            IF( ABS(xs_mat(i)%sigma_scat(1,g,gp)) > 2.24E-16_d_t .AND. most_thermal>g ) &
              most_thermal=g
          END DO
        END DO
      END IF
      ! Compute tsigs
      DO gp=1,egmax
        xs_mat(i)%tsigs(gp)=0.0_d_t
        DO g=1,egmax
          xs_mat(i)%tsigs(gp)=xs_mat(i)%tsigs(gp)+xs_mat(i)%sigma_scat(1,g,gp)
        END DO
      END DO
    ENDDO
    ! set neven
    neven=1_li+scatt_ord+(scatt_ord+1_li)*scatt_ord/2_li
    ! set scat_mult
    ALLOCATE(scat_mult(0:scatt_ord,0:scatt_ord))
    scat_mult=0.0_d_t
    IF(scat_mult_flag.EQ.1_li) THEN
      DO j=0,scatt_ord
        DO i=0,j
          IF(i.EQ.0_li) THEN
            scat_mult(j,i)=1.0_d_t/REAL(2_li*j+1_li,d_t)
          ELSE
            scat_mult(j,i)=2.0_d_t/REAL(2_li*j+1_li,d_t)
          END IF
        END DO
      END DO
    ELSE
      DO j=0,scatt_ord
        DO i=0,j
          IF(i.EQ.0_li) THEN
            scat_mult(j,i)=1.0_d_t
          ELSE
            scat_mult(j,i)=2.0_d_t
          END IF
        END DO
      END DO
    END IF

    mat_id_min=9999999
    mat_id_max=-9999999
    !find max and min id values
    DO i=1,num_mat
      IF(xs_mat(i)%mat_id .GE. mat_id_max)mat_id_max=xs_mat(i)%mat_id
      IF(xs_mat(i)%mat_id .LE. mat_id_min)mat_id_min=xs_mat(i)%mat_id
    ENDDO
    !allocate and assign pointer values
    !everywhere else will be 0
    ALLOCATE(material_ids(mat_id_min:mat_id_max))
    material_ids=0
    DO i=1,num_mat
      material_ids(xs_mat(i)%mat_id)=i
    ENDDO
  END SUBROUTINE read_xs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE xs_read_legacy_v0()

    INTEGER(kind=li) :: alloc_stat, e1, order, eg_to, eg_from,m

    READ(local_unit,*) num_mat

    ! Allocate cross-section arrays and check if enough memory is available
    ALLOCATE(xs_mat(num_mat),&
          eg_bounds(egmax+1),&
          stat=alloc_stat)
    IF(alloc_stat /= 0) CALL raise_fatal_error("*** Not enough memory ***")

    DO m=1,num_mat
      ALLOCATE(xs_mat(m)%chi(egmax),xs_mat(m)%sigma_f(egmax),xs_mat(m)%nu(egmax), &
        xs_mat(m)%sigma_t(egmax),xs_mat(m)%tsigs(egmax),xs_mat(m)%nusig_f(egmax), &
        xs_mat(m)%sigma_scat(xs_ord+1,egmax,egmax),stat=alloc_stat)
      IF(alloc_stat /= 0) CALL raise_fatal_error("*** Not enough memory ***")
      !set everything to 0
      xs_mat(m)%chi(:)=0
      xs_mat(m)%sigma_f(:)=0
      xs_mat(m)%nu(:)=0
      xs_mat(m)%sigma_t(:)=0
      xs_mat(m)%tsigs(:)=0
      xs_mat(m)%nusig_f(:)=0
      xs_mat(m)%sigma_scat(:,:,:)=0
    ENDDO
    xs_mat(:)%mat_id=0
    eg_bounds(:)=0
    xs_mat(:)%fissile=.TRUE.

    ! Read material total & scattering cross-section
    DO m=1, num_mat
      READ(local_unit,*) xs_mat(m)%mat_id
      WRITE(xs_mat(m)%mat_name,'(A,I0)')'mat_',xs_mat(m)%mat_id
      IF(multiplying .NE. 0)THEN
        READ(local_unit,*) (xs_mat(m)%chi(e1),e1=1,egmax)
        READ(local_unit,*) (eg_bounds(e1),e1=1,egmax)
        eg_bounds(egmax+1)=0.0_d_t
        READ(local_unit,*) (xs_mat(m)%sigma_f(e1),e1=1,egmax)
        IF(MAXVAL(xs_mat(m)%sigma_f(:)) .LE. 0.0 .AND. MINVAL(xs_mat(m)%sigma_f(:)) .GE. 0.0) &
          xs_mat(m)%fissile=.FALSE.
        READ(local_unit,*) (xs_mat(m)%nu(e1),e1=1,egmax)
        IF(MAXVAL(xs_mat(m)%nu(:)) .LE. 0.0 .AND. MINVAL(xs_mat(m)%nu(:)) .GE. 0.0) &
          xs_mat(m)%fissile=.FALSE.
      ELSE
        xs_mat(m)%fissile=.FALSE.
        DO e1 = 1, egmax
          xs_mat(m)%chi(e1) = zero
          eg_bounds(e1) = zero
          xs_mat(m)%sigma_f(e1) = zero
          xs_mat(m)%nu(e1) = zero
        END DO
        eg_bounds(egmax+1) = zero
      END IF
      READ(local_unit,*) (xs_mat(m)%sigma_t(e1),e1=1,egmax)
      ! Initialize scattering matrix
      DO order=1, xs_ord+1
        DO eg_to=1,egmax
          DO eg_from=1,egmax
            xs_mat(m)%sigma_scat(order,eg_to,eg_from)=zero
          END DO
        END DO
      END DO
      ! Read scattering matrix, note: thermal groups are separated from fast
      ! groups but old cross section format remains valid
      IF(upscattering.EQ.0) THEN
        DO order=1, xs_ord+1
          DO eg_to=1,egmax
            READ(local_unit,*) (xs_mat(m)%sigma_scat(order,eg_to,eg_from),eg_from=1,eg_to)
          END DO
        END DO
      ELSE
        DO order=1, xs_ord+1
          DO eg_to=1,egmax
            READ(local_unit,*) (xs_mat(m)%sigma_scat(order,eg_to,eg_from),eg_from=1,egmax)
          END DO
        END DO
      END IF
      xs_mat(m)%nusig_f(:)=xs_mat(m)%nu(:)*xs_mat(m)%sigma_f(:)
    END DO
  END SUBROUTINE xs_read_legacy_v0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!currently supports max of 1000 groups
  SUBROUTINE xs_read_current()
    CHARACTER(10000), ALLOCATABLE :: words(:)
    INTEGER :: nwords,i,g,alloc_stat,gp,j

    ALLOCATE(words(1000))

    num_mat=0
    egmax=0
    xs_ord=0
    !get the specification data
    DO
      CALL get_next_line(words,nwords)
      words(1)=TRIM(ADJUSTL(words(1)))
      IF(words(1) .EQ. 'THOR_XS_V1')THEN
        EXIT
      ENDIF
    ENDDO

    !get the material, group, and order amount
    words(2)=TRIM(ADJUSTL(words(2)))
    words(3)=TRIM(ADJUSTL(words(3)))
    words(4)=TRIM(ADJUSTL(words(4)))
    READ(words(2),*)num_mat
    READ(words(3),*)egmax
    READ(words(4),*)xs_ord
    ! Allocate cross-section arrays and check if enough memory is available
    ALLOCATE(xs_mat(num_mat),eg_bounds(egmax+1),stat=alloc_stat)
    IF(alloc_stat /= 0) CALL raise_fatal_error("*** Not enough memory ***")

    DO i=1,num_mat
      ALLOCATE(xs_mat(i)%chi(egmax),xs_mat(i)%sigma_f(egmax),xs_mat(i)%nu(egmax), &
        xs_mat(i)%sigma_t(egmax),xs_mat(i)%tsigs(egmax),xs_mat(i)%nusig_f(egmax), &
        xs_mat(i)%sigma_scat(xs_ord+1,egmax,egmax),stat=alloc_stat)
      IF(alloc_stat /= 0) CALL raise_fatal_error("*** Not enough memory ***")
      !set everything to 0
      xs_mat(i)%chi(:)=0
      xs_mat(i)%sigma_f(:)=0
      xs_mat(i)%nu(:)=0
      xs_mat(i)%sigma_t(:)=0
      xs_mat(i)%tsigs(:)=0
      xs_mat(i)%nusig_f(:)=0
      xs_mat(i)%sigma_scat(:,:,:)=0
    ENDDO
    xs_mat(:)%mat_id=0
    xs_mat(:)%mat_name=''
    eg_bounds(:)=0
    xs_mat(:)%fissile=.TRUE.
    !set everything to 0
    xs_mat(:)%mat_id=0
    eg_bounds(:)=0

    !read in energy bounds if present
    CALL get_next_line(words,nwords)
    words(1)=TRIM(ADJUSTL(lowercase(words(1))))
    IF(words(1) .NE. 'id')THEN
      IF(nwords .NE. egmax)CALL raise_fatal_error('bad amount of energy data on line in xs file')
      DO g=1,egmax
        READ(words(g),*)eg_bounds(g)
      ENDDO
      eg_bounds(egmax+1)=0.0
    ENDIF
    REWIND(local_unit)

    !get all material datas
    DO i=1,num_mat
      DO
        !get next line
        CALL get_next_line(words,nwords)
        words(1)=TRIM(ADJUSTL(lowercase(words(1))))
        IF(words(1) .EQ. 'id')THEN
          IF(nwords .LT. 2)CALL raise_fatal_error('no id number in xs file')
          READ(words(2),*)xs_mat(i)%mat_id
          words(4)=TRIM(ADJUSTL(words(4)))
          !get or assign mat name
          IF(words(4) .NE. '')THEN
            xs_mat(i)%mat_name=TRIM(words(4))
          ELSE
            WRITE(xs_mat(i)%mat_name,'(A,I0)')'mat_',xs_mat(i)%mat_id
          ENDIF
          !read in the fission spectrum
          CALL get_next_line(words,nwords)
          IF(nwords .NE. egmax)THEN
            CALL raise_fatal_error('bad amount of xs data on chi line for mat '//TRIM(STR(xs_mat(i)%mat_id)))
          ENDIF
          DO g=1,egmax
            READ(words(g),*)xs_mat(i)%chi(g)
          ENDDO
          !read in SigmaF
          CALL get_next_line(words,nwords)
          IF(nwords .NE. egmax)THEN
            CALL raise_fatal_error('bad amount of xs data on SigmaF line for mat '//TRIM(STR(xs_mat(i)%mat_id)))
          ENDIF
          DO g=1,egmax
            READ(words(g),*)xs_mat(i)%sigma_f(g)
          ENDDO
          IF(MAXVAL(xs_mat(i)%sigma_f(:)) .LE. 0.0 .AND. MINVAL(xs_mat(i)%sigma_f(:)) .GE. 0.0) &
            xs_mat(i)%fissile=.FALSE.
          !read in nu
          CALL get_next_line(words,nwords)
          IF(nwords .NE. egmax)THEN
            CALL raise_fatal_error('bad amount of xs data on nu line for mat '//TRIM(STR(xs_mat(i)%mat_id)))
          ENDIF
          DO g=1,egmax
            READ(words(g),*)xs_mat(i)%nu(g)
          ENDDO
          IF(MAXVAL(xs_mat(i)%nu(:)) .LE. 0.0 .AND. MINVAL(xs_mat(i)%nu(:)) .GE. 0.0) &
            xs_mat(i)%fissile=.FALSE.
          !read in total/transport xs
          CALL get_next_line(words,nwords)
          IF(nwords .NE. egmax)THEN
            CALL raise_fatal_error('bad amount of xs data on total xs line for mat '//TRIM(STR(xs_mat(i)%mat_id)))
          ENDIF
          DO g=1,egmax
            READ(words(g),*)xs_mat(i)%sigma_t(g)
          ENDDO
          !read in scattering xs format gp->g
          DO j=1,xs_ord+1
            DO g=1,egmax
              CALL get_next_line(words,nwords)
              IF(nwords .NE. egmax)THEN
                CALL raise_fatal_error('bad amount of xs data on sigmas line for mat '//TRIM(STR(xs_mat(i)%mat_id)))
              ENDIF
              DO gp=1,egmax
                READ(words(gp),*)xs_mat(i)%sigma_scat(j,g,gp)
              ENDDO
            ENDDO
          ENDDO
          EXIT
        ENDIF
      ENDDO
      xs_mat(i)%nusig_f(:)=xs_mat(i)%nu(:)*xs_mat(i)%sigma_f(:)
    ENDDO
    DEALLOCATE(words)
  END SUBROUTINE xs_read_current

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_next_line(words,nwords)
    INTEGER,INTENT(OUT) :: nwords
    CHARACTER(10000),INTENT(OUT) :: words(:)
    CHARACTER(10000) :: line
    INTEGER :: ios
    DO
      READ(local_unit,'(A10000)',IOSTAT=ios)line
      IF( ios .NE. 0)CALL raise_fatal_error('end of xs file was reached before all data/materials were found')
      line=TRIM(ADJUSTL(line))
      !finding uncommented line that isn't empty
      IF(line(1:1) .NE. '!' .AND. line .NE. '')THEN
        !ignore commented portions of line
        CALL parse(line,'!',words,nwords)
        line=TRIM(ADJUSTL(words(1)))
        CALL parse(line,' ',words,nwords)
        EXIT
      ENDIF
    ENDDO
  ENDSUBROUTINE get_next_line


END MODULE read_cross_section_module
