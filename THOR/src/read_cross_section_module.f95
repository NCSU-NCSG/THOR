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
  USE global_variables
  USE termination_module
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

    INTEGER(kind=li) :: rank,ios,mpi_err
    CHARACTER(200) :: xs_format

    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)
    local_unit=100+rank
    ! Open and read mesh file

    OPEN(unit=local_unit,file=TRIM(cross_section_filename),status='old',action='read',IOSTAT=ios)
    IF(ios .NE. 0)THEN
      WRITE(*,*)'error opening ',TRIM(cross_section_filename)
      STOP 'fatal error'
    ENDIF

    ! Set most_thermal to no upscattering

    most_thermal=egmax+1

    DO
      READ(local_unit,*,IOSTAT=ios)xs_format
      IF(ios .NE. 0)THEN
        !legacy original THOR xs format with no format indicater
        REWIND(local_unit)
        CALL xs_read_legacyv0()
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
  END SUBROUTINE read_xs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE xs_read_legacyv0()

    INTEGER(kind=li) :: alloc_stat, e1, order, eg_to, eg_from,l,m

    READ(local_unit,*) num_mat

    ! Allocate cross-section arrays and check if enough memory is available

    ALLOCATE(xs_mat(0:num_mat),chi(0:num_mat,egmax),&
          eg_bounds(0:num_mat,egmax+1),fiss(0:num_mat,egmax),&
          nu(0:num_mat,egmax),sigma_t(0:num_mat,egmax),tsigs(0:num_mat,egmax),&
          sigma_scat(0:num_mat,xs_ord+1,egmax,egmax),&
          stat=alloc_stat)
    IF(alloc_stat /= 0) CALL stop_thor(2_li)

    ! Read material total & scattering cross-section

    DO m=1, num_mat
      READ(local_unit,*) xs_mat(m)%mat
      IF(multiplying .NE. 0)THEN
        READ(local_unit,*) (chi(xs_mat(m)%mat,e1)%xs,e1=1,egmax)
        READ(local_unit,*) (eg_bounds(xs_mat(m)%mat,e1)%xs,e1=1,egmax)
        eg_bounds(xs_mat(m)%mat,egmax+1)%xs=0.0_d_t
        READ(local_unit,*) (fiss(xs_mat(m)%mat,e1)%xs,e1=1,egmax)
        READ(local_unit,*) (nu(xs_mat(m)%mat,e1)%xs,e1=1,egmax)
      ELSE
        DO e1 = 1, egmax
          chi(xs_mat(m)%mat,e1)%xs = zero
          eg_bounds(xs_mat(m)%mat,e1)%xs = zero
          fiss(xs_mat(m)%mat,e1)%xs = zero
          nu(xs_mat(m)%mat,e1)%xs = zero
        END DO
        eg_bounds(xs_mat(m)%mat,egmax+1)%xs = zero
      END IF
      READ(local_unit,*) (sigma_t(xs_mat(m)%mat,e1)%xs,e1=1,egmax)

      ! Initialize scattering matrix

      DO order=1, xs_ord+1
        DO eg_to=1,egmax
          DO eg_from=1,egmax
            sigma_scat(xs_mat(m)%mat,order,eg_to,eg_from)%xs=zero
          END DO
        END DO
      END DO

      ! Read scattering matrix, note: thermal groups are separated from fast
      ! groups but old cross section format remains valid

      IF(upscattering.EQ.0) THEN
        DO order=1, xs_ord+1
          DO eg_to=1,egmax
            READ(local_unit,*) (sigma_scat(xs_mat(m)%mat,order,&
                  eg_to,eg_from)%xs,eg_from=1,eg_to)
          END DO
        END DO
      ELSE
        DO order=1, xs_ord+1
          DO eg_to=1,egmax
            READ(local_unit,*) (sigma_scat(xs_mat(m)%mat,order,&
                  eg_to,eg_from)%xs,eg_from=1,egmax)
          END DO
        END DO
      END IF

      ! Determine most_thermal

      IF(upscattering .NE. 0) THEN
        DO eg_to=1,egmax
          DO eg_from=eg_to+1,egmax
            IF( ABS(sigma_scat(xs_mat(m)%mat,1,eg_to,eg_from)%xs) > 2.24E-16_d_t .AND. most_thermal>eg_to ) most_thermal=eg_to
          END DO
        END DO
      END IF

      ! Compute tsigs
      DO eg_from=1,egmax
        tsigs(xs_mat(m)%mat,eg_from)%xs=0.0_d_t
        DO eg_to=1,egmax
          tsigs(xs_mat(m)%mat,eg_from)%xs=tsigs(xs_mat(m)%mat,eg_from)%xs+sigma_scat(xs_mat(m)%mat,1,eg_to,eg_from)%xs
        END DO
      END DO

    END DO

    ! set neven

    neven=1_li+scatt_ord+(scatt_ord+1_li)*scatt_ord/2_li

    ! set scat_mult

    ALLOCATE(scat_mult(0:scatt_ord,0:scatt_ord))
    scat_mult=0.0_d_t
    IF(scat_mult_flag.EQ.1_li) THEN
      DO l=0,scatt_ord
        DO m=0,l
          IF(m.EQ.0_li) THEN
            scat_mult(l,m)=1.0_d_t/REAL(2_li*l+1_li,d_t)
          ELSE
            scat_mult(l,m)=2.0_d_t/REAL(2_li*l+1_li,d_t)
          END IF
        END DO
      END DO
    ELSE
      DO l=0,scatt_ord
        DO m=0,l
          IF(m.EQ.0_li) THEN
            scat_mult(l,m)=1.0_d_t
          ELSE
            scat_mult(l,m)=2.0_d_t
          END IF
        END DO
      END DO
    END IF

  END SUBROUTINE xs_read_legacyv0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!currently supports max of 1000 groups
  SUBROUTINE xs_read_current()
    CHARACTER(10000) :: line
    CHARACTER(10000) :: words(1000)
    INTEGER :: nwords,ios,mat_num,i,g,alloc_stat,gp,j

    num_mat=0
    egmax=0
    xs_ord=0
    !get the specification data
    DO
      READ(local_unit,'(A10000)')line
      line=TRIM(ADJUSTL(line))
      CALL parse(line,' ',words,nwords)
      words(1)=TRIM(ADJUSTL(words(1)))
      IF(words(1) .EQ. 'THOR_XS_V1')THEN
        !get the material, group, and order amount
        words(2)=TRIM(ADJUSTL(words(2)))
        words(3)=TRIM(ADJUSTL(words(3)))
        words(4)=TRIM(ADJUSTL(words(4)))
        READ(words(2),*)num_mat
        READ(words(3),*)egmax
        READ(words(4),*)xs_ord
        REWIND(local_unit)
        CALL SLEEP(5)
        EXIT
      ENDIF
    ENDDO

    ! Allocate cross-section arrays and check if enough memory is available

    ALLOCATE(xs_mat(0:num_mat),chi(0:num_mat,egmax),&
          eg_bounds(0:num_mat,egmax+1),fiss(0:num_mat,egmax),&
          nu(0:num_mat,egmax),sigma_t(0:num_mat,egmax),tsigs(0:num_mat,egmax),&
          sigma_scat(0:num_mat,xs_ord+1,egmax,egmax),mat_name(0:num_mat),&
          stat=alloc_stat)
    IF(alloc_stat /= 0) CALL stop_thor(2_li)
    !set everything to 0
    xs_mat(:)%mat=0
    chi(:,:)%xs=0
    eg_bounds(:,:)%xs=0
    fiss(:,:)%xs=0
    nu(:,:)%xs=0
    sigma_t(:,:)%xs=0
    tsigs(:,:)%xs=0
    sigma_scat(:,:,:,:)%xs=0
    mat_name(:)=''

    DO i=1,num_mat
      DO
        READ(local_unit,'(A10000)',IOSTAT=ios)line
        IF(ios .NE. 0)STOP 'end of xs file was reached before all materials were found'
        line=TRIM(ADJUSTL(line))
        !uncommented line
        IF(line(1:1) .NE. '!' .AND. line .NE. '')THEN
          !ignore commented portions of line
          CALL parse(line,'!',words,nwords)
          line=TRIM(ADJUSTL(words(1)))
          !separate into data
          CALL parse(line,' ',words,nwords)
          words(1)=TRIM(ADJUSTL(words(1)))
          IF(words(1) .EQ. 'id')THEN
            READ(words(2),*)mat_num
            xs_mat(mat_num)%mat=mat_num
            words(4)=TRIM(ADJUSTL(words(4)))
            !get or assign mat name
            IF(words(4) .NE. '')THEN
              mat_name(mat_num)=words(4)
            ELSE
              WRITE(mat_name(mat_num),'(A,I0)')'mat_',mat_num
            ENDIF
            !read in the fission spectrum
            DO
              READ(local_unit,'(A10000)',IOSTAT=ios)line
              IF(ios .NE. 0)STOP 'end of xs file was reached before all data was found'
              line=TRIM(ADJUSTL(line))
              !uncommented line
              IF(line(1:1) .NE. '!' .AND. line .NE. '')THEN
                !ignore commented portions of line
                CALL parse(line,'!',words,nwords)
                line=TRIM(ADJUSTL(words(1)))
                !separate into data
                CALL parse(line,' ',words,nwords)
                IF(nwords .NE. egmax)STOP 'bad amount of xs data on line'
                DO g=1,egmax
                  READ(words(g),*)chi(mat_num,g)%xs
                ENDDO
                EXIT
              ENDIF
            ENDDO
            !read in the energy bounds
            DO
              READ(local_unit,'(A10000)',IOSTAT=ios)line
              IF(ios .NE. 0)STOP 'end of xs file was reached before all data was found'
              line=TRIM(ADJUSTL(line))
              !uncommented line
              IF(line(1:1) .NE. '!' .AND. line .NE. '')THEN
                !ignore commented portions of line
                CALL parse(line,'!',words,nwords)
                line=TRIM(ADJUSTL(words(1)))
                !separate into data
                CALL parse(line,' ',words,nwords)
                IF(nwords .NE. egmax)STOP 'bad amount of xs data on line'
                DO g=1,egmax
                  READ(words(g),*)eg_bounds(mat_num,g)%xs
                ENDDO
                eg_bounds(mat_num,egmax+1)%xs=0.0
                EXIT
              ENDIF
            ENDDO
            !read in SigmaF
            DO
              READ(local_unit,'(A10000)',IOSTAT=ios)line
              IF(ios .NE. 0)STOP 'end of xs file was reached before all data was found'
              line=TRIM(ADJUSTL(line))
              !uncommented line
              IF(line(1:1) .NE. '!' .AND. line .NE. '')THEN
                !ignore commented portions of line
                CALL parse(line,'!',words,nwords)
                line=TRIM(ADJUSTL(words(1)))
                !separate into data
                CALL parse(line,' ',words,nwords)
                IF(nwords .NE. egmax)STOP 'bad amount of xs data on line'
                DO g=1,egmax
                  READ(words(g),*)fiss(mat_num,g)%xs
                ENDDO
                EXIT
              ENDIF
            ENDDO
            !read in nu
            DO
              READ(local_unit,'(A10000)',IOSTAT=ios)line
              IF(ios .NE. 0)STOP 'end of xs file was reached before all data was found'
              line=TRIM(ADJUSTL(line))
              !uncommented line
              IF(line(1:1) .NE. '!' .AND. line .NE. '')THEN
                !ignore commented portions of line
                CALL parse(line,'!',words,nwords)
                line=TRIM(ADJUSTL(words(1)))
                !separate into data
                CALL parse(line,' ',words,nwords)
                IF(nwords .NE. egmax)STOP 'bad amount of xs data on line'
                DO g=1,egmax
                  READ(words(g),*)nu(mat_num,g)%xs
                ENDDO
                EXIT
              ENDIF
            ENDDO
            !read in total/transport xs
            DO
              READ(local_unit,'(A10000)',IOSTAT=ios)line
              IF(ios .NE. 0)STOP 'end of xs file was reached before all data was found'
              line=TRIM(ADJUSTL(line))
              !uncommented line
              IF(line(1:1) .NE. '!' .AND. line .NE. '')THEN
                !ignore commented portions of line
                CALL parse(line,'!',words,nwords)
                line=TRIM(ADJUSTL(words(1)))
                !separate into data
                CALL parse(line,' ',words,nwords)
                IF(nwords .NE. egmax)STOP 'bad amount of xs data on line'
                DO g=1,egmax
                  READ(words(g),*)sigma_t(mat_num,g)%xs
                ENDDO
                EXIT
              ENDIF
            ENDDO
            !read in scattering xs format gp->g
            DO j=1,xs_ord+1
              DO g=1,egmax
                DO
                  READ(local_unit,'(A10000)',IOSTAT=ios)line
                  IF(ios .NE. 0)STOP 'end of xs file was reached before all data was found'
                  line=TRIM(ADJUSTL(line))
                  !uncommented line
                  IF(line(1:1) .NE. '!' .AND. line .NE. '')THEN
                    !ignore commented portions of line
                    CALL parse(line,'!',words,nwords)
                    line=TRIM(ADJUSTL(words(1)))
                    !separate into data
                    CALL parse(line,' ',words,nwords)
                    IF(nwords .NE. egmax)STOP 'bad amount of xs data on line'
                    DO gp=1,egmax
                      READ(words(gp),*)sigma_scat(mat_num,j,g,gp)%xs
                    ENDDO
                    EXIT
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
            EXIT
          ENDIF
        ENDIF
      ENDDO
    ENDDO

    ! Determine most_thermal
    DO i=1,num_mat
      write(*,*)i,xs_mat(i)%mat
      IF(upscattering .NE. 0) THEN
        DO g=1,egmax
          DO gp=g+1,egmax
            IF( ABS(sigma_scat(xs_mat(i)%mat,1,g,gp)%xs) > 2.24E-16_d_t .AND. most_thermal>g ) &
              most_thermal=g
          END DO
        END DO
      END IF

      ! Compute tsigs
      DO gp=1,egmax
        tsigs(xs_mat(i)%mat,gp)%xs=0.0_d_t
        DO g=1,egmax
          tsigs(xs_mat(i)%mat,gp)%xs=tsigs(xs_mat(i)%mat,gp)%xs+sigma_scat(xs_mat(i)%mat,1,g,gp)%xs
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

  END SUBROUTINE xs_read_current

  ! SUBROUTINE get_next_line()
  ! ENDSUBROUTINE


END MODULE read_cross_section_module
