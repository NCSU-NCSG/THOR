MODULE dump_inguess_module
  !***********************************************************************
  !
  ! This module contains subroutines to create and read restart(dump) files
  !
  !***********************************************************************

  ! Use derived-type modules

  USE types
  USE parameter_types
  USE filename_types
  USE vector_types
  USE cross_section_types
  USE geometry_types
  USE angle_types
  USE multindex_types
  USE global_variables
  USE error_module

  IMPLICIT NONE

CONTAINS


  SUBROUTINE inguess_eig(flux,keff,flag)
    !*********************************************************************
    !
    ! Subroutine inguess_eig reads a file that is created by either
    ! dump_PI or dump_jfnk and uses it as initial guess
    !
    !*********************************************************************

    ! Pass arguments

    INTEGER(kind=li) :: flag
    REAL(kind=d_t) :: keff
    REAL(kind=d_t) :: flux(num_moments_v,namom,num_cells,egmax,niter)

    ! local variables

    INTEGER(kind=li) :: eg,i,l,p
    INTEGER(kind=li) :: egmax_t,num_cells_t,namom_t,num_moments_v_t
    LOGICAL :: existence

    ! open file
    INQUIRE(file=TRIM(inguess_file), exist=existence)
    IF (existence .EQV. .TRUE.) THEN
      OPEN(UNIT = 66, FILE = inguess_file, STATUS = "OLD", ACTION = "READ",&
            form='unformatted')
      ! read inguess file
      READ(66) egmax_t,num_cells_t,namom_t,num_moments_v_t
      ! check if data in file and data from stdin are identical
      IF(egmax_t .NE. egmax .OR. num_cells_t .NE. num_cells .OR.    &
            num_moments_v_t .NE. num_moments_v .OR. namom_t .NE. namom ) THEN
        CALL raise_fatal_error("Reading initial guess file failed")
      END IF
      DO eg =1,egmax
        DO i=1,num_cells
          DO p=1,namom
            DO l=1,num_moments_v
              READ(66) flux(l,p,i,eg,flag)
            END DO
          END DO
        END DO
      END DO
      READ(66) keff
      ! close file

      CLOSE(unit=66)
      WRITE(stdout_unit,*)
      WRITE(stdout_unit,*) 'Initial guess file ',TRIM(inguess_file),' was successfully read'
      WRITE(log_unit,*)
      WRITE(log_unit,*) 'Initial guess file ',TRIM(inguess_file),' was successfully read'

    ELSE
      WRITE(stdout_unit,*)
      WRITE(stdout_unit,*) "Attempted to read initial guess file ",TRIM(inguess_file)," but could not find it."
      WRITE(stdout_unit,*)
      WRITE(log_unit,*)
      WRITE(log_unit,*) "Attempted to read initial guess file ",TRIM(inguess_file)," but could not find it."
      WRITE(log_unit,*)
      flux(:,:,:,:,niter)=zero
      DO eg =1,egmax
        DO i=1,num_cells
          flux(1,1,i,eg,niter)=one
        END DO
      END DO
      keff = 1.0_d_t
    END IF

  END SUBROUTINE inguess_eig

  SUBROUTINE dump_jfnk(flux,keff)
    !*********************************************************************
    !
    ! Subroutine dump_jfnk creates restart file after each Newton iteration
    !
    !*********************************************************************

    ! Pass arguments

    REAL(kind=d_t) :: keff
    REAL(kind=d_t) :: flux(num_moments_v,namom,num_cells,egmax,niter)

    ! local variables

    INTEGER(kind=li) :: eg,i,l,p
    LOGICAL :: existence

    ! open file
    INQUIRE(file=TRIM(dump_file), exist=existence)
    IF (existence .EQV. .TRUE.) THEN
      OPEN(UNIT = 56, FILE = dump_file, STATUS = "OLD", ACTION = "WRITE",&
            form='unformatted')
    ELSE
      OPEN(UNIT = 56, FILE = dump_file, STATUS = "NEW", ACTION ="WRITE",&
            form='unformatted')
    END IF

    ! write dump file
    WRITE(56) egmax,num_cells,namom,num_moments_v
    DO eg =1,egmax
      DO i=1,num_cells
        DO p=1,namom
          DO l=1,num_moments_v
            WRITE(56) flux(l,p,i,eg,1)
          END DO
        END DO
      END DO
    END DO
    WRITE(56) keff

    ! close file

    CLOSE(unit=56)

  END SUBROUTINE dump_jfnk

  SUBROUTINE dump_PI(flux,keff)
    !*********************************************************************
    !
    ! Subroutine dump_PI creates restart filse after each power iteration
    !
    !*********************************************************************

    ! Pass arguments

    REAL(kind=d_t) :: keff
    REAL(kind=d_t) :: flux(num_moments_v,namom,num_cells,egmax,niter)

    ! local variables

    INTEGER(kind=li) :: eg,i,l,p
    LOGICAL :: existence

    ! open file
    INQUIRE(file=TRIM(dump_file), exist=existence)
    IF (existence .EQV. .TRUE.) THEN
      OPEN(UNIT = 56, FILE = dump_file, STATUS = "OLD", ACTION = "WRITE",&
            form='unformatted')
    ELSE
      OPEN(UNIT = 56, FILE = dump_file, STATUS = "NEW", ACTION ="WRITE",&
            form='unformatted')
    END IF

    ! write dump file
    WRITE(56) egmax,num_cells,namom,num_moments_v
    DO eg =1,egmax
      DO i=1,num_cells
        DO p=1,namom
          DO l=1,num_moments_v
            WRITE(56) flux(l,p,i,eg,niter)
          END DO
        END DO
      END DO
    END DO
    WRITE(56) keff

    ! close file

    CLOSE(unit=56)

  END SUBROUTINE dump_PI

END MODULE dump_inguess_module
