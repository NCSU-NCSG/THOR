MODULE read_cross_section_module
  !***********************************************************************
  !
  ! Read cross-section module  contains all subroutines needed to read
  ! cross-section file
  !
  !***********************************************************************

  USE types
  USE parameter_types
  USE filename_types
  USE cross_section_types
  USE global_variables
  USE termination_module

  IMPLICIT NONE

CONTAINS

  SUBROUTINE read_xs
    !*********************************************************************
    !
    ! Subroutine reads cross-section in 'unique' ahot format
    ! (could be adapted for other formats, of course)
    !
    !*********************************************************************

    INTEGER(kind=li) :: alloc_stat, e1, order, eg_to, eg_from,l,m

    ! Open and read mesh file

    OPEN(unit=11,file=TRIM(cross_section_filename),status='old',action='read')

    READ(11,*) num_mat

    ! Allocate cross-section arrays and check if enough memory is available

    ALLOCATE(xs_mat(0:num_mat),chi(0:num_mat,egmax),&
          eg_bounds(0:num_mat,egmax+1),fiss(0:num_mat,egmax),&
          nu(0:num_mat,egmax),sigma_t(0:num_mat,egmax),tsigs(0:num_mat,egmax),&
          sigma_scat(0:num_mat,xs_ord+1,egmax,egmax),&
          stat=alloc_stat)
    IF(alloc_stat /= 0) CALL stop_thor(2_li)

    ! Set most_thermal to no upscattering

    most_thermal=egmax+1

    ! Read material total & scattering cross-section

    DO m=1, num_mat
      READ(11,*) xs_mat(m)%mat
      IF(multiplying .NE. 0)THEN
        READ(11,*) (chi(xs_mat(m)%mat,e1)%xs,e1=1,egmax)
        READ(11,*) (eg_bounds(xs_mat(m)%mat,e1)%xs,e1=1,egmax)
        eg_bounds(xs_mat(m)%mat,egmax+1)%xs=0.0_d_t
        READ(11,*) (fiss(xs_mat(m)%mat,e1)%xs,e1=1,egmax)
        READ(11,*) (nu(xs_mat(m)%mat,e1)%xs,e1=1,egmax)
      ELSE
        DO e1 = 1, egmax
          chi(xs_mat(m)%mat,e1)%xs = zero
          eg_bounds(xs_mat(m)%mat,e1)%xs = zero
          fiss(xs_mat(m)%mat,e1)%xs = zero
          nu(xs_mat(m)%mat,e1)%xs = zero
        END DO
        eg_bounds(xs_mat(m)%mat,egmax+1)%xs = zero
      END IF
      READ(11,*) (sigma_t(xs_mat(m)%mat,e1)%xs,e1=1,egmax)

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
            READ(11,*) (sigma_scat(xs_mat(m)%mat,order,&
                  eg_to,eg_from)%xs,eg_from=1,eg_to)
          END DO
        END DO
      ELSE
        DO order=1, xs_ord+1
          DO eg_to=1,egmax
            READ(11,*) (sigma_scat(xs_mat(m)%mat,order,&
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

    ! Close mesh file

    CLOSE(11)

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

  END SUBROUTINE read_xs


END MODULE read_cross_section_module
