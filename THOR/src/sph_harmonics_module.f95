MODULE sph_harmonics_module
  !*********************************************************************
  !
  ! Module spherical harmonics constructs real spherical harmonic
  ! expansion coefficient
  !
  !*********************************************************************
  ! User derived-type modules

  USE types
  USE parameter_types
  USE vector_types
  USE angle_types
  USE multindex_types

  ! Use modules that pertain setting up problem

  USE termination_module

  IMPLICIT NONE

CONTAINS

  SUBROUTINE spherical_harmonics(sord,nang,quad,max_p,Y)
    !*********************************************************************
    !
    ! Subroutine spherical harmonics constructs real spherical harmonic
    ! expansion coefficient
    !
    !*********************************************************************
    ! Pass input parameters

    INTEGER(kind=li), INTENT(in) :: sord, nang

    ! Pass angular derived type used globally

    TYPE(ordinate), DIMENSION(nang), INTENT(in) :: &
          quad

    ! Define temporary variables

    INTEGER(kind=li) :: q, octant,l,m,indx
    INTEGER(kind=li), INTENT(inout) :: max_p
    REAL(kind=d_t), DIMENSION(nang,8,max_p) :: Y
    REAL(kind=d_t) :: mu,eta,xi
    TYPE(vector) :: omega

    ! Do even contributions
    DO l=0,sord
      DO m=0,l
        indx=1_li + m + (l+1_li)*l/2_li
        DO octant=1,8
          DO q=1, nang
            IF(octant == 1)THEN
              mu  = quad(q)%mu%x1
              eta = quad(q)%mu%x2
              xi  = quad(q)%mu%x3
            ELSEIF(octant == 2)THEN
              mu  =-quad(q)%mu%x1
              eta = quad(q)%mu%x2
              xi  = quad(q)%mu%x3
            ELSEIF(octant == 3)THEN
              mu  =-quad(q)%mu%x1
              eta =-quad(q)%mu%x2
              xi  = quad(q)%mu%x3
            ELSEIF(octant == 4)THEN
              mu  = quad(q)%mu%x1
              eta =-quad(q)%mu%x2
              xi  = quad(q)%mu%x3
            ELSEIF(octant == 5)THEN
              mu  = quad(q)%mu%x1
              eta = quad(q)%mu%x2
              xi  = -quad(q)%mu%x3
            ELSEIF(octant == 6)THEN
              mu  =-quad(q)%mu%x1
              eta = quad(q)%mu%x2
              xi  = -quad(q)%mu%x3
            ELSEIF(octant == 7)THEN
              mu  =-quad(q)%mu%x1
              eta =-quad(q)%mu%x2
              xi  = -quad(q)%mu%x3
            ELSE
              mu  = quad(q)%mu%x1
              eta =-quad(q)%mu%x2
              xi  = -quad(q)%mu%x3
            END IF
            Y(q,octant,indx)=y_e(l,m,mu,eta,xi)
          END DO
        END DO
      END DO
    END DO
    ! Odd contributions
    DO l=1,sord
      DO m=1,l
        indx= neven + m  + (l-1_li)*l/2_li
        DO octant=1,8
          DO q=1, nang
            IF(octant == 1)THEN
              mu  = quad(q)%mu%x1
              eta = quad(q)%mu%x2
              xi  = quad(q)%mu%x3
            ELSEIF(octant == 2)THEN
              mu  =-quad(q)%mu%x1
              eta = quad(q)%mu%x2
              xi  = quad(q)%mu%x3
            ELSEIF(octant == 3)THEN
              mu  =-quad(q)%mu%x1
              eta =-quad(q)%mu%x2
              xi  = quad(q)%mu%x3
            ELSEIF(octant == 4)THEN
              mu  = quad(q)%mu%x1
              eta =-quad(q)%mu%x2
              xi  = quad(q)%mu%x3
            ELSEIF(octant == 5)THEN
              mu  = quad(q)%mu%x1
              eta = quad(q)%mu%x2
              xi  = -quad(q)%mu%x3
            ELSEIF(octant == 6)THEN
              mu  =-quad(q)%mu%x1
              eta = quad(q)%mu%x2
              xi  = -quad(q)%mu%x3
            ELSEIF(octant == 7)THEN
              mu  =-quad(q)%mu%x1
              eta =-quad(q)%mu%x2
              xi  = -quad(q)%mu%x3
            ELSE
              mu  = quad(q)%mu%x1
              eta =-quad(q)%mu%x2
              xi  = -quad(q)%mu%x3
            END IF
            Y(q,octant,indx)=y_o(l,m,mu,eta,xi)
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE spherical_harmonics

  FUNCTION y_e(l,m,mu,eta,xi)
    INTEGER(kind=li), INTENT(in) :: l,m
    REAL(kind=d_t), INTENT(in) :: mu,eta,xi
    REAL(kind=d_t) :: y_e
    REAL(kind=d_t) :: phi

    !     phi = acos(eta/sqrt(1.0_d_t-mu**2))
    phi = ATAN2(eta,mu)
    y_e = SQRT(REAL(2*l+1,d_t)*factorial_d_t(l-m)/factorial_d_t(l+m))*plegendre(l,m,xi)*&
          COS(REAL(m,d_t)*phi)*(-1.0_d_t)**m

  END FUNCTION y_e

  FUNCTION y_o(l,m,mu,eta,xi)
    INTEGER(kind=li), INTENT(in) :: l,m
    REAL(kind=d_t), INTENT(in) :: mu,eta,xi
    REAL(kind=d_t) :: y_o
    REAL(kind=d_t) :: phi

    !     phi = acos(eta/sqrt(1.0_d_t-mu**2))
    phi = ATAN2(eta,mu)
    y_o = SQRT(REAL(2*l+1,d_t)*factorial_d_t(l-m)/factorial_d_t(l+m))*plegendre(l,m,xi)*&
          SIN(REAL(m,d_t)*phi)*(-1.0_d_t)**m

  END FUNCTION y_o

  ! adopted from numerical recipes plegendre.h
  FUNCTION plegendre(l, m, x)
    INTEGER(kind=li), INTENT(in) :: l,m
    REAL(kind=d_t), INTENT(in) :: x
    REAL(kind=d_t) :: plegendre

    INTEGER(kind=li) :: i, j, ll
    REAL(kind=d_t) :: fact, oldfact, pll, pmm, pmmp1, omx2, prod

    REAL(kind=d_t), PARAMETER :: PI = 3.141592653589793_d_t

    IF (m < 0 .or. m > l .or. abs(x) > 1.0_d_t) THEN
      CALL stop_thor(1000_li, "Bad arguments in routine plegendre")
    END IF

    pmm = 1.0_d_t
    IF (m > 0) THEN
      omx2 = (1.0_d_t - x) * (1.0_d_t + x)
      fact = 1.0_d_t
      DO i = 1, m
        pmm = pmm * omx2 * fact / (fact + 1.0_d_t)
        fact = fact + 2.0_d_t
      END DO
    END IF

    pmm = sqrt((2.0_d_t * m + 1.0_d_t) * pmm / (4.0_d_t * PI));
    IF (m == 1) pmm = -pmm

    IF (l == m) THEN
      plegendre = pmm
    ELSE
      pmmp1 = x * sqrt(2.0_d_t * m + 3.0_d_t) * pmm;
      IF (l == m + 1) THEN
        plegendre = pmmp1
      ELSE
        oldfact=sqrt(2.0_d_t * m + 3.0_d_t);
        DO ll = m + 2, l
          fact = sqrt((4.0_d_t * ll * ll - 1.0_d_t)/(ll * ll - m * m));
          pll = (x * pmmp1 - pmm / oldfact) * fact;
          oldfact = fact
          pmm = pmmp1
          pmmp1 = pll
        END DO
        plegendre = pll
      END IF
    END IF

    prod = 1.0_d_t;
    DO j = l - m + 1, l + m
      prod = prod * j
    END DO
    plegendre = plegendre * sqrt(4.0_d_t * PI * prod / (2.0_d_t * l + 1.0_d_t))
  END FUNCTION

END MODULE sph_harmonics_module
