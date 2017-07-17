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

    y_e = SQRT(REAL(2*l+1,d_t)*factorial_d_t(l-m)/factorial_d_t(l+m))*aLegendreP(l,m,xi)*&
          COS(REAL(m,d_t)*phi)*(-1.0_d_t)**m

  END FUNCTION y_e

  FUNCTION y_o(l,m,mu,eta,xi)
    INTEGER(kind=li), INTENT(in) :: l,m
    REAL(kind=d_t), INTENT(in) :: mu,eta,xi
    REAL(kind=d_t) :: y_o
    REAL(kind=d_t) :: phi

    !     phi = acos(eta/sqrt(1.0_d_t-mu**2))
    phi = ATAN2(eta,mu)

    y_o = SQRT(REAL(2*l+1,d_t)*factorial_d_t(l-m)/factorial_d_t(l+m))*aLegendreP(l,m,xi)*&
          SIN(REAL(m,d_t)*phi)*(-1.0_d_t)**m

  END FUNCTION y_o

  FUNCTION aLegendreP(l,m,x)
    INTEGER(kind=li), INTENT(in) :: l,m
    REAL(kind=d_t), INTENT(in) :: x
    REAL(kind=d_t) :: aLegendreP

    IF(l == 0_li) THEN
      IF(m == 0_li) THEN
        aLegendreP=1.0_d_t
      ELSE
        CALL stop_thor(32_li)
      END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ELSE IF(l == 1_li) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF(m == -1_li) THEN
        aLegendreP=0.5_d_t*SQRT(1.0_d_t-x*x)
      ELSE IF (m == 0_li) THEN
        aLegendreP=x
      ELSE IF (m == 1_li) THEN
        aLegendreP=-SQRT(1.0_d_t-x*x)
      ELSE
        CALL stop_thor(32_li)
      END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ELSE IF(l == 2_li) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF(m == -2_li) THEN
        aLegendreP=0.125_d_t*(1.0_d_t-x*x)
      ELSE IF(m == -1_li) THEN
        aLegendreP=0.5_d_t*x*SQRT(1.0_d_t-x*x)
      ELSE IF(m == 0_li) THEN
        aLegendreP=0.5_d_t*(3.0_d_t*x*x-1.0_d_t)
      ELSE IF(m == 1_li) THEN
        aLegendreP=-3.0_d_t*x*SQRT(1.0_d_t-x*x)
      ELSE IF(m == 2_li) THEN
        aLegendreP=-3.0_d_t*(x*x-1.0_d_t)
      ELSE
        CALL stop_thor(32_li)
      END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ELSE IF(l == 3_li) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF(m == -3_li) THEN
        aLegendreP=0.02083333333333333_d_t*SQRT((1.0_d_t-x*x)**3)
      ELSE IF(m == -2_li) THEN
        aLegendreP=-0.125_d_t*x*(x*x-1.0_d_t)
      ELSE IF(m == -1_li) THEN
        aLegendreP=0.125_d_t*SQRT(1.0_d_t-x*x)*(5.0_d_t*x*x-1.0_d_t)
      ELSE IF(m == 0_li) THEN
        aLegendreP=0.5_d_t*(-3.0_d_t*x+5.0_d_t*x*x*x)
      ELSE IF(m == 1_li) THEN
        aLegendreP=-1.5_d_t*SQRT(1.0_d_t-x*x)*(5.0_d_t*x*x-1.0_d_t)
      ELSE IF(m == 2_li) THEN
        aLegendreP=-15.0_d_t*x*(x*x-1.0_d_t)
      ELSE IF(m == 3_li) THEN
        aLegendreP=-15.0_d_t**SQRT((1.0_d_t-x*x)**3)
      ELSE
        CALL stop_thor(32_li)
      END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ELSE IF(l == 4_li) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF(m == -4_li) THEN
        aLegendreP=0.002604166666666667_d_t*(-1.0_d_t+x*x)**2
      ELSE IF(m == -3_li) THEN
        aLegendreP=0.02083333333333333_d_t*x*SQRT((1.0_d_t-x*x)**3)
      ELSE IF(m == -2_li) THEN
        aLegendreP=-0.02083333333333333_d_t*(-1.0_d_t+x*x)*(-1.0_d_t+7.0_d_t*x*x)
      ELSE IF(m == -1_li) THEN
        aLegendreP=0.125_d_t*SQRT(1.0_d_t-x*x)*(-3.0_d_t*x+7.0_d_t*x*x*x)
      ELSE IF(m == 0_li) THEN
        aLegendreP=0.125_d_t*(3.0_d_t-30.0_d_t*x*x+35.0_d_t*x*x*x*x)
      ELSE IF(m == 1_li) THEN
        aLegendreP=-2.5_d_t*SQRT(1.0_d_t-x*x)*(-3.0_d_t*x+7.0_d_t*x*x*x)
      ELSE IF(m == 2_li) THEN
        aLegendreP=-7.5_d_t*(-1.0_d_t+x*x)*(-1.0_d_t+7.0_d_t*x*x)
      ELSE IF(m == 3_li) THEN
        aLegendreP=-105.0_d_t*x*SQRT((1.0_d_t-x*x)**3)
      ELSE IF(m == 4_li) THEN
        aLegendreP=105.0_d_t*(-1.0_d_t+x*x)**2
      ELSE
        CALL stop_thor(32_li)
      END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ELSE IF(l == 5_li) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF(m == -5_li) THEN
        aLegendreP=0.0002604166666666667_d_t*SQRT((1.0_d_t-x*x)**5)
      ELSE IF(m == -4_li) THEN
        aLegendreP=0.002604166666666667_d_t*x*(-1.0_d_t+x*x)**2
      ELSE IF(m == -3_li) THEN
        aLegendreP=-0.002604166666666667_d_t*SQRT(1.0_d_t-x*x)*(-1.0_d_t+x*x)*(-1.0_d_t+9.0_d_t*x*x)
      ELSE IF(m == -2_li) THEN
        aLegendreP=-0.0625_d_t*(-1.0_d_t+x*x)*(-x+3.0_d_t*x*x*x)
      ELSE IF(m == -1_li) THEN
        aLegendreP=0.0625_d_t*SQRT(1.0_d_t-x*x)*(1.0_d_t-14.0_d_t*x*x+21.0_d_t*x*x*x*x)
      ELSE IF(m == 0_li) THEN
        aLegendreP=0.125_d_t*(15.0_d_t*x-70.0_d_t*x*x*x+63.0_d_t*x*x*x*x*x)
      ELSE IF(m == 1_li) THEN
        aLegendreP=-1.875_d_t*SQRT(1.0_d_t-x*x)*(1.0_d_t-14.0_d_t*x*x+21.0_d_t*x*x*x*x)
      ELSE IF(m == 2_li) THEN
        aLegendreP=-52.5_d_t*(-1.0_d_t+x*x)*(-x+3.0_d_t*x*x*x)
      ELSE IF(m == 3_li) THEN
        aLegendreP=52.5_d_t*SQRT(1.0_d_t-x*x)*(-1.0_d_t+x*x)*(-1.0_d_t+9.0_d_t*x*x)
      ELSE IF(m == 4_li) THEN
        aLegendreP=945.0_d_t*x*(-1.0_d_t+x*x)**2
      ELSE IF(m == 5_li) THEN
        aLegendreP=-945.0_d_t*SQRT((1.0_d_t-x*x)**5)
      ELSE
        CALL stop_thor(32_li)
      END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ELSE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL stop_thor(33_li)
    END IF

  END FUNCTION aLegendreP

  FUNCTION LegendreP(l,x)
    INTEGER(kind=li), INTENT(in) :: l
    REAL(kind=d_t), INTENT(in) :: x
    REAL(kind=d_t) :: LegendreP

    IF(l == 0_li)THEN
      LegendreP=1
    ELSEIF(l == 1_li)THEN
      LegendreP=x
    ELSEIF(l == 2_li)THEN
      LegendreP=0.5_d_t*(3.0_d_t*x*x-1)
    ELSEIF(l == 3_li)THEN
      LegendreP=0.5_d_t*(5.0_d_t*x*x*x-3.0_d_t*x)
    ELSEIF(l == 4_li)THEN
      LegendreP=0.125_d_t*(35.0_d_t*x*x*x*x-30.0_d_t*x*x+3.0_d_t)
    ELSEIF(l == 5_li)THEN
      LegendreP=0.125_d_t*(63.0_d_t*x*x*x*x*x-70.0_d_t*x*x*x+&
            15.0_d_t*x)
    ELSEIF(l == 6_li)THEN
      LegendreP=0.0625_d_t*(231.0_d_t*x*x*x*x*x*x-315.0_d_t*x*x*x*x+&
            105.0_d_t*x*x-5.0_d_t)
    ELSEIF(l == 7_li)THEN
      LegendreP=0.0625_d_t*(429.0_d_t*x*x*x*x*x*x*x-&
            693.0_d_t*x*x*x*x*x+315.0_d_t*x*x*x-35.0_d_t*x)
    ELSE
      CALL stop_thor(23_li)
    END IF

  END FUNCTION LegendreP

END MODULE sph_harmonics_module
