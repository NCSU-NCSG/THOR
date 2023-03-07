!THOR is licensed under the MIT License.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!> @brief Vector types.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
MODULE vector_types
  !***********************************************************************
  ! This is a general derived-type definition for vectors and their
  ! operations
  !***********************************************************************

  USE types

  IMPLICIT NONE

  TYPE vector
    REAL(kind=d_t) :: x1, x2, x3
  END TYPE vector

  INTERFACE OPERATOR(+)
    MODULE PROCEDURE plus
  END INTERFACE

  INTERFACE OPERATOR(-)
    MODULE PROCEDURE minus
  END INTERFACE

  INTERFACE OPERATOR(*)
    MODULE PROCEDURE mult
  END INTERFACE

  INTERFACE OPERATOR(*)
    MODULE PROCEDURE matrix_mult
  END INTERFACE

  INTERFACE ASSIGNMENT(=)
    MODULE PROCEDURE assign_r
  END INTERFACE

  INTERFACE OPERATOR(.dot.)
    MODULE PROCEDURE vdot
  END INTERFACE

  INTERFACE OPERATOR(.cross.)
    MODULE PROCEDURE vcross
  END INTERFACE

  INTERFACE abs
    MODULE PROCEDURE vmag
  END INTERFACE

  INTERFACE OPERATOR(.proj.)
    MODULE PROCEDURE vprojection
  END INTERFACE

CONTAINS

  PURE FUNCTION plus(v1,v2)
    TYPE(vector) :: plus
    TYPE(vector), INTENT(in) :: v1,v2

    plus%x1 = v1%x1 + v2%x1
    plus%x2 = v1%x2 + v2%x2
    plus%x3 = v1%x3 + v2%x3

  END FUNCTION plus

  PURE FUNCTION minus(v1, v2)
    TYPE(vector) :: minus
    TYPE(vector), INTENT(in) :: v1,v2

    minus%x1 = v1%x1 - v2%x1
    minus%x2 = v1%x2 - v2%x2
    minus%x3 = v1%x3 - v2%x3

  END FUNCTION minus

  PURE FUNCTION mult(r, vec)
    TYPE(vector) :: mult
    TYPE(vector), INTENT(in) :: vec
    REAL(kind=d_t),INTENT(in) :: r

    mult%x1 = r * vec%x1
    mult%x2 = r * vec%x2
    mult%x3 = r * vec%x3

  END FUNCTION mult

  PURE FUNCTION matrix_mult(m, v1)
    TYPE(vector) :: matrix_mult
    TYPE(vector), INTENT(in) :: v1
    REAL(kind=d_t), DIMENSION(3,3), INTENT(in) :: m

    matrix_mult%x1 = m(1,1) * v1%x1 + m(1,2) * v1%x2 + m(1,3) * v1%x3
    matrix_mult%x2 = m(2,1) * v1%x1 + m(2,2) * v1%x2 + m(2,3) * v1%x3
    matrix_mult%x3 = m(3,1) * v1%x1 + m(3,2) * v1%x2 + m(3,3) * v1%x3

  END FUNCTION matrix_mult

  PURE FUNCTION vdot(v1, v2)
    REAL(kind=d_t) :: vdot
    TYPE(vector), INTENT(in) :: v1,v2

    vdot = v1%x1 * v2%x1 + v1%x2 * v2%x2 + v1%x3 * v2%x3

  END FUNCTION vdot

  PURE FUNCTION vcross(v1, v2)
    TYPE(vector) :: vcross
    TYPE(vector), INTENT(in) :: v1,v2

    vcross%x1 = v1%x2 * v2%x3 - v1%x3 * v2%x2
    vcross%x2 = v1%x3 * v2%x1 - v1%x1 * v2%x3
    vcross%x3 = v1%x1 * v2%x2 - v1%x2 * v2%x1

  END FUNCTION vcross

  PURE FUNCTION vmag_sqr(vec)
    TYPE(vector), INTENT(in) :: vec
    REAL(kind=d_t) :: vmag_sqr

    vmag_sqr = vdot(vec,vec)

  END FUNCTION vmag_sqr

  PURE SUBROUTINE assign_r(V, Ain)
    TYPE(vector), INTENT(inout) :: v
    REAL, DIMENSION(3), INTENT(in) :: Ain

    V%x1 = Ain(1)
    V%x2 = Ain(2)
    V%x3 = Ain(3)

  END SUBROUTINE assign_r

  PURE FUNCTION vmag(vec)
    TYPE(vector), INTENT(in) :: vec
    REAL(kind=d_t) :: vmag

    vmag =  SQRT(vmag_sqr(vec))
  END FUNCTION vmag

  PURE FUNCTION vprojection(from_v, onto_v)
    TYPE(vector), INTENT(in) :: from_v, onto_v
    TYPE(vector) :: vprojection

    vprojection = (vdot(from_v, onto_v)/vmag_sqr(onto_v)) * onto_v
  END FUNCTION vprojection

  ! Project from_v onto the plane normal to normal_v
  ! length(normal_v) must be 1 for this to work properly

  PURE FUNCTION perp_projection(from_v, normal_v)
    TYPE(vector), INTENT(in) :: from_v, normal_v
    TYPE(vector) :: perp_projection

    perp_projection = from_v - vdot(from_v,normal_v)*normal_v
  END FUNCTION perp_projection

  PURE FUNCTION normalize(v)
    TYPE(vector),INTENT(in) :: v
    TYPE(vector) :: normalize

    normalize = (1.0_d_t/vmag(v)) * v

  END FUNCTION normalize

END MODULE vector_types
