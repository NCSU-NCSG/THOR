module vector_types
!***********************************************************************
! This is a general derived-type definition for vectors and their
! operations
!***********************************************************************

  use types

  implicit none

  type vector
    real(kind=d_t) :: x1, x2, x3 
  end type vector

  interface operator(+)
     module procedure plus
  end interface

  interface operator(-)
     module procedure minus
  end interface

  interface operator(*) 
     module procedure mult
  end interface

  interface operator(*) 
     module procedure matrix_mult
  end interface

  interface assignment(=)
     module procedure assign_r
  end interface

  interface operator(.dot.)
     module procedure vdot
  end interface

  interface operator(.cross.)
     module procedure vcross
  end interface

  interface abs
     module procedure vmag
  end interface

  interface operator(.proj.)
     module procedure vprojection
  end interface

contains
  
  pure function plus(v1,v2) 
    type(vector) :: plus
    type(vector), intent(in) :: v1,v2
    
    plus%x1 = v1%x1 + v2%x1
    plus%x2 = v1%x2 + v2%x2
    plus%x3 = v1%x3 + v2%x3

  end function plus

  pure function minus(v1, v2)
    type(vector) :: minus
    type(vector), intent(in) :: v1,v2
    
    minus%x1 = v1%x1 - v2%x1
    minus%x2 = v1%x2 - v2%x2
    minus%x3 = v1%x3 - v2%x3

  end function minus

  pure function mult(r, vec)
    type(vector) :: mult
    type(vector), intent(in) :: vec
    real(kind=d_t),intent(in) :: r
    
    mult%x1 = r * vec%x1
    mult%x2 = r * vec%x2
    mult%x3 = r * vec%x3
    
  end function mult

  pure function matrix_mult(m, v1)
    type(vector) :: matrix_mult
    type(vector), intent(in) :: v1
    real(kind=d_t), dimension(3,3), intent(in) :: m
    
    matrix_mult%x1 = m(1,1) * v1%x1 + m(1,2) * v1%x2 + m(1,3) * v1%x3
    matrix_mult%x2 = m(2,1) * v1%x1 + m(2,2) * v1%x2 + m(2,3) * v1%x3
    matrix_mult%x3 = m(3,1) * v1%x1 + m(3,2) * v1%x2 + m(3,3) * v1%x3
    
  end function matrix_mult

  pure function vdot(v1, v2)
    real(kind=d_t) :: vdot
    type(vector), intent(in) :: v1,v2
    
    vdot = v1%x1 * v2%x1 + v1%x2 * v2%x2 + v1%x3 * v2%x3

  end function vdot

  pure function vcross(v1, v2)
    type(vector) :: vcross
    type(vector), intent(in) :: v1,v2

    vcross%x1 = v1%x2 * v2%x3 - v1%x3 * v2%x2
    vcross%x2 = v1%x3 * v2%x1 - v1%x1 * v2%x3
    vcross%x3 = v1%x1 * v2%x2 - v1%x2 * v2%x1
    
  end function vcross

  pure function vmag_sqr(vec) 
    type(vector), intent(in) :: vec
    real(kind=d_t) :: vmag_sqr
    
    vmag_sqr = vdot(vec,vec)

  end function vmag_sqr

  pure subroutine assign_r(V, Ain)
    type(vector), intent(inout) :: v
    real, dimension(3), intent(in) :: Ain
    
    V%x1 = Ain(1)
    V%x2 = Ain(2)
    V%x3 = Ain(3)

  end subroutine assign_r

  pure function vmag(vec)
    type(vector), intent(in) :: vec
    real(kind=d_t) :: vmag
    
    vmag =  sqrt(vmag_sqr(vec))
  end function vmag

  pure function vprojection(from_v, onto_v)
    type(vector), intent(in) :: from_v, onto_v
    type(vector) :: vprojection
    
    vprojection = (vdot(from_v, onto_v)/vmag_sqr(onto_v)) * onto_v
  end function vprojection

  ! Project from_v onto the plane normal to normal_v
  ! length(normal_v) must be 1 for this to work properly

  pure function perp_projection(from_v, normal_v)
    type(vector), intent(in) :: from_v, normal_v
    type(vector) :: perp_projection
    
    perp_projection = from_v - vdot(from_v,normal_v)*normal_v
  end function perp_projection

  pure function normalize(v)
    type(vector),intent(in) :: v
    type(vector) :: normalize
    
    normalize = (1.0_d_t/vmag(v)) * v

  end function normalize

end module vector_types

