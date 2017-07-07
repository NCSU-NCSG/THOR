module angle_types
!*********************************************************************** 
! This is hardwired to be discrete ordinates
!***********************************************************************

  use types
  use vector_types

  implicit none

! Discrete Ordinate quadrature is retrived from precomputed database
  
  type ordinate
     type(vector) :: mu
     real(kind=d_t) :: wt
  end type ordinate

end module angle_types

