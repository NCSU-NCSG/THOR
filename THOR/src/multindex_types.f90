module multindex_types
!***********************************************************************
! These are defined for double precision
!***********************************************************************

  use types
  use vector_types

  implicit none

! Multi indices are pre-calculated for increased efficiency

  type indices_v
     integer(kind=li) :: i1, i2, i3
  end type indices_v

  type indices_f
     integer(kind=li) :: i1, i2
  end type indices_f

end module multindex_types

