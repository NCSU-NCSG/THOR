module types
!***********************************************************************
! This derived type sets the precision of other declared variables 
!***********************************************************************

  implicit none

  integer, parameter :: li = selected_int_kind(8)
  integer, parameter :: d_t = selected_real_kind(15,307)

end module types

