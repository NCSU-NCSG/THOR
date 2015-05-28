module cross_section_types
!*********************************************************************** 
! Cross-section derived type
!***********************************************************************

  use types

  implicit none

! Define total and scattering cross-section

  type cross_section_mat
     integer(kind=li) :: mat
  end type cross_section_mat
  
  type cross_section
     real(kind=d_t) :: xs
  end type cross_section

end module cross_section_types

