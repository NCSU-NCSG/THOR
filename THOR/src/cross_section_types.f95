MODULE cross_section_types
  !***********************************************************************
  ! Cross-section derived type
  !***********************************************************************

  USE types

  IMPLICIT NONE

  ! Define total and scattering cross-section

  TYPE cross_section_mat
    INTEGER(kind=li) :: mat
  END TYPE cross_section_mat

  TYPE cross_section
    REAL(kind=d_t) :: xs
  END TYPE cross_section

END MODULE cross_section_types
