MODULE multindex_types
  !***********************************************************************
  ! These are defined for double precision
  !***********************************************************************

  USE types
  USE vector_types

  IMPLICIT NONE

  ! Multi indices are pre-calculated for increased efficiency

  TYPE indices_v
    INTEGER(kind=li) :: i1, i2, i3
  END TYPE indices_v

  TYPE indices_f
    INTEGER(kind=li) :: i1, i2
  END TYPE indices_f

END MODULE multindex_types
