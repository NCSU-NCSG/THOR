MODULE cross_section_types
  !***********************************************************************
  ! Cross-section derived type
  !***********************************************************************

  USE types

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: cross_section_mat
  PUBLIC :: cross_section
  PUBLIC :: xs_material_type

  INTEGER,PARAMETER :: name_size=64

  ! Define total and scattering cross-section

  TYPE cross_section_mat
    INTEGER(kind=li) :: mat
  END TYPE cross_section_mat

  TYPE cross_section
    REAL(kind=d_t) :: xs
  END TYPE cross_section

  !cross section material type
  TYPE :: xs_material_type
    !material name
    CHARACTER(name_size) :: mat_name
    !fission spectrum
    CHARACTER(name_size), DIMENSION(:), ALLOCATABLE :: chi
  ENDTYPE xs_material_type

END MODULE cross_section_types
