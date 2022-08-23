MODULE cross_section_types
  !***********************************************************************
  ! Cross-section derived type
  !***********************************************************************

  USE types

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: xs_material_type

  INTEGER,PARAMETER :: name_size=64

  !cross section material type
  TYPE :: xs_material_type
    !material name
    CHARACTER(name_size) :: mat_name
    !material ID
    INTEGER(kind=li) :: mat_id
    !indicator if material is fissile
    LOGICAL :: fissile=.TRUE.
    !fission spectrum
    REAL(kind=d_t), DIMENSION(:), ALLOCATABLE :: chi
    !fission cross section SigmaF
    REAL(kind=d_t), DIMENSION(:), ALLOCATABLE :: sigma_f
    !fission production nu
    REAL(kind=d_t), DIMENSION(:), ALLOCATABLE :: nu
    !fission production cross section nuSigmaF
    REAL(kind=d_t), DIMENSION(:), ALLOCATABLE :: nusig_f
    !total or transport cross section sigma_t
    REAL(kind=d_t), DIMENSION(:), ALLOCATABLE :: sigma_t
    !total scattering cross section for reaction rates
    REAL(kind=d_t), DIMENSION(:), ALLOCATABLE :: tsigs
    !scattering matrix
    REAL(kind=d_t), DIMENSION(:,:,:), ALLOCATABLE :: sigma_scat
  ENDTYPE xs_material_type

END MODULE cross_section_types
