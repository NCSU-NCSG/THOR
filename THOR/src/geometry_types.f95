MODULE geometry_types
  !***********************************************************************
  ! Derived types for geometry (algorithm was built around these!)
  !***********************************************************************

  USE types
  USE vector_types

  IMPLICIT NONE

  TYPE vertex
    TYPE(vector) :: v
  END TYPE vertex

  TYPE cell
    INTEGER(kind=li), DIMENSION(0:3) :: R
    INTEGER(kind=li), DIMENSION(0:3) :: face
    INTEGER(kind=li) :: reg
    INTEGER(kind=li) :: src
    REAL(kind=d_t) :: volume
  END TYPE cell

  TYPE boundary_cell
    INTEGER(kind=li) :: cell, face, bc, ptr
  END TYPE boundary_cell

  TYPE boundary_condition
    INTEGER(kind=li), DIMENSION(0:3) :: bc
  END TYPE boundary_condition

  TYPE list
    INTEGER(kind=li) :: cell, face
  END TYPE list

  TYPE, EXTENDS(list) :: SDD_list_entry
    INTEGER(kind=li) :: proc
  END TYPE SDD_list_entry

  TYPE linked_list
    INTEGER(kind=li) :: cell_id
    TYPE(linked_list), POINTER :: next
    TYPE(linked_list), POINTER :: prev
  END TYPE linked_list

  TYPE upstream_face_le
    INTEGER(kind=li) :: cell,face
    TYPE(upstream_face_le),POINTER :: next
  END TYPE upstream_face_le

  TYPE cycle_breaker
    INTEGER(kind=li),ALLOCATABLE :: faces(:)
    INTEGER(kind=li),ALLOCATABLE :: cells(:)
    REAL(kind=d_t)  ,ALLOCATABLE :: face_fluxes(:,:)
  END TYPE cycle_breaker

END MODULE geometry_types
