module geometry_types
!***********************************************************************
! Derived types for geometry (algorithm was built around these!)
!***********************************************************************

  use types
  use vector_types

  implicit none

  type vertex
     type(vector) :: v
  end type vertex

  type cell
     integer(kind=li), dimension(0:3) :: R
     integer(kind=li), dimension(0:3) :: face
     integer(kind=li) :: reg
     integer(kind=li) :: src
     real(kind=d_t) :: volume
  end type cell

  type boundary_cell
     integer(kind=li) :: cell, face, bc, ptr
  end type boundary_cell

  type boundary_condition
     integer(kind=li), dimension(0:3) :: bc
  end type boundary_condition

  type list
     integer(kind=li) :: cell, face
  end type list

  type linked_list
     integer(kind=li) :: cell_id
     type(linked_list), pointer :: next
     type(linked_list), pointer :: prev
  end type linked_list

  type upstream_face_le
     integer(kind=li) :: cell,face
     type(upstream_face_le),pointer :: next
  end type upstream_face_le

  type cycle_breaker
     integer(kind=li),allocatable :: faces(:)
     integer(kind=li),allocatable :: cells(:)
     real(kind=d_t)  ,allocatable :: face_fluxes(:,:)
  end type

end module geometry_types
