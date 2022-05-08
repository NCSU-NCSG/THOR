!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! THOR MESH GENERATION UTILITY
!   Globals Module:
!
!>    This module contains variables and functions used to store and manipulate
!!    data common to both the input mesh and output mesh.
!> @author Nicholas F. Herring
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
MODULE globals
  IMPLICIT NONE

  !input mesh filename
  CHARACTER(200) :: mesh_infile

  !output mesh filename
  CHARACTER(200) :: mesh_outfile

  !number of vertices
  INTEGER :: num_verts

  !number of tets
  INTEGER :: num_tets

  !number of boundary condition faces
  INTEGER :: num_bcf

  !vertex data
  REAL(8), ALLOCATABLE :: vertex(:,:)

  !element data
  INTEGER, ALLOCATABLE :: element(:,:)

  !element region tags
  INTEGER, ALLOCATABLE :: el_tag(:)

  !adjacency list
  INTEGER, ALLOCATABLE :: adj_list(:,:)

  !boundary conditions data
  INTEGER, ALLOCATABLE :: bc_data(:,:)

  !boundary conditions on each side ordered -x, +x, -y, +y, -z, +z
  INTEGER :: side_bc(6)=0

  !side flatnesses, ordered same as BCs.
  LOGICAL :: side_flat(6)=.FALSE.

  !pi
  REAL(8),PARAMETER :: PI=4.D0*DATAN(1.D0)
CONTAINS

END MODULE globals
