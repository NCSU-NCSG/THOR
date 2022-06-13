!globals module
MODULE globals
  IMPLICIT NONE

  !number of materials
  INTEGER :: nummats
  !number of energy groups
  INTEGER :: numgroups
  !level of anisotropy
  INTEGER :: levelanis
  !xs input filename
  CHARACTER(64) :: xsin
  !xs output filename
  CHARACTER(80) :: xsout
  !xs output format
  CHARACTER(64) :: outformat
  !energy group structure
  REAL(8),ALLOCATABLE :: eg_struc(:)
  !cross section variables
  REAL(8),ALLOCATABLE :: chi(:,:),sigmaf(:,:),nuf(:,:),sigmat(:,:),sigmas(:,:,:,:)
CONTAINS

END MODULE globals
