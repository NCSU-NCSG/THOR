!globals module
MODULE globals
  USE precisions
  IMPLICIT NONE

  REAL(kr8),ALLOCATABLE :: flux(:,:),volume(:),resp_func(:,:)

  INTEGER(ki4) :: num_cells,num_groups

  REAL(kr8) :: resp_value

  CHARACTER(64) :: response_type

  !input filename
  CHARACTER(64) :: response_inp

CONTAINS

END MODULE globals
