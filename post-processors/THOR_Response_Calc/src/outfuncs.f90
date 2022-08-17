!output functions
MODULE outfuncs
  USE globals
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: outputresp
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !outputs response
  SUBROUTINE outputresp()
    WRITE(*,'(A,ES16.8)')'Response: ',resp_value
  ENDSUBROUTINE outputresp
END MODULE outfuncs
