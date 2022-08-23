!-------------------------------------------------------------------------------
! THOR_Response_calc to calculate a response function from THOR outputs
!> @author Nicholas F. Herring
!-------------------------------------------------------------------------------
PROGRAM THOR_Response_calc
  USE globals
  USE infuncs
  USE outfuncs
  USE resp_mod
  IMPLICIT NONE

  !read in input
  CALL readinp()
  !calculate the response function
  CALL calcresp()
  !output response result
  CALL outputresp()
END PROGRAM THOR_Response_calc
