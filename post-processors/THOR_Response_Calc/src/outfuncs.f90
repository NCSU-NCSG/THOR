!THOR is licensed under the MIT License.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!> @brief Output module.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
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

    OPEN(UNIT=30,FILE=TRIM(ADJUSTL(response_inp))//'_response.out',ACTION='WRITE',STATUS='REPLACE')
    WRITE(30,'(A,ES24.16)')'Computed Response: ',resp_value
    CLOSE(30)
  ENDSUBROUTINE outputresp
END MODULE outfuncs
