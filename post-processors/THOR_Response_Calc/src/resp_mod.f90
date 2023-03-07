!THOR is licensed under the MIT License.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!> @brief Response calculation.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
MODULE resp_mod
  USE globals
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: calcresp
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Gets command line arguments and read input file
  SUBROUTINE calcresp()
    INTEGER(ki4) :: i,g

    resp_value=0.0D0
    DO i=1,num_cells
      DO g=1,num_groups
        resp_value=resp_value+volume(i)*flux(i,g)*resp_func(i,g)
      ENDDO
    ENDDO
  ENDSUBROUTINE calcresp
END MODULE resp_mod
