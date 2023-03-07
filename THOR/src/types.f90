!THOR is licensed under the MIT License.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!> @brief General types.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
MODULE types
  !***********************************************************************
  ! This derived type sets the precision of other declared variables
  !***********************************************************************

  IMPLICIT NONE

  INTEGER, PARAMETER :: li = SELECTED_INT_KIND(8)
  INTEGER, PARAMETER :: d_t = SELECTED_REAL_KIND(15,307)

END MODULE types
