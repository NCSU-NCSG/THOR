!THOR is licensed under the MIT License.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!> @brief Angle types.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
MODULE angle_types
  !***********************************************************************
  ! This is hardwired to be discrete ordinates
  !***********************************************************************

  USE types
  USE vector_types

  IMPLICIT NONE

  ! Discrete Ordinate quadrature is retrived from precomputed database

  TYPE ordinate
    TYPE(vector) :: mu
    REAL(kind=d_t) :: wt
  END TYPE ordinate

END MODULE angle_types
