!THOR is licensed under the MIT License.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!> @brief Source types.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
MODULE source_types
  !***********************************************************************
  ! Cross-section derived type
  !***********************************************************************

  USE types

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: source_type

  INTEGER,PARAMETER :: name_size=64

  !source type
  TYPE :: source_type
    !source id
    INTEGER(kind=li) :: src_id
    !strength of the source by spatial moment, angular moment, and group
    REAL(kind=d_t), DIMENSION(:,:,:), ALLOCATABLE :: mom
  ENDTYPE source_type

END MODULE source_types
