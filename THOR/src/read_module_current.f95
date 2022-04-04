!***********************************************************************
!
! The module contains all legacy subroutines for old yaml_input versions. These
! subroutines should generally not be changed to garauntee backwards
! compatibility.
!
!***********************************************************************
MODULE read_module_current

  IMPLICIT NONE
  PRIVATE
  !
  ! List of public members
  PUBLIC :: inputfile_read

!> A list of valid card names for this block.
CHARACTER(100),DIMENSION(1),PARAMETER :: &
    cardnames=(/'type'/)

CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE inputfile_read(local_unit)
    !Command Line Arg Parsing
    INTEGER, INTENT(IN) :: local_unit
    STOP 'inputfile_read not yet complete'
  END SUBROUTINE inputfile_read

END MODULE read_module_current