!-------------------------------------------------------------------------------
! THOR XS converter to convert different XS formats back and forth
!> @author Nicholas F. Herring
!-------------------------------------------------------------------------------
PROGRAM thor_xs_converter
  USE globals
  USE infuncs
  USE outfuncs
  IMPLICIT NONE

  !read in command line arguments
  CALL readcl()
  WRITE(*,'(A)')'Reading in xs data'
  !read in xs file
  CALL readxs()
  WRITE(*,'(A)')'Outputting xs data'
  !output xs file
  CALL outputxs()

  WRITE(*,'(A)')'XS Converter completed! No detected errors in conversion.'
END PROGRAM thor_xs_converter
