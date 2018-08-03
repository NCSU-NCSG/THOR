MODULE read_source_module
  !***********************************************************************
  !
  ! Read source module contains all subroutines needed to read source file
  !
  !***********************************************************************
  USE types
  USE parameter_types
  USE filename_types
  USE multindex_types
  USE global_variables
  USE termination_module

  IMPLICIT NONE

CONTAINS

  SUBROUTINE read_src
    !*********************************************************************
    !
    ! Subroutine reads source in 'unique' ahot format
    ! (could be adapted for other formats, of course)
    !
    !*********************************************************************

    ! Declare temporary variables

    INTEGER(kind=li) :: alloc_stat, eg, m, l

    ! Open and read source file

    OPEN(unit=10,file=TRIM(source_filename),status='unknown',action='read')

    ! Read source strength and moments from file

    READ(10,*) num_src_mat

    ALLOCATE(src_mat(num_src_mat),stat=alloc_stat)
    IF(alloc_stat /= 0) CALL stop_thor(2_li)

    ALLOCATE(src_str(0:num_src_mat-1,egmax),stat=alloc_stat)
    IF(alloc_stat /= 0) CALL stop_thor(2_li)

    ALLOCATE(src_m(num_moments_v,0:num_src_mat-1,egmax),stat=alloc_stat)
    IF(alloc_stat /= 0) CALL stop_thor(2_li)

    DO m=1, num_src_mat
      READ(10,*) src_mat(m)
      IF (src_mat(m) .ge. num_src_mat) CALL stop_thor(-1_li,&
                                       "source region ID is 0-indexed and must be < # source regions")
      DO eg=1, egmax
        READ(10,*) src_str(src_mat(m),eg)
        READ(10,*) (src_m(l,src_mat(m),eg), l = 1, num_moments_v)
      END DO
    END DO

    ! Close mesh file

    CLOSE(10)

  END SUBROUTINE read_src

END MODULE read_source_module
