MODULE read_inflow_module
  !***********************************************************************
  !
  ! Read source module contains all subroutines needed to read boundary
  ! source file
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

  SUBROUTINE read_finflow
    !*********************************************************************
    !
    ! Subroutine reads source in 'unique' ahot format
    ! (could be adapted for other formats, of course)
    !
    !*********************************************************************

    ! Declare temporary variables

    INTEGER(kind=li) :: alloc_stat, eg, q, octant, m, cell, f, face
    REAL(kind=d_t)   :: dmy

    ! Depending on page_iflw set eg_iflw

    IF(page_iflw .EQ. 1) THEN
      eg_iflw=1_li
    ELSE
      eg_iflw=egmax
    END IF

    ! allocate binflx: indices followig column major. Fastest running index
    ! is fixed boundary faces then octants then angles, then groups

    ALLOCATE(binflx(num_moments_f,fside_cells,8,nangle,eg_iflw),stat=alloc_stat)
    IF(alloc_stat /= 0) CALL stop_thor(.FALSE.,"*** Not enough memory ***")

    ! if page_iflw == 0 then read binflow in full, otherwise just open a file

    IF(page_iflw.EQ.0) THEN

      ! Open and read source file

      OPEN(unit=10,file=TRIM(finflow_filename),status='unknown',action='read')

      ! Read source strength and moments from file

      DO eg=1,egmax
        DO q=1,nangle
          DO octant=1,8
            DO f=1,fside_cells
              READ(10,*) face    ! this is the face number in the b_cells array
              IF( b_cells(face)%bc .NE. 2 ) THEN
                CALL stop_thor(.FALSE.,"A fixed inflow flux value is placed in a boundary face that is not declared fixed inflow.")
              END IF
              face = b_cells(face)%ptr   ! make sure that bc are stored in the same order as in fb_cells
              READ(10,*) (binflx(m,face,octant,q,eg),m=1,num_moments_f)
            END DO
          END DO
        END DO
      END DO

      ! Close inflow file

      CLOSE(10)

    ELSE

      OPEN(unit=97,file=TRIM(finflow_filename),status='unknown',action='read')
      DO eg=1,egmax
        DO q=1,nangle
          DO octant=1,8
            DO f=1,fside_cells
              READ(97,*) face    ! this is the face number in the b_cells array
              IF( b_cells(face)%bc .NE. 2 ) THEN
                CALL stop_thor(.FALSE.,"A fixed inflow flux value is placed in a boundary face that is not declared fixed inflow.")
              END IF
              READ(97,*) (dmy,m=1,num_moments_f)
            END DO
          END DO
        END DO
      END DO
      REWIND(unit=97)

    END IF

  END SUBROUTINE read_finflow

END MODULE read_inflow_module
