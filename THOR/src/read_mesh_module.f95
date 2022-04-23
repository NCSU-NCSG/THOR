MODULE read_mesh_module
  !***********************************************************************
  !
  ! Read mesh module contains all subroutines needed for reading
  ! tetrahedral mesh file
  !
  !***********************************************************************
  USE mpi
  USE types
  USE parameter_types
  USE filename_types
  USE vector_types
  USE geometry_types
  USE global_variables
  USE error_module

  IMPLICIT NONE

CONTAINS

  SUBROUTINE read_tetmesh
    !*********************************************************************
    !
    ! Subroutine reads tetrahedral mesh from THOR mesh file
    !
    !*********************************************************************

    ! Define temporary variables

    INTEGER(kind=li) :: alloc_stat, i, j, k, dummy
    INTEGER(kind=li) :: kv,kr,kf
    INTEGER(kind=li), DIMENSION(:), ALLOCATABLE:: cell_temp, face, &
          side_cells_tmp
    REAL(kind=d_t) :: rdummy
    INTEGER ::rank,mpi_err, localunit,temp_minreg,temp_maxreg
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)
    localunit = rank+100

    ! Set temp_minreg,temp_maxreg

    temp_minreg= 100000_li
    temp_maxreg=-1_li

    ! Open and read mesh file

    OPEN(unit=localunit,file=TRIM(mesh_filename),status='old',action='read')

    ! Read general parameters

    READ(localunit,*) num_vert
    READ(localunit,*) num_cells
    READ(localunit,*) num_cell_blk
    READ(localunit,*) num_side_sets

    ! Allocate vertices type dimension

    ALLOCATE(vertices(num_vert),stat=alloc_stat)
    IF(alloc_stat /= 0) CALL raise_fatal_error("*** Not enough memory ***")

    ! Read vertices

    DO i=1, num_vert
      READ(localunit,*) dummy, vertices(i)%v
    END DO

    ! Allocate cell type dimension

    ALLOCATE(cells(num_cells),stat=alloc_stat)
    IF(alloc_stat /= 0) CALL raise_fatal_error("*** Not enough memory ***")

    ! Read cell block property list (region and source composition)
    DO i=1, num_cells
      READ(localunit,*) dummy, cells(i)%reg, cells(i)%src
      IF(cells(i)%reg>temp_maxreg) temp_maxreg=cells(i)%reg
      IF(cells(i)%reg<temp_minreg) temp_minreg=cells(i)%reg
    END DO

    IF(minreg .EQ. 100000_li .AND. maxreg .EQ. -1_li)THEN
      !region map was not found in input file, so assign here
      maxreg=temp_maxreg
      minreg=temp_minreg
      ALLOCATE(reg2mat(minreg:maxreg))
      DO i=minreg,maxreg
        reg2mat(i)=i
      ENDDO
    ELSE
      !region map was found in input file, so check that they match
      IF(temp_maxreg .GT. maxreg .OR. temp_minreg .LT. minreg)THEN
        WRITE(*,*)'Error: Min index in input file',minreg,'min index in mesh',temp_minreg
        WRITE(*,*)'Error: Max index in input file',maxreg,'max index in mesh',temp_maxreg
        CALL raise_fatal_error('Region map in input file and regions in mesh file do not match! Be &
          & sure to index from 1 if using legacy regionmap!')
      ENDIF
    ENDIF

    ! Read vertices to cell mapping

    DO i=1, num_cells
      READ(localunit,*) dummy, cells(i)%R(0), cells(i)%R(1), &
            cells(i)%R(2), cells(i)%R(3)
    END DO

    ! Go through side_sets once to count side_cells

    side_cells =0
    vside_cells=0
    rside_cells=0
    fside_cells=0

    ALLOCATE(side_cells_tmp(num_side_sets),stat=alloc_stat)
    IF(alloc_stat /= 0) CALL raise_fatal_error("*** Not enough memory ***")

    DO i=1, num_side_sets

      READ(localunit,*) side_cells_tmp(i)

      DO j=1, side_cells_tmp(i)
        READ(localunit,*) dummy,dummy,k
        IF      (k.EQ.0) THEN
          vside_cells=vside_cells+1
        ELSE IF (k.EQ.1) THEN
          rside_cells=rside_cells+1
        ELSE IF (k.EQ.2) THEN
          fside_cells=fside_cells+1
        END IF
      END DO

      side_cells=side_cells+side_cells_tmp(i)

    END DO

    ! If there are fixed inflow boundary faces then finflow flag must be ==1

    IF(fside_cells>0 .AND. finflow .EQ. 0) THEN
      CALL raise_fatal_error("Some boundary faces feature a fixed inflow flux but the  &
        & inflow flag is set to false. Execution terminates.")
    END IF

    ! Close file, reopen and read everything up to side_sets into dummies

    CLOSE(unit=localunit)

    OPEN(unit=localunit,file=TRIM(mesh_filename),status='old',action='read')

    READ(localunit,*) dummy
    READ(localunit,*) dummy
    READ(localunit,*) dummy
    READ(localunit,*) dummy

    DO i=1, num_vert
      READ(localunit,*) dummy, rdummy,rdummy,rdummy
    END DO

    DO i=1, num_cells
      READ(localunit,*) dummy, dummy, dummy
    END DO

    DO i=1, num_cells
      READ(localunit,*) dummy, dummy,dummy,dummy,dummy
    END DO

    ! Allocate b_cells and read from side sets

    ALLOCATE(b_cells(side_cells),stat=alloc_stat)
    IF(alloc_stat /= 0) CALL raise_fatal_error("*** Not enough memory ***")

    IF(vside_cells .GT. 0) THEN
      ALLOCATE(vb_cells(vside_cells),stat=alloc_stat)
      IF(alloc_stat /= 0) CALL raise_fatal_error("*** Not enough memory ***")
    END IF

    IF(rside_cells .GT. 0) THEN
      ALLOCATE(rb_cells(rside_cells),stat=alloc_stat)
      IF(alloc_stat /= 0) CALL raise_fatal_error("*** Not enough memory ***")
    END IF

    IF(fside_cells .GT. 0) THEN
      ALLOCATE(fb_cells(fside_cells),stat=alloc_stat)
      IF(alloc_stat /= 0) CALL raise_fatal_error("*** Not enough memory ***")
    END IF

    k=1 ; kv=1 ; kr=1 ; kf=1

    DO i=1, num_side_sets
      READ(localunit,*) dummy
      DO j=1, side_cells_tmp(i)
        READ(localunit,*) b_cells(k)%cell,b_cells(k)%face,b_cells(k)%bc
        ! work on boundary data
        IF      (b_cells(k)%bc.EQ.0) THEN
          vb_cells(kv)%cell=b_cells(k)%cell
          vb_cells(kv)%face=b_cells(k)%face
          vb_cells(kv)%bc  =0
          vb_cells(kv)%ptr =k                 ! ptr points to the position in the b_cells array
          b_cells(k)%ptr   =kv                ! ptr points from kth position in b_cells to kv position in vb_cells array
          kv=kv+1
        ELSE IF (b_cells(k)%bc.EQ.1) THEN
          rb_cells(kr)%cell=b_cells(k)%cell
          rb_cells(kr)%face=b_cells(k)%face
          rb_cells(kr)%bc  =1
          rb_cells(kr)%ptr =k
          b_cells(k)%ptr   =kr
          kr=kr+1
        ELSE IF (b_cells(k)%bc.EQ.2) THEN
          fb_cells(kf)%cell=b_cells(k)%cell
          fb_cells(kf)%face=b_cells(k)%face
          fb_cells(kf)%bc  =2
          fb_cells(kf)%ptr =k
          b_cells(k)%ptr   =kf
          kf=kf+1
        END IF
        k=k+1
      END DO
    END DO

    DEALLOCATE(side_cells_tmp)

    ! Read adjacency list

    READ(localunit,*) adjacent_cells
    ALLOCATE(adjacency_list(adjacent_cells,0:3),&
          cell_temp(adjacent_cells),face(adjacent_cells),&
          stat=alloc_stat)
    IF(alloc_stat /= 0) CALL raise_fatal_error("*** Not enough memory ***")

    DO i=1, adjacent_cells
      READ(localunit,*) cell_temp(i), face(i),&
            adjacency_list(cell_temp(i),face(i))%cell, &
            adjacency_list(cell_temp(i),face(i))%face
    END DO

    DEALLOCATE(cell_temp,face)

    ! Close mesh file

    CLOSE(localunit)

  END SUBROUTINE read_tetmesh

END MODULE read_mesh_module
