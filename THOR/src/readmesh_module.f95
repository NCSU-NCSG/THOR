MODULE readmesh_module
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
  USE termination_module
  USE SDD_global_variables

  IMPLICIT NONE

CONTAINS

  SUBROUTINE read_tetmesh
    !*********************************************************************
    !
    ! Subroutine reads tetrahedral mesh from CUBIT Exodus II file
    !
    !*********************************************************************

    ! Define temporary variables

    INTEGER(kind=li) :: alloc_stat, i, j, k, dummy,max_cell=0
    INTEGER(kind=li) :: kv,kr,kf,ks
    INTEGER(kind=li), DIMENSION(:), ALLOCATABLE:: cell_temp, face, &
          side_cells_tmp
    REAL(kind=d_t) :: rdummy
    CHARACTER::comment_check
    CHARACTER(50):: temp_str
    INTEGER ::rank=0,mpi_err, localunit, neigh_cell, neigh_proc, neigh_face

    !! ADD REMOVED - OCT 2019
    !CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)
    rank = 0
    localunit = rank+100
    !! ADD REMOVED - OCT 2019

    dummy = INDEX(mesh_filename,".", BACK=.TRUE.)
    READ(mesh_filename( (dummy+1) :),*) rank
    rank = rank - 1

    ! Set minreg,maxreg

    minreg= 100000_li
    maxreg=-1_li

    ! Open and read mesh file

    OPEN(unit=localunit,file=TRIM(mesh_filename),status='old',action='read')


    ! Read general parameters

    READ(localunit,*) num_vert, tot_num_vert_G
    READ(localunit,*) num_cells, tot_num_cells_G
    READ(localunit,*) num_cell_blk
    READ(localunit,*) num_side_sets

    !Allocate global to local / local to global

    IF ((ITMM.EQ.2).OR.(ITMM.EQ.3)) THEN
        ALLOCATE(SDD_cells_g2l_G(tot_num_cells_G))
        ALLOCATE(SDD_cells_l2g_G(num_cells))
        ALLOCATE(SDD_vert_g2l_G(tot_num_vert_G))
        ALLOCATE(SDD_vert_l2g_G(num_vert))
        SDD_cells_l2g_G = 0
        SDD_cells_g2l_G = -1
        SDD_vert_l2g_G = 0
        SDD_vert_g2l_G = -1
    END IF


    !Check for header line
    READ(localunit,*) comment_check
    IF (comment_check.NE.'!') BACKSPACE(localunit)

    ! Allocate vertices type dimension

    IF (ITMM.NE.1 .AND. ITMM.NE.0) THEN
        ALLOCATE(vertices(num_vert),stat=alloc_stat)
        IF(alloc_stat /= 0) CALL stop_thor(2_li)
    END IF

    ! Read vertices

    DO i=1, num_vert
      READ(localunit,*) SDD_vert_l2g_G(i), vertices(i)%v
      SDD_vert_g2l_G(SDD_vert_l2g_G(i))=i
    END DO

    ! Allocate cell type dimension

    IF (ITMM.NE.1 .AND. ITMM.NE.0) THEN
        ALLOCATE(cells(num_cells),stat=alloc_stat)
        IF(alloc_stat /= 0) CALL stop_thor(2_li)
    END IF

    !Check for header line
    READ(localunit,*) comment_check
    IF (comment_check.NE.'!') BACKSPACE(localunit)

    ! Read cell block property list (region and source composition)
    DO i=1, num_cells
      READ(localunit,*) SDD_cells_l2g_G(i), cells(i)%reg, cells(i)%src
      SDD_cells_g2l_G(SDD_cells_l2g_G(i))=i
      IF (SDD_cells_l2g_G(i).GT.max_cell) max_cell=SDD_cells_l2g_G(i)
      IF(cells(i)%reg>maxreg) maxreg=cells(i)%reg
      IF(cells(i)%reg<minreg) minreg=cells(i)%reg
    END DO

    ! Print *, num_cells, tot_num_cells_G
    ! PRINT *, "L2G: ",SDD_cells_l2g_G
    ! PRINT *, "G2L: ",SDD_cells_g2l_G

    !Check for header line
    READ(localunit,*) comment_check
    IF (comment_check.NE.'!') BACKSPACE(localunit)

    ! Read vertices to cell mapping (Element List)

    DO i=1, num_cells
      READ(localunit,*) dummy, cells(i)%R(0), cells(i)%R(1), &
            cells(i)%R(2), cells(i)%R(3)
            DO j=0,3
              cells(i)%R(j)=SDD_vert_g2l_G(cells(i)%R(j))
            END DO
    END DO

    !Check for header line
    READ(localunit,*) comment_check
    IF (comment_check.NE.'!') BACKSPACE(localunit)

    ! Go through side_sets once to count side_cells

    side_cells =0
    vside_cells=0
    rside_cells=0
    fside_cells=0
    SDD_side_cells = 0


    ALLOCATE(side_cells_tmp(num_side_sets),stat=alloc_stat)
    IF(alloc_stat /= 0) CALL stop_thor(2_li)

    IF ((ITMM.EQ.2).OR.(ITMM.EQ.3))N_side_SDbound=0
    DO i=1, num_side_sets

      READ(localunit,*) side_cells_tmp(i)

	  IF ((ITMM.EQ.2).OR.(ITMM.EQ.3)) N_side_SDbound=N_side_SDbound+side_cells_tmp(i)

      DO j=1, side_cells_tmp(i)
        READ(localunit,*) dummy,dummy,k
        IF      ((k.EQ.0).OR.(ITMM.EQ.2).OR.(ITMM.EQ.3)) THEN
          vside_cells=vside_cells+1
        ELSE IF (k.EQ.1) THEN
          rside_cells=rside_cells+1
        ELSE IF (k.EQ.2) THEN
          fside_cells=fside_cells+1
        ELSE IF (k .EQ. -1) THEN
          SDD_side_cells = SDD_side_cells + 1
          !vside_cells=vside_cells+1
        END IF
      END DO

      side_cells=side_cells+side_cells_tmp(i)

    END DO

    ! If there are fixed inflow boundary faces then finflow flag must be ==1

    IF(fside_cells>0 .AND. finflow .EQ. 0) THEN
      CALL stop_thor(12_li)
    END IF

    ! Close file, reopen and read everything up to side_sets into dummies

    CLOSE(unit=localunit)

    OPEN(unit=localunit,file=TRIM(mesh_filename),status='old',action='read')

    READ(localunit,*) dummy
    READ(localunit,*) dummy
    READ(localunit,*) dummy
    READ(localunit,*) dummy

    !Check for header line
    READ(localunit,*) comment_check
    IF (comment_check.NE.'!') BACKSPACE(localunit)

    DO i=1, num_vert
      READ(localunit,*) dummy, rdummy,rdummy,rdummy
    END DO

    !Check for header line
    READ(localunit,*) comment_check
    IF (comment_check.NE.'!') BACKSPACE(localunit)

    DO i=1, num_cells
      READ(localunit,*) dummy, dummy, dummy
    END DO

    !Check for header line
    READ(localunit,*) comment_check
    IF (comment_check.NE.'!') BACKSPACE(localunit)

    DO i=1, num_cells
      READ(localunit,*) dummy, dummy,dummy,dummy,dummy
    END DO

    ! Allocate b_cells and read from side sets
    !First, it we are redoing these BC after GFIC,
    !we need to deallocate, so they can be reallocated
    IF ((ITMM.EQ.1).OR.(ITMM.EQ.0)) DEALLOCATE(b_cells,vb_cells)

    ALLOCATE(b_cells(side_cells),stat=alloc_stat)
    IF((alloc_stat /= 0).AND.(ITMM.NE.1) .AND. ITMM.NE.0) CALL stop_thor(2_li)

    IF(vside_cells .GT. 0) THEN
      !ALLOCATE(vb_cells(vside_cells+SDD_side_cells),stat=alloc_stat)
      ALLOCATE(vb_cells(vside_cells),stat=alloc_stat)
      IF((alloc_stat /= 0).AND.(ITMM.NE.1) .AND. ITMM.NE.0) CALL stop_thor(2_li)
    END IF

    IF(rside_cells .GT. 0) THEN
      ALLOCATE(rb_cells(rside_cells),stat=alloc_stat)
      IF((alloc_stat /= 0).AND.(ITMM.NE.1) .AND. ITMM.NE.0) CALL stop_thor(2_li)
    END IF

    IF(fside_cells .GT. 0) THEN
      ALLOCATE(fb_cells(fside_cells),stat=alloc_stat)
      IF((alloc_stat /= 0).AND.(ITMM.NE.1) .AND. ITMM.NE.0) CALL stop_thor(2_li)
    END IF

    IF(SDD_side_cells .GT. 0) THEN
      ALLOCATE(SDDb_cells(SDD_side_cells),stat=alloc_stat)
      IF((alloc_stat /= 0).AND.(ITMM.NE.1) .AND. ITMM.NE.0) CALL stop_thor(2_li)
    END IF

    k=1 ; kv=1 ; kr=1 ; kf=1; ks=1

    !Check for header line
    READ(localunit,*) comment_check
    IF (comment_check.NE.'!') BACKSPACE(localunit)

    DO i=1, num_side_sets
      READ(localunit,*) dummy
      DO j=1, side_cells_tmp(i)
        READ(localunit,*) b_cells(k)%cell,b_cells(k)%face,b_cells(k)%bc
        ! Convert global cell numbers to local for BC enumeration
        b_cells(k)%cell = SDD_cells_g2l_G(b_cells(k)%cell)
        ! work on boundary data
        IF      ((b_cells(k)%bc.EQ.0).OR.(ITMM.EQ.2).OR.(ITMM.EQ.3)) THEN !If this is an ITMM construction, we impose vacuum BCs
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
        ELSE IF (b_cells(k)%bc.EQ.-1) THEN !SDD Side Cells will go here
          !print *, b_cells(k)%cell, b_cells(k)%face
          SDDb_cells(ks)%cell=b_cells(k)%cell
          SDDb_cells(ks)%face=b_cells(k)%face
          SDDb_cells(ks)%bc  =-1
          SDDb_cells(ks)%ptr =k                 ! ptr points to the position in the b_cells array
          b_cells(k)%ptr   =ks                ! ptr points from kth position in b_cells to kv position in vb_cells array
          ks=ks+1
          !vb_cells(kv)%cell=b_cells(k)%cell
          !vb_cells(kv)%face=b_cells(k)%face
          !vb_cells(kv)%bc  =0
          !vb_cells(kv)%ptr =k                 ! ptr points to the position in the b_cells array
          !b_cells(k)%ptr   =kv                ! ptr points from kth position in b_cells to kv position in vb_cells array
          !kv=kv+1
        END IF
        k=k+1
      END DO
    END DO

    DEALLOCATE(side_cells_tmp)

    IF (ITMM.EQ.1 .OR. ITMM.EQ.0) CLOSE(UNIT=localunit)
    IF (ITMM.EQ.1 .OR. ITMM.EQ.0) RETURN

    !Check for header line
    READ(localunit,*) comment_check
    IF (comment_check.NE.'!') BACKSPACE(localunit)

    ! Read adjacency list

    READ(localunit,*) adjacent_cells
    ALLOCATE(adjacency_list(adjacent_cells,0:3),&
          cell_temp(adjacent_cells),face(adjacent_cells),&
          stat=alloc_stat)
    IF(alloc_stat /= 0) CALL stop_thor(2_li)
    ALLOCATE(SDD_adjacency_list(max_cell,0:3), stat=alloc_stat)
    IF(alloc_stat /= 0) CALL stop_thor(2_li)
    DO i=1, adjacent_cells
      READ(localunit, *) cell_temp(i), face(i), neigh_cell, neigh_face, neigh_proc
      SDD_adjacency_list(cell_temp(i),face(i))%cell = neigh_cell
      SDD_adjacency_list(cell_temp(i),face(i))%face = neigh_face
      SDD_adjacency_list(cell_temp(i),face(i))%proc = neigh_proc

      IF (neigh_proc .EQ. (rank + 1)) THEN
        adjacency_list(SDD_cells_g2l_G(cell_temp(i)),face(i))%cell = SDD_cells_g2l_G(neigh_cell)
        adjacency_list(SDD_cells_g2l_G(cell_temp(i)),face(i))%face = neigh_face
      ELSE
        adjacency_list(SDD_cells_g2l_G(cell_temp(i)),face(i))%cell = 0
        adjacency_list(SDD_cells_g2l_G(cell_temp(i)),face(i))%face = 0
      END IF
    END DO

    DEALLOCATE(cell_temp,face)

    ! Close mesh file

    CLOSE(localunit)

  END SUBROUTINE read_tetmesh

END MODULE readmesh_module
