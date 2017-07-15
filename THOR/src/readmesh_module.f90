module readmesh_module
!***********************************************************************
!
! Read mesh module contains all subroutines needed for reading 
! tetrahedral mesh file
!
!***********************************************************************
  use mpi
  use types
  use parameter_types
  use filename_types
  use vector_types
  use geometry_types
  use global_variables
  use termination_module

  implicit none

contains

  subroutine read_tetmesh
  !*********************************************************************
  !
  ! Subroutine reads tetrahedral mesh from CUBIT Exodus II file
  !
  !*********************************************************************

  ! Define temporary variables

    integer(kind=li) :: alloc_stat, i, j, k, dummy
    integer(kind=li) :: kv,kr,kf
    integer(kind=li), dimension(:), allocatable:: cell_temp, face, &
         side_cells_tmp
    real(kind=d_t) :: rdummy
    integer ::rank,mpi_err, localunit
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)
    localunit = rank+100

  ! Set minreg,maxreg

    minreg= 100000_li
    maxreg=-1_li

  ! Open and read mesh file 

    open(unit=localunit,file=trim(mesh_filename),status='old',action='read')

  ! Read general parameters

    read(localunit,*) num_vert
    read(localunit,*) num_cells
    read(localunit,*) num_cell_blk
    read(localunit,*) num_side_sets

  ! Allocate vertices type dimension

    allocate(vertices(num_vert),stat=alloc_stat)
    if(alloc_stat /= 0) call stop_thor(2_li)

  ! Read vertices

    do i=1, num_vert
       read(localunit,*) dummy, vertices(i)%v
    end do

  ! Allocate cell type dimension

    allocate(cells(num_cells),stat=alloc_stat)
    if(alloc_stat /= 0) call stop_thor(2_li)

  ! Read cell block property list (region and source composition)

    do i=1, num_cells
       read(localunit,*) dummy, cells(i)%reg, cells(i)%src
       if(cells(i)%reg>maxreg) maxreg=cells(i)%reg
       if(cells(i)%reg<minreg) minreg=cells(i)%reg
    end do

  ! Read vertices to cell mapping

    do i=1, num_cells
       read(localunit,*) dummy, cells(i)%R(0), cells(i)%R(1), &
            cells(i)%R(2), cells(i)%R(3)
    end do

  ! Go through side_sets once to count side_cells

    side_cells =0
    vside_cells=0
    rside_cells=0
    fside_cells=0

    allocate(side_cells_tmp(num_side_sets),stat=alloc_stat)
    if(alloc_stat /= 0) call stop_thor(2_li)
    
    do i=1, num_side_sets
       
       read(localunit,*) side_cells_tmp(i)
      
       do j=1, side_cells_tmp(i)
          read(localunit,*) dummy,dummy,k
          if      (k.eq.0) then 
             vside_cells=vside_cells+1
          else if (k.eq.1) then
             rside_cells=rside_cells+1
          else if (k.eq.2) then
             fside_cells=fside_cells+1
          end if   
       end do
       
       side_cells=side_cells+side_cells_tmp(i)

    end do

  ! If there are fixed inflow boundary faces then finflow flag must be ==1

    if(fside_cells>0 .and. finflow .eq. 0) then
       call stop_thor(12_li)
    end if

  ! Close file, reopen and read everything up to side_sets into dummies

    close(unit=localunit)
 
    open(unit=localunit,file=trim(mesh_filename),status='old',action='read')

    read(localunit,*) dummy
    read(localunit,*) dummy
    read(localunit,*) dummy
    read(localunit,*) dummy

    do i=1, num_vert
       read(localunit,*) dummy, rdummy,rdummy,rdummy
    end do

    do i=1, num_cells
       read(localunit,*) dummy, dummy, dummy
    end do

    do i=1, num_cells
       read(localunit,*) dummy, dummy,dummy,dummy,dummy
    end do

  ! Allocate b_cells and read from side sets  
    
    allocate(b_cells(side_cells),stat=alloc_stat)
    if(alloc_stat /= 0) call stop_thor(2_li)

    if(vside_cells .gt. 0) then
      allocate(vb_cells(vside_cells),stat=alloc_stat)
      if(alloc_stat /= 0) call stop_thor(2_li)
    end if

    if(rside_cells .gt. 0) then
      allocate(rb_cells(rside_cells),stat=alloc_stat)
      if(alloc_stat /= 0) call stop_thor(2_li)
    end if

    if(fside_cells .gt. 0) then
      allocate(fb_cells(fside_cells),stat=alloc_stat)
      if(alloc_stat /= 0) call stop_thor(2_li)
    end if

    k=1 ; kv=1 ; kr=1 ; kf=1

    do i=1, num_side_sets
       read(localunit,*) dummy
       do j=1, side_cells_tmp(i)
          read(localunit,*) b_cells(k)%cell,b_cells(k)%face,b_cells(k)%bc 
          ! work on boundary data
          if      (b_cells(k)%bc.eq.0) then
             vb_cells(kv)%cell=b_cells(k)%cell 
             vb_cells(kv)%face=b_cells(k)%face 
             vb_cells(kv)%bc  =0 
             vb_cells(kv)%ptr =k                 ! ptr points to the position in the b_cells array 
             b_cells(k)%ptr   =kv                ! ptr points from kth position in b_cells to kv position in vb_cells array
             kv=kv+1      
          else if (b_cells(k)%bc.eq.1) then
             rb_cells(kr)%cell=b_cells(k)%cell 
             rb_cells(kr)%face=b_cells(k)%face 
             rb_cells(kr)%bc  =1
             rb_cells(kr)%ptr =k
             b_cells(k)%ptr   =kr
             kr=kr+1 
          else if (b_cells(k)%bc.eq.2) then
             fb_cells(kf)%cell=b_cells(k)%cell 
             fb_cells(kf)%face=b_cells(k)%face 
             fb_cells(kf)%bc  =2 
             fb_cells(kf)%ptr =k
             b_cells(k)%ptr   =kf
             kf=kf+1
          end if 
          k=k+1
       end do
    end do

    deallocate(side_cells_tmp)
    
    ! Read adjacency list

    read(localunit,*) adjacent_cells
    allocate(adjacency_list(adjacent_cells,0:3),&
         cell_temp(adjacent_cells),face(adjacent_cells),&
         stat=alloc_stat)
    if(alloc_stat /= 0) call stop_thor(2_li)

    do i=1, adjacent_cells
       read(localunit,*) cell_temp(i), face(i),&
            adjacency_list(cell_temp(i),face(i))%cell, &
            adjacency_list(cell_temp(i),face(i))%face
    end do

    deallocate(cell_temp,face)

  ! Close mesh file

    close(localunit)

  end subroutine read_tetmesh

end module readmesh_module

