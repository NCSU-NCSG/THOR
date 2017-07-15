program main

!*****************************************************************************80
!
!! MAIN is the main program for TET_MESH_TET_NEIGHBORS.
!
!  Discussion:
!
!    TET_MESH_TET_NEIGHBORS manages the tet mesh neighbor calculation.
!
!    A tet mesh of order 4 or order 10 may be used.
!
!  Usage:
!
!    tet_mesh_tet_neighbors prefix
!
!    where
!
!    * prefix_elements.txt contains the element definitions;
!    * prefix_element_neighbors.txt contains the neighbor information.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 August 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer ( kind = 4 ) arg_num
  integer ( kind = 4 ) dim_num
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iarg
  integer ( kind = 4 ) iargc
  character ( len = 255 ) :: element_filename = ' '
  integer ( kind = 4 ) node_num
  character ( len = 255 ) :: neighbor_filename = ' '
  character ( len = 255 ) prefix
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: tetra_neighbor
  integer ( kind = 4 ), allocatable, dimension ( :, : ) :: tetra_node
  integer ( kind = 4 ) tetra_num
  integer ( kind = 4 ) tetra_order

  write ( *, '(a)' ) ' '
  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TET_MESH_TET_NEIGHBORS:'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Read a tet mesh dataset of TETRA_NUM'
  write ( *, '(a)' ) '  tetrahedrons, using 4 or 10 nodes.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Compute tet mesh neighbors and write this'
  write ( *, '(a)' ) '  information to a file'
!
!  Get the number of command line arguments.
!
  arg_num = iargc ( )
!
!  Argument 1 is the common file prefix.
!
  if ( 1 <= arg_num ) then

    iarg = 1
    call getarg ( iarg, prefix )

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TET_MESH_TET_NEIGHBORS:'
    write ( *, '(a)' ) '  Please enter the filename prefix:'

    read ( *, '(a)' ) prefix

  end if
!
!  Create the filenames.
!
  element_filename = trim ( prefix ) // '_elements.txt'
  neighbor_filename = trim ( prefix ) // '_element_neighbors.txt'
!
!  Read the tetra data.
!
  call i4mat_header_read ( element_filename, tetra_order, &
    tetra_num )

  if ( tetra_order /= 4 .and. tetra_order /= 10 ) then
    write ( *, * ) ' '
    write ( *, '(a)' ) 'TET_MESH_TET_NEIGHBORS - Fatal error!'
    write ( *, '(a)' ) '  Data is not for a 4 or 10 node tet mesh.'
    stop
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read the header of "' &
    // trim ( element_filename ) //'".'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Tetrahedron order = ', tetra_order
  write ( *, '(a,i8)' ) '  Number of tetras  = ', tetra_num

  allocate ( tetra_node(1:tetra_order,1:tetra_num) )
  allocate ( tetra_neighbor(1:4,1:tetra_num) )

  call i4mat_data_read ( element_filename, tetra_order, &
    tetra_num, tetra_node )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Read the data in "' &
    // trim ( element_filename ) //'".'

  call i4mat_transpose_print_some ( tetra_order, tetra_num, &
    tetra_node, 1, 1, tetra_order, 5, '  First 5 tetrahedrons:' )
!
!  Compute the neighbor information.
!
  call tet_mesh_neighbor_tets ( tetra_order, tetra_num, tetra_node, &
    tetra_neighbor )

  call i4mat_transpose_print_some ( 4, tetra_num, &
    tetra_neighbor, 1, 1, 4, 5, '  First 5 neighbor sets:' )
!
!  Write the neighbor information to a file.
!
  call i4mat_write ( neighbor_filename, 4, tetra_num, tetra_neighbor )

  write ( *, '(a)' ) '  Wrote the tetrahedron neighbor information to "' &
    // trim ( neighbor_filename ) //'".'
!
!  Free memory.
!
  deallocate ( tetra_neighbor )
  deallocate ( tetra_node )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TET_MESH_TET_NEIGHBORS:'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
