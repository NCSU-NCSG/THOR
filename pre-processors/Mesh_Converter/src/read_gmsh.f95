!-------------------------------------------------------------------------------
!> This module contains the functionality necessary to ingest a file in the
!! Gmesh 3.0 format (.msh)
!> @author Nicholas F. Herring
!-------------------------------------------------------------------------------
MODULE read_gmsh
  USE globals
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_gmsh_file
CONTAINS
  SUBROUTINE read_gmsh_file()
    !temp variables for non stored data or data in transition
    CHARACTER(300) :: temp_str

    !error variable
    INTEGER :: ios

    OPEN(unit=20,FILE=mesh_infile,ACTION='READ',STATUS='OLD')

    READ(20,*)temp_str
    IF(temp_str .NE. '$MeshFormat')STOP 'Not a gmesh file'
    READ(20,*)temp_str
    IF(temp_str .NE. '4.1')STOP 'gmesh format is wrong, should be 4.1'

    !find the nodes
    REWIND(20)
    DO
      READ(20,*,iostat=ios)temp_str
      IF(temp_str .EQ. '$Nodes')THEN
        !and read them
        CALL read_nodes()
        EXIT
      ENDIF
      IF(ios .NE. 0)STOP 'did not find Nodes in gmesh file'
    ENDDO

    !find the elements
    REWIND(20)
    DO
      READ(20,*,iostat=ios)temp_str
      IF(temp_str .EQ. '$Elements')THEN
        !and read them
        CALL read_elements()
        EXIT
      ENDIF
      IF(ios .NE. 0)STOP 'did not find Elements in gmesh file'
    ENDDO
    CLOSE (20)
  ENDSUBROUTINE read_gmsh_file

  SUBROUTINE read_nodes()
    INTEGER :: num_entities,minnode,maxnode
    INTEGER :: temp_int,loc_num_nodes,i,j,node_indx
    !get the nodes info
    READ(20,*)num_entities,num_verts,minnode,maxnode
    IF(minnode .NE. 1 .OR. maxnode .NE. num_verts)STOP 'nodes should be indexed 1 to amount'
    ALLOCATE(vertex(num_verts,3))
    vertex=0.0
    !read in all node data
    node_indx=1
    DO i=1,num_entities
      READ(20,*)temp_int,temp_int,temp_int,loc_num_nodes
      !get past the indices
      DO j=1,loc_num_nodes
        READ(20,*)
      ENDDO
      !read in the node datas
      DO j=1,loc_num_nodes
        READ(20,*)vertex(node_indx,:)
        node_indx=node_indx+1
      ENDDO
    ENDDO
  ENDSUBROUTINE read_nodes

  SUBROUTINE read_elements()
    INTEGER :: num_entities,numel,minel,maxel
    INTEGER :: temp_int,loc_num_el,i,j,el_indx,el_dim,ent_tag
    INTEGER, ALLOCATABLE :: temp_array(:,:)
    !get the elements info
    READ(20,*)num_entities,numel,minel,maxel
    !last value of array is tag (region/material)
    ALLOCATE(temp_array(numel,5))
    IF(minel .NE. 1 .OR. maxel .NE. numel)STOP 'elements should be indexed 1 to amount'
    !read in all node data
    num_tets=0
    el_indx=0
    DO i=1,num_entities
      READ(20,*)el_dim,ent_tag,temp_int,loc_num_el
      !actually counting the tets
      IF(el_dim .EQ. 3)THEN
        num_tets=num_tets+loc_num_el
        !get data since this is a tet set
        DO j=1,loc_num_el
          el_indx=el_indx+1
          READ(20,*)temp_int,temp_array(el_indx,1:4)
          temp_array(el_indx,5)=ent_tag
        ENDDO
      ELSE
        !get past element data, just counting tets right now
        DO j=1,loc_num_el
          READ(20,*)
        ENDDO
      ENDIF
    ENDDO
    !allocate elements and element tags
    ALLOCATE(element(num_tets,4),el_tag(num_tets))
    element=0
    el_tag=0
    DO i=1,num_tets
      element(i,:)=temp_array(i,1:4)
      el_tag(i)=temp_array(i,5)
    ENDDO
    DEALLOCATE(temp_array)
  ENDSUBROUTINE read_elements
END MODULE read_gmsh
