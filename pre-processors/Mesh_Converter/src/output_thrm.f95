!-------------------------------------------------------------------------------
!> This module contains the functionality necessary to output a file in the
!! Thor_mesh format (.thrm)
!> @author Nicholas F. Herring
!-------------------------------------------------------------------------------
MODULE output_thrm
  USE globals
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: output_gmsh_file
CONTAINS

  SUBROUTINE output_gmsh_file()
    INTEGER :: i

    OPEN(UNIT=30,FILE=TRIM(ADJUSTL(mesh_infile))//'_out.thrm',ACTION='WRITE',STATUS='REPLACE')

    !print out base data
    WRITE(30,'(I0)')num_verts
    WRITE(30,'(I0)')num_tets
    WRITE(30,'(I0)')1
    WRITE(30,'(I0)')1
    !print out vertices
    DO i=1,num_verts
      WRITE(30,'(I0,3ES24.16)')i,vertex(i,:)
    ENDDO
    !print out tet regions and source. Assume each region has its own source. User can change this later
    DO i=1,num_tets
      WRITE(30,'(I0,A,I0,A,I0)')i,' ',el_tag(i),' ',el_tag(i)
    ENDDO
    !print out tet composition
    DO i=1,num_tets
      WRITE(30,'(I0,A,I0,A,I0,A,I0,A,I0)')i,' ',element(i,1),' ',element(i,2),' ',element(i,3),' ' &
        ,element(i,4)
    ENDDO
    !print out boundary condiitons
    WRITE(30,'(I0)')num_bcf
    DO i=1,num_bcf
      WRITE(30,'(I0,A,I0,A,I0)')bc_data(i,1),' ',bc_data(i,2),' ',bc_data(i,3)
    ENDDO
    !print out adjacency list
    WRITE(30,'(I0)')num_tets*4
    DO i=1,num_tets*4
      WRITE(30,'(I0,A,I0,A,I0,A,I0)')adj_list(i,1),' ',adj_list(i,2),' ',adj_list(i,3),' ',adj_list(i,4)
    ENDDO

    CLOSE(30)
  ENDSUBROUTINE output_gmsh_file
END MODULE output_thrm
