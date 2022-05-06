!-------------------------------------------------------------------------------
!> This module contains the functionality necessary to create boundary
!! conditions and adjacency info for a given set of elements
!> @author Nicholas F. Herring
!-------------------------------------------------------------------------------
MODULE boundary_conditions
    USE globals
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: adjacency_calc

    INTEGER, ALLOCATABLE :: tbound_cond(:,:)
CONTAINS

  SUBROUTINE adjacency_calc()
    INTEGER :: i,j,og_face(3),comp_face(3),adj_idx

    DO i=1,num_tets
      element(i,:)=orderedverts(element(i,:))
    ENDDO

    ALLOCATE(adj_list(num_tets*4,4),tbound_cond(num_tets*4,3))
    adj_list=0
    !loop over all tets
    adj_idx=0
    num_bcf=0
    DO i=1,num_tets
      !first face
      og_face=(/element(i,2),element(i,3),element(i,4)/)
      CALL find_adj(og_face,i,0,adj_idx)
      !second face
      og_face=(/element(i,1),element(i,3),element(i,4)/)
      CALL find_adj(og_face,i,1,adj_idx)
      !third face
      og_face=(/element(i,1),element(i,2),element(i,4)/)
      CALL find_adj(og_face,i,2,adj_idx)
      !fourth face
      og_face=(/element(i,1),element(i,2),element(i,3)/)
      CALL find_adj(og_face,i,3,adj_idx)
    ENDDO
    ALLOCATE(bc_data(num_bcf,3))
    DO i=1,num_bcf
      bc_data(i,:)=tbound_cond(i,:)
    ENDDO
    DEALLOCATE(tbound_cond)
  ENDSUBROUTINE adjacency_calc

  SUBROUTINE find_adj(face,el_idx,faceid,adj_idx)
    INTEGER,INTENT(IN) :: face(3)
    INTEGER,INTENT(IN) :: el_idx
    INTEGER,INTENT(IN) :: faceid
    INTEGER,INTENT(INOUT) :: adj_idx
    INTEGER :: j,comp_face(3)
    LOGICAL :: match

    match=.FALSE.
    DO j=1,num_tets
      !compare for first face
      comp_face=(/element(j,2),element(j,3),element(j,4)/)
      CALL check_face(face,comp_face,el_idx,j,faceid,0,adj_idx,match)
      IF(match)EXIT
      !compare for second face
      comp_face=(/element(j,1),element(j,3),element(j,4)/)
      CALL check_face(face,comp_face,el_idx,j,faceid,1,adj_idx,match)
      IF(match)EXIT
      !compare for third face
      comp_face=(/element(j,1),element(j,2),element(j,4)/)
      CALL check_face(face,comp_face,el_idx,j,faceid,2,adj_idx,match)
      IF(match)EXIT
      !compare for fourth face
      comp_face=(/element(j,1),element(j,2),element(j,3)/)
      CALL check_face(face,comp_face,el_idx,j,faceid,3,adj_idx,match)
      IF(match)EXIT
    ENDDO
    !if we go through the whole thing and don't exit
    IF(j .EQ. num_tets+1)THEN
      !we didn't find a matching face so it's a boundary condition
      adj_idx=adj_idx+1
      adj_list(adj_idx,1)=el_idx
      adj_list(adj_idx,2)=faceid
      adj_list(adj_idx,3)=0
      adj_list(adj_idx,4)=0
      num_bcf=num_bcf+1
      tbound_cond(num_bcf,1)=el_idx
      tbound_cond(num_bcf,2)=faceid
      IF(MINVAL(side_bc(:)) .EQ. MAXVAL(side_bc(:)))THEN
        tbound_cond(num_bcf,3)=0
      ElSE
        CALL determine_side()
      ENDIF
    ENDIF
  ENDSUBROUTINE find_adj

  SUBROUTINE check_face(face1,face2,el_idx1,el_idx2,faceid1,faceid2,adj_idx,match)
    INTEGER,INTENT(IN) :: face1(3)
    INTEGER,INTENT(IN) :: face2(3)
    INTEGER,INTENT(IN) :: el_idx1
    INTEGER,INTENT(IN) :: el_idx2
    INTEGER,INTENT(IN) :: faceid1
    INTEGER,INTENT(IN) :: faceid2
    INTEGER,INTENT(INOUT) :: adj_idx
    LOGICAL,INTENT(INOUT) :: match
    match=.FALSE.
    !check to see if the faces match
    IF(face1(1) .EQ. face2(1) .AND. face1(2) .EQ. face2(2) &
        .AND. face1(3) .EQ. face2(3) .AND. el_idx1 .NE. el_idx2)THEN
      adj_idx=adj_idx+1
      adj_list(adj_idx,1)=el_idx1
      adj_list(adj_idx,2)=faceid1
      adj_list(adj_idx,3)=el_idx2
      adj_list(adj_idx,4)=faceid2
      match=.TRUE.
    ENDIF
  ENDSUBROUTINE check_face

  SUBROUTINE determine_side()
    STOP 'need to determine the side that this boundary face is on'
  ENDSUBROUTINE determine_side

  FUNCTION orderedverts(verts)
    INTEGER, INTENT(IN) :: verts(4)
    INTEGER :: orderedverts(4)
    INTEGER :: i,temp_vert,changes
    orderedverts=verts
    !bubble sort algorithm, pretty cheap for only 4 elements
    DO
      changes=0
      DO i=1,3
        IF(orderedverts(i) .GE. orderedverts(i+1))THEN
          temp_vert=orderedverts(i)
          orderedverts(i)=orderedverts(i+1)
          orderedverts(i+1)=temp_vert
          changes=changes+1
        ENDIF
      ENDDO
      IF(changes .EQ. 0)EXIT
    ENDDO
  ENDFUNCTION
END MODULE boundary_conditions
