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
    INTEGER, ALLOCATABLE :: bc_side(:)
CONTAINS

  SUBROUTINE adjacency_calc()
    INTEGER :: i,j,og_face(3),comp_face(3),adj_idx

    DO i=1,num_tets
      element(i,:)=orderedverts(element(i,:))
    ENDDO

    ALLOCATE(adj_list(num_tets*4,4),tbound_cond(num_tets*4,3))
    adj_list=0
    tbound_cond=0
    !loop over all tets
    adj_idx=0
    num_bcf=0
    prog=0
    WRITE(*,'(A)',ADVANCE='NO')'Progress:'
    DO i=1,num_tets
      IF(MOD(i,CEILING(num_tets*1.0/(max_prog-1.0))) .EQ. 0)THEN
        WRITE(*,'(A)',ADVANCE='NO')'*'
        prog=prog+1
      ENDIF
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
    ALLOCATE(bc_data(num_bcf,3),bc_side(num_bcf))
    bc_data=0
    bc_side=0
    DO i=1,num_bcf
      bc_data(i,:)=tbound_cond(i,:)
    ENDDO

    bc_side=0
    !determine which sides are flat
    IF(MINVAL(side_bc) .EQ. MAXVAL(side_bc))THEN
      bc_data(:,3)=side_bc(1)
    ELSE
      CALL det_side_flatness()
    ENDIF
    DO i=prog,max_prog
      WRITE(*,'(A)',ADVANCE='NO')'*'
    ENDDO
    WRITE(*,*)
    DEALLOCATE(tbound_cond,bc_side)
  ENDSUBROUTINE adjacency_calc

  SUBROUTINE find_adj(face,el_idx,faceid,adj_idx)
    INTEGER,INTENT(IN) :: face(3)
    INTEGER,INTENT(IN) :: el_idx
    INTEGER,INTENT(IN) :: faceid
    INTEGER,INTENT(INOUT) :: adj_idx
    INTEGER :: j,comp_face(3),side_id
    LOGICAL :: match,flat

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
      !assign bc to 0 for now (this will change later
      tbound_cond(num_bcf,3)=0
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

  SUBROUTINE det_side_flatness()
    INTEGER :: i,j,el_id,vert_id,face_idx,loc_max,loc_min
    REAL(8) :: face_point(3,3),ext_point(3),norm_vec(3),lambda,offset

    !sides are assumed to be flat, and if they are not this is set to false.
    side_flat=.TRUE.
    DO i=1,num_bcf
      el_id=bc_data(i,1)
      face_idx=0
      !assign extruded and face boundary points
      DO j=1,4
        IF(bc_data(i,2) .NE. j-1)THEN
          face_idx=face_idx+1
          face_point(face_idx,:)=vertex(element(el_id,j),:)
        ELSE
          ext_point(:)=vertex(element(el_id,j),:)
        ENDIF
      ENDDO
      !get the outward going unit normal vector for the tet for this face
      norm_vec=cross(face_point(2,:)-face_point(1,:), face_point(3,:)-face_point(1,:))
      offset=face_point(1,1)*norm_vec(1)+face_point(1,2)*norm_vec(2)+face_point(1,3)*norm_vec(3)
      lambda=(offset-norm_vec(1)*ext_point(1)-norm_vec(2)*ext_point(2)-norm_vec(3)*ext_point(3)) &
        /(norm_vec(1)**2+norm_vec(2)**2+norm_vec(3)**2)
      norm_vec=norm_vec*lambda
      norm_vec=norm_vec/(SQRT(norm_vec(1)**2+norm_vec(2)**2+norm_vec(3)**2))

      !figure out which side of the problem this bc faces
      IF(ABS(MAXVAL(norm_vec)) .GT. ABS(MINVAL(norm_vec)))THEN
        !face on a positive side, assign the side for that face
        IF(MAXLOC(norm_vec,1) .EQ. 1)THEN
          bc_side(i)=2
        ELSEIF(MAXLOC(norm_vec,1) .EQ. 2)THEN
          bc_side(i)=4
        ELSE
          bc_side(i)=6
        ENDIF
        IF(ABS(MAXVAL(norm_vec)-1.0) .LE. 1.0E-14)THEN
          !face is flat, do nothing
        ElSE
          !face is not flat, it only takes 1 for the side to not be flat
          side_flat(bc_side(i))=.FALSE.
        ENDIF
      ELSE
        !face on a negative side, assign the side for that face
        IF(MINLOC(norm_vec,1) .EQ. 1)THEN
          bc_side(i)=1
        ELSEIF(MINLOC(norm_vec,1) .EQ. 2)THEN
          bc_side(i)=3
        ELSE
          bc_side(i)=5
        ENDIF
        IF(ABS(MINVAL(norm_vec)+1.0) .LE. 1.0E-14)THEN
          !face is flat, do nothing
        ElSE
          !face is not flat, it only takes 1 for the side to not be flat
          side_flat(bc_side(i))=.FALSE.
        ENDIF
      ENDIF
    ENDDO
    DO i=1,6
      IF(.NOT. side_flat(i) .AND. side_bc(i) .EQ. 1)THEN
        WRITE(*,*)'ERROR: side ',i,' is reflective but not flat'
        STOP
      ENDIF
    ENDDO
    DO i=1,num_bcf
      bc_data(i,3)=side_bc(bc_side(i))
    ENDDO
  ENDSUBROUTINE det_side_flatness

  FUNCTION cross(a, b)
    REAL(8) :: cross(3)
    REAL(8), INTENT(IN) :: a(3), b(3)

    cross(1) = a(2) * b(3) - a(3) * b(2)
    cross(2) = a(3) * b(1) - a(1) * b(3)
    cross(3) = a(1) * b(2) - a(2) * b(1)
  END FUNCTION cross

  FUNCTION orderedverts(verts)
    INTEGER :: orderedverts(4)
    INTEGER, INTENT(IN) :: verts(4)
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
