!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! THOR MESH GENERATION UTILITY
!   Boundary Conditions Module:
!
!> This module contains the functionality necessary to create boundary
!! conditions for a given set of elements and
!
!> @author Nicholas Herring
!> @version 1.0
!> @date Mar, 2018
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
MODULE boundary_conditions
    IMPLICIT NONE
CONTAINS

    !-------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------
    !> Get faces for every single tet, make sure they are in order from low to high
    !! in this local list. i.e. 3 1 2 is represented always as 1 2 3. This prevents
    !! redundancy in the adjancency checking protocol.
    !!
    !! @param faces Faces of the elements
    !! @param elements Info on elements
    !! @param numelements Number of elements
    !! @param numproc Number of threads
    !-------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------
    SUBROUTINE getfaces(faces,elements,numelements,numproc)
        !> Number of elements
        INTEGER,INTENT(IN)::numelements
        !> Number of threads
        INTEGER,INTENT(IN)::numproc
        !> Info on elements
        INTEGER,INTENT(IN)::elements(numelements,4)
        !> Faces of the elements
        INTEGER,INTENT(OUT)::faces(numelements,4,3)
        !> Temporary integer used to order faces
        INTEGER::tempint=0
        !> Loop control
        INTEGER::i,j

        !loop over all elements, this is not expensive but can be parallel anyway
        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,tempint) NUM_THREADS(numproc)
        !$OMP DO
        DO i=1,numelements
            !loop over all faces for a given element
            DO j=1,4
                !gets data for nodes of a given face of an element
                faces(i,j,1)=elements(i,1+FLOOR(j/4.))
                faces(i,j,2)=elements(i,2+FLOOR(j/3.))
                faces(i,j,3)=elements(i,3+FLOOR((j+1)/3.))
                !orders face nodes properly for the given element, could write a sorting algorithm, but there are only three nodes per face so it's just not worth it.
                IF((faces(i,j,2) .LT. faces(i,j,1)) .AND. (faces(i,j,2) .LT. faces(i,j,3)))THEN
                    tempint=faces(i,j,1)
                    faces(i,j,1)=faces(i,j,2)
                    faces(i,j,2)=tempint
                ELSE IF((faces(i,j,3) .LT. faces(i,j,1)) .AND. (faces(i,j,3) .LT. faces(i,j,2)))THEN
                    tempint=faces(i,j,1)
                    faces(i,j,1)=faces(i,j,3)
                    faces(i,j,3)=tempint
                END IF
                IF(faces(i,j,3) .LT. faces(i,j,2))THEN
                    tempint=faces(i,j,2)
                    faces(i,j,2)=faces(i,j,3)
                    faces(i,j,3)=tempint
                END IF
            END DO
        END DO
        !$OMP END DO
        !$OMP END PARALLEL
    END SUBROUTINE getfaces

    !-------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------
    !> Find the boundary faces for a general set of faces for a general set of
    !! elements, see
    !!
    !! @param faces Faces for each element, contains their nodes
    !! @param numel Number of elements
    !! @param boundfaces Boundary faces nodes
    !! @param totbound Total number of boundary faces
    !! @param numproc Number of processors
    !-------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------
    SUBROUTINE findboundaries(faces,numel,boundfaces,totbound,numproc)
        IMPLICIT NONE
        !> Number of elements
        INTEGER,INTENT(IN)::numel
        !> Number of processors
        INTEGER,INTENT(IN)::numproc
        !> Faces for each element, contains their nodes
        INTEGER,INTENT(IN)::faces(numel,4,3)
        !> Total number of boundary faces
        INTEGER,INTENT(OUT)::totbound
        !> Boundary faces nodes
        INTEGER,ALLOCATABLE,INTENT(OUT)::boundfaces(:,:)
        !> Adjacency count for each face
        INTEGER::fadjcount(numel,4)
        !> Loop control variables
        INTEGER::i,j,k

        !get adjacencies for the faces
        CALL getadjacencies(faces,numel,fadjcount,numproc)

        !at this point any boundary faces will have an adjacency count of 0, and all other faces will have adjacency of 1
        !here we tally up the boundary faces, and as a sanity check inform the user of how many faces are boundaries vs how many faces there are total
        totbound=0
        DO i=1,numel
            DO k=1,4
                IF(fadjcount(i,k) .EQ. 0)THEN
                    totbound=totbound+1
                END IF
            END DO
        END DO
        WRITE(*,'(2(I0,A))')totbound,' boundary faces, ',4*numel,' faces total'

        !assign all boundary faces, only assign if the adjacency count for that face is 0
        j=0
        ALLOCATE(boundfaces(totbound,3))
        DO i=1,numel
            DO k=1,4
                IF(fadjcount(i,k) .EQ. 0)THEN
                    j=j+1
                    boundfaces(j,1)=faces(i,k,1)
                    boundfaces(j,2)=faces(i,k,2)
                    boundfaces(j,3)=faces(i,k,3)
                END IF
            END DO
        END DO
    END SUBROUTINE findboundaries

    !-------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------
    !> Calculates the adjacency count for each face, this is the expensive part,
    !! optimization attempts should focus here.
    !!
    !! @param faces Faces for each element, contains their nodes
    !! @param numel Number of elements
    !! @param fadjcount Adjacency count for each face
    !! @param numproc Number of processors
    !-------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------
    SUBROUTINE getadjacencies(faces,numel,fadjcount,numproc)
        IMPLICIT NONE
        !> Number of elements
        INTEGER,INTENT(IN)::numel
        !> Number of processors
        INTEGER,INTENT(IN)::numproc
        !> Faces for each element, contains their nodes
        INTEGER,INTENT(IN)::faces(numel,4,3)
        !> Adjacency count for each face
        INTEGER,INTENT(OUT)::fadjcount(numel,4)
        !> Timer variables
        INTEGER(8)::time1,time2,clock_rate
        REAL(8)::timetaken
        !> Loop control variables
        INTEGER::i,j,k,l

        !preassign all adjacencies for each elements face as -1
        fadjcount=-1

        ! when an adjacency hit is registered, the adjacency count will increment by 1.
        !One of the hits will be the face originally pulled from, thus the -1 pre-assignment
        ! loop over all faces for the checked face
        CALL SYSTEM_CLOCK(time1)
        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,k,j,l) NUM_THREADS(numproc)
        !$OMP DO
        DO i=1,numel
            DO k=1,4
                !loop over all faces to check against the face
                DO j=i,numel
                    IF(fadjcount(i,k) .GT. 0)EXIT
                    DO l=1,4
                        IF(fadjcount(i,k) .GT. 0)EXIT
                        !check against all nodes of a face, if those same nodes are all found on another face it has an adjacency and is not a boundary face
                        IF((faces(i,k,1) .EQ. faces(j,l,1)) .AND. (faces(i,k,2) .EQ. &
                        faces(j,l,2)) .AND. (faces(i,k,3) .EQ. faces(j,l,3)))THEN
                            fadjcount(i,k)=fadjcount(i,k)+1
                        END IF
                    END DO
                END DO
                IF(i .LE. 1)EXIT
                DO j=i-1,1,-1
                    IF(fadjcount(i,k) .GT. 0)EXIT
                    DO l=1,4
                        IF(fadjcount(i,k) .GT. 0)EXIT
                        !check against all nodes of a face, if those same nodes are all found on another face it has an adjacency and is not a boundary face
                        IF((faces(i,k,1) .EQ. faces(j,l,1)) .AND. (faces(i,k,2) .EQ. &
                        faces(j,l,2)) .AND. (faces(i,k,3) .EQ. faces(j,l,3)))THEN
                            fadjcount(i,k)=fadjcount(i,k)+1
                        END IF
                    END DO
                END DO
            END DO
        END DO
        !$OMP END DO
        !$OMP END PARALLEL
        CALL SYSTEM_CLOCK(time2,clock_rate)
        timetaken=(time2*1.0_8-time1*1.0_8)/(clock_rate*1.0_8)

        !output time taken because this can be a fairly expensive code for large problems
        WRITE(*,'(A,ES14.8,A)')'Time taken calculating adjacencies is ', timetaken , ' seconds.'
    END SUBROUTINE getadjacencies
END MODULE boundary_conditions
