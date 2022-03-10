MODULE vtk_module

  IMPLICIT NONE

  CONTAINS
    SUBROUTINE generate_vtk(meshfile, partfile, vtkfile)
    IMPLICIT NONE

    CHARACTER(*), INTENT(IN)::meshfile,partfile,vtkfile
    INTEGER::Nnodes,Nelem,temp,i,j, num_side_sets, BC_cell, BC_type
    REAL(8),ALLOCATABLE,DIMENSION(:,:)::nodexyz
    INTEGER,ALLOCATABLE,DIMENSION(:)::part
    INTEGER,ALLOCATABLE,DIMENSION(:,:)::elemnodes
    INTEGER, ALLOCATABLE::src(:), region(:), BC(:), BC_count(:)
    CHARACTER:: com_check

    !Read in mesh file
    OPEN(UNIT=31,FILE=meshfile,STATUS='OLD')
    READ(31,*)Nnodes
    READ(31,*)Nelem
    READ(31,*)temp
    READ(31,*)num_side_sets
    !Element coordinates
    ALLOCATE(nodexyz(Nnodes,3))
    READ(31,*)com_check
    IF (com_check.NE.'!') BACKSPACE(31)
    DO i=1,Nnodes
        READ(31,*)temp,nodexyz(i,1),nodexyz(i,2),nodexyz(i,3)
    END DO
    !Read Region & SRC
    ALLOCATE(src(Nelem), region(Nelem))
    READ(31,*)com_check
    IF (com_check.NE.'!') BACKSPACE(31)
    DO i=1,Nelem
        READ(31,*)temp, region(i), src(i)
    END DO
    !Read element node list
    ALLOCATE(elemnodes(Nelem,4))
    READ(31,*)com_check
    IF (com_check.NE.'!') BACKSPACE(31)
    DO i=1,Nelem
        READ(31,*)temp,elemnodes(i,1),elemnodes(i,2),elemnodes(i,3),elemnodes(i,4)
    END DO
    !Read in all BC data
    ALLOCATE(BC_count(num_side_sets))
    READ(31,*)com_check
    IF (com_check.NE.'!') BACKSPACE(31)
    DO i = 1, num_side_sets
      READ(31,*) BC_count(i)
      DO j = 1, BC_count(i)
        READ(31,*)
      END DO
    END DO
    DO i = 1, SUM(BC_count) + num_side_sets
      BACKSPACE(31)
    END DO
    ALLOCATE(BC(Nelem))
    DO i = 1, num_side_sets
      READ(31,*) temp
      DO j = 1, BC_count(i)
        READ(31,*) BC_cell, temp, BC_type
        BC(BC_cell) = BC_type
      END DO
    END DO
    CLOSE(31)

    !Read the partition file
    ALLOCATE(part(Nelem))
    OPEN(UNIT=32,FILE=partfile,STATUS='OLD')
    DO i=1,Nelem
        READ(32,*)part(i)
    END DO
    CLOSE(32)

    !Write out vtk file
    OPEN(UNIT=33,FILE=vtkfile,STATUS='REPLACE')
    WRITE(33,'(A)')'# vtk DataFile Version 2.0'
    WRITE(33,'(A)')meshfile
    WRITE(33,'(A)')'ASCII'
    WRITE(33,'(A)')'DATASET UNSTRUCTURED_GRID'
    WRITE(33,'(A,I0,A)')'POINTS ',Nnodes,' float'
    !Node coordinates
    DO i=1,Nnodes
        WRITE(33,*)nodexyz(i,1),nodexyz(i,2),nodexyz(i,3)
    END DO
    WRITE(33,'(A,I0,A,I0)')'CELLS ',Nelem,' ',Nelem*5
    !Element node list
    DO i=1,Nelem
        WRITE(33,'(A,4(I10,X))')'4 ',elemnodes(i,1)-1,elemnodes(i,2)-1,elemnodes(i,3)-1,elemnodes(i,4)-1
    END DO
    WRITE(33,'(A,I0)')'CELL_TYPES ',Nelem
    !Write a ton of tens
    DO i=1,Nelem
        WRITE(33,'(A)',ADVANCE='NO')'10 '
        IF (MOD(i,25) .EQ. 0) WRITE(33,*)
    END DO
    IF (MOD(i,25) .NE. 0)WRITE(33,*)
    WRITE(33,'(A,I0)')'CELL_DATA ',Nelem
    WRITE(33,'(A)')'SCALARS partitioning int 1'
    WRITE(33,'(A)')'LOOKUP_TABLE partitioning'
    !Write the partition list
    DO i=1,Nelem
        WRITE(33,'(I0)')part(i)
    END DO
    WRITE(33,'(A)')'SCALARS region int 1'
    WRITE(33,'(A)')'LOOKUP_TABLE region'
    !Write the partition list
    DO i=1,Nelem
        WRITE(33,'(I0)')region(i)
    END DO
    WRITE(33,'(A)')'SCALARS source int 1'
    WRITE(33,'(A)')'LOOKUP_TABLE source'
    !Write the partition list
    DO i=1,Nelem
        WRITE(33,'(I0)')src(i)
    END DO
    WRITE(33,'(A)')'SCALARS BCs int 1'
    WRITE(33,'(A)')'LOOKUP_TABLE BCs'
    !Write the partition list
    DO i=1,Nelem
        WRITE(33,'(I0)')BC(i)
    END DO
    CLOSE(33)

  END SUBROUTINE

END MODULE
