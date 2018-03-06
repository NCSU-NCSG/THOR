!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! THOR MESH GENERATION UTILITY
!   UNV Module:
!
!> This module contains the functionality necessary to ingest a file in the
!! UNV format (.unv)
!
!> @author Nicholas Herring
!> @version 1.0
!> @date Mar, 2018
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
MODULE unv
    USE globals
    USE boundary_conditions
    IMPLICIT NONE
CONTAINS

    !-------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------
    !> Extracts elements, node, bc, and block_id data from the unv file, then
    !! converts them to proper format for the globals used in the mesh creator
    !-------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------
    SUBROUTINE ingestUNV()
        !> Number of processors to use
        INTEGER::numproc
        !> All of the faces
        INTEGER,ALLOCATABLE::faces(:,:,:)

        !read in the unv file
        CALL readinunv()

        !set number of processors to 1 for non parallel runs, prompt for number of processors
        numproc=1
        !$ WRITE(*,'(A)')"How many processors do you wish to use for boundary addition?"
        !$ READ(*,*)numproc

        !allocate the faces array
        ALLOCATE(faces(element_count,4,3))

        ! get all the faces
        CALL getfaces(faces,element_list,element_count,numproc)

        !find the boundaries
        CALL findboundaries(faces,element_count,bc_list,bc_count,numproc)

        !there is no side set information in a unv file, because there is no boundary information
        ALLOCATE(side_set(bc_count))
        side_set=0
    END SUBROUTINE ingestUNV

    !-------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------
    !> Extracts elements, node, and material data from the unv file
    !-------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------
    SUBROUTINE readinunv()
        !> Delimiter to find various pieces of the mesh in the unv file
        CHARACTER(64)::delim
        !> Temporary character variable
        CHARACTER(64)::tempcharacter
        !> loop control
        INTEGER::i,j
        !> Conversion factor for unv nodes
        REAL(8)::conversion

        !open input mesh file
        OPEN(UNIT = in_unit, FILE = TRIM(in_file), ACTION = 'READ', STATUS = 'OLD', IOSTAT = io_status, IOMSG=tempcharacter)
        IF(io_status .NE. 0)THEN
            WRITE(*,*)tempcharacter
            STOP
        END IF

        !get to data card
        delim="164"
        CALL finddelim(in_unit,delim)

        !get conversion factor, code outputs units of cm
        READ(in_unit,*)tempcharacter
        SELECT CASE(tempcharacter)
            CASE('1Meter', '3Meter')
                conversion=100.
            CASE('2Foot', '4Foot')
                conversion=30.48
            CASE('5mm', '8mm', '10mm')
                conversion=0.1
            CASE('6cm')
                conversion=1
            CASE('7Inch')
                conversion=2.54
            CASE DEFAULT
                STOP 'Not a valid set of units'
        END SELECT

        !get to nodes card
        delim="2411"
        CALL finddelim(in_unit,delim)

        !count the number of nodes until -1 delimiter
        delim="-1"
        CALL counter(in_unit,delim,node_count)

        !fix count for those extra lines in format
        node_count=(node_count-1)/2

        !get to elements card
        delim="2412"
        CALL finddelim(in_unit,delim)

        !count number of elements until -1 delimiter
        delim="-1"
        CALL counter(in_unit,delim,element_count)

        !fix count for those extra lines in format
        element_count=(element_count-1)/2

        !allocate nodes and elements arrays
        ALLOCATE(node_list(node_count,3),element_list(element_count,4),block_id(element_count))

        !get to nodes card
        delim="2411"
        CALL finddelim(in_unit,delim)

        !read in all nodes data, skipping extra lines
        DO i=1,node_count
            READ(in_unit,*)tempcharacter
            READ(in_unit,*)node_list(i,:)
            !convert those nodes to cm
            DO j=1,3
                node_list(i,j)=node_list(i,j)*conversion
            END DO
        END DO

        !get to elements card
        delim="2412"
        CALL finddelim(in_unit,delim)

        !read in all elements data, extra lines contain relevant data so don't skip them!
        DO i=1,element_count
            READ(in_unit,*)tempcharacter,tempcharacter,tempcharacter,block_id(i)
            READ(in_unit,*)element_list(i,:)
        END DO

        !close input file
        CLOSE(UNIT=in_unit,IOSTAT=io_status)
        IF(io_status .NE. 0)THEN
            WRITE(*,*)'Could not close input file ', in_file
            STOP
        END IF
    END SUBROUTINE readinunv

    !-------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------
    !> Find delimiter value in a unv file
    !!
    !! @param filenum File id number
    !! @param delim Delimiter string
    !-------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------
    SUBROUTINE finddelim(filenum,delim)
        !> File id number
        INTEGER,INTENT(IN)::filenum
        !> Delimiter string
        CHARACTER(64),INTENT(IN)::delim
        !> Temporary character variables
        CHARACTER(64)::tc1="",tc2=""
        !> File io status
        INTEGER::ios
        !> Temporary character variable
        CHARACTER(64)::tempcharacter

        !go to beginning of file
        REWIND(filenum)

        !find delimiter, if not found give error, this is tricky for unv files as the
        !delimiter may appear somewhere and be incorrect, need to also check that it
        !is proceeded by two -1 values
        tempcharacter=""
        DO WHILE(tempcharacter .NE. delim .OR. tc1 .NE. "-1" .OR. tc2 .NE. "-1")
            tc1=tc2
            tc2=tempcharacter
            READ(filenum,*,IOSTAT=ios)tempcharacter
            IF(ios .NE. 0)THEN
                WRITE(*,'(2A)') "End of file reached without finding delimiter: ", delim
                STOP "Error: delimiter not found!"
            END IF
        END DO
    END SUBROUTINE finddelim

    !-------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------
    !> Count the number of values in the section of the unv file
    !!
    !! @param filenum File id number
    !! @param delim Delimiter string
    !! @param count Number counted
    !-------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------
    SUBROUTINE counter(filenum,delim,count)
        !> File id number
        INTEGER,INTENT(IN)::filenum
        !> Delimiter string
        CHARACTER(64),INTENT(IN)::delim
        !> Number counted
        INTEGER,INTENT(OUT)::count
        !> Temporary character variable
        CHARACTER(64)::tempcharacter

        !count number until delimiter
        count=0
        tempcharacter=""
        DO WHILE(tempcharacter .NE. delim)
            READ(filenum,*)tempcharacter
            count=count+1
        END DO
    END SUBROUTINE counter
END MODULE unv
