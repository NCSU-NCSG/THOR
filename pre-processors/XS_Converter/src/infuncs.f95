!input functions
MODULE infuncs
  USE globals
  USE stringmod
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: readcl,readxs
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Gets command line arguments
  SUBROUTINE readcl()
    INTEGER::arg_count
    arg_count=COMMAND_ARGUMENT_COUNT()

    IF(arg_count .GT. 2)STOP 'Only two argument variables are allowed for this program! &
        &"xsin" "outformat"'

    !either use or prompt for input file name
    IF(arg_count .GE. 1)THEN
        CALL GET_COMMAND_ARGUMENT(1, xsin)
    ELSE
        WRITE(*,'(A)')'Input cross sections filename?'
        WRITE(*,'(A)',ADVANCE='NO')'> '
        READ(*,*)xsin
    END IF
    xsin=TRIM(ADJUSTL(xsin))

    !either use or prompt for output file name
    IF(arg_count .GE. 2)THEN
        CALL GET_COMMAND_ARGUMENT(2, outformat)
    ELSE
        WRITE(*,'(A)')'Output cross sections format? Available formats below:'
        WRITE(*,'(A)')'THOR'
        WRITE(*,'(A)')'OpenMC'
        WRITE(*,'(A)',ADVANCE='NO')'> '
        READ(*,*)outformat
    END IF
    outformat=TRIM(ADJUSTL(lowercase(outformat)))

    xsout=TRIM(xsin)//'_'//TRIM(outformat)//'.out'
  ENDSUBROUTINE readcl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Get xs format and call the xs reader based on format
  SUBROUTINE readxs()
    INTEGER :: ios,nwords
    CHARACTER(1000) :: tchar1,informat
    CHARACTER(100) :: words(200)

    !open xsin file
    OPEN(UNIT=22,FILE=xsin,STATUS='OLD',ACTION='READ',IOSTAT=ios,IOMSG=tchar1)
    IF(ios .NE. 0)THEN
        WRITE(*,'(A)')tchar1
        STOP
    ENDIF

    !find the input format
    informat=''
    DO
      READ(22,*,IOSTAT=ios)tchar1
      SELECTCASE(tchar1)
        CASE('VERSION')
          BACKSPACE(22)
          READ(22,'(A20000)')tchar1
          CALL parse(tchar1,'=',words,nwords)
          tchar1=TRIM(ADJUSTL(words(2)))
          CALL parse(tchar1,"'",words,nwords)
          tchar1=TRIM(ADJUSTL(words(2)))
          CALL parse(tchar1," ",words,nwords)
          words(1)=TRIM(ADJUSTL(lowercase(words(1))))
          tchar1=TRIM(ADJUSTL(words(2)))
          IF(words(1) .EQ. 'serpent')THEN
            SELECTCASE(tchar1(1:1))
              CASE('1')
                informat='serp_gen_v1'
                EXIT
              CASE('2')
                informat='serp_gen_v2'
                EXIT
            ENDSELECT
          ENDIF
        CASE('THOR_XS_V1')
          informat='thor_v1'
          EXIT
        CASE DEFAULT
      ENDSELECT
      IF(ios .NE. 0)STOP 'No valid format identifier found'
    ENDDO

    !call the appropriate XS reader
    REWIND(22)
    SELECTCASE(informat)
      CASE('serp_gen_v1')
        CALL read_serp_v1()
      CASE('serp_gen_v2')
        CALL read_serp_v2()
      CASE('thor_v1')
        CALL read_thor_v1()
      CASE  DEFAULT
        WRITE(*,'(3A)')'ERROR: ',TRIM(informat),' not a known input xs format.'
        STOP 'Fatal error'
    ENDSELECT
    !close the input file
    CLOSE(22)
  ENDSUBROUTINE readxs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE read_serp_v1()
    STOP 'read_serp_v1 not complete'
  ENDSUBROUTINE read_serp_v1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE read_serp_v2()
    CHARACTER(1000) :: tchar1
    INTEGER :: ios,nwords,m
    CHARACTER(100) :: words(2000)

    numgroups=0
    nummats=0
    levelanis=7
    !find the number of materials and number of energy groups
    DO WHILE(ios .EQ. 0)
      READ(22,*,IOSTAT=ios)tchar1
      !get the number of materials
      IF(tchar1 .EQ. "GC_UNIVERSE_NAME")nummats=nummats+1
      !get the number of energy groups
      IF(tchar1 .EQ. "MACRO_NG")THEN
        IF(numgroups .EQ. 0)THEN
          BACKSPACE(22)
          READ(22,'(A20000)')tchar1
          CALL parse(tchar1,"=",words,nwords)
          READ(words(2),*)numgroups
          ALLOCATE(eg_struc(numgroups+1))
          eg_struc=0.0
        ENDIF
      ENDIF
      !get the energy group structure
      IF(tchar1 .EQ. "MACRO_E")THEN
        IF(ALLOCATED(eg_struc))THEN
          BACKSPACE(22)
          READ(22,'(A20000)')tchar1
          CALL parse(tchar1,"=",words,nwords)
          tchar1=words(2)
          CALL parse(tchar1,"[",words,nwords)
          READ(words(2),*)eg_struc(:)
        ENDIF
      ENDIF
    ENDDO
    REWIND(22)

    ALLOCATE(chi(nummats,numgroups),sigmaf(nummats,numgroups),nuf(nummats,numgroups))
    ALLOCATE(sigmat(nummats,numgroups),sigmas(nummats,levelanis+1,numgroups,numgroups))

    m=0
    DO
      READ(22,*)tchar1
      !found the next material
      IF(tchar1 .EQ. 'GC_UNIVERSE_NAME')THEN
        m=m+1

        !total xs comes first
        CALL getserpv2xsdata('INF_TOT',sigmat(m,:))
        !then fission xs
        CALL getserpv2xsdata('INF_FISS',sigmaf(m,:))
        !then nuf
        CALL getserpv2xsdata('INF_NUBAR',nuf(m,:))
        !then chi
        CALL getserpv2xsdata('INF_CHIT',chi(m,:))
        !then sigmass values for each scattering order
        CALL getserpv2scatdata('INF_SP0',sigmas(m,1,:,:))
        CALL getserpv2scatdata('INF_SP1',sigmas(m,2,:,:))
        CALL getserpv2scatdata('INF_SP2',sigmas(m,3,:,:))
        CALL getserpv2scatdata('INF_SP3',sigmas(m,4,:,:))
        CALL getserpv2scatdata('INF_SP4',sigmas(m,5,:,:))
        CALL getserpv2scatdata('INF_SP5',sigmas(m,6,:,:))
        CALL getserpv2scatdata('INF_SP6',sigmas(m,7,:,:))
        CALL getserpv2scatdata('INF_SP7',sigmas(m,8,:,:))
      ENDIF
      IF(m .GE. nummats)EXIT
    ENDDO
  ENDSUBROUTINE read_serp_v2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE getserpv2xsdata(xsID,xsvec)
    CHARACTER(*),INTENT(IN) :: xsID
    REAL(8) :: xsvec(numgroups)
    CHARACTER(20000) :: tchar1
    INTEGER :: nwords,g
    CHARACTER(100) :: words(2000)

    !get the xs data
    CALL find_line(xsID)
    READ(22,'(A20000)')tchar1
    CALL parse(tchar1,'=',words,nwords)
    tchar1=words(2)
    CALL parse(tchar1,'[',words,nwords)
    tchar1=TRIM(ADJUSTL(words(2)))
    CALL parse(tchar1,' ',words,nwords)
    DO g=1,numgroups
      READ(words(g*2-1),*)xsvec(g)
    ENDDO
  ENDSUBROUTINE getserpv2xsdata

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE getserpv2scatdata(xsID,xsvec)
    CHARACTER(*),INTENT(IN) :: xsID
    REAL(8) :: xsvec(numgroups,numgroups)
    CHARACTER(20000) :: tchar1
    INTEGER :: nwords,g,gp
    CHARACTER(100) :: words(2000)

    !get the xs data
    CALL find_line(xsID)
    READ(22,'(A20000)')tchar1
    CALL parse(tchar1,'=',words,nwords)
    tchar1=words(2)
    CALL parse(tchar1,'[',words,nwords)
    tchar1=TRIM(ADJUSTL(words(2)))
    CALL parse(tchar1,' ',words,nwords)
    DO g=1,numgroups
      DO gp=1,numgroups
        READ(words((g-1)*numgroups*2+gp*2-1),*)xsvec(gp,g)
      ENDDO
    ENDDO
  ENDSUBROUTINE getserpv2scatdata

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE read_thor_v1()
    INTEGER :: nwords,g,m,gp,l
    CHARACTER(100) :: words(2000)
    DO
      CALL get_thor_line(words,nwords)
      words(1)=TRIM(ADJUSTL(words(1)))
      IF(words(1) .EQ. 'THOR_XS_V1')THEN
        EXIT
      ENDIF
    ENDDO
    READ(words(2),*)nummats
    READ(words(3),*)numgroups
    READ(words(4),*)levelanis

    ALLOCATE(chi(nummats,numgroups),sigmaf(nummats,numgroups),nuf(nummats,numgroups))
    ALLOCATE(sigmat(nummats,numgroups),sigmas(nummats,levelanis+1,numgroups,numgroups))
    ALLOCATE(eg_struc(numgroups))
    eg_struc=0.0

    CALL get_thor_line(words,nwords)
    words(1)=TRIM(ADJUSTL(lowercase(words(1))))
    IF(words(1) .NE. 'id')THEN
      DO g=1,numgroups
        READ(words(g),*)eg_struc(g)
      ENDDO
    ENDIF
    REWIND(22)
    !get all material datas
    DO m=1,nummats
      DO
        !get next line
        CALL get_thor_line(words,nwords)
        words(1)=TRIM(ADJUSTL(lowercase(words(1))))
        IF(words(1) .EQ. 'id')THEN
          !read in the fission spectrum
          CALL get_thor_line(words,nwords)
          DO g=1,numgroups
            READ(words(g),*)chi(m,g)
          ENDDO
          !read in SigmaF
          CALL get_thor_line(words,nwords)
          DO g=1,numgroups
            READ(words(g),*)sigmaf(m,g)
          ENDDO
          !read in nu
          CALL get_thor_line(words,nwords)
          DO g=1,numgroups
            READ(words(g),*)nuf(m,g)
          ENDDO
          !read in total/transport xs
          CALL get_thor_line(words,nwords)
          DO g=1,numgroups
            READ(words(g),*)sigmat(m,g)
          ENDDO
          !read in scattering xs format gp->g
          DO l=1,levelanis+1
            DO g=1,numgroups
              CALL get_thor_line(words,nwords)
              DO gp=1,numgroups
                READ(words(gp),*)sigmas(m,l,g,gp)
              ENDDO
            ENDDO
          ENDDO
          EXIT
        ENDIF
      ENDDO
    ENDDO
  ENDSUBROUTINE read_thor_v1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_thor_line(words,nwords)
    INTEGER,INTENT(OUT) :: nwords
    CHARACTER(100),INTENT(OUT) :: words(2000)
    CHARACTER(20000) :: tchar1
    INTEGER :: ios
    DO
      READ(22,'(A10000)',IOSTAT=ios)tchar1
      IF( ios .NE. 0)STOP 'end of xs file was reached before all data/materials were found'
      tchar1=TRIM(ADJUSTL(tchar1))
      !finding uncommented line that isn't empty
      IF(tchar1(1:1) .NE. '!' .AND. tchar1 .NE. '')THEN
        !ignore commented portions of line
        CALL parse(tchar1,'!',words,nwords)
        tchar1=TRIM(ADJUSTL(words(1)))
        CALL parse(tchar1,' ',words,nwords)
        EXIT
      ENDIF
    ENDDO
  ENDSUBROUTINE get_thor_line

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE find_line(linechar)
    CHARACTER(*),INTENT(IN) :: linechar
    CHARACTER(80) :: tchar1

    DO
      READ(22,*)tchar1
      IF(TRIM(ADJUSTL(tchar1)) .EQ. TRIM(ADJUSTL(linechar)))EXIT
    ENDDO
    BACKSPACE(22)
  ENDSUBROUTINE find_line
END MODULE infuncs
