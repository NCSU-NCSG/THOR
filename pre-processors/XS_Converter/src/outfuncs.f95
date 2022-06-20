!output functions
MODULE outfuncs
  USE globals
  USE HDF5
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: outputxs
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !based on xs output format, call the xs out-putter
  SUBROUTINE outputxs()
    SELECTCASE(outformat)
      CASE('thor')
        CALL out_thor()
      CASE('mcnp')
        CALL out_mcnp()
      CASE('openmc')
        CALL out_openmc()
      CASE DEFAULT
        STOP 'bad output format'
    ENDSELECT
  ENDSUBROUTINE outputxs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE out_thor()
    INTEGER :: ios,m,gp,l
    CHARACTER(64) :: tchar1

    !open xsout file
    OPEN(UNIT=32,FILE=xsout,STATUS='REPLACE',ACTION='WRITE',IOSTAT=ios,IOMSG=tchar1)
    IF(ios .NE. 0)THEN
        WRITE(*,'(A)')tchar1
        STOP
    ENDIF

    !output the xs characteristics
    WRITE(32,'(A,I0,A,I0,A,I0)')'THOR_XS_V1 ',nummats,' ',numgroups,' ',levelanis
    !output the energy group structure
    WRITE(tchar1,'(10000ES20.12)')eg_struc(1:numgroups)
    WRITE(32,'(A)')TRIM(ADJUSTL(tchar1))
    !output the xs data
    DO m=1,nummats
      WRITE(32,'(A,I0)')'id ',m
      WRITE(tchar1,'(10000ES20.12)')chi(m,:)
      WRITE(32,'(A)')TRIM(ADJUSTL(tchar1))
      WRITE(tchar1,'(10000ES20.12)')sigmaf(m,:)
      WRITE(32,'(A)')TRIM(ADJUSTL(tchar1))
      WRITE(tchar1,'(10000ES20.12)')nuf(m,:)
      WRITE(32,'(A)')TRIM(ADJUSTL(tchar1))
      WRITE(tchar1,'(10000ES20.12)')sigmat(m,:)
      WRITE(32,'(A)')TRIM(ADJUSTL(tchar1))
      DO l=1,levelanis+1
        DO gp=1,numgroups
          WRITE(32,'(10000ES20.12)')sigmas(m,l,gp,:)
        ENDDO
      ENDDO
    ENDDO
    !close the output file
    CLOSE(32)
  ENDSUBROUTINE out_thor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE out_mcnp()
    STOP 'out_mcnp not yet complete'
  ENDSUBROUTINE out_mcnp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE out_openmc()
    INTEGER(4) :: err1
    INTEGER(8) :: h5fileid
    xsout=TRIM(ADJUSTL(xsout))//'.hdf5'

    h5fileid=32

    !initialize hdf5
    CALL h5open_f(err1)
    IF(err1 .NE. 0)STOP 'error opening output file'

    !create h5 file
    CALL h5fcreate_f(xsout,H5F_ACC_TRUNC_F,h5fileid,err1)
    IF(err1 .NE. 0)STOP 'error opening output file'

    STOP 'out_openmc not yet complete'
  ENDSUBROUTINE out_openmc
END MODULE outfuncs
