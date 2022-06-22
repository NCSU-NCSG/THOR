!output functions
MODULE outfuncs
  USE globals
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
      WRITE(tchar1,'(10000ES16.8)')chi(m,:)
      WRITE(32,'(A)')TRIM(ADJUSTL(tchar1))
      WRITE(tchar1,'(10000ES16.8)')sigmaf(m,:)
      WRITE(32,'(A)')TRIM(ADJUSTL(tchar1))
      WRITE(tchar1,'(10000ES16.8)')nuf(m,:)
      WRITE(32,'(A)')TRIM(ADJUSTL(tchar1))
      WRITE(tchar1,'(10000ES16.8)')sigmat(m,:)
      WRITE(32,'(A)')TRIM(ADJUSTL(tchar1))
      DO l=1,levelanis+1
        DO gp=1,numgroups
          WRITE(32,'(10000ES16.8)')sigmas(m,l,gp,:)
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
    INTEGER :: ios,g,m,gp,l
    CHARACTER(64) :: tchar1
    xsout=TRIM(xsout)//'.py'
    !open xsout file
    OPEN(UNIT=32,FILE=xsout,STATUS='REPLACE',ACTION='WRITE',IOSTAT=ios,IOMSG=tchar1)
    IF(ios .NE. 0)THEN
        WRITE(*,'(A)')tchar1
        STOP
    ENDIF
    WRITE(32,'(A)')'import numpy as np'
    WRITE(32,'(A)')'import openmc'
    WRITE(32,'(A)')
    WRITE(32,'(A)')'#INSTRUCTIONS: place these cross sections at the top of your OpenMC python script'
    WRITE(32,'(A)')'#OR build your script from this baseline'
    WRITE(32,'(A)')'#If you wish to reduce the number of angular moments, do so manually in the script'
    WRITE(32,'(A)')'#OR alter the XS input file'
    WRITE(32,'(A)')'####################################################################################################'
    WRITE(32,'(A)')'###################################-Beginning of Cross Sections-####################################'
    WRITE(32,'(A)')'####################################################################################################'
    WRITE(32,'(A)')'#Group structure:'
    WRITE(32,'(A,ES15.8)',ADVANCE='NO')'groups = openmc.mgxs.EnergyGroups(',0.0D0
    DO g=numgroups,1,-1
      IF(MOD(g,6) .EQ. 0)WRITE(32,'(A)')'                                  '
      WRITE(32,'(A,ES15.8)',ADVANCE='NO')',',eg_struc(g)
    ENDDO
    WRITE(32,'(A)')')'
    DO m=1,nummats
      WRITE(32,'(A,I0,A)')'#Data for Material ',m,':'
      WRITE(32,'(A,I0,A,I0,A)')"mat",m,"_xsdat = openmc.XSdata('mat_",m,"', groups)"
      WRITE(32,'(A,I0,A,I0)')"mat",m,"_xsdat.order = ",levelanis
      !total xs
      CALL print_xs_openmc(m,'total',sigmat(m,:))
      CALL print_xs_openmc(m,'absorption',sigmaa(m,:))
      CALL print_xs_openmc(m,'fission',sigmaf(m,:))
      CALL print_xs_openmc(m,'nu_fission',nuf(m,:)*sigmaf(m,:))
      CALL print_xs_openmc(m,'chi',chi(m,:))
      !print the scattering matrix, a bit more involved...
      WRITE(32,'(A)')'scatter_matrix = np.array(\'
      WRITE(32,'(A)',ADVANCE='NO')'    ['
      DO l=1,levelanis+1
        IF(l .NE. 1)WRITE(32,'(A)',ADVANCE='NO')'     '
        WRITE(32,'(A)',ADVANCE='NO')'['
        DO gp=1,numgroups
          IF(gp .NE. 1)WRITE(32,'(A)',ADVANCE='NO')'      '
          WRITE(32,'(A)',ADVANCE='NO')'['
          DO g=1,numgroups
            IF(g .NE. 1)WRITE(32,'(A)',ADVANCE='NO')','
            WRITE(32,'(ES15.8)',ADVANCE='NO')sigmas(m,l,gp,g)
          ENDDO
          WRITE(32,'(A)',ADVANCE='NO')']'
          IF(gp .NE. numgroups)WRITE(32,'(A)')','
        ENDDO
        WRITE(32,'(A)',ADVANCE='NO')']'
        IF(l .NE. levelanis+1)WRITE(32,'(A)')','
      ENDDO
      WRITE(32,'(A)')'])'
      WRITE(32,'(A)')'scatter_matrix = np.transpose(scatter_matrix)'
      WRITE(32,'(A,I0,A)')"mat",m,"_xsdat.set_scatter_matrix(scatter_matrix, temperature=294.)"
    ENDDO
    WRITE(32,'(A)')'#Create the cross sections hdf5 file:'
    WRITE(32,'(A)')'mg_cross_sections_file = openmc.MGXSLibrary(groups)'
    DO m=1,nummats
      WRITE(32,'(A,I0,A)')'mg_cross_sections_file.add_xsdata(mat',m,'_xsdat)'
    ENDDO
    WRITE(32,'(A)')"mg_cross_sections_file.export_to_hdf5('cross_sections.h5')"
    WRITE(32,'(A)')'#Assign each cross section to a separate material:'
    WRITE(32,'(A)')'materials = {}'
    WRITE(32,'(A)',ADVANCE='NO')"for xs in ["
    DO m=1,nummats
      IF(m .NE. 1)WRITE(32,'(A)',ADVANCE='NO')','
      WRITE(32,'(A,I0,A)',ADVANCE='NO')"'mat_",m,"'"
    ENDDO
    WRITE(32,'(A)')']:'
    WRITE(32,'(A)')"    materials[xs] = openmc.Material(name=xs)"
    WRITE(32,'(A)')"    materials[xs].set_density('macro', 1.)"
    WRITE(32,'(A)')"    materials[xs].add_macroscopic(xs)"
    WRITE(32,'(A)')'#Create the materials file for this specification:'
    WRITE(32,'(A)')"materials_file = openmc.Materials(materials.values())"
    WRITE(32,'(A)')"materials_file.cross_sections = 'cross_sections.h5'"
    WRITE(32,'(A)')"materials_file.export_to_xml()"
    WRITE(32,'(A)')'####################################################################################################'
    WRITE(32,'(A)')'######################################-End of Cross Sections-#######################################'
    WRITE(32,'(A)')'####################################################################################################'
    WRITE(32,'(A)')"#The following creates the settings and sets the energy mode to multi-group:"
    WRITE(32,'(A)')"settings_file = openmc.Settings()"
    WRITE(32,'(A)')"settings_file.energy_mode = 'multi-group'"
    WRITE(32,'(A)')"#The user should complete the OpenMC input below by specifying geometry, settings, etc."
  ENDSUBROUTINE out_openmc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE print_xs_openmc(matnum,xstype,xsarr)
    INTEGER,INTENT(IN) :: matnum
    CHARACTER(*),INTENT(IN) :: xstype
    REAL(8),INTENT(IN) :: xsarr(*)
    INTEGER :: g

    WRITE(32,'(A,I0,3A,ES15.8)',ADVANCE='NO')'mat',matnum,'_xsdat.set_',TRIM(ADJUSTL(xstype)),'([',xsarr(1)
    DO g=2,numgroups
      IF(MOD(g,6) .EQ. 0)WRITE(32,'(A)')'                                  '
      WRITE(32,'(A,ES15.8)',ADVANCE='NO')',',xsarr(g)
    ENDDO
    WRITE(32,'(A)')'], temperature=294.)'
  ENDSUBROUTINE print_xs_openmc
END MODULE outfuncs
