MODULE check_input
  !***********************************************************************
  ! This module contains variables and subroutines that enable writing
  ! a *.vtk file that can be used in VisIt to visualize mesh and geometry.
  !***********************************************************************

  ! User derived-type modules

  USE types
  USE parameter_types
  USE filename_types
  USE vector_types
  USE cross_section_types
  USE geometry_types
  USE angle_types
  USE multindex_types
  USE global_variables

  IMPLICIT NONE

CONTAINS

  SUBROUTINE plot_source
    !*********************************************************************
    !
    ! Subroutine prints a *.vtk file containing source information
    !
    !*********************************************************************

    ! local variables

    INTEGER(kind=li) :: i,l,eg,src_indx

    OPEN(unit=10,file=TRIM(vtk_src_filename),status='unknown',action='write')

    WRITE(10,'(a26)') '# vtk DataFile Version 3.0'
    WRITE(10,'(a72)') jobname
    WRITE(10,'(a5)') 'ASCII'
    WRITE(10,'(a25)') 'DATASET UNSTRUCTURED_GRID'
    WRITE(10,'(a6,1x,i12,1x,a5)') 'POINTS',num_vert,'float'

    DO i=1, num_vert
      WRITE(10,'(3(1x,es12.5))') vertices(i)%v%x1,vertices(i)%v%x2,&
            vertices(i)%v%x3
    END DO

    WRITE(10,*)
    WRITE(10,'(a5,1x,i12,1x,i12)') 'CELLS',num_cells,num_cells+&
          4*num_cells

    DO i=1, num_cells
      WRITE(10,'(5(i12,1x))') 4,cells(i)%R(0)-1,cells(i)%R(1)-1,&
            cells(i)%R(2)-1,cells(i)%R(3)-1
    END DO

    WRITE(10,*)
    WRITE(10,'(a11,1x,i12)') 'CELL_TYPES',num_cells

    DO i=1, num_cells
      WRITE(10,'(i12)') 10
    END DO

    WRITE(10,*)
    WRITE(10,'(a9,1x,i12)') 'CELL_DATA', num_cells
    WRITE(10,'(a5,1x,a16,1x,i12)') 'FIELD','Neutronics_Edits',egmax

    l=1
    DO eg=1, egmax
      WRITE(10,'(i12,1x,i12,1x,i12,1x,a5)') eg,1,num_cells,'float'
      DO i=1, num_cells
        src_indx=source_ids(cells(i)%src)
        WRITE(10,'(es12.5)')ext_src(src_indx)%mom(1,1,eg)
      END DO
    END DO

    CLOSE(10)

  END SUBROUTINE plot_source

  SUBROUTINE plot_region
    !*********************************************************************
    !
    ! Subroutine prints a *.vtk file containing cell material (and region) IDs
    !
    !*********************************************************************

    ! local variables

    INTEGER(kind=li) :: i,l

    OPEN(unit=10,file=vtk_reg_filename,status='unknown',action='write')

    WRITE(10,'(a26)') '# vtk DataFile Version 3.0'
    WRITE(10,'(a72)') jobname
    WRITE(10,'(a5)') 'ASCII'
    WRITE(10,'(a25)') 'DATASET UNSTRUCTURED_GRID'
    WRITE(10,'(a6,1x,i12,1x,a5)') 'POINTS',num_vert,'float'

    DO i=1, num_vert
      WRITE(10,'(3(1x,es12.5))') vertices(i)%v%x1,vertices(i)%v%x2,&
            vertices(i)%v%x3
    END DO

    WRITE(10,*)
    WRITE(10,'(a5,1x,i12,1x,i12)') 'CELLS',num_cells,num_cells+&
          4*num_cells

    DO i=1, num_cells
      WRITE(10,'(5(i12,1x))') 4,cells(i)%R(0)-1,cells(i)%R(1)-1,&
            cells(i)%R(2)-1,cells(i)%R(3)-1
    END DO

    WRITE(10,*)
    WRITE(10,'(a11,1x,i12)') 'CELL_TYPES',num_cells

    DO i=1, num_cells
      WRITE(10,'(i12)') 10
    END DO

    WRITE(10,*)
    WRITE(10,'(a9,1x,i12)') 'CELL_DATA', num_cells
    WRITE(10,'(a5,1x,a16,1x,i12)') 'FIELD','Neutronics_Edits',1

    l=1
    WRITE(10,'(i12,1x,i12,1x,i12,1x,a5)') 1,1,num_cells,'float'
    DO i=1, num_cells
      WRITE(10,'(es12.3)') REAL(cells(i)%reg,d_t)
    END DO

    CLOSE(10)

  END SUBROUTINE plot_region

  SUBROUTINE plot_material
    !*********************************************************************
    !
    ! Subroutine prints a *.vtk file containing cell material (and region) IDs
    !
    !*********************************************************************

    ! local variables

    INTEGER(kind=li) :: i,l

    OPEN(unit=10,file=vtk_mat_filename,status='unknown',action='write')

    WRITE(10,'(a26)') '# vtk DataFile Version 3.0'
    WRITE(10,'(a72)') jobname
    WRITE(10,'(a5)') 'ASCII'
    WRITE(10,'(a25)') 'DATASET UNSTRUCTURED_GRID'
    WRITE(10,'(a6,1x,i12,1x,a5)') 'POINTS',num_vert,'float'

    DO i=1, num_vert
      WRITE(10,'(3(1x,es12.5))') vertices(i)%v%x1,vertices(i)%v%x2,&
            vertices(i)%v%x3
    END DO

    WRITE(10,*)
    WRITE(10,'(a5,1x,i12,1x,i12)') 'CELLS',num_cells,num_cells+&
          4*num_cells

    DO i=1, num_cells
      WRITE(10,'(5(i12,1x))') 4,cells(i)%R(0)-1,cells(i)%R(1)-1,&
            cells(i)%R(2)-1,cells(i)%R(3)-1
    END DO

    WRITE(10,*)
    WRITE(10,'(a11,1x,i12)') 'CELL_TYPES',num_cells

    DO i=1, num_cells
      WRITE(10,'(i12)') 10
    END DO

    WRITE(10,*)
    WRITE(10,'(a9,1x,i12)') 'CELL_DATA', num_cells
    WRITE(10,'(a5,1x,a16,1x,i12)') 'FIELD','Neutronics_Edits',1

    l=1
    WRITE(10,'(i12,1x,i12,1x,i12,1x,a5)') 1,1,num_cells,'float'
    DO i=1, num_cells
      WRITE(10,'(es12.3)') REAL(reg2mat(cells(i)%reg),d_t)
    END DO

    CLOSE(10)

  END SUBROUTINE plot_material

  SUBROUTINE print_xs
    !*********************************************************************
    !
    ! Subroutine prints cross section sets
    !
    !*********************************************************************

    ! Local variables

    INTEGER :: g,gp,m,order
    REAL(kind=d_t) :: sigs(egmax)

    ! Print cross sections

    WRITE(6,*)
    WRITE(6,*) '------------------------------------------------------------------'
    WRITE(6,*) '--------------------- Echoing Cross Sections ---------------------'
    WRITE(6,*) '------------------------------------------------------------------'
    WRITE(6,*)
    WRITE(6,101) 'Number of materials:             ',num_mat
    WRITE(6,101) 'Number of groups:                ',egmax
    WRITE(6,101) 'Scattering expansion order read: ',xs_ord
    IF ( most_thermal>egmax ) THEN
      WRITE(6,*) 'No upscattering present.'
    ELSE
      WRITE(6,101) 'Most thermal group:              ',most_thermal
    END IF
    DO m=1,num_mat
      WRITE(6,102) 'Material ',TRIM(xs_mat(m)%mat_name),' with ID ',xs_mat(m)%mat_id
      WRITE(6,103) 'Group','SigT','SigF','SigS','nu*SigF','chi'
      DO g=1,egmax
        sigs=0.0_d_t
        DO gp=1,egmax
          sigs(g)=sigs(g)+xs_mat(m)%sigma_scat(1,gp,g)
        END DO
        WRITE(6,104) g,xs_mat(m)%sigma_t(g),xs_mat(m)%sigma_f(g), &
              sigs(g),xs_mat(m)%sigma_f(g)*xs_mat(m)%nu(g),&
              xs_mat(m)%chi(g)
      END DO
      WRITE(6,*) 'Scattering Matrix, from -> columns, to -> row'
      DO order=1, xs_ord+1
        WRITE(6,101) 'Scattering order: ',order-1
        DO g=1,egmax
          WRITE(6,105) (xs_mat(m)%sigma_scat(order,g,gp),gp=1,egmax)
        END DO
      END DO
    END DO
    WRITE(6,*)
101 FORMAT(1X,A,I8)
102 FORMAT(1X,A,A,A,I8)
103 FORMAT(1X,A9,5A15)
104 FORMAT(1X,I9,5ES15.4)
105 FORMAT(1X,12ES15.4)
  END SUBROUTINE print_xs

END MODULE check_input
