module check_input
!***********************************************************************
! This module contains variables and subroutines that enable writing 
! a *.vtk file that can be used in VisIt to visualize mesh and geometry.
!***********************************************************************

! User derived-type modules

  use types
  use parameter_types
  use filename_types
  use vector_types
  use cross_section_types
  use geometry_types
  use angle_types
  use multindex_types
  use global_variables

  implicit none

  contains

  subroutine plot_source
  !*********************************************************************
  !
  ! Subroutine prints a *.vtk file containing source information
  !
  !*********************************************************************

       ! local variables  

       integer(kind=li) :: i,l,eg

       open(unit=10,file=trim(vtk_src_filename),status='unknown',action='write')

       write(10,'(a26)') '# vtk DataFile Version 3.0' 
       write(10,'(a72)') jobname
       write(10,'(a5)') 'ASCII'
       write(10,'(a25)') 'DATASET UNSTRUCTURED_GRID'
       write(10,'(a6,1x,i12,1x,a5)') 'POINTS',num_vert,'float'
       
       do i=1, num_vert
          write(10,'(3(1x,es12.5))') vertices(i)%v%x1,vertices(i)%v%x2,&
               vertices(i)%v%x3
       end do
       
       write(10,*)
       write(10,'(a5,1x,i12,1x,i12)') 'CELLS',num_cells,num_cells+&
            4*num_cells
       
       do i=1, num_cells
          write(10,'(5(i12,1x))') 4,cells(i)%R(0)-1,cells(i)%R(1)-1,&
               cells(i)%R(2)-1,cells(i)%R(3)-1
       end do
       
       write(10,*)
       write(10,'(a11,1x,i12)') 'CELL_TYPES',num_cells
       
       do i=1, num_cells
          write(10,'(i12)') 10
       end do
       
       write(10,*)
       write(10,'(a9,1x,i12)') 'CELL_DATA', num_cells
       write(10,'(a5,1x,a16,1x,i12)') 'FIELD','Neutronics_Edits',egmax
       
       l=1
       do eg=1, egmax
          write(10,'(i12,1x,i12,1x,i12,1x,a5)') eg,1,num_cells,'float'
          do i=1, num_cells     
             write(10,'(es12.5)') &
                  src_str(cells(i)%src,eg)*src_m(1,cells(i)%src,eg)  
          end do
       end do

       close(10)

  end subroutine
  
  subroutine plot_region
  !*********************************************************************
  !
  ! Subroutine prints a *.vtk file containing cell material (and region) IDs
  !
  !*********************************************************************

  ! local variables  

    integer(kind=li) :: i,l

    open(unit=10,file=vtk_reg_filename,status='unknown',action='write')

    write(10,'(a26)') '# vtk DataFile Version 3.0' 
    write(10,'(a72)') jobname
    write(10,'(a5)') 'ASCII'
    write(10,'(a25)') 'DATASET UNSTRUCTURED_GRID'
    write(10,'(a6,1x,i12,1x,a5)') 'POINTS',num_vert,'float'
       
    do i=1, num_vert
      write(10,'(3(1x,es12.5))') vertices(i)%v%x1,vertices(i)%v%x2,&
            vertices(i)%v%x3
    end do
       
    write(10,*)
    write(10,'(a5,1x,i12,1x,i12)') 'CELLS',num_cells,num_cells+&
          4*num_cells
       
    do i=1, num_cells
      write(10,'(5(i12,1x))') 4,cells(i)%R(0)-1,cells(i)%R(1)-1,&
            cells(i)%R(2)-1,cells(i)%R(3)-1
    end do
       
    write(10,*)
    write(10,'(a11,1x,i12)') 'CELL_TYPES',num_cells
       
    do i=1, num_cells
      write(10,'(i12)') 10
    end do
       
    write(10,*)
    write(10,'(a9,1x,i12)') 'CELL_DATA', num_cells
    write(10,'(a5,1x,a16,1x,i12)') 'FIELD','Neutronics_Edits',1
       
    l=1
    write(10,'(i12,1x,i12,1x,i12,1x,a5)') 1,1,num_cells,'float'
    do i=1, num_cells     
      write(10,'(es12.3)') real(cells(i)%reg,d_t)
    end do

    close(10)

  end subroutine

  subroutine plot_material
  !*********************************************************************
  !
  ! Subroutine prints a *.vtk file containing cell material (and region) IDs
  !
  !*********************************************************************

  ! local variables  

    integer(kind=li) :: i,l

    open(unit=10,file=vtk_mat_filename,status='unknown',action='write')

    write(10,'(a26)') '# vtk DataFile Version 3.0' 
    write(10,'(a72)') jobname
    write(10,'(a5)') 'ASCII'
    write(10,'(a25)') 'DATASET UNSTRUCTURED_GRID'
    write(10,'(a6,1x,i12,1x,a5)') 'POINTS',num_vert,'float'
       
    do i=1, num_vert
      write(10,'(3(1x,es12.5))') vertices(i)%v%x1,vertices(i)%v%x2,&
            vertices(i)%v%x3
    end do
       
    write(10,*)
    write(10,'(a5,1x,i12,1x,i12)') 'CELLS',num_cells,num_cells+&
          4*num_cells
       
    do i=1, num_cells
      write(10,'(5(i12,1x))') 4,cells(i)%R(0)-1,cells(i)%R(1)-1,&
            cells(i)%R(2)-1,cells(i)%R(3)-1
    end do
       
    write(10,*)
    write(10,'(a11,1x,i12)') 'CELL_TYPES',num_cells
       
    do i=1, num_cells
      write(10,'(i12)') 10
    end do
       
    write(10,*)
    write(10,'(a9,1x,i12)') 'CELL_DATA', num_cells
    write(10,'(a5,1x,a16,1x,i12)') 'FIELD','Neutronics_Edits',1
       
    l=1
    write(10,'(i12,1x,i12,1x,i12,1x,a5)') 1,1,num_cells,'float'
    do i=1, num_cells     
     write(10,'(es12.3)') real(reg2mat(cells(i)%reg),d_t)
    end do

    close(10)

  end subroutine

  subroutine print_xs
  !*********************************************************************
  !
  ! Subroutine prints cross section sets
  !
  !*********************************************************************

  ! Local variables
    
    integer :: g,gp,m,l,order 
    real(kind=d_t) :: sigs(egmax)
 
  ! Print cross sections

    write(6,*) 
    write(6,*) '------------------------------------------------------------------'
    write(6,*) '--------------------- Echoing Cross Sections ---------------------'  
    write(6,*) '------------------------------------------------------------------'
    write(6,*)
    write(6,101) 'Number of materials:             ',num_mat
    write(6,101) 'Number of groups:                ',egmax
    write(6,101) 'Scattering expansion order read: ',xs_ord
    if ( most_thermal>egmax ) then
      write(6,*) 'No upscattering present.' 
    else
      write(6,101) 'Most thermal group:              ',most_thermal
    end if
    do m=1,num_mat
       write(6,102) 'Material ',m,' with ID ',xs_mat(m)%mat
       write(6,103) 'Group','SigT','SigF','SigS','nu*SigF','chi'
       do g=1,egmax
         sigs=0.0_d_t
         do gp=1,egmax
            sigs(g)=sigs(g)+sigma_scat(xs_mat(m)%mat,1,gp,g)%xs
         end do
         write(6,104) g,sigma_t(xs_mat(m)%mat,g)%xs,fiss(xs_mat(m)%mat,g)%xs, &
                      sigs(g),fiss(xs_mat(m)%mat,g)%xs*nu(xs_mat(m)%mat,g)%xs,&
                      chi(xs_mat(m)%mat,g)%xs
       end do
       write(6,*) 'Scattering Matrix, from -> columns, to -> row'
       do order=1, xs_ord+1
         write(6,101) 'Scattering order: ',order 
         do g=1,egmax
            write(6,105) (sigma_scat(xs_mat(m)%mat,order,g,gp)%xs,gp=1,egmax)
         end do
       end do
    end do
    write(6,*) 
    101 FORMAT(1X,A,I8)
    102 FORMAT(1X,A,I8,A,I8)
    103 FORMAT(1X,A9,5A15)
    104 FORMAT(1X,I9,5ES15.4)
    105 FORMAT(1X,12ES15.4)
  end subroutine  

end module
