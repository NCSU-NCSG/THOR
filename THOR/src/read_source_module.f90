module read_source_module
!***********************************************************************
!
! Read source module contains all subroutines needed to read source file
!
!***********************************************************************
  use types
  use parameter_types
  use filename_types
  use multindex_types
  use global_variables
  use termination_module

  implicit none

contains

  subroutine read_src
  !*********************************************************************
  !
  ! Subroutine reads source in 'unique' ahot format
  ! (could be adapted for other formats, of course)
  !
  !*********************************************************************

  ! Declare temporary variables

    integer(kind=li) :: alloc_stat, eg, m, l

  ! Open and read source file 

    open(unit=10,file=trim(source_filename),status='unknown',action='read')

  ! Read source strength and moments from file

    read(10,*) num_src_mat

    allocate(src_mat(num_src_mat),stat=alloc_stat)
    if(alloc_stat /= 0) call stop_thor(2_li)

    allocate(src_str(0:num_src_mat-1,egmax),stat=alloc_stat)
    if(alloc_stat /= 0) call stop_thor(2_li)

    allocate(src_m(num_moments_v,0:num_src_mat-1,egmax),stat=alloc_stat)
    if(alloc_stat /= 0) call stop_thor(2_li)

    do m=1, num_src_mat
       read(10,*) src_mat(m)
       do eg=1, egmax
          read(10,*) src_str(src_mat(m),eg)
          do l=1, num_moments_v
             read(10,*) src_m(l,src_mat(m),eg)  
          end do
       end do
    end do

  ! Close mesh file

    close(10)
 
  end subroutine read_src

end module read_source_module

