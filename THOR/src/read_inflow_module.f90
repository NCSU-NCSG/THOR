module read_inflow_module
!***********************************************************************
!
! Read source module contains all subroutines needed to read boundary 
! source file
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

  subroutine read_finflow
  !*********************************************************************
  !
  ! Subroutine reads source in 'unique' ahot format
  ! (could be adapted for other formats, of course)
  !
  !*********************************************************************

  ! Declare temporary variables

    integer(kind=li) :: alloc_stat, eg, q, octant, m, i, cell, f, face
    real(kind=d_t)   :: dmy

  ! Depending on page_iflw set eg_iflw

    if(page_iflw .eq. 1) then
       eg_iflw=1_li
    else 
       eg_iflw=egmax
    end if 

  ! allocate binflx: indices followig column major. Fastest running index
  ! is fixed boundary faces then octants then angles, then groups

    allocate(binflx(num_moments_f,fside_cells,8,nangle,eg_iflw),stat=alloc_stat)
    if(alloc_stat /= 0) call stop_thor(2_li)

  ! if page_iflw == 0 then read binflow in full, otherwise just open a file

    if(page_iflw.eq.0) then
  
    ! Open and read source file 

      open(unit=10,file=trim(finflow_filename),status='unknown',action='read')

    ! Read source strength and moments from file

      do eg=1,egmax
        do q=1,nangle
          do octant=1,8
            do f=1,fside_cells 
               read(10,*) face    ! this is the face number in the b_cells array
               if( b_cells(face)%bc .ne. 2 ) then
                  call stop_thor(16_li) 
               end if
               face = b_cells(face)%ptr   ! make sure that bc are stored in the same order as in fb_cells
               read(10,*) (binflx(m,face,octant,q,eg),m=1,num_moments_f)            
            end do     
          end do
        end do
      end do

    ! Close inflow file
 
      close(10)

    else

      open(unit=97,file=trim(finflow_filename),status='unknown',action='read')
      do eg=1,egmax
        do q=1,nangle
          do octant=1,8
            do f=1,fside_cells
               read(97,*) face    ! this is the face number in the b_cells array
               if( b_cells(face)%bc .ne. 2 ) then
                  call stop_thor(16_li)
               end if
               read(97,*) (dmy,m=1,num_moments_f)
            end do
          end do
        end do
      end do
      rewind(unit=97)
   
    end if

  end subroutine read_finflow

end module read_inflow_module

