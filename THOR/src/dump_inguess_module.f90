module dump_inguess_module
!***********************************************************************
!
! This module contains subroutines to create and read restart(dump) files
!
!***********************************************************************

! Use derived-type modules

  use types
  use parameter_types
  use filename_types
  use vector_types
  use cross_section_types
  use geometry_types
  use angle_types
  use multindex_types
  use global_variables
  use termination_module

  implicit none

contains


  subroutine inguess_eig(flux,keff,flag)
  !*********************************************************************
  !
  ! Subroutine inguess_eig reads a file that is created by either 
  ! dump_PI or dump_jfnk and uses it as initial guess
  !
  !*********************************************************************
 
  ! Pass arguments

    integer(kind=li) :: flag
    real(kind=d_t) :: keff
    real(kind=d_t) :: flux(num_moments_v,namom,num_cells,egmax,niter)

  ! local variables

    integer(kind=li) :: eg,i,l,p,q
    integer(kind=li) :: egmax_t,num_cells_t,namom_t,num_moments_v_t
    logical :: existence

  ! open file 
    inquire(file=trim(inguess_file), exist=existence)
    if (existence .eqv. .TRUE.) THEN
       open (UNIT = 66, FILE = inguess_file, STATUS = "OLD", ACTION = "READ",&
             form='unformatted')
      ! read inguess file
       read(66) egmax_t,num_cells_t,namom_t,num_moments_v_t
      ! check if data in file and data from stdin are identical
       if(egmax_t .ne. egmax .or. num_cells_t .ne. num_cells .or.    &
          num_moments_v_t .ne. num_moments_v .or. namom_t .ne. namom ) then
          call stop_thor(24_li)
       end if
       do eg =1,egmax
         do i=1,num_cells
           do p=1,namom
             do l=1,num_moments_v 
                read(66) flux(l,p,i,eg,flag) 
             end do
           end do
         end do
       end do
       read(66) keff
      ! close file

        close(unit=66) 
        write(6,*)
        write(6,*) 'Initial guess file ',trim(inguess_file),' was successfully read'  
 
    else
       write(6,*)
       write(6,*) "Attempted to read initial guess file ",trim(inguess_file)," but could not find it."
       write(6,*)
       flux(:,:,:,:,niter)=zero
       do eg =1,egmax
         do i=1,num_cells
           flux(1,1,i,eg,niter)=one
         end do
       end do
       keff = 1.0_d_t
    end if

  end subroutine

  subroutine dump_jfnk(flux,keff)
  !*********************************************************************
  !
  ! Subroutine dump_jfnk creates restart file after each Newton iteration
  !
  !*********************************************************************
  
  ! Pass arguments
    
    real(kind=d_t) :: keff
    real(kind=d_t) :: flux(num_moments_v,namom,num_cells,egmax,niter)
 
  ! local variables

    integer(kind=li) :: eg,i,l,p,q
    logical :: existence

  ! open file 
    inquire(file=trim(dump_file), exist=existence)
    if (existence .eqv. .TRUE.) THEN
       open (UNIT = 56, FILE = dump_file, STATUS = "OLD", ACTION = "WRITE",&
             form='unformatted')
    else
       open (UNIT = 56, FILE = dump_file, STATUS = "NEW", ACTION ="WRITE",& 
       form='unformatted')
    end if

  ! write dump file
      write(56) egmax,num_cells,namom,num_moments_v
       do eg =1,egmax
         do i=1,num_cells
           do p=1,namom
             do l=1,num_moments_v
                write(56) flux(l,p,i,eg,1)
             end do
           end do
         end do
       end do
      write(56) keff

  ! close file
    
    close(unit=56)
 
  end subroutine 

  subroutine dump_PI(flux,keff)
  !*********************************************************************
  !
  ! Subroutine dump_PI creates restart filse after each power iteration
  !
  !*********************************************************************

  ! Pass arguments

    real(kind=d_t) :: keff
    real(kind=d_t) :: flux(num_moments_v,namom,num_cells,egmax,niter)

  ! local variables

    integer(kind=li) :: eg,i,l,p,q
    logical :: existence

  ! open file 
    inquire(file=trim(dump_file), exist=existence)
    if (existence .eqv. .TRUE.) THEN
       open (UNIT = 56, FILE = dump_file, STATUS = "OLD", ACTION = "WRITE",&
             form='unformatted')
    else
       open (UNIT = 56, FILE = dump_file, STATUS = "NEW", ACTION ="WRITE",&
       form='unformatted')
    end if

  ! write dump file
      write(56) egmax,num_cells,namom,num_moments_v
      do eg =1,egmax
        do i=1,num_cells
          do p=1,namom
            do l=1,num_moments_v
               write(56) flux(l,p,i,eg,niter) 
            end do
          end do   
        end do
      end do
      write(56) keff
   
  ! close file
    
    close(unit=56)

  end subroutine
 
end module 
