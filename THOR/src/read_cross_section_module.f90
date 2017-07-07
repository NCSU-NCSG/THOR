module read_cross_section_module
!***********************************************************************
!
! Read cross-section module  contains all subroutines needed to read
! cross-section file
!
!***********************************************************************

  use types
  use parameter_types
  use filename_types
  use cross_section_types
  use global_variables
  use termination_module

  implicit none

contains

  subroutine read_xs
  !*********************************************************************
  !
  ! Subroutine reads cross-section in 'unique' ahot format
  ! (could be adapted for other formats, of course)
  !
  !*********************************************************************

    integer(kind=li) :: alloc_stat, e1, order, eg_to, eg_from,l,m

  ! Open and read mesh file 

    open(unit=11,file=trim(cross_section_filename),status='old',action='read')

    read(11,*) num_mat    

  ! Allocate cross-section arrays and check if enough memory is available

    allocate(xs_mat(0:num_mat),chi(0:num_mat,egmax),&
         eg_bounds(0:num_mat,egmax+1),fiss(0:num_mat,egmax),&
         nu(0:num_mat,egmax),sigma_t(0:num_mat,egmax),tsigs(0:num_mat,egmax),&
         sigma_scat(0:num_mat,xs_ord+1,egmax,egmax),&
         stat=alloc_stat)
    if(alloc_stat /= 0) call stop_thor(2_li) 

  ! Set most_thermal to no upscattering

    most_thermal=egmax+1

  ! Read material total & scattering cross-section

    do m=1, num_mat
       read(11,*) xs_mat(m)%mat
       if(multiplying .ne. 0)then
          read(11,*) (chi(xs_mat(m)%mat,e1)%xs,e1=1,egmax)
          read(11,*) (eg_bounds(xs_mat(m)%mat,e1)%xs,e1=1,egmax)
          eg_bounds(xs_mat(m)%mat,egmax+1)%xs=0.0_d_t
          read(11,*) (fiss(xs_mat(m)%mat,e1)%xs,e1=1,egmax)
          read(11,*) (nu(xs_mat(m)%mat,e1)%xs,e1=1,egmax)
       end if
       read(11,*) (sigma_t(xs_mat(m)%mat,e1)%xs,e1=1,egmax)

    ! Initialize scattering matrix
   
       do order=1, xs_ord+1
          do eg_to=1,egmax
             do eg_from=1,egmax
               sigma_scat(xs_mat(m)%mat,order,eg_to,eg_from)%xs=zero 
             end do
          end do
       end do

    ! Read scattering matrix, note: thermal groups are separated from fast
    ! groups but old cross section format remains valid  

       if(upscattering.eq.0) then
          do order=1, xs_ord+1
             do eg_to=1,egmax
                read(11,*) (sigma_scat(xs_mat(m)%mat,order,&
                  eg_to,eg_from)%xs,eg_from=1,eg_to)
             end do
          end do
       else
          do order=1, xs_ord+1
            do eg_to=1,egmax
               read(11,*) (sigma_scat(xs_mat(m)%mat,order,&
                    eg_to,eg_from)%xs,eg_from=1,egmax)
            end do 
          end do
       end if 

    ! Determine most_thermal
 
      if(upscattering .ne. 0) then
         do eg_to=1,egmax
            do eg_from=eg_to+1,egmax
               if( abs(sigma_scat(xs_mat(m)%mat,1,eg_to,eg_from)%xs) > 2.24E-16_d_t .and. most_thermal>eg_to ) most_thermal=eg_to
            end do
         end do
      end if 

    ! Compute tsigs
      do eg_from=1,egmax
         tsigs(xs_mat(m)%mat,eg_from)%xs=0.0_d_t
         do eg_to=1,egmax
            tsigs(xs_mat(m)%mat,eg_from)%xs=tsigs(xs_mat(m)%mat,eg_from)%xs+sigma_scat(xs_mat(m)%mat,1,eg_to,eg_from)%xs
         end do
      end do
       
    end do

  ! Close mesh file

    close(11)

  ! set neven

    neven=1_li+scatt_ord+(scatt_ord+1_li)*scatt_ord/2_li

  ! set scat_mult

    allocate(scat_mult(0:scatt_ord,0:scatt_ord))
    scat_mult=0.0_d_t
    if(scat_mult_flag.eq.1_li) then
       do l=0,scatt_ord
          do m=0,l
             if(m.eq.0_li) then
                scat_mult(l,m)=1.0_d_t/real(2_li*l+1_li,d_t)
             else
                scat_mult(l,m)=2.0_d_t/real(2_li*l+1_li,d_t)
             end if
          end do
       end do
    else
       do l=0,scatt_ord
          do m=0,l   
             if(m.eq.0_li) then
                scat_mult(l,m)=1.0_d_t
             else
                scat_mult(l,m)=2.0_d_t
             end if 
          end do
       end do
    end if

  end subroutine read_xs


end module read_cross_section_module

