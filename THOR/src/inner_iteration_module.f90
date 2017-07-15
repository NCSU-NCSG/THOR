module inner_iteration_module
!***********************************************************************
!
! Inner iteration module calls sweep subroutine and updates distributed
! source
!
!***********************************************************************

! User derived-type modules
  
  use mpi
  use types
  use parameter_types
  use filename_types
  use vector_types
  use cross_section_types
  use geometry_types
  use angle_types
  use multindex_types
  use global_variables

! Use modules that pertain setting up problem

  use sweep_module
  use termination_module

  implicit none

contains

  !> This subroutine performs and inner iteration (compare Alg. 1 in primer)
  subroutine inner_iteration(eg,sc_flux,q_external,LL,U,Lf,Uf,rs,reflected_flux,prnt)
  !*********************************************************************
  !
  ! Subroutine inner iteration calls sweep subroutine to compute all 
  ! fluxes average cell scalar
  !
  !*********************************************************************

  ! Inner Iterations executed for group eg

    integer(kind=li), intent(in) :: eg

  ! Declare angular and scalar flux types used globally

    real(kind=d_t) :: sc_flux(num_moments_v,namom,num_cells)

  ! Declare source types

    real(kind=d_t) :: q_external(num_moments_v,namom,num_cells)

  ! Declare pre-computed matrices

    real(kind=d_t), dimension(num_moments_v,num_moments_v) :: LL, U
    real(kind=d_t), dimension(num_moments_f,num_moments_f) :: Lf, Uf
   
  ! Define reflected flux array

    integer(kind=li) :: rs
    real(kind=d_t),dimension(num_moments_f,rs,8,nangle) :: reflected_flux

  ! print convergence monitor flag

    logical :: prnt 

  ! Define temporary variables

    integer(kind=li) :: l, i, q, octant, order, alloc_stat,&
                        n,m, k, indx  
    real(kind=d_t) :: error,te,ts

  ! Define self-scatter source

    real(kind=d_t) :: self_scatter(num_moments_v,namom,num_cells)

  ! Define old scalar flux
 
    real(kind=d_t) :: sc_flux_old (num_moments_v,namom,num_cells)
    
  ! Define MPI environment and get rank
    integer ::rank,mpi_err, localunit
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)



  ! Set self_scatter and old scalar flux to 0
   
    self_scatter = zero
    sc_flux_old  = zero

  ! write header for convergence monitor
    if (rank .eq. 0) then
      if(prnt) write(6,102) '  grp  itn       error        time'  
    end if
  ! Begin inner iteration

    do inner=1, max_inner

  ! start timer

      call cpu_time(ts) 

  ! Recompute self-scattering source and add external source

      do i=1,num_cells
         ! even contribution
         do l=0,scatt_ord
           do m=0,l
             indx=1_li+m+(l+1_li)*l/2_li
             do k=1, num_moments_v
              self_scatter(k,indx,i) = scat_mult(l,m)                             *&
                               sigma_scat(reg2mat(cells(i)%reg),l+1,eg,eg)%xs     *&
                               dens_fact(cells(i)%reg)                            *&
                               sc_flux(k,indx,i) + q_external(k,indx,i)
             end do
           end do
         end do
         ! odd contribution
         do l=1,scatt_ord
           do m=1,l
             indx=neven+m+(l-1_li)*l/2_li
             do k=1, num_moments_v
              self_scatter(k,indx,i) = scat_mult(l,m) *&
                               sigma_scat(reg2mat(cells(i)%reg),l+1,eg,eg)%xs     *&
                               dens_fact(cells(i)%reg)                            *&
                               sc_flux(k,indx,i) + q_external(k,indx,i)
             end do
           end do
         end do

      end do

  ! Call sweep algorithm
  
      sc_flux=zero
      call sweep(eg,sc_flux,self_scatter,rs,reflected_flux,LL,U,Lf,Uf)
  ! Computer error

       max_error(eg)=zero

       do i=1, num_cells
          if(abs(sc_flux(1,1,i)) > 1e-12)then
             error=abs( (  sc_flux(1,1,i) - sc_flux_old(1,1,i) )/&
                           sc_flux(1,1,i) )
          else
             error=abs(    sc_flux(1,1,i) - sc_flux_old(1,1,i) )
          end if
          if(error >= max_error(eg))then
             max_error(eg)=error
          end if
       end do

  ! Test convergence of current-group scalar flux

       if(inner > 1 .and. max_error(eg) < inner_conv)then
          go to 10
       end if

       if(inner < max_inner)then
         do i=1, num_cells
           do n=1,namom
             do l=1, num_moments_v
               sc_flux_old(l,n,i)=sc_flux(l,n,i)
             end do
           end do
          end do
       end if

    ! stop timer

      call cpu_time(te) 

    ! write convergence monitor
      if (rank .eq. 0) then
        if(prnt) then
          write(6,101) eg,inner,max_error(eg),te-ts,' % '
          flush(6)
        end if
        if(prnt .and. print_conv.eq.1) then 
          write(21,101) eg,inner,max_error(eg),te-ts,' % '
          flush(21)
        end if
        101 FORMAT (1X,2I5,2ES12.4,A)
        102 FORMAT (1X,A)
      end if
    end do

    inner=inner-1_li

10  continue

    tot_nInners=tot_nInners+inner

  end subroutine inner_iteration

end module inner_iteration_module

