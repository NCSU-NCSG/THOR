module sweep_module
!***********************************************************************
!
! Sweep module selects discrete ordinate directions and performs sweep
! over the eight Cartesian octants. Also contains subroutines that
! determine surface normals, incoming/outgoing declarations, and 
! initiate infinite queue over all cells. Finally, this module prepares
! all input data for the transport kernel to calculate angular fluxes.  
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

  use cell_splitting_module

  implicit none

! Define angular flux data type

  type angular_flux
     real(kind=d_t), dimension(:,:)  , allocatable :: vol_flux
     real(kind=d_t), dimension(:,:,:), allocatable :: face_flux
  end type angular_flux

contains

  !> This subroutine performs a transport sweep, i.e. inversion of (Omega * nabla + sigt_g)
  !> It explicitly loops over all angular directions. For each angular direction a directionally
  !> dependent source is computed (scattering + fission + fixed source). Then you do the mesh sweep
  !> on that source. The mesh sweep is done in subroutine queue. It is called queue because on an
  !> unstructured mesh the mesh sweep path is computed using a Breadth First Search algorithm that
  !> queues up all the elements in permissible order. Note: A lot of the stuff done here is for 
  !> reflective boundary conditions. These are implicit boundary conditions that need to be stored
  !> to ensure convergence. Angular fluxes, in general, are not stored in THOR!  
  subroutine sweep(eg,sc_flux,src,rs,reflected_flux,LL,U,Lf,Uf)
  !*********************************************************************
  !
  ! Subroutine sweep loops over all angular directions and call 
  ! transport solver
  !
  !*********************************************************************

  ! Pass group number
    
    integer(kind=li), intent(in) :: eg,rs

  ! Pass scalar flux
 
    real(kind=d_t) :: sc_flux(num_moments_v,namom,num_cells)

  ! Pass source
   
    real(kind=d_t) :: src(num_moments_v,namom,num_cells)

  ! Pass reflected flux

    real(kind=d_t), dimension(num_moments_f,rs,8,nangle) :: reflected_flux

  ! Pass pre-computed matrices

    real(kind=d_t), dimension(num_moments_v,num_moments_v) :: LL, U
    real(kind=d_t), dimension(num_moments_f,num_moments_f) :: Lf, Uf

  ! Declare local angular flux

    type(angular_flux) :: an_flux

  ! Declare face_known array=> now pre-filled in sweep and then passed to queue 

    integer(kind=1), dimension(0:3,num_cells) :: face_known

  ! Define temporary variables

    integer(kind=li) :: alloc_stat, q, oct, octt, octant, i, f, l, rcell,j,k,m,indx,nk,nf
    type(vector) :: omega

  ! Define source along direction q

    real(kind=d_t) :: dir_src(num_moments_v,num_cells)
 
  ! Temporary variables

    integer(kind=li) :: tpe, mate,giflx, parallel_i
    
    integer ::rank,mpi_err, localunit, num_p, optimal_tasks
    
    
  ! Define sc_flux parallel recieve buffer & reflected buffer

    real(kind=d_t) :: sc_flux_buffer (num_moments_v,namom,num_cells)    
    real(kind=d_t) :: reflected_buffer(num_moments_f, rs, 8, nangle)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_p, mpi_err)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)
    optimal_tasks = ceiling((nangle*8.0)/(num_p)) 
    sc_flux_buffer =zero
    reflected_buffer= zero
  ! Allocate angular flux

    allocate(an_flux%vol_flux(num_moments_v,num_cells)     ,stat=alloc_stat)
    if(alloc_stat /= 0) call stop_thor(2_li)
    allocate(an_flux%face_flux(num_moments_f,0:3,num_cells),stat=alloc_stat)
    if(alloc_stat /= 0) call stop_thor(2_li)
    
  ! Begin parallel loop over quadrature octant

    do parallel_i = 1, optimal_tasks
      k = parallel_map_l2g(parallel_i, rank+1)
      if (k .eq. 0) exit
      oct = mod(k,8)+1
      q = ceiling(k/8.0)
          octant= oct !ordering(oct)
          if(octant == 1)then
             omega%x1=quadrature(q)%mu%x1
             omega%x2=quadrature(q)%mu%x2
             omega%x3=quadrature(q)%mu%x3
          elseif(octant == 2)then
             omega%x1=-1*quadrature(q)%mu%x1
             omega%x2=quadrature(q)%mu%x2
             omega%x3=quadrature(q)%mu%x3
          elseif(octant == 3)then
             omega%x1=-1*quadrature(q)%mu%x1
             omega%x2=-1*quadrature(q)%mu%x2
             omega%x3=quadrature(q)%mu%x3
          elseif(octant == 4)then
             omega%x1=quadrature(q)%mu%x1
             omega%x2=-1*quadrature(q)%mu%x2
             omega%x3=quadrature(q)%mu%x3
          elseif(octant == 5)then
             omega%x1=quadrature(q)%mu%x1
             omega%x2=quadrature(q)%mu%x2
             omega%x3=-1*quadrature(q)%mu%x3
          elseif(octant == 6)then
             omega%x1=-1*quadrature(q)%mu%x1
             omega%x2=quadrature(q)%mu%x2
             omega%x3=-1*quadrature(q)%mu%x3
          elseif(octant == 7)then
             omega%x1=-1*quadrature(q)%mu%x1
             omega%x2=-1*quadrature(q)%mu%x2
             omega%x3=-1*quadrature(q)%mu%x3
          else
             omega%x1=quadrature(q)%mu%x1
             omega%x2=-1*quadrature(q)%mu%x2
             omega%x3=-1*quadrature(q)%mu%x3
          end if

          ! Zero angular flux and directed source

          dir_src=0.0_d_t             
          an_flux%vol_flux  = 0.0_d_t
          an_flux%face_flux = 0.0_d_t
          face_known = 0

          ! Set boundary faces and mark as known

            ! 1. Vacuum boundary conditions
            do i=1,vside_cells
               k=vb_cells(i)%cell
               f=vb_cells(i)%face
               if ( (omega .dot. outward_normal(k,f)) < 0.0_d_t )then 
                 face_known(f,k) = 1
                 ! no nneed to assign 0 to face flux, already done before !!
               end if
            end do   

            ! 2. Reflective boundary conditions
            
            do i=1,rside_cells
               k=rb_cells(i)%cell
               f=rb_cells(i)%face
               if ( (omega .dot. outward_normal(k,f)) < 0.0_d_t )then 
                 face_known(f,k) = 1
                 ! reflected flux(oct) stores outflow for octant oct
                 ! so we need to find the mate of oct to get the right inflow
                 tpe=refl_face_tpe(i)  
                 if      (tpe .eq. 1_li .or. tpe .eq. -1_li) then
                    mate = mu_mate(octant)
                 else if (tpe .eq. 2_li .or. tpe .eq. -2_li) then
                    mate = eta_mate(octant)
                 else if (tpe .eq. 3_li .or. tpe .eq. -3_li) then
                    mate = xi_mate(octant)
                 end if
                 an_flux%face_flux(:,f,k) = reflected_flux(:,i,mate,q)
               end if
            end do
            
            
            ! 3. Fixed inflow boundary conditions
           
            if(page_iflw.eq.1_li) then
              giflx=1_li
            else
              giflx=eg
            end if
            
            do i=1,fside_cells
               k=fb_cells(i)%cell
               f=fb_cells(i)%face
               if ( (omega .dot. outward_normal(k,f)) < 0.0_d_t )then   
                 face_known(f,k) = 1
                 an_flux%face_flux(:,f,k) = binflx(:,i,octant,q,giflx)
               end if
            end do 

            ! 4. Faces that are assumed known for breaking cycles
 
            is_cycle=0 
            do i=1, neldep(octant,q)
               k=eldep(octant,q,eg)%cells(i)
               f=eldep(octant,q,eg)%faces(i)
               is_cycle(f,k)   = 1 
               face_known(f,k) = 1
               an_flux%face_flux(:,f,k) = eldep(octant,q,eg)%face_fluxes(:,i) 
            end do 

          ! Compute directed source - so far only isotropic 

          do i=1,num_cells
            ! even contributions
            do l=0,scatt_ord
              do m=0,l
                 indx=1_li+m+(l+1_li)*l/2_li
                 do k=1,num_moments_v
                   dir_src(k,i)=dir_src(k,i)+Ysh(q,octant,indx)*src(k,indx,i) 
                 end do
              end do
            end do
            ! odd contributions
            do l=1,scatt_ord
              do m=1,l
                 indx=neven+m+(l-1_li)*l/2_li
                 do k=1,num_moments_v
                   dir_src(k,i)=dir_src(k,i)+Ysh(q,octant,indx)*src(k,indx,i) 
                 end do
              end do
            end do 
          end do

          ! If page_sweep == 1 read sweep_path from file

          if(page_sweep.eq.1) then 
             sweep_path=0
             read(99,rec=8*(q-1)+octant) sweep_path(:,1,1)
          end if

          ! Call queue for mesh sweep 
          
          if (sweep_tpe .eq.1) then
            call queue(eg,q,octant,dir_src,an_flux,LL,U,Lf,Uf,omega,face_known)
          else  
            call stop_thor(9_li)
          end if
          
          ! Increment scalar flux using 

            call angular_moments(octant,q,sc_flux,an_flux)

          ! update reflected flux array, loop through r           
            do i=1,rside_cells
               k=rb_cells(i)%cell
               f=rb_cells(i)%face
               if ( (omega .dot. outward_normal(k,f)) > 0.0_d_t )then
                 ! store outflow on reflective faces 
                 reflected_buffer(:,i,octant,q) = an_flux%face_flux(:,f,k)
               end if
            end do
          ! update the faces that are assumed known for cycles
            do i=1, neldep(octant,q)
               k=eldep(octant,q,eg)%cells(i)
               f=eldep(octant,q,eg)%faces(i)
               nk=adjacency_list(k,f)%cell 
               nf=adjacency_list(k,f)%face
               eldep(octant,q,eg)%face_fluxes(:,i)=an_flux%face_flux(:,nf,nk)
            end do
    end do
    reflected_flux=reflected_buffer
  ! Deallocation
    call MPI_AllREDUCE(reflected_flux, reflected_buffer, num_moments_f*rs*8*nangle, &
                                                   MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpi_err)
    reflected_flux=reflected_buffer
    call MPI_AllREDUCE(sc_flux, sc_flux_buffer, num_moments_v*namom*num_cells, &
                                                   MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpi_err)
    sc_flux = sc_flux_buffer
    deallocate(an_flux%vol_flux,an_flux%face_flux)

  end subroutine sweep

!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

  !> This subroutine performs a mesh sweep. It uses the sweep_path variable to obtain
  !> the sweep order.
  subroutine queue(eg,q,octant,qm,an_flux,LL,U,Lf,Uf,omega,face_known)
  !*********************************************************************
  !
  ! Subroutine queue uses the pre-computed sweep paths to perform a 
  ! mesh sweep.
  !
  !*********************************************************************

  ! Pass group and angular indices

    integer(kind=li), intent(in) :: eg, q, octant

  ! Pass qm
 
    real(kind=d_t) :: qm(num_moments_v,num_cells)

  ! Pass angular flux

    type(angular_flux), intent(inout) :: an_flux

  ! Pass pre-computed matrices

    real(kind=d_t), dimension(num_moments_v,num_moments_v) :: LL, U
    real(kind=d_t), dimension(num_moments_f,num_moments_f) :: Lf, Uf

  ! Pass logical array face known

    integer(kind=1), dimension(0:3,num_cells) :: face_known

  ! Pass omega

    type(vector) :: omega 

  ! Local variables

    integer(kind=li) :: alloc_stat, i,j, oct, f, l, case
    real(kind=d_t) :: sigmat
    logical :: all_incoming_known
    type(vector) :: n0, n1, n2, n3
    integer(kind=li) :: soct,sq
    
    integer ::rank,mpi_err, num_p
    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_p, mpi_err)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)

  ! Set soct and sq depending on whether page_sweep=0,1

    if      (page_sweep .eq. 0) then
       soct=octant;sq=q
    else if (page_sweep .eq. 1) then
       soct=1_li;sq=1_li
    end if

  ! Initiate loop over all cells
    
    do j=1,num_cells

       i=sweep_path(j,soct,sq)
       
       sigmat=dens_fact(cells(i)%reg)*sigma_t(reg2mat(cells(i)%reg),eg)%xs

       n0=outward_normal(i,0)
       n1=outward_normal(i,1)
       n2=outward_normal(i,2)
       n3=outward_normal(i,3)
       
       call cell_orientation(omega,n0,n1,n2,n3,case)
       
       call cell_splitting(sigmat,qm(:,i),an_flux%vol_flux,an_flux%face_flux,   &
                           octant,LL,U,Lf,Uf,omega,i,face_known,case,n0,n1,n2,n3)
              
    end do
    

  end subroutine queue

!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

  
  !> This subroutine orientation determines whether cell's faces are
  !> incoming or outgoing and assigns them a 'case' for the transport
  !> calculation
  subroutine cell_orientation(omega,n0,n1,n2,n3,case)
  !*********************************************************************
  !
  ! Subroutine cell orientation determines whether cell's faces are
  ! incoming or outgoing and assigns them a 'case' for the transport
  ! calculation
  !
  !*********************************************************************
  ! Define temporary variables

    integer(kind=li) :: incoming, outgoing
    integer(kind=li), intent(inout) :: case
    type(vector), intent(inout) :: omega
    type(vector), intent(in) :: n0, n1, n2, n3

  ! Determine which case cell corresponds with respect to current ordinate

    incoming=0
    outgoing=0

    if((omega .dot. n0) < 0.0)then
       incoming=incoming+1
    elseif((omega .dot. n0) == 0.0)then
    else
       outgoing=outgoing+1
    end if

    if((omega .dot. n1) < 0.0)then
       incoming=incoming+1
    elseif((omega .dot. n1) == 0.0)then
    else
       outgoing=outgoing+1
    end if

    if((omega .dot. n2) < 0.0)then
       incoming=incoming+1
    elseif((omega .dot. n2) == 0.0)then
    else
       outgoing=outgoing+1
    end if

    if((omega .dot. n3) < 0.0)then
       incoming=incoming+1
    elseif((omega .dot. n3) == 0.0)then
    else
       outgoing=outgoing+1
    end if

    if(incoming == 3 .and. outgoing == 1)then
       case=1
    elseif(incoming == 1 .and. outgoing == 3)then
       case=2
    elseif(incoming == 2 .and. outgoing == 2)then
       case=3
    elseif(incoming == 2 .and. outgoing == 1)then
       case=4
    elseif(incoming == 1 .and. outgoing == 2)then
       case=5
    elseif(incoming == 1 .and. outgoing == 1)then
       case=6
    else
       call stop_thor(10_li)
    end if

  end subroutine cell_orientation

!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!

  subroutine check_incoming(face,case,all_incoming_known)
  !*********************************************************************
  !
  ! Subroutine check incoming determines if cell is ready for transport
  ! calculation (all incoming angular fluxes are known)
  !
  !*********************************************************************
  ! Define geometry variables

    logical, dimension(0:3), intent(in) :: face
    logical, intent(inout) :: all_incoming_known
    integer(kind=li), intent(inout) :: case

  ! Define temporary variables

    integer(kind=li) :: known_faces, f, l

    known_faces=0

  ! Check if incoming face moments to all computed orders for cell are known

  ! Case 1 (3 incoming, 1 outgoing)

    if(case == 1)then
       do f=0, 3
          l=1
          if(face(f) .eqv. .true.)then
             known_faces=known_faces+1
          end if
       end do

       if(known_faces == 3)then
          all_incoming_known=.true.
       else
          all_incoming_known=.false.
       end if


  ! Case 2 (1 incoming, 3 outgoing)

    elseif(case == 2)then
       do f=0, 3
          l=1
          if(face(f) .eqv. .true.)then
             known_faces=known_faces+1
          end if
       end do

       if(known_faces == 1)then
          all_incoming_known=.true.
       else
          all_incoming_known=.false.
       end if

  ! Case 3 (2 incoming, 2 outgoing)

    elseif(case == 3)then
       do f=0, 3
          l=1
          if(face(f) .eqv. .true.)then
             known_faces=known_faces+1
          end if
       end do

       if(known_faces == 2)then
          all_incoming_known=.true.
       else
          all_incoming_known=.false.
       end if

  ! Case 4 (2 incoming, 1 outgoing)

    elseif(case == 4)then
       do f=0, 3
          l=1
          if(face(f) .eqv. .true.)then
             known_faces=known_faces+1
          end if
       end do

       if(known_faces == 2)then
          all_incoming_known=.true.
       else
          all_incoming_known=.false.
       end if

  ! Case 5 (1 incoming, 2 outgoing)

    elseif(case == 5)then
       do f=0, 3
          l=1
          if(face(f) .eqv. .true.)then
             known_faces=known_faces+1
          end if
       end do

       if(known_faces == 1)then
          all_incoming_known=.true.
       else
          all_incoming_known=.false.
       end if

   ! Case 6 (1 incoming, 1 outgoing)

    elseif(case == 6)then
       do f=0, 3
          l=1
          if(face(f) .eqv. .true.)then
             known_faces=known_faces+1
          end if
       end do

       if(known_faces == 1)then
          all_incoming_known=.true.
       else
          all_incoming_known=.false.
       end if

    else
       call stop_thor(11_li)
    end if
    
  end subroutine check_incoming

!-----------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------!
 
  !> This subroutine evaluates angular flux moments phi_lm from the angular fluxes. 
  !> Note: The angular flux moments are accumulated on-the-fly from angular fluxes so that
  !> angular fluxes do not need to be stored.  
  subroutine angular_moments(octant,q,sc_flux,an_flux)
  !*********************************************************************
  !
  ! Subroutine angular moments computes the real spherical harmonics
  ! angular moments of the angular flux
  !
  !*********************************************************************
  ! Pass input parameters

    integer(kind=li), intent(in) :: q,octant

  ! Declare angular and scalar fluxes

    real(kind=d_t) :: sc_flux(num_moments_v,namom,num_cells)
    type(angular_flux), intent(in) :: an_flux 

  ! Define temporary variables

    integer(kind=li) :: i, l, order, k, indx, m

  ! Compute scalar flux angular moments, right now only isotropic
  ! flux is accumulated

    do i=1, num_cells
      ! even contributions  
      do l=0,scatt_ord
        do m=0,l
          do k=1, num_moments_v
             indx=1_li+m+(l+1_li)*l/2_li
             sc_flux(k,indx,i)=sc_flux(k,indx,i)+&
                  (1.0_d_t/8.0_d_t)*quadrature(q)%wt*Ysh(q,octant,indx)*an_flux%vol_flux(k,i)
          end do
        end do
      end do
      ! odd contributions  
      do l=1,scatt_ord
        do m=1,l
          do k=1, num_moments_v
             indx=neven+m+(l-1_li)*l/2_li
             sc_flux(k,indx,i)=sc_flux(k,indx,i)+&
                  (1.0_d_t/8.0_d_t)*quadrature(q)%wt*Ysh(q,octant,indx)*an_flux%vol_flux(k,i)
          end do
        end do
      end do 
    end do

  end subroutine angular_moments

end module sweep_module

