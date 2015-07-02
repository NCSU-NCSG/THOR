module setup_module
!***********************************************************************
!
! Setup module contains all subroutines to setup problem
!
!***********************************************************************

! Use derived-type modules

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

  use read_module
  use quadrature_module
  use sweep_module
  use termination_module

  implicit none

contains

  !> This subroutine reads the input and prepares data for execution
  subroutine setup
  !*********************************************************************
  !
  ! Subroutine setup calls routines to read and prepare input for
  ! execution
  !
  !*********************************************************************

  ! Local variables

    integer(kind=li) :: i,alloc_stat
    integer(kind=li) :: q,oct,summe
    type(vector) :: v0, v1, v2, v3, n0, n1, n2, n3
    integer(kind=li) :: reflective(3,2)
    real(kind=d_t) :: ts,te
    real(kind=d_t), dimension(3,3) :: J, J_inv
    real(kind=d_t)                 :: det_J,tot_vol
    logical :: existence
    integer ::rank,mpi_err, localunit
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)

  ! Read input

    call read

  ! Assign local angles for parallel analysis

    call assign_work

  ! Allocate and initialize density factors and region volume

    allocate(dens_fact(minreg:maxreg));dens_fact=1.0_d_t
    allocate(reg_vol  (minreg:maxreg));reg_vol  =0.0_d_t

    do i=1, num_cells
       v0=vertices(cells(i)%R(0))%v
       v1=vertices(cells(i)%R(1))%v
       v2=vertices(cells(i)%R(2))%v
       v3=vertices(cells(i)%R(3))%v
       call cell_jacobian(v0,v1,v2,v3,J,det_J,J_inv)
       cells(i)%volume=(1.0_d_t/6.0_d_t)*det_J
       reg_vol(cells(i)%reg)=reg_vol(cells(i)%reg)+cells(i)%volume
    end do

    if(dfact_opt.ne.0) then
      call read_density_factors
      if(dfact_opt.eq.1) then
        do i=minreg,maxreg
          if(dens_fact(i)<0.0_d_t .or. reg_vol(i)<2.24d-12) then
            dens_fact(i)=1.0_d_t
          else
            dens_fact(i)=dens_fact(i)/reg_vol(i)
          end if
        end do
      else
        do i=minreg,maxreg
          if(dens_fact(i)<0.0_d_t) dens_fact(i)=1.0_d_t
        end do
      end if
    end if

    call write_reg_info

  ! Set comvergence flag to 0 => no convergence

    conv_flag=0_li

  ! Set iteration counters to 0

    tot_nInners = 0_li
    outer       = 0_li
    tot_kit     = 0_li
    nit         = 0_li

  ! Set number of angular flux moments to be computed

    namom=(scatt_ord+1)**2

  ! Set number of retained flux iterates

    niter=3_li

  ! Estimate memory consumption
    if (rank .eq. 0) then
      call memory_estimate(mem_req)
      write(6,*) "--------------------------------------------------------"
      write(6,*) "   Memory estimate   "
      write(6,*) "--------------------------------------------------------"
      write(6,102) 'Memory estimate: ',mem_req,' MB'
      write(6,*)
      102 FORMAT(1X,A18,ES12.4,A3)
    end if
  ! Pre-compute the cell outward normal vectors

    allocate(outward_normal(num_cells,0:3),stat=alloc_stat)
    if(alloc_stat /=0) call stop_thor(2_li)

  ! Generate outward normal vectors for all cells for performance

    do i=1, num_cells
       v0=vertices(cells(i)%R(0))%v
       v1=vertices(cells(i)%R(1))%v
       v2=vertices(cells(i)%R(2))%v
       v3=vertices(cells(i)%R(3))%v
       call generate_normals(v0,v1,v2,v3,n0,n1,n2,n3)
       outward_normal(i,0)=n0
       outward_normal(i,1)=n1
       outward_normal(i,2)=n2
       outward_normal(i,3)=n3
    end do

    if (rside_cells > 0) then

    ! Check reflective boundary faces: They need to be on a face that is
    ! perpendicuar to either x,y or z-axis.

      call check_reflective

    ! Classify reflective boundary faces by axis they are orthogonal to:
    ! +1 = outward normal along positive x axis, -1 = outward normal along negative x-axis
    ! +2 = outward normal along positive y axis, -2 = outward normal along negative y-axis
    ! +3 = outward normal along positive z axis, -3 = outward normal along negative z-axis

      call classify_reflective

    end if

    ! Set up the variables to memorize rings in the mesh sweep

    allocate(neldep(8,nangle),stat=alloc_stat)
    if(alloc_stat /=0) call stop_thor(1_li)
    neldep=0_li
    allocate(eldep(8,nangle,egmax) ,stat=alloc_stat)
    if(alloc_stat /=0) call stop_thor(1_li)
    allocate( is_cycle(0:3,num_cells) )

    ! Pre-compute sweep if sweep_tpe = 1

    if (sweep_tpe .eq. 1) then

       ! allocate sweep path
         if      (page_sweep .eq. 0) then
           allocate( sweep_path(num_cells,8,nangle),stat=alloc_stat)
           if(alloc_stat /=0) call stop_thor(1_li)
         else if (page_sweep .eq. 1) then
           allocate( sweep_path(num_cells,1,1),stat=alloc_stat)
           if(alloc_stat /=0) call stop_thor(1_li)
           ! open scratch file that contains sweep path
           inquire(file="sweep_path.pg",exist=existence)
           if( existence .eqv. .true.) then
              open(unit=99,file='sweep_path.pg',status='replace',form='unformatted',access='direct',recl=4*num_cells)
           else
              open(unit=99,file='sweep_path.pg',status='new'    ,form='unformatted',access='direct',recl=4*num_cells)
           end if
         end if
         sweep_path=0_li

       ! pre-compute sweep
         if (rank .eq. 0) then
           write(6,*) "--------------------------------------------------------"
           write(6,*) "   Precomputing sweep path   "
           write(6,*) "--------------------------------------------------------"
         end if
!FIX ME - pre sweep efficiency and threading
         call cpu_time(ts)
         call pre_sweep
         call cpu_time(te)
         if (rank .eq. 0) then
           write(6,101) '-- Finished precomputing sweep path. Time (sec.) ', te-ts
           write(6,*)
         end if
    end if

    ! Determine which boundary faces are reflective and save in reflective

    reflective=0_li
    do i=1,rside_cells
       q=abs(refl_face_tpe(i))
       if(refl_face_tpe(i)>0) then
          reflective(q,2)=1
       else
          reflective(q,1)=1
       end if
    end do

    ! Determine ordering of octants

    call set_ordering(reflective)

    ! Make sure that if JFNK is selected reflective BC in opposite faces is not
    ! permitted

    if(problem .eq. 1 .and. eig_switch .eq.1 .and. rd_method > 1) then
      do i=1,3
         if(reflective(i,1) .eq. 1 .and. reflective(i,2).eq.1 ) call stop_thor(31_li)
      end do
    end if

   101 FORMAT(1X,A,ES12.4)

  ! set variables that should eventually be moved to input

    exmax = 0.5_d_t
    extol = 0.05_d_t

  end subroutine setup

!------------------------------------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------------------------------------!

  subroutine write_reg_info
  !*********************************************************************
  !
  ! Write region information
  !
  !*********************************************************************

  ! Local variables

    integer(kind=li) :: i

  ! Write region information
    if (rank .eq. 0) then
      write(6,*)
      write(6,*) "--------------------------------------------------------"
      write(6,*) "   Region Information  "
      write(6,*) "--------------------------------------------------------"
      write(6,*)
      write(6,*) "     Region ID","   Material ID","  Region Volume"," Density Factor"
      do i = minreg,maxreg
         write(6,104) i,reg2mat(i),reg_vol(i),dens_fact(i)
      end do
    end if

    104 FORMAT(1X,I14,I14,3ES15.4)

  end subroutine

!------------------------------------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------------------------------------!

  subroutine read_density_factors
  !*********************************************************************
  !
  ! Reads density factor file
  !
  !*********************************************************************

    ! Local variables
    logical :: ex
    integer(kind=li) :: i

    ! Open and read file, then close, terminate THOR if file does not exist

    inquire(file=dens_fact_filename,exist=ex)
    if(ex) then
      open(file=dens_fact_filename,unit=51)
      read(51,*) (dens_fact(i),i=minreg,maxreg)
      close(unit=51)
    else
      call stop_thor(34_li)
      return
    end if
  end subroutine
!------------------------------------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------------------------------------!

  subroutine set_ordering(reflective)
  !*********************************************************************
  !
  ! Subroutine computes ordering such that sweep starts from vacuum BC
  !
  !*********************************************************************

    ! --- Arguments
    integer(kind=li) :: reflective(3,2)
    integer(kind=li) :: xd = 0,yd = 0,zd=0

    ! --- Set ordering
    !  -- x-dir
    if      (reflective(1,1) .eq. 0_li .and. reflective(1,2) .eq. 0_li) then
       xd=1_li
    else if (reflective(1,1) .eq. 1_li .and. reflective(1,2) .eq. 0_li) then
       xd=-1_li
    else if (reflective(1,1) .eq. 0_li .and. reflective(1,2) .eq. 1_li) then
       xd=1_li
    else if (reflective(1,1) .eq. 1_li .and. reflective(1,2) .eq. 1_li) then
       xd=1_li
    end if
    !  -- y-dir
    if      (reflective(2,1) .eq. 0_li .and. reflective(2,2) .eq. 0_li) then
       yd=1_li
    else if (reflective(2,1) .eq. 1_li .and. reflective(2,2) .eq. 0_li) then
       yd=-1_li
    else if (reflective(2,1) .eq. 0_li .and. reflective(2,2) .eq. 1_li) then
       yd=1_li
    else if (reflective(2,1) .eq. 1_li .and. reflective(2,2) .eq. 1_li) then
       yd=1_li
    end if
    !  -- z-dir
    if      (reflective(3,1) .eq. 0_li .and. reflective(3,2) .eq. 0_li) then
       zd=1_li
    else if (reflective(3,1) .eq. 1_li .and. reflective(3,2) .eq. 0_li) then
       zd=-1_li
    else if (reflective(3,1) .eq. 0_li .and. reflective(3,2) .eq. 1_li) then
       zd=1_li
    else if (reflective(3,1) .eq. 1_li .and. reflective(3,2) .eq. 1_li) then
       zd=1_li
    end if

    !  -- Chose octant ordering
    if      (xd .eq.  1_li .and. yd .eq.  1_li .and. zd .eq.  1_li) then
       ordering =(/1,2,3,4,5,6,7,8/)
    else if (xd .eq.  1_li .and. yd .eq.  1_li .and. zd .eq. -1_li) then
       ordering =(/5,6,7,8,1,2,3,4/)
    else if (xd .eq.  1_li .and. yd .eq. -1_li .and. zd .eq.  1_li) then
       ordering =(/4,3,8,7,2,1,5,4/)
    else if (xd .eq. -1_li .and. yd .eq.  1_li .and. zd .eq.  1_li) then
       ordering =(/2,3,6,7,1,4,5,8/)
    else if (xd .eq.  1_li .and. yd .eq. -1_li .and. zd .eq. -1_li) then
       ordering =(/8,7,6,5,4,3,2,1/)
    else if (xd .eq. -1_li .and. yd .eq.  1_li .and. zd .eq. -1_li) then
       ordering =(/6,7,8,5,2,3,4,1/)
    else if (xd .eq. -1_li .and. yd .eq. -1_li .and. zd .eq.  1_li) then
       ordering =(/3,4,1,2,7,8,5,6/)
    else if (xd .eq. -1_li .and. yd .eq. -1_li .and. zd .eq. -1_li) then
       ordering =(/7,6,8,5,3,2,4,1/)
    end if

  end subroutine

!------------------------------------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------------------------------------!

  !> This subroutine pre-computes the sweep path and saves it in
  !> variable sweep_path
  subroutine pre_sweep
  !*********************************************************************
  !
  ! Subroutine pre-sweep pre-computes the sweep path and saves it in
  ! variable sweep_path
  !
  !*********************************************************************

  ! Define temporary variables

    integer(kind=li) :: alloc_stat, q, oct, octt, octant, i, f, l, rcell,j,k
    type(vector) :: omega

  ! Temporary variables

    integer(kind=li) :: tpe, mate
    real(kind=d_t) :: te,ts

    integer ::rank,mpi_err, localunit, num_p, optimal_tasks
    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_p, mpi_err)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)

  ! Begin parallel loop over quadrature octant
    optimal_tasks = ceiling((nangle*8.0)/(num_p))

    do i=1, optimal_tasks

      k = parallel_map_l2g(i, rank+1)
      if (k .eq. 0) exit
      oct = mod(k,8)+1
      q = ceiling(k/8.0)

!FIX ME - Octant ordering is lost here
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

      ! Start timer

      call cpu_time(ts)

      call pre_queue(q,octant,omega)

      ! End timer and print message

      call cpu_time(te)
!FIX ME - Only prints angles belonging to root process
      if (rank .eq. 0) then
        write(6,101) ' - Computation of sweep path for octant: ',octant,' angle: ',q,' . Ex. Time(sec.): ',te-ts
      end if


    end do

101 FORMAT(1X,A,I4,A,I4,A,ES12.4)
  end subroutine

!------------------------------------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------------------------------------!

  !> This subroutine performs a BFS of the mesh along direction omega.
  !> It deals with rings by severing the connection across one face
  !> within the ring and saving the cell/face #.
  subroutine pre_queue(q,octant,omega)
  !*********************************************************************
  !
  ! This subroutine performs a BFS of the mesh along direction omega.
  ! It deals with rings by severing the connection across one face
  ! within the ring and saving the cell/face #.
  !
  !*********************************************************************

    ! --- Arguments
    integer(kind=li) :: q,octant
    type(vector) :: omega

    ! --- Local variables
    integer(kind=li) :: queue(num_cells)
    integer(kind=1)  :: idegr(num_cells)
    integer(kind=1)  :: odegr(num_cells)
    integer(kind=li) :: qmarker,iqmarker,fmarker,ifmarker,finished
    integer(kind=li) :: g,i,j,neigh,pvt=0,cell,face
    integer(kind=1)  :: solved(num_cells)
    integer(kind=li) :: dadj(0:3,num_cells)
    type(vector)     :: onor
    type(upstream_face_le),pointer :: cycle_list,current_el,next

    ! --- Zero sweep_path if page_sweep==1
    if(page_sweep.eq.1) sweep_path(:,1,1)=0

    ! --- Prepare cycle list
    nullify(cycle_list);allocate(cycle_list)
    nullify(cycle_list%next)
    current_el => cycle_list

    ! --- Constructing the directed adjacency list
!-------------------------------
! Problem w/ omega.dot.onor ==0 not resolved
!-------------------------------
    dadj=0_li
    do cell=1,num_cells
       do face=0,3
          ! -- get current normal
          onor=outward_normal(cell,face)
          ! -- decide if inflow or outflow
          if      ( (omega .dot. onor) > 0.0_d_t) then   ! - outflow
             dadj(face,cell)=adjacency_list(cell,face)%cell
          else if ( (omega .dot. onor) < 0.0_d_t) then   ! - inflow
             dadj(face,cell)=-adjacency_list(cell,face)%cell
          end if
       end do
    end do

    ! --- Start algorithm by setting indegree and outdegree:
    !     indegree : # upstream dependencies
    !     outdegree: # downstream dependencies

    idegr=0 ; odegr=0

    do cell=1,num_cells
      do face=0,3
        if      (dadj(face,cell)>0) then
           odegr(cell)=odegr(cell)+1
        else if (dadj(face,cell)<0) then
           idegr(cell)=idegr(cell)+1
        end if
      end do
    end do

    ! --- Start main loop
    finished=0
    fmarker =0
    ifmarker=num_cells+1
    qmarker =0
    iqmarker=num_cells+1
    solved  =0

    do while(finished<num_cells)

      ! - Fill queue with all tets that initially have no upstream
      !   dependencies
      do i=1,num_cells
         if(idegr(i) .eq. 0 .and. solved(i) .eq. 0) then
            qmarker=qmarker+1
            queue(qmarker)=i
         end if
      end do

      ! -- Downstream trim
      do while(fmarker .lt. qmarker)

        ! - Get current tet
        finished=finished+1
        fmarker = fmarker+1
        cell=queue(fmarker)

        ! - list cell in order
        if      ( page_sweep .eq. 0) then
          sweep_path(fmarker,octant,q)=cell
        else if ( page_sweep .eq. 1) then
          sweep_path(fmarker,1,1)=cell
        end if

        solved(cell)  =1
        ! - Decrements the upstream dependencies of the downstream neighbors
        do j=0,3
           neigh=dadj(j,cell)
           if(neigh>0) then
              idegr(neigh)=idegr(neigh)-1
              if(idegr(neigh).eq.0 .and. solved(neigh) .eq. 0) then
                qmarker=qmarker+1
                queue(qmarker)=neigh
              end if
           end if
        end do

      end do

      ! -- If downstream trim did not touch all cells, then there are cycles
      if(finished<num_cells) then

      ! - Fill queue reversely with all tets that initially have no downstream
      !   dependencies
        do i=1,num_cells
          if(odegr(i) .eq. 0 .and. solved(i) .eq. 0) then
            iqmarker=iqmarker-1
            queue(iqmarker)=i
          end if
        end do

      ! -- Upstream trim
        do while(ifmarker .gt. iqmarker)

          ! - Get current tet
          finished=finished+1
          ifmarker = ifmarker-1
          cell=queue(ifmarker)

          ! - list cell in order
          if   ( page_sweep .eq. 0) then
            sweep_path(ifmarker,octant,q)=cell
          else if ( page_sweep .eq. 1) then
            sweep_path(ifmarker,1,1)=cell
          end if

          solved(cell)  =1
          ! - Decrements the upstream dependencies of the downstream neighbors
          do j=0,3
            neigh=dadj(j,cell)
            if(neigh<0) then
              odegr(abs(neigh))=odegr(abs(neigh))-1
              if(odegr(abs(neigh)).eq.0 .and. solved(abs(neigh)) .eq. 0) then
                iqmarker=iqmarker-1
                queue(iqmarker)=abs(neigh)
              end if
            end if
          end do

        end do

      end if

      if(finished<num_cells) then

      ! -- Find element with only one upstream dependence or use
      !    and element with two upstream dependencies
        j=3;i=0
        do while(j>1 .and. i.lt.num_cells)
           i=i+1
           if(solved(i) .eq. 0) then
             if(idegr(i)<j) then
                j=idegr(i)
                pvt=i
             end if
           end if
        end do

      ! -- List pvt in order and increment finished counter
      !
        finished=finished+1
        fmarker=fmarker+1
        qmarker=qmarker+1
        if   ( page_sweep .eq. 0) then
          sweep_path(fmarker,octant,q)=pvt
        else if ( page_sweep .eq. 1) then
          sweep_path(fmarker,1,1)=pvt
        end if
        solved(pvt)  = 1

      ! -- decrement upstream and downstream dependencies
        do j=0,3
          neigh=dadj(j,pvt);i=abs(neigh)
          if(i > 0) then
            if      ( neigh<0 .and. solved(i) .eq. 0 ) then
              allocate(current_el%next);nullify(current_el%next%next)
              current_el%cell = pvt
              current_el%face = j
              current_el => current_el%next
              neldep(octant,q)=neldep(octant,q)+1_li
              odegr(i)=odegr(i)-1
            else if ( neigh>0 .and. solved(i) .eq. 0 ) then
              idegr(i)=idegr(i)-1
            end if
          end if
        end do

      end if

    end do
  ! --- End main loop

  ! --- Copy over the eliminated dependencies into
    j=max(neldep(octant,q),1)
    if(j>1) write(6,*) '  -- Cycles detected. Number: ',j
    do g=1,egmax
      allocate( eldep(octant,q,g)%cells(j) )
      allocate( eldep(octant,q,g)%faces(j) )
      allocate( eldep(octant,q,g)%face_fluxes(num_moments_f,j) )
    end do
    i=0
    current_el=>cycle_list
    do while( associated(current_el%next) )
      i=i+1_li
      do g=1,egmax
        eldep(octant,q,g)%cells(i)=current_el%cell
        eldep(octant,q,g)%faces(i)=current_el%face
        eldep(octant,q,g)%face_fluxes(:,i)=0.0_d_t
      end do
      current_el => current_el%next
    end do

  ! --- Deallocate linked list

    current_el=>cycle_list
    do while( associated(current_el) )
       next=>current_el%next
       deallocate(current_el)
       nullify(current_el)
       current_el=>next
    end do

  ! --- If page_sweep == 1 then write path to disk

    if(page_sweep.eq.1) then
         write(99,rec=8*(q-1)+octant) sweep_path(:,1,1)
    end if
    101 FORMAT(1X,I12)

  end subroutine

!------------------------------------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------------------------------------!

  subroutine classify_reflective

    ! local variables

      integer(kind=li) :: f,face,c,alloc_stat
      type(vector)     :: normal
      real(kind=d_t)   :: d1,d2,d3,magn
      real(kind=d_t)   :: x1,x2,x3
      real(kind=d_t)   :: tol=1e-6_d_t

    ! allocate refl_face_tpe

      allocate( refl_face_tpe(rside_cells) , stat=alloc_stat)
      if(alloc_stat /=0) call stop_thor(2_li)

    ! run through all boundary faces and check their type

      do face=1,rside_cells
        c = rb_cells(face)%cell
        f = rb_cells(face)%face
        normal=outward_normal(c,f)
        magn=sqrt(normal%x1**2+normal%x2**2+normal%x3**2)
        x1=abs( normal%x1 ) / magn
        x2=abs( normal%x2 ) / magn
        x3=abs( normal%x3 ) / magn
        d1 = sqrt( (x1-one)**2 + (x2)**2 + (x3)**2 )
        d2 = sqrt( (x1)**2 + (x2-one)**2 + (x3)**2 )
        d3 = sqrt( (x1)**2 + (x2)**2 + (x3-one)**2 )
        if      (d1 < tol .and. normal%x1 > zero) then
           refl_face_tpe(face)=1_li
        else if (d1 < tol .and. normal%x1 < zero) then
           refl_face_tpe(face)=-1_li
        else if (d2 < tol .and. normal%x2 > zero) then
           refl_face_tpe(face)=2_li
        else if (d2 < tol .and. normal%x2 < zero) then
           refl_face_tpe(face)=-2_li
        else if (d3 < tol .and. normal%x3 > zero) then
           refl_face_tpe(face)=3_li
        else if (d3 < tol .and. normal%x3 < zero) then
           refl_face_tpe(face)=-3_li
        else
           call stop_thor(14_li)
        end if
      end do

  end subroutine

!------------------------------------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------------------------------------!

  subroutine check_reflective
  !*********************************************************************
  ! This subroutine ensures that reflective boundary faces are perpendicular
  ! to one of the coordinate axis
  !*********************************************************************

  ! Local variables

    integer(kind=li) :: face,f,c
    real(kind=d_t)   :: x1,x2,x3
    real(kind=d_t)   :: d1,d2,d3
    real(kind=d_t)   :: magn
    real(kind=d_t)   :: tol=1e-6_d_t
    type(vector)     :: normal

    do face=1,rside_cells
       c = rb_cells(face)%cell
       f = rb_cells(face)%face
       normal=outward_normal(c,f)
       magn=sqrt(normal%x1**2+normal%x2**2+normal%x3**2)
       x1=abs( normal%x1 ) / magn
       x2=abs( normal%x2 ) / magn
       x3=abs( normal%x3 ) / magn
       d1 = sqrt( (x1-one)**2 + (x2)**2 + (x3)**2 )
       d2 = sqrt( (x1)**2 + (x2-one)**2 + (x3)**2 )
       d3 = sqrt( (x1)**2 + (x2)**2 + (x3-one)**2 )
       if(d1 > tol .and. d2 > tol .and. d3 > tol) then
          write(6,101) 'Reflective boundary conditions on boundary face #, cell # ',c,f
          write(6,*)   'located on a non-permissible boundary face.'
          call stop_thor(-1_li)
          101 FORMAT(1X,A,2I8)
       end if
    end do

  end subroutine

!------------------------------------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------------------------------------!

  subroutine generate_normals(v0,v1,v2,v3,n0,n1,n2,n3)
  !*********************************************************************
  !
  ! Subroutine generate normals creates normal vectors for later use
  !
  !*********************************************************************
  ! Pass cell vertices

    type(vector), intent(in)  :: v0, v1, v2, v3
    type(vector), intent(out) :: n0, n1, n2, n3

  ! Local variables

    type(vector) :: dmy

  ! Compute vector normal to each face

    n0=(v2-v1) .cross. (v3-v1)
    n1=(v3-v0) .cross. (v2-v0)
    n2=(v1-v0) .cross. (v3-v0)
    n3=(v2-v0) .cross. (v1-v0)

  ! make sure the normal is an outward normal: compute face midpoint, then
  ! compute vector from point outside of face to midpoint and dot with outer
  ! normal. Results must be positive.

    ! a. n0
    dmy%x1 = third * (v1%x1+v2%x1+v3%x1) - v0%x1
    dmy%x2 = third * (v1%x2+v2%x2+v3%x2) - v0%x2
    dmy%x3 = third * (v1%x3+v2%x3+v3%x3) - v0%x3
    if( (n0 .dot. dmy) < zero) then
      n0%x1 = -n0%x1
      n0%x2 = -n0%x2
      n0%x3 = -n0%x3
    end if
    ! b. n1
    dmy%x1 = third * (v0%x1+v2%x1+v3%x1) - v1%x1
    dmy%x2 = third * (v0%x2+v2%x2+v3%x2) - v1%x2
    dmy%x3 = third * (v0%x3+v2%x3+v3%x3) - v1%x3
    if( (n1 .dot. dmy) < zero) then
      n1%x1 = -n1%x1
      n1%x2 = -n1%x2
      n1%x3 = -n1%x3
    end if
    ! b. n2
    dmy%x1 = third * (v1%x1+v0%x1+v3%x1) - v2%x1
    dmy%x2 = third * (v1%x2+v0%x2+v3%x2) - v2%x2
    dmy%x3 = third * (v1%x3+v0%x3+v3%x3) - v2%x3
    if( (n2 .dot. dmy) < zero) then
      n2%x1 = -n2%x1
      n2%x2 = -n2%x2
      n2%x3 = -n2%x3
    end if
    ! b. n3
    dmy%x1 = third * (v1%x1+v2%x1+v0%x1) - v3%x1
    dmy%x2 = third * (v1%x2+v2%x2+v0%x2) - v3%x2
    dmy%x3 = third * (v1%x3+v2%x3+v0%x3) - v3%x3
    if( (n3 .dot. dmy) < zero) then
      n3%x1 = -n3%x1
      n3%x2 = -n3%x2
      n3%x3 = -n3%x3
    end if

  end subroutine generate_normals

!------------------------------------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------------------------------------!

  subroutine memory_estimate(mem)
  !*********************************************************************
  !
  ! Computes an estimate of the consumption of memory based on the
  ! largest arrays used in THOR
  !
  !*********************************************************************

    ! --- Return value
    real(kind=d_t) :: mem

    ! --- Local variables
    integer(kind=li) :: nvar

    ! --- Initialize mem
    mem=0.0d0

    ! --- Contributions for all types
    mem = mem +&
          4.0d0*8.0d0*real(num_cells,d_t)                                               +&  ! adj list
          8.0d0*2.0d0*real(namom,d_t)*real(num_cells,d_t)*real(egmax*num_moments_v,d_t) +&  ! flux
          8.0d0*3.0d0*real(namom,d_t)*real(num_cells,d_t)*real(egmax*num_moments_v,d_t) +&  ! flux + src
          8.0d0*real(num_cells,d_t)*real(num_moments_v,d_t)                             +&
          8.0d0*4.0d0*real(num_cells,d_t)*real(num_moments_f,d_t)                       +&  ! angular flux
          8.0d0*real(num_cells,d_t)*real(num_moments_v,d_t)                                 ! directed source

    ! --- Reflected flux
    if(page_refl.eq.0_li) then
      mem = mem +&
          8.0d0*8.0d0*real(nangle,d_t)*real(rside_cells,d_t)*real(egmax*num_moments_v,d_t)        ! reflected flux
    else
      mem = mem +&
          8.0d0*8.0d0*real(nangle,d_t)*real(rside_cells,d_t)*real(num_moments_v,d_t)              ! reflected flux
    end if

    ! --- Sweep path

    if(page_sweep.eq.0) then
       mem = mem +&
             4.0d0*8.0d0*real(nangle,d_t)*real(num_cells,d_t)
    else
       mem = mem + 4.0d0*real(num_cells,d_t)
    end if

    ! --- Go through execution modii
    if(problem.eq.0) then
       if(page_iflw.eq.0_li) then
         mem = mem +&
               8.0d0*real(nangle,d_t)*real(fside_cells,d_t)*real(egmax*num_moments_f,d_t)  ! binflx
       else
         mem = mem +&
               8.0d0*real(nangle,d_t)*real(fside_cells,d_t)*real(num_moments_f,d_t)        ! binflx
       end if
    else
       if(eig_switch.eq.0) then
          mem = mem +&
                8.0d0*2.0d0*real(num_cells,d_t)*real(num_moments_v,d_t)                  ! fiss_src
       else
          nvar = egmax*num_moments_v*namom*num_cells
          mem = mem +&
                8.0d0 * 3.0d0 * real(nvar,d_t) +&                                        ! residual,du and dflx
                8.0d0 * real( (nvar+4)*(rd_restart+2)+(rd_restart+1)*rd_restart/2,d_t)   ! work
       end if
    end if

    ! --- Divide by 1E6 to get MB
    mem = mem * 1.0d-6

  end subroutine

  subroutine assign_work

    implicit none

    integer:: i, p, k=0
    integer ::rank,mpi_err, localunit, num_p, optimal_tasks
    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_p, mpi_err)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)

    optimal_tasks = ceiling((nangle*8.0)/(num_p))

    allocate(parallel_map_g2l(8*nangle,2) , parallel_map_l2g(optimal_tasks, num_p))

    do i=1, optimal_tasks
      do p = 1, num_p
        k=k+1
        if (k .gt. 8*nangle) then
          parallel_map_l2g(i,p) = 0
          cycle
        end if
        parallel_map_g2l(k,1) = p
        parallel_map_g2l(k,2) = i
        parallel_map_l2g(i,p) = k

      end do
    end do

  end subroutine

!------------------------------------------------------------------------------------------------------------------------!
! The end
!------------------------------------------------------------------------------------------------------------------------!
end module setup_module
