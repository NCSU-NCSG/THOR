MODULE setup_module
  !***********************************************************************
  !
  ! Setup module contains all subroutines to setup problem
  !
  !***********************************************************************

  ! Use derived-type modules

  USE mpi
  USE types
  USE parameter_types
  USE filename_types
  USE vector_types
  USE cross_section_types
  USE geometry_types
  USE angle_types
  USE multindex_types
  USE globals

  ! Use modules that pertain setting up problem

  USE read_module
  USE quadrature_module
  USE sweep_module
  USE error_module

  IMPLICIT NONE

CONTAINS

  !> This subroutine reads the input and prepares data for execution
  SUBROUTINE setup
    !*********************************************************************
    !
    ! Subroutine setup calls routines to read and prepare input for
    ! execution
    !
    !*********************************************************************

    ! Local variables

    INTEGER(kind=li) :: i,alloc_stat
    INTEGER(kind=li) :: q
    TYPE(vector) :: v0, v1, v2, v3, n0, n1, n2, n3
    INTEGER(kind=li) :: reflective(3,2)
    REAL(kind=d_t) :: ts,te
    REAL(kind=d_t), DIMENSION(3,3) :: J, J_inv
    REAL(kind=d_t)                 :: det_J
    LOGICAL :: existence
    INTEGER :: rank, num_p, mpi_err
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, num_p, mpi_err)

    ! Read input

    CALL READ

    ! Assign local angles for parallel analysis

    CALL assign_work

    ! Allocate and initialize density factors and region volume

    ALLOCATE(dens_fact(minreg:maxreg));dens_fact=1.0_d_t
    ALLOCATE(reg_vol  (minreg:maxreg));reg_vol  =0.0_d_t

    DO i=1, num_cells
      v0=vertices(cells(i)%R(0))%v
      v1=vertices(cells(i)%R(1))%v
      v2=vertices(cells(i)%R(2))%v
      v3=vertices(cells(i)%R(3))%v
      CALL cell_jacobian(v0,v1,v2,v3,J,det_J,J_inv)
      cells(i)%volume=(1.0_d_t/6.0_d_t)*det_J
      reg_vol(cells(i)%reg)=reg_vol(cells(i)%reg)+cells(i)%volume
    END DO

    IF(dfact_opt.NE.0) THEN
      CALL read_density_factors
      IF(dfact_opt.EQ.1) THEN
        DO i=minreg,maxreg
          IF(dens_fact(i)<0.0_d_t .OR. reg_vol(i)<2.24d-12) THEN
            dens_fact(i)=1.0_d_t
          ELSE
            dens_fact(i)=dens_fact(i)/reg_vol(i)
          END IF
        END DO
      ELSE
        DO i=minreg,maxreg
          IF(dens_fact(i)<0.0_d_t) dens_fact(i)=1.0_d_t
        END DO
      END IF
    END IF

    CALL write_reg_info

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
    IF (eig_switch .EQ. 1) THEN
      niter=2_li
    ELSE
      niter=3_li
    END IF

    ! Pre-compute the cell outward normal vectors

    ALLOCATE(outward_normal(num_cells,0:3),stat=alloc_stat)
    IF(alloc_stat /=0) CALL raise_fatal_error("*** Not enough memory ***")

    ! Generate outward normal vectors for all cells for performance

    DO i=1, num_cells
      v0=vertices(cells(i)%R(0))%v
      v1=vertices(cells(i)%R(1))%v
      v2=vertices(cells(i)%R(2))%v
      v3=vertices(cells(i)%R(3))%v
      CALL generate_normals(v0,v1,v2,v3,n0,n1,n2,n3)
      outward_normal(i,0)=n0
      outward_normal(i,1)=n1
      outward_normal(i,2)=n2
      outward_normal(i,3)=n3
    END DO

    IF (rside_cells > 0) THEN

      ! Check reflective boundary faces: They need to be on a face that is
      ! perpendicuar to either x,y or z-axis.

      CALL check_reflective

      ! Classify reflective boundary faces by axis they are orthogonal to:
      ! +1 = outward normal along positive x axis, -1 = outward normal along negative x-axis
      ! +2 = outward normal along positive y axis, -2 = outward normal along negative y-axis
      ! +3 = outward normal along positive z axis, -3 = outward normal along negative z-axis

      CALL classify_reflective

    END IF

    ! Determine which boundary faces are reflective and save in reflective

    reflective=0_li
    DO i=1,rside_cells
      q=ABS(refl_face_tpe(i))
      IF(refl_face_tpe(i)>0) THEN
        reflective(q,2)=1
      ELSE
        reflective(q,1)=1
      END IF
    END DO

    ! Determine ordering of octants

    CALL set_ordering(reflective)

    ! Make sure that if JFNK is selected reflective BC in opposite faces is not
    ! permitted
    IF(problem .EQ. 1 .AND. eig_switch .EQ. 1 .AND. rd_method > 1) THEN
      DO i = 1, 3
        IF(reflective(i,1) .EQ. 1 .AND. reflective(i,2).EQ.1 ) CALL raise_fatal_error("JFNK option &
          & in combination with reflective boundary conditions on opposite faces is not permitted")
      END DO
    END IF

    ! Error out if we have reflecting BCs, JFNK, and nproc > 1
    IF (eig_switch .EQ. 1 .and. ABS(SUM(reflective)) .GT. 1.0D-16 .and. num_p > 1) THEN
      CALL raise_fatal_error( "Parallel execution of JFNK with reflecting BC is currently not supported in THOR.")
      !TODO: This should be analyzed and figured out if we can do something different, or if we can just set it to 1 thread
    END IF

101 FORMAT(A,ES12.4)

    ! Set up the variables to memorize rings in the mesh sweep

    ALLOCATE(neldep(8,nangle),stat=alloc_stat)
    IF(alloc_stat /=0) CALL raise_fatal_error('*** Not enough memory ***')
    neldep=0_li
    ALLOCATE(eldep(8,nangle,egmax) ,stat=alloc_stat)
    IF(alloc_stat /=0) CALL raise_fatal_error('*** Not enough memory ***')
    ALLOCATE( is_cycle(0:3,num_cells) )

    ! Pre-compute sweep if sweep_tpe = 1

    IF (sweep_tpe .EQ. 1) THEN

      ! allocate sweep path
      IF      (page_sweep .EQ. 0) THEN
        ALLOCATE( sweep_path(num_cells,8,nangle),stat=alloc_stat)
        IF(alloc_stat /=0) CALL raise_fatal_error('*** Not enough memory ***')
      ELSE IF (page_sweep .EQ. 1) THEN
        ALLOCATE( sweep_path(num_cells,1,1),stat=alloc_stat)
        IF(alloc_stat /=0) CALL raise_fatal_error('*** Not enough memory ***')
        ! open scratch file that contains sweep path
        INQUIRE(file="sweep_path.pg",exist=existence)
        IF( existence .EQV. .TRUE.) THEN
          OPEN(unit=99,file='sweep_path.pg',status='replace',form='unformatted',access='direct',recl=4*num_cells)
        ELSE
          OPEN(unit=99,file='sweep_path.pg',status='new'    ,form='unformatted',access='direct',recl=4*num_cells)
        END IF
      END IF
      sweep_path=0_li

      ! pre-compute sweep
      IF (rank .EQ. 0) THEN
        CALL printlog("--------------------------------------------------------")
        CALL printlog("   Precomputing sweep path   ")
        CALL printlog("--------------------------------------------------------")
      END IF
      !FIX ME - pre sweep efficiency and threading
      CALL CPU_TIME(ts)
      CALL pre_sweep
      CALL CPU_TIME(te)
      IF (rank .EQ. 0) THEN
        WRITE(amsg,101)'-- Finished precomputing sweep path. Time (sec.) ', te-ts
        CALL printlog(amsg)
        CALL printlog('')
      END IF
    END IF

    ! set variables that should eventually be moved to input

    exmax = 0.5_d_t
    extol = 0.05_d_t

  END SUBROUTINE setup

  !------------------------------------------------------------------------------------------------------------------------!
  !------------------------------------------------------------------------------------------------------------------------!

  SUBROUTINE write_reg_info
    !*********************************************************************
    !
    ! Write region information
    !
    !*********************************************************************

    ! Local variables

    INTEGER(kind=li) :: i

    ! Write region information
    IF (rank .EQ. 0) THEN
      CALL printlog('')
      CALL printlog("--------------------------------------------------------")
      CALL printlog("   Region Information  ")
      CALL printlog("--------------------------------------------------------")
      CALL printlog('')
      CALL printlog("     Region ID   Material ID  Region Volume Density Factor  Material Name")
      DO i = minreg,maxreg
        WRITE(amsg,104) i,reg2mat(i),reg_vol(i),dens_fact(i),'  ',TRIM(xs_mat(i)%mat_name)
        CALL printlog(amsg)
      END DO
    END IF

104 FORMAT(I14,I14,2ES15.4,2A)

  END SUBROUTINE write_reg_info

  !------------------------------------------------------------------------------------------------------------------------!
  !------------------------------------------------------------------------------------------------------------------------!

  SUBROUTINE read_density_factors
    !*********************************************************************
    !
    ! Reads density factor file
    !
    !*********************************************************************

    ! Local variables
    LOGICAL :: ex
    INTEGER(kind=li) :: i

    ! Open and read file, then close, terminate THOR if file does not exist

    INQUIRE(file=dens_fact_filename,exist=ex)
    IF(ex) THEN
      OPEN(file=dens_fact_filename,unit=51)
      READ(51,*) (dens_fact(i),i=minreg,maxreg)
      CLOSE(unit=51)
    ELSE
      CALL raise_fatal_error("Density factors were requested but referenced file was not found.")
      RETURN
    END IF
  END SUBROUTINE read_density_factors
  !------------------------------------------------------------------------------------------------------------------------!
  !------------------------------------------------------------------------------------------------------------------------!

  SUBROUTINE set_ordering(reflective)
    !*********************************************************************
    !
    ! Subroutine computes ordering such that sweep starts from vacuum BC
    !
    !*********************************************************************

    ! --- Arguments
    INTEGER(kind=li) :: reflective(3,2)

    ! --- Local variables
    INTEGER(kind=li) :: xd = 0, yd = 0, zd=0
    INTEGER(kind=li) :: parallel_i, k, oct, permutation(8), q
    INTEGER :: rank, num_p, mpi_err

    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, num_p, mpi_err)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)

    ! --- Set ordering
    !  -- x-dir
    IF      (reflective(1,1) .EQ. 0_li .AND. reflective(1,2) .EQ. 0_li) THEN
      xd = 1_li
    ELSE IF (reflective(1,1) .EQ. 1_li .AND. reflective(1,2) .EQ. 0_li) THEN
      xd = -1_li
    ELSE IF (reflective(1,1) .EQ. 0_li .AND. reflective(1,2) .EQ. 1_li) THEN
      xd = 1_li
    ELSE IF (reflective(1,1) .EQ. 1_li .AND. reflective(1,2) .EQ. 1_li) THEN
      xd = 1_li
    END IF
    !  -- y-dir
    IF      (reflective(2,1) .EQ. 0_li .AND. reflective(2,2) .EQ. 0_li) THEN
      yd = 1_li
    ELSE IF (reflective(2,1) .EQ. 1_li .AND. reflective(2,2) .EQ. 0_li) THEN
      yd = -1_li
    ELSE IF (reflective(2,1) .EQ. 0_li .AND. reflective(2,2) .EQ. 1_li) THEN
      yd = 1_li
    ELSE IF (reflective(2,1) .EQ. 1_li .AND. reflective(2,2) .EQ. 1_li) THEN
      yd = 1_li
    END IF
    !  -- z-dir
    IF      (reflective(3,1) .EQ. 0_li .AND. reflective(3,2) .EQ. 0_li) THEN
      zd = 1_li
    ELSE IF (reflective(3,1) .EQ. 1_li .AND. reflective(3,2) .EQ. 0_li) THEN
      zd = -1_li
    ELSE IF (reflective(3,1) .EQ. 0_li .AND. reflective(3,2) .EQ. 1_li) THEN
      zd = 1_li
    ELSE IF (reflective(3,1) .EQ. 1_li .AND. reflective(3,2) .EQ. 1_li) THEN
      zd = 1_li
    END IF

    !  -- Chose octant ordering
    ! TODO: check the ordering...might not be right
    IF      (xd .EQ.  1_li .AND. yd .EQ.  1_li .AND. zd .EQ.  1_li) THEN
      ordering = (/1,2,3,4,5,6,7,8/)
    ELSE IF (xd .EQ.  1_li .AND. yd .EQ.  1_li .AND. zd .EQ. -1_li) THEN
      ordering = (/5,6,7,8,1,2,3,4/)
    ELSE IF (xd .EQ.  1_li .AND. yd .EQ. -1_li .AND. zd .EQ.  1_li) THEN
      ordering = (/4,3,1,8,7,5,2,6/)
    ELSE IF (xd .EQ. -1_li .AND. yd .EQ.  1_li .AND. zd .EQ.  1_li) THEN
      ordering = (/2,3,6,7,1,4,5,8/)
    ELSE IF (xd .EQ.  1_li .AND. yd .EQ. -1_li .AND. zd .EQ. -1_li) THEN
      ordering = (/8,7,6,5,4,3,2,1/)
    ELSE IF (xd .EQ. -1_li .AND. yd .EQ.  1_li .AND. zd .EQ. -1_li) THEN
      ordering = (/6,7,8,5,2,3,4,1/)
    ELSE IF (xd .EQ. -1_li .AND. yd .EQ. -1_li .AND. zd .EQ.  1_li) THEN
      ordering = (/3,4,1,2,7,8,5,6/)
    ELSE IF (xd .EQ. -1_li .AND. yd .EQ. -1_li .AND. zd .EQ. -1_li) THEN
      ordering = (/7,6,8,5,3,2,4,1/)
    END IF

    ! get the ordering of oct in the loop implemented in sweep_module
    octants_to_sweep = 9_li
    permutation = (/1, 2, 3, 4, 5, 6, 7, 8/)
    q = one
    DO parallel_i = 1, CEILING((nangle * 8.0) / (num_p))
      k = parallel_map_l2g(parallel_i, rank + 1)
      IF (k .EQ. 0) EXIT
      oct = MOD(k, 8) + 1
      IF (indexOf(oct, octants_to_sweep) == -1_li) THEN
        octants_to_sweep(q) = oct
        q = q + 1
      END IF
    END DO
    ordered_octants_to_sweep = octants_to_sweep
    CALL quickSortInteger(ordered_octants_to_sweep, permutation)
  END SUBROUTINE set_ordering

  !------------------------------------------------------------------------------------------------------------------------!
  !------------------------------------------------------------------------------------------------------------------------!

  !> This subroutine pre-computes the sweep path and saves it in
  !> variable sweep_path
  SUBROUTINE pre_sweep
    !*********************************************************************
    !
    ! Subroutine pre-sweep pre-computes the sweep path and saves it in
    ! variable sweep_path
    !
    !*********************************************************************

    ! Define temporary variables

    INTEGER(kind=li) :: q, oct, octant,&
                        i, k, index
    TYPE(vector) :: omega

    ! Temporary variables

    REAL(kind=d_t) :: te,ts

    INTEGER ::rank,mpi_err, num_p, optimal_tasks
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, num_p, mpi_err)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)

    ! Begin parallel loop over quadrature octant
    optimal_tasks = CEILING((nangle*8.0)/(num_p))

    DO i=1, optimal_tasks

      k = parallel_map_l2g(i, rank+1)
      IF (k .EQ. 0) EXIT
      oct = MOD(k,8)+1
      q = CEILING(k/8.0)

      index = indexOf(oct, octants_to_sweep)
      IF (index > 8_li .or. index < 1) CALL raise_fatal_error( "Sweep order index out of bounds")
      octant = ordering(ordered_octants_to_sweep(index))

      IF(octant == 1)THEN
        omega%x1=quadrature(q)%mu%x1
        omega%x2=quadrature(q)%mu%x2
        omega%x3=quadrature(q)%mu%x3
      ELSEIF(octant == 2)THEN
        omega%x1=-1*quadrature(q)%mu%x1
        omega%x2=quadrature(q)%mu%x2
        omega%x3=quadrature(q)%mu%x3
      ELSEIF(octant == 3)THEN
        omega%x1=-1*quadrature(q)%mu%x1
        omega%x2=-1*quadrature(q)%mu%x2
        omega%x3=quadrature(q)%mu%x3
      ELSEIF(octant == 4)THEN
        omega%x1=quadrature(q)%mu%x1
        omega%x2=-1*quadrature(q)%mu%x2
        omega%x3=quadrature(q)%mu%x3
      ELSEIF(octant == 5)THEN
        omega%x1=quadrature(q)%mu%x1
        omega%x2=quadrature(q)%mu%x2
        omega%x3=-1*quadrature(q)%mu%x3
      ELSEIF(octant == 6)THEN
        omega%x1=-1*quadrature(q)%mu%x1
        omega%x2=quadrature(q)%mu%x2
        omega%x3=-1*quadrature(q)%mu%x3
      ELSEIF(octant == 7)THEN
        omega%x1=-1*quadrature(q)%mu%x1
        omega%x2=-1*quadrature(q)%mu%x2
        omega%x3=-1*quadrature(q)%mu%x3
      ELSE
        omega%x1=quadrature(q)%mu%x1
        omega%x2=-1*quadrature(q)%mu%x2
        omega%x3=-1*quadrature(q)%mu%x3
      END IF

      ! Start timer

      CALL CPU_TIME(ts)

      CALL pre_queue(q,octant,omega)

      ! End timer and print message

      CALL CPU_TIME(te)
      !FIX ME - Only prints angles belonging to root process
      IF (rank .EQ. 0) THEN
        WRITE(amsg,101) ' - Computation of sweep path for octant: ',octant,' angle: ',q,' . Ex. Time(sec.): ',te-ts
        CALL printlog(amsg)
      END IF


    END DO

101 FORMAT(A,I4,A,I4,A,ES12.4)
  END SUBROUTINE pre_sweep

  !------------------------------------------------------------------------------------------------------------------------!
  !------------------------------------------------------------------------------------------------------------------------!

  !> This subroutine performs a BFS of the mesh along direction omega.
  !> It deals with rings by severing the connection across one face
  !> within the ring and saving the cell/face #.
  SUBROUTINE pre_queue(q,octant,omega)
    !*********************************************************************
    !
    ! This subroutine performs a BFS of the mesh along direction omega.
    ! It deals with rings by severing the connection across one face
    ! within the ring and saving the cell/face #.
    !
    !*********************************************************************

    ! --- Arguments
    INTEGER(kind=li) :: q,octant
    TYPE(vector) :: omega

    ! --- Local variables
    INTEGER(kind=li) :: queue(num_cells)
    INTEGER(kind=1)  :: idegr(num_cells)
    INTEGER(kind=1)  :: odegr(num_cells)
    INTEGER(kind=li) :: qmarker,iqmarker,fmarker,ifmarker,finished
    INTEGER(kind=li) :: g,i,j,neigh,pvt=0,cell,face
    INTEGER(kind=1)  :: solved(num_cells)
    INTEGER(kind=li) :: dadj(0:3,num_cells)
    TYPE(vector)     :: onor
    TYPE(upstream_face_le),POINTER :: cycle_list,current_el,next

    ! --- Zero sweep_path if page_sweep==1
    IF(page_sweep.EQ.1) sweep_path(:,1,1)=0

    ! --- Prepare cycle list
    NULLIFY(cycle_list);ALLOCATE(cycle_list)
    NULLIFY(cycle_list%next)
    current_el => cycle_list

    ! --- Constructing the directed adjacency list
    !-------------------------------
    ! Problem w/ omega.dot.onor ==0 not resolved
    !-------------------------------
    dadj=0_li
    DO cell=1,num_cells
      DO face=0,3
        ! -- get current normal
        onor=outward_normal(cell,face)
        ! -- decide if inflow or outflow
        IF      ( (omega .dot. onor) > 0.0_d_t) THEN   ! - outflow
          dadj(face,cell)=adjacency_list(cell,face)%cell
        ELSE IF ( (omega .dot. onor) < 0.0_d_t) THEN   ! - inflow
          dadj(face,cell)=-adjacency_list(cell,face)%cell
        END IF
      END DO
    END DO

    ! --- Start algorithm by setting indegree and outdegree:
    !     indegree : # upstream dependencies
    !     outdegree: # downstream dependencies

    idegr=0 ; odegr=0

    DO cell=1,num_cells
      DO face=0,3
        IF      (dadj(face,cell)>0) THEN
          odegr(cell)=odegr(cell)+INT(1,1)
        ELSE IF (dadj(face,cell)<0) THEN
          idegr(cell)=idegr(cell)+INT(1,1)
        END IF
      END DO
    END DO

    ! --- Start main loop
    finished=0
    fmarker =0
    ifmarker=num_cells+1
    qmarker =0
    iqmarker=num_cells+1
    solved  =0

    DO WHILE(finished<num_cells)

      ! - Fill queue with all tets that initially have no upstream
      !   dependencies
      DO i=1,num_cells
        IF(idegr(i) .EQ. 0 .AND. solved(i) .EQ. 0) THEN
          qmarker=qmarker+1
          queue(qmarker)=i
        END IF
      END DO

      ! -- Downstream trim
      DO WHILE(fmarker .LT. qmarker)

        ! - Get current tet
        finished=finished+1
        fmarker = fmarker+1
        cell=queue(fmarker)

        ! - list cell in order
        IF      ( page_sweep .EQ. 0) THEN
          sweep_path(fmarker,octant,q)=cell
        ELSE IF ( page_sweep .EQ. 1) THEN
          sweep_path(fmarker,1,1)=cell
        END IF

        solved(cell)  =1
        ! - Decrements the upstream dependencies of the downstream neighbors
        DO j=0,3
          neigh=dadj(j,cell)
          IF(neigh>0) THEN
            idegr(neigh)=idegr(neigh)-INT(1,1)
            IF(idegr(neigh).EQ.0 .AND. solved(neigh) .EQ. 0) THEN
              qmarker=qmarker+1
              queue(qmarker)=neigh
            END IF
          END IF
        END DO

      END DO

      ! -- If downstream trim did not touch all cells, then there are cycles
      IF(finished<num_cells) THEN

        ! - Fill queue reversely with all tets that initially have no downstream
        !   dependencies
        DO i=1,num_cells
          IF(odegr(i) .EQ. 0 .AND. solved(i) .EQ. 0) THEN
            iqmarker=iqmarker-1
            queue(iqmarker)=i
          END IF
        END DO

        ! -- Upstream trim
        DO WHILE(ifmarker .GT. iqmarker)

          ! - Get current tet
          finished=finished+1
          ifmarker = ifmarker-1
          cell=queue(ifmarker)

          ! - list cell in order
          IF   ( page_sweep .EQ. 0) THEN
            sweep_path(ifmarker,octant,q)=cell
          ELSE IF ( page_sweep .EQ. 1) THEN
            sweep_path(ifmarker,1,1)=cell
          END IF

          solved(cell)  =1
          ! - Decrements the upstream dependencies of the downstream neighbors
          DO j=0,3
            neigh=dadj(j,cell)
            IF(neigh<0) THEN
              odegr(ABS(neigh))=odegr(ABS(neigh))-INT(1,1)
              IF(odegr(ABS(neigh)).EQ.0 .AND. solved(ABS(neigh)) .EQ. 0) THEN
                iqmarker=iqmarker-1
                queue(iqmarker)=ABS(neigh)
              END IF
            END IF
          END DO

        END DO

      END IF

      IF(finished<num_cells) THEN

        ! -- Find element with only one upstream dependence or use
        !    and element with two upstream dependencies
        j=3;i=0
        DO WHILE(j>1 .AND. i.LT.num_cells)
          i=i+1
          IF(solved(i) .EQ. 0) THEN
            IF(idegr(i)<j) THEN
              j=idegr(i)
              pvt=i
            END IF
          END IF
        END DO

        ! -- List pvt in order and increment finished counter
        !
        finished=finished+1
        fmarker=fmarker+1
        qmarker=qmarker+1
        IF   ( page_sweep .EQ. 0) THEN
          sweep_path(fmarker,octant,q)=pvt
        ELSE IF ( page_sweep .EQ. 1) THEN
          sweep_path(fmarker,1,1)=pvt
        END IF
        solved(pvt)  = 1

        ! -- decrement upstream and downstream dependencies
        DO j=0,3
          neigh=dadj(j,pvt);i=ABS(neigh)
          IF(i > 0) THEN
            IF      ( neigh<0 .AND. solved(i) .EQ. 0 ) THEN
              ALLOCATE(current_el%next);NULLIFY(current_el%next%next)
              current_el%cell = pvt
              current_el%face = j
              current_el => current_el%next
              neldep(octant,q)=neldep(octant,q)+1_li
              odegr(i)=odegr(i)-INT(1,1)
            ELSE IF ( neigh>0 .AND. solved(i) .EQ. 0 ) THEN
              idegr(i)=idegr(i)-INT(1,1)
            END IF
          END IF
        END DO

      END IF

    END DO
    ! --- End main loop

    ! --- Copy over the eliminated dependencies into
    j=MAX(neldep(octant,q),1)
    IF(j>1) WRITE(amsg,'(A,I0)') '  -- Cycles detected. Number: ',j
    IF(j>1) CALL printlog(amsg)
    DO g=1,egmax
      ALLOCATE( eldep(octant,q,g)%cells(j) )
      ALLOCATE( eldep(octant,q,g)%faces(j) )
      ALLOCATE( eldep(octant,q,g)%face_fluxes(num_moments_f,j) )
    END DO
    i=0
    current_el=>cycle_list
    DO WHILE( ASSOCIATED(current_el%next) )
      i=i+1_li
      DO g=1,egmax
        eldep(octant,q,g)%cells(i)=current_el%cell
        eldep(octant,q,g)%faces(i)=current_el%face
        eldep(octant,q,g)%face_fluxes(:,i)=0.0_d_t
      END DO
      current_el => current_el%next
    END DO

    ! --- Deallocate linked list

    current_el=>cycle_list
    DO WHILE( ASSOCIATED(current_el) )
      next=>current_el%next
      DEALLOCATE(current_el)
      NULLIFY(current_el)
      current_el=>next
    END DO

    ! --- If page_sweep == 1 then write path to disk

    IF(page_sweep.EQ.1) THEN
      WRITE(99,rec=8*(q-1)+octant) sweep_path(:,1,1)
    END IF

  END SUBROUTINE pre_queue

  !------------------------------------------------------------------------------------------------------------------------!
  !------------------------------------------------------------------------------------------------------------------------!

  SUBROUTINE classify_reflective

    ! local variables

    INTEGER(kind=li) :: f,face,c,alloc_stat
    TYPE(vector)     :: normal
    REAL(kind=d_t)   :: d1,d2,d3,magn
    REAL(kind=d_t)   :: x1,x2,x3
    REAL(kind=d_t)   :: tol=1e-6_d_t

    ! allocate refl_face_tpe

    ALLOCATE( refl_face_tpe(rside_cells) , stat=alloc_stat)
    IF(alloc_stat /=0) CALL raise_fatal_error("*** Not enough memory ***")

    ! run through all boundary faces and check their type

    DO face=1,rside_cells
      c = rb_cells(face)%cell
      f = rb_cells(face)%face
      normal=outward_normal(c,f)
      magn=SQRT(normal%x1**2+normal%x2**2+normal%x3**2)
      x1=ABS( normal%x1 ) / magn
      x2=ABS( normal%x2 ) / magn
      x3=ABS( normal%x3 ) / magn
      d1 = SQRT( (x1-one)**2 + (x2)**2 + (x3)**2 )
      d2 = SQRT( (x1)**2 + (x2-one)**2 + (x3)**2 )
      d3 = SQRT( (x1)**2 + (x2)**2 + (x3-one)**2 )
      IF      (d1 < tol .AND. normal%x1 > zero) THEN
        refl_face_tpe(face)=1_li
      ELSE IF (d1 < tol .AND. normal%x1 < zero) THEN
        refl_face_tpe(face)=-1_li
      ELSE IF (d2 < tol .AND. normal%x2 > zero) THEN
        refl_face_tpe(face)=2_li
      ELSE IF (d2 < tol .AND. normal%x2 < zero) THEN
        refl_face_tpe(face)=-2_li
      ELSE IF (d3 < tol .AND. normal%x3 > zero) THEN
        refl_face_tpe(face)=3_li
      ELSE IF (d3 < tol .AND. normal%x3 < zero) THEN
        refl_face_tpe(face)=-3_li
      ELSE
        CALL raise_fatal_error("Reflective boundary face could not be determined.")
      END IF
    END DO

  END SUBROUTINE classify_reflective

  !------------------------------------------------------------------------------------------------------------------------!
  !------------------------------------------------------------------------------------------------------------------------!

  SUBROUTINE check_reflective
    !*********************************************************************
    ! This subroutine ensures that reflective boundary faces are perpendicular
    ! to one of the coordinate axis
    !*********************************************************************

    ! Local variables

    INTEGER(kind=li) :: face,f,c
    REAL(kind=d_t)   :: x1,x2,x3
    REAL(kind=d_t)   :: d1,d2,d3
    REAL(kind=d_t)   :: magn
    REAL(kind=d_t)   :: tol=1e-6_d_t
    TYPE(vector)     :: normal

    DO face=1,rside_cells
      c = rb_cells(face)%cell
      f = rb_cells(face)%face
      normal=outward_normal(c,f)
      magn=SQRT(normal%x1**2+normal%x2**2+normal%x3**2)
      x1=ABS( normal%x1 ) / magn
      x2=ABS( normal%x2 ) / magn
      x3=ABS( normal%x3 ) / magn
      d1 = SQRT( (x1-one)**2 + (x2)**2 + (x3)**2 )
      d2 = SQRT( (x1)**2 + (x2-one)**2 + (x3)**2 )
      d3 = SQRT( (x1)**2 + (x2)**2 + (x3-one)**2 )
      IF(d1 > tol .AND. d2 > tol .AND. d3 > tol) THEN
        CALL raise_fatal_error('Reflective boundary conditions on boundary face #, cell # ' &
          //TRIM(STR(c))//', '//TRIM(STR(f))//' located on a non-permissible boundary face.')
      END IF
    END DO

  END SUBROUTINE check_reflective

  !------------------------------------------------------------------------------------------------------------------------!
  !------------------------------------------------------------------------------------------------------------------------!

  SUBROUTINE generate_normals(v0,v1,v2,v3,n0,n1,n2,n3)
    !*********************************************************************
    !
    ! Subroutine generate normals creates normal vectors for later use
    !
    !*********************************************************************
    ! Pass cell vertices

    TYPE(vector), INTENT(in)  :: v0, v1, v2, v3
    TYPE(vector), INTENT(out) :: n0, n1, n2, n3

    ! Local variables

    TYPE(vector) :: dmy

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
    IF( (n0 .dot. dmy) < zero) THEN
      n0%x1 = -n0%x1
      n0%x2 = -n0%x2
      n0%x3 = -n0%x3
    END IF
    ! b. n1
    dmy%x1 = third * (v0%x1+v2%x1+v3%x1) - v1%x1
    dmy%x2 = third * (v0%x2+v2%x2+v3%x2) - v1%x2
    dmy%x3 = third * (v0%x3+v2%x3+v3%x3) - v1%x3
    IF( (n1 .dot. dmy) < zero) THEN
      n1%x1 = -n1%x1
      n1%x2 = -n1%x2
      n1%x3 = -n1%x3
    END IF
    ! b. n2
    dmy%x1 = third * (v1%x1+v0%x1+v3%x1) - v2%x1
    dmy%x2 = third * (v1%x2+v0%x2+v3%x2) - v2%x2
    dmy%x3 = third * (v1%x3+v0%x3+v3%x3) - v2%x3
    IF( (n2 .dot. dmy) < zero) THEN
      n2%x1 = -n2%x1
      n2%x2 = -n2%x2
      n2%x3 = -n2%x3
    END IF
    ! b. n3
    dmy%x1 = third * (v1%x1+v2%x1+v0%x1) - v3%x1
    dmy%x2 = third * (v1%x2+v2%x2+v0%x2) - v3%x2
    dmy%x3 = third * (v1%x3+v2%x3+v0%x3) - v3%x3
    IF( (n3 .dot. dmy) < zero) THEN
      n3%x1 = -n3%x1
      n3%x2 = -n3%x2
      n3%x3 = -n3%x3
    END IF

  END SUBROUTINE generate_normals

  !------------------------------------------------------------------------------------------------------------------------!
  !------------------------------------------------------------------------------------------------------------------------!

  SUBROUTINE assign_work

    IMPLICIT NONE

    INTEGER:: i, p, k=0
    INTEGER ::rank,mpi_err, num_p, optimal_tasks
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, num_p, mpi_err)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)

    optimal_tasks = CEILING((nangle*8.0)/(num_p))

    ALLOCATE(parallel_map_g2l(8*nangle,2) , parallel_map_l2g(optimal_tasks, num_p))

    DO i=1, optimal_tasks
      DO p = 1, num_p
        k=k+1
        IF (k .GT. 8*nangle) THEN
          parallel_map_l2g(i,p) = 0
          CYCLE
        END IF
        parallel_map_g2l(k,1) = p
        parallel_map_g2l(k,2) = i
        parallel_map_l2g(i,p) = k

      END DO
    END DO

  END SUBROUTINE assign_work

  !------------------------------------------------------------------------------------------------------------------------!
  ! The end
  !------------------------------------------------------------------------------------------------------------------------!
END MODULE setup_module
