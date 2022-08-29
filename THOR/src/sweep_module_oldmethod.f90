MODULE sweep_module
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

  USE cell_splitting_module

  IMPLICIT NONE

  ! Define angular flux data type

  TYPE angular_flux
    REAL(kind=d_t), DIMENSION(:,:)  , ALLOCATABLE :: vol_flux
    REAL(kind=d_t), DIMENSION(:,:,:), ALLOCATABLE :: face_flux
  END TYPE angular_flux

CONTAINS

  !> This subroutine performs a transport sweep, i.e. inversion of (Omega * nabla + sigt_g)
  !> It explicitly loops over all angular directions. For each angular direction a directionally
  !> dependent source is computed (scattering + fission + fixed source). Then you do the mesh sweep
  !> on that source. The mesh sweep is done in subroutine queue. It is called queue because on an
  !> unstructured mesh the mesh sweep path is computed using a Breadth First Search algorithm that
  !> queues up all the elements in permissible order. Note: A lot of the stuff done here is for
  !> reflective boundary conditions. These are implicit boundary conditions that need to be stored
  !> to ensure convergence. Angular fluxes, in general, are not stored in THOR!
  SUBROUTINE sweep(eg,sc_flux,src,rs,reflected_flux,LL,U,Lf,Uf)
    !*********************************************************************
    !
    ! Subroutine sweep loops over all angular directions and call
    ! transport solver
    !
    !*********************************************************************

    ! Pass group number

    INTEGER(kind=li), INTENT(in) :: eg,rs

    ! Pass scalar flux

    REAL(kind=d_t) :: sc_flux(num_moments_v,namom,num_cells)

    ! Pass source

    REAL(kind=d_t) :: src(num_moments_v,namom,num_cells)

    ! Pass reflected flux

    REAL(kind=d_t), DIMENSION(num_moments_f,rs,8,nangle) :: reflected_flux

    ! Pass pre-computed matrices

    REAL(kind=d_t), DIMENSION(num_moments_v,num_moments_v) :: LL, U
    REAL(kind=d_t), DIMENSION(num_moments_f,num_moments_f) :: Lf, Uf

    ! Declare local angular flux

    TYPE(angular_flux) :: an_flux

    ! Declare face_known array=> now pre-filled in sweep and then passed to queue

    INTEGER(kind=1), DIMENSION(0:3,num_cells) :: face_known

    ! Define temporary variables

    INTEGER(kind=li) :: alloc_stat, q, oct, octant, i, f, l,&
                        k, m, indx, nk, nf
    TYPE(vector) :: omega

    ! Define source along direction q

    REAL(kind=d_t) :: dir_src(num_moments_v,num_cells)

    ! Temporary variables

    INTEGER(kind=li) :: tpe, mate,giflx, parallel_i

    INTEGER ::rank,mpi_err, num_p, optimal_tasks


    ! Define sc_flux parallel recieve buffer & reflected buffer

    REAL(kind=d_t) :: sc_flux_buffer(num_moments_v,namom,num_cells)
    REAL(kind=d_t) :: reflected_buffer(num_moments_f, rs, 8, nangle)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, num_p, mpi_err)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)
    optimal_tasks = CEILING((nangle * 8.0) / (num_p))
    sc_flux_buffer =zero
    reflected_buffer= zero
    ! Allocate angular flux

    ALLOCATE(an_flux%vol_flux(num_moments_v,num_cells)     ,stat=alloc_stat)
    IF(alloc_stat /= 0) CALL raise_fatal_error("*** Not enough memory ***")
    ALLOCATE(an_flux%face_flux(num_moments_f,0:3,num_cells),stat=alloc_stat)
    IF(alloc_stat /= 0) CALL raise_fatal_error("*** Not enough memory ***")

    ! Begin parallel loop over quadrature octant
    DO parallel_i = 1, optimal_tasks
      k = parallel_map_l2g(parallel_i, rank+1)
      IF (k .EQ. 0) EXIT
      oct = MOD(k,8)+1
      q = CEILING(k/8.0)
      octant= oct !ordering(oct)
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

      ! Zero angular flux and directed source

      dir_src=0.0_d_t
      an_flux%vol_flux  = 0.0_d_t
      an_flux%face_flux = 0.0_d_t
      face_known = 0

      ! Set boundary faces and mark as known

      ! 1. Vacuum boundary conditions
      DO i=1,vside_cells
        k=vb_cells(i)%cell
        f=vb_cells(i)%face
        IF ( (omega .dot. outward_normal(k,f)) < 0.0_d_t )THEN
          face_known(f,k) = 1
          ! no nneed to assign 0 to face flux, already done before !!
        END IF
      END DO

      ! 2. Reflective boundary conditions

      DO i=1,rside_cells
        k=rb_cells(i)%cell
        f=rb_cells(i)%face
        IF ( (omega .dot. outward_normal(k,f)) < 0.0_d_t )THEN
          face_known(f,k) = 1
          ! reflected flux(oct) stores outflow for octant oct
          ! so we need to find the mate of oct to get the right inflow
          tpe=refl_face_tpe(i)
          IF      (tpe .EQ. 1_li .OR. tpe .EQ. -1_li) THEN
            mate = mu_mate(octant)
          ELSE IF (tpe .EQ. 2_li .OR. tpe .EQ. -2_li) THEN
            mate = eta_mate(octant)
          ELSE IF (tpe .EQ. 3_li .OR. tpe .EQ. -3_li) THEN
            mate = xi_mate(octant)
          END IF
          an_flux%face_flux(:,f,k) = reflected_flux(:,i,mate,q)
        END IF
      END DO


      ! 3. Fixed inflow boundary conditions

      IF(page_iflw.EQ.1_li) THEN
        giflx=1_li
      ELSE
        giflx=eg
      END IF

      DO i=1,fside_cells
        k=fb_cells(i)%cell
        f=fb_cells(i)%face
        IF ( (omega .dot. outward_normal(k,f)) < 0.0_d_t )THEN
          face_known(f,k) = 1
          an_flux%face_flux(:,f,k) = binflx(:,i,octant,q,giflx)
        END IF
      END DO

      ! 4. Faces that are assumed known for breaking cycles

      is_cycle=0
      DO i=1, neldep(octant,q)
        k=eldep(octant,q,eg)%cells(i)
        f=eldep(octant,q,eg)%faces(i)
        is_cycle(f,k)   = 1
        face_known(f,k) = 1
        an_flux%face_flux(:,f,k) = eldep(octant,q,eg)%face_fluxes(:,i)
      END DO

      ! Compute directed source - so far only isotropic

      DO i=1,num_cells
        ! even contributions
        DO l=0,scatt_ord
          DO m=0,l
            indx=1_li+m+(l+1_li)*l/2_li
            DO k=1,num_moments_v
              dir_src(k,i)=dir_src(k,i)+Ysh(q,octant,indx)*src(k,indx,i)
            END DO
          END DO
        END DO
        ! odd contributions
        DO l=1,scatt_ord
          DO m=1,l
            indx=neven+m+(l-1_li)*l/2_li
            DO k=1,num_moments_v
              dir_src(k,i)=dir_src(k,i)+Ysh(q,octant,indx)*src(k,indx,i)
            END DO
          END DO
        END DO
      END DO

      ! If page_sweep == 1 read sweep_path from file

      IF(page_sweep.EQ.1) THEN
        sweep_path=0
        READ(99,rec=8*(q-1)+octant) sweep_path(:,1,1)
      END IF

      ! Call queue for mesh sweep

      IF (sweep_tpe .EQ.1) THEN
        CALL queue(eg,q,octant,dir_src,an_flux,LL,U,Lf,Uf,omega,face_known)
      ELSE
        CALL raise_fatal_error("Sweep type specification not recognized. Execution terminates.")
      END IF

      ! Increment scalar flux using

      CALL angular_moments(octant,q,sc_flux,an_flux)

      ! update reflected flux array, loop through r
      DO i=1,rside_cells
        k=rb_cells(i)%cell
        f=rb_cells(i)%face
        IF ( (omega .dot. outward_normal(k,f)) > 0.0_d_t )THEN
          ! store outflow on reflective faces
          reflected_buffer(:,i,octant,q) = an_flux%face_flux(:,f,k)
        END IF
      END DO
      ! update the faces that are assumed known for cycles
      DO i=1, neldep(octant,q)
        k=eldep(octant,q,eg)%cells(i)
        f=eldep(octant,q,eg)%faces(i)
        nk=adjacency_list(k,f)%cell
        nf=adjacency_list(k,f)%face
        eldep(octant,q,eg)%face_fluxes(:,i)=an_flux%face_flux(:,nf,nk)
      END DO
    END DO
    reflected_flux=reflected_buffer
    ! Deallocation
    CALL MPI_AllREDUCE(reflected_flux, reflected_buffer, num_moments_f*rs*8*nangle, &
          MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpi_err)
    reflected_flux=reflected_buffer
    CALL MPI_AllREDUCE(sc_flux, sc_flux_buffer, num_moments_v*namom*num_cells, &
          MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpi_err)
    sc_flux = sc_flux_buffer
    DEALLOCATE(an_flux%vol_flux,an_flux%face_flux)

  END SUBROUTINE sweep

  !-----------------------------------------------------------------------------------------------!
  !-----------------------------------------------------------------------------------------------!

  !> This subroutine performs a mesh sweep. It uses the sweep_path variable to obtain
  !> the sweep order.
  SUBROUTINE queue(eg,q,octant,qm,an_flux,LL,U,Lf,Uf,omega,face_known)
    !*********************************************************************
    !
    ! Subroutine queue uses the pre-computed sweep paths to perform a
    ! mesh sweep.
    !
    !*********************************************************************

    ! Pass group and angular indices

    INTEGER(kind=li), INTENT(in) :: eg, q, octant

    ! Pass qm

    REAL(kind=d_t) :: qm(num_moments_v,num_cells)

    ! Pass angular flux

    TYPE(angular_flux), INTENT(inout) :: an_flux

    ! Pass pre-computed matrices

    REAL(kind=d_t), DIMENSION(num_moments_v,num_moments_v) :: LL, U
    REAL(kind=d_t), DIMENSION(num_moments_f,num_moments_f) :: Lf, Uf

    ! Pass logical array face known

    INTEGER(kind=1), DIMENSION(0:3,num_cells) :: face_known

    ! Pass omega

    TYPE(vector) :: omega

    ! Local variables

    INTEGER(kind=li) :: i,j, CASE
    REAL(kind=d_t) :: sigmat
    TYPE(vector) :: n0, n1, n2, n3
    INTEGER(kind=li) :: soct,sq,mat_indx

    INTEGER ::rank,mpi_err, num_p
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, num_p, mpi_err)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)

    ! Set soct and sq depending on whether page_sweep=0,1

    IF      (page_sweep .EQ. 0) THEN
      soct=octant;sq=q
    ELSE IF (page_sweep .EQ. 1) THEN
      soct=1_li;sq=1_li
    END IF

    ! Initiate loop over all cells

    DO j=1,num_cells

      i=sweep_path(j,soct,sq)

      mat_indx=material_ids(reg2mat(cells(i)%reg))
      sigmat=dens_fact(cells(i)%reg)*xs_mat(mat_indx)%sigma_t(eg)

      n0=outward_normal(i,0)
      n1=outward_normal(i,1)
      n2=outward_normal(i,2)
      n3=outward_normal(i,3)

      CALL cell_orientation(omega,n0,n1,n2,n3,CASE)

      CALL cell_splitting(sigmat,qm(:,i),an_flux%vol_flux,an_flux%face_flux,   &
            LL,U,Lf,Uf,omega,i,face_known,CASE,n0,n1,n2,n3)

    END DO


  END SUBROUTINE queue

  !-----------------------------------------------------------------------------------------------!
  !-----------------------------------------------------------------------------------------------!


  !> This subroutine orientation determines whether cell's faces are
  !> incoming or outgoing and assigns them a 'case' for the transport
  !> calculation
  SUBROUTINE cell_orientation(omega,n0,n1,n2,n3,CASE)
    !*********************************************************************
    !
    ! Subroutine cell orientation determines whether cell's faces are
    ! incoming or outgoing and assigns them a 'case' for the transport
    ! calculation
    !
    !*********************************************************************
    ! Define temporary variables

    INTEGER(kind=li) :: incoming, outgoing
    INTEGER(kind=li), INTENT(inout) :: CASE
    TYPE(vector), INTENT(inout) :: omega
    TYPE(vector), INTENT(in) :: n0, n1, n2, n3

    ! Determine which case cell corresponds with respect to current ordinate

    incoming=0
    outgoing=0

    IF((omega .dot. n0) < 0.0)THEN
      incoming=incoming+1
    ELSEIF(ABS(omega .dot. n0) .LE. 1.0D-16)THEN
    ELSE
      outgoing=outgoing+1
    END IF

    IF((omega .dot. n1) < 0.0)THEN
      incoming=incoming+1
    ELSEIF(ABS(omega .dot. n1) .LE. 1.0D-16)THEN
    ELSE
      outgoing=outgoing+1
    END IF

    IF((omega .dot. n2) < 0.0)THEN
      incoming=incoming+1
    ELSEIF(ABS(omega .dot. n2) .LE. 1.0D-16)THEN
    ELSE
      outgoing=outgoing+1
    END IF

    IF((omega .dot. n3) < 0.0)THEN
      incoming=incoming+1
    ELSEIF(ABS(omega .dot. n3) .LE. 1.0D-16)THEN
    ELSE
      outgoing=outgoing+1
    END IF

    IF(incoming == 3 .AND. outgoing == 1)THEN
      CASE=1
    ELSEIF(incoming == 1 .AND. outgoing == 3)THEN
      CASE=2
    ELSEIF(incoming == 2 .AND. outgoing == 2)THEN
      CASE=3
    ELSEIF(incoming == 2 .AND. outgoing == 1)THEN
      CASE=4
    ELSEIF(incoming == 1 .AND. outgoing == 2)THEN
      CASE=5
    ELSEIF(incoming == 1 .AND. outgoing == 1)THEN
      CASE=6
    ELSE
      CALL raise_fatal_error("Unacceptable splitting in cell")
    END IF

  END SUBROUTINE cell_orientation

  !-----------------------------------------------------------------------------------------------!
  !-----------------------------------------------------------------------------------------------!

  SUBROUTINE check_incoming(face,CASE,all_incoming_known)
    !*********************************************************************
    !
    ! Subroutine check incoming determines if cell is ready for transport
    ! calculation (all incoming angular fluxes are known)
    !
    !*********************************************************************
    ! Define geometry variables

    LOGICAL, DIMENSION(0:3), INTENT(in) :: face
    LOGICAL, INTENT(inout) :: all_incoming_known
    INTEGER(kind=li), INTENT(inout) :: CASE

    ! Define temporary variables

    INTEGER(kind=li) :: known_faces, f, l

    known_faces=0

    ! Check if incoming face moments to all computed orders for cell are known

    ! Case 1 (3 incoming, 1 outgoing)

    IF(CASE == 1)THEN
      DO f=0, 3
        l=1
        IF(face(f) .EQV. .TRUE.)THEN
          known_faces=known_faces+1
        END IF
      END DO

      IF(known_faces == 3)THEN
        all_incoming_known=.TRUE.
      ELSE
        all_incoming_known=.FALSE.
      END IF


      ! Case 2 (1 incoming, 3 outgoing)

    ELSEIF(CASE == 2)THEN
      DO f=0, 3
        l=1
        IF(face(f) .EQV. .TRUE.)THEN
          known_faces=known_faces+1
        END IF
      END DO

      IF(known_faces == 1)THEN
        all_incoming_known=.TRUE.
      ELSE
        all_incoming_known=.FALSE.
      END IF

      ! Case 3 (2 incoming, 2 outgoing)

    ELSEIF(CASE == 3)THEN
      DO f=0, 3
        l=1
        IF(face(f) .EQV. .TRUE.)THEN
          known_faces=known_faces+1
        END IF
      END DO

      IF(known_faces == 2)THEN
        all_incoming_known=.TRUE.
      ELSE
        all_incoming_known=.FALSE.
      END IF

      ! Case 4 (2 incoming, 1 outgoing)

    ELSEIF(CASE == 4)THEN
      DO f=0, 3
        l=1
        IF(face(f) .EQV. .TRUE.)THEN
          known_faces=known_faces+1
        END IF
      END DO

      IF(known_faces == 2)THEN
        all_incoming_known=.TRUE.
      ELSE
        all_incoming_known=.FALSE.
      END IF

      ! Case 5 (1 incoming, 2 outgoing)

    ELSEIF(CASE == 5)THEN
      DO f=0, 3
        l=1
        IF(face(f) .EQV. .TRUE.)THEN
          known_faces=known_faces+1
        END IF
      END DO

      IF(known_faces == 1)THEN
        all_incoming_known=.TRUE.
      ELSE
        all_incoming_known=.FALSE.
      END IF

      ! Case 6 (1 incoming, 1 outgoing)

    ELSEIF(CASE == 6)THEN
      DO f=0, 3
        l=1
        IF(face(f) .EQV. .TRUE.)THEN
          known_faces=known_faces+1
        END IF
      END DO

      IF(known_faces == 1)THEN
        all_incoming_known=.TRUE.
      ELSE
        all_incoming_known=.FALSE.
      END IF

    ELSE
      CALL raise_fatal_error("Unacceptable case from cell splitting in cell")
    END IF

  END SUBROUTINE check_incoming

  !-----------------------------------------------------------------------------------------------!
  !-----------------------------------------------------------------------------------------------!

  !> This subroutine evaluates angular flux moments phi_lm from the angular fluxes.
  !> Note: The angular flux moments are accumulated on-the-fly from angular fluxes so that
  !> angular fluxes do not need to be stored.
  SUBROUTINE angular_moments(octant,q,sc_flux,an_flux)
    !*********************************************************************
    !
    ! Subroutine angular moments computes the real spherical harmonics
    ! angular moments of the angular flux
    !
    !*********************************************************************
    ! Pass input parameters

    INTEGER(kind=li), INTENT(in) :: q,octant

    ! Declare angular and scalar fluxes

    REAL(kind=d_t) :: sc_flux(num_moments_v,namom,num_cells)
    TYPE(angular_flux), INTENT(in) :: an_flux

    ! Define temporary variables

    INTEGER(kind=li) :: i, l, k, indx, m

    ! Compute scalar flux angular moments, right now only isotropic
    ! flux is accumulated

    DO i=1, num_cells
      ! even contributions
      DO l=0,scatt_ord
        DO m=0,l
          DO k=1, num_moments_v
            indx=1_li+m+(l+1_li)*l/2_li
            sc_flux(k,indx,i)=sc_flux(k,indx,i)+&
                  (1.0_d_t/8.0_d_t)*quadrature(q)%wt*Ysh(q,octant,indx)*an_flux%vol_flux(k,i)
          END DO
        END DO
      END DO
      ! odd contributions
      DO l=1,scatt_ord
        DO m=1,l
          DO k=1, num_moments_v
            indx=neven+m+(l-1_li)*l/2_li
            sc_flux(k,indx,i)=sc_flux(k,indx,i)+&
                  (1.0_d_t/8.0_d_t)*quadrature(q)%wt*Ysh(q,octant,indx)*an_flux%vol_flux(k,i)
          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE angular_moments

END MODULE sweep_module
