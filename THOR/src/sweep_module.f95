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
  USE global_variables
  USE SDD_global_variables

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

    INTEGER(kind=li) :: alloc_stat, q, oct, octant, i, f, l, rcell,&
                        j, k, m, indx, nk, nf
    TYPE(vector) :: omega

    ! Define source along direction q

    REAL(kind=d_t) :: dir_src(num_moments_v,num_cells), reallocate_start,reallocate_end

    ! Temporary variables

    INTEGER(kind=li) :: tpe, mate,giflx, parallel_i,kpsi_counter

    INTEGER ::rank,mpi_err, localunit, num_p, optimal_tasks, ind, calcK


    ! Define sc_flux parallel recieve buffer & reflected buffer

    REAL(kind=d_t) :: sc_flux_buffer(num_moments_v,namom,num_cells)
    REAL(kind=d_t) :: reflected_buffer(num_moments_f, rs, 8, nangle)

    !! ADD REMOVED - OCT 2019
    !CALL MPI_COMM_SIZE(MPI_COMM_WORLD, num_p, mpi_err)
    !CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)

    num_p = 1
    rank = 0
    !! ADD REMOVED - OCT 2019
    optimal_tasks = CEILING((nangle * 8.0) / (num_p))
    sc_flux_buffer =zero
    reflected_buffer= zero
    ! Allocate angular flux

    ALLOCATE(an_flux%vol_flux(num_moments_v,num_cells)     ,stat=alloc_stat)
    IF(alloc_stat /= 0 .AND. ITMM.NE.0) CALL stop_thor(2_li)
    ALLOCATE(an_flux%face_flux(num_moments_f,0:3,num_cells),stat=alloc_stat)
    IF(alloc_stat /= 0 .AND. ITMM.NE.0) CALL stop_thor(2_li)

    !If we are constructing ITMM matrices, initialize indout
    IF ((ITMM.EQ.2).OR.(ITMM.EQ.3)) THEN
    	indout=0
    END IF

    ! Begin parallel loop over quadrature octant
    DO parallel_i = 1, optimal_tasks
      k = parallel_map_l2g(parallel_i, rank+1)
      IF (k .EQ. 0) EXIT
      oct = MOD(k, 8) + 1
      q = CEILING(k / 8.0)
      indx = indexOf(oct, octants_to_sweep)
      IF (indx > 8_li .or. indx < 1) CALL stop_thor(-1_li, "Sweep order index out of bounds")
      octant = ordering(ordered_octants_to_sweep(indx))
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

      !If we're constructing K ITMM matrices, we have to rezero sc_flux with each direction
      IF (((ITMM.EQ.2) .OR. (ITMM .EQ. 3)).AND.(inner.GT.num_cells)) sc_flux=zero

      ! Set boundary faces and mark as known

      ! 1. Vacuum boundary conditions
      DO i=1,vside_cells
        k=vb_cells(i)%cell
        f=vb_cells(i)%face
        IF ( (omega .dot. outward_normal(k,f)) < 0.0_d_t )THEN
          face_known(f,k) = 1
          ! no nneed to assign 0 to face flux, already done before !!
          ! Unless this is ITMM construction
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
          !tpe=refl_face_tpe(i)
          !IF      (tpe .EQ. 1_li .OR. tpe .EQ. -1_li) THEN
          !  mate = mu_mate(octant)
          !ELSE IF (tpe .EQ. 2_li .OR. tpe .EQ. -2_li) THEN
          !  mate = eta_mate(octant)
          !ELSE IF (tpe .EQ. 3_li .OR. tpe .EQ. -3_li) THEN
          !  mate = xi_mate(octant)
          !END IF
          !an_flux%face_flux(:,f,k) = reflected_flux(:,i,mate,q)
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
        IF ( (omega .dot. outward_normal(k,f)) < 0.0_d_t ) THEN
          face_known(f,k) = 1
          !an_flux%face_flux(:,f,k) = binflx(:,i,octant,q,giflx)
        END IF
      END DO

      ! 4. Faces that are assumed known for breaking cycles

      is_cycle=0
      DO i=1, neldep(octant,q)
        k=eldep(octant,q,eg)%cells(i)
        f=eldep(octant,q,eg)%faces(i)
        is_cycle(f,k)   = 1
        face_known(f,k) = 1
        !an_flux%face_flux(:,f,k) = eldep(octant,q,eg)%face_fluxes(:,i)
      END DO

      ! 5. Faces that are known on an SDD boundary
      DO i=1,SDD_side_cells
        k=SDDb_cells(i)%cell
        f=SDDb_cells(i)%face
        IF ( (omega .dot. outward_normal(k,f)) < 0.0_d_t )THEN
          face_known(f,k) = 1
          ! no nneed to assign 0 to face flux, already done before !!
          ! Unless this is ITMM construction
        END IF
      END DO

      !If this is IPBJ solution, use psi_in to set angular fluxes on incoming boundaries of the SD
      IF ((ITMM.EQ.0).AND.(inner.GT.1)) CALL IPBJ_unpack(eg,octant,q,an_flux)

      !If we're constructing K ITMM matrices, check if the current angular flux is incoming on the current boundary face
      IF ((ITMM.EQ.2).AND.(inner.GT.num_cells)) THEN
      	k=vb_cells(inner-num_cells)%cell
        f=vb_cells(inner-num_cells)%face
        !Incoming
        IF ( (omega .dot. outward_normal(k,f)) < 0.0_d_t ) THEN
        	indin=indin+1
        	an_flux%face_flux(:,f,k) = 1.0_d_t
        	calcK=1
       		ITMMKindex(indin,1,2)=k
       		ITMMKindex(indin,2,2)=f
       		ITMMKindex(indin,3,2)=octant
       		ITMMKindex(indin,4,2)=q
        !Outgoing or parallel
        ELSE
            !IF( (omega .dot. outward_normal(k,f)) .EQ. 0.0_d_t ) &
            !& WRITE(*,*)'%%%WARNING: Parallel angle on sub-domain boundary, incoming K prescription.'
        	!Skip this angle
        	calcK=0
        	GOTO 86
        END IF
      END IF

      !If we are performing IPBJ setup, we still need to construct ITMMKindex
      IF ((ITMM.EQ.3).AND.(inner.GT.num_cells)) THEN
        k=vb_cells(inner-num_cells)%cell
        f=vb_cells(inner-num_cells)%face
        !Incoming
        IF ( (omega .dot. outward_normal(k,f)) < 0.0_d_t ) THEN
        	indin=indin+1
       		ITMMKindex(indin,1,2)=k
       		ITMMKindex(indin,2,2)=f
       		ITMMKindex(indin,3,2)=octant
       		ITMMKindex(indin,4,2)=q
        ELSE
            !IF( (omega .dot. outward_normal(k,f)) .EQ. 0.0_d_t ) &
            !& WRITE(*,*)'%%%WARNING: Parallel angle on sub-domain boundary, pre-process instructions.'
   		END IF
   		calcK=0
   		GOTO 86
      ELSE IF (ITMM.EQ.3) THEN
        calcK=0
        GOTO 87
      END IF



      ! Compute directed source - so far only isotropic

      IF (ITMM.NE.3) THEN
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
      END IF

      ! If page_sweep == 1 read sweep_path from file

      IF(page_sweep.EQ.1) THEN
        sweep_path=0
        READ(99,rec=8*(q-1)+octant) sweep_path(:,1,1)
      END IF

      ! Call queue for mesh sweep

      IF (sweep_tpe .EQ.1) THEN
        CALL queue(eg,q,octant,dir_src,an_flux,LL,U,Lf,Uf,omega,face_known)
      ELSE
        CALL stop_thor(9_li)
      END IF

      ! Increment scalar flux using

      CALL angular_moments(octant,q,sc_flux,an_flux)

      !IF we are constructing ITMM matrices with unit scalar flux, store Jpsi values
      IF ((ITMM.EQ.2).AND.(inner.LE.num_cells)) THEN
      	DO i=1,vside_cells
        	k=vb_cells(i)%cell
        	f=vb_cells(i)%face
        	IF ( (omega .dot. outward_normal(k,f)) > 0.0_d_t ) THEN
          		indout=indout+1
          		Jpsi(indout,inner,eg)=an_flux%face_flux(1,f,k)
        	ELSEIF ( (omega .dot. outward_normal(k,f)) .EQ. 0.0_d_t ) THEN
        		!CALL stop_thor(35_li)
        		!WRITE(*,*)'%%%WARNING: Parallel angle on sub-domain boundary, Jpsi storage.'
        	END IF
      	END DO
      END IF

      !If we are constructing ITMM matrices with unit incoming flux, store Kphi and Kpsi values
      IF ((ITMM.EQ.2).AND.(inner.GT.num_cells)) THEN
      	!Kphi
      	DO ind=1,num_cells
      		Kphi(ind,indin,eg)=sc_flux(1,1,ind)
      	END DO
      	!Kpsi
86      DO i=1,vside_cells
        	k=vb_cells(i)%cell
        	f=vb_cells(i)%face
        	IF ( (omega .dot. outward_normal(k,f)) > 0.0_d_t)  THEN
          		indout=indout+1
          		IF (inner.EQ.num_cells+1) THEN
          		    ITMMKindex(indout,1,1)=k
        		    ITMMKindex(indout,2,1)=f
        		    ITMMKindex(indout,3,1)=octant
        		    ITMMKindex(indout,4,1)=q
        	    END IF
          		IF (calcK.EQ.1 .AND. (an_flux%face_flux(1,f,k).NE.0.0d0)) THEN
          		    !Kpsi(indout,indin)=an_flux%face_flux(1,f,k)
          		    ind_k=ind_k+1
          		    
          		    !If this is the first energy group, check to see if Kpsi needs reallocation
          		    IF (eg.EQ.1) THEN
          		        !Initial allocation
          		        IF (ind_k.EQ.1) THEN
          		            ALLOCATE(KpsiIndexes(2,Kpsi_reallocate),KpsiElements(Kpsi_reallocate,1))

      		            !Reallocation
      		            ELSEIF (MOD(ind_k,Kpsi_reallocate).EQ.1) THEN
      		                reallocate_start=MPI_WTIME()
      		                !Allocate temp vectors
      		                ALLOCATE(KpsiIndexes_temp(2,ind_k-1),KpsiElements_temp(ind_k-1,1))
      		                !Store currents values in temp vectors
      		                DO Kpsi_counter=1,ind_k-1
  		                        KpsiIndexes_temp(:,Kpsi_counter)=KpsiIndexes(:,Kpsi_counter)
  		                        KpsiElements_temp(Kpsi_counter,:)=KpsiElements(Kpsi_counter,:)
      		                END DO
      		                !Reallocate permanent vectors
      		                DEALLOCATE(KpsiIndexes,KpsiElements)
      		                ALLOCATE(KpsiIndexes(2,ind_k-1+Kpsi_reallocate),KpsiElements(ind_k-1+Kpsi_reallocate,1))
      		                KpsiIndexes=0
      		                KpsiElements=0.0d0

      		                !Return values to permanent vectors
      		                DO Kpsi_counter=1,ind_k-1
  		                        KpsiIndexes(:,Kpsi_counter)=KpsiIndexes_temp(:,Kpsi_counter)
  		                        KpsiElements(Kpsi_counter,:)=KpsiElements_temp(Kpsi_counter,:)
      		                END DO
      		                DEALLOCATE(KpsiIndexes_temp,KpsiElements_temp)
      		                Kpsi_reallocate=Kpsi_reallocate*2
      		                reallocate_end=MPI_WTIME()
      		                IF (PBJrank.EQ.0) WRITE(*,'(A,F16.8,A)')'NOTICE: Kpsi sparse matrix storage was reallocated: ', &
      		                & reallocate_end-reallocate_start,' seconds used'
  		                END IF
          		    END IF
          		    
          		    !Store the Kpsi values
          		    KpsiIndexes(1,ind_k)=indout
          		    KpsiIndexes(2,ind_k)=indin
          		    KpsiElements(ind_k,eg)=an_flux%face_flux(1,f,k)
      		    END IF
        	ELSEIF ( (omega .dot. outward_normal(k,f)) .EQ. 0.0_d_t ) THEN
        		!CALL stop_thor(35_li)
        		!WRITE(*,*)'%%%WARNING: Parallel angle on sub-domain boundary, K matrix storage.'
        	END IF
      	END DO
      END IF



      ! update reflected flux array, loop through r
      DO i=1,rside_cells
        k=rb_cells(i)%cell
        f=rb_cells(i)%face
        IF ( (omega .dot. outward_normal(k,f)) > 0.0_d_t )THEN
          ! store outflow on reflective faces
          reflected_buffer(:,i,octant,q) = an_flux%face_flux(:,f,k)

          ! this is necessary to have immediate updates of potentially
          ! relevant reflected fluxes
          reflected_flux(:,i,octant,q) = an_flux%face_flux(:,f,k)
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

      !If this is IPBJ, pack the portion of psi_out coresponding to this angle/octant combo
      IF (ITMM.EQ.0) CALL IPBJ_pack(eg,octant,q,an_flux)

    END DO
    !if ((pbjrank.eq.0).and.(itmm.eq.0).and.(inner.eq.2)) stop

87  CONTINUE

    !! ADD REMOVED - OCT 2019
    reflected_flux=reflected_buffer

    ! Deallocation
    !CALL MPI_AllREDUCE(reflected_flux, reflected_buffer, num_moments_f*rs*8*nangle, &
    !      MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpi_err)
    reflected_flux=reflected_buffer
    !CALL MPI_AllREDUCE(sc_flux, sc_flux_buffer, num_moments_v*namom*num_cells, &
    !      MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpi_err)
    !sc_flux = sc_flux_buffer
    DEALLOCATE(an_flux%vol_flux,an_flux%face_flux)
    !! ADD REMOVED - OCT 2019

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

    INTEGER(kind=li) :: alloc_stat, i,j, oct, f, l, CASE
    REAL(kind=d_t) :: sigmat
    LOGICAL :: all_incoming_known
    TYPE(vector) :: n0, n1, n2, n3
    INTEGER(kind=li) :: soct,sq

    INTEGER ::rank,mpi_err, num_p


    ! Set soct and sq depending on whether page_sweep=0,1

    IF      (page_sweep .EQ. 0) THEN
      soct=octant;sq=q
    ELSE IF (page_sweep .EQ. 1) THEN
      soct=1_li;sq=1_li
    END IF

    ! Initiate loop over all cells

    DO j=1,num_cells

      i=sweep_path(j,soct,sq)

      sigmat=dens_fact(cells(i)%reg)*sigma_t(reg2mat(cells(i)%reg),eg)%xs

      n0=outward_normal(i,0)
      n1=outward_normal(i,1)
      n2=outward_normal(i,2)
      n3=outward_normal(i,3)

      CALL cell_orientation(omega,n0,n1,n2,n3,CASE)

      CALL cell_splitting(sigmat,qm(:,i),an_flux%vol_flux,an_flux%face_flux,   &
            octant,LL,U,Lf,Uf,omega,i,face_known,CASE,n0,n1,n2,n3)

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
    ELSEIF((omega .dot. n0) == 0.0)THEN
    ELSE
      outgoing=outgoing+1
    END IF

    IF((omega .dot. n1) < 0.0)THEN
      incoming=incoming+1
    ELSEIF((omega .dot. n1) == 0.0)THEN
    ELSE
      outgoing=outgoing+1
    END IF

    IF((omega .dot. n2) < 0.0)THEN
      incoming=incoming+1
    ELSEIF((omega .dot. n2) == 0.0)THEN
    ELSE
      outgoing=outgoing+1
    END IF

    IF((omega .dot. n3) < 0.0)THEN
      incoming=incoming+1
    ELSEIF((omega .dot. n3) == 0.0)THEN
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
      CALL stop_thor(10_li)
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
      CALL stop_thor(11_li)
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

    INTEGER(kind=li) :: i, l, order, k, indx, m

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

!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------


  SUBROUTINE IPBJ_pack(eg,octant,q,an_flux)

    !This subroutine creates the psi_out vector from the psi values when using IPBJ

    INTEGER,INTENT(IN)::eg,octant,q
    INTEGER::i,k,f,octant_temp,q_temp
    TYPE(angular_flux),INTENT(IN) :: an_flux

    DO i=1,4*nangle*N_side_SDbound

        !Get the location in phase space for this vector entry
        k=ITMMKindex(i,1,1)
        f=ITMMKindex(i,2,1)
        octant_temp=ITMMKindex(i,3,1)
        q_temp=ITMMKindex(i,4,1)

        !Only pack this value if it is the correct q and octant
        IF ((q.EQ.q_temp).AND.(octant.EQ.octant_temp)) THEN
          psiout(i,eg)=an_flux%face_flux(1,f,k)
        END IF

    END DO



  END SUBROUTINE IPBJ_pack


!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------


  SUBROUTINE IPBJ_unpack(eg,octant,q,an_flux)

    !This subroutine creates the proper an_flux values using the psi_in vector

    INTEGER,INTENT(IN)::eg,octant,q
    INTEGER::i,k,f,octant_temp,q_temp,ierr
    TYPE(angular_flux),INTENT(INOUT) :: an_flux

    DO i=1,4*nangle*N_side_SDbound

        !Get the location in phase space for this vector entry
        k=ITMMKindex(i,1,2)
        f=ITMMKindex(i,2,2)
        octant_temp=ITMMKindex(i,3,2)
        q_temp=ITMMKindex(i,4,2)

        !Only unpack this value if it is the correct q and octant
        IF ((q.EQ.q_temp).AND.(octant.EQ.octant_temp)) THEN
          an_flux%face_flux(1,f,k)=psiin(i,eg)
        END IF

    END DO
        !if (pbjrank.eq.0) stop
!if (inner.eq.2)write(*,*)pbjrank,'pre',sum(an_flux%face_flux)
!stop
  END SUBROUTINE IPBJ_unpack


END MODULE sweep_module
