MODULE inner_iteration_module
  !***********************************************************************
  !
  ! Inner iteration module calls sweep subroutine and updates distributed
  ! source
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

  ! Use modules that pertain setting up problem

  USE sweep_module
  USE termination_module
  USE ITMM_module
  USE SDD_global_variables

  IMPLICIT NONE

CONTAINS

  !> This subroutine performs and inner iteration (compare Alg. 1 in primer)
  SUBROUTINE inner_iteration(eg,sc_flux,q_external,LL,U,Lf,Uf,rs,reflected_flux,prnt)
    !*********************************************************************
    !
    ! Subroutine inner iteration calls sweep subroutine to compute all
    ! fluxes average cell scalar
    !
    !*********************************************************************

    ! Inner Iterations executed for group eg

    INTEGER(kind=li), INTENT(in) :: eg

    ! Declare angular and scalar flux types used globally

    REAL(kind=d_t) :: sc_flux(num_moments_v,namom,num_cells)

    ! Declare source types

    REAL(kind=d_t) :: q_external(num_moments_v,namom,num_cells)

    ! Declare pre-computed matrices

    REAL(kind=d_t), DIMENSION(num_moments_v,num_moments_v) :: LL, U
    REAL(kind=d_t), DIMENSION(num_moments_f,num_moments_f) :: Lf, Uf

    ! Define reflected flux array

    INTEGER(kind=li) :: rs
    REAL(kind=d_t),DIMENSION(num_moments_f,rs,8,nangle) :: reflected_flux

    ! print convergence monitor flag

    LOGICAL :: prnt

    ! Define temporary variables

    INTEGER(kind=li) :: l, i, q, octant, order, alloc_stat,&
          n,m, k, indx, ind, j, Kpsi_counter
    REAL(kind=d_t) :: error,te,ts

    ! Define self-scatter source

    REAL(kind=d_t) :: self_scatter(num_moments_v,namom,num_cells)

    ! Define old scalar flux

    REAL(kind=d_t) :: sc_flux_old (num_moments_v,namom,num_cells)

    ! Define MPI environment and get rank
    INTEGER ::rank,mpi_err, localunit

    REAL(kind=d_t) :: per_inner_start_time, per_inner_end_time, per_comm_start_time
    !! ADD REMOVED - OCT 2019
    !CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)
    rank = 0
    !! ADD REMOVED - OCT 2019


    ! Set self_scatter and old scalar flux to 0

    self_scatter = zero
    sc_flux_old  = zero

    ! write header for convergence monitor
    IF (PBJrank .EQ. 0 .AND. ITMM .NE. 2 .AND. ITMM.NE.3) THEN
      IF(prnt) WRITE(6,102) '  grp  itn       error        time'
    END IF
    IF (PBJrank .EQ. 0 .AND. ITMM .EQ. 2) THEN
      IF(prnt) WRITE(6,'(A,I5,A,I5)') '  Beginning GFIC for group: ', eg, " of: ", egmax
    END IF
    IF (PBJrank .EQ. 0 .AND. ITMM .EQ. 3) THEN
      IF(prnt) WRITE(6,'(A,I5,A,I5)') '  Beginning IPBJ pre-process for group: ', eg, " of: ", egmax
    END IF
    ! Begin inner iteration
    IF ((ITMM.EQ.1).OR.(ITMM.EQ.0)) THEN
  	  psiin=0.0d0
    END IF

    DO inner=1, max_inner

        per_inner_start_time=MPI_WTIME()

        IF (PBJrank .EQ. 0 .AND. ITMM .EQ. 2) THEN
            IF (inner .EQ. 1) WRITE(6,*)  "Begin J matrices"
            IF (inner .EQ. num_cells + 1) WRITE(6,*) "Begin K matrices"
        END IF

      ! start timer

      CALL CPU_TIME(ts)

      !If we're constructing ITMM matrices, zero out sc_flux and set one value to 1
      IF ((ITMM.EQ.2).OR.(ITMM.EQ.3)) THEN
      	sc_flux=zero
      	!J construction
      	IF (inner.LE.num_cells) THEN
      		sc_flux(:,:,inner)=one
      	!K construction
      	ELSE
      		IF (inner.EQ.num_cells+1) THEN
      		    indin=0
      		    ind_k=0
  		    END IF
      	END IF
      END IF

      ! Recompute self-scattering source and add external source
      IF (ITMM.NE.3) THEN
          DO i=1,num_cells
            ! even contribution
            DO l=0,scatt_ord
              DO m=0,l
                indx=1_li+m+(l+1_li)*l/2_li
                DO k=1, num_moments_v
                  self_scatter(k,indx,i) = scat_mult(l,m)                             *&
                        sigma_scat(reg2mat(cells(i)%reg),l+1,eg,eg)%xs     *&
                        dens_fact(cells(i)%reg)                            *&
                        sc_flux(k,indx,i) + q_external(k,indx,i)
                END DO
              END DO
            END DO
            ! odd contribution
            DO l=1,scatt_ord
            DO m=1,l
                indx=neven+m+(l-1_li)*l/2_li
                DO k=1, num_moments_v
                self_scatter(k,indx,i) = scat_mult(l,m) *&
                        sigma_scat(reg2mat(cells(i)%reg),l+1,eg,eg)%xs     *&
                        dens_fact(cells(i)%reg)                            *&
                        sc_flux(k,indx,i) + q_external(k,indx,i)
                END DO
            END DO
            END DO

        END DO
      END IF

      ! Call sweep algorithm or ITMM solve
      IF (ITMM.NE.1) THEN
      	sc_flux=zero
      	CALL sweep(eg,sc_flux,self_scatter,rs,reflected_flux,LL,U,Lf,Uf)
      ELSE
      	CALL ITMMsolve(eg,sc_flux,q_external)
      END IF

      ! Computer error

      !If we are performing GFIC and have a unit scalar flux, store column of Jphi
      IF ((ITMM.EQ.2).AND.(inner.LE.num_cells)) THEN
      	DO ind=1,num_cells
      		Jphi(ind,inner,eg)=sc_flux(1,1,ind)
      	END DO
      END IF
      
      
      !If this is the first energy group, determine final size of Kpsi vectors
      IF ((ITMM.EQ.2) .AND. (eg.EQ.1) .AND. (inner.EQ.max_inner)) THEN
      
        !Allocate temp vectors
        ALLOCATE(KpsiIndexes_temp(2,ind_k),KpsiElements_temp(ind_k,1))
        IF (PBJrank.EQ.0) WRITE(*,'(2(A,I0),A,F5.1,A)')'INFO: Kpsi is ',ind_k,' elements. Maximum allocated in construction: ', &
                                                & Kpsi_reallocate,' (', (ind_k/(1.0*kpsi_reallocate))*100,'%)'
        !Store currents values in temp vectors
        DO Kpsi_counter=1,ind_k
            KpsiIndexes_temp(:,Kpsi_counter)=KpsiIndexes(:,Kpsi_counter)
            KpsiElements_temp(Kpsi_counter,:)=KpsiElements(Kpsi_counter,:)
        END DO
        !Reallocate permanent vectors
        DEALLOCATE(KpsiIndexes,KpsiElements)
        ALLOCATE(KpsiIndexes(2,ind_k),KpsiElements(ind_k,egmax))
        KpsiElements=0.0d0
        !Return values to permanent vectors
        DO Kpsi_counter=1,ind_k
            KpsiIndexes(:,Kpsi_counter)=KpsiIndexes_temp(:,Kpsi_counter)
            KpsiElements(Kpsi_counter,1)=KpsiElements_temp(Kpsi_counter,1)
        END DO
        DEALLOCATE(KpsiIndexes_temp,KpsiElements_temp)
        nonzero=ind_k
      
      END IF
    
      !Transition Kpsi to sparse storage
      !IF ((ITMM.EQ.2).AND.(inner.EQ.max_inner)) THEN

      	!If this is the first group, allocate
      	!IF (eg.EQ.1) THEN

      		!Determine number of nonzero values
      		!nonzero=0
      		!DO i=1,4*nangle*N_side_SDbound
      			!DO j=1,4*nangle*N_side_SDbound
      				!IF (ABS(Kpsi(i,j)).GT.0.0d0) nonzero=nonzero+1
      			!END DO
      		!END DO

      		!Allocate
      		!ALLOCATE(KpsiElements(nonzero,egmax),KpsiIndexes(2,nonzero))

      		!Get indexes
      		!nonzero=0
      		!DO i=1,4*nangle*N_side_SDbound
      			!DO j=1,4*nangle*N_side_SDbound
      				!IF (ABS(Kpsi(i,j)).GT.0.0d0) THEN
      					!nonzero=nonzero+1
      					!KpsiIndexes(1,nonzero)=i
      					!KpsiIndexes(2,nonzero)=j
      				!END IF
      			!END DO
      		!END DO

      	!END IF

      	!Pull nonzero elements from Kspi to KpsiElements
      	!DO i=1,nonzero
      		!KpsiElements(i,eg)=Kpsi(KpsiIndexes(1,i),KpsiIndexes(2,i))
      	!END DO

      	!Reset Kpsi to zero
      	!Kpsi=zero

      !END IF
      
      

      max_error(eg)=zero

      IF ((ITMM.EQ.0).OR.(ITMM.EQ.1)) THEN

          DO i=1, num_cells
            IF(ABS(sc_flux(1,1,i)) > 1e-12)THEN
              error=ABS( (  sc_flux(1,1,i) - sc_flux_old(1,1,i) )/&
                    sc_flux(1,1,i) )
            ELSE
              error=ABS(    sc_flux(1,1,i) - sc_flux_old(1,1,i) )
            END IF
            IF(error >= max_error(eg))THEN
              max_error(eg)=error
            END IF
          END DO

      END IF

      !Get the global maximum error
      IF ((ITMM .EQ. 0) .OR. (ITMM .EQ. 1)) THEN
        comm_start_time=MPI_WTIME()
        CALL MPI_ALLREDUCE(max_error(eg),max_error_g,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,mpi_err)
        comm_time=comm_time+MPI_WTIME()-comm_start_time
      END IF

      ! Test convergence of current-group scalar flux

      IF(inner > 1 .AND. max_error_g < inner_conv .AND. ((ITMM .EQ. 1) .OR. (ITMM .EQ. 0)))THEN
        go to 10
      END IF

      IF ((ITMM.EQ.0).OR.(ITMM.EQ.1)) THEN

          IF(inner < max_inner)THEN
            DO i=1, num_cells
              DO n=1,namom
                DO l=1, num_moments_v
                  sc_flux_old(l,n,i)=sc_flux(l,n,i)
                END DO
              END DO
            END DO
          END IF

      END IF

      ! stop timer

      CALL CPU_TIME(te)

      ! write convergence monitor
      IF (PBJrank .EQ. 0 .AND. ((ITMM .EQ. 0) .OR. (ITMM .EQ. 1))) THEN
        IF(prnt) THEN
          WRITE(6,101) eg,inner,max_error(eg),te-ts,' % '
          flush(6)
        END IF
        IF(prnt .AND. print_conv.EQ.1) THEN
          WRITE(21,101) eg,inner,max_error(eg),te-ts,' % '
          flush(21)
        END IF
101     FORMAT (1X,2I5,2ES12.4,A)
102     FORMAT (1X,A)
      END IF

      per_comm_start_time=MPI_WTIME()

      !ITMM (or IPBJ) Communication
      IF ((ITMM.EQ.1).OR.(ITMM.EQ.0)) THEN
        
        comm_start_time=MPI_WTIME()
        !Pack the messages
        CALL ITMM_pack(eg)

        !Communicate with adjacent processors
        CALL ITMM_comm(eg)

        !Unpack the messages
        CALL ITMM_unpack(eg)
        comm_time=comm_time+MPI_WTIME()-comm_start_time

      END IF

      per_inner_end_time=MPI_WTIME()

      IF (tot_nInners+inner.LE.num_dbg_records) per_comm_time(tot_nInners+inner) = &
        & per_inner_end_time-per_comm_start_time
      IF (tot_nInners+inner.LE.num_dbg_records) per_iter_time(tot_nInners+inner)= &
        & per_inner_end_time-per_inner_start_time

    END DO

    inner=inner-1_li

10  CONTINUE

    per_inner_end_time=MPI_WTIME()

    IF (tot_nInners+inner.LE.num_dbg_records) per_comm_time(tot_nInners+inner) = &
      & per_inner_end_time-per_comm_start_time
    IF (tot_nInners+inner.LE.num_dbg_records) per_iter_time(tot_nInners+inner)= &
      & per_inner_end_time-per_inner_start_time

    tot_nInners=tot_nInners+inner




  END SUBROUTINE inner_iteration

END MODULE inner_iteration_module
