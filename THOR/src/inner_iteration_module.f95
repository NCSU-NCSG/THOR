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
  USE globals

  ! Use modules that pertain setting up problem

  USE sweep_module
  USE error_module

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

    INTEGER(kind=li) :: l, i,&
          n,m, k, indx, mat_indx
    REAL(kind=d_t) :: error,te,ts

    ! Define self-scatter source

    REAL(kind=d_t) :: self_scatter(num_moments_v,namom,num_cells)

    ! Define old scalar flux

    REAL(kind=d_t) :: sc_flux_old (num_moments_v,namom,num_cells)

    ! Define MPI environment and get rank
    INTEGER ::rank,mpi_err
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)



    ! Set self_scatter and old scalar flux to 0

    self_scatter = zero
    sc_flux_old  = zero

    ! write header for convergence monitor
    IF (rank .EQ. 0) THEN
      IF(prnt)CALL printlog('   grp  itn       error        time')
    END IF
    ! Begin inner iteration

    DO inner=1, max_inner

      ! start timer

      CALL CPU_TIME(ts)

      ! Recompute self-scattering source and add external source

      DO i=1,num_cells
        mat_indx=material_ids(reg2mat(cells(i)%reg))
        ! even contribution
        DO l=0,scatt_ord
          DO m=0,l
            indx=1_li+m+(l+1_li)*l/2_li
            DO k=1, num_moments_v
              self_scatter(k,indx,i) = scat_mult(l,m)                             *&
                    xs_mat(mat_indx)%sigma_scat(l+1,eg,eg)*dens_fact(cells(i)%reg)*&
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
                    xs_mat(mat_indx)%sigma_scat(l+1,eg,eg)*dens_fact(cells(i)%reg)*&
                    sc_flux(k,indx,i) + q_external(k,indx,i)
            END DO
          END DO
        END DO

      END DO

      ! Call sweep algorithm

      sc_flux=zero
      CALL sweep(eg,sc_flux,self_scatter,rs,reflected_flux,LL,U,Lf,Uf)
      ! Computer error

      max_error(eg)=zero

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

      ! Test convergence of current-group scalar flux

      IF(inner > 1 .AND. max_error(eg) < inner_conv)THEN
        EXIT
      END IF

      IF(inner < max_inner)THEN
        DO i=1, num_cells
          DO n=1,namom
            DO l=1, num_moments_v
              sc_flux_old(l,n,i)=sc_flux(l,n,i)
            END DO
          END DO
        END DO
      END IF

      ! stop timer

      CALL CPU_TIME(te)

      ! write convergence monitor
      IF (rank .EQ. 0) THEN
        IF(prnt) THEN
          WRITE(amsg,101) eg,inner,max_error(eg),te-ts,' % '
          CALL printlog(amsg)
          flush(6)
        END IF
        IF(prnt .AND. print_conv.EQ.1) THEN
          WRITE(21,101) eg,inner,max_error(eg),te-ts,' % '
          flush(21)
        END IF
101     FORMAT (1X,2I5,2ES12.4,A)
      END IF
    END DO

    inner=inner-1_li

    tot_nInners=tot_nInners+inner

  END SUBROUTINE inner_iteration

END MODULE inner_iteration_module
