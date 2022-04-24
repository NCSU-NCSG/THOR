MODULE solver_module
  !***********************************************************************
  !
  ! Solver module initializes fluxes and sources, then determines type of
  ! problem and calls outer iteration
  !
  !***********************************************************************

  ! User derived-type modules

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

  USE ahotc_matrix_module
  USE sph_harmonics_module
  USE outer_iteration_module
  USE jfnk_module
  USE error_module

  IMPLICIT NONE

CONTAINS

  !> Prepares some data and allocates memory, then hands the solve over to
  !> the outer_iteration_ext subroutine
  SUBROUTINE solver_ext(flux)
    !*********************************************************************
    !
    ! Subroutine solver calls outer iteration subroutine
    !
    !*********************************************************************

    ! Pass scalar flux

    REAL(kind=d_t) :: flux(num_moments_v,namom,num_cells,egmax,niter)

    ! Define temporary variables

    INTEGER(kind=li) :: alloc_stat
    REAL(kind=d_t), DIMENSION(:,:), ALLOCATABLE :: M, LL, U, &
          Mf, Lf, Uf

    ! Allocate and initialize all necessary derived types

    ALLOCATE(M(num_moments_v,num_moments_v),&
          LL(num_moments_v,num_moments_v),&
          U(num_moments_v,num_moments_v),&
          Mf(num_moments_f,num_moments_f),&
          Lf(num_moments_f,num_moments_f),&
          Uf(num_moments_f,num_moments_f),&
          stat=alloc_stat);M=zero;&
          LL=zero;U=zero;Mf=zero;Lf=zero;Uf=zero;
    IF(alloc_stat /= 0) CALL raise_fatal_error("*** Not enough memory ***")


    ALLOCATE(Ysh(nangle,8,namom),stat=alloc_stat)
    IF(alloc_stat /=0) CALL raise_fatal_error("*** Not enough memory ***")

    ! Pre-compute and apply LU decomposition to 'mass matrices'

    CALL ahotc_matrix(num_moments_v,num_moments_f,index_v,&
          index_f,M,LL,U,Mf,Lf,Uf)

    ! Pre-compute the spherical harmonics basis functions

    CALL spherical_harmonics(scatt_ord,nangle,quadrature,namom,Ysh)

    ! Initate outer iteration
    IF (rank .EQ. 0) THEN
      WRITE(6,*) '-- Commencing fixed source computation.'
    END IF
    CALL outer_iteration_ext(flux,LL,U,Lf,Uf)

    ! Deallocate temporary arrays

    DEALLOCATE(M,Mf,LL,U,Lf,Uf)

  END SUBROUTINE solver_ext

  !-----------------------------------------------------------------------------------------
  !-----------------------------------------------------------------------------------------

  !> Prepares some data and allocates memory, then hands the solve over to
  !> the outer_iteration_eig subroutine or the JFNK solver. JFNK is a different
  !> solution method developed by Dan Gill.
  SUBROUTINE solver_eig(flux,keff)
    !*********************************************************************
    !
    ! Subroutine solver calls outer iteration subroutine
    !
    !*********************************************************************



    ! Pass eigenvalue and scalar flux

    REAL(kind=d_t), INTENT(inout) :: keff
    REAL(kind=d_t) :: flux(num_moments_v,namom,num_cells,egmax,niter)

    ! Define temporary variables

    INTEGER(kind=li) :: alloc_stat
    REAL(kind=d_t), DIMENSION(:,:), ALLOCATABLE :: M, LL, U, &
          Mf, Lf, Uf

    ! Allocate and initialize all necessary derived types

    ALLOCATE(M(num_moments_v,num_moments_v),&
          LL(num_moments_v,num_moments_v),&
          U(num_moments_v,num_moments_v),&
          Mf(num_moments_f,num_moments_f),&
          Lf(num_moments_f,num_moments_f),&
          Uf(num_moments_f,num_moments_f),&
          stat=alloc_stat);M=zero;&
          LL=zero;U=zero;Mf=zero;Lf=zero;Uf=zero;
    IF(alloc_stat /= 0) CALL raise_fatal_error("*** Not enough memory ***")

    ALLOCATE(Ysh(nangle,8,namom),stat=alloc_stat)
    IF(alloc_stat /=0) CALL raise_fatal_error("*** Not enough memory ***")

    ! Pre-compute and apply LU decomposition to 'mass matrices'

    CALL ahotc_matrix(num_moments_v,num_moments_f,index_v,&
          index_f,M,LL,U,Mf,Lf,Uf)

    ! Pre-compute the spherical harmonics basis functions

    CALL spherical_harmonics(scatt_ord,nangle,quadrature,namom,Ysh)

    ! Initate outer iteration
    IF (eig_switch .EQ. 0) THEN
      IF (rank .EQ. 0) THEN
        WRITE(6,*) '-- Commencing eigenvalue computation.'
        WRITE(6,*) '-- A power iteration computation is executed.'
      END IF
      CALL outer_iteration_eig(flux,keff,LL,U,Lf,Uf)
    ELSE
      IF (rank .EQ. 0) THEN
        WRITE(6,*) '-- Commencing eigenvalue computation.'
        WRITE(6,*) '-- A JFNK computation is executed.'
      END IF
      CALL set_jfnk

      CALL do_jfnk(flux,keff,LL,U,Lf,Uf)

      CALL clean_jfnk

    END IF

    ! Deallocate temporary arrays

    DEALLOCATE(M,Mf,LL,U,Lf,Uf)

  END SUBROUTINE solver_eig

  !-----------------------------------------------------------------------------------------
  ! End
  !-----------------------------------------------------------------------------------------

END MODULE solver_module
