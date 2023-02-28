MODULE cell_splitting_module
  !***********************************************************************
  !
  ! Cell splitting module divides each tetrahedral cell into subcells
  ! along the discrete ordinate of interest. Six possible configurations
  ! account for all the possible combinations of incoming/outgoing faces.
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
  USE globals
  USE stringmod

  ! Use modules that pertain setting up problem
  ! FIXME: You need to find a better place to store cell_jacobian
  !        it's defined in all transport kernels...and then used in cell_splitting
  !
  USE general_utility_module
  USE transport_kernel_module_sc, ONLY: transport_kernel_sc
  USE transport_kernel_module_lc, ONLY: transport_kernel_lc
  USE transport_kernel_module_cce, ONLY: transport_kernel_cce
  USE error_module

  IMPLICIT NONE

CONTAINS

  !> THOR is based in the AHOT-C methods. The AHOT-C method splits tetrahedra
  !> into characteristic tetrahedra and solves the characteristic equations
  !> on each of these. This is the entry point for the single mesh-cell solver.
  SUBROUTINE cell_splitting(sigmat,qm_moments,vol_flux,face_flux,  &
        LL,U,Lf,Uf,omega,i,face_known,CASE,n0,n1,n2,n3)
    !*********************************************************************
    !
    ! Subroutine transport kernel divides tetrahedra into 'CTs' and
    ! calculates outgoing face angular fluxes based on incoming ones
    !
    !*********************************************************************

    ! Pass arguments

    REAL(kind=d_t), INTENT(in)                             :: sigmat
    REAL(kind=d_t), DIMENSION(num_moments_v), INTENT(in)   :: qm_moments
    REAL(kind=d_t), DIMENSION(num_moments_v,num_cells)     :: vol_flux
    REAL(kind=d_t), DIMENSION(num_moments_f,0:3,num_cells) :: face_flux
    INTEGER(kind=li), INTENT(in)                           :: i, CASE
    REAL(kind=d_t), DIMENSION(num_moments_v,num_moments_v) :: LL,U
    REAL(kind=d_t), DIMENSION(num_moments_f,num_moments_f) :: Lf,Uf
    TYPE(vector), INTENT(in)                               :: omega, n0, n1, n2,n3
    INTEGER(kind=1), DIMENSION(0:3,num_cells)              :: face_known

    ! Define temporary variables

    INTEGER(kind=li)                            :: alloc_stat, subcells
    INTEGER(kind=li), DIMENSION(:), ALLOCATABLE :: incoming_face,outgoing_face
    INTEGER(kind=li), DIMENSION(:), ALLOCATABLE :: bc
    REAL(kind=d_t)                              :: t, det_J
    REAL(kind=d_t), DIMENSION(3,3)              :: J, J_inv
    REAL(kind=d_t), DIMENSION(:,:), ALLOCATABLE :: upstream_moments
    INTEGER(kind=1)                             :: face0, face1, face2, face3
    TYPE(vector)                                :: omega_local, v0, v1, v2, v3
    TYPE(vector), DIMENSION(:), ALLOCATABLE     :: r0, r1, r2, r3

    ! Call Jacobian subroutine to generate cell transformation

    v0=vertices(cells(i)%R(0))%v
    v1=vertices(cells(i)%R(1))%v
    v2=vertices(cells(i)%R(2))%v
    v3=vertices(cells(i)%R(3))%v

    CALL cell_jacobian(v0,v1,v2,v3,J,det_J,J_inv)

    ! Compute local coordinate omega_local

    omega_local=J_inv*omega

    ! Split cell into subcells depending on cases

    face0=face_known(0,i)
    face1=face_known(1,i)
    face2=face_known(2,i)
    face3=face_known(3,i)

    ! CASE 1: 3 incoming, 1 outgoing
    IF(CASE == 1)THEN
      subcells=3

      ALLOCATE(r0(subcells),r1(subcells),r2(subcells),&
            r3(subcells),incoming_face(subcells),&
            outgoing_face(subcells),upstream_moments(num_moments_f,subcells),&
            bc(subcells),stat=alloc_stat)
      IF(alloc_stat /= 0) CALL raise_fatal_error("*** Not enough memory ***")

      CALL case1(v0,v1,v2,v3,omega,r0,r1,r2,r3,face0,face1,face2,&
            face3,incoming_face,outgoing_face,omega_local,t)

      ! CASE 2: 1 incoming, 3 outgoing
    ELSEIF(CASE == 2)THEN
      subcells=3

      ALLOCATE(r0(subcells),r1(subcells),r2(subcells),&
            r3(subcells),incoming_face(subcells),&
            outgoing_face(subcells),upstream_moments(num_moments_f,subcells),&
            bc(subcells),stat=alloc_stat)
      IF(alloc_stat /= 0) CALL raise_fatal_error("*** Not enough memory ***")

      CALL case2(v0,v1,v2,v3,omega,r0,r1,r2,r3,face0,face1,face2,&
            face3,incoming_face,outgoing_face,omega_local,t)

      ! CASE 3: 2 incoming, 2 outgoing
    ELSEIF(CASE == 3)THEN
      subcells=4

      ALLOCATE(r0(subcells),r1(subcells),r2(subcells),&
            r3(subcells),incoming_face(subcells),&
            outgoing_face(subcells),upstream_moments(num_moments_f,subcells),&
            bc(subcells),stat=alloc_stat)
      IF(alloc_stat /= 0) CALL raise_fatal_error("*** Not enough memory ***")

      CALL case3(v0,v1,v2,v3,r0,r1,r2,r3,face0,face1,face2,face3,&
            incoming_face,outgoing_face,J,omega_local,t)

      ! CASE 4: 2 incoming, 1 outgoing
    ELSEIF(CASE == 4)THEN
      subcells=2

      ALLOCATE(r0(subcells),r1(subcells),r2(subcells),&
            r3(subcells),incoming_face(subcells),&
            outgoing_face(subcells),upstream_moments(num_moments_f,subcells),&
            bc(subcells),stat=alloc_stat)
      IF(alloc_stat /= 0) CALL raise_fatal_error("*** Not enough memory ***")

      CALL case4(v0,v1,v2,v3,omega,n0,n1,n2,n3,r0,r1,r2,r3,face0,&
            face1,face2,face3,incoming_face,outgoing_face,omega_local,t)

      ! CASE 5: 1 incoming, 2 outgoing
    ELSEIF(CASE == 5)THEN
      subcells=2

      ALLOCATE(r0(subcells),r1(subcells),r2(subcells),&
            r3(subcells),incoming_face(subcells),&
            outgoing_face(subcells),upstream_moments(num_moments_f,subcells),&
            bc(subcells),stat=alloc_stat)
      IF(alloc_stat /= 0) CALL raise_fatal_error("*** Not enough memory ***")

      CALL case5(v0,v1,v2,v3,omega,n0,n1,n2,n3,r0,r1,r2,r3,face0,&
            face1,face2,face3,incoming_face,outgoing_face,omega_local,t)

      ! CASE 6: 1 incoming, 1 outgoing
    ELSEIF(CASE == 6)THEN
      subcells=1

      ALLOCATE(r0(subcells),r1(subcells),r2(subcells),&
            r3(subcells),incoming_face(subcells),&
            outgoing_face(subcells),upstream_moments(num_moments_f,subcells),&
            bc(subcells),stat=alloc_stat)
      IF(alloc_stat /= 0) CALL raise_fatal_error("*** Not enough memory ***")

      CALL case6(v0,v1,v2,v3,omega,n0,n1,n2,n3,r0,r1,r2,r3,face0,&
            face1,face2,face3,incoming_face,outgoing_face,omega_local,t)

      ! CASE UNKNOWN: failure
    ELSE
      CALL raise_fatal_error("Unacceptable case from cell splitting in cell"//TRIM(STR(i)))
    END IF

    CALL upstream_mom(adjacent_cells,num_moments_f,     &
          adjacency_list,face_flux,i,subcells,incoming_face,&
          upstream_moments)

    IF (space_ord == 0) THEN
      CALL transport_kernel_SC(sigmat,qm_moments,vol_flux(:,i),                &
            face_flux(:,:,i),i,                  &
            t,subcells,r0,r1,r2,r3,face_known,incoming_face,&
            outgoing_face,upstream_moments)
    ELSE IF (space_ord == -1) THEN
      CALL transport_kernel_LC(sigmat,qm_moments,vol_flux(:,i),                &
            face_flux(:,:,i),i,                  &
            t,subcells,r0,r1,r2,r3,face_known,incoming_face,&
            outgoing_face,J,J_inv,upstream_moments)
    ELSE
      CALL transport_kernel_CCE(sigmat,qm_moments,vol_flux(:,i),                &
            face_flux(:,:,i),LL,U,Lf,Uf,i,                  &
            t,subcells,r0,r1,r2,r3,face_known,incoming_face,&
            outgoing_face,J,J_inv,upstream_moments)

    END IF

    DEALLOCATE(r0,r1,r2,r3,upstream_moments)

  END SUBROUTINE cell_splitting

  SUBROUTINE upstream_mom(nadj,numf,adj_list,face_flux,i,&
        subcells,incoming_face,upstream_moments)
    !*********************************************************************
    !
    ! Subroutine upstream_mom determines upstream moments
    !
    !*********************************************************************
    ! Pass geometry parameters

    INTEGER(kind=li), INTENT(in) :: nadj

    ! Declare index size

    INTEGER(kind=li), INTENT(in) :: numf

    ! Pass geometry derived type

    TYPE(list), DIMENSION(nadj,0:3), INTENT(in) :: adj_list

    ! Declare angular and scalar flux types used globally

    REAL(kind=d_t) :: face_flux(num_moments_f,0:3,num_cells)

    ! Define temporary variables

    INTEGER(kind=li)                                        :: ii, up_cell, up_face, l
    INTEGER(kind=li), INTENT(in)                            :: i, subcells
    INTEGER(kind=li), DIMENSION(subcells), INTENT(in)       :: incoming_face
    REAL(kind=d_t), DIMENSION(numf,subcells), INTENT(inout) :: upstream_moments

    ! Assign upstream values and boundary conditions


    DO ii=1, subcells
      up_cell=adj_list(i,incoming_face(ii))%cell
      up_face=adj_list(i,incoming_face(ii))%face
      IF(adj_list(i,incoming_face(ii))%cell .NE. 0 .AND. is_cycle(incoming_face(ii),i).EQ.0 )THEN
        DO l=1, numf
          upstream_moments(l,ii)=face_flux(l,up_face,up_cell)
        END DO
      ELSE
        DO l=1, numf
          upstream_moments(l,ii)=face_flux(l,incoming_face(ii),i)
        END DO
      END IF
    END DO

  END SUBROUTINE upstream_mom

  SUBROUTINE case1(v0,v1,v2,v3,omega,r0,r1,r2,r3,face0,face1,face2,&
        face3,incoming_face,outgoing_face,omega_local,t)
    !*********************************************************************
    !
    ! Subroutine case1 generates local coordinate subcells vectors
    !
    !*********************************************************************
    ! Define temporary variables

    INTEGER(kind=li), DIMENSION(3), INTENT(inout) :: incoming_face,&
          outgoing_face
    REAL(kind=d_t), INTENT(out) :: t
    INTEGER(kind=1), INTENT(in) :: face0,face1,face2,face3
    TYPE(vector) :: Rin, Rout
    TYPE(vector), INTENT(in) :: omega, omega_local, v0, v1, v2, v3
    TYPE(vector), DIMENSION(3),INTENT(inout) :: r0, r1, r2, r3

    ! Split cells and define subcells vertices r0, r1, r2, & r3

    IF(face0 .EQ. 0)THEN
      t=1.0_d_t/(omega_local%x1)
      Rin=v0
      Rout=Rin+t*omega
      ! Subcell 1
      incoming_face(1)=3
      outgoing_face(1)=0
      r0(1)=v1
      r1(1)=v2
      r2(1)=Rin
      r3(1)=Rout
      ! Subcell 2
      incoming_face(2)=1
      outgoing_face(2)=0
      r0(2)=v2
      r1(2)=v3
      r2(2)=Rin
      r3(2)=Rout
      ! Subcell 3
      incoming_face(3)=2
      outgoing_face(3)=0
      r0(3)=v3
      r1(3)=v1
      r2(3)=Rin
      r3(3)=Rout
    ELSEIF(face1 .EQ. 0)THEN
      t=1.0_d_t/(omega_local%x2-omega_local%x1)
      Rin=v1
      Rout=Rin+t*omega
      ! Subcell 1
      incoming_face(1)=2
      outgoing_face(1)=1
      r0(1)=v0
      r1(1)=v3
      r2(1)=Rin
      r3(1)=Rout
      ! Subcell 2
      incoming_face(2)=0
      outgoing_face(2)=1
      r0(2)=v3
      r1(2)=v2
      r2(2)=Rin
      r3(2)=Rout
      ! Subcell 3
      incoming_face(3)=3
      outgoing_face(3)=1
      r0(3)=v2
      r1(3)=v0
      r2(3)=Rin
      r3(3)=Rout
    ELSEIF(face2 .EQ. 0)THEN
      t=1.0_d_t/(omega_local%x3-omega_local%x2)
      Rin=v2
      Rout=Rin+t*omega
      ! Subcell 1
      incoming_face(1)=3
      outgoing_face(1)=2
      r0(1)=v0
      r1(1)=v1
      r2(1)=Rin
      r3(1)=Rout
      ! Subcell 2
      incoming_face(2)=0
      outgoing_face(2)=2
      r0(2)=v1
      r1(2)=v3
      r2(2)=Rin
      r3(2)=Rout
      ! Subcell 3
      incoming_face(3)=1
      outgoing_face(3)=2
      r0(3)=v3
      r1(3)=v0
      r2(3)=Rin
      r3(3)=Rout
    ELSEIF(face3 .EQ. 0)THEN
      t=1.0_d_t/(-1.0_d_t*omega_local%x3)
      Rin=v3
      Rout=Rin+t*omega
      ! Subcell 1
      incoming_face(1)=1
      outgoing_face(1)=3
      r0(1)=v0
      r1(1)=v2
      r2(1)=Rin
      r3(1)=Rout
      ! Subcell 2
      incoming_face(2)=0
      outgoing_face(2)=3
      r0(2)=v2
      r1(2)=v1
      r2(2)=Rin
      r3(2)=Rout
      ! Subcell 3
      incoming_face(3)=2
      outgoing_face(3)=3
      r0(3)=v1
      r1(3)=v0
      r2(3)=Rin
      r3(3)=Rout
    ELSE
      CALL raise_fatal_error("Cell splitting case 1 failed for cell")
    END IF

  END SUBROUTINE case1

  SUBROUTINE case2(v0,v1,v2,v3,omega,r0,r1,r2,r3,face0,face1,face2,&
        face3,incoming_face,outgoing_face,omega_local,t)
    !*********************************************************************
    !
    ! Subroutine case2 generates local coordinate subcell vectors
    !
    !*********************************************************************
    ! Define temporary variables

    INTEGER(kind=li), DIMENSION(3), INTENT(inout) :: incoming_face,&
          outgoing_face
    REAL(kind=d_t), INTENT(out) :: t
    INTEGER(kind=1), INTENT(in) :: face0,face1,face2,face3
    TYPE(vector) :: Rin, Rout
    TYPE(vector), INTENT(in) :: omega, omega_local, v0, v1, v2, v3
    TYPE(vector), DIMENSION(3),INTENT(inout) :: r0, r1, r2, r3

    ! Split cell and define subcell vertices r0, r1, r2, & r3

    IF(face0 .EQ. 1)THEN
      t=1.0_d_t/(-1.0_d_t*omega_local%x1)
      Rout=v0
      Rin=Rout-t*omega
      ! Subcell 1
      incoming_face(1)=0
      outgoing_face(1)=2
      r0(1)=v1
      r1(1)=v3
      r2(1)=Rin
      r3(1)=Rout
      ! Subcell 2
      incoming_face(2)=0
      outgoing_face(2)=1
      r0(2)=v3
      r1(2)=v2
      r2(2)=Rin
      r3(2)=Rout
      ! Subcell 3
      incoming_face(3)=0
      outgoing_face(3)=3
      r0(3)=v2
      r1(3)=v1
      r2(3)=Rin
      r3(3)=Rout
    ELSEIF(face1 .EQ. 1)THEN
      t=1.0_d_t/(omega_local%x1-omega_local%x2)
      Rout=v1
      Rin=Rout-t*omega
      ! Subcell 1
      incoming_face(1)=1
      outgoing_face(1)=3
      r0(1)=v0
      r1(1)=v2
      r2(1)=Rin
      r3(1)=Rout
      ! Subcell 2
      incoming_face(2)=1
      outgoing_face(2)=0
      r0(2)=v2
      r1(2)=v3
      r2(2)=Rin
      r3(2)=Rout
      ! Subcell 3
      incoming_face(3)=1
      outgoing_face(3)=2
      r0(3)=v3
      r1(3)=v0
      r2(3)=Rin
      r3(3)=Rout
    ELSEIF(face2 .EQ. 1)THEN
      t=1.0_d_t/(omega_local%x2-omega_local%x3)
      Rout=v2
      Rin=Rout-t*omega
      ! Subcell 1
      incoming_face(1)=2
      outgoing_face(1)=1
      r0(1)=v0
      r1(1)=v3
      r2(1)=Rin
      r3(1)=Rout
      ! Subcell 2
      incoming_face(2)=2
      outgoing_face(2)=0
      r0(2)=v3
      r1(2)=v1
      r2(2)=Rin
      r3(2)=Rout
      ! Subcell 3
      incoming_face(3)=2
      outgoing_face(3)=3
      r0(3)=v1
      r1(3)=v0
      r2(3)=Rin
      r3(3)=Rout
    ELSEIF(face3 .EQ. 1)THEN
      t=1.0_d_t/(omega_local%x3)
      Rout=v3
      Rin=Rout-t*omega
      ! Subcell 1
      incoming_face(1)=3
      outgoing_face(1)=2
      r0(1)=v0
      r1(1)=v1
      r2(1)=Rin
      r3(1)=Rout
      ! Subcell 2
      incoming_face(2)=3
      outgoing_face(2)=0
      r0(2)=v1
      r1(2)=v2
      r2(2)=Rin
      r3(2)=Rout
      ! Subcell 3
      incoming_face(3)=3
      outgoing_face(3)=1
      r0(3)=v2
      r1(3)=v0
      r2(3)=Rin
      r3(3)=Rout
    ELSE
      CALL raise_fatal_error("Cell splitting case 2 failed for cell")
    END IF

  END SUBROUTINE case2

  SUBROUTINE case3(v0,v1,v2,v3,r0,r1,r2,r3,face0,face1,face2,face3,&
        incoming_face,outgoing_face,J,omega_local,t)
    !*********************************************************************
    !
    ! Subroutine case3 generates local coordinate subcell vectors
    !
    !*********************************************************************
    ! Define temporary variables

    INTEGER(kind=li), DIMENSION(4), INTENT(inout) :: incoming_face,&
          outgoing_face
    REAL(kind=d_t), INTENT(out) :: t
    REAL(kind=d_t), DIMENSION(3,3), INTENT(in) :: J
    INTEGER(kind=1), INTENT(in) :: face0,face1,face2,face3
    TYPE(vector) :: Rin, Rout, Rin_local, Rout_local
    TYPE(vector), INTENT(in) :: omega_local, v0, v1, v2, v3
    TYPE(vector), DIMENSION(4),INTENT(inout) :: r0, r1, r2, r3

    ! Split cell and define subcell vertices r0, r1, r2, & r3

    IF((face2 .EQ. 1) .AND. (face3 .EQ. 1))THEN
      t=1.0_d_t/(omega_local%x2)
      ! Assign local coordinate incoming & outgoing vectors
      Rin_local%x1=1-t*omega_local%x1
      Rin_local%x2=zero
      Rin_local%x3=zero
      Rout_local%x1=one
      Rout_local%x2=one
      Rout_local%x3=t*omega_local%x3
      Rin=v0+J*Rin_local
      Rout=v0+J*Rout_local
      ! Subcell 1
      incoming_face(1)=2
      outgoing_face(1)=1
      r0(1)=v0
      r1(1)=v3
      r2(1)=Rin
      r3(1)=Rout
      ! Subcell 2
      incoming_face(2)=2
      outgoing_face(2)=0
      r0(2)=v3
      r1(2)=v1
      r2(2)=Rin
      r3(2)=Rout
      ! Subcell 3
      incoming_face(3)=3
      outgoing_face(3)=0
      r0(3)=v1
      r1(3)=v2
      r2(3)=Rin
      r3(3)=Rout
      ! Subcell 4
      incoming_face(4)=3
      outgoing_face(4)=1
      r0(4)=v2
      r1(4)=v0
      r2(4)=Rin
      r3(4)=Rout
    ELSEIF((face1 .EQ. 1) .AND. (face3 .EQ. 1))THEN
      t=1.0_d_t/(omega_local%x1-omega_local%x2+omega_local%x3)
      ! Assign local coordinate incoming & outgoing vectors
      Rin_local%x1=1-t*omega_local%x1
      Rin_local%x2=1-t*omega_local%x1
      Rin_local%x3=zero
      Rout_local%x1=one
      Rout_local%x2=t*omega_local%x3
      Rout_local%x3=t*omega_local%x3
      Rin=v0+J*Rin_local
      Rout=v0+J*Rout_local
      ! Subcell 1
      incoming_face(1)=3
      outgoing_face(1)=2
      r0(1)=v0
      r1(1)=v1
      r2(1)=Rin
      r3(1)=Rout
      ! Subcell 2
      incoming_face(2)=3
      outgoing_face(2)=0
      r0(2)=v1
      r1(2)=v2
      r2(2)=Rin
      r3(2)=Rout
      ! Subcell 3
      incoming_face(3)=1
      outgoing_face(3)=0
      r0(3)=v2
      r1(3)=v3
      r2(3)=Rin
      r3(3)=Rout
      ! Subcell 4
      incoming_face(4)=1
      outgoing_face(4)=2
      r0(4)=v3
      r1(4)=v0
      r2(4)=Rin
      r3(4)=Rout
    ELSEIF((face1 .EQ. 1) .AND. (face2 .EQ. 1))THEN
      t=1.0_d_t/(omega_local%x1-omega_local%x3)
      ! Assign local coordinate incoming & outgoing vectors
      Rin_local%x1=-t*omega_local%x3
      Rin_local%x2=-t*omega_local%x3
      Rin_local%x3=-t*omega_local%x3
      Rout_local%x1=one
      Rout_local%x2=t*omega_local%x2-t*omega_local%x3
      Rout_local%x3=zero
      Rin=v0+J*Rin_local
      Rout=v0+J*Rout_local
      ! Subcell 1
      incoming_face(1)=1
      outgoing_face(1)=3
      r0(1)=v0
      r1(1)=v2
      r2(1)=Rin
      r3(1)=Rout
      ! Subcell 2
      incoming_face(2)=1
      outgoing_face(2)=0
      r0(2)=v2
      r1(2)=v3
      r2(2)=Rin
      r3(2)=Rout
      ! Subcell 3
      incoming_face(3)=2
      outgoing_face(3)=0
      r0(3)=v3
      r1(3)=v1
      r2(3)=Rin
      r3(3)=Rout
      ! Subcell 4
      incoming_face(4)=2
      outgoing_face(4)=3
      r0(4)=v1
      r1(4)=v0
      r2(4)=Rin
      r3(4)=Rout
    ELSEIF((face0 .EQ. 1) .AND. (face3 .EQ. 1))THEN
      t=1.0_d_t/(omega_local%x3-omega_local%x1)
      ! Assign local coordinate incoming & outgoing vectors
      Rin_local%x1=one
      Rin_local%x2=t*omega_local%x3-t*omega_local%x2
      Rin_local%x3=zero
      Rout_local%x1=t*omega_local%x3
      Rout_local%x2=t*omega_local%x3
      Rout_local%x3=t*omega_local%x3
      Rin=v0+J*Rin_local
      Rout=v0+J*Rout_local
      ! Subcell 1
      incoming_face(1)=3
      outgoing_face(1)=2
      r0(1)=v0
      r1(1)=v1
      r2(1)=Rin
      r3(1)=Rout
      ! Subcell 2
      incoming_face(2)=0
      outgoing_face(2)=2
      r0(2)=v1
      r1(2)=v3
      r2(2)=Rin
      r3(2)=Rout
      ! Subcell 3
      incoming_face(3)=0
      outgoing_face(3)=1
      r0(3)=v3
      r1(3)=v2
      r2(3)=Rin
      r3(3)=Rout
      ! Subcell 4
      incoming_face(4)=3
      outgoing_face(4)=1
      r0(4)=v2
      r1(4)=v0
      r2(4)=Rin
      r3(4)=Rout
    ELSEIF((face0 .EQ. 1) .AND. (face2 .EQ. 1))THEN
      t=1.0_d_t/(-omega_local%x1+omega_local%x2-omega_local%x3)
      ! Assign local coordinate incoming & outgoing vectors
      Rin_local%x1=one
      Rin_local%x2=-t*omega_local%x3
      Rin_local%x3=-t*omega_local%x3
      Rout_local%x1=one+t*omega_local%x1
      Rout_local%x2=one+t*omega_local%x1
      Rout_local%x3=zero
      Rin=v0+J*Rin_local
      Rout=v0+J*Rout_local
      ! Subcell 1
      incoming_face(1)=2
      outgoing_face(1)=1
      r0(1)=v0
      r1(1)=v3
      r2(1)=Rin
      r3(1)=Rout
      ! Subcell 2
      incoming_face(2)=0
      outgoing_face(2)=1
      r0(2)=v3
      r1(2)=v2
      r2(2)=Rin
      r3(2)=Rout
      ! Subcell 3
      incoming_face(3)=0
      outgoing_face(3)=3
      r0(3)=v2
      r1(3)=v1
      r2(3)=Rin
      r3(3)=Rout
      ! Subcell 4
      incoming_face(4)=2
      outgoing_face(4)=3
      r0(4)=v1
      r1(4)=v0
      r2(4)=Rin
      r3(4)=Rout
    ELSEIF((face0 .EQ. 1) .AND. (face1 .EQ. 1))THEN
      t=1.0_d_t/(-omega_local%x2)
      ! Assign local coordinate incoming & outgoing vectors
      Rin_local%x1=one
      Rin_local%x2=one
      Rin_local%x3=-t*omega_local%x3
      Rout_local%x1=one+t*omega_local%x1
      Rout_local%x2=zero
      Rout_local%x3=zero
      Rin=v0+J*Rin_local
      Rout=v0+J*Rout_local
      ! Subcell 1
      incoming_face(1)=1
      outgoing_face(1)=3
      r0(1)=v0
      r1(1)=v2
      r2(1)=Rin
      r3(1)=Rout
      ! Subcell 2
      incoming_face(2)=0
      outgoing_face(2)=3
      r0(2)=v2
      r1(2)=v1
      r2(2)=Rin
      r3(2)=Rout
      ! Subcell 3
      incoming_face(3)=0
      outgoing_face(3)=2
      r0(3)=v1
      r1(3)=v3
      r2(3)=Rin
      r3(3)=Rout
      ! Subcell 4
      incoming_face(4)=1
      outgoing_face(4)=2
      r0(4)=v3
      r1(4)=v0
      r2(4)=Rin
      r3(4)=Rout
    ELSE
      CALL raise_fatal_error("Cell splitting case 3 failed for cell")
    END IF

  END SUBROUTINE case3

  SUBROUTINE case4(v0,v1,v2,v3,omega,n0,n1,n2,n3,r0,r1,r2,r3,face0,&
        face1,face2,face3,incoming_face,outgoing_face,omega_local,t)
    !*********************************************************************
    !
    ! Subroutine case4 generates local coordinate subcell vectors
    !
    !*********************************************************************
    ! Define temporary variables

    INTEGER(kind=li), DIMENSION(2), INTENT(inout) :: incoming_face,&
          outgoing_face
    REAL(kind=d_t), INTENT(out) :: t
    INTEGER(kind=1), INTENT(in) :: face0,face1,face2,face3
    TYPE(vector) :: Rin, Rout
    TYPE(vector), INTENT(in) :: omega, omega_local, v0, v1, v2, v3, &
          n0,n1,n2,n3
    TYPE(vector), DIMENSION(2),INTENT(inout) :: r0, r1, r2, r3

    ! Split cell and define subcell vertices r0, r1, r2, & r3

    IF((face0 .EQ. 0) .AND. ABS(omega .dot. n0) .GT. 0.0)THEN
      t=1.0_d_t/(omega_local%x1)
      Rin=v0
      Rout=Rin+t*omega
      IF((face2 .EQ. 1) .AND. (face3 .EQ. 1))THEN
        ! Subcell 1
        incoming_face(1)=2
        outgoing_face(1)=0
        r0(1)=v3
        r1(1)=v1
        r2(1)=Rin
        r3(1)=Rout
        ! Subcell 2
        incoming_face(2)=3
        outgoing_face(2)=0
        r0(2)=v1
        r1(2)=v2
        r2(2)=Rin
        r3(2)=Rout
      ELSEIF((face1 .EQ. 1) .AND. (face3 .EQ. 1))THEN
        ! Subcell 1
        incoming_face(1)=3
        outgoing_face(1)=0
        r0(1)=v1
        r1(1)=v2
        r2(1)=Rin
        r3(1)=Rout
        ! Subcell 2
        incoming_face(2)=1
        outgoing_face(2)=0
        r0(2)=v2
        r1(2)=v3
        r2(2)=Rin
        r3(2)=Rout
      ELSEIF((face1 .EQ. 1) .AND. (face2 .EQ. 1))THEN
        ! Subcell 1
        incoming_face(1)=1
        outgoing_face(1)=0
        r0(1)=v2
        r1(1)=v3
        r2(1)=Rin
        r3(1)=Rout
        ! Subcell 2
        incoming_face(2)=2
        outgoing_face(2)=0
        r0(2)=v3
        r1(2)=v1
        r2(2)=Rin
        r3(2)=Rout
      ENDIF
    ELSEIF((face1 .EQ. 0) .AND. ABS(omega .dot. n1) .GT. zero)THEN
      t=1.0_d_t/(omega_local%x2-omega_local%x1)
      Rin=v1
      Rout=Rin+t*omega
      IF((face2 .EQ. 1) .AND. (face3 .EQ. 1))THEN
        ! Subcell 1
        incoming_face(1)=3
        outgoing_face(1)=1
        r0(1)=v2
        r1(1)=v0
        r2(1)=Rin
        r3(1)=Rout
        ! Subcell 2
        incoming_face(2)=2
        outgoing_face(2)=1
        r0(2)=v0
        r1(2)=v3
        r2(2)=Rin
        r3(2)=Rout
      ELSEIF((face0 .EQ. 1) .AND. (face3 .EQ. 1))THEN
        ! Subcell 1
        incoming_face(1)=0
        outgoing_face(1)=1
        r0(1)=v3
        r1(1)=v2
        r2(1)=Rin
        r3(1)=Rout
        ! Subcell 2
        incoming_face(2)=3
        outgoing_face(2)=1
        r0(2)=v2
        r1(2)=v0
        r2(2)=Rin
        r3(2)=Rout
      ELSEIF((face0 .EQ. 1) .AND. (face2 .EQ. 1))THEN
        ! Subcell 1
        incoming_face(1)=2
        outgoing_face(1)=1
        r0(1)=v0
        r1(1)=v3
        r2(1)=Rin
        r3(1)=Rout
        ! Subcell 2
        incoming_face(2)=0
        outgoing_face(2)=1
        r0(2)=v3
        r1(2)=v2
        r2(2)=Rin
        r3(2)=Rout
      ENDIF
    ELSEIF((face2 .EQ. 0) .AND. ABS(omega .dot. n2) .GT. 0.0)THEN
      t=1.0_d_t/(omega_local%x3-omega_local%x2)
      Rin=v2
      Rout=Rin+t*omega
      IF((face1 .EQ. 1) .AND. (face3 .EQ. 1))THEN
        ! Subcell 1
        incoming_face(1)=1
        outgoing_face(1)=2
        r0(1)=v3
        r1(1)=v0
        r2(1)=Rin
        r3(1)=Rout
        ! Subcell 2
        incoming_face(2)=3
        outgoing_face(2)=2
        r0(2)=v0
        r1(2)=v1
        r2(2)=Rin
        r3(2)=Rout
      ELSEIF((face0 .EQ. 1) .AND. (face3 .EQ. 1))THEN
        ! Subcell 1
        incoming_face(1)=3
        outgoing_face(1)=2
        r0(1)=v0
        r1(1)=v1
        r2(1)=Rin
        r3(1)=Rout
        ! Subcell 2
        incoming_face(2)=0
        outgoing_face(2)=2
        r0(2)=v1
        r1(2)=v3
        r2(2)=Rin
        r3(2)=Rout
      ELSEIF((face0 .EQ. 1) .AND. (face1 .EQ. 1))THEN
        ! Subcell 1
        incoming_face(1)=0
        outgoing_face(1)=2
        r0(1)=v1
        r1(1)=v3
        r2(1)=Rin
        r3(1)=Rout
        ! Subcell 2
        incoming_face(2)=1
        outgoing_face(2)=2
        r0(2)=v3
        r1(2)=v0
        r2(2)=Rin
        r3(2)=Rout
      ENDIF
    ELSEIF((face3 .EQ. 0) .AND. ABS(omega .dot. n3) .GT. 0.0)THEN
      t=1.0_d_t/(-1.0_d_t*omega_local%x3)
      Rin=v3
      Rout=Rin+t*omega
      IF((face1 .EQ. 1) .AND. (face2 .EQ. 1))THEN
        ! Subcell 1
        incoming_face(1)=2
        outgoing_face(1)=3
        r0(1)=v1
        r1(1)=v0
        r2(1)=Rin
        r3(1)=Rout
        ! Subcell 2
        incoming_face(2)=1
        outgoing_face(2)=3
        r0(2)=v0
        r1(2)=v2
        r2(2)=Rin
        r3(2)=Rout
      ELSEIF((face0 .EQ. 1) .AND. (face2 .EQ. 1))THEN
        ! Subcell 1
        incoming_face(1)=0
        outgoing_face(1)=3
        r0(1)=v2
        r1(1)=v1
        r2(1)=Rin
        r3(1)=Rout
        ! Subcell 2
        incoming_face(2)=2
        outgoing_face(2)=3
        r0(2)=v1
        r1(2)=v0
        r2(2)=Rin
        r3(2)=Rout
      ELSEIF((face0 .EQ. 1) .AND. (face1 .EQ. 1))THEN
        ! Subcell 1
        incoming_face(1)=1
        outgoing_face(1)=3
        r0(1)=v0
        r1(1)=v2
        r2(1)=Rin
        r3(1)=Rout
        ! Subcell 2
        incoming_face(2)=0
        outgoing_face(2)=3
        r0(2)=v2
        r1(2)=v1
        r2(2)=Rin
        r3(2)=Rout
      ENDIF
    ELSE
      CALL raise_fatal_error("Cell splitting case 4 failed for cell")
    END IF

  END SUBROUTINE case4

  SUBROUTINE case5(v0,v1,v2,v3,omega,n0,n1,n2,n3,r0,r1,r2,r3,face0,&
        face1,face2,face3,incoming_face,outgoing_face,omega_local,t)
    !*********************************************************************
    !
    ! Subroutine case5 generates local coordinate subcell vectors
    !
    !*********************************************************************
    ! Define temporary variables
    INTEGER(kind=li), DIMENSION(2), INTENT(inout) :: incoming_face,&
          outgoing_face
    REAL(kind=d_t), INTENT(out) :: t
    INTEGER(kind=1), INTENT(in) :: face0,face1,face2,face3
    TYPE(vector) :: Rin, Rout
    TYPE(vector), INTENT(in) :: omega, omega_local, v0, v1, v2, v3, &
          n0,n1,n2,n3
    TYPE(vector), DIMENSION(2),INTENT(inout) :: r0, r1, r2, r3

    ! Split cell and define subcell vertices r0, r1, r2, & r3

    IF(face0 .EQ. 1)THEN
      t=1.0_d_t/(-1.0_d_t*omega_local%x1)
      Rout=v0
      Rin=Rout-t*omega
      IF(ABS(omega .dot. n1) .LE. 1.0D-16)THEN
        ! Subcell 1
        incoming_face(1)=0
        outgoing_face(1)=3
        r0(1)=v2
        r1(1)=v1
        r2(1)=Rin
        r3(1)=Rout
        ! Subcell 2
        incoming_face(2)=0
        outgoing_face(2)=2
        r0(2)=v1
        r1(2)=v3
        r2(2)=Rin
        r3(2)=Rout
      ELSEIF(ABS(omega .dot. n2) .LE. 1.0D-16)THEN
        ! Subcell 1
        incoming_face(1)=0
        outgoing_face(1)=1
        r0(1)=v3
        r1(1)=v2
        r2(1)=Rin
        r3(1)=Rout
        ! Subcell 2
        incoming_face(2)=0
        outgoing_face(2)=3
        r0(2)=v2
        r1(2)=v1
        r2(2)=Rin
        r3(2)=Rout
      ELSEIF(ABS(omega .dot. n3) .LE. 1.0D-16)THEN
        ! Subcell 1
        incoming_face(1)=0
        outgoing_face(1)=2
        r0(1)=v1
        r1(1)=v3
        r2(1)=Rin
        r3(1)=Rout
        ! Subcell 2
        incoming_face(2)=0
        outgoing_face(2)=1
        r0(2)=v3
        r1(2)=v2
        r2(2)=Rin
        r3(2)=Rout
      ENDIF
    ELSEIF(face1 .EQ. 1)THEN
      t=1.0_d_t/(omega_local%x1-omega_local%x2)
      Rout=v1
      Rin=Rout-t*omega
      IF(ABS(omega .dot. n0) .LE. 1.0D-16)THEN
        ! Subcell 1
        incoming_face(1)=1
        outgoing_face(1)=2
        r0(1)=v3
        r1(1)=v0
        r2(1)=Rin
        r3(1)=Rout
        ! Subcell 2
        incoming_face(2)=1
        outgoing_face(2)=3
        r0(2)=v0
        r1(2)=v2
        r2(2)=Rin
        r3(2)=Rout
      ELSEIF(ABS(omega .dot. n2) .LE. 1.0D-16)THEN
        ! Subcell 1
        incoming_face(1)=1
        outgoing_face(1)=3
        r0(1)=v0
        r1(1)=v2
        r2(1)=Rin
        r3(1)=Rout
        ! Subcell 2
        incoming_face(2)=1
        outgoing_face(2)=0
        r0(2)=v2
        r1(2)=v3
        r2(2)=Rin
        r3(2)=Rout
      ELSEIF(ABS(omega .dot. n3) .LE. 1.0D-16)THEN
        ! Subcell 1
        incoming_face(1)=1
        outgoing_face(1)=0
        r0(1)=v2
        r1(1)=v3
        r2(1)=Rin
        r3(1)=Rout
        ! Subcell 2
        incoming_face(2)=1
        outgoing_face(2)=2
        r0(2)=v3
        r1(2)=v0
        r2(2)=Rin
        r3(2)=Rout
      ENDIF
    ELSEIF(face2 .EQ. 1)THEN
      t=1.0_d_t/(omega_local%x2-omega_local%x3)
      Rout=v2
      Rin=Rout-t*omega
      IF(ABS(omega .dot. n0) .LE. 1.0D-16)THEN
        ! Subcell 1
        incoming_face(1)=2
        outgoing_face(1)=3
        r0(1)=v1
        r1(1)=v0
        r2(1)=Rin
        r3(1)=Rout
        ! Subcell 2
        incoming_face(2)=2
        outgoing_face(2)=1
        r0(2)=v0
        r1(2)=v3
        r2(2)=Rin
        r3(2)=Rout
      ELSEIF(ABS(omega .dot. n1) .LE. 1.0D-16)THEN
        ! Subcell 1
        incoming_face(1)=2
        outgoing_face(1)=0
        r0(1)=v3
        r1(1)=v1
        r2(1)=Rin
        r3(1)=Rout
        ! Subcell 2
        incoming_face(2)=2
        outgoing_face(2)=3
        r0(2)=v1
        r1(2)=v0
        r2(2)=Rin
        r3(2)=Rout
      ELSEIF(ABS(omega .dot. n3) .LE. 1.0D-16)THEN
        ! Subcell 1
        incoming_face(1)=2
        outgoing_face(1)=1
        r0(1)=v0
        r1(1)=v3
        r2(1)=Rin
        r3(1)=Rout
        ! Subcell 2
        incoming_face(2)=2
        outgoing_face(2)=0
        r0(2)=v3
        r1(2)=v1
        r2(2)=Rin
        r3(2)=Rout
      ENDIF
    ELSEIF(face3 .EQ. 1)THEN
      t=1.0_d_t/(omega_local%x3)
      Rout=v3
      Rin=Rout-t*omega
      IF(ABS(omega .dot. n0) .LE. 1.0D-16)THEN
        ! Subcell 1
        incoming_face(1)=3
        outgoing_face(1)=1
        r0(1)=v2
        r1(1)=v0
        r2(1)=Rin
        r3(1)=Rout
        ! Subcell 2
        incoming_face(2)=3
        outgoing_face(2)=2
        r0(2)=v0
        r1(2)=v1
        r2(2)=Rin
        r3(2)=Rout
      ELSEIF(ABS(omega .dot. n1) .LE. 1.0D-16)THEN
        ! Subcell 1
        incoming_face(1)=3
        outgoing_face(1)=2
        r0(1)=v0
        r1(1)=v1
        r2(1)=Rin
        r3(1)=Rout
        ! Subcell 2
        incoming_face(2)=3
        outgoing_face(2)=0
        r0(2)=v1
        r1(2)=v2
        r2(2)=Rin
        r3(2)=Rout
      ELSEIF(ABS(omega .dot. n2) .LE. 1.0D-16)THEN
        ! Subcell 1
        incoming_face(1)=3
        outgoing_face(1)=0
        r0(1)=v1
        r1(1)=v2
        r2(1)=Rin
        r3(1)=Rout
        ! Subcell 2
        incoming_face(2)=3
        outgoing_face(2)=1
        r0(2)=v2
        r1(2)=v0
        r2(2)=Rin
        r3(2)=Rout
      ENDIF
    ELSE
      CALL raise_fatal_error("Cell splitting case 5 failed for cell")
    END IF

  END SUBROUTINE case5

  SUBROUTINE case6(v0,v1,v2,v3,omega,n0,n1,n2,n3,r0,r1,r2,r3,face0,&
        face1,face2,face3,incoming_face,outgoing_face,omega_local,t)
    !*********************************************************************
    !
    ! Subroutine case6 generates local coordinate subcell vectors
    !
    !*********************************************************************
    ! Define temporary variables

    INTEGER(kind=li), DIMENSION(1), INTENT(inout) :: incoming_face,&
          outgoing_face
    REAL(kind=d_t), INTENT(out) :: t
    INTEGER(kind=1), INTENT(in) :: face0,face1,face2,face3
    TYPE(vector) :: Rin, Rout
    TYPE(vector), INTENT(in) :: omega, omega_local, v0, v1, v2, v3, &
          n0,n1,n2,n3
    TYPE(vector), DIMENSION(1), INTENT(inout) :: r0, r1, r2, r3

    ! Split cell and define subcell vertices r0, r1, r2, & r3

    IF((face0 .EQ. 0) .AND. ABS(omega .dot. n0) .GT. 0.0)THEN
      t=1.0_d_t/(omega_local%x1)
      Rin=v0
      Rout=Rin+t*omega
      IF(face1 .EQ. 1)THEN
        ! Subcell 1
        incoming_face(1)=1
        outgoing_face(1)=0
        r0(1)=v2
        r1(1)=v3
        r2(1)=Rin
        r3(1)=Rout
      ELSEIF(face2 .EQ. 1)THEN
        ! Subcell 1
        incoming_face(1)=2
        outgoing_face(1)=0
        r0(1)=v3
        r1(1)=v1
        r2(1)=Rin
        r3(1)=Rout
      ELSEIF(face3 .EQ. 1)THEN
        ! Subcell 1
        incoming_face(1)=3
        outgoing_face(1)=0
        r0(1)=v1
        r1(1)=v2
        r2(1)=Rin
        r3(1)=Rout
      ENDIF
    ELSEIF((face1 .EQ. 0) .AND. ABS(omega .dot. n1) .GT. 0.0)THEN
      t=1.0_d_t/(omega_local%x2-omega_local%x1)
      Rin=v1
      Rout=Rin+t*omega
      IF(face0 .EQ. 1)THEN
        ! Subcell 1
        incoming_face(1)=0
        outgoing_face(1)=1
        r0(1)=v3
        r1(1)=v2
        r2(1)=Rin
        r3(1)=Rout
      ELSEIF(face2 .EQ. 1)THEN
        ! Subcell 1
        incoming_face(1)=2
        outgoing_face(1)=1
        r0(1)=v0
        r1(1)=v3
        r2(1)=Rin
        r3(1)=Rout
      ELSEIF(face3 .EQ. 1)THEN
        ! Subcell 1
        incoming_face(1)=3
        outgoing_face(1)=1
        r0(1)=v2
        r1(1)=v0
        r2(1)=Rin
        r3(1)=Rout
      ENDIF
    ELSEIF((face2 .EQ. 0) .AND. ABS(omega .dot. n2) .GT. 0.0)THEN
      t=1.0_d_t/(omega_local%x3-omega_local%x2)
      Rin=v2
      Rout=Rin+t*omega
      IF(face0 .EQ. 1)THEN
        ! Subcell 1
        incoming_face(1)=0
        outgoing_face(1)=2
        r0(1)=v1
        r1(1)=v3
        r2(1)=Rin
        r3(1)=Rout
      ELSEIF(face1 .EQ. 1)THEN
        ! Subcell 1
        incoming_face(1)=1
        outgoing_face(1)=2
        r0(1)=v3
        r1(1)=v0
        r2(1)=Rin
        r3(1)=Rout
      ELSEIF(face3 .EQ. 1)THEN
        ! Subcell 1
        incoming_face(1)=3
        outgoing_face(1)=2
        r0(1)=v0
        r1(1)=v1
        r2(1)=Rin
        r3(1)=Rout
      ENDIF
    ELSEIF((face3 .EQ. 0) .AND. ABS(omega .dot. n3) .GT. 0.0)THEN
      t=1.0_d_t/(-1.0_d_t*omega_local%x3)
      Rin=v3
      Rout=Rin+t*omega
      IF(face0 .EQ. 1)THEN
        ! Subcell 1
        incoming_face(1)=0
        outgoing_face(1)=3
        r0(1)=v2
        r1(1)=v1
        r2(1)=Rin
        r3(1)=Rout
      ELSEIF(face1 .EQ. 1)THEN
        ! Subcell 1
        incoming_face(1)=1
        outgoing_face(1)=3
        r0(1)=v0
        r1(1)=v2
        r2(1)=Rin
        r3(1)=Rout
      ELSEIF(face2 .EQ. 1)THEN
        ! Subcell 1
        incoming_face(1)=2
        outgoing_face(1)=3
        r0(1)=v1
        r1(1)=v0
        r2(1)=Rin
        r3(1)=Rout
      ENDIF
    ELSE
      CALL raise_fatal_error("Cell splitting case 6 failed for cell")
    END IF

  END SUBROUTINE case6

END MODULE cell_splitting_module
