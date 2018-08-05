MODULE transport_kernel_module_CCE
  !***********************************************************************
  !
  ! Transport kernel module performs the transport calculation in which
  ! a tetrahedra is divided into 'canonical' or 'characteristic'
  ! tetrahedra and the outgoing angular flux is computed based on the
  ! incoming angular flux by integrating over the faces and assuming
  ! a characteristic flux shape and an arbitrary polynomial expansion
  ! of the source and flux functions.
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
  USE termination_module
  USE general_utility_module

  IMPLICIT NONE

CONTAINS

  SUBROUTINE transport_kernel_CCE(sigmat,q_moments,vol_moment,face_moment,       &
        LL,U,Lf,Uf,i,t,subcells,r0,r1,r2,r3,           &
        face_known,incoming_face,outgoing_face,J,J_inv,&
        upstream_moments)
    !*********************************************************************
    !
    ! Subroutine transport cell computes outgoing angular face and
    ! volume moments based on incoming angular face fluxes and source
    ! moments
    !
    !*********************************************************************

    ! Pass arguments

    REAL(kind=d_t), INTENT(in) :: sigmat, t
    REAL(kind=d_t), DIMENSION(num_moments_v), INTENT(in) :: q_moments
    REAL(kind=d_t), DIMENSION(num_moments_v), INTENT(inout) :: vol_moment
    REAL(kind=d_t), DIMENSION(num_moments_f,0:3), INTENT(inout) :: face_moment
    REAL(kind=d_t), DIMENSION(num_moments_v,num_moments_v) :: LL, U
    REAL(kind=d_t), DIMENSION(num_moments_f,num_moments_f) :: Lf, Uf
    INTEGER(kind=li), INTENT(in) :: i, subcells
    TYPE(vector), DIMENSION(subcells), INTENT(in) :: r0, r1, r2, r3
    INTEGER(kind=1), DIMENSION(0:3,num_cells), INTENT(inout) :: face_known
    INTEGER(kind=li), DIMENSION(subcells), INTENT(in) :: incoming_face, outgoing_face
    REAL(kind=d_t), DIMENSION(3,3), INTENT(in) :: J, J_inv
    REAL(kind=d_t), DIMENSION(num_moments_f,subcells), INTENT(in) :: upstream_moments


    ! Define temporary variables

    INTEGER(kind=li) :: alloc_stat, ii, l, f, adjcnt_cell
    REAL(kind=d_t) :: det_Jup, det_Js, e
    REAL(kind=d_t), DIMENSION(2) :: af
    REAL(kind=d_t), DIMENSION(3) :: a_temp, a
    REAL(kind=d_t), DIMENSION(2,2) :: bf
    REAL(kind=d_t), DIMENSION(3,2) :: JFup, JFdown, JF, Jsf
    REAL(kind=d_t), DIMENSION(2,3) :: JFup_inv, JF_inv
    REAL(kind=d_t), DIMENSION(3,3) :: Jup, Jup_inv, Js, Js_inv, b
    REAL(kind=d_t), DIMENSION(num_moments_v) :: q_expansion
    REAL(kind=d_t), DIMENSION(:), ALLOCATABLE :: area, volume, &
          subcell_upstream_moments, transformed_moments, &
          transformed_flux, incoming_flux, outgoing_moments,cell_source, &
          subcell_source, subcell_flux, projected_temp, proj_moments, &
          proj_flux_moments, face_angular_mom, face_cell_temp, y
    REAL(kind=d_t), DIMENSION(:,:), ALLOCATABLE :: &
          subcell_source_moments, projected_moments,&
          face_angular, projected_flux
    TYPE(vector) :: R0down, R0up, R0F, v0, v1, v2, v3
    TYPE(vector), DIMENSION(0:3) :: v

    ! Allocate variables

    ALLOCATE(area(subcells),volume(subcells),&
          subcell_upstream_moments(num_moments_f),&
          incoming_flux(num_moments_f),&
          transformed_moments(num_moments_f),&
          transformed_flux(num_moments_f),&
          cell_source(num_moments_v),&
          outgoing_moments(num_moments_f),&
          subcell_source(num_moments_v),&
          subcell_flux(num_moments_v),&
          projected_moments(subcells,num_moments_f),&
          projected_temp(num_moments_f),&
          proj_moments(num_moments_f),&
          projected_flux(subcells,num_moments_v),&
          proj_flux_moments(num_moments_v),&
          face_angular_mom(num_moments_f),&
          face_cell_temp(num_moments_v),&
          face_angular(0:3,num_moments_v),&
          subcell_source_moments(subcells,num_moments_v),&
          y(num_moments_v),&
          stat=alloc_stat);area=0.0_d_t;volume=0.0_d_t;&
          subcell_upstream_moments=0.0_d_t;incoming_flux=0.0_d_t;&
          transformed_flux=0.0_d_t;transformed_moments=0.0_d_t;&
          cell_source=0.0_d_t;subcell_source=0.0_d_t;&
          outgoing_moments=0.0_d_t;projected_temp=0.0_d_t;&
          proj_moments=0.0_d_t;projected_moments=0.0_d_t;&
          face_angular_mom=0.0_d_t;face_angular=0.0_d_t;&
          face_cell_temp=0.0_d_t;subcell_source_moments=0.0_d_t;y=0.0_d_t
    IF(alloc_stat /= 0) CALL stop_thor(2_li)

    ! Compute optical thickness

    e=sigmat*t

    ! Loop over subcells and transport through each

    DO ii=1, subcells

      ! Transform upstream cell outgoing face moments into incoming face moments

      IF(adjacency_list(i,incoming_face(ii))%cell /= 0)THEN
        adjcnt_cell=adjacency_list(i,incoming_face(ii))%cell
        v0=vertices(cells(adjcnt_cell)%R(0))%v
        v1=vertices(cells(adjcnt_cell)%R(1))%v
        v2=vertices(cells(adjcnt_cell)%R(2))%v
        v3=vertices(cells(adjcnt_cell)%R(3))%v
        CALL cell_jacobian(v0,v1,v2,v3,Jup,det_Jup,Jup_inv)
        IF(adjacency_list(i,incoming_face(ii))%face == 0)THEN
          JFup(1,1)=Jup(1,2)
          JFup(2,1)=Jup(2,2)
          JFup(3,1)=Jup(3,2)
          JFup(1,2)=Jup(1,3)
          JFup(2,2)=Jup(2,3)
          JFup(3,2)=Jup(3,3)
          R0up=vertices(cells(adjacency_list(i,&
                incoming_face(ii))%cell)%R(1))%v
        ELSEIF(adjacency_list(i,incoming_face(ii))%face == 1)THEN
          JFup(1,1)=Jup(1,1)+Jup(1,2)
          JFup(2,1)=Jup(2,1)+Jup(2,2)
          JFup(3,1)=Jup(3,1)+Jup(3,2)
          JFup(1,2)=Jup(1,3)
          JFup(2,2)=Jup(2,3)
          JFup(3,2)=Jup(3,3)
          R0up=vertices(cells(adjacency_list(i,&
                incoming_face(ii))%cell)%R(0))%v
        ELSEIF(adjacency_list(i,incoming_face(ii))%face == 2)THEN
          JFup(1,1)=Jup(1,1)
          JFup(2,1)=Jup(2,1)
          JFup(3,1)=Jup(3,1)
          JFup(1,2)=Jup(1,2)+Jup(1,3)
          JFup(2,2)=Jup(2,2)+Jup(2,3)
          JFup(3,2)=Jup(3,2)+Jup(3,3)
          R0up=vertices(cells(adjacency_list(i,&
                incoming_face(ii))%cell)%R(0))%v
        ELSEIF(adjacency_list(i,incoming_face(ii))%face == 3)THEN
          JFup(1,1)=Jup(1,1)
          JFup(2,1)=Jup(2,1)
          JFup(3,1)=Jup(3,1)
          JFup(1,2)=Jup(1,2)
          JFup(2,2)=Jup(2,2)
          JFup(3,2)=Jup(3,2)
          R0up=vertices(cells(adjacency_list(i,&
                incoming_face(ii))%cell)%R(0))%v
        ELSE
          CALL stop_thor(3_li)
        END IF

        IF(incoming_face(ii) == 0)THEN
          JFdown(1,1)=J(1,2)
          JFdown(2,1)=J(2,2)
          JFdown(3,1)=J(3,2)
          JFdown(1,2)=J(1,3)
          JFdown(2,2)=J(2,3)
          JFdown(3,2)=J(3,3)
          R0down=vertices(cells(i)%R(1))%v
        ELSEIF(incoming_face(ii) == 1)THEN
          JFdown(1,1)=J(1,1)+J(1,2)
          JFdown(2,1)=J(2,1)+J(2,2)
          JFdown(3,1)=J(3,1)+J(3,2)
          JFdown(1,2)=J(1,3)
          JFdown(2,2)=J(2,3)
          JFdown(3,2)=J(3,3)
          R0down=vertices(cells(i)%R(0))%v
        ELSEIF(incoming_face(ii) == 2)THEN
          JFdown(1,1)=J(1,1)
          JFdown(2,1)=J(2,1)
          JFdown(3,1)=J(3,1)
          JFdown(1,2)=J(1,2)+J(1,3)
          JFdown(2,2)=J(2,2)+J(2,3)
          JFdown(3,2)=J(3,2)+J(3,3)
          R0down=vertices(cells(i)%R(0))%v
        ELSEIF(incoming_face(ii) == 3)THEN
          JFdown(1,1)=J(1,1)
          JFdown(2,1)=J(2,1)
          JFdown(3,1)=J(3,1)
          JFdown(1,2)=J(1,2)
          JFdown(2,2)=J(2,2)
          JFdown(3,2)=J(3,2)
          R0down=vertices(cells(i)%R(0))%v
        ELSE
          CALL stop_thor(4_li)
        END IF

        CALL invert_face_jacobian(JFup,JFup_inv)

        a_temp(1)=R0down%x1-R0up%x1
        a_temp(2)=R0down%x2-R0up%x2
        a_temp(3)=R0down%x3-R0up%x3
        af=MATMUL(JFup_inv,a_temp)
        bf=MATMUL(JFup_inv,JFdown)

        DO l=1, num_moments_f
          subcell_upstream_moments(l)=upstream_moments(l,ii)
        END DO

        CALL transform_incoming_moments(Lf,Uf,af,bf,subcell_upstream_moments,&
              transformed_moments,transformed_flux)

      ELSE
        DO l=1, num_moments_f
          transformed_moments(l)=upstream_moments(l,ii)
        END DO
        CALL transform_boundary_moments(Lf,Uf,transformed_moments,transformed_flux)
      END IF

      DO l=1, num_moments_f
        face_moment(l,incoming_face(ii))=transformed_moments(l)
      END DO

      ! Transform incoming cell face moments into subcell face moments

      v(0)=r0(ii)

      v(1)=r1(ii)

      v(2)=r2(ii)

      v(3)=r3(ii)

      CALL subcell_jacobian(v,Js,det_Js,Js_inv)

      volume(ii)=(1.0_d_t/6.0_d_t)*det_Js

      Jsf(1,1)=Js(1,1)
      Jsf(2,1)=Js(2,1)
      Jsf(3,1)=Js(3,1)
      Jsf(1,2)=Js(1,2)
      Jsf(2,2)=Js(2,2)
      Jsf(3,2)=Js(3,2)

      IF(incoming_face(ii) == 0)THEN
        JF(1,1)=J(1,2)
        JF(2,1)=J(2,2)
        JF(3,1)=J(3,2)
        JF(1,2)=J(1,3)
        JF(2,2)=J(2,3)
        JF(3,2)=J(3,3)
        R0F=vertices(cells(i)%R(1))%v
      ELSEIF(incoming_face(ii) == 1)THEN
        JF(1,1)=J(1,1)+J(1,2)
        JF(2,1)=J(2,1)+J(2,2)
        JF(3,1)=J(3,1)+J(3,2)
        JF(1,2)=J(1,3)
        JF(2,2)=J(2,3)
        JF(3,2)=J(3,3)
        R0F=vertices(cells(i)%R(0))%v
      ELSEIF(incoming_face(ii) == 2)THEN
        JF(1,1)=J(1,1)
        JF(2,1)=J(2,1)
        JF(3,1)=J(3,1)
        JF(1,2)=J(1,2)+J(1,3)
        JF(2,2)=J(2,2)+J(2,3)
        JF(3,2)=J(3,2)+J(3,3)
        R0F=vertices(cells(i)%R(0))%v
      ELSEIF(incoming_face(ii) == 3)THEN
        JF(1,1)=J(1,1)
        JF(2,1)=J(2,1)
        JF(3,1)=J(3,1)
        JF(1,2)=J(1,2)
        JF(2,2)=J(2,2)
        JF(3,2)=J(3,2)
        R0F=vertices(cells(i)%R(0))%v
      ELSE
        CALL stop_thor(5_li)
      END IF

      CALL invert_face_jacobian(JF,JF_inv)

      a_temp(1)=r0(ii)%x1-R0F%x1
      a_temp(2)=r0(ii)%x2-R0F%x2
      a_temp(3)=r0(ii)%x3-R0F%x3
      af=MATMUL(JF_inv,a_temp)
      bf=MATMUL(JF_inv,Jsf)

      CALL incoming_cell_subcell_project(Lf,Uf,af,bf,transformed_flux,incoming_flux)

      ! Compute cell source expansion coefficients based on source moments

      CALL cell_source_expansion(q_moments,LL,U,q_expansion)

      ! Project cell source moment into subcell moments in subcell system

      a_temp(1)=r0(ii)%x1-vertices(cells(i)%R(0))%v%x1
      a_temp(2)=r0(ii)%x2-vertices(cells(i)%R(0))%v%x2
      a_temp(3)=r0(ii)%x3-vertices(cells(i)%R(0))%v%x3
      a=MATMUL(J_inv,a_temp)
      b=MATMUL(J_inv,Js)

      DO l=1, num_moments_v
        cell_source(l)=q_expansion(l)
      END DO

      CALL source_projection(LL,U,a,b,cell_source,subcell_source)

      ! Compute outgoing face moments in subcell with characteristic relation

      ! Project outgoing subcell moments into cell outgoing moments

      Jsf(1,1)=Js(1,1)
      Jsf(2,1)=Js(2,1)
      Jsf(3,1)=Js(3,1)
      Jsf(1,2)=Js(1,2)+Js(1,3)
      Jsf(2,2)=Js(2,2)+Js(2,3)
      Jsf(3,2)=Js(3,2)+Js(3,3)

      IF(outgoing_face(ii) == 0)THEN
        JF(1,1)=J(1,2)
        JF(2,1)=J(2,2)
        JF(3,1)=J(3,2)
        JF(1,2)=J(1,3)
        JF(2,2)=J(2,3)
        JF(3,2)=J(3,3)
        R0F=vertices(cells(i)%R(1))%v
      ELSEIF(outgoing_face(ii) == 1)THEN
        JF(1,1)=J(1,1)+J(1,2)
        JF(2,1)=J(2,1)+J(2,2)
        JF(3,1)=J(3,1)+J(3,2)
        JF(1,2)=J(1,3)
        JF(2,2)=J(2,3)
        JF(3,2)=J(3,3)
        R0F=vertices(cells(i)%R(0))%v
      ELSEIF(outgoing_face(ii) == 2)THEN
        JF(1,1)=J(1,1)
        JF(2,1)=J(2,1)
        JF(3,1)=J(3,1)
        JF(1,2)=J(1,2)+J(1,3)
        JF(2,2)=J(2,2)+J(2,3)
        JF(3,2)=J(3,2)+J(3,3)
        R0F=vertices(cells(i)%R(0))%v
      ELSEIF(outgoing_face(ii) == 3)THEN
        JF(1,1)=J(1,1)
        JF(2,1)=J(2,1)
        JF(3,1)=J(3,1)
        JF(1,2)=J(1,2)
        JF(2,2)=J(2,2)
        JF(3,2)=J(3,2)
        R0F=vertices(cells(i)%R(0))%v
      ELSE
        CALL stop_thor(6_li)
      END IF

      CALL invert_face_jacobian(JF,JF_inv)

      a_temp(1)=r0(ii)%x1-R0F%x1
      a_temp(2)=r0(ii)%x2-R0F%x2
      a_temp(3)=r0(ii)%x3-R0F%x3
      af=MATMUL(JF_inv,a_temp)
      bf=MATMUL(JF_inv,Jsf)

      a_temp(1)=r0(ii)%x1-vertices(cells(i)%R(0))%v%x1
      a_temp(2)=r0(ii)%x2-vertices(cells(i)%R(0))%v%x2
      a_temp(3)=r0(ii)%x3-vertices(cells(i)%R(0))%v%x3
      a=MATMUL(J_inv,a_temp)
      b=MATMUL(J_inv,Js)

      CALL characteristic_solver(i,t,e,a,b,af,bf,incoming_flux,&
            subcell_source,outgoing_moments,subcell_flux)

      DO l=1, num_moments_f
        projected_moments(ii,l)=outgoing_moments(l)
      END DO

      ! Project subcell flux moments into cell flux moments

      DO l=1, num_moments_v
        projected_flux(ii,l)= subcell_flux(l)
      END DO
    END DO

    ! Compute area weighted outgoing face angular moments and update
    ! face_known array

    DO ii=1, subcells
      area(ii)=(0.5_d_t)*vmag((r1(ii)-r0(ii)) .cross. (r3(ii)-r1(ii)))
      DO l=1, num_moments_f

        face_moment(l,outgoing_face(ii))=&
              face_moment(l,outgoing_face(ii))+(area(ii)/((0.5_d_t)*&
              vmag(outward_normal(i,outgoing_face(ii)))))*&
              projected_moments(ii,l)
      END DO
      IF(adjacency_list(i,outgoing_face(ii))%cell /= 0)THEN
        face_known(adjacency_list(i,outgoing_face(ii))%face,&
              adjacency_list(i,outgoing_face(ii))%cell)=1
      END IF
    END DO

    ! Transform face flux into face moments

    DO f=0, 3
      DO l=1, num_moments_f
        face_angular_mom(l)=face_moment(l,f)
      END DO

      CALL face_moment_transformation(Lf,Uf,f,face_angular_mom,face_cell_temp)

      DO l=1, num_moments_v
        face_angular(f,l)=face_cell_temp(l)
      END DO
    END DO

    ! Compute volume weighted angular flux moments over cell

    DO l=1, num_moments_v
      DO ii=1, subcells
        y(l)=y(l)+(volume(ii)/cells(i)%volume)*projected_flux(ii,l)
      END DO
    END DO

    DO l=1, num_moments_v
      vol_moment(l)=y(l)
    END DO

    DEALLOCATE(area,volume,subcell_upstream_moments,incoming_flux,&
          transformed_moments,transformed_flux,cell_source,&
          subcell_source,outgoing_moments,projected_moments,&
          projected_temp,proj_moments,projected_flux,&
          proj_flux_moments,face_angular,subcell_flux,&
          subcell_source_moments)

  END SUBROUTINE transport_kernel_CCE

  SUBROUTINE transform_incoming_moments(Lf,Uf,af,bf,upstream_moments,transformed_moments,transformed_flux)
    !*********************************************************************
    !
    ! Subroutine transform incoming moments transforms the upstream cell
    ! face outgoing moments into downstream incoming cell face moments
    !
    !*********************************************************************

    ! Define variables

    INTEGER(kind=li) :: alloc_stat, l, q, i1, i2, iup1, iup2, m11, &
          m12, m21, m22
    REAL(kind=d_t), DIMENSION(num_moments_f), INTENT(in) :: &
          upstream_moments
    REAL(kind=d_t), DIMENSION(num_moments_f), INTENT(out) :: &
          transformed_moments, transformed_flux
    REAL(kind=d_t), DIMENSION(2), INTENT(in) :: af
    REAL(kind=d_t), DIMENSION(2,2), INTENT(in) :: bf
    REAL(kind=d_t), DIMENSION(num_moments_f) :: upstream_flux
    REAL(kind=d_t), DIMENSION(num_moments_f,num_moments_f), &
          INTENT(in) :: Lf, Uf
    REAL(kind=d_t), DIMENSION(:,:), ALLOCATABLE :: T

    ! Allocate transformation matrix

    ALLOCATE(T(num_moments_f,num_moments_f),stat=alloc_stat);&
          T=0.0_d_t
    IF(alloc_stat /= 0) CALL stop_thor(2_li)

    ! Create transformation matrix

    DO q=1, num_moments_f
      i1=index_f(q)%i1
      i2=index_f(q)%i2
      DO l=1, num_moments_f
        iup1=index_f(l)%i1
        iup2=index_f(l)%i2
        DO m11=0, iup1
          DO m12=0, iup1-m11
            DO m21=0, iup2
              DO m22=0, iup2-m21
                T(q,l)=T(q,l)+2.0_d_t*(factorial_d_t(iup1)*&
                      factorial_d_t(iup2)*af(1)**(iup1-m11-m12)*&
                      af(2)**(iup2-m21-m22)*bf(1,1)**(m11)*&
                      bf(1,2)**(m12)*bf(2,1)**(m21)*bf(2,2)**(m22))&
                      /(factorial_d_t(m11)*factorial_d_t(m12)*&
                      factorial_d_t(m21)*factorial_d_t(m22)*&
                      factorial_d_t(iup1-m11-m12)*&
                      factorial_d_t(iup2-m21-m22)*&
                      REAL((i1+m11+m21+i2+m12+m22+2.0_d_t)*&
                      (i2+m12+m22+1.0_d_t)))
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

    ! Solve for the upstream expansion coefficients

    CALL back_substitution(num_moments_f,upstream_moments,Lf,Uf,&
          upstream_flux)

    ! Solve for the downstream expansion coefficients

    transformed_moments=MATMUL(T,upstream_flux)

    CALL back_substitution(num_moments_f,transformed_moments,Lf,Uf,&
          transformed_flux)

    DEALLOCATE(T)

  END SUBROUTINE transform_incoming_moments

  SUBROUTINE transform_boundary_moments(Lf,Uf,transformed_moments,transformed_flux)
    !*********************************************************************
    !
    ! Subroutine transform incoming moments transforms the upstream cell
    ! face outgoing moments into downstream incoming cell face moments
    !
    !*********************************************************************

    ! Define temporary variables

    REAL(kind=d_t), DIMENSION(num_moments_f), INTENT(in) :: &
          transformed_moments
    REAL(kind=d_t), DIMENSION(num_moments_f), INTENT(out) :: &
          transformed_flux
    REAL(kind=d_t), DIMENSION(num_moments_f,num_moments_f), &
          INTENT(in) :: Lf, Uf

    CALL back_substitution(num_moments_f,transformed_moments,Lf,Uf,&
          transformed_flux)

  END SUBROUTINE transform_boundary_moments

  SUBROUTINE incoming_cell_subcell_project(Lf,Uf,af,bf,transformed_flux,incoming_flux)
    !*********************************************************************
    !
    ! Subroutine incoming cell to subcell moments projects the incoming
    ! cell face moments into incoming subcell face moments
    !
    !*********************************************************************

    ! Define variables

    INTEGER(kind=li) :: alloc_stat, l, q, i1, i2, m11, m12, m21, m22
    REAL(kind=d_t), DIMENSION(num_moments_f), INTENT(in) :: &
          transformed_flux
    REAL(kind=d_t), DIMENSION(num_moments_f), INTENT(out) :: &
          incoming_flux
    REAL(kind=d_t), DIMENSION(2), INTENT(in) :: af
    REAL(kind=d_t), DIMENSION(2,2), INTENT(in) :: bf
    REAL(kind=d_t), DIMENSION(num_moments_f) :: incoming_moments
    REAL(kind=d_t), DIMENSION(num_moments_f,num_moments_f) :: Lf, Uf
    REAL(kind=d_t), DIMENSION(:,:), ALLOCATABLE :: p

    ! Allocate projection matrix

    ALLOCATE(p(num_moments_f,num_moments_f),stat=alloc_stat);&
          p=0.0_d_t
    IF(alloc_stat /= 0) CALL stop_thor(2_li)

    ! Create projection matrix

    DO q=1, num_moments_f
      DO l=1, num_moments_f
        i1=index_f(l)%i1
        i2=index_f(l)%i2
        DO m11=0, i1
          DO m12=0, i1-m11
            DO m21=0, i2
              DO m22=0, i2-m21
                p(q,l)=p(q,l)+2.0_d_t*(factorial_d_t(i1)*&
                      factorial_d_t(i2)*af(1)**(i1-m11-m12)*&
                      af(2)**(i2-m21-m22)*bf(1,1)**(m11)*&
                      bf(1,2)**(m12)*bf(2,1)**(m21)*bf(2,2)**(m22))&
                      /(factorial_d_t(m11)*factorial_d_t(m12)*&
                      factorial_d_t(m21)*factorial_d_t(m22)*&
                      factorial_d_t(i1-m11-m12)*&
                      factorial_d_t(i2-m21-m22)*&
                      REAL((index_f(q)%i1+m11+m21+&
                      index_f(q)%i2+m12+m22+2.0_d_t)*&
                      (index_f(q)%i2+m12+m22+1.0_d_t)))
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

    incoming_moments=MATMUL(p,transformed_flux)

    CALL back_substitution(num_moments_f,incoming_moments,Lf,Uf,&
          incoming_flux)

    DEALLOCATE(p)

  END SUBROUTINE incoming_cell_subcell_project

  SUBROUTINE cell_source_expansion(q_moments,LL,U,q_expansion)
    !*********************************************************************
    !
    ! Subroutine cell source expansion solves for the source expansion
    ! coefficients based on the distributed source moments
    !
    !*********************************************************************

    ! Pass source derived type

    REAL(kind=d_t), DIMENSION(num_moments_v), INTENT(in) :: q_moments

    ! Define temporary variables

    INTEGER(kind=li) ::  l
    REAL(kind=d_t), DIMENSION(num_moments_v) :: y, x
    REAL(kind=d_t), DIMENSION(num_moments_v), INTENT(inout) :: &
          q_expansion
    REAL(kind=d_t), DIMENSION(num_moments_v,num_moments_v),&
          INTENT(in) :: LL, U

    DO l=1, num_moments_v
      y(l)=q_moments(l)
    END DO

    CALL back_substitution(num_moments_v,y,LL,U,x)

    DO l=1, num_moments_v
      q_expansion(l)=x(l)
    END DO

  END SUBROUTINE cell_source_expansion

  SUBROUTINE source_projection(LL,U,a,b,cell_source,subcell_source)
    !*********************************************************************
    !
    ! Subroutine source projection projects cell source moments into
    ! subcell source moments
    !
    !*********************************************************************

    ! Define variables

    INTEGER(kind=li) ::  alloc_stat, l, q, i1, i2, i3, m11, m12, m13, &
          m21, m22, m23, m31, m32, m33
    REAL(kind=d_t) :: fact, a1temp, a2temp, a3temp
    REAL(kind=d_t), DIMENSION(3), INTENT(in) :: a
    REAL(kind=d_t), DIMENSION(3,3), INTENT(in) :: b
    REAL(kind=d_t), DIMENSION(num_moments_v) :: subcell_moments
    REAL(kind=d_t), DIMENSION(num_moments_v), INTENT(in) :: &
          cell_source
    REAL(kind=d_t), DIMENSION(num_moments_v), INTENT(out) :: &
          subcell_source
    REAL(kind=d_t), DIMENSION(num_moments_v,num_moments_v),&
          INTENT(in) :: LL, U
    REAL(kind=d_t), DIMENSION(:,:), ALLOCATABLE :: P

    ! Allocate projection matrix

    ALLOCATE(P(num_moments_v,num_moments_v),stat=alloc_stat);&
          P=0.0_d_t;
    IF(alloc_stat /= 0) CALL stop_thor(2_li)

    ! Create projection matrix

    DO q=1, num_moments_v
      DO l=1, num_moments_v
        i1=index_v(l)%i1
        i2=index_v(l)%i2
        i3=index_v(l)%i3
        fact=factorial_d_t(i1)*factorial_d_t(i2)*factorial_d_t(i3)
        DO m11=0, i1
          DO m12=0, i1-m11
            DO m13=0, i1-m11-m12
              a1temp=a(1)**(i1-m11-m12-m13)
              DO m21=0, i2
                DO m22=0, i2-m21
                  DO m23=0, i2-m21-m22
                    a2temp=a(2)**(i2-m21-m22-m23)
                    DO m31=0, i3
                      DO m32=0, i3-m31
                        DO m33=0, i3-m31-m32
                          a3temp=a(3)**(i3-m31-m32-m33)
                          P(q,l)=P(q,l)+&
                                6.0_d_t*(fact*&
                                a1temp*a2temp*a3temp*&
                                b(1,1)**(m11)*b(1,2)**(m12)*&
                                b(1,3)**(m13)*b(2,1)**(m21)*&
                                b(2,2)**(m22)*b(2,3)**(m23)*&
                                b(3,1)**(m31)*b(3,2)**(m32)*&
                                b(3,3)**(m33))/&
                                (factorial_d_t(m11)*&
                                factorial_d_t(m12)*&
                                factorial_d_t(m13)*&
                                factorial_d_t(m21)*&
                                factorial_d_t(m22)*&
                                factorial_d_t(m23)*&
                                factorial_d_t(m31)*&
                                factorial_d_t(m32)*&
                                factorial_d_t(m33)*&
                                factorial_d_t(i1-m11-m12-m13)*&
                                factorial_d_t(i2-m21-m22-m23)*&
                                factorial_d_t(i3-m31-m32-m33)*&
                                REAL(&
                                (index_v(q)%i1+m11+m21+m31+&
                                index_v(q)%i2+m12+m22+m32+&
                                index_v(q)%i3+m13+m23+m33+&
                                3.0_d_t)*&
                                (index_v(q)%i2+m12+m22+m32+&
                                index_v(q)%i3+m13+m23+m33+&
                                2.0_d_t)*&
                                (index_v(q)%i3+m13+m23+m33+&
                                1.0_d_t)))
                        END DO
                      END DO
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

    subcell_moments=MATMUL(P,cell_source)

    CALL back_substitution(num_moments_v,subcell_moments,LL,U,&
          subcell_source)

    DEALLOCATE(P)

  END SUBROUTINE source_projection

  SUBROUTINE characteristic_solver(i,t,e,a,b,af,bf,incoming_flux,subcell_source,&
        outgoing_moments,subcell_flux)
    !*********************************************************************
    !
    ! Subroutine characteristic solver computes the outgoing face
    ! angular flux based on incoming and source subcell moments
    !
    !*********************************************************************

    ! Define variables

    INTEGER(kind=li), INTENT(in) :: i
    INTEGER(kind=li) ::  alloc_stat, l, q, i1, i2, i3, m11, m12, m13, &
          m21, m22, m23, m31, m32, m33, g1, g2, g3, g3p
    REAL(kind=d_t) :: e1, e2, e3, e4, e5, e6, fact, af1temp, af2temp, &
          a1temp, a2temp, a3temp
    REAL(kind=d_t), INTENT(in) :: t, e
    REAL(kind=d_t), DIMENSION(2), INTENT(in) :: af
    REAL(kind=d_t), DIMENSION(2,2), INTENT(in) :: bf
    REAL(kind=d_t), DIMENSION(3), INTENT(in) :: a
    REAL(kind=d_t), DIMENSION(3,3), INTENT(in) :: b
    REAL(kind=d_t), DIMENSION(num_moments_f), INTENT(in) :: &
          incoming_flux
    REAL(kind=d_t), DIMENSION(num_moments_v), INTENT(in) :: &
          subcell_source
    REAL(kind=d_t), DIMENSION(num_moments_v), INTENT(out) :: &
          subcell_flux
    REAL(kind=d_t), DIMENSION(num_moments_f), INTENT(out) :: &
          outgoing_moments
    REAL(kind=d_t), DIMENSION(:,:), ALLOCATABLE :: FF, FV, F, V

    ! Allocate face and volume matrix

    ALLOCATE(FF(num_moments_f,num_moments_f),&
          FV(num_moments_f,num_moments_v),F(num_moments_v,num_moments_f),&
          V(num_moments_v,num_moments_v),stat=alloc_stat);&
          FF=0.0_d_t;FV=0.0_d_t;F=0.0_d_t;V=0.0_d_t
    IF(alloc_stat /= 0) CALL stop_thor(2_li)

    ! Pre-compute optimisers

    e1=e
    e2=e*e1
    e3=e*e2
    e4=e*e3
    e5=e*e4
    e6=e*e5

    ! Create face integral matrix (outgoing face moments)

    DO q=1, num_moments_f
      i1=index_f(q)%i1
      i2=index_f(q)%i2
      fact=factorial_d_t(i1)*factorial_d_t(i2)
      DO l=1, num_moments_f
        DO m11=0, i1
          DO m12=0, i1-m11
            af1temp=af(1)**(i1-m11-m12)
            DO m21=0, i2
              DO m22=0, i2-m21
                af2temp=af(2)**(i2-m21-m22)
                g1=index_f(l)%i1+m11+m21
                g2=index_f(l)%i2+m12+m22
                FF(q,l)=FF(q,l)+((fact*af1temp*af2temp*&
                      bf(1,1)**(m11)*bf(1,2)**(m12)*&
                      bf(2,1)**(m21)*bf(2,2)**(m22))/(&
                      factorial_d_t(m11)*factorial_d_t(m12)*&
                      factorial_d_t(m21)*factorial_d_t(m22)*&
                      factorial_d_t(i1-m11-m12)*&
                      factorial_d_t(i2-m21-m22)))*&
                      face_moment1(g1,g2,e)
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

    ! Create volume integral matrix (outgoing face moments)

    DO q=1, num_moments_f
      i1=index_f(q)%i1
      i2=index_f(q)%i2
      fact=factorial_d_t(i1)*factorial_d_t(i2)
      DO l=1, num_moments_v
        DO m11=0, i1
          DO m12=0, i1-m11
            af1temp=af(1)**(i1-m11-m12)
            DO m21=0, i2
              DO m22=0, i2-m21
                af2temp=af(2)**(i2-m21-m22)
                g1=index_v(l)%i1+m11+m21
                g2=index_v(l)%i2+m12+m22
                g3=index_v(l)%i3
                FV(q,l)=FV(q,l)+((fact*af1temp*af2temp*&
                      bf(1,1)**(m11)*bf(1,2)**(m12)*&
                      bf(2,1)**(m21)*bf(2,2)**(m22))/(&
                      factorial_d_t(m11)*factorial_d_t(m12)*&
                      factorial_d_t(m21)*factorial_d_t(m22)*&
                      factorial_d_t(i1-m11-m12)*&
                      factorial_d_t(i2-m21-m22)))*&
                      face_moment2(g1,g2,g3,t,e)
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

    ! Create face integral matrix (subcell flux moments)

    DO q=1, num_moments_v
      i1=index_v(q)%i1
      i2=index_v(q)%i2
      i3=index_v(q)%i3
      fact=factorial_d_t(i1)*factorial_d_t(i2)*factorial_d_t(i3)
      DO l=1, num_moments_f
        DO m11=0, i1
          DO m12=0, i1-m11
            DO m13=0, i1-m11-m12
              a1temp=a(1)**(i1-m11-m12-m13)
              DO m21=0, i2
                DO m22=0, i2-m21
                  DO m23=0, i2-m21-m22
                    a2temp=a(2)**(i2-m21-m22-m23)
                    DO m31=0, i3
                      DO m32=0, i3-m31
                        DO m33=0, i3-m31-m32
                          a3temp=a(3)**(i3-m31-m32-m33)
                          g1=index_f(l)%i1+m11+m21+m31
                          g2=index_f(l)%i2+m12+m22+m32
                          g3=m13+m23+m33
                          F(q,l)=F(q,l)+&
                                ((fact*a1temp*a2temp*a3temp*&
                                b(1,1)**(m11)*b(1,2)**(m12)*&
                                b(1,3)**(m13)*b(2,1)**(m21)*&
                                b(2,2)**(m22)*b(2,3)**(m23)*&
                                b(3,1)**(m31)*b(3,2)**(m32)*&
                                b(3,3)**(m33))/&
                                (factorial_d_t(m11)*&
                                factorial_d_t(m12)*&
                                factorial_d_t(m13)*&
                                factorial_d_t(m21)*&
                                factorial_d_t(m22)*&
                                factorial_d_t(m23)*&
                                factorial_d_t(m31)*&
                                factorial_d_t(m32)*&
                                factorial_d_t(m33)*&
                                factorial_d_t(i1-m11-m12-m13)*&
                                factorial_d_t(i2-m21-m22-m23)*&
                                factorial_d_t(i3-m31-m32-m33)))*&
                                volume_moment1(g1,g2,g3,e)
                        END DO
                      END DO
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

    ! Create volume integral matrix (subcell flux moments)

    DO q=1, num_moments_v
      i1=index_v(q)%i1
      i2=index_v(q)%i2
      i3=index_v(q)%i3
      fact=factorial_d_t(i1)*factorial_d_t(i2)*factorial_d_t(i3)
      DO l=1, num_moments_v
        DO m11=0, i1
          DO m12=0, i1-m11
            DO m13=0, i1-m11-m12
              a1temp=a(1)**(i1-m11-m12-m13)
              DO m21=0, i2
                DO m22=0, i2-m21
                  DO m23=0, i2-m21-m22
                    a2temp=a(2)**(i2-m21-m22-m23)
                    DO m31=0, i3
                      DO m32=0, i3-m31
                        DO m33=0, i3-m31-m32
                          a3temp=a(3)**(i3-m31-m32-m33)
                          g1=index_v(l)%i1+m11+m21+m31
                          g2=index_v(l)%i2+m12+m22+m32
                          g3=m13+m23+m33
                          g3p=index_v(l)%i3
                          V(q,l)=V(q,l)+&
                                ((fact*a1temp*a2temp*a3temp*&
                                b(1,1)**(m11)*b(1,2)**(m12)*&
                                b(1,3)**(m13)*b(2,1)**(m21)*&
                                b(2,2)**(m22)*b(2,3)**(m23)*&
                                b(3,1)**(m31)*b(3,2)**(m32)*&
                                b(3,3)**(m33))/&
                                (factorial_d_t(m11)*&
                                factorial_d_t(m12)*&
                                factorial_d_t(m13)*&
                                factorial_d_t(m21)*&
                                factorial_d_t(m22)*&
                                factorial_d_t(m23)*&
                                factorial_d_t(m31)*&
                                factorial_d_t(m32)*&
                                factorial_d_t(m33)*&
                                factorial_d_t(i1-m11-m12-m13)*&
                                factorial_d_t(i2-m21-m22-m23)*&
                                factorial_d_t(i3-m31-m32-&
                                m33)))*&
                                volume_moment2(g1,g2,g3,g3p,t,e)
                        END DO
                      END DO
                    END DO
                  END DO
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

    outgoing_moments=MATMUL(FF,incoming_flux)+MATMUL(FV,subcell_source)

    subcell_flux=MATMUL(F,incoming_flux)+MATMUL(V,subcell_source)

    DEALLOCATE(FF,FV,F,V)

  END SUBROUTINE characteristic_solver

  SUBROUTINE face_moment_transformation(Lf,Uf,f,face_angular_mom,face_cell_temp)
    !*********************************************************************
    !
    ! Subroutine face moment transformation transforms the face angular
    ! flux moments from their face coordinate systems into the cell
    ! volume coordinate system through matrix multiplication
    !
    !*********************************************************************

    ! Define variables

    INTEGER(kind=li), INTENT(in) :: f
    INTEGER(kind=li) :: alloc_stat, l, q, i1F, i2F
    REAL(kind=d_t), DIMENSION(num_moments_f), INTENT(in) :: &
          face_angular_mom
    REAL(kind=d_t), DIMENSION(num_moments_v), INTENT(out) :: &
          face_cell_temp
    REAL(kind=d_t), DIMENSION(num_moments_f) :: x
    REAL(kind=d_t), DIMENSION(num_moments_f,num_moments_f) :: Lf, Uf
    REAL(kind=d_t), DIMENSION(:,:), ALLOCATABLE :: TF

    ! Allocate projection matrix

    ALLOCATE(TF(num_moments_v,num_moments_f),stat=alloc_stat);&
          TF=0.0_d_t
    IF(alloc_stat /= 0) CALL stop_thor(2_li)

    ! Create projection matrix

    DO q=1, num_moments_v
      IF(f == 0)THEN
        i1F=index_v(q)%i2
        i2F=index_v(q)%i3
      ELSEIF(f == 1)THEN
        i1F=index_v(q)%i1+index_v(q)%i2
        i2F=index_v(q)%i3
      ELSEIF(f == 2)THEN
        i1F=index_v(q)%i1
        i2F=index_v(q)%i2+index_v(q)%i3
      ELSEIF(f == 3)THEN
        i1F=index_v(q)%i1
        i2F=index_v(q)%i2
      ELSE
        CALL stop_thor(7_li)
      ENDIF
      DO l=1, num_moments_f
        TF(q,l)=2.0_d_t/REAL((index_f(l)%i1+i1F+index_f(l)%i2+&
              i2F+2.0_d_t)*(index_f(l)%i2+i2F+1.0_d_t))
        IF(f == 3 .AND. index_v(q)%i3 /= 0)THEN
          TF(q,l)=0.0_d_t
        END IF
      END DO
    END DO

    ! Solve for expansion coefficient on face

    CALL back_substitution(num_moments_f,face_angular_mom,Lf,Uf,x)

    face_cell_temp=MATMUL(TF,x)

    DEALLOCATE(TF)

  END SUBROUTINE face_moment_transformation

  FUNCTION face_moment1(gam1,gam2,e)
    !*********************************************************************
    !
    ! Subroutine calculates face moments of characteristic integral
    ! using a Taylor expansion
    !
    !*********************************************************************
    USE types
    IMPLICIT NONE

    ! Pass input parameters
    INTEGER(kind=li), INTENT(in) :: gam1, gam2
    REAL(kind=d_t), INTENT(in) :: e

    ! Define temporary variables
    INTEGER(kind=li) :: n, max_n, ggam
    REAL(kind=d_t) :: r, new_term, face_moment1, rerror
    REAL(kind=d_t), PARAMETER :: eps=1.0e-8

    max_n=1000
    ggam=gam1+gam2
    face_moment1=0.0_d_t
    rerror=1.0_d_t

    DO n=0, max_n
      IF(n == 0)THEN
        r=2.0_d_t
      ELSE
        r=r*(-e)/n
      END IF
      new_term=r/((ggam+n+2_li)*(gam2+n+1_li))
      face_moment1=face_moment1+new_term
      rerror=ABS(new_term/face_moment1)
      IF(rerror > eps)THEN
      ELSE
        go to 10
      END IF
    END DO

10 END FUNCTION face_moment1

  FUNCTION face_moment2(gam1,gam2,gam3,t,e)
    !*********************************************************************
    !
    ! Subroutine calculates volume moments of characteristic integral
    ! using a Taylor expansion
    !
    !*********************************************************************
    USE types
    IMPLICIT NONE
    INTEGER(kind=li), INTENT(in) :: gam1, gam2, gam3
    REAL(kind=d_t), INTENT(in) :: t, e
    INTEGER(kind=li) :: n, max_n, m, max_m, ggam3, ggam2, ggam1
    REAL(kind=d_t) :: r1, r2, new_term1, new_term2, face_moment2, &
          rerror1, rerror2
    REAL(kind=d_t), PARAMETER :: eps=1.0e-8

    max_n=1000
    max_m=1000
    ggam3=gam3+gam2+gam1
    ggam2=gam3+gam2
    ggam1=gam3

    r1=0.0_d_t
    r2=0.0_d_t
    new_term1=0.0_d_t
    new_term2=0.0_d_t
    face_moment2=0.0_d_t
    rerror1=1.0_d_t
    rerror2=1.0_d_t

    DO n=0, max_n
      IF(n == 0)THEN
        r1=1.0_d_t
      ELSE
        r1=r1*e/n
      END IF
      DO m=0, max_m
        IF(m == 0)THEN
          r2=2.0_d_t*t
        ELSE
          r2=r2*(-e)/m
        END IF
        new_term1=r1*r2/((ggam3+n+m+3_li)*(ggam2+n+m+2_li)*&
              (ggam1+n+1_li))
        new_term2=new_term2+new_term1
        IF(new_term2 /= 0)THEN
          rerror1=ABS(new_term1/new_term2)
        END IF
        IF(rerror1 > eps)THEN
        ELSE
          go to 10
        END IF
      END DO

10    CONTINUE

      face_moment2=face_moment2+new_term2
      rerror2=ABS(new_term2/face_moment2)
      IF(rerror2 > eps)THEN
        new_term2=0.0_d_t
      ELSE
        go to 11
      END IF
    END DO

11 END FUNCTION face_moment2

  FUNCTION volume_moment1(gam1,gam2,gam3,e)
    !*********************************************************************
    !
    ! Subroutine calculates face moments of characteristic integral
    ! using a Taylor expansion
    !
    !*********************************************************************
    USE types
    IMPLICIT NONE

    ! Pass input parameters
    INTEGER(kind=li), INTENT(in) :: gam1, gam2, gam3
    REAL(kind=d_t), INTENT(in) :: e

    ! Define temporary variables
    INTEGER(kind=li) :: n, max_n, ggam3, ggam2
    REAL(kind=d_t) :: r, new_term, volume_moment1, rerror
    REAL(kind=d_t), PARAMETER :: eps=1.0e-8

    max_n=1000
    ggam3=gam3+gam2+gam1
    ggam2=gam3+gam2
    volume_moment1=0.0_d_t
    rerror=1.0_d_t

    DO n=0, max_n
      IF(n == 0)THEN
        r=6.0_d_t
      ELSE
        r=r*(-e)/n
      END IF
      new_term=r/((ggam3+n+3_li)*(ggam2+n+2_li)*(gam3+n+1_li))
      volume_moment1=volume_moment1+new_term
      rerror=ABS(new_term/volume_moment1)
      IF(rerror > eps)THEN
      ELSE
        go to 10
      END IF
    END DO

10 END FUNCTION volume_moment1

  FUNCTION volume_moment2(gam1,gam2,gam3,gam3p,t,e)
    !*********************************************************************
    !
    ! Subroutine calculates volume moments of characteristic integral
    ! using a Taylor expansion
    !
    !*********************************************************************
    USE types
    IMPLICIT NONE
    INTEGER(kind=li), INTENT(in) :: gam1, gam2, gam3, gam3p
    REAL(kind=d_t), INTENT(in) :: t, e
    INTEGER(kind=li) :: n, max_n, m, max_m, ggam3, ggam2, ggam1
    REAL(kind=d_t) :: r1, r2, new_term1, new_term2, volume_moment2, &
          rerror1, rerror2
    REAL(kind=d_t), PARAMETER :: eps=1.0e-8

    max_n=1000
    max_m=1000
    ggam3=gam3p+gam3+gam2+gam1
    ggam2=gam3p+gam3+gam2
    ggam1=gam3p+gam3

    r1=0.0_d_t
    r2=0.0_d_t
    new_term1=0.0_d_t
    new_term2=0.0_d_t
    volume_moment2=0.0_d_t
    rerror1=1.0_d_t
    rerror2=1.0_d_t

    DO n=0, max_n
      IF(n == 0)THEN
        r1=1.0_d_t
      ELSE
        r1=r1*e/n
      END IF
      DO m=0, max_m
        IF(m == 0)THEN
          r2=6.0_d_t*t
        ELSE
          r2=r2*(-e)/m
        END IF
        new_term1=r1*r2/((ggam3+n+m+4_li)*(ggam2+n+m+3_li)*&
              (ggam1+n+m+2_li)*(gam3p+n+1_li))
        new_term2=new_term2+new_term1
        IF(new_term2 /= 0)THEN
          rerror1=ABS(new_term1/new_term2)
        END IF
        IF(rerror1 > eps)THEN
        ELSE
          go to 10
        END IF
      END DO

10    CONTINUE

      volume_moment2=volume_moment2+new_term2
      rerror2=ABS(new_term2/volume_moment2)
      IF(rerror2 > eps)THEN
        new_term2=0.0_d_t
      ELSE
        go to 11
      END IF
    END DO

11 END FUNCTION volume_moment2

END MODULE transport_kernel_module_CCE
