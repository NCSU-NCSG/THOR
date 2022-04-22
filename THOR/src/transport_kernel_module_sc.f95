MODULE transport_kernel_module_SC
  !***********************************************************************
  !
  ! Transport kernel module performs the transport calculation
  ! based on the step characteristic (LC) discretization.
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

  SUBROUTINE transport_kernel_SC(sigmat,q_moments,vol_moment,face_moment, &
        i,t,subcells,r0,r1,r2,r3,     &
        face_known,incoming_face,outgoing_face,&
        upstream_moments)
    !*********************************************************************
    !
    ! Subroutine transport cell computes outgoing angular face and
    ! volume moments based on incoming angular face fluxes and source
    ! moments
    !
    !*********************************************************************

    ! Pass source derived type

    REAL(kind=d_t),    INTENT(in)    :: q_moments(num_moments_v)

    ! Declare angular and scalar flux types used globally

    REAL(kind=d_t), INTENT(inout)    :: vol_moment(num_moments_v)
    REAL(kind=d_t), INTENT(inout)    :: face_moment(num_moments_v,0:3)

    ! Define temporary variables

    INTEGER(kind=li)                 :: ii
    INTEGER(kind=li)                 :: f
    INTEGER(kind=li),  INTENT(in)    :: i
    INTEGER(kind=li),  INTENT(in)    :: subcells
    INTEGER(kind=li),  INTENT(in)    :: incoming_face(subcells)
    INTEGER(kind=li),  INTENT(in)    :: outgoing_face(subcells)
    REAL(kind=d_t)                   :: det_Js
    REAL(kind=d_t)                   :: e
    REAL(kind=d_t)                   :: Js(3,3)
    REAL(kind=d_t)                   :: Js_inv(3,3)
    REAL(kind=d_t)                   :: q_expansion
    REAL(kind=d_t)                   :: area(subcells)
    REAL(kind=d_t)                   :: volume(subcells)
    REAL(kind=d_t)                   :: cell_source
    REAL(kind=d_t)                   :: subcell_source
    REAL(kind=d_t)                   :: subcell_flux
    REAL(kind=d_t)                   :: incoming_flux
    REAL(kind=d_t)                   :: transformed_flux
    REAL(kind=d_t)                   :: outgoing_moments
    REAL(kind=d_t)                   :: projected_moments(subcells)
    REAL(kind=d_t)                   :: face_angular(0:3)
    REAL(kind=d_t)                   :: projected_flux(subcells)
    REAL(kind=d_t),    INTENT(in)    :: upstream_moments(num_moments_v,subcells)
    REAL(kind=d_t),    INTENT(in)    :: sigmat
    REAL(kind=d_t),    INTENT(in)    :: t
    INTEGER(kind=1),   INTENT(inout) :: face_known(0:3,num_cells)
    TYPE(vector)                     :: v(0:3)
    TYPE(vector),      INTENT(in)    :: r0(subcells)
    TYPE(vector),      INTENT(in)    :: r1(subcells)
    TYPE(vector),      INTENT(in)    :: r2(subcells)
    TYPE(vector),      INTENT(in)    :: r3(subcells)

    ! Compute optical thickness

    e=sigmat*t

    ! Loop over subcells and transport through each

    DO ii=1, subcells

      ! Transform upstream cell outgoing face moments into incoming face moments

      transformed_flux=upstream_moments(1,ii)
      face_moment(1,incoming_face(ii))=transformed_flux

      ! Transform incoming cell face moments into subcell face moments

      v(0)=r0(ii)
      v(1)=r1(ii)
      v(2)=r2(ii)
      v(3)=r3(ii)

      CALL subcell_jacobian(v,Js,det_Js,Js_inv)

      volume(ii)=(1.0_d_t/6.0_d_t)*det_Js

      incoming_flux=transformed_flux

      ! Compute cell source expansion coefficients based on source moments

      q_expansion=q_moments(1)

      ! Project cell source moment into subcell moments in subcell system

      cell_source=q_expansion

      subcell_source=cell_source

      ! Compute outgoing face moments in subcell with characteristic relation

      CALL characteristic_solver(t,e,incoming_flux,subcell_source,&
            outgoing_moments,subcell_flux)

      projected_moments(ii)=outgoing_moments

      ! Project subcell flux moments into cell flux moments

      projected_flux(ii)=subcell_flux

    END DO

    ! Compute area weighted outgoing face angular moments and update
    ! face_known array

    DO ii=1, subcells
      area(ii)=(0.5_d_t)*vmag((r1(ii)-r0(ii)) .cross. (r3(ii)-r1(ii)))

      face_moment(1,outgoing_face(ii))=&
            face_moment(1,outgoing_face(ii))+(area(ii)/((0.5_d_t)*&
            vmag(outward_normal(i,outgoing_face(ii)))))*&
            projected_moments(ii)

      IF(adjacency_list(i,outgoing_face(ii))%cell /= 0)THEN
        face_known(adjacency_list(i,outgoing_face(ii))%face,&
              adjacency_list(i,outgoing_face(ii))%cell)=1
      END IF
    END DO

    ! Transform face flux into face moments

    DO f=0, 3
      face_angular(f)=face_moment(1,f)
    END DO

    ! Compute volume weighted angular flux moments over cell

    DO ii=1, subcells
      vol_moment(1)=vol_moment(1)+&
            (volume(ii)/cells(i)%volume)*projected_flux(ii)
    END DO

  END SUBROUTINE transport_kernel_SC

  SUBROUTINE characteristic_solver(t,e,incoming_flux,subcell_source,&
        outgoing_moments,subcell_flux)
    !*********************************************************************
    !
    ! Subroutine characteristic solver computes the outgoing face
    ! angular flux based on incoming and source subcell moments
    !
    !*********************************************************************

    IMPLICIT NONE

    ! Define temporary variables

    REAL(kind=d_t), INTENT(in)  :: t
    REAL(kind=d_t), INTENT(in)  :: e
    REAL(kind=d_t), INTENT(in)  :: incoming_flux
    REAL(kind=d_t), INTENT(in)  :: subcell_source
    REAL(kind=d_t), INTENT(out) :: subcell_flux
    REAL(kind=d_t), INTENT(out) :: outgoing_moments
    REAL(kind=d_t)              :: e2
    REAL(kind=d_t)              :: e3
    REAL(kind=d_t)              :: e4
    REAL(kind=d_t)              :: e5
    REAL(kind=d_t)              :: e6
    REAL(kind=d_t)              :: FF
    REAL(kind=d_t)              :: FV
    REAL(kind=d_t)              :: F
    REAL(kind=d_t)              :: V

    e2=e*e
    e3=e*e2
    e4=e*e3
    e5=e*e4
    e6=e*e5

    IF(e > 0.25_d_t)THEN
      FF=2.0/e-2.0/e2+2.0*EXP(-e)/e2
      FV=1.0/e-2.0/e2+2.0/e3-2.0*EXP(-e)/e3
      F=3.0/e-6.0/e2+6.0/e3-6.0*EXP(-e)/e3
      V=1.0/e-3.0/e2+6.0/e3-6.0/e4+6.0*EXP(-e)/e4
    ELSE
      FF=1.0-0.333333*e+0.0833333*e2-0.0166667*e3+&
            0.00277778*e4-0.000396825*e5+0.0000496032*e6
      FV=0.333333-0.0833333*e+0.0166667*e2-0.00277778*e3+&
            0.000396825*e4-0.0000496032*e5+(5.51146e-6)*e6
      F=1.0-0.25*e+0.05*e2-0.00833333*e3+0.00119048*e4-&
            0.00014881*e5+0.0000165344*e6
      V=0.25-0.05*e+0.00833333*e2-0.00119048*e3+0.00014881*e4-&
            0.0000165344*e5+(1.65344e-6)*e6
    END IF

    outgoing_moments=FF*incoming_flux+t*FV*subcell_source

    subcell_flux=F*incoming_flux+t*V*subcell_source

  END SUBROUTINE characteristic_solver

END MODULE transport_kernel_module_SC
