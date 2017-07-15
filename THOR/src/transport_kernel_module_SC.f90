module transport_kernel_module_SC
!***********************************************************************
!
! Transport kernel module performs the transport calculation
! based on the step characteristic (LC) discretization.
!
!***********************************************************************

! User derived-type modules

  use types
  use parameter_types
  use filename_types
  use vector_types
  use cross_section_types
  use geometry_types
  use angle_types
  use multindex_types
  use global_variables
  use termination_module
  use general_utility_module

  implicit none

contains

  subroutine transport_kernel_SC(sigmat,q_moments,vol_moment,face_moment, &
                                 LL,U,Lf,Uf,i,t,subcells,r0,r1,r2,r3,     &
                                 face_known,incoming_face,outgoing_face,J,&
                                 J_inv,upstream_moments)
  !*********************************************************************
  !
  ! Subroutine transport cell computes outgoing angular face and 
  ! volume moments based on incoming angular face fluxes and source 
  ! moments
  !
  !*********************************************************************

  ! Pass source derived type

    real(kind=d_t),    intent(in)    :: q_moments(num_moments_v)

  ! Declare angular and scalar flux types used globally

    real(kind=d_t), intent(inout)    :: vol_moment(num_moments_v)
    real(kind=d_t), intent(inout)    :: face_moment(num_moments_v,0:3) 

  ! Declare the matrices LL,U,Lf,Uf

    real(kind=d_t), dimension(num_moments_v,num_moments_v) :: LL, U
    real(kind=d_t), dimension(num_moments_f,num_moments_f) :: Lf, Uf

  ! Define temporary variables

    integer(kind=li)                 :: ii
    integer(kind=li)                 :: l
    integer(kind=li)                 :: f
    integer(kind=li),  intent(in)    :: i
    integer(kind=li),  intent(in)    :: subcells 
    integer(kind=li),  intent(in)    :: incoming_face(subcells)
    integer(kind=li),  intent(in)    :: outgoing_face(subcells)
    real(kind=d_t)                   :: det_Jup
    real(kind=d_t)                   :: det_Js
    real(kind=d_t)                   :: e
    real(kind=d_t)                   :: af(2)
    real(kind=d_t)                   :: a_temp(3)
    real(kind=d_t)                   :: a(3)
    real(kind=d_t)                   :: bf(2,2)
    real(kind=d_t)                   :: JFup(3,2)
    real(kind=d_t)                   :: JFdown(3,2)
    real(kind=d_t)                   :: JF(3,2)
    real(kind=d_t)                   :: Jsf(3,2)
    real(kind=d_t)                   :: JFup_inv(2,3)
    real(kind=d_t)                   :: JF_inv(2,3)
    real(kind=d_t)                   :: Jup(3,3)
    real(kind=d_t)                   :: Jup_inv(3,3)
    real(kind=d_t)                   :: Js(3,3)
    real(kind=d_t)                   :: Js_inv(3,3)
    real(kind=d_t)                   :: b(3,3)
    real(kind=d_t),    intent(in)    :: J(3,3)
    real(kind=d_t),    intent(in)    :: J_inv(3,3)
    real(kind=d_t)                   :: q_expansion
    real(kind=d_t)                   :: area(subcells)
    real(kind=d_t)                   :: volume(subcells)
    real(kind=d_t)                   :: cell_source
    real(kind=d_t)                   :: subcell_source
    real(kind=d_t)                   :: subcell_flux
    real(kind=d_t)                   :: proj_flux_moments
    real(kind=d_t)                   :: face_cell_mom
    real(kind=d_t)                   :: face_cell_temp
    real(kind=d_t)                   :: y
    real(kind=d_t)                   :: subcell_upstream_moments
    real(kind=d_t)                   :: incoming_flux
    real(kind=d_t)                   :: outgoing_flux
    real(kind=d_t)                   :: transformed_moments
    real(kind=d_t)                   :: transformed_flux
    real(kind=d_t)                   :: outgoing_moments
    real(kind=d_t)                   :: projected_temp
    real(kind=d_t)                   :: proj_moments
    real(kind=d_t)                   :: face_angular_mom
    real(kind=d_t)                   :: subcell_source_moments(subcells)
    real(kind=d_t)                   :: projected_moments(subcells)
    real(kind=d_t)                   :: face_angular(0:3)
    real(kind=d_t)                   :: projected_flux(subcells)
    real(kind=d_t),    intent(in)    :: upstream_moments(num_moments_v,subcells) 
    real(kind=d_t),    intent(in)    :: sigmat
    real(kind=d_t),    intent(in)    :: t
    integer(kind=1),   intent(inout) :: face_known(0:3,num_cells)
    type(vector)                     :: R0down 
    type(vector)                     :: R0up 
    type(vector)                     :: R0F 
    type(vector)                     :: v0
    type(vector)                     :: v1
    type(vector)                     :: v2
    type(vector)                     :: v3
    type(vector)                     :: v(0:3)
    type(vector),      intent(in)    :: r0(subcells)
    type(vector),      intent(in)    :: r1(subcells)
    type(vector),      intent(in)    :: r2(subcells)
    type(vector),      intent(in)    :: r3(subcells)

  ! Compute optical thickness

    e=sigmat*t

  ! Loop over subcells and transport through each

    do ii=1, subcells

  ! Transform upstream cell outgoing face moments into incoming face moments

       transformed_flux=upstream_moments(1,ii)
       face_moment(1,incoming_face(ii))=transformed_flux

  ! Transform incoming cell face moments into subcell face moments

       v(0)=r0(ii)
       v(1)=r1(ii)
       v(2)=r2(ii)
       v(3)=r3(ii)
       
       call subcell_jacobian(v,Js,det_Js,Js_inv)
       
       volume(ii)=(1.0_d_t/6.0_d_t)*det_Js

       incoming_flux=transformed_flux

  ! Compute cell source expansion coefficients based on source moments

       q_expansion=q_moments(1)

  ! Project cell source moment into subcell moments in subcell system

       cell_source=q_expansion

       subcell_source=cell_source

  ! Compute outgoing face moments in subcell with characteristic relation

       call characteristic_solver(t,e,incoming_flux,subcell_source,&
            outgoing_moments,subcell_flux)

       projected_moments(ii)=outgoing_moments
 
  ! Project subcell flux moments into cell flux moments

       projected_flux(ii)=subcell_flux

    end do

  ! Compute area weighted outgoing face angular moments and update 
  ! face_known array

    do ii=1, subcells
       area(ii)=(0.5_d_t)*vmag((r1(ii)-r0(ii)) .cross. (r3(ii)-r1(ii)))

       face_moment(1,outgoing_face(ii))=&
            face_moment(1,outgoing_face(ii))+(area(ii)/((0.5_d_t)*&
            vmag(outward_normal(i,outgoing_face(ii)))))*&
            projected_moments(ii)

       if(adjacency_list(i,outgoing_face(ii))%cell /= 0)then
          face_known(adjacency_list(i,outgoing_face(ii))%face,&
               adjacency_list(i,outgoing_face(ii))%cell)=1
       end if
    end do

  ! Transform face flux into face moments

    do f=0, 3
       face_angular(f)=face_moment(1,f)
    end do

  ! Compute volume weighted angular flux moments over cell

    do ii=1, subcells
       vol_moment(1)=vol_moment(1)+&
            (volume(ii)/cells(i)%volume)*projected_flux(ii)
    end do

  end subroutine

  subroutine characteristic_solver(t,e,incoming_flux,subcell_source,&
       outgoing_moments,subcell_flux)
  !*********************************************************************
  !
  ! Subroutine characteristic solver computes the outgoing face 
  ! angular flux based on incoming and source subcell moments
  !
  !*********************************************************************

    implicit none

  ! Define temporary variables

    real(kind=d_t), intent(in)  :: t
    real(kind=d_t), intent(in)  :: e
    real(kind=d_t), intent(in)  :: incoming_flux
    real(kind=d_t), intent(in)  :: subcell_source
    real(kind=d_t), intent(out) :: subcell_flux 
    real(kind=d_t), intent(out) :: outgoing_moments
    real(kind=d_t)              :: e2
    real(kind=d_t)              :: e3
    real(kind=d_t)              :: e4
    real(kind=d_t)              :: e5
    real(kind=d_t)              :: e6
    real(kind=d_t)              :: FF
    real(kind=d_t)              :: FV
    real(kind=d_t)              :: F
    real(kind=d_t)              :: V

    e2=e*e
    e3=e*e2
    e4=e*e3
    e5=e*e4
    e6=e*e5

    if(e > 0.25_d_t)then
       FF=2.0/e-2.0/e2+2.0*exp(-e)/e2
       FV=1.0/e-2.0/e2+2.0/e3-2.0*exp(-e)/e3
       F=3.0/e-6.0/e2+6.0/e3-6.0*exp(-e)/e3
       V=1.0/e-3.0/e2+6.0/e3-6.0/e4+6.0*exp(-e)/e4
    else
       FF=1.0-0.333333*e+0.0833333*e2-0.0166667*e3+&
          0.00277778*e4-0.000396825*e5+0.0000496032*e6
       FV=0.333333-0.0833333*e+0.0166667*e2-0.00277778*e3+&
          0.000396825*e4-0.0000496032*e5+(5.51146e-6)*e6
       F=1.0-0.25*e+0.05*e2-0.00833333*e3+0.00119048*e4-&
         0.00014881*e5+0.0000165344*e6
       V=0.25-0.05*e+0.00833333*e2-0.00119048*e3+0.00014881*e4-& 
            0.0000165344*e5+(1.65344e-6)*e6
    end if

    outgoing_moments=FF*incoming_flux+t*FV*subcell_source

    subcell_flux=F*incoming_flux+t*V*subcell_source

  end subroutine characteristic_solver

end module transport_kernel_module_SC

