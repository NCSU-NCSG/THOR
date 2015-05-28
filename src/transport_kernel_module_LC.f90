module transport_kernel_module_LC
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
  
  subroutine transport_kernel_LC(sigmat,q_moments,vol_moment,face_moment, &
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

  ! Pass arguments

    real(kind=d_t), intent(in) :: sigmat, t
    real(kind=d_t), dimension(num_moments_v), intent(in) :: q_moments
    real(kind=d_t), dimension(num_moments_v), intent(inout) :: vol_moment
    real(kind=d_t), dimension(num_moments_f,0:3), intent(inout) :: face_moment
    real(kind=d_t), dimension(num_moments_v,num_moments_v) :: LL, U
    real(kind=d_t), dimension(num_moments_f,num_moments_f) :: Lf, Uf
    integer(kind=li), intent(in) :: i, subcells
    type(vector), dimension(subcells), intent(in) :: r0, r1, r2, r3
    integer(kind=1), dimension(0:3,num_cells), intent(inout) :: face_known
    integer(kind=li), dimension(subcells), intent(in) :: incoming_face, outgoing_face
    real(kind=d_t), dimension(3,3), intent(in) :: J, J_inv
    real(kind=d_t), dimension(num_moments_f,subcells), intent(in) :: upstream_moments
    

  ! Define local variables

    integer(kind=li) :: alloc_stat, ii, l, f, adjcnt_cell
    real(kind=d_t) :: det_Jup, det_Js, e
    real(kind=d_t), dimension(2) :: af
    real(kind=d_t), dimension(3) :: a_temp, a
    real(kind=d_t), dimension(2,2) :: bf
    real(kind=d_t), dimension(3,2) :: JFup, JFdown, JF, Jsf
    real(kind=d_t), dimension(2,3) :: JFup_inv, JF_inv
    real(kind=d_t), dimension(3,3) :: Jup, Jup_inv, Js, Js_inv, b
    real(kind=d_t), dimension(num_moments_v) :: q_expansion
    type(vector) :: R0down, R0up, R0F, v0, v1, v2, v3
    type(vector), dimension(0:3) :: v
    real(kind=d_t)  :: area(subcells)
    real(kind=d_t)  :: volume(subcells)
    real(kind=d_t)  :: subcell_upstream_moments(num_moments_f)
    real(kind=d_t)  :: incoming_flux(num_moments_f)
    real(kind=d_t)  :: transformed_moments(num_moments_f)
    real(kind=d_t)  :: transformed_flux(num_moments_f)
    real(kind=d_t)  :: cell_source(num_moments_v)
    real(kind=d_t)  :: outgoing_moments(num_moments_f)
    real(kind=d_t)  :: subcell_source(num_moments_v)
    real(kind=d_t)  :: subcell_flux(num_moments_v)
    real(kind=d_t)  :: proj_moments(num_moments_f)
    real(kind=d_t)  :: projected_moments(subcells,num_moments_f)
    real(kind=d_t)  :: projected_temp(num_moments_f)
    real(kind=d_t)  :: projected_flux(subcells,num_moments_v)
    real(kind=d_t)  :: proj_flux_moments(num_moments_v)
    real(kind=d_t)  :: face_angular_mom(num_moments_f)
    real(kind=d_t)  :: face_cell_temp(num_moments_v)
    real(kind=d_t)  :: face_angular(0:3,num_moments_v)
    real(kind=d_t)  :: subcell_source_moments(subcells,num_moments_v)
    real(kind=d_t)  :: y(num_moments_v)

  ! Initialize

    area=0.0_d_t
    volume=0.0_d_t
    subcell_upstream_moments=0.0_d_t
    incoming_flux=0.0_d_t
    transformed_moments=0.0_d_t
    transformed_flux=0.0_d_t
    cell_source=0.0_d_t
    outgoing_moments=0.0_d_t
    subcell_source=0.0_d_t
    subcell_flux=0.0_d_t
    proj_moments=0.0_d_t
    projected_moments=0.0_d_t
    projected_temp=0.0_d_t
    projected_flux=0.0_d_t
    proj_flux_moments=0.0_d_t
    face_angular_mom=0.0_d_t
    face_cell_temp=0.0_d_t
    face_angular=0.0_d_t
    subcell_source_moments=0.0_d_t
    y=0.0_d_t

  ! Compute optical thickness

    e=sigmat*t

  ! Loop over subcells and transport through each

    do ii=1, subcells

  ! Transform upstream cell outgoing face moments into incoming face moments

       if(adjacency_list(i,incoming_face(ii))%cell /= 0)then
          adjcnt_cell=adjacency_list(i,incoming_face(ii))%cell
          v0=vertices(cells(adjcnt_cell)%R(0))%v
          v1=vertices(cells(adjcnt_cell)%R(1))%v
          v2=vertices(cells(adjcnt_cell)%R(2))%v
          v3=vertices(cells(adjcnt_cell)%R(3))%v
          call cell_jacobian(v0,v1,v2,v3,Jup,det_Jup,Jup_inv)
          if(adjacency_list(i,incoming_face(ii))%face == 0)then
             JFup(1,1)=Jup(1,2)
             JFup(2,1)=Jup(2,2)
             JFup(3,1)=Jup(3,2)
             JFup(1,2)=Jup(1,3)
             JFup(2,2)=Jup(2,3)
             JFup(3,2)=Jup(3,3)
             R0up=vertices(cells(adjacency_list(i,&
                  incoming_face(ii))%cell)%R(1))%v
          elseif(adjacency_list(i,incoming_face(ii))%face == 1)then
             JFup(1,1)=Jup(1,1)+Jup(1,2)
             JFup(2,1)=Jup(2,1)+Jup(2,2)
             JFup(3,1)=Jup(3,1)+Jup(3,2)
             JFup(1,2)=Jup(1,3)
             JFup(2,2)=Jup(2,3)
             JFup(3,2)=Jup(3,3)
             R0up=vertices(cells(adjacency_list(i,&
                  incoming_face(ii))%cell)%R(0))%v
          elseif(adjacency_list(i,incoming_face(ii))%face == 2)then
             JFup(1,1)=Jup(1,1)
             JFup(2,1)=Jup(2,1)
             JFup(3,1)=Jup(3,1)
             JFup(1,2)=Jup(1,2)+Jup(1,3)
             JFup(2,2)=Jup(2,2)+Jup(2,3)
             JFup(3,2)=Jup(3,2)+Jup(3,3)
             R0up=vertices(cells(adjacency_list(i,&
                  incoming_face(ii))%cell)%R(0))%v
          elseif(adjacency_list(i,incoming_face(ii))%face == 3)then
             JFup(1,1)=Jup(1,1)
             JFup(2,1)=Jup(2,1)
             JFup(3,1)=Jup(3,1)
             JFup(1,2)=Jup(1,2)
             JFup(2,2)=Jup(2,2)
             JFup(3,2)=Jup(3,2)
             R0up=vertices(cells(adjacency_list(i,&
                  incoming_face(ii))%cell)%R(0))%v
          else
             call stop_thor(3_li) 
          end if

          if(incoming_face(ii) == 0)then
             JFdown(1,1)=J(1,2)
             JFdown(2,1)=J(2,2)
             JFdown(3,1)=J(3,2)
             JFdown(1,2)=J(1,3)
             JFdown(2,2)=J(2,3)
             JFdown(3,2)=J(3,3)
             R0down=vertices(cells(i)%R(1))%v
          elseif(incoming_face(ii) == 1)then
             JFdown(1,1)=J(1,1)+J(1,2)
             JFdown(2,1)=J(2,1)+J(2,2)
             JFdown(3,1)=J(3,1)+J(3,2)
             JFdown(1,2)=J(1,3)
             JFdown(2,2)=J(2,3)
             JFdown(3,2)=J(3,3)
             R0down=vertices(cells(i)%R(0))%v
          elseif(incoming_face(ii) == 2)then
             JFdown(1,1)=J(1,1)
             JFdown(2,1)=J(2,1)
             JFdown(3,1)=J(3,1)
             JFdown(1,2)=J(1,2)+J(1,3)
             JFdown(2,2)=J(2,2)+J(2,3)
             JFdown(3,2)=J(3,2)+J(3,3)
             R0down=vertices(cells(i)%R(0))%v
          elseif(incoming_face(ii) == 3)then
             JFdown(1,1)=J(1,1)
             JFdown(2,1)=J(2,1)
             JFdown(3,1)=J(3,1)
             JFdown(1,2)=J(1,2)
             JFdown(2,2)=J(2,2)
             JFdown(3,2)=J(3,2)
             R0down=vertices(cells(i)%R(0))%v
          else
             call stop_thor(4_li)
          end if

          call invert_face_jacobian(JFup,JFup_inv)


          a_temp(1)=R0down%x1-R0up%x1
          a_temp(2)=R0down%x2-R0up%x2
          a_temp(3)=R0down%x3-R0up%x3
          af=MATMUL(JFup_inv,a_temp)
          bf=MATMUL(JFup_inv,JFdown)

          do l=1, num_moments_f
             subcell_upstream_moments(l)=upstream_moments(l,ii)
          end do

          call transform_incoming_moments(Lf,Uf,af,bf,subcell_upstream_moments,&
                                          transformed_moments,transformed_flux )

       else
          do l=1, num_moments_f
             transformed_moments(l)=upstream_moments(l,ii)
          end do
          call transform_boundary_moments(Lf,Uf,transformed_moments,transformed_flux)
       end if

       do l=1, num_moments_f
          face_moment(l,incoming_face(ii))=transformed_moments(l)
       end do

  ! Transform incoming cell face moments into subcell face moments

       v(0)=r0(ii)

       v(1)=r1(ii)

       v(2)=r2(ii)

       v(3)=r3(ii)

       call subcell_jacobian(v,Js,det_Js,Js_inv)

       volume(ii)=(1.0_d_t/6.0_d_t)*det_Js
       
       Jsf(1,1)=Js(1,1)
       Jsf(2,1)=Js(2,1)
       Jsf(3,1)=Js(3,1)
       Jsf(1,2)=Js(1,2)
       Jsf(2,2)=Js(2,2)
       Jsf(3,2)=Js(3,2)
       
       if(incoming_face(ii) == 0)then
          JF(1,1)=J(1,2)
          JF(2,1)=J(2,2)
          JF(3,1)=J(3,2)
          JF(1,2)=J(1,3)
          JF(2,2)=J(2,3)
          JF(3,2)=J(3,3)
          R0F=vertices(cells(i)%R(1))%v
       elseif(incoming_face(ii) == 1)then
          JF(1,1)=J(1,1)+J(1,2)
          JF(2,1)=J(2,1)+J(2,2)
          JF(3,1)=J(3,1)+J(3,2)
          JF(1,2)=J(1,3)
          JF(2,2)=J(2,3)
          JF(3,2)=J(3,3)
          R0F=vertices(cells(i)%R(0))%v
       elseif(incoming_face(ii) == 2)then
          JF(1,1)=J(1,1)
          JF(2,1)=J(2,1)
          JF(3,1)=J(3,1)
          JF(1,2)=J(1,2)+J(1,3)
          JF(2,2)=J(2,2)+J(2,3)
          JF(3,2)=J(3,2)+J(3,3)
          R0F=vertices(cells(i)%R(0))%v
       elseif(incoming_face(ii) == 3)then
          JF(1,1)=J(1,1)
          JF(2,1)=J(2,1)
          JF(3,1)=J(3,1)
          JF(1,2)=J(1,2)
          JF(2,2)=J(2,2)
          JF(3,2)=J(3,2)
          R0F=vertices(cells(i)%R(0))%v
       else
          call stop_thor(5_li)
       end if

       call invert_face_jacobian(JF,JF_inv)

       a_temp(1)=r0(ii)%x1-R0F%x1
       a_temp(2)=r0(ii)%x2-R0F%x2
       a_temp(3)=r0(ii)%x3-R0F%x3
       af=MATMUL(JF_inv,a_temp)
       bf=MATMUL(JF_inv,Jsf)

       call incoming_cell_subcell_project(Lf,Uf,af,bf,transformed_flux,incoming_flux)

  ! Compute cell source expansion coefficients based on source moments

       call cell_source_expansion(q_moments,LL,U,q_expansion)
       
  ! Project cell source moment into subcell moments in subcell system

       a_temp(1)=r0(ii)%x1-vertices(cells(i)%R(0))%v%x1
       a_temp(2)=r0(ii)%x2-vertices(cells(i)%R(0))%v%x2
       a_temp(3)=r0(ii)%x3-vertices(cells(i)%R(0))%v%x3
       a=MATMUL(J_inv,a_temp)
       b=MATMUL(J_inv,Js)

       do l=1, num_moments_v
          cell_source(l)=q_expansion(l)
       end do

       call source_projection(LL,U,a,b,cell_source,subcell_source)

  ! Compute outgoing face moments in subcell with characteristic relation

  ! Project outgoing subcell moments into cell outgoing moments

       Jsf(1,1)=Js(1,1)
       Jsf(2,1)=Js(2,1)
       Jsf(3,1)=Js(3,1)
       Jsf(1,2)=Js(1,2)+Js(1,3)
       Jsf(2,2)=Js(2,2)+Js(2,3)
       Jsf(3,2)=Js(3,2)+Js(3,3)

       if(outgoing_face(ii) == 0)then
          JF(1,1)=J(1,2)
          JF(2,1)=J(2,2)
          JF(3,1)=J(3,2)
          JF(1,2)=J(1,3)
          JF(2,2)=J(2,3)
          JF(3,2)=J(3,3)
          R0F=vertices(cells(i)%R(1))%v
       elseif(outgoing_face(ii) == 1)then
          JF(1,1)=J(1,1)+J(1,2)
          JF(2,1)=J(2,1)+J(2,2)
          JF(3,1)=J(3,1)+J(3,2)
          JF(1,2)=J(1,3)
          JF(2,2)=J(2,3)
          JF(3,2)=J(3,3)
          R0F=vertices(cells(i)%R(0))%v
       elseif(outgoing_face(ii) == 2)then
          JF(1,1)=J(1,1)
          JF(2,1)=J(2,1)
          JF(3,1)=J(3,1)
          JF(1,2)=J(1,2)+J(1,3)
          JF(2,2)=J(2,2)+J(2,3)
          JF(3,2)=J(3,2)+J(3,3)
          R0F=vertices(cells(i)%R(0))%v
       elseif(outgoing_face(ii) == 3)then
          JF(1,1)=J(1,1)
          JF(2,1)=J(2,1)
          JF(3,1)=J(3,1)
          JF(1,2)=J(1,2)
          JF(2,2)=J(2,2)
          JF(3,2)=J(3,2)
          R0F=vertices(cells(i)%R(0))%v
       else
          call stop_thor(6_li)
       end if

       call invert_face_jacobian(JF,JF_inv)

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

       call characteristic_solver(i,t,e,a,b,af,bf,incoming_flux,&
                    subcell_source,outgoing_moments,subcell_flux)

       do l=1, num_moments_f
          projected_moments(ii,l)=outgoing_moments(l)
       end do
 
  ! Project subcell flux moments into cell flux moments

       do l=1, num_moments_v
          projected_flux(ii,l)= subcell_flux(l)
       end do
    end do

  ! Compute area weighted outgoing face angular moments and update 
  ! face_known array

    do ii=1, subcells
       area(ii)=(0.5_d_t)*vmag((r1(ii)-r0(ii)) .cross. (r3(ii)-r1(ii))) 
       do l=1, num_moments_f

          face_moment(l,outgoing_face(ii))=&
               face_moment(l,outgoing_face(ii))+(area(ii)/((0.5_d_t)*&
               vmag(outward_normal(i,outgoing_face(ii)))))*&
               projected_moments(ii,l)
       end do
       if(adjacency_list(i,outgoing_face(ii))%cell /= 0)then
          face_known(adjacency_list(i,outgoing_face(ii))%face,&
               adjacency_list(i,outgoing_face(ii))%cell)=1
       end if
    end do

  ! Transform face flux into face moments

    do f=0, 3
       do l=1, num_moments_f
          face_angular_mom(l)=face_moment(l,f)
       end do
       
       call face_moment_transformation(Lf,Uf,f,face_angular_mom,face_cell_temp)
       
       do l=1, num_moments_v
          face_angular(f,l)=face_cell_temp(l)
       end do
    end do

  ! Compute volume weighted angular flux moments over cell

    do l=1, num_moments_v
       do ii=1, subcells
          y(l)=y(l)+(volume(ii)/cells(i)%volume)*projected_flux(ii,l)
       end do
    end do

    do l=1, num_moments_v
       vol_moment(l)=y(l)
    end do

  end subroutine transport_kernel_LC

  subroutine transform_incoming_moments(Lf,Uf,af,bf,upstream_moments,transformed_moments,&
                                        transformed_flux)
  !*********************************************************************
  !
  ! Subroutine transform incoming moments transforms the upstream cell
  ! face outgoing moments into downstream incoming cell face moments
  !
  !*********************************************************************

  ! Define variables

    integer(kind=li) :: alloc_stat, l, q, i1, i2, iup1, iup2, m11, &
         m12, m21, m22
    real(kind=d_t), dimension(num_moments_f), intent(in) :: &
         upstream_moments
    real(kind=d_t), dimension(num_moments_f), intent(out) :: &
         transformed_moments, transformed_flux
    real(kind=d_t), dimension(2), intent(in) :: af
    real(kind=d_t), dimension(2,2), intent(in) :: bf
    real(kind=d_t), dimension(num_moments_f) :: upstream_flux
    real(kind=d_t), dimension(num_moments_f,num_moments_f), &
         intent(in) :: Lf, Uf
    real(kind=d_t), dimension(3,3) :: T

  ! Create transformation matrix

    call project_face_moments(af,bf,T)

  ! Solve for the upstream expansion coefficients

!!    call back_substitution(num_moments_f,upstream_moments,Lf,Uf,&
!!         upstream_flux)
    call back_substitution(num_moments_f,upstream_moments,upstream_flux)

    transformed_moments=MATMUL(T,upstream_flux)

!!    call back_substitution(num_moments_f,transformed_moments,Lf,Uf,&
!!         transformed_flux)
    call back_substitution(num_moments_f,transformed_moments,transformed_flux)

  end subroutine transform_incoming_moments

  subroutine transform_boundary_moments(Lf,Uf,transformed_moments,transformed_flux)
  !*********************************************************************
  !
  ! Subroutine transform incoming moments transforms the upstream cell
  ! face outgoing moments into downstream incoming cell face moments
  !
  !*********************************************************************

  ! Define variables

    real(kind=d_t), dimension(num_moments_f), intent(in) :: &
         transformed_moments
    real(kind=d_t), dimension(num_moments_f), intent(out) :: &
         transformed_flux
    real(kind=d_t), dimension(num_moments_f,num_moments_f), &
         intent(in) :: Lf, Uf

!!    call back_substitution(num_moments_f,transformed_moments,Lf,Uf,&
!!         transformed_flux)
    call back_substitution(num_moments_f,transformed_moments,transformed_flux)

  end subroutine transform_boundary_moments

  subroutine incoming_cell_subcell_project(Lf,Uf,af,bf,transformed_flux,incoming_flux)
  !*********************************************************************
  !
  ! Subroutine incoming cell to subcell moments projects the incoming
  ! cell face moments into incoming subcell face moments
  !
  !*********************************************************************

  ! Define variables

    integer(kind=li) :: alloc_stat, l, q, i1, i2, m11, m12, m21, m22
    real(kind=d_t), dimension(num_moments_f), intent(in) :: &
         transformed_flux
    real(kind=d_t), dimension(num_moments_f), intent(out) :: &
         incoming_flux
    real(kind=d_t), dimension(2), intent(in) :: af
    real(kind=d_t), dimension(2,2), intent(in) :: bf
    real(kind=d_t), dimension(num_moments_f) :: incoming_moments
    real(kind=d_t), dimension(num_moments_f,num_moments_f) :: Lf, Uf
    real(kind=d_t), dimension(3,3)  :: p

  ! Create projection matrix

    call project_face_moments(af,bf,p)    

    incoming_moments=MATMUL(p,transformed_flux)

!!    call back_substitution(num_moments_f,incoming_moments,Lf,Uf,&
!!         incoming_flux)
    call back_substitution(num_moments_f,incoming_moments,incoming_flux)

  end subroutine incoming_cell_subcell_project

  subroutine cell_source_expansion(q_moments,LL,U,q_expansion)
  !*********************************************************************
  !
  ! Subroutine cell source expansion solves for the source expansion 
  ! coefficients based on the distributed source moments
  !
  !*********************************************************************

  ! Pass source derived type

    real(kind=d_t), dimension(num_moments_v), intent(in) :: q_moments

  ! Define temporary variables

    integer(kind=li) ::  l
    real(kind=d_t), dimension(num_moments_v) :: y, x
    real(kind=d_t), dimension(num_moments_v), intent(inout) :: &
         q_expansion
    real(kind=d_t), dimension(num_moments_v,num_moments_v),&
         intent(in) :: LL, U

    do l=1, num_moments_v
       y(l)=q_moments(l)
    end do

!!    call back_substitution(num_moments_v,y,LL,U,x)
    call back_substitution(num_moments_v,y,x)

    do l=1, num_moments_v
       q_expansion(l)=x(l)
    end do

  end subroutine cell_source_expansion

  subroutine source_projection(LL,U,a,b,cell_source,subcell_source)
  !*********************************************************************
  !
  ! Subroutine source projection projects cell source moments into
  ! subcell source moments
  !
  !*********************************************************************

  ! Define variables

    integer(kind=li) ::  alloc_stat, l, q, i1, i2, i3, m11, m12, m13, &
         m21, m22, m23, m31, m32, m33
    real(kind=d_t) :: fact, a1temp, a2temp, a3temp
    real(kind=d_t), dimension(3), intent(in) :: a
    real(kind=d_t), dimension(3,3), intent(in) :: b
    real(kind=d_t), dimension(num_moments_v) :: subcell_moments
    real(kind=d_t), dimension(num_moments_v), intent(in) :: &
         cell_source
    real(kind=d_t), dimension(num_moments_v), intent(out) :: &
         subcell_source
    real(kind=d_t), dimension(num_moments_v,num_moments_v),&
         intent(in) :: LL, U
    real(kind=d_t), dimension(4,4) :: P

  ! Create projection matrix

    call project_moments(a,b,P)

    subcell_moments=MATMUL(P,cell_source)

!!    call back_substitution(num_moments_v,subcell_moments,LL,U,&
!!         subcell_source)
    call back_substitution(num_moments_v,subcell_moments,subcell_source)

  end subroutine source_projection

  subroutine characteristic_solver(i,t,e,a,b,af,bf,incoming_flux,subcell_source,&
                                   outgoing_moments,subcell_flux)
  !*********************************************************************
  !
  ! Subroutine characteristic solver computes the outgoing face 
  ! angular flux based on incoming and source subcell moments
  !
  !*********************************************************************

  ! Define variables

    integer(kind=li), intent(in) :: i
    integer(kind=li) ::  alloc_stat, l, q, i1, i2, i3, m11, m12, m13, &
         m21, m22, m23, m31, m32, m33, g1, g2, g3, g3p
    real(kind=d_t) :: fact, af1temp, af2temp, &
         a1temp, a2temp, a3temp
    real(kind=d_t), intent(in) :: t, e
    real(kind=d_t), dimension(2), intent(in) :: af
    real(kind=d_t), dimension(2,2), intent(in) :: bf
    real(kind=d_t), dimension(3), intent(in) :: a
    real(kind=d_t), dimension(3,3), intent(in) :: b
    real(kind=d_t), dimension(num_moments_f), intent(in) :: &
         incoming_flux
    real(kind=d_t), dimension(num_moments_v), intent(in) :: &
         subcell_source
    real(kind=d_t), dimension(num_moments_v), intent(out) :: &
         subcell_flux 
    real(kind=d_t), dimension(num_moments_f), intent(out) :: &
         outgoing_moments
    real(kind=d_t) :: CT_face_moments(3)
    real(kind=d_t) :: CT_moments(4)
    real(kind=d_t) :: PF(3,3)
    real(kind=d_t) :: PV(4,4)
    real(kind=d_t) :: FF(3,3)
    real(kind=d_t) :: FV(3,4)
    real(kind=d_t) :: F(4,3)
    real(kind=d_t) :: V(4,4)
    real(kind=d_t) :: e1,e2,e3,e4,e5,e6,expe,expne

  ! Pre-compute optimisers

    e1=e
    e2=e*e1
    e3=e*e2
    e4=e*e3
    e5=e*e4
    e6=e*e5

    expe=exp(e)
    expne=exp(-e)

  ! Create face integral matrix (outgoing face moments)

    if(e < 0.1_d_t)then
       FF(1,1)=1-0.3333333333333333*e1+0.08333333333333333*e2-0.016666666666666666*e3+&
            0.002777777777777778*e4
       FF(1,2)=0.6666666666666666-0.25*e1+0.06666666666666667*e2-0.013888888888888888*e3+&
            0.002380952380952381*e4
       FF(1,3)=0.3333333333333333-0.16666666666666666*e1+0.05*e2-0.011111111111111112*e3+&
            0.001984126984126984*e4

       FF(2,2)=0.5-0.2*e1+0.05555555555555555*e2-0.011904761904761904*e3+0.0020833333333333333*e4
       FF(2,3)=0.25-0.13333333333333333*e1+0.041666666666666664*e2-0.009523809523809525*e3+&
            0.001736111111111111*e4

       FF(3,3)=0.16666666666666666-0.1*e1+0.03333333333333333*e2-0.007936507936507936*e3+&
            0.001488095238095238*e4
    else
       FF(1,1)=2.0_d_t*(-1.0_d_t+expne+e)/e2
       FF(1,2)=2.0_d_t*(1.0_d_t/(2.0_d_t*e)+(1.0_d_t-expe+e)*expne/e3)
       FF(1,3)=(2.0_d_t*(2.0_d_t+expe*(-2_d_t+e)+e))*expne/e3

       FF(2,2)=2.0_d_t*(1.0_d_t/(3.0_d_t*e)+(2.0_d_t-2.0_d_t*expe+2.0_d_t*e+e2)*expne/e4)
       FF(2,3)=(expe*(-6.0_d_t+e2)+2.0_d_t*(3.0_d_t+3.0_d_t*e+e2))*expne/e4

       FF(3,3)=(2.0_d_t*(6.0_d_t+2.0_d_t*expe*(-3.0_d_t+e)+4.0_d_t*e+e2))*expne/e4
    end if

    FF(2,1)=FF(1,2)
    FF(3,1)=FF(1,3)
    FF(3,2)=FF(2,3)

  ! Create volume integral matrix (outgoing face moments)

    if(e < 0.1_d_t)then
       FV(1,1)=t*(0.3333333333333333-0.08333333333333333*e1+0.016666666666666666*e2-&
            0.002777777777777778*e3+0.0003968253968253968*e4)
       FV(1,2)=t*(0.25-0.06666666666666667*e1+0.013888888888888888*e2-0.002380952380952381*e3+&
            0.00034722222222222224*e4)
       FV(1,3)=t*(0.16666666666666666-0.05*e1+0.011111111111111112*e2-0.001984126984126984*e3+&
            0.00029761904761904765*e4)
       FV(1,4)=t*(0.08333333333333333-0.016666666666666666*e1+0.002777777777777778*e2-&
            0.0003968253968253968*e3+0.0000496031746031746*e4)

       FV(2,2)=t*(0.2-0.05555555555555555*e1+0.011904761904761904*e2-0.0020833333333333333*e3+&
            0.00030864197530864197*e4)
       FV(2,3)=t*(0.13333333333333333-0.041666666666666664*e1+0.009523809523809525*e2-&
            0.001736111111111111*e3+0.00026455026455026457*e4)
       FV(2,4)=t*(0.06666666666666667-0.013888888888888888*e1+0.002380952380952381*e2-&
            0.00034722222222222224*e3+0.00004409171075837743*e4)

       FV(3,3)=t*(0.1-0.03333333333333333*e1+0.007936507936507936*e2-0.001488095238095238*e3+&
            0.0002314814814814815*e4)
       FV(3,4)=t*(0.05-0.011111111111111112*e1+0.001984126984126984*e2-0.00029761904761904765*e3+&
            0.000038580246913580246*e4)
    else
       FV(1,1)=(t*(2.0_d_t-2.0_d_t*expne-2.0_d_t*e+e2))/e3
       FV(1,2)=(t*(6.0_d_t-3.0_d_t*e2+2.0_d_t*e3-(6.0_d_t*(1.0_d_t+e))*expne))/(3.0_d_t*e4)
       FV(1,3)=(t*(-6.0_d_t*(2.0_d_t+e)+expe*(12.0_d_t-6.0_d_t*e+e3)))*expne/(3.0_d_t*e4)
       FV(1,4)=(t*(-6.0_d_t+6.0_d_t*expne+e*(6.0_d_t+(-3.0_d_t+e)*e)))/(3.0_d_t*e4)

       FV(2,2)=2.0_d_t*t*(-1.0_d_t/(3.0_d_t*e2)+1.0_d_t/(4.0_d_t*e)+&
            (2.0_d_t-(2.0_d_t+2.0_d_t*e+e2)*expne)/e5)
       FV(2,3)=(t*(-8.0_d_t*(3.0_d_t+3.0_d_t*e+e2)+expe*(24.0_d_t-4.0_d_t*e2+e4)))/&
             (4.0_d_t*expe*e5)
       FV(2,4)=(t*(-24.0_d_t+12.0_d_t*e2-8.0_d_t*e3+3.0_d_t*e4+(24.0_d_t*(1.0_d_t+e))/&
             expe))/(12.0_d_t*e5)

       FV(3,3)=(t*(-12.0_d_t*(6.0_d_t+4.0_d_t*e+e2)+expe*(72.0_d_t-24.0_d_t*e+e4)))/&
            (6.0_d_t*expe*e5)
       FV(3,4)=(t*((12.0_d_t*(2.0_d_t+e))*expne+(-2.0_d_t+e)*(12.0_d_t+e3)))/(6.0_d_t*e5)
    end if

    FV(2,1)=FV(1,2)
    FV(3,1)=FV(1,3)
    FV(3,2)=FV(2,3)

  ! Create face integral matrix (subcell flux moments)

    if(e < 0.1_d_t)then
       F(1,1)=1-0.25*e1+0.05*e2-0.008333333333333333*e3+0.0011904761904761906*e4
       F(1,2)=0.75-0.2*e1+0.041666666666666664*e2-0.007142857142857143*e3+0.0010416666666666667*e4
       F(1,3)=0.5-0.15*e1+0.03333333333333333*e2-0.005952380952380952*e3+0.0008928571428571428*e4

       F(2,2)=0.6-0.16666666666666666*e1+0.03571428571428571*e2-0.00625*e3+0.000925925925925926*e4
       F(2,3)=0.4-0.125*e1+0.02857142857142857*e2-0.005208333333333333*e3+0.0007936507936507937*e4

       F(3,3)=0.3-0.1*e1+0.023809523809523808*e2-0.004464285714285714*e3+0.0006944444444444445*e4

       F(4,1)=0.25-0.1*e1+0.025*e2-0.004761904761904762*e3+0.000744047619047619*e4
       F(4,2)=0.2-0.08333333333333333*e1+0.02142857142857143*e2-0.004166666666666667*e3+&
            0.0006613756613756613*e4
       F(4,3)=0.15-0.06666666666666667*e1+0.017857142857142856*e2-0.0035714285714285713*e3+&
            0.0005787037037037037*e4
    else
       F(1,1)=(3.0_d_t*(2.0_d_t-2.0_d_t*expne-2.0_d_t*e+e2))/e3
       F(1,2)=(6.0_d_t-3.0_d_t*e2+2.0_d_t*e3-(6.0_d_t*(1.0_d_t+e))*expne)/e4
       F(1,3)=(-6.0_d_t*(2.0_d_t+e)+expe*(12.0_d_t-6.0_d_t*e+e3))*expne/e4

       F(2,2)=6.0_d_t*(-1.0_d_t/(3.0_d_t*e2)+1.0_d_t/(4.0_d_t*e)+(2.0_d_t-&
            (2.0_d_t+2.0_d_t*e+e2)*expne)/e5)
       F(2,3)=(3.0_d_t*(-8.0_d_t*(3.0_d_t+3.0_d_t*e+e2)+expe*(24.0_d_t-4.0_d_t*e2+e4)))/&
            (4.0_d_t*expe*e5)

       F(3,3)=(-12.0_d_t*(6.0_d_t+4.0_d_t*e+e2)+expe*(72.0_d_t-24.0_d_t*e+e4))/&
            (2.0_d_t*expe*e5)

       F(4,1)=(3.0_d_t*(-2.0_d_t*(3.0_d_t+e)+expe*(6.0_d_t-4.0_d_t*e+e2)))*expne/e4
       F(4,2)=(2.0_d_t*(-3.0_d_t*(2.0_d_t+e)**2+expe*(12.0_d_t-3.0_d_t*e2+e3)))*expne/e5
       F(4,3)=(-6.0_d_t*(8.0_d_t+5.0_d_t*e+e2)+expe*(48.0_d_t-18.0_d_t*e+e3))*expne/e5
    end if

    F(2,1)=F(1,2)
    F(3,1)=F(1,3)
    F(3,2)=F(2,3)

  ! Create volume integral matrix (subcell flux moments)

    if(e < 0.1_d_t)then
       V(1,1)=t*(0.25-0.05*e1+0.008333333333333333*e2-0.0011904761904761906*e3+0.00014880952380952382*e4)
       V(1,2)=t*(0.2-0.041666666666666664*e1+0.007142857142857143*e2-0.0010416666666666667*e3+&
            0.00013227513227513228*e4)
       V(1,3)=t*(0.15-0.03333333333333333*e1+0.005952380952380952*e2-0.0008928571428571428*e3+&
            0.00011574074074074075*e4)
       V(1,4)=t*(0.05-0.008333333333333333*e1+0.0011904761904761906*e2-0.00014880952380952382*e3+&
            0.000016534391534391536*e4)

       V(2,2)=t*(0.16666666666666666-0.03571428571428571*e1+0.00625*e2-0.000925925925925926*e3+&
            0.00011904761904761905*e4)
       V(2,3)=t*(0.125-0.02857142857142857*e1+0.005208333333333333*e2-0.0007936507936507937*e3+&
            0.00010416666666666667*e4)
       V(2,4)=t*(0.041666666666666664-0.007142857142857143*e1+0.0010416666666666667*e2-&
            0.00013227513227513228*e3+0.000014880952380952381*e4)

       V(3,3)=t*(0.1-0.023809523809523808*e1+0.004464285714285714*e2-0.0006944444444444445*e3+&
            0.00009259259259259259*e4)
       V(3,4)=t*(0.03333333333333333-0.005952380952380952*e1+0.0008928571428571428*e2-&
            0.00011574074074074075*e3+0.000013227513227513228*e4)

       V(4,1)=t*(0.1-0.025*e1+0.004761904761904762*e2-0.000744047619047619*e3+&
            0.0000992063492063492*e4)
       V(4,2)=t*(0.08333333333333333-0.02142857142857143*e1+0.004166666666666667*e2-&
            0.0006613756613756613*e3+0.00008928571428571429*e4)
       V(4,3)=t*(0.06666666666666667-0.017857142857142856*e1+0.0035714285714285713*e2-&
            0.0005787037037037037*e3+0.00007936507936507937*e4)
       V(4,4)=t*(0.025-0.004761904761904762*e1+0.000744047619047619*e2-0.0000992063492063492*e3+&
            0.000011574074074074073*e4)
    else
       V(1,1)=(t*(-6.0_d_t+6.0_d_t*expne+e*(6.0_d_t+(-3.0_d_t+e)*e)))/e4
       V(1,2)=(t*(-24.0_d_t+12.0_d_t*e2-8.0_d_t*e3+3.0_d_t*e4+(24.0_d_t*(1.0_d_t+e))*expne))/&
           (4.0_d_t*e5)
       V(1,3)=(t*((12.0_d_t*(2.0_d_t+e))*expne+(-2.0_d_t+e)*(12.0_d_t+e3)))/(2.0_d_t*e5)
       V(1,4)=(t*(24.0_d_t-24.0_d_t*expne+e*(-24.0_d_t+e*(12.0_d_t+(-4.0_d_t+e)*e))))/&
           (4.0_d_t*e5)

       V(2,2)=(t*(-120.0_d_t+20.0_d_t*e3-15.0_d_t*e4+6.0_d_t*e5+(60.0_d_t*(2.0_d_t+2.0_d_t*e+e2))/&
           expe))/(10.0_d_t*e6)
       V(2,3)=(t*(-360.0_d_t+60.0_d_t*e2-15.0_d_t*e4+8.0_d_t*e5+(120.0_d_t*(3.0_d_t+e*(3.0_d_t+e)))/&
           expe))/(20.0_d_t*e6)
       V(2,4)=(t*(120.0_d_t-(120.0_d_t*(1+e))*expne+e2*(-60.0_d_t+e*(40.0_d_t+e*(-15.0_d_t+&
           4.0_d_t*e)))))/(20.0_d_t*e6)

       V(3,3)=(t*(120.0_d_t*(-3.0_d_t+e)+e4*(-5.0_d_t+3.0_d_t*e)+(60.0_d_t*(6.0_d_t+e*(4.0_d_t+e)))/&
           expe))/(10.0_d_t*e6)
       V(3,4)=(t*(-120.0_d_t*(-2.0_d_t+e)-(120.0_d_t*(2.0_d_t+e))*expne+e3*(20.0_d_t+e*(-10.0_d_t+&
           3.0_d_t*e))))/(20.0_d_t*e6)

       V(4,1)=(t*(-72.0_d_t+(24.0_d_t*(3.0_d_t+e))*expne+e*(48.0_d_t-12.0_d_t*e+e3)))/(4.0_d_t*e5)
       V(4,2)=(t*(-120.0_d_t+(30.0_d_t*(2.0_d_t+e)**2)*expne+e2*(30.0_d_t-10.0_d_t*e+e3)))/&
           (5.0_d_t*e6) 
       V(4,3)=(t*(-960.0_d_t+360.0_d_t*e-20.0_d_t*e3+3.0_d_t*e5+(120.0_d_t*(8.0_d_t+e*(5.0_d_t+e)))/&
           expe))/(20.0_d_t*e6) 
       V(4,4)=(t*(360.0_d_t-(120.0_d_t*(3.0_d_t+e))*expne+e*(-240.0_d_t+60.0_d_t*e-5.0_d_t*e3+&
           2.0_d_t*e4)))/(20.0_d_t*e6)
    end if

    V(2,1)=V(1,2)
    V(3,1)=V(1,3)
    V(3,2)=V(2,3)

  ! Create face projection matrix 

    call project_face_moments(af,bf,PF)

  ! Create volume projection matrix 

    call project_moments(a,b,PV)

  ! Compute outgoing face and volume moment contribution from CT to Cell 

    CT_face_moments=MATMUL(FF,incoming_flux)+MATMUL(FV,subcell_source)

    outgoing_moments(1)=CT_face_moments(1)
    outgoing_moments(2)=af(1)*CT_face_moments(1)+bf(1,1)*CT_face_moments(2)+bf(1,2)*CT_face_moments(3)
    outgoing_moments(3)=af(2)*CT_face_moments(1)+bf(2,1)*CT_face_moments(2)+bf(2,2)*CT_face_moments(3)

    CT_moments=MATMUL(F,incoming_flux)+MATMUL(V,subcell_source)

    subcell_flux(1)=CT_moments(1)
    subcell_flux(2)=a(1)*CT_moments(1)+b(1,1)*CT_moments(2)+b(1,2)*CT_moments(3)+b(1,3)*CT_moments(4)
    subcell_flux(3)=a(2)*CT_moments(1)+b(2,1)*CT_moments(2)+b(2,2)*CT_moments(3)+b(2,3)*CT_moments(4)
    subcell_flux(4)=a(3)*CT_moments(1)+b(3,1)*CT_moments(2)+b(3,2)*CT_moments(3)+b(3,3)*CT_moments(4)

  end subroutine characteristic_solver

  subroutine face_moment_transformation(Lf,Uf,f,face_angular_mom,face_cell_temp)
  !*********************************************************************
  !
  ! Subroutine face moment transformation transforms the face angular
  ! flux moments from their face coordinate systems into the cell
  ! volume coordinate system through matrix multiplication
  !
  !*********************************************************************

  ! Define variables

    integer(kind=li) :: f
    integer(kind=li) :: alloc_stat, l, q, i1F, i2F
    real(kind=d_t), dimension(num_moments_f), intent(in) :: &
         face_angular_mom
    real(kind=d_t), dimension(num_moments_v), intent(out) :: &
         face_cell_temp
    real(kind=d_t), dimension(num_moments_f) :: x
    real(kind=d_t), dimension(num_moments_f,num_moments_f) :: Lf, Uf
    real(kind=d_t), dimension(4,3) :: TF

  ! Create projection matrix

    do q=1, num_moments_v
       if(f == 0)then
          i1F=index_v(q)%i2
          i2F=index_v(q)%i3
       elseif(f == 1)then
          i1F=index_v(q)%i1+index_v(q)%i2
          i2F=index_v(q)%i3
       elseif(f == 2)then
          i1F=index_v(q)%i1
          i2F=index_v(q)%i2+index_v(q)%i3
       elseif(f == 3)then
          i1F=index_v(q)%i1
          i2F=index_v(q)%i2
       else
          call stop_thor(7_li)
       endif
       do l=1, num_moments_f
          TF(q,l)=2.0_d_t/REAL((index_f(l)%i1+i1F+index_f(l)%i2+&
               i2F+2.0_d_t)*(index_f(l)%i2+i2F+1.0_d_t))
          if(f == 3 .and. index_v(q)%i3 /= 0)then
             TF(q,l)=0.0_d_t
          end if
       end do
    end do
    
  ! Solve for expansion coefficient on face

!!    call back_substitution(num_moments_f,face_angular_mom,Lf,Uf,x)
    call back_substitution(num_moments_f,face_angular_mom,x)
    
    face_cell_temp=MATMUL(TF,x)

  end subroutine face_moment_transformation

  subroutine back_substitution(n,xin,xout)
  !*********************************************************************
  !
  ! Subroutine back_substitutions performs back substitution based on
  ! either 3x3 or 4x4 matrix 
  !
  !*********************************************************************

    implicit none

  ! Define variables

    integer(kind=li), intent(in)  :: n
    real(kind=d_t),   intent(in)  :: xin(n)
    real(kind=d_t),   intent(out) :: xout(n)

  ! Determine array size

    if(n == 3)then
       xout(1)=9.0_d_t*xin(1)-12.0_d_t*xin(2)
       xout(2)=-12.0_d_t*xin(1)+24.0_d_t*xin(2)-12.0_d_t*xin(3)
       xout(3)=-12.0_d_t*xin(2)+24.0_d_t*xin(3)
    elseif(n == 4)then
       xout(1)=16.0_d_t*xin(1)-20.0_d_t*xin(2)
       xout(2)=-20.0_d_t*xin(1)+40.0_d_t*xin(2)-20.0_d_t*xin(3)
       xout(3)=-20.0_d_t*xin(2)+40.0_d_t*xin(3)-20.0_d_t*xin(4)
       xout(4)=-20.0_d_t*xin(3)+40.0_d_t*xin(4)
    end if

  end subroutine back_substitution

!!  subroutine back_substitution(n,b,L,U,x)
!!  !**********************************************************************
!!  !
!!  ! Subroutine back substitution solves for unknown after LU
!!  !
!!  !**********************************************************************
!!  ! Pass input parameters
!!
!!    integer(kind=li), intent(in) :: n
!!    real(kind=d_t), dimension(n), intent(in) :: b
!!    real(kind=d_t), dimension(n,n), intent(in) :: L, U
!!    real(kind=d_t), dimension(n), intent(out) :: x
!!
!!  ! Declare spatial order moments index
!!
!!    integer(kind=li) :: i, j
!!    real(kind=d_t), dimension(n) :: y
!!
!!    do i=1, n
!!       y(i)=b(i)
!!       do j=1, i-1
!!          y(i)=y(i)-L(i,j)*y(j)
!!       end do
!!    end do
!!
!!    do i=n, 1, -1
!!       x(i)=y(i)/U(i,i)
!!       do j=i+1, n
!!          x(i)=x(i)-U(i,j)*x(j)/U(i,i)
!!       end do
!!    end do
!!
!!  end subroutine back_substitution

  subroutine project_face_moments(af,bf,M)
  !*********************************************************************
  !
  ! Subroutine project face moments build generic transformation  
  ! matrix 
  !
  !*********************************************************************

    implicit none

  ! Define temporary variables

    real(kind=d_t), intent(in)  :: af(2)
    real(kind=d_t), intent(in)  :: bf(2,2)
    real(kind=d_t), intent(out) :: M(3,3)

  ! Create transformation matrix

    M(1,1)=1.0_d_t
    M(1,2)=0.3333333333333333_d_t*(3.0_d_t*af(1)+2.0_d_t*bf(1,1)+bf(1,2))
    M(1,3)=0.3333333333333333_d_t*(3.0_d_t*af(2)+2.0_d_t*bf(2,1)+bf(2,2))

    M(2,1)=0.6666666666666666_d_t
    M(2,2)=0.08333333333333333_d_t*(8.0_d_t*af(1)+6.0_d_t*bf(1,1)+3.0_d_t*bf(1,2))
    M(2,3)=0.08333333333333333_d_t*(8.0_d_t*af(2)+6.0_d_t*bf(2,1)+3.0_d_t*bf(2,2))

    M(3,1)=0.3333333333333333_d_t
    M(3,2)=0.08333333333333333_d_t*(4.0_d_t*af(1)+3.0_d_t*bf(1,1)+2.0_d_t*bf(1,2))
    M(3,3)=0.08333333333333333_d_t*(4.0_d_t*af(2)+3.0_d_t*bf(2,1)+2.0_d_t*bf(2,2))

  end subroutine project_face_moments

  subroutine project_moments(a,b,M)
  !*********************************************************************
  !
  ! Subroutine project moments build generic transformation matrix 
  !
  !*********************************************************************

    implicit none

  ! Define temporary variables

    real(kind=d_t), intent(in)  :: a(3)
    real(kind=d_t), intent(in)  :: b(3,3)
    real(kind=d_t)              :: M(4,4)

  ! Create projection matrix

    M(1,1)=1.0_d_t
    M(1,2)=0.25_d_t*(4.0_d_t*a(1)+3.0_d_t*b(1,1)+2.0_d_t*b(1,2)+b(1,3))
    M(1,3)=0.25_d_t*(4.0_d_t*a(2)+3.0_d_t*b(2,1)+2.0_d_t*b(2,2)+b(2,3))
    M(1,4)=0.25_d_t*(4.0_d_t*a(3)+3.0_d_t*b(3,1)+2.0_d_t*b(3,2)+b(3,3))

    M(2,1)=0.75_d_t
    M(2,2)=0.05_d_t*(15.0*a(1)+4.0_d_t*(3.0_d_t*b(1,1)+2.0_d_t*b(1,2)+b(1,3)))
    M(2,3)=0.05_d_t*(15.0*a(2)+4.0_d_t*(3.0_d_t*b(2,1)+2.0_d_t*b(2,2)+b(2,3)))
    M(2,4)=0.05_d_t*(15.0*a(3)+4.0_d_t*(3.0_d_t*b(3,1)+2.0_d_t*b(3,2)+b(3,3)))

    M(3,1)=0.5_d_t
    M(3,2)=0.05_d_t*(10.0_d_t*a(1)+8.0_d_t*b(1,1)+6.0_d_t*b(1,2)+3.0_d_t*b(1,3))
    M(3,3)=0.05_d_t*(10.0_d_t*a(2)+8.0_d_t*b(2,1)+6.0_d_t*b(2,2)+3.0_d_t*b(2,3))
    M(3,4)=0.05_d_t*(10.0_d_t*a(3)+8.0_d_t*b(3,1)+6.0_d_t*b(3,2)+3.0_d_t*b(3,3))

    M(4,1)=0.25_d_t
    M(4,2)=0.05_d_t*(5.0_d_t*a(1)+4.0_d_t*b(1,1)+3.0_d_t*b(1,2)+2.0_d_t*b(1,3))
    M(4,3)=0.05_d_t*(5.0_d_t*a(2)+4.0_d_t*b(2,1)+3.0_d_t*b(2,2)+2.0_d_t*b(2,3))
    M(4,4)=0.05_d_t*(5.0_d_t*a(3)+4.0_d_t*b(3,1)+3.0_d_t*b(3,2)+2.0_d_t*b(3,3))

  end subroutine project_moments

end module transport_kernel_module_LC

