module transport_kernel_module_CCE
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
  
  subroutine transport_kernel_CCE(sigmat,q_moments,vol_moment,face_moment,       &
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
  

  ! Define temporary variables

    integer(kind=li) :: alloc_stat, ii, l, f, adjcnt_cell
    real(kind=d_t) :: det_Jup, det_Js, e
    real(kind=d_t), dimension(2) :: af
    real(kind=d_t), dimension(3) :: a_temp, a
    real(kind=d_t), dimension(2,2) :: bf
    real(kind=d_t), dimension(3,2) :: JFup, JFdown, JF, Jsf
    real(kind=d_t), dimension(2,3) :: JFup_inv, JF_inv
    real(kind=d_t), dimension(3,3) :: Jup, Jup_inv, Js, Js_inv, b
    real(kind=d_t), dimension(num_moments_v) :: q_expansion
    real(kind=d_t), dimension(:), allocatable :: area, volume, &
         subcell_upstream_moments, transformed_moments, &
         transformed_flux, incoming_flux, outgoing_moments,cell_source, &
         subcell_source, subcell_flux, projected_temp, proj_moments, &
         proj_flux_moments, face_angular_mom, face_cell_temp, y
    real(kind=d_t), dimension(:,:), allocatable :: &
         subcell_source_moments, projected_moments,&
         face_angular, projected_flux
    type(vector) :: R0down, R0up, R0F, v0, v1, v2, v3
    type(vector), dimension(0:3) :: v

  ! Allocate variables

    allocate(area(subcells),volume(subcells),&
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
    if(alloc_stat /= 0) call stop_thor(2_li)

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
                                          transformed_moments,transformed_flux)

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

    deallocate(area,volume,subcell_upstream_moments,incoming_flux,&
         transformed_moments,transformed_flux,cell_source,&
         subcell_source,outgoing_moments,projected_moments,&
         projected_temp,proj_moments,projected_flux,&
         proj_flux_moments,face_angular,subcell_flux,&
         subcell_source_moments)
    
  end subroutine transport_kernel_CCE

  subroutine transform_incoming_moments(Lf,Uf,af,bf,upstream_moments,transformed_moments,transformed_flux)
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
    real(kind=d_t), dimension(:,:), allocatable :: T

  ! Allocate transformation matrix 

    allocate(T(num_moments_f,num_moments_f),stat=alloc_stat);&
         T=0.0_d_t
    if(alloc_stat /= 0) call stop_thor(2_li)

  ! Create transformation matrix

    do q=1, num_moments_f
       i1=index_f(q)%i1
       i2=index_f(q)%i2
       do l=1, num_moments_f
          iup1=index_f(l)%i1
          iup2=index_f(l)%i2
          do m11=0, iup1
             do m12=0, iup1-m11
                do m21=0, iup2
                   do m22=0, iup2-m21
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
                   end do
                end do
             end do
          end do
       end do
    end do

  ! Solve for the upstream expansion coefficients

    call back_substitution(num_moments_f,upstream_moments,Lf,Uf,&
         upstream_flux)

  ! Solve for the downstream expansion coefficients

    transformed_moments=MATMUL(T,upstream_flux)

    call back_substitution(num_moments_f,transformed_moments,Lf,Uf,&
         transformed_flux)

    deallocate(T)

  end subroutine transform_incoming_moments

  subroutine transform_boundary_moments(Lf,Uf,transformed_moments,transformed_flux)
  !*********************************************************************
  !
  ! Subroutine transform incoming moments transforms the upstream cell
  ! face outgoing moments into downstream incoming cell face moments
  !
  !*********************************************************************

  ! Define temporary variables

    real(kind=d_t), dimension(num_moments_f), intent(in) :: &
         transformed_moments
    real(kind=d_t), dimension(num_moments_f), intent(out) :: &
         transformed_flux
    real(kind=d_t), dimension(num_moments_f,num_moments_f), &
         intent(in) :: Lf, Uf

    call back_substitution(num_moments_f,transformed_moments,Lf,Uf,&
         transformed_flux)

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
    real(kind=d_t), dimension(:,:), allocatable :: p

  ! Allocate projection matrix

    allocate(p(num_moments_f,num_moments_f),stat=alloc_stat);&
         p=0.0_d_t
    if(alloc_stat /= 0) call stop_thor(2_li)

  ! Create projection matrix
    
    do q=1, num_moments_f
       do l=1, num_moments_f
          i1=index_f(l)%i1
          i2=index_f(l)%i2
          do m11=0, i1
             do m12=0, i1-m11
                do m21=0, i2
                   do m22=0, i2-m21
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
                   end do
                end do
             end do
          end do
       end do
    end do

    incoming_moments=MATMUL(p,transformed_flux)

    call back_substitution(num_moments_f,incoming_moments,Lf,Uf,&
         incoming_flux)

    deallocate(p)

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

    call back_substitution(num_moments_v,y,LL,U,x)

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
    real(kind=d_t), dimension(:,:), allocatable :: P

  ! Allocate projection matrix

    allocate(P(num_moments_v,num_moments_v),stat=alloc_stat);&
         P=0.0_d_t;
    if(alloc_stat /= 0) call stop_thor(2_li)

  ! Create projection matrix

    do q=1, num_moments_v
       do l=1, num_moments_v
          i1=index_v(l)%i1
          i2=index_v(l)%i2
          i3=index_v(l)%i3
          fact=factorial_d_t(i1)*factorial_d_t(i2)*factorial_d_t(i3)
          do m11=0, i1
             do m12=0, i1-m11
                do m13=0, i1-m11-m12
                   a1temp=a(1)**(i1-m11-m12-m13)
                   do m21=0, i2
                      do m22=0, i2-m21
                         do m23=0, i2-m21-m22
                            a2temp=a(2)**(i2-m21-m22-m23)
                            do m31=0, i3
                               do m32=0, i3-m31
                                  do m33=0, i3-m31-m32
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
                                  end do
                               end do
                            end do
                         end do
                      end do
                   end do
                end do
             end do
          end do
       end do
    end do

    subcell_moments=MATMUL(P,cell_source)

    call back_substitution(num_moments_v,subcell_moments,LL,U,&
         subcell_source)

    deallocate(P)

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
    real(kind=d_t) :: e1, e2, e3, e4, e5, e6, fact, af1temp, af2temp, &
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
    real(kind=d_t), dimension(:,:), allocatable :: FF, FV, F, V

  ! Allocate face and volume matrix

    allocate(FF(num_moments_f,num_moments_f),&
         FV(num_moments_f,num_moments_v),F(num_moments_v,num_moments_f),&
         V(num_moments_v,num_moments_v),stat=alloc_stat);&
         FF=0.0_d_t;FV=0.0_d_t;F=0.0_d_t;V=0.0_d_t
    if(alloc_stat /= 0) call stop_thor(2_li)

  ! Pre-compute optimisers

    e1=e
    e2=e*e1
    e3=e*e2
    e4=e*e3
    e5=e*e4
    e6=e*e5

  ! Create face integral matrix (outgoing face moments)

    do q=1, num_moments_f
       i1=index_f(q)%i1
       i2=index_f(q)%i2
       fact=factorial_d_t(i1)*factorial_d_t(i2)
       do l=1, num_moments_f
          do m11=0, i1
             do m12=0, i1-m11
                af1temp=af(1)**(i1-m11-m12)
                do m21=0, i2
                   do m22=0, i2-m21
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
                   end do
                end do
             end do
          end do
       end do
    end do

  ! Create volume integral matrix (outgoing face moments)

    do q=1, num_moments_f
       i1=index_f(q)%i1
       i2=index_f(q)%i2
       fact=factorial_d_t(i1)*factorial_d_t(i2)
       do l=1, num_moments_v
          do m11=0, i1
             do m12=0, i1-m11
                af1temp=af(1)**(i1-m11-m12)
                do m21=0, i2
                   do m22=0, i2-m21
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
                   end do
                end do
             end do
          end do
       end do
    end do

  ! Create face integral matrix (subcell flux moments)

    do q=1, num_moments_v
       i1=index_v(q)%i1
       i2=index_v(q)%i2
       i3=index_v(q)%i3
       fact=factorial_d_t(i1)*factorial_d_t(i2)*factorial_d_t(i3)
       do l=1, num_moments_f
          do m11=0, i1
             do m12=0, i1-m11
                do m13=0, i1-m11-m12
                   a1temp=a(1)**(i1-m11-m12-m13)
                   do m21=0, i2
                      do m22=0, i2-m21
                         do m23=0, i2-m21-m22
                            a2temp=a(2)**(i2-m21-m22-m23)
                            do m31=0, i3
                               do m32=0, i3-m31
                                  do m33=0, i3-m31-m32
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
                                  end do
                               end do
                            end do
                         end do
                      end do
                   end do
                end do
             end do
          end do
       end do
    end do

  ! Create volume integral matrix (subcell flux moments)

    do q=1, num_moments_v
       i1=index_v(q)%i1
       i2=index_v(q)%i2
       i3=index_v(q)%i3
       fact=factorial_d_t(i1)*factorial_d_t(i2)*factorial_d_t(i3)
       do l=1, num_moments_v
          do m11=0, i1
             do m12=0, i1-m11
                do m13=0, i1-m11-m12
                   a1temp=a(1)**(i1-m11-m12-m13)
                   do m21=0, i2
                      do m22=0, i2-m21
                         do m23=0, i2-m21-m22
                            a2temp=a(2)**(i2-m21-m22-m23)
                            do m31=0, i3
                               do m32=0, i3-m31
                                  do m33=0, i3-m31-m32
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
                                  end do
                               end do
                            end do
                         end do
                      end do
                   end do
                end do
             end do
          end do
       end do
    end do

    outgoing_moments=MATMUL(FF,incoming_flux)+MATMUL(FV,subcell_source)

    subcell_flux=MATMUL(F,incoming_flux)+MATMUL(V,subcell_source)

    deallocate(FF,FV,F,V)

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

    integer(kind=li), intent(in) :: f
    integer(kind=li) :: alloc_stat, l, q, i1F, i2F
    real(kind=d_t), dimension(num_moments_f), intent(in) :: &
         face_angular_mom
    real(kind=d_t), dimension(num_moments_v), intent(out) :: &
         face_cell_temp
    real(kind=d_t), dimension(num_moments_f) :: x
    real(kind=d_t), dimension(num_moments_f,num_moments_f) :: Lf, Uf
    real(kind=d_t), dimension(:,:), allocatable :: TF

  ! Allocate projection matrix

    allocate(TF(num_moments_v,num_moments_f),stat=alloc_stat);&
         TF=0.0_d_t
    if(alloc_stat /= 0) call stop_thor(2_li)

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

    call back_substitution(num_moments_f,face_angular_mom,Lf,Uf,x)
    
    face_cell_temp=MATMUL(TF,x)

    deallocate(TF)
    
  end subroutine face_moment_transformation

  subroutine back_substitution(n,b,L,U,x)
  !**********************************************************************
  !
  ! Subroutine back substitution solves for unknown after LU
  !
  !**********************************************************************
  ! Pass input parameters

    integer(kind=li), intent(in) :: n
    real(kind=d_t), dimension(n), intent(in) :: b
    real(kind=d_t), dimension(n,n), intent(in) :: L, U
    real(kind=d_t), dimension(n), intent(out) :: x

  ! Declare spatial order moments index

    integer(kind=li) :: i, j
    real(kind=d_t), dimension(n) :: y

    do i=1, n
       y(i)=b(i)
       do j=1, i-1
          y(i)=y(i)-L(i,j)*y(j)
       end do
    end do

    do i=n, 1, -1
       x(i)=y(i)/U(i,i)
       do j=i+1, n
          x(i)=x(i)-U(i,j)*x(j)/U(i,i)
       end do
    end do

  end subroutine back_substitution

  function face_moment1(gam1,gam2,e)
  !*********************************************************************
  !
  ! Subroutine calculates face moments of characteristic integral
  ! using a Taylor expansion
  !
  !*********************************************************************
    use types
    implicit none

  ! Pass input parameters
    integer(kind=li), intent(in) :: gam1, gam2
    real(kind=d_t), intent(in) :: e

  ! Define temporary variables
    integer(kind=li) :: n, max_n, ggam
    real(kind=d_t) :: r, new_term, face_moment1, rerror
    real(kind=d_t), parameter :: eps=1.0e-8

    max_n=1000
    ggam=gam1+gam2
    face_moment1=0.0_d_t
    rerror=1.0_d_t

    do n=0, max_n
       if(n == 0)then
          r=2.0_d_t
       else
          r=r*(-e)/n
       end if
       new_term=r/((ggam+n+2_li)*(gam2+n+1_li))
       face_moment1=face_moment1+new_term
       rerror=abs(new_term/face_moment1)
       if(rerror > eps)then
       else
          go to 10
       end if
    end do
    
10 end function face_moment1

  function face_moment2(gam1,gam2,gam3,t,e)
  !*********************************************************************
  !
  ! Subroutine calculates volume moments of characteristic integral
  ! using a Taylor expansion
  !
  !*********************************************************************
    use types
    implicit none
    integer(kind=li), intent(in) :: gam1, gam2, gam3
    real(kind=d_t), intent(in) :: t, e
    integer(kind=li) :: n, max_n, m, max_m, ggam3, ggam2, ggam1
    real(kind=d_t) :: r1, r2, new_term1, new_term2, face_moment2, &
         rerror1, rerror2
    real(kind=d_t), parameter :: eps=1.0e-8
    
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

    do n=0, max_n
       if(n == 0)then
          r1=1.0_d_t
       else
          r1=r1*e/n
       end if
       do m=0, max_m
          if(m == 0)then
             r2=2.0_d_t*t
          else
             r2=r2*(-e)/m
          end if
          new_term1=r1*r2/((ggam3+n+m+3_li)*(ggam2+n+m+2_li)*&
               (ggam1+n+1_li))
          new_term2=new_term2+new_term1
          if(new_term2 /= 0)then
             rerror1=abs(new_term1/new_term2)
          end if
          if(rerror1 > eps)then
          else
             go to 10
          end if
       end do

10     continue

       face_moment2=face_moment2+new_term2
       rerror2=abs(new_term2/face_moment2)
       if(rerror2 > eps)then
          new_term2=0.0_d_t
       else
          go to 11
       end if
    end do

11 end function face_moment2

  function volume_moment1(gam1,gam2,gam3,e)
  !*********************************************************************
  !
  ! Subroutine calculates face moments of characteristic integral
  ! using a Taylor expansion
  !
  !*********************************************************************
    use types
    implicit none

  ! Pass input parameters
    integer(kind=li), intent(in) :: gam1, gam2, gam3
    real(kind=d_t), intent(in) :: e

  ! Define temporary variables
    integer(kind=li) :: n, max_n, ggam3, ggam2
    real(kind=d_t) :: r, new_term, volume_moment1, rerror
    real(kind=d_t), parameter :: eps=1.0e-8

    max_n=1000
    ggam3=gam3+gam2+gam1
    ggam2=gam3+gam2
    volume_moment1=0.0_d_t
    rerror=1.0_d_t

    do n=0, max_n
       if(n == 0)then
          r=6.0_d_t
       else
          r=r*(-e)/n
       end if
       new_term=r/((ggam3+n+3_li)*(ggam2+n+2_li)*(gam3+n+1_li))
       volume_moment1=volume_moment1+new_term
       rerror=abs(new_term/volume_moment1)
       if(rerror > eps)then
       else
          go to 10
       end if
    end do
    
10 end function volume_moment1

  function volume_moment2(gam1,gam2,gam3,gam3p,t,e)
  !*********************************************************************
  !
  ! Subroutine calculates volume moments of characteristic integral
  ! using a Taylor expansion
  !
  !*********************************************************************
    use types
    implicit none
    integer(kind=li), intent(in) :: gam1, gam2, gam3, gam3p
    real(kind=d_t), intent(in) :: t, e
    integer(kind=li) :: n, max_n, m, max_m, ggam3, ggam2, ggam1
    real(kind=d_t) :: r1, r2, new_term1, new_term2, volume_moment2, &
         rerror1, rerror2
    real(kind=d_t), parameter :: eps=1.0e-8
    
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

    do n=0, max_n
       if(n == 0)then
          r1=1.0_d_t
       else
          r1=r1*e/n
       end if
       do m=0, max_m
          if(m == 0)then
             r2=6.0_d_t*t
          else
             r2=r2*(-e)/m
          end if
          new_term1=r1*r2/((ggam3+n+m+4_li)*(ggam2+n+m+3_li)*&
               (ggam1+n+m+2_li)*(gam3p+n+1_li))
          new_term2=new_term2+new_term1
          if(new_term2 /= 0)then
             rerror1=abs(new_term1/new_term2)
          end if
          if(rerror1 > eps)then
          else
             go to 10
          end if
       end do

10     continue

       volume_moment2=volume_moment2+new_term2
       rerror2=abs(new_term2/volume_moment2)
       if(rerror2 > eps)then
          new_term2=0.0_d_t
       else
          go to 11
       end if
    end do

11 end function volume_moment2

end module transport_kernel_module_CCE

