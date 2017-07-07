module cell_splitting_module
!***********************************************************************
!
! Cell splitting module divides each tetrahedral cell into subcells
! along the discrete ordinate of interest. Six possible configurations
! account for all the possible combinations of incoming/outgoing faces.
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

! Use modules that pertain setting up problem
  ! FIXME: You need to find a better place to store cell_jacobian
  !        it's defined in all transport kernels...and then used in cell_splitting
  !
  use general_utility_module  
  use transport_kernel_module_SC, only: transport_kernel_SC
  use transport_kernel_module_LC, only: transport_kernel_LC
  use transport_kernel_module_CCE, only: transport_kernel_CCE
  use termination_module

  implicit none

contains

  !> THOR is based in the AHOT-C methods. The AHOT-C method splits tetrahedra
  !> into characteristic tetrahedra and solves the characteristic equations
  !> on each of these. This is the entry point for the single mesh-cell solver.
  subroutine cell_splitting(sigmat,qm_moments,vol_flux,face_flux,octant,  &
                            LL,U,Lf,Uf,omega,i,face_known,case,n0,n1,n2,n3) 
  !*********************************************************************
  !
  ! Subroutine transport kernel divides tetrahedra into 'CTs' and 
  ! calculates outgoing face angular fluxes based on incoming ones
  !
  !*********************************************************************

  ! Pass arguments

    real(kind=d_t), intent(in)                             :: sigmat
    real(kind=d_t), dimension(num_moments_v), intent(in)   :: qm_moments
    real(kind=d_t), dimension(num_moments_v,num_cells)     :: vol_flux
    real(kind=d_t), dimension(num_moments_f,0:3,num_cells) :: face_flux
    integer(kind=li), intent(in)                           :: octant, i, case
    real(kind=d_t), dimension(num_moments_v,num_moments_v) :: LL,U
    real(kind=d_t), dimension(num_moments_f,num_moments_f) :: Lf,Uf
    type(vector), intent(in)                               :: omega, n0, n1, n2,n3
    integer(kind=1), dimension(0:3,num_cells)              :: face_known

  ! Define temporary variables

    integer(kind=li)                            :: alloc_stat, subcells
    integer(kind=li), dimension(:), allocatable :: incoming_face,outgoing_face
    integer(kind=li), dimension(:), allocatable :: bc
    real(kind=d_t)                              :: t, det_J
    real(kind=d_t), dimension(3,3)              :: J, J_inv
    real(kind=d_t), dimension(:,:), allocatable :: upstream_moments
    integer(kind=1)                             :: face0, face1, face2, face3
    type(vector)                                :: omega_local, v0, v1, v2, v3
    type(vector), dimension(:), allocatable     :: r0, r1, r2, r3

  ! Call Jacobian subroutine to generate cell transformation

    v0=vertices(cells(i)%R(0))%v
    v1=vertices(cells(i)%R(1))%v
    v2=vertices(cells(i)%R(2))%v
    v3=vertices(cells(i)%R(3))%v

    call cell_jacobian(v0,v1,v2,v3,J,det_J,J_inv)

  ! Compute local coordinate omega_local

    omega_local=J_inv*omega

  ! Split cell into subcells depending on cases

    face0=face_known(0,i)
    face1=face_known(1,i)
    face2=face_known(2,i)
    face3=face_known(3,i)

  ! CASE 1: 3 incoming, 1 outgoing
    if(case == 1)then
       subcells=3
       
       allocate(r0(subcells),r1(subcells),r2(subcells),&
            r3(subcells),incoming_face(subcells),&
            outgoing_face(subcells),upstream_moments(num_moments_f,subcells),&
            bc(subcells),stat=alloc_stat)
       if(alloc_stat /= 0) call stop_thor(2_li)
      
       call case1(v0,v1,v2,v3,omega,r0,r1,r2,r3,face0,face1,face2,&
            face3,incoming_face,outgoing_face,omega_local,t)

  ! CASE 2: 1 incoming, 3 outgoing
    elseif(case == 2)then
       subcells=3

       allocate(r0(subcells),r1(subcells),r2(subcells),&
            r3(subcells),incoming_face(subcells),&
            outgoing_face(subcells),upstream_moments(num_moments_f,subcells),&
            bc(subcells),stat=alloc_stat)
       if(alloc_stat /= 0) call stop_thor(2_li)

       call case2(v0,v1,v2,v3,omega,r0,r1,r2,r3,face0,face1,face2,&
            face3,incoming_face,outgoing_face,omega_local,t)

  ! CASE 3: 2 incoming, 2 outgoing
    elseif(case == 3)then
       subcells=4
       
       allocate(r0(subcells),r1(subcells),r2(subcells),&
            r3(subcells),incoming_face(subcells),&
            outgoing_face(subcells),upstream_moments(num_moments_f,subcells),&
            bc(subcells),stat=alloc_stat)
       if(alloc_stat /= 0) call stop_thor(2_li)
       
       call case3(v0,v1,v2,v3,r0,r1,r2,r3,face0,face1,face2,face3,&
            incoming_face,outgoing_face,J,omega_local,t)

  ! CASE 4: 2 incoming, 1 outgoing
    elseif(case == 4)then
       subcells=2

       allocate(r0(subcells),r1(subcells),r2(subcells),&
            r3(subcells),incoming_face(subcells),&
            outgoing_face(subcells),upstream_moments(num_moments_f,subcells),&
            bc(subcells),stat=alloc_stat)
       if(alloc_stat /= 0) call stop_thor(2_li)

       call case4(v0,v1,v2,v3,omega,n0,n1,n2,n3,r0,r1,r2,r3,face0,&
            face1,face2,face3,incoming_face,outgoing_face,omega_local,t)

  ! CASE 5: 1 incoming, 2 outgoing
    elseif(case == 5)then
       subcells=2

       allocate(r0(subcells),r1(subcells),r2(subcells),&
            r3(subcells),incoming_face(subcells),&
            outgoing_face(subcells),upstream_moments(num_moments_f,subcells),&
            bc(subcells),stat=alloc_stat)
       if(alloc_stat /= 0) call stop_thor(2_li)

       call case5(v0,v1,v2,v3,omega,n0,n1,n2,n3,r0,r1,r2,r3,face0,&
            face1,face2,face3,incoming_face,outgoing_face,omega_local,t)
       
  ! CASE 6: 1 incoming, 1 outgoing
    elseif(case == 6)then
       subcells=1
        
       allocate(r0(subcells),r1(subcells),r2(subcells),&
            r3(subcells),incoming_face(subcells),&
            outgoing_face(subcells),upstream_moments(num_moments_f,subcells),&
            bc(subcells),stat=alloc_stat)
       if(alloc_stat /= 0) call stop_thor(2_li)

       call case6(v0,v1,v2,v3,omega,n0,n1,n2,n3,r0,r1,r2,r3,face0,&
            face1,face2,face3,incoming_face,outgoing_face,omega_local,t)

  ! CASE UNKNOWN: failure
    else
       write(6,*) "Unacceptable case from cell splitting in cell", i
       call stop_thor(-1_li)
    end if

    call upstream_mom(num_cells,adjacent_cells,num_moments_f,     &
         adjacency_list,face_flux,octant,i,subcells,incoming_face,&
         n0,n1,n2,n3,upstream_moments)

    if (space_ord == 0) then
      call transport_kernel_SC(sigmat,qm_moments,vol_flux(:,i),                &
                               face_flux(:,:,i),LL,U,Lf,Uf,i,                  &
                               t,subcells,r0,r1,r2,r3,face_known,incoming_face,&
                               outgoing_face,J,J_inv,upstream_moments)
    else if (space_ord == -1) then
      call transport_kernel_LC(sigmat,qm_moments,vol_flux(:,i),                &
                               face_flux(:,:,i),LL,U,Lf,Uf,i,                  &
                               t,subcells,r0,r1,r2,r3,face_known,incoming_face,&
                               outgoing_face,J,J_inv,upstream_moments)
    else
      call transport_kernel_CCE(sigmat,qm_moments,vol_flux(:,i),                &
                                face_flux(:,:,i),LL,U,Lf,Uf,i,                  &
                                t,subcells,r0,r1,r2,r3,face_known,incoming_face,&
                                outgoing_face,J,J_inv,upstream_moments)
    
    end if

    deallocate(r0,r1,r2,r3,upstream_moments)

  end subroutine cell_splitting

  subroutine upstream_mom(nc,nadj,numf,adj_list,face_flux,octant,i,&
                          subcells,incoming_face,n0,n1,n2,n3,upstream_moments)
  !*********************************************************************
  !
  ! Subroutine upstream_mom determines upstream moments
  !
  !*********************************************************************
  ! Pass geometry parameters

    integer(kind=li), intent(in) :: nc, nadj

  ! Declare index size

    integer(kind=li), intent(in) :: numf
    
  ! Pass geometry derived type
    
    type(list), dimension(nadj,0:3), intent(in) :: adj_list

  ! Declare angular and scalar flux types used globally

    real(kind=d_t) :: face_flux(num_moments_f,0:3,num_cells)

  ! Define temporary variables

    integer(kind=li)                                        :: ii, up_cell, up_face, l, qoct
    integer(kind=li), intent(in)                            :: octant, i, subcells
    integer(kind=li), dimension(subcells), intent(in)       :: incoming_face
    real(kind=d_t), dimension(numf,subcells), intent(inout) :: upstream_moments
    type(vector), intent(in)                                :: n0, n1, n2, n3

  ! Assign upstream values and boundary conditions


    do ii=1, subcells
       up_cell=adj_list(i,incoming_face(ii))%cell
       up_face=adj_list(i,incoming_face(ii))%face
       if(adj_list(i,incoming_face(ii))%cell .ne. 0 .and. is_cycle(incoming_face(ii),i).eq.0 )then
         do l=1, numf
           upstream_moments(l,ii)=face_flux(l,up_face,up_cell)
         end do
       else 
         do l=1, numf
           upstream_moments(l,ii)=face_flux(l,incoming_face(ii),i)
         end do 
      end if
    end do

  end subroutine upstream_mom

  subroutine case1(v0,v1,v2,v3,omega,r0,r1,r2,r3,face0,face1,face2,&
       face3,incoming_face,outgoing_face,omega_local,t)
  !*********************************************************************
  !
  ! Subroutine case1 generates local coordinate subcells vectors
  !
  !*********************************************************************
  ! Define temporary variables

    integer(kind=li), dimension(3), intent(inout) :: incoming_face,&
         outgoing_face
    real(kind=d_t), intent(out) :: t
    integer(kind=1), intent(in) :: face0,face1,face2,face3
    type(vector) :: Rin, Rout
    type(vector), intent(in) :: omega, omega_local, v0, v1, v2, v3
    type(vector), dimension(3),intent(inout) :: r0, r1, r2, r3

  ! Split cells and define subcells vertices r0, r1, r2, & r3

    if(face0 .eq. 0)then
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
    elseif(face1 .eq. 0)then
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
    elseif(face2 .eq. 0)then
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
    elseif(face3 .eq. 0)then
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
    else
       write(6,*) "Cell splitting case 1 failed for cell"
    end if

  end subroutine case1

  subroutine case2(v0,v1,v2,v3,omega,r0,r1,r2,r3,face0,face1,face2,&
       face3,incoming_face,outgoing_face,omega_local,t)
  !*********************************************************************
  !
  ! Subroutine case2 generates local coordinate subcell vectors
  !
  !*********************************************************************
  ! Define temporary variables

    integer(kind=li), dimension(3), intent(inout) :: incoming_face,&
         outgoing_face
    real(kind=d_t), intent(out) :: t
    integer(kind=1), intent(in) :: face0,face1,face2,face3
    type(vector) :: Rin, Rout
    type(vector), intent(in) :: omega, omega_local, v0, v1, v2, v3
    type(vector), dimension(3),intent(inout) :: r0, r1, r2, r3

  ! Split cell and define subcell vertices r0, r1, r2, & r3

    if(face0 .eq. 1)then
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
    elseif(face1 .eq. 1)then
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
    elseif(face2 .eq. 1)then
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
    elseif(face3 .eq. 1)then
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
    else
       write(6,*) "Cell splitting case 2 failed for cell"
    end if

  end subroutine case2

  subroutine case3(v0,v1,v2,v3,r0,r1,r2,r3,face0,face1,face2,face3,&
       incoming_face,outgoing_face,J,omega_local,t)
  !*********************************************************************
  !
  ! Subroutine case3 generates local coordinate subcell vectors
  !
  !*********************************************************************
  ! Define temporary variables

    integer(kind=li), dimension(4), intent(inout) :: incoming_face,&
         outgoing_face
    real(kind=d_t), intent(out) :: t
    real(kind=d_t), dimension(3,3), intent(in) :: J
    integer(kind=1), intent(in) :: face0,face1,face2,face3
    type(vector) :: Rin, Rout, Rin_local, Rout_local
    type(vector), intent(in) :: omega_local, v0, v1, v2, v3
    type(vector), dimension(4),intent(inout) :: r0, r1, r2, r3

  ! Split cell and define subcell vertices r0, r1, r2, & r3

    if((face2 .eq. 1) .and. (face3 .eq. 1))then
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
    elseif((face1 .eq. 1) .and. (face3 .eq. 1))then
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
    elseif((face1 .eq. 1) .and. (face2 .eq. 1))then
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
    elseif((face0 .eq. 1) .and. (face3 .eq. 1))then
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
    elseif((face0 .eq. 1) .and. (face2 .eq. 1))then
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
    elseif((face0 .eq. 1) .and. (face1 .eq. 1))then
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
    else
       write(6,*) "Cell splitting case 3 failed for cell"
    end if

  end subroutine case3

  subroutine case4(v0,v1,v2,v3,omega,n0,n1,n2,n3,r0,r1,r2,r3,face0,&
       face1,face2,face3,incoming_face,outgoing_face,omega_local,t)
  !*********************************************************************
  !
  ! Subroutine case4 generates local coordinate subcell vectors
  !
  !*********************************************************************
  ! Define temporary variables

    integer(kind=li), dimension(2), intent(inout) :: incoming_face,&
         outgoing_face
    real(kind=d_t), intent(out) :: t
    integer(kind=1), intent(in) :: face0,face1,face2,face3
    type(vector) :: Rin, Rout
    type(vector), intent(in) :: omega, omega_local, v0, v1, v2, v3, &
         n0,n1,n2,n3
    type(vector), dimension(2),intent(inout) :: r0, r1, r2, r3

  ! Split cell and define subcell vertices r0, r1, r2, & r3

    if((face0 .eq. 0) .and. (omega .dot. n0) /= 0.0)then
       t=1.0_d_t/(omega_local%x1)
       Rin=v0
       Rout=Rin+t*omega
       if((face2 .eq. 1) .and. (face3 .eq. 1))then
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
       elseif((face1 .eq. 1) .and. (face3 .eq. 1))then
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
       elseif((face1 .eq. 1) .and. (face2 .eq. 1))then
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
       endif
    elseif((face1 .eq. 0) .and. (omega .dot. n1) /= zero)then
       t=1.0_d_t/(omega_local%x2-omega_local%x1)
       Rin=v1
       Rout=Rin+t*omega
       if((face2 .eq. 1) .and. (face3 .eq. 1))then
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
       elseif((face0 .eq. 1) .and. (face3 .eq. 1))then
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
       elseif((face0 .eq. 1) .and. (face2 .eq. 1))then
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
       endif
    elseif((face2 .eq. 0) .and. (omega .dot. n2) /= 0.0)then
       t=1.0_d_t/(omega_local%x3-omega_local%x2)
       Rin=v2
       Rout=Rin+t*omega
       if((face1 .eq. 1) .and. (face3 .eq. 1))then
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
       elseif((face0 .eq. 1) .and. (face3 .eq. 1))then
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
       elseif((face0 .eq. 1) .and. (face1 .eq. 1))then
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
       endif
    elseif((face3 .eq. 0) .and. (omega .dot. n3) /= 0.0)then
       t=1.0_d_t/(-1.0_d_t*omega_local%x3)
       Rin=v3
       Rout=Rin+t*omega
       if((face1 .eq. 1) .and. (face2 .eq. 1))then
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
       elseif((face0 .eq. 1) .and. (face2 .eq. 1))then
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
       elseif((face0 .eq. 1) .and. (face1 .eq. 1))then
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
       endif
    else
       write(6,*) "Cell splitting case 4 failed for cell"
    end if

  end subroutine case4

  subroutine case5(v0,v1,v2,v3,omega,n0,n1,n2,n3,r0,r1,r2,r3,face0,&
       face1,face2,face3,incoming_face,outgoing_face,omega_local,t)
  !*********************************************************************
  !
  ! Subroutine case5 generates local coordinate subcell vectors
  !
  !*********************************************************************
  ! Define temporary variables
    integer(kind=li), dimension(2), intent(inout) :: incoming_face,&
         outgoing_face
    real(kind=d_t), intent(out) :: t
    integer(kind=1), intent(in) :: face0,face1,face2,face3
    type(vector) :: Rin, Rout
    type(vector), intent(in) :: omega, omega_local, v0, v1, v2, v3, &
         n0,n1,n2,n3
    type(vector), dimension(2),intent(inout) :: r0, r1, r2, r3

  ! Split cell and define subcell vertices r0, r1, r2, & r3

    if(face0 .eq. 1)then
       t=1.0_d_t/(-1.0_d_t*omega_local%x1)
       Rout=v0
       Rin=Rout-t*omega
       if((omega .dot. n1) == 0.0)then
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
       elseif((omega .dot. n2) == 0.0)then
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
       elseif((omega .dot. n3) == 0.0)then
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
       endif
    elseif(face1 .eq. 1)then
       t=1.0_d_t/(omega_local%x1-omega_local%x2)
       Rout=v1
       Rin=Rout-t*omega
       if((omega .dot. n0) == 0.0)then
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
       elseif((omega .dot. n2) == 0.0)then
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
       elseif((omega .dot. n3) == 0.0)then
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
       endif
    elseif(face2 .eq. 1)then
       t=1.0_d_t/(omega_local%x2-omega_local%x3)
       Rout=v2
       Rin=Rout-t*omega
       if((omega .dot. n0) == 0.0)then
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
       elseif((omega .dot. n1) == 0.0)then
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
       elseif((omega .dot. n3) == 0.0)then
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
       endif
    elseif(face3 .eq. 1)then
       t=1.0_d_t/(omega_local%x3)
       Rout=v3
       Rin=Rout-t*omega
       if((omega .dot. n0) == 0.0)then
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
       elseif((omega .dot. n1) == 0.0)then
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
       elseif((omega .dot. n2) == 0.0)then
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
       endif
    else
       write(6,*) "Cell splitting case 5 failed for cell"
    end if

  end subroutine case5

  subroutine case6(v0,v1,v2,v3,omega,n0,n1,n2,n3,r0,r1,r2,r3,face0,&
       face1,face2,face3,incoming_face,outgoing_face,omega_local,t)
  !*********************************************************************
  !
  ! Subroutine case6 generates local coordinate subcell vectors
  !
  !*********************************************************************
  ! Define temporary variables

    integer(kind=li), dimension(1), intent(inout) :: incoming_face,&
         outgoing_face
    real(kind=d_t), intent(out) :: t
    integer(kind=1), intent(in) :: face0,face1,face2,face3
    type(vector) :: Rin, Rout
    type(vector), intent(in) :: omega, omega_local, v0, v1, v2, v3, &
         n0,n1,n2,n3
    type(vector), dimension(1), intent(inout) :: r0, r1, r2, r3

  ! Split cell and define subcell vertices r0, r1, r2, & r3

    if((face0 .eq. 0) .and. (omega .dot. n0) /= 0.0)then
       t=1.0_d_t/(omega_local%x1)
       Rin=v0
       Rout=Rin+t*omega
       if(face1 .eq. 1)then
      ! Subcell 1
          incoming_face(1)=1
          outgoing_face(1)=0
          r0(1)=v2
          r1(1)=v3
          r2(1)=Rin
          r3(1)=Rout
       elseif(face2 .eq. 1)then
      ! Subcell 1
          incoming_face(1)=2
          outgoing_face(1)=0
          r0(1)=v3
          r1(1)=v1
          r2(1)=Rin
          r3(1)=Rout
       elseif(face3 .eq. 1)then
      ! Subcell 1
          incoming_face(1)=3
          outgoing_face(1)=0
          r0(1)=v1
          r1(1)=v2
          r2(1)=Rin
          r3(1)=Rout
       endif
    elseif((face1 .eq. 0) .and. (omega .dot. n1) /= 0.0)then
       t=1.0_d_t/(omega_local%x2-omega_local%x1)
       Rin=v1
       Rout=Rin+t*omega
       if(face0 .eq. 1)then
      ! Subcell 1
          incoming_face(1)=0
          outgoing_face(1)=1
          r0(1)=v3
          r1(1)=v2
          r2(1)=Rin
          r3(1)=Rout
       elseif(face2 .eq. 1)then
      ! Subcell 1
          incoming_face(1)=2
          outgoing_face(1)=1
          r0(1)=v0
          r1(1)=v3
          r2(1)=Rin
          r3(1)=Rout
       elseif(face3 .eq. 1)then
      ! Subcell 1
          incoming_face(1)=3
          outgoing_face(1)=1
          r0(1)=v2
          r1(1)=v0
          r2(1)=Rin
          r3(1)=Rout
       endif
    elseif((face2 .eq. 0) .and. (omega .dot. n2) /= 0.0)then
       t=1.0_d_t/(omega_local%x3-omega_local%x2)
       Rin=v2
       Rout=Rin+t*omega
       if(face0 .eq. 1)then
      ! Subcell 1
          incoming_face(1)=0
          outgoing_face(1)=2
          r0(1)=v1
          r1(1)=v3
          r2(1)=Rin
          r3(1)=Rout
       elseif(face1 .eq. 1)then
      ! Subcell 1
          incoming_face(1)=1
          outgoing_face(1)=2
          r0(1)=v3
          r1(1)=v0
          r2(1)=Rin
          r3(1)=Rout
       elseif(face3 .eq. 1)then
      ! Subcell 1
          incoming_face(1)=3
          outgoing_face(1)=2
          r0(1)=v0
          r1(1)=v1
          r2(1)=Rin
          r3(1)=Rout
       endif
    elseif((face3 .eq. 0) .and. (omega .dot. n3) /= 0.0)then
       t=1.0_d_t/(-1.0_d_t*omega_local%x3)
       Rin=v3
       Rout=Rin+t*omega
       if(face0 .eq. 1)then
      ! Subcell 1
          incoming_face(1)=0
          outgoing_face(1)=3
          r0(1)=v2
          r1(1)=v1
          r2(1)=Rin
          r3(1)=Rout
       elseif(face1 .eq. 1)then
      ! Subcell 1
          incoming_face(1)=1
          outgoing_face(1)=3
          r0(1)=v0
          r1(1)=v2
          r2(1)=Rin
          r3(1)=Rout
       elseif(face2 .eq. 1)then
      ! Subcell 1
          incoming_face(1)=2
          outgoing_face(1)=3
          r0(1)=v1
          r1(1)=v0
          r2(1)=Rin
          r3(1)=Rout
       endif
    else
       write(6,*) "Cell splitting case 6 failed for cell"
    end if

  end subroutine case6

end module cell_splitting_module

