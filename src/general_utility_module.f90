module general_utility_module
!***********************************************************************
!
! Contains general utility functions used throughout THOR. 
! This module should not depend on THOR modules except basic type definitions
! and the utilities provided should not manipulate global variables.
!
!***********************************************************************

  use precision
  use types
  use vector_types 
  implicit none

contains

  subroutine cell_jacobian(v0,v1,v2,v3,J,det_J,J_inv)
  !*********************************************************************
  !
  ! Subroutine jacobian generates the Jacobian matrix and its inverse
  !
  !*********************************************************************
  ! Pass cell vertices

    type(vector), intent(in) :: v0, v1, v2, v3

  ! Define temporary variables

    real(kind=d_t), intent(out) :: det_J
    real(kind=d_t), dimension(3,3), intent(out) :: J, J_inv

  ! Construct Jacobian matrix

    J(1,1)=v1%x1-v0%x1
    J(2,1)=v1%x2-v0%x2
    J(3,1)=v1%x3-v0%x3
    
    J(1,2)=v2%x1-v1%x1
    J(2,2)=v2%x2-v1%x2
    J(3,2)=v2%x3-v1%x3
    
    J(1,3)=v3%x1-v2%x1
    J(2,3)=v3%x2-v2%x2
    J(3,3)=v3%x3-v2%x3

  ! Compute determinant and check if matrix is nonsingular

    det_J=     J(1,1)*J(2,2)*J(3,3)-J(1,1)*J(2,3)*J(3,2)+&
               J(1,2)*J(2,3)*J(3,1)-J(1,2)*J(2,1)*J(3,3)+&
               J(1,3)*J(2,1)*J(3,2)-J(1,3)*J(2,2)*J(3,1) 

  ! Construct inverse Jacobian matrix

    J_inv(1,1)=J(2,2)*J(3,3)*(1.0_d_t/det_J)-&
         J(2,3)*J(3,2)*(1.0_d_t/det_J)
    J_inv(1,2)=J(1,3)*J(3,2)*(1.0_d_t/det_J)-&
         J(1,2)*J(3,3)*(1.0_d_t/det_J)
    J_inv(1,3)=J(1,2)*J(2,3)*(1.0_d_t/det_J)-&
         J(1,3)*J(2,2)*(1.0_d_t/det_J)

    J_inv(2,1)=J(2,3)*J(3,1)*(1.0_d_t/det_J)-&
         J(2,1)*J(3,3)*(1.0_d_t/det_J)
    J_inv(2,2)=J(1,1)*J(3,3)*(1.0_d_t/det_J)-&
         J(1,3)*J(3,1)*(1.0_d_t/det_J)
    J_inv(2,3)=J(1,3)*J(2,1)*(1.0_d_t/det_J)-&
         J(1,1)*J(2,3)*(1.0_d_t/det_J)

    J_inv(3,1)=J(2,1)*J(3,2)*(1.0_d_t/det_J)-&
         J(2,2)*J(3,1)*(1.0_d_t/det_J)
    J_inv(3,2)=J(1,2)*J(3,1)*(1.0_d_t/det_J)-&
         J(1,1)*J(3,2)*(1.0_d_t/det_J)
    J_inv(3,3)=J(1,1)*J(2,2)*(1.0_d_t/det_J)-&
         J(1,2)*J(2,1)*(1.0_d_t/det_J)

    ! return abs(det_J)

      det_J=abs(det_J)

  end subroutine cell_jacobian

  subroutine invert_face_jacobian(Jf,Jf_inv)
  !*********************************************************************
  !
  ! Subroutine invert face jacobian invert 3x2 matrix into 2x3 inverse
  !
  !*********************************************************************

  ! Define temporary variables
  !  integer(kind=li) :: q,l
    real(kind=d_t) :: det_M
    real(kind=d_t), dimension(3,2), intent(in) :: Jf
    real(kind=d_t), dimension(2,3), intent(out) :: Jf_inv
    real(kind=d_t), dimension(2,3) :: JfT
    real(kind=d_t), dimension(2,2) :: M, M_inv
  !  real(kind=d_t), dimension(2,2) :: Identity

  ! Construct matrix transpose of Jf

    JfT(1,1)=Jf(1,1)
    JfT(1,2)=Jf(2,1)
    JfT(1,3)=Jf(3,1)
    JfT(2,1)=Jf(1,2)
    JfT(2,2)=Jf(2,2)
    JfT(2,3)=Jf(3,2)

  ! Construct temporary matrix A

    M=MATMUL(JfT,Jf)

  ! Invert temporary matrix A

    det_M=M(1,1)*M(2,2)-M(1,2)*M(2,1)

    M_inv(1,1)=M(2,2)*(1.0_d_t/det_M)
    M_inv(1,2)=-M(1,2)*(1.0_d_t/det_M)
    M_inv(2,1)=-M(2,1)*(1.0_d_t/det_M)
    M_inv(2,2)=M(1,1)*(1.0_d_t/det_M)

  ! Compute left pseudo-inverse
    
    Jf_inv=MATMUL(M_inv,JfT)

  end subroutine invert_face_jacobian

  subroutine subcell_jacobian(v,Js,det_Js,Js_inv)
  !*********************************************************************
  !
  ! Subroutine subcell jacobian generates the subcell Jacobian 
  ! matrix and its inverse
  !
  !*********************************************************************

  ! Define temporary variables

    integer(kind=li) :: i, j
    real(kind=d_t), intent(out) :: det_Js
  !  real(kind=d_t), dimension(3,3) :: Identity
    type(vector), dimension(0:3), intent(in) :: v
    real(kind=d_t), dimension(3,3),intent(out) :: Js, Js_inv

  ! Construct Jacobian matrix

    do j=1, 3
       Js(1,j)=v(j)%x1-v(j-1)%x1
       Js(2,j)=v(j)%x2-v(j-1)%x2
       Js(3,j)=v(j)%x3-v(j-1)%x3
    end do

  ! Compute determinant and check if matrix is nonsingular

    det_Js=     Js(1,1)*Js(2,2)*Js(3,3)-Js(1,1)*Js(2,3)*Js(3,2)+&
                Js(1,2)*Js(2,3)*Js(3,1)-Js(1,2)*Js(2,1)*Js(3,3)+&
                Js(1,3)*Js(2,1)*Js(3,2)-Js(1,3)*Js(2,2)*Js(3,1) 

  ! Construct inverse Jacobian matrix

    Js_inv(1,1)=Js(2,2)*Js(3,3)*(1.0_d_t/det_Js)-&
         Js(2,3)*Js(3,2)*(1.0_d_t/det_Js)
    Js_inv(1,2)=Js(1,3)*Js(3,2)*(1.0_d_t/det_Js)-&
         Js(1,2)*Js(3,3)*(1.0_d_t/det_Js)
    Js_inv(1,3)=Js(1,2)*Js(2,3)*(1.0_d_t/det_Js)-&
         Js(1,3)*Js(2,2)*(1.0_d_t/det_Js)

    Js_inv(2,1)=Js(2,3)*Js(3,1)*(1.0_d_t/det_Js)-&
         Js(2,1)*Js(3,3)*(1.0_d_t/det_Js)
    Js_inv(2,2)=Js(1,1)*Js(3,3)*(1.0_d_t/det_Js)-&
         Js(1,3)*Js(3,1)*(1.0_d_t/det_Js)
    Js_inv(2,3)=Js(1,3)*Js(2,1)*(1.0_d_t/det_Js)-&
         Js(1,1)*Js(2,3)*(1.0_d_t/det_Js)

    Js_inv(3,1)=Js(2,1)*Js(3,2)*(1.0_d_t/det_Js)-&
         Js(2,2)*Js(3,1)*(1.0_d_t/det_Js)
    Js_inv(3,2)=Js(1,2)*Js(3,1)*(1.0_d_t/det_Js)-&
         Js(1,1)*Js(3,2)*(1.0_d_t/det_Js)
    Js_inv(3,3)=Js(1,1)*Js(2,2)*(1.0_d_t/det_Js)-&
         Js(1,2)*Js(2,1)*(1.0_d_t/det_Js)

  ! return abs(det_Js)

    det_Js=abs(det_Js)

  end subroutine subcell_jacobian

  subroutine lu_decomposition(n,A,L,U)
  !**********************************************************************
  !
  ! Subroutine lu decomposition decomposes L and U matrices from A
  !
  !**********************************************************************
  ! Pass input parameters

    integer(kind=li), intent(in) :: n
    real(kind=d_t), dimension(n,n), intent(in) :: A
    real(kind=d_t), dimension(n,n), intent(out) :: &
         L, U

  ! Declare spatial order moments index

    integer(kind=li) :: i, j, k

    do j=1, n
       U(1,j)=A(1,j)
    end do

    do i=1, n
       k=i
       L(i,k)=1.0_d_t
    end do

    do k=1, n
       do j=k, n
          U(k,j)=A(k,j)
          do i=1, k-1
             U(k,j)=U(k,j)-L(k,i)*U(i,j)
          end do
       end do
       do i=k+1, n
          L(i,k)=A(i,k)/U(k,k)
          do j=1, k-1
             L(i,k)=L(i,k)-(L(i,j)*U(j,k))/U(k,k)
          end do
       end do
    end do

  end subroutine lu_decomposition

end module general_utility_module 

