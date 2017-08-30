MODULE general_utility_module
  !***********************************************************************
  !
  ! Contains general utility functions used throughout THOR.
  ! This module should not depend on THOR modules except basic type definitions
  ! and the utilities provided should not manipulate global variables.
  !
  !***********************************************************************

  USE PRECISION
  USE types
  USE vector_types
  IMPLICIT NONE

CONTAINS

  SUBROUTINE cell_jacobian(v0,v1,v2,v3,J,det_J,J_inv)
    !*********************************************************************
    !
    ! Subroutine jacobian generates the Jacobian matrix and its inverse
    !
    !*********************************************************************
    ! Pass cell vertices

    TYPE(vector), INTENT(in) :: v0, v1, v2, v3

    ! Define temporary variables

    REAL(kind=d_t), INTENT(out) :: det_J
    REAL(kind=d_t), DIMENSION(3,3), INTENT(out) :: J, J_inv

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

    det_J=ABS(det_J)

  END SUBROUTINE cell_jacobian

  SUBROUTINE invert_face_jacobian(Jf,Jf_inv)
    !*********************************************************************
    !
    ! Subroutine invert face jacobian invert 3x2 matrix into 2x3 inverse
    !
    !*********************************************************************

    ! Define temporary variables
    !  integer(kind=li) :: q,l
    REAL(kind=d_t) :: det_M
    REAL(kind=d_t), DIMENSION(3,2), INTENT(in) :: Jf
    REAL(kind=d_t), DIMENSION(2,3), INTENT(out) :: Jf_inv
    REAL(kind=d_t), DIMENSION(2,3) :: JfT
    REAL(kind=d_t), DIMENSION(2,2) :: M, M_inv
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

  END SUBROUTINE invert_face_jacobian

  SUBROUTINE subcell_jacobian(v,Js,det_Js,Js_inv)
    !*********************************************************************
    !
    ! Subroutine subcell jacobian generates the subcell Jacobian
    ! matrix and its inverse
    !
    !*********************************************************************

    ! Define temporary variables

    INTEGER(kind=li) :: i, j
    REAL(kind=d_t), INTENT(out) :: det_Js
    !  real(kind=d_t), dimension(3,3) :: Identity
    TYPE(vector), DIMENSION(0:3), INTENT(in) :: v
    REAL(kind=d_t), DIMENSION(3,3),INTENT(out) :: Js, Js_inv

    ! Construct Jacobian matrix

    DO j=1, 3
      Js(1,j)=v(j)%x1-v(j-1)%x1
      Js(2,j)=v(j)%x2-v(j-1)%x2
      Js(3,j)=v(j)%x3-v(j-1)%x3
    END DO

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

    det_Js=ABS(det_Js)

  END SUBROUTINE subcell_jacobian

  SUBROUTINE lu_decomposition(n,A,L,U)
    !**********************************************************************
    !
    ! Subroutine lu decomposition decomposes L and U matrices from A
    !
    !**********************************************************************
    ! Pass input parameters

    INTEGER(kind=li), INTENT(in) :: n
    REAL(kind=d_t), DIMENSION(n,n), INTENT(in) :: A
    REAL(kind=d_t), DIMENSION(n,n), INTENT(out) :: &
          L, U

    ! Declare spatial order moments index

    INTEGER(kind=li) :: i, j, k

    DO j=1, n
      U(1,j)=A(1,j)
    END DO

    DO i=1, n
      k=i
      L(i,k)=1.0_d_t
    END DO

    DO k=1, n
      DO j=k, n
        U(k,j)=A(k,j)
        DO i=1, k-1
          U(k,j)=U(k,j)-L(k,i)*U(i,j)
        END DO
      END DO
      DO i=k+1, n
        L(i,k)=A(i,k)/U(k,k)
        DO j=1, k-1
          L(i,k)=L(i,k)-(L(i,j)*U(j,k))/U(k,k)
        END DO
      END DO
    END DO

  END SUBROUTINE lu_decomposition

  SUBROUTINE back_substitution(n,b,L,U,x)
    !**********************************************************************
    !
    ! Subroutine back substitution solves for unknown after LU
    !
    !**********************************************************************
    ! Pass input parameters

    INTEGER(kind=li), INTENT(in) :: n
    REAL(kind=d_t), DIMENSION(n), INTENT(in) :: b
    REAL(kind=d_t), DIMENSION(n,n), INTENT(in) :: L, U
    REAL(kind=d_t), DIMENSION(n), INTENT(out) :: x

    ! Declare spatial order moments index

    INTEGER(kind=li) :: i, j
    REAL(kind=d_t), DIMENSION(n) :: y

    DO i=1, n
      y(i)=b(i)
      DO j=1, i-1
        y(i)=y(i)-L(i,j)*y(j)
      END DO
    END DO

    DO i=n, 1, -1
      x(i)=y(i)/U(i,i)
      DO j=i+1, n
        x(i)=x(i)-U(i,j)*x(j)/U(i,i)
      END DO
    END DO

  END SUBROUTINE back_substitution

  SUBROUTINE linearSolve(n, mat, rhs, sol)

    ! pass variables
    INTEGER(kind=li), INTENT(in) :: n
    REAL(kind=d_t), INTENT(in) :: mat(n,n), rhs(n)
    REAL(kind=d_t), INTENT(inout) :: sol(n)

    ! local variables
    REAL(kind=8) :: A(n,n), b(n, 1)
    INTEGER :: info, ipiv(n)

    ! copy over some data
    A = mat
    b(:, 1) = rhs

    ! call lapack linear solver package
    CALL dgesv(n, 1, A, n, ipiv, b, n, info)

    sol = b(:, 1)

  END SUBROUTINE linearSolve

  SUBROUTINE barycentricCoordinates(p, v0, v1, v2, v3, lambda)

    ! pass arguments
    TYPE(vector), INTENT(in) :: p, v0, v1, v2, v3
    TYPE(vector), INTENT(inout) :: lambda

    ! local arguments to make math interface simpler
    REAL(kind=d_t) :: matrix(3,3), rhs(3), x(3)

    ! set the local variables
    rhs(1) = p%x1 - v3%x1
    rhs(2) = p%x2 - v3%x2
    rhs(3) = p%x3 - v3%x3

    matrix(1, 1) = v0%x1 - v3%x1
    matrix(2, 1) = v0%x2 - v3%x2
    matrix(3, 1) = v0%x3 - v3%x3

    matrix(1, 2) = v1%x1 - v3%x1
    matrix(2, 2) = v1%x2 - v3%x2
    matrix(3, 2) = v1%x3 - v3%x3

    matrix(1, 3) = v2%x1 - v3%x1
    matrix(2, 3) = v2%x2 - v3%x2
    matrix(3, 3) = v2%x3 - v3%x3

    CALL linearSolve(3, matrix, rhs, x)

    lambda%x1 = x(1)
    lambda%x2 = x(2)
    lambda%x3 = x(3)
  END SUBROUTINE barycentricCoordinates

END MODULE general_utility_module
