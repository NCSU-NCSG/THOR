module ahotc_matrix_module
!***********************************************************************
!
! Ahotc_matrix module generates and performs LU decomposition on the 
! mass matrices
!
!***********************************************************************

! User derived-type modules

  use types
  use parameter_types
  use multindex_types

  implicit none

contains

  subroutine ahotc_matrix(numv,numf,indv,indf,tM,tLL,tU,tMf,tLf,tUf)
  !*********************************************************************
  !
  ! Subroutine ahotc_matrix computes M, LL, and U
  !
  !*********************************************************************
  ! Declare index size

    integer(kind=li), intent(in) :: numv, numf

  ! Declare spatial order moments index

    type(indices_v), dimension(numv), intent(in) :: indv
    type(indices_f), dimension(numf), intent(in) :: indf

  ! Define temporary variables

    real(kind=d_t), dimension(numv,numv), intent(out) :: tM, tLL, tU
    real(kind=d_t), dimension(numf,numf), intent(out) :: tMf,tLf, tUf

    call mass_matrix(numv,numf,indv,indf,tM,tMf)

    call lu_decomposition(numv,tM,tLL,tU)

    call lu_decomposition(numf,tMf,tLf,tUf)

  end subroutine ahotc_matrix

  subroutine mass_matrix(numv,numf,indv,indf,tM,tMf)
  !**********************************************************************
  !
  ! Subroutine mass matrix generates mass matrix
  !
  !**********************************************************************
  ! Pass input parameters

    integer(kind=li), intent(in) :: numv, numf
    real(kind=d_t), dimension(numv,numv),intent(out) :: tM
    real(kind=d_t), dimension(numf,numf),intent(out) :: tMf

  ! Declare spatial order moments index
    
    integer(kind=li) :: q, l
    type(indices_v), dimension(numv), intent(in) :: indv
    type(indices_f), dimension(numf), intent(in) :: indf

    do q=1, numv
       do l=1, numv
          tM(q,l)=(6.0_d_t)/((indv(q)%i1+indv(l)%i1+indv(q)%i2+&
               indv(l)%i2+indv(q)%i3+indv(l)%i3+3.0_d_t)*&
               (indv(q)%i2+indv(l)%i2+indv(q)%i3+indv(l)%i3+&
               2.0_d_t)*(indv(q)%i3+indv(l)%i3+1.0_d_t))
       end do
    end do
    
    do q=1, numf
       do l=1, numf
          tMf(q,l)=(2.0_d_t)/&
               ((indf(q)%i1+indf(l)%i1+indf(q)%i2+&
               indf(l)%i2+2.0_d_t)*(indf(q)%i2+indf(l)%i2+&
               1.0_d_t))
       end do
    end do

  end subroutine mass_matrix

  subroutine lu_decomposition(n,A,L,U)
  !**********************************************************************
  !
  ! Subroutine lu decomposition decomposes L and U matrices from A
  !
  !**********************************************************************
  ! Pass input parameters

    integer(kind=li), intent(in) :: n
    real(kind=d_t), dimension(n,n), intent(in)  :: A
    real(kind=d_t), dimension(n,n), intent(out) :: L, U

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

end module ahotc_matrix_module

