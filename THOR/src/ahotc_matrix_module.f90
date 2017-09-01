MODULE ahotc_matrix_module
  !***********************************************************************
  !
  ! Ahotc_matrix module generates and performs LU decomposition on the
  ! mass matrices
  !
  !***********************************************************************

  ! User derived-type modules

  USE types
  USE parameter_types
  USE multindex_types
  USE general_utility_module

  IMPLICIT NONE

CONTAINS

  SUBROUTINE ahotc_matrix(numv,numf,indv,indf,tM,tLL,tU,tMf,tLf,tUf)
    !*********************************************************************
    !
    ! Subroutine ahotc_matrix computes M, LL, and U
    !
    !*********************************************************************
    ! Declare index size

    INTEGER(kind=li), INTENT(in) :: numv, numf

    ! Declare spatial order moments index

    TYPE(indices_v), DIMENSION(numv), INTENT(in) :: indv
    TYPE(indices_f), DIMENSION(numf), INTENT(in) :: indf

    ! Define temporary variables

    REAL(kind=d_t), DIMENSION(numv,numv), INTENT(out) :: tM, tLL, tU
    REAL(kind=d_t), DIMENSION(numf,numf), INTENT(out) :: tMf,tLf, tUf

    CALL mass_matrix(numv,numf,indv,indf,tM,tMf)

    CALL lu_decomposition(numv,tM,tLL,tU)

    CALL lu_decomposition(numf,tMf,tLf,tUf)

  END SUBROUTINE ahotc_matrix

  SUBROUTINE mass_matrix(numv,numf,indv,indf,tM,tMf)
    !**********************************************************************
    !
    ! Subroutine mass matrix generates mass matrix
    !
    !**********************************************************************
    ! Pass input parameters

    INTEGER(kind=li), INTENT(in) :: numv, numf
    REAL(kind=d_t), DIMENSION(numv,numv),INTENT(out) :: tM
    REAL(kind=d_t), DIMENSION(numf,numf),INTENT(out) :: tMf

    ! Declare spatial order moments index

    INTEGER(kind=li) :: q, l
    TYPE(indices_v), DIMENSION(numv), INTENT(in) :: indv
    TYPE(indices_f), DIMENSION(numf), INTENT(in) :: indf

    DO q=1, numv
      DO l=1, numv
        tM(q,l)=(6.0_d_t)/((indv(q)%i1+indv(l)%i1+indv(q)%i2+&
              indv(l)%i2+indv(q)%i3+indv(l)%i3+3.0_d_t)*&
              (indv(q)%i2+indv(l)%i2+indv(q)%i3+indv(l)%i3+&
              2.0_d_t)*(indv(q)%i3+indv(l)%i3+1.0_d_t))
      END DO
    END DO

    DO q=1, numf
      DO l=1, numf
        tMf(q,l)=(2.0_d_t)/&
              ((indf(q)%i1+indf(l)%i1+indf(q)%i2+&
              indf(l)%i2+2.0_d_t)*(indf(q)%i2+indf(l)%i2+&
              1.0_d_t))
      END DO
    END DO

  END SUBROUTINE mass_matrix

END MODULE ahotc_matrix_module
