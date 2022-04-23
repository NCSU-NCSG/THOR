MODULE quadrature_module
  !***********************************************************************
  !
  ! Quadrature module contains subroutines used to generate discrete
  ! ordinates quadrature. Currently only selects from pre-computed list
  !
  !***********************************************************************

  ! User derived-type modules

  USE types
  USE parameter_types
  USE vector_types
  USE angle_types
  USE global_variables
  USE error_module

  IMPLICIT NONE

CONTAINS

  SUBROUTINE quad_gen
    !*********************************************************************
    !
    ! Selects desired quadrature from list
    !
    !*********************************************************************

    INTEGER(kind=li) :: alloc_stat
    REAL(kind=d_t) :: wsum
    INTEGER(kind=li) :: i

    ! Select quadrature from below (LQ only)
    IF(quad_tpe==1) THEN

      ! Allocate quadrature and check if enough memory is available

      IF(MOD(quad_ord,2).NE. 0 .OR. quad_ord>16) THEN
        IF (rank .EQ. 0) THEN
          WRITE(6,*) 'Level symmetric quadrature is selected. The order must be even &
            & and less than 16.'
        END IF
      END IF
      nangle=quad_ord*(quad_ord+2)/8
      ALLOCATE(quadrature(quad_ord*(quad_ord+2)/8),stat=alloc_stat)
      IF(alloc_stat /= 0) CALL raise_fatal_error("*** Not enough memory ***")

      SELECT CASE(quad_ord)
      CASE(2)
        quadrature(1)%mu%x1=0.577350269190_d_t
        quadrature(1)%mu%x2=0.577350269190_d_t
        quadrature(1)%mu%x3=0.577350269190_d_t
        quadrature(1)%wt=1.000000000000_d_t
      CASE(4)
        quadrature(1)%mu%x1=0.350021174582_d_t
        quadrature(1)%mu%x2=0.350021174582_d_t
        quadrature(1)%mu%x3=0.868890300722_d_t
        quadrature(1)%wt=0.333333333333_d_t
        quadrature(2)%mu%x1=0.350021174582_d_t
        quadrature(2)%mu%x2=0.868890300722_d_t
        quadrature(2)%mu%x3=0.350021174582_d_t
        quadrature(2)%wt=0.333333333333_d_t
        quadrature(3)%mu%x1=0.868890300722_d_t
        quadrature(3)%mu%x2=0.350021174582_d_t
        quadrature(3)%mu%x3=0.350021174582_d_t
        quadrature(3)%wt=0.333333333333_d_t
      CASE(6)
        quadrature(1)%mu%x1=0.266635401517_d_t
        quadrature(1)%mu%x2=0.266635401517_d_t
        quadrature(1)%mu%x3=0.926180935517_d_t
        quadrature(1)%wt=0.176126130863_d_t
        quadrature(2)%mu%x1=0.266635401517_d_t
        quadrature(2)%mu%x2=0.926180935517_d_t
        quadrature(2)%mu%x3=0.266635401517_d_t
        quadrature(2)%wt=0.176126130863_d_t
        quadrature(3)%mu%x1=0.926180935517_d_t
        quadrature(3)%mu%x2=0.266635401517_d_t
        quadrature(3)%mu%x3=0.266635401517_d_t
        quadrature(3)%wt=0.176126130863_d_t
        quadrature(4)%mu%x1=0.681507726537_d_t
        quadrature(4)%mu%x2=0.681507726537_d_t
        quadrature(4)%mu%x3=0.266635401517_d_t
        quadrature(4)%wt=0.157207202470_d_t
        quadrature(5)%mu%x1=0.681507726537_d_t
        quadrature(5)%mu%x2=0.266635401517_d_t
        quadrature(5)%mu%x3=0.681507726537_d_t
        quadrature(5)%wt=0.157207202470_d_t
        quadrature(6)%mu%x1=0.266635401517_d_t
        quadrature(6)%mu%x2=0.681507726537_d_t
        quadrature(6)%mu%x3=0.681507726537_d_t
        quadrature(6)%wt=0.157207202470_d_t
      CASE(8)
        quadrature(1)%mu%x1=0.218217890236_d_t
        quadrature(1)%mu%x2=0.218217890236_d_t
        quadrature(1)%mu%x3=0.951189731211_d_t
        quadrature(1)%wt=0.120987654321_d_t
        quadrature(2)%mu%x1=0.218217890236_d_t
        quadrature(2)%mu%x2=0.951189731211_d_t
        quadrature(2)%mu%x3=0.218217890236_d_t
        quadrature(2)%wt=0.120987654321_d_t
        quadrature(3)%mu%x1=0.951189731211_d_t
        quadrature(3)%mu%x2=0.218217890236_d_t
        quadrature(3)%mu%x3=0.218217890236_d_t
        quadrature(3)%wt=0.120987654321_d_t
        quadrature(4)%mu%x1=0.577350269190_d_t
        quadrature(4)%mu%x2=0.786795792469_d_t
        quadrature(4)%mu%x3=0.218217890236_d_t
        quadrature(4)%wt=0.0907407407405_d_t
        quadrature(5)%mu%x1=0.786795792469_d_t
        quadrature(5)%mu%x2=0.577350269190_d_t
        quadrature(5)%mu%x3=0.218217890236_d_t
        quadrature(5)%wt=0.0907407407405_d_t
        quadrature(6)%mu%x1=0.786795792469_d_t
        quadrature(6)%mu%x2=0.218217890236_d_t
        quadrature(6)%mu%x3=0.577350269190_d_t
        quadrature(6)%wt=0.0907407407405_d_t
        quadrature(7)%mu%x1=0.577350269190_d_t
        quadrature(7)%mu%x2=0.218217890236_d_t
        quadrature(7)%mu%x3=0.786795792469_d_t
        quadrature(7)%wt=0.0907407407405_d_t
        quadrature(8)%mu%x1=0.218217890236_d_t
        quadrature(8)%mu%x2=0.786795792469_d_t
        quadrature(8)%mu%x3=0.577350269190_d_t
        quadrature(8)%wt=0.0907407407405_d_t
        quadrature(9)%mu%x1=0.218217890236_d_t
        quadrature(9)%mu%x2=0.577350269190_d_t
        quadrature(9)%mu%x3=0.786795792469_d_t
        quadrature(9)%wt=0.0907407407405_d_t
        quadrature(10)%mu%x1=0.577350269190_d_t
        quadrature(10)%mu%x2=0.577350269190_d_t
        quadrature(10)%mu%x3=0.577350269190_d_t
        quadrature(10)%wt=0.092592592594_d_t
      CASE(12)
        quadrature(1)%mu%x1=0.167212652823_d_t
        quadrature(1)%mu%x2=0.167212652823_d_t
        quadrature(1)%mu%x3=0.971637719251_d_t
        quadrature(1)%wt=0.0707625899701_d_t
        quadrature(2)%mu%x1=0.167212652823_d_t
        quadrature(2)%mu%x2=0.971637719251_d_t
        quadrature(2)%mu%x3=0.167212652823_d_t
        quadrature(2)%wt=0.0707625899701_d_t
        quadrature(3)%mu%x1=0.971637719251_d_t
        quadrature(3)%mu%x2=0.167212652823_d_t
        quadrature(3)%mu%x3=0.167212652823_d_t
        quadrature(3)%wt=0.0707625899701_d_t
        quadrature(4)%mu%x1=0.167212652823_d_t
        quadrature(4)%mu%x2=0.459547634643_d_t
        quadrature(4)%mu%x3=0.872270543026_d_t
        quadrature(4)%wt=0.055881101565_d_t
        quadrature(5)%mu%x1=0.459547634643_d_t
        quadrature(5)%mu%x2=0.167212652823_d_t
        quadrature(5)%mu%x3=0.872270543026_d_t
        quadrature(5)%wt=0.055881101565_d_t
        quadrature(6)%mu%x1=0.167212652823_d_t
        quadrature(6)%mu%x2=0.872270543026_d_t
        quadrature(6)%mu%x3=0.459547634643_d_t
        quadrature(6)%wt=0.055881101565_d_t
        quadrature(7)%mu%x1=0.459547634643_d_t
        quadrature(7)%mu%x2=0.872270543026_d_t
        quadrature(7)%mu%x3=0.167212652823_d_t
        quadrature(7)%wt=0.055881101565_d_t
        quadrature(8)%mu%x1=0.872270543026_d_t
        quadrature(8)%mu%x2=0.167212652823_d_t
        quadrature(8)%mu%x3=0.459547634643_d_t
        quadrature(8)%wt=0.055881101565_d_t
        quadrature(9)%mu%x1=0.872270543026_d_t
        quadrature(9)%mu%x2=0.459547634643_d_t
        quadrature(9)%mu%x3=0.167212652823_d_t
        quadrature(9)%wt=0.055881101565_d_t
        quadrature(10)%mu%x1=0.167212652823_d_t
        quadrature(10)%mu%x2=0.628019096642_d_t
        quadrature(10)%mu%x3=0.760021014834_d_t
        quadrature(10)%wt=0.037337673759_d_t
        quadrature(11)%mu%x1=0.628019096642_d_t
        quadrature(11)%mu%x2=0.167212652823_d_t
        quadrature(11)%mu%x3=0.760021014834_d_t
        quadrature(11)%wt=0.037337673759_d_t
        quadrature(12)%mu%x1=0.167212652823_d_t
        quadrature(12)%mu%x2=0.760021014834_d_t
        quadrature(12)%mu%x3=0.628019096642_d_t
        quadrature(12)%wt=0.037337673759_d_t
        quadrature(13)%mu%x1=0.628019096642_d_t
        quadrature(13)%mu%x2=0.760021014834_d_t
        quadrature(13)%mu%x3=0.167212652823_d_t
        quadrature(13)%wt=0.037337673759_d_t
        quadrature(14)%mu%x1=0.760021014834_d_t
        quadrature(14)%mu%x2=0.167212652823_d_t
        quadrature(14)%mu%x3=0.628019096642_d_t
        quadrature(14)%wt=0.037337673759_d_t
        quadrature(15)%mu%x1=0.760021014834_d_t
        quadrature(15)%mu%x2=0.628019096642_d_t
        quadrature(15)%mu%x3=0.167212652823_d_t
        quadrature(15)%wt=0.037337673759_d_t
        quadrature(16)%mu%x1=0.459547634643_d_t
        quadrature(16)%mu%x2=0.459547634643_d_t
        quadrature(16)%mu%x3=0.760021014834_d_t
        quadrature(16)%wt=0.050281901060_d_t
        quadrature(17)%mu%x1=0.459547634643_d_t
        quadrature(17)%mu%x2=0.760021014834_d_t
        quadrature(17)%mu%x3=0.459547634643_d_t
        quadrature(17)%wt=0.050281901060_d_t
        quadrature(18)%mu%x1=0.760021014834_d_t
        quadrature(18)%mu%x2=0.459547634643_d_t
        quadrature(18)%mu%x3=0.459547634643_d_t
        quadrature(18)%wt=0.050281901060_d_t
        quadrature(19)%mu%x1=0.459547634643_d_t
        quadrature(19)%mu%x2=0.628019096642_d_t
        quadrature(19)%mu%x3=0.628019096642_d_t
        quadrature(19)%wt=0.025851291656_d_t
        quadrature(20)%mu%x1=0.628019096642_d_t
        quadrature(20)%mu%x2=0.459547634643_d_t
        quadrature(20)%mu%x3=0.628019096642_d_t
        quadrature(20)%wt=0.025851291656_d_t
        quadrature(21)%mu%x1=0.628019096642_d_t
        quadrature(21)%mu%x2=0.628019096642_d_t
        quadrature(21)%mu%x3=0.459547634643_d_t
        quadrature(21)%wt=0.025851291656_d_t
      CASE(16)
        quadrature(1)%mu%x1=0.138956875068_d_t
        quadrature(1)%mu%x2=0.138956875068_d_t
        quadrature(1)%mu%x3=0.980500879012_d_t
        quadrature(1)%wt=0.048987239158_d_t
        quadrature(2)%mu%x1=0.138956875068_d_t
        quadrature(2)%mu%x2=0.980500879012_d_t
        quadrature(2)%mu%x3=0.138956875068_d_t
        quadrature(2)%wt=0.048987239158_d_t
        quadrature(3)%mu%x1=0.138956875068_d_t
        quadrature(3)%mu%x2=0.980500879012_d_t
        quadrature(3)%mu%x3=0.138956875068_d_t
        quadrature(3)%wt=0.048987239158_d_t
        quadrature(4)%mu%x1=0.138956875068_d_t
        quadrature(4)%mu%x2=0.392289261445_d_t
        quadrature(4)%mu%x3=0.909285500944_d_t
        quadrature(4)%wt=0.041329597870_d_t
        quadrature(5)%mu%x1=0.392289261445_d_t
        quadrature(5)%mu%x2=0.138956875068_d_t
        quadrature(5)%mu%x3=0.909285500944_d_t
        quadrature(5)%wt=0.041329597870_d_t
        quadrature(6)%mu%x1=0.138956875068_d_t
        quadrature(6)%mu%x2=0.909285500944_d_t
        quadrature(6)%mu%x3=0.392289261445_d_t
        quadrature(6)%wt=0.041329597870_d_t
        quadrature(7)%mu%x1=0.392289261445_d_t
        quadrature(7)%mu%x2=0.909285500944_d_t
        quadrature(7)%mu%x3=0.138956875068_d_t
        quadrature(7)%wt=0.041329597870_d_t
        quadrature(8)%mu%x1=0.909285500944_d_t
        quadrature(8)%mu%x2=0.138956875068_d_t
        quadrature(8)%mu%x3=0.392289261445_d_t
        quadrature(8)%wt=0.041329597870_d_t
        quadrature(9)%mu%x1=0.909285500944_d_t
        quadrature(9)%mu%x2=0.392289261445_d_t
        quadrature(9)%mu%x3=0.138956875068_d_t
        quadrature(9)%wt=0.041329597870_d_t
        quadrature(10)%mu%x1=0.138956875068_d_t
        quadrature(10)%mu%x2=0.537096561301_d_t
        quadrature(10)%mu%x3=0.831996556910_d_t
        quadrature(10)%wt=0.021232626483_d_t
        quadrature(11)%mu%x1=0.537096561301_d_t
        quadrature(11)%mu%x2=0.138956875068_d_t
        quadrature(11)%mu%x3=0.831996556910_d_t
        quadrature(11)%wt=0.021232626483_d_t
        quadrature(12)%mu%x1=0.138956875068_d_t
        quadrature(12)%mu%x2=0.831996556910_d_t
        quadrature(12)%mu%x3=0.537096561301_d_t
        quadrature(12)%wt=0.021232626483_d_t
        quadrature(13)%mu%x1=0.537096561301_d_t
        quadrature(13)%mu%x2=0.831996556910_d_t
        quadrature(13)%mu%x3=0.138956875068_d_t
        quadrature(13)%wt=0.021232626483_d_t
        quadrature(14)%mu%x1=0.831996556910_d_t
        quadrature(14)%mu%x2=0.138956875068_d_t
        quadrature(14)%mu%x3=0.537096561301_d_t
        quadrature(14)%wt=0.021232626483_d_t
        quadrature(15)%mu%x1=0.831996556910_d_t
        quadrature(15)%mu%x2=0.537096561301_d_t
        quadrature(15)%mu%x3=0.138956875068_d_t
        quadrature(15)%wt=0.021232626483_d_t
        quadrature(16)%mu%x1=0.138956875068_d_t
        quadrature(16)%mu%x2=0.650426450629_d_t
        quadrature(16)%mu%x3=0.746750573615_d_t
        quadrature(16)%wt=0.025620650037_d_t
        quadrature(17)%mu%x1=0.746750573615_d_t
        quadrature(17)%mu%x2=0.138956875068_d_t
        quadrature(17)%mu%x3=0.746750573615_d_t
        quadrature(17)%wt=0.025620650037_d_t
        quadrature(18)%mu%x1=0.138956875068_d_t
        quadrature(18)%mu%x2=0.746750573615_d_t
        quadrature(18)%mu%x3=0.746750573615_d_t
        quadrature(18)%wt=0.025620650037_d_t
        quadrature(19)%mu%x1=0.746750573615_d_t
        quadrature(19)%mu%x2=0.746750573615_d_t
        quadrature(19)%mu%x3=0.138956875068_d_t
        quadrature(19)%wt=0.025620650037_d_t
        quadrature(20)%mu%x1=0.746750573615_d_t
        quadrature(20)%mu%x2=0.138956875068_d_t
        quadrature(20)%mu%x3=0.746750573615_d_t
        quadrature(20)%wt=0.025620650037_d_t
        quadrature(21)%mu%x1=0.746750573615_d_t
        quadrature(21)%mu%x2=0.746750573615_d_t
        quadrature(21)%mu%x3=0.138956875068_d_t
        quadrature(21)%wt=0.025620650037_d_t
        quadrature(22)%mu%x1=0.392289261445_d_t
        quadrature(22)%mu%x2=0.392289261445_d_t
        quadrature(22)%mu%x3=0.831996556910_d_t
        quadrature(22)%wt=0.036048589308_d_t
        quadrature(23)%mu%x1=0.392289261445_d_t
        quadrature(23)%mu%x2=0.831996556910_d_t
        quadrature(23)%mu%x3=0.392289261445_d_t
        quadrature(23)%wt=0.036048589308_d_t
        quadrature(24)%mu%x1=0.831996556910_d_t
        quadrature(24)%mu%x2=0.392289261445_d_t
        quadrature(24)%mu%x3=0.392289261445_d_t
        quadrature(24)%wt=0.036048589308_d_t
        quadrature(25)%mu%x1=0.392289261445_d_t
        quadrature(25)%mu%x2=0.537096561301_d_t
        quadrature(25)%mu%x3=0.746750573615_d_t
        quadrature(25)%wt=0.014458930523_d_t
        quadrature(26)%mu%x1=0.537096561301_d_t
        quadrature(26)%mu%x2=0.392289261445_d_t
        quadrature(26)%mu%x3=0.746750573615_d_t
        quadrature(26)%wt=0.014458930523_d_t
        quadrature(27)%mu%x1=0.392289261445_d_t
        quadrature(27)%mu%x2=0.746750573615_d_t
        quadrature(27)%mu%x3=0.537096561301_d_t
        quadrature(27)%wt=0.014458930523_d_t
        quadrature(28)%mu%x1=0.537096561301_d_t
        quadrature(28)%mu%x2=0.746750573615_d_t
        quadrature(28)%mu%x3=0.392289261445_d_t
        quadrature(28)%wt=0.014458930523_d_t
        quadrature(29)%mu%x1=0.746750573615_d_t
        quadrature(29)%mu%x2=0.392289261445_d_t
        quadrature(29)%mu%x3=0.537096561301_d_t
        quadrature(29)%wt=0.014458930523_d_t
        quadrature(30)%mu%x1=0.746750573615_d_t
        quadrature(30)%mu%x2=0.537096561301_d_t
        quadrature(30)%mu%x3=0.392289261445_d_t
        quadrature(30)%wt=0.014458930523_d_t
        quadrature(31)%mu%x1=0.392289261445_d_t
        quadrature(31)%mu%x2=0.650426450629_d_t
        quadrature(31)%mu%x3=0.650426450629_d_t
        quadrature(31)%wt=0.034495788690_d_t
        quadrature(32)%mu%x1=0.650426450629_d_t
        quadrature(32)%mu%x2=0.392289261445_d_t
        quadrature(32)%mu%x3=0.650426450629_d_t
        quadrature(32)%wt=0.034495788690_d_t
        quadrature(33)%mu%x1=0.650426450629_d_t
        quadrature(33)%mu%x2=0.650426450629_d_t
        quadrature(33)%mu%x3=0.392289261445_d_t
        quadrature(33)%wt=0.034495788690_d_t
        quadrature(34)%mu%x1=0.537096561301_d_t
        quadrature(34)%mu%x2=0.537096561301_d_t
        quadrature(34)%mu%x3=0.650426450629_d_t
        quadrature(34)%wt=0.008518106351_d_t
        quadrature(35)%mu%x1=0.537096561301_d_t
        quadrature(35)%mu%x2=0.650426450629_d_t
        quadrature(35)%mu%x3=0.537096561301_d_t
        quadrature(35)%wt=0.008518106351_d_t
        quadrature(36)%mu%x1=0.650426450629_d_t
        quadrature(36)%mu%x2=0.537096561301_d_t
        quadrature(36)%mu%x3=0.537096561301_d_t
        quadrature(36)%wt=0.008518106351_d_t
      END SELECT

      wsum=zero
      DO i=1,nangle
        wsum=wsum+quadrature(i)%wt
      END DO
      DO i=1,nangle
        quadrature(i)%wt=quadrature(i)%wt/wsum
      END DO
    ELSE IF(quad_tpe==2) THEN

      IF(MOD(quad_ord,2).NE. 0 .OR. quad_ord>26) THEN
        CALL raise_fatal_error("SLC quadrature only available for orders 1, 2, 3 and 5")
      END IF
      ! Square Legendre-Chebychev quadrature

      nangle=quad_ord*(quad_ord+2_li)/8
      ALLOCATE(quadrature(nangle),stat=alloc_stat)
      IF(alloc_stat /= 0) CALL raise_fatal_error("*** Not enough memory ***")

      SELECT CASE(quad_ord)
      CASE(2)
        quadrature( 1 )%mu%x1 = 0.5773502691896257645091488_d_t
        quadrature( 1 )%mu%x2 = 0.5773502691896257645091488_d_t
        quadrature( 1 )%mu%x3 = 0.5773502691896257645091488_d_t
        quadrature( 1 )%wt = 1.0_d_t
      CASE(4)
        quadrature( 1 )%mu%x1 = 0.35988785622265202697900193_d_t
        quadrature( 1 )%mu%x2 = 0.8688461434261049924200219_d_t
        quadrature( 1 )%mu%x3 = 0.3399810435848562648026658_d_t
        quadrature( 1 )%wt = 0.326072577431273071313468_d_t
        quadrature( 2 )%mu%x1 = 0.8688461434261049924200219_d_t
        quadrature( 2 )%mu%x2 = 0.3598878562226520269790019_d_t
        quadrature( 2 )%mu%x3 = 0.3399810435848562648026658_d_t
        quadrature( 2 )%wt = 0.326072577431273071313468_d_t
        quadrature( 3 )%mu%x1 = 0.359474792477991975479987_d_t
        quadrature( 3 )%mu%x2 = 0.359474792477991975479987_d_t
        quadrature( 3 )%mu%x3 = 0.8611363115940525752239465_d_t
        quadrature( 3 )%wt = 0.347854845137453857373064_d_t
      CASE(6)
        quadrature( 1 )%mu%x1 = 0.25134259601688142911001209_d_t
        quadrature( 1 )%mu%x2 = 0.93802333844126041890357974_d_t
        quadrature( 1 )%mu%x3 = 0.2386191860831969086305017_d_t
        quadrature( 1 )%wt = 0.155971311524230349129957_d_t
        quadrature( 2 )%mu%x1 = 0.68668074242437898979356765_d_t
        quadrature( 2 )%mu%x2 = 0.6866807424243789897935676_d_t
        quadrature( 2 )%mu%x3 = 0.2386191860831969086305017_d_t
        quadrature( 2 )%wt = 0.155971311524230349129957_d_t
        quadrature( 3 )%mu%x1 = 0.93802333844126041890357974_d_t
        quadrature( 3 )%mu%x2 = 0.2513425960168814291100121_d_t
        quadrature( 3 )%mu%x3 = 0.2386191860831969086305017_d_t
        quadrature( 3 )%wt = 0.155971311524230349129957_d_t
        quadrature( 4 )%mu%x1 = 0.2870896484226266406169828_d_t
        quadrature( 4 )%mu%x2 = 0.6930957228388288698287195_d_t
        quadrature( 4 )%mu%x3 = 0.6612093864662645136613996_d_t
        quadrature( 4 )%wt = 0.18038078652406930378492_d_t
        quadrature( 5 )%mu%x1 = 0.6930957228388288698287195_d_t
        quadrature( 5 )%mu%x2 = 0.287089648422626640616983_d_t
        quadrature( 5 )%mu%x3 = 0.6612093864662645136613996_d_t
        quadrature( 5 )%wt = 0.18038078652406930378492_d_t
        quadrature( 6 )%mu%x1 = 0.255441387681927592007806_d_t
        quadrature( 6 )%mu%x2 = 0.255441387681927592007806_d_t
        quadrature( 6 )%mu%x3 = 0.9324695142031520278123016_d_t
        quadrature( 6 )%wt = 0.1713244923791703450403_d_t
      CASE(8)
        quadrature( 1 )%mu%x1 = 0.191780011462650977226973348_d_t
        quadrature( 1 )%mu%x2 = 0.96414322542653309175199058_d_t
        quadrature( 1 )%mu%x3 = 0.1834346424956498049394761_d_t
        quadrature( 1 )%wt = 0.090670945844590495741288_d_t
        quadrature( 2 )%mu%x1 = 0.54614326613289737776850015_d_t
        quadrature( 2 )%mu%x2 = 0.8173611593354460189712852_d_t
        quadrature( 2 )%mu%x3 = 0.1834346424956498049394761_d_t
        quadrature( 2 )%wt = 0.090670945844590495741288_d_t
        quadrature( 3 )%mu%x1 = 0.8173611593354460189712852_d_t
        quadrature( 3 )%mu%x2 = 0.54614326613289737776850015_d_t
        quadrature( 3 )%mu%x3 = 0.1834346424956498049394761_d_t
        quadrature( 3 )%wt = 0.090670945844590495741288_d_t
        quadrature( 4 )%mu%x1 = 0.96414322542653309175199058_d_t
        quadrature( 4 )%mu%x2 = 0.1917800114626509772269733_d_t
        quadrature( 4 )%mu%x3 = 0.1834346424956498049394761_d_t
        quadrature( 4 )%wt = 0.090670945844590495741288_d_t
        quadrature( 5 )%mu%x1 = 0.22019640583286784184941713_d_t
        quadrature( 5 )%mu%x2 = 0.8217841742123186716931386_d_t
        quadrature( 5 )%mu%x3 = 0.525532409916328985817739_d_t
        quadrature( 5 )%wt = 0.10456888195929576244599_d_t
        quadrature( 6 )%mu%x1 = 0.6015877683794508298437215_d_t
        quadrature( 6 )%mu%x2 = 0.6015877683794508298437215_d_t
        quadrature( 6 )%mu%x3 = 0.525532409916328985817739_d_t
        quadrature( 6 )%wt = 0.10456888195929576244599_d_t
        quadrature( 7 )%mu%x1 = 0.8217841742123186716931386_d_t
        quadrature( 7 )%mu%x2 = 0.220196405832867841849417_d_t
        quadrature( 7 )%mu%x3 = 0.525532409916328985817739_d_t
        quadrature( 7 )%wt = 0.10456888195929576244599_d_t
        quadrature( 8 )%mu%x1 = 0.2313011996193396979048773_d_t
        quadrature( 8 )%mu%x2 = 0.5584104931141764685268045_d_t
        quadrature( 8 )%mu%x3 = 0.7966664774136267395915539_d_t
        quadrature( 8 )%wt = 0.1111905172266872352722_d_t
        quadrature( 9 )%mu%x1 = 0.558410493114176468526805_d_t
        quadrature( 9 )%mu%x2 = 0.231301199619339697904877_d_t
        quadrature( 9 )%mu%x3 = 0.7966664774136267395915539_d_t
        quadrature( 9 )%wt = 0.1111905172266872352722_d_t
        quadrature( 10 )%mu%x1 = 0.197285822485982594739779_d_t
        quadrature( 10 )%mu%x2 = 0.197285822485982594739779_d_t
        quadrature( 10 )%mu%x3 = 0.9602898564975362316835609_d_t
        quadrature( 10 )%wt = 0.1012285362903762591525_d_t
      CASE(10)
        quadrature( 1 )%mu%x1 = 0.15469117853983970604829247_d_t
        quadrature( 1 )%mu%x2 = 0.97668166281278046592604053_d_t
        quadrature( 1 )%mu%x3 = 0.148874338981631210884826_d_t
        quadrature( 1 )%wt = 0.0591048449429505740347786_d_t
        quadrature( 2 )%mu%x1 = 0.44893128526722284995566994_d_t
        quadrature( 2 )%mu%x2 = 0.88107725671538119766654488_d_t
        quadrature( 2 )%mu%x3 = 0.148874338981631210884826_d_t
        quadrature( 2 )%wt = 0.0591048449429505740347786_d_t
        quadrature( 3 )%mu%x1 = 0.69922686990446182426345806_d_t
        quadrature( 3 )%mu%x2 = 0.69922686990446182426345806_d_t
        quadrature( 3 )%mu%x3 = 0.148874338981631210884826_d_t
        quadrature( 3 )%wt = 0.0591048449429505740347786_d_t
        quadrature( 4 )%mu%x1 = 0.88107725671538119766654488_d_t
        quadrature( 4 )%mu%x2 = 0.44893128526722284995566994_d_t
        quadrature( 4 )%mu%x3 = 0.148874338981631210884826_d_t
        quadrature( 4 )%wt = 0.0591048449429505740347786_d_t
        quadrature( 5 )%mu%x1 = 0.97668166281278046592604053_d_t
        quadrature( 5 )%mu%x2 = 0.1546911785398397060482925_d_t
        quadrature( 5 )%mu%x3 = 0.148874338981631210884826_d_t
        quadrature( 5 )%wt = 0.0591048449429505740347786_d_t
        quadrature( 6 )%mu%x1 = 0.17581615504536439384202461_d_t
        quadrature( 6 )%mu%x2 = 0.883887499613281598548222_d_t
        quadrature( 6 )%mu%x3 = 0.4333953941292471907992659_d_t
        quadrature( 6 )%wt = 0.06731667982749908877281_d_t
        quadrature( 7 )%mu%x1 = 0.5006820493078507269226395_d_t
        quadrature( 7 )%mu%x2 = 0.7493236402572959148727155_d_t
        quadrature( 7 )%mu%x3 = 0.4333953941292471907992659_d_t
        quadrature( 7 )%wt = 0.06731667982749908877281_d_t
        quadrature( 8 )%mu%x1 = 0.7493236402572959148727155_d_t
        quadrature( 8 )%mu%x2 = 0.5006820493078507269226395_d_t
        quadrature( 8 )%mu%x3 = 0.4333953941292471907992659_d_t
        quadrature( 8 )%wt = 0.06731667982749908877281_d_t
        quadrature( 9 )%mu%x1 = 0.883887499613281598548222_d_t
        quadrature( 9 )%mu%x2 = 0.175816155045364393842025_d_t
        quadrature( 9 )%mu%x3 = 0.4333953941292471907992659_d_t
        quadrature( 9 )%wt = 0.06731667982749908877281_d_t
        quadrature( 10 )%mu%x1 = 0.1899108686923032762948602_d_t
        quadrature( 10 )%mu%x2 = 0.7087570108692174579781006_d_t
        quadrature( 10 )%mu%x3 = 0.6794095682990244062343274_d_t
        quadrature( 10 )%wt = 0.0730287875053273479985_d_t
        quadrature( 11 )%mu%x1 = 0.5188461421769141816832404_d_t
        quadrature( 11 )%mu%x2 = 0.51884614217691418168324_d_t
        quadrature( 11 )%mu%x3 = 0.6794095682990244062343274_d_t
        quadrature( 11 )%wt = 0.0730287875053273479985_d_t
        quadrature( 12 )%mu%x1 = 0.7087570108692174579781006_d_t
        quadrature( 12 )%mu%x2 = 0.18991086869230327629486_d_t
        quadrature( 12 )%mu%x3 = 0.6794095682990244062343274_d_t
        quadrature( 12 )%wt = 0.0730287875053273479985_d_t
        quadrature( 13 )%mu%x1 = 0.1919779684697349464999435_d_t
        quadrature( 13 )%mu%x2 = 0.463475815156468523900573_d_t
        quadrature( 13 )%mu%x3 = 0.8650633666889845107320967_d_t
        quadrature( 13 )%wt = 0.0747256745752902965729_d_t
        quadrature( 14 )%mu%x1 = 0.463475815156468523900573_d_t
        quadrature( 14 )%mu%x2 = 0.191977968469734946499943_d_t
        quadrature( 14 )%mu%x3 = 0.8650633666889845107320967_d_t
        quadrature( 14 )%wt = 0.0747256745752902965729_d_t
        quadrature( 15 )%mu%x1 = 0.160477527572603428321382_d_t
        quadrature( 15 )%mu%x2 = 0.160477527572603428321382_d_t
        quadrature( 15 )%mu%x3 = 0.973906528517171720077964_d_t
        quadrature( 15 )%wt = 0.0666713443086881375936_d_t
      CASE(12)
        quadrature( 1 )%mu%x1 = 0.129498599586659577006058655_d_t
        quadrature( 1 )%mu%x2 = 0.98363952040251694892577123_d_t
        quadrature( 1 )%mu%x3 = 0.1252334085114689154724414_d_t
        quadrature( 1 )%wt = 0.0415245076355671308334271_d_t
        quadrature( 2 )%mu%x1 = 0.379670683204702274644847654_d_t
        quadrature( 2 )%mu%x2 = 0.91660611262825110569866883_d_t
        quadrature( 2 )%mu%x3 = 0.1252334085114689154724414_d_t
        quadrature( 2 )%wt = 0.0415245076355671308334271_d_t
        quadrature( 3 )%mu%x1 = 0.603968837197814674280923572_d_t
        quadrature( 3 )%mu%x2 = 0.78710751304159152869261018_d_t
        quadrature( 3 )%mu%x3 = 0.1252334085114689154724414_d_t
        quadrature( 3 )%wt = 0.0415245076355671308334271_d_t
        quadrature( 4 )%mu%x1 = 0.78710751304159152869261018_d_t
        quadrature( 4 )%mu%x2 = 0.60396883719781467428092357_d_t
        quadrature( 4 )%mu%x3 = 0.1252334085114689154724414_d_t
        quadrature( 4 )%wt = 0.0415245076355671308334271_d_t
        quadrature( 5 )%mu%x1 = 0.91660611262825110569866883_d_t
        quadrature( 5 )%mu%x2 = 0.37967068320470227464484765_d_t
        quadrature( 5 )%mu%x3 = 0.1252334085114689154724414_d_t
        quadrature( 5 )%wt = 0.0415245076355671308334271_d_t
        quadrature( 6 )%mu%x1 = 0.98363952040251694892577123_d_t
        quadrature( 6 )%mu%x2 = 0.1294985995866595770060587_d_t
        quadrature( 6 )%mu%x3 = 0.1252334085114689154724414_d_t
        quadrature( 6 )%wt = 0.0415245076355671308334271_d_t
        quadrature( 7 )%mu%x1 = 0.14546722962429073158076143_d_t
        quadrature( 7 )%mu%x2 = 0.9184439413759479071797386_d_t
        quadrature( 7 )%mu%x3 = 0.3678314989981801937526915_d_t
        quadrature( 7 )%wt = 0.04669850730767096175217_d_t
        quadrature( 8 )%mu%x1 = 0.42216234290746094077632518_d_t
        quadrature( 8 )%mu%x2 = 0.8285402492188506385915938_d_t
        quadrature( 8 )%mu%x3 = 0.3678314989981801937526915_d_t
        quadrature( 8 )%wt = 0.04669850730767096175217_d_t
        quadrature( 9 )%mu%x1 = 0.6575332646888489409452314_d_t
        quadrature( 9 )%mu%x2 = 0.6575332646888489409452314_d_t
        quadrature( 9 )%mu%x3 = 0.3678314989981801937526915_d_t
        quadrature( 9 )%wt = 0.04669850730767096175217_d_t
        quadrature( 10 )%mu%x1 = 0.8285402492188506385915938_d_t
        quadrature( 10 )%mu%x2 = 0.4221623429074609407763252_d_t
        quadrature( 10 )%mu%x3 = 0.3678314989981801937526915_d_t
        quadrature( 10 )%wt = 0.04669850730767096175217_d_t
        quadrature( 11 )%mu%x1 = 0.9184439413759479071797386_d_t
        quadrature( 11 )%mu%x2 = 0.145467229624290731580761_d_t
        quadrature( 11 )%mu%x3 = 0.3678314989981801937526915_d_t
        quadrature( 11 )%wt = 0.04669850730767096175217_d_t
        quadrature( 12 )%mu%x1 = 0.15789758121964557206310076_d_t
        quadrature( 12 )%mu%x2 = 0.7938047457766728219797501_d_t
        quadrature( 12 )%mu%x3 = 0.5873179542866174472967024_d_t
        quadrature( 12 )%wt = 0.0507918566807664804373_d_t
        quadrature( 13 )%mu%x1 = 0.4496542682633837369147182_d_t
        quadrature( 13 )%mu%x2 = 0.6729551690901138115304159_d_t
        quadrature( 13 )%mu%x3 = 0.5873179542866174472967024_d_t
        quadrature( 13 )%wt = 0.0507918566807664804373_d_t
        quadrature( 14 )%mu%x1 = 0.6729551690901138115304159_d_t
        quadrature( 14 )%mu%x2 = 0.449654268263383736914718_d_t
        quadrature( 14 )%mu%x3 = 0.5873179542866174472967024_d_t
        quadrature( 14 )%wt = 0.0507918566807664804373_d_t
        quadrature( 15 )%mu%x1 = 0.7938047457766728219797501_d_t
        quadrature( 15 )%mu%x2 = 0.157897581219645572063101_d_t
        quadrature( 15 )%mu%x3 = 0.5873179542866174472967024_d_t
        quadrature( 15 )%wt = 0.0507918566807664804373_d_t
        quadrature( 16 )%mu%x1 = 0.165168303853986527475871_d_t
        quadrature( 16 )%mu%x2 = 0.616416501783052127706427_d_t
        quadrature( 16 )%mu%x3 = 0.7699026741943046870368938_d_t
        quadrature( 16 )%wt = 0.0533594428477820754449_d_t
        quadrature( 17 )%mu%x1 = 0.4512481979290656002305559_d_t
        quadrature( 17 )%mu%x2 = 0.451248197929065600230556_d_t
        quadrature( 17 )%mu%x3 = 0.7699026741943046870368938_d_t
        quadrature( 17 )%wt = 0.0533594428477820754449_d_t
        quadrature( 18 )%mu%x1 = 0.6164165017830521277064269_d_t
        quadrature( 18 )%mu%x2 = 0.165168303853986527475871_d_t
        quadrature( 18 )%mu%x3 = 0.7699026741943046870368938_d_t
        quadrature( 18 )%wt = 0.0533594428477820754449_d_t
        quadrature( 19 )%mu%x1 = 0.1635146734385509748740351_d_t
        quadrature( 19 )%mu%x2 = 0.394759342262357452235565_d_t
        quadrature( 19 )%mu%x3 = 0.9041172563704748566784659_d_t
        quadrature( 19 )%wt = 0.05346966299765921548_d_t
        quadrature( 20 )%mu%x1 = 0.394759342262357452235565_d_t
        quadrature( 20 )%mu%x2 = 0.163514673438550974874035_d_t
        quadrature( 20 )%mu%x3 = 0.9041172563704748566784659_d_t
        quadrature( 20 )%wt = 0.05346966299765921548_d_t
        quadrature( 21 )%mu%x1 = 0.135164198841960802123231_d_t
        quadrature( 21 )%mu%x2 = 0.13516419884196080212323_d_t
        quadrature( 21 )%mu%x3 = 0.9815606342467192506905491_d_t
        quadrature( 21 )%wt = 0.0471753363865118271946_d_t
      END SELECT

      !  else if(quad_tpe==2) then
      !
      !    if(quad_ord < 1 .or. quad_ord >  5 .or. quad_ord .eq. 4) then
      !       call raise_fatal_error("SLC quadrature only available for orders 1, 2, 3 and 5")
      !    end if
      !
      !  ! Square Legendre-Chebychev quadrature
      !
      !    nangle=quad_ord**2
      !    allocate(quadrature(nangle),stat=alloc_stat)
      !    if(alloc_stat /= 0) call raise_fatal_error("*** Not enough memory ***")
      !
      !    select case(quad_ord)
      !      case(1)
      !        quadrature(1)%mu%x1= 5.773500D-01
      !        quadrature(1)%mu%x2= 5.773500D-01
      !        quadrature(1)%wt= 1.000000D+00
      !      case(2)
      !        quadrature(1)%mu%x1= 3.399810D-01
      !        quadrature(1)%mu%x2= 8.688460D-01
      !        quadrature(1)%wt= 3.260724D-01
      !        quadrature(2)%mu%x1= 8.611360D-01
      !        quadrature(2)%mu%x2= 4.696760D-01
      !        quadrature(2)%wt= 1.739276D-01
      !        quadrature(3)%mu%x1= 3.399810D-01
      !        quadrature(3)%mu%x2= 3.598880D-01
      !        quadrature(3)%wt= 3.260724D-01
      !        quadrature(4)%mu%x1= 8.611360D-01
      !        quadrature(4)%mu%x2= 1.945460D-01
      !        quadrature(4)%wt= 1.739276D-01
      !      case(3)
      !        quadrature(1)%mu%x1= 6.612090D-01
      !        quadrature(1)%mu%x2= 7.246390D-01
      !        quadrature(1)%wt= 1.202540D-01
      !        quadrature(2)%mu%x1= 2.386190D-01
      !        quadrature(2)%mu%x2= 6.866810D-01
      !        quadrature(2)%wt= 1.559712D-01
      !        quadrature(3)%mu%x1= 6.612090D-01
      !        quadrature(3)%mu%x2= 5.304730D-01
      !        quadrature(3)%wt= 1.202540D-01
      !        quadrature(4)%mu%x1= 9.324690D-01
      !        quadrature(4)%mu%x2= 3.489390D-01
      !        quadrature(4)%wt= 5.710800D-02
      !        quadrature(5)%mu%x1= 9.324690D-01
      !        quadrature(5)%mu%x2= 2.554410D-01
      !        quadrature(5)%wt= 5.710800D-02
      !        quadrature(6)%mu%x1= 2.386190D-01
      !        quadrature(6)%mu%x2= 2.513430D-01
      !        quadrature(6)%wt= 1.559712D-01
      !        quadrature(7)%mu%x1= 6.612090D-01
      !        quadrature(7)%mu%x2= 1.941660D-01
      !        quadrature(7)%wt= 1.202540D-01
      !        quadrature(8)%mu%x1= 9.324690D-01
      !        quadrature(8)%mu%x2= 9.349800D-02
      !        quadrature(8)%wt= 5.710800D-02
      !      case(5)
      !        quadrature(1)%mu%x1= 1.488740D-01
      !        quadrature(1)%mu%x2= 9.766820D-01
      !        quadrature(1)%wt= 5.910480D-02
      !        quadrature(2)%mu%x1= 4.333950D-01
      !        quadrature(2)%mu%x2= 8.901090D-01
      !        quadrature(2)%wt= 5.385320D-02
      !        quadrature(3)%mu%x1= 1.488740D-01
      !        quadrature(3)%mu%x2= 8.810770D-01
      !        quadrature(3)%wt= 5.910480D-02
      !        quadrature(4)%mu%x1= 4.333950D-01
      !        quadrature(4)%mu%x2= 8.029780D-01
      !        quadrature(4)%wt= 5.385320D-02
      !        quadrature(5)%mu%x1= 6.794100D-01
      !        quadrature(5)%mu%x2= 7.247260D-01
      !        quadrature(5)%wt= 4.381720D-02
      !        quadrature(6)%mu%x1= 1.488740D-01
      !        quadrature(6)%mu%x2= 6.992270D-01
      !        quadrature(6)%wt= 5.910480D-02
      !        quadrature(7)%mu%x1= 6.794100D-01
      !        quadrature(7)%mu%x2= 6.537840D-01
      !        quadrature(7)%wt= 4.381720D-02
      !        quadrature(8)%mu%x1= 4.333950D-01
      !        quadrature(8)%mu%x2= 6.372470D-01
      !        quadrature(8)%wt= 5.385320D-02
      !        quadrature(9)%mu%x1= 6.794100D-01
      !        quadrature(9)%mu%x2= 5.188460D-01
      !        quadrature(9)%wt= 4.381720D-02
      !        quadrature(10)%mu%x1= 8.650630D-01
      !        quadrature(10)%mu%x2= 4.954860D-01
      !        quadrature(10)%wt= 2.989028D-02
      !        quadrature(11)%mu%x1= 1.488740D-01
      !        quadrature(11)%mu%x2= 4.489310D-01
      !        quadrature(11)%wt= 5.910480D-02
      !        quadrature(12)%mu%x1= 8.650630D-01
      !        quadrature(12)%mu%x2= 4.469850D-01
      !        quadrature(12)%wt= 2.989028D-02
      !        quadrature(13)%mu%x1= 4.333950D-01
      !        quadrature(13)%mu%x2= 4.091380D-01
      !        quadrature(13)%wt= 5.385320D-02
      !        quadrature(14)%mu%x1= 8.650630D-01
      !        quadrature(14)%mu%x2= 3.547290D-01
      !        quadrature(14)%wt= 2.989028D-02
      !        quadrature(15)%mu%x1= 6.794100D-01
      !        quadrature(15)%mu%x2= 3.331200D-01
      !        quadrature(15)%wt= 4.381720D-02
      !        quadrature(16)%mu%x1= 8.650630D-01
      !        quadrature(16)%mu%x2= 2.277500D-01
      !        quadrature(16)%wt= 2.989028D-02
      !        quadrature(17)%mu%x1= 9.739070D-01
      !        quadrature(17)%mu%x2= 2.241550D-01
      !        quadrature(17)%wt= 1.333428D-02
      !        quadrature(18)%mu%x1= 9.739070D-01
      !        quadrature(18)%mu%x2= 2.022130D-01
      !        quadrature(18)%wt= 1.333428D-02
      !        quadrature(19)%mu%x1= 9.739070D-01
      !        quadrature(19)%mu%x2= 1.604770D-01
      !        quadrature(19)%wt= 1.333428D-02
      !        quadrature(20)%mu%x1= 1.488740D-01
      !        quadrature(20)%mu%x2= 1.546910D-01
      !        quadrature(20)%wt= 5.910480D-02
      !        quadrature(21)%mu%x1= 4.333950D-01
      !        quadrature(21)%mu%x2= 1.409790D-01
      !        quadrature(21)%wt= 5.385320D-02
      !        quadrature(22)%mu%x1= 6.794100D-01
      !        quadrature(22)%mu%x2= 1.147850D-01
      !        quadrature(22)%wt= 4.381720D-02
      !        quadrature(23)%mu%x1= 9.739070D-01
      !        quadrature(23)%mu%x2= 1.030330D-01
      !        quadrature(23)%wt= 1.333428D-02
      !        quadrature(24)%mu%x1= 8.650630D-01
      !        quadrature(24)%mu%x2= 7.847730D-02
      !        quadrature(24)%wt= 2.989028D-02
      !        quadrature(25)%mu%x1= 9.739070D-01
      !        quadrature(25)%mu%x2= 3.550270D-02
      !        quadrature(25)%wt= 1.333428D-02
      !    end select

      ! Compute xi and renormalize weights

      !        wsum=zero
      !        do i=1,nangle
      !          wsum=wsum+quadrature(i)%wt
      !          quadrature(i)%mu%x3=sqrt(one - quadrature(i)%mu%x1**2 - &
      !                                         quadrature(i)%mu%x2**2 )
      !        end do
      !        do i=1,nangle
      !           quadrature(i)%wt=quadrature(i)%wt/wsum
      !        end do

    ELSE IF(quad_tpe==3) THEN

    END IF
  END SUBROUTINE quad_gen

END MODULE quadrature_module
