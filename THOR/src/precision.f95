MODULE PRECISION

  ! Real kinds

  INTEGER, PARAMETER :: kr4 = SELECTED_REAL_KIND(6,37)       ! single precision real
  INTEGER, PARAMETER :: kr8 = SELECTED_REAL_KIND(15,307)     ! double precision real

  ! Integer kinds

  INTEGER, PARAMETER :: ki4 = SELECTED_INT_KIND(9)           ! single precision integer
  INTEGER, PARAMETER :: ki8 = SELECTED_INT_KIND(18)          ! double precision integer

  !Complex kinds

  INTEGER, PARAMETER :: kc4 = kr4                            ! single precision complex
  INTEGER, PARAMETER :: kc8 = kr8                            ! double precision complex

END MODULE PRECISION
