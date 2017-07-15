! Copyright (c) 2005-2010, 2012-2013, Andrew Hang Chen and contributors,
! All rights reserved.
! Licensed under the 3-clause BSD license.

module stringmod_test
  use fruit
  implicit none

contains

!===============================================================================
!Subroutine > test_uppercase
!===============================================================================
!> Tests the stringmod function "uppercase"
!> Passes lower, upper, numeric, & quoted data
  subroutine test_uppercase
    use precision
    use strings, only: uppercase

    character(50):: give = "abcdefghijklmnopqrstuvwxyz0123456789""a""A"
    character(50):: get, expct = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789""a""A"

    get = uppercase(give)
    call assert_equals(expct, get)

  end subroutine

!===============================================================================
!Subroutine > test_is_digit
!===============================================================================
!> Tests the stringmod function "test_is_digit"
!> Passes digit and non digit data
  subroutine test_is_digit
    use precision
    use strings, only: is_digit
    character:: in1 = '4', in2 ='d'

    call assert_true(is_digit(in1))
    call assert_false(is_digit(in2))

  end subroutine

!===============================================================================
!Subroutine > test_is_letter
!===============================================================================
!> Tests the stringmod function "test_is_letter"
!> Passes letter and digit data
  subroutine test_is_letter
    use precision
    use strings, only: is_letter
    character:: in1 = '4', in2 ='d', in3='A'

    call assert_true(is_letter(in3))
    call assert_true(is_letter(in2))
    call assert_false(is_letter(in1))

  end subroutine

!===============================================================================
!Subroutine test_removebksl
!===============================================================================
!> Tests the stringmod function "test_removebksl"
!> Passes a string containing a combination of text, '\\',  '\'
  subroutine test_removebksl
    use precision
    use strings, only: removebksl
    character(50):: give = "   \a\b\\c  "
    character(50):: expct = "ab\c"

    call removebksl(give)
    call assert_equals(expct, give)
  end subroutine

end module stringmod_test
