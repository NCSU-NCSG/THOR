! Copyright (c) 2005-2010, 2012-2013, Andrew Hang Chen and contributors,
! All rights reserved.
! Licensed under the 3-clause BSD license.

module stringmod_test
  use fruit
  implicit none

contains                          !fortran 95 limits subroutine name to 31 char.
  subroutine test_uppercase
    use precision
    use strings, only: uppercase
    character(50):: give = "abcdefghijklmnopqrstuvwxyz0123456789"
    character(50):: get, expct = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"

    get = uppercase(give)
    call assert_equals(get, expct)

  end subroutine 
end module stringmod_test

