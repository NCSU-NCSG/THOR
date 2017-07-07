!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! algorithms module
!     Provides general algorithms that are simply not included in standard
!     FORTRAN intrinsics.
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

MODULE integer_array_tools
IMPLICIT NONE

PRIVATE
PUBLIC quickSortInteger, numUniqueEntries, isSorted, hasUniqueEntries, &
       uniqueEntries, mapIndexOf, indexOf

CONTAINS

  !
  ! Adopted from https://rosettacode.org/wiki/Sorting_algorithms/Quicksort#Fortran
  ! Sorts the array values in ascending order and returns order, where order(i)
  ! is the original index of the item that got sorted to position i.
  !
  ! :values: (input, output) array of integers to be sorted
  ! :order: (output) permutation array, see description
  !
  RECURSIVE SUBROUTINE quickSortInteger(values, order, map_int)

    INTEGER, DIMENSION (:), INTENT(INOUT) :: values
    INTEGER, DIMENSION (:), INTENT(OUT) :: order
    INTEGER, DIMENSION (:), INTENT(INOUT), OPTIONAL  :: map_int

    ! LOCAL VARIABLES
    INTEGER :: nA
    INTEGER :: left, right
    REAL :: random
    INTEGER :: pivot
    INTEGER :: temp_value
    INTEGER ::  temp_map_int_value
    INTEGER :: temp_order
    INTEGER :: marker
    LOGICAL :: map_int_present
    INTEGER, POINTER :: null_pass => NULL()

    map_int_present = PRESENT(map_int)
    null_pass => NULL()

    nA = SIZE(values)

    IF (nA .NE. SIZE(order)) THEN
      STOP 'Array to be sorted and ordering array must be of size N'
    END IF

    IF(map_int_present) THEN
      IF (na .NE. SIZE(map_int)) THEN
        STOP 'To apply Quicksort to a map, both the key and data arrays must be of size N'
      END IF
    END IF

    IF (nA > 1) THEN

      CALL random_number(random)
      pivot = values(int(random*REAL(nA-1))+1)
      left = 0
      right = nA + 1

      DO WHILE (left < right)
        right = right - 1
        DO WHILE (values(right) > pivot)
          right = right - 1
        END DO
        left = left + 1
        DO WHILE (values(left) < pivot)
          left = left + 1
        END DO
        IF (left < right) THEN

          !Sort the primary array
          temp_value = values(left)
          values(left) = values(right)
          values(right) = temp_value

          !If present, sort the optional map_int array
          IF (map_int_present) THEN
            temp_map_int_value = map_int(left)
            map_int(left) = map_int(right)
            map_int(right) = temp_map_int_value
          END IF

          !Update the permutation array
          temp_order = order(left)
          order(left) = order(right)
          order(right) = temp_order
        END IF
      END DO

      IF (left == right) THEN
        marker = left + 1
      ELSE
        marker = left
      END IF

      !Recurse for no optional arrays
      IF (.NOT. map_int_present) THEN
        CALL quickSortInteger(values(:marker-1),order(:marker-1))
        CALL quickSortInteger(values(marker:),order(marker:))
      END IF

      !Recurse for optional int array
      IF (map_int_present) THEN
        CALL quickSortInteger(values(:marker-1),order(:marker-1),map_int(:marker-1))
        CALL quickSortInteger(values(marker:),order(marker:), map_int(marker:))
      END IF

    END IF

  END SUBROUTINE quickSortInteger

  !
  ! Assigns 1, 2, 3, ..., N to array of size N
  !
  ! :array: (input, output) the array to be filled with 1:N
  SUBROUTINE range(array)

    INTEGER, DIMENSION (:), INTENT(INOUT)  :: array
    INTEGER :: s
    INTEGER :: j

    s = SIZE(array)
    DO j = 1, s
      array(j) = j
    END DO

  END SUBROUTINE range

  !
  ! Determines the number of unique integers in array values
  ! :values: the input array of integers
  ! :sorted: whether the values array is sorted or not, if omitted => not sorted
  ! :returns: numUniqueEntries number of unique values in the array
  !
  INTEGER FUNCTION numUniqueEntries(values, sorted)

    INTEGER, DIMENSION (:), INTENT(IN)  :: values
    LOGICAL, OPTIONAL :: sorted

    INTEGER :: s
    LOGICAL :: is_sorted
    INTEGER :: j
    INTEGER, DIMENSION(:), ALLOCATABLE :: order
    INTEGER, DIMENSION(:), ALLOCATABLE :: ordered_values

    IF (.NOT. PRESENT(sorted)) THEN
      is_sorted = .FALSE.
    ELSE
      is_sorted = sorted
    END IF

    ! prepare temporary arrays
    s = SIZE(values)
    ALLOCATE(order(s))
    ALLOCATE(ordered_values(s))
    ordered_values = values

    ! determine num_unique
    CALL range(order)
    IF (.NOT. is_sorted) CALL quickSortInteger(ordered_values, order)

    ! determine unique values, note first one is unique so we start with 1
    numUniqueEntries = 1
    DO j = 1, s - 1
      IF (ordered_values(j + 1) .ne. ordered_values(j)) &
        numUniqueEntries = numUniqueEntries + 1
    END DO
    DEALLOCATE(order, ordered_values)
    RETURN
  END FUNCTION numUniqueEntries

  !
  ! Returns the unique integers in array values
  ! :values: the input array of integers
  ! :unique_values: the unique values in array values; note! must have size
  !                 num_unique_entries
  ! :sorted: whether the values array is sorted or not, if omitted => not sorted
  !
  SUBROUTINE uniqueEntries(values, unique_values, sorted)

    INTEGER, DIMENSION (:), INTENT(IN)  :: values
    INTEGER, DIMENSION (:), INTENT(OUT)  :: unique_values
    LOGICAL, OPTIONAL :: sorted

    INTEGER :: s
    LOGICAL :: is_sorted
    INTEGER :: j
    INTEGER :: counter
    INTEGER, DIMENSION(:), ALLOCATABLE :: order
    INTEGER, DIMENSION(:), ALLOCATABLE :: ordered_values

    IF (.NOT. PRESENT(sorted)) THEN
      is_sorted = .FALSE.
    ELSE
      is_sorted = sorted
    END IF

    ! prepare temporary arrays
    s = SIZE(values)
    ALLOCATE(order(s))
    ALLOCATE(ordered_values(s))
    ordered_values = values

    ! determine num_unique
    CALL range(order)
    IF (.NOT. is_sorted) CALL quickSortInteger(ordered_values, order)

    ! determine unique values, note first one is unique so we start with 1
    counter = 1
    unique_values(counter) = ordered_values(1)
    DO j = 1, s - 1
      IF (ordered_values(j + 1) .ne. ordered_values(j)) THEN
        counter = counter + 1
        unique_values(counter) = ordered_values(j + 1)
      END IF
    END DO
    DEALLOCATE(order, ordered_values)
    RETURN
  END SUBROUTINE uniqueEntries

  !
  ! isSorted checks if an array is sorted in ascending order
  ! :array: the input array of integers
  ! :returns: logical true if sorted, otherwise false
  !
  LOGICAL FUNCTION isSorted(array)

    INTEGER, DIMENSION (:), INTENT(IN)  :: array

    INTEGER :: s
    INTEGER :: j

    s = SIZE(array)
    isSorted = .TRUE.
    DO j = 1, s - 1
      IF (array(j + 1) .LT. array(j)) THEN
        isSorted = .FALSE.
        RETURN
      END IF
    END DO
  END FUNCTION isSorted

  !
  ! hasUniqueEntries checks if an array has unique entries
  ! :array: the input array of integers
  ! :sorted: optional argument, set to true if array is sorted already
  ! :returns: logical true if array has all entries only once, false otherwise
  !
  LOGICAL FUNCTION hasUniqueEntries(array, sorted)

    INTEGER, DIMENSION (:), INTENT(IN)  :: array
    LOGICAL, OPTIONAL :: sorted

    INTEGER :: s
    LOGICAL :: is_sorted
    INTEGER :: j
    INTEGER, DIMENSION(:), ALLOCATABLE :: order
    INTEGER, DIMENSION(:), ALLOCATABLE :: ordered_array

    IF (.NOT. PRESENT(sorted)) THEN
      is_sorted = .FALSE.
    ELSE
      is_sorted = sorted
    END IF

    s = SIZE(array)
    ALLOCATE(ordered_array(s))
    ordered_array = array
    IF (.NOT. is_sorted) THEN
      ALLOCATE(order(s))
      CALL range(order)
      CALL quickSortInteger(ordered_array, order)
    END IF

    hasUniqueEntries = .TRUE.
    DO j = 1, s - 1
      IF (ordered_array(j + 1) .EQ. ordered_array(j)) THEN
        hasUniqueEntries = .FALSE.
        DEALLOCATE(ordered_array)
        IF (.NOT. is_sorted) DEALLOCATE(order)
        RETURN
      END IF
    END DO

    DEALLOCATE(ordered_array)
    IF (.NOT. is_sorted) DEALLOCATE(order)

  END FUNCTION hasUniqueEntries

  !Returns the index of the target value in the provided array
  !Does not require sorted or unique data
  !If target occurs multiple times, first occurence will be returned
  !A return of -1 implies item not found
  INTEGER FUNCTION indexOf(key, array)
    INTEGER, INTENT(IN) :: key
    INTEGER, INTENT(IN) :: array(:)

    INTEGER :: length
    INTEGER :: i

    length = SIZE(array)
    indexOf = -1
    DO i = 1, length
      IF (array(i) .EQ. key) THEN
        indexOf = i
        EXIT
      END IF
    END DO
  END FUNCTION

  !Returns the index of the target key in the provded array
  !Binary search assumes sorted list
  !Behavior is undefined for arrays with repeated entries of key
  ! -1 means not found
  INTEGER FUNCTION mapIndexOf(key, array)
    INTEGER, INTENT(IN) :: array(:)
    INTEGER, INTENT(IN) :: key

    INTEGER :: start, finish, center, delta, val

    !Assume not found
    mapIndexOf = -1

    !Set initial conditions
    start = 1
    finish = SIZE(array)
    delta = finish - start
    center = (start + finish) / 2
    DO WHILE (delta .GE. 0)
      val = array(center)
      IF (val .EQ. key) THEN
        mapIndexOf = center
        EXIT
      END IF

      IF (val .LT. key) start = center + 1
      IF (val .GT. key) finish = center -1

      delta = finish - start
      center = (start + finish) / 2
    END DO
  END FUNCTION mapIndexOf

END MODULE integer_array_tools
