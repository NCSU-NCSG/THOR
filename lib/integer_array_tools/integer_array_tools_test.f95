PROGRAM TEST
  USE integer_array_tools
  IMPLICIT NONE
  INTEGER :: A(5) = (/5, 3, 2, 1, 4/), i
  INTEGER :: P(5) = (/1, 2, 3, 4, 5/), P2(10)=0

  INTEGER :: B(5) = (/1, 3, 4, 5, 2/)

  INTEGER, ALLOCATABLE:: C(:), D(:)

  CALL quickSortInteger(A,P)
  WRITE(*,'(A,5(I0,X))') 'Expected: ', 1 , 2, 3 ,4 ,5
  WRITE(*, '(A,5(I0,X))')  'Recieved: ', A

  P = (/1,2,3,4,5/)
  A = (/5, 3, 2, 1, 4/)

  CALL quickSortInteger(A, P, B)
  WRITE(*,*)
  WRITE(*, '(A,3(5(I0,X), 3X))') 'Expected: ', 1,2,3,4,5,4,3,2,5,1,5,4,3,2,1
  WRITE(*, '(A,3(5(I0,X), 3X))')  'Recieved: ', A, P,B


  ALLOCATE(C(10), D(10))
  DO i = 1, 10
    D(i) = i
    P2(i) = i
    C(11-i) = i
  END DO
  CALL quickSortInteger(C, P2, D)

  WRITE(*, '(A,3(10(I0,X), 3X))') 'Expected: ', 1,2,3,4,5,6,7,8,9,10,10,9,8,7,6,5,4,3,2,1,10,9,8,7,6,5,4,3,2,1
  WRITE(*, '(A,3(10(I0,X), 3X))')  'Recieved: ', C, P2 ,D

  WRITE(*,*) 'Expected: ', 3
  WRITE(*,*) 'Recieved: ',mapIndexOf(3, C)
  WRITE(*,*) 'Expected: ', 7
  WRITE(*,*) 'Recieved: ', mapIndexOf(7, C)
  WRITE(*,*) 'Expected: ', 1
  WRITE(*,*) 'Recieved: ', mapIndexOf(1, C)
  WRITE(*,*) 'Expected: ', 10
  WRITE(*,*) 'Recieved: ', mapIndexOf(10, C)
  WRITE(*,*) 'Expected: ', -1
  WRITE(*,*) 'Recieved: ', mapIndexOf(0, C)
  WRITE(*,*) 'Expected: ', -1
  WRITE(*,*) 'Recieved: ', mapIndexOf(11, C)


END PROGRAM
