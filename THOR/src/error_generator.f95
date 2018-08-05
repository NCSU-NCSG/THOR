MODULE ERROR_GENERATOR
CONTAINS
  SUBROUTINE genError(unit, message, err_code)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: unit
    INTEGER, INTENT(IN), OPTIONAL:: err_code
    CHARACTER(*), INTENT(IN):: message

    WRITE(unit, '(A)') message
    IF (PRESENT(err_code)) WRITE(unit, '(A,I0)') "Error Code: ", err_code

    STOP

  END SUBROUTINE genError
END MODULE ERROR_GENERATOR
