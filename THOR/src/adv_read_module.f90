
!===============================================================================
!===============================================================================
!===============================================================================
!::A generic input flag parser subroutine is provided at the bottom of the file
!::To insert a new command line flag parser, copy the generic and perform:
!
!1: update the comment header
!2: update the subroutine name and subroutine end statement
!3: Add a description and default example to verbose == 1
!4: Add a description and all possible inputs to verbose == 2
!5: Update invalid parameter error
!6: Provide parsing logic
!7: Provide a default action
!8: Provide a sanity check action (Validator))
!
!To add the parser to the read sequence, add a function call in the main read routine
!===============================================================================
!===============================================================================
!===============================================================================
MODULE adv_read_module
  USE error_generator
  USE global_variables
  USE parameter_types
  IMPLICIT NONE

CONTAINS
  SUBROUTINE adv_read(local_unit)

    !Command Line Arg Parsing
    INTEGER, INTENT(IN) :: local_unit
    !INTEGER :: arg_count = 0
    INTEGER :: indx
    !CHARACTER(50):: arg

    !Control
    !INTEGER :: i = 1
    INTEGER :: verbose=0
    INTEGER :: ioerr=0
    CHARACTER(50) :: in_file = ""
    !CHARACTER(50) :: out_file = ""

    !Input File Parsing
    CHARACTER(200) :: read_str
    CHARACTER(200) :: key_string
    CHARACTER(200) :: data_string
    LOGICAL :: set_defaults
    LOGICAL :: sanity_check
    WRITE(*,*) "UNIT IS: ", local_unit


    ! !Fetch command line argument count
    ! arg_count = COMMAND_ARGUMENT_COUNT()
    !
    ! !Parse args in any order
    ! !NOTE:: Probably will not work if multiple args share a
    ! !       substring - i.e <-i> <-iter> would collide
    ! DO WHILE (i .LE. arg_count)
    !   CALL GETARG(i, arg)
    !
    !   IF(INDEX(arg,'-i') .NE. 0) THEN
    !     i=i+1
    !     CALL GETARG(i, arg)
    !     IF (arg(1:1) .EQ. '-') STOP 'Input File Name Must Follow -i'
    !     in_file = TRIM(arg)
    !
    !   ELSE IF(INDEX(arg,'-o') .NE. 0) THEN
    !     i=i+1
    !     CALL GETARG(i, arg)
    !     IF (arg(1:1) .EQ. '-') THEN
    !       STOP 'Output File Name Must Follow -o'
    !     END IF
    !     out_file = TRIM(arg)
    !
    !   ELSE IF(INDEX(arg,'--verbose=') .NE. 0) THEN
    !     indx = LNBLNK(arg)
    !     read_str = arg(indx : indx)
    !     READ(read_str, *) verbose
    !     IF (verbose .NE. 1 .AND. verbose .NE. 2) THEN
    !       CALL genError(0, "--verbose=n must be specified with n=1 or n=2 ")
    !     END IF
    !   END IF
    !
    !   i=i+1
    ! END DO
    !
    ! !Check for critical arguments
    ! IF (arg_count .EQ. 0 .OR. in_file .EQ. "" .OR. out_file .EQ. "") THEN
    !   WRITE(*,*) "BAD INVOCATION - try: ./a.out -i in.file -o out.file <--verbose=1/2>"
    !   STOP
    ! END IF
    !
    ! !Print out invocation summary
    ! WRITE(*,*) "INPUT SPECIFICATION"
    ! WRITE(*,*) "IN FILE: ", in_file
    ! WRITE(*,*) "OUT FILE: ", out_file
    ! WRITE(*,*) "VERBOSITY: ", verbose
    !
    ! IF (verbose .EQ. 1) THEN
    !   OPEN(UNIT=14, FILE="sample_input.yaml", ACTION="WRITE", STATUS="REPLACE", IOSTAT=ioerr)
    !   IF (ioerr .NE. 0) CALL genError(0, "Error opening verbose outout file")
    !   WRITE(14, *) '#Sample THOR input file with all values set to default'
    !   WRITE(14, *) '#Provide as -i argument to run a THOR'
    !   WRITE(14, *) '#/////////////////////////////////////////////////////'
    !   WRITE(14, *)
    ! END IF
    ! IF (verbose .EQ. 2) THEN
    !   OPEN(UNIT=14, FILE="complete_input.yaml", ACTION="WRITE", STATUS="REPLACE", IOSTAT=ioerr)
    !   IF (ioerr .NE. 0) CALL genError(0, "Error opening verbose outout file")
    !   WRITE(14,*) '#Full THOR input specification.'
    !   WRITE(14,*) '#Select one argument per key to create in input file'
    !   WRITE(14, *) '#/////////////////////////////////////////////////////'
    !   WRITE(14, *)
    ! END IF
    !OPEN(UNIT=local_unit, FILE=in_file, ACTION="READ", STATUS="OLD", IOSTAT=ioerr)
    !IF (ioerr .NE. 0) CALL genError(0, "ERROR READING INPUT FILE")

    !Until file finished
    set_defaults = .TRUE.
    sanity_check = .FALSE.
    DO WHILE (ioerr .EQ. 0 .OR. sanity_check)
      IF (.NOT. set_defaults .AND. .NOT. sanity_check) THEN
        !Consume a line and shift left
        READ(local_unit, '(A)', IOSTAT=ioerr) read_str
        WRITE(*,*) read_str
        IF (ioerr .NE. 0) THEN
          sanity_check = .TRUE.
          CYCLE
        END IF
        read_str = ADJUSTL(read_str)
        !TODO: Force upper || lower case

        !Ignore comments
        IF (TRIM(read_str) .EQ. '---') CYCLE
        IF (read_str(1:1) .EQ. '#') CYCLE

        !Split on ':' delimiter
        indx = INDEX(read_str,':')
        key_string  = read_str(1:indx-1)
        data_string = ADJUSTL(read_str(indx+2:LNBLNK(read_str)))
      END IF

      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "execute"))   &
            CALL input_flag_execute(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "type"))   &
            CALL input_flag_type(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "lambda"))   &
            CALL input_flag_lambda(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "inflow"))   &
            CALL input_flag_inflow(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "piacc"))   &
            CALL input_flag_piacc(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "sweep"))   &
            CALL input_flag_sweep(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "page_sweep"))   &
            CALL input_flag_page_sweep(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "page_refl"))   &
            CALL input_flag_page_refl(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "page_inflow"))   &
            CALL input_flag_page_inflow(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "maxouter"))   &
            CALL input_flag_maxouter(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "maxinner"))   &
            CALL input_flag_maxinner(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "innerconv"))   &
            CALL input_flag_innerconv(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "outerconv"))   &
            CALL input_flag_outerconv(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "kconv"))   &
            CALL input_flag_kconv(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "keigsolve"))   &
            CALL input_flag_keigsolve(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "jfnk_krsze"))   &
            CALL input_flag_jfnk_krsze(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "jfnk_maxkr"))   &
            CALL input_flag_jfnk_maxkr(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "jfnk_method"))   &
            CALL input_flag_jfnk_method(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "initial_guess"))   &
            CALL input_flag_initial_guess(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "save_restart"))   &
            CALL input_flag_save_restart(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "ipiter"))   &
            CALL input_flag_ipiter(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "print_conv"))   &
            CALL input_flag_print_conv(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "density_factor"))   &
            CALL input_flag_density_factor(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "source_file"))   &
            CALL input_flag_source_file(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "inflow_file"))   &
            CALL input_flag_inflow_file(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "xs_file"))   &
            CALL input_flag_xs_file(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "mesh_file"))   &
            CALL input_flag_mesh_file(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "flux_file"))   &
            CALL input_flag_flux_file(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "vtk_flux_file"))   &
            CALL input_flag_vtk_flux_file(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "quad_file"))   &
            CALL input_flag_quad_file(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "restart_file"))   &
            CALL input_flag_restart_file(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "inguess_file"))   &
            CALL input_flag_inguess_file(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "vtk_flux_file"))   &
            CALL input_flag_vtk_flux_file(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "vtk_reg_file"))   &
            CALL input_flag_vtk_reg_file(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "vtk_src_file"))   &
            CALL input_flag_vtk_src_file(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "density_factor_file"))   &
            CALL input_flag_density_factor_file(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "print_xs"))   &
            CALL input_flag_print_xs(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "vtk_flux"))   &
            CALL input_flag_vtk_flux(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "vtk_reg"))   &
            CALL input_flag_vtk_reg(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "vtk_mat"))   &
            CALL input_flag_vtk_mat(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "vtk_src"))   &
            CALL input_flag_vtk_src(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "qdorder"))   &
            CALL input_flag_qdorder(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "qdtype"))   &
            CALL input_flag_qdtype(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "ngroups"))   &
            CALL input_flag_ngroups(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "pnorder"))   &
            CALL input_flag_pnorder(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "pnread"))   &
            CALL input_flag_pnread(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "upscattering"))   &
            CALL input_flag_upscattering(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "multiplying"))   &
            CALL input_flag_multiplying(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "scatt_mult_included"))   &
            CALL input_flag_scatt_mult_included(data_string, set_defaults, sanity_check, verbose)



      set_defaults = .FALSE.
      sanity_check = .FALSE.
    END DO

  END SUBROUTINE adv_read

  !===============================================================================
  ! [Execute] flag !==============================================================
  !===============================================================================
  SUBROUTINE input_flag_execute(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    !DATA_TYPE:: data_value
    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "# Allows for execution of input validation"
        WRITE(14, *) "execution: yes #Perform solve as described in input"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "# Allows for execution of input validation"
        WRITE(14, *) "execution: yes #Perform solve as described in input [DEFAULT]"
        WRITE(14, *) "execution: no  #Stop after parsing input"
        WRITE(14, *)
      END IF
      execution = 1
    ELSE IF(sanity_check) THEN
      IF (execution .NE. 1 .AND. execution .NE. 0) &
            CALL genError(0, '[execution] failed input validation')
    ELSE
      IF (data_string .EQ. "yes") THEN
        execution = 1
      ELSE IF (data_string .EQ. "no") THEN
        execution = 0
      ELSE
        CALL genError(0, "Invalid parameter for [execute] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_execute

  !===============================================================================
  ! [type] flag !=================================================================
  !===============================================================================
  SUBROUTINE input_flag_type(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    !DATA_TYPE:: data_value
    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "# Selects problem type"
        WRITE(14, *) "type: keig"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "# Selects problem type"
        WRITE(14, *) "type: keig #Establishes a k-eigenvalue search problem [DEFAULT]"
        WRITE(14, *) "type: fsrc #Establishes an external source problem"
        WRITE(14, *)
      END IF
      problem = 1
    ELSE IF(sanity_check) THEN
      IF (problem .NE. 1 .AND. problem .NE. 0) &
            CALL genError(0, '[type] failed input validation')
    ELSE
      IF (data_string .EQ. "keig") THEN
        problem = 1
      ELSE IF (data_string .EQ. "fsrc") THEN
        problem = 0
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [type] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_type

  !===============================================================================
  ! [lambda] flag !===============================================================
  !===============================================================================
  SUBROUTINE input_flag_lambda(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose, data_int
    WRITE(*,*) data_string

    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "# Sets the spatial order of the problem. Accepts any integer n>=0"
        WRITE(14, *) "lambda: 0 # Zero order expansion"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "# Sets the spatial order of the problem. Accepts any integer n>=0"
        WRITE(14, *) "lambda: 0 # Zero order expansion [DEFAULT]"
        WRITE(14, *)
      END IF
      space_ord = 0
    ELSE IF(sanity_check) THEN
      IF (space_ord .LT. 0) &
            CALL genError(0, "Invalid parameter for [lambda] input flag")
    ELSE
      READ(data_string, *) data_int
      IF (data_int .GE. 0) THEN
        space_ord = data_int
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [lambda] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_lambda

  !===============================================================================
  ! [inflow] flag !===============================================================
  !===============================================================================
  SUBROUTINE input_flag_inflow(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    !DATA_TYPE:: data_value
    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "# Determines wheter an boundary inflow file is present"
        WRITE(14, *) "inflow: yes"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "# Determines wheter an boundary inflow file is present"
        WRITE(14, *) "inflow: yes # [DEFAULT]"
        WRITE(14, *) "inflow: no"
        WRITE(14, *)
      END IF
      !Set variable to default value here
    ELSE IF(sanity_check) THEN
      !Validate user choice here
    ELSE
      IF (data_string .EQ. "yes") THEN
        !Set variable appropriately
      ELSE IF (data_string .EQ. "no") THEN
        !Set variable appropriately
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [inflow] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_inflow

  !===============================================================================
  ! [piacc] flag !================================================================
  !===============================================================================
  SUBROUTINE input_flag_piacc(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    !DATA_TYPE:: data_value
    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "# Sets the power iteration acceleration mode"
        WRITE(14, *) "piacc: none"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "# Sets the power iteration acceleration mode"
        WRITE(14, *) "piacc: none #no PI acceleration [DEFAULT]"
        WRITE(14, *) "piacc: errmode #Fission Source Extrapolation"
        WRITE(14, *) "piacc: chebychev #Chebychev acceleration"
        WRITE(14, *)
      END IF
      !Set variable to default value here
    ELSE IF(sanity_check) THEN
      !Validate user choice here
    ELSE
      IF (data_string .EQ. "yes") THEN
        !Set variable appropriately
      ELSE IF (data_string .EQ. "no") THEN
        !Set variable appropriately
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [piacc] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_piacc

  !===============================================================================
  ! [sweep] flag !================================================================
  !===============================================================================
  SUBROUTINE input_flag_sweep(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    !DATA_TYPE:: data_value
    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "# Sets whether to use a precomputed sweep order or to calculate on the fly"
        WRITE(14, *) "sweep: precomp #Use a precomputed path"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "# Sets whether to use a precomputed sweep order or to calculate on the fly"
        WRITE(14, *) "sweep: precomp #Use a precomputed sweep path"
        WRITE(14, *)
      END IF
      sweep_tpe = 1
    ELSE IF(sanity_check) THEN
      !Validate user choice here
    ELSE
      IF (data_string .EQ. "precomp") THEN
        sweep_tpe = 1
      ELSE
        CALL genError(0, "Invalid parameter for [sweep] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_sweep

  !===============================================================================
  ! [page_sweep] flag !===========================================================
  !===============================================================================
  SUBROUTINE input_flag_page_sweep(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    !DATA_TYPE:: data_value
    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "# Flush sweep path to file between each inner iteration"
        WRITE(14, *) "page_sweep: yes"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "# Flush sweep path to file between each inner iteration"
        WRITE(14, *) "page_sweep: no # [DEFAULT]"
        WRITE(14, *) "page_sweep: yes "
        WRITE(14, *)
      END IF
      !Set variable to default value here
    ELSE IF(sanity_check) THEN
      !Validate user choice here
    ELSE
      IF (data_string .EQ. "yes") THEN
        !Set variable appropriately
      ELSE IF (data_string .EQ. "no") THEN
        !Set variable appropriately
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [page_sweep] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_page_sweep

  !===============================================================================
  ! [page_refl] flag !============================================================
  !===============================================================================
  SUBROUTINE input_flag_page_refl(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    !DATA_TYPE:: data_value
    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "# Flush reflected bc data to file between each inner iteration"
        WRITE(14, *) "page_refl: save "
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "# Flush reflected bc data to file between each inner iteration"
        WRITE(14, *) "page_refl: save  # Save in memory [DEFAULT]"
        WRITE(14, *) "page_refl: page  # Page out to file"
        WRITE(14, *) "page_refl: inner # Page every inner"
        WRITE(14, *)
      END IF
      page_refl = 0
    ELSE IF(sanity_check) THEN
      !Validate user choice here
    ELSE
      IF (data_string .EQ. "save") THEN
        page_refl = 0
      ELSE IF (data_string .EQ. "page") THEN
        page_refl = 1
      ELSE IF (data_string .EQ. "inner") THEN
        page_refl = 2
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [page_refl] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_page_refl

  !===============================================================================
  ! [page_inflow] flag !==========================================================
  !===============================================================================
  SUBROUTINE input_flag_page_inflow(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    !DATA_TYPE:: data_value
    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "# Flush boundary inflow data to file between each inner iteration"
        WRITE(14, *) "page_inflow: no"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "# Flush boundary inflow data to file between each inner iteration"
        WRITE(14, *) "page_inflow: no # [DEFAULT]"
        WRITE(14, *) "page_inflow: yes"
        WRITE(14, *)
      END IF
      !Set variable to default value here
    ELSE IF(sanity_check) THEN
      !Validate user choice here
    ELSE
      IF (data_string .EQ. "yes") THEN
        !Set variable appropriately
      ELSE IF (data_string .EQ. "no") THEN
        !Set variable appropriately
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [page_inflow] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_page_inflow

  !===============================================================================
  ! [maxouter] flag !=============================================================
  !===============================================================================
  SUBROUTINE input_flag_maxouter(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose, data_int

    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "# Sets the maximum allowed number of outer iterations. Accepts all N > 0"
        WRITE(14, *) "maxouter: 5"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "# Sets the maximum allowed number of outer iterations. Accepts all N > 0"
        WRITE(14, *) "maxouter: 5 #[DEFAULT]"
        WRITE(14, *)
      END IF
      max_outer = 5
    ELSE IF(sanity_check) THEN
      !Validate user choice here
    ELSE
      READ(data_string, *) data_int
      IF (data_int .GT. 0) THEN
        max_outer = data_int
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [maxouter] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_maxouter

  !===============================================================================
  ! [maxinner] flag !=============================================================
  !===============================================================================
  SUBROUTINE input_flag_maxinner(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose, data_int

    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "# Sets the maximum allowed number of inner iterations. Accepts all N > 0"
        WRITE(14, *) "maxinner: 10"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) !Description
        WRITE(14, *) "# Sets the maximum allowed number of inner iterations. Accepts all N > 0"
        WRITE(14, *) "maxinner: 10 #[DEFAULT]"
        WRITE(14, *)
      END IF
      max_inner = 10
    ELSE IF(sanity_check) THEN
      !Validate user choice here
    ELSE
      READ(data_string, *) data_int
      IF (data_int .GT. 0) THEN
        max_inner = data_int
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [maxinner] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_maxinner

  !===============================================================================
  ! [innerconv] flag !============================================================
  !===============================================================================
  SUBROUTINE input_flag_innerconv(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    REAL(kind=d_t) :: data_real
    !DATA_TYPE:: data_value
    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "# Sets the inner iteration convergence threshold. Accepts any real value."
        WRITE(14, *) "innerconv: 1E-4"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "# Sets the inner iteration convergence threshold. Accepts any real value."
        WRITE(14, *) "innerconv: 1E-4 #[DEFAULT]"
        WRITE(14, *)
      END IF
      inner_conv= 1.0E-4_d_t
    ELSE IF(sanity_check) THEN
      !Validate user choice here
    ELSE
      READ(data_string , *) data_real
      IF (data_real .GT. 0) THEN
        inner_conv = data_real
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [innerconv] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_innerconv

  !===============================================================================
  ! [outerconv] flag !============================================================
  !===============================================================================
  SUBROUTINE input_flag_outerconv(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    REAL(kind=d_t) :: data_real

    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "# Sets the outer iteration convergence threshold. Accepts any real value."
        WRITE(14, *) "outerconv: 1E-3"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "# Sets the outer iteration convergence threshold. Accepts any real value."
        WRITE(14, *) "outerconv: 1E-3 #[DEFAULT]"
        WRITE(14, *)
      END IF
      outer_conv= 1.0E-3_d_t
    ELSE IF(sanity_check) THEN
      !Validate user choice here
    ELSE
      READ(data_string , *) data_real
      IF (data_real .GT. 0) THEN
        outer_conv = data_real
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [outerconv] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_outerconv

  !===============================================================================
  ! [kconv] flag !================================================================
  !===============================================================================
  SUBROUTINE input_flag_kconv(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    !DATA_TYPE:: data_value
    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "# Sets the k-eigenvalue iteration convergence threshold. Accepts any real value."
        WRITE(14, *) "kconv: 1E-4"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "# Sets the k-eigenvalue iteration convergence threshold. Accepts any real value."
        WRITE(14, *) "kconv: 1E-4 #[DEFAULT]"
        WRITE(14, *)
      END IF
      !Set variable to default value here
    ELSE IF(sanity_check) THEN
      !Validate user choice here
    ELSE
      IF (data_string .EQ. "yes") THEN
        !Set variable appropriately
      ELSE IF (data_string .EQ. "no") THEN
        !Set variable appropriately
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [kconv] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_kconv

  !===============================================================================
  ! [keigsolve] flag !============================================================
  !===============================================================================
  SUBROUTINE input_flag_keigsolve(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    !DATA_TYPE:: data_value
    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "# Selects the method by which to solve k-eigenvalue searches"
        WRITE(14, *) "keigsolver: pi #Power Iteration mode"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "# Selects the method by which to solve k-eigenvalue searches"
        WRITE(14, *) "keigsolver: pi #Power Iteration mode [DEFAULT]"
        WRITE(14, *) "keigsolver: jfnk #Jacobi-Free Newton Krylov mode"
        WRITE(14, *)
      END IF
      !Set variable to default value here
    ELSE IF(sanity_check) THEN
      !Validate user choice here
    ELSE
      IF (data_string .EQ. "yes") THEN
        !Set variable appropriately
      ELSE IF (data_string .EQ. "no") THEN
        !Set variable appropriately
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [keigsolver] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_keigsolve

  !===============================================================================
  ! [jfnk_krsze] flag !===========================================================
  !===============================================================================
  SUBROUTINE input_flag_jfnk_krsze(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    !DATA_TYPE:: data_value
    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "# The maximum number of solution vectors to be stored within the krylov system. Accepts integers > 0"
        WRITE(14, *) "jfnk_krsze: 25"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "# The maximum number of solution vectors to be stored within the krylov system. Accepts integers > 0"
        WRITE(14, *) "jfnk_krsze: 25 # [DEFAULT]"
        WRITE(14, *)
      END IF
      !Set variable to default value here
    ELSE IF(sanity_check) THEN
      !Validate user choice here
    ELSE
      IF (data_string .EQ. "yes") THEN
        !Set variable appropriately
      ELSE IF (data_string .EQ. "no") THEN
        !Set variable appropriately
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [jfnk_krsze] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_jfnk_krsze

  !===============================================================================
  ! [jfnk_maxkr] flag !===========================================================
  !===============================================================================
  SUBROUTINE input_flag_jfnk_maxkr(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    !DATA_TYPE:: data_value
    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "The maximum number of krylov iterations per non-linear step. Accepts integers > 0"
        WRITE(14, *) "jfnk_maxkr: 250"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "The maximum number of krylov iterations per non-linear step. Accepts integers > 0"
        WRITE(14, *) "jfnk_maxkr: 250 # [DEFAULT]"
        WRITE(14, *)
      END IF
      !Set variable to default value here
    ELSE IF(sanity_check) THEN
      !Validate user choice here
    ELSE
      IF (data_string .EQ. "yes") THEN
        !Set variable appropriately
      ELSE IF (data_string .EQ. "no") THEN
        !Set variable appropriately
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [jfnk_maxkr] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_jfnk_maxkr

  !===============================================================================
  ! [jfnk_method] flag !==========================================================
  !===============================================================================
  SUBROUTINE input_flag_jfnk_method(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    !DATA_TYPE:: data_value
    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "#Choose JFNK solver formulation"
        WRITE(14, *) "jfnk_method: flat # non-linear function evaluation is a sweep on constructed source"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "#Choose JFNK solver formulation"
        WRITE(14, *) "jfnk_method: flat # non-linear function evaluation is a sweep on constructed source [DEFAULT]"
        WRITE(14, *) "jfnk_method: outer # non-linear function evaluation is a single outer iterations"
        WRITE(14, *) "jfnk_method: flat_wds # non-linear function evaluation is a sweep on constructed" &
              // " source with updated downscattering"
        WRITE(14, *)
      END IF
      !Set variable to default value here
    ELSE IF(sanity_check) THEN
      !Validate user choice here
    ELSE
      IF (data_string .EQ. "yes") THEN
        !Set variable appropriately
      ELSE IF (data_string .EQ. "no") THEN
        !Set variable appropriately
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [jfnk_method] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_jfnk_method

  !===============================================================================
  ! [initial_guess] flag !========================================================
  !===============================================================================
  SUBROUTINE input_flag_initial_guess(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    !DATA_TYPE:: data_value
    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "# Specifies whether an initial guess file is provided"
        WRITE(14, *) "initial_guess: no #No initial guess file provided"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "# Specifies whether an initial guess file is provided"
        WRITE(14, *) "initial_guess: no #[#DEFAULT]"
        WRITE(14, *) "initial_guess: yes"
        WRITE(14, *)
      END IF
      !Set variable to default value here
    ELSE IF(sanity_check) THEN
      !Validate user choice here
    ELSE
      IF (data_string .EQ. "yes") THEN
        !Set variable appropriately
      ELSE IF (data_string .EQ. "no") THEN
        !Set variable appropriately
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [initial_guess] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_initial_guess

  !===============================================================================
  ! [save_restart] flag !=========================================================
  !===============================================================================
  SUBROUTINE input_flag_save_restart(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    !DATA_TYPE:: data_value
    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "# Determines whether a restart data file should be output after each outer iteration"
        WRITE(14, *) "save_restart: yes"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "# Determines whether a restart data file should be output after each outer iteration"
        WRITE(14, *) "save_restart: yes #[DEFAULT]"
        WRITE(14, *) "save_restart: no"
        WRITE(14, *)
      END IF
      !Set variable to default value here
    ELSE IF(sanity_check) THEN
      !Validate user choice here
    ELSE
      IF (data_string .EQ. "yes") THEN
        !Set variable appropriately
      ELSE IF (data_string .EQ. "no") THEN
        !Set variable appropriately
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [save_restart] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_save_restart

  !===============================================================================
  ! [ipiter] flag !===============================================================
  !===============================================================================
  SUBROUTINE input_flag_ipiter(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    !DATA_TYPE:: data_value
    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "# Defines the number of initial power iterations during a jfnk solve. Accepts integers >= 0"
        WRITE(14, *) "ipiter: 0"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "# Defines the number of initial power iterations during a jfnk solve. Accepts integers >= 0"
        WRITE(14, *) "ipiter: 0 # [DEFAULT]"
        WRITE(14, *)
      END IF
      !Set variable to default value here
    ELSE IF(sanity_check) THEN
      !Validate user choice here
    ELSE
      IF (data_string .EQ. "yes") THEN
        !Set variable appropriately
      ELSE IF (data_string .EQ. "no") THEN
        !Set variable appropriately
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [ipiter] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_ipiter

  !===============================================================================
  ! [print_conv] flag !===========================================================
  !===============================================================================
  SUBROUTINE input_flag_print_conv(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    !DATA_TYPE:: data_value
    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "# Determines whether to print a running convergence history to file during execution"
        WRITE(14, *) "print_conv: no"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "# Determines whether to print a running convergence history to file during execution"
        WRITE(14, *) "print_conv: no #[DEFAULT]"
        WRITE(14, *) "print_conv: yes"
        WRITE(14, *)
      END IF
      !Set variable to default value here
    ELSE IF(sanity_check) THEN
      !Validate user choice here
    ELSE
      IF (data_string .EQ. "yes") THEN
        !Set variable appropriately
      ELSE IF (data_string .EQ. "no") THEN
        !Set variable appropriately
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [generic] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_print_conv

  !===============================================================================
  ! [density_factor] flag !=======================================================
  !===============================================================================
  SUBROUTINE input_flag_density_factor(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    !DATA_TYPE:: data_value
    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) '# Determines whether density factors should be applied'
        WRITE(14, *) "density_factor: no"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) '# Determines whether density factors should be applied'
        WRITE(14, *) "density_factor: no #[DEFAULT]"
        WRITE(14, *) "density_factor: byvolume #Assign density factors by region volume"
        WRITE(14, *) "density_factor: fromfile #Assign density factors based on input file"
        WRITE(14, *)
      END IF
      !Set variable to default value here
    ELSE IF(sanity_check) THEN
      !Validate user choice here
    ELSE
      IF (data_string .EQ. "yes") THEN
        !Set variable appropriately
      ELSE IF (data_string .EQ. "no") THEN
        !Set variable appropriately
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [density_factor] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_density_factor

  !===============================================================================
  ! [source_file] flag !==========================================================
  !===============================================================================
  SUBROUTINE input_flag_source_file(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    LOGICAL :: data_logical

    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "# Specifies the name of the source input file. Accepts a string of length < 80"
        WRITE(14, *) "source_file: file.src"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "# Specifies the name of the source input file. Accepts a string of length < 80"
        WRITE(14, *) "source_file: file.src #[DEFAULT]"
        WRITE(14, *)
      END IF
      source_filename = "file.src"
    ELSE IF(sanity_check) THEN
      !Validate user choice here
    ELSE
      INQUIRE(FILE = data_string, EXIST = data_logical)
      IF(data_logical) THEN
        source_filename = data_string
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [source_file] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_source_file

  !===============================================================================
  ! [inflow_file] flag !==========================================================
  !===============================================================================
  SUBROUTINE input_flag_inflow_file(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    !DATA_TYPE:: data_value
    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "# Defines the filename for the boundary inflow input file"
        WRITE(14, *) "inflow_file: file.bc"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "# Defines the filename for the boundary inflow input file"
        WRITE(14, *) "inflow_file: file.bc # [DEFAULT]"
        WRITE(14, *)
      END IF
      !Set variable to default value here
    ELSE IF(sanity_check) THEN
      !Validate user choice here
    ELSE
      IF (data_string .EQ. "yes") THEN
        !Set variable appropriately
      ELSE IF (data_string .EQ. "no") THEN
        !Set variable appropriately
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [inflow_file] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_inflow_file

  !===============================================================================
  ! [xs_file] flag !==============================================================
  !===============================================================================
  SUBROUTINE input_flag_xs_file(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    LOGICAL :: data_logical
    !DATA_TYPE:: data_value
    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "# Specifies the name of the cross-section input file. Accepts any string of length < 80"
        WRITE(14, *) "xs_file: file.xs"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "# Specifies the name of the cross-section input file. Accepts any string of length < 80"
        WRITE(14, *) "xs_file: file.xs #[DEFAULT]"
        WRITE(14, *)
      END IF
      cross_section_filename = "file.xs"
    ELSE IF(sanity_check) THEN
      !Validate user choice here
    ELSE
      INQUIRE(FILE = data_string, EXIST = data_logical)
      IF(data_logical) THEN
        cross_section_filename = data_string
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [xs_file] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_xs_file

  !===============================================================================
  ! [mesh_file] flag !============================================================
  !===============================================================================
  SUBROUTINE input_flag_mesh_file(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    LOGICAL :: data_logical

    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "# Specifies the name of the mesh input file. Accepts any string of length < 80"
        WRITE(14, *) "mesh_file: file.thrm"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "# Specifies the name of the mesh input file. Accepts any string of length < 80"
        WRITE(14, *) "mesh_file: file.xs #[DEFAULT]"
        WRITE(14, *)
      END IF
      mesh_filename = "file.mesh"
    ELSE IF(sanity_check) THEN
      !Validate user choice here
    ELSE
      INQUIRE(FILE = data_string, EXIST = data_logical)
      IF(data_logical) THEN
        mesh_filename = data_string
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [mesh_file] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_mesh_file

  !===============================================================================
  ! [flux_file] flag !============================================================
  !===============================================================================
  SUBROUTINE input_flag_flux_file(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    LOGICAL :: data_logical

    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "#File name for scalar flux output dump"
        WRITE(14, *) "flux_file: file.flux"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "#File name for flux output dump"
        WRITE(14, *) "flux_file: file.flux # [DEFAULT]"
        WRITE(14, *)
      END IF
      flux_filename = "file.flux"
    ELSE IF(sanity_check) THEN
      !Validate user choice here
    ELSE
      INQUIRE(FILE = data_string, EXIST = data_logical)
      IF(data_logical) THEN
        flux_filename = data_string
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [flux_file] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_flux_file

  !===============================================================================
  ! [vtk_flux_file] flag !========================================================
  !===============================================================================
  SUBROUTINE input_flag_vtk_flux_file(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    !DATA_TYPE:: data_value
    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "#Filename for the .vtk formatted flux output file"
        WRITE(14, *) "vtk_flux_file: flux.vtk"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "#Filename for the .vtk formatted flux output file"
        WRITE(14, *) "vtk_flux_file: flux.vtk # [DEFAULT]"
        WRITE(14, *)
      END IF
      !Set variable to default value here
    ELSE IF(sanity_check) THEN
      !Validate user choice here
    ELSE
      IF (data_string .EQ. "yes") THEN
        !Set variable appropriately
      ELSE IF (data_string .EQ. "no") THEN
        !Set variable appropriately
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [vtk_flux_file] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_vtk_flux_file

  !===============================================================================
  ! [quad_file] flag !============================================================
  !===============================================================================
  SUBROUTINE input_flag_quad_file(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    !DATA_TYPE:: data_value
    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "# Specifies the name of the quadrature input file. Accepts any string of length < 80"
        WRITE(14, *) "quad_file: file.quad"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "# Specifies the name of the quadrature input file. Accepts any string of length < 80"
        WRITE(14, *) "quad_file: file.quad #[DEFAULT]"
        WRITE(14, *)
      END IF
      !Set variable to default value here
    ELSE IF(sanity_check) THEN
      !Validate user choice here
    ELSE
      IF (data_string .EQ. "yes") THEN
        !Set variable appropriately
      ELSE IF (data_string .EQ. "no") THEN
        !Set variable appropriately
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [quad_file] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_quad_file

  !===============================================================================
  ! [restart_file] flag !=========================================================
  !===============================================================================
  SUBROUTINE input_flag_restart_file(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    !DATA_TYPE:: data_value
    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "# Specifies the name of the restart file. Accepts any string of length < 80"
        WRITE(14, *) "restart_file: file.restart"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "# Specifies the name of the restart file. Accepts any string of length < 80"
        WRITE(14, *) "restart_file: file.restart"
        WRITE(14, *)
      END IF
      !Set variable to default value here
    ELSE IF(sanity_check) THEN
      !Validate user choice here
    ELSE
      IF (data_string .EQ. "yes") THEN
        !Set variable appropriately
      ELSE IF (data_string .EQ. "no") THEN
        !Set variable appropriately
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [restart_file] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_restart_file

  !===============================================================================
  ! [inguess_file] flag !=========================================================
  !===============================================================================
  SUBROUTINE input_flag_inguess_file(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    !DATA_TYPE:: data_value
    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "# Filename for initial guess file (also accepts restart files)"
        WRITE(14, *) "inguess_file: inguess"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "# Filename for initial guess file (also accepts restart files)"
        WRITE(14, *) "inguess_file: inguess # [DEFAULT]"
        WRITE(14, *)
      END IF
      !Set variable to default value here
    ELSE IF(sanity_check) THEN
      !Validate user choice here
    ELSE
      IF (data_string .EQ. "yes") THEN
        !Set variable appropriately
      ELSE IF (data_string .EQ. "no") THEN
        !Set variable appropriately
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [inguess_file] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_inguess_file

  !===============================================================================
  ! [vtk_mat_file] flag !=========================================================
  !===============================================================================
  SUBROUTINE input_flag_vtk_mat_file(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    !DATA_TYPE:: data_value
    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "#Filename for the .vtk formatted material map output file"
        WRITE(14, *) "vtk_mat_file: mat.vtk "
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "#Filename for the .vtk formatted flux output file"
        WRITE(14, *) "vtk_mat_file: mat.vtk # [DEFAULT]"
        WRITE(14, *)
      END IF
      !Set variable to default value here
    ELSE IF(sanity_check) THEN
      !Validate user choice here
    ELSE
      IF (data_string .EQ. "yes") THEN
        !Set variable appropriately
      ELSE IF (data_string .EQ. "no") THEN
        !Set variable appropriately
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [vtk_mat_file] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_vtk_mat_file

  !===============================================================================
  ! [vtk_reg_file] flag !=========================================================
  !===============================================================================
  SUBROUTINE input_flag_vtk_reg_file(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    !DATA_TYPE:: data_value
    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "#Filename for the .vtk formatted region specification output file"
        WRITE(14, *) "vtk_reg_file: reg.vtk"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "#Filename for the .vtk formatted region specification output file"
        WRITE(14, *) "vtk_reg_file: reg.vtk # [DEFAULT]"
        WRITE(14, *)
      END IF
      !Set variable to default value here
    ELSE IF(sanity_check) THEN
      !Validate user choice here
    ELSE
      IF (data_string .EQ. "yes") THEN
        !Set variable appropriately
      ELSE IF (data_string .EQ. "no") THEN
        !Set variable appropriately
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [vtk_reg_file] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_vtk_reg_file

  !===============================================================================
  ! [vtk_src_file] flag !=========================================================
  !===============================================================================
  SUBROUTINE input_flag_vtk_src_file(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    !DATA_TYPE:: data_value
    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "#Filename for the .vtk formatted source specification output file"
        WRITE(14, *) "vtk_src_file: src.vtk"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "#Filename for the .vtk formatted source specification output file"
        WRITE(14, *) "vtk_src_file: src.vtk"
        WRITE(14, *)
      END IF
      !Set variable to default value here
    ELSE IF(sanity_check) THEN
      !Validate user choice here
    ELSE
      IF (data_string .EQ. "yes") THEN
        !Set variable appropriately
      ELSE IF (data_string .EQ. "no") THEN
        !Set variable appropriately
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [vtk_src_file] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_vtk_src_file

  !===============================================================================
  ! [density_factor_file] flag !==================================================
  !===============================================================================
  SUBROUTINE input_flag_density_factor_file(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    !DATA_TYPE:: data_value
    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "# Specifies the name of the density factor input file. Accepts any string of length < 80"
        WRITE(14, *) "density_factor_file: density_factors.dat"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "# Specifies the name of the density factor input file. Accepts any string of length < 80"
        WRITE(14, *) "density_factor_file: density_factors.dat #[DEFAULT]"
        WRITE(14, *)
      END IF
      !Set variable to default value here
    ELSE IF(sanity_check) THEN
      !Validate user choice here
    ELSE
      IF (data_string .EQ. "yes") THEN
        !Set variable appropriately
      ELSE IF (data_string .EQ. "no") THEN
        !Set variable appropriately
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [generic] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_density_factor_file

  !===============================================================================
  ! [print_xs] flag !=============================================================
  !===============================================================================
  SUBROUTINE input_flag_print_xs(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    !DATA_TYPE:: data_value
    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "# Determine whether to print a cross-section summary before execution"
        WRITE(14, *) "print_xs: yes"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "# Determine whether to print a cross-section summary before execution"
        WRITE(14, *) "print_xs: yes #[DEFAULT]"
        WRITE(14, *) "print_xs: no"
        WRITE(14, *)
      END IF
      print_xs_flag = 1
    ELSE IF(sanity_check) THEN
      !Validate user choice here
    ELSE
      IF (data_string .EQ. "yes") THEN
        print_xs_flag = 1
      ELSE IF (data_string .EQ. "no") THEN
        print_xs_flag = 0
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [print_xs] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_print_xs

  !===============================================================================
  ! [vtk_flux] flag !=============================================================
  !===============================================================================
  SUBROUTINE input_flag_vtk_flux(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    !DATA_TYPE:: data_value
    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "# Determines wheter to output a .vtk formatted flux file after execution"
        WRITE(14, *) "vtk_flux: no"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "# Determines wheter to output a .vtk formatted flux file after execution"
        WRITE(14, *) "vtk_flux: no # [DEFAULT]"
        WRITE(14, *) "vtk_flux: yes"
      END IF
      !Set variable to default value here
    ELSE IF(sanity_check) THEN
      !Validate user choice here
    ELSE
      IF (data_string .EQ. "yes") THEN
        !Set variable appropriately
      ELSE IF (data_string .EQ. "no") THEN
        !Set variable appropriately
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [vtk_flux] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_vtk_flux

  !===============================================================================
  ! [vtk_reg] flag !==============================================================
  !===============================================================================
  SUBROUTINE input_flag_vtk_reg(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    !DATA_TYPE:: data_value
    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "# Determines wheter to output a .vtk formatted region mapping file after execution"
        WRITE(14, *) "vtk_reg: no"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "# Determines wheter to output a .vtk formatted region mapping file after execution"
        WRITE(14, *) "vtk_reg: no # [DEFAULT]"
        WRITE(14, *) "vtk_reg: yes"
        WRITE(14, *)
      END IF
      !Set variable to default value here
    ELSE IF(sanity_check) THEN
      !Validate user choice here
    ELSE
      IF (data_string .EQ. "yes") THEN
        !Set variable appropriately
      ELSE IF (data_string .EQ. "no") THEN
        !Set variable appropriately
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [vtk_reg] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_vtk_reg

  !===============================================================================
  ! [vtk_mat] flag !==============================================================
  !===============================================================================
  SUBROUTINE input_flag_vtk_mat(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    !DATA_TYPE:: data_value
    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "# Determines wheter to output a .vtk formatted material mapping file after execution"
        WRITE(14, *) "vtk_mat: no"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "# Determines wheter to output a .vtk formatted material mapping file after execution"
        WRITE(14, *) "vtk_mat: no # [DEFAULT]"
        WRITE(14, *) "vtk_mat: yes"
        WRITE(14, *)
      END IF
      !Set variable to default value here
    ELSE IF(sanity_check) THEN
      !Validate user choice here
    ELSE
      IF (data_string .EQ. "yes") THEN
        !Set variable appropriately
      ELSE IF (data_string .EQ. "no") THEN
        !Set variable appropriately
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [vtk_mat] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_vtk_mat

  !===============================================================================
  ! [vtk_src] flag !==============================================================
  !===============================================================================
  SUBROUTINE input_flag_vtk_src(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    !DATA_TYPE:: data_value
    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "# Determines wheter to output a .vtk formatted source mapping file after execution"
        WRITE(14, *) "vtk_src: no"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "# Determines wheter to output a .vtk formatted source mapping file after execution"
        WRITE(14, *) "vtk_src: no # [DEFAULT]"
        WRITE(14, *) "vtk_src: yes"
        WRITE(14, *)
      END IF
      !Set variable to default value here
    ELSE IF(sanity_check) THEN
      !Validate user choice here
    ELSE
      IF (data_string .EQ. "yes") THEN
        !Set variable appropriately
      ELSE IF (data_string .EQ. "no") THEN
        !Set variable appropriately
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [vtk_src] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_vtk_src

  !===============================================================================
  ! [qdorder] flag !==============================================================
  !===============================================================================
  SUBROUTINE input_flag_qdorder(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose, data_int

    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "# Sets the quadrature order. Accepts any even integer >= 2"
        WRITE(14, *) "qdorder: 4"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "# Sets the quadrature order. Accepts any even integer >= 2"
        WRITE(14, *) "qdorder: 4 #[DEFAULT]"
        WRITE(14, *)
      END IF
      !Set variable to default value here
    ELSE IF(sanity_check) THEN
      !Validate user choice here
    ELSE
      READ(data_string, *) data_int
      IF (data_int .GE. 2) THEN
        quad_ord = data_int
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [qdorder] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_qdorder

  !===============================================================================
  ! [qdtype] flag !===============================================================
  !===============================================================================
  SUBROUTINE input_flag_qdtype(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    !DATA_TYPE:: data_value
    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "# Determines the quadrature type"
        WRITE(14, *) "qdtype: levelsym # Level-symmetric quadrature set"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "# Determines the quadrature type"
        WRITE(14, *) "qdtype: levelsym # Level-symmetric quadrature set [DEFAULT]"
        WRITE(14, *) "qdtype: legcheb # Legendre-Chebychev quadrature set"
        WRITE(14, *) "qdtype: fromfile # read in quadrature set from file"
        WRITE(14, *)
      END IF
      quad_tpe = 1
    ELSE IF(sanity_check) THEN
      !Validate user choice here
    ELSE
      IF (data_string .EQ. "levelsym") THEN
        quad_tpe = 1
      ELSE IF (data_string .EQ. "legcheb") THEN
        quad_tpe = 2
      ELSE IF (data_string .EQ. "fromfile") THEN
        quad_tpe = 3
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [qdtype] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_qdtype

  !===============================================================================
  ! [ngroups] flag !==============================================================
  !===============================================================================
  SUBROUTINE input_flag_ngroups(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose, data_int

    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "# Sets the number of energy groups. Accepts integers > 0"
        WRITE(14, *) "ngroups: 1"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "# Sets the number of energy groups. Accepts integers > 0"
        WRITE(14, *) "ngroups: 1 #[DEFAULT]"
        WRITE(14, *)
      END IF
      egmax = 1
    ELSE IF(sanity_check) THEN
      !Validate user choice here
    ELSE
      READ(data_string, *) data_int
      IF (data_int .GT. 0) THEN
        egmax = data_int
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [ngroups] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_ngroups

  !===============================================================================
  ! [pnorder] flag !==============================================================
  !===============================================================================
  SUBROUTINE input_flag_pnorder(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    !DATA_TYPE:: data_value
    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "# Sets the order of the scattering term expansion. Accepts integers >= 0"
        WRITE(14, *) "pnorder: 0"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "# Sets the order of the scattering term expansion. Accepts integers >= 0"
        WRITE(14, *) "pnorder: 0 #[DEFAULT]"
        WRITE(14, *)
      END IF
      !Set variable to default value here
    ELSE IF(sanity_check) THEN
      !Validate user choice here
    ELSE
      IF (data_string .EQ. "yes") THEN
        !Set variable appropriately
      ELSE IF (data_string .EQ. "no") THEN
        !Set variable appropriately
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [pnorder] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_pnorder

  !===============================================================================
  ! [pnread] flag !===============================================================
  !===============================================================================
  SUBROUTINE input_flag_pnread(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    !DATA_TYPE:: data_value
    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "# Sets the order of the cross-section term expansion. Accepts integers >= 0"
        WRITE(14, *) "pnread: 0 #[DEFAULT]"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "# Sets the order of the cross-section term expansion. Accepts integers >= 0"
        WRITE(14, *) "pnread: 0 #[DEFAULT]"
        WRITE(14, *)
      END IF
      !Set variable to default value here
    ELSE IF(sanity_check) THEN
      !Validate user choice here
    ELSE
      IF (data_string .EQ. "yes") THEN
        !Set variable appropriately
      ELSE IF (data_string .EQ. "no") THEN
        !Set variable appropriately
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [pnread] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_pnread

  !===============================================================================
  ! [upscattering] flag !=========================================================
  !===============================================================================
  SUBROUTINE input_flag_upscattering(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    !DATA_TYPE:: data_value
    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "# Determines whether THOR will calculate the upscattering contribution"
        WRITE(14, *) "upscattering: yes"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "# Determines whether THOR will calculate the upscattering contribution"
        WRITE(14, *) "upscattering: yes #[DEFAULT]"
        WRITE(14, *) "upscattering: no"
        WRITE(14, *)
      END IF
      upscattering = 1
    ELSE IF(sanity_check) THEN
      !Validate user choice here
    ELSE
      IF (data_string .EQ. "yes") THEN
        upscattering = 1
      ELSE IF (data_string .EQ. "no") THEN
        upscattering = 0
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [upscattering] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_upscattering

  !===============================================================================
  ! [multiplying] flag !==========================================================
  !===============================================================================
  SUBROUTINE input_flag_multiplying(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    !DATA_TYPE:: data_value
    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) '# Determines whether THOR will calculate the fission multiplication contribution'
        WRITE(14, *) 'multiplying: yes'
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) '# Determines whether THOR will calculate the fission multiplication contribution'
        WRITE(14, *) 'multiplying: yes #[DEFAULT]'
        WRITE(14, *) 'multiplying: no'
        WRITE(14, *)
      END IF
      multiplying  = 1
    ELSE IF(sanity_check) THEN
      !Validate user choice here
    ELSE
      IF (data_string .EQ. "yes") THEN
        multiplying  = 1
      ELSE IF (data_string .EQ. "no") THEN
        multiplying  = 0
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [multiplying] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_multiplying

  !===============================================================================
  ! [scatt_mult_included] flag !==================================================
  !===============================================================================
  SUBROUTINE input_flag_scatt_mult_included(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    !DATA_TYPE:: data_value
    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "Do the provided scattering cross-section moments include a (2*l+1) factor?"
        WRITE(14, *) "scatt_mult_included: no"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "Do the provided scattering cross-section moments include a (2*l+1) factor?"
        WRITE(14, *) "scatt_mult_included: no # [DEFAULT]"
        WRITE(14, *) "scatt_mult_included: yes"
        WRITE(14, *)
      END IF
      !Set variable to default value here
    ELSE IF(sanity_check) THEN
      !Validate user choice here
    ELSE
      IF (data_string .EQ. "yes") THEN
        !Set variable appropriately
      ELSE IF (data_string .EQ. "no") THEN
        !Set variable appropriately
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [scatt_mult_included] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_scatt_mult_included

  !===============================================================================
  ! [GENERIC] flag !==============================================================
  !===============================================================================
  SUBROUTINE input_flag_generic(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    !DATA_TYPE:: data_value
    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) !Description
        WRITE(14, *) !Default Value Description
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) !Description
        WRITE(14, *) !Default Value Description and [DEFAULT]
        WRITE(14, *) !Other option Description
        WRITE(14, *)
      END IF
      !Set variable to default value here
    ELSE IF(sanity_check) THEN
      !Validate user choice here
    ELSE
      IF (data_string .EQ. "yes") THEN
        !Set variable appropriately
      ELSE IF (data_string .EQ. "no") THEN
        !Set variable appropriately
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [generic] input flag")
      END IF
    END IF
  END SUBROUTINE input_flag_generic

END MODULE adv_read_module
