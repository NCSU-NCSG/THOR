!***********************************************************************
!
! The module contains all legacy subroutines for old yaml_input versions. These
! subroutines should generally not be changed to garauntee backwards
! compatibility.
!
!***********************************************************************
MODULE read_module_legacy
  USE mpi
  USE stringmod
  USE error_generator
  USE global_variables
  USE parameter_types
  IMPLICIT NONE

CONTAINS
  SUBROUTINE yaml_read(local_unit)

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
    !   WRITE(14, *) '#Sample THOR yaml_input file with all values set to default'
    !   WRITE(14, *) '#Provide as -i argument to run a THOR'
    !   WRITE(14, *) '#/////////////////////////////////////////////////////'
    !   WRITE(14, *)
    ! END IF
    ! IF (verbose .EQ. 2) THEN
    !   OPEN(UNIT=14, FILE="complete_input.yaml", ACTION="WRITE", STATUS="REPLACE", IOSTAT=ioerr)
    !   IF (ioerr .NE. 0) CALL genError(0, "Error opening verbose outout file")
    !   WRITE(14,*) '#Full THOR yaml_input specification.'
    !   WRITE(14,*) '#Select one argument per key to create in yaml_input file'
    !   WRITE(14, *) '#/////////////////////////////////////////////////////'
    !   WRITE(14, *)
    ! END IF
    !OPEN(UNIT=local_unit, FILE=in_file, ACTION="READ", STATUS="OLD", IOSTAT=ioerr)
    !IF (ioerr .NE. 0) CALL genError(0, "ERROR READING yaml_input FILE")

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
            CALL yaml_input_flag_execute(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "type"))   &
            CALL yaml_input_flag_type(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "lambda"))   &
            CALL yaml_input_flag_lambda(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "inflow"))   &
            CALL yaml_input_flag_inflow(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "piacc"))   &
            CALL yaml_input_flag_piacc(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "sweep"))   &
            CALL yaml_input_flag_sweep(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "page_sweep"))   &
            CALL yaml_input_flag_page_sweep(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "page_refl"))   &
            CALL yaml_input_flag_page_refl(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "page_inflow"))   &
            CALL yaml_input_flag_page_inflow(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "maxouter"))   &
            CALL yaml_input_flag_maxouter(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "maxinner"))   &
            CALL yaml_input_flag_maxinner(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "innerconv"))   &
            CALL yaml_input_flag_innerconv(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "outerconv"))   &
            CALL yaml_input_flag_outerconv(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "kconv"))   &
            CALL yaml_input_flag_kconv(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "keigsolve"))   &
            CALL yaml_input_flag_keigsolve(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "jfnk_krsze"))   &
            CALL yaml_input_flag_jfnk_krsze(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "jfnk_maxkr"))   &
            CALL yaml_input_flag_jfnk_maxkr(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "jfnk_method"))   &
            CALL yaml_input_flag_jfnk_method(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "initial_guess"))   &
            CALL yaml_input_flag_initial_guess(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "save_restart"))   &
            CALL yaml_input_flag_save_restart(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "ipiter"))   &
            CALL yaml_input_flag_ipiter(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "print_conv"))   &
            CALL yaml_input_flag_print_conv(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "density_factor"))   &
            CALL yaml_input_flag_density_factor(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "source_file"))   &
            CALL yaml_input_flag_source_file(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "inflow_file"))   &
            CALL yaml_input_flag_inflow_file(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "xs_file"))   &
            CALL yaml_input_flag_xs_file(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "mesh_file"))   &
            CALL yaml_input_flag_mesh_file(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "flux_file"))   &
            CALL yaml_input_flag_flux_file(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "vtk_flux_file"))   &
            CALL yaml_input_flag_vtk_flux_file(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "quad_file"))   &
            CALL yaml_input_flag_quad_file(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "restart_file"))   &
            CALL yaml_input_flag_restart_file(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "inguess_file"))   &
            CALL yaml_input_flag_inguess_file(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "vtk_flux_file"))   &
            CALL yaml_input_flag_vtk_flux_file(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "vtk_reg_file"))   &
            CALL yaml_input_flag_vtk_reg_file(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "vtk_src_file"))   &
            CALL yaml_input_flag_vtk_src_file(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "density_factor_file"))   &
            CALL yaml_input_flag_density_factor_file(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "print_xs"))   &
            CALL yaml_input_flag_print_xs(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "vtk_flux"))   &
            CALL yaml_input_flag_vtk_flux(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "vtk_reg"))   &
            CALL yaml_input_flag_vtk_reg(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "vtk_mat"))   &
            CALL yaml_input_flag_vtk_mat(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "vtk_src"))   &
            CALL yaml_input_flag_vtk_src(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "qdorder"))   &
            CALL yaml_input_flag_qdorder(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "qdtype"))   &
            CALL yaml_input_flag_qdtype(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "ngroups"))   &
            CALL yaml_input_flag_ngroups(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "pnorder"))   &
            CALL yaml_input_flag_pnorder(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "pnread"))   &
            CALL yaml_input_flag_pnread(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "upscattering"))   &
            CALL yaml_input_flag_upscattering(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "multiplying"))   &
            CALL yaml_input_flag_multiplying(data_string, set_defaults, sanity_check, verbose)
      !///////////////////////////////////////////////////////////////////////////
      IF (set_defaults .OR. sanity_check .OR. (key_string .EQ. "scatt_mult_included"))   &
            CALL yaml_input_flag_scatt_mult_included(data_string, set_defaults, sanity_check, verbose)



      set_defaults = .FALSE.
      sanity_check = .FALSE.
    END DO

  END SUBROUTINE yaml_read

  !===============================================================================
  ! [Execute] flag !==============================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_execute(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    !DATA_TYPE:: data_value
    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "# Allows for execution of yaml_input validation"
        WRITE(14, *) "execution: yes #Perform solve as described in yaml_input"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "# Allows for execution of yaml_input validation"
        WRITE(14, *) "execution: yes #Perform solve as described in yaml_input [DEFAULT]"
        WRITE(14, *) "execution: no  #Stop after parsing yaml_input"
        WRITE(14, *)
      END IF
      execution = 1
    ELSE IF(sanity_check) THEN
      IF (execution .NE. 1 .AND. execution .NE. 0) &
            CALL genError(0, '[execution] failed yaml_input validation')
    ELSE
      IF (data_string .EQ. "yes") THEN
        execution = 1
      ELSE IF (data_string .EQ. "no") THEN
        execution = 0
      ELSE
        CALL genError(0, "Invalid parameter for [execute] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_execute

  !===============================================================================
  ! [type] flag !=================================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_type(data_string, set_defaults, sanity_check, verbose)
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
            CALL genError(0, '[type] failed yaml_input validation')
    ELSE
      IF (data_string .EQ. "keig") THEN
        problem = 1
      ELSE IF (data_string .EQ. "fsrc") THEN
        problem = 0
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [type] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_type

  !===============================================================================
  ! [lambda] flag !===============================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_lambda(data_string, set_defaults, sanity_check, verbose)
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
            CALL genError(0, "Invalid parameter for [lambda] yaml_input flag")
    ELSE
      READ(data_string, *) data_int
      IF (data_int .GE. 0) THEN
        space_ord = data_int
      ELSE
        !Ammend error message
        CALL genError(0, "Invalid parameter for [lambda] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_lambda

  !===============================================================================
  ! [inflow] flag !===============================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_inflow(data_string, set_defaults, sanity_check, verbose)
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
        CALL genError(0, "Invalid parameter for [inflow] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_inflow

  !===============================================================================
  ! [piacc] flag !================================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_piacc(data_string, set_defaults, sanity_check, verbose)
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
        CALL genError(0, "Invalid parameter for [piacc] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_piacc

  !===============================================================================
  ! [sweep] flag !================================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_sweep(data_string, set_defaults, sanity_check, verbose)
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
        CALL genError(0, "Invalid parameter for [sweep] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_sweep

  !===============================================================================
  ! [page_sweep] flag !===========================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_page_sweep(data_string, set_defaults, sanity_check, verbose)
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
        CALL genError(0, "Invalid parameter for [page_sweep] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_page_sweep

  !===============================================================================
  ! [page_refl] flag !============================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_page_refl(data_string, set_defaults, sanity_check, verbose)
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
        CALL genError(0, "Invalid parameter for [page_refl] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_page_refl

  !===============================================================================
  ! [page_inflow] flag !==========================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_page_inflow(data_string, set_defaults, sanity_check, verbose)
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
        CALL genError(0, "Invalid parameter for [page_inflow] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_page_inflow

  !===============================================================================
  ! [maxouter] flag !=============================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_maxouter(data_string, set_defaults, sanity_check, verbose)
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
        CALL genError(0, "Invalid parameter for [maxouter] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_maxouter

  !===============================================================================
  ! [maxinner] flag !=============================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_maxinner(data_string, set_defaults, sanity_check, verbose)
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
        CALL genError(0, "Invalid parameter for [maxinner] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_maxinner

  !===============================================================================
  ! [innerconv] flag !============================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_innerconv(data_string, set_defaults, sanity_check, verbose)
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
        CALL genError(0, "Invalid parameter for [innerconv] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_innerconv

  !===============================================================================
  ! [outerconv] flag !============================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_outerconv(data_string, set_defaults, sanity_check, verbose)
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
        CALL genError(0, "Invalid parameter for [outerconv] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_outerconv

  !===============================================================================
  ! [kconv] flag !================================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_kconv(data_string, set_defaults, sanity_check, verbose)
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
        CALL genError(0, "Invalid parameter for [kconv] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_kconv

  !===============================================================================
  ! [keigsolve] flag !============================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_keigsolve(data_string, set_defaults, sanity_check, verbose)
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
        CALL genError(0, "Invalid parameter for [keigsolver] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_keigsolve

  !===============================================================================
  ! [jfnk_krsze] flag !===========================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_jfnk_krsze(data_string, set_defaults, sanity_check, verbose)
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
        CALL genError(0, "Invalid parameter for [jfnk_krsze] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_jfnk_krsze

  !===============================================================================
  ! [jfnk_maxkr] flag !===========================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_jfnk_maxkr(data_string, set_defaults, sanity_check, verbose)
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
        CALL genError(0, "Invalid parameter for [jfnk_maxkr] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_jfnk_maxkr

  !===============================================================================
  ! [jfnk_method] flag !==========================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_jfnk_method(data_string, set_defaults, sanity_check, verbose)
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
        CALL genError(0, "Invalid parameter for [jfnk_method] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_jfnk_method

  !===============================================================================
  ! [initial_guess] flag !========================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_initial_guess(data_string, set_defaults, sanity_check, verbose)
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
        CALL genError(0, "Invalid parameter for [initial_guess] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_initial_guess

  !===============================================================================
  ! [save_restart] flag !=========================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_save_restart(data_string, set_defaults, sanity_check, verbose)
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
        CALL genError(0, "Invalid parameter for [save_restart] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_save_restart

  !===============================================================================
  ! [ipiter] flag !===============================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_ipiter(data_string, set_defaults, sanity_check, verbose)
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
        CALL genError(0, "Invalid parameter for [ipiter] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_ipiter

  !===============================================================================
  ! [print_conv] flag !===========================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_print_conv(data_string, set_defaults, sanity_check, verbose)
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
        CALL genError(0, "Invalid parameter for [generic] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_print_conv

  !===============================================================================
  ! [density_factor] flag !=======================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_density_factor(data_string, set_defaults, sanity_check, verbose)
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
        WRITE(14, *) "density_factor: fromfile #Assign density factors based on yaml_input file"
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
        CALL genError(0, "Invalid parameter for [density_factor] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_density_factor

  !===============================================================================
  ! [source_file] flag !==========================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_source_file(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    LOGICAL :: data_logical

    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "# Specifies the name of the source yaml_input file. Accepts a string of length < 80"
        WRITE(14, *) "source_file: file.src"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "# Specifies the name of the source yaml_input file. Accepts a string of length < 80"
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
        CALL genError(0, "Invalid parameter for [source_file] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_source_file

  !===============================================================================
  ! [inflow_file] flag !==========================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_inflow_file(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    !DATA_TYPE:: data_value
    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "# Defines the filename for the boundary inflow yaml_input file"
        WRITE(14, *) "inflow_file: file.bc"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "# Defines the filename for the boundary inflow yaml_input file"
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
        CALL genError(0, "Invalid parameter for [inflow_file] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_inflow_file

  !===============================================================================
  ! [xs_file] flag !==============================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_xs_file(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    LOGICAL :: data_logical
    !DATA_TYPE:: data_value
    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "# Specifies the name of the cross-section yaml_input file. Accepts any string of length < 80"
        WRITE(14, *) "xs_file: file.xs"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "# Specifies the name of the cross-section yaml_input file. Accepts any string of length < 80"
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
        CALL genError(0, "Invalid parameter for [xs_file] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_xs_file

  !===============================================================================
  ! [mesh_file] flag !============================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_mesh_file(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    LOGICAL :: data_logical

    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "# Specifies the name of the mesh yaml_input file. Accepts any string of length < 80"
        WRITE(14, *) "mesh_file: file.thrm"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "# Specifies the name of the mesh yaml_input file. Accepts any string of length < 80"
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
        CALL genError(0, "Invalid parameter for [mesh_file] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_mesh_file

  !===============================================================================
  ! [flux_file] flag !============================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_flux_file(data_string, set_defaults, sanity_check, verbose)
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
        CALL genError(0, "Invalid parameter for [flux_file] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_flux_file

  !===============================================================================
  ! [vtk_flux_file] flag !========================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_vtk_flux_file(data_string, set_defaults, sanity_check, verbose)
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
        CALL genError(0, "Invalid parameter for [vtk_flux_file] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_vtk_flux_file

  !===============================================================================
  ! [quad_file] flag !============================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_quad_file(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    !DATA_TYPE:: data_value
    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "# Specifies the name of the quadrature yaml_input file. Accepts any string of length < 80"
        WRITE(14, *) "quad_file: file.quad"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "# Specifies the name of the quadrature yaml_input file. Accepts any string of length < 80"
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
        CALL genError(0, "Invalid parameter for [quad_file] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_quad_file

  !===============================================================================
  ! [restart_file] flag !=========================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_restart_file(data_string, set_defaults, sanity_check, verbose)
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
        CALL genError(0, "Invalid parameter for [restart_file] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_restart_file

  !===============================================================================
  ! [inguess_file] flag !=========================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_inguess_file(data_string, set_defaults, sanity_check, verbose)
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
        CALL genError(0, "Invalid parameter for [inguess_file] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_inguess_file

  !===============================================================================
  ! [vtk_mat_file] flag !=========================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_vtk_mat_file(data_string, set_defaults, sanity_check, verbose)
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
        CALL genError(0, "Invalid parameter for [vtk_mat_file] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_vtk_mat_file

  !===============================================================================
  ! [vtk_reg_file] flag !=========================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_vtk_reg_file(data_string, set_defaults, sanity_check, verbose)
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
        CALL genError(0, "Invalid parameter for [vtk_reg_file] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_vtk_reg_file

  !===============================================================================
  ! [vtk_src_file] flag !=========================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_vtk_src_file(data_string, set_defaults, sanity_check, verbose)
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
        CALL genError(0, "Invalid parameter for [vtk_src_file] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_vtk_src_file

  !===============================================================================
  ! [density_factor_file] flag !==================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_density_factor_file(data_string, set_defaults, sanity_check, verbose)
    IMPLICIT NONE

    CHARACTER(*):: data_string
    LOGICAL:: set_defaults, sanity_check
    INTEGER:: verbose
    !DATA_TYPE:: data_value
    IF(set_defaults) THEN
      IF (verbose .EQ. 1) THEN
        WRITE(14, *) "# Specifies the name of the density factor yaml_input file. Accepts any string of length < 80"
        WRITE(14, *) "density_factor_file: density_factors.dat"
        WRITE(14, *)
      ELSE IF (verbose .EQ. 2) THEN
        WRITE(14, *) "# Specifies the name of the density factor yaml_input file. Accepts any string of length < 80"
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
        CALL genError(0, "Invalid parameter for [generic] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_density_factor_file

  !===============================================================================
  ! [print_xs] flag !=============================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_print_xs(data_string, set_defaults, sanity_check, verbose)
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
        CALL genError(0, "Invalid parameter for [print_xs] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_print_xs

  !===============================================================================
  ! [vtk_flux] flag !=============================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_vtk_flux(data_string, set_defaults, sanity_check, verbose)
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
        CALL genError(0, "Invalid parameter for [vtk_flux] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_vtk_flux

  !===============================================================================
  ! [vtk_reg] flag !==============================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_vtk_reg(data_string, set_defaults, sanity_check, verbose)
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
        CALL genError(0, "Invalid parameter for [vtk_reg] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_vtk_reg

  !===============================================================================
  ! [vtk_mat] flag !==============================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_vtk_mat(data_string, set_defaults, sanity_check, verbose)
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
        CALL genError(0, "Invalid parameter for [vtk_mat] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_vtk_mat

  !===============================================================================
  ! [vtk_src] flag !==============================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_vtk_src(data_string, set_defaults, sanity_check, verbose)
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
        CALL genError(0, "Invalid parameter for [vtk_src] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_vtk_src

  !===============================================================================
  ! [qdorder] flag !==============================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_qdorder(data_string, set_defaults, sanity_check, verbose)
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
        CALL genError(0, "Invalid parameter for [qdorder] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_qdorder

  !===============================================================================
  ! [qdtype] flag !===============================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_qdtype(data_string, set_defaults, sanity_check, verbose)
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
        CALL genError(0, "Invalid parameter for [qdtype] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_qdtype

  !===============================================================================
  ! [ngroups] flag !==============================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_ngroups(data_string, set_defaults, sanity_check, verbose)
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
        CALL genError(0, "Invalid parameter for [ngroups] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_ngroups

  !===============================================================================
  ! [pnorder] flag !==============================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_pnorder(data_string, set_defaults, sanity_check, verbose)
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
        CALL genError(0, "Invalid parameter for [pnorder] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_pnorder

  !===============================================================================
  ! [pnread] flag !===============================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_pnread(data_string, set_defaults, sanity_check, verbose)
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
        CALL genError(0, "Invalid parameter for [pnread] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_pnread

  !===============================================================================
  ! [upscattering] flag !=========================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_upscattering(data_string, set_defaults, sanity_check, verbose)
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
        CALL genError(0, "Invalid parameter for [upscattering] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_upscattering

  !===============================================================================
  ! [multiplying] flag !==========================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_multiplying(data_string, set_defaults, sanity_check, verbose)
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
        CALL genError(0, "Invalid parameter for [multiplying] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_multiplying

  !===============================================================================
  ! [scatt_mult_included] flag !==================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_scatt_mult_included(data_string, set_defaults, sanity_check, verbose)
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
        CALL genError(0, "Invalid parameter for [scatt_mult_included] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_scatt_mult_included

  !===============================================================================
  ! [GENERIC] flag !==============================================================
  !===============================================================================
  SUBROUTINE yaml_input_flag_generic(data_string, set_defaults, sanity_check, verbose)
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
        CALL genError(0, "Invalid parameter for [generic] yaml_input flag")
      END IF
    END IF
  END SUBROUTINE yaml_input_flag_generic

  SUBROUTINE legacyv1_read(localunit,regmap)
    INTEGER,INTENT(IN) :: localunit
    CHARACTER(100000),INTENT(OUT) :: regmap
    !***********************************************************************
    !
    ! This subroutine reads the standard input file
    !
    !***********************************************************************

    ! local variables

    CHARACTER(100) :: buffer, fname, tchar
    LOGICAL :: done
    INTEGER :: i, rank,mpi_err,ios,legacyv
    CALL GET_COMMAND_ARGUMENT(1,fname)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)
    ! read title
    READ(localunit,*) jobname

    ! main read loop

    done = .FALSE.
    DO WHILE(done .EQV. .FALSE.)

      READ(localunit,101,END=999) buffer
      IF     ( INDEX( lowercase(buffer) ,'start') > 0 .AND. INDEX( lowercase(buffer) ,'problem')>0 ) THEN
        CALL legacyv1_read_problem
      ELSE IF( INDEX( lowercase(buffer) ,'start') > 0 .AND. INDEX( lowercase(buffer) ,'inout')>0   ) THEN
        CALL legacyv1_read_inout
      ELSE IF( INDEX( lowercase(buffer) ,'start') > 0 .AND. INDEX( lowercase(buffer) ,'cross_sections')>0   ) THEN
        CALL legacyv1_read_cross_sections
      ELSE IF( INDEX( lowercase(buffer) ,'start') > 0 .AND. INDEX( lowercase(buffer) ,'quadrature')>0   ) THEN
        CALL legacyv1_read_quadrature_field
      ELSE IF( INDEX( lowercase(buffer) ,'start') > 0 .AND. INDEX( lowercase(buffer) ,'postprocess')>0   ) THEN
        CALL legacyv1_read_postprocess_field
      ELSE IF( INDEX( lowercase(buffer) ,'start') > 0 .AND. INDEX( lowercase(buffer) ,'regionmap')>0   ) THEN
        CALL legacyv1_read_regionmap_field(regmap)
      ELSE IF( INDEX( lowercase(buffer) ,'end') > 0   .AND. INDEX( lowercase(buffer) ,'file')   >0 ) THEN
        done=.TRUE.
      END IF

    END DO

999 CONTINUE

101 FORMAT(A100)

  END SUBROUTINE legacyv1_read

  SUBROUTINE legacyv1_read_regionmap_field(regmap)

    ! reads the regionmap into array regmap

    ! Arguments

    CHARACTER(100000) :: regmap

    ! local variables

    CHARACTER(100) :: line, fname
    INTEGER :: l,lr
    INTEGER :: i, rank,mpi_err, localunit
    CALL GET_COMMAND_ARGUMENT(1,fname)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)
    localunit = rank+100

    ! read input line by line

    regmap=""
    DO WHILE(.TRUE.)
      READ(localunit,101) line
      IF ( INDEX( lowercase(line) ,'end') > 0 ) THEN
        RETURN
      ELSE
        l      = LEN(TRIM(regmap))
        lr     = LEN(TRIM(line))
        regmap(l+1:l+1+lr) = TRIM(line)
      END IF
    END DO

101 FORMAT(A100)

  END SUBROUTINE legacyv1_read_regionmap_field

  SUBROUTINE legacyv1_read_quadrature_field

    ! local variables
    INTEGER :: nwords,ntmp,i,nwwords,ios
    CHARACTER(100) :: buffer, fname
    CHARACTER(100) :: words(100),wwords(2)
    INTEGER :: rank,mpi_err, localunit
    CALL GET_COMMAND_ARGUMENT(1,fname)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)
    localunit = rank+100

    ! read loop over inout block
    DO WHILE(.TRUE.)
      READ(localunit,101) buffer
      IF ( INDEX( lowercase(buffer) ,'end') > 0 ) THEN
        RETURN
      ELSE
        CALL parse(buffer,";",words,nwords)
        DO i=1,nwords
          CALL parse(words(i),"=",wwords,nwwords)
          IF(nwwords .EQ. 2) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'qdtype'
            IF( TRIM(lowercase(wwords(1))) .EQ. 'qdtype' ) THEN
              wwords(2)=TRIM(lowercase(wwords(2)))
              IF      ( wwords(2) .EQ. 'levelsym') THEN
                quad_tpe=1
              ELSE IF ( wwords(2) .EQ. 'legcheb') THEN
                quad_tpe=2
              ELSE IF ( wwords(2) .EQ. 'fromfile') THEN
                quad_tpe=3
              ELSE
                WRITE(6,*) 'Error. This is not a valid quadrature type -- ',wwords(2),' --'
                WRITE(6,*) 'Execution will terminate.'
                STOP
              END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'qdorder'
            ELSE IF( TRIM(lowercase(wwords(1))) .EQ. 'qdorder' ) THEN
              wwords(2)=TRIM(lowercase(wwords(2)))
              READ(wwords(2),'(i10)',iostat=ios) quad_ord
              IF(ios.NE.0 ) THEN
                WRITE(6,*) 'Invalid quadrature order -- ',TRIM(wwords(2)),' --'
                STOP
              END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > default if keyword is unknown
            ELSE
              WRITE(6,*) 'Unknown keyword in quadrature specification -- ',TRIM(wwords(1)),' --'
              WRITE(6,*) 'Execution will terminate.'
              STOP
            END IF
          ELSE
            WRITE(6,*) 'Error while reading cross section specification'
            WRITE(6,*) 'Do not understand entry: ',TRIM(words(i))
            WRITE(6,*) 'Execution will terminate.'
            STOP
          END IF
        END DO
      END IF
    END DO

101 FORMAT(A100)

  END SUBROUTINE legacyv1_read_quadrature_field

  SUBROUTINE legacyv1_read_postprocess_field

    ! local variables
    INTEGER :: nwords, ntmp, i, l, j, nwwords, ios, nwwwords
    CHARACTER(1000) :: buffer, fname
    CHARACTER(1000) :: words(100), wwords(2), wwwords(100)
    INTEGER :: rank, mpi_err, localunit
    CALL GET_COMMAND_ARGUMENT(1,fname)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)
    localunit = rank+100

    ! read loop over inout block
    DO WHILE(.TRUE.)
      READ(localunit,101) buffer
      IF ( INDEX( lowercase(buffer) ,'end') > 0 ) THEN
        RETURN
      ELSE
        CALL parse(buffer,";",words,nwords)
        DO i=1,nwords
          CALL parse(words(i),"=",wwords,nwwords)
          IF(nwwords .EQ. 2) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'qdtype'
            IF( TRIM(lowercase(wwords(1))) .EQ. 'cartesian_map' ) THEN
              wwords(2)=TRIM(lowercase(wwords(2)))
              glob_do_cartesian_mesh = .TRUE.
              ! wwords must be an array with of length 9
              CALL parse(wwords(2), " ", wwwords, nwwwords)
              IF (nwwwords .NE. 9) THEN
                WRITE(6,*) 'Following cartesian map nine entries are required; Found: ',&
                      TRIM(wwords(2)),' has ', nwwwords, ' entries.'
              END IF
              glob_cmap_min_x = string_to_real(wwwords(1), 'Conversion to cartesian map xmin failed')
              glob_cmap_max_x = string_to_real(wwwords(2), 'Conversion to cartesian map xmax failed')
              IF (ABS(glob_cmap_max_x - glob_cmap_min_x) < small_real) THEN
                WRITE(6, *) "cartesian_map xmin and xmax are too close to each other"
              END IF
              glob_cmap_nx = string_to_int(wwwords(3), 'Conversion to cartesian map nx failed', 1)
              glob_cmap_min_y = string_to_real(wwwords(4), 'Conversion to cartesian map ymin failed')
              glob_cmap_max_y = string_to_real(wwwords(5), 'Conversion to cartesian map ymax failed')
              IF (ABS(glob_cmap_max_y - glob_cmap_min_y) < small_real) THEN
                WRITE(6, *) "cartesian_map xmin and xmax are too close to each other"
              END IF
              glob_cmap_ny = string_to_int(wwwords(6), 'Conversion to cartesian map ny failed', 1)
              glob_cmap_min_z = string_to_real(wwwords(7), 'Conversion to cartesian map zmin failed')
              glob_cmap_max_z = string_to_real(wwwords(8), 'Conversion to cartesian map zmax failed')
              IF (ABS(glob_cmap_max_z - glob_cmap_min_z) < small_real) THEN
                WRITE(6, *) "cartesian_map zmin and zmax are too close to each other"
              END IF
              glob_cmap_nz = string_to_int(wwwords(9), 'Conversion to cartesian map nz failed', 1)
            ELSE IF ( TRIM(lowercase(wwords(1))) .EQ. 'point_value_locations' ) THEN
              wwords(2)=TRIM(lowercase(wwords(2)))
              CALL parse(wwords(2), " ", wwwords, nwwwords)
              ! must be divisible by 3
              IF (modulo(nwwwords, 3) .ne. 0) THEN
                WRITE(6,*) 'point_value_locations number of entries must be divisible by 3; Found: ',&
                      TRIM(wwords(2)),' has ', nwwwords, ' entries.'
              ELSE
                number_point_flux_locations = nwwwords / 3
                ALLOCATE(point_flux_locations(number_point_flux_locations, 3))
                DO l = 1, number_point_flux_locations
                  DO j = 1, 3
                    point_flux_locations(l, j) = string_to_real(wwwords((l - 1) * 3 + j),&
                      'Conversion to point flux location failed')
                  END DO
                END DO
              END IF
            ELSE
              WRITE(6,*) 'Unknown keyword in postprocess specification -- ',TRIM(wwords(1)),' --'
              WRITE(6,*) 'Execution will terminate.'
              STOP
            END IF
          ELSE
            WRITE(6,*) 'Error while reading postprocess specification'
            WRITE(6,*) 'Do not understand entry: ',TRIM(words(i))
            WRITE(6,*) 'Execution will terminate.'
            STOP
          END IF
        END DO
      END IF
    END DO

101 FORMAT(A1000)

  END SUBROUTINE legacyv1_read_postprocess_field

  SUBROUTINE legacyv1_read_cross_sections

    ! local variables
    INTEGER :: nwords,ntmp,nwwords,ios
    CHARACTER(100) :: buffer, fname
    CHARACTER(100) :: words(100),wwords(2)
    INTEGER :: i, rank,mpi_err, localunit
    CALL GET_COMMAND_ARGUMENT(1,fname)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)
    localunit = rank+100

    ! read loop over inout block
    DO WHILE(.TRUE.)
      READ(localunit,101) buffer
      IF ( INDEX( lowercase(buffer) ,'end') > 0 ) THEN
        RETURN
      ELSE
        CALL parse(buffer,";",words,nwords)
        DO i=1,nwords
          CALL parse(words(i),"=",wwords,nwwords)
          IF(nwwords .EQ. 2) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'ngroups'
            IF( TRIM(lowercase(wwords(1))) .EQ. 'ngroups' ) THEN
              wwords(2)=TRIM(lowercase(wwords(2)))
              READ(wwords(2),'(i10)',iostat=ios) egmax
              IF(ios.NE.0 .OR. egmax<1) THEN
                WRITE(6,*) 'Invalid number of energy groups in cross section specification -- ',TRIM(wwords(2)),' --'
                STOP
              END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'pnorder'
            ELSE IF( TRIM(lowercase(wwords(1))) .EQ. 'pnorder' ) THEN
              wwords(2)=TRIM(lowercase(wwords(2)))
              READ(wwords(2),'(i10)',iostat=ios) scatt_ord
              IF(ios.NE.0 .OR. scatt_ord < 0) THEN
                WRITE(6,*) 'Invalid scattering expansion in cross section specification -- ',TRIM(wwords(2)),' --'
                STOP
              END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword pnread'
            ELSE IF( TRIM(lowercase(wwords(1))) .EQ. 'pnread' ) THEN
              wwords(2)=TRIM(lowercase(wwords(2)))
              READ(wwords(2),'(i10)',iostat=ios) xs_ord
              IF(ios.NE.0 .OR. xs_ord < 0) THEN
                WRITE(6,*) 'Invalid cross section expansion in cross section specification -- ',TRIM(wwords(2)),' --'
                STOP
              END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'upscattering'
            ELSE IF( TRIM(lowercase(wwords(1))) .EQ. 'upscattering' ) THEN
              wwords(2)=TRIM(lowercase(wwords(2)))
              IF      ( wwords(2) .EQ. 'yes') THEN
                upscattering=1
              ELSE IF ( wwords(2) .EQ. 'no') THEN
                upscattering=0
              ELSE
                WRITE(6,*) 'Error. This is not a valid upscattering flag -- ',wwords(2),' --'
                WRITE(6,*) 'Execution will terminate.'
                STOP
              END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword !'multiplying'
            ELSE IF( TRIM(lowercase(wwords(1))) .EQ. 'multiplying' ) THEN
              wwords(2)=TRIM(lowercase(wwords(2)))
              IF      ( wwords(2) .EQ. 'yes') THEN
                multiplying=1
              ELSE IF ( wwords(2) .EQ. 'no') THEN
                multiplying=0
              ELSE
                WRITE(6,*) 'Error. This is not a valid multiplying flag -- ',wwords(2),' --'
                WRITE(6,*) 'Execution will terminate.'
                STOP
              END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'scatt_mult_included'
            ELSE IF( TRIM(lowercase(wwords(1))) .EQ. 'scatt_mult_included' ) THEN
              wwords(2)=TRIM(lowercase(wwords(2)))
              IF      ( wwords(2) .EQ. 'yes') THEN
                scat_mult_flag=1
              ELSE IF ( wwords(2) .EQ. 'no') THEN
                scat_mult_flag=0
              ELSE
                WRITE(6,*) 'Error. This is not a valid scattering multiplier flag -- ',wwords(2),' --'
                WRITE(6,*) 'Execution will terminate.'
                STOP
              END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > default if keyword is unknown
            ELSE
              WRITE(6,*) 'Unknown keyword in cross section specification -- ',TRIM(wwords(1)),' --'
              WRITE(6,*) 'Execution will terminate.'
              STOP
            END IF
          ELSE
            WRITE(6,*) 'Error while reading cross section specification'
            WRITE(6,*) 'Do not understand entry: ',TRIM(words(i))
            WRITE(6,*) 'Execution will terminate.'
            STOP
          END IF
        END DO
      END IF
    END DO

101 FORMAT(A100)

  END SUBROUTINE legacyv1_read_cross_sections

  SUBROUTINE legacyv1_read_inout

    ! local variables
    INTEGER :: nwords,ntmp,i,nwwords,ios
    CHARACTER(100) :: buffer, fname
    CHARACTER(100) :: words(100),wwords(2)
    INTEGER :: rank,mpi_err, localunit
    CALL GET_COMMAND_ARGUMENT(1,fname)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)
    localunit = rank+100

    ! read loop over inout block
    DO WHILE(.TRUE.)
      READ(localunit,101) buffer
      IF ( INDEX( lowercase(buffer) ,'end') > 0 ) THEN
        RETURN
      ELSE
        CALL parse(buffer,";",words,nwords)
        DO i=1,nwords
          CALL parse(words(i),"=",wwords,nwwords)
          IF(nwwords .EQ. 2) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'meshi_file'
            IF( TRIM(lowercase(wwords(1))) .EQ. 'mesh_file' ) THEN
              mesh_filename=TRIM(wwords(2))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'inflow_file'
            ELSE IF( TRIM(lowercase(wwords(1))) .EQ. 'inflow_file' ) THEN
              finflow_filename=TRIM(wwords(2))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'source_file'
            ELSE IF( TRIM(lowercase(wwords(1))) .EQ. 'source_file' ) THEN
              source_filename=TRIM(wwords(2))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'flux_file'
            ELSE IF( TRIM(lowercase(wwords(1))) .EQ. 'flux_file' ) THEN
              flux_filename=TRIM(wwords(2))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'xs_file'
            ELSE IF( TRIM(lowercase(wwords(1))) .EQ. 'xs_file' ) THEN
              cross_section_filename=TRIM(wwords(2))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'quad_file'
            ELSE IF( TRIM(lowercase(wwords(1))) .EQ. 'density_factor_file' ) THEN
              dens_fact_filename=TRIM(wwords(2))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'quad_file'
            ELSE IF( TRIM(lowercase(wwords(1))) .EQ. 'quad_file' ) THEN
              quad_file=TRIM(wwords(2))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'vtk'
            ELSE IF( TRIM(lowercase(wwords(1))) .EQ. 'vtk' ) THEN
              wwords(2)=TRIM(lowercase(wwords(2)))
              IF( INDEX(wwords(2),'flux') > 0 ) vtk_flux_output=1
              IF( INDEX(wwords(2),'mat'  ) > 0 ) vtk_mat_output=1
              IF( INDEX(wwords(2),'reg') > 0 ) vtk_reg_output=1
              IF( INDEX(wwords(2),'src') > 0 ) vtk_src_output=1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'vtk_flux_file'
            ELSE IF( TRIM(lowercase(wwords(1))) .EQ. 'vtk_flux_file' ) THEN
              vtk_flux_filename=TRIM(wwords(2))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'vtk_mat_file'
            ELSE IF( TRIM(lowercase(wwords(1))) .EQ. 'vtk_mat_file' ) THEN
              vtk_mat_filename=TRIM(wwords(2))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'vtk_reg_file'
            ELSE IF( TRIM(lowercase(wwords(1))) .EQ. 'vtk_reg_file' ) THEN
              vtk_reg_filename=TRIM(wwords(2))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'vtk_src_file'
            ELSE IF( TRIM(lowercase(wwords(1))) .EQ. 'vtk_src_file' ) THEN
              vtk_src_filename=TRIM(wwords(2))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'restart_file'
            ELSE IF( TRIM(lowercase(wwords(1))) .EQ. 'restart_file' ) THEN
              dump_file=TRIM(wwords(2))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'inguess_file'
            ELSE IF( TRIM(lowercase(wwords(1))) .EQ. 'inguess_file' ) THEN
              inguess_file=TRIM(wwords(2))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'inguess_file'
            ELSE IF( TRIM(lowercase(wwords(1))) .EQ. 'cartesian_map_file' ) THEN
              cartesian_map_filename = TRIM(wwords(2))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'print_xs'
            ELSE IF( TRIM(lowercase(wwords(1))) .EQ. 'print_xs' ) THEN
              wwords(2)=TRIM(lowercase(wwords(2)))
              IF      ( wwords(2) .EQ. 'yes') THEN
                print_xs_flag=1
              ELSE IF ( wwords(2) .EQ. 'no') THEN
                print_xs_flag=0
              ELSE
                WRITE(6,*) 'Error. This is not a valid cross section print option (yes/no) -- ',wwords(2),' --'
                WRITE(6,*) 'Execution will terminate.'
                STOP
              END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > default if keyword is unknown
            ELSE
              WRITE(6,*) 'Unknown keyword in inout field -- ',TRIM(wwords(1)),' --'
              WRITE(6,*) 'Execution will terminate.'
              STOP
            END IF
          ELSE
            WRITE(6,*) 'Error while reading inout specification'
            WRITE(6,*) 'Do not understand entry: ',TRIM(words(i))
            WRITE(6,*) 'Execution will terminate.'
            STOP
          END IF
        END DO
      END IF
    END DO

101 FORMAT(A100)

  END SUBROUTINE legacyv1_read_inout

  SUBROUTINE legacyv1_read_problem

    ! local variables
    INTEGER :: nwords,ntmp,i,nwwords,ios
    CHARACTER(100) :: buffer, fname
    CHARACTER(100) :: words(100),wwords(2)
    INTEGER :: rank,mpi_err, localunit
    CALL GET_COMMAND_ARGUMENT(1,fname)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)
    localunit = rank+100

    ! read loop over problem block
    DO WHILE(.TRUE.)
      READ(localunit,101) buffer
      IF ( INDEX( lowercase(buffer) ,'end') > 0 ) THEN
        RETURN
      ELSE
        CALL parse(buffer,";",words,nwords)
        DO i=1,nwords
          CALL parse(words(i),"=",wwords,nwwords)
          IF(nwwords .EQ. 2) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'type'
            IF( TRIM(lowercase(wwords(1))) .EQ. 'type' ) THEN
              wwords(2)=TRIM(lowercase(wwords(2)))
              IF      ( wwords(2) .EQ. 'keig') THEN
                problem=1
              ELSE IF ( wwords(2) .EQ. 'fsrc') THEN
                problem=0
              ELSE
                WRITE(6,*) 'Error. This is not a valid problem type -- ',wwords(2),' --'
                WRITE(6,*) 'Execution will terminate.'
                STOP
              END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'keigsolver'
            ELSE IF( TRIM(lowercase(wwords(1))) .EQ. 'keigsolver' ) THEN
              wwords(2)=TRIM(lowercase(wwords(2)))
              IF      ( wwords(2) .EQ. 'pi') THEN
                eig_switch=0
              ELSE IF ( wwords(2) .EQ. 'jfnk') THEN
                eig_switch=1
              ELSE
                WRITE(6,*) 'Error. This is not a valid eigenvalue solver -- ',wwords(2),' --'
                WRITE(6,*) 'Execution will terminate.'
                STOP
              END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'lambda'
            ELSE IF( TRIM(lowercase(wwords(1))) .EQ. 'lambda' ) THEN
              wwords(2)=TRIM(lowercase(wwords(2)))
              READ(wwords(2),'(i10)',iostat=ios) space_ord
              IF(ios.NE.0) THEN
                WRITE(6,*) 'Invalid spatial order in problem specification -- ',TRIM(wwords(2)),' --'
                STOP
              END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'inflow'
            ELSE IF( TRIM(lowercase(wwords(1))) .EQ. 'inflow' ) THEN
              wwords(2)=TRIM(lowercase(wwords(2)))
              IF      ( wwords(2) .EQ. 'yes') THEN
                finflow=1
              ELSE IF ( wwords(2) .EQ. 'no') THEN
                finflow=0
              ELSE
                WRITE(6,*) 'Error. This is not a valid inflow flag -- ',wwords(2),' --'
                WRITE(6,*) 'Execution will terminate.'
                STOP
              END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'PIacc'
            ELSE IF( TRIM(lowercase(wwords(1))) .EQ. 'piacc' ) THEN
              wwords(2)=TRIM(lowercase(wwords(2)))
              IF      ( wwords(2) .EQ. 'errmode') THEN
                outer_acc=2
              ELSE IF ( wwords(2) .EQ. 'none') THEN
                outer_acc=1
              ELSE
                WRITE(6,*) 'Error. This is not a valid acceleration option for PI -- ',wwords(2),' --'
                WRITE(6,*) 'Execution will terminate.'
                STOP
              END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'sweep'
            ELSE IF( TRIM(lowercase(wwords(1))) .EQ. 'sweep' ) THEN
              wwords(2)=TRIM(lowercase(wwords(2)))
              IF      ( wwords(2) .EQ. 'precomp') THEN
                sweep_tpe=1
              ELSE
                WRITE(6,*) 'Error. This is not a valid mesh sweep type -- ',wwords(2),' --'
                WRITE(6,*) 'Execution will terminate.'
                STOP
              END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'sweep_page'
            ELSE IF( TRIM(lowercase(wwords(1))) .EQ. 'page_sweep' ) THEN
              wwords(2)=TRIM(lowercase(wwords(2)))
              IF      ( wwords(2) .EQ. 'yes') THEN
                page_sweep=1
              ELSE IF ( wwords(2) .EQ. 'no') THEN
                page_sweep=0
              ELSE
                WRITE(6,*) 'Error. This is not a valid sweep page option -- ',wwords(2),' --'
                WRITE(6,*) 'Execution will terminate.'
                STOP
              END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'page_refl'
            ELSE IF( TRIM(lowercase(wwords(1))) .EQ. 'page_refl' ) THEN
              wwords(2)=TRIM(lowercase(wwords(2)))
              IF      ( wwords(2) .EQ. 'page') THEN
                page_refl=1
              ELSE IF ( wwords(2) .EQ. 'save') THEN
                page_refl=0
              ELSE IF ( wwords(2) .EQ. 'inner') THEN
                page_refl=2
              ELSE
                WRITE(6,*) 'Error. This is not a valid page reflective BC option -- ',wwords(2),' --'
                WRITE(6,*) 'Execution will terminate.'
                STOP
              END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'page_iflw'
            ELSE IF( TRIM(lowercase(wwords(1))) .EQ. 'page_inflow' ) THEN
              wwords(2)=TRIM(lowercase(wwords(2)))
              IF      ( wwords(2) .EQ. 'bygroup') THEN
                page_iflw=1
              ELSE IF ( wwords(2) .EQ. 'all') THEN
                page_iflw=0
              ELSE
                WRITE(6,*) 'Error. This is not a valid page inflow option -- ',wwords(2),' --'
                WRITE(6,*) 'Execution will terminate.'
                STOP
              END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'kconv'
            ELSE IF( TRIM(lowercase(wwords(1))) .EQ. 'kconv' ) THEN
              wwords(2)=TRIM(lowercase(wwords(2)))
              READ(wwords(2),*,iostat=ios) k_conv
              IF(ios.NE.0) THEN
                WRITE(6,*) 'Invalid stopping criterion for keff in problem specification -- ',TRIM(wwords(2)),' --'
                STOP
              END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'innerconv'
            ELSE IF( TRIM(lowercase(wwords(1))) .EQ. 'innerconv' ) THEN
              wwords(2)=TRIM(lowercase(wwords(2)))
              READ(wwords(2),*,iostat=ios) inner_conv
              IF(ios.NE.0) THEN
                WRITE(6,*) 'Invalid stopping criterion for inner iterations in problem specification -- ',TRIM(wwords(2)),' --'
                STOP
              END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'outerconv'
            ELSE IF( TRIM(lowercase(wwords(1))) .EQ. 'outerconv' ) THEN
              wwords(2)=TRIM(lowercase(wwords(2)))
              READ(wwords(2),*,iostat=ios) outer_conv
              IF(ios.NE.0) THEN
                WRITE(6,*) 'Invalid stopping criterion for outer iterations in problem specification -- ',TRIM(wwords(2)),' --'
                STOP
              END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'maxinner'
            ELSE IF( TRIM(lowercase(wwords(1))) .EQ. 'maxinner' ) THEN
              wwords(2)=TRIM(lowercase(wwords(2)))
              READ(wwords(2),'(i10)',iostat=ios) max_inner
              IF(ios.NE.0 .OR. max_inner<1 ) THEN
                WRITE(6,*) 'Invalid maximum number of inner iteration in problem specification -- ',TRIM(wwords(2)),' --'
                STOP
              END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'maxouter'
            ELSE IF( TRIM(lowercase(wwords(1))) .EQ. 'maxouter' ) THEN
              wwords(2)=TRIM(lowercase(wwords(2)))
              READ(wwords(2),'(i10)',iostat=ios) max_outer
              IF(ios.NE.0 .OR. max_outer<1 ) THEN
                WRITE(6,*) 'Invalid maximum number of outer iteration in problem specification -- ',TRIM(wwords(2)),' --'
                STOP
              END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'jfnk_krsze'
            ELSE IF( TRIM(lowercase(wwords(1))) .EQ. 'jfnk_krsze' ) THEN
              wwords(2)=TRIM(lowercase(wwords(2)))
              READ(wwords(2),'(i10)',iostat=ios) rd_restart
              IF(ios.NE.0 .OR. rd_restart<1 ) THEN
                WRITE(6,*) 'Invalid Krylov subspace size for JFNK in problem specification -- ',TRIM(wwords(2)),' --'
                STOP
              END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'jfnk_maxkr'
            ELSE IF( TRIM(lowercase(wwords(1))) .EQ. 'jfnk_maxkr' ) THEN
              wwords(2)=TRIM(lowercase(wwords(2)))
              READ(wwords(2),'(i10)',iostat=ios) rd_max_kit
              IF(ios.NE.0 .OR. rd_max_kit < 1 ) THEN
                WRITE(6,*) 'Invalid maximum number of Krylov iterations for JFNK in problem specification -- '&
                      ,TRIM(wwords(2)),' --'
                STOP
              END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'jfnk_method'
            ELSE IF( TRIM(lowercase(wwords(1))) .EQ. 'jfnk_method' ) THEN
              wwords(2)=TRIM(lowercase(wwords(2)))
              IF      ( wwords(2) .EQ. 'outer') THEN
                rd_method=1
              ELSE IF ( wwords(2) .EQ. 'flat') THEN
                rd_method=2
              ELSE IF ( wwords(2) .EQ. 'flat_wds') THEN
                rd_method=3
              ELSE
                WRITE(6,*) 'Error. This is not a valid jfnk solution method -- ',wwords(2),' --'
                WRITE(6,*) 'Execution will terminate.'
                STOP
              END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'initial_guess'
            ELSE IF( TRIM(lowercase(wwords(1))) .EQ. 'initial_guess' ) THEN
              wwords(2)=TRIM(lowercase(wwords(2)))
              IF      ( wwords(2) .EQ. 'yes') THEN
                inguess_flag=1
              ELSE IF ( wwords(2) .EQ. 'no') THEN
                inguess_flag=0
              ELSE
                WRITE(6,*) 'Error. This is not a valid initial guess option (yes/no) -- ',wwords(2),' --'
                WRITE(6,*) 'Execution will terminate.'
                STOP
              END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'save_restart'
            ELSE IF( TRIM(lowercase(wwords(1))) .EQ. 'save_restart' ) THEN
              wwords(2)=TRIM(lowercase(wwords(2)))
              IF      ( wwords(2) .EQ. 'yes') THEN
                dump_flag=1
              ELSE IF ( wwords(2) .EQ. 'no') THEN
                dump_flag=0
              ELSE
                WRITE(6,*) 'Error. This is not a valid restart file option (yes/no) -- ',wwords(2),' --'
                WRITE(6,*) 'Execution will terminate.'
                STOP
              END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword !'ipiter'
            ELSE IF( TRIM(lowercase(wwords(1))) .EQ. 'ipiter' ) THEN
              wwords(2)=TRIM(lowercase(wwords(2)))
              READ(wwords(2),'(i10)',iostat=ios) ipow
              IF(ios.NE.0 ) THEN
                WRITE(6,*) 'Invalid number of initial power iterations -- '&
                      ,TRIM(wwords(2)),' --'
                STOP
              END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'print_conv'
            ELSE IF( TRIM(lowercase(wwords(1))) .EQ. 'print_conv' ) THEN
              wwords(2)=TRIM(lowercase(wwords(2)))
              IF      ( wwords(2) .EQ. 'yes') THEN
                print_conv=1
              ELSE IF ( wwords(2) .EQ. 'no') THEN
                print_conv=0
              ELSE
                WRITE(6,*) 'Error. This is not a valid execution option (yes/no) -- ',wwords(2),' --'
                WRITE(6,*) 'Execution will terminate.'
                STOP
              END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'density factor option'
            ELSE IF( TRIM(lowercase(wwords(1))) .EQ. 'density_factor' ) THEN
              wwords(2)=TRIM(lowercase(wwords(2)))
              IF      ( wwords(2) .EQ. 'no') THEN
                dfact_opt = 0
              ELSE
                IF      ( wwords(2) .EQ. 'byvolume') THEN
                  dfact_opt = 1
                ELSE IF ( wwords(2) .EQ. 'fromfile') THEN
                  dfact_opt = 2
                ELSE
                  WRITE(6,*) 'Error. This is not a valid density factor option &
                        (no/byvolume/fromfile) -- ',wwords(2),' --'
                  WRITE(6,*) 'Execution will terminate.'
                  STOP
                END IF
              END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'execution'
            ELSE IF( TRIM(lowercase(wwords(1))) .EQ. 'execution' ) THEN
              wwords(2)=TRIM(lowercase(wwords(2)))
              IF      ( wwords(2) .EQ. 'yes') THEN
                execution=1
              ELSE IF ( wwords(2) .EQ. 'no') THEN
                execution=0
              ELSE
                WRITE(6,*) 'Error. This is not a valid execution option (yes/no) -- ',wwords(2),' --'
                WRITE(6,*) 'Execution will terminate.'
                STOP
              END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > default if keyword is unknown
            ELSE
              WRITE(6,*) 'Unknown keyword in problem specification -- ',TRIM(wwords(1)),' --'
              WRITE(6,*) 'Execution will terminate.'
              STOP
            END IF
          ELSE
            WRITE(6,*) 'Error while reading problem specification'
            WRITE(6,*) 'Do not understand entry: ',TRIM(words(i))
            WRITE(6,*) 'Execution will terminate.'
            STOP
          END IF
        END DO
      END IF

    END DO

101 FORMAT(A100)
  END SUBROUTINE legacyv1_read_problem

  REAL(kind=d_t) FUNCTION string_to_real(string, msg)
    CHARACTER(100), INTENT(in) :: string
    CHARACTER(100), INTENT(in) :: msg

    INTEGER(kind=li) :: ios

    READ(string, *, iostat=ios) string_to_real
    IF(ios .NE. 0) THEN
      WRITE(6,*) TRIM(msg), TRIM(string)
      STOP
    END IF
  END FUNCTION string_to_real

  INTEGER(kind=li) FUNCTION string_to_int(string, msg, min_int)
    CHARACTER(100), INTENT(in) :: string
    CHARACTER(100), INTENT(in) :: msg
    INTEGER(kind=li), OPTIONAL, INTENT(inout) :: min_int

    INTEGER(kind=li) :: ios

    IF (.NOT. PRESENT(min_int)) min_int = -glob_max_int
    READ(string, '(i10)', iostat=ios) string_to_int
    IF(ios .NE. 0 .OR. string_to_int < min_int) THEN
      WRITE(6,*) TRIM(msg), TRIM(string)
      STOP
    END IF
  END FUNCTION string_to_int

END MODULE read_module_legacy