MODULE read_module
  !***********************************************************************
  !
  ! The module io contains subroutines and variables pertinent to input and
  ! output.
  !
  !***********************************************************************

  ! User derived-type modules

  USE mpi
  USE stringmod
  USE types
  USE parameter_types
  USE filename_types
  USE vector_types
  USE cross_section_types
  USE geometry_types
  USE angle_types
  USE multindex_types
  USE globals

  ! Use modules that pertain setting up problem

  USE read_cross_section_module
  USE read_mesh_module
  USE read_source_module
  USE read_inflow_module
  USE quadrature_module
  USE check_input
  USE error_module
  USE read_inp_module_legacy
  USE read_inp_module

  IMPLICIT NONE

CONTAINS

  SUBROUTINE READ
    !*********************************************************************
    !
    ! Subroutine read calls subroutines to read input and mesh file
    !
    !*********************************************************************

    ! Declare temporary variable

    INTEGER(kind=li) :: alloc_stat
    INTEGER ::rank,mpi_err
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)

    ! Set the defaults

    CALL  set_default

    ! Call read_input to read input file

    CALL read_input

    ! Call read_mesh to read tetrahedral mesh file

    CALL read_tetmesh

    ! Check the input

    CALL check_standard_input

    ! Echo input
    IF (rank.EQ.0) THEN
      CALL echo_input
    END IF

    ! Call read_xs to read cross-section file

    CALL read_xs

    ! Print cross sections if desired
    IF(rank .EQ. 0) THEN
      IF(print_xs_flag .EQ. 1 ) CALL print_xs
      IF(vtk_mat_output .EQ. 1) CALL plot_material
      IF(vtk_reg_output .EQ. 1) CALL plot_region
    END IF

    ! If execution is not desired then stop here

    IF(execution .EQ. 0) CALL raise_fatal_error('*** Not enough memory ***')

    !set the namom here so we know for the source file
    namom=(scatt_ord+1)**2

    ! Call generate_num_moments and generate_multi_index to create indices
    CALL generate_num_moments(space_ord,num_moments_v,num_moments_f)

    ALLOCATE(index_v(num_moments_v),index_f(num_moments_f),&
          stat=alloc_stat)
    IF(alloc_stat /= 0) CALL raise_fatal_error("*** Not enough memory ***")

    CALL generate_multi_index(space_ord,num_moments_v,num_moments_f,&
          index_v,index_f)

    ! Call read_src to read external source file

    IF(problem == 0)THEN
      CALL read_src
    END IF

    ! Call quad_gen to generate quadrature

    CALL quad_gen

    ! Call read_finflow to read fixed inflow bc

    IF(finflow /= 0 .AND. problem .EQ. 0)THEN
      CALL read_finflow
    END IF

  END SUBROUTINE READ

  SUBROUTINE generate_num_moments(spord,numv,numf)
    !*********************************************************************
    !
    ! Subroutine generates the number of total moments
    !
    !*********************************************************************

    ! Pass input parameters

    INTEGER(kind=li), INTENT(in) :: spord
    INTEGER(kind=li), INTENT(out) :: numv, numf

    ! Define temporary variables
    INTEGER(kind=li) :: i, j, k, l

    l=1
    IF(spord <= 0)THEN
      DO k=0, ABS(spord)
        DO j=0, ABS(spord)
          DO i=0, ABS(spord)
            IF(i+j+k <= ABS(spord))THEN
              l=l+1
            END IF
          END DO
        END DO
      END DO
      numv=l-1
      l=1
      DO j=0, ABS(spord)
        DO i=0, ABS(spord)
          IF(i+j <= ABS(spord))THEN
            l=l+1
          END IF
        END DO
      END DO
      numf=l-1
    ELSE
      DO k=0, spord
        DO j=0, spord
          DO i=0, spord
            IF(i <= spord .OR. &
                  j <= spord .OR. &
                  k <= spord)THEN
              l=l+1
            END IF
          END DO
        END DO
      END DO
      numv=l-1
      l=1
      DO j=0, 2*spord
        DO i=0, 2*spord
          IF(i <= 2*spord .OR. &
                j <= 2*spord)THEN
            l=l+1
          END IF
        END DO
      END DO
      numf=l-1
    END IF

  END SUBROUTINE generate_num_moments

  SUBROUTINE generate_multi_index(spord,numv,numf,indv,indf)
    !*********************************************************************
    !
    ! Subroutine generates the multi-index
    !
    !*********************************************************************

    ! Pass input parameters

    INTEGER(kind=li), INTENT(in) :: spord, numv,&
          numf

    ! Pass index derived type
    TYPE(indices_v), DIMENSION(numv), INTENT(inout) :: indv
    TYPE(indices_f), DIMENSION(numf), INTENT(inout) :: indf

    ! Define temporary variables
    INTEGER(kind=li) :: i, j, k, l

    l=1

    IF(spord <= 0)THEN
      DO i=0, ABS(spord)
        indv(l)%i1=i
        indv(l)%i2=0
        indv(l)%i3=0
        l=l+1
      END DO
      DO j=1, ABS(spord)
        indv(l)%i1=0
        indv(l)%i2=j
        indv(l)%i3=0
        l=l+1
      END DO
      DO k=1, ABS(spord)
        indv(l)%i1=0
        indv(l)%i2=0
        indv(l)%i3=k
        l=l+1
      END DO
      DO j=1, ABS(spord)
        DO i=1, ABS(spord)
          IF(i+j <= ABS(spord))THEN
            indv(l)%i1=i
            indv(l)%i2=j
            indv(l)%i3=0
            l=l+1
          END IF
        END DO
      END DO
      DO k=1, ABS(spord)
        DO i=1, ABS(spord)
          IF(i+k <= ABS(spord))THEN
            indv(l)%i1=i
            indv(l)%i2=0
            indv(l)%i3=k
            l=l+1
          END IF
        END DO
      END DO
      DO k=1, ABS(spord)
        DO j=1, ABS(spord)
          IF(j+k <= ABS(spord))THEN
            indv(l)%i1=0
            indv(l)%i2=j
            indv(l)%i3=k
            l=l+1
          END IF
        END DO
      END DO
      DO k=1, ABS(spord)
        DO j=1, ABS(spord)
          DO i=1, ABS(spord)
            IF(i+j+k <= ABS(spord))THEN
              indv(l)%i1=i
              indv(l)%i2=j
              indv(l)%i3=k
              l=l+1
            END IF
          END DO
        END DO
      END DO
      IF(l == numv+1)THEN
      ELSE
        CALL raise_fatal_error("Counter fail to correctly count number of moments!")
      END IF
      l=1
      DO i=0, ABS(spord)
        indf(l)%i1=i
        indf(l)%i2=0
        l=l+1
      END DO
      DO j=1, ABS(spord)
        indf(l)%i1=0
        indf(l)%i2=j
        l=l+1
      END DO
      DO j=1, ABS(spord)
        DO i=1, ABS(spord)
          IF(i+j <= ABS(spord))THEN
            indf(l)%i1=i
            indf(l)%i2=j
            l=l+1
          END IF
        END DO
      END DO
      IF(l == numf+1)THEN
      ELSE
        CALL raise_fatal_error("Counter fail to correctly count number of moments!")
      END IF
    ELSE
      DO i=0, spord
        indv(l)%i1=i
        indv(l)%i2=0
        indv(l)%i3=0
        l=l+1
      END DO
      DO j=1, spord
        indv(l)%i1=0
        indv(l)%i2=j
        indv(l)%i3=0
        l=l+1
      END DO
      DO k=1, spord
        indv(l)%i1=0
        indv(l)%i2=0
        indv(l)%i3=k
        l=l+1
      END DO
      DO j=1, spord
        DO i=1, spord
          IF(i <= spord .OR. &
                j <= spord)THEN
            indv(l)%i1=i
            indv(l)%i2=j
            indv(l)%i3=0
            l=l+1
          END IF
        END DO
      END DO
      DO k=1, spord
        DO i=1, spord
          IF(i <= spord .OR. &
                k <= spord)THEN
            indv(l)%i1=i
            indv(l)%i2=0
            indv(l)%i3=k
            l=l+1
          END IF
        END DO
      END DO
      DO k=1, spord
        DO j=1, spord
          IF(j <= spord .OR. &
                k <= spord)THEN
            indv(l)%i1=0
            indv(l)%i2=j
            indv(l)%i3=k
            l=l+1
          END IF
        END DO
      END DO
      DO k=1, spord
        DO j=1, spord
          DO i=1, spord
            IF(i <= spord .OR. &
                  j <= spord .OR. &
                  k <= spord)THEN
              indv(l)%i1=i
              indv(l)%i2=j
              indv(l)%i3=k
              l=l+1
            END IF
          END DO
        END DO
      END DO
      IF(l == numv+1)THEN
      ELSE
        CALL raise_fatal_error("Counter fail to correctly count number of moments!")
      END IF
      l=1
      DO i=0, 2*spord
        indf(l)%i1=i
        indf(l)%i2=0
        l=l+1
      END DO
      DO j=1, 2*spord
        indf(l)%i1=0
        indf(l)%i2=j
        l=l+1
      END DO
      DO j=1, 2*spord
        DO i=1, 2*spord
          IF(i <= 2*spord .OR. &
                j <= 2*spord)THEN
            indf(l)%i1=i
            indf(l)%i2=j
            l=l+1
          END IF
        END DO
      END DO
      IF(l == numf+1)THEN
      ELSE
        CALL raise_fatal_error("Counter fail to correctly count number of moments!")
      END IF
    END IF

  END SUBROUTINE generate_multi_index

  SUBROUTINE set_default

    ! Set the input default

    problem                     = 1
    eig_switch                  = 0
    space_ord                   = 0
    finflow                     = 0
    finflow_filename            = "finflow.dat"
    outer_acc                   = 1
    page_sweep                  = 0
    page_refl                   = 0
    page_iflw                   = 0
    k_conv                      = 1.0E-4_d_t
    inner_conv                  = 1.0E-4_d_t
    outer_conv                  = 1.0E-3_d_t
    max_inner                   = 10
    max_outer                   = 100
    rd_restart                  = 25
    rd_max_kit                  = 250
    rd_method                   = 2
    inguess_flag                = 0
    inguess_file                = "initial_guess.dat"
    dump_flag                   = 0
    dump_file                   = "restart.out"
    ipow                        = 0
    print_conv                  = 0
    dfact_opt                   = 0
    dens_fact_filename          = "density_factor.dat"
    execution                   = 1
    mesh_filename               = "mesh.thrm"
    source_filename             = "source.dat"
    flux_filename               = "flux.out"
    cross_section_filename      = "xs.dat"
    vtk_flux_output             = 0
    vtk_flux_filename           = "flux.vtk"
    vtk_mat_output              = 0
    vtk_mat_filename            = "mat.vtk"
    vtk_reg_output              = 0
    vtk_reg_filename            = "reg.vtk"
    vtk_src_output              = 0
    vtk_src_filename            = "src.vtk"
    glob_do_cartesian_mesh      = .FALSE.
    cartesian_map_filename      = "cartesian_map.out"
    print_xs_flag               = 0
    egmax                       = 1
    scatt_ord                   = 0
    xs_ord                      = 0
    upscattering                = 1
    multiplying                 = 1
    scat_mult_flag              = 0
    quad_tpe                    = 1
    quad_file                   = "quad.dat"
    quad_ord                    = 4
    number_point_flux_locations = 0_li
    !only 1 sweep type right now and so not an input
    sweep_tpe                   = 1

  END SUBROUTINE set_default

  SUBROUTINE check_standard_input

    ! Check the input for errors

    ! problem must be 0 or 1
    IF(problem > 1 .OR. problem < 0)THEN
      CALL raise_fatal_error("Please select problem type 1 (eigenvalue search) or &
        & 2 (external source) to perform calculation")
    ENDIF

    ! jfnk parameters
    IF(eig_switch==1) THEN
      IF (rd_restart<0) THEN
        CALL raise_fatal_error("Number of Krylov Iterations between restarts must be greater than 0")
      END IF
      IF (rd_max_kit<rd_restart) THEN
        CALL raise_fatal_error("Maximum number of krylov iterations must be greater than &
          & number of iterations between restarts.")
      END IF
      IF (rd_method<1 .OR. rd_method >3) THEN
        CALL raise_fatal_error("Method has to be 1 (Outer iteration with lagged upscattering),  &
          & 2(Flat iteration), or 3(Flat iteration with updated downscattering).")
      END IF
    END IF

    ! in case problem == 1, no inflow file allowed

    IF(problem .EQ. 1 .AND. finflow .EQ. 1) THEN
      CALL raise_fatal_error("Inflow file only allowed for fixed source problems")
    END IF

    ! in case problem == 1, no vtk source file allowed

    IF(problem .EQ. 1 .AND. vtk_src_output .EQ. 1) THEN
      CALL raise_fatal_error("VTK source file only allowed for fixed source problems")
    END IF


  END SUBROUTINE check_standard_input

  SUBROUTINE echo_input

    ! Echo out problem specifications

    IF (rank .EQ. 0) THEN
      CALL printlog('')
      CALL printlog("--------------------------------------------------------")
      CALL printlog("   Input Summary  ")
      CALL printlog("--------------------------------------------------------")

      IF(problem == 1)THEN
        IF ( eig_switch == 0 ) THEN
          CALL printlog("Eigenvalue calculation using PI selected")
        ELSE
          CALL printlog("Eigenvalue calculation using JFNK selected")
          IF(rd_method == 1) THEN
            CALL printlog("Method: F(u) is evaluated using one outer iteration with lagged upscattering.")
          ELSE IF (rd_method == 2) THEN
            CALL printlog("Method: F(u) is evaluated using flat iteration scheme.")
          ELSE
            CALL printlog("Method: F(u) is evaluated using flat iteration scheme wit updated downscattering.")
          END IF
        END IF
      ELSE IF(problem == 0)THEN
        CALL printlog("External source calculation selected")
      END IF

      ! Print rest of input

      WRITE(amsg,'(2A)')   "Problem Title:                                              ", jobname
      CALL printlog(amsg)
      WRITE(amsg,101) "Spatial order:                                              ", space_ord
      CALL printlog(amsg)
      IF      (sweep_tpe .EQ. 1) THEN
        CALL printlog('Precomputed mesh sweep is used.')
      END IF
      IF      (outer_acc .EQ. 1 .AND. problem .EQ. 1 .AND. eig_switch .EQ. 0) THEN
        CALL printlog('Power iterations are not accelerated')
      ELSE IF (outer_acc .EQ. 2 .AND. problem .EQ. 1 .AND. eig_switch .EQ. 0) THEN
        CALL printlog('Error mode extrapolation used for accelerating power iterations')
      END IF
      WRITE(amsg,101) "Scattering order:                                           ", scatt_ord
      CALL printlog(amsg)
      IF(quad_tpe == 1) THEN
        WRITE(amsg,101) "Level symmetric quadrature of order:                        ", quad_ord
        CALL printlog(amsg)
      ELSE IF(quad_tpe == 2) THEN
        WRITE(amsg,101) "Square Legendre-Chebychev quadrature of order:              ", quad_ord
        CALL printlog(amsg)
      ELSE IF(quad_tpe == 3) THEN
        WRITE(amsg,101) "Quadrature read from file. #Angles/octant:                  ",nangle
        CALL printlog(amsg)
        WRITE(amsg,'(2A)')"Quadrature file name:                                       ",quad_file
        CALL printlog(amsg)
      END IF
      WRITE(amsg,101) "Cross-section order:                                        ", xs_ord
      CALL printlog(amsg)
      WRITE(amsg,101) "Energy groups:                                              ", egmax
      CALL printlog(amsg)
      IF (problem == 0 .OR. (problem == 1 .AND. eig_switch == 0 ) ) THEN
        WRITE(amsg,101) "Maximum number of outer iterations:                         ", max_outer
        CALL printlog(amsg)
        WRITE(amsg,101) "Maximum number of inner iterations:                         ", max_inner
        CALL printlog(amsg)
        WRITE(amsg,102)   "Inner convergence criteria:                                 ", &
          inner_conv
        CALL printlog(amsg)
        IF(problem==1 .AND. eig_switch==0) THEN
          WRITE(amsg,102)   "Eigenvalue convergence &
            & criteria:                            ", k_conv
          CALL printlog(amsg)
        ENDIF
        WRITE(amsg,102)   "Outer convergence criteria:                                 ", outer_conv
        CALL printlog(amsg)
      ELSE
        WRITE(amsg,101) "Maximum number of newton iterations:                        ", max_outer
        CALL printlog(amsg)
        WRITE(amsg,101) "Maximum number of inner iterations(if used):                ", max_inner
        CALL printlog(amsg)
        WRITE(amsg,102) " Newton convergence criteria:                                  ", outer_conv
        CALL printlog(amsg)
        WRITE(amsg,102) "Inner convergence criteria (if used):                         ", inner_conv
        CALL printlog(amsg)
        WRITE(amsg,101) "Number of krylov iterations between restarts:               ",rd_restart
        CALL printlog(amsg)
        WRITE(amsg,101) "Maximum number of krylov iterations per newton iteration:   ",rd_max_kit
        CALL printlog(amsg)
      END IF
      IF(problem == 0)THEN
        WRITE(amsg,'(2A)') "File containing external source:                            ", source_filename
        CALL printlog(amsg)
      ENDIF

      IF(finflow /= 0 .AND. problem == 0)THEN
        WRITE(amsg,'(2A)') "File containing fixed inflow boundary conditions:           ", finflow_filename
        CALL printlog(amsg)
      END IF

      WRITE(amsg,'(2A)') "File containing cross-sections:                             ", cross_section_filename
      CALL printlog(amsg)
      WRITE(amsg,'(2A)') "File containing mesh:                                       ", mesh_filename
      CALL printlog(amsg)
      WRITE(amsg,'(2A)') "Flux output file:                                           ", flux_filename
      CALL printlog(amsg)

      IF(vtk_flux_output /= 0)THEN
        WRITE(amsg,'(2A)') "VTK-format flux output file:                                ", vtk_flux_filename
      CALL printlog(amsg)
      END IF
      IF(vtk_mat_output /= 0)THEN
        WRITE(amsg,'(2A)') "VTK-format material output file:                            ", vtk_mat_filename
      CALL printlog(amsg)
      END IF
      IF(vtk_reg_output /= 0)THEN
        WRITE(amsg,'(2A)') "VTK-format region output file:                              ", vtk_reg_filename
      CALL printlog(amsg)
      END IF
      IF(vtk_src_output /= 0)THEN
        WRITE(amsg,'(2A)') "VTK-format region output file:                              ", vtk_src_filename
      CALL printlog(amsg)
      END IF

      IF(inguess_flag /= 0)THEN
        WRITE(amsg,'(2A)') "Initial guess read from file:                               ", inguess_file
      CALL printlog(amsg)
      END IF

      IF(dump_flag /= 0)THEN
        WRITE(amsg,'(2A)') "Restart file:                                               ", dump_file
      CALL printlog(amsg)
      END IF
    END IF
    ! Formats

101 FORMAT(A60,I5)
102 FORMAT(A60,ES12.4)

  END SUBROUTINE echo_input

  SUBROUTINE read_input
    !***********************************************************************
    !
    ! This subroutine reads the standard input file
    !
    !***********************************************************************

    ! local variables

    CHARACTER(100) :: fname, tchar
    INTEGER :: rank,mpi_err,localunit,ios,legacy_v
    CALL GET_COMMAND_ARGUMENT(1,fname)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)
    localunit = rank+100

    OPEN(unit = localunit, file = fname , status = 'old', action = 'read',IOSTAT=ios)
    jobname=TRIM(ADJUSTL(fname))
    IF(ios .NE. 0)THEN
      CALL raise_fatal_error('error opening '//TRIM(fname))
    ENDIF
    IF(rank .EQ. 0)THEN
      CALL printlog("<><><><><><><><>"//TRIM(fname))
    ENDIF
    dump_file=TRIM(fname)//'_restart.out'
    flux_filename=TRIM(fname)//'_flux.out'
    vtk_flux_filename=TRIM(fname)//'_flux.vtk'
    vtk_mat_filename=TRIM(fname)//'_mat.vtk'
    vtk_reg_filename=TRIM(fname)//'_reg.vtk'
    vtk_src_filename=TRIM(fname)//'_src.vtk'
    cartesian_map_filename=TRIM(fname)//'_cartesian_map.out'
    converge_filename=TRIM(fname)//'_conv.convergence'

    legacy_v=-9999
    !determine if the input is yaml
    IF(fname(LEN(TRIM(fname))-4: LEN(TRIM(fname))) .EQ. ".yaml")THEN
      legacy_v=-1
    ELSE
      !determine if the input is legacy
      DO
        READ(localunit,'(A100)',IOSTAT=ios) tchar
        !legacy version 2
        IF(trim(adjustl(tchar)) .EQ. "start problem")legacy_v=0
        IF(ios .NE. 0)EXIT
      ENDDO
      REWIND(localunit)
    ENDIF

    SELECTCASE(legacy_v)
      CASE(-1)
        !yaml read case
        CALL yaml_read(localunit)
      CASE(0)
        !original THOR input version, for backwards compatibility
        CALL legacy_v0_read(localunit)
      CASE(-9999)
        !current THOR input format
        CALL inputfile_read(localunit)
      CASE DEFAULT
        CALL raise_fatal_error('invalid input format')
    ENDSELECT

    CLOSE(localunit)

  END SUBROUTINE read_input

END MODULE read_module
