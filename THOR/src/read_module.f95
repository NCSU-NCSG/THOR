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
  USE global_variables

  ! Use modules that pertain setting up problem

  USE read_cross_section_module
  USE readmesh_module
  USE read_source_module
  USE read_inflow_module
  USE quadrature_module
  USE check_input
  USE termination_module
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
    INTEGER ::rank,mpi_err, localunit
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)

    ! Set the defaults

    CALL  set_default

    ! Call read_input to read input file

    CALL read_input

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

    IF(execution .EQ. 0) CALL stop_thor(1_li)

    ! Call generate_num_moments and generate_multi_index to create indices

    CALL generate_num_moments(space_ord,num_moments_v,num_moments_f)

    ALLOCATE(index_v(num_moments_v),index_f(num_moments_f),&
          stat=alloc_stat)
    IF(alloc_stat /= 0) CALL stop_thor(2_li)

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
        CALL stop_thor(21_li)
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
        CALL stop_thor(21_li)
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
        CALL stop_thor(21_li)
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
        CALL stop_thor(21_li)
      END IF
    END IF

  END SUBROUTINE generate_multi_index

  SUBROUTINE set_default

    ! Set the input default

    execution               = 1
    problem                 = 1
    space_ord               = 0
    finflow                 = 0
    outer_acc               = 1
    sweep_tpe               = 1
    page_sweep              = 0
    page_refl               = 0
    page_iflw               = 0
    max_outer               = 5
    max_inner               = 10
    inner_conv              = 1.0E-4_d_t
    outer_conv              = 1.0E-3_d_t
    k_conv                  = 1.0E-4_d_t
    eig_switch              = 0
    rd_restart              = 25
    rd_max_kit              = 250
    rd_method               = 2
    inguess_flag            = 0
    dump_flag               = 0
    ipow                    = 0
    print_conv              = 0
    dfact_opt               = 0

    source_filename         = "source.dat"
    finflow_filename        = "finflow.dat"
    cross_section_filename  = "xs.dat"
    mesh_filename           = "mesh.thrm"
    flux_filename           = "flux.out"
    vtk_flux_filename       = "flux.vtk"
    quad_file               = "quad.dat"
    dump_file               = "restart.out"
    inguess_file            = "initial_guess.dat"
    vtk_mat_filename        = "mat.vtk"
    vtk_reg_filename        = "reg.vtk"
    vtk_src_filename        = "src.vtk"
    dens_fact_filename      = "density_factor.dat"
    print_xs_flag           = 0
    vtk_flux_output         = 0
    vtk_reg_output          = 0
    vtk_mat_output          = 0
    vtk_src_output          = 0

    quad_ord                = 4
    quad_tpe                = 1

    egmax                   = 1
    scatt_ord               = 0
    xs_ord                  = 0
    upscattering            = 1
    multiplying             = 1
    scat_mult_flag          = 0

    glob_do_cartesian_mesh = .FALSE.
    cartesian_map_filename = "cartesian_map.out"

    number_point_flux_locations = 0_li

  END SUBROUTINE set_default

  SUBROUTINE check_standard_input

    ! Check the input for errors

    ! problem must be 0 or 1
    IF(problem > 1 .OR. problem < 0)THEN
      CALL stop_thor(8_li)
    ENDIF

    ! jfnk parameters
    IF(eig_switch==1) THEN
      IF (rd_restart<0) THEN
        CALL stop_thor(18_li)
      END IF
      IF (rd_max_kit<rd_restart) THEN
        CALL stop_thor(19_li)
      END IF
      IF (rd_method<1 .OR. rd_method >3) THEN
        CALL stop_thor(20_li)
      END IF
    END IF

    ! in case problem == 1, no inflow file allowed

    IF(problem .EQ. 1 .AND. finflow .EQ. 1) THEN
      CALL stop_thor(29_li)
    END IF

    ! in case problem == 1, no vtk source file allowed

    IF(problem .EQ. 1 .AND. vtk_src_output .EQ. 1) THEN
      CALL stop_thor(30_li)
    END IF


  END SUBROUTINE check_standard_input

  SUBROUTINE echo_input

    ! Echo out problem specifications

    INTEGER :: i

    IF (rank .EQ. 0) THEN
      WRITE(6,*)
      WRITE(6,*) "--------------------------------------------------------"
      WRITE(6,*) "   Input Summary  "
      WRITE(6,*) "--------------------------------------------------------"

      IF(problem == 1)THEN
        IF ( eig_switch == 0 ) THEN
          WRITE(6,*) "Eigenvalue calculation using PI selected"
        ELSE
          WRITE(6,*) "Eigenvalue calculation using JFNK selected"
          IF(rd_method == 1) THEN
            WRITE(6,*) "Method: F(u) is evaluated using one outer iteration with lagged upscattering."
          ELSE IF (rd_method == 2) THEN
            WRITE(6,*) "Method: F(u) is evaluated using flat iteration scheme."
          ELSE
            WRITE(6,*) "Method: F(u) is evaluated using flat iteration scheme wit updated downscattering."
          END IF
        END IF
      ELSE IF(problem == 0)THEN
        WRITE(6,*) "External source calculation selected"
      END IF

      ! Print rest of input

      WRITE(6,*)   "Problem Title:                                              ", jobname
      WRITE(6,101) "Spatial order:                                              ", space_ord
      IF      (sweep_tpe .EQ. 1) THEN
        WRITE(6,*) 'Precomputed mesh sweep is used.'
      END IF
      IF      (outer_acc .EQ. 1 .AND. problem .EQ. 1 .AND. eig_switch .EQ. 0) THEN
        WRITE(6,*) 'Power iterations are not accelerated'
      ELSE IF (outer_acc .EQ. 2 .AND. problem .EQ. 1 .AND. eig_switch .EQ. 0) THEN
        WRITE(6,*) 'Error mode extrapolation used for accelerating power iterations'
      END IF
      WRITE(6,101) "Scattering order:                                           ", scatt_ord
      IF(quad_tpe == 1) THEN
        WRITE(6,101) "Level symmetric quadrature of order:                        ", quad_ord
      ELSE IF(quad_tpe == 2) THEN
        WRITE(6,101) "Square Legendre-Chebychev quadrature of order:              ", quad_ord
      ELSE IF(quad_tpe == 3) THEN
        WRITE(6,101) "Quadrature read from file. #Angles/octant:                  ",nangle
        WRITE(6,*)   "Quadrature file name:                                       ",quad_file
      END IF
      WRITE(6,101) "Cross-section order:                                        ", xs_ord
      WRITE(6,101) "Energy groups:                                              ", egmax
      IF (problem == 0 .OR. (problem == 1 .AND. eig_switch == 0 ) ) THEN
        WRITE(6,101) "Maximum number of outer iterations:                         ", max_outer
        WRITE(6,101) "Maximum number of inner iterations:                         ", max_inner
        WRITE(6,102)   " Inner convergence criteria:                                 ", inner_conv
        IF(problem==1 .AND. eig_switch==0) WRITE(6,102)   " Eigenvalue convergence criteria:                            ", k_conv
        WRITE(6,102)   " Outer convergence criteria:                                 ", outer_conv
      ELSE
        WRITE(6,101) "Maximum number of newton iterations:                        ", max_outer
        WRITE(6,101) "Maximum number of inner iterations(if used):                ", max_inner
        WRITE(6,102) " Newton convergence criteria:                                  ", outer_conv
        WRITE(6,102) " Inner convergence criteria (if used):                         ", inner_conv
        WRITE(6,101) "Number of krylov iterations between restarts:               ",rd_restart
        WRITE(6,101) "Maximum number of krylov iterations per newton iteration:   ",rd_max_kit
      END IF
      IF(problem == 0)THEN
        WRITE(6,*) "File containing external source:                            ", source_filename
      ENDIF

      IF(finflow /= 0 .AND. problem == 0)THEN
        WRITE(6,*) "File containing fixed inflow boundary conditions:           ", finflow_filename
      END IF

      WRITE(6,*) "File containing cross-sections:                             ", cross_section_filename
      WRITE(6,*) "File containing mesh:                                       ", mesh_filename
      WRITE(6,*) "Flux output file:                                           ", flux_filename

      IF(vtk_flux_output /= 0)THEN
        WRITE(6,*) "VTK-format flux output file:                                ", vtk_flux_filename
      END IF
      IF(vtk_mat_output /= 0)THEN
        WRITE(6,*) "VTK-format material output file:                            ", vtk_mat_filename
      END IF
      IF(vtk_reg_output /= 0)THEN
        WRITE(6,*) "VTK-format region output file:                              ", vtk_reg_filename
      END IF
      IF(vtk_src_output /= 0)THEN
        WRITE(6,*) "VTK-format region output file:                              ", vtk_src_filename
      END IF

      IF(inguess_flag /= 0)THEN
        WRITE(6,*) "Initial guess read from file:                               ", inguess_file
      END IF

      IF(dump_flag /= 0)THEN
        WRITE(6,*) "Restart file:                                               ", dump_file
      END IF
    END IF
    ! Formats

101 FORMAT(1X,A60,I5)
102 FORMAT(A60,ES12.4)
202 FORMAT(1X,A60,ES12.4)
103 FORMAT(1X,I5,ES12.4)
104 FORMAT(1X,I14,A,I14)

  END SUBROUTINE echo_input

  SUBROUTINE read_input
    !***********************************************************************
    !
    ! This subroutine reads the standard input file
    !
    !***********************************************************************

    ! local variables

    CHARACTER(100) :: fname, tchar
    CHARACTER(100000) :: regmap
    INTEGER :: i,rank,mpi_err,localunit,ios,legacyv
    CALL GET_COMMAND_ARGUMENT(1,fname)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)
    localunit = rank+100

    OPEN(unit = localunit, file = fname , status = 'old', action = 'read')
    IF(rank .EQ. 0)THEN
      WRITE(*,*) "<><><><><><><><>", fname
    ENDIF

    legacyv=0
    !determine if the input is yaml
    IF(fname(LEN(TRIM(fname))-4: LEN(TRIM(fname))) .EQ. ".yaml")THEN
      legacyv=-1
    ELSE
      !determine if the input is legacy
      DO
        READ(localunit,'(A100)',IOSTAT=ios) tchar
        !legacy version 2
        IF(trim(adjustl(tchar)) .EQ. "start problem")legacyv=1
        IF(ios .NE. 0)EXIT
      ENDDO
      REWIND(localunit)
    ENDIF

    SELECTCASE(legacyv)
      CASE(-1)
        !yaml read case
        CALL yaml_read(localunit)
      CASE(1)
        !original THOR input version, for backwards compatibility
        CALL legacyv1_read(localunit,regmap)
      CASE(0)
        !current THOR input format
        CALL inputfile_read(localunit)
      CASE DEFAULT
        STOP 'invalid input format'
    ENDSELECT

    CLOSE(localunit)

    ! Call read_mesh to read tetrahedral mesh file
    CALL read_tetmesh

    ! Read reg2mat
    ALLOCATE(reg2mat(minreg:maxreg))
    IF(fname(LEN(TRIM(fname))-4: LEN(TRIM(fname))) .EQ. ".yaml") THEN
      reg2mat(0)=1
    ELSE
      READ(regmap,*) (reg2mat(i),i=minreg,maxreg)
    END IF

  END SUBROUTINE read_input

END MODULE read_module
