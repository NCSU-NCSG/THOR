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
  USE adv_read_module

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

    execution = 1
    problem   = 1
    space_ord = 0
    finflow   = 0
    outer_acc = 1
    sweep_tpe = 1
    page_sweep= 0
    page_refl = 0
    page_iflw = 0
    max_outer = 5
    max_inner = 10
    inner_conv= 1.0E-4_d_t
    outer_conv= 1.0E-3_d_t
    k_conv    = 1.0E-4_d_t
    eig_switch= 0
    rd_restart= 25
    rd_max_kit= 250
    rd_method = 2
    inguess_flag =0
    dump_flag    =0
    ipow         =0
    print_conv   =0
    dfact_opt    =0

    source_filename        = "file.src"
    finflow_filename       = "file.bc"
    cross_section_filename = "file.xs"
    mesh_filename          = "file.mesh"
    flux_filename          = "file.flux"
    vtk_flux_filename      = "flux.vtk"
    quad_file              = "file.quad"
    dump_file              = "restart"
    inguess_file           = "inguess"
    vtk_mat_filename       = "mat.vtk"
    vtk_reg_filename       = "reg.vtk"
    vtk_src_filename       = "src.vtk"
    dens_fact_filename     = "density_factors.dat"
    print_xs_flag          = 1
    vtk_flux_output        = 0
    vtk_reg_output         = 0
    vtk_mat_output         = 0
    vtk_src_output         = 0

    quad_ord = 4
    quad_tpe = 1

    egmax = 1
    scatt_ord = 0
    xs_ord    = 0
    upscattering = 1
    multiplying  = 1
    scat_mult_flag=0

    glob_do_cartesian_mesh = .FALSE.
    cartesian_map_filename = "cartesian_map_output.dat"

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

    CHARACTER(100) :: buffer, fname
    CHARACTER(100000) :: regmap
    LOGICAL :: done
    INTEGER :: i, rank,mpi_err, localunit
    CALL GET_COMMAND_ARGUMENT(1,fname)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)
    localunit = rank+100

    OPEN(unit = localunit, file = fname , status = 'old', action = 'read')
    WRITE(*,*) "<><><><><><><><>", fname(LEN(TRIM(fname))-4: LEN(TRIM(fname)))

    IF(fname(LEN(TRIM(fname))-4: LEN(TRIM(fname))) .EQ. ".yaml") THEN
      WRITE(*,*) "ADV READ"
      CALL adv_read(localunit)
    ELSE
      ! read title


      READ(localunit,*) jobname

      ! main read loop

      done = .FALSE.
      DO WHILE(done .EQV. .FALSE.)

        READ(localunit,101,END=999) buffer
        IF     ( INDEX( lowercase(buffer) ,'start') > 0 .AND. INDEX( lowercase(buffer) ,'problem')>0 ) THEN
          CALL read_problem
        ELSE IF( INDEX( lowercase(buffer) ,'start') > 0 .AND. INDEX( lowercase(buffer) ,'inout')>0   ) THEN
          CALL read_inout
        ELSE IF( INDEX( lowercase(buffer) ,'start') > 0 .AND. INDEX( lowercase(buffer) ,'cross_sections')>0   ) THEN
          CALL read_cross_sections
        ELSE IF( INDEX( lowercase(buffer) ,'start') > 0 .AND. INDEX( lowercase(buffer) ,'quadrature')>0   ) THEN
          CALL read_quadrature_field
        ELSE IF( INDEX( lowercase(buffer) ,'start') > 0 .AND. INDEX( lowercase(buffer) ,'postprocess')>0   ) THEN
          CALL read_postprocess_field
        ELSE IF( INDEX( lowercase(buffer) ,'start') > 0 .AND. INDEX( lowercase(buffer) ,'regionmap')>0   ) THEN
          CALL read_regionmap_field(regmap)
        ELSE IF( INDEX( lowercase(buffer) ,'end') > 0   .AND. INDEX( lowercase(buffer) ,'file')   >0 ) THEN
          done=.TRUE.
        END IF

      END DO
    END IF

999 CONTINUE

    ! Call read_mesh to read tetrahedral mesh file

    CALL read_tetmesh

    ! Read reg2mat
    ALLOCATE(reg2mat(minreg:maxreg))
    IF(fname(LEN(TRIM(fname))-4: LEN(TRIM(fname))) .EQ. ".yaml") THEN
      reg2mat(0)=1
    ELSE
      READ(regmap,*) (reg2mat(i),i=minreg,maxreg)
    END IF

101 FORMAT(A100)

  END SUBROUTINE read_input

  SUBROUTINE read_regionmap_field(regmap)

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

  END SUBROUTINE read_regionmap_field

  SUBROUTINE read_quadrature_field

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

  END SUBROUTINE read_quadrature_field

  SUBROUTINE read_postprocess_field

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

  END SUBROUTINE read_postprocess_field

  SUBROUTINE read_cross_sections

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

  END SUBROUTINE read_cross_sections

  SUBROUTINE read_inout

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

  END SUBROUTINE read_inout

  SUBROUTINE read_problem

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
  END SUBROUTINE read_problem

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

END MODULE read_module
