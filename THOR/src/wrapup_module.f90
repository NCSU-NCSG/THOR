MODULE wrapup_module
  !***********************************************************************
  !
  ! Wraup module contains all subroutines to finish up problem and echo
  ! output
  !
  !***********************************************************************

  ! Use derived-type modules

  USE types
  USE parameter_types
  USE filename_types
  USE vector_types
  USE cross_section_types
  USE geometry_types
  USE angle_types
  USE globals
  USE error_module
  USE general_utility_module

  IMPLICIT NONE

CONTAINS

  SUBROUTINE wrapup(flux, keff, unit_number, suffix, is_final)
    !**********************************************************************
    !
    ! Subroutine wrapup calls routines that print output into file(s)
    !
    !**********************************************************************

    ! Pass eigenvalue

    REAL(kind=d_t), INTENT(in) :: keff

    ! Pass scalar flux

    REAL(kind=d_t) :: flux(num_moments_v,namom,num_cells,egmax,niter)

    ! Pass optional arguments

    INTEGER(kind=li) :: unit_number
    CHARACTER(*) :: suffix
    LOGICAL :: is_final

    ! Declare temporary variables

    INTEGER(kind=li) :: i, j, eg, l, region, ix, iy, iz,mat_indx
    REAL(kind=d_t) :: reac_rates(4,minreg:maxreg,egmax+1)
    REAL(kind=d_t) :: reg_volume(minreg:maxreg)
    REAL(kind=d_t), ALLOCATABLE :: cartesian_map(:,:,:,:)
    REAL(kind=d_t) :: centroid(3)
    TYPE(vector) :: vertex
    REAL(kind=d_t) :: delx, dely, delz, cartesian_vol
    TYPE(vector) :: v0, v1, v2, v3, point, barycentric
    INTEGER(kind=li) :: point_flux_location_element_indices(number_point_flux_locations)

    ! Print runtime
    IF (rank .EQ. 0) THEN
      WRITE(unit_number,*)
      WRITE(unit_number,'(A)') "--------------------------------------------------------"
      WRITE(unit_number,'(A)') "   Execution Summary   "
      WRITE(unit_number,'(A)') "--------------------------------------------------------"
      IF (conv_flag==0) THEN
        WRITE(unit_number,'(A)') "Warning! Execution finished without satisfying all stopping criteria. Warning!"
      ELSE IF(conv_flag==1) THEN
        WRITE(unit_number,'(A)') "Execution finished successfully. All stopping criteria satisfied."
      END IF
      WRITE(unit_number,202) "Runtime (seconds):                                         ", finish-start
      IF (problem == 0 .OR. (problem == 1 .AND. eig_switch == 0 ) ) THEN
        WRITE(unit_number,101) "Number of outer iterations:                                   ", outer
        WRITE(unit_number,101) "Number of inner iterations:                                   ", tot_nInners
      ELSE
        WRITE(unit_number,101) "Number of newton iterations:                                ", nit
        IF (rd_method==1) THEN
          WRITE(unit_number,101) "Number of inner iterations:                                 ", tot_nInners
        END IF
        WRITE(unit_number,101)   "Total number of krylov iterations:                          ", tot_kit
        WRITE(unit_number,*)
      END IF

      IF(problem == 1)THEN
        WRITE(unit_number,102) "Final eigenvalue:                                          ", keff
      END IF
      IF(problem==0) THEN
        WRITE(unit_number,102) "Maximum scalar group flux error (outer iteration):         ", max_outer_error
        WRITE(unit_number,'(A)') "Maximum scalar flux error by group:"
        DO eg=1, egmax
          WRITE(unit_number,103) eg,max_error(eg)
        END DO
      ELSE IF(problem==1 .AND. eig_switch == 1) THEN
        WRITE(unit_number,'(A)') "Maximum residual by group:"
        DO eg=1, egmax
          WRITE(unit_number,103) eg,max_error(eg)
        END DO
      ELSE IF(problem==1 .AND. eig_switch == 0) THEN
        WRITE(unit_number,102) "Eigenvalue error:                                          ", k_error
        WRITE(unit_number,102) "Maximum fission source error:                              ", f_error
        WRITE(unit_number,102) "Maximum scalar group flux error:                           ", max_outer_error
      END IF
    END IF

    ! Formats

101 FORMAT(A,I0)
102 FORMAT(A,ES25.16)
202 FORMAT(A,ES25.16)
103 FORMAT(I5,ES25.16)

    ! Open angular flux file and write volume and face flux in file

    IF(vtk_flux_output /= 0 .AND. rank .EQ. 0)THEN
      OPEN(unit=10,file=TRIM(vtk_flux_filename)//TRIM(suffix),status='unknown',action='write')

      WRITE(10,'(a26)') '# vtk DataFile Version 3.0'
      WRITE(10,'(a72)') jobname
      WRITE(10,'(a5)') 'ASCII'
      WRITE(10,'(a25)') 'DATASET UNSTRUCTURED_GRID'
      WRITE(10,'(a6,1x,i12,1x,a5)') 'POINTS',num_vert,'float'

      DO i=1, num_vert
        WRITE(10,'(3(1x,es12.5))') vertices(i)%v%x1,vertices(i)%v%x2,&
              vertices(i)%v%x3
      END DO

      WRITE(10,*)
      WRITE(10,'(a5,1x,i12,1x,i12)') 'CELLS',num_cells,num_cells+&
            4*num_cells

      DO i=1, num_cells
        WRITE(10,'(5(i12,1x))') 4,cells(i)%R(0)-1,cells(i)%R(1)-1,&
              cells(i)%R(2)-1,cells(i)%R(3)-1
      END DO

      WRITE(10,*)
      WRITE(10,'(a11,1x,i12)') 'CELL_TYPES',num_cells

      DO i=1, num_cells
        WRITE(10,'(i12)') 10
      END DO

      WRITE(10,*)
      WRITE(10,'(a9,1x,i12)') 'CELL_DATA', num_cells
      WRITE(10,'(a5,1x,a16,1x,i12)') 'FIELD','Neutronics_Edits',egmax

      l=1
      DO eg=1, egmax
        WRITE(10,'(i12,1x,i12,1x,i12,1x,a5)') eg,1,num_cells,'float'
        DO i=1, num_cells
          WRITE(10,'(es12.5)') flux(1,1,i,eg,niter)
        END DO
      END DO

      CLOSE(10)

    END IF

    ! Open flux file and write volume and cell flux in file
    IF (rank .EQ. 0) THEN
      OPEN(unit=20,file=TRIM(flux_filename)//TRIM(suffix),status='unknown',action='write')

      l=1
      WRITE(20,'(I0)') num_cells
      DO i=1, num_cells
        WRITE(20,'(ES24.16)', ADVANCE='NO') cells(i)%volume
        DO eg=1, egmax
          WRITE(20,'(ES24.16)', ADVANCE='NO') flux(1,1,i,eg,niter)
        END DO
        WRITE(20,*)
      END DO

      CLOSE(20)

      WRITE(unit_number,'(A)') "A flux file was been succesfully written!"
    END IF

    ! Compute region averaged fluxes and reaction rates

    IF (rank .EQ. 0) THEN
      WRITE(unit_number,*)
      WRITE(unit_number,'(A)') "--------------------------------------------------------"
      WRITE(unit_number,'(A)') "   Region averaged reaction rates  "
      WRITE(unit_number,'(A)') "--------------------------------------------------------"
      WRITE(unit_number,*)
      reg_volume=0.0_d_t
      reac_rates=0.0_d_t
      DO i=1,num_cells
        mat_indx=material_ids(reg2mat(cells(i)%reg))
        DO eg=1,egmax
          IF(eg.EQ.1) reg_volume(cells(i)%reg)=reg_volume(cells(i)%reg)+ cells(i)%volume
          reac_rates(1,cells(i)%reg,eg)=reac_rates(1,cells(i)%reg,eg)  + cells(i)%volume * &
                flux(1,1,i,eg,niter)
          reac_rates(2,cells(i)%reg,eg)=reac_rates(2,cells(i)%reg,eg)  + cells(i)%volume  &
            * xs_mat(mat_indx)%sigma_f(eg)* &
                flux(1,1,i,eg,niter)
          reac_rates(3,cells(i)%reg,eg)=reac_rates(3,cells(i)%reg,eg)  + cells(i)%volume *  &
                (xs_mat(mat_indx)%sigma_t(eg) &
                - xs_mat(mat_indx)%tsigs(eg) )         * &
                flux(1,1,i,eg,niter)
          reac_rates(4,cells(i)%reg,eg)=reac_rates(4,cells(i)%reg,eg)  + cells(i)%volume *  &
                xs_mat(mat_indx)%nusig_f(eg)* &
                flux(1,1,i,eg,niter)

        END DO
      END DO
      DO eg=1,egmax
        DO region=minreg,maxreg
          reac_rates(1,region,egmax+1)=reac_rates(1,region,egmax+1)+reac_rates(1,region,eg)
          reac_rates(2,region,egmax+1)=reac_rates(2,region,egmax+1)+reac_rates(2,region,eg)
          reac_rates(3,region,egmax+1)=reac_rates(3,region,egmax+1)+reac_rates(3,region,eg)
          reac_rates(4,region,egmax+1)=reac_rates(4,region,egmax+1)+reac_rates(4,region,eg)
        END DO
      END DO
      DO eg=1,egmax+1
        DO region=minreg,maxreg
          reac_rates(1,region,eg)=reac_rates(1,region,eg)/reg_volume(region)
          reac_rates(2,region,eg)=reac_rates(2,region,eg)/reg_volume(region)
          reac_rates(3,region,eg)=reac_rates(3,region,eg)/reg_volume(region)
          reac_rates(4,region,eg)=reac_rates(4,region,eg)/reg_volume(region)
        END DO
      END DO
      DO region=minreg,maxreg
        WRITE(unit_number,501) '-- Region --',region,' -- Material -- ', &
          TRIM(xs_mat(material_ids(reg2mat(region)))%mat_name),' Volume = ',reg_volume(region)
        WRITE(unit_number,502) '   Group          Flux       Fission    Absorption      Fiss Src'
        DO eg=1,egmax
          WRITE(unit_number,503) eg,reac_rates(1,region,eg),reac_rates(2,region,eg),&
                reac_rates(3,region,eg),reac_rates(4,region,eg)
        END DO
        WRITE(unit_number,504) '   Total',reac_rates(1,region,egmax+1),reac_rates(2,region,egmax+1),&
              reac_rates(3,region,egmax+1),reac_rates(4,region,egmax+1)
        WRITE(unit_number,*)
      END DO
    END IF

    ! if desired print cartesian flux map to file

    IF (rank .EQ. 0 .AND. glob_do_cartesian_mesh .AND. is_final) THEN

      ! compute reaction rates for each cartesian cell
      ! computes flux(1), total(2), absorption(3), total scattering (4),
      ! fission(5), fission production(6)

      ! allocate the array holding the data
      ALLOCATE(cartesian_map(6, glob_cmap_nx, glob_cmap_ny, glob_cmap_nz))

      ! compute the spacing
      delx = (glob_cmap_max_x - glob_cmap_min_x) / REAL(glob_cmap_nx)
      dely = (glob_cmap_max_y - glob_cmap_min_y) / REAL(glob_cmap_ny)
      delz = (glob_cmap_max_z - glob_cmap_min_z) / REAL(glob_cmap_nz)
      cartesian_vol = delx * dely * delz

      ! make sure it's properly initialized
      cartesian_map = zero

      ! computation of the averaged reaction rates (only energy integrated!)
      DO eg = 1, egmax
        DO i = 1, num_cells
          mat_indx=material_ids(reg2mat(cells(i)%reg))

          ! this is the current cell with index i, first we need to find out
          ! which x/y/z cartesian this cell belongs to
          centroid = zero

          DO l = 0, 3
            vertex = vertices(cells(i)%R(l))%v
            centroid(1) = centroid(1) + vertex%x1
            centroid(2) = centroid(2) + vertex%x2
            centroid(3) = centroid(3) + vertex%x3
          END DO
          centroid = centroid * fourth

          ! get the correct entry in the cartesian_map array but note that
          ! numbering starts from 1 so we need to add one
          ix = FLOOR((centroid(1) - glob_cmap_min_x) / delx) + 1
          iy = FLOOR((centroid(2) - glob_cmap_min_y) / dely) + 1
          iz = FLOOR((centroid(3) - glob_cmap_min_z) / delz) + 1

          ! we have to make sure that the x/y/z coordinates are permissible
          IF (ix .GT. 0 .AND. ix .LE. glob_cmap_nx .AND.&
                iy .GT. 0 .AND. iy .LE. glob_cmap_ny .AND.&
                iz .GT. 0 .AND. iz .LE. glob_cmap_nz) THEN

            ! Now go through all reaction types and accumulate into cartesian_map
            ! 1. the scalar flux
            cartesian_map(1, ix, iy, iz) = cartesian_map(1, ix, iy, iz) + &
                  cells(i)%volume * flux(1, 1, i, eg, niter) / cartesian_vol
            ! 2. total interaction rate
            cartesian_map(2, ix, iy, iz) = cartesian_map(2, ix, iy, iz) + &
                  cells(i)%volume * flux(1, 1, i, eg, niter) * &
                  xs_mat(mat_indx)%sigma_t(eg) / cartesian_vol
            ! 3. absorption rate
            cartesian_map(3, ix, iy, iz) = cartesian_map(3, ix, iy, iz) + &
                  cells(i)%volume * flux(1, 1, i, eg, niter) * &
                  (xs_mat(mat_indx)%sigma_t(eg) &
                  - xs_mat(mat_indx)%tsigs(eg)) / &
                  cartesian_vol
            ! 4. total scattering rate
            cartesian_map(4, ix, iy, iz) = cartesian_map(4, ix, iy, iz) + &
                  cells(i)%volume * flux(1, 1, i, eg, niter) * &
                  xs_mat(mat_indx)%tsigs(eg) / cartesian_vol
            ! 5. fission rate
            cartesian_map(5, ix, iy, iz) = cartesian_map(5, ix, iy, iz) + &
                  cells(i)%volume * flux(1, 1, i, eg, niter) * &
                  xs_mat(mat_indx)%sigma_f(eg) / cartesian_vol
            ! 6. fission source rate
            cartesian_map(6, ix, iy, iz) = cartesian_map(6, ix, iy, iz) + &
                  cells(i)%volume * flux(1, 1, i, eg, niter) * &
                  xs_mat(mat_indx)%nusig_f(eg) / &
                  cartesian_vol
          END IF
        END DO
      END DO

      ! write it to file
      OPEN(unit=30, file = TRIM(cartesian_map_filename), status='unknown', action='write')
      WRITE(30, '(A)') "     ix    iy    iz             flux             total&
            &        absorption        scattering           fission&
            &       fission src"
      DO iz = 1, glob_cmap_nz
        DO iy = 1, glob_cmap_ny
          DO ix = 1, glob_cmap_nx
            WRITE(30, 499) ix, iy, iz, cartesian_map(1, ix, iy, iz),&
                  cartesian_map(2, ix, iy, iz), cartesian_map(3, ix, iy, iz),&
                  cartesian_map(4, ix, iy, iz), cartesian_map(5, ix, iy, iz),&
                  cartesian_map(6, ix, iy, iz)
          END DO
        END DO
      END DO

      CLOSE(30)

      ! wrap up
      DEALLOCATE(cartesian_map)
    END IF

    ! fluxes at point locations
    IF (rank .EQ. 0 .AND. number_point_flux_locations .gt. 0 .AND. is_final) THEN
      point_flux_location_element_indices = 0_li
      ! find the right tets
      DO i = 1, num_cells
        v0 = vertices(cells(i)%R(0))%v
        v1 = vertices(cells(i)%R(1))%v
        v2 = vertices(cells(i)%R(2))%v
        v3 = vertices(cells(i)%R(3))%v

        DO j = 1, number_point_flux_locations
          point%x1 = point_flux_locations(j, 1)
          point%x2 = point_flux_locations(j, 2)
          point%x3 = point_flux_locations(j, 3)

          ! point_flux_location_element_indices
          CALL barycentricCoordinates(point, v0, v1, v2, v3, barycentric)

          ! is the point in the tetrahedron?
          IF (barycentric%x1 .ge. 0.0_d_t .and. barycentric%x1 .le. 1.0_d_t .and. &
              barycentric%x2 .ge. 0.0_d_t .and. barycentric%x2 .le. 1.0_d_t .and. &
              barycentric%x3 .ge. 0.0_d_t .and. barycentric%x3 .le. 1.0_d_t .and. &
              1.0_d_t - barycentric%x1 - barycentric%x2 - barycentric%x3 .ge. 0.0_d_t .and. &
              1.0_d_t - barycentric%x1 - barycentric%x2 - barycentric%x3 .le. 1.0_d_t) THEN
              point_flux_location_element_indices(j) = i
          END IF
        END DO
      END DO
      ! print the group averages fluxes
      WRITE(unit_number, *)
      WRITE(unit_number, '(A)') "--------------------------------------------------------"
      WRITE(unit_number, '(A)') "   Point flux values  "
      WRITE(unit_number, '(A)') "--------------------------------------------------------"
      WRITE(unit_number, *)
      DO j = 1, number_point_flux_locations
        IF (point_flux_location_element_indices(j) .GT. 0_li) THEN
          WRITE(unit_number, 510) "Point: ", point_flux_locations(j, :)
          DO eg=1, egmax
            WRITE(unit_number, 511) "Group: ", eg, flux(1, 1, point_flux_location_element_indices(j), eg, niter)
          END DO
        ELSE
          WRITE(unit_number, 510) "Warning. Point: ", point_flux_locations(j, :), " was not found."
        CALL raise_warning("A flux point was not found.")
        END IF
      END DO
    END IF



    ! write a csv output file containing all region averaged information
    OPEN(unit=20, file=TRIM(jobname)//TRIM('_out.csv'), status='unknown', action='write')
    IF(problem == 1)WRITE(20,'(A,F20.16)')'k-eff Eigenvalue: ',keff
    WRITE(20, 502, ADVANCE = "NO") "Region, Material, "
    DO eg = 1, egmax - 1
      WRITE(20, 505, ADVANCE = "NO") " flux g = ", eg, ","
    END DO
    WRITE(20, 506) " flux g = ", eg

    DO region=minreg,maxreg
      WRITE(20, 509, ADVANCE = "NO") region, ", ", &
        TRIM(xs_mat(material_ids(reg2mat(region)))%mat_name),','
      DO eg= 1, egmax - 1
        WRITE(20, 507, ADVANCE = "NO") reac_rates(1, region, eg), ","
      END DO
      WRITE(20, 508) reac_rates(1, region, egmax)
    END DO
    CLOSE(20)

499 FORMAT(3I6,6ES18.8)
501 FORMAT(A,I4,3A,ES14.6)
502 FORMAT(A)
503 FORMAT(I8,4ES14.6)
504 FORMAT(A8,4ES14.6)
505 FORMAT(A9,I3,A1)
506 FORMAT(A9,I3)
507 FORMAT(ES25.16,A1)
508 FORMAT(ES25.16)
509 FORMAT(I3,3A)
510 FORMAT(A,3ES14.6)
511 FORMAT(A,I3,ES14.6)
  END SUBROUTINE wrapup

END MODULE wrapup_module
