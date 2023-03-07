!THOR is licensed under the MIT License.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!> @brief Filename types.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
MODULE filename_types
  !***********************************************************************
  ! These are internal file names, actual names are specified in input
  !***********************************************************************

  IMPLICIT NONE

  INTEGER, PARAMETER :: length = 100
  CHARACTER(length) jobname, source_filename, finflow_filename, &
        cross_section_filename, mesh_filename, flux_filename,    &
        vtk_flux_filename,quad_file,dump_file,inguess_file,      &
        vtk_mat_filename,vtk_reg_filename,vtk_src_filename,      &
        dens_fact_filename, cartesian_map_filename, converge_filename

END MODULE filename_types
