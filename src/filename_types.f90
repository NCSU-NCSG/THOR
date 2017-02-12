module filename_types
!***********************************************************************
! These are internal file names, actual names are specified in input
!***********************************************************************

  implicit none

  integer, parameter :: length = 80
  character(length) jobname, source_filename, finflow_filename, &
       cross_section_filename, mesh_filename, flux_filename,    &
       vtk_flux_filename,quad_file,dump_file,inguess_file,      &
       vtk_mat_filename,vtk_reg_filename,vtk_src_filename,      &
       dens_fact_filename, cartesian_map_filename

end module filename_types
