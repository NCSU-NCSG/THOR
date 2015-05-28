#! /bin/bash
#

# get sparskit ready
if [ -f sparskit/SPARSKIT2/libskit.a ]; then
    rm sparskit/SPARSKIT2/libskit.a
fi
cd sparskit/SPARSKIT2
make
cd -
mv sparskit/SPARSKIT2/libskit.a .

execute="yes"
if [ $1 = "ifort" ]; then
  echo Using ifort as compiler. 
  cmple1="ifort -c -O3"
  cmple2="ifort -c -O3 -r8"
  cmple3="ifort -O3 -o"
elif [ $1 = "ifort_d" ]; then
  echo Using ifort with debug flags as compiler.
  cmple1="ifort -c -g -traceback -check all -fpe0"
  cmple2="ifort -c -g -traceback -check all -fpe0 -r8"
  cmple3="ifort -o"
elif [ $1 = "ifort_idbc" ]; then
  echo Using ifort for debugging with idbc.
  cmple1="ifort -c -debug -heap-arrays"
  cmple2="ifort -c -debug -heap-arrays"
  cmple3="ifort -debug  -heap-arrays  -o"
elif [ $1 = "ifort_h" ]; then
  echo Using ifort. Subroutine memory is allocated in heap, not on stack.
  cmple1="ifort -c -O3 -heap-arrays"
  cmple2="ifort -c -O3 -r8 -heap-arrays"
  cmple3="ifort -O3 -heap-arrays  -o"
elif [ $1 = "gfortran" ]; then
  echo Using gfortran as compiler.
  cmple1="gfortran -c -O3"
  cmple2="gfortran -c -O3 -fdefault-real-8"
  cmple3="gfortran -O3 -o"
elif [ $1 = "gfortran_d" ]; then
  echo Using gfortran with debugging flags as compiler.
  cmple1="gfortran -c -g -O0"
  cmple2="gfortran -c -g -O0 -fdefault-real-8"
  cmple3="gfortran -O0 -o"
elif [ $1 = "gfortran_p" ]; then
  echo Using gfortran as compiler. Setting profiling flags.
  cmple1="gfortran -c -pg"
  cmple2="gfortran -c -pg -fdefault-real-8"
  cmple3="gfortran -pg -o"
elif [ $1 = "mpifort" ]; then
  echo Using mpifort as compiler.
  cmple1="mpifort -c -O3 -fdefault-real-8"
  cmple2="mpifort -c -O3 -fdefault-real-8"
  cmple3="mpifort -O3 -o"
else
   echo Compiler unknown.
   execute="no"
fi

if [ $execute = "yes" ]; then
  $cmple1 precmod.f90
  $cmple1 stringmod.f90
  $cmple1 types.f90
  $cmple1 parameter_types.f90
  $cmple1 filename_types.f90
  $cmple1 vector_types.f90
  $cmple1 general_utility_module.f90
  $cmple1 cross_section_types.f90
  $cmple1 geometry_types.f90
  $cmple1 angle_types.f90
  $cmple1 multindex_types.f90
  $cmple1 global_variables.f90
  $cmple1 termination_module.f90
  $cmple1 input_check.f90
  $cmple1 quadrature_module.f90
  $cmple1 read_cross_section_module.f90
  $cmple1 readmesh_module.f90
  $cmple1 read_source_module.f90
  $cmple1 read_inflow_module.f90
  $cmple1 read_module.f90
  $cmple1 dump_inguess_module.f90
  $cmple1 ahotc_matrix_module.f90
  $cmple1 sph_harmonics_module.f90
  $cmple1 transport_kernel_module_CCE.f90
  $cmple2 transport_kernel_module_SC.f90
  $cmple2 transport_kernel_module_LC.f90
  $cmple1 cell_splitting_module.f90
  $cmple1 sweep_module.f90
  $cmple1 setup_module.f90
  $cmple1 inner_iteration_module.f90
  $cmple1 wrapup_module.f90
  $cmple1 outer_iteration_module.f90
  $cmple1 sweep_module.f90
  $cmple1 jfnk_module.f90 
  $cmple1 solver_module.f90
  $cmple1 execution_module.f90
  $cmple1 distdot.f90  
  $cmple1 main.f90
  $cmple3 thor-1.0.exe  *.o libskit.a 
  mv thor-1.0.exe ../

#  rm transport_kernel_module.f90 read_module.f90
  rm *.o *.mod
  rm libskit.a
fi
#
