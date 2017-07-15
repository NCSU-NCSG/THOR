#Test THOR Makefile using gcc automatic rule development (-M -cpp)

# List of all variables

#compiler
	CC = mpif90

#options
	OPT1 = -c -O3 -g
	OPT2 = -fdefault-real-8 
	OPT3 = -O3 -g -o

# Setting up file names
	SRC = $(wildcard *.f90)
	# If there are any source files in the src folder that 
	# need to be ignored for compilation, 
	# add them here (integrate cube a b ..., $SRC))
	SRC1 = $(filter-out integrate_cube.f90, $(SRC))
	OBJ = $(SRC1:.f90=.o)
	ARC = libskit.a

#Executable Name
	EXE = thor-1.0.exe

#Begin Global Builds

all: $(EXE)

$(EXE): $(OBJ) $(ARC)
	$(CC) $(OPT3) $(EXE)  $(OBJ) $(ARC) 
	@echo $(CC) $(OPT3) $(EXE) '*.o' $(ARC) 
	@mv ./$(EXE) ../ 

clean: 
	rm *.mod *.o *.a

libskit.a: 
	cd ./sparskit/SPARSKIT2; make
	mv sparskit/SPARSKIT2/libskit.a .

#Automatic Items Shown Below
#-------------------------------------------------------------------------------
ahotc_matrix_module.o ahotc_matrix_module.mod: ahotc_matrix_module.f90 \
 types.mod parameter_types.mod multindex_types.mod
	$(CC) $(OPT1) ahotc_matrix_module.f90
#-------------------------------------------------------------------------------
angle_types.o angle_types.mod: angle_types.f90 types.mod vector_types.mod
	$(CC) $(OPT1) angle_types.f90
#-------------------------------------------------------------------------------
cell_splitting_module.o cell_splitting_module.mod: \
 cell_splitting_module.f90 types.mod parameter_types.mod \
 filename_types.mod vector_types.mod cross_section_types.mod \
 geometry_types.mod angle_types.mod multindex_types.mod \
 global_variables.mod general_utility_module.mod \
 transport_kernel_module_sc.mod transport_kernel_module_lc.mod \
 transport_kernel_module_cce.mod termination_module.mod
	$(CC) $(OPT1) cell_splitting_module.f90
#-------------------------------------------------------------------------------
cross_section_types.o cross_section_types.mod: cross_section_types.f90 \
 types.mod
	$(CC) $(OPT1) cross_section_types.f90
#-------------------------------------------------------------------------------
distdot.o: distdot.f90
	$(CC) $(OPT1) distdot.f90
#-------------------------------------------------------------------------------
dump_inguess_module.o dump_inguess_module.mod: dump_inguess_module.f90 \
 types.mod parameter_types.mod filename_types.mod vector_types.mod \
 cross_section_types.mod geometry_types.mod angle_types.mod \
 multindex_types.mod global_variables.mod termination_module.mod
	$(CC) $(OPT1) dump_inguess_module.f90
#-------------------------------------------------------------------------------
execution_module.o execution_module.mod: execution_module.f90 types.mod \
 parameter_types.mod filename_types.mod vector_types.mod \
 cross_section_types.mod geometry_types.mod angle_types.mod \
 multindex_types.mod global_variables.mod solver_module.mod
	$(CC) $(OPT1) execution_module.f90
#-------------------------------------------------------------------------------
filename_types.o filename_types.mod: filename_types.f90
	$(CC) $(OPT1) filename_types.f90
#-------------------------------------------------------------------------------
general_utility_module.o general_utility_module.mod: \
 general_utility_module.f90 precision.mod types.mod vector_types.mod
	$(CC) $(OPT1) general_utility_module.f90
#-------------------------------------------------------------------------------
geometry_types.o geometry_types.mod: geometry_types.f90 types.mod \
 vector_types.mod
	$(CC) $(OPT1) geometry_types.f90
#-------------------------------------------------------------------------------
global_variables.o global_variables.mod: global_variables.f90 types.mod \
 parameter_types.mod filename_types.mod vector_types.mod \
 cross_section_types.mod geometry_types.mod angle_types.mod \
 multindex_types.mod
	$(CC) $(OPT1) global_variables.f90
#-------------------------------------------------------------------------------
inner_iteration_module.o inner_iteration_module.mod: \
 inner_iteration_module.f90  types.mod \
 parameter_types.mod filename_types.mod vector_types.mod \
 cross_section_types.mod geometry_types.mod angle_types.mod \
 multindex_types.mod global_variables.mod sweep_module.mod \
 termination_module.mod
	$(CC) $(OPT1) inner_iteration_module.f90
#-------------------------------------------------------------------------------
input_check.o check_input.mod: input_check.f90 types.mod \
 parameter_types.mod filename_types.mod vector_types.mod \
 cross_section_types.mod geometry_types.mod angle_types.mod \
 multindex_types.mod global_variables.mod
	$(CC) $(OPT1) input_check.f90
#-------------------------------------------------------------------------------

jfnk_module.o jfnk_module.mod: jfnk_module.f90 types.mod \
 parameter_types.mod filename_types.mod vector_types.mod \
 cross_section_types.mod geometry_types.mod angle_types.mod \
 multindex_types.mod global_variables.mod inner_iteration_module.mod \
 dump_inguess_module.mod
	$(CC) $(OPT1) jfnk_module.f90
#-------------------------------------------------------------------------------
main.o: main.f90  types.mod parameter_types.mod \
 filename_types.mod vector_types.mod cross_section_types.mod \
 geometry_types.mod angle_types.mod multindex_types.mod \
 global_variables.mod termination_module.mod setup_module.mod \
 execution_module.mod wrapup_module.mod  \
 global_variables.mod
	$(CC) $(OPT1) main.f90
#-------------------------------------------------------------------------------
multindex_types.o multindex_types.mod: multindex_types.f90 types.mod \
 vector_types.mod
	$(CC) $(OPT1) multindex_types.f90
#-------------------------------------------------------------------------------
outer_iteration_module.o outer_iteration_module.mod: \
 outer_iteration_module.f90 types.mod parameter_types.mod \
 filename_types.mod vector_types.mod cross_section_types.mod \
 geometry_types.mod angle_types.mod multindex_types.mod \
 global_variables.mod wrapup_module.mod termination_module.mod \
 inner_iteration_module.mod dump_inguess_module.mod
	$(CC) $(OPT1) outer_iteration_module.f90
#-------------------------------------------------------------------------------
parameter_types.o parameter_types.mod: parameter_types.f90 types.mod
	$(CC) $(OPT1) parameter_types.f90
#-------------------------------------------------------------------------------
precmod.o precision.mod: precmod.f90
	$(CC) $(OPT1) precmod.f90
#-------------------------------------------------------------------------------
quadrature_module.o quadrature_module.mod: quadrature_module.f90 \
 types.mod parameter_types.mod vector_types.mod angle_types.mod \
 global_variables.mod termination_module.mod
	$(CC) $(OPT1) quadrature_module.f90
#-------------------------------------------------------------------------------
read_cross_section_module.o read_cross_section_module.mod: \
 read_cross_section_module.f90 types.mod parameter_types.mod \
 filename_types.mod cross_section_types.mod global_variables.mod \
 termination_module.mod
	$(CC) $(OPT1) read_cross_section_module.f90
#-------------------------------------------------------------------------------
read_inflow_module.o read_inflow_module.mod: read_inflow_module.f90 \
 types.mod parameter_types.mod filename_types.mod multindex_types.mod \
 global_variables.mod termination_module.mod
	$(CC) $(OPT1) read_inflow_module.f90
#-------------------------------------------------------------------------------
read_module.o read_module.mod: read_module.f90  \
 strings.mod types.mod parameter_types.mod filename_types.mod \
 vector_types.mod cross_section_types.mod geometry_types.mod \
 angle_types.mod multindex_types.mod global_variables.mod \
 read_cross_section_module.mod readmesh_module.mod read_source_module.mod \
 read_inflow_module.mod quadrature_module.mod check_input.mod \
 termination_module.mod
	$(CC) $(OPT1) read_module.f90
#-------------------------------------------------------------------------------
read_source_module.o read_source_module.mod: read_source_module.f90 \
 types.mod parameter_types.mod filename_types.mod multindex_types.mod \
 global_variables.mod termination_module.mod
	$(CC) $(OPT1) read_source_module.f90
#-------------------------------------------------------------------------------
readmesh_module.o readmesh_module.mod: readmesh_module.f90 \
  types.mod parameter_types.mod filename_types.mod \
 vector_types.mod geometry_types.mod global_variables.mod \
 termination_module.mod
	$(CC) $(OPT1) readmesh_module.f90
#-------------------------------------------------------------------------------
setup_module.o setup_module.mod: setup_module.f90  \
 types.mod parameter_types.mod filename_types.mod vector_types.mod \
 cross_section_types.mod geometry_types.mod angle_types.mod \
 multindex_types.mod global_variables.mod read_module.mod \
 quadrature_module.mod sweep_module.mod termination_module.mod
	$(CC) $(OPT1) setup_module.f90
#-------------------------------------------------------------------------------
solver_module.o solver_module.mod: solver_module.f90 types.mod \
 parameter_types.mod filename_types.mod vector_types.mod \
 cross_section_types.mod geometry_types.mod angle_types.mod \
 multindex_types.mod global_variables.mod ahotc_matrix_module.mod \
 sph_harmonics_module.mod outer_iteration_module.mod jfnk_module.mod \
 termination_module.mod
	$(CC) $(OPT1) solver_module.f90
#-------------------------------------------------------------------------------
sph_harmonics_module.o sph_harmonics_module.mod: sph_harmonics_module.f90 \
 types.mod parameter_types.mod vector_types.mod angle_types.mod \
 multindex_types.mod termination_module.mod
	$(CC) $(OPT1) sph_harmonics_module.f90
#-------------------------------------------------------------------------------
stringmod.o strings.mod: stringmod.f90 precision.mod
	$(CC) $(OPT1) stringmod.f90
#-------------------------------------------------------------------------------
sweep_module.o sweep_module.mod: sweep_module.f90  \
 types.mod parameter_types.mod filename_types.mod vector_types.mod \
 cross_section_types.mod geometry_types.mod angle_types.mod \
 multindex_types.mod global_variables.mod cell_splitting_module.mod
	$(CC) $(OPT1) sweep_module.f90
#-------------------------------------------------------------------------------
termination_module.o termination_module.mod: termination_module.f90 \
 global_variables.mod
	$(CC) $(OPT1) termination_module.f90
#-------------------------------------------------------------------------------
transport_kernel_module_CCE.o transport_kernel_module_cce.mod: \
 transport_kernel_module_CCE.f90 types.mod parameter_types.mod \
 filename_types.mod vector_types.mod cross_section_types.mod \
 geometry_types.mod angle_types.mod multindex_types.mod \
 global_variables.mod termination_module.mod general_utility_module.mod \
 types.mod types.mod types.mod types.mod
	$(CC) $(OPT1) transport_kernel_module_CCE.f90
#-------------------------------------------------------------------------------
transport_kernel_module_LC.o transport_kernel_module_lc.mod: \
    transport_kernel_module_LC.f90 types.mod parameter_types.mod \
    filename_types.mod vector_types.mod cross_section_types.mod \
    geometry_types.mod angle_types.mod multindex_types.mod \
    global_variables.mod termination_module.mod general_utility_module.mod
	$(CC) $(OPT1) $(OPT2) transport_kernel_module_LC.f90

#-------------------------------------------------------------------------------

transport_kernel_module_SC.o transport_kernel_module_sc.mod: \
    transport_kernel_module_SC.f90 types.mod parameter_types.mod \
    filename_types.mod vector_types.mod cross_section_types.mod \
    geometry_types.mod angle_types.mod multindex_types.mod \
    global_variables.mod termination_module.mod general_utility_module.mod
	$(CC) $(OPT1) $(OPT2) transport_kernel_module_SC.f90
#-------------------------------------------------------------------------------
types.o types.mod: types.f90
	$(CC) $(OPT1) types.f90
#-------------------------------------------------------------------------------
vector_types.o vector_types.mod: vector_types.f90 types.mod
	$(CC) $(OPT1) vector_types.f90
#-------------------------------------------------------------------------------
wrapup_module.o wrapup_module.mod: wrapup_module.f90 types.mod \
 parameter_types.mod filename_types.mod vector_types.mod \
 cross_section_types.mod geometry_types.mod angle_types.mod \
 global_variables.mod termination_module.mod
	$(CC) $(OPT1) wrapup_module.f90
#-------------------------------------------------------------------------------
