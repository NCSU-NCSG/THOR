#-----------------------------------------------------------------------------------
#Compile Options
#-----------------------------------------------------------------------------------
FC := mpifort
EXE := thor-1.0
EXT := .exe

all:      TYPE :=
debug:    TYPE := _Debug
catch:    TYPE := _Catch
coverage: TYPE := _Coverage

all:      OPT := -O3
debug:    OPT := -O0 -g -fbacktrace -ffpe-trap=zero,overflow,underflow
catch:    OPT := -O0 -g -Wall -W -Wsurprising -Werror
coverage: OPT := -O0 -g --coverage

#-----------------------------------------------------------------------------------
#Directory Paths
#-----------------------------------------------------------------------------------
CONTRIB := ../../contrib
LIB     := ../../lib

#-----------------------------------------------------------------------------------
#Add source files as necessary
#If the files are not compiled using the generic rules, add commands under their
#rule declaration. Add these items to FILTER
#-----------------------------------------------------------------------------------
SRC :=  \
        ahotc_matrix_module.f95						\
        angle_types.f95										\
        cell_splitting_module.f95					\
        cross_section_types.f95						\
        source_types.f95						\
        dump_inguess_module.f95						\
        error_generator.f95								\
        execution_module.f95							\
        filename_types.f95								\
        general_utility_module.f95				\
        geometry_types.f95								\
        global_variables.f95							\
        inner_iteration_module.f95				\
        input_check.f95										\
        jfnk_module.f95										\
        main.f95													\
        multindex_types.f95								\
        outer_iteration_module.f95				\
        parameter_types.f95								\
        precision.f95											\
        quadrature_module.f95							\
        read_cross_section_module.f95			\
        read_inflow_module.f95						\
        read_inp_module_legacy.f95				\
        read_inp_module.f95 	        		\
        read_module.f95										\
        read_source_module.f95						\
        read_mesh_module.f95								\
        setup_module.f95									\
        solver_module.f95									\
        sph_harmonics_module.f95					\
        sweep_module.f95									\
        termination_module.f95						\
        transport_kernel_module_cce.f95		\
        transport_kernel_module_lc.f95		\
        transport_kernel_module_sc.f95		\
        types.f95													\
        vector_types.f95									\
        wrapup_module.f95									\


OBJ := $(SRC:.f95=.o)
MOD := $(OBJ:.o=.mod)

FILTER := \
          transport_kernel_module_lc.f95		\
          transport_kernel_module_sc.f95		\

OBJ_FILTER := $(FILTER:.f95=.o)
MOD_FILTER := $(FILTER:.f95=.mod)

#-----------------------------------------------------------------------------------
#Add simple compile contrib files
#-----------------------------------------------------------------------------------
CONT_SRC := \
            $(CONTRIB)/distdot/distdot.f95      \
            $(CONTRIB)/stringmod/stringmod.f95  \

CONT_OBJ := $(patsubst %.f95,%.o,$(notdir $(CONT_SRC)))

CONT_MOD := $(CONT_OBJ:.o=.mod)

#-----------------------------------------------------------------------------------
#Add lib files
#-----------------------------------------------------------------------------------
LIB_SRC :=  \
            $(LIB)/integer_array_tools/integer_array_tools.f95 \

LIB_OBJ := $(patsubst %.f95,%.o,$(notdir $(LIB_SRC)))

LIB_MOD := $(LIB_OBJ:.o=.mod)

#-----------------------------------------------------------------------------------
#Complex external dependencies. Each of these requires a custom build rule.
#Intended for dependencies with their own makefiles or other build system
#-----------------------------------------------------------------------------------
COMP_DEP := \
            libskit.a \
            liblapack.a\
            librefblas.a\
            libtmglib.a

#-----------------------------------------------------------------------------------
#Phony targets for cleaning and building
#-----------------------------------------------------------------------------------
.PHONY: all debug catch coverage clean reset

print-%  : ; @echo $* = $($*)

all: $(EXE)

debug: $(EXE)

catch: $(EXE)

coverage: $(EXE)

#Intended to clean up compilation artifacts but leave executable & coverage
clean:
	rm -f $(OBJ) $(CONT_OBJ) $(LIB_OBJ)
	rm -f $(MOD) $(CONT_MOD) $(LIB_MOD)
	rm -f $(COMP_DEP)

#Intended to reset directory to fresh state with no exe or artifacts
reset: clean
	rm -f *.gcno *.gcda *.o *.mod
	rm -f -r $(EXE)*.dSYM
	rm -f ../$(EXE)*

#-----------------------------------------------------------------------------------
#Generics for source files
#-----------------------------------------------------------------------------------
$(filter-out $(OBJ_FILTER), $(OBJ)): %.o:	%.f95
	$(FC) -c $(OPT) $<

$(filter-out $(MOD_FILTER), $(MOD)):	%.mod:	%.f95
	$(FC) -c $(OPT) $<

$(EXE): $(OBJ) $(CONT_OBJ) $(LIB_OBJ) $(COMP_DEP)
	$(FC) -o $@$(TYPE)$(EXT) $(OPT) $(OBJ) $(CONT_OBJ) $(LIB_OBJ) $(COMP_DEP)
	mv ./$(EXE)$(TYPE)$(EXT) ../

#-----------------------------------------------------------------------------------
#Generics for contrib files
#-----------------------------------------------------------------------------------
$(CONT_OBJ): %.o: $(filter %.f95, $(CONT_SRC))
	$(FC) -c $(OPT) $^

$(CONT_MOD):	%.mod:	$(filter %.f95, $(CONT_SRC))
	$(FC) -c $(OPT) $^

#-----------------------------------------------------------------------------------
#Generics for lib files
#-----------------------------------------------------------------------------------
$(LIB_OBJ):	%.o: $(filter %.f95, $(LIB_SRC))
	$(FC) -c $(OPT) $^

$(LIB_MOD): %.mod:	$(filter %.f95, $(LIB_SRC))
	$(FC) -c $(OPT) $^

#-----------------------------------------------------------------------------------
#Rules for entries in COMP_DEP. Don't forget to add them to make clean / reset
#-----------------------------------------------------------------------------------
libskit.a:
	cd $(CONTRIB)/SPARSKIT2; make
	cp $(CONTRIB)/SPARSKIT2/libskit.a .

liblapack.a:
	cp $(CONTRIB)/libs/liblapack.a .

librefblas.a:
	cp $(CONTRIB)/libs/librefblas.a .

libtmglib.a:
	cp $(CONTRIB)/libs/libtmglib.a .

#-----------------------------------------------------------------------------------
#Dependency List
#Use [gfortran -M -cpp *.f95] repeatedly until clean compile to update rules below
#-----------------------------------------------------------------------------------
ahotc_matrix_module.o ahotc_matrix_module.mod: ahotc_matrix_module.f95 \
  types.mod parameter_types.mod multindex_types.mod general_utility_module.mod

angle_types.o angle_types.mod: angle_types.f95 types.mod vector_types.mod

cell_splitting_module.o cell_splitting_module.mod: \
  cell_splitting_module.f95 types.mod parameter_types.mod \
  filename_types.mod vector_types.mod cross_section_types.mod \
  geometry_types.mod angle_types.mod multindex_types.mod \
  global_variables.mod general_utility_module.mod \
  transport_kernel_module_sc.mod transport_kernel_module_lc.mod \
  transport_kernel_module_cce.mod termination_module.mod stringmod.mod

cross_section_types.o cross_section_types.mod: cross_section_types.f95 \
  types.mod

dump_inguess_module.o dump_inguess_module.mod: dump_inguess_module.f95 \
  types.mod parameter_types.mod filename_types.mod vector_types.mod \
  cross_section_types.mod geometry_types.mod angle_types.mod \
  multindex_types.mod global_variables.mod termination_module.mod

error_generator.o error_generator.mod: error_generator.f95

execution_module.o execution_module.mod: execution_module.f95 types.mod \
  parameter_types.mod filename_types.mod vector_types.mod \
  cross_section_types.mod geometry_types.mod angle_types.mod \
  multindex_types.mod global_variables.mod solver_module.mod

filename_types.o filename_types.mod: filename_types.f95

general_utility_module.o general_utility_module.mod: \
  general_utility_module.f95 precision.mod types.mod vector_types.mod

geometry_types.o geometry_types.mod: geometry_types.f95 types.mod \
  vector_types.mod

global_variables.o global_variables.mod: global_variables.f95 types.mod \
  parameter_types.mod filename_types.mod vector_types.mod \
  cross_section_types.mod source_types.mod geometry_types.mod angle_types.mod \
  multindex_types.mod integer_array_tools.mod

inner_iteration_module.o inner_iteration_module.mod: \
  inner_iteration_module.f95  types.mod \
  parameter_types.mod filename_types.mod vector_types.mod \
  cross_section_types.mod geometry_types.mod angle_types.mod \
  multindex_types.mod global_variables.mod sweep_module.mod \
  termination_module.mod

input_check.o check_input.mod: input_check.f95 types.mod \
  parameter_types.mod filename_types.mod vector_types.mod \
  cross_section_types.mod geometry_types.mod angle_types.mod \
  multindex_types.mod global_variables.mod

jfnk_module.o jfnk_module.mod: jfnk_module.f95 types.mod \
  parameter_types.mod filename_types.mod vector_types.mod \
  cross_section_types.mod geometry_types.mod angle_types.mod \
  multindex_types.mod global_variables.mod inner_iteration_module.mod \
  dump_inguess_module.mod outer_iteration_module.mod

main.o: main.f95  types.mod parameter_types.mod \
  filename_types.mod vector_types.mod cross_section_types.mod \
  geometry_types.mod angle_types.mod multindex_types.mod \
  global_variables.mod termination_module.mod setup_module.mod \
  execution_module.mod wrapup_module.mod  \
  global_variables.mod

multindex_types.o multindex_types.mod: multindex_types.f95 types.mod \
  vector_types.mod

outer_iteration_module.o outer_iteration_module.mod: \
  outer_iteration_module.f95 types.mod parameter_types.mod \
  filename_types.mod vector_types.mod cross_section_types.mod \
  geometry_types.mod angle_types.mod multindex_types.mod \
  global_variables.mod wrapup_module.mod termination_module.mod \
  inner_iteration_module.mod dump_inguess_module.mod

parameter_types.o parameter_types.mod: parameter_types.f95 types.mod

precision.o precision.mod: precision.f95

quadrature_module.o quadrature_module.mod: quadrature_module.f95 \
  types.mod parameter_types.mod vector_types.mod angle_types.mod \
  global_variables.mod termination_module.mod

read_cross_section_module.o read_cross_section_module.mod: \
  read_cross_section_module.f95 types.mod parameter_types.mod \
  filename_types.mod cross_section_types.mod global_variables.mod \
  termination_module.mod stringmod.mod

read_inflow_module.o read_inflow_module.mod: read_inflow_module.f95 \
  types.mod parameter_types.mod filename_types.mod multindex_types.mod \
  global_variables.mod termination_module.mod

read_inp_module.o read_inp_module.mod: read_inp_module.f95  \
  stringmod.mod global_variables.mod

read_inp_module_legacy.o read_inp_module_legacy.mod: read_inp_module_legacy.f95 \
  stringmod.mod error_generator.mod global_variables.mod parameter_types.mod

read_mesh_module.o read_mesh_module.mod: read_mesh_module.f95 \
  types.mod parameter_types.mod filename_types.mod \
  vector_types.mod geometry_types.mod global_variables.mod \
  termination_module.mod

read_module.o read_module.mod: read_module.f95  \
  stringmod.mod types.mod parameter_types.mod filename_types.mod \
  vector_types.mod cross_section_types.mod geometry_types.mod \
  angle_types.mod multindex_types.mod global_variables.mod \
  read_cross_section_module.mod read_mesh_module.mod read_source_module.mod \
  read_inflow_module.mod quadrature_module.mod check_input.mod \
  termination_module.mod read_inp_module_legacy.mod read_inp_module.mod

read_source_module.o read_source_module.mod: read_source_module.f95 \
  types.mod parameter_types.mod filename_types.mod multindex_types.mod \
  global_variables.mod termination_module.mod stringmod.mod

setup_module.o setup_module.mod: setup_module.f95  \
  types.mod parameter_types.mod filename_types.mod vector_types.mod \
  cross_section_types.mod geometry_types.mod angle_types.mod \
  multindex_types.mod global_variables.mod read_module.mod \
  quadrature_module.mod sweep_module.mod termination_module.mod

source_types.o source_types.mod: source_types.f95 \
  types.mod

solver_module.o solver_module.mod: solver_module.f95 types.mod \
  parameter_types.mod filename_types.mod vector_types.mod \
  cross_section_types.mod geometry_types.mod angle_types.mod \
  multindex_types.mod global_variables.mod ahotc_matrix_module.mod \
  sph_harmonics_module.mod outer_iteration_module.mod jfnk_module.mod \
  termination_module.mod

sph_harmonics_module.o sph_harmonics_module.mod: sph_harmonics_module.f95 \
  types.mod parameter_types.mod vector_types.mod angle_types.mod \
  multindex_types.mod termination_module.mod global_variables.mod

sweep_module.o sweep_module.mod: sweep_module.f95  \
  types.mod parameter_types.mod filename_types.mod vector_types.mod \
  cross_section_types.mod geometry_types.mod angle_types.mod \
  multindex_types.mod global_variables.mod cell_splitting_module.mod

termination_module.o termination_module.mod: termination_module.f95 \
  global_variables.mod

transport_kernel_module_cce.o transport_kernel_module_cce.mod: \
  transport_kernel_module_cce.f95 types.mod parameter_types.mod \
  filename_types.mod vector_types.mod cross_section_types.mod \
  geometry_types.mod angle_types.mod multindex_types.mod \
  global_variables.mod termination_module.mod general_utility_module.mod \
  types.mod types.mod types.mod types.mod

transport_kernel_module_lc.o transport_kernel_module_lc.mod: \
  transport_kernel_module_lc.f95 types.mod parameter_types.mod \
  filename_types.mod vector_types.mod cross_section_types.mod \
  geometry_types.mod angle_types.mod multindex_types.mod \
  global_variables.mod termination_module.mod general_utility_module.mod
		$(FC) -c $(OPT) -fdefault-real-8 $<

transport_kernel_module_sc.o transport_kernel_module_sc.mod: \
  transport_kernel_module_sc.f95 types.mod parameter_types.mod \
  filename_types.mod vector_types.mod cross_section_types.mod \
  geometry_types.mod angle_types.mod multindex_types.mod \
  global_variables.mod termination_module.mod general_utility_module.mod
		$(FC) -c $(OPT) -fdefault-real-8 $<

types.o types.mod: types.f95

vector_types.o vector_types.mod: vector_types.f95 types.mod

wrapup_module.o wrapup_module.mod: wrapup_module.f95 types.mod \
  parameter_types.mod filename_types.mod vector_types.mod \
  cross_section_types.mod geometry_types.mod angle_types.mod \
  global_variables.mod termination_module.mod general_utility_module.mod
