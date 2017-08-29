#-----------------------------------------------------------------------------------
#Compile Options
#-----------------------------------------------------------------------------------
FC := mpifort
EXE := thor-1.0
EXT := .exe

all:      TYPE :=
debug:    TYPE := _Debug
coverage: TYPE := _Coverage

all:      OPT := -O3
debug:    OPT := -O0 -g -fbacktrace -ffpe-trap=zero,overflow,underflow
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
        adv_read_module.f90								\
        ahotc_matrix_module.f90						\
        angle_types.f90										\
        cell_splitting_module.f90					\
        cross_section_types.f90						\
        dump_inguess_module.f90						\
        error_generator.f90								\
        execution_module.f90							\
        filename_types.f90								\
        general_utility_module.f90				\
        geometry_types.f90								\
        global_variables.f90							\
        inner_iteration_module.f90				\
        input_check.f90										\
        jfnk_module.f90										\
        main.f90													\
        multindex_types.f90								\
        outer_iteration_module.f90				\
        parameter_types.f90								\
        precision.f90												\
        quadrature_module.f90							\
        read_cross_section_module.f90			\
        read_inflow_module.f90						\
        read_module.f90										\
        read_source_module.f90						\
        readmesh_module.f90								\
        setup_module.f90									\
        solver_module.f90									\
        sph_harmonics_module.f90					\
        sweep_module.f90									\
        termination_module.f90						\
        transport_kernel_module_cce.f90		\
        transport_kernel_module_lc.f90		\
        transport_kernel_module_sc.f90		\
        types.f90													\
        vector_types.f90									\
        wrapup_module.f90									\


OBJ := $(SRC:.f90=.o)
MOD := $(OBJ:.o=.mod)

FILTER := \
          transport_kernel_module_lc.f90		\
          transport_kernel_module_sc.f90		\

OBJ_FILTER := $(FILTER:.f90=.o)
MOD_FILTER := $(FILTER:.f90=.mod)

#-----------------------------------------------------------------------------------
#Add simple compile contrib files
#-----------------------------------------------------------------------------------
CONT_SRC := \
            $(CONTRIB)/distdot/distdot.f90      \
            $(CONTRIB)/stringmod/stringmod.f90  \

CONT_OBJ := $(patsubst %.f90,%.o,$(notdir $(CONT_SRC)))

CONT_MOD := $(CONT_OBJ:.o=.mod)

#-----------------------------------------------------------------------------------
#Add lib files
#-----------------------------------------------------------------------------------
LIB_SRC :=  \

LIB_OBJ := $(patsubst %.f90,%.o,$(notdir $(LIB_SRC)))

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
.PHONY: all debug coverage clean reset

all: $(EXE)

debug: $(EXE)

coverage: $(EXE)

#Intended to clean up compilation artifacts but leave executable & coverage
clean:
	rm -f $(OBJ) $(CONT_OBJ) $(LIB_OBJ)
	rm -f $(MOD) $(CONT_MOD) $(LIB_MOD)
	rm -f $(COMP_DEP)

#Intended to reset directory to fresh state with no exe or artifacts
reset: clean
	rm -f *.gcno *.gcda
	rm -f -r $(EXE)*.dSYM
	rm -f ../$(EXE)*

#-----------------------------------------------------------------------------------
#Generics for source files
#-----------------------------------------------------------------------------------
$(filter-out $(OBJ_FILTER), $(OBJ)): %.o:	%.f90
	$(FC) -c $(OPT) $<

$(filter-out $(MOD_FILTER), $(MOD)):	%.mod:	%.f90
	$(FC) -c $(OPT) $<

$(EXE): $(OBJ) $(CONT_OBJ) $(LIB_OBJ) $(COMP_DEP)
	$(FC) -o $@$(TYPE)$(EXT) $(OPT) $(OBJ) $(CONT_OBJ) $(LIB_OBJ) $(COMP_DEP)
	mv ./$(EXE)$(TYPE)$(EXT) ../

#-----------------------------------------------------------------------------------
#Generics for contrib files
#-----------------------------------------------------------------------------------
$(CONT_OBJ): %.o: $(filter %.f90, $(CONT_SRC))
	$(FC) -c $(OPT) $^

$(CONT_MOD):	%.mod:	$(filter %.f90, $(CONT_SRC))
	$(FC) -c $(OPT) $^

#-----------------------------------------------------------------------------------
#Generics for lib files
#-----------------------------------------------------------------------------------
$(LIB_OBJ):	%.o: $(filter %.f90, $(LIB_SRC))
	$(FC) -c $(OPT) $^

$(LIB_MOD): %.mod:	$(filter %.f90, $(LIB_SRC))
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
#Use [gfortran -M -cpp *.f90] repeatedly until clean compile to update rules below
#-----------------------------------------------------------------------------------
adv_read_module.o adv_read_module.mod: adv_read_module.f90 \
  error_generator.mod global_variables.mod parameter_types.mod

ahotc_matrix_module.o ahotc_matrix_module.mod: ahotc_matrix_module.f90 \
  types.mod parameter_types.mod multindex_types.mod

angle_types.o angle_types.mod: angle_types.f90 types.mod vector_types.mod

cell_splitting_module.o cell_splitting_module.mod: \
  cell_splitting_module.f90 types.mod parameter_types.mod \
  filename_types.mod vector_types.mod cross_section_types.mod \
  geometry_types.mod angle_types.mod multindex_types.mod \
  global_variables.mod general_utility_module.mod \
  transport_kernel_module_sc.mod transport_kernel_module_lc.mod \
  transport_kernel_module_cce.mod termination_module.mod

cross_section_types.o cross_section_types.mod: cross_section_types.f90 \
  types.mod

dump_inguess_module.o dump_inguess_module.mod: dump_inguess_module.f90 \
  types.mod parameter_types.mod filename_types.mod vector_types.mod \
  cross_section_types.mod geometry_types.mod angle_types.mod \
  multindex_types.mod global_variables.mod termination_module.mod

error_generator.o error_generator.mod: error_generator.f90

execution_module.o execution_module.mod: execution_module.f90 types.mod \
  parameter_types.mod filename_types.mod vector_types.mod \
  cross_section_types.mod geometry_types.mod angle_types.mod \
  multindex_types.mod global_variables.mod solver_module.mod

filename_types.o filename_types.mod: filename_types.f90

general_utility_module.o general_utility_module.mod: \
  general_utility_module.f90 precision.mod types.mod vector_types.mod

geometry_types.o geometry_types.mod: geometry_types.f90 types.mod \
  vector_types.mod

global_variables.o global_variables.mod: global_variables.f90 types.mod \
  parameter_types.mod filename_types.mod vector_types.mod \
  cross_section_types.mod geometry_types.mod angle_types.mod \
  multindex_types.mod

inner_iteration_module.o inner_iteration_module.mod: \
  inner_iteration_module.f90  types.mod \
  parameter_types.mod filename_types.mod vector_types.mod \
  cross_section_types.mod geometry_types.mod angle_types.mod \
  multindex_types.mod global_variables.mod sweep_module.mod \
  termination_module.mod

input_check.o check_input.mod: input_check.f90 types.mod \
  parameter_types.mod filename_types.mod vector_types.mod \
  cross_section_types.mod geometry_types.mod angle_types.mod \
  multindex_types.mod global_variables.mod

jfnk_module.o jfnk_module.mod: jfnk_module.f90 types.mod \
  parameter_types.mod filename_types.mod vector_types.mod \
  cross_section_types.mod geometry_types.mod angle_types.mod \
  multindex_types.mod global_variables.mod inner_iteration_module.mod \
  dump_inguess_module.mod

main.o: main.f90  types.mod parameter_types.mod \
  filename_types.mod vector_types.mod cross_section_types.mod \
  geometry_types.mod angle_types.mod multindex_types.mod \
  global_variables.mod termination_module.mod setup_module.mod \
  execution_module.mod wrapup_module.mod  \
  global_variables.mod

multindex_types.o multindex_types.mod: multindex_types.f90 types.mod \
  vector_types.mod

outer_iteration_module.o outer_iteration_module.mod: \
  outer_iteration_module.f90 types.mod parameter_types.mod \
  filename_types.mod vector_types.mod cross_section_types.mod \
  geometry_types.mod angle_types.mod multindex_types.mod \
  global_variables.mod wrapup_module.mod termination_module.mod \
  inner_iteration_module.mod dump_inguess_module.mod

parameter_types.o parameter_types.mod: parameter_types.f90 types.mod

precmod.o precision.mod: precision.f90

quadrature_module.o quadrature_module.mod: quadrature_module.f90 \
  types.mod parameter_types.mod vector_types.mod angle_types.mod \
  global_variables.mod termination_module.mod

read_cross_section_module.o read_cross_section_module.mod: \
  read_cross_section_module.f90 types.mod parameter_types.mod \
  filename_types.mod cross_section_types.mod global_variables.mod \
  termination_module.mod

read_inflow_module.o read_inflow_module.mod: read_inflow_module.f90 \
  types.mod parameter_types.mod filename_types.mod multindex_types.mod \
  global_variables.mod termination_module.mod

read_module.o read_module.mod: read_module.f90  \
  stringmod.mod types.mod parameter_types.mod filename_types.mod \
  vector_types.mod cross_section_types.mod geometry_types.mod \
  angle_types.mod multindex_types.mod global_variables.mod \
  read_cross_section_module.mod readmesh_module.mod read_source_module.mod \
  read_inflow_module.mod quadrature_module.mod check_input.mod \
  termination_module.mod adv_read_module.mod

read_source_module.o read_source_module.mod: read_source_module.f90 \
  types.mod parameter_types.mod filename_types.mod multindex_types.mod \
  global_variables.mod termination_module.mod

readmesh_module.o readmesh_module.mod: readmesh_module.f90 \
  types.mod parameter_types.mod filename_types.mod \
  vector_types.mod geometry_types.mod global_variables.mod \
  termination_module.mod

setup_module.o setup_module.mod: setup_module.f90  \
  types.mod parameter_types.mod filename_types.mod vector_types.mod \
  cross_section_types.mod geometry_types.mod angle_types.mod \
  multindex_types.mod global_variables.mod read_module.mod \
  quadrature_module.mod sweep_module.mod termination_module.mod

solver_module.o solver_module.mod: solver_module.f90 types.mod \
  parameter_types.mod filename_types.mod vector_types.mod \
  cross_section_types.mod geometry_types.mod angle_types.mod \
  multindex_types.mod global_variables.mod ahotc_matrix_module.mod \
  sph_harmonics_module.mod outer_iteration_module.mod jfnk_module.mod \
  termination_module.mod

sph_harmonics_module.o sph_harmonics_module.mod: sph_harmonics_module.f90 \
  types.mod parameter_types.mod vector_types.mod angle_types.mod \
  multindex_types.mod termination_module.mod

sweep_module.o sweep_module.mod: sweep_module.f90  \
  types.mod parameter_types.mod filename_types.mod vector_types.mod \
  cross_section_types.mod geometry_types.mod angle_types.mod \
  multindex_types.mod global_variables.mod cell_splitting_module.mod

termination_module.o termination_module.mod: termination_module.f90 \
  global_variables.mod

transport_kernel_module_cce.o transport_kernel_module_cce.mod: \
  transport_kernel_module_cce.f90 types.mod parameter_types.mod \
  filename_types.mod vector_types.mod cross_section_types.mod \
  geometry_types.mod angle_types.mod multindex_types.mod \
  global_variables.mod termination_module.mod general_utility_module.mod \
  types.mod types.mod types.mod types.mod

transport_kernel_module_lc.o transport_kernel_module_lc.mod: \
  transport_kernel_module_lc.f90 types.mod parameter_types.mod \
  filename_types.mod vector_types.mod cross_section_types.mod \
  geometry_types.mod angle_types.mod multindex_types.mod \
  global_variables.mod termination_module.mod general_utility_module.mod
		$(FC) -c $(OPT) -fdefault-real-8 $<

transport_kernel_module_sc.o transport_kernel_module_sc.mod: \
  transport_kernel_module_sc.f90 types.mod parameter_types.mod \
  filename_types.mod vector_types.mod cross_section_types.mod \
  geometry_types.mod angle_types.mod multindex_types.mod \
  global_variables.mod termination_module.mod general_utility_module.mod
		$(FC) -c $(OPT) -fdefault-real-8 $<

types.o types.mod: types.f90

vector_types.o vector_types.mod: vector_types.f90 types.mod

wrapup_module.o wrapup_module.mod: wrapup_module.f90 types.mod \
  parameter_types.mod filename_types.mod vector_types.mod \
  cross_section_types.mod geometry_types.mod angle_types.mod \
  global_variables.mod termination_module.mod
