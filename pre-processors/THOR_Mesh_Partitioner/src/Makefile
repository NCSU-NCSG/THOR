#-----------------------------------------------------------------------------------
#Compile Options
#-----------------------------------------------------------------------------------
FC := gfortran
EXE := Thor_Mesh_Partitioner
EXT := .exe

all:			TYPE :=
debug:    TYPE := _Debug
coverage: TYPE := _Coverage

all:      OPT := -O3 -fopenmp
debug:    OPT := -O0 -g -fbacktrace -ffpe-trap=zero,overflow,underflow
coverage: OPT := -O0 -g --coverage

#-----------------------------------------------------------------------------------
#Directory Paths
#-----------------------------------------------------------------------------------
CONTRIB :=
LIB     :=

#-----------------------------------------------------------------------------------
#Add source files as necessary
#If the files are not comiled using the genieric rules, add commands under their
#rule declaration. Add these items to FILTER
#-----------------------------------------------------------------------------------
SRC := global_variables.f95					\
			 mesh_splitter.f95						\
			 metis_module.f95							\
			 output_module.f95						\
			 reconstruction_module.f95		\
			 statistics_module.f95 				\
			 vtk_module.f95								\


OBJ := $(SRC:.f95=.o)
MOD := $(OBJ:.o=.mod)

FILTER := \

OBJ_FILTER := $(FILTER:.f95=.o)
MOD_FILTER := $(FILTER:.f95=.mod)

#-----------------------------------------------------------------------------------
#Add simple compile contrib files
#-----------------------------------------------------------------------------------
CONT_SRC :=

CONT_OBJ := $(patsubst %.f95,%.o,$(notdir $(CONT_SRC)))

CONT_MOD := $(CONT_OBJ:.o=.mod)

#-----------------------------------------------------------------------------------
#Add lib files
#-----------------------------------------------------------------------------------
LIB_SRC :=

LIB_OBJ := $(patsubst %.f95,%.o,$(notdir $(LIB_SRC)))

LIB_MOD := $(LIB_OBJ:.o=.mod)

#-----------------------------------------------------------------------------------
#Complex external dependencies. Each of these requires a custom build rule.
#Intended for depenencies with thier own makefiles or other build system
#-----------------------------------------------------------------------------------
COMP_DEP :=

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
$(filter-out $(OBJ_FILTER), $(OBJ)): %.o:	%.f95
	$(FC) -c $(OPT) $<

$(filter-out $(MOD_FILTER), $(MOD)):	%.mod:	%.f95
	$(FC) -c $(OPT) $<

$(EXE): $(OBJ) $(CONT_OBJ) $(LIB_OBJ)
	$(FC) -o $@$(TYPE)$(EXT) $(OPT) $(OBJ) $(CONT_OBJ) $(LIB_OBJ)
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

#-----------------------------------------------------------------------------------
#Dependency List
#Use [gfortran -M -cpp *.f95] repeatedly until clean compile to update rules below
#-----------------------------------------------------------------------------------
global_variables.o global_variables.mod: global_variables.f95

mesh_splitter.o: mesh_splitter.f95 global_variables.mod metis_module.mod \
 vtk_module.mod reconstruction_module.mod statistics_module.mod \
 output_module.mod global_variables.mod global_variables.mod \
 global_variables.mod global_variables.mod

metis_module.o metis_module.mod: metis_module.f95 global_variables.mod

output_module.o output_module.mod: output_module.f95 global_variables.mod

reconstruction_module.o reconstruction_module.mod: \
 reconstruction_module.f95 global_variables.mod

statistics_module.o statistics_module.mod: statistics_module.f95 \
 global_variables.mod

vtk_module.o vtk_module.mod: vtk_module.f95
