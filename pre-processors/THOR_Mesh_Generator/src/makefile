#-----------------------------------------------------------------------------------
#Compile Options
#-----------------------------------------------------------------------------------
FC := gfortran
EXE := Thor_Mesh_Generator
EXT := .exe

all:			TYPE :=
MPparallel:    TYPE := _MP
debug:    TYPE := _Debug
coverage: TYPE := _Coverage

all:      OPT := -O3
MPparallel: OPT := -O3 -fopenmp
debug:    OPT := -O0 -g
coverage: OPT := -O0 -g --coverage

#-----------------------------------------------------------------------------------
#Directory Paths
#-----------------------------------------------------------------------------------
CONTRIB := ../../../contrib
LIB     := ../../../lib

#-----------------------------------------------------------------------------------
#Add source files as necessary
#If the files are not comiled using the genieric rules, add commands under their
#rule declaration. Add these items to FILTER
#-----------------------------------------------------------------------------------
SRC := \
				mesh_converter.f95 	      	\
				globals.f95        					\
				thor_mesh.f95      					\
				unv.f95      					\
				boundary_conditions.f95      					\
				gmesh.f95

OBJ := $(SRC:.f95=.o)
MOD := $(OBJ:.o=.mod)

FILTER := \

OBJ_FILTER := $(FILTER:.f95=.o)
MOD_FILTER := $(FILTER:.f95=.mod)

#-----------------------------------------------------------------------------------
#Add simple compile contrib files
#-----------------------------------------------------------------------------------
CONT_SRC :=\
						$(CONTRIB)/tet_mesh_tet_neighbors/tet_mesh_tet_neighbors.f95

CONT_OBJ := $(patsubst %.f95,%.o,$(notdir $(CONT_SRC)))

CONT_MOD := $(CONT_OBJ:.o=.mod)

#-----------------------------------------------------------------------------------
#Add lib files
#-----------------------------------------------------------------------------------
LIB_SRC := \
						$(LIB)/integer_array_tools/integer_array_tools.f95

LIB_OBJ := $(patsubst %.f95,%.o,$(notdir $(LIB_SRC)))

LIB_MOD := $(LIB_OBJ:.o=.mod)

#-----------------------------------------------------------------------------------
#Complex external dependencies. Each of these requires a custom build rule.
#Intended for depenencies with thier own makefiles or other build system
#-----------------------------------------------------------------------------------
COMP_DEP := \
#Phony targets for cleaning and building

#-----------------------------------------------------------------------------------
#Phony targets for cleaning and building
#-----------------------------------------------------------------------------------
.PHONY: all clean reset

all: $(EXE)

MPparallel: $(EXE)

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
boundary_conditions.o boundary_conditions.mod: boundary_conditions.f95

globals.o globals.mod: globals.f95 integer_array_tools.mod \
 tet_mesh_tet_neighbors.mod

gmesh.o gmesh.mod: gmesh.f95 globals.mod

mesh_converter.o: mesh_converter.f95 globals.mod gmesh.mod unv.mod \
 thor_mesh.mod integer_array_tools.mod

thor_mesh.o thor_mesh.mod: thor_mesh.f95 globals.mod \
 tet_mesh_tet_neighbors.mod

unv.o unv.mod: unv.f95 globals.mod boundary_conditions.mod
