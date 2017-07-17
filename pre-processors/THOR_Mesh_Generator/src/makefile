
#Compile Options
FC := gfortran
					EXE := Thor_Mesh_Generator.exe
debug:    EXE := Thor_Mesh_Generator_Debug.exe
coverage: EXE := Thor_Mesh_Generator_Coverage.exe

all:      OPT := -O3
debug:    OPT := -O0 -g
coverage: OPT := -O0 -g --coverage

#Directory Paths
CONTRIB := ../../../contrib
LIB     := ../../../lib

#Add source files as necessary
SRC := \
				mesh_converter.f95 	      	\
				globals.f95        					\
				thor_mesh.f95      					\
				gmesh.f95

OBJ := $(SRC:.f95=.o)
MOD := $(OBJ:.o=.mod)

#Add contrib files
CONT_SRC :=\
						$(CONTRIB)/tet_mesh_tet_neighbors/tet_mesh_tet_neighbors.f95

CONT_OBJ := $(patsubst %.f95,%.o,$(notdir $(CONT_SRC)))

CONT_MOD := $(CONT_OBJ:.o=.mod)

#Add lib files
LIB_SRC := \
						$(LIB)/integer_array_tools/integer_array_tools.f95

LIB_OBJ := $(patsubst %.f95,%.o,$(notdir $(LIB_SRC)))

LIB_MOD := $(LIB_OBJ:.o=.mod)

#Phony targets for cleaning and building
.PHONY: all clean reset

all: $(EXE)

debug: $(EXE)

coverage: $(EXE)

clean:
	rm -f $(OBJ) $(MOD) $(CONT_OBJ) $(CONT_MOD) $(LIB_OBJ) $(LIB_MOD)

reset:
	rm -f $(OBJ) $(CONT_OBJ) $(LIB_OBJ)
	rm -f $(MOD) $(CONT_MOD) $(LIB_MOD)
	rm -f *.gcno *.gcda
	rm -f -r $(EXE)*.dSYM
	rm -f $(EXE)*

#Generics for source files
$(OBJ): %.o:	%.f95
	$(FC) -c $(OPT) $<

$(MOD):	%.mod:	%.f95
	$(FC) -c $(OPT) $<

$(EXE): $(OBJ) $(CONT_OBJ) $(LIB_OBJ)
	$(FC) -o $@ $(OPT) $(OBJ) $(CONT_OBJ) $(LIB_OBJ)

#Generics for contrib files
$(CONT_OBJ): %.o: $(filter %.f95, $(CONT_SRC))
	$(FC) -c $(OPT) $<

$(CONT_MOD):	%.mod:	$(filter %.f95, $(CONT_SRC))
	$(FC) -c $(OPT) $<

#Generics for lib files
$(LIB_OBJ):	%.o: $(filter %.f95, $(LIB_SRC))
	$(FC) -c $(OPT) $<

$(LIB_MOD): %.mod:	$(filter %.f95, $(LIB_SRC))
	$(FC) -c $(OPT) $<

#Dependency List
#Use [gfortran -M -cpp *.f95] repeatedly until clean compile to update rules below
globals.o globals.mod: globals.f95 tet_mesh_tet_neighbors.mod integer_array_tools.mod
gmesh.o gmesh.mod: gmesh.f95 globals.mod
mesh_converter.o: mesh_converter.f95 globals.mod gmesh.mod thor_mesh.mod integer_array_tools.mod
thor_mesh.o thor_mesh.mod: thor_mesh.f95 tet_mesh_tet_neighbors.mod globals.mod
