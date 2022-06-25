#-----------------------------------------------------------------------------------
#Compile Options
#-----------------------------------------------------------------------------------
FC := gfortran
EXE := XS_Converter
EXT := .exe

all:			TYPE :=
MPparallel:    TYPE := _MP
debug:    TYPE := _Debug
catch:    TYPE := _Catch
coverage: TYPE := _Coverage

all:      OPT := -O3
MPparallel: OPT := -O3 -fopenmp
debug:    OPT := -O0 -g -fbacktrace -ffpe-trap=zero,overflow,underflow
catch:    OPT := -O0 -g -Wall -W -Wsurprising -Werror
coverage: OPT := -O0 -g --coverage

#-----------------------------------------------------------------------------------
#Directory Paths
#-----------------------------------------------------------------------------------
CONTRIB := ../../../contrib

#-----------------------------------------------------------------------------------
#Add source files as necessary
#If the files are not compiled using the generic rules, add commands under their
#rule declaration. Add these items to FILTER
#-----------------------------------------------------------------------------------
SRC := \
				main.f90 \
				globals.f90 \
        precision.f90 \
				outfuncs.f90 \
				infuncs.f90

OBJ := $(SRC:.f90=.o)
MOD := $(OBJ:.o=.mod)

FILTER := \

OBJ_FILTER := $(FILTER:.f90=.o)
MOD_FILTER := $(FILTER:.f90=.mod)

#h5 datafiles inclusion
LIB=usr/lib/x86_64-linux-gnu/hdf5/serial
INC=usr/include/hdf5/serial

#-----------------------------------------------------------------------------------
#Add simple compile contrib files
#-----------------------------------------------------------------------------------
CONT_SRC := \
            $(CONTRIB)/stringmod/stringmod.f90

CONT_OBJ := $(patsubst %.f90,%.o,$(notdir $(CONT_SRC)))

CONT_MOD := $(CONT_OBJ:.o=.mod)

#-----------------------------------------------------------------------------------
#Phony targets for cleaning and building
#-----------------------------------------------------------------------------------
.PHONY: all clean reset

all: $(EXE)

MPparallel: $(EXE)

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
	rm -f *.gcno *.gcda
	rm -f -r $(EXE)*.dSYM
	rm -f ../$(EXE)*

#-----------------------------------------------------------------------------------
#Generics for source files
#-----------------------------------------------------------------------------------
$(filter-out $(OBJ_FILTER), $(OBJ)): %.o:	%.f90
	$(FC) -c $(OPT) $< -I/$(INC) -L/$(LIB) -lhdf5_fortran -lhdf5

$(filter-out $(MOD_FILTER), $(MOD)):	%.mod:	%.f90
	$(FC) -c $(OPT) $< -I/$(INC) -L/$(LIB) -lhdf5_fortran -lhdf5

$(EXE): $(OBJ) $(CONT_OBJ) $(LIB_OBJ)
	$(FC) -o $@$(TYPE)$(EXT) $(OPT) $(OBJ) $(CONT_OBJ) $(LIB_OBJ) -I/$(INC) -L/$(LIB) -lhdf5_fortran -lhdf5
	mv ./$(EXE)$(TYPE)$(EXT) ../

#-----------------------------------------------------------------------------------
#Generics for contrib files
#-----------------------------------------------------------------------------------
$(CONT_OBJ): %.o: $(filter %.f90, $(CONT_SRC))
	$(FC) -c $(OPT) $^ -I/$(INC) -L/$(LIB) -lhdf5_fortran -lhdf5

$(CONT_MOD):	%.mod:	$(filter %.f90, $(CONT_SRC))
	$(FC) -c $(OPT) $^ -I/$(INC) -L/$(LIB) -lhdf5_fortran -lhdf5

#-----------------------------------------------------------------------------------
#Generics for lib files
#-----------------------------------------------------------------------------------
$(LIB_OBJ):	%.o: $(filter %.f90, $(LIB_SRC))
	$(FC) -c $(OPT) $^ -I/$(INC) -L/$(LIB) -lhdf5_fortran -lhdf5

$(LIB_MOD): %.mod:	$(filter %.f90, $(LIB_SRC))
	$(FC) -c $(OPT) $^ -I/$(INC) -L/$(LIB) -lhdf5_fortran -lhdf5


#-----------------------------------------------------------------------------------
#Rules for entries in COMP_DEP. Don't forget to add them to make clean / reset
#-----------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------
#Dependency List
#Use [gfortran -M -cpp *.f90] repeatedly until clean compile to update rules below
#-----------------------------------------------------------------------------------
globals.o globals.mod: globals.f90 precision.mod
infuncs.o infuncs.mod: infuncs.f90 globals.mod stringmod.mod
main.o: main.f90 globals.mod infuncs.mod outfuncs.mod
outfuncs.o outfuncs.mod: outfuncs.f90 globals.mod
precision.o precision.mod: precision.f90