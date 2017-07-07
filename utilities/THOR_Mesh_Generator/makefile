
FC := gfortran

EXE := Thor_Mesh_Generator

#Add source files as necessary
SRC := \
      geompack.f95       \
      tet_mesh_tet_neighbors.f95 \
			integer_array_tools.f95 \
			mesh_converter.f95 \
			globals.f95        \
			thor_mesh.f95      \
			gmsh.f95

OBJ := $(SRC:.f95=.o)
MOD := $(OBJ:.o=.mod)

.PHONY: all clean reset

all: $(EXE)

clean:
	rm -f $(OBJ) $(MOD)

reset:
	rm -f $(OBJ)
	rm -f $(MOD)
	rm -f $(EXE)

%.o:	%.f95
	$(FC) -c $<

%.mod:	%.f95
	$(FC) -c $<

$(EXE): $(OBJ)
	$(FC) -o $@ $(OBJ)

#Use [gfortran -M -cpp *.f95] repeatedly until clean compile to update rules below
geompack.o geompack.mod: geompack.f95 globals.mod
tet_mesh_tet_neighbors.o: tet_mesh_tet_neighbors.f95
integer_array_tools.o integer_array_tools.mod: integer_array_tools.f95
globals.o globals.mod: globals.f95 integer_array_tools.mod
gmsh.o gmsh.mod: gmsh.f95 globals.mod
mesh_converter.o: mesh_converter.f95 globals.mod gmsh.mod thor_mesh.mod integer_array_tools.mod
thor_mesh.o thor_mesh.mod: thor_mesh.f95 globals.mod
