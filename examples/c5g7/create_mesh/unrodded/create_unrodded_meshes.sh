#!/bin/bash

# r1 files
../../../yak/contrib/instant/instant_mesh_generator-opt c5g7_instant_r1.xml
mpiexec -np 12 ~/projects/make_C5G7_meshes_THOR_study/yak/yak-opt -i extrude_c5g7_unrodded_simple_a1.i Mesh/file=c5g7-2d_r1.e UserObjects/dump_mesh/filename=raw_mesh.dat
python TriPrismToTet.py raw_mesh.dat c5g7_r1a1.thor 4

mpiexec -np 12 ~/projects/make_C5G7_meshes_THOR_study/yak/yak-opt -i extrude_c5g7_unrodded_simple_a2.i Mesh/file=c5g7-2d_r1.e UserObjects/dump_mesh/filename=raw_mesh.dat
python TriPrismToTet.py raw_mesh.dat c5g7_r1a2.thor 8

mpiexec -np 12 ~/projects/make_C5G7_meshes_THOR_study/yak/yak-opt -i extrude_c5g7_unrodded_simple_a3.i Mesh/file=c5g7-2d_r1.e UserObjects/dump_mesh/filename=raw_mesh.dat
python TriPrismToTet.py raw_mesh.dat c5g7_r1a3.thor 16

# r2 files
../../../yak/contrib/instant/instant_mesh_generator-opt c5g7_instant_r2.xml
mpiexec -np 12 ~/projects/make_C5G7_meshes_THOR_study/yak/yak-opt -i extrude_c5g7_unrodded_simple_a1.i Mesh/file=c5g7-2d_r2.e UserObjects/dump_mesh/filename=raw_mesh.dat
python TriPrismToTet.py raw_mesh.dat c5g7_r2a1.thor 4

mpiexec -np 12 ~/projects/make_C5G7_meshes_THOR_study/yak/yak-opt -i extrude_c5g7_unrodded_simple_a2.i Mesh/file=c5g7-2d_r2.e UserObjects/dump_mesh/filename=raw_mesh.dat
python TriPrismToTet.py raw_mesh.dat c5g7_r2a2.thor 8

mpiexec -np 12 ~/projects/make_C5G7_meshes_THOR_study/yak/yak-opt -i extrude_c5g7_unrodded_simple_a3.i Mesh/file=c5g7-2d_r2.e UserObjects/dump_mesh/filename=raw_mesh.dat
python TriPrismToTet.py raw_mesh.dat c5g7_r2a3.thor 16

# r3 files
../../../yak/contrib/instant/instant_mesh_generator-opt c5g7_instant_r3.xml
mpiexec -np 12 ~/projects/make_C5G7_meshes_THOR_study/yak/yak-opt -i extrude_c5g7_unrodded_simple_a1.i Mesh/file=c5g7-2d_r3.e UserObjects/dump_mesh/filename=raw_mesh.dat
python TriPrismToTet.py raw_mesh.dat c5g7_r3a1.thor 4

mpiexec -np 12 ~/projects/make_C5G7_meshes_THOR_study/yak/yak-opt -i extrude_c5g7_unrodded_simple_a2.i Mesh/file=c5g7-2d_r3.e UserObjects/dump_mesh/filename=raw_mesh.dat
python TriPrismToTet.py raw_mesh.dat c5g7_r3a2.thor 8

mpiexec -np 12 ~/projects/make_C5G7_meshes_THOR_study/yak/yak-opt -i extrude_c5g7_unrodded_simple_a3.i Mesh/file=c5g7-2d_r3.e UserObjects/dump_mesh/filename=raw_mesh.dat
python TriPrismToTet.py raw_mesh.dat c5g7_r3a3.thor 16

rm raw_mesh.dat
