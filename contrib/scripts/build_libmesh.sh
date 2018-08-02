#!/usr/bin/env bash

ulimit -l unlimited
ulimit -n 10240

if [ -d "../libmesh" ]; then
  rm -rf ../libmesh
fi

cd ../
git clone git@github.com:libMesh/libmesh.git
cd libmesh
mkdir build
cd build
../configure --with-methods="opt dbg devel"
make -j $1

if grep "THOR_LIBMESH_DIRECTORY" $HOME/.bash_profile > /dev/null
then
  echo THOR_LIBMESH_DIRECTORY already defined in bash_profile
else
  echo "export THOR_LIBMESH_DIRECTORY=$PWD" >> $HOME/.bash_profile
fi
