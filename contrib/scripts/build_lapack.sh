#!/bin/bash

if [ -d "../libs" ]; then
  cd ../libs
  cd -
else
  mkdir ../libs
fi

if [ -d "../lapack" ]; then
  cp make.inc ../lapack
  cd ../lapack
  make clean
  make
  mv *a ../libs
else
  echo "Error contrib/lapack does not exist"
fi
