#!/bin/bash

codeexec="../../../thor-1.0.exe"

for folder in */ ; do
  cd $folder
  pwd
  for file in *.inp ; do
    echo $file
    mpiexec -np 24 $codeexec $file | tee $file".out"
  done
  cd ../
done
