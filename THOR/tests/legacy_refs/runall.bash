#!/bin/bash

codeexec="/mnt/c/Users/nfherrin/Documents/THOR/THOR/thor-1.0.exe"

for folder in */ ; do
  cd $folder
  pwd
  for file in *.inp ; do
    echo $file
    $codeexec $file > $file".out"
  done
  cd ../
done
