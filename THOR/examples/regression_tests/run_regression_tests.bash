#!/bin/bash

numprocs=$1

re='^[0-9]+$'
if ! [[ $numprocs =~ $re ]] ; then
   numprocs=1
fi

codeexec="../../../thor-1.0.exe"

bash ./rmresults.bash

for folder in */ ; do
  cd $folder
  pwd
  for file in *.inp ; do
    echo $file
    mpiexec -np $numprocs $codeexec $file
  done
  cd ../
done

python3 compare_results.py