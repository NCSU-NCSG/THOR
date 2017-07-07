#!/bin/bash
#mpiexec -np 4 ../../thor-1.0.exe  cube_eig_pi.i 1 > cube_eig_pi_ua.out
#mpiexec -np 4 ../../thor-1.0.exe  cube_eig_pi.i 2 > cube_eig_pi_fs.out
#mpiexec -np 4 ../../thor-1.0.exe  cube_eig_pi.i 3 > cube_eig_pi_ch.out

#mpiexec -np 4 ../../thor-1.0.exe  cube_eig_jfnk-1-I2.i > cube_eig_jfnk-1-I2.out
mpiexec -np 4 ../../thor-1.0.exe  cube_eig_jfnk-1-I4.i > cube_eig_jfnk-1-I4.out
mpiexec -np 4 ../../thor-1.0.exe  cube_eig_jfnk-2.i > cube_eig_jfnk-2.out
