#!/bin/bash

diff -q ./all_reflective_6_group_analytical/six_group_keig_analytical_convergence.gold ./all_reflective_6_group_analytical/six_group_keig_analytical_conv.convergence
diff -q ./all_reflective_6_group_analytical/six_group_keig_analytical_csv.gold         ./all_reflective_6_group_analytical/six_group_keig_analytical_out.csv

diff -q ./demonstrate_cartesian_map/cube_cartesian_map.gold      ./demonstrate_cartesian_map/cube_cartesian_map.dat
diff -q ./demonstrate_cartesian_map/cube_eig_pi_convergence.gold ./demonstrate_cartesian_map/cube_eig_pi_conv.convergence
diff -q ./demonstrate_cartesian_map/cube_eig_pi_csv.gold         ./demonstrate_cartesian_map/cube_eig_pi_out.csv

diff -q ./godiva/godiva_c_convergence.gold ./godiva/godiva_c_conv.convergence
diff -q ./godiva/godiva_c_csv.gold         ./godiva/godiva_c_out.csv
diff -q ./godiva/godiva_o_convergence.gold ./godiva/godiva_o_conv.convergence
diff -q ./godiva/godiva_o_csv.gold         ./godiva/godiva_o_out.csv

diff -q ./homogeneous_cube_keig_1G/keig_lc_convergence.gold   ./homogeneous_cube_keig_1G/keig_lc_conv.convergence
diff -q ./homogeneous_cube_keig_1G/keig_lc_csv.gold           ./homogeneous_cube_keig_1G/keig_lc_out.csv
diff -q ./homogeneous_cube_keig_1G/keig_sc_convergence.gold   ./homogeneous_cube_keig_1G/keig_sc_conv.convergence
diff -q ./homogeneous_cube_keig_1G/keig_sc_csv.gold           ./homogeneous_cube_keig_1G/keig_sc_out.csv

diff -q ./homogeneous_cube_source_1G/fixed_source_convergence.gold ./homogeneous_cube_source_1G/fixed_source_conv.convergence
diff -q ./homogeneous_cube_source_1G/fixed_source_csv.gold         ./homogeneous_cube_source_1G/fixed_source_out.csv

diff -q ./kobayashi-1/kobayashi-LC_convergence.gold ./kobayashi-1/kobayashi-LC_conv.convergence
diff -q ./kobayashi-1/kobayashi-LC_csv.gold         ./kobayashi-1/kobayashi-LC_out.csv
diff -q ./kobayashi-1/kobayashi-SC_convergence.gold ./kobayashi-1/kobayashi-SC_conv.convergence
diff -q ./kobayashi-1/kobayashi-SC_csv.gold         ./kobayashi-1/kobayashi-SC_out.csv

diff -q ./kobayashi-3/kobayashi-LC_convergence.gold ./kobayashi-3/kobayashi-LC_conv.convergence
diff -q ./kobayashi-3/kobayashi-LC_csv.gold         ./kobayashi-3/kobayashi-LC_out.csv
diff -q ./kobayashi-3/kobayashi-SC_convergence.gold ./kobayashi-3/kobayashi-SC_conv.convergence
diff -q ./kobayashi-3/kobayashi-SC_csv.gold         ./kobayashi-3/kobayashi-SC_out.csv

diff -q ./multigroup_source/multigroup_source_convergence.gold      ./multigroup_source/multigroup_source_conv.convergence
diff -q ./multigroup_source/multigroup_source_csv.gold              ./multigroup_source/multigroup_source_out.csv
diff -q ./multigroup_source/multigroup_source_fine_convergence.gold ./multigroup_source/multigroup_source_fine_conv.convergence
diff -q ./multigroup_source/multigroup_source_fine_csv.gold         ./multigroup_source/multigroup_source_fine_out.csv

diff -q ./noshield_snap_2group/noshield_2group_convergence.gold ./noshield_snap_2group/noshield_2group_conv.convergence
diff -q ./noshield_snap_2group/noshield_2group_csv.gold         ./noshield_snap_2group/noshield_2group_out.csv

diff -q ./noshield_snap_32group/noshield_32group_convergence.gold ./noshield_snap_32group/noshield_32group_conv.convergence
diff -q ./noshield_snap_32group/noshield_32group_csv.gold         ./noshield_snap_32group/noshield_32group_out.csv

diff -q ./power_iterations/godiva_PI_S4_convergence.gold ./power_iterations/godiva_PI_S4_conv.convergence
diff -q ./power_iterations/godiva_PI_S4_csv.gold         ./power_iterations/godiva_PI_S4_out.csv

diff -q ./simple_cube_eigenvalue/cube_eig_pi_convergence.gold        ./simple_cube_eigenvalue/cube_eig_pi_conv.convergence
diff -q ./simple_cube_eigenvalue/cube_eig_pi_csv.gold                ./simple_cube_eigenvalue/cube_eig_pi_out.csv

diff -q ./takeda-IV/takeda_convergence.gold      ./takeda-IV/takeda_conv.convergence
diff -q ./takeda-IV/takeda_csv.gold              ./takeda-IV/takeda_out.csv

diff -q ./two_region_two_group_P1_keig/2r2gP1_LC_convergence.gold ./two_region_two_group_P1_keig/2r2gP1_LC_conv.convergence
diff -q ./two_region_two_group_P1_keig/2r2gP1_LC_csv.gold         ./two_region_two_group_P1_keig/2r2gP1_LC_out.csv
diff -q ./two_region_two_group_P1_keig/2r2gP1_SC_convergence.gold ./two_region_two_group_P1_keig/2r2gP1_SC_conv.convergence
diff -q ./two_region_two_group_P1_keig/2r2gP1_SC_csv.gold         ./two_region_two_group_P1_keig/2r2gP1_SC_out.csv
