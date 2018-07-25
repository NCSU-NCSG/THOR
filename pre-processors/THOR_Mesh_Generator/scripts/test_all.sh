message1=">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
message2="<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n\n"

printf $message1
echo Test 1: basic_cube_mesh_test
cd ../test/basic_cube_mesh_test
  pwd
  bash ./run_test.sh
  cd -
printf $message2

printf $message1
echo Test 2: homogeneous_domain
cd ../test/homogeneous_domain
  pwd
  bash ./run_test.sh
  cd -
printf $message2

printf $message1
echo Test 3: bad_invocation_no_-i
cd ../test/bad_invocation_no_-i
  pwd
  bash ./run_test.sh
cd -
printf $message2

printf $message1
echo Test 4: bad_invocation_wrong_arg
cd ../test/bad_invocation_wrong_arg
  pwd
  bash ./run_test.sh
cd -
printf $message2

printf $message1
echo Test 5: bad_invocation_no_input_file
cd ../test/bad_invocation_no_input_file
  pwd
  bash ./run_test.sh
cd -
printf $message2

printf $message1
echo Test 6: bad_invocation_unsupported_file_type_unknown
cd ../test/bad_invocation_unsupported_file_type_unknown
  pwd
  bash ./run_test.sh
cd -
printf $message2

printf $message1
echo Test 7: bad_invocation_unsupported_file_type_exodus
cd ../test/bad_invocation_unsupported_file_type_exodus
  pwd
  bash ./run_test.sh
cd -
printf $message2

printf $message1
echo Test 8: bad_gmesh_non_tet_element
cd ../test/bad_gmesh_non_tet_element
  pwd
  bash ./run_test.sh
cd -
printf $message2

printf $message1
echo Test 9: bad_gmesh_non_tri_face
cd ../test/bad_gmesh_non_tri_face
  pwd
  bash ./run_test.sh
cd -
printf $message2

printf $message1
echo Test 10: bad_gmesh_no_elements_block
cd ../test/bad_gmesh_no_elements_block
  pwd
  bash ./run_test.sh
cd -
printf $message2

printf $message1
echo Test 11: bad_gmesh_no_nodes_block
cd ../test/bad_gmesh_no_nodes_block
  pwd
  bash ./run_test.sh
cd -
printf $message2

printf $message1
echo Test 12: bad_gmesh_no_format_block
cd ../test/bad_gmesh_no_format_block
  pwd
  bash ./run_test.sh
cd -
printf $message2

printf $message1
echo Test 13: unv_sphere_in_shell_in_box
cd ../test/unv_sphere_in_shell_in_box
  pwd
  bash ./run_test.sh
cd -
printf $message2

printf $message1
echo Test 14: bad_unv_sphere_in_shell_in_box_no_2411
cd ../test/bad_unv_sphere_in_shell_in_box_no_2411
  pwd
  bash ./run_test.sh
cd -
printf $message2

printf $message1
echo Test 15: convert_old_to_new_THOR
cd ../test/convert_old_to_new_THOR
  pwd
  bash ./run_test.sh
cd -
printf $message2
