cd ../test/basic_cube_mesh_test
  pwd
  bash ./run_test.sh
  cd -

cd ../test/homogeneous_domain
  pwd
  bash ./run_test.sh
  cd -

cd ../test/bad_invocation_no_-i
  pwd
  bash ./run_test.sh
cd -

cd ../test/bad_invocation_wrong_arg
  pwd
  bash ./run_test.sh
cd -

cd ../test/bad_invocation_no_input_file
  pwd
  bash ./run_test.sh
cd -

cd ../test/bad_invocation_unsupported_file_type_unknown
  pwd
  bash ./run_test.sh
cd -

cd ../test/bad_invocation_unsupported_file_type_exodus
  pwd
  bash ./run_test.sh
cd -

cd ../test/bad_gmesh_non_tet_element
  pwd
  bash ./run_test.sh
cd -

cd ../test/bad_gmesh_non_tri_face
  pwd
  bash ./run_test.sh
cd -

cd ../test/bad_gmesh_no_elements_block
  pwd
  bash ./run_test.sh
cd -

cd ../test/bad_gmesh_no_nodes_block
  pwd
  bash ./run_test.sh
cd -

cd ../test/bad_gmesh_no_format_block
  pwd
  bash ./run_test.sh
cd -

cd ../test/unv_sphere_in_shell_in_box
  pwd
  bash ./run_test.sh
cd -

cd ../test/bad_unv_sphere_in_shell_in_box_no_2411
  pwd
  bash ./run_test.sh
cd -
