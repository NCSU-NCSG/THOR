rm -f unvtest.thrm
../../Thor_Mesh_Generator_MP.exe -i ./unvtest.in
echo "Diff should be zero. Confirm below."
diff ref.thrm unvtest.thrm