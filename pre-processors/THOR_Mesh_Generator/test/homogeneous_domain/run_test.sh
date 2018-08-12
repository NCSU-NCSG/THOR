rm homogeneous.thrm homogeneous_from_exodus.thrm 

../../Thor_Mesh_Generator_MP.exe -i ./homogeneous.in
diff homogeneous.thrm homogeneous_ref.thrm
echo "Diff should be zero. Confirm below."
echo "Diff end"
../../Thor_Mesh_Generator_MP.exe -i ./homogeneous_from_exodus.in
diff homogeneous_from_exodus.thrm homogeneous_from_exodus_ref.thrm
echo "Diff should be zero. Confirm below."
echo "Diff end"

rm homogeneous.thrm homogeneous_from_exodus.thrm

../../Thor_Mesh_Generator_MP.exe -i ./homogeneous.in -r 1
diff homogeneous.thrm homogeneous_ref_r1.thrm
echo "Diff should be zero. Confirm below."
echo "Diff end"
../../Thor_Mesh_Generator_MP.exe -i ./homogeneous_from_exodus.in -r 1
diff homogeneous_from_exodus.thrm homogeneous_from_exodus_ref_r1.thrm
echo "Diff should be zero. Confirm below."
echo "Diff end"
