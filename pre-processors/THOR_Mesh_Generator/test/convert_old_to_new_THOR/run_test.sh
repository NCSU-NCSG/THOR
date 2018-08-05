rm -f takeda.thrm
../../Thor_Mesh_Generator_MP.exe -i takeda_conversion.in
echo "Diff should be zero. Confirm below."
diff ref.thrm takeda.thrm
