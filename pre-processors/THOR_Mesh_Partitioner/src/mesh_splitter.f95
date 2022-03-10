PROGRAM MESH_SPLIT

  USE global_variables
  USE metis_module
  USE vtk_module
  USE reconstruction_module
  USE statistics_module
  USE output_module

  IMPLICIT NONE

  !! Start execution timer and clear the screen
  CALL CPU_TIME(t_start)
  CALL SYSTEM('clear')

  !!Check for the correct invocation
  temp_int = COMMAND_ARGUMENT_COUNT()
  IF ((temp_int .NE. 4)) STOP 'Provide root file name from [root name].thrm, # procs, -s or -n for stats, and seed'

  !! Get provided parameters and convert #procs to an integer
  CALL GET_COMMAND_ARGUMENT (1, root_name)
  CALL GET_COMMAND_ARGUMENT (2, num_procs_str)
  CALL GET_COMMAND_ARGUMENT (3, stats)
  CALL GET_COMMAND_ARGUMENT (4, seed)
  
  READ(num_procs_str,*) num_procs
  WRITE(*,'(/,A,I0)') "Target Processor Count: ",num_procs

  !! Generate necessary file structure
  CALL file_structure_check()

  !! Build the input THRM Name
  WRITE(thrm_file,"(A)") TRIM(thrm_path)//TRIM(root_name)//".thrm"
  PRINT *, thrm_file

  !! Read in the mesh element list to prepare for METIS run
  CALL get_mesh_element_list()

  !! Write the metis input file
  CALL generate_metis_input()

  !! Call Metis
  CALL execute_metis()

  !! Generate VTK file name and call
  WRITE(vtk_file, "(A)") "./vtk/"//TRIM(root_name)//"-"//TRIM(num_procs_str)//".vtu"

  CALL generate_vtk(TRIM(thrm_file), TRIM(metis_elem_out_file), TRIM(vtk_file))

  !! End timer for all operations prior to reconstruction
  CALL CPU_TIME(t_end1)

  !! Combine METIS output and oringal .THRM file to create fully specified
  !!  partition files
  CALL reconstruct_meshes()

  !! Caluclate basic statistics for the new meshes
  IF (TRIM(stats).EQ.'-s') CALL calculate_stats()

  !! Write out the new <p> separate mesh files
  CALL write_files()

  !! Generate a unified tarball
  CALL pack_output()

  !! Print summary of execution for user
  CALL execution_summary(6)

END PROGRAM

SUBROUTINE file_structure_check()
  USE global_variables
  IMPLICIT NONE

  INQUIRE(FILE = "./meshes", EXIST = temp_log)
  IF (.NOT. temp_log) CALL SYSTEM("mkdir ./meshes")
  INQUIRE(FILE = "./metis_in", EXIST = temp_log)
  IF (.NOT. temp_log) CALL SYSTEM("mkdir ./metis_in")
  INQUIRE(FILE = "./metis_out", EXIST = temp_log)
  IF (.NOT. temp_log) CALL SYSTEM("mkdir ./metis_out")
  INQUIRE(FILE = "./vtk", EXIST = temp_log)
  IF (.NOT. temp_log) CALL SYSTEM("mkdir ./vtk")
  INQUIRE(FILE = "./part_meshes", EXIST = temp_log)
  IF (.NOT. temp_log) CALL SYSTEM("mkdir ./part_meshes")
  INQUIRE(FILE = "./final_output", EXIST = temp_log)
  IF (.NOT. temp_log) CALL SYSTEM("mkdir ./final_output")
END SUBROUTINE file_structure_check

SUBROUTINE get_mesh_element_list()
  USE global_variables
  IMPLICIT NONE
  
  !!Strip out relevant THM data to make METIS input
  OPEN(UNIT = 21, FILE = thrm_file, ACTION = "READ", STATUS= "OLD")
  READ(21,*) num_node
  READ(21,*) num_elem
  ALLOCATE(elem_list(num_elem,4))
  READ(21,*) num_cell_blocks
  READ(21,*) num_side_sets
  READ(21,*)com_check
  IF (com_check.NE.'!') BACKSPACE(21)
  DO i = 1, num_node
    READ(21,*)
  END DO
  READ(21,*)com_check
  IF (com_check.NE.'!') BACKSPACE(21)
  DO i = 1, num_elem
    READ(21,*)
  END DO
  READ(21,*)com_check
  IF (com_check.NE.'!') BACKSPACE(21)
  DO i = 1, num_elem
    READ(21,*) temp_int, elem_list(i, :)
  END DO
  CLOSE(21)

END SUBROUTINE get_mesh_element_list

SUBROUTINE pack_output()
  use global_variables
  IMPLICIT NONE
  CALL SYSTEM("mkdir ./final_output/"//TRIM(root_name)//"-"//TRIM(num_procs_str))
  OPEN(UNIT=95, &
  & FILE='./final_output/'//TRIM(root_name)//"-"//TRIM(num_procs_str)//'/stats.txt'& 
  & ,STATUS='REPLACE',ACTION='WRITE',iostat=ierr)
  IF (ierr.NE.0) THEN
    WRITE(*,*)'FATAL ERROR: Could not create statistics file.'
    STOP
  END IF
  CALL execution_summary(95)
  CALL SYSTEM("cp -r ./part_meshes ./final_output/"//TRIM(root_name)//"-"//TRIM(num_procs_str)//"/")
  CALL SYSTEM("cp ./meshes/"//TRIM(root_name)//".thrm"//" ./final_output/"//TRIM(root_name)//"-"//TRIM(num_procs_str))
  CALL SYSTEM("cp ./vtk/"//TRIM(root_name)//"-"//TRIM(num_procs_str)// &
             & ".vtu ./final_output/"//TRIM(root_name)//"-"//TRIM(num_procs_str))
  CALL SYSTEM("cd ./final_output; tar -cf ./"//TRIM(root_name)//"-"//TRIM(num_procs_str) &
             & //".tar " //TRIM(root_name)//"-"//TRIM(num_procs_str))

  CALL SYSTEM("rm -r ./metis_in ./metis_out ./part_meshes")
  CALL SYSTEM("rm -r ./final_output/"//TRIM(root_name)//"-"//TRIM(num_procs_str))
  WRITE(*,'(/)')
END SUBROUTINE pack_output

SUBROUTINE execution_summary(filenum)
  USE global_variables
  IMPLICIT NONE
  
  INTEGER,INTENT(IN)::filenum
    
    IF (TRIM(stats).EQ.'-s') THEN
        IF (filenum.EQ.6) THEN
            IF (contig) THEN
                WRITE(filenum, '(2A)')        "CONTIGUOUS                           : ", ACHAR(027)//"[4;42mTRUE"//ACHAR(027)//"[0m"
            ELSE
                WRITE(filenum, '(2A)')       "CONTIGUOUS                           : ", ACHAR(027)//"[4;41mFALSE"//ACHAR(027)//"[0m"
            END IF
        ELSE
            IF (contig) THEN
                WRITE(filenum, '(2A)')        "CONTIGUOUS                           : ", "TRUE"
            ELSE
                WRITE(filenum, '(2A)')        "CONTIGUOUS                           : ", "FALSE"
            END IF
        END IF
        WRITE(filenum,"(A, F8.2, 2I8)") "# ELEM PER PROCESSOR (AVG/MIN/MAX)   :    ", load_AVG, load_MIN, load_MAX
        WRITE(filenum,"(A, 3F8.5)")     "EXT FACES PER ELEMENT (AVG/MIN/MAX)  :    ", &
        &   ext_int_ratio_AVG, ext_int_ratio_MIN, ext_int_ratio_MAX
        WRITE(filenum,'(A,I0)')         "# OF ELEMENTS WITH 3 BOUNDARY FACES  :    ", num_spikes
        WRITE(filenum,"(A, F8.2, 2I8)") "# OF NEIGHBORS (AVG/MIN/MAX)         :    ", neigh_AVG, neigh_MIN, neigh_MAX
    ELSE
        IF (filenum.EQ.6) THEN
            WRITE(filenum, '(2A)')"CONTIGUOUS                           : ", ACHAR(027)//"[4;43mNOT EVALUATED"//ACHAR(027)//"[0m"
        ELSE 
            WRITE(filenum, '(2A)')"CONTIGUOUS                           : ", "NOT EVALUATED"
        END IF
        WRITE(filenum,*)'Statistics not generated. Use -s option to generate statistics.'
    END IF
    CALL CPU_TIME(t_end2)
    WRITE(filenum, '(/)')
    WRITE(filenum,"(A, F16.8)") "Mesh Decomposition Time      (s) : ", t_end1 - t_start
    WRITE(filenum,"(A, F16.8)") "Mesh File Construction Time  (s) : ", t_end2 - t_end1
    WRITE(filenum,"(A, F16.8)") "Total Time (Excliding METIS) (s) : ", t_end2 - t_start

    mem_usage = mem_usage + SIZEOF(node_proc_map)
    mem_usage = mem_usage + SIZEOF(elem_list)
    mem_usage = mem_usage + SIZEOF(elem_proc_indexes)
    mem_usage = mem_usage + SIZEOF(node_proc_indexes)
    mem_usage = mem_usage + SIZEOF(elem_node_proc_map)
    mem_usage = mem_usage + SIZEOF(element_proc_map)
    mem_usage = mem_usage + SIZEOF(src_proc_indexes)
    mem_usage = mem_usage + SIZEOF(adj_proc_map)
    mem_usage = mem_usage + SIZEOF(adj_proc_indexes)
    mem_usage = mem_usage + SIZEOF(BC_proc_map)
    mem_usage = mem_usage + SIZEOF(BC_proc_indexes)
    mem_usage = mem_usage + SIZEOF(proc_map)

    WRITE(filenum,*) "Allocatable Memory Used: ", (FLOAT(mem_usage)/2.0**20), " MBytes"
END SUBROUTINE execution_summary
