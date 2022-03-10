MODULE metis_module

USE global_variables

IMPLICIT NONE

CONTAINS

  SUBROUTINE generate_metis_input()
    IMPLICIT NONE
    !!Make METIS File
    WRITE(metis_in_file,'(A)')  TRIM(metis_in_path)//(TRIM(root_name)//".met")
    PRINT *, metis_in_file
    OPEN(UNIT = 22, FILE = metis_in_file, ACTION="WRITE", STATUS="REPLACE")
    WRITE(22,*) num_elem
    DO i = 1, num_elem
      WRITE(22,*) elem_list(i, :)
    END DO
  END SUBROUTINE


  SUBROUTINE execute_metis()
    IMPLICIT NONE
    !!Run METIS
    WRITE(metis_elem_out_file, "(A)") "./metis_out/"//TRIM(root_name)//".met.epart."//TRIM(num_procs_str)
    IF (num_procs .GT. 1) THEN
      WRITE(*, '(2/)')
      CALL SYSTEM("mpmetis "//metis_in_file//" "//num_procs_str // " -contig -ncommon=3 -seed="//TRIM(seed))
      CALL SYSTEM("mv ./metis_in/*part* ./metis_out/")
      WRITE(*, '(2/)')

      PRINT *, metis_elem_out_file
      WRITE(metis_node_out_file, "(A)") "./metis_out/"//TRIM(root_name)//".met.npart."//TRIM(num_procs_str)
      PRINT *, metis_node_out_file
    ELSE IF (num_procs .EQ. 1) THEN
      OPEN(UNIT=31, FILE = metis_elem_out_file, ACTION = "WRITE", STATUS = "REPLACE")
      DO i = 1, num_elem
        WRITE(31, '(I0)') 0
      END DO
      CLOSE(31)
    ELSE
      STOP 'INVALID PROC COUNT'
    END IF
    WRITE(*, '(2/)')
  END SUBROUTINE

END MODULE
