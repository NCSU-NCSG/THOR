MODULE statistics_module

  USE global_variables

  IMPLICIT NONE

CONTAINS

  SUBROUTINE calculate_stats()

    IMPLICIT NONE
    !!Generate Mesh Statistics
    !! Load Balance
    !Average / min / max of #elements
    load_AVG = (SUM(elem_proc_indexes)-num_procs)/REAL(num_procs)
    load_MIN = MINVAL(elem_proc_indexes) - 1
    load_MAX = MAXVAL(elem_proc_indexes) - 1
    !!External vs. internal ratio
    !Average / min / max of #external faces per element)
    ! BC_proc_indexes = BC_proc_indexes - 1
    ! elem_proc_indexes = elem_proc_indexes - 1
    IF (num_side_sets .GT. 1) WRITE(*,'(A,A)') ACHAR(027)//"[4;43mWARN::"//ACHAR(027)//"[0m", &
          & "BC metrics do not support meshes with more than 1 side set"
    ext_int_ratio_AVG = SUM(REAL(BC_proc_indexes(:,1)-1) / REAL(elem_proc_indexes-1)) / REAL(num_procs)
    ext_int_ratio_MIN = MINVAL(REAL(BC_proc_indexes(:,1) - 1) / REAL(elem_proc_indexes - 1))
    ext_int_ratio_MAX = MAXVAL(REAL(BC_proc_indexes(:,1) - 1) / REAL(elem_proc_indexes - 1))
    !! # with three external faces
    ! Just count the occurances in the BC list
    !! Contiguous (no cells with 4 external neighbors & #cells>1)
    ! Same check, but # = 4, not 3
      contig = .TRUE.
      num_spikes = 0
      !$omp parallel do default(shared) firstprivate(temp_int,j) private(i)
    DO i = 1, num_procs
      temp_int = 1
      DO j = 1, elem_proc_indexes(i)-1
        temp_int = COUNT(BC_proc_map(i,:,:,1) .EQ. element_proc_map(i,j,1))
        IF (temp_int .EQ. 3) num_spikes = num_spikes + 1
        IF (temp_int .EQ. 4) contig = .FALSE.
      END DO
    END DO
    !$omp end parallel do
    !! # neighbors
    ! AVG / min / max of # of different processors in adj. list - 1
    ALLOCATE(neigh_count(num_procs))
    neigh_count = 0
    !$omp parallel do default(shared) firstprivate(j,temp_int) private(i)
    DO i = 1, num_procs
      DO j = 1, num_procs
        IF (i .EQ. j) CYCLE
        temp_int = COUNT(adj_proc_map(i,:,5) .EQ. j)
        IF (temp_int .GT. 0) neigh_count(i) = neigh_count(i) + 1
      END DO
    END DO
    !$omp end parallel do
    neigh_AVG = SUM(neigh_count) / num_procs
    neigh_MIN = MINVAL(neigh_count)
    neigh_MAX = MAXVAL(neigh_count)
  END SUBROUTINE


END MODULE
