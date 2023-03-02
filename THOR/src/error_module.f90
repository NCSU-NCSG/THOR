MODULE error_module
  !***********************************************************************
  ! This module contains subroutines for terminating the execution of
  ! THOR either successfully or unsuccessfully.
  !***********************************************************************
  USE globals
  USE mpi
  USE stringmod
  IMPLICIT NONE
  PRIVATE
  !
  PUBLIC :: raise_fatal_error, thor_success, raise_warning

CONTAINS

  !------------------------------------------------------------------------------------------------------------!
  !------------------------------------------------------------------------------------------------------------!

  SUBROUTINE cleanup_thor
    !**********************************************************************
    !
    ! Deallocates global allocatable variables
    !
    !**********************************************************************

    ! src types
    IF ( ALLOCATED( ext_src   ) ) DEALLOCATE( ext_src)
    IF ( ALLOCATED( source_ids   ) ) DEALLOCATE( source_ids)

    ! cross section types

    IF( ALLOCATED(xs_mat) )     DEALLOCATE(xs_mat)
    IF( ALLOCATED(eg_bounds) )  DEALLOCATE(eg_bounds)
    IF( ALLOCATED(scat_mult)  ) DEALLOCATE(scat_mult)
    IF( ALLOCATED(material_ids)  ) DEALLOCATE(material_ids)

    ! fixed inflow flux

    IF(ALLOCATED( binflx )) DEALLOCATE( binflx )

    ! mesh data types

    IF( ALLOCATED(vertices) )       DEALLOCATE(vertices)
    IF( ALLOCATED(cells) )          DEALLOCATE(cells)
    IF( ALLOCATED(b_cells) )        DEALLOCATE(b_cells)
    IF( ALLOCATED(rb_cells) )       DEALLOCATE(rb_cells)
    IF( ALLOCATED(vb_cells) )       DEALLOCATE(vb_cells)
    IF( ALLOCATED(fb_cells) )       DEALLOCATE(fb_cells)
    IF( ALLOCATED(adjacency_list) ) DEALLOCATE(adjacency_list)
    IF( ALLOCATED(quadrature) )     DEALLOCATE(quadrature)
    IF( ALLOCATED(index_v) )        DEALLOCATE(index_v)
    IF( ALLOCATED(index_f) )        DEALLOCATE(index_f)
    IF( ALLOCATED(max_error) )      DEALLOCATE(max_error)
    IF( ALLOCATED(Ysh) )            DEALLOCATE(Ysh)
    IF( ALLOCATED(outward_normal))  DEALLOCATE(outward_normal)
    IF( ALLOCATED(refl_face_tpe) )  DEALLOCATE(refl_face_tpe)
    IF( ALLOCATED(sweep_path) )     DEALLOCATE(sweep_path)
    IF( ALLOCATED(reg2mat))         DEALLOCATE(reg2mat)
    IF( ALLOCATED(dens_fact) )      DEALLOCATE(dens_fact)
    IF( ALLOCATED(reg_vol) )        DEALLOCATE(reg_vol)

    ! variables related to rings

    IF( ALLOCATED(neldep) )   DEALLOCATE(neldep)
    IF( ALLOCATED( eldep) )   DEALLOCATE( eldep)

    !cycle variables
    IF( ALLOCATED( is_cycle) )   DEALLOCATE( is_cycle)

    !parallel variables
    IF( ALLOCATED( parallel_map_g2l) )   DEALLOCATE( parallel_map_g2l)
    IF( ALLOCATED( parallel_map_l2g) )   DEALLOCATE( parallel_map_l2g)

    ! post-processing variables
    IF( ALLOCATED(point_flux_locations) ) DEALLOCATE(point_flux_locations)

    ! Close files that might be open

    IF(page_sweep.EQ.1_li) CLOSE(unit=99)
    IF(page_iflw .EQ.1_li) CLOSE(unit=97)

  END SUBROUTINE cleanup_thor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE raise_fatal_error(message)
    !**********************************************************************
    !
    ! Terminates THOR and prints message according to scode
    !
    !**********************************************************************

    ! Pass argument
    CHARACTER(*), OPTIONAL,INTENT(IN) :: message
    INTEGER :: mpi_err

    ! Print message
    IF(rank .EQ. 0)THEN
      CALL printlog("*****************************************************************************")
      CALL printlog("*****************************************************************************")
      CALL printlog("*****************************************************************************")
      IF (PRESENT(message)) THEN
        CALL printlog('FATAL ERROR!')
        CALL printlog('ERROR: '//TRIM(ADJUSTL(message)))
      ELSE
        CALL printlog('FATAL ERROR!')
        CALL printlog('...')
        CALL printlog('No error message given')
      ENDIF
      CALL printlog('>> THOR encountered a fatal error!')
      CALL printlog('>> Execution of THOR terminated UNsuccessfully!')
      CALL printlog("*****************************************************************************")
      CALL printlog("*****************************************************************************")
      CALL printlog("*****************************************************************************")
    ENDIF
    CALL MPI_BARRIER(MPI_COMM_WORLD, mpi_err)
    CALL MPI_FINALIZE(mpi_err)


    ! Cleanup and terminate

    CALL cleanup_thor
    STOP

  END SUBROUTINE raise_fatal_error

  SUBROUTINE thor_success(message)
    !**********************************************************************
    !
    ! Terminates THOR and prints message according to scode
    !
    !**********************************************************************

    ! Pass argument
    CHARACTER(*), OPTIONAL,INTENT(IN) :: message

    IF(rank .EQ. 0)THEN
      IF (PRESENT(message)) THEN
        CALL printlog(TRIM(ADJUSTL(message)))
      END IF
      CALL printlog('')
      CALL printlog("--------------------------------------------------------")
      CALL printlog("   Execution of THOR completed successfully  ")
      IF(num_warnings .GT. 0)THEN
        CALL raise_warning(TRIM(str(num_warnings))//" warning(s) were found during execution.")
        CALL raise_warning("Review of log is recommended to avoid unexpected behavior.")
      ENDIF
      CALL printlog("--------------------------------------------------------")
      CALL printlog('')
    ENDIF


    ! Cleanup and terminate

    CALL cleanup_thor
    STOP

  ENDSUBROUTINE thor_success

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE raise_warning(message)
    !**********************************************************************
    !
    ! Throws a warning for THOR
    !
    !**********************************************************************
    CHARACTER(*), INTENT(IN):: message

    IF(rank .EQ. 0)CALL printlog('WARNING: '//TRIM(ADJUSTL(message)))
    num_warnings=num_warnings+1

  END SUBROUTINE raise_warning
END MODULE error_module
