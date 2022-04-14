MODULE termination_module
  !***********************************************************************
  ! This module contains subroutines for terminating the execution of
  ! THOR either successfully or unsuccessfully.
  !***********************************************************************
  USE global_variables
  IMPLICIT NONE

CONTAINS

  !------------------------------------------------------------------------------------------------------------!
  !------------------------------------------------------------------------------------------------------------!

  SUBROUTINE cleanup
    !**********************************************************************
    !
    ! Deallocates global allocatable variables
    !
    !**********************************************************************

    ! local variables

    INTEGER(kind=li) :: m,eg,q,octant,f

    ! src types
    IF ( ALLOCATED( ext_src   ) ) DEALLOCATE( ext_src)

    ! cross section types

    IF( ALLOCATED(xs_mat) )     DEALLOCATE(xs_mat)
    IF( ALLOCATED(eg_bounds) )  DEALLOCATE(eg_bounds)
    IF( ALLOCATED(scat_mult)  ) DEALLOCATE(scat_mult)

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

    ! post-processing variables
    IF( ALLOCATED(point_flux_locations) ) DEALLOCATE(point_flux_locations)

    ! Close files that might be open

    IF(page_sweep.EQ.1_li) CLOSE(unit=99)
    IF(page_iflw .EQ.1_li) CLOSE(unit=97)

  END SUBROUTINE cleanup

  !------------------------------------------------------------------------------------------------------------!
  !------------------------------------------------------------------------------------------------------------!

  SUBROUTINE stop_thor(scode, message)
    !**********************************************************************
    !
    ! Terminates THOR and prints message according to scode
    !
    !**********************************************************************

    ! Pass argument

    INTEGER(kind=li) :: scode
    CHARACTER(*), OPTIONAL:: message

    ! Print message
    SELECT CASE(scode)

    CASE(1_li)
      IF (rank .EQ. 0) THEN
        WRITE(6,*)
        WRITE(6,*) "--------------------------------------------------------"
        WRITE(6,*) "   Execution of THOR completed successfully  "
        WRITE(6,*) "--------------------------------------------------------"
        WRITE(6,*)
      END IF
    CASE(2_li)
      WRITE(6,*) '*** Not enough memory ***'
      WRITE(6,*) '>> Execution terminates unsuccessfully!'
    CASE(3_li)
      WRITE(6,*) "Upstream face moment determination failed."
      WRITE(6,*) '>> Execution terminates unsuccessfully!'
    CASE(4_li)
      WRITE(6,*) "Downstream face moment determination failed."
      WRITE(6,*) '>> Execution terminates unsuccessfully!'
    CASE(5_li)
      WRITE(6,*) "Incoming face moments determination failed."
      WRITE(6,*) '>> Execution terminates unsuccessfully!'
    CASE(6_li)
      WRITE(6,*) "Outgoing face moments determination failed."
      WRITE(6,*) '>> Execution terminates unsuccessfully!'
    CASE(7_li)
      WRITE(6,*) 'Face moments transformation failed.'
      WRITE(6,*) '>> Execution terminates unsuccessfully!'
    CASE(8_li)
      WRITE(6,*) "Please select problem type 1 (eigenvalue search)",&
            "or 2 (external source) to perform calculation"
      WRITE(6,*) '>> Execution terminates unsuccessfully!'
    CASE(9_li)
      WRITE(6,*) 'Sweep type specification not recognized. Execution terminates.'
      WRITE(6,*) '>> Execution terminates unsuccessfully!'
    CASE(10_li)
      WRITE(6,*) "Unacceptable splitting in cell"
      WRITE(6,*) '>> Execution terminates unsuccessfully!'
    CASE(11_li)
      WRITE(6,*) "Unacceptable case from cell splitting in cell"
      WRITE(6,*) '>> Execution terminates unsuccessfully!'
    CASE(12_li)
      WRITE(6,*) 'Some boundary faces feature a fixed inflow flux but the inflow flag is set to false. Execution terminates.'
      WRITE(6,*) '>> Execution terminates unsuccessfully!'
    CASE(13_li)
      WRITE(6,*) 'Precomputed sweep path inconsistent. Execution terminates'
      WRITE(6,*) '>> Execution terminates unsuccessfully!'
    CASE(14_li)
      WRITE(6,*)  'Reflective boundary face could not be determined.'
      WRITE(6,*) '>> Execution terminates unsuccessfully!'
    CASE(15_li)
      WRITE(6,*) 'SLC quadrature only available for orders 1, 2, 3 and 5'
      WRITE(6,*) '>> Execution terminates unsuccessfully!'
    CASE(16_li)
      WRITE(6,*) 'A fixed inflow flux value is placed in a boundary face that is not declared fixed inflow.'
      WRITE(6,*) '>> Execution terminates unsuccessfully!'
    CASE(17_li)
      WRITE(6,*) "Only lambda=-1 is allowed!"
      WRITE(6,*) '>> Execution terminates unsuccessfully!'
    CASE(18_li)
      WRITE(6,*) "Number of Krylov Iterations between restarts must be greater than 0"
      WRITE(6,*) '>> Execution terminates unsuccessfully!'
    CASE(19_li)
      WRITE(6,*) "Maximum number of krylov iterations must be greater than number of",&
            "iterations between restarts."
      WRITE(6,*) '>> Execution terminates unsuccessfully!'
    CASE(20_li)
      WRITE(6,*) "Method has to be 1 (Outer iteration with lagged upscattering)",&
            " 2(Flat iteration) or 3(Flat iteration with updated downscattering)."
      WRITE(6,*) '>> Execution terminates unsuccessfully!'
    CASE(21_li)
      WRITE(6,*) 'Counter fail to correctly count number of moments!'
      WRITE(6,*) '>> Execution terminates unsuccessfully!'
    CASE(22_li)
      WRITE(6,*) 'Higher Order Spherical Harmonic Not Available!'
      WRITE(6,*) '>> Execution terminates unsuccessfully!'
    CASE(23_li)
      WRITE(6,*) 'Legendre Scattering Order Limited to P7'
      WRITE(6,*) '>> Execution terminates unsuccessfully!'
    CASE(24_li)
      WRITE(6,*) 'Reading initial guess file failed'
      WRITE(6,*) '>> Execution terminates unsuccessfully!'
    CASE(25_li)
      WRITE(6,*) "JFNK module requires niter to be equal to 2. This is a coding mistake!"
      WRITE(6,*) '>> Execution terminates unsuccessfully!'
    CASE(26_li)
      WRITE(6,*) 'Method needs to be 1, 2 or 3.'
      WRITE(6,*) '>> Execution terminates unsuccessfully!'
    CASE(27_li)
      WRITE(6,*) 'GMRES error'
      WRITE(6,*) '>> Execution terminates unsuccessfully!'
    CASE(28_li)
      WRITE(6,*) 'Only lambda=0 is allowed'
      WRITE(6,*) '>> Execution terminates unsuccessfully!'
    CASE(29_li)
      WRITE(6,*) 'Inflow file only allowed for fixed source problems'
      WRITE(6,*) '>> Execution terminates unsuccessfully!'
    CASE(30_li)
      WRITE(6,*) 'VTK source file only allowed for fixed source problems'
      WRITE(6,*) '>> Execution terminates unsuccessfully!'
    CASE(31_li)
      WRITE(6,*) 'JFNK option in combination with reflective boundary conditions on opposite faces is not permitted'
      WRITE(6,*) '>> Execution terminates unsuccessfully!'
    CASE(32_li)
      WRITE(6,*) 'Call to associated Legendre Polynomials features impossible combination of l and m.'
      WRITE(6,*) 'This is a coding mistake. Contact the developers.'
      WRITE(6,*) '>> Execution terminates unsuccessfully!'
    CASE(33_li)
      WRITE(6,*) 'Associated Legendre Polynomial order limited to 5.'
      WRITE(6,*) '>> Execution terminates unsuccessfully!'
    CASE(34_li)
      WRITE(6,*) 'Density factors were requested but referenced file was not found.'
      WRITE(6,*) '>> Execution terminates unsuccessfully!'
    CASE default
      IF (PRESENT(message)) THEN
        WRITE(6,*) message
      ELSE
        WRITE(6,*) "No Error message given - Format should be <ERR #> <Messsage>"
      END IF
    END SELECT


    ! Cleanup and terminate

    CALL cleanup
    STOP

  END SUBROUTINE stop_thor

END MODULE termination_module
