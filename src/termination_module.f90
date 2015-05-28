module termination_module
!***********************************************************************
! This module contains subroutines for terminating the execution of 
! THOR either successfully or unsuccessfully.
!***********************************************************************
use global_variables
implicit none

contains

!------------------------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------------------------!

  subroutine cleanup
  !**********************************************************************
  !
  ! Deallocates global allocatable variables 
  ! 
  !**********************************************************************
  
    ! local variables
  
      integer(kind=li) :: m,eg,q,octant,f    
  
    ! src types
     
      if ( allocated( src_m   ) ) deallocate( src_m)
      if ( allocated( src_mat ) ) deallocate( src_mat )  
      if ( allocated( src_str ) ) deallocate( src_str )

    ! cross section types
    
      if( allocated(xs_mat) )     deallocate(xs_mat) 
      if( allocated(chi) )        deallocate(chi) 
      if( allocated(eg_bounds) )  deallocate(eg_bounds)
      if( allocated(fiss) )       deallocate(fiss)
      if( allocated(nu) )         deallocate(nu)
      if( allocated(sigma_t) )    deallocate(sigma_t)
      if( allocated(sigma_scat) ) deallocate(sigma_scat)        
      if( allocated(scat_mult)  ) deallocate(scat_mult)
      if( allocated(tsigs) )      deallocate(tsigs)

    ! fixed inflow flux

      if(allocated( binflx )) deallocate( binflx )

    ! mesh data types
 
      if( allocated(vertices) )       deallocate(vertices) 
      if( allocated(cells) )          deallocate(cells) 
      if( allocated(b_cells) )        deallocate(b_cells) 
      if( allocated(rb_cells) )       deallocate(rb_cells) 
      if( allocated(vb_cells) )       deallocate(vb_cells) 
      if( allocated(fb_cells) )       deallocate(fb_cells) 
      if( allocated(adjacency_list) ) deallocate(adjacency_list)
      if( allocated(quadrature) )     deallocate(quadrature)
      if( allocated(index_v) )        deallocate(index_v)     
      if( allocated(index_f) )        deallocate(index_f)     
      if( allocated(max_error) )      deallocate(max_error)
      if( allocated(Ysh) )            deallocate(Ysh) 
      if( allocated(outward_normal))  deallocate(outward_normal)  
      if( allocated(refl_face_tpe) )  deallocate(refl_face_tpe)
      if( allocated(sweep_path) )     deallocate(sweep_path)
      if( allocated(reg2mat))         deallocate(reg2mat)   
      if( allocated(dens_fact) )      deallocate(dens_fact) 
      if( allocated(reg_vol) )        deallocate(reg_vol) 

    ! variables related to rings

      if( allocated(neldep) )   deallocate(neldep)      
      if( allocated( eldep) )   deallocate( eldep)      

    ! Close files that might be open

      if(page_sweep.eq.1_li) close(unit=99)
      if(page_iflw .eq.1_li) close(unit=97)

  end subroutine cleanup

!------------------------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------------------------!

  subroutine stop_thor(scode, message)
  !**********************************************************************
  ! 
  ! Terminates THOR and prints message according to scode
  !
  !**********************************************************************

  ! Pass argument
     
    integer(kind=li) :: scode
    character(250), optional:: message

  ! Print message
    select case(scode)

      case(1_li)
            if (rank .eq. 0) then 
             write(6,*)
             write(6,*) "--------------------------------------------------------"
             write(6,*) "   Execution of THOR completed successfully  "
             write(6,*) "--------------------------------------------------------"
             write(6,*) 
            end if
      case(2_li) 
             write(6,*) '*** Not enough memory ***'
             write(6,*) '>> Execution terminates unsuccessfully!' 
      case(3_li)
             write(6,*) "Upstream face moment determination failed."
             write(6,*) '>> Execution terminates unsuccessfully!' 
      case(4_li)
             write(6,*) "Downstream face moment determination failed."
             write(6,*) '>> Execution terminates unsuccessfully!' 
      case(5_li)
          write(6,*) "Incoming face moments determination failed."
          write(6,*) '>> Execution terminates unsuccessfully!' 
      case(6_li)
          write(6,*) "Outgoing face moments determination failed."
          write(6,*) '>> Execution terminates unsuccessfully!' 
      case(7_li)
          write(6,*) 'Face moments transformation failed.'
          write(6,*) '>> Execution terminates unsuccessfully!' 
      case(8_li)
          write(6,*) "Please select problem type 1 (eigenvalue search)",&
                      "or 2 (external source) to perform calculation"
          write(6,*) '>> Execution terminates unsuccessfully!' 
      case(9_li)
          write(6,*) 'Sweep type specification not recognized. Execution terminates.'
          write(6,*) '>> Execution terminates unsuccessfully!' 
      case(10_li)
          write(6,*) "Unacceptable splitting in cell"
          write(6,*) '>> Execution terminates unsuccessfully!' 
      case(11_li)
          write(6,*) "Unacceptable case from cell splitting in cell"
          write(6,*) '>> Execution terminates unsuccessfully!' 
      case(12_li)
          write(6,*) 'Some boundary faces feature a fixed inflow flux but the inflow flag is set to false. Execution terminates.' 
          write(6,*) '>> Execution terminates unsuccessfully!' 
      case(13_li)
          write(6,*) 'Precomputed sweep path inconsistent. Execution terminates' 
          write(6,*) '>> Execution terminates unsuccessfully!' 
      case(14_li)
          write(6,*)  'Reflective boundary face could not be determined.' 
          write(6,*) '>> Execution terminates unsuccessfully!' 
      case(15_li)
          write(6,*) 'SLC quadrature only available for orders 1, 2, 3 and 5'
          write(6,*) '>> Execution terminates unsuccessfully!' 
      case(16_li)
          write(6,*) 'A fixed inflow flux value is placed in a boundary face that is not declared fixed inflow.'
          write(6,*) '>> Execution terminates unsuccessfully!' 
      case(17_li)
          write(6,*) "Only lambda=-1 is allowed!"
          write(6,*) '>> Execution terminates unsuccessfully!' 
      case(18_li)
          write(6,*) "Number of Krylov Iterations between restarts must be greater than 0" 
          write(6,*) '>> Execution terminates unsuccessfully!' 
      case(19_li)
          write(6,*) "Maximum number of krylov iterations must be greater than number of",&
                     "iterations between restarts." 
          write(6,*) '>> Execution terminates unsuccessfully!' 
      case(20_li)
          write(6,*) "Method has to be 1 (Outer iteration with lagged upscattering)",&
                     " 2(Flat iteration) or 3(Flat iteration with updated downscattering)." 
          write(6,*) '>> Execution terminates unsuccessfully!' 
      case(21_li)
          write(6,*) 'Counter fail to correctly count number of moments!'
          write(6,*) '>> Execution terminates unsuccessfully!' 
      case(22_li)
          write(6,*) 'Higher Order Spherical Harmonic Not Available!'  
          write(6,*) '>> Execution terminates unsuccessfully!' 
      case(23_li)
          write(6,*) 'Legendre Scattering Order Limited to P7'
          write(6,*) '>> Execution terminates unsuccessfully!' 
      case(24_li)
          write(6,*) 'Reading initial guess file failed' 
          write(6,*) '>> Execution terminates unsuccessfully!' 
      case(25_li)
          write(6,*) "JFNK module requires niter to be equal to 2. This is a coding mistake!"    
          write(6,*) '>> Execution terminates unsuccessfully!' 
      case(26_li)
          write(6,*) 'Method needs to be 1, 2 or 3.'
          write(6,*) '>> Execution terminates unsuccessfully!' 
      case(27_li)
          write(6,*) 'GMRES error'  
          write(6,*) '>> Execution terminates unsuccessfully!' 
      case(28_li)
          write(6,*) 'Only lambda=0 is allowed'  
          write(6,*) '>> Execution terminates unsuccessfully!' 
      case(29_li)
          write(6,*) 'Inflow file only allowed for fixed source problems'
          write(6,*) '>> Execution terminates unsuccessfully!' 
      case(30_li)
          write(6,*) 'VTK source file only allowed for fixed source problems'
          write(6,*) '>> Execution terminates unsuccessfully!' 
      case(31_li)
          write(6,*) 'JFNK option in combination with reflective boundary conditions on opposite faces is not permitted'
          write(6,*) '>> Execution terminates unsuccessfully!' 
      case(32_li)
          write(6,*) 'Call to associated Legendre Polynomials features impossible combination of l and m.'
          write(6,*) 'This is a coding mistake. Contact the developers.'
          write(6,*) '>> Execution terminates unsuccessfully!' 
      case(33_li)
          write(6,*) 'Associated Legendre Polynomial order limited to 5.'
          write(6,*) '>> Execution terminates unsuccessfully!'
      case(34_li)
          write(6,*) 'Density factors were requested but referenced file was not found.'
          write(6,*) '>> Execution terminates unsuccessfully!'
      case default
        if (present(message)) then
          write(6,*) message
        else
          write(6,*) "No Error message given"
        end if
    end select


  ! Cleanup and terminate

    call cleanup
    stop

  end subroutine stop_thor

end module
