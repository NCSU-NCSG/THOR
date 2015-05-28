module wrapup_module
!***********************************************************************
!
! Wraup module contains all subroutines to finish up problem and echo
! output
!
!***********************************************************************

! Use derived-type modules

  use types
  use parameter_types
  use filename_types
  use vector_types
  use cross_section_types
  use geometry_types
  use angle_types
  use global_variables
  use termination_module
 
  implicit none

contains

  subroutine wrapup(flux,keff,unit_number,suffix)
  !**********************************************************************
  !
  ! Subroutine wrapup calls routines that print output into file(s)
  !
  !**********************************************************************

  ! Pass eigenvalue

    real(kind=d_t), intent(in) :: keff

  ! Pass scalar flux 

    real(kind=d_t) :: flux(num_moments_v,namom,num_cells,egmax,niter)

  ! Pass optional arguments
    
    integer(kind=li) :: unit_number
    character(100)   :: suffix 

  ! Declare temporary variables

    integer(kind=li) :: i, eg, l, region
    real(kind=d_t)   :: reac_rates(4,minreg:maxreg,egmax+1) 
    real(kind=d_t)   :: reg_volume(minreg:maxreg) 

  ! Print runtime
    if (rank .eq. 0) then
      write(unit_number,*)
      write(unit_number,*) "--------------------------------------------------------"
      write(unit_number,*) "   Execution Summary   "
      write(unit_number,*) "--------------------------------------------------------"
      if (conv_flag==0) then
         write(unit_number,*) "Warning! Execution finished without satisfying all stopping criteria. Warning!"
      else if(conv_flag==1) then
         write(unit_number,*) "Execution finished successfully. All stopping criteria satisfied."
      end if
      write(unit_number,202) "Runtime (seconds):                                          ", finish-start
      if (problem == 0 .or. (problem == 1 .and. eig_switch == 0 ) ) then
        write(unit_number,101) "Number of outer iterations:                                 ", outer
        write(unit_number,101) "Number of inner iterations:                                 ", tot_nInners
      else
        write(unit_number,101) "Number of newton iterations:                                ", nit
        if (rd_method==1) then
          write(unit_number,101) "Number of inner iterations:                                 ", tot_nInners   
        end if
        write(unit_number,101)   "Total number of krylov iterations:                          ", tot_kit
        write(unit_number,*) 
      end if 
        
      if(problem == 1)then
         write(unit_number,102) "Final eigenvalue:                                          ", keff
      end if
      if(problem==0) then
        write(unit_number,*) "Maximum scalar group flux error (outer iteration):         ", max_outer_error  
        write(unit_number,*) "Maximum scalar flux error by group:"
        do eg=1, egmax
           write(unit_number,103) eg,max_error(eg)
        end do
      else if(problem==1 .and. eig_switch == 1) then
        write(unit_number,*) "Maximum residual by group:"
        do eg=1, egmax
           write(unit_number,103) eg,max_error(eg)
        end do
      else if(problem==1 .and. eig_switch == 0) then
        write(unit_number,*) "Eigenvalue error:                                           ", k_error
        write(unit_number,*) "Maximum fission source error:                               ", f_error   
        write(unit_number,*) "Maximum scalar group flux error:                            ", max_outer_error  
      end if
    end if
  ! Formats
   
    101 FORMAT(1X,A60,I5)
    102 FORMAT(A60,ES25.16) 
    202 FORMAT(1X,A60,ES25.16) 
    103 FORMAT(1X,I5,ES25.16)
  ! Open angular flux file and write volume and face flux in file

    if(vtk_flux_output /= 0 .and. rank .eq. 0)then
       open(unit=10,file=trim(vtk_flux_filename)//trim(suffix),status='unknown',action='write')
 
       write(10,'(a26)') '# vtk DataFile Version 3.0' 
       write(10,'(a72)') jobname
       write(10,'(a5)') 'ASCII'
       write(10,'(a25)') 'DATASET UNSTRUCTURED_GRID'
       write(10,'(a6,1x,i12,1x,a5)') 'POINTS',num_vert,'float'
       
       do i=1, num_vert
          write(10,'(3(1x,es12.5))') vertices(i)%v%x1,vertices(i)%v%x2,&
               vertices(i)%v%x3
       end do
       
       write(10,*)
       write(10,'(a5,1x,i12,1x,i12)') 'CELLS',num_cells,num_cells+&
            4*num_cells
       
       do i=1, num_cells
          write(10,'(5(i12,1x))') 4,cells(i)%R(0)-1,cells(i)%R(1)-1,&
               cells(i)%R(2)-1,cells(i)%R(3)-1
       end do
       
       write(10,*)
       write(10,'(a11,1x,i12)') 'CELL_TYPES',num_cells
       
       do i=1, num_cells
          write(10,'(i12)') 10
       end do
       
       write(10,*)
       write(10,'(a9,1x,i12)') 'CELL_DATA', num_cells
       write(10,'(a5,1x,a16,1x,i12)') 'FIELD','Neutronics_Edits',egmax
       
       l=1
       do eg=1, egmax
          write(10,'(i12,1x,i12,1x,i12,1x,a5)') eg,1,num_cells,'float'
          do i=1, num_cells     
             write(10,'(es12.5)') flux(1,1,i,eg,niter) 
          end do
       end do

       close(10)

    end if

  ! Open flux file and write volume and cell flux in file
    if (rank .eq. 0) then
      open(unit=20,file=trim(flux_filename)//trim(suffix),status='unknown',action='write')

      l=1
      write(20,*) num_cells
      do eg=1, egmax
         do i=1, num_cells     
            write(20,*) cells(i)%volume,flux(1,1,i,eg,niter)
         end do
      end do
      
      close(20)

      write(unit_number,*) "A flux file was been succesfully written!"
    end if
    
  ! Compute region averaged fluxes and reaction rates
    
    if (rank .eq. 0) then
      write(unit_number,*)
      write(unit_number,*) "--------------------------------------------------------"
      write(unit_number,*) "   Region averaged reaction rates  "
      write(unit_number,*) "   (not corrected by )  "
      write(unit_number,*) "--------------------------------------------------------"
      write(unit_number,*) 
      reg_volume=0.0_d_t
      reac_rates=0.0_d_t
      do eg=1,egmax
        do i=1,num_cells
          if(eg.eq.1) reg_volume(cells(i)%reg)=reg_volume(cells(i)%reg)+ cells(i)%volume
          reac_rates(1,cells(i)%reg,eg)=reac_rates(1,cells(i)%reg,eg)  + cells(i)%volume * &
                                        flux(1,1,i,eg,niter)   
          reac_rates(2,cells(i)%reg,eg)=reac_rates(2,cells(i)%reg,eg)  + cells(i)%volume * fiss(reg2mat(cells(i)%reg),eg)%xs * &
                                        flux(1,1,i,eg,niter)   
          reac_rates(3,cells(i)%reg,eg)=reac_rates(3,cells(i)%reg,eg)  + cells(i)%volume *  & 
                                        (sigma_t(reg2mat(cells(i)%reg),eg)%xs - tsigs(reg2mat(cells(i)%reg),eg)%xs )         * &
                                        flux(1,1,i,eg,niter)   
          reac_rates(4,cells(i)%reg,eg)=reac_rates(4,cells(i)%reg,eg)  + cells(i)%volume *  &
                                        nu(reg2mat(cells(i)%reg),eg)%xs*fiss(reg2mat(cells(i)%reg),eg)%xs                    * &
                                        flux(1,1,i,eg,niter)   

        end do
      end do
      do eg=1,egmax
         do region=minreg,maxreg
            reac_rates(1,region,egmax+1)=reac_rates(1,region,egmax+1)+reac_rates(1,region,eg)
            reac_rates(2,region,egmax+1)=reac_rates(2,region,egmax+1)+reac_rates(2,region,eg)
            reac_rates(3,region,egmax+1)=reac_rates(3,region,egmax+1)+reac_rates(3,region,eg)
            reac_rates(4,region,egmax+1)=reac_rates(4,region,egmax+1)+reac_rates(4,region,eg)
         end do
      end do    
      do eg=1,egmax+1
         do region=minreg,maxreg
            reac_rates(1,region,eg)=reac_rates(1,region,eg)/reg_volume(region)
            reac_rates(2,region,eg)=reac_rates(2,region,eg)/reg_volume(region)
            reac_rates(3,region,eg)=reac_rates(3,region,eg)/reg_volume(region)
            reac_rates(4,region,eg)=reac_rates(4,region,eg)/reg_volume(region)
         end do
      end do
      do region=minreg,maxreg
         write(unit_number,501) '-- Region --',region,' Volume= ',reg_volume(region)
         write(unit_number,*) 
         write(unit_number,502) '   Group          Flux       Fission    Absorption      Fiss Src'
         do eg=1,egmax
           write(unit_number,503) eg,reac_rates(1,region,eg),reac_rates(2,region,eg),&
                           reac_rates(3,region,eg),reac_rates(4,region,eg)
         end do
         write(unit_number,504) 'Total   ',reac_rates(1,region,egmax+1),reac_rates(2,region,egmax+1),&
                                 reac_rates(3,region,egmax+1),reac_rates(4,region,egmax+1)   
      end do 
    end if
    
    501 FORMAT(1X,A,I4,A,ES14.6)
    502 FORMAT(1X,A)
    503 FORMAT(1X,I8,4ES14.6)
    504 FORMAT(1X,A8,4ES14.6)
  end subroutine wrapup

  
end module wrapup_module

