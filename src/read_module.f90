module read_module
!***********************************************************************
!
! The module io contains subroutines and variables pertinent to input and
! output.
!
!***********************************************************************

! User derived-type modules

  use mpi
  use strings
  use types
  use parameter_types
  use filename_types
  use vector_types
  use cross_section_types
  use geometry_types
  use angle_types
  use multindex_types
  use global_variables

! Use modules that pertain setting up problem

  use read_cross_section_module
  use readmesh_module
  use read_source_module
  use read_inflow_module
  use quadrature_module
  use check_input
  use termination_module

  implicit none

contains

  subroutine read
  !*********************************************************************
  !
  ! Subroutine read calls subroutines to read input and mesh file
  !
  !*********************************************************************

  ! Declare temporary variable

    integer(kind=li) :: alloc_stat
    integer ::rank,mpi_err, localunit
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)

  ! Set the defaults

    call  set_default

  ! Call read_input to read input file

    call read_input

  ! Check the input

    call check_standard_input

  ! Echo input
    if (rank.eq.0) then
      call  echo_input
    end if

  ! Call read_xs to read cross-section file

    call read_xs

 ! Print cross sections if desired
   if(rank .eq. 0) then
     if(print_xs_flag .eq. 1 ) call print_xs
     if(vtk_mat_output .eq. 1) call plot_material
     if(vtk_reg_output .eq. 1) call plot_region
    end if

 ! If execution is not desired then stop here

   if(execution .eq. 0) call stop_thor(1_li)

 ! Call generate_num_moments and generate_multi_index to create indices

    call generate_num_moments(space_ord,num_moments_v,num_moments_f)

    allocate(index_v(num_moments_v),index_f(num_moments_f),&
         stat=alloc_stat)
    if(alloc_stat /= 0) call stop_thor(2_li)

    call generate_multi_index(space_ord,num_moments_v,num_moments_f,&
         index_v,index_f)

  ! Call read_src to read external source file

    if(problem == 0)then
       call read_src
    end if

  ! Call quad_gen to generate quadrature

    call quad_gen

  ! Call read_finflow to read fixed inflow bc

    if(finflow /= 0 .and. problem .eq. 0)then
       call read_finflow
    end if

  end subroutine read

  subroutine generate_num_moments(spord,numv,numf)
  !*********************************************************************
  !
  ! Subroutine generates the number of total moments
  !
  !*********************************************************************

  ! Pass input parameters

    integer(kind=li), intent(in) :: spord
    integer(kind=li), intent(out) :: numv, numf

  ! Define temporary variables
    integer(kind=li) :: i, j, k, l

    l=1
    if(spord <= 0)then
       do k=0, abs(spord)
          do j=0, abs(spord)
             do i=0, abs(spord)
                if(i+j+k <= abs(spord))then
                   l=l+1
                end if
             end do
          end do
       end do
       numv=l-1
       l=1
       do j=0, abs(spord)
          do i=0, abs(spord)
             if(i+j <= abs(spord))then
                l=l+1
             end if
          end do
       end do
       numf=l-1
    else
       do k=0, spord
          do j=0, spord
             do i=0, spord
                if(i <= spord .or. &
                     j <= spord .or. &
                     k <= spord)then
                   l=l+1
                end if
             end do
          end do
       end do
       numv=l-1
       l=1
       do j=0, 2*spord
          do i=0, 2*spord
             if(i <= 2*spord .or. &
                  j <= 2*spord)then
                l=l+1
             end if
          end do
       end do
       numf=l-1
    end if

  end subroutine generate_num_moments

  subroutine generate_multi_index(spord,numv,numf,indv,indf)
  !*********************************************************************
  !
  ! Subroutine generates the multi-index
  !
  !*********************************************************************

  ! Pass input parameters

    integer(kind=li), intent(in) :: spord, numv,&
         numf

  ! Pass index derived type
    type(indices_v), dimension(numv), intent(inout) :: indv
    type(indices_f), dimension(numf), intent(inout) :: indf

  ! Define temporary variables
    integer(kind=li) :: i, j, k, l

    l=1

    if(spord <= 0)then
       do i=0, abs(spord)
          indv(l)%i1=i
          indv(l)%i2=0
          indv(l)%i3=0
          l=l+1
       end do
       do j=1, abs(spord)
          indv(l)%i1=0
          indv(l)%i2=j
          indv(l)%i3=0
          l=l+1
       end do
       do k=1, abs(spord)
          indv(l)%i1=0
          indv(l)%i2=0
          indv(l)%i3=k
          l=l+1
       end do
       do j=1, abs(spord)
          do i=1, abs(spord)
             if(i+j <= abs(spord))then
                indv(l)%i1=i
                indv(l)%i2=j
                indv(l)%i3=0
                l=l+1
             end if
          end do
       end do
       do k=1, abs(spord)
          do i=1, abs(spord)
             if(i+k <= abs(spord))then
                indv(l)%i1=i
                indv(l)%i2=0
                indv(l)%i3=k
                l=l+1
             end if
          end do
       end do
       do k=1, abs(spord)
          do j=1, abs(spord)
             if(j+k <= abs(spord))then
                indv(l)%i1=0
                indv(l)%i2=j
                indv(l)%i3=k
                l=l+1
             end if
          end do
       end do
       do k=1, abs(spord)
          do j=1, abs(spord)
             do i=1, abs(spord)
                if(i+j+k <= abs(spord))then
                   indv(l)%i1=i
                   indv(l)%i2=j
                   indv(l)%i3=k
                   l=l+1
                end if
             end do
          end do
       end do
       if(l == numv+1)then
       else
         call stop_thor(21_li)
       end if
       l=1
       do i=0, abs(spord)
          indf(l)%i1=i
          indf(l)%i2=0
          l=l+1
       end do
       do j=1, abs(spord)
          indf(l)%i1=0
          indf(l)%i2=j
          l=l+1
       end do
       do j=1, abs(spord)
          do i=1, abs(spord)
             if(i+j <= abs(spord))then
                indf(l)%i1=i
                indf(l)%i2=j
                l=l+1
             end if
          end do
       end do
       if(l == numf+1)then
       else
          call stop_thor(21_li)
       end if
    else
       do i=0, spord
          indv(l)%i1=i
          indv(l)%i2=0
          indv(l)%i3=0
          l=l+1
       end do
       do j=1, spord
          indv(l)%i1=0
          indv(l)%i2=j
          indv(l)%i3=0
          l=l+1
       end do
       do k=1, spord
          indv(l)%i1=0
          indv(l)%i2=0
          indv(l)%i3=k
          l=l+1
       end do
       do j=1, spord
          do i=1, spord
             if(i <= spord .or. &
                  j <= spord)then
                indv(l)%i1=i
                indv(l)%i2=j
                indv(l)%i3=0
                l=l+1
             end if
          end do
       end do
       do k=1, spord
          do i=1, spord
             if(i <= spord .or. &
                  k <= spord)then
                indv(l)%i1=i
                indv(l)%i2=0
                indv(l)%i3=k
                l=l+1
             end if
          end do
       end do
       do k=1, spord
          do j=1, spord
             if(j <= spord .or. &
                  k <= spord)then
                indv(l)%i1=0
                indv(l)%i2=j
                indv(l)%i3=k
                l=l+1
             end if
          end do
       end do
       do k=1, spord
          do j=1, spord
             do i=1, spord
                if(i <= spord .or. &
                     j <= spord .or. &
                     k <= spord)then
                   indv(l)%i1=i
                   indv(l)%i2=j
                   indv(l)%i3=k
                   l=l+1
                end if
             end do
          end do
       end do
       if(l == numv+1)then
       else
          call stop_thor(21_li)
       end if
       l=1
       do i=0, 2*spord
          indf(l)%i1=i
          indf(l)%i2=0
          l=l+1
       end do
       do j=1, 2*spord
          indf(l)%i1=0
          indf(l)%i2=j
          l=l+1
       end do
       do j=1, 2*spord
          do i=1, 2*spord
             if(i <= 2*spord .or. &
                  j <= 2*spord)then
                indf(l)%i1=i
                indf(l)%i2=j
                l=l+1
             end if
          end do
       end do
       if(l == numf+1)then
       else
          call stop_thor(21_li)
       end if
    end if

  end subroutine generate_multi_index

subroutine set_default

  ! Set the input default

    execution = 1
    problem   = 1
    space_ord = 0
    finflow   = 0
    outer_acc = 1
    sweep_tpe = 1
    page_sweep= 0
    page_refl = 0
    page_iflw = 0
    max_outer = 5
    max_inner = 10
    inner_conv= 1.0E-4_d_t
    outer_conv= 1.0E-3_d_t
    k_conv    = 1.0E-4_d_t
    eig_switch= 0
    rd_restart= 25
    rd_max_kit= 250
    rd_method = 2
    inguess_flag =0
    dump_flag    =0
    ipow         =0
    print_conv   =0
    dfact_opt    =0

    source_filename        = "file.src"
    finflow_filename       = "file.bc"
    cross_section_filename = "file.xs"
    mesh_filename          = "file.mesh"
    flux_filename          = "file.flux"
    vtk_flux_filename      = "flux.vtk"
    quad_file              = "file.quad"
    dump_file              = "restart"
    inguess_file           = "inguess"
    vtk_mat_filename       = "mat.vtk"
    vtk_reg_filename       = "reg.vtk"
    vtk_src_filename       = "src.vtk"
    dens_fact_filename     = "density_factors.dat"
    print_xs_flag          = 1
    vtk_flux_output        = 0
    vtk_reg_output         = 0
    vtk_mat_output         = 0
    vtk_src_output         = 0

    quad_ord = 4
    quad_tpe = 1

    egmax = 1
    scatt_ord = 0
    xs_ord    = 0
    upscattering = 1
    multiplying  = 1
    scat_mult_flag=0

    glob_do_cartesian_mesh = .false.

end subroutine

subroutine check_standard_input

  ! Check the input for errors

    ! problem must be 0 or 1
    if(problem > 1 .or. problem < 0)then
      call stop_thor(8_li)
    endif

    ! jfnk parameters
    if(eig_switch==1) then
      if (rd_restart<0) then
        call stop_thor(18_li)
      end if
      if (rd_max_kit<rd_restart) then
        call stop_thor(19_li)
      end if
      if (rd_method<1 .or. rd_method >3) then
        call stop_thor(20_li)
      end if
    end if

    ! in case problem == 1, no inflow file allowed

    if(problem .eq. 1 .and. finflow .eq. 1) then
       call stop_thor(29_li)
    end if

    ! in case problem == 1, no vtk source file allowed

    if(problem .eq. 1 .and. vtk_src_output .eq. 1) then
       call stop_thor(30_li)
    end if


end subroutine

subroutine echo_input

  ! Echo out problem specifications

    integer :: i

    if (rank .eq. 0) then
      write(6,*)
      write(6,*) "--------------------------------------------------------"
      write(6,*) "   Input Summary  "
      write(6,*) "--------------------------------------------------------"

      if(problem == 1)then
         if ( eig_switch == 0 ) then
           write(6,*) "Eigenvalue calculation using PI selected"
         else
           write(6,*) "Eigenvalue calculation using JFNK selected"
           if(rd_method == 1) then
             write(6,*) "Method: F(u) is evaluated using one outer iteration with lagged upscattering."
           else if (rd_method == 2) then
             write(6,*) "Method: F(u) is evaluated using flat iteration scheme."
           else
             write(6,*) "Method: F(u) is evaluated using flat iteration scheme wit updated downscattering."
           end if
         end if
      else if(problem == 0)then
         write(6,*) "External source calculation selected"
      end if

    ! Print rest of input

      write(6,*)   "Problem Title:                                              ", jobname
      write(6,101) "Spatial order:                                              ", space_ord
      if      (sweep_tpe .eq. 1) then
         write(6,*) 'Precomputed mesh sweep is used.'
      end if
      if      (outer_acc .eq. 1 .and. problem .eq. 1 .and. eig_switch .eq. 0) then
         write(6,*) 'Power iterations are not accelerated'
      else if (outer_acc .eq. 2 .and. problem .eq. 1 .and. eig_switch .eq. 0) then
         write(6,*) 'Error mode extrapolation used for accelerating power iterations'
      end if
      write(6,101) "Scattering order:                                           ", scatt_ord
      if(quad_tpe == 1) then
        write(6,101) "Level symmetric quadrature of order:                        ", quad_ord
      else if(quad_tpe == 2) then
        write(6,101) "Square Legendre-Chebychev quadrature of order:              ", quad_ord
      else if(quad_tpe == 3) then
        write(6,101) "Quadrature read from file. #Angles/octant:                  ",nangle
        write(6,*)   "Quadrature file name:                                       ",quad_file
      end if
      write(6,101) "Cross-section order:                                        ", xs_ord
      write(6,101) "Energy groups:                                              ", egmax
      if (problem == 0 .or. (problem == 1 .and. eig_switch == 0 ) ) then
        write(6,101) "Maximum number of outer iterations:                         ", max_outer
        write(6,101) "Maximum number of inner iterations:                         ", max_inner
        write(6,102)   " Inner convergence criteria:                                 ", inner_conv
        if(problem==1 .and. eig_switch==0) write(6,102)   " Eigenvalue convergence criteria:                            ", k_conv
        write(6,102)   " Outer convergence criteria:                                 ", outer_conv
      else
        write(6,101) "Maximum number of newton iterations:                        ", max_outer
        write(6,101) "Maximum number of inner iterations(if used):                ", max_inner
        write(6,102) " Newton convergence criteria:                                  ", outer_conv
        write(6,102) " Inner convergence criteria (if used):                         ", inner_conv
        write(6,101) "Number of krylov iterations between restarts:               ",rd_restart
        write(6,101) "Maximum number of krylov iterations per newton iteration:   ",rd_max_kit
      end if
      if(problem == 0)then
         write(6,*) "File containing external source:                            ", source_filename
      endif

      if(finflow /= 0 .and. problem == 0)then
         write(6,*) "File containing fixed inflow boundary conditions:           ", finflow_filename
      end if

      write(6,*) "File containing cross-sections:                             ", cross_section_filename
      write(6,*) "File containing mesh:                                       ", mesh_filename
      write(6,*) "Flux output file:                                           ", flux_filename

      if(vtk_flux_output /= 0)then
         write(6,*) "VTK-format flux output file:                                ", vtk_flux_filename
      end if
      if(vtk_mat_output /= 0)then
         write(6,*) "VTK-format material output file:                            ", vtk_mat_filename
      end if
      if(vtk_reg_output /= 0)then
         write(6,*) "VTK-format region output file:                              ", vtk_reg_filename
      end if
      if(vtk_src_output /= 0)then
         write(6,*) "VTK-format region output file:                              ", vtk_src_filename
      end if

      if(inguess_flag /= 0)then
         write(6,*) "Initial guess read from file:                               ", inguess_file
      end if

      if(dump_flag /= 0)then
         write(6,*) "Restart file:                                               ", dump_file
      end if
    end if
    ! Formats

    101 FORMAT(1X,A60,I5)
    102 FORMAT(A60,ES12.4)
    202 FORMAT(1X,A60,ES12.4)
    103 FORMAT(1X,I5,ES12.4)
    104 FORMAT(1X,I14,A,I14)

end subroutine

subroutine read_input
!***********************************************************************
!
! This subroutine reads the standard input file
!
!***********************************************************************

 ! local variables

   character(100) :: buffer, fname
   character(100000) :: regmap
   logical :: done
   integer :: i, rank,mpi_err, localunit
   call GET_COMMAND_ARGUMENT(1,fname)
   call MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)
   localunit = rank+100

   open(unit = localunit, file = fname , status = 'old', action = 'read')

 ! read title


    read(localunit,*) jobname

 ! main read loop

    done = .false.
    do while(done .eqv. .false.)

      read(localunit,101,end=999) buffer
      if     ( index( lowercase(buffer) ,'start') > 0 .and. index( lowercase(buffer) ,'problem')>0 ) then
          call read_problem
      else if( index( lowercase(buffer) ,'start') > 0 .and. index( lowercase(buffer) ,'inout')>0   ) then
          call read_inout
      else if( index( lowercase(buffer) ,'start') > 0 .and. index( lowercase(buffer) ,'cross_sections')>0   ) then
          call read_cross_sections
      else if( index( lowercase(buffer) ,'start') > 0 .and. index( lowercase(buffer) ,'quadrature')>0   ) then
          call read_quadrature_field
      else if( index( lowercase(buffer) ,'start') > 0 .and. index( lowercase(buffer) ,'postprocess')>0   ) then
          call read_postprocess_field
      else if( index( lowercase(buffer) ,'start') > 0 .and. index( lowercase(buffer) ,'regionmap')>0   ) then
          call read_regionmap_field(regmap)
      else if( index( lowercase(buffer) ,'end') > 0   .and. index( lowercase(buffer) ,'file')   >0 ) then
         done=.true.
      end if

    end do

999 continue

 ! Call read_mesh to read tetrahedral mesh file

   call read_tetmesh

 ! Read reg2mat

    allocate(reg2mat(minreg:maxreg))
    read(regmap,*) (reg2mat(i),i=minreg,maxreg)

101 FORMAT(A100)

end subroutine

subroutine read_regionmap_field(regmap)

  ! reads the regionmap into array regmap

  ! Arguments

    character(100000) :: regmap

  ! local variables

    character(100) :: line, fname
    integer :: l,lr
    integer :: i, rank,mpi_err, localunit
    call GET_COMMAND_ARGUMENT(1,fname)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)
    localunit = rank+100

  ! read input line by line

    regmap=""
    do while(.true.)
      read(localunit,101) line
      if ( index( lowercase(line) ,'end') > 0 ) then
        return
      else
        l      = len(trim(regmap))
        lr     = len(trim(line))
        regmap(l+1:l+1+lr) = trim(line)
      end if
    end do

101 FORMAT(A100)

end subroutine

subroutine read_quadrature_field

  ! local variables
  integer :: nwords,ntmp,i,nwwords,ios
  character(100) :: buffer, fname
  character(100) :: words(100),wwords(2)
  integer :: rank,mpi_err, localunit
  call GET_COMMAND_ARGUMENT(1,fname)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)
  localunit = rank+100

  ! read loop over inout block
  do while(.true.)
    read(localunit,101) buffer
    if ( index( lowercase(buffer) ,'end') > 0 ) then
       return
     else
       call parse(buffer,";",words,nwords)
       do i=1,nwords
          call parse(words(i),"=",wwords,nwwords)
          if(nwwords .eq. 2) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'qdtype'
            if( trim(lowercase(wwords(1))) .eq. 'qdtype' ) then
               wwords(2)=trim(lowercase(wwords(2)))
               if      ( wwords(2) .eq. 'levelsym') then
                  quad_tpe=1
               else if ( wwords(2) .eq. 'legcheb') then
                  quad_tpe=2
               else if ( wwords(2) .eq. 'fromfile') then
                  quad_tpe=3
               else
                  write(6,*) 'Error. This is not a valid quadrature type -- ',wwords(2),' --'
                  write(6,*) 'Execution will terminate.'
                  stop
               end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'qdorder'
            else if( trim(lowercase(wwords(1))) .eq. 'qdorder' ) then
              wwords(2)=trim(lowercase(wwords(2)))
              read(wwords(2),'(i10)',iostat=ios) quad_ord
              if(ios.ne.0 ) then
                write(6,*) 'Invalid quadrature order -- ',trim(wwords(2)),' --'
                stop
              end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > default if keyword is unknown
            else
               write(6,*) 'Unknown keyword in quadrature specification -- ',trim(wwords(1)),' --'
               write(6,*) 'Execution will terminate.'
               stop
            end if
          else
             write(6,*) 'Error while reading cross section specification'
             write(6,*) 'Do not understand entry: ',trim(words(i))
             write(6,*) 'Execution will terminate.'
             stop
          end if
       end do
     end if
  end do

101 FORMAT(A100)

end subroutine

subroutine read_postprocess_field

  ! local variables
  integer :: nwords, ntmp, i, nwwords, ios, nwwwords
  character(100) :: buffer, fname
  character(100) :: words(100), wwords(2), wwwords(100)
  integer :: rank, mpi_err, localunit
  call GET_COMMAND_ARGUMENT(1,fname)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)
  localunit = rank+100

  ! read loop over inout block
  do while(.true.)
    read(localunit,101) buffer
    if ( index( lowercase(buffer) ,'end') > 0 ) then
       return
    else
      call parse(buffer,";",words,nwords)
      do i=1,nwords
        call parse(words(i),"=",wwords,nwwords)
        if(nwwords .eq. 2) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'qdtype'
          if( trim(lowercase(wwords(1))) .eq. 'cartesian_map' ) then
             wwords(2)=trim(lowercase(wwords(2)))
             glob_do_cartesian_mesh = .true.
             ! wwords must be an array with of length 6
             call parse(wwords(2), " ", wwwords, nwwwords)
             if (nwwwords .ne. 9) then
               write(6,*) 'Following cartesian map nine entries are required; Found: ',&
                           trim(wwords(2)),' has ', nwwwords, ' entries.'
             end if
             glob_cmap_min_x = string_to_real(wwwords(1), 'Conversion to cartesian map xmin failed')
             glob_cmap_max_x = string_to_real(wwwords(1), 'Conversion to cartesian map xmax failed')
             glob_cmap_nx = string_to_int(wwwords(3), 'Conversion to cartesian map nx failed', 1)
             glob_cmap_min_y = string_to_real(wwwords(4), 'Conversion to cartesian map ymin failed')
             glob_cmap_max_y = string_to_real(wwwords(5), 'Conversion to cartesian map ymax failed')
             glob_cmap_ny = string_to_int(wwwords(6), 'Conversion to cartesian map ny failed', 1)
             glob_cmap_min_z = string_to_real(wwwords(7), 'Conversion to cartesian map zmin failed')
             glob_cmap_max_z = string_to_real(wwwords(8), 'Conversion to cartesian map zmax failed')
             glob_cmap_nz = string_to_int(wwwords(9), 'Conversion to cartesian map nz failed', 1)
          else
             write(6,*) 'Unknown keyword in postprocess specification -- ',trim(wwords(1)),' --'
             write(6,*) 'Execution will terminate.'
             stop
          end if
        else
          write(6,*) 'Error while reading postprocess specification'
          write(6,*) 'Do not understand entry: ',trim(words(i))
          write(6,*) 'Execution will terminate.'
          stop
        end if
      end do
    end if
  end do

101 FORMAT(A100)

end subroutine

subroutine read_cross_sections

  ! local variables
  integer :: nwords,ntmp,nwwords,ios
  character(100) :: buffer, fname
  character(100) :: words(100),wwords(2)
  integer :: i, rank,mpi_err, localunit
  call GET_COMMAND_ARGUMENT(1,fname)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)
  localunit = rank+100

  ! read loop over inout block
  do while(.true.)
    read(localunit,101) buffer
    if ( index( lowercase(buffer) ,'end') > 0 ) then
       return
     else
       call parse(buffer,";",words,nwords)
       do i=1,nwords
          call parse(words(i),"=",wwords,nwwords)
          if(nwwords .eq. 2) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'ngroups'
            if( trim(lowercase(wwords(1))) .eq. 'ngroups' ) then
              wwords(2)=trim(lowercase(wwords(2)))
               read(wwords(2),'(i10)',iostat=ios) egmax
               if(ios.ne.0 .or. egmax<1) then
                  write(6,*) 'Invalid number of energy groups in cross section specification -- ',trim(wwords(2)),' --'
                  stop
               end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'pnorder'
            else if( trim(lowercase(wwords(1))) .eq. 'pnorder' ) then
              wwords(2)=trim(lowercase(wwords(2)))
               read(wwords(2),'(i10)',iostat=ios) scatt_ord
               if(ios.ne.0 .or. scatt_ord < 0) then
                  write(6,*) 'Invalid scattering expansion in cross section specification -- ',trim(wwords(2)),' --'
                  stop
               end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword pnread'
            else if( trim(lowercase(wwords(1))) .eq. 'pnread' ) then
              wwords(2)=trim(lowercase(wwords(2)))
               read(wwords(2),'(i10)',iostat=ios) xs_ord
               if(ios.ne.0 .or. xs_ord < 0) then
                  write(6,*) 'Invalid cross section expansion in cross section specification -- ',trim(wwords(2)),' --'
                  stop
               end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'upscattering'
            else if( trim(lowercase(wwords(1))) .eq. 'upscattering' ) then
               wwords(2)=trim(lowercase(wwords(2)))
               if      ( wwords(2) .eq. 'yes') then
                  upscattering=1
               else if ( wwords(2) .eq. 'no') then
                  upscattering=0
               else
                  write(6,*) 'Error. This is not a valid upscattering flag -- ',wwords(2),' --'
                  write(6,*) 'Execution will terminate.'
                  stop
               end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword !'multiplying'
            else if( trim(lowercase(wwords(1))) .eq. 'multiplying' ) then
               wwords(2)=trim(lowercase(wwords(2)))
               if      ( wwords(2) .eq. 'yes') then
                  multiplying=1
               else if ( wwords(2) .eq. 'no') then
                  multiplying=0
               else
                  write(6,*) 'Error. This is not a valid multiplying flag -- ',wwords(2),' --'
                  write(6,*) 'Execution will terminate.'
                  stop
               end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'scatt_mult_included'
            else if( trim(lowercase(wwords(1))) .eq. 'scatt_mult_included' ) then
               wwords(2)=trim(lowercase(wwords(2)))
               if      ( wwords(2) .eq. 'yes') then
                  scat_mult_flag=1
               else if ( wwords(2) .eq. 'no') then
                  scat_mult_flag=0
               else
                  write(6,*) 'Error. This is not a valid scattering multiplier flag -- ',wwords(2),' --'
                  write(6,*) 'Execution will terminate.'
                  stop
               end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > default if keyword is unknown
            else
               write(6,*) 'Unknown keyword in cross section specification -- ',trim(wwords(1)),' --'
               write(6,*) 'Execution will terminate.'
               stop
            end if
          else
             write(6,*) 'Error while reading cross section specification'
             write(6,*) 'Do not understand entry: ',trim(words(i))
             write(6,*) 'Execution will terminate.'
             stop
          end if
       end do
     end if
  end do

101 FORMAT(A100)

end subroutine

subroutine read_inout

  ! local variables
  integer :: nwords,ntmp,i,nwwords,ios
  character(100) :: buffer, fname
  character(100) :: words(100),wwords(2)
  integer :: rank,mpi_err, localunit
  call GET_COMMAND_ARGUMENT(1,fname)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)
  localunit = rank+100

  ! read loop over inout block
  do while(.true.)
    read(localunit,101) buffer
    if ( index( lowercase(buffer) ,'end') > 0 ) then
       return
     else
       call parse(buffer,";",words,nwords)
       do i=1,nwords
          call parse(words(i),"=",wwords,nwwords)
          if(nwwords .eq. 2) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'meshi_file'
            if( trim(lowercase(wwords(1))) .eq. 'mesh_file' ) then
               mesh_filename=trim(wwords(2))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'inflow_file'
            else if( trim(lowercase(wwords(1))) .eq. 'inflow_file' ) then
               finflow_filename=trim(wwords(2))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'source_file'
            else if( trim(lowercase(wwords(1))) .eq. 'source_file' ) then
               source_filename=trim(wwords(2))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'flux_file'
            else if( trim(lowercase(wwords(1))) .eq. 'flux_file' ) then
               flux_filename=trim(wwords(2))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'xs_file'
            else if( trim(lowercase(wwords(1))) .eq. 'xs_file' ) then
               cross_section_filename=trim(wwords(2))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'quad_file'
            else if( trim(lowercase(wwords(1))) .eq. 'density_factor_file' ) then
               dens_fact_filename=trim(wwords(2))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'quad_file'
            else if( trim(lowercase(wwords(1))) .eq. 'quad_file' ) then
               quad_file=trim(wwords(2))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'vtk'
            else if( trim(lowercase(wwords(1))) .eq. 'vtk' ) then
               wwords(2)=trim(lowercase(wwords(2)))
               if( index(wwords(2),'flux') > 0 ) vtk_flux_output=1
               if( index(wwords(2),'mat'  ) > 0 ) vtk_mat_output=1
               if( index(wwords(2),'reg') > 0 ) vtk_reg_output=1
               if( index(wwords(2),'src') > 0 ) vtk_src_output=1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'vtk_flux_file'
            else if( trim(lowercase(wwords(1))) .eq. 'vtk_flux_file' ) then
               vtk_flux_filename=trim(wwords(2))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'vtk_mat_file'
            else if( trim(lowercase(wwords(1))) .eq. 'vtk_mat_file' ) then
               vtk_mat_filename=trim(wwords(2))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'vtk_reg_file'
            else if( trim(lowercase(wwords(1))) .eq. 'vtk_reg_file' ) then
               vtk_reg_filename=trim(wwords(2))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'vtk_src_file'
            else if( trim(lowercase(wwords(1))) .eq. 'vtk_src_file' ) then
               vtk_src_filename=trim(wwords(2))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'restart_file'
            else if( trim(lowercase(wwords(1))) .eq. 'restart_file' ) then
               dump_file=trim(wwords(2))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'inguess_file'
            else if( trim(lowercase(wwords(1))) .eq. 'inguess_file' ) then
               inguess_file=trim(wwords(2))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'print_xs'
            else if( trim(lowercase(wwords(1))) .eq. 'print_xs' ) then
               wwords(2)=trim(lowercase(wwords(2)))
               if      ( wwords(2) .eq. 'yes') then
                  print_xs_flag=1
               else if ( wwords(2) .eq. 'no') then
                  print_xs_flag=0
               else
                  write(6,*) 'Error. This is not a valid cross section print option (yes/no) -- ',wwords(2),' --'
                  write(6,*) 'Execution will terminate.'
                  stop
               end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > default if keyword is unknown
            else
               write(6,*) 'Unknown keyword in inout field -- ',trim(wwords(1)),' --'
               write(6,*) 'Execution will terminate.'
               stop
            end if
          else
             write(6,*) 'Error while reading inout specification'
             write(6,*) 'Do not understand entry: ',trim(words(i))
             write(6,*) 'Execution will terminate.'
             stop
          end if
       end do
     end if
  end do

101 FORMAT(A100)

end subroutine

subroutine read_problem

  ! local variables
  integer :: nwords,ntmp,i,nwwords,ios
  character(100) :: buffer, fname
  character(100) :: words(100),wwords(2)
  integer :: rank,mpi_err, localunit
  call GET_COMMAND_ARGUMENT(1,fname)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)
  localunit = rank+100

  ! read loop over problem block
  do while(.true.)
     read(localunit,101) buffer
     if ( index( lowercase(buffer) ,'end') > 0 ) then
       return
     else
       call parse(buffer,";",words,nwords)
       do i=1,nwords
          call parse(words(i),"=",wwords,nwwords)
          if(nwwords .eq. 2) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'type'
            if( trim(lowercase(wwords(1))) .eq. 'type' ) then
               wwords(2)=trim(lowercase(wwords(2)))
               if      ( wwords(2) .eq. 'keig') then
                  problem=1
               else if ( wwords(2) .eq. 'fsrc') then
                  problem=0
               else
                  write(6,*) 'Error. This is not a valid problem type -- ',wwords(2),' --'
                  write(6,*) 'Execution will terminate.'
                  stop
               end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'keigsolver'
            else if( trim(lowercase(wwords(1))) .eq. 'keigsolver' ) then
               wwords(2)=trim(lowercase(wwords(2)))
               if      ( wwords(2) .eq. 'pi') then
                  eig_switch=0
               else if ( wwords(2) .eq. 'jfnk') then
                  eig_switch=1
               else
                  write(6,*) 'Error. This is not a valid eigenvalue solver -- ',wwords(2),' --'
                  write(6,*) 'Execution will terminate.'
                  stop
               end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'lambda'
            else if( trim(lowercase(wwords(1))) .eq. 'lambda' ) then
               wwords(2)=trim(lowercase(wwords(2)))
               read(wwords(2),'(i10)',iostat=ios) space_ord
               if(ios.ne.0) then
                  write(6,*) 'Invalid spatial order in problem specification -- ',trim(wwords(2)),' --'
                  stop
               end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'inflow'
            else if( trim(lowercase(wwords(1))) .eq. 'inflow' ) then
               wwords(2)=trim(lowercase(wwords(2)))
               if      ( wwords(2) .eq. 'yes') then
                  finflow=1
               else if ( wwords(2) .eq. 'no') then
                  finflow=0
               else
                  write(6,*) 'Error. This is not a valid inflow flag -- ',wwords(2),' --'
                  write(6,*) 'Execution will terminate.'
                  stop
               end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'PIacc'
            else if( trim(lowercase(wwords(1))) .eq. 'piacc' ) then
               wwords(2)=trim(lowercase(wwords(2)))
               if      ( wwords(2) .eq. 'errmode') then
                  outer_acc=2
               else if ( wwords(2) .eq. 'none') then
                  outer_acc=1
               else
                  write(6,*) 'Error. This is not a valid acceleration option for PI -- ',wwords(2),' --'
                  write(6,*) 'Execution will terminate.'
                  stop
               end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'sweep'
            else if( trim(lowercase(wwords(1))) .eq. 'sweep' ) then
               wwords(2)=trim(lowercase(wwords(2)))
               if      ( wwords(2) .eq. 'precomp') then
                  sweep_tpe=1
               else
                  write(6,*) 'Error. This is not a valid mesh sweep type -- ',wwords(2),' --'
                  write(6,*) 'Execution will terminate.'
                  stop
               end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'sweep_page'
            else if( trim(lowercase(wwords(1))) .eq. 'page_sweep' ) then
               wwords(2)=trim(lowercase(wwords(2)))
               if      ( wwords(2) .eq. 'yes') then
                  page_sweep=1
               else if ( wwords(2) .eq. 'no') then
                  page_sweep=0
               else
                  write(6,*) 'Error. This is not a valid sweep page option -- ',wwords(2),' --'
                  write(6,*) 'Execution will terminate.'
                  stop
               end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'page_refl'
            else if( trim(lowercase(wwords(1))) .eq. 'page_refl' ) then
               wwords(2)=trim(lowercase(wwords(2)))
               if      ( wwords(2) .eq. 'page') then
                  page_refl=1
               else if ( wwords(2) .eq. 'save') then
                  page_refl=0
               else if ( wwords(2) .eq. 'inner') then
                  page_refl=2
               else
                  write(6,*) 'Error. This is not a valid page reflective BC option -- ',wwords(2),' --'
                  write(6,*) 'Execution will terminate.'
                  stop
               end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'page_iflw'
            else if( trim(lowercase(wwords(1))) .eq. 'page_inflow' ) then
               wwords(2)=trim(lowercase(wwords(2)))
               if      ( wwords(2) .eq. 'bygroup') then
                  page_iflw=1
               else if ( wwords(2) .eq. 'all') then
                  page_iflw=0
               else
                  write(6,*) 'Error. This is not a valid page inflow option -- ',wwords(2),' --'
                  write(6,*) 'Execution will terminate.'
                  stop
               end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'kconv'
            else if( trim(lowercase(wwords(1))) .eq. 'kconv' ) then
               wwords(2)=trim(lowercase(wwords(2)))
               read(wwords(2),*,iostat=ios) k_conv
               if(ios.ne.0) then
                  write(6,*) 'Invalid stopping criterion for keff in problem specification -- ',trim(wwords(2)),' --'
                  stop
               end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'innerconv'
            else if( trim(lowercase(wwords(1))) .eq. 'innerconv' ) then
               wwords(2)=trim(lowercase(wwords(2)))
               read(wwords(2),*,iostat=ios) inner_conv
               if(ios.ne.0) then
                  write(6,*) 'Invalid stopping criterion for inner iterations in problem specification -- ',trim(wwords(2)),' --'
                  stop
               end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'outerconv'
            else if( trim(lowercase(wwords(1))) .eq. 'outerconv' ) then
               wwords(2)=trim(lowercase(wwords(2)))
               read(wwords(2),*,iostat=ios) outer_conv
               if(ios.ne.0) then
                  write(6,*) 'Invalid stopping criterion for outer iterations in problem specification -- ',trim(wwords(2)),' --'
                  stop
               end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'maxinner'
            else if( trim(lowercase(wwords(1))) .eq. 'maxinner' ) then
              wwords(2)=trim(lowercase(wwords(2)))
               read(wwords(2),'(i10)',iostat=ios) max_inner
               if(ios.ne.0 .or. max_inner<1 ) then
                  write(6,*) 'Invalid maximum number of inner iteration in problem specification -- ',trim(wwords(2)),' --'
                  stop
               end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'maxouter'
            else if( trim(lowercase(wwords(1))) .eq. 'maxouter' ) then
              wwords(2)=trim(lowercase(wwords(2)))
               read(wwords(2),'(i10)',iostat=ios) max_outer
               if(ios.ne.0 .or. max_outer<1 ) then
                  write(6,*) 'Invalid maximum number of outer iteration in problem specification -- ',trim(wwords(2)),' --'
                  stop
               end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'jfnk_krsze'
            else if( trim(lowercase(wwords(1))) .eq. 'jfnk_krsze' ) then
              wwords(2)=trim(lowercase(wwords(2)))
               read(wwords(2),'(i10)',iostat=ios) rd_restart
               if(ios.ne.0 .or. rd_restart<1 ) then
                  write(6,*) 'Invalid Krylov subspace size for JFNK in problem specification -- ',trim(wwords(2)),' --'
                  stop
               end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'jfnk_maxkr'
            else if( trim(lowercase(wwords(1))) .eq. 'jfnk_maxkr' ) then
              wwords(2)=trim(lowercase(wwords(2)))
               read(wwords(2),'(i10)',iostat=ios) rd_max_kit
               if(ios.ne.0 .or. rd_max_kit < 1 ) then
                  write(6,*) 'Invalid maximum number of Krylov iterations for JFNK in problem specification -- '&
                             ,trim(wwords(2)),' --'
                  stop
               end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'jfnk_method'
            else if( trim(lowercase(wwords(1))) .eq. 'jfnk_method' ) then
               wwords(2)=trim(lowercase(wwords(2)))
               if      ( wwords(2) .eq. 'outer') then
                  rd_method=1
               else if ( wwords(2) .eq. 'flat') then
                  rd_method=2
               else if ( wwords(2) .eq. 'flat_wds') then
                  rd_method=3
               else
                  write(6,*) 'Error. This is not a valid jfnk solution method -- ',wwords(2),' --'
                  write(6,*) 'Execution will terminate.'
                  stop
               end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'initial_guess'
            else if( trim(lowercase(wwords(1))) .eq. 'initial_guess' ) then
               wwords(2)=trim(lowercase(wwords(2)))
               if      ( wwords(2) .eq. 'yes') then
                  inguess_flag=1
               else if ( wwords(2) .eq. 'no') then
                  inguess_flag=0
               else
                  write(6,*) 'Error. This is not a valid initial guess option (yes/no) -- ',wwords(2),' --'
                  write(6,*) 'Execution will terminate.'
                  stop
               end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'save_restart'
            else if( trim(lowercase(wwords(1))) .eq. 'save_restart' ) then
               wwords(2)=trim(lowercase(wwords(2)))
               if      ( wwords(2) .eq. 'yes') then
                  dump_flag=1
               else if ( wwords(2) .eq. 'no') then
                  dump_flag=0
               else
                  write(6,*) 'Error. This is not a valid restart file option (yes/no) -- ',wwords(2),' --'
                  write(6,*) 'Execution will terminate.'
                  stop
               end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword !'ipiter'
            else if( trim(lowercase(wwords(1))) .eq. 'ipiter' ) then
              wwords(2)=trim(lowercase(wwords(2)))
               read(wwords(2),'(i10)',iostat=ios) ipow
               if(ios.ne.0 ) then
                  write(6,*) 'Invalid number of initial power iterations -- '&
                             ,trim(wwords(2)),' --'
                  stop
               end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'print_conv'
             else if( trim(lowercase(wwords(1))) .eq. 'print_conv' ) then
               wwords(2)=trim(lowercase(wwords(2)))
               if      ( wwords(2) .eq. 'yes') then
                  print_conv=1
               else if ( wwords(2) .eq. 'no') then
                  print_conv=0
               else
                  write(6,*) 'Error. This is not a valid execution option (yes/no) -- ',wwords(2),' --'
                  write(6,*) 'Execution will terminate.'
                  stop
               end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'density factor option'
             else if( trim(lowercase(wwords(1))) .eq. 'density_factor' ) then
               wwords(2)=trim(lowercase(wwords(2)))
               if      ( wwords(2) .eq. 'no') then
                  dfact_opt = 0
               else
                 if      ( wwords(2) .eq. 'byvolume') then
                   dfact_opt = 1
                 else if ( wwords(2) .eq. 'fromfile') then
                   dfact_opt = 2
                 else
                   write(6,*) 'Error. This is not a valid density factor option &
                              (no/byvolume/fromfile) -- ',wwords(2),' --'
                   write(6,*) 'Execution will terminate.'
                   stop
                 end if
               end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > keyword 'execution'
            else if( trim(lowercase(wwords(1))) .eq. 'execution' ) then
               wwords(2)=trim(lowercase(wwords(2)))
               if      ( wwords(2) .eq. 'yes') then
                  execution=1
               else if ( wwords(2) .eq. 'no') then
                  execution=0
               else
                  write(6,*) 'Error. This is not a valid execution option (yes/no) -- ',wwords(2),' --'
                  write(6,*) 'Execution will terminate.'
                  stop
               end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! > default if keyword is unknown
            else
               write(6,*) 'Unknown keyword in problem specification -- ',trim(wwords(1)),' --'
               write(6,*) 'Execution will terminate.'
               stop
            end if
          else
            write(6,*) 'Error while reading problem specification'
            write(6,*) 'Do not understand entry: ',trim(words(i))
            write(6,*) 'Execution will terminate.'
            stop
          end if
       end do
     end if

  end do

101 FORMAT(A100)
end subroutine

real(kind=d_t) function string_to_real(string, msg)
  character(100), intent(in) :: string
  character(100), intent(in) :: msg

  integer(kind=li) :: ios

  read(string, *, iostat=ios) string_to_real
  if(ios .ne. 0) then
     write(6,*) trim(msg), trim(string)
     stop
  end if
end function

integer(kind=li) function string_to_int(string, msg, min_int)
  character(100), intent(in) :: string
  character(100), intent(in) :: msg
  integer(kind=li), optional, intent(inout) :: min_int

  integer(kind=li) :: ios

  if (.not. present(min_int)) min_int = -glob_max_int
  read(string, '(i10)', iostat=ios) string_to_int
  if(ios .ne. 0 .or. string_to_int < min_int) then
     write(6,*) trim(msg), trim(string)
     stop
  end if
end function

end module
