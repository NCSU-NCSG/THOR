!***********************************************************************
!
! The module contains all legacy subroutines for old yaml_input versions. These
! subroutines should generally not be changed to garauntee backwards
! compatibility.
!
!***********************************************************************
MODULE read_module_current

  USE mpi
  USE stringmod
  USE global_variables

  IMPLICIT NONE
  PRIVATE
  !
  ! List of public members
  PUBLIC :: inputfile_read

!> The maximum length of a cardname
INTEGER,PARAMETER :: MAX_CARDNAME_LEN=32
!> The number of cards we have
INTEGER,PARAMETER :: num_cards=48
!> The maximum length of a line in the input file
INTEGER,PARAMETER :: ll_max=200
!(also need to change in interface if changed)
!> The maximum number of params on a given line or inputs for a param
INTEGER,PARAMETER :: lp_max=100
!(also need to change in interface if changed)

TYPE :: cardType
  !character name of the card
  CHARACTER(MAX_CARDNAME_LEN) :: cname
  !logical to tell if card has already appeared
  LOGICAL :: used=.FALSE.
  !readin procedure, unique for each card
  PROCEDURE(prototype_wordarg),POINTER :: getcard => NULL()
ENDTYPE cardType

TYPE(cardType) :: cards(num_cards)

!> Simple abstract interface for an object method subroutine with word array argument
ABSTRACT INTERFACE
  SUBROUTINE prototype_wordarg(thisCard,twords)
    IMPORT :: cardType
    CLASS(cardType),INTENT(INOUT) :: thisCard
    CHARACTER(200),INTENT(IN) :: twords(100)
  ENDSUBROUTINE prototype_wordarg
ENDINTERFACE

CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE inputfile_read(local_unit)
    INTEGER, INTENT(IN) :: local_unit
    CHARACTER(ll_max) :: tchar
    CHARACTER(ll_max) :: words(lp_max),wwords(lp_max)
    INTEGER :: ios,rank,mpi_err,nwords,nwwords,iw,iww,ic
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)
    IF(local_unit .NE. 100+rank)STOP 'rank mismatch in input file reading'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !setup the card objects, add one to here if you need a new one, make sure to change the number
    !of cards param as well
    cards(1)%cname='type'
    cards(1)%getcard => get_type
    cards(2)%cname='keigsolver'
    cards(2)%getcard => get_keigsolver
    cards(3)%cname='lambda'
    cards(3)%getcard => get_lambda
    cards(4)%cname='inflow'
    cards(4)%getcard => get_inflow
    cards(5)%cname='piacc'
    cards(5)%getcard => get_piacc
    cards(6)%cname='page_sweep'
    cards(6)%getcard => get_page_sweep
    cards(7)%cname='page_refl'
    cards(7)%getcard => get_page_refl
    cards(8)%cname='page_iflw'
    cards(8)%getcard => get_page_iflw
    cards(9)%cname='kconv'
    cards(9)%getcard => get_kconv
    cards(10)%cname='innerconv'
    cards(10)%getcard => get_innerconv
    cards(11)%cname='outerconv'
    cards(11)%getcard => get_outerconv
    cards(12)%cname='maxinner'
    cards(12)%getcard => get_maxinner
    cards(13)%cname='maxouter'
    cards(13)%getcard => get_maxouter
    cards(14)%cname='jfnk_krsze'
    cards(14)%getcard => get_jfnk_krsze
    cards(15)%cname='jfnk_maxkr'
    cards(15)%getcard => get_jfnk_maxkr
    cards(16)%cname='jfnk_method'
    cards(16)%getcard => get_jfnk_method
    cards(17)%cname='initial_guess'
    cards(17)%getcard => get_initial_guess
    cards(18)%cname='save_restart'
    cards(18)%getcard => get_save_restart
    cards(19)%cname='ipiter'
    cards(19)%getcard => get_ipiter
    cards(20)%cname='print_conv'
    cards(20)%getcard => get_print_conv
    cards(21)%cname='density_factor'
    cards(21)%getcard => get_density_factor
    cards(22)%cname='execution'
    cards(22)%getcard => get_execution
    cards(23)%cname='mesh_file'
    cards(23)%getcard => get_mesh_file
    cards(24)%cname='inflow_file'
    cards(24)%getcard => get_inflow_file
    cards(25)%cname='source_file'
    cards(25)%getcard => get_source_file
    cards(26)%cname='flux_file'
    cards(26)%getcard => get_flux_file
    cards(27)%cname='xs_file'
    cards(27)%getcard => get_xs_file
    cards(28)%cname='density_factor_file'
    cards(28)%getcard => get_density_factor_file
    cards(29)%cname='quad_file'
    cards(29)%getcard => get_quad_file
    cards(30)%cname='vtk'
    cards(30)%getcard => get_vtk
    cards(31)%cname='vtk_flux_file'
    cards(31)%getcard => get_vtk_flux_file
    cards(32)%cname='vtk_mat_file'
    cards(32)%getcard => get_vtk_mat_file
    cards(33)%cname='vtk_reg_file'
    cards(33)%getcard => get_vtk_reg_file
    cards(34)%cname='vtk_src_file'
    cards(34)%getcard => get_vtk_src_file
    cards(35)%cname='restart_file'
    cards(35)%getcard => get_restart_file
    cards(36)%cname='inguess_file'
    cards(36)%getcard => get_inguess_file
    cards(37)%cname='cartesian_map_file'
    cards(37)%getcard => get_cartesian_map_file
    cards(38)%cname='print_xs'
    cards(38)%getcard => get_print_xs
    cards(39)%cname='ngroups'
    cards(39)%getcard => get_ngroups
    cards(40)%cname='pnorder'
    cards(40)%getcard => get_pnorder
    cards(41)%cname='pnread'
    cards(41)%getcard => get_pnread
    cards(42)%cname='upscattering'
    cards(42)%getcard => get_upscattering
    cards(43)%cname='multiplying'
    cards(43)%getcard => get_multiplying
    cards(44)%cname='scatt_mult_included'
    cards(44)%getcard => get_scatt_mult_included
    cards(45)%cname='qdtype'
    cards(45)%getcard => get_qdtype
    cards(46)%cname='qdorder'
    cards(46)%getcard => get_qdorder
    cards(47)%cname='cartesian_map'
    cards(47)%getcard => get_cartesian_map
    cards(48)%cname='point_value_locations'
    cards(48)%getcard => get_point_value_locations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    DO
      READ(local_unit,'(A200)',IOSTAT=ios) tchar
      tchar=TRIM(ADJUSTL(tchar))
      !ignore blank and commented lines
      IF(tchar .NE. '' .AND. tchar(1:1) .NE. '!')THEN
        CALL parse(tchar,";",words,nwords)
        !loop over inputs on a line
        DO iw=1,nwords
          words(iw)=TRIM(ADJUSTL(words(iw)))
          CALL parse(words(iw)," ",wwords,nwwords)
          wwords(1)=TRIM(ADJUSTL(lowercase(wwords(1))))
          !find the card it belongs to
          DO ic=1,num_cards
            IF(wwords(1) .EQ. cards(ic)%cname)THEN
              IF(cards(ic)%used)STOP 'duplicate params'
              CALL cards(ic)%getcard(wwords)
              cards(ic)%used=.TRUE.
              EXIT
            ENDIF
          ENDDO
          IF(ic .GE. num_cards+1)THEN
            IF(rank .EQ. 0)WRITE(*,'(2A)')'STOPPING: bad input, unrecognized card: ',wwords(1)
            STOP
          ENDIF
        ENDDO
      ENDIF
      IF(ios .NE. 0)EXIT
    ENDDO
    STOP 'inputfile_read not yet complete'
  END SUBROUTINE inputfile_read

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_type(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(IN) :: wwords(lp_max)
    CHARACTER(ll_max) :: tchar1
    tchar1=TRIM(lowercase(wwords(2)))
    IF(tchar1 .EQ. 'keig') THEN
      problem=1
    ELSEIF(tchar1 .EQ. 'fsrc') THEN
      problem=0
    ELSE
      WRITE(6,*) 'Error. This is not a valid problem type -- ',wwords(2),' --'
      WRITE(6,*) 'Execution will terminate.'
      STOP
    ENDIF
  END SUBROUTINE get_type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_keigsolver(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(IN) :: wwords(lp_max)
    CHARACTER(ll_max) :: tchar1
    tchar1=TRIM(lowercase(wwords(2)))
    IF(tchar1 .EQ. 'pi') THEN
      eig_switch=0
    ELSEIF(tchar1 .EQ. 'jfnk') THEN
      eig_switch=1
    ELSE
      WRITE(6,*) 'Error. This is not a valid eigenvalue solver -- ',wwords(2),' --'
      WRITE(6,*) 'Execution will terminate.'
      STOP
    ENDIF
  END SUBROUTINE get_keigsolver

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_lambda(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(IN) :: wwords(lp_max)
    INTEGER :: ios
    CHARACTER(ll_max) :: tchar1
    tchar1=TRIM(lowercase(wwords(2)))
    READ(tchar1,'(i10)',iostat=ios) space_ord
    IF(ios.NE.0) THEN
      WRITE(6,*) 'Invalid spatial order in problem specification -- ',TRIM(wwords(2)),' --'
      STOP
    ENDIF
  END SUBROUTINE get_lambda

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_inflow(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(IN) :: wwords(lp_max)
    CHARACTER(ll_max) :: tchar1
    tchar1=TRIM(lowercase(wwords(2)))
    IF(tchar1.EQ. 'yes') THEN
      finflow=1
    ELSEIF(tchar1 .EQ. 'no') THEN
      finflow=0
    ELSE
      WRITE(6,*) 'Error. This is not a valid inflow flag -- ',wwords(2),' --'
      WRITE(6,*) 'Execution will terminate.'
      STOP
    ENDIF
  END SUBROUTINE get_inflow

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_piacc(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(IN) :: wwords(lp_max)
    CHARACTER(ll_max) :: tchar1
    tchar1=TRIM(lowercase(wwords(2)))
    IF(tchar1 .EQ. 'errmode') THEN
      outer_acc=2
    ELSEIF(tchar1 .EQ. 'none') THEN
      outer_acc=1
    ELSE
      WRITE(6,*) 'Error. This is not a valid acceleration option for PI -- ',wwords(2),' --'
      WRITE(6,*) 'Execution will terminate.'
      STOP
    ENDIF
  END SUBROUTINE get_piacc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_page_sweep(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(IN) :: wwords(lp_max)
    CHARACTER(ll_max) :: tchar1
    tchar1=TRIM(lowercase(wwords(2)))
    IF(tchar1 .EQ. 'yes') THEN
      page_sweep=1
    ELSEIF(tchar1 .EQ. 'no') THEN
      page_sweep=0
    ELSE
      WRITE(6,*) 'Error. This is not a valid sweep page option -- ',wwords(2),' --'
      WRITE(6,*) 'Execution will terminate.'
      STOP
    ENDIF
  END SUBROUTINE get_page_sweep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_page_refl(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(IN) :: wwords(lp_max)
    STOP 'thisget not yet complete'
  END SUBROUTINE get_page_refl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_page_iflw(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(IN) :: wwords(lp_max)
    STOP 'thisget not yet complete'
  END SUBROUTINE get_page_iflw

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_kconv(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(IN) :: wwords(lp_max)
    STOP 'thisget not yet complete'
  END SUBROUTINE get_kconv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_innerconv(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(IN) :: wwords(lp_max)
    STOP 'thisget not yet complete'
  END SUBROUTINE get_innerconv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_outerconv(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(IN) :: wwords(lp_max)
    STOP 'thisget not yet complete'
  END SUBROUTINE get_outerconv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_maxinner(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(IN) :: wwords(lp_max)
    STOP 'thisget not yet complete'
  END SUBROUTINE get_maxinner

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_maxouter(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(IN) :: wwords(lp_max)
    STOP 'thisget not yet complete'
  END SUBROUTINE get_maxouter

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_jfnk_krsze(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(IN) :: wwords(lp_max)
    STOP 'thisget not yet complete'
  END SUBROUTINE get_jfnk_krsze

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_jfnk_maxkr(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(IN) :: wwords(lp_max)
    STOP 'thisget not yet complete'
  END SUBROUTINE get_jfnk_maxkr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_jfnk_method(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(IN) :: wwords(lp_max)
    STOP 'thisget not yet complete'
  END SUBROUTINE get_jfnk_method

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_initial_guess(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(IN) :: wwords(lp_max)
    STOP 'thisget not yet complete'
  END SUBROUTINE get_initial_guess

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_save_restart(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(IN) :: wwords(lp_max)
    STOP 'thisget not yet complete'
  END SUBROUTINE get_save_restart

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_ipiter(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(IN) :: wwords(lp_max)
    STOP 'thisget not yet complete'
  END SUBROUTINE get_ipiter

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_print_conv(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(IN) :: wwords(lp_max)
    STOP 'thisget not yet complete'
  END SUBROUTINE get_print_conv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_density_factor(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(IN) :: wwords(lp_max)
    STOP 'thisget not yet complete'
  END SUBROUTINE get_density_factor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_execution(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(IN) :: wwords(lp_max)
    STOP 'thisget not yet complete'
  END SUBROUTINE get_execution

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_mesh_file(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(IN) :: wwords(lp_max)
    STOP 'thisget not yet complete'
  END SUBROUTINE get_mesh_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_inflow_file(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(IN) :: wwords(lp_max)
    STOP 'thisget not yet complete'
  END SUBROUTINE get_inflow_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_source_file(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(IN) :: wwords(lp_max)
    STOP 'thisget not yet complete'
  END SUBROUTINE get_source_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_flux_file(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(IN) :: wwords(lp_max)
    STOP 'thisget not yet complete'
  END SUBROUTINE get_flux_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_xs_file(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(IN) :: wwords(lp_max)
    STOP 'thisget not yet complete'
  END SUBROUTINE get_xs_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_density_factor_file(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(IN) :: wwords(lp_max)
    STOP 'thisget not yet complete'
  END SUBROUTINE get_density_factor_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_quad_file(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(IN) :: wwords(lp_max)
    STOP 'thisget not yet complete'
  END SUBROUTINE get_quad_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_vtk(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(IN) :: wwords(lp_max)
    STOP 'thisget not yet complete'
  END SUBROUTINE get_vtk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_vtk_flux_file(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(IN) :: wwords(lp_max)
    STOP 'thisget not yet complete'
  END SUBROUTINE get_vtk_flux_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_vtk_mat_file(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(IN) :: wwords(lp_max)
    STOP 'thisget not yet complete'
  END SUBROUTINE get_vtk_mat_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_vtk_reg_file(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(IN) :: wwords(lp_max)
    STOP 'thisget not yet complete'
  END SUBROUTINE get_vtk_reg_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_vtk_src_file(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(IN) :: wwords(lp_max)
    STOP 'thisget not yet complete'
  END SUBROUTINE get_vtk_src_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_restart_file(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(IN) :: wwords(lp_max)
    STOP 'thisget not yet complete'
  END SUBROUTINE get_restart_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_inguess_file(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(IN) :: wwords(lp_max)
    STOP 'thisget not yet complete'
  END SUBROUTINE get_inguess_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_cartesian_map_file(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(IN) :: wwords(lp_max)
    STOP 'thisget not yet complete'
  END SUBROUTINE get_cartesian_map_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_print_xs(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(IN) :: wwords(lp_max)
    STOP 'thisget not yet complete'
  END SUBROUTINE get_print_xs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_ngroups(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(IN) :: wwords(lp_max)
    STOP 'thisget not yet complete'
  END SUBROUTINE get_ngroups

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_pnorder(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(IN) :: wwords(lp_max)
    STOP 'thisget not yet complete'
  END SUBROUTINE get_pnorder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_pnread(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(IN) :: wwords(lp_max)
    STOP 'thisget not yet complete'
  END SUBROUTINE get_pnread

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_upscattering(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(IN) :: wwords(lp_max)
    STOP 'thisget not yet complete'
  END SUBROUTINE get_upscattering

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_multiplying(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(IN) :: wwords(lp_max)
    STOP 'thisget not yet complete'
  END SUBROUTINE get_multiplying

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_scatt_mult_included(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(IN) :: wwords(lp_max)
    STOP 'thisget not yet complete'
  END SUBROUTINE get_scatt_mult_included

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_qdtype(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(IN) :: wwords(lp_max)
    STOP 'thisget not yet complete'
  END SUBROUTINE get_qdtype

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_qdorder(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(IN) :: wwords(lp_max)
    STOP 'thisget not yet complete'
  END SUBROUTINE get_qdorder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_cartesian_map(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(IN) :: wwords(lp_max)
    STOP 'thisget not yet complete'
  END SUBROUTINE get_cartesian_map

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_point_value_locations(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(IN) :: wwords(lp_max)
    STOP 'thisget not yet complete'
  END SUBROUTINE get_point_value_locations

END MODULE read_module_current