!***********************************************************************
!
! The module contains all current subroutines for input file reading
!
!***********************************************************************
MODULE read_inp_module

  USE mpi
  USE stringmod
  USE globals
  USE error_module

  IMPLICIT NONE
  PRIVATE
  !
  ! List of public members
  PUBLIC :: inputfile_read

!> The maximum length of a cardname
INTEGER,PARAMETER :: MAX_CARDNAME_LEN=32
!> The number of cards we have
INTEGER,PARAMETER :: num_cards=43
!> The number of cards for deprecated info
INTEGER,PARAMETER :: num_dep_cards=5
!> The maximum length of a line in the input file
INTEGER,PARAMETER :: ll_max=200
!(need to change in interface IFchanged)
!> The maximum number of params on a given line or inputs for a param
INTEGER,PARAMETER :: lp_max=100
!(also need to change in interface IFchanged)
!> The local unit for the input file
INTEGER :: local_unit

TYPE :: cardType
  !character name of the card
  CHARACTER(MAX_CARDNAME_LEN) :: cname
  !logical to tell IFcard has already appeared
  LOGICAL :: used=.FALSE.
  !readin procedure, unique for each card
  PROCEDURE(prototype_wordarg),POINTER :: getcard => NULL()
  !card argument
  CHARACTER(ll_max*lp_max) :: carg
  !card subject matter
  CHARACTER(MAX_CARDNAME_LEN) :: csub
ENDTYPE cardType

TYPE(cardType) :: cards(num_cards)

!> Simple abstract interface for an object method subroutine with word array argument
ABSTRACT INTERFACE
  SUBROUTINE prototype_wordarg(thisCard,twords)
    IMPORT :: cardType
    CLASS(cardType),INTENT(INOUT) :: thisCard
    CHARACTER(200),INTENT(INOUT) :: twords(100)
  ENDSUBROUTINE prototype_wordarg
ENDINTERFACE

CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE inputfile_read(localunit)
    INTEGER, INTENT(IN) :: localunit
    CHARACTER(ll_max) :: tchar
    CHARACTER(ll_max) :: words(lp_max),wwords(lp_max)
    INTEGER :: ios,rank,mpi_err,nwords,nwwords,iw,ic,card_indx
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)
    IF(localunit .NE. 100+rank)CALL raise_fatal_error('rank mismatch in input file reading')
    local_unit=localunit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !setup the card objects, add one to here IFyou need a new one, make sure to change the total
    !number of cards param as well
    cards(:)%carg=''
    cards(:)%csub=''
    !type card
    card_indx=1
    cards(card_indx)%cname='type'
    cards(card_indx)%carg='keig'
    cards(card_indx)%getcard => get_type
    !keigsolver card
    card_indx=2
    cards(card_indx)%cname='keigsolver'
    cards(card_indx)%carg='pi'
    cards(card_indx)%getcard => get_keigsolver
    cards(card_indx)%csub='keig'
    !lambda card
    card_indx=3
    cards(card_indx)%cname='lambda'
    WRITE(cards(card_indx)%carg,'(I0)')space_ord
    cards(card_indx)%getcard => get_lambda
    !inflow card
    card_indx=4
    cards(card_indx)%cname='inflow'
    cards(card_indx)%carg='no'
    cards(card_indx)%getcard => get_inflow
    cards(card_indx)%csub='srcprob'
    !piacc card
    card_indx=5
    cards(card_indx)%cname='piacc'
    cards(card_indx)%carg='no'
    cards(card_indx)%getcard => get_piacc
    cards(card_indx)%csub='poweriter'
    !page_sweep card
    card_indx=6
    cards(card_indx)%cname='page_sweep'
    cards(card_indx)%carg='no'
    cards(card_indx)%getcard => get_page_sweep
    !page_refl card
    card_indx=7
    cards(card_indx)%cname='page_refl'
    cards(card_indx)%carg='save'
    cards(card_indx)%getcard => get_page_refl
    !page_iflw card
    card_indx=8
    cards(card_indx)%cname='page_iflw'
    cards(card_indx)%carg='all'
    cards(card_indx)%getcard => get_page_iflw
    cards(card_indx)%csub='srcprob'
    !kconv card
    card_indx=9
    cards(card_indx)%cname='kconv'
    WRITE(cards(card_indx)%carg,'(ES12.4)')k_conv
    cards(card_indx)%getcard => get_kconv
    cards(card_indx)%csub='keig'
    !innerconv card
    card_indx=10
    cards(card_indx)%cname='innerconv'
    WRITE(cards(card_indx)%carg,'(ES12.4)')inner_conv
    cards(card_indx)%getcard => get_innerconv
    !outerconv card
    card_indx=11
    cards(card_indx)%cname='outerconv'
    WRITE(cards(card_indx)%carg,'(ES12.4)')outer_conv
    cards(card_indx)%getcard => get_outerconv
    !maxinner card
    card_indx=12
    cards(card_indx)%cname='maxinner'
    WRITE(cards(card_indx)%carg,'(I0)')max_inner
    cards(card_indx)%getcard => get_maxinner
    !maxouter card
    card_indx=13
    cards(card_indx)%cname='maxouter'
    WRITE(cards(card_indx)%carg,'(I0)')max_outer
    cards(card_indx)%getcard => get_maxouter
    !jfnk_krsze card
    card_indx=14
    cards(card_indx)%cname='jfnk_krsze'
    WRITE(cards(card_indx)%carg,'(I0)')rd_restart
    cards(card_indx)%getcard => get_jfnk_krsze
    cards(card_indx)%csub='jfnk'
    !jfnk_maxkr card
    card_indx=15
    cards(card_indx)%cname='jfnk_maxkr'
    WRITE(cards(card_indx)%carg,'(I0)')rd_max_kit
    cards(card_indx)%getcard => get_jfnk_maxkr
    cards(card_indx)%csub='jfnk'
    !jfnk_method card
    card_indx=16
    cards(card_indx)%cname='jfnk_method'
    cards(card_indx)%carg='flat'
    cards(card_indx)%getcard => get_jfnk_method
    cards(card_indx)%csub='jfnk'
    !initial_guess card
    card_indx=17
    cards(card_indx)%cname='initial_guess'
    cards(card_indx)%carg='no'
    cards(card_indx)%getcard => get_initial_guess
    !restart_out card
    card_indx=18
    cards(card_indx)%cname='restart_out'
    cards(card_indx)%carg='no'
    cards(card_indx)%getcard => get_restart_out
    !ipiter card
    card_indx=19
    cards(card_indx)%cname='ipiter'
    WRITE(cards(card_indx)%carg,'(I0)')ipow
    cards(card_indx)%getcard => get_ipiter
    cards(card_indx)%csub='jfnk'
    !print_conv card
    card_indx=20
    cards(card_indx)%cname='print_conv'
    cards(card_indx)%carg='no'
    cards(card_indx)%getcard => get_print_conv
    !density_factor card
    card_indx=21
    cards(card_indx)%cname='density_factor'
    cards(card_indx)%carg='no'
    cards(card_indx)%getcard => get_density_factor
    !execution card
    card_indx=22
    cards(card_indx)%cname='execution'
    cards(card_indx)%carg='yes'
    cards(card_indx)%getcard => get_execution
    !mesh card
    card_indx=23
    cards(card_indx)%cname='mesh'
    cards(card_indx)%carg=mesh_filename
    cards(card_indx)%getcard => get_mesh
    !source card
    card_indx=24
    cards(card_indx)%cname='source'
    cards(card_indx)%carg=source_filename
    cards(card_indx)%getcard => get_source
    cards(card_indx)%csub='srcprob'
    !flux_out card
    card_indx=25
    cards(card_indx)%cname='flux_out'
    cards(card_indx)%carg=flux_filename
    cards(card_indx)%getcard => get_flux_out
    !xs card
    card_indx=26
    cards(card_indx)%cname='xs'
    cards(card_indx)%carg=cross_section_filename
    cards(card_indx)%getcard => get_xs
    !vtk_flux_out card
    card_indx=27
    cards(card_indx)%cname='vtk_flux_out'
    cards(card_indx)%carg='no'
    cards(card_indx)%getcard => get_vtk_flux_out
    !vtk_mat_out card
    card_indx=28
    cards(card_indx)%cname='vtk_mat_out'
    cards(card_indx)%carg='no'
    cards(card_indx)%getcard => get_vtk_mat_out
    !vtk_reg_out card
    card_indx=29
    cards(card_indx)%cname='vtk_reg_out'
    cards(card_indx)%carg='no'
    cards(card_indx)%getcard => get_vtk_reg_out
    !vtk_src_out card
    card_indx=30
    cards(card_indx)%cname='vtk_src_out'
    cards(card_indx)%carg='no'
    cards(card_indx)%getcard => get_vtk_src_out
    !cartesian_map_out card
    card_indx=31
    cards(card_indx)%cname='cartesian_map_out'
    cards(card_indx)%carg='no'
    cards(card_indx)%getcard => get_cartesian_map_out
    !print_xs card
    card_indx=32
    cards(card_indx)%cname='print_xs'
    cards(card_indx)%carg='no'
    cards(card_indx)%getcard => get_print_xs
    !pnorder card
    card_indx=33
    cards(card_indx)%cname='pnorder'
    WRITE(cards(card_indx)%carg,'(I0)')scatt_ord
    cards(card_indx)%getcard => get_pnorder
    !qdtype card
    card_indx=34
    cards(card_indx)%cname='qdtype'
    cards(card_indx)%carg='levelsym'
    cards(card_indx)%getcard => get_qdtype
    !qdorder card
    card_indx=35
    cards(card_indx)%cname='qdorder'
    WRITE(cards(card_indx)%carg,'(I0)')quad_ord
    cards(card_indx)%getcard => get_qdorder
    !cartesian_map card
    card_indx=36
    cards(card_indx)%cname='cartesian_map'
    cards(card_indx)%carg='no'
    cards(card_indx)%getcard => get_cartesian_map
    !point_value_locations card
    card_indx=37
    cards(card_indx)%cname='point_value_locations'
    cards(card_indx)%carg='no'
    cards(card_indx)%getcard => get_point_value_locations
    !region_map card
    card_indx=38
    cards(card_indx)%cname='region_map'
    cards(card_indx)%carg='no'
    cards(card_indx)%getcard => get_region_map
    !ngroups card
    card_indx=39
    cards(card_indx)%cname='ngroups'
    WRITE(cards(card_indx)%carg,'(I0)')egmax
    cards(card_indx)%getcard => get_ngroups
    cards(card_indx)%csub='legacyxs'
    !pnread card
    card_indx=40
    cards(card_indx)%cname='pnread'
    WRITE(cards(card_indx)%carg,'(I0)')xs_ord
    cards(card_indx)%getcard => get_pnread
    cards(card_indx)%csub='legacyxs'
    !upscattering card
    card_indx=41
    cards(card_indx)%cname='upscattering'
    cards(card_indx)%carg='yes'
    cards(card_indx)%getcard => get_upscattering
    cards(card_indx)%csub='legacyxs'
    !multiplying card
    card_indx=42
    cards(card_indx)%cname='multiplying'
    cards(card_indx)%carg='yes'
    cards(card_indx)%getcard => get_multiplying
    cards(card_indx)%csub='legacyxs'
    !scatt_mult_included card
    card_indx=43
    cards(card_indx)%cname='scatt_mult_included'
    cards(card_indx)%carg='yes'
    cards(card_indx)%getcard => get_scatt_mult_included
    cards(card_indx)%csub='legacyxs'
    !end of input cards
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    minreg= 100000_li
    maxreg=-1_li

    DO
      CALL get_next_line(tchar,ios)
      IF(ios .NE. 0)EXIT
      CALL parse(tchar,";",words,nwords)
      !loop over inputs on a line
      DO iw=1,nwords
        words(iw)=TRIM(ADJUSTL(words(iw)))
        CALL parse(words(iw)," ",wwords,nwwords)
        wwords(1)=TRIM(ADJUSTL(lowercase(wwords(1))))
        !find the card it belongs to
        DO ic=1,num_cards
          IF(wwords(1) .EQ. cards(ic)%cname)THEN
            IF(cards(ic)%used)THEN
              CALL raise_fatal_error('duplicate params: '//TRIM(cards(ic)%cname))
            ENDIF
            CALL cards(ic)%getcard(wwords)
            cards(ic)%used=.TRUE.
            EXIT
          ENDIF
        ENDDO
        IF(ic .GE. num_cards+1)THEN
          CALL raise_fatal_error('bad input, unrecognized card: '//TRIM(wwords(1)))
        ENDIF
      ENDDO
    ENDDO
    IF(rank .EQ. 0)CALL echo_equiv_inp()
  END SUBROUTINE inputfile_read

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_type(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)

    this_card%carg=TRIM(lowercase(wwords(2)))
    IF(this_card%carg .EQ. 'keig') THEN
      problem=1
    ELSEIF(this_card%carg .EQ. 'fsrc') THEN
      problem=0
    ELSE
      CALL raise_fatal_error('This is not a valid problem type -- '//TRIM(this_card%carg)//' --')
    ENDIF
  END SUBROUTINE get_type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_keigsolver(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)

    this_card%carg=TRIM(lowercase(wwords(2)))
    IF(this_card%carg .EQ. 'pi') THEN
      eig_switch=0
    ELSEIF(this_card%carg .EQ. 'jfnk') THEN
      eig_switch=1
    ELSE
      CALL raise_fatal_error('This is not a valid eigenvalue solver -- '//TRIM(this_card%carg)//' --')
    ENDIF
  END SUBROUTINE get_keigsolver

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_lambda(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    INTEGER :: ios

    this_card%carg=TRIM(lowercase(wwords(2)))
    READ(this_card%carg,'(i10)',iostat=ios) space_ord
    IF(ios.NE.0) THEN
      CALL raise_fatal_error('Invalid spatial order in problem specification -- '//TRIM(this_card%carg)//' --')
    ENDIF
  END SUBROUTINE get_lambda

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_inflow(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)

    this_card%carg=TRIM(wwords(2))
    finflow=1
    IF(lowercase(this_card%carg) .EQ. 'yes') THEN
      !do nothing, default filename
      this_card%carg=finflow_filename
    ELSEIF(lowercase(this_card%carg) .EQ. 'no' .OR. lowercase(this_card%carg) .EQ. 'none') THEN
      finflow=0
    ELSE
      finflow_filename=TRIM(this_card%carg)
    ENDIF
  END SUBROUTINE get_inflow

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_piacc(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)

    this_card%carg=TRIM(lowercase(wwords(2)))
    IF(this_card%carg .EQ. 'errmode') THEN
      outer_acc=2
    ELSEIF(this_card%carg .EQ. 'no' .OR. wwords(2) .EQ. 'none') THEN
      outer_acc=1
    ELSE
      CALL raise_fatal_error('This is not a valid acceleration option for PI -- '//TRIM(this_card%carg)//' --')
    ENDIF
  END SUBROUTINE get_piacc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_page_sweep(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)

    this_card%carg=TRIM(lowercase(wwords(2)))
    IF(this_card%carg .EQ. 'yes') THEN
      page_sweep=1
    ELSEIF(this_card%carg .EQ. 'no' .OR. this_card%carg .EQ. 'none') THEN
      page_sweep=0
    ELSE
      CALL raise_fatal_error('This is not a valid sweep page option -- '//TRIM(this_card%carg)//' --')
    ENDIF
  END SUBROUTINE get_page_sweep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_page_refl(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)

    this_card%carg=TRIM(lowercase(wwords(2)))
    IF(this_card%carg .EQ. 'page') THEN
      page_refl=1
    ELSEIF(this_card%carg .EQ. 'save') THEN
      page_refl=0
    ELSEIF(this_card%carg .EQ. 'inner') THEN
      page_refl=2
    ELSE
      CALL raise_fatal_error('This is not a valid page refl -- '//TRIM(this_card%carg)//' --')
    ENDIF
  END SUBROUTINE get_page_refl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_page_iflw(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)

    this_card%carg=TRIM(lowercase(wwords(2)))
    IF(this_card%carg .EQ. 'bygroup') THEN
      page_iflw=1
    ELSEIF(this_card%carg .EQ. 'all') THEN
      page_iflw=0
    ELSE
      CALL raise_fatal_error('This is not a valid page inflow option -- '//TRIM(this_card%carg)//' --')
    END IF
  END SUBROUTINE get_page_iflw

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_kconv(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    INTEGER :: ios

    this_card%carg=TRIM(lowercase(wwords(2)))
    READ(this_card%carg,*,iostat=ios) k_conv
    IF(ios.NE.0) THEN
      CALL raise_fatal_error('Invalid stopping criterion for keff in problem specification -- '//TRIM(this_card%carg)//' --')
    END IF
  END SUBROUTINE get_kconv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_innerconv(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    INTEGER :: ios

    this_card%carg=TRIM(lowercase(wwords(2)))
    READ(this_card%carg,*,iostat=ios) inner_conv
    IF(ios.NE.0) THEN
      CALL raise_fatal_error('Invalid stopping criterion for inner iterations in problem &
        & specification -- '//TRIM(this_card%carg)//' --')
    END IF
  END SUBROUTINE get_innerconv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_outerconv(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    INTEGER :: ios

    this_card%carg=TRIM(lowercase(wwords(2)))
    READ(this_card%carg,*,iostat=ios) outer_conv
    IF(ios.NE.0) THEN
      CALL raise_fatal_error('Invalid stopping criterion for outer iterations in problem &
        & specification -- '//TRIM(this_card%carg)//' --')
    END IF
  END SUBROUTINE get_outerconv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_maxinner(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    INTEGER :: ios

    this_card%carg=TRIM(lowercase(wwords(2)))
    READ(this_card%carg,'(i10)',iostat=ios) max_inner
    IF(ios.NE.0 .OR. max_inner<1 ) THEN
      CALL raise_fatal_error('Invalid maximum number of inner iteration in problem specification -- ' &
        //TRIM(this_card%carg)//' --')
    END IF
  END SUBROUTINE get_maxinner

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_maxouter(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    INTEGER :: ios

    this_card%carg=TRIM(lowercase(wwords(2)))
    READ(this_card%carg,'(i10)',iostat=ios) max_outer
    IF(ios.NE.0 .OR. max_outer<1 ) THEN
      CALL raise_fatal_error('Invalid maximum number of outer iteration in problem specification -- '//TRIM(this_card%carg)//' --')
    END IF
  END SUBROUTINE get_maxouter

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_jfnk_krsze(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    INTEGER :: ios

    this_card%carg=TRIM(lowercase(wwords(2)))
    READ(this_card%carg,'(i10)',iostat=ios) rd_restart
    IF(ios.NE.0 .OR. rd_restart<1 ) THEN
      CALL raise_fatal_error('Invalid Krylov subspace size for JFNK in problem specification -- '//TRIM(this_card%carg)//' --')
    END IF
  END SUBROUTINE get_jfnk_krsze

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_jfnk_maxkr(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    INTEGER :: ios

    this_card%carg=TRIM(lowercase(wwords(2)))
    READ(this_card%carg,'(i10)',iostat=ios) rd_max_kit
    IF(ios.NE.0 .OR. rd_max_kit < 1 ) THEN
      CALL raise_fatal_error('Invalid maximum number of Krylov iterations for JFNK in problem specification -- '&
            //TRIM(this_card%carg)//' --')
    ENDIF
  END SUBROUTINE get_jfnk_maxkr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_jfnk_method(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)

    this_card%carg=TRIM(lowercase(wwords(2)))
    IF(this_card%carg .EQ. 'outer') THEN
      rd_method=1
    ELSEIF(this_card%carg .EQ. 'flat') THEN
      rd_method=2
    ELSEIF(this_card%carg .EQ. 'flat_wds') THEN
      rd_method=3
    ELSE
      CALL raise_fatal_error('This is not a valid jfnk solution method -- '//TRIM(this_card%carg)//' --')
    END IF
  END SUBROUTINE get_jfnk_method

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_initial_guess(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)

    this_card%carg=TRIM(wwords(2))
    inguess_flag=1
    IF(lowercase(this_card%carg) .EQ. 'yes') THEN
      !do nothing, default
      this_card%carg=inguess_file
    ELSEIF(lowercase(this_card%carg) .EQ. 'no' .OR. lowercase(this_card%carg) .EQ. 'none') THEN
      inguess_flag=0
    ELSE
      inguess_file=TRIM(this_card%carg)
    END IF
  END SUBROUTINE get_initial_guess

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_restart_out(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)

    this_card%carg=TRIM(wwords(2))
    dump_flag=1
    IF(lowercase(this_card%carg) .EQ. 'yes') THEN
      !do nothing, default
      this_card%carg=dump_file
    ELSEIF(lowercase(this_card%carg) .EQ. 'no' .OR. lowercase(this_card%carg) .EQ. 'none') THEN
      dump_flag=0
    ELSE
      dump_file=TRIM(this_card%carg)
    END IF
  END SUBROUTINE get_restart_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_ipiter(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    INTEGER :: ios

    this_card%carg=TRIM(lowercase(wwords(2)))
    READ(this_card%carg,'(i10)',iostat=ios) ipow
    IF(ios.NE.0 ) THEN
      CALL raise_fatal_error('Invalid number of initial power iterations -- '&
            //TRIM(this_card%carg)//' --')
    END IF
  END SUBROUTINE get_ipiter

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_print_conv(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)

    this_card%carg=TRIM(lowercase(wwords(2)))
    IF(this_card%carg .EQ. 'yes') THEN
      print_conv=1
    ELSEIF(this_card%carg .EQ. 'no' .OR. wwords(2) .EQ. 'none') THEN
      print_conv=0
    ELSE
      CALL raise_fatal_error('This is not a valid execution option (yes/no) -- '//TRIM(this_card%carg)//' --')
    END IF
  END SUBROUTINE get_print_conv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_density_factor(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)

    this_card%carg=TRIM(wwords(2))
    IF(lowercase(this_card%carg) .EQ. 'no' .OR. lowercase(this_card%carg) .EQ. 'none') THEN
      dfact_opt = 0
    ELSE
      dfact_opt = 99
      dens_fact_filename=TRIM(wwords(2))
    END IF
  END SUBROUTINE get_density_factor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_execution(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)

    this_card%carg=TRIM(lowercase(wwords(2)))
    IF(this_card%carg .EQ. 'yes') THEN
      execution=1
    ELSEIF(this_card%carg .EQ. 'no' .OR. this_card%carg .EQ. 'none') THEN
      execution=0
    ELSE
      CALL raise_fatal_error('This is not a valid execution option (yes/no) -- '//TRIM(this_card%carg)//' --')
    END IF
  END SUBROUTINE get_execution

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_mesh(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)

    this_card%carg=TRIM(wwords(2))
    mesh_filename=TRIM(this_card%carg)
  END SUBROUTINE get_mesh

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_source(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)

    this_card%carg=TRIM(wwords(2))
    source_filename=TRIM(this_card%carg)
  END SUBROUTINE get_source

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_flux_out(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)

    this_card%carg=TRIM(wwords(2))
    flux_filename=TRIM(wwords(2))
  END SUBROUTINE get_flux_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_xs(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)

    this_card%carg=TRIM(wwords(2))
    cross_section_filename=TRIM(this_card%carg)
  END SUBROUTINE get_xs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_vtk_flux_out(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)

    vtk_flux_output=1
    this_card%carg=TRIM(wwords(2))
    IF(lowercase(this_card%carg) .EQ. 'no' .OR. lowercase(this_card%carg) .EQ. 'none')THEN
      vtk_flux_output=0
    ELSEIF(lowercase(this_card%carg) .EQ. 'yes')THEN
      !do nothing, default
      this_card%carg=vtk_flux_filename
    ELSE
      vtk_flux_filename=TRIM(this_card%carg)
    ENDIF
  END SUBROUTINE get_vtk_flux_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_vtk_mat_out(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)

    vtk_mat_output=1
    this_card%carg=TRIM(wwords(2))
    IF(lowercase(this_card%carg) .EQ. 'no' .OR. lowercase(this_card%carg) .EQ. 'none')THEN
      vtk_mat_output=0
    ELSEIF(lowercase(this_card%carg) .EQ. 'yes')THEN
      !do nothing, default
      this_card%carg=vtk_mat_filename
    ELSE
      vtk_mat_filename=TRIM(this_card%carg)
    ENDIF
  END SUBROUTINE get_vtk_mat_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_vtk_reg_out(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)

    vtk_reg_output=1
    this_card%carg=TRIM(wwords(2))
    IF(lowercase(this_card%carg) .EQ. 'no' .OR. lowercase(this_card%carg) .EQ. 'none')THEN
      vtk_reg_output=0
    ELSEIF(lowercase(this_card%carg) .EQ. 'yes')THEN
      !do nothing, default
      this_card%carg=vtk_reg_filename
    ELSE
      vtk_reg_filename=TRIM(this_card%carg)
    ENDIF
  END SUBROUTINE get_vtk_reg_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_vtk_src_out(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)

    vtk_src_output=1
    this_card%carg=TRIM(wwords(2))
    IF(lowercase(this_card%carg) .EQ. 'no' .OR. lowercase(this_card%carg) .EQ. 'none')THEN
      vtk_src_output=0
    ELSEIF(lowercase(this_card%carg) .EQ. 'yes')THEN
      !do nothing, default
      this_card%carg=vtk_src_filename
    ELSE
      vtk_src_filename=TRIM(this_card%carg)
    ENDIF
  END SUBROUTINE get_vtk_src_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_cartesian_map_out(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)

    this_card%carg=TRIM(wwords(2))
    cartesian_map_filename=TRIM(this_card%carg)
  END SUBROUTINE get_cartesian_map_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_print_xs(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)

    this_card%carg=TRIM(lowercase(wwords(2)))
    IF(this_card%carg .EQ. 'yes') THEN
      print_xs_flag=1
    ELSEIF(this_card%carg .EQ. 'no' .OR. this_card%carg .EQ. 'none') THEN
      print_xs_flag=0
    ELSE
      CALL raise_fatal_error('This is not a valid cross section print option (yes/no) -- '//TRIM(this_card%carg)//' --')
    END IF
  END SUBROUTINE get_print_xs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_ngroups(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    INTEGER :: ios

    this_card%carg=TRIM(lowercase(wwords(2)))
    READ(this_card%carg,'(i10)',iostat=ios) egmax
    IF(ios.NE.0 .OR. egmax<1) THEN
      CALL raise_fatal_error('Invalid number of energy groups in cross section specification -- '//TRIM(this_card%carg)//' --')
    END IF
  END SUBROUTINE get_ngroups

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_pnorder(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    INTEGER :: ios

    this_card%carg=TRIM(lowercase(wwords(2)))
    READ(wwords(2),'(i10)',iostat=ios) scatt_ord
    IF(ios.NE.0 .OR. scatt_ord < 0) THEN
      CALL raise_fatal_error('Invalid scattering expansion in cross section specification -- '//TRIM(this_card%carg)//' --')
    END IF
  END SUBROUTINE get_pnorder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_pnread(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    INTEGER :: ios

    this_card%carg=TRIM(lowercase(wwords(2)))
    READ(this_card%carg,'(i10)',iostat=ios) xs_ord
    IF(ios.NE.0 .OR. xs_ord < 0) THEN
      CALL raise_fatal_error('Invalid cross section expansion in cross section specification -- '//TRIM(this_card%carg)//' --')
    END IF
  END SUBROUTINE get_pnread

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_upscattering(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)

    this_card%carg=TRIM(lowercase(wwords(2)))
    IF(this_card%carg .EQ. 'yes') THEN
      upscattering=1
    ELSEIF(this_card%carg .EQ. 'no' .OR. this_card%carg .EQ. 'none') THEN
      upscattering=0
    ELSE
      CALL raise_fatal_error('This is not a valid upscattering flag -- '//TRIM(this_card%carg)//' --')
    END IF
  END SUBROUTINE get_upscattering

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_multiplying(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)

    this_card%carg=TRIM(lowercase(wwords(2)))
    IF(this_card%carg .EQ. 'yes') THEN
      multiplying=1
    ELSEIF(this_card%carg .EQ. 'no' .OR. this_card%carg .EQ. 'none') THEN
      multiplying=0
    ELSE
      CALL raise_fatal_error('This is not a valid multiplying flag -- '//TRIM(this_card%carg)//' --')
    END IF
  END SUBROUTINE get_multiplying

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_scatt_mult_included(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)

    this_card%carg=TRIM(lowercase(wwords(2)))
    IF(this_card%carg .EQ. 'yes') THEN
      scat_mult_flag=1
    ELSEIF(this_card%carg .EQ. 'no' .OR. this_card%carg .EQ. 'none') THEN
      scat_mult_flag=0
    ELSE
      CALL raise_fatal_error('This is not a valid scattering multiplier flag -- '//TRIM(this_card%carg)//' --')
    END IF
  END SUBROUTINE get_scatt_mult_included

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_qdtype(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)

    this_card%carg=TRIM(wwords(2))
    IF(lowercase(this_card%carg) .EQ. 'levelsym') THEN
      quad_tpe=1
    ELSEIF(lowercase(this_card%carg) .EQ. 'legcheb') THEN
      quad_tpe=2
    ELSE
      quad_tpe=3
      quad_file=TRIM(this_card%carg)
      CALL raise_fatal_error('Quadrature file given. THOR does not currently support user provided &
        & user provied quadrature files -- '//TRIM(this_card%carg)//' --')
    END IF
  END SUBROUTINE get_qdtype

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_qdorder(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    INTEGER :: ios

    this_card%carg=TRIM(lowercase(wwords(2)))
    READ(this_card%carg,'(i10)',iostat=ios) quad_ord
    IF(ios.NE.0 ) THEN
      CALL raise_fatal_error('Invalid quadrature order -- '//TRIM(this_card%carg)//' --')
    END IF
  END SUBROUTINE get_qdorder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_cartesian_map(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    INTEGER :: nwwwords,i,minint
    CHARACTER(ll_max) :: wwwords(lp_max),msg

    IF(lowercase(wwords(2)) .NE. 'none' .AND. lowercase(wwords(2)) .NE. 'no')THEN
      minint=1
      !get the cartesian map array
      this_card%carg=''
      DO i=2,lp_max
        wwwords(i-1)=TRIM(ADJUSTL(lowercase(wwords(i))))
        IF(TRIM(ADJUSTL(wwords(i))) .EQ. '')EXIT
        nwwwords=i-1
        this_card%carg=TRIM(this_card%carg)//' '//TRIM(wwwords(i-1))
      ENDDO
      glob_do_cartesian_mesh = .TRUE.
      ! wwords must be an array with of length 9
      IF(nwwwords .NE. 9) THEN
        WRITE(amsg,'(3A,I0,A)') 'Following cartesian map nine entries are required; Found: ',&
              TRIM(wwords(2)),' has ', nwwwords, ' entries.'
        CALL printlog(amsg)
      END IF
      msg='Conversion to cartesian map xmin failed'
      glob_cmap_min_x = string_to_real(wwwords(1), msg)
      msg='Conversion to cartesian map xmax failed'
      glob_cmap_max_x = string_to_real(wwwords(2), msg)
      IF(ABS(glob_cmap_max_x - glob_cmap_min_x) < small_real) THEN
        CALL printlog("cartesian_map xmin and xmax are too close to each other")
      END IF
      msg='Conversion to cartesian map nx failed'
      glob_cmap_nx = string_to_int(wwwords(3), msg, minint)
      msg='Conversion to cartesian map ymin failed'
      glob_cmap_min_y = string_to_real(wwwords(4), msg)
      msg='Conversion to cartesian map ymax failed'
      glob_cmap_max_y = string_to_real(wwwords(5), msg)
      IF(ABS(glob_cmap_max_y - glob_cmap_min_y) < small_real) THEN
        CALL printlog("cartesian_map xmin and xmax are too close to each other")
      END IF
      msg='Conversion to cartesian map ny failed'
      glob_cmap_ny = string_to_int(wwwords(6), msg, minint)
      msg='Conversion to cartesian map zmin failed'
      glob_cmap_min_z = string_to_real(wwwords(7), msg)
      msg='Conversion to cartesian map zmax failed'
      glob_cmap_max_z = string_to_real(wwwords(8), msg)
      IF(ABS(glob_cmap_max_z - glob_cmap_min_z) < small_real) THEN
        CALL printlog("cartesian_map zmin and zmax are too close to each other")
      END IF
      msg='Conversion to cartesian map nz failed'
      glob_cmap_nz = string_to_int(wwwords(9), msg, minint)
    ENDIF
  END SUBROUTINE get_cartesian_map

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_point_value_locations(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    INTEGER :: nwwwords,j,l
    CHARACTER(ll_max) :: wwwords(lp_max),msg

    IF(lowercase(wwords(2)) .NE. 'none' .AND. lowercase(wwords(2)) .NE. 'no')THEN
      msg='Conversion to point flux location failed'
      !get the point value locations array
      this_card%carg=''
      DO j=2,lp_max
        wwwords(j-1)=TRIM(ADJUSTL(lowercase(wwords(j))))
        IF(TRIM(ADJUSTL(wwords(j))) .EQ. '')EXIT
        nwwwords=j-1
        this_card%carg=TRIM(this_card%carg)//' '//TRIM(wwwords(j-1))
      ENDDO
      ! must be divisible by 3
      IF(modulo(nwwwords, 3) .ne. 0) THEN
        WRITE(amsg,'(3A,I0,A)') 'point_value_locations number of entries must be divisible by 3; Found: ',&
              TRIM(wwords(2)),' has ', nwwwords, ' entries.'
        CALL printlog(amsg)
      ELSE
        number_point_flux_locations = nwwwords / 3
        ALLOCATE(point_flux_locations(number_point_flux_locations, 3))
        DO l = 1, number_point_flux_locations
          DO j = 1, 3
            point_flux_locations(l, j) = string_to_real(wwwords((l - 1) * 3 + j),msg)
          END DO
        END DO
      END IF
    ENDIF
  END SUBROUTINE get_point_value_locations

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_region_map(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    CHARACTER(100000) :: regmap
    CHARACTER(200) :: line
    CHARACTER(MAX_CARDNAME_LEN) :: words(200)
    INTEGER :: rank,mpi_err,local_unit,i,l,lr,nwords,ios
    INTEGER,ALLOCATABLE :: tempintarray(:)

    IF(lowercase(wwords(2)) .NE. 'none' .AND. lowercase(wwords(2)) .NE. 'no')THEN
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)
      local_unit=rank+100

      regmap=""
      !get any region map info on the same line
      this_card%carg=''
      DO i=2,lp_max
        wwords(i)=ADJUSTL(TRIM(wwords(i)))
        IF(wwords(i) .EQ. "")EXIT
        l      = LEN(TRIM(regmap))
        lr     = LEN(TRIM(wwords(i)))
        regmap(l+2:l+2+lr) = TRIM(wwords(i))
        this_card%carg=TRIM(this_card%carg)//' '//TRIM(wwords(i))
      ENDDO
      !get any region map info on subsequent lines
      D1: DO
        READ(local_unit,'(A200)',IOSTAT=ios)line
        IF(ios .NE. 0)EXIT
        line=TRIM(ADJUSTL(line))
        IF(line .NE. '' .AND. line(1:1) .NE. '!')THEN
          CALL parse(line," ",words,nwords)
          words(1)=TRIM(ADJUSTL(words(1)))
          !make sure we're not onto the next param
          DO i=1,num_cards
            IF(words(1).EQ. cards(i)%cname)THEN
              BACKSPACE(local_unit)
              EXIT D1
            ENDIF
          ENDDO
          !IFit didn't exit, then it's legit input
          DO i=1,nwords
            words(i)=TRIM(ADJUSTL(words(i)))
            l      = LEN(TRIM(regmap))
            lr     = LEN(TRIM(words(i)))
            regmap(l+2:l+2+lr) = TRIM(words(i))
            this_card%carg=TRIM(this_card%carg)//' '//TRIM(wwords(i))
          ENDDO
        ENDIF
      ENDDO D1
      regmap=TRIM(ADJUSTL(regmap))
      CALL parse(regmap," ",words,nwords)
      ALLOCATE(tempintarray(nwords))
      READ(regmap,*)(tempintarray(i),i=1,nwords)
      !odd indexed are the actual region indeces
      DO i=1,nwords,2
        IF(tempintarray(i) .GE. maxreg)maxreg=tempintarray(i)
        IF(tempintarray(i) .LE. minreg)minreg=tempintarray(i)
      ENDDO
      IF(ABS(nwords-2*(maxreg-minreg+1)) .GT. 0)CALL raise_fatal_error("region map bounds and number of entries don't match")
      ALLOCATE(reg2mat(minreg:maxreg))
      !assign region mapping
      DO i=1,nwords,2
        reg2mat(tempintarray(i))=tempintarray(i+1)
      ENDDO
    ENDIF
  END SUBROUTINE get_region_map

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER(kind=li) FUNCTION string_to_int(string, msg, min_int)
    CHARACTER(lp_max), INTENT(in) :: string
    CHARACTER(100), INTENT(in) :: msg
    INTEGER(kind=li), OPTIONAL, INTENT(inout) :: min_int

    INTEGER(kind=li) :: ios

    IF(.NOT. PRESENT(min_int)) min_int = -glob_max_int
    READ(string, '(i10)', iostat=ios) string_to_int
    IF(ios .NE. 0 .OR. string_to_int < min_int) THEN
      CALL raise_fatal_error('string_to_int failed: '//TRIM(msg)//' '//TRIM(string))
    END IF
  END FUNCTION string_to_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL(kind=d_t) FUNCTION string_to_real(string, msg)
    CHARACTER(ll_max), INTENT(in) :: string
    CHARACTER(ll_max), INTENT(in) :: msg

    INTEGER(kind=li) :: ios

    READ(string, *, iostat=ios) string_to_real
    IF(ios .NE. 0) THEN
      CALL raise_fatal_error('string_to_real failed: '//TRIM(msg)//' '//TRIM(string))
    END IF
  END FUNCTION string_to_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_next_line(line,ios)
    CHARACTER(ll_max),INTENT(OUT) :: line
    INTEGER,INTENT(OUT) :: ios
    CHARACTER(ll_max) :: words(100)
    INTEGER :: nwords
    DO
      READ(local_unit,'(A10000)',IOSTAT=ios)line
      IF(ios .NE. 0)EXIT
      line=TRIM(ADJUSTL(line))
      !finding uncommented line that isn't empty
      IF(line(1:1) .NE. '!' .AND. line .NE. '')THEN
        !ignore commented portions of line
        CALL parse(line,'!',words,nwords)
        line=TRIM(ADJUSTL(words(1)))
        EXIT
      ENDIF
    ENDDO
  ENDSUBROUTINE get_next_line

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE echo_equiv_inp
    INTEGER :: i,strlg_max,cnamelg_max

    strlg_max=0
    cnamelg_max=0
    DO i=1,num_cards
      cards(i)%carg=TRIM(ADJUSTL(cards(i)%carg))
      strlg_max=MAX(strlg_max,LEN_TRIM(cards(i)%carg))
      cards(i)%cname=TRIM(ADJUSTL(cards(i)%cname))
      cnamelg_max=MAX(cnamelg_max,LEN_TRIM(cards(i)%cname))
    ENDDO
    strlg_max=MIN(strlg_max,MAX_CARDNAME_LEN)
    CALL printlog('')
    CALL printlog('***********************************************************************')
    CALL printlog('*******************Echoing verbose equivalent input********************')
    CALL printlog('***********************************************************************')
    DO i=1,num_cards
      IF(cards(i)%used .OR. i .LE. num_cards-num_dep_cards)THEN
        cards(i)%carg=TRIM(ADJUSTL(cards(i)%carg))
        !echo cards and data
        IF(LEN_TRIM(cards(i)%carg) .LE. MAX_CARDNAME_LEN)THEN
          WRITE(amsg,'(4A)')cards(i)%cname(1:cnamelg_max),' ',cards(i)%carg(1:strlg_max),' !'
          CALL printlog(amsg,ADVANCING=.FALSE.)
        ELSE
          WRITE(amsg,'(4A)')cards(i)%cname(1:cnamelg_max),' ',TRIM(cards(i)%carg),' !'
          CALL printlog(amsg,ADVANCING=.FALSE.)
        ENDIF
        !echo if card is provided
        IF(cards(i)%used)THEN
          CALL printlog('card provided',ADVANCING=.FALSE.)
        ELSE
          CALL printlog('card NOT provided (default assumed)',ADVANCING=.FALSE.)
        ENDIF
        !print info on ignored input cards
        SELECT CASE(cards(i)%csub)
          CASE('keig')
            IF(problem .NE. 1)THEN
              CALL printlog(': ignored, not a keig problem',ADVANCING=.FALSE.)
            ENDIF
          CASE('poweriter')
            IF(eig_switch .NE. 0)THEN
              CALL printlog(': ignored, not using power iterations',ADVANCING=.FALSE.)
            ENDIF
          CASE('jfnk')
            IF(eig_switch .NE. 1)THEN
              CALL printlog(': ignored, not using jfnk',ADVANCING=.FALSE.)
            ENDIF
          CASE('legacyxs')
            CALL printlog(': only used with legaxy xs',ADVANCING=.FALSE.)
          CASE('srcprob')
            IF(problem .NE. 0)THEN
              CALL printlog(': ignored, not a fixed source problem',ADVANCING=.FALSE.)
            ENDIF
          CASE DEFAULT
        ENDSELECT
        CALL printlog('')
      ENDIF
    ENDDO
    CALL printlog('***********************************************************************')
    CALL printlog('******************Equivalent complete input finished*******************')
    CALL printlog('***********************************************************************')
  ENDSUBROUTINE echo_equiv_inp

END MODULE read_inp_module