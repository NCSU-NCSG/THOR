!***********************************************************************
!
! The module contains all current subroutines for input file reading
!
!***********************************************************************
MODULE read_inp_module

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
INTEGER,PARAMETER :: num_cards=43
!> The maximum length of a line in the input file
INTEGER,PARAMETER :: ll_max=200
!(also need to change in interface if changed)
!> The maximum number of params on a given line or inputs for a param
INTEGER,PARAMETER :: lp_max=100
!(also need to change in interface if changed)
!> The local unit for the input file
INTEGER :: local_unit

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
    CHARACTER(200),INTENT(INOUT) :: twords(100)
  ENDSUBROUTINE prototype_wordarg
ENDINTERFACE

CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE inputfile_read(localunit)
    INTEGER, INTENT(IN) :: localunit
    CHARACTER(ll_max) :: tchar
    CHARACTER(ll_max) :: words(lp_max),wwords(lp_max)
    INTEGER :: ios,rank,mpi_err,nwords,nwwords,iw,iww,ic
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)
    IF(localunit .NE. 100+rank)STOP 'rank mismatch in input file reading'
    local_unit=localunit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !setup the card objects, add one to here if you need a new one, make sure to change the total
    !number of cards param as well
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
    cards(18)%cname='restart_out'
    cards(18)%getcard => get_restart_out
    cards(19)%cname='ipiter'
    cards(19)%getcard => get_ipiter
    cards(20)%cname='print_conv'
    cards(20)%getcard => get_print_conv
    cards(21)%cname='density_factor'
    cards(21)%getcard => get_density_factor
    cards(22)%cname='execution'
    cards(22)%getcard => get_execution
    cards(23)%cname='mesh'
    cards(23)%getcard => get_mesh
    cards(24)%cname='source'
    cards(24)%getcard => get_source
    cards(25)%cname='flux_out'
    cards(25)%getcard => get_flux_out
    cards(26)%cname='xs'
    cards(26)%getcard => get_xs
    cards(27)%cname='vtk_flux_out'
    cards(27)%getcard => get_vtk_flux_out
    cards(28)%cname='vtk_mat_out'
    cards(28)%getcard => get_vtk_mat_out
    cards(29)%cname='vtk_reg_out'
    cards(29)%getcard => get_vtk_reg_out
    cards(30)%cname='vtk_src_out'
    cards(30)%getcard => get_vtk_src_out
    cards(31)%cname='cartesian_map_out'
    cards(31)%getcard => get_cartesian_map_out
    cards(32)%cname='print_xs'
    cards(32)%getcard => get_print_xs
    cards(33)%cname='ngroups'
    cards(33)%getcard => get_ngroups
    cards(34)%cname='pnorder'
    cards(34)%getcard => get_pnorder
    cards(35)%cname='pnread'
    cards(35)%getcard => get_pnread
    cards(36)%cname='upscattering'
    cards(36)%getcard => get_upscattering
    cards(37)%cname='multiplying'
    cards(37)%getcard => get_multiplying
    cards(38)%cname='scatt_mult_included'
    cards(38)%getcard => get_scatt_mult_included
    cards(39)%cname='qdtype'
    cards(39)%getcard => get_qdtype
    cards(40)%cname='qdorder'
    cards(40)%getcard => get_qdorder
    cards(41)%cname='cartesian_map'
    cards(41)%getcard => get_cartesian_map
    cards(42)%cname='point_value_locations'
    cards(42)%getcard => get_point_value_locations
    cards(42)%cname='region_map'
    cards(42)%getcard => get_region_map
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
              WRITE(*,*) 'duplicate params: ',cards(ic)%cname
              STOP
            ENDIF
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
    ENDDO
  END SUBROUTINE inputfile_read

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_type(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    wwords(2)=TRIM(lowercase(wwords(2)))
    IF(wwords(2) .EQ. 'keig') THEN
      problem=1
    ELSEIF(wwords(2) .EQ. 'fsrc') THEN
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
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    wwords(2)=TRIM(lowercase(wwords(2)))
    IF(wwords(2) .EQ. 'pi') THEN
      eig_switch=0
    ELSEIF(wwords(2) .EQ. 'jfnk') THEN
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
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    INTEGER :: ios
    wwords(2)=TRIM(lowercase(wwords(2)))
    READ(wwords(2),'(i10)',iostat=ios) space_ord
    IF(ios.NE.0) THEN
      WRITE(6,*) 'Invalid spatial order in problem specification -- ',TRIM(wwords(2)),' --'
      STOP
    ENDIF
  END SUBROUTINE get_lambda

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_inflow(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    wwords(2)=TRIM(wwords(2))
    finflow=1
    IF(lowercase(wwords(2)) .EQ. 'yes') THEN
      !do nothing, default filename
    ELSEIF(lowercase(wwords(2)) .EQ. 'no' .OR. lowercase(wwords(2)) .EQ. 'none') THEN
      finflow=0
    ELSE
      finflow_filename=wwords(2)
    ENDIF
  END SUBROUTINE get_inflow

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_piacc(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    wwords(2)=TRIM(lowercase(wwords(2)))
    IF(wwords(2) .EQ. 'errmode') THEN
      outer_acc=2
    ELSEIF(wwords(2) .EQ. 'no' .OR. wwords(2) .EQ. 'none') THEN
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
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    wwords(2)=TRIM(lowercase(wwords(2)))
    IF(wwords(2) .EQ. 'yes') THEN
      page_sweep=1
    ELSEIF(wwords(2) .EQ. 'no' .OR. wwords(2) .EQ. 'none') THEN
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
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    wwords(2)=TRIM(lowercase(wwords(2)))
    IF(wwords(2) .EQ. 'page') THEN
      page_refl=1
    ELSEIF(wwords(2) .EQ. 'save') THEN
      page_refl=0
    ELSEIF(wwords(2) .EQ. 'inner') THEN
      page_refl=2
    ELSE
      WRITE(6,*) 'Error. This is not a valid page reflective BC option -- ',wwords(2),' --'
      WRITE(6,*) 'Execution will terminate.'
      STOP
    ENDIF
  END SUBROUTINE get_page_refl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_page_iflw(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    wwords(2)=TRIM(lowercase(wwords(2)))
    IF(wwords(2) .EQ. 'bygroup') THEN
      page_iflw=1
    ELSEIF(wwords(2) .EQ. 'all') THEN
      page_iflw=0
    ELSE
      WRITE(6,*) 'Error. This is not a valid page inflow option -- ',wwords(2),' --'
      WRITE(6,*) 'Execution will terminate.'
      STOP
    END IF
  END SUBROUTINE get_page_iflw

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_kconv(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    INTEGER :: ios
    wwords(2)=TRIM(lowercase(wwords(2)))
    READ(wwords(2),*,iostat=ios) k_conv
    IF(ios.NE.0) THEN
      WRITE(6,*) 'Invalid stopping criterion for keff in problem specification -- ',TRIM(wwords(2)),' --'
      STOP
    END IF
  END SUBROUTINE get_kconv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_innerconv(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    INTEGER :: ios
    wwords(2)=TRIM(lowercase(wwords(2)))
    READ(wwords(2),*,iostat=ios) inner_conv
    IF(ios.NE.0) THEN
      WRITE(6,*) 'Invalid stopping criterion for inner iterations in problem specification -- ',TRIM(wwords(2)),' --'
      STOP
    END IF
  END SUBROUTINE get_innerconv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_outerconv(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    INTEGER :: ios
    wwords(2)=TRIM(lowercase(wwords(2)))
    READ(wwords(2),*,iostat=ios) outer_conv
    IF(ios.NE.0) THEN
      WRITE(6,*) 'Invalid stopping criterion for outer iterations in problem specification -- ',TRIM(wwords(2)),' --'
      STOP
    END IF
  END SUBROUTINE get_outerconv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_maxinner(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    INTEGER :: ios
    wwords(2)=TRIM(lowercase(wwords(2)))
    READ(wwords(2),'(i10)',iostat=ios) max_inner
    IF(ios.NE.0 .OR. max_inner<1 ) THEN
      WRITE(6,*) 'Invalid maximum number of inner iteration in problem specification -- ',TRIM(wwords(2)),' --'
      STOP
    END IF
  END SUBROUTINE get_maxinner

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_maxouter(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    INTEGER :: ios
    wwords(2)=TRIM(lowercase(wwords(2)))
    READ(wwords(2),'(i10)',iostat=ios) max_outer
    IF(ios.NE.0 .OR. max_outer<1 ) THEN
      WRITE(6,*) 'Invalid maximum number of outer iteration in problem specification -- ',TRIM(wwords(2)),' --'
      STOP
    END IF
  END SUBROUTINE get_maxouter

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_jfnk_krsze(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    INTEGER :: ios
    wwords(2)=TRIM(lowercase(wwords(2)))
    READ(wwords(2),'(i10)',iostat=ios) rd_restart
    IF(ios.NE.0 .OR. rd_restart<1 ) THEN
      WRITE(6,*) 'Invalid Krylov subspace size for JFNK in problem specification -- ',TRIM(wwords(2)),' --'
      STOP
    END IF
  END SUBROUTINE get_jfnk_krsze

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_jfnk_maxkr(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    INTEGER :: ios
    wwords(2)=TRIM(lowercase(wwords(2)))
    READ(wwords(2),'(i10)',iostat=ios) rd_max_kit
    IF(ios.NE.0 .OR. rd_max_kit < 1 ) THEN
      WRITE(6,*) 'Invalid maximum number of Krylov iterations for JFNK in problem specification -- '&
            ,TRIM(wwords(2)),' --'
      STOP
    ENDIF
  END SUBROUTINE get_jfnk_maxkr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_jfnk_method(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    wwords(2)=TRIM(lowercase(wwords(2)))
    IF(wwords(2) .EQ. 'outer') THEN
      rd_method=1
    ELSEIF(wwords(2) .EQ. 'flat') THEN
      rd_method=2
    ELSEIF(wwords(2) .EQ. 'flat_wds') THEN
      rd_method=3
    ELSE
      WRITE(6,*) 'Error. This is not a valid jfnk solution method -- ',wwords(2),' --'
      WRITE(6,*) 'Execution will terminate.'
      STOP
    END IF
  END SUBROUTINE get_jfnk_method

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_initial_guess(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    wwords(2)=TRIM(wwords(2))
    inguess_flag=1
    IF      ( lowercase(wwords(2)) .EQ. 'yes') THEN
      !do nothing, default
    ELSE IF ( lowercase(wwords(2)) .EQ. 'no' .OR. lowercase(wwords(2)) .EQ. 'none') THEN
      inguess_flag=0
    ELSE
      inguess_file=wwords(2)
    END IF
  END SUBROUTINE get_initial_guess

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_restart_out(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    wwords(2)=TRIM(wwords(2))
    dump_flag=1
    IF      ( lowercase(wwords(2)) .EQ. 'yes') THEN
      !do nothing, default
    ELSE IF ( lowercase(wwords(2)) .EQ. 'no' .OR. lowercase(wwords(2)) .EQ. 'none') THEN
      dump_flag=0
    ELSE
      dump_file=wwords(2)
    END IF
  END SUBROUTINE get_restart_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_ipiter(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    INTEGER :: ios
    wwords(2)=TRIM(lowercase(wwords(2)))
    READ(wwords(2),'(i10)',iostat=ios) ipow
    IF(ios.NE.0 ) THEN
      WRITE(6,*) 'Invalid number of initial power iterations -- '&
            ,TRIM(wwords(2)),' --'
      STOP
    END IF
  END SUBROUTINE get_ipiter

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_print_conv(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    wwords(2)=TRIM(lowercase(wwords(2)))
    IF      ( wwords(2) .EQ. 'yes') THEN
      print_conv=1
    ELSE IF ( wwords(2) .EQ. 'no' .OR. wwords(2) .EQ. 'none') THEN
      print_conv=0
    ELSE
      WRITE(6,*) 'Error. This is not a valid execution option (yes/no) -- ',wwords(2),' --'
      WRITE(6,*) 'Execution will terminate.'
      STOP
    END IF
  END SUBROUTINE get_print_conv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_density_factor(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    wwords(2)=TRIM(lowercase(wwords(2)))
    IF      ( wwords(2) .EQ. 'no' .OR. wwords(2) .EQ. 'none') THEN
      dfact_opt = 0
    ELSE
      IF      ( wwords(2) .EQ. 'byvolume') THEN
        dfact_opt = 1
      ELSE IF ( wwords(2) .EQ. 'fromfile') THEN
        dfact_opt = 2
      ELSE
        WRITE(6,*) 'Error. This is not a valid density factor option &
              (no/byvolume/fromfile) -- ',wwords(2),' --'
        WRITE(6,*) 'Execution will terminate.'
        STOP
      END IF
    END IF
    wwords(3)=TRIM(wwords(3))
    IF(dfact_opt .NE. 0)THEN
      IF(lowercase(wwords(3)) .NE. '')THEN
        dens_fact_filename=wwords(3)
      ENDIF
    ENDIF
  END SUBROUTINE get_density_factor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_execution(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    wwords(2)=TRIM(lowercase(wwords(2)))
    IF      ( wwords(2) .EQ. 'yes') THEN
      execution=1
    ELSE IF ( wwords(2) .EQ. 'no' .OR. wwords(2) .EQ. 'none') THEN
      execution=0
    ELSE
      WRITE(6,*) 'Error. This is not a valid execution option (yes/no) -- ',wwords(2),' --'
      WRITE(6,*) 'Execution will terminate.'
      STOP
    END IF
  END SUBROUTINE get_execution

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_mesh(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    mesh_filename=TRIM(wwords(2))
  END SUBROUTINE get_mesh

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_source(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    wwords(2)=TRIM(wwords(2))
    IF      ( lowercase(wwords(2)) .EQ. 'yes') THEN
      !do nothing, default
    ELSE IF ( lowercase(wwords(2)) .EQ. 'no' .OR. lowercase(wwords(2)) .EQ. 'none') THEN
    ELSE
      source_filename=wwords(2)
    ENDIF
  END SUBROUTINE get_source

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_flux_out(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    wwords(2)=TRIM(wwords(2))
    IF      ( lowercase(wwords(2)) .EQ. 'yes') THEN
      !do nothing, default
    ELSE IF ( lowercase(wwords(2)) .EQ. 'no' .OR. lowercase(wwords(2)) .EQ. 'none') THEN
    ELSE
      flux_filename=wwords(2)
    ENDIF
  END SUBROUTINE get_flux_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_xs(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    cross_section_filename=TRIM(wwords(2))
  END SUBROUTINE get_xs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_vtk_flux_out(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    vtk_flux_output=1
    wwords(2)=TRIM(wwords(2))
    IF(lowercase(wwords(2)) .EQ. 'no' .OR. lowercase(wwords(2)) .EQ. 'none')THEN
      vtk_flux_output=0
    ELSEIF(lowercase(wwords(2)) .EQ. 'yes')THEN
      !do nothing, default
    ELSE
      vtk_flux_filename=wwords(2)
    ENDIF
  END SUBROUTINE get_vtk_flux_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_vtk_mat_out(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    vtk_mat_output=1
    wwords(2)=TRIM(wwords(2))
    IF(lowercase(wwords(2)) .EQ. 'no' .OR. lowercase(wwords(2)) .EQ. 'none')THEN
      vtk_mat_output=0
    ELSEIF(lowercase(wwords(2)) .EQ. 'yes')THEN
      !do nothing, default
    ELSE
      vtk_mat_filename=wwords(2)
    ENDIF
  END SUBROUTINE get_vtk_mat_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_vtk_reg_out(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    vtk_reg_output=1
    wwords(2)=TRIM(wwords(2))
    IF(lowercase(wwords(2)) .EQ. 'no' .OR. lowercase(wwords(2)) .EQ. 'none')THEN
      vtk_reg_output=0
    ELSEIF(lowercase(wwords(2)) .EQ. 'yes')THEN
      !do nothing, default
    ELSE
      vtk_reg_filename=wwords(2)
    ENDIF
  END SUBROUTINE get_vtk_reg_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_vtk_src_out(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    vtk_src_output=1
    wwords(2)=TRIM(wwords(2))
    IF(lowercase(wwords(2)) .EQ. 'no' .OR. lowercase(wwords(2)) .EQ. 'none')THEN
      vtk_src_output=0
    ELSEIF(lowercase(wwords(2)) .EQ. 'yes')THEN
      !do nothing, default
    ELSE
      vtk_src_filename=wwords(2)
    ENDIF
  END SUBROUTINE get_vtk_src_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_cartesian_map_out(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    wwords(2)=TRIM(wwords(2))
    IF      ( lowercase(wwords(2)) .EQ. 'yes') THEN
      !do nothing, default
    ELSE IF ( lowercase(wwords(2)) .EQ. 'no' .OR. lowercase(wwords(2)) .EQ. 'none') THEN
      glob_do_cartesian_mesh = .FALSE.
    ELSE
      cartesian_map_filename=wwords(2)
    ENDIF
  END SUBROUTINE get_cartesian_map_out

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_print_xs(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    wwords(2)=TRIM(lowercase(wwords(2)))
    IF      ( wwords(2) .EQ. 'yes') THEN
      print_xs_flag=1
    ELSE IF ( wwords(2) .EQ. 'no' .OR. wwords(2) .EQ. 'none') THEN
      print_xs_flag=0
    ELSE
      WRITE(6,*) 'Error. This is not a valid cross section print option (yes/no) -- ',wwords(2),' --'
      WRITE(6,*) 'Execution will terminate.'
      STOP
    END IF
  END SUBROUTINE get_print_xs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_ngroups(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    INTEGER :: ios
    wwords(2)=TRIM(lowercase(wwords(2)))
    READ(wwords(2),'(i10)',iostat=ios) egmax
    IF(ios.NE.0 .OR. egmax<1) THEN
      WRITE(6,*) 'Invalid number of energy groups in cross section specification -- ',TRIM(wwords(2)),' --'
      STOP
    END IF
  END SUBROUTINE get_ngroups

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_pnorder(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    INTEGER :: ios
    wwords(2)=TRIM(lowercase(wwords(2)))
    READ(wwords(2),'(i10)',iostat=ios) scatt_ord
    IF(ios.NE.0 .OR. scatt_ord < 0) THEN
      WRITE(6,*) 'Invalid scattering expansion in cross section specification -- ',TRIM(wwords(2)),' --'
      STOP
    END IF
  END SUBROUTINE get_pnorder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_pnread(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    INTEGER :: ios
    wwords(2)=TRIM(lowercase(wwords(2)))
    READ(wwords(2),'(i10)',iostat=ios) xs_ord
    IF(ios.NE.0 .OR. xs_ord < 0) THEN
      WRITE(6,*) 'Invalid cross section expansion in cross section specification -- ',TRIM(wwords(2)),' --'
      STOP
    END IF
  END SUBROUTINE get_pnread

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_upscattering(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    wwords(2)=TRIM(lowercase(wwords(2)))
    IF      ( wwords(2) .EQ. 'yes') THEN
      upscattering=1
    ELSE IF ( wwords(2) .EQ. 'no' .OR. wwords(2) .EQ. 'none') THEN
      upscattering=0
    ELSE
      WRITE(6,*) 'Error. This is not a valid upscattering flag -- ',wwords(2),' --'
      WRITE(6,*) 'Execution will terminate.'
      STOP
    END IF
  END SUBROUTINE get_upscattering

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_multiplying(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    wwords(2)=TRIM(lowercase(wwords(2)))
    IF      ( wwords(2) .EQ. 'yes') THEN
      multiplying=1
    ELSE IF ( wwords(2) .EQ. 'no' .OR. wwords(2) .EQ. 'none') THEN
      multiplying=0
    ELSE
      WRITE(6,*) 'Error. This is not a valid multiplying flag -- ',wwords(2),' --'
      WRITE(6,*) 'Execution will terminate.'
      STOP
    END IF
  END SUBROUTINE get_multiplying

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_scatt_mult_included(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    wwords(2)=TRIM(lowercase(wwords(2)))
    IF      ( wwords(2) .EQ. 'yes') THEN
      scat_mult_flag=1
    ELSE IF ( wwords(2) .EQ. 'no' .OR. wwords(2) .EQ. 'none') THEN
      scat_mult_flag=0
    ELSE
      WRITE(6,*) 'Error. This is not a valid scattering multiplier flag -- ',wwords(2),' --'
      WRITE(6,*) 'Execution will terminate.'
      STOP
    END IF
  END SUBROUTINE get_scatt_mult_included

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_qdtype(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    wwords(2)=TRIM(wwords(2))
    IF      ( lowercase(wwords(2)) .EQ. 'levelsym') THEN
      quad_tpe=1
    ELSE IF ( lowercase(wwords(2)) .EQ. 'legcheb') THEN
      quad_tpe=2
    ELSE
      quad_tpe=3
      quad_file=wwords(2)
    END IF
  END SUBROUTINE get_qdtype

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_qdorder(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    INTEGER :: ios
    wwords(2)=TRIM(lowercase(wwords(2)))
    READ(wwords(2),'(i10)',iostat=ios) quad_ord
    IF(ios.NE.0 ) THEN
      WRITE(6,*) 'Invalid quadrature order -- ',TRIM(wwords(2)),' --'
      STOP
    END IF
  END SUBROUTINE get_qdorder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_cartesian_map(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    INTEGER :: nwwwords
    CHARACTER(ll_max) :: wwwords(lp_max)
    wwords(2)=TRIM(lowercase(wwords(2)))
    glob_do_cartesian_mesh = .TRUE.
    ! wwords must be an array with of length 9
    CALL parse(wwords(2), " ", wwwords, nwwwords)
    IF (nwwwords .NE. 9) THEN
      WRITE(6,*) 'Following cartesian map nine entries are required; Found: ',&
            TRIM(wwords(2)),' has ', nwwwords, ' entries.'
    END IF
    glob_cmap_min_x = string_to_real(wwwords(1), 'Conversion to cartesian map xmin failed')
    glob_cmap_max_x = string_to_real(wwwords(2), 'Conversion to cartesian map xmax failed')
    IF (ABS(glob_cmap_max_x - glob_cmap_min_x) < small_real) THEN
      WRITE(6, *) "cartesian_map xmin and xmax are too close to each other"
    END IF
    glob_cmap_nx = string_to_int(wwwords(3), 'Conversion to cartesian map nx failed', 1)
    glob_cmap_min_y = string_to_real(wwwords(4), 'Conversion to cartesian map ymin failed')
    glob_cmap_max_y = string_to_real(wwwords(5), 'Conversion to cartesian map ymax failed')
    IF (ABS(glob_cmap_max_y - glob_cmap_min_y) < small_real) THEN
      WRITE(6, *) "cartesian_map xmin and xmax are too close to each other"
    END IF
    glob_cmap_ny = string_to_int(wwwords(6), 'Conversion to cartesian map ny failed', 1)
    glob_cmap_min_z = string_to_real(wwwords(7), 'Conversion to cartesian map zmin failed')
    glob_cmap_max_z = string_to_real(wwwords(8), 'Conversion to cartesian map zmax failed')
    IF (ABS(glob_cmap_max_z - glob_cmap_min_z) < small_real) THEN
      WRITE(6, *) "cartesian_map zmin and zmax are too close to each other"
    END IF
    glob_cmap_nz = string_to_int(wwwords(9), 'Conversion to cartesian map nz failed', 1)
  END SUBROUTINE get_cartesian_map

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_point_value_locations(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    INTEGER :: nwwwords,j,l
    CHARACTER(ll_max) :: wwwords(lp_max)
    wwords(2)=TRIM(lowercase(wwords(2)))
    CALL parse(wwords(2), " ", wwwords, nwwwords)
    ! must be divisible by 3
    IF (modulo(nwwwords, 3) .ne. 0) THEN
      WRITE(6,*) 'point_value_locations number of entries must be divisible by 3; Found: ',&
            TRIM(wwords(2)),' has ', nwwwords, ' entries.'
    ELSE
      number_point_flux_locations = nwwwords / 3
      ALLOCATE(point_flux_locations(number_point_flux_locations, 3))
      DO l = 1, number_point_flux_locations
        DO j = 1, 3
          point_flux_locations(l, j) = string_to_real(wwwords((l - 1) * 3 + j),&
            'Conversion to point flux location failed')
        END DO
      END DO
    END IF
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
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)
    local_unit=rank+100

    regmap=""
    !get any region map info on the same line
    DO i=2,lp_max
      wwords(i)=ADJUSTL(TRIM(wwords(i)))
      IF(wwords(i) .EQ. "")EXIT
      l      = LEN(TRIM(regmap))
      lr     = LEN(TRIM(wwords(i)))
      regmap(l+2:l+2+lr) = TRIM(wwords(i))
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
        !if it didn't exit, then it's legit input
        DO i=1,nwords
          words(i)=TRIM(ADJUSTL(words(i)))
          l      = LEN(TRIM(regmap))
          lr     = LEN(TRIM(words(i)))
          regmap(l+2:l+2+lr) = TRIM(words(i))
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
    IF(ABS(nwords-2*(maxreg-minreg+1)) .GT. 0) &
      STOP "region map bounds and number of entries don't match"
    ALLOCATE(reg2mat(minreg:maxreg))
    !assign region mapping
    DO i=1,nwords,2
      reg2mat(tempintarray(i))=tempintarray(i+1)
    ENDDO
  END SUBROUTINE get_region_map

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER(kind=li) FUNCTION string_to_int(string, msg, min_int)
    CHARACTER(lp_max), INTENT(in) :: string
    CHARACTER(100), INTENT(in) :: msg
    INTEGER(kind=li), OPTIONAL, INTENT(inout) :: min_int

    INTEGER(kind=li) :: ios

    IF (.NOT. PRESENT(min_int)) min_int = -glob_max_int
    READ(string, '(i10)', iostat=ios) string_to_int
    IF(ios .NE. 0 .OR. string_to_int < min_int) THEN
      WRITE(6,*) TRIM(msg), TRIM(string)
      STOP
    END IF
  END FUNCTION string_to_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL(kind=d_t) FUNCTION string_to_real(string, msg)
    CHARACTER(lp_max), INTENT(in) :: string
    CHARACTER(lp_max), INTENT(in) :: msg

    INTEGER(kind=li) :: ios

    READ(string, *, iostat=ios) string_to_real
    IF(ios .NE. 0) THEN
      WRITE(6,*) TRIM(msg), TRIM(string)
      STOP
    END IF
  END FUNCTION string_to_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_next_line(line,ios)
    CHARACTER(200),INTENT(OUT) :: line
    INTEGER,INTENT(OUT) :: ios
    CHARACTER(200) :: words(100)
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
  ENDSUBROUTINE

END MODULE read_inp_module