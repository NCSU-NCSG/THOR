!***********************************************************************
!
! The module contains all legacy subroutines for old yaml_input versions. These
! subroutines should generally not be changed to garauntee backwards
! compatibility.
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
INTEGER,PARAMETER :: num_cards=42
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
    CHARACTER(200),INTENT(INOUT) :: twords(100)
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    DO
      READ(local_unit,'(A200)',IOSTAT=ios) tchar
      IF(ios .NE. 0)EXIT
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
      ENDIF
    ENDDO
    STOP 'inputfile_read not yet complete'
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
      finflow_filename='finflow.txt'
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
    ELSEIF(wwords(2) .EQ. 'none') THEN
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
      inguess_file='initial_guess.txt'
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
      dump_file='restart.out'
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
      ELSE
        dens_fact_filename='density_factor.txt'
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
      source_filename='source.txt'
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
      flux_filename='flux.out'
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
      vtk_flux_filename='vtk_flux.out'
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
      vtk_mat_filename='vtk_mat.out'
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
      vtk_reg_filename='vtk_reg.out'
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
      vtk_src_filename='vtk_src.out'
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
      cartesian_map_filename='cartesian_map.out'
    ELSE IF ( lowercase(wwords(2)) .EQ. 'no' .OR. lowercase(wwords(2)) .EQ. 'none') THEN
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
    STOP 'get_cartesian_map not yet complete'
  END SUBROUTINE get_cartesian_map

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_point_value_locations(this_card,wwords)
    CLASS(cardType),INTENT(INOUT) :: this_card
    CHARACTER(ll_max),INTENT(INOUT) :: wwords(lp_max)
    STOP 'get_point_value_locations not yet complete'
  END SUBROUTINE get_point_value_locations

END MODULE read_inp_module