!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/mc [1.0] *
!!
!!   NAME
!!     mc_control
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     MJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE mc_control
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dbl
  USE parser,                          ONLY: get_next,&
                                             p_error,&
                                             parser_end,&
                                             parser_init,&
                                             read_line,&
                                             search_label,&
                                             stop_parser,&
                                             test_next
  USE physcon,                         ONLY: boltzmann,joule
  USE string_utilities,                ONLY: uppercase

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_mc_section, mc_simulation_parameters_type


!!***
!****************************************************************************

  TYPE mc_simulation_parameters_type
     INTEGER :: nstep
     INTEGER :: iprint
     INTEGER :: iuptraion
     INTEGER :: iuptrans
     INTEGER :: iblock
     INTEGER :: nbox
     INTEGER :: nchain
     INTEGER :: nunits
     INTEGER :: nunits_tot
     INTEGER :: nmolty
     INTEGER, DIMENSION ( : ), POINTER :: moltyp
     INTEGER, DIMENSION ( : ), POINTER :: ininch
     INTEGER, DIMENSION ( : ), POINTER :: inix
     INTEGER, DIMENSION ( : ), POINTER :: iniy
     INTEGER, DIMENSION ( : ), POINTER :: iniz
     INTEGER :: iratv
     INTEGER :: iratp
     CHARACTER ( LEN = 3 ) :: lstop
     CHARACTER ( LEN = 20 ) :: ensemble
     CHARACTER ( LEN = 50 ) :: program
     REAL ( dbl ) :: rmion
     REAL ( dbl ) :: rmrot
     REAL ( dbl ) :: rmtrans
     REAL ( dbl ) :: temp
     REAL ( dbl ) :: pmtraion
     REAL ( dbl ) :: pmtrans
     REAL ( dbl ) :: BETA
  END TYPE mc_simulation_parameters_type

!!***
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** mc/read_mc_section [1.0] *
!!
!!   NAME
!!     read_mc_section
!!
!!   SYNOPSIS
!!     Subroutine read_mc_section(mc_par, mcpar, mcio)
!!       Implicit None
!!       Type(simulation_parameters_type), Intent (OUT):: mc_par
!!       Type(global_environment_type), Intent (IN):: mcpar
!!       Type(mcio_parameters_type), Intent (OUT):: mcio
!!     End Subroutine read_mc_section
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     MJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   NOTES
!!I---------------------------------------------------------------------------I
!!I SECTION: &mc ... &end                                                     I
!!I                                                                           I
!!I nsteps:       n                       total number of steps                I
!!I temperature: T                       temperature [K]                      I
!!I                                                                           I
!!I---------------------------------------------------------------------------I
!!*** ************************************************************************

SUBROUTINE read_mc_section ( mc_par, mcpar )

  IMPLICIT NONE

! Arguments
  TYPE ( mc_simulation_parameters_type ), INTENT ( OUT ) :: mc_par
  TYPE ( global_environment_type ), INTENT ( IN ) :: mcpar

! Locals
  INTEGER :: ierror, ilen, msglen, ia, ie, iw, allgrp, source
  INTEGER :: ibox,imolty
  CHARACTER ( LEN = 20 ) :: string
  CHARACTER ( LEN = 5 ) :: label

!------------------------------------------------------------------------------

!..defaults


  iw = mcpar % scr
  WRITE ( iw, * )
!..filenames

!..parse the input section
  label = '&MC'
  CALL parser_init(mcpar % input_file_name,mcpar)
  CALL search_label(label,ierror,ignore_case=.TRUE.)
  IF (ierror /= 0 ) THEN
     WRITE ( iw, *) ierror
     IF( mcpar % ionode ) &
        WRITE ( iw, '( a )' ) ' No input section &MC found '
  ELSE
     CALL read_line
     DO WHILE ( test_next() /= 'X' )
        ilen = 9
        CALL get_next ( string, ilen )
        CALL uppercase ( string )

        SELECT CASE (string)
        CASE DEFAULT
           CALL p_error()
           CALL stop_parser ( 'read_mc_section', 'unknown option' )

        CASE ( 'NSTEP' )
           CALL get_next ( mc_par % nstep )
        CASE ( 'LSTOP' )
           CALL get_next ( mc_par % lstop, ilen )
        CASE ( 'IUPTRAION' )
           CALL get_next ( mc_par % iuptraion )
        CASE ( 'IUPTRANS' )
           CALL get_next ( mc_par % iuptrans )
        CASE ( 'PROGRAM' )
           CALL get_next ( mc_par % ensemble, ilen )
        CASE ( 'TEMP' )
           CALL get_next ( mc_par % temp )
! figure out what beta is in atomic units (1/kT)
           mc_par % BETA = 1 / mc_par%temp / boltzmann * joule
        CASE ( 'IPRINT' )
           CALL get_next ( mc_par % iprint )
        CASE ( 'IBLOCK' )
           CALL get_next ( mc_par % iblock )
        CASE ( 'NBOX' )
           CALL get_next ( mc_par % nbox )
        CASE ( 'PMTRAION' )
           CALL get_next ( mc_par % pmtraion )
        CASE ( 'PMTRANS' )
           CALL get_next ( mc_par % pmtrans )
        CASE ( 'RMION' )
           CALL get_next ( mc_par % rmion )
        CASE ( 'RMROT' )
           CALL get_next ( mc_par % rmrot )
        CASE ( 'RMTRANS' )
           CALL get_next ( mc_par % rmtrans )
        CASE ( 'NCHAIN' )
           CALL get_next ( mc_par % nchain )
        CASE ( 'NUNITS' )
           CALL get_next ( mc_par % nunits )
        CASE ( 'NMOLTY' )
           CALL get_next ( mc_par % nmolty )
        CASE ( 'MOLTYP' )
           ALLOCATE(mc_par % moltyp(1:mc_par%nmolty))
           DO imolty=1,mc_par%nmolty
               CALL get_next ( mc_par % moltyp(imolty) )
           ENDDO
        CASE ( 'ININCH' )
           ALLOCATE(mc_par % ininch(1:mc_par%nbox))
           DO ibox=1,mc_par%nbox
               CALL get_next ( mc_par % ininch(ibox) )
           ENDDO
        CASE ( 'INIX' )
           ALLOCATE(mc_par % inix(1:mc_par%nbox))
           DO ibox=1,mc_par%nbox
               CALL get_next ( mc_par % inix(ibox) )
           ENDDO
        CASE ( 'INIY' )
           ALLOCATE(mc_par % iniy(1:mc_par%nbox))
           DO ibox=1,mc_par%nbox
               CALL get_next ( mc_par % iniy(ibox) )
           ENDDO
        CASE ( 'INIZ' )
           ALLOCATE(mc_par % iniz(1:mc_par%nbox))
           DO ibox=1,mc_par%nbox
               CALL get_next ( mc_par % iniz(ibox) )
           ENDDO
        CASE ( 'IRATV' )
           CALL get_next ( mc_par % iratv )
        CASE ( 'IRATP' )
           CALL get_next ( mc_par % iratp )
        END SELECT
        CALL read_line
     END DO

  END IF
  CALL parser_end
!..end of parsing the input section

! find the total number of units in the system
  mc_par%nunits_tot=mc_par%nunits*mc_par%nchain

! find out how many steps we're doing
  IF (mc_par%lstop .EQ. 'no') THEN
      mc_par%nstep=mc_par%nstep*mc_par%nunits_tot
  ENDIF

!..write some information to output
  IF (mcpar % ionode .AND. mcpar % print_level>=0) THEN
     WRITE ( iw, '( A )' ) ' MC| Monte Carlo Protocol '
     WRITE ( iw, '( A,T71,I10 )' ) ' MC| total number of steps ', &
          mc_par % nstep
     WRITE ( iw, '( A,T71,I10 )' ) ' MC| iprint ', &
          mc_par % iprint
     WRITE ( iw, '( A,T71,I10 )' ) ' MC| iblock ', &
          mc_par % iblock
     WRITE ( iw, '( A,T71,I10 )' ) ' MC| nbox ', &
          mc_par % nbox
     WRITE ( iw, '( A,T71,I10 )' ) ' MC| number of molecules', &
          mc_par % nchain
     WRITE ( iw, '( A,A,T71,I10 )' ) ' MC| number of interactions ',&
              'sites per molecule', mc_par % nunits
     WRITE ( iw, '( A,A,T71,I10 )' ) ' MC| total number of inter',&
              'action sites', mc_par % nunits_tot
     WRITE ( iw, '( A,T71,F10.3 )' ) ' MC| pmtraion ', &
          mc_par % pmtraion
     WRITE ( iw, '( A,T71,F10.3 )' ) ' MC| pmtrans ', &
          mc_par % pmtrans
     WRITE ( iw, '( A,T71,I10 )' ) ' MC| iuptraion ', &
          mc_par % iuptraion
     WRITE ( iw, '( A,T71,I10 )' ) ' MC| iuptrans ', &
          mc_par % iuptrans
     WRITE ( iw, '( A,T71,A20 )' ) ' MC| ensemble ', &
          mc_par % ensemble
     WRITE ( iw, '( A,T71,F10.3 )' ) ' MC| rmion ', &
          mc_par % rmion
     WRITE ( iw, '( A,T71,F10.3 )' ) ' MC| rmtrans ', &
          mc_par % rmtrans
     WRITE ( iw, '( A,T71,F10.3 )' ) ' MC| rmrot ', &
          mc_par % rmrot
     WRITE ( iw, '( A,T71,I10 )' ) ' MC| nmolty ', &
          mc_par % nmolty
     DO imolty=1,mc_par%nmolty
         WRITE( iw, '(A,I3,A,I2)' ) 'We have ',mc_par % moltyp(imolty), &
             ' molecules of type ', imolty
     ENDDO
     DO ibox=1,mc_par%nbox
         WRITE( iw, '(A,I1,A,I3)' ) 'The initial number of chains in box ',ibox, &
             ' is ', mc_par % ininch(ibox)
         WRITE( iw, '(2A,I2,A,3I2,2A)') 'The number of grid pts. in the initial ',&
             'simulation box ',ibox, ' are ',mc_par % inix(ibox),mc_par%iniy(ibox),&
             mc_par%iniz(ibox),' in the x, y, and z directions.'
     ENDDO
     WRITE ( iw, '( A,T71,I10 )' ) ' MC| iratv ', &
          mc_par % iratv
     WRITE ( iw, '( A,T71,I10 )' ) ' MC| iratp ', &
          mc_par % iratp
     WRITE ( iw, '( A,T71,F10.2 )' ) ' MC| Temperature [K] ', &
          mc_par % temp
  END IF

END SUBROUTINE read_mc_section

!*****************************************************************************

END MODULE mc_control

!******************************************************************************
