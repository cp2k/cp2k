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
  USE convert_units,                   ONLY: convert_to_cp2k_units
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
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
  USE timings,                         ONLY: timeset,&
                                             timestop


  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_mc_section, mc_simulation_parameters_type


!!***
!****************************************************************************

  TYPE mc_simulation_parameters_type
     INTEGER :: nstep
     INTEGER :: iprint
     INTEGER :: iupvolume
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
     CHARACTER ( LEN = 200 ) :: box2_file
     CHARACTER ( LEN = 50 ) :: program
     REAL ( dp ) :: rmvolume
     REAL ( dp ) :: rmion
     REAL ( dp ) :: rmrot
     REAL ( dp ) :: rmtrans
     REAL ( dp ) :: temp
     REAL ( dp ) :: pressure
     REAL ( dp ) :: pmvolume
     REAL ( dp ) :: pmtraion
     REAL ( dp ) :: pmtrans
     REAL ( dp ) :: BETA
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
!!     Subroutine read_mc_section(mc_par, globenv)
!!       Implicit None
!!       Type(simulation_parameters_type), Intent (OUT):: mc_par
!!       Type(global_environment_type), Intent (IN):: globenv
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

SUBROUTINE read_mc_section ( mc_par, globenv )

  IMPLICIT NONE

! Arguments
  TYPE ( mc_simulation_parameters_type ), INTENT ( OUT ) :: mc_par
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv

!   *** Local parameters ***

  CHARACTER(LEN=*), PARAMETER :: routineN = "read_mc_section"

! Locals
  INTEGER :: ierror, ilen, msglen, ia, ie, iw, allgrp, source
  INTEGER :: ibox,imolty,handle
  CHARACTER ( LEN = 20 ) :: string
  CHARACTER ( LEN = 5 ) :: label
  REAL( KIND = dp ) :: length

!------------------------------------------------------------------------------

! begin the timing of the subroutine
  CALL timeset(routineN,'I','',handle)

!..defaults

  iw = globenv % scr
  WRITE ( iw, * )
!..filenames
!..parse the input section
  label = '&MC'
  CALL parser_init(globenv % input_file_name,globenv)
  CALL search_label(label,ierror,ignore_case=.TRUE.)
  IF (ierror /= 0 ) THEN
     WRITE ( iw, *) ierror
     IF( globenv % ionode ) &
        WRITE ( iw, '( a )' ) ' No input section &MC found '
  ELSE
     CALL read_line
     DO WHILE ( test_next() /= 'X' )
        ilen = 20
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
        CASE ( 'IUPVOLUME' )
           CALL get_next ( mc_par % iupvolume )
        CASE ( 'IUPTRANS' )
           CALL get_next ( mc_par % iuptrans )
        CASE ( 'PROGRAM' )
           CALL get_next ( mc_par % ensemble, ilen )
        CASE ( 'BOX2_FILE' )
           CALL get_next ( mc_par % box2_file, ilen )
        CASE ( 'PRESSURE' )
           CALL get_next ( mc_par % pressure )
! convert from bar to a.u.
           CALL convert_to_cp2k_units(pressure=mc_par%pressure)
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
        CASE ( 'PMVOLUME' )
           CALL get_next ( mc_par % pmvolume )
        CASE ( 'PMTRAION' )
           CALL get_next ( mc_par % pmtraion )
        CASE ( 'PMTRANS' )
           CALL get_next ( mc_par % pmtrans )
        CASE ( 'RMVOLUME' )
           CALL get_next ( mc_par % rmvolume )
! convert from angstrom to a.u.
           length=mc_par%rmvolume**0.333333333333333d0
           CALL convert_to_cp2k_units(unit_char='ANGSTROM',&
                            length=length)
           mc_par%rmvolume=length**3
        CASE ( 'RMION' )
           CALL get_next ( mc_par % rmion )
! convert from angstrom to a.u.
           CALL convert_to_cp2k_units(unit_char='ANGSTROM',&
                            length=mc_par%rmion)
        CASE ( 'RMROT' )
           CALL get_next ( mc_par % rmrot )
        CASE ( 'RMTRANS' )
           CALL get_next ( mc_par % rmtrans )
! convert from angstrom to a.u.
           CALL convert_to_cp2k_units(unit_char='ANGSTROM',&
                            length=mc_par%rmtrans)
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
!  IF (mc_par%lstop .EQ. 'no') THEN
!      mc_par%nstep=mc_par%nstep*mc_par%nunits_tot
!  ENDIF

!..write some information to output
  IF (globenv % ionode .AND. globenv % print_level>=0) THEN
     WRITE ( iw, '( A )' ) ' MC| Monte Carlo Protocol '
     WRITE ( iw, '( A,T71,I10 )' ) ' MC| total number of steps ', &
          mc_par % nstep
!     WRITE ( iw, '( A,T71,I10 )' ) ' MC| iprint ', &
!          mc_par % iprint
!     WRITE ( iw, '( A,T71,I10 )' ) ' MC| iblock ', &
!          mc_par % iblock
!     WRITE ( iw, '( A,T71,I10 )' ) ' MC| nbox ', &
!          mc_par % nbox
     WRITE ( iw, '( A,T71,I10 )' ) ' MC| number of molecules', &
          mc_par % nchain
     WRITE ( iw, '( A,A,T71,I10 )' ) ' MC| number of interactions ',&
              'sites per molecule', mc_par % nunits
     WRITE ( iw, '( A,A,T71,I10 )' ) ' MC| total number of inter',&
              'action sites', mc_par % nunits_tot
     WRITE ( iw, '( A,T71,F10.3 )' ) ' MC| pmvolume ', &
          mc_par % pmvolume
     WRITE ( iw, '( A,T71,F10.3 )' ) ' MC| pmtraion ', &
          mc_par % pmtraion
     WRITE ( iw, '( A,T71,F10.3 )' ) ' MC| pmtrans ', &
          mc_par % pmtrans
     WRITE ( iw, '( A,T71,I10 )' ) ' MC| iupvolume ', &
          mc_par % iupvolume
     WRITE ( iw, '( A,T71,I10 )' ) ' MC| iuptrans ', &
          mc_par % iuptrans
     WRITE ( iw, '( A,T68,A20 )' ) ' MC| ensemble ', &
          mc_par % ensemble
     IF (mc_par%ensemble .EQ. 'GEMC-NVT') THEN
         WRITE ( iw, '( A,T68,A20 )' ) ' MC| Box 2 file', &
              mc_par % box2_file
     ENDIF
     WRITE ( iw, '( A,T71,F10.3 )' ) ' MC| rmvolume [cubic a.u.]', &
          mc_par % rmvolume
     WRITE ( iw, '( A,T71,F10.3 )' ) ' MC| rmion [a.u.]', &
          mc_par % rmion
     WRITE ( iw, '( A,T71,F10.3 )' ) ' MC| rmtrans [a.u.]', &
          mc_par % rmtrans
     WRITE ( iw, '( A,T71,F10.3 )' ) ' MC| rmrot [radians]', &
          mc_par % rmrot
!     DO imolty=1,mc_par%nmolty
!         WRITE( iw, '(A,I3,A,I2)' ) 'We have ',mc_par % moltyp(imolty), &
!             ' molecules of type ', imolty
!     ENDDO
!     DO ibox=1,mc_par%nbox
!         WRITE( iw, '(A,I1,A,I3)' ) 'The initial number of chains in box ',ibox, &
!             ' is ', mc_par % ininch(ibox)
!         WRITE( iw, '(2A,I2,A,3I2,2A)') 'The number of grid pts. in the initial ',&
!             'simulation box ',ibox, ' are ',mc_par % inix(ibox),mc_par%iniy(ibox),&
!             mc_par%iniz(ibox),' in the x, y, and z directions.'
!     ENDDO
!     WRITE ( iw, '( A,T71,I10 )' ) ' MC| iratv ', &
!          mc_par % iratv
!     WRITE ( iw, '( A,T71,I10 )' ) ' MC| iratp ', &
!          mc_par % iratp
     WRITE ( iw, '( A,T71,F10.2 )' ) ' MC| Temperature [K] ', &
          mc_par % temp
     WRITE ( iw, '( A,T71,E10.5 )' ) ' MC| Pressure [a.u.] ', &
          mc_par % pressure
  END IF

! end the timing
  CALL timestop(0.0d0,handle)


END SUBROUTINE read_mc_section

!*****************************************************************************

END MODULE mc_control

!******************************************************************************
