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
  USE cp_error_handling,               ONLY: cp_error_type
  USE cp_subsystem_types,              ONLY: cp_subsystem_p_type,&
                                             cp_subsys_get
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
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
  USE particle_list_types,             ONLY: particle_list_type
  USE physcon,                         ONLY: boltzmann,joule
  USE simulation_cell,                 ONLY: cell_type
  USE string_utilities,                ONLY: uppercase
  USE timings,                         ONLY: timeset,&
                                             timestop


  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_mc_section, mc_simulation_parameters_type,&
            write_mc_restart,read_mc_restart


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
     INTEGER :: nmoves
     INTEGER :: crd
     INTEGER :: ene
     INTEGER :: rm
     INTEGER :: cl
     INTEGER :: data
     INTEGER :: nstart
     LOGICAL :: lclassical
     LOGICAL :: restart
     CHARACTER ( LEN = 3 ) :: lstop
     CHARACTER ( LEN = 20 ) :: ensemble
     CHARACTER ( LEN = 200 ) :: restart_file
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
     REAL ( dp ) :: epsilon
     REAL ( dp ) :: sigma
     REAL ( dp ) :: charge
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
! units to write various data to
  mc_par % crd=75
  mc_par % rm=76
  mc_par % data=77
  mc_par % ene=78
  mc_par % cl=79
! force field parameters
  mc_par % sigma = 5.48d0 ! a.u., 2.9 angstroms
  mc_par % epsilon = 0.000247  ! Hartree, 78 Kelvin
  mc_par % charge = 1.122d-19  ! Coulombs, -0.7 e
! starting step
  mc_par % nstart = 0

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
        CASE ( 'NMOVES' )
          CALL get_next ( mc_par % nmoves )
        CASE ( 'CRD' )
          CALL get_next ( mc_par % crd )
        CASE ( 'DATA' )
          CALL get_next ( mc_par % data )
        CASE ( 'ENE' )
          CALL get_next ( mc_par % ene )
        CASE ( 'CL' )
          CALL get_next ( mc_par % cl )
        CASE ( 'RM' )
          CALL get_next ( mc_par % rm )
        CASE ( 'LSTOP' )
           CALL get_next ( mc_par % lstop, ilen )
        CASE ( 'LCLASSICAL' )
           CALL get_next ( string, ilen )
           CALL uppercase ( string )
           IF ( string(1:2)=='NO' ) mc_par % lclassical = .FALSE.
           IF ( string(1:3)=='YES' ) mc_par % lclassical = .TRUE.
        CASE ( 'RESTART' )
           CALL get_next ( string, ilen )
           CALL uppercase ( string )
           IF ( string(1:2)=='NO' ) mc_par % restart = .FALSE.
           IF ( string(1:3)=='YES' ) mc_par % restart = .TRUE.
        CASE ( 'IUPVOLUME' )
           CALL get_next ( mc_par % iupvolume )
        CASE ( 'IUPTRANS' )
           CALL get_next ( mc_par % iuptrans )
        CASE ( 'PROGRAM' )
           CALL get_next ( mc_par % ensemble, ilen )
        CASE ( 'RESTART_FILE' )
           CALL get_next ( mc_par % restart_file, ilen )
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
        END SELECT
        CALL read_line
     END DO

  END IF
  CALL parser_end
!..end of parsing the input section

! problems
      IF ((mc_par % nmoves .EQ. 1) .AND. mc_par % lclassical ) THEN
         STOP 'If lclassical=.TRUE., nmoves must be greater than 1.'
      ENDIF

! find out how many steps we're doing
!  IF (mc_par%lstop .EQ. 'no') THEN
!      mc_par%nstep=mc_par%nstep*mc_par%nunits_tot
!  ENDIF

!..write some information to output
  IF (globenv % ionode .AND. globenv % print_level>=0) THEN
     WRITE ( iw, '( A )' ) ' MC| Monte Carlo Protocol '
     WRITE ( iw, '( A,T71,I10 )' ) ' MC| total number of steps ', &
          mc_par % nstep
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
     WRITE ( iw, '( A,A,T78,I3 )' ) ' MC| Unit to write',&
              ' coordinates to', mc_par % crd
     WRITE ( iw, '( A,A,T78,I3 )' ) ' MC| Unit to write',&
              ' energies to', mc_par % ene
     WRITE ( iw, '( A,A,T78,I3 )' ) ' MC| Unit to write',&
              ' maximum move values to', mc_par % rm
     WRITE ( iw, '( A,A,T78,I3 )' ) ' MC| Unit to write',&
              ' molecule data to', mc_par % data
     WRITE ( iw, '( A,A,T78,I3 )' ) ' MC| Unit to write',&
              ' cell lengths to', mc_par % cl
     WRITE ( iw, '( A,T68,A )' ) ' MC| Name of restart file',&
                   mc_par % restart_file
     WRITE ( iw, '( A,T71,F10.3 )' ) ' MC| rmvolume [cubic a.u.]', &
          mc_par % rmvolume
     WRITE ( iw, '( A,T71,F10.3 )' ) ' MC| rmion [a.u.]', &
          mc_par % rmion
     WRITE ( iw, '( A,T71,F10.3 )' ) ' MC| rmtrans [a.u.]', &
          mc_par % rmtrans
     WRITE ( iw, '( A,T71,F10.3 )' ) ' MC| rmrot [radians]', &
          mc_par % rmrot
     WRITE ( iw, '( A,T71,F10.2 )' ) ' MC| Temperature [K] ', &
          mc_par % temp
     WRITE ( iw, '( A,T71,E10.5 )' ) ' MC| Pressure [a.u.] ', &
          mc_par % pressure
     IF ( mc_par % restart ) THEN
          WRITE ( iw, '(A,A)') ' MC| Initial data will be read from a',&
               ' restart file.'
     ENDIF
     IF ( mc_par % lclassical ) THEN
          WRITE ( iw, '(A,A)') ' MC| The moves will be biased with a',&
          ' classical force field.'
     ELSE
          WRITE ( iw, '(A,A)') ' MC| The moves will not be biased',&
          ' with a classical force field.'
     ENDIF
     IF (mc_par%nmoves .EQ. 1) THEN
          WRITE(iw,'(A,A)') ' MC| A Quickstep energy calculation ',&
          'will be done at every step.'
     ELSE
          WRITE( iw, '(A,I4,A,A)' ) ' MC| ',mc_par%nmoves,&
          ' moves will be attempted ',&
          'before a Quickstep energy calculation takes place.'
     ENDIF
  END IF

! end the timing
  CALL timestop(0.0d0,handle)


END SUBROUTINE read_mc_section

SUBROUTINE write_mc_restart ( nnstep, mc_par, force_env )

! write a bunch of data to a file that we can restart from

  IMPLICIT NONE

! Arguments
  TYPE ( mc_simulation_parameters_type ), POINTER :: mc_par
  TYPE ( force_env_type ), POINTER :: force_env
  INTEGER,INTENT(IN) :: nnstep

!   *** Local parameters ***

  CHARACTER(LEN=*), PARAMETER :: routineN = "write_mc_restart"

! Locals
      INTEGER :: ipart,handle
      TYPE(cell_type),POINTER :: cell
      TYPE(cp_subsystem_p_type),DIMENSION(:),POINTER :: subsys
      TYPE(particle_list_type), POINTER :: particles
      TYPE(cp_error_type)      :: error
!------------------------------------------------------------------------------

! begin the timing of the subroutine
  CALL timeset(routineN,'I','',handle)

! open the file and write some simulation parameters
  OPEN(UNIT=1,FILE=mc_par%restart_file)
  WRITE(1,*) nnstep
  WRITE(1,*) mc_par%temp,mc_par%nunits_tot
  WRITE(1,*) mc_par%ensemble

! get the cell length and coordinates
  CALL force_env_get(force_env,cell=cell,subsys=subsys,&
                         error=error) 
  CALL cp_subsys_get(subsys(1)%subsys, &
             particles=particles, error=error)

  WRITE(1,'(F8.4)') cell%hmat(1,1)*0.529177d0  ! in angstoms
  WRITE(1,*)
  DO ipart=1,mc_par%nunits_tot
      WRITE(1,'(1X,3(F15.10,3X))') &
                 particles%els(ipart)%r(1:3)*0.529177d0
  ENDDO
  CLOSE(1)

! end the timing
  CALL timestop(0.0d0,handle)


END SUBROUTINE write_mc_restart

!*****************************************************************************

SUBROUTINE read_mc_restart ( mc_par, force_env ,iw )

! read a bunch of data from a restart file

  IMPLICIT NONE

! Arguments
  TYPE ( mc_simulation_parameters_type ), POINTER :: mc_par
  TYPE ( force_env_type ), POINTER :: force_env
  INTEGER,INTENT(IN) :: iw

!   *** Local parameters ***

  CHARACTER(LEN=*), PARAMETER :: routineN = "read_mc_restart"

! Locals
      INTEGER :: ipart,handle,nunits_tot
      TYPE(cell_type),POINTER :: cell
      TYPE(cp_subsystem_p_type),DIMENSION(:),POINTER :: subsys
      TYPE(particle_list_type), POINTER :: particles
      TYPE(cp_error_type)      :: error
      REAL(KIND=dp),DIMENSION(1:3) :: r
      REAL(KIND=dp) :: temp,box_length
      CHARACTER(LEN=20) :: ensemble
!------------------------------------------------------------------------------

! begin the timing of the subroutine
  CALL timeset(routineN,'I','',handle)

! open the file and write some simulation parameters
  OPEN(UNIT=1,FILE=mc_par%restart_file)
  READ(1,*) mc_par%nstart
  READ(1,*) temp,nunits_tot
  READ(1,*) ensemble

! do some checking
  IF(ABS(temp - mc_par%temp) .GT. 0.01d0) THEN
      WRITE(iw,*) 'The temperature in the restart file is ',&
            'not the same as the input file.'
      WRITE(iw,*) 'Input file temp =',mc_par%temp
      WRITE(iw,*) 'Restart file temp =',temp
      STOP
  ENDIF
  IF(nunits_tot .NE. mc_par%nunits_tot) THEN
      WRITE(iw,*) 'The total number of units in the restart file is ',&
            'not the same as the input file.'
      WRITE(iw,*) 'Input file units =',mc_par%nunits_tot
      WRITE(iw,*) 'Restart file units =',nunits_tot
      STOP
  ENDIF
  IF(ensemble .NE. mc_par%ensemble) THEN
      WRITE(iw,*) 'The ensemble in the restart file is ',&
            'not the same as the input file.'
      WRITE(iw,*) 'Input file ensemble =',mc_par%ensemble
      WRITE(iw,*) 'Restart file ensemble =',ensemble
      STOP
  ENDIF

! get the cell length and coordinates
  CALL force_env_get(force_env,cell=cell,subsys=subsys,&
                         error=error) 
  CALL cp_subsys_get(subsys(1)%subsys, &
             particles=particles, error=error)

  READ(1,*)  box_length ! in angstoms
  IF((ABS(box_length - cell%hmat(1,1)*0.529177d0) .GT. 0.0001d0) &
               .AND. (mc_par%pmvolume .EQ. 0.0d0)) THEN
      WRITE(iw,*) 'The cell length in the restart file is ',&
            'not the same as the input file.'
      WRITE(iw,*) 'Input file cell length =',cell%hmat(1,1)*0.529177d0
      WRITE(iw,*) 'Restart file cell length =',box_length
      STOP
  ENDIF
  READ(1,*)
! read in the particle coordinates
  DO ipart=1,mc_par%nunits_tot
      READ(1,*) r(1:3)
      particles%els(ipart)%r(1:3)=r(1:3)/0.529177d0
  ENDDO
  CLOSE(1)

! end the timing
  CALL timestop(0.0d0,handle)


END SUBROUTINE read_mc_restart

!*****************************************************************************

END MODULE mc_control

!******************************************************************************
