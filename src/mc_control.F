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
  USE cp_subsystem_types,              ONLY: cp_subsys_get,&
                                             cp_subsystem_p_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: pi
  USE parser,                          ONLY: get_next,&
                                             p_error,&
                                             parser_end,&
                                             parser_init,&
                                             read_line,&
                                             search_label,&
                                             stop_parser,&
                                             test_next
  USE particle_list_types,             ONLY: particle_list_type
  USE physcon,                         ONLY: boltzmann,&
                                             joule,&
                                             angstrom
  USE qs_parser,                       ONLY: close_file,&
                                             open_file
  USE simulation_cell,                 ONLY: cell_type
  USE string_utilities,                ONLY: uppercase,&
                                             xstring
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE 
  PUBLIC :: read_mc_section, mc_simulation_parameters_type,&
            write_mc_restart,read_mc_restart


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
     INTEGER :: nvolmoves
     INTEGER :: crd
     INTEGER :: ene
     INTEGER :: rm
     INTEGER :: cl
     INTEGER :: data
     INTEGER :: nstart
     INTEGER :: source
     INTEGER :: group
     LOGICAL :: lclassical
     LOGICAL :: ionode
     LOGICAL :: restart
     CHARACTER ( LEN = 3 ) :: lstop
     CHARACTER ( LEN = 20 ) :: ensemble
     CHARACTER ( LEN = 30 ) :: restart_file
     CHARACTER ( LEN = 30 ) :: coords_file
     CHARACTER ( LEN = 30 ) :: energy_file
     CHARACTER ( LEN = 30 ) :: displacement_file
     CHARACTER ( LEN = 30 ) :: cell_file
     CHARACTER ( LEN = 30 ) :: data_file
     CHARACTER ( LEN = 200 ) :: box2_file
     CHARACTER ( LEN = 50 ) :: program
     REAL ( dp ) :: rmvolume
     REAL ( dp ) :: rmbond
     REAL ( dp ) :: rmangle
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


    TYPE(mc_simulation_parameters_type), &
      INTENT(OUT)                            :: mc_par
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=*), PARAMETER     :: routineN = "read_mc_section"

    CHARACTER(LEN=20)                        :: string
    CHARACTER(LEN=5)                         :: label
    INTEGER                                  :: handle,ierror,ilen, iw,&
                                                ia,ie
    REAL(KIND=dp)                            :: length

!------------------------------------------------------------------------------
! begin the timing of the subroutine

  CALL timeset(routineN,'I','',handle)

! set the ionode and mepos
  mc_par % ionode = globenv % ionode
  mc_par % group = globenv % group
  mc_par % source = globenv % source

!..defaults
  mc_par % nstart = 0
  mc_par % nvolmoves = 4

  iw = globenv % scr
  WRITE ( iw, * )

!..filenames
  CALL xstring ( mc_par % program, ia, ie )
  mc_par%coords_file = mc_par % program(ia:ie) // '.coordinates'
  mc_par%energy_file = mc_par % program(ia:ie) // '.energy'
  mc_par%cell_file = mc_par % program(ia:ie) // '.cell'
  mc_par%displacement_file= mc_par % program(ia:ie)&
      // '.max_displacements'
  mc_par%data_file = mc_par % program(ia:ie) // '.data'

!..parse the input section
  label = '&MC'
  CALL parser_init(globenv % input_file_name,globenv)
  CALL search_label(label,ierror,ignore_case=.TRUE.)
  IF (ierror /= 0 ) THEN
     WRITE ( iw, *) ierror
     IF( mc_par % ionode ) &
        WRITE ( iw, '( a )' ) ' No input section &MC found '
  ELSE
     CALL read_line
     DO WHILE ( test_next() /= 'X' )
        ilen = 40
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
        CASE ( 'NVOLMOVES' )
          CALL get_next ( mc_par % nvolmoves )
        CASE ( 'LSTOP' )
           CALL get_next ( mc_par % lstop, ilen )
        CASE ( 'LCLASSICAL' )
           CALL get_next ( string, ilen )
           CALL uppercase ( string )
           IF ( string(1:2)=='NO' ) THEN
              mc_par % lclassical = .FALSE.
              mc_par % nvolmoves = 1
           ENDIF
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
        CASE ( 'COORDINATE_FILE' )
           CALL get_next ( mc_par % coords_file, ilen )
        CASE ( 'ENERGY_FILE' )
           CALL get_next ( mc_par % energy_file, ilen )
        CASE ( 'DATA_FILE' )
           CALL get_next ( mc_par % data_file, ilen )
        CASE ( 'CELL_FILE' )
           CALL get_next ( mc_par % cell_file, ilen )
        CASE ( 'MAX_DISP_FILE' )
           CALL get_next ( mc_par % displacement_file, ilen )
        CASE ( 'BOX2_FILE' )
           CALL get_next ( mc_par % box2_file, ilen )
        CASE ( 'PRESSURE' )
           CALL get_next ( mc_par % pressure )
        CASE ( 'TEMP' )
           CALL get_next ( mc_par % temp )
        CASE ( 'PMVOLUME' )
           CALL get_next ( mc_par % pmvolume )
        CASE ( 'PMTRAION' )
           CALL get_next ( mc_par % pmtraion )
        CASE ( 'PMTRANS' )
           CALL get_next ( mc_par % pmtrans )
        CASE ( 'RMVOLUME' )
           CALL get_next ( mc_par % rmvolume )
        CASE ( 'RMBOND' )
           CALL get_next ( mc_par % rmbond )
        CASE ( 'RMANGLE' )
           CALL get_next ( mc_par % rmangle )
        CASE ( 'RMROT' )
           CALL get_next ( mc_par % rmrot )
        CASE ( 'RMTRANS' )
           CALL get_next ( mc_par % rmtrans )
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

!..write some information to output
  IF (mc_par % ionode .AND. globenv % print_level>=0) THEN
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
     WRITE ( iw, '( A,T58,A20 )' ) ' MC| ensemble ', &
          mc_par % ensemble
     IF (mc_par%ensemble .EQ. 'GEMC-NVT') THEN
         WRITE ( iw, '( A,T58,A)' ) ' MC| Box 2 file', &
                TRIM(mc_par % box2_file)
     ENDIF
     WRITE ( iw, '( A,T58,A )' ) ' MC| Name of restart file:',&
                   TRIM(mc_par % restart_file)
     WRITE ( iw, '( A,T44,A )' ) ' MC| Name of output coordinate file:',&
                   TRIM(mc_par % coords_file)
     WRITE ( iw, '( A,T44,A )' ) ' MC| Name of output data file:',&
                   TRIM(mc_par % data_file)
     WRITE ( iw, '( A,T44,A )' ) ' MC| Name of output energy file:',&
                   TRIM(mc_par % energy_file)
     WRITE ( iw, '( A,T44,A )' ) ' MC| Name of output cell file:',&
                   TRIM(mc_par % cell_file)
     WRITE ( iw, '( A,A,T44,A )' ) ' MC| Name of output',&
                   ' displacement file:',&
                   TRIM(mc_par % displacement_file)
     WRITE ( iw, '( A,T71,F10.3 )' ) ' MC| rmvolume [cubic angstroms]', &
          mc_par % rmvolume
     WRITE ( iw, '( A,T71,F10.3 )' ) ' MC| rmbond [angstroms]', &
          mc_par % rmbond
     WRITE ( iw, '( A,T71,F10.3 )' ) ' MC| rmangle [degrees]', &
          mc_par % rmangle
     WRITE ( iw, '( A,T71,F10.3 )' ) ' MC| rmtrans [angstroms]', &
          mc_par % rmtrans
     WRITE ( iw, '( A,T71,F10.3 )' ) ' MC| rmrot [degrees]', &
          mc_par % rmrot
     WRITE ( iw, '( A,T71,F10.2 )' ) ' MC| Temperature [K] ', &
          mc_par % temp
     WRITE ( iw, '( A,T71,F10.5 )' ) ' MC| Pressure [bar] ', &
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
          'before a Quickstep energy calculation'
          WRITE( iw, '(A)' ) ' MC|      takes place.'
     ENDIF
     IF (mc_par%pmvolume .GT. 0.0d0 ) THEN
          WRITE( iw, '(A,I4,A,A)') ' MC| ',mc_par%nvolmoves,&
          ' volume moves will be attempted ',&
          'before a Quickstep energy calculation'
          WRITE( iw, '(A)' ) ' MC|      takes place.'
     ENDIF
  END IF

! figure out what beta (1/kT) is in atomic units (1/Hartree)
      mc_par % BETA = 1 / mc_par%temp / boltzmann * joule
! convert from degrees to radians
      mc_par%rmrot = mc_par%rmrot/180.0d0*pi
! convert from bar to a.u.
      CALL convert_to_cp2k_units(pressure=mc_par%pressure)
! convert from angstrom to a.u.
      CALL convert_to_cp2k_units(unit_char='ANGSTROM',&
                     length=mc_par%rmtrans)
! convert from degrees to radians
      mc_par%rmangle = mc_par%rmangle/180.0d0*pi
! convert from angstrom to a.u.
      CALL convert_to_cp2k_units(unit_char='ANGSTROM',&
                    length=mc_par%rmbond)
! convert from angstrom to a.u.
      CALL convert_to_cp2k_units(unit_char='ANGSTROM',&
                            length=mc_par%rmvolume,l_power=3)

! end the timing
  CALL timestop(0.0d0,handle)


END SUBROUTINE read_mc_section

SUBROUTINE write_mc_restart ( nnstep, mc_par, force_env )

! write a bunch of data to a file that we can restart from


    INTEGER, INTENT(IN)                      :: nnstep
    TYPE(mc_simulation_parameters_type), &
      POINTER                                :: mc_par
    TYPE(force_env_type), POINTER            :: force_env

    CHARACTER(LEN=*), PARAMETER        :: routineN = "write_mc_restart"

    INTEGER                                  :: handle, ipart, unit
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type)                      :: error
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(particle_list_type), POINTER        :: particles

!------------------------------------------------------------------------------
! begin the timing of the subroutine

  CALL timeset(routineN,'I','',handle)

! open the file and write some simulation parameters
  CALL open_file(file_name=mc_par%restart_file,unit_number=unit,&
                 file_action='WRITE',file_status='OLD')
  WRITE(unit,*) nnstep
  WRITE(unit,*) mc_par%temp,mc_par%nunits_tot
  WRITE(unit,*) mc_par%ensemble

! get the cell length and coordinates
  CALL force_env_get(force_env,cell=cell,subsys=subsys,&
                         error=error) 
  CALL cp_subsys_get(subsys(1)%subsys, &
             particles=particles, error=error)

  WRITE(unit,'(F10.6)') cell%hmat(1,1)*angstrom  ! in angstroms
  WRITE(unit,*)
  DO ipart=1,mc_par%nunits_tot
      WRITE(unit,'(1X,3(F15.10,3X))') &
                 particles%els(ipart)%r(1:3)*angstrom
  ENDDO
  CLOSE(unit)

! end the timing
  CALL timestop(0.0d0,handle)


END SUBROUTINE write_mc_restart

!*****************************************************************************

SUBROUTINE read_mc_restart ( mc_par, force_env ,iw )

! read a bunch of data from a restart file


    TYPE(mc_simulation_parameters_type), &
      POINTER                                :: mc_par
    TYPE(force_env_type), POINTER            :: force_env
    INTEGER, INTENT(IN)                      :: iw

    CHARACTER(LEN=*), PARAMETER          :: routineN = "read_mc_restart"

    CHARACTER(LEN=20)                        :: ensemble
    INTEGER                              :: handle, ipart, nunits_tot,&
                                                unit
    REAL(KIND=dp)                            :: box_length, temp
    REAL(KIND=dp), DIMENSION(1:3)            :: r
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type)                      :: error
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(particle_list_type), POINTER        :: particles

!------------------------------------------------------------------------------
! begin the timing of the subroutine

  CALL timeset(routineN,'I','',handle)

! open the file and write some simulation parameters
  CALL open_file(file_name=mc_par%restart_file,unit_number=unit)
!  OPEN(UNIT=1,FILE=mc_par%restart_file)
  READ(unit,*) mc_par%nstart
  READ(unit,*) temp,nunits_tot
  READ(unit,*) ensemble

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

  READ(unit,*)  box_length ! in angstroms
  IF(ABS(box_length - cell%hmat(1,1)*angstrom) .GT. 0.0001d0 ) THEN
      WRITE(iw,*) 'The cell length in the restart file is ',&
            'not the same as the input file.  Change the restart file.'
      WRITE(iw,*) 'Input file cell length =',cell%hmat(1,1)*angstrom
      WRITE(iw,*) 'Restart file cell length =',box_length
      STOP
  ENDIF

  READ(unit,*)
! read in the particle coordinates
  DO ipart=1,mc_par%nunits_tot
      READ(unit,*) r(1:3)
      particles%els(ipart)%r(1:3)=r(1:3)/angstrom
  ENDDO
  CLOSE(unit)

! end the timing
  CALL timestop(0.0d0,handle)


END SUBROUTINE read_mc_restart


!---------------------------------
!      FUNCTION RANDOM() RESULT(rand)

!      IMPLICIT NONE
!      REAL(KIND=dp) :: rand 
!      RANDOM = RCARRY()
!      RCARRY=5.0d0
!      rand = RCARRY()


!      END FUNCTION RANDOM
!---------------------------------
!      INTEGER FUNCTION RANDx(ISEED)
!----------------------------------------------------------------------C
!  Random number generator, fast and rough, machine independent.
!  Returns an uniformly distributed deviate in the 0 to 1 interval.
!  This random number generator is portable, machine-independent and
!  reproducible, for any machine with at least 32 bits / real number.
!  REF: Press, Flannery, Teukolsky, Vetterling, Numerical Recipes (1986)
!----------------------------------------------------------------------C
!      IMPLICIT NONE
!      INTEGER IA , IC , ISEED , M1
!      double precision RANDx , RM
!      PARAMETER (M1=714025,IA=1366,IC=150889,RM=1.D+0/M1)
!
!      ISEED = MOD(IA*ISEED+IC,M1)
!      RANDx= ISEED*RM
!      IF ( RANDx.LT.0.D+0 ) STOP '*** Random number is negative ***'
!c
!      RETURN
!      END
!
!      SUBROUTINE RANSET(ISEED)
!      IMPLICIT NONE
!      INTEGER ISEED
!
!      CALL RSTART(ISEED)
!      RETURN
!      END
!
!      SUBROUTINE RSTART(ISEEDA)
!C----------------------------------------------------------------------C
!C       Initialize Marsaglia list of 24 random numbers.
!C----------------------------------------------------------------------C
!      IMPLICIT NONE
!      double precision Carry , ran , RANDx , Seed
!      INTEGER i , I24 , Iseed , ISEEDA , J24
!      COMMON /RANDM/ Seed(24) , Carry , I24 , J24 , Iseed

!      I24 = 24
!      J24 = 10
!      Carry = 0.D+0
!      Iseed = ISEEDA
!c
!c       get rid of initial correlations in rand by throwing
!c       away the first 100 random numbers generated.
!c
!      DO 10 i = 1 , 100
!        ran = RANDx(Iseed)
!   10 CONTINUE
!c
!c       initialize the 24 elements of seed
!c

!      DO 20 i = 1 , 24
!        Seed(i) = RANDx(Iseed)
!   20 CONTINUE
!
!      RETURN
!      END


!      FUNCTION RCARRY() RESULT(rcarryout)
!----------------------------------------------------------------------C
!       Random number generator from Marsaglia.
!----------------------------------------------------------------------C
!      IMPLICIT NONE
!      REAL(KIND=dp) ::  Carry , RCARRY , Seed , TWOM24 , TWOP24 , uni
!      REAL(KIND=dp),PARAMETER :: TWOP24=16777216.0D0,TWOM24=1.0D0/TWOP24
!      INTEGER       ::  I24 , Iseed , J24

!      COMMON /RANDM/ Seed(24) , Carry , I24 , J24 , Iseed

!       f.james Comp. Phys. Comm. 60, 329  (1990)
!       algorithm by G. Marsaglia and A. Zaman
!       base b = 2**24  lags r=24 and s=10

!      uni = Seed(J24) - Seed(I24) - Carry
!      IF ( uni.LT.0.D+0 ) THEN
!        uni = uni + 1.D+0
!        Carry = TWOM24
!      ELSE
!        Carry = 0.D+0
!      ENDIF
!      Seed(I24) = uni
!      I24 = I24 - 1
!      IF ( I24.EQ.0 ) I24 = 24
!      J24 = J24 - 1
!      IF ( J24.EQ.0 ) J24 = 24
!      rcarryout = uni

!      END FUNCTION RCARRY


!*****************************************************************************

END MODULE mc_control

!******************************************************************************
