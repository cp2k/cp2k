!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** mc_types.F *********************************************************
!!
!!   NAME
!!     mc_types
!!
!!   FUNCTION
!!     holds all the structure types needed for Monte Carlo, except
!!     the mc_environment_type
!!
!!   AUTHOR
!!     MJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE mc_types
  USE convert_units,                   ONLY: convert_to_cp2k_units
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
  USE physcon,                         ONLY: angstrom,&
                                             boltzmann,&
                                             joule
  USE string_utilities,                ONLY: xstring,&
                                             uppercase
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE 
! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "mc_types"

 PUBLIC :: mc_simulation_parameters_type,&
            mc_simulation_parameters_p_type,&
            mc_averages_type,mc_averages_p_type,&
            mc_moves_type, mc_moves_p_type,accattempt,&
            get_mc_par,set_mc_par,read_mc_section


  TYPE mc_simulation_parameters_type
      PRIVATE
      INTEGER :: nstep
      INTEGER :: iupvolume
      INTEGER :: iuptrans
      INTEGER :: nbox
      INTEGER :: nchain
      INTEGER :: nunits
      INTEGER :: nunits_tot
      INTEGER :: nmoves
      INTEGER :: nvolmoves
      INTEGER :: nswapmoves
      INTEGER :: rm
      INTEGER :: cl
      INTEGER :: diff
      INTEGER :: nstart
      INTEGER :: source
      INTEGER :: group
      LOGICAL :: lclassical
      LOGICAL :: lclassical_run
      LOGICAL :: ionode
      LOGICAL :: lrestart
      LOGICAL :: lstop
      CHARACTER ( LEN = 20 ) :: ensemble
      CHARACTER ( LEN = 40 ) :: restart_file
      CHARACTER ( LEN = 40 ) :: molecules_file
      CHARACTER ( LEN = 40 ) :: moves_file
      CHARACTER ( LEN = 40 ) :: coords_file
      CHARACTER ( LEN = 40 ) :: energy_file
      CHARACTER ( LEN = 40 ) :: displacement_file
      CHARACTER ( LEN = 40 ) :: cell_file
      CHARACTER ( LEN = 40 ) :: dat_file
      CHARACTER ( LEN = 40 ) :: data_file
      CHARACTER ( LEN = 200 ) :: box2_file
      CHARACTER ( LEN = 50 ) :: program
      REAL ( dp ) :: rmvolume
      REAL ( dp ) :: rmbond
      REAL ( dp ) :: rmangle
      REAL ( dp ) :: rmrot
      REAL ( dp ) :: rmtrans
      REAL ( dp ) :: temp
      REAL ( dp ) :: pressure
      REAL ( dp ) :: pmswap
      REAL ( dp ) :: pmvolume
      REAL ( dp ) :: pmtraion
      REAL ( dp ) :: pmtrans
      REAL ( dp ) :: BETA
      REAL ( dp ) :: cutoff
      REAL ( dp ) :: inner_cutoff
  END TYPE mc_simulation_parameters_type

  TYPE mc_simulation_parameters_p_type
      TYPE (mc_simulation_parameters_type),POINTER :: mc_par
  END TYPE mc_simulation_parameters_p_type

  TYPE mc_averages_type
      REAL(KIND = dp) :: ave_energy
      REAL(KIND = dp) :: ave_energy_squared
      REAL(KIND = dp) :: ave_cell_length
      REAL(KIND = dp) :: bond1
      REAL(KIND = dp) :: bond2
      REAL(KIND = dp) :: angle
      REAL(KIND = dp) :: molecules
  END TYPE mc_averages_type

  TYPE mc_averages_p_type
      TYPE (mc_averages_type),POINTER :: averages
  END TYPE mc_averages_p_type

  TYPE mc_moves_type
      TYPE ( accattempt ), POINTER :: classical_bond
      TYPE ( accattempt ), POINTER :: classical_angle
      TYPE ( accattempt ), POINTER :: classical_trans
      TYPE ( accattempt ), POINTER :: classical_rot
      TYPE ( accattempt ), POINTER :: classical_volume
      TYPE ( accattempt ), POINTER :: bond
      TYPE ( accattempt ), POINTER :: angle
      TYPE ( accattempt ), POINTER :: trans
      TYPE ( accattempt ), POINTER :: rot
      TYPE ( accattempt ), POINTER :: swap
      TYPE ( accattempt ), POINTER :: volume
      TYPE ( accattempt ), POINTER :: Quickstep
      REAL(KIND = dp) :: classical_energy,trans_dis,qtrans_dis
      integer :: empty,grown
  END TYPE mc_moves_type

  TYPE accattempt
      INTEGER :: successes
      INTEGER :: qsuccesses
      INTEGER :: attempts
  END TYPE accattempt

  TYPE mc_moves_p_type
      TYPE(mc_moves_type), POINTER :: moves
  END TYPE mc_moves_p_type

CONTAINS
!****************************************************************************

 SUBROUTINE get_mc_par ( mc_par, nstep, iuptrans, iupvolume, nchain,&
       nunits,nunits_tot,nmoves,nvolmoves,nswapmoves,rm,cl,diff,nstart,&
       source,group,lclassical,ionode,lrestart,lstop,rmvolume,rmbond,rmangle,&
       rmrot,rmtrans,temp,pressure,BETA,pmswap,pmvolume,pmtraion,pmtrans,&
       ensemble,program,restart_file,molecules_file,moves_file,coords_file,&
       energy_file,displacement_file,cell_file,dat_file,data_file,box2_file,&
       lclassical_run,cutoff,inner_cutoff)

! *****************************************************************************
!!****s* mc_types.F/get_mc_par
!!
!!   NAME
!!     get_mc_par
!!
!!   FUNCTION
!!     accesses the private elements of the mc_parameters_type
!!
!!   ARGUMENTS
!!     - mc_par: the structure mc parameters you want 
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************

    TYPE(mc_simulation_parameters_type), POINTER       :: mc_par
    INTEGER,INTENT(OUT),OPTIONAL :: nstep,iuptrans,iupvolume,nchain,nunits,&
      nunits_tot,nmoves,nvolmoves,nswapmoves,rm,cl,diff,nstart,source,&
      group

    LOGICAL,INTENT(OUT),OPTIONAL     :: lclassical,ionode,lrestart,lstop,&
      lclassical_run

    REAL(KIND=dp),INTENT(OUT),OPTIONAL   :: rmvolume,rmbond,rmangle,rmrot,&
      rmtrans,temp,pressure,BETA,pmswap,pmvolume,pmtraion,pmtrans,cutoff,&
      inner_cutoff

    CHARACTER(LEN=*),INTENT(OUT),OPTIONAL :: ensemble,program,restart_file,&
       molecules_file,moves_file,coords_file,&
       energy_file,displacement_file,cell_file,dat_file,data_file,box2_file

    CHARACTER(LEN=*), PARAMETER :: routine_name = "get_mc_par", &
      routineP = module_name//':'//routine_name

!   ---------------------------------------------------------------------------

    IF ( PRESENT ( nstep ) ) nstep = mc_par % nstep
    IF ( PRESENT ( iuptrans ) ) iuptrans = mc_par % iuptrans
    IF ( PRESENT ( iupvolume ) ) iupvolume = mc_par % iupvolume
    IF ( PRESENT ( nchain ) ) nchain = mc_par % nchain
    IF ( PRESENT ( nunits ) ) nunits = mc_par % nunits
    IF ( PRESENT ( nunits_tot ) ) nunits_tot = mc_par % nunits_tot
    IF ( PRESENT ( nmoves ) ) nmoves = mc_par % nmoves
    IF ( PRESENT ( nvolmoves ) ) nvolmoves = mc_par % nvolmoves
    IF ( PRESENT ( nswapmoves ) ) nswapmoves = mc_par % nswapmoves
    IF ( PRESENT ( rm ) ) rm = mc_par % rm
    IF ( PRESENT ( cl ) ) cl = mc_par % cl
    IF ( PRESENT ( diff ) ) diff = mc_par % diff
    IF ( PRESENT ( nstart ) ) nstart = mc_par % nstart
    IF ( PRESENT ( source ) ) source = mc_par % source
    IF ( PRESENT ( group ) ) group = mc_par % group

    IF ( PRESENT ( lclassical ) ) lclassical = mc_par % lclassical
    IF ( PRESENT ( lclassical_run ) ) lclassical_run = &
      mc_par % lclassical_run
    IF ( PRESENT ( ionode ) ) ionode = mc_par % ionode
    IF ( PRESENT ( lrestart ) ) lrestart = mc_par % lrestart
    IF ( PRESENT ( lstop ) ) lstop = mc_par % lstop

    IF ( PRESENT ( cutoff ) ) cutoff = mc_par % cutoff
    IF ( PRESENT ( inner_cutoff ) ) inner_cutoff = mc_par % inner_cutoff
    IF ( PRESENT ( rmvolume ) ) rmvolume = mc_par % rmvolume
    IF ( PRESENT ( rmbond ) ) rmbond = mc_par % rmbond
    IF ( PRESENT ( rmangle ) ) rmangle = mc_par % rmangle
    IF ( PRESENT ( rmrot ) ) rmrot = mc_par % rmrot
    IF ( PRESENT ( rmtrans ) ) rmtrans = mc_par % rmtrans
    IF ( PRESENT ( temp ) ) temp = mc_par % temp
    IF ( PRESENT ( pressure ) ) pressure = mc_par % pressure
    IF ( PRESENT ( BETA ) ) BETA = mc_par % BETA
    IF ( PRESENT ( pmswap ) ) pmswap = mc_par % pmswap
    IF ( PRESENT ( pmvolume) ) pmvolume = mc_par % pmvolume
    IF ( PRESENT ( pmtraion ) ) pmtraion = mc_par % pmtraion
    IF ( PRESENT ( pmtrans ) ) pmtrans = mc_par % pmtrans

    IF ( PRESENT ( ensemble ) ) ensemble = mc_par % ensemble
    IF ( PRESENT ( program ) ) program = mc_par % program
    IF ( PRESENT ( restart_file ) ) restart_file = mc_par % restart_file
    IF ( PRESENT ( moves_file ) ) moves_file = mc_par % moves_file
    IF ( PRESENT ( coords_file ) ) coords_file = mc_par % coords_file
    IF ( PRESENT ( molecules_file ) ) molecules_file = mc_par % molecules_file
    IF ( PRESENT ( energy_file ) ) energy_file = mc_par % energy_file
    IF ( PRESENT ( displacement_file ) ) displacement_file = mc_par % displacement_file
    IF ( PRESENT ( cell_file ) ) cell_file = mc_par % cell_file
    IF ( PRESENT ( dat_file ) ) dat_file = mc_par % dat_file
    IF ( PRESENT ( data_file ) ) data_file = mc_par % data_file
    IF ( PRESENT ( box2_file ) ) box2_file = mc_par % box2_file


  END SUBROUTINE get_mc_par

 SUBROUTINE set_mc_par ( mc_par, nstep, iuptrans, iupvolume, nchain,&
       nunits,nunits_tot,nmoves,nvolmoves,nswapmoves,rm,cl,diff,nstart,&
       source,group,lclassical,ionode,lrestart,lstop,rmvolume,rmbond,rmangle,&
       rmrot,rmtrans,temp,pressure,BETA,pmswap,pmvolume,pmtraion,pmtrans,&
       ensemble,program,restart_file,molecules_file,moves_file,coords_file,&
       energy_file,displacement_file,cell_file,dat_file,data_file,box2_file,&
       lclassical_run)

! *****************************************************************************
!!****s* mc_types.F/get_mc_par
!!
!!   NAME
!!     set_mc_par
!!
!!   FUNCTION
!!     changes the private elements of the mc_parameters_type
!!
!!   ARGUMENTS
!!     - mc_par: the structure mc parameters you want 
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************

    TYPE(mc_simulation_parameters_type), POINTER       :: mc_par
    INTEGER,INTENT(IN),OPTIONAL :: nstep,iuptrans,iupvolume,nchain,nunits,&
      nunits_tot,nmoves,nvolmoves,nswapmoves,rm,cl,diff,nstart,source,&
      group

    LOGICAL,INTENT(IN),OPTIONAL     :: lclassical,ionode,lrestart,lstop,&
      lclassical_run

    REAL(KIND=dp),INTENT(IN),OPTIONAL   :: rmvolume,rmbond,rmangle,rmrot,&
      rmtrans,temp,pressure,BETA,pmswap,pmvolume,pmtraion,pmtrans

    CHARACTER(LEN=*),INTENT(IN),OPTIONAL :: ensemble,program,restart_file,&
       molecules_file,moves_file,coords_file,&
       energy_file,displacement_file,cell_file,dat_file,data_file,box2_file

    CHARACTER(LEN=*), PARAMETER :: routine_name = "set_mc_par", &
      routineP = module_name//':'//routine_name

!   ---------------------------------------------------------------------------

! I've commented out all the ones that are read in in read_mc_section and
! shouldn't be changed in the course of a simulation

!    IF ( PRESENT ( nstep ) ) mc_par % nstep = nstep
!    IF ( PRESENT ( iuptrans ) ) mc_par % iuptrans = iuptrans
!    IF ( PRESENT ( iupvolume ) ) mc_par % iupvolume = iupvolume
    IF ( PRESENT ( nchain ) ) mc_par % nchain = nchain
    IF ( PRESENT ( nunits ) ) mc_par % nunits = nunits
    IF ( PRESENT ( nunits_tot ) ) mc_par % nunits_tot = nunits_tot
!    IF ( PRESENT ( nmoves ) ) mc_par % nmoves = nmoves
!    IF ( PRESENT ( nvolmoves ) ) mc_par % nvolmoves = nvolmoves
!    IF ( PRESENT ( nswapmoves ) ) mc_par % nswapmoves = nswapmoves
    IF ( PRESENT ( rm ) ) mc_par % rm = rm
    IF ( PRESENT ( cl ) ) mc_par % cl = cl
    IF ( PRESENT ( diff ) ) mc_par % diff = diff
    IF ( PRESENT ( nstart ) ) mc_par % nstart = nstart
!    IF ( PRESENT ( source ) ) mc_par % source = source
!    IF ( PRESENT ( group ) ) mc_par % group = group

!    IF ( PRESENT ( lclassical ) ) mc_par % lclassical = lclassical
!    IF ( PRESENT ( ionode ) ) mc_par % ionode = ionode
!    IF ( PRESENT ( lrestart ) ) mc_par % lrestart = lrestart
!    IF ( PRESENT ( lstop ) ) mc_par % lstop = lstop

    IF ( PRESENT ( rmvolume ) ) mc_par % rmvolume = rmvolume
    IF ( PRESENT ( rmbond ) ) mc_par % rmbond = rmbond
    IF ( PRESENT ( rmangle ) ) mc_par % rmangle = rmangle
    IF ( PRESENT ( rmrot ) ) mc_par % rmrot = rmrot
    IF ( PRESENT ( rmtrans ) ) mc_par % rmtrans = rmtrans
!    IF ( PRESENT ( temp ) ) mc_par % temp = temp
!    IF ( PRESENT ( pressure ) ) mc_par % pressure = pressure
!    IF ( PRESENT ( BETA ) ) mc_par % BETA = BETA
!    IF ( PRESENT ( pmswap ) ) mc_par % pmswap = pmswap
!    IF ( PRESENT ( pmvolume) ) mc_par % pmvolume = pmvolume
!    IF ( PRESENT ( pmtraion ) ) mc_par % pmtraion = pmtraion
!    IF ( PRESENT ( pmtrans ) ) mc_par % pmtrans = pmtrans

!    IF ( PRESENT ( ensemble ) ) mc_par % ensemble = ensemble
    IF ( PRESENT ( program ) ) mc_par % program = program
!    IF ( PRESENT ( restart_file ) ) mc_par % restart_file = restart_file
!    IF ( PRESENT ( moves_file ) ) mc_par % moves_file = moves_file
!    IF ( PRESENT ( coords_file ) ) mc_par % coords_file = coords_file
!    IF ( PRESENT ( molecules_file ) ) mc_par % molecules_file = molecules_file
!    IF ( PRESENT ( energy_file ) ) mc_par % energy_file = energy_file
!    IF ( PRESENT ( displacement_file ) ) mc_par % displacement_file = displacement_file
!    IF ( PRESENT ( cell_file ) ) mc_par % cell_file = cell_file
!    IF ( PRESENT ( dat_file ) ) mc_par % dat_file = dat_file
!    IF ( PRESENT ( data_file ) ) mc_par % data_file = data_file
!    IF ( PRESENT ( box2_file ) ) mc_par % box2_file = box2_file


  END SUBROUTINE set_mc_par

!******************************************************************************
!!****** mc/read_mc_section [1.0] *
!!
!!   NAME
!!     read_mc_section
!!
!!   FUNCTION
!!     reads in the Monte Carlo simulation parameters from an input file
!!
!!   ARGUMENTS
!!     - mc_par: the structure that will store the parameters
!!     - globenv: the global environment for the simulation
!!
!!   AUTHOR
!!     MJM
!!
!!   NOTES
!!I---------------------------------------------------------------------------I
!!I SECTION: &MC ... &END                                                     I
!!I                                                                           I
!!I program:     (traditional,gemc)       type of simulation
!!I nstep:                                number of steps
!!I lstop:       (yes,no)                 if no, peform nstep cycles
!!I lclassical:  (yes,no)                 pre-bias moves with classical 
!!                                        potential
!!I lclassical_run:  (yes,no)             using a classical potential instead
!!                                        of Quickstep
!!I nmoves:                               number of moves between Quickstep
!!                                        calculations
!!I rmbond:                               maximum displacement for bond
!!                                        changes (in angstroms)
!!I rmangle:                              maximum displacement for angle
!!                                        changes (in degrees)
!!I rmtrans:                              maximum displacement for molecule
!!                                        translations (in angstroms)
!!I rmrot:                                maximum displacement for molecule
!!                                        rotations (in degrees)
!!I rmvolume:                             maximum displacement for volume
!!                                        changes (in cubic angstroms)
!!I iuptrans:                             every iuptrans step, update maximum
!!                                        displacements for angle,bond,rot,
!!                                        trans
!!I iupvolume:                            every iupvolume step, update maximum
!!                                        volume displacement
!!I cutoff:                               the cutoff for the classical 
!!                                        potential
!!I inner_cutoff:                         a smaller cutoff for the classical 
!!                                        potential
!!I pmswap:                               probability of doing a swap move
!!I pmvolume:                             probability of doing a volume move
!!I pmtraion:                             probability of doing a conformational
!!                                        move
!!I pmtrans:                              probability of doing a translation
!!I temp:                                 temperature (in Kelvin)
!!I pressure:                             pressure (in bar)
!!I lrestart:     (yes,no)                read initial coordinates from a
!!                                        restart file
!!I restart_file:                         name of the restart file
!!I box2_file:                            name of the input file that contains
!!                                        the coords for the second box in
!!                                        a GEMC simulation
!!I---------------------------------------------------------------------------I
!!*****************************************************************************

SUBROUTINE read_mc_section ( mc_par, globenv )


    TYPE(mc_simulation_parameters_type), &
      INTENT(OUT)                            :: mc_par
    TYPE(global_environment_type), &
      INTENT(IN)                             :: globenv

    CHARACTER(LEN=*), PARAMETER      &
      :: routineN = "read_mc_section"

    CHARACTER(LEN=20)                        :: string
    CHARACTER(LEN=5)                         :: label
    INTEGER                                  :: handle, ia, ie, ierror, ilen, &
                                                iw, stop_num
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
  mc_par % nvolmoves = 1
  mc_par % nswapmoves = 1
  mc_par % lclassical_run = .FALSE.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  mc_par % cutoff = 7.8d0 / angstrom
  mc_par % inner_cutoff = 4.0d0 / angstrom

  iw = globenv % scr
  IF(globenv%ionode) WRITE ( iw, * )

!..filenames
  mc_par % program = globenv % input_file_name
  CALL xstring ( mc_par % program, ia, ie )
  mc_par%coords_file = mc_par % program(ia:ie) // '.coordinates'
  mc_par%molecules_file = mc_par % program(ia:ie) // '.molecules'
  mc_par%moves_file = mc_par % program(ia:ie) // '.moves'
  mc_par%energy_file = mc_par % program(ia:ie) // '.energy'
  mc_par%cell_file = mc_par % program(ia:ie) // '.cell'
  mc_par%displacement_file= mc_par % program(ia:ie)&
      // '.max_displacements'
  mc_par%data_file = mc_par % program(ia:ie) // '.data'
  stop_num=ie-3
  mc_par%dat_file = mc_par % program(ia:stop_num) // 'dat'

!..parse the input section
  label = '&MC'
  CALL parser_init(globenv % input_file_name,globenv)
  CALL search_label(label,ierror,ignore_case=.TRUE.)
  IF (ierror /= 0 ) THEN
     IF( mc_par % ionode ) THEN
        WRITE ( iw, *) ierror
        WRITE ( iw, '( a )' ) ' No input section &MC found '
     ENDIF
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
        CASE ( 'NSWAPMOVES' )
          CALL get_next ( mc_par % nswapmoves )
        CASE ( 'LSTOP' )
           CALL get_next ( string, ilen )
           CALL uppercase ( string )
           IF ( string(1:2)=='NO' ) mc_par % lstop = .FALSE.
           IF ( string(1:3)=='YES' ) mc_par % lstop = .TRUE.
        CASE ( 'LCLASSICAL_RUN' )
           CALL get_next ( string, ilen )
           CALL uppercase ( string )
           IF ( string(1:3)=='YES' ) mc_par % lclassical_run = .TRUE.
        CASE ( 'LCLASSICAL' )
           CALL get_next ( string, ilen )
           CALL uppercase ( string )
           IF ( string(1:2)=='NO' ) THEN
              mc_par % lclassical = .FALSE.
              mc_par % nvolmoves = 1
              mc_par % nswapmoves = 1
           ENDIF
           IF ( string(1:3)=='YES' ) mc_par % lclassical = .TRUE.
        CASE ( 'RESTART' )
           CALL get_next ( string, ilen )
           CALL uppercase ( string )
           IF ( string(1:2)=='NO' ) mc_par % lrestart = .FALSE.
           IF ( string(1:3)=='YES' ) mc_par % lrestart = .TRUE.
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
        CASE ( 'PMSWAP' )
           CALL get_next ( mc_par % pmswap )
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
!      IF ((mc_par % nmoves .EQ. 1) .AND. mc_par % lclassical ) THEN
!         STOP 'If lclassical=.TRUE., nmoves must be greater than 1.'
!      ENDIF

!..write some information to output
  IF (mc_par % ionode .AND. globenv % print_level>=0) THEN
     WRITE ( iw, '( A )' ) ' MC| Monte Carlo Protocol '
     WRITE ( iw, '( A,T71,I10 )' ) ' MC| total number of steps ', &
          mc_par % nstep
     WRITE ( iw, '( A,T71,F10.3 )' ) ' MC| pmswap ', &
          mc_par % pmswap
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
     WRITE ( iw, '( A,T44,A )' ) ' MC| Name of output molecules file:',&
                   TRIM(mc_par %molecules_file)
     WRITE ( iw, '( A,T44,A )' ) ' MC| Name of output moves file:',&
                   TRIM(mc_par % moves_file)
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
     IF ( mc_par % lrestart ) THEN
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
     IF ( mc_par % lclassical_run ) THEN
          WRITE ( iw, '(A)') ' MC| This will be a classical simulation.'
     ELSE
          WRITE ( iw, '(A)') ' MC| This will be a Quickstep simulation.'
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
     IF (mc_par%pmswap .GT. 0.0d0 ) THEN
          WRITE( iw, '(A,I4,A,A)') ' MC| ',mc_par%nswapmoves,&
          ' classical swap insertions will be attempted ',&
          'per Quickstep swap move'
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

END MODULE mc_types

!******************************************************************************
