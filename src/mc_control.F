!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** mc_control.F *********************************************************
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
  USE kinds,                           ONLY: dp,&
                                             dp_size
  USE mathconstants,                   ONLY: pi
  USE mc_misc,                         ONLY: mc_make_dat_file
  USE mc_types,                        ONLY: mc_simulation_parameters_type
  USE parser,                          ONLY: get_next,&
                                             p_error,&
                                             parser_end,&
                                             parser_init,&
                                             read_line,&
                                             search_label,&
                                             stop_parser,&
                                             test_next
  USE particle_list_types,             ONLY: particle_list_type
  USE physcon,                         ONLY: angstrom,&
                                             boltzmann,&
                                             joule
  USE qs_parser,                       ONLY: open_file
  USE simulation_cell,                 ONLY: cell_type
  USE string_utilities,                ONLY: uppercase,&
                                             xstring
  USE termination,                     ONLY: stop_memory
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

  PRIVATE 
! *** Global parameters ***

  CHARACTER(LEN=*), PARAMETER :: module_name = "mc_control"

  PUBLIC :: read_mc_section, &
            write_mc_restart,read_mc_restart


!****************************************************************************

CONTAINS

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
!!I pmswap:                               probability of doing a swap move
!!I pmvolume:                             probability of doing a volume move
!!I pmtraion:                             probability of doing a conformational
!!                                        move
!!I pmtrans:                              probability of doing a translation
!!I temp:                                 temperature (in Kelvin)
!!I pressure:                             pressure (in bar)
!!I restart:     (yes,no)                 read initial coordinates from a
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

    CHARACTER(LEN=*), PARAMETER              :: routineN = "read_mc_section"

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

  iw = globenv % scr
  IF(globenv%ionode) WRITE ( iw, * )

!..filenames
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

! *****************************************************************************
!!****s* mc_control.F/write_mc_restart
!!
!!   NAME
!!     write_mc_restart
!!
!!   FUNCTION
!!     writes the coordinates of the current step to a file that can
!!     be read in at the start of the next simulation
!!
!!   ARGUMENTS
!!     - mc_par: the mc parameters for the force env
!!     - force_env: the force environment to write the coords from
!!     - nnstep: how many steps the simulation has run
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************
SUBROUTINE write_mc_restart ( nnstep, mc_par, force_env )

    INTEGER, INTENT(IN)                      :: nnstep
    TYPE(mc_simulation_parameters_type), &
      POINTER                                :: mc_par
    TYPE(force_env_type), POINTER            :: force_env

    CHARACTER(LEN=*), PARAMETER              :: routineN = "write_mc_restart"

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
                 file_action='WRITE',file_status='REPLACE')
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
  DO ipart=1,mc_par%nunits_tot,mc_par%nunits
      WRITE(unit,'(1X,A,1X,3(F15.10,3X),A)') &
                 'O ',particles%els(ipart)%r(1:3)*angstrom,' H2O'
      WRITE(unit,'(1X,A,1X,3(F15.10,3X),A)') &
                 'H ',particles%els(ipart+1)%r(1:3)*angstrom,' H2O'
      WRITE(unit,'(1X,A,1X,3(F15.10,3X),A)') &
                 'H ',particles%els(ipart+2)%r(1:3)*angstrom,' H2O'
  ENDDO
  CLOSE(unit)

! end the timing
  CALL timestop(0.0d0,handle)


END SUBROUTINE write_mc_restart

! *****************************************************************************
!!****s* mc_control.F/read_mc_restart
!!
!!   NAME
!!     read_mc_restart
!!
!!   FUNCTION
!!     reads the input coordinates of the simulation from a file written above
!!
!!   ARGUMENTS
!!     - mc_par: the mc parameters for the force env
!!     - force_env: the force environment to write the coords from
!!     - iw: the unit to write an error message to, in case current
!!           simulation parameters don't match what's in the restart file
!!
!!   AUTHOR
!!     MJM
!!
!!*****************************************************************************
SUBROUTINE read_mc_restart ( mc_par, force_env ,iw )

    TYPE(mc_simulation_parameters_type), &
      POINTER                                :: mc_par
    TYPE(force_env_type), POINTER            :: force_env
    INTEGER, INTENT(IN)                      :: iw

    CHARACTER(LEN=*), PARAMETER              :: routineN = "read_mc_restart"

    CHARACTER(LEN=20)                        :: chardum, ensemble
    INTEGER                                  :: handle, ipart, istat, &
                                                nunits_tot, unit
    REAL(KIND=dp)                            :: box_length, temp
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: r
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_error_type)                      :: error
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: subsys
    TYPE(particle_list_type), POINTER        :: particles

!------------------------------------------------------------------------------
! begin the timing of the subroutine

      CALL timeset(routineN,'I','',handle)

! open the file and read some simulation parameters
      CALL open_file(file_name=mc_par%restart_file,unit_number=unit,&
         file_action='READ',file_status='OLD')

      READ(unit,*) mc_par%nstart
      READ(unit,*) temp,nunits_tot
      READ(unit,*) ensemble

! allocate the array holding the coordiantes
      ALLOCATE (r(3,nunits_tot),STAT=istat)
      IF (istat /= 0) CALL stop_memory(routineN,module_name,__LINE__,&
         "r",3*nunits_tot*dp_size)

! do some checking
      IF(ABS(temp - mc_par%temp) .GT. 0.01d0) THEN
         IF(mc_par%ionode) THEN
            WRITE(iw,*) 'The temperature in the restart file is ',&
               'not the same as the input file.'
            WRITE(iw,*) 'Input file temp =',mc_par%temp
            WRITE(iw,*) 'Restart file temp =',temp
         ENDIF
         STOP
      ENDIF
      IF(nunits_tot .NE. mc_par%nunits_tot) THEN
         IF(mc_par%ionode) THEN
            WRITE(iw,*) 'The total number of units in the restart ',&
            'file is not the same as the input file.'
            WRITE(iw,*) 'Input file units =',mc_par%nunits_tot
            WRITE(iw,*) 'Restart file units =',nunits_tot
         ENDIF
         mc_par%nunits_tot=nunits_tot
         mc_par%nchain=mc_par%nunits_tot/mc_par%nunits
!     STOP
      ENDIF
      IF(ensemble .NE. mc_par%ensemble) THEN
         IF(mc_par%ionode) THEN
            WRITE(iw,*) 'The ensemble in the restart file is ',&
               'not the same as the input file.'
            WRITE(iw,*) 'Input file ensemble =',mc_par%ensemble
            WRITE(iw,*) 'Restart file ensemble =',ensemble
         ENDIF
         STOP
      ENDIF

! get the cell length and coordinates
      CALL force_env_get(force_env,cell=cell,subsys=subsys,&
         error=error) 
      CALL cp_subsys_get(subsys(1)%subsys, &
         particles=particles, error=error)

      READ(unit,*)  box_length  ! in angstroms
      IF(ABS(box_length - cell%hmat(1,1)*angstrom) .GT. 0.0001d0 ) THEN
         IF (mc_par%ionode) THEN
            WRITE(iw,*) 'The cell length in the restart file is ',&
               'not the same as the input file.'
            WRITE(iw,*) 'Input file cell length =',cell%hmat(1,1)*angstrom
            WRITE(iw,*) 'Restart file cell length =',box_length
         ENDIF
!     STOP
      ENDIF
      
      READ(unit,*)
! read in the particle coordinates
      DO ipart=1,mc_par%nunits_tot
         READ(unit,*) chardum,r(1:3,ipart),chardum
         r(1:3,ipart)=r(1:3,ipart)/angstrom
      ENDDO
      CLOSE(unit)
      
! write the dat file so we can make a new force_env
      IF(mc_par%ionode) &
         CALL mc_make_dat_file(r(:,:),mc_par%nunits_tot,&
         box_length,mc_par%dat_file)

! end the timing
      CALL timestop(0.0d0,handle)

! deallcoate
      DEALLOCATE(r)

END SUBROUTINE read_mc_restart


END MODULE mc_control

!******************************************************************************
