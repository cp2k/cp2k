!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/md [1.0] *
!!
!!   NAME
!!     md
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE md
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE convert_units,                   ONLY: convert_to_cp2k_units
  USE cp2k_units,                      ONLY: get_cp2k_units,&
                                             unit_convert_type
  USE global_types,                    ONLY: global_environment_type
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_bcast
  USE parser,                          ONLY: get_next,&
                                             p_error,&
                                             parser_end,&
                                             parser_init,&
                                             read_line,&
                                             search_label,&
                                             stop_parser,&
                                             test_next
  USE particle_types,                  ONLY: particle_type
  USE string_utilities,                ONLY: uppercase,&
                                             xstring
  USE util,                            ONLY: gasdev

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: read_md_section, simulation_parameters_type, &
            initialize_velocities, mdio_parameters_type, &
            qs_md_parameters_type, simpar_units_type


!!****s* md/qs_md_parameters_type [1.0] *
!!
!!   NAME
!!     qs_md_parameters_type
!!
!!   FUNCTION
!!     parameters for the quickstep md
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - neighbor_lists_from_scratch: if the neighbor_lists should be build
!!       from scratch
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     03.2003 created [fawzi]
!!
!!   SOURCE
!***************************************************************************
  TYPE qs_md_parameters_type
     LOGICAL :: neighbor_lists_from_scratch
  END TYPE qs_md_parameters_type
!!***
!****************************************************************************

  TYPE simpar_units_type 
      CHARACTER ( LEN = 8 ) :: time, tau,            &
                               verlet_skin, tau_cell,&
                               press
  END TYPE 

  TYPE simulation_parameters_type
     INTEGER :: nsteps
     REAL (KIND=dp) :: walltime ! target maximum wall time of a simulation (alternative to nsteps)
     REAL (KIND=dp) :: dt
     REAL (KIND=dp) :: temp_ext
     REAL (KIND=dp) :: temp_baro
     REAL (KIND=dp) :: temp_tol
     REAL (KIND=dp) :: p_ext
     REAL (KIND=dp) :: v0
     REAL (KIND=dp) :: e0
     REAL (KIND=dp) :: v_shock
     REAL (KIND=dp) :: p0
     CHARACTER ( LEN = 20 ) :: ensemble
     CHARACTER ( LEN = 20 ) :: read_type
     CHARACTER ( LEN = 20 ) :: output_format
     CHARACTER ( LEN = 20 ) :: PROGRAM
     LOGICAL :: constraint
     LOGICAL :: pv_numerical
     LOGICAL :: pv_availability
     INTEGER :: nfree
     INTEGER :: nc
     INTEGER :: nyosh
     INTEGER :: nhclen
     REAL (KIND=dp) :: tau_nhc
     REAL (KIND=dp) :: tau_cell
     REAL (KIND=dp), POINTER, DIMENSION ( : ) :: dt_yosh
     REAL (KIND=dp) :: shake_tol
     TYPE ( qs_md_parameters_type ) :: qs_md_param
     TYPE ( simpar_units_type ) :: units
     LOGICAL :: fix_atoms
  END TYPE simulation_parameters_type

  TYPE mdio_parameters_type
     INTEGER :: icrd, icoef, ivel, iptens, iener, itemp, idump, iscreen, ipot
     CHARACTER ( LEN = 50 ) :: crd_file_name
     CHARACTER ( LEN = 50 ) :: vel_file_name
     CHARACTER ( LEN = 50 ) :: ptens_file_name
     CHARACTER ( LEN = 50 ) :: pot_file_name
     CHARACTER ( LEN = 50 ) :: ener_file_name
     CHARACTER ( LEN = 50 ) :: temp_file_name
     CHARACTER ( LEN = 50 ) :: dump_file_name
     CHARACTER ( LEN = 50 ) :: coef_file_name
  END TYPE mdio_parameters_type

! *** Global parameters ***
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='md'

!!***
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** md/read_md_section [1.0] *
!!
!!   NAME
!!     read_md_section
!!
!!   SYNOPSIS
!!     Subroutine read_md_section(simpar, mdpar, mdio)
!!       Implicit None
!!       Type(simulation_parameters_type), Intent (OUT):: simpar
!!       Type(global_environment_type), Intent (IN):: mdpar
!!       Type(mdio_parameters_type), Intent (OUT):: mdio
!!     End Subroutine read_md_section
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   NOTES
!!I---------------------------------------------------------------------------I
!!I SECTION: &md ... &end                                                     I
!!I                                                                           I
!!I output_format[xmol,atomic...]        format for the output of pos,vel,    I
!!I                                      forces, default is XMOL.             I
!!I ensemble:    [nve,nvt,isokin]        ensemble type                        I
!!I steps:       n                       total number of steps                I
!!I timestep:    dt                      time step [fs]                       I
!!I temperature: T                       temperature [K]                      I
!!I temperature: T_BARO                  temperature [K]                      I
!!I temp_tol: T                          temperature tolerance [K]            I
!!I restart:     [pos,posvel]            start type                           I
!!I constraints: [on,off]                constraints                          I
!!I pv_numerical: [true, false]          pv_numerical                         I
!!I pv_availability [true, false]        pv_availability                      I
!!I shake:       eps                     convergence for Shake algorithm      I
!!I verlet_skin  vs                      buffer size for Verlet lists         I
!!I cells        x y z                   Number of cells in each direction    I
!!I Uniaxial_shock        &              parameters for uniaxial shock        I
!!I   v_shock: v_shock  &                 velocity of shock front           I
!!I   pressure : p0       &                 Initial pressure                  I
!!I   energy   : e0       &                 Initial energy                    I
!!I   volume   : v0       &                 Initial compression               I
!!I Nose_parameter        &              parameter of Nose-Hoover-chains      I
!!I   length:    nhclen   &                 length of thermostat              I
!!I   Yoshida:   nyosh    &                 order of Yoshida integrator       I
!!I   mts:       nc       &                 number of multiple time steps     I
!!I   timecon:   tau_nhc                    NHC time constant                 I
!!I Barostat_parameter                      parameter of barostat             I
!!I   pressure: P         &                 pressure [bar]                    I
!!I   timecon:  tau_cell  &                 barostat time constant            I
!!I qs_md_parameter &                                                         I
!!I   neighbor_lists_from_scratch: [true, false]                              I
!!! fix_atoms                            fix atoms during an MD run
!!I FILES                           &                                         I
!!I   filebody    "filename"        &                                         I
!!I   coordinates "filename.crd"    &                                         I
!!I   coefficient "filename.coef"   &                                         I
!!I   velocities  "filename.vel"    &                                         I
!!I   pressure    "filename.ptens"  &                                         I
!!I   energies    "filename.ener"   &                                         I
!!I   temperature "filename.temp"   &                                         I
!!I   dump        "filename.restart"                                          I
!!I PRINT                           &                                         I
!!I   screen       iscreen          &                                         I
!!I   files        ifiles           &                                         I
!!I   coordinates  icrd             &                                         I
!!I   coefficents  icoef            &                                         I
!!I   velocities   ivel             &                                         I
!!I   pressure     iptens           &                                         I
!!I   energies     iener            &                                         I
!!I   temperature  itemp            &                                         I
!!I   dump         idump                                                      I
!!I                                                                           I
!!I                                                                           I
!!I---------------------------------------------------------------------------I
!!*** ************************************************************************

SUBROUTINE read_md_section ( simpar, mdpar, mdio )


    TYPE(simulation_parameters_type), &
      INTENT(OUT)                            :: simpar
    TYPE(global_environment_type), &
      INTENT(IN)                             :: mdpar
    TYPE(mdio_parameters_type), INTENT(OUT)  :: mdio

    CHARACTER(LEN=20)                        :: string
    CHARACTER(LEN=5)                         :: label
    INTEGER                                  :: ia, ie, ierror, ilen, iw

!------------------------------------------------------------------------------
!..defaults

  simpar % units % time = 'FS'
  simpar % units % press = 'BAR'
  simpar % units % tau = 'FS'
  simpar % units % tau_cell = 'FS'
  simpar % units % verlet_skin = 'ANGSTROMS'
  simpar % output_format = 'XMOL'
  simpar % ensemble = 'NVE'
  simpar % nsteps = 0
  simpar % walltime = -1.0_dp ! default: don't care
  simpar % dt = 0.0_dp
  simpar % temp_ext = 0.0_dp
  simpar % temp_baro = 0.0_dp
  simpar % temp_tol = 0.0_dp
  simpar % p_ext = 0.0_dp
  simpar % read_type = 'POS'
  simpar % constraint = .FALSE.
  simpar % pv_numerical=.FALSE.
  simpar % pv_availability=.FALSE.
  simpar % shake_tol = 1.0E-6_dp
  simpar % nhclen = 3
  simpar % nc = 2
  simpar % nyosh = 3
  simpar % tau_nhc = 1000.0_dp
  simpar % tau_cell = 1000.0_dp
  simpar % qs_md_param % neighbor_lists_from_scratch = .TRUE.
  simpar % fix_atoms = .FALSE.
  simpar % p0 = 0._dp
  simpar % v0 = 0._dp
  simpar % e0 = 0._dp
  simpar % v_shock = 0.0_dp

  iw = mdpar % scr

!..filenames
  CALL xstring ( mdpar % project_name, ia, ie )
  mdio % crd_file_name = mdpar % project_name(ia:ie) // '.crd'
  mdio % vel_file_name = mdpar % project_name(ia:ie) // '.vel'
  mdio % ptens_file_name = mdpar % project_name(ia:ie) // '.ptens'
  mdio % ener_file_name = mdpar % project_name(ia:ie) // '.ener'
  mdio % pot_file_name = mdpar % project_name(ia:ie) // '.pot'
  mdio % temp_file_name = mdpar % project_name(ia:ie) // '.temp'
  mdio % dump_file_name = mdpar % project_name(ia:ie) // '.restart'
  mdio % coef_file_name = mdpar % project_name(ia:ie) // '.coef'

!..print frequency
  mdio % icrd = 1
  mdio % icoef = 1
  mdio % ivel = 1
  mdio % iptens = 1
  mdio % iener = 1
  mdio % ipot = 1
  mdio % itemp = 1
  mdio % idump = 1
  mdio % iscreen = 1

!..parse the input section
  label = '&MD'
  CALL parser_init(mdpar % input_file_name,mdpar)
  CALL search_label(label,ierror,ignore_case=.TRUE.)
  IF (ierror /= 0 ) THEN
     IF( mdpar % ionode ) &
        WRITE ( iw, '( a )' ) ' No input section &MD found '
  ELSE
     CALL read_line
     DO WHILE ( test_next() /= 'X' )
        ilen = 6
        CALL get_next ( string, ilen )
        CALL uppercase ( string )

        SELECT CASE ( string )
        CASE DEFAULT
           CALL p_error()
           CALL stop_parser ( 'read_md_section', 'unknown option' )

        CASE ( 'FIX_AT' )
           simpar % fix_atoms = .TRUE.

        CASE ( 'OUTPUT' )
           ilen = 0
           CALL get_next ( string, ilen )
           CALL uppercase ( string )
           simpar % output_format = string ( 1:ilen )

        CASE ( 'ENSEMB' )
           ilen = 0
           CALL get_next ( string, ilen )
           CALL uppercase ( string )
           simpar % ensemble = string ( 1:ilen )

        CASE ( 'STEPS' )
           CALL get_next ( simpar % nsteps )

        CASE ( 'WALLTI' ) ! WALLTIME
           CALL get_next ( simpar % walltime )
        CASE ( 'TIMEST' )
           CALL get_next ( simpar % dt )
           IF ( test_next () == 'C' ) THEN
             CALL get_next ( string, ilen )
             CALL uppercase ( string )
             SELECT CASE ( string )
             CASE DEFAULT
               simpar % units % time = 'FS'  
             CASE ('FS','AU')
               simpar % units % time = string
             END SELECT
           ELSE
              simpar % units % time = 'FS'
           ENDIF
           CALL uppercase ( simpar % units % time )

        CASE ( 'TEMPER' )
           CALL get_next ( simpar % temp_ext )

        CASE ( 'TEMP_T' )
           CALL get_next ( simpar % temp_tol )

        CASE ( 'RESTAR' )
           ilen = 0
           CALL get_next ( string, ilen )
           CALL uppercase ( string )
           simpar % read_type = string ( 1:ilen )

        CASE ( 'PV_AVA' )
           ilen = 0
           CALL get_next ( string, ilen )
           CALL uppercase ( string )
           IF ( string(1:1)=='T' ) simpar % pv_availability = .TRUE.
           IF ( string(1:1)=='F' ) simpar % pv_availability = .FALSE.

        CASE ( 'PV_NUM' )
           ilen = 0
           CALL get_next ( string, ilen )
           CALL uppercase ( string )
           IF ( string(1:1)=='T' ) simpar % pv_numerical=.TRUE.
           IF ( string(1:1)=='F' ) simpar % pv_numerical=.FALSE.

        CASE ( 'CONSTR' )
           ilen = 0
           CALL get_next ( string, ilen )
           CALL uppercase ( string )
           IF ( string(1:2)=='ON' ) simpar % constraint = .TRUE.
           IF ( string(1:3)=='OFF' ) simpar % constraint = .FALSE.

        CASE ( 'SHAKE' )
           CALL get_next ( simpar % shake_tol )
        CASE ( 'UNIAXI' )
           DO WHILE ( test_next() == 'C' )

              ilen = 6
              CALL get_next ( string, ilen )
              CALL uppercase ( string )

              SELECT CASE ( string )
              CASE DEFAULT
                 CALL p_error()
                 CALL stop_parser ( 'read_md_section%unixial_shock', &
                      'unknown suboption' )
              CASE ( 'PRESSU' )
                 CALL get_next ( simpar % p0 )
                 IF ( test_next () == 'C' ) THEN
                   CALL get_next ( simpar % units % press, ilen )
                 ELSE
                    simpar % units % press = 'BAR'
                 ENDIF
                 CALL uppercase ( simpar % units % press )
              CASE ( 'ENERGY' )
                 CALL get_next ( simpar % e0 )
              CASE ( 'VOLUME' )
                 CALL get_next ( simpar % v0 )
              CASE ( 'TEMP_B' )
                 CALL get_next ( simpar % temp_baro )
                 simpar % ensemble = "NPH_UNIAXIAL_DAMPED"
              CASE ( 'TIMECO' )
                 CALL get_next ( simpar % tau_cell )
                 IF ( test_next () == 'C' ) THEN
                   CALL get_next ( string, ilen )
                   CALL uppercase ( string )
                   SELECT CASE ( string )
                   CASE DEFAULT
                    simpar % units % tau_cell = '1/CM'
                   CASE ('FS','AU')
                     simpar % units % tau_cell = string
                   END SELECT
                 ELSE
                    simpar % units % tau_cell = '1/CM'
                 ENDIF
                 CALL uppercase ( simpar % units % tau_cell )
              CASE ( 'V_SHOC' )
                 CALL get_next ( simpar % v_shock )
!!! input in m/sec
                 simpar % v_shock = simpar % v_shock * 1.e10_dp
                 simpar % v_shock = simpar % v_shock * 1.e-15_dp
              END SELECT
           END DO
        CASE ( 'NOSE_P' )
           DO WHILE ( test_next() == 'C' )

              ilen = 6
              CALL get_next ( string, ilen )
              CALL uppercase ( string )

              SELECT CASE ( string )
              CASE DEFAULT
                 CALL p_error()
                 CALL stop_parser ( 'read_md_section%nose', &
                      'unknown suboption' )
              CASE ( 'LENGTH' )
                 CALL get_next ( simpar % nhclen )
              CASE ( 'YOSHID' )
                 CALL get_next ( simpar % nyosh )
              CASE ( 'TIMECO' )
                 CALL get_next ( simpar % tau_nhc )
                 IF ( test_next () == 'C' ) THEN
                   CALL get_next ( string, ilen )
                   CALL uppercase ( string )
                   SELECT CASE ( string )
                   CASE DEFAULT
                    simpar % units % tau = '1/CM'
                   CASE ('FS','AU')
                     simpar % units % tau = string
                   END SELECT
                 ELSE
                    simpar % units % tau = '1/CM'
                 ENDIF
                 CALL uppercase ( simpar % units % tau )
              CASE ( 'MTS' )
                 CALL get_next ( simpar % nc )
              END SELECT
           END DO

        CASE ( 'BAROST' )
           DO WHILE ( test_next() == 'C' )

              ilen = 6
              CALL get_next ( string, ilen )
              CALL uppercase ( string )

              SELECT CASE ( string )
              CASE DEFAULT
                 CALL p_error()
                 CALL stop_parser ( 'read_md_section%barost', &
                      'unknown suboption' )
              CASE ( 'TIMECO' )
                 CALL get_next ( simpar % tau_cell )
                 IF ( test_next () == 'C' ) THEN
                   CALL get_next ( string, ilen )
                   CALL uppercase ( string )
                   SELECT CASE ( string )
                   CASE DEFAULT
                    simpar % units % tau_cell = '1/CM'
                   CASE ('FS','AU')
                     simpar % units % tau_cell = string
                   END SELECT
                 ELSE
                    simpar % units % tau_cell = '1/CM'
                 ENDIF
                 CALL uppercase ( simpar % units % tau_cell )
              CASE ( 'PRESSU' )
                 CALL get_next ( simpar % p_ext )
                 IF ( test_next () == 'C' ) THEN
                   CALL get_next ( simpar % units % press, ilen )
                 ELSE
                    simpar % units % press = 'BAR'
                 ENDIF
                 CALL uppercase ( simpar % units % press )
              END SELECT
           END DO

        CASE ( 'FILES' )
           DO WHILE ( test_next() == 'C' )

              ilen = 6
              CALL get_next ( string, ilen )
              CALL uppercase ( string )

              SELECT CASE ( string )
              CASE DEFAULT
                 CALL p_error()
                 CALL stop_parser ( 'read_md_section%files', &
                      'unknown suboption' )
              CASE ( 'FILEBO' )
                 ilen = 0
                 CALL get_next ( string, ilen )
                 mdio % crd_file_name = string(1:ilen) // '.crd'
                 mdio % vel_file_name = string(1:ilen) // '.vel'
                 mdio % ptens_file_name = string(1:ilen) // '.ptens'
                 mdio % ener_file_name = string(1:ilen) // '.ener'
                 mdio % pot_file_name = string(1:ilen) // '.pot'
                 mdio % temp_file_name = string(1:ilen) // '.temp'
                 mdio % dump_file_name = string(1:ilen) // '.restart'
                 mdio % coef_file_name = string(1:ilen) // '.coef'
              CASE ( 'COORDI' )
                 ilen = 0
                 CALL get_next ( mdio % crd_file_name, ilen )
              CASE ( 'ENERGI' )
                 ilen = 0
                 CALL get_next ( mdio % ener_file_name, ilen )
              CASE ( 'POTENT' )
                 ilen = 0
                 CALL get_next ( mdio % pot_file_name, ilen )
              CASE ( 'PRESSU' )
                 ilen = 0
                 CALL get_next ( mdio % ptens_file_name, ilen )
              CASE ( 'VELOCI' )
                 ilen = 0
                 CALL get_next ( mdio % vel_file_name, ilen )
              CASE ( 'TEMPER' )
                 ilen = 0
                 CALL get_next ( mdio % temp_file_name, ilen )
              CASE ( 'DUMP' )
                 ilen = 0
                 CALL get_next ( mdio % dump_file_name, ilen )
              CASE ( 'COEF' )
                 ilen = 0
                 CALL get_next ( mdio % coef_file_name, ilen )
              END SELECT
           END DO

        CASE ( 'PRINT' )
           DO WHILE ( test_next() == 'C' )

              ilen = 6
              CALL get_next ( string, ilen )
              CALL uppercase ( string )

              SELECT CASE ( string )
              CASE DEFAULT
                 CALL p_error()
                 CALL stop_parser ( 'read_md_section%print', &
                      'unknown suboption' )

              CASE ( 'SCREEN' )
                 CALL get_next ( mdio % iscreen )

              CASE ( 'FILES' )
                 CALL get_next ( ia )
                 mdio % icrd = ia
                 mdio % ivel = ia
                 mdio % iener = ia
                 mdio % ipot = ia
                 mdio % iptens = ia
                 mdio % itemp = ia
                 mdio % idump = ia
                 mdio % icoef = ia
              CASE ( 'COORDI' )
                 CALL get_next ( mdio % icrd )
              CASE ( 'ENERGI' )
                 CALL get_next ( mdio % iener )
              CASE ( 'POTENT' )
                 CALL get_next ( mdio % ipot )
              CASE ( 'PRESSU' )
                 CALL get_next ( mdio % iptens )
              CASE ( 'VELOCI' )
                 CALL get_next ( mdio % ivel )
              CASE ( 'TEMPER' )
                 CALL get_next ( mdio % itemp )
              CASE ( 'DUMP' )
                 CALL get_next ( mdio % idump )
              CASE ( 'COEF' )
                 CALL get_next ( mdio % icoef )
              END SELECT
           END DO

        CASE ( 'QS_MD_' )
           DO WHILE ( test_next() == 'C' )

              ilen = 9
              CALL get_next ( string, ilen )
              CALL uppercase ( string )

              SELECT CASE ( string )
              CASE DEFAULT
                 CALL p_error()
                 CALL stop_parser ( 'read_md_section%qs_md_param', &
                      'unknown suboption' )
              CASE ( 'NEIGHB' )
                 ilen = 5
                 CALL get_next ( string, ilen )
                 CALL uppercase ( string )
                 
                 SELECT CASE ( string )
                 CASE DEFAULT
                    CALL p_error()
                    CALL stop_parser ( 'read_md_section%qs_md_param', &
                         'neighbor_lists_from_scratch must be either true '//&
                         'or false' )
                 CASE ( 'TRUE')
                    simpar % qs_md_param % neighbor_lists_from_scratch = &
                         .TRUE.
                 CASE ( 'FALSE')
                    simpar % qs_md_param % neighbor_lists_from_scratch = &
                         .FALSE.
                 END SELECT
              END SELECT
           END DO
           
        END SELECT
        CALL read_line
     END DO

  END IF
  CALL parser_end
!..end of parsing the input section

!..write some information to output
  IF (mdpar % ionode .AND. mdpar % print_level>=0) THEN
     WRITE ( iw, '( A )' ) ' MD| Molecular Dynamics Protocol '
     WRITE ( iw, '( A,T61,A )' ) ' MD| Output Format ', &
          ADJUSTR ( simpar % output_format )
     WRITE ( iw, '( A,T61,A )' ) ' MD| Ensemble Type ', &
          ADJUSTR ( simpar % ensemble )
     WRITE ( iw, '( A,T61,A )' ) ' MD| Restart Type ', &
          ADJUSTR ( simpar % read_type )
     WRITE ( iw, '( A,T71,I10 )' ) ' MD| Number of Time Steps ', &
          simpar % nsteps
     WRITE ( iw, '( A,T71,F10.2 )' ) ' MD| Target Walltime [sec]', &
          simpar % walltime
     WRITE ( iw, '( A,A2,A,T71,F10.2 )' ) ' MD| Time Step [',simpar % units % time,'] ', simpar % dt
     WRITE ( iw, '( A,T71,F10.2 )' ) ' MD| Temperature [K] ', &
          simpar % temp_ext
     WRITE ( iw, '( A,T71,F10.2 )' ) ' MD| Temperature tolerance [K] ', &
          simpar % temp_tol
     IF ( ( simpar % ensemble ( 1:3 ) == 'NVT' ) &
          .OR. (simpar % ensemble ( 1:3 ) == 'NPT' ) ) THEN
        WRITE ( iw, '( A )' ) ' MD| Nose-Hoover-Chain parameters '
        WRITE ( iw, '( A,T77,I4 )' ) ' MD| Nose-Hoover-Chain length ', &
             simpar % nhclen
        WRITE ( iw, '( A,A4,A,T71,F10.4 )' ) &
             ' MD| Nose-Hoover-Chain time constant [',simpar % units % tau,'] ', simpar % tau_nhc
        WRITE ( iw, '( A,T77,I4 )' ) ' MD| Order of Yoshida integrator ', &
             simpar % nyosh
        WRITE ( iw, '( A,T77,I4 )' ) ' MD| Number of multiple time steps', &
             simpar % nc
     END IF
     IF ( simpar % ensemble ( 1:3 ) == 'NPT' ) THEN
        WRITE ( iw, '( A,A3,A, T71, F10.2 )' ) &
             ' MD| Pressure [',simpar % units % press,'] ', simpar % p_ext
        WRITE ( iw, '( A,A4,A, T71, F10.4 )' ) &
             ' MD| Barostat time constant [',simpar % units % tau_cell,'] ', simpar % tau_cell
     END IF
     IF ( simpar % ensemble ( 1:6 ) == 'ISOKIN' ) THEN
        WRITE ( iw, '( A )' ) ' MD| Simulation in the isokinetic ensemble'
     END IF
     IF ( simpar % constraint ) THEN
        WRITE ( iw, '( A )' ) ' MD| Constraints activated '
        WRITE ( iw, '( A,T71,G10.4 )' ) ' MD| Tolerance for shake ', &
             simpar % shake_tol
     END IF

     WRITE ( iw, '( A,T63,i10,A )' ) ' MD| Print MD information every', &
          mdio % iscreen, ' step(s)'
     WRITE ( iw, '( A,T20,A,T71,A )' ) ' MD| File type', &
          'Print frequency[steps]', 'File names'
     WRITE ( iw, '( A,T25,i5,T31,A )' ) ' MD| Coordinates', mdio % icrd, &
          ADJUSTR ( mdio % crd_file_name ( 1:50 ) )
     WRITE ( iw, '( A,T25,i5,T31,A )' ) ' MD| Velocities', mdio % ivel, &
          ADJUSTR ( mdio % vel_file_name ( 1:50 ) )
     WRITE ( iw, '( A,T25,i5,T31,A )' ) ' MD| Energies', mdio % iener, &
          ADJUSTR ( mdio % ener_file_name ( 1:50 ) )
     WRITE ( iw, '( A,T25,i5,T31,A )' ) ' MD| Potential Energies',     &
          mdio % ipot, ADJUSTR ( mdio % pot_file_name ( 1:50 ) )
     WRITE ( iw, '( A,T25,i5,T31,A )' ) ' MD| Temperature', mdio % itemp, &
          ADJUSTR ( mdio % temp_file_name ( 1:50 ) )
     WRITE ( iw, '( A,T25,i5,T31,A )' ) ' MD| Pressure Tensor', mdio % iptens, &
          ADJUSTR ( mdio % ptens_file_name ( 1:50 ) )
     WRITE ( iw, '( A,T25,i5,T31,A )' ) ' MD| Dump', mdio % idump, &
          ADJUSTR ( mdio % dump_file_name ( 1:50 ) )
     WRITE (iw,*)
  END IF

  CALL convert_to_cp2k_units ( simpar % units % time, time = simpar % dt ) 
  CALL convert_to_cp2k_units ( temperature = simpar % temp_ext ) 
  CALL convert_to_cp2k_units ( temperature = simpar % temp_baro ) 
  CALL convert_to_cp2k_units ( temperature = simpar % temp_tol ) 
  CALL convert_to_cp2k_units ( simpar % units % tau, tau_nhc = simpar % tau_nhc ) 
  CALL convert_to_cp2k_units ( simpar % units % tau_cell, tau_nhc = simpar % tau_cell ) 
  CALL convert_to_cp2k_units ( simpar % units % press, pressure = simpar % p_ext ) 
  CALL convert_to_cp2k_units ( simpar % units % press, pressure = simpar % p0 ) 
  CALL convert_to_cp2k_units ( "ANGSTROM", length = simpar % v_shock ) 
  CALL convert_to_cp2k_units ( "FS", time = simpar % v_shock, t_power = -1 ) 

END SUBROUTINE read_md_section

!******************************************************************************
!!****** md/initialize_velocities [1.0] *
!!
!!   NAME
!!     initialize_velocities
!!
!!   SYNOPSIS
!!     Subroutine initialize_velocities(simpar, part, mdpar)
!!       Implicit None
!!       Type(simulation_parameters_type), Intent (IN):: simpar
!!       Type(particle_type), Intent (INOUT), Dimension(:):: part
!!       Type(global_environment_type), Intent (INOUT):: mdpar
!!     End Subroutine initialize_velocities
!!
!!   FUNCTION
!!     Initializes the velocities to the Maxwellian distribution
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!*** *************************************************************************

SUBROUTINE initialize_velocities ( simpar, part, mdpar )


    TYPE(simulation_parameters_type), &
      INTENT(INOUT)                          :: simpar
    TYPE(particle_type), DIMENSION(:), &
      INTENT(INOUT)                          :: part
    TYPE(global_environment_type), &
      INTENT(INOUT)                          :: mdpar

    INTEGER                                  :: i, iw, natoms
    REAL(KIND=dp)                                :: akin, denom, mass, sc, temp, mass_tot
    REAL(KIND=dp), DIMENSION(3)                  :: vcom, vcom_new
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(unit_convert_type), POINTER         :: unit

!------------------------------------------------------------------------------

  iw = mdpar % scr
  IF (mdpar % ionode .AND. mdpar % print_level>1) THEN
     WRITE ( iw, '( A, A, A )' ) ' **********************', &
          ' begin of velocity initialization ', '***********************'
  END IF

! Initializing parameters
  vcom = 0.0_dp
  mass_tot = 0.0_dp
  vcom_new = 0.0_dp
  denom = 0.0_dp
  natoms = SIZE ( part )

  CALL get_cp2k_units ( unit ) 

! initializing velocities
  IF (mdpar % ionode) THEN

! computing the COM velocities
     DO i = 1, natoms
        atomic_kind => part ( i ) % atomic_kind
        CALL get_atomic_kind ( atomic_kind = atomic_kind, mass = mass )
        mass_tot = mass_tot + mass
        IF (part(i)%is_fixed) THEN
          part ( i ) % v ( 1 ) = 0.0_dp
          part ( i ) % v ( 2 ) = 0.0_dp
          part ( i ) % v ( 3 ) = 0.0_dp
        ELSE
          part ( i ) % v ( 1 ) = gasdev ( mdpar % idum )
          part ( i ) % v ( 2 ) = gasdev ( mdpar % idum )
          part ( i ) % v ( 3 ) = gasdev ( mdpar % idum )
          vcom ( 1 ) = vcom ( 1 ) + mass * part ( i ) % v ( 1 )
          vcom ( 2 ) = vcom ( 2 ) + mass * part ( i ) % v ( 2 )
          vcom ( 3 ) = vcom ( 3 ) + mass * part ( i ) % v ( 3 )
          denom = denom + mass
        END IF
     END DO
     vcom = vcom / denom

     akin = 0.0_dp
     DO i = 1, natoms
        atomic_kind => part ( i ) % atomic_kind
        CALL get_atomic_kind ( atomic_kind = atomic_kind, mass = mass )
! first subtract the COM
        IF (.NOT.part(i)%is_fixed) THEN
          part ( i ) % v ( : ) = part ( i ) % v ( : ) - vcom ( : )
! now compute the temp
          akin = akin + 0.5_dp *  &
           ( mass * part ( i ) % v ( 1 ) * part ( i ) % v ( 1 ) )
          akin = akin + 0.5_dp *  &
           ( mass * part ( i ) % v ( 2 ) * part ( i ) % v ( 2 ) )
          akin = akin + 0.5_dp *  &
           ( mass * part ( i ) % v ( 3 ) * part ( i ) % v ( 3 ) )
        END IF
     END DO
     temp = 2.0_dp * akin / REAL ( simpar % nfree,KIND=dp)
 
! scale velocities to get the correct initial temperature and
     sc = SQRT ( simpar % temp_ext / temp )
     DO i = 1, natoms
        atomic_kind => part ( i ) % atomic_kind
        CALL get_atomic_kind ( atomic_kind = atomic_kind, mass = mass )
        part ( i ) % v ( 1 ) = sc * part ( i ) % v ( 1 )
        part ( i ) % v ( 2 ) = sc * part ( i ) % v ( 2 )
        part ( i ) % v ( 3 ) = sc * part ( i ) % v ( 3 )
     END DO

! Moving the velocities so that the COM velocity is zero
     akin = 0.0_dp
     vcom_new = 0.0_dp
     DO i = 1, natoms
        atomic_kind => part ( i ) % atomic_kind
        CALL get_atomic_kind ( atomic_kind = atomic_kind, mass = mass )
        vcom_new ( : ) = vcom_new ( : ) + mass * part ( i ) % v ( : )
        akin = akin + 0.5_dp *  &
         ( mass * part ( i ) % v ( 1 ) * part ( i ) % v ( 1 ) )
        akin = akin + 0.5_dp *  &
         ( mass * part ( i ) % v ( 2 ) * part ( i ) % v ( 2 ) )
        akin = akin + 0.5_dp *  &
         ( mass * part ( i ) % v ( 3 ) * part ( i ) % v ( 3 ) )
     END DO
  END IF


#if defined(__parallel)
  CALL mp_bcast ( mass_tot, mdpar % source, mdpar % group )
  DO i = 1, natoms
     CALL mp_bcast ( part ( i ) % v, mdpar % source, mdpar % group )
  END DO
#endif

  IF ( ( simpar % read_type == 'INIT' ) .AND. &
        mdpar % ionode .AND. mdpar % print_level > 1 ) THEN
     temp = 2.0_dp * akin / REAL ( simpar % nfree,KIND=dp)
     WRITE ( iw, '( A, T61, F18.2, A2 )' ) &
          ' Initial Temperature ', temp / unit % kelvin, " K"
     WRITE ( iw, '( A, T61, F20.12 )' ) &
          ' Centre of mass velocity in direction x:', vcom_new ( 1 )
     WRITE ( iw, '( A, T61, F20.12 )' ) &
          ' Centre of mass velocity in direction y:', vcom_new ( 2 )
     WRITE ( iw, '( A, T61, F20.12 )' ) &
          ' Centre of mass velocity in direction z:', vcom_new ( 3 )
     WRITE ( iw, '( A, A, A, / )' ) ' ***********************', &
          ' end of velocity initialization ', '************************'
  END IF

  IF ( simpar % ensemble (1:3) == "NPH" ) &
  simpar % v_shock = simpar % v_shock * SQRT ( mass_tot )
END SUBROUTINE initialize_velocities



END MODULE md

!******************************************************************************
