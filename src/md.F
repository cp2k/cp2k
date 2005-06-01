!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****** cp2k/md [1.0] *
!!
!!   NAME
!!     md
!!
!!   FUNCTION
!!
!!   AUTHOR  INTEGER, PARAMETER, PUBLIC, DIMENSION(9) ::
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE md
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_logger_type,&
                                             cp_warning_level
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
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_get,&
                                             section_vals_val_get
  USE input_cp2k,                      ONLY: parsed_cp2k_input
  USE input_constants,                 ONLY: nve_ensemble,&
                                             nvt_ensemble,&
                                             npt_i_ensemble,&
                                             npt_f_ensemble,&
                                             nph_ensemble,&
                                             nph_uniaxial_ensemble,&
                                             nph_uniaxial_damped_ensemble,&
                                             isokin_ensemble, &
                                             reftraj_ensemble,&
                                             nhc_global,&
                                             nhc_massive,&
                                             ensemble_name_types,&
                                             ensemble_name_index,&
                                             md_init_start,&
                                             md_pos_restart,&
                                             md_pos_vel_restart,&
                                             md_all_restart,&
                                             restart_name_types,&
                                             restart_name_index,&
                                             dump_xyz,&
                                             dump_crd,&
                                             dump_atomic,&
                                             dump_xmol,&
                                             output_name_types,&
                                             output_name_index

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
     REAL (KIND=dp) :: f_annealing
     INTEGER        :: ensemble
     INTEGER        :: read_type
     INTEGER        :: output_format
     CHARACTER ( LEN = 20 ) :: PROGRAM
     LOGICAL :: constraint
     LOGICAL :: annealing
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
!!     Subroutine read_md_section(simpar, globenv, mdio)
!!       Implicit None
!!       Type(simulation_parameters_type), Intent (OUT):: simpar
!!       Type(global_environment_type), Intent (IN):: globenv
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

  SUBROUTINE read_md_section ( simpar, globenv, mdio, error)
    IMPLICIT NONE
    ! Arguments
    TYPE(simulation_parameters_type), &
         INTENT(OUT)                            :: simpar
    TYPE(global_environment_type), &
         INTENT(IN)                             :: globenv
    TYPE(mdio_parameters_type), INTENT(OUT)     :: mdio
    TYPE(cp_error_type), INTENT(inout)            :: error
    ! Local Variables
    CHARACTER(LEN=20)                        :: string
    CHARACTER(LEN=5)                         :: label
    INTEGER                                  :: ia, ie, ierror, ILEN, iw
    TYPE(section_vals_type), POINTER         :: md_section
    TYPE(unit_convert_type), POINTER         :: unit

    NULLIFY(unit, md_section)
    CALL set_md_defaults(simpar, globenv, mdio)
    iw = globenv % scr

    md_section => section_vals_get_subs_vals(globenv%input_file,&
       "MOTION%MD_NEW",error=error)
    IF (parsed_cp2k_input(md_section,check_this_section=.TRUE.,error=error)) THEN
       CALL read_md_new_section(simpar, globenv, mdio, md_section, error)
    ELSE
!
!TODEL_START_TL
!
!..parse the input section
       label = '&MD'
       CALL parser_init(globenv % input_file_name,globenv)
       CALL search_label(label,ierror,ignore_case=.TRUE.)
       IF (ierror /= 0 ) THEN
          IF( globenv % ionode ) &
               WRITE ( iw, '( a )' ) ' No input section &MD found '
       ELSE
          CALL read_line
          DO WHILE ( test_next() /= 'X' )
             ILEN = 6
             CALL get_next ( string, ILEN )
             CALL uppercase ( string )
             SELECT CASE ( string )
             CASE DEFAULT
                CALL p_error()
                CALL stop_parser ( 'read_md_section', 'unknown option' )
                
             CASE ( 'ANNEAL' )
                simpar % annealing = .TRUE.
                CALL get_next ( simpar % f_annealing )
                
             CASE ( 'OUTPUT' )
                ILEN = 0
                CALL get_next ( string, ILEN )
                CALL uppercase ( string )
                SELECT CASE (string ( 1:ILEN ))
                CASE ("XYZ")
                   simpar % output_format = dump_xyz
                CASE ("CRD")
                   simpar % output_format = dump_crd
                CASE ("ATOMIC")
                   simpar % output_format = dump_atomic
                CASE ("XMOL")
                   simpar % output_format = dump_xmol
                CASE DEFAULT
                   call p_error
                END SELECT
             CASE ( 'ENSEMB' )
                ILEN = 0
                CALL get_next ( string, ILEN )
                CALL uppercase ( string )
                SELECT CASE(string( 1:ILEN ))
                CASE ("NVE")
                   simpar % ensemble = nve_ensemble
                CASE ("NVT")
                   simpar % ensemble = nvt_ensemble
                CASE ("NPT_I")
                   simpar % ensemble = npt_i_ensemble
                 CASE ("NPT_F")
                   simpar % ensemble = npt_f_ensemble
                CASE ("NPH")
                   simpar % ensemble = nph_ensemble
                CASE ("NPH_UNIAXIAL")
                   simpar % ensemble = nph_uniaxial_ensemble            
                CASE ("NPH_UNIAXIAL_DUMPED")
                   simpar % ensemble = nph_uniaxial_damped_ensemble
                CASE ("ISOKIN")
                   simpar % ensemble = isokin_ensemble
                CASE ("REFTRAJ")
                   simpar % ensemble = reftraj_ensemble
                CASE DEFAULT 
                   CALL p_error()
                END SELECT
             CASE ( 'STEPS' )
                CALL get_next ( simpar % nsteps )
                
             CASE ( 'WALLTI' ) ! WALLTIME
                CALL get_next ( simpar % walltime )
             CASE ( 'TIMEST' )
                CALL get_next ( simpar % dt )
                IF ( test_next () == 'C' ) THEN
                   CALL get_next ( string, ILEN )
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
                ILEN = 0
                CALL get_next ( string, ILEN )
                CALL uppercase ( string )
                SELECT CASE(string ( 1:ILEN ))
                CASE ("INIT")
                   simpar % read_type = md_init_start
                CASE ("POS")
                   simpar % read_type = md_pos_restart
                CASE ("POSVEL")
                   simpar % read_type = md_pos_vel_restart
                CASE ("ALL")
                   simpar % read_type = md_all_restart
                CASE DEFAULT
                   call p_error()
                END SELECT
             CASE ( 'PV_AVA' )
                ILEN = 0
                CALL get_next ( string, ILEN )
                CALL uppercase ( string )
                IF ( string(1:1)=='T' ) simpar % pv_availability = .TRUE.
                IF ( string(1:1)=='F' ) simpar % pv_availability = .FALSE.
                
             CASE ( 'PV_NUM' )
                ILEN = 0
                CALL get_next ( string, ILEN )
                CALL uppercase ( string )
                IF ( string(1:1)=='T' ) simpar % pv_numerical=.TRUE.
                IF ( string(1:1)=='F' ) simpar % pv_numerical=.FALSE.
                
             CASE ( 'CONSTR' )
                ILEN = 0
                CALL get_next ( string, ILEN )
                CALL uppercase ( string )
                IF ( string(1:2)=='ON' ) simpar % constraint = .TRUE.
                IF ( string(1:3)=='OFF' ) simpar % constraint = .FALSE.
                
             CASE ( 'SHAKE' )
                CALL get_next ( simpar % shake_tol )
             CASE ( 'UNIAXI' )
                DO WHILE ( test_next() == 'C' )
                   
                   ILEN = 6
                   CALL get_next ( string, ILEN )
                   CALL uppercase ( string )
                   
                   SELECT CASE ( string )
                   CASE DEFAULT
                      CALL p_error()
                      CALL stop_parser ( 'read_md_section%unixial_shock', &
                           'unknown suboption' )
                   CASE ( 'PRESSU' )
                      CALL get_next ( simpar % p0 )
                      IF ( test_next () == 'C' ) THEN
                         CALL get_next ( simpar % units % press, ILEN )
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
                      simpar % ensemble = nph_uniaxial_damped_ensemble
                   CASE ( 'TIMECO' )
                      CALL get_next ( simpar % tau_cell )
                      IF ( test_next () == 'C' ) THEN
                         CALL get_next ( string, ILEN )
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
                   
                   ILEN = 6
                   CALL get_next ( string, ILEN )
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
                         CALL get_next ( string, ILEN )
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
                   
                   ILEN = 6
                   CALL get_next ( string, ILEN )
                   CALL uppercase ( string )
                   
                   SELECT CASE ( string )
                   CASE DEFAULT
                      CALL p_error()
                      CALL stop_parser ( 'read_md_section%barost', &
                           'unknown suboption' )
                   CASE ( 'TIMECO' )
                      CALL get_next ( simpar % tau_cell )
                      IF ( test_next () == 'C' ) THEN
                         CALL get_next ( string, ILEN )
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
                         CALL get_next ( simpar % units % press, ILEN )
                      ELSE
                         simpar % units % press = 'BAR'
                      ENDIF
                      CALL uppercase ( simpar % units % press )
                   END SELECT
                END DO
                
             CASE ( 'FILES' )
                DO WHILE ( test_next() == 'C' )
                   
                   ILEN = 6
                   CALL get_next ( string, ILEN )
                   CALL uppercase ( string )
                   
                   SELECT CASE ( string )
                   CASE DEFAULT
                      CALL p_error()
                      CALL stop_parser ( 'read_md_section%files', &
                           'unknown suboption' )
                   CASE ( 'FILEBO' )
                      ILEN = 0
                      CALL get_next ( string, ILEN )
                      mdio % crd_file_name = string(1:ILEN) // '.crd'
                      mdio % vel_file_name = string(1:ILEN) // '.vel'
                      mdio % ptens_file_name = string(1:ILEN) // '.ptens'
                      mdio % ener_file_name = string(1:ILEN) // '.ener'
                      mdio % pot_file_name = string(1:ILEN) // '.pot'
                      mdio % temp_file_name = string(1:ILEN) // '.temp'
                      mdio % dump_file_name = string(1:ILEN) // '.restart'
                      mdio % coef_file_name = string(1:ILEN) // '.coef'
                   CASE ( 'COORDI' )
                      ILEN = 0
                      CALL get_next ( mdio % crd_file_name, ILEN )
                   CASE ( 'ENERGI' )
                      ILEN = 0
                      CALL get_next ( mdio % ener_file_name, ILEN )
                   CASE ( 'POTENT' )
                      ILEN = 0
                      CALL get_next ( mdio % pot_file_name, ILEN )
                   CASE ( 'PRESSU' )
                      ILEN = 0
                      CALL get_next ( mdio % ptens_file_name, ILEN )
                   CASE ( 'VELOCI' )
                      ILEN = 0
                      CALL get_next ( mdio % vel_file_name, ILEN )
                   CASE ( 'TEMPER' )
                      ILEN = 0
                      CALL get_next ( mdio % temp_file_name, ILEN )
                   CASE ( 'DUMP' )
                      ILEN = 0
                      CALL get_next ( mdio % dump_file_name, ILEN )
                   CASE ( 'COEF' )
                      ILEN = 0
                      CALL get_next ( mdio % coef_file_name, ILEN )
                   END SELECT
                END DO
                
             CASE ( 'PRINT' )
                DO WHILE ( test_next() == 'C' )

                   ILEN = 6
                   CALL get_next ( string, ILEN )
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
                   
                   ILEN = 9
                   CALL get_next ( string, ILEN )
                   CALL uppercase ( string )
                   
                   SELECT CASE ( string )
                   CASE DEFAULT
                      CALL p_error()
                      CALL stop_parser ( 'read_md_section%qs_md_param', &
                           'unknown suboption' )
                   CASE ( 'NEIGHB' )
                      ILEN = 5
                      CALL get_next ( string, ILEN )
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
!
!TODEL_END_TL
!
    END IF
!
!TODEL_START_TL -- The conversion is intrinsic to the new input style.. not needed anymore
!
    IF (.NOT.parsed_cp2k_input(md_section,check_this_section=.TRUE.,error=error)) THEN
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
    END IF    
!
!TODEL_END_TL
!
!..write some information to output
    IF (globenv % ionode .AND. globenv % print_level>=0) THEN
       CALL get_cp2k_units ( unit )
       WRITE ( iw, '( A )' ) ' MD| Molecular Dynamics Protocol '
       WRITE ( iw, '( A,T61,A )' ) ' MD| Output Format ', &
            ADJUSTR ( output_name_types(simpar % output_format)(1:20))
       WRITE ( iw, '( A,T61,A )' ) ' MD| Ensemble Type ', & 
            ADJUSTR ( ensemble_name_types(simpar % ensemble)(1:20))
       WRITE ( iw, '( A,T61,A )' ) ' MD| Restart Type ', &
            ADJUSTR ( restart_name_types(simpar % read_type)(1:20))
       WRITE ( iw, '( A,T71,I10 )' ) ' MD| Number of Time Steps ', &
            simpar % nsteps
       WRITE ( iw, '( A,T71,F10.2 )' ) ' MD| Target Walltime [sec]', &
            simpar % walltime
       WRITE ( iw, '( A,A2,A,T71,F10.2 )' ) ' MD| Time Step [','fs','] ',&
            simpar % dt /unit%fs
       WRITE ( iw, '( A,T71,F10.2 )' ) ' MD| Temperature [K] ', &
            simpar % temp_ext/ unit % kelvin
       WRITE ( iw, '( A,T71,F10.2 )' ) ' MD| Temperature tolerance [K] ', &
            simpar % temp_tol/ unit % kelvin
       IF ( simpar % annealing ) &
            WRITE ( iw, '( A,T71,F10.2 )' ) ' MD| Annealing ion factor      ', &
            simpar % f_annealing
       IF ( simpar % ensemble == nvt_ensemble .OR. simpar % ensemble == npt_i_ensemble &
            .OR. simpar % ensemble == npt_f_ensemble ) THEN
          WRITE ( iw, '( A )' ) ' MD| Nose-Hoover-Chain parameters '
          WRITE ( iw, '( A,T77,I4 )' ) ' MD| Nose-Hoover-Chain length ', &
               simpar % nhclen
          WRITE ( iw, '( A,A4,A,T71,F10.4 )' ) &
               ' MD| Nose-Hoover-Chain time constant [', 'fs' ,'] ',&
               simpar % tau_nhc/unit%fs
          WRITE ( iw, '( A,T77,I4 )' ) ' MD| Order of Yoshida integrator ', &
               simpar % nyosh
          WRITE ( iw, '( A,T77,I4 )' ) ' MD| Number of multiple time steps', &
               simpar % nc
       END IF
       IF (  simpar % ensemble == npt_i_ensemble .OR. simpar % ensemble == npt_f_ensemble ) THEN
          WRITE ( iw, '( A,A3,A, T71, F10.2 )' ) &
               ' MD| Pressure [','Bar','] ', simpar % p_ext/unit%bar
          WRITE ( iw, '( A,A4,A, T71, F10.4 )' ) &
               ' MD| Barostat time constant [','fs','] ', &
               simpar % tau_cell/unit%fs
       END IF
       IF (  simpar % ensemble == isokin_ensemble ) THEN
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

END SUBROUTINE read_md_section

!******************************************************************************
!!****** md/read_md_new_section [1.0] *
!!
!!   NAME
!!     read_md_new_section
!!
!!
!!   FUNCTION
!!     Parses the MD input section 
!!
!!   AUTHOR
!!     teo
!!
!!   MODIFICATION HISTORY
!!  
!!
!!*** *************************************************************************  
SUBROUTINE read_md_new_section(simpar, globenv, mdio, md_section, error)
  IMPLICIT NONE
  ! Arguments
  TYPE(simulation_parameters_type), &
       INTENT(INOUT)                          :: simpar
  TYPE(global_environment_type), &
       INTENT(IN)                             :: globenv
  TYPE(mdio_parameters_type), INTENT(OUT)     :: mdio
  TYPE(cp_error_type), INTENT(inout)          :: error
  TYPE(section_vals_type), POINTER            :: md_section
  ! Local Variables
  CHARACTER(len=*), PARAMETER :: routineN = 'read_md_section', &
       routineP = moduleN//':'//routineN
  TYPE(section_vals_type), POINTER            :: tmp_section
  LOGICAL :: explicit, failure
  INTEGER :: n_rep

  failure = .FALSE.
  IF (.NOT.failure) THEN
     CALL section_vals_val_get(md_section,"ANNEALING",r_val=simpar%f_annealing,error=error)
     IF (simpar%f_annealing /= 1.0_dp) simpar % annealing = .TRUE.
     CALL section_vals_val_get(md_section,"ENSEMBLE",i_val=simpar%ensemble,error=error)
     CALL section_vals_val_get(md_section,"STEPS",i_val=simpar%nsteps,error=error)
     CALL section_vals_val_get(md_section,"TIMESTEP",r_val=simpar%dt,error=error)
     CALL section_vals_val_get(md_section,"TEMPERATURE",r_val=simpar%temp_ext,error=error)
     CALL section_vals_val_get(md_section,"TEMP_TOL",r_val=simpar%temp_tol,error=error)
     CALL section_vals_val_get(md_section,"RESTART",i_val=simpar%read_type,error=error)
     CALL section_vals_val_get(md_section,"PV_AVA",l_val=simpar%pv_availability,error=error)
     CALL section_vals_val_get(md_section,"PV_NUM",l_val=simpar%pv_numerical,error=error)

     tmp_section => section_vals_get_subs_vals(globenv%input_file,"GLOBAL",error=error)
     CALL section_vals_val_get(tmp_section,"OUTPUT",i_val=simpar%output_format,error=error)
     CALL section_vals_val_get(tmp_section,"WALLTIME",r_val=simpar%walltime,error=error)

     tmp_section => section_vals_get_subs_vals(globenv%input_file,"FORCE_EVAL%SUBSYS%CONSTRAINT",error=error)
     CALL section_vals_get(tmp_section,explicit=explicit,error=error)
     IF (explicit) THEN 
        simpar%constraint = .TRUE.
        CALL section_vals_val_get(tmp_section,"SHAKE",r_val=simpar%shake_tol,error=error)
     END IF

     IF (simpar%ensemble == nph_uniaxial_ensemble.OR.simpar%ensemble == nph_uniaxial_damped_ensemble) THEN
        tmp_section => section_vals_get_subs_vals(md_section,"UNIAXIAL",error=error)
        CALL section_vals_val_get(tmp_section,"PRESSURE",r_val=simpar%p0,error=error)
        CALL section_vals_val_get(tmp_section,"ENERGY",r_val=simpar%e0,error=error)
        CALL section_vals_val_get(tmp_section,"VOLUME",r_val=simpar%v0,error=error)
        IF (simpar%ensemble == nph_uniaxial_damped_ensemble) THEN
           CALL section_vals_val_get(tmp_section,"TEMP_B",r_val=simpar%temp_baro,n_rep_val=n_rep,error=error)
           IF (n_rep == 0) CPPostcondition(.FALSE.,cp_failure_level,routineP,error,failure)
        END IF
        CALL section_vals_val_get(tmp_section,"TIMECO",r_val=simpar%tau_cell,error=error)
        CALL section_vals_val_get(tmp_section,"VSHOC",r_val=simpar%v_shock,error=error)
     END IF

     tmp_section => section_vals_get_subs_vals(md_section,"BAROSTAT",error=error)
     CALL section_vals_get(tmp_section,explicit=explicit,error=error)
     IF (explicit) THEN
        CALL section_vals_val_get(tmp_section,"PRESSURE",r_val=simpar%p_ext,error=error)
        CALL section_vals_val_get(tmp_section,"TIMECO",r_val=simpar%tau_cell,error=error)
     END IF
     
     tmp_section => section_vals_get_subs_vals(md_section,"NOSE",error=error)
     CALL section_vals_get(tmp_section,explicit=explicit,error=error)
     IF (explicit) THEN
        CALL section_vals_val_get(tmp_section,"LENGTH",i_val=simpar%nhclen,error=error)
        CALL section_vals_val_get(tmp_section,"YOSHIDA",i_val=simpar%nyosh,error=error)
        CALL section_vals_val_get(tmp_section,"TIMECON",r_val=simpar%tau_nhc,error=error)
        CALL section_vals_val_get(tmp_section,"MTS",i_val=simpar%nc,error=error)
     END IF
     
     tmp_section => section_vals_get_subs_vals(md_section,"QS_MD",error=error)
     CALL section_vals_get(tmp_section,explicit=explicit,error=error)
     IF (explicit) THEN
        CALL section_vals_val_get(tmp_section,"NEIGHBOR_LISTS_FROM_SCRATCH",&
             l_val=simpar%qs_md_param%neighbor_lists_from_scratch, error=error)
     END IF
  END IF

END SUBROUTINE read_md_new_section

!******************************************************************************
!!****** md/set_md_defaults [1.0] *
!!
!!   NAME
!!     set_md_defaults
!!
!!
!!   FUNCTION
!!     Initializes the simpar and mdio environment
!!
!!   AUTHOR
!!     teo
!!
!!   MODIFICATION HISTORY
!!  
!!
!!*** *************************************************************************  
  SUBROUTINE set_md_defaults(simpar, globenv, mdio)
    IMPLICIT NONE
    ! Arguments
    TYPE(simulation_parameters_type), &
         INTENT(OUT)                            :: simpar
    TYPE(global_environment_type), &
         INTENT(IN)                             :: globenv
    TYPE(mdio_parameters_type), INTENT(OUT)     :: mdio    
    ! Local Variables
    INTEGER                                  :: ia, ie

    simpar % units % time = 'FS'
    simpar % units % press = 'BAR'
    simpar % units % tau = 'FS'
    simpar % units % tau_cell = 'FS'
    simpar % units % verlet_skin = 'ANGSTROMS'
    simpar % output_format = dump_xmol
    simpar % ensemble = nve_ensemble
    simpar % nsteps = 0
    simpar % walltime = -1.0_dp ! default: don't care
    simpar % dt = 0.0_dp
    simpar % temp_ext = 0.0_dp
    simpar % temp_baro = 0.0_dp
    simpar % temp_tol = 0.0_dp
    simpar % p_ext = 0.0_dp
    simpar % f_annealing = 1.0_dp
    simpar % read_type = md_pos_restart
    simpar % constraint = .FALSE.
    simpar % annealing = .FALSE.
    simpar % pv_numerical=.FALSE.
    simpar % pv_availability=.FALSE.
    simpar % shake_tol = 1.0E-6_dp
    simpar % nhclen = 3
    simpar % nc = 2
    simpar % nyosh = 3
    simpar % tau_nhc = 1000.0_dp
    simpar % tau_cell = 1000.0_dp
    simpar % qs_md_param % neighbor_lists_from_scratch = .TRUE.
    simpar % p0 = 0._dp
    simpar % v0 = 0._dp
    simpar % e0 = 0._dp
    simpar % v_shock = 0.0_dp
    
    !..filenames
    CALL xstring ( globenv % project_name, ia, ie )
    mdio % crd_file_name   = globenv % project_name(ia:ie) // '.crd'
    mdio % vel_file_name   = globenv % project_name(ia:ie) // '.vel'
    mdio % ptens_file_name = globenv % project_name(ia:ie) // '.ptens'
    mdio % ener_file_name  = globenv % project_name(ia:ie) // '.ener'
    mdio % pot_file_name   = globenv % project_name(ia:ie) // '.pot'
    mdio % temp_file_name  = globenv % project_name(ia:ie) // '.temp'
    mdio % dump_file_name  = globenv % project_name(ia:ie) // '.restart'
    mdio % coef_file_name  = globenv % project_name(ia:ie) // '.coef'

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
    
  END SUBROUTINE set_md_defaults





!******************************************************************************
!!****** md/initialize_velocities [1.0] *
!!
!!   NAME
!!     initialize_velocities
!!
!!   SYNOPSIS
!!     Subroutine initialize_velocities(simpar, part, globenv)
!!       Implicit None
!!       Type(simulation_parameters_type), Intent (IN):: simpar
!!       Type(particle_type), Intent (INOUT), Dimension(:):: part
!!       Type(global_environment_type), Intent (INOUT):: globenv
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

SUBROUTINE initialize_velocities ( simpar, part, globenv )


    TYPE(simulation_parameters_type), &
      INTENT(INOUT)                          :: simpar
    TYPE(particle_type), DIMENSION(:), &
      INTENT(INOUT)                          :: part
    TYPE(global_environment_type), &
      INTENT(INOUT)                          :: globenv

    INTEGER                                  :: i, iw, natoms
    REAL(KIND=dp)                                :: akin, denom, mass, sc, temp, mass_tot
    REAL(KIND=dp), DIMENSION(3)                  :: vcom, vcom_new
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(unit_convert_type), POINTER         :: unit

!------------------------------------------------------------------------------

  iw = globenv % scr
  IF (globenv % ionode .AND. globenv % print_level>1) THEN
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
  IF (globenv % ionode) THEN

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
          part ( i ) % v ( 1 ) = gasdev ( globenv % idum )
          part ( i ) % v ( 2 ) = gasdev ( globenv % idum )
          part ( i ) % v ( 3 ) = gasdev ( globenv % idum )
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
     temp = 0.0_dp
     sc   = 0.0_dp
     IF (simpar % nfree /= 0) THEN 
        temp = 2.0_dp * akin / REAL ( simpar % nfree,KIND=dp)
        ! scale velocities to get the correct initial temperature and
        sc = SQRT ( simpar % temp_ext / temp )
     END IF
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
  CALL mp_bcast ( mass_tot, globenv % source, globenv % group )
  DO i = 1, natoms
     CALL mp_bcast ( part ( i ) % v, globenv % source, globenv % group )
  END DO
#endif

  IF ( (simpar % read_type == md_init_start ) .AND. &
        globenv % ionode .AND. globenv % print_level > 1 ) THEN
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

  IF ( simpar % ensemble == nph_ensemble .OR.&
       simpar % ensemble == nph_uniaxial_ensemble .OR.&
       simpar % ensemble == nph_uniaxial_damped_ensemble ) &
  simpar % v_shock = simpar % v_shock * SQRT ( mass_tot )
END SUBROUTINE initialize_velocities



END MODULE md

!******************************************************************************
