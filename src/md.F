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
  
  USE global_types, ONLY : global_environment_type
  USE kinds, ONLY : dbl
  USE message_passing, ONLY : mp_bcast
  USE parser, ONLY : parser_init, parser_end, read_line, test_next, &
       cfield, p_error, get_real, get_int, stop_parser
  USE particle_types, ONLY : particle_type
  USE string_utilities, ONLY : uppercase, xstring
  USE util, ONLY : gasdev
  
  IMPLICIT NONE
  
  PRIVATE
  PUBLIC :: read_md_section, simulation_parameters_type, &  
       initialize_velocities, thermodynamic_type, mdio_parameters_type
  
  TYPE simulation_parameters_type
     INTEGER :: nsteps
     REAL ( dbl ) :: dt
     REAL ( dbl ) :: temp_ext
     REAL ( dbl ) :: p_ext
     CHARACTER ( LEN = 20 ) :: ensemble
     CHARACTER ( LEN = 20 ) :: read_type
     CHARACTER ( LEN = 20 ) :: program
     LOGICAL :: constraint
     INTEGER :: nfree
     INTEGER :: nc
     INTEGER :: nyosh
     INTEGER :: nhclen
     INTEGER, DIMENSION ( 3 ) :: n_cell
     REAL ( dbl ) :: tau_nhc
     REAL ( dbl ) :: tau_cell
     REAL ( dbl ), POINTER, DIMENSION ( : ) :: dt_yosh
     REAL ( dbl ) :: shake_tol
     REAL ( dbl ) :: verlet_skin
  END TYPE simulation_parameters_type
  
  TYPE thermodynamic_type
     REAL ( dbl ), DIMENSION ( 3, 3 ) :: ptens, pv_kin, pv, pv_const
     REAL ( dbl ) :: kin, pot, gspace, e_self, e_neut, e_bonded
  END TYPE thermodynamic_type
  
  TYPE mdio_parameters_type
     INTEGER :: icrd, ivel, iptens, iener, itemp, idump, iscreen
     CHARACTER ( LEN = 50 ) :: crd_file_name
     CHARACTER ( LEN = 50 ) :: vel_file_name
     CHARACTER ( LEN = 50 ) :: ptens_file_name
     CHARACTER ( LEN = 50 ) :: ener_file_name
     CHARACTER ( LEN = 50 ) :: temp_file_name
     CHARACTER ( LEN = 50 ) :: dump_file_name
  END TYPE mdio_parameters_type
  
!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** md/read_md_section [1.0] *
!!
!!   NAME
!!     read_md_section
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
!!I ensemble:    [nve,nvt]               ensemble type                        I
!!I steps:       n                       total number of steps                I
!!I timestep:    dt                      time step [fs]                       I
!!I temperature: T                       temperature [K]                      I
!!I restart:     [pos,posvel]            start type                           I
!!I constraints: [on,off]                constraints                          I
!!I shake:       eps                     convergence for Shake algorithm      I
!!I verlet_skin  vs                      buffer size for Verlet lists         I
!!I cells        x y z                   Number of cells in each direction    I
!!I Nose_parameter        &              parameter of Nose-Hoover-chains      I
!!I   length:    nhclen   &                 length of thermostat              I
!!I   Yoshida:   nyosh    &                 order of Yoshida integrator       I
!!I   mts:       nc       &                 number of multiple time steps     I
!!I   timecon:   tau_nhc                    NHC time constant                 I
!!I Barostat_parameter                      parameter of barostat             I
!!I   pressure: P         &                 pressure [bar]                    I
!!I   timecon:  tau_cell  &                 barostat time constant            I
!!I FILES                           &                                         I
!!I   filebody    "filename"        &                                         I
!!I   coordinates "filename.crd"    &                                         I
!!I   velocities  "filename.vel"    &                                         I
!!I   pressure    "filename.ptens"  &                                         I
!!I   energies    "filename.ener"   &                                         I
!!I   temperature "filename.temp"   &                                         I
!!I   dump        "filename.dum"                                              I
!!I PRINT                           &                                         I
!!I   screen       iscreen          &                                         I
!!I   files        ifiles           &                                         I
!!I   coordinates  icrd             &                                         I
!!I   velocities   ivel             &                                         I
!!I   pressure     iptens           &                                         I
!!I   energies     iener            &                                         I
!!I   temperature  itemp            &                                         I
!!I   dump         idump                                                      I
!!I                                                                           I
!!I                                                                           I
!!I---------------------------------------------------------------------------I
!!   SOURCE
!******************************************************************************

SUBROUTINE read_md_section ( simpar, mdpar, mdio )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( simulation_parameters_type ), INTENT ( OUT ) :: simpar
  TYPE ( global_environment_type ), INTENT ( IN ) :: mdpar
  TYPE ( mdio_parameters_type ), INTENT ( OUT ) :: mdio
  
! Locals
  INTEGER :: ierror, ilen, msglen, ia, ie, iw, allgrp, source
  CHARACTER ( LEN = 20 ) :: string
  CHARACTER ( LEN = 5 ) :: label
  
!------------------------------------------------------------------------------
  
!..defaults
  simpar % ensemble = 'NVE'
  simpar % nsteps = 0
  simpar % dt = 0.0_dbl
  simpar % temp_ext = 0.0_dbl
  simpar % p_ext = 0.0_dbl
  simpar % read_type = 'POS'
  simpar % constraint = .FALSE.
  simpar % shake_tol = 1.0E-6_dbl
  simpar % nhclen = 1
  simpar % nc = 1
  simpar % n_cell = 1
  simpar % nyosh = 1
  simpar % tau_nhc = 1000.0_dbl
  simpar % tau_cell = 1000.0_dbl
  simpar % verlet_skin = 1.0_dbl
  
  iw = mdpar % scr
     
!..filenames
  CALL xstring ( mdpar % project_name, ia, ie )
  mdio % crd_file_name = mdpar % project_name(ia:ie) // '.crd'
  mdio % vel_file_name = mdpar % project_name(ia:ie) // '.vel'
  mdio % ptens_file_name = mdpar % project_name(ia:ie) // '.ptens'
  mdio % ener_file_name = mdpar % project_name(ia:ie) // '.ener'
  mdio % temp_file_name = mdpar % project_name(ia:ie) // '.temp'
  mdio % dump_file_name = mdpar % project_name(ia:ie) // '.dump'
     
!..print frequency
  mdio % icrd = 1
  mdio % ivel = 1
  mdio % iptens = 1
  mdio % iener = 1
  mdio % itemp = 1
  mdio % idump = 1
  mdio % iscreen = 1
     
!..parse the input section
  label = '&MD'
  CALL parser_init(mdpar % input_file_name,label,ierror,mdpar)
  IF (ierror /= 0 ) THEN
     IF( mdpar % ionode ) &
        WRITE ( iw, '( a )' ) ' No input section &MD found '
  ELSE
     CALL read_line
     DO WHILE ( test_next() /= 'X' )
        ilen = 6
        CALL cfield ( string, ilen )
        CALL uppercase ( string )
        
        SELECT CASE ( string )
        CASE DEFAULT
           CALL p_error()
           CALL stop_parser ( 'read_md_section', 'unknown option' )
           
        CASE ( 'ENSEMB' )
           ilen = 0
           CALL cfield ( string, ilen )
           CALL uppercase ( string )
           simpar % ensemble = string ( 1:ilen )
           
        CASE ( 'STEPS' )
           simpar % nsteps = get_int()
           
        CASE ( 'VERLET' )
           simpar % verlet_skin = get_real()
           
        CASE ( 'CELL' )
           simpar % n_cell(1) = get_int()
           simpar % n_cell(2) = get_int()
           simpar % n_cell(3) = get_int()
           
        CASE ( 'TIMEST' )
           simpar % dt = get_real()
           
        CASE ( 'TEMPER' )
           simpar % temp_ext = get_real()
           
        CASE ( 'RESTAR' )
           ilen = 0
           CALL cfield ( string, ilen )
           CALL uppercase ( string )
           simpar % read_type = string ( 1:ilen )
           
        CASE ( 'CONSTR' )
           ilen = 0
           CALL cfield ( string, ilen )
           CALL uppercase ( string )
           IF ( string(1:2)=='ON' ) simpar % constraint = .TRUE.
           IF ( string(1:3)=='OFF' ) simpar % constraint = .FALSE.
           
        CASE ( 'SHAKE' )
           simpar % shake_tol = get_real()
           
        CASE ( 'NOSE_P' )
           DO WHILE ( test_next() == 'C' )
              
              ilen = 6
              CALL cfield ( string, ilen )
              CALL uppercase ( string )
              
              SELECT CASE ( string )
              CASE DEFAULT
                 CALL p_error()
                 CALL stop_parser ( 'read_md_section%nose', &
                      'unknown suboption' )
              CASE ( 'LENGTH' )
                 simpar % nhclen = get_int()
              CASE ( 'YOSHID' )
                 simpar % nyosh = get_int()
              CASE ( 'TIMECO' )
                 simpar % tau_nhc = get_real()
              CASE ( 'MTS' )
                 simpar % nc = get_int()
              END SELECT
           END DO
           
        CASE ( 'BAROST' )
           DO WHILE ( test_next() == 'C' )
              
              ilen = 6
              CALL cfield ( string, ilen )
              CALL uppercase ( string )
              
              SELECT CASE ( string )
              CASE DEFAULT
                 CALL p_error()
                 CALL stop_parser ( 'read_md_section%barost', &
                      'unknown suboption' )
              CASE ( 'TIMECO' )
                 simpar % tau_cell = get_real()
              CASE ( 'PRESSU' )
                 simpar % p_ext = get_real()
              END SELECT
           END DO
           
        CASE ( 'FILES' )
           DO WHILE ( test_next() == 'C' )
              
              ilen = 6
              CALL cfield ( string, ilen )
              CALL uppercase ( string )
              
              SELECT CASE ( string )
              CASE DEFAULT
                 CALL p_error()
                 CALL stop_parser ( 'read_md_section%files', &
                      'unknown suboption' )
              CASE ( 'FILEBO' )
                 ilen = 0
                 CALL cfield ( string, ilen )
                 mdio % crd_file_name = string(1:ilen) // '.crd'
                 mdio % vel_file_name = string(1:ilen) // '.vel'
                 mdio % ptens_file_name = string(1:ilen) // '.ptens'
                 mdio % ener_file_name = string(1:ilen) // '.ener'
                 mdio % temp_file_name = string(1:ilen) // '.temp'
                 mdio % dump_file_name = string(1:ilen) // '.dump'
              CASE ( 'COORDI' )
                 ilen = 0
                 CALL cfield ( mdio % crd_file_name, ilen )
              CASE ( 'ENERGI' )
                 ilen = 0
                 CALL cfield ( mdio % ener_file_name, ilen )
              CASE ( 'PRESSU' )
                 ilen = 0
                 CALL cfield ( mdio % ptens_file_name, ilen )
              CASE ( 'VELOCI' )
                 ilen = 0
                 CALL cfield ( mdio % vel_file_name, ilen )
              CASE ( 'TEMPER' )
                 ilen = 0
                 CALL cfield ( mdio % temp_file_name, ilen )
              CASE ( 'DUMP' )
                 ilen = 0
                 CALL cfield ( mdio % dump_file_name, ilen )
              END SELECT
           END DO
           
        CASE ( 'PRINT' )
           DO WHILE ( test_next() == 'C' )
              
              ilen = 6
              CALL cfield ( string, ilen )
              CALL uppercase ( string )
              
              SELECT CASE ( string )
              CASE DEFAULT
                 CALL p_error()
                 CALL stop_parser ( 'read_md_section%print', &
                      'unknown suboption' )
                 
              CASE ( 'SCREEN' )
                 mdio % iscreen = get_int()
                 
              CASE ( 'FILES' )
                 ia = get_int()
                 mdio % icrd = ia
                 mdio % ivel = ia
                 mdio % iener = ia
                 mdio % iptens = ia
                 mdio % itemp = ia
                 mdio % idump = ia
              CASE ( 'COORDI' )
                 mdio % icrd = get_int()
              CASE ( 'ENERGI' )
                 mdio % iener = get_int()
              CASE ( 'PRESSU' )
                 mdio % iptens = get_int()
              CASE ( 'VELOCI' )
                 mdio % ivel = get_int()
              CASE ( 'TEMPER' )
                 mdio % itemp = get_int()
              CASE ( 'DUMP' )
                 mdio % idump = get_int()
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
     WRITE ( iw, '( A,T61,A )' ) ' MD| Ensemble Type ', &
          ADJUSTR ( simpar % ensemble )
     WRITE ( iw, '( A,T61,A )' ) ' MD| Restart Type ', &
          ADJUSTR ( simpar % read_type )
     WRITE ( iw, '( A,T71,I10 )' ) ' MD| Number of Time Steps ', &
          simpar % nsteps
     WRITE ( iw, '( A,T71,F10.2 )' ) ' MD| Time Step [fs] ', simpar % dt
     WRITE ( iw, '( A,T71,F10.2 )' ) ' MD| Temperature [K] ', &
          simpar % temp_ext
     WRITE ( iw, '( A,T71,F10.2 )' ) &
          ' MD| Verlet buffer skin [Angstrom] ', simpar % verlet_skin
     IF ( ( simpar % ensemble ( 1:3 ) == 'NVT' ) &
          .OR. (simpar % ensemble ( 1:3 ) == 'NPT' ) ) THEN
        WRITE ( iw, '( A )' ) ' MD| Nose-Hoover-Chain parameters '
        WRITE ( iw, '( A,T77,I4 )' ) ' MD| Nose-Hoover-Chain length ', &
             simpar % nhclen
        WRITE ( iw, '( A,T71,F10.4 )' ) &
             ' MD| Nose-Hoover-Chain time constant [fs] ', simpar % tau_nhc
        WRITE ( iw, '( A,T77,I4 )' ) ' MD| Order of Yoshida integrator ', &
             simpar % nyosh
        WRITE ( iw, '( A,T77,I4 )' ) ' MD| Number of multiple time steps', &
             simpar % nc
     END IF
     IF ( simpar % ensemble ( 1:3 ) == 'NPT' ) THEN
        WRITE ( iw, '( A, T71, F10.2 )' ) &
             ' MD| Pressure [bar] ', simpar % p_ext
        WRITE ( iw, '( A, T71, F10.4 )' ) &
             ' MD| Barostat time constant [fs] ', simpar % tau_cell
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
     WRITE ( iw, '( A,T25,i5,T31,A )' ) ' MD| Temperature', mdio % itemp, &
          ADJUSTR ( mdio % temp_file_name ( 1:50 ) )
     WRITE ( iw, '( A,T25,i5,T31,A )' ) ' MD| Pressure Tensor', mdio % iptens, &
          ADJUSTR ( mdio % ptens_file_name ( 1:50 ) )
     WRITE ( iw, '( A,T25,i5,T31,A )' ) ' MD| Dump', mdio % idump, &
          ADJUSTR ( mdio % dump_file_name ( 1:50 ) )
     WRITE (iw,*)
  END IF
  
END SUBROUTINE read_md_section

!!*****
!******************************************************************************
!!****** md/initialize_velocities [1.0] *
!!
!!   NAME
!!     initialize_velocities
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
!!   SOURCE
!******************************************************************************

SUBROUTINE initialize_velocities ( simpar, part, mdpar )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( simulation_parameters_type ), INTENT ( IN ) :: simpar
  TYPE ( particle_type ), INTENT ( INOUT ), DIMENSION ( : ) :: part
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: mdpar
  
! Locals
  INTEGER :: i, id, natoms, iw
  REAL ( dbl ), DIMENSION ( 3 ) :: vcom
  REAL ( dbl ), DIMENSION ( 3 ) :: sum
  REAL ( dbl ) :: denom, sc, mass
  REAL ( dbl ) :: akin, temp
  LOGICAL :: domdec
  
!------------------------------------------------------------------------------
  
  iw = mdpar % scr
  IF (mdpar % ionode .AND. mdpar % print_level>1) THEN
     WRITE ( iw, '( A, A, A )' ) ' **********************', &
          ' begin of velocity initialization ', '***********************'
  END IF
  
! Initializing parameters
  sum = 0.0_dbl
  vcom = 0.0_dbl
  denom = 0.0_dbl
  natoms = SIZE ( part )
  
! initializing velocities
  IF (mdpar % ionode) THEN
     DO i = 1, natoms
        part(i) % v(1) = gasdev(mdpar % idum)
        part(i) % v(2) = gasdev(mdpar % idum)
        part(i) % v(3) = gasdev(mdpar % idum)
     END DO
     akin = 0.0_dbl
     DO i = 1, natoms
        mass = part(i) % prop % mass
        akin = akin + 0.5_dbl*(mass*part(i) % v(1)*part(i) % v(1))
        akin = akin + 0.5_dbl*(mass*part(i) % v(2)*part(i) % v(2))
        akin = akin + 0.5_dbl*(mass*part(i) % v(3)*part(i) % v(3))
     END DO
     
     temp = 2.0_dbl * akin / REAL ( simpar % nfree, dbl )
     
! scale velocities to get the correct initial temperature
     sc = SQRT ( simpar % temp_ext / temp )
     DO i = 1, natoms
        part(i) % v(1) = sc*part(i) % v(1)
        part(i) % v(2) = sc*part(i) % v(2)
        part(i) % v(3) = sc*part(i) % v(3)
     END DO
     
! Computing the COM velocities
     DO i = 1, natoms
        mass = part(i) % prop % mass
        vcom(1) = vcom(1) + mass*part(i) % v(1)
        vcom(2) = vcom(2) + mass*part(i) % v(2)
        vcom(3) = vcom(3) + mass*part(i) % v(3)
        denom = denom + mass
     END DO
     
! Moving the velocities so that the COM velocity is zero
     vcom = vcom / denom
     DO i = 1, natoms
        mass = part(i) % prop % mass
        part(i) % v ( : ) = part(i) % v ( : ) - vcom ( : )
        sum ( : ) = sum ( : ) + mass*part(i) % v ( : )
     END DO
  END IF
  
#if defined(__parallel)
  DO i = 1, natoms
     CALL mp_bcast(part(i) % v,mdpar % source,mdpar % group)
  END DO
#endif
  
  IF ( mdpar % ionode .AND. mdpar % print_level > 1 ) THEN
     WRITE ( iw, '( A, T69, E12.6 )' ) &
          ' Centre of mass velocity in direction x:', sum ( 1 )
     WRITE ( iw, '( A, T69, E12.6 )' ) &
          ' Centre of mass velocity in direction y:', sum ( 2 )
     WRITE ( iw, '( A, T69, E12.6 )' ) &
          ' Centre of mass velocity in direction z:', sum ( 3 )
     WRITE ( iw, '( A, A, A, / )' ) ' ***********************', &
          ' end of velocity initialization ', '************************'
  END IF
  
END SUBROUTINE initialize_velocities

!!*****
!******************************************************************************

END MODULE md

!******************************************************************************
