!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/cp2k_units [1.0] *
!!
!!   NAME
!!     cp2k_units
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

MODULE cp2k_units
  USE ewald_parameters_types, ONLY : ewald_parameters_type
  USE kinds, ONLY : dbl
  USE mathconstants, ONLY : pi
  USE md, ONLY : simulation_parameters_type
  USE molecule_types, ONLY : intra_parameters_type
  USE pair_potential, ONLY : potentialparm_type
  USE particle_types, ONLY : particle_type
  USE physcon, ONLY : e_charge, angstrom, a_bohr, e_mass, joule, a_mass, &
       boltzmann, massunit, bohr, evolt, h_bar, permittivity
  USE simulation_cell, ONLY : cell_type
  USE termination, ONLY : stop_program

  PRIVATE
  PUBLIC :: set_cp2k_units, unit_convert_type, get_cp2k_units

  TYPE unit_convert_type
     PRIVATE
     REAL ( dbl ) :: mass, fs, kcal, kelvin, bar, wavenum, angstrom
     CHARACTER ( LEN = 10 ) :: e_label, l_label, vol_label,  &
                               temp_label, pres_label,       &
                               angl_label
  END TYPE unit_convert_type
  
  TYPE ( unit_convert_type ), POINTER :: units

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** units/init_cp2k_units [1.0] *
!!
!!   NAME
!!     init_cp2k_units
!!
!!   FUNCTION
!!     routine which defines internal units and conversion factors
!!     internal units of CP2K are in atomic_units
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE init_cp2k_units ( )

  IMPLICIT NONE

!------------------------------------------------------------------------------

  ALLOCATE ( units, stat = isos )
  IF ( isos < 0 ) CALL stop_program ( )

  units % temp_label = '[K]'
  units % pres_label = '[bar]'
  units % angl_label = '[deg]'
  units % l_label = '[bohr]'
  units % vol_label = '[bohr^3]'
  units % e_label = '[hartree]'

! amu to mass AU
  units % mass = massunit
! fs to time AU
  units % fs = joule / 1.e15_dbl /  h_bar
! 1/cm to time AU
  units % wavenum = 1._dbl / 100._dbl / c_light * h_bar / joule
! kcal/mole to energy AU
  units % kcal = 1._dbl / kcalmol
! kelvin to energy AU
  units % kelvin = boltzmann / joule
! angstrom to length AU
  units % angstrom = bohr
! bar to pressure AU
  units % bar = joule / ( a_bohr ** 3 * 1.0E+05_dbl )
! permittivity to AU
  units % eps0 = joule * a_bohr / e_charge ** 2

END SUBROUTINE init_cp2k_units

!******************************************************************************
!!****** units/get_cp2k_units [1.0] *
!!
!!   NAME
!!     get_cp2k_units
!!
!!   FUNCTION
!!     routine which defines internal units and conversion factors
!!     internal units of CP2K are in atomic_units
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE get_cp2k_units ( unit_struc )
  IMPLICIT NONE

  TYPE ( unit_convert_type ), POINTER :: unit_struc
  unit_struc => units

END SUBROUTINE get_cp2k_units

!******************************************************************************
!!****** units/time_convert [1.0] *
!!
!!   NAME
!!     time_convert
!!
!!   FUNCTION
!!     Converts time in AU to the internal units of CP2K
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE time_convert ( time )
  IMPLICIT NONE
  REAL ( dbl ), INTENT ( INOUT ) :: time
  time = time / units % fs * 1.e-15_dbl / ( a_bohr * SQRT ( e_mass / joule ) )
END SUBROUTINE time_convert 

!******************************************************************************
!!****** cp2k_units/convert_to_cp2k_units [1.0] *
!!
!!   NAME
!!     convert_to_cp2k_units
!!
!!   FUNCTION
!!     Converts user defined input into AU
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE convert_to_cp2k_units ( unit_char, simpar, part, box, potparm, &
                                   intra_param, ewald_param )

  IMPLICIT NONE

! Arguments
  TYPE (unit_convert_type ), INTENT ( IN ), OPTIONAL :: unit_char
  TYPE (simulation_parameters_type ), INTENT ( INOUT ), OPTIONAL :: simpar
  TYPE (particle_type ), INTENT ( INOUT ), OPTIONAL :: part ( : )
  TYPE (cell_type ), INTENT ( INOUT ), OPTIONAL :: box
  TYPE (intra_parameters_type ), INTENT ( INOUT ), OPTIONAL :: intra_param
  TYPE (ewald_parameters_type ), INTENT ( INOUT ), OPTIONAL :: ewald_param
  TYPE (potentialparm_type ), INTENT ( INOUT ), OPTIONAL :: potparm ( :, : )

! Locals
  REAL ( dbl ) :: conv2k
  INTEGER :: i

!------------------------------------------------------------------------------

! convert hmat
  IF ( PRESENT ( box ) ) THEN
   IF ( box % unit_of_length_name == "ANGSTROMS" ) &
   box % hmat = box % hmat * units % angstrom
  ENDIF

  IF ( PRESENT ( simpar ) ) THEN
! convert time
     IF ( simpar % unit == 'FS' ) &
     simpar%dt = simpar % dt * units % fs

     CALL time_convert ( simpar % dt )
! convert temperature
     simpar%temp_ext = simpar % temp_ext * units % kelvin
! convert temperature of coefficents
     simpar%temp_coef = simpar % temp_coef * units % kelvin
! convert pressure
     simpar % p_ext = simpar % p_ext / units % bar
! convert verlet skin
     IF ( simpar % unit % verlet_skin == 'ANGSTROM' ) &
     simpar % verlet_skin = simpar % verlet_skin * units % angstrom

!!!!
! convert omega_p (wp) of PIMD
!     simpar%pimd_params%wp = simpar%pimd_params%wp * units%econv / units % tconv
!!!!
     
! converting extended system time constants
     SELECT CASE (simpar%ensemble)
     CASE ( 'NVT')
        IF ( simpar % units % tau == "1/CM" ) &
        simpar % tau_nhc = simpar % tau_nhc * units % wavenum
        IF ( simpar % units % tau == "FS" ) &
        simpar % tau_nhc = simpar % tau_nhc * units % fs

        CALL time_convert ( simpar % tau_nhc )
     CASE ( 'NPT_I')
        IF ( simpar % units % tau_cell == "1/CM" ) &
        simpar % tau_cell = simpar % tau_cell * units % wavenum
        IF ( simpar % units % tau_cell == "FS" ) &
        simpar % tau_cell = simpar % tau_cell * units % fs

        IF ( simpar % units % tau == "1/CM" ) &
        simpar % tau_nhc = simpar % tau_nhc * units % wavenum
        IF ( simpar % units % tau == "FS" ) &
        simpar % tau_nhc = simpar % tau_nhc * units % fs

        CALL time_convert ( simpar % tau_cell )
        CALL time_convert ( simpar % tau )
     CASE ( 'NPT_F')
        IF ( simpar % units % tau_cell == "1/CM" ) &
        simpar % tau_cell = simpar % tau_cell * units % wavenum
        IF ( simpar % units % tau_cell == "FS" ) &
        simpar % tau_cell = simpar % tau_cell * units % fs

        IF ( simpar % units % tau == "1/CM" ) &
        simpar % tau_nhc = simpar % tau_nhc * units % wavenum
        IF ( simpar % units % tau == "FS" ) &
        simpar % tau_nhc = simpar % tau_nhc * units % fs

        CALL time_convert ( simpar % tau_cell )
        CALL time_convert ( simpar % tau )
     END SELECT
     IF ( simpar % ensemble_coef == 'NVT' ) THEN
        IF ( simpar % units % tau == "1/CM" ) &
        simpar % tau_nhc_coef = simpar % tau_nhc_coef * units % wavenum
        IF ( simpar % units % tau == "FS" ) &
        simpar % tau_nhc_coef = simpar % tau_nhc_coef * units % fs

        CALL time_convert ( simpar % tau_nhc_coef )
     END IF
  END IF

  IF ( PRESENT ( part ) ) THEN
! convert position
     IF ( unit_char == "ANGSTROM" ) THEN
       DO i = 1, SIZE ( part )
          part(i) %r(1) = part(i) %r(1) * units % angstrom
          part(i) %r(2) = part(i) %r(2) * units % angstrom
          part(i) %r(3) = part(i) %r(3) * units % angstrom
       END DO
     ENDIF
  END IF

  IF ( PRESENT ( potparm ) ) THEN
! convert LJ parameters

     IF ( unit_char == "KCAL" )
        potparm ( :, : ) % lj % epsilon = potparm ( :, : ) % lj % epsilon * units % kcal
        potparm ( :, : ) % lj % sigma6 = potparm ( :, : ) % lj % sigma6 * ( units % angstrom ) ** 6
        potparm ( :, : ) % lj % sigma12 = potparm ( :, : ) % lj % sigma12 * ( units % angstrom ) ** 12
     ELSEIF ( unit_char == "KELVIN" )
        potparm ( :, : ) % lj % epsilon = potparm ( :, : ) % lj % epsilon * units % kelvin
        potparm ( :, : ) % lj % sigma6 = potparm ( :, : ) % lj % sigma6 * ( units % angstrom ) ** 6
        potparm ( :, : ) % lj % sigma12 = potparm ( :, : ) % lj % sigma12 * ( units % angstrom ) ** 12
     ENDIF 

! convert Williams parameters
     IF ( unit_char == "KCAL" )
       potparm ( :, : ) %willis%a = potparm ( :, : ) %willis%a*units%kcal
       potparm ( :, : ) %willis%b = potparm ( :, : ) %willis%b/units%angstrom
       potparm ( :, : ) %willis%c = potparm ( :, : ) %willis%c*units%kcal &
                                                      units%angstrom**6
     ELSEIF ( unit_char == "KELVIN" )
       potparm ( :, : ) %willis%a = potparm ( :, : ) %willis%a*units%kelvin
       potparm ( :, : ) %willis%b = potparm ( :, : ) %willis%b/units%angstrom
       potparm ( :, : ) %willis%c = potparm ( :, : ) %willis%c*units%kelvin &
                                                      units%angstrom**6
     ENDIF

! convert Cutoff distances
     IF ( unit_char /= "AU" ) &
     potparm ( :, : ) %rcutsq = potparm ( :, : ) % rcutsq * ( units % angstrom ) ** 2
  END IF
!
! note: e_cutoff_coul and energy_cutoff do not need to be converted
! because ewald_correction and the initialization of the energy_cutoff
! are called later and the cut-off will be computed in the appropriate units
!
! convert electrostatic parameters
  IF ( PRESENT ( ewald_param ) ) THEN
     ewald_param % eps0 = permittivity * units % eps0
     IF ( ewald_param % unit == "ANGSTROM" ) &
     ewald_param % alpha = ewald_param % alpha / units % angstrom
  END IF

! converting bond, bends and torsion and constraint parameters
  IF ( PRESENT ( intra_param ) ) THEN
    IF ( ASSOCIATED ( intra_param % bond_param )) THEN
       IF ( unit_char == "KCAL" ) THEN
         intra_param % bond_param ( : ) % r0 = intra_param % bond_param ( : ) % r0 &
                                               *units%angstrom
         intra_param % bond_param ( : ) % k = intra_param % bond_param ( : ) % k &
                                               *units%kcal/units%angstrom/units%angstrom
       ELSEIF ( unit_char == "KELVIN" ) THEN
         intra_param % bond_param ( : ) % r0 = intra_param % bond_param ( : ) % r0 &
                                               *units%angstrom
         intra_param % bond_param ( : ) % k = intra_param % bond_param ( : ) % k &
                                               *units%kelvin/units%angtrom/units%angstrom
       ENDIF
    END IF

    IF ( ASSOCIATED ( intra_param % bend_param )) THEN
       IF ( unit_char == "KCAL" ) THEN
         intra_param % bend_param( : ) % k &
         = intra_param % bend_param( : ) % k * units%kcal
       ELSEIF ( unit_char == "KCAL" ) THEN
         intra_param % bend_param( : ) % k &
         = intra_param % bend_param( : ) % k * units%kelvin
       ENDIF
    END IF

    IF ( ASSOCIATED ( intra_param % torsion_param ) ) THEN
       IF ( unit_char == "KCAL" ) THEN
         DO i = 0, SIZE ( intra_param % torsion_param(1) % a ) - 1
            intra_param % torsion_param( : ) % a(i) =  &
            intra_param % torsion_param( : ) % a(i)*units%kcal
         END DO
       ELSEIF ( unit_char == "KELVIN" ) THEN
         DO i = 0, SIZE ( intra_param % torsion_param(1) % a ) - 1
            intra_param % torsion_param( : ) % a(i) =  &
            intra_param % torsion_param( : ) % a(i)*units%kelvin
         END DO
       ENDIF
    END IF

    IF ( ASSOCIATED ( intra_param % constraint_distance ) ) THEN
       IF ( ( unit_char == "KELVIN" ) .OR. ( unit_char == "KCAL" ) ) &
       intra_param % constraint_distance ( : ) =  &
       intra_param % constraint_distance ( : ) * units%angstrom
    END IF
  ENDIF

END SUBROUTINE convert_to_cp2k_units

!******************************************************************************

END MODULE cp2k_units

!******************************************************************************
