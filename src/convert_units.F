!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/convert_units [1.0] *
!!
!!   NAME
!!     convert_units
!!
!!   FUNCTION
!!
!!   AUTHOR
!!
!!   MODIFICATION HISTORY
!!     JGH (26-01-2002) : force field parameters are now stored in tables
!!
!!   SOURCE
!******************************************************************************

MODULE convert_units

  USE cp2k_units,                      ONLY: get_cp2k_units,&
                                             unit_convert_type
  USE kinds,                           ONLY: dbl
  USE pair_potential_types,            ONLY: pair_potential_type
  USE particle_types,                  ONLY: particle_type
  USE simulation_cell,                 ONLY: cell_type

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: convert_to_cp2k_units

!!*****
!-----------------------------------------------------------------------------!

CONTAINS

!******************************************************************************
!!****** convert_units/convert_to_cp2k_units [1.0] *
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
SUBROUTINE convert_to_cp2k_units ( unit_char, time, temperature, pressure, length,  &
                                   tau_nhc, mass, angle, eps0, l_power, part, box, &
                                   energy, ff_nonbonded )


    CHARACTER(len=*), INTENT(IN), OPTIONAL   :: unit_char
    REAL(dbl), INTENT(INOUT), OPTIONAL       :: time, temperature, pressure, &
                                                length, tau_nhc, mass, eps0, energy,angle
    INTEGER, INTENT(IN), OPTIONAL            :: l_power
    TYPE(particle_type), INTENT(INOUT), &
      OPTIONAL                               :: part( : )
    TYPE(cell_type), INTENT(INOUT), OPTIONAL :: box
    TYPE(pair_potential_type), &
      POINTER,   OPTIONAL                    :: ff_nonbonded ( : )
    INTEGER                                  :: i, ll
    TYPE(unit_convert_type), POINTER         :: units

!------------------------------------------------------------------------------

  CALL get_cp2k_units ( units )

  ll = 1 
  IF ( PRESENT ( l_power ) ) ll = l_power 

! convert hmat
  IF ( PRESENT ( box ) ) THEN
   IF ( box % unit_of_length_name == 'ANGSTROMS' ) &
   box % hmat = box % hmat * units % angstrom
  ENDIF

  IF ( PRESENT ( angle ) ) THEN
! convert time
     IF ( unit_char == 'DEGREE' ) &
     angle = angle * units % degree
  ENDIF

  IF ( PRESENT ( time ) ) THEN
! convert time
     IF ( unit_char == 'FS' ) &
     time = time * units % fs
  ENDIF

  IF ( PRESENT ( temperature ) ) THEN
! convert temperature
     temperature = temperature * units % kelvin
  ENDIF

  IF ( PRESENT ( pressure ) ) THEN
! convert pressure
     pressure =  pressure / units % bar
  ENDIF 

  IF ( PRESENT ( length ) ) THEN
     IF ( unit_char == 'ANGSTROM' ) &
     length = length * units % angstrom ** ll
     IF ( unit_char == 'BOHR' ) &
     length = length * units % bohr ** ll
  ENDIF

  IF ( PRESENT ( mass ) ) THEN
   mass = mass * units % mass
  ENDIF

  IF ( PRESENT ( eps0 ) ) THEN
   eps0 = units % eps0
  ENDIF

  IF ( PRESENT ( tau_nhc ) ) THEN
     IF ( unit_char == '1/CM' ) &
     tau_nhc =  units % wavenum / tau_nhc
     IF ( unit_char == 'FS' ) &
     tau_nhc = tau_nhc * units % fs
  ENDIF 

!!!!
! convert omega_p (KIND = dp) of PIMD
!     simpar%pimd_params%wp = simpar%pimd_params%wp * units%econv / units % tconv
!!!!
     
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

  IF ( PRESENT ( ff_nonbonded ) ) THEN
   IF(ASSOCIATED(ff_nonbonded)) THEN
! convert LJ parameters

     IF ( unit_char == 'KCAL' ) THEN
        ff_nonbonded ( : ) % lj % epsilon = ff_nonbonded ( : ) % lj % epsilon * units % kcal
        ff_nonbonded ( : ) % lj % sigma6 = ff_nonbonded ( : ) % lj % sigma6 * ( units % angstrom ) ** 6
        ff_nonbonded ( : ) % lj % sigma12 = ff_nonbonded ( : ) % lj % sigma12 * ( units % angstrom ) ** 12
     ELSEIF ( unit_char == 'KELVIN' ) THEN
        ff_nonbonded ( : ) % lj % epsilon = ff_nonbonded ( : ) % lj % epsilon * units % kelvin
        ff_nonbonded ( : ) % lj % sigma6 = ff_nonbonded ( : ) % lj % sigma6 * ( units % angstrom ) ** 6
        ff_nonbonded ( : ) % lj % sigma12 = ff_nonbonded ( : ) % lj % sigma12 * ( units % angstrom ) ** 12
     ENDIF 

! convert Williams parameters
     IF ( unit_char == 'KCAL' ) THEN
       ff_nonbonded ( : ) %willis%a = ff_nonbonded ( : ) %willis%a*units%kcal
       ff_nonbonded ( : ) %willis%b = ff_nonbonded ( : ) %willis%b/units%angstrom
       ff_nonbonded ( : ) %willis%c = ff_nonbonded ( : ) %willis%c*units%kcal* &
                                                      units%angstrom**6
     ELSEIF ( unit_char == 'KELVIN' ) THEN
       ff_nonbonded ( : ) %willis%a = ff_nonbonded ( : ) %willis%a*units%kelvin
       ff_nonbonded ( : ) %willis%b = ff_nonbonded ( : ) %willis%b/units%angstrom
       ff_nonbonded ( : ) %willis%c = ff_nonbonded ( : ) %willis%c*units%kelvin* &
                                                      units%angstrom**6
     ENDIF

! convert Cutoff distances
     IF ( ( unit_char /= 'ATOMIC' ) .OR. ( unit_char /= 'AU' ) ) &
     ff_nonbonded ( : ) % rcutsq = ff_nonbonded ( : ) % rcutsq * ( units % angstrom ) ** 2
    END IF
  END IF

! converting bond, bends and torsion and constraint parameters
  IF ( PRESENT ( energy ) ) THEN
     IF ( unit_char == 'KCAL' ) THEN
       energy  = energy * units%kcal
     ELSEIF ( unit_char == 'KELVIN' ) THEN
       energy  = energy * units%kelvin
     ENDIF
  END IF


END SUBROUTINE convert_to_cp2k_units

!******************************************************************************

END MODULE convert_units

!******************************************************************************
