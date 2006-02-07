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

  USE cell_types,                      ONLY: cell_type
  USE cp2k_units,                      ONLY: get_cp2k_units,&
                                             unit_convert_type
  USE kinds,                           ONLY: dp
  USE particle_types,                  ONLY: particle_type

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
                                   tau_nhc, mass, angle, eps0, l_power, t_power, part, box, &
                                   energy )


    CHARACTER(len=*), INTENT(IN), OPTIONAL   :: unit_char
    REAL(KIND=dp), INTENT(INOUT), OPTIONAL   :: time, temperature, pressure, &
                                                length, tau_nhc, mass, angle, &
                                                eps0
    INTEGER, INTENT(IN), OPTIONAL            :: l_power, t_power
    TYPE(particle_type), INTENT(INOUT), &
      OPTIONAL                               :: part( : )
    TYPE(cell_type), INTENT(INOUT), OPTIONAL :: box
    REAL(KIND=dp), INTENT(INOUT), OPTIONAL   :: energy

    INTEGER                                  :: i, ll, tt
    TYPE(unit_convert_type), POINTER         :: units

!------------------------------------------------------------------------------

  CALL get_cp2k_units ( units )

  ll = 1 
  tt = 1
  IF ( PRESENT ( l_power ) ) ll = l_power 
  IF ( PRESENT ( t_power ) ) tt = t_power 

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
     time = time * units % fs ** tt
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

! converting bond, bends and torsion and constraint parameters
  IF ( PRESENT ( energy ) ) THEN
     IF ( unit_char == 'KCAL' ) THEN
       energy  = energy * units%kcal
     ELSEIF ( unit_char == 'KELVIN' ) THEN
       energy  = energy * units%kelvin
     ELSEIF ( unit_char == 'EV' ) THEN
       energy  = energy * units%ev
     ENDIF
  END IF


END SUBROUTINE convert_to_cp2k_units

!******************************************************************************

END MODULE convert_units

!******************************************************************************
