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

  USE ewald_parameters_types, ONLY : ewald_parameters_type
  USE kinds, ONLY : dbl
  USE md, ONLY : simulation_parameters_type
  USE molecule_types, ONLY : intra_parameters_type
  USE pair_potential, ONLY : potentialparm_type
  USE particle_types, ONLY : particle_type
  USE physcon, ONLY : joule, evolt, boltzmann
  USE simulation_cell, ONLY : cell_type
  USE unit, ONLY : unit_convert_type

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: convert

!!*****
!-----------------------------------------------------------------------------!

CONTAINS

!-----------------------------------------------------------------------------!
!
! routine to convert input variables from external units to internal units
SUBROUTINE convert ( units, simpar, part, box, potparm, &
     intra_param, ewald_param )

  IMPLICIT NONE

! Arguments
  TYPE (unit_convert_type ), INTENT ( IN ) :: units
  TYPE (simulation_parameters_type ), INTENT ( INOUT ), OPTIONAL :: simpar
  TYPE (particle_type ), INTENT ( INOUT ), OPTIONAL :: part ( : )
  TYPE (cell_type ), INTENT ( INOUT ), OPTIONAL :: box
  TYPE (intra_parameters_type ), INTENT ( INOUT ), OPTIONAL :: intra_param
  TYPE (ewald_parameters_type ), INTENT ( INOUT ), OPTIONAL :: ewald_param
  TYPE (potentialparm_type ), INTENT ( INOUT ), OPTIONAL :: potparm ( :, : )
!  TYPE (particle_prop_type ), INTENT ( INOUT ), OPTIONAL :: pstat ( : )

! Locals
  REAL ( dbl ), PARAMETER :: permittivity = 8.854187817E-12_dbl
  REAL ( dbl ) :: conv2k
  INTEGER :: i

!------------------------------------------------------------------------------

! convert hmat
  IF ( PRESENT ( box)) box%hmat = box%hmat*units%lconv

  IF ( PRESENT ( simpar)) THEN
! convert time
     simpar%dt = simpar%dt*units%tconv
! convert temperature
     simpar%temp_ext = simpar%temp_ext*units%econv
! convert temperature of coefficents
     simpar%temp_coef = simpar%temp_coef*units%econv
! convert pressure
     simpar%p_ext = simpar%p_ext/units%pconv
! convert verlet skin
     simpar%verlet_skin = simpar%verlet_skin*units%lconv

! convert omega_p (wp) of PIMD
     simpar%pimd_params%wp = simpar%pimd_params%wp * units%econv / units%tconv
     
! converting extended system time constants
     SELECT CASE (simpar%ensemble)
     CASE ( 'NVT')
        simpar%tau_nhc = simpar%tau_nhc*units%tconv
     CASE ( 'NPT_I')
        simpar%tau_cell = simpar%tau_cell*units%tconv
        simpar%tau_nhc = simpar%tau_nhc*units%tconv
     CASE ( 'NPT_F')
        simpar%tau_cell = simpar%tau_cell*units%tconv
        simpar%tau_nhc = simpar%tau_nhc*units%tconv
     END SELECT
     IF (simpar % ensemble_coef =='NVT') THEN
        simpar%tau_nhc_coef = simpar%tau_nhc_coef*units%tconv
     END IF
  END IF

!! gt units to be fixed !!
!!  IF ( PRESENT ( pstat)) THEN
! convert mass
!!     DO i = 1, size(pstat)
!!        pstat(i) %mass = pstat(i) %mass*units%mconv
!!     END DO
!!  END IF

  IF ( PRESENT ( part)) THEN
! convert position
     DO i = 1, size(part)
        part(i) %r(1) = part(i) %r(1)*units%lconv
        part(i) %r(2) = part(i) %r(2)*units%lconv
        part(i) %r(3) = part(i) %r(3)*units%lconv
     END DO
  END IF

  IF ( PRESENT ( potparm)) THEN
! convert LJ parameters
     potparm ( :, : ) %lj%epsilon = potparm ( :, : ) %lj%epsilon*units%econv
     potparm ( :, : ) %lj%sigma6 = potparm ( :, : ) %lj%sigma6*(units%lconv)**6
     potparm ( :, : ) %lj%sigma12 = potparm ( :, : ) %lj%sigma12*(units%lconv)**12
! convert Williams parameters
     potparm ( :, : ) %willis%a = potparm ( :, : ) %willis%a*units%econv
     potparm ( :, : ) %willis%b = potparm ( :, : ) %willis%b/units%lconv
     potparm ( :, : ) %willis%c = potparm ( :, : ) %willis%c*units%econv* &
          units%lconv**6
! convert Goodwin parameters ( ev -> K has to be included )
     conv2k = joule/(evolt*boltzmann)
     potparm(:,:) %goodwin%vr0 = potparm(:,:) %goodwin%vr0*units%econv* &
          conv2k
     potparm(:,:) %goodwin%d = potparm(:,:) %goodwin%d*units%lconv
     potparm(:,:) %goodwin%dc = potparm(:,:) %goodwin%dc*units%lconv
! convert Cutoff distances
     potparm ( :, : ) %rcutsq = potparm ( :, : ) %rcutsq*(units%lconv)**2
  END IF
! convert hopping matrix elements
! orbital energies
!       DO i = 1, size(tbatom)
!         tbatom(i) %orbital_energy(:,:) = tbatom(i) %orbital_energy(:,:)* &
!           units%econv*conv_2k
!       END DO
! convert Goodwin parameters
!       tbhop(:,:) %rcutsq = tbhop(:,:) %rcutsq*(units%lconv)**2
!       DO i = 1, size(tbhop,1)
!         DO j = 1, size(tbhop,2)
!           tbhop(i,j) %gsp(:) %er0 = tbhop(i,j) %gsp(:) %er0*units%econv* &
!             conv_2k
!           tbhop(i,j) %gsp(:) %r0 = tbhop(i,j) %gsp(:) %r0*units%lconv
!           tbhop(i,j) %gsp(:) %rt = tbhop(i,j) %gsp(:) %rt*units%lconv
!         END DO
!       END DO
!
! note: e_cutoff_coul and energy_cutoff do not need to be converted
! because ewald_correction and the initialization of the energy_cutoff
! are called later and the cut-off will be computed in the appropriate units
!
! convert electrostatic parameters
  IF ( PRESENT ( ewald_param)) THEN
     ewald_param % eps0 = permittivity*units%eps0
     ewald_param % alpha = ewald_param % alpha/units % lconv
  END IF

! converting bond, bends and torsion and constraint parameters
  IF ( PRESENT ( intra_param ) ) THEN
    IF ( ASSOCIATED ( intra_param % bond_param )) THEN
       intra_param % bond_param ( : ) % r0 = intra_param % bond_param ( : ) % r0 &
            *units%lconv
       intra_param % bond_param ( : ) % k = intra_param % bond_param ( : ) % k &
            * units%econv/units%lconv/units%lconv
    END IF

    IF ( ASSOCIATED ( intra_param % bend_param )) THEN
       intra_param % bend_param( : ) % k &
            = intra_param % bend_param( : ) % k * units%econv
    END IF

    IF ( ASSOCIATED ( intra_param % torsion_param )) THEN
       DO i = 0, SIZE ( intra_param % torsion_param(1) % a ) - 1
          intra_param % torsion_param( : ) % a(i) =  &
               intra_param % torsion_param( : ) % a(i)*units%econv
       END DO
    END IF

    IF ( ASSOCIATED ( intra_param % constraint_distance )) THEN
       intra_param % constraint_distance ( : ) &
            = intra_param % constraint_distance ( : ) * units%lconv
    END IF
  ENDIF

END SUBROUTINE convert

!******************************************************************************

END MODULE convert_units

!******************************************************************************
