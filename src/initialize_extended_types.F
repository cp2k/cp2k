!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/initialize_extended_types [1.0] *
!!
!!   NAME
!!     initialize_extended_types
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     cjm, FEB 20 2001:  added subroutine assign_extended_parameters
!!     cjm, MAY 03 2001:  reorganized and added separtate routines for
!!                        nhc_part, nhc_baro, nhc_ao, npt
!!   SOURCE
!******************************************************************************

MODULE initialize_extended_types

  USE atoms_input, ONLY : system_type
  USE coefficient_types, ONLY : coeff_type
  USE global_types, ONLY : global_environment_type
  USE kinds, ONLY : dbl
  USE md, ONLY : simulation_parameters_type
  USE molecule_types, ONLY : molecule_structure_type, molecule_type
  USE message_passing, ONLY : mp_sum
  USE nhc, ONLY : nhc_particle_setup, nhc_barostat_setup, nhc_pol_setup
  USE nose, ONLY : initialize_nhc, lnhc_parameters_type, initialize_npt, &
                   npt_info_type
  USE simulation_cell, ONLY : cell_type
  USE termination, ONLY : stop_program

  IMPLICIT NONE

  PUBLIC :: initialize_nhc_part, initialize_nhc_baro, initialize_npt_type, &
            initialize_nhc_pol, assign_extended_parameters
  PRIVATE

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** initialize_extended_types/initialize_npt_type [1.0] *
!!
!!   NAME
!!     initialize_npt_type
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!    
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE initialize_npt_type ( simpar, globenv, npt_info, box )

  IMPLICIT NONE

! Arguments
  TYPE ( simulation_parameters_type ), INTENT ( INOUT ) :: simpar
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv
  TYPE ( npt_info_type ), DIMENSION ( :, : ), POINTER :: npt_info
  TYPE ( cell_type ), INTENT (inout) :: box
! locals
  INTEGER :: isos

!------------------------------------------------------------------------------

! first allocating the npt_info_type
  SELECT CASE ( simpar % ensemble )
  CASE ( 'NPT_I')
     ALLOCATE ( npt_info ( 1, 1 ), STAT = isos )
     IF ( isos /= 0 ) CALL stop_program ( 'initial_npt', &
          'failed to allocate npt_info')

  CASE ( 'NPT_F')
     ALLOCATE ( npt_info ( 3, 3 ), STAT = isos )
     IF ( isos /= 0 ) CALL stop_program ( 'initial_npt', &
          'failed to allocate npt_info')

  CASE  DEFAULT
     ALLOCATE ( npt_info ( 0, 0 ), STAT = isos )
     IF ( isos /= 0 ) CALL stop_program ( 'initial_npt', &
          'failed to allocate npt_info')
  END SELECT

  IF ( simpar % read_type /= 'ALL' ) &
  CALL initialize_npt ( npt_info, simpar % tau_cell, simpar % temp_ext,  &
                     box % deth, simpar % nfree, simpar % ensemble, globenv )

END SUBROUTINE initialize_npt_type

!******************************************************************************
!!****** initialize_extended_types/initialize_nhc_pol [1.0] *
!!
!!   NAME
!!     initialize_nhc_pol
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!    
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE initialize_nhc_pol ( simpar, coeff, globenv, nhc )

  IMPLICIT NONE

! Arguments
  TYPE ( coeff_type ), INTENT ( IN ) :: coeff
  TYPE ( simulation_parameters_type ), INTENT ( INOUT ) :: simpar
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv
  TYPE ( lnhc_parameters_type ),  INTENT ( INOUT ) :: nhc
! locals
  INTEGER :: isos

!------------------------------------------------------------------------------

! fire up the thermostats
  CALL nhc_pol_setup ( simpar, coeff, nhc, globenv )

  IF ( simpar % pol_read_type /= 'ALL' ) THEN

! Initializing thermostat forces and velocities for the Nose-Hoover
! Chain variables
      IF ( nhc % nhc_len /= 0 ) CALL initialize_nhc ( nhc, simpar % pol_tau_nhc,  &
                                      simpar % pol_temp, globenv )
  END IF

END SUBROUTINE initialize_nhc_pol

!******************************************************************************
!!****** initialize_extended_types/initialize_nhc_baro [1.0] *
!!
!!   NAME
!!     initialize_nhc_baro
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!    
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE initialize_nhc_baro ( simpar, globenv, nhc )

  IMPLICIT NONE

! Arguments
  TYPE ( simulation_parameters_type ), INTENT ( INOUT ) :: simpar
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv
  TYPE ( lnhc_parameters_type ),  INTENT ( INOUT ) :: nhc
! locals
  INTEGER :: isos

!------------------------------------------------------------------------------

! fire up the thermostats, if NPT
  CALL nhc_barostat_setup ( simpar, nhc, globenv )

  IF ( simpar % read_type /= 'ALL' ) THEN

! Initializing thermostat forces and velocities for the Nose-Hoover
! Chain variables
      IF ( nhc % nhc_len /= 0 ) &
       CALL initialize_nhc ( nhc, simpar % tau_nhc,  &
                       simpar % temp_ext, globenv )
  END IF

END SUBROUTINE initialize_nhc_baro

!******************************************************************************
!!****** initialize_extended_types/initialize_nhc_part [1.0] *
!!
!!   NAME
!!     initialize_nhc_part
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE initialize_nhc_part ( box, simpar, molecule, mol_setup, &
     globenv, nhc )

  IMPLICIT NONE

! Arguments
  TYPE ( simulation_parameters_type ), INTENT ( INOUT ) :: simpar
  TYPE ( molecule_structure_type ), DIMENSION ( : ), INTENT ( INOUT ) :: &
       molecule
  TYPE ( molecule_type ), DIMENSION ( : ), INTENT ( INOUT ) :: mol_setup
  TYPE ( cell_type ), INTENT ( IN ) :: box
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv
  TYPE ( lnhc_parameters_type ),  INTENT ( INOUT ) :: nhc
! locals
  INTEGER :: isos

!------------------------------------------------------------------------------

! Initialize simpar%nfree to dim (because at this point no
! assumptions have been made about thermostatting options
! and COM is ALWAYS zeroed in this code.  If there are
! multiple thermostats etc...simpar % nfree will be set
! to zero in nhc_particle_setup.f)

  simpar % nfree = 6 - SUM ( box % perd ( 1:3 ) )

! fire up the thermostats, if not NVE
  CALL nhc_particle_setup ( simpar, molecule, mol_setup, nhc, globenv )

  IF ( simpar % read_type /= 'ALL' ) THEN

! Initializing thermostat forces and velocities for the Nose-Hoover
! Chain variables
     IF ( nhc % nhc_len /= 0 ) &
       CALL initialize_nhc ( nhc, simpar % tau_nhc,  &
                        simpar % temp_ext, globenv )

  ENDIF

! compute the total number of degrees of freedom for temperature


  CALL get_nfree ( simpar, molecule, globenv )

END SUBROUTINE initialize_nhc_part

!!*****
!******************************************************************************
!!****** initialize_extended_types/assign_extended_parameters [1.0] *
!!
!!   NAME
!!     assign_extended_parameters
!!
!!   FUNCTION
!!    Assigns extended parameters from the restart file.
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************
SUBROUTINE assign_extended_parameters ( deth, ensemble, ainp, npt,  &
                                        nhc_part, nhc_baro )
  IMPLICIT NONE
  CHARACTER ( len = * ), INTENT ( IN ) :: ensemble
  TYPE ( system_type ), INTENT ( IN ) :: ainp
  TYPE ( lnhc_parameters_type ), INTENT ( OUT ) :: nhc_baro
  TYPE ( lnhc_parameters_type ), INTENT ( OUT ) :: nhc_part
  TYPE ( npt_info_type ), DIMENSION ( :, : ), INTENT ( OUT ) :: npt
  REAL ( dbl ), INTENT ( IN ) :: deth
! locals
  INTEGER :: size1, size2, size3, size4, i, j
! assign nhc_part % eta
  size1 = SIZE ( nhc_part % nvt, 1 )
  size2 = SIZE ( nhc_part % nvt, 2 )
  size3 = SIZE ( ainp % eta_part, 1 )
  size4 = SIZE ( ainp % eta_part, 2 )
  IF  ( ( size1 /= 0 ) .AND. ( size2 /= 0 ) ) THEN
    IF ( ( size1 == size3 ) .AND. ( size2 == size4 ) ) THEN
      nhc_part % nvt ( :, : ) % eta = ainp % eta_part ( :, : )
    ELSE
      CALL stop_program ( 'assign_extended_parameters', 'memory mismatch' )
    ENDIF
  ENDIF
! assign nhc_part % v
  size3 = SIZE ( ainp % veta_part, 1 )
  size4 = SIZE ( ainp % veta_part, 2 )
  IF ( ( size1 /= 0 ) .AND. ( size2 /= 0 ) ) THEN
    IF ( ( size1 == size3 ) .AND. ( size2 == size4 ) ) THEN
      nhc_part % nvt ( :, : ) % v = ainp % veta_part ( :, : )
    ELSE
      CALL stop_program ( 'assign_extended_parameters', 'memory mismatch' )
    ENDIF
  ENDIF
! assign nhc_part % mass
  size3 = SIZE ( ainp % mass_eta_part, 1 )
  size4 = SIZE ( ainp % mass_eta_part, 2 )
  IF  ( ( size1 /= 0 ) .AND. ( size2 /= 0 ) ) THEN
    IF ( ( size1 == size3 ) .AND. ( size2 == size4 ) ) THEN
      nhc_part % nvt ( :, : ) % mass = ainp % mass_eta_part ( :, : )
    ELSE
      CALL stop_program ( 'assign_extended_parameters', 'memory mismatch' )
    ENDIF
  ENDIF
! assign the forces
  DO i = 1, size2
     DO j = 2, size1
        nhc_part % nvt ( j, i ) % f = nhc_part % nvt ( j-1, i ) % mass * &
                                      nhc_part % nvt ( j-1, i ) % v *  &
                                      nhc_part % nvt ( j-1, i ) % v -  &
                                      nhc_part % nvt ( j, i ) % nkt
        nhc_part % nvt ( j, i ) % f = nhc_part % nvt ( j, i ) % f &
                                     /nhc_part % nvt ( j, i ) % mass
     END DO
  END DO
! assign nhc_baro % eta
  size1 = SIZE ( nhc_baro % nvt, 1 )
  size2 = SIZE ( nhc_baro % nvt, 2 )
  size3 = SIZE ( ainp % eta_baro, 1 )
  size4 = SIZE ( ainp % eta_baro, 2 )
  IF ( ( size1 /= 0 ) .AND. ( size2 /= 0 ) ) THEN
    IF ( ( size1 == size3 ) .AND. ( size2 == size4 ) ) THEN
      nhc_baro % nvt ( :, : ) % eta = ainp % eta_baro ( :, : )
    ELSE
      CALL stop_program ( 'assign_extended_parameters', 'memory mismatch' )
    ENDIF
  ENDIF
! assign nhc_part % v
  size3 = SIZE ( ainp % veta_baro, 1 )
  size4 = SIZE ( ainp % veta_baro, 2 )
  IF ( ( size1 /= 0 ) .AND. ( size2 /= 0 ) ) THEN
    IF ( ( size1 == size3 ) .AND. ( size2 == size4 ) ) THEN
      nhc_baro % nvt ( :, : ) % v = ainp % veta_baro ( :, : )
    ELSE
      CALL stop_program ( 'assign_extended_parameters', 'memory mismatch' )
    ENDIF
  ENDIF
! assign nhc_baro % mass
  size3 = SIZE ( ainp % mass_eta_baro, 1 )
  size4 = SIZE ( ainp % mass_eta_baro, 2 )
  IF ( ( size1 /= 0 ) .AND. ( size2 /= 0 ) ) THEN
    IF ( ( size1 == size3 ) .AND. ( size2 == size4 ) ) THEN
      nhc_baro % nvt ( :, : ) % mass = ainp % mass_eta_baro ( :, : )
    ELSE
      CALL stop_program ( 'assign_extended_parameters', 'memory mismatch' )
    ENDIF
  ENDIF
! assign the forces
  DO i = 1, size2
     DO j = 2, size1
        nhc_baro % nvt ( j, i ) % f = nhc_baro % nvt ( j-1, i ) % mass * &
                                      nhc_baro % nvt ( j-1, i ) % v *  &
                                      nhc_baro % nvt ( j-1, i ) % v -  &
                                      nhc_baro % nvt ( j, i ) % nkt
        nhc_baro % nvt ( j, i ) % f = nhc_baro % nvt ( j, i ) % f &
                                     /nhc_baro % nvt ( j, i ) % mass
     END DO
  END DO

! assign npt % v
  size1 = SIZE ( npt, 1 )
  size2 = SIZE ( npt, 2 )
  size3 = SIZE ( ainp % veps, 1 )
  size4 = SIZE ( ainp % veps, 2 )
  IF ( ( size1 /= 0 ) .AND. ( size2 /= 0 ) ) THEN
    IF ( ( size1 == size3 ) .AND. ( size2 == size4 ) ) THEN
      npt ( :, : ) % v = ainp % veps ( :, : )
    ELSE
      CALL stop_program ( 'assign_extended_parameters', 'memory mismatch' )
    ENDIF
  ENDIF
! assign npt % mass
  size3 = SIZE ( ainp % mass_eps, 1 )
  size4 = SIZE ( ainp % mass_eps, 2 )
  IF ( ( size1 /= 0 ) .AND. ( size2 /= 0 ) ) THEN
    IF ( ( size1 == size3 ) .AND. ( size2 == size4 ) ) THEN
      npt ( :, : ) % mass = ainp % mass_eps ( :, : )
    ELSE
      CALL stop_program ( 'assign_extended_parameters', 'memory mismatch' )
    ENDIF
  ENDIF
! assign npt % eps
  IF ( ensemble == 'NPT_I' ) npt ( :, : ) % eps = LOG ( deth ) / 3.0_dbl


END SUBROUTINE assign_extended_parameters

!!*****
!******************************************************************************
!!****** initialize_extended_types/get_nfree [1.0] *
!!
!!   NAME
!!     get_nfree
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     29-Nov-00 (JGH) correct counting of DOF if constraints are off
!!     18-Dec-00 (JGH) bug fix in initialization of nconst
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE get_nfree ( simpar, molecule, globenv )

  IMPLICIT NONE

! Arguments
  TYPE ( simulation_parameters_type ), INTENT ( INOUT ) :: simpar
  TYPE ( molecule_structure_type ), INTENT ( IN ), DIMENSION ( : ) :: molecule
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv

! Locals
  INTEGER :: imol, nconst, natoms, nfree, iw

!------------------------------------------------------------------------------

  natoms = 0
  nconst = 0
  DO imol = 1, SIZE ( molecule )
     natoms = natoms + molecule ( imol ) % natoms_mol
  END DO
  IF ( simpar % constraint ) THEN
    DO imol = 1, SIZE ( molecule )
       nconst = nconst + molecule ( imol ) % ndistconst_mol
       nconst = nconst + 3 * molecule ( imol ) % ng3x3const_mol
    END DO
  END IF

#if defined(__parallel)
  CALL mp_sum(nconst,globenv % group)
  CALL mp_sum(natoms,globenv % group)
#endif

! definining simpar%nfree.  If multiple thermostats present
! it is initialized to 0.  It one thermostat or no thermostats
! it is initialized to 3
  nfree = 3 * natoms - nconst - simpar % nfree

  IF ( globenv % ionode .AND. globenv % print_level > 2 ) THEN
     iw = globenv % scr
     WRITE ( iw, '( /, A )' ) ' Calculation of degrees of freedom'
     WRITE ( iw, '( T48, A, T71, I10 )' ) '       Number of atoms:', natoms
     WRITE ( iw, '( T48, A, T71, I10 )' ) ' Number of constraints:', nconst
     WRITE ( iw, '( T34, A, T71, I10 )' ) &
          ' Invariants(translation + rotations):', simpar % nfree
     WRITE ( iw, '( T48, A, T71, I10 )' ) '    Degrees of freedom:', nfree
     WRITE ( iw, '()' )
  END IF
  simpar % nfree = nfree
END SUBROUTINE get_nfree

!!*****
!******************************************************************************

END MODULE initialize_extended_types

!******************************************************************************
