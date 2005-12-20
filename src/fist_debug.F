!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/fist_debug [1.0] *
!!
!!   NAME
!!     fist_debug
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!    CJM 2/01.  Passed in the variable box_ref
!!
!!   SOURCE
!******************************************************************************

MODULE fist_debug
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_assertion_failed
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_logger_type,&
                                             cp_warning_level
  USE atomic_kind_types,                        ONLY : atomic_kind_type
  USE distribution_1d_types,                    ONLY : distribution_1d_type
  USE ewald_environment_types,                  ONLY : ewald_environment_type,&
                                                       ewald_env_get
  USE ewald_pw_types,                           ONLY : ewald_pw_type
  USE fist_environment_types,                   ONLY : fist_environment_type,&
                                                       get_fist_env
  USE fist_force,                               ONLY : fist_force_control,&
                                                       debug_variables_type
  USE fist_force_numer,                         ONLY : force_bond_numer,&
                                                       force_bend_numer,&
                                                       force_torsion_numer,&
                                                       force_imptors_numer,&
                                                       force_nonbond_numer,&
                                                       force_onefour_numer,&
                                                       force_ub_numer,&
                                                       ewald_energy_numer,&
                                                       ewald_force_numer,&
                                                       pv_bond_numer,&
                                                       pv_bend_numer,&
                                                       pv_torsion_numer,&
                                                       pv_imptors_numer,&
                                                       pv_nonbond_numer,&
                                                       pv_onefour_numer,&
                                                       pv_ub_numer,&
                                                       pv_g_numer,&
                                                       pme_force_numer
  USE fist_nonbond_env_types,                   ONLY : fist_nonbond_env_type,&
                                                       fist_nonbond_env_get
  USE force_env_types,                          ONLY : force_env_type
  USE global_types,                             ONLY : global_environment_type
  USE kinds,                                    ONLY : dp,default_string_length
  USE message_passing,                          ONLY :  mp_sum
  USE molecule_kind_types,                      ONLY : molecule_kind_type
  USE molecule_types_new,                       ONLY : molecule_type
  USE pair_potential_types,                     ONLY : pair_potential_pp_type
  USE particle_types,                           ONLY : particle_type
  USE simulation_cell,                          ONLY : cell_type
  USE termination,                              ONLY : stop_memory, stop_program
  USE virial_types,                             ONLY : virial_type
  USE input_constants,                 ONLY: do_ewald_none,&
                                             do_ewald_ewald,&
                                             do_ewald_pme,&
                                             do_ewald_spme
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: debug_control

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** fist_debug/debug_control [1.0] *
!!
!!   NAME
!!     debug_control
!!
!!   FUNCTION
!!     Routine to execute numerical tests for debugging
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     Changed Formats (JGH) 22-Nov-2000
!!     Harald Forbert (Dec-2000): Changes for multiple linked lists
!!     Torsions added (DG) 05-Dec-2000
!!     Additional output information (DG) 05-Dec-2000
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE  debug_control ( force_env, globenv, error)

  IMPLICIT NONE
! Arguments
  TYPE ( force_env_type ), POINTER :: force_env
  TYPE ( global_environment_type ), pointer :: globenv
  TYPE(cp_error_type), INTENT(inout), &
       OPTIONAL                               :: error
! Locals
  INTEGER ::  iatom, iw, iwf, ir, nbond, nbend, ntorsion, group
  INTEGER, DIMENSION ( 2 ) :: dum
  INTEGER :: i, natoms, nnodes, isos, n, nlocal_particles
  INTEGER :: ewald_type
  TYPE (pair_potential_pp_type), POINTER :: potparm
  TYPE ( particle_type ), POINTER :: particle_set ( : )
  TYPE(distribution_1d_type), POINTER      :: local_molecules, &
                                              local_particles
  TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
  TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
  TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
  TYPE ( fist_nonbond_env_type ), POINTER    :: fist_nonbond_env
  TYPE ( cell_type ), POINTER                :: cell
  TYPE ( ewald_pw_type ),   POINTER          :: ewald_pw
  TYPE ( ewald_environment_type ),   POINTER :: ewald_env
  TYPE ( debug_variables_type ) :: dbg
  TYPE ( fist_environment_type ), POINTER :: fist_env
  REAL (KIND=dp) :: delta, energy_numer, cutoff, ecut
  REAL (KIND=dp) :: e_numer, pv_test ( 3, 3 )
  REAL (KIND=dp) :: err1, numer, denom1, vec ( 3 ), e_bc, e_real, energy_tot
  REAL (KIND=dp) :: denom2, err2
  REAL (KIND=dp), DIMENSION ( :, : ), ALLOCATABLE :: f_numer
  REAL (KIND=dp), DIMENSION ( :, : ), ALLOCATABLE :: rel, diff
  LOGICAL, PARAMETER :: box_change = .FALSE.
  LOGICAL :: ionode

  ionode = globenv % ionode
  group  = globenv % group
  fist_env => force_env % fist_env
! associcating local pointers
  CALL get_fist_env ( fist_env, ewald_pw = ewald_pw, ewald_env = ewald_env, &
                      local_particles = local_particles, particle_set = particle_set, &
                      atomic_kind_set = atomic_kind_set, molecule_set = molecule_set, &
                      local_molecules = local_molecules, &
                      molecule_kind_set = molecule_kind_set, &
                      fist_nonbond_env = fist_nonbond_env, cell = cell )


  CALL fist_nonbond_env_get ( fist_nonbond_env, potparm = potparm )
  CALL ewald_env_get ( ewald_env, ewald_type = ewald_type )
  isos = 0

  natoms = SIZE ( particle_set )

!------------------------------------------------------------------------------
  iw = globenv % scr
  iwf = 7
  ir = 5
  IF ( ionode ) THEN
    WRITE ( iw, '( A,A )' ) ' **************************************', &
                          '*****************************************'
    WRITE ( iw, '( 16x,A )' ) '                   BEGIN FIST DEBUG'
    WRITE ( iw, '( A,A )' ) ' **************************************', &
                          '*****************************************'
  ENDIF
! allocating all arrays involved with debug
  ALLOCATE ( diff ( 3, natoms ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( "debug_control", "diff", 3 * natoms )
  ALLOCATE ( rel ( 3, natoms ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( "debug_control", "rel", 3 * natoms )
  ALLOCATE ( f_numer ( 3, natoms ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( "debug_control", "f_numer", 3 * natoms )
  ALLOCATE ( dbg % f_nonbond ( 3, natoms ), STAT = isos )
  IF ( isos /= 0 ) &
       CALL stop_memory ( "debug_control", "f_nonbond", 3 * natoms )
  ALLOCATE ( dbg % f_g ( 3, natoms ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( "debug_control", "f_g", 3 * natoms )
  ALLOCATE ( dbg % f_bond ( 3, natoms ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( "debug_control", "f_bond", 3 * natoms )
  ALLOCATE ( dbg % f_bend ( 3, natoms ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( "debug_control", "f_bend", 3 * natoms )
  ALLOCATE ( dbg % f_torsion ( 3, natoms ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( "debug_control", "f_torsion", 3 * natoms )
  ALLOCATE ( dbg % f_imptors ( 3, natoms ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( "debug_control", "f_imptors", 3 * natoms )
  ALLOCATE ( dbg % f_ub ( 3, natoms ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( "debug_control", "f_ub", 3 * natoms )
  ALLOCATE ( dbg % f_onef ( 3, natoms ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( "debug_control", "f_onef", 3 * natoms )
  
   dbg % f_nonbond          = 0.0_dp
   dbg % f_g                = 0.0_dp
   dbg % f_bond             = 0.0_dp
   dbg % f_bend             = 0.0_dp
   dbg % f_torsion          = 0.0_dp
   dbg % f_imptors          = 0.0_dp
   dbg % f_onef             = 0.0_dp
   dbg % f_ub               = 0.0_dp
   dbg % pot_nonbond        = 0.0_dp
   dbg % pot_g              = 0.0_dp
   dbg % pot_bond           = 0.0_dp
   dbg % pot_bend           = 0.0_dp
   dbg % pot_torsion        = 0.0_dp
   dbg % pot_imptors        = 0.0_dp
   dbg % pot_onef           = 0.0_dp
   dbg % pot_urey_bradley   = 0.0_dp

! Call to analytical force routine:
  CALL fist_force_control ( fist_env, force_env%virial, globenv, &
                            box_change=.FALSE., debug = dbg, error=error )

  CALL mp_sum ( dbg % f_bond, group )
  CALL mp_sum ( dbg % f_bend, group )
  CALL mp_sum ( dbg % f_torsion, group )
  CALL mp_sum ( dbg % f_imptors, group )
  CALL mp_sum ( dbg % f_onef, group )
  CALL mp_sum ( dbg % f_ub, group )
  
! defining delta
  delta = 1.0E-5_dp
  IF ( ionode ) &
  WRITE ( iw, '( A, T61, G20.12 )' ) ' DELTA = ', delta

! Debug real-space non-bonded interactions

  IF ( ionode ) THEN
    WRITE ( iw, '( A,A )' ) ' **************************************', &
                            '*****************************************'
    WRITE ( iw, '( 16x,A )' ) 'NON-BONDED CONTRIBUTIONS'
    WRITE ( iw, '( A,A )' ) ' **************************************', &
                            '*****************************************'
  END IF
  f_numer = 0.0_dp
  energy_numer = 0.0_dp
  pv_test = 0.0_dp

  CALL force_nonbond_numer ( atomic_kind_set, particle_set, local_particles, &
                             ewald_env, cell, potparm, delta, f_numer, &
                             energy_numer, globenv )
  CALL mp_sum( f_numer, group )
  CALL mp_sum( energy_numer, group )
  IF ( ionode )  THEN
     WRITE ( iw, '( A, T61, G20.12 )' ) &
          ' E ANA NON-BONDED = ', dbg % pot_nonbond
     WRITE ( iw, '( A, T61, G20.12 )' ) &
          ' E DBG NON-BONDED = ', energy_numer

! write out all particle numbers and forces
     WRITE ( iw, '( /,A,A )' ) &
          ' THE FORCES ON ALL PARTICLES ', &
          '(FILE=non_bonded)'
     OPEN (UNIT = iwf, &
     FILE = 'non_bonded', STATUS = 'REPLACE', ACTION = 'WRITE')
     DO iatom = 1, natoms
      WRITE ( iwf, '( A,T71,I10 )' ) ' ATOM NUMBER ', iatom
      WRITE ( iwf, '( A,T21,3G20.12 )' ) &
      ' F ANA NON-BONDED =', dbg % f_nonbond(:,iatom)
      WRITE ( iwf, '( A,T21,3G20.12,/ )' ) &
      ' F NUM NON-BONDED =', f_numer(:,iatom)
     END DO
     CLOSE (iwf)
! compute the absolute value of the differences in the forces
     diff = ABS ( ABS (dbg % f_nonbond) - ABS (f_numer) )
! find the maximum difference and the relative and absolute errors.
     IF ( .NOT. ANY ( dbg % f_nonbond == 0.0_dp ) )  rel = diff / ABS (dbg % f_nonbond)

! find the maximum difference and the relative and absolute errors
     WRITE ( iw, '( A, T61, G20.12 )' ) &
          ' MAXIMUM ABSOLUTE ERROR = ', MAXVAL ( diff )
     IF ( .NOT. ANY ( dbg % f_nonbond == 0.0_dp ) )  THEN
       WRITE ( iw, '( A, T61, G20.12 )' ) &
          ' MAXIMUM RELATIVE ERROR = ', MAXVAL ( rel )
     ENDIF

! write out the particle number and forces of
! the max absolute and relative error
     dum = MAXLOC ( diff )
     WRITE ( iw, '( A, T71, I10 )' ) &
          ' PARTICLE WITH MAX ABSOLUTE ERROR IS ', dum(2)
     WRITE ( iw, '( A, T21, 3G20.12 )' ) &
          ' F ANA NON-BONDED =', dbg % f_nonbond ( :, dum ( 2 ) )
     WRITE ( iw, '( A,T21,3G20.12 )' ) &
          ' F NUM NON-BONDED =', f_numer ( :, dum ( 2 ) )

     IF ( .NOT. ANY ( dbg % f_nonbond == 0.0_dp ) )  THEN
       dum = MAXLOC ( rel )
       WRITE ( iw, '( A,T71,I10 )' ) &
          ' PARTICLE WITH MAX RELATIVE ERROR IS ', dum(2)
       WRITE ( iw, '( A,T21,3G20.12 )' ) &
          ' F ANA NON-BONDED =', dbg % f_nonbond ( :, dum ( 2 ) )
       WRITE ( iw, '( A,T21,3G20.12 )' ) &
          ' F NUM NON-BONDED =', f_numer ( :, dum ( 2 ) )
     ENDIF
  END IF
! Debug non-bonded virial

  CALL pv_nonbond_numer ( fist_nonbond_env, ewald_env, atomic_kind_set, &
                          particle_set, cell, pv_test, delta )

  CALL mp_sum( pv_test, group )
 
  IF ( ionode ) THEN
! write out the virials
    WRITE ( iw, '( A )' ) ' PV ANA NON-BONDED'
    DO i = 1, 3
      WRITE ( iw, '( T21,3G20.12 )' ) dbg % pv_nonbond(i,:)
    END DO
    WRITE ( iw, '( A )' ) ' PV NUM NON-BONDED'
    DO i = 1, 3
     WRITE ( iw, '( T21,3G20.12 )' ) pv_test(i,:)
    END DO
    WRITE ( iw, '( A,A )' ) ' **************************************', &
                            '*****************************************'
    WRITE ( iw, '( A,A )' ) ' **************************************', &
                            '*****************************************'
    WRITE ( iw, '( 16x,A )' ) ' '
    WRITE ( iw, '( A,A )' ) ' **************************************', &
                            '*****************************************'
    WRITE ( iw, '( 16x,A )' ) 'BONDED CONTRIBUTIONS'
    WRITE ( iw, '( A,A )' ) ' **************************************', &
                            '*****************************************'
  ENDIF

! Debug bonds
  f_numer = 0.0_dp
  energy_numer = 0.0_dp
  CALL force_bond_numer ( molecule_set, molecule_kind_set, local_molecules, &
                          particle_set, f_numer, energy_numer, delta )
  CALL mp_sum( f_numer, group )
  CALL mp_sum( energy_numer, group )
  IF ( ionode ) THEN
     WRITE ( iw, '( A, T61, G20.12 )' ) &
       ' E ANA BOND = ', dbg % pot_bond
     WRITE ( iw, '( A, T61, G20.12 )' ) &
       ' E DBG BOND = ', energy_numer

! write out all particle numbers and forces
     WRITE ( iw, '( /,A )' ) &
          ' THE FORCES ON ALL PARTICLES (FILE=bonds)'
     OPEN (UNIT = iwf, &
           FILE = 'bonds', STATUS = 'REPLACE', ACTION = 'WRITE')
     DO iatom = 1, natoms
        WRITE ( iwf, '( A,T71,I10 )' ) ' ATOM NUMBER ', iatom
        WRITE ( iwf, '( A,T21,3G20.12 )' ) &
             ' F ANA BONDS =', dbg % f_bond(:,iatom)
        WRITE ( iwf, '( A,T21,3G20.12,/ )' ) &
             ' F NUM BONDS =', f_numer(:,iatom)
     END DO
     CLOSE (iwf)

! compute the absolute value of the differences in the forces
     diff = ABS ( ABS (dbg % f_bond) - ABS (f_numer) )
! find the maximum difference and the relative and absolute errors.
     IF ( .NOT. ANY ( dbg % f_bond == 0.0_dp ) )  rel = diff / ABS (dbg % f_bond)

! find the maximum difference and the relative and absolute errors
     WRITE ( iw, '( A,T61,G20.12 )' ) &
          ' MAXIMUM ABSOLUTE ERROR = ', MAXVAL ( diff )
     IF ( .NOT. ANY ( dbg % f_bond == 0.0_dp ) ) THEN
       WRITE ( iw, '( A,T61,G20.12 )' ) &
           ' MAXIMUM RELATIVE ERROR = ', MAXVAL ( rel )
     ENDIF

! write out the particle number and forces of
! the max absolute and relative error
     dum = MAXLOC ( diff )
     WRITE ( iw, '( A,T71,I10 )' ) &
          ' PARTICLE WITH MAX ABSOLUTE ERROR IS ', dum(2)
     WRITE ( iw, '( A,T21,3G20.12 )' ) &
          ' F ANA BOND =', dbg % f_bond(:,dum(2))
     WRITE ( iw, '( A,T21,3G20.12 )' ) &
          ' F NUM BOND =', f_numer(:,dum(2))

     IF ( .NOT. ANY ( dbg % f_bond == 0.0_dp ) ) THEN
       dum = MAXLOC ( rel )
       WRITE ( iw, '( A,T71,I10 )' ) &
          ' PARTICLE WITH MAX RELATIVE ERROR IS ', dum(2)
       WRITE ( iw, '( A,T21,3G20.12 )' ) &
          ' F ANA BOND =', dbg % f_bond(:,dum(2))
       WRITE ( iw, '( A,T21,3G20.12 )' ) &
          ' F NUM BOND =', f_numer(:,dum(2))
     END IF
   ENDIF
! Debug bond virial

   pv_test = 0.0_dp
   CALL pv_bond_numer ( molecule_set, molecule_kind_set, local_molecules, &
                        particle_set, cell, pv_test, delta )
   CALL mp_sum ( pv_test, group )
   IF ( ionode ) THEN
! write out the virials
     WRITE ( iw, '( A )' ) ' PV ANA BOND'
     DO i = 1, 3
        WRITE ( iw, '( T21,3G20.12 )' ) dbg % pv_bond(i,:)
     END DO
     WRITE ( iw, '( A )' ) ' PV NUM BOND'
     DO i = 1, 3
        WRITE ( iw, '( T21, 3G20.12 )' ) pv_test ( i, : )
     END DO
     WRITE ( iw, '( A,A )' ) ' **************************************', &
                          '*****************************************'
     WRITE ( iw, '( A,A )' ) ' **************************************', &
                          '*****************************************'
     WRITE ( iw, '( 16x,A )' ) ' '
     WRITE ( iw, '( A,A )' ) ' **************************************', &
                          '*****************************************'
     WRITE ( iw, '( 16x,A )' ) 'UREY-BRADLEY CONTRIBUTIONS'
     WRITE ( iw, '( A,A )' ) ' **************************************', &
                          '*****************************************'
   ENDIF

! Debug urey-bradley
  f_numer = 0.0_dp
  energy_numer = 0.0_dp
  pv_test = 0.0_dp
  CALL force_ub_numer ( molecule_set, molecule_kind_set, local_molecules, &
                        particle_set, f_numer, energy_numer, delta )
  
  CALL mp_sum( f_numer, group )
  CALL mp_sum( energy_numer, group )
  IF ( ionode ) THEN
    WRITE ( iw, '( A, T61, G20.12 )' ) &
       ' E ANA UREY-BRADLEY = ', dbg % pot_urey_bradley
    WRITE ( iw, '( A, T61, G20.12 )' ) &
       ' E DBG UREY-BRADLEY = ', energy_numer

! write out all particle numbers and forces
    WRITE ( iw, '( /,A )' ) &
       'THE FORCES ON ALL PARTICLES (FILE=urey-bradley)'
    OPEN (UNIT = iwf, &
        FILE = 'urey-bradley', STATUS = 'REPLACE', ACTION = 'WRITE')
    DO iatom = 1, natoms
        WRITE ( iwf, '( A,T71,I10 )' ) ' ATOM NUMBER ', iatom
        WRITE ( iwf, '( A,T21,3G20.12 )' ) &
             ' F ANA UREY-BRADLEY =', dbg % f_ub(:,iatom)
        WRITE ( iwf, '( A,T21,3G20.12,/ )' ) &
             ' F NUM UREY-BRADLEY =', f_numer(:,iatom)
    END DO
    CLOSE (iwf)

! compute the absolute value of the differences in the forces
    diff = ABS ( ABS (dbg % f_ub) - ABS (f_numer) )
! find the maximum difference and the relative and absolute errors.
    IF ( .NOT. ANY ( dbg % f_ub == 0.0_dp ) )  rel = diff / ABS (dbg % f_ub)

! find the maximum difference and the relative and absolute errors
    WRITE ( iw, '( A,T61,G20.12 )' ) &
       ' MAXIMUM ABSOLUTE ERROR = ', MAXVAL ( diff )
    IF ( .NOT. ANY ( dbg % f_ub == 0.0_dp ) ) THEN
     WRITE ( iw, '( A,T61,G20.12 )' ) &
    ' MAXIMUM RELATIVE ERROR = ', MAXVAL ( rel )
    ENDIF

! write out the particle number and forces of
! the max absolute and relative error
     dum = MAXLOC ( diff )
    WRITE ( iw, '( A,T71,I10 )' ) &
       ' PARTICLE WITH MAX ABSOLUTE ERROR IS ', dum(2)
    WRITE ( iw, '( A,T21,3G20.12 )' ) &
       ' F ANA UREY-BRADLEY =', dbg % f_ub(:,dum(2))
    WRITE ( iw, '( A,T21,3G20.12 )' ) &
       ' F NUM UREY-BRADLEY =', f_numer(:,dum(2))

    IF ( .NOT. ANY ( dbg % f_ub == 0.0_dp ) ) THEN
     dum = MAXLOC ( rel )
      WRITE ( iw, '( A,T71,I10 )' ) &
         ' PARTICLE WITH MAX RELATIVE ERROR IS ', dum(2)
      WRITE ( iw, '( A,T21,3G20.12 )' ) &
       ' F ANA UREY-BRADLEY =', dbg % f_ub(:,dum(2))
      WRITE ( iw, '( A,T21,3G20.12 )' ) &
         ' F NUM UREY-BRADLEY =', f_numer(:,dum(2))
    END IF
  ENDIF
! Debug ub virial
  CALL pv_ub_numer ( molecule_set, molecule_kind_set, local_molecules, &
                     particle_set, cell, pv_test, delta )

  CALL mp_sum ( pv_test, group )
  IF ( ionode ) THEN
! write out the virials
   WRITE ( iw, '( A )' ) ' PV ANA UREY-BRADLEY'
   DO i = 1, 3
     WRITE ( iw, '( T21,3G20.12 )' ) dbg % pv_ub(i,:)
   END DO
   WRITE ( iw, '( A )' ) ' PV NUM UREY-BRADLEY'
   DO i = 1, 3
     WRITE ( iw, '( T21, 3G20.12 )' ) pv_test ( i, : )
   END DO
   WRITE ( iw, '( A,A )' ) ' **************************************', &
                          '*****************************************'
   WRITE ( iw, '( A,A )' ) ' **************************************', &
                          '*****************************************'
   WRITE ( iw, '( 16x,A )' ) ' '
   WRITE ( iw, '( A,A )' ) ' **************************************', &
                          '*****************************************'
   WRITE ( iw, '( 16x,A )' ) 'ONE-FOUR CONTRIBUTIONS'
   WRITE ( iw, '( A,A )' ) ' **************************************', &
                           '*****************************************'
  ENDIF
! Debug onefours
  f_numer = 0.0_dp
  energy_numer = 0.0_dp
  pv_test = 0.0_dp
  CALL force_onefour_numer ( molecule_set, molecule_kind_set, local_molecules, &
                              particle_set, f_numer, energy_numer, delta )
  CALL mp_sum( f_numer, group )
  CALL mp_sum( energy_numer, group )
  IF ( ionode ) THEN
    WRITE ( iw, '( A, T61, G20.12 )' ) &
       ' E ANA ONE-FOUR = ', dbg % pot_onef
    WRITE ( iw, '( A, T61, G20.12 )' ) &
       ' E DBG ONE-FOUR = ', energy_numer

! write out all particle numbers and forces
    WRITE ( iw, '( /,A )' ) &
      'THE FORCES ON ALL PARTICLES (FILE=onefours)'
    OPEN (UNIT = iwf, &
        FILE = 'onefours', STATUS = 'REPLACE', ACTION = 'WRITE')
    DO iatom = 1, natoms
      WRITE ( iwf, '( A,T71,I10 )' ) ' ATOM NUMBER ', iatom
      WRITE ( iwf, '( A,T21,3G20.12 )' ) &
         ' F ANA ONE-FOUR =', dbg % f_onef(:,iatom)
      WRITE ( iwf, '( A,T21,3G20.12,/ )' ) &
         ' F NUM ONE-FOUR =', f_numer(:,iatom)
    END DO
    CLOSE (iwf)

! compute the absolute value of the differences in the forces
    diff = ABS ( ABS (dbg % f_onef) - ABS (f_numer) )
! find the maximum difference and the relative and absolute errors.
    IF ( .NOT. ANY ( dbg % f_onef == 0.0_dp ) )  rel = diff / ABS (dbg % f_onef)

! find the maximum difference and the relative and absolute errors
    WRITE ( iw, '( A,T61,G20.12 )' ) &
       ' MAXIMUM ABSOLUTE ERROR = ', MAXVAL ( diff )
    IF ( .NOT. ANY ( dbg % f_onef == 0.0_dp ) ) THEN
      WRITE ( iw, '( A,T61,G20.12 )' ) &
        ' MAXIMUM RELATIVE ERROR = ', MAXVAL ( rel )
    ENDIF

! write out the particle number and forces of
! the max absolute and relative error
  dum = MAXLOC ( diff )
    WRITE ( iw, '( A,T71,I10 )' ) &
       ' PARTICLE WITH MAX ABSOLUTE ERROR IS ', dum(2)
    WRITE ( iw, '( A,T21,3G20.12 )' ) &
       ' F ANA ONE-FOUR =', dbg % f_onef(:,dum(2))
    WRITE ( iw, '( A,T21,3G20.12 )' ) &
       ' F NUM ONE-FOUR =', f_numer(:,dum(2))

    IF ( .NOT. ANY ( dbg % f_ub == 0.0_dp ) ) THEN
      dum = MAXLOC ( rel )
      WRITE ( iw, '( A,T71,I10 )' ) &
       ' PARTICLE WITH MAX RELATIVE ERROR IS ', dum(2)
      WRITE ( iw, '( A,T21,3G20.12 )' ) &
       ' F ANA ONE-FOUR =', dbg % f_onef(:,dum(2))
      WRITE ( iw, '( A,T21,3G20.12 )' ) &
       ' F NUM ONE-FOUR =', f_numer(:,dum(2))
    END IF
  ENDIF
! Debug ub virial
  CALL pv_onefour_numer ( molecule_set, molecule_kind_set, local_molecules, &
                          particle_set, cell, pv_test, delta )
  CALL mp_sum ( pv_test, group )
  IF ( ionode ) THEN
! write out the virials
    WRITE ( iw, '( A )' ) ' PV ANA ONE-FOUR'
    DO i = 1, 3
       WRITE ( iw, '( T21,3G20.12 )' ) dbg % pv_onef(i,:)
    END DO
    WRITE ( iw, '( A )' ) ' PV NUM ONE-FOUR'
    DO i = 1, 3
       WRITE ( iw, '( T21, 3G20.12 )' ) pv_test ( i, : )
    END DO
    WRITE ( iw, '( A,A )' ) ' **************************************', &
                            '*****************************************'
    WRITE ( iw, '( A,A )' ) ' **************************************', &
                          '*****************************************'
    WRITE ( iw, '( 16x,A )' ) ' '
    WRITE ( iw, '( A,A )' ) ' **************************************', &
                            '*****************************************'
    WRITE ( iw, '( 16x,A )' ) 'BEND CONTRIBUTIONS'
    WRITE ( iw, '( A,A )' ) ' **************************************', &
                            '*****************************************'
  ENDIF
! Debug bends
  f_numer = 0.0_dp
  energy_numer = 0.0_dp
  pv_test = 0.0_dp
  CALL force_bend_numer ( molecule_set, molecule_kind_set, local_molecules, &
                          particle_set, f_numer, energy_numer, delta )
  CALL mp_sum( f_numer, group )
  CALL mp_sum( energy_numer, group )

  IF ( ionode ) THEN
    WRITE ( iw, '( A, T61, G20.12 )' ) &
       ' E ANA BEND = ', dbg % pot_bend
    WRITE ( iw, '( A, T61, G20.12 )' ) &
       ' E DBG BEND = ', energy_numer

! write out all particle numbers and forces
    WRITE ( iw, '( /,A )' ) &
       ' THE FORCES ON ALL PARTICLES (FILE=bends)'
    OPEN (UNIT = iwf, FILE = 'bends', STATUS = 'REPLACE', ACTION = 'WRITE')
    DO iatom = 1, natoms
       WRITE ( iwf, '( A,T71,I10 )' ) ' ATOM NUMBER ', iatom
       WRITE ( iwf, '( A,T21,3G20.12 )' ) &
          ' F ANA BENDS =', dbg % f_bend(:,iatom)
       WRITE ( iwf, '( A,T21,3G20.12,/ )' ) &
          ' F NUM BENDS =', f_numer(:,iatom)
    END DO
    CLOSE (iwf)

! compute the absolute value of the differences in the forces
    diff = ABS ( ABS (dbg % f_bend) - ABS (f_numer) )
    IF ( .NOT. ANY ( dbg % f_bend == 0.0_dp ) ) rel = diff / ABS (dbg % f_bend)

! find the maximum difference and the relative and absolute errors
    WRITE ( iw, '( A,T61,G20.12 )' ) &
       ' MAXIMUM ABSOLUTE ERROR = ', MAXVAL ( diff )
    IF ( .NOT. ANY ( dbg % f_bend == 0.0_dp ) ) THEN
      WRITE ( iw, '( A,T61,G20.12 )' ) &
        ' MAXIMUM RELATIVE ERROR = ', MAXVAL ( rel )
    ENDIF

! write out the particle number and forces of
! the max absolute and relative error
    dum = MAXLOC ( diff )
    WRITE ( iw, '( A,T71,I10 )' ) &
       ' PARTICLE WITH MAX ABSOLUTE ERROR IS ', dum(2)
    WRITE ( iw, '( A,T21,3G20.12 )' ) &
       ' F ANA BEND =', dbg % f_bend(:,dum(2))
    WRITE ( iw, '( A,T21,3G20.12 )' ) &
             ' F NUM BEND =', f_numer(:,dum(2))

    IF ( .NOT. ANY ( dbg % f_bend == 0.0_dp ) ) THEN
      dum = MAXLOC ( rel )
      WRITE ( iw, '( A,T71,I10 )' ) &
       ' PARTICLE WITH MAX RELATIVE ERROR IS ', dum(2)
      WRITE ( iw, '( A,T21,3G20.12 )' ) &
       ' F ANA BEND =', dbg % f_bend(:,dum(2))
      WRITE ( iw, '( A,T21,3G20.12 )' ) &
       ' F NUM BEND =', f_numer(:,dum(2))
    ENDIF
  END IF
! Debug bend virial
  CALL pv_bend_numer ( molecule_set, molecule_kind_set, local_molecules, &
                       particle_set, cell, pv_test, delta )
  CALL mp_sum ( pv_test, group )
 
  IF ( ionode ) THEN
! write out the virials
    WRITE ( iw, '( A )' ) ' PV ANA BEND'
    DO i = 1, 3
     WRITE ( iw, '( T21, 3G20.12 )' ) dbg % pv_bend(i,:)
    END DO
    WRITE ( iw, '( A )' ) ' PV NUM BEND'
    DO i = 1, 3
       WRITE ( iw, '( T21, 3G20.12 )' ) pv_test(i,:)
    END DO

    WRITE ( iw, '( A,A )' ) ' **************************************', &
                            '*****************************************'
    WRITE ( iw, '( A,A )' ) ' **************************************', &
                            '*****************************************'
    WRITE ( iw, '( 16x,A )' ) ' '
    WRITE ( iw, '( A,A )' ) ' **************************************', &
                            '*****************************************'
    WRITE ( iw, '( 16x,A )' ) 'PROPER TORSION CONTRIBUTIONS'
    WRITE ( iw, '( A,A )' ) ' **************************************', &
                            '*****************************************'
  ENDIF
! Debug torsions
  f_numer = 0.0_dp
  energy_numer = 0.0_dp
  pv_test = 0.0_dp
  CALL force_torsion_numer ( molecule_set, molecule_kind_set, local_molecules, &
                            particle_set, f_numer, energy_numer, delta )
  CALL mp_sum( f_numer, group )
  CALL mp_sum( energy_numer, group )
  
  IF ( ionode ) THEN
    WRITE ( iw, '( A, T61, G20.12 )' ) &
       ' E ANA TORSION = ', dbg % pot_torsion
    WRITE ( iw, '( A, T61, G20.12 )' ) &
       ' E DBG TORSION = ', energy_numer

! write out all particle numbers and forces
    WRITE ( iw, '( /,A,A )' ) &
    ' THE FORCES ON ALL PARTICLES ',&
    '(FILE=torsions)'
    OPEN (UNIT = iwf, &
        FILE = 'torsions', STATUS = 'REPLACE', ACTION = 'WRITE')
    DO iatom = 1, natoms
      WRITE ( iwf, '( A,T71,I10 )' ) ' ATOM NUMBER ', iatom
      WRITE ( iwf, '( A,T21,3G20.12 )' ) ' F ANA TORSION =', &
       dbg % f_torsion(:,iatom)
      WRITE ( iwf, '( A,T21,3G20.12,/ )' ) ' F NUM TORSION =', &
       f_numer(:,iatom)
    END DO
    CLOSE (iwf)

! compute the absolute value of the differences in the forces
    diff = ABS ( ABS (dbg % f_torsion) - ABS (f_numer) )
    IF ( .NOT. ANY ( dbg % f_torsion ( :, : ) == 0.0_dp ) )  &
    rel = diff / ABS (dbg % f_torsion)

! find the maximum difference and the relative and absolute errors
    WRITE ( iw, '( A,T61,G20.12 )' ) &
       ' MAXIMUM ABSOLUTE ERROR = ', MAXVAL ( diff )
    IF ( .NOT. ANY ( dbg % f_torsion == 0.0_dp ) ) THEN
      WRITE ( iw, '( A,T61,G20.12 )' ) &
       ' MAXIMUM RELATIVE ERROR = ', MAXVAL ( rel )
    ENDIF

! write out the particle number and forces of
! the max absolute and relative error
    dum = MAXLOC ( diff )
    WRITE ( iw, '( A,T71,I10 )' ) &
       ' PARTICLE WITH MAX ABSOLUTE ERROR IS ', dum(2)
    WRITE ( iw, '( A,T21,3G20.12 )' ) &
       ' F ANA TORSION =', dbg % f_torsion(:,dum(2))
    WRITE ( iw, '( A,T21,3G20.12 )' ) &
       ' F NUM TORSION =', f_numer(:,dum(2))

    IF ( .NOT. ANY ( dbg % f_torsion == 0.0_dp ) ) THEN
      dum = MAXLOC ( rel )
      WRITE ( iw, '( A,T71,I10 )' ) &
         ' PARTICLE WITH MAX RELATIVE ERROR IS ', dum(2)
      WRITE ( iw, '( A,T21,3G20.12 )' ) &
       ' F ANA TORSION =', dbg % f_torsion(:,dum(2))
      WRITE ( iw, '( A,T21,3G20.12 )' ) &
       ' F NUM TORSION =', f_numer(:,dum(2))
    ENDIF
  END IF
! Debug torsion virial
  CALL pv_torsion_numer ( molecule_set, molecule_kind_set, local_molecules, &
                          particle_set, cell, pv_test, delta )
  CALL mp_sum ( pv_test, group )

  IF ( ionode ) THEN
! write out the virials
    WRITE ( iw, '( A )' ) ' PV ANA TORSION'
    DO i = 1, 3
       WRITE ( iw, '( T21, 3G20.12 )' ) dbg % pv_torsion(i,:)
    END DO
    WRITE ( iw, '( A )' ) ' PV NUM TORSION'
    DO i = 1, 3
       WRITE ( iw, '( T21, 3G20.12 )' ) pv_test(i,:)
    END DO
    WRITE ( iw, '( A,A )' ) ' **************************************', &
                            '*****************************************'
    WRITE ( iw, '( A,A )' ) ' **************************************', &
                            '*****************************************'
    WRITE ( iw, '( 16x,A )' ) ' '
    WRITE ( iw, '( A,A )' ) ' **************************************', &
                          '*****************************************'
    WRITE ( iw, '( 16x,A )' ) 'IMPROPER TORSION CONTRIBUTIONS'
    WRITE ( iw, '( A,A )' ) ' **************************************', &
                            '*****************************************'
  ENDIF

! Debug improper torsions
  f_numer = 0.0_dp
  energy_numer = 0.0_dp
  pv_test = 0.0_dp
  CALL force_imptors_numer ( molecule_set, molecule_kind_set, local_molecules, &
                                  particle_set, f_numer, energy_numer, delta )
  CALL mp_sum( f_numer, group )
  CALL mp_sum( energy_numer, group )
  IF ( ionode ) THEN
    WRITE ( iw, '( A, T61, G20.12 )' ) &
       ' E ANA TORSION = ', dbg % pot_imptors
    WRITE ( iw, '( A, T61, G20.12 )' ) &
       ' E DBG TORSION = ', energy_numer

! write out all particle numbers and forces
    WRITE ( iw, '( /,A,A )' ) &
       ' THE FORCES ON ALL PARTICLES ',&
       '(FILE=imptors)'
    OPEN (UNIT = iwf, &
        FILE = 'imptors', STATUS = 'REPLACE', ACTION = 'WRITE')
    DO iatom = 1, natoms
     WRITE ( iwf, '( A,T71,I10 )' ) ' ATOM NUMBER ', iatom
     WRITE ( iwf, '( A,T21,3G20.12 )' ) ' F ANA IMP TORSION =', &
          dbg % f_imptors(:,iatom)
     WRITE ( iwf, '( A,T21,3G20.12,/ )' ) ' F NUM IMP TORSION =', &
          f_numer(:,iatom)
    END DO
    CLOSE (iwf)

! compute the absolute value of the differences in the forces
    diff = ABS ( ABS (dbg % f_imptors) - ABS (f_numer) )
    IF ( .NOT. ANY ( dbg % f_imptors ( :, : ) == 0.0_dp ) ) &
      rel = diff / ABS (dbg % f_imptors)

! find the maximum difference and the relative and absolute errors
    WRITE ( iw, '( A,T61,G20.12 )' ) &
       ' MAXIMUM ABSOLUTE ERROR = ', MAXVAL ( diff )
    IF ( .NOT. ANY ( dbg % f_imptors == 0.0_dp ) ) THEN
      WRITE ( iw, '( A,T61,G20.12 )' ) &
       ' MAXIMUM RELATIVE ERROR = ', MAXVAL ( rel )
    ENDIF

! write out the particle number and forces of
! the max absolute and relative error
    dum = MAXLOC ( diff )
    WRITE ( iw, '( A,T71,I10 )' ) &
       ' PARTICLE WITH MAX ABSOLUTE ERROR IS ', dum(2)
    WRITE ( iw, '( A,T21,3G20.12 )' ) &
       ' F ANA IMP TORSION =', dbg % f_imptors(:,dum(2))
    WRITE ( iw, '( A,T21,3G20.12 )' ) &
       ' F NUM IMP TORSION =', f_numer(:,dum(2))

    IF ( .NOT. ANY ( dbg % f_imptors == 0.0_dp ) ) THEN
      dum = MAXLOC ( rel )
      WRITE ( iw, '( A,T71,I10 )' ) &
       ' PARTICLE WITH MAX RELATIVE ERROR IS ', dum(2)
          WRITE ( iw, '( A,T21,3G20.12 )' ) &
       ' F ANA IMP TORSION =', dbg % f_imptors(:,dum(2))
      WRITE ( iw, '( A,T21,3G20.12 )' ) &
       ' F NUM IMP TORSION =', f_numer(:,dum(2))
    ENDIF
  END IF
! Debug torsion virial
  CALL pv_imptors_numer ( molecule_set, molecule_kind_set, local_molecules, &
                          particle_set, cell, pv_test, delta )
  CALL mp_sum ( pv_test, group )

  IF ( ionode ) THEN
! write out the virials
    WRITE ( iw, '( A )' ) ' PV ANA IMP TORSION'
    DO i = 1, 3
     WRITE ( iw, '( T21, 3G20.12 )' ) dbg % pv_imptors (i,:)
    END DO
    WRITE ( iw, '( A )' ) ' PV NUM IMP TORSION'
    DO i = 1, 3
     WRITE ( iw, '( T21, 3G20.12 )' ) pv_test(i,:)
    END DO
    WRITE ( iw, '( A,A )' ) ' **************************************', &
                          '*****************************************'
    WRITE ( iw, '( A,A )' ) ' **************************************', &
                          '*****************************************'
    WRITE ( iw, '( 16x,A )' ) ' '

  END IF
! Debug g-space
  CALL ewald_env_get ( ewald_env, ewald_type = ewald_type )
  IF ( ionode ) THEN
    IF ( ewald_type /= do_ewald_none ) THEN
      WRITE ( iw, '( A,A )' ) ' **************************************', &
                              '*****************************************'
      WRITE ( iw, '( 16x,A )' ) 'G-SPACE CONTRIBUTIONS'
      WRITE ( iw, '( A,A )' ) ' **************************************', &
                              '*****************************************'
    ENDIF
  ENDIF
  SELECT CASE ( ewald_type )
  CASE ( do_ewald_ewald )
    energy_numer = 0.0_dp
    CALL ewald_energy_numer ( ewald_env, ewald_pw, atomic_kind_set, &
                              particle_set, local_particles, energy_numer )
    IF ( ionode ) &
    WRITE ( iw, '( A, T61, G20.12 )' ) ' G-SPACE DBG ENERGY = ', &
         energy_numer
  CASE ( do_ewald_pme )
    IF ( ionode ) &
    WRITE ( iw, '( A, T61, G20.12 )' ) ' G-SPACE ANA ENERGY = ', &
         dbg % pot_g
  CASE ( do_ewald_spme )
    CALL stop_program ( "fist_debug", " NO SPME energy debug routines " )
  END SELECT

! Now forces
  SELECT CASE ( ewald_type )
  CASE ( do_ewald_ewald )
    f_numer = 0.0_dp
    CALL ewald_force_numer ( ewald_env, ewald_pw, atomic_kind_set,  &
                             particle_set, local_particles, delta, f_numer )
    CALL mp_sum ( f_numer, group )    
    CALL mp_sum ( dbg % f_g, group )    
  CASE ( do_ewald_pme )
    f_numer = 0.0_dp
    CALL pme_force_numer ( ewald_env, ewald_pw, particle_set, delta, f_numer )
  CASE ( do_ewald_spme )
    CALL stop_program ( "fist_debug", " NO SPME debug routines " )
  END SELECT
! write out all particle numbers and forces
  IF ( ionode ) THEN
   WRITE ( iw, '( /,A )' ) &
   ' THE FORCES ON ALL PARTICLES (FILE=gspace)'
   OPEN (UNIT = iwf, FILE = 'gspace', STATUS = 'REPLACE', ACTION = 'WRITE')
   DO iatom = 1, natoms
     WRITE ( iwf, '( A,T71,I10 )' ) ' ATOM NUMBER ', iatom
     WRITE ( iwf, '( A,T21,3G20.12 )' ) &
   ' F ANA GSPACE =', dbg % f_g(:,iatom)
     WRITE ( iwf, '( A,T21,3G20.12,/ )' ) &
   ' F NUM GSPACE =', f_numer(:,iatom)
   END DO
   CLOSE (iwf)

! compute the absolute value of the differences in the forces
   diff = ABS ( ABS (dbg % f_g) - ABS (f_numer) )
   IF ( .NOT. ANY ( dbg % f_g == 0.0_dp ) ) rel = diff / ABS (dbg % f_g)

! find the maximum difference and the relative and absolute errors
   WRITE ( iw, '( A,T61,G20.12 )' ) 'MAXIMUM ABSOLUTE ERROR = ', &
       MAXVAL (diff)
   IF ( .NOT. ANY ( dbg % f_g == 0.0_dp ) ) THEN 
     WRITE ( iw, '( A,T61,G20.12 )' ) 'MAXIMUM RELATIVE ERROR = ', &
       MAXVAL (rel)
   ENDIF

! write out the particle number and forces of the max absolute
! and relative error
   dum = MAXLOC ( diff )
   WRITE ( iw, '( A,T71,I10 )' ) &
        ' PARTICLE WITH MAX ABSOLUTE ERROR IS ', dum(2)
   WRITE ( iw, '( A,T21,3G20.12 )' ) ' F ANA G-SPACE =', &
        dbg % f_g(:,dum(2))
   WRITE ( iw, '( A,T21,3G20.12 )' ) ' F NUM G-SPACE =', &
        f_numer(:,dum(2))

   IF ( .NOT. ANY ( dbg % f_g == 0.0_dp ) ) THEN 
     dum = MAXLOC ( rel )
     WRITE ( iw, '( A,T71,I10 )' ) &
        'PARTICLE WITH MAX RELATIVE ERROR IS ', dum(2)
     WRITE ( iw, '( A,T21,3G20.12 )' ) ' F ANA G-SPACE =', &
        dbg % f_g(:,dum(2))
     WRITE ( iw, '( A,T21,3G20.12 )' ) ' F NUM G-SPACE =', &
        f_numer(:,dum(2))
   END IF
  ENDIF

! Debug g-space virial
  SELECT CASE ( ewald_type )
  CASE ( do_ewald_ewald )
    pv_test = 0.0_dp
    CALL pv_g_numer ( ewald_env, ewald_pw, atomic_kind_set, &
                       particle_set, local_particles, delta, pv_test )
  CASE ( do_ewald_pme )
    CALL stop_program ( "fist_debug", " NO PME debug routines " )
  CASE ( do_ewald_spme )
    CALL stop_program ( "fist_debug", " NO SPME debug routines " )
  END SELECT

! write out the virials
  IF ( ionode ) THEN
   WRITE ( iw, '( A )' ) ' PV NUM G-SPACE'
   DO i = 1, 3
    WRITE ( iw, '( T21,3G20.12 )' ) pv_test(i,:)
   END DO
   WRITE ( iw, '( A )' ) ' PV ANA G-SPACE'
   DO i = 1, 3
     WRITE ( iw, '( T21,3G20.12 )' ) dbg % pv_g(i,:)
   END DO
  ENDIF
! deallocating all arrays involved with debug
  DEALLOCATE ( diff, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( "debug_control", "diff", 3 * natoms )
  DEALLOCATE ( rel, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( "debug_control", "rel", 3 * natoms )
  DEALLOCATE ( f_numer, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( "debug_control", "f_numer", 3 * natoms )
  DEALLOCATE ( dbg % f_nonbond, STAT = isos )
  IF ( isos /= 0 ) &
       CALL stop_memory ( "debug_control", "f_nonbond", 3 * natoms )
  DEALLOCATE ( dbg % f_g, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( "debug_control", "f_g", 3 * natoms )
  DEALLOCATE ( dbg % f_bond, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( "debug_control", "f_bond", 3 * natoms )
  DEALLOCATE ( dbg % f_bend, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( "debug_control", "f_bend", 3 * natoms )
  DEALLOCATE ( dbg % f_torsion, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( "debug_control", "f_torsion", 3 * natoms )
  DEALLOCATE ( dbg % f_imptors, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( "debug_control", "f_imptors", 3 * natoms )
  DEALLOCATE ( dbg % f_ub, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( "debug_control", "f_ub", 3 * natoms )
  DEALLOCATE ( dbg % f_onef, STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( "debug_control", "f_onef", 3 * natoms )
  IF ( ionode ) THEN
    WRITE ( iw, '( )' )
    WRITE ( iw, '( A,A )' ) ' **************************************', &
                            '*****************************************'
    WRITE ( iw, '( 16x, A )' ) &
         '*                     END FIST DEBUG                   *'
    WRITE ( iw, '( A,A )' ) ' **************************************', &
                            '*****************************************'
  ENDIF
END SUBROUTINE debug_control

!!*****
!******************************************************************************

END MODULE fist_debug

!******************************************************************************
