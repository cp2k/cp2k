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

  USE atomic_kind_types,                        ONLY : atomic_kind_type
  USE distribution_1d_types,                    ONLY : distribution_1d_type
  USE ewald_environment_types,                  ONLY : ewald_environment_type, &
                                                       ewald_env_get
  USE ewald_pw_types,                           ONLY : ewald_pw_type
  USE fist_environment_types,                   ONLY: fist_environment_type, &
                                                      get_fist_env
  USE fist_force,                               ONLY : fist_force_control,         &
                                                       debug_variables_type
  USE fist_force_numer,                          ONLY : force_bond_numer,      &
                                                        force_bend_numer,      &
                                                        force_torsion_numer,   &
                                                        force_imptors_numer,   &
                                                        force_nonbond_numer,   &
                                                        force_onefour_numer,   &
                                                        force_ub_numer,   &
                                                        ewald_energy_numer,    &
                                                        ewald_force_numer,     &
                                                        pv_bond_numer,         &
                                                        pv_bend_numer,         &
                                                        pv_torsion_numer,      &
                                                        pv_imptors_numer,      &
                                                        pv_nonbond_numer,      &
                                                        pv_onefour_numer,      &
                                                        pv_ub_numer,      &
                                                        pv_g_numer,            &
                                                        pme_force_numer
  USE fist_nonbond_env_types,          ONLY : fist_nonbond_env_type, fist_nonbond_env_get
  USE force_env_types,                          ONLY: force_env_type
  USE global_types,                             ONLY : global_environment_type
  USE kinds,                                    ONLY : dbl, default_string_length
  USE molecule_kind_types,                      ONLY:  molecule_kind_type
  USE molecule_types_new,                       ONLY: molecule_type
  USE pair_potential_types,                     ONLY: pair_potential_type
  USE particle_types,                           ONLY : particle_type
  USE simulation_cell,                          ONLY : cell_type
  USE splines,                                  ONLY: spline_environment_type
  USE termination,                              ONLY : stop_memory, stop_program
  USE virial_types,                             ONLY : virial_type

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

SUBROUTINE  debug_control ( force_env, globenv )

  IMPLICIT NONE
! Arguments
  TYPE ( force_env_type ), POINTER :: force_env
  TYPE ( global_environment_type ), INTENT ( IN ) :: globenv

! Locals
  INTEGER :: iflag, iatom, iw, iwf, ir, nbond, nbend, ntorsion
  INTEGER, DIMENSION ( 2 ) :: dum
  INTEGER :: i, natoms, nnodes, isos, n, nlocal_particles
  CHARACTER ( len = default_string_length ) :: ewald_type
  TYPE (pair_potential_type), POINTER, DIMENSION ( :, : ) :: potparm
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
  TYPE ( spline_environment_type ), POINTER :: spline_env
  REAL ( dbl ) :: delta, energy_numer, cutoff, ecut
  REAL ( dbl ) :: e_numer, pv_test ( 3, 3 )
  REAL ( dbl ) :: err1, numer, denom1, vec ( 3 ), e_bc, e_real, energy_tot
  REAL ( dbl ) :: denom2, err2
  REAL ( dbl ), DIMENSION ( :, : ), ALLOCATABLE :: f_numer
  REAL ( dbl ), DIMENSION ( :, : ), ALLOCATABLE :: rel, diff
  LOGICAL, PARAMETER :: box_change = .FALSE.

  fist_env => force_env % fist_env
! associcating local pointers
  CALL get_fist_env ( fist_env, ewald_pw = ewald_pw, ewald_env = ewald_env, &
                      local_particles = local_particles, particle_set = particle_set, &
                      atomic_kind_set = atomic_kind_set, molecule_set = molecule_set, &
                      local_molecules = local_molecules, &
                      molecule_kind_set = molecule_kind_set, &
                      fist_nonbond_env = fist_nonbond_env, cell = cell )


  CALL fist_nonbond_env_get ( fist_nonbond_env, spline_env = spline_env )
  potparm => spline_env % potparm
  CALL ewald_env_get ( ewald_env, ewald_type = ewald_type )
  isos = 0

  natoms = SIZE ( particle_set )

!------------------------------------------------------------------------------
  iw = globenv % scr
  iwf = 7
  ir = 5
  WRITE ( iw, '( 16x,A )' ) &
       '***************************************************'
  WRITE ( iw, '( 16x,A )' ) '                   BEGIN FIST DEBUG'
  WRITE ( iw, '( 16x,A )' ) &
       '***************************************************'

! allocating all arrays involved with debug
  IF ( .NOT. ALLOCATED ( diff ) ) ALLOCATE ( diff ( 3, natoms ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( "debug_control", "diff", 3 * natoms )
  IF ( .NOT. ALLOCATED ( rel ) ) ALLOCATE ( rel ( 3, natoms ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( "debug_control", "rel", 3 * natoms )
  IF ( .NOT. ALLOCATED ( f_numer ) ) &
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
  
   dbg % f_onef = 0._dbl
   dbg % f_ub = 0._dbl
   dbg % f_imptors = 0._dbl
   dbg % f_torsion = 0._dbl
   dbg % f_bond= 0._dbl
   dbg % f_bend= 0._dbl
   dbg % f_g= 0._dbl

! Call to analytical force routine:
  CALL fist_force_control ( fist_env, force_env%virial, globenv, &
                            box_change=.FALSE., debug = dbg )

! Debug real-space non-bonded interactions
  WRITE ( iw, '( A )' ) &
       ' DO YOU WANT TO DEBUG YOUR REAL-SPACE NON-BONDED INTERACTIONS (1=yes)?'
  READ ( ir, * ) iflag
  IF ( iflag == 1 ) THEN

! get numerical force for the non-bonded interactions
     WRITE ( iw, '( A )' ) ' ENTER A DELTA LESS THAN 1'
     READ ( ir, * ) delta
     IF ( delta >= 1.0_dbl .OR. delta <= 0.0_dbl ) THEN
        delta = 1.0E-5_dbl
        WRITE ( iw, '( A, T61, G20.12 )' ) &
             ' DELTA (changed to default) = ', delta
     ELSE
        WRITE ( iw, '( A, T61, G20.12 )' ) ' DELTA = ', delta
     END IF
     CALL force_nonbond_numer ( atomic_kind_set, particle_set, local_particles, &
                                ewald_env, cell, potparm, delta, f_numer, energy_numer, & 
                                globenv )
     WRITE ( iw, '( A, T61, G20.12 )' ) &
          ' E ANA NON-BONDED = ', dbg % pot_nonbond
     WRITE ( iw, '( A, T61, G20.12 )' ) &
          ' E DBG NON-BONDED = ', energy_numer

! write out all particle numbers and forces
     WRITE ( iw, '( /,A,A )' ) &
          ' DO YOU WANT TO GET THE FORCES ON ALL PARTICLES ', &
          '(FILE=non_bonded, 1=yes)?'
     READ ( ir, * ) iflag
     IF ( iflag == 1 ) THEN
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
     END IF

! compute the absolute value of the differences in the forces
     diff = ABS ( ABS (dbg % f_nonbond) - ABS (f_numer) )
! find the maximum difference and the relative and absolute errors.
    IF ( .NOT. ANY ( dbg % f_nonbond == 0._dbl ) )  rel = diff / ABS (dbg % f_nonbond)

! find the maximum difference and the relative and absolute errors
     WRITE ( iw, '( A, T61, G20.12 )' ) &
          ' MAXIMUM ABSOLUTE ERROR = ', MAXVAL ( diff )
     IF ( .NOT. ANY ( dbg % f_nonbond == 0._dbl ) )  THEN
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

     IF ( .NOT. ANY ( dbg % f_nonbond == 0._dbl ) )  THEN
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
  WRITE ( iw, '( /,A )' ) &
       ' DO YOU WANT TO DEBUG YOUR NON-BONDED VIRIAL (1=yes)?'
  READ ( ir, * ) iflag
  IF ( iflag == 1 ) THEN

! get numerical virial
     WRITE ( iw, '( A )' ) ' ENTER A DELTA LESS THAN 1'
     READ ( ir, * ) delta
     IF ( delta >= 1.0_dbl .OR. delta <= 0.0_dbl ) THEN
        delta = 1.0E-5_dbl
        WRITE ( iw, '( A, T61, G20.12 )' ) &
             ' DELTA (changed to default) = ', delta
     ELSE
        WRITE ( iw, '( A, T61, G20.12 )' ) ' DELTA = ', delta
     END IF
     CALL pv_nonbond_numer ( fist_nonbond_env, ewald_env, atomic_kind_set, &
                             particle_set, cell, pv_test, delta )

! write out the virials
     WRITE ( iw, '( A )' ) ' PV ANA NON-BONDED'
     DO i = 1, 3
        WRITE ( iw, '( T21,3G20.12 )' ) dbg % pv_nonbond(i,:)
     END DO
     WRITE ( iw, '( A )' ) ' PV NUM NON-BONDED'
     DO i = 1, 3
        WRITE ( iw, '( T21,3G20.12 )' ) pv_test(i,:)
     END DO
  END IF

! Debug bonds
     WRITE ( iw, '( /,A )' ) &
          ' DO YOU WANT TO DEBUG YOUR BONDS (1=yes)?'
     READ ( ir, * ) iflag
     IF ( iflag == 1 ) THEN
        WRITE ( iw, '( A )' ) ' ENTER A DELTA LESS THAN 1'
        READ ( ir, * ) delta
        IF ( delta >= 1.0_dbl .OR. delta <= 0.0_dbl ) THEN
           delta = 1.0E-5_dbl
           WRITE ( iw, '( A, T61, G20.12 )' ) &
                ' DELTA (changed to default) = ', delta
        ELSE
           WRITE ( iw, '( A, T61, G20.12 )' ) ' DELTA = ', delta
        END IF
        CALL force_bond_numer ( molecule_set, molecule_kind_set, local_molecules, &
                                particle_set, f_numer, energy_numer, delta )
        WRITE ( iw, '( A, T61, G20.12 )' ) &
             ' E ANA BOND = ', dbg % pot_bond
        WRITE ( iw, '( A, T61, G20.12 )' ) &
             ' E DBG BOND = ', energy_numer

! write out all particle numbers and forces
     WRITE ( iw, '( /,A )' ) &
          ' DO YOU WANT TO GET THE FORCES ON ALL PARTICLES (FILE=bonds, 1=yes)?'
     READ ( ir, * ) iflag
     IF ( iflag == 1 ) THEN
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
     END IF

! compute the absolute value of the differences in the forces
        diff = ABS ( ABS (dbg % f_bond) - ABS (f_numer) )
! find the maximum difference and the relative and absolute errors.
        IF ( .NOT. ANY ( dbg % f_bond == 0._dbl ) )  rel = diff / ABS (dbg % f_bond)

! find the maximum difference and the relative and absolute errors
        WRITE ( iw, '( A,T61,G20.12 )' ) &
             ' MAXIMUM ABSOLUTE ERROR = ', MAXVAL ( diff )
        IF ( .NOT. ANY ( dbg % f_bond == 0._dbl ) ) THEN
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

        IF ( .NOT. ANY ( dbg % f_bond == 0._dbl ) ) THEN
          dum = MAXLOC ( rel )
          WRITE ( iw, '( A,T71,I10 )' ) &
             ' PARTICLE WITH MAX RELATIVE ERROR IS ', dum(2)
          WRITE ( iw, '( A,T21,3G20.12 )' ) &
             ' F ANA BOND =', dbg % f_bond(:,dum(2))
          WRITE ( iw, '( A,T21,3G20.12 )' ) &
             ' F NUM BOND =', f_numer(:,dum(2))
        END IF
     END IF

! Debug bond virial
     WRITE ( iw, '( /,A )' ) &
          ' DO YOU WANT TO DEBUG YOUR BOND VIRIAL (1=yes)?'
     READ ( ir, * ) iflag
     IF ( iflag == 1 ) THEN

! get numerical virial
        WRITE ( iw, '( A )' ) ' ENTER A DELTA LESS THAN 1'
        READ ( ir, * ) delta
        IF ( delta >= 1 .OR. delta <= 0 ) THEN
           delta = 1.0E-5_dbl
           WRITE ( iw, '( A, T61, G20.12 )' ) &
                ' DELTA (changed to default) = ', delta
        ELSE
           WRITE ( iw, '( A, T61, G20.12 )' ) ' DELTA = ', delta
        END IF
        CALL pv_bond_numer ( molecule_set, molecule_kind_set, local_molecules, &
                             particle_set, cell, pv_test, delta )

! write out the virials
        WRITE ( iw, '( A )' ) ' PV ANA BOND'
        DO i = 1, 3
           WRITE ( iw, '( T21,3G20.12 )' ) dbg % pv_bond(i,:)
        END DO
        WRITE ( iw, '( A )' ) ' PV NUM BOND'
        DO i = 1, 3
           WRITE ( iw, '( T21, 3G20.12 )' ) pv_test ( i, : )
        END DO
     END IF

! Debug urey-bradley
     WRITE ( iw, '( /,A )' ) &
          ' DO YOU WANT TO DEBUG YOUR UREY-BRADLEYS (1=yes)?'
     READ ( ir, * ) iflag
     IF ( iflag == 1 ) THEN
        WRITE ( iw, '( A )' ) ' ENTER A DELTA LESS THAN 1'
        READ ( ir, * ) delta
        IF ( delta >= 1.0_dbl .OR. delta <= 0.0_dbl ) THEN
           delta = 1.0E-5_dbl
           WRITE ( iw, '( A, T61, G20.12 )' ) &
                ' DELTA (changed to default) = ', delta
        ELSE
           WRITE ( iw, '( A, T61, G20.12 )' ) ' DELTA = ', delta
        END IF
        CALL force_ub_numer ( molecule_set, molecule_kind_set, local_molecules, &
                                particle_set, f_numer, energy_numer, delta )
        WRITE ( iw, '( A, T61, G20.12 )' ) &
             ' E ANA UREY-BRADLEY = ', dbg % pot_urey_bradley
        WRITE ( iw, '( A, T61, G20.12 )' ) &
             ' E DBG UREY-BRADLEY = ', energy_numer

! write out all particle numbers and forces
     WRITE ( iw, '( /,A )' ) &
          ' DO YOU WANT TO GET THE FORCES ON ALL PARTICLES (FILE=urey-bradley, 1=yes)?'
     READ ( ir, * ) iflag
     IF ( iflag == 1 ) THEN
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
     END IF

! compute the absolute value of the differences in the forces
        diff = ABS ( ABS (dbg % f_ub) - ABS (f_numer) )
! find the maximum difference and the relative and absolute errors.
        IF ( .NOT. ANY ( dbg % f_ub == 0._dbl ) )  rel = diff / ABS (dbg % f_ub)

! find the maximum difference and the relative and absolute errors
        WRITE ( iw, '( A,T61,G20.12 )' ) &
             ' MAXIMUM ABSOLUTE ERROR = ', MAXVAL ( diff )
        IF ( .NOT. ANY ( dbg % f_ub == 0._dbl ) ) THEN
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

        IF ( .NOT. ANY ( dbg % f_ub == 0._dbl ) ) THEN
          dum = MAXLOC ( rel )
          WRITE ( iw, '( A,T71,I10 )' ) &
             ' PARTICLE WITH MAX RELATIVE ERROR IS ', dum(2)
          WRITE ( iw, '( A,T21,3G20.12 )' ) &
             ' F ANA UREY-BRADLEY =', dbg % f_ub(:,dum(2))
          WRITE ( iw, '( A,T21,3G20.12 )' ) &
             ' F NUM UREY-BRADLEY =', f_numer(:,dum(2))
        END IF
     END IF

! Debug ub virial
     WRITE ( iw, '( /,A )' ) &
          ' DO YOU WANT TO DEBUG YOUR UREY-BRADLEY VIRIAL (1=yes)?'
     READ ( ir, * ) iflag
     IF ( iflag == 1 ) THEN

! get numerical virial
        WRITE ( iw, '( A )' ) ' ENTER A DELTA LESS THAN 1'
        READ ( ir, * ) delta
        IF ( delta >= 1 .OR. delta <= 0 ) THEN
           delta = 1.0E-5_dbl
           WRITE ( iw, '( A, T61, G20.12 )' ) &
                ' DELTA (changed to default) = ', delta
        ELSE
           WRITE ( iw, '( A, T61, G20.12 )' ) ' DELTA = ', delta
        END IF
        CALL pv_ub_numer ( molecule_set, molecule_kind_set, local_molecules, &
                             particle_set, cell, pv_test, delta )

! write out the virials
        WRITE ( iw, '( A )' ) ' PV ANA UREY-BRADLEY'
        DO i = 1, 3
           WRITE ( iw, '( T21,3G20.12 )' ) dbg % pv_ub(i,:)
        END DO
        WRITE ( iw, '( A )' ) ' PV NUM UREY-BRADLEY'
        DO i = 1, 3
           WRITE ( iw, '( T21, 3G20.12 )' ) pv_test ( i, : )
        END DO
     END IF

! Debug onefours
     WRITE ( iw, '( /,A )' ) &
          ' DO YOU WANT TO DEBUG YOUR ONE-FOURS (1=yes)?'
     READ ( ir, * ) iflag
     IF ( iflag == 1 ) THEN
        WRITE ( iw, '( A )' ) ' ENTER A DELTA LESS THAN 1'
        READ ( ir, * ) delta
        IF ( delta >= 1.0_dbl .OR. delta <= 0.0_dbl ) THEN
           delta = 1.0E-5_dbl
           WRITE ( iw, '( A, T61, G20.12 )' ) &
                ' DELTA (changed to default) = ', delta
        ELSE
           WRITE ( iw, '( A, T61, G20.12 )' ) ' DELTA = ', delta
        END IF
        CALL force_onefour_numer ( molecule_set, molecule_kind_set, local_molecules, &
                              particle_set, f_numer, energy_numer, delta )
        WRITE ( iw, '( A, T61, G20.12 )' ) &
             ' E ANA ONE-FOUR = ', dbg % pot_onef
        WRITE ( iw, '( A, T61, G20.12 )' ) &
             ' E DBG ONE-FOUR = ', energy_numer

! write out all particle numbers and forces
     WRITE ( iw, '( /,A )' ) &
          ' DO YOU WANT TO GET THE FORCES ON ALL PARTICLES (FILE=onefours, 1=yes)?'
     READ ( ir, * ) iflag
     IF ( iflag == 1 ) THEN
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
     END IF

! compute the absolute value of the differences in the forces
        diff = ABS ( ABS (dbg % f_onef) - ABS (f_numer) )
! find the maximum difference and the relative and absolute errors.
        IF ( .NOT. ANY ( dbg % f_onef == 0._dbl ) )  rel = diff / ABS (dbg % f_onef)

! find the maximum difference and the relative and absolute errors
        WRITE ( iw, '( A,T61,G20.12 )' ) &
             ' MAXIMUM ABSOLUTE ERROR = ', MAXVAL ( diff )
        IF ( .NOT. ANY ( dbg % f_onef == 0._dbl ) ) THEN
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

        IF ( .NOT. ANY ( dbg % f_ub == 0._dbl ) ) THEN
          dum = MAXLOC ( rel )
          WRITE ( iw, '( A,T71,I10 )' ) &
             ' PARTICLE WITH MAX RELATIVE ERROR IS ', dum(2)
          WRITE ( iw, '( A,T21,3G20.12 )' ) &
             ' F ANA ONE-FOUR =', dbg % f_onef(:,dum(2))
          WRITE ( iw, '( A,T21,3G20.12 )' ) &
             ' F NUM ONE-FOUR =', f_numer(:,dum(2))
        END IF
     END IF

! Debug ub virial
     WRITE ( iw, '( /,A )' ) &
          ' DO YOU WANT TO DEBUG YOUR ONE-FOUR VIRIAL (1=yes)?'
     READ ( ir, * ) iflag
     IF ( iflag == 1 ) THEN

! get numerical virial
        WRITE ( iw, '( A )' ) ' ENTER A DELTA LESS THAN 1'
        READ ( ir, * ) delta
        IF ( delta >= 1 .OR. delta <= 0 ) THEN
           delta = 1.0E-5_dbl
           WRITE ( iw, '( A, T61, G20.12 )' ) &
                ' DELTA (changed to default) = ', delta
        ELSE
           WRITE ( iw, '( A, T61, G20.12 )' ) ' DELTA = ', delta
        END IF
        CALL pv_onefour_numer ( molecule_set, molecule_kind_set, local_molecules, &
                                particle_set, cell, pv_test, delta )

! write out the virials
        WRITE ( iw, '( A )' ) ' PV ANA ONE-FOUR'
        DO i = 1, 3
           WRITE ( iw, '( T21,3G20.12 )' ) dbg % pv_onef(i,:)
        END DO
        WRITE ( iw, '( A )' ) ' PV NUM ONE-FOUR'
        DO i = 1, 3
           WRITE ( iw, '( T21, 3G20.12 )' ) pv_test ( i, : )
        END DO
     END IF

! Debug bends
     WRITE ( iw, '( /,A )' ) &
          ' DO YOU WANT TO DEBUG YOUR BENDS (1=yes)?'
     READ ( ir, * ) iflag
     IF ( iflag == 1 ) THEN
        WRITE ( iw, '( A )' ) ' ENTER A DELTA LESS THAN 1'
        READ ( ir, * ) delta
        IF ( delta >= 1.0_dbl .OR. delta <= 0.0_dbl ) THEN
           delta = 1.0E-5_dbl
           WRITE ( iw, '( A, T61, G20.12 )' ) &
                ' DELTA (changed to default) = ', delta
        ELSE
           WRITE ( iw, '( A, T61, G20.12 )' ) ' DELTA = ', delta
        END IF
        CALL force_bend_numer ( molecule_set, molecule_kind_set, local_molecules, &
                                particle_set, f_numer, energy_numer, delta )
        WRITE ( iw, '( A, T61, G20.12 )' ) &
             ' E ANA BEND = ', dbg % pot_bend
        WRITE ( iw, '( A, T61, G20.12 )' ) &
             ' E DBG BEND = ', energy_numer

! write out all particle numbers and forces
     WRITE ( iw, '( /,A )' ) &
          ' DO YOU WANT TO GET THE FORCES ON ALL PARTICLES (FILE=bends, 1=yes)?'
     READ ( ir, * ) iflag
     IF ( iflag == 1 ) THEN
        OPEN (UNIT = iwf, FILE = 'bends', STATUS = 'REPLACE', ACTION = 'WRITE')
        DO iatom = 1, natoms
        WRITE ( iwf, '( A,T71,I10 )' ) ' ATOM NUMBER ', iatom
        WRITE ( iwf, '( A,T21,3G20.12 )' ) &
             ' F ANA BENDS =', dbg % f_bend(:,iatom)
        WRITE ( iwf, '( A,T21,3G20.12,/ )' ) &
             ' F NUM BENDS =', f_numer(:,iatom)
        END DO
        CLOSE (iwf)
     END IF

! compute the absolute value of the differences in the forces
        diff = ABS ( ABS (dbg % f_bend) - ABS (f_numer) )
        IF ( .NOT. ANY ( dbg % f_bend == 0._dbl ) ) rel = diff / ABS (dbg % f_bend)

! find the maximum difference and the relative and absolute errors
        WRITE ( iw, '( A,T61,G20.12 )' ) &
             ' MAXIMUM ABSOLUTE ERROR = ', MAXVAL ( diff )
        IF ( .NOT. ANY ( dbg % f_bend == 0._dbl ) ) THEN
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

        IF ( .NOT. ANY ( dbg % f_bend == 0._dbl ) ) THEN
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
     WRITE ( iw, '( /,A )' ) &
          ' DO YOU WANT TO DEBUG YOUR BEND VIRIAL (1=yes)?'
     READ ( ir, * ) iflag
     IF ( iflag == 1 ) THEN

! get numerical virial
        WRITE ( iw, '( A )' ) ' ENTER A DELTA LESS THAN 1'
        READ ( ir, * ) delta
        IF ( delta >= 1.0_dbl .OR. delta <= 0.0_dbl ) THEN
           delta = 1.0E-5_dbl
           WRITE ( iw, '( A, T61, G20.12 )' ) &
                ' DELTA (changed to default) = ', delta
        ELSE
           WRITE ( iw, '( A, T61, G20.12 )' ) ' DELTA = ', delta
        END IF
        CALL pv_bend_numer ( molecule_set, molecule_kind_set, local_molecules, &
                             particle_set, cell, pv_test, delta )

! write out the virials
        WRITE ( iw, '( A )' ) ' PV ANA BEND'
        DO i = 1, 3
           WRITE ( iw, '( T21, 3G20.12 )' ) dbg % pv_bend(i,:)
        END DO
        WRITE ( iw, '( A )' ) ' PV NUM BEND'
        DO i = 1, 3
           WRITE ( iw, '( T21, 3G20.12 )' ) pv_test(i,:)
        END DO
     END IF

! Debug torsions
     WRITE ( iw, '( /,A )' ) &
          ' DO YOU WANT TO DEBUG YOUR TORSIONS (1=yes)?'
     READ ( ir, * ) iflag
     IF ( iflag == 1 ) THEN
        WRITE ( iw, '( A )' ) ' ENTER A DELTA LESS THAN 1'
        READ ( ir, * ) delta
        IF ( delta >= 1.0_dbl .OR. delta <= 0.0_dbl ) THEN
           delta = 1.0E-5_dbl
           WRITE ( iw, '( A, T61, G20.12 )' ) &
                ' DELTA (changed to default) = ', delta
        ELSE
           WRITE ( iw, '( A, T61, G20.12 )' ) ' DELTA = ', delta
        END IF
        CALL force_torsion_numer ( molecule_set, molecule_kind_set, local_molecules, &
                                  particle_set, f_numer, energy_numer, delta )
        WRITE ( iw, '( A, T61, G20.12 )' ) &
             ' E ANA TORSION = ', dbg % pot_torsion
        WRITE ( iw, '( A, T61, G20.12 )' ) &
             ' E DBG TORSION = ', energy_numer

! write out all particle numbers and forces
     WRITE ( iw, '( /,A,A )' ) &
          ' DO YOU WANT TO GET THE FORCES ON ALL PARTICLES ',&
          '(FILE=torsions, 1=yes)?'
     READ ( ir, * ) iflag
     IF ( iflag == 1 ) THEN
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
     END IF

! compute the absolute value of the differences in the forces
        diff = ABS ( ABS (dbg % f_torsion) - ABS (f_numer) )
        IF ( .NOT. ANY ( dbg % f_torsion ( :, : ) == 0._dbl ) ) rel = diff / ABS (dbg % f_torsion)

! find the maximum difference and the relative and absolute errors
        WRITE ( iw, '( A,T61,G20.12 )' ) &
             ' MAXIMUM ABSOLUTE ERROR = ', MAXVAL ( diff )
        IF ( .NOT. ANY ( dbg % f_torsion == 0._dbl ) ) THEN
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

        IF ( .NOT. ANY ( dbg % f_torsion == 0._dbl ) ) THEN
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
     WRITE ( iw, '( /,A )' ) &
          ' DO YOU WANT TO DEBUG YOUR TORSION VIRIAL (1=yes)?'
     READ ( ir, * ) iflag
     IF ( iflag == 1 ) THEN

! get numerical virial
        WRITE ( iw, '( A )' ) ' ENTER A DELTA LESS THAN 1'
        READ ( ir, * ) delta
        IF ( delta >= 1.0_dbl .OR. delta <= 0.0_dbl ) THEN
           delta = 1.0E-5_dbl
           WRITE ( iw, '( A, T61, G20.12 )' ) &
                ' DELTA (changed to default) = ', delta
        ELSE
           WRITE ( iw, '( A, T61, G20.12 )' ) ' DELTA = ', delta
        END IF
        CALL pv_torsion_numer ( molecule_set, molecule_kind_set, local_molecules, &
                                particle_set, cell, pv_test, delta )

! write out the virials
        WRITE ( iw, '( A )' ) ' PV ANA TORSION'
        DO i = 1, 3
           WRITE ( iw, '( T21, 3G20.12 )' ) dbg % pv_torsion(i,:)
        END DO
        WRITE ( iw, '( A )' ) ' PV NUM TORSION'
        DO i = 1, 3
           WRITE ( iw, '( T21, 3G20.12 )' ) pv_test(i,:)
        END DO
     END IF

! Debug improper torsions
     WRITE ( iw, '( /,A )' ) &
          ' DO YOU WANT TO DEBUG YOUR IMPROPER TORSIONS (1=yes)?'
     READ ( ir, * ) iflag
     IF ( iflag == 1 ) THEN
        WRITE ( iw, '( A )' ) ' ENTER A DELTA LESS THAN 1'
        READ ( ir, * ) delta
        IF ( delta >= 1.0_dbl .OR. delta <= 0.0_dbl ) THEN
           delta = 1.0E-5_dbl
           WRITE ( iw, '( A, T61, G20.12 )' ) &
                ' DELTA (changed to default) = ', delta
        ELSE
           WRITE ( iw, '( A, T61, G20.12 )' ) ' DELTA = ', delta
        END IF
        CALL force_imptors_numer ( molecule_set, molecule_kind_set, local_molecules, &
                                  particle_set, f_numer, energy_numer, delta )
        WRITE ( iw, '( A, T61, G20.12 )' ) &
             ' E ANA TORSION = ', dbg % pot_imptors
        WRITE ( iw, '( A, T61, G20.12 )' ) &
             ' E DBG TORSION = ', energy_numer

! write out all particle numbers and forces
     WRITE ( iw, '( /,A,A )' ) &
          ' DO YOU WANT TO GET THE FORCES ON ALL PARTICLES ',&
          '(FILE=imptors, 1=yes)?'
     READ ( ir, * ) iflag
     IF ( iflag == 1 ) THEN
        OPEN (UNIT = iwf, &
              FILE = 'imptors', STATUS = 'REPLACE', ACTION = 'WRITE')
        DO iatom = 1, natoms
        WRITE ( iwf, '( A,T71,I10 )' ) ' ATOM NUMBER ', iatom
        WRITE ( iwf, '( A,T21,3G20.12 )' ) ' F ANA IMP TORSION =', &
             dbg % f_torsion(:,iatom)
        WRITE ( iwf, '( A,T21,3G20.12,/ )' ) ' F NUM IMP TORSION =', &
             f_numer(:,iatom)
        END DO
        CLOSE (iwf)
     END IF

! compute the absolute value of the differences in the forces
        diff = ABS ( ABS (dbg % f_imptors) - ABS (f_numer) )
        IF ( .NOT. ANY ( dbg % f_imptors ( :, : ) == 0._dbl ) ) &
            rel = diff / ABS (dbg % f_imptors)

! find the maximum difference and the relative and absolute errors
        WRITE ( iw, '( A,T61,G20.12 )' ) &
             ' MAXIMUM ABSOLUTE ERROR = ', MAXVAL ( diff )
        IF ( .NOT. ANY ( dbg % f_imptors == 0._dbl ) ) THEN
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

        IF ( .NOT. ANY ( dbg % f_imptors == 0._dbl ) ) THEN
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
     WRITE ( iw, '( /,A )' ) &
          ' DO YOU WANT TO DEBUG YOUR IMP TORSION VIRIAL (1=yes)?'
     READ ( ir, * ) iflag
     IF ( iflag == 1 ) THEN

! get numerical virial
        WRITE ( iw, '( A )' ) ' ENTER A DELTA LESS THAN 1'
        READ ( ir, * ) delta
        IF ( delta >= 1.0_dbl .OR. delta <= 0.0_dbl ) THEN
           delta = 1.0E-5_dbl
           WRITE ( iw, '( A, T61, G20.12 )' ) &
                ' DELTA (changed to default) = ', delta
        ELSE
           WRITE ( iw, '( A, T61, G20.12 )' ) ' DELTA = ', delta
        END IF
        CALL pv_imptors_numer ( molecule_set, molecule_kind_set, local_molecules, &
                                particle_set, cell, pv_test, delta )

! write out the virials
        WRITE ( iw, '( A )' ) ' PV ANA IMP TORSION'
        DO i = 1, 3
           WRITE ( iw, '( T21, 3G20.12 )' ) dbg % pv_imptors (i,:)
        END DO
        WRITE ( iw, '( A )' ) ' PV NUM IMP TORSION'
        DO i = 1, 3
           WRITE ( iw, '( T21, 3G20.12 )' ) pv_test(i,:)
        END DO
     END IF

! Debug g-space
 CALL ewald_env_get ( ewald_env, ewald_type = ewald_type )
 IF ( ewald_type /= 'NONE' ) THEN
    WRITE ( iw, '( A )' ) &
         ' DO YOU WANT TO DEBUG YOUR G-SPACE INTERACTIONS (1=yes)?'
    READ ( ir, * ) iflag
    IF ( iflag == 1 ) THEN
       WRITE ( iw, '( A )' ) &
            ' DO YOU WANT TO DEBUG YOUR G-SPACE ENERGIES (1=yes)?'
       READ ( ir, * ) iflag
       IF ( iflag == 1 ) THEN
          SELECT CASE ( ewald_type )
          CASE ( "EWALD" )
            CALL ewald_energy_numer ( ewald_env, ewald_pw, atomic_kind_set, &
                                      particle_set, local_particles, energy_numer )
          CASE ( "PME" )
            CALL stop_program ( "fist_debug", " NO PME energy debug routines " )
          CASE ( "SPME" )
            CALL stop_program ( "fist_debug", " NO SPME energy debug routines " )
          END SELECT
          WRITE ( iw, '( A, T61, G20.12 )' ) ' G-SPACE ANA ENERGY = ', &
               dbg % pot_g
          WRITE ( iw, '( A, T61, G20.12 )' ) ' G-SPACE DBG ENERGY = ', &
               energy_numer
       END IF
       WRITE ( iw, '( A )' ) &
            ' DO YOU WANT TO DEBUG YOUR G-SPACE FORCES (1=yes)?'
       READ ( ir, * ) iflag
       IF ( iflag == 1 ) THEN
          WRITE ( iw, '( A )' ) ' ENTER A DELTA LESS THAN 1'
          READ ( ir, * ) delta
          IF ( delta >= 1.0_dbl .OR. delta <= 0.0_dbl ) THEN
             delta = 1.0E-5_dbl
             WRITE ( iw, '( A, T61, G20.12 )' ) &
                  ' DELTA (changed to default) = ', delta
          ELSE
             WRITE ( iw, '( A, T61, G20.12 )' ) ' DELTA = ', delta
          END IF
          SELECT CASE ( ewald_type )
          CASE ( "EWALD" )
            CALL ewald_force_numer ( ewald_env, ewald_pw, atomic_kind_set,  &
                                      particle_set, local_particles, delta, f_numer )
          CASE ( "PME" )
            CALL pme_force_numer ( ewald_env, ewald_pw, particle_set, delta, f_numer )
          CASE ( "SPME" )
            CALL stop_program ( "fist_debug", " NO SPME debug routines " )
          END SELECT
! write out all particle numbers and forces
          WRITE ( iw, '( /,A )' ) &
          ' DO YOU WANT TO GET THE FORCES ON ALL PARTICLES (FILE=gspace, 1=yes)?'
          READ ( ir, * ) iflag
          IF ( iflag == 1 ) THEN
             OPEN (UNIT = iwf, FILE = 'gspace', STATUS = 'REPLACE', ACTION = 'WRITE')
             DO iatom = 1, natoms
             WRITE ( iwf, '( A,T71,I10 )' ) ' ATOM NUMBER ', iatom
             WRITE ( iwf, '( A,T21,3G20.12 )' ) &
                  ' F ANA GSPACE =', dbg % f_g(:,iatom)
             WRITE ( iwf, '( A,T21,3G20.12,/ )' ) &
                  ' F NUM GSPACE =', f_numer(:,iatom)
             END DO
             CLOSE (iwf)
          END IF

! compute the absolute value of the differences in the forces

! compute the absolute value of the differences in the forces
          diff = ABS ( ABS (dbg % f_g) - ABS (f_numer) )
          IF ( .NOT. ANY ( dbg % f_g == 0._dbl ) ) rel = diff / ABS (dbg % f_g)

! find the maximum difference and the relative and absolute errors
          WRITE ( iw, '( A,T61,G20.12 )' ) 'MAXIMUM ABSOLUTE ERROR = ', &
               MAXVAL (diff)
          IF ( .NOT. ANY ( dbg % f_g == 0._dbl ) ) THEN 
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

          IF ( .NOT. ANY ( dbg % f_g == 0._dbl ) ) THEN 
            dum = MAXLOC ( rel )
            WRITE ( iw, '( A,T71,I10 )' ) &
               'PARTICLE WITH MAX RELATIVE ERROR IS ', dum(2)
            WRITE ( iw, '( A,T21,3G20.12 )' ) ' F ANA G-SPACE =', &
               dbg % f_g(:,dum(2))
            WRITE ( iw, '( A,T21,3G20.12 )' ) ' F NUM G-SPACE =', &
               f_numer(:,dum(2))
          END IF
       END IF

! Debug g-space virial
       WRITE ( iw, '( A )' ) &
            ' DO YOU WANT TO DEBUG YOUR G-SPACE VIRIAL (1=yes)?'
       READ ( ir, * ) iflag
       IF ( iflag == 1 ) THEN

! get numerical virial
          WRITE ( iw, '( A )' ) ' ENTER A DELTA LESS THAN 1'
          READ ( ir, * ) delta
          IF ( delta >= 1.0_dbl .OR. delta <= 0.0_dbl ) THEN
             delta = 1.0E-5_dbl
             WRITE ( iw, '( A,T71,F10.6 )' ) &
                  ' DELTA (changed to default) = ', delta
          ELSE
             WRITE ( iw, '( A,T71,F10.6 )' ) ' DELTA = ', delta
          END IF
          SELECT CASE ( ewald_type )
          CASE ( "EWALD" )
            CALL pv_g_numer ( ewald_env, ewald_pw, atomic_kind_set, &
                               particle_set, local_particles, delta, pv_test )
          CASE ( "PME" )
            CALL stop_program ( "fist_debug", " NO PME debug routines " )
          CASE ( "SPME" )
            CALL stop_program ( "fist_debug", " NO SPME debug routines " )
          END SELECT

! write out the virials
          WRITE ( iw, '( A )' ) ' PV NUM G-SPACE'
          DO i = 1, 3
             WRITE ( iw, '( T21,3G20.12 )' ) pv_test(i,:)
          END DO
          WRITE ( iw, '( A )' ) ' PV ANA G-SPACE'
          DO i = 1, 3
             WRITE ( iw, '( T21,3G20.12 )' ) dbg % pv_g(i,:)
          END DO
       END IF
    END IF
 END IF

  WRITE ( iw, '( )' )
  WRITE ( iw, '( 16x, A )' ) &
       '***************************************************'
  WRITE ( iw, '( 16x, A )' ) &
       '*                     END FIST DEBUG                   *'
  WRITE ( iw, '( 16x, A )' ) &
       '***************************************************'
END SUBROUTINE debug_control

!!*****
!******************************************************************************

END MODULE fist_debug

!******************************************************************************
