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
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE fist_debug
  
  USE ewald_parameters_types, ONLY : ewald_parameters_type
  USE fist_force_numer, ONLY : force_bond_numer, force_bend_numer, &
       force_torsion_numer, force_nonbond_numer, force_recip_numer, &
       pv_bond_numer, pv_bend_numer, pv_torsion_numer, ptens_numer, &
       pvg_numer, potential_g_numer, de_g_numer, energy_recip_numer
  USE fist_force, ONLY : force_control, debug_variables_type
  USE global_types, ONLY : global_environment_type
  USE kinds, ONLY : dbl
  USE linklists, ONLY : bonds, bends, torsions, dist_constraints, &
       g3x3_constraints
  USE md, ONLY : thermodynamic_type
  USE molecule_types, ONLY : molecule_structure_type, particle_node_type
  USE particle_types, ONLY : particle_type
  USE pair_potential, ONLY : potentialparm_type, ener_coul
  USE pw_grid_types, ONLY : pw_grid_type, HALFSPACE
  USE pw_grids, ONLY : pw_find_cutoff, pw_grid_setup
  USE simulation_cell, ONLY : cell_type
  USE stop_program, ONLY : stop_memory, stop_prg
  USE linklist_types, ONLY : linklist_internal_data_type
  
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

SUBROUTINE debug_control ( globenv, ewald_param, part, pnode, molecule, box, &
     thermo, potparm, ensemble, lldata )
  
  IMPLICIT NONE
  
! Arguments
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv
  TYPE ( ewald_parameters_type ), INTENT ( INOUT ) :: ewald_param
  TYPE ( particle_type ), INTENT ( INOUT ), DIMENSION ( : ) :: part
  TYPE ( particle_node_type ), INTENT ( INOUT ), DIMENSION ( : ) :: pnode
  TYPE ( potentialparm_type ), INTENT ( IN ), DIMENSION ( :, : ) :: potparm
  TYPE ( molecule_structure_type ), INTENT ( IN ), DIMENSION ( : ) :: molecule
  TYPE ( cell_type ), INTENT ( INOUT ) :: box
  TYPE ( thermodynamic_type ), INTENT ( INOUT ) :: thermo
  CHARACTER ( LEN = * ), INTENT ( IN ) :: ensemble
  TYPE ( linklist_internal_data_type ) :: lldata
  
! Locals
  TYPE ( debug_variables_type ) :: dbg
  TYPE ( pw_grid_type ) :: ewald_grid
  INTEGER :: iflag, i, natoms, isos, iatom, iw, iwf, ir, npts_s(3), gmax
  INTEGER, DIMENSION ( 2 ) :: dum
  REAL ( dbl ) :: delta, energy_numer, cutoff
  REAL ( dbl ) :: e_numer, pv_test ( 3, 3 )
  REAL ( dbl ) :: err1, numer, denom1, vec ( 3 ), e_bc, e_real, energy_tot
  REAL ( dbl ) :: denom2, err2
  REAL ( dbl ), DIMENSION ( :, : ), ALLOCATABLE :: f_bc, f_real
  REAL ( dbl ), DIMENSION ( :, : ), ALLOCATABLE :: f_numer
  REAL ( dbl ), DIMENSION ( :, : ), ALLOCATABLE :: rel, diff
  
!------------------------------------------------------------------------------
  
  iw = globenv % scr
  iwf = 7
  ir = 5
  WRITE ( iw, '( 16x,A )' ) &
       '***************************************************'
  WRITE ( iw, '( 16x,A )' ) '                   BEGIN DEBUG'
  WRITE ( iw, '( 16x,A )' ) &
       '***************************************************'
  
! allocating all arrays involved with debug
  natoms = SIZE ( pnode )
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
  
! Call to analytical force routine: 
  CALL force_control ( molecule, pnode, part, box, thermo, & 
       potparm, ewald_param, ensemble, globenv, lldata, dbg )
  
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
     CALL force_nonbond_numer ( ewald_param, pnode, box, potparm,  &
          delta, f_numer, energy_numer )
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
     rel = diff / ABS (dbg % f_nonbond)
     
! find the maximum difference and the relative and absolute errors
     WRITE ( iw, '( A, T61, G20.12 )' ) &
          ' MAXIMUM ABSOLUTE ERROR = ', MAXVAL ( diff )
     WRITE ( iw, '( A, T61, G20.12 )' ) &
          ' MAXIMUM RELATIVE ERROR = ', MAXVAL ( rel )
     
! write out the particle number and forces of 
! the max absolute and relative error         
     dum = MAXLOC ( diff )
     WRITE ( iw, '( A, T71, I10 )' ) &
          ' PARTICLE WITH MAX ABSOLUTE ERROR IS ', dum(2)
     WRITE ( iw, '( A, T21, 3G20.12 )' ) &
          ' F ANA NON-BONDED =', dbg % f_nonbond ( :, dum ( 2 ) )
     WRITE ( iw, '( A,T21,3G20.12 )' ) &
          ' F NUM NON-BONDED =', f_numer ( :, dum ( 2 ) )
     
     dum = MAXLOC ( rel )
     WRITE ( iw, '( A,T71,I10 )' ) &
          ' PARTICLE WITH MAX RELATIVE ERROR IS ', dum(2)
     WRITE ( iw, '( A,T21,3G20.12 )' ) &
          ' F ANA NON-BONDED =', dbg % f_nonbond ( :, dum ( 2 ) )
     WRITE ( iw, '( A,T21,3G20.12 )' ) &
          ' F NUM NON-BONDED =', f_numer ( :, dum ( 2 ) )
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
     CALL ptens_numer ( ewald_param, pnode, box, potparm, pv_test, delta )

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
  IF ( bonds ) THEN
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
        CALL force_bond_numer ( molecule, f_numer, energy_numer, delta )
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
        rel = diff / ABS (dbg % f_bond)
        
! find the maximum difference and the relative and absolute errors
        WRITE ( iw, '( A,T61,G20.12 )' ) &
             ' MAXIMUM ABSOLUTE ERROR = ', MAXVAL ( diff )
        WRITE ( iw, '( A,T61,G20.12 )' ) &
             ' MAXIMUM RELATIVE ERROR = ', MAXVAL ( rel )
        
! write out the particle number and forces of 
! the max absolute and relative error         
        dum = MAXLOC ( diff )
        WRITE ( iw, '( A,T71,I10 )' ) &
             ' PARTICLE WITH MAX ABSOLUTE ERROR IS ', dum(2)
        WRITE ( iw, '( A,T21,3G20.12 )' ) &
             ' F ANA BOND =', dbg % f_bond(:,dum(2))
        WRITE ( iw, '( A,T21,3G20.12 )' ) &
             ' F NUM BOND =', f_numer(:,dum(2))
        
        dum = MAXLOC ( rel )
        WRITE ( iw, '( A,T71,I10 )' ) &
             ' PARTICLE WITH MAX RELATIVE ERROR IS ', dum(2)
        WRITE ( iw, '( A,T21,3G20.12 )' ) &
             ' F ANA BOND =', dbg % f_bond(:,dum(2))
        WRITE ( iw, '( A,T21,3G20.12 )' ) &
             ' F NUM BOND =', f_numer(:,dum(2))
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
        CALL pv_bond_numer ( molecule, pv_test, box, pnode, delta )
        
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
  END IF
  
! Debug bends
  IF ( bends ) THEN
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
        CALL force_bend_numer ( molecule, f_numer, energy_numer, delta )
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
        rel = diff / ABS (dbg % f_bend)
        
! find the maximum difference and the relative and absolute errors
        WRITE ( iw, '( A,T61,G20.12 )' ) &
             ' MAXIMUM ABSOLUTE ERROR = ', MAXVAL ( diff )
        WRITE ( iw, '( A,T61,G20.12 )' ) &
             ' MAXIMUM RELATIVE ERROR = ', MAXVAL ( rel )
        
! write out the particle number and forces of 
! the max absolute and relative error         
        dum = MAXLOC ( diff )
        WRITE ( iw, '( A,T71,I10 )' ) &
             ' PARTICLE WITH MAX ABSOLUTE ERROR IS ', dum(2)
        WRITE ( iw, '( A,T21,3G20.12 )' ) &
             ' F ANA BEND =', dbg % f_bend(:,dum(2))
        WRITE ( iw, '( A,T21,3G20.12 )' ) &
             ' F NUM BEND =', f_numer(:,dum(2))
        
        dum = MAXLOC ( rel )
        WRITE ( iw, '( A,T71,I10 )' ) &
             ' PARTICLE WITH MAX RELATIVE ERROR IS ', dum(2)
        WRITE ( iw, '( A,T21,3G20.12 )' ) &
             ' F ANA BEND =', dbg % f_bend(:,dum(2))
        WRITE ( iw, '( A,T21,3G20.12 )' ) &
             ' F NUM BEND =', f_numer(:,dum(2))
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
        CALL pv_bend_numer ( molecule, pv_test, box, pnode, delta )
        
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
  END IF
  
! Debug torsions
  IF ( torsions ) THEN
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
        CALL force_torsion_numer ( molecule, f_numer, energy_numer, delta )
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
        rel = diff / ABS (dbg % f_torsion)
        
! find the maximum difference and the relative and absolute errors
        WRITE ( iw, '( A,T61,G20.12 )' ) &
             ' MAXIMUM ABSOLUTE ERROR = ', MAXVAL ( diff )
        WRITE ( iw, '( A,T61,G20.12 )' ) &
             ' MAXIMUM RELATIVE ERROR = ', MAXVAL ( rel )
        
! write out the particle number and forces of 
! the max absolute and relative error         
        dum = MAXLOC ( diff )
        WRITE ( iw, '( A,T71,I10 )' ) &
             ' PARTICLE WITH MAX ABSOLUTE ERROR IS ', dum(2)
        WRITE ( iw, '( A,T21,3G20.12 )' ) &
             ' F ANA TORSION =', dbg % f_torsion(:,dum(2))
        WRITE ( iw, '( A,T21,3G20.12 )' ) &
             ' F NUM TORSION =', f_numer(:,dum(2))
        
        dum = MAXLOC ( rel )
        WRITE ( iw, '( A,T71,I10 )' ) &
             ' PARTICLE WITH MAX RELATIVE ERROR IS ', dum(2)
        WRITE ( iw, '( A,T21,3G20.12 )' ) &
             ' F ANA TORSION =', dbg % f_torsion(:,dum(2))
        WRITE ( iw, '( A,T21,3G20.12 )' ) &
             ' F NUM TORSION =', f_numer(:,dum(2))
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
        CALL pv_torsion_numer ( molecule, pv_test, box, pnode, delta )
        
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
  END IF
  
! Debug g-space 
 IF ( ewald_param % ewald_type /= 'NONE' ) THEN
    WRITE ( iw, '( A )' ) &
         ' DO YOU WANT TO DEBUG YOUR G-SPACE INTERACTIONS (1=yes)?'
    READ ( ir, * ) iflag
    IF ( iflag == 1 ) THEN
       gmax = ewald_param % gmax
       ewald_grid % bounds ( 1, : ) = -gmax / 2
       ewald_grid % bounds ( 2, : ) = +gmax / 2
       npts_s = (/ gmax, gmax, gmax /)
       ewald_grid % grid_span = HALFSPACE
        
       CALL pw_find_cutoff ( npts_s, box, cutoff )
        
       CALL pw_grid_setup( box, ewald_grid, cutoff)

       WRITE ( iw, '( A )' ) &
            ' DO YOU WANT TO DEBUG YOUR G-SPACE ENERGIES (1=yes)?'
       READ ( ir, * ) iflag
       IF ( iflag == 1 ) THEN
          WRITE ( iw, '( A, T71, I10 )' ) ' TOTAL NUMBER OF G-VECTORS = ', &
               ewald_grid % ngpts_cut
          WRITE ( iw, '( A, T71, F10.4 )' ) ' ALPHA= ', &
               ewald_param % alpha
          CALL energy_recip_numer ( ewald_param, pnode, box,  &
	          ewald_grid, energy_numer, ewald_grid % ngpts_cut )
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
          
          CALL force_recip_numer ( ewald_param, pnode, box, ewald_grid, &
	                           delta, f_numer )
          
! compute the absolute value of the differences in the forces
          diff = ABS ( ABS (dbg % f_g) - ABS (f_numer) )
          rel = diff / ABS (dbg % f_g)
          
! find the maximum difference and the relative and absolute errors
          WRITE ( iw, '( A,T61,G20.12 )' ) 'MAXIMUM ABSOLUTE ERROR = ', &
               MAXVAL (diff)
          WRITE ( iw, '( A,T61,G20.12 )' ) 'MAXIMUM RELATIVE ERROR = ', &                       
               MAXVAL (rel)
         
! write out the particle number and forces of the max absolute  
! and relative error         
          dum = MAXLOC ( diff )
          WRITE ( iw, '( A,T71,I10 )' ) &
               ' PARTICLE WITH MAX ABSOLUTE ERROR IS ', dum(2)
          WRITE ( iw, '( A,T21,3G20.12 )' ) ' F ANA G-SPACE =', &
               dbg % f_g(:,dum(2))
          WRITE ( iw, '( A,T21,3G20.12 )' ) ' F NUM G-SPACE =', &
               f_numer(:,dum(2))
          
          dum = MAXLOC ( rel )
          WRITE ( iw, '( A,T71,I10 )' ) &
               'PARTICLE WITH MAX RELATIVE ERROR IS ', dum(2)
          WRITE ( iw, '( A,T21,3G20.12 )' ) ' F ANA G-SPACE =', &
               dbg % f_g(:,dum(2))
          WRITE ( iw, '( A,T21,3G20.12 )' ) ' F NUM G-SPACE =', &
               f_numer(:,dum(2))
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
          CALL pvg_numer ( ewald_param, pnode, box, &
	                   ewald_grid, pv_test, delta )
          
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
       '*                     END DEBUG                   *'
  WRITE ( iw, '( 16x, A )' ) &
       '***************************************************'
END SUBROUTINE debug_control

!!*****
!******************************************************************************

END MODULE fist_debug
