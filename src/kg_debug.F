!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/kg_debug [1.0] *
!!
!!   NAME
!!     kg_debug
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     GT   02-Oct-2002  
!!
!!   MODIFICATION HISTORY
!!
!!   SOURCE
!******************************************************************************

MODULE kg_debug

  USE cube_utils, ONLY : cube_info_type
  USE gaussian_gridlevels, ONLY: gridlevel_info_type
  USE kg_environment_types, ONLY : kg_environment_type
  USE kg_force_numer, ONLY : fpart_gspace_numer,&
                             fpart_xc_numer, fpart_core_rspace_numer
  USE kg_force, ONLY : force_control, debug_variables_type, get_grids
  USE global_types, ONLY : global_environment_type
  USE kinds, ONLY : dbl
  USE l_utils, ONLY : l_info_type
  USE linklists, ONLY : bonds, bends, torsions, dist_constraints, &
       g3x3_constraints
  USE molecule_types, ONLY : molecule_structure_type, particle_node_type
  USE particle_types, ONLY : particle_type
  USE pw_grid_types, ONLY : pw_grid_type
  USE simulation_cell, ONLY : cell_type
  USE termination, ONLY : stop_memory, stop_program

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: debug_control

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** kg_debug/debug_control [1.0] *
!!
!!   NAME
!!     debug_control
!!
!!   FUNCTION
!!     Routine to execute numerical tests for debugging
!!
!!   AUTHOR
!!     gt 02-Oct-2002
!!
!!   MODIFICATION HISTORY
!!   SOURCE
!******************************************************************************

SUBROUTINE  debug_control ( globenv, kg_env )

  IMPLICIT NONE

! Arguments
  TYPE ( kg_environment_type ), INTENT ( INOUT ) :: kg_env
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv

! Locals

  TYPE ( cell_type ), POINTER :: box
  TYPE ( cell_type ), POINTER  :: box_ref
  TYPE ( cube_info_type ), POINTER :: cube_info ( : ) 
  TYPE ( debug_variables_type ) :: dbg
  TYPE ( gridlevel_info_type ), POINTER :: gridlevel_info
  TYPE ( l_info_type ), POINTER    :: l_info
  TYPE ( molecule_structure_type ), POINTER :: molecule ( : )
  TYPE ( particle_type ), POINTER  :: part ( : )
  TYPE ( particle_node_type ), POINTER  :: pnode ( : )
  TYPE ( pw_grid_type ), POINTER, DIMENSION ( : ) :: pw_grid
  INTEGER :: iflag, i, natoms, isos, iatom, iw, iwf, ir, ref_grid, gmax
  INTEGER, DIMENSION ( 2 ) :: dum
  REAL ( dbl ) :: delta, energy_numer
  REAL ( dbl ) :: err1, numer, denom1, vec ( 3 ), e_bc, e_real, energy_tot
  REAL ( dbl ) :: denom2, err2
  REAL ( dbl ), DIMENSION ( :, : ), ALLOCATABLE :: f_bc, f_real
  REAL ( dbl ), DIMENSION ( :, : ), ALLOCATABLE :: f_numer
  REAL ( dbl ), DIMENSION ( :, : ), ALLOCATABLE :: rel, diff
  LOGICAL, PARAMETER :: box_change = .FALSE.

!------------------------------------------------------------------------------

! Assigning local pointers
  pnode => kg_env % pnode
  part => kg_env % part
  molecule => kg_env % molecule
  box_ref => kg_env % box_ref
  box => kg_env % box

  iw = globenv % scr
  iwf = 7
  ir = 5
  WRITE ( iw, '( 16x,A )' ) &
       '***************************************************'
  WRITE ( iw, '( 16x,A )' ) '                   BEGIN KIM-GORDON DEBUG'
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
  ALLOCATE ( dbg % f_hartree ( 3, natoms ), STAT = isos )
  IF ( isos /= 0 ) &
       CALL stop_memory ( "debug_control", "f_hartree", 3 * natoms )
  ALLOCATE ( dbg % f_bc ( 3, natoms ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( "debug_control", "f_bc", 3 * natoms )
  ALLOCATE ( dbg % f_xc ( 3, natoms ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( "debug_control", "f_xc", 3 * natoms )
  ALLOCATE ( dbg % f_rspace_core ( 3, natoms ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( "debug_control", "f_rspace_core", 3 * natoms )

! Call to analytical force routine:
  CALL force_control (  kg_env, globenv, dbg )

  WRITE ( iw, '( A )' ) &
       ' DO YOU WANT TO DEBUG YOUR G-SPACE INTERACTIONS (1=yes)?'
  READ ( ir, * ) iflag
  IF ( iflag == 1 ) THEN
! get numerical force for the total gspace contribution on the particles
    WRITE ( iw, '(A)' ) ' ENTER A DELTA LESS THAN 1'
    READ ( ir, * ) delta
    IF ( delta >= 1 .OR. delta <= 0 ) THEN
      delta = 1.E-5
      WRITE ( iw, '(A,T71,F10.6)' ) &
        ' DELTA ( changed to default ) = ', delta
    ELSE
      WRITE ( iw, '(A,T71,F10.6)' ) ' DELTA = ', delta
    END IF

    CALL get_grids ( pw_grids = pw_grid, &              
                      gridlevel_info = gridlevel_info, &
                      l_info = l_info, cube_info = cube_info, &
                      ref_grid = ref_grid ) 

    CALL fpart_gspace_numer ( delta, kg_env, pw_grid, ref_grid, gridlevel_info, &
                              l_info, cube_info, f_numer, energy_numer, globenv )

    WRITE ( iw, '( A, T61, G20.12 )' ) ' G-SPACE ANA ENERGY = ', &
               dbg % pot_hartree
    WRITE ( iw, '( A, T61, G20.12 )' ) ' G-SPACE DBG ENERGY = ', &
               energy_numer

! write out all particle numbers and forces
    WRITE ( iw, '( /,A )' ) &
     ' DO YOU WANT TO GET THE FORCES ON ALL PARTICLES (FILE=gspace, 1=yes)?'
    READ ( ir, * ) iflag
    IF ( iflag == 1 ) THEN
       OPEN (UNIT = iwf, FILE = 'gspace', STATUS = 'REPLACE', ACTION = 'WRITE')
       DO iatom = 1, natoms
         WRITE ( iwf, '( A,T71,I10 )' ) ' ATOM NUMBER ', iatom
         WRITE ( iwf, '( A,T21,3G20.12 )' ) &
            ' F ANA GSPACE =', dbg % f_hartree(:,iatom)
         WRITE ( iwf, '( A,T21,3G20.12,/ )' ) &
         ' F NUM GSPACE =', f_numer(:,iatom)
       END DO
       CLOSE (iwf)
    END IF

! compute the absolute value of the differences in the forces

! compute the absolute value of the differences in the forces
    diff = ABS ( ABS (dbg % f_hartree) - ABS (f_numer) )
    IF ( .NOT. ANY ( dbg % f_hartree == 0._dbl ) ) rel = diff / ABS (dbg % f_hartree)

! find the maximum difference and the relative and absolute errors
    WRITE ( iw, '( A,T61,G20.12 )' ) 'MAXIMUM ABSOLUTE ERROR = ', &
        MAXVAL (diff)
    IF ( .NOT. ANY ( dbg % f_hartree == 0._dbl ) ) THEN 
      WRITE ( iw, '( A,T61,G20.12 )' ) 'MAXIMUM RELATIVE ERROR = ', &
        MAXVAL (rel)
    ENDIF

! write out the particle number and forces of the max absolute
! and relative error
    dum = MAXLOC ( diff )
    WRITE ( iw, '( A,T71,I10 )' ) &
        ' PARTICLE WITH MAX ABSOLUTE ERROR IS ', dum(2)
    WRITE ( iw, '( A,T21,3G20.12 )' ) ' F ANA G-SPACE =', &
        dbg % f_hartree(:,dum(2))
    WRITE ( iw, '( A,T21,3G20.12 )' ) ' F NUM G-SPACE =', &
        f_numer(:,dum(2))

    IF ( .NOT. ANY ( dbg % f_hartree == 0._dbl ) ) THEN 
      dum = MAXLOC ( rel )
      WRITE ( iw, '( A,T71,I10 )' ) &
        'PARTICLE WITH MAX RELATIVE ERROR IS ', dum(2)
      WRITE ( iw, '( A,T21,3G20.12 )' ) ' F ANA G-SPACE =', &
        dbg % f_hartree(:,dum(2))
      WRITE ( iw, '( A,T21,3G20.12 )' ) ' F NUM G-SPACE =', &
        f_numer(:,dum(2))
    END IF
  END IF

  WRITE ( iw, '( A )' ) &
       ' DO YOU WANT TO DEBUG YOUR BONDED_CORRECT INTERACTIONS (1=yes)?'
!  READ ( ir, * ) iflag
  iflag = 5
  IF ( iflag == 1 ) THEN
! get numerical force for the total bc contribution on the particles
    WRITE ( iw, '(A)' ) ' ENTER A DELTA LESS THAN 1'
    READ ( ir, * ) delta
    IF ( delta >= 1 .OR. delta <= 0 ) THEN
      delta = 1.E-5
      WRITE ( iw, '(A,T71,F10.6)' ) &
        ' DELTA ( changed to default ) = ', delta
    ELSE
      WRITE ( iw, '(A,T71,F10.6)' ) ' DELTA = ', delta
    END IF

    CALL get_grids ( pw_grids = pw_grid, &              
                      gridlevel_info = gridlevel_info, &
                      l_info = l_info, cube_info = cube_info, &
                      ref_grid = ref_grid ) 

!    CALL fpart_bc_numer ( delta, kg_env, pw_grid, ref_grid, gridlevel_info, &
!                              l_info, cube_info, f_numer, energy_numer, globenv )

    WRITE ( iw, '( A, T61, G20.12 )' ) ' BOND_CORR ANA ENERGY = ', &
               dbg % pot_bond_corr
    WRITE ( iw, '( A, T61, G20.12 )' ) ' BOND_CORR DBG ENERGY = ', &
               energy_numer

! write out all particle numbers and forces
    WRITE ( iw, '( /,A )' ) &
     ' DO YOU WANT TO GET THE FORCES ON ALL PARTICLES (FILE=part_bc, 1=yes)?'
    READ ( ir, * ) iflag
    IF ( iflag == 1 ) THEN
       OPEN (UNIT = iwf, FILE = 'part_bc', STATUS = 'REPLACE', ACTION = 'WRITE')
       DO iatom = 1, natoms
         WRITE ( iwf, '( A,T71,I10 )' ) ' ATOM NUMBER ', iatom
         WRITE ( iwf, '( A,T21,3G20.12 )' ) &
            ' F ANA THOMAS-FERMI =', dbg % f_bc(:,iatom)
         WRITE ( iwf, '( A,T21,3G20.12,/ )' ) &
         ' F NUM THOMAS-FERMI =', f_numer(:,iatom)
       END DO
       CLOSE (iwf)
    END IF

! compute the absolute value of the differences in the forces

! compute the absolute value of the differences in the forces
    diff = ABS ( ABS (dbg % f_bc) - ABS (f_numer) )
    IF ( .NOT. ANY ( dbg % f_bc == 0._dbl ) ) rel = diff / ABS (dbg % f_bc)

! find the maximum difference and the relative and absolute errors
    WRITE ( iw, '( A,T61,G20.12 )' ) 'MAXIMUM ABSOLUTE ERROR = ', &
         MAXVAL (diff)
    IF ( .NOT. ANY ( dbg % f_bc == 0._dbl ) ) THEN 
      WRITE ( iw, '( A,T61,G20.12 )' ) 'MAXIMUM RELATIVE ERROR = ', &
         MAXVAL (rel)
    ENDIF

! write out the particle number and forces of the max absolute
! and relative error
    dum = MAXLOC ( diff )
    WRITE ( iw, '( A,T71,I10 )' ) &
         ' PARTICLE WITH MAX ABSOLUTE ERROR IS ', dum(2)
    WRITE ( iw, '( A,T21,3G20.12 )' ) ' F ANA BOND_CORR =', &
         dbg % f_bc(:,dum(2))
    WRITE ( iw, '( A,T21,3G20.12 )' ) ' F NUM BOND_CORR =', &
         f_numer(:,dum(2))

    IF ( .NOT. ANY ( dbg % f_bc == 0._dbl ) ) THEN 
      dum = MAXLOC ( rel )
      WRITE ( iw, '( A,T71,I10 )' ) &
         'PARTICLE WITH MAX RELATIVE ERROR IS ', dum(2)
      WRITE ( iw, '( A,T21,3G20.12 )' ) ' F ANA BOND_CORR =', &
         dbg % f_bc(:,dum(2))
      WRITE ( iw, '( A,T21,3G20.12 )' ) ' F NUM BOND_CORR =', &
         f_numer(:,dum(2))
    END IF
  END IF

  WRITE ( iw, '( A )' ) &
       ' DO YOU WANT TO DEBUG YOUR XC INTERACTIONS (1=yes)?'
  READ ( ir, * ) iflag
  IF ( iflag == 1 ) THEN
! get numerical force for the total xc contribution on the particles
    WRITE ( iw, '(A)' ) ' ENTER A DELTA LESS THAN 1'
    READ ( ir, * ) delta
    IF ( delta >= 1 .OR. delta <= 0 ) THEN
      delta = 1.E-5
      WRITE ( iw, '(A,T71,F10.6)' ) &
        ' DELTA ( changed to default ) = ', delta
    ELSE
      WRITE ( iw, '(A,T71,F10.6)' ) ' DELTA = ', delta
    END IF

    CALL get_grids ( pw_grids = pw_grid, &              
                      gridlevel_info = gridlevel_info, &
                      l_info = l_info, cube_info = cube_info, &
                      ref_grid = ref_grid ) 

    CALL fpart_xc_numer ( delta, kg_env, pw_grid, ref_grid, gridlevel_info, &
                              l_info, cube_info, f_numer, energy_numer, globenv )

    WRITE ( iw, '( A, T61, G20.12 )' ) ' XC ANA ENERGY = ', &
               dbg % pot_xc
    WRITE ( iw, '( A, T61, G20.12 )' ) ' XC DBG ENERGY = ', &
               energy_numer

! write out all particle numbers and forces
    WRITE ( iw, '( /,A )' ) &
     ' DO YOU WANT TO GET THE FORCES ON ALL PARTICLES (FILE=part_xc, 1=yes)?'
    READ ( ir, * ) iflag
    IF ( iflag == 1 ) THEN
       OPEN (UNIT = iwf, FILE = 'part_xc', STATUS = 'REPLACE', ACTION = 'WRITE')
       DO iatom = 1, natoms
         WRITE ( iwf, '( A,T71,I10 )' ) ' ATOM NUMBER ', iatom
         WRITE ( iwf, '( A,T21,3G20.12 )' ) &
            ' F ANA XC =', dbg % f_xc(:,iatom)
         WRITE ( iwf, '( A,T21,3G20.12,/ )' ) &
         ' F NUM XC =', f_numer(:,iatom)
       END DO
       CLOSE (iwf)
    END IF

! compute the absolute value of the differences in the forces

! compute the absolute value of the differences in the forces
    diff = ABS ( ABS (dbg % f_xc) - ABS (f_numer) )
    IF ( .NOT. ANY ( dbg % f_xc == 0._dbl ) ) rel = diff / ABS (dbg % f_xc)

! find the maximum difference and the relative and absolute errors
    WRITE ( iw, '( A,T61,G20.12 )' ) 'MAXIMUM ABSOLUTE ERROR = ', &
         MAXVAL (diff)
    IF ( .NOT. ANY ( dbg % f_xc == 0._dbl ) ) THEN 
      WRITE ( iw, '( A,T61,G20.12 )' ) 'MAXIMUM RELATIVE ERROR = ', &
         MAXVAL (rel)
    ENDIF

! write out the particle number and forces of the max absolute
! and relative error
    dum = MAXLOC ( diff )
    WRITE ( iw, '( A,T71,I10 )' ) &
         ' PARTICLE WITH MAX ABSOLUTE ERROR IS ', dum(2)
    WRITE ( iw, '( A,T21,3G20.12 )' ) ' F ANA XC =', &
         dbg % f_xc(:,dum(2))
    WRITE ( iw, '( A,T21,3G20.12 )' ) ' F NUM XC =', &
         f_numer(:,dum(2))

    IF ( .NOT. ANY ( dbg % f_xc == 0._dbl ) ) THEN 
      dum = MAXLOC ( rel )
      WRITE ( iw, '( A,T71,I10 )' ) &
         'PARTICLE WITH MAX RELATIVE ERROR IS ', dum(2)
      WRITE ( iw, '( A,T21,3G20.12 )' ) ' F ANA XC =', &
         dbg % f_xc(:,dum(2))
      WRITE ( iw, '( A,T21,3G20.12 )' ) ' F NUM XC =', &
         f_numer(:,dum(2))
    END IF
  END IF

  WRITE ( iw, '( A )' ) &
       ' DO YOU WANT TO DEBUG CORE CHARGE RSPACE INTERACTIONS (1=yes)?'
  READ ( ir, * ) iflag
  IF ( iflag == 1 ) THEN
! get numerical force for the total contribution on the particles
    WRITE ( iw, '(A)' ) ' ENTER A DELTA LESS THAN 1'
    READ ( ir, * ) delta
    IF ( delta >= 1 .OR. delta <= 0 ) THEN
      delta = 1.E-5
      WRITE ( iw, '(A,T71,F10.6)' ) &
        ' DELTA ( changed to default ) = ', delta
    ELSE
      WRITE ( iw, '(A,T71,F10.6)' ) ' DELTA = ', delta
    END IF

    CALL fpart_core_rspace_numer ( delta, kg_env, f_numer, &
                                   energy_numer, globenv )

    WRITE ( iw, '( A, T61, G20.12 )' ) ' RSPACE ANA ENERGY = ', &
               dbg % pot_rspace_core
    WRITE ( iw, '( A, T61, G20.12 )' ) ' RSPACE DBG ENERGY = ', &
               energy_numer

! write out all particle numbers and forces
    WRITE ( iw, '( /,A )' ) &
     ' DO YOU WANT TO GET THE FORCES ON ALL PARTICLES (FILE=part_rspace, 1=yes)?'
    READ ( ir, * ) iflag
    IF ( iflag == 1 ) THEN
       OPEN (UNIT = iwf, FILE = 'part_rspace', STATUS = 'REPLACE', ACTION = 'WRITE')
       DO iatom = 1, natoms
         WRITE ( iwf, '( A,T71,I10 )' ) ' ATOM NUMBER ', iatom
         WRITE ( iwf, '( A,T21,3G20.12 )' ) &
            ' F ANA RSPACE CORE =', dbg % f_rspace_core(:,iatom)
         WRITE ( iwf, '( A,T21,3G20.12,/ )' ) &
         ' F NUM RSPACE CORE =', f_numer(:,iatom)
       END DO
       CLOSE (iwf)
    END IF

! compute the absolute value of the differences in the forces

! compute the absolute value of the differences in the forces
    diff = ABS ( ABS (dbg % f_rspace_core) - ABS (f_numer) )
    IF ( .NOT. ANY ( dbg % f_rspace_core == 0._dbl ) ) rel = diff / ABS (dbg % f_rspace_core)

! find the maximum difference and the relative and absolute errors
    WRITE ( iw, '( A,T61,G20.12 )' ) 'MAXIMUM ABSOLUTE ERROR = ', &
         MAXVAL (diff)
    IF ( .NOT. ANY ( dbg % f_rspace_core == 0._dbl ) ) THEN 
      WRITE ( iw, '( A,T61,G20.12 )' ) 'MAXIMUM RELATIVE ERROR = ', &
         MAXVAL (rel)
    ENDIF

! write out the particle number and forces of the max absolute
! and relative error
    dum = MAXLOC ( diff )
    WRITE ( iw, '( A,T71,I10 )' ) &
         ' PARTICLE WITH MAX ABSOLUTE ERROR IS ', dum(2)
    WRITE ( iw, '( A,T21,3G20.12 )' ) ' F ANA RSPACE CORE =', &
         dbg % f_rspace_core(:,dum(2))
    WRITE ( iw, '( A,T21,3G20.12 )' ) ' F NUM RSPACE CORE =', &
         f_numer(:,dum(2))

    IF ( .NOT. ANY ( dbg % f_rspace_core == 0._dbl ) ) THEN 
      dum = MAXLOC ( rel )
      WRITE ( iw, '( A,T71,I10 )' ) &
         'PARTICLE WITH MAX RELATIVE ERROR IS ', dum(2)
      WRITE ( iw, '( A,T21,3G20.12 )' ) ' F ANA RSPACE CORE =', &
         dbg % f_rspace_core(:,dum(2))
      WRITE ( iw, '( A,T21,3G20.12 )' ) ' F NUM RSPACE CORE =', &
         f_numer(:,dum(2))
    END IF
  END IF

  END SUBROUTINE debug_control

!!*****
!******************************************************************************

END MODULE kg_debug

!******************************************************************************
