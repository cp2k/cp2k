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

  USE coefficient_types, ONLY : coeff_type
  USE cube_utils, ONLY : cube_info_type
  Use cp_control_types, ONLY : dft_control_type
  USE gaussian_gridlevels, ONLY: gridlevel_info_type
  USE kg_environment_types, ONLY : get_kg_env,kg_environment_type
  USE kg_force_numer, ONLY : fpart_gspace_numer, fpart_bc_numer,&
                             fpart_xc_numer, fpart_core_rspace_numer,&
                             fpart_p_XC_numer, fpart_p_gspace_numer,&
                             fcoef_p_XC_numer, fcoef_p_gspace_numer,&
                             fpart_ppl_numer, fpart_p_ppl_numer, &
                             fcoef_p_ppl_numer
  USE kg_force, ONLY : kg_force_control !!!, debug_variables_type, get_grids
  USE global_types, ONLY : global_environment_type
  USE kinds, ONLY : dbl
  USE l_utils, ONLY : l_info_type
  USE linklists, ONLY : bonds, bends, torsions, dist_constraints, &
       g3x3_constraints
  USE molecule_types, ONLY : molecule_structure_type, particle_node_type
  USE particle_types, ONLY : particle_type
  USE pw_grid_types, ONLY : pw_grid_type
  USE pw_pool_types, ONLY : pw_pool_type 
  USE simulation_cell, ONLY : cell_type
  USE termination, ONLY : stop_memory, stop_program

  IMPLICIT NONE

! debug variables
! gt this is to change! Here only for build purposes
   TYPE debug_variables_type
     REAL ( dbl ) :: pot_bond_corr, pot_hartree, pot_xc, &
                     pot_rspace_core, pot_ppl, pot_ppl2
     REAL ( dbl ), DIMENSION ( :, : ), POINTER :: f_hartree, &
           f_rspace_core, f_bc, f_xc, f_ppl
     REAL ( dbl ), DIMENSION ( : ), POINTER :: fc_hartree,&
                                               fc_ppl, fc_xc
  END TYPE debug_variables_type

  PRIVATE
  PUBLIC :: kg_debug_control, debug_variables_type

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** kg_debug/kg_debug_control [1.0] *
!!
!!   NAME
!!     kg_debug_control
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

SUBROUTINE  kg_debug_control ( globenv, kg_env )

  IMPLICIT NONE

! Arguments
  TYPE ( kg_environment_type ), POINTER   :: kg_env
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv

! Locals

  TYPE ( cell_type ), POINTER :: cell
  TYPE ( coeff_type ), POINTER :: coef_pos
  TYPE ( coeff_type ), POINTER :: coef_vel
  TYPE ( coeff_type ), POINTER :: coef_force
  TYPE ( dft_control_type), POINTER:: dft_control  
  TYPE ( cube_info_type ), POINTER :: cube_info ( : ) 
!dbg gt to be changed
   TYPE ( debug_variables_type ) :: dbg
!dbg gt to be done
  TYPE ( gridlevel_info_type ), POINTER :: gridlevel_info
  TYPE ( l_info_type ), POINTER    :: l_info
  TYPE ( particle_type ), POINTER  :: particle_set ( : )
  TYPE ( pw_grid_type ), POINTER, DIMENSION ( : ) :: pw_grid
  TYPE ( pw_pool_type ), POINTER :: pw_pool
  INTEGER :: iflag, i, natoms, isos, iatom, iw, iwf, ir, ref_grid, gmax
  INTEGER :: NCOEF, icoef
  INTEGER, DIMENSION ( 2 ) :: dum
  INTEGER, DIMENSION ( 1 ) :: dum_c
  REAL ( dbl ) :: delta, energy_numer
  REAL ( dbl ) :: err1, numer, denom1, vec ( 3 ), e_bc, e_real, energy_tot
  REAL ( dbl ) :: denom2, err2
  REAL ( dbl ), DIMENSION ( :, : ), ALLOCATABLE :: f_bc, f_real
  REAL ( dbl ), DIMENSION ( :, : ), ALLOCATABLE :: f_numer
  REAL ( dbl ), DIMENSION ( :, : ), ALLOCATABLE :: rel, diff
  REAL ( dbl ), DIMENSION (:), ALLOCATABLE :: fc_numer, rel_c, diff_c
  LOGICAL :: kgpol 

!------------------------------------------------------------------------------

! Assigning local pointers
  NULLIFY (dft_control, cell, particle_set)
  CALL get_kg_env (kg_env=kg_env, particle_set=particle_set, cell=cell,&
                   dft_control=dft_control)

  kgpol = dft_control % kg_control % polarization
  iw = globenv % scr
  iwf = 7
  ir = 5
  WRITE ( iw, '( 16x,A )' ) &
       '***************************************************'
  WRITE ( iw, '( 16x,A )' ) '                   BEGIN KIM-GORDON DEBUG'
  WRITE ( iw, '( 16x,A )' ) &
       '***************************************************'

! allocating all arrays involved with debug
  natoms = SIZE ( particle_set )
  IF ( .NOT. ALLOCATED ( diff ) ) ALLOCATE ( diff ( 3, natoms ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( "kg_debug_control", "diff", 3 * natoms )
  IF ( .NOT. ALLOCATED ( rel ) ) ALLOCATE ( rel ( 3, natoms ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( "kg_debug_control", "rel", 3 * natoms )
  IF ( .NOT. ALLOCATED ( f_numer ) ) &
       ALLOCATE ( f_numer ( 3, natoms ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( "kg_debug_control", "f_numer", 3 * natoms )
  ALLOCATE ( dbg % f_hartree ( 3, natoms ), STAT = isos )
  IF ( isos /= 0 ) &
       CALL stop_memory ( "kg_debug_control", "f_hartree", 3 * natoms )
  ALLOCATE ( dbg % f_bc ( 3, natoms ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( "kg_debug_control", "f_bc", 3 * natoms )
  ALLOCATE ( dbg % f_ppl ( 3, natoms ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( "kg_debug_control", "f_ppl", 3 * natoms )
  ALLOCATE ( dbg % f_xc ( 3, natoms ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( "kg_debug_control", "f_xc", 3 * natoms )
  ALLOCATE ( dbg % f_rspace_core ( 3, natoms ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( "kg_debug_control", "f_rspace_core", 3 * natoms )

  IF (kgpol) THEN
! gt this part is to be done!
!    coef_pos => kg_env% coef_pos (1)
!    coef_force => kg_env% coef_force (1)
    ncoef= size(coef_pos%ao%cr)
    IF ( .NOT. ALLOCATED ( diff_c ) ) ALLOCATE ( diff_c ( ncoef ), STAT = isos )
    IF ( isos /= 0 ) CALL stop_memory ( "kg_debug_control", "diff_c", ncoef )
    IF ( .NOT. ALLOCATED ( rel_c ) ) ALLOCATE ( rel_c ( ncoef ), STAT = isos )
    IF ( isos /= 0 ) CALL stop_memory ( "kg_debug_control", "rel_c",  ncoef )
    IF ( .NOT. ALLOCATED ( fc_numer ) ) &
        ALLOCATE ( fc_numer (  ncoef ), STAT = isos )
    IF ( isos /= 0 ) CALL stop_memory ( "kg_debug_control", "fc_numer", ncoef )
    ALLOCATE ( dbg % fc_hartree ( ncoef ), STAT = isos )
    IF ( isos /= 0 ) &
        CALL stop_memory ( "kg_debug_control", "fc_hartree", ncoef )
    ALLOCATE ( dbg % fc_xc ( ncoef ), STAT = isos )
    IF ( isos /= 0 ) CALL stop_memory ( "kg_debug_control", "fc_xc", ncoef )
    ALLOCATE ( dbg % fc_ppl ( ncoef ), STAT = isos )
    IF ( isos /= 0 ) CALL stop_memory ( "kg_debug_control", "fc_ppl", ncoef )
  END IF

! Call to analytical force routine:
  CALL kg_force_control (  kg_env, globenv )

  WRITE ( iw, '( A )' ) &
       ' DO YOU WANT TO DEBUG YOUR G-SPACE INTERACTIONS (1=yes)?'
  READ ( ir, * ) iflag
  IF ( iflag == 1 ) THEN
! get numerical force for the total gspace contribution on the particles
    WRITE ( iw, '(A)' ) ' ENTER A DELTA LESS THAN 1'
    READ ( ir, * ) delta
    IF ( delta >= 1 .OR. delta <= 0 ) THEN
      delta = 1.D-5
      WRITE ( iw, '(A,T71,F10.6)' ) &
        ' DELTA ( changed to default ) = ', delta
    ELSE
      WRITE ( iw, '(A,T71,F10.6)' ) ' DELTA = ', delta
    END IF

!    CALL get_grids ( pw_grids = pw_grid, &              
!                      gridlevel_info = gridlevel_info, &
!                      l_info = l_info, cube_info = cube_info, &
!                      ref_grid = ref_grid ) 
    IF (kgpol) THEN
      CALL fpart_p_gspace_numer ( delta, kg_env, pw_grid, ref_grid, gridlevel_info, &
                                l_info, cube_info, f_numer, energy_numer, globenv )

    ELSE
      CALL fpart_gspace_numer ( delta, kg_env, pw_grid, ref_grid, gridlevel_info, &
                                l_info, cube_info, f_numer, energy_numer, globenv )
    END IF

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
  READ ( ir, * ) iflag
  IF ( iflag == 1 ) THEN
! get numerical force for the total bc contribution on the particles
    WRITE ( iw, '(A)' ) ' ENTER A DELTA LESS THAN 1'
    READ ( ir, * ) delta
    IF ( delta >= 1 .OR. delta <= 0 ) THEN
      delta = 1.D-5
      WRITE ( iw, '(A,T71,F10.6)' ) &
        ' DELTA ( changed to default ) = ', delta
    ELSE
      WRITE ( iw, '(A,T71,F10.6)' ) ' DELTA = ', delta
    END IF

    CALL fpart_bc_numer ( delta, kg_env, f_numer, energy_numer, globenv) 

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
      delta = 1.D-5
      WRITE ( iw, '(A,T71,F10.6)' ) &
        ' DELTA ( changed to default ) = ', delta
    ELSE
      WRITE ( iw, '(A,T71,F10.6)' ) ' DELTA = ', delta
    END IF

!    CALL get_grids ( pw_grids = pw_grid, pw_pool=pw_pool, &              
!                      gridlevel_info = gridlevel_info, &
!                      l_info = l_info, cube_info = cube_info, &
!                      ref_grid = ref_grid ) 
    IF (kgpol) THEN

      CALL fpart_p_xc_numer ( delta, kg_env, pw_grid, pw_pool, ref_grid, gridlevel_info, &
                                l_info, cube_info, f_numer, energy_numer, globenv )
    ELSE
      CALL fpart_xc_numer ( delta, kg_env, pw_grid, pw_pool, ref_grid, gridlevel_info, &
                                l_info, cube_info, f_numer, energy_numer, globenv )
    END IF

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
      delta = 1.D-5
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

  WRITE ( iw, '( A )' ) &
       ' DO YOU WANT TO DEBUG LOCAL PP INTERACTIONS (1=yes)?'
  READ ( ir, * ) iflag
  IF ( iflag == 1 ) THEN
! get numerical force for the total contribution on the particles
    WRITE ( iw, '(A)' ) ' ENTER A DELTA LESS THAN 1'
    READ ( ir, * ) delta
    IF ( delta >= 1 .OR. delta <= 0 ) THEN
      delta = 1.D-5
      WRITE ( iw, '(A,T71,F10.6)' ) &
        ' DELTA ( changed to default ) = ', delta
    ELSE
      WRITE ( iw, '(A,T71,F10.6)' ) ' DELTA = ', delta
    END IF

    IF(kgpol) THEN
      CALL fpart_p_ppl_numer ( delta, kg_env, f_numer, &
                                   energy_numer, globenv )
      WRITE ( iw, '( A, T61, G20.12 )' ) 'PPL ANA ENERGY = ', &
               dbg % pot_ppl + dbg % pot_ppl2

    ELSE
      CALL fpart_ppl_numer ( delta, kg_env, f_numer, &
                                   energy_numer, globenv )
      WRITE ( iw, '( A, T61, G20.12 )' ) 'PPL ANA ENERGY = ', &
               dbg % pot_ppl
    END IF

    WRITE ( iw, '( A, T61, G20.12 )' ) 'PPL DBG ENERGY = ', &
               energy_numer

! write out all particle numbers and forces
    WRITE ( iw, '( /,A )' ) &
     ' DO YOU WANT TO GET THE FORCES ON ALL PARTICLES (FILE=part_ppl, 1=yes)?'
    READ ( ir, * ) iflag
    IF ( iflag == 1 ) THEN
       OPEN (UNIT = iwf, FILE = 'part_ppl', STATUS = 'REPLACE', ACTION = 'WRITE')
       DO iatom = 1, natoms
         WRITE ( iwf, '( A,T71,I10 )' ) ' ATOM NUMBER ', iatom
         WRITE ( iwf, '( A,T21,3G20.12 )' ) &
            ' F ANA PPL =', dbg % f_ppl(:,iatom)
         WRITE ( iwf, '( A,T21,3G20.12,/ )' ) &
         ' F NUM PPL =', f_numer(:,iatom)
       END DO
       CLOSE (iwf)
    END IF

! compute the absolute value of the differences in the forces

! compute the absolute value of the differences in the forces
    diff = ABS ( ABS (dbg % f_ppl) - ABS (f_numer) )
    IF ( .NOT. ANY ( dbg % f_ppl == 0._dbl ) ) rel = diff / ABS (dbg % f_ppl)

! find the maximum difference and the relative and absolute errors
    WRITE ( iw, '( A,T61,G20.12 )' ) 'MAXIMUM ABSOLUTE ERROR = ', &
         MAXVAL (diff)
    IF ( .NOT. ANY ( dbg % f_ppl == 0._dbl ) ) THEN 
      WRITE ( iw, '( A,T61,G20.12 )' ) 'MAXIMUM RELATIVE ERROR = ', &
         MAXVAL (rel)
    ENDIF

! write out the particle number and forces of the max absolute
! and relative error
    dum = MAXLOC ( diff )
    WRITE ( iw, '( A,T71,I10 )' ) &
         ' PARTICLE WITH MAX ABSOLUTE ERROR IS ', dum(2)
    WRITE ( iw, '( A,T21,3G20.12 )' ) ' F ANA PPL =', &
         dbg % f_ppl(:,dum(2))
    WRITE ( iw, '( A,T21,3G20.12 )' ) ' F NUM PPL =', &
         f_numer(:,dum(2))

    IF ( .NOT. ANY ( dbg % f_ppl == 0._dbl ) ) THEN 
      dum = MAXLOC ( rel )
      WRITE ( iw, '( A,T71,I10 )' ) &
         'PARTICLE WITH MAX RELATIVE ERROR IS ', dum(2)
      WRITE ( iw, '( A,T21,3G20.12 )' ) ' F ANA PPL =', &
         dbg % f_ppl(:,dum(2))
      WRITE ( iw, '( A,T21,3G20.12 )' ) ' F NUM PPL =', &
         f_numer(:,dum(2))
    END IF
  END IF


  IF (kgpol) THEN
   WRITE ( iw, '( A )' ) &
         ' DO YOU WANT TO DEBUG YOUR G-SPACE INTERACTIONS on COEFS (1=yes)?'
   READ ( ir, * ) iflag
   IF ( iflag == 1 ) THEN
! get numerical force for the total gspace contribution on the COEFS
    WRITE ( iw, '(A)' ) ' ENTER A DELTA LESS THAN 1'
    READ ( ir, * ) delta
    IF ( delta >= 1 .OR. delta <= 0 ) THEN
      delta = 1.D-5
      WRITE ( iw, '(A,T71,F10.6)' ) &
        ' DELTA ( changed to default ) = ', delta
    ELSE
      WRITE ( iw, '(A,T71,F10.6)' ) ' DELTA = ', delta
    END IF

!    CALL get_grids ( pw_grids = pw_grid, &
!                      gridlevel_info = gridlevel_info, &
!                      l_info = l_info, cube_info = cube_info, &
!                      ref_grid = ref_grid )
    CALL fcoef_p_gspace_numer ( delta, kg_env, pw_grid, ref_grid, gridlevel_info, &
                                l_info, cube_info, fc_numer, energy_numer, globenv )

    WRITE ( iw, '( A, T61, G20.12 )' ) ' G-SPACE ANA ENERGY = ', &
               dbg % pot_hartree
    WRITE ( iw, '( A, T61, G20.12 )' ) ' G-SPACE DBG ENERGY = ', &
               energy_numer

! write out all particle numbers and forces
    WRITE ( iw, '( /,A )' ) &
     ' DO YOU WANT TO GET THE FORCES ON ALL coefs     (FILE=coef_gspace, 1=yes)?'
    READ ( ir, * ) iflag
    IF ( iflag == 1 ) THEN
       OPEN (UNIT = iwf, FILE = 'coef_gspace', STATUS = 'REPLACE', ACTION = 'WRITE')
       DO icoef = 1, ncoef
         WRITE ( iwf, '( A,T71,I10 )' ) ' ATOM NUMBER ', icoef
         WRITE ( iwf, '( A,T21,G20.12 )' ) &
            ' F ANA GSPACE =', dbg % fc_hartree(icoef)
         WRITE ( iwf, '( A,T21,G20.12,/ )' ) &
         ' F NUM GSPACE =', fc_numer(icoef)
       END DO
       CLOSE (iwf)
    END IF

! compute the absolute value of the differences in the forces
! compute the absolute value of the differences in the forces
    diff_c = ABS ( ABS (dbg % fc_hartree) - ABS (fc_numer) )
    IF ( .NOT. ANY ( dbg % fc_hartree == 0._dbl ) ) rel_c = diff_c / ABS (dbg % fc_hartree)

! find the maximum difference and the relative and absolute errors
    WRITE ( iw, '( A,T61,G20.12 )' ) 'MAXIMUM ABSOLUTE ERROR = ', &
        MAXVAL (diff_c)
    IF ( .NOT. ANY ( dbg % fc_hartree == 0._dbl ) ) THEN
      WRITE ( iw, '( A,T61,G20.12 )' ) 'MAXIMUM RELATIVE ERROR = ', &
        MAXVAL (rel_c)
    ENDIF

! write out the particle number and forces of the max absolute
! and relative error
    dum_c = MAXLOC ( diff_c )
    WRITE ( iw, '( A,T71,I10 )' ) &
        ' PARTICLE WITH MAX ABSOLUTE ERROR IS ', dum_c
    WRITE ( iw, '( A,T21,3G20.12 )' ) ' F ANA G-SPACE =', &
        dbg % fc_hartree(dum_c)
    WRITE ( iw, '( A,T21,3G20.12 )' ) ' F NUM G-SPACE =', &
        fc_numer(dum_c)

    IF ( .NOT. ANY ( dbg % fc_hartree == 0._dbl ) ) THEN
      dum_c = MAXLOC ( rel_c )
      WRITE ( iw, '( A,T71,I10 )' ) &
        'PARTICLE WITH MAX RELATIVE ERROR IS ', dum_c
      WRITE ( iw, '( A,T21,3G20.12 )' ) ' F ANA G-SPACE =', &
        dbg % fc_hartree(dum_c)
      WRITE ( iw, '( A,T21,3G20.12 )' ) ' F NUM G-SPACE =', &
        fc_numer(dum_c)
    END IF

   END IF

   WRITE ( iw, '( A )' ) &
       ' DO YOU WANT TO DEBUG YOUR XC INTERACTIONS ON COEF (1=yes)?'
   READ ( ir, * ) iflag
   IF ( iflag == 1 ) THEN
! get numerical force for the total xc contribution on the COEFS
    WRITE ( iw, '(A)' ) ' ENTER A DELTA LESS THAN 1'
    READ ( ir, * ) delta
    IF ( delta >= 1 .OR. delta <= 0 ) THEN
      delta = 1.D-5
      WRITE ( iw, '(A,T71,F10.6)' ) &
        ' DELTA ( changed to default ) = ', delta
    ELSE
      WRITE ( iw, '(A,T71,F10.6)' ) ' DELTA = ', delta
    END IF

!    CALL get_grids ( pw_grids = pw_grid, pw_pool = pw_pool, & 
!                      gridlevel_info = gridlevel_info, &
!                      l_info = l_info, cube_info = cube_info, &
!                      ref_grid = ref_grid )

    CALL fcoef_p_xc_numer ( delta, kg_env, pw_grid, pw_pool, ref_grid, gridlevel_info, &
                                l_info, cube_info, fc_numer, energy_numer, globenv )

    WRITE ( iw, '( A, T61, G20.12 )' ) ' XC ANA ENERGY = ', &
               dbg % pot_xc
    WRITE ( iw, '( A, T61, G20.12 )' ) ' XC DBG ENERGY = ', &
               energy_numer

! write out all particle numbers and forces
    WRITE ( iw, '( /,A )' ) &
     ' DO YOU WANT TO GET THE FORCES ON ALL COEfS (FILE=coef_xc, 1=yes)?'
    READ ( ir, * ) iflag
    IF ( iflag == 1 ) THEN
       OPEN (UNIT = iwf, FILE = 'coef_xc', STATUS = 'REPLACE', ACTION = 'WRITE')
       DO icoef = 1, ncoef 
         WRITE ( iwf, '( A,T71,I10 )' ) ' coef NUMBER ', icoef
         WRITE ( iwf, '( A,T21,3G20.12 )' ) &
            ' F ANA XC =', dbg % fc_xc(icoef)
         WRITE ( iwf, '( A,T21,3G20.12,/ )' ) &
         ' F NUM XC =', fc_numer(icoef)
       END DO
       CLOSE (iwf)
    END IF
! compute the absolute value of the differences in the forces
    diff_c = ABS ( ABS (dbg % fc_xc) - ABS (fc_numer) )
    IF ( .NOT. ANY ( dbg % fc_xc == 0._dbl ) ) rel_c = diff_c / ABS (dbg % fc_xc)

! find the maximum difference and the relative and absolute errors
    WRITE ( iw, '( A,T61,G20.12 )' ) 'MAXIMUM ABSOLUTE ERROR = ', &
        MAXVAL (diff_c)
    IF ( .NOT. ANY ( dbg % fc_xc == 0._dbl ) ) THEN
      WRITE ( iw, '( A,T61,G20.12 )' ) 'MAXIMUM RELATIVE ERROR = ', &
        MAXVAL (rel_c)
    ENDIF

! write out the particle number and forces of the max absolute
! and relative error
    dum_c = MAXLOC ( diff_c )
    WRITE ( iw, '( A,T71,I10 )' ) &
        ' PARTICLE WITH MAX ABSOLUTE ERROR IS ', dum_c
    WRITE ( iw, '( A,T21,3G20.12 )' ) ' F ANA XC =', &
        dbg % fc_xc(dum_c)
    WRITE ( iw, '( A,T21,3G20.12 )' ) ' F NUM XC =', &
        fc_numer(dum_c)

    IF ( .NOT. ANY ( dbg % fc_xc == 0._dbl ) ) THEN
      dum_c = MAXLOC ( rel_c )
      WRITE ( iw, '( A,T71,I10 )' ) &
        'PARTICLE WITH MAX RELATIVE ERROR IS ', dum_c
      WRITE ( iw, '( A,T21,3G20.12 )' ) ' F ANA XC =', &
        dbg % fc_xc(dum_c)
      WRITE ( iw, '( A,T21,3G20.12 )' ) ' F NUM XC =', &
        fc_numer(dum_c)
    END IF

   END IF


   WRITE ( iw, '( A )' ) &
       ' DO YOU WANT TO DEBUG YOUR PPL INTERACTIONS ON COEF (1=yes)?'
   READ ( ir, * ) iflag
   IF ( iflag == 1 ) THEN
! get numerical force for the total ppl contribution on the COEFS
    WRITE ( iw, '(A)' ) ' ENTER A DELTA LESS THAN 1'
    READ ( ir, * ) delta
    IF ( delta >= 1 .OR. delta <= 0 ) THEN
      delta = 1.D-5
      WRITE ( iw, '(A,T71,F10.6)' ) &
        ' DELTA ( changed to default ) = ', delta
    ELSE
      WRITE ( iw, '(A,T71,F10.6)' ) ' DELTA = ', delta
    END IF

    CALL fcoef_p_ppl_numer ( delta, kg_env, &
                                 fc_numer, energy_numer, globenv )

    WRITE ( iw, '( A, T61, G20.12 )' ) ' PPL ANA ENERGY = ', &
               dbg % pot_ppl2
    WRITE ( iw, '( A, T61, G20.12 )' ) ' PPL DBG ENERGY = ', &
               energy_numer

! write out all particle numbers and forces
    WRITE ( iw, '( /,A )' ) &
     ' DO YOU WANT TO GET THE FORCES ON ALL COEfS (FILE=coef_ppl, 1=yes)?'
    READ ( ir, * ) iflag
    IF ( iflag == 1 ) THEN
       OPEN (UNIT = iwf, FILE = 'coef_ppl', STATUS = 'REPLACE', ACTION = 'WRITE')
       DO icoef = 1, ncoef 
         WRITE ( iwf, '( A,T71,I10 )' ) ' coef NUMBER ', icoef
         WRITE ( iwf, '( A,T21,3G20.12 )' ) &
            ' F ANA ppl =', dbg % fc_ppl(icoef)
         WRITE ( iwf, '( A,T21,3G20.12,/ )' ) &
         ' F NUM ppl =', fc_numer(icoef)
       END DO
       CLOSE (iwf)
    END IF
! compute the absolute value of the differences in the forces
    diff_c = ABS ( ABS (dbg % fc_ppl) - ABS (fc_numer) )
    IF ( .NOT. ANY ( dbg % fc_ppl == 0._dbl ) ) rel_c = diff_c / ABS (dbg % fc_ppl)

! find the maximum difference and the relative and absolute errors
    WRITE ( iw, '( A,T61,G20.12 )' ) 'MAXIMUM ABSOLUTE ERROR = ', &
        MAXVAL (diff_c)
    IF ( .NOT. ANY ( dbg % fc_ppl == 0._dbl ) ) THEN
      WRITE ( iw, '( A,T61,G20.12 )' ) 'MAXIMUM RELATIVE ERROR = ', &
        MAXVAL (rel_c)
    ENDIF

! write out the particle number and forces of the max absolute
! and relative error
    dum_c = MAXLOC ( diff_c )
    WRITE ( iw, '( A,T71,I10 )' ) &
        ' PARTICLE WITH MAX ABSOLUTE ERROR IS ', dum_c
    WRITE ( iw, '( A,T21,3G20.12 )' ) ' F ANA PPL =', &
        dbg % fc_ppl(dum_c)
    WRITE ( iw, '( A,T21,3G20.12 )' ) ' F NUM PPL =', &
        fc_numer(dum_c)

    IF ( .NOT. ANY ( dbg % fc_ppl == 0._dbl ) ) THEN
      dum_c = MAXLOC ( rel_c )
      WRITE ( iw, '( A,T71,I10 )' ) &
        'PARTICLE WITH MAX RELATIVE ERROR IS ', dum_c
      WRITE ( iw, '( A,T21,3G20.12 )' ) ' F ANA PPL =', &
        dbg % fc_ppl(dum_c)
      WRITE ( iw, '( A,T21,3G20.12 )' ) ' F NUM PPL =', &
        fc_numer(dum_c)
    END IF

   END IF


  END IF


  END SUBROUTINE kg_debug_control

!!*****
!******************************************************************************

END MODULE kg_debug

!******************************************************************************
