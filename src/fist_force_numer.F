!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2005 CP2K developers group                           !
!-----------------------------------------------------------------------------!
!!****** cp2k/fist_force_numer [1.0] *
!!
!!   NAME
!!     fist_force_numer
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     Torsions added (DG) 05-Dec-2000
!!     Variable names changed (DG) 05-Dec-2000
!!
!!   SOURCE
!******************************************************************************

MODULE fist_force_numer

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE cell_types,                      ONLY: &
       cell_create, cell_release, cell_retain, cell_type, get_cell, get_hinv, &
       pbc, plane_distance, real_to_scaled, scaled_to_real
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE ewald_environment_types,         ONLY: ewald_env_get,&
                                             ewald_environment_type
  USE ewald_pw_types,                  ONLY: ewald_pw_get,&
                                             ewald_pw_type
  USE fist_neighbor_list_types,        ONLY: &
       first_list, first_node, get_neighbor_list, get_neighbor_list_set, &
       get_neighbor_node, neighbor_list_set_p_type, neighbor_list_set_type, &
       neighbor_list_type, neighbor_node_type, next
  USE fist_nonbond_env_types,          ONLY: fist_nonbond_env_get,&
                                             fist_nonbond_env_type
  USE global_types,                    ONLY: global_environment_type
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp,&
                                             dp_size,&
                                             int_size
  USE message_passing,                 ONLY: mp_sum
  USE mol_force,                       ONLY: force_bends,&
                                             force_bonds,&
                                             force_imp_torsions,&
                                             force_torsions
  USE molecule_kind_types,             ONLY: bend_type,&
                                             bond_type,&
                                             get_molecule_kind,&
                                             impr_type,&
                                             molecule_kind_type,&
                                             onfo_type,&
                                             torsion_type,&
                                             ub_type
  USE molecule_types_new,              ONLY: get_molecule,&
                                             molecule_type
  USE pair_potential,                  ONLY: potential_f,&
                                             potential_s
  USE pair_potential_types,            ONLY: pair_potential_pp_type
  USE particle_types,                  ONLY: particle_type
  USE pme,                             ONLY: pme_evaluate
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_grids,                        ONLY: pw_grid_change
  USE pw_pool_types,                   ONLY: pw_pool_type
  USE qs_neighbor_lists,               ONLY: local_atoms_type
  USE splines,                         ONLY: spline_data_type
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"  

  PRIVATE
  PUBLIC :: force_bond_numer, force_bend_numer, force_torsion_numer, &
            force_nonbond_numer, ewald_force_numer, force_imptors_numer
  PUBLIC :: pv_bond_numer, pv_bend_numer, pv_torsion_numer, pv_nonbond_numer
  PUBLIC :: pv_g_numer, force_g_numer, ewald_energy_numer
  PUBLIC :: pv_ub_numer, force_ub_numer
  PUBLIC :: pv_onefour_numer, force_onefour_numer
  PUBLIC :: pme_force_numer, pv_imptors_numer

! *** Global parameters (in this module) ***

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'fist_force_numer'

! *** Public subroutines ***


!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** fist_force_numer/force_bond_numer [1.0] *
!!
!!   NAME
!!     force_bond_numer
!!
!!   FUNCTION
!!     Computes numerically the forces from the bonds
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE force_bond_numer ( molecule_set, molecule_kind_set, local_molecules, &
                              particle_set, f_bond_numer, pot_bond, numerical_shift )


    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: f_bond_numer
    REAL(KIND=dp), INTENT(INOUT)             :: pot_bond
    REAL(KIND=dp), INTENT(IN)                :: numerical_shift

    INTEGER                                  :: first_atom, i, ibond, idim, &
                                                ikind, imol, index_a, &
                                                index_b, nbonds, nkind, &
                                                nmol_per_kind
    REAL(KIND=dp)                            :: e0, em, ep, fscalar
    REAL(KIND=dp), DIMENSION(3)              :: delta, rij, rij_minus_delta, &
                                                rij_plus_delta
    TYPE(bond_type), POINTER                 :: bond_list( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule

!------------------------------------------------------------------------------

  nkind = SIZE ( molecule_kind_set )
  delta = 0.0_dp
  MOL: DO ikind = 1, nkind
     nmol_per_kind = local_molecules % n_el ( ikind )
     DO imol = 1, nmol_per_kind
       i = local_molecules % list ( ikind ) % array ( imol )
       molecule => molecule_set ( i )
       molecule_kind => molecule % molecule_kind
       CALL get_molecule_kind ( molecule_kind, nbond = nbonds,  &
                                bond_list = bond_list )
     BOND: DO ibond = 1, nbonds
        CALL get_molecule ( molecule, first_atom=first_atom )
        index_a = bond_list ( ibond ) % a + first_atom - 1
        index_b = bond_list ( ibond ) % b + first_atom - 1
        rij = particle_set ( index_a ) % r - particle_set ( index_b ) % r
! get the energy
        CALL force_bonds ( rij, bond_list ( ibond ) % bond_kind % r0, &
             bond_list ( ibond ) % bond_kind  % k, e0, fscalar )
        pot_bond = pot_bond + e0

! get the forces
        DO idim = 1, 3
           delta(idim) = numerical_shift
           rij_minus_delta = rij - delta
           CALL force_bonds ( rij_minus_delta, bond_list ( ibond ) % bond_kind % r0, &
             bond_list ( ibond ) % bond_kind  % k, em, fscalar )
           rij_plus_delta = rij + delta
           CALL force_bonds ( rij_plus_delta, bond_list ( ibond ) % bond_kind % r0, &
                bond_list ( ibond ) % bond_kind % k, ep, fscalar )
           f_bond_numer(idim,index_a) = f_bond_numer(idim,index_a) + (em-ep)
           f_bond_numer(idim,index_b) = f_bond_numer(idim,index_b) - (em-ep)
           delta = 0.0_dp
        END DO

     END DO BOND
     END DO
    END DO MOL

    f_bond_numer = f_bond_numer / ( 2.0_dp * numerical_shift )

END SUBROUTINE force_bond_numer
!******************************************************************************
!!****** fist_force_numer/force_onefour_numer [1.0] *
!!
!!   NAME
!!     force_onefour_numer
!!
!!   FUNCTION
!!     Computes numerically the forces from the onefours
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE force_onefour_numer ( molecule_set, molecule_kind_set, local_molecules, &
                                particle_set, f_numer, pot, numerical_shift )


    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: f_numer
    REAL(KIND=dp), INTENT(INOUT)             :: pot
    REAL(KIND=dp), INTENT(IN)                :: numerical_shift

    INTEGER :: akind, bkind, first_atom, i, ibond, idim, ikind, imol, &
      index_a, index_b, nkind, nmol_per_kind, nonfo
    REAL(KIND=dp)                            :: e0, em, ep, fscalar, rab2, &
                                                rab2_minus_delta, &
                                                rab2_plus_delta
    REAL(KIND=dp), DIMENSION(3)              :: delta, rij
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule
    TYPE(onfo_type), POINTER                 :: onfo_list( : )
    TYPE(spline_data_type), POINTER          :: spline_data

!------------------------------------------------------------------------------
! initializing the forces

  nkind = SIZE ( molecule_kind_set )
  delta = 0.0_dp
  MOL: DO ikind = 1, nkind
     nmol_per_kind = local_molecules % n_el ( ikind )
     DO imol = 1, nmol_per_kind
       i = local_molecules % list ( ikind ) % array ( imol )
       molecule => molecule_set ( i )
       molecule_kind => molecule % molecule_kind
       CALL get_molecule_kind ( molecule_kind, nonfo = nonfo,  &
                               onfo_list = onfo_list )
     BOND: DO ibond = 1, nonfo
        CALL get_molecule ( molecule, first_atom=first_atom )
        index_a = onfo_list ( ibond ) % a + first_atom - 1
        atomic_kind => particle_set ( index_a ) % atomic_kind
        CALL get_atomic_kind ( atomic_kind, kind_number =  akind )
        index_b = onfo_list ( ibond ) % b + first_atom - 1
        atomic_kind => particle_set ( index_b ) % atomic_kind
        CALL get_atomic_kind ( atomic_kind, kind_number =  bkind )
        rij = particle_set ( index_b ) % r - particle_set ( index_a ) % r
! get the energy
        rab2 = DOT_PRODUCT ( rij, rij )
        spline_data => onfo_list ( ibond ) % onfo_kind % spline_data
        CALL potential_s ( spline_data, rab2, e0, fscalar )
        pot = pot + e0

! get the forces
        DO idim = 1, 3
           delta(idim) = numerical_shift
           rab2_minus_delta = DOT_PRODUCT ( rij - delta, rij - delta )
           CALL potential_s ( spline_data, rab2_minus_delta, em, fscalar )
           rab2_plus_delta = DOT_PRODUCT ( rij + delta, rij + delta )
           CALL potential_s ( spline_data, rab2_plus_delta, ep, fscalar )
           f_numer(idim,index_a) = f_numer(idim,index_a) + (ep-em)
           f_numer(idim,index_b) = f_numer(idim,index_b) - (ep-em)
           delta = 0.0_dp
        END DO

     END DO BOND
     END DO
    END DO MOL

    f_numer = f_numer / ( 2.0_dp * numerical_shift )

END SUBROUTINE force_onefour_numer

!******************************************************************************
!!****** fist_force_numer/force_ub_numer [1.0] *
!!
!!   NAME
!!     force_bond_numer
!!
!!   FUNCTION
!!     Computes numerically the forces from the urey bradley
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE force_ub_numer ( molecule_set, molecule_kind_set, local_molecules, &
                            particle_set, f_ub_numer, pot_ub, numerical_shift )


    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: f_ub_numer
    REAL(KIND=dp), INTENT(INOUT)             :: pot_ub
    REAL(KIND=dp), INTENT(IN)                :: numerical_shift

    INTEGER                                  :: first_atom, i, idim, ikind, &
                                                imol, index_a, index_b, iub, &
                                                nkind, nmol_per_kind, nub
    REAL(KIND=dp)                            :: e0, em, ep, fscalar
    REAL(KIND=dp), DIMENSION(3)              :: delta, rij, rij_minus_delta, &
                                                rij_plus_delta
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule
    TYPE(ub_type), POINTER                   :: ub_list( : )

!------------------------------------------------------------------------------
! initializing the forces

  nkind = SIZE ( molecule_kind_set )
  delta = 0.0_dp
  MOL: DO ikind = 1, nkind
     nmol_per_kind = local_molecules % n_el ( ikind )
     DO imol = 1, nmol_per_kind
       i = local_molecules % list ( ikind ) % array ( imol )
       molecule => molecule_set ( i )
       molecule_kind => molecule % molecule_kind
       CALL get_molecule_kind ( molecule_kind, nub = nub,  &
                                ub_list = ub_list )
     BOND: DO iub = 1, nub
        CALL get_molecule ( molecule, first_atom=first_atom )
        index_a = ub_list ( iub ) % a + first_atom - 1
        index_b = ub_list ( iub ) % c + first_atom - 1
        rij = particle_set ( index_a ) % r - particle_set ( index_b ) % r
! get the energy
        CALL force_bonds ( rij, ub_list ( iub ) % ub_kind % r0, &
             ub_list ( iub ) % ub_kind  % k, e0, fscalar )
        pot_ub = pot_ub + e0

! get the forces
        DO idim = 1, 3
           delta(idim) = numerical_shift
           rij_minus_delta = rij - delta
           CALL force_bonds ( rij_minus_delta, ub_list ( iub ) % ub_kind % r0, &
             ub_list ( iub ) % ub_kind  % k, em, fscalar )
           rij_plus_delta = rij + delta
           CALL force_bonds ( rij_plus_delta, ub_list ( iub ) % ub_kind % r0, &
                ub_list ( iub ) % ub_kind % k, ep, fscalar )
           f_ub_numer(idim,index_a) = f_ub_numer(idim,index_a) + (em-ep)
           f_ub_numer(idim,index_b) = f_ub_numer(idim,index_b) - (em-ep)
           delta = 0.0_dp
        END DO

     END DO BOND
     END DO
    END DO MOL

    f_ub_numer = f_ub_numer / ( 2.0_dp * numerical_shift )

END SUBROUTINE force_ub_numer

!!*****
!******************************************************************************
!!****** fist_force_numer/force_bend_numer [1.0] *
!!
!!   NAME
!!     force_bend_numer
!!
!!   FUNCTION
!!     Computes numerically the forces from the bends
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE force_bend_numer ( molecule_set, molecule_kind_set, local_molecules, &
                              particle_set, f_bend_numer, pot_bend, numerical_shift )


    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: f_bend_numer
    REAL(KIND=dp), INTENT(INOUT)             :: pot_bend
    REAL(KIND=dp), INTENT(IN)                :: numerical_shift

    INTEGER                                  :: first_atom, i, ibend, idim, &
                                                ikind, imol, index_a, &
                                                index_b, index_c, nbends, &
                                                nkind, nmol_per_kind
    REAL(KIND=dp)                            :: d12, d12m, d12p, d32, d32m, &
                                                d32p, dist, e0, em, ep, &
                                                fscalar, id12, id12m, id12p, &
                                                id32, id32m, id32p, theta
    REAL(KIND=dp), DIMENSION(3)              :: b12, b12_minus_delta, &
                                                b12_plus_delta, b32, &
                                                b32_minus_delta, &
                                                b32_plus_delta, delta, g1, &
                                                g2, g3
    TYPE(bend_type), POINTER                 :: bend_list( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule

!------------------------------------------------------------------------------
! initializing the forces

  nkind = SIZE ( molecule_kind_set )
  MOL: DO ikind = 1, nkind
     nmol_per_kind = local_molecules % n_el ( ikind )
     DO imol = 1, nmol_per_kind
       i = local_molecules % list ( ikind ) % array ( imol )
       molecule => molecule_set ( i )
       molecule_kind => molecule % molecule_kind
       CALL get_molecule_kind ( molecule_kind, nbend = nbends, &
                                bend_list = bend_list )
BEND: DO ibend = 1, nbends
        CALL get_molecule ( molecule, first_atom=first_atom )
        index_a = bend_list ( ibend ) % a + first_atom - 1
        index_b = bend_list ( ibend ) % b + first_atom - 1
        index_c = bend_list ( ibend ) % c + first_atom - 1
        b12 = particle_set ( index_a ) % r - particle_set ( index_b ) % r
        b32 = particle_set ( index_c ) % r - particle_set ( index_b ) % r
        d12 = SQRT(DOT_PRODUCT(b12,b12))
        id12 = 1.0_dp/d12
        d32 = SQRT(DOT_PRODUCT(b32,b32))
        id32 = 1.0_dp/d32
        dist = DOT_PRODUCT(b12,b32)
        theta = ACOS(dist*id12*id32)
! get the energy
        CALL force_bends ( b12, b32, d12, d32, id12, id32, dist, theta, &
                           bend_list ( ibend ) % bend_kind % theta0,  &
                           bend_list ( ibend ) % bend_kind % k, &
                           g1, g2, g3, e0, fscalar )
        pot_bend = pot_bend + e0

! get the force on atom1 and atom2
        delta = 0.0_dp
        DO idim = 1, 3
           delta(idim) = numerical_shift
           b12_minus_delta = b12 - delta
           d12m = SQRT(DOT_PRODUCT(b12_minus_delta,b12_minus_delta))
           id12m = 1.0_dp / d12m
           dist = DOT_PRODUCT(b12_minus_delta,b32)
           theta = ACOS(dist*id12m*id32)
           CALL force_bends ( b12_minus_delta, b32, d12m, d32, id12m, id32, dist, theta, &
                              bend_list ( ibend ) % bend_kind % theta0,  &
                              bend_list ( ibend ) % bend_kind % k, &
                              g1, g2, g3, em, fscalar )
           b12_plus_delta = b12 + delta
           d12p = SQRT(DOT_PRODUCT(b12_plus_delta,b12_plus_delta))
           id12p = 1.0_dp/d12p
           dist = DOT_PRODUCT(b12_plus_delta,b32)
           theta = ACOS(dist*id12p*id32)
           CALL force_bends ( b12_plus_delta, b32, d12p, d32, id12p, id32, dist, theta, &
                              bend_list ( ibend ) % bend_kind % theta0,  &
                              bend_list ( ibend ) % bend_kind % k, &
                              g1, g2, g3, ep, fscalar )
           f_bend_numer(idim,index_a) = f_bend_numer(idim,index_a) + (em-ep)
           f_bend_numer(idim,index_b) = f_bend_numer(idim,index_b) - (em-ep)
           delta = 0.0_dp
        END DO

! get the force on atom3 and atom2
        delta = 0.0_dp
        DO idim = 1, 3
           delta(idim) = numerical_shift
           b32_minus_delta = b32 - delta
           d32m = SQRT(DOT_PRODUCT(b32_minus_delta,b32_minus_delta))
           id32m = 1.0_dp/d32m
           dist = DOT_PRODUCT(b32_minus_delta,b12)
           theta = ACOS(dist*id32m*id12)
           CALL force_bends(b12,b32_minus_delta,d12,d32m,id12,id32m,dist,theta, &
                            bend_list ( ibend ) % bend_kind % theta0,  &
                            bend_list ( ibend ) % bend_kind % k, &
                            g1, g2, g3, em, fscalar )
           b32_plus_delta = b32 + delta
           d32p = SQRT(DOT_PRODUCT(b32_plus_delta,b32_plus_delta))
           id32p = 1.0_dp/d32p
           dist = DOT_PRODUCT(b32_plus_delta,b12)
           theta = ACOS(dist*id32p*id12)
           CALL force_bends ( b12, b32_plus_delta, d12, d32p, id12, id32p, dist, theta, &
                              bend_list ( ibend ) % bend_kind % theta0,  &
                              bend_list ( ibend ) % bend_kind % k, &
                              g1, g2, g3, ep, fscalar )
           f_bend_numer(idim,index_c) = f_bend_numer(idim,index_c) + (em-ep)
           f_bend_numer(idim,index_b) = f_bend_numer(idim,index_b) - (em-ep)
           delta = 0.0_dp
        END DO

     END DO BEND
    END DO
  END DO MOL

     f_bend_numer = f_bend_numer / ( 2.0_dp * numerical_shift )

END SUBROUTINE force_bend_numer

!!*****
!******************************************************************************
!!****** fist_force_numer/force_torsion_numer [1.0] *
!!
!!   NAME
!!     force_torsion_numer
!!
!!   FUNCTION
!!     Computes numerically the forces from the torsions
!!
!!   AUTHOR
!!     DG
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE force_torsion_numer ( molecule_set, molecule_kind_set, local_molecules, &
                              particle_set, f_torsion_numer, pot_torsion, numerical_shift )


    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: f_torsion_numer
    REAL(KIND=dp), INTENT(INOUT)             :: pot_torsion
    REAL(KIND=dp), INTENT(IN)                :: numerical_shift

    INTEGER :: first_atom, i, idim, ikind, imol, imul, index_a, index_b, &
      index_c, index_d, itorsion, nkind, nmol_per_kind, ntorsions
    REAL(KIND=dp) :: dist1, dist1m, dist1p, dist2, dist2m, dist2p, e0, em, &
      ep, fscalar, is32, is32m, is32p, iw1, iw1m, iw1p, iw2, iw2m, iw2p, s32, &
      s32m, s32p, w1, w1m, w1p, w2, w2m, w2p
    REAL(KIND=dp), DIMENSION(3) :: delta, gt1, gt2, gt3, gt4, t12, &
      t12_minus_delta, t12_plus_delta, t32, t32_minus_delta, t32_plus_delta, &
      t43, t43_minus_delta, t43_plus_delta, u1, u1m, u1p, u2, u2m, u2p
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule
    TYPE(torsion_type), POINTER              :: torsion_list( : )

!------------------------------------------------------------------------------
! initializing the forces

  nkind = SIZE ( molecule_kind_set )
  MOL: DO ikind = 1, nkind
     nmol_per_kind = local_molecules % n_el ( ikind )
     DO imol = 1, nmol_per_kind
       i = local_molecules % list ( ikind ) % array ( imol )
       molecule => molecule_set ( i )
       molecule_kind => molecule % molecule_kind
       CALL get_molecule_kind ( molecule_kind, ntorsion = ntorsions, &
                                torsion_list = torsion_list )

     TORSION: DO itorsion = 1, ntorsions
        CALL get_molecule ( molecule, first_atom=first_atom )
        index_a = torsion_list ( itorsion ) % a + first_atom - 1
        index_b = torsion_list ( itorsion ) % b + first_atom - 1
        index_c = torsion_list ( itorsion ) % c + first_atom - 1
        index_d = torsion_list ( itorsion ) % d + first_atom - 1
        t12 = particle_set ( index_a ) % r - particle_set ( index_b ) % r
        t32 = particle_set ( index_c ) % r - particle_set ( index_b ) % r
        t43 = particle_set ( index_d ) % r - particle_set ( index_c ) % r
        s32 = SQRT(DOT_PRODUCT(t32,t32))
        is32 = 1.0_dp/s32
        dist1 = DOT_PRODUCT(t12,t32)
        dist2 = DOT_PRODUCT(t43,t32)
        u1 = t12 - dist1*t32*is32**2
        u2 = t43 - dist2*t32*is32**2
        w1 = SQRT(DOT_PRODUCT(u1,u1))
        w2 = SQRT(DOT_PRODUCT(u2,u2))
        iw1 = 1.0_dp/w1
        iw2 = 1.0_dp/w2
        ! get the energy
       DO imul = 1, torsion_list(itorsion) % torsion_kind % nmul
        CALL force_torsions(is32,dist1,dist2,u1,u2,iw1,iw2, &
                            torsion_list(itorsion)%torsion_kind%a, &
                            torsion_list(itorsion)%torsion_kind%k(imul), &
                            torsion_list(itorsion)%torsion_kind%cosphi0(imul), &
                            torsion_list(itorsion)%torsion_kind%m(imul), &
                            torsion_list(itorsion)%torsion_kind%name, &
                            gt1,gt2,gt3,gt4,e0,fscalar)
        pot_torsion = pot_torsion + e0

! get the force on atom1 and atom2
        delta = 0.0_dp
        DO idim = 1, 3
           delta(idim) = numerical_shift
           t12_minus_delta = t12 - delta
           s32 = SQRT(DOT_PRODUCT(t32,t32))
           is32 = 1.0_dp/s32
           dist1m = DOT_PRODUCT(t12_minus_delta,t32)
           dist2 = DOT_PRODUCT(t43,t32)
           u1m = t12_minus_delta - dist1m*t32*is32**2
           u2 = t43 - dist2*t32*is32**2
           w1m = SQRT(DOT_PRODUCT(u1m,u1m))
           w2 = SQRT(DOT_PRODUCT(u2,u2))
           iw1m = 1.0_dp/w1m
           iw2 = 1.0_dp/w2
           CALL force_torsions(is32,dist1m,dist2,u1m,u2,iw1m,iw2, &
                            torsion_list(itorsion)%torsion_kind%a, &
                            torsion_list(itorsion)%torsion_kind%k(imul), &
                            torsion_list(itorsion)%torsion_kind%cosphi0(imul), &
                            torsion_list(itorsion)%torsion_kind%m(imul), &
                            torsion_list(itorsion)%torsion_kind%name, &
                            gt1,gt2,gt3,gt4,em,fscalar)
           t12_plus_delta = t12 + delta
           s32 = SQRT(DOT_PRODUCT(t32,t32))
           is32 = 1.0_dp/s32
           dist1p = DOT_PRODUCT(t12_plus_delta,t32)
           dist2 = DOT_PRODUCT(t43,t32)
           u1p = t12_plus_delta - dist1p*t32*is32**2
           u2 = t43 - dist2*t32*is32**2
           w1p = SQRT(DOT_PRODUCT(u1p,u1p))
           w2 = SQRT(DOT_PRODUCT(u2,u2))
           iw1p = 1.0_dp/w1p
           iw2 = 1.0_dp/w2
           CALL force_torsions(is32,dist1p,dist2,u1p,u2,iw1p,iw2, &
                            torsion_list(itorsion)%torsion_kind%a, &
                            torsion_list(itorsion)%torsion_kind%k(imul), &
                            torsion_list(itorsion)%torsion_kind%cosphi0(imul), &
                            torsion_list(itorsion)%torsion_kind%m(imul), &
                            torsion_list(itorsion)%torsion_kind%name, &
                            gt1,gt2,gt3,gt4,ep,fscalar)
           f_torsion_numer(idim,index_a) = f_torsion_numer(idim,index_a) + (em-ep)
           f_torsion_numer(idim,index_b) = f_torsion_numer(idim,index_b) - (em-ep)
           delta = 0.0_dp
        END DO

! get the force on atom3 and atom2
        delta = 0.0_dp
        DO idim = 1, 3
           delta(idim) = numerical_shift
           t32_minus_delta = t32 - delta
           s32m = SQRT(DOT_PRODUCT(t32_minus_delta,t32_minus_delta))
           is32m = 1.0_dp/s32m
           dist1m = DOT_PRODUCT(t12,t32_minus_delta)
           dist2m = DOT_PRODUCT(t43,t32_minus_delta)
           u1m = t12 - dist1m*t32_minus_delta*is32m**2
           u2m = t43 - dist2m*t32_minus_delta*is32m**2
           w1m = SQRT(DOT_PRODUCT(u1m,u1m))
           w2m = SQRT(DOT_PRODUCT(u2m,u2m))
           iw1m = 1.0_dp/w1m
           iw2m = 1.0_dp/w2m
           CALL force_torsions(is32m,dist1m,dist2m,u1m,u2m,iw1m,iw2m, &
                            torsion_list(itorsion)%torsion_kind%a, &
                            torsion_list(itorsion)%torsion_kind%k(imul), &
                            torsion_list(itorsion)%torsion_kind%cosphi0(imul), &
                            torsion_list(itorsion)%torsion_kind%m(imul), &
                            torsion_list(itorsion)%torsion_kind%name, &
                             gt1,gt2,gt3,gt4,em,fscalar)
           t32_plus_delta = t32 + delta
           s32p = SQRT(DOT_PRODUCT(t32_plus_delta,t32_plus_delta))
           is32p = 1.0_dp/s32p
           dist1p = DOT_PRODUCT(t12,t32_plus_delta)
           dist2p = DOT_PRODUCT(t43,t32_plus_delta)
           u1p = t12 - dist1p*t32_plus_delta*is32p**2
           u2p = t43 - dist2p*t32_plus_delta*is32p**2
           w1p = SQRT(DOT_PRODUCT(u1p,u1p))
           w2p = SQRT(DOT_PRODUCT(u2p,u2p))
           iw1p = 1.0_dp/w1p
           iw2p = 1.0_dp/w2p
           CALL force_torsions(is32p,dist1p,dist2p,u1p,u2p,iw1p,iw2p, &
                            torsion_list(itorsion)%torsion_kind%a, &
                            torsion_list(itorsion)%torsion_kind%k(imul), &
                            torsion_list(itorsion)%torsion_kind%cosphi0(imul), &
                            torsion_list(itorsion)%torsion_kind%m(imul), &
                            torsion_list(itorsion)%torsion_kind%name, &
                            gt1,gt2,gt3,gt4,ep,fscalar)
           f_torsion_numer(idim,index_c) = f_torsion_numer(idim,index_c) + (em-ep)
           f_torsion_numer(idim,index_b) = f_torsion_numer(idim,index_b) - (em-ep)
           delta = 0.0_dp
        END DO

! get the force on atom4 and atom3
        delta = 0.0_dp
        DO idim = 1, 3
           delta(idim) = numerical_shift
           t43_minus_delta = t43 - delta
           s32 = SQRT(DOT_PRODUCT(t32,t32))
           is32 = 1.0_dp/s32
           dist1 = DOT_PRODUCT(t12,t32)
           dist2m = DOT_PRODUCT(t43_minus_delta,t32)
           u1 = t12 - dist1*t32*is32**2
           u2m = t43_minus_delta - dist2m*t32*is32**2
           w1 = SQRT(DOT_PRODUCT(u1,u1))
           w2m = SQRT(DOT_PRODUCT(u2m,u2m))
           iw1 = 1.0_dp/w1
           iw2m = 1.0_dp/w2m
           CALL force_torsions(is32,dist1,dist2m,u1,u2m,iw1,iw2m, &
                            torsion_list(itorsion)%torsion_kind%a, &
                            torsion_list(itorsion)%torsion_kind%k(imul), &
                            torsion_list(itorsion)%torsion_kind%cosphi0(imul), &
                            torsion_list(itorsion)%torsion_kind%m(imul), &
                            torsion_list(itorsion)%torsion_kind%name, &
                            gt1,gt2,gt3,gt4,em,fscalar)
           t43_plus_delta = t43 + delta
           s32 = SQRT(DOT_PRODUCT(t32,t32))
           is32 = 1.0_dp/s32
           dist1 = DOT_PRODUCT(t12,t32)
           dist2p = DOT_PRODUCT(t43_plus_delta,t32)
           u1 = t12 - dist1*t32*is32**2
           u2p = t43_plus_delta - dist2p*t32*is32**2
           w1 = SQRT(DOT_PRODUCT(u1,u1))
           w2p = SQRT(DOT_PRODUCT(u2p,u2p))
           iw1 = 1.0_dp/w1
           iw2p = 1.0_dp/w2p
           CALL force_torsions(is32,dist1,dist2p,u1,u2p,iw1,iw2p, &
                            torsion_list(itorsion)%torsion_kind%a, &
                            torsion_list(itorsion)%torsion_kind%k(imul), &
                            torsion_list(itorsion)%torsion_kind%cosphi0(imul), &
                            torsion_list(itorsion)%torsion_kind%m(imul), &
                            torsion_list(itorsion)%torsion_kind%name, &
                            gt1,gt2,gt3,gt4,ep,fscalar)
           f_torsion_numer(idim,index_d) = f_torsion_numer(idim,index_d) + (em-ep)
           f_torsion_numer(idim,index_c) = f_torsion_numer(idim,index_c) - (em-ep)
           delta = 0.0_dp
        END DO
       END DO

     END DO TORSION
   END DO
  END DO MOL

  f_torsion_numer = f_torsion_numer / ( 2.0_dp * numerical_shift )

END SUBROUTINE force_torsion_numer
!!*****
!******************************************************************************
!!****** fist_force_numer/force_imptors_numer [1.0] *
!!
!!   NAME
!!     force_imptors_numer
!!
!!   FUNCTION
!!     Computes numerically the forces from the torsions
!!
!!   AUTHOR
!!     DG
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE force_imptors_numer ( molecule_set, molecule_kind_set, local_molecules, &
                                 particle_set, f_numer, pot, numerical_shift )


    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: f_numer
    REAL(KIND=dp), INTENT(INOUT)             :: pot
    REAL(KIND=dp), INTENT(IN)                :: numerical_shift

    INTEGER :: first_atom, i, idim, ikind, imol, index_a, index_b, index_c, &
      index_d, itorsion, nimpr, nkind, nmol_per_kind
    REAL(KIND=dp) :: dist1, dist1m, dist1p, dist2, dist2m, dist2p, e0, em, &
      ep, fscalar, is32, is32_m, is32_p, ism, ism_m, ism_p, isn, isn_m, &
      isn_p, s32, s32_m, s32_p, sm, sm_m, sm_p, sn, sn_m, sn_p
    REAL(KIND=dp), DIMENSION(3) :: delta, gt1, gt2, gt3, gt4, t12, t12_m, &
      t12_p, t32, t32_m, t32_p, t34, t34_m, t34_p, tm, tm_m, tm_p, tn, tn_m, &
      tn_p
    TYPE(impr_type), POINTER                 :: impr_list( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule

!------------------------------------------------------------------------------

  nkind = SIZE ( molecule_kind_set )
  MOL: DO ikind = 1, nkind
     nmol_per_kind = local_molecules % n_el ( ikind )
     DO imol = 1, nmol_per_kind
       i = local_molecules % list ( ikind ) % array ( imol )
       molecule => molecule_set ( i )
       molecule_kind => molecule % molecule_kind
       CALL get_molecule_kind ( molecule_kind, nimpr = nimpr, &
                                impr_list = impr_list )

     TORSION: DO itorsion = 1, nimpr
        CALL get_molecule ( molecule, first_atom=first_atom )
        index_a = impr_list ( itorsion ) % a + first_atom - 1
        index_b = impr_list ( itorsion ) % b + first_atom - 1
        index_c = impr_list ( itorsion ) % c + first_atom - 1
        index_d = impr_list ( itorsion ) % d + first_atom - 1
        t12 = particle_set ( index_a ) % r - particle_set ( index_b ) % r
        t32 = particle_set ( index_c ) % r - particle_set ( index_b ) % r
        t34 = particle_set ( index_c ) % r - particle_set ( index_d ) % r
! t12 x t32
        tm(1) = t12(2)*t32(3)-t32(2)*t12(3)
        tm(2) = -t12(1)*t32(3)+t32(1)*t12(3)
        tm(3) = t12(1)*t32(2)-t32(1)*t12(2)
! t32 x t34
        tn(1) = t32(2)*t34(3)-t34(2)*t32(3)
        tn(2) = -t32(1)*t34(3)+t34(1)*t32(3)
        tn(3) = t32(1)*t34(2)-t34(1)*t32(2)
        sm = SQRT ( DOT_PRODUCT(tm,tm) )
        ism = 1.0_dp/sm
        sn = SQRT ( DOT_PRODUCT(tn,tn) )
        isn = 1.0_dp/sn
        s32 = SQRT ( DOT_PRODUCT(t32,t32) )
        is32 = 1.0_dp/s32
        dist1 = DOT_PRODUCT(t12,t32)
        dist2 = DOT_PRODUCT(t34,t32)
        ! get the energy
        CALL force_imp_torsions(s32, is32,ism,isn,dist1,dist2,tm,tn, t12, t32, &
                            impr_list ( itorsion ) % impr_kind % k, &
                            impr_list ( itorsion ) % impr_kind % phi0, &
                            gt1,gt2,gt3,gt4,e0,fscalar)
        pot = pot + e0

! get the force on atom1 and atom2
        delta = 0.0_dp
        DO idim = 1, 3
           delta(idim) = numerical_shift
           t12_m = t12 - delta
! t12 x t32
           tm_m(1) = t12_m(2)*t32(3)-t32(2)*t12_m(3)
           tm_m(2) = -t12_m(1)*t32(3)+t32(1)*t12_m(3)
           tm_m(3) = t12_m(1)*t32(2)-t32(1)*t12_m(2)
           sm_m = SQRT ( DOT_PRODUCT(tm_m,tm_m) )
           ism_m = 1.0_dp/sm_m
           dist1m = DOT_PRODUCT(t12_m,t32)
           CALL force_imp_torsions(s32, is32,ism_m,isn,dist1m,dist2,tm_m,tn,t12_m,t32, &
                                   impr_list ( itorsion ) % impr_kind % k, &
                                   impr_list ( itorsion ) % impr_kind % phi0, &
                                   gt1,gt2,gt3,gt4,em,fscalar)
           t12_p = t12 + delta
           tm_p(1) = t12_p(2)*t32(3)-t32(2)*t12_p(3)
           tm_p(2) = -t12_p(1)*t32(3)+t32(1)*t12_p(3)
           tm_p(3) = t12_p(1)*t32(2)-t32(1)*t12_p(2)
           sm_p = SQRT ( DOT_PRODUCT(tm_p,tm_p) )
           ism_p = 1.0_dp/sm_p
           dist1p = DOT_PRODUCT(t12_p,t32)
           CALL force_imp_torsions(s32, is32,ism_p,isn,dist1p,dist2,tm_p,tn,t12_p,t32, &
                                   impr_list ( itorsion ) % impr_kind % k, &
                                   impr_list ( itorsion ) % impr_kind % phi0, &
                                   gt1,gt2,gt3,gt4,ep,fscalar)

           f_numer(idim,index_a) = f_numer(idim,index_a) + (em-ep)
           f_numer(idim,index_b) = f_numer(idim,index_b) - (em-ep)
           delta = 0.0_dp
        END DO

! get the force on atom3 and atom2
        delta = 0.0_dp
        DO idim = 1, 3
           delta(idim) = numerical_shift
           t32_m = t32 - delta
! t12 x t32
           tm_m(1) = t12(2)*t32_m(3)-t32_m(2)*t12(3)
           tm_m(2) = -t12(1)*t32_m(3)+t32_m(1)*t12(3)
           tm_m(3) = t12(1)*t32_m(2)-t32_m(1)*t12(2)
! t32 x t34
           tn_m(1) = t32_m(2)*t34(3)-t34(2)*t32_m(3)
           tn_m(2) = -t32_m(1)*t34(3)+t34(1)*t32_m(3)
           tn_m(3) = t32_m(1)*t34(2)-t34(1)*t32_m(2)
           sn_m = SQRT ( DOT_PRODUCT(tn_m,tn_m) )
           isn_m= 1.0_dp/sn_m
           s32_m = SQRT ( DOT_PRODUCT(t32_m,t32_m) )
           is32_m = 1.0_dp/s32_m
           sm_m = SQRT ( DOT_PRODUCT(tm_m,tm_m) )
           ism_m = 1.0_dp/sm_m
           dist1m = DOT_PRODUCT(t12,t32_m)
           dist2m = DOT_PRODUCT(t34,t32_m)
           CALL force_imp_torsions(s32_m, is32_m,ism_m,isn_m,dist1m,dist2m,tm_m,tn_m, t12, t32_m, &
                                   impr_list ( itorsion ) % impr_kind % k, &
                                   impr_list ( itorsion ) % impr_kind % phi0, &
                                   gt1,gt2,gt3,gt4,em,fscalar)
           t32_p = t32 + delta
! t12 x t32
           tm_p(1) = t12(2)*t32_p(3)-t32_p(2)*t12(3)
           tm_p(2) = -t12(1)*t32_p(3)+t32_p(1)*t12(3)
           tm_p(3) = t12(1)*t32_p(2)-t32_p(1)*t12(2)
! t32 x t34
           tn_p(1) = t32_p(2)*t34(3)-t34(2)*t32_p(3)
           tn_p(2) = -t32_p(1)*t34(3)+t34(1)*t32_p(3)
           tn_p(3) = t32_p(1)*t34(2)-t34(1)*t32_p(2)
           sn_p = SQRT ( DOT_PRODUCT(tn_p,tn_p) )
           isn_p= 1.0_dp/sn_p
           s32_p = SQRT ( DOT_PRODUCT(t32_p,t32_p) )
           is32_p = 1.0_dp/s32_p
           sm_p = SQRT ( DOT_PRODUCT(tm_p,tm_p) )
           ism_p = 1.0_dp/sm_p
           dist1p = DOT_PRODUCT(t12,t32_p)
           dist2p = DOT_PRODUCT(t34,t32_p)
           CALL force_imp_torsions(s32_p, is32_p,ism_p,isn_p,dist1p,dist2p,tm_p,tn_p, t12, t32_p, &
                                   impr_list ( itorsion ) % impr_kind % k, &
                                   impr_list ( itorsion ) % impr_kind % phi0, &
                                   gt1,gt2,gt3,gt4,ep,fscalar)
           f_numer(idim,index_c) = f_numer(idim,index_c) + (em-ep)
           f_numer(idim,index_b) = f_numer(idim,index_b) - (em-ep)
           delta = 0.0_dp
        END DO

! get the force on atom4 and atom3
        delta = 0.0_dp
        DO idim = 1, 3
           delta(idim) = numerical_shift
           t34_m = t34 - delta
! t32 x t34
           tn_m(1) = t32(2)*t34_m(3)-t34_m(2)*t32(3)
           tn_m(2) = -t32(1)*t34_m(3)+t34_m(1)*t32(3)
           tn_m(3) = t32(1)*t34_m(2)-t34_m(1)*t32(2)
           sn_m = SQRT ( DOT_PRODUCT(tn_m,tn_m) )
           isn_m= 1.0_dp/sn_m
           dist2m = DOT_PRODUCT(t34_m,t32)
           CALL force_imp_torsions(s32, is32,ism,isn_m,dist1,dist2m,tm,tn_m,t12,t32, &
                                   impr_list ( itorsion ) % impr_kind % k, &
                                   impr_list ( itorsion ) % impr_kind % phi0, &
                                   gt1,gt2,gt3,gt4,em,fscalar)
           t34_p = t34 + delta
! t32 x t34
           tn_p(1) = t32(2)*t34_p(3)-t34_p(2)*t32(3)
           tn_p(2) = -t32(1)*t34_p(3)+t34_p(1)*t32(3)
           tn_p(3) = t32(1)*t34_p(2)-t34_p(1)*t32(2)
           sn_p = SQRT ( DOT_PRODUCT(tn_p,tn_p) )
           isn_p= 1.0_dp/sn_p
           dist2p = DOT_PRODUCT(t34_p,t32)
           CALL force_imp_torsions(s32, is32,ism,isn_p,dist1,dist2p,tm,tn_p,t12,t32, &
                                   impr_list ( itorsion ) % impr_kind % k, &
                                   impr_list ( itorsion ) % impr_kind % phi0, &
                                   gt1,gt2,gt3,gt4,ep,fscalar)
           f_numer(idim,index_c) = f_numer(idim,index_c) + (em-ep)
           f_numer(idim,index_d) = f_numer(idim,index_d) - (em-ep)
           delta = 0.0_dp
        END DO

     END DO TORSION
   END DO
  END DO MOL

  f_numer = f_numer / ( 2.0_dp * numerical_shift )

END SUBROUTINE force_imptors_numer

!!*****
!******************************************************************************
!!****** fist_force_numer/ewald_force_numer [1.0] *
!!
!!   NAME
!!     ewald_force_numer
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

SUBROUTINE ewald_force_numer ( ewald_env, ewald_pw, atomic_kind_set, particle_set,  &
                               local_particles, numerical_shift, f_numer )


    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_set( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_particles
    REAL(KIND=dp), INTENT(IN)                :: numerical_shift
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: f_numer

    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:)                           :: sum_igr
    INTEGER :: group, idim, iparticle, iparticle_kind, iparticle_local, isos, &
      ngtot, nparticle_kind, nparticle_local
    REAL(KIND=dp)                            :: alpha, em, ep, epsilon0, q
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: gauss, r_delta
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: g
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cell_type), POINTER                 :: cell
    TYPE(pw_grid_type), POINTER              :: pw_grid
    TYPE(pw_pool_type), POINTER              :: pw_pool

!------------------------------------------------------------------------------

  nparticle_kind = SIZE ( atomic_kind_set )
  CALL ewald_pw_get ( ewald_pw, pw_big_pool = pw_pool, cell = cell )
  pw_grid => pw_pool % pw_grid
  g => pw_grid % g
  ngtot = pw_grid % ngpts_cut

! allocating
  ALLOCATE (sum_igr(ngtot),STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'ewald_force_numer', 'sum_igr', ngtot )
  ALLOCATE (gauss(ngtot),STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'ewald_force_numer', 'gauss', ngtot )
  ALLOCATE (r_delta(3),STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'ewald_force_numer', 'r_delta', 3 )

! defining alpha and epsilon0
  CALL ewald_env_get ( ewald_env, alpha = alpha, eps0 = epsilon0, group = group )

! first initialize the arrays gauss and sum_igr
  CALL structure_factor ( atomic_kind_set, particle_set, local_particles, &
                          sum_igr, gauss, alpha, g, ngtot, group )

! computing the numerical force on each atom
  DO iparticle_kind = 1, nparticle_kind
    atomic_kind => atomic_kind_set(iparticle_kind)
    CALL get_atomic_kind(atomic_kind=atomic_kind, qeff=q)
    nparticle_local = local_particles%n_el(iparticle_kind)
    DO iparticle_local=1,nparticle_local
       iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
       r_delta ( : ) = particle_set ( iparticle ) % r ( : )
       DO idim = 1, 3
          r_delta(idim) = r_delta(idim) + numerical_shift
          CALL force_g_numer(ep, sum_igr, gauss, particle_set ( iparticle ) % r ,  &
                             r_delta, q, g, ngtot)
          r_delta(idim) = r_delta(idim) - 2.0_dp*numerical_shift
          CALL force_g_numer(em, sum_igr, gauss, particle_set ( iparticle ) % r ,  &
                             r_delta, q, g, ngtot)
          f_numer(idim,iparticle) = f_numer ( idim, iparticle ) + em - ep  
          r_delta(idim) = r_delta(idim) + numerical_shift
       END DO
    END DO
  END DO
  f_numer = f_numer /  2.0_dp / numerical_shift / epsilon0 / pw_grid % vol

  DEALLOCATE (r_delta,STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'ewald_force_numer','r_delta')
  DEALLOCATE (sum_igr,STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'ewald_force_numer','sum_igr')
  DEALLOCATE (gauss,STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'ewald_force_numer', 'gauss' )

END SUBROUTINE ewald_force_numer

!!*****
!******************************************************************************
!!****** fist_force_numer/pv_bond_numer [1.0] *
!!
!!   NAME
!!     pv_bond_numer
!!
!!   FUNCTION
!!     Computes numerically the pressure tensor from the bonds
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pv_bond_numer ( molecule_set, molecule_kind_set, local_molecules, &
                           particle_set, cell, pv_test, delta )


    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: pv_test
    REAL(KIND=dp), INTENT(IN)                :: delta

    INTEGER                                  :: i, ii, isos, j, jj, k, natoms
    REAL(KIND=dp)                            :: idelta, vm, vp
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: s
    REAL(KIND=dp), DIMENSION(3, 3)           :: dvdh, hmat
    TYPE(cell_type), POINTER                 :: cell_local

!------------------------------------------------------------------------------

  NULLIFY ( cell_local )
 
  CALL cell_create ( cell_local, hmat = cell % hmat, periodic = cell % perd )
  CALL cell_retain ( cell_local )

! assigning the local variables
  natoms = SIZE(particle_set)
  ALLOCATE (s(3,natoms),STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'pv_bond_numer', 's', 3 * natoms )
  DO i = 1, natoms
     s(:,i) = MATMUL ( cell_local % h_inv, particle_set ( i ) % r ( : ) )
  END DO

! Defining the increments
  idelta = 1.0_dp/(2.0_dp*delta)
  dvdh = 0.0_dp
  LOOP1: DO ii = 1, 3
     LOOP2: DO jj = 1, 3
        cell_local % hmat(ii,jj) = cell_local % hmat(ii,jj) + delta

! get inverse
        CALL get_hinv(cell_local)

! get coords from perturbed hmat
        DO i = 1, natoms
           particle_set ( i ) % r ( : ) = MATMUL(cell_local % hmat,s(:,i))
        END DO
        CALL potential_bond ( molecule_set, molecule_kind_set, local_molecules, &
                              particle_set, vp )

! tweak again
        cell_local % hmat(ii,jj) = cell_local % hmat(ii,jj) - 2.0_dp*delta

! get inverse
        CALL get_hinv ( cell_local )

! get coords from perturbed hmat
        DO i = 1, natoms
           particle_set ( i ) % r ( : ) = MATMUL(cell_local % hmat,s(:,i))
        END DO
        CALL potential_bond ( molecule_set, molecule_kind_set, local_molecules, &
                              particle_set, vm )

! calculate the derivative
        dvdh(ii,jj) = (vm-vp)*idelta
! resetting the hmat
        cell_local % hmat(ii,jj) = cell_local % hmat(ii,jj) + delta
        CALL get_hinv(cell_local)
! resetting particle_set
        DO i = 1, natoms
           particle_set ( i ) % r = MATMUL(cell_local % hmat,s(:,i))
        END DO
     END DO LOOP2
  END DO LOOP1

! computing pv_test
  DO i = 1, 3
     DO j = 1, 3
        DO k = 1, 3
           pv_test(i,j) = pv_test(i,j) + dvdh(i,k)*cell_local % hmat(j,k)
        END DO
     END DO
  END DO
  
  WRITE ( *, * ) 'PVTEST BOND'
  DO i =1,3
   WRITE ( *, * ) ( pv_test  ( i,j ), j = 1,3 ) 
  ENDDO

  DEALLOCATE (s,STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'pv_bond_numer', 's' )

  CALL cell_release ( cell_local )

END SUBROUTINE pv_bond_numer

!!*****
!******************************************************************************
!!****** fist_force_numer/pv_onefour_numer [1.0] *
!!
!!   NAME
!!     pv_onefour_numer
!!
!!   FUNCTION
!!     Computes numerically the pressure tensor from the onefours
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pv_onefour_numer ( molecule_set, molecule_kind_set, local_molecules, &
                           particle_set, cell, pv_test, delta )


    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: pv_test
    REAL(KIND=dp), INTENT(IN)                :: delta

    INTEGER                                  :: i, ii, isos, j, jj, k, natoms
    REAL(KIND=dp)                            :: idelta, vm, vp
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: s
    REAL(KIND=dp), DIMENSION(3, 3)           :: dvdh, hmat
    TYPE(cell_type), POINTER                 :: cell_local

!------------------------------------------------------------------------------

  NULLIFY ( cell_local )
 
  CALL cell_create ( cell_local, hmat = cell % hmat, periodic = cell % perd )
  CALL cell_retain ( cell_local )

! assigning the local variables
  natoms = SIZE(particle_set)
  ALLOCATE (s(3,natoms),STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'pv_bond_numer', 's', 3 * natoms )
  DO i = 1, natoms
     s(:,i) = MATMUL ( cell_local % h_inv, particle_set ( i ) % r ( : ) )
  END DO

! Defining the increments
  idelta = 1.0_dp/(2.0_dp*delta)
  dvdh = 0.0_dp
  LOOP1: DO ii = 1, 3
     LOOP2: DO jj = 1, 3
        cell_local % hmat(ii,jj) = cell_local % hmat(ii,jj) + delta

! get inverse
        CALL get_hinv(cell_local)

! get coords from perturbed hmat
        DO i = 1, natoms
           particle_set ( i ) % r ( : ) = MATMUL(cell_local % hmat,s(:,i))
        END DO
        CALL potential_onefour ( molecule_set, molecule_kind_set, local_molecules, &
                              particle_set, vp )

! tweak again
        cell_local % hmat(ii,jj) = cell_local % hmat(ii,jj) - 2.0_dp*delta

! get inverse
        CALL get_hinv ( cell_local )

! get coords from perturbed hmat
        DO i = 1, natoms
           particle_set ( i ) % r ( : ) = MATMUL(cell_local % hmat,s(:,i))
        END DO
        CALL potential_onefour ( molecule_set, molecule_kind_set, local_molecules, &
                                 particle_set, vm )

! calculate the derivative
        dvdh(ii,jj) = (vm-vp)*idelta
! resetting the hmat
        cell_local % hmat(ii,jj) = cell_local % hmat(ii,jj) + delta
        CALL get_hinv(cell_local)
! resetting particle_set
        DO i = 1, natoms
           particle_set ( i ) % r = MATMUL(cell_local % hmat,s(:,i))
        END DO
     END DO LOOP2
  END DO LOOP1

! computing pv_test
  DO i = 1, 3
     DO j = 1, 3
        DO k = 1, 3
           pv_test(i,j) = pv_test(i,j) + dvdh(i,k)*cell_local % hmat(j,k)
        END DO
     END DO
  END DO

  DEALLOCATE (s,STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'pv_onefour_numer', 's' )

  CALL cell_release ( cell_local )

END SUBROUTINE pv_onefour_numer

!!*****
!******************************************************************************
!!****** fist_force_numer/pv_ub_numer [1.0] *
!!
!!   NAME
!!     pv_ub_numer
!!
!!   FUNCTION
!!     Computes numerically the pressure tensor from the urey-bradley
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pv_ub_numer ( molecule_set, molecule_kind_set, local_molecules, &
                           particle_set, cell, pv_test, delta )


    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: pv_test
    REAL(KIND=dp), INTENT(IN)                :: delta

    INTEGER                                  :: i, ii, isos, j, jj, k, natoms
    REAL(KIND=dp)                            :: idelta, vm, vp
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: s
    REAL(KIND=dp), DIMENSION(3, 3)           :: dvdh, hmat
    TYPE(cell_type), POINTER                 :: cell_local

!------------------------------------------------------------------------------

  NULLIFY ( cell_local )
 
  CALL cell_create ( cell_local, hmat = cell % hmat, periodic = cell % perd )
  CALL cell_retain ( cell_local )

! assigning the local variables
  natoms = SIZE(particle_set)
  ALLOCATE (s(3,natoms),STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'pv_bond_numer', 's', 3 * natoms )
  DO i = 1, natoms
     s(:,i) = MATMUL ( cell_local % h_inv, particle_set ( i ) % r ( : ) )
  END DO

! Defining the increments
  idelta = 1.0_dp/(2.0_dp*delta)
  dvdh = 0.0_dp
  LOOP1: DO ii = 1, 3
     LOOP2: DO jj = 1, 3
        cell_local % hmat(ii,jj) = cell_local % hmat(ii,jj) + delta

! get inverse
        CALL get_hinv(cell_local)

! get coords from perturbed hmat
        DO i = 1, natoms
           particle_set ( i ) % r ( : ) = MATMUL(cell_local % hmat,s(:,i))
        END DO
        CALL potential_ub ( molecule_set, molecule_kind_set, local_molecules, &
                              particle_set, vp )

! tweak again
        cell_local % hmat(ii,jj) = cell_local % hmat(ii,jj) - 2.0_dp*delta

! get inverse
        CALL get_hinv ( cell_local )

! get coords from perturbed hmat
        DO i = 1, natoms
           particle_set ( i ) % r ( : ) = MATMUL(cell_local % hmat,s(:,i))
        END DO
        CALL potential_ub ( molecule_set, molecule_kind_set, local_molecules, &
                              particle_set, vm )

! calculate the derivative
        dvdh(ii,jj) = (vm-vp)*idelta
! resetting the hmat
        cell_local % hmat(ii,jj) = cell_local % hmat(ii,jj) + delta
        CALL get_hinv(cell_local)
! resetting particle_set
        DO i = 1, natoms
           particle_set ( i ) % r = MATMUL(cell_local % hmat,s(:,i))
        END DO
     END DO LOOP2
  END DO LOOP1

! computing pv_test
  DO i = 1, 3
     DO j = 1, 3
        DO k = 1, 3
           pv_test(i,j) = pv_test(i,j) + dvdh(i,k)*cell_local % hmat(j,k)
        END DO
     END DO
  END DO

  DEALLOCATE (s,STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'pv_ub_numer', 's' )

  CALL cell_release ( cell_local )

END SUBROUTINE pv_ub_numer

!!*****
!******************************************************************************
!!****** fist_force_numer/pv_bend_numer [1.0] *
!!
!!   NAME
!!     pv_bend_numer
!!
!!   FUNCTION
!!     Computes numerically the pressure tensor from the bends
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pv_bend_numer ( molecule_set, molecule_kind_set, local_molecules, &
                           particle_set, cell, pv_test, delta )


    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: pv_test
    REAL(KIND=dp), INTENT(IN)                :: delta

    INTEGER                                  :: i, ii, isos, j, jj, k, natoms
    REAL(KIND=dp)                            :: idelta, vm, vp
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: s
    REAL(KIND=dp), DIMENSION(3, 3)           :: dvdh
    TYPE(cell_type), POINTER                 :: cell_local

!------------------------------------------------------------------------------

  NULLIFY ( cell_local )
 
  CALL cell_create ( cell_local, hmat = cell % hmat, periodic = cell % perd )
  CALL cell_retain ( cell_local )

! assigning the local variables
  natoms = SIZE(particle_set)
  ALLOCATE (s(3,natoms),STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'pv_bend_numer', 's', 3 * natoms )
  DO i = 1, natoms
     s(:,i) = MATMUL(cell_local % h_inv,particle_set(i) % r)
  END DO


! Defining the increments
  idelta = 1.0_dp/(2.0_dp*delta)
  dvdh = 0.0_dp
  LOOP1: DO ii = 1, 3
     LOOP2: DO jj = 1, 3
        cell_local % hmat(ii,jj) = cell_local % hmat(ii,jj) + delta
! get inverse
        CALL get_hinv(cell_local)
! get coords from perturbed hmat
        DO i = 1, natoms
           particle_set ( i ) % r = MATMUL(cell_local % hmat,s(:,i))
        END DO
        CALL potential_bend ( molecule_set, molecule_kind_set, local_molecules, &
                              particle_set, vp )
! tweak again
        cell_local % hmat(ii,jj) = cell_local % hmat(ii,jj) - 2.0_dp*delta
! get inverse
        CALL get_hinv(cell_local)
! get coords from perturbed hmat
        DO i = 1, natoms
           particle_set ( i ) % r = MATMUL(cell_local % hmat,s(:,i))
        END DO
        CALL potential_bend ( molecule_set, molecule_kind_set, local_molecules, &
                              particle_set, vm )
! calculate the derivative
        dvdh(ii,jj) = (vm-vp)*idelta
! resetting the hmat
        cell_local % hmat(ii,jj) = cell_local % hmat(ii,jj) + delta
        CALL get_hinv(cell_local)
! resetting the hmat
        DO i = 1, natoms
           particle_set ( i ) % r = MATMUL(cell_local % hmat,s(:,i))
        END DO
     END DO LOOP2
  END DO LOOP1

! computing pv_test
  DO i = 1, 3
     DO j = 1, 3
        DO k = 1, 3
           pv_test(i,j) = pv_test(i,j) + dvdh(i,k)*cell_local % hmat(j,k)
        END DO
     END DO
  END DO

  DEALLOCATE (s,STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'pv_bend_numer', 's' )

  CALL cell_release ( cell_local )
  
END SUBROUTINE pv_bend_numer

!!*****
!******************************************************************************
!!****** fist_force_numer/pv_torsion_numer [1.0] *
!!
!!   NAME
!!     pv_torsion_numer
!!
!!   FUNCTION
!!     Computes numerically the pressure tensor from the torsions
!!
!!   AUTHOR
!!     DG
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pv_torsion_numer ( molecule_set, molecule_kind_set, local_molecules, &
                              particle_set, cell, pv_test, delta )


    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: pv_test
    REAL(KIND=dp), INTENT(IN)                :: delta

    INTEGER                                  :: i, ii, isos, j, jj, k, natoms
    REAL(KIND=dp)                            :: idelta, vm, vp
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: s
    REAL(KIND=dp), DIMENSION(3, 3)           :: dvdh
    TYPE(cell_type), POINTER                 :: cell_local

!------------------------------------------------------------------------------

  NULLIFY ( cell_local )
 
  CALL cell_create ( cell_local, hmat = cell % hmat, periodic = cell % perd )
  CALL cell_retain ( cell_local )

! assigning the local variables
  natoms = SIZE(particle_set)
  ALLOCATE (s(3,natoms),STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'pv_torsion_numer', 's', 3 * natoms )
  DO i = 1, natoms
     s(:,i) = MATMUL(cell_local % h_inv,particle_set(i)% r)
  END DO


! Defining the increments
  idelta = 1.0_dp/(2.0_dp*delta)
  dvdh = 0.0_dp
  LOOP1: DO ii = 1, 3
     LOOP2: DO jj = 1, 3
        cell_local % hmat(ii,jj) = cell_local % hmat(ii,jj) + delta
! get inverse
        CALL get_hinv(cell_local)
! get coords from perturbed hmat
        DO i = 1, natoms
           particle_set (i) % r = MATMUL(cell_local % hmat,s(:,i))
        END DO
        CALL potential_torsion ( molecule_set, molecule_kind_set, local_molecules, &
                              particle_set, vp )
! tweak again
        cell_local % hmat(ii,jj) = cell_local % hmat(ii,jj) - 2.0_dp*delta
! get inverse
        CALL get_hinv( cell_local )
! get coords from perturbed hmat
        DO i = 1, natoms
           particle_set(i)%r = MATMUL(cell_local % hmat,s(:,i))
        END DO
        CALL potential_torsion ( molecule_set, molecule_kind_set, local_molecules, &
                              particle_set, vm )
! calculate the derivative
        dvdh(ii,jj) = (vm-vp)*idelta
! resetting the hmat and particle_set
        cell_local % hmat(ii,jj) = cell_local % hmat(ii,jj) + delta
        CALL get_hinv(cell_local)
        DO i = 1, natoms
           particle_set(i)%r = MATMUL(cell_local % hmat,s(:,i))
        END DO
     END DO LOOP2
  END DO LOOP1

! computing pv_test
  DO i = 1, 3
     DO j = 1, 3
        DO k = 1, 3
           pv_test(i,j) = pv_test(i,j) + dvdh(i,k)*cell_local % hmat(j,k)
        END DO
     END DO
  END DO

  DEALLOCATE (s,STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'pv_torsion_numer', 's' )

  CALL cell_release ( cell_local )

END SUBROUTINE pv_torsion_numer
!!*****
!******************************************************************************
!!****** fist_force_numer/pv_imptors_numer [1.0] *
!!
!!   NAME
!!     pv_imptors_numer
!!
!!   FUNCTION
!!     Computes numerically the pressure tensor from the torsions
!!
!!   AUTHOR
!!     DG
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pv_imptors_numer ( molecule_set, molecule_kind_set, local_molecules, &
                              particle_set, cell, pv_test, delta )


    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: pv_test
    REAL(KIND=dp), INTENT(IN)                :: delta

    INTEGER                                  :: i, ii, isos, j, jj, k, natoms
    REAL(KIND=dp)                            :: idelta, vm, vp
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: s
    REAL(KIND=dp), DIMENSION(3, 3)           :: dvdh
    TYPE(cell_type), POINTER                 :: cell_local

!------------------------------------------------------------------------------

  NULLIFY ( cell_local )
 
  CALL cell_create ( cell_local, hmat = cell % hmat, periodic = cell % perd )
  CALL cell_retain ( cell_local )

! assigning the local variables
  natoms = SIZE(particle_set)
  ALLOCATE (s(3,natoms),STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'pv_imptors_numer', 's', 3 * natoms )
  DO i = 1, natoms
     s(:,i) = MATMUL(cell_local % h_inv,particle_set(i)% r)
  END DO


! Defining the increments
  idelta = 1.0_dp/(2.0_dp*delta)
  dvdh = 0.0_dp
  LOOP1: DO ii = 1, 3
     LOOP2: DO jj = 1, 3
        cell_local % hmat(ii,jj) = cell_local % hmat(ii,jj) + delta
! get inverse
        CALL get_hinv(cell_local)
! get coords from perturbed hmat
        DO i = 1, natoms
           particle_set (i) % r = MATMUL(cell_local % hmat,s(:,i))
        END DO
        CALL potential_imptors ( molecule_set, molecule_kind_set, local_molecules, &
                              particle_set, vp )
! tweak again
        cell_local % hmat(ii,jj) = cell_local % hmat(ii,jj) - 2.0_dp*delta
! get inverse
        CALL get_hinv( cell_local )
! get coords from perturbed hmat
        DO i = 1, natoms
           particle_set(i)%r = MATMUL(cell_local % hmat,s(:,i))
        END DO
        CALL potential_imptors ( molecule_set, molecule_kind_set, local_molecules, &
                                particle_set, vm )
! calculate the derivative
        dvdh(ii,jj) = (vm-vp)*idelta
! resetting the hmat and particle_set
        cell_local % hmat(ii,jj) = cell_local % hmat(ii,jj) + delta
        CALL get_hinv(cell_local)
        DO i = 1, natoms
           particle_set(i)%r = MATMUL(cell_local % hmat,s(:,i))
        END DO
     END DO LOOP2
  END DO LOOP1

! computing pv_test
  DO i = 1, 3
     DO j = 1, 3
        DO k = 1, 3
           pv_test(i,j) = pv_test(i,j) + dvdh(i,k)*cell_local % hmat(j,k)
        END DO
     END DO
  END DO

  DEALLOCATE (s,STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'pv_imptors_numer', 's' )

  CALL cell_release ( cell_local )

END SUBROUTINE pv_imptors_numer

!!*****
!******************************************************************************
!!****** fist_force_numer/pv_nonbond_numer [1.0] *
!!
!!   NAME
!!     ptens_numer
!!
!!   FUNCTION
!!     Computes the numerical pressure tensor
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pv_nonbond_numer ( fist_nonbond_env, ewald_env, atomic_kind_set, &
                              particle_set, cell, pv_test, delta )


    TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_set( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), DIMENSION(3, 3), &
      INTENT(INOUT)                          :: pv_test
    REAL(KIND=dp), INTENT(IN)                :: delta

    INTEGER                                  :: i, ii, isos, j, jj, k, natoms
    REAL(KIND=dp)                            :: idelta, vm, vnb, vp
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: s
    REAL(KIND=dp), DIMENSION(3, 3)           :: dvdh
    TYPE(cell_type), POINTER                 :: cell_local

!------------------------------------------------------------------------------

  NULLIFY ( cell_local )
 
  CALL cell_create ( cell_local, hmat = cell % hmat, periodic = cell % perd )
  CALL cell_retain ( cell_local )

! assigning the local variables
  natoms = SIZE(particle_set)
  ALLOCATE (s(3,natoms),STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'ptens_numer', 's', 3 * natoms )

  DO i = 1, natoms
     s(:,i) = MATMUL(cell_local % h_inv, particle_set ( i ) % r )
  END DO


! Defining the increments
  idelta = 1.0_dp / (2.0_dp * delta )
  dvdh = 0.0_dp
  LOOP1: DO ii = 1, 3
     LOOP2: DO jj = 1, 3
        cell_local % hmat(ii,jj) = cell_local % hmat(ii,jj) + delta

! get inverse
        CALL get_hinv ( cell_local )

! get coords from perturbed hmat
        DO i = 1, natoms
           particle_set ( i ) % r ( : ) = MATMUL ( cell_local % hmat, s(:,i) )
!           particle_set ( i ) % rho = 0._dp
!           particle_set ( i ) % f_embedding = 0._dp
        END DO
!        CALL density_nonbond_numer
!        CALL energy_manybody_numer
        CALL potential_nonbond ( fist_nonbond_env, particle_set, atomic_kind_set, &
                                 ewald_env, cell_local, vnb )
        vp = vnb

! tweak again
        cell_local % hmat(ii,jj) = cell_local % hmat(ii,jj) - 2.0_dp*delta
        CALL get_hinv(cell_local)
        DO i = 1, natoms
           particle_set (i) % r = MATMUL(cell_local % hmat,s(:,i))
!           particle_set ( i ) % rho = 0._dp
!           particle_set ( i ) % f_embedding = 0._dp
        END DO
!        CALL density_nonbond_numer
!        CALL energy_manybody_numer
        CALL potential_nonbond ( fist_nonbond_env, particle_set, atomic_kind_set, &
                                 ewald_env, cell_local, vnb )
        vm = vnb

! calculate the derivative
        dvdh(ii,jj) = (vm-vp)*idelta

! resetting the hmat
        cell_local% hmat(ii,jj) = cell_local % hmat(ii,jj) + delta
        CALL get_hinv(cell_local)
        DO i = 1, natoms
           particle_set (i) % r = MATMUL(cell_local % hmat,s(:,i))
        END DO
     END DO LOOP2
  END DO LOOP1

! computing pv_test
  DO i = 1, 3
     DO j = 1, 3
        DO k = 1, 3
           pv_test(i,j) = pv_test(i,j) + dvdh(i,k)*cell_local % hmat(j,k)
        END DO
     END DO
  END DO

  DEALLOCATE (s,STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory( 'ptens_numer', 's' )

  CALL cell_release ( cell_local )

END SUBROUTINE pv_nonbond_numer

!!*****
!******************************************************************************
!!****** fist_force_numer/pv_g_numer [1.0] *
!!
!!   NAME
!!     pv_g_numer
!!
!!   FUNCTION
!!     Computes the numerical pressure tensor
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE pv_g_numer ( ewald_env, ewald_pw, atomic_kind_set, particle_set,  &
                        local_particles, delta, pv_test )


    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_set( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_particles
    REAL(KIND=dp), INTENT(IN)                :: delta
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: pv_test

    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:)                           :: sum_igr
    INTEGER                                  :: group, i, ig, ii, isos, j, &
                                                jj, k, natoms, ngtot
    REAL(KIND=dp)                            :: alpha, epsilon0, idelta, vm, &
                                                vp
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: gauss
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: s
    REAL(KIND=dp), DIMENSION(3, 3)           :: dvdh
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: g
    TYPE(cell_type), POINTER                 :: cell, cell_local
    TYPE(pw_grid_type), POINTER              :: pw_grid
    TYPE(pw_pool_type), POINTER              :: pw_pool

!------------------------------------------------------------------------------

  CALL ewald_pw_get ( ewald_pw, pw_big_pool = pw_pool, cell = cell )
  pw_grid => pw_pool % pw_grid
  ngtot = pw_grid % ngpts_cut

  NULLIFY ( cell_local )
 
  CALL cell_create ( cell_local, hmat = cell % hmat, periodic = cell % perd )
  CALL cell_retain ( cell_local )

! allocating
  natoms = SIZE ( particle_set )
  ALLOCATE ( s ( 3, natoms ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'pvg_numer', 's', 3 * natoms )
  ALLOCATE (sum_igr(ngtot),STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'pvg_numer', 'sum_igr', ngtot )
  ALLOCATE (gauss(ngtot),STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'pvg_numer', 'gauss', ngtot )

! assigning the local variables
  DO i = 1, natoms
     s(:,i) = MATMUL ( cell_local % h_inv, particle_set ( i ) % r )
  END DO

! defining alpha and epsilon0 (see ewald.f for details)
  CALL ewald_env_get ( ewald_env, alpha = alpha, eps0 = epsilon0, group = group )

! Initializing pv_test
  pv_test = 0.0_dp

! Defining the increments
  idelta = 1.0_dp/(2.0_dp*delta)
  dvdh = 0.0_dp
  LOOP1: DO ii = 1, 3
     LOOP2: DO jj = 1, 3
        cell_local % hmat(ii,jj) = cell_local % hmat(ii,jj) + delta

! get inverse
        CALL get_hinv(cell_local)
        CALL pw_grid_change ( cell_local, pw_grid )
        g => pw_grid % g
        
! recompute new partical postions from old scaled coords.
        DO i = 1, natoms
           particle_set ( i ) % r = MATMUL ( cell_local % hmat, s(:,i) )
        END DO

! compute the g-vectors
        CALL structure_factor ( atomic_kind_set, particle_set, local_particles, &
                                sum_igr, gauss, alpha, g, ngtot, group )
        vp = 0.0_dp
        DO ig = 1, ngtot
          vp = vp + gauss(ig)* &
                   REAL ( sum_igr ( ig ) *  &
                   CONJG( sum_igr ( ig ) ),KIND=dp)
        END DO

        vp = vp /  epsilon0 / pw_grid % vol

! tweak again
        cell_local % hmat(ii,jj) = cell_local % hmat(ii,jj) - 2.0_dp * delta

! get inverse
        CALL get_hinv ( cell_local )
        CALL pw_grid_change ( cell_local, pw_grid )
        g => pw_grid % g

! recompute new partical postions from old scaled coords.
        DO i = 1, natoms
           particle_set ( i ) % r = MATMUL ( cell_local % hmat, s(:,i) )
        END DO

! compute the g-vectors

        CALL structure_factor ( atomic_kind_set, particle_set, local_particles, &
                                sum_igr, gauss, alpha, g, ngtot, group )
        vm = 0.0_dp
        DO ig = 1, ngtot
          vm = vm + gauss(ig)* &
                   REAL ( sum_igr ( ig ) *  &
                   CONJG( sum_igr ( ig ) ),KIND=dp)
        END DO

        vm = vm / epsilon0 / pw_grid % vol

! calculate the derivative
        dvdh ( ii, jj ) = ( vm - vp ) * idelta

! resetting the hmat
        cell_local % hmat(ii,jj) = cell_local % hmat(ii,jj) + delta
        CALL pw_grid_change ( cell_local, pw_grid )

! resetting particle_set
        DO i = 1, natoms
           particle_set ( i ) % r = MATMUL ( cell_local % hmat, s(:,i) )
        END DO
     END DO LOOP2
  END DO LOOP1

! computing pv_test
  DO i = 1, 3
     DO j = 1, 3
        DO k = 1, 3
           pv_test(i,j) = pv_test(i,j) + dvdh(i,k)*cell_local % hmat(j,k)
        END DO
     END DO
  END DO

  DEALLOCATE (s,STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory( 'pvg_numer', 's' )
  DEALLOCATE (sum_igr,STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory( 'pvg_numer', 'sum_igr' )
  DEALLOCATE (gauss,STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory( 'pvg_numer', 'gauss' )

  CALL cell_release ( cell_local )

END SUBROUTINE pv_g_numer

!!*****
!******************************************************************************
!!****** fist_force_numer/structure_factor [1.0] *
!!
!!   NAME
!!     potential_g_numer
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     CJM
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE structure_factor ( atomic_kind_set, particle_set, local_particles, &
                              sum_igr, gauss, alpha, glocal, ngtot, group)

    TYPE(atomic_kind_type), POINTER          :: atomic_kind_set( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_particles
    COMPLEX(KIND=dp), DIMENSION(:), &
      INTENT(OUT)                            :: sum_igr
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: gauss
    REAL(KIND=dp), INTENT(IN)                :: alpha
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: glocal
    INTEGER, INTENT(IN)                      :: ngtot, group

    INTEGER                                  :: ig, iparticle, &
                                                iparticle_kind, &
                                                iparticle_local, &
                                                nparticle_kind, &
                                                nparticle_local
    REAL(KIND=dp)                            :: gdotr, gsq, q
    TYPE(atomic_kind_type), POINTER          :: atomic_kind

!------------------------------------------------------------------------------

  nparticle_kind = SIZE ( atomic_kind_set )
! numerical energy
  sum_igr = CMPLX ( 0.0_dp, 0.0_dp,KIND=dp)
  DO ig = 1, ngtot
     gsq = DOT_PRODUCT(glocal(:,ig),glocal(:,ig))
     IF ( gsq <= 1.0E-10_dp ) CYCLE
     gauss(ig) = EXP(-gsq*0.25_dp/alpha/alpha)/gsq
     DO iparticle_kind = 1, nparticle_kind
       atomic_kind => atomic_kind_set(iparticle_kind)
       CALL get_atomic_kind(atomic_kind=atomic_kind, qeff=q)
       nparticle_local = local_particles%n_el(iparticle_kind)
       DO iparticle_local=1,nparticle_local
          iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
          gdotr = DOT_PRODUCT(particle_set ( iparticle ) % r ( : ), glocal(:,ig))
          sum_igr(ig) = sum_igr(ig) + q*CMPLX(COS(gdotr),SIN(gdotr),dp)
       END DO
     END DO

  END DO
  CALL mp_sum ( sum_igr, group )
  END SUBROUTINE structure_factor

!!*****
!******************************************************************************
!!****** fist_force_numer/force_g_numer [1.0] *
!!
!!   NAME
!!     de_g_numer
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

SUBROUTINE force_g_numer(energy,sum_igr,gauss,r,r_delta,charge,glocal, &
                         igtot)


    REAL(KIND=dp), INTENT(INOUT)             :: energy
    COMPLEX(KIND=dp), DIMENSION(:), &
      INTENT(IN)                             :: sum_igr
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: gauss, r, r_delta
    REAL(KIND=dp), INTENT(IN)                :: charge
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: glocal
    INTEGER, INTENT(IN)                      :: igtot

    COMPLEX(KIND=dp)                         :: sum
    INTEGER                                  :: ig
    REAL(KIND=dp)                            :: gdotr, gdotr_delta, gsq

!------------------------------------------------------------------------------

  sum = CMPLX( 0.0_dp, 0.0_dp,KIND=dp)
  energy = 0.0_dp
  DO ig = 1, igtot
     gsq = DOT_PRODUCT( glocal(:,ig), glocal(:,ig))
     IF ( gsq <= 1.0E-10_dp ) CYCLE
! compute g.r and g.(r+delta)
     gdotr = DOT_PRODUCT ( r ( : ), glocal ( :, ig ) )
     gdotr_delta = DOT_PRODUCT ( r_delta ( : ), glocal ( :, ig ) )
! subtract off exp(ig.r) and add exp(ig.(r+delta))
     sum = sum_igr ( ig ) - charge * &
           CMPLX ( COS( gdotr ), SIN( gdotr ),KIND=dp) + &
           charge * CMPLX( COS( gdotr_delta ), SIN( gdotr_delta ),KIND=dp)
! recompute energy
     energy = energy + gauss(ig) * REAL ( sum * CONJG ( sum ),KIND=dp)
  END DO

END SUBROUTINE force_g_numer

!!*****
!******************************************************************************
!!****** fist_force_numer/pme_force_numer [1.0] *
!!
!!   NAME
!!     pme_force_numer
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

SUBROUTINE pme_force_numer ( ewald_env, ewald_pw, particle_set, delta, f_part ,&
     force_env_section, error)
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    REAL(KIND=dp), INTENT(IN)                :: delta
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: f_part
    TYPE(section_vals_type), POINTER         :: force_env_section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    INTEGER                                  :: i, id, natoms
    REAL(KIND=dp)                            :: eminus, eplus

  natoms = SIZE ( particle_set ) 
  DO i = 1, natoms 
    DO id = 1, 3
      particle_set ( i ) % r ( id ) = particle_set ( i ) % r ( id ) + delta 
      CALL pme_evaluate ( ewald_env, ewald_pw, particle_set, eplus,&
           force_env_section=force_env_section,error=error)
      particle_set ( i ) % r ( id ) = particle_set ( i ) % r ( id ) - 2.0_dp * delta 
      CALL pme_evaluate ( ewald_env, ewald_pw, particle_set, eminus,&
           force_env_section=force_env_section,error=error )
      f_part ( id, i ) = eminus - eplus
      particle_set ( i ) % r ( id ) = particle_set ( i ) % r ( id ) + delta 
    ENDDO
  ENDDO
  f_part = f_part / 2.0_dp / delta
   
END SUBROUTINE pme_force_numer

!!*****
!******************************************************************************
!!****** fist_force_numer/ewald_energy_numer [1.0] *
!!
!!   NAME
!!     ewald_energy_numer
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
 SUBROUTINE ewald_energy_numer ( ewald_env, ewald_pw, atomic_kind_set, &
                                 particle_set, local_particles, energy_numer )


    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(atomic_kind_type), POINTER          :: atomic_kind_set( : )
    TYPE(particle_type), POINTER             :: particle_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_particles
    REAL(KIND=dp), INTENT(INOUT)             :: energy_numer

    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:)                           :: sum_igr
    INTEGER                                  :: group, ig, isos, ngtot
    REAL(KIND=dp)                            :: alpha, epsilon0
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: gauss
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: g
    TYPE(cell_type), POINTER                 :: cell
    TYPE(pw_grid_type), POINTER              :: pw_grid
    TYPE(pw_pool_type), POINTER              :: pw_pool

!------------------------------------------------------------------------------

  CALL ewald_pw_get ( ewald_pw, pw_big_pool = pw_pool, cell = cell )
  pw_grid => pw_pool % pw_grid
  g => pw_grid % g
  ngtot = pw_grid % ngpts_cut

! allocating
  ALLOCATE ( sum_igr ( ngtot ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'ewald_energy_numer', 'sum_igr', ngtot )
  ALLOCATE ( gauss ( ngtot ), STAT = isos )
  IF ( isos /= 0 ) CALL stop_memory ( 'ewald_energy_numer', 'gauss', ngtot )

! defining alpha and epsilon0
  CALL ewald_env_get ( ewald_env, alpha = alpha, eps0 = epsilon0, group = group )

! computing the energy
  CALL structure_factor ( atomic_kind_set, particle_set, local_particles, &
                          sum_igr, gauss, alpha, g, ngtot, group )

  
! computing the potential energy
  energy_numer = 0.0_dp
  DO ig = 1, ngtot
     energy_numer = energy_numer + gauss(ig)* REAL ( sum_igr ( ig ) *  &
                                  CONJG( sum_igr ( ig ) ),KIND=dp)
  END DO


  energy_numer = energy_numer / epsilon0 / pw_grid % vol 
  DEALLOCATE (sum_igr,STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'ewald_energy_numer', 'sum_igr' )
  DEALLOCATE (gauss,STAT=isos)
  IF ( isos /= 0 ) CALL stop_memory ( 'ewald_energy_numer', 'gauss' )

END SUBROUTINE ewald_energy_numer

!!*****
!******************************************************************************
!!****** fist_force_numer/potential_bond [1.0] *
!!
!!   NAME
!!     potential_bond
!!
!!   FUNCTION
!!     Computes analytically the potential from the bonds
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE potential_bond ( molecule_set, molecule_kind_set, local_molecules, &
                           particle_set, energy )


    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: energy

    INTEGER                                  :: first_atom, i, ibond, ikind, &
                                                imol, index_a, index_b, &
                                                nbonds, nkind, nmol_per_kind
    REAL(KIND=dp)                            :: e0, fscalar
    REAL(KIND=dp), DIMENSION(3)              :: rij
    TYPE(bond_type), POINTER                 :: bond_list( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule

!------------------------------------------------------------------------------

  nkind = SIZE ( molecule_kind_set )
  energy = 0._dp
  e0 = 0.0_dp
  fscalar = 0.0_dp
  MOL: DO ikind = 1, nkind
     nmol_per_kind = local_molecules % n_el ( ikind )
     DO imol = 1, nmol_per_kind
       i = local_molecules % list ( ikind ) % array ( imol )
       molecule => molecule_set ( i )
       molecule_kind => molecule % molecule_kind
       CALL get_molecule_kind ( molecule_kind, nbond = nbonds, &
                                bond_list = bond_list )
     BOND: DO ibond = 1, nbonds
        CALL get_molecule ( molecule, first_atom=first_atom )
        index_a = bond_list ( ibond ) % a + first_atom - 1
        index_b = bond_list ( ibond ) % b + first_atom - 1
        rij = particle_set ( index_a ) % r - particle_set ( index_b ) % r
        CALL force_bonds(rij,bond_list ( ibond ) % bond_kind % r0, &
                         bond_list ( ibond ) % bond_kind % k, &
                         e0,fscalar)
        energy = energy + e0
     END DO BOND
    END DO
  END DO MOL

END SUBROUTINE potential_bond
!!*****
!******************************************************************************
!!****** fist_force_numer/potential_onefour [1.0] *
!!
!!   NAME
!!     potential_bond
!!
!!   FUNCTION
!!     Computes analytically the potential from the onefours
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE potential_onefour ( molecule_set, molecule_kind_set, local_molecules, &
                               particle_set, energy )


    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: energy

    INTEGER                                  :: akind, bkind, first_atom, i, &
                                                ibond, ikind, imol, index_a, &
                                                index_b, nkind, &
                                                nmol_per_kind, nonfo
    REAL(KIND=dp)                            :: e0, fscalar, rab2
    REAL(KIND=dp), DIMENSION(3)              :: rij
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule
    TYPE(onfo_type), POINTER                 :: onfo_list( : )
    TYPE(spline_data_type), POINTER          :: spline_data

!------------------------------------------------------------------------------

  nkind = SIZE ( molecule_kind_set )
  energy = 0.0_dp
  e0 = 0.0_dp
  fscalar = 0.0_dp
  MOL: DO ikind = 1, nkind
     nmol_per_kind = local_molecules % n_el ( ikind )
     DO imol = 1, nmol_per_kind
       i = local_molecules % list ( ikind ) % array ( imol )
       molecule => molecule_set ( i )
       molecule_kind => molecule % molecule_kind
       CALL get_molecule_kind ( molecule_kind, nonfo = nonfo, &
                                onfo_list = onfo_list )
     BOND: DO ibond = 1, nonfo
        CALL get_molecule ( molecule, first_atom=first_atom )
        index_a = onfo_list ( ibond ) % a + first_atom - 1
        atomic_kind => particle_set ( index_a ) % atomic_kind
        CALL get_atomic_kind ( atomic_kind, kind_number =  akind )
        index_b = onfo_list ( ibond ) % b + first_atom - 1
        atomic_kind => particle_set ( index_b ) % atomic_kind
        CALL get_atomic_kind ( atomic_kind, kind_number =  bkind )
        rij = particle_set ( index_b ) % r - particle_set ( index_a ) % r
        rab2 = DOT_PRODUCT ( rij, rij )
        spline_data => onfo_list ( ibond ) % onfo_kind % spline_data
        CALL potential_s ( spline_data, rab2, e0, fscalar )
        energy = energy + e0
     END DO BOND
    END DO
  END DO MOL

END SUBROUTINE potential_onefour

!!*****
!******************************************************************************
!!****** fist_force_numer/potential_ub [1.0] *
!!
!!   NAME
!!     potential_ub
!!
!!   FUNCTION
!!     Computes analytically the potential from the urey-bradley
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE potential_ub ( molecule_set, molecule_kind_set, local_molecules, &
                           particle_set, energy )


    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: energy

    INTEGER                                  :: first_atom, i, ikind, imol, &
                                                index_a, index_b, iub, nkind, &
                                                nmol_per_kind, nub
    REAL(KIND=dp)                            :: e0, fscalar
    REAL(KIND=dp), DIMENSION(3)              :: rij
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule
    TYPE(ub_type), POINTER                   :: ub_list( : )

!------------------------------------------------------------------------------

  nkind = SIZE ( molecule_kind_set )
  energy = 0.0_dp
  e0 = 0.0_dp
  fscalar = 0.0_dp
  MOL: DO ikind = 1, nkind
     nmol_per_kind = local_molecules % n_el ( ikind )
     DO imol = 1, nmol_per_kind
       i = local_molecules % list ( ikind ) % array ( imol )
       molecule => molecule_set ( i )
       molecule_kind => molecule % molecule_kind
       CALL get_molecule_kind ( molecule_kind, nub = nub, &
                                ub_list = ub_list )
     BOND: DO iub = 1, nub
        CALL get_molecule ( molecule, first_atom=first_atom )
        index_a = ub_list ( iub ) % a + first_atom - 1
        index_b = ub_list ( iub ) % c + first_atom - 1
        rij = particle_set ( index_a ) % r - particle_set ( index_b ) % r
        CALL force_bonds(rij,ub_list ( iub ) % ub_kind % r0, &
                         ub_list ( iub ) % ub_kind % k, &
                         e0,fscalar)
        energy = energy + e0
     END DO BOND
    END DO
  END DO MOL

END SUBROUTINE potential_ub

!!*****
!******************************************************************************
!!****** fist_force_numer/potential_bend [1.0] *
!!
!!   NAME
!!     potential_bend
!!
!!   FUNCTION
!!     Computes analytically the potential from the bends
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE potential_bend ( molecule_set, molecule_kind_set, local_molecules, &
                           particle_set, energy )


    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: energy

    INTEGER                                  :: first_atom, i, ibend, ikind, &
                                                imol, index_a, index_b, &
                                                index_c, nbends, nkind, &
                                                nmol_per_kind
    REAL(KIND=dp)                            :: d12, d32, dist, e0, fscalar, &
                                                id12, id32, theta
    REAL(KIND=dp), DIMENSION(3)              :: b12, b32, g1, g2, g3
    TYPE(bend_type), POINTER                 :: bend_list( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule

!------------------------------------------------------------------------------

  nkind = SIZE ( molecule_kind_set )
  energy = 0.0_dp
  e0 = 0.0_dp
  fscalar = 0.0_dp
  MOL: DO ikind = 1, nkind
     nmol_per_kind = local_molecules % n_el ( ikind )
     DO imol = 1, nmol_per_kind
       i = local_molecules % list ( ikind ) % array ( imol )
       molecule => molecule_set ( i )
       molecule_kind => molecule % molecule_kind
       CALL get_molecule_kind ( molecule_kind, nbend = nbends, &
                                bend_list = bend_list )
     BEND: DO ibend = 1, nbends
        CALL get_molecule ( molecule, first_atom=first_atom )
        index_a = bend_list ( ibend ) % a + first_atom - 1
        index_b = bend_list ( ibend ) % b + first_atom - 1
        index_c = bend_list ( ibend ) % c + first_atom - 1
        b12 = particle_set ( index_a ) % r - particle_set ( index_b ) % r
        b32 = particle_set ( index_c ) % r - particle_set ( index_b ) % r
        d12 = SQRT(DOT_PRODUCT(b12,b12))
        id12 = 1.0_dp/d12
        d32 = SQRT(DOT_PRODUCT(b32,b32))
        id32 = 1.0_dp/d32
        dist = DOT_PRODUCT(b12,b32)
        theta = ACOS(dist*id12*id32)
        CALL force_bends(b12,b32,d12,d32,id12,id32,dist,theta, &
                         bend_list ( ibend ) % bend_kind % theta0, &
                         bend_list ( ibend ) % bend_kind % k, &
                         g1,g2,g3,e0,fscalar)
        energy = energy + e0
     END DO BEND
    END DO
  END DO MOL

END SUBROUTINE potential_bend

!!*****
!******************************************************************************
!!****** fist_force_numer/potential_torsion [1.0] *
!!
!!   NAME
!!     potential_torsion
!!
!!   FUNCTION
!!     Computes analytically the potential from the torsions
!!
!!   AUTHOR
!!     DG
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE potential_torsion ( molecule_set, molecule_kind_set, local_molecules, &
                               particle_set,energy )


    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: energy

    INTEGER :: first_atom, i, ikind, imol, imul, index_a, index_b, index_c, &
      index_d, itorsion, nkind, nmol_per_kind, ntorsions
    REAL(KIND=dp)                            :: dist1, dist2, e0, fscalar, &
                                                is32, iw1, iw2, s32, w1, w2
    REAL(KIND=dp), DIMENSION(3)              :: gt1, gt2, gt3, gt4, t12, t32, &
                                                t43, u1, u2
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule
    TYPE(torsion_type), POINTER              :: torsion_list( : )

!------------------------------------------------------------------------------

  nkind = SIZE ( molecule_kind_set )
  energy = 0.0_dp
  e0 = 0.0_dp
  fscalar = 0.0_dp
  MOL: DO ikind = 1, nkind
     nmol_per_kind = local_molecules % n_el ( ikind )
     DO imol = 1, nmol_per_kind
       i = local_molecules % list ( ikind ) % array ( imol )
       molecule => molecule_set ( i )
       molecule_kind => molecule % molecule_kind
       CALL get_molecule_kind ( molecule_kind, ntorsion = ntorsions, &
                                torsion_list = torsion_list )
     TORSION: DO itorsion = 1, ntorsions
        CALL get_molecule ( molecule, first_atom=first_atom )
        index_a = torsion_list ( itorsion ) % a + first_atom - 1
        index_b = torsion_list ( itorsion ) % b + first_atom - 1
        index_c = torsion_list ( itorsion ) % c + first_atom - 1
        index_d = torsion_list ( itorsion ) % d + first_atom - 1
        t12 = particle_set ( index_a ) % r - particle_set ( index_b ) % r
        t32 = particle_set ( index_c ) % r - particle_set ( index_b ) % r
        t43 = particle_set ( index_d ) % r - particle_set ( index_c ) % r
        s32 = SQRT(DOT_PRODUCT(t32,t32))
        is32 = 1.0_dp/s32
        dist1 = DOT_PRODUCT(t12,t32)
        dist2 = DOT_PRODUCT(t43,t32)
        u1 = t12 - dist1*t32*is32**2
        u2 = t43 - dist2*t32*is32**2
        w1 = SQRT(DOT_PRODUCT(u1,u1))
        w2 = SQRT(DOT_PRODUCT(u2,u2))
        iw1 = 1.0_dp/w1
        iw2 = 1.0_dp/w2
       DO imul = 1, torsion_list(itorsion) % torsion_kind % nmul
        CALL force_torsions(is32,dist1,dist2,u1,u2,iw1,iw2, &
                            torsion_list(itorsion)%torsion_kind%a, &
                            torsion_list(itorsion)%torsion_kind%k(imul), &
                            torsion_list(itorsion)%torsion_kind%cosphi0(imul), &
                            torsion_list(itorsion)%torsion_kind%m(imul), &
                            torsion_list(itorsion)%torsion_kind%name, &
                            gt1,gt2,gt3,gt4,e0,fscalar)
        energy = energy + e0
       END DO
     END DO TORSION
   END DO
  END DO MOL

END SUBROUTINE potential_torsion
!!*****
!******************************************************************************
!!****** fist_force_numer/potential_imptors [1.0] *
!!
!!   NAME
!!     potential_imptors
!!
!!   FUNCTION
!!     Computes analytically the potential from the torsions
!!
!!   AUTHOR
!!     DG
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE potential_imptors ( molecule_set, molecule_kind_set, local_molecules, &
                               particle_set,energy )


    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(particle_type), POINTER             :: particle_set( : )
    REAL(KIND=dp), INTENT(INOUT)             :: energy

    INTEGER :: first_atom, i, ikind, imol, index_a, index_b, index_c, &
      index_d, itorsion, nimpr, nkind, nmol_per_kind
    REAL(KIND=dp)                            :: dist1, dist2, e0, fscalar, &
                                                is32, ism, isn, s32, sm, sn
    REAL(KIND=dp), DIMENSION(3)              :: gt1, gt2, gt3, gt4, t12, t32, &
                                                t34, tm, tn
    TYPE(impr_type), POINTER                 :: impr_list( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule

!------------------------------------------------------------------------------

  nkind = SIZE ( molecule_kind_set )
  e0 = 0.0_dp
  fscalar = 0.0_dp
  energy = 0.0_dp
  MOL: DO ikind = 1, nkind
     nmol_per_kind = local_molecules % n_el ( ikind )
     DO imol = 1, nmol_per_kind
       i = local_molecules % list ( ikind ) % array ( imol )
       molecule => molecule_set ( i )
       molecule_kind => molecule % molecule_kind
       CALL get_molecule_kind ( molecule_kind, nimpr = nimpr, &
                                impr_list = impr_list )
     TORSION: DO itorsion = 1, nimpr
        CALL get_molecule ( molecule, first_atom=first_atom )
        index_a = impr_list ( itorsion ) % a + first_atom - 1
        index_b = impr_list ( itorsion ) % b + first_atom - 1
        index_c = impr_list ( itorsion ) % c + first_atom - 1
        index_d = impr_list ( itorsion ) % d + first_atom - 1
        t12 = particle_set ( index_a ) % r - particle_set ( index_b ) % r
        t32 = particle_set ( index_c ) % r - particle_set ( index_b ) % r
        t34 = particle_set ( index_c ) % r - particle_set ( index_d ) % r
! t12 x t32
        tm(1) = t12(2)*t32(3)-t32(2)*t12(3)
        tm(2) = -t12(1)*t32(3)+t32(1)*t12(3)
        tm(3) = t12(1)*t32(2)-t32(1)*t12(2)
! t32 x t34
        tn(1) = t32(2)*t34(3)-t34(2)*t32(3)
        tn(2) = -t32(1)*t34(3)+t34(1)*t32(3)
        tn(3) = t32(1)*t34(2)-t34(1)*t32(2)
        sm = SQRT ( DOT_PRODUCT(tm,tm) )
        ism = 1.0_dp/sm
        sn = SQRT ( DOT_PRODUCT(tn,tn) )
        isn = 1.0_dp/sn
        s32 = SQRT ( DOT_PRODUCT(t32,t32) )
        is32 = 1.0_dp/s32
        dist1 = DOT_PRODUCT(t12,t32)
        dist2 = DOT_PRODUCT(t34,t32)
        CALL force_imp_torsions(s32, is32,ism,isn,dist1,dist2,tm,tn, t12, t32, &
                            impr_list ( itorsion ) % impr_kind % k, &
                            impr_list ( itorsion ) % impr_kind % phi0, &
                            gt1,gt2,gt3,gt4,e0,fscalar)
        energy = energy + e0
     END DO TORSION
   END DO
  END DO MOL

END SUBROUTINE potential_imptors

!!*****
!******************************************************************************
!!****** fist_force_numer/potential_nonbond [1.0] *
!!
!!   NAME
!!     getv_ptens
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

SUBROUTINE potential_nonbond ( fist_nonbond_env, particle_set, atomic_kind_set, &
                               ewald_env, cell, vnb )


    TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), INTENT(INOUT)             :: vnb

    INTEGER                                  :: ab, atom_a, atom_b, ikind, &
                                                ilist, inode, jkind, nkinds, &
                                                nlist, nnode
    REAL(KIND=dp)                            :: energy, qi, qj, rab2, rab2_max
    REAL(KIND=dp), DIMENSION(3)              :: dra, drb, rab
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: r_last_update
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: nonbonded
    TYPE(neighbor_list_set_type), POINTER    :: neighbor_list_set
    TYPE(neighbor_list_type), POINTER        :: neighbor_list
    TYPE(neighbor_node_type), POINTER        :: neighbor_node
    TYPE(pair_potential_pp_type), POINTER    :: potparm

!------------------------------------------------------------------------------

  CALL fist_nonbond_env_get ( fist_nonbond_env, nonbonded = nonbonded, &
                              r_last_update = r_last_update,            &
                              potparm = potparm )
!
! initializing the potential energy
!
  vnb = 0.0_dp

  nkinds = SIZE ( potparm%pot, 1 )

!
! starting the force loop
!
    DO ikind=1,nkinds
      DO jkind=1,nkinds

        ab = ikind + nkinds*(jkind - 1)

        neighbor_list_set => nonbonded(ab)%neighbor_list_set

        IF (.NOT.ASSOCIATED(neighbor_list_set)) CYCLE
! obtain charges
        atomic_kind => atomic_kind_set ( ikind )
        CALL get_atomic_kind ( atomic_kind, qeff = qi )
        atomic_kind => atomic_kind_set ( jkind )
        CALL get_atomic_kind ( atomic_kind, qeff = qj )

        rab2_max = potparm %pot( ikind, jkind )%pot % rcutsq 

        CALL get_neighbor_list_set(neighbor_list_set=neighbor_list_set,&
                                   nlist=nlist)

        neighbor_list => first_list(neighbor_list_set)

        DO ilist=1,nlist

          CALL get_neighbor_list(neighbor_list=neighbor_list,&
                                 atom=atom_a,&
                                 nnode=nnode)

          dra(:) = particle_set(atom_a)%r(:) - r_last_update(:,atom_a)

          neighbor_node => first_node(neighbor_list)

          DO inode=1,nnode

            CALL get_neighbor_node(neighbor_node=neighbor_node,&
                                   neighbor=atom_b,&
                                   r=rab)

            drb(:) = particle_set(atom_b)%r(:) - r_last_update(:,atom_b)
            rab(:) = rab(:) - dra(:) + drb(:)
            rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)

            IF (rab2 <= rab2_max) THEN
              IF (qi==0.AND.qj==0) THEN
                CALL potential_f ( rab2, potparm, qi, qj, &
                                   ikind, jkind, energy )
              ELSE
                CALL potential_f ( rab2, potparm, qi, qj, &
                          ikind, jkind, energy, ewald_env )
              ENDIF
!
! summing up the potential energy
!
                 vnb = vnb + energy
            END IF

            neighbor_node => next(neighbor_node)

          END DO

          neighbor_list => next(neighbor_list)

        END DO

      END DO

    END DO

END SUBROUTINE potential_nonbond

! *****************************************************************************
  SUBROUTINE force_nonbond_numer (atomic_kind_set,particle_set,&
                                  local_particles, ewald_env, cell, &
                                  potparm, delta, f_numer,    &
                                  energy_numer, globenv)

!   Purpose: Build all the required neighbor lists for FIST.

!   History: - Creation (19.11.2002,MK)

!   ***************************************************************************

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(cell_type), POINTER                 :: cell
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    REAL(KIND=dp), INTENT(IN)                :: delta
    REAL(KIND=dp), INTENT(INOUT)             :: f_numer( :, : ), energy_numer
    TYPE(global_environment_type), POINTER   :: globenv

    CHARACTER(LEN=*), PARAMETER :: routineN = 'force_nonbond_numer'

    INTEGER :: atom_a, handle, iatom, iatom_local, ikind, istat, maxatom, &
      maxatom_local, natom_a, natom_local_a, nkind
    REAL(KIND=dp)                            :: qeff
    REAL(KIND=dp), DIMENSION(3)              :: r_pbc
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(local_atoms_type), ALLOCATABLE, &
      DIMENSION(:)                           :: atom

!   ---------------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             maxatom=maxatom)
! Zero 
    energy_numer = 0.0_dp
    f_numer = 0.0_dp

!   *** Allocate work storage ***

    nkind = SIZE(atomic_kind_set)

    ALLOCATE (atom(nkind),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "atom",nkind*int_size)

!   *** Initialize the local data structures ***

    maxatom_local = 0

    DO ikind=1,nkind

      atomic_kind => atomic_kind_set(ikind)

      NULLIFY (atom(ikind)%list)
      NULLIFY (atom(ikind)%list_local_a_index)
      NULLIFY (atom(ikind)%list_local_b_index)
      NULLIFY (atom(ikind)%r_pbc)
      NULLIFY (atom(ikind)%s_pbc)
      atom ( ikind ) % qeff = 0.0_dp

      CALL get_atomic_kind(atomic_kind=atomic_kind,&
                           atom_list=atom(ikind)%list, &
                           qeff = qeff)

      natom_a = SIZE(atom(ikind)%list)
      natom_local_a = local_particles%n_el(ikind)

      maxatom_local = MAX(maxatom_local,natom_local_a)

        IF (natom_local_a > 0) THEN

          ALLOCATE (atom(ikind)%list_local_a_index(natom_local_a),STAT=istat)
          IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                           "atom(ikind)%list_local_a_index",&
                                           natom_local_a*int_size)

!         *** Build index vector for mapping ***

          DO iatom_local=1,natom_local_a
            atom_a = local_particles%list(ikind)%array(iatom_local)
            DO iatom=1,natom_a
              IF (atom_a == atom(ikind)%list(iatom)) THEN
                atom(ikind)%list_local_a_index(iatom_local) = iatom
                EXIT
              END IF
            END DO
          END DO

        END IF

        ALLOCATE (atom(ikind)%r_pbc(3,natom_a),STAT=istat)
        IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                         "atom(ikind)%r_pbc",3*natom_a*dp_size)

        ALLOCATE (atom(ikind)%s_pbc(3,natom_a),STAT=istat)
        IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                         "atom(ikind)%s_pbc",3*natom_a*dp_size)

!       *** Calculate PBC coordinates ***

        DO iatom=1,natom_a
          atom_a = atom(ikind)%list(iatom)
          r_pbc(:) = pbc(particle_set(atom_a)%r(:),cell )
          atom(ikind)%r_pbc(:,iatom) = r_pbc(:)
          atom(ikind)%s_pbc(:,iatom) = real_to_scaled(r_pbc(:),cell )
        END DO

        atom ( ikind ) % qeff = qeff

    END DO

!   *** evaluate the forces ***
    CALL force_nonbond_evaluate ( particle_set, atom, ewald_env, cell, &
                                  potparm, energy_numer, f_numer, delta,   &
                                  maxatom_local )

    f_numer = f_numer/(2.0_dp*delta)
!   *** Release work storage ***

    DO ikind=1,nkind
      NULLIFY (atom(ikind)%list)
      IF (ASSOCIATED(atom(ikind)%list_local_a_index)) THEN
        DEALLOCATE (atom(ikind)%list_local_a_index,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                         "atom(ikind)%list_local_a_index")
      END IF
      IF (ASSOCIATED(atom(ikind)%r_pbc)) THEN
        DEALLOCATE (atom(ikind)%r_pbc,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                         "atom(ikind)%r_pbc")
      END IF
      IF (ASSOCIATED(atom(ikind)%s_pbc)) THEN
        DEALLOCATE (atom(ikind)%s_pbc,STAT=istat)
        IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                         "atom(ikind)%s_pbc")
      END IF
    END DO

    DEALLOCATE (atom,STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                     "atom")

    CALL timestop(0.0_dp,handle)


  END SUBROUTINE force_nonbond_numer

! *****************************************************************************

  SUBROUTINE force_nonbond_evaluate ( particle_set, atom, ewald_env, cell, &
                                      potparm, energy_numer, f_numer,    &
                                      numerical_shift, maxatom_local )

!   Purpose: Build simple pair neighbor lists.

!   History: - Creation (20.03.2002,MK)

!   ***************************************************************************

    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(local_atoms_type), DIMENSION(:), &
      INTENT(IN)                             :: atom
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(cell_type), POINTER                 :: cell
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    REAL(KIND=dp), INTENT(OUT)               :: energy_numer
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(OUT)                            :: f_numer
    REAL(KIND=dp), INTENT(IN)                :: numerical_shift
    INTEGER, INTENT(IN)                      :: maxatom_local

    CHARACTER(LEN=*), PARAMETER :: routineN = 'force_nonbond_evaluate'

    INTEGER :: atom_a, atom_b, handle, i, iatom, iatom_local, icell, ikind, &
      jatom_local, jcell, jkind, kcell, natom_local_a, natom_local_b, nkind
    INTEGER, DIMENSION(3)                    :: cell_b, ncell, periodic
    INTEGER, DIMENSION(:), POINTER           :: exclusion_list
    LOGICAL                                  :: exclude_flag, include_ab
    REAL(KIND=dp)                            :: qi, qj, rab2, rab2_max
    REAL(KIND=dp), DIMENSION(3)              :: r, rab, rb, sab, sab_max, sb, &
                                                sb_max, sb_min, sb_pbc

!   ---------------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)
   
    f_numer=0.0_dp
    energy_numer=0.0_dp

!   *** Allocate and initialize the new neighbor list structure ***

    nkind = SIZE(atom)

!   *** Load informations about the simulation cell ***

    CALL get_cell(cell=cell,&
                  periodic=periodic)

!   *** Loop over all atomic kind pairs ***

    DO ikind=1,nkind

      IF (.NOT.ASSOCIATED(atom(ikind)%list_local_a_index)) CYCLE

      natom_local_a = SIZE(atom(ikind)%list_local_a_index)
      qi = atom(ikind)%qeff

      DO jkind=1,nkind

        natom_local_b = SIZE(atom(jkind)%list)
        qj = atom(jkind)%qeff

!       *** Calculate the square of the maximum interaction distance ***


        rab2_max = potparm %pot( ikind, jkind )%pot % rcutsq 
        r(:) = SQRT ( rab2_max )
        sab_max(1) = r(1)/plane_distance(1,0,0,cell)
        sab_max(2) = r(2)/plane_distance(0,1,0,cell)
        sab_max(3) = r(3)/plane_distance(0,0,1,cell)

        ncell(:) = (INT(sab_max(:)) + 1)*periodic(:)

!       *** No subcells will be considered in the debug for safety check****

!         *** Case 1: No subcell grid is needed ***

          DO iatom_local=1,natom_local_a
            iatom = atom(ikind)%list_local_a_index(iatom_local)
            atom_a = atom(ikind)%list(iatom)
          END DO

          DO jatom_local=1,natom_local_b

            atom_b = atom(jkind)%list(jatom_local)
            sb_pbc(:) = atom(jkind)%s_pbc(:,jatom_local)

            loop1_kcell: DO kcell=-ncell(3),ncell(3)

              sb(3) = sb_pbc(3) + REAL(kcell,dp)
              sb_min(3) = sb(3) - sab_max(3)
              IF (sb_min(3) >= 0.5_dp) EXIT loop1_kcell
              sb_max(3) = sb(3) + sab_max(3)
              IF (sb_max(3) < -0.5_dp) CYCLE loop1_kcell
              cell_b(3) = kcell

              loop1_jcell: DO jcell=-ncell(2),ncell(2)

                sb(2) = sb_pbc(2) + REAL(jcell,dp)
                sb_min(2) = sb(2) - sab_max(2)
                IF (sb_min(2) >= 0.5_dp) EXIT loop1_jcell
                sb_max(2) = sb(2) + sab_max(2)
                IF (sb_max(2) < -0.5_dp) CYCLE loop1_jcell
                cell_b(2) = jcell

                loop1_icell: DO icell=-ncell(1),ncell(1)

                  sb(1) = sb_pbc(1) + REAL(icell,dp)
                  sb_min(1) = sb(1) - sab_max(1)
                  IF (sb_min(1) >= 0.5_dp) EXIT loop1_icell
                  sb_max(1) = sb(1) + sab_max(1)
                  IF (sb_max(1) < -0.5_dp) CYCLE loop1_icell
                  cell_b(1) = icell

                  rb(:) = scaled_to_real(sb(:),cell)

                  DO iatom_local=1,natom_local_a
                    iatom = atom(ikind)%list_local_a_index(iatom_local)
                    atom_a = atom(ikind)%list(iatom)
                    IF (atom_a > atom_b) THEN
                      include_ab = (MODULO(atom_a + atom_b,2) /= 0)
                    ELSE
                      include_ab = (MODULO(atom_a + atom_b,2) == 0)
                    END IF
                    IF (include_ab) THEN
                      rab(:) = rb(:) - atom(ikind)%r_pbc(:,iatom)
                      rab2 = rab(1)*rab(1) + rab(2)*rab(2) + rab(3)*rab(3)
                      IF (rab2 < rab2_max ) THEN
                        IF (rab2 > 1.0E-6_dp) THEN
                          sab(:) = real_to_scaled(rab(:),cell)
                          exclude_flag = .FALSE.
                          exclusion_list =>  &
                          particle_set ( atom_a ) % lpi % exclusion ( jkind ) % list_vdw
                          IF ( ASSOCIATED ( exclusion_list ) ) THEN
                            DO iatom=1,SIZE(exclusion_list)
                              IF (exclusion_list(iatom) == 0) EXIT
                              IF (exclusion_list(iatom) == atom_b) &
                              exclude_flag = .TRUE. 
                            END DO
                          END IF
                          IF (.NOT.exclude_flag ) THEN
                            SELECT CASE ( TRIM ( potparm %pot( ikind, jkind )%pot % type ) )
                            CASE ( "EAM" )
                              CALL stop_program ( "force_nonbond_numer", &
                                "DEBUG results for EAM are in regtest" )
                            CASE DEFAULT
                              CALL finite_difference_evaluate ( ewald_env,     &
                                                              potparm, rab,    &
                                                              ikind, jkind, &
                                                              atom_a, atom_b, &
                                                              qi, qj, &
                                                              numerical_shift, &
                                                              f_numer,     &
                                                              energy_numer )
                            END SELECT
                          END IF
                        END IF
                      END IF
                    END IF
                  END DO

                END DO loop1_icell
              END DO loop1_jcell
            END DO loop1_kcell

          END DO

      END DO
    END DO

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE force_nonbond_evaluate
!******************************************************************************

SUBROUTINE finite_difference_evaluate ( ewald_env, potparm, rab,              &
                                        ikind, jkind, atom_a, atom_b, qi, qj, &
                                        numerical_shift, f_numer, energy_numer )


    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(pair_potential_pp_type), POINTER    :: potparm
    REAL(KIND=dp), DIMENSION(3), INTENT(IN)  :: rab
    INTEGER, INTENT(IN)                      :: ikind, jkind, atom_a, atom_b
    REAL(KIND=dp), INTENT(IN)                :: qi, qj, numerical_shift
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: f_numer
    REAL(KIND=dp), INTENT(INOUT)             :: energy_numer

    INTEGER                                  :: id
    REAL(KIND=dp)                            :: energy, energy_minus, &
                                                energy_plus, rab2, &
                                                rij_minus_delta_sq, &
                                                rij_plus_delta_sq
    REAL(KIND=dp), DIMENSION(3)              :: delta

  delta = 0.0_dp
  rab2 = DOT_PRODUCT(rab,rab)
  DIM_LOOP: DO id = 1, 3
     delta(id) = numerical_shift
     rij_minus_delta_sq = DOT_PRODUCT(rab-delta,rab-delta)
     rij_plus_delta_sq = DOT_PRODUCT(rab+delta,rab+delta)
     IF (qi==0.AND.qj==0) THEN
        CALL potential_f(rij_minus_delta_sq,potparm,qi,qj, &
             ikind, jkind, energy_minus )
        CALL potential_f(rij_plus_delta_sq,potparm,qi,qj, &
             ikind, jkind, energy_plus )
     ELSE
        CALL potential_f(rij_minus_delta_sq,potparm,qi,qj, &
             ikind, jkind, energy_minus, ewald_env )
        CALL potential_f(rij_plus_delta_sq,potparm,qi,qj, &
             ikind, jkind, energy_plus, ewald_env )
     ENDIF
     f_numer(id,atom_a) = f_numer(id,atom_a) + (energy_plus-energy_minus)
     f_numer(id,atom_b) = f_numer(id,atom_b) - (energy_plus-energy_minus)
     delta = 0.0_dp
  END DO DIM_LOOP
  ! compute energy
  IF (qi==0.AND.qj==0) THEN
     CALL potential_f ( rab2, potparm, qi, qj, &
          ikind, jkind, energy )
  ELSE
     CALL potential_f ( rab2, potparm, qi, qj, &
          ikind, jkind, energy, ewald_env )
  ENDIF
  energy_numer = energy_numer + energy
  
END SUBROUTINE finite_difference_evaluate

!!*****
!******************************************************************************

END MODULE fist_force_numer

!******************************************************************************
