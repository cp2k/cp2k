!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/fist_intra_force [1.0] *
!!
!!   NAME
!!     fist_intra_force
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

MODULE fist_intra_force

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE kinds,                           ONLY: dbl
  USE mol_force,                       ONLY: force_bends,&
                                             force_bonds,&
                                             force_torsions,&
                                             force_imp_torsions,&
                                             get_pv_bend,&
                                             get_pv_bond,&
                                             get_pv_torsion
  USE molecule_kind_types,             ONLY: get_molecule_kind,&
                                             molecule_kind_type,&
                                             bond_type, &
                                             bend_type, &
                                             ub_type, &
                                             impr_type, &
                                             torsion_type
  USE molecule_types_new,              ONLY: molecule_type, & 
                                             get_molecule
  USE particle_types,                  ONLY: particle_type
  USE timings,                         ONLY: timeset,&
                                             timestop

  IMPLICIT NONE

 PRIVATE
  PUBLIC :: force_intra_control

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** fist_intra_force/force_intra_control [1.0] *
!!
!!   NAME
!!     force_intra_control
!!
!!   FUNCTION
!!     Computes the the intramolecular energies, forces, and pressure tensors
!!
!!   AUTHOR
!!     CJM
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!!******************************************************************************

SUBROUTINE force_intra_control ( molecule_set, molecule_kind_set, local_molecules, &
                               particle_set, pot_bond, pot_bend, pot_urey_bradley, &
                               pot_torsion, pot_imp_torsion, pv_bond, pv_bend, pv_urey_bradley, &
                               pv_torsion, pv_imp_torsion, f_bond, f_bend, f_torsion )

! parses up intramolecular properties


  TYPE(distribution_1d_type), POINTER      :: local_molecules
  TYPE(molecule_type), POINTER             :: molecule_set( : )
  TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
  TYPE(particle_type), POINTER             :: particle_set ( : )
  REAL(dbl), INTENT(INOUT)                 :: pot_bond, pot_bend, &
                                                pot_torsion, pot_urey_bradley, &
                                                pot_imp_torsion

  REAL(dbl), DIMENSION(:, :), &
    INTENT(INOUT)                          :: pv_bond, pv_bend, pv_torsion, &
                                              pv_urey_bradley, pv_imp_torsion
  REAL(dbl), DIMENSION(:, :), &
    INTENT(OUT), OPTIONAL                  :: f_bond, f_bend, f_torsion
 
  INTEGER                                  :: handle, ibend, ibond, imol, &
                                             itorsion, ikind, nkind, index_a,&
                                             index_b, index_c, index_d, i, &
                                             nbonds, nbends, ntorsions, first_atom, &
                                             nmol_per_kind, nub, nimptors
  REAL(dbl)                                :: d13, d12, d32, dist, dist1, dist2, &
                                             energy, fscalar, id12, id32, &
                                             is32, iw1, iw2, s32, theta, &
                                                w1, w2
  REAL(dbl), DIMENSION(3)                  :: b13, b12, b32, g1, g2, g3, gt1, &
                                                gt2, gt3, gt4, rij, t12, t32, &
                                                t43, u1, u2
  TYPE ( bond_type ), POINTER              :: bond_list ( : )
  TYPE ( bend_type ), POINTER              :: bend_list ( : )
  TYPE ( ub_type ), POINTER                :: ub_list ( : )
  TYPE ( impr_type ), POINTER              :: impr_list ( : )
  TYPE ( torsion_type ), POINTER           :: torsion_list ( : )
  TYPE ( molecule_type ), POINTER          :: molecule
  TYPE ( molecule_kind_type ), POINTER      :: molecule_kind

!------------------------------------------------------------------------------

  CALL timeset ( 'FORCE_INTRA_CONTROL','I',' ',handle )

  IF ( PRESENT ( f_bond)) f_bond = 0._dbl
  IF ( PRESENT ( f_bend)) f_bend = 0._dbl
  IF ( PRESENT ( f_torsion)) f_torsion = 0._dbl
  pot_bond =  0._dbl
  pot_bend =  0._dbl
  pot_urey_bradley =  0._dbl
  pot_torsion =  0._dbl
  pot_imp_torsion =  0._dbl

  nkind = SIZE ( molecule_kind_set )
  MOL: DO ikind = 1, nkind
     nmol_per_kind = local_molecules % n_el ( ikind )
     DO imol = 1, nmol_per_kind
       i = local_molecules % list ( ikind ) % array ( imol )
       molecule => molecule_set ( i )
       molecule_kind => molecule % molecule_kind
       CALL get_molecule_kind ( molecule_kind, nbend = nbends, nbond = nbonds,  &
                                nimpr = nimptors, nub = nub, ntorsion = ntorsions, &
                                bond_list = bond_list, bend_list = bend_list, &
                                torsion_list = torsion_list )

     BOND: DO ibond = 1, nbonds
        CALL get_molecule ( molecule, first_atom=first_atom )
        index_a = bond_list ( ibond ) % a + first_atom - 1
        index_b = bond_list ( ibond ) % b + first_atom - 1
        rij = particle_set ( index_a ) % r - particle_set ( index_b ) % r
        CALL force_bonds(rij,bond_list ( ibond ) % bond_kind % r0, &
                         bond_list ( ibond ) % bond_kind % k, &
                         energy,fscalar)
        pot_bond = pot_bond + energy
        particle_set ( index_a )% f = particle_set ( index_a ) % f &
                                      -rij * fscalar
        particle_set ( index_b )% f = particle_set ( index_b ) % f &
                                      +rij * fscalar

! computing the pressure tensor
        CALL get_pv_bond ( -rij * fscalar, rij * fscalar, &
                           particle_set ( index_a ) % r,  &
                           particle_set ( index_b ) % r,  &
                           pv_bond )

! the contribution from the bonds. ONLY FOR DEBUG
!        IF ( PRESENT ( f_bond ) ) THEN
!           f_bond(:,llbond%index(1)) = f_bond(:,llbond%index(1)) - &
!                rij*fscalar
!           f_bond(:,llbond%index(2)) = f_bond(:,llbond%index(2)) + &
!                rij*fscalar
!        END IF

     END DO BOND

     UREY_BRADLEY: DO ibend = 1, nub
        CALL get_molecule ( molecule, first_atom=first_atom )
        index_a = ub_list ( ibend ) % a + first_atom - 1
        index_b = ub_list ( ibend ) % b + first_atom - 1
        index_c = ub_list ( ibend ) % c + first_atom - 1
        b13 = particle_set ( index_a ) % r - particle_set ( index_c ) % r
        d13 = sqrt(dot_product(b13,b13))
        CALL force_bonds (b13,  ub_list ( ibend ) % ub_kind % r0, &
                          ub_list ( ibend ) % ub_kind % k, energy, fscalar)
        pot_bend = pot_bend + energy
        particle_set ( index_a )% f = particle_set ( index_a ) % f &
                                      - b13 * fscalar
        particle_set ( index_c )% f = particle_set ( index_c ) % f &
                                      + b13 * fscalar

! computing the pressure tensor
        CALL get_pv_bond ( -b13 * fscalar, b13 * fscalar, &
                           particle_set ( index_a ) % r,  &
                           particle_set ( index_c ) % r,  &
                           pv_urey_bradley )

! the contribution from the bonds. ONLY FOR DEBUG
!        IF ( PRESENT ( f_bond ) ) THEN
!           f_bond(:,llbond%index(1)) = f_bond(:,llbond%index(1)) - &
!                rij*fscalar
!           f_bond(:,llbond%index(2)) = f_bond(:,llbond%index(2)) + &
!                rij*fscalar
!        END IF

     END DO UREY_BRADLEY


     BEND: DO ibend = 1, nbends
        CALL get_molecule ( molecule, first_atom=first_atom )
        index_a = bend_list ( ibend ) % a + first_atom - 1
        index_b = bend_list ( ibend ) % b + first_atom - 1
        index_c = bend_list ( ibend ) % c + first_atom - 1
        b12 = particle_set ( index_a ) % r - particle_set ( index_b ) % r
        b32 = particle_set ( index_c ) % r - particle_set ( index_b ) % r
        d12 = sqrt(dot_product(b12,b12))
        id12 = 1._dbl/d12
        d32 = sqrt(dot_product(b32,b32))
        id32 = 1._dbl/d32
        dist = dot_product(b12,b32)
        theta = acos(dist*id12*id32)
        CALL force_bends(b12,b32,d12,d32,id12,id32,dist,theta, &
                         bend_list ( ibend ) % bend_kind % theta0, &
                         bend_list ( ibend ) % bend_kind % k, &
                         g1,g2,g3,energy,fscalar)
        pot_bend = pot_bend + energy
        particle_set ( index_a )% f = particle_set ( index_a ) % f &
                                      + g1 * fscalar
        particle_set ( index_b )% f = particle_set ( index_b ) % f &
                                      + g2 * fscalar
        particle_set ( index_c )% f = particle_set ( index_c ) % f &
                                      + g3 * fscalar

! computing the pressure tensor
        CALL get_pv_bend(fscalar*g1,fscalar*g2,fscalar*g3, &
                         particle_set ( index_a ) % r,    & 
                         particle_set ( index_b ) % r,    & 
                         particle_set ( index_c ) % r,    & 
                         pv_bend)

! the contribution from the bends. ONLY FOR DEBUG
!        IF ( PRESENT ( f_bend)) THEN
!           f_bend(:,llbend%index(1)) = f_bend(:,llbend%index(1)) + &
!                fscalar*g1
!           f_bend(:,llbend%index(2)) = f_bend(:,llbend%index(2)) + &
!                fscalar*g2
!           f_bend(:,llbend%index(3)) = f_bend(:,llbend%index(3)) + &
!                fscalar*g3
!        END IF

     END DO BEND

     TORSION: DO itorsion = 1, ntorsions
        CALL get_molecule ( molecule, first_atom=first_atom )
        index_a = torsion_list ( itorsion ) % a + first_atom - 1
        index_b = torsion_list ( itorsion ) % b + first_atom - 1
        index_c = torsion_list ( itorsion ) % c + first_atom - 1
        index_d = torsion_list ( itorsion ) % d + first_atom - 1
        t12 = particle_set ( index_a ) % r - particle_set ( index_b ) % r
        t32 = particle_set ( index_c ) % r - particle_set ( index_b ) % r
        t43 = particle_set ( index_d ) % r - particle_set ( index_c ) % r
        s32 = sqrt(dot_product(t32,t32))
        is32 = 1._dbl/s32
        dist1 = dot_product(t12,t32)
        dist2 = dot_product(t43,t32)
        u1 = t12 - dist1*t32*is32**2
        u2 = t43 - dist2*t32*is32**2
        w1 = sqrt(dot_product(u1,u1))
        w2 = sqrt(dot_product(u2,u2))
        iw1 = 1._dbl/w1
        iw2 = 1._dbl/w2
        CALL force_torsions(is32,dist1,dist2,u1,u2,iw1,iw2, &
                            torsion_list ( itorsion ) % torsion_kind % a, &
                            gt1,gt2,gt3,gt4,energy,fscalar)
        pot_torsion = pot_torsion + energy
        particle_set ( index_a )% f = particle_set ( index_a ) % f &
                                      + gt1 * fscalar
        particle_set ( index_b )% f = particle_set ( index_b ) % f &
                                      + gt2 * fscalar
        particle_set ( index_c )% f = particle_set ( index_c ) % f &
                                      + gt3 * fscalar
        particle_set ( index_d )% f = particle_set ( index_d ) % f &
                                      + gt3 * fscalar

! computing the pressure tensor
        CALL get_pv_torsion(fscalar*gt1,fscalar*gt2,fscalar*gt3,fscalar*gt4, &
                         particle_set ( index_a ) % r,    & 
                         particle_set ( index_b ) % r,    & 
                         particle_set ( index_c ) % r,    & 
                         particle_set ( index_d ) % r,    & 
                         pv_torsion)

! the contribution from the torsions. ONLY FOR DEBUG
!        IF ( PRESENT ( f_torsion)) THEN
!           f_torsion(:,lltorsion%index(1)) = f_torsion(:,lltorsion%index(1)) + &
!                fscalar*gt1
!           f_torsion(:,lltorsion%index(2)) = f_torsion(:,lltorsion%index(2)) + &
!                fscalar*gt2
!           f_torsion(:,lltorsion%index(3)) = f_torsion(:,lltorsion%index(3)) + &
!                fscalar*gt3
!           f_torsion(:,lltorsion%index(4)) = f_torsion(:,lltorsion%index(4)) + &
!                fscalar*gt4
!        END IF

     END DO TORSION

     IMP_TORSION: DO itorsion = 1, nimptors
        CALL get_molecule ( molecule, first_atom=first_atom )
        index_a = torsion_list ( itorsion ) % a + first_atom - 1
        index_b = torsion_list ( itorsion ) % b + first_atom - 1
        index_c = torsion_list ( itorsion ) % c + first_atom - 1
        index_d = torsion_list ( itorsion ) % d + first_atom - 1
        t12 = particle_set ( index_a ) % r - particle_set ( index_b ) % r
        t32 = particle_set ( index_c ) % r - particle_set ( index_b ) % r
        t43 = particle_set ( index_d ) % r - particle_set ( index_c ) % r
        s32 = sqrt(dot_product(t32,t32))
        is32 = 1._dbl/s32
        dist1 = dot_product(t12,t32)
        dist2 = dot_product(t43,t32)
        u1 = t12 - dist1*t32*is32**2
        u2 = t43 - dist2*t32*is32**2
        w1 = sqrt(dot_product(u1,u1))
        w2 = sqrt(dot_product(u2,u2))
        iw1 = 1._dbl/w1
        iw2 = 1._dbl/w2
        CALL force_imp_torsions(is32,dist1,dist2,u1,u2,iw1,iw2, &
                            torsion_list ( itorsion ) % torsion_kind % k, &
                            torsion_list ( itorsion ) % torsion_kind % phi0, &
                            gt1,gt2,gt3,gt4,energy,fscalar)
        pot_imp_torsion = pot_imp_torsion + energy
        particle_set ( index_a )% f = particle_set ( index_a ) % f &
                                      + gt1 * fscalar
        particle_set ( index_b )% f = particle_set ( index_b ) % f &
                                      + gt2 * fscalar
        particle_set ( index_c )% f = particle_set ( index_c ) % f &
                                      + gt3 * fscalar
        particle_set ( index_d )% f = particle_set ( index_d ) % f &
                                      + gt3 * fscalar

! computing the pressure tensor
        CALL get_pv_torsion(fscalar*gt1,fscalar*gt2,fscalar*gt3,fscalar*gt4, &
                         particle_set ( index_a ) % r,    & 
                         particle_set ( index_b ) % r,    & 
                         particle_set ( index_c ) % r,    & 
                         particle_set ( index_d ) % r,    & 
                         pv_imp_torsion)

! the contribution from the torsions. ONLY FOR DEBUG
!        IF ( PRESENT ( f_torsion)) THEN
!           f_torsion(:,lltorsion%index(1)) = f_torsion(:,lltorsion%index(1)) + &
!                fscalar*gt1
!           f_torsion(:,lltorsion%index(2)) = f_torsion(:,lltorsion%index(2)) + &
!                fscalar*gt2
!           f_torsion(:,lltorsion%index(3)) = f_torsion(:,lltorsion%index(3)) + &
!                fscalar*gt3
!           f_torsion(:,lltorsion%index(4)) = f_torsion(:,lltorsion%index(4)) + &
!                fscalar*gt4
!        END IF
     END DO IMP_TORSION
    END DO
  END DO MOL

  CALL timestop ( 0._dbl, handle )

END SUBROUTINE force_intra_control

!!*****
!******************************************************************************

END MODULE fist_intra_force

!******************************************************************************
