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

!  USE kinds,                           ONLY: dbl
!  USE mol_force,                       ONLY: force_bends,&
!                                             force_bonds,&
!                                             force_torsions,&
!                                             get_pv_bend,&
!                                             get_pv_bond,&
!                                             get_pv_torsion
!  USE timings,                         ONLY: timeset,&
!                                             timestop
!
!  IMPLICIT NONE
!
!  PRIVATE
!  PUBLIC :: force_intra_control
!
!!!*****
!!******************************************************************************
!
!CONTAINS
!
!!******************************************************************************
!!!****** fist_intra_force/force_intra_control [1.0] *
!!!
!!!   NAME
!!!     force_intra_control
!!!
!!!   FUNCTION
!!!     Computes the the intramolecular energies, forces, and pressure tensors
!!!
!!!   AUTHOR
!!!     CJM
!!!
!!!   MODIFICATION HISTORY
!!!     none
!!!
!!!   SOURCE
!!******************************************************************************
!
!SUBROUINTE force_intra_control ( molecule_set, molecule_kind_set, local_molecules, &
!                               particle_set, pot_bond, pot_bend, pot_torsion, &
!                               pv_bond, pv_bend, pv_torsion, f_bond, f_bend, f_torsion )
!
!! parses up intramolecular properties
!
!
!    REAL(dbl), INTENT(INOUT)                 :: pot_bond, pot_bend, &
!                                                pot_torsion
!    REAL(dbl), DIMENSION(:, :), &
!      INTENT(INOUT)                          :: pv_bond, pv_bend, pv_torsion
!    REAL(dbl), DIMENSION(:, :), &
!      INTENT(OUT), OPTIONAL                  :: f_bond, f_bend, f_torsion
!
!    INTEGER                                  :: handle, ibend, ibond, imol, &
!                                                itorsion
!    REAL(dbl)                                :: d12, d32, dist, dist1, dist2, &
!                                                energy, fscalar, id12, id32, &
!                                                is32, iw1, iw2, s32, theta, &
!                                                w1, w2
!    REAL(dbl), DIMENSION(3)                  :: b12, b32, g1, g2, g3, gt1, &
!                                                gt2, gt3, gt4, rij, t12, t32, &
!                                                t43, u1, u2
!
!!------------------------------------------------------------------------------
!
!  CALL timeset ( 'FORCE_INTRA_CONTROL','I',' ',handle )
!
!  IF ( PRESENT ( f_bond)) f_bond = 0._dbl
!  IF ( PRESENT ( f_bend)) f_bend = 0._dbl
!  IF ( PRESENT ( f_torsion)) f_torsion = 0._dbl
!  pot_bond =  0._dbl
!  pot_bend =  0._dbl
!  pot_torsion =  0._dbl
!
!  nkind = SIZE ( molecule_kind_set )
!  MOL: DO ikind = 1, nkind
!     nmol_per_kind = local_molecules % n_el ( ikind )
!     DO imol = 1, nmol_per_kind
!       i = local_molecules % list ( ikind ) % array ( imol )
!       molecule => molecule_set ( i )
!       molecule_kind => molecule % molecule_kind
!       CALL get_molecule_kind ( molecule_kind, nbends = nbends,  &
!                                nbonds = nbonds, ntorsions = ntorsions, &
!                                bond_list = bond_list, bend_list = bend_list, &
!                                torsion_list = torsion_list )
!
!     BOND: DO ibond = 1, nbonds
!        CALL get_molecule ( molecule, first_atom=first_atom, lbonds=lbonds )
!        index_a = lbonds ( ibond ) % a + first_atom - 1
!        index_b = lbonds ( ibond ) % b + first_atom - 1
!        rij = particle_set ( index_a ) % r - particle_set ( index_b ) % r
!        CALL force_bonds(rij,llbond%bond_param%r0,llbond%bond_param%k, &
!             energy,fscalar)
!        pot_bond = pot_bond + energy
!        llbond%p1%f = llbond%p1%f - rij*fscalar
!        llbond%p2%f = llbond%p2%f + rij*fscalar
!
!! computing the pressure tensor
!        CALL get_pv_bond ( -rij * fscalar, rij * fscalar, &
!             llbond%p1%r, llbond%p2%r, pv_bond )
!
!! the contribution from the bonds. ONLY FOR DEBUG
!        IF ( PRESENT ( f_bond)) THEN
!           f_bond(:,llbond%index(1)) = f_bond(:,llbond%index(1)) - &
!                rij*fscalar
!           f_bond(:,llbond%index(2)) = f_bond(:,llbond%index(2)) + &
!                rij*fscalar
!        END IF
!
!        llbond => llbond%next
!     END DO BOND
!
!     BEND: DO ibend = 1, molecule(imol) %nbends_mol
!        b12 = llbend%p1%r - llbend%p2%r
!        b32 = llbend%p3%r - llbend%p2%r
!        d12 = sqrt(dot_product(b12,b12))
!        id12 = 1._dbl/d12
!        d32 = sqrt(dot_product(b32,b32))
!        id32 = 1._dbl/d32
!        dist = dot_product(b12,b32)
!        theta = acos(dist*id12*id32)
!        CALL force_bends(b12,b32,d12,d32,id12,id32,dist,theta, &
!             llbend%bend_param%theta0,llbend%bend_param%k,g1,g2,g3,energy, &
!             fscalar)
!        pot_bend = pot_bend + energy
!        llbend%p1%f = llbend%p1%f + fscalar*g1
!        llbend%p2%f = llbend%p2%f + fscalar*g2
!        llbend%p3%f = llbend%p3%f + fscalar*g3
!
!! computing the pressure tensor
!        CALL get_pv_bend(fscalar*g1,fscalar*g2,fscalar*g3,llbend%p1%r, &
!             llbend%p2%r,llbend%p3%r,pv_bend)
!
!! the contribution from the bends. ONLY FOR DEBUG
!        IF ( PRESENT ( f_bend)) THEN
!           f_bend(:,llbend%index(1)) = f_bend(:,llbend%index(1)) + &
!                fscalar*g1
!           f_bend(:,llbend%index(2)) = f_bend(:,llbend%index(2)) + &
!                fscalar*g2
!           f_bend(:,llbend%index(3)) = f_bend(:,llbend%index(3)) + &
!                fscalar*g3
!        END IF
!
!        llbend => llbend%next
!     END DO BEND
!
!     TORSION: DO itorsion = 1, molecule(imol) %ntorsions_mol
!        t12 = lltorsion%p1%r - lltorsion%p2%r
!        t32 = lltorsion%p3%r - lltorsion%p2%r
!        t43 = lltorsion%p4%r - lltorsion%p3%r
!        s32 = sqrt(dot_product(t32,t32))
!        is32 = 1._dbl/s32
!        dist1 = dot_product(t12,t32)
!        dist2 = dot_product(t43,t32)
!        u1 = t12 - dist1*t32*is32**2
!        u2 = t43 - dist2*t32*is32**2
!        w1 = sqrt(dot_product(u1,u1))
!        w2 = sqrt(dot_product(u2,u2))
!        iw1 = 1._dbl/w1
!        iw2 = 1._dbl/w2
!        CALL force_torsions(is32,dist1,dist2,u1,u2,iw1,iw2, &
!             lltorsion%torsion_param%a,gt1,gt2,gt3,gt4,energy,fscalar)
!        pot_torsion = pot_torsion + energy
!        lltorsion%p1%f = lltorsion%p1%f + fscalar*gt1
!        lltorsion%p2%f = lltorsion%p2%f + fscalar*gt2
!        lltorsion%p3%f = lltorsion%p3%f + fscalar*gt3
!        lltorsion%p4%f = lltorsion%p4%f + fscalar*gt4
!
!! computing the pressure tensor
!        CALL get_pv_torsion(fscalar*gt1,fscalar*gt2,fscalar*gt3,fscalar*gt4, &
!             lltorsion%p1%r,lltorsion%p2%r,lltorsion%p3%r,lltorsion%p4%r, &
!             pv_torsion)
!
!! the contribution from the torsions. ONLY FOR DEBUG
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
!
!        lltorsion => lltorsion%next
!     END DO TORSION
!
!  END DO MOL
!
!  CALL timestop ( 0._dbl, handle )
!
!END SUBROUTINE force_intra_control
!
!!!*****
!!******************************************************************************

END MODULE fist_intra_force

!******************************************************************************
