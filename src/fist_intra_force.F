!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      Torsions added (DG) 05-Dec-2000
!>      Variable names changed (DG) 05-Dec-2000
!> \author CJM
! *****************************************************************************
MODULE fist_intra_force

  USE atprop_types,                    ONLY: atprop_type
  USE cell_types,                      ONLY: cell_type,&
                                             pbc
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE f77_blas
  USE kinds,                           ONLY: dp
  USE mol_force,                       ONLY: force_bends,&
                                             force_bonds,&
                                             force_imp_torsions,&
                                             force_opbends,&
                                             force_torsions,&
                                             get_pv_bend,&
                                             get_pv_bond,&
                                             get_pv_torsion
  USE molecule_kind_types,             ONLY: &
       bend_type, bond_type, get_molecule_kind, impr_type, &
       molecule_kind_type, opbend_type, shell_type, torsion_type, ub_type
  USE molecule_types_new,              ONLY: get_molecule,&
                                             molecule_type
  USE particle_types,                  ONLY: particle_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'fist_intra_force'
  PUBLIC :: force_intra_control

CONTAINS

! *****************************************************************************
!> \brief Computes the the intramolecular energies, forces, and pressure tensors
!> \param molecule_set ...
!> \param molecule_kind_set ...
!> \param local_molecules ...
!> \param particle_set ...
!> \param shell_particle_set ...
!> \param core_particle_set ...
!> \param pot_bond ...
!> \param pot_bend ...
!> \param pot_urey_bradley ...
!> \param pot_torsion ...
!> \param pot_imp_torsion ...
!> \param pot_opbend ...
!> \param pot_shell ...
!> \param pv_bond ...
!> \param pv_bend ...
!> \param pv_urey_bradley ...
!> \param pv_torsion ...
!> \param pv_imp_torsion ...
!> \param pv_opbend ...
!> \param f_bond ...
!> \param f_bend ...
!> \param f_torsion ...
!> \param f_ub ...
!> \param f_imptor ...
!> \param f_opbend ...
!> \param cell ...
!> \param use_virial ...
!> \param atprop_env ...
!> \param error ...
!> \par History
!>      none
!> \author CJM
! *****************************************************************************
  SUBROUTINE force_intra_control(molecule_set, molecule_kind_set, &
       local_molecules, particle_set, shell_particle_set, core_particle_set, &
       pot_bond, pot_bend, pot_urey_bradley, pot_torsion, pot_imp_torsion, &
       pot_opbend, pot_shell, pv_bond, pv_bend, pv_urey_bradley, pv_torsion, &
       pv_imp_torsion, pv_opbend, f_bond, f_bend, f_torsion, f_ub, &
       f_imptor, f_opbend, cell, use_virial, atprop_env, error)

    TYPE(molecule_type), POINTER             :: molecule_set( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind_set( : )
    TYPE(distribution_1d_type), POINTER      :: local_molecules
    TYPE(particle_type), POINTER             :: particle_set( : ), &
                                                shell_particle_set( : ), &
                                                core_particle_set( : )
    REAL(KIND=dp), INTENT(INOUT) :: pot_bond, pot_bend, pot_urey_bradley, &
      pot_torsion, pot_imp_torsion, pot_opbend, pot_shell
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT)                          :: pv_bond, pv_bend, &
                                                pv_urey_bradley, pv_torsion, &
                                                pv_imp_torsion, pv_opbend
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT), OPTIONAL                :: f_bond, f_bend, f_torsion, &
                                                f_ub, f_imptor, f_opbend
    TYPE(cell_type), POINTER                 :: cell
    LOGICAL, INTENT(IN)                      :: use_virial
    TYPE(atprop_type), POINTER               :: atprop_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'force_intra_control', &
      routineP = moduleN//':'//routineN

    INTEGER :: first_atom, handle, i, ibend, ibond, ikind, imol, imul, &
      index_a, index_b, index_c, index_d, iopbend, ishell, itorsion, nbends, &
      nbonds, nimptors, nkind, nmol_per_kind, nopbends, nshell, ntorsions, nub
    LOGICAL                                  :: atener, atstress
    REAL(KIND=dp) :: d12, d32, dist, dist1, dist2, energy, fscalar, id12, &
      id32, is32, ism, isn, k2_spring, k4_spring, r2, s32, sm, sn, theta
    REAL(KIND=dp), DIMENSION(3)              :: b12, b32, g1, g2, g3, gt1, &
                                                gt2, gt3, gt4, k1, k2, k3, &
                                                k4, rij, t12, t32, t34, t41, &
                                                t42, t43, tm, tn
    REAL(KIND=dp), DIMENSION(:), POINTER     :: ener_a
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: pv_a
    TYPE(bend_type), POINTER                 :: bend_list( : )
    TYPE(bond_type), POINTER                 :: bond_list( : )
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(impr_type), POINTER                 :: impr_list( : )
    TYPE(molecule_kind_type), POINTER        :: molecule_kind
    TYPE(molecule_type), POINTER             :: molecule
    TYPE(opbend_type), POINTER               :: opbend_list( : )
    TYPE(shell_type), POINTER                :: shell_list( : )
    TYPE(torsion_type), POINTER              :: torsion_list( : )
    TYPE(ub_type), POINTER                   :: ub_list( : )

    CALL timeset ( routineN, handle )
    NULLIFY(logger)
    logger => cp_error_get_logger(error)

    IF ( PRESENT ( f_bond)) f_bond = 0.0_dp
    IF ( PRESENT ( f_bend)) f_bend = 0.0_dp
    IF ( PRESENT ( f_torsion)) f_torsion = 0.0_dp
    IF ( PRESENT ( f_imptor)) f_imptor = 0.0_dp
    IF ( PRESENT ( f_ub)) f_ub = 0.0_dp

    pot_bond =  0.0_dp
    pot_bend =  0.0_dp
    pot_urey_bradley =  0.0_dp
    pot_torsion =  0.0_dp
    pot_imp_torsion =  0.0_dp
    pot_opbend =0.0_dp
    pot_shell =  0.0_dp

    atener = atprop_env%energy
    IF (atener) ener_a => atprop_env%atener
    atstress = atprop_env%stress
    IF (atstress) pv_a => atprop_env%atstress

    nkind = SIZE ( molecule_kind_set )
    MOL: DO ikind = 1, nkind
       nmol_per_kind = local_molecules % n_el ( ikind )

       DO imol = 1, nmol_per_kind
          i = local_molecules % list ( ikind ) % array ( imol )
          molecule => molecule_set ( i )
          molecule_kind => molecule % molecule_kind
          CALL get_molecule_kind ( molecule_kind, nbend = nbends, nbond = nbonds,  &
                                   nimpr = nimptors, nub = nub, ntorsion = ntorsions, &
                                   nopbend = nopbends, nshell=nshell, &
                                   bond_list = bond_list, ub_list = ub_list, &
                                   bend_list = bend_list, torsion_list = torsion_list, &
                                   impr_list = impr_list, opbend_list = opbend_list, &
                                   shell_list=shell_list )

          CALL get_molecule ( molecule, first_atom=first_atom )

          BOND: DO ibond = 1, nbonds
             index_a = bond_list ( ibond ) % a + first_atom - 1
             index_b = bond_list ( ibond ) % b + first_atom - 1
             rij = particle_set ( index_a ) % r - particle_set ( index_b ) % r
             rij = pbc(rij,cell)
             CALL force_bonds(bond_list(ibond)%bond_kind%id_type, rij, &
                              bond_list ( ibond ) % bond_kind % r0, &
                              bond_list ( ibond ) % bond_kind % k, &
                              bond_list ( ibond ) % bond_kind % cs, &
                              energy, fscalar)
             pot_bond = pot_bond + energy
             IF (atener) THEN
                ener_a(index_a) = ener_a(index_a) + 0.5_dp*energy
                ener_a(index_b) = ener_a(index_b) + 0.5_dp*energy
             END IF

             particle_set ( index_a )% f(1) = particle_set ( index_a ) % f(1)  -rij(1) * fscalar
             particle_set ( index_a )% f(2) = particle_set ( index_a ) % f(2)  -rij(2) * fscalar
             particle_set ( index_a )% f(3) = particle_set ( index_a ) % f(3)  -rij(3) * fscalar
             particle_set ( index_b )% f(1) = particle_set ( index_b ) % f(1)  +rij(1) * fscalar
             particle_set ( index_b )% f(2) = particle_set ( index_b ) % f(2)  +rij(2) * fscalar
             particle_set ( index_b )% f(3) = particle_set ( index_b ) % f(3)  +rij(3) * fscalar

             ! computing the pressure tensor
             k2 =  - rij * fscalar
             IF (use_virial) CALL get_pv_bond  ( k2, rij, pv_bond )
             IF (atstress) THEN
               CALL get_pv_bond  ( 0.5_dp*k2, rij, pv_a(:,:,index_a) )
               CALL get_pv_bond  ( 0.5_dp*k2, rij, pv_a(:,:,index_b) )
             END IF

             ! the contribution from the bonds. ONLY FOR DEBUG
             IF ( PRESENT ( f_bond ) ) THEN
                f_bond(1,index_a) = f_bond(1,index_a) - rij(1)*fscalar
                f_bond(2,index_a) = f_bond(2,index_a) - rij(2)*fscalar
                f_bond(3,index_a) = f_bond(3,index_a) - rij(3)*fscalar
                f_bond(1,index_b) = f_bond(1,index_b) + rij(1)*fscalar
                f_bond(2,index_b) = f_bond(2,index_b) + rij(2)*fscalar
                f_bond(3,index_b) = f_bond(3,index_b) + rij(3)*fscalar
             END IF

          END DO BOND

          SHELL: DO ishell=1,nshell
             index_a = shell_list(ishell)%a + first_atom - 1
             index_b = particle_set(index_a)%shell_index
             rij = core_particle_set(index_b)%r - shell_particle_set(index_b)%r
             rij = pbc(rij,cell)
             k2_spring = shell_list(ishell)%shell_kind%k2_spring
             k4_spring = shell_list(ishell)%shell_kind%k4_spring
             r2  = DOT_PRODUCT(rij,rij)
             energy  = 0.5_dp*(k2_spring + k4_spring*r2/12.0_dp)*r2
             fscalar = k2_spring + k4_spring*r2/6.0_dp
             pot_shell = pot_shell + energy
             IF (atener) THEN
                ener_a(index_a) = ener_a(index_a) + energy
             END IF
             core_particle_set(index_b)%f(1) = core_particle_set(index_b)%f(1) - rij(1)*fscalar
             core_particle_set(index_b)%f(2) = core_particle_set(index_b)%f(2) - rij(2)*fscalar
             core_particle_set(index_b)%f(3) = core_particle_set(index_b)%f(3) - rij(3)*fscalar
             shell_particle_set(index_b)%f(1) = shell_particle_set(index_b)%f(1) + rij(1)*fscalar
             shell_particle_set(index_b)%f(2) = shell_particle_set(index_b)%f(2) + rij(2)*fscalar
             shell_particle_set(index_b)%f(3) = shell_particle_set(index_b)%f(3) + rij(3)*fscalar
             ! Compute the pressure tensor, if requested
             IF (use_virial) THEN
                k1 = -rij*fscalar
                CALL get_pv_bond(k1,rij,pv_bond)
             END IF
             IF (atstress) THEN
                CALL get_pv_bond(k1,rij,pv_a(:,:,index_a))
             END IF
          END DO SHELL

          UREY_BRADLEY: DO ibend = 1, nub
             index_a = ub_list ( ibend ) % a + first_atom - 1
             index_b = ub_list ( ibend ) % c + first_atom - 1
             rij = particle_set ( index_a ) % r - particle_set ( index_b ) % r
             rij = pbc(rij,cell)
             CALL force_bonds (ub_list(ibend)%ub_kind%id_type, rij, &
                               ub_list(ibend)%ub_kind%r0, &
                               ub_list(ibend)%ub_kind%k, 0.0_dp, energy, fscalar)
             pot_urey_bradley = pot_urey_bradley + energy
             IF (atener) THEN
                ener_a(index_a) = ener_a(index_a) + 0.5_dp*energy
                ener_a(index_b) = ener_a(index_b) + 0.5_dp*energy
             END IF
             particle_set ( index_a )% f(1) = particle_set ( index_a ) % f(1)  - rij(1) * fscalar
             particle_set ( index_a )% f(2) = particle_set ( index_a ) % f(2)  - rij(2) * fscalar
             particle_set ( index_a )% f(3) = particle_set ( index_a ) % f(3)  - rij(3) * fscalar
             particle_set ( index_b )% f(1) = particle_set ( index_b ) % f(1)  + rij(1) * fscalar
             particle_set ( index_b )% f(2) = particle_set ( index_b ) % f(2)  + rij(2) * fscalar
             particle_set ( index_b )% f(3) = particle_set ( index_b ) % f(3)  + rij(3) * fscalar

             ! computing the pressure tensor
             k2 =  - rij * fscalar
             IF (use_virial) CALL get_pv_bond ( k2,rij,pv_urey_bradley )
             IF (atstress) THEN
               CALL get_pv_bond  ( 0.5_dp*k2, rij, pv_a(:,:,index_a) )
               CALL get_pv_bond  ( 0.5_dp*k2, rij, pv_a(:,:,index_b) )
             END IF

             ! the contribution from the ub. ONLY FOR DEBUG
             IF ( PRESENT ( f_ub ) ) THEN
                f_ub(:,index_a) = f_ub(:,index_a) - rij*fscalar
                f_ub(:,index_b) = f_ub(:,index_b) + rij*fscalar
             END IF

          END DO UREY_BRADLEY

          BEND: DO ibend = 1, nbends
             index_a = bend_list ( ibend ) % a + first_atom - 1
             index_b = bend_list ( ibend ) % b + first_atom - 1
             index_c = bend_list ( ibend ) % c + first_atom - 1
             b12 = particle_set ( index_a ) % r - particle_set ( index_b ) % r
             b32 = particle_set ( index_c ) % r - particle_set ( index_b ) % r
             b12 = pbc(b12,cell)
             b32 = pbc(b32,cell)
             d12 = SQRT(DOT_PRODUCT(b12,b12))
             id12 = 1.0_dp/d12
             d32 = SQRT(DOT_PRODUCT(b32,b32))
             id32 = 1.0_dp/d32
             dist = DOT_PRODUCT(b12,b32)
             theta = (dist*id12*id32)
             IF(theta<-1.0_dp) theta=-1.0_dp
             IF(theta>+1.0_dp) theta=+1.0_dp
             theta = ACOS(theta)
             CALL force_bends( bend_list(ibend)%bend_kind%id_type, &
                               b12,b32,d12,d32,id12,id32,dist,theta, &
                               bend_list(ibend)%bend_kind%theta0, &
                               bend_list(ibend)%bend_kind%k, &
                               bend_list(ibend)%bend_kind%cb, &
                               bend_list(ibend)%bend_kind%r012, &
                               bend_list(ibend)%bend_kind%r032, &
                               bend_list(ibend)%bend_kind%kbs12, &
                               bend_list(ibend)%bend_kind%kbs32, &
                               bend_list(ibend)%bend_kind%kss, &
                               g1,g2,g3,energy,fscalar)
             pot_bend = pot_bend + energy
             IF (atener) THEN
                ener_a(index_a) = ener_a(index_a) + energy/3._dp
                ener_a(index_b) = ener_a(index_b) + energy/3._dp
                ener_a(index_c) = ener_a(index_c) + energy/3._dp
             END IF
             particle_set ( index_a )% f(1) = particle_set ( index_a ) % f(1) + g1(1) * fscalar
             particle_set ( index_a )% f(2) = particle_set ( index_a ) % f(2) + g1(2) * fscalar
             particle_set ( index_a )% f(3) = particle_set ( index_a ) % f(3) + g1(3) * fscalar
             particle_set ( index_b )% f(1) = particle_set ( index_b ) % f(1) + g2(1) * fscalar
             particle_set ( index_b )% f(2) = particle_set ( index_b ) % f(2) + g2(2) * fscalar
             particle_set ( index_b )% f(3) = particle_set ( index_b ) % f(3) + g2(3) * fscalar
             particle_set ( index_c )% f(1) = particle_set ( index_c ) % f(1) + g3(1) * fscalar
             particle_set ( index_c )% f(2) = particle_set ( index_c ) % f(2) + g3(2) * fscalar
             particle_set ( index_c )% f(3) = particle_set ( index_c ) % f(3) + g3(3) * fscalar

             ! computing the pressure tensor
             k1 = fscalar*g1
             k3 = fscalar*g3
             IF (use_virial) CALL get_pv_bend(k1,k3,b12,b32,pv_bend)
             IF (atstress) THEN
               k1 = fscalar*g1/3._dp
               k3 = fscalar*g3/3._dp
               CALL get_pv_bend(k1,k3,b12,b32,pv_a(:,:,index_a))
               CALL get_pv_bend(k1,k3,b12,b32,pv_a(:,:,index_b))
               CALL get_pv_bend(k1,k3,b12,b32,pv_a(:,:,index_c))
             END IF

             ! the contribution from the bends. ONLY FOR DEBUG
             IF ( PRESENT ( f_bend)) THEN
                f_bend(:,index_a) = f_bend(:,index_a) + fscalar*g1
                f_bend(:,index_b) = f_bend(:,index_b) + fscalar*g2
                f_bend(:,index_c) = f_bend(:,index_c) + fscalar*g3
             END IF
          END DO BEND

          TORSION: DO itorsion = 1, ntorsions
             index_a = torsion_list ( itorsion ) % a + first_atom - 1
             index_b = torsion_list ( itorsion ) % b + first_atom - 1
             index_c = torsion_list ( itorsion ) % c + first_atom - 1
             index_d = torsion_list ( itorsion ) % d + first_atom - 1
             t12 = particle_set ( index_a ) % r - particle_set ( index_b ) % r
             t32 = particle_set ( index_c ) % r - particle_set ( index_b ) % r
             t34 = particle_set ( index_c ) % r - particle_set ( index_d ) % r
             t43 = particle_set ( index_d ) % r - particle_set ( index_c ) % r
             t12 = pbc(t12,cell)
             t32 = pbc(t32,cell)
             t34 = pbc(t34,cell)
             t43 = pbc(t43,cell)
             ! t12 x t32
             tm(1) =  t12(2)*t32(3)-t32(2)*t12(3)
             tm(2) = -t12(1)*t32(3)+t32(1)*t12(3)
             tm(3) =  t12(1)*t32(2)-t32(1)*t12(2)
             ! t32 x t34
             tn(1) =  t32(2)*t34(3)-t34(2)*t32(3)
             tn(2) = -t32(1)*t34(3)+t34(1)*t32(3)
             tn(3) =  t32(1)*t34(2)-t34(1)*t32(2)
             sm = SQRT ( DOT_PRODUCT(tm,tm) )
             ism = 1.0_dp/sm
             sn = SQRT ( DOT_PRODUCT(tn,tn) )
             isn = 1.0_dp/sn
             s32 = SQRT ( DOT_PRODUCT(t32,t32) )
             is32 = 1.0_dp/s32
             dist1 = DOT_PRODUCT(t12,t32)
             dist2 = DOT_PRODUCT(t34,t32)
             DO imul = 1, torsion_list(itorsion)%torsion_kind%nmul
                CALL force_torsions(torsion_list(itorsion)%torsion_kind%id_type, &
                                    s32, is32,ism,isn,dist1,dist2,tm,tn, t12, t32,&
                                    torsion_list(itorsion)%torsion_kind%k(imul), &
                                    torsion_list(itorsion)%torsion_kind%phi0(imul), &
                                    torsion_list(itorsion)%torsion_kind%m(imul), &
                                    gt1,gt2,gt3,gt4,energy,fscalar)
                pot_torsion = pot_torsion + energy
                IF (atener) THEN
                   ener_a(index_a) = ener_a(index_a) + energy*0.25_dp
                   ener_a(index_b) = ener_a(index_b) + energy*0.25_dp
                   ener_a(index_c) = ener_a(index_c) + energy*0.25_dp
                   ener_a(index_d) = ener_a(index_d) + energy*0.25_dp
                END IF
                particle_set ( index_a )% f(1) = particle_set ( index_a ) % f(1) + gt1(1) * fscalar
                particle_set ( index_a )% f(2) = particle_set ( index_a ) % f(2) + gt1(2) * fscalar
                particle_set ( index_a )% f(3) = particle_set ( index_a ) % f(3) + gt1(3) * fscalar
                particle_set ( index_b )% f(1) = particle_set ( index_b ) % f(1) + gt2(1) * fscalar
                particle_set ( index_b )% f(2) = particle_set ( index_b ) % f(2) + gt2(2) * fscalar
                particle_set ( index_b )% f(3) = particle_set ( index_b ) % f(3) + gt2(3) * fscalar
                particle_set ( index_c )% f(1) = particle_set ( index_c ) % f(1) + gt3(1) * fscalar
                particle_set ( index_c )% f(2) = particle_set ( index_c ) % f(2) + gt3(2) * fscalar
                particle_set ( index_c )% f(3) = particle_set ( index_c ) % f(3) + gt3(3) * fscalar
                particle_set ( index_d )% f(1) = particle_set ( index_d ) % f(1) + gt4(1) * fscalar
                particle_set ( index_d )% f(2) = particle_set ( index_d ) % f(2) + gt4(2) * fscalar
                particle_set ( index_d )% f(3) = particle_set ( index_d ) % f(3) + gt4(3) * fscalar

                ! computing the pressure tensor
                k1 = fscalar*gt1
                k3 = fscalar*gt3
                k4 = fscalar*gt4
                IF (use_virial) CALL get_pv_torsion(k1,k3,k4,t12,t32,t43,pv_torsion)
                IF (atstress) THEN
                   k1 = 0.25_dp*fscalar*gt1
                   k3 = 0.25_dp*fscalar*gt3
                   k4 = 0.25_dp*fscalar*gt4
                   CALL get_pv_torsion(k1,k3,k4,t12,t32,t43,pv_a(:,:,index_a))
                   CALL get_pv_torsion(k1,k3,k4,t12,t32,t43,pv_a(:,:,index_b))
                   CALL get_pv_torsion(k1,k3,k4,t12,t32,t43,pv_a(:,:,index_c))
                   CALL get_pv_torsion(k1,k3,k4,t12,t32,t43,pv_a(:,:,index_d))
                END IF

                ! the contribution from the torsions. ONLY FOR DEBUG
                IF ( PRESENT ( f_torsion)) THEN
                   f_torsion(:,index_a) = f_torsion(:,index_a) + fscalar*gt1
                   f_torsion(:,index_b) = f_torsion(:,index_b) + fscalar*gt2
                   f_torsion(:,index_c) = f_torsion(:,index_c) + fscalar*gt3
                   f_torsion(:,index_d) = f_torsion(:,index_d) + fscalar*gt4
                END IF
             END DO
          END DO TORSION

          IMP_TORSION: DO itorsion = 1, nimptors
             index_a = impr_list ( itorsion ) % a + first_atom - 1
             index_b = impr_list ( itorsion ) % b + first_atom - 1
             index_c = impr_list ( itorsion ) % c + first_atom - 1
             index_d = impr_list ( itorsion ) % d + first_atom - 1
             t12 = particle_set ( index_a ) % r - particle_set ( index_b ) % r
             t32 = particle_set ( index_c ) % r - particle_set ( index_b ) % r
             t34 = particle_set ( index_c ) % r - particle_set ( index_d ) % r
             t43 = particle_set ( index_d ) % r - particle_set ( index_c ) % r
             t12 = pbc(t12,cell)
             t32 = pbc(t32,cell)
             t34 = pbc(t34,cell)
             t43 = pbc(t43,cell)
             ! t12 x t32
             tm(1) =  t12(2)*t32(3)-t32(2)*t12(3)
             tm(2) = -t12(1)*t32(3)+t32(1)*t12(3)
             tm(3) =  t12(1)*t32(2)-t32(1)*t12(2)
             ! t32 x t34
             tn(1) =  t32(2)*t34(3)-t34(2)*t32(3)
             tn(2) = -t32(1)*t34(3)+t34(1)*t32(3)
             tn(3) =  t32(1)*t34(2)-t34(1)*t32(2)
             sm = SQRT ( DOT_PRODUCT(tm,tm) )
             ism = 1.0_dp/sm
             sn = SQRT ( DOT_PRODUCT(tn,tn) )
             isn = 1.0_dp/sn
             s32 = SQRT ( DOT_PRODUCT(t32,t32) )
             is32 = 1.0_dp/s32
             dist1 = DOT_PRODUCT(t12,t32)
             dist2 = DOT_PRODUCT(t34,t32)
             CALL force_imp_torsions( impr_list(itorsion)%impr_kind%id_type, &
                                      s32, is32,ism,isn,dist1,dist2,tm,tn, t12, t32, &
                                      impr_list(itorsion)%impr_kind%k, &
                                      impr_list(itorsion)%impr_kind%phi0, &
                                      gt1,gt2,gt3,gt4,energy,fscalar)
             pot_imp_torsion = pot_imp_torsion + energy
             IF (atener) THEN
                ener_a(index_a) = ener_a(index_a) + energy*0.25_dp
                ener_a(index_b) = ener_a(index_b) + energy*0.25_dp
                ener_a(index_c) = ener_a(index_c) + energy*0.25_dp
                ener_a(index_d) = ener_a(index_d) + energy*0.25_dp
             END IF
             particle_set ( index_a )% f(1) = particle_set ( index_a ) % f(1) + gt1(1) * fscalar
             particle_set ( index_a )% f(2) = particle_set ( index_a ) % f(2) + gt1(2) * fscalar
             particle_set ( index_a )% f(3) = particle_set ( index_a ) % f(3) + gt1(3) * fscalar
             particle_set ( index_b )% f(1) = particle_set ( index_b ) % f(1) + gt2(1) * fscalar
             particle_set ( index_b )% f(2) = particle_set ( index_b ) % f(2) + gt2(2) * fscalar
             particle_set ( index_b )% f(3) = particle_set ( index_b ) % f(3) + gt2(3) * fscalar
             particle_set ( index_c )% f(1) = particle_set ( index_c ) % f(1) + gt3(1) * fscalar
             particle_set ( index_c )% f(2) = particle_set ( index_c ) % f(2) + gt3(2) * fscalar
             particle_set ( index_c )% f(3) = particle_set ( index_c ) % f(3) + gt3(3) * fscalar
             particle_set ( index_d )% f(1) = particle_set ( index_d ) % f(1) + gt4(1) * fscalar
             particle_set ( index_d )% f(2) = particle_set ( index_d ) % f(2) + gt4(2) * fscalar
             particle_set ( index_d )% f(3) = particle_set ( index_d ) % f(3) + gt4(3) * fscalar

             ! computing the pressure tensor
             k1 = fscalar*gt1
             k3 = fscalar*gt3
             k4 = fscalar*gt4
             IF (use_virial) CALL get_pv_torsion(k1,k3,k4,t12,t32,t43,pv_imp_torsion)
             IF (atstress) THEN
                k1 = 0.25_dp*fscalar*gt1
                k3 = 0.25_dp*fscalar*gt3
                k4 = 0.25_dp*fscalar*gt4
                CALL get_pv_torsion(k1,k3,k4,t12,t32,t43,pv_a(:,:,index_a))
                CALL get_pv_torsion(k1,k3,k4,t12,t32,t43,pv_a(:,:,index_b))
                CALL get_pv_torsion(k1,k3,k4,t12,t32,t43,pv_a(:,:,index_c))
                CALL get_pv_torsion(k1,k3,k4,t12,t32,t43,pv_a(:,:,index_d))
             END IF

             ! the contribution from the torsions. ONLY FOR DEBUG
             IF ( PRESENT ( f_imptor)) THEN
                f_imptor(:,index_a) = f_imptor(:,index_a) + fscalar*gt1
                f_imptor(:,index_b) = f_imptor(:,index_b) + fscalar*gt2
                f_imptor(:,index_c) = f_imptor(:,index_c) + fscalar*gt3
                f_imptor(:,index_d) = f_imptor(:,index_d) + fscalar*gt4
             END IF
          END DO IMP_TORSION

          OPBEND: DO iopbend = 1, nopbends
             index_a = opbend_list ( iopbend ) % a + first_atom - 1
             index_b = opbend_list ( iopbend ) % b + first_atom - 1
             index_c = opbend_list ( iopbend ) % c + first_atom - 1
             index_d = opbend_list ( iopbend ) % d + first_atom - 1

             t12 = particle_set ( index_a ) % r - particle_set ( index_b ) % r
             t32 = particle_set ( index_c ) % r - particle_set ( index_b ) % r
             t34 = particle_set ( index_c ) % r - particle_set ( index_d ) % r
             t43 = particle_set ( index_d ) % r - particle_set ( index_c ) % r
             t41 = particle_set ( index_d ) % r - particle_set ( index_a ) % r
             t42 = pbc(t41+t12,cell)
             t12 = pbc(t12,cell)
             t32 = pbc(t32,cell)
             t41 = pbc(t41,cell)
             t43 = pbc(t43,cell)
             ! tm = t32 x t12
             tm(1) =  t32(2)*t12(3)-t12(2)*t32(3)
             tm(2) = -t32(1)*t12(3)+t12(1)*t32(3)
             tm(3) =  t32(1)*t12(2)-t12(1)*t32(2)
             ! tn = - t41 x tm
             tn(1) = -( t41(2)*tm(3)-tm(2)*t41(3))
             tn(2) = -(-t41(1)*tm(3)+tm(1)*t41(3))
             tn(3) = -( t41(1)*tm(2)-tm(1)*t41(2))
             sm = SQRT ( DOT_PRODUCT(tm,tm) )
             ism = 1.0_dp/sm
             sn = SQRT ( DOT_PRODUCT(tn,tn) )
             isn = 1.0_dp/sn
             s32 = SQRT ( DOT_PRODUCT(t32,t32) )
             is32 = 1.0_dp/s32
             CALL force_opbends( opbend_list(iopbend)%opbend_kind%id_type, &
                                      s32, is32,ism,isn,tm,tn, t41, t42, t43, &
                                      opbend_list(iopbend)%opbend_kind%k, &
                                      opbend_list(iopbend)%opbend_kind%phi0, &
                                      gt1,gt2,gt3,gt4,energy,fscalar)
             pot_opbend = pot_opbend + energy
             IF (atener) THEN
                ener_a(index_a) = ener_a(index_a) + energy*0.25_dp
                ener_a(index_b) = ener_a(index_b) + energy*0.25_dp
                ener_a(index_c) = ener_a(index_c) + energy*0.25_dp
                ener_a(index_d) = ener_a(index_d) + energy*0.25_dp
             END IF
             particle_set ( index_a )% f(1) = particle_set ( index_a ) % f(1) + gt1(1) * fscalar
             particle_set ( index_a )% f(2) = particle_set ( index_a ) % f(2) + gt1(2) * fscalar
             particle_set ( index_a )% f(3) = particle_set ( index_a ) % f(3) + gt1(3) * fscalar
             particle_set ( index_b )% f(1) = particle_set ( index_b ) % f(1) + gt2(1) * fscalar
             particle_set ( index_b )% f(2) = particle_set ( index_b ) % f(2) + gt2(2) * fscalar
             particle_set ( index_b )% f(3) = particle_set ( index_b ) % f(3) + gt2(3) * fscalar
             particle_set ( index_c )% f(1) = particle_set ( index_c ) % f(1) + gt3(1) * fscalar
             particle_set ( index_c )% f(2) = particle_set ( index_c ) % f(2) + gt3(2) * fscalar
             particle_set ( index_c )% f(3) = particle_set ( index_c ) % f(3) + gt3(3) * fscalar
             particle_set ( index_d )% f(1) = particle_set ( index_d ) % f(1) + gt4(1) * fscalar
             particle_set ( index_d )% f(2) = particle_set ( index_d ) % f(2) + gt4(2) * fscalar
             particle_set ( index_d )% f(3) = particle_set ( index_d ) % f(3) + gt4(3) * fscalar

             ! computing the pressure tensor
             k1 = fscalar*gt1
             k3 = fscalar*gt3
             k4 = fscalar*gt4

             IF (use_virial) CALL get_pv_torsion(k1,k3,k4,t12,t32,t43,pv_opbend)
             IF (atstress) THEN
                k1 = 0.25_dp*fscalar*gt1
                k3 = 0.25_dp*fscalar*gt3
                k4 = 0.25_dp*fscalar*gt4
                CALL get_pv_torsion(k1,k3,k4,t12,t32,t43,pv_a(:,:,index_a))
                CALL get_pv_torsion(k1,k3,k4,t12,t32,t43,pv_a(:,:,index_b))
                CALL get_pv_torsion(k1,k3,k4,t12,t32,t43,pv_a(:,:,index_c))
                CALL get_pv_torsion(k1,k3,k4,t12,t32,t43,pv_a(:,:,index_d))
             END IF

             ! the contribution from the opbends. ONLY FOR DEBUG
             IF ( PRESENT ( f_opbend)) THEN
                f_opbend(:,index_a) = f_opbend(:,index_a) + fscalar*gt1
                f_opbend(:,index_b) = f_opbend(:,index_b) + fscalar*gt2
                f_opbend(:,index_c) = f_opbend(:,index_c) + fscalar*gt3
                f_opbend(:,index_d) = f_opbend(:,index_d) + fscalar*gt4
             END IF
          END DO OPBEND
       END DO
    END DO MOL

    CALL timestop(handle)

  END SUBROUTINE force_intra_control

END MODULE fist_intra_force

