!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of dispersion using pair potentials
!> \author Johann Pototschnig
! **************************************************************************************************
MODULE qs_dispersion_d4

   USE atomic_kind_types, ONLY: atomic_kind_type, &
                                get_atomic_kind, &
                                get_atomic_kind_set
   USE cell_types, ONLY: cell_type, &
                         plane_distance, &
                         pbc, &
                         get_cell
   USE qs_environment_types, ONLY: get_qs_env, &
                                   qs_environment_type
   USE qs_force_types, ONLY: qs_force_type
   USE qs_kind_types, ONLY: get_qs_kind, &
                            qs_kind_type, &
                            set_qs_kind
   USE qs_neighbor_list_types, ONLY: get_iterator_info, &
                                     neighbor_list_iterate, &
                                     neighbor_list_iterator_create, &
                                     neighbor_list_iterator_p_type, &
                                     neighbor_list_iterator_release, &
                                     neighbor_list_set_p_type
   USE virial_types, ONLY: virial_type
   USE kinds, ONLY: dp
   USE particle_types, ONLY: particle_type
   USE qs_dispersion_types, ONLY: qs_dispersion_type
   USE qs_dispersion_utils, ONLY: cellhash
   USE message_passing, ONLY: mp_para_env_type

#if defined(__DFTD4)
!&<
   USE dftd4,                           ONLY: d4_model, &
                                              damping_param, &
                                              get_dispersion, &
                                              get_rational_damping, &
                                              new, &
                                              new_d4_model, &
                                              realspace_cutoff, &
                                              structure_type, &
                                              rational_damping_param, &
                                              get_coordination_number, &
                                              get_lattice_points
   USE dftd4_charge,                    ONLY: get_charges
!&>
#endif
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_dispersion_d4'

   PUBLIC :: calculate_dispersion_d4_pairpot

! **************************************************************************************************

CONTAINS

#if defined(__DFTD4)
! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param dispersion_env ...
!> \param evdw ...
!> \param calculate_forces ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE calculate_dispersion_d4_pairpot(qs_env, dispersion_env, evdw, calculate_forces, unit_nr)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_dispersion_type), INTENT(IN), POINTER      :: dispersion_env
      REAL(KIND=dp), INTENT(INOUT)                       :: evdw
      LOGICAL, INTENT(IN)                                :: calculate_forces
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_dispersion_d4_pairpot'
      LOGICAL, PARAMETER                                 :: debug = .FALSE.

      INTEGER                                            :: atoma, handle, i, iatom, ikind, mref, &
                                                            natom
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, atomtype, kind_of
      INTEGER, DIMENSION(3)                              :: periodic
      LOGICAL                                            :: grad, use_virial
      LOGICAL, DIMENSION(3)                              :: lperiod
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: gradient, xyz
      REAL(KIND=dp), DIMENSION(3, 3)                     :: stress
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(virial_type), POINTER                         :: virial

      CLASS(damping_param), ALLOCATABLE                  :: param
      TYPE(d4_model)                                     :: disp
      TYPE(structure_type)                               :: mol
      TYPE(realspace_cutoff)                             :: cutoff
      REAL(dp), ALLOCATABLE :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
      REAL(dp), ALLOCATABLE :: q(:), dqdr(:, :, :), dqdL(:, :, :)
      REAL(dp), ALLOCATABLE :: gwvec(:, :), gwdcn(:, :), gwdq(:, :)
      REAL(dp), ALLOCATABLE :: c6(:, :), dc6dcn(:, :), dc6dq(:, :)
      REAL(dp), ALLOCATABLE :: dEdcn(:), dEdq(:), energies(:)
      REAL(dp), ALLOCATABLE :: lattr(:, :)

      MARK_USED(unit_nr)

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, particle_set=particle_set, atomic_kind_set=atomic_kind_set, &
                      cell=cell, force=force, virial=virial, para_env=para_env)
      CALL get_atomic_kind_set(atomic_kind_set, atom_of_kind=atom_of_kind, kind_of=kind_of)

      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)

      !get information about particles
      natom = SIZE(particle_set)
      ALLOCATE (xyz(3, natom), atomtype(natom))
      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set)
      DO iatom = 1, natom
         xyz(:, iatom) = particle_set(iatom)%r(:)
         ikind = kind_of(iatom)
         CALL get_qs_kind(qs_kind_set(ikind), zatom=atomtype(iatom))
      END DO

      !get information about cell / lattice
      CALL get_cell(cell=cell, periodic=periodic)
      lperiod(1) = periodic(1) == 1
      lperiod(2) = periodic(2) == 1
      lperiod(3) = periodic(3) == 1

      !prepare for the call to the dispersion function
      CALL new(mol, atomtype, xyz, lattice=cell%hmat, periodic=lperiod)
      CALL new_d4_model(disp, mol)
      CALL get_rational_damping(dispersion_env%ref_functional, param, s9=dispersion_env%s9)

      SELECT TYPE (param)
      TYPE is (rational_damping_param)
         dispersion_env%s6 = param%s6
         dispersion_env%s8 = param%s8
         dispersion_env%a1 = param%a1
         dispersion_env%a2 = param%a2
         dispersion_env%alp = param%alp
      END SELECT

      ! Coordination number cutoff
      cutoff%cn = dispersion_env%rc_cn
      ! Two-body interaction cutoff
      cutoff%disp2 = dispersion_env%rc_d4*2._dp
      ! Three-body interaction cutoff
      cutoff%disp3 = dispersion_env%rc_disp*2._dp

      IF (calculate_forces .OR. use_virial .OR. dispersion_env%verbose) THEN
         grad = .TRUE.
      ELSE
         grad = .FALSE.
      END IF

      IF (dispersion_env%d4_reference_code) THEN

         !> Wrapper to handle the evaluation of dispersion energy and derivatives
         IF (.NOT. dispersion_env%doabc) THEN
            CPABORT("Using D4_REFERENCE_CODE enforces calculation of C9 term.")
         END IF
         IF (grad) THEN
            ALLOCATE (gradient(3, natom))
            CALL get_dispersion(mol, disp, param, cutoff, evdw, gradient, stress)
            IF (use_virial) THEN
               virial%pv_virial = virial%pv_virial - stress/para_env%num_pe
            END IF
            DO iatom = 1, natom
               ikind = kind_of(iatom)
               atoma = atom_of_kind(iatom)
               force(ikind)%dispersion(:, atoma) = &
                  force(ikind)%dispersion(:, atoma) + gradient(:, iatom)/para_env%num_pe
            END DO
            DEALLOCATE (gradient)
         ELSE
            CALL get_dispersion(mol, disp, param, cutoff, evdw)
         END IF
         !dispersion energy is computed by every MPI process
         evdw = evdw/para_env%num_pe

      ELSE

         mref = MAXVAL(disp%ref)

         IF (grad) ALLOCATE (gradient(3, natom))

         ! Coordination numbers
         ALLOCATE (cn(natom))
         IF (grad) ALLOCATE (dcndr(3, natom, natom), dcndL(3, 3, natom))
         CALL get_lattice_points(mol%periodic, mol%lattice, cutoff%cn, lattr)
         CALL get_coordination_number(mol, lattr, cutoff%cn, disp%rcov, disp%en, &
            & cn, dcndr, dcndL)

         ! EEQ charges
         ALLOCATE (q(natom))
         IF (grad) ALLOCATE (dqdr(3, natom, natom), dqdL(3, 3, natom))
         CALL get_charges(mol, q, dqdr, dqdL)

         ! Weights for C6 calculation
         ALLOCATE (gwvec(mref, natom))
         IF (grad) ALLOCATE (gwdcn(mref, natom), gwdq(mref, natom))
         CALL disp%weight_references(mol, cn, q, gwvec, gwdcn, gwdq)

         ! C6 calculation
         ALLOCATE (c6(mol%nat, mol%nat))
         IF (grad) ALLOCATE (dc6dcn(natom, natom), dc6dq(natom, natom))
         CALL disp%get_atomic_c6(mol, gwvec, gwdcn, gwdq, c6, dc6dcn, dc6dq)

         ALLOCATE (energies(natom))
         energies(:) = 0.0_dp
         IF (grad) THEN
            ALLOCATE (dEdcn(natom), dEdq(natom))
            dEdcn(:) = 0.0_dp
            dEdq(:) = 0.0_dp
            gradient(:, :) = 0.0_dp
            stress(:, :) = 0.0_dp
         END IF
         IF (debug) THEN
            CALL get_lattice_points(mol%periodic, mol%lattice, cutoff%disp2, lattr)
            CALL param%get_dispersion2(mol, lattr, cutoff%disp2, disp%r4r2, &
               & c6, dc6dcn, dc6dq, energies, dEdcn, dEdq, gradient, stress)
         ELSE
            CALL dispersion_2b(qs_env, dispersion_env, cutoff%disp2, disp%r4r2, &
                               c6, dc6dcn, dc6dq, energies, dEdcn, dEdq, &
                               grad, gradient, stress)
         END IF
         IF (grad) THEN
            DO i = 1, 3
               gradient(i, :) = gradient(i, :) + MATMUL(dqdr(i, :, :), dEdq(:))
               stress(i, :) = stress(i, :) + MATMUL(dqdL(i, :, :), dEdq(:))
            END DO
         END IF

         IF (dispersion_env%doabc) THEN
            q(:) = 0.0_dp
            CALL disp%weight_references(mol, cn, q, gwvec, gwdcn, gwdq)
            CALL disp%get_atomic_c6(mol, gwvec, gwdcn, gwdq, c6, dc6dcn, dc6dq)

            IF (debug) THEN
               CALL get_lattice_points(mol%periodic, mol%lattice, cutoff%disp3, lattr)
               CALL param%get_dispersion3(mol, lattr, cutoff%disp3, disp%r4r2, &
                  & c6, dc6dcn, dc6dq, energies, dEdcn, dEdq, gradient, stress)
            ELSE
               CALL dispersion_3b(qs_env, dispersion_env, cutoff%disp3, disp%r4r2, &
                                  c6, dc6dcn, dc6dq, energies, dEdcn, dEdq, &
                                  grad, gradient, stress)
            END IF
         END IF
         IF (grad) THEN
            DO i = 1, 3
               gradient(i, :) = gradient(i, :) + MATMUL(dcndr(i, :, :), dEdcn(:))
               stress(i, :) = stress(i, :) + MATMUL(dcndL(i, :, :), dEdcn(:))
            END DO
         END IF

         evdw = SUM(energies)
         IF (debug) evdw = evdw/para_env%num_pe

         IF (use_virial) THEN
            IF (debug) THEN
               virial%pv_virial = virial%pv_virial - stress/para_env%num_pe
            ELSE
               virial%pv_virial = virial%pv_virial - stress
            END IF
         END IF
         IF (grad) THEN
            DO iatom = 1, natom
               ikind = kind_of(iatom)
               atoma = atom_of_kind(iatom)
               force(ikind)%dispersion(:, atoma) = &
                  force(ikind)%dispersion(:, atoma) + gradient(:, iatom)
            END DO
            DEALLOCATE (gradient)
         END IF

      END IF

      DEALLOCATE (xyz, atomtype)

      CALL timestop(handle)

   END SUBROUTINE calculate_dispersion_d4_pairpot

! **************************************************************************************************

#else

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param dispersion_env ...
!> \param evdw ...
!> \param calculate_forces ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE calculate_dispersion_d4_pairpot(qs_env, dispersion_env, evdw, calculate_forces, unit_nr)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_dispersion_type), INTENT(IN), POINTER      :: dispersion_env
      REAL(KIND=dp), INTENT(INOUT)                       :: evdw
      LOGICAL, INTENT(IN)                                :: calculate_forces
      INTEGER, INTENT(IN)                                :: unit_nr

      MARK_USED(qs_env)
      MARK_USED(dispersion_env)
      MARK_USED(evdw)
      MARK_USED(calculate_forces)
      MARK_USED(unit_nr)

      CPABORT("CP2K build without DFTD4")

   END SUBROUTINE calculate_dispersion_d4_pairpot

#endif

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param dispersion_env ...
!> \param cutoff ...
!> \param r4r2 ...
!> \param c6 ...
!> \param dc6dcn ...
!> \param dc6dq ...
!> \param energies ...
!> \param dEdcn ...
!> \param dEdq ...
!> \param calculate_forces ...
!> \param gradient ...
!> \param stress ...
! **************************************************************************************************
   SUBROUTINE dispersion_2b(qs_env, dispersion_env, cutoff, r4r2, &
                            c6, dc6dcn, dc6dq, energies, dEdcn, dEdq, &
                            calculate_forces, gradient, stress)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_dispersion_type), POINTER                  :: dispersion_env
      REAL(KIND=dp), INTENT(IN)                          :: cutoff
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: r4r2
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: c6, dc6dcn, dc6dq
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: energies, dEdcn, dEdq
      LOGICAL, INTENT(IN)                                :: calculate_forces
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: gradient, stress

      INTEGER                                            :: iatom, ikind, jatom, jkind, mepos, num_pe
      REAL(KINd=dp)                                      :: a1, a2, c6ij, d6, d8, dE, dr, dr2, &
                                                            edisp, fac, gdisp, r0ij, rrij, s6, s8, &
                                                            t6, t8
      REAL(KINd=dp), DIMENSION(3)                        :: dG, rij
      REAL(KINd=dp), DIMENSION(3, 3)                     :: dS
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_vdw

      MARK_USED(qs_env)

      a1 = dispersion_env%a1
      a2 = dispersion_env%a2
      s6 = dispersion_env%s6
      s8 = dispersion_env%s8

      sab_vdw => dispersion_env%sab_vdw

      num_pe = 1
      CALL neighbor_list_iterator_create(nl_iterator, sab_vdw, nthread=num_pe)

      mepos = 0
      DO WHILE (neighbor_list_iterate(nl_iterator, mepos=mepos) == 0)
         CALL get_iterator_info(nl_iterator, mepos=mepos, ikind=ikind, jkind=jkind, &
                                iatom=iatom, jatom=jatom, r=rij)
         ! vdW potential
         dr2 = SUM(rij(:)**2)
         dr = SQRT(dr2)
         IF (dr <= cutoff .AND. dr > 0.001_dp) THEN
            rrij = 3._dp*r4r2(ikind)*r4r2(jkind)
            r0ij = a1*SQRT(rrij) + a2
            c6ij = c6(jatom, iatom)
            fac = 1._dp
            IF (iatom == jatom) fac = 0.5_dp
            t6 = 1.0_dp/(dr2**3 + r0ij**6)
            t8 = 1.0_dp/(dr2**4 + r0ij**8)

            edisp = (s6*t6 + s8*rrij*t8)*fac
            dE = -c6ij*edisp
            energies(iatom) = energies(iatom) + dE*0.5_dp
            energies(jatom) = energies(jatom) + dE*0.5_dp

            IF (calculate_forces) THEN
               d6 = -6.0_dp*dr2**2*t6**2
               d8 = -8.0_dp*dr2**3*t8**2
               gdisp = (s6*d6 + s8*rrij*d8)*fac
               dG(:) = -c6ij*gdisp*rij(:)
               gradient(:, iatom) = gradient(:, iatom) - dG
               gradient(:, jatom) = gradient(:, jatom) + dG
               dS(:, :) = SPREAD(dG, 1, 3)*SPREAD(rij, 2, 3)
               dEdcn(iatom) = dEdcn(iatom) - dc6dcn(iatom, jatom)*edisp
               dEdq(iatom) = dEdq(iatom) - dc6dq(iatom, jatom)*edisp
               dEdcn(jatom) = dEdcn(jatom) - dc6dcn(jatom, iatom)*edisp
               dEdq(jatom) = dEdq(jatom) - dc6dq(jatom, iatom)*edisp
               stress(:, :) = stress + dS
            END IF
         END IF
      END DO

      CALL neighbor_list_iterator_release(nl_iterator)

   END SUBROUTINE dispersion_2b

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param dispersion_env ...
!> \param cutoff ...
!> \param r4r2 ...
!> \param c6 ...
!> \param dc6dcn ...
!> \param dc6dq ...
!> \param energies ...
!> \param dEdcn ...
!> \param dEdq ...
!> \param calculate_forces ...
!> \param gradient ...
!> \param stress ...
! **************************************************************************************************
   SUBROUTINE dispersion_3b(qs_env, dispersion_env, cutoff, r4r2, &
                            c6, dc6dcn, dc6dq, energies, dEdcn, dEdq, &
                            calculate_forces, gradient, stress)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_dispersion_type), POINTER                  :: dispersion_env
      REAL(KIND=dp), INTENT(IN)                          :: cutoff
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: r4r2
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: c6, dc6dcn, dc6dq
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: energies, dEdcn, dEdq
      LOGICAL, INTENT(IN)                                :: calculate_forces
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: gradient, stress

      INTEGER                                            :: hashb, hashc, iatom, icx, icy, icz, &
                                                            ikind, jatom, jkind, katom, kkind, &
                                                            kstart, mepos, natom, num_pe
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of
      INTEGER, DIMENSION(3)                              :: cell_b, cell_c, ncell, periodic
      LOGICAL                                            :: is000
      REAL(KINd=dp)                                      :: a1, a2, alp, ang, c6ij, c6ik, c6jk, c9, &
                                                            cutoff2, dang, dE, dfdmp, fac, fac0, &
                                                            fdmp, r0, r0ij, r0ik, r0jk, r1, r2, &
                                                            r2ij, r2ik, r2jk, r3, r5, rr, s6, s8, &
                                                            s9
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: rcpbc
      REAL(KINd=dp), DIMENSION(3)                        :: dGij, dGik, dGjk, ra, rb, rb0, rc, rc0, &
                                                            rij, sab_max, vij, vik, vjk
      REAL(KINd=dp), DIMENSION(3, 3)                     :: dS
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_vdw
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL get_qs_env(qs_env=qs_env, natom=natom, cell=cell, &
                      atomic_kind_set=atomic_kind_set, particle_set=particle_set)

      ALLOCATE (rcpbc(3, natom))
      DO iatom = 1, natom
         rcpbc(:, iatom) = pbc(particle_set(iatom)%r(:), cell)
      END DO
      CALL get_atomic_kind_set(atomic_kind_set, kind_of=kind_of)

      a1 = dispersion_env%a1
      a2 = dispersion_env%a2
      s6 = dispersion_env%s6
      s8 = dispersion_env%s8
      s9 = dispersion_env%s9
      alp = dispersion_env%alp

      cutoff2 = cutoff**2
      CALL get_cell(cell=cell, periodic=periodic)
      sab_max(1) = cutoff/plane_distance(1, 0, 0, cell)
      sab_max(2) = cutoff/plane_distance(0, 1, 0, cell)
      sab_max(3) = cutoff/plane_distance(0, 0, 1, cell)
      ncell(:) = (INT(sab_max(:)) + 1)*periodic(:)

      sab_vdw => dispersion_env%sab_vdw

      num_pe = 1
      CALL neighbor_list_iterator_create(nl_iterator, sab_vdw, nthread=num_pe)

      mepos = 0
      DO WHILE (neighbor_list_iterate(nl_iterator, mepos=mepos) == 0)
         CALL get_iterator_info(nl_iterator, mepos=mepos, ikind=ikind, jkind=jkind, iatom=iatom, jatom=jatom, r=rij)

         r2ij = SUM(rij(:)**2)
         c6ij = c6(jatom, iatom)
         r0ij = a1*SQRT(3._dp*r4r2(jkind)*r4r2(ikind)) + a2
         IF (r2ij <= cutoff2 .AND. r2ij > EPSILON(1._dp)) THEN
            CALL get_iterator_info(nl_iterator, cell=cell_b)
            hashb = cellhash(cell_b, ncell)
            is000 = (ALL(cell_b == 0))
            rb0(:) = MATMUL(cell%hmat, cell_b)
            ra(:) = pbc(particle_set(iatom)%r(:), cell)
            rb(:) = pbc(particle_set(jatom)%r(:), cell) + rb0
            DO icx = -ncell(1), ncell(1)
               DO icy = -ncell(2), ncell(2)
                  DO icz = -ncell(3), ncell(3)
                     cell_c(1) = icx
                     cell_c(2) = icy
                     cell_c(3) = icz
                     hashc = cellhash(cell_c, ncell)
                     IF (is000 .AND. (ALL(cell_c == 0))) THEN
                        ! CASE 1: all atoms in (000), use only ordered triples
                        kstart = MAX(jatom + 1, iatom + 1)
                        fac0 = 1.0_dp
                     ELSE IF (is000) THEN
                        ! CASE 2: AB in (000), C in other cell
                        !         This case covers also all instances with BC in same
                        !         cell not (000)
                        kstart = 1
                        fac0 = 1.0_dp
                     ELSE
                        ! These are case 2 again, cycle
                        IF (hashc == hashb) CYCLE
                        IF (ALL(cell_c == 0)) CYCLE
                        ! CASE 3: A in (000) and B and C in different cells
                        kstart = 1
                        fac0 = 1.0_dp/3.0_dp
                     END IF
                     !
                     rc0 = MATMUL(cell%hmat, cell_c)
                     !
                     DO katom = kstart, natom
                        kkind = kind_of(katom)
                        c6ik = c6(katom, iatom)
                        c6jk = c6(katom, jatom)
                        c9 = -s9*SQRT(ABS(c6ij*c6ik*c6jk))
                        r0ik = a1*SQRT(3._dp*r4r2(kkind)*r4r2(ikind)) + a2
                        r0jk = a1*SQRT(3._dp*r4r2(kkind)*r4r2(jkind)) + a2
                        r0 = r0ij*r0ik*r0jk
                        rc(:) = rcpbc(:, katom) + rc0(:)
                        r2jk = SUM((rb - rc)**2)
                        IF (r2jk >= cutoff2 .OR. r2jk < EPSILON(1._dp)) CYCLE
                        r2ik = SUM((rc - ra)**2)
                        IF (r2ik >= cutoff2 .OR. r2ik < EPSILON(1._dp)) CYCLE
                        r2 = r2ij*r2jk*r2ik
                        r1 = SQRT(r2)
                        r3 = r2*r1
                        r5 = r3*r2
                        fdmp = 1.0_dp/(1.0_dp + 6.0_dp*(r0/r1)**(alp/3.0_dp))
                        ang = 0.375_dp*(r2ij + r2jk - r2ik)*(r2ij - r2jk + r2ik)* &
                              (-r2ij + r2jk + r2ik)/r5 + 1.0_dp/r3

                        ! avoid double counting!
                        fac = 1._dp
                        IF (iatom == jatom .OR. iatom == katom .OR. jatom == katom) fac = 0.5_dp
                        IF (iatom == jatom .AND. iatom == katom) fac = 1._dp/3._dp
                        fac = fac*fac0

                        rr = ang*fdmp
                        dE = rr*c9*fac
                        energies(iatom) = energies(iatom) - dE/3._dp
                        energies(jatom) = energies(jatom) - dE/3._dp
                        energies(katom) = energies(katom) - dE/3._dp

                        IF (calculate_forces) THEN
                           vij = rb - ra
                           vjk = rc - rb
                           vik = rc - ra

                           dfdmp = -2.0_dp*alp*(r0/r1)**(alp/3.0_dp)*fdmp**2

                           ! d/drij
                           dang = -0.375_dp*(r2ij**3 + r2ij**2*(r2jk + r2ik) &
                                             + r2ij*(3.0_dp*r2jk**2 + 2.0_dp*r2jk*r2ik &
                                                     + 3.0_dp*r2ik**2) &
                                             - 5.0_dp*(r2jk - r2ik)**2*(r2jk + r2ik))/r5
                           dGij(:) = c9*(-dang*fdmp + ang*dfdmp)/r2ij*vij

                           ! d/drik
                           dang = -0.375_dp*(r2ik**3 + r2ik**2*(r2jk + r2ij) &
                                             + r2ik*(3.0_dp*r2jk**2 + 2.0_dp*r2jk*r2ij &
                                                     + 3.0_dp*r2ij**2) &
                                             - 5.0_dp*(r2jk - r2ij)**2*(r2jk + r2ij))/r5
                           dGik(:) = c9*(-dang*fdmp + ang*dfdmp)/r2ik*vik

                           ! d/drjk
                           dang = -0.375_dp*(r2jk**3 + r2jk**2*(r2ik + r2ij) &
                                             + r2jk*(3.0_dp*r2ik**2 + 2.0_dp*r2ik*r2ij &
                                                     + 3.0_dp*r2ij**2) &
                                             - 5.0_dp*(r2ik - r2ij)**2*(r2ik + r2ij))/r5
                           dGjk(:) = c9*(-dang*fdmp + ang*dfdmp)/r2jk*vjk

                           gradient(:, iatom) = gradient(:, iatom) - dGij - dGik
                           gradient(:, jatom) = gradient(:, jatom) + dGij - dGjk
                           gradient(:, katom) = gradient(:, katom) + dGik + dGjk

                           dS(:, :) = SPREAD(dGij, 1, 3)*SPREAD(vij, 2, 3) &
                                      + SPREAD(dGik, 1, 3)*SPREAD(vik, 2, 3) &
                                      + SPREAD(dGjk, 1, 3)*SPREAD(vjk, 2, 3)

                           stress(:, :) = stress + dS*fac

                           dEdcn(iatom) = dEdcn(iatom) - dE*0.5_dp &
                                         & *(dc6dcn(iatom, jatom)/c6ij + dc6dcn(iatom, katom)/c6ik)
                           dEdcn(jatom) = dEdcn(jatom) - dE*0.5_dp &
                                         & *(dc6dcn(jatom, iatom)/c6ij + dc6dcn(jatom, katom)/c6jk)
                           dEdcn(katom) = dEdcn(katom) - dE*0.5_dp &
                                         & *(dc6dcn(katom, iatom)/c6ik + dc6dcn(katom, jatom)/c6jk)

                           dEdq(iatom) = dEdq(iatom) - dE*0.5_dp &
                                        & *(dc6dq(iatom, jatom)/c6ij + dc6dq(iatom, katom)/c6ik)
                           dEdq(jatom) = dEdq(jatom) - dE*0.5_dp &
                                        & *(dc6dq(jatom, iatom)/c6ij + dc6dq(jatom, katom)/c6jk)
                           dEdq(katom) = dEdq(katom) - dE*0.5_dp &
                                        & *(dc6dq(katom, iatom)/c6ik + dc6dq(katom, jatom)/c6jk)

                        END IF

                     END DO
                  END DO
               END DO
            END DO
         END IF
      END DO

      CALL neighbor_list_iterator_release(nl_iterator)

      DEALLOCATE (rcpbc)

   END SUBROUTINE dispersion_3b

END MODULE qs_dispersion_d4
