!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of gCP pair potentials
!> \author JGH
! **************************************************************************************************
MODULE qs_gcp_method
   USE ai_overlap,                      ONLY: overlap_ab
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE atprop_types,                    ONLY: atprop_array_init,&
                                              atprop_type
   USE cell_types,                      ONLY: cell_type
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_sum
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: kcalmol
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_gcp_types,                    ONLY: qs_gcp_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
   USE virial_methods,                  ONLY: virial_pair_force
   USE virial_types,                    ONLY: virial_type
#include "base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_gcp_method'

   PUBLIC :: calculate_gcp_pairpot

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param gcp_env ...
!> \param energy ...
!> \param calculate_forces ...
!> \note
!> \note energy_correction_type: also add gcp_env and egcp to the type
!> \note
! **************************************************************************************************
   SUBROUTINE calculate_gcp_pairpot(qs_env, gcp_env, energy, calculate_forces)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_gcp_type), POINTER                         :: gcp_env
      REAL(KIND=dp), INTENT(OUT)                         :: energy
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_gcp_pairpot'

      INTEGER                                            :: atom_a, atom_b, handle, i, iatom, ikind, &
                                                            jatom, jkind, mepos, natom, nkind, &
                                                            nsto, unit_nr
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, kind_of, ngcpat
      LOGICAL                                            :: atenergy, atstress, use_virial, verbose
      REAL(KIND=dp)                                      :: eama, eamb, egcp, expab, fac, fda, fdb, &
                                                            gnorm, nvirta, nvirtb, rcc, sint, sqa, &
                                                            sqb
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: egcpat
      REAL(KIND=dp), DIMENSION(3)                        :: dsint, fdij, rij
      REAL(KIND=dp), DIMENSION(3, 3)                     :: dvirial
      REAL(KIND=dp), DIMENSION(6)                        :: cla, clb, rcut, zeta, zetb
      REAL(KIND=dp), DIMENSION(6, 6)                     :: sab
      REAL(KIND=dp), DIMENSION(6, 6, 3)                  :: dab
      REAL(KIND=dp), DIMENSION(:), POINTER               :: atener
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: atstr
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atprop_type), POINTER                         :: atprop
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_gcp
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(virial_type), POINTER                         :: virial

      energy = 0._dp
      IF (.NOT. gcp_env%do_gcp) RETURN

      CALL timeset(routineN, handle)

      NULLIFY (atomic_kind_set, qs_kind_set, particle_set, sab_gcp)

      CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set, &
                      cell=cell, virial=virial, para_env=para_env, atprop=atprop)
      nkind = SIZE(atomic_kind_set)
      NULLIFY (particle_set)
      CALL get_qs_env(qs_env=qs_env, particle_set=particle_set)
      natom = SIZE(particle_set)

      verbose = gcp_env%verbose
      IF (verbose) THEN
         unit_nr = cp_logger_get_default_io_unit()
      ELSE
         unit_nr = -1
      END IF
      ! atomic energy and stress arrays
      atenergy = atprop%energy
      IF (atenergy) THEN
         CALL atprop_array_init(atprop%ategcp, natom)
         atener => atprop%ategcp
      END IF
      atstress = atprop%stress
      IF (atstress) THEN
         atstr => atprop%atstress
      END IF

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, *)
         WRITE (unit_nr, *) " Pair potential geometrical counterpoise (gCP) calculation"
         WRITE (unit_nr, *)
         WRITE (unit_nr, "(T15,A,T74,F7.4)") " Gloabal Parameters:     sigma = ", gcp_env%sigma, &
            "                         alpha = ", gcp_env%alpha, &
            "                         beta  = ", gcp_env%beta, &
            "                         eta   = ", gcp_env%eta
         WRITE (unit_nr, *)
         WRITE (unit_nr, "(T31,4(A5,10X))") " kind", "nvirt", "Emiss", " asto"
         DO ikind = 1, nkind
            WRITE (unit_nr, "(T31,i5,F15.1,F15.4,F15.4)") ikind, gcp_env%gcp_kind(ikind)%nbvirt, &
               gcp_env%gcp_kind(ikind)%eamiss, gcp_env%gcp_kind(ikind)%asto
         END DO
         WRITE (unit_nr, *)
      END IF

      IF (calculate_forces) THEN
         NULLIFY (force)
         CALL get_qs_env(qs_env=qs_env, force=force)
         ALLOCATE (atom_of_kind(natom), kind_of(natom))
         CALL get_atomic_kind_set(atomic_kind_set, atom_of_kind=atom_of_kind, kind_of=kind_of)
         use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)
         IF (use_virial) dvirial = virial%pv_virial
      END IF

      ! include all integrals in the list
      rcut = 1.e6_dp

      egcp = 0.0_dp
      IF (verbose) THEN
         ALLOCATE (egcpat(natom), ngcpat(natom))
         egcpat = 0.0_dp
         ngcpat = 0
      END IF

      nsto = 6
      DO ikind = 1, nkind
         CPASSERT(nsto == SIZE(gcp_env%gcp_kind(jkind)%al))
      END DO

      sab_gcp => gcp_env%sab_gcp
      CALL neighbor_list_iterator_create(nl_iterator, sab_gcp)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)

         CALL get_iterator_info(nl_iterator, mepos=mepos, ikind=ikind, jkind=jkind, iatom=iatom, jatom=jatom, r=rij)

         rcc = SQRT(rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3))
         IF (rcc > 1.e-6_dp) THEN
            fac = 1._dp
            IF (iatom == jatom) fac = 0.5_dp
            nvirta = gcp_env%gcp_kind(ikind)%nbvirt
            nvirtb = gcp_env%gcp_kind(jkind)%nbvirt
            eama = gcp_env%gcp_kind(ikind)%eamiss
            eamb = gcp_env%gcp_kind(jkind)%eamiss
            expab = EXP(-gcp_env%alpha*rcc**gcp_env%beta)
            zeta(1:nsto) = gcp_env%gcp_kind(ikind)%al(1:nsto)
            zetb(1:nsto) = gcp_env%gcp_kind(jkind)%al(1:nsto)
            cla(1:nsto) = gcp_env%gcp_kind(ikind)%cl(1:nsto)
            clb(1:nsto) = gcp_env%gcp_kind(jkind)%cl(1:nsto)
            IF (calculate_forces) THEN
               CALL overlap_ab(0, 0, nsto, rcut, zeta, 0, 0, nsto, rcut, zetb, rij, sab, dab)
               DO i = 1, 3
                  dsint(i) = SUM(cla*MATMUL(dab(:, :, i), clb))
               END DO
            ELSE
               CALL overlap_ab(0, 0, nsto, rcut, zeta, 0, 0, nsto, rcut, zetb, rij, sab)
            END IF
            sint = SUM(cla*MATMUL(sab, clb))
            IF (sint < 1.e-16_dp) CYCLE
            sqa = SQRT(sint*nvirta)
            sqb = SQRT(sint*nvirtb)
            IF (sqb > 1.e-12_dp) THEN
               fda = gcp_env%sigma*eama*expab/sqb
            ELSE
               fda = 0.0_dp
            END IF
            IF (sqa > 1.e-12_dp) THEN
               fdb = gcp_env%sigma*eamb*expab/sqa
            ELSE
               fdb = 0.0_dp
            END IF
            egcp = egcp + fac*(fda + fdb)
            IF (verbose) THEN
               egcpat(iatom) = egcpat(iatom) + fac*fda
               egcpat(jatom) = egcpat(jatom) + fac*fdb
               ngcpat(iatom) = ngcpat(iatom) + 1
               ngcpat(jatom) = ngcpat(jatom) + 1
            END IF
            IF (calculate_forces) THEN
               fdij = -fac*(fda + fdb)*(gcp_env%alpha*gcp_env%beta*rcc**(gcp_env%beta - 1.0_dp)*rij(1:3)/rcc)
               IF (sqa > 1.e-12_dp) THEN
                  fdij = fdij + 0.25_dp*fac*fdb/(sqa*sqa)*dsint(1:3)
               END IF
               IF (sqb > 1.e-12_dp) THEN
                  fdij = fdij + 0.25_dp*fac*fda/(sqb*sqb)*dsint(1:3)
               END IF
               atom_a = atom_of_kind(iatom)
               atom_b = atom_of_kind(jatom)
               force(ikind)%gcp(:, atom_a) = force(ikind)%gcp(:, atom_a) - fdij(:)
               force(jkind)%gcp(:, atom_b) = force(jkind)%gcp(:, atom_b) + fdij(:)
               IF (use_virial) THEN
                  CALL virial_pair_force(virial%pv_virial, -1._dp, fdij, rij)
               END IF
               IF (atstress) THEN
                  CALL virial_pair_force(atstr(:, :, iatom), -0.5_dp, fdij, rij)
                  CALL virial_pair_force(atstr(:, :, jatom), -0.5_dp, fdij, rij)
               END IF
            END IF
            IF (atenergy) THEN
               atener(iatom) = atener(iatom) + fda*fac
               atener(jatom) = atener(jatom) + fdb*fac
            END IF
         END IF
      END DO

      CALL neighbor_list_iterator_release(nl_iterator)

      ! set gCP energy
      CALL mp_sum(egcp, para_env%group)
      energy = egcp
      IF (verbose) THEN
         CALL mp_sum(egcpat, para_env%group)
         CALL mp_sum(ngcpat, para_env%group)
      END IF

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, "(T15,A,T61,F20.10)") " Total gCP energy [au]     :", egcp
         WRITE (unit_nr, "(T15,A,T61,F20.10)") " Total gCP energy [kcal]   :", egcp*kcalmol
         WRITE (unit_nr, *)
         WRITE (unit_nr, "(T19,A)") " gCP atomic energy contributions"
         WRITE (unit_nr, "(T19,A,T60,A20)") " #             sites", "      BSSE [kcal/mol]"
         DO i = 1, natom
            WRITE (unit_nr, "(12X,I8,10X,I8,T61,F20.10)") i, ngcpat(i), egcpat(i)*kcalmol
         END DO
      END IF
      IF (calculate_forces) THEN
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, *) " gCP Forces         "
            WRITE (unit_nr, *) " Atom   Kind                            Forces    "
         END IF
         gnorm = 0._dp
         DO iatom = 1, natom
            ikind = kind_of(iatom)
            atom_a = atom_of_kind(iatom)
            fdij(1:3) = force(ikind)%gcp(:, atom_a)
            CALL mp_sum(fdij, para_env%group)
            gnorm = gnorm + SUM(ABS(fdij))
            IF (unit_nr > 0) WRITE (unit_nr, "(i5,i7,3F20.14)") iatom, ikind, fdij
         END DO
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, *)
            WRITE (unit_nr, *) " |G| = ", gnorm
            WRITE (unit_nr, *)
         END IF
         IF (use_virial) THEN
            dvirial = virial%pv_virial - dvirial
            CALL mp_sum(dvirial, para_env%group)
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, *) " Stress Tensor (gCP)"
               WRITE (unit_nr, "(3G20.12)") dvirial
               WRITE (unit_nr, *) "  Tr(P)/3 :  ", (dvirial(1, 1) + dvirial(2, 2) + dvirial(3, 3))/3._dp
               WRITE (unit_nr, *)
            END IF
         END IF
      END IF
      IF (verbose) THEN
         DEALLOCATE (egcpat, ngcpat)
      END IF

      CALL timestop(handle)

   END SUBROUTINE calculate_gcp_pairpot

! **************************************************************************************************

END MODULE qs_gcp_method
