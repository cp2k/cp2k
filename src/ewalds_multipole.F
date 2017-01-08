!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Treats the electrostatic for multipoles (up to quadrupoles)
!> \author Teodoro Laino [tlaino] - 12.2007 - University of Zurich
!> inclusion of optional electric field damping for the polarizable atoms
!> Rodolphe Vuilleumier and Mathieu Salanne - 12.2009
! **************************************************************************************************
MODULE ewalds_multipole
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE bibliography,                    ONLY: Aguado2003,&
                                              Laino2008,&
                                              cite_reference
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE damping_dipole_types,            ONLY: damping_type,&
                                              no_damping,&
                                              tang_toennies
   USE dg_rho0_types,                   ONLY: dg_rho0_type
   USE dg_types,                        ONLY: dg_get,&
                                              dg_type
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE erf_fn,                          ONLY: erf,&
                                              erfc
   USE ewald_environment_types,         ONLY: ewald_env_get,&
                                              ewald_environment_type
   USE ewald_pw_types,                  ONLY: ewald_pw_get,&
                                              ewald_pw_type
   USE fist_neighbor_list_control,      ONLY: list_control
   USE fist_neighbor_list_types,        ONLY: fist_neighbor_type,&
                                              neighbor_kind_pairs_type
   USE fist_nonbond_env_types,          ONLY: fist_nonbond_env_get,&
                                              fist_nonbond_env_type,&
                                              pos_type
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: fourpi,&
                                              oorootpi,&
                                              pi,&
                                              sqrthalf
   USE mathlib,                         ONLY: matvec_3x3
   USE message_passing,                 ONLY: mp_sum
   USE parallel_rng_types,              ONLY: UNIFORM,&
                                              create_rng_stream,&
                                              delete_rng_stream,&
                                              random_numbers,&
                                              rng_stream_type
   USE particle_types,                  ONLY: particle_type
   USE pw_grid_types,                   ONLY: pw_grid_type
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE structure_factor_types,          ONLY: structure_factor_type
   USE structure_factors,               ONLY: structure_factor_allocate,&
                                              structure_factor_deallocate,&
                                              structure_factor_evaluate
#include "./base/base_uses.f90"

#:include "ewalds_multipole_sr.fypp"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .FALSE.
   LOGICAL, PRIVATE, PARAMETER :: debug_r_space = .FALSE.
   LOGICAL, PRIVATE, PARAMETER :: debug_g_space = .FALSE.
   LOGICAL, PRIVATE, PARAMETER :: debug_e_field = .FALSE.
   LOGICAL, PRIVATE, PARAMETER :: debug_e_field_en = .FALSE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ewalds_multipole'

   PUBLIC :: ewald_multipole_evaluate

CONTAINS

! **************************************************************************************************
!> \brief  Computes the potential and the force for a lattice sum of multipoles (up to quadrupole)
!> \param ewald_env ...
!> \param ewald_pw ...
!> \param nonbond_env ...
!> \param cell ...
!> \param particle_set ...
!> \param local_particles ...
!> \param energy_local ...
!> \param energy_glob ...
!> \param e_neut ...
!> \param e_self ...
!> \param task ...
!> \param do_correction_bonded ...
!> \param do_forces ...
!> \param do_stress ...
!> \param do_efield ...
!> \param radii ...
!> \param charges ...
!> \param dipoles ...
!> \param quadrupoles ...
!> \param forces_local ...
!> \param forces_glob ...
!> \param pv_local ...
!> \param pv_glob ...
!> \param efield0 ...
!> \param efield1 ...
!> \param efield2 ...
!> \param iw ...
!> \param do_debug ...
!> \param atomic_kind_set ...
!> \param mm_section ...
!> \par    Note
!>         atomic_kind_set and mm_section are between the arguments only
!>         for debug purpose (therefore optional) and can be avoided when this
!>         function is called in other part of the program
!> \par    Note
!>         When a gaussian multipole is used instead of point multipole, i.e.
!>         when radii(i)>0, the electrostatic fields (efield0, efield1, efield2)
!>         become derivatives of the electrostatic potential energy towards
!>         these gaussian multipoles.
!> \author Teodoro Laino [tlaino] - 12.2007 - University of Zurich
! **************************************************************************************************
   RECURSIVE SUBROUTINE ewald_multipole_evaluate(ewald_env, ewald_pw, nonbond_env, &
                                                 cell, particle_set, local_particles, energy_local, energy_glob, e_neut, e_self, &
                                                 task, do_correction_bonded, do_forces, do_stress, &
                                                 do_efield, radii, charges, dipoles, &
                                                 quadrupoles, forces_local, forces_glob, pv_local, pv_glob, efield0, efield1, &
                                                 efield2, iw, do_debug, atomic_kind_set, mm_section)
      TYPE(ewald_environment_type), POINTER              :: ewald_env
      TYPE(ewald_pw_type), POINTER                       :: ewald_pw
      TYPE(fist_nonbond_env_type), POINTER               :: nonbond_env
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), POINTER                       :: particle_set(:)
      TYPE(distribution_1d_type), POINTER                :: local_particles
      REAL(KIND=dp), INTENT(INOUT)                       :: energy_local, energy_glob
      REAL(KIND=dp), INTENT(OUT)                         :: e_neut, e_self
      LOGICAL, DIMENSION(3), INTENT(IN)                  :: task
      LOGICAL, INTENT(IN)                                :: do_correction_bonded, do_forces, &
                                                            do_stress, do_efield
      REAL(KIND=dp), DIMENSION(:), OPTIONAL, POINTER     :: radii, charges
      REAL(KIND=dp), DIMENSION(:, :), OPTIONAL, POINTER  :: dipoles
      REAL(KIND=dp), DIMENSION(:, :, :), OPTIONAL, &
         POINTER                                         :: quadrupoles
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT), &
         OPTIONAL                                        :: forces_local, forces_glob, pv_local, &
                                                            pv_glob
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT), OPTIONAL :: efield0
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT), &
         OPTIONAL                                        :: efield1, efield2
      INTEGER, INTENT(IN)                                :: iw
      LOGICAL, INTENT(IN)                                :: do_debug
      TYPE(atomic_kind_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: atomic_kind_set
      TYPE(section_vals_type), OPTIONAL, POINTER         :: mm_section

      CHARACTER(len=*), PARAMETER :: routineN = 'ewald_multipole_evaluate', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: group, handle, i, j, size1, size2
      LOGICAL                                            :: check_debug, check_efield, check_forces, &
                                                            do_task(3)
      LOGICAL, DIMENSION(3, 3)                           :: my_task
      REAL(KIND=dp)                                      :: e_bonded, e_bonded_t, e_rspace, &
                                                            e_rspace_t, energy_glob_t
      REAL(KIND=dp), DIMENSION(:), POINTER               :: efield0_lr, efield0_sr
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: efield1_lr, efield1_sr, efield2_lr, &
                                                            efield2_sr

      CALL cite_reference(Aguado2003)
      CALL cite_reference(Laino2008)
      CALL timeset(routineN, handle)
      CPASSERT(ASSOCIATED(nonbond_env))
      check_debug = (debug_this_module .OR. debug_r_space .OR. debug_g_space .OR. debug_e_field .OR. debug_e_field_en) &
                    .EQV. debug_this_module
      CPASSERT(check_debug)
      check_forces = do_forces .EQV. (PRESENT(forces_local) .AND. PRESENT(forces_glob))
      CPASSERT(check_forces)
      check_efield = do_efield .EQV. (PRESENT(efield0) .OR. PRESENT(efield1) .OR. PRESENT(efield2))
      CPASSERT(check_efield)
      ! Debugging this module
      IF (debug_this_module .AND. do_debug) THEN
         ! Debug specifically real space part
         IF (debug_r_space) THEN
            CALL debug_ewald_multipoles(ewald_env, ewald_pw, nonbond_env, cell, &
                                        particle_set, local_particles, iw, debug_r_space)
            CPABORT("Debug Multipole Requested:  Real Part!")
         END IF
         ! Debug electric fields and gradients as pure derivatives
         IF (debug_e_field) THEN
            CPASSERT(PRESENT(atomic_kind_set))
            CPASSERT(PRESENT(mm_section))
            CALL debug_ewald_multipoles_fields(ewald_env, ewald_pw, nonbond_env, &
                                               cell, particle_set, local_particles, radii, charges, dipoles, &
                                               quadrupoles, task, iw, atomic_kind_set, mm_section)
            CPABORT("Debug Multipole Requested:  POT+EFIELDS+GRAD!")
         END IF
         ! Debug the potential, electric fields and electric fields gradient in oder
         ! to retrieve the correct energy
         IF (debug_e_field_en) THEN
            CALL debug_ewald_multipoles_fields2(ewald_env, ewald_pw, nonbond_env, &
                                                cell, particle_set, local_particles, radii, charges, dipoles, &
                                                quadrupoles, task, iw)
            CPABORT("Debug Multipole Requested:  POT+EFIELDS+GRAD to give the correct energy!!")
         END IF
      END IF

      ! Setup the tasks (needed to skip useless parts in the real-space term)
      do_task = task
      DO i = 1, 3
         IF (do_task(i)) THEN
            SELECT CASE (i)
            CASE (1)
               do_task(1) = ANY(charges /= 0.0_dp)
            CASE (2)
               do_task(2) = ANY(dipoles /= 0.0_dp)
            CASE (3)
               do_task(3) = ANY(quadrupoles /= 0.0_dp)
            END SELECT
         END IF
      END DO
      DO i = 1, 3
         DO j = i, 3
            my_task(j, i) = do_task(i) .AND. do_task(j)
            my_task(i, j) = my_task(j, i)
         END DO
      END DO

      ! Allocate arrays for the evaluation of the potential, fields and electrostatic field gradients
      NULLIFY (efield0_sr, efield0_lr, efield1_sr, efield1_lr, efield2_sr, efield2_lr)
      IF (do_efield) THEN
         IF (PRESENT(efield0)) THEN
            size1 = SIZE(efield0)
            ALLOCATE (efield0_sr(size1))
            ALLOCATE (efield0_lr(size1))
            efield0_sr = 0.0_dp
            efield0_lr = 0.0_dp
         END IF
         IF (PRESENT(efield1)) THEN
            size1 = SIZE(efield1, 1)
            size2 = SIZE(efield1, 2)
            ALLOCATE (efield1_sr(size1, size2))
            ALLOCATE (efield1_lr(size1, size2))
            efield1_sr = 0.0_dp
            efield1_lr = 0.0_dp
         END IF
         IF (PRESENT(efield2)) THEN
            size1 = SIZE(efield2, 1)
            size2 = SIZE(efield2, 2)
            ALLOCATE (efield2_sr(size1, size2))
            ALLOCATE (efield2_lr(size1, size2))
            efield2_sr = 0.0_dp
            efield2_lr = 0.0_dp
         END IF
      END IF

      e_rspace = 0.0_dp
      e_bonded = 0.0_dp
      IF ((.NOT. debug_g_space) .AND. (nonbond_env%do_nonbonded)) THEN
         ! Compute the Real Space (Short-Range) part of the Ewald sum.
         ! This contribution is only added when the nonbonded flag in the input
         ! is set, because these contributions depend. the neighborlists.
         CALL ewald_multipole_SR(nonbond_env, ewald_env, atomic_kind_set, &
                                 particle_set, cell, e_rspace, my_task, &
                                 do_forces, do_efield, do_stress, radii, charges, dipoles, quadrupoles, &
                                 forces_glob, pv_glob, efield0_sr, efield1_sr, efield2_sr)
         energy_glob = energy_glob+e_rspace

         IF (do_correction_bonded) THEN
            ! The corrections for bonded interactions are stored in the Real Space
            ! (Short-Range) part of the fields array.
            CALL ewald_multipole_bonded(nonbond_env, particle_set, ewald_env, &
                                        cell, e_bonded, my_task, do_forces, do_efield, do_stress, &
                                        charges, dipoles, quadrupoles, forces_glob, pv_glob, &
                                        efield0_sr, efield1_sr, efield2_sr)
            energy_glob = energy_glob+e_bonded
         END IF
      END IF

      e_neut = 0.0_dp
      e_self = 0.0_dp
      energy_local = 0.0_dp
      IF (.NOT. debug_r_space) THEN
         ! Compute the Reciprocal Space (Long-Range) part of the Ewald sum
         CALL ewald_multipole_LR(ewald_env, ewald_pw, cell, particle_set, &
                                 local_particles, energy_local, my_task, do_forces, do_efield, do_stress, &
                                 charges, dipoles, quadrupoles, forces_local, pv_local, efield0_lr, efield1_lr, &
                                 efield2_lr)

         ! Self-Interactions corrections
         CALL ewald_multipole_self(ewald_env, cell, local_particles, e_self, &
                                   e_neut, my_task, do_efield, radii, charges, dipoles, quadrupoles, &
                                   efield0_lr, efield1_lr, efield2_lr)
      END IF

      ! Sumup energy contributions for possible IO
      CALL ewald_env_get(ewald_env, group=group)
      energy_glob_t = energy_glob
      e_rspace_t = e_rspace
      e_bonded_t = e_bonded
      CALL mp_sum(energy_glob_t, group)
      CALL mp_sum(e_rspace_t, group)
      CALL mp_sum(e_bonded_t, group)
      ! Print some info about energetics
      CALL ewald_multipole_print(iw, energy_local, e_rspace_t, e_bonded_t, e_self, e_neut)

      ! Gather the components of the potential, fields and electrostatic field gradients
      IF (do_efield) THEN
         IF (PRESENT(efield0)) THEN
            efield0 = efield0_sr+efield0_lr
            CALL mp_sum(efield0, group)
            DEALLOCATE (efield0_sr)
            DEALLOCATE (efield0_lr)
         END IF
         IF (PRESENT(efield1)) THEN
            efield1 = efield1_sr+efield1_lr
            CALL mp_sum(efield1, group)
            DEALLOCATE (efield1_sr)
            DEALLOCATE (efield1_lr)
         END IF
         IF (PRESENT(efield2)) THEN
            efield2 = efield2_sr+efield2_lr
            CALL mp_sum(efield2, group)
            DEALLOCATE (efield2_sr)
            DEALLOCATE (efield2_lr)
         END IF
      END IF
      CALL timestop(handle)
   END SUBROUTINE ewald_multipole_evaluate

! **************************************************************************************************
!> \brief computes the potential and the force for a lattice sum of multipoles
!>      up to quadrupole - Short Range (Real Space) Term
!> \param nonbond_env ...
!> \param ewald_env ...
!> \param atomic_kind_set ...
!> \param particle_set ...
!> \param cell ...
!> \param energy ...
!> \param task ...
!> \param do_forces ...
!> \param do_efield ...
!> \param do_stress ...
!> \param radii ...
!> \param charges ...
!> \param dipoles ...
!> \param quadrupoles ...
!> \param forces ...
!> \param pv ...
!> \param efield0 ...
!> \param efield1 ...
!> \param efield2 ...
!> \author Teodoro Laino [tlaino] - 12.2007 - University of Zurich
! **************************************************************************************************
   SUBROUTINE ewald_multipole_SR(nonbond_env, ewald_env, atomic_kind_set, &
                                 particle_set, cell, energy, task, &
                                 do_forces, do_efield, do_stress, radii, charges, dipoles, quadrupoles, &
                                 forces, pv, efield0, efield1, efield2)
      TYPE(fist_nonbond_env_type), POINTER               :: nonbond_env
      TYPE(ewald_environment_type), POINTER              :: ewald_env
      TYPE(atomic_kind_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: atomic_kind_set
      TYPE(particle_type), POINTER                       :: particle_set(:)
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), INTENT(INOUT)                       :: energy
      LOGICAL, DIMENSION(3, 3), INTENT(IN)               :: task
      LOGICAL, INTENT(IN)                                :: do_forces, do_efield, do_stress
      REAL(KIND=dp), DIMENSION(:), OPTIONAL, POINTER     :: radii, charges
      REAL(KIND=dp), DIMENSION(:, :), OPTIONAL, POINTER  :: dipoles
      REAL(KIND=dp), DIMENSION(:, :, :), OPTIONAL, &
         POINTER                                         :: quadrupoles
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT), &
         OPTIONAL                                        :: forces, pv
      REAL(KIND=dp), DIMENSION(:), POINTER               :: efield0
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: efield1, efield2

      CHARACTER(len=*), PARAMETER :: routineN = 'ewald_multipole_SR', &
         routineP = moduleN//':'//routineN

      INTEGER :: a, atom_a, atom_b, b, c, d, e, handle, i, iend, igrp, ikind, ilist, ipair, &
         istart, itype_ij, itype_ji, jkind, k, kind_a, kind_b, kk, nkdamp_ij, nkdamp_ji, nkinds, &
         npairs
      INTEGER, DIMENSION(:, :), POINTER                  :: list
      LOGICAL                                            :: do_efield0, do_efield1, do_efield2, &
                                                            force_eval
      REAL(KIND=dp) :: alpha, beta, ch_i, ch_j, dampa_ij, dampa_ji, dampaexpi, dampaexpj, &
         dampfac_ij, dampfac_ji, dampfuncdiffi, dampfuncdiffj, dampfunci, dampfuncj, dampsumfi, &
         dampsumfj, ef0_i, ef0_j, eloc, fac, fac_ij, factorial, ir, irab2, ptens11, ptens12, &
         ptens13, ptens21, ptens22, ptens23, ptens31, ptens32, ptens33, r, rab2, rab2_max, radius, &
         rcut, tij, tmp, tmp1, tmp11, tmp12, tmp13, tmp2, tmp21, tmp22, tmp23, tmp31, tmp32, &
         tmp33, tmp_ij, tmp_ji, xf
      REAL(KIND=dp), DIMENSION(0:5)                      :: f
      REAL(KIND=dp), DIMENSION(3)                        :: cell_v, cvi, damptij_a, damptji_a, dp_i, &
                                                            dp_j, ef1_i, ef1_j, fr, rab, tij_a
      REAL(KIND=dp), DIMENSION(3, 3)                     :: damptij_ab, damptji_ab, ef2_i, ef2_j, &
                                                            qp_i, qp_j, tij_ab
      REAL(KIND=dp), DIMENSION(3, 3, 3)                  :: tij_abc
      REAL(KIND=dp), DIMENSION(3, 3, 3, 3)               :: tij_abcd
      REAL(KIND=dp), DIMENSION(3, 3, 3, 3, 3)            :: tij_abcde
      TYPE(damping_type)                                 :: damping_ij, damping_ji
      TYPE(fist_neighbor_type), POINTER                  :: nonbonded
      TYPE(neighbor_kind_pairs_type), POINTER            :: neighbor_kind_pair
      TYPE(pos_type), DIMENSION(:), POINTER              :: r_last_update, r_last_update_pbc

      CALL timeset(routineN, handle)
      NULLIFY (nonbonded, r_last_update, r_last_update_pbc)
      do_efield0 = do_efield .AND. ASSOCIATED(efield0)
      do_efield1 = do_efield .AND. ASSOCIATED(efield1)
      do_efield2 = do_efield .AND. ASSOCIATED(efield2)
      IF (do_stress) THEN
         ptens11 = 0.0_dp; ptens12 = 0.0_dp; ptens13 = 0.0_dp
         ptens21 = 0.0_dp; ptens22 = 0.0_dp; ptens23 = 0.0_dp
         ptens31 = 0.0_dp; ptens32 = 0.0_dp; ptens33 = 0.0_dp
      END IF
      ! Get nonbond_env info
      CALL fist_nonbond_env_get(nonbond_env, nonbonded=nonbonded, natom_types=nkinds, &
                                r_last_update=r_last_update, r_last_update_pbc=r_last_update_pbc)
      CALL ewald_env_get(ewald_env, alpha=alpha, rcut=rcut)
      rab2_max = rcut**2
      IF (debug_r_space) THEN
         rab2_max = HUGE(0.0_dp)
      END IF
      ! Starting the force loop
      Lists: DO ilist = 1, nonbonded%nlists
         neighbor_kind_pair => nonbonded%neighbor_kind_pairs(ilist)
         npairs = neighbor_kind_pair%npairs
         IF (npairs == 0) CYCLE
         list => neighbor_kind_pair%list
         cvi = neighbor_kind_pair%cell_vector
         CALL matvec_3x3(cell_v, cell%hmat, cvi)
         Kind_Group_Loop: DO igrp = 1, neighbor_kind_pair%ngrp_kind
            istart = neighbor_kind_pair%grp_kind_start(igrp)
            iend = neighbor_kind_pair%grp_kind_end(igrp)
            ikind = neighbor_kind_pair%ij_kind(1, igrp)
            jkind = neighbor_kind_pair%ij_kind(2, igrp)

            itype_ij = no_damping
            nkdamp_ij = 1
            dampa_ij = 1.0_dp
            dampfac_ij = 0.0_dp

            itype_ji = no_damping
            nkdamp_ji = 1
            dampa_ji = 1.0_dp
            dampfac_ji = 0.0_dp
            IF (PRESENT(atomic_kind_set)) THEN
               IF (ASSOCIATED(atomic_kind_set(jkind)%damping)) THEN
                  damping_ij = atomic_kind_set(jkind)%damping%damp(ikind)
                  itype_ij = damping_ij%itype
                  nkdamp_ij = damping_ij%order
                  dampa_ij = damping_ij%bij
                  dampfac_ij = damping_ij%cij
               END IF

               IF (ASSOCIATED(atomic_kind_set(ikind)%damping)) THEN
                  damping_ji = atomic_kind_set(ikind)%damping%damp(jkind)
                  itype_ji = damping_ji%itype
                  nkdamp_ji = damping_ji%order
                  dampa_ji = damping_ji%bij
                  dampfac_ji = damping_ji%cij
               END IF
            END IF

            Pairs: DO ipair = istart, iend
               IF (ipair <= neighbor_kind_pair%nscale) THEN
                  ! scale the electrostatic interaction if needed
                  ! (most often scaled to zero)
                  fac_ij = neighbor_kind_pair%ei_scale(ipair)
                  IF (fac_ij <= 0) CYCLE
               ELSE
                  fac_ij = 1.0_dp
               END IF
               atom_a = list(1, ipair)
               atom_b = list(2, ipair)
               kind_a = particle_set(atom_a)%atomic_kind%kind_number
               kind_b = particle_set(atom_b)%atomic_kind%kind_number
               IF (atom_a == atom_b) fac_ij = 0.5_dp
               rab = r_last_update_pbc(atom_b)%r-r_last_update_pbc(atom_a)%r
               rab = rab+cell_v
               rab2 = rab(1)**2+rab(2)**2+rab(3)**2
               IF (rab2 <= rab2_max) THEN
                  IF (PRESENT(radii)) THEN
                     radius = SQRT(radii(atom_a)*radii(atom_a)+radii(atom_b)*radii(atom_b))
                  ELSE
                     radius = 0.0_dp
                  END IF
                  IF (radius > 0.0_dp) THEN
                     beta = sqrthalf/radius
$: ewalds_multipole_sr_macro(mode="SCREENED_COULOMB_GAUSS", damping=False, store_energy=True, store_forces=True)
                  ELSE
$: ewalds_multipole_sr_macro(mode="SCREENED_COULOMB_ERFC", damping=True, store_energy=True, store_forces=True )
                  END IF
               END IF
            END DO Pairs
         END DO Kind_Group_Loop
      END DO Lists
      IF (do_stress) THEN
         pv(1, 1) = pv(1, 1)+ptens11
         pv(1, 2) = pv(1, 2)+(ptens12+ptens21)*0.5_dp
         pv(1, 3) = pv(1, 3)+(ptens13+ptens31)*0.5_dp
         pv(2, 1) = pv(1, 2)
         pv(2, 2) = pv(2, 2)+ptens22
         pv(2, 3) = pv(2, 3)+(ptens23+ptens32)*0.5_dp
         pv(3, 1) = pv(1, 3)
         pv(3, 2) = pv(2, 3)
         pv(3, 3) = pv(3, 3)+ptens33
      END IF

      CALL timestop(handle)
   END SUBROUTINE ewald_multipole_SR

! **************************************************************************************************
!> \brief computes the bonded correction for the potential and the force for a
!>        lattice sum of multipoles up to quadrupole
!> \param nonbond_env ...
!> \param particle_set ...
!> \param ewald_env ...
!> \param cell ...
!> \param energy ...
!> \param task ...
!> \param do_forces ...
!> \param do_efield ...
!> \param do_stress ...
!> \param charges ...
!> \param dipoles ...
!> \param quadrupoles ...
!> \param forces ...
!> \param pv ...
!> \param efield0 ...
!> \param efield1 ...
!> \param efield2 ...
!> \author Teodoro Laino [tlaino] - 05.2009
! **************************************************************************************************
   SUBROUTINE ewald_multipole_bonded(nonbond_env, particle_set, ewald_env, &
                                     cell, energy, task, do_forces, do_efield, do_stress, charges, &
                                     dipoles, quadrupoles, forces, pv, efield0, efield1, efield2)

      TYPE(fist_nonbond_env_type), POINTER               :: nonbond_env
      TYPE(particle_type), POINTER                       :: particle_set(:)
      TYPE(ewald_environment_type), POINTER              :: ewald_env
      TYPE(cell_type), POINTER                           :: cell
      REAL(KIND=dp), INTENT(INOUT)                       :: energy
      LOGICAL, DIMENSION(3, 3), INTENT(IN)               :: task
      LOGICAL, INTENT(IN)                                :: do_forces, do_efield, do_stress
      REAL(KIND=dp), DIMENSION(:), OPTIONAL, POINTER     :: charges
      REAL(KIND=dp), DIMENSION(:, :), OPTIONAL, POINTER  :: dipoles
      REAL(KIND=dp), DIMENSION(:, :, :), OPTIONAL, &
         POINTER                                         :: quadrupoles
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT), &
         OPTIONAL                                        :: forces, pv
      REAL(KIND=dp), DIMENSION(:), POINTER               :: efield0
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: efield1, efield2

      CHARACTER(len=*), PARAMETER :: routineN = 'ewald_multipole_bonded', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: a, atom_a, atom_b, b, c, d, e, handle, &
                                                            i, iend, igrp, ilist, ipair, istart, &
                                                            k, nscale
      INTEGER, DIMENSION(:, :), POINTER                  :: list
      LOGICAL                                            :: do_efield0, do_efield1, do_efield2, &
                                                            force_eval
      REAL(KIND=dp) :: alpha, ch_i, ch_j, ef0_i, ef0_j, eloc, fac, fac_ij, ir, irab2, ptens11, &
         ptens12, ptens13, ptens21, ptens22, ptens23, ptens31, ptens32, ptens33, r, rab2, tij, &
         tmp, tmp1, tmp11, tmp12, tmp13, tmp2, tmp21, tmp22, tmp23, tmp31, tmp32, tmp33, tmp_ij, &
         tmp_ji
      REAL(KIND=dp), DIMENSION(0:5)                      :: f
      REAL(KIND=dp), DIMENSION(3)                        :: dp_i, dp_j, ef1_i, ef1_j, fr, rab, tij_a
      REAL(KIND=dp), DIMENSION(3, 3)                     :: ef2_i, ef2_j, qp_i, qp_j, tij_ab
      REAL(KIND=dp), DIMENSION(3, 3, 3)                  :: tij_abc
      REAL(KIND=dp), DIMENSION(3, 3, 3, 3)               :: tij_abcd
      REAL(KIND=dp), DIMENSION(3, 3, 3, 3, 3)            :: tij_abcde
      TYPE(fist_neighbor_type), POINTER                  :: nonbonded
      TYPE(neighbor_kind_pairs_type), POINTER            :: neighbor_kind_pair

      CALL timeset(routineN, handle)
      do_efield0 = do_efield .AND. ASSOCIATED(efield0)
      do_efield1 = do_efield .AND. ASSOCIATED(efield1)
      do_efield2 = do_efield .AND. ASSOCIATED(efield2)
      IF (do_stress) THEN
         ptens11 = 0.0_dp; ptens12 = 0.0_dp; ptens13 = 0.0_dp
         ptens21 = 0.0_dp; ptens22 = 0.0_dp; ptens23 = 0.0_dp
         ptens31 = 0.0_dp; ptens32 = 0.0_dp; ptens33 = 0.0_dp
      END IF
      CALL ewald_env_get(ewald_env, alpha=alpha)
      CALL fist_nonbond_env_get(nonbond_env, nonbonded=nonbonded)

      ! Starting the force loop
      Lists: DO ilist = 1, nonbonded%nlists
         neighbor_kind_pair => nonbonded%neighbor_kind_pairs(ilist)
         nscale = neighbor_kind_pair%nscale
         IF (nscale == 0) CYCLE
         list => neighbor_kind_pair%list
         Kind_Group_Loop: DO igrp = 1, neighbor_kind_pair%ngrp_kind
            istart = neighbor_kind_pair%grp_kind_start(igrp)
            IF (istart > nscale) CYCLE
            iend = MIN(neighbor_kind_pair%grp_kind_end(igrp), nscale)
            Pairs: DO ipair = istart, iend
               ! only use pairs that are (partially) excluded for electrostatics
               fac_ij = -1.0_dp+neighbor_kind_pair%ei_scale(ipair)
               IF (fac_ij >= 0) CYCLE

               atom_a = list(1, ipair)
               atom_b = list(2, ipair)

               rab = particle_set(atom_b)%r-particle_set(atom_a)%r
               rab = pbc(rab, cell)
               rab2 = rab(1)**2+rab(2)**2+rab(3)**2
$: ewalds_multipole_sr_macro(mode="SCREENED_COULOMB_ERF", store_energy=True, store_forces=True)
            END DO Pairs
         END DO Kind_Group_Loop
      END DO Lists
      IF (do_stress) THEN
         pv(1, 1) = pv(1, 1)+ptens11
         pv(1, 2) = pv(1, 2)+(ptens12+ptens21)*0.5_dp
         pv(1, 3) = pv(1, 3)+(ptens13+ptens31)*0.5_dp
         pv(2, 1) = pv(1, 2)
         pv(2, 2) = pv(2, 2)+ptens22
         pv(2, 3) = pv(2, 3)+(ptens23+ptens32)*0.5_dp
         pv(3, 1) = pv(1, 3)
         pv(3, 2) = pv(2, 3)
         pv(3, 3) = pv(3, 3)+ptens33
      END IF

      CALL timestop(handle)
   END SUBROUTINE ewald_multipole_bonded

! **************************************************************************************************
!> \brief computes the potential and the force for a lattice sum of multipoles
!>      up to quadrupole - Long Range (Reciprocal Space) Term
!> \param ewald_env ...
!> \param ewald_pw ...
!> \param cell ...
!> \param particle_set ...
!> \param local_particles ...
!> \param energy ...
!> \param task ...
!> \param do_forces ...
!> \param do_efield ...
!> \param do_stress ...
!> \param charges ...
!> \param dipoles ...
!> \param quadrupoles ...
!> \param forces ...
!> \param pv ...
!> \param efield0 ...
!> \param efield1 ...
!> \param efield2 ...
!> \author Teodoro Laino [tlaino] - 12.2007 - University of Zurich
! **************************************************************************************************
   SUBROUTINE ewald_multipole_LR(ewald_env, ewald_pw, cell, particle_set, &
                                 local_particles, energy, task, do_forces, do_efield, do_stress, &
                                 charges, dipoles, quadrupoles, forces, pv, efield0, efield1, efield2)
      TYPE(ewald_environment_type), POINTER              :: ewald_env
      TYPE(ewald_pw_type), POINTER                       :: ewald_pw
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), POINTER                       :: particle_set(:)
      TYPE(distribution_1d_type), POINTER                :: local_particles
      REAL(KIND=dp), INTENT(INOUT)                       :: energy
      LOGICAL, DIMENSION(3, 3), INTENT(IN)               :: task
      LOGICAL, INTENT(IN)                                :: do_forces, do_efield, do_stress
      REAL(KIND=dp), DIMENSION(:), OPTIONAL, POINTER     :: charges
      REAL(KIND=dp), DIMENSION(:, :), OPTIONAL, POINTER  :: dipoles
      REAL(KIND=dp), DIMENSION(:, :, :), OPTIONAL, &
         POINTER                                         :: quadrupoles
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT), &
         OPTIONAL                                        :: forces, pv
      REAL(KIND=dp), DIMENSION(:), POINTER               :: efield0
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: efield1, efield2

      CHARACTER(len=*), PARAMETER :: routineN = 'ewald_multipole_LR', &
         routineP = moduleN//':'//routineN

      COMPLEX(KIND=dp)                                   :: atm_factor, atm_factor_st(3), cnjg_fac, &
                                                            fac, summe_tmp
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: summe_ef
      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:, :)     :: summe_st
      INTEGER :: gpt, group, handle, iparticle, iparticle_kind, iparticle_local, lp, mp, nnodes, &
         node, np, nparticle_kind, nparticle_local
      INTEGER, DIMENSION(:, :), POINTER                  :: bds
      LOGICAL                                            :: do_efield0, do_efield1, do_efield2
      REAL(KIND=dp)                                      :: alpha, denom, dipole_t(3), f0, factor, &
                                                            four_alpha_sq, gauss, pref, q_t, tmp, &
                                                            trq_t
      REAL(KIND=dp), DIMENSION(3)                        :: tmp_v, vec
      REAL(KIND=dp), DIMENSION(3, 3)                     :: pv_tmp
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: rho0
      TYPE(dg_rho0_type), POINTER                        :: dg_rho0
      TYPE(dg_type), POINTER                             :: dg
      TYPE(pw_grid_type), POINTER                        :: pw_grid
      TYPE(pw_pool_type), POINTER                        :: pw_pool
      TYPE(structure_factor_type)                        :: exp_igr

      CALL timeset(routineN, handle)
      do_efield0 = do_efield .AND. ASSOCIATED(efield0)
      do_efield1 = do_efield .AND. ASSOCIATED(efield1)
      do_efield2 = do_efield .AND. ASSOCIATED(efield2)

      ! Gathering data from the ewald environment
      CALL ewald_env_get(ewald_env, alpha=alpha, group=group)
      CALL ewald_pw_get(ewald_pw, pw_big_pool=pw_pool, dg=dg)
      CALL dg_get(dg, dg_rho0=dg_rho0)
      rho0 => dg_rho0%density%pw%cr3d
      pw_grid => pw_pool%pw_grid
      bds => pw_grid%bounds

      ! Allocation of working arrays
      nparticle_kind = SIZE(local_particles%n_el)
      nnodes = 0
      DO iparticle_kind = 1, nparticle_kind
         nnodes = nnodes+local_particles%n_el(iparticle_kind)
      ENDDO
      CALL structure_factor_allocate(pw_grid%bounds, nnodes, exp_igr)

      ALLOCATE (summe_ef(1:pw_grid%ngpts_cut))
      summe_ef = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      ! Stress Tensor
      IF (do_stress) THEN
         pv_tmp = 0.0_dp
         ALLOCATE (summe_st(3, 1:pw_grid%ngpts_cut))
         summe_st = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      END IF

      ! Defining four_alpha_sq
      four_alpha_sq = 4.0_dp*alpha**2
      dipole_t = 0.0_dp
      q_t = 0.0_dp
      trq_t = 0.0_dp
      ! Zero node count
      node = 0
      DO iparticle_kind = 1, nparticle_kind
         nparticle_local = local_particles%n_el(iparticle_kind)
         DO iparticle_local = 1, nparticle_local
            node = node+1
            iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
            CALL matvec_3x3(vec, cell%h_inv, particle_set(iparticle)%r)
            CALL structure_factor_evaluate(vec, exp_igr%lb, &
                                           exp_igr%ex(:, node), exp_igr%ey(:, node), exp_igr%ez(:, node))

            ! Computing the total charge, dipole and quadrupole trace (if any)
            IF (ANY(task(1, :))) THEN
               q_t = q_t+charges(iparticle)
            END IF
            IF (ANY(task(2, :))) THEN
               dipole_t = dipole_t+dipoles(:, iparticle)
            END IF
            IF (ANY(task(3, :))) THEN
               trq_t = trq_t+quadrupoles(1, 1, iparticle)+ &
                       quadrupoles(2, 2, iparticle)+ &
                       quadrupoles(3, 3, iparticle)
            END IF
         END DO
      END DO

      ! Looping over the positive g-vectors
      DO gpt = 1, pw_grid%ngpts_cut_local
         lp = pw_grid%mapl%pos(pw_grid%g_hat(1, gpt))
         mp = pw_grid%mapm%pos(pw_grid%g_hat(2, gpt))
         np = pw_grid%mapn%pos(pw_grid%g_hat(3, gpt))

         lp = lp+bds(1, 1)
         mp = mp+bds(1, 2)
         np = np+bds(1, 3)

         ! Initializing sum to be used in the energy and force
         node = 0
         DO iparticle_kind = 1, nparticle_kind
            nparticle_local = local_particles%n_el(iparticle_kind)
            DO iparticle_local = 1, nparticle_local
               node = node+1
               iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
               ! Density for energy and forces
               CALL get_atom_factor(atm_factor, pw_grid, gpt, iparticle, task, charges, &
                                    dipoles, quadrupoles)
               summe_tmp = exp_igr%ex(lp, node)*exp_igr%ey(mp, node)*exp_igr%ez(np, node)
               summe_ef(gpt) = summe_ef(gpt)+atm_factor*summe_tmp

               ! Precompute pseudo-density for stress tensor calculation
               IF (do_stress) THEN
                  CALL get_atom_factor_stress(atm_factor_st, pw_grid, gpt, iparticle, task, &
                                              dipoles, quadrupoles)
                  summe_st(1:3, gpt) = summe_st(1:3, gpt)+atm_factor_st(1:3)*summe_tmp
               END IF
            END DO
         END DO
      END DO
      CALL mp_sum(q_t, group)
      CALL mp_sum(dipole_t, group)
      CALL mp_sum(trq_t, group)
      CALL mp_sum(summe_ef, group)
      IF (do_stress) CALL mp_sum(summe_st, group)

      ! Looping over the positive g-vectors
      DO gpt = 1, pw_grid%ngpts_cut_local
         ! computing the potential energy
         lp = pw_grid%mapl%pos(pw_grid%g_hat(1, gpt))
         mp = pw_grid%mapm%pos(pw_grid%g_hat(2, gpt))
         np = pw_grid%mapn%pos(pw_grid%g_hat(3, gpt))

         lp = lp+bds(1, 1)
         mp = mp+bds(1, 2)
         np = np+bds(1, 3)

         IF (pw_grid%gsq(gpt) == 0.0_dp) THEN
            ! G=0 vector for dipole-dipole and charge-quadrupole
            energy = energy+(1.0_dp/6.0_dp)*DOT_PRODUCT(dipole_t, dipole_t) &
                     -(1.0_dp/9.0_dp)*q_t*trq_t
            ! Stress tensor
            IF (do_stress) THEN
               pv_tmp(1, 1) = pv_tmp(1, 1)+(1.0_dp/6.0_dp)*DOT_PRODUCT(dipole_t, dipole_t)
               pv_tmp(2, 2) = pv_tmp(2, 2)+(1.0_dp/6.0_dp)*DOT_PRODUCT(dipole_t, dipole_t)
               pv_tmp(3, 3) = pv_tmp(3, 3)+(1.0_dp/6.0_dp)*DOT_PRODUCT(dipole_t, dipole_t)
            END IF
            ! Corrections for G=0 to potential, field and field gradient
            IF (do_efield .AND. (debug_e_field_en .OR. (.NOT. debug_this_module))) THEN
               ! This term is important and may give problems if one is debugging
               ! VS finite differences since it comes from a residual integral in
               ! the complex plane (cannot be reproduced with finite differences)
               node = 0
               DO iparticle_kind = 1, nparticle_kind
                  nparticle_local = local_particles%n_el(iparticle_kind)
                  DO iparticle_local = 1, nparticle_local
                     node = node+1
                     iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)

                     ! Potential
                     IF (do_efield0) THEN
                        efield0(iparticle) = efield0(iparticle)
                     END IF
                     ! Electrostatic field
                     IF (do_efield1) THEN
                        efield1(1:3, iparticle) = efield1(1:3, iparticle)-(1.0_dp/6.0_dp)*dipole_t(1:3)
                     END IF
                     ! Electrostatic field gradients
                     IF (do_efield2) THEN
                        efield2(1, iparticle) = efield2(1, iparticle)-(1.0_dp/(18.0_dp))*q_t
                        efield2(5, iparticle) = efield2(5, iparticle)-(1.0_dp/(18.0_dp))*q_t
                        efield2(9, iparticle) = efield2(9, iparticle)-(1.0_dp/(18.0_dp))*q_t
                     END IF
                  END DO
               END DO
            END IF
            CYCLE
         END IF
         gauss = (rho0(lp, mp, np)*pw_grid%vol)**2/pw_grid%gsq(gpt)
         factor = gauss*REAL(summe_ef(gpt)*CONJG(summe_ef(gpt)), KIND=dp)
         energy = energy+factor

         IF (do_forces .OR. do_efield) THEN
            node = 0
            DO iparticle_kind = 1, nparticle_kind
               nparticle_local = local_particles%n_el(iparticle_kind)
               DO iparticle_local = 1, nparticle_local
                  node = node+1
                  iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
                  fac = exp_igr%ex(lp, node)*exp_igr%ey(mp, node)*exp_igr%ez(np, node)
                  cnjg_fac = CONJG(fac)

                  ! Forces
                  IF (do_forces) THEN
                     CALL get_atom_factor(atm_factor, pw_grid, gpt, iparticle, task, charges, &
                                          dipoles, quadrupoles)

                     tmp = gauss*AIMAG(summe_ef(gpt)*(cnjg_fac*CONJG(atm_factor)))
                     forces(1, node) = forces(1, node)+tmp*pw_grid%g(1, gpt)
                     forces(2, node) = forces(2, node)+tmp*pw_grid%g(2, gpt)
                     forces(3, node) = forces(3, node)+tmp*pw_grid%g(3, gpt)
                  END IF

                  ! Electric field
                  IF (do_efield) THEN
                     ! Potential
                     IF (do_efield0) THEN
                        efield0(iparticle) = efield0(iparticle)+gauss*REAL(fac*CONJG(summe_ef(gpt)), KIND=dp)
                     END IF
                     ! Electric field
                     IF (do_efield1) THEN
                        tmp = AIMAG(fac*CONJG(summe_ef(gpt)))*gauss
                        efield1(1, iparticle) = efield1(1, iparticle)-tmp*pw_grid%g(1, gpt)
                        efield1(2, iparticle) = efield1(2, iparticle)-tmp*pw_grid%g(2, gpt)
                        efield1(3, iparticle) = efield1(3, iparticle)-tmp*pw_grid%g(3, gpt)
                     END IF
                     ! Electric field gradient
                     IF (do_efield2) THEN
                        tmp_v(1) = REAL(fac*CONJG(summe_ef(gpt)), KIND=dp)*pw_grid%g(1, gpt)*gauss
                        tmp_v(2) = REAL(fac*CONJG(summe_ef(gpt)), KIND=dp)*pw_grid%g(2, gpt)*gauss
                        tmp_v(3) = REAL(fac*CONJG(summe_ef(gpt)), KIND=dp)*pw_grid%g(3, gpt)*gauss

                        efield2(1, iparticle) = efield2(1, iparticle)+tmp_v(1)*pw_grid%g(1, gpt)
                        efield2(2, iparticle) = efield2(2, iparticle)+tmp_v(1)*pw_grid%g(2, gpt)
                        efield2(3, iparticle) = efield2(3, iparticle)+tmp_v(1)*pw_grid%g(3, gpt)
                        efield2(4, iparticle) = efield2(4, iparticle)+tmp_v(2)*pw_grid%g(1, gpt)
                        efield2(5, iparticle) = efield2(5, iparticle)+tmp_v(2)*pw_grid%g(2, gpt)
                        efield2(6, iparticle) = efield2(6, iparticle)+tmp_v(2)*pw_grid%g(3, gpt)
                        efield2(7, iparticle) = efield2(7, iparticle)+tmp_v(3)*pw_grid%g(1, gpt)
                        efield2(8, iparticle) = efield2(8, iparticle)+tmp_v(3)*pw_grid%g(2, gpt)
                        efield2(9, iparticle) = efield2(9, iparticle)+tmp_v(3)*pw_grid%g(3, gpt)
                     END IF
                  END IF
               END DO
            END DO
         END IF

         ! Compute the virial P*V
         IF (do_stress) THEN
            ! The Stress Tensor can be decomposed in two main components.
            ! The first one is just a normal ewald component for reciprocal space
            denom = 1.0_dp/four_alpha_sq+1.0_dp/pw_grid%gsq(gpt)
            pv_tmp(1, 1) = pv_tmp(1, 1)+factor*(1.0_dp-2.0_dp*pw_grid%g(1, gpt)*pw_grid%g(1, gpt)*denom)
            pv_tmp(1, 2) = pv_tmp(1, 2)-factor*(2.0_dp*pw_grid%g(1, gpt)*pw_grid%g(2, gpt)*denom)
            pv_tmp(1, 3) = pv_tmp(1, 3)-factor*(2.0_dp*pw_grid%g(1, gpt)*pw_grid%g(3, gpt)*denom)
            pv_tmp(2, 1) = pv_tmp(2, 1)-factor*(2.0_dp*pw_grid%g(2, gpt)*pw_grid%g(1, gpt)*denom)
            pv_tmp(2, 2) = pv_tmp(2, 2)+factor*(1.0_dp-2.0_dp*pw_grid%g(2, gpt)*pw_grid%g(2, gpt)*denom)
            pv_tmp(2, 3) = pv_tmp(2, 3)-factor*(2.0_dp*pw_grid%g(2, gpt)*pw_grid%g(3, gpt)*denom)
            pv_tmp(3, 1) = pv_tmp(3, 1)-factor*(2.0_dp*pw_grid%g(3, gpt)*pw_grid%g(1, gpt)*denom)
            pv_tmp(3, 2) = pv_tmp(3, 2)-factor*(2.0_dp*pw_grid%g(3, gpt)*pw_grid%g(2, gpt)*denom)
            pv_tmp(3, 3) = pv_tmp(3, 3)+factor*(1.0_dp-2.0_dp*pw_grid%g(3, gpt)*pw_grid%g(3, gpt)*denom)
            ! The second one can be written in the following way
            f0 = 2.0_dp*gauss
            pv_tmp(1, 1) = pv_tmp(1, 1)+f0*pw_grid%g(1, gpt)*REAL(summe_st(1, gpt)*CONJG(summe_ef(gpt)), KIND=dp)
            pv_tmp(1, 2) = pv_tmp(1, 2)+f0*pw_grid%g(1, gpt)*REAL(summe_st(2, gpt)*CONJG(summe_ef(gpt)), KIND=dp)
            pv_tmp(1, 3) = pv_tmp(1, 3)+f0*pw_grid%g(1, gpt)*REAL(summe_st(3, gpt)*CONJG(summe_ef(gpt)), KIND=dp)
            pv_tmp(2, 1) = pv_tmp(2, 1)+f0*pw_grid%g(2, gpt)*REAL(summe_st(1, gpt)*CONJG(summe_ef(gpt)), KIND=dp)
            pv_tmp(2, 2) = pv_tmp(2, 2)+f0*pw_grid%g(2, gpt)*REAL(summe_st(2, gpt)*CONJG(summe_ef(gpt)), KIND=dp)
            pv_tmp(2, 3) = pv_tmp(2, 3)+f0*pw_grid%g(2, gpt)*REAL(summe_st(3, gpt)*CONJG(summe_ef(gpt)), KIND=dp)
            pv_tmp(3, 1) = pv_tmp(3, 1)+f0*pw_grid%g(3, gpt)*REAL(summe_st(1, gpt)*CONJG(summe_ef(gpt)), KIND=dp)
            pv_tmp(3, 2) = pv_tmp(3, 2)+f0*pw_grid%g(3, gpt)*REAL(summe_st(2, gpt)*CONJG(summe_ef(gpt)), KIND=dp)
            pv_tmp(3, 3) = pv_tmp(3, 3)+f0*pw_grid%g(3, gpt)*REAL(summe_st(3, gpt)*CONJG(summe_ef(gpt)), KIND=dp)
         END IF
      END DO
      pref = fourpi/pw_grid%vol
      energy = energy*pref

      CALL structure_factor_deallocate(exp_igr)
      DEALLOCATE (summe_ef)
      IF (do_stress) THEN
         pv_tmp = pv_tmp*pref
         ! Symmetrize the tensor
         pv(1, 1) = pv(1, 1)+pv_tmp(1, 1)
         pv(1, 2) = pv(1, 2)+(pv_tmp(1, 2)+pv_tmp(2, 1))*0.5_dp
         pv(1, 3) = pv(1, 3)+(pv_tmp(1, 3)+pv_tmp(3, 1))*0.5_dp
         pv(2, 1) = pv(1, 2)
         pv(2, 2) = pv(2, 2)+pv_tmp(2, 2)
         pv(2, 3) = pv(2, 3)+(pv_tmp(2, 3)+pv_tmp(3, 2))*0.5_dp
         pv(3, 1) = pv(1, 3)
         pv(3, 2) = pv(2, 3)
         pv(3, 3) = pv(3, 3)+pv_tmp(3, 3)
         DEALLOCATE (summe_st)
      END IF
      IF (do_forces) THEN
         forces = 2.0_dp*forces*pref
      END IF
      IF (do_efield0) THEN
         efield0 = 2.0_dp*efield0*pref
      END IF
      IF (do_efield1) THEN
         efield1 = 2.0_dp*efield1*pref
      END IF
      IF (do_efield2) THEN
         efield2 = 2.0_dp*efield2*pref
      END IF
      CALL timestop(handle)

   END SUBROUTINE ewald_multipole_LR

! **************************************************************************************************
!> \brief Computes the atom factor including charge, dipole and quadrupole
!> \param atm_factor ...
!> \param pw_grid ...
!> \param gpt ...
!> \param iparticle ...
!> \param task ...
!> \param charges ...
!> \param dipoles ...
!> \param quadrupoles ...
!> \par History
!>      none
!> \author Teodoro Laino [tlaino] - 12.2007 - University of Zurich
! **************************************************************************************************
   SUBROUTINE get_atom_factor(atm_factor, pw_grid, gpt, iparticle, task, charges, &
                              dipoles, quadrupoles)
      COMPLEX(KIND=dp), INTENT(OUT)                      :: atm_factor
      TYPE(pw_grid_type), POINTER                        :: pw_grid
      INTEGER, INTENT(IN)                                :: gpt
      INTEGER                                            :: iparticle
      LOGICAL, DIMENSION(3, 3), INTENT(IN)               :: task
      REAL(KIND=dp), DIMENSION(:), OPTIONAL, POINTER     :: charges
      REAL(KIND=dp), DIMENSION(:, :), OPTIONAL, POINTER  :: dipoles
      REAL(KIND=dp), DIMENSION(:, :, :), OPTIONAL, &
         POINTER                                         :: quadrupoles

      CHARACTER(len=*), PARAMETER :: routineN = 'get_atom_factor', &
         routineP = moduleN//':'//routineN

      COMPLEX(KIND=dp)                                   :: tmp
      INTEGER                                            :: i, j

      atm_factor = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      IF (task(1, 1)) THEN
         ! Charge
         atm_factor = atm_factor+charges(iparticle)
      END IF
      IF (task(2, 2)) THEN
         ! Dipole
         tmp = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         DO i = 1, 3
            tmp = tmp+dipoles(i, iparticle)*pw_grid%g(i, gpt)
         END DO
         atm_factor = atm_factor+tmp*CMPLX(0.0_dp, -1.0_dp, KIND=dp)
      END IF
      IF (task(3, 3)) THEN
         ! Quadrupole
         tmp = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
         DO i = 1, 3
            DO j = 1, 3
               tmp = tmp+quadrupoles(j, i, iparticle)*pw_grid%g(j, gpt)*pw_grid%g(i, gpt)
            END DO
         END DO
         atm_factor = atm_factor-1.0_dp/3.0_dp*tmp
      END IF

   END SUBROUTINE get_atom_factor

! **************************************************************************************************
!> \brief Computes the atom factor including charge, dipole and quadrupole
!> \param atm_factor ...
!> \param pw_grid ...
!> \param gpt ...
!> \param iparticle ...
!> \param task ...
!> \param dipoles ...
!> \param quadrupoles ...
!> \par History
!>      none
!> \author Teodoro Laino [tlaino] - 12.2007 - University of Zurich
! **************************************************************************************************
   SUBROUTINE get_atom_factor_stress(atm_factor, pw_grid, gpt, iparticle, task, &
                                     dipoles, quadrupoles)
      COMPLEX(KIND=dp), INTENT(OUT)                      :: atm_factor(3)
      TYPE(pw_grid_type), POINTER                        :: pw_grid
      INTEGER, INTENT(IN)                                :: gpt
      INTEGER                                            :: iparticle
      LOGICAL, DIMENSION(3, 3), INTENT(IN)               :: task
      REAL(KIND=dp), DIMENSION(:, :), OPTIONAL, POINTER  :: dipoles
      REAL(KIND=dp), DIMENSION(:, :, :), OPTIONAL, &
         POINTER                                         :: quadrupoles

      CHARACTER(len=*), PARAMETER :: routineN = 'get_atom_factor_stress', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i

      atm_factor = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      IF (ANY(task(2, :))) THEN
         ! Dipole
         atm_factor = dipoles(:, iparticle)*CMPLX(0.0_dp, -1.0_dp, KIND=dp)
      END IF
      IF (ANY(task(3, :))) THEN
         ! Quadrupole
         DO i = 1, 3
            atm_factor(1) = atm_factor(1)-1.0_dp/3.0_dp* &
                            (quadrupoles(1, i, iparticle)*pw_grid%g(i, gpt)+ &
                             quadrupoles(i, 1, iparticle)*pw_grid%g(i, gpt))
            atm_factor(2) = atm_factor(2)-1.0_dp/3.0_dp* &
                            (quadrupoles(2, i, iparticle)*pw_grid%g(i, gpt)+ &
                             quadrupoles(i, 2, iparticle)*pw_grid%g(i, gpt))
            atm_factor(3) = atm_factor(3)-1.0_dp/3.0_dp* &
                            (quadrupoles(3, i, iparticle)*pw_grid%g(i, gpt)+ &
                             quadrupoles(i, 3, iparticle)*pw_grid%g(i, gpt))
         END DO
      END IF

   END SUBROUTINE get_atom_factor_stress

! **************************************************************************************************
!> \brief Computes the self interaction from g-space and the neutralizing background
!>     when using multipoles
!> \param ewald_env ...
!> \param cell ...
!> \param local_particles ...
!> \param e_self ...
!> \param e_neut ...
!> \param task ...
!> \param do_efield ...
!> \param radii ...
!> \param charges ...
!> \param dipoles ...
!> \param quadrupoles ...
!> \param efield0 ...
!> \param efield1 ...
!> \param efield2 ...
!> \author Teodoro Laino [tlaino] - University of Zurich - 12.2007
! **************************************************************************************************
   SUBROUTINE ewald_multipole_self(ewald_env, cell, local_particles, e_self, &
                                   e_neut, task, do_efield, radii, charges, dipoles, quadrupoles, efield0, &
                                   efield1, efield2)
      TYPE(ewald_environment_type), POINTER              :: ewald_env
      TYPE(cell_type), INTENT(IN)                        :: cell
      TYPE(distribution_1d_type), POINTER                :: local_particles
      REAL(KIND=dp), INTENT(OUT)                         :: e_self, e_neut
      LOGICAL, DIMENSION(3, 3), INTENT(IN)               :: task
      LOGICAL, INTENT(IN)                                :: do_efield
      REAL(KIND=dp), DIMENSION(:), OPTIONAL, POINTER     :: radii, charges
      REAL(KIND=dp), DIMENSION(:, :), OPTIONAL, POINTER  :: dipoles
      REAL(KIND=dp), DIMENSION(:, :, :), OPTIONAL, &
         POINTER                                         :: quadrupoles
      REAL(KIND=dp), DIMENSION(:), POINTER               :: efield0
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: efield1, efield2

      CHARACTER(len=*), PARAMETER :: routineN = 'ewald_multipole_self', &
         routineP = moduleN//':'//routineN
      REAL(KIND=dp), PARAMETER                           :: f23 = 2.0_dp/3.0_dp, &
                                                            f415 = 4.0_dp/15.0_dp

      INTEGER                                            :: ewald_type, group, i, iparticle, &
                                                            iparticle_kind, iparticle_local, j, &
                                                            nparticle_local
      LOGICAL                                            :: do_efield0, do_efield1, do_efield2, &
                                                            lradii
      REAL(KIND=dp)                                      :: alpha, ch_qu_self, ch_qu_self_tmp, &
                                                            dipole_self, fac1, fac2, fac3, fac4, &
                                                            q, q_neutg, q_self, q_sum, qu_qu_self, &
                                                            radius

      CALL ewald_env_get(ewald_env, ewald_type=ewald_type, alpha=alpha, &
                         group=group)

      do_efield0 = do_efield .AND. ASSOCIATED(efield0)
      do_efield1 = do_efield .AND. ASSOCIATED(efield1)
      do_efield2 = do_efield .AND. ASSOCIATED(efield2)
      q_self = 0.0_dp
      q_sum = 0.0_dp
      dipole_self = 0.0_dp
      ch_qu_self = 0.0_dp
      qu_qu_self = 0.0_dp
      fac1 = 2.0_dp*alpha*oorootpi
      fac2 = 6.0_dp*(f23**2)*(alpha**3)*oorootpi
      fac3 = (2.0_dp*oorootpi)*f23*alpha**3
      fac4 = (4.0_dp*oorootpi)*f415*alpha**5
      lradii = PRESENT(radii)
      radius = 0.0_dp
      q_neutg = 0.0_dp
      DO iparticle_kind = 1, SIZE(local_particles%n_el)
         nparticle_local = local_particles%n_el(iparticle_kind)
         DO iparticle_local = 1, nparticle_local
            iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
            IF (ANY(task(1, :))) THEN
               ! Charge - Charge
               q = charges(iparticle)
               IF (lradii) radius = radii(iparticle)
               IF (radius > 0) THEN
                  q_neutg = q_neutg+2.0_dp*q*radius**2
               END IF
               q_self = q_self+q*q
               q_sum = q_sum+q
               ! Potential
               IF (do_efield0) THEN
                  efield0(iparticle) = efield0(iparticle)-q*fac1
               END IF

               IF (task(1, 3)) THEN
                  ! Charge - Quadrupole
                  ch_qu_self_tmp = 0.0_dp
                  DO i = 1, 3
                     ch_qu_self_tmp = ch_qu_self_tmp+quadrupoles(i, i, iparticle)*q
                  END DO
                  ch_qu_self = ch_qu_self+ch_qu_self_tmp
                  ! Electric Field Gradient
                  IF (do_efield2) THEN
                     efield2(1, iparticle) = efield2(1, iparticle)+fac2*q
                     efield2(5, iparticle) = efield2(5, iparticle)+fac2*q
                     efield2(9, iparticle) = efield2(9, iparticle)+fac2*q
                  END IF
               END IF
            END IF
            IF (ANY(task(2, :))) THEN
               ! Dipole - Dipole
               DO i = 1, 3
                  dipole_self = dipole_self+dipoles(i, iparticle)**2
               END DO
               ! Electric Field
               IF (do_efield1) THEN
                  efield1(1, iparticle) = efield1(1, iparticle)+fac3*dipoles(1, iparticle)
                  efield1(2, iparticle) = efield1(2, iparticle)+fac3*dipoles(2, iparticle)
                  efield1(3, iparticle) = efield1(3, iparticle)+fac3*dipoles(3, iparticle)
               END IF
            END IF
            IF (ANY(task(3, :))) THEN
               ! Quadrupole - Quadrupole
               DO i = 1, 3
                  DO j = 1, 3
                     qu_qu_self = qu_qu_self+quadrupoles(j, i, iparticle)**2
                  END DO
               END DO
               ! Electric Field Gradient
               IF (do_efield2) THEN
                  efield2(1, iparticle) = efield2(1, iparticle)+fac4*quadrupoles(1, 1, iparticle)
                  efield2(2, iparticle) = efield2(2, iparticle)+fac4*quadrupoles(2, 1, iparticle)
                  efield2(3, iparticle) = efield2(3, iparticle)+fac4*quadrupoles(3, 1, iparticle)
                  efield2(4, iparticle) = efield2(4, iparticle)+fac4*quadrupoles(1, 2, iparticle)
                  efield2(5, iparticle) = efield2(5, iparticle)+fac4*quadrupoles(2, 2, iparticle)
                  efield2(6, iparticle) = efield2(6, iparticle)+fac4*quadrupoles(3, 2, iparticle)
                  efield2(7, iparticle) = efield2(7, iparticle)+fac4*quadrupoles(1, 3, iparticle)
                  efield2(8, iparticle) = efield2(8, iparticle)+fac4*quadrupoles(2, 3, iparticle)
                  efield2(9, iparticle) = efield2(9, iparticle)+fac4*quadrupoles(3, 3, iparticle)
               END IF
            END IF
         END DO
      END DO

      CALL mp_sum(q_neutg, group)
      CALL mp_sum(q_self, group)
      CALL mp_sum(q_sum, group)
      CALL mp_sum(dipole_self, group)
      CALL mp_sum(ch_qu_self, group)
      CALL mp_sum(qu_qu_self, group)

      e_self = -(q_self+f23*(dipole_self-f23*ch_qu_self+f415*qu_qu_self*alpha**2)*alpha**2)*alpha*oorootpi
      fac1 = pi/(2.0_dp*cell%deth)
      e_neut = -q_sum*fac1*(q_sum/alpha**2-q_neutg)

      ! Correcting Potential for the neutralizing background charge
      DO iparticle_kind = 1, SIZE(local_particles%n_el)
         nparticle_local = local_particles%n_el(iparticle_kind)
         DO iparticle_local = 1, nparticle_local
            iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
            IF (ANY(task(1, :))) THEN
               ! Potential energy
               IF (do_efield0) THEN
                  efield0(iparticle) = efield0(iparticle)-q_sum*2.0_dp*fac1/alpha**2
                  IF (lradii) radius = radii(iparticle)
                  IF (radius > 0) THEN
                     q = charges(iparticle)
                     efield0(iparticle) = efield0(iparticle)+fac1*radius**2*(q_sum+q)
                  END IF
               END IF
            END IF
         END DO
      END DO
   END SUBROUTINE ewald_multipole_self

! **************************************************************************************************
!> \brief ...
!> \param iw ...
!> \param e_gspace ...
!> \param e_rspace ...
!> \param e_bonded ...
!> \param e_self ...
!> \param e_neut ...
!> \author Teodoro Laino [tlaino] - University of Zurich - 12.2007
! **************************************************************************************************
   SUBROUTINE ewald_multipole_print(iw, e_gspace, e_rspace, e_bonded, e_self, e_neut)

      INTEGER, INTENT(IN)                                :: iw
      REAL(KIND=dp), INTENT(IN)                          :: e_gspace, e_rspace, e_bonded, e_self, &
                                                            e_neut

      CHARACTER(len=*), PARAMETER :: routineN = 'ewald_multipole_print', &
         routineP = moduleN//':'//routineN

      IF (iw > 0) THEN
         WRITE (iw, '( A, A )') ' *********************************', &
            '**********************************************'
         WRITE (iw, '( A, A, T35, A, T56, E25.15 )') ' INITIAL GSPACE ENERGY', &
            '[hartree]', '= ', e_gspace
         WRITE (iw, '( A, A, T35, A, T56, E25.15 )') ' INITIAL RSPACE ENERGY', &
            '[hartree]', '= ', e_rspace
         WRITE (iw, '( A, A, T35, A, T56, E25.15 )') ' BONDED CORRECTION', &
            '[hartree]', '= ', e_bonded
         WRITE (iw, '( A, A, T35, A, T56, E25.15 )') ' SELF ENERGY CORRECTION', &
            '[hartree]', '= ', e_self
         WRITE (iw, '( A, A, T35, A, T56, E25.15 )') ' NEUTRALIZ. BCKGR. ENERGY', &
            '[hartree]', '= ', e_neut
         WRITE (iw, '( A, A, T35, A, T56, E25.15 )') ' TOTAL ELECTROSTATIC EN.', &
            '[hartree]', '= ', e_rspace+e_bonded+e_gspace+e_self+e_neut
         WRITE (iw, '( A, A )') ' *********************************', &
            '**********************************************'
      END IF
   END SUBROUTINE ewald_multipole_print

! **************************************************************************************************
!> \brief  Debug routines for multipoles
!> \param ewald_env ...
!> \param ewald_pw ...
!> \param nonbond_env ...
!> \param cell ...
!> \param particle_set ...
!> \param local_particles ...
!> \param iw ...
!> \param debug_r_space ...
!> \date   05.2008
!> \author Teodoro Laino [tlaino] - University of Zurich - 05.2008
! **************************************************************************************************
   SUBROUTINE debug_ewald_multipoles(ewald_env, ewald_pw, nonbond_env, cell, &
                                     particle_set, local_particles, iw, debug_r_space)
      TYPE charge_mono_type
         REAL(KIND=dp), DIMENSION(:), &
            POINTER                          :: charge
         REAL(KIND=dp), DIMENSION(:, :), &
            POINTER                          :: pos
      END TYPE charge_mono_type
      TYPE multi_charge_type
         TYPE(charge_mono_type), DIMENSION(:), &
            POINTER                          :: charge_typ
      END TYPE multi_charge_type
      TYPE(ewald_environment_type), POINTER    :: ewald_env
      TYPE(ewald_pw_type), POINTER             :: ewald_pw
      TYPE(fist_nonbond_env_type), POINTER     :: nonbond_env
      TYPE(cell_type), POINTER                 :: cell
      TYPE(particle_type), DIMENSION(:), &
         POINTER                                :: particle_set
      TYPE(distribution_1d_type), POINTER      :: local_particles
      INTEGER, INTENT(IN)                      :: iw
      LOGICAL, INTENT(IN)                      :: debug_r_space

      CHARACTER(len=*), PARAMETER :: routineN = 'debug_ewald_multipoles', &
                                     routineP = "ewalds_multipole_debug"//':'//routineN

      INTEGER                                  :: nparticles
      LOGICAL, DIMENSION(3)                    :: task
      REAL(KIND=dp)                            :: e_neut, e_self, g_energy, &
                                                  r_energy, debug_energy
      REAL(KIND=dp), POINTER, DIMENSION(:)     :: charges
      REAL(KIND=dp), POINTER, &
         DIMENSION(:, :)                     :: dipoles, g_forces, g_pv, &
                                                r_forces, r_pv, e_field1, &
                                                e_field2
      REAL(KIND=dp), POINTER, &
         DIMENSION(:, :, :)                  :: quadrupoles
      TYPE(rng_stream_type), POINTER           :: random_stream
      TYPE(multi_charge_type), DIMENSION(:), &
         POINTER                             :: multipoles

      NULLIFY (random_stream, multipoles, charges, dipoles, g_forces, g_pv, &
               r_forces, r_pv, e_field1, e_field2)
      CALL create_rng_stream(random_stream, name="DEBUG_EWALD_MULTIPOLE", &
                             distribution_type=UNIFORM)
      ! check:  charge - charge
      task = .FALSE.
      nparticles = SIZE(particle_set)

      ! Allocate charges, dipoles, quadrupoles
      ALLOCATE (charges(nparticles))
      ALLOCATE (dipoles(3, nparticles))
      ALLOCATE (quadrupoles(3, 3, nparticles))

      ! Allocate arrays for forces
      ALLOCATE (r_forces(3, nparticles))
      ALLOCATE (g_forces(3, nparticles))
      ALLOCATE (e_field1(3, nparticles))
      ALLOCATE (e_field2(3, nparticles))
      ALLOCATE (g_pv(3, 3))
      ALLOCATE (r_pv(3, 3))

      ! Debug CHARGES-CHARGES
      task(1) = .TRUE.
      charges = 0.0_dp
      dipoles = 0.0_dp
      quadrupoles = 0.0_dp
      r_forces = 0.0_dp
      g_forces = 0.0_dp
      e_field1 = 0.0_dp
      e_field2 = 0.0_dp
      g_pv = 0.0_dp
      r_pv = 0.0_dp
      g_energy = 0.0_dp
      r_energy = 0.0_dp
      e_neut = 0.0_dp
      e_self = 0.0_dp

      CALL create_multi_type(multipoles, nparticles, 1, nparticles/2, "CHARGE", echarge=-1.0_dp, &
                             random_stream=random_stream, charges=charges)
      CALL create_multi_type(multipoles, nparticles, nparticles/2+1, nparticles, "CHARGE", echarge=1.0_dp, &
                             random_stream=random_stream, charges=charges)
      CALL debug_ewald_multipole_low(particle_set, cell, nonbond_env, multipoles, debug_energy, &
                                     debug_r_space)

      WRITE (iw, *) "DEBUG ENERGY (CHARGE-CHARGE): ", debug_energy
      CALL ewald_multipole_evaluate(ewald_env, ewald_pw, nonbond_env, cell, &
                                    particle_set, local_particles, g_energy, r_energy, e_neut, e_self, &
                                    task, do_correction_bonded=.FALSE., do_forces=.TRUE., do_stress=.TRUE., do_efield=.FALSE., &
                                    charges=charges, dipoles=dipoles, quadrupoles=quadrupoles, forces_local=g_forces, &
                                    forces_glob=r_forces, pv_local=g_pv, pv_glob=r_pv, iw=iw, do_debug=.FALSE.)
      CALL release_multi_type(multipoles)

      ! Debug CHARGES-DIPOLES
      task(1) = .TRUE.
      task(2) = .TRUE.
      charges = 0.0_dp
      dipoles = 0.0_dp
      quadrupoles = 0.0_dp
      r_forces = 0.0_dp
      g_forces = 0.0_dp
      e_field1 = 0.0_dp
      e_field2 = 0.0_dp
      g_pv = 0.0_dp
      r_pv = 0.0_dp
      g_energy = 0.0_dp
      r_energy = 0.0_dp
      e_neut = 0.0_dp
      e_self = 0.0_dp

      CALL create_multi_type(multipoles, nparticles, 1, nparticles/2, "CHARGE", echarge=-1.0_dp, &
                             random_stream=random_stream, charges=charges)
      CALL create_multi_type(multipoles, nparticles, nparticles/2+1, nparticles, "DIPOLE", echarge=0.5_dp, &
                             random_stream=random_stream, dipoles=dipoles)
      WRITE (iw, '("CHARGES",F15.9)') charges
      WRITE (iw, '("DIPOLES",3F15.9)') dipoles
      CALL debug_ewald_multipole_low(particle_set, cell, nonbond_env, multipoles, debug_energy, &
                                     debug_r_space)

      WRITE (iw, *) "DEBUG ENERGY (CHARGE-DIPOLE): ", debug_energy
      CALL ewald_multipole_evaluate(ewald_env, ewald_pw, nonbond_env, cell, &
                                    particle_set, local_particles, g_energy, r_energy, e_neut, e_self, &
                                    task, do_correction_bonded=.FALSE., do_forces=.TRUE., do_stress=.TRUE., do_efield=.FALSE., &
                                    charges=charges, dipoles=dipoles, quadrupoles=quadrupoles, forces_local=g_forces, &
                                    forces_glob=r_forces, pv_local=g_pv, pv_glob=r_pv, iw=iw, do_debug=.FALSE.)
      CALL release_multi_type(multipoles)

      ! Debug DIPOLES-DIPOLES
      task(2) = .TRUE.
      charges = 0.0_dp
      dipoles = 0.0_dp
      quadrupoles = 0.0_dp
      r_forces = 0.0_dp
      g_forces = 0.0_dp
      e_field1 = 0.0_dp
      e_field2 = 0.0_dp
      g_pv = 0.0_dp
      r_pv = 0.0_dp
      g_energy = 0.0_dp
      r_energy = 0.0_dp
      e_neut = 0.0_dp
      e_self = 0.0_dp

      CALL create_multi_type(multipoles, nparticles, 1, nparticles/2, "DIPOLE", echarge=10000.0_dp, &
                             random_stream=random_stream, dipoles=dipoles)
      CALL create_multi_type(multipoles, nparticles, nparticles/2+1, nparticles, "DIPOLE", echarge=20000._dp, &
                             random_stream=random_stream, dipoles=dipoles)
      WRITE (iw, '("DIPOLES",3F15.9)') dipoles
      CALL debug_ewald_multipole_low(particle_set, cell, nonbond_env, multipoles, debug_energy, &
                                     debug_r_space)

      WRITE (iw, *) "DEBUG ENERGY (DIPOLE-DIPOLE): ", debug_energy
      CALL ewald_multipole_evaluate(ewald_env, ewald_pw, nonbond_env, cell, &
                                    particle_set, local_particles, g_energy, r_energy, e_neut, e_self, &
                                    task, do_correction_bonded=.FALSE., do_forces=.TRUE., do_stress=.TRUE., do_efield=.FALSE., &
                                    charges=charges, dipoles=dipoles, quadrupoles=quadrupoles, forces_local=g_forces, &
                                    forces_glob=r_forces, pv_local=g_pv, pv_glob=r_pv, iw=iw, do_debug=.FALSE.)
      CALL release_multi_type(multipoles)

      ! Debug CHARGES-QUADRUPOLES
      task(1) = .TRUE.
      task(3) = .TRUE.
      charges = 0.0_dp
      dipoles = 0.0_dp
      quadrupoles = 0.0_dp
      r_forces = 0.0_dp
      g_forces = 0.0_dp
      e_field1 = 0.0_dp
      e_field2 = 0.0_dp
      g_pv = 0.0_dp
      r_pv = 0.0_dp
      g_energy = 0.0_dp
      r_energy = 0.0_dp
      e_neut = 0.0_dp
      e_self = 0.0_dp

      CALL create_multi_type(multipoles, nparticles, 1, nparticles/2, "CHARGE", echarge=-1.0_dp, &
                             random_stream=random_stream, charges=charges)
      CALL create_multi_type(multipoles, nparticles, nparticles/2+1, nparticles, "QUADRUPOLE", echarge=10.0_dp, &
                             random_stream=random_stream, quadrupoles=quadrupoles)
      WRITE (iw, '("CHARGES",F15.9)') charges
      WRITE (iw, '("QUADRUPOLES",9F15.9)') quadrupoles
      CALL debug_ewald_multipole_low(particle_set, cell, nonbond_env, multipoles, debug_energy, &
                                     debug_r_space)

      WRITE (iw, *) "DEBUG ENERGY (CHARGE-QUADRUPOLE): ", debug_energy
      CALL ewald_multipole_evaluate(ewald_env, ewald_pw, nonbond_env, cell, &
                                    particle_set, local_particles, g_energy, r_energy, e_neut, e_self, &
                                    task, do_correction_bonded=.FALSE., do_forces=.TRUE., do_stress=.TRUE., do_efield=.FALSE., &
                                    charges=charges, dipoles=dipoles, quadrupoles=quadrupoles, forces_local=g_forces, &
                                    forces_glob=r_forces, pv_local=g_pv, pv_glob=r_pv, iw=iw, do_debug=.FALSE.)
      CALL release_multi_type(multipoles)

      ! Debug DIPOLES-QUADRUPOLES
      task(2) = .TRUE.
      task(3) = .TRUE.
      charges = 0.0_dp
      dipoles = 0.0_dp
      quadrupoles = 0.0_dp
      r_forces = 0.0_dp
      g_forces = 0.0_dp
      e_field1 = 0.0_dp
      e_field2 = 0.0_dp
      g_pv = 0.0_dp
      r_pv = 0.0_dp
      g_energy = 0.0_dp
      r_energy = 0.0_dp
      e_neut = 0.0_dp
      e_self = 0.0_dp

      CALL create_multi_type(multipoles, nparticles, 1, nparticles/2, "DIPOLE", echarge=10000.0_dp, &
                             random_stream=random_stream, dipoles=dipoles)
      CALL create_multi_type(multipoles, nparticles, nparticles/2+1, nparticles, "QUADRUPOLE", echarge=10000.0_dp, &
                             random_stream=random_stream, quadrupoles=quadrupoles)
      WRITE (iw, '("DIPOLES",3F15.9)') dipoles
      WRITE (iw, '("QUADRUPOLES",9F15.9)') quadrupoles
      CALL debug_ewald_multipole_low(particle_set, cell, nonbond_env, multipoles, debug_energy, &
                                     debug_r_space)

      WRITE (iw, *) "DEBUG ENERGY (DIPOLE-QUADRUPOLE): ", debug_energy
      CALL ewald_multipole_evaluate(ewald_env, ewald_pw, nonbond_env, cell, &
                                    particle_set, local_particles, g_energy, r_energy, e_neut, e_self, &
                                    task, do_correction_bonded=.FALSE., do_forces=.TRUE., do_stress=.TRUE., do_efield=.FALSE., &
                                    charges=charges, dipoles=dipoles, quadrupoles=quadrupoles, forces_local=g_forces, &
                                    forces_glob=r_forces, pv_local=g_pv, pv_glob=r_pv, iw=iw, do_debug=.FALSE.)
      CALL release_multi_type(multipoles)

      ! Debug QUADRUPOLES-QUADRUPOLES
      task(3) = .TRUE.
      charges = 0.0_dp
      dipoles = 0.0_dp
      quadrupoles = 0.0_dp
      r_forces = 0.0_dp
      g_forces = 0.0_dp
      e_field1 = 0.0_dp
      e_field2 = 0.0_dp
      g_pv = 0.0_dp
      r_pv = 0.0_dp
      g_energy = 0.0_dp
      r_energy = 0.0_dp
      e_neut = 0.0_dp
      e_self = 0.0_dp

      CALL create_multi_type(multipoles, nparticles, 1, nparticles/2, "QUADRUPOLE", echarge=-20000.0_dp, &
                             random_stream=random_stream, quadrupoles=quadrupoles)
      CALL create_multi_type(multipoles, nparticles, nparticles/2+1, nparticles, "QUADRUPOLE", echarge=10000.0_dp, &
                             random_stream=random_stream, quadrupoles=quadrupoles)
      WRITE (iw, '("QUADRUPOLES",9F15.9)') quadrupoles
      CALL debug_ewald_multipole_low(particle_set, cell, nonbond_env, multipoles, debug_energy, &
                                     debug_r_space)

      WRITE (iw, *) "DEBUG ENERGY (QUADRUPOLE-QUADRUPOLE): ", debug_energy
      CALL ewald_multipole_evaluate(ewald_env, ewald_pw, nonbond_env, cell, &
                                    particle_set, local_particles, g_energy, r_energy, e_neut, e_self, &
                                    task, do_correction_bonded=.FALSE., do_forces=.TRUE., do_stress=.TRUE., do_efield=.FALSE., &
                                    charges=charges, dipoles=dipoles, quadrupoles=quadrupoles, forces_local=g_forces, &
                                    forces_glob=r_forces, pv_local=g_pv, pv_glob=r_pv, iw=iw, do_debug=.FALSE.)
      CALL release_multi_type(multipoles)

      DEALLOCATE (charges)
      DEALLOCATE (dipoles)
      DEALLOCATE (quadrupoles)
      DEALLOCATE (r_forces)
      DEALLOCATE (g_forces)
      DEALLOCATE (e_field1)
      DEALLOCATE (e_field2)
      DEALLOCATE (g_pv)
      DEALLOCATE (r_pv)
      CALL delete_rng_stream(random_stream)

   CONTAINS
! **************************************************************************************************
!> \brief  Debug routines for multipoles - low level - charge interactions
!> \param particle_set ...
!> \param cell ...
!> \param nonbond_env ...
!> \param multipoles ...
!> \param energy ...
!> \param debug_r_space ...
!> \date   05.2008
!> \author Teodoro Laino [tlaino] - University of Zurich - 05.2008
! **************************************************************************************************
      SUBROUTINE debug_ewald_multipole_low(particle_set, cell, nonbond_env, multipoles, &
                                           energy, debug_r_space)
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(fist_nonbond_env_type), POINTER               :: nonbond_env
      TYPE(multi_charge_type), DIMENSION(:), POINTER     :: multipoles
      REAL(KIND=dp), INTENT(OUT)                         :: energy
      LOGICAL, INTENT(IN)                                :: debug_r_space

      CHARACTER(len=*), PARAMETER :: routineN = 'debug_ewald_multipole_low', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: atom_a, atom_b, icell, iend, igrp, &
                                                            ikind, ilist, ipair, istart, jcell, &
                                                            jkind, k, k1, kcell, l, l1, ncells, &
                                                            nkinds, npairs
      INTEGER, DIMENSION(:, :), POINTER                  :: list
      REAL(KIND=dp)                                      :: fac_ij, q, r, rab2, rab2_max
      REAL(KIND=dp), DIMENSION(3)                        :: cell_v, cvi, rab, rab0, rm
      TYPE(fist_neighbor_type), POINTER                  :: nonbonded
      TYPE(neighbor_kind_pairs_type), POINTER            :: neighbor_kind_pair
      TYPE(pos_type), DIMENSION(:), POINTER              :: r_last_update, r_last_update_pbc

         energy = 0.0_dp
         CALL fist_nonbond_env_get(nonbond_env, nonbonded=nonbonded, natom_types=nkinds, &
                                   r_last_update=r_last_update, r_last_update_pbc=r_last_update_pbc)
         rab2_max = HUGE(0.0_dp)
         IF (debug_r_space) THEN
            ! This debugs the real space part of the multipole Ewald summation scheme
            ! Starting the force loop
            Lists: DO ilist = 1, nonbonded%nlists
               neighbor_kind_pair => nonbonded%neighbor_kind_pairs(ilist)
               npairs = neighbor_kind_pair%npairs
               IF (npairs == 0) CYCLE
               list => neighbor_kind_pair%list
               cvi = neighbor_kind_pair%cell_vector
               CALL matvec_3x3(cell_v, cell%hmat, cvi)
               Kind_Group_Loop: DO igrp = 1, neighbor_kind_pair%ngrp_kind
                  istart = neighbor_kind_pair%grp_kind_start(igrp)
                  iend = neighbor_kind_pair%grp_kind_end(igrp)
                  ikind = neighbor_kind_pair%ij_kind(1, igrp)
                  jkind = neighbor_kind_pair%ij_kind(2, igrp)
                  Pairs: DO ipair = istart, iend
                     fac_ij = 1.0_dp
                     atom_a = list(1, ipair)
                     atom_b = list(2, ipair)
                     IF (atom_a == atom_b) fac_ij = 0.5_dp
                     rab = r_last_update_pbc(atom_b)%r-r_last_update_pbc(atom_a)%r
                     rab = rab+cell_v
                     rab2 = rab(1)**2+rab(2)**2+rab(3)**2
                     IF (rab2 <= rab2_max) THEN

                        DO k = 1, SIZE(multipoles(atom_a)%charge_typ)
                           DO k1 = 1, SIZE(multipoles(atom_a)%charge_typ(k)%charge)

                              DO l = 1, SIZE(multipoles(atom_b)%charge_typ)
                                 DO l1 = 1, SIZE(multipoles(atom_b)%charge_typ(l)%charge)

                                    rm = rab+multipoles(atom_b)%charge_typ(l)%pos(:, l1)-multipoles(atom_a)%charge_typ(k)%pos(:, k1)
                                    r = SQRT(DOT_PRODUCT(rm, rm))
                                    q = multipoles(atom_b)%charge_typ(l)%charge(l1)*multipoles(atom_a)%charge_typ(k)%charge(k1)
                                    energy = energy+q/r*fac_ij
                                 END DO
                              END DO

                           END DO
                        END DO

                     END IF
                  END DO Pairs
               END DO Kind_Group_Loop
            END DO Lists
         ELSE
            ncells = 6
            !Debugs the sum of real + space terms.. (Charge-Charge and Charge-Dipole should be anyway wrong but
            !all the other terms should be correct)
            DO atom_a = 1, SIZE(particle_set)
            DO atom_b = atom_a, SIZE(particle_set)
               fac_ij = 1.0_dp
               IF (atom_a == atom_b) fac_ij = 0.5_dp
               rab0 = r_last_update_pbc(atom_b)%r-r_last_update_pbc(atom_a)%r
               ! Loop over cells
               DO icell = -ncells, ncells
               DO jcell = -ncells, ncells
               DO kcell = -ncells, ncells
                  cell_v = MATMUL(cell%hmat, REAL((/icell, jcell, kcell/), KIND=dp))
                  IF (ALL(cell_v == 0.0_dp) .AND. (atom_a == atom_b)) CYCLE
                  rab = rab0+cell_v
                  rab2 = rab(1)**2+rab(2)**2+rab(3)**2
                  IF (rab2 <= rab2_max) THEN

                     DO k = 1, SIZE(multipoles(atom_a)%charge_typ)
                        DO k1 = 1, SIZE(multipoles(atom_a)%charge_typ(k)%charge)

                           DO l = 1, SIZE(multipoles(atom_b)%charge_typ)
                              DO l1 = 1, SIZE(multipoles(atom_b)%charge_typ(l)%charge)

                                 rm = rab+multipoles(atom_b)%charge_typ(l)%pos(:, l1)-multipoles(atom_a)%charge_typ(k)%pos(:, k1)
                                 r = SQRT(DOT_PRODUCT(rm, rm))
                                 q = multipoles(atom_b)%charge_typ(l)%charge(l1)*multipoles(atom_a)%charge_typ(k)%charge(k1)
                                 energy = energy+q/r*fac_ij
                              END DO
                           END DO

                        END DO
                     END DO

                  END IF
               END DO
               END DO
               END DO
            END DO
            END DO
         END IF
      END SUBROUTINE debug_ewald_multipole_low

! **************************************************************************************************
!> \brief  create multi_type for multipoles
!> \param multipoles ...
!> \param idim ...
!> \param istart ...
!> \param iend ...
!> \param label ...
!> \param echarge ...
!> \param random_stream ...
!> \param charges ...
!> \param dipoles ...
!> \param quadrupoles ...
!> \date   05.2008
!> \author Teodoro Laino [tlaino] - University of Zurich - 05.2008
! **************************************************************************************************
      SUBROUTINE create_multi_type(multipoles, idim, istart, iend, label, echarge, &
                                   random_stream, charges, dipoles, quadrupoles)
      TYPE(multi_charge_type), DIMENSION(:), POINTER     :: multipoles
      INTEGER, INTENT(IN)                                :: idim, istart, iend
      CHARACTER(LEN=*), INTENT(IN)                       :: label
      REAL(KIND=dp), INTENT(IN)                          :: echarge
      TYPE(rng_stream_type), POINTER                     :: random_stream
      REAL(KIND=dp), DIMENSION(:), OPTIONAL, POINTER     :: charges
      REAL(KIND=dp), DIMENSION(:, :), OPTIONAL, POINTER  :: dipoles
      REAL(KIND=dp), DIMENSION(:, :, :), OPTIONAL, &
         POINTER                                         :: quadrupoles

      CHARACTER(len=*), PARAMETER :: routineN = 'create_multi_type', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, isize, k, l, m
      REAL(KIND=dp)                                      :: dx, r2, rvec(3), rvec1(3), rvec2(3)

         IF (ASSOCIATED(multipoles)) THEN
            CPASSERT(SIZE(multipoles) == idim)
         ELSE
            ALLOCATE (multipoles(idim))
            DO i = 1, idim
               NULLIFY (multipoles(i)%charge_typ)
            END DO
         END IF
         DO i = istart, iend
            IF (ASSOCIATED(multipoles(i)%charge_typ)) THEN
               ! make a copy of the array and enlarge the array type by 1
               isize = SIZE(multipoles(i)%charge_typ)+1
            ELSE
               isize = 1
            END IF
            CALL reallocate_charge_type(multipoles(i)%charge_typ, 1, isize)
            SELECT CASE (label)
            CASE ("CHARGE")
               CPASSERT(PRESENT(charges))
               CPASSERT(ASSOCIATED(charges))
               ALLOCATE (multipoles(i)%charge_typ(isize)%charge(1))
               ALLOCATE (multipoles(i)%charge_typ(isize)%pos(3, 1))

               multipoles(i)%charge_typ(isize)%charge(1) = echarge
               multipoles(i)%charge_typ(isize)%pos(1:3, 1) = 0.0_dp
               charges(i) = charges(i)+echarge
            CASE ("DIPOLE")
               dx = 1.0E-4_dp
               CPASSERT(PRESENT(dipoles))
               CPASSERT(ASSOCIATED(dipoles))
               ALLOCATE (multipoles(i)%charge_typ(isize)%charge(2))
               ALLOCATE (multipoles(i)%charge_typ(isize)%pos(3, 2))
               CALL random_numbers(rvec, random_stream)
               rvec = rvec/(2.0_dp*SQRT(DOT_PRODUCT(rvec, rvec)))*dx
               multipoles(i)%charge_typ(isize)%charge(1) = echarge
               multipoles(i)%charge_typ(isize)%pos(1:3, 1) = rvec
               multipoles(i)%charge_typ(isize)%charge(2) = -echarge
               multipoles(i)%charge_typ(isize)%pos(1:3, 2) = -rvec

               dipoles(:, i) = dipoles(:, i)+2.0_dp*echarge*rvec
            CASE ("QUADRUPOLE")
               dx = 1.0E-2_dp
               CPASSERT(PRESENT(quadrupoles))
               CPASSERT(ASSOCIATED(quadrupoles))
               ALLOCATE (multipoles(i)%charge_typ(isize)%charge(4))
               ALLOCATE (multipoles(i)%charge_typ(isize)%pos(3, 4))
               CALL random_numbers(rvec1, random_stream)
               CALL random_numbers(rvec2, random_stream)
               rvec1 = rvec1/SQRT(DOT_PRODUCT(rvec1, rvec1))
               rvec2 = rvec2-DOT_PRODUCT(rvec2, rvec1)*rvec1
               rvec2 = rvec2/SQRT(DOT_PRODUCT(rvec2, rvec2))
               !
               rvec1 = rvec1/2.0_dp*dx
               rvec2 = rvec2/2.0_dp*dx
               !       + (4)  ^      - (1)
               !              |rvec2
               !              |
               !              0------> rvec1
               !
               !
               !       - (3)         + (2)
               multipoles(i)%charge_typ(isize)%charge(1) = -echarge
               multipoles(i)%charge_typ(isize)%pos(1:3, 1) = rvec1+rvec2
               multipoles(i)%charge_typ(isize)%charge(2) = echarge
               multipoles(i)%charge_typ(isize)%pos(1:3, 2) = rvec1-rvec2
               multipoles(i)%charge_typ(isize)%charge(3) = -echarge
               multipoles(i)%charge_typ(isize)%pos(1:3, 3) = -rvec1-rvec2
               multipoles(i)%charge_typ(isize)%charge(4) = echarge
               multipoles(i)%charge_typ(isize)%pos(1:3, 4) = -rvec1+rvec2

               DO k = 1, 4
                  r2 = DOT_PRODUCT(multipoles(i)%charge_typ(isize)%pos(:, k), multipoles(i)%charge_typ(isize)%pos(:, k))
                  DO l = 1, 3
                     DO m = 1, 3
                        quadrupoles(m, l, i) = quadrupoles(m, l, i)+3.0_dp*0.5_dp*multipoles(i)%charge_typ(isize)%charge(k)* &
                                               multipoles(i)%charge_typ(isize)%pos(l, k)* &
                                               multipoles(i)%charge_typ(isize)%pos(m, k)
                        IF (m == l) quadrupoles(m, l, i) = quadrupoles(m, l, i)-0.5_dp*multipoles(i)%charge_typ(isize)%charge(k)*r2
                     END DO
                  END DO
               END DO

            END SELECT
         END DO
      END SUBROUTINE create_multi_type

! **************************************************************************************************
!> \brief  release multi_type for multipoles
!> \param multipoles ...
!> \date   05.2008
!> \author Teodoro Laino [tlaino] - University of Zurich - 05.2008
! **************************************************************************************************
      SUBROUTINE release_multi_type(multipoles)
      TYPE(multi_charge_type), DIMENSION(:), POINTER     :: multipoles

      CHARACTER(len=*), PARAMETER :: routineN = 'release_multi_type', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, j

         IF (ASSOCIATED(multipoles)) THEN
            DO i = 1, SIZE(multipoles)
               DO j = 1, SIZE(multipoles(i)%charge_typ)
                  DEALLOCATE (multipoles(i)%charge_typ(j)%charge)
                  DEALLOCATE (multipoles(i)%charge_typ(j)%pos)
               END DO
               DEALLOCATE (multipoles(i)%charge_typ)
            END DO
         END IF
      END SUBROUTINE release_multi_type

! **************************************************************************************************
!> \brief  reallocates multi_type for multipoles
!> \param charge_typ ...
!> \param istart ...
!> \param iend ...
!> \date   05.2008
!> \author Teodoro Laino [tlaino] - University of Zurich - 05.2008
! **************************************************************************************************
      SUBROUTINE reallocate_charge_type(charge_typ, istart, iend)
      TYPE(charge_mono_type), DIMENSION(:), POINTER      :: charge_typ
      INTEGER, INTENT(IN)                                :: istart, iend

      CHARACTER(len=*), PARAMETER :: routineN = 'reallocate_charge_type', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, isize, j, jsize, jsize1, jsize2
      TYPE(charge_mono_type), DIMENSION(:), POINTER      :: charge_typ_bk

         IF (ASSOCIATED(charge_typ)) THEN
            isize = SIZE(charge_typ)
            ALLOCATE (charge_typ_bk(1:isize))
            DO j = 1, isize
               jsize = SIZE(charge_typ(j)%charge)
               ALLOCATE (charge_typ_bk(j)%charge(jsize))
               jsize1 = SIZE(charge_typ(j)%pos, 1)
               jsize2 = SIZE(charge_typ(j)%pos, 2)
               ALLOCATE (charge_typ_bk(j)%pos(jsize1, jsize2))
               charge_typ_bk(j)%pos = charge_typ(j)%pos
               charge_typ_bk(j)%charge = charge_typ(j)%charge
            END DO
            DO j = 1, SIZE(charge_typ)
               DEALLOCATE (charge_typ(j)%charge)
               DEALLOCATE (charge_typ(j)%pos)
            END DO
            DEALLOCATE (charge_typ)
            ! Reallocate
            ALLOCATE (charge_typ_bk(istart:iend))
            DO i = istart, isize
               jsize = SIZE(charge_typ_bk(j)%charge)
               ALLOCATE (charge_typ(j)%charge(jsize))
               jsize1 = SIZE(charge_typ_bk(j)%pos, 1)
               jsize2 = SIZE(charge_typ_bk(j)%pos, 2)
               ALLOCATE (charge_typ(j)%pos(jsize1, jsize2))
               charge_typ(j)%pos = charge_typ_bk(j)%pos
               charge_typ(j)%charge = charge_typ_bk(j)%charge
            END DO
            DO j = 1, SIZE(charge_typ_bk)
               DEALLOCATE (charge_typ_bk(j)%charge)
               DEALLOCATE (charge_typ_bk(j)%pos)
            END DO
            DEALLOCATE (charge_typ_bk)
         ELSE
            ALLOCATE (charge_typ(istart:iend))
         END IF

      END SUBROUTINE reallocate_charge_type

   END SUBROUTINE debug_ewald_multipoles

! **************************************************************************************************
!> \brief  Routine to debug potential, field and electric field gradients
!> \param ewald_env ...
!> \param ewald_pw ...
!> \param nonbond_env ...
!> \param cell ...
!> \param particle_set ...
!> \param local_particles ...
!> \param radii ...
!> \param charges ...
!> \param dipoles ...
!> \param quadrupoles ...
!> \param task ...
!> \param iw ...
!> \param atomic_kind_set ...
!> \param mm_section ...
!> \date   05.2008
!> \author Teodoro Laino [tlaino] - University of Zurich - 05.2008
! **************************************************************************************************
   SUBROUTINE debug_ewald_multipoles_fields(ewald_env, ewald_pw, nonbond_env, cell, &
                                            particle_set, local_particles, radii, charges, dipoles, quadrupoles, task, iw, &
                                            atomic_kind_set, mm_section)
      TYPE(ewald_environment_type), POINTER              :: ewald_env
      TYPE(ewald_pw_type), POINTER                       :: ewald_pw
      TYPE(fist_nonbond_env_type), POINTER               :: nonbond_env
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), POINTER                       :: particle_set(:)
      TYPE(distribution_1d_type), POINTER                :: local_particles
      REAL(KIND=dp), DIMENSION(:), OPTIONAL, POINTER     :: radii, charges
      REAL(KIND=dp), DIMENSION(:, :), OPTIONAL, POINTER  :: dipoles
      REAL(KIND=dp), DIMENSION(:, :, :), OPTIONAL, &
         POINTER                                         :: quadrupoles
      LOGICAL, DIMENSION(3), INTENT(IN)                  :: task
      INTEGER, INTENT(IN)                                :: iw
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind_set(:)
      TYPE(section_vals_type), POINTER                   :: mm_section

      CHARACTER(len=*), PARAMETER :: routineN = 'debug_ewald_multipoles_fields', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, iparticle_kind, j, k, &
                                                            nparticle_local, nparticles
      REAL(KIND=dp) :: coord(3), dq, e_neut, e_self, efield1n(3), efield2n(3, 3), ene(2), &
         energy_glob, energy_local, enev(3, 2), o_tot_ene, pot, pv_glob(3, 3), pv_local(3, 3), &
         tot_ene
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: efield1, efield2, forces_glob, &
                                                            forces_local
      REAL(KIND=dp), DIMENSION(:), POINTER               :: efield0, lcharges
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(particle_type), DIMENSION(:), POINTER         :: core_particle_set, shell_particle_set

      NULLIFY (lcharges, shell_particle_set, core_particle_set)
      NULLIFY (logger)
      logger => cp_get_default_logger()

      nparticles = SIZE(particle_set)
      nparticle_local = 0
      DO iparticle_kind = 1, SIZE(local_particles%n_el)
         nparticle_local = nparticle_local+local_particles%n_el(iparticle_kind)
      END DO
      ALLOCATE (lcharges(nparticles))
      ALLOCATE (forces_glob(3, nparticles))
      ALLOCATE (forces_local(3, nparticle_local))
      ALLOCATE (efield0(nparticles))
      ALLOCATE (efield1(3, nparticles))
      ALLOCATE (efield2(9, nparticles))
      forces_glob = 0.0_dp
      forces_local = 0.0_dp
      efield0 = 0.0_dp
      efield1 = 0.0_dp
      efield2 = 0.0_dp
      pv_local = 0.0_dp
      pv_glob = 0.0_dp
      energy_glob = 0.0_dp
      energy_local = 0.0_dp
      e_neut = 0.0_dp
      e_self = 0.0_dp
      CALL ewald_multipole_evaluate(ewald_env, ewald_pw, nonbond_env, cell, particle_set, &
                                    local_particles, energy_local, energy_glob, e_neut, e_self, task, .FALSE., .TRUE., .TRUE., &
                                    .TRUE., radii, charges, dipoles, quadrupoles, forces_local, forces_glob, pv_local, pv_glob, &
                                    efield0, efield1, efield2, iw, do_debug=.FALSE.)
      o_tot_ene = energy_local+energy_glob+e_neut+e_self
      WRITE (iw, *) "TOTAL ENERGY :: ========>", o_tot_ene
      ! Debug Potential
      dq = 0.001_dp
      tot_ene = 0.0_dp
      DO i = 1, nparticles
         DO k = 1, 2
            lcharges = charges
            lcharges(i) = charges(i)+(-1.0_dp)**k*dq
            forces_glob = 0.0_dp
            forces_local = 0.0_dp
            pv_local = 0.0_dp
            pv_glob = 0.0_dp
            energy_glob = 0.0_dp
            energy_local = 0.0_dp
            e_neut = 0.0_dp
            e_self = 0.0_dp
            CALL ewald_multipole_evaluate(ewald_env, ewald_pw, nonbond_env, cell, particle_set, &
                                          local_particles, energy_local, energy_glob, e_neut, e_self, &
                                          task, .FALSE., .FALSE., .FALSE., .FALSE., radii, &
                                          lcharges, dipoles, quadrupoles, iw=iw, do_debug=.FALSE.)
            ene(k) = energy_local+energy_glob+e_neut+e_self
         END DO
         pot = (ene(2)-ene(1))/(2.0_dp*dq)
         WRITE (iw, '(A,I8,3(A,F15.9))') "POTENTIAL FOR ATOM: ", i, " NUMERICAL: ", pot, " ANALYTICAL: ", efield0(i), &
            " ERROR: ", pot-efield0(i)
         tot_ene = tot_ene+0.5_dp*efield0(i)*charges(i)
      END DO
      WRITE (iw, *) "ENERGIES: ", o_tot_ene, tot_ene, o_tot_ene-tot_ene
      WRITE (iw, '(/,/,/)')
      ! Debug Field
      dq = 0.001_dp
      DO i = 1, nparticles
         coord = particle_set(i)%r
         DO j = 1, 3
            DO k = 1, 2
               particle_set(i)%r(j) = coord(j)+(-1.0_dp)**k*dq

               ! Rebuild neighbor lists
               CALL list_control(atomic_kind_set, particle_set, local_particles, &
                                 cell, nonbond_env, logger%para_env, mm_section, &
                                 shell_particle_set, core_particle_set)

               forces_glob = 0.0_dp
               forces_local = 0.0_dp
               pv_local = 0.0_dp
               pv_glob = 0.0_dp
               energy_glob = 0.0_dp
               energy_local = 0.0_dp
               e_neut = 0.0_dp
               e_self = 0.0_dp
               efield0 = 0.0_dp
               CALL ewald_multipole_evaluate(ewald_env, ewald_pw, nonbond_env, cell, particle_set, &
                                             local_particles, energy_local, energy_glob, e_neut, e_self, &
                                             task, .FALSE., .TRUE., .TRUE., .TRUE., radii, &
                                             charges, dipoles, quadrupoles, forces_local, forces_glob, &
                                             pv_local, pv_glob, efield0, iw=iw, do_debug=.FALSE.)
               ene(k) = efield0(i)
               particle_set(i)%r(j) = coord(j)
            END DO
            efield1n(j) = -(ene(2)-ene(1))/(2.0_dp*dq)
         END DO
         WRITE (iw, '(/,A,I8)') "FIELD FOR ATOM: ", i
         WRITE (iw, '(A,3F15.9)') " NUMERICAL: ", efield1n, " ANALYTICAL: ", efield1(:, i), &
            " ERROR: ", efield1n-efield1(:, i)
         IF (task(2)) THEN
            tot_ene = tot_ene-0.5_dp*DOT_PRODUCT(efield1(:, i), dipoles(:, i))
         END IF
      END DO
      WRITE (iw, *) "ENERGIES: ", o_tot_ene, tot_ene, o_tot_ene-tot_ene

      ! Debug Field Gradient
      dq = 0.0001_dp
      DO i = 1, nparticles
         coord = particle_set(i)%r
         DO j = 1, 3
            DO k = 1, 2
               particle_set(i)%r(j) = coord(j)+(-1.0_dp)**k*dq

               ! Rebuild neighbor lists
               CALL list_control(atomic_kind_set, particle_set, local_particles, &
                                 cell, nonbond_env, logger%para_env, mm_section, &
                                 shell_particle_set, core_particle_set)

               forces_glob = 0.0_dp
               forces_local = 0.0_dp
               pv_local = 0.0_dp
               pv_glob = 0.0_dp
               energy_glob = 0.0_dp
               energy_local = 0.0_dp
               e_neut = 0.0_dp
               e_self = 0.0_dp
               efield1 = 0.0_dp
               CALL ewald_multipole_evaluate(ewald_env, ewald_pw, nonbond_env, cell, particle_set, &
                                             local_particles, energy_local, energy_glob, e_neut, e_self, &
                                             task, .FALSE., .TRUE., .TRUE., .TRUE., radii, &
                                             charges, dipoles, quadrupoles, forces_local, forces_glob, &
                                             pv_local, pv_glob, efield1=efield1, iw=iw, do_debug=.FALSE.)
               enev(:, k) = efield1(:, i)
               particle_set(i)%r(j) = coord(j)
            END DO
            efield2n(:, j) = (enev(:, 2)-enev(:, 1))/(2.0_dp*dq)
         END DO
         WRITE (iw, '(/,A,I8)') "FIELD GRADIENT FOR ATOM: ", i
         WRITE (iw, '(A,9F15.9)') " NUMERICAL:  ", efield2n, &
            " ANALYTICAL: ", efield2(:, i), &
            " ERROR:      ", RESHAPE(efield2n, (/9/))-efield2(:, i)
      END DO
   END SUBROUTINE debug_ewald_multipoles_fields

! **************************************************************************************************
!> \brief  Routine to debug potential, field and electric field gradients
!> \param ewald_env ...
!> \param ewald_pw ...
!> \param nonbond_env ...
!> \param cell ...
!> \param particle_set ...
!> \param local_particles ...
!> \param radii ...
!> \param charges ...
!> \param dipoles ...
!> \param quadrupoles ...
!> \param task ...
!> \param iw ...
!> \date   05.2008
!> \author Teodoro Laino [tlaino] - University of Zurich - 05.2008
! **************************************************************************************************
   SUBROUTINE debug_ewald_multipoles_fields2(ewald_env, ewald_pw, nonbond_env, cell, &
                                             particle_set, local_particles, radii, charges, dipoles, quadrupoles, task, iw)
      TYPE(ewald_environment_type), POINTER              :: ewald_env
      TYPE(ewald_pw_type), POINTER                       :: ewald_pw
      TYPE(fist_nonbond_env_type), POINTER               :: nonbond_env
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), POINTER                       :: particle_set(:)
      TYPE(distribution_1d_type), POINTER                :: local_particles
      REAL(KIND=dp), DIMENSION(:), OPTIONAL, POINTER     :: radii, charges
      REAL(KIND=dp), DIMENSION(:, :), OPTIONAL, POINTER  :: dipoles
      REAL(KIND=dp), DIMENSION(:, :, :), OPTIONAL, &
         POINTER                                         :: quadrupoles
      LOGICAL, DIMENSION(3), INTENT(IN)                  :: task
      INTEGER, INTENT(IN)                                :: iw

      CHARACTER(len=*), PARAMETER :: routineN = 'debug_ewald_multipoles_fields2', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: i, ind, iparticle_kind, j, k, &
                                                            nparticle_local, nparticles
      REAL(KIND=dp)                                      :: e_neut, e_self, energy_glob, &
                                                            energy_local, o_tot_ene, prod, &
                                                            pv_glob(3, 3), pv_local(3, 3), tot_ene
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: efield1, efield2, forces_glob, &
                                                            forces_local
      REAL(KIND=dp), DIMENSION(:), POINTER               :: efield0
      TYPE(cp_logger_type), POINTER                      :: logger

      NULLIFY (logger)
      logger => cp_get_default_logger()

      nparticles = SIZE(particle_set)
      nparticle_local = 0
      DO iparticle_kind = 1, SIZE(local_particles%n_el)
         nparticle_local = nparticle_local+local_particles%n_el(iparticle_kind)
      END DO
      ALLOCATE (forces_glob(3, nparticles))
      ALLOCATE (forces_local(3, nparticle_local))
      ALLOCATE (efield0(nparticles))
      ALLOCATE (efield1(3, nparticles))
      ALLOCATE (efield2(9, nparticles))
      forces_glob = 0.0_dp
      forces_local = 0.0_dp
      efield0 = 0.0_dp
      efield1 = 0.0_dp
      efield2 = 0.0_dp
      pv_local = 0.0_dp
      pv_glob = 0.0_dp
      energy_glob = 0.0_dp
      energy_local = 0.0_dp
      e_neut = 0.0_dp
      e_self = 0.0_dp
      CALL ewald_multipole_evaluate(ewald_env, ewald_pw, nonbond_env, cell, particle_set, &
                                    local_particles, energy_local, energy_glob, e_neut, e_self, task, .FALSE., .TRUE., .TRUE., &
                                    .TRUE., radii, charges, dipoles, quadrupoles, forces_local, forces_glob, pv_local, pv_glob, &
                                    efield0, efield1, efield2, iw, do_debug=.FALSE.)
      o_tot_ene = energy_local+energy_glob+e_neut+e_self
      WRITE (iw, *) "TOTAL ENERGY :: ========>", o_tot_ene

      ! Debug Potential
      tot_ene = 0.0_dp
      IF (task(1)) THEN
         DO i = 1, nparticles
            tot_ene = tot_ene+0.5_dp*efield0(i)*charges(i)
         END DO
         WRITE (iw, *) "ENERGIES: ", o_tot_ene, tot_ene, o_tot_ene-tot_ene
         WRITE (iw, '(/,/,/)')
      END IF

      ! Debug Field
      IF (task(2)) THEN
         DO i = 1, nparticles
            tot_ene = tot_ene-0.5_dp*DOT_PRODUCT(efield1(:, i), dipoles(:, i))
         END DO
         WRITE (iw, *) "ENERGIES: ", o_tot_ene, tot_ene, o_tot_ene-tot_ene
         WRITE (iw, '(/,/,/)')
      END IF

      ! Debug Field Gradient
      IF (task(3)) THEN
         DO i = 1, nparticles
            ind = 0
            prod = 0.0_dp
            DO j = 1, 3
               DO k = 1, 3
                  ind = ind+1
                  prod = prod+efield2(ind, i)*quadrupoles(j, k, i)
               END DO
            END DO
            tot_ene = tot_ene-0.5_dp*(1.0_dp/3.0_dp)*prod
         END DO
         WRITE (iw, *) "ENERGIES: ", o_tot_ene, tot_ene, o_tot_ene-tot_ene
         WRITE (iw, '(/,/,/)')
      END IF

   END SUBROUTINE debug_ewald_multipoles_fields2

END MODULE ewalds_multipole
