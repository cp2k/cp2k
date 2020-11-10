!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Build up the plane wave density by collocating the primitive Gaussian
!>      functions (pgf).
!> \par History
!>      Joost VandeVondele (02.2002)
!>            1) rewrote collocate_pgf for increased accuracy and speed
!>            2) collocate_core hack for PGI compiler
!>            3) added multiple grid feature
!>            4) new way to go over the grid
!>      Joost VandeVondele (05.2002)
!>            1) prelim. introduction of the real space grid type
!>      JGH [30.08.02] multigrid arrays independent from potential
!>      JGH [17.07.03] distributed real space code
!>      JGH [23.11.03] refactoring and new loop ordering
!>      JGH [04.12.03] OpneMP parallelization of main loops
!>      Joost VandeVondele (12.2003)
!>           1) modified to compute tau
!>      Joost removed incremental build feature
!>      Joost introduced map consistent
!>      Rewrote grid integration/collocation routines, [Joost VandeVondele,03.2007]
!> \author Matthias Krack (03.04.2001)
! **************************************************************************************************
MODULE qs_integrate_potential_single
   USE ao_util,                         ONLY: exp_radius_very_extended
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE atprop_types,                    ONLY: atprop_array_init,&
                                              atprop_type
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              pbc
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_get_block_p,&
                                              dbcsr_type
   USE external_potential_types,        ONLY: get_potential,&
                                              gth_potential_type
   USE gaussian_gridlevels,             ONLY: gaussian_gridlevel,&
                                              gridlevel_info_type
   USE grid_api,                        ONLY: integrate_pgf_product
   USE kinds,                           ONLY: dp
   USE lri_environment_types,           ONLY: lri_kind_type
   USE memory_utilities,                ONLY: reallocate
   USE message_passing,                 ONLY: mp_sum
   USE orbital_pointers,                ONLY: coset,&
                                              ncoset
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_types,                        ONLY: pw_p_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE realspace_grid_types,            ONLY: pw2rs,&
                                              realspace_grid_p_type,&
                                              realspace_grid_type,&
                                              rs_grid_release,&
                                              rs_grid_retain,&
                                              rs_grid_zero,&
                                              rs_pw_transfer
   USE rs_pw_interface,                 ONLY: potential_pw2rs
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .FALSE.

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_integrate_potential_single'

! *** Public subroutines ***
! *** Don't include this routines directly, use the interface to
! *** qs_integrate_potential

   PUBLIC :: integrate_v_rspace_one_center, &
             integrate_v_rspace_diagonal, &
             integrate_v_core_rspace, &
             integrate_ppl_rspace, &
             integrate_rho_nlcc

CONTAINS

! **************************************************************************************************
!> \brief computes the forces/virial due to the local pseudopotential
!> \param rho_rspace ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE integrate_ppl_rspace(rho_rspace, qs_env)
      TYPE(pw_p_type), INTENT(INOUT)                     :: rho_rspace
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'integrate_ppl_rspace'

      INTEGER                                            :: atom_a, handle, iatom, ikind, j, lppl, &
                                                            n, natom_of_kind, ni, npme
      INTEGER, DIMENSION(:), POINTER                     :: atom_list, cores
      LOGICAL                                            :: use_virial
      REAL(KIND=dp)                                      :: alpha, eps_rho_rspace, radius
      REAL(KIND=dp), DIMENSION(3)                        :: force_a, force_b, ra
      REAL(KIND=dp), DIMENSION(3, 3)                     :: my_virial_a, my_virial_b
      REAL(KIND=dp), DIMENSION(:), POINTER               :: cexp_ppl
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: hab, pab
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gth_potential_type), POINTER                  :: gth_potential
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(realspace_grid_type), POINTER                 :: rs_v
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      NULLIFY (pw_env, cores)

      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env)
      CALL pw_env_get(pw_env=pw_env, auxbas_rs_grid=rs_v)
      CALL rs_grid_retain(rs_v)

      CALL rs_pw_transfer(rs_v, rho_rspace%pw, pw2rs)

      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      cell=cell, &
                      dft_control=dft_control, &
                      particle_set=particle_set, &
                      pw_env=pw_env, &
                      force=force, virial=virial)

      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)

      eps_rho_rspace = dft_control%qs_control%eps_rho_rspace

      DO ikind = 1, SIZE(atomic_kind_set)

         CALL get_atomic_kind(atomic_kind_set(ikind), natom=natom_of_kind, atom_list=atom_list)
         CALL get_qs_kind(qs_kind_set(ikind), gth_potential=gth_potential)

         IF (.NOT. ASSOCIATED(gth_potential)) CYCLE
         CALL get_potential(potential=gth_potential, alpha_ppl=alpha, nexp_ppl=lppl, cexp_ppl=cexp_ppl)

         IF (lppl <= 0) CYCLE

         ni = ncoset(2*lppl - 2)
         ALLOCATE (hab(ni, 1), pab(ni, 1))
         pab = 0._dp

         CALL reallocate(cores, 1, natom_of_kind)
         npme = 0
         cores = 0

         ! prepare core function
         DO j = 1, lppl
            SELECT CASE (j)
            CASE (1)
               pab(1, 1) = cexp_ppl(1)
            CASE (2)
               n = coset(2, 0, 0)
               pab(n, 1) = cexp_ppl(2)
               n = coset(0, 2, 0)
               pab(n, 1) = cexp_ppl(2)
               n = coset(0, 0, 2)
               pab(n, 1) = cexp_ppl(2)
            CASE (3)
               n = coset(4, 0, 0)
               pab(n, 1) = cexp_ppl(3)
               n = coset(0, 4, 0)
               pab(n, 1) = cexp_ppl(3)
               n = coset(0, 0, 4)
               pab(n, 1) = cexp_ppl(3)
               n = coset(2, 2, 0)
               pab(n, 1) = 2._dp*cexp_ppl(3)
               n = coset(2, 0, 2)
               pab(n, 1) = 2._dp*cexp_ppl(3)
               n = coset(0, 2, 2)
               pab(n, 1) = 2._dp*cexp_ppl(3)
            CASE (4)
               n = coset(6, 0, 0)
               pab(n, 1) = cexp_ppl(4)
               n = coset(0, 6, 0)
               pab(n, 1) = cexp_ppl(4)
               n = coset(0, 0, 6)
               pab(n, 1) = cexp_ppl(4)
               n = coset(4, 2, 0)
               pab(n, 1) = 3._dp*cexp_ppl(4)
               n = coset(4, 0, 2)
               pab(n, 1) = 3._dp*cexp_ppl(4)
               n = coset(2, 4, 0)
               pab(n, 1) = 3._dp*cexp_ppl(4)
               n = coset(2, 0, 4)
               pab(n, 1) = 3._dp*cexp_ppl(4)
               n = coset(0, 4, 2)
               pab(n, 1) = 3._dp*cexp_ppl(4)
               n = coset(0, 2, 4)
               pab(n, 1) = 3._dp*cexp_ppl(4)
               n = coset(2, 2, 2)
               pab(n, 1) = 6._dp*cexp_ppl(4)
            CASE DEFAULT
               CPABORT("")
            END SELECT
         END DO

         DO iatom = 1, natom_of_kind
            atom_a = atom_list(iatom)
            ra(:) = pbc(particle_set(atom_a)%r, cell)
            IF (rs_v%desc%parallel .AND. .NOT. rs_v%desc%distributed) THEN
               ! replicated realspace grid, split the atoms up between procs
               IF (MODULO(iatom, rs_v%desc%group_size) == rs_v%desc%my_pos) THEN
                  npme = npme + 1
                  cores(npme) = iatom
               END IF
            ELSE
               npme = npme + 1
               cores(npme) = iatom
            END IF
         END DO

         DO j = 1, npme

            iatom = cores(j)
            atom_a = atom_list(iatom)
            ra(:) = pbc(particle_set(atom_a)%r, cell)
            hab(:, 1) = 0.0_dp
            force_a(:) = 0.0_dp
            force_b(:) = 0.0_dp
            IF (use_virial) THEN
               my_virial_a = 0.0_dp
               my_virial_b = 0.0_dp
            END IF
            ni = 2*lppl - 2

            radius = exp_radius_very_extended(la_min=0, la_max=ni, lb_min=0, lb_max=0, &
                                              ra=ra, rb=ra, rp=ra, &
                                              zetp=alpha, eps=eps_rho_rspace, &
                                              pab=pab, o1=0, o2=0, &  ! without map_consistent
                                              prefactor=1.0_dp, cutoff=1.0_dp)

            CALL integrate_pgf_product(ni, alpha, 0, &
                                       0, 0.0_dp, 0, ra, (/0.0_dp, 0.0_dp, 0.0_dp/), &
                                       rs_v, cell, pw_env%cube_info(1), hab, pab=pab, o1=0, o2=0, &
                                       radius=radius, &
                                       calculate_forces=.TRUE., force_a=force_a, &
                                       force_b=force_b, use_virial=use_virial, my_virial_a=my_virial_a, &
                                       my_virial_b=my_virial_b, use_subpatch=.TRUE., subpatch_pattern=0)

            force(ikind)%gth_ppl(:, iatom) = &
               force(ikind)%gth_ppl(:, iatom) + force_a(:)*rho_rspace%pw%pw_grid%dvol

            IF (use_virial) THEN
               virial%pv_ppl = virial%pv_ppl + my_virial_a*rho_rspace%pw%pw_grid%dvol
               virial%pv_virial = virial%pv_virial + my_virial_a*rho_rspace%pw%pw_grid%dvol
               CPABORT("Virial not debuged for CORE_PPL")
            END IF
         END DO

         DEALLOCATE (hab, pab)

      END DO

      CALL rs_grid_release(rs_v)

      DEALLOCATE (cores)

      CALL timestop(handle)

   END SUBROUTINE integrate_ppl_rspace

! **************************************************************************************************
!> \brief computes the forces/virial due to the nlcc pseudopotential
!> \param rho_rspace ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE integrate_rho_nlcc(rho_rspace, qs_env)
      TYPE(pw_p_type), INTENT(INOUT)                     :: rho_rspace
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'integrate_rho_nlcc'

      INTEGER                                            :: atom_a, handle, iatom, iexp_nlcc, ikind, &
                                                            ithread, j, n, natom, nc, nexp_nlcc, &
                                                            ni, npme, nthread
      INTEGER, DIMENSION(:), POINTER                     :: atom_list, cores, nct_nlcc
      LOGICAL                                            :: nlcc, use_virial
      REAL(KIND=dp)                                      :: alpha, eps_rho_rspace, radius
      REAL(KIND=dp), DIMENSION(3)                        :: force_a, force_b, ra
      REAL(KIND=dp), DIMENSION(3, 3)                     :: my_virial_a, my_virial_b
      REAL(KIND=dp), DIMENSION(:), POINTER               :: alpha_nlcc
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: cval_nlcc, hab, pab
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gth_potential_type), POINTER                  :: gth_potential
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(realspace_grid_type), POINTER                 :: rs_v
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      NULLIFY (pw_env, cores)

      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env)
      CALL pw_env_get(pw_env=pw_env, auxbas_rs_grid=rs_v)
      CALL rs_grid_retain(rs_v)

      CALL rs_pw_transfer(rs_v, rho_rspace%pw, pw2rs)

      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      cell=cell, &
                      dft_control=dft_control, &
                      particle_set=particle_set, &
                      pw_env=pw_env, &
                      force=force, virial=virial)

      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)

      eps_rho_rspace = dft_control%qs_control%eps_rho_rspace

      DO ikind = 1, SIZE(atomic_kind_set)

         CALL get_atomic_kind(atomic_kind_set(ikind), natom=natom, atom_list=atom_list)
         CALL get_qs_kind(qs_kind_set(ikind), gth_potential=gth_potential)

         IF (.NOT. ASSOCIATED(gth_potential)) CYCLE
         CALL get_potential(potential=gth_potential, nlcc_present=nlcc, nexp_nlcc=nexp_nlcc, &
                            alpha_nlcc=alpha_nlcc, nct_nlcc=nct_nlcc, cval_nlcc=cval_nlcc)

         IF (.NOT. nlcc) CYCLE

         DO iexp_nlcc = 1, nexp_nlcc

            alpha = alpha_nlcc(iexp_nlcc)
            nc = nct_nlcc(iexp_nlcc)

            ni = ncoset(2*nc - 2)

            nthread = 1
            ithread = 0

            ALLOCATE (hab(ni, 1), pab(ni, 1))
            pab = 0._dp

            CALL reallocate(cores, 1, natom)
            npme = 0
            cores = 0

            ! prepare core function
            DO j = 1, nc
               SELECT CASE (j)
               CASE (1)
                  pab(1, 1) = cval_nlcc(1, iexp_nlcc)
               CASE (2)
                  n = coset(2, 0, 0)
                  pab(n, 1) = cval_nlcc(2, iexp_nlcc)/alpha**2
                  n = coset(0, 2, 0)
                  pab(n, 1) = cval_nlcc(2, iexp_nlcc)/alpha**2
                  n = coset(0, 0, 2)
                  pab(n, 1) = cval_nlcc(2, iexp_nlcc)/alpha**2
               CASE (3)
                  n = coset(4, 0, 0)
                  pab(n, 1) = cval_nlcc(3, iexp_nlcc)/alpha**4
                  n = coset(0, 4, 0)
                  pab(n, 1) = cval_nlcc(3, iexp_nlcc)/alpha**4
                  n = coset(0, 0, 4)
                  pab(n, 1) = cval_nlcc(3, iexp_nlcc)/alpha**4
                  n = coset(2, 2, 0)
                  pab(n, 1) = 2._dp*cval_nlcc(3, iexp_nlcc)/alpha**4
                  n = coset(2, 0, 2)
                  pab(n, 1) = 2._dp*cval_nlcc(3, iexp_nlcc)/alpha**4
                  n = coset(0, 2, 2)
                  pab(n, 1) = 2._dp*cval_nlcc(3, iexp_nlcc)/alpha**4
               CASE (4)
                  n = coset(6, 0, 0)
                  pab(n, 1) = cval_nlcc(4, iexp_nlcc)/alpha**6
                  n = coset(0, 6, 0)
                  pab(n, 1) = cval_nlcc(4, iexp_nlcc)/alpha**6
                  n = coset(0, 0, 6)
                  pab(n, 1) = cval_nlcc(4, iexp_nlcc)/alpha**6
                  n = coset(4, 2, 0)
                  pab(n, 1) = 3._dp*cval_nlcc(4, iexp_nlcc)/alpha**6
                  n = coset(4, 0, 2)
                  pab(n, 1) = 3._dp*cval_nlcc(4, iexp_nlcc)/alpha**6
                  n = coset(2, 4, 0)
                  pab(n, 1) = 3._dp*cval_nlcc(4, iexp_nlcc)/alpha**6
                  n = coset(2, 0, 4)
                  pab(n, 1) = 3._dp*cval_nlcc(4, iexp_nlcc)/alpha**6
                  n = coset(0, 4, 2)
                  pab(n, 1) = 3._dp*cval_nlcc(4, iexp_nlcc)/alpha**6
                  n = coset(0, 2, 4)
                  pab(n, 1) = 3._dp*cval_nlcc(4, iexp_nlcc)/alpha**6
                  n = coset(2, 2, 2)
                  pab(n, 1) = 6._dp*cval_nlcc(4, iexp_nlcc)/alpha**6
               CASE DEFAULT
                  CPABORT("")
               END SELECT
            END DO
            IF (dft_control%nspins == 2) pab = pab*0.5_dp

            DO iatom = 1, natom
               atom_a = atom_list(iatom)
               ra(:) = pbc(particle_set(atom_a)%r, cell)
               IF (rs_v%desc%parallel .AND. .NOT. rs_v%desc%distributed) THEN
                  ! replicated realspace grid, split the atoms up between procs
                  IF (MODULO(iatom, rs_v%desc%group_size) == rs_v%desc%my_pos) THEN
                     npme = npme + 1
                     cores(npme) = iatom
                  END IF
               ELSE
                  npme = npme + 1
                  cores(npme) = iatom
               END IF
            END DO

            DO j = 1, npme

               iatom = cores(j)
               atom_a = atom_list(iatom)
               ra(:) = pbc(particle_set(atom_a)%r, cell)
               hab(:, 1) = 0.0_dp
               force_a(:) = 0.0_dp
               force_b(:) = 0.0_dp
               IF (use_virial) THEN
                  my_virial_a = 0.0_dp
                  my_virial_b = 0.0_dp
               END IF
               ni = 2*nc - 2

               radius = exp_radius_very_extended(la_min=0, la_max=ni, lb_min=0, lb_max=0, &
                                                 ra=ra, rb=ra, rp=ra, &
                                                 zetp=1/(2*alpha**2), eps=eps_rho_rspace, &
                                                 pab=pab, o1=0, o2=0, &  ! without map_consistent
                                                 prefactor=1.0_dp, cutoff=1.0_dp)

               CALL integrate_pgf_product(ni, 1/(2*alpha**2), 0, &
                                          0, 0.0_dp, 0, ra, (/0.0_dp, 0.0_dp, 0.0_dp/), &
                                          rs_v, cell, pw_env%cube_info(1), hab, pab=pab, o1=0, o2=0, &
                                          radius=radius, &
                                          calculate_forces=.TRUE., force_a=force_a, &
                                          force_b=force_b, use_virial=use_virial, my_virial_a=my_virial_a, &
                                          my_virial_b=my_virial_b, use_subpatch=.TRUE., subpatch_pattern=0)

               force(ikind)%gth_nlcc(:, iatom) = &
                  force(ikind)%gth_nlcc(:, iatom) + force_a(:)*rho_rspace%pw%pw_grid%dvol

               IF (use_virial) THEN
                  virial%pv_nlcc = virial%pv_nlcc + my_virial_a*rho_rspace%pw%pw_grid%dvol
                  virial%pv_virial = virial%pv_virial + my_virial_a*rho_rspace%pw%pw_grid%dvol
               END IF
            END DO

            DEALLOCATE (hab, pab)

         END DO

      END DO

      CALL rs_grid_release(rs_v)

      DEALLOCATE (cores)

      CALL timestop(handle)

   END SUBROUTINE integrate_rho_nlcc

! **************************************************************************************************
!> \brief computes the forces/virial due to the ionic cores with a potential on
!>      grid
!> \param v_rspace ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE integrate_v_core_rspace(v_rspace, qs_env)
      TYPE(pw_p_type), INTENT(INOUT)                     :: v_rspace
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'integrate_v_core_rspace'

      INTEGER                                            :: atom_a, handle, iatom, ikind, j, natom, &
                                                            natom_of_kind, npme
      INTEGER, DIMENSION(:), POINTER                     :: atom_list, cores
      LOGICAL                                            :: paw_atom, skip_fcore, use_virial
      REAL(KIND=dp)                                      :: alpha_core_charge, ccore_charge, &
                                                            eps_rho_rspace, radius
      REAL(KIND=dp), DIMENSION(3)                        :: force_a, force_b, ra
      REAL(KIND=dp), DIMENSION(3, 3)                     :: my_virial_a, my_virial_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: hab, pab
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atprop_type), POINTER                         :: atprop
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(realspace_grid_type), POINTER                 :: rs_v
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)
      NULLIFY (virial, force, atprop, dft_control)

      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control)

      !If gapw, check for gpw kinds
      skip_fcore = .FALSE.
      IF (dft_control%qs_control%gapw) THEN
         IF (.NOT. dft_control%qs_control%gapw_control%nopaw_as_gpw) skip_fcore = .TRUE.
      END IF

      IF (.NOT. skip_fcore) THEN
         NULLIFY (pw_env)
         ALLOCATE (cores(1))
         ALLOCATE (hab(1, 1))
         ALLOCATE (pab(1, 1))

         CALL get_qs_env(qs_env=qs_env, pw_env=pw_env)
         CALL pw_env_get(pw_env=pw_env, auxbas_rs_grid=rs_v)
         CALL rs_grid_retain(rs_v)

         CALL rs_pw_transfer(rs_v, v_rspace%pw, pw2rs)

         CALL get_qs_env(qs_env=qs_env, &
                         atomic_kind_set=atomic_kind_set, &
                         qs_kind_set=qs_kind_set, &
                         cell=cell, &
                         dft_control=dft_control, &
                         particle_set=particle_set, &
                         pw_env=pw_env, &
                         force=force, &
                         virial=virial, &
                         atprop=atprop)

         ! atomic energy contributions
         IF (ASSOCIATED(atprop)) THEN
            natom = SIZE(particle_set)
            CALL atprop_array_init(atprop%ateb, natom)
         END IF

         use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)

         eps_rho_rspace = dft_control%qs_control%eps_rho_rspace

         DO ikind = 1, SIZE(atomic_kind_set)

            CALL get_atomic_kind(atomic_kind_set(ikind), natom=natom_of_kind, atom_list=atom_list)
            CALL get_qs_kind(qs_kind_set(ikind), paw_atom=paw_atom, &
                             alpha_core_charge=alpha_core_charge, &
                             ccore_charge=ccore_charge)

            IF (paw_atom) CYCLE

            pab(1, 1) = -ccore_charge
            IF (alpha_core_charge == 0.0_dp .OR. pab(1, 1) == 0.0_dp) CYCLE

            CALL reallocate(cores, 1, natom_of_kind)
            npme = 0
            cores = 0

            DO iatom = 1, natom_of_kind
               atom_a = atom_list(iatom)
               ra(:) = pbc(particle_set(atom_a)%r, cell)
               IF (rs_v%desc%parallel .AND. .NOT. rs_v%desc%distributed) THEN
                  ! replicated realspace grid, split the atoms up between procs
                  IF (MODULO(iatom, rs_v%desc%group_size) == rs_v%desc%my_pos) THEN
                     npme = npme + 1
                     cores(npme) = iatom
                  END IF
               ELSE
                  npme = npme + 1
                  cores(npme) = iatom
               END IF
            END DO

            DO j = 1, npme

               iatom = cores(j)
               atom_a = atom_list(iatom)
               ra(:) = pbc(particle_set(atom_a)%r, cell)
               hab(1, 1) = 0.0_dp
               force_a(:) = 0.0_dp
               force_b(:) = 0.0_dp
               IF (use_virial) THEN
                  my_virial_a = 0.0_dp
                  my_virial_b = 0.0_dp
               END IF

               radius = exp_radius_very_extended(la_min=0, la_max=0, lb_min=0, lb_max=0, &
                                                 ra=ra, rb=ra, rp=ra, &
                                                 zetp=alpha_core_charge, eps=eps_rho_rspace, &
                                                 pab=pab, o1=0, o2=0, &  ! without map_consistent
                                                 prefactor=1.0_dp, cutoff=1.0_dp)

               CALL integrate_pgf_product(0, alpha_core_charge, 0, &
                                          0, 0.0_dp, 0, ra, (/0.0_dp, 0.0_dp, 0.0_dp/), &
                                          rs_v, cell, pw_env%cube_info(1), hab, pab=pab, o1=0, o2=0, &
                                          radius=radius, &
                                          calculate_forces=.TRUE., force_a=force_a, &
                                          force_b=force_b, use_virial=use_virial, my_virial_a=my_virial_a, &
                                          my_virial_b=my_virial_b, use_subpatch=.TRUE., subpatch_pattern=0)

               IF (ASSOCIATED(force)) THEN
                  force(ikind)%rho_core(:, iatom) = force(ikind)%rho_core(:, iatom) + force_a(:)
               END IF
               IF (use_virial) THEN
                  virial%pv_ehartree = virial%pv_ehartree + my_virial_a
                  virial%pv_virial = virial%pv_virial + my_virial_a
               END IF
               IF (ASSOCIATED(atprop)) THEN
                  atprop%ateb(atom_a) = atprop%ateb(atom_a) + 0.5_dp*hab(1, 1)*pab(1, 1)
               END IF

            END DO

         END DO

         CALL rs_grid_release(rs_v)

         DEALLOCATE (hab, pab, cores)

      END IF

      CALL timestop(handle)

   END SUBROUTINE integrate_v_core_rspace

! **************************************************************************************************
!> \brief computes integrals of product of v_rspace times a one-center function
!>        required for LRIGPW
!> \param v_rspace ...
!> \param qs_env ...
!> \param int_res ...
!> \param calculate_forces ...
!> \param basis_type ...
!> \param atomlist ...
!> \author Dorothea Golze
! **************************************************************************************************
   SUBROUTINE integrate_v_rspace_one_center(v_rspace, qs_env, int_res, &
                                            calculate_forces, basis_type, atomlist)
      TYPE(pw_p_type), INTENT(IN)                        :: v_rspace
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(lri_kind_type), DIMENSION(:), POINTER         :: int_res
      LOGICAL, INTENT(IN)                                :: calculate_forces
      CHARACTER(len=*), INTENT(IN)                       :: basis_type
      INTEGER, DIMENSION(:), OPTIONAL                    :: atomlist

      CHARACTER(len=*), PARAMETER :: routineN = 'integrate_v_rspace_one_center'

      INTEGER :: atom_a, dir, group_size, handle, i, iatom, igrid_level, ikind, ipgf, iset, m1, &
         maxco, maxsgf_set, my_pos, na1, natom_of_kind, ncoa, nkind, nseta, offset, sgfa
      INTEGER, DIMENSION(3)                              :: lb, location, tp, ub
      INTEGER, DIMENSION(:), POINTER                     :: atom_list, la_max, la_min, npgfa, &
                                                            nsgf_seta
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa
      LOGICAL                                            :: use_virial
      LOGICAL, ALLOCATABLE, DIMENSION(:)                 :: map_it
      REAL(KIND=dp)                                      :: eps_rho_rspace, radius
      REAL(KIND=dp), DIMENSION(3)                        :: force_a, force_b, ra
      REAL(KIND=dp), DIMENSION(3, 3)                     :: my_virial_a, my_virial_b
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: hab, pab, rpgfa, sphi_a, work_f, work_i, &
                                                            zeta
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gridlevel_info_type), POINTER                 :: gridlevel_info
      TYPE(gto_basis_set_type), POINTER                  :: lri_basis_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(realspace_grid_p_type), DIMENSION(:), POINTER :: rs_v
      TYPE(realspace_grid_type), POINTER                 :: rs_grid
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      NULLIFY (atomic_kind_set, qs_kind_set, atom_list, cell, dft_control, &
               first_sgfa, gridlevel_info, hab, la_max, la_min, lri_basis_set, &
               npgfa, nsgf_seta, pab, particle_set, pw_env, rpgfa, &
               rs_grid, rs_v, virial, set_radius_a, sphi_a, work_f, &
               work_i, zeta)

      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env)

      CALL pw_env_get(pw_env, rs_grids=rs_v)
      DO i = 1, SIZE(rs_v)
         CALL rs_grid_retain(rs_v(i)%rs_grid)
         CALL rs_grid_zero(rs_v(i)%rs_grid)
      END DO

      gridlevel_info => pw_env%gridlevel_info

      CALL potential_pw2rs(rs_v, v_rspace, pw_env)

      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      cell=cell, &
                      dft_control=dft_control, &
                      nkind=nkind, &
                      particle_set=particle_set, &
                      pw_env=pw_env, &
                      virial=virial)

      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)

      eps_rho_rspace = dft_control%qs_control%eps_rho_rspace

      offset = 0
      my_pos = v_rspace%pw%pw_grid%para%my_pos
      group_size = v_rspace%pw%pw_grid%para%group_size

      DO ikind = 1, nkind

         CALL get_atomic_kind(atomic_kind_set(ikind), natom=natom_of_kind, atom_list=atom_list)
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=lri_basis_set, basis_type=basis_type)
         CALL get_gto_basis_set(gto_basis_set=lri_basis_set, &
                                first_sgf=first_sgfa, &
                                lmax=la_max, &
                                lmin=la_min, &
                                maxco=maxco, &
                                maxsgf_set=maxsgf_set, &
                                npgf=npgfa, &
                                nset=nseta, &
                                nsgf_set=nsgf_seta, &
                                pgf_radius=rpgfa, &
                                set_radius=set_radius_a, &
                                sphi=sphi_a, &
                                zet=zeta)

         ALLOCATE (hab(maxco, 1), pab(maxco, 1))
         hab = 0._dp
         pab(:, 1) = 0._dp

         DO iatom = 1, natom_of_kind

            atom_a = atom_list(iatom)
            IF (PRESENT(atomlist)) THEN
               IF (atomlist(atom_a) == 0) CYCLE
            END IF
            ra(:) = pbc(particle_set(atom_a)%r, cell)
            force_a(:) = 0._dp
            force_b(:) = 0._dp
            my_virial_a(:, :) = 0._dp
            my_virial_b(:, :) = 0._dp

            m1 = MAXVAL(npgfa(1:nseta))
            ALLOCATE (map_it(m1))

            DO iset = 1, nseta
               !
               map_it = .FALSE.
               DO ipgf = 1, npgfa(iset)
                  igrid_level = gaussian_gridlevel(gridlevel_info, zeta(ipgf, iset))
                  rs_grid => rs_v(igrid_level)%rs_grid
                  IF (.NOT. ALL(rs_grid%desc%perd == 1)) THEN
                     DO dir = 1, 3
                        ! bounds of local grid (i.e. removing the 'wings'), if periodic
                        tp(dir) = FLOOR(DOT_PRODUCT(cell%h_inv(dir, :), ra)*rs_grid%desc%npts(dir))
                        tp(dir) = MODULO(tp(dir), rs_grid%desc%npts(dir))
                        IF (rs_grid%desc%perd(dir) .NE. 1) THEN
                           lb(dir) = rs_grid%lb_local(dir) + rs_grid%desc%border
                           ub(dir) = rs_grid%ub_local(dir) - rs_grid%desc%border
                        ELSE
                           lb(dir) = rs_grid%lb_local(dir)
                           ub(dir) = rs_grid%ub_local(dir)
                        END IF
                        ! distributed grid, only map if it is local to the grid
                        location(dir) = tp(dir) + rs_grid%desc%lb(dir)
                     END DO
                     IF (lb(1) <= location(1) .AND. location(1) <= ub(1) .AND. &
                         lb(2) <= location(2) .AND. location(2) <= ub(2) .AND. &
                         lb(3) <= location(3) .AND. location(3) <= ub(3)) THEN
                        map_it(ipgf) = .TRUE.
                     END IF
                  ELSE
                     ! not distributed, just a round-robin distribution over the full set of CPUs
                     IF (MODULO(offset, group_size) == my_pos) map_it(ipgf) = .TRUE.
                  END IF
               END DO
               offset = offset + 1
               !
               IF (ANY(map_it(1:npgfa(iset)))) THEN
                  sgfa = first_sgfa(1, iset)
                  ncoa = npgfa(iset)*ncoset(la_max(iset))
                  hab(:, 1) = 0._dp
                  ALLOCATE (work_i(nsgf_seta(iset), 1))
                  work_i = 0.0_dp

                  ! get fit coefficients for forces
                  IF (calculate_forces) THEN
                     m1 = sgfa + nsgf_seta(iset) - 1
                     ALLOCATE (work_f(nsgf_seta(iset), 1))
                     work_f(1:nsgf_seta(iset), 1) = int_res(ikind)%acoef(iatom, sgfa:m1)
                     CALL dgemm("N", "N", ncoa, 1, nsgf_seta(iset), 1.0_dp, sphi_a(1, sgfa), &
                                SIZE(sphi_a, 1), work_f(1, 1), SIZE(work_f, 1), 0.0_dp, pab(1, 1), &
                                SIZE(pab, 1))
                     DEALLOCATE (work_f)
                  END IF

                  DO ipgf = 1, npgfa(iset)
                     na1 = (ipgf - 1)*ncoset(la_max(iset))
                     igrid_level = gaussian_gridlevel(gridlevel_info, zeta(ipgf, iset))
                     rs_grid => rs_v(igrid_level)%rs_grid

                     radius = exp_radius_very_extended(la_min=la_min(iset), la_max=la_max(iset), &
                                                       lb_min=0, lb_max=0, ra=ra, rb=ra, rp=ra, &
                                                       zetp=zeta(ipgf, iset), eps=eps_rho_rspace, &
                                                       prefactor=1.0_dp, cutoff=1.0_dp)

                     IF (map_it(ipgf)) THEN
                        IF (.NOT. calculate_forces) THEN
                           CALL integrate_pgf_product(la_max=la_max(iset), &
                                                      zeta=zeta(ipgf, iset), la_min=la_min(iset), &
                                                      lb_max=0, zetb=0.0_dp, lb_min=0, &
                                                      ra=ra, rab=(/0.0_dp, 0.0_dp, 0.0_dp/), &
                                                      rsgrid=rs_grid, cell=cell, &
                                                      cube_info=pw_env%cube_info(igrid_level), &
                                                      hab=hab, o1=na1, o2=0, radius=radius, &
                                                      calculate_forces=calculate_forces)
                        ELSE
                           CALL integrate_pgf_product(la_max=la_max(iset), &
                                                      zeta=zeta(ipgf, iset), la_min=la_min(iset), &
                                                      lb_max=0, zetb=0.0_dp, lb_min=0, &
                                                      ra=ra, rab=(/0.0_dp, 0.0_dp, 0.0_dp/), &
                                                      rsgrid=rs_grid, cell=cell, &
                                                      cube_info=pw_env%cube_info(igrid_level), &
                                                      hab=hab, pab=pab, o1=na1, o2=0, radius=radius, &
                                                      calculate_forces=calculate_forces, &
                                                      force_a=force_a, force_b=force_b, &
                                                      use_virial=use_virial, &
                                                      my_virial_a=my_virial_a, my_virial_b=my_virial_b)
                        END IF
                     END IF
                  END DO
                  ! contract hab
                  CALL dgemm("T", "N", nsgf_seta(iset), 1, ncoa, 1.0_dp, sphi_a(1, sgfa), &
                             SIZE(sphi_a, 1), hab(1, 1), SIZE(hab, 1), 0.0_dp, work_i(1, 1), SIZE(work_i, 1))

                  int_res(ikind)%v_int(iatom, sgfa:sgfa - 1 + nsgf_seta(iset)) = &
                     int_res(ikind)%v_int(iatom, sgfa:sgfa - 1 + nsgf_seta(iset)) + work_i(1:nsgf_seta(iset), 1)
                  DEALLOCATE (work_i)
               END IF
            END DO
            !
            IF (calculate_forces) THEN
               int_res(ikind)%v_dfdr(iatom, :) = int_res(ikind)%v_dfdr(iatom, :) + force_a(:)
               IF (use_virial) THEN
                  virial%pv_lrigpw = virial%pv_lrigpw + my_virial_a
                  virial%pv_virial = virial%pv_virial + my_virial_a
               END IF
            END IF

            DEALLOCATE (map_it)

         END DO

         DEALLOCATE (hab, pab)
      END DO

      DO i = 1, SIZE(rs_v)
         CALL rs_grid_release(rs_v(i)%rs_grid)
      END DO

      CALL timestop(handle)

   END SUBROUTINE integrate_v_rspace_one_center

! **************************************************************************************************
!> \brief computes integrals of product of v_rspace times the diagonal block basis functions
!>        required for LRIGPW with exact 1c terms
!> \param v_rspace ...
!> \param ksmat ...
!> \param pmat ...
!> \param qs_env ...
!> \param calculate_forces ...
!> \param basis_type ...
!> \author JGH
! **************************************************************************************************
   SUBROUTINE integrate_v_rspace_diagonal(v_rspace, ksmat, pmat, qs_env, calculate_forces, basis_type)
      TYPE(pw_p_type), INTENT(IN)                        :: v_rspace
      TYPE(dbcsr_type), INTENT(INOUT)                    :: ksmat, pmat
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: calculate_forces
      CHARACTER(len=*), INTENT(IN)                       :: basis_type

      CHARACTER(len=*), PARAMETER :: routineN = 'integrate_v_rspace_diagonal'

      INTEGER :: atom_a, dir, group_size, handle, i, iatom, igrid_level, ikind, ipgf, iset, jpgf, &
         jset, m1, maxco, maxsgf_set, my_pos, na1, na2, natom_of_kind, nb1, nb2, ncoa, ncob, &
         nkind, nseta, nsgfa, offset, sgfa, sgfb
      INTEGER, DIMENSION(3)                              :: lb, location, tp, ub
      INTEGER, DIMENSION(:), POINTER                     :: atom_list, la_max, la_min, npgfa, &
                                                            nsgf_seta
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa
      LOGICAL                                            :: found, use_virial
      LOGICAL, ALLOCATABLE, DIMENSION(:, :)              :: map_it2
      REAL(KIND=dp)                                      :: eps_rho_rspace, radius, zetp
      REAL(KIND=dp), DIMENSION(3)                        :: force_a, force_b, ra
      REAL(KIND=dp), DIMENSION(3, 3)                     :: my_virial_a, my_virial_b
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_a
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: h_block, hab, hmat, p_block, pab, pblk, &
                                                            rpgfa, sphi_a, work, zeta
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gridlevel_info_type), POINTER                 :: gridlevel_info
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(realspace_grid_p_type), DIMENSION(:), POINTER :: rs_v
      TYPE(realspace_grid_type), POINTER                 :: rs_grid
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env)
      CALL pw_env_get(pw_env, rs_grids=rs_v)
      DO i = 1, SIZE(rs_v)
         CALL rs_grid_retain(rs_v(i)%rs_grid)
      END DO
      CALL potential_pw2rs(rs_v, v_rspace, pw_env)

      gridlevel_info => pw_env%gridlevel_info

      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      cell=cell, &
                      dft_control=dft_control, &
                      nkind=nkind, &
                      particle_set=particle_set, &
                      force=force, &
                      virial=virial, &
                      para_env=para_env)

      eps_rho_rspace = dft_control%qs_control%eps_rho_rspace
      use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)

      offset = 0
      my_pos = v_rspace%pw%pw_grid%para%my_pos
      group_size = v_rspace%pw%pw_grid%para%group_size

      DO ikind = 1, nkind

         CALL get_atomic_kind(atomic_kind_set(ikind), natom=natom_of_kind, atom_list=atom_list)
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set, basis_type=basis_type)
         CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                                lmax=la_max, lmin=la_min, maxco=maxco, maxsgf_set=maxsgf_set, &
                                npgf=npgfa, nset=nseta, nsgf_set=nsgf_seta, nsgf=nsgfa, &
                                first_sgf=first_sgfa, pgf_radius=rpgfa, set_radius=set_radius_a, &
                                sphi=sphi_a, zet=zeta)

         ALLOCATE (hab(maxco, maxco), work(maxco, maxsgf_set), hmat(nsgfa, nsgfa))
         IF (calculate_forces) ALLOCATE (pab(maxco, maxco), pblk(nsgfa, nsgfa))

         DO iatom = 1, natom_of_kind
            atom_a = atom_list(iatom)
            ra(:) = pbc(particle_set(atom_a)%r, cell)
            hmat = 0.0_dp
            IF (calculate_forces) THEN
               CALL dbcsr_get_block_p(matrix=pmat, row=atom_a, col=atom_a, BLOCK=p_block, found=found)
               IF (found) THEN
                  pblk(1:nsgfa, 1:nsgfa) = p_block(1:nsgfa, 1:nsgfa)
               ELSE
                  pblk = 0.0_dp
               END IF
               CALL mp_sum(pblk, para_env%group)
               force_a(:) = 0._dp
               force_b(:) = 0._dp
               IF (use_virial) THEN
                  my_virial_a = 0.0_dp
                  my_virial_b = 0.0_dp
               END IF
            END IF
            m1 = MAXVAL(npgfa(1:nseta))
            ALLOCATE (map_it2(m1, m1))
            DO iset = 1, nseta
               sgfa = first_sgfa(1, iset)
               ncoa = npgfa(iset)*ncoset(la_max(iset))
               DO jset = 1, nseta
                  sgfb = first_sgfa(1, jset)
                  ncob = npgfa(jset)*ncoset(la_max(jset))
                  !
                  map_it2 = .FALSE.
                  DO ipgf = 1, npgfa(iset)
                     DO jpgf = 1, npgfa(jset)
                        zetp = zeta(ipgf, iset) + zeta(jpgf, jset)
                        igrid_level = gaussian_gridlevel(gridlevel_info, zetp)
                        rs_grid => rs_v(igrid_level)%rs_grid
                        IF (.NOT. ALL(rs_grid%desc%perd == 1)) THEN
                           DO dir = 1, 3
                              ! bounds of local grid (i.e. removing the 'wings'), if periodic
                              tp(dir) = FLOOR(DOT_PRODUCT(cell%h_inv(dir, :), ra)*rs_grid%desc%npts(dir))
                              tp(dir) = MODULO(tp(dir), rs_grid%desc%npts(dir))
                              IF (rs_grid%desc%perd(dir) .NE. 1) THEN
                                 lb(dir) = rs_grid%lb_local(dir) + rs_grid%desc%border
                                 ub(dir) = rs_grid%ub_local(dir) - rs_grid%desc%border
                              ELSE
                                 lb(dir) = rs_grid%lb_local(dir)
                                 ub(dir) = rs_grid%ub_local(dir)
                              END IF
                              ! distributed grid, only map if it is local to the grid
                              location(dir) = tp(dir) + rs_grid%desc%lb(dir)
                           END DO
                           IF (lb(1) <= location(1) .AND. location(1) <= ub(1) .AND. &
                               lb(2) <= location(2) .AND. location(2) <= ub(2) .AND. &
                               lb(3) <= location(3) .AND. location(3) <= ub(3)) THEN
                              map_it2(ipgf, jpgf) = .TRUE.
                           END IF
                        ELSE
                           ! not distributed, just a round-robin distribution over the full set of CPUs
                           IF (MODULO(offset, group_size) == my_pos) map_it2(ipgf, jpgf) = .TRUE.
                        END IF
                     END DO
                  END DO
                  offset = offset + 1
                  !
                  IF (ANY(map_it2(1:npgfa(iset), 1:npgfa(jset)))) THEN
                     hab(:, :) = 0._dp
                     IF (calculate_forces) THEN
                        CALL dgemm("N", "N", ncoa, nsgf_seta(jset), nsgf_seta(iset), &
                                   1.0_dp, sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                                   pblk(sgfa, sgfb), SIZE(pblk, 1), &
                                   0.0_dp, work(1, 1), SIZE(work, 1))
                        CALL dgemm("N", "T", ncoa, ncob, nsgf_seta(jset), &
                                   1.0_dp, work(1, 1), SIZE(work, 1), &
                                   sphi_a(1, sgfb), SIZE(sphi_a, 1), &
                                   0.0_dp, pab(1, 1), SIZE(pab, 1))
                     END IF

                     DO ipgf = 1, npgfa(iset)
                        na1 = (ipgf - 1)*ncoset(la_max(iset))
                        na2 = ipgf*ncoset(la_max(iset))
                        DO jpgf = 1, npgfa(jset)
                           nb1 = (jpgf - 1)*ncoset(la_max(jset))
                           nb2 = jpgf*ncoset(la_max(jset))
                           zetp = zeta(ipgf, iset) + zeta(jpgf, jset)
                           igrid_level = gaussian_gridlevel(gridlevel_info, zetp)
                           rs_grid => rs_v(igrid_level)%rs_grid

                           radius = exp_radius_very_extended(la_min=la_min(iset), la_max=la_max(iset), &
                                                             lb_min=la_min(jset), lb_max=la_max(jset), &
                                                             ra=ra, rb=ra, rp=ra, &
                                                             zetp=zetp, eps=eps_rho_rspace, &
                                                             prefactor=1.0_dp, cutoff=1.0_dp)

                           IF (map_it2(ipgf, jpgf)) THEN
                              IF (calculate_forces) THEN
                                 CALL integrate_pgf_product( &
                                    la_max(iset), zeta(ipgf, iset), la_min(iset), &
                                    la_max(jset), zeta(jpgf, jset), la_min(jset), &
                                    ra=ra, rab=(/0.0_dp, 0.0_dp, 0.0_dp/), &
                                    rsgrid=rs_v(igrid_level)%rs_grid, cell=cell, &
                                    cube_info=pw_env%cube_info(igrid_level), &
                                    hab=hab, pab=pab, o1=na1, o2=nb1, &
                                    radius=radius, &
                                    calculate_forces=.TRUE., &
                                    force_a=force_a, force_b=force_b, &
                                    use_virial=use_virial, my_virial_a=my_virial_a, my_virial_b=my_virial_b)
                              ELSE
                                 CALL integrate_pgf_product( &
                                    la_max(iset), zeta(ipgf, iset), la_min(iset), &
                                    la_max(jset), zeta(jpgf, jset), la_min(jset), &
                                    ra=ra, rab=(/0.0_dp, 0.0_dp, 0.0_dp/), &
                                    rsgrid=rs_v(igrid_level)%rs_grid, cell=cell, &
                                    cube_info=pw_env%cube_info(igrid_level), &
                                    hab=hab, o1=na1, o2=nb1, &
                                    radius=radius, &
                                    calculate_forces=.FALSE.)
                              END IF
                           END IF
                        END DO
                     END DO
                     ! contract hab
                     CALL dgemm("N", "N", ncoa, nsgf_seta(jset), ncob, &
                                1.0_dp, hab(1, 1), SIZE(hab, 1), &
                                sphi_a(1, sgfb), SIZE(sphi_a, 1), &
                                0.0_dp, work(1, 1), SIZE(work, 1))
                     CALL dgemm("T", "N", nsgf_seta(iset), nsgf_seta(jset), ncoa, &
                                1.0_dp, sphi_a(1, sgfa), SIZE(sphi_a, 1), &
                                work(1, 1), SIZE(work, 1), &
                                1.0_dp, hmat(sgfa, sgfb), SIZE(hmat, 1))
                  END IF
               END DO
            END DO
            DEALLOCATE (map_it2)
            ! update KS matrix
            CALL mp_sum(hmat, para_env%group)
            CALL dbcsr_get_block_p(matrix=ksmat, row=atom_a, col=atom_a, BLOCK=h_block, found=found)
            IF (found) THEN
               h_block(1:nsgfa, 1:nsgfa) = h_block(1:nsgfa, 1:nsgfa) + hmat(1:nsgfa, 1:nsgfa)
            END IF
            IF (calculate_forces) THEN
               force(ikind)%rho_elec(:, iatom) = force(ikind)%rho_elec(:, iatom) + 2.0_dp*force_a(:)
               IF (use_virial) THEN
                  IF (use_virial .AND. calculate_forces) THEN
                     virial%pv_lrigpw = virial%pv_lrigpw + 2.0_dp*my_virial_a
                     virial%pv_virial = virial%pv_virial + 2.0_dp*my_virial_a
                  END IF
               END IF
            END IF
         END DO
         DEALLOCATE (hab, work, hmat)
         IF (calculate_forces) DEALLOCATE (pab, pblk)
      END DO

      DO i = 1, SIZE(rs_v)
         CALL rs_grid_release(rs_v(i)%rs_grid)
      END DO

      CALL timestop(handle)

   END SUBROUTINE integrate_v_rspace_diagonal

END MODULE qs_integrate_potential_single
