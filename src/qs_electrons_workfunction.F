!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Work function and Fermi level readback from the converged SCF, using the planar
!>        average of the Hartree potential over a user-defined region.
!> \author Nanchen DONGFANG
! **************************************************************************************************
MODULE qs_electrons_workfunction
   USE IEEE_ARITHMETIC,                 ONLY: IEEE_IS_FINITE
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type
   USE cp_fm_types,                     ONLY: cp_fm_type
   USE force_env_types,                 ONLY: force_env_get,&
                                              force_env_type
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: kpoint_env_type,&
                                              kpoint_type
   USE message_passing,                 ONLY: mp_para_env_type
   USE physcon,                         ONLY: bohr
   USE pw_types,                        ONLY: pw_r3d_rs_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_mo_methods,                   ONLY: calculate_subspace_eigenvalues
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: qs_electrons_compute_workfunction

   REAL(KIND=dp), PARAMETER, PRIVATE :: HARTREE_TO_EV = 27.211386245988_dp

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_electrons_workfunction'

CONTAINS

! **************************************************************************************************
!> \brief Average the Hartree potential over a slab perpendicular to a lattice axis.
!> \param v_hartree_rspace Hartree potential on the real-space grid
!> \param para_env parallel environment of the grid
!> \param axis axis (1=x, 2=y, 3=z) along which the potential is planar-averaged
!> \param minimum_A lower bound of the averaging slab, Angstrom
!> \param maximum_A upper bound of the averaging slab, Angstrom
!> \param v_avg_hartree the averaged potential, Hartree
! **************************************************************************************************
   SUBROUTINE planar_hartree_average(v_hartree_rspace, para_env, axis, minimum_A, maximum_A, &
                                     v_avg_hartree)
      TYPE(pw_r3d_rs_type), POINTER                      :: v_hartree_rspace
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: axis
      REAL(KIND=dp), INTENT(IN)                          :: minimum_A, maximum_A
      REAL(KIND=dp), INTENT(OUT)                         :: v_avg_hartree

      INTEGER                                            :: grid_index, ix, iy, iz
      INTEGER, DIMENSION(2, 3)                           :: bounds_global, bounds_local
      REAL(KIND=dp) :: cell_hmat(3, 3), coordinate_bohr, cross_product(3), face_norm, height_bohr, &
         local_count, local_sum, maximum_bohr, minimum_bohr, udvol

      IF (.NOT. ASSOCIATED(para_env)) THEN
         CPABORT("Planar-potential average could not access para_env.")
      END IF
      IF (.NOT. ASSOCIATED(v_hartree_rspace)) THEN
         CPABORT("Planar-potential average could not access the Hartree potential grid.")
      END IF
      IF (.NOT. ASSOCIATED(v_hartree_rspace%array)) THEN
         CPABORT("Planar-potential average found an unallocated Hartree potential grid.")
      END IF
      IF (.NOT. ASSOCIATED(v_hartree_rspace%pw_grid)) THEN
         CPABORT("Planar-potential average found a Hartree potential grid without metadata.")
      END IF
      IF (v_hartree_rspace%pw_grid%dvol <= 0.0_dp) THEN
         CPABORT("Planar-potential average found an invalid Hartree potential grid volume.")
      END IF

      minimum_bohr = minimum_A*bohr
      maximum_bohr = maximum_A*bohr
      IF (axis < 1 .OR. axis > 3) THEN
         CPABORT("Planar-potential average received an invalid averaging axis.")
      END IF
      IF (maximum_bohr <= minimum_bohr) THEN
         CPABORT("Planar-potential average requires potential_average_max_A > min_A.")
      END IF

      cell_hmat(:, 1) = v_hartree_rspace%pw_grid%dh(:, 1)* &
                        REAL(v_hartree_rspace%pw_grid%npts(1), KIND=dp)
      cell_hmat(:, 2) = v_hartree_rspace%pw_grid%dh(:, 2)* &
                        REAL(v_hartree_rspace%pw_grid%npts(2), KIND=dp)
      cell_hmat(:, 3) = v_hartree_rspace%pw_grid%dh(:, 3)* &
                        REAL(v_hartree_rspace%pw_grid%npts(3), KIND=dp)
      SELECT CASE (axis)
      CASE (1)
         cross_product(1) = cell_hmat(2, 2)*cell_hmat(3, 3) - cell_hmat(3, 2)*cell_hmat(2, 3)
         cross_product(2) = cell_hmat(3, 2)*cell_hmat(1, 3) - cell_hmat(1, 2)*cell_hmat(3, 3)
         cross_product(3) = cell_hmat(1, 2)*cell_hmat(2, 3) - cell_hmat(2, 2)*cell_hmat(1, 3)
      CASE (2)
         cross_product(1) = cell_hmat(2, 3)*cell_hmat(3, 1) - cell_hmat(3, 3)*cell_hmat(2, 1)
         cross_product(2) = cell_hmat(3, 3)*cell_hmat(1, 1) - cell_hmat(1, 3)*cell_hmat(3, 1)
         cross_product(3) = cell_hmat(1, 3)*cell_hmat(2, 1) - cell_hmat(2, 3)*cell_hmat(1, 1)
      CASE (3)
         cross_product(1) = cell_hmat(2, 1)*cell_hmat(3, 2) - cell_hmat(3, 1)*cell_hmat(2, 2)
         cross_product(2) = cell_hmat(3, 1)*cell_hmat(1, 2) - cell_hmat(1, 1)*cell_hmat(3, 2)
         cross_product(3) = cell_hmat(1, 1)*cell_hmat(2, 2) - cell_hmat(2, 1)*cell_hmat(1, 2)
      END SELECT
      face_norm = SQRT(SUM(cross_product**2))
      IF (.NOT. IEEE_IS_FINITE(face_norm) .OR. face_norm <= 0.0_dp) THEN
         CPABORT("Planar-potential average found a singular cell face.")
      END IF
      height_bohr = v_hartree_rspace%pw_grid%vol/face_norm
      IF (.NOT. IEEE_IS_FINITE(height_bohr) .OR. height_bohr <= 0.0_dp) THEN
         CPABORT("Planar-potential average found an invalid cell height.")
      END IF
      IF (minimum_bohr < -1.0e-10_dp .OR. maximum_bohr > height_bohr + 1.0e-10_dp) THEN
         CPABORT("CP2K potential_average_region lies outside the cell")
      END IF

      udvol = 1.0_dp/v_hartree_rspace%pw_grid%dvol
      local_sum = 0.0_dp
      local_count = 0.0_dp

      bounds_global = v_hartree_rspace%pw_grid%bounds
      bounds_local = v_hartree_rspace%pw_grid%bounds_local

      DO iz = bounds_local(1, 3), bounds_local(2, 3)
         DO iy = bounds_local(1, 2), bounds_local(2, 2)
            DO ix = bounds_local(1, 1), bounds_local(2, 1)
               SELECT CASE (axis)
               CASE (1)
                  grid_index = ix - bounds_global(1, 1)
               CASE (2)
                  grid_index = iy - bounds_global(1, 2)
               CASE (3)
                  grid_index = iz - bounds_global(1, 3)
               END SELECT
               coordinate_bohr = REAL(grid_index, KIND=dp)/ &
                                 REAL(v_hartree_rspace%pw_grid%npts(axis), KIND=dp)*height_bohr
               IF (coordinate_bohr >= minimum_bohr .AND. coordinate_bohr <= maximum_bohr) THEN
                  local_sum = local_sum + v_hartree_rspace%array(ix, iy, iz)*udvol
                  local_count = local_count + 1.0_dp
               END IF
            END DO
         END DO
      END DO

      CALL para_env%sum(local_sum)
      CALL para_env%sum(local_count)
      IF (local_count <= 0.0_dp) THEN
         CPABORT("Planar-potential average found no grid points in the requested region.")
      END IF

      v_avg_hartree = local_sum/local_count
   END SUBROUTINE planar_hartree_average

! **************************************************************************************************
!> \brief Compute the work function (vacuum Hartree level minus Fermi level, in eV)
!>      from the planar average of the Hartree potential of this force_env's own
!>      converged SCF and its Fermi level.
!> \param force_env force_env whose qs_env holds the converged SCF result
!> \param axis axis (1=x, 2=y, 3=z) along which the potential is planar-averaged
!> \param minimum_A lower bound of the averaging slab, Angstrom
!> \param maximum_A upper bound of the averaging slab, Angstrom
!> \param fermi_level_eV returned Fermi level, eV
!> \param workfunction_eV returned work function, eV
! **************************************************************************************************
   SUBROUTINE qs_electrons_compute_workfunction(force_env, axis, minimum_A, maximum_A, &
                                                fermi_level_eV, workfunction_eV)
      TYPE(force_env_type), POINTER                      :: force_env
      INTEGER, INTENT(IN)                                :: axis
      REAL(KIND=dp), INTENT(IN)                          :: minimum_A, maximum_A
      REAL(KIND=dp), INTENT(OUT)                         :: fermi_level_eV, workfunction_eV

      CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_electrons_compute_workfunction'

      INTEGER                                            :: handle, homo, ik, ispin, nfound, nmo, &
                                                            nspin
      LOGICAL                                            :: do_kpoints
      REAL(KIND=dp)                                      :: fermi_hartree, v_avg_hartree
      REAL(KIND=dp), DIMENSION(:), POINTER               :: mo_eigenvalues
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(pw_r3d_rs_type), POINTER                      :: v_hartree_rspace
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CALL timeset(routineN, handle)

      IF (axis < 1 .OR. axis > 3) THEN
         CPABORT("workfunction axis must be x, y, or z.")
      END IF
      IF (maximum_A <= minimum_A) THEN
         CPABORT("workfunction region requires maximum_A > minimum_A.")
      END IF

      NULLIFY (qs_env)
      CALL force_env_get(force_env, qs_env=qs_env)
      IF (.NOT. ASSOCIATED(qs_env)) THEN
         CPABORT("workfunction update requires a QS force environment.")
      END IF

      NULLIFY (energy, mos, dft_control, matrix_ks, mo_coeff, mo_eigenvalues, para_env, &
               v_hartree_rspace, kpoints, kp)
      do_kpoints = .FALSE.
      CALL get_qs_env(qs_env, energy=energy, mos=mos, dft_control=dft_control, &
                      para_env=para_env, v_hartree_rspace=v_hartree_rspace, &
                      do_kpoints=do_kpoints, kpoints=kpoints)
      IF (.NOT. ASSOCIATED(energy)) THEN
         CPABORT("workfunction update could not access QS energy data.")
      END IF

      ! matrix_ks is deliberately NOT in the call above.  ks_env stores it in a
      ! kpoint_transitional_type, and get_ks_env resolves that through
      ! get_1d_pointer, which does
      !
      !     IF (ASSOCIATED(t%ptr_2d)) CPABORT("Method not implemented for k-points")
      !
      ! (kpoint_transitional.F).  Under k-points the 2-D pointer is the one that
      ! is set, so merely ASKING for matrix_ks aborts -- before any branch below
      ! can decide not to use it.  Only the OT fallback needs it, and that lives
      ! on the non-k-point side, so it is fetched there.
      IF (.NOT. do_kpoints) THEN
         CALL get_qs_env(qs_env, matrix_ks=matrix_ks)
      END IF

      IF (do_kpoints) THEN
         ! Tier 0: k-point run (including GAMMA-only, which is how OT reaches
         ! ADDED_MOS + SMEARING: qs_scf_initialization refuses OT with added
         ! MOs unless complex k points and OT%ENERGIES are both active).
         ! Under k-points NONE of the tiers below work.  qs_env%mos is only a
         ! template -- kpoint_env_create sizes the k-point MO sets from it and
         ! then restores its nmo, its coefficients are never filled, and
         ! nothing writes mu into it.  energy%efermi is likewise never
         ! assigned on the k-point path.  So tier 1 reads zero, tier 2 reads
         ! an untouched 0.0, and tier 3 would diagonalize the KS matrix in a
         ! subspace of uninitialized orbitals and return a meaningless number
         ! WITHOUT aborting -- the silent failure that makes this branch
         ! necessary rather than merely convenient.
         !
         ! The real chemical potential lives on the k-point MO sets, where
         ! kpoint_set_mo_occupation stores it as mos(ispin)%mu.  That value is
         ! already global in k space: the routine sums all eigenvalues over
         ! para_env_inter_kp and solves for a single mu against the weighted
         ! k-point sum, so every local k-point of every k-point group carries
         ! the same number.  Reading the first one visible here is therefore
         ! exact, not an approximation -- no re-weighting by wkp is needed.
         IF (.NOT. ASSOCIATED(kpoints)) THEN
            CPABORT("Fermi level: k-point run without a k-point environment.")
         END IF
         fermi_hartree = -HUGE(0.0_dp)
         nfound = 0
         IF (ASSOCIATED(kpoints%kp_env)) THEN
            DO ik = 1, SIZE(kpoints%kp_env)
               kp => kpoints%kp_env(ik)%kpoint_env
               IF (.NOT. ASSOCIATED(kp)) CYCLE
               IF (.NOT. ASSOCIATED(kp%mos)) CYCLE
               DO ispin = 1, SIZE(kp%mos)
                  IF (.NOT. IEEE_IS_FINITE(kp%mos(ispin)%mu)) CYCLE
                  ! Open shell with a fixed magnetic moment gives the two spin
                  ! channels separate mu.  Take the HIGHER one, matching the OT
                  ! fallback below: it is the level an electron leaves first,
                  ! whereas the mean sits between the channels and biases e_F
                  ! whenever they differ.
                  fermi_hartree = MAX(fermi_hartree, kp%mos(ispin)%mu)
                  nfound = nfound + 1
               END DO
            END DO
         END IF
         ! kp_env is allocated only on the ranks of a k-point group that owns
         ! at least one local k-point, so the result has to be spread over the
         ! whole force_env before anyone uses it.
         CALL para_env%sum(nfound)
         CALL para_env%max(fermi_hartree)
         IF (nfound <= 0) THEN
            CPABORT("Fermi level unavailable from the k-point MO sets.")
         END IF
      ELSE
         ! Tier 1: standard diagonalization sets energy%efermi directly.
         fermi_hartree = energy%efermi
         IF (fermi_hartree == 0.0_dp) THEN
            ! Tier 2: fall back to the average MO chemical potential (e.g. smearing).
            IF (.NOT. ASSOCIATED(mos)) THEN
               CPABORT("Fermi level unset; no MO data available.")
            END IF
            IF (SIZE(mos) <= 0) CPABORT("MO array is empty.")
            nspin = SIZE(mos)
            fermi_hartree = 0.0_dp
            DO ispin = 1, nspin
               IF (.NOT. IEEE_IS_FINITE(mos(ispin)%mu)) THEN
                  CPABORT("MO chemical potential not finite.")
               END IF
               fermi_hartree = fermi_hartree + mos(ispin)%mu
            END DO
            fermi_hartree = fermi_hartree/REAL(nspin, KIND=dp)
         END IF

         IF (ABS(fermi_hartree) <= 10.0_dp*EPSILON(fermi_hartree)) THEN
            ! Tier 3: OT does not populate mo_set%mu or energy%efermi. Compute the
            ! occupied-subspace Ritz eigenvalues from the converged KS matrix and use
            ! the HOMO as the zero-temperature electronic chemical potential.
            ! Open shell: the chemical potential is the HIGHER of the two spin
            ! HOMOs (the level an electron leaves first), not their mean, which
            ! sits between them and biases e_F whenever the channels differ.
            IF (.NOT. ASSOCIATED(mos) .OR. .NOT. ASSOCIATED(matrix_ks)) THEN
               CPABORT("OT fallback needs MO and KS-matrix data.")
            END IF
            IF (SIZE(matrix_ks) < SIZE(mos)) THEN
               CPABORT("OT fallback found too few KS matrices.")
            END IF
            IF (ASSOCIATED(dft_control)) THEN
               IF (dft_control%do_admm) THEN
                  CPABORT("OT fallback does not support ADMM.")
               END IF
            END IF
            nspin = SIZE(mos)
            fermi_hartree = -HUGE(0.0_dp)
            DO ispin = 1, nspin
               NULLIFY (mo_coeff, mo_eigenvalues)
               homo = 0
               nmo = 0
               CALL get_mo_set(mos(ispin), mo_coeff=mo_coeff, &
                               eigenvalues=mo_eigenvalues, homo=homo, nmo=nmo)
               IF (.NOT. ASSOCIATED(mo_coeff) .OR. .NOT. ASSOCIATED(mo_eigenvalues)) THEN
                  CPABORT("OT fallback could not access orbitals.")
               END IF
               IF (homo <= 0 .OR. homo > nmo .OR. homo > SIZE(mo_eigenvalues)) THEN
                  CPABORT("OT fallback found an invalid HOMO index.")
               END IF
               CALL calculate_subspace_eigenvalues(orbitals=mo_coeff, &
                                                   ks_matrix=matrix_ks(ispin)%matrix, &
                                                   evals_arg=mo_eigenvalues, &
                                                   do_rotation=.FALSE.)
               IF (.NOT. IEEE_IS_FINITE(mo_eigenvalues(homo))) THEN
                  CPABORT("OT HOMO energy not finite.")
               END IF
               fermi_hartree = MAX(fermi_hartree, mo_eigenvalues(homo))
            END DO
         END IF

         IF (ABS(fermi_hartree) <= 10.0_dp*EPSILON(fermi_hartree)) THEN
            CPABORT("Fermi level unset after OT fallback.")
         END IF
      END IF

      IF (.NOT. IEEE_IS_FINITE(fermi_hartree)) THEN
         CPABORT("Fermi level not finite.")
      END IF
      fermi_level_eV = fermi_hartree*HARTREE_TO_EV

      ! The average is taken from this force_env's own Hartree potential, so that the two
      ! endpoints of a charge-mixing run each contribute their own vacuum level.
      CALL planar_hartree_average(v_hartree_rspace, para_env, axis, minimum_A, maximum_A, &
                                  v_avg_hartree)
      workfunction_eV = v_avg_hartree*HARTREE_TO_EV - fermi_level_eV

      CALL timestop(handle)

   END SUBROUTINE qs_electrons_compute_workfunction

END MODULE qs_electrons_workfunction
