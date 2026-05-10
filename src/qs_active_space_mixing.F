!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Dense density mixing for active-space embedding.
! **************************************************************************************************
MODULE qs_active_space_mixing
   USE cp_fm_types,                     ONLY: cp_fm_get_element,&
                                              cp_fm_set_element,&
                                              cp_fm_type
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: invert_matrix
   USE qs_active_space_types,           ONLY: active_space_type
   USE qs_density_mixing_types,         ONLY: &
        broyden_mixing_nr, direct_mixing_nr, gspace_mixing_nr, mixing_storage_create, &
        mixing_storage_type, modified_broyden_mixing_nr, multisecant_mixing_nr, no_mixing_nr, &
        pulay_mixing_nr
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: active_space_mixing_label
   PUBLIC :: initialize_active_space_mixing
   PUBLIC :: update_active_density

CONTAINS

! **************************************************************************************************
!> \brief Initialize the density mixer used in the self-consistent active-space embedding loop.
!> \param active_space_env active space environment
!> \param as_input ACTIVE_SPACE input section
! **************************************************************************************************
   SUBROUTINE initialize_active_space_mixing(active_space_env, as_input)
      TYPE(active_space_type), POINTER                   :: active_space_env
      TYPE(section_vals_type), POINTER                   :: as_input

      LOGICAL                                            :: do_mixing, legacy_alpha_explicit, &
                                                            mixing_alpha_explicit, mixing_explicit
      REAL(KIND=dp)                                      :: legacy_alpha
      TYPE(section_vals_type), POINTER                   :: mixing_section

      NULLIFY (mixing_section)

      CALL section_vals_val_get(as_input, "ALPHA", r_val=legacy_alpha, &
                                explicit=legacy_alpha_explicit)
      IF (legacy_alpha < 0.0_dp .OR. legacy_alpha > 1.0_dp) THEN
         CPABORT("Specify an active-space damping factor between 0 and 1.")
      END IF

      mixing_section => section_vals_get_subs_vals(as_input, "MIXING")
      CALL section_vals_get(mixing_section, explicit=mixing_explicit)
      CALL section_vals_val_get(mixing_section, "_SECTION_PARAMETERS_", l_val=do_mixing)

      active_space_env%as_mixing_dim = 0
      active_space_env%as_mixing_iter = 0

      IF (.NOT. do_mixing) THEN
         active_space_env%as_mixing_method = no_mixing_nr
         active_space_env%alpha = 1.0_dp
         RETURN
      END IF

      CALL section_vals_val_get(mixing_section, "METHOD", &
                                i_val=active_space_env%as_mixing_method)

      SELECT CASE (active_space_env%as_mixing_method)
      CASE (no_mixing_nr)
         active_space_env%alpha = 1.0_dp
         RETURN
      CASE (direct_mixing_nr, pulay_mixing_nr, broyden_mixing_nr, modified_broyden_mixing_nr)
         CONTINUE
      CASE (gspace_mixing_nr)
         CALL cp_abort(__LOCATION__, &
                       "ACTIVE_SPACE%MIXING%METHOD KERKER_MIXING is not supported. "// &
                       "The active-space density lives in the active MO subspace, "// &
                       "not in G-space.")
      CASE (multisecant_mixing_nr)
         CPABORT("ACTIVE_SPACE%MIXING%METHOD MULTISECANT_MIXING is not yet supported.")
      CASE DEFAULT
         CPABORT("Unknown ACTIVE_SPACE%MIXING%METHOD.")
      END SELECT

      IF (ASSOCIATED(active_space_env%as_mixing_store)) THEN
         CPABORT("Active-space mixing storage already initialized.")
      END IF
      ALLOCATE (active_space_env%as_mixing_store)
      CALL mixing_storage_create(active_space_env%as_mixing_store, mixing_section, &
                                 active_space_env%as_mixing_method, ecut=0.0_dp)

      CALL section_vals_val_get(mixing_section, "ALPHA", explicit=mixing_alpha_explicit)
      IF ((legacy_alpha_explicit .AND. (.NOT. mixing_alpha_explicit)) .OR. &
          ((.NOT. mixing_explicit) .AND. (.NOT. mixing_alpha_explicit))) THEN
         active_space_env%as_mixing_store%alpha = legacy_alpha
      END IF
      IF (active_space_env%as_mixing_store%alpha < 0.0_dp .OR. &
          active_space_env%as_mixing_store%alpha > 1.0_dp) THEN
         CPABORT("Specify an active-space mixing ALPHA between 0 and 1.")
      END IF
      IF (active_space_env%as_mixing_store%nbuffer < 1 .AND. &
          active_space_env%as_mixing_method /= direct_mixing_nr) THEN
         CPABORT("ACTIVE_SPACE%MIXING%NBUFFER has to be positive.")
      END IF

      active_space_env%alpha = active_space_env%as_mixing_store%alpha

   END SUBROUTINE initialize_active_space_mixing

! **************************************************************************************************
!> \brief Return the current active-space mixer label for iteration output.
!> \param active_space_env active space environment
!> \return short mixer label
! **************************************************************************************************
   FUNCTION active_space_mixing_label(active_space_env) RESULT(label)
      TYPE(active_space_type), POINTER                   :: active_space_env
      CHARACTER(len=15)                                  :: label

      SELECT CASE (active_space_env%as_mixing_method)
      CASE (direct_mixing_nr)
         label = "P_Mix"
      CASE (pulay_mixing_nr)
         label = "Pulay"
      CASE (broyden_mixing_nr)
         label = "Broy."
      CASE (modified_broyden_mixing_nr)
         label = "MBroy"
      CASE DEFAULT
         label = "NoMix"
      END SELECT
      IF (ASSOCIATED(active_space_env%as_mixing_store)) THEN
         IF (active_space_env%as_mixing_iter > 0 .AND. &
             LEN_TRIM(active_space_env%as_mixing_store%iter_method) > 0) THEN
            label = active_space_env%as_mixing_store%iter_method
         END IF
      END IF

   END FUNCTION active_space_mixing_label

! **************************************************************************************************
!> \brief Release dense active-space mixing history buffers.
!> \param active_space_env active space environment
! **************************************************************************************************
   SUBROUTINE release_active_mixing_history(active_space_env)
      TYPE(active_space_type), POINTER                   :: active_space_env

      IF (ASSOCIATED(active_space_env%as_mix_r_old)) THEN
         DEALLOCATE (active_space_env%as_mix_r_old)
      END IF
      IF (ASSOCIATED(active_space_env%as_mix_weight)) THEN
         DEALLOCATE (active_space_env%as_mix_weight)
      END IF
      IF (ASSOCIATED(active_space_env%as_mix_x_old)) THEN
         DEALLOCATE (active_space_env%as_mix_x_old)
      END IF
      IF (ASSOCIATED(active_space_env%as_mix_r_buffer)) THEN
         DEALLOCATE (active_space_env%as_mix_r_buffer)
      END IF
      IF (ASSOCIATED(active_space_env%as_mix_x_buffer)) THEN
         DEALLOCATE (active_space_env%as_mix_x_buffer)
      END IF
      active_space_env%as_mixing_dim = 0

   END SUBROUTINE release_active_mixing_history

! **************************************************************************************************
!> \brief Ensure dense active-space mixing history buffers are allocated.
!> \param active_space_env active space environment
!> \param ndim length of the flattened density vector
! **************************************************************************************************
   SUBROUTINE ensure_active_mixing_history(active_space_env, ndim)
      TYPE(active_space_type), POINTER                   :: active_space_env
      INTEGER, INTENT(IN)                                :: ndim

      INTEGER                                            :: nbuffer
      TYPE(mixing_storage_type), POINTER                 :: mixing_store

      IF (.NOT. ASSOCIATED(active_space_env%as_mixing_store)) RETURN
      IF (active_space_env%as_mixing_method == direct_mixing_nr) RETURN

      mixing_store => active_space_env%as_mixing_store
      nbuffer = mixing_store%nbuffer
      IF (nbuffer < 1) CPABORT("ACTIVE_SPACE%MIXING%NBUFFER has to be positive.")

      IF (active_space_env%as_mixing_dim == ndim .AND. &
          ASSOCIATED(active_space_env%as_mix_r_buffer)) RETURN

      CALL release_active_mixing_history(active_space_env)

      active_space_env%as_mixing_dim = ndim
      ALLOCATE (active_space_env%as_mix_r_old(ndim))
      ALLOCATE (active_space_env%as_mix_weight(nbuffer))
      ALLOCATE (active_space_env%as_mix_x_old(ndim))
      ALLOCATE (active_space_env%as_mix_r_buffer(nbuffer, ndim))
      ALLOCATE (active_space_env%as_mix_x_buffer(nbuffer, ndim))

      active_space_env%as_mix_r_old = 0.0_dp
      active_space_env%as_mix_weight = 1.0_dp
      active_space_env%as_mix_x_old = 0.0_dp
      active_space_env%as_mix_r_buffer = 0.0_dp
      active_space_env%as_mix_x_buffer = 0.0_dp
      mixing_store%ncall = 0

   END SUBROUTINE ensure_active_mixing_history

! **************************************************************************************************
!> \brief Apply a direct dense active-space density mix.
!> \param p_old current active-space density vector
!> \param p_solver active-space density vector returned by the external solver
!> \param alpha mixing damping
!> \param p_mixed mixed active-space density vector
! **************************************************************************************************
   SUBROUTINE active_space_direct_mix(p_old, p_solver, alpha, p_mixed)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: p_old, p_solver
      REAL(KIND=dp), INTENT(IN)                          :: alpha
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: p_mixed

      p_mixed = p_old + alpha*(p_solver - p_old)

   END SUBROUTINE active_space_direct_mix

! **************************************************************************************************
!> \brief Apply Pulay mixing to the dense active-space density vector.
!> \param active_space_env active space environment
!> \param p_old current active-space density vector
!> \param p_solver active-space density vector returned by the external solver
!> \param p_mixed mixed active-space density vector
! **************************************************************************************************
   SUBROUTINE active_space_pulay_mixing(active_space_env, p_old, p_solver, p_mixed)
      TYPE(active_space_type), POINTER                   :: active_space_env
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: p_old, p_solver
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: p_mixed

      INTEGER                                            :: i, ib, ibb, j, nb, nbuffer, ndim
      REAL(KIND=dp)                                      :: inv_err, norm_c_inv, res_norm
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: alpha_c, p_diis
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: c, c_inv
      TYPE(mixing_storage_type), POINTER                 :: mixing_store

      ndim = SIZE(p_old)
      CALL ensure_active_mixing_history(active_space_env, ndim)
      mixing_store => active_space_env%as_mixing_store
      nbuffer = mixing_store%nbuffer

      ib = MODULO(mixing_store%ncall, nbuffer) + 1
      mixing_store%ncall = mixing_store%ncall + 1
      nb = MIN(mixing_store%ncall, nbuffer)
      ibb = MODULO(mixing_store%ncall, nbuffer) + 1

      active_space_env%as_mix_x_buffer(ib, :) = p_old
      active_space_env%as_mix_r_buffer(ib, :) = p_solver - p_old
      res_norm = SQRT(DOT_PRODUCT(active_space_env%as_mix_r_buffer(ib, :), &
                                  active_space_env%as_mix_r_buffer(ib, :)))

      IF (nb == 1 .OR. res_norm < 1.E-14_dp) THEN
         CALL active_space_direct_mix(p_old, p_solver, mixing_store%alpha, p_mixed)
      ELSE
         ALLOCATE (c(nb, nb))
         ALLOCATE (c_inv(nb, nb))
         ALLOCATE (alpha_c(nb))
         ALLOCATE (p_diis(ndim))

         c(:, :) = 0.0_dp
         DO i = 1, nb
            DO j = i, nb
               c(j, i) = DOT_PRODUCT(active_space_env%as_mix_r_buffer(i, :), &
                                     active_space_env%as_mix_r_buffer(j, :))
               c(i, j) = c(j, i)
            END DO
         END DO

         CALL invert_matrix(c, c_inv, inv_err, improve=.TRUE.)
         norm_c_inv = SUM(c_inv)
         IF (ABS(norm_c_inv) < 1.E-14_dp) THEN
            CALL active_space_direct_mix(p_old, p_solver, mixing_store%alpha, p_mixed)
         ELSE
            DO i = 1, nb
               alpha_c(i) = SUM(c_inv(:, i))/norm_c_inv
            END DO

            p_diis(:) = 0.0_dp
            DO i = 1, nb
               p_diis(:) = p_diis(:) + alpha_c(i)*(active_space_env%as_mix_x_buffer(i, :) + &
                                                   mixing_store%pulay_beta*active_space_env%as_mix_r_buffer(i, :))
            END DO
            IF (mixing_store%pulay_alpha > 0.0_dp) THEN
               p_mixed = mixing_store%pulay_alpha*p_solver + &
                         (1.0_dp - mixing_store%pulay_alpha)*p_diis
            ELSE
               p_mixed = p_diis
            END IF
         END IF

         DEALLOCATE (alpha_c)
         DEALLOCATE (p_diis)
         DEALLOCATE (c)
         DEALLOCATE (c_inv)
      END IF

      active_space_env%as_mix_x_buffer(ibb, :) = p_mixed
      mixing_store%iter_method = "Pulay"

   END SUBROUTINE active_space_pulay_mixing

! **************************************************************************************************
!> \brief Apply original or modified Broyden mixing to the dense active-space density vector.
!> \param active_space_env active space environment
!> \param p_old current active-space density vector
!> \param p_solver active-space density vector returned by the external solver
!> \param p_mixed mixed active-space density vector
!> \param modified use dynamic residual weights of modified Broyden
! **************************************************************************************************
   SUBROUTINE active_space_broyden_mixing(active_space_env, p_old, p_solver, p_mixed, modified)
      TYPE(active_space_type), POINTER                   :: active_space_env
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: p_old, p_solver
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: p_mixed
      LOGICAL, INTENT(IN)                                :: modified

      INTEGER                                            :: ib, j, k, nb, nbuffer, ndim
      LOGICAL                                            :: can_update
      REAL(KIND=dp)                                      :: delta_norm, inv_err, res_norm, weight
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: c, g, p_res
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: a, b
      TYPE(mixing_storage_type), POINTER                 :: mixing_store

      ndim = SIZE(p_old)
      CALL ensure_active_mixing_history(active_space_env, ndim)
      mixing_store => active_space_env%as_mixing_store
      nbuffer = mixing_store%nbuffer

      ALLOCATE (p_res(ndim))
      p_res(:) = p_solver(:) - p_old(:)
      res_norm = SQRT(DOT_PRODUCT(p_res, p_res))

      mixing_store%ncall = mixing_store%ncall + 1
      IF (mixing_store%ncall == 1) THEN
         CALL active_space_direct_mix(p_old, p_solver, mixing_store%alpha, p_mixed)
         active_space_env%as_mix_x_old = p_old
         active_space_env%as_mix_r_old = p_res
         IF (modified) THEN
            mixing_store%iter_method = "MBroy"
         ELSE
            mixing_store%iter_method = "Broy."
         END IF
         DEALLOCATE (p_res)
         RETURN
      END IF

      nb = MIN(mixing_store%ncall - 1, nbuffer)
      ib = MODULO(mixing_store%ncall - 2, nbuffer) + 1

      active_space_env%as_mix_r_buffer(ib, :) = p_res - active_space_env%as_mix_r_old
      active_space_env%as_mix_x_buffer(ib, :) = p_old - active_space_env%as_mix_x_old
      delta_norm = SQRT(DOT_PRODUCT(active_space_env%as_mix_r_buffer(ib, :), &
                                    active_space_env%as_mix_r_buffer(ib, :)))
      can_update = res_norm > 1.E-14_dp .AND. delta_norm > 1.E-14_dp

      IF (can_update) THEN
         active_space_env%as_mix_r_buffer(ib, :) = active_space_env%as_mix_r_buffer(ib, :)/delta_norm
         active_space_env%as_mix_x_buffer(ib, :) = active_space_env%as_mix_x_buffer(ib, :)/delta_norm + &
                                                   mixing_store%alpha*active_space_env%as_mix_r_buffer(ib, :)

         IF (modified) THEN
            IF (res_norm > (mixing_store%wc/mixing_store%wmax)) THEN
               active_space_env%as_mix_weight(ib) = mixing_store%wc/res_norm
            ELSE
               active_space_env%as_mix_weight(ib) = mixing_store%wmax
            END IF
            active_space_env%as_mix_weight(ib) = MAX(1.0_dp, active_space_env%as_mix_weight(ib))
         END IF

         ALLOCATE (a(nb, nb))
         ALLOCATE (b(nb, nb))
         ALLOCATE (c(nb))
         ALLOCATE (g(nb))

         a(:, :) = 0.0_dp
         c(:) = 0.0_dp
         DO j = 1, nb
            DO k = j, nb
               a(k, j) = DOT_PRODUCT(active_space_env%as_mix_r_buffer(j, :), &
                                     active_space_env%as_mix_r_buffer(k, :))
               a(j, k) = a(k, j)
            END DO
         END DO

         DO j = 1, nb
            c(j) = DOT_PRODUCT(active_space_env%as_mix_r_buffer(j, :), p_res)
            IF (modified) THEN
               c(j) = active_space_env%as_mix_weight(j)*c(j)
               DO k = 1, nb
                  a(k, j) = active_space_env%as_mix_weight(k)* &
                            active_space_env%as_mix_weight(j)*a(k, j)
               END DO
               a(j, j) = mixing_store%broy_w0*mixing_store%broy_w0 + a(j, j)
            ELSE
               a(j, j) = mixing_store%broy_w0 + a(j, j)
            END IF
         END DO

         CALL invert_matrix(a, b, inv_err)
         g(:) = 0.0_dp
         DO j = 1, nb
            DO k = 1, nb
               g(j) = g(j) + b(k, j)*c(k)
            END DO
         END DO

         CALL active_space_direct_mix(p_old, p_solver, mixing_store%alpha, p_mixed)
         DO j = 1, nb
            weight = 1.0_dp
            IF (modified) weight = active_space_env%as_mix_weight(j)
            p_mixed = p_mixed - weight*g(j)*active_space_env%as_mix_x_buffer(j, :)
         END DO

         DEALLOCATE (a)
         DEALLOCATE (b)
         DEALLOCATE (c)
         DEALLOCATE (g)
      ELSE
         active_space_env%as_mix_r_buffer(ib, :) = 0.0_dp
         active_space_env%as_mix_x_buffer(ib, :) = 0.0_dp
         CALL active_space_direct_mix(p_old, p_solver, mixing_store%alpha, p_mixed)
      END IF

      active_space_env%as_mix_x_old = p_old
      active_space_env%as_mix_r_old = p_res
      IF (modified) THEN
         mixing_store%iter_method = "MBroy"
      ELSE
         mixing_store%iter_method = "Broy."
      END IF

      DEALLOCATE (p_res)

   END SUBROUTINE active_space_broyden_mixing

! **************************************************************************************************
!> \brief Mix the dense active-space density vector.
!> \param active_space_env active space environment
!> \param p_old current active-space density vector
!> \param p_solver active-space density vector returned by the external solver
!> \param p_mixed mixed active-space density vector
! **************************************************************************************************
   SUBROUTINE mix_active_density_vector(active_space_env, p_old, p_solver, p_mixed)
      TYPE(active_space_type), POINTER                   :: active_space_env
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: p_old, p_solver
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: p_mixed

      INTEGER                                            :: active_mix_iter
      TYPE(mixing_storage_type), POINTER                 :: mixing_store

      IF (active_space_env%as_mixing_method == no_mixing_nr .OR. &
          .NOT. ASSOCIATED(active_space_env%as_mixing_store)) THEN
         p_mixed = p_solver
         RETURN
      END IF

      mixing_store => active_space_env%as_mixing_store
      active_space_env%as_mixing_iter = active_space_env%as_mixing_iter + 1

      IF (active_space_env%as_mixing_iter <= mixing_store%nskip_mixing) THEN
         p_mixed = p_solver
         mixing_store%iter_method = "NoMix"
         RETURN
      END IF

      active_mix_iter = active_space_env%as_mixing_iter - mixing_store%nskip_mixing
      IF (active_space_env%as_mixing_method == direct_mixing_nr .OR. &
          active_mix_iter <= mixing_store%n_simple_mix) THEN
         CALL active_space_direct_mix(p_old, p_solver, mixing_store%alpha, p_mixed)
         mixing_store%iter_method = "P_Mix"
         RETURN
      END IF

      SELECT CASE (active_space_env%as_mixing_method)
      CASE (pulay_mixing_nr)
         CALL active_space_pulay_mixing(active_space_env, p_old, p_solver, p_mixed)
      CASE (broyden_mixing_nr)
         CALL active_space_broyden_mixing(active_space_env, p_old, p_solver, p_mixed, .FALSE.)
      CASE (modified_broyden_mixing_nr)
         CALL active_space_broyden_mixing(active_space_env, p_old, p_solver, p_mixed, .TRUE.)
      CASE DEFAULT
         CPABORT("Unsupported ACTIVE_SPACE%MIXING%METHOD.")
      END SELECT

   END SUBROUTINE mix_active_density_vector

! **************************************************************************************************
!> \brief Update active space density matrix from Fortran arrays
!> \param p_act_mo_a alpha density matrix in active space MO basis
!> \param active_space_env active space environment
!> \param p_act_mo_b beta density matrix in active space MO basis
!> \author Vladimir Rybkin
! **************************************************************************************************
   SUBROUTINE update_active_density(p_act_mo_a, active_space_env, p_act_mo_b)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: p_act_mo_a
      TYPE(active_space_type), POINTER                   :: active_space_env
      REAL(KIND=dp), DIMENSION(:), INTENT(IN), OPTIONAL  :: p_act_mo_b

      INTEGER                                            :: i1, i2, idx, ispin, m1, m2, nact2, &
                                                            nmo_active, nspins
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: p_mixed, p_old, p_solver
      TYPE(cp_fm_type), POINTER                          :: p_active

      nmo_active = active_space_env%nmo_active
      nact2 = nmo_active*nmo_active
      nspins = active_space_env%nspins
      CPASSERT(SIZE(p_act_mo_a) == nact2)
      IF (nspins == 2) THEN
         IF (.NOT. PRESENT(p_act_mo_b)) CPABORT("Missing beta active-space density.")
         CPASSERT(SIZE(p_act_mo_b) == nact2)
      END IF

      ALLOCATE (p_mixed(nspins*nact2))
      ALLOCATE (p_old(nspins*nact2))
      ALLOCATE (p_solver(nspins*nact2))

      p_solver(1:nact2) = p_act_mo_a
      IF (nspins == 2) THEN
         p_solver(nact2 + 1:2*nact2) = p_act_mo_b
      END IF

      idx = 0
      DO ispin = 1, nspins
         p_active => active_space_env%p_active(ispin)
         DO i1 = 1, nmo_active
            m1 = active_space_env%active_orbitals(i1, ispin)
            DO i2 = 1, nmo_active
               idx = idx + 1
               m2 = active_space_env%active_orbitals(i2, ispin)
               CALL cp_fm_get_element(p_active, m1, m2, p_old(idx))
            END DO
         END DO
      END DO

      CALL mix_active_density_vector(active_space_env, p_old, p_solver, p_mixed)

      idx = 0
      DO ispin = 1, nspins
         p_active => active_space_env%p_active(ispin)
         DO i1 = 1, nmo_active
            m1 = active_space_env%active_orbitals(i1, ispin)
            DO i2 = 1, nmo_active
               idx = idx + 1
               m2 = active_space_env%active_orbitals(i2, ispin)
               CALL cp_fm_set_element(p_active, m1, m2, p_mixed(idx))
            END DO
         END DO
      END DO

      DEALLOCATE (p_mixed)
      DEALLOCATE (p_old)
      DEALLOCATE (p_solver)

   END SUBROUTINE update_active_density

END MODULE qs_active_space_mixing
