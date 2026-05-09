!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
MODULE qs_charge_mixing

#if defined(__TBLITE)
   USE mctc_env, ONLY: error_type
   USE tblite_scc_mixer, ONLY: new_cp2k_tblite_mixer
#endif
   USE input_constants, ONLY: tblite_mixer_damping_default, &
                              tblite_mixer_iterations_default, &
                              tblite_mixer_max_weight_default, &
                              tblite_mixer_min_weight_default, &
                              tblite_mixer_omega0_default, &
                              tblite_mixer_weight_factor_default, &
                              tblite_scc_mixer_auto, &
                              tblite_scc_mixer_cp2k, &
                              tblite_scc_mixer_none, &
                              tblite_scc_mixer_tblite
   USE kinds, ONLY: dp
   USE mathlib, ONLY: get_pseudo_inverse_svd
   USE message_passing, ONLY: mp_para_env_type
   USE qs_density_mixing_types, ONLY: broyden_mixing_nr, &
                                      gspace_mixing_nr, &
                                      mixing_storage_type, &
                                      modified_broyden_mixing_nr, &
                                      multisecant_mixing_nr, &
                                      pulay_mixing_nr
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_charge_mixing'

   PUBLIC :: charge_mixing, charge_mixing_scc_error

   REAL(KIND=dp), PARAMETER, PRIVATE :: tblite_scc_pconv = 2.0E-5_dp

CONTAINS

! **************************************************************************************************
!> \brief  Driver for TB SCC variable mixing, calls the requested method.
!> \param mixing_method ...
!> \param mixing_store ...
!> \param charges ...
!> \param para_env ...
!> \param iter_count ...
!> \param scc_mixer ...
!> \param tblite_mixer_iterations ...
!> \param tblite_mixer_damping ...
!> \param tblite_mixer_memory ...
!> \param tblite_mixer_omega0 ...
!> \param tblite_mixer_min_weight ...
!> \param tblite_mixer_max_weight ...
!> \param tblite_mixer_weight_factor ...
!> \par History
!> \author JGH
! **************************************************************************************************
   SUBROUTINE charge_mixing(mixing_method, mixing_store, charges, para_env, iter_count, &
                            scc_mixer, tblite_mixer_iterations, tblite_mixer_damping, &
                            tblite_mixer_memory, tblite_mixer_omega0, tblite_mixer_min_weight, &
                            tblite_mixer_max_weight, tblite_mixer_weight_factor)
      INTEGER, INTENT(IN)                                :: mixing_method
      TYPE(mixing_storage_type), POINTER                 :: mixing_store
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: charges
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: iter_count
      INTEGER, INTENT(IN), OPTIONAL                      :: scc_mixer, tblite_mixer_iterations
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: tblite_mixer_damping
      INTEGER, INTENT(IN), OPTIONAL                      :: tblite_mixer_memory
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: tblite_mixer_omega0, &
                                                            tblite_mixer_min_weight, &
                                                            tblite_mixer_max_weight, &
                                                            tblite_mixer_weight_factor

      CHARACTER(len=*), PARAMETER                        :: routineN = 'charge_mixing'

      INTEGER                                            :: effective_scc_mixer, handle, ia, ii, &
                                                            imin, inow, mixer_iterations, &
                                                            mixer_memory, nbuffer, ns, nvec
      REAL(dp)                                           :: alpha, mixer_damping, mixer_max_weight, &
                                                            mixer_min_weight, mixer_omega0, &
                                                            mixer_weight_factor

      CALL timeset(routineN, handle)

      effective_scc_mixer = tblite_scc_mixer_cp2k
      IF (PRESENT(scc_mixer)) effective_scc_mixer = scc_mixer
      IF (ASSOCIATED(mixing_store)) mixing_store%tb_scc_mixer_error = 0.0_dp

      SELECT CASE (effective_scc_mixer)
      CASE (tblite_scc_mixer_auto, tblite_scc_mixer_cp2k)
         ! Use the regular CP2K SCC-variable mixing path below.
      CASE (tblite_scc_mixer_tblite)
         CPASSERT(ASSOCIATED(mixing_store))
#if defined(__TBLITE)
         mixer_damping = tblite_mixer_damping_default
         IF (PRESENT(tblite_mixer_damping)) mixer_damping = tblite_mixer_damping
         IF (mixer_damping <= 0.0_dp) CPABORT("tblite SCC mixer DAMPING must be positive")
         mixer_omega0 = tblite_mixer_omega0_default
         IF (PRESENT(tblite_mixer_omega0)) mixer_omega0 = tblite_mixer_omega0
         IF (mixer_omega0 <= 0.0_dp) CPABORT("tblite SCC mixer OMEGA0 must be positive")
         mixer_min_weight = tblite_mixer_min_weight_default
         IF (PRESENT(tblite_mixer_min_weight)) mixer_min_weight = tblite_mixer_min_weight
         IF (mixer_min_weight <= 0.0_dp) CPABORT("tblite SCC mixer MIN_WEIGHT must be positive")
         mixer_max_weight = tblite_mixer_max_weight_default
         IF (PRESENT(tblite_mixer_max_weight)) mixer_max_weight = tblite_mixer_max_weight
         IF (mixer_max_weight <= 0.0_dp) CPABORT("tblite SCC mixer MAX_WEIGHT must be positive")
         IF (mixer_max_weight < mixer_min_weight) &
            CPABORT("tblite SCC mixer MAX_WEIGHT must not be smaller than MIN_WEIGHT")
         mixer_weight_factor = tblite_mixer_weight_factor_default
         IF (PRESENT(tblite_mixer_weight_factor)) mixer_weight_factor = tblite_mixer_weight_factor
         IF (mixer_weight_factor <= 0.0_dp) CPABORT("tblite SCC mixer WEIGHT_FACTOR must be positive")
         mixer_iterations = tblite_mixer_iterations_default
         IF (PRESENT(tblite_mixer_iterations)) mixer_iterations = tblite_mixer_iterations
         IF (mixer_iterations < 1) CPABORT("tblite SCC mixer ITERATIONS must be positive")
         IF (iter_count > mixer_iterations) CPABORT("tblite SCC mixer exceeded ITERATIONS")
         mixer_memory = MAX(1, mixing_store%nbuffer)
         IF (PRESENT(tblite_mixer_memory)) mixer_memory = tblite_mixer_memory
         IF (mixer_memory < 1) CPABORT("tblite SCC mixer MEMORY must be positive")
         CALL tblite_charge_mixing(mixing_store, charges, para_env, iter_count, &
                                   mixer_damping, mixer_memory, mixer_omega0, mixer_min_weight, &
                                   mixer_max_weight, mixer_weight_factor)
         CALL timestop(handle)
         RETURN
#else
         IF (iter_count == 1) THEN
            CALL cp_warn(__LOCATION__, &
                         "SCC_MIXER TBLITE requested but CP2K was built without tblite; "// &
                         "falling back to the CP2K SCC mixer.")
         END IF
#endif
      CASE (tblite_scc_mixer_none)
         IF (ASSOCIATED(mixing_store)) mixing_store%iter_method = "NoMix"
         CALL timestop(handle)
         RETURN
      CASE DEFAULT
         CPABORT("Unknown SCC mixer for TB charge mixing")
      END SELECT

      IF (mixing_method >= gspace_mixing_nr) THEN
         CPASSERT(ASSOCIATED(mixing_store))
         mixing_store%ncall = mixing_store%ncall + 1
         ns = SIZE(charges, 2)
         IF (ns > mixing_store%max_shell) THEN
            CPABORT("Mixing storage too small for TB SCC variables")
         END IF
         alpha = mixing_store%alpha
         nbuffer = mixing_store%nbuffer
         inow = MOD(mixing_store%ncall - 1, nbuffer) + 1
         imin = inow - 1
         IF (imin == 0) imin = nbuffer
         IF (mixing_store%ncall > nbuffer) THEN
            nvec = nbuffer
         ELSE
            nvec = mixing_store%ncall - 1
         END IF
         IF (mixing_store%ncall > 1) THEN
            ! store in/out charge difference
            DO ia = 1, mixing_store%nat_local
               ii = mixing_store%atlist(ia)
               mixing_store%dacharge(ia, 1:ns, imin) = mixing_store%acharge(ia, 1:ns, imin) - charges(ii, 1:ns)
            END DO
         END IF
         IF ((iter_count == 1) .OR. (iter_count + 1 <= mixing_store%nskip_mixing)) THEN
            ! skip mixing
            mixing_store%iter_method = "NoMix"
         ELSEIF (((iter_count + 1 - mixing_store%nskip_mixing) <= mixing_store%n_simple_mix) .OR. (nvec == 1)) THEN
            CALL mix_charges_only(mixing_store, charges, alpha, imin, ns, para_env)
            mixing_store%iter_method = "Mixing"
         ELSEIF (mixing_method == gspace_mixing_nr) THEN
            CPABORT("Kerker method not available for Charge Mixing")
         ELSEIF (mixing_method == pulay_mixing_nr) THEN
            CPABORT("Pulay method not available for Charge Mixing")
         ELSEIF (mixing_method == broyden_mixing_nr) THEN
            CALL broyden_mixing(mixing_store, charges, imin, nvec, ns, para_env)
            mixing_store%iter_method = "Broy."
         ELSEIF (mixing_method == modified_broyden_mixing_nr) THEN
            CPABORT("Modified Broyden mixing is only available for DFT density mixing")
         ELSEIF (mixing_method == multisecant_mixing_nr) THEN
            CPABORT("Multisecant_mixing method not available for Charge Mixing")
         END IF

         ! store new 'input' charges
         DO ia = 1, mixing_store%nat_local
            ii = mixing_store%atlist(ia)
            mixing_store%acharge(ia, 1:ns, inow) = charges(ii, 1:ns)
         END DO

      END IF

      CALL timestop(handle)

   END SUBROUTINE charge_mixing

! **************************************************************************************************
!> \brief Return the CP2K-side tblite SCC-mixer residual on the CP2K EPS_SCF scale.
!> \param mixing_store ...
!> \param eps_scf ...
!> \return ...
! **************************************************************************************************
   FUNCTION charge_mixing_scc_error(mixing_store, eps_scf) RESULT(mixer_error)
      TYPE(mixing_storage_type), POINTER                 :: mixing_store
      REAL(KIND=dp), INTENT(IN)                          :: eps_scf
      REAL(KIND=dp)                                      :: mixer_error

      mixer_error = 0.0_dp
      IF (.NOT. ASSOCIATED(mixing_store)) RETURN
      IF (mixing_store%tb_scc_mixer_step <= 1) RETURN

      IF (eps_scf > 0.0_dp) THEN
         mixer_error = eps_scf*mixing_store%tb_scc_mixer_error/tblite_scc_pconv
      ELSE
         mixer_error = mixing_store%tb_scc_mixer_error
      END IF

   END FUNCTION charge_mixing_scc_error

! **************************************************************************************************
!> \brief TBLite modified-Broyden mixing for a complete TB SCC-variable vector.
!> \param mixing_store ...
!> \param charges ...
!> \param para_env ...
!> \param iter_count ...
!> \param damping ...
!> \param memory ...
!> \param omega0 ...
!> \param min_weight ...
!> \param max_weight ...
!> \param weight_factor ...
! **************************************************************************************************
   SUBROUTINE tblite_charge_mixing(mixing_store, charges, para_env, iter_count, damping, memory, omega0, &
                                   min_weight, max_weight, weight_factor)
      TYPE(mixing_storage_type), POINTER                 :: mixing_store
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: charges
      TYPE(mp_para_env_type), POINTER                    :: para_env
      INTEGER, INTENT(IN)                                :: iter_count, memory
      REAL(KIND=dp), INTENT(IN)                          :: damping, max_weight, min_weight, omega0, &
                                                            weight_factor

#if defined(__TBLITE)
      TYPE(error_type), ALLOCATABLE                      :: error
#endif
      INTEGER                                            :: natom, ndim, ns
      LOGICAL                                            :: on_source, reset_mixer
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: qvec

      natom = SIZE(charges, 1)
      ns = SIZE(charges, 2)
      ndim = natom*ns
      ALLOCATE (qvec(ndim))
      qvec(:) = RESHAPE(charges, [ndim])
      on_source = para_env%mepos == para_env%source
      reset_mixer = (iter_count == 1) .OR. (mixing_store%tb_scc_mixer_step == 0) .OR. &
                    (mixing_store%tb_scc_mixer_natom /= natom) .OR. &
                    (mixing_store%tb_scc_mixer_ns /= ns) .OR. &
                    (mixing_store%tb_scc_mixer_memory /= memory)
      mixing_store%tb_scc_mixer_error = 0.0_dp

#if defined(__TBLITE)
      IF (reset_mixer) THEN
         IF (ALLOCATED(mixing_store%tb_scc_mixer)) DEALLOCATE (mixing_store%tb_scc_mixer)
         IF (on_source) THEN
            CALL new_cp2k_tblite_mixer(mixing_store%tb_scc_mixer, memory, ndim, damping, omega0, &
                                       min_weight, max_weight, weight_factor)
            CALL mixing_store%tb_scc_mixer%set(qvec)
         END IF
         mixing_store%tb_scc_mixer_natom = natom
         mixing_store%tb_scc_mixer_ns = ns
         mixing_store%tb_scc_mixer_memory = memory
         mixing_store%tb_scc_mixer_step = 1
         mixing_store%iter_method = "NoMix"
         CALL para_env%bcast(qvec)
         charges = RESHAPE(qvec, SHAPE(charges))
         RETURN
      END IF

      IF (on_source) THEN
         CPASSERT(ALLOCATED(mixing_store%tb_scc_mixer))
         CALL mixing_store%tb_scc_mixer%diff(qvec)
         mixing_store%tb_scc_mixer_error = REAL(mixing_store%tb_scc_mixer%get_error(), KIND=dp)
         CALL mixing_store%tb_scc_mixer%next(error)
         IF (ALLOCATED(error)) CPABORT("tblite SCC mixer failed")
         CALL mixing_store%tb_scc_mixer%get(qvec)
      END IF
      CALL para_env%bcast(qvec)
      CALL para_env%bcast(mixing_store%tb_scc_mixer_error)
      charges = RESHAPE(qvec, SHAPE(charges))
      mixing_store%tb_scc_mixer_step = mixing_store%tb_scc_mixer_step + 1
      mixing_store%iter_method = "TBLITE"
#else
      MARK_USED(mixing_store)
      MARK_USED(charges)
      MARK_USED(para_env)
      MARK_USED(iter_count)
      MARK_USED(damping)
      MARK_USED(memory)
      MARK_USED(omega0)
      MARK_USED(min_weight)
      MARK_USED(max_weight)
      MARK_USED(weight_factor)
      CPABORT("SCC_MIXER TBLITE requires CP2K to be built with tblite")
#endif

   END SUBROUTINE tblite_charge_mixing

! **************************************************************************************************
!> \brief Simple charge mixing
!> \param mixing_store ...
!> \param charges ...
!> \param alpha ...
!> \param imin ...
!> \param ns ...
!> \param para_env ...
!> \author JGH
! **************************************************************************************************
   SUBROUTINE mix_charges_only(mixing_store, charges, alpha, imin, ns, para_env)
      TYPE(mixing_storage_type), POINTER                 :: mixing_store
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: charges
      REAL(KIND=dp), INTENT(IN)                          :: alpha
      INTEGER, INTENT(IN)                                :: imin, ns
      TYPE(mp_para_env_type), POINTER                    :: para_env

      INTEGER                                            :: ia, ii

      charges = 0.0_dp

      DO ia = 1, mixing_store%nat_local
         ii = mixing_store%atlist(ia)
         charges(ii, 1:ns) = alpha*mixing_store%dacharge(ia, 1:ns, imin) - mixing_store%acharge(ia, 1:ns, imin)
      END DO

      CALL para_env%sum(charges)

   END SUBROUTINE mix_charges_only

! **************************************************************************************************
!> \brief Broyden charge mixing
!> \param mixing_store ...
!> \param charges ...
!> \param inow ...
!> \param nvec ...
!> \param ns ...
!> \param para_env ...
!> \author JGH
! **************************************************************************************************
   SUBROUTINE broyden_mixing(mixing_store, charges, inow, nvec, ns, para_env)
      TYPE(mixing_storage_type), POINTER                 :: mixing_store
      REAL(KIND=dp), DIMENSION(:, :), INTENT(INOUT)      :: charges
      INTEGER, INTENT(IN)                                :: inow, nvec, ns
      TYPE(mp_para_env_type), POINTER                    :: para_env

      INTEGER                                            :: i, ia, ii, imin, j, nbuffer, nv
      REAL(KIND=dp)                                      :: alpha, broy_w0, rskip, wdf, wprod
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: cvec, gammab
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: amat, beta
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: dq_last, dq_now, q_last, q_now

      CPASSERT(nvec > 1)

      nbuffer = mixing_store%nbuffer
      alpha = mixing_store%alpha
      imin = inow - 1
      IF (imin == 0) imin = nvec
      nv = nvec - 1

      ! charge vectors
      q_now => mixing_store%acharge(:, :, inow)
      q_last => mixing_store%acharge(:, :, imin)
      dq_now => mixing_store%dacharge(:, :, inow)
      dq_last => mixing_store%dacharge(:, :, imin)

      IF (nvec == nbuffer) THEN
         ! reshuffel Broyden storage n->n-1
         DO i = 1, nv - 1
            mixing_store%wbroy(i) = mixing_store%wbroy(i + 1)
            mixing_store%dfbroy(:, :, i) = mixing_store%dfbroy(:, :, i + 1)
            mixing_store%ubroy(:, :, i) = mixing_store%ubroy(:, :, i + 1)
         END DO
         DO i = 1, nv - 1
            DO j = 1, nv - 1
               mixing_store%abroy(i, j) = mixing_store%abroy(i + 1, j + 1)
            END DO
         END DO
      END IF

      broy_w0 = mixing_store%broy_w0
      mixing_store%wbroy(nv) = 1.0_dp

      ! dfbroy
      mixing_store%dfbroy(:, :, nv) = 0.0_dp
      mixing_store%dfbroy(:, 1:ns, nv) = dq_now(:, 1:ns) - dq_last(:, 1:ns)
      wdf = SUM(mixing_store%dfbroy(:, 1:ns, nv)**2)
      CALL para_env%sum(wdf)
      wdf = 1.0_dp/SQRT(wdf)
      mixing_store%dfbroy(:, 1:ns, nv) = wdf*mixing_store%dfbroy(:, 1:ns, nv)

      ! abroy matrix
      DO i = 1, nv
         wprod = SUM(mixing_store%dfbroy(:, 1:ns, i)*mixing_store%dfbroy(:, 1:ns, nv))
         CALL para_env%sum(wprod)
         mixing_store%abroy(i, nv) = wprod
         mixing_store%abroy(nv, i) = wprod
      END DO

      ! broyden matrices
      ALLOCATE (amat(nv, nv), beta(nv, nv), cvec(nv), gammab(nv))
      DO i = 1, nv
         wprod = SUM(mixing_store%dfbroy(:, 1:ns, i)*dq_now(:, 1:ns))
         CALL para_env%sum(wprod)
         cvec(i) = mixing_store%wbroy(i)*wprod
      END DO

      DO i = 1, nv
         DO j = 1, nv
            beta(j, i) = mixing_store%wbroy(j)*mixing_store%wbroy(i)*mixing_store%abroy(j, i)
         END DO
         beta(i, i) = beta(i, i) + broy_w0
      END DO

      rskip = 1.e-12_dp
      CALL get_pseudo_inverse_svd(beta, amat, rskip)
      gammab(1:nv) = MATMUL(cvec(1:nv), amat(1:nv, 1:nv))

      ! build ubroy
      mixing_store%ubroy(:, :, nv) = 0.0_dp
      mixing_store%ubroy(:, 1:ns, nv) = alpha*mixing_store%dfbroy(:, 1:ns, nv) + &
                                        wdf*(q_now(:, 1:ns) - q_last(:, 1:ns))

      charges = 0.0_dp
      DO ia = 1, mixing_store%nat_local
         ii = mixing_store%atlist(ia)
         charges(ii, 1:ns) = q_now(ia, 1:ns) + alpha*dq_now(ia, 1:ns)
      END DO
      DO i = 1, nv
         DO ia = 1, mixing_store%nat_local
            ii = mixing_store%atlist(ia)
            charges(ii, 1:ns) = charges(ii, 1:ns) - mixing_store%wbroy(i)*gammab(i)*mixing_store%ubroy(ia, 1:ns, i)
         END DO
      END DO
      CALL para_env%sum(charges)

      DEALLOCATE (amat, beta, cvec, gammab)

   END SUBROUTINE broyden_mixing

END MODULE qs_charge_mixing
