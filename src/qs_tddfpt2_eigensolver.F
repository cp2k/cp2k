!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

MODULE qs_tddfpt2_eigensolver
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: tddfpt2_control_type
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply
   USE cp_fm_basic_linalg,              ONLY: cp_fm_contracted_trace,&
                                              cp_fm_scale,&
                                              cp_fm_scale_and_add,&
                                              cp_fm_trace
   USE cp_fm_diag,                      ONLY: choose_eigv_solver
   USE cp_fm_pool_types,                ONLY: fm_pool_create_fm,&
                                              fm_pool_give_back_fm
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: &
        cp_fm_copy_general, cp_fm_create, cp_fm_get_info, cp_fm_get_submatrix, cp_fm_maxabsval, &
        cp_fm_p_type, cp_fm_release, cp_fm_set_all, cp_fm_set_submatrix, cp_fm_to_fm, cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_log_handling,                 ONLY: cp_logger_type
   USE cp_output_handling,              ONLY: cp_iterate
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_get_info,&
                                              dbcsr_p_type,&
                                              dbcsr_type
   USE input_constants,                 ONLY: tddfpt_kernel_full,&
                                              tddfpt_kernel_stda
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE machine,                         ONLY: m_flush,&
                                              m_walltime
   USE memory_utilities,                ONLY: reallocate
   USE physcon,                         ONLY: evolt
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_scf_methods,                  ONLY: eigensolver
   USE qs_tddfpt2_fhxc,                 ONLY: fhxc_kernel,&
                                              stda_kernel
   USE qs_tddfpt2_operators,            ONLY: tddfpt_apply_energy_diff,&
                                              tddfpt_apply_hfx
   USE qs_tddfpt2_restart,              ONLY: tddfpt_write_restart
   USE qs_tddfpt2_subgroups,            ONLY: tddfpt_subgroup_env_type
   USE qs_tddfpt2_types,                ONLY: full_kernel_env_type,&
                                              tddfpt_ground_state_mos,&
                                              tddfpt_kernel_env_type,&
                                              tddfpt_work_matrices
   USE qs_tddfpt2_utils,                ONLY: tddfpt_total_number_of_states
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddfpt2_eigensolver'

   LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .FALSE.
   ! number of first derivative components (3: d/dx, d/dy, d/dz)
   INTEGER, PARAMETER, PRIVATE          :: nderivs = 3
   INTEGER, PARAMETER, PRIVATE          :: maxspins = 2

   PUBLIC :: tddfpt_davidson_solver, tddfpt_orthogonalize_psi1_psi0, &
             tddfpt_orthonormalize_psi1_psi1

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief Make TDDFPT trial vectors orthogonal to all occupied molecular orbitals.
!> \param evects            trial vectors distributed across all processors (modified on exit)
!> \param S_C0_C0T          matrix product S * C_0 * C_0^T, where C_0 is the ground state
!>                          wave function for each spin expressed in atomic basis set,
!>                          and S is the corresponding overlap matrix
!> \par History
!>    * 05.2016 created [Sergey Chulkov]
!>    * 05.2019 use a temporary work matrix [JHU]
!> \note  Based on the subroutine p_preortho() which was created by Thomas Chassaing on 09.2002.
!>        Should be useless when ground state MOs are computed with extremely high accuracy,
!>        as all virtual orbitals are already orthogonal to the occupied ones by design.
!>        However, when the norm of residual vectors is relatively small (e.g. less then SCF_EPS),
!>        new Krylov's vectors seem to be random and should be orthogonalised even with respect to
!>        the occupied MOs.
! **************************************************************************************************
   SUBROUTINE tddfpt_orthogonalize_psi1_psi0(evects, S_C0_C0T)
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: evects
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(in)       :: S_C0_C0T

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_orthogonalize_psi1_psi0', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, ivect, nactive, nao, &
                                                            nspins, nvects
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct
      TYPE(cp_fm_type), POINTER                          :: evortho

      CALL timeset(routineN, handle)

      nspins = SIZE(evects, 1)
      nvects = SIZE(evects, 2)

      IF (nvects > 0) THEN
         DO ispin = 1, nspins
            CALL cp_fm_get_info(matrix=evects(ispin, 1)%matrix, matrix_struct=matrix_struct, &
                                nrow_global=nao, ncol_global=nactive)
            NULLIFY (evortho)
            CALL cp_fm_create(evortho, matrix_struct)
            DO ivect = 1, nvects
               ! evortho: C0 * C0^T * S * C1 == (S * C0 * C0^T)^T * C1
               CALL cp_gemm('T', 'N', nao, nactive, nao, 1.0_dp, S_C0_C0T(ispin)%matrix, &
                            evects(ispin, ivect)%matrix, 0.0_dp, evortho)
               CALL cp_fm_scale_and_add(1.0_dp, evects(ispin, ivect)%matrix, -1.0_dp, evortho)
            END DO
            CALL cp_fm_release(evortho)
         END DO
      END IF

      CALL timestop(handle)

   END SUBROUTINE tddfpt_orthogonalize_psi1_psi0

! **************************************************************************************************
!> \brief Check that orthogonalised TDDFPT trial vectors remain orthogonal to
!>        occupied molecular orbitals.
!> \param evects    trial vectors
!> \param S_C0      matrix product S * C_0, where C_0 is the ground state wave function
!>                  for each spin in atomic basis set, and S is the corresponding overlap matrix
!> \param max_norm  the largest possible overlap between the ground state and
!>                  excited state wave functions
!> \return true if trial vectors are non-orthogonal to occupied molecular orbitals
!> \par History
!>    * 07.2016 created [Sergey Chulkov]
!>    * 05.2019 use temporary work matrices [JHU]
! **************************************************************************************************
   FUNCTION tddfpt_is_nonorthogonal_psi1_psi0(evects, S_C0, max_norm) &
      RESULT(is_nonortho)
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: evects
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(in)       :: S_C0
      REAL(kind=dp), INTENT(in)                          :: max_norm
      LOGICAL                                            :: is_nonortho

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_is_nonorthogonal_psi1_psi0', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, ivect, nactive, nao, &
                                                            nocc, nspins, nvects
      REAL(kind=dp)                                      :: maxabs_val
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct, matrix_struct_tmp
      TYPE(cp_fm_type), POINTER                          :: aortho

      CALL timeset(routineN, handle)

      nspins = SIZE(evects, 1)
      nvects = SIZE(evects, 2)

      is_nonortho = .FALSE.

      loop: DO ispin = 1, nspins
         CALL cp_fm_get_info(matrix=S_C0(ispin)%matrix, ncol_global=nocc)
         CALL cp_fm_get_info(matrix=evects(ispin, 1)%matrix, matrix_struct=matrix_struct, &
                             nrow_global=nao, ncol_global=nactive)
         CALL cp_fm_struct_create(matrix_struct_tmp, nrow_global=nocc, &
                                  ncol_global=nactive, template_fmstruct=matrix_struct)
         NULLIFY (aortho)
         CALL cp_fm_create(aortho, matrix_struct_tmp)
         CALL cp_fm_struct_release(matrix_struct_tmp)
         DO ivect = 1, nvects
            ! aortho = S_0^T * S * C_1
            CALL cp_gemm('T', 'N', nocc, nactive, nao, 1.0_dp, S_C0(ispin)%matrix, &
                         evects(ispin, ivect)%matrix, 0.0_dp, aortho)
            CALL cp_fm_maxabsval(aortho, maxabs_val)
            is_nonortho = maxabs_val > max_norm
            IF (is_nonortho) THEN
               CALL cp_fm_release(aortho)
               EXIT loop
            END IF
         END DO
         CALL cp_fm_release(aortho)
      END DO loop

      CALL timestop(handle)

   END FUNCTION tddfpt_is_nonorthogonal_psi1_psi0

! **************************************************************************************************
!> \brief Make new TDDFPT trial vectors orthonormal to all previous trial vectors.
!> \param evects      trial vectors (modified on exit)
!> \param nvects_new  number of new trial vectors to orthogonalise
!> \param S_evects    set of matrices to store matrix product S * evects (modified on exit)
!> \param matrix_s    overlap matrix
!> \par History
!>    * 05.2016 created [Sergey Chulkov]
!>    * 02.2017 caching the matrix product S * evects [Sergey Chulkov]
!> \note \parblock
!>       Based on the subroutines reorthogonalize() and normalize() which were originally created
!>       by Thomas Chassaing on 03.2003.
!>
!>       In order to orthogonalise a trial vector C3 = evects(:,3) with respect to previously
!>       orthogonalised vectors C1 = evects(:,1) and C2 = evects(:,2) we need to compute the
!>       quantity C3'' using the following formulae:
!>          C3'  = C3  - Tr(C3^T  * S * C1) * C1,
!>          C3'' = C3' - Tr(C3'^T * S * C2) * C2,
!>       which can be expanded as:
!>          C3'' = C3 - Tr(C3^T  * S * C1) * C1 - Tr(C3^T * S * C2) * C2 +
!>                 Tr(C3^T * S * C1) * Tr(C2^T * S * C1) * C2 .
!>       In case of unlimited float-point precision, the last term in above expression is exactly 0,
!>       due to orthogonality condition between C1 and C2. In this case the expression could be
!>       simplified as (taking into account the identity: Tr(A * S * B) = Tr(B * S * A)):
!>          C3'' = C3 - Tr(C1^T  * S * C3) * C1 - Tr(C2^T * S * C3) * C2 ,
!>       which means we do not need the variable S_evects to keep the matrix products S * Ci .
!>
!>       In reality, however, we deal with limited float-point precision arithmetic meaning that
!>       the trace Tr(C2^T * S * C1) is close to 0 but does not equal to 0 exactly. The term
!>          Tr(C3^T * S * C1) * Tr(C2^T * S * C1) * C2
!>       can not be ignored anymore. Ignorance of this term will lead to numerical instability
!>       when the trace Tr(C3^T * S * C1) is large enough.
!>       \endparblock
! **************************************************************************************************
   SUBROUTINE tddfpt_orthonormalize_psi1_psi1(evects, nvects_new, S_evects, matrix_s)
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: evects
      INTEGER, INTENT(in)                                :: nvects_new
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: S_evects
      TYPE(dbcsr_type), POINTER                          :: matrix_s

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_orthonormalize_psi1_psi1', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, ivect, jvect, nspins, &
                                                            nvects_old, nvects_total
      INTEGER, DIMENSION(maxspins)                       :: nactive
      REAL(kind=dp)                                      :: norm
      REAL(kind=dp), DIMENSION(maxspins)                 :: weights

      CALL timeset(routineN, handle)

      nspins = SIZE(evects, 1)
      nvects_total = SIZE(evects, 2)
      nvects_old = nvects_total-nvects_new

      IF (debug_this_module) THEN
         CPASSERT(SIZE(S_evects, 1) == nspins)
         CPASSERT(SIZE(S_evects, 2) == nvects_total)
         CPASSERT(nvects_old >= 0)
      END IF

      DO ispin = 1, nspins
         CALL cp_fm_get_info(matrix=evects(ispin, 1)%matrix, ncol_global=nactive(ispin))
      END DO

      DO jvect = nvects_old+1, nvects_total
         ! <psi1_i | psi1_j>
         DO ivect = 1, jvect-1
            CALL cp_fm_trace(evects(:, jvect), S_evects(:, ivect), weights(1:nspins), accurate=.FALSE.)
            norm = SUM(weights(1:nspins))

            DO ispin = 1, nspins
               CALL cp_fm_scale_and_add(1.0_dp, evects(ispin, jvect)%matrix, -norm, evects(ispin, ivect)%matrix)
            END DO
         END DO

         ! <psi1_j | psi1_j>
         DO ispin = 1, nspins
            CALL cp_dbcsr_sm_fm_multiply(matrix_s, evects(ispin, jvect)%matrix, S_evects(ispin, jvect)%matrix, &
                                         ncol=nactive(ispin), alpha=1.0_dp, beta=0.0_dp)
         END DO

         CALL cp_fm_trace(evects(:, jvect), S_evects(:, jvect), weights(1:nspins), accurate=.FALSE.)

         norm = SUM(weights(1:nspins))
         norm = 1.0_dp/SQRT(norm)

         DO ispin = 1, nspins
            CALL cp_fm_scale(norm, evects(ispin, jvect)%matrix)
            CALL cp_fm_scale(norm, S_evects(ispin, jvect)%matrix)
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE tddfpt_orthonormalize_psi1_psi1

! **************************************************************************************************
!> \brief Compute action matrix-vector products.
!> \param Aop_evects            action of TDDFPT operator on trial vectors (modified on exit)
!> \param evects                TDDFPT trial vectors
!> \param S_evects              cached matrix product S * evects where S is the overlap matrix
!>                              in primary basis set
!> \param gs_mos                molecular orbitals optimised for the ground state
!> \param tddfpt_control        control section for tddfpt
!> \param do_hfx                flag that activates computation of exact-exchange terms
!> \param matrix_ks             Kohn-Sham matrix
!> \param qs_env                Quickstep environment
!> \param kernel_env            kernel environment
!> \param sub_env               parallel (sub)group environment
!> \param work_matrices         collection of work matrices (modified on exit)
!> \par History
!>    * 06.2016 created [Sergey Chulkov]
!>    * 03.2017 refactored [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_compute_Aop_evects(Aop_evects, evects, S_evects, gs_mos, tddfpt_control, &
                                        do_hfx, matrix_ks, qs_env, kernel_env, &
                                        sub_env, work_matrices)
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: Aop_evects, evects, S_evects
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control
      LOGICAL, INTENT(in)                                :: do_hfx
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(in)       :: matrix_ks
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddfpt_kernel_env_type), INTENT(in)           :: kernel_env
      TYPE(tddfpt_subgroup_env_type), INTENT(in)         :: sub_env
      TYPE(tddfpt_work_matrices), INTENT(inout)          :: work_matrices

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_compute_Aop_evects', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, ivect, nspins, nvects
      INTEGER, DIMENSION(maxspins)                       :: nmo_occ
      LOGICAL                                            :: do_admm, is_rks_triplets, is_stda
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(full_kernel_env_type), POINTER                :: full_kernel_env, kernel_env_admm_aux

      CALL timeset(routineN, handle)

      full_kernel_env => kernel_env%full_kernel

      nspins = SIZE(evects, 1)
      nvects = SIZE(evects, 2)
      do_admm = ASSOCIATED(sub_env%admm_A)
      IF (do_admm) THEN
         kernel_env_admm_aux => kernel_env%admm_kernel
      ELSE
         NULLIFY (kernel_env_admm_aux)
      END IF
      is_rks_triplets = tddfpt_control%rks_triplets

      IF (debug_this_module) THEN
         CPASSERT(nspins > 0)
         CPASSERT(SIZE(Aop_evects, 1) == nspins)
         CPASSERT(SIZE(Aop_evects, 2) == nvects)
         CPASSERT(SIZE(S_evects, 1) == nspins)
         CPASSERT(SIZE(S_evects, 2) == nvects)
         CPASSERT(SIZE(gs_mos) == nspins)
      END IF

      DO ispin = 1, nspins
         nmo_occ(ispin) = SIZE(gs_mos(ispin)%evals_occ)
      END DO

      IF (nvects > 0) THEN
         CALL cp_fm_get_info(evects(1, 1)%matrix, para_env=para_env)
         IF (ALLOCATED(work_matrices%evects_sub)) THEN
            DO ivect = 1, nvects
               DO ispin = 1, nspins
                  CALL cp_fm_copy_general(evects(ispin, ivect)%matrix, &
                                          work_matrices%evects_sub(ispin, ivect)%matrix, para_env)
               END DO
            END DO
         END IF

         IF (tddfpt_control%kernel == tddfpt_kernel_full) THEN
            ! full TDDFPT kernel
            CALL fhxc_kernel(Aop_evects, evects, is_rks_triplets, do_hfx, qs_env, &
                             full_kernel_env, kernel_env_admm_aux, sub_env, work_matrices)
            is_stda = .FALSE.
         ELSE IF (tddfpt_control%kernel == tddfpt_kernel_stda) THEN
            ! sTDA kernel
            CALL stda_kernel(Aop_evects, evects, is_rks_triplets, qs_env, tddfpt_control%stda_control, &
                             kernel_env%stda_kernel, sub_env, work_matrices)
            is_stda = .TRUE.
         ELSE
            CPABORT("Kernel type not implemented")
         END IF

         IF (ALLOCATED(work_matrices%evects_sub)) THEN
            DO ivect = 1, nvects
               DO ispin = 1, nspins
                  CALL cp_fm_copy_general(work_matrices%Aop_evects_sub(ispin, ivect)%matrix, &
                                          Aop_evects(ispin, ivect)%matrix, para_env)
               END DO
            END DO
         END IF

         ! orbital energy difference term
         CALL tddfpt_apply_energy_diff(Aop_evects=Aop_evects, evects=evects, S_evects=S_evects, &
                                       gs_mos=gs_mos, matrix_ks=matrix_ks)

         IF (do_hfx .AND. .NOT. is_stda) THEN
            CALL tddfpt_apply_hfx(Aop_evects=Aop_evects, evects=evects, gs_mos=gs_mos, do_admm=do_admm, &
                                  qs_env=qs_env, work_rho_ia_ao=work_matrices%hfx_rho_ao, &
                                  work_hmat=work_matrices%hfx_hmat, wfm_rho_orb=work_matrices%hfx_fm_ao_ao)
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE tddfpt_compute_Aop_evects

! **************************************************************************************************
!> \brief Solve eigenproblem for the reduced action matrix and find new Ritz eigenvectors and
!>        eigenvalues.
!> \param ritz_vects       Ritz eigenvectors (initialised on exit)
!> \param Aop_ritz         approximate action of TDDFPT operator on Ritz vectors (initialised on exit)
!> \param evals            Ritz eigenvalues (initialised on exit)
!> \param krylov_vects     Krylov's vectors
!> \param Aop_krylov       action of TDDFPT operator on Krylov's vectors
!> \param Atilde ...
!> \param nkvo ...
!> \param nkvn ...
!> \par History
!>    * 06.2016 created [Sergey Chulkov]
!>    * 03.2017 altered prototype, OpenMP parallelisation [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_compute_ritz_vects(ritz_vects, Aop_ritz, evals, krylov_vects, Aop_krylov, &
                                        Atilde, nkvo, nkvn)
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: ritz_vects, Aop_ritz
      REAL(kind=dp), DIMENSION(:), INTENT(out)           :: evals
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: krylov_vects, Aop_krylov
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: Atilde
      INTEGER, INTENT(IN)                                :: nkvo, nkvn

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_compute_ritz_vects', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ikv, irv, ispin, nkvs, nrvs, &
                                                            nspins
      REAL(kind=dp)                                      :: act
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: at12, at21, at22, evects_Atilde
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_global
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: Atilde_fm, evects_Atilde_fm

      CALL timeset(routineN, handle)

      nspins = SIZE(krylov_vects, 1)
      nkvs = SIZE(krylov_vects, 2)
      nrvs = SIZE(ritz_vects, 2)
      CPASSERT(nkvs == nkvo+nkvn)

      CALL cp_fm_get_info(krylov_vects(1, 1)%matrix, context=blacs_env_global)

      CALL cp_fm_struct_create(fm_struct, nrow_global=nkvs, ncol_global=nkvs, context=blacs_env_global)
      CALL cp_fm_create(Atilde_fm, fm_struct)
      CALL cp_fm_create(evects_Atilde_fm, fm_struct)
      CALL cp_fm_struct_release(fm_struct)

      ! *** compute upper-diagonal reduced action matrix ***
      CALL reallocate(Atilde, 1, nkvs, 1, nkvs)
      ! TO DO: the subroutine 'cp_fm_contracted_trace' will compute all elements of
      ! the matrix 'Atilde', however only upper-triangular elements are actually needed
      !
      IF (nkvo == 0) THEN
         CALL cp_fm_contracted_trace(Aop_krylov(:, 1:nkvs), krylov_vects(:, 1:nkvs), &
                                     Atilde(1:nkvs, 1:nkvs), accurate=.FALSE.)
      ELSE
         ALLOCATE (at12(nkvn, nkvo), at21(nkvo, nkvn), at22(nkvn, nkvn))
         CALL cp_fm_contracted_trace(Aop_krylov(:, nkvo+1:nkvs), krylov_vects(:, 1:nkvo), &
                                     at12, accurate=.FALSE.)
         Atilde(nkvo+1:nkvs, 1:nkvo) = at12(1:nkvn, 1:nkvo)
         CALL cp_fm_contracted_trace(Aop_krylov(:, 1:nkvo), krylov_vects(:, nkvo+1:nkvs), &
                                     at21, accurate=.FALSE.)
         Atilde(1:nkvo, nkvo+1:nkvs) = at21(1:nkvo, 1:nkvn)
         CALL cp_fm_contracted_trace(Aop_krylov(:, nkvo+1:nkvs), krylov_vects(:, nkvo+1:nkvs), &
                                     at22, accurate=.FALSE.)
         Atilde(nkvo+1:nkvs, nkvo+1:nkvs) = at22(1:nkvn, 1:nkvn)
         DEALLOCATE (at12, at21, at22)
      END IF
      Atilde = 0.5_dp*(Atilde+TRANSPOSE(Atilde))
      CALL cp_fm_set_submatrix(Atilde_fm, Atilde)

      ! *** solve an eigenproblem for the reduced matrix ***
      CALL choose_eigv_solver(Atilde_fm, evects_Atilde_fm, evals(1:nkvs))

      ALLOCATE (evects_Atilde(nkvs, nrvs))
      CALL cp_fm_get_submatrix(evects_Atilde_fm, evects_Atilde, start_row=1, start_col=1, n_rows=nkvs, n_cols=nrvs)
      CALL cp_fm_release(evects_Atilde_fm)
      CALL cp_fm_release(Atilde_fm)

!$OMP PARALLEL DO DEFAULT(NONE), &
!$OMP             PRIVATE(act, ikv, irv, ispin), &
!$OMP             SHARED(Aop_krylov, Aop_ritz, krylov_vects, evects_Atilde, nkvs, nrvs, nspins, ritz_vects)
      DO irv = 1, nrvs
         DO ispin = 1, nspins
            CALL cp_fm_set_all(ritz_vects(ispin, irv)%matrix, 0.0_dp)
            CALL cp_fm_set_all(Aop_ritz(ispin, irv)%matrix, 0.0_dp)
         END DO

         DO ikv = 1, nkvs
            act = evects_Atilde(ikv, irv)
            DO ispin = 1, nspins
               CALL cp_fm_scale_and_add(1.0_dp, ritz_vects(ispin, irv)%matrix, &
                                        act, krylov_vects(ispin, ikv)%matrix)
               CALL cp_fm_scale_and_add(1.0_dp, Aop_ritz(ispin, irv)%matrix, &
                                        act, Aop_krylov(ispin, ikv)%matrix)
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO

      DEALLOCATE (evects_Atilde)

      CALL timestop(handle)

   END SUBROUTINE tddfpt_compute_ritz_vects

! **************************************************************************************************
!> \brief Expand Krylov space by computing residual vectors.
!> \param residual_vects          residual vectors (modified on exit)
!> \param evals                   Ritz eigenvalues (modified on exit)
!> \param ritz_vects              Ritz eigenvectors
!> \param Aop_ritz                approximate action of TDDFPT operator on Ritz vectors
!> \param gs_mos                  molecular orbitals optimised for the ground state
!> \param matrix_s                overlap matrix
!> \par History
!>    * 06.2016 created [Sergey Chulkov]
!>    * 03.2017 refactored to achieve significant performance gain [Sergey Chulkov]
! **************************************************************************************************
   SUBROUTINE tddfpt_compute_residual_vects(residual_vects, evals, ritz_vects, Aop_ritz, gs_mos, &
                                            matrix_s)
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: residual_vects
      REAL(kind=dp), DIMENSION(:), INTENT(in)            :: evals
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(in)    :: ritz_vects, Aop_ritz
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      TYPE(dbcsr_type), POINTER                          :: matrix_s

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_compute_residual_vects', &
         routineP = moduleN//':'//routineN
      REAL(kind=dp), PARAMETER :: eref_scale = 0.99_dp, threshold = 16.0_dp*EPSILON(1.0_dp)

      INTEGER                                            :: handle, icol_local, irow_local, irv, &
                                                            ispin, nao, ncols_local, nrows_local, &
                                                            nrvs, nspins
      INTEGER, DIMENSION(:), POINTER                     :: col_indices_local, row_indices_local
      INTEGER, DIMENSION(maxspins)                       :: nactive, nmo_virt
      REAL(kind=dp)                                      :: e_occ_plus_lambda, eref, lambda
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: weights_ldata
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: awork, vomat
      TYPE(cp_fm_struct_type), POINTER                   :: ao_mo_struct, virt_mo_struct

      CALL timeset(routineN, handle)

      nspins = SIZE(residual_vects, 1)
      nrvs = SIZE(residual_vects, 2)

      IF (nrvs > 0) THEN
         CALL dbcsr_get_info(matrix_s, nfullrows_total=nao)
         ALLOCATE (awork(nspins), vomat(nspins))
         DO ispin = 1, nspins
            nmo_virt(ispin) = SIZE(gs_mos(ispin)%evals_virt)
            !
            CALL cp_fm_get_info(matrix=ritz_vects(ispin, 1)%matrix, matrix_struct=ao_mo_struct, &
                                ncol_global=nactive(ispin))
            NULLIFY (awork(ispin)%matrix)
            CALL cp_fm_create(awork(ispin)%matrix, ao_mo_struct)
            !
            CALL cp_fm_struct_create(virt_mo_struct, nrow_global=nmo_virt(ispin), &
                                     ncol_global=nactive(ispin), template_fmstruct=ao_mo_struct)
            NULLIFY (vomat(ispin)%matrix)
            CALL cp_fm_create(vomat(ispin)%matrix, virt_mo_struct)
            CALL cp_fm_struct_release(virt_mo_struct)
         END DO

         ! *** actually compute residual vectors ***
         DO irv = 1, nrvs
            lambda = evals(irv)

            DO ispin = 1, nspins
               CALL cp_fm_get_info(vomat(ispin)%matrix, nrow_local=nrows_local, &
                                   ncol_local=ncols_local, row_indices=row_indices_local, &
                                   col_indices=col_indices_local, local_data=weights_ldata)

               ! awork := Ab(ispin, irv) - evals(irv) b(ispin, irv), where 'b' is a Ritz vector
               CALL cp_dbcsr_sm_fm_multiply(matrix_s, ritz_vects(ispin, irv)%matrix, awork(ispin)%matrix, &
                                            ncol=nactive(ispin), alpha=-lambda, beta=0.0_dp)
               CALL cp_fm_scale_and_add(1.0_dp, awork(ispin)%matrix, 1.0_dp, Aop_ritz(ispin, irv)%matrix)
               !
               CALL cp_gemm('T', 'N', nmo_virt(ispin), nactive(ispin), nao, 1.0_dp, gs_mos(ispin)%mos_virt, &
                            awork(ispin)%matrix, 0.0_dp, vomat(ispin)%matrix)

               DO icol_local = 1, ncols_local
                  e_occ_plus_lambda = gs_mos(ispin)%evals_occ(col_indices_local(icol_local))+lambda

                  DO irow_local = 1, nrows_local
                     eref = gs_mos(ispin)%evals_virt(row_indices_local(irow_local))-e_occ_plus_lambda

                     ! eref = e_virt - e_occ - lambda = e_virt - e_occ - (eref_scale*lambda + (1-eref_scale)*lambda);
                     ! eref_new = e_virt - e_occ - eref_scale*lambda = eref + (1 - eref_scale)*lambda
                     IF (ABS(eref) < threshold) &
                        eref = eref+(1.0_dp-eref_scale)*lambda

                     weights_ldata(irow_local, icol_local) = weights_ldata(irow_local, icol_local)/eref
                  END DO
               END DO

               CALL cp_gemm('N', 'N', nao, nactive(ispin), nmo_virt(ispin), 1.0_dp, gs_mos(ispin)%mos_virt, &
                            vomat(ispin)%matrix, 0.0_dp, residual_vects(ispin, irv)%matrix)
            END DO
         END DO
         !
         DO ispin = 1, nspins
            CALL cp_fm_release(awork(ispin)%matrix)
            CALL cp_fm_release(vomat(ispin)%matrix)
         END DO
         DEALLOCATE (awork, vomat)
      END IF

      CALL timestop(handle)

   END SUBROUTINE tddfpt_compute_residual_vects

! **************************************************************************************************
!> \brief Perform Davidson iterations.
!> \param evects                TDDFPT trial vectors (modified on exit)
!> \param evals                 TDDFPT eigenvalues (modified on exit)
!> \param S_evects              cached matrix product S * evects (modified on exit)
!> \param gs_mos                molecular orbitals optimised for the ground state
!> \param do_hfx                flag that activates computation of exact-exchange terms
!> \param tddfpt_control        TDDFPT control parameters
!> \param matrix_ks             Kohn-Sham matrix
!> \param qs_env                Quickstep environment
!> \param kernel_env            kernel environment
!> \param sub_env               parallel (sub)group environment
!> \param logger                CP2K logger
!> \param iter_unit             I/O unit to write basic iteration information
!> \param energy_unit           I/O unit to write detailed energy information
!> \param tddfpt_print_section  TDDFPT print input section (need to write TDDFPT restart files)
!> \param work_matrices         collection of work matrices (modified on exit)
!> \return energy convergence achieved (in Hartree)
!> \par History
!>    * 03.2017 code related to Davidson eigensolver has been moved here from the main subroutine
!>              tddfpt() [Sergey Chulkov]
!> \note Based on the subroutines apply_op() and iterative_solver() originally created by
!>       Thomas Chassaing in 2002.
! **************************************************************************************************
   FUNCTION tddfpt_davidson_solver(evects, evals, S_evects, gs_mos, do_hfx, tddfpt_control, &
                                   matrix_ks, qs_env, kernel_env, &
                                   sub_env, logger, iter_unit, energy_unit, &
                                   tddfpt_print_section, work_matrices) RESULT(conv)
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(inout) :: evects
      REAL(kind=dp), DIMENSION(:), INTENT(inout)         :: evals
      TYPE(cp_fm_p_type), DIMENSION(:, :), INTENT(inout) :: S_evects
      TYPE(tddfpt_ground_state_mos), DIMENSION(:), &
         INTENT(in)                                      :: gs_mos
      LOGICAL, INTENT(in)                                :: do_hfx
      TYPE(tddfpt2_control_type), POINTER                :: tddfpt_control
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddfpt_kernel_env_type), INTENT(in)           :: kernel_env
      TYPE(tddfpt_subgroup_env_type), INTENT(in)         :: sub_env
      TYPE(cp_logger_type), POINTER                      :: logger
      INTEGER, INTENT(in)                                :: iter_unit, energy_unit
      TYPE(section_vals_type), POINTER                   :: tddfpt_print_section
      TYPE(tddfpt_work_matrices), INTENT(inout)          :: work_matrices
      REAL(kind=dp)                                      :: conv

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_davidson_solver', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, istate, iter, &
                                                            max_krylov_vects, nspins, nstates, &
                                                            nstates_conv, nvects_exist, nvects_new
      INTEGER(kind=int_8)                                :: nstates_total
      LOGICAL                                            :: is_nonortho
      REAL(kind=dp)                                      :: t1, t2
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: evals_last
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: Atilde
      TYPE(cp_fm_p_type), ALLOCATABLE, DIMENSION(:, :)   :: Aop_krylov, Aop_ritz, krylov_vects, &
                                                            S_krylov
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s

      CALL timeset(routineN, handle)

      nspins = SIZE(gs_mos)
      nstates = tddfpt_control%nstates
      nstates_total = tddfpt_total_number_of_states(gs_mos)

      IF (debug_this_module) THEN
         CPASSERT(SIZE(evects, 1) == nspins)
         CPASSERT(SIZE(evects, 2) == nstates)
         CPASSERT(SIZE(evals) == nstates)
      END IF

      CALL get_qs_env(qs_env, matrix_s=matrix_s)

      ! adjust the number of Krylov vectors
      max_krylov_vects = tddfpt_control%nkvs
      IF (max_krylov_vects < nstates) max_krylov_vects = nstates
      IF (INT(max_krylov_vects, kind=int_8) > nstates_total) max_krylov_vects = INT(nstates_total)

      ALLOCATE (Aop_ritz(nspins, nstates))
      DO istate = 1, nstates
         DO ispin = 1, nspins
            NULLIFY (Aop_ritz(ispin, istate)%matrix)
            CALL fm_pool_create_fm(work_matrices%fm_pool_ao_mo_occ(ispin)%pool, Aop_ritz(ispin, istate)%matrix)
         END DO
      END DO

      ALLOCATE (evals_last(max_krylov_vects))
      ALLOCATE (Aop_krylov(nspins, max_krylov_vects), krylov_vects(nspins, max_krylov_vects), &
                S_krylov(nspins, max_krylov_vects))
      DO istate = 1, max_krylov_vects
         DO ispin = 1, nspins
            NULLIFY (Aop_krylov(ispin, istate)%matrix, krylov_vects(ispin, istate)%matrix, &
                     S_krylov(ispin, istate)%matrix)
         END DO
      END DO

      DO istate = 1, nstates
         DO ispin = 1, nspins
            CALL fm_pool_create_fm(work_matrices%fm_pool_ao_mo_occ(ispin)%pool, krylov_vects(ispin, istate)%matrix)
            CALL cp_fm_to_fm(evects(ispin, istate)%matrix, krylov_vects(ispin, istate)%matrix)

            CALL fm_pool_create_fm(work_matrices%fm_pool_ao_mo_occ(ispin)%pool, S_krylov(ispin, istate)%matrix)
            CALL cp_fm_to_fm(S_evects(ispin, istate)%matrix, S_krylov(ispin, istate)%matrix)

            CALL fm_pool_create_fm(work_matrices%fm_pool_ao_mo_occ(ispin)%pool, Aop_krylov(ispin, istate)%matrix)
         END DO
      END DO

      nvects_exist = 0
      nvects_new = nstates

      t1 = m_walltime()

      ALLOCATE (Atilde(1, 1))

      DO
         ! davidson iteration
         CALL cp_iterate(logger%iter_info, iter_nr_out=iter)

         CALL tddfpt_compute_Aop_evects(Aop_evects=Aop_krylov(:, nvects_exist+1:nvects_exist+nvects_new), &
                                        evects=krylov_vects(:, nvects_exist+1:nvects_exist+nvects_new), &
                                        S_evects=S_krylov(:, nvects_exist+1:nvects_exist+nvects_new), &
                                        gs_mos=gs_mos, tddfpt_control=tddfpt_control, &
                                        do_hfx=do_hfx, matrix_ks=matrix_ks, &
                                        qs_env=qs_env, kernel_env=kernel_env, &
                                        sub_env=sub_env, &
                                        work_matrices=work_matrices)

         CALL tddfpt_compute_ritz_vects(ritz_vects=evects, Aop_ritz=Aop_ritz, &
                                        evals=evals_last(1:nvects_exist+nvects_new), &
                                        krylov_vects=krylov_vects(:, 1:nvects_exist+nvects_new), &
                                        Aop_krylov=Aop_krylov(:, 1:nvects_exist+nvects_new), &
                                        Atilde=Atilde, nkvo=nvects_exist, nkvn=nvects_new)

         CALL tddfpt_write_restart(evects=evects, evals=evals_last(1:nstates), gs_mos=gs_mos, &
                                   logger=logger, tddfpt_print_section=tddfpt_print_section)

         conv = MAXVAL(ABS(evals_last(1:nstates)-evals(1:nstates)))

         nvects_exist = nvects_exist+nvects_new
         IF (nvects_exist+nvects_new > max_krylov_vects) &
            nvects_new = max_krylov_vects-nvects_exist
         IF (iter >= tddfpt_control%niters) nvects_new = 0

         IF (conv > tddfpt_control%conv .AND. nvects_new > 0) THEN
            ! compute residual vectors for the next iteration
            DO istate = 1, nvects_new
               DO ispin = 1, nspins
                  NULLIFY (Aop_krylov(ispin, nvects_exist+istate)%matrix, &
                           krylov_vects(ispin, nvects_exist+istate)%matrix, &
                           S_krylov(ispin, nvects_exist+istate)%matrix)
                  CALL fm_pool_create_fm(work_matrices%fm_pool_ao_mo_occ(ispin)%pool, &
                                         krylov_vects(ispin, nvects_exist+istate)%matrix)
                  CALL fm_pool_create_fm(work_matrices%fm_pool_ao_mo_occ(ispin)%pool, &
                                         S_krylov(ispin, nvects_exist+istate)%matrix)
                  CALL fm_pool_create_fm(work_matrices%fm_pool_ao_mo_occ(ispin)%pool, &
                                         Aop_krylov(ispin, nvects_exist+istate)%matrix)
               END DO
            END DO

            CALL tddfpt_compute_residual_vects(residual_vects=krylov_vects(:, nvects_exist+1:nvects_exist+nvects_new), &
                                               evals=evals_last(1:nvects_new), &
                                               ritz_vects=evects(:, 1:nvects_new), Aop_ritz=Aop_ritz(:, 1:nvects_new), &
                                               gs_mos=gs_mos, matrix_s=matrix_s(1)%matrix)

            CALL tddfpt_orthogonalize_psi1_psi0(krylov_vects(:, nvects_exist+1:nvects_exist+nvects_new), &
                                                work_matrices%S_C0_C0T)

            CALL tddfpt_orthonormalize_psi1_psi1(krylov_vects(:, 1:nvects_exist+nvects_new), nvects_new, &
                                                 S_krylov(:, 1:nvects_exist+nvects_new), matrix_s(1)%matrix)

            is_nonortho = tddfpt_is_nonorthogonal_psi1_psi0(krylov_vects(:, nvects_exist+1:nvects_exist+nvects_new), &
                                                            work_matrices%S_C0, tddfpt_control%orthogonal_eps)
         ELSE
            ! convergence or the maximum number of Krylov vectors have been achieved
            nvects_new = 0
            is_nonortho = .FALSE.
         END IF

         t2 = m_walltime()
         IF (energy_unit > 0) THEN
            WRITE (energy_unit, '(/,4X,A,T14,A,T36,A)') "State", "Exc. energy (eV)", "Convergence (eV)"
            DO istate = 1, nstates
               WRITE (energy_unit, '(1X,I8,T12,F14.7,T38,ES11.4)') istate, &
                  evals_last(istate)*evolt, (evals_last(istate)-evals(istate))*evolt
            END DO
            WRITE (energy_unit, *)
            CALL m_flush(energy_unit)
         END IF

         IF (iter_unit > 0) THEN
            nstates_conv = 0
            DO istate = 1, nstates
               IF (ABS(evals_last(istate)-evals(istate)) <= tddfpt_control%conv) &
                  nstates_conv = nstates_conv+1
            END DO

            WRITE (iter_unit, '(T7,I8,T24,F7.1,T40,ES11.4,T66,I8)') iter, t2-t1, conv, nstates_conv
            CALL m_flush(iter_unit)
         END IF

         t1 = t2
         evals(1:nstates) = evals_last(1:nstates)

         ! nvects_new == 0 if iter >= tddfpt_control%niters
         IF (nvects_new == 0 .OR. is_nonortho) THEN
            ! restart Davidson iterations
            CALL tddfpt_orthogonalize_psi1_psi0(evects, work_matrices%S_C0_C0T)
            CALL tddfpt_orthonormalize_psi1_psi1(evects, nstates, S_evects, matrix_s(1)%matrix)

            EXIT
         END IF
      END DO

      DEALLOCATE (Atilde)

      DO istate = nvects_exist+nvects_new, 1, -1
         DO ispin = nspins, 1, -1
            CALL fm_pool_give_back_fm(work_matrices%fm_pool_ao_mo_occ(ispin)%pool, Aop_krylov(ispin, istate)%matrix)
            CALL fm_pool_give_back_fm(work_matrices%fm_pool_ao_mo_occ(ispin)%pool, S_krylov(ispin, istate)%matrix)
            CALL fm_pool_give_back_fm(work_matrices%fm_pool_ao_mo_occ(ispin)%pool, krylov_vects(ispin, istate)%matrix)
         END DO
      END DO
      DEALLOCATE (Aop_krylov, krylov_vects, S_krylov)
      DEALLOCATE (evals_last)

      DO istate = nstates, 1, -1
         DO ispin = nspins, 1, -1
            CALL fm_pool_give_back_fm(work_matrices%fm_pool_ao_mo_occ(ispin)%pool, Aop_ritz(ispin, istate)%matrix)
         END DO
      END DO
      DEALLOCATE (Aop_ritz)

      CALL timestop(handle)

   END FUNCTION tddfpt_davidson_solver

END MODULE qs_tddfpt2_eigensolver
