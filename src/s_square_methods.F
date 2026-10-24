!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Methods related to (\cal S)^2 (i.e. spin)
!> \par History
!>      03.2006 copied compute_s_square from qs_scf_post (Joost VandeVondele)
!>      08.2021 revised (Matthias Krack, MK)
!> \author Joost VandeVondele
! **************************************************************************************************
MODULE s_square_methods

   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_schur_product_cc
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_info,&
                                              cp_cfm_release,&
                                              cp_cfm_to_fm,&
                                              cp_cfm_type
   USE cp_control_types,                ONLY: s2_restraint_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply
   USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                              cp_fm_row_scale,&
                                              cp_fm_schur_product
   USE cp_fm_pool_types,                ONLY: cp_fm_pool_type,&
                                              fm_pool_create_fm,&
                                              fm_pool_give_back_fm
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_sumall,&
                                              cp_fm_type
   USE input_constants,                 ONLY: do_s2_constraint,&
                                              do_s2_restraint
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: get_kpoint_env,&
                                              get_kpoint_info,&
                                              kpoint_env_type,&
                                              kpoint_type
   USE mathconstants,                   ONLY: one,&
                                              z_one,&
                                              z_zero,&
                                              zero
   USE message_passing,                 ONLY: mp_para_env_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE qs_kpoint_operators,             ONLY: kpoint_operator_context_create,&
                                              kpoint_operator_context_release,&
                                              kpoint_operator_context_type,&
                                              kpoint_operator_finish,&
                                              kpoint_operator_prepare,&
                                              kpoint_operator_start,&
                                              kpoint_slot_s,&
                                              kpoint_spin_free
   USE qs_matrix_pools,                 ONLY: mpools_get,&
                                              qs_matrix_pools_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              has_uniform_occupation,&
                                              mo_set_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   ! Global parameters

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 's_square_methods'

   PUBLIC :: compute_s_square, compute_s_square_kp
   PUBLIC :: s2_restraint

CONTAINS

! **************************************************************************************************
!> \brief Compute the expectation value <(\cal S)^2> of the single determinant defined by the
!>        spin up (alpha) and spin down (beta) orbitals
!> \param mos [in] MO set with all MO information including the alpha and beta MO coefficients
!> \param matrix_s [in] AO overlap matrix S (do not mix with the spin operator (\cal S))
!> \param s_square [out] <(\cal S)^2> including potential spin contaminations
!> \param s_square_ideal [out] Ideal value for <(\cal S)^2> without any spin contaminations
!> \par History
!>      07.2004 created (Joost VandeVondele)
!>      08.2021 revised (Matthias Krack, MK)
!>      10.2026 adapted and generalized to partial occupation and k-points
!> \note
!>      See Eqs. 2.271 and 2.272 in Modern Quantum Chemistry by A. Szabo and N. S. Ostlund
! **************************************************************************************************
   SUBROUTINE compute_s_square(mos, matrix_s, s_square, s_square_ideal)

      TYPE(mo_set_type), DIMENSION(:), INTENT(IN)        :: mos
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      REAL(KIND=dp), INTENT(OUT)                         :: s_square, s_square_ideal

      CHARACTER(len=*), PARAMETER                        :: routineN = 'compute_s_square'

      INTEGER                                            :: handle, nalpha, nao, nao_beta, nbeta
      REAL(KIND=dp)                                      :: nela, nelb, s2
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:), POINTER   :: occa, occb
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: asb, asb2, scb
      TYPE(cp_fm_type), POINTER                          :: c_alpha, c_beta

      CALL timeset(routineN, handle)

      CALL get_mo_set(mo_set=mos(1), mo_coeff=c_alpha, nmo=nalpha, nao=nao)
      CALL get_mo_set(mo_set=mos(2), mo_coeff=c_beta, nmo=nbeta, nao=nao_beta)
      CPASSERT(nao == nao_beta)
      CALL get_mo_set(mo_set=mos(1), n_el_f=nela, occupation_numbers=occa)
      CALL get_mo_set(mo_set=mos(2), n_el_f=nelb, occupation_numbers=occb)
      ! Eq. 2.272 in Modern Quantum Chemistry by A. Szabo and N. S. Ostlund
      s_square_ideal = REAL((nela - nelb)*(nela - nelb + 2), KIND=dp)/4.0_dp
      ! Create MO overlap matrix
      CALL cp_fm_get_info(c_beta, matrix_struct=fm_struct)
      ! Prepare C(alpha)^T*S*C(beta)
      CALL cp_fm_create(asb, fm_struct, nrow=nalpha, ncol=nbeta)
      CALL cp_fm_create(asb2, fm_struct, nrow=nalpha, ncol=nbeta)
      CALL cp_fm_create(scb, fm_struct)
      ! Compute S*C(beta)
      CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, c_beta, scb, nbeta)
      ! Compute C(alpha)^T*S*C(beta)
      CALL parallel_gemm('T', 'N', nalpha, nbeta, nao, one, c_alpha, scb, zero, asb)
      ! square the matrix elements
      CALL cp_fm_schur_product(asb, asb, asb2)
      ! occupation number weighting
      CALL cp_fm_row_scale(asb2, occa)
      CALL cp_fm_column_scale(asb2, occb)
      s2 = cp_fm_sumall(asb2)
      s_square = s_square_ideal + nelb - s2
      CALL cp_fm_release(scb)
      CALL cp_fm_release(asb)
      CALL cp_fm_release(asb2)

      CALL timestop(handle)

   END SUBROUTINE compute_s_square

! **************************************************************************************************
!> \brief Compute the expectation value <(\cal S)^2> of the single determinant defined by the
!>        spin up (alpha) and spin down (beta) orbitals
!> \param mosin ...
!> \param kpoints ...
!> \param matrix_s ...
!> \param s_square [out] <(\cal S)^2> including potential spin contaminations
!> \param s_square_ideal [out] Ideal value for <(\cal S)^2> without any spin contaminations
!> \par History
!>      10.2026 adapted and generalized to partial occupation and k-points
! **************************************************************************************************
   SUBROUTINE compute_s_square_kp(mosin, kpoints, matrix_s, s_square, s_square_ideal)

      TYPE(mo_set_type), DIMENSION(:), INTENT(IN)        :: mosin
      TYPE(kpoint_type), INTENT(IN), POINTER             :: kpoints
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s
      REAL(KIND=dp), INTENT(OUT)                         :: s_square, s_square_ideal

      CHARACTER(len=*), PARAMETER :: routineN = 'compute_s_square_kp'

      INTEGER                                            :: handle, ikp, kplocal, kpmax, nalpha, &
                                                            nao, nbeta
      INTEGER, DIMENSION(2)                              :: kp_range
      INTEGER, DIMENSION(:, :), POINTER                  :: kp_dist
      LOGICAL                                            :: wfn_real
      REAL(KIND=dp)                                      :: nela, nelb, s2, wkp
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:), POINTER   :: occa, occb
      TYPE(cp_cfm_type)                                  :: casb, casb2, cscb, csmat
      TYPE(cp_cfm_type), POINTER                         :: c_alpha, c_beta
      TYPE(cp_fm_pool_type), POINTER                     :: ao_ao_fm_pools
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct, matrix_struct
      TYPE(cp_fm_type)                                   :: asb, asb2, fmglobal, fmlocal, rsmat, scb
      TYPE(cp_fm_type), POINTER                          :: moc, r_alpha, r_beta
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(kpoint_operator_context_type)                 :: op_ctx
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(qs_matrix_pools_type), POINTER                :: mpools

      CALL timeset(routineN, handle)

      CALL get_kpoint_info(kpoints, use_real_wfn=wfn_real, kp_range=kp_range, &
                           kp_dist=kp_dist, sab_nl=sab_nl)
      CPASSERT(ASSOCIATED(sab_nl))

      kplocal = kp_range(2) - kp_range(1) + 1
      kpmax = MAXVAL(kp_dist(2, :) - kp_dist(1, :) + 1)

      ! fm pools to be used within a kpoint group
      CALL get_kpoint_info(kpoints, mpools=mpools)
      CALL mpools_get(mpools, ao_ao_fm_pools=ao_ao_fm_pools)
      CALL fm_pool_create_fm(ao_ao_fm_pools, fmlocal)
      CALL cp_fm_get_info(fmlocal, matrix_struct=matrix_struct)
      IF (wfn_real) THEN
         CALL cp_fm_create(rsmat, matrix_struct)
      ELSE
         CALL cp_cfm_create(csmat, matrix_struct)
      END IF
      CALL fm_pool_give_back_fm(ao_ao_fm_pools, fmlocal)

      ! fm matrix_structure for the global environment
      CALL get_mo_set(mo_set=mosin(1), mo_coeff=moc, nao=nao)
      CALL cp_fm_create(fmglobal, moc%matrix_struct, nrow=nao, ncol=nao)
      ! The overlap is spin free: the routine drives every instance under the
      ! single spin-free key.
      CALL kpoint_operator_context_create(op_ctx, kpoints, fmglobal, 1)
      CALL kpoint_operator_prepare(op_ctx, matrix_s, kpoint_slot_s, matrix_row=kpoint_spin_free)
      DO ikp = 1, kpmax
         CALL kpoint_operator_start(op_ctx, ikp, kpoint_spin_free, matrix_s, kpoint_slot_s, &
                                    use_grid=.FALSE., matrix_row=kpoint_spin_free)
      END DO
      CALL cp_fm_release(fmglobal)

      s2 = 0.0_dp
      DO ikp = 1, kpmax
         IF (wfn_real) THEN
            CALL kpoint_operator_finish(op_ctx, ikp, kpoint_spin_free, kpoint_slot_s, fm_re=rsmat)
         ELSE
            CALL kpoint_operator_finish(op_ctx, ikp, kpoint_spin_free, kpoint_slot_s, cmat=csmat)
         END IF

         IF (ikp > kplocal) CYCLE

         kp => kpoints%kp_env(ikp)%kpoint_env
         CALL get_kpoint_env(kp, wkp=wkp, mos=mos)
         CALL get_mo_set(mo_set=mos(1), occupation_numbers=occa)
         CALL get_mo_set(mo_set=mos(2), occupation_numbers=occb)

         IF (wfn_real) THEN
            CALL get_mo_set(mo_set=mos(1), mo_coeff=r_alpha, nmo=nalpha)
            CALL get_mo_set(mo_set=mos(2), mo_coeff=r_beta, nmo=nbeta)
            CALL cp_fm_get_info(r_beta, matrix_struct=fm_struct)
            CALL cp_fm_create(scb, fm_struct)
            CALL cp_fm_create(asb, fm_struct, nrow=nalpha, ncol=nbeta)
            CALL cp_fm_create(asb2, fm_struct, nrow=nalpha, ncol=nbeta)
            CALL parallel_gemm('N', 'N', nao, nbeta, nao, one, rsmat, r_beta, zero, scb)
            CALL parallel_gemm('T', 'N', nalpha, nbeta, nao, one, r_alpha, scb, zero, asb)
            CALL cp_fm_schur_product(asb, asb, asb2)
            CALL cp_fm_release(scb)
            CALL cp_fm_release(asb)
         ELSE
            CALL get_mo_set(mo_set=mos(1), cmo_coeff=c_alpha, nmo=nalpha)
            CALL get_mo_set(mo_set=mos(2), cmo_coeff=c_beta, nmo=nbeta)
            CALL cp_cfm_get_info(c_beta, matrix_struct=fm_struct)
            CALL cp_cfm_create(cscb, fm_struct)
            CALL cp_cfm_create(casb, fm_struct, nrow=nalpha, ncol=nbeta)
            CALL cp_cfm_create(casb2, fm_struct, nrow=nalpha, ncol=nbeta)
            CALL cp_fm_create(asb2, fm_struct, nrow=nalpha, ncol=nbeta)
            CALL parallel_gemm('N', 'N', nao, nbeta, nao, z_one, csmat, c_beta, z_zero, cscb)
            CALL parallel_gemm('C', 'N', nalpha, nbeta, nao, z_one, c_alpha, cscb, z_zero, casb)
            CALL cp_cfm_schur_product_cc(casb, casb, casb2)
            CALL cp_cfm_to_fm(casb2, asb2)
            CALL cp_cfm_release(cscb)
            CALL cp_cfm_release(casb)
            CALL cp_cfm_release(casb2)
         END IF
         CALL cp_fm_row_scale(asb2, occa)
         CALL cp_fm_column_scale(asb2, occb)
         s2 = s2 + wkp*cp_fm_sumall(asb2)
         CALL cp_fm_release(asb2)

      END DO

      !sum s2 over kp group
      CALL kpoints%para_env_inter_kp%sum(s2)

      ikp = 1
      IF (ikp <= kplocal) THEN
         kp => kpoints%kp_env(ikp)%kpoint_env
         CALL get_kpoint_env(kp, mos=mos)
         CALL get_mo_set(mo_set=mos(1), occupation_numbers=occa)
         CALL get_mo_set(mo_set=mos(2), occupation_numbers=occb)
         nela = SUM(occa)
         nelb = SUM(occb)
      ELSE
         nela = 0
         nelb = 0
      END IF
      !
      CALL kpoints%para_env_kp%max(nela)
      CALL kpoints%para_env_kp%max(nelb)
      !
      s_square_ideal = REAL((nela - nelb)*(nela - nelb + 2), KIND=dp)/4.0_dp
      s_square = s_square_ideal + nelb - s2

      CALL kpoint_operator_context_release(op_ctx)

      IF (wfn_real) THEN
         CALL cp_fm_release(rsmat)
      ELSE
         CALL cp_cfm_release(csmat)
      END IF

      CALL timestop(handle)

   END SUBROUTINE compute_s_square_kp

! **************************************************************************************************
!> \brief restrains/constrains the value of s2 in a calculation
!> \param mos input
!> \param matrix_s input
!> \param mo_derivs inout if present, add the derivative of s_square wrt mos to mo_derivs
!> \param energy ...
!> \param s2_restraint_control ...
!> \param just_energy ...
!> \par History
!>      07.2004 created [ Joost VandeVondele ]
! **************************************************************************************************
   SUBROUTINE s2_restraint(mos, matrix_s, mo_derivs, energy, &
                           s2_restraint_control, just_energy)

      TYPE(mo_set_type), DIMENSION(:), INTENT(IN)        :: mos
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(cp_fm_type), DIMENSION(:), INTENT(INOUT)      :: mo_derivs
      REAL(kind=dp)                                      :: energy
      TYPE(s2_restraint_type), POINTER                   :: s2_restraint_control
      LOGICAL                                            :: just_energy

      CHARACTER(len=*), PARAMETER                        :: routineN = 's2_restraint'

      INTEGER                                            :: handle
      REAL(kind=dp)                                      :: s_square, s_square_ideal

      CALL timeset(routineN, handle)

      SELECT CASE (s2_restraint_control%functional_form)
      CASE (do_s2_constraint)
         IF (just_energy) THEN
            CALL compute_s2_energyforce(mos, matrix_s, s_square, s_square_ideal)
         ELSE
            CALL compute_s2_energyforce(mos, matrix_s, s_square, s_square_ideal, &
                                        mo_derivs, s2_restraint_control%strength)
         END IF
         energy = s2_restraint_control%strength*(s_square - s2_restraint_control%target)
         s2_restraint_control%s2_order_p = s_square
      CASE (do_s2_restraint) ! not yet implemented
         CPABORT("do_s2_restraint for <S**2> not yet implemented")
      CASE DEFAULT
         CPABORT("Unknown <S**2> constraint/restraint form")
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE s2_restraint

! **************************************************************************************************
!> \brief Compute the expectation value <(\cal S)^2> of the single determinant defined by the
!>        spin up (alpha) and spin down (beta) orbitals
!> \param mos [in] MO set with all MO information including the alpha and beta MO coefficients
!> \param matrix_s [in] AO overlap matrix S (do not mix with the spin operator (\cal S))
!> \param s_square [out] <(\cal S)^2> including potential spin contaminations
!> \param s_square_ideal [out] Ideal value for <(\cal S)^2> without any spin contaminations
!> \param mo_derivs [inout] If present, add the derivative of s_square wrt the MOs to mo_derivs
!> \param strength [in] Strength for constraining or restraining (\cal S)^2
!> \par History
!>      07.2004 created (Joost VandeVondele)
!>      08.2021 revised (Matthias Krack, MK)
!> \note
!>      See Eqs. 2.271 and 2.272 in Modern Quantum Chemistry by A. Szabo and N. S. Ostlund
! **************************************************************************************************
   SUBROUTINE compute_s2_energyforce(mos, matrix_s, s_square, s_square_ideal, mo_derivs, strength)

      TYPE(mo_set_type), DIMENSION(:), INTENT(IN)        :: mos
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      REAL(KIND=dp), INTENT(OUT)                         :: s_square, s_square_ideal
      TYPE(cp_fm_type), DIMENSION(:), INTENT(INOUT), &
         OPTIONAL                                        :: mo_derivs
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: strength

      CHARACTER(len=*), PARAMETER :: routineN = 'compute_s2_energyforce'

      INTEGER                                            :: handle, i, j, nalpha, nao, nao_beta, &
                                                            nbeta, ncol_local, nmo_alpha, &
                                                            nmo_beta, nrow_local
      LOGICAL                                            :: has_uniform_occupation_alpha, &
                                                            has_uniform_occupation_beta
      REAL(KIND=dp)                                      :: s2
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         POINTER                                         :: local_data
      TYPE(cp_blacs_env_type), POINTER                   :: context
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_tmp
      TYPE(cp_fm_type)                                   :: catscb, sca, scb
      TYPE(cp_fm_type), POINTER                          :: c_alpha, c_beta
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      NULLIFY (context)
      NULLIFY (fm_struct_tmp)
      NULLIFY (local_data)
      NULLIFY (para_env)

      SELECT CASE (SIZE(mos))
      CASE (1)
         ! Spin restricted case, i.e. nothing to do
         s_square = 0.0_dp
         s_square_ideal = 0.0_dp
         ! Restraining or constraining (\cal S) does not make sense
         CPASSERT(PRESENT(mo_derivs) .OR. PRESENT(strength))
      CASE (2)
         CALL get_mo_set(mo_set=mos(1), mo_coeff=c_alpha, homo=nalpha, nmo=nmo_alpha, nao=nao)
         CALL get_mo_set(mo_set=mos(2), mo_coeff=c_beta, homo=nbeta, nmo=nmo_beta, nao=nao_beta)
         CPASSERT(nao == nao_beta)
         has_uniform_occupation_alpha = has_uniform_occupation(mo_set=mos(1), last_mo=nalpha)
         has_uniform_occupation_beta = has_uniform_occupation(mo_set=mos(2), last_mo=nbeta)
         ! This makes only sense if we have uniform occupations for the alpha and beta spin MOs while
         ! ignoring potentially added MOs with zero occupation
         IF (has_uniform_occupation_alpha .AND. has_uniform_occupation_beta) THEN
            ! Eq. 2.272 in Modern Quantum Chemistry by A. Szabo and N. S. Ostlund
            s_square_ideal = REAL((nalpha - nbeta)*(nalpha - nbeta + 2), KIND=dp)/4.0_dp
            ! Create overlap matrix
            CALL cp_fm_get_info(c_alpha, para_env=para_env, context=context)
            CALL cp_fm_struct_create(fm_struct_tmp, para_env=para_env, context=context, &
                                     nrow_global=nalpha, ncol_global=nbeta)
            ! Prepare C(alpha)^T*S*C(beta)
            CALL cp_fm_create(catscb, fm_struct_tmp, name="C(alpha)^T*S*C(beta)")
            CALL cp_fm_struct_release(fm_struct_tmp)
            ! Create S*C(beta)
            CALL cp_fm_get_info(c_beta, matrix_struct=fm_struct_tmp)
            CALL cp_fm_create(scb, fm_struct_tmp, name="S*C(beta)", set_zero=.TRUE.)
            ! Compute S*C(beta)
            CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, c_beta, scb, nbeta)
            ! Compute C(alpha)^T*S*C(beta)
            CALL parallel_gemm('T', 'N', nalpha, nbeta, nao, 1.0_dp, c_alpha, scb, 0.0_dp, catscb)
            ! Eq. 2.271 in Modern Quantum Chemistry by A. Szabo and N. S. Ostlund
            CALL cp_fm_get_info(catscb, local_data=local_data, nrow_local=nrow_local, ncol_local=ncol_local)
            s2 = 0.0_dp
            DO j = 1, ncol_local
               DO i = 1, nrow_local
                  s2 = s2 + local_data(i, j)**2
               END DO
            END DO
            CALL para_env%sum(s2)
            s_square = s_square_ideal + nbeta - s2
            ! Only needed for constraining or restraining (\cal S)
            IF (PRESENT(mo_derivs)) THEN
               CPASSERT(SIZE(mo_derivs, 1) == 2)
               ! This gets really wrong for fractional occupations
               CALL get_mo_set(mo_set=mos(1), uniform_occupation=has_uniform_occupation_alpha)
               CPASSERT(has_uniform_occupation_alpha)
               CALL get_mo_set(mo_set=mos(2), uniform_occupation=has_uniform_occupation_beta)
               CPASSERT(has_uniform_occupation_beta)
               ! Add -strength*S*C(beta)*(C(alpha)^T*S*C(beta))^T to the alpha MO derivatives
               CALL parallel_gemm('N', 'T', nao, nalpha, nbeta, -strength, scb, catscb, 1.0_dp, mo_derivs(1))
               ! Create S*C(alpha)
               CALL cp_fm_get_info(c_alpha, matrix_struct=fm_struct_tmp)
               CALL cp_fm_create(sca, fm_struct_tmp, name="S*C(alpha)", set_zero=.TRUE.)
               ! Compute S*C(alpha)
               CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, c_alpha, sca, nalpha)
               ! Add -strength*S*C(alpha)*C(alpha)^T*S*C(beta) to the beta MO derivatives
               CALL parallel_gemm('N', 'N', nao, nbeta, nalpha, -strength, sca, catscb, 1.0_dp, mo_derivs(2))
               CALL cp_fm_release(sca)
            END IF
            CALL cp_fm_release(scb)
            CALL cp_fm_release(catscb)
         ELSE
            IF (.NOT. has_uniform_occupation_alpha) THEN
               CPHINT("The alpha orbitals (MOs) have a non-uniform occupation")
            END IF
            IF (.NOT. has_uniform_occupation_beta) THEN
               CPHINT("The beta orbitals (MOs) have a non-uniform occupation")
            END IF
            CPHINT("Skipping S**2 calculation")
         END IF
      CASE DEFAULT
         ! We should never get here
         CPABORT("Alpha, beta, what else?")
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE compute_s2_energyforce

END MODULE s_square_methods
