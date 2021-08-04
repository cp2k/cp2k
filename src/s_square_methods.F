!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
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
   USE cp_control_types,                ONLY: s2_restraint_type
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_p_type
   USE input_constants,                 ONLY: do_s2_constraint,&
                                              do_s2_restraint
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_sum
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              has_uniform_occupation,&
                                              mo_set_p_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   ! Global parameters

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 's_square_methods'

   PUBLIC :: compute_s_square, s2_restraint

CONTAINS

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
   SUBROUTINE compute_s_square(mos, matrix_s, s_square, s_square_ideal, mo_derivs, strength)

      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      REAL(KIND=dp), INTENT(OUT)                         :: s_square, s_square_ideal
      TYPE(cp_fm_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: mo_derivs
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: strength

      CHARACTER(len=*), PARAMETER                        :: routineN = 'compute_s_square'

      INTEGER                                            :: handle, i, j, nalpha, nao, nao_beta, &
                                                            nbeta, ncol_local, nmo_alpha, &
                                                            nmo_beta, nrow_local
      LOGICAL                                            :: has_uniform_occupation_alpha, &
                                                            has_uniform_occupation_beta
      REAL(KIND=dp)                                      :: s2
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: local_data
      TYPE(cp_blacs_env_type), POINTER                   :: context
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_tmp
      TYPE(cp_fm_type), POINTER                          :: c_alpha, c_beta, catscb, sca, scb
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      NULLIFY (catscb)
      NULLIFY (c_alpha)
      NULLIFY (c_beta)
      NULLIFY (context)
      NULLIFY (fm_struct_tmp)
      NULLIFY (local_data)
      NULLIFY (para_env)
      NULLIFY (sca)
      NULLIFY (scb)

      SELECT CASE (SIZE(mos))
      CASE (1)
         ! Spin restricted case, i.e. nothing to do
         s_square = 0.0_dp
         s_square_ideal = 0.0_dp
         ! Restraining or constraining (\cal S) does not make sense
         CPASSERT(PRESENT(mo_derivs) .OR. PRESENT(strength))
      CASE (2)
         CALL get_mo_set(mo_set=mos(1)%mo_set, mo_coeff=c_alpha, homo=nalpha, nmo=nmo_alpha, nao=nao)
         CALL get_mo_set(mo_set=mos(2)%mo_set, mo_coeff=c_beta, homo=nbeta, nmo=nmo_beta, nao=nao_beta)
         CPASSERT(nao == nao_beta)
         has_uniform_occupation_alpha = has_uniform_occupation(mo_set=mos(1)%mo_set, last_mo=nalpha)
         has_uniform_occupation_beta = has_uniform_occupation(mo_set=mos(2)%mo_set, last_mo=nbeta)
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
            CALL cp_fm_create(scb, fm_struct_tmp, name="S*C(beta)")
            ! Compute S*C(beta)
            CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, c_beta, scb, nbeta)
            ! Compute C(alpha)^T*S*C(beta)
            CALL cp_gemm('T', 'N', nalpha, nbeta, nao, 1.0_dp, c_alpha, scb, 0.0_dp, catscb)
            ! Eq. 2.271 in Modern Quantum Chemistry by A. Szabo and N. S. Ostlund
            CALL cp_fm_get_info(catscb, local_data=local_data, nrow_local=nrow_local, ncol_local=ncol_local)
            s2 = 0.0_dp
            DO j = 1, ncol_local
               DO i = 1, nrow_local
                  s2 = s2 + local_data(i, j)**2
               END DO
            END DO
            CALL mp_sum(s2, para_env%group)
            s_square = s_square_ideal + nbeta - s2
            ! Only needed for constraining or restraining (\cal S)
            IF (PRESENT(mo_derivs)) THEN
               CPASSERT(SIZE(mo_derivs, 1) == 2)
               ! This gets really wrong for fractional occupations
               CALL get_mo_set(mo_set=mos(1)%mo_set, uniform_occupation=has_uniform_occupation_alpha)
               CPASSERT(has_uniform_occupation_alpha)
               CALL get_mo_set(mo_set=mos(2)%mo_set, uniform_occupation=has_uniform_occupation_beta)
               CPASSERT(has_uniform_occupation_beta)
               ! Add -strength*S*C(beta)*(C(alpha)^T*S*C(beta))^T to the alpha MO derivatives
               CALL cp_gemm('N', 'T', nao, nalpha, nbeta, -strength, scb, catscb, 1.0_dp, mo_derivs(1)%matrix)
               ! Create S*C(alpha)
               CALL cp_fm_get_info(c_alpha, matrix_struct=fm_struct_tmp)
               CALL cp_fm_create(sca, fm_struct_tmp, name="S*C(alpha)")
               ! Compute S*C(alpha)
               CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, c_alpha, sca, nalpha)
               ! Add -strength*S*C(alpha)*C(alpha)^T*S*C(beta) to the beta MO derivatives
               CALL cp_gemm('N', 'N', nao, nbeta, nalpha, -strength, sca, catscb, 1.0_dp, mo_derivs(2)%matrix)
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
         CPABORT("Alpha, beta, what else ?")
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE compute_s_square

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

      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: mo_derivs
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
            CALL compute_s_square(mos, matrix_s, s_square, s_square_ideal)
         ELSE
            CALL compute_s_square(mos, matrix_s, s_square, s_square_ideal, &
                                  mo_derivs, s2_restraint_control%strength)
         END IF
         energy = s2_restraint_control%strength*(s_square - s2_restraint_control%target)
         s2_restraint_control%s2_order_p = s_square
      CASE (do_s2_restraint) ! not yet implemented
         CPABORT("")
      CASE DEFAULT
         CPABORT("")
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE s2_restraint

END MODULE s_square_methods
