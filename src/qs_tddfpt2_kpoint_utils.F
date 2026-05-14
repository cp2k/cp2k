!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE qs_tddfpt2_kpoint_utils
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_create,&
                                              dbcsr_deallocate_matrix,&
                                              dbcsr_get_info,&
                                              dbcsr_p_type,&
                                              dbcsr_set,&
                                              dbcsr_type,&
                                              dbcsr_type_antisymmetric,&
                                              dbcsr_type_symmetric
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_get_submatrix,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_set_submatrix,&
                                              cp_fm_type
   USE kinds,                           ONLY: dp
   USE kpoint_methods,                  ONLY: kpoint_density_transform,&
                                              rskp_transform
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_env_p_type,&
                                              kpoint_env_type,&
                                              kpoint_type
   USE message_passing,                 ONLY: mp_para_env_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_rho_methods,                  ONLY: qs_rho_rebuild,&
                                              qs_rho_update_rho
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddfpt2_kpoint_utils'
   INTEGER, PARAMETER, PRIVATE          :: maxspins = 2

   TYPE, PUBLIC :: tddfpt_kpoint_transition_type
      INTEGER       :: ikp = 0
      INTEGER       :: ispin = 0
      INTEGER       :: iocc = 0
      INTEGER       :: ivirt = 0
      REAL(kind=dp) :: weight = 0.0_dp
      REAL(kind=dp) :: energy = 0.0_dp
   END TYPE tddfpt_kpoint_transition_type

   TYPE, PUBLIC :: tddfpt_kpoint_response_type
      INTEGER                                      :: ntrans = 0
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)     :: coeff_re
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)     :: coeff_im
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)     :: sigma_re
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)     :: sigma_im
   END TYPE tddfpt_kpoint_response_type

   PUBLIC :: tddfpt_kpoint_apply_energy_diff, tddfpt_kpoint_create_matrix_set, &
             tddfpt_kpoint_hmat_from_rs, tddfpt_kpoint_hmat_to_response, &
             tddfpt_kpoint_init_transitions, tddfpt_kpoint_release_matrix_set, &
             tddfpt_kpoint_response_add_scaled_sigma, &
             tddfpt_kpoint_response_copy_scaled_coeff, tddfpt_kpoint_response_create, &
             tddfpt_kpoint_response_dot, tddfpt_kpoint_response_norm, &
             tddfpt_kpoint_response_normalize, tddfpt_kpoint_response_project_out, &
             tddfpt_kpoint_response_release, tddfpt_kpoint_response_set_unit, &
             tddfpt_kpoint_response_to_pmat, tddfpt_kpoint_response_update_rho, &
             tddfpt_kpoint_response_zero, &
             tddfpt_kpoint_transition_offset

CONTAINS

! **************************************************************************************************
!> \brief Allocate a local k-point AO dense matrix set with the same distribution as kp%pmat.
!> \param qs_env ...
!> \param matrix_set ...
! **************************************************************************************************
   SUBROUTINE tddfpt_kpoint_create_matrix_set(qs_env, matrix_set)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), ALLOCATABLE, &
         DIMENSION(:, :, :), INTENT(out)                 :: matrix_set

      CHARACTER(len=*), PARAMETER :: routineN = 'tddfpt_kpoint_create_matrix_set'

      INTEGER                                            :: handle, ic, ikp_local, ispin, nkp_local, &
                                                            nspins
      INTEGER, DIMENSION(2)                              :: kp_range
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_env_p_type), DIMENSION(:), POINTER     :: kp_env
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(kpoint_type), POINTER                         :: kpoints

      CALL timeset(routineN, handle)

      NULLIFY (dft_control, kp, kp_env, kpoints, matrix_struct)
      CALL get_qs_env(qs_env, dft_control=dft_control, kpoints=kpoints)
      CPASSERT(ASSOCIATED(dft_control))
      CPASSERT(ASSOCIATED(kpoints))
      CALL get_kpoint_info(kpoints, kp_range=kp_range, kp_env=kp_env)
      CPASSERT(ASSOCIATED(kp_env))

      nspins = dft_control%nspins
      nkp_local = MAX(0, kp_range(2) - kp_range(1) + 1)
      ALLOCATE (matrix_set(nkp_local, 2, nspins))

      DO ikp_local = 1, nkp_local
         kp => kp_env(ikp_local)%kpoint_env
         CPASSERT(ASSOCIATED(kp))
         DO ic = 1, 2
            DO ispin = 1, nspins
               CALL cp_fm_get_info(kp%pmat(ic, ispin), matrix_struct=matrix_struct)
               CALL cp_fm_create(matrix_set(ikp_local, ic, ispin), matrix_struct)
               CALL cp_fm_set_all(matrix_set(ikp_local, ic, ispin), 0.0_dp)
            END DO
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE tddfpt_kpoint_create_matrix_set

! **************************************************************************************************
!> \brief Release a local k-point AO dense matrix set.
!> \param matrix_set ...
! **************************************************************************************************
   SUBROUTINE tddfpt_kpoint_release_matrix_set(matrix_set)
      TYPE(cp_fm_type), ALLOCATABLE, &
         DIMENSION(:, :, :), INTENT(inout)               :: matrix_set

      INTEGER                                            :: ic, ikp_local, ispin

      IF (ALLOCATED(matrix_set)) THEN
         DO ikp_local = 1, SIZE(matrix_set, 1)
            DO ic = 1, SIZE(matrix_set, 2)
               DO ispin = 1, SIZE(matrix_set, 3)
                  CALL cp_fm_release(matrix_set(ikp_local, ic, ispin))
               END DO
            END DO
         END DO
         DEALLOCATE (matrix_set)
      END IF

   END SUBROUTINE tddfpt_kpoint_release_matrix_set

! **************************************************************************************************
!> \brief Allocate a k-point TDDFPT response vector in transition space.
!> \param response ...
!> \param ntrans ...
! **************************************************************************************************
   SUBROUTINE tddfpt_kpoint_response_create(response, ntrans)
      TYPE(tddfpt_kpoint_response_type), INTENT(out)     :: response
      INTEGER, INTENT(in)                                :: ntrans

      CPASSERT(ntrans >= 0)

      response%ntrans = ntrans
      ALLOCATE (response%coeff_re(ntrans), response%coeff_im(ntrans), &
                response%sigma_re(ntrans), response%sigma_im(ntrans))
      CALL tddfpt_kpoint_response_zero(response)

   END SUBROUTINE tddfpt_kpoint_response_create

! **************************************************************************************************
!> \brief Release a k-point TDDFPT response vector.
!> \param response ...
! **************************************************************************************************
   SUBROUTINE tddfpt_kpoint_response_release(response)
      TYPE(tddfpt_kpoint_response_type), INTENT(inout)   :: response

      IF (ALLOCATED(response%coeff_re)) DEALLOCATE (response%coeff_re)
      IF (ALLOCATED(response%coeff_im)) DEALLOCATE (response%coeff_im)
      IF (ALLOCATED(response%sigma_re)) DEALLOCATE (response%sigma_re)
      IF (ALLOCATED(response%sigma_im)) DEALLOCATE (response%sigma_im)
      response%ntrans = 0

   END SUBROUTINE tddfpt_kpoint_response_release

! **************************************************************************************************
!> \brief Reset coefficients and operator result of a k-point TDDFPT response vector.
!> \param response ...
! **************************************************************************************************
   SUBROUTINE tddfpt_kpoint_response_zero(response)
      TYPE(tddfpt_kpoint_response_type), INTENT(inout)   :: response

      IF (ALLOCATED(response%coeff_re)) response%coeff_re = 0.0_dp
      IF (ALLOCATED(response%coeff_im)) response%coeff_im = 0.0_dp
      IF (ALLOCATED(response%sigma_re)) response%sigma_re = 0.0_dp
      IF (ALLOCATED(response%sigma_im)) response%sigma_im = 0.0_dp

   END SUBROUTINE tddfpt_kpoint_response_zero

! **************************************************************************************************
!> \brief Set a single real or imaginary unit transition in a k-point TDDFPT response vector.
!> \param response ...
!> \param itrans ...
!> \param imaginary ...
! **************************************************************************************************
   SUBROUTINE tddfpt_kpoint_response_set_unit(response, itrans, imaginary)
      TYPE(tddfpt_kpoint_response_type), INTENT(inout)   :: response
      INTEGER, INTENT(in)                                :: itrans
      LOGICAL, INTENT(in), OPTIONAL                      :: imaginary

      LOGICAL                                            :: use_imaginary

      CPASSERT(itrans >= 1)
      CPASSERT(itrans <= response%ntrans)

      use_imaginary = .FALSE.
      IF (PRESENT(imaginary)) use_imaginary = imaginary

      CALL tddfpt_kpoint_response_zero(response)
      IF (use_imaginary) THEN
         response%coeff_im(itrans) = 1.0_dp
      ELSE
         response%coeff_re(itrans) = 1.0_dp
      END IF

   END SUBROUTINE tddfpt_kpoint_response_set_unit

! **************************************************************************************************
!> \brief Apply the independent-particle energy-difference term in transition space.
!> \param response ...
!> \param transitions ...
! **************************************************************************************************
   SUBROUTINE tddfpt_kpoint_apply_energy_diff(response, transitions)
      TYPE(tddfpt_kpoint_response_type), INTENT(inout)   :: response
      TYPE(tddfpt_kpoint_transition_type), &
         DIMENSION(:), INTENT(in)                        :: transitions

      CPASSERT(response%ntrans == SIZE(transitions))

      response%sigma_re(:) = response%sigma_re(:) + transitions(:)%energy*response%coeff_re(:)
      response%sigma_im(:) = response%sigma_im(:) + transitions(:)%energy*response%coeff_im(:)

   END SUBROUTINE tddfpt_kpoint_apply_energy_diff

! **************************************************************************************************
!> \brief Weighted complex dot product between response coefficients.
!> \param left ...
!> \param right ...
!> \param transitions ...
!> \return ...
! **************************************************************************************************
   COMPLEX(kind=dp) FUNCTION tddfpt_kpoint_response_dot(left, right, transitions) RESULT(dot_product)
      TYPE(tddfpt_kpoint_response_type), INTENT(in)      :: left, right
      TYPE(tddfpt_kpoint_transition_type), &
         DIMENSION(:), INTENT(in)                        :: transitions

      REAL(kind=dp)                                      :: dot_im, dot_re

      CPASSERT(left%ntrans == right%ntrans)
      CPASSERT(left%ntrans == SIZE(transitions))

      dot_re = SUM(transitions(:)%weight*(left%coeff_re(:)*right%coeff_re(:) + &
                                          left%coeff_im(:)*right%coeff_im(:)))
      dot_im = SUM(transitions(:)%weight*(left%coeff_re(:)*right%coeff_im(:) - &
                                          left%coeff_im(:)*right%coeff_re(:)))
      dot_product = CMPLX(dot_re, dot_im, KIND=dp)

   END FUNCTION tddfpt_kpoint_response_dot

! **************************************************************************************************
!> \brief Weighted norm of a k-point TDDFPT response vector.
!> \param response ...
!> \param transitions ...
!> \return ...
! **************************************************************************************************
   REAL(kind=dp) FUNCTION tddfpt_kpoint_response_norm(response, transitions) RESULT(response_norm)
      TYPE(tddfpt_kpoint_response_type), INTENT(in)      :: response
      TYPE(tddfpt_kpoint_transition_type), &
         DIMENSION(:), INTENT(in)                        :: transitions

      response_norm = SQRT(MAX(0.0_dp, REAL(tddfpt_kpoint_response_dot(response, response, transitions), dp)))

   END FUNCTION tddfpt_kpoint_response_norm

! **************************************************************************************************
!> \brief Normalize a k-point TDDFPT response vector with the weighted S(k) metric.
!> \param response ...
!> \param transitions ...
! **************************************************************************************************
   SUBROUTINE tddfpt_kpoint_response_normalize(response, transitions)
      TYPE(tddfpt_kpoint_response_type), INTENT(inout)   :: response
      TYPE(tddfpt_kpoint_transition_type), &
         DIMENSION(:), INTENT(in)                        :: transitions

      REAL(kind=dp)                                      :: response_norm

      response_norm = tddfpt_kpoint_response_norm(response, transitions)
      IF (response_norm <= 0.0_dp) CPABORT("Cannot normalize empty k-point TDDFPT response")
      response%coeff_re(:) = response%coeff_re(:)/response_norm
      response%coeff_im(:) = response%coeff_im(:)/response_norm

   END SUBROUTINE tddfpt_kpoint_response_normalize

! **************************************************************************************************
!> \brief Project one normalized k-point TDDFPT response vector out of another.
!> \param response ...
!> \param basis ...
!> \param transitions ...
! **************************************************************************************************
   SUBROUTINE tddfpt_kpoint_response_project_out(response, basis, transitions)
      TYPE(tddfpt_kpoint_response_type), INTENT(inout)   :: response
      TYPE(tddfpt_kpoint_response_type), INTENT(in)      :: basis
      TYPE(tddfpt_kpoint_transition_type), &
         DIMENSION(:), INTENT(in)                        :: transitions

      COMPLEX(kind=dp)                                   :: coeff

      coeff = tddfpt_kpoint_response_dot(basis, response, transitions)
      response%coeff_re(:) = response%coeff_re(:) - &
                             REAL(coeff, dp)*basis%coeff_re(:) + AIMAG(coeff)*basis%coeff_im(:)
      response%coeff_im(:) = response%coeff_im(:) - &
                             REAL(coeff, dp)*basis%coeff_im(:) - AIMAG(coeff)*basis%coeff_re(:)

   END SUBROUTINE tddfpt_kpoint_response_project_out

! **************************************************************************************************
!> \brief Copy response coefficients with a k-weight metric scaling.
!> \param dest ...
!> \param source ...
!> \param transitions ...
!> \param weight_power ...
! **************************************************************************************************
   SUBROUTINE tddfpt_kpoint_response_copy_scaled_coeff(dest, source, transitions, weight_power)
      TYPE(tddfpt_kpoint_response_type), INTENT(inout)   :: dest
      TYPE(tddfpt_kpoint_response_type), INTENT(in)      :: source
      TYPE(tddfpt_kpoint_transition_type), &
         DIMENSION(:), INTENT(in)                        :: transitions
      REAL(kind=dp), INTENT(in)                          :: weight_power

      INTEGER                                            :: itrans
      REAL(kind=dp)                                      :: scale

      CPASSERT(dest%ntrans == source%ntrans)
      CPASSERT(dest%ntrans == SIZE(transitions))

      CALL tddfpt_kpoint_response_zero(dest)
      DO itrans = 1, dest%ntrans
         CPASSERT(transitions(itrans)%weight > 0.0_dp)
         scale = transitions(itrans)%weight**weight_power
         dest%coeff_re(itrans) = scale*source%coeff_re(itrans)
         dest%coeff_im(itrans) = scale*source%coeff_im(itrans)
      END DO

   END SUBROUTINE tddfpt_kpoint_response_copy_scaled_coeff

! **************************************************************************************************
!> \brief Accumulate an operator result with a k-weight metric scaling.
!> \param dest ...
!> \param source ...
!> \param transitions ...
!> \param weight_power ...
! **************************************************************************************************
   SUBROUTINE tddfpt_kpoint_response_add_scaled_sigma(dest, source, transitions, weight_power)
      TYPE(tddfpt_kpoint_response_type), INTENT(inout)   :: dest
      TYPE(tddfpt_kpoint_response_type), INTENT(in)      :: source
      TYPE(tddfpt_kpoint_transition_type), &
         DIMENSION(:), INTENT(in)                        :: transitions
      REAL(kind=dp), INTENT(in)                          :: weight_power

      INTEGER                                            :: itrans
      REAL(kind=dp)                                      :: scale

      CPASSERT(dest%ntrans == source%ntrans)
      CPASSERT(dest%ntrans == SIZE(transitions))

      DO itrans = 1, dest%ntrans
         CPASSERT(transitions(itrans)%weight > 0.0_dp)
         scale = transitions(itrans)%weight**weight_power
         dest%sigma_re(itrans) = dest%sigma_re(itrans) + scale*source%sigma_re(itrans)
         dest%sigma_im(itrans) = dest%sigma_im(itrans) + scale*source%sigma_im(itrans)
      END DO

   END SUBROUTINE tddfpt_kpoint_response_add_scaled_sigma

! **************************************************************************************************
!> \brief Build AO transition-density matrices from a k-point TDDFPT response vector.
!> \param qs_env ...
!> \param response ...
!> \param transitions ...
!> \param pmat_ext ...
! **************************************************************************************************
   SUBROUTINE tddfpt_kpoint_response_to_pmat(qs_env, response, transitions, pmat_ext)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddfpt_kpoint_response_type), INTENT(in)      :: response
      TYPE(tddfpt_kpoint_transition_type), &
         DIMENSION(:), INTENT(in)                        :: transitions
      TYPE(cp_fm_type), DIMENSION(:, :, :), &
         INTENT(inout)                                   :: pmat_ext

      CHARACTER(len=*), PARAMETER :: routineN = 'tddfpt_kpoint_response_to_pmat'

      INTEGER                                            :: handle, ikp, ikp_local, ispin, itrans, &
                                                            nao, nspins
      INTEGER, DIMENSION(2)                              :: kp_range
      LOGICAL                                            :: my_kpgrp
      REAL(kind=dp)                                      :: xi, xr
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: delta_im, delta_re, density_im, &
                                                            density_re, occ_im, occ_re, virt_im, &
                                                            virt_re
      TYPE(cp_fm_type), POINTER                          :: mo_coeff_im, mo_coeff_re
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_env_p_type), DIMENSION(:), POINTER     :: kp_env
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mo_set_type), DIMENSION(:, :), POINTER        :: mos_kp

      CALL timeset(routineN, handle)

      NULLIFY (dft_control, kp, kp_env, kpoints, mo_coeff_im, mo_coeff_re, mos_kp)
      CALL get_qs_env(qs_env, dft_control=dft_control, kpoints=kpoints)
      CPASSERT(ASSOCIATED(dft_control))
      CPASSERT(ASSOCIATED(kpoints))
      CPASSERT(response%ntrans == SIZE(transitions))
      CPASSERT(SIZE(pmat_ext, 2) >= 2)

      CALL get_kpoint_info(kpoints, kp_range=kp_range, kp_env=kp_env)
      CPASSERT(ASSOCIATED(kp_env))

      nspins = dft_control%nspins
      CPASSERT(SIZE(pmat_ext, 3) >= nspins)

      DO ikp_local = 1, SIZE(pmat_ext, 1)
         DO ispin = 1, nspins
            CALL cp_fm_set_all(pmat_ext(ikp_local, 1, ispin), 0.0_dp)
            CALL cp_fm_set_all(pmat_ext(ikp_local, 2, ispin), 0.0_dp)
         END DO
      END DO

      DO ikp = 1, kpoints%nkp
         my_kpgrp = (ikp >= kp_range(1) .AND. ikp <= kp_range(2))
         IF (.NOT. my_kpgrp) CYCLE
         ikp_local = ikp - kp_range(1) + 1
         kp => kp_env(ikp_local)%kpoint_env
         mos_kp => kp%mos
         CPASSERT(ASSOCIATED(mos_kp))

         DO ispin = 1, nspins
            CALL get_mo_set(mos_kp(1, ispin), mo_coeff=mo_coeff_re, nao=nao)
            CALL get_mo_set(mos_kp(2, ispin), mo_coeff=mo_coeff_im)
            CPASSERT(ASSOCIATED(mo_coeff_re))
            CPASSERT(ASSOCIATED(mo_coeff_im))

            ALLOCATE (density_re(nao, nao), density_im(nao, nao), delta_re(nao, 1), &
                      delta_im(nao, 1), occ_re(nao, 1), occ_im(nao, 1), &
                      virt_re(nao, 1), virt_im(nao, 1))
            density_re = 0.0_dp
            density_im = 0.0_dp

            DO itrans = 1, response%ntrans
               IF (transitions(itrans)%ikp /= ikp) CYCLE
               IF (transitions(itrans)%ispin /= ispin) CYCLE
               xr = response%coeff_re(itrans)
               xi = response%coeff_im(itrans)
               IF (xr == 0.0_dp .AND. xi == 0.0_dp) CYCLE

               CALL cp_fm_get_submatrix(mo_coeff_re, occ_re, start_row=1, &
                                        start_col=transitions(itrans)%iocc, n_rows=nao, n_cols=1)
               CALL cp_fm_get_submatrix(mo_coeff_im, occ_im, start_row=1, &
                                        start_col=transitions(itrans)%iocc, n_rows=nao, n_cols=1)
               CALL cp_fm_get_submatrix(mo_coeff_re, virt_re, start_row=1, &
                                        start_col=transitions(itrans)%ivirt, n_rows=nao, n_cols=1)
               CALL cp_fm_get_submatrix(mo_coeff_im, virt_im, start_row=1, &
                                        start_col=transitions(itrans)%ivirt, n_rows=nao, n_cols=1)

               delta_re(:, 1) = xr*virt_re(:, 1) - xi*virt_im(:, 1)
               delta_im(:, 1) = xr*virt_im(:, 1) + xi*virt_re(:, 1)

               density_re = density_re + 0.5_dp*(MATMUL(occ_re, TRANSPOSE(delta_re)) + &
                                                 MATMUL(occ_im, TRANSPOSE(delta_im)) + &
                                                 MATMUL(delta_re, TRANSPOSE(occ_re)) + &
                                                 MATMUL(delta_im, TRANSPOSE(occ_im)))
               density_im = density_im + 0.5_dp*(MATMUL(occ_im, TRANSPOSE(delta_re)) - &
                                                 MATMUL(occ_re, TRANSPOSE(delta_im)) + &
                                                 MATMUL(delta_im, TRANSPOSE(occ_re)) - &
                                                 MATMUL(delta_re, TRANSPOSE(occ_im)))
            END DO

            CALL cp_fm_set_submatrix(pmat_ext(ikp_local, 1, ispin), density_re, &
                                     start_row=1, start_col=1, n_rows=nao, n_cols=nao)
            CALL cp_fm_set_submatrix(pmat_ext(ikp_local, 2, ispin), density_im, &
                                     start_row=1, start_col=1, n_rows=nao, n_cols=nao)
            DEALLOCATE (density_re, density_im, delta_re, delta_im, occ_re, occ_im, virt_re, virt_im)
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE tddfpt_kpoint_response_to_pmat

! **************************************************************************************************
!> \brief Build and collocate a real-space response density from a k-point response vector.
!> \param qs_env ...
!> \param response ...
!> \param transitions ...
!> \param rho1 ...
!> \param pmat_ext ...
! **************************************************************************************************
   SUBROUTINE tddfpt_kpoint_response_update_rho(qs_env, response, transitions, rho1, pmat_ext)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddfpt_kpoint_response_type), INTENT(in)      :: response
      TYPE(tddfpt_kpoint_transition_type), &
         DIMENSION(:), INTENT(in)                        :: transitions
      TYPE(qs_rho_type), POINTER                         :: rho1
      TYPE(cp_fm_type), DIMENSION(:, :, :), &
         INTENT(inout)                                   :: pmat_ext

      CHARACTER(len=*), PARAMETER :: routineN = 'tddfpt_kpoint_response_update_rho'

      INTEGER                                            :: handle, ic, nao
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_all
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fmwork
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s_kp, rho1_ao
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_kp

      CALL timeset(routineN, handle)

      NULLIFY (blacs_env_all, fm_struct, kpoints, matrix_s_kp, para_env, rho1_ao, sab_kp)
      CPASSERT(ASSOCIATED(rho1))
      CALL get_qs_env(qs_env, kpoints=kpoints, matrix_s_kp=matrix_s_kp)
      CPASSERT(ASSOCIATED(kpoints))
      CPASSERT(ASSOCIATED(matrix_s_kp))
      CALL get_kpoint_info(kpoints, blacs_env_all=blacs_env_all, para_env=para_env, sab_nl=sab_kp)
      CPASSERT(ASSOCIATED(blacs_env_all))
      CPASSERT(ASSOCIATED(para_env))
      CPASSERT(ASSOCIATED(sab_kp))

      CALL tddfpt_kpoint_response_to_pmat(qs_env, response, transitions, pmat_ext)

      CALL qs_rho_rebuild(rho1, qs_env, rebuild_ao=.FALSE., rebuild_grids=.TRUE.)
      CALL qs_rho_get(rho1, rho_ao_kp=rho1_ao)
      CPASSERT(ASSOCIATED(rho1_ao))

      CALL dbcsr_get_info(matrix_s_kp(1, 1)%matrix, nfullrows_total=nao)
      ALLOCATE (fmwork(2))
      CALL cp_fm_struct_create(fm_struct, nrow_global=nao, ncol_global=nao, &
                               para_env=para_env, context=blacs_env_all)
      DO ic = 1, 2
         CALL cp_fm_create(fmwork(ic), fm_struct)
      END DO
      CALL cp_fm_struct_release(fm_struct)

      CALL kpoint_density_transform(kpoints, rho1_ao, .FALSE., matrix_s_kp(1, 1)%matrix, &
                                    sab_kp, fmwork, pmat_ext=pmat_ext)
      CALL qs_rho_update_rho(rho1, qs_env)

      DO ic = 1, 2
         CALL cp_fm_release(fmwork(ic))
      END DO
      DEALLOCATE (fmwork)

      CALL timestop(handle)

   END SUBROUTINE tddfpt_kpoint_response_update_rho

! **************************************************************************************************
!> \brief Fourier-transform real-space AO operator matrices to local k-point dense matrices.
!> \param qs_env ...
!> \param hmat_rs ...
!> \param hmat_ext ...
! **************************************************************************************************
   SUBROUTINE tddfpt_kpoint_hmat_from_rs(qs_env, hmat_rs, hmat_ext)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: hmat_rs
      TYPE(cp_fm_type), DIMENSION(:, :, :), &
         INTENT(inout)                                   :: hmat_ext

      CHARACTER(len=*), PARAMETER :: routineN = 'tddfpt_kpoint_hmat_from_rs'

      INTEGER                                            :: handle, ikp, ikp_local, ispin, nspins
      INTEGER, DIMENSION(2)                              :: kp_range
      TYPE(dbcsr_type), POINTER                          :: cmatrix, rmatrix
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_kp

      CALL timeset(routineN, handle)

      NULLIFY (cmatrix, dft_control, kpoints, rmatrix, sab_kp)
      CALL get_qs_env(qs_env, dft_control=dft_control, kpoints=kpoints)
      CPASSERT(ASSOCIATED(dft_control))
      CPASSERT(ASSOCIATED(kpoints))
      CPASSERT(ASSOCIATED(hmat_rs))
      CALL get_kpoint_info(kpoints, kp_range=kp_range, sab_nl=sab_kp)
      CPASSERT(ASSOCIATED(sab_kp))

      nspins = dft_control%nspins
      CPASSERT(SIZE(hmat_ext, 2) >= 2)
      CPASSERT(SIZE(hmat_ext, 3) >= nspins)

      ALLOCATE (rmatrix, cmatrix)
      CALL dbcsr_create(rmatrix, template=hmat_rs(1, 1)%matrix, matrix_type=dbcsr_type_symmetric)
      CALL dbcsr_create(cmatrix, template=hmat_rs(1, 1)%matrix, matrix_type=dbcsr_type_antisymmetric)
      CALL cp_dbcsr_alloc_block_from_nbl(rmatrix, sab_kp)
      CALL cp_dbcsr_alloc_block_from_nbl(cmatrix, sab_kp)

      DO ikp = 1, kpoints%nkp
         IF (ikp < kp_range(1) .OR. ikp > kp_range(2)) CYCLE
         ikp_local = ikp - kp_range(1) + 1
         DO ispin = 1, nspins
            CALL dbcsr_set(rmatrix, 0.0_dp)
            CALL dbcsr_set(cmatrix, 0.0_dp)
            CALL rskp_transform(rmatrix=rmatrix, cmatrix=cmatrix, rsmat=hmat_rs, &
                                ispin=ispin, xkp=kpoints%xkp(:, ikp), &
                                cell_to_index=kpoints%cell_to_index, sab_nl=sab_kp)
            CALL copy_dbcsr_to_fm(rmatrix, hmat_ext(ikp_local, 1, ispin))
            CALL copy_dbcsr_to_fm(cmatrix, hmat_ext(ikp_local, 2, ispin))
         END DO
      END DO

      CALL dbcsr_deallocate_matrix(rmatrix)
      CALL dbcsr_deallocate_matrix(cmatrix)

      CALL timestop(handle)

   END SUBROUTINE tddfpt_kpoint_hmat_from_rs

! **************************************************************************************************
!> \brief Project k-point AO operator matrices onto the TDDFPT transition space.
!> \param qs_env ...
!> \param transitions ...
!> \param hmat_ext ...
!> \param response ...
! **************************************************************************************************
   SUBROUTINE tddfpt_kpoint_hmat_to_response(qs_env, transitions, hmat_ext, response)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddfpt_kpoint_transition_type), &
         DIMENSION(:), INTENT(in)                        :: transitions
      TYPE(cp_fm_type), DIMENSION(:, :, :), INTENT(in)   :: hmat_ext
      TYPE(tddfpt_kpoint_response_type), INTENT(inout)   :: response

      CHARACTER(len=*), PARAMETER :: routineN = 'tddfpt_kpoint_hmat_to_response'

      INTEGER                                            :: handle, ikp, ikp_local, ispin, itrans, &
                                                            nao, nspins
      INTEGER, DIMENSION(2)                              :: kp_range
      LOGICAL                                            :: my_kpgrp
      REAL(kind=dp)                                      :: sigma_im, sigma_re
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: work_im, work_re
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: hmat_im, hmat_re, occ_im, occ_re, &
                                                            virt_im, virt_re
      TYPE(cp_fm_type), POINTER                          :: mo_coeff_im, mo_coeff_re
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_env_p_type), DIMENSION(:), POINTER     :: kp_env
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mo_set_type), DIMENSION(:, :), POINTER        :: mos_kp
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      NULLIFY (dft_control, kp, kp_env, kpoints, mo_coeff_im, mo_coeff_re, mos_kp, para_env)
      CALL get_qs_env(qs_env, dft_control=dft_control, kpoints=kpoints)
      CPASSERT(ASSOCIATED(dft_control))
      CPASSERT(ASSOCIATED(kpoints))
      CPASSERT(response%ntrans == SIZE(transitions))
      CPASSERT(SIZE(hmat_ext, 2) >= 2)

      CALL get_kpoint_info(kpoints, kp_range=kp_range, kp_env=kp_env, para_env=para_env)
      CPASSERT(ASSOCIATED(kp_env))
      CPASSERT(ASSOCIATED(para_env))

      nspins = dft_control%nspins
      CPASSERT(SIZE(hmat_ext, 3) >= nspins)

      DO ikp = 1, kpoints%nkp
         my_kpgrp = (ikp >= kp_range(1) .AND. ikp <= kp_range(2))
         IF (.NOT. my_kpgrp) CYCLE
         ikp_local = ikp - kp_range(1) + 1
         kp => kp_env(ikp_local)%kpoint_env
         mos_kp => kp%mos
         CPASSERT(ASSOCIATED(mos_kp))

         DO ispin = 1, nspins
            CALL get_mo_set(mos_kp(1, ispin), mo_coeff=mo_coeff_re, nao=nao)
            CALL get_mo_set(mos_kp(2, ispin), mo_coeff=mo_coeff_im)
            CPASSERT(ASSOCIATED(mo_coeff_re))
            CPASSERT(ASSOCIATED(mo_coeff_im))

            ALLOCATE (hmat_re(nao, nao), hmat_im(nao, nao), occ_re(nao, 1), &
                      occ_im(nao, 1), virt_re(nao, 1), virt_im(nao, 1), &
                      work_re(nao), work_im(nao))
            CALL cp_fm_get_submatrix(hmat_ext(ikp_local, 1, ispin), hmat_re, &
                                     start_row=1, start_col=1, n_rows=nao, n_cols=nao)
            CALL cp_fm_get_submatrix(hmat_ext(ikp_local, 2, ispin), hmat_im, &
                                     start_row=1, start_col=1, n_rows=nao, n_cols=nao)

            DO itrans = 1, response%ntrans
               IF (transitions(itrans)%ikp /= ikp) CYCLE
               IF (transitions(itrans)%ispin /= ispin) CYCLE

               CALL cp_fm_get_submatrix(mo_coeff_re, occ_re, start_row=1, &
                                        start_col=transitions(itrans)%iocc, n_rows=nao, n_cols=1)
               CALL cp_fm_get_submatrix(mo_coeff_im, occ_im, start_row=1, &
                                        start_col=transitions(itrans)%iocc, n_rows=nao, n_cols=1)
               CALL cp_fm_get_submatrix(mo_coeff_re, virt_re, start_row=1, &
                                        start_col=transitions(itrans)%ivirt, n_rows=nao, n_cols=1)
               CALL cp_fm_get_submatrix(mo_coeff_im, virt_im, start_row=1, &
                                        start_col=transitions(itrans)%ivirt, n_rows=nao, n_cols=1)

               work_re(:) = MATMUL(hmat_re, occ_re(:, 1)) - MATMUL(hmat_im, occ_im(:, 1))
               work_im(:) = MATMUL(hmat_re, occ_im(:, 1)) + MATMUL(hmat_im, occ_re(:, 1))
               sigma_re = SUM(virt_re(:, 1)*work_re(:) + virt_im(:, 1)*work_im(:))
               sigma_im = SUM(virt_re(:, 1)*work_im(:) - virt_im(:, 1)*work_re(:))
               response%sigma_re(itrans) = response%sigma_re(itrans) + sigma_re
               response%sigma_im(itrans) = response%sigma_im(itrans) + sigma_im
            END DO
            DEALLOCATE (hmat_re, hmat_im, occ_re, occ_im, virt_re, virt_im, work_re, work_im)
         END DO
      END DO

      CALL para_env%sum(response%sigma_re)
      CALL para_env%sum(response%sigma_im)

      CALL timestop(handle)

   END SUBROUTINE tddfpt_kpoint_hmat_to_response

! **************************************************************************************************
!> \brief Return the zero-based transition offset of a spin block within a single k point.
!> \param ispin ...
!> \param homo_spin ...
!> \param nvirt_spin ...
!> \return ...
! **************************************************************************************************
   INTEGER FUNCTION tddfpt_kpoint_transition_offset(ispin, homo_spin, nvirt_spin) RESULT(offset)
      INTEGER, INTENT(in)                                :: ispin
      INTEGER, DIMENSION(:), INTENT(in)                  :: homo_spin, nvirt_spin

      INTEGER                                            :: jspin

      offset = 0
      DO jspin = 1, ispin - 1
         offset = offset + homo_spin(jspin)*nvirt_spin(jspin)
      END DO

   END FUNCTION tddfpt_kpoint_transition_offset

! **************************************************************************************************
!> \brief Build the k/spin/occupied/virtual transition map for k-point TDDFPT.
!> \param qs_env ...
!> \param transitions ...
!> \param ntrans_kpoint ...
!> \param nmo_spin ...
!> \param homo_spin ...
!> \param nvirt_spin ...
!> \param nao_spin ...
! **************************************************************************************************
   SUBROUTINE tddfpt_kpoint_init_transitions(qs_env, transitions, ntrans_kpoint, &
                                             nmo_spin, homo_spin, nvirt_spin, nao_spin)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddfpt_kpoint_transition_type), ALLOCATABLE, &
         DIMENSION(:), INTENT(out)                       :: transitions
      INTEGER, INTENT(out)                               :: ntrans_kpoint
      INTEGER, DIMENSION(:), INTENT(out)                 :: nmo_spin, homo_spin, nvirt_spin, nao_spin

      CHARACTER(len=*), PARAMETER :: routineN = 'tddfpt_kpoint_init_transitions'

      INTEGER                                            :: handle, ikp, ikp_local, iocc, ispin, &
                                                            ivirt, nkp, nkp_local, nspins, offset, &
                                                            trans_index
      INTEGER, DIMENSION(2)                              :: kp_range
      LOGICAL                                            :: my_kpgrp
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues_kp
      REAL(kind=dp), DIMENSION(:), POINTER               :: eigenvalues, wkp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_env_p_type), DIMENSION(:), POINTER     :: kp_env
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mo_set_type), DIMENSION(:, :), POINTER        :: mos_kp
      TYPE(mp_para_env_type), POINTER                    :: para_env, para_env_kp

      CALL timeset(routineN, handle)

      NULLIFY (dft_control, eigenvalues, kp, kp_env, kpoints, mos_kp, para_env, para_env_kp, wkp)
      CALL get_qs_env(qs_env, dft_control=dft_control, kpoints=kpoints)
      CPASSERT(ASSOCIATED(dft_control))
      CPASSERT(ASSOCIATED(kpoints))

      CALL get_kpoint_info(kpoints, nkp=nkp, kp_range=kp_range, kp_env=kp_env, &
                           para_env=para_env, para_env_kp=para_env_kp, wkp=wkp)
      CPASSERT(ASSOCIATED(kp_env))
      CPASSERT(ASSOCIATED(para_env))
      CPASSERT(ASSOCIATED(para_env_kp))
      CPASSERT(ASSOCIATED(wkp))

      nspins = dft_control%nspins
      CPASSERT(SIZE(nmo_spin) >= nspins)
      CPASSERT(SIZE(homo_spin) >= nspins)
      CPASSERT(SIZE(nvirt_spin) >= nspins)
      CPASSERT(SIZE(nao_spin) >= nspins)

      nmo_spin = 0
      homo_spin = 0
      nvirt_spin = 0
      nao_spin = 0

      nkp_local = MAX(0, kp_range(2) - kp_range(1) + 1)
      IF (nkp_local > 0) THEN
         kp => kp_env(1)%kpoint_env
         mos_kp => kp%mos
         CPASSERT(ASSOCIATED(mos_kp))
         CPASSERT(SIZE(mos_kp, 2) == nspins)
         DO ispin = 1, nspins
            CALL get_mo_set(mos_kp(1, ispin), nmo=nmo_spin(ispin), &
                            homo=homo_spin(ispin), nao=nao_spin(ispin))
         END DO
      END IF
      CALL para_env%max(nmo_spin)
      CALL para_env%max(homo_spin)
      CALL para_env%max(nao_spin)

      ntrans_kpoint = 0
      DO ispin = 1, nspins
         nvirt_spin(ispin) = nmo_spin(ispin) - homo_spin(ispin)
         IF (homo_spin(ispin) <= 0 .OR. nvirt_spin(ispin) <= 0) &
            CPABORT("Missing occupied or unoccupied MO for k-point TDDFPT")
         ntrans_kpoint = ntrans_kpoint + homo_spin(ispin)*nvirt_spin(ispin)
      END DO
      IF (nkp*ntrans_kpoint <= 0) CPABORT("No k-point TDDFPT transitions available")

      ALLOCATE (transitions(nkp*ntrans_kpoint))

      DO ikp = 1, nkp
         my_kpgrp = (ikp >= kp_range(1) .AND. ikp <= kp_range(2))
         IF (my_kpgrp) THEN
            ikp_local = ikp - kp_range(1) + 1
            kp => kp_env(ikp_local)%kpoint_env
            mos_kp => kp%mos
         ELSE
            NULLIFY (kp, mos_kp)
         END IF

         DO ispin = 1, nspins
            ALLOCATE (eigenvalues_kp(nmo_spin(ispin)))
            eigenvalues_kp = 0.0_dp

            IF (my_kpgrp) THEN
               CALL get_mo_set(mos_kp(1, ispin), eigenvalues=eigenvalues)
               CPASSERT(ASSOCIATED(eigenvalues))
               IF (para_env_kp%is_source()) &
                  eigenvalues_kp(1:nmo_spin(ispin)) = eigenvalues(1:nmo_spin(ispin))
            END IF
            CALL para_env%sum(eigenvalues_kp)

            offset = tddfpt_kpoint_transition_offset(ispin, homo_spin, nvirt_spin)
            DO iocc = 1, homo_spin(ispin)
               DO ivirt = homo_spin(ispin) + 1, nmo_spin(ispin)
                  trans_index = (ikp - 1)*ntrans_kpoint + offset + &
                                (iocc - 1)*nvirt_spin(ispin) + ivirt - homo_spin(ispin)
                  transitions(trans_index)%ikp = ikp
                  transitions(trans_index)%ispin = ispin
                  transitions(trans_index)%iocc = iocc
                  transitions(trans_index)%ivirt = ivirt
                  transitions(trans_index)%weight = wkp(ikp)
                  transitions(trans_index)%energy = eigenvalues_kp(ivirt) - eigenvalues_kp(iocc)
               END DO
            END DO

            DEALLOCATE (eigenvalues_kp)
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE tddfpt_kpoint_init_transitions

END MODULE qs_tddfpt2_kpoint_utils
