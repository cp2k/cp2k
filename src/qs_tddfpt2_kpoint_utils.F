!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE qs_tddfpt2_kpoint_utils
   USE cp_control_types,                ONLY: dft_control_type
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_env_p_type,&
                                              kpoint_env_type,&
                                              kpoint_type
   USE message_passing,                 ONLY: mp_para_env_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
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

   PUBLIC :: tddfpt_kpoint_apply_energy_diff, tddfpt_kpoint_init_transitions, &
             tddfpt_kpoint_response_create, tddfpt_kpoint_response_dot, tddfpt_kpoint_response_norm, &
             tddfpt_kpoint_response_normalize, tddfpt_kpoint_response_project_out, &
             tddfpt_kpoint_response_release, tddfpt_kpoint_response_set_unit, &
             tddfpt_kpoint_response_zero, tddfpt_kpoint_transition_offset

CONTAINS

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
!> \brief Set a single real unit transition in a k-point TDDFPT response vector.
!> \param response ...
!> \param itrans ...
! **************************************************************************************************
   SUBROUTINE tddfpt_kpoint_response_set_unit(response, itrans)
      TYPE(tddfpt_kpoint_response_type), INTENT(inout)   :: response
      INTEGER, INTENT(in)                                :: itrans

      CPASSERT(itrans >= 1)
      CPASSERT(itrans <= response%ntrans)

      CALL tddfpt_kpoint_response_zero(response)
      response%coeff_re(itrans) = 1.0_dp

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
