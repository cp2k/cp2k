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

   PUBLIC :: tddfpt_kpoint_init_transitions, tddfpt_kpoint_transition_offset

CONTAINS

! **************************************************************************************************
!> \brief Return the zero-based transition offset of a spin block within a single k point.
!> \param ispin ...
!> \param homo_spin ...
!> \param nvirt_spin ...
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
