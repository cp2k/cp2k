!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculate Energy Decomposition analysis
!> \par History
!>      07.2023 created [JGH]
!> \author JGH
! **************************************************************************************************
MODULE ed_analysis
   USE bibliography,                    ONLY: Eriksen2020,&
                                              cite_reference
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE iao_analysis,                    ONLY: iao_wfn_analysis
   USE iao_types,                       ONLY: iao_env_type,&
                                              iao_set_default
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   TYPE block_type
      REAL(KIND=dp), DIMENSION(:, :), ALLOCATABLE  :: mat
   END TYPE block_type

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ed_analysis'

   PUBLIC ::  edmf_analysis

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param input_section ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE edmf_analysis(qs_env, input_section, unit_nr)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: input_section
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(len=*), PARAMETER                        :: routineN = 'edmf_analysis'

      INTEGER                                            :: handle, ispin, nimages, nspin
      LOGICAL                                            :: explicit
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: c_iao_coef, c_loc_coef
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(iao_env_type)                                 :: iao_env
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos

      CALL section_vals_get(input_section, explicit=explicit)
      IF (.NOT. explicit) RETURN

      CALL timeset(routineN, handle)

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(/,T2,A)') '!-----------------------------------------------------------------------------!'
         WRITE (UNIT=unit_nr, FMT="(T26,A)") "ENERGY DECOMPOSITION ANALYSIS"
         WRITE (UNIT=unit_nr, FMT="(T22,A)") "Janus J Eriksen, JCP 153 214109 (2020)"
         WRITE (unit_nr, '(T2,A)') '!-----------------------------------------------------------------------------!'
      END IF
      CALL cite_reference(Eriksen2020)

      ! k-points?
      CALL get_qs_env(qs_env, dft_control=dft_control)
      nimages = dft_control%nimages
      IF (nimages > 1) THEN
         IF (unit_nr > 0) THEN
            WRITE (UNIT=unit_nr, FMT="(T2,A)") &
               "K-Points: MAO's determined and analyzed using Gamma-Point only."
         END IF
      END IF

      CALL get_qs_env(qs_env, mos=mos)
      nspin = dft_control%nspins
      ALLOCATE (c_loc_coef(nspin))

      CALL iao_set_default(iao_env)
      iao_env%do_iao = .TRUE.
      iao_env%do_bondorbitals = .TRUE.
      DO ispin = 1, nspin
         CALL get_mo_set(mos(ispin), mo_coeff=mo_coeff)
         CALL cp_fm_create(c_loc_coef(ispin), mo_coeff%matrix_struct)
      END DO
      CALL iao_wfn_analysis(qs_env, iao_env, unit_nr, c_iao_coef, c_loc_coef)
      CALL cp_fm_release(c_iao_coef)
      CALL cp_fm_release(c_loc_coef)

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(/,T2,A)') &
            '!----------------------------END OF ED ANALYSIS-------------------------------!'
      END IF

      CALL timestop(handle)

   END SUBROUTINE edmf_analysis

END MODULE ed_analysis
