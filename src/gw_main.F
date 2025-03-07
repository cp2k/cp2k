!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief
!> \author Jan Wilhelm
!> \date 07.2023
! **************************************************************************************************
MODULE gw_main
   USE gw_large_cell_Gamma,             ONLY: gw_calc_large_cell_Gamma
   USE gw_small_cell_full_kp,           ONLY: gw_calc_small_cell_full_kp
   USE gw_utils,                        ONLY: create_and_init_bs_env_for_gw
   USE input_constants,                 ONLY: large_cell_Gamma,&
                                              small_cell_full_kp
   USE input_section_types,             ONLY: section_vals_type
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE qs_environment_types,            ONLY: qs_environment_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_main'

   PUBLIC :: gw

CONTAINS

! **************************************************************************************************
!> \brief Perform GW band structure calculation
!> \param qs_env ...
!> \param bs_env ...
!> \param post_scf_bandstructure_section ...
!> \par History
!>    * 07.2023 created [Jan Wilhelm]
! **************************************************************************************************
   SUBROUTINE gw(qs_env, bs_env, post_scf_bandstructure_section)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(section_vals_type), POINTER                   :: post_scf_bandstructure_section

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'gw'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL create_and_init_bs_env_for_gw(qs_env, bs_env, post_scf_bandstructure_section)

      SELECT CASE (bs_env%small_cell_full_kp_or_large_cell_Gamma)

      CASE (small_cell_full_kp)

         CALL gw_calc_small_cell_full_kp(qs_env, bs_env)

      CASE (large_cell_Gamma)

         CALL gw_calc_large_cell_Gamma(qs_env, bs_env)

      END SELECT

      CALL timestop(handle)

   END SUBROUTINE gw

END MODULE gw_main

