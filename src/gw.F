!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief
!> \author Jan Wilhelm
!> \date 07.2023
! **************************************************************************************************
MODULE gw
   USE gw_ri_rs_large_cell_Gamma,       ONLY: gw_calc_ri_rs_large_cell_Gamma
   USE gw_ri_rs_non_periodic,           ONLY: gw_calc_ri_rs_non_periodic
   USE gw_tensor_large_cell_Gamma,      ONLY: gw_calc_tensor_large_cell_Gamma
   USE gw_tensor_small_cell_full_kp,    ONLY: gw_calc_tensor_small_cell_full_kp
   USE gw_utils,                        ONLY: create_and_init_bs_env_for_gw
   USE input_constants,                 ONLY: ri_rs_large_cell_Gamma,&
                                              ri_rs_non_periodic,&
                                              tensor_large_cell_Gamma,&
                                              tensor_small_cell_full_kp
   USE input_section_types,             ONLY: section_vals_type
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE qs_environment_types,            ONLY: qs_environment_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw'

   PUBLIC :: run_gw

CONTAINS

! **************************************************************************************************
!> \brief Perform GW band structure calculation
!> \param qs_env ...
!> \param bs_env ...
!> \param post_scf_bandstructure_section ...
!> \par History
!>    * 04.2026 Added GW RI-RS module [Ritaj Tyagi]
!>    * 07.2023 created [Jan Wilhelm]
! **************************************************************************************************
   SUBROUTINE run_gw(qs_env, bs_env, post_scf_bandstructure_section)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(section_vals_type), POINTER                   :: post_scf_bandstructure_section

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'run_gw'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL create_and_init_bs_env_for_gw(qs_env, bs_env, post_scf_bandstructure_section)

      SELECT CASE (bs_env%gw_implementation)

      CASE (tensor_small_cell_full_kp)

         CALL gw_calc_tensor_small_cell_full_kp(qs_env, bs_env)

      CASE (tensor_large_cell_Gamma)

         CALL gw_calc_tensor_large_cell_Gamma(qs_env, bs_env)

      CASE (ri_rs_large_cell_Gamma)

         CALL gw_calc_ri_rs_large_cell_Gamma(qs_env, bs_env)

      CASE (ri_rs_non_periodic)

         CALL gw_calc_ri_rs_non_periodic(qs_env, bs_env)

      CASE DEFAULT

         CPABORT("Unknown GW implementation.")

      END SELECT

      CALL timestop(handle)

   END SUBROUTINE run_gw

END MODULE gw
