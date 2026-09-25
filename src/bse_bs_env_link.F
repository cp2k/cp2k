!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Bethe-Salpeter equation on the RI-AO GW of PROPERTIES%BANDSTRUCTURE%GW: reads
!>        PROPERTIES%BSE, builds the window-sized RI slabs and the effective inverse dielectric
!>        matrix from bs_env, and hands them to the solver core through its windowed seam.
!>        Non-periodic, Gamma point only. Reads no real-time propagation input.
!> \author Max Graml
! **************************************************************************************************
MODULE bse_bs_env_link
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE input_section_types,             ONLY: section_vals_type
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'bse_bs_env_link'

   PUBLIC :: bse_from_bs_env

CONTAINS

! **************************************************************************************************
!> \brief Entry from the properties driver after the bandstructure block
!> \param qs_env ...
!> \param bse_section PROPERTIES%BSE
! **************************************************************************************************
   SUBROUTINE bse_from_bs_env(qs_env, bse_section)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: bse_section

      CHARACTER(LEN=*), PARAMETER :: routineN = 'bse_from_bs_env'

      INTEGER                                            :: handle, unit_nr
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      unit_nr = cp_logger_get_default_io_unit(logger)
      NULLIFY (bs_env)
      CALL get_qs_env(qs_env, bs_env=bs_env)
      CPASSERT(ASSOCIATED(bs_env))
      ! the requirement reader at GW init refused RI_RS and periodic runs already
      CPASSERT(bs_env%do_gw .AND. .NOT. bs_env%do_periodic .AND. .NOT. bs_env%do_gw_ri_rs)
      CPASSERT(bs_env%keep_W_freq_zero .AND. bs_env%keep_3c_nl)

      CPABORT("PROPERTIES%BSE: adapter not yet available (lands with the next commit)")
      MARK_USED(bse_section)
      MARK_USED(unit_nr)

      CALL timestop(handle)

   END SUBROUTINE bse_from_bs_env

END MODULE bse_bs_env_link
