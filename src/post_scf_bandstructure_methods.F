!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE post_scf_bandstructure_methods
   USE gw_main,                         ONLY: gw
   USE input_section_types,             ONLY: section_vals_type
   USE post_scf_bandstructure_utils,    ONLY: create_and_init_bs_env,&
                                              dos_pdos_ldos,&
                                              soc
   USE qs_environment_types,            ONLY: qs_environment_type
   USE qs_scf,                          ONLY: scf
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'post_scf_bandstructure_methods'

   PUBLIC :: post_scf_bandstructure

CONTAINS

! **************************************************************************************************
!> \brief Perform post-SCF band structure calculations from higher level methods
!> \param qs_env  Quickstep environment
!> \param post_scf_bandstructure_section ...
!> \par History
!>    * 07.2023 created [Jan Wilhelm]
! **************************************************************************************************
   SUBROUTINE post_scf_bandstructure(qs_env, post_scf_bandstructure_section)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: post_scf_bandstructure_section

      CHARACTER(LEN=*), PARAMETER :: routineN = 'post_scf_bandstructure'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! general setup of post SCF bandstructure calculation
      CALL create_and_init_bs_env(qs_env, qs_env%bs_env, post_scf_bandstructure_section)

      ! shifts of eigenvalues/bandstructure due to spin-orbit coupling from pseudopotentials
      IF (qs_env%bs_env%do_soc) THEN
         CALL soc(qs_env, qs_env%bs_env)
      END IF

      ! GW calculation for eigenvalues/bandstructure for molecules and periodic systems
      IF (qs_env%bs_env%do_gw) THEN
         CALL gw(qs_env, qs_env%bs_env, post_scf_bandstructure_section)
      END IF

      ! density of states (DOS), projected DOS, local DOS for DFT, DFT+SOC, G0W0, G0W0+SOC
      CALL dos_pdos_ldos(qs_env, qs_env%bs_env)

      CALL timestop(handle)

   END SUBROUTINE post_scf_bandstructure

END MODULE post_scf_bandstructure_methods
