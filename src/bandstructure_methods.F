!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE bandstructure_methods
   USE gw_methods,                      ONLY: gw
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type
   USE qs_environment_types,            ONLY: qs_environment_type
   USE qs_scf,                          ONLY: scf
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'bandstructure_methods'

   PUBLIC :: post_SCF_bandstructure

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief Perform post-SCF band structure calculations from higher level methods
!> \param qs_env  Quickstep environment
!> \par History
!>    * 07.2023 created [Jan Wilhelm]
! **************************************************************************************************
   SUBROUTINE post_SCF_bandstructure(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'post_SCF_bandstructure'

      INTEGER                                            :: handle
      LOGICAL                                            :: do_gw
      TYPE(section_vals_type), POINTER                   :: gw_input_section

      CALL timeset(routineN, handle)

      ! GW calculation
      NULLIFY (gw_input_section)
      gw_input_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%BANDSTRUCTURE%GW")
      CALL section_vals_get(gw_input_section, explicit=do_gw)
      IF (do_gw) THEN
         CALL gw(qs_env, qs_env%gw_env, gw_input_section)
      END IF

      ! TODO: SOC calculation

      ! TODO: local bandgap calculation

      CALL timestop(handle)

   END SUBROUTINE post_SCF_bandstructure

END MODULE bandstructure_methods
