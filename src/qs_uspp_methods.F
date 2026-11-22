!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE qs_uspp_methods
   USE cp_dbcsr_api,                    ONLY: dbcsr_add,&
                                              dbcsr_set,&
                                              dbcsr_type
   USE kinds,                           ONLY: dp
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_uspp_types,                   ONLY: qs_uspp_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: qs_uspp_init_ao_matrix
CONTAINS
! **************************************************************************************************
!> \brief Retain inverse-mediated AO blocks at Gamma. Complete Bloch operators use native group tiles.
!> \param ks_env Prepared inverse allocation pattern, if present.
!> \param matrix New AO matrix, cleared after adding the pattern.
! **************************************************************************************************
   SUBROUTINE qs_uspp_init_ao_matrix(ks_env, matrix)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix

      LOGICAL                                            :: do_kpoints
      TYPE(qs_uspp_type), POINTER                        :: uspp

      CALL get_ks_env(ks_env, uspp=uspp, do_kpoints=do_kpoints)
      IF (.NOT. ASSOCIATED(uspp) .OR. do_kpoints) RETURN
      CPASSERT(ASSOCIATED(uspp%ao_pattern))
      CALL dbcsr_add(matrix, uspp%ao_pattern, 1.0_dp, 1.0_dp)
      CALL dbcsr_set(matrix, 0.0_dp)
   END SUBROUTINE qs_uspp_init_ao_matrix
END MODULE qs_uspp_methods
