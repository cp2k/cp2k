!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Shared radial evaluator tables for native integral kernels.
! **************************************************************************************************
MODULE r12_fundamentals
   USE cp_files,                        ONLY: close_file,&
                                              discover_file,&
                                              open_file
   USE gamma,                           ONLY: deallocate_md_ftable,&
                                              init_md_ftable
   USE input_constants,                 ONLY: do_potential_mix_cl_trunc,&
                                              do_potential_truncated
   USE integral_library_types,          ONLY: coulomb_operator_type
   USE kinds,                           ONLY: default_path_length
   USE message_passing,                 ONLY: mp_para_env_type
   USE t_c_g0,                          ONLY: free_truncated_table => free_C0,&
                                              get_lmax_init,&
                                              init_truncated_table => init
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: r12_fundamentals_init, r12_fundamentals_cleanup

   CHARACTER(default_path_length) :: truncated_filename = ""

CONTAINS

! **************************************************************************************************
!> \brief Prepare shared radial evaluator tables before entering a threaded evaluation region.
!> \param m_max maximum radial order including any derivative increment
!> \param potential_parameter physical radial operator and table filename
!> \param para_env communicator; table initialization is collective on this communicator
!> \note Call outside OpenMP regions. Per-thread integral contexts borrow these read-only tables.
! **************************************************************************************************
   SUBROUTINE r12_fundamentals_init(m_max, potential_parameter, para_env)
      INTEGER, INTENT(IN)                                :: m_max
      TYPE(coulomb_operator_type), INTENT(IN)            :: potential_parameter
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env

      INTEGER                                            :: table_order, unit_id
      LOGICAL                                            :: truncated

      CALL init_md_ftable(m_max)
      truncated = .FALSE.
      SELECT CASE (potential_parameter%potential_type)
      CASE (do_potential_truncated, do_potential_mix_cl_trunc)
         truncated = m_max > get_lmax_init() .OR. potential_parameter%filename /= truncated_filename
      END SELECT
      IF (.NOT. truncated) RETURN

      unit_id = -1
      IF (para_env%is_source()) THEN
         CALL open_file(unit_number=unit_id, file_name=discover_file(potential_parameter%filename))
      END IF
      table_order = MAX(m_max, get_lmax_init())
      CALL init_truncated_table(table_order, unit_id, para_env%mepos, para_env)
      truncated_filename = potential_parameter%filename
      IF (para_env%is_source()) CALL close_file(unit_id)
   END SUBROUTINE r12_fundamentals_init

! **************************************************************************************************
!> \brief Release shared radial evaluator tables at process finalization, after their consumers.
! **************************************************************************************************
   SUBROUTINE r12_fundamentals_cleanup()
      CALL free_truncated_table()
      CALL deallocate_md_ftable()
      truncated_filename = ""
   END SUBROUTINE r12_fundamentals_cleanup

END MODULE r12_fundamentals
