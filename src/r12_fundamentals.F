!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Shared integral tables and external radial fundamentals for TShPSC.
!>        Standard Libint operators generate their fundamentals inside Engine.
! **************************************************************************************************
MODULE r12_fundamentals
   USE cp_files,                        ONLY: close_file,&
                                              discover_file,&
                                              open_file
   USE gamma,                           ONLY: deallocate_md_ftable,&
                                              init_md_ftable
   USE input_constants,                 ONLY: do_potential_TShPSC,&
                                              do_potential_mix_cl_trunc,&
                                              do_potential_truncated
   USE integral_library_types,          ONLY: coulomb_operator_type
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE mathconstants,                   ONLY: pi
   USE message_passing,                 ONLY: mp_para_env_type
   USE t_c_g0,                          ONLY: free_truncated_table => free_C0,&
                                              get_lmax_init,&
                                              init_truncated_table => init
   USE t_sh_p_s_c,                      ONLY: free_tshpsc_table => free_C0,&
                                              init_tshpsc_table => init,&
                                              trunc_CS_poly_n20
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: compute_tshpsc_fundamentals, r12_fundamentals_init, r12_fundamentals_cleanup

   INTEGER :: tshpsc_lmax = -1
   CHARACTER(default_path_length) :: truncated_filename = "", tshpsc_filename = ""

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
      LOGICAL                                            :: truncated, tshpsc

      CALL init_md_ftable(m_max)
      truncated = .FALSE.
      tshpsc = .FALSE.
      SELECT CASE (potential_parameter%potential_type)
      CASE (do_potential_truncated, do_potential_mix_cl_trunc)
         truncated = m_max > get_lmax_init() .OR. potential_parameter%filename /= truncated_filename
      CASE (do_potential_TShPSC)
         tshpsc = m_max > tshpsc_lmax .OR. potential_parameter%filename /= tshpsc_filename
      END SELECT
      IF (.NOT. truncated .AND. .NOT. tshpsc) RETURN

      unit_id = -1
      IF (para_env%is_source()) THEN
         CALL open_file(unit_number=unit_id, file_name=discover_file(potential_parameter%filename))
      END IF
      IF (truncated) THEN
         table_order = MAX(m_max, get_lmax_init())
         CALL init_truncated_table(table_order, unit_id, para_env%mepos, para_env)
         truncated_filename = potential_parameter%filename
      ELSE
         tshpsc_lmax = MAX(m_max, tshpsc_lmax)
         CALL init_tshpsc_table(tshpsc_lmax, unit_id, para_env%mepos, para_env)
         tshpsc_filename = potential_parameter%filename
      END IF
      IF (para_env%is_source()) CALL close_file(unit_id)
   END SUBROUTINE r12_fundamentals_init

! **************************************************************************************************
!> \brief Release shared radial evaluator tables at process finalization, after their consumers.
! **************************************************************************************************
   SUBROUTINE r12_fundamentals_cleanup()
      CALL free_truncated_table()
      CALL free_tshpsc_table()
      CALL deallocate_md_ftable()
      tshpsc_lmax = -1
      truncated_filename = ""
      tshpsc_filename = ""
   END SUBROUTINE r12_fundamentals_cleanup

! **************************************************************************************************
!> \brief Compute external TShPSC radial seeds, retaining the legacy evaluation order.
!> \param F coefficients of orders zero through m_max, excluding normalization
!> \param normalization common multiplicative normalization
!> \param m_max maximum radial order
!> \param rho reduced Gaussian exponent
!> \param T reduced squared product-center distance
!> \param cutoff_radius TShPSC cutoff radius
! **************************************************************************************************
   SUBROUTINE compute_tshpsc_fundamentals(F, normalization, m_max, rho, T, cutoff_radius)
      REAL(dp), INTENT(OUT)                              :: normalization
      INTEGER, INTENT(IN)                                :: m_max
      REAL(dp), INTENT(OUT)                              :: F(m_max + 1)
      REAL(dp), INTENT(IN)                               :: rho, T, cutoff_radius

      REAL(dp)                                           :: R

      normalization = 2.0_dp*pi/rho
      R = cutoff_radius*SQRT(rho)
      F = 0.0_dp
      CALL trunc_CS_poly_n20(F, R, T, m_max)
   END SUBROUTINE compute_tshpsc_fundamentals

END MODULE r12_fundamentals
