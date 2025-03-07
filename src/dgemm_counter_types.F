!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Counters to determine the performance of parallel DGEMMs
!> \par History
!>       2022.05 created [Mauro Del Ben]
!> \author FS, Refactored from mp2_types
! **************************************************************************************************
MODULE dgemm_counter_types
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: m_walltime
   USE message_passing,                 ONLY: mp_para_env_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dgemm_counter_types'

   PUBLIC :: dgemm_counter_type, &
             dgemm_counter_init, &
             dgemm_counter_start, &
             dgemm_counter_stop, &
             dgemm_counter_write

   TYPE dgemm_counter_type
      PRIVATE
      REAL(KIND=dp) :: flop_rate = 0.0_dp, t_start = 0.0_dp
      INTEGER :: num_dgemm_call = 0
      INTEGER :: unit_nr = -1
      LOGICAL :: print_info = .FALSE.
   END TYPE

CONTAINS

! **************************************************************************************************
!> \brief Initialize a dgemm_counter
!> \param dgemm_counter ...
!> \param unit_nr ...
!> \param print_info ...
! **************************************************************************************************
   ELEMENTAL SUBROUTINE dgemm_counter_init(dgemm_counter, unit_nr, print_info)
      TYPE(dgemm_counter_type), INTENT(OUT)              :: dgemm_counter
      INTEGER, INTENT(IN)                                :: unit_nr
      LOGICAL, INTENT(IN)                                :: print_info

      dgemm_counter%unit_nr = unit_nr
      dgemm_counter%print_info = print_info

   END SUBROUTINE dgemm_counter_init

! **************************************************************************************************
!> \brief start timer of the counter
!> \param dgemm_counter ...
! **************************************************************************************************
   SUBROUTINE dgemm_counter_start(dgemm_counter)
      TYPE(dgemm_counter_type), INTENT(INOUT)            :: dgemm_counter

      dgemm_counter%t_start = m_walltime()

   END SUBROUTINE dgemm_counter_start

! **************************************************************************************************
!> \brief stop timer of the counter and provide matrix sizes
!> \param dgemm_counter ...
!> \param size1 ...
!> \param size2 ...
!> \param size3 ...
! **************************************************************************************************
   SUBROUTINE dgemm_counter_stop(dgemm_counter, size1, size2, size3)
      TYPE(dgemm_counter_type), INTENT(INOUT)            :: dgemm_counter
      INTEGER, INTENT(IN)                                :: size1, size2, size3

      REAL(KIND=dp)                                      :: flop_rate, t_end

      t_end = m_walltime()
      flop_rate = 2.0_dp*REAL(size1, dp)*REAL(size2, dp)*REAL(size3, dp)/MAX(0.001_dp, t_end - dgemm_counter%t_start)
      dgemm_counter%num_dgemm_call = dgemm_counter%num_dgemm_call + 1
      IF (dgemm_counter%unit_nr > 0 .AND. dgemm_counter%print_info) THEN
         WRITE (UNIT=dgemm_counter%unit_nr, FMT="(T3,A,I10)") &
            "PERFORMANCE| DGEMM call #", dgemm_counter%num_dgemm_call
         WRITE (UNIT=dgemm_counter%unit_nr, FMT="(T3,A,I12,A,I12,A,I12)") &
            "PERFORMANCE| DGEMM size (M x N x K) = ", size1, " x ", size2, " x ", size3
         WRITE (UNIT=dgemm_counter%unit_nr, FMT="(T3,A,F15.5)") &
            "PERFORMANCE| DGEMM time (s) = ", t_end - dgemm_counter%t_start
      END IF
      dgemm_counter%flop_rate = dgemm_counter%flop_rate + flop_rate

   END SUBROUTINE dgemm_counter_stop

! **************************************************************************************************
!> \brief calculate and print flop rates
!> \param dgemm_counter ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE dgemm_counter_write(dgemm_counter, para_env)
      TYPE(dgemm_counter_type), INTENT(INOUT)            :: dgemm_counter
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env

      dgemm_counter%flop_rate = dgemm_counter%flop_rate/REAL(MAX(dgemm_counter%num_dgemm_call, 1), dp)/1.0E9_dp
      ! Average over all ranks
      CALL para_env%sum(dgemm_counter%flop_rate)
      dgemm_counter%flop_rate = dgemm_counter%flop_rate/(REAL(para_env%num_pe, dp)*REAL(para_env%num_pe, dp))
      IF (dgemm_counter%unit_nr > 0) WRITE (UNIT=dgemm_counter%unit_nr, FMT="(T3,A,T66,F15.2)") &
         "PERFORMANCE| Average DGEMM flop rate (Gflops / MPI rank):", dgemm_counter%flop_rate

   END SUBROUTINE dgemm_counter_write

END MODULE dgemm_counter_types
