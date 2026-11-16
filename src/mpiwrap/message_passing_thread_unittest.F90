!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Exercises point-to-point MPI wrappers concurrently from multiple OpenMP threads.
! **************************************************************************************************
PROGRAM message_passing_thread_unittest
   USE OMP_LIB,                         ONLY: omp_get_num_threads,&
                                              omp_get_thread_num,&
                                              omp_set_dynamic,&
                                              omp_set_num_threads
   USE message_passing,                 ONLY: mp_collect_timings,&
                                              mp_comm_type,&
                                              mp_world_finalize,&
                                              mp_world_init
   USE timings,                         ONLY: add_timer_env,&
                                              rm_timer_env,&
                                              timings_register_hooks
#include "../base/base_uses.f90"

   IMPLICIT NONE

   TYPE(mp_comm_type)                                  :: world

   CALL mp_world_init(world)

   CALL timings_register_hooks()
   CALL add_timer_env()
   mp_collect_timings = .TRUE.

   CALL omp_set_dynamic(.FALSE.)
   CALL omp_set_num_threads(2)
!$OMP PARALLEL DEFAULT(NONE) SHARED(world)
   CALL exercise_concurrent_point_to_point(world)
!$OMP END PARALLEL

   mp_collect_timings = .FALSE.
   CALL rm_timer_env()
   CALL mp_world_finalize()

CONTAINS

! **************************************************************************************************
!> \brief Exercise concurrent point-to-point communication on a shared communicator.
!> \param comm communicator shared by all OpenMP threads
! **************************************************************************************************
   SUBROUTINE exercise_concurrent_point_to_point(comm)
      TYPE(mp_comm_type), INTENT(IN)                     :: comm

      INTEGER, PARAMETER                                 :: num_iterations = 128

      INTEGER                                            :: dest, expected, iter, nthreads, &
                                                            recv_value, send_value, source, tag, &
                                                            thread_id

      thread_id = omp_get_thread_num()
      nthreads = omp_get_num_threads()
      CPASSERT(nthreads == 2)
      CPASSERT(thread_id >= 0 .AND. thread_id < nthreads)

      dest = MODULO(comm%mepos + 1, comm%num_pe)
      source = MODULO(comm%mepos - 1, comm%num_pe)

      DO iter = 1, num_iterations
         tag = (iter - 1)*nthreads + thread_id
         send_value = (comm%num_pe*(iter - 1) + comm%mepos)*nthreads + thread_id
         CALL comm%sendrecv(send_value, dest, recv_value, source, tag=tag)
         expected = (comm%num_pe*(iter - 1) + source)*nthreads + thread_id
         CPASSERT(recv_value == expected)
      END DO
   END SUBROUTINE exercise_concurrent_point_to_point

END PROGRAM message_passing_thread_unittest
