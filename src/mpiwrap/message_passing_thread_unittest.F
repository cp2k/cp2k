!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Exercises MPI wrappers concurrently from multiple OpenMP threads.
! **************************************************************************************************
PROGRAM message_passing_thread_unittest
   USE OMP_LIB,                         ONLY: omp_get_num_threads,&
                                              omp_get_thread_num,&
                                              omp_set_dynamic,&
                                              omp_set_num_threads
   USE message_passing,                 ONLY: mp_collect_timings,&
                                              mp_comm_type,&
                                              mp_thread_comm_pool_type,&
                                              mp_world_finalize,&
                                              mp_world_init
   USE timings,                         ONLY: add_timer_env,&
                                              rm_timer_env,&
                                              timings_register_hooks
#include "../base/base_uses.f90"

   IMPLICIT NONE

   TYPE(mp_comm_type)                                  :: world
   TYPE(mp_comm_type), POINTER                         :: thread_comm => NULL()
   TYPE(mp_thread_comm_pool_type)                      :: thread_comms
   INTEGER                                             :: nthreads

   CALL mp_world_init(world)

   CALL timings_register_hooks()
   CALL add_timer_env()
   mp_collect_timings = .TRUE.

   CALL omp_set_dynamic(.FALSE.)
   CALL omp_set_num_threads(2)
   nthreads = 0
!$OMP PARALLEL DEFAULT(NONE) SHARED(nthreads)
!$OMP SINGLE
   nthreads = omp_get_num_threads()
!$OMP END SINGLE
!$OMP END PARALLEL
   CALL thread_comms%init(world, num_lanes=nthreads)
   CPASSERT(thread_comms%is_initialized())
   CPASSERT(thread_comms%get_num_lanes() == nthreads)

!$OMP PARALLEL DEFAULT(NONE) NUM_THREADS(nthreads) &
!$OMP    PRIVATE(thread_comm) SHARED(thread_comms, world)
   CALL thread_comms%get_comm(thread_comm)
   CALL exercise_thread_comm(thread_comm, world)
!$OMP END PARALLEL

   CALL thread_comms%free()
   CPASSERT(.NOT. thread_comms%is_initialized())

   mp_collect_timings = .FALSE.
   CALL rm_timer_env()
   CALL mp_world_finalize()

CONTAINS

! **************************************************************************************************
!> \brief Exercise the communicator owned by the calling OpenMP thread.
!> \param thread_comm ...
!> \param world ...
! **************************************************************************************************
   SUBROUTINE exercise_thread_comm(thread_comm, world)
      TYPE(mp_comm_type), INTENT(IN)                     :: thread_comm, world

      INTEGER, PARAMETER                                 :: num_iterations = 128

      INTEGER                                            :: dest, expected, iter, nthreads, &
                                                            recv_value, send_value, source, &
                                                            thread_id, value

      thread_id = omp_get_thread_num()
      nthreads = omp_get_num_threads()
      CPASSERT(thread_id >= 0 .AND. thread_id < nthreads)

      dest = MODULO(world%mepos + 1, world%num_pe)
      source = MODULO(world%mepos - 1, world%num_pe)

      DO iter = 1, num_iterations
         value = world%mepos + 1 + 10*thread_id
         CALL thread_comm%sum(value)
         expected = world%num_pe*(world%num_pe + 1)/2 + 10*thread_id*world%num_pe
         CPASSERT(value == expected)

         value = -1
         IF (world%is_source()) value = iter + thread_id
         CALL thread_comm%bcast(value)
         CPASSERT(value == iter + thread_id)

         send_value = (world%num_pe*(iter - 1) + world%mepos)*nthreads + thread_id
         CALL thread_comm%sendrecv(send_value, dest, recv_value, source, tag=iter)
         expected = (world%num_pe*(iter - 1) + source)*nthreads + thread_id
         CPASSERT(recv_value == expected)
      END DO
   END SUBROUTINE exercise_thread_comm

END PROGRAM message_passing_thread_unittest
