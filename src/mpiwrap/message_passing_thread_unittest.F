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
   USE message_passing,                 ONLY: mp_collect_timings,&
                                              mp_comm_type,&
                                              mp_has_thread_multiple,&
                                              mp_world_finalize,&
                                              mp_world_init
   USE timings,                         ONLY: add_timer_env,&
                                              rm_timer_env,&
                                              timings_register_hooks

!$ USE OMP_LIB, ONLY: omp_get_num_threads, omp_get_thread_num, omp_set_dynamic, omp_set_num_threads
#include "../base/base_uses.f90"

   IMPLICIT NONE

   TYPE(mp_comm_type)                                  :: world
#if defined(__parallel)
   INTEGER                                             :: i, nthreads
   TYPE(mp_comm_type), ALLOCATABLE, DIMENSION(:)       :: thread_comm
#endif

   CALL mp_world_init(world)

#if defined(__parallel)
   CPASSERT(mp_has_thread_multiple())

   CALL timings_register_hooks()
   CALL add_timer_env()
   mp_collect_timings = .TRUE.

!$ CALL omp_set_dynamic(.FALSE.)
!$ CALL omp_set_num_threads(2)
   nthreads = 1
!$OMP PARALLEL DEFAULT(NONE) SHARED(nthreads)
!$OMP SINGLE
!$ nthreads = omp_get_num_threads()
!$OMP END SINGLE
!$OMP END PARALLEL

   ALLOCATE (thread_comm(0:nthreads - 1))
   DO i = 0, nthreads - 1
      CALL thread_comm(i)%from_dup(world)
   END DO

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP    SHARED(nthreads, thread_comm, world)
   CALL exercise_thread_comm()
!$OMP END PARALLEL

   DO i = 0, nthreads - 1
      CALL thread_comm(i)%free()
   END DO
   DEALLOCATE (thread_comm)

   mp_collect_timings = .FALSE.
   CALL rm_timer_env()
#endif

   CALL mp_world_finalize()

CONTAINS

#if defined(__parallel)
! **************************************************************************************************
!> \brief Exercise the communicator owned by the calling OpenMP thread.
! **************************************************************************************************
   SUBROUTINE exercise_thread_comm()
      INTEGER, PARAMETER                                 :: num_iterations = 128

      INTEGER                                            :: dest, expected, iter, recv_value, &
                                                            send_value, source, thread_id, value

      thread_id = 0
!$    thread_id = omp_get_thread_num()
      CPASSERT(thread_id >= 0 .AND. thread_id < nthreads)

      dest = MODULO(world%mepos + 1, world%num_pe)
      source = MODULO(world%mepos - 1, world%num_pe)

      DO iter = 1, num_iterations
         value = world%mepos + 1 + 10*thread_id
         CALL thread_comm(thread_id)%sum(value)
         expected = world%num_pe*(world%num_pe + 1)/2 + 10*thread_id*world%num_pe
         CPASSERT(value == expected)

         value = -1
         IF (world%is_source()) value = iter + thread_id
         CALL thread_comm(thread_id)%bcast(value)
         CPASSERT(value == iter + thread_id)

         send_value = (world%num_pe*(iter - 1) + world%mepos)*nthreads + thread_id
         CALL thread_comm(thread_id)%sendrecv(send_value, dest, recv_value, source, tag=iter)
         expected = (world%num_pe*(iter - 1) + source)*nthreads + thread_id
         CPASSERT(recv_value == expected)
      END DO
   END SUBROUTINE exercise_thread_comm
#endif

END PROGRAM message_passing_thread_unittest
