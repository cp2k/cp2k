!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Interface to the message passing library MPI
!> \par History
!>      JGH (02-Jan-2001): New error handling
!>                         Performance tools
!>      JGH (14-Jan-2001): New routines mp_comm_compare, mp_cart_coords,
!>                                      mp_rank_compare, mp_alltoall
!>      JGH (06-Feb-2001): New routines mp_comm_free
!>      JGH (22-Mar-2001): New routines mp_comm_dup
!>      fawzi (04-NOV-2004): storable performance info (for f77 interface)
!>      Wrapper routine for mpi_gatherv added (22.12.2005,MK)
!>      JGH (13-Feb-2006): Flexible precision
!>      JGH (15-Feb-2006): single precision mp_alltoall
!> \author JGH
! **************************************************************************************************
MODULE mp_perf_test
   USE kinds, ONLY: dp
   USE message_passing, ONLY: mp_comm_type
   ! some benchmarking code
#include "../base/base_uses.f90"

#if defined(__parallel)
#if defined(__MPI_F08)
   USE mpi_f08, ONLY: mpi_wtime
#else
   USE mpi, ONLY: mpi_wtime
#endif
#endif

   PRIVATE

   PUBLIC :: mpi_perf_test

CONTAINS

! **************************************************************************************************
!> \brief Tests the MPI library
!> \param comm the relevant, initialized communicator
!> \param npow number of sizes to test, 10**1 .. 10**npow
!> \param output_unit where to direct output
!> \par History
!>      JGH  6-Feb-2001 : Test and performance code
!> \author JGH  1-JAN-2001
!> \note
!>      quickly adapted benchmark code, will only work on an even number of CPUs.
! **************************************************************************************************
   SUBROUTINE mpi_perf_test(comm, npow, output_unit)
      CLASS(mp_comm_type), INTENT(IN) :: comm
      INTEGER, INTENT(IN)                      :: npow, output_unit

#if defined(__parallel)

      INTEGER :: I, itask, itests, J, jtask, left, nbufmax, &
                 ncount, Ngrid, Nloc, nprocs, Ntot, partner, right, taskid, tag, source
      INTEGER, ALLOCATABLE, DIMENSION(:)       :: rcount, rdispl, scount, sdispl
      LOGICAL                                  :: ionode
      REAL(KIND=dp)                            :: maxdiff, t1, &
                                                  t2, t3, t4, t5
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: buffer1, buffer2, buffer3, &
                                                  lgrid, lgrid2, lgrid3
      REAL(KIND=dp), ALLOCATABLE, &
         DIMENSION(:, :)                        :: grid, grid2, grid3, &
                                                   send_timings, send_timings2
      REAL(KIND=dp), PARAMETER :: threshold = 1.0E-8_dp

      ! set system sizes !
      ngrid = 10**npow

      taskid = comm%mepos
      nprocs = comm%num_pe
      ionode = comm%is_source()
      IF (ionode .AND. output_unit > 0) THEN
         WRITE (output_unit, *) "Running with ", nprocs
         WRITE (output_unit, *) "running messages with npow = ", npow
         WRITE (output_unit, *) "use MPI X in the input for larger (e.g. 6) of smaller (e.g. 3) messages"
         IF (MODULO(nprocs, 2) .NE. 0) WRITE (output_unit, *) "Testing only with an even number of tasks"
      END IF

      IF (MODULO(nprocs, 2) .NE. 0) RETURN

      ! equal loads
      Nloc = Ngrid/nprocs
      Ntot = Nprocs*Nloc
      nbufmax = 10**npow
      !
      ALLOCATE (rcount(nprocs))
      ALLOCATE (scount(nprocs))
      ALLOCATE (sdispl(nprocs))
      ALLOCATE (rdispl(nprocs))
      ALLOCATE (buffer1(nbufmax))
      ALLOCATE (buffer2(nbufmax))
      ALLOCATE (buffer3(nbufmax))
      ALLOCATE (grid(Nloc, Nprocs))
      ALLOCATE (grid2(Nloc, Nprocs))
      ALLOCATE (grid3(Nloc, Nprocs))
      ALLOCATE (lgrid(Nloc))
      ALLOCATE (lgrid2(Nloc))
      ALLOCATE (lgrid3(Nloc))
      ALLOCATE (send_timings(0:nprocs - 1, 0:nprocs - 1))
      ALLOCATE (send_timings2(0:nprocs - 1, 0:nprocs - 1))
      buffer1 = 0.0_dp
      buffer2 = 0.0_dp
      buffer3 = 0.0_dp
      ! timings
      send_timings = 0.0_dp
      send_timings2 = 0.0_dp
      ! -------------------------------------------------------------------------------------------
      ! ------------------------------ some in memory tests                   ---------------------
      ! -------------------------------------------------------------------------------------------
      CALL comm%sync()
      IF (ionode .AND. output_unit > 0) THEN
         WRITE (output_unit, *) "Testing in memory copies just 1 CPU "
         WRITE (output_unit, *) "  could tell something about the motherboard / cache / compiler "
      END IF
      DO i = 1, npow
         ncount = 10**i
         t2 = 0.0E0_dp
         IF (ncount .GT. nbufmax) CPABORT("")
         DO j = 1, 3**(npow - i)
            CALL comm%sync()
            t1 = MPI_WTIME()
            buffer2(1:ncount) = buffer1(1:ncount)
            t2 = t2 + MPI_WTIME() - t1 + threshold
         END DO
         CALL comm%max(t2, 0)
         IF (ionode .AND. output_unit > 0) THEN
            WRITE (output_unit, '(I9,A,F12.4,A)') 8*ncount, " Bytes ", (3**(npow - i))*ncount*8.0E-6_dp/t2, " MB/s"
         END IF
      END DO
      ! -------------------------------------------------------------------------------------------
      ! ------------------------------ some in memory tests                   ---------------------
      ! -------------------------------------------------------------------------------------------
      CALL comm%sync()
      IF (ionode .AND. output_unit > 0) THEN
         WRITE (output_unit, *) "Testing in memory copies all cpus"
         WRITE (output_unit, *) "  is the memory bandwidth affected on an SMP machine ?"
      END IF
      DO i = 1, npow
         ncount = 10**i
         t2 = 0.0E0_dp
         IF (ncount .GT. nbufmax) CPABORT("")
         DO j = 1, 3**(npow - i)
            CALL comm%sync()
            t1 = MPI_WTIME()
            buffer2(1:ncount) = buffer1(1:ncount)
            t2 = t2 + MPI_WTIME() - t1 + threshold
         END DO
         CALL comm%max(t2, 0)
         IF (ionode .AND. output_unit > 0) THEN
            WRITE (output_unit, '(I9,A,F12.4,A)') 8*ncount, " Bytes ", (3**(npow - i))*ncount*8.0E-6_dp/t2, " MB/s"
         END IF
      END DO
      ! -------------------------------------------------------------------------------------------
      ! ------------------------------ first test point to point communication ---------------------
      ! -------------------------------------------------------------------------------------------
      CALL comm%sync()
      IF (ionode .AND. output_unit > 0) THEN
         WRITE (output_unit, *) "Testing truly point to point communication (i with j only)"
         WRITE (output_unit, *) "  is there some different connection between i j (e.g. shared memory comm)"
      END IF
      ncount = 10**npow
      IF (ionode .AND. output_unit > 0) WRITE (output_unit, *) "For messages of ", ncount*8, " bytes"
      IF (ncount .GT. nbufmax) CPABORT("")
      DO itask = 0, nprocs - 1
         DO jtask = itask + 1, nprocs - 1
            CALL comm%sync()
            t1 = MPI_WTIME()
            IF (taskid .EQ. itask) THEN
               CALL comm%send(buffer1, jtask, itask*jtask)
            END IF
            IF (taskid .EQ. jtask) THEN
               source = itask
               tag = itask*jtask
               CALL comm%recv(buffer1, source, tag)
            END IF
            send_timings(itask, jtask) = MPI_WTIME() - t1 + threshold
         END DO
      END DO
      CALL comm%max(send_timings, 0)
      IF (ionode .AND. output_unit > 0) THEN
         DO itask = 0, nprocs - 1
            DO jtask = itask + 1, nprocs - 1
               WRITE (output_unit, '(I4,I4,F12.4,A)') itask, jtask, ncount*8.0E-6_dp/send_timings(itask, jtask), " MB/s"
            END DO
         END DO
      END IF
      CALL comm%sync()
      ! -------------------------------------------------------------------------------------------
      ! ------------------------------ second test point to point communication -------------------
      ! -------------------------------------------------------------------------------------------
      IF (ionode .AND. output_unit > 0) THEN
         WRITE (output_unit, *) "Testing all nearby point to point communication (0,1)(2,3)..."
         WRITE (output_unit, *) "    these could / should all be on the same shared memory node "
      END IF
      DO i = 1, npow
         ncount = 10**i
         t2 = 0.0E0_dp
         IF (ncount .GT. nbufmax) CPABORT("")
         DO j = 1, 3**(npow - i)
            CALL comm%sync()
            t1 = MPI_WTIME()
            IF (MODULO(taskid, 2) == 0) THEN
               CALL comm%send(buffer1, taskid + 1, 0)
            ELSE
               source = taskid - 1
               tag = 0
               CALL comm%recv(buffer1, source, tag)
            END IF
            t2 = t2 + MPI_WTIME() - t1 + threshold
         END DO
         CALL comm%max(t2, 0)
         IF (ionode .AND. output_unit > 0) THEN
            WRITE (output_unit, '(I9,A,F12.4,A)') 8*ncount, " Bytes ", (3**(npow - i))*ncount*8.0E-6_dp/t2, " MB/s"
         END IF
      END DO
      CALL comm%sync()
      ! -------------------------------------------------------------------------------------------
      ! ------------------------------ third test point to point communication -------------------
      ! -------------------------------------------------------------------------------------------
      IF (ionode .AND. output_unit > 0) THEN
         WRITE (output_unit, *) "Testing all far point to point communication (0,nprocs/2),(1,nprocs/2+1),.."
         WRITE (output_unit, *) "    these could all be going over the network, and stress it a lot"
      END IF
      DO i = 1, npow
         ncount = 10**i
         t2 = 0.0E0_dp
         IF (ncount .GT. nbufmax) CPABORT("")
         DO j = 1, 3**(npow - i)
            CALL comm%sync()
            t1 = MPI_WTIME()
            ! first half with partner
            IF (taskid .LT. nprocs/2) THEN
               CALL comm%send(buffer1, taskid + nprocs/2, 0)
            ELSE
               source = taskid - nprocs/2
               tag = 0
               CALL comm%recv(buffer1, source, tag)
            END IF
            t2 = t2 + MPI_WTIME() - t1 + threshold
         END DO
         CALL comm%max(t2, 0)
         IF (ionode .AND. output_unit > 0) THEN
            WRITE (output_unit, '(I9,A,F12.4,A)') 8*ncount, " Bytes ", (3**(npow - i))*ncount*8.0E-6_dp/t2, " MB/s"
         END IF
      END DO
      ! -------------------------------------------------------------------------------------------
      ! ------------------------------ test root to all broadcast               -------------------
      ! -------------------------------------------------------------------------------------------
      CALL comm%sync()
      IF (ionode .AND. output_unit > 0) THEN
         WRITE (output_unit, *) "Testing root to all broadcast "
         WRITE (output_unit, *) "    using trees at least ? "
      END IF
      DO i = 1, npow
         ncount = 10**i
         t2 = 0.0E0_dp
         IF (ncount .GT. nbufmax) CPABORT("")
         DO j = 1, 3**(npow - i)
            CALL comm%sync()
            t1 = MPI_WTIME()
            CALL comm%bcast(buffer1, 0)
            t2 = t2 + MPI_WTIME() - t1 + threshold
         END DO
         CALL comm%max(t2, 0)
         IF (ionode .AND. output_unit > 0) THEN
            WRITE (output_unit, '(I9,A,F12.4,A)') 8*ncount, " Bytes ", (3**(npow - i))*ncount*8.0E-6_dp/t2, " MB/s"
         END IF
      END DO
      ! -------------------------------------------------------------------------------------------
      ! ------------------------------ test parallel sum like behavior                -------------------
      ! -------------------------------------------------------------------------------------------
      CALL comm%sync()
      IF (ionode .AND. output_unit > 0) WRITE (output_unit, *) "Test global summation (mpi_allreduce) "
      DO i = 1, npow
         ncount = 10**i
         t2 = 0.0E0_dp
         IF (ncount .GT. nbufmax) CPABORT("")
         DO j = 1, 3**(npow - i)
            buffer2(:) = buffer1
            CALL comm%sync()
            t1 = MPI_WTIME()
            CALL comm%sum(buffer2)
            t2 = t2 + MPI_WTIME() - t1 + threshold
         END DO
         CALL comm%max(t2, 0)
         IF (ionode .AND. output_unit > 0) THEN
            WRITE (output_unit, '(I9,A,F12.4,A)') 8*ncount, " Bytes ", (3**(npow - i))*ncount*8.0E-6_dp/t2, " MB/s"
         END IF
      END DO
      ! -------------------------------------------------------------------------------------------
      ! ------------------------------ test all to all communication            -------------------
      ! -------------------------------------------------------------------------------------------
      CALL comm%sync()
      IF (ionode .AND. output_unit > 0) THEN
         WRITE (output_unit, *) "Test all to all communication (mpi_alltoallv)"
         WRITE (output_unit, *) "    mpi/network getting confused ? "
      END IF
      DO i = 1, npow
         ncount = 10**i
         t2 = 0.0E0_dp
         IF (ncount .GT. nbufmax) CPABORT("")
         scount = ncount/nprocs
         rcount = ncount/nprocs
         DO j = 1, nprocs
            sdispl(j) = (j - 1)*(ncount/nprocs)
            rdispl(j) = (j - 1)*(ncount/nprocs)
         END DO
         DO j = 1, 3**(npow - i)
            CALL comm%sync()
            t1 = MPI_WTIME()
            CALL comm%alltoall(buffer1, scount, sdispl, buffer2, rcount, rdispl)
            t2 = t2 + MPI_WTIME() - t1 + threshold
         END DO
         CALL comm%max(t2, 0)
         IF (ionode .AND. output_unit > 0) THEN
            WRITE (output_unit, '(I9,A,F12.4,A)') 8*(ncount/nprocs)*nprocs, " Bytes ", &
               (3**(npow - i))*(ncount/nprocs)*nprocs*8.0E-6_dp/t2, " MB/s"
         END IF
      END DO

      ! -------------------------------------------------------------------------------------------
      ! ------------------------------ other stuff                            ---------------------
      ! -------------------------------------------------------------------------------------------
      IF (ionode .AND. output_unit > 0) THEN
         WRITE (output_unit, *) " Clean tests completed "
         WRITE (output_unit, *) " Testing MPI_REDUCE scatter"
      END IF
      rcount = Nloc
      DO itests = 1, 3
         IF (ionode .AND. output_unit > 0) &
            WRITE (output_unit, *) "------------------------------- test ", itests, " ------------------------"
         ! *** reference ***
         DO j = 1, Nprocs
            DO i = 1, Nloc
               grid(i, j) = MODULO(i*j*taskid, itests)
            END DO
         END DO
         t1 = MPI_WTIME()
         CALL comm%mp_sum_scatter_dv(grid, lgrid, rcount)
         t2 = MPI_WTIME() - t1 + threshold
         CALL comm%max(t2)
         IF (ionode .AND. output_unit > 0) WRITE (output_unit, *) "MPI_REDUCE_SCATTER    ", t2
         ! *** simple shift ***
         DO j = 1, Nprocs
            DO i = 1, Nloc
               grid2(i, j) = MODULO(i*j*taskid, itests)
            END DO
         END DO
         t3 = MPI_WTIME()
         lgrid2 = 0.0E0_dp
         DO i = 1, Nprocs
            lgrid2(:) = lgrid2 + grid(:, MODULO(taskid - i, Nprocs) + 1)
            IF (i .EQ. nprocs) EXIT
            CALL comm%shift(lgrid2, 1)
         END DO
         t4 = MPI_WTIME() - t3 + threshold
         CALL comm%max(t4)
         maxdiff = MAXVAL(ABS(lgrid2 - lgrid))
         CALL comm%max(maxdiff)
         IF (ionode .AND. output_unit > 0) WRITE (output_unit, *) "MPI_SENDRECV_REPLACE  ", t4, maxdiff
         ! *** involved shift ****
         IF (MODULO(nprocs, 2) /= 0) CPABORT("")
         DO j = 1, Nprocs
            DO i = 1, Nloc
               grid3(i, j) = MODULO(i*j*taskid, itests)
            END DO
         END DO
         t3 = MPI_WTIME()
         ! first sum the grid in pairs (0,1),(2,3) should be within an LPAR and fast XXXXXXXXX
         ! 0 will only need parts 0,2,4,... correctly summed
         ! 1 will only need parts 1,3,5,... correctly summed
         ! *** could nicely be generalised ****
         IF (MODULO(taskid, 2) == 0) THEN
            partner = taskid + 1
            DO i = 1, Nprocs, 2 ! sum the full grid with the partner
               CALL comm%sendrecv(grid3(:, i + 1), partner, lgrid3, partner, 17)
               grid3(:, i) = grid3(:, i) + lgrid3(:)
            END DO
         ELSE
            partner = taskid - 1
            DO i = 1, Nprocs, 2
               CALL comm%sendrecv(grid3(:, i), partner, lgrid3, partner, 17)
               grid3(:, i + 1) = grid3(:, i + 1) + lgrid3(:)
            END DO
         END IF
         t4 = MPI_WTIME() - t3 + threshold
         ! now send a given buffer from 1 to 3 to 5 .. adding the right part of the data
         ! since we've summed an lgrid does only need to pass by even or odd tasks
         left = MODULO(taskid - 2, Nprocs)
         right = MODULO(taskid + 2, Nprocs)
         t3 = MPI_WTIME()
         lgrid3 = 0.0E0_dp
         DO i = 1, Nprocs, 2
            lgrid3(:) = lgrid3 + grid3(:, MODULO(taskid - i - 1, Nprocs) + 1)
            IF (i .EQ. nprocs - 1) EXIT
            CALL comm%shift(lgrid3, 2)
         END DO
         t5 = MPI_WTIME() - t3 + threshold
         CALL comm%max(t4)
         CALL comm%max(t5)
         maxdiff = MAXVAL(ABS(lgrid3 - lgrid))
         CALL comm%max(maxdiff)
         IF (ionode .AND. output_unit > 0) WRITE (output_unit, *) "INVOLVED SHIFT        ", t4 + t5, "(", t4, ",", t5, ")", maxdiff
      END DO
      DEALLOCATE (rcount)
      DEALLOCATE (scount)
      DEALLOCATE (sdispl)
      DEALLOCATE (rdispl)
      DEALLOCATE (buffer1)
      DEALLOCATE (buffer2)
      DEALLOCATE (buffer3)
      DEALLOCATE (grid)
      DEALLOCATE (grid2)
      DEALLOCATE (grid3)
      DEALLOCATE (lgrid)
      DEALLOCATE (lgrid2)
      DEALLOCATE (lgrid3)
      DEALLOCATE (send_timings)
      DEALLOCATE (send_timings2)
#else
      MARK_USED(comm)
      MARK_USED(npow)
      IF (output_unit > 0) WRITE (output_unit, *) "No MPI tests for a serial program"
#endif
   END SUBROUTINE mpi_perf_test

END MODULE mp_perf_test
