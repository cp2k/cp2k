!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Defines all routines to deal with the performance of MPI routines
! **************************************************************************************************
MODULE mp_perf_env
   ! performance gathering
   USE kinds,                           ONLY: dp
#include "../base/base_uses.f90"

   PRIVATE

   PUBLIC :: mp_perf_env_type
   PUBLIC :: mp_perf_env_retain, mp_perf_env_release
   PUBLIC :: add_mp_perf_env, rm_mp_perf_env, get_mp_perf_env, describe_mp_perf_env
   PUBLIC :: add_perf

   TYPE mp_perf_type
      CHARACTER(LEN=20) :: name = ""
      INTEGER :: count = 0
      REAL(KIND=dp) :: msg_size = 0.0_dp
   END TYPE mp_perf_type

   INTEGER, PARAMETER :: MAX_PERF = 28

! **************************************************************************************************
   TYPE mp_perf_env_type
      PRIVATE
      INTEGER :: ref_count = -1
      TYPE(mp_perf_type), DIMENSION(MAX_PERF) :: mp_perfs = mp_perf_type()
   CONTAINS
      PROCEDURE, PUBLIC, PASS(perf_env), NON_OVERRIDABLE :: retain => mp_perf_env_retain
   END TYPE mp_perf_env_type

! **************************************************************************************************
   TYPE mp_perf_env_p_type
      TYPE(mp_perf_env_type), POINTER         :: mp_perf_env => Null()
   END TYPE mp_perf_env_p_type

   ! introduce a stack of mp_perfs, first index is the stack pointer, for convenience is replacing
   INTEGER, PARAMETER :: max_stack_size = 10
   INTEGER            :: stack_pointer = 0
   TYPE(mp_perf_env_p_type), DIMENSION(max_stack_size), SAVE :: mp_perf_stack

   CHARACTER(LEN=20), PARAMETER :: sname(MAX_PERF) = &
                                   (/"MP_Group            ", "MP_Bcast            ", "MP_Allreduce        ", &
                                     "MP_Gather           ", "MP_Sync             ", "MP_Alltoall         ", &
                                     "MP_SendRecv         ", "MP_ISendRecv        ", "MP_Wait             ", &
                                     "MP_comm_split       ", "MP_ISend            ", "MP_IRecv            ", &
                                     "MP_Send             ", "MP_Recv             ", "MP_Memory           ", &
                                     "MP_Put              ", "MP_Get              ", "MP_Fence            ", &
                                     "MP_Win_Lock         ", "MP_Win_Create       ", "MP_Win_Free         ", &
                                     "MP_IBcast           ", "MP_IAllreduce       ", "MP_IScatter         ", &
                                     "MP_RGet             ", "MP_Isync            ", "MP_Read_All         ", &
                                     "MP_Write_All        "/)

CONTAINS

! **************************************************************************************************
!> \brief start and stop the performance indicators
!>      for every call to start there has to be (exactly) one call to stop
!> \param perf_env ...
!> \par History
!>      2.2004 created [Joost VandeVondele]
!> \note
!>      can be used to measure performance of a sub-part of a program.
!>      timings measured here will not show up in the outer start/stops
!>      Doesn't need a fresh communicator
! **************************************************************************************************
   SUBROUTINE add_mp_perf_env(perf_env)
      TYPE(mp_perf_env_type), OPTIONAL, POINTER          :: perf_env

      stack_pointer = stack_pointer + 1
      IF (stack_pointer > max_stack_size) THEN
         CPABORT("stack_pointer too large : message_passing @ add_mp_perf_env")
      END IF
      NULLIFY (mp_perf_stack(stack_pointer)%mp_perf_env)
      IF (PRESENT(perf_env)) THEN
         mp_perf_stack(stack_pointer)%mp_perf_env => perf_env
         IF (ASSOCIATED(perf_env)) CALL mp_perf_env_retain(perf_env)
      END IF
      IF (.NOT. ASSOCIATED(mp_perf_stack(stack_pointer)%mp_perf_env)) THEN
         CALL mp_perf_env_create(mp_perf_stack(stack_pointer)%mp_perf_env)
      END IF
   END SUBROUTINE add_mp_perf_env

! **************************************************************************************************
!> \brief ...
!> \param perf_env ...
! **************************************************************************************************
   SUBROUTINE mp_perf_env_create(perf_env)
      TYPE(mp_perf_env_type), OPTIONAL, POINTER          :: perf_env

      INTEGER                                            :: i

      NULLIFY (perf_env)
      ALLOCATE (perf_env)
      perf_env%ref_count = 1
      DO i = 1, MAX_PERF
         perf_env%mp_perfs(i)%name = sname(i)
      END DO

   END SUBROUTINE mp_perf_env_create

! **************************************************************************************************
!> \brief ...
!> \param perf_env ...
! **************************************************************************************************
   SUBROUTINE mp_perf_env_release(perf_env)
      TYPE(mp_perf_env_type), POINTER                    :: perf_env

      IF (ASSOCIATED(perf_env)) THEN
         IF (perf_env%ref_count < 1) THEN
            CPABORT("invalid ref_count: message_passing @ mp_perf_env_release")
         END IF
         perf_env%ref_count = perf_env%ref_count - 1
         IF (perf_env%ref_count == 0) THEN
            DEALLOCATE (perf_env)
         END IF
      END IF
      NULLIFY (perf_env)
   END SUBROUTINE mp_perf_env_release

! **************************************************************************************************
!> \brief ...
!> \param perf_env ...
! **************************************************************************************************
   ELEMENTAL SUBROUTINE mp_perf_env_retain(perf_env)
      CLASS(mp_perf_env_type), INTENT(INOUT)                    :: perf_env

      perf_env%ref_count = perf_env%ref_count + 1
   END SUBROUTINE mp_perf_env_retain

!.. reports the performance counters for the MPI run
! **************************************************************************************************
!> \brief ...
!> \param perf_env ...
!> \param iw ...
! **************************************************************************************************
   SUBROUTINE mp_perf_env_describe(perf_env, iw)
      TYPE(mp_perf_env_type), INTENT(IN)       :: perf_env
      INTEGER, INTENT(IN)                      :: iw

#if defined(__parallel)
      INTEGER                                  :: i
      REAL(KIND=dp)                            :: vol
#endif

      IF (perf_env%ref_count < 1) THEN
         CPABORT("invalid perf_env%ref_count : message_passing @ mp_perf_env_describe")
      END IF
#if defined(__parallel)
      IF (iw > 0) THEN
         WRITE (iw, '( /, 1X, 79("-") )')
         WRITE (iw, '( " -", 77X, "-" )')
         WRITE (iw, '( " -", 24X, A, 24X, "-" )') ' MESSAGE PASSING PERFORMANCE '
         WRITE (iw, '( " -", 77X, "-" )')
         WRITE (iw, '( 1X, 79("-"), / )')
         WRITE (iw, '( A, A, A )') ' ROUTINE', '             CALLS ', &
            '     AVE VOLUME [Bytes]'
         DO i = 1, MAX_PERF

            IF (perf_env%mp_perfs(i)%count > 0) THEN
               vol = perf_env%mp_perfs(i)%msg_size/REAL(perf_env%mp_perfs(i)%count, KIND=dp)
               IF (vol < 1.0_dp) THEN
                  WRITE (iw, '(1X,A15,T17,I10)') &
                     ADJUSTL(perf_env%mp_perfs(i)%name), perf_env%mp_perfs(i)%count
               ELSE
                  WRITE (iw, '(1X,A15,T17,I10,T40,F11.0)') &
                     ADJUSTL(perf_env%mp_perfs(i)%name), perf_env%mp_perfs(i)%count, &
                     vol
               END IF
            END IF

         END DO
         WRITE (iw, '( 1X, 79("-"), / )')
      END IF
#else
      MARK_USED(iw)
#endif
   END SUBROUTINE mp_perf_env_describe

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE rm_mp_perf_env()
      IF (stack_pointer < 1) THEN
         CPABORT("no perf_env in the stack : message_passing @ rm_mp_perf_env")
      END IF
      CALL mp_perf_env_release(mp_perf_stack(stack_pointer)%mp_perf_env)
      stack_pointer = stack_pointer - 1
   END SUBROUTINE rm_mp_perf_env

! **************************************************************************************************
!> \brief ...
!> \return ...
! **************************************************************************************************
   FUNCTION get_mp_perf_env() RESULT(res)
      TYPE(mp_perf_env_type), POINTER                    :: res

      IF (stack_pointer < 1) THEN
         CPABORT("no perf_env in the stack : message_passing @ get_mp_perf_env")
      END IF
      res => mp_perf_stack(stack_pointer)%mp_perf_env
   END FUNCTION get_mp_perf_env

! **************************************************************************************************
!> \brief ...
!> \param scr ...
! **************************************************************************************************
   SUBROUTINE describe_mp_perf_env(scr)
      INTEGER, INTENT(in)                                :: scr

      TYPE(mp_perf_env_type), POINTER                    :: perf_env

      perf_env => get_mp_perf_env()
      CALL mp_perf_env_describe(perf_env, scr)
   END SUBROUTINE describe_mp_perf_env

! **************************************************************************************************
!> \brief adds the performance informations of one call
!> \param perf_id ...
!> \param count ...
!> \param msg_size ...
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE add_perf(perf_id, count, msg_size)
      INTEGER, INTENT(in)                      :: perf_id
      INTEGER, INTENT(in), OPTIONAL            :: count
      INTEGER, INTENT(in), OPTIONAL            :: msg_size

#if defined(__parallel)
      TYPE(mp_perf_type), POINTER              :: mp_perf

      IF (.NOT. ASSOCIATED(mp_perf_stack(stack_pointer)%mp_perf_env)) RETURN

      mp_perf => mp_perf_stack(stack_pointer)%mp_perf_env%mp_perfs(perf_id)
      IF (PRESENT(count)) THEN
         mp_perf%count = mp_perf%count + count
      END IF
      IF (PRESENT(msg_size)) THEN
         mp_perf%msg_size = mp_perf%msg_size + REAL(msg_size, dp)
      END IF
#else
      MARK_USED(perf_id)
      MARK_USED(count)
      MARK_USED(msg_size)
#endif

   END SUBROUTINE add_perf

END MODULE mp_perf_env
