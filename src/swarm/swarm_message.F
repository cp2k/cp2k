!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Swarm-message, a convenient data-container for with build-in serialization.
!> \author Ole Schuett
! **************************************************************************************************
MODULE swarm_message

   USE cp_parser_methods, ONLY: parser_get_next_line
   USE cp_parser_types, ONLY: cp_parser_type
   USE kinds, ONLY: default_string_length, &
                    int_4, &
                    int_8, &
                    real_4, &
                    real_8
   USE message_passing, ONLY: mp_comm_type
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'swarm_message'

   TYPE swarm_message_type
      PRIVATE
      TYPE(message_entry_type), POINTER :: root => Null()
   END TYPE swarm_message_type

   INTEGER, PARAMETER  :: key_length = 20

   TYPE message_entry_type
      CHARACTER(LEN=key_length)                      :: key = ""
      TYPE(message_entry_type), POINTER              :: next => Null()
      CHARACTER(LEN=default_string_length), POINTER  :: value_str => Null()
      INTEGER(KIND=int_4), POINTER                   :: value_i4 => Null()
      INTEGER(KIND=int_8), POINTER                   :: value_i8 => Null()
      REAL(KIND=real_4), POINTER                     :: value_r4 => Null()
      REAL(KIND=real_8), POINTER                     :: value_r8 => Null()
      INTEGER(KIND=int_4), DIMENSION(:), POINTER     :: value_1d_i4 => Null()
      INTEGER(KIND=int_8), DIMENSION(:), POINTER     :: value_1d_i8 => Null()
      REAL(KIND=real_4), DIMENSION(:), POINTER       :: value_1d_r4 => Null()
      REAL(KIND=real_8), DIMENSION(:), POINTER       :: value_1d_r8 => Null()
   END TYPE message_entry_type

! **************************************************************************************************
!> \brief Adds an entry from a swarm-message.
!> \author Ole Schuett
! **************************************************************************************************
   INTERFACE swarm_message_add
      MODULE PROCEDURE swarm_message_add_str
      MODULE PROCEDURE swarm_message_add_i4, swarm_message_add_i8
      MODULE PROCEDURE swarm_message_add_r4, swarm_message_add_r8
      MODULE PROCEDURE swarm_message_add_1d_i4, swarm_message_add_1d_i8
      MODULE PROCEDURE swarm_message_add_1d_r4, swarm_message_add_1d_r8
   END INTERFACE swarm_message_add

! **************************************************************************************************
!> \brief Returns an entry from a swarm-message.
!> \author Ole Schuett
! **************************************************************************************************
   INTERFACE swarm_message_get
      MODULE PROCEDURE swarm_message_get_str
      MODULE PROCEDURE swarm_message_get_i4, swarm_message_get_i8
      MODULE PROCEDURE swarm_message_get_r4, swarm_message_get_r8
      MODULE PROCEDURE swarm_message_get_1d_i4, swarm_message_get_1d_i8
      MODULE PROCEDURE swarm_message_get_1d_r4, swarm_message_get_1d_r8
   END INTERFACE swarm_message_get

   PUBLIC :: swarm_message_type, swarm_message_add, swarm_message_get
   PUBLIC :: swarm_message_mpi_send, swarm_message_mpi_recv, swarm_message_mpi_bcast
   PUBLIC :: swarm_message_file_write, swarm_message_file_read
   PUBLIC :: swarm_message_haskey, swarm_message_equal
   PUBLIC :: swarm_message_free

CONTAINS

! **************************************************************************************************
!> \brief Returns the number of entries contained in a swarm-message.
!> \param msg ...
!> \return ...
!> \author Ole Schuett
! **************************************************************************************************
   FUNCTION swarm_message_length(msg) RESULT(l)
      TYPE(swarm_message_type), INTENT(IN)               :: msg
      INTEGER                                            :: l

      TYPE(message_entry_type), POINTER                  :: curr_entry

      l = 0
      curr_entry => msg%root
      DO WHILE (ASSOCIATED(curr_entry))
         l = l + 1
         curr_entry => curr_entry%next
      END DO
   END FUNCTION swarm_message_length

! **************************************************************************************************
!> \brief Checks if a swarm-message contains an entry with the given key.
!> \param msg ...
!> \param key ...
!> \return ...
!> \author Ole Schuett
! **************************************************************************************************
   FUNCTION swarm_message_haskey(msg, key) RESULT(res)
      TYPE(swarm_message_type), INTENT(IN)               :: msg
      CHARACTER(LEN=*), INTENT(IN)                       :: key
      LOGICAL                                            :: res

      TYPE(message_entry_type), POINTER                  :: curr_entry

      res = .FALSE.
      curr_entry => msg%root
      DO WHILE (ASSOCIATED(curr_entry))
         IF (TRIM(curr_entry%key) == TRIM(key)) THEN
            res = .TRUE.
            EXIT
         END IF
         curr_entry => curr_entry%next
      END DO
   END FUNCTION swarm_message_haskey

! **************************************************************************************************
!> \brief Deallocates all entries contained in a swarm-message.
!> \param msg ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE swarm_message_free(msg)
      TYPE(swarm_message_type), INTENT(INOUT)            :: msg

      TYPE(message_entry_type), POINTER                  :: ENTRY, old_entry

      ENTRY => msg%root
      DO WHILE (ASSOCIATED(ENTRY))
         IF (ASSOCIATED(entry%value_str)) DEALLOCATE (entry%value_str)
         IF (ASSOCIATED(entry%value_i4)) DEALLOCATE (entry%value_i4)
         IF (ASSOCIATED(entry%value_i8)) DEALLOCATE (entry%value_i8)
         IF (ASSOCIATED(entry%value_r4)) DEALLOCATE (entry%value_r4)
         IF (ASSOCIATED(entry%value_r8)) DEALLOCATE (entry%value_r8)
         IF (ASSOCIATED(entry%value_1d_i4)) DEALLOCATE (entry%value_1d_i4)
         IF (ASSOCIATED(entry%value_1d_i8)) DEALLOCATE (entry%value_1d_i8)
         IF (ASSOCIATED(entry%value_1d_r4)) DEALLOCATE (entry%value_1d_r4)
         IF (ASSOCIATED(entry%value_1d_r8)) DEALLOCATE (entry%value_1d_r8)
         old_entry => ENTRY
         ENTRY => entry%next
         DEALLOCATE (old_entry)
      END DO

      NULLIFY (msg%root)

      CPASSERT(swarm_message_length(msg) == 0)
   END SUBROUTINE swarm_message_free

! **************************************************************************************************
!> \brief Checks if two swarm-messages are equal
!> \param msg1 ...
!> \param msg2 ...
!> \return ...
!> \author Ole Schuett
! **************************************************************************************************
   FUNCTION swarm_message_equal(msg1, msg2) RESULT(res)
      TYPE(swarm_message_type), INTENT(IN)               :: msg1, msg2
      LOGICAL                                            :: res

      res = swarm_message_equal_oneway(msg1, msg2) .AND. &
            swarm_message_equal_oneway(msg2, msg1)

   END FUNCTION swarm_message_equal

! **************************************************************************************************
!> \brief Sends a swarm message via MPI.
!> \param msg ...
!> \param group ...
!> \param dest ...
!> \param tag ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE swarm_message_mpi_send(msg, group, dest, tag)
      TYPE(swarm_message_type), INTENT(IN)               :: msg
      CLASS(mp_comm_type), INTENT(IN) :: group
      INTEGER, INTENT(IN)                                :: dest, tag

      TYPE(message_entry_type), POINTER                  :: curr_entry

      CALL group%send(swarm_message_length(msg), dest, tag)
      curr_entry => msg%root
      DO WHILE (ASSOCIATED(curr_entry))
         CALL swarm_message_entry_mpi_send(curr_entry, group, dest, tag)
         curr_entry => curr_entry%next
      END DO
   END SUBROUTINE swarm_message_mpi_send

! **************************************************************************************************
!> \brief Receives a swarm message via MPI.
!> \param msg ...
!> \param group ...
!> \param src ...
!> \param tag ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE swarm_message_mpi_recv(msg, group, src, tag)
      TYPE(swarm_message_type), INTENT(INOUT)            :: msg
      CLASS(mp_comm_type), INTENT(IN)                                :: group
      INTEGER, INTENT(INOUT)                             :: src, tag

      INTEGER                                            :: i, length
      TYPE(message_entry_type), POINTER                  :: new_entry

      IF (ASSOCIATED(msg%root)) CPABORT("message not empty")
      CALL group%recv(length, src, tag)
      DO i = 1, length
         ALLOCATE (new_entry)
         CALL swarm_message_entry_mpi_recv(new_entry, group, src, tag)
         new_entry%next => msg%root
         msg%root => new_entry
      END DO

   END SUBROUTINE swarm_message_mpi_recv

! **************************************************************************************************
!> \brief Broadcasts a swarm message via MPI.
!> \param msg ...
!> \param src ...
!> \param group ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE swarm_message_mpi_bcast(msg, src, group)
      TYPE(swarm_message_type), INTENT(INOUT)            :: msg
      INTEGER, INTENT(IN)                                :: src
      CLASS(mp_comm_type), INTENT(IN) :: group

      INTEGER                                            :: i, length
      TYPE(message_entry_type), POINTER                  :: curr_entry

      ASSOCIATE (mepos => group%mepos)

         IF (mepos /= src .AND. ASSOCIATED(msg%root)) CPABORT("message not empty")
         length = swarm_message_length(msg)
         CALL group%bcast(length, src)

         IF (mepos == src) curr_entry => msg%root

         DO i = 1, length
            IF (mepos /= src) ALLOCATE (curr_entry)

            CALL swarm_message_entry_mpi_bcast(curr_entry, src, group, mepos)

            IF (mepos == src) THEN
               curr_entry => curr_entry%next
            ELSE
               curr_entry%next => msg%root
               msg%root => curr_entry
            END IF
         END DO
      END ASSOCIATE

   END SUBROUTINE swarm_message_mpi_bcast

! **************************************************************************************************
!> \brief Write a swarm-message to a given file / unit.
!> \param msg ...
!> \param unit ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE swarm_message_file_write(msg, unit)
      TYPE(swarm_message_type), INTENT(IN)               :: msg
      INTEGER, INTENT(IN)                                :: unit

      INTEGER                                            :: handle
      TYPE(message_entry_type), POINTER                  :: curr_entry

      IF (unit <= 0) RETURN

      CALL timeset("swarm_message_file_write", handle)
      WRITE (unit, "(A)") "BEGIN SWARM_MESSAGE"
      WRITE (unit, "(A,I10)") "msg_length: ", swarm_message_length(msg)

      curr_entry => msg%root
      DO WHILE (ASSOCIATED(curr_entry))
         CALL swarm_message_entry_file_write(curr_entry, unit)
         curr_entry => curr_entry%next
      END DO

      WRITE (unit, "(A)") "END SWARM_MESSAGE"
      WRITE (unit, "()")
      CALL timestop(handle)
   END SUBROUTINE swarm_message_file_write

! **************************************************************************************************
!> \brief Reads a swarm-message from a given file / unit.
!> \param msg ...
!> \param parser ...
!> \param at_end ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE swarm_message_file_read(msg, parser, at_end)
      TYPE(swarm_message_type), INTENT(OUT)              :: msg
      TYPE(cp_parser_type), INTENT(INOUT)                :: parser
      LOGICAL, INTENT(INOUT)                             :: at_end

      INTEGER                                            :: handle

      CALL timeset("swarm_message_file_read", handle)
      CALL swarm_message_file_read_low(msg, parser, at_end)
      CALL timestop(handle)
   END SUBROUTINE swarm_message_file_read

! **************************************************************************************************
!> \brief Helper routine, does the actual work of swarm_message_file_read().
!> \param msg ...
!> \param parser ...
!> \param at_end ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE swarm_message_file_read_low(msg, parser, at_end)
      TYPE(swarm_message_type), INTENT(OUT)              :: msg
      TYPE(cp_parser_type), INTENT(INOUT)                :: parser
      LOGICAL, INTENT(INOUT)                             :: at_end

      CHARACTER(LEN=20)                                  :: label
      INTEGER                                            :: i, length
      TYPE(message_entry_type), POINTER                  :: new_entry

      CALL parser_get_next_line(parser, 1, at_end)
      at_end = at_end .OR. LEN_TRIM(parser%input_line(1:10)) == 0
      IF (at_end) RETURN
      CPASSERT(TRIM(parser%input_line(1:20)) == "BEGIN SWARM_MESSAGE")

      CALL parser_get_next_line(parser, 1, at_end)
      at_end = at_end .OR. LEN_TRIM(parser%input_line(1:10)) == 0
      IF (at_end) RETURN
      READ (parser%input_line(1:40), *) label, length
      CPASSERT(TRIM(label) == "msg_length:")

      DO i = 1, length
         ALLOCATE (new_entry)
         CALL swarm_message_entry_file_read(new_entry, parser, at_end)
         new_entry%next => msg%root
         msg%root => new_entry
      END DO

      CALL parser_get_next_line(parser, 1, at_end)
      at_end = at_end .OR. LEN_TRIM(parser%input_line(1:10)) == 0
      IF (at_end) RETURN
      CPASSERT(TRIM(parser%input_line(1:20)) == "END SWARM_MESSAGE")

   END SUBROUTINE swarm_message_file_read_low

! **************************************************************************************************
!> \brief Helper routine for swarm_message_equal
!> \param msg1 ...
!> \param msg2 ...
!> \return ...
!> \author Ole Schuett
! **************************************************************************************************
   FUNCTION swarm_message_equal_oneway(msg1, msg2) RESULT(res)
      TYPE(swarm_message_type), INTENT(IN)               :: msg1, msg2
      LOGICAL                                            :: res

      LOGICAL                                            :: found
      TYPE(message_entry_type), POINTER                  :: entry1, entry2

      res = .FALSE.

      !loop over entries of msg1
      entry1 => msg1%root
      DO WHILE (ASSOCIATED(entry1))

         ! finding matching entry in msg2
         entry2 => msg2%root
         found = .FALSE.
         DO WHILE (ASSOCIATED(entry2))
            IF (TRIM(entry2%key) == TRIM(entry1%key)) THEN
               found = .TRUE.
               EXIT
            END IF
            entry2 => entry2%next
         END DO
         IF (.NOT. found) RETURN

         !compare the two entries
         IF (ASSOCIATED(entry1%value_str)) THEN
            IF (.NOT. ASSOCIATED(entry2%value_str)) RETURN
            IF (TRIM(entry1%value_str) /= TRIM(entry2%value_str)) RETURN

         ELSE IF (ASSOCIATED(entry1%value_i4)) THEN
            IF (.NOT. ASSOCIATED(entry2%value_i4)) RETURN
            IF (entry1%value_i4 /= entry2%value_i4) RETURN

         ELSE IF (ASSOCIATED(entry1%value_i8)) THEN
            IF (.NOT. ASSOCIATED(entry2%value_i8)) RETURN
            IF (entry1%value_i8 /= entry2%value_i8) RETURN

         ELSE IF (ASSOCIATED(entry1%value_r4)) THEN
            IF (.NOT. ASSOCIATED(entry2%value_r4)) RETURN
            IF (ABS(entry1%value_r4 - entry2%value_r4) > 1e-5) RETURN

         ELSE IF (ASSOCIATED(entry1%value_r8)) THEN
            IF (.NOT. ASSOCIATED(entry2%value_r8)) RETURN
            IF (ABS(entry1%value_r8 - entry2%value_r8) > 1e-10) RETURN

         ELSE IF (ASSOCIATED(entry1%value_1d_i4)) THEN
            IF (.NOT. ASSOCIATED(entry2%value_1d_i4)) RETURN
            IF (ANY(entry1%value_1d_i4 /= entry2%value_1d_i4)) RETURN

         ELSE IF (ASSOCIATED(entry1%value_1d_i8)) THEN
            IF (.NOT. ASSOCIATED(entry2%value_1d_i8)) RETURN
            IF (ANY(entry1%value_1d_i8 /= entry2%value_1d_i8)) RETURN

         ELSE IF (ASSOCIATED(entry1%value_1d_r4)) THEN
            IF (.NOT. ASSOCIATED(entry2%value_1d_r4)) RETURN
            IF (ANY(ABS(entry1%value_1d_r4 - entry2%value_1d_r4) > 1e-5)) RETURN

         ELSE IF (ASSOCIATED(entry1%value_1d_r8)) THEN
            IF (.NOT. ASSOCIATED(entry2%value_1d_r8)) RETURN
            IF (ANY(ABS(entry1%value_1d_r8 - entry2%value_1d_r8) > 1e-10)) RETURN
         ELSE
            CPABORT("no value ASSOCIATED")
         END IF

         entry1 => entry1%next
      END DO

      ! if we reach this point no differences were found
      res = .TRUE.
   END FUNCTION swarm_message_equal_oneway

! **************************************************************************************************
!> \brief Helper routine for swarm_message_mpi_send.
!> \param ENTRY ...
!> \param group ...
!> \param dest ...
!> \param tag ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE swarm_message_entry_mpi_send(ENTRY, group, dest, tag)
      TYPE(message_entry_type), INTENT(IN)               :: ENTRY
      CLASS(mp_comm_type), INTENT(IN) :: group
      INTEGER, INTENT(IN)                                :: dest, tag

      INTEGER, DIMENSION(default_string_length)          :: value_str_arr
      INTEGER, DIMENSION(key_length)                     :: key_arr

      key_arr = str2iarr(entry%key)
      CALL group%send(key_arr, dest, tag)

      IF (ASSOCIATED(entry%value_i4)) THEN
         CALL group%send(1, dest, tag)
         CALL group%send(entry%value_i4, dest, tag)

      ELSE IF (ASSOCIATED(entry%value_i8)) THEN
         CALL group%send(2, dest, tag)
         CALL group%send(entry%value_i8, dest, tag)

      ELSE IF (ASSOCIATED(entry%value_r4)) THEN
         CALL group%send(3, dest, tag)
         CALL group%send(entry%value_r4, dest, tag)

      ELSE IF (ASSOCIATED(entry%value_r8)) THEN
         CALL group%send(4, dest, tag)
         CALL group%send(entry%value_r8, dest, tag)

      ELSE IF (ASSOCIATED(entry%value_1d_i4)) THEN
         CALL group%send(5, dest, tag)
         CALL group%send(SIZE(entry%value_1d_i4), dest, tag)
         CALL group%send(entry%value_1d_i4, dest, tag)

      ELSE IF (ASSOCIATED(entry%value_1d_i8)) THEN
         CALL group%send(6, dest, tag)
         CALL group%send(SIZE(entry%value_1d_i8), dest, tag)
         CALL group%send(entry%value_1d_i8, dest, tag)

      ELSE IF (ASSOCIATED(entry%value_1d_r4)) THEN
         CALL group%send(7, dest, tag)
         CALL group%send(SIZE(entry%value_1d_r4), dest, tag)
         CALL group%send(entry%value_1d_r4, dest, tag)

      ELSE IF (ASSOCIATED(entry%value_1d_r8)) THEN
         CALL group%send(8, dest, tag)
         CALL group%send(SIZE(entry%value_1d_r8), dest, tag)
         CALL group%send(entry%value_1d_r8, dest, tag)

      ELSE IF (ASSOCIATED(entry%value_str)) THEN
         CALL group%send(9, dest, tag)
         value_str_arr = str2iarr(entry%value_str)
         CALL group%send(value_str_arr, dest, tag)
      ELSE
         CPABORT("no value ASSOCIATED")
      END IF
   END SUBROUTINE swarm_message_entry_mpi_send

! **************************************************************************************************
!> \brief Helper routine for swarm_message_mpi_recv.
!> \param ENTRY ...
!> \param group ...
!> \param src ...
!> \param tag ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE swarm_message_entry_mpi_recv(ENTRY, group, src, tag)
      TYPE(message_entry_type), INTENT(INOUT)            :: ENTRY
      CLASS(mp_comm_type), INTENT(IN)                                :: group
      INTEGER, INTENT(INOUT)                             :: src, tag

      INTEGER                                            :: datatype, s
      INTEGER, DIMENSION(default_string_length)          :: value_str_arr
      INTEGER, DIMENSION(key_length)                     :: key_arr

      CALL group%recv(key_arr, src, tag)
      entry%key = iarr2str(key_arr)

      CALL group%recv(datatype, src, tag)

      SELECT CASE (datatype)
      CASE (1)
         ALLOCATE (entry%value_i4)
         CALL group%recv(entry%value_i4, src, tag)
      CASE (2)
         ALLOCATE (entry%value_i8)
         CALL group%recv(entry%value_i8, src, tag)
      CASE (3)
         ALLOCATE (entry%value_r4)
         CALL group%recv(entry%value_r4, src, tag)
      CASE (4)
         ALLOCATE (entry%value_r8)
         CALL group%recv(entry%value_r8, src, tag)
      CASE (5)
         CALL group%recv(s, src, tag)
         ALLOCATE (entry%value_1d_i4(s))
         CALL group%recv(entry%value_1d_i4, src, tag)
      CASE (6)
         CALL group%recv(s, src, tag)
         ALLOCATE (entry%value_1d_i8(s))
         CALL group%recv(entry%value_1d_i8, src, tag)
      CASE (7)
         CALL group%recv(s, src, tag)
         ALLOCATE (entry%value_1d_r4(s))
         CALL group%recv(entry%value_1d_r4, src, tag)
      CASE (8)
         CALL group%recv(s, src, tag)
         ALLOCATE (entry%value_1d_r8(s))
         CALL group%recv(entry%value_1d_r8, src, tag)
      CASE (9)
         ALLOCATE (entry%value_str)
         CALL group%recv(value_str_arr, src, tag)
         entry%value_str = iarr2str(value_str_arr)
      CASE DEFAULT
         CPABORT("unknown datatype")
      END SELECT
   END SUBROUTINE swarm_message_entry_mpi_recv

! **************************************************************************************************
!> \brief Helper routine for swarm_message_mpi_bcast.
!> \param ENTRY ...
!> \param src ...
!> \param group ...
!> \param mepos ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE swarm_message_entry_mpi_bcast(ENTRY, src, group, mepos)
      TYPE(message_entry_type), INTENT(INOUT)            :: ENTRY
      INTEGER, INTENT(IN)                                :: src, mepos
      CLASS(mp_comm_type), INTENT(IN) :: group

      INTEGER                                            :: datasize, datatype
      INTEGER, DIMENSION(default_string_length)          :: value_str_arr
      INTEGER, DIMENSION(key_length)                     :: key_arr

      IF (src == mepos) key_arr = str2iarr(entry%key)
      CALL group%bcast(key_arr, src)
      IF (src /= mepos) entry%key = iarr2str(key_arr)

      IF (src == mepos) THEN
         datasize = 1
         IF (ASSOCIATED(entry%value_i4)) THEN
            datatype = 1
         ELSE IF (ASSOCIATED(entry%value_i8)) THEN
            datatype = 2
         ELSE IF (ASSOCIATED(entry%value_r4)) THEN
            datatype = 3
         ELSE IF (ASSOCIATED(entry%value_r8)) THEN
            datatype = 4
         ELSE IF (ASSOCIATED(entry%value_1d_i4)) THEN
            datatype = 5
            datasize = SIZE(entry%value_1d_i4)
         ELSE IF (ASSOCIATED(entry%value_1d_i8)) THEN
            datatype = 6
            datasize = SIZE(entry%value_1d_i8)
         ELSE IF (ASSOCIATED(entry%value_1d_r4)) THEN
            datatype = 7
            datasize = SIZE(entry%value_1d_r4)
         ELSE IF (ASSOCIATED(entry%value_1d_r8)) THEN
            datatype = 8
            datasize = SIZE(entry%value_1d_r8)
         ELSE IF (ASSOCIATED(entry%value_str)) THEN
            datatype = 9
         ELSE
            CPABORT("no value ASSOCIATED")
         END IF
      END IF
      CALL group%bcast(datatype, src)
      CALL group%bcast(datasize, src)

      SELECT CASE (datatype)
      CASE (1)
         IF (src /= mepos) ALLOCATE (entry%value_i4)
         CALL group%bcast(entry%value_i4, src)
      CASE (2)
         IF (src /= mepos) ALLOCATE (entry%value_i8)
         CALL group%bcast(entry%value_i8, src)
      CASE (3)
         IF (src /= mepos) ALLOCATE (entry%value_r4)
         CALL group%bcast(entry%value_r4, src)
      CASE (4)
         IF (src /= mepos) ALLOCATE (entry%value_r8)
         CALL group%bcast(entry%value_r8, src)
      CASE (5)
         IF (src /= mepos) ALLOCATE (entry%value_1d_i4(datasize))
         CALL group%bcast(entry%value_1d_i4, src)
      CASE (6)
         IF (src /= mepos) ALLOCATE (entry%value_1d_i8(datasize))
         CALL group%bcast(entry%value_1d_i8, src)
      CASE (7)
         IF (src /= mepos) ALLOCATE (entry%value_1d_r4(datasize))
         CALL group%bcast(entry%value_1d_r4, src)
      CASE (8)
         IF (src /= mepos) ALLOCATE (entry%value_1d_r8(datasize))
         CALL group%bcast(entry%value_1d_r8, src)
      CASE (9)
         IF (src == mepos) value_str_arr = str2iarr(entry%value_str)
         CALL group%bcast(value_str_arr, src)
         IF (src /= mepos) THEN
            ALLOCATE (entry%value_str)
            entry%value_str = iarr2str(value_str_arr)
         END IF
      CASE DEFAULT
         CPABORT("unknown datatype")
      END SELECT

   END SUBROUTINE swarm_message_entry_mpi_bcast

! **************************************************************************************************
!> \brief Helper routine for swarm_message_file_write.
!> \param ENTRY ...
!> \param unit ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE swarm_message_entry_file_write(ENTRY, unit)
      TYPE(message_entry_type), INTENT(IN)               :: ENTRY
      INTEGER, INTENT(IN)                                :: unit

      INTEGER                                            :: i

      WRITE (unit, "(A,A)") "key: ", entry%key
      IF (ASSOCIATED(entry%value_i4)) THEN
         WRITE (unit, "(A)") "datatype: i4"
         WRITE (unit, "(A,I10)") "value: ", entry%value_i4

      ELSE IF (ASSOCIATED(entry%value_i8)) THEN
         WRITE (unit, "(A)") "datatype: i8"
         WRITE (unit, "(A,I20)") "value: ", entry%value_i8

      ELSE IF (ASSOCIATED(entry%value_r4)) THEN
         WRITE (unit, "(A)") "datatype: r4"
         WRITE (unit, "(A,E30.20)") "value: ", entry%value_r4

      ELSE IF (ASSOCIATED(entry%value_r8)) THEN
         WRITE (unit, "(A)") "datatype: r8"
         WRITE (unit, "(A,E30.20)") "value: ", entry%value_r8

      ELSE IF (ASSOCIATED(entry%value_str)) THEN
         WRITE (unit, "(A)") "datatype: str"
         WRITE (unit, "(A,A)") "value: ", entry%value_str

      ELSE IF (ASSOCIATED(entry%value_1d_i4)) THEN
         WRITE (unit, "(A)") "datatype: 1d_i4"
         WRITE (unit, "(A,I10)") "size: ", SIZE(entry%value_1d_i4)
         DO i = 1, SIZE(entry%value_1d_i4)
            WRITE (unit, *) entry%value_1d_i4(i)
         END DO

      ELSE IF (ASSOCIATED(entry%value_1d_i8)) THEN
         WRITE (unit, "(A)") "datatype: 1d_i8"
         WRITE (unit, "(A,I20)") "size: ", SIZE(entry%value_1d_i8)
         DO i = 1, SIZE(entry%value_1d_i8)
            WRITE (unit, *) entry%value_1d_i8(i)
         END DO

      ELSE IF (ASSOCIATED(entry%value_1d_r4)) THEN
         WRITE (unit, "(A)") "datatype: 1d_r4"
         WRITE (unit, "(A,I8)") "size: ", SIZE(entry%value_1d_r4)
         DO i = 1, SIZE(entry%value_1d_r4)
            WRITE (unit, "(1X,E30.20)") entry%value_1d_r4(i)
         END DO

      ELSE IF (ASSOCIATED(entry%value_1d_r8)) THEN
         WRITE (unit, "(A)") "datatype: 1d_r8"
         WRITE (unit, "(A,I8)") "size: ", SIZE(entry%value_1d_r8)
         DO i = 1, SIZE(entry%value_1d_r8)
            WRITE (unit, "(1X,E30.20)") entry%value_1d_r8(i)
         END DO

      ELSE
         CPABORT("no value ASSOCIATED")
      END IF
   END SUBROUTINE swarm_message_entry_file_write

! **************************************************************************************************
!> \brief Helper routine for swarm_message_file_read.
!> \param ENTRY ...
!> \param parser ...
!> \param at_end ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE swarm_message_entry_file_read(ENTRY, parser, at_end)
      TYPE(message_entry_type), INTENT(INOUT)            :: ENTRY
      TYPE(cp_parser_type), INTENT(INOUT)                :: parser
      LOGICAL, INTENT(INOUT)                             :: at_end

      CHARACTER(LEN=15)                                  :: datatype, label
      INTEGER                                            :: arr_size, i
      LOGICAL                                            :: is_scalar

      CALL parser_get_next_line(parser, 1, at_end)
      at_end = at_end .OR. LEN_TRIM(parser%input_line(1:10)) == 0
      IF (at_end) RETURN
      READ (parser%input_line(1:key_length + 10), *) label, entry%key
      CPASSERT(TRIM(label) == "key:")

      CALL parser_get_next_line(parser, 1, at_end)
      at_end = at_end .OR. LEN_TRIM(parser%input_line(1:10)) == 0
      IF (at_end) RETURN
      READ (parser%input_line(1:30), *) label, datatype
      CPASSERT(TRIM(label) == "datatype:")

      CALL parser_get_next_line(parser, 1, at_end)
      at_end = at_end .OR. LEN_TRIM(parser%input_line(1:10)) == 0
      IF (at_end) RETURN

      is_scalar = .TRUE.
      SELECT CASE (TRIM(datatype))
      CASE ("i4")
         ALLOCATE (entry%value_i4)
         READ (parser%input_line(1:40), *) label, entry%value_i4
      CASE ("i8")
         ALLOCATE (entry%value_i8)
         READ (parser%input_line(1:40), *) label, entry%value_i8
      CASE ("r4")
         ALLOCATE (entry%value_r4)
         READ (parser%input_line(1:40), *) label, entry%value_r4
      CASE ("r8")
         ALLOCATE (entry%value_r8)
         READ (parser%input_line(1:40), *) label, entry%value_r8
      CASE ("str")
         ALLOCATE (entry%value_str)
         READ (parser%input_line(1:40), *) label, entry%value_str
      CASE DEFAULT
         is_scalar = .FALSE.
      END SELECT

      IF (is_scalar) THEN
         CPASSERT(TRIM(label) == "value:")
         RETURN
      END IF

      ! musst be an array-datatype
      READ (parser%input_line(1:30), *) label, arr_size
      CPASSERT(TRIM(label) == "size:")

      SELECT CASE (TRIM(datatype))
      CASE ("1d_i4")
         ALLOCATE (entry%value_1d_i4(arr_size))
      CASE ("1d_i8")
         ALLOCATE (entry%value_1d_i8(arr_size))
      CASE ("1d_r4")
         ALLOCATE (entry%value_1d_r4(arr_size))
      CASE ("1d_r8")
         ALLOCATE (entry%value_1d_r8(arr_size))
      CASE DEFAULT
         CPABORT("unknown datatype")
      END SELECT

      DO i = 1, arr_size
         CALL parser_get_next_line(parser, 1, at_end)
         at_end = at_end .OR. LEN_TRIM(parser%input_line(1:10)) == 0
         IF (at_end) RETURN

         !Numbers were written with at most 31 characters.
         SELECT CASE (TRIM(datatype))
         CASE ("1d_i4")
            READ (parser%input_line(1:31), *) entry%value_1d_i4(i)
         CASE ("1d_i8")
            READ (parser%input_line(1:31), *) entry%value_1d_i8(i)
         CASE ("1d_r4")
            READ (parser%input_line(1:31), *) entry%value_1d_r4(i)
         CASE ("1d_r8")
            READ (parser%input_line(1:31), *) entry%value_1d_r8(i)
         CASE DEFAULT
            CPABORT("swarm_message_entry_file_read: unknown datatype")
         END SELECT
      END DO

   END SUBROUTINE swarm_message_entry_file_read

! **************************************************************************************************
!> \brief Helper routine, converts a string into an integer-array
!> \param str ...
!> \return ...
!> \author Ole Schuett
! **************************************************************************************************
   PURE FUNCTION str2iarr(str) RESULT(arr)
      CHARACTER(LEN=*), INTENT(IN)                       :: str
      INTEGER, DIMENSION(LEN(str))                       :: arr

      INTEGER                                            :: i

      DO i = 1, LEN(str)
         arr(i) = ICHAR(str(i:i))
      END DO
   END FUNCTION str2iarr

! **************************************************************************************************
!> \brief Helper routine, converts an integer-array into a string
!> \param arr ...
!> \return ...
!> \author Ole Schuett
! **************************************************************************************************
   PURE FUNCTION iarr2str(arr) RESULT(str)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: arr
      CHARACTER(LEN=SIZE(arr))                           :: str

      INTEGER                                            :: i

      DO i = 1, SIZE(arr)
         str(i:i) = CHAR(arr(i))
      END DO
   END FUNCTION iarr2str

   #:set instances = {'str'   : 'CHARACTER(LEN=*)', &
      'i4'    : 'INTEGER(KIND=int_4)', &
      'i8'    : 'INTEGER(KIND=int_8)', &
      'r4'    : 'REAL(KIND=real_4)', &
      'r8'    : 'REAL(KIND=real_8)', &
      '1d_i4' : 'INTEGER(KIND=int_4), DIMENSION(:)', &
      '1d_i8' : 'INTEGER(KIND=int_8), DIMENSION(:)', &
      '1d_r4' : 'REAL(KIND=real_4), DIMENSION(:)', &
      '1d_r8' : 'REAL(KIND=real_8), DIMENSION(:)' }

   #:for label, type in instances.items()

! **************************************************************************************************
!> \brief Addes an entry from a swarm-message.
!> \param msg ...
!> \param key ...
!> \param value ...
!> \author Ole Schuett
! **************************************************************************************************
      SUBROUTINE swarm_message_add_${label}$ (msg, key, value)
         TYPE(swarm_message_type), INTENT(INOUT)   :: msg
         CHARACTER(LEN=*), INTENT(IN)              :: key
         ${type}$, INTENT(IN)                      :: value

         TYPE(message_entry_type), POINTER :: new_entry

         IF (swarm_message_haskey(msg, key)) &
            CPABORT("swarm_message_add_${label}$: key already exists: "//TRIM(key))

         ALLOCATE (new_entry)
         new_entry%key = key

         #:if label.startswith("1d_")
            ALLOCATE (new_entry%value_${label}$ (SIZE(value)))
         #:else
            ALLOCATE (new_entry%value_${label}$)
         #:endif

         new_entry%value_${label}$ = value

         !WRITE (*,*) "swarm_message_add_${label}$: key=",key, " value=",new_entry%value_${label}$

         IF (.NOT. ASSOCIATED(msg%root)) THEN
            msg%root => new_entry
         ELSE
            new_entry%next => msg%root
            msg%root => new_entry
         END IF

      END SUBROUTINE swarm_message_add_${label}$

! **************************************************************************************************
!> \brief Returns an entry from a swarm-message.
!> \param msg ...
!> \param key ...
!> \param value ...
!> \author Ole Schuett
! **************************************************************************************************
      SUBROUTINE swarm_message_get_${label}$ (msg, key, value)
         TYPE(swarm_message_type), INTENT(IN)  :: msg
         CHARACTER(LEN=*), INTENT(IN)          :: key

         #:if label=="str"
            CHARACTER(LEN=default_string_length)  :: value
         #:elif label.startswith("1d_")
            ${type}$, POINTER                     :: value
         #:else
            ${type}$, INTENT(OUT)                 :: value
         #:endif

         TYPE(message_entry_type), POINTER :: curr_entry
         !WRITE (*,*) "swarm_message_get_${label}$: key=",key

         #:if label.startswith("1d_")
            IF (ASSOCIATED(value)) CPABORT("swarm_message_get_${label}$: value already associated")
         #:endif

         curr_entry => msg%root
         DO WHILE (ASSOCIATED(curr_entry))
            IF (TRIM(curr_entry%key) == TRIM(key)) THEN
               IF (.NOT. ASSOCIATED(curr_entry%value_${label}$)) &
                  CPABORT("swarm_message_get_${label}$: value not associated key: "//TRIM(key))
               #:if label.startswith("1d_")
                  ALLOCATE (value(SIZE(curr_entry%value_${label}$)))
               #:endif
               value = curr_entry%value_${label}$
               !WRITE (*,*) "swarm_message_get_${label}$: value=",value
               RETURN
            END IF
            curr_entry => curr_entry%next
         END DO
         CPABORT("swarm_message_get: key not found: "//TRIM(key))
      END SUBROUTINE swarm_message_get_${label}$

   #:endfor

END MODULE swarm_message

