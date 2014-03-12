!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Swarm-message, a convenient data-container for with build-in serialization.
!> \author Ole Schuett
! *****************************************************************************
MODULE swarm_message

  USE cp_parser_methods,               ONLY: parser_get_next_line
  USE cp_parser_types,                 ONLY: cp_parser_type
  USE kinds,                           ONLY: default_string_length,&
                                             int_4,&
                                             int_8,&
                                             real_4,&
                                             real_8
  USE message_passing,                 ONLY: mp_abort,&
                                             mp_bcast,&
                                             mp_environ,&
                                             mp_recv,&
                                             mp_send
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

 IMPLICIT NONE
 PRIVATE

 TYPE swarm_message_type
   PRIVATE
   TYPE(message_entry_type), POINTER :: root => Null()
 END TYPE swarm_message_type

 INTEGER, PARAMETER  :: key_length      = 20

 TYPE message_entry_type
   CHARACTER(LEN=key_length)                        :: key
   TYPE(message_entry_type),            POINTER     :: next        => Null()
   CHARACTER(LEN=default_string_length),POINTER     :: value_str   => Null()
   INTEGER(KIND=int_4),                 POINTER     :: value_i4    => Null()
   INTEGER(KIND=int_8),                 POINTER     :: value_i8    => Null()
   REAL(KIND=real_4),                   POINTER     :: value_r4    => Null()
   REAL(KIND=real_8),                   POINTER     :: value_r8    => Null()
   INTEGER(KIND=int_4), DIMENSION(:),   POINTER     :: value_i4_1d => Null()
   INTEGER(KIND=int_8), DIMENSION(:),   POINTER     :: value_i8_1d => Null()
   REAL(KIND=real_4),   DIMENSION(:),   POINTER     :: value_r4_1d => Null()
   REAL(KIND=real_8),   DIMENSION(:),   POINTER     :: value_r8_1d => Null()
 END TYPE message_entry_type


! *****************************************************************************
!> \brief Adds an entry from a swarm-message.
!> \author Ole Schuett
! *****************************************************************************
 INTERFACE swarm_message_add
    MODULE PROCEDURE swarm_message_add_str
    MODULE PROCEDURE swarm_message_add_i4,    swarm_message_add_i8
    MODULE PROCEDURE swarm_message_add_r4,    swarm_message_add_r8
    MODULE PROCEDURE swarm_message_add_i4_1d, swarm_message_add_i8_1d
    MODULE PROCEDURE swarm_message_add_r4_1d, swarm_message_add_r8_1d
 END INTERFACE swarm_message_add


! *****************************************************************************
!> \brief Returns an entry from a swarm-message.
!> \author Ole Schuett
! *****************************************************************************
 INTERFACE swarm_message_get
    MODULE PROCEDURE swarm_message_get_str
    MODULE PROCEDURE swarm_message_get_i4,    swarm_message_get_i8
    MODULE PROCEDURE swarm_message_get_r4,    swarm_message_get_r8
    MODULE PROCEDURE swarm_message_get_i4_1d, swarm_message_get_i8_1d
    MODULE PROCEDURE swarm_message_get_r4_1d, swarm_message_get_r8_1d
 END INTERFACE swarm_message_get


 PUBLIC :: swarm_message_type, swarm_message_add, swarm_message_get
 PUBLIC :: swarm_message_mpi_send, swarm_message_mpi_recv, swarm_message_mpi_bcast
 PUBLIC :: swarm_message_file_write, swarm_message_file_read
 PUBLIC :: swarm_message_length, swarm_message_haskey, swarm_message_equal
 PUBLIC :: swarm_message_free


 CONTAINS


! *****************************************************************************
!> \brief Returns the number of entries contained in a swarm-message.
!> \author Ole Schuett
! *****************************************************************************
  FUNCTION swarm_message_length(msg) RESULT(l)
    TYPE(swarm_message_type), INTENT(IN)     :: msg
    INTEGER                                  :: l

    TYPE(message_entry_type), POINTER        :: curr_entry

    l = 0
    curr_entry => msg%root
    DO WHILE(ASSOCIATED(curr_entry))
      l = l + 1
      curr_entry => curr_entry%next
    END DO
  END FUNCTION swarm_message_length


! *****************************************************************************
!> \brief Checks if a swarm-message contains an entry with the given key.
!> \author Ole Schuett
! *****************************************************************************
  FUNCTION swarm_message_haskey(msg, key) RESULT(res)
    TYPE(swarm_message_type), INTENT(IN)     :: msg
    CHARACTER(LEN=*), INTENT(IN)             :: key
    LOGICAL                                  :: res

    TYPE(message_entry_type), POINTER        :: curr_entry

    res = .FALSE.
    curr_entry => msg%root
    DO WHILE(ASSOCIATED(curr_entry))
      IF(TRIM(curr_entry%key) == TRIM(key)) THEN
         res = .TRUE.
         EXIT
      END IF
      curr_entry => curr_entry%next
    END DO
  END FUNCTION swarm_message_haskey


! *****************************************************************************
!> \brief Deallocates all entries contained in a swarm-message.
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE swarm_message_free(msg)
    TYPE(swarm_message_type), INTENT(INOUT)  :: msg

    TYPE(message_entry_type), POINTER        :: ENTRY, old_entry

     ENTRY => msg%root
     DO WHILE(ASSOCIATED(ENTRY))
       IF(ASSOCIATED(entry%value_str))      DEALLOCATE(entry%value_str)
       IF(ASSOCIATED(entry%value_i4))       DEALLOCATE(entry%value_i4)
       IF(ASSOCIATED(entry%value_i8))       DEALLOCATE(entry%value_i8)
       IF(ASSOCIATED(entry%value_r4))       DEALLOCATE(entry%value_r4)
       IF(ASSOCIATED(entry%value_r8))       DEALLOCATE(entry%value_r8)
       IF(ASSOCIATED(entry%value_i4_1d))    DEALLOCATE(entry%value_i4_1d)
       IF(ASSOCIATED(entry%value_i8_1d))    DEALLOCATE(entry%value_i8_1d)
       IF(ASSOCIATED(entry%value_r4_1d))    DEALLOCATE(entry%value_r4_1d)
       IF(ASSOCIATED(entry%value_r8_1d))    DEALLOCATE(entry%value_r8_1d)
       old_entry => ENTRY
       ENTRY => entry%next
       DEALLOCATE(old_entry)
     END DO

     NULLIFY(msg%root)

     IF(swarm_message_length(msg) > 0) STOP "swarm_message_free failed"
  END SUBROUTINE swarm_message_free


! *****************************************************************************
!> \brief Checks if two swarm-messages are equal
!> \author Ole Schuett
! *****************************************************************************
  FUNCTION swarm_message_equal(msg1, msg2) RESULT(res)
    TYPE(swarm_message_type), INTENT(IN)     :: msg1, msg2
    LOGICAL                                  :: res

     res =  swarm_message_equal_oneway(msg1, msg2) .AND. &
            swarm_message_equal_oneway(msg2, msg1)

  END FUNCTION swarm_message_equal


! *****************************************************************************
!> \brief Sends a swarm message via MPI.
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE swarm_message_mpi_send(msg, group, dest, tag)
    TYPE(swarm_message_type), INTENT(IN)     :: msg
    INTEGER, INTENT(IN)                      :: group, dest, tag

    TYPE(message_entry_type), POINTER        :: curr_entry

    CALL mp_send(swarm_message_length(msg), dest, tag, group)
    curr_entry => msg%root
    DO WHILE(ASSOCIATED(curr_entry))
      CALL swarm_message_entry_mpi_send(curr_entry, group, dest, tag)
      curr_entry => curr_entry%next
    END DO
  END SUBROUTINE swarm_message_mpi_send


! *****************************************************************************
!> \brief Receives a swarm message via MPI.
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE swarm_message_mpi_recv(msg, group, src, tag)
    TYPE(swarm_message_type), INTENT(INOUT)  :: msg
    INTEGER, INTENT(IN)                      :: group
    INTEGER, INTENT(INOUT)                   :: src, tag

    INTEGER                                  :: i, length
    TYPE(message_entry_type), POINTER        :: new_entry

    IF(ASSOCIATED(msg%root)) STOP "swarm_message_mpi_recv: message not empty"
    CALL mp_recv(length, src, tag, group)
    DO i=1, length
       ALLOCATE(new_entry)
       CALL swarm_message_entry_mpi_recv(new_entry, group, src, tag)
       new_entry%next => msg%root
       msg%root => new_entry
    END DO

  END SUBROUTINE swarm_message_mpi_recv


! *****************************************************************************
!> \brief Broadcasts a swarm message via MPI.
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE swarm_message_mpi_bcast(msg, src, group)
    TYPE(swarm_message_type), INTENT(INOUT)  :: msg
    INTEGER, INTENT(IN)                      :: src, group

    INTEGER                                  :: i, length, mepos, num_pe
    TYPE(message_entry_type), POINTER        :: curr_entry

    CALL mp_environ(num_pe, mepos, group)

    IF(mepos/=src .AND. ASSOCIATED(msg%root)) STOP "swarm_message_mpi_bcast: message not empty"
    length = swarm_message_length(msg)
    CALL mp_bcast(length, src, group)

    IF(mepos==src) curr_entry => msg%root

    DO i=1, length
       IF(mepos/=src) ALLOCATE(curr_entry)

       CALL swarm_message_entry_mpi_bcast(curr_entry, src, group, mepos)

       IF(mepos==src) THEN
          curr_entry => curr_entry%next
       ELSE
          curr_entry%next => msg%root
          msg%root => curr_entry
       END IF
    END DO

  END SUBROUTINE swarm_message_mpi_bcast



! *****************************************************************************
!> \brief Write a swarm-message to a given file / unit.
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE swarm_message_file_write(msg, unit)
    TYPE(swarm_message_type), INTENT(IN)     :: msg
    INTEGER, INTENT(IN)                      :: unit

    INTEGER                                  :: handle
    TYPE(message_entry_type), POINTER        :: curr_entry

    IF(unit <= 0) RETURN

    CALL timeset("swarm_message_file_write", handle)
    WRITE(unit,"(A)") "BEGIN SWARM_MESSAGE"
    WRITE(unit,"(A,I10)") "msg_length: ", swarm_message_length(msg)

    curr_entry => msg%root
    DO WHILE(ASSOCIATED(curr_entry))
      CALL swarm_message_entry_file_write(curr_entry, unit)
      curr_entry => curr_entry%next
    END DO

    WRITE(unit,"(A)") "END SWARM_MESSAGE"
    WRITE(unit,"()")
    CALL timestop(handle)
  END SUBROUTINE swarm_message_file_write


! *****************************************************************************
!> \brief Reads a swarm-message from a given file / unit.
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE swarm_message_file_read(msg, parser, at_end, error)
    TYPE(swarm_message_type), INTENT(OUT)    :: msg
    TYPE(cp_parser_type), POINTER            :: parser
    LOGICAL, INTENT(INOUT)                   :: at_end
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    INTEGER                                  :: handle

    CALL timeset("swarm_message_file_read", handle)
    CALL swarm_message_file_read_low(msg, parser, at_end, error)
    CALL timestop(handle)
  END SUBROUTINE swarm_message_file_read


! *****************************************************************************
!> \brief Helper routine, does the actual work of swarm_message_file_read().
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE swarm_message_file_read_low(msg, parser, at_end, error)
    TYPE(swarm_message_type), INTENT(OUT)    :: msg
    TYPE(cp_parser_type), POINTER            :: parser
    LOGICAL, INTENT(INOUT)                   :: at_end
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=20)                        :: label
    INTEGER                                  :: i, length
    TYPE(message_entry_type), POINTER        :: new_entry

!IF(ASSOCIATED(msg%root)) STOP "swarm_message_file_read: message not empty"

    CALL parser_get_next_line(parser, 1, at_end, error)
    at_end = at_end .OR. LEN_TRIM(parser%input_line(1:10))==0
    IF(at_end) RETURN
    IF(TRIM(parser%input_line(1:20))/="BEGIN SWARM_MESSAGE") STOP "swarm_message_file_read failed (1)"

    CALL parser_get_next_line(parser, 1, at_end, error)
    at_end = at_end .OR. LEN_TRIM(parser%input_line(1:10))==0
    IF(at_end) RETURN
    READ (parser%input_line(1:40),*) label, length
    IF(TRIM(label)/="msg_length:") STOP "swarm_message_file_read failed (2)"

    DO i=1, length
       ALLOCATE(new_entry)
       CALL swarm_message_entry_file_read(new_entry, parser, at_end, error)
       new_entry%next => msg%root
       msg%root => new_entry
    END DO

    CALL parser_get_next_line(parser, 1, at_end, error)
    at_end = at_end .OR. LEN_TRIM(parser%input_line(1:10))==0
    IF(at_end) RETURN
    IF(TRIM(parser%input_line(1:20))/="END SWARM_MESSAGE") STOP "swarm_message_file_read failed (3)"

  END SUBROUTINE swarm_message_file_read_low


! *****************************************************************************
!> \brief Helper routine for swarm_message_equal
!> \author Ole Schuett
! *****************************************************************************
  FUNCTION swarm_message_equal_oneway(msg1, msg2) RESULT(res)
    TYPE(swarm_message_type), INTENT(IN)     :: msg1, msg2
    LOGICAL                                  :: res

    LOGICAL                                  :: found
    TYPE(message_entry_type), POINTER        :: entry1, entry2

     res = .FALSE.

     !loop over entries of msg1
     entry1 => msg1%root
     DO WHILE(ASSOCIATED(entry1))

        ! finding matching entry in msg2
        entry2 => msg2%root
        found = .FALSE.
        DO WHILE(ASSOCIATED(entry2))
           IF(TRIM(entry2%key) == TRIM(entry1%key)) THEN
              found = .TRUE.
              EXIT
           END IF
           entry2 => entry2%next
        END DO
        IF(.NOT. found) RETURN

        !compare the two entries
        IF(ASSOCIATED(entry1%value_str)) THEN
           IF(.NOT.ASSOCIATED(entry2%value_str)) RETURN
           IF(TRIM(entry1%value_str) /= TRIM(entry2%value_str)) RETURN

        ELSE IF(ASSOCIATED(entry1%value_i4)) THEN
           IF(.NOT. ASSOCIATED(entry2%value_i4)) RETURN
           IF(entry1%value_i4 /= entry2%value_i4) RETURN

        ELSE IF(ASSOCIATED(entry1%value_i8)) THEN
           IF(.NOT. ASSOCIATED(entry2%value_i8)) RETURN
           IF(entry1%value_i8 /= entry2%value_i8) RETURN

        ELSE IF(ASSOCIATED(entry1%value_r4)) THEN
           IF(.NOT. ASSOCIATED(entry2%value_r4)) RETURN
           IF(ABS(entry1%value_r4-entry2%value_r4)>1e-5) RETURN

        ELSE IF(ASSOCIATED(entry1%value_r8)) THEN
           IF(.NOT. ASSOCIATED(entry2%value_r8)) RETURN
           IF(ABS(entry1%value_r8-entry2%value_r8)>1e-10) RETURN

        ELSE IF(ASSOCIATED(entry1%value_i4_1d)) THEN
           IF(.NOT. ASSOCIATED(entry2%value_i4_1d)) RETURN
           IF(ANY(entry1%value_i4_1d /= entry2%value_i4_1d)) RETURN

        ELSE IF(ASSOCIATED(entry1%value_i8_1d)) THEN
           IF(.NOT. ASSOCIATED(entry2%value_i8_1d)) RETURN
           IF(ANY(entry1%value_i8_1d /= entry2%value_i8_1d)) RETURN

        ELSE IF(ASSOCIATED(entry1%value_r4_1d)) THEN
           IF(.NOT. ASSOCIATED(entry2%value_r4_1d)) RETURN
           IF(ANY(ABS(entry1%value_r4_1d-entry2%value_r4_1d)>1e-5)) RETURN

        ELSE IF(ASSOCIATED(entry1%value_r8_1d)) THEN
           IF(.NOT. ASSOCIATED(entry2%value_r8_1d)) RETURN
           IF(ANY(ABS(entry1%value_r8_1d-entry2%value_r8_1d)>1e-10)) RETURN
        ELSE
           STOP "swarm_message_equal_oneway: no value ASSOCIATED"
        END IF

        entry1 => entry1%next
     END DO

     ! if we reach this point no differences were found
     res = .TRUE.
  END FUNCTION swarm_message_equal_oneway


! *****************************************************************************
!> \brief Helper routine for swarm_message_mpi_send.
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE swarm_message_entry_mpi_send(ENTRY, group, dest, tag)
    TYPE(message_entry_type), INTENT(IN)     :: ENTRY
    INTEGER, INTENT(IN)                      :: group, dest, tag

    INTEGER, &
      DIMENSION(default_string_length)       :: value_str_arr
    INTEGER, DIMENSION(key_length)           :: key_arr

    key_arr = str2iarr(entry%key)
    CALL mp_send(key_arr, dest, tag, group)

    IF(ASSOCIATED(entry%value_i4)) THEN
       CALL mp_send(1, dest, tag, group)
       CALL mp_send(entry%value_i4, dest, tag, group)

    ELSE IF(ASSOCIATED(entry%value_i8)) THEN
       CALL mp_send(2, dest, tag, group)
       CALL mp_send(entry%value_i8, dest, tag, group)

    ELSE IF(ASSOCIATED(entry%value_r4)) THEN
       CALL mp_send(3, dest, tag, group)
       CALL mp_send(entry%value_r4, dest, tag, group)

    ELSE IF(ASSOCIATED(entry%value_r8)) THEN
       CALL mp_send(4, dest, tag, group)
       CALL mp_send(entry%value_r8, dest, tag, group)

    ELSE IF(ASSOCIATED(entry%value_i4_1d)) THEN
       CALL mp_send(5, dest, tag, group)
       CALL mp_send(SIZE(entry%value_i4_1d), dest, tag, group)
       CALL mp_send(entry%value_i4_1d, dest, tag, group)

    ELSE IF(ASSOCIATED(entry%value_i8_1d)) THEN
       CALL mp_send(6, dest, tag, group)
       CALL mp_send(SIZE(entry%value_i8_1d), dest, tag, group)
       CALL mp_send(entry%value_i8_1d, dest, tag, group)

    ELSE IF(ASSOCIATED(entry%value_r4_1d)) THEN
       CALL mp_send(7, dest, tag, group)
       CALL mp_send(SIZE(entry%value_r4_1d), dest, tag, group)
       CALL mp_send(entry%value_r4_1d, dest, tag, group)

    ELSE IF(ASSOCIATED(entry%value_r8_1d)) THEN
       CALL mp_send(8, dest, tag, group)
       CALL mp_send(SIZE(entry%value_r8_1d), dest, tag, group)
       CALL mp_send(entry%value_r8_1d, dest, tag, group)

    ELSE IF(ASSOCIATED(entry%value_str)) THEN
       CALL mp_send(9, dest, tag, group)
       value_str_arr = str2iarr(entry%value_str)
       CALL mp_send(value_str_arr, dest, tag, group)
    ELSE
       CALL mp_abort()
       STOP "swarm_message_entry_mpi_send: no value ASSOCIATED"
    END IF
  END SUBROUTINE swarm_message_entry_mpi_send


! *****************************************************************************
!> \brief Helper routine for swarm_message_mpi_recv.
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE swarm_message_entry_mpi_recv(ENTRY, group, src, tag)
    TYPE(message_entry_type), INTENT(INOUT)  :: ENTRY
    INTEGER, INTENT(IN)                      :: group
    INTEGER, INTENT(INOUT)                   :: src, tag

    INTEGER                                  :: datatype, s
    INTEGER, &
      DIMENSION(default_string_length)       :: value_str_arr
    INTEGER, DIMENSION(key_length)           :: key_arr

    CALL mp_recv(key_arr, src, tag, group)
    entry%key = iarr2str(key_arr)

    CALL mp_recv(datatype, src, tag, group)

    SELECT CASE(datatype)
    CASE(1)
       ALLOCATE(entry%value_i4)
       CALL mp_recv(entry%value_i4, src, tag, group)
    CASE(2)
       ALLOCATE(entry%value_i8)
       CALL mp_recv(entry%value_i8, src, tag, group)
    CASE(3)
       ALLOCATE(entry%value_r4)
       CALL mp_recv(entry%value_r4, src, tag, group)
    CASE(4)
       ALLOCATE(entry%value_r8)
       CALL mp_recv(entry%value_r8, src, tag, group)

    CASE(5)
       CALL mp_recv(s, src, tag, group)
       ALLOCATE(entry%value_i4_1d(s))
       CALL mp_recv(entry%value_i4_1d, src, tag, group)
    CASE(6)
       CALL mp_recv(s, src, tag, group)
       ALLOCATE(entry%value_i8_1d(s))
       CALL mp_recv(entry%value_i8_1d, src, tag, group)
    CASE(7)
       CALL mp_recv(s, src, tag, group)
       ALLOCATE(entry%value_r4_1d(s))
       CALL mp_recv(entry%value_r4_1d, src, tag, group)
    CASE(8)
       CALL mp_recv(s, src, tag, group)
       ALLOCATE(entry%value_r8_1d(s))
       CALL mp_recv(entry%value_r8_1d, src, tag, group)
    CASE(9)
       ALLOCATE(entry%value_str)
       CALL mp_recv(value_str_arr, src, tag, group)
       entry%value_str = iarr2str(value_str_arr)
    CASE DEFAULT
       STOP "swarm_message_entry_mpi_recv: unkown datatype"
    END SELECT
  END SUBROUTINE swarm_message_entry_mpi_recv


! *****************************************************************************
!> \brief Helper routine for swarm_message_mpi_bcast.
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE swarm_message_entry_mpi_bcast(ENTRY, src, group, mepos)
    TYPE(message_entry_type), INTENT(INOUT)  :: ENTRY
    INTEGER, INTENT(IN)                      :: src, group, mepos

    INTEGER                                  :: datasize, datatype
    INTEGER, &
      DIMENSION(default_string_length)       :: value_str_arr
    INTEGER, DIMENSION(key_length)           :: key_arr

    IF(src==mepos) key_arr = str2iarr(entry%key)
    CALL mp_bcast(key_arr, src, group)
    IF(src/=mepos) entry%key = iarr2str(key_arr)

    IF(src==mepos) THEN
       datasize = 1
       IF(ASSOCIATED(entry%value_i4)) THEN
          datatype = 1
       ELSE IF(ASSOCIATED(entry%value_i8)) THEN
          datatype = 2
       ELSE IF(ASSOCIATED(entry%value_r4)) THEN
          datatype = 3
       ELSE IF(ASSOCIATED(entry%value_r8)) THEN
          datatype = 4
       ELSE IF(ASSOCIATED(entry%value_i4_1d)) THEN
          datatype = 5
          datasize = SIZE(entry%value_i4_1d)
       ELSE IF(ASSOCIATED(entry%value_i8_1d)) THEN
          datatype = 6
          datasize = SIZE(entry%value_i8_1d)
       ELSE IF(ASSOCIATED(entry%value_r4_1d)) THEN
          datatype = 7
          datasize = SIZE(entry%value_r4_1d)
       ELSE IF(ASSOCIATED(entry%value_r8_1d)) THEN
          datatype = 8
          datasize = SIZE(entry%value_r8_1d)
       ELSE IF(ASSOCIATED(entry%value_str)) THEN
          datatype = 9
       ELSE
          STOP "swarm_message_entry_mpi_bcast: no value ASSOCIATED"
       END IF
    END IF
    CALL mp_bcast(datatype, src, group)
    CALL mp_bcast(datasize, src, group)

    SELECT CASE(datatype)
    CASE(1)
       IF(src/=mepos) ALLOCATE(entry%value_i4)
       CALL mp_bcast(entry%value_i4, src, group)
    CASE(2)
       IF(src/=mepos) ALLOCATE(entry%value_i8)
       CALL mp_bcast(entry%value_i8, src, group)
    CASE(3)
       IF(src/=mepos) ALLOCATE(entry%value_r4)
       CALL mp_bcast(entry%value_r4, src, group)
    CASE(4)
       IF(src/=mepos) ALLOCATE(entry%value_r8)
       CALL mp_bcast(entry%value_r8, src, group)
    CASE(5)
       IF(src/=mepos) ALLOCATE(entry%value_i4_1d(datasize))
       CALL mp_bcast(entry%value_i4_1d, src, group)
    CASE(6)
       IF(src/=mepos) ALLOCATE(entry%value_i8_1d(datasize))
       CALL mp_bcast(entry%value_i8_1d, src, group)
    CASE(7)
       IF(src/=mepos) ALLOCATE(entry%value_r4_1d(datasize))
       CALL mp_bcast(entry%value_r4_1d, src, group)
    CASE(8)
       IF(src/=mepos) ALLOCATE(entry%value_r8_1d(datasize))
       CALL mp_bcast(entry%value_r8_1d, src, group)
    CASE(9)
       IF(src==mepos) value_str_arr = str2iarr(entry%value_str)
       CALL mp_bcast(value_str_arr, src, group)
       IF(src/=mepos) THEN
          ALLOCATE(entry%value_str)
          entry%value_str = iarr2str(value_str_arr)
       END IF
    CASE DEFAULT
       STOP "swarm_message_entry_mpi_bcast: unkown datatype"
    END SELECT

  END SUBROUTINE swarm_message_entry_mpi_bcast


! *****************************************************************************
!> \brief Helper routine for swarm_message_file_write.
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE swarm_message_entry_file_write(ENTRY, unit)
    TYPE(message_entry_type), INTENT(IN)     :: ENTRY
    INTEGER, INTENT(IN)                      :: unit

    INTEGER                                  :: i

    WRITE(unit,"(A,A)") "key: ",  entry%key
    IF(ASSOCIATED(entry%value_i4)) THEN
       WRITE(unit,"(A)") "datatype: i4"
       WRITE(unit,"(A,I10)") "value: ", entry%value_i4

    ELSE IF(ASSOCIATED(entry%value_i8)) THEN
       WRITE(unit,"(A)") "datatype: i8"
       WRITE(unit,"(A,I20)") "value: ", entry%value_i8

    ELSE IF(ASSOCIATED(entry%value_r4)) THEN
       WRITE(unit,"(A)") "datatype: r4"
       WRITE(unit,"(A,E30.20)") "value: ", entry%value_r4

    ELSE IF(ASSOCIATED(entry%value_r8)) THEN
       WRITE(unit,"(A)") "datatype: r8"
       WRITE(unit,"(A,E30.20)") "value: ", entry%value_r8

    ELSE IF(ASSOCIATED(entry%value_str)) THEN
       WRITE(unit,"(A)") "datatype: str"
       WRITE(unit,"(A,A)") "value: ", entry%value_str

    ELSE IF(ASSOCIATED(entry%value_i4_1d)) THEN
       WRITE(unit,"(A)") "datatype: i4_1d"
       WRITE(unit,"(A,I10)") "size: ", SIZE(entry%value_i4_1d)
       DO i=1, SIZE(entry%value_i4_1d)
         WRITE(unit,*) entry%value_i4_1d(i)
       END DO

    ELSE IF(ASSOCIATED(entry%value_i8_1d)) THEN
       WRITE(unit,"(A)") "datatype: i8_1d"
       WRITE(unit,"(A,I20)") "size: ", SIZE(entry%value_i8_1d)
       DO i=1, SIZE(entry%value_i8_1d)
         WRITE(unit,*) entry%value_i8_1d(i)
       END DO

    ELSE IF(ASSOCIATED(entry%value_r4_1d)) THEN
       WRITE(unit,"(A)") "datatype: r4_1d"
       WRITE(unit,"(A,I8)") "size: ", SIZE(entry%value_r4_1d)
       DO i=1, SIZE(entry%value_r4_1d)
         WRITE(unit,"(1X,E30.20)") entry%value_r4_1d(i)
       END DO

    ELSE IF(ASSOCIATED(entry%value_r8_1d)) THEN
       WRITE(unit,"(A)") "datatype: r8_1d"
       WRITE(unit,"(A,I8)") "size: ", SIZE(entry%value_r8_1d)
       DO i=1, SIZE(entry%value_r8_1d)
         WRITE(unit,"(1X,E30.20)") entry%value_r8_1d(i)
       END DO

    ELSE
       CALL mp_abort()
       STOP "swarm_message_entry_file_write: no value ASSOCIATED"
    END IF
  END SUBROUTINE swarm_message_entry_file_write


! *****************************************************************************
!> \brief Helper routine for swarm_message_file_read.
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE swarm_message_entry_file_read(ENTRY, parser, at_end, error)
    TYPE(message_entry_type), INTENT(INOUT)  :: ENTRY
    TYPE(cp_parser_type), POINTER            :: parser
    LOGICAL, INTENT(INOUT)                   :: at_end
    TYPE(cp_error_type)                      :: error

    CHARACTER(LEN=15)                        :: datatype, label
    INTEGER                                  :: arr_size, i
    LOGICAL                                  :: is_scalar

    CALL parser_get_next_line(parser, 1, at_end, error)
    at_end = at_end .OR. LEN_TRIM(parser%input_line(1:10))==0
    IF(at_end) RETURN
    READ (parser%input_line(1:key_length+10),*) label, entry%key
    IF(TRIM(label)/="key:") STOP "swarm_message_entry_file_read failed (1)"

    CALL parser_get_next_line(parser, 1, at_end, error)
    at_end = at_end .OR. LEN_TRIM(parser%input_line(1:10))==0
    IF(at_end) RETURN
    READ (parser%input_line(1:30),*) label, datatype
    IF(TRIM(label)/="datatype:") STOP "swarm_message_entry_file_read failed (2)"

    CALL parser_get_next_line(parser, 1, at_end, error)
    at_end = at_end .OR. LEN_TRIM(parser%input_line(1:10))==0
    IF(at_end) RETURN

    is_scalar=.TRUE.
    SELECT CASE (TRIM(datatype))
      CASE("i4")
         ALLOCATE(entry%value_i4)
         READ (parser%input_line(1:40),*) label, entry%value_i4
      CASE("i8")
         ALLOCATE(entry%value_i8)
         READ (parser%input_line(1:40),*) label, entry%value_i8
      CASE("r4")
         ALLOCATE(entry%value_r4)
         READ (parser%input_line(1:40),*) label, entry%value_r4
      CASE("r8")
         ALLOCATE(entry%value_r8)
         READ (parser%input_line(1:40),*) label, entry%value_r8
      CASE("str")
         ALLOCATE(entry%value_str)
         READ (parser%input_line(1:40),*) label, entry%value_str
      CASE DEFAULT
         is_scalar = .FALSE.
    END SELECT

    IF(is_scalar) THEN
       IF(TRIM(label)/="value:") STOP "swarm_message_entry_file_read failed (3)"
       RETURN
    END IF

    ! musst be an array-datatype
    READ (parser%input_line(1:30),*) label, arr_size
    IF(TRIM(label)/="size:") STOP "swarm_message_entry_file_read failed (4)"

    SELECT CASE (TRIM(datatype))
      CASE("i4_1d")
        ALLOCATE(entry%value_i4_1d(arr_size))
      CASE("i8_1d")
       ALLOCATE(entry%value_i8_1d(arr_size))
      CASE("r4_1d")
        ALLOCATE(entry%value_r4_1d(arr_size))
      CASE("r8_1d")
        ALLOCATE(entry%value_r8_1d(arr_size))
      CASE DEFAULT
       STOP "swarm_message_entry_file_read: unkown datatype"
    END SELECT

    DO i=1, arr_size
      CALL parser_get_next_line(parser, 1, at_end, error)
      at_end = at_end .OR. LEN_TRIM(parser%input_line(1:10))==0
      IF(at_end) RETURN

      !Numbers were written with at most 31 characters.
      SELECT CASE (TRIM(datatype))
        CASE("i4_1d")
          READ (parser%input_line(1:31),*) entry%value_i4_1d(i)
        CASE("i8_1d")
          READ (parser%input_line(1:31),*) entry%value_i8_1d(i)
        CASE("r4_1d")
          READ (parser%input_line(1:31),*) entry%value_r4_1d(i)
        CASE("r8_1d")
          READ (parser%input_line(1:31),*) entry%value_r8_1d(i)
        CASE DEFAULT
          STOP "swarm_message_entry_file_read: unkown datatype"
      END SELECT
    END DO


  END SUBROUTINE swarm_message_entry_file_read


! *****************************************************************************
!> \brief Helper routine, converts a string into an integer-array
!> \author Ole Schuett
! *****************************************************************************
  PURE FUNCTION str2iarr(str) RESULT(arr)
    CHARACTER(LEN=*), INTENT(IN)             :: str
    INTEGER, DIMENSION(LEN(str))             :: arr

    INTEGER                                  :: i

    DO i=1,LEN(str)
      arr(i) = ICHAR(str(i:i))
    ENDDO
  END FUNCTION str2iarr


! *****************************************************************************
!> \brief Helper routine, converts an integer-array into a string
!> \author Ole Schuett
! *****************************************************************************
  PURE FUNCTION iarr2str(arr) RESULT(str)
    INTEGER, DIMENSION(:), INTENT(IN)        :: arr
    CHARACTER(LEN=SIZE(arr))                 :: str

    INTEGER                                  :: i

    DO i=1,SIZE(arr)
      str(i:i) = CHAR(arr(i))
    ENDDO
  END FUNCTION iarr2str


#include "swarm_message_str.f90"
#include "swarm_message_i4.f90"
#include "swarm_message_i8.f90"
#include "swarm_message_r4.f90"
#include "swarm_message_r8.f90"
#include "swarm_message_i4_1d.f90"
#include "swarm_message_i8_1d.f90"
#include "swarm_message_r4_1d.f90"
#include "swarm_message_r8_1d.f90"


END MODULE swarm_message


