!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2003  CP2K developers group
!-----------------------------------------------------------------------------!
!!****** cp2k/string_utilities [1.0] *
!!
!!   NAME
!!     string_utilities
!!
!!   FUNCTION
!!     Utilities for string manipulations
!!
!!   AUTHOR
!!     MK & JGH
!!
!!   MODIFICATION HISTORY
!!     Adapted compress and uppercase for use in CP2K (JGH)
!!     string_to_integer and integer_to_string added (06.02.2001, MK)
!!     Cleaned (04.01.2004,MK)
!!
!!   SOURCE
!******************************************************************************

MODULE string_utilities

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ascii_to_string,&
            compress,&
            integer_to_string,&
            make_tuple,&
            str_comp,&
            string_to_ascii,&
            uppercase,&
            xstring,str_search

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE ascii_to_string(nascii,string)

!   Purpose: Convert a sequence of integer numbers (ASCII code) to a string.
!            Blanks are inserted for invalid ASCII code numbers.

!   History: - Creation (19.10.2000,MK)

!   ***************************************************************************

    INTEGER, DIMENSION(:), INTENT(IN)        :: nascii
    CHARACTER(LEN=*), INTENT(OUT)            :: string

    INTEGER                                  :: i

!   ---------------------------------------------------------------------------

    string = ""

    DO i=1,MIN(LEN(string),SIZE(nascii))
      IF ((nascii(i) >= 0).AND.(nascii(i) <= 127)) THEN
        string(i:i) = CHAR(nascii(i))
      ELSE
        string(i:i) = " "
      END IF
    END DO

  END SUBROUTINE ascii_to_string

! *****************************************************************************

  SUBROUTINE compress(string,full)

!   Purpose: Eliminate multiple space characters in a string.
!            If full is .TRUE., then all spaces are eliminated.

!   History: - Creation (23.06.1998,MK)

!   ***************************************************************************

    CHARACTER(LEN=*), INTENT(INOUT)          :: string
    LOGICAL, INTENT(IN), OPTIONAL            :: full

    INTEGER                                  :: i, z
    LOGICAL                                  :: remove_all

!   ---------------------------------------------------------------------------

    IF (PRESENT(full)) THEN
      remove_all = full
    ELSE
      remove_all = .FALSE.
    END IF

    z = 1

    DO i=1,LEN_TRIM(string)
      IF ((z == 1).OR.remove_all) THEN
        IF (string(i:i) /= " ") THEN
          string(z:z) = string(i:i)
          z = z + 1
        END IF
      ELSE
        IF ((string(i:i) /= " ").OR.(string(z-1:z-1) /= " ")) THEN
          string(z:z) = string(i:i)
          z = z + 1
        END IF
      END IF
    END DO

    string(z:) = ""

  END SUBROUTINE compress

! *****************************************************************************

  SUBROUTINE integer_to_string(inumber,string)

!   Purpose: Converts an integer number to a string.
!            The WRITE statement will return an error message, if the number of
!            digits of the integer number is larger the than the length of the
!            supplied string.

!   History: - Creation (05.01.2004,MK)

!   ***************************************************************************

    INTEGER, INTENT(IN)                      :: inumber
    CHARACTER(LEN=*), INTENT(OUT)            :: string

!   ---------------------------------------------------------------------------

    string = ""
    WRITE (UNIT=string,FMT=*) inumber
    string = ADJUSTL(string)

  END SUBROUTINE integer_to_string

! *****************************************************************************

  SUBROUTINE string_to_ascii(string,nascii)

!   Purpose: Convert a string to sequence of integer numbers.

!   History: - Creation (19.10.2000,MK)

!   ***************************************************************************

    CHARACTER(LEN=*), INTENT(IN)             :: string
    INTEGER, DIMENSION(:), INTENT(OUT)       :: nascii

    INTEGER                                  :: i

!   ---------------------------------------------------------------------------

    nascii(:) = 0

    DO i=1,MIN(LEN(string),SIZE(nascii))
      nascii(i) = ICHAR(string(i:i))
    END DO

  END SUBROUTINE string_to_ascii

! *****************************************************************************

  SUBROUTINE uppercase(string)

!   Purpose: Convert all lower case characters in a string to upper case.

!   History: - Creation (22.06.1998,MK)

!   ***************************************************************************

    CHARACTER(LEN=*), INTENT(INOUT)          :: string

    INTEGER                                  :: i, iascii

!   ---------------------------------------------------------------------------

    DO i=1,LEN_TRIM(string)
      iascii = ICHAR(string(i:i))
      IF ((iascii >= 97).AND.(iascii <= 122)) THEN
        string(i:i) = CHAR(iascii - 32)
      END IF
    END DO

  END SUBROUTINE uppercase

! *****************************************************************************

SUBROUTINE xstring(string,ia,ib)

    CHARACTER(LEN=*), INTENT(IN)             :: string
    INTEGER, INTENT(OUT)                     :: ia, ib

!------------------------------------------------------------------------------

  ia = 1
  ib = LEN_TRIM(string)
  IF (ib>0) THEN
     DO WHILE (string(ia:ia)==' ')
        ia = ia + 1
     END DO
  END IF

END SUBROUTINE xstring

!******************************************************************************

SUBROUTINE make_tuple(int,nt,na,tuple)

    INTEGER, INTENT(IN)                      :: INT( :, : ), nt, na
    CHARACTER(LEN=*), INTENT(OUT)            :: tuple( : )

    INTEGER                                  :: i, nm

  nm = MAXVAL(INT(1:nt,1:na))
  SELECT CASE (nt)
  CASE DEFAULT
     STOP 'make_tuple: case not programmed'
  CASE (1)
     IF (nm<100) THEN
        DO i = 1, na
           WRITE (tuple(i),'(A,I2,A )' ) '[', INT(1,i), ']'
        END DO
     ELSE
        DO i = 1, na
           WRITE (tuple(i),'(A,I4,A )' ) '[', INT(1,i), ']'
        END DO
     END IF
  CASE (2)
     IF (nm<100) THEN
        DO i = 1, na
           WRITE (tuple(i),'(A,I2,A,I2,A )' ) '[', INT(1,i), '-', INT(2,i), &
                ']'
        END DO
     ELSE
        DO i = 1, na
           WRITE (tuple(i),'(A,I4,A,I4,A )' ) '[', INT(1,i), '-', INT(2,i), &
                ']'
        END DO
     END IF
  CASE (3)
     IF (nm<100) THEN
        DO i = 1, na
           WRITE (tuple(i),'(A,I2,A,I2,A,I2,A )' ) '[', INT(1,i), '-', &
                INT(2,i), '-', INT(3,i), ']'
        END DO
     ELSE
        DO i = 1, na
           WRITE (tuple(i),'(A,I4,A,I4,A,I4,A )' ) '[', INT(1,i), '-', &
                INT(2,i), '-', INT(3,i), ']'
        END DO
     END IF
  CASE (4)
     IF (nm<100) THEN
        DO i = 1, na
           WRITE (tuple(i),'(A,I2,A,I2,A,I2,A,I2,A )' ) '[', INT(1,i), '-', &
                INT(2,i), '-', INT(3,i), '-', INT(4,i), ']'
        END DO
     ELSE
        DO i = 1, na
           WRITE (tuple(i),'(A,I4,A,I4,A,I4,A,I4,A )' ) '[', INT(1,i), '-', &
                INT(2,i), '-', INT(3,i), '-', INT(4,i), ']'
        END DO
     END IF
  END SELECT

END SUBROUTINE make_tuple

!******************************************************************************

FUNCTION str_comp(str1,str2) RESULT (equal)

    CHARACTER(LEN=*), INTENT(IN)             :: str1, str2
    LOGICAL                                  :: equal

    INTEGER                                  :: i1, i2, j1, j2

  i1 = 0
  i2 = 0
  j1 = 0
  j2 = 0
  CALL xstring(str1,i1,i2)
  CALL xstring(str2,j1,j2)
  equal = (str1(i1:i2)==str2(j1:j2))
END FUNCTION str_comp

!******************************************************************************

FUNCTION str_search(str1,n,str2) RESULT (pos)
    CHARACTER(LEN=*), INTENT(IN)             :: str1( : )
    INTEGER, INTENT(IN)                      :: n
    CHARACTER(LEN=*), INTENT(IN)             :: str2
    INTEGER                                  :: pos

    INTEGER                                  :: i

  pos = 0
  DO i = 1, n
     IF (str_comp(str1(i),str2)) THEN
        pos = i
        EXIT
     END IF
  END DO
END FUNCTION str_search

!******************************************************************************

END MODULE string_utilities
