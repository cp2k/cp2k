!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/string_utilities [1.0] *
!!
!!   NAME
!!     string_utilities
!!
!!   FUNCTION
!!     Utilities for string manipulations.
!!
!!   AUTHOR
!!     Matthias Krack & JGH
!!
!!   MODIFICATION HISTORY
!!     Adapted compress and uppercase for use in CP2K (JGH)
!!     string_to_integer and integer_to_string added (06.02.2001, MK)
!!
!!   SOURCE
!******************************************************************************

MODULE string_utilities

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: compress,uppercase,xstring,str_search,str_comp,make_tuple,&
            integer_to_string,string_to_integer


!!*****
! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE compress(string)

!   Purpose: Eliminate multiple space characters in a string.

!   History: - Creation (23.06.1998, Matthias Krack)

!   ***************************************************************************

    CHARACTER(LEN=*), INTENT(INOUT) :: string

!   *** Local variables ***

    INTEGER :: i,z

!   ---------------------------------------------------------------------------

    z = 1

    DO i=1,LEN_TRIM(string)
      IF (z == 1) THEN
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

  SUBROUTINE integer_to_string(nascii,string)

!   Purpose: Convert a sequence of integer numbers to a string.

!   History: - Creation (19.10.2000, Matthias Krack)

!   ***************************************************************************

    CHARACTER(LEN=*), INTENT(OUT)     :: string
    INTEGER, DIMENSION(:), INTENT(IN) :: nascii

!   *** Local variables ***

    INTEGER :: i

!   ---------------------------------------------------------------------------

    string = ""

    DO i=1,MIN(LEN(string),SIZE(nascii))
      IF ((nascii(i) >= 0).AND.(nascii(i) <= 127)) THEN
        string(i:i) = CHAR(nascii(i))
      ELSE
        string(i:i) = " "
      END IF
    END DO

  END SUBROUTINE integer_to_string

! *****************************************************************************

  SUBROUTINE string_to_integer(string,nascii)

!   Purpose: Convert a string to sequence of integer numbers.

!   History: - Creation (19.10.2000, Matthias Krack)

!   ***************************************************************************

    CHARACTER(LEN=*), INTENT(IN)       :: string
    INTEGER, DIMENSION(:), INTENT(OUT) :: nascii

!   *** Local variables ***

    INTEGER :: i

!   ---------------------------------------------------------------------------

    nascii(:) = 0

    DO i=1,MIN(LEN(string),SIZE(nascii))
      nascii(i) = ICHAR(string(i:i))
    END DO

  END SUBROUTINE string_to_integer

! *****************************************************************************

  SUBROUTINE uppercase(string)

!   Purpose: Convert all lower case characters in a string to upper case.

!   History: - Creation (22.06.1998, Matthias Krack)

!   ***************************************************************************

    CHARACTER(LEN=*), INTENT(INOUT) :: string

!   *** Local variables ***

    INTEGER :: i,iascii

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

  IMPLICIT NONE

  CHARACTER ( LEN = * ), INTENT ( IN ) :: string
  INTEGER, INTENT ( OUT ) :: ia, ib

!------------------------------------------------------------------------------

  ia = 1
  ib = len_trim(string)
  IF (ib>0) THEN
     DO WHILE (string(ia:ia)==' ')
        ia = ia + 1
     END DO
  END IF

END SUBROUTINE xstring

!******************************************************************************

SUBROUTINE make_tuple(int,nt,na,tuple)
  IMPLICIT NONE
  INTEGER, INTENT ( IN ) :: int ( :, : ), nt, na
  INTEGER :: i, nm
  CHARACTER ( LEN = * ), INTENT ( OUT ) :: tuple ( : )

  nm = maxval(int(1:nt,1:na))
  SELECT CASE (nt)
  CASE DEFAULT
     STOP 'make_tuple: case not programmed'
  CASE (1)
     IF (nm<100) THEN
        DO i = 1, na
           WRITE (tuple(i),'(A,I2,A )' ) '[', int(1,i), ']'
        END DO
     ELSE
        DO i = 1, na
           WRITE (tuple(i),'(A,I4,A )' ) '[', int(1,i), ']'
        END DO
     END IF
  CASE (2)
     IF (nm<100) THEN
        DO i = 1, na
           WRITE (tuple(i),'(A,I2,A,I2,A )' ) '[', int(1,i), '-', int(2,i), &
                ']'
        END DO
     ELSE
        DO i = 1, na
           WRITE (tuple(i),'(A,I4,A,I4,A )' ) '[', int(1,i), '-', int(2,i), &
                ']'
        END DO
     END IF
  CASE (3)
     IF (nm<100) THEN
        DO i = 1, na
           WRITE (tuple(i),'(A,I2,A,I2,A,I2,A )' ) '[', int(1,i), '-', &
                int(2,i), '-', int(3,i), ']'
        END DO
     ELSE
        DO i = 1, na
           WRITE (tuple(i),'(A,I4,A,I4,A,I4,A )' ) '[', int(1,i), '-', &
                int(2,i), '-', int(3,i), ']'
        END DO
     END IF
  CASE (4)
     IF (nm<100) THEN
        DO i = 1, na
           WRITE (tuple(i),'(A,I2,A,I2,A,I2,A,I2,A )' ) '[', int(1,i), '-', &
                int(2,i), '-', int(3,i), '-', int(4,i), ']'
        END DO
     ELSE
        DO i = 1, na
           WRITE (tuple(i),'(A,I4,A,I4,A,I4,A,I4,A )' ) '[', int(1,i), '-', &
                int(2,i), '-', int(3,i), '-', int(4,i), ']'
        END DO
     END IF
  END SELECT

END SUBROUTINE make_tuple

!******************************************************************************

FUNCTION str_comp(str1,str2) RESULT (equal)
  IMPLICIT NONE
  CHARACTER ( LEN = * ), INTENT ( IN ) :: str1, str2
  LOGICAL equal
  INTEGER :: i1, i2, j1, j2

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
  CHARACTER ( LEN = * ), INTENT ( IN ) :: str1 ( : ), str2
  INTEGER, INTENT ( IN ) :: n
  INTEGER :: pos, i

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
