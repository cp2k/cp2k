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
  USE kinds,                           ONLY: dbl
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: compress,uppercase,xstring,str_search,str_comp,make_tuple,&
            integer_to_string,string_to_integer, compress_all

  INTERFACE integer_to_string
     MODULE PROCEDURE integer_to_string_n, integer_to_string_1
  END INTERFACE

  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN="string_utilities"

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

  SUBROUTINE integer_to_string_n(nascii,string)

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

  END SUBROUTINE integer_to_string_n

! *****************************************************************************

  SUBROUTINE integer_to_string_1(zahl,string)

!   Purpose: Convert an integer to a string.

!   History: - Creation (29.4.2001, JGH)

!   ***************************************************************************

    CHARACTER(LEN=*), INTENT(OUT)     :: string
    INTEGER, INTENT(IN) :: zahl

!   *** Local variables ***

    INTEGER :: i, z, n
    INTEGER,ALLOCATABLE,DIMENSION(:) :: nascii

!   ---------------------------------------------------------------------------

    string = ""
    n = INT ( LOG ( REAL ( zahl, dbl ) ) / LOG ( 10._dbl ) ) + 1
    ALLOCATE ( nascii ( n ) )
    z = zahl
    DO i = 1, n
      nascii ( i ) = z / 10**(n-i)
      z = z - nascii ( i )*10**(n-i)
      nascii ( i ) = nascii ( i ) + ICHAR("0")
    END DO
    CALL integer_to_string_n(nascii,string)
    DEALLOCATE ( nascii )

  END SUBROUTINE integer_to_string_1

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
  ib = LEN_TRIM(string)
  IF (ib>0) THEN
     DO WHILE (string(ia:ia)==' ')
        ia = ia + 1
     END DO
  END IF

END SUBROUTINE xstring

!******************************************************************************

SUBROUTINE make_tuple(int,nt,na,tuple)
  IMPLICIT NONE
  INTEGER, INTENT ( IN ) :: INT ( :, : ), nt, na
  INTEGER :: i, nm
  CHARACTER ( LEN = * ), INTENT ( OUT ) :: tuple ( : )

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

!!****f* string_utilities/compress_all [1.0] *
!!
!!   NAME
!!     compress_all
!!
!!   SYNOPSIS
!!     Subroutine compress_all(string)
!!       Character(Len=*), Intent (INOUT):: string
!!     End Subroutine compress_all
!!
!!   FUNCTION
!!     removes all the blancs from the given string, putting them at end
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE compress_all(string)
  CHARACTER(len=*), INTENT(inout) :: string
  
  CHARACTER(len=*), PARAMETER :: routineN='compress_all',&
        routineP=moduleN//':'//routineN
  INTEGER :: i, pos2

  pos2=1
  DO i=1,LEN_TRIM(string)
     IF (string(i:i)/=' ') THEN
        string(pos2:pos2)=string(i:i)
        pos2=pos2+1
     END IF
  END DO
  IF (pos2<LEN(string)) THEN
     string(pos2:LEN(string))=""
  END IF

END SUBROUTINE compress_all
!***************************************************************************

END MODULE string_utilities
