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

  USE glob_matching,                   ONLY: pattern_match=>string_match

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=1), PARAMETER :: backslash = '\\'
  CHARACTER(len=1), PARAMETER :: star      = '*'
  CHARACTER(len=1), PARAMETER :: question  = '?'


  PUBLIC :: ascii_to_string,&
            compress,&
            integer_to_string,&
            make_tuple,&
            str_comp,&
            string_to_ascii,&
            substitute_special_xml_tokens,&
            uppercase,&
            xstring,str_search, s2a, &
            pattern_match,typo_match,backslash,star,question

  INTERFACE s2a
    MODULE PROCEDURE s2a_1,s2a_2,s2a_3,s2a_4,s2a_5,s2a_6, &
                     s2a_7,s2a_8,s2a_9,s2a_10,s2a_11,s2a_12, s2a_13, s2a_14  ! should be clear how to add more
  END INTERFACE 


! *****************************************************************************

CONTAINS
!!****f* string_utilities/typo_match *
!!
!!   NAME
!!    typo_match
!!
!!   FUNCTION
!!    returns true if typo_string equals string apart from a few typos. It is case sensitive, apart from typos.
!! returns true if there is any string match for two 'stars'/'questions' walking the typo_string and string
!!
!!   NOTES
!!    could maybe be made a bit smarter
!!
!!   MODIFICATION HISTORY
!!     02.2006 created [Joost VandeVondele]
!!
!!   SOURCE
!!*** **********************************************************************
FUNCTION typo_match(string,typo_string) RESULT(match)
    CHARACTER(LEN=*), INTENT(IN)             :: string, typo_string
    LOGICAL                                  :: match

    CHARACTER(LEN=1)                         :: kind
    CHARACTER(LEN=LEN(string))               :: tmp2
    CHARACTER(LEN=LEN(typo_string))          :: tmp
    INTEGER                                  :: i, j

    match=.FALSE.
    IF (LEN_TRIM(typo_string).LE.4) THEN
       kind=question 
    ELSE
       kind=star
    ENDIF
    DO i=1,LEN_TRIM(typo_string)
     DO j=i,LEN_TRIM(typo_string)
      tmp=typo_string
      tmp(i:i)=kind
      tmp(j:j)=kind
      IF (i==j .AND. LEN_TRIM(typo_string)>1 ) tmp(i:i)=star
   !dbg   IF (pattern_match(string=string,pattern=tmp)) write(6,*) "->"//string//"<-","->"//tmp//"<-"
      match=match.OR.pattern_match(string=string,pattern=tmp)
     ENDDO
    ENDDO
    IF (LEN_TRIM(string).LE.4) THEN
       kind=question 
    ELSE
       kind=star
    ENDIF
    DO i=1,LEN_TRIM(string)
     DO j=i,LEN_TRIM(string)
      tmp2=string
      tmp2(i:i)=kind
      tmp2(j:j)=kind
      IF (i==j .AND. LEN_TRIM(typo_string)>1 ) tmp(i:i)=star
   !dbg   IF (pattern_match(string=typo_string,pattern=tmp2)) write(6,*) "->"//typo_string//"<-","->"//TRIM(tmp2)//"<-"
      match=match.OR.pattern_match(string=typo_string,pattern=tmp2)
     ENDDO
    ENDDO

END FUNCTION typo_match

!!****f* string_utilities/s2a *
!!
!!   NAME
!!    s2a_x
!!
!!   FUNCTION
!!    converts a bunch of strings of different length to an array of strings of the same length
!!   NOTES
!!    can be used instead of the illegal (/"12","1234"/) generating s2a("12","1234").EQ.(/"12  ","1234"/)
!!
!!   MODIFICATION HISTORY
!!     11.2004 created [Joost VandeVondele ]
!!
!!   SOURCE
!!*** **********************************************************************
PURE FUNCTION s2a_1(s1) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN)             :: s1
    CHARACTER(LEN=1000), DIMENSION(1)        :: a

  a(1)=s1
END FUNCTION
PURE FUNCTION s2a_2(s1,s2) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN)             :: s1, s2
    CHARACTER(LEN=1000), DIMENSION(2)        :: a

  a(1)=s1; a(2)=s2
END FUNCTION
PURE FUNCTION s2a_3(s1,s2,s3) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN)             :: s1, s2, s3
    CHARACTER(LEN=1000), DIMENSION(3)        :: a

  a(1)=s1; a(2)=s2; a(3)=s3
END FUNCTION
PURE FUNCTION s2a_4(s1,s2,s3,s4) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN)             :: s1, s2, s3, s4
    CHARACTER(LEN=1000), DIMENSION(4)        :: a

  a(1)=s1; a(2)=s2; a(3)=s3; a(4)=s4
END FUNCTION
PURE FUNCTION s2a_5(s1,s2,s3,s4,s5) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN)             :: s1, s2, s3, s4, s5
    CHARACTER(LEN=1000), DIMENSION(5)        :: a

  a(1)=s1; a(2)=s2; a(3)=s3; a(4)=s4; a(5)=s5
END FUNCTION
PURE FUNCTION s2a_6(s1,s2,s3,s4,s5,s6) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN)             :: s1, s2, s3, s4, s5, s6
    CHARACTER(LEN=1000), DIMENSION(6)        :: a

  a(1)=s1; a(2)=s2; a(3)=s3; a(4)=s4; a(5)=s5; a(6)=s6
END FUNCTION
PURE FUNCTION s2a_7(s1,s2,s3,s4,s5,s6,s7) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN)             :: s1, s2, s3, s4, s5, s6, s7
    CHARACTER(LEN=1000), DIMENSION(7)        :: a

  a(1)=s1; a(2)=s2; a(3)=s3; a(4)=s4; a(5)=s5; a(6)=s6; a(7)=s7
END FUNCTION
PURE FUNCTION s2a_8(s1,s2,s3,s4,s5,s6,s7,s8) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN)             :: s1, s2, s3, s4, s5, s6, s7, s8
    CHARACTER(LEN=1000), DIMENSION(8)        :: a

  a(1)=s1; a(2)=s2; a(3)=s3; a(4)=s4; a(5)=s5; a(6)=s6; a(7)=s7
  a(8)=s8
END FUNCTION
PURE FUNCTION s2a_9(s1,s2,s3,s4,s5,s6,s7,s8,s9) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN)             :: s1, s2, s3, s4, s5, s6, s7, &
                                                s8, s9
    CHARACTER(LEN=1000), DIMENSION(9)        :: a

  a(1)=s1; a(2)=s2; a(3)=s3; a(4)=s4; a(5)=s5; a(6)=s6; a(7)=s7
  a(8)=s8; a(9)=s9
END FUNCTION
PURE FUNCTION s2a_10(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN)             :: s1, s2, s3, s4, s5, s6, s7, &
                                                s8, s9, s10
    CHARACTER(LEN=1000), DIMENSION(10)       :: a

  a(1)=s1; a(2)=s2; a(3)=s3; a(4)=s4; a(5)=s5; a(6)=s6; a(7)=s7
  a(8)=s8; a(9)=s9; a(10)=s10
END FUNCTION
PURE FUNCTION s2a_11(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN)             :: s1, s2, s3, s4, s5, s6, s7, &
                                                s8, s9, s10, s11
    CHARACTER(LEN=1000), DIMENSION(11)       :: a

  a(1)=s1; a(2)=s2; a(3)=s3; a(4)=s4; a(5)=s5; a(6)=s6; a(7)=s7
  a(8)=s8; a(9)=s9; a(10)=s10; a(11)=s11
END FUNCTION
PURE FUNCTION s2a_12(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN)             :: s1, s2, s3, s4, s5, s6, s7, &
                                                s8, s9, s10, s11, s12
    CHARACTER(LEN=1000), DIMENSION(12)       :: a

  a(1)=s1; a(2)=s2; a(3)=s3; a(4)=s4; a(5)=s5; a(6)=s6; a(7)=s7
  a(8)=s8; a(9)=s9; a(10)=s10; a(11)=s11; a(12)=s12
END FUNCTION
PURE FUNCTION s2a_13(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN)             :: s1, s2, s3, s4, s5, s6, s7, &
                                                s8, s9, s10, s11, s12, s13
    CHARACTER(LEN=1000), DIMENSION(13)       :: a

  a(1)=s1; a(2)=s2; a(3)=s3; a(4)=s4; a(5)=s5; a(6)=s6; a(7)=s7
  a(8)=s8; a(9)=s9; a(10)=s10; a(11)=s11; a(12)=s12; a(13)=s13
END FUNCTION
PURE FUNCTION s2a_14(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14) RESULT(a)
    CHARACTER(LEN=*), INTENT(IN)             :: s1, s2, s3, s4, s5, s6, s7, &
                                                s8, s9, s10, s11, s12, s13, &
                                                s14
    CHARACTER(LEN=1000), DIMENSION(14)       :: a

  a(1)=s1; a(2)=s2; a(3)=s3; a(4)=s4; a(5)=s5; a(6)=s6; a(7)=s7
  a(8)=s8; a(9)=s9; a(10)=s10; a(11)=s11; a(12)=s12; a(13)=s13; a(14)=s14
END FUNCTION
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

  SUBROUTINE substitute_special_xml_tokens(inp_string,out_string,ltu)

!   Purpose: Substitute special XML tokens like "<" or ">" in inp_string.
!            Optionally convert also all lowercase characters to uppercase, if
!            ltu is true.

!   History: - Creation (10.03.2005,MK)

!   ***************************************************************************

    CHARACTER(LEN=*), INTENT(IN)             :: inp_string
    CHARACTER(LEN=*), INTENT(OUT)            :: out_string
    LOGICAL, INTENT(IN), OPTIONAL            :: ltu

    CHARACTER(LEN=LEN(inp_string))           :: string
    INTEGER                                  :: i, j

!   ---------------------------------------------------------------------------

    string = inp_string
    out_string = ""

    IF (PRESENT(ltu)) THEN
      IF (ltu) CALL uppercase(string)
    END IF

    j = 1
    DO i=1,LEN_TRIM(string)
      IF (string(i:i) == "<") THEN
        out_string(j:j+3) = "&lt;"
        j = j + 4
      ELSE IF (string(i:i) == ">") THEN
        out_string(j:j+3) = "&gt;"
        j = j + 4
      ELSE IF (string(i:i) == "&") THEN
        out_string(j:j+4) = "&amp;"
        j = j + 5
      ELSE IF (string(i:i) == """") THEN
        out_string(j:j+5) = "&quot;"
        j = j + 6
      ELSE
        out_string(j:j) = string(i:i)
        j = j + 1
      END IF
    END DO

  END SUBROUTINE substitute_special_xml_tokens

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
