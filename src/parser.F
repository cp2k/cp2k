!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****** cp2k/parser [1.0] *
!!
!!   NAME
!!     parser
!!
!!   FUNCTION
!!     Provides all the functions needed to parse an input section
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     JGH (28.11.2001) : New overloaded call get_next, replacing interface
!!                        for cfield, get_int, get_real
!!                      : delete get_array, it was never used
!!
!!   NOTES
!!----------------------------------------------------------------------------
!!          F R E E - F O R M A T  I N P U T  R O U T I N E S
!!
!! Subroutine  READ_LINE : reads an input line, physical end of lines
!!                         can be masked by a & character, no further
!!                         characters are alowed to the right
!! Function    TEST_NEXT() : tests contains of next input field, return values
!!                           are : E (blank); C (text); N (number);
!!                                 K (bracket); S (! comment); X (&END);
!!                                 O all other
!! Subroutine  GET_NEXT(...) : overloaded function of CFIELD, GET_INT, GET_REAL
!!
!! Subroutine  CFIELD(string,length) : extracts the next field from the
!!                                     input line, at most length charcters
!!                                     are returned in string, if length is 0
!!                                     the full field is returned, possible
!!                                     field seperators are: blank, comma,
!!                                     colon, semicolon, bracket, equal sign
!! Subroutine  GET_INT(int) : read an integer from the input line
!! Subroutine  GET_REAL(real) : read an real value from the input line
!!
!! Subroutine  P_ERROR : writes out current input line and current position
!!                       within the line, does not stop code
!!----------------------------------------------------------------------------
!!
!!   SOURCE
!-----------------------------------------------------------------------------

MODULE parser

  USE kinds, ONLY : dbl

  USE global_types, ONLY : global_environment_type
  USE mathconstants, ONLY : zero, one
  USE string_utilities, ONLY : uppercase, compress
  USE termination, ONLY : stop_program
  USE util, ONLY : get_unit
  USE message_passing, ONLY : mp_bcast, mp_sync

  PRIVATE
  PUBLIC :: parser_init, parser_end, read_line, test_next, get_next, &
       p_error, stop_parser, search_label

  INTEGER, PARAMETER :: max_line = 20
  CHARACTER ( LEN = 132 ) :: inp_line(max_line)
  INTEGER :: num_inp_line
  INTEGER :: p_line, p_char
  INTEGER :: iunit
  CHARACTER ( LEN = 1 ), PARAMETER :: newline = '&'
  CHARACTER ( LEN = 1 ), PARAMETER :: blank = ' '

  TYPE(global_environment_type) :: parenv
  LOGICAL :: parallel_io, ionode
  CHARACTER ( LEN = 1 ), PARAMETER :: comment_default='!'
  CHARACTER ( LEN = 1 ) :: comment_sign = comment_default
  CHARACTER ( LEN = 4 ), PARAMETER :: end_section_default='&END'
  CHARACTER ( LEN = 10 ) :: end_section = end_section_default
  INTEGER :: end_section_length

  INTERFACE get_next
    MODULE PROCEDURE get_real, get_int, cfield
  END INTERFACE

  INTERFACE read_line
    MODULE PROCEDURE read_line_s, read_line_p
  END INTERFACE

!******************************************************************************

CONTAINS

!******************************************************************************

SUBROUTINE parser_init ( file_name, globenv, section_end, comment )
  IMPLICIT NONE
  CHARACTER ( LEN = * ), INTENT ( IN ) :: file_name
  TYPE(global_environment_type), INTENT ( IN ), OPTIONAL :: globenv
  CHARACTER ( LEN = * ), INTENT ( IN ), OPTIONAL :: section_end
  CHARACTER ( LEN = 1 ), INTENT ( IN ), OPTIONAL :: comment
  CHARACTER ( LEN = 132 ) :: line
  INTEGER :: ios
  LOGICAL :: exists

!------------------------------------------------------------------------------

! setup parameters for parallel parser
  IF ( PRESENT(globenv) ) THEN
    ionode = globenv%ionode
    parallel_io = .true.
    parenv = globenv
  ELSE
    ionode = .true.
    parallel_io = .false.
  ENDIF

  IF ( PRESENT(section_end) ) THEN
    end_section = section_end
    end_section_length = LEN_TRIM ( section_end )
  ELSE
    end_section = end_section_default
    end_section_length = 4
  ENDIF

  IF ( PRESENT(comment) ) THEN
    comment_sign = comment
  ELSE
    comment_sign = comment_default
  ENDIF



  IF(ionode) THEN

!..find a not associated unit
    iunit = get_unit()
    line = blank

!..check the file name
    INQUIRE ( FILE = file_name, EXIST = exists )
    IF ( .NOT. exists ) THEN
       CALL stop_parser ( 'parser_init', 'file does not exist' )
    ELSE
       OPEN (unit=iunit,file=file_name)
    END IF
  END IF

END SUBROUTINE parser_init

!******************************************************************************

SUBROUTINE parser_end

  IMPLICIT NONE

  IF(ionode) CLOSE ( UNIT = iunit )

END SUBROUTINE parser_end

!******************************************************************************

SUBROUTINE search_label ( label, ierr, ignore_case )

  IMPLICIT NONE

  INTEGER, INTENT ( OUT ) :: ierr
  CHARACTER ( LEN = * ), INTENT ( INOUT ) :: label
  LOGICAL, INTENT ( IN ), OPTIONAL :: ignore_case

  CHARACTER ( LEN = 132 ) :: line
  INTEGER :: ios
  LOGICAL :: exists, igc

!------------------------------------------------------------------------------

  ierr = 0

  IF(ionode) THEN

    igc = .FALSE.
    IF ( PRESENT ( ignore_case ) ) igc = ignore_case

    IF ( igc ) CALL uppercase ( label )

    REWIND ( iunit )

    DO
       READ ( iunit, iostat = ios, fmt = '( a )' ) line
       IF ( ios > 0 ) THEN
          CALL stop_parser ( 'search_label','read error')
       ELSE IF ( ios < 0 ) THEN
          ierr = 2
          EXIT
       ELSE
          IF ( igc ) CALL uppercase ( line )
          IF ( index ( line, label ) /= 0 ) EXIT
       END IF
    END DO
  END IF

  IF ( parallel_io ) call mp_bcast ( ierr, parenv%source, parenv%group )
  
END SUBROUTINE search_label

!******************************************************************************
SUBROUTINE read_line_s ( line )
  IMPLICIT NONE
  CHARACTER ( LEN = * ) :: line
  INTEGER :: ios

!------------------------------------------------------------------------------

  IF ( ionode ) THEN

    READ (iunit,iostat=ios,fmt='(a)' ) line
    IF ( ios /= 0 ) CALL stop_parser ( 'read_line','error while reading file')

  END IF

!..distribute the line to all processors (for parallel parser)
  IF(parallel_io) THEN
    CALL mp_bcast(line,parenv%source,parenv%group)
  END IF

END SUBROUTINE read_line_s

!******************************************************************************

RECURSIVE SUBROUTINE read_line_p
  IMPLICIT NONE
  INTEGER :: llen, i, j, ii, ios
  CHARACTER ( LEN = 132 ) :: line
  CHARACTER ( LEN = 1 ) :: last
  LOGICAL :: comment

!------------------------------------------------------------------------------

  IF(ionode) THEN

!..skip all empty lines
    llen = 0
    DO WHILE (llen==0)
       READ (iunit,iostat=ios,fmt='(a )' ) line
       IF ( ios /= 0 ) CALL stop_parser ( 'read_line','error while reading file')
       llen = MAX(len_trim(line),0) ! again the xlf problem (see qs_parser.F)
    END DO

!..now read the input line
    last = line(llen:llen)
    num_inp_line = 1
    inp_line(num_inp_line) = line

    DO WHILE (last==newline)
       IF (num_inp_line>=max_line) CALL stop_parser ( 'read_line','max_line')
       READ (iunit,iostat=ios,fmt='(a )' ) line
       IF ( ios /= 0 ) CALL stop_parser ( 'read_line','error while reading file')
       llen = len_trim(line)
       last = line(llen:llen)
       num_inp_line = num_inp_line + 1
       inp_line(num_inp_line) = line
    END DO

!..delete newline characters
    DO i = 1, num_inp_line - 1
       llen = len_trim(inp_line(i))
       inp_line(i) (llen:llen) = blank
    END DO

!..convert all delimiter characters to blanks
    DO i = 1, num_inp_line
       llen = len_trim(inp_line(i))
       DO j = 1, llen
          SELECT CASE (inp_line(i)(j:j))
          CASE DEFAULT
          CASE ( '=')
             inp_line(i) (j:j) = blank
          CASE ( ',')
             inp_line(i) (j:j) = blank
          CASE ( ':')
             inp_line(i) (j:j) = blank
          CASE ( ';')
             inp_line(i) (j:j) = blank
          END SELECT
       END DO
    END DO

!..delete comment sections
    DO i = 1, num_inp_line
       line = blank
       llen = len_trim(inp_line(i))
       comment = .FALSE.
       DO j = 1, llen
          IF (inp_line(i)(j:j)==comment_sign .AND. .NOT. comment) THEN
             comment = .TRUE.
             inp_line(i)(j:j)= blank
          ELSE IF (inp_line(i)(j:j)==comment_sign .AND. comment) THEN
             comment = .FALSE.
             inp_line(i)(j:j)= blank
          ELSE IF (comment) THEN
             inp_line(i)(j:j)= blank
          END IF
       END DO
       CALL compress ( inp_line ( i ) )
    END DO

!..delete empty lines
    j = 0
    DO i = 1, num_inp_line
       llen = len_trim(inp_line(i))
       IF (llen>0) THEN
          j = j + 1
          IF (j/=i) inp_line(j) = inp_line(i)
       END IF
    END DO

    num_inp_line = j

!..add line seperators
    DO i = 1, num_inp_line - 1
       llen = len_trim(inp_line(i))
       inp_line(i) (llen+2:llen+2) = '$'
    END DO
    IF (num_inp_line>0) THEN
       llen = len_trim(inp_line(num_inp_line))
       inp_line(num_inp_line) (llen+2:llen+2) = '%'
    END IF

  END IF

!..distribute the input line to all processors (for parallel parser)
  IF(parallel_io) THEN
    CALL mp_bcast(num_inp_line,parenv%source,parenv%group)
    DO i = 1, num_inp_line
      CALL mp_bcast(inp_line(i),parenv%source,parenv%group)
    END DO
  END IF

!..initialize pointers
  p_line = 1
  p_char = 1

  IF ( num_inp_line == 0 ) CALL read_line_p

END SUBROUTINE read_line_p

!******************************************************************************

FUNCTION test_next()

  IMPLICIT NONE

! Locals
  CHARACTER ( LEN = 1 ) :: test_next
  INTEGER :: n

!------------------------------------------------------------------------------

  SELECT CASE (inp_line(p_line)(p_char:p_char))
  CASE (blank)
     test_next = 'E'
  CASE ( 'A':'Z')
     test_next = 'C'
  CASE ( 'a':'z')
     test_next = 'C'
  CASE ( '0':'9')
     test_next = 'N'
  CASE ( '+')
     test_next = 'N'
  CASE ( '-')
     test_next = 'N'
  CASE ( '.')
     test_next = 'N'
  CASE ( '[')
     test_next = 'K'
  CASE ( '!')
     test_next = 'S'
  CASE ( '%')
     test_next = 'P'
  CASE DEFAULT
     test_next = 'O'
  END SELECT
  n = end_section_length
  IF (inp_line(p_line)(p_char:p_char+n-1)==end_section(1:n)) test_next = 'X'
END FUNCTION test_next

!******************************************************************************

SUBROUTINE cfield ( string, length )

  IMPLICIT NONE

! Arguments
  CHARACTER ( LEN = * ), INTENT ( INOUT ) :: string
  INTEGER, INTENT ( INOUT ) :: length

! Locals
  CHARACTER ( LEN = 1 ) :: char1
  INTEGER :: l, ll, i, ii

!------------------------------------------------------------------------------

  string = blank
  IF (length==0) length = 999
  l = len(inp_line(p_line)(p_char:))
  ll = length
  IF (l<length) ll = l
  DO i = 1, ll
     char1 = inp_line(p_line) (p_char+i-1:p_char+i-1)
     IF (char1==blank) EXIT
     IF (char1==',') EXIT
     IF (char1==';') EXIT
     IF (char1==':') EXIT
     IF (char1=='=') EXIT
     IF (char1=='[') EXIT
     IF (char1=='!') EXIT
     string(i:i) = char1
  END DO
  l = len_trim(inp_line(p_line))
  DO i = p_char, l
     IF (inp_line(p_line)(i:i)==blank .OR. inp_line(p_line)(i:i)==',' &
          .OR. inp_line(p_line)(i:i)==':' .OR. inp_line(p_line)(i:i)==';' &
          .OR. inp_line(p_line)(i:i)=='=') THEN
        ii = i
        p_char = ii + 1
        EXIT
     ELSE IF (inp_line(p_line)(i:i)=='[' .OR. inp_line(p_line)(i:i)=='$' &
          .OR. inp_line(p_line)(i:i)=='%') THEN
        ii = i
        p_char = ii
        EXIT
     END IF
     ii = i
  END DO

  IF (inp_line(p_line)(p_char:p_char)=='$' .AND. p_char==l) THEN
     p_line = p_line + 1
     p_char = 1
  END IF

  IF (length==999) length = len_trim(string)

END SUBROUTINE cfield

!******************************************************************************

SUBROUTINE p_error ( ioutunit )

  IMPLICIT NONE

! Arguments
  INTEGER, INTENT ( IN ), OPTIONAL :: ioutunit

! Locals
  INTEGER :: i, llen, iout

!------------------------------------------------------------------------------

  IF ( PRESENT ( ioutunit ) ) THEN
     iout = ioutunit
  ELSEIF ( parallel_io ) THEN
     iout = parenv%scr
  ELSE
     iout = 6
  END IF
  IF ( parallel_io ) call mp_sync(parenv%group)
  IF ( ionode ) then
    WRITE (iout,'(a )' ) ' Error while parsing input file '
    WRITE (iout,'(a )' ) ' The last line read was : '
    DO i = 1, num_inp_line
       llen = len_trim(inp_line(i))
       WRITE (iout,'(1x,a )' ) inp_line(i) (1:llen)
    END DO
    WRITE (iout,'(a )' ) ' The current position in the line is:'
    llen = len_trim(inp_line(p_line))
    WRITE (iout,'(1x,a )' ) inp_line(p_line) (p_char:llen)
  END IF
  IF ( parallel_io ) call mp_sync(parenv%group)

END SUBROUTINE p_error

!******************************************************************************

SUBROUTINE get_int(int)
  IMPLICIT NONE

! Return value
  INTEGER, INTENT(OUT) :: int

! Localse
  REAL ( dbl ) :: rival

!------------------------------------------------------------------------------

  CALL get_real(rival)
  int = NINT ( rival )
  IF ( ABS ( int - rival ) > EPSILON ( rival ) ) THEN
     CALL p_error()
     CALL stop_parser ( 'get_int','not an integer field')
  END IF

END SUBROUTINE get_int

!******************************************************************************

SUBROUTINE get_real(reval)

  IMPLICIT NONE

!..most general format :  I((sign)_._E(sign)_/(sign)_._E(sign)_)

! Return value
  REAL ( dbl ), INTENT ( OUT ) :: reval

! Locals
  REAL ( dbl ) :: xx, v1, v2
  REAL ( dbl ), SAVE :: value
  LOGICAL, SAVE :: repetition = .FALSE.
  INTEGER, SAVE :: rep_num=0
  CHARACTER ( LEN = 256 ) :: string
  INTEGER :: ie, llen, state, ii
  REAL ( dbl ) :: s1, s2, e1, e2, se1, se2, x1, x2, y1, y2

!------------------------------------------------------------------------------

  IF (repetition) THEN
     reval = value
     rep_num = rep_num - 1
     IF (rep_num==0) repetition = .FALSE.
  ELSE
     s1 = one
     s2 = one
     e1 = zero
     e2 = zero
     se1 = one
     se2 = one
     x1 = zero
     x2 = one
     y1 = zero
     y2 = zero
     IF (test_next()/='N') THEN
        CALL p_error
        CALL stop_parser ( 'get_real','not an numeric field')
     END IF
     llen = 0
     CALL cfield(string,llen)
     ie = scan(string,'(')
     IF (ie /= 0 ) THEN
        xx = zero
        CALL get_num(string(1:ie),xx,llen)
        rep_num = nint(xx)
        IF (rep_num<=0) THEN
           CALL p_error
           CALL stop_parser ( 'get_real','non valid repetition number')
        END IF
        IF (rep_num>1) repetition = .TRUE.
        rep_num = rep_num - 1
        string(1:ie) = blank
        ie = scan(string,' )' )
        string(ie:) = blank
        string = adjustl(string)
     END IF
     state = 0
     llen = len_trim(string)
     ii = 1
     DO WHILE (ii<=llen)
        SELECT CASE (string(ii:ii))
        CASE ( '-')
           IF (state==0) THEN
              s1 = -one
              state = 1
              ii = ii + 1
           ELSE IF (state==5) THEN
              se1 = -one
              state = 6
              ii = ii + 1
           ELSE IF (state==8) THEN
              s2 = -one
              state = 9
              ii = ii + 1
           ELSE IF (state==13) THEN
              se2 = -one
              state = 14
              ii = ii + 1
           ELSE
              CALL p_error ()
              CALL stop_parser ( 'get_real','not a valid number')
           END IF
        CASE ( '+')
           IF (state==0) THEN
              state = 1
              ii = ii + 1
           ELSE IF (state==5) THEN
              state = 6
              ii = ii + 1
           ELSE IF (state==8) THEN
              state = 9
              ii = ii + 1
           ELSE IF (state==13) THEN
              state = 14
              ii = ii + 1
           ELSE
              CALL p_error ()
              CALL stop_parser ( 'get_real','not a valid number')
           END IF
        CASE ( '.')
           IF (state<=2) THEN
              state = 3
              ii = ii + 1
           ELSE IF (state>7 .AND. state<=10) THEN
              state = 11
              ii = ii + 1
           ELSE
              CALL p_error ()
              CALL stop_parser ( 'get_real','not a valid number')
           END IF
        CASE ( 'E', 'e', 'D', 'd' )
           IF (state>=2 .AND. state<=4) THEN
              state = 5
              ii = ii + 1
           ELSE IF (state>=10 .AND. state<=12) THEN
              state = 13
              ii = ii + 1
           ELSE
              CALL p_error ()
              CALL stop_parser ( 'get_real','not a valid number')
           END IF
        CASE ( '/')
           x2 = zero
           IF (state>=2 .AND. state<=4) THEN
              state = 8
              ii = ii + 1
           ELSE IF (state==7) THEN
              state = 8
              ii = ii + 1
           ELSE
              CALL p_error ()
              CALL stop_parser ( 'get_real','not a valid number')
           END IF
        CASE ( '0':'9')
           IF (state==0 .OR. state==1) THEN
              CALL get_num(string(ii:),x1,ie)
              state = 2
              ii = ii + ie
           ELSE IF (state==3) THEN
              CALL get_num(string(ii:),y1,ie)
              y1 = y1*10._dbl**(-ie)
              state = 4
              ii = ii + ie
           ELSE IF (state==5 .OR. state==6) THEN
              CALL get_num(string(ii:),e1,ie)
              state = 7
              ii = ii + ie
           ELSE IF (state==8 .OR. state==9) THEN
              CALL get_num(string(ii:),x2,ie)
              state = 10
              ii = ii + ie
           ELSE IF (state==11) THEN
              CALL get_num(string(ii:),y2,ie)
              y2 = y2*10._dbl**(-ie)
              state = 12
              ii = ii + ie
           ELSE IF (state==13 .OR. state==14) THEN
              CALL get_num(string(ii:),e2,ie)
              state = 15
              ii = llen + 1
           ELSE
              CALL p_error ()
              CALL stop_parser ( 'get_real','not a valid number')
           END IF
        CASE DEFAULT
           ii = ii + 1
        END SELECT
     END DO

     v1 = s1*(x1+y1)*10._dbl**(se1*e1)
     v2 = s2*(x2+y2)*10._dbl**(se2*e2)
     IF (abs(v2)<epsilon(v2)) THEN
        CALL p_error ()
        CALL stop_parser ( 'get_real','divisor too small')
     END IF
     value = v1/v2
  END IF
  reval = value

END SUBROUTINE get_real

!******************************************************************************

SUBROUTINE get_num ( string, result, rlen )

  IMPLICIT NONE

! Arguments
  REAL ( dbl ), INTENT ( OUT ) :: result
  CHARACTER ( LEN = * ), INTENT ( IN ) :: string
  INTEGER, INTENT ( OUT ) :: rlen

! Locals
  REAL ( dbl ) :: sign, new
  INTEGER :: istart, llen, i

  REAL ( dbl ), PARAMETER :: ten = 10.0_dbl

!------------------------------------------------------------------------------

  result = zero
  llen = len_trim(string)
  IF (llen==0) THEN
     CALL p_error ()
     CALL stop_parser ( 'get_num','no number')
  END IF
  IF (string(1:1)=='-') THEN
     sign = -one
     istart = 2
     rlen = 1
  ELSE IF (string(1:1)=='+') THEN
     sign = one
     istart = 2
     rlen = 1
  ELSE
     sign = one
     istart = 1
     rlen = 0
  END IF

  DO i = istart, llen
     SELECT CASE (string(i:i))
     CASE ( '0')
        new = 0._dbl
        rlen = rlen + 1
     CASE ( '1')
        new = 1._dbl
        rlen = rlen + 1
     CASE ( '2')
        new = 2._dbl
        rlen = rlen + 1
     CASE ( '3')
        new = 3._dbl
        rlen = rlen + 1
     CASE ( '4')
        new = 4._dbl
        rlen = rlen + 1
     CASE ( '5')
        new = 5._dbl
        rlen = rlen + 1
     CASE ( '6')
        new = 6._dbl
        rlen = rlen + 1
     CASE ( '7')
        new = 7._dbl
        rlen = rlen + 1
     CASE ( '8')
        new = 8._dbl
        rlen = rlen + 1
     CASE ( '9')
        new = 9._dbl
        rlen = rlen + 1
     CASE DEFAULT
        EXIT
     END SELECT
     result = ten*result + new
  END DO

  result = sign * result

END SUBROUTINE get_num

!******************************************************************************

SUBROUTINE afield ( str )

! Arguments
  CHARACTER ( LEN = * ), INTENT ( OUT ) :: str

! Locals
  INTEGER :: i, l, ii

!------------------------------------------------------------------------------

  str = blank

  l = len(inp_line(p_line)(p_char:))
  DO i = 1, l
     str(i:i) = inp_line(p_line) (p_char+i-1:p_char+i-1)
     IF (str(i:i)==']') EXIT
  END DO

  l = len_trim(str)
  IF (str(l:l)/=']') THEN
     CALL p_error ()
     CALL stop_parser ( 'afield','not an array field')
  END IF
  p_char = p_char + l

  l = len_trim(inp_line(p_line))
  DO i = p_char, l
     IF (inp_line(p_line)(i:i)==blank .OR. inp_line(p_line)(i:i)==',' &
          .OR. inp_line(p_line)(i:i)==':' .OR. inp_line(p_line)(i:i)==';' &
          .OR. inp_line(p_line)(i:i)=='=') THEN
        ii = i
        p_char = ii + 1
        EXIT
     ELSE IF (inp_line(p_line)(i:i)=='[' .OR. inp_line(p_line)(i:i)=='$' &
          .OR. inp_line(p_line)(i:i)=='%') THEN
        ii = i
        p_char = ii
        EXIT
     END IF
     ii = i
  END DO

  IF (inp_line(p_line)(p_char:p_char)=='$' .AND. p_char==l) THEN
     p_line = p_line + 1
     p_char = 1
  END IF

END SUBROUTINE afield

!******************************************************************************

SUBROUTINE stop_parser ( routine, commentary )

  IMPLICIT NONE

! Arguments
  CHARACTER ( LEN = * ), INTENT ( IN ) :: routine, commentary

  IF ( parallel_io) THEN
    IF ( ionode ) CALL stop_program( routine, commentary )
    CALL mp_sync ( parenv % group )
  ELSE
    CALL stop_program( routine, commentary )
  END IF

END SUBROUTINE stop_parser

!******************************************************************************

END MODULE parser

!******************************************************************************
!!*****
