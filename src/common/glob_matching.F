!
!  Imported from Arjen Markus' flibs project
!  http://flibs.sourceforge.net/
!
! globmatch.f90 --
!     Match strings according to (simplified) glob patterns
!
!     The pattern matching is limited to literals, * and ?
!     (character classes are not supported). A backslash escapes
!     any character.
!
! *****************************************************************************
MODULE glob_matching
  USE f77_blas

    IMPLICIT NONE

    CHARACTER(len=1), PARAMETER, PRIVATE :: backslash = '\\'
    CHARACTER(len=1), PARAMETER, PRIVATE :: star      = '*'
    CHARACTER(len=1), PARAMETER, PRIVATE :: question  = '?'

  PRIVATE
  PUBLIC :: string_match

CONTAINS

! string_match --
!     Tries to match the given string with the pattern
! Arguments:
!     string     String to be examined
!     pattern    Glob pattern to be used for the matching
! Result:
!     .true. if the entire string matches the pattern, .false.
!     otherwise
! Note:
!     Trailing blanks are ignored
!
! *****************************************************************************
!> \brief ...
!> \param string ...
!> \param pattern ...
!> \retval match ...
! *****************************************************************************
RECURSIVE FUNCTION string_match( string, pattern ) RESULT(match)
    CHARACTER(len=*), INTENT(in)             :: string, pattern
    LOGICAL                                  :: match

    CHARACTER(len=LEN(pattern))              :: literal
    INTEGER                                  :: k, ll, method, p, ptrim, &
                                                start, strim

    match  = .FALSE.
    method = 0
    ptrim  = LEN_TRIM( pattern )
    strim  = LEN_TRIM( string )
    p      = 1
    ll     = 0
    start  = 1

    !
    ! Split off a piece of the pattern
    !
    DO WHILE ( p <= ptrim )
        SELECT CASE ( pattern(p:p) )
            CASE( star )
                IF ( ll .NE. 0 ) EXIT
                method = 1
            CASE( question )
                IF ( ll .NE. 0 ) EXIT
                method = 2
                start  = start + 1
            CASE( backslash )
                p  = p + 1
                ll = ll + 1
                literal(ll:ll) = pattern(p:p)
            CASE default
                ll = ll + 1
                literal(ll:ll) = pattern(p:p)
        END SELECT

        p = p + 1
    ENDDO

    !
    ! Now look for the literal string (if any!)
    !
    IF ( method == 0 ) THEN
        !
        ! We are at the end of the pattern, and of the string?
        !
        IF ( strim == 0 .AND. ptrim == 0 ) THEN
            match = .TRUE.
        ELSE
            !
            ! The string matches a literal part?
            !
            IF ( ll > 0 ) THEN
                IF ( string(start:MIN(strim,start+ll-1)) == literal(1:ll) ) THEN
                    start = start + ll
                    match = string_match( string(start:), pattern(p:) )
                ENDIF
            ENDIF
        ENDIF
    ENDIF

    IF ( method == 1 ) THEN
        !
        ! Scan the whole of the remaining string ...
        !
        IF ( ll == 0 ) THEN
            match = .TRUE.
        ELSE
            DO WHILE ( start <= strim )
                k     = INDEX( string(start:), literal(1:ll) )
                IF ( k > 0 ) THEN
                    start = start + k + ll - 1
                    match = string_match( string(start:), pattern(p:) )
                    IF ( match ) THEN
                        EXIT
                    ENDIF
                ENDIF

                start = start + 1
            ENDDO
        ENDIF
    ENDIF

    IF ( method == 2 .AND. ll > 0 ) THEN
        !
        ! Scan the whole of the remaining string ...
        !
        IF ( string(start:MIN(strim,start+ll-1)) == literal(1:ll) ) THEN
            match = string_match( string(start+ll:), pattern(p:) )
        ENDIF
    ENDIF
    RETURN
END FUNCTION string_match

END MODULE glob_matching
