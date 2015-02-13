! xmlparse.f90 - Simple, limited XML parser in Fortran
!
! $Id: xmlparse.f90,v 1.14 2007/12/07 10:10:19 arjenmarkus Exp $
!
! Arjen Markus
!
! General information:
! The module reads XML files by:
! - Identifying the tag and all attributes and data belonging
!   to the tag.
! - Returning to the calling subprogram to let it take care of
!   the tag, attributes and data.
! - If the tag is actually an ending tag, then this is flagged
!   too.
! - Handling all the data is left to the calling subprogram,
!   the module merely facilitates in the parsing.
!
! Note:
! The module in its current version has a number of limitations:
! - It does not handle escape sequences (like &gt. to signify
!   a ">" sign)
! - It does not handle tags with attributes that are spread
!   over more than one line
! - The maximum length of a line is 1000 characters
! - It may report too many lines of data (empty lines)
! - No DOM support nor support for an object tree
! - It is probably not very robust in detecting malformed XML files
!
! Some questions:
! - What to do with leading blanks?
!
! Update - several ideas:
! - Introduce at least two options (via xml_options):
!   - ignore_whitespace  - remove leading blanks and leading and trailing
!                          empty lines from the PCDATA
!   - no_data_truncation - consider truncation of data (more
!                          attributes or lines of character data than
!                          can be stored) a read error
! - Introduce convenience functions and subroutines:
!   - xml_ok()           - all is well, reading can continue
!   - xml_data_trunc()   - was there truncation of the data?
!   - xml_find_attrib()  - find an attribute by name
!
! Further ideas:
!   - simple checking via a table: parent, tag, id, min, max
!
MODULE xml_parser

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: XML_PARSE,xml_close,xml_get,xml_ok,xml_open

   INTEGER, PARAMETER :: XML_BUFFER_LENGTH = 1000

   ! Define the data type that holds the parser information

   TYPE XML_PARSE
      INTEGER          :: lun                ! LU-number of the XML-file
      INTEGER          :: level              ! Indentation level (output)
      INTEGER          :: lineno             ! Line in file
      LOGICAL          :: ignore_whitespace  ! Ignore leading blanks etc.
      LOGICAL          :: no_data_truncation ! Do not allow data truncation
      LOGICAL          :: too_many_attribs   ! More attributes than could be stored?
      LOGICAL          :: too_many_data      ! More lines of data than could be stored?
      LOGICAL          :: eof                ! End of file?
      LOGICAL          :: error              ! Invalid XML file or other error?
      CHARACTER(len=XML_BUFFER_LENGTH) :: line  ! Buffer
   END TYPE XML_PARSE

   ! Global options

   INTEGER, PARAMETER    :: XML_STDOUT       = -1
   INTEGER, PRIVATE      :: report_lun_      = XML_STDOUT
   LOGICAL, PRIVATE      :: report_errors_   = .FALSE.
   LOGICAL, PRIVATE      :: report_details_  = .FALSE.

   ! Global data (the ampersand must come first)

   CHARACTER(len=10), DIMENSION(2,3), SAVE, PRIVATE :: entities = &
      RESHAPE( (/ '&    ', '&amp;', &
                  '>    ', '&gt; ',  &
                  '<    ', '&lt; ' /), (/2,3/) )

   ! Auxiliary routines - private

   PRIVATE               :: xml_compress_
   PRIVATE               :: xml_put_open_tag_
   PRIVATE               :: xml_put_element_
   PRIVATE               :: xml_put_close_tag_
   PRIVATE               :: xml_replace_entities_
   !
   ! Interfaces to reporting routines
   !
   PRIVATE               :: xml_report_details_int_
   PRIVATE               :: xml_report_details_string_
   PRIVATE               :: xml_report_errors_int_
   PRIVATE               :: xml_report_errors_string_

   INTERFACE xml_report_details
      MODULE PROCEDURE xml_report_details_int_
      MODULE PROCEDURE xml_report_details_string_
   END INTERFACE
   INTERFACE xml_report_errors
      MODULE PROCEDURE xml_report_errors_int_
      MODULE PROCEDURE xml_report_errors_string_
      MODULE PROCEDURE xml_report_errors_extern_
   END INTERFACE

CONTAINS

! xml_report_details_int_ --
!    Routine to write a text with an integer value
! Arguments:
!    text        Text to be written
!    int         Integer value to be added
!
! *****************************************************************************
!> \brief ...
!> \param text ...
!> \param int ...
! *****************************************************************************
SUBROUTINE xml_report_details_int_( text, int )
    CHARACTER(len=*), INTENT(in)             :: text
    INTEGER, INTENT(in)                      :: int

   IF ( report_details_ ) THEN
      IF ( report_lun_ .EQ. XML_STDOUT ) THEN
         WRITE(*,*) TRIM(text), int
      ELSE
         WRITE(report_lun_,*) TRIM(text), int
      ENDIF
   ENDIF
END SUBROUTINE xml_report_details_int_

! xml_report_details_string_ --
!    Routine to write a text with a string value
! Arguments:
!    text        Text to be written
!    string      String to be added
!
! *****************************************************************************
!> \brief ...
!> \param text ...
!> \param string ...
! *****************************************************************************
SUBROUTINE xml_report_details_string_( text, string )
    CHARACTER(len=*), INTENT(in)             :: text, string

   IF ( report_details_ ) THEN
      IF ( report_lun_ .EQ. XML_STDOUT ) THEN
         WRITE(*,*) TRIM(text), ' ', TRIM(string)
      ELSE
         WRITE(report_lun_,*) TRIM(text), ' ', TRIM(string)
      ENDIF
   ENDIF
END SUBROUTINE xml_report_details_string_


! xml_report_errors_string_ --
!    Routine to write an error message text with an integer value
! Arguments:
!    text        Text to be written
!    int         Integer value to be added
!    lineno      Line number in the file
!
! *****************************************************************************
!> \brief ...
!> \param text ...
!> \param int ...
!> \param lineno ...
! *****************************************************************************
SUBROUTINE xml_report_errors_int_( text, int, lineno )
    CHARACTER(len=*), INTENT(in)             :: text
    INTEGER, INTENT(in)                      :: int
    INTEGER, INTENT(in), OPTIONAL            :: lineno

   IF ( report_errors_ .OR. report_details_ ) THEN
      IF ( report_lun_ .EQ. XML_STDOUT ) THEN
         WRITE(*,*) TRIM(text), int
         IF ( PRESENT(lineno) ) THEN
            WRITE(*,*) '   At or near line', lineno
         ENDIF
      ELSE
         WRITE(report_lun_,*) TRIM(text), int
         IF ( PRESENT(lineno) ) THEN
            WRITE(report_lun_,*) '   At or near line', lineno
         ENDIF
      ENDIF
   ENDIF
END SUBROUTINE xml_report_errors_int_

! xml_report_errors_string_ --
!    Routine to write an error message text with a string value
! Arguments:
!    text        Text to be written
!    string      String to be added
!    lineno      Line number in the file
!
! *****************************************************************************
!> \brief ...
!> \param text ...
!> \param string ...
!> \param lineno ...
! *****************************************************************************
SUBROUTINE xml_report_errors_string_( text, string, lineno )
    CHARACTER(len=*), INTENT(in)             :: text, string
    INTEGER, INTENT(in), OPTIONAL            :: lineno

   IF ( report_errors_ .OR. report_details_ ) THEN
      IF ( report_lun_ .EQ. XML_STDOUT ) THEN
         WRITE(*,*) TRIM(text), ' ', TRIM(string)
         IF ( PRESENT(lineno) ) THEN
            WRITE(*,*) '   At or near line', lineno
         ENDIF
      ELSE
         WRITE(report_lun_,*) TRIM(text), ' ', TRIM(string)
         IF ( PRESENT(lineno) ) THEN
            WRITE(report_lun_,*) '   At or near line', lineno
         ENDIF
      ENDIF
   ENDIF
END SUBROUTINE xml_report_errors_string_

! xml_report_errors_extern_ --
!    Routine to write an error message text with a string value
! Arguments:
!    info        Structure holding information on the XML-file
!    text        Text to be written
! Note:
!    This routine is meant for use by routines outside
!    this module
!
! *****************************************************************************
!> \brief ...
!> \param info ...
!> \param text ...
! *****************************************************************************
SUBROUTINE xml_report_errors_extern_( info, text )
    TYPE(XML_PARSE), INTENT(in)              :: info
    CHARACTER(len=*), INTENT(in)             :: text

   IF ( report_lun_ .EQ. XML_STDOUT ) THEN
      WRITE(*,*) TRIM(text), ' - at or near line', info%lineno
   ELSE
      WRITE(report_lun_,*) TRIM(text), ' - at or near line', info%lineno
   ENDIF
END SUBROUTINE xml_report_errors_extern_

! xml_open --
!    Routine to open an XML file for reading or writing
! Arguments:
!    info        Structure holding information on the XML-file
!    fname       Name of the file
!    mustread    The file will be read (.true.) or written (.false.)
!
! *****************************************************************************
!> \brief ...
!> \param info ...
!> \param fname ...
!> \param mustread ...
! *****************************************************************************
SUBROUTINE xml_open( info, fname, mustread )
    TYPE(XML_PARSE), INTENT(out)             :: info
    CHARACTER(len=*), INTENT(in)             :: fname
    LOGICAL, INTENT(in)                      :: mustread

    INTEGER                                  :: i, ierr, k, kend
    LOGICAL                                  :: exists, opend

   info%lun = 10
   info%ignore_whitespace  = .FALSE.
   info%no_data_truncation = .FALSE.
   info%too_many_attribs   = .FALSE.
   info%too_many_data      = .FALSE.
   info%eof                = .FALSE.
   info%error              = .FALSE.
   info%level              = -1
   info%lineno             =  0

   DO i = 10,99
      INQUIRE( unit = i, opened = opend )
      IF ( .NOT. opend ) THEN
         info%lun = i
         INQUIRE( file = fname, exist = exists )
         IF ( .NOT. exists .AND. mustread ) THEN
            CALL xml_report_errors( 'XML_OPEN: file does not exist:', TRIM(fname))
            info%lun   = -1
            info%error = .TRUE.
         ELSE
            OPEN( unit = info%lun, file = fname )
            CALL xml_report_details( 'XML_OPEN: opened file ', TRIM(fname) )
            CALL xml_report_details( 'at LU-number: ', info%lun )
         ENDIF
         EXIT
      ENDIF
   ENDDO
   IF ( .NOT. info%error .AND. mustread ) THEN
      k = 1
      DO WHILE ( k .GE. 1 )
         READ( info%lun, '(a)', iostat = ierr ) info%line
         IF ( ierr .EQ. 0 ) THEN
            info%line = ADJUSTL(  info%line )
            k         = INDEX( info%line, '<?' )
            !
            ! Assume (for now at least) that <?xml ... ?> appears on a single line!
            !
            IF ( k .GE. 1 ) THEN
               kend = INDEX( info%line, '?>' )
               IF ( kend .LE. 0 ) THEN
                  CALL xml_report_errors( 'XML_OPEN: error reading file with LU-number: ', info%lun )
                  CALL xml_report_errors( 'Line starting with "<?xml" should end with "?>"', ' ' )
                  info%error = .TRUE.
                  EXIT
               ENDIF
            ENDIF
         ELSE
            CALL xml_report_errors( 'XML_OPEN: error reading file with LU-number: ', info%lun )
            CALL xml_report_errors( 'Possibly no line starting with "<?xml"', ' ' )
            CALL xml_close( info )
            info%error = .TRUE.
            EXIT
         ENDIF
      ENDDO
   ENDIF
   IF ( .NOT. info%error .AND. .NOT. mustread ) THEN
       WRITE( info%lun, '(a)' ) '<?xml version="1.0"?>'
   ENDIF
END SUBROUTINE xml_open

! xml_close --
!    Routine to close an XML file
! Arguments:
!    info        Structure holding information on the XML-file
!
! *****************************************************************************
!> \brief ...
!> \param info ...
! *****************************************************************************
SUBROUTINE xml_close( info )
    TYPE(XML_PARSE), INTENT(inout)           :: info

   CLOSE( info%lun )

   !
   ! Only clean up the LU-number, so that the calling program
   ! can examine the last condition
   !
   CALL xml_report_details( 'XML_CLOSE: Closing file with LU-number ', info%lun )
   info%lun              = -1
END SUBROUTINE xml_close

! xml_get --
!    Routine to get the next bit of information from an XML file
! Arguments:
!    info        Structure holding information on the XML-file
!    tag         Tag that was encountered
!    endtag      Whether the end of the element was encountered
!    attribs     List of attribute-value pairs
!    no_attribs  Number of pairs in the list
!    data        Lines of character data found
!    no_data     Number of lines of character data
!
! *****************************************************************************
!> \brief ...
!> \param info ...
!> \param tag ...
!> \param endtag ...
!> \param attribs ...
!> \param no_attribs ...
!> \param DATA ...
!> \param no_data ...
! *****************************************************************************
SUBROUTINE xml_get( info, tag, endtag, attribs, no_attribs, &
                   DATA, no_data )
    TYPE(XML_PARSE), INTENT(inout)           :: info
    CHARACTER(len=*), INTENT(out)            :: tag
    LOGICAL, INTENT(out)                     :: endtag
    CHARACTER(len=*), DIMENSION(:, :), &
      INTENT(out)                            :: attribs
    INTEGER, INTENT(out)                     :: no_attribs
    CHARACTER(len=*), DIMENSION(:), &
      INTENT(out)                            :: DATA
    INTEGER, INTENT(out)                     :: no_data

    CHARACTER(len=XML_BUFFER_LENGTH)         :: nextline
    INTEGER                                  :: idxat, idxdat, ierr, kend, &
                                                keq, kfirst, ksecond, kspace
    LOGICAL                                  :: close_bracket, comment_tag

!
! Initialise the output
!

   endtag     = .FALSE.
   no_attribs = 0
   no_data    = 0

   info%too_many_attribs = .FALSE.
   info%too_many_data    = .FALSE.

   IF ( info%lun .LT. 0 ) THEN
      CALL xml_report_details( 'XML_GET on closed file ', ' ' )
      RETURN
   ENDIF

   !
   ! From the previous call or the call to xmlopen we have
   ! the line that we need to parse already in memory:
   ! <tag attrib1="..." attrib2="..." />
   !
   comment_tag   = .FALSE.
   close_bracket = .FALSE.
   kspace        = INDEX( info%line, ' ' )
   kend          = INDEX( info%line, '>' )
   DO WHILE ( kend .LE. 0 )
      READ( info%lun, '(a)', iostat = ierr ) nextline
      info%lineno = info%lineno + 1

      IF ( ierr .EQ. 0 ) THEN
         info%line = TRIM(info%line) // ' ' // ADJUSTL(nextline)
      ELSE
         info%error = .TRUE.
         CALL xml_report_errors( 'XML_GET - end of tag not found ', &
                 '(buffer too small?)', info%lineno )
         CALL xml_close( info )
         RETURN
      ENDIF
      kend = INDEX( info%line, '>' )
   ENDDO
   IF ( kend .GT. kspace ) THEN
      kend = kspace
   ELSE
      close_bracket = .TRUE.
   ENDIF

   !
   ! Check for the end of an ordianry tag and of
   ! a comment tag
   !
   IF ( info%line(1:3) .EQ. '-->' ) THEN
      endtag = .TRUE.
      tag    = info%line(4:kend-1)
   ELSE IF ( info%line(1:2) .EQ. '</' ) THEN
      endtag = .TRUE.
      tag    = info%line(3:kend-1)
   ELSE
      IF ( info%line(1:1) .EQ. '<' ) THEN
         tag    = info%line(2:kend-1)
         CALL xml_report_details( 'XML_GET - tag found: ', TRIM(tag) )
      ELSE
         kend   = 0 ! Beginning of data!
      ENDIF
   ENDIF

   info%line = ADJUSTL( info%line(kend+1:) )

   idxat     = 0
   idxdat    = 0

   IF ( tag(1:3) .EQ. '!--' ) comment_tag = .TRUE.

   DO WHILE ( info%line .NE. ' ' .AND. .NOT. close_bracket .AND. .NOT. comment_tag )

      keq  = INDEX( info%line, '=' )
      kend = INDEX( info%line, '>' )
      IF ( keq .GT. kend ) keq = 0 ! Guard against multiple tags
                                   ! with attributes on one line

      !
      ! No attributes any more?
      !
      IF ( keq .LT. 1 ) THEN
         kend = INDEX( info%line, '/>' )
         IF ( kend .GE. 1 ) THEN
            kend   = kend + 1 ! To go beyond the ">" character
            endtag = .TRUE.
         ELSE
            kend = INDEX( info%line, '>' )
            IF ( kend .LT. 1 ) THEN
               CALL xml_report_errors( 'XML_GET - wrong ending of tag ', &
                       TRIM(info%line), info%lineno  )
               info%error = .TRUE. ! Wrong ending of line!
               CALL xml_close( info )
               RETURN
            ELSE
               close_bracket = .TRUE.
            ENDIF
         ENDIF
         IF ( kend .GE. 1 ) THEN
            info%line = ADJUSTL( info%line(kend+1:) )
         ENDIF
         EXIT
      ENDIF

      idxat = idxat + 1
      IF ( idxat .LE. SIZE(attribs,2) ) THEN
         no_attribs = idxat
         attribs(1,idxat) = ADJUSTL(info%line(1:keq-1)) ! Use adjustl() to avoid
                                                        ! multiple spaces, etc
         info%line = ADJUSTL( info%line(keq+1:) )

         !
         ! We have almost found the start of the attribute's value
         !
         kfirst  = INDEX( info%line, '"' )
         IF ( kfirst .LT. 1 ) THEN
            CALL xml_report_errors( 'XML_GET - malformed attribute-value pair: ', &
                    TRIM(info%line), info%lineno  )
            info%error = .TRUE. ! Wrong form of attribute-value pair
            CALL xml_close( info )
            RETURN
         ENDIF

         ksecond = INDEX( info%line(kfirst+1:), '"' ) + kfirst
         IF ( ksecond .LT. 1 ) THEN
            CALL xml_report_errors( 'XML_GET - malformed attribute-value pair: ', &
                    TRIM(info%line), info%lineno  )
            info%error = .TRUE. ! Wrong form of attribute-value pair
            CALL xml_close( info )
            RETURN
         ENDIF

         attribs(2,idxat) = info%line(kfirst+1:ksecond-1)
         info%line = ADJUSTL( info%line(ksecond+1:) )
      ENDIF

      IF ( idxat .GT. SIZE(attribs,2) ) THEN
         CALL xml_report_errors( 'XML_GET - more attributes than could be stored: ', &
                 TRIM(info%line), info%lineno  )
         info%too_many_attribs = .TRUE.
         info%line             = ' '
         EXIT
      ENDIF
   ENDDO

   !
   ! Now read the data associated with the current tag
   ! - all the way to the next "<" character
   !
   ! To do: reduce the number of data lines - empty ones
   ! at the end should not count.
   !
   DO
      IF ( comment_tag ) THEN
         kend   = INDEX( info%line, '-->' )
      ELSE
         kend   = INDEX( info%line, '<' )
      ENDIF
      idxdat = idxdat + 1
      IF ( idxdat .LE. SIZE(DATA) ) THEN
         no_data = idxdat
         IF ( kend .GE. 1 ) THEN
            DATA(idxdat) = info%line(1:kend-1)
            info%line    = info%line(kend:)
         ELSE
            DATA(idxdat) = info%line
         ENDIF
      ELSE
         CALL xml_report_errors( 'XML_GET - more data lines than could be stored: ', &
                 TRIM(info%line), info%lineno  )
         info%too_many_data = .TRUE.
         EXIT
      ENDIF

      !
      ! No more data? Otherwise, read on
      !
      IF ( kend .GE. 1 ) THEN
         EXIT
      ELSE
         READ( info%lun, '(a)', iostat = ierr ) info%line
         info%lineno = info%lineno + 1

         IF ( ierr .LT. 0 ) THEN
            CALL xml_report_details( 'XML_GET - end of file found - LU-number: ', &
                    info%lun )
            info%eof = .TRUE.
         ELSEIF ( ierr .GT. 0 ) THEN
            CALL xml_report_errors( 'XML_GET - error reading file with LU-number ', &
                    info%lun, info%lineno  )
            info%error = .TRUE.
         ENDIF
         IF ( ierr .NE. 0 ) THEN
            EXIT
         ENDIF
      ENDIF
   ENDDO

   !
   ! Compress the data?
   !
   IF ( info%ignore_whitespace ) THEN
      CALL xml_compress_( DATA, no_data )
   ENDIF

   !
   ! Replace the entities, if any
   !
   CALL xml_replace_entities_( DATA, no_data )

   CALL xml_report_details( 'XML_GET - number of attributes: ', no_attribs )
   CALL xml_report_details( 'XML_GET - number of data lines: ', no_data    )

END SUBROUTINE xml_get

! xml_put --
!    Routine to write a tag with the associated data to an XML file
! Arguments:
!    info        Structure holding information on the XML-file
!    tag         Tag that was encountered
!    endtag      Whether the end of the element was encountered
!    attribs     List of attribute-value pairs
!    no_attribs  Number of pairs in the list
!    data        Lines of character data found
!    no_data     Number of lines of character data
!    type        Type of action:
!                open - just the opening tag with attributes
!                elem - complete element
!                close - just the closing tag
!
! *****************************************************************************
!> \brief ...
!> \param info ...
!> \param tag ...
!> \param attribs ...
!> \param no_attribs ...
!> \param DATA ...
!> \param no_data ...
!> \param TYPE ...
! *****************************************************************************
SUBROUTINE xml_put(info, tag, attribs, no_attribs, &
                   DATA, no_data, TYPE)

    TYPE(XML_PARSE), INTENT(inout)           :: info
    CHARACTER(len=*), INTENT(in)             :: tag
    CHARACTER(len=*), DIMENSION(:, :), &
      INTENT(in)                             :: attribs
    INTEGER, INTENT(in)                      :: no_attribs
    CHARACTER(len=*), DIMENSION(:), &
      INTENT(in)                             :: DATA
    INTEGER, INTENT(in)                      :: no_data
    CHARACTER(len=*)                         :: TYPE

   SELECT CASE(TYPE)
      CASE('open')
         CALL xml_put_open_tag_(info, tag, attribs, no_attribs, &
                  DATA, no_data)
      CASE('elem')
         CALL xml_put_element_(info, tag, attribs, no_attribs, &
                  DATA, no_data)
      CASE('close')
         CALL xml_put_close_tag_(info, tag, attribs, no_attribs, &
                  DATA, no_data)
   END SELECT

END SUBROUTINE xml_put

! xml_put_open_tag_ --
!    Routine to write the opening tag with the attributes
! Arguments:
!    info        Structure holding information on the XML-file
!    tag         Tag that was encountered
!    endtag      Whether the end of the element was encountered
!    attribs     List of attribute-value pairs
!    no_attribs  Number of pairs in the list
!    data        Lines of character data found
!    no_data     Number of lines of character data
!
! *****************************************************************************
!> \brief ...
!> \param info ...
!> \param tag ...
!> \param attribs ...
!> \param no_attribs ...
!> \param DATA ...
!> \param no_data ...
! *****************************************************************************
SUBROUTINE xml_put_open_tag_(info, tag, attribs, no_attribs, &
                             DATA, no_data)

    TYPE(XML_PARSE), INTENT(inout)           :: info
    CHARACTER(len=*), INTENT(in)             :: tag
    CHARACTER(len=*), DIMENSION(:, :), &
      INTENT(in)                             :: attribs
    INTEGER, INTENT(in)                      :: no_attribs
    CHARACTER(len=*), DIMENSION(:), &
      INTENT(in)                             :: DATA
    INTEGER, INTENT(in)                      :: no_data

    CHARACTER(len=300), PARAMETER            :: indent = ' '

    INTEGER                                  :: i

   WRITE( info%lun, '(3a)', advance = 'no' ) &
         indent(1:3*info%level), '<', ADJUSTL(tag)
   DO i=1,no_attribs
      IF (attribs(2,i).NE.'') THEN
         WRITE( info%lun, '(5a)', advance = 'no' ) &
            ' ',TRIM(attribs(1,i)),'="', TRIM(attribs(2,i)),'"'
      ENDIF
   ENDDO
   WRITE( info%lun, '(a)' ) '>'
   info%level = info%level + 1

END SUBROUTINE xml_put_open_tag_

! xml_put_element_ --
!    Routine to write the complete element
! Arguments:
!    info        Structure holding information on the XML-file
!    tag         Tag that was encountered
!    endtag      Whether the end of the element was encountered
!    attribs     List of attribute-value pairs
!    no_attribs  Number of pairs in the list
!    data        Lines of character data found
!    no_data     Number of lines of character data
!
! *****************************************************************************
!> \brief ...
!> \param info ...
!> \param tag ...
!> \param attribs ...
!> \param no_attribs ...
!> \param DATA ...
!> \param no_data ...
! *****************************************************************************
SUBROUTINE xml_put_element_(info, tag, attribs, no_attribs, &
                            DATA, no_data)

    TYPE(XML_PARSE), INTENT(inout)           :: info
    CHARACTER(len=*), INTENT(in)             :: tag
    CHARACTER(len=*), DIMENSION(:, :), &
      INTENT(in)                             :: attribs
    INTEGER, INTENT(in)                      :: no_attribs
    CHARACTER(len=*), DIMENSION(:), &
      INTENT(in)                             :: DATA
    INTEGER, INTENT(in)                      :: no_data

    CHARACTER(len=300), PARAMETER            :: indent = ' '

    CHARACTER(len=1)                         :: aa
    INTEGER                                  :: i, ii
    LOGICAL                                  :: logic

   IF ( (no_attribs.eq.0 .AND. no_data.eq.0) ) THEN
      RETURN
   ELSE
      logic = .TRUE.
      DO ii = 1,no_attribs
         logic = logic .AND. (attribs(2,ii).EQ.'')
      ENDDO
      DO ii = 1,no_data
         logic = logic .AND. (DATA(ii).EQ.'')
      ENDDO
      IF ( logic ) THEN
         RETURN
      ELSE
         WRITE( info%lun, '(3a)', advance = 'no' ) &
             indent(1:3*info%level), '<', ADJUSTL(tag)
         DO i = 1,no_attribs
            IF (attribs(2,i).NE.'') THEN
               WRITE( info%lun, '(5a)', advance = 'no' ) &
                   ' ',TRIM(attribs(1,i)),'="', TRIM(attribs(2,i)),'"'
            ENDIF
         ENDDO
         IF ( no_attribs.gt.0 .AND. no_data.eq.0 ) THEN
            aa='a'
         ELSEIF ( (no_attribs.gt.0 .AND. no_data.gt.0) .OR. &
                  (no_attribs.eq.0 .AND. no_data.gt.0) ) THEN
            aa='b'
         ELSE
            WRITE(*,*) no_attribs, no_data
         ENDIF
      ENDIF
   ENDIF

   SELECT CASE(aa)
      CASE('a')
         WRITE( info%lun, '(a)' ) '/>'
      CASE('b')
         WRITE( info%lun, '(a)',advance='no' ) '>'
         WRITE( info%lun, '(2a)', advance='no') &
               ( ' ', TRIM(DATA(i)), i=1,no_data )
         WRITE( info%lun, '(4a)' ) ' ','</', tag, '>'
    END SELECT

END SUBROUTINE xml_put_element_

! xml_put_close_tag_ --
!    Routine to write the closing tag
! Arguments:
!    info        Structure holding information on the XML-file
!    tag         Tag that was encountered
!    endtag      Whether the end of the element was encountered
!    attribs     List of attribute-value pairs
!    no_attribs  Number of pairs in the list
!    data        Lines of character data found
!    no_data     Number of lines of character data
!
! *****************************************************************************
!> \brief ...
!> \param info ...
!> \param tag ...
!> \param attribs ...
!> \param no_attribs ...
!> \param DATA ...
!> \param no_data ...
! *****************************************************************************
SUBROUTINE xml_put_close_tag_(info, tag, attribs, no_attribs, &
                     DATA, no_data)

    TYPE(XML_PARSE), INTENT(inout)           :: info
    CHARACTER(len=*), INTENT(in)             :: tag
    CHARACTER(len=*), DIMENSION(:, :), &
      INTENT(in)                             :: attribs
    INTEGER, INTENT(in)                      :: no_attribs
    CHARACTER(len=*), DIMENSION(:), &
      INTENT(in)                             :: DATA
    INTEGER, INTENT(in)                      :: no_data

    CHARACTER(len=300), PARAMETER            :: indent = ' '

   info%level=info%level-1
   WRITE( info%lun, '(4a)' ) &
        indent(1:3*info%level), '</', ADJUSTL(tag), '>'

END SUBROUTINE xml_put_close_tag_

! xml_compress_ --
!    Routine to remove empty lines from the character data
! Arguments:
!    data        Lines of character data found
!    no_data     (Nett) number of lines of character data
!
! *****************************************************************************
!> \brief ...
!> \param DATA ...
!> \param no_data ...
! *****************************************************************************
SUBROUTINE xml_compress_( DATA, no_data )
    CHARACTER(len=*), DIMENSION(:), &
      INTENT(inout)                          :: DATA
    INTEGER, INTENT(inout)                   :: no_data

    INTEGER                                  :: i, j
    LOGICAL                                  :: empty

   j     = 0
   empty = .TRUE.
   DO i = 1,no_data
      IF ( LEN_TRIM(DATA(i)) .NE. 0 .OR. .NOT. empty ) THEN
         j       = j + 1
         DATA(j) = ADJUSTL(DATA(i))
         empty = .FALSE.
      ENDIF
   ENDDO

   no_data = j

   DO i = no_data,1,-1
      IF ( LEN_TRIM(DATA(i)) .NE. 0 ) THEN
         EXIT
      ELSE
         no_data = no_data - 1
      ENDIF
   ENDDO

END SUBROUTINE xml_compress_

! xml_replace_entities_ --
!    Routine to replace entities such as &gt; by their
!    proper character representation
! Arguments:
!    data        Lines of character data found
!    no_data     (Nett) number of lines of character data
!
! *****************************************************************************
!> \brief ...
!> \param DATA ...
!> \param no_data ...
! *****************************************************************************
SUBROUTINE xml_replace_entities_( DATA, no_data )
    CHARACTER(len=*), DIMENSION(:), &
      INTENT(inout)                          :: DATA
    INTEGER, INTENT(inout)                   :: no_data

    INTEGER                                  :: i, j, j2, k, pos
    LOGICAL                                  :: found

   DO i = 1,no_data
      j = 1
      DO
         DO k = 1,SIZE(entities,2)
            found = .FALSE.
            pos   = INDEX( DATA(i)(j:), TRIM(entities(2,k)) )
            IF ( pos .GT. 0 ) THEN
               found = .TRUE.
               j     = j + pos - 1
               j2    = j + LEN_TRIM(entities(2,k))
               DATA(i)(j:) = TRIM(entities(1,k)) // DATA(i)(j2:)
               j     = j2
            ENDIF
         ENDDO
         IF ( .NOT. found ) EXIT
      ENDDO
   ENDDO

END SUBROUTINE xml_replace_entities_

! xml_options --
!    Routine to handle the parser options
! Arguments:
!    info                Structure holding information on the XML-file
!    ignore_whitespace   Ignore whitespace (leading blanks, empty lines) or not
!    no_data_truncation  Consider truncation of strings an error or not
!    report_lun          LU-number for reporting information
!    report_errors       Write messages about errors or not
!    report_details      Write messages about all kinds of actions or not
!
! *****************************************************************************
!> \brief ...
!> \param info ...
!> \param ignore_whitespace ...
!> \param no_data_truncation ...
!> \param report_lun ...
!> \param report_errors ...
!> \param report_details ...
! *****************************************************************************
SUBROUTINE xml_options( info, ignore_whitespace, no_data_truncation, &
                        report_lun, report_errors, &
                        report_details )
    TYPE(XML_PARSE), INTENT(inout)           :: info
    LOGICAL, INTENT(in), OPTIONAL            :: ignore_whitespace, &
                                                no_data_truncation
    INTEGER, INTENT(in), OPTIONAL            :: report_lun
    LOGICAL, INTENT(in), OPTIONAL            :: report_errors, report_details

   IF ( PRESENT(ignore_whitespace) ) THEN
      info%ignore_whitespace = ignore_whitespace
   ENDIF
   IF ( PRESENT(no_data_truncation) ) THEN
      info%no_data_truncation = no_data_truncation
   ENDIF
   IF ( PRESENT(report_lun) ) THEN
      report_lun_ = report_lun
   ENDIF
   IF ( PRESENT(report_errors) ) THEN
      report_errors_ = report_errors
   ENDIF
   IF ( PRESENT(report_details) ) THEN
      report_details_ = report_details
   ENDIF
END SUBROUTINE xml_options

! xml_ok --
!    Function that returns whether all was okay or not
! Arguments:
!    info                Structure holding information on the XML-file
! Returns:
!    .true. if there was no error, .false. otherwise
!
! *****************************************************************************
!> \brief ...
!> \param info ...
!> \retval xml_ok ...
! *****************************************************************************
LOGICAL FUNCTION xml_ok( info )
   TYPE(XML_PARSE),  INTENT(in)               :: info

   xml_ok = info%eof .OR. info%error .OR. &
            ( info%no_data_truncation .AND.    &
                 ( info%too_many_attribs .OR. info%too_many_data ) )
   xml_ok = .NOT. xml_ok
END FUNCTION xml_ok

! xml_error --
!    Function that returns whether there was an error
! Arguments:
!    info                Structure holding information on the XML-file
! Returns:
!    .true. if there was an error, .false. if there was none
!
! *****************************************************************************
!> \brief ...
!> \param info ...
!> \retval xml_error ...
! *****************************************************************************
LOGICAL FUNCTION xml_error( info )
   TYPE(XML_PARSE),  INTENT(in)               :: info

   xml_error = info%error .OR. &
            ( info%no_data_truncation .AND.    &
                 ( info%too_many_attribs .OR. info%too_many_data ) )
END FUNCTION xml_error

! xml_data_trunc --
!    Function that returns whether data were truncated or not
! Arguments:
!    info                Structure holding information on the XML-file
! Returns:
!    .true. if data were truncated, .false. otherwise
!
! *****************************************************************************
!> \brief ...
!> \param info ...
!> \retval xml_data_trunc ...
! *****************************************************************************
LOGICAL FUNCTION xml_data_trunc( info )
   TYPE(XML_PARSE),  INTENT(in)               :: info

   xml_data_trunc = info%too_many_attribs .OR. info%too_many_data
END FUNCTION xml_data_trunc

! *****************************************************************************
!> \brief ...
!> \param attribs ...
!> \param no_attribs ...
!> \param name ...
!> \param value ...
!> \retval xml_find_attrib ...
! *****************************************************************************
INTEGER FUNCTION xml_find_attrib( attribs, no_attribs, name, value )
   CHARACTER(len=*), DIMENSION(:,:)  :: attribs
   INTEGER                           :: no_attribs
   CHARACTER(len=*)                  :: name
   CHARACTER(len=*)                  :: value

   INTEGER :: i

   xml_find_attrib = -1
   DO i = 1,no_attribs
      IF ( name .EQ. attribs(1,i) ) THEN
         value           = attribs(2,i)
         xml_find_attrib = i
         EXIT
      ENDIF
   ENDDO

END FUNCTION xml_find_attrib

! xml_process --
!    Routine to read the XML file as a whole and distribute processing
!    the contents over three user-defined subroutines
! Arguments:
!    filename            Name of the file to process
!    attribs             Array for holding the attributes
!    data                Array for holding the character data
!    startfunc           Subroutine to handle the start of elements
!    datafunc            Subroutine to handle the character data
!    endfunc             Subroutine to handle the end of elements
!    error               Indicates if there was an error or not
! Note:
!    The routine is declared recursive to allow inclusion of XML files
!    (common with XSD schemas). This extends to the auxiliary routines.
!
RECURSIVE &
! *****************************************************************************
!> \brief ...
!> \param filename ...
!> \param attribs ...
!> \param DATA ...
!> \param startfunc ...
!> \param datafunc ...
!> \param endfunc ...
!> \param lunrep ...
!> \param error ...
! *****************************************************************************
SUBROUTINE xml_process( filename, attribs, DATA, startfunc, datafunc, endfunc, lunrep, error )
    CHARACTER(len=*)                         :: filename
    CHARACTER(len=*), DIMENSION(:, :)        :: attribs
    CHARACTER(len=*), DIMENSION(:)           :: DATA
   INTERFACE
      RECURSIVE SUBROUTINE startfunc( tag, attribs, error )
         CHARACTER(len=*)                  :: tag
         CHARACTER(len=*), DIMENSION(:,:)  :: attribs
         LOGICAL                           :: error
      END SUBROUTINE
   END INTERFACE
   INTERFACE
      RECURSIVE SUBROUTINE datafunc( tag, DATA, error )
         CHARACTER(len=*)                  :: tag
         CHARACTER(len=*), DIMENSION(:)    :: DATA
         LOGICAL                           :: error
      END SUBROUTINE
   END INTERFACE
   INTERFACE
      RECURSIVE SUBROUTINE endfunc( tag, error )
         CHARACTER(len=*)                  :: tag
         LOGICAL                           :: error
      END SUBROUTINE
   END INTERFACE
    INTEGER                                  :: lunrep
    LOGICAL                                  :: error

    CHARACTER(len=80)                        :: tag
    INTEGER                                  :: noattribs, nodata
    LOGICAL                                  :: endtag
    TYPE(XML_PARSE)                          :: info

   CALL xml_options( info, report_lun = lunrep, report_details = .FALSE. )
   CALL xml_open( info, filename, .TRUE. )

   error = .FALSE.
   DO
      CALL xml_get( info, tag, endtag, attribs, noattribs, DATA, nodata )
      IF ( .NOT. xml_ok(info) ) THEN
         EXIT
      ENDIF

      IF ( xml_error(info) ) THEN
         WRITE(lunrep,*) 'Error reading XML file!'
         error = .TRUE.
         EXIT
      ENDIF

      IF ( .NOT. endtag .OR. noattribs .NE. 0 ) THEN
         CALL startfunc( tag, attribs(:,1:noattribs), error )
         IF ( error ) EXIT

         CALL datafunc( tag, DATA(1:nodata), error )
         IF ( error ) EXIT
      ENDIF

      IF ( endtag ) THEN
         CALL endfunc( tag, error )
         IF ( error ) EXIT
      ENDIF
   ENDDO
   CALL xml_close( info )
END SUBROUTINE xml_process

END MODULE xml_parser
