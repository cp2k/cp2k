!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief provides a uniform framework to add references to CP2K
!>      cite and output these
!> \note
!>      references need to be input using the ISI citation format, because it is
!>      uniform, easy to parse, and can be exported for example from web of science
!>      furthermore, it can be easily converted to and from using the bibutils tools
!>      a collection of easy to use conversion programs that can be found at
!>      http://www.scripps.edu/~cdputnam/software/bibutils/
!>      by Chris Putnam
!>
!>      see thebibliography.F on how to add references easily
!> \par History
!>      08.2007 [Joost VandeVondele]
!> \author Joost VandeVondele
! **************************************************************************************************
MODULE reference_manager
   USE kinds,                           ONLY: default_string_length
   USE message_passing,                 ONLY: mp_para_env_type
   USE string_utilities,                ONLY: substitute_special_xml_tokens,&
                                              uppercase
   USE util,                            ONLY: sort
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PUBLIC :: cite_reference
   PUBLIC :: collect_citations_from_ranks
   PUBLIC :: print_cited_references
   PUBLIC :: export_references_as_xml

   PUBLIC :: add_reference          ! use this one only in bibliography.F
   PUBLIC :: remove_all_references  ! use only in f77_interface.F
   PUBLIC :: get_citation_key       ! a string key describing the reference (e.g. Kohn1965b)

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'reference_manager'

   ! maximum number of reference that can be added
   INTEGER, PARAMETER :: max_reference = 1024

   ! storage of a reference
   INTEGER, PARAMETER :: doi_length = 128
   INTEGER, PARAMETER :: ISI_length = 128

   ! the way we store a reference, should remain fully private
   TYPE reference_type
      PRIVATE
      ! the reference in a format as returned by the web of science
      CHARACTER(LEN=ISI_length), DIMENSION(:), POINTER :: ISI_record => NULL()
      ! the doi only, i.e. without "https://doi.org/"
      CHARACTER(LEN=doi_length)                        :: DOI = ""
      ! has this reference been cited in the program run
      LOGICAL                                          :: is_cited = .FALSE.
      ! this is a citation key for output in the reference lists
      CHARACTER(LEN=ISI_length)                        :: citation_key = ""
   END TYPE reference_type

   ! useful to build arrays
   TYPE reference_p_type
      TYPE(reference_type), POINTER :: ref => NULL()
   END TYPE

   ! the bibliography
   INTEGER, SAVE :: nbib = 0
   TYPE(reference_p_type), DIMENSION(max_reference) :: thebib

CONTAINS

! **************************************************************************************************
!> \brief marks a given reference as cited.
!> \param key citation key as returned from add_reference
!> \par History
!>      XX.2007 created [ ]
! **************************************************************************************************
   SUBROUTINE cite_reference(key)
      INTEGER, INTENT(IN)                                :: key

      IF (key < 1 .OR. key > max_reference) CPABORT("citation key out of range")

      ! set as cited
      thebib(key)%ref%is_cited = .TRUE.

   END SUBROUTINE

! **************************************************************************************************
!> \brief Checks for each reference if any mpi-rank has marked it for citation.
!> \param para_env ...
!> \par History
!>      12.2013 created [Ole Schuett]
! **************************************************************************************************
   SUBROUTINE collect_citations_from_ranks(para_env)
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env

      INTEGER                                            :: i, t

      DO i = 1, nbib
         t = 0
         IF (thebib(i)%ref%is_cited) t = 1
         CALL para_env%max(t)
         thebib(i)%ref%is_cited = (t == 1)
      END DO

   END SUBROUTINE collect_citations_from_ranks

! **************************************************************************************************
!> \brief add a reference to the bibliography
!> \param key output, this handle is needed to cite this reference later
!> \param ISI_record ...
!> \param DOI ...
!> \par History
!>      08.2007 created [Joost VandeVondele]
!> \note
!>      - see bibliography.F for it use.
!>      - the ISI record is space sensitive, in particular the first three characters need to be blank
!>        or contain a key indicating the record type. See the header of this file for tools
!>        that can convert e.g. bibtex or endnote files to the ISI format
!>      - DOI: provide the DOI without a link. The link will be automatically created as needed.
! **************************************************************************************************
   SUBROUTINE add_reference(key, ISI_record, DOI)
      INTEGER, INTENT(OUT)                               :: key
      CHARACTER(LEN=*), DIMENSION(:), INTENT(IN)         :: ISI_record
      CHARACTER(LEN=*), INTENT(IN)                       :: DOI

      CHARACTER                                          :: tmp
      CHARACTER(LEN=ISI_length)                          :: author, citation_key, key_a, key_b, year
      INTEGER                                            :: commaloc, i, ires, line, match, mylen, &
                                                            nlines

      IF (nbib + 1 > max_reference) CPABORT("increase max_reference")
      nbib = nbib + 1
      key = nbib

      ! initialize reference to zero
      ALLOCATE (thebib(key)%ref)
      NULLIFY (thebib(key)%ref%ISI_record)
      thebib(key)%ref%DOI = ""
      thebib(key)%ref%is_cited = .FALSE.

      ! Assign DOI
      thebib(key)%ref%DOI = DOI

      ! Assign ISI_record
      nlines = SIZE(ISI_record, 1)
      ALLOCATE (thebib(key)%ref%ISI_record(nlines))
      thebib(key)%ref%ISI_record = ISI_record

      ! construct a citation_key
      line = 1
      author = get_next_author(thebib(key)%ref%ISI_record, line)
      commaloc = INDEX(author, ',')
      IF (commaloc .GT. 0) author = author(1:commaloc - 1)
      CPASSERT(LEN_TRIM(author) > 0)
      year = get_year(thebib(key)%ref%ISI_record)
      CPASSERT(LEN_TRIM(year) == 4)
      citation_key = TRIM(author)//TRIM(year)

      ! avoid special characters in names, just remove them
      mylen = LEN_TRIM(citation_key)
      ires = 0
      DO I = 1, mylen
       IF (INDEX("0123456789thequickbrownfoxjumpsoverthelazydogTHEQUICKBROWNFOXJUMPSOVERTHELAZYDOG", citation_key(i:i)) .NE. 0) THEN
            ires = ires + 1
            tmp = citation_key(i:i)
            citation_key(ires:ires) = tmp
         END IF
      END DO
      citation_key(ires + 1:) = ""
      CPASSERT(LEN_TRIM(citation_key) > 4) ! At least one character of the author should be left.

      ! avoid duplicates, search through the list for matches (case-insensitive)
      mylen = LEN_TRIM(citation_key)
      key_a = citation_key(1:mylen)
      CALL uppercase(key_a)
      match = 0
      DO I = 1, nbib - 1
         key_b = thebib(I)%ref%citation_key(1:mylen)
         CALL uppercase(key_b)
         IF (key_a == key_b) match = match + 1
      END DO
      IF (match > 0) citation_key = citation_key(1:mylen)//CHAR(ICHAR('a') + match)

      ! finally store it
      thebib(key)%ref%citation_key = citation_key

   END SUBROUTINE add_reference

! **************************************************************************************************
!> \brief deallocate the bibliography
!> \par History
!>      08.2007 Joost VandeVondele [ ]
! **************************************************************************************************
   SUBROUTINE remove_all_references()
      INTEGER                                            :: i

      DO i = 1, nbib
         IF (ASSOCIATED(thebib(i)%ref%ISI_record)) DEALLOCATE (thebib(i)%ref%ISI_record)
         thebib(i)%ref%DOI = ""

         DEALLOCATE (thebib(i)%ref)
      END DO
   END SUBROUTINE remove_all_references

! **************************************************************************************************
!> \brief printout of all cited references in the journal format sorted by publication year
!> \param unit ...
!> \par History
!>      08.2007 Joost VandeVondele
!>      07.2024 Ole Schuett
! **************************************************************************************************
   SUBROUTINE print_cited_references(unit)
      INTEGER, INTENT(IN)                                :: unit

      INTEGER                                            :: i
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: irank, ival

      ALLOCATE (ival(nbib), irank(nbib))

      ! we'll sort the references wrt to the publication year
      ! the most recent first, publications without a year get last
      DO i = 1, nbib
         irank(i) = i
         ival(i) = -get_epoch(thebib(i)%ref%ISI_record)
      END DO
      CALL sort(ival, nbib, irank)

      DO i = 1, nbib
         IF (thebib(irank(i))%ref%is_cited) THEN
            CALL print_reference_journal(key=irank(i), unit=unit)
            WRITE (unit, '(A)') ""
         END IF
      END DO

   END SUBROUTINE print_cited_references

! **************************************************************************************************
!> \brief prints a reference in a journal style citation format,
!>      adding also a DOI link, which is convenient
!> \param key ...
!> \param unit ...
!> \par History
!>      08.2007 created [Joost VandeVondele]
! **************************************************************************************************
   SUBROUTINE print_reference_journal(key, unit)
      INTEGER, INTENT(IN)                                :: key, unit

      CHARACTER(LEN=4*ISI_length)                        :: journal
      CHARACTER(LEN=ISI_length)                          :: author, title
      INTEGER                                            :: iauthor, ipos_line, ititle, line

! write the author list

      WRITE (unit, '(T2,A)', ADVANCE="NO") ""
      line = 1; iauthor = 0; ipos_line = 2
      author = get_next_author(thebib(key)%ref%ISI_record, line)
      DO WHILE (author .NE. "")
         iauthor = iauthor + 1
         IF (ipos_line + LEN_TRIM(author) > 71) THEN
            WRITE (unit, '(A)') ";"
            WRITE (unit, '(T2,A)', ADVANCE="NO") ""
            ipos_line = 2
         ELSE
            IF (iauthor .NE. 1) WRITE (unit, '(A)', ADVANCE="NO") "; "
            ipos_line = ipos_line + 2
         END IF
         WRITE (unit, '(A)', ADVANCE="NO") TRIM(author)
         ipos_line = ipos_line + LEN_TRIM(author)
         author = get_next_author(thebib(key)%ref%ISI_record, line)
      END DO
      IF (iauthor > 0) THEN
         WRITE (unit, '(A)', ADVANCE="NO") ". "
         ipos_line = ipos_line + 2
      END IF

      ! Journal, volume (issue), pages (year).
      journal = TRIM(get_source(thebib(key)%ref%ISI_record))
      IF (get_volume(thebib(key)%ref%ISI_record) .NE. "") THEN
         journal = TRIM(journal)//", "//get_volume(thebib(key)%ref%ISI_record)
         IF (get_issue(thebib(key)%ref%ISI_record) .NE. "") THEN
            journal = TRIM(journal)//" ("//TRIM(get_issue(thebib(key)%ref%ISI_record))//")"
         END IF
      END IF
      journal = TRIM(journal)//", "//get_pages(thebib(key)%ref%ISI_record)
      IF (get_year(thebib(key)%ref%ISI_record) .NE. "") THEN
         journal = TRIM(journal)//" ("//TRIM(get_year(thebib(key)%ref%ISI_record))//")."
      END IF
      IF (ipos_line + LEN_TRIM(journal) > 71) THEN
         WRITE (unit, '(A)') ""
         WRITE (unit, '(T2,A)', ADVANCE="NO") ""
         ipos_line = 2
      END IF
      IF (ipos_line + LEN_TRIM(journal) > 71) THEN
         WRITE (unit, '(A)') TRIM(journal(1:69))
         WRITE (unit, '(A)', ADVANCE="NO") TRIM(journal(69:))
      ELSE
         WRITE (unit, '(A)', ADVANCE="NO") TRIM(journal)
      END IF

      WRITE (unit, '(T2,A)') ""
      ! Title
      line = 1; ititle = 0
      title = get_next_title(thebib(key)%ref%ISI_record, line)
      DO WHILE (title .NE. "")
         ititle = ititle + 1
         IF (ititle .NE. 1) WRITE (unit, '(A)') ""
         WRITE (unit, '(T2,A)', ADVANCE="NO") TRIM(title)
         title = get_next_title(thebib(key)%ref%ISI_record, line)
      END DO
      IF (ititle > 0) WRITE (unit, '(A)') "."

      ! DOI
      IF (thebib(key)%ref%DOI .NE. "") THEN
         WRITE (unit, '(T2,A)') "https://doi.org/"//TRIM(thebib(key)%ref%DOI)
      END IF

   END SUBROUTINE print_reference_journal

! **************************************************************************************************
!> \brief Exports all references as XML.
!> \param unit ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE export_references_as_xml(unit)
      INTEGER, INTENT(IN)                                :: unit

      CHARACTER(LEN=ISI_length)                          :: author, title
      INTEGER                                            :: i, line

      DO i = 1, nbib
         WRITE (unit, '(T2,A)') '<REFERENCE key="'//TRIM(thebib(i)%ref%citation_key)//'">'

         ! Authors
         line = 1
         DO WHILE (.TRUE.)
            author = get_next_author(thebib(i)%ref%ISI_record, line)
            IF (LEN_TRIM(author) == 0) EXIT
            WRITE (unit, '(T3,A)') '<AUTHOR>'//TRIM(author)//'</AUTHOR>'
         END DO

         ! DOI, journal, volume (issue), pages (year).
         WRITE (unit, '(T3,A)') '<DOI>'//TRIM(substitute_special_xml_tokens(thebib(i)%ref%DOI))//'</DOI>'
         WRITE (unit, '(T3,A)') '<SOURCE>'//TRIM(get_source(thebib(i)%ref%ISI_record))//'</SOURCE>'
         WRITE (unit, '(T3,A)') '<VOLUME>'//TRIM(get_volume(thebib(i)%ref%ISI_record))//'</VOLUME>'
         WRITE (unit, '(T3,A)') '<ISSUE>'//TRIM(get_issue(thebib(i)%ref%ISI_record))//'</ISSUE>'
         WRITE (unit, '(T3,A)') '<PAGES>'//TRIM(get_pages(thebib(i)%ref%ISI_record))//'</PAGES>'
         WRITE (unit, '(T3,A)') '<YEAR>'//TRIM(get_year(thebib(i)%ref%ISI_record))//'</YEAR>'
         WRITE (unit, '(T3,A)') '<MONTH>'//TRIM(get_month(thebib(i)%ref%ISI_record))//'</MONTH>'
         WRITE (unit, '(T3,A)') '<DAY>'//TRIM(ADJUSTL(get_day(thebib(i)%ref%ISI_record)))//'</DAY>'

         ! Title
         line = 1
         title = get_next_title(thebib(i)%ref%ISI_record, line)
         WRITE (unit, '(T3,A)', ADVANCE='NO') '<TITLE>'//TRIM(title)
         DO WHILE (.TRUE.)
            title = ADJUSTL(get_next_title(thebib(i)%ref%ISI_record, line))
            IF (LEN_TRIM(title) == 0) EXIT
            WRITE (unit, '(1X,A)', ADVANCE="NO") TRIM(title)
         END DO
         WRITE (unit, '(A)') '</TITLE>'

         WRITE (unit, '(T2,A)') '</REFERENCE>'
      END DO

   END SUBROUTINE export_references_as_xml

! **************************************************************************************************
!> \brief returns the corresponding fields from an ISI record.
!>       returns an empty string if the field can not be found
!>       iline_start should be initialized to 1 to obtain the first matching entry
!>       on return it is updated, so that successive calls give successive fields
!> \param ISI_record ...
!> \param iline_start ...
!> \return ...
!> \par History
!>      08.2007 created [Joost VandeVondele]
! **************************************************************************************************
   FUNCTION get_next_author(ISI_record, iline_start) RESULT(res)
      CHARACTER(LEN=ISI_length), DIMENSION(:), &
         INTENT(IN)                                      :: ISI_record
      INTEGER, INTENT(INOUT)                             :: iline_start
      CHARACTER(LEN=ISI_length)                          :: res

      INTEGER                                            :: I, N
      LOGICAL                                            :: in_au_section

      res = ""
      in_au_section = .FALSE.
      N = SIZE(ISI_record, 1)
      IF (iline_start > N) RETURN
      line_loop: DO I = 1, N
         IF (ISI_record(I) (1:3) == "AU ") in_au_section = .TRUE.
         IF (in_au_section .AND. (ISI_record(I) (1:3) /= "AU " .AND. ISI_record(I) (1:3) /= "   ")) in_au_section = .FALSE.
         IF (in_au_section) THEN
            IF (I >= iline_start) THEN
               iline_start = I + 1
               res = ISI_record(I) (4:)
               EXIT line_loop
            END IF
         END IF
      END DO line_loop

      ! We might want to fixup the initials, adding a dot after each of them

   END FUNCTION get_next_author

! **************************************************************************************************
!> \brief ...
!> \param ISI_record ...
!> \param iline_start ...
!> \return ...
! **************************************************************************************************
   FUNCTION get_next_title(ISI_record, iline_start) RESULT(res)
      CHARACTER(LEN=ISI_length), DIMENSION(:), &
         INTENT(IN)                                      :: ISI_record
      INTEGER, INTENT(INOUT)                             :: iline_start
      CHARACTER(LEN=ISI_length)                          :: res

      INTEGER                                            :: I, N
      LOGICAL                                            :: in_ti_section

      res = ""

      in_ti_section = .FALSE.
      N = SIZE(ISI_record, 1)
      IF (iline_start > N) RETURN
      line_loop: DO I = 1, N
         IF (ISI_record(I) (1:3) == "TI ") in_ti_section = .TRUE.
         IF (in_ti_section .AND. (ISI_record(I) (1:3) /= "TI " .AND. ISI_record(I) (1:3) /= "   ")) in_ti_section = .FALSE.
         IF (in_ti_section) THEN
            IF (I >= iline_start) THEN
               iline_start = I + 1
               res = ISI_record(I) (4:)
               EXIT line_loop
            END IF
         END IF
      END DO line_loop

   END FUNCTION get_next_title

! **************************************************************************************************
!> \brief ...
!> \param ISI_record ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION get_source(ISI_record) RESULT(res)
      CHARACTER(LEN=ISI_length), DIMENSION(:), &
         INTENT(IN)                                      :: ISI_record
      CHARACTER(LEN=4*ISI_length)                        :: res

      INTEGER                                            :: I, J, N

      N = SIZE(ISI_record, 1)
      res = ""
      DO I = 1, N
         IF (ISI_record(I) (1:3) == "SO ") THEN
            res = ISI_record(I) (4:)
            DO J = I + 1, N
               IF (ISI_record(J) (1:3) == "   ") THEN
                  res = TRIM(res)//" "//ISI_record(J) (4:)
               ELSE
                  EXIT
               END IF
            END DO
            EXIT
         END IF
      END DO
   END FUNCTION get_source

! **************************************************************************************************
!> \brief ...
!> \param ISI_record ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION get_year(ISI_record) RESULT(res)
      CHARACTER(LEN=ISI_length), DIMENSION(:), &
         INTENT(IN)                                      :: ISI_record
      CHARACTER(LEN=ISI_length)                          :: res

      INTEGER                                            :: I, N

      N = SIZE(ISI_record, 1)
      res = ""
      DO I = 1, N
         IF (ISI_record(I) (1:3) == "PY ") res = ISI_record(I) (4:)
      END DO
   END FUNCTION get_year

! **************************************************************************************************
!> \brief ...
!> \param ISI_record ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION get_month(ISI_record) RESULT(res)
      CHARACTER(LEN=ISI_length), DIMENSION(:), &
         INTENT(IN)                                      :: ISI_record
      CHARACTER(LEN=ISI_length)                          :: res

      INTEGER                                            :: I, N

      N = SIZE(ISI_record, 1)
      res = ""
      DO I = 1, N
         IF (ISI_record(I) (1:3) == "PD ") res = ISI_record(I) (4:6)
      END DO
   END FUNCTION get_month

! **************************************************************************************************
!> \brief ...
!> \param ISI_record ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION get_day(ISI_record) RESULT(res)
      CHARACTER(LEN=ISI_length), DIMENSION(:), &
         INTENT(IN)                                      :: ISI_record
      CHARACTER(LEN=ISI_length)                          :: res

      INTEGER                                            :: D, I, N

      N = SIZE(ISI_record, 1)
      res = ""
      DO I = 1, N
         IF (ISI_record(I) (1:3) == "PD ") res = ISI_record(I) (7:)
      END DO
      ! PD can be e.g. OCT-NOV or OCT or OCT 27
      ! if res can't be read as an integer, it is not a day, and we bail out
      READ (res, *, ERR=998, END=998) D
      ! if the day is not in the expected range, we assume it is a parse error
      IF (D < 0 .OR. D > 31) res = ""
      RETURN
998   CONTINUE
      res = ""
   END FUNCTION get_day

! **************************************************************************************************
!> \brief ...
!> \param ISI_record ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION get_volume(ISI_record) RESULT(res)
      CHARACTER(LEN=ISI_length), DIMENSION(:), &
         INTENT(IN)                                      :: ISI_record
      CHARACTER(LEN=ISI_length)                          :: res

      INTEGER                                            :: I, N

      N = SIZE(ISI_record, 1)
      res = ""
      DO I = 1, N
         IF (ISI_record(I) (1:3) == "VL ") res = ISI_record(I) (4:)
      END DO
   END FUNCTION get_volume

! **************************************************************************************************
!> \brief ...
!> \param ISI_record ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION get_issue(ISI_record) RESULT(res)
      CHARACTER(LEN=ISI_length), DIMENSION(:), &
         INTENT(IN)                                      :: ISI_record
      CHARACTER(LEN=ISI_length)                          :: res

      INTEGER                                            :: I, N

      N = SIZE(ISI_record, 1)
      res = ""
      DO I = 1, N
         IF (ISI_record(I) (1:3) == "IS ") res = ISI_record(I) (4:)
      END DO
   END FUNCTION get_issue

! **************************************************************************************************
!> \brief ...
!> \param ISI_record ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION get_pages(ISI_record) RESULT(res)
      CHARACTER(LEN=ISI_length), DIMENSION(:), &
         INTENT(IN)                                      :: ISI_record
      CHARACTER(LEN=ISI_length)                          :: res

      CHARACTER(LEN=ISI_length)                          :: ar, bp, ep
      INTEGER                                            :: I, N

      N = SIZE(ISI_record, 1)
      res = ""
      bp = ""
      ep = ""
      ar = ""

      DO I = 1, N
         IF (ISI_record(I) (1:3) == "BP ") bp = ISI_record(I) (4:)
         IF (ISI_record(I) (1:3) == "EP ") ep = ISI_record(I) (4:)
         IF (ISI_record(I) (1:3) == "AR ") ar = ISI_record(I) (4:)
      END DO
      IF (bp .NE. "") THEN
         res = bp
         IF (ep .NE. "") res = TRIM(res)//"-"//ep
      END IF
      IF (res .EQ. "" .AND. ar .NE. "") res = ar
   END FUNCTION get_pages

! **************************************************************************************************
!> \brief ...
!> \param key ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION get_citation_key(key) RESULT(res)
      INTEGER, INTENT(IN)                                :: key
      CHARACTER(LEN=default_string_length)               :: res

      res = thebib(key)%ref%citation_key
   END FUNCTION get_citation_key

!
! This returns something epoch like, but can only be used to order the records
! missing years, months, days are implied zero(1900)
!
! **************************************************************************************************
!> \brief ...
!> \param ISI_record ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION get_epoch(ISI_record) RESULT(res)
      CHARACTER(LEN=ISI_length), DIMENSION(:), &
         INTENT(IN)                                      :: ISI_record
      INTEGER                                            :: res

      CHARACTER(LEN=ISI_length)                          :: tmp
      INTEGER                                            :: day, istat, month, year

! read year

      tmp = get_year(ISI_record)
      READ (tmp, *, IOSTAT=istat) year
      IF (istat .NE. 0) year = 1900

      ! read day
      tmp = get_day(ISI_record)
      READ (tmp, *, IOSTAT=istat) day
      IF (istat .NE. 0) day = 0

      ! read month
      tmp = get_month(ISI_record)
      SELECT CASE (tmp)
      CASE ("JAN")
         month = 1
      CASE ("FEB")
         month = 2
      CASE ("MAR")
         month = 3
      CASE ("APR")
         month = 4
      CASE ("MAY")
         month = 5
      CASE ("JUN")
         month = 6
      CASE ("JUL")
         month = 7
      CASE ("AUG")
         month = 8
      CASE ("SEP")
         month = 9
      CASE ("OCT")
         month = 10
      CASE ("NOV")
         month = 11
      CASE ("DEC")
         month = 12
      CASE DEFAULT
         month = 0
      END SELECT

      res = day + 31*month + 12*31*(year - 1900)

   END FUNCTION get_epoch

END MODULE reference_manager
