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
!>      07.2024 [Ole Schuett]
!> \author Joost VandeVondele
! **************************************************************************************************
MODULE reference_manager
   USE kinds,                           ONLY: default_string_length
   USE message_passing,                 ONLY: mp_para_env_type
   USE string_utilities,                ONLY: integer_to_string,&
                                              substitute_special_xml_tokens,&
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

   TYPE reference_type
      PRIVATE
      CHARACTER(LEN=default_string_length), DIMENSION(:), ALLOCATABLE :: authors
      CHARACTER(LEN=:), ALLOCATABLE                                   :: title
      CHARACTER(LEN=:), ALLOCATABLE                                   :: source
      CHARACTER(LEN=:), ALLOCATABLE                                   :: volume
      CHARACTER(LEN=:), ALLOCATABLE                                   :: pages
      INTEGER                                                         :: year = 0
      CHARACTER(LEN=:), ALLOCATABLE                                   :: doi
      ! has this reference been cited in the program run
      LOGICAL                                                         :: is_cited = .FALSE.
      ! this is a citation key for output in the reference lists
      CHARACTER(LEN=default_string_length)                            :: citation_key = ""
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
!> \param authors ...
!> \param title ...
!> \param source ...
!> \param volume ...
!> \param pages ...
!> \param year ...
!> \param doi ...
!> \par History
!>      08.2007 created [Joost VandeVondele]
!>      07.2024 complete rewrite [Ole Schuett]
!> \note
!>      - see bibliography.F for it use.
! **************************************************************************************************
   SUBROUTINE add_reference(key, authors, title, source, volume, pages, year, doi)
      INTEGER, INTENT(OUT)                               :: key
      CHARACTER(LEN=*), DIMENSION(:), INTENT(IN)         :: authors
      CHARACTER(LEN=*), INTENT(IN)                       :: title, source
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: volume, pages
      INTEGER, INTENT(IN)                                :: year
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL             :: doi

      CHARACTER                                          :: tmp
      CHARACTER(LEN=default_string_length)               :: author, citation_key, key_a, key_b
      INTEGER                                            :: i, ires, match, mylen, periodloc

      IF (nbib + 1 > max_reference) CPABORT("increase max_reference")
      nbib = nbib + 1
      key = nbib

      ALLOCATE (thebib(key)%ref)

      ! Copy authors.
      ALLOCATE (thebib(key)%ref%authors(SIZE(authors)))
      DO i = 1, SIZE(authors)
         CPASSERT(LEN_TRIM(authors(i)) <= default_string_length)
         thebib(key)%ref%authors(i) = authors(i)
      END DO

      ! Copy mandatory attributes.
      thebib(key)%ref%title = TRIM(title)
      thebib(key)%ref%source = TRIM(source)
      thebib(key)%ref%year = year

      ! Copy optional attributes.
      IF (PRESENT(volume)) THEN
         thebib(key)%ref%volume = TRIM(volume)
      END IF
      IF (PRESENT(pages)) THEN
         thebib(key)%ref%pages = TRIM(pages)
      END IF
      IF (PRESENT(doi)) THEN
         thebib(key)%ref%doi = TRIM(doi)
      END IF

      ! construct a citation_key
      author = authors(1)
      periodloc = INDEX(author, '.', back=.TRUE.)
      IF (periodloc > 0) author = author(periodloc + 1:)
      CPASSERT(LEN_TRIM(author) > 0)
      WRITE (citation_key, '(A,I4)') TRIM(author), year

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
         ival(i) = -thebib(i)%ref%year
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

      CHARACTER(LEN=:), ALLOCATABLE                      :: text
      CHARACTER(LEN=default_string_length)               :: year_str
      INTEGER                                            :: iauthor

      ! Authors
      text = thebib(key)%ref%authors(1)
      DO iauthor = 2, SIZE(thebib(key)%ref%authors)
         text = TRIM(text)//", "//thebib(key)%ref%authors(iauthor)
      END DO
      CALL write_long_text(TRIM(text)//".", unit)

      ! Journal, volume, pages (year).
      text = thebib(key)%ref%source
      IF (ALLOCATED(thebib(key)%ref%volume)) THEN
         text = text//" "//thebib(key)%ref%volume
      END IF
      IF (ALLOCATED(thebib(key)%ref%pages)) THEN
         text = TRIM(text)//", "//thebib(key)%ref%pages
      END IF
      IF (thebib(key)%ref%year > 0) THEN
         CALL integer_to_string(thebib(key)%ref%year, year_str)
         text = TRIM(text)//" ("//TRIM(year_str)//")"
      END IF
      CALL write_long_text(TRIM(text)//".", unit)

      ! Title
      CALL write_long_text(thebib(key)%ref%title//".", unit)

      ! DOI
      IF (ALLOCATED(thebib(key)%ref%doi)) THEN
         WRITE (unit, '(T2,A)') "https://doi.org/"//TRIM(thebib(key)%ref%doi)
      END IF

   END SUBROUTINE print_reference_journal

! **************************************************************************************************
!> \brief Exports all references as XML.
!> \param unit ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE export_references_as_xml(unit)
      INTEGER, INTENT(IN)                                :: unit

      INTEGER                                            :: i, j

      DO i = 1, nbib
         WRITE (unit, '(T2,A)') '<REFERENCE key="'//TRIM(thebib(i)%ref%citation_key)//'">'

         ! Authors
         DO j = 1, SIZE(thebib(i)%ref%authors)
            WRITE (unit, '(T3,A)') '<AUTHOR>'//TRIM(thebib(i)%ref%authors(j))//'</AUTHOR>'
         END DO

         ! Title and source.
         WRITE (unit, '(T3,A)') '<TITLE>'//thebib(i)%ref%title//'</TITLE>'
         WRITE (unit, '(T3,A)') '<SOURCE>'//thebib(i)%ref%source//'</SOURCE>'

         ! DOI, volume, pages, year, month.
         IF (ALLOCATED(thebib(i)%ref%doi)) &
            WRITE (unit, '(T3,A)') '<DOI>'//TRIM(substitute_special_xml_tokens(thebib(i)%ref%doi))//'</DOI>'
         IF (ALLOCATED(thebib(i)%ref%volume)) &
            WRITE (unit, '(T3,A)') '<VOLUME>'//thebib(i)%ref%volume//'</VOLUME>'
         IF (ALLOCATED(thebib(i)%ref%pages)) &
            WRITE (unit, '(T3,A)') '<PAGES>'//thebib(i)%ref%pages//'</PAGES>'
         IF (thebib(i)%ref%year > 0) &
            WRITE (unit, '(T3,A,I4.4,A)') '<YEAR>', thebib(i)%ref%year, '</YEAR>'
         WRITE (unit, '(T2,A)') '</REFERENCE>'
      END DO

   END SUBROUTINE export_references_as_xml

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

! **************************************************************************************************
!> \brief Helper routine for print_reference_journal()
!> \param text ...
!> \param unit ...
!> \return ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE write_long_text(text, unit)
      CHARACTER(LEN=*), INTENT(IN)                       :: text
      INTEGER, INTENT(IN)                                :: unit

      INTEGER                                            :: a, b

      a = 1; b = -1
      DO WHILE (b < LEN(text))
         b = next_linebreak(text, pos=a, rowlen=78)
         WRITE (unit, '(T2,A)') text(a:b)
         a = b + 1
      END DO
   END SUBROUTINE write_long_text

! **************************************************************************************************
!> \brief Helper routine for write_long_text()
!> \param text ...
!> \param pos ...
!> \param rowlen ...
!> \return ...
!> \author Ole Schuett
! **************************************************************************************************
   FUNCTION next_linebreak(text, pos, rowlen) RESULT(ibreak)
      CHARACTER(LEN=*), INTENT(IN)                       :: text
      INTEGER, INTENT(IN)                                :: pos, rowlen
      INTEGER                                            :: ibreak

      INTEGER                                            :: i, n

      n = LEN_TRIM(text)
      IF (n - pos <= rowlen) THEN
         ibreak = n ! remaining text shorter than line
      ELSE
         i = INDEX(text(pos + 1:pos + 1 + rowlen), " ", BACK=.TRUE.)
         IF (i == 0) THEN
            ibreak = pos + rowlen - 1 ! no space found, break mid-word
         ELSE
            ibreak = pos + i ! break at space closest to rowlen
         END IF
      END IF
   END FUNCTION next_linebreak

END MODULE reference_manager
