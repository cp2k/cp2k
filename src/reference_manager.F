!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2007 - 2008  CP2K developers group                                 !
!-----------------------------------------------------------------------------!
!!****h* cp2k/reference_manager *
!!
!!   NAME
!!     reference_manager
!!
!!   FUNCTION
!!     provides a uniform framework to add references to CP2K
!!     cite and output these
!!
!!   NOTES
!!     references need to be input using the ISI citation format, because it is
!!     uniform, easy to parse, and can be exported for example from web of science
!!     furthermore, it can be easily converted to and from using the bibutils tools
!!     a collection of easy to use conversion programs that can be found at
!!     http://www.scripps.edu/~cdputnam/software/bibutils/
!!     by Chris Putnam
!!
!!     see thebibliography.F on how to add references easily
!!     
!!
!!   AUTHOR
!!     Joost VandeVondele
!!     
!!
!!   MODIFICATION HISTORY
!!     08.2007 [Joost VandeVondele]
!!
!!   SOURCE
!****************************************************************************
MODULE reference_manager
  USE kinds,                           ONLY: default_string_length
  USE util,                            ONLY: sort

  IMPLICIT NONE

  PUBLIC ::  print_reference, print_all_references, cite_reference

  INTEGER, PUBLIC, PARAMETER :: print_format_isi=101, &
                                print_format_journal=102, &
                                print_format_html=103

  PRIVATE

  ! maximum number of reference that can be added
  INTEGER, PARAMETER :: max_reference=1024

  ! storage of a reference
  INTEGER, PARAMETER :: doi_length=128
  INTEGER, PARAMETER :: ISI_length=128

  ! the way we store a reference, should remain fully private
  TYPE reference_type
       PRIVATE
       ! the reference in a format as returned by the web of science
       CHARACTER(LEN=ISI_length), DIMENSION(:), POINTER :: ISI_record
       ! the doi only, i.e. without "http://dx.doi.org/"
       CHARACTER(LEN=doi_length)                        :: DOI
       ! has this reference been cited in the program run
       LOGICAL                                          :: is_cited
       ! this is a citation key for output in the reference lists
       CHARACTER(LEN=ISI_length)                        :: citation_key
  END TYPE reference_type

  ! useful to build arrays
  TYPE reference_p_type
     TYPE(reference_type), POINTER :: ref
  END TYPE

  ! thebibliography
  INTEGER, SAVE :: nbib=0
  TYPE(reference_p_type), DIMENSION(max_reference) :: thebib

  PUBLIC :: add_reference, & ! use this one only in bibliography.F
            remove_all_references, & ! use only in f77_interface.F
            get_citation_key         ! a string key describing the reference (e.g. Kohn1965b) 

CONTAINS
!!****f* reference_manager/cite_reference *
!!
!!   NAME
!!     cite_reference
!!
!!   FUNCTION
!!     marks a given reference as cited. 
!!
!!   NOTES
!!
!!   INPUTS
!!    - citation key as returned from add_reference
!!
!!   MODIFICATION HISTORY
!!     XX.2007 created [ ]
!!
!!*** **********************************************************************
  SUBROUTINE cite_reference(key) 
    INTEGER, INTENT(IN)                      :: key

    IF (key<1 .OR. key>max_reference) STOP "citation key out of range" 

    ! set as cited
    thebib(key)%ref%is_cited=.TRUE.

  END SUBROUTINE

!!****f* reference_manager/add_reference *
!!
!!   NAME
!!     add_reference
!!
!!   FUNCTION
!!     add a reference to the bibliography
!!
!!   NOTES
!!     - see bibliography.F for it use.
!!     - the ISI record is space sensitive, in particular the first three characters need to be blank
!!       or contain a key indicating the record type. See the header of this file for tools
!!       that can convert e.g. bibtex or endnote files to the ISI format
!!     - DOI: provide the DOI without a link. The link will be automatically created as needed.
!!
!!   INPUTS
!!    - key : output, this handle is needed to cite this reference later
!!
!!   MODIFICATION HISTORY
!!     08.2007 created [Joost VandeVondele]
!!
!!*** **********************************************************************
  SUBROUTINE add_reference(key,ISI_record,DOI)
    INTEGER, INTENT(OUT)                     :: key
    CHARACTER(LEN=*), DIMENSION(:)           :: ISI_record
    CHARACTER(LEN=*)                         :: DOI

    CHARACTER(LEN=ISI_length)                :: author, citation_key, year
    INTEGER                                  :: commaloc, i, ires, line, &
                                                match, mylen, nlines

    IF (nbib+1>max_reference) STOP "increase max_reference" 
    nbib=nbib+1 
    key =nbib

    ! initialize reference to zero
    ALLOCATE(thebib(key)%ref)
    NULLIFY(thebib(key)%ref%ISI_record)
    thebib(key)%ref%DOI=""
    thebib(key)%ref%is_cited=.FALSE.

    ! Assign DOI
    thebib(key)%ref%DOI=DOI

    ! Assign ISI_record
    nlines=SIZE(ISI_record,1)
    ALLOCATE(thebib(key)%ref%ISI_record(nlines))
    thebib(key)%ref%ISI_record=ISI_record

    ! construct a citation_key
    line=1 
    author=get_next_author(thebib(key)%ref%ISI_record,line)
    commaloc=INDEX(author,',')
    IF (commaloc.GT.0) author=author(1:commaloc-1)
    year=get_year(thebib(key)%ref%ISI_record)
    citation_key=TRIM(author)//TRIM(year)
    IF (citation_key=="") citation_key="unknown"

    ! avoid special characters in names, just remove them
    mylen=LEN_TRIM(citation_key)
    ires=0
    DO I=1,mylen
       IF (INDEX("0123456789thequickbrownfoxjumpsoverthelazydogTHEQUICKBROWNFOXJUMPSOVERTHELAZYDOG",citation_key(i:i)).NE.0) THEN
          ires=ires+1
          citation_key(ires:ires)=citation_key(i:i)
       ENDIF
    ENDDO
    citation_key(ires+1:)=""

    ! avoid duplicates, search through the list for matches
    mylen=LEN_TRIM(citation_key)
    match=0
    DO I=1,nbib-1
       IF (thebib(I)%ref%citation_key(1:mylen)==citation_key(1:mylen)) match=match+1
    ENDDO
    IF (match>0) citation_key=citation_key(1:mylen)//CHAR(ICHAR('a')+match)

    ! finally store it 
    thebib(key)%ref%citation_key=citation_key

  END SUBROUTINE add_reference

!!****f* reference_manager/remove_all_references *
!!
!!   NAME
!!     remove_all_references
!!
!!   FUNCTION
!!     deallocate the bibliography
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     08.2007 Joost VandeVondele [ ]
!!
!!*** **********************************************************************
  SUBROUTINE remove_all_references()
    INTEGER                                  :: i

    DO i=1,nbib        
       IF (ASSOCIATED(thebib(i)%ref%ISI_record)) DEALLOCATE(thebib(i)%ref%ISI_record)
       thebib(i)%ref%DOI=""

       DEALLOCATE(thebib(i)%ref)
    ENDDO
  END SUBROUTINE remove_all_references
!****f* reference_manager/print_all_references *
!!
!!   NAME
!!     print_all_references
!!
!!   FUNCTION
!!     printout of all references in a specific format
!!     optionally printing only those that are actually cited
!!     during program execution
!!
!!   NOTES
!!
!!   INPUTS
!!     cited_only: print only those marked as cited
!!     sorted: sort entries most recent first according to the date, 
!!             otherways sort with respect to key
!!     format: see module parameters print_format_XXXXXXXXX
!!     list: optionally, output a sub-list only
!!
!!   MODIFICATION HISTORY
!!     08.2007 Joost VandeVondele [ ]
!!
!!*** **********************************************************************
  SUBROUTINE print_all_references(cited_only,sorted,FORMAT,unit,list)
    LOGICAL, INTENT(IN)                      :: cited_only, sorted
    INTEGER, INTENT(IN)                      :: FORMAT, unit
    INTEGER, DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: list

    INTEGER                                  :: I, irecord, nref
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: indx, irank, ival

! we'll sort the references wrt to the publication year
! the most recent first, publications without a year get last

    IF (PRESENT(list)) THEN
       nref=SIZE(list)
    ELSE
       nref=nbib
    ENDIF

    ALLOCATE(ival(nref))
    ALLOCATE(irank(nref))
    ALLOCATE(indx(nref))

    IF (PRESENT(list)) THEN
       indx=list
    ELSE
       DO I=1,nref
          indx(I)=I
       ENDDO
    ENDIF

    DO I=1,nref
       irank(I)=I
    ENDDO

    IF (sorted) THEN
        DO I=1,nref
           ival(I)=-get_epoch(thebib(indx(I))%ref%ISI_record)
        ENDDO
    ELSE
        DO I=1,nref
           ival(I)=indx(I)
        ENDDO
    ENDIF
    CALL sort(ival,nref,irank)

    IF (FORMAT.EQ.print_format_html) THEN
       WRITE(unit,'(A)') '<TABLE border="1">'
    ENDIF
    DO I=1,nref
       irecord=indx(irank(I))
       IF (.NOT. cited_only .OR. thebib(irecord)%ref%is_cited) THEN
          SELECT CASE(FORMAT)
          CASE(print_format_isi)
          CASE(print_format_journal)
            WRITE(unit,'(A)') ""
          CASE(print_format_html)
            WRITE(unit,'(A)') "<TR><TD>"//'['//TRIM(thebib(irecord)%ref%citation_key)//']'//"</TD><TD>"
          CASE DEFAULT
            STOP "print_reference: wrong format"
          END SELECT

          CALL print_reference(irecord,FORMAT,unit)

          SELECT CASE(FORMAT)
          CASE(print_format_isi)
          CASE(print_format_journal)
            WRITE(unit,'(A)') ""
          CASE(print_format_html)
            WRITE(unit,'(A)') '</TD></TR>'
          CASE DEFAULT
            STOP "print_reference: wrong format"
          END SELECT
       ENDIF
    ENDDO
    IF (FORMAT.EQ.print_format_html) THEN
       WRITE(unit,'(A)') "</TABLE>"
    ENDIF

  END SUBROUTINE print_all_references
!****f* reference_manager/print_reference *
!!
!!   NAME
!!     print_reference
!!
!!   FUNCTION
!!     printout of a specified reference to a given unit in a selectable format
!!
!!   NOTES
!!
!!   INPUTS
!!     key: as returned from add_reference
!!     format: see module parameters print_format_XXXXXXXXX
!!
!!   MODIFICATION HISTORY
!!     08.2007 Joost VandeVondele [ ]
!!
!!*** **********************************************************************
  SUBROUTINE print_reference(key,FORMAT,unit)
    INTEGER, INTENT(IN)                      :: key, FORMAT, unit

    INTEGER                                  :: I

    IF (key<1 .OR. key>max_reference) STOP "citation key out of range" 

    SELECT CASE(FORMAT)
    CASE(print_format_isi)
      DO I=1,SIZE(thebib(key)%ref%ISI_record)
          WRITE(unit,'(T2,A)') TRIM(thebib(key)%ref%ISI_record(I))
      ENDDO
    CASE(print_format_journal)
      CALL print_reference_journal(key,unit)
    CASE(print_format_html)
      CALL print_reference_html(key,unit)
    CASE DEFAULT
      STOP "print_reference: wrong format"
    END SELECT
  END SUBROUTINE print_reference
!!****f* reference_manager/print_reference_journal *
!!
!!   NAME
!!     print_reference_journal
!!
!!   FUNCTION
!!     prints a reference in a journal style citation format,
!!     adding also a DOI link, which is convenient
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     08.2007 created [Joost VandeVondele]
!!
!!*** **********************************************************************
  SUBROUTINE print_reference_journal(key,unit)
    INTEGER, INTENT(IN)                      :: key, unit

    CHARACTER(LEN=ISI_length)                :: author, journal, title
    INTEGER                                  :: iauthor, ipos_line, ititle, &
                                                line

! write the author list

    WRITE(unit,'(T2,A)',ADVANCE="NO") ""
    line=1 ; iauthor=0 ; ipos_line=2
    author=get_next_author(thebib(key)%ref%ISI_record,line)
    DO WHILE (author.NE."")
       iauthor=iauthor+1
       IF (ipos_line+LEN_TRIM(author)>71) THEN
          WRITE(unit,'(A)') ";"
          WRITE(unit,'(T2,A)',ADVANCE="NO") ""
          ipos_line=2
       ELSE
          IF (iauthor.NE.1) WRITE(unit,'(A)',ADVANCE="NO") "; "
          ipos_line=ipos_line+2
       ENDIF
       WRITE(unit,'(A)',ADVANCE="NO") TRIM(author)
       ipos_line=ipos_line+LEN_TRIM(author)
       author=get_next_author(thebib(key)%ref%ISI_record,line)
    ENDDO
    IF (iauthor>0) THEN
        WRITE(unit,'(A)',ADVANCE="NO") ". "
        ipos_line=ipos_line+2
    ENDIF

    ! Journal, volume (issue), pages (year).
    journal=TRIM(get_source(thebib(key)%ref%ISI_record))//","
    journal=TRIM(journal)//" "//get_volume(thebib(key)%ref%ISI_record)
    IF (get_issue(thebib(key)%ref%ISI_record).NE."") THEN
        journal=TRIM(journal)//" ("//TRIM(get_issue(thebib(key)%ref%ISI_record))//")"
    ENDIF
    journal=TRIM(journal)//", "//get_pages(thebib(key)%ref%ISI_record)
    IF (get_year(thebib(key)%ref%ISI_record).NE."") THEN
        journal=TRIM(journal)//" ("//TRIM(get_year(thebib(key)%ref%ISI_record))//")."
    ENDIF
    IF (ipos_line+LEN_TRIM(journal)>71) THEN
       WRITE(unit,'(A)') "" 
       WRITE(unit,'(T2,A)',ADVANCE="NO") ""
       ipos_line=2
    ENDIF
    WRITE(unit,'(A)',ADVANCE="NO") TRIM(journal)

    WRITE(unit,'(T2,A)') ""
    ! Title
    line=1 ; ititle=0 
    title=get_next_title(thebib(key)%ref%ISI_record,line)
    DO WHILE (title.NE."")
       ititle=ititle+1
       IF (ititle.NE.1) WRITE(unit,'(A)') ""
       WRITE(unit,'(T2,A)', ADVANCE="NO") TRIM(title)
       title=get_next_title(thebib(key)%ref%ISI_record,line)
    ENDDO
    IF (ititle>0) WRITE(unit,'(A)') "."

    ! DOI
    IF (thebib(key)%ref%DOI .NE. "") THEN
      WRITE(unit,'(T2,A)') "http://dx.doi.org/"//TRIM(thebib(key)%ref%DOI)
    ENDIF
    
  END SUBROUTINE print_reference_journal
!!****f* reference_manager/print_reference_html *
!!
!!   NAME
!!     print_reference_html
!!
!!   FUNCTION
!!     prints a reference in a journal style citation format,
!!     adding 'beautifying' html tags, and a link to the journal
!!     using the DOI
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     08.2007 created [Joost VandeVondele]
!!
!!*** **********************************************************************
  SUBROUTINE print_reference_html(key,unit)
    INTEGER, INTENT(IN)                      :: key, unit

    CHARACTER(LEN=ISI_length)                :: author, journal, title
    INTEGER                                  :: iauthor, ititle, line

! write the author list

    WRITE(unit,'(T2,A,I0,A)',ADVANCE="NO") '<A NAME="reference_',key,'">'
    line=1 ; iauthor=0 
    author=get_next_author(thebib(key)%ref%ISI_record,line)
    DO WHILE (author.NE."")
       iauthor=iauthor+1
       IF (iauthor.NE.1) WRITE(unit,'(A)',ADVANCE="NO") "; "
       WRITE(unit,'(A)',ADVANCE="NO") TRIM(author)
       author=get_next_author(thebib(key)%ref%ISI_record,line)
    ENDDO
    IF (iauthor>0) WRITE(unit,'(A)') ".<br>"

    ! DOI
    IF (thebib(key)%ref%DOI .NE. "") THEN
      WRITE(unit,'(T2,A)',ADVANCE="NO") '<A HREF="http://dx.doi.org/'//TRIM(thebib(key)%ref%DOI)//'">'
    ENDIF
    ! Journal, volume (issue), pages (year).
    journal=TRIM(get_source(thebib(key)%ref%ISI_record))//","
    journal=TRIM(journal)//" "//get_volume(thebib(key)%ref%ISI_record)
    IF (get_issue(thebib(key)%ref%ISI_record).NE."") THEN
        journal=TRIM(journal)//" ("//TRIM(get_issue(thebib(key)%ref%ISI_record))//")"
    ENDIF
    journal=TRIM(journal)//", "//get_pages(thebib(key)%ref%ISI_record)
    IF (get_year(thebib(key)%ref%ISI_record).NE."") THEN
        journal=TRIM(journal)//" ("//TRIM(get_year(thebib(key)%ref%ISI_record))//")."
    ENDIF
    WRITE(unit,'(A)',ADVANCE="NO") TRIM(journal)
    IF (thebib(key)%ref%DOI .NE. "") THEN
      WRITE(unit,'(A)',ADVANCE="NO") '</A>'
    ENDIF
    WRITE(unit,'(A)') "</A><br>"

    ! Title
    line=1 ; ititle=0 
    title=get_next_title(thebib(key)%ref%ISI_record,line)
    DO WHILE (title.NE."")
       ititle=ititle+1
       IF (ititle.NE.1) WRITE(unit,'(A)') ""
       WRITE(unit,'(T2,A)', ADVANCE="NO") TRIM(title)
       title=get_next_title(thebib(key)%ref%ISI_record,line)
    ENDDO
    IF (ititle>0) WRITE(unit,'(A)') "."

    
  END SUBROUTINE print_reference_html

!!****f* reference_manager/XXXX *
!!
!!   NAME
!!      get_next_author, get_next_title, get_source, get_year, get_pages,
!!      get_volume, get_issue
!!
!!   FUNCTION
!!      returns the corresponding fields from an ISI record.
!!      returns an empty string if the field can not be found
!!      iline_start should be initialized to 1 to obtain the first matching entry
!!      on return it is updated, so that succesive calls give succesive fields
!!
!!   NOTES
!!
!!   INPUTS
!!
!!   MODIFICATION HISTORY
!!     08.2007 created [Joost VandeVondele]
!!
!!*** **********************************************************************
  FUNCTION get_next_author(ISI_record,iline_start) RESULT(res)
    CHARACTER(LEN=ISI_length), &
      DIMENSION(:), POINTER                  :: ISI_record
    INTEGER, INTENT(INOUT)                   :: iline_start
    CHARACTER(LEN=ISI_length)                :: res

    INTEGER                                  :: I, N
    LOGICAL                                  :: in_au_section

    res=""
    in_au_section=.FALSE.
    N=SIZE(ISI_record,1)
    IF (iline_start>N) RETURN
    line_loop: DO I=1,N
       IF (ISI_record(I)(1:3)=="AU ") in_au_section=.TRUE.
       IF (in_au_section .AND. (ISI_record(I)(1:3)/="AU " .AND. ISI_record(I)(1:3)/="   ")) in_au_section=.FALSE.
       IF (in_au_section) THEN
          IF (I>=iline_start) THEN
             iline_start=I+1
             res=ISI_record(I)(4:)
             EXIT line_loop
          ENDIF
       ENDIF
    ENDDO line_loop
 
    ! We might want to fixup the initials, adding a dot after each of them

  END FUNCTION get_next_author

  FUNCTION get_next_title(ISI_record,iline_start) RESULT(res)
    CHARACTER(LEN=ISI_length), &
      DIMENSION(:), POINTER                  :: ISI_record
    INTEGER, INTENT(INOUT)                   :: iline_start
    CHARACTER(LEN=ISI_length)                :: res

    INTEGER                                  :: I, N
    LOGICAL                                  :: in_ti_section

    res=""

    in_ti_section=.FALSE.
    N=SIZE(ISI_record,1)
    IF (iline_start>N) RETURN
    line_loop: DO I=1,N
       IF (ISI_record(I)(1:3)=="TI ") in_ti_section=.TRUE.
       IF (in_ti_section .AND. (ISI_record(I)(1:3)/="TI " .AND. ISI_record(I)(1:3)/="   ")) in_ti_section=.FALSE.
       IF (in_ti_section) THEN
          IF (I>=iline_start) THEN
             iline_start=I+1
             res=ISI_record(I)(4:)
             EXIT line_loop
          ENDIF
       ENDIF
    ENDDO line_loop

  END FUNCTION get_next_title

  FUNCTION get_source(ISI_record) RESULT(res)
    CHARACTER(LEN=ISI_length), &
      DIMENSION(:), POINTER                  :: ISI_record
    CHARACTER(LEN=ISI_length)                :: res

    INTEGER                                  :: I, N

    N=SIZE(ISI_record,1)
    res=""
    DO I=1,N
       IF (ISI_record(I)(1:3)=="SO ") res=ISI_record(I)(4:)
    ENDDO
  END FUNCTION get_source

  FUNCTION get_year(ISI_record) RESULT(res)
    CHARACTER(LEN=ISI_length), &
      DIMENSION(:), POINTER                  :: ISI_record
    CHARACTER(LEN=ISI_length)                :: res

    INTEGER                                  :: I, N

    N=SIZE(ISI_record,1)
    res=""
    DO I=1,N
       IF (ISI_record(I)(1:3)=="PY ") res=ISI_record(I)(4:)
    ENDDO
  END FUNCTION get_year

  FUNCTION get_month(ISI_record) RESULT(res)
    CHARACTER(LEN=ISI_length), &
      DIMENSION(:), POINTER                  :: ISI_record
    CHARACTER(LEN=ISI_length)                :: res

    INTEGER                                  :: I, N

    N=SIZE(ISI_record,1)
    res=""
    DO I=1,N
       IF (ISI_record(I)(1:3)=="PD ") res=ISI_record(I)(4:6)
    ENDDO
  END FUNCTION get_month

  FUNCTION get_day(ISI_record) RESULT(res)
    CHARACTER(LEN=ISI_length), &
      DIMENSION(:), POINTER                  :: ISI_record
    CHARACTER(LEN=ISI_length)                :: res

    INTEGER                                  :: D, I, N

    N=SIZE(ISI_record,1)
    res=""
    DO I=1,N
       IF (ISI_record(I)(1:3)=="PD ") res=ISI_record(I)(7:)
    ENDDO
    ! PD can be e.g. OCT-NOV or OCT or OCT 27
    ! if res can't be read as an integer, it is not a day, and we bail out
    READ(res,*,ERR=998,END=998) D
    RETURN
998 CONTINUE
    res=""
  END FUNCTION get_day

  FUNCTION get_volume(ISI_record) RESULT(res)
    CHARACTER(LEN=ISI_length), &
      DIMENSION(:), POINTER                  :: ISI_record
    CHARACTER(LEN=ISI_length)                :: res

    INTEGER                                  :: I, N

    N=SIZE(ISI_record,1)
    res=""
    DO I=1,N
       IF (ISI_record(I)(1:3)=="VL ") res=ISI_record(I)(4:)
    ENDDO
  END FUNCTION get_volume

  FUNCTION get_issue(ISI_record) RESULT(res)
    CHARACTER(LEN=ISI_length), &
      DIMENSION(:), POINTER                  :: ISI_record
    CHARACTER(LEN=ISI_length)                :: res

    INTEGER                                  :: I, N

    N=SIZE(ISI_record,1)
    res=""
    DO I=1,N
       IF (ISI_record(I)(1:3)=="IS ") res=ISI_record(I)(4:)
    ENDDO
  END FUNCTION get_issue

  FUNCTION get_pages(ISI_record) RESULT(res)
    CHARACTER(LEN=ISI_length), &
      DIMENSION(:), POINTER                  :: ISI_record
    CHARACTER(LEN=ISI_length)                :: res

    CHARACTER(LEN=ISI_length)                :: ar, bp, ep
    INTEGER                                  :: I, N

    N=SIZE(ISI_record,1)
    res=""
    bp=""
    ep=""
    ar=""

    DO I=1,N
       IF (ISI_record(I)(1:3)=="BP ") bp=ISI_record(I)(4:)
       IF (ISI_record(I)(1:3)=="EP ") ep=ISI_record(I)(4:)
       IF (ISI_record(I)(1:3)=="AR ") ar=ISI_record(I)(4:)
    ENDDO
    IF (bp.NE."") THEN
       res=bp
       IF (ep.NE."") res=TRIM(res)//"-"//ep
    ENDIF
    IF (res.EQ."".AND.ar.NE."") res=ar
  END FUNCTION get_pages

  FUNCTION get_citation_key(key) RESULT(res)
    INTEGER                                  :: key
    CHARACTER(LEN=default_string_length)     :: res

    res=thebib(key)%ref%citation_key
  END FUNCTION get_citation_key

!
! This returns something epoch like, but can only be used to order the records
! missing years, months, days are implied zero(1900)
!
  FUNCTION get_epoch(ISI_record) RESULT(res)
    CHARACTER(LEN=ISI_length), &
      DIMENSION(:), POINTER                  :: ISI_record
    INTEGER                                  :: res

    CHARACTER(LEN=ISI_length)                :: tmp
    INTEGER                                  :: day, istat, month, year

! read year

    tmp=get_year(ISI_record)
    READ(tmp,*,IOSTAT=istat) year
    IF (istat.NE.0) year=1900

    ! read day
    tmp=get_day(ISI_record)
    READ(tmp,*,IOSTAT=istat) day
    IF (istat.NE.0) day=0

    ! read month
    tmp=get_month(ISI_record)
    SELECT CASE(tmp)
    CASE("JAN")
      month=1
    CASE("FEB")
      month=2
    CASE("MAR")
      month=3
    CASE("APR")
      month=4
    CASE("MAY")
      month=5
    CASE("JUN")
      month=6
    CASE("JUL")
      month=7
    CASE("AUG")
      month=8
    CASE("SEP")
      month=9
    CASE("OCT")
      month=10
    CASE("NOV")
      month=11
    CASE("DEC")
      month=12
    CASE DEFAULT
      month=0
    END SELECT

    res=day+31*month+12*31*(year-1900)

  END FUNCTION get_epoch

END MODULE reference_manager
