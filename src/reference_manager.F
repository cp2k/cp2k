!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2007  CP2K developers group                                 !
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
  USE util, ONLY: sort
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
  END TYPE reference_type

  ! useful to build arrays
  TYPE reference_p_type
     TYPE(reference_type), POINTER :: ref
  END TYPE

  ! thebibliography
  INTEGER, SAVE :: nbib=0
  TYPE(reference_p_type), DIMENSION(max_reference) :: thebib

  PUBLIC :: add_reference, & ! use this one only in bibliography.F
            remove_all_references ! use only in f77_interface.F

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
    INTEGER, INTENT(IN) :: key

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
    INTEGER, INTENT(OUT) :: key
    CHARACTER(LEN=*), DIMENSION(:)  :: ISI_record
    CHARACTER(LEN=*)                :: DOI

    INTEGER :: nlines

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
    INTEGER :: i

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
  SUBROUTINE print_all_references(cited_only,sorted,format,unit,list)
    LOGICAL, INTENT(IN) :: cited_only,sorted
    INTEGER, INTENT(IN) :: format
    INTEGER, INTENT(IN) :: unit
    INTEGER, INTENT(IN), DIMENSION(:), OPTIONAL :: list


    INTEGER :: I, tmp, irecord, nref
    INTEGER, DIMENSION(:), ALLOCATABLE :: ival,irank,indx
    CHARACTER(LEN=ISI_length)                        :: year

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
           year=get_year(thebib(indx(I))%ref%ISI_record) 
           ival(I)=+1
           read(year,*,ERR=999,END=999) tmp
           ival(I)=-tmp
999        continue ! error reading
        ENDDO
    ELSE
        DO I=1,nref
           ival(I)=indx(I)
        ENDDO
    ENDIF
    CALL sort(ival,nref,irank)

    IF (format.EQ.print_format_html) THEN
       write(unit,'(A)') "<ol>"
    ENDIF
    DO I=1,nref
       irecord=indx(irank(I))
       IF (.NOT. cited_only .OR. thebib(irecord)%ref%is_cited) THEN
          SELECT CASE(format)
          CASE(print_format_isi)
          CASE(print_format_journal)
            write(unit,'(A)') ""
          CASE(print_format_html)
            write(unit,'(A)') "<li>"
          CASE DEFAULT
            STOP "print_reference: wrong format"
          END SELECT
          CALL print_reference(irecord,format,unit)
       ENDIF
    ENDDO
    IF (format.EQ.print_format_html) THEN
       write(unit,'(A)') "</ol>"
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
  SUBROUTINE print_reference(key,format,unit)
    INTEGER, INTENT(IN) :: key
    INTEGER, INTENT(IN) :: format
    INTEGER, INTENT(IN) :: unit

    INTEGER :: I

    IF (key<1 .OR. key>max_reference) STOP "citation key out of range" 

    SELECT CASE(format)
    CASE(print_format_isi)
      DO I=1,SIZE(thebib(key)%ref%ISI_record)
          write(unit,'(T2,A)') TRIM(thebib(key)%ref%ISI_record(I))
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
    INTEGER, INTENT(IN) :: key
    INTEGER, INTENT(IN) :: unit

    CHARACTER(LEN=ISI_length)                        :: author,title,journal
    INTEGER :: line, ipos_line, iauthor, ititle

    ! write the author list
    write(unit,'(T2,A)',ADVANCE="NO") ""
    line=1 ; iauthor=0 ; ipos_line=2
    author=get_next_author(thebib(key)%ref%ISI_record,line)
    DO WHILE (author.NE."")
       iauthor=iauthor+1
       IF (ipos_line+LEN_TRIM(author)>71) THEN
          write(unit,'(A)') ";"
          write(unit,'(T2,A)',ADVANCE="NO") ""
          ipos_line=2
       ELSE
          IF (iauthor.NE.1) write(unit,'(A)',ADVANCE="NO") "; "
          ipos_line=ipos_line+2
       ENDIF
       write(unit,'(A)',ADVANCE="NO") TRIM(author)
       ipos_line=ipos_line+LEN_TRIM(author)
       author=get_next_author(thebib(key)%ref%ISI_record,line)
    ENDDO
    IF (iauthor>0) THEN
        write(unit,'(A)',ADVANCE="NO") ". "
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
       write(unit,'(A)') "" 
       write(unit,'(T2,A)',ADVANCE="NO") ""
       ipos_line=2
    ENDIF
    write(unit,'(A)',ADVANCE="NO") TRIM(journal)

    write(unit,'(T2,A)') ""
    ! Title
    line=1 ; ititle=0 
    title=get_next_title(thebib(key)%ref%ISI_record,line)
    DO WHILE (title.NE."")
       ititle=ititle+1
       IF (ititle.NE.1) write(unit,'(A)') ""
       write(unit,'(T2,A)', ADVANCE="NO") TRIM(title)
       title=get_next_title(thebib(key)%ref%ISI_record,line)
    ENDDO
    IF (ititle>0) write(unit,'(A)') "."

    ! DOI
    IF (thebib(key)%ref%DOI .NE. "") THEN
      write(unit,'(T2,A)') "http://dx.doi.org/"//TRIM(thebib(key)%ref%DOI)
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
    INTEGER, INTENT(IN) :: key
    INTEGER, INTENT(IN) :: unit

    CHARACTER(LEN=ISI_length)                        :: author,title,journal
    INTEGER :: line, ipos_line, iauthor, ititle

    ! write the author list
    write(unit,'(T2,A,I0,A)',ADVANCE="NO") '<A NAME="reference_',key,'">'
    line=1 ; iauthor=0 ; ipos_line=2
    author=get_next_author(thebib(key)%ref%ISI_record,line)
    DO WHILE (author.NE."")
       iauthor=iauthor+1
       IF (iauthor.NE.1) write(unit,'(A)',ADVANCE="NO") "; "
       write(unit,'(A)',ADVANCE="NO") TRIM(author)
       author=get_next_author(thebib(key)%ref%ISI_record,line)
    ENDDO
    IF (iauthor>0) write(unit,'(A)') ".<br>"

    ! DOI
    IF (thebib(key)%ref%DOI .NE. "") THEN
      write(unit,'(T2,A)',ADVANCE="NO") '<A HREF="http://dx.doi.org/'//TRIM(thebib(key)%ref%DOI)//'">'
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
    write(unit,'(A)',ADVANCE="NO") TRIM(journal)
    IF (thebib(key)%ref%DOI .NE. "") THEN
      write(unit,'(A)',ADVANCE="NO") '</A>'
    ENDIF
    write(unit,'(A)') "</A><br>"

    ! Title
    line=1 ; ititle=0 
    title=get_next_title(thebib(key)%ref%ISI_record,line)
    DO WHILE (title.NE."")
       ititle=ititle+1
       IF (ititle.NE.1) write(unit,'(A)') ""
       write(unit,'(T2,A)', ADVANCE="NO") TRIM(title)
       title=get_next_title(thebib(key)%ref%ISI_record,line)
    ENDDO
    IF (ititle>0) write(unit,'(A)') "."

    
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
    CHARACTER(LEN=ISI_length)                        :: res
    CHARACTER(LEN=ISI_length), DIMENSION(:), POINTER :: ISI_record
    INTEGER, INTENT(INOUT)                           :: iline_start

    INTEGER :: I,N
    LOGICAL :: in_au_section

    res=""
    in_au_section=.FALSE.
    N=SIZE(ISI_record,1)
    IF (iline_start>N) return
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
    CHARACTER(LEN=ISI_length)                        :: res
    CHARACTER(LEN=ISI_length), DIMENSION(:), POINTER :: ISI_record
    INTEGER, INTENT(INOUT)                           :: iline_start

    INTEGER :: I,N
    LOGICAL :: in_ti_section
    res=""

    in_ti_section=.FALSE.
    N=SIZE(ISI_record,1)
    IF (iline_start>N) return
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
    CHARACTER(LEN=ISI_length)                        :: res
    CHARACTER(LEN=ISI_length), DIMENSION(:), POINTER :: ISI_record

    INTEGER :: I,N
    N=SIZE(ISI_record,1)
    res=""
    DO I=1,N
       IF (ISI_record(I)(1:3)=="SO ") res=ISI_record(I)(4:)
    ENDDO
  END FUNCTION get_source

  FUNCTION get_year(ISI_record) RESULT(res)
    CHARACTER(LEN=ISI_length)                        :: res
    CHARACTER(LEN=ISI_length), DIMENSION(:), POINTER :: ISI_record

    INTEGER :: I,N
    N=SIZE(ISI_record,1)
    res=""
    DO I=1,N
       IF (ISI_record(I)(1:3)=="PY ") res=ISI_record(I)(4:)
    ENDDO
  END FUNCTION get_year

  FUNCTION get_volume(ISI_record) RESULT(res)
    CHARACTER(LEN=ISI_length)                        :: res
    CHARACTER(LEN=ISI_length), DIMENSION(:), POINTER :: ISI_record

    INTEGER :: I,N
    N=SIZE(ISI_record,1)
    res=""
    DO I=1,N
       IF (ISI_record(I)(1:3)=="VL ") res=ISI_record(I)(4:)
    ENDDO
  END FUNCTION get_volume

  FUNCTION get_issue(ISI_record) RESULT(res)
    CHARACTER(LEN=ISI_length)                        :: res
    CHARACTER(LEN=ISI_length), DIMENSION(:), POINTER :: ISI_record

    INTEGER :: I,N
    N=SIZE(ISI_record,1)
    res=""
    DO I=1,N
       IF (ISI_record(I)(1:3)=="IS ") res=ISI_record(I)(4:)
    ENDDO
  END FUNCTION get_issue

  FUNCTION get_pages(ISI_record) RESULT(res)
    CHARACTER(LEN=ISI_length)                        :: res
    CHARACTER(LEN=ISI_length), DIMENSION(:), POINTER :: ISI_record

    INTEGER :: I,N
    CHARACTER(LEN=ISI_length)                        :: bp,ep,ar
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
    IF (res.NE."".AND.ar.NE."") res=ar
  END FUNCTION get_pages

END MODULE reference_manager
