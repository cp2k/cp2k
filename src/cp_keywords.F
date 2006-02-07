!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2004  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!!****s* cp2k/cp_keywords [1.0] *
!!
!!   NAME
!!     cp_keywords
!!
!!   FUNCTION
!!     Defines a structure for defining keywords in the input, listing them
!!     and providing basic help.
!!
!!   AUTHOR
!!     Joost VandeVondele
!!
!!   HISTORY
!!     Creation [05.2004]
!!
!!   SOURCE
!******************************************************************************

MODULE cp_keywords
  IMPLICIT NONE
  TYPE keyword_type
       CHARACTER(LEN=80)  :: keyword          ! name of the keyword  : "MAXSCF"
       CHARACTER(LEN=100) :: usage            ! how to use it        : "MAXSCF 10"
       CHARACTER(LEN=400) :: description      ! what does it do ?    : "MAXSCF : determines the maximum number of steps in an SCF run"
       CHARACTER(LEN=80)  :: variant1         ! alternate names of keywords
       CHARACTER(LEN=80)  :: variant2         ! alternate names of keywords
       CHARACTER(LEN=80)  :: variant3         ! alternate names of keywords
  END TYPE

  PRIVATE
  PUBLIC :: list_keywords, find_keyword, keyword_type

CONTAINS

! print a list of of keywords 
  SUBROUTINE list_keywords(keywords,unit)  
    TYPE(keyword_type), DIMENSION(:)         :: keywords
    INTEGER                                  :: unit

    INTEGER                                  :: i

     WRITE(unit,*)
     WRITE(unit,*) " +++ KNOWN KEYWORDS FOR THIS SECTION +++ "
     DO I=1,SIZE(keywords)
        WRITE(unit,*) "                                        ---"//TRIM(keywords(I)%keyword)//"---"
        WRITE(unit,*) "usage       : "//TRIM(keywords(I)%usage)
        WRITE(unit,*) "description : "//TRIM(keywords(I)%description)
        IF (TRIM(keywords(I)%variant1).NE."" .OR. TRIM(keywords(I)%variant2).NE."" .OR. TRIM(keywords(I)%variant3).NE."") THEN
            WRITE(unit,*) "variants    : ",TRIM(keywords(I)%variant1)," ",TRIM(keywords(I)%variant2)," ",TRIM(keywords(I)%variant3)
        ENDIF
     ENDDO
     WRITE(unit,*) " +++ END KNOWN KEYWORDS FOR THIS SECTION +++"
     WRITE(unit,*)
  END SUBROUTINE list_keywords

! return the location of the matching entry
  FUNCTION find_keyword(keywords,key)
    TYPE(keyword_type), DIMENSION(:)         :: keywords
    CHARACTER(LEN=*)                         :: key
    INTEGER                                  :: find_keyword

    INTEGER                                  :: I

     find_keyword=-1
     DO I=1,SIZE(keywords)
        IF (TRIM(keywords(I)%keyword).EQ.TRIM(key))  find_keyword=I
        IF (TRIM(keywords(I)%variant1).EQ.TRIM(key)) find_keyword=I
        IF (TRIM(keywords(I)%variant2).EQ.TRIM(key)) find_keyword=I
        IF (TRIM(keywords(I)%variant3).EQ.TRIM(key)) find_keyword=I
     ENDDO
 END FUNCTION find_keyword

END MODULE  cp_keywords
