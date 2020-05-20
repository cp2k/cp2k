!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief represents keywords in an input
!> \par History
!>      06.2004 created, based on Joost cp_keywords proposal [fawzi]
!> \author fawzi
! **************************************************************************************************
MODULE input_keyword_types
   USE cp_units,                        ONLY: cp_unit_create,&
                                              cp_unit_desc,&
                                              cp_unit_release,&
                                              cp_unit_type
   USE input_enumeration_types,         ONLY: enum_create,&
                                              enum_release,&
                                              enum_retain,&
                                              enumeration_type
   USE input_val_types,                 ONLY: &
        char_t, enum_t, integer_t, lchar_t, logical_t, no_t, real_t, val_create, val_release, &
        val_retain, val_type, val_write, val_write_internal
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE print_messages,                  ONLY: print_message
   USE reference_manager,               ONLY: get_citation_key
   USE string_utilities,                ONLY: a2s,&
                                              compress,&
                                              substitute_special_xml_tokens,&
                                              typo_match,&
                                              uppercase
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_keyword_types'

   INTEGER, PARAMETER, PUBLIC :: usage_string_length = default_string_length*2
   INTEGER, SAVE, PRIVATE :: last_keyword_id = 0

   PUBLIC :: keyword_p_type, keyword_type, keyword_create, keyword_retain, &
             keyword_release, keyword_get, keyword_describe, &
             write_keyword_xml, keyword_typo_match

! **************************************************************************************************
!> \brief represent a pointer to a keyword (to make arrays of pointers)
!> \param keyword the pointer to the keyword
!> \author fawzi
! **************************************************************************************************
   TYPE keyword_p_type
      TYPE(keyword_type), POINTER :: keyword
   END TYPE keyword_p_type

! **************************************************************************************************
!> \brief represent a keyword in the input
!> \param names the names of the current keyword (at least one should be
!>        present) for example "MAXSCF"
!> \param location is where in the source code (file and line) the keyword is created
!> \param usage how to use it "MAXSCF 10"
!> \param description what does it do: "MAXSCF : determines the maximum
!>        number of steps in an SCF run"
!> \param deprecation_notice show this warning that the keyword is deprecated
!> \param citations references to literature associated with this keyword
!> \param type_of_var the type of keyword (controls how it is parsed)
!>        it can be one of: no_parse_t,logical_t, integer_t, real_t,
!>        char_t
!> \param n_var number of values that should be parsed (-1=unknown)
!> \param repeats if the keyword can be present more than once in the
!>        section
!> \param removed to trigger a CPABORT when encountered while parsing the input
!> \param enum enumeration that defines the mapping between integers and
!>        strings
!> \param unit the default unit this keyword is read in (to automatically
!>        convert to the internal cp2k units during parsing)
!> \param default_value the default value for the keyword
!> \param lone_keyword_value value to be used in presence of the keyword
!>        without any parameter
!> \note
!>      I have expressely avoided a format string for the type of keywords:
!>      they should easily map to basic types of fortran, if you need more
!>      information use a subsection. [fawzi]
!> \author Joost & fawzi
! **************************************************************************************************
   TYPE keyword_type
      INTEGER :: ref_count, id_nr
      CHARACTER(LEN=default_string_length), DIMENSION(:), POINTER:: names
      CHARACTER(LEN=usage_string_length) :: location
      CHARACTER(LEN=usage_string_length) :: usage
      CHARACTER, DIMENSION(:), POINTER :: description => null()
      CHARACTER(LEN=:), ALLOCATABLE :: deprecation_notice
      INTEGER, POINTER, DIMENSION(:) :: citations
      INTEGER :: type_of_var, n_var
      LOGICAL :: repeats, removed
      TYPE(enumeration_type), POINTER :: enum
      TYPE(cp_unit_type), POINTER :: unit
      TYPE(val_type), POINTER :: default_value
      TYPE(val_type), POINTER :: lone_keyword_value
   END TYPE keyword_type

CONTAINS

! **************************************************************************************************
!> \brief creates a keyword object
!> \param keyword the keyword object to be created
!> \param location from where in the source code keyword_create() is called
!> \param name the name of the keyword
!> \param description ...
!> \param usage ...
!> \param type_of_var ...
!> \param n_var ...
!> \param repeats ...
!> \param variants ...
!> \param default_val ...
!> \param default_l_val ...
!> \param default_r_val ...
!> \param default_lc_val ...
!> \param default_c_val ...
!> \param default_i_val ...
!> \param default_l_vals ...
!> \param default_r_vals ...
!> \param default_c_vals ...
!> \param default_i_vals ...
!> \param lone_keyword_val ...
!> \param lone_keyword_l_val ...
!> \param lone_keyword_r_val ...
!> \param lone_keyword_c_val ...
!> \param lone_keyword_i_val ...
!> \param lone_keyword_l_vals ...
!> \param lone_keyword_r_vals ...
!> \param lone_keyword_c_vals ...
!> \param lone_keyword_i_vals ...
!> \param enum_c_vals ...
!> \param enum_i_vals ...
!> \param enum ...
!> \param enum_strict ...
!> \param enum_desc ...
!> \param unit_str ...
!> \param citations ...
!> \param deprecation_notice ...
!> \param removed ...
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE keyword_create(keyword, location, name, description, usage, type_of_var, &
                             n_var, repeats, variants, default_val, &
                             default_l_val, default_r_val, default_lc_val, default_c_val, default_i_val, &
                             default_l_vals, default_r_vals, default_c_vals, default_i_vals, &
                             lone_keyword_val, lone_keyword_l_val, lone_keyword_r_val, lone_keyword_c_val, &
                             lone_keyword_i_val, lone_keyword_l_vals, lone_keyword_r_vals, &
                             lone_keyword_c_vals, lone_keyword_i_vals, enum_c_vals, enum_i_vals, &
                             enum, enum_strict, enum_desc, unit_str, citations, deprecation_notice, removed)
      TYPE(keyword_type), POINTER                        :: keyword
      CHARACTER(len=*), INTENT(in)                       :: location, name, description
      CHARACTER(len=*), INTENT(in), OPTIONAL             :: usage
      INTEGER, INTENT(in), OPTIONAL                      :: type_of_var, n_var
      LOGICAL, INTENT(in), OPTIONAL                      :: repeats
      CHARACTER(len=*), DIMENSION(:), INTENT(in), &
         OPTIONAL                                        :: variants
      TYPE(val_type), OPTIONAL, POINTER                  :: default_val
      LOGICAL, INTENT(in), OPTIONAL                      :: default_l_val
      REAL(KIND=DP), INTENT(in), OPTIONAL                :: default_r_val
      CHARACTER(len=*), INTENT(in), OPTIONAL             :: default_lc_val, default_c_val
      INTEGER, INTENT(in), OPTIONAL                      :: default_i_val
      LOGICAL, DIMENSION(:), INTENT(in), OPTIONAL        :: default_l_vals
      REAL(KIND=DP), DIMENSION(:), INTENT(in), OPTIONAL  :: default_r_vals
      CHARACTER(len=*), DIMENSION(:), INTENT(in), &
         OPTIONAL                                        :: default_c_vals
      INTEGER, DIMENSION(:), INTENT(in), OPTIONAL        :: default_i_vals
      TYPE(val_type), OPTIONAL, POINTER                  :: lone_keyword_val
      LOGICAL, INTENT(in), OPTIONAL                      :: lone_keyword_l_val
      REAL(KIND=DP), INTENT(in), OPTIONAL                :: lone_keyword_r_val
      CHARACTER(len=*), INTENT(in), OPTIONAL             :: lone_keyword_c_val
      INTEGER, INTENT(in), OPTIONAL                      :: lone_keyword_i_val
      LOGICAL, DIMENSION(:), INTENT(in), OPTIONAL        :: lone_keyword_l_vals
      REAL(KIND=DP), DIMENSION(:), INTENT(in), OPTIONAL  :: lone_keyword_r_vals
      CHARACTER(len=*), DIMENSION(:), INTENT(in), &
         OPTIONAL                                        :: lone_keyword_c_vals
      INTEGER, DIMENSION(:), INTENT(in), OPTIONAL        :: lone_keyword_i_vals
      CHARACTER(len=*), DIMENSION(:), INTENT(in), &
         OPTIONAL                                        :: enum_c_vals
      INTEGER, DIMENSION(:), INTENT(in), OPTIONAL        :: enum_i_vals
      TYPE(enumeration_type), OPTIONAL, POINTER          :: enum
      LOGICAL, INTENT(in), OPTIONAL                      :: enum_strict
      CHARACTER(len=*), DIMENSION(:), INTENT(in), &
         OPTIONAL                                        :: enum_desc
      CHARACTER(len=*), INTENT(in), OPTIONAL             :: unit_str
      INTEGER, DIMENSION(:), INTENT(in), OPTIONAL        :: citations
      CHARACTER(len=*), INTENT(in), OPTIONAL             :: deprecation_notice
      LOGICAL, INTENT(in), OPTIONAL                      :: removed

      CHARACTER(len=*), PARAMETER :: routineN = 'keyword_create', routineP = moduleN//':'//routineN

      INTEGER                                            :: i, n
      LOGICAL                                            :: check

      CPASSERT(.NOT. ASSOCIATED(keyword))
      ALLOCATE (keyword)
      keyword%ref_count = 1
      last_keyword_id = last_keyword_id + 1
      keyword%id_nr = last_keyword_id
      NULLIFY (keyword%unit)
      keyword%location = location
      keyword%removed = .FALSE.

      IF (PRESENT(variants)) THEN
         ALLOCATE (keyword%names(SIZE(variants) + 1))
         keyword%names(1) = name
         DO i = 1, SIZE(variants)
            keyword%names(i + 1) = variants(i)
         END DO
      ELSE
         ALLOCATE (keyword%names(1))
         keyword%names(1) = name
      END IF
      DO i = 1, SIZE(keyword%names)
         CALL uppercase(keyword%names(i))
      END DO

      IF (PRESENT(usage)) THEN
         CPASSERT(LEN_TRIM(usage) <= LEN(keyword%usage))
         keyword%usage = usage
      ELSE
         keyword%usage = ""
      END IF

      n = LEN_TRIM(description)
      ALLOCATE (keyword%description(n))
      DO i = 1, n
         keyword%description(i) = description(i:i)
      END DO

      IF (PRESENT(citations)) THEN
         ALLOCATE (keyword%citations(SIZE(citations, 1)))
         keyword%citations = citations
      ELSE
         NULLIFY (keyword%citations)
      ENDIF

      keyword%repeats = .FALSE.
      IF (PRESENT(repeats)) keyword%repeats = repeats

      NULLIFY (keyword%enum)
      IF (PRESENT(enum)) THEN
         keyword%enum => enum
         IF (ASSOCIATED(enum)) CALL enum_retain(enum)
      END IF
      IF (PRESENT(enum_i_vals)) THEN
         CPASSERT(PRESENT(enum_c_vals))
         CPASSERT(.NOT. ASSOCIATED(keyword%enum))
         CALL enum_create(keyword%enum, c_vals=enum_c_vals, i_vals=enum_i_vals, &
                          desc=enum_desc, strict=enum_strict)
      ELSE
         CPASSERT(.NOT. PRESENT(enum_c_vals))
      END IF

      NULLIFY (keyword%default_value, keyword%lone_keyword_value)
      IF (PRESENT(default_val)) THEN
         IF (PRESENT(default_l_val) .OR. PRESENT(default_l_vals) .OR. &
             PRESENT(default_i_val) .OR. PRESENT(default_i_vals) .OR. &
             PRESENT(default_r_val) .OR. PRESENT(default_r_vals) .OR. &
             PRESENT(default_c_val) .OR. PRESENT(default_c_vals)) &
            CPABORT("you should pass either default_val or a default value, not both")
         keyword%default_value => default_val
         IF (ASSOCIATED(default_val%enum)) THEN
            IF (ASSOCIATED(keyword%enum)) THEN
               CPASSERT(keyword%enum%id_nr == default_val%enum%id_nr)
            ELSE
               keyword%enum => default_val%enum
               CALL enum_retain(keyword%enum)
            END IF
         ELSE
            CPASSERT(.NOT. ASSOCIATED(keyword%enum))
         END IF
         CALL val_retain(default_val)
      END IF
      IF (.NOT. ASSOCIATED(keyword%default_value)) THEN
         CALL val_create(keyword%default_value, l_val=default_l_val, &
                         l_vals=default_l_vals, i_val=default_i_val, i_vals=default_i_vals, &
                         r_val=default_r_val, r_vals=default_r_vals, c_val=default_c_val, &
                         c_vals=default_c_vals, lc_val=default_lc_val, enum=keyword%enum)
      END IF

      keyword%type_of_var = keyword%default_value%type_of_var
      IF (keyword%default_value%type_of_var == no_t) THEN
         CALL val_release(keyword%default_value)
      END IF

      IF (keyword%type_of_var == no_t) THEN
         IF (PRESENT(type_of_var)) THEN
            keyword%type_of_var = type_of_var
         ELSE
            CALL cp_warn(__LOCATION__, &
                         "keyword "//TRIM(keyword%names(1))// &
                         " assumed undefined type by default")
         END IF
      ELSE IF (PRESENT(type_of_var)) THEN
         IF (keyword%type_of_var /= type_of_var) &
            CALL cp_abort(__LOCATION__, &
                          "keyword "//TRIM(keyword%names(1))// &
                          " has a type different from the type of the default_value")
         keyword%type_of_var = type_of_var
      END IF

      IF (keyword%type_of_var == no_t) THEN
         CALL val_create(keyword%default_value)
      END IF

      IF (PRESENT(lone_keyword_val)) THEN
         IF (PRESENT(lone_keyword_l_val) .OR. PRESENT(lone_keyword_l_vals) .OR. &
             PRESENT(lone_keyword_i_val) .OR. PRESENT(lone_keyword_i_vals) .OR. &
             PRESENT(lone_keyword_r_val) .OR. PRESENT(lone_keyword_r_vals) .OR. &
             PRESENT(lone_keyword_c_val) .OR. PRESENT(lone_keyword_c_vals)) &
            CALL cp_abort(__LOCATION__, &
                          "you should pass either lone_keyword_val or a lone_keyword value, not both")
         keyword%lone_keyword_value => lone_keyword_val
         CALL val_retain(lone_keyword_val)
         IF (ASSOCIATED(lone_keyword_val%enum)) THEN
            IF (ASSOCIATED(keyword%enum)) THEN
               IF (keyword%enum%id_nr /= lone_keyword_val%enum%id_nr) &
                  CPABORT("keyword%enum%id_nr==lone_keyword_val%enum%id_nr")
            ELSE
               IF (ASSOCIATED(keyword%lone_keyword_value)) THEN
                  CPABORT(".NOT. ASSOCIATED(keyword%lone_keyword_value)")
               END IF
               keyword%enum => lone_keyword_val%enum
               CALL enum_retain(keyword%enum)
            END IF
         ELSE
            CPASSERT(.NOT. ASSOCIATED(keyword%enum))
         END IF
      END IF
      IF (.NOT. ASSOCIATED(keyword%lone_keyword_value)) THEN
         CALL val_create(keyword%lone_keyword_value, l_val=lone_keyword_l_val, &
                         l_vals=lone_keyword_l_vals, i_val=lone_keyword_i_val, i_vals=lone_keyword_i_vals, &
                         r_val=lone_keyword_r_val, r_vals=lone_keyword_r_vals, c_val=lone_keyword_c_val, &
                         c_vals=lone_keyword_c_vals, enum=keyword%enum)
      END IF
      IF (ASSOCIATED(keyword%lone_keyword_value)) THEN
         IF (keyword%lone_keyword_value%type_of_var == no_t) THEN
            CALL val_release(keyword%lone_keyword_value)
         ELSE
            IF (keyword%lone_keyword_value%type_of_var /= keyword%type_of_var) &
               CPABORT("lone_keyword_value type incompatible with keyword type")
            ! lc_val cannot have lone_keyword_value!
            IF (keyword%type_of_var == enum_t) THEN
               IF (keyword%enum%strict) THEN
                  check = .FALSE.
                  DO i = 1, SIZE(keyword%enum%i_vals)
                     check = check .OR. (keyword%default_value%i_val(1) == keyword%enum%i_vals(i))
                  END DO
                  IF (.NOT. check) &
                     CPABORT("default value not in enumeration : "//keyword%names(1))
               ENDIF
            ENDIF
         END IF
      END IF

      keyword%n_var = 1
      IF (ASSOCIATED(keyword%default_value)) THEN
         SELECT CASE (keyword%default_value%type_of_var)
         CASE (logical_t)
            keyword%n_var = SIZE(keyword%default_value%l_val)
         CASE (integer_t)
            keyword%n_var = SIZE(keyword%default_value%i_val)
         CASE (enum_t)
            IF (keyword%enum%strict) THEN
               check = .FALSE.
               DO i = 1, SIZE(keyword%enum%i_vals)
                  check = check .OR. (keyword%default_value%i_val(1) == keyword%enum%i_vals(i))
               END DO
               IF (.NOT. check) &
                  CPABORT("default value not in enumeration : "//keyword%names(1))
            ENDIF
            keyword%n_var = SIZE(keyword%default_value%i_val)
         CASE (real_t)
            keyword%n_var = SIZE(keyword%default_value%r_val)
         CASE (char_t)
            keyword%n_var = SIZE(keyword%default_value%c_val)
         CASE (lchar_t)
            keyword%n_var = 1
         CASE (no_t)
            keyword%n_var = 0
         CASE default
            CPABORT("")
         END SELECT
      END IF
      IF (PRESENT(n_var)) keyword%n_var = n_var
      IF (keyword%type_of_var == lchar_t .AND. keyword%n_var /= 1) &
         CPABORT("arrays of lchar_t not supported : "//keyword%names(1))

      IF (PRESENT(unit_str)) THEN
         CALL cp_unit_create(keyword%unit, unit_str)
      END IF

      IF (PRESENT(deprecation_notice)) THEN
         keyword%deprecation_notice = TRIM(deprecation_notice)
      END IF

      IF (PRESENT(removed)) THEN
         keyword%removed = removed
      END IF
   END SUBROUTINE keyword_create

! **************************************************************************************************
!> \brief retains the given keyword (see doc/ReferenceCounting.html)
!> \param keyword the keyword to retain
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE keyword_retain(keyword)
      TYPE(keyword_type), POINTER                        :: keyword

      CHARACTER(len=*), PARAMETER :: routineN = 'keyword_retain', routineP = moduleN//':'//routineN

      CPASSERT(ASSOCIATED(keyword))
      CPASSERT(keyword%ref_count > 0)
      keyword%ref_count = keyword%ref_count + 1
   END SUBROUTINE keyword_retain

! **************************************************************************************************
!> \brief releases the given keyword (see doc/ReferenceCounting.html)
!> \param keyword the keyword to release
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE keyword_release(keyword)
      TYPE(keyword_type), POINTER                        :: keyword

      CHARACTER(len=*), PARAMETER :: routineN = 'keyword_release', &
         routineP = moduleN//':'//routineN

      IF (ASSOCIATED(keyword)) THEN
         CPASSERT(keyword%ref_count > 0)
         keyword%ref_count = keyword%ref_count - 1
         IF (keyword%ref_count == 0) THEN
            DEALLOCATE (keyword%names)
            DEALLOCATE (keyword%description)
            CALL val_release(keyword%default_value)
            CALL val_release(keyword%lone_keyword_value)
            CALL enum_release(keyword%enum)
            CALL cp_unit_release(keyword%unit)
            IF (ASSOCIATED(keyword%citations)) THEN
               DEALLOCATE (keyword%citations)
            ENDIF
            DEALLOCATE (keyword)
         END IF
      END IF
      NULLIFY (keyword)
   END SUBROUTINE keyword_release

! **************************************************************************************************
!> \brief ...
!> \param keyword ...
!> \param names ...
!> \param usage ...
!> \param description ...
!> \param type_of_var ...
!> \param n_var ...
!> \param default_value ...
!> \param lone_keyword_value ...
!> \param repeats ...
!> \param enum ...
!> \param citations ...
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE keyword_get(keyword, names, usage, description, type_of_var, n_var, &
                          default_value, lone_keyword_value, repeats, enum, citations)
      TYPE(keyword_type), POINTER                        :: keyword
      CHARACTER(len=default_string_length), &
         DIMENSION(:), OPTIONAL, POINTER                 :: names
      CHARACTER(len=*), INTENT(out), OPTIONAL            :: usage, description
      INTEGER, INTENT(out), OPTIONAL                     :: type_of_var, n_var
      TYPE(val_type), OPTIONAL, POINTER                  :: default_value, lone_keyword_value
      LOGICAL, INTENT(out), OPTIONAL                     :: repeats
      TYPE(enumeration_type), OPTIONAL, POINTER          :: enum
      INTEGER, DIMENSION(:), OPTIONAL, POINTER           :: citations

      CHARACTER(len=*), PARAMETER :: routineN = 'keyword_get', routineP = moduleN//':'//routineN

      CPASSERT(ASSOCIATED(keyword))
      CPASSERT(keyword%ref_count > 0)
      IF (PRESENT(names)) names => keyword%names
      IF (PRESENT(usage)) usage = keyword%usage
      IF (PRESENT(description)) description = a2s(keyword%description)
      IF (PRESENT(type_of_var)) type_of_var = keyword%type_of_var
      IF (PRESENT(n_var)) n_var = keyword%n_var
      IF (PRESENT(repeats)) repeats = keyword%repeats
      IF (PRESENT(default_value)) default_value => keyword%default_value
      IF (PRESENT(lone_keyword_value)) lone_keyword_value => keyword%lone_keyword_value
      IF (PRESENT(enum)) enum => keyword%enum
      IF (PRESENT(citations)) citations => keyword%citations
   END SUBROUTINE keyword_get

! **************************************************************************************************
!> \brief writes out a description of the keyword
!> \param keyword the keyword to describe
!> \param unit_nr the unit to write to
!> \param level the description level (0 no description, 1 name
!>        2: +usage, 3: +variants+description+default_value+repeats
!>        4: +type_of_var)
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE keyword_describe(keyword, unit_nr, level)
      TYPE(keyword_type), POINTER                        :: keyword
      INTEGER, INTENT(in)                                :: unit_nr, level

      CHARACTER(len=*), PARAMETER :: routineN = 'keyword_describe', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=default_string_length)               :: c_string
      INTEGER                                            :: i, l

      CPASSERT(ASSOCIATED(keyword))
      CPASSERT(keyword%ref_count > 0)
      IF (level > 0 .AND. (unit_nr > 0)) THEN
         WRITE (unit_nr, "(a,a,a)") "                           ---", &
            TRIM(keyword%names(1)), "---"
         IF (level > 1) THEN
            WRITE (unit_nr, "(a,a)") "usage         : ", TRIM(keyword%usage)
         END IF
         IF (level > 2) THEN
            WRITE (unit_nr, "(a)") "description   : "
            CALL print_message(TRIM(a2s(keyword%description)), unit_nr, 0, 0, 0)
            IF (level > 3) THEN
               SELECT CASE (keyword%type_of_var)
               CASE (logical_t)
                  IF (keyword%n_var == -1) THEN
                     WRITE (unit_nr, "('  A list of logicals is expected')")
                  ELSE IF (keyword%n_var == 1) THEN
                     WRITE (unit_nr, "('  A logical is expected')")
                  ELSE
                     WRITE (unit_nr, "(i6,'  logicals are expected')") keyword%n_var
                  END IF
                  WRITE (unit_nr, "('  (T,TRUE,YES,ON) and (F,FALSE,NO,OFF) are synonyms')")
               CASE (integer_t)
                  IF (keyword%n_var == -1) THEN
                     WRITE (unit_nr, "('  A list of integers is expected')")
                  ELSE IF (keyword%n_var == 1) THEN
                     WRITE (unit_nr, "('  An integer is expected')")
                  ELSE
                     WRITE (unit_nr, "(i6,' integers are expected')") keyword%n_var
                  END IF
               CASE (real_t)
                  IF (keyword%n_var == -1) THEN
                     WRITE (unit_nr, "('  A list of reals is expected')")
                  ELSE IF (keyword%n_var == 1) THEN
                     WRITE (unit_nr, "('  A real is expected')")
                  ELSE
                     WRITE (unit_nr, "(i6,' reals are expected')") keyword%n_var
                  END IF
                  IF (ASSOCIATED(keyword%unit)) THEN
                     c_string = cp_unit_desc(keyword%unit, accept_undefined=.TRUE.)
                     WRITE (unit_nr, "('the default unit of measure is ',a)") &
                        TRIM(c_string)
                  END IF
               CASE (char_t)
                  IF (keyword%n_var == -1) THEN
                     WRITE (unit_nr, "('  A list of words is expected')")
                  ELSE IF (keyword%n_var == 1) THEN
                     WRITE (unit_nr, "('  A word is expected')")
                  ELSE
                     WRITE (unit_nr, "(i6,' words are expected')") keyword%n_var
                  END IF
               CASE (lchar_t)
                  WRITE (unit_nr, "('  A string is expected')")
               CASE (enum_t)
                  IF (keyword%n_var == -1) THEN
                     WRITE (unit_nr, "('  A list of keywords is expected')")
                  ELSE IF (keyword%n_var == 1) THEN
                     WRITE (unit_nr, "('  A keyword is expected')")
                  ELSE
                     WRITE (unit_nr, "(i6,' keywords are expected')") keyword%n_var
                  END IF
               CASE (no_t)
                  WRITE (unit_nr, "('  Non-standard type.')")
               CASE default
                  CPABORT("")
               END SELECT
            END IF
            IF (keyword%type_of_var == enum_t) THEN
               IF (level > 3) THEN
                  WRITE (unit_nr, "('  valid keywords:')")
                  DO i = 1, SIZE(keyword%enum%c_vals)
                     c_string = keyword%enum%c_vals(i)
                     IF (LEN_TRIM(a2s(keyword%enum%desc(i)%chars)) > 0) THEN
                        WRITE (unit_nr, "('  - ',a,' : ',a,'.')") &
                           TRIM(c_string), TRIM(a2s(keyword%enum%desc(i)%chars))
                     ELSE
                        WRITE (unit_nr, "('  - ',a)") TRIM(c_string)
                     END IF
                  END DO
               ELSE
                  WRITE (unit_nr, "('  valid keywords:')", advance='NO')
                  l = 17
                  DO i = 1, SIZE(keyword%enum%c_vals)
                     c_string = keyword%enum%c_vals(i)
                     IF (l + LEN_TRIM(c_string) > 72 .AND. l > 14) THEN
                        WRITE (unit_nr, "(/,'    ')", advance='NO')
                        l = 4
                     END IF
                     WRITE (unit_nr, "(' ',a)", advance='NO') TRIM(c_string)
                     l = LEN_TRIM(c_string) + 3
                  END DO
                  WRITE (unit_nr, "()")
               END IF
               IF (.NOT. keyword%enum%strict) THEN
                  WRITE (unit_nr, "('     other integer values are also accepted.')")
               END IF
            END IF
            IF (ASSOCIATED(keyword%default_value) .AND. keyword%type_of_var /= no_t) THEN
               WRITE (unit_nr, "('default_value : ')", advance="NO")
               CALL val_write(keyword%default_value, unit_nr=unit_nr)
            END IF
            IF (ASSOCIATED(keyword%lone_keyword_value) .AND. keyword%type_of_var /= no_t) THEN
               WRITE (unit_nr, "('lone_keyword  : ')", advance="NO")
               CALL val_write(keyword%lone_keyword_value, unit_nr=unit_nr)
            END IF
            IF (keyword%repeats) THEN
               WRITE (unit_nr, "(' and it can be repeated more than once')", advance="NO")
            END IF
            WRITE (unit_nr, "()")
            IF (SIZE(keyword%names) > 1) THEN
               WRITE (unit_nr, "(a)", advance="NO") "variants    : "
               DO i = 2, SIZE(keyword%names)
                  WRITE (unit_nr, "(a,' ')", advance="NO") keyword%names(i)
               END DO
               WRITE (unit_nr, "()")
            ENDIF
         END IF
      END IF
   END SUBROUTINE keyword_describe

! **************************************************************************************************
!> \brief Prints a description of a keyword in XML format
!> \param keyword The keyword to describe
!> \param level ...
!> \param unit_number Number of the output unit
!> \author Matthias Krack
! **************************************************************************************************
   SUBROUTINE write_keyword_xml(keyword, level, unit_number)

      TYPE(keyword_type), POINTER                        :: keyword
      INTEGER, INTENT(IN)                                :: level, unit_number

      CHARACTER(LEN=*), PARAMETER :: routineN = 'write_keyword_xml', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=1000)                                :: string
      CHARACTER(LEN=3)                                   :: removed, repeats
      CHARACTER(LEN=8)                                   :: short_string
      INTEGER                                            :: i, l0, l1, l2, l3, l4

      CPASSERT(ASSOCIATED(keyword))
      CPASSERT(keyword%ref_count > 0)

      ! Indentation for current level, next level, etc.

      l0 = level
      l1 = level + 1
      l2 = level + 2
      l3 = level + 3
      l4 = level + 4

      IF (keyword%repeats) THEN
         repeats = "yes"
      ELSE
         repeats = "no "
      END IF

      IF (keyword%removed) THEN
         removed = "yes"
      ELSE
         removed = "no "
      END IF

      ! Write (special) keyword element

      IF (keyword%names(1) == "_SECTION_PARAMETERS_") THEN
         WRITE (UNIT=unit_number, FMT="(A)") &
            REPEAT(" ", l0)//"<SECTION_PARAMETERS repeats="""//TRIM(repeats)// &
            """ removed="""//TRIM(removed)//""">", &
            REPEAT(" ", l1)//"<NAME type=""default"">SECTION_PARAMETERS</NAME>"
      ELSE IF (keyword%names(1) == "_DEFAULT_KEYWORD_") THEN
         WRITE (UNIT=unit_number, FMT="(A)") &
            REPEAT(" ", l0)//"<DEFAULT_KEYWORD repeats="""//TRIM(repeats)//""">", &
            REPEAT(" ", l1)//"<NAME type=""default"">DEFAULT_KEYWORD</NAME>"
      ELSE
         WRITE (UNIT=unit_number, FMT="(A)") &
            REPEAT(" ", l0)//"<KEYWORD repeats="""//TRIM(repeats)// &
            """ removed="""//TRIM(removed)//""">", &
            REPEAT(" ", l1)//"<NAME type=""default"">"// &
            TRIM(keyword%names(1))//"</NAME>"
      END IF

      DO i = 2, SIZE(keyword%names)
         WRITE (UNIT=unit_number, FMT="(A)") &
            REPEAT(" ", l1)//"<NAME type=""alias"">"// &
            TRIM(keyword%names(i))//"</NAME>"
      END DO

      SELECT CASE (keyword%type_of_var)
      CASE (logical_t)
         WRITE (UNIT=unit_number, FMT="(A)") &
            REPEAT(" ", l1)//"<DATA_TYPE kind=""logical"">"
      CASE (integer_t)
         WRITE (UNIT=unit_number, FMT="(A)") &
            REPEAT(" ", l1)//"<DATA_TYPE kind=""integer"">"
      CASE (real_t)
         WRITE (UNIT=unit_number, FMT="(A)") &
            REPEAT(" ", l1)//"<DATA_TYPE kind=""real"">"
      CASE (char_t)
         WRITE (UNIT=unit_number, FMT="(A)") &
            REPEAT(" ", l1)//"<DATA_TYPE kind=""word"">"
      CASE (lchar_t)
         WRITE (UNIT=unit_number, FMT="(A)") &
            REPEAT(" ", l1)//"<DATA_TYPE kind=""string"">"
      CASE (enum_t)
         WRITE (UNIT=unit_number, FMT="(A)") &
            REPEAT(" ", l1)//"<DATA_TYPE kind=""keyword"">"
         IF (keyword%enum%strict) THEN
            WRITE (UNIT=unit_number, FMT="(A)") &
               REPEAT(" ", l2)//"<ENUMERATION strict=""yes"">"
         ELSE
            WRITE (UNIT=unit_number, FMT="(A)") &
               REPEAT(" ", l2)//"<ENUMERATION strict=""no"">"
         END IF
         DO i = 1, SIZE(keyword%enum%c_vals)
            WRITE (UNIT=unit_number, FMT="(A)") &
               REPEAT(" ", l3)//"<ITEM>", &
               REPEAT(" ", l4)//"<NAME>"// &
               TRIM(ADJUSTL(substitute_special_xml_tokens(keyword%enum%c_vals(i))))//"</NAME>", &
               REPEAT(" ", l4)//"<DESCRIPTION>"// &
               TRIM(ADJUSTL(substitute_special_xml_tokens(a2s(keyword%enum%desc(i)%chars)))) &
               //"</DESCRIPTION>", REPEAT(" ", l3)//"</ITEM>"
         END DO
         WRITE (UNIT=unit_number, FMT="(A)") REPEAT(" ", l2)//"</ENUMERATION>"
      CASE (no_t)
         WRITE (UNIT=unit_number, FMT="(A)") &
            REPEAT(" ", l1)//"<DATA_TYPE kind=""non-standard type"">"
      CASE DEFAULT
         CPABORT("")
      END SELECT

      short_string = ""
      WRITE (UNIT=short_string, FMT="(I8)") keyword%n_var
      WRITE (UNIT=unit_number, FMT="(A)") &
         REPEAT(" ", l2)//"<N_VAR>"//TRIM(ADJUSTL(short_string))//"</N_VAR>", &
         REPEAT(" ", l1)//"</DATA_TYPE>"

      WRITE (UNIT=unit_number, FMT="(A)") REPEAT(" ", l1)//"<USAGE>"// &
         TRIM(substitute_special_xml_tokens(keyword%usage)) &
         //"</USAGE>"

      WRITE (UNIT=unit_number, FMT="(A)") REPEAT(" ", l1)//"<DESCRIPTION>"// &
         TRIM(substitute_special_xml_tokens(a2s(keyword%description))) &
         //"</DESCRIPTION>"

      IF (ALLOCATED(keyword%deprecation_notice)) &
         WRITE (UNIT=unit_number, FMT="(A)") REPEAT(" ", l1)//"<DEPRECATION_NOTICE>"// &
         TRIM(substitute_special_xml_tokens(keyword%deprecation_notice)) &
         //"</DEPRECATION_NOTICE>"

      IF (ASSOCIATED(keyword%default_value) .AND. &
          (keyword%type_of_var /= no_t)) THEN
         IF (ASSOCIATED(keyword%unit)) THEN
            CALL val_write_internal(val=keyword%default_value, &
                                    string=string, &
                                    unit=keyword%unit)
         ELSE
            CALL val_write_internal(val=keyword%default_value, &
                                    string=string)
         END IF
         CALL compress(string)
         WRITE (UNIT=unit_number, FMT="(A)") &
            REPEAT(" ", l1)//"<DEFAULT_VALUE>"// &
            TRIM(ADJUSTL(substitute_special_xml_tokens(string)))//"</DEFAULT_VALUE>"
      END IF

      IF (ASSOCIATED(keyword%unit)) THEN
         string = cp_unit_desc(keyword%unit, accept_undefined=.TRUE.)
         WRITE (UNIT=unit_number, FMT="(A)") &
            REPEAT(" ", l1)//"<DEFAULT_UNIT>"// &
            TRIM(ADJUSTL(string))//"</DEFAULT_UNIT>"
      END IF

      IF (ASSOCIATED(keyword%lone_keyword_value) .AND. &
          (keyword%type_of_var /= no_t)) THEN
         CALL val_write_internal(val=keyword%lone_keyword_value, &
                                 string=string)
         WRITE (UNIT=unit_number, FMT="(A)") &
            REPEAT(" ", l1)//"<LONE_KEYWORD_VALUE>"// &
            TRIM(ADJUSTL(substitute_special_xml_tokens(string)))//"</LONE_KEYWORD_VALUE>"
      END IF

      IF (ASSOCIATED(keyword%citations)) THEN
         DO i = 1, SIZE(keyword%citations, 1)
            short_string = ""
            WRITE (UNIT=short_string, FMT="(I8)") keyword%citations(i)
            WRITE (UNIT=unit_number, FMT="(A)") &
               REPEAT(" ", l1)//"<REFERENCE>", &
               REPEAT(" ", l2)//"<NAME>"//TRIM(get_citation_key(keyword%citations(i)))//"</NAME>", &
               REPEAT(" ", l2)//"<NUMBER>"//TRIM(ADJUSTL(short_string))//"</NUMBER>", &
               REPEAT(" ", l1)//"</REFERENCE>"
         END DO
      END IF

      WRITE (UNIT=unit_number, FMT="(A)") &
         REPEAT(" ", l1)//"<LOCATION>"//TRIM(keyword%location)//"</LOCATION>"

      ! Close (special) keyword section

      IF (keyword%names(1) == "_SECTION_PARAMETERS_") THEN
         WRITE (UNIT=unit_number, FMT="(A)") &
            REPEAT(" ", l0)//"</SECTION_PARAMETERS>"
      ELSE IF (keyword%names(1) == "_DEFAULT_KEYWORD_") THEN
         WRITE (UNIT=unit_number, FMT="(A)") &
            REPEAT(" ", l0)//"</DEFAULT_KEYWORD>"
      ELSE
         WRITE (UNIT=unit_number, FMT="(A)") &
            REPEAT(" ", l0)//"</KEYWORD>"
      END IF

   END SUBROUTINE write_keyword_xml

! **************************************************************************************************
!> \brief ...
!> \param keyword ...
!> \param unknown_string ...
!> \param location_string ...
!> \param matching_rank ...
!> \param matching_string ...
!> \param bonus ...
! **************************************************************************************************
   SUBROUTINE keyword_typo_match(keyword, unknown_string, location_string, matching_rank, matching_string, bonus)

      TYPE(keyword_type), POINTER                        :: keyword
      CHARACTER(LEN=*)                                   :: unknown_string, location_string
      INTEGER, DIMENSION(:), INTENT(INOUT)               :: matching_rank
      CHARACTER(LEN=*), DIMENSION(:), INTENT(INOUT)      :: matching_string
      INTEGER, INTENT(IN)                                :: bonus

      CHARACTER(LEN=*), PARAMETER :: routineN = 'keyword_typo_match', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=LEN(matching_string(1)))             :: line
      INTEGER                                            :: i, imatch, imax, irank, j, k

      CPASSERT(ASSOCIATED(keyword))
      CPASSERT(keyword%ref_count > 0)

      DO i = 1, SIZE(keyword%names)
         imatch = typo_match(TRIM(keyword%names(i)), TRIM(unknown_string))
         IF (imatch > 0) THEN
            imatch = imatch + bonus
            WRITE (line, '(T2,A)') " keyword "//TRIM(keyword%names(i))//" in section "//TRIM(location_string)
            imax = SIZE(matching_rank, 1)
            irank = imax + 1
            DO k = imax, 1, -1
               IF (imatch > matching_rank(k)) irank = k
            ENDDO
            IF (irank <= imax) THEN
               matching_rank(irank + 1:imax) = matching_rank(irank:imax - 1)
               matching_string(irank + 1:imax) = matching_string(irank:imax - 1)
               matching_rank(irank) = imatch
               matching_string(irank) = line
            ENDIF
         END IF

         IF (keyword%type_of_var == enum_t) THEN
            DO j = 1, SIZE(keyword%enum%c_vals)
               imatch = typo_match(TRIM(keyword%enum%c_vals(j)), TRIM(unknown_string))
               IF (imatch > 0) THEN
                  imatch = imatch + bonus
                  WRITE (line, '(T2,A)') " enum "//TRIM(keyword%enum%c_vals(j))// &
                     " in section "//TRIM(location_string)// &
                     " for keyword "//TRIM(keyword%names(i))
                  imax = SIZE(matching_rank, 1)
                  irank = imax + 1
                  DO k = imax, 1, -1
                     IF (imatch > matching_rank(k)) irank = k
                  ENDDO
                  IF (irank <= imax) THEN
                     matching_rank(irank + 1:imax) = matching_rank(irank:imax - 1)
                     matching_string(irank + 1:imax) = matching_string(irank:imax - 1)
                     matching_rank(irank) = imatch
                     matching_string(irank) = line
                  ENDIF
               END IF
            END DO
         END IF
      ENDDO

   END SUBROUTINE keyword_typo_match

END MODULE input_keyword_types
