MODULE method_specifications

! Purpose: Calculation method specification.

! History: - Creation (13.06.2000, Matthias Krack)

! *****************************************************************************

! Index:

! SUBROUTINE read_method_specifications(start_section,end_section)
! SUBROUTINE write_method_specifications(lunit)

! *****************************************************************************

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=8) :: method = "GAPW"
  INTEGER          :: maxder = 0
  LOGICAL          :: allchem = .FALSE.,&
                      gapw = .TRUE.,&
                      gpw = .TRUE.,&
                      oks = .FALSE.,&
                      periodic = .TRUE.

! *** Public variables ***

  PUBLIC :: allchem,&
            gapw,&
            gpw,&
            maxder,&
            method,&
            oks,&
            periodic

! *** Public subroutines ***

  PUBLIC :: read_method_specifications,&
            write_method_specifications

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE read_method_specifications(start_section,end_section,globenv)

!   Purpose: Read the method specifications.

!   History: - Creation (13.06.2000, Matthias Krack)

!   ***************************************************************************

    USE global_types,     ONLY: global_environment_type
    USE input_utilities,  ONLY: finish_input_session,&
                                read_object,&
                                search,&
                                start_input_session
    USE string_utilities, ONLY: uppercase
    USE termination,      ONLY: stop_program

    TYPE(global_environment_type), INTENT(IN) :: globenv
    CHARACTER(LEN=*), INTENT(IN)              :: start_section,end_section

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE read_method_specifications (MODULE method_specifications)"

!   *** Local variables ***

    CHARACTER(LEN=40) :: keyword
    LOGICAL           :: found

!   ---------------------------------------------------------------------------

!   *** Load the default values ***

    allchem = .FALSE.
    gapw = .TRUE.
    gpw = .TRUE.
    maxder = 0
    method = "GAPW"
    oks = .FALSE.
    periodic = .TRUE.

    CALL start_input_session(globenv%input_file_name,globenv)

!   *** Search for the input section with the cell parameters ***

    CALL search(start_section,.TRUE.,found)

    IF (found) THEN

      DO

        CALL read_object(keyword,skip_lines=1)

        CALL uppercase(keyword)

        SELECT CASE (TRIM(keyword))
        CASE ("METHOD")
          CALL read_object(method)
          CALL uppercase(method)
          SELECT CASE (TRIM(method))
          CASE ("ALLCHEM")
            allchem = .TRUE.
            gpw = .FALSE.
            gapw = .FALSE.
            periodic = .FALSE.
          CASE ("GAPW")
            allchem = .FALSE.
            gpw = .TRUE.
            gapw = .TRUE.
          CASE ("GPW")
            allchem = .FALSE.
            gpw = .TRUE.
            gapw = .FALSE.
          CASE DEFAULT
            CALL stop_program(routine,&
                              "Invalid keyword <"//TRIM(method)//&
                              "> specified for METHOD in the input "//&
                              "section <"//TRIM(start_section)//">")
          END SELECT
        CASE ("CALCULATE_FORCES","FORCES","FORCE")
          maxder = 1
        CASE ("NO_PERIODICITY","NOPER","NON_PERIODIC","NONPER")
          periodic = .FALSE.
        CASE DEFAULT
          IF (keyword == end_section) THEN
            EXIT
          ELSE
            CALL stop_program(routine,&
                              "Invalid keyword <"//TRIM(keyword)//&
                              "> found in the input section <"//&
                              TRIM(start_section)//">")
          END IF
        END SELECT

      END DO

    END IF

    CALL finish_input_session

  END SUBROUTINE read_method_specifications

! *****************************************************************************

  SUBROUTINE write_method_specifications(lunit)

!   Purpose: Write the method specifications to the logical unit number
!            "lunit".

!   ***************************************************************************

    INTEGER, INTENT(IN) :: lunit

!   ---------------------------------------------------------------------------

    WRITE (lunit,"(/,/,T2,A)") "METHOD SPECIFICATIONS"

    WRITE (lunit, "(/,T3,A,/)") "Method: "//TRIM(method)

    IF (maxder > 0) THEN
      WRITE (lunit,"(T3,A)") "Force calculation:    YES"
    ELSE
      WRITE (lunit,"(T3,A)") "Force calculation:    NO"
    END IF

    IF (periodic) THEN
      WRITE (lunit,"(T3,A)") "Periodic calculation: YES"
    ELSE
      WRITE (lunit,"(T3,A)") "Periodic calculation: NO"
    END IF

  END SUBROUTINE write_method_specifications

! *****************************************************************************

END MODULE method_specifications
