MODULE method_specifications

! Purpose: Calculation method specification.

! History: - Creation (13.06.2000, Matthias Krack)

! *****************************************************************************

! Index:

! SUBROUTINE read_method_specifications(start_section,end_section)
! SUBROUTINE write_method_specifications(lunit)

! *****************************************************************************

  USE functionals, ONLY: xc_type

  IMPLICIT NONE

  PRIVATE

  TYPE(xc_type)    :: xc_info
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
            periodic,&
            xc_info

! *** Public subroutines ***

  PUBLIC :: read_method_specifications,&
            write_method_specifications,&
            get_maxder

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
                                test_object,&
                                start_input_session
    USE string_utilities, ONLY: uppercase
    USE termination,      ONLY: stop_program

    TYPE(global_environment_type), INTENT(IN) :: globenv
    CHARACTER(LEN=*), INTENT(IN)              :: start_section,end_section

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE read_method_specifications (MODULE method_specifications)"

!   *** Local variables ***

    CHARACTER(LEN=60) :: keyword,name,test_result
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
    xc_info%x_functional = "Pade"
    xc_info%c_functional = "Pade"
    xc_info%gradient_functional = .FALSE.

    CALL start_input_session(globenv%input_file_name,globenv)

!   *** Search for the input section with the cell parameters ***

    CALL search(start_section,.TRUE.,found)

    IF (found) THEN

      CALL test_object(test_result)

      IF (test_result /= "end of line") THEN
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
                            "Invalid method keyword <"//TRIM(method)//&
                            "> specified in the input section <"//&
                            TRIM(start_section)//">")
        END SELECT
      END IF

      DO

        CALL read_object(keyword,skip_lines=1)

        CALL uppercase(keyword)

        SELECT CASE (TRIM(keyword))
        CASE ("CALCULATE_FORCES","FORCES","FORCE")
          maxder = 1
        CASE ("NO_PERIODICITY","NOPER","NON_PERIODIC","NONPER")
          periodic = .FALSE.
        CASE ("EXCHANGE-FUNCTIONAL","X-FUNCTIONAL","X-FUN")
          CALL read_object(name)
          CALL uppercase(name)
          SELECT CASE (name)
          CASE ("BECKE88","B88")
            xc_info%x_functional = "Becke88"
            xc_info%gradient_functional = .TRUE.
          CASE ("NONE")
            xc_info%x_functional = "None"
          CASE ("PBE")
            xc_info%x_functional = "PBE"
            xc_info%gradient_functional = .TRUE.
          CASE ("PERDEW86","P86")
            xc_info%x_functional = "Perdew86"
            xc_info%gradient_functional = .TRUE.
          CASE ("SLATER","XALPHA")
            xc_info%x_functional = "Slater"
          CASE ("VWN")
            xc_info%x_functional = "VWN"
          CASE DEFAULT
            CALL stop_program(routine,&
                              "Invalid functional name <"//TRIM(name)//&
                              "> specified for the exchange "//&
                              " functional in the input section <"//&
                              TRIM(start_section)//">")
          END SELECT
        CASE ("CORRELATION-FUNCTIONAL","C-FUNCTIONAL","C-FUN")
          CALL read_object(name)
          CALL uppercase(name)
          SELECT CASE (name)
          CASE ("LEE-YANG-PARR","LYP")
            xc_info%c_functional = "LYP"
            xc_info%gradient_functional = .TRUE.
          CASE ("NONE")
            xc_info%c_functional = "None"
          CASE ("PBE")
            xc_info%c_functional = "PBE"
            xc_info%gradient_functional = .TRUE.
          CASE ("PERDEW86","P86")
            xc_info%c_functional = "Perdew86"
            xc_info%gradient_functional = .TRUE.
          CASE ("VWN")
            xc_info%c_functional = "VWN"
          CASE DEFAULT
            CALL stop_program(routine,&
                              "Invalid functional name <"//TRIM(name)//&
                              "> specified for the correlation "//&
                              " functional in the input section <"//&
                              TRIM(start_section)//">")
          END SELECT
        CASE ("EXCHANGE-CORRELATION-FUNCTIONAL","XC-FUNCTIONAL","FUNCTIONAL")
          CALL read_object(name)
          CALL uppercase(name)
          SELECT CASE (name)
          CASE ("BLYP")
            xc_info%x_functional = "Becke88"
            xc_info%c_functional = "LYP"
            xc_info%gradient_functional = .TRUE.
          CASE ("BP")
            xc_info%x_functional = "Becke88"
            xc_info%c_functional = "Perdew86"
            xc_info%gradient_functional = .TRUE.
          CASE ("HCTH/93","HCTH-93","HCTH93")
            xc_info%x_functional = "HCTH/93"
            xc_info%c_functional = "HCTH/93"
            xc_info%gradient_functional = .TRUE.
          CASE ("HCTH/120","HCTH-120","HCTH120")
            xc_info%x_functional = "HCTH/120"
            xc_info%c_functional = "HCTH/120"
            xc_info%gradient_functional = .TRUE.
          CASE ("HCTH/147","HCTH-147","HCTH147")
            xc_info%x_functional = "HCTH/147"
            xc_info%c_functional = "HCTH/147"
            xc_info%gradient_functional = .TRUE.
          CASE ("HCTH/407","HCTH-407","HCTH407")
            xc_info%x_functional = "HCTH/407"
            xc_info%c_functional = "HCTH/407"
            xc_info%gradient_functional = .TRUE.
          CASE ("NONE")
            xc_info%x_functional = "None"
            xc_info%c_functional = "None"
          CASE ("PADE")
            xc_info%x_functional = "Pade"
            xc_info%c_functional = "Pade"
          CASE ("PBE")
            xc_info%x_functional = "PBE"
            xc_info%c_functional = "PBE"
            xc_info%gradient_functional = .TRUE.
          CASE ("PERDEW86","P86")
            xc_info%x_functional = "Perdew86"
            xc_info%c_functional = "Perdew86"
            xc_info%gradient_functional = .TRUE.
          CASE ("SLATER","XALPHA")
            xc_info%x_functional = "Slater"
            xc_info%c_functional = "None"
          CASE ("VWN")
            xc_info%x_functional = "VWN"
            xc_info%c_functional = "VWN"
          CASE DEFAULT
            CALL stop_program(routine,&
                              "Invalid functional name <"//TRIM(name)//&
                              "> specified for the exchange-correlation "//&
                              " functional in the input section <"//&
                              TRIM(start_section)//">")
          END SELECT
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

    WRITE (lunit, "(/,T3,A,/,T3,A)")&
      "Exchange functional:    "//TRIM(xc_info%x_functional),&
      "Correlation functional: "//TRIM(xc_info%c_functional)

  END SUBROUTINE write_method_specifications

! *****************************************************************************
  SUBROUTINE get_maxder(maxder_index)
    IMPLICIT NONE
    INTEGER, INTENT(in), OPTIONAL :: maxder_index  ! maxder
    IF (PRESENT(maxder_index)) maxder = maxder_index
  END SUBROUTINE get_maxder
! *****************************************************************************

END MODULE method_specifications
