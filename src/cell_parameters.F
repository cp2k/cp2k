!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/cell_parameters [1.0] *
!!
!!   NAME
!!     cell_parameters
!!
!!   FUNCTION
!!     Definition of the periodic cell parameters.
!!
!!   AUTHOR
!!     Matthias Krack (12.05.2000)
!!
!!   MODIFICATION HISTORY
!!     1) Modified to call init_cell (Joost VandeVondele, 22.02.02)
!!
!!   SOURCE
!******************************************************************************

MODULE cell_parameters

! *****************************************************************************

! Index:

! SUBROUTINE read_cell_parameters(start_section,end_section)
! SUBROUTINE write_cell_parameters(lunit)

! FUNCTION pbc(r) RESULT(r_pbc)

! *****************************************************************************

  USE kinds, ONLY: wp => dp

  USE simulation_cell, ONLY: cell_type,init_cell,get_hinv

  IMPLICIT NONE

  PRIVATE

  TYPE(cell_type)  :: cell
  CHARACTER(LEN=8) :: unit_of_length_name = "BOHR"
  REAL(wp)         :: cutoff,omega = 0.0_wp,&
                      unit_of_length = 1.0_wp
  LOGICAL          :: scaled_coordinates = .FALSE.

  REAL(wp), DIMENSION(3) :: abc = (/0.0_wp,0.0_wp,0.0_wp/)

! *** Public variables ***

  PUBLIC :: abc,&
            cell,&
            cutoff,&
            omega,&
            scaled_coordinates,&
            unit_of_length,&
            unit_of_length_name

! *** Public subroutines ***

  PUBLIC :: pbc,&
            read_cell_parameters,&
            write_cell_parameters

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE read_cell_parameters(start_section,end_section,globenv)

!   Purpose: Read the periodic cell parameters.

!   History: - Creation (19.05.2000, Matthias Krack)
!   modified to call the init_cell routine (22.02.2002, Joost VandeVondele)

!   ***************************************************************************

    USE global_types,          ONLY: global_environment_type
    USE input_utilities,       ONLY: finish_input_session,&
                                     read_object,&
                                     search,&
                                     start_input_session
    USE method_specifications, ONLY: gpw
    USE physcon,               ONLY: bohr
    USE string_utilities,      ONLY: uppercase
    USE termination,           ONLY: stop_program

    TYPE(global_environment_type), INTENT(IN) :: globenv
    CHARACTER(LEN=*), INTENT(IN)              :: start_section,end_section

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE read_cell_parameters (MODULE cell_parameters)"

!   *** Local variables ***

    CHARACTER(LEN=40) :: keyword,string
    LOGICAL           :: found
    real(wp)          :: hmat(3,3)
    integer           :: periodic(3)
    CHARACTER ( LEN = 80 ) :: poisson_solver

!   ---------------------------------------------------------------------------

!   *** Load the default values ***

    unit_of_length_name = "BOHR"
    unit_of_length = 1.0_wp
    scaled_coordinates = .FALSE.
    abc = (/0.0_wp,0.0_wp,0.0_wp/)
    periodic = (/1,1,1/)
    poisson_solver = "ANALYTIC"

    CALL start_input_session(globenv%input_file_name,globenv)

!   *** Search for the input section with the cell parameters ***

    CALL search(start_section,.TRUE.,found)

    IF (found) THEN

      DO

        CALL read_object(keyword,skip_lines=1)

        CALL uppercase(keyword)

        SELECT CASE (TRIM(keyword))
        CASE ("ISOLATED")
          periodic=(/0,0,0/)
        CASE ("ABC")
          CALL read_object(abc(1))
          CALL read_object(abc(2))
          CALL read_object(abc(3))
        CASE ("CUTOFF")
          CALL read_object(cutoff)
          cutoff = 0.5_wp*cutoff
        CASE ("UNIT")
          CALL read_object(string)
          SELECT CASE (TRIM(string))
          CASE ("ANGSTROM")
            unit_of_length = bohr
            unit_of_length_name = "ANGSTROM"
            scaled_coordinates = .FALSE.
          CASE ("BOHR")
            unit_of_length = 1.0_wp
            unit_of_length_name = "BOHR"
            scaled_coordinates = .FALSE.
          CASE ("SCALED_ANGSTROM")
            unit_of_length = bohr
            unit_of_length_name = "ANGSTROM"
            scaled_coordinates = .TRUE.
          CASE ("SCALED_BOHR")
            unit_of_length = 1.0_wp
            unit_of_length_name = "BOHR"
            scaled_coordinates = .TRUE.
          CASE DEFAULT
            CALL stop_program(routine,&
                              "Invalid name for the unit of length <"//&
                              TRIM(string)//&
                              "> found in the input section <"//&
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

      abc(:) = abc(:)*unit_of_length
      omega = PRODUCT(abc(:))

    ELSE

!     *** Stop program, if the end of file is reached ***

      IF (gpw) THEN
        CALL stop_program(routine,&
                          "No input section <"//TRIM(start_section)//&
                          "> found in the input file <"//&
                          TRIM(globenv%input_file_name)//">")
      END IF

    END IF

    CALL finish_input_session

!   *** Check the cell parameters ***

    IF (gpw.OR.scaled_coordinates) THEN
      IF ((abc(1) <= 0.0_wp).OR.&
          (abc(2) <= 0.0_wp).OR.&
          (abc(3) <= 0.0_wp)) THEN
        CALL stop_program(routine,"Invalid cell dimensions were specified")
      END IF
    END IF

    hmat(:,:) = 0.0_wp
    hmat(1,1) = abc(1)
    hmat(2,2) = abc(2)
    hmat(3,3) = abc(3)

    CALL init_cell(cell,hmat,periodic,poisson_solver)

  END SUBROUTINE read_cell_parameters

! *****************************************************************************

  SUBROUTINE write_cell_parameters(lunit)

!   Purpose: Write the periodic cell parameters to the logical unit number
!            "lunit".

!   History: - Creation (02.06.2000, Matthias Krack)

!   ***************************************************************************

    INTEGER, INTENT(IN) :: lunit

!   ---------------------------------------------------------------------------

    WRITE (UNIT=lunit,FMT="(/,/,T2,A,/)") "CELL PARAMETERS"

    WRITE (UNIT=lunit,FMT="(T3,A,T26,3F10.3,/,/,T3,A,T31,F15.3,A)")&
      "Cell size in "//TRIM(unit_of_length_name)//":",&
      abc(:)/unit_of_length,&
      "Cell volume in "//TRIM(unit_of_length_name)//"**3:",&
      omega/unit_of_length**3

  END SUBROUTINE write_cell_parameters

! *****************************************************************************

  FUNCTION pbc(r) RESULT(r_pbc)

!   Purpose: Apply the periodic boundary condition to the position vector r.

!   History: - Creation (22.05.2000, Matthias Krack)

!   ***************************************************************************

    USE method_specifications, ONLY: periodic

    REAL(wp), DIMENSION(3) :: r

    REAL(wp), DIMENSION(3) :: r_pbc

!   ---------------------------------------------------------------------------

    IF (periodic) THEN
      r_pbc(:) = r(:) - abc(:)*ANINT(r(:)/abc(:))
    ELSE
      r_pbc(:) = r(:)
    END IF

  END FUNCTION pbc

! *****************************************************************************

END MODULE cell_parameters
