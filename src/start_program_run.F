!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 1999  MPI fuer Festkoerperforschung, Stuttgart              !
!-----------------------------------------------------------------------------!
!!****** cp2k/start_program_run [1.0] *
!!
!!   NAME
!!     start_program_run
!!
!!   FUNCTION
!!      Perform the program startup.
!!
!!   AUTHOR
!!     Matthias Krack (13.07.2000)
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

MODULE start_program_run

! *****************************************************************************

! Index:

! SUBROUTINE start_qs_program_run(globenv)
! SUBROUTINE read_io_parameters(start_section,end_section,globenv)
! SUBROUTINE read_print_keys(start_section,end_section,globenv)
! SUBROUTINE write_program_banner(lunit)

! *****************************************************************************

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: start_qs_program_run

! *****************************************************************************

CONTAINS

! *****************************************************************************

  SUBROUTINE start_qs_program_run(globenv)

!   Purpose: Start a Quickstep program run.

!   History: - Creation (13.07.2000, Matthias Krack)

!   ***************************************************************************

    USE global_types,       ONLY: global_environment_type
    USE input_utilities,    ONLY: get_special_input_character

    TYPE(global_environment_type), INTENT(INOUT) :: globenv

!   *** Local variables ***

    CHARACTER(LEN=4) :: end_section
    CHARACTER(LEN=1) :: section_character

!   ---------------------------------------------------------------------------

    CALL get_special_input_character(section=section_character)

    end_section = section_character//"END"

!   *** Read the input section with the input/output parameters ***

    CALL read_io_parameters(section_character//"IO",end_section,globenv)

!   *** Read the input section with the print keys ***

    CALL read_print_keys(section_character//"PRINT",end_section,globenv)

!   *** Print program banner (copyright and version number) ***

    IF (globenv%ionode.AND.globenv%print%program_banner) THEN
      CALL write_program_banner(globenv%scr)
    END IF

  END SUBROUTINE start_qs_program_run

! *****************************************************************************

  SUBROUTINE read_io_parameters(start_section,end_section,globenv)

!   Purpose: Read the input/output parameters (file names).

!   History: - Creation (13.07.2000, Matthias Krack)

!   ***************************************************************************

    USE global_types,     ONLY: global_environment_type
    USE input_utilities,  ONLY: finish_input_session,&
                                read_object,&
                                search,&
                                start_input_session
    USE string_utilities, ONLY: uppercase
    USE termination,      ONLY: stop_program

    TYPE(global_environment_type), INTENT(INOUT) :: globenv
    CHARACTER(LEN=*), INTENT(IN)                 :: start_section,end_section

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE read_io_parameters (MODULE start_program_run)"

!   *** Local variables ***

    CHARACTER(LEN=40) :: keyword
    LOGICAL           :: found

!   ---------------------------------------------------------------------------

    CALL start_input_session(globenv%input_file_name,globenv)

!   *** Search for the requested input section ***

    CALL search(start_section,.TRUE.,found)

    IF (found) THEN

      DO

        CALL read_object(keyword,skip_lines=1)

        CALL uppercase(keyword)

        SELECT CASE (TRIM(keyword))
        CASE ("BASIS_SET_FILE_NAME","BASIS_SET_FILE","BASIS_SET","BSF")
          CALL read_object(globenv%basis_set_file_name)
        CASE ("POTENTIAL_FILE_NAME","POTENTIAL_FILE","POTENTIAL","PPF")
          CALL read_object(globenv%potential_file_name)
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

  END SUBROUTINE read_io_parameters

! *****************************************************************************

  SUBROUTINE read_print_keys(start_section,end_section,globenv)

!   Purpose: Read the print keys.

!   History: - Creation (13.07.2000, Matthias Krack)

!   ***************************************************************************

    USE global_types,     ONLY: global_environment_type
    USE input_utilities,  ONLY: finish_input_session,&
                                read_object,&
                                search,&
                                start_input_session,&
                                test_object
    USE print_keys,       ONLY: DEBUG,HIGH,init_print_keys,LOW,MEDIUM,SILENT
    USE string_utilities, ONLY: uppercase
    USE termination,      ONLY: stop_program

    TYPE(global_environment_type), INTENT(INOUT) :: globenv
    CHARACTER(LEN=*), INTENT(IN)                 :: start_section,end_section

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE read_print_keys (MODULE start_program_run)"

!   *** Local variables ***

    CHARACTER(LEN=40) :: keyword,test_result
    INTEGER           :: ipos
    LOGICAL           :: found,print_request

!   ---------------------------------------------------------------------------

    CALL start_input_session(globenv%input_file_name,globenv)

!   *** Search for the input section with the print keys ***

    CALL search(start_section,.TRUE.,found)

    IF (found) THEN

      CALL test_object(test_result)

      IF (test_result == "end of line") THEN
        globenv%print%level = 1
      ELSE IF (test_result == "integer") THEN
        CALL read_object(globenv%print%level)
      ELSE
        CALL read_object(keyword)
        CALL uppercase(keyword)
        SELECT CASE (TRIM(keyword))
        CASE ("SILENT")
          globenv%print%level = SILENT
        CASE ("LOW","DEFAULT")
          globenv%print%level = LOW
        CASE ("MEDIUM")
          globenv%print%level = MEDIUM
        CASE ("HIGH")
          globenv%print%level = HIGH
        CASE ("DEBUG","ALL","FULL")
          globenv%print%level = DEBUG
        CASE DEFAULT
          CALL stop_program(routine,&
                            "Invalid keyword <"//TRIM(keyword)//&
                            "> found in the input section <"//&
                            TRIM(start_section)//&
                            "> for the default print level")
        END SELECT
      END IF

      CALL init_print_keys(globenv%print,globenv%print%level)

      DO

        CALL read_object(keyword,skip_lines=1)

        CALL uppercase(keyword)

        IF (keyword(1:3) == "NO_") THEN
          ipos = 4
          print_request = .FALSE.
        ELSE
          ipos = 1
          print_request = .TRUE.
        END IF

        SELECT CASE (TRIM(keyword(ipos:)))
        CASE ("ANGLES")
          globenv%print%angles = print_request
        CASE ("ATOMIC_COORDINATES","COORDINATES","COORD")
          globenv%print%atomic_coordinates = print_request
        CASE ("BASIC_DATA_TYPES")
          globenv%print%basic_data_types = print_request
        CASE ("BASIS_SETS","BASIS_SET","BASIS")
          globenv%print%basis_set = print_request
        CASE ("CARTESIAN_OPERATOR_MATRICES","CARTESIAN_MATRICES")
          globenv%print%cartesian_om = print_request
        CASE ("CELL_PARAMETERS","CELL")
          globenv%print%cell_parameters = print_request
        CASE ("CHARGE_DENSITY_MATRIX","G_MATRIX")
          globenv%print%charge_density_matrix = print_request
        CASE ("CORE_HAMILTONIAN_MATRIX","H_MATRIX")
          globenv%print%core_hamiltonian_matrix = print_request
        CASE ("CORE_CHARGE_RADII","CORE_RADII")
          globenv%print%core_charge_radii = print_request
        CASE ("DENSITY_MATRIX","P_MATRIX")
          globenv%print%density_matrix = print_request
        CASE ("DIIS_INFORMATION","DIIS_INFO","DIIS")
          globenv%print%diis_information = print_request
        CASE ("DISTRIBUTION")
          globenv%print%distribution = print_request
        CASE ("EACH_SCF_STEP")
          globenv%print%each_scf_step = print_request
        CASE ("HARTREE_MATRIX")
          globenv%print%hartree_matrix = print_request
        CASE ("INTERACTION_THRESHOLDS","THRESHOLDS")
          globenv%print%interaction_thresholds = print_request
        CASE ("INTERATOMIC_DISTANCES","DISTANCES")
          globenv%print%interatomic_distances = print_request
        CASE ("KIND_INFORMATION","KIND_INFO","KINDS")
          globenv%print%basis_set = print_request
          globenv%print%potential = print_request
        CASE ("KIND_RADII")
          globenv%print%kind_radii = print_request
        CASE ("KINETIC_ENERGY_MATRIX","KINETIC_ENERGY","KINETIC","T_MATRIX")
          globenv%print%kinetic_energy_matrix = print_request
        CASE ("KOHN_SHAM_MATRIX","KS_MATRIX")
          globenv%print%kohn_sham_matrix = print_request
        CASE ("MEMORY")
          globenv%print%memory = print_request
        CASE ("METHOD_SPECIFICATIONS","METHOD")
          globenv%print%method_specifications = print_request
        CASE ("MO_EIGENVALUES")
          globenv%print%mo_eigenvalues = print_request
        CASE ("MO_EIGENVECTORS","MOS")
          globenv%print%mo_eigenvalues = print_request
          globenv%print%mo_eigenvectors = print_request
          globenv%print%mo_occupation_numbers = print_request
        CASE ("MO_OCCUPATION_NUMBERS")
          globenv%print%mo_eigenvalues = print_request
          globenv%print%mo_occupation_numbers = print_request
        CASE ("NEIGHBOR_LISTS")
          globenv%print%nl_overlap = print_request
          globenv%print%nl_verfc = print_request
          globenv%print%nl_vppl = print_request
          globenv%print%nl_vppnl = print_request
        CASE ("NL_OVERLAP")
          globenv%print%nl_overlap = print_request
        CASE ("NL_VERFC")
          globenv%print%nl_verfc = print_request
        CASE ("NL_VPPL")
          globenv%print%nl_vppl = print_request
        CASE ("NL_VPPNL")
          globenv%print%nl_vppnl = print_request
        CASE ("ORTHOGONALISATION_MATRIX","ORTHO_MATRIX")
          globenv%print%ortho_matrix = print_request
        CASE ("OVERLAP_MATRIX","OVERLAP","S_MATRIX")
          globenv%print%overlap_matrix = print_request
        CASE ("PGF_RADII")
          globenv%print%pgf_radii = print_request
        CASE ("PHYSICAL_CONSTANTS","PHYSCON")
          globenv%print%physcon = print_request
        CASE ("POTENTIALS","POTENTIAL")
          globenv%print%potential = print_request
        CASE ("PPL_RADII")
          globenv%print%ppl_radii = print_request
        CASE ("PPNL_RADII")
          globenv%print%ppnl_radii = print_request
        CASE ("PROGRAM_BANNER","BANNER")
          globenv%print%program_banner = print_request
        CASE ("PROGRAM_RUN_INFORMATION","PROGRAM_RUN_INFO","RUN_INFO")
          globenv%print%program_run_information = print_request
        CASE ("PROJECTORS")
          globenv%print%projectors = print_request
        CASE ("RADII")
          globenv%print%core_charge_radii = print_request
          globenv%print%kind_radii = print_request
          globenv%print%ppl_radii = print_request
          globenv%print%ppnl_radii = print_request
          globenv%print%set_radii = print_request
          globenv%print%pgf_radii = print_request
        CASE ("PW_GRID_INFORMATION","PW_GRID_INFO","GRID_INFO")
          globenv%print%pw_grid_information = print_request
        CASE ("SCF")
          globenv%print%scf = print_request
        CASE ("SCF_ENERGIES")
          globenv%print%scf = print_request
          globenv%print%scf_energies = print_request
        CASE ("SET_RADII")
          globenv%print%set_radii = print_request
        CASE ("SPHERICAL_HARMONICS")
          globenv%print%spherical_harmonics = print_request
        CASE ("TIMING_INFORMATION","TIMING")
          globenv%print%timing_information = print_request
        CASE ("TITLE")
          globenv%print%title = print_request
        CASE ("TOTAL_NUMBERS")
          globenv%print%total_numbers = print_request
        CASE ("TOTAL_DENSITIES")
          globenv%print%total_densities = print_request
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

  END SUBROUTINE read_print_keys

! *****************************************************************************

  SUBROUTINE write_program_banner(lunit)

!   Purpose: Write the program banner with the copyright and the program
!            version information to the logical unit number "lunit".

!   History: - Creation (15.06.2000, Matthias Krack)

!   ***************************************************************************

    INTEGER, INTENT(IN) :: lunit

!   ---------------------------------------------------------------------------

    WRITE (lunit,"(T2,A78)")&
      "******************************************************************************",&
      "******************************************************************************",&
      "**                                                                          **",&
      "**     #####                         ##              ##                     **",&
      "**    ##   ##            ##          ##              ##                     **",&
      "**   ##     ##                       ##            ######                   **",&
      "**   ##     ##  ##   ##  ##   #####  ##  ##   ####   ##    #####    #####   **",&
      "**   ##     ##  ##   ##  ##  ##      ## ##   ##      ##   ##   ##  ##   ##  **",&
      "**   ##  ## ##  ##   ##  ##  ##      ####     ###    ##   ######   ######   **",&
      "**    ##  ###   ##   ##  ##  ##      ## ##      ##   ##   ##       ##       **",&
      "**     #######   #####   ##   #####  ##  ##  ####    ##    #####   ##       **",&
      "**           ##                                                    ##       **",&
      "**                                                                          **",&
      "**                                                ... make the atoms dance  **",&
      "**                                                                          **",&
      "**                         Version 3.0 (May 2001)                           **",&
      "**                                                                          **",&
      "**     Copyright (C) by MPI fuer Festkoerperforschung, Stuttgart (2000)     **",&
      "**                                                                          **",&
      "******************************************************************************",&
      "******************************************************************************"

  END SUBROUTINE write_program_banner

! *****************************************************************************

END MODULE start_program_run
