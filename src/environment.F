!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2002  CP2K developers group
!-----------------------------------------------------------------------------!
!!****** cp2k/environment [1.0] *
!!
!!   NAME
!!     environment
!!
!!   FUNCTION
!!     Sets up and terminates the global environment variables
!!
!!   AUTHOR
!!     JGH,MK
!!
!!   MODIFICATION HISTORY
!!     - init_mathcon is now called (27.02.2001, MK)
!!     - Merged with Quickstep MODULE start_program_run (17.01.2002, MK)
!!     - Compile information added (16.01.2002,MK)
!!
!!   SOURCE
!******************************************************************************

MODULE environment

  USE global_types,     ONLY: global_environment_type,&
                              init_print_keys,&
                              DEBUG,HIGH,LOW,MEDIUM,SILENT
  USE machine,          ONLY: m_datum,&
                              m_hostnm,&
                              m_getcwd,&
                              m_getlog,&
                              m_getuid,&
                              m_getpid,&
                              m_getarg
  USE mathconstants,    ONLY: init_mathcon,zero
  USE physcon,          ONLY: init_physcon
  USE qs_parser,        ONLY: finish_parser,&
                              read_object,&
                              search_object,&
                              start_parser,&
                              stop_parser,&
                              test_object
  USE string_utilities, ONLY: uppercase
  USE termination,      ONLY: set_error_unit,&
                              stop_program
  USE timings,          ONLY: timeprint,&
                              timeset,&
                              timestop
  USE util,             ONLY: ran2

  IMPLICIT NONE

  PRIVATE

#if defined(__COMPILE_ARCH)
  CHARACTER(LEN=*), PARAMETER :: compile_arch = __COMPILE_ARCH
#endif
#if defined(__COMPILE_DATE)
  CHARACTER(LEN=*), PARAMETER :: compile_date = __COMPILE_DATE
#endif
#if defined(__COMPILE_HOST)
  CHARACTER(LEN=*), PARAMETER :: compile_host = __COMPILE_HOST
#endif
  CHARACTER(LEN=60) :: curdir = " " !! current directory name
  CHARACTER(LEN=32) :: hname = " "  !! name of the computer
  CHARACTER(LEN=10) :: user = " "   !! user name
  INTEGER           :: handle
  INTEGER           :: my_pid = 0   !! process ID
  INTEGER           :: my_uid = 0   !! user ID

! *** Public subroutines ***

  PUBLIC :: initialisation,&
            trailer

!!*****
!******************************************************************************

CONTAINS

!******************************************************************************
!!****** environment/initialisation [1.0] *
!!
!!   NAME
!!     initialisation
!!
!!   FUNCTION
!!     Initialises global environment variables
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     JGH (28.11.2001) : default for pp_library_path
!!     - print keys added (17.01.2002, MK)
!!
!!   SOURCE
!******************************************************************************

  SUBROUTINE initialisation(globenv)

    TYPE(global_environment_type), INTENT(INOUT) :: globenv

!   *** Local variables ***

    CHARACTER(LEN=26) :: datx
    CHARACTER(LEN=4)  :: end_section
    CHARACTER(LEN=1)  :: section_character
    INTEGER           :: i1,i2,initial_random,l,output_unit
    LOGICAL           :: ionode

!   ---------------------------------------------------------------------------

!   *** Set the default logical output unit number ***

    globenv%scr = 6
    CALL set_error_unit(globenv%scr)

    globenv%pp_library_path     = "./"

    output_unit = globenv%scr
    ionode = globenv%ionode

    IF (ionode) CALL set_environment(globenv)

!   *** Initialise timing ***

    CALL timeset ("CP2K","I"," ",handle)

!   *** Initialise random number generator ***

    globenv%idum = -1
    initial_random = ran2(globenv%idum)

!   *** Initialise mathematical constants ***

    CALL init_mathcon()

!   *** Initialise physical constants ***

    CALL init_physcon()

!   *** Read the input section with the print keys ***

    CALL read_print_keys(globenv)

    IF (ionode.AND.globenv%print%program_banner) THEN
      CALL m_datum(datx)
      WRITE (UNIT=output_unit,&
             FMT="(A,T55,A26,/,A,T49,A32,/,A,T71,A10,/,A,T71,I10)")&
        '  **** **** ******  **  PROGRAM STARTED AT ', ADJUSTR ( datx ), &
        ' ***** ** ***  *** **   PROGRAM STARTED ON ', ADJUSTR ( hname ), &
        ' **    ****   ******    PROGRAM STARTED BY ', ADJUSTR ( user ), &
        ' ***** **    ** ** **   PROGRAM PROCESS ID ', my_pid
      l = LEN_TRIM(curdir)
      IF ( l <= 38 ) THEN
        WRITE (UNIT=output_unit,FMT="(A,T43,A38)") &
          '  **** **  *******  **  PROGRAM STARTED IN ', &
          ADJUSTR(curdir(1:l))
      ELSE
        WRITE (UNIT=output_unit,FMT="(A)") &
          '  **** **  *******  **  PROGRAM STARTED IN '
        WRITE (UNIT=output_unit,FMT="(T21,A60)") ADJUSTR(curdir)
      END IF
      WRITE (UNIT=output_unit,FMT=*)
#if defined(__COMPILE_DATE)
      WRITE (UNIT=*,FMT="(T2,A,T31,A50)")&
        "CP2K| Program compiled at",&
        ADJUSTR(compile_date(1:MIN(50,LEN(compile_date))))
#endif
#if defined(__COMPILE_HOST)
      WRITE (UNIT=*,FMT="(T2,A,T31,A50)")&
        "CP2K| Program compiled on",&
        ADJUSTR(compile_host(1:MIN(50,LEN(compile_host))))
#endif
#if defined(__COMPILE_ARCH)
      WRITE (UNIT=*,FMT="(T2,A,T31,A50)")&
        "CP2K| Program compiled for",&
        ADJUSTR(compile_arch(1:MIN(50,LEN(compile_arch))))
#endif
    END IF

!   *** Read the input section with the input/output parameters ***

    CALL read_io_parameters(globenv)

  END SUBROUTINE initialisation

! *****************************************************************************

  SUBROUTINE read_io_parameters(globenv)

!   Purpose: Read the input/output parameters (file names).

!   History: - Creation (13.07.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(INOUT) :: globenv

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE read_io_parameters (MODULE environment)"

!   *** Local variables ***

    CHARACTER(LEN=60) :: keyword
    LOGICAL           :: section_found

!   ---------------------------------------------------------------------------

    globenv%basis_set_file_name = "BASIS_SET"
    globenv%potential_file_name = "POTENTIAL"
    globenv%restart_file_name = "RESTART"

    CALL start_parser(file_name=globenv%input_file_name,&
                      globenv=globenv,&
                      start_section_label="IO",&
                      section_found=section_found)

    DO WHILE (section_found.AND.(test_object(newline=.TRUE.) /= "EOS"))
      CALL read_object(keyword,lower_to_upper=.TRUE.,newline=.TRUE.)
      SELECT CASE (TRIM(keyword))
      CASE ("BASIS_SET_FILE_NAME","BASIS_SET_FILE","BASIS_SET","BSF")
        CALL read_object(globenv%basis_set_file_name)
      CASE ("POTENTIAL_FILE_NAME","POTENTIAL_FILE","POTENTIAL","PPF")
        CALL read_object(globenv%potential_file_name)
      CASE ("RESTART_FILE_NAME","RESTART_FILE")
        CALL read_object(globenv%restart_file_name)
      END SELECT
    END DO

    CALL finish_parser()

  END SUBROUTINE read_io_parameters

! *****************************************************************************

  SUBROUTINE read_print_keys(globenv)

!   Purpose: Read the print keys.

!   History: - Creation (13.07.2000, Matthias Krack)

!   ***************************************************************************

    TYPE(global_environment_type), INTENT(INOUT) :: globenv

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routine =&
      "SUBROUTINE read_print_keys (MODULE environment)"

!   *** Local variables ***

    CHARACTER(LEN=60) :: keyword
    CHARACTER(LEN=3)  :: test_result
    INTEGER           :: ipos1,ipos2
    LOGICAL           :: print_request,section_found

!   ---------------------------------------------------------------------------

!   *** Initialise print keys (the default print level is LOW) ***

    globenv%print%level = LOW

!   *** Modify the default precision of flexible format output ***

    globenv%print%precision = 0

    CALL start_parser(file_name=globenv%input_file_name,&
                      globenv=globenv,&
                      start_section_label="PRINT",&
                      section_found=section_found)

!   *** Check for print level specification ***

    IF (section_found) THEN
      test_result = test_object()
      IF (test_result == "EOL") THEN
        globenv%print%level = LOW
      ELSE IF (test_result == "INT") THEN
        CALL read_object(globenv%print%level)
      ELSE IF (test_result == "STR") THEN
        CALL read_object(keyword,lower_to_upper=.TRUE.)
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
          CALL stop_parser(routine,&
                           "Invalid keyword <"//TRIM(keyword)//&
                           "> found for the default print level")
        END SELECT
      END IF
    END IF

    CALL init_print_keys(globenv%print,globenv%print%level)

    DO WHILE (section_found.AND.(test_object(newline=.TRUE.) /= "EOS"))

      CALL read_object(keyword,lower_to_upper=.TRUE.,newline=.TRUE.)

      IF (keyword(1:3) == "NO_") THEN
        ipos1= 4
        print_request = .FALSE.
      ELSE
        ipos1= 1
        print_request = .TRUE.
      END IF

      ipos2 = LEN_TRIM(keyword)

      SELECT CASE (keyword(ipos1:ipos2))
      CASE ("ANGLES")
        globenv%print%angles = print_request
      CASE ("ATOMIC_COORDINATES","COORDINATES","COORD")
        globenv%print%atomic_coordinates = print_request
      CASE ("BASIC_DATA_TYPES")
        globenv%print%basic_data_types = print_request
      CASE ("BASIS_SETS","BASIS_SET","BASIS")
        globenv%print%basis_set = print_request
      CASE ("BLACS_INFORMATION","BLACS_INFO")
        globenv%print%blacs_info = print_request
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
      CASE ("CPUTIME")
        globenv%print%cputime = print_request
      CASE ("DENSITY_MATRIX","P_MATRIX")
        globenv%print%density_matrix = print_request
      CASE ("DENSITY_MATRIX_MAGNITUDE")
        globenv%print%density_matrix_magnitude = print_request
      CASE ("OVERLAP_MATRIX_MAGNITUDE")
        globenv%print%overlap_matrix_magnitude = print_request
      CASE ("DERIVATIVES")
        globenv%print%derivatives = print_request
      CASE ("DIIS_INFORMATION","DIIS_INFO","DIIS")
        globenv%print%diis_information = print_request
      CASE ("DISTRIBUTION")
        globenv%print%distribution = print_request
      CASE ("EACH_SCF_STEP")
        globenv%print%each_scf_step = print_request
      CASE ("FORCES")
        globenv%print%forces = print_request
      CASE ("HARTREE_MATRIX")
        globenv%print%hartree_matrix = print_request
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
        globenv%print%sab_orb_neighbor_lists = print_request
        globenv%print%sac_ppl_neighbor_lists = print_request
        globenv%print%sca_ppnl_neighbor_lists = print_request
        globenv%print%scb_ppnl_neighbor_lists = print_request
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
      CASE ("PRECISION")
        CALL read_object(globenv%print%precision)
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
      CASE ("SAB_ORB_NEIGHBOR_LIST")
        globenv%print%sab_orb_neighbor_lists = print_request
      CASE ("SAC_PPL_NEIGHBOR_LIST")
        globenv%print%sac_ppl_neighbor_lists = print_request
      CASE ("SCA_PPNL_NEIGHBOR_LIST")
        globenv%print%sca_ppnl_neighbor_lists = print_request
      CASE ("SCB_PPNL_NEIGHBOR_LIST")
        globenv%print%scb_ppnl_neighbor_lists = print_request
      CASE ("SCF")
        globenv%print%scf = print_request
      CASE ("SCF_ENERGIES")
        globenv%print%scf = print_request
        globenv%print%scf_energies = print_request
      CASE ("SET_RADII")
        globenv%print%set_radii = print_request
      CASE ("SPHERICAL_HARMONICS")
        globenv%print%spherical_harmonics = print_request
      CASE ("TITLE")
        globenv%print%title = print_request
      CASE ("TOTAL_NUMBERS")
        globenv%print%total_numbers = print_request
      CASE ("TOTAL_DENSITIES")
        globenv%print%total_densities = print_request
      CASE ("W_MATRIX","WMATRIX")
        globenv%print%w_matrix = print_request
      CASE ("CUBE","CUBES")
        globenv%print%e_density_cube = print_request
        globenv%print%v_hartree_cube = print_request
        globenv%print%lumos = print_request
        globenv%print%homos = print_request
      CASE ("E_DENSITY_CUBE")
        globenv%print%e_density_cube = print_request
      CASE ("V_HARTREE_CUBE")
        globenv%print%v_hartree_cube = print_request
      CASE ("LUMO","LUMOS")
        globenv%print%lumos = print_request
        IF (test_object().eq."INT") THEN
           CALL read_object(globenv%print%nlumos)
        ENDIF
        if (globenv%print%nlumos.le.0) globenv%print%nlumos = 4
      CASE ("HOMO","HOMOS")
        globenv%print%homos = print_request
        IF (test_object().eq."INT") THEN
           CALL read_object(globenv%print%nhomos)
        ENDIF
        if (globenv%print%nhomos.le.0) globenv%print%nhomos = 4
      CASE ("TIMECUMUL")
        globenv%print%timecumul = print_request
      END SELECT

    END DO

    CALL finish_parser()

  END SUBROUTINE read_print_keys

!!******
!******************************************************************************
!!****** environment/set_environment [1.0] *
!!
!!   NAME
!!     set_environment
!!
!!   FUNCTION
!!     Calls for UNIX environment variabls
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE set_environment ( globenv )

  IMPLICIT NONE

! Arguments
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv

! Locals
  INTEGER :: narg

  INTEGER :: iargc

!------------------------------------------------------------------------------

  CALL m_hostnm ( hname )
  CALL m_getcwd ( curdir )
  CALL m_getlog ( user )
  CALL m_getuid ( my_uid )
  CALL m_getpid ( my_pid )

  narg = iargc()

  SELECT CASE ( narg )
  CASE (0)
    globenv%input_file_name = "input"
  CASE DEFAULT
    CALL m_getarg (1,globenv%input_file_name)
  END SELECT

  END SUBROUTINE set_environment

!!*****
!******************************************************************************
!!****** environment/trailer [1.0] *
!!
!!   NAME
!!     trailer
!!
!!   FUNCTION
!!     Writes final timings and banner
!!
!!   AUTHOR
!!     JGH
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!   SOURCE
!******************************************************************************

SUBROUTINE trailer ( globenv )

  IMPLICIT NONE

! Arguments
  TYPE ( global_environment_type ), INTENT ( INOUT ) :: globenv ! check intent

! Locals
  CHARACTER ( LEN = 26 ) :: datx
  INTEGER :: iw, l

!------------------------------------------------------------------------------

  CALL timestop ( zero, handle )

  IF ( globenv % ionode ) THEN

     iw = globenv % scr

     CALL timeprint ( iw, globenv % print_level )

     WRITE ( iw, '( )' )

     CALL m_datum ( datx )

     IF (globenv%print%program_banner) THEN
       WRITE ( iw, &
          '( A, T55, A26, /, A, T49, A32, /, A, T71, A10, /, A, T71, I10 )' ) &
          '  **** **** ******  **  PROGRAM ENDED AT   ', ADJUSTR ( datx ), &
          ' ***** ** ***  *** **   PROGRAM RAN ON     ', ADJUSTR ( hname ), &
          ' **    ****   ******    PROGRAM RAN BY     ', ADJUSTR ( user ), &
          ' ***** **    ** ** **   PROGRAM PROCESS ID ', my_pid

       l = LEN_TRIM ( curdir )
       IF ( l <= 38 ) THEN
          WRITE ( iw, '( A, T43, A38 )' ) &
               '  **** **  *******  **  PROGRAM RAN IN     ', &
               ADJUSTR ( curdir ( 1 : l ) )
       ELSE
          WRITE ( iw, '( A )' ) &
            '  **** **  *******  **  PROGRAM RAN IN     '
          WRITE ( iw, '( T21, A60 )' ) ADJUSTR ( curdir )
       END IF
      END IF
    END IF

  END SUBROUTINE trailer

!!*****
!******************************************************************************

END MODULE environment
