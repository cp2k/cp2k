!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations
!   Copyright (C) 2001 - 2004  CP2K developers group
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
!!     - init_mathcon is now called (27.02.2001,MK)
!!     - Merged with Quickstep MODULE start_program_run (17.01.2002,MK)
!!     - Compile information added (16.01.2002,MK)
!!     - Merged with MODULE cp2k_input, some rearrangements (30.10.2002,MK)
!!
!!******
!******************************************************************************

MODULE environment
  USE cp2k_units,                      ONLY: init_cp2k_units
  USE cp_log_handling,                 ONLY: cp_init_default_logger, &
                                             cp_destroy_default_logger
  USE cp_para_env,                     ONLY: cp_para_env_release, & 
                                             cp_para_env_create
  USE fft_tools,                       ONLY: fft3d,&
                                             init_fft
  USE global_types,                    ONLY: BLACS_GRID_COL,&
                                             BLACS_GRID_ROW,&
                                             BLACS_GRID_SQUARE,&
                                             DEBUG,&
                                             HIGH,&
                                             LOW,&
                                             MEDIUM,&
                                             SILENT,&
                                             WANNIER_CRAZY,&
                                             WANNIER_DIRECT,&
                                             WANNIER_JACOBI,&
                                             WANNIER_NONE,&
                                             global_environment_type,&
                                             init_print_keys
  USE kinds,                           ONLY: dp,&
                                             print_kind_info
  USE machine,                         ONLY: m_datum,&
                                             m_getarg,&
                                             m_getcwd,&
                                             m_getlog,&
                                             m_getpid,&
                                             m_hostnm,&
                                             m_iargc
  USE mathconstants,                   ONLY: init_mathcon
  USE message_passing,                 ONLY: mp_bcast, &
                                             mp_perf_end, &
                                             mp_perf_start, &
                                             mp_environ
  USE physcon,                         ONLY: init_physcon,&
                                             write_physcon
  USE qs_parser,                       ONLY: finish_parser,&
                                             open_file,&
                                             read_object,&
                                             start_parser,&
                                             stop_parser,&
                                             test_object
  USE termination,                     ONLY: set_error_unit,&
                                             stop_program
  USE timings,                         ONLY: timeprint,&
                                             timeset,&
                                             timestop,&
                                             trace_debug,&
                                             init_timer_env,&
                                             destroy_timer_env,&
                                             use_HPM
  USE util,                            ONLY: ran2

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
#if defined(__COMPILE_LASTCVS)
  CHARACTER(LEN=*), PARAMETER :: compile_lastcvs = __COMPILE_LASTCVS
#else
  CHARACTER(LEN=*), PARAMETER :: compile_lastcvs = "/No CVS/Entry in src dir//"
#endif
  CHARACTER(LEN=*), PARAMETER :: module_name = "environment"

  INTEGER :: handle = 0

! *** Public subroutines ***

  PUBLIC :: cp2k_finalize,&
            cp2k_init

!******************************************************************************

CONTAINS

!******************************************************************************
!!****f* environment/cp2k_init [1.0] *
!!
!!   NAME
!!     cp2k_init
!!
!!   SYNOPSIS
!!     Subroutine cp2k_init(globenv)
!!       Type(global_environment_type), Intent (INOUT):: globenv
!!     End Subroutine cp2k_init
!!
!!   FUNCTION
!!     Initializes a CP2K run (setting of the global environment variables)
!!
!!   AUTHOR
!!     JGH,MK
!!
!!   MODIFICATION HISTORY
!!     JGH (28.11.2001) : default for pp_library_path
!!     - print keys added (17.01.2002, MK)
!!     - merged with cp2k_input (30.10.2002,MK)
!!
!!******
!******************************************************************************

  SUBROUTINE cp2k_init(globenv)

    TYPE(global_environment_type), INTENT(INOUT) :: globenv

!   *** Local variables ***

    CHARACTER(LEN=60) :: cwd
    CHARACTER(LEN=30) :: host_name,user_name
    CHARACTER(LEN=26) :: datx
    INTEGER           :: initial_random,l,narg,output_unit,pid,strstart,strend
    LOGICAL           :: ionode
#if defined(__HPM)
#include "f_hpm.h"    
#endif

!   ---------------------------------------------------------------------------

    CALL m_getpid(pid)
    CALL m_getlog(user_name)
    CALL m_hostnm(host_name)

    ! set up the necessary varialbes in globenv
    CALL mp_environ ( globenv % num_pe, globenv % mepos, globenv % group )
    globenv % source = 0
    globenv % ionode = (globenv % mepos==globenv % source)
    nullify(globenv%para_env)
    CALL cp_para_env_create (globenv%para_env, group=globenv%group,&
        owns_group=.false.,source=globenv%source)

    ! only on ionode, otherwise we need a local logger for just a single line
    ! if this info is needed on the screen, send it to ionode for output
    IF (globenv % ionode) &
        WRITE (UNIT=globenv%scr,FMT="(T2,A,T71,I10)")&
          TRIM(user_name)//"@"//TRIM(host_name)//" has created process number",pid

    ! create a timer_env
    CALL init_timer_env()

    ! message passing performance
    CALL mp_perf_start( )

!   *** Set the default logical output and error unit number ***

    CALL set_error_unit(globenv%scr)

!   *** Init the default logger 

    globenv%print%ref_count=1 ! ugly
    CALL cp_init_default_logger(para_env=globenv%para_env,&
         print_keys=globenv%print, &
         default_global_unit_nr=globenv%scr, &
         close_global_unit_on_dealloc=.FALSE.)

    ionode = globenv%ionode

    IF (ionode) THEN

      CALL m_datum(datx)
      CALL m_getcwd(cwd)

    END IF

    CALL mp_bcast(globenv%input_file_name,globenv%source,globenv%group)

!   *** Initialize timing ***

    CALL timeset ("CP2K","I","",handle)

!   *** Initialize random number generator ***

    globenv%idum = -1
    initial_random = ran2(globenv%idum)

!   *** Initialize mathematical constants ***

    CALL init_mathcon()

!   *** Initialize physical constants ***

    CALL init_physcon()
    CALL init_cp2k_units()

!   *** Read the input section with the print keys ***

    CALL read_print_section(globenv)

    output_unit=globenv%scr

    IF (ionode.AND.globenv%print%program_banner) THEN
      WRITE (UNIT=output_unit,&
             FMT="(/,A,T55,A26,/,A,T51,A30,/,A,T51,A30,/,A,T71,I10)")&
        "  **** **** ******  **  PROGRAM STARTED AT ",ADJUSTR(datx),&
        " ***** ** ***  *** **   PROGRAM STARTED ON ",ADJUSTR(host_name),&
        " **    ****   ******    PROGRAM STARTED BY ",ADJUSTR(user_name),&
        " ***** **    ** ** **   PROGRAM PROCESS ID ",pid
      l = LEN_TRIM(cwd)
      IF (l <= 38) THEN
        WRITE (UNIT=output_unit,FMT="(A,T43,A38,/)")&
          "  **** **  *******  **  PROGRAM STARTED IN ",ADJUSTR(cwd(1:l))
      ELSE
        WRITE (UNIT=output_unit,FMT="(A)")&
          "  **** **  *******  **  PROGRAM STARTED IN "
        WRITE (UNIT=output_unit,FMT="(T21,A60,/)") ADJUSTR(cwd)
      END IF
#if defined(__COMPILE_DATE)
      WRITE (UNIT=output_unit,FMT="(T2,A,T31,A50)")&
        "CP2K| Program compiled at",&
        ADJUSTR(compile_date(1:MIN(50,LEN(compile_date))))
#endif
#if defined(__COMPILE_HOST)
      WRITE (UNIT=output_unit,FMT="(T2,A,T31,A50)")&
        "CP2K| Program compiled on",&
        ADJUSTR(compile_host(1:MIN(50,LEN(compile_host))))
#endif
#if defined(__COMPILE_ARCH)
      WRITE (UNIT=output_unit,FMT="(T2,A,T31,A50)")&
        "CP2K| Program compiled for",&
        ADJUSTR(compile_arch(1:MIN(50,LEN(compile_arch))))
#endif
      strend=MAX(1,LEN(compile_lastcvs)-2)
      strstart=MAX(2,strend-50)
      WRITE (UNIT=output_unit,FMT="(T2,A,T31,A50)")&
        "CP2K| Last CVS entry",&
        ADJUSTR(compile_lastcvs(strstart:strend))
    END IF


!   *** Read the input/output section ***

    CALL read_io_section(globenv)

#if defined(__HPM)
    IF (globenv%print%HPM) THEN
        CALL f_hpminit(globenv%para_env%mepos,"CP2K") 
        use_HPM = .TRUE.
    ENDIF
#endif

!   *** Read the CP2K section ***

    CALL read_cp2k_section(globenv)

    IF (ionode) THEN
      IF (globenv%print%basic_data_types) CALL print_kind_info(output_unit)
      IF (globenv%print%physcon) CALL write_physcon(output_unit)
    END IF

  END SUBROUTINE cp2k_init

! *****************************************************************************

!******************************************************************************
!!****f* environment/read_print_section [1.0] *
!!
!!   NAME
!!     read_print_section
!!
!!   FUNCTION
!!     Read the input section with the print keys.
!!
!!   AUTHOR
!!     MK
!!
!!   MODIFICATION HISTORY
!!     Creation (13.07.2000,MK)
!!
!!******
!******************************************************************************
  SUBROUTINE read_print_section(globenv)

    TYPE(global_environment_type), &
      INTENT(INOUT)                          :: globenv

    CHARACTER(LEN=*), PARAMETER :: routine_name = "read_print_section", &
      start_section_label = "PRINT"

    CHARACTER(LEN=3)                         :: test_result
    CHARACTER(LEN=60)                        :: keyword
    INTEGER                                  :: ipos1, ipos2
    LOGICAL                                  :: print_request, section_found

!   ---------------------------------------------------------------------------
!   *** Initialize print keys (the default print level is LOW) ***

    globenv%print%level = LOW

!   *** Modify the default precision of flexible format output ***

    globenv%print%precision = 0

    CALL start_parser(file_name=globenv%input_file_name,&
                      globenv=globenv,&
                      start_section_label=start_section_label,&
                      section_found=section_found)

!   *** Check for print level specification ***

    IF (section_found) THEN
      test_result = test_object()
      IF (test_result == "EOL") THEN
        globenv%print%level = LOW
      ELSE IF (test_result == "INT") THEN
        CALL read_object(globenv%print%level)
        IF (globenv%print%level < SILENT) globenv%print%level = SILENT
        IF (globenv%print%level > DEBUG) globenv%print%level = DEBUG
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
        CASE ("DEBUG","ALL")
          globenv%print%level = DEBUG
        CASE DEFAULT
          CALL stop_parser(routine_name,module_name,__LINE__,&
                           "Invalid definition <"//TRIM(keyword)//&
                           "> found for the default print level")
        END SELECT
      END IF
    END IF

    CALL init_print_keys(globenv%print,globenv%print%level)

    IF (section_found) THEN
      DO
        test_result = test_object(newline=.TRUE.)
        IF (test_result == "EOS") THEN
          EXIT
        ELSE IF (test_result == "EOF") THEN
          CALL stop_parser(routine_name,module_name,__LINE__,"EOF")
        ELSE
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
          CASE ("CORE_HAMILTONIAN_MATRIX","H_MATRIX","H-MATRIX")
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
            globenv%print%kind_information = print_request
            globenv%print%basis_set = print_request
            globenv%print%potential = print_request
          CASE ("KIND_RADII")
            globenv%print%kind_radii = print_request
          CASE ("KINETIC_ENERGY_MATRIX","T_MATRIX","T-MATRIX")
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
            globenv%print%sac_ppnl_neighbor_lists = print_request
            globenv%print%sbc_ppnl_neighbor_lists = print_request
          CASE ("OPTICAL_CONDUCTIVITY")
            globenv%print%optical_conductivity = print_request
          CASE ("ORTHOGONALISATION_MATRIX","ORTHO_MATRIX")
            globenv%print%ortho_matrix = print_request
          CASE ("OVERLAP_MATRIX","S_MATRIX","S-MATRIX")
            globenv%print%overlap_matrix = print_request
          CASE ("PGF_RADII")
            globenv%print%pgf_radii = print_request
          CASE ("PHYSICAL_CONSTANTS","PHYSCON")
            globenv%print%physcon = print_request
          CASE ("POTENTIALS","POTENTIAL")
            globenv%print%potential = print_request
          CASE ("POTENTIAL_ENERGY_MATRIX","V_MATRIX","V-MATRIX")
            globenv%print%potential_energy_matrix = print_request
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
          CASE ("SAC_PPNL_NEIGHBOR_LIST")
            globenv%print%sac_ppnl_neighbor_lists = print_request
          CASE ("SBC_PPNL_NEIGHBOR_LIST")
            globenv%print%sbc_ppnl_neighbor_lists = print_request
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
          CASE ("W-MATRIX","WMATRIX")
            globenv%print%w_matrix = print_request
!MK
          CASE ("CUBE","CUBES")
            globenv%print%cubes = print_request
          CASE ("E_DENSITY_CUBE","DENSITY")
            globenv%print%e_density_cube = print_request
          CASE ("V_HARTREE_CUBE","HARTREE")
            globenv%print%v_hartree_cube = print_request
          CASE ("LUMO","LUMOS")
            globenv%print%lumos = print_request
            IF (test_object() == "INT") CALL read_object(globenv%print%nlumos)
            IF (globenv%print%nlumos <= 0) globenv%print%nlumos = 4
          CASE ("HOMO","HOMOS")
            globenv%print%homos = print_request
            IF (test_object() == "INT") CALL read_object(globenv%print%nhomos)
            IF (globenv%print%nhomos <= 0) globenv%print%nhomos = 4
          CASE ("LOCALIZATION")
            globenv%print%localization = print_request
            IF (test_object() == "STR") THEN
               CALL read_object(keyword,lower_to_upper=.TRUE.)
               SELECT CASE(TRIM(keyword))
               CASE("JACOBI")
                  globenv%print%localization_method = WANNIER_JACOBI
               CASE("CRAZY")
                  globenv%print%localization_method = WANNIER_CRAZY
               CASE("DIRECT")
                  globenv%print%localization_method = WANNIER_DIRECT
               CASE("NONE")
                  globenv%print%localization_method = WANNIER_NONE
               CASE DEFAULT
                  CALL stop_parser(routine_name,module_name,__LINE__,&
                        "INVALID_KEYWORD")
               END SELECT
            ENDIF
          CASE ("MOLECULAR_STATES")
            globenv%print%molecular_states = print_request
          CASE ("TIMECUMUL")
            globenv%print%timecumul = print_request
          CASE ("HPM")
            globenv%print%HPM = print_request
!MK
          CASE ("EP_MATRIXES")
             globenv%print%ep_matrixes = print_request
          CASE ("EP_RHO1_CUBE")
             globenv%print%ep_rho1_cube = print_request
          CASE ("EP_OPERATOR")
             globenv%print%ep_operator = print_request
          CASE ("MULLIKEN","MULLIKEN_ANALYSIS")
             globenv%print%mulliken_analysis = print_request
          CASE ("ORTHO","MO_ORTHONORMALITY")
             globenv%print%mo_orthonormality = print_request
          CASE DEFAULT
            CALL stop_parser(routine_name,module_name,__LINE__,&
                        "INVALID_KEYWORD")
          END SELECT
        END IF
      END DO
    END IF

    CALL finish_parser()

  END SUBROUTINE read_print_section

!******************************************************************************
!!****f* environment/read_io_section [1.0] *
!!
!!   NAME
!!     read_io_section
!!
!!   FUNCTION
!!     Read the input/output section (path and file names).
!!
!!   AUTHOR
!!     MK
!!
!!   MODIFICATION HISTORY
!!     Creation (13.07.2000,MK)
!!
!!******
!******************************************************************************
  SUBROUTINE read_io_section(globenv)

    TYPE(global_environment_type), &
      INTENT(INOUT)                          :: globenv

    CHARACTER(LEN=*), PARAMETER :: routine_name = "read_io_section", &
      start_section_label = "IO"

    CHARACTER(LEN=3)                         :: test_result
    CHARACTER(LEN=60)                        :: keyword
    INTEGER                                  :: output_unit, unit_number
    LOGICAL                                  :: section_found

!   ---------------------------------------------------------------------------
!   *** Load the default settings ***

    globenv%pp_library_path = "./"

    globenv%basis_set_file_name = "BASIS_SET"
    globenv%coord_file_name = globenv%input_file_name
    globenv%error_file_name = "(stdout)"
    globenv%output_file_name = "(stdout)"
    globenv%potential_file_name = "POTENTIAL"
    globenv%restart_file_name = "RESTART"

    CALL start_parser(file_name=globenv%input_file_name,&
                      globenv=globenv,&
                      start_section_label=start_section_label,&
                      section_found=section_found)

    IF (section_found) THEN
      DO
        test_result = test_object()
        IF (test_result /= "EOL") THEN
          CALL read_object(keyword)
          CALL stop_parser(routine_name,module_name,__LINE__,&
                      "INVALID_KEYWORD")
        END IF
        test_result = test_object(newline=.TRUE.)
        IF (test_result == "EOS") THEN
          EXIT
        ELSE IF (test_result == "EOF") THEN
          CALL stop_parser(routine_name,module_name,__LINE__,"EOF")
        ELSE
          CALL read_object(keyword,lower_to_upper=.TRUE.,newline=.TRUE.)
          SELECT CASE (TRIM(keyword))
          CASE ("BASIS_SET_FILE_NAME","BASIS_SET_FILE")
            CALL read_object(globenv%basis_set_file_name)
          CASE ("COORD_FILE_NAME","COORD_FILE")
            CALL read_object(globenv%coord_file_name)
          CASE ("ERROR_FILE_NAME","ERROR_FILE")
            CALL read_object(globenv%error_file_name)
            CALL open_file(file_name=globenv%error_file_name,&
                           file_action="WRITE",&
                           file_form="FORMATTED",&
                           file_status="REPLACE",&
                           unit_number=unit_number)
            CALL set_error_unit(unit_number)
          CASE ("OUTPUT_FILE_NAME","OUTPUT_FILE")
            CALL read_object(globenv%output_file_name)
            CALL open_file(file_name=globenv%output_file_name,&
                           file_action="WRITE",&
                           file_form="FORMATTED",&
                           file_status="REPLACE",&
                           unit_number=globenv%scr)
          CASE ("POTENTIAL_FILE_NAME","POTENTIAL_FILE")
            CALL read_object(globenv%potential_file_name)
          CASE ("PP_LIBRARY_PATH")
            CALL read_object(globenv%pp_library_path)
          CASE ("RESTART_FILE_NAME","RESTART_FILE")
            CALL read_object(globenv%restart_file_name)
          CASE DEFAULT
            CALL stop_parser(routine_name,module_name,__LINE__,&
                      "INVALID_KEYWORD")
          END SELECT
        END IF
      END DO
    END IF

    CALL finish_parser()

    IF (globenv%ionode.AND.(globenv%print%level > SILENT)) THEN
      output_unit = globenv%scr
      WRITE (UNIT=output_unit,FMT="(/,(T2,A,T41,A))")&
        start_section_label//"| Input file name",&
        ADJUSTR(globenv%input_file_name(:40)),&
        start_section_label//"| Output file name",&
        ADJUSTR(globenv%output_file_name(:40)),&
        start_section_label//"| Error file name",&
        ADJUSTR(globenv%error_file_name(:40)),&
        start_section_label//"| Library path",&
        ADJUSTR(globenv%pp_library_path(:40)),&
        start_section_label//"| Basis set file name",&
        ADJUSTR(globenv%basis_set_file_name(:40)),&
        start_section_label//"| Potential file name",&
        ADJUSTR(globenv%potential_file_name(:40)),&
        start_section_label//"| Restart file name",&
        ADJUSTR(globenv%restart_file_name(:40)),&
        start_section_label//"| Coordinate file name",&
        ADJUSTR(globenv%coord_file_name(:40))
    END IF

  END SUBROUTINE read_io_section

!******************************************************************************
!!****f* environment/read_cp2k_section [1.0] *
!!
!!   NAME
!!     read_cp2k_section
!!
!!   FUNCTION
!!
!!   AUTHOR
!!     JGH,MK
!!
!!   MODIFICATION HISTORY
!!     2-Dec-2000 (JGH) added default fft library
!!
!!******
!******************************************************************************

  SUBROUTINE read_cp2k_section(globenv)

    TYPE(global_environment_type), &
      INTENT(INOUT)                          :: globenv

    CHARACTER(LEN=*), PARAMETER :: routine_name = "read_cp2k_section", &
      start_section_label = "CP2K"

    CHARACTER(LEN=3)                         :: test_result
    CHARACTER(LEN=60)                        :: keyvalue, keyword
    COMPLEX(KIND=dp), DIMENSION(4, 4, 4)     :: zz
    INTEGER                                  :: num_threads, &
                                                omp_get_num_threads, &
                                                output_unit, stat
    INTEGER, DIMENSION(3)                    :: n
    LOGICAL                                  :: section_found

!   ---------------------------------------------------------------------------
!   *** Load the default settings ***

    globenv%blacs_grid_layout = BLACS_GRID_SQUARE
    globenv%default_fft_library = "FFTSG"
    globenv%print_level = 0
    globenv%program_name = " "
    globenv%project_name = "unspecified"
    globenv%run_type = "default"

    CALL start_parser(file_name=globenv%input_file_name,&
                      globenv=globenv,&
                      start_section_label=start_section_label,&
                      section_found=section_found)

    IF (section_found) THEN
      DO
        test_result = test_object()
        IF (test_result /= "EOL") THEN
          CALL read_object(keyword)
          CALL stop_parser(routine_name,module_name,__LINE__,&
                       "INVALID_KEYWORD")
        END IF
        test_result = test_object(newline=.TRUE.)
        IF (test_result == "EOS") THEN
          EXIT
        ELSE IF (test_result == "EOF") THEN
          CALL stop_parser(routine_name,module_name,__LINE__,"EOF")
        ELSE
          CALL read_object(keyword,lower_to_upper=.TRUE.,newline=.TRUE.)
          SELECT CASE (TRIM(keyword))
          CASE ("BLACS_GRID")
            CALL read_object(keyvalue,lower_to_upper=.TRUE.)
            SELECT CASE(TRIM(keyvalue))
            CASE ("SQUARE")
                 globenv%blacs_grid_layout = BLACS_GRID_SQUARE
            CASE ("ROW")
                 globenv%blacs_grid_layout = BLACS_GRID_ROW
            CASE ("COL","COLUMN")
                 globenv%blacs_grid_layout = BLACS_GRID_COL
            CASE DEFAULT
              CALL stop_parser(routine_name,module_name,__LINE__,&
                       "INVALID_KEYWORD")
            END SELECT
          CASE ("FFT_LIBRARY","FFT_LIB","FFTLIB")
            CALL read_object(globenv%default_fft_library,lower_to_upper=.TRUE.)
          CASE ("IOLEVEL","PRINT_LEVEL")
            CALL read_object(globenv%print_level)
          CASE ("PROGRAM_NAME","PROGRAM")
            CALL read_object(globenv%program_name,lower_to_upper=.TRUE.)
          CASE ("PROJECT_NAME","PROJECT")
            CALL read_object(globenv%project_name)
          CASE ("RUN_TYPE")
            CALL read_object(globenv%run_type,lower_to_upper=.TRUE.)
          CASE ("TRACE")
            CALL trace_debug("start",groupid=globenv%group)
          CASE DEFAULT
            CALL stop_parser(routine_name,module_name,__LINE__,&
                       "INVALID_KEYWORD")
          END SELECT
        END IF
      END DO
    ELSE
      CALL stop_parser(routine_name,module_name,__LINE__,"SECTION_NOT_FOUND")
    END IF

    CALL finish_parser()

!   *** Check and unify arguments ***

    SELECT CASE (globenv%program_name)
    CASE ("ATOM")
    CASE ("FARMING")
    CASE ("FIST")
    CASE ("KIM-GORDON","KG")
      globenv%program_name = "KG"
    CASE ("QUICKSTEP","QS")
      globenv%program_name = "QS"
    CASE ("EP")
    CASE ("TRAJANA")
    CASE ("TEST")
    CASE (" ")
      CALL stop_program(routine_name,module_name,__LINE__,&
                        "No PROGRAM_NAME was specified in the "//&
                        TRIM(start_section_label)//" section",globenv)
    CASE DEFAULT
      CALL stop_program(routine_name,module_name,__LINE__,&
                        "Invalid program name <"//&
                        TRIM(globenv%program_name)//"> found",globenv)
    END SELECT

    SELECT CASE (globenv%run_type)
    CASE ("GEOMETRY_OPTIMIZATION","GEO_OPT")
      globenv%run_type = "GEOMETRY OPTIMIZATION"
    CASE ("MONTE_CARLO","MC")
      globenv%run_type = "MONTE CARLO"
    CASE ("MOLECULAR_DYNAMICS","MD")
      globenv%run_type = "MOLECULAR DYNAMICS"
    CASE ("WAVEFUNCTION_OPTIMIZATION","WFN_OPT")
      globenv%run_type = "WAVEFUNCTION OPTIMIZATION"
    CASE ("DEBUG_FORCES","DEBUG")
      globenv%run_type = "DEBUG"
    CASE ("default")
       SELECT CASE (globenv%program_name)
       CASE ("ATOM")
          globenv%run_type = "NONE"
       CASE ("FARMING")
          globenv%run_type = "NONE"
       CASE ("FIST")
          globenv%run_type = "MOLECULAR DYNAMICS"
       CASE ("KIM-GORDON","KG")
          globenv%run_type = "MOLECULAR DYNAMICS"
       CASE ("QUICKSTEP","QS")
          globenv%run_type = "WAVEFUNCTION OPTIMIZATION"
       CASE ("TEST")
          globenv%run_type = "NONE"
       CASE ("TRAJANA")
          globenv%run_type = "ANALYSIS"
       CASE DEFAULT
         CALL stop_program(routine_name,module_name,__LINE__,&
                           "Invalid program name <"//&
                           TRIM(globenv%program_name)//"> found",globenv)
       END SELECT
    CASE (" ")
      CALL stop_program(routine_name,module_name,__LINE__,&
                        "No RUN_TYPE was specified in the "//&
                        TRIM(start_section_label)//" section",globenv)
    CASE DEFAULT
      CALL stop_program(routine_name,module_name,__LINE__,&
                        "Invalid run type <"//&
                        TRIM(globenv%run_type)//"> found",globenv)
    END SELECT

    SELECT CASE(globenv%program_name)
    CASE ("ATOM")
       IF (globenv%run_type/="NONE") THEN
          CALL stop_program(routine_name,module_name,__LINE__,&
               "ATOM program supports only NONE as run type",globenv)
       END IF
    CASE ("FARMING")
       IF (globenv%run_type/="NONE") THEN
          CALL stop_program(routine_name,module_name,__LINE__,&
               "ATOM program supports only NONE as run type",globenv)
       END IF
    CASE ("TEST")
       IF (globenv%run_type/="NONE") THEN
          CALL stop_program(routine_name,module_name,__LINE__,&
               "TEST program supports only NONE as run type",globenv)
       END IF
    END SELECT

!   *** Initialize FFT library ***
    CALL init_fft(fftlib=TRIM(globenv%default_fft_library))

!   *** Check for FFT library ***
    n(:) = 4
    zz(:,:,:) = 0.0_dp
    CALL fft3d(1,n,zz,status=stat)
    IF (stat /= 0) THEN
      IF (globenv%ionode) THEN
         write(globenv%scr,*) &
           " WARNING : FFT library "//TRIM(globenv%default_fft_library)//&
           " is not available "
         write(globenv%scr,*) "           Trying FFTSG as a default "
      ENDIF
      globenv%default_fft_library="FFTSG"
      CALL init_fft(fftlib=TRIM(globenv%default_fft_library))
      CALL fft3d(1,n,zz,status=stat) 
      IF (stat /= 0) THEN
           CALL stop_program(routine_name,module_name,__LINE__,&
                             "The default FFT library <"//&
                             TRIM(globenv%default_fft_library)//&
                             "> is not available",globenv)
      ENDIF
    END IF

    num_threads=1
!$omp parallel
!$ num_threads = omp_get_num_threads()
!$omp end parallel

    IF (globenv%ionode.AND.(globenv%print%level > SILENT)) THEN
      output_unit = globenv%scr
      WRITE (UNIT=output_unit,FMT="(T2,A,T41,A)")&
        start_section_label//"| Program name",&
        ADJUSTR(globenv%program_name(:40)),&
        start_section_label//"| Project name",&
        ADJUSTR(globenv%project_name(:40)),&
        start_section_label//"| Default FFT library",&
        ADJUSTR(globenv%default_fft_library(:40)),&
        start_section_label//"| Run type",&
        ADJUSTR(globenv%run_type(:40))
      WRITE (UNIT=output_unit,FMT="(T2,A,T75,I6)")&
        start_section_label//"| Global print level",globenv%print_level,&
        start_section_label//"| Total number of message passing processes",&
                             globenv%num_pe,&
        start_section_label//"| Number of threads for this process",&
                             num_threads,&
        start_section_label//"| This output is from process",globenv%mepos
      WRITE (UNIT=output_unit,FMT='()')
    END IF

  END SUBROUTINE read_cp2k_section

!******************************************************************************
!!****f* environment/cp2k_finalize [1.0] *
!!
!!   NAME
!!     cp2k_finalize
!!
!!   SYNOPSIS
!!     Subroutine cp2k_finalize(globenv)
!!       Type(global_environment_type), Intent (INOUT):: globenv
!!     End Subroutine cp2k_finalize
!!
!!   FUNCTION
!!     Writes final timings and banner for CP2K
!!
!!   AUTHOR
!!     JGH,MK
!!
!!   MODIFICATION HISTORY
!!     none
!!
!!******
!******************************************************************************

  SUBROUTINE cp2k_finalize(globenv)

    TYPE(global_environment_type), INTENT(INOUT) :: globenv

!   *** Local variables ***

    CHARACTER(LEN=60) :: cwd
    CHARACTER(LEN=30) :: host_name,user_name
    CHARACTER(LEN=26) :: datx
    INTEGER           :: l,output_unit,pid
#if defined(__HPM)
#include "f_hpm.h"
#endif

!   ---------------------------------------------------------------------------

    ! write message passing performance info
    CALL mp_perf_end ( globenv%group,globenv%scr )

#if defined(__HPM)
    IF (globenv%print%HPM) THEN
        CALL f_hpmterminate(globenv%group) 
        use_HPM=.FALSE.
    ENDIF
#endif

    CALL timestop(0.0_dp,handle) ! corresponding the "CP2K" in cp2k_init

    IF (globenv%ionode) THEN
      output_unit = globenv%scr
      CALL timeprint(output_unit,globenv%print_level)
    END IF

    CALL destroy_timer_env()

    IF (globenv%ionode) THEN
      output_unit = globenv%scr
      CALL m_datum(datx)
      CALL m_getcwd(cwd)
      CALL m_getlog(user_name)
      CALL m_getpid(pid)
      CALL m_hostnm(host_name)

      IF (globenv%ionode.AND.globenv%print%program_banner) THEN
        WRITE (UNIT=output_unit,&
               FMT="(A,T55,A26,/,A,T51,A30,/,A,T51,A30,/,A,T71,I10)")&
         "  **** **** ******  **  PROGRAM ENDED AT   ",ADJUSTR(datx),&
         " ***** ** ***  *** **   PROGRAM RAN ON     ",ADJUSTR(host_name),&
         " **    ****   ******    PROGRAM RAN BY     ",ADJUSTR(user_name),&
         " ***** **    ** ** **   PROGRAM PROCESS ID ",pid
        l = LEN_TRIM(cwd)
        IF (l <= 38) THEN
          WRITE (UNIT=output_unit,FMT="(A,T43,A38,/)")&
            "  **** **  *******  **  PROGRAM STARTED IN ",ADJUSTR(cwd(1:l))
        ELSE
          WRITE (UNIT=output_unit,FMT="(A)")&
            "  **** **  *******  **  PROGRAM STARTED IN "
          WRITE (UNIT=output_unit,FMT="(T21,A60,/)") ADJUSTR(cwd)
        END IF
      END IF

    END IF

!   *** Release message passing environment ***
    CALL cp_destroy_default_logger()
    CALL cp_para_env_release(globenv%para_env)

  END SUBROUTINE cp2k_finalize

! *****************************************************************************

END MODULE environment

! *****************************************************************************
