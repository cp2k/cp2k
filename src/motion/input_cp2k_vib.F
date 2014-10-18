!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief builds the input structure for the VIBRATIONAL_ANALYSIS module
!> \par History
!>      01.2008 [tlaino] Teodoro Laino - University of Zurich
!>                       Creating an own module for vibrational analysis
!> \author [tlaino]
! *****************************************************************************
MODULE input_cp2k_vib
  USE cp_output_handling,              ONLY: add_last_numeric,&
                                             cp_print_key_section_create,&
                                             debug_print_level,&
                                             low_print_level,&
                                             medium_print_level,&
                                             silent_print_level
  USE f77_blas
  USE input_constants,                 ONLY: do_rep_blocked,&
                                             do_rep_interleaved,&
                                             ms_guess_atomic,&
                                             ms_guess_bfgs,&
                                             ms_guess_molden,&
                                             ms_guess_restart,&
                                             ms_guess_restart_vec
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_add_subsection,&
                                             section_create,&
                                             section_release,&
                                             section_type
  USE input_val_types,                 ONLY: integer_t,&
                                             real_t
  USE kinds,                           ONLY: dp
  USE string_utilities,                ONLY: s2a
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_vib'

  PUBLIC :: create_vib_section
CONTAINS

! *****************************************************************************
!> \brief Creates the exteranal restart section
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author tlaino
! *****************************************************************************
  SUBROUTINE create_vib_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_vib_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="VIBRATIONAL_ANALYSIS",&
            description="Section to setup parameters to perform a Normal Modes analysis.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            supported_feature=.TRUE.,error=error)
       NULLIFY(keyword, subsection)

       CALL keyword_create(keyword, name="DX",&
            description="Specify the increment to be used to construct the HESSIAN with "//&
            "finite difference method",&
            default_r_val=1.0E-2_dp, unit_str="bohr", supported_feature=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NPROC_REP",&
            description="Specify the number of processors to be used per replica "//&
            "environment (for parallel runs). "//&
            "In case of mode selective calculations more than one replica will start"//&
            " a block Davidson algorithm to track more than only one frequency",&
            default_i_val=1, supported_feature=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PROC_DIST_TYPE",&
            description="Specify the topology of the mapping of processors into replicas.",&
            usage="PROC_DIST_TYPE (INTERLEAVED|BLOCKED)",&
            enum_c_vals=s2a("INTERLEAVED",&
                            "BLOCKED"),&
            enum_desc=s2a( "Interleaved distribution",&
                           "Blocked distribution"),&
            enum_i_vals=(/do_rep_interleaved,do_rep_blocked/),&
            default_i_val=do_rep_blocked, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="FULLY_PERIODIC",&
            description="Avoids to clean rotations from the Hessian matrix.",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., supported_feature=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="INTENSITIES",&
            description="Calculation of the IR-Intensities. Calculation of dipols has to be specified explicitly ",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., supported_feature=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_mode_selective_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_print_vib_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)
    END IF
  END SUBROUTINE create_vib_section

! *****************************************************************************
!> \brief Create the print section for VIB
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Teodoro Laino [tlaino] - 10.2008
! *****************************************************************************
  SUBROUTINE create_print_vib_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_print_vib_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="PRINT",&
            description="Section controlling the print information during a vibrational "//&
            "analysis.",n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            supported_feature=.TRUE.,error=error)
       NULLIFY(keyword, print_key)

       CALL cp_print_key_section_create(print_key,"BANNER",&
            description="Controls the printing of the vibrational analysis banner", &
            print_level=low_print_level, common_iter_levels=1,&
            filename="__STD_OUT__",error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"PROGRAM_RUN_INFO",&
            description="Controls the printing basic info about the vibrational method", &
            print_level=medium_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
            supported_feature=.TRUE., error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"MOLDEN_VIB",&
            description="Controls the printing for visualization in molden format", &
            print_level=low_print_level,add_last=add_last_numeric,filename="VIBRATIONS",&
            supported_feature=.TRUE., error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"ROTATIONAL_INFO",&
            description="Controls the printing basic info during the cleaning of the "//&
            "rotational degrees of freedom.", &
            print_level=debug_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
            supported_feature=.TRUE., error=error)
       ! Print_key keywords
       CALL keyword_create(keyword, name="COORDINATES",&
            description="Prints atomic coordinates after rotation",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., supported_feature=.TRUE.,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

    END IF
  END SUBROUTINE create_print_vib_section

! *****************************************************************************
!> \brief Create the input section for MODE selective
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fschiff
! *****************************************************************************
  SUBROUTINE create_mode_selective_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'create_mode_selective_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, subsection

    failure=.FALSE.

    NULLIFY(keyword, subsection, print_key)
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="MODE_SELECTIVE",&
            description="All parameters needed for to run a mode selective vibrational analysis",&
            n_keywords=5, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            error=error)

       CALL keyword_create(keyword,name="FREQUENCY",&
            description="value close to the expected value of the frequency for to look for. "//&
            "If the block Davidson algorithm is applied, the nrep closest frequencies are tracked. ",&
            usage="FREQUENCY {REAL}", default_r_val=-1._dp,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="RANGE",&
            description="Track modes in a given range of frequencies. "//&
            "No warranty that the set of frequencies is complete.",&
            usage="RANGE {REAL} {REAL}",&
            n_var=-1,type_of_var=real_t,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ATOMS",&
            description="Specifies the list of atoms which should be displaced for the Initial guess",&
            usage="ATOMS {integer} {integer} .. {integer}", required=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="EPS_MAX_VAL",&
            description="Convergence criterium for the davidson algorithm. Specifies the maximal value in the "//&
            "residuum vectors ",&
            usage="EPS_MAX_VAL {REAL}", default_r_val=5.0E-7_dp,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="EPS_NORM",&
            description="Convergence criterium for the davidson algorithm. Specifies the maximal value of the norm "//&
            "of the residuum vectors ",&
            usage="EPS_NORM {REAL}", default_r_val=2.0E-6_dp,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="INITIAL_GUESS",&
            description="The type of initial guess for the normal modes",&
            usage="INITIAL_GUESS BFGS_HESS",&
            default_i_val=ms_guess_atomic,&
            enum_c_vals=s2a("BFGS_HESS", "ATOMIC", "RESTART","RESTART_VEC","MOLDEN_RESTART"),&
            enum_desc=s2a("get the first displacement vector out of the BFGS approximate Hessian", &
                          "use random displacements for a set of atoms specified",&
                          "use data from MS_RESTART as initial guess",&
                          "use a vector from MS_RESTART, useful if you want to increase accurcy by changing functionals or basis",&
                          "use the .mol file of a former run, to restart a vector"//&
                          "(similar to Restart_vec, but a different file FORMAT is used)"),&
            enum_i_vals=(/ms_guess_bfgs,ms_guess_atomic,ms_guess_restart,ms_guess_restart_vec,ms_guess_molden/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RESTART_FILE_NAME",&
            description="Specifies the name of the file used to create the restarted vectors",&
            usage="RESTART_FILE_NAME {filename}",&
            default_lc_val="",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_involved_atoms_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL section_create(subsection, name="PRINT",&
            description="Controls the printing mode selective vibrational analysis",&
            n_keywords=0, n_subsections=1, repeats=.TRUE., required=.FALSE.,&
            error=error)

       CALL cp_print_key_section_create(print_key,"MS_RESTART",&
            description="Controls the printing of the Mode Selective Restart file.", &
            print_level=silent_print_level, common_iter_levels=1,  &
            add_last=add_last_numeric, filename="",error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE create_mode_selective_section

! *****************************************************************************
!> \brief Create the input section for Ivolved_atoms keyword in mode selective
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fschiff
! *****************************************************************************
  SUBROUTINE create_involved_atoms_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'create_involved_atoms_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    NULLIFY(keyword)
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="INVOLVED_ATOMS",&
            description="All parameters needed for the tracking of modes dominated by the motion of selected atoms ",&
            n_keywords=2, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       CALL keyword_create(keyword,name="RANGE",&
            description=" Specifies the range of wavenumbers in which the modes related to the ATOMS have to be tracked. "//&
            " If not specified frequencies >400cm-1 will be used to avoid tracking of translational or rotational modes",&
            usage="RANGE {REAL} {REAL}",&
            n_var=-1,type_of_var=real_t,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="INVOLVED_ATOMS",&
            description="Specifies the list of atoms on which the tracked eigenvector should have the highest value "//&
            "similar to looking for the vibration of a set of atoms",&
            usage="INVOLVED_ATOMS {integer} {integer} .. {integer}", required=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)



    END IF
  END SUBROUTINE create_involved_atoms_section

END MODULE input_cp2k_vib
