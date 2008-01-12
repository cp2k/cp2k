!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!!****h* cp2k/input_cp2k_vib [1.0] *
!!
!!   NAME
!!     input_cp2k_vib
!!
!!   FUNCTION
!!     builds the input structure for the VIBRATIONAL_ANALYSIS module
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     [tlaino]
!!
!!   MODIFICATION HISTORY
!!     01.2008 [tlaino] Teodoro Laino - University of Zurich
!!                      Creating an own module for vibrational analysis
!!
!****************************************************************************
MODULE input_cp2k_vib
  USE cp_output_handling,              ONLY: cp_print_key_section_create
  USE f77_blas
  USE input_constants
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
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
  
  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_vib'
  
  PUBLIC :: create_vib_section
!!***
!****************************************************************************
CONTAINS

!!****f* input_cp2k_vib/create_vib_section [1.0] *
!!
!!   NAME
!!     create_vib_section
!!
!!   FUNCTION
!!     Creates the exteranal restart section
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - section: the section to create
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     tlaino
!!
!!*** **********************************************************************
  SUBROUTINE create_vib_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_vib_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="VIBRATIONAL_ANALYSIS",&
            description="Section to setup parameters to perform a Normal Modes analysis.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            supported_feature=.TRUE.,error=error)
       NULLIFY(keyword, print_key,subsection)
       
       CALL keyword_create(keyword, name="DX",&
            description="Specify the increment to be used to construct the HESSIAN with "//&
            "finite difference method",&
            default_r_val=1.0E-2_dp, supported_feature=.TRUE.,error=error)
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

       CALL cp_print_key_section_create(print_key,"PROGRAM_RUN_INFO",&
            description="Controls the printing basic info about the vibrational method", &
            print_level=medium_print_level,each=(/1/),add_last=(/1/),filename="__STD_OUT__",&
            supported_feature=.TRUE., error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"MOLDEN_VIB",&
            description="Controls the printing for visualization in molden format", &
            print_level=low_print_level,each=(/1/),add_last=(/1/),filename="VIBRATIONS",&
            supported_feature=.TRUE., error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)
       
       CALL create_mode_selective_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE create_vib_section
!***************************************************************************

!!****f* input_cp2k/create_mode_selective_section [1.0] *
!!
!!   NAME
!!     create_mode_selective_section
!!
!!   FUNCTION
!!     Create the input section for MODE selective
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - section: the section to create
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fschiff
!!
!!*** **********************************************************************
  SUBROUTINE create_mode_selective_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'create_mode_selective_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    NULLIFY(keyword)
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="MODE_SELECTIVE",&
            description="All parameters needed for to run a mode selective vibrational analysis",&
            n_keywords=0, n_subsections=5, repeats=.FALSE., required=.TRUE.,&
            error=error)

       CALL keyword_create(keyword,name="FREQUENCY",&
            description="value close to the expected value of the frequency for to look for. "//&
            "If the block Davidson algorithm is applied, the nrep closest frequencies are tracked. ",&
            usage="FREQUENCY <REAL>", default_r_val=-1._dp,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="RANGE",&
            description="Track modes in a given range of frequencies. "//&
            "No warranty that the set of frequencies is complete.",&
            usage="FREQUENCY <REAL> <REAL>",&
            n_var=-1,type_of_var=real_t,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ATOMS",&
            description="Specifies the list of atoms which should be mainly involved in the selected frequency",&
            usage="ATOMS {integer} {integer} .. {integer}", required=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="EPS_MAX_VAL",&
            description="Convergence criterium for the davidson algorithm. Specifies the maximal value in the "//&
            "residuum vectors ",&
            usage="EPS_MAX_VAL <REAL>", default_r_val=5.0E-7_dp,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="EPS_NORM",&
            description="Convergence criterium for the davidson algorithm. Specifies the maximal value of the norm "//&
            "of the residuum vectors ",&
            usage="EPS_MAX_VAL <REAL>", default_r_val=2.0E-6_dp,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="INITIAL_GUESS",&
            description="The type of initial guess for the normal modes",&
            usage="INITIAL_GUESS BFGS_HESS",&
            default_i_val=ms_guess_atomic,&
            enum_c_vals=s2a("BFGS_HESS", "ATOMIC", "RESTART"),&
            enum_desc=s2a("get the first displacement vector out of the BFGS approximate Hessian", &
                          "use random displacements for a set of atoms specified",&
                          "use data from MS_RESTART as initial guess"),&
            enum_i_vals=(/ms_guess_bfgs,ms_guess_atomic,ms_guess_restart/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PRECONDITIONER",&
            description="The preconditioner for the davidson algorithm (do not use yet they have to improve first )",&
            usage="PRECONDITIONER NONE",&
            default_i_val=no_ms_precond,&
            enum_c_vals=s2a("NONE","BFGS"),&
            enum_desc=s2a("no preconditioner to be applied", &
                          "BFGS approximate Hessian used for preconditioning"),&
            enum_i_vals=(/no_ms_precond,bfgs_ms_precond/),&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_mode_selective_section

END MODULE input_cp2k_vib
