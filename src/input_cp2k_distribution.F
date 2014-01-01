!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief function that builds the distribution section of the input
!> \par History
!>      04.2007 created
!> \author Joost VandeVondele
! *****************************************************************************
MODULE input_cp2k_distribution
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
  USE kinds,                           ONLY: dp
  USE string_utilities,                ONLY: s2a
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.FALSE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_distribution'

  PUBLIC :: create_distribution_section

CONTAINS

! *****************************************************************************
!> \brief Creates the distribution section
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE create_distribution_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_distribution_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="DISTRIBUTION",&
            description="can be used used to tune the parallel distribution of the data",&
            n_keywords=2, n_subsections=2, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword, subsection)

       CALL keyword_create(keyword,name="COST_MODEL",&
            description="The cost model that needs to be minimized ",&
            usage="COST_MODEL BLOCK_COUNT",&
            enum_c_vals=s2a("BLOCK_COUNT","BLOCK_SURFACE", "BLOCK_LMAX"),&
            enum_i_vals=(/model_block_count, model_block_surface,model_block_lmax/),&
            enum_desc=s2a("the number of blocks",&
                          "the number of blocks weighted by the number elements per block",&
                          "the number of blocks weighted by the sum of the lmax"), &
            default_i_val=model_block_count,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="2D_MOLECULAR_DISTRIBUTION",&
            description="Distribute the atoms so that atoms belonging to a given molecule"//&
                        " are on the same CPU for the 2D distribution. This might give rise to a"//&
                        " worse distribution but reduces memory needs of finding the optimal distribution.",&
            usage="2D_MOLECULAR_DISTRIBUTION TRUE",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SKIP_OPTIMIZATION",&
            description="Do not optimize the distribution, go for something very simple."//&
                        " Might be useful if the optimization, which scales quadratically in system size, is too expensive.",&
            usage="SKIP_OPTIMIZATION TRUE",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="BASIC_OPTIMIZATION",&
            description="Creates a distribution based on a few heuristics using only minimal memory "//&
                        "and CPU time.",&
            usage="BASIC_OPTIMIZATION TRUE",&
            default_l_val=.TRUE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="BASIC_SPATIAL_OPTIMIZATION",&
            description="Creates a distribution with spatial info, using only minimal memory "//&
                        "and CPU time.",&
            usage="BASIC_SPATIAL_OPTIMIZATION TRUE",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SYMMETRIC",&
            description="Take the symmetry of the distribution_2d into account.",&
            usage="SYMMETRIC TRUE",&
            default_l_val=.TRUE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_2d_mc_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE create_distribution_section


! *****************************************************************************
!> \brief provides options for the distribute_2d_monte_carlo
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE create_2d_mc_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_2d_mc_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="2D_MC",&
            description="provides parameters for the 2D Monte Carlo procedure",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="MOVES",&
            description="Total number of Monte Carlo moves to use in the optimization",&
            usage="MOVES 1000000 ",&
            n_var=1,default_i_val=400000,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="CYCLES",&
            description="How many annealing cycles to use",&
            usage="CYCLES 10", default_i_val=10, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="REDUCTION_STEPS",&
            description="How many temperature reduction steps have to be performed per cycle",&
            usage="REDUCTION_STEPS 5", default_i_val=5, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="REDUCTION_FACTOR",&
            description="How to scale the temperature at each reduction step",&
            usage="REDUCTION_FACTOR 0.8", default_r_val=0.8_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="TEMPERATURE",&
            description="Highest temperature to be used (in units of the maximum cost of a non-zero 2D block)",&
            usage="TEMPERATURE 0.5", default_r_val=0.5_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="TOLERANCE",&
            description="If the found distribution is within a fraction of the theoretical optimal one, the optimization stops",&
            usage="TOLERANCE 0.01", default_r_val=0.02_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="SWAP_PROBABILITY",&
            description="Probability of a col/row assignment swap",&
            usage="SWAP_PROBABILITY 0.8", default_r_val=0.9_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="OPT_TARGET",&
            description="The aim of the MC optimization ",&
            usage="OPT_TARGET PERFECT_BALANCE",&
            enum_c_vals=s2a("SMALLEST_MAXIMUM","PERFECT_BALANCE"),&
            enum_i_vals=(/opt_target_minimum, opt_target_median/),&
            enum_desc=s2a("minimize the maximum cost",&
                          "minimize the deviance from the average"), &
            default_i_val=opt_target_minimum,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       NULLIFY(print_key)
       CALL cp_print_key_section_create(print_key,"PROGRAM_RUN_INFO",&
            description="Controls the printing of basic information during the run", &
            print_level=low_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section, print_key, error=error)
       CALL section_release(print_key,error=error)

    END IF

  END SUBROUTINE create_2d_mc_section

END MODULE input_cp2k_distribution
