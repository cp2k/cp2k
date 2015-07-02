!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief function that builds the distribution section of the input
!> \par History
!>      04.2007 created
!> \author Joost VandeVondele
! *****************************************************************************
MODULE input_cp2k_distribution
  
  USE input_constants,                 ONLY: model_block_count,&
                                             model_block_lmax,&
                                             model_block_surface
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_create,&
                                             section_type
  USE string_utilities,                ONLY: s2a
#include "./common/cp_common_uses.f90"

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

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="DISTRIBUTION",&
            description="can be used used to tune the parallel distribution of the data",&
            n_keywords=2, n_subsections=2, repeats=.FALSE., &
            error=error)

       NULLIFY(keyword)

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

       CALL keyword_create(keyword, name="BASIC_CLUSTER_OPTIMIZATION",&
            description="Creates a distribution with spatial info, using recursively KMEANS clustering. ",&
            usage="BASIC_CLUSTER_OPTIMIZATION TRUE",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SYMMETRIC",&
            description="Take the symmetry of the distribution_2d into account.",&
            usage="SYMMETRIC TRUE",&
            default_l_val=.TRUE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_distribution_section

END MODULE input_cp2k_distribution
