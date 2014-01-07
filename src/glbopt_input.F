!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief performs global geometry optimization
!> \par History
!> \author Ole
! *****************************************************************************
MODULE glbopt_input
  USE bibliography,                    ONLY: Goedecker2004
  USE cp_output_handling,              ONLY: cp_print_key_section_create
  USE input_constants,                 ONLY: glbopt_do_mincrawl,&
                                             glbopt_do_minhop,&
                                             low_print_level
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
#include "cp_common_uses.h"

 IMPLICIT NONE
 PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'glbopt_input'

 PUBLIC :: glbopt_declare_input

  CONTAINS



!! *****************************************************************************
!! *****************************************************************************
  SUBROUTINE glbopt_declare_input(swarm_section,error)
    TYPE(section_type), POINTER              :: swarm_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: glbopt_section, printkey

    NULLIFY(glbopt_section, keyword, printkey)

    CALL section_create(glbopt_section,name="GLOBAL_OPT",&
         description="Section to control global geometry optimizations.",&
         repeats=.FALSE., supported_feature=.TRUE.,error=error)

    CALL keyword_create(keyword, name="METHOD",&
         description="Methods to use for optimization.", required=.TRUE.,&
         default_i_val=glbopt_do_minhop,&
         enum_c_vals=s2a("MINIMA_HOPPING","MINIMA_CRAWLING"),&
         enum_desc=s2a("Runs Minima-Hopping algorithm.",&
                       "Runs Minima-Crawling algorithm."),&
         enum_i_vals=(/glbopt_do_minhop, glbopt_do_mincrawl/),&
         error=error)
    CALL section_add_keyword(glbopt_section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="E_TARGET",&
         description="Target Energy, the optimization will quit once a lower potential energy is reached.",&
         default_r_val=-1*HUGE(1.0_dp),type_of_var=real_t,unit_str="hartree",error=error)
    CALL section_add_keyword(glbopt_section, keyword, error=error)
    CALL keyword_release(keyword, error=error)

    CALL keyword_create(keyword, name="MD_BUMPS_MAX",&
         description="Number of bumps in potential energy after which MD runs ends.",&
         type_of_var=integer_t,default_i_val=3,error=error)
    CALL section_add_keyword(glbopt_section, keyword, error=error)
    CALL keyword_release(keyword, error=error)

    CALL keyword_create(keyword, name="BUMP_STEPS_UPWARDS",&
         description="Number of MD steps with potential energy increases required for a bump.",&
         type_of_var=integer_t,default_i_val=2,error=error)
    CALL section_add_keyword(glbopt_section, keyword, error=error)
    CALL keyword_release(keyword, error=error)

    CALL keyword_create(keyword, name="BUMP_STEPS_DOWNWARDS",&
         description="Number of MD steps with potential energy decreases required for a bump.",&
         type_of_var=integer_t,default_i_val=2,error=error)
    CALL section_add_keyword(glbopt_section, keyword, error=error)
    CALL keyword_release(keyword, error=error)

    CALL keyword_create(keyword, name="FRAGMENTATION_THRESHOLD",&
         description="Threshold for atom distance used for detecting fragmentation of clusters.", &
         default_r_val=2.0_dp, unit_str="angstrom",type_of_var=real_t,error=error)
    CALL section_add_keyword(glbopt_section, keyword, error=error)
    CALL keyword_release(keyword, error=error)

    !CALL keyword_create(keyword, name="MD_ADAPTIVE_TIMESTEP",&
    !     description="Make MD timesteps longer for lower temperatures.", &
    !     default_r_val=0.0_dp, type_of_var=real_t,error=error)
    !CALL section_add_keyword(glbopt_section, keyword, error=error)
    !CALL keyword_release(keyword, error=error)

    CALL cp_print_key_section_create(printkey,"PROGRESS_TRAJECTORY",&
            description="Printkey to controll the writting of the progress trajectory. "//&
            "This trajectory contains the minima, which are lower in energy than the by then lowerest.",&
            print_level=low_print_level, common_iter_levels=1,&
            filename="",unit_str="angstrom",error=error)
    CALL section_add_subsection(glbopt_section,printkey,error=error)
    CALL section_release(printkey,error=error)

    CALL history_declare_input(glbopt_section, error=error)
    CALL minhop_declare_input(glbopt_section, error=error)
    CALL mincrawl_declare_input(glbopt_section, error=error)

    CALL section_add_subsection(swarm_section, glbopt_section,error=error)
    CALL section_release(glbopt_section,error=error)
  END SUBROUTINE glbopt_declare_input

!! *****************************************************************************
!! *****************************************************************************
  SUBROUTINE history_declare_input(glbopt_section,error)
    TYPE(section_type), POINTER              :: glbopt_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: history_section

    NULLIFY(history_section, keyword)

    CALL section_create(history_section,name="HISTORY",&
         description="Section controling the history of visited minima and "//&
         "how minima are recognized at a later point.",&
         repeats=.FALSE., supported_feature=.TRUE.,error=error)

    CALL keyword_create(keyword, name="ENERGY_PRECISION",&
         description="If the difference of two energies is below this threshold "//&
         "they are considert equal.", &
         default_r_val=1.0e-5_dp, type_of_var=real_t,error=error)
    CALL section_add_keyword(history_section, keyword, error=error)
    CALL keyword_release(keyword, error=error)

    CALL keyword_create(keyword, name="FINGERPRINT_PRECISION",&
         description="If the euclidean distance of two fingerprints is below "//&
         "this threshold, they are considert equal.", &
         default_r_val=1.0e-2_dp, type_of_var=real_t,error=error)
    CALL section_add_keyword(history_section, keyword, error=error)
    CALL keyword_release(keyword, error=error)

    CALL section_add_subsection(glbopt_section,history_section,error=error)
    CALL section_release(history_section,error=error)
  END SUBROUTINE history_declare_input


!! *****************************************************************************
!! *****************************************************************************
  SUBROUTINE minhop_declare_input(glbopt_section,error)
    TYPE(section_type), POINTER              :: glbopt_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: minhop_section

    NULLIFY(minhop_section, keyword)

    CALL section_create(minhop_section,name="MINIMA_HOPPING",&
         description="Section controlling the Minima Hopping method.",&
         citations=(/Goedecker2004/),&
         repeats=.FALSE., supported_feature=.TRUE.,error=error)

    CALL keyword_create(keyword, name="BETA_1",&
         description="Factor used to increase temperature when escape failed, "//&
         "should be larger than 1.", &
         default_r_val=1.1_dp, type_of_var=real_t,error=error)
    CALL section_add_keyword(minhop_section, keyword, error=error)
    CALL keyword_release(keyword, error=error)

    CALL keyword_create(keyword, name="BETA_2",&
         description="Factor used to increase temperature when escape found "//&
         "known minima, should be larger than 1.", &
         default_r_val=1.1_dp, type_of_var=real_t,error=error)
    CALL section_add_keyword(minhop_section, keyword, error=error)
    CALL keyword_release(keyword, error=error)

    CALL keyword_create(keyword, name="BETA_3",&
         description="Factor used to decrease temperature when escape succeeded"//&
         ", should be smaller than 1.", &
         default_r_val=1.0/1.1_dp, type_of_var=real_t,error=error)
    CALL section_add_keyword(minhop_section, keyword, error=error)
    CALL keyword_release(keyword, error=error)

    CALL keyword_create(keyword, name="ALPHA_1",&
         description="Factor used to decrease acceptance energy, when minima was accepted"//&
         ", should be smaller than 1.", &
         default_r_val=0.98_dp,type_of_var=real_t,error=error)
    CALL section_add_keyword(minhop_section, keyword, error=error)
    CALL keyword_release(keyword, error=error)

    CALL keyword_create(keyword, name="ALPHA_2",&
         description="Factor used to increase acceptance energy, when minima was rejected"//&
         ", should be larger than 1.", &
         default_r_val=1.0/0.98_dp,type_of_var=real_t,error=error)
    CALL section_add_keyword(minhop_section, keyword, error=error)
    CALL keyword_release(keyword, error=error)

    CALL keyword_create(keyword, name="E_ACCEPT_INIT",&
         description="Initial value of acceptance energy",  &
         default_r_val=0.005_dp, type_of_var=real_t,unit_str="hartree",error=error)
    CALL section_add_keyword(minhop_section, keyword, error=error)
    CALL keyword_release(keyword, error=error)

    CALL keyword_create(keyword, name="TEMPERATURE_INIT",&
         description="Initially temperature in Kelvin",  &
         default_r_val=100.0_dp,type_of_var=real_t,error=error)
    CALL section_add_keyword(minhop_section, keyword, error=error)
    CALL keyword_release(keyword, error=error)

    CALL keyword_create(keyword, name="SHARE_HISTORY",&
         description="If set all worker will use a single share history of visited minima.",&
         default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,supported_feature=.TRUE.,error=error)
    CALL section_add_keyword(minhop_section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL section_add_subsection(glbopt_section,minhop_section,error=error)
    CALL section_release(minhop_section,error=error)
  END SUBROUTINE minhop_declare_input


!! *****************************************************************************
!! *****************************************************************************
  SUBROUTINE mincrawl_declare_input(glbopt_section,error)
    TYPE(section_type), POINTER              :: glbopt_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: mincrawl_section, printkey

    NULLIFY(mincrawl_section, keyword, printkey)

    CALL section_create(mincrawl_section,name="MINIMA_CRAWLING",&
         description="Section controls Minima Crawling run.",&
         repeats=.FALSE., supported_feature=.TRUE.,error=error)

    CALL keyword_create(keyword, name="TEMPSTEP_BASE",&
         description="Base used to calculate temperature steps base**n", &
         default_r_val=1.1_dp, type_of_var=real_t,error=error)
    CALL section_add_keyword(mincrawl_section, keyword, error=error)
    CALL keyword_release(keyword, error=error)

    CALL keyword_create(keyword, name="TEMPSTEP_MAX",&
         description="Maximum number of temperature steps.", &
         default_i_val=100, type_of_var=integer_t,error=error)
    CALL section_add_keyword(mincrawl_section, keyword, error=error)
    CALL keyword_release(keyword, error=error)

    CALL keyword_create(keyword, name="TEMPDIST_UPDATE_WIDTH",&
         description="Width of gaussian used to update temperature distribution.", &
         default_r_val=2.0_dp, type_of_var=real_t,error=error)
    CALL section_add_keyword(mincrawl_section, keyword, error=error)
    CALL keyword_release(keyword, error=error)

    CALL keyword_create(keyword, name="TEMPDIST_UPDATE_HEIGHT",&
         description="Height of gaussian used to update temperature distribution.", &
         default_r_val=0.1_dp, type_of_var=real_t,error=error)
    CALL section_add_keyword(mincrawl_section, keyword, error=error)
    CALL keyword_release(keyword, error=error)

    CALL keyword_create(keyword, name="TEMPERATURE_INIT",&
         description="Initial temperature in Kelvin",  &
         default_r_val=100.0_dp,type_of_var=real_t,error=error)
    CALL section_add_keyword(mincrawl_section, keyword, error=error)
    CALL keyword_release(keyword, error=error)

    CALL keyword_create(keyword, name="TEMPDIST_INIT_WIDTH",&
         description="Initial width of temperature distribution.", &
         default_r_val=5.0_dp, type_of_var=real_t,error=error)
    CALL section_add_keyword(mincrawl_section, keyword, error=error)
    CALL keyword_release(keyword, error=error)

    CALL keyword_create(keyword, name="WORKER_PER_MINIMA",&
         description="Maximum number of active workers per Minima.", &
         default_i_val=3, type_of_var=integer_t,error=error)
    CALL section_add_keyword(mincrawl_section, keyword, error=error)
    CALL keyword_release(keyword, error=error)

    CALL keyword_create(keyword, name="ESCAPE_HISTORY_LENGTH",&
         description="Number of escapes averaged for scoring of minima.", &
         default_i_val=10, type_of_var=integer_t,error=error)
    CALL section_add_keyword(mincrawl_section, keyword, error=error)
    CALL keyword_release(keyword, error=error)

    CALL cp_print_key_section_create(printkey,"MINIMA_TRAJECTORY",&
            description="Printkey to controll the writting of the minima trajectory. "//&
            "This trajectory contains all encountered local minima.",&
            print_level=low_print_level, common_iter_levels=1,&
            filename="",unit_str="angstrom",error=error)
    CALL section_add_subsection(mincrawl_section,printkey,error=error)
    CALL section_release(printkey,error=error)

    CALL section_add_subsection(glbopt_section,mincrawl_section,error=error)
    CALL section_release(mincrawl_section,error=error)
  END SUBROUTINE mincrawl_declare_input

END MODULE glbopt_input

