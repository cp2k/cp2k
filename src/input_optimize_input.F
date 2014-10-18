!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief builds the input structure for optimize_input
!> \par History
!>      09.2010 created [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
MODULE input_optimize_input
  USE cp_output_handling,              ONLY: cp_print_key_section_create,&
                                             low_print_level
  USE f77_blas
  USE input_constants,                 ONLY: opt_force_matching
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_add_subsection,&
                                             section_create,&
                                             section_release,&
                                             section_type
  USE input_val_types,                 ONLY: char_t,&
                                             real_t
  USE kinds,                           ONLY: dp
  USE string_utilities,                ONLY: s2a
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_optimize_input'
  PUBLIC :: create_optimize_input_section

CONTAINS

! *****************************************************************************
!> \brief creates the optimize_input section
!> \param section ...
!> \param error ...
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE create_optimize_input_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'create_optimize_input_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: sub_section, subsubsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    CALL section_create(section,name="OPTIMIZE_INPUT",&
         description="describes an input optimization job, in which parameters in input files get optimized.",&
         repeats=.FALSE., required=.TRUE., supported_feature=.TRUE.,error=error)
    NULLIFY(keyword)

    CALL keyword_create(keyword, name="METHOD",&
         description="What kind of input optimization to perform.",&
         usage="METHOD FORCE_MATCHING",&
         enum_c_vals=s2a("FORCE_MATCHING"),&
         enum_desc=s2a("Perform a force matching minimization."),&
         enum_i_vals=(/opt_force_matching/), &
         default_i_val=opt_force_matching,&
         required=.TRUE.,&
         error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="ACCURACY",&
         description="Final accuracy requested in optimization (RHOEND)",&
         usage="ACCURACY 0.00001",&
         default_r_val=1.e-5_dp,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="STEP_SIZE",&
         description="Initial step size for search algorithm (RHOBEG)",&
         usage="STEP_SIZE 0.005",&
         default_r_val=0.05_dp,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="MAX_FUN",&
         description="Maximum number of function evaluations",&
         usage="MAX_FUN 1000",&
         default_i_val=5000,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="ITER_START_VAL",&
         description="Used for restarting, starting value of the iteration",&
         usage="ITER_START_VAL 0",&
         default_i_val=0,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    CALL keyword_create(keyword, name="RANDOMIZE_VARIABLES",&
         description="Percentage randomization of the free variables applied initially",&
         usage="RANDOMIZE_VARIABLES 20",&
         default_r_val=0.00_dp,error=error)
    CALL section_add_keyword(section,keyword,error=error)
    CALL keyword_release(keyword,error=error)

    !
    ! variables section
    !

    NULLIFY(sub_section)
    CALL section_create(sub_section,name="VARIABLE",&
         description="Defines initial values for variables and their labels",&
         n_subsections=0, repeats=.TRUE., required=.TRUE.,&
         error=error)

       CALL keyword_create(keyword, name="VALUE",&
            description="Initial value of the variable",&
            usage="VALUE 0.0",&
            required=.TRUE.,&
            type_of_var=real_t, unit_str="internal_cp2k",error=error)
       CALL section_add_keyword(sub_section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="FIXED",&
            description="Is this variable fixed or should it be optimized.",&
            usage="FIXED",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(sub_section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="LABEL",&
            description="The label used in the input file, i.e. ${LABEL} will be replaced by the VALUE specified.",&
            usage="LABEL PRM01",&
            required=.TRUE.,&
            type_of_var=char_t,&
            error=error)
       CALL section_add_keyword(sub_section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    CALL section_add_subsection(section, sub_section, error=error)
    CALL section_release(sub_section,error=error)


    !
    ! force matching sub sectiong
    !

    NULLIFY(sub_section)
    CALL section_create(sub_section,name="FORCE_MATCHING",&
         description="Specify the force matching input.",&
         repeats=.TRUE., required=.TRUE., supported_feature=.TRUE.,error=error)

         CALL keyword_create(keyword,name="OPTIMIZE_FILE_NAME",&
              description="the filename of the input file which contains the parameters to be optimized",&
              usage="OPTIMIZE_FILE_NAME my_input.inp",&
              default_lc_val="",&
              required=.TRUE.,error=error)
         CALL section_add_keyword(sub_section,keyword,error=error)
         CALL keyword_release(keyword,error=error)

         CALL keyword_create(keyword,name="REF_TRAJ_FILE_NAME",&
              description="the filename of the reference coordinates.",&
              usage="REF_TRAJ_FILE_NAME pos.xyz",&
              default_lc_val="",&
              error=error)
         CALL section_add_keyword(sub_section,keyword,error=error)
         CALL keyword_release(keyword,error=error)

         CALL keyword_create(keyword,name="REF_FORCE_FILE_NAME",&
              description="the filename of the reference forces, should also contain the energy",&
              usage="REF_FORCE_FILE_NAME frc.xyz",&
              default_lc_val="",&
              error=error)
         CALL section_add_keyword(sub_section,keyword,error=error)
         CALL keyword_release(keyword,error=error)

         CALL keyword_create(keyword,name="REF_CELL_FILE_NAME",&
              description="the filename of the reference cell",&
              usage="REF_CELL_FILE_NAME project.cell",&
              default_lc_val="",&
              error=error)
         CALL section_add_keyword(sub_section,keyword,error=error)
         CALL keyword_release(keyword,error=error)

         CALL keyword_create(keyword, name="GROUP_SIZE",&
              description="Gives the preferred size of a working group, "//&
              "groups will always be equal or larger than this size."//&
              "Usually this should take the number of cores per socket into account for good performance.",&
              usage="group_size 2", default_i_val=6, supported_feature=.TRUE.,error=error)
         CALL section_add_keyword(sub_section,keyword,error=error)
         CALL keyword_release(keyword,error=error)

         CALL keyword_create(keyword, name="FRAME_START",&
              description="starting frame to be used from the reference trajectory",&
              usage="FRAME_START 1", default_i_val=1, error=error)
         CALL section_add_keyword(sub_section,keyword,error=error)
         CALL keyword_release(keyword,error=error)

         CALL keyword_create(keyword, name="FRAME_STOP",&
              description="final frame to be used from the reference trajectory (all=-1)",&
              usage="FRAME_STOP -1", default_i_val=-1, error=error)
         CALL section_add_keyword(sub_section,keyword,error=error)
         CALL keyword_release(keyword,error=error)

         CALL keyword_create(keyword, name="FRAME_STRIDE",&
              description="stride when using the reference trajectory",&
              usage="FRAME_STRIDE 1", default_i_val=1, error=error)
         CALL section_add_keyword(sub_section,keyword,error=error)
         CALL keyword_release(keyword,error=error)

         CALL keyword_create(keyword, name="FRAME_COUNT",&
              description="Use at most FRAME_COUNT frames from the reference trajectory, "//&
                          "adjusting the stride to have them as fas apart as possible (all=-1).",&
              usage="FRAME_COUNT 100", default_i_val=-1, error=error)
         CALL section_add_keyword(sub_section,keyword,error=error)
         CALL keyword_release(keyword,error=error)

         CALL keyword_create(keyword, name="ENERGY_WEIGHT",&
              description="Relative weight of the energy RMSD vs the force RMSD",&
              usage="ENERGY_WEIGHT 0.1", default_r_val=0.1_dp, error=error)
         CALL section_add_keyword(sub_section,keyword,error=error)
         CALL keyword_release(keyword,error=error)

         CALL keyword_create(keyword, name="SHIFT_AVERAGE",&
              description="Shift averages of the energies before computing energy RMSD.",&
              usage="SHIFT_AVERAGE", default_l_val=.FALSE.,lone_keyword_l_val=.TRUE., error=error)
         CALL section_add_keyword(sub_section,keyword,error=error)
         CALL keyword_release(keyword,error=error)

         CALL keyword_create(keyword, name="SHIFT_QM",&
              description="Shift of the reference energies applied before computing energy RMSD.",&
              usage="SHIFT_QM -17.0", default_r_val=0.0_dp, error=error)
         CALL section_add_keyword(sub_section,keyword,error=error)
         CALL keyword_release(keyword,error=error)

         CALL keyword_create(keyword, name="SHIFT_MM",&
              description="Shift of the fit energies applied before computing energy RMSD.",&
              usage="SHIFT_MM 0.0", default_r_val=0.0_dp, error=error)
         CALL section_add_keyword(sub_section,keyword,error=error)
         CALL keyword_release(keyword,error=error)

         NULLIFY(subsubsection)
         CALL cp_print_key_section_create(subsubsection,"COMPARE_ENERGIES",&
              description="A comparison of energies between fit and reference",&
              print_level=low_print_level,filename="compare_energies",common_iter_levels=1,&
              error=error)
         CALL section_add_subsection(sub_section,subsubsection,error=error)
         CALL section_release(subsubsection,error=error)

         NULLIFY(subsubsection)
         CALL cp_print_key_section_create(subsubsection,"COMPARE_FORCES",&
              description="A comparison of forces between fit and reference",&
              print_level=low_print_level,filename="compare_forces",common_iter_levels=1,&
              error=error)
         CALL section_add_subsection(sub_section,subsubsection,error=error)
         CALL section_release(subsubsection,error=error)

    CALL section_add_subsection(section, sub_section, error=error)
    CALL section_release(sub_section,error=error)

    NULLIFY(subsubsection)
    CALL cp_print_key_section_create(subsubsection,"HISTORY",&
         description="writes a history of the function value and parameters",&
         print_level=low_print_level,filename="history",common_iter_levels=1,&
         error=error)
    CALL section_add_subsection(section,subsubsection,error=error)
    CALL section_release(subsubsection,error=error)

    CALL cp_print_key_section_create(subsubsection,"RESTART",&
         description="writes an input file that can be used to restart ",&
         print_level=low_print_level,filename="optimize",common_iter_levels=1,&
         error=error)
    CALL keyword_create(keyword, name="BACKUP_COPIES",&
          description="Specifies the maximum index of backup copies.",&
          usage="BACKUP_COPIES {int}",&
          default_i_val=3, error=error)
    CALL section_add_keyword(subsubsection,keyword,error=error)
    CALL keyword_release(keyword,error=error)
    CALL section_add_subsection(section,subsubsection,error=error)
    CALL section_release(subsubsection,error=error)


  END SUBROUTINE create_optimize_input_section

END MODULE input_optimize_input

