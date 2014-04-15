!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      06.2013 split input_cp2k into smaller modules 
!> \author Mandes
! *****************************************************************************
MODULE input_cp2k_tmc
  USE input_cp2k_mm,                   ONLY: create_CHARGE_section
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_add_subsection,&
                                             section_create,&
                                             section_release,&
                                             section_type
  USE input_val_types,                 ONLY: char_t,&
                                             integer_t,&
                                             real_t
  USE kinds,                           ONLY: dp
  USE tmc_stati,                       ONLY: tmc_default_dot_file_name,&
                                             tmc_default_unspecified_name
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_tmc'

  PUBLIC :: create_TMC_section

CONTAINS

! *****************************************************************************
!> \brief creates the TreeMonteCarlo subsection
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes
! *****************************************************************************
  SUBROUTINE create_TMC_section(section, error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_TMC_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    IF (.NOT. failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="TMC",&
            description="A parallelized MC algorithm, presampling the configurations, espacially the Markov chain elements",&
            n_keywords=1, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword, subsection)

       ! size of a group calculating the exact energy 
       CALL keyword_create(keyword=keyword,&
            name="GROUP_ENERGY_SIZE",&
            description="amount of CPUs per group for energy calculation",&
            usage="GROUP_ENERGY_SIZE {INTEGER}",&
            default_i_val=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! size of a group calculating the exact energy 
       CALL keyword_create(keyword=keyword,&
            name="GROUP_ENERGY_NR",&
            description="amount of groups for exact energy calculation",&
            usage="GROUP_ENERGY_NR {INTEGER}",&
            default_i_val=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! size of a group calculating the exact energy 
       CALL keyword_create(keyword=keyword,&
            name="GROUP_CC_SIZE",&
            description="amount of of CPUs per group for configurational change",&
            usage="GROUP_CC_SIZE {INTEGER}",&
            default_i_val=0, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! size of a group calculating the exact energy 
       CALL keyword_create(keyword=keyword,&
            name="GROUP_ANLYSIS_NR",&
            description="amount of groups (cores) for analysing the configurations",&
            usage="GROUP_ANALYSIS_NR {INTEGER}",&
            default_i_val=1, lone_keyword_i_val=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,&
            name="NUM_MC_ELEM",&
            description="the minimum Markov Chain elements, to be sampled",&
            usage="NUM_MC_ELEM {INTEGER}",&
            default_i_val=-1, required=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! the start value for the random number generator 
       CALL keyword_create(keyword=keyword,&
            name="RND_DETERMINISTIC",&
            description="the initialisation number for the random number generator",&
            usage="RND_INIT {INTEGER}",&
            default_i_val=-1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,&
            name="TASK_TYPE",&
            description="Select specialized types. Selectable:"//&
                        "IDEAL_GAS (constant configuration energy E=0.0), ",&
            usage="TASK_TYPE {OPTION}",&
            default_c_val="", lone_keyword_c_val=tmc_default_unspecified_name,&
            required=.FALSE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,&
            name="NR_TEMPERATURE",&
            description="the number of different temperature for parallel tempering",&
            usage="NR_TEMP {INTEGER}",&
            default_i_val=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, &
            name="TEMPERATURE",&
            description="one temperatue OR for parallel tempering: "//&
                        "Tmin Tmax or each temperature T1 T2 T3 ..."//&
                        " If every single temperature is specified, "//&
                        "do NOT use keyword NR_TEMPERATURE",&
            usage="TEMPERATURE {REAL} |OR| TEMPERATURE {REAL} {REAL} ...",&
            default_r_vals=(/330.0_dp/), &
            n_var=-1,type_of_var=real_t,required=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,&
            name="NUM_MV_ELEM_IN_CELL",&
            description="the number of elements (atoms or molecules) "//&
                        "moves in cell or sub box."//&
                        "if 0 all elements are moved once in a MC move",&
            usage="NUM_MV_ELEM_IN_CELL {INTEGER}",&
            default_i_val=0, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

!       CALL keyword_create(keyword=keyword,&
!            name="NR_NMC_STEPS",&
!            description="the number of Nested Mont Carlo moves with in one MC move "//&
!                        "should be huge enough to reach euilibrium state", &
!            usage="NR_NMC_STEPS {INTEGER}",&
!            default_i_val=-1, error=error)
!       CALL section_add_keyword(section,keyword,error=error)
!       CALL keyword_release(keyword,error=error)

       ! the moves MOVE_TYPE on exact potential
       CALL create_TMC_move_type_section(subsection, error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       ! the moves MOVE_TYPE on approx potential
       CALL create_TMC_NMC_move_type_section(subsection, error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL keyword_create(keyword, &
            name="SUB_BOX",&
            description="specifies the size ot the sub box. "//&
            "Standard moves only within subbox of random position, "//&
            "to compensate the potential difference of the approximate potential.", &
            usage="SUB_BOX {REAL} {REAL} {REAL} OR SUB_BOX {REAL} for cubic",&
            default_r_vals=(/-1.0_dp/), &
            n_var=-1,type_of_var=real_t,required=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,&
            name="PRESSURE",&
            description="enables NPT calculation with specified constant pressure [bar]", &
            usage="PRESSURE {REAL}",&
            default_r_val=-1.0_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,&
            name="VOLUME_ISOTROPIC",&
            description="volume move is equal in each direction", &
            usage="VOLUME_ISOTROPIC {LOGICAL}",&
            default_l_val=.TRUE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,&
            name="MOVE_CENTER_OF_MASS",&
            description="Moves the center of mass of defined molecules (in volume moves)", &
            usage="MOVE_CENTER_OF_MASS {LOGICAL}",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, &
            name="ESIMATE_ACC_PROB",&
            description="set the estimation of the acceptance probability using run time information of the energy", &
            usage="ESIMATE_ACC_PROB {LOGICAL}",&
            default_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,&
            name="SPECULATIVE_CANCELING",&
            description="enables or disables the speculative canceling. If we have further knowledge of "//&
                        "acceptance probabilities using parent acceptance or the estimated energy.",&
            usage="SPECULATIVE_CANCELING {LOGICAL}",&
            default_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,&
            name="USE_SCF_ENERGY_INFO",&
            description="enables or disables the usage of SCF energy information for "//&
                        " estimating the acceptance probability. ",&
            usage="USE_SCF_ENERGY_INFO {LOGICAL}",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,&
            name="RESULT_LIST_IN_MEMORY",&
            description="enables the storing of the whole Markov Chain", &
            usage="RESULT_LIST_IN_MEMORY {LOGICAL}",&
            default_l_val=.FALSE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,&
            name="INFO_OUT_STEP_SIZE",&
            description="the number the amount of calculated configurations between to output printings.",&
            usage="INFO_OUT_STEP_SIZE {INTEGER}",&
            default_i_val=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,&
            name="RESTART_IN",&
            description="if existing use the last restart file", &
            usage="RESTART or RSTART {FILENAME}",&
            default_c_val="", lone_keyword_c_val=tmc_default_unspecified_name,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,&
            name="RESTART_OUT",&
            description="Defines the frequency to write restart files. "//&
                        "If no frequency is specified (lone keyword) "//&
                        "the restart file is written at the end (only). "//&
                        "If the value is 0, no restart file is written at all. "//&
                        "The frequency specifies is related "//&
                        "to the calculated Markov chain elements",&
            usage="RESTART or RESTART {INTEGER}",&
            default_i_val=-1, lone_keyword_i_val=-9,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ENERGY_FILE_NAME",&
            description="input file name for the exact potential energy calculation.",&
            usage="ENERGY_FILE_NAME {filename}",&
            default_c_val="", required=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

!       CALL keyword_create(keyword, name="NMC_FILE_NAME",&
!            description="input file name for the approximate potential for Nested Monte Carlo.",&
!            usage="NMC_FILE_NAME {filename}",&
!            default_c_val="",lone_keyword_c_val=tmc_default_unspecified_name,&
!            required=.FALSE., error=error)
!       CALL section_add_keyword(section,keyword,error=error)
!       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,&
            name="PRINT_ONLY_ACC",&
            description="printing only accepted elements of the Markov Chain.", &
            usage="PRINT_ONLY_ACC {LOGICAL}",&
            default_l_val=.TRUE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,&
            name="PRINT_COORDS",&
            description="printing coordinates of the Markov Chain elements", &
            usage="PRINT_COORDS {LOGICAL}",&
            default_l_val=.TRUE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,&
            name="PRINT_FORCES",&
            description="printing forces of the Markov Chain elements", &
            usage="PRINT_FORCES {LOGICAL}",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,&
            name="PRINT_DIPOLE",&
            description="enables the calculation and printing the exact cell dipoles"//&
            " (only for QS methods)", &
            usage="PRINT_DIPOLE {LOGICAL}",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,&
            name="PRINT_CELL",&
            description="printing the cell vectors of the Markov Chain elements", &
            usage="PRINT_CELL {LOGICAL}",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,&
            name="PRINT_ENERGIES",&
            description="printing the different calculated energies (approximated, scf and exact)", &
            usage="PRINT_ENERGIES {LOGICAL}",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DOT_TREE",&
            description="file name for GrapgViz dot file",&
            usage="DOT_TREE {filename}",&
            default_c_val="", lone_keyword_c_val=tmc_default_dot_file_name,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ALL_CONF_FILE_NAME",&
            description="file name for printing every single calculated configuration (e.g. for fitting).",&
            usage="ALL_CONF_FILE_NAME {filename}",&
            default_lc_val="",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,&
            name="PRINT_TEST_OUTPUT",&
            description="printing different values for regtest comparison", &
            usage="PRINT_TEST_OUTPUT {LOGICAL}",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_TMC_ana_section(subsection, error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_TMC_ana_files(subsection, error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)
    END IF

  END SUBROUTINE create_TMC_section

! *****************************************************************************
!> \brief creates the TMC section to select the move types performed within the
!>        NMC (on approximate potential)
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes
! *****************************************************************************
  SUBROUTINE create_TMC_NMC_move_type_section(section, error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'create_TMC_NMC_move_type_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.
    NULLIFY(subsection)

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
      CALL section_create(section,name="NMC_MOVES",&
            description="This section specifies the TMC move types, "//&
                        "which are performed within the nested Monte Carlo (NMC). "//&
                        "thus using the approximate potential.",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.FALSE.,&
            error=error)

      NULLIFY(keyword)

      CALL keyword_create(keyword=keyword,&
           name="NR_NMC_STEPS",&
           description="the number of Nested Mont Carlo moves with in one MC move "//&
                       "should be huge enough to reach euilibrium state", &
           usage="NR_NMC_STEPS {INTEGER}",&
           default_i_val=-1, error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword, name="NMC_FILE_NAME",&
           description="input file name for the approximate potential for Nested Monte Carlo.",&
           usage="NMC_FILE_NAME {filename}",&
           default_c_val="",lone_keyword_c_val=tmc_default_unspecified_name,&
           required=.FALSE., error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword, name="PROB",&
           description="Defines the probability of the NMC move "//&
                       "(considering the ration between the selected other moves)"//&
                       "the probabilities of the move types in the NMC section "//&
                       "defines only the weight within the NMC steps",&
           usage="PROB {real}", required=.FALSE., type_of_var=real_t,&
           default_r_val=1.0_dp, error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword, name="INIT_ACC_PROB",&
           description="Defines the initial probability of accepting the move. ",&
           usage="INIT_ACC_PROB {real}", required=.FALSE., type_of_var=real_t,&
           default_r_val=0.5_dp, n_var=1, error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      ! the moves types
      CALL create_TMC_move_type_section(subsection, error)
      CALL section_add_subsection(section, subsection, error=error)
      CALL section_release(subsection,error=error)
    END IF
  END SUBROUTINE create_TMC_NMC_move_type_section

! *****************************************************************************
!> \brief creates the TMC section to select the move types
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes
! *****************************************************************************
  SUBROUTINE create_TMC_move_type_section(section, error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_TMC_move_type_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="MOVE_TYPE",&
            description="This section specifies the TMC move type, "//&
                        "and its properties. "//&
                        "Selectable types are: "//&
                        "ATOM_TRANS atom translation, "//&
                        "MOL_TRANS molecule translation, "//&
                        "MOL_ROT molecule rotation, "//&
                        "PROT_REORDER proton reordering, "//&
                        "PT_SWAP Parallel Tempering swap, "//&
                        "VOL_MOVE volume change, "//&
                        "ATOM_SWAP swaps two atoms of different type.",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="_SECTION_PARAMETERS_",&
            description="The name of the move type described in this section.",&
            usage="ATOM_TRANS", default_c_val="DEFAULT", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SIZE",&
            description="Defines the size of the move:"//&
                        "ATOM_TRANS [A], "//&
                        "MOL_TRANS [A], "//&
                        "MOL_ROT [degree], "//&
                        "PROT_REORDER [], "//&
                        "VOL_MOVE [A], "//&
                        "ATOM_SWAP",&
            usage="SIZE {real}", required=.FALSE., type_of_var=real_t,&
            default_r_val=-1.0_dp, n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PROB",&
            description="Defines the probability of the move "//&
                        "(considering the ration between the selected moves)",&
            usage="PROB {real}", required=.FALSE., type_of_var=real_t,&
            default_r_val=1.0_dp, n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="INIT_ACC_PROB",&
            description="Defines the initial probability of accepting the move. ",&
            usage="INIT_ACC_PROB {real}", required=.FALSE., type_of_var=real_t,&
            default_r_val=0.23_dp, n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ATOMS",&
            description="Defines the atomic kinds involved in the move. "//&
                        "Up to now only used for the atom swap.",&
            usage="ATOMS {KIND1} {KIND2} ... ", required=.FALSE., type_of_var=char_t,&
            n_var=-1, repeats=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
    

  END SUBROUTINE create_TMC_move_type_section

! *****************************************************************************
!> \brief creates the TreeMonteCarlo subsection
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes
! *****************************************************************************
  SUBROUTINE create_TMC_ana_section(section, error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_TMC_ana_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="TMC_ANALYSIS",&
            description="Analysing the Markov Chain elments with the specified methods",&
            n_keywords=1, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            error=error)

       CALL create_TMC_ana_kinds(section=section, error=error)
    END IF

  END SUBROUTINE create_TMC_ana_section

! *****************************************************************************
!> \brief creates the TreeMonteCarlo subsection
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes
! *****************************************************************************
  SUBROUTINE create_TMC_ana_files(section, error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_TMC_ana_files', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    NULLIFY(keyword)
    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="TMC_ANALYSIS_FILES",&
            description="Analysing the Markov Chain elments with the specified methods, "//&
                        "reading form default or specified files",&
            n_keywords=1, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            error=error)

       CALL create_TMC_ana_kinds(section=section, error=error)

       CALL keyword_create(keyword=keyword,&
            name="NR_TEMPERATURE",&
            description="the number of different temperature for parallel tempering",&
            usage="NR_TEMP {INTEGER}",&
            default_i_val=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, &
            name="TEMPERATURE",&
            description="one temperatue OR for parallel tempering: "//&
                        "Tmin Tmax or each temperature T1 T2 T3 ..."//&
                        " If every single temperature is specified, "//&
                        "do NOT use keyword NR_TEMPERATURE",&
            usage="TEMPERATURE {REAL} |OR| TEMPERATURE {REAL} {REAL} ...",&
            default_r_vals=(/330.0_dp/), &
            n_var=-1,type_of_var=real_t,required=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, &
            name="DIRECTORIES",&
            description="Analysing multiple directories, "//&
                        "created by standard parallel MC (e.g. using TMC farming ",&
            usage="DIRECTORIES {DIR1/} {DIR2/} ...",&
            default_c_vals=(/"./"/), &
            n_var=-1,type_of_var=char_t,required=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="FORCE_ENV_FILE",&
            description="input file name for force env, "//&
                        "to get initial configuration with dimensions and cell",&
            usage="FORCE_ENV_FILE {filename}",&
            default_c_val="",lone_keyword_c_val="",&
            required=.FALSE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="POSITION_FILE",&
            description="file name for analysing the position file",&
            usage="POSITION_FILE {filename}",&
            default_c_val="",lone_keyword_c_val="",&
            required=.FALSE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="CELL_FILE",&
            description="file name for analysing the cell file",&
            usage="CELL_FILE {filename}",&
            default_c_val="",lone_keyword_c_val="",&
            required=.FALSE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DIPOLE_FILE",&
            description="file name for analysing the dipole file",&
            usage="DIPOLE_FILE {filename}",&
            default_c_val="",lone_keyword_c_val="",&
            required=.FALSE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,&
            name="START_ELEM",&
            description="start analysis at element with number #",&
            usage="START_ELEM {INTEGER}",&
            default_i_val=-1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword=keyword,&
            name="END_ELEM",&
            description="end analysis at element with number #",&
            usage="END_ELEM {INTEGER}",&
            default_i_val=-1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
    END IF

  END SUBROUTINE create_TMC_ana_files

! *****************************************************************************
!> \brief creates the TreeMonteCarlo subsection
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Mandes
! *****************************************************************************
  SUBROUTINE create_TMC_ana_kinds(section, error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_TMC_ana_kinds', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    NULLIFY(subsection, keyword)

    failure=.FALSE.

    IF (.NOT. failure) THEN
      CALL keyword_create(keyword=keyword,&
           name="RESTART",&
           description="Enables/disables the reading and writing of "//&
                       " analysis restart files ", &
           usage="RESTART {LOGICAL}",&
           default_l_val=.TRUE., lone_keyword_l_val=.TRUE., required=.FALSE., &
           error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword, name="PREFIX_ANA_FILES",&
            description="specifies a prefix for all analysis files.",&
            usage="ANA_FILES_PREFIX {prefix}",&
            default_c_val="", required=.FALSE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword=keyword,&
           name="DENSITY",&
           description="Mass density in the simulation cell, or if specified in sub cubes",&
           usage="DENSITY or DENSITY {INTEGER} {INTEGER} {INTEGER}",&
           default_i_vals=(/1/), required=.FALSE., &
           n_var=-1,type_of_var=integer_t,error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword=keyword,&
           name="G_R",&
           description="Radial Distribution Function for each pair of atoms"//&
                       "using the amount of specified bins within MAX(cell_lenght)/2",&
           usage="G_R or G_R {INTEGER}",&
           default_i_val=-1, lone_keyword_i_val=-1, required=.FALSE., &
           error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword=keyword,&
           name="CLASSICAL_DIPOLE_MOMENTS",&
           description="calculates the classical dipole Moment. "//&
                       "Following flag specifies if they should be written. "//&
                       "Class. Dip. Mom. are also used to unfold the exact dipole moment. ",&
           usage="CLASSICAL_DIPOLE_MOMENTS or CLASSICAL_DIPOLE_MOMENTS {LOGICAL}",&
           default_l_val=.TRUE., lone_keyword_l_val=.TRUE., required=.FALSE., &
           error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      ! for calculating the classical dipole moment we need charges
      CALL create_CHARGE_section(subsection, error=error)
      CALL section_add_subsection(section, subsection, error=error)
      CALL section_release(subsection,error=error)

      CALL keyword_create(keyword=keyword,&
           name="DIPOLE_ANALYSIS",&
           description="Enables analysis of dipole moments, espacially dielectric constant. "//&
                       "An additional type can be specified, e.g. analyzing ice structures."//&
                       "using SYM_XYZ also dipoles (-x,y,z) .. .. (-x,-y,z).... (-x,-y-z) "//&
                       "are regarded, only use it if this configurations have "//&
                       "all the same energy. ", &
           usage="DIPOLE_ANALYSIS or DIPOLE_ANALYSIS {type}",&
           default_c_val="", lone_keyword_c_val=tmc_default_unspecified_name,&
           required=.FALSE., error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword=keyword,&
           name="DEVIATION",&
           description="Calculates the deviation of the position from the last configuration",&
           usage="DEVIATION {LOGICAL}",&
           default_l_val=.FALSE., lone_keyword_l_val=.TRUE., required=.FALSE., &
           error=error)
      CALL section_add_keyword(section,keyword,error=error)
      CALL keyword_release(keyword,error=error)
    END IF
  END SUBROUTINE create_TMC_ana_kinds

END MODULE input_cp2k_tmc
