!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      - taken out of input_cp2k_motion
!> \author teo & fawzi
! *****************************************************************************
MODULE input_cp2k_colvar
  USE bibliography,                    ONLY: Branduardi2007
  USE colvar_types,                    ONLY: &
       do_clv_fix_point, do_clv_geo_center, do_clv_x, do_clv_xy, do_clv_xyz, &
       do_clv_xz, do_clv_y, do_clv_yz, do_clv_z, plane_def_atoms, &
       plane_def_vec
  USE cp_output_handling,              ONLY: add_last_numeric,&
                                             cp_print_key_section_create,&
                                             high_print_level,&
                                             low_print_level
  USE input_constants,                 ONLY: gaussian,&
                                             numerical,&
                                             rmsd_all,&
                                             rmsd_list,&
                                             rmsd_weightlist
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
                                             lchar_t,&
                                             real_t
  USE kinds,                           ONLY: dp
  USE string_utilities,                ONLY: s2a
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.FALSE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_colvar'

PUBLIC :: create_colvar_section,&
          create_colvar_xyz_d_section,&
          create_colvar_xyz_od_section

CONTAINS

! *****************************************************************************
!> \brief creates the colvar section
!> \param section the section to be created
!> \param skip_recursive_colvar ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  RECURSIVE SUBROUTINE create_colvar_section(section,skip_recursive_colvar,error)
    TYPE(section_type), POINTER              :: section
    LOGICAL, OPTIONAL                        :: skip_recursive_colvar
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_colvar_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure, skip
    TYPE(section_type), POINTER              :: print_key, subsection

    failure=.FALSE.
    skip=.FALSE.
    IF(PRESENT(skip_recursive_colvar))skip=skip_recursive_colvar
    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="COLVAR",&
            description="This section specifies the nature of the collective variables.",&
            n_keywords=1, n_subsections=1, repeats=.TRUE., &
            error=error)
       NULLIFY(subsection, print_key)

       CALL create_colvar_var_section(subsection=subsection,&
            section=section,skip_recursive_colvar=skip,error=error)

       CALL section_create(subsection,name="PRINT",&
            description="Controls the printing of the colvar specifications",&
            n_keywords=0, n_subsections=1, repeats=.TRUE., &
            error=error)
       NULLIFY(print_key)
       CALL cp_print_key_section_create(print_key,"PROGRAM_RUN_INFO",&
            description="Controls the printing of basic information during colvar setup.", &
            print_level=low_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_clv_info_section(subsection, error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)
    END IF

  END SUBROUTINE create_colvar_section

! *****************************************************************************
!> \brief Create the restart section for colvar restraints
!>      This section will be only used for restraint restarts.
!>      Constraints are handled automatically
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Teodoro Laino 08.2006
! *****************************************************************************
  SUBROUTINE create_clv_info_section(section, error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_clv_info_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       NULLIFY( keyword)
       CALL section_create(section,name="COLVAR_FUNC_INFO",&
            description="Specify further data possibly used by colvars, depending "//&
            "on the starting geometry, for computing the functions value.",&
            n_subsections=0, repeats=.FALSE., &
            error=error)

       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="Colvar function data."//&
            " The order is an internal order. So if you decide to edit/modify/add these values by hand"//&
            " you should know very well what you are doing.!",repeats=.TRUE.,&
            usage="{Real} ...", type_of_var=real_t, n_var=-1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_clv_info_section

! *****************************************************************************
!> \brief creates the collective variables for the colvar section
!> \param subsection ...
!> \param section the section to be created
!> \param skip_recursive_colvar ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  RECURSIVE SUBROUTINE create_colvar_var_section(subsection, section, skip_recursive_colvar, error)
    TYPE(section_type), POINTER              :: subsection, section
    LOGICAL, INTENT(IN)                      :: skip_recursive_colvar
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_colvar_var_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(subsection),cp_failure_level,routineP,error,failure)
       CPPrecondition(ASSOCIATED(section),cp_failure_level,routineP,error,failure)

       CALL create_colvar_dist_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_colvar_angle_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_colvar_torsion_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_colvar_coord_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_colvar_pop_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_colvar_gyr_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_colvar_d_pl_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_colvar_a_pl_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_colvar_rot_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_colvar_dfunct_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_colvar_qparm_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_colvar_hydronium_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_colvar_rmsd_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_colvar_xyz_d_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_colvar_xyz_od_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_colvar_u_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_colvar_wc_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_colvar_hbp_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_colvar_ring_puckering_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_colvar_cond_dist_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       IF(.NOT. skip_recursive_colvar)THEN
          CALL create_colvar_rpath_section(subsection,error)
          CALL section_add_subsection(section, subsection, error=error)
          CALL section_release(subsection,error=error)

          CALL create_colvar_dpath_section(subsection,error)
          CALL section_add_subsection(section, subsection, error=error)
          CALL section_release(subsection,error=error)

          CALL create_colvar_comb_section(subsection,error)
          CALL section_add_subsection(section, subsection, error=error)
          CALL section_release(subsection,error=error)
       ENDIF

    END IF
  END SUBROUTINE create_colvar_var_section

! *****************************************************************************
!> \brief collective variables specifying coordination
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_colvar_coord_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_colvar_coord_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="coordination",&
            description="Section to define the coordination number as a collective variable.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., &
            error=error)

       NULLIFY(subsection, keyword)

       CALL keyword_create(keyword, name="ATOMS_FROM",&
            variants=(/"POINTS_FROM"/),&
            description="Specify indexes of atoms/points building the coordination variable. ",&
            usage="ATOMS_FROM {integer} {integer} ..",repeats=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ATOMS_TO",&
            variants=(/"POINTS_TO"/),&
            description="Specify indexes of atoms/points building the coordination variable. ",&
            usage="ATOMS_TO {integer} {integer} ..",repeats=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ATOMS_TO_B",&
            variants=(/"POINTS_TO_B"/),&
            description="For the CV given by the multiplication of two coorination numbers,"//&
            " here specify indexes of the third set of atoms/points. ",&
            usage="ATOMS_TO_B {integer} {integer} ..",repeats=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="KINDS_FROM",&
            description="Specify alternatively kinds of atoms building the coordination variable.",&
            usage="KINDS_FROM {CHAR} {CHAR} ..",repeats=.TRUE.,&
            n_var=-1, type_of_var=char_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="KINDS_TO",&
            description="Specify alternatively kinds of atoms building the coordination variable.",&
            usage="KINDS_TO {CHAR} {CHAR} ..",repeats=.TRUE.,&
            n_var=-1, type_of_var=char_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="KINDS_TO_B",&
            description="For the CV given by the multiplication of two coorination numbers,"//&
           " here specify alternatively kinds of atoms building the coordination variable.",&
            usage="KINDS_TO_B {CHAR} {CHAR} ..",repeats=.TRUE.,&
            n_var=-1, type_of_var=char_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! Must be present in each colvar and handled properly
       CALL create_point_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL keyword_create(keyword, name="R0",&
            variants=(/"R_0"/),&
            description="Specify the R0 parameter in the coordination function.",&
            usage="R0 {real}",default_r_val=3.0_dp,&
            unit_str="bohr",n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NN",&
            variants=(/"EXPON_NUMERATOR"/),&
            description="Sets the value of the numerator of the exponential factor"//&
            "in the coordination FUNCTION.",&
            usage="NN {integer}",default_i_val=6,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ND",&
            variants=(/"EXPON_DENOMINATOR"/),&
            description="Sets the value of the denominator of the exponential factor"//&
            "in the coordination FUNCTION.",&
            usage="ND {integer}",default_i_val=12,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="R0_B",&
            variants=(/"R_0_B"/),&
            description="For the CV given by the multiplication of two coorination numbers,"//&
            " specify the R0 parameter in the second coordination function.",&
            usage="R0_B {real}",default_r_val=3.0_dp,&
            unit_str="bohr",n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NN_B",&
            variants=(/"EXPON_NUMERATOR_B"/),&
            description="For the CV given by the multiplication of two coorination numbers,"//&
            "Sets the value of the numerator of the exponential factor"//&
            "in the coordination FUNCTION.",&
            usage="NN_B {integer}",default_i_val=6,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ND_B",&
            variants=(/"EXPON_DENOMINATOR_B"/),&
            description="For the CV given by the multiplication of two coorination numbers,"//&
            "Sets the value of the denominator of the exponential factor"//&
            "in the coordination FUNCTION.",&
            usage="ND_B {integer}",default_i_val=12,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_colvar_coord_section



! *****************************************************************************
!> \brief ...
!> \param section ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE create_colvar_cond_dist_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'create_colvar_cond_dist_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="CONDITIONED_DISTANCE",&
            description="Section to define the conditioned distance as a collective variable.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., &
            error=error)

       NULLIFY(subsection, keyword)

       CALL keyword_create(keyword, name="ATOMS_DISTANCE",&
            description="Specify indexes of atoms/points from which the distance is computed. ",&
            usage="ATOMS_DISTANCE {integer} {integer} ..",repeats=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)


       CALL keyword_create(keyword, name="ATOMS_FROM",&
            variants=(/"POINTS_FROM"/),&
            description="Specify indexes of atoms/points building the coordination variable. ",&
            usage="ATOMS_FROM {integer} {integer} ..",repeats=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ATOMS_TO",&
            variants=(/"POINTS_TO"/),&
            description="Specify indexes of atoms/points building the coordination variable. ",&
            usage="ATOMS_TO {integer} {integer} ..",repeats=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)


       CALL keyword_create(keyword, name="KINDS_FROM",&
            description="Specify alternatively kinds of atoms building the coordination variable.",&
            usage="KINDS_FROM {CHAR} {CHAR} ..",repeats=.TRUE.,&
            n_var=-1, type_of_var=char_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="KINDS_TO",&
            description="Specify alternatively kinds of atoms building the coordination variable.",&
            usage="KINDS_TO {CHAR} {CHAR} ..",repeats=.TRUE.,&
            n_var=-1, type_of_var=char_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! Must be present in each colvar and handled properly
       CALL create_point_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL keyword_create(keyword, name="R0",&
            variants=(/"R_0"/),&
            description="Specify the R0 parameter in the coordination function.",&
            usage="R0 {real}",default_r_val=3.0_dp,&
            unit_str="bohr",n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NN",&
            variants=(/"EXPON_NUMERATOR"/),&
            description="Sets the value of the numerator of the exponential factor"//&
            "in the coordination FUNCTION.",&
            usage="NN {integer}",default_i_val=6,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ND",&
            variants=(/"EXPON_DENOMINATOR"/),&
            description="Sets the value of the denominator of the exponential factor"//&
            "in the coordination FUNCTION.",&
            usage="ND {integer}",default_i_val=12,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="LAMBDA",&
            description="Specify the lambda parameter at the exponent of the conditioned distance function.",&
            usage="R0 {real}",default_r_val=3.0_dp,&
            unit_str="bohr",n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)


    END IF

  END SUBROUTINE create_colvar_cond_dist_section




! *****************************************************************************
!> \brief collective variables specifying population of a specie based on coordination
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \date  01.2009
!> \author Fabio Sterpone
! *****************************************************************************
  SUBROUTINE create_colvar_pop_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_colvar_pop_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="population",&
            description="Section to define the population of specie as a collective variable.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., &
            error=error)

       NULLIFY(subsection, keyword)

       CALL keyword_create(keyword, name="ATOMS_FROM",&
            variants=(/"POINTS_FROM"/),&
            description="Specify indexes of atoms/points building the coordination variable. ",&
            usage="ATOMS_FROM {integer} {integer} ..",repeats=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ATOMS_TO",&
            variants=(/"POINTS_TO"/),&
            description="Specify indexes of atoms/points building the coordination variable. ",&
            usage="ATOMS_TO {integer} {integer} ..",repeats=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="KINDS_FROM",&
            description="Specify alternatively kinds of atoms building the coordination variable.",&
            usage="KINDS_FROM {CHAR} {CHAR} ..",repeats=.TRUE.,&
            n_var=-1, type_of_var=char_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="KINDS_TO",&
            description="Specify alternatively kinds of atoms building the coordination variable.",&
            usage="KINDS_TO {CHAR} {CHAR} ..",repeats=.TRUE.,&
            n_var=-1, type_of_var=char_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! Must be present in each colvar and handled properly
       CALL create_point_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL keyword_create(keyword, name="R0",&
            variants=(/"R_0"/),&
            description="Specify the R0 parameter in the coordination function.",&
            usage="R0 {real}",default_r_val=3.0_dp,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NN",&
            variants=(/"EXPON_NUMERATOR"/),&
            description="Sets the value of the numerator of the exponential factor"//&
            "in the coordination FUNCTION.",&
            usage="NN {integer}",default_i_val=6,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ND",&
            variants=(/"EXPON_DENOMINATOR"/),&
            description="Sets the value of the denominator of the exponential factor"//&
            "in the coordination FUNCTION.",&
            usage="ND {integer}",default_i_val=12,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)


       CALL keyword_create(keyword, name="n0",&
            variants=(/"n_0"/),&
            description="Specify the n0 parameter that sets the coordination of the species.",&
            usage="n0 {integer}",default_i_val=4,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

      CALL keyword_create(keyword, name="SIGMA",&
            description="Specify the gaussian width of used to build the population istogram.",&
            usage="SIGMA {real}",default_r_val=0.5_dp,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
    END IF

  END SUBROUTINE create_colvar_pop_section


! *****************************************************************************
!> \brief ...
!> \param section ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE create_colvar_gyr_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_colvar_gyr_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       CALL section_create(section,name="GYRATION_RADIUS",&
            description="Section to define the gyration radius as a collective variable.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., &
            error=error)

       NULLIFY(subsection, keyword)

       CALL keyword_create(keyword, name="ATOMS",&
            variants=(/"POINTS"/),&
            description="Specify indexes of atoms/points defyining the gyration radius variable. ",&
            usage="ATOMS {integer} {integer} ..",repeats=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="KINDS",&
            description="Specify alternatively kinds of atoms defining the gyration radius.",&
            usage="KINDS {CHAR} {CHAR} ..",repeats=.TRUE.,&
            n_var=-1, type_of_var=char_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! Must be present in each colvar and handled properly
       CALL create_point_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

     END IF

  END SUBROUTINE create_colvar_gyr_section

! *****************************************************************************
!> \brief collective variables specifying torsion
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_colvar_dfunct_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_colvar_dfunct_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="DISTANCE_FUNCTION",&
            description="Section to define functions between two distances as collective variables."//&
            " The function is defined as d1+coeff*d2",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., &
            error=error)

       NULLIFY(keyword, subsection)

       CALL keyword_create(keyword, name="ATOMS",&
            variants=(/"POINTS"/),&
            description="Specifies the indexes of atoms/points for the two bonds d1=(1-2) d2=(3-4).",&
            usage="ATOMS {integer} {integer} {integer} {integer}",&
            n_var=4, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="COEFFICIENT",&
            description="Specifies the coefficient in the function for the constraint."//&
            " -1.0 has to be used for distance difference, 1.0 for distance addition",&
            usage="COEFFICIENT {real}",&
            type_of_var=real_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PBC",&
            description="Whether periodic boundary conditions should be applied on the "//&
            "atomic position before computing the colvar or not.",&
            usage="PBC",&
            default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! Must be present in each colvar and handled properly
       CALL create_point_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_colvar_dfunct_section

! *****************************************************************************
!> \brief collective variables specifying torsion
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_colvar_torsion_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'create_colvar_torsion_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="torsion",&
            description="Section to define the torsion as a collective variables.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., &
            error=error)

       NULLIFY(keyword, subsection)

       CALL keyword_create(keyword, name="ATOMS",&
            variants=(/"POINTS"/),&
            description="Specifies the indexes of atoms/points defining the torsion.",&
            usage="ATOMS {integer} {integer} {integer} {integer}",&
            n_var=4, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! Must be present in each colvar and handled properly
       CALL create_point_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_colvar_torsion_section

! *****************************************************************************
!> \brief collective variables specifying torsion
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_colvar_rot_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_colvar_rot_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="bond_rotation",&
            description="Section to define the rotation of a bond/line with respect to"//&
                 "another bond/line",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., &
            error=error)

       NULLIFY(keyword, subsection)

       CALL keyword_create(keyword, name="P1_BOND1",&
            description="Specifies the index of atom/point defining the first point"//&
            "of the first bond/line.",&
            usage="P1_BOND1 {integer}",&
            n_var=1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="P2_BOND1",&
            description="Specifies the index of atom/point defining the second point"//&
            "of the first bond/line.",&
            usage="P2_BOND1 {integer}",&
            n_var=1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="P1_BOND2",&
            description="Specifies the index of atom/point defining the first point"//&
            "of the second bond/line.",&
            usage="P1_BOND2 {integer}",&
            n_var=1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="P2_BOND2",&
            description="Specifies the index of atom/point defining the second point"//&
            "of the second bond/line.",&
            usage="P2_BOND2 {integer}",&
            n_var=1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

        ! Must be present in each colvar and handled properly
       CALL create_point_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_colvar_rot_section

! *****************************************************************************
!> \brief collective variables specifying angles
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_colvar_angle_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_colvar_angle_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="angle",&
            description="Section to define the angle as a collective variables.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., &
            error=error)
       NULLIFY(keyword, subsection)

       CALL keyword_create(keyword, name="ATOMS",&
            variants=(/"POINTS"/),&
            description="Specifies the indexes of atoms/points defining the angle.",&
            usage="ATOMS {integer} {integer} {integer}",&
            n_var=3, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! Must be present in each colvar and handled properly
       CALL create_point_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_colvar_angle_section

! *****************************************************************************
!> \brief creates the colvar section regarded to the collective variables dist
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_colvar_dist_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_colvar_dist_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="distance",&
            description="Section to define the distance as a collective variables.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., &
            error=error)
       NULLIFY(keyword, subsection)

       CALL keyword_create(keyword, name="ATOMS",&
            variants=(/"POINTS"/),&
            description="Specifies the indexes of atoms/points defining the distance.",&
            usage="ATOMS {integer} {integer}",&
            n_var=2, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="AXIS",&
            description="Define the axes along which the colvar should be evaluated",&
            usage="AXIS (XYZ | X | Y | Z | XY| XZ | YZ)",&
            enum_c_vals=s2a( "XYZ","X", "Y", "Z", "XY", "XZ", "YZ"),&
            enum_i_vals=(/ do_clv_xyz, do_clv_x, do_clv_y,do_clv_z, do_clv_xy, do_clv_xz, do_clv_yz/),&
            default_i_val=do_clv_xyz, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! Must be present in each colvar and handled properly
       CALL create_point_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_colvar_dist_section

! *****************************************************************************
!> \brief creates the colvar section regarded to the collective variables dist
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_colvar_xyz_d_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_colvar_xyz_d_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="XYZ_DIAG",&
            description="Section to define the distance of an atom from its starting "//&
            "position ((X-X(0))^2+(Y-Y(0))^2+(Z-Z(0))^2) or part of its components as a collective variable."//&
            "If absolute_position is specified, instead the CV is represented by the "//&
            "instantaneous position of the atom (only available for X, Y or Z components).",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., &
            error=error)
       NULLIFY(keyword, subsection)

       CALL keyword_create(keyword, name="ATOM",&
            variants=(/"POINT"/),&
            description="Specifies the index of the atom/point.",&
            usage="ATOM {integer}",&
            n_var=1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="COMPONENT",&
            description="Define the component of the position vector which will be used "//&
            "as a colvar.",&
            usage="AXIS (XYZ | X | Y | Z | XY| XZ | YZ)",&
            enum_c_vals=s2a( "XYZ","X", "Y", "Z", "XY", "XZ", "YZ"),&
            enum_i_vals=(/ do_clv_xyz, do_clv_x, do_clv_y,do_clv_z, do_clv_xy, do_clv_xz, do_clv_yz/),&
            default_i_val=do_clv_xyz, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PBC",&
            description="Whether periodic boundary conditions should be applied on the "//&
            "atomic position before computing the colvar or not.",&
            usage="PBC",&
            default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ABSOLUTE_POSITION",&
            description="If enabled, the absolute position of the atoms will be used. ",&
            usage="ABSOLUTE_POSITION",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! Must be present in each colvar and handled properly
       CALL create_point_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_colvar_xyz_d_section

! *****************************************************************************
!> \brief creates the colvar section regarded to the collective variables dist
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_colvar_xyz_od_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_colvar_xyz_od_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="XYZ_OUTERDIAG",&
            description="Section to define the cross term (XA-XA(0))*(XB-XB(0))+(XA-XA(0))*(YB-YB(0))"//&
            " or part of its components as a collective variable. The final term is given by the product "//&
            " of the components of A with the components of B.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., &
            error=error)
       NULLIFY(keyword, subsection)

       CALL keyword_create(keyword, name="ATOMS",&
            variants=(/"POINTS"/),&
            description="Specifies the index of the atoms/points A and B.",&
            usage="ATOMS {integer} {integer}",&
            n_var=2, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="COMPONENT_A",&
            description="Define the component of the position vector which will be used "//&
            "as a colvar for atom A.",&
            usage="AXIS (XYZ | X | Y | Z | XY| XZ | YZ)",&
            enum_c_vals=s2a( "XYZ","X", "Y", "Z", "XY", "XZ", "YZ"),&
            enum_i_vals=(/ do_clv_xyz, do_clv_x, do_clv_y,do_clv_z, do_clv_xy, do_clv_xz, do_clv_yz/),&
            default_i_val=do_clv_xyz, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="COMPONENT_B",&
            description="Define the component of the position vector which will be used "//&
            "as a colvar for atom B.",&
            usage="AXIS (XYZ | X | Y | Z | XY| XZ | YZ)",&
            enum_c_vals=s2a( "XYZ","X", "Y", "Z", "XY", "XZ", "YZ"),&
            enum_i_vals=(/ do_clv_xyz, do_clv_x, do_clv_y,do_clv_z, do_clv_xy, do_clv_xz, do_clv_yz/),&
            default_i_val=do_clv_xyz, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PBC",&
            description="Whether periodic boundary conditions should be applied on the "//&
            "atomic position before computing the colvar or not.",&
            usage="PBC",&
            default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! Must be present in each colvar and handled properly
       CALL create_point_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_colvar_xyz_od_section

! *****************************************************************************
!> \brief energy as collective variable
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Sebastiano Caravati
! *****************************************************************************
  SUBROUTINE create_colvar_u_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_colvar_u_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="u",&
            description="Section to define the energy as a generalized collective variable.",&
            n_keywords=0, n_subsections=0, repeats=.FALSE., &
            error=error)

       NULLIFY(subsection, keyword)
       CALL section_create(subsection,name="MIXED",&
            description="This section allows to use any function of the energy subsystems "//&
                        " in a mixed_env calculation as a collective variable.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., &
            error=error)

       CALL keyword_create(keyword, name="ENERGY_FUNCTION",&
            description="Specifies the functional form of the collective variable in mathematical notation.",&
            usage="ENERGY_FUNCTION (E1+E2-LOG(E1/E2))",  type_of_var=lchar_t,&
            n_var=1, error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="VARIABLES",&
            description="Defines the variables of the functional form. To allow an efficient"//&
            " mapping the order of the energy variables will be considered identical to the"//&
            " order of the force_eval in the force_eval_order list.",&
            usage="VARIABLES x",  type_of_var=char_t,&
            n_var=-1, error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PARAMETERS",&
            description="Defines the parameters of the functional form",&
            usage="PARAMETERS a b D",  type_of_var=char_t,&
            n_var=-1, repeats=.TRUE., error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="VALUES",&
            description="Defines the values of parameter of the functional form",&
            usage="VALUES ",  type_of_var=real_t,&
            n_var=-1, repeats=.TRUE., unit_str="internal_cp2k", error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="UNITS",&
            description="Optionally, allows to define valid CP2K unit strings for each parameter value. "//&
                        "It is assumed that the corresponding parameter value is specified in this unit.",&
            usage="UNITS angstrom eV*angstrom^-1 angstrom^1 K",  type_of_var=char_t,&
            n_var=-1, repeats=.TRUE., error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DX",&
            description="Parameter used for computing the derivative with the Ridders' method.",&
            usage="DX <REAL>", default_r_val=0.1_dp, unit_str="bohr", error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ERROR_LIMIT",&
            description="Checks that the error in computing the derivative is not larger than "//&
            "the value set; in case error is larger a warning message is printed.",&
            usage="ERROR_LIMIT <REAL>", default_r_val=1.0E-12_dp, error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)
    END IF

  END SUBROUTINE create_colvar_u_section

! *****************************************************************************
!> \brief creates the colvar section regarded to the collective variables distance
!>      of a point from a plane
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_colvar_d_pl_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_colvar_d_pl_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="distance_point_plane",&
            description="Section to define the distance of a point from a plane "//&
            "as a collective variables.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., &
            error=error)
       NULLIFY(keyword, subsection)

       CALL keyword_create(keyword, name="PBC",&
            description="Whether periodic boundary conditions should be applied on the "//&
            "atomic position before computing the colvar or not.",&
            usage="PBC",&
            default_l_val=.TRUE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ATOMS_PLANE",&
            variants=(/"POINTS_PLANE"/),&
            description="Specifies the indexes of atoms/points defining the plane.",&
            usage="ATOMS_PLANE <INTEGER> <INTEGER> <INTEGER>",&
            n_var=3, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ATOM_POINT",&
            variants=(/"POINT_POINT"/),&
            description="Specifies the atom/point index defining the point.",&
            usage="ATOM_POINT <INTEGER>",&
            n_var=1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! Must be present in each colvar and handled properly
       CALL create_point_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_colvar_d_pl_section

! *****************************************************************************
!> \brief creates the colvar section regarded to the collective variables
!>      angles betweem two planes
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_colvar_a_pl_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_colvar_a_pl_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="angle_plane_plane",&
            description="This section defines the angle between two planes "//&
            "as a collective variables.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., &
            error=error)
       NULLIFY(keyword, subsection)

       CALL section_create(subsection,name="PLANE",&
            description="This section defines the plane. When using this colvar, "//&
            "two plane section must be defined!",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., &
            error=error)

       CALL keyword_create(keyword, name="DEF_TYPE",&
            description="Specify how the plane is defined: either by 3 atoms or by a fixed normal "//&
            "vector. At least one plane must be defined through atoms.",&
            usage="DEF_TYPE ATOMS",&
            default_i_val=plane_def_atoms,&
            enum_c_vals=s2a("ATOMS", "VECTOR"),&
            enum_desc=s2a("Plane defined by the position of 3 atoms",&
            "Plane defined by a fixed normal vector"),&
            enum_i_vals=(/plane_def_atoms,plane_def_vec/),error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)


       CALL keyword_create(keyword, name="ATOMS",&
            description="Specifies the indexes of 3 atoms/points defining the plane.",&
            usage="ATOMS <INTEGER> <INTEGER> <INTEGER>",&
            n_var=3, type_of_var=integer_t, error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NORMAL_VECTOR",&
            description="Alternatively to 3 atoms/points one can define one of the two, "//&
            "planes by defining its NORMAL vector.",&
            usage="NORMAL_VECTOR 0.0 1.0 0.0",&
            n_var=3, type_of_var=real_t, error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       ! Must be present in each colvar and handled properly
       CALL create_point_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)
    END IF
  END SUBROUTINE create_colvar_a_pl_section

! *****************************************************************************
!> \brief create a geometrical point as a function of several atom coordinates
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_point_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_point_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="POINT",&
            description="Enables the possibility to use geometrical centers instead of single atoms"//&
            " to define colvars",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., &
            error=error)

       NULLIFY(keyword)

      CALL keyword_create(keyword, name="TYPE",&
            description="Chooses the type of geometrical point",&
            usage="type (GEO_CENTER|FIX_POINT)",&
            enum_c_vals=s2a( "GEO_CENTER","FIX_POINT"),&
            enum_desc=s2a("Conmputes the geometrical center of the listed atoms",&
                          "Defines a fixed point in space"),&
            enum_i_vals=(/ do_clv_geo_center, do_clv_fix_point /),&
            default_i_val=do_clv_geo_center, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ATOMS",&
            description="Specifies the indexes of atoms defining the geometrical center",&
            usage="ATOMS {integer} {integer} {integer} {integer}",&
            n_var=-1, type_of_var=integer_t, repeats=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="WEIGHTS",&
            description="Specifies the weights for a weighted geometrical center. Default is 1/natoms for every atom",&
            usage="WEIGHTS {real} {real} {real} {real}",&
            n_var=-1, type_of_var=real_t, repeats=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="XYZ",&
            description="Specifies the xyz of the fixed point (if the case)",&
            usage="XYZ {real} {real} {real}",&
            n_var=3, type_of_var=real_t, unit_str="bohr",&
            repeats=.FALSE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_point_section

! *****************************************************************************
!> \brief collective variables specifying torsion
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_colvar_qparm_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_colvar_qparm_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="qparm",&
            description="Section to define the Q parameter (crystalline order parameter) as a collective variable.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., &
            error=error)

       NULLIFY(keyword, subsection)

       CALL keyword_create(keyword, name="ATOMS_FROM",&
            variants=(/"POINTS_FROM"/),&
            description="Specify indexes of atoms/points building the coordination variable. ",&
            usage="ATOMS_FROM {integer} {integer} ..",repeats=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ATOMS_TO",&
            variants=(/"POINTS_TO"/),&
            description="Specify indexes of atoms/points building the coordination variable. ",&
            usage="ATOMS_TO {integer} {integer} ..",repeats=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RCUT",&
            description="Specifies the distance cutoff for neighbors.",&
            usage="RCUT {real}",&
            n_var=1, unit_str="angstrom",type_of_var=real_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="L",&
            description="Specifies the L spherical harmonics from Ylm.",&
            usage="L {integer}",&
            n_var=1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ALPHA",&
            description="Specifies the width of the Fermi-Dirac style smearing around RCUT.",&
            usage="ALPHA {real}",unit_str="angstrom^-1",default_r_val=0.0_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! Must be present in each colvar and handled properly
       CALL create_point_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_colvar_qparm_section

! *****************************************************************************
!> \brief collective variables specifying hydronium solvation
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Marcel Baer
! *****************************************************************************
  SUBROUTINE create_colvar_hydronium_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'create_colvar_hydronium_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="hydronium",&
            description="Section to define the formation of a hydronium as a collective variable.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., &
            error=error)

       NULLIFY(keyword, subsection)

       CALL keyword_create(keyword, name="OXYGENS",&
            description="Specify indexes of atoms building the coordination variable. ",&
            usage="OXYGENS {integer} {integer} ..", repeats=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="HYDROGENS",&
            description="Specify indexes of atoms building the coordination variable. ",&
            usage="HYDROGENS {integer} {integer} ..", repeats=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ROO",&
            variants=(/"R_OO"/),&
            description="Specify the ROO parameter in the coordination function.",&
            usage="ROO {real}",default_r_val=3.0_dp,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="pNO",&
            variants=(/"EXPON_NUMERATORA"/),&
            description="Sets the value of the numerator of the exponential factor"//&
            "in the coordination FUNCTION.",&
            usage="pNO {integer}",default_i_val=6,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="qNO",&
            variants=(/"EXPON_DENOMINATORA"/),&
            description="Sets the value of the denominator of the exponential factor"//&
            "in the coordination FUNCTION.",&
            usage="qNO {integer}",default_i_val=12,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ROH",&
            variants=(/"R_OH"/),&
            description="Specify the ROH parameter in the coordination function.",&
            usage="ROH {real}",default_r_val=3.0_dp,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="pNH",&
            variants=(/"EXPON_NUMERATORB"/),&
            description="Sets the value of the numerator of the exponential factor"//&
            "in the coordination FUNCTION.",&
            usage="pNH {integer}",default_i_val=6,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="qNH",&
            variants=(/"EXPON_DENOMINATORB"/),&
            description="Sets the value of the denominator of the exponential factor"//&
            "in the coordination FUNCTION.",&
            usage="qNH {integer}",default_i_val=12,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NH",&
            variants=(/"NHtest"/),&
            description="Specify the NH parameter in the hydronium function.",&
            usage="NH {real}",default_r_val=3.0_dp,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="p",&
            variants=(/"EXPON_NUMERATOR"/),&
            description="Sets the value of the numerator of the exponential factor"//&
            "in the coordination FUNCTION.",&
            usage="p {integer}",default_i_val=8,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="q",&
            variants=(/"EXPON_DENOMINATOR"/),&
            description="Sets the value of the denominator of the exponential factor"//&
            "in the coordination FUNCTION.",&
            usage="q {integer}",default_i_val=16,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="LAMBDA",&
            variants=(/"LAMBDA"/),&
            description="Specify the LAMBDA parameter in the hydronium function.",&
            usage="LAMBDA {real}",default_r_val=10.0_dp,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! Must be present in each colvar and handled properly
       CALL create_point_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE create_colvar_hydronium_section

! *****************************************************************************
!> \brief ...
!> \param section ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE create_colvar_rmsd_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_colvar_rmsd_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection, subsubsection

    failure=.FALSE.

    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="rmsd",&
            description="Section to define a CV as function of RMSD computed with respect to"//&
            " given reference configurations. For 2 configurations the colvar is equal to:"//&
            " ss = (RMSDA-RMSDB)/(RMSDA+RMSDB), while if only 1 configuration is given, then the"//&
            " colvar is just the RMSD from that frame.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., &
            error=error)

       NULLIFY(keyword, subsection, subsubsection)
       CALL keyword_create(keyword, name="SUBSET_TYPE",&
            description="Define the subsytem used to compute the RMSD",&
            usage="SUBSET_TYPE ALL",&
            enum_c_vals=s2a( "ALL","LIST", "WEIGHT_LIST"),&
            enum_i_vals=(/ rmsd_all, rmsd_list, rmsd_weightlist/),&
            default_i_val=rmsd_all, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ALIGN_FRAMES",&
            description="Whether the reference frames should be aligned to minimize the RMSD",&
            usage="ALIGN_FRAME",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ATOMS",&
            description="Specify indexes of atoms building the subset. ",&
            usage="ATOMS {integer} {integer} ..",repeats=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="WEIGHTS",&
            description="Specify weights of atoms building the subset. ",&
            usage="weightS {real} {real} ..",repeats=.TRUE.,&
            n_var=-1, type_of_var=real_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_create(subsection,name="FRAME",&
            description="Specify coordinates of the frame (number of frames can be either 1 or 2)",&
            repeats=.TRUE.,  error=error)

       CALL keyword_create(keyword, name="COORD_FILE_NAME",&
            description="Name of the xyz file with coordinates (alternative to &COORD section)",&
            usage="COORD_FILE_NAME <CHAR>",&
            default_lc_val="",error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_coord_section_cv(subsubsection,"RMSD",error=error)
       CALL section_add_subsection(subsection,subsubsection,error=error)
       CALL section_release(subsubsection,error=error)

       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF

  END  SUBROUTINE create_colvar_rmsd_section

! *****************************************************************************
!> \brief collective variables specifying the space orthogonal to the reaction path
!>        in the space spanned by the involved collective coordinates
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author fschiff
! *****************************************************************************
  SUBROUTINE create_colvar_rpath_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_colvar_rpath_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="REACTION_PATH",&
            description="Section defining a one dimensional reaction path in an Q-dimensional space of colvars. "//&
            "Constraining this colvar, allows to sample the space orthogonal to the reaction path, "//&
            "both in the Q-dimensional colvar and 3N-Q remaining coordinates. "//&
            "For the details of the function see cited literature.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., &
            citations=(/Branduardi2007/), &
            error=error)

       CALL keywords_colvar_path(section,error)
    END IF
  END SUBROUTINE create_colvar_rpath_section

! *****************************************************************************
!> \brief Distance from reaction path
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author 01.2010
! *****************************************************************************
  SUBROUTINE create_colvar_dpath_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_colvar_dpath_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.
    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="DISTANCE_FROM_PATH",&
            description="Section defining the distance from a one dimensional reaction "//&
            "path in an Q-dimensional space of colvars. "//&
            "Constraining this colvar, allows to sample the space equidistant to the reaction path, "//&
            "both in the Q-dimensional colvar and 3N-Q remaining coordinates. "//&
            "For the details of the function see cited literature.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., &
            citations=(/Branduardi2007/), &
            error=error)

       CALL keywords_colvar_path(section,error)
    END IF
  END SUBROUTINE create_colvar_dpath_section

! *****************************************************************************
!> \brief Section describinf keywords for both reaction path and distance from reaction path
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author 01.2010
! *****************************************************************************
  SUBROUTINE keywords_colvar_path(section,error)

    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'keywords_colvar_path', &
      routineP = moduleN//':'//routineN

    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, subsection, &
                                                subsubsection

       NULLIFY(keyword, subsection, subsubsection, print_key)
       CALL create_colvar_section(subsection,skip_recursive_colvar=.TRUE.,error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL keyword_create(keyword, name="DISTANCES_RMSD",&
            description=" ",&
            usage="DISTANCES_RMSD T", &
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RMSD",&
            description=" ",&
            usage="RMSD T", &
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)


       CALL keyword_create(keyword, name="SUBSET_TYPE",&
            description="Define the subsytem used to compute the RMSD",&
            usage="SUBSET_TYPE ALL",&
            enum_c_vals=s2a( "ALL","LIST"),&
            enum_i_vals=(/ rmsd_all, rmsd_list/),&
            default_i_val=rmsd_all, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ALIGN_FRAMES",&
            description="Whether the reference frames should be aligned to minimize the RMSD",&
            usage="ALIGN_FRAME",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ATOMS",&
            description="Specify indexes of atoms building the subset. ",&
            usage="ATOMS {integer} {integer} ..",repeats=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_create(subsection,name="FRAME",&
            description="Specify coordinates of the frame",&
            repeats=.TRUE.,  error=error)

       CALL keyword_create(keyword, name="COORD_FILE_NAME",&
            description="Name of the xyz file with coordinates (alternative to &COORD section)",&
            usage="COORD_FILE_NAME <CHAR>",&
            default_lc_val="",error=error)
       CALL section_add_keyword(subsection,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_coord_section_cv(subsubsection,"RMSD",error=error)
       CALL section_add_subsection(subsection,subsubsection,error=error)
       CALL section_release(subsubsection,error=error)

       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL keyword_create(keyword, name="FUNCTION",&
            description="Specifies the ith element of the vector valued function that defines the reaction path. "//&
            "This keyword needs to repeat exactly Q times, and the order must match the order of the colvars. "//&
            "The VARIABLE (e.g. T) which parametrises the curve can be used as the target of a constraint.",&
            usage="FUNCTION (sin(T+2)+2*T)", type_of_var=lchar_t,&
            n_var=1, default_lc_val="0",repeats=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="VARIABLE",&
            description="Specifies the name of the variable that parametrises the FUNCTION "//&
            "defining the reaction path.",&
            usage="VARIABLE T", type_of_var=char_t,&
            n_var=1, repeats=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="LAMBDA",&
            description="Specifies the exponent of the Gaussian used in the integral representation of the colvar."//&
            "The shape of the space orthogonal to the reaction path is defined by this choice. "//&
            "In the limit of large values, it is given by the plane orthogonal to the path."//&
            "In practice, modest values are required for stable numerical integration.",&
            usage="LAMBDA {real}",&
            type_of_var=real_t,default_r_val=5.0_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="STEP_SIZE",&
            description="Step size in the numerical integration, "//&
            "a few thousand points are common, and the proper number also depends on LAMBDA.",&
            usage="STEP_SIZE {real}",&
            type_of_var=real_t,default_r_val=0.01_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="RANGE",&
            description="The range of VARIABLE used for the parametrisation.",&
            usage="RANGE <REAL> <REAL>",&
            n_var=2,type_of_var=real_t,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL cp_print_key_section_create(print_key,name="MAP",&
            description="Activating this print key will print once a file with the values of the FUNCTION on a grid "//&
            "of COLVAR values in a specified range. "//&
            "GRID_SPACING and RANGE for every COLVAR has to be specified again in the same order as they are in the input.", &
            print_level=high_print_level,filename="PATH",&
            error=error)

       CALL keyword_create(keyword,name="RANGE",&
            description="The range of of the grid of the COLVAR.",&
            usage="RANGE <REAL> <REAL>",&
            n_var=2,type_of_var=real_t,repeats=.TRUE.,&
            error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="GRID_SPACING",&
            description="Distance between two gridpoints for the grid on the COLVAR",&
            usage="STEP_SIZE {real}",repeats=.TRUE.,&
            type_of_var=real_t,default_r_val=0.01_dp, error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

  END SUBROUTINE keywords_colvar_path

! *****************************************************************************
!> \brief Colvar allowing a combination of COLVARS
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Teodoro Laino [tlaino] - 12.2008
! *****************************************************************************
  SUBROUTINE create_colvar_comb_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_colvar_comb_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="COMBINE_COLVAR",&
            description="Allows the possibility to combine several COLVARs into one COLVAR "//&
            "with a generic function.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., &
            error=error)

       NULLIFY(keyword, subsection)
       CALL create_colvar_section(subsection,skip_recursive_colvar=.TRUE.,error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL keyword_create(keyword, name="FUNCTION",&
            description="Specifies the function used to combine different COLVARs into one.",&
! *****************************************************************************
!> \brief ...
!> \param CV1^2 ...
!> \param CV2^2 ...
!> \param " ...
!> \param type_of_var=lchar_t ...
!> \param n_var=1 ...
!> \param error=error ...
! *****************************************************************************
            usage="FUNCTION SQRT(CV1^2+CV2^2)", type_of_var=lchar_t,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="VARIABLES",&
            description="Specifies the name of the variable that parametrises the FUNCTION "//&
            "defining how COLVARS should be combined. The matching follows the same order of the "//&
            "COLVARS definition in the input file.",&
            usage="VARIABLE CV1 CV2 CV3", type_of_var=char_t,n_var=-1, repeats=.FALSE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="PARAMETERS",&
            description="Defines the parameters of the functional form",&
            usage="PARAMETERS a b D",  type_of_var=char_t,&
            n_var=-1, repeats=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="VALUES",&
            description="Defines the values of  parameter of the functional form",&
            usage="VALUES ",  type_of_var=real_t,&
            n_var=-1, repeats=.TRUE., unit_str="internal_cp2k", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="DX",&
            description="Parameter used for computing the derivative of the combination "//&
            "of COLVARs with the Ridders' method.",&
            usage="DX <REAL>", default_r_val=0.1_dp, unit_str="bohr", error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ERROR_LIMIT",&
            description="Checks that the error in computing the derivative is not larger than "//&
            "the value set; in case error is larger a warning message is printed.",&
            usage="ERROR_LIMIT <REAL>", default_r_val=1.0E-12_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_colvar_comb_section

! *****************************************************************************
!> \brief Creates the coord section
!> \param section the section to create
!> \param name ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author teo
! *****************************************************************************
  SUBROUTINE create_coord_section_cv(section,name,error)
    TYPE(section_type), POINTER              :: section
    CHARACTER(LEN=*), INTENT(IN)             :: name
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_coord_section_cv', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="coord",&
            description="The positions for "//TRIM(name)//" used for restart",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., &
            error=error)
       NULLIFY(keyword)

       CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_",&
            description="Specify positions of the system",repeats=.TRUE.,&
            usage="{Real} ...", type_of_var=real_t, n_var=-1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF
  END SUBROUTINE create_coord_section_cv

! *****************************************************************************
!> \brief collective variables specifying h bonds
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author alin m elena
! *****************************************************************************
  SUBROUTINE create_colvar_wc_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_colvar_wc_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="wc",&
            description="Section to define the hbond wannier centre as a collective variables.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., &
            error=error)
       NULLIFY(keyword, subsection)

       CALL keyword_create(keyword, name="RCUT",&
            description="Parameter used for computing the cutoff radius for searching "//&
            "the wannier centres around an atom",&
            usage="RCUT <REAL>", default_r_val=0.529177208590000_dp, unit_str="angstrom",&
            type_of_var=real_t,repeats=.FALSE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ATOMS",&
            variants=(/"POINTS"/),&
            description="Specifies the indexes of atoms/points defining the bond (Od, H, Oa).",&
            usage="ATOMS {integer} {integer} {integer}",&
            n_var=3, type_of_var=integer_t,repeats=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! Must be present in each colvar and handled properly
       CALL create_point_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_colvar_wc_section


   ! *****************************************************************************
!> \brief collective variables specifying h bonds= wire
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author alin m elena
! *****************************************************************************
  SUBROUTINE create_colvar_hbp_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_colvar_hbp_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="hbp",&
            description="Section to define the hbond wannier centre as a collective variables.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., &
            error=error)
       NULLIFY(keyword, subsection)

       CALL keyword_create(keyword, name="RCUT",&
            description="Parameter used for computing the cutoff radius for searching "//&
            "the wannier centres around an atom",&
            usage="RCUT <REAL>", default_r_val=0.529177208590000_dp, unit_str="angstrom",&
            type_of_var=real_t,repeats=.FALSE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SHIFT",&
            description="Parameter used for shifting each term in the sum ", &
            usage="SHIFT <REAL>", default_r_val=0.5_dp,&
            type_of_var=real_t,repeats=.FALSE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="NPOINTS",&
            description="The number of points in the path",&
            usage="NPOINTS {integer}",default_i_val=-1,&
            n_var=1,type_of_var=integer_t,repeats=.FALSE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ATOMS",&
            variants=(/"POINTS"/),&
            description="Specifies the indexes of atoms/points defining the bond (Od, H, Oa).",&
            usage="ATOMS {integer} {integer} {integer}",&
            n_var=3, type_of_var=integer_t,repeats=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! Must be present in each colvar and handled properly
       CALL create_point_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_colvar_hbp_section

! *****************************************************************************
!> \brief collective variables specifying ring puckering     
!> \brief D. Cremer and J.A. Pople, JACS 97 1354 (1975)
!> \param section the section to be created
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Marcel Baer
! *****************************************************************************
  SUBROUTINE create_colvar_ring_puckering_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'create_colvar_ring_puckering_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="RING_PUCKERING",&
            description="Section to define general ring puckering collective variables.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., &
            error=error)

       NULLIFY(keyword, subsection)

       CALL keyword_create(keyword, name="ATOMS",&
            variants=(/"POINTS"/),&
            description="Specifies the indexes of atoms/points defining the ring."//&
            "At least 4 Atoms are needed.",&
            usage="ATOMS {integer} {integer} {integer} ..",&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="COORDINATE",&
            description="Indicate the coordinate to be used. Follow the Cremer-Pople definition for a N ring."//&
            "0 is the total puckering variable Q,"//&
            "2..[N/2] are puckering coordinates."//&
            "-2..-[N/2-1] are puckering angles.",&
            usage="COORDINATE {integer}",default_i_val=0,&
            n_var=1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! Must be present in each colvar and handled properly
       CALL create_point_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE create_colvar_ring_puckering_section

! *****************************************************************************


END MODULE input_cp2k_colvar
