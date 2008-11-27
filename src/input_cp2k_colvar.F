!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2008  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      - taken out of input_cp2k_motion
!> \author teo & fawzi
! *****************************************************************************
MODULE input_cp2k_colvar
  USE bibliography,                    ONLY: Branduardi2007
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
  USE input_val_types,                 ONLY: char_t,&
                                             integer_t,&
                                             lchar_t,&
                                             real_t
  USE kinds,                           ONLY: dp
  USE string_utilities,                ONLY: s2a
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_colvar'

PUBLIC :: create_colvar_section

CONTAINS

! *****************************************************************************
!> \brief creates the colvar section
!> \param section the section to be created
!> \param error variable to control error logging, stopping,... 
!>        see module cp_error_handling 
!> \author teo
! *****************************************************************************
  RECURSIVE SUBROUTINE create_colvar_section(section,skip_reaction_path,error)
    TYPE(section_type), POINTER              :: section
    LOGICAL, OPTIONAL                        :: skip_reaction_path
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_colvar_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure, skip
    TYPE(section_type), POINTER              :: print_key, subsection

    failure=.FALSE.
    skip=.FALSE.
    IF(PRESENT(skip_reaction_path))skip=skip_reaction_path
    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="COLVAR",&
            description="This section specifies the nature of the collective variables.",&
            n_keywords=1, n_subsections=1, repeats=.TRUE., required=.TRUE.,&
            error=error)
       NULLIFY(subsection, print_key)

       CALL create_colvar_var_section(subsection=subsection,&
            section=section,skip_reaction_path=skip,error=error)

       CALL section_create(subsection,name="PRINT",&
            description="Controls the printing of the colvar specifications",&
            n_keywords=0, n_subsections=1, repeats=.TRUE., required=.FALSE.,&
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

    END IF

  END SUBROUTINE create_colvar_section

! *****************************************************************************
!> \brief creates the collective variables for the colvar section
!> \param section the section to be created
!> \param error variable to control error logging, stopping,... 
!>        see module cp_error_handling 
!> \author teo
! *****************************************************************************
  RECURSIVE SUBROUTINE create_colvar_var_section(subsection, section, skip_reaction_path, error)
    TYPE(section_type), POINTER              :: subsection, section
    LOGICAL, INTENT(IN)                      :: skip_reaction_path
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

       CALL create_colvar_d_pl_section(subsection,error)
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

       IF(.NOT. skip_reaction_path)THEN
          CALL create_colvar_reaction_path_section(subsection,error)
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
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)
        
       NULLIFY(subsection, keyword)

       CALL keyword_create(keyword, name="ATOMS_FROM",&
            variants=(/"POINTS_FROM"/),&
            description="Specify indexes of atoms/points building the coordination variable. ",&
            usage="ATOMS_FROM {integer} {integer} ..",required=.TRUE.,repeats=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ATOMS_TO",&
            variants=(/"POINTS_TO"/),&
            description="Specify indexes of atoms/points building the coordination variable. ",&
            usage="ATOMS_TO {integer} {integer} ..",required=.TRUE.,repeats=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ATOMS_TO_B",&
            variants=(/"POINTS_TO_B"/),&
            description="For the CV given by the multiplication of two coorination numbers,"//&
            " here specify indexes of the third set of atoms/points. ",&
            usage="ATOMS_TO_B {integer} {integer} ..",required=.TRUE.,repeats=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="KINDS_FROM",&
            description="Specify alternatively kinds of atoms building the coordination variable.",&
            usage="KINDS_FROM {CHAR} {CHAR} ..",required=.FALSE.,repeats=.TRUE.,&
            n_var=-1, type_of_var=char_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="KINDS_TO",&
            description="Specify alternatively kinds of atoms building the coordination variable.",&
            usage="KINDS_TO {CHAR} {CHAR} ..",required=.FALSE.,repeats=.TRUE.,&
            n_var=-1, type_of_var=char_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="KINDS_TO_B",&
            description="For the CV given by the multiplication of two coorination numbers,"//&
           " here specify alternatively kinds of atoms building the coordination variable.",&
            usage="KINDS_TO_B {CHAR} {CHAR} ..",required=.FALSE.,repeats=.TRUE.,&
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
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)
        
       NULLIFY(keyword, subsection)

       CALL keyword_create(keyword, name="ATOMS",&
            variants=(/"POINTS"/),&
            description="Specifies the indexes of atoms/points for the two bonds d1=(1-2) d2=(3-4).",&
            usage="ATOMS {integer} {integer} {integer} {integer}",required=.TRUE.,&
            n_var=4, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="COEFFICIENT",&
            description="Specifies the coefficient in the function for the constraint."//&
            " -1.0 has to be used for distance difference, 1.0 for distance addition",&
            usage="COEFFICIENT {real}",required=.TRUE.,&
            type_of_var=real_t, error=error)
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
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)
        
       NULLIFY(keyword, subsection)

       CALL keyword_create(keyword, name="ATOMS",&
            variants=(/"POINTS"/),&
            description="Specifies the indexes of atoms/points defining the torsion.",&
            usage="ATOMS {integer} {integer} {integer} {integer}",required=.TRUE.,&
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
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)
        
       NULLIFY(keyword, subsection)

       CALL keyword_create(keyword, name="P1_BOND1",&
            description="Specifies the index of atom/point defining the first point"//&
            "of the first bond/line.",&
            usage="P1_BOND1 {integer}",required=.TRUE.,&
            n_var=1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="P2_BOND1",&
            description="Specifies the index of atom/point defining the second point"//&
            "of the first bond/line.",&
            usage="P2_BOND1 {integer}",required=.TRUE.,&
            n_var=1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="P1_BOND2",&
            description="Specifies the index of atom/point defining the first point"//&
            "of the second bond/line.",&
            usage="P1_BOND2 {integer}",required=.TRUE.,&
            n_var=1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="P2_BOND2",&
            description="Specifies the index of atom/point defining the second point"//&
            "of the second bond/line.",&
            usage="P2_BOND2 {integer}",required=.TRUE.,&
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
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)
       NULLIFY(keyword, subsection)

       CALL keyword_create(keyword, name="ATOMS",&
            variants=(/"POINTS"/),&
            description="Specifies the indexes of atoms/points defining the angle.",&
            usage="ATOMS {integer} {integer} {integer}",required=.TRUE.,&
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
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)
       NULLIFY(keyword, subsection)

       CALL keyword_create(keyword, name="ATOMS",&
            variants=(/"POINTS"/),&
            description="Specifies the indexes of atoms/points defining the distance.",&
            usage="ATOMS {integer} {integer}",required=.TRUE.,&
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
            description="Section to define the distance of a point from a plane"//&
            "as a collective variables.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)
       NULLIFY(keyword, subsection)

       CALL keyword_create(keyword, name="ATOMS_PLANE",&
            variants=(/"POINTS_PLANE"/),&
            description="Specifies the indexes of atoms/points defining the plane.",&
            usage="ATOMS_PLANE <INTEGER> <INTEGER> <INTEGER>",required=.TRUE.,&
            n_var=3, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ATOM_POINT",&
            variants=(/"POINT_POINT"/),&
            description="Specifies the atom/point index defining the point.",&
            usage="ATOM_POINT <INTEGER>",required=.TRUE.,&
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
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
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
            usage="ATOMS {integer} {integer} {integer} {integer}",required=.TRUE.,&
            n_var=-1, type_of_var=integer_t, repeats=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="WEIGHTS",&
            description="Specifies the weights for a weighted geometrical center. Default is 1/natoms for every atom",&
            usage="WEIGHTS {real} {real} {real} {real}",required=.FALSE.,&
            n_var=-1, type_of_var=real_t, repeats=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="XYZ",&
            description="Specifies the xyz of the fixed point (if the case)",&
            usage="XYZ {real} {real} {real}",required=.TRUE.,&
            n_var=3, type_of_var=real_t, repeats=.FALSE., error=error)
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
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)
        
       NULLIFY(keyword, subsection)

       CALL keyword_create(keyword, name="ATOMS_FROM",&
            variants=(/"POINTS_FROM"/),&
            description="Specify indexes of atoms/points building the coordination variable. ",&
            usage="ATOMS_FROM {integer} {integer} ..",required=.TRUE.,repeats=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ATOMS_TO",&
            variants=(/"POINTS_TO"/),&
            description="Specify indexes of atoms/points building the coordination variable. ",&
            usage="ATOMS_TO {integer} {integer} ..",required=.TRUE.,repeats=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RCUT",&
            description="Specifies the distance cutoff for neighbors.",&
            usage="RCUT {real}",required=.TRUE.,&
            n_var=1, unit_str="angstrom",type_of_var=real_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="L",&
            description="Specifies the L spherical harmonics from Ylm.",&
            usage="L {integer}",required=.TRUE.,&
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
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)
        
       NULLIFY(keyword, subsection)

       CALL keyword_create(keyword, name="OXYGENS",&
            description="Specify indexes of atoms building the coordination variable. ",&
            usage="OXYGENS {integer} {integer} ..",required=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="HYDROGENS",&
            description="Specify indexes of atoms building the coordination variable. ",&
            usage="HYDROGENS {integer} {integer} ..",required=.TRUE.,&
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
!> \brief collective variables specifying the space orthogonal to the reaction path
!>        in the space spanned by the involved collective coordinates
!> \param section the section to be created
!> \param error variable to control error logging, stopping,... 
!>        see module cp_error_handling 
!> \author fschiff
! *****************************************************************************
  SUBROUTINE create_colvar_reaction_path_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'create_colvar_reaction_path_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key, subsection

    failure=.FALSE.

    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="REACTION_PATH",&
            description="Section defining a one dimensional reaction path in an Q-dimensional space of colvars. "//&
            "Constraining this colvar, allows to sample the space orthogonal to the reaction path, "//&
            "both in the Q-dimensional colvar and 3N-Q remaining coordinates. "//&
            "For the details of the function see cited literature.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            citations=(/Branduardi2007/), &
            error=error)
        
       NULLIFY(keyword, subsection,print_key)
       CALL create_colvar_section(subsection,skip_reaction_path=.TRUE.,error=error)
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
            usage="LAMBDA {real}",required=.TRUE.,&
            type_of_var=real_t,default_r_val=5.0_dp, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="STEP_SIZE",&
            description="Step size in the numerical integration, "//&
                        "a few thousand points are common, and the proper number also depends on LAMBDA.",&
            usage="STEP_SIZE {real}",required=.TRUE.,&
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
            print_level=high_print_level,filename="FUNCTION_MAP",&
            error=error)

       CALL keyword_create(keyword,name="RANGE",&
            description="The range of of the grid of the COLVAR.",&
            usage="RANGE <REAL> <REAL>",&
            n_var=2,type_of_var=real_t,repeats=.TRUE.,required=.TRUE.,&
            error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="GRID_SPACING",&
            description="Distance between two gridpoints for the grid on the COLVAR",&
            usage="STEP_SIZE {real}",required=.TRUE.,repeats=.TRUE.,&
            type_of_var=real_t,default_r_val=0.01_dp, error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)
    END IF
  END SUBROUTINE create_colvar_reaction_path_section

END MODULE input_cp2k_colvar
