!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2007  CP2K developers group                          !
!-----------------------------------------------------------------------------!

!!****h* cp2k/input_cp2k_colvar *
!!
!!   NAME
!!     input_cp2k_colvar
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     teo & fawzi
!!
!!   MODIFICATION HISTORY
!!     - taken out of input_cp2k_motion
!!
!!   SOURCE
!****************************************************************************
MODULE input_cp2k_colvar
  USE cp_output_handling,              ONLY: cp_print_key_section_create
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
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_colvar'

PUBLIC :: create_colvar_section

CONTAINS

!!****f* input_cp2k/create_colvar_section *
!!
!!   NAME
!!     create_colvar_section
!!
!!   FUNCTION
!!     creates the colvar section
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - section: the section to be created
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo
!!
!!*** **********************************************************************
  SUBROUTINE create_colvar_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_colvar_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(section_type), POINTER              :: print_key, subsection

    failure=.FALSE.
    
    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="COLVAR",&
            description="This section specifies the nature of the collective variables.",&
            n_keywords=1, n_subsections=1, repeats=.TRUE., required=.TRUE.,&
            error=error)
       NULLIFY(subsection, print_key)

       CALL create_colvar_var_section(subsection, section, error)

       CALL section_create(subsection,name="print",&
            description="Controls the printing of the colvar specifications",&
            n_keywords=0, n_subsections=1, repeats=.TRUE., required=.FALSE.,&
            error=error)
       NULLIFY(print_key)
       CALL cp_print_key_section_create(print_key,"program_run_info",&
            description="Controls the printing of basic information during colvar setup.", &
            print_level=low_print_level,each=(/1/),add_last=(/1/),filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(subsection,print_key,error=error)
       CALL section_release(print_key,error=error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_colvar_section

!!****f* input_cp2k/create_colvar_var_section *
!!
!!   NAME
!!     create_colvar_var_section
!!
!!   FUNCTION
!!     creates the collective variables for the colvar section
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - section: the section to be created
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo
!!
!!*** **********************************************************************
  SUBROUTINE create_colvar_var_section(subsection, section, error)
    TYPE(section_type), POINTER              :: subsection, section
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

       CALL create_colvar_ddiff_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_colvar_qparm_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_colvar_hydronium_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE create_colvar_var_section

!!****f* input_cp2k/create_colvar_coord_section *
!!
!!   NAME
!!     create_colvar_coord_section
!!
!!   FUNCTION
!!     collective variables specifying coordination
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - section: the section to be created
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo
!!
!!*** **********************************************************************
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
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
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

    END IF

  END SUBROUTINE create_colvar_coord_section

!!****f* input_cp2k/create_colvar_ddiff_section *
!!
!!   NAME
!!     create_colvar_ddiff_section
!!
!!   FUNCTION
!!     collective variables specifying torsion
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - section: the section to be created
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo
!!
!!*** **********************************************************************
  SUBROUTINE create_colvar_ddiff_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_colvar_ddiff_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="distance_difference",&
            description="Section to define difference between two distances as collective variables.",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)
        
       NULLIFY(keyword, subsection)

       CALL keyword_create(keyword, name="ATOMS",&
            variants=(/"POINTS"/),&
            description="Specifies the indexes of atoms/points for the two bonds d1=(1-2) d2=(3-4)."//&
            " The difference is defined as d1-d2",&
            usage="ATOMS {integer} {integer} {integer} {integer}",required=.TRUE.,&
            n_var=4, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! Must be present in each colvar and handled properly
       CALL create_point_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_colvar_ddiff_section

!!****f* input_cp2k/create_colvar_torsion_section *
!!
!!   NAME
!!     create_colvar_torsion_section
!!
!!   FUNCTION
!!     collective variables specifying torsion
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - section: the section to be created
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo
!!
!!*** **********************************************************************
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
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
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

!!****f* input_cp2k/create_colvar_rot_section *
!!
!!   NAME
!!     create_colvar_rot_section
!!
!!   FUNCTION
!!     collective variables specifying torsion
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - section: the section to be created
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo
!!
!!*** **********************************************************************
  SUBROUTINE create_colvar_rot_section(section,error,description)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: description

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
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
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

!!****f* input_cp2k/create_colvar_angle_section *
!!
!!   NAME
!!     create_colvar_angle_section
!!
!!   FUNCTION
!!     collective variables specifying angles
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - section: the section to be created
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo
!!
!!*** **********************************************************************
  SUBROUTINE create_colvar_angle_section(section,error,description)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: description

    CHARACTER(len=*), PARAMETER :: routineN = 'create_colvar_angle_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.
    
    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       IF (PRESENT(description)) THEN
          CALL section_create(section,name="angle",&
               description=description,&
               n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
               error=error)
       ELSE
          CALL section_create(section,name="angle",&
               description="Section to define the angle as a collective variables.",&
               n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
               error=error)
       END IF

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
!****************************************************************************


!!****f* input_cp2k/create_colvar_dist_section *
!!
!!   NAME
!!     create_colvar_dist_section
!!
!!   FUNCTION
!!     creates the colvar section regarded to the collective variables dist
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - section: the section to be created
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo
!!
!!*** **********************************************************************
  SUBROUTINE create_colvar_dist_section(section,error,description)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: description

    CHARACTER(len=*), PARAMETER :: routineN = 'create_colvar_dist_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.
    
    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       IF (PRESENT(description)) THEN
          CALL section_create(section,name="distance",&
               description=description,&
               n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
               error=error)
       ELSE
          CALL section_create(section,name="distance",&
               description="Section to define the distance as a collective variables.",&
               n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
               error=error)
       END IF

       NULLIFY(keyword, subsection)

       CALL keyword_create(keyword, name="ATOMS",&
            variants=(/"POINTS"/),&
            description="Specifies the indexes of atoms/points defining the distance.",&
            usage="ATOMS {integer} {integer}",required=.TRUE.,&
            n_var=2, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       ! Must be present in each colvar and handled properly
       CALL create_point_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

    END IF

  END SUBROUTINE create_colvar_dist_section

!!****f* input_cp2k/create_colvar_d_pl_section *
!!
!!   NAME
!!     create_colvar_d_pl_section
!!
!!   FUNCTION
!!     creates the colvar section regarded to the collective variables distance 
!!     of a point from a plane
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - section: the section to be created
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo
!!
!!*** **********************************************************************
  SUBROUTINE create_colvar_d_pl_section(section,error,description)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: description

    CHARACTER(len=*), PARAMETER :: routineN = 'create_colvar_d_pl_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.
    
    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       IF (PRESENT(description)) THEN
          CALL section_create(section,name="distance_point_plane",&
               description=description,&
               n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
               error=error)
       ELSE
          CALL section_create(section,name="distance_point_plane",&
               description="Section to define the distance of a point from a plane"//&
               "as a collective variables.",&
               n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
               error=error)
       END IF

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


!!****f* input_cp2k/create_point_section *
!!
!!   NAME
!!     create_point_section
!!
!!   FUNCTION
!!     create a geometrical point as a function of several atom coordinates
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - section: the section to be created
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo
!!
!!*** **********************************************************************
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
            n_var=-1, type_of_var=real_t, repeats=.FALSE., error=error)
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

!!****f* input_cp2k/create_colvar_qparm_section *
!!
!!   NAME
!!     create_colvar_qparm_section
!!
!!   FUNCTION
!!     collective variables specifying torsion
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - section: the section to be created
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     teo
!!
!!*** **********************************************************************
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
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
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

!!****f* input_cp2k/create_colvar_hudronium_section *
!!
!!   NAME
!!     create_colvar_hyrdonium_section
!!
!!   FUNCTION
!!     collective variables specifying hydronium solvation 
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - section: the section to be created
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Marcel Baer
!!
!!*** **********************************************************************
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
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
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

END MODULE input_cp2k_colvar
