!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2006  CP2K developers group                                 !
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
  use kinds, only: dp
  use cp_para_types, only: cp_para_env_type
  use timings, only: timeset, timestop
  USE input_section_types, ONLY: section_type, section_create, &
       section_add_keyword, section_add_subsection, section_release,&
       section_vals_type,&
       section_vals_get_subs_vals, section_vals_get
  USE input_constants
  USE input_keyword_types, ONLY: keyword_type, keyword_create,&
       keyword_release
  USE input_val_types, ONLY: val_type, no_t, logical_t, integer_t, real_t, char_t,&
       enum_t
  USE cp_output_handling, only: cp_print_key_section_create
  USE cp_units, only: cp_unit_to_cp2k
  USE string_utilities, ONLY: s2a
#include "cp_common_uses.h"
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='input_cp2k_colvar'

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
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_colvar_section',&
         routineP=moduleN//':'//routineN
    TYPE(section_type), POINTER :: subsection
    TYPE(keyword_type), pointer :: keyword    

    failure=.FALSE.
    
    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="colvar",&
            description="This section specify the nature of the collective variables.",&
            n_keywords=1, n_subsections=1, repeats=.TRUE., required=.TRUE.,&
            error=error)
        
       NULLIFY(subsection,keyword)

       CALL keyword_create(keyword, name="LAMBDA",&
            description="Specifies the lambda parameter of the collective variable in the"//&
            " extended lagrangian scheme.",&
            usage="LAMBDA {real}",type_of_var=real_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="MASS",&
            description="Specifies the mass parameter of the collective variable in the"//&
            " extended lagrangian scheme.",&
            usage="MASS {real}",unit_str='amu',&
            default_r_val=0.0_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SCALE",&
            description="Specifies the scale factor for the following collective variable.",&
            usage="SCALE {real}",&
            default_r_val=1.0_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="WALL_PLUS",&
            variants=(/"WALLP"/),&
            description="Activates the reflective wall on the upper limit of the collective variable.",&
            usage="WALL_PLUS {real}",type_of_var=real_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="WALL_MINUS",&
            variants=(/"WALLM"/),&
            description="Activates the reflective wall on the lower limit of the collective variable.",&
            usage="WALL_MINUS {real}",type_of_var=real_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SUBSYS",&
            description="Specifies the subsystem on which to apply the collective variables. Default 1.",&
            usage="SUBSYS {integer}",&
            default_i_val=1,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)       
    
       CALL create_colvar_var_section(subsection, section, error)
       
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
    TYPE(section_type), POINTER :: section, subsection
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_colvar_var_section',&
         routineP=moduleN//':'//routineN

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

    END IF
  END SUBROUTINE create_colvar_var_section

!!****f* input_cp2k/create_colvar_coord_section *
!!
!!   NAME
!!     create_colvar_coord_section
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
  SUBROUTINE create_colvar_coord_section(section,error)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_colvar_coord_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), pointer :: keyword    

    failure=.FALSE.
    
    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       CALL section_create(section,name="coordination",&
            description="Section to define the coordination as a collective variables.",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)
        
       NULLIFY(keyword)

       CALL keyword_create(keyword, name="ATOMS_FROM",&
            description="Specify indexes of atoms building the coordination variable. ",&
            usage="ATOMS_FROM {integer} {integer} ..",required=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="ATOMS_TO",&
            description="Specify indexes of atoms building the coordination variable. ",&
            usage="ATOMS_TO {integer} {integer} ..",required=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="KINDS_FROM",&
            description="Specify kinds of atoms building the coordination variable. ",&
            usage="KINDS_FROM {integer} {integer} ..",required=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="KINDS_TO",&
            description="Specify kinds of atoms building the coordination variable. ",&
            usage="KINDS_TO {integer} {integer} ..",required=.TRUE.,&
            n_var=-1, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

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


!!****f* input_cp2k/create_colvar_torsion_section *
!!
!!   NAME
!!     create_colvar_torsion_section
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
  SUBROUTINE create_colvar_torsion_section(section,error,description)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: description
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_colvar_torsion_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), pointer :: keyword    

    failure=.FALSE.

    IF (.NOT.failure) THEN
       CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
       IF (PRESENT(description)) THEN
          CALL section_create(section,name="torsion",&
               description=description,&
               n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
               error=error)
       ELSE
          CALL section_create(section,name="torsion",&
               description="Section to define the torsion as a collective variables.",&
               n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
               error=error)
       END IF
        
       NULLIFY(keyword)

       CALL keyword_create(keyword, name="ATOMS",&
            description="Specifies the atom indexes defining the torsion.",&
            usage="ATOMS {integer} {integer} {integer} {integer}",required=.TRUE.,&
            n_var=4, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_colvar_torsion_section

!!****f* input_cp2k/create_colvar_angle_section *
!!
!!   NAME
!!     create_colvar_angle_section
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
  SUBROUTINE create_colvar_angle_section(section,error,description)
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: description
    
    
    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_colvar_angle_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), pointer :: keyword    

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

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="ATOMS",&
            description="Specifies the atom indexes defining the angle.",&
            usage="ATOMS {integer} {integer} {integer}",required=.TRUE.,&
            n_var=3, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

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
    TYPE(section_type), POINTER :: section
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: description

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='create_colvar_dist_section',&
         routineP=moduleN//':'//routineN
    TYPE(keyword_type), pointer :: keyword    

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

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="ATOMS",&
            description="Specifies the atom indexes defining the distance.",&
            usage="ATOMS {integer} {integer}",required=.TRUE.,&
            n_var=2, type_of_var=integer_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_colvar_dist_section
END MODULE input_cp2k_colvar
