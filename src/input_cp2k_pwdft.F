!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Creates the PW section of the input
!> \par History
!>      07.2018 created
!> \author JHU
! **************************************************************************************************
MODULE input_cp2k_pwdft
   USE input_keyword_types,             ONLY: keyword_create,&
                                              keyword_release,&
                                              keyword_type
   USE input_section_types,             ONLY: section_add_keyword,&
                                              section_add_subsection,&
                                              section_create,&
                                              section_release,&
                                              section_type
   USE input_val_types,                 ONLY: lchar_t
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_pwdft'

   PUBLIC :: create_pwdft_section

CONTAINS

! **************************************************************************************************
!> \brief Create the input section for PW calculations using SIRIUS
!> \param section the section to create
!> \par History
!>      07.2018 created
!> \author JHU
! **************************************************************************************************
   SUBROUTINE create_pwdft_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_pwdft_section', &
         routineP = moduleN//':'//routineN

      TYPE(section_type), POINTER                        :: subsection

!   ------------------------------------------------------------------------

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, name="PW_DFT", &
                          description="This section contains all information to run an "// &
                          "SIRIUS PW calculation.", &
                          n_keywords=1, n_subsections=1, repeats=.FALSE.)

      NULLIFY (subsection)
      CALL create_control_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_parameter_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_solver_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_mixer_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_pwdft_section

! **************************************************************************************************
!> \brief input section for PWDFT control
!> \param section will contain the CONTROL section
!> \author JHU
! **************************************************************************************************
   SUBROUTINE create_control_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_control_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, "CONTROL", &
                          description="Control parameters for SIRIUS in list format.", &
                          n_keywords=1, repeats=.FALSE.)

      NULLIFY (keyword)
      CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_", &
                          description="List of parameters:"// &
                          '<p><tt>"PARAMETER" : value,<tt><p>', &
                          repeats=.TRUE., usage="{{String}}", &
                          type_of_var=lchar_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_control_section

! **************************************************************************************************
!> \brief input section for PWDFT parameters
!> \param section will contain the PARAMETER section
!> \author JHU
! **************************************************************************************************
   SUBROUTINE create_parameter_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_parameter_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, "PARAMETER", &
                          description="Parameters for SIRIUS in list format.", &
                          n_keywords=1, repeats=.FALSE.)

      NULLIFY (keyword)
      CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_", &
                          description="List of parameters:"// &
                          '<p><tt>"PARAMETER" : value,<tt><p>', &
                          repeats=.TRUE., usage="{{String}}", &
                          type_of_var=lchar_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_parameter_section

! **************************************************************************************************
!> \brief input section for PWDFT solvers
!> \param section will contain the SOLVER section
!> \author JHU
! **************************************************************************************************
   SUBROUTINE create_solver_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_solver_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, "SOLVER", &
                          description="Solver parameters for SIRIUS in list format.", &
                          n_keywords=1, repeats=.FALSE.)

      NULLIFY (keyword)
      CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_", &
                          description="List of parameters:"// &
                          '<p><tt>"PARAMETER" : value,<tt><p>', &
                          repeats=.TRUE., usage="{{String}}", &
                          type_of_var=lchar_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_solver_section

! **************************************************************************************************
!> \brief input section for PWDFT density mixing
!> \param section will contain the MIXER section
!> \author JHU
! **************************************************************************************************
   SUBROUTINE create_mixer_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_mixer_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, "MIXER", &
                          description="Mixer parameters for SIRIUS in list format.", &
                          n_keywords=1, repeats=.FALSE.)

      NULLIFY (keyword)
      CALL keyword_create(keyword, name="_DEFAULT_KEYWORD_", &
                          description="List of parameters:"// &
                          '<p><tt>"PARAMETER" : value,<tt><p>', &
                          repeats=.TRUE., usage="{{String}}", &
                          type_of_var=lchar_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_mixer_section

! **************************************************************************************************

END MODULE input_cp2k_pwdft
