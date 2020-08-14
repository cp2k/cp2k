!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!
!
! **************************************************************************************************
! > \brief Creates the PW section of the input
! > \par History
! >      07.2018 created
! > \author JHU
! **************************************************************************************************

MODULE input_cp2k_pwdft
#if defined(__SIRIUS)
   USE SIRIUS, ONLY: &
      sirius_option_get_description_usage, sirius_option_get_double, sirius_option_get_int, &
      sirius_option_get_length, sirius_option_get_logical, sirius_option_get_name_and_type, &
      sirius_option_get_number_of_possible_values, sirius_option_get_string, &
      sirius_option_string_get_value
#endif
   USE input_keyword_types, ONLY: keyword_create, &
                                  keyword_release, &
                                  keyword_type
   USE input_section_types, ONLY: section_add_keyword, &
                                  section_add_subsection, &
                                  section_create, &
                                  section_release, &
                                  section_type
   USE input_val_types, ONLY: char_t, &
                              integer_t, &
                              lchar_t, &
                              logical_t, &
                              real_t
   USE cp_output_handling, ONLY: add_last_numeric, &
                                 cp_print_key_section_create, &
                                 debug_print_level, &
                                 high_print_level, &
                                 low_print_level, &
                                 medium_print_level, &
                                 silent_print_level
   USE kinds, ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_pwdft'

   PUBLIC :: create_pwdft_section

CONTAINS

#if defined(__SIRIUS)
! **************************************************************************************************
!> \brief Create the input section for PW calculations using SIRIUS
!> \param section the section to create
!> \par History
!>      07.2018 created
!> \author JHU
! **************************************************************************************************
   SUBROUTINE create_pwdft_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_pwdft_section'

      TYPE(section_type), POINTER                        :: subsection

!     ------------------------------------------------------------------------

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="PW_DFT", &
                          description="DFT calculation using plane waves basis can be set in this section. "// &
                          "The backend called SIRIUS, computes the basic properties of the system, "// &
                          "such as ground state, forces and stresses tensors which can be used by "// &
                          "cp2k afterwards. The engine has all these features build-in, support of "// &
                          "pseudo-potentials and full-potentials, spin-orbit coupling, collinear and "// &
                          "non collinear magnetism, Hubbard correction, all exchange functionals "// &
                          "supported by libxc and Van der Waals corrections (libvdwxc).")

      NULLIFY (subsection)
      CALL create_sirius_section(subsection, 'control')
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_sirius_section(subsection, 'parameters')
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_sirius_section(subsection, 'mixer')
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_sirius_section(subsection, 'iterative_solver')
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_print_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_pwdft_section

! **************************************************************************************************
!> \brief input section for PWDFT control
!> \param section will contain the CONTROL section
!> \param section_name ...
!> \author JHU
! **************************************************************************************************
   SUBROUTINE create_sirius_section(section, section_name)
      TYPE(section_type), POINTER                        :: section
      CHARACTER(len=*), INTENT(in)                       :: section_name

      CHARACTER(len=*), PARAMETER :: routineN = 'create_sirius_section'

      INTEGER                                            :: length

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL sirius_option_get_length(TRIM(ADJUSTL(section_name))//CHAR(0), length)

      CALL section_create(section, __LOCATION__, &
                          name=TRIM(ADJUSTL(section_name)), &
                          description=TRIM(section_name)//" section", &
                          n_subsections=0, &
                          n_keywords=length, &
                          repeats=.FALSE.)

      CALL fill_in_section(section, TRIM(ADJUSTL(section_name))//CHAR(0))
   END SUBROUTINE create_sirius_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
!> \param section_name ...
! **************************************************************************************************
   SUBROUTINE fill_in_section(section, section_name)
      TYPE(section_type), POINTER                        :: section
      CHARACTER(len=*), INTENT(in)                       :: section_name

      CHARACTER(len=*), PARAMETER                        :: routineN = 'fill_in_section'

      CHARACTER(len=128)                                 :: name, name1, possible_values(1:16)
      CHARACTER(len=512)                                 :: default_string_val, description, usage
      INTEGER                                            :: ctype, dummy_i, enum_i_val(1:16), i, j, &
                                                            length, num_possible_values, vec_length
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ivec
      LOGICAL                                            :: lvecl(1:16)
      LOGICAL(1)                                         :: dummy_l
      LOGICAL(1), ALLOCATABLE, DIMENSION(:)              :: lvec
      REAL(kind=dp)                                      :: dummy_r
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: rvec
      TYPE(keyword_type), POINTER                        :: keyword

      ALLOCATE (ivec(1:16))
      ALLOCATE (rvec(1:16))
      ALLOCATE (lvec(1:16))
      CALL sirius_option_get_length(section_name, length)
      DO i = 0, length - 1
         NULLIFY (keyword)
         name = CHAR(0)
         ! return a non null terminated string. Stupid fortran does not understand the \0 terminated string when comparing things
         CALL sirius_option_get_name_and_type(section_name, i, name, ctype)

!     do not invert these two lines
         name1 = TRIM(ADJUSTL(name))

!     we need to null char since SIRIUS interface is basically C
         name = TRIM(ADJUSTL(name))//CHAR(0)
         ! I exclude these two keywords because one of them is for debugging
         ! purpose the other one is useless and replaced by a dedicated in cp2k
         !
         ! Moreover xc_functionals would need a special treatment.
         IF ((name1 /= 'xc_functionals') .AND. (name1 /= 'memory_usage')) THEN
            description = CHAR(0)
            usage = CHAR(0)
            CALL sirius_option_get_description_usage(section_name, name, description, usage)
            SELECT CASE (ctype)
            CASE (1)
               CALL sirius_option_get_int(section_name, name, dummy_i, vec_length)
               CALL keyword_create(keyword, __LOCATION__, &
                                   name=TRIM(name1), &
                                   description=TRIM(ADJUSTL(description)), &
                                   usage=TRIM(ADJUSTL(usage)), &
                                   repeats=.FALSE., &
                                   default_i_val=dummy_i)
            CASE (11)
               CALL sirius_option_get_int(section_name, name, ivec(1), vec_length)
               CALL keyword_create(keyword, __LOCATION__, &
                                   name=name1, &
                                   description=TRIM(ADJUSTL(description)), &
                                   usage=TRIM(ADJUSTL(usage)), &
                                   repeats=.FALSE., &
                                   n_var=vec_length, &
                                   default_i_vals=ivec(1:vec_length))
            CASE (2)
               CALL sirius_option_get_double(section_name, name, dummy_r, vec_length)
               CALL keyword_create(keyword, __LOCATION__, &
                                   name=name1, &
                                   description=TRIM(ADJUSTL(description)), &
                                   usage=TRIM(ADJUSTL(usage)), &
                                   repeats=.FALSE., &
                                   default_r_val=dummy_r)
            CASE (12)
               CALL sirius_option_get_double(section_name, name, rvec(1), vec_length)
               CALL keyword_create(keyword, __LOCATION__, &
                                   name=name1, &
                                   description=TRIM(ADJUSTL(description)), &
                                   usage=TRIM(ADJUSTL(usage)), &
                                   repeats=.FALSE., &
                                   n_var=vec_length, &
                                   default_r_vals=rvec(1:vec_length))
            CASE (3)
               CALL sirius_option_get_logical(section_name, name, dummy_l, vec_length)
               IF (dummy_l) THEN
                  CALL keyword_create(keyword, __LOCATION__, &
                                      name=name1, &
                                      description=TRIM(ADJUSTL(description)), &
                                      usage=TRIM(ADJUSTL(usage)), &
                                      repeats=.FALSE., &
                                      default_l_val=.TRUE., &
                                      lone_keyword_l_val=.TRUE.)
               ELSE
                  CALL keyword_create(keyword, __LOCATION__, &
                                      name=name1, &
                                      description=TRIM(ADJUSTL(description)), &
                                      usage=TRIM(ADJUSTL(usage)), &
                                      repeats=.FALSE., &
                                      default_l_val=.FALSE., &
                                      lone_keyword_l_val=.TRUE.)
               ENDIF
            CASE (13)
               CALL sirius_option_get_logical(section_name, name, lvec(1), vec_length)
               DO j = 1, vec_length
                  lvecl(j) = lvec(j)
               ENDDO
               CALL keyword_create(keyword, __LOCATION__, &
                                   name=name1, &
                                   description=TRIM(ADJUSTL(description)), &
                                   usage=TRIM(ADJUSTL(usage)), &
                                   repeats=.FALSE., &
                                   n_var=vec_length, &
                                   default_l_vals=lvecl(1:vec_length))
            CASE (4)
               !     string need a special treatment because the parameters can only have dedicated values
               default_string_val = CHAR(0)
               CALL sirius_option_get_string(section_name, name, default_string_val)
               default_string_val = TRIM(ADJUSTL(default_string_val))
               CALL sirius_option_get_number_of_possible_values(section_name, name, num_possible_values)
               IF (num_possible_values > 0) THEN
                  DO j = 0, num_possible_values - 1
                     possible_values(j + 1) = CHAR(0)
                     CALL sirius_option_string_get_value(section_name, name, j, possible_values(j + 1))
                     enum_i_val(j + 1) = j
                  END DO
                  CALL keyword_create(keyword, __LOCATION__, &
                                      name=name1, &
                                      description=TRIM(ADJUSTL(description)), &
                                      usage=TRIM(ADJUSTL(usage)), &
                                      repeats=.FALSE., &
                                      enum_i_vals=enum_i_val(1:num_possible_values), &
                                      enum_c_vals=possible_values(1:num_possible_values), &
                                      default_i_val=0)
               ELSE
                  CALL keyword_create(keyword, __LOCATION__, &
                                      name=name1, &
                                      description=TRIM(ADJUSTL(description)), &
                                      usage=TRIM(ADJUSTL(usage)), &
                                      type_of_var=char_t, &
                                      default_c_val="", &
                                      repeats=.FALSE.)
               END IF
            CASE default
            END SELECT
            CALL section_add_keyword(section, keyword)
            CALL keyword_release(keyword)
         END IF
      END DO
   END SUBROUTINE fill_in_section

! **************************************************************************************************
!> \brief Create the print section for sirius
!> \param section the section to create
!> \author jgh
! **************************************************************************************************
   SUBROUTINE create_print_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_print_section'

      TYPE(section_type), POINTER                        :: print_key

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="PRINT", &
                          description="Section of possible print options in PW_DFT code.", &
                          n_keywords=0, n_subsections=1, repeats=.FALSE.)

      NULLIFY (print_key)
      CALL create_dos_section(print_key)
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

   END SUBROUTINE create_print_section

! **************************************************************************************************
!> \brief ...
!> \param print_key ...
! **************************************************************************************************
   SUBROUTINE create_dos_section(print_key)

      TYPE(section_type), POINTER                        :: print_key

      CHARACTER(len=*), PARAMETER :: routineN = 'create_dos_section'

      TYPE(keyword_type), POINTER                        :: keyword

      NULLIFY (keyword)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "DOS", &
                                       description="Print Density of States (DOS) (only available states from SCF)", &
                                       print_level=debug_print_level, common_iter_levels=1, filename="")

      CALL keyword_create(keyword, __LOCATION__, name="APPEND", &
                          description="Append the DOS obtained at different iterations to the output file."// &
                          "By default the file is overwritten", &
                          usage="APPEND", default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DELTA_E", &
                          description="Histogramm energy spacing.", &
                          usage="DELTA_E 0.0005", type_of_var=real_t, default_r_val=0.001_dp)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_dos_section

#else
! **************************************************************************************************
!> \brief ...
!> \param section ...
! **************************************************************************************************
   SUBROUTINE create_pwdft_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_pwdft_section', &
         routineP = moduleN//':'//routineN

      CPASSERT(.NOT. ASSOCIATED(section))

      CALL section_create(section, __LOCATION__, name="PW_DFT", &
                          description="This section contains all information to run an "// &
                          "SIRIUS PW calculation.", &
                          n_subsections=0, &
                          repeats=.FALSE.)

   END SUBROUTINE create_pwdft_section

#endif

END MODULE input_cp2k_pwdft
