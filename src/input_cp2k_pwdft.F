!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
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
   USE ISO_C_BINDING, ONLY: C_LOC
   USE SIRIUS, ONLY: &
      sirius_option_get, &
      sirius_option_get_section_length, sirius_option_get_info, &
      SIRIUS_INTEGER_TYPE, SIRIUS_NUMBER_TYPE, SIRIUS_STRING_TYPE, &
      SIRIUS_LOGICAL_TYPE, SIRIUS_ARRAY_TYPE, SIRIUS_INTEGER_ARRAY_TYPE, SIRIUS_LOGICAL_ARRAY_TYPE, &
      SIRIUS_NUMBER_ARRAY_TYPE, SIRIUS_STRING_ARRAY_TYPE, string_f2c
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

      CALL create_sirius_section(subsection, 'settings')
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_sirius_section(subsection, 'mixer')
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_sirius_section(subsection, 'iterative_solver')
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      !
      ! uncomment these lines when nlcg is officialy supported in cp2k
      !

      !      CALL create_sirius_section(subsection, 'nlcg')
      !      CALL section_add_subsection(section, subsection)
      !      CALL section_release(subsection)

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

      INTEGER                                            :: length

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL sirius_option_get_section_length(TRIM(ADJUSTL(section_name)), length)

      CALL section_create(section, __LOCATION__, &
                          name=TRIM(ADJUSTL(section_name)), &
                          description=TRIM(section_name)//" section", &
                          n_subsections=0, &
                          n_keywords=length, &
                          repeats=.FALSE.)

      CALL fill_in_section(section, TRIM(ADJUSTL(section_name)))
   END SUBROUTINE create_sirius_section

! **************************************************************************************************
!> \brief ...
!> \param section ...
!> \param section_name ...
! **************************************************************************************************
   SUBROUTINE fill_in_section(section, section_name)
      TYPE(section_type), POINTER                        :: section
      CHARACTER(len=*), INTENT(in)                       :: section_name

      CHARACTER(len=128)                                 :: name
      CHARACTER(len=128), TARGET                         :: possible_values(1:16)
      CHARACTER(len=4096)                                :: description, usage
      INTEGER                                            :: ctype, enum_i_val(1:16), enum_length, i, &
                                                            j, length, num_possible_values
      INTEGER, ALLOCATABLE, DIMENSION(:), TARGET         :: ivec
      INTEGER, TARGET                                    :: dummy_i
      LOGICAL                                            :: lvecl(1:16)
      LOGICAL(4), ALLOCATABLE, DIMENSION(:), TARGET      :: lvec
      LOGICAL(4), TARGET                                 :: dummy_l
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:), TARGET   :: rvec
      REAL(kind=dp), TARGET                              :: dummy_r
      TYPE(keyword_type), POINTER                        :: keyword

      ALLOCATE (ivec(1:16))
      ALLOCATE (rvec(1:16))
      ALLOCATE (lvec(1:16))
      CALL sirius_option_get_section_length(section_name, length)
      DO i = 1, length
         NULLIFY (keyword)
         name = ''
         description = ''
         usage = ''
         CALL sirius_option_get_info(section_name, &
                                     i, &
                                     name, &
                                     128, &
                                     ctype, &
                                     num_possible_values, &
                                     enum_length, &
                                     description, &
                                     4096, &
                                     usage, &
                                     4096)
         ! description and usage are ignored here
         ! it is a minor inconvenience from the api.

         name = TRIM(ADJUSTL(name))
         ! I exclude these three keywords because one of them is for debugging
         ! purpose the other are replaced by a dedicated call in cp2k
         !
         ! Moreover xc_functionals would need a special treatment.

         IF ((name /= 'xc_functionals') .AND. (name /= 'memory_usage') .AND. (name /= 'vk')) THEN
            !     we need to null char since SIRIUS interface is basically C
            SELECT CASE (ctype)
            CASE (SIRIUS_INTEGER_TYPE)
               CALL sirius_option_get(section_name, name, ctype, C_LOC(dummy_i))
               CALL keyword_create(keyword, __LOCATION__, &
                                   name=TRIM(name), &
                                   description=TRIM(ADJUSTL(description)), &
                                   !                                   usage=TRIM(ADJUSTL(usage)), &
                                   type_of_var=integer_t, &
                                   repeats=.FALSE., &
                                   default_i_val=dummy_i)
               CALL section_add_keyword(section, keyword)
               CALL keyword_release(keyword)
            CASE (SIRIUS_NUMBER_TYPE)
               CALL sirius_option_get(section_name, name, ctype, C_LOC(dummy_r))
               CALL keyword_create(keyword, __LOCATION__, &
                                   name=name, &
                                   description=TRIM(ADJUSTL(description)), &
                                   !                                   usage=TRIM(ADJUSTL(usage)), &
                                   type_of_var=real_t, &
                                   repeats=.FALSE., &
                                   default_r_val=dummy_r)
               CALL section_add_keyword(section, keyword)
               CALL keyword_release(keyword)
            CASE (SIRIUS_LOGICAL_TYPE)
               CALL sirius_option_get(section_name, name, ctype, C_LOC(dummy_l))
               IF (dummy_l) THEN
                  CALL keyword_create(keyword, __LOCATION__, &
                                      name=name, &
                                      description=TRIM(ADJUSTL(description)), &
                                      !                                      usage=TRIM(ADJUSTL(usage)), &
                                      type_of_var=logical_t, &
                                      repeats=.FALSE., &
                                      default_l_val=.TRUE., &
                                      lone_keyword_l_val=.TRUE.)
               ELSE
                  CALL keyword_create(keyword, __LOCATION__, &
                                      name=name, &
                                      description=TRIM(ADJUSTL(description)), &
                                      !                                      usage=TRIM(ADJUSTL(usage)), &
                                      type_of_var=logical_t, &
                                      repeats=.FALSE., &
                                      default_l_val=.FALSE., &
                                      lone_keyword_l_val=.TRUE.)
               END IF
               CALL section_add_keyword(section, keyword)
               CALL keyword_release(keyword)
            CASE (SIRIUS_STRING_TYPE)
               IF (enum_length >= 1) THEN
                  DO j = 1, enum_length
                     possible_values(j) = ''
                     CALL sirius_option_get(section_name, name, ctype, C_LOC(possible_values(j)), max_length=128, enum_idx=j)
                     enum_i_val(j) = j
                     possible_values(j) = TRIM(ADJUSTL(possible_values(j)))
                  END DO

                  IF (enum_length > 1) THEN
                     CALL keyword_create(keyword, __LOCATION__, &
                                         name=name, &
                                         description=TRIM(ADJUSTL(description)), &
                                         !                                      usage=TRIM(ADJUSTL(usage)), &
                                         repeats=.FALSE., &
                                         enum_i_vals=enum_i_val(1:enum_length), &
                                         enum_c_vals=possible_values(1:enum_length), &
                                         default_i_val=1)
                  ELSE
                     CALL keyword_create(keyword, __LOCATION__, &
                                         name=name, &
                                         description=TRIM(ADJUSTL(description)), &
                                         !                                      usage=TRIM(ADJUSTL(usage)), &
                                         type_of_var=char_t, &
                                         default_c_val=possible_values(1), &
                                         repeats=.FALSE.)
                  END IF
               ELSE
                  CALL keyword_create(keyword, __LOCATION__, &
                                      name=name, &
                                      description=TRIM(ADJUSTL(description)), &
                                      !                                      usage=TRIM(ADJUSTL(usage)), &
                                      type_of_var=char_t, &
                                      default_c_val='', &
                                      repeats=.FALSE.)
               END IF
               CALL section_add_keyword(section, keyword)
               CALL keyword_release(keyword)
            CASE (SIRIUS_INTEGER_ARRAY_TYPE)
               CALL sirius_option_get(section_name, name, ctype, C_LOC(ivec(1)), max_length=16)

               IF (num_possible_values .EQ. 0) THEN
                  CALL keyword_create(keyword, __LOCATION__, &
                                      name=name, &
                                      description=TRIM(ADJUSTL(description)), &
                                      type_of_var=integer_t, &
                                      n_var=-1, &
                                      repeats=.FALSE.)
               ELSE
                  CALL keyword_create(keyword, __LOCATION__, &
                                      name=name, &
                                      description=TRIM(ADJUSTL(description)), &
                                      type_of_var=integer_t, &
                                      repeats=.FALSE., &
                                      n_var=num_possible_values, &
                                      default_i_vals=ivec(1:num_possible_values))
               END IF
               CALL section_add_keyword(section, keyword)
               CALL keyword_release(keyword)
            CASE (SIRIUS_LOGICAL_ARRAY_TYPE)
               CALL sirius_option_get(section_name, name, ctype, C_LOC(lvec(1)), max_length=16)
               DO j = 1, num_possible_values
                  lvecl(j) = lvec(j)
               END DO
               IF (num_possible_values > 0) THEN
                  CALL keyword_create(keyword, __LOCATION__, &
                                      name=name, &
                                      description=TRIM(ADJUSTL(description)), &
                                      !usage=TRIM(ADJUSTL(usage)), &
                                      type_of_var=logical_t, &
                                      repeats=.FALSE., &
                                      n_var=num_possible_values, &
                                      default_l_vals=lvecl(1:num_possible_values))
               ELSE
                  CALL keyword_create(keyword, __LOCATION__, &
                                      name=name, &
                                      description=TRIM(ADJUSTL(description)), &
                                      !usage=TRIM(ADJUSTL(usage)), &
                                      type_of_var=logical_t, &
                                      repeats=.FALSE., &
                                      n_var=-1)
               END IF
               CALL section_add_keyword(section, keyword)
               CALL keyword_release(keyword)
            CASE (SIRIUS_NUMBER_ARRAY_TYPE)
               CALL sirius_option_get(section_name, name, ctype, C_LOC(rvec(1)), max_length=16)

               IF (num_possible_values .EQ. 0) THEN
                  CALL keyword_create(keyword, __LOCATION__, &
                                      name=name, &
                                      description=TRIM(ADJUSTL(description)), &
                                      !                                   usage=TRIM(ADJUSTL(usage)), &
                                      type_of_var=real_t, &
                                      repeats=.FALSE., &
                                      n_var=-1)
               ELSE
                  CALL keyword_create(keyword, __LOCATION__, &
                                      name=name, &
                                      description=TRIM(ADJUSTL(description)), &
                                      !     usage=TRIM(ADJUSTL(usage)), &
                                      type_of_var=real_t, &
                                      repeats=.FALSE., &
                                      n_var=num_possible_values, &
                                      default_r_vals=rvec(1:num_possible_values))
               END IF
               CALL section_add_keyword(section, keyword)
               CALL keyword_release(keyword)
            CASE default
            END SELECT
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
