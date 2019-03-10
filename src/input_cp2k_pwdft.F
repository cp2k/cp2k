!--------------------------------------------------------------------------------------------------!
!     CP2K: A general program to perform molecular dynamics simulations                              !
!     Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

!     **************************************************************************************************
!     > \brief Creates the PW section of the input
!     > \par History
!     >      07.2018 created
!     > \author JHU
!     **************************************************************************************************
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
      !USE mpi,                             ONLY: MPI_COMM_WORLD
      USE SIRIUS, ONLY: &
           sirius_option_set_string, sirius_option_set_logical, sirius_option_set_double, &
           sirius_option_set_int, sirius_option_get_string, sirius_option_get_logical, &
           sirius_option_get_double, sirius_option_get_int, sirius_option_get_name_and_type, &
           sirius_option_get_length, sirius_create_context, sirius_option_get_number_of_possible_values, &
           sirius_option_string_get_value, sirius_option_get_description_usage
      USE SIRIUS_INTERFACE, ONLY : cp_sirius_init, cp_sirius_finalize
      USE ISO_C_BINDING, ONLY: &
      C_ASSOCIATED, C_BOOL, C_DOUBLE, C_FUNPTR, C_F_POINTER, C_INT, C_NULL_FUNPTR, C_NULL_PTR, &
      C_PTR
      USE kinds, ONLY: dp, default_string_length

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
      CHARACTER(len=32)                                  :: section_name

!     ------------------------------------------------------------------------

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, name="PW_DFT", &
      description="This section contains all information to run an "// &
      "SIRIUS PW calculation.", &
      n_subsections=4, &
      repeats=.FALSE.)

      NULLIFY (subsection)
      section_name= ''
      section_name = 'control'
      CALL create_sirius_section(subsection, section_name)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)
      section_name= ''
      section_name = 'parameters'

      CALL create_sirius_section(subsection, section_name)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)
      section_name = ''
      section_name = 'mixer'

      CALL create_sirius_section(subsection, section_name)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)
      section_name = ''
      section_name = 'iterative_solver'

      CALL create_sirius_section(subsection, section_name)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      END SUBROUTINE create_pwdft_section

! **************************************************************************************************
!> \brief input section for PWDFT control
!> \param section will contain the CONTROL section
!> \author JHU
! **************************************************************************************************
      SUBROUTINE create_sirius_section(section, section_name)
      TYPE(section_type), POINTER                        :: section
      CHARACTER(len=32), intent(in)                      :: section_name
      CHARACTER(len=*), PARAMETER :: routineN = 'create_sirius_section', &
      routineP = moduleN//':'//routineN
      integer                                            :: length
      CPASSERT(.NOT. ASSOCIATED(section))
      CALL sirius_option_get_length(TRIM(adjustl(section_name))//char(0), length)
      CALL section_create(section, &
      name=TRIM(ADJUSTL(section_name)), &
      description=TRIM(section_name)//" section", &
      n_subsections=0, &
      n_keywords=length, &
      repeats=.FALSE.)

      CALL fill_in_section(section, TRIM(adjustl(section_name))//char(0))
    END SUBROUTINE create_sirius_section


! **************************************************************************************************
      SUBROUTINE fill_in_section(section, section_name)
      TYPE(section_type), POINTER                        :: section
      CHARACTER(len=32), intent(in)                      :: section_name
      CHARACTER(len=*), PARAMETER :: routineN = 'fill_in_section', &
           routineP = moduleN//':'//routineN


      TYPE(keyword_type), POINTER                        :: keyword
      INTEGER                                            :: length, i, ctype, num_possible_values, vec_length, j
      CHARACTER(len=512)                                 :: description, usage, default_string_val
      CHARACTER(len=128)                                 :: name
      CHARACTER(len=128)                                 :: name1
      CHARACTER(len=128)                                 :: possible_values(1:16)
      INTEGER                                            :: dummy_i
      REAL(kind=dp)                                      :: dummy_r
      LOGICAL(1)                                         :: dummy_l
      INTEGER                                            :: enum_i_val(1:16)
      INTEGER, DIMENSION(:), allocatable                 :: ivec
      LOGICAL(1), DIMENSION(:), allocatable              :: lvec
      LOGICAL                                            :: lvecl(1:16)
      REAL(kind=dp), DIMENSION(:), allocatable           :: rvec


      allocate(ivec(1:16))
      allocate(rvec(1:16))
      allocate(lvec(1:16))
#if defined(__SIRIUS)
      CALL sirius_option_get_length(section_name, length)
      DO i = 0, length - 1
         NULLIFY(keyword)
         name = char(0)
      ! return a non null terminated string. Stupid fortran does not understand the \0 terminated string when comparing things
      CALL sirius_option_get_name_and_type(section_name, i, name, ctype)

!     do not invert these two lines
         name1 = TRIM(ADJUSTL(name))

!     we need to null char since SIRIUS interface is basically C
         name = TRIM(ADJUSTL(name))//char(0)
         description = char(0)
         usage = char(0)
         CALL sirius_option_get_description_usage(section_name, name, description, usage)
         select case(ctype)
      case (1)
         call sirius_option_get_int(section_name, name, dummy_i, vec_length);
         call keyword_create(keyword, &
         name=TRIM(name1), &
         description=TRIM(ADJUSTL(description)), &
         usage=TRIM(ADJUSTL(usage)), &
         repeats=.FALSE., &
         default_i_val=dummy_i)
      case (11)
         call sirius_option_get_int(section_name, name, ivec(1), vec_length);
         call keyword_create(keyword, &
         name=name1, &
         description=TRIM(ADJUSTL(description)), &
         usage=TRIM(ADJUSTL(usage)), &
         repeats=.FALSE., &
         n_var=vec_length, &
         default_i_vals=ivec(1:vec_length))
      case (2)
         call sirius_option_get_double(section_name, name, dummy_r, vec_length);
         call keyword_create(keyword, &
         name=name1, &
         description=TRIM(ADJUSTL(description)), &
         usage=TRIM(ADJUSTL(usage)), &
         repeats=.FALSE., &
         default_r_val=dummy_r)
      case (12)
         call sirius_option_get_double(section_name, name, rvec(1), vec_length);
         call keyword_create(keyword, &
         name=name1, &
         description=TRIM(ADJUSTL(description)), &
         usage=TRIM(ADJUSTL(usage)), &
         repeats=.FALSE., &
         n_var=vec_length, &
         default_r_vals=rvec(1:vec_length))
      case (3)
         call sirius_option_get_logical(section_name, name, dummy_l, vec_length);
         if (dummy_l) then
            call keyword_create(keyword, &
            name=name1, &
                 description=TRIM(ADJUSTL(description)), &
                 usage=TRIM(ADJUSTL(usage)), &
                 repeats=.FALSE., &
                 default_l_val=.true., &
                 lone_keyword_l_val=.true.)
         else
            call keyword_create(keyword, &
            name=name1, &
            description=TRIM(ADJUSTL(description)), &
            usage=TRIM(ADJUSTL(usage)), &
            repeats=.FALSE., &
            default_l_val=.false., &
            lone_keyword_l_val=.true.)
         endif
      case (13)
         CALL sirius_option_get_logical(section_name, name, lvec(1), vec_length);
         do j = 1, vec_length
            lvecl(j) = lvec(j)
         enddo
         call keyword_create(keyword, &
         name=name1, &
         description=TRIM(ADJUSTL(description)), &
         usage=TRIM(ADJUSTL(usage)), &
         repeats=.FALSE., &
         n_var=vec_length, &
         default_l_vals=lvecl(1:vec_length))
      case (4)
         !     string need a special treatment because the parameters can only have dedicated values
         default_string_val = char(0)
         call sirius_option_get_string(section_name, name, default_string_val)
         default_string_val = trim(adjustl(default_string_val))
         call sirius_option_get_number_of_possible_values(section_name, name, num_possible_values)
         if (num_possible_values > 0) then
            do j = 0, num_possible_values - 1
               possible_values(j + 1) = char(0)
               CALL sirius_option_string_get_value(section_name, name, j, possible_values(j + 1))
               enum_i_val(j+1) = j
            end do
            call keyword_create(keyword, &
            name=name1, &
            description=TRIM(ADJUSTL(description)), &
            usage=TRIM(ADJUSTL(usage)), &
            repeats=.FALSE., &
            enum_i_vals=enum_i_val(1:num_possible_values), &
            enum_c_vals=possible_values(1:num_possible_values), &
            default_i_val=0)
         else
            call keyword_create(keyword, &
            name=name1, &
            description=TRIM(ADJUSTL(description)), &
            usage=TRIM(ADJUSTL(usage)), &
            repeats=.FALSE.)
         end if
      case default
      end select
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)
      END DO
#endif
      END SUBROUTINE fill_in_section
      END MODULE input_cp2k_pwdft
