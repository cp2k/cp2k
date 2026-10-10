!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Creates and parses an input section to define a generic function
!> \note
!>      This module handles the input definition and retrieval, while the task of actually parsing
!>      and evaluating a function is handled by common/fparser module. The general workflow is:
!>      1. In the module that creates input sections, use the subroutine to add the subsection
!>         &GENERIC_FUNCTION with keywords FUNCTION, PARAMETERS, UNITS, VALUES and (optionally)
!>         VARIABLES. If a list of symbols has been specified (like input_variables=["X", "Y", "Z"])
!>         then they will be reserved as the variables and the VARIABLES keyword will not appear.
!>         If no such list is given then the VARIABLES keyword requesting input of symbol(s) for
!>         the variable(s) will appear, with no restriction on the number of symbols.
!>            NULLIFY(subsection)
!>            CALL create_function_section(subsection, "GENERIC_FUNCTION", "FUNCTION", &
!>                                         "This subsection defines a generic function.")
!>            ! Add extra keywords here if needed
!>            ! NULLIFY (keyword)
!>            ! CALL keyword_create(keyword, __LOCATION__, ...)
!>            ! CALL section_add_keyword(print_key, keyword)
!>            ! CALL keyword_release(keyword)
!>            ! ...
!>            ! Done adding extra keywords
!>            CALL section_add_subsection(section, subsection)
!>            CALL section_release(subsection)
!>      2. Define a type that stores information about the function.
!>            TYPE my_func_type
!>               CHARACTER(LEN=default_path_length) :: function = ""
!>               CHARACTER(LEN=default_string_length), DIMENSION(:), POINTER :: parameters => NULL()
!>               REAL(KIND=dp), DIMENSION(:), POINTER :: values => NULL()
!>            END TYPE my_func_type
!>      3. In the module that parses input, use the subroutine to get and save the definitions. The
!>         usage of input_variables should be consistent with above. Note that the variables and the
!>         parameters are all saved in my_func%parameters and my_func%values; with n_var variables
!>         and n_par parameters, both my_func%parameters and my_func%values will have a dimension of
!>         (n_var + n_par), variables first.
!>            TYPE(my_func_type), POINTER :: my_func
!>            NULLIFY (my_func)
!>            CALL parse_function_section(gen_section, "FUNCTION", my_func%function, &
!>                                        my_func%parameters, my_func%values)
!>      4. Use the subroutines from the fparser module to evaluate the function.
!>            CALL initf(1)
!>            CALL parsef(1, TRIM(my_func%function), my_func%parameters)
!>            ! Substitute new values for the variables if needed
!>            ! my_func%values(1) = ...
!>            ! Done substituting
!>            result = evalf(1, my_func%values)
!>            CPASSERT(EvalErrType <= 0)
!>            CALL finalizef()
!> \par History
!>      10.2026 created
!> \author Zilong He
! **************************************************************************************************
MODULE input_cp2k_function
   USE cp_units,                        ONLY: cp_unit_to_cp2k
   USE fparser,                         ONLY: docf
   USE input_keyword_types,             ONLY: keyword_create,&
                                              keyword_release,&
                                              keyword_type
   USE input_section_types,             ONLY: section_add_keyword,&
                                              section_create,&
                                              section_type,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE input_val_types,                 ONLY: char_t,&
                                              lchar_t,&
                                              real_t
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE memory_utilities,                ONLY: reallocate
   USE string_utilities,                ONLY: compress,&
                                              integer_to_string,&
                                              newline
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_function'

   CHARACTER(len=*), PARAMETER, PRIVATE :: name_kw_parameters = "PARAMETERS", &
                                           name_kw_units = "UNITS", &
                                           name_kw_values = "VALUES", &
                                           name_kw_variables = "VARIABLES"

   PUBLIC :: create_function_section, parse_function_section

CONTAINS

! **************************************************************************************************
!> \brief Create the input section for a generic function
!> \param gen_section the section to be created for the generic function
!> \param section_name the name of the subsection for the generic function
!> \param func_name the name of the keyword for the function
!> \param description a custom description for the subsection
!> \param repeats whether the section can be repeated
!> \param input_variables a list of symbols that will be parsed as variables
!> \par History
!>      10.2026 created from scratch matching argument names with parse_function_section
! **************************************************************************************************
   SUBROUTINE create_function_section(gen_section, section_name, func_name, description, &
                                      repeats, input_variables)
      TYPE(section_type), POINTER                        :: gen_section
      CHARACTER(LEN=*), INTENT(IN)                       :: section_name, func_name, description
      LOGICAL, INTENT(IN), OPTIONAL                      :: repeats
      CHARACTER(LEN=*), DIMENSION(:), OPTIONAL           :: input_variables

      CHARACTER(LEN=*), PARAMETER :: default_values_unit_str = "internal_cp2k"

      CHARACTER(LEN=:), ALLOCATABLE                      :: input_variables_info, my_description
      CHARACTER(LEN=default_string_length)               :: my_size_default, my_size_var
      INTEGER                                            :: i, mydim
      LOGICAL                                            :: has_kw_variables, my_repeats
      TYPE(keyword_type), POINTER                        :: keyword

      my_repeats = .FALSE.
      IF (PRESENT(repeats)) my_repeats = repeats
      has_kw_variables = .FALSE.
      my_description = TRIM(description)
      mydim = 0
      CALL integer_to_string(default_string_length, my_size_default)

      IF (PRESENT(input_variables)) THEN
         has_kw_variables = .TRUE.
         mydim = SIZE(input_variables)
         CALL integer_to_string(mydim, my_size_var)
         input_variables_info = "The following "//TRIM(my_size_var)//" symbol(s) will be reserved "// &
                                "as the variable(s) to be substituted for internal evaluation: "
         DO i = 1, SIZE(input_variables)
            input_variables_info = input_variables_info//"`"//TRIM(input_variables(i))//"`"
            IF (i < SIZE(input_variables)) THEN
               input_variables_info = input_variables_info//", "
            ELSE
               input_variables_info = input_variables_info//". "
            END IF
         END DO
         input_variables_info = input_variables_info//"They cannot be used as parameters."
         my_description = my_description//newline//newline//input_variables_info
      END IF

      NULLIFY (keyword)
      CPASSERT(.NOT. ASSOCIATED(gen_section))

      CALL section_create(gen_section, __LOCATION__, name=TRIM(section_name), &
                          description=my_description, &
                          n_keywords=0, n_subsections=0, repeats=my_repeats)

      CALL keyword_create(keyword, __LOCATION__, name=TRIM(func_name), &
                          description=docf(), &
                          usage=TRIM(func_name)//" A*X+B*(2/(1+EXP(-0.001*Y))-1)", &
                          type_of_var=lchar_t, n_var=1)
      CALL section_add_keyword(gen_section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name=name_kw_parameters, &
                          description="Defines the symbols for the parameter(s) in the functional "// &
                          "form, delimited by whitespaces. Symbols are case sensitive and each is "// &
                          "at most "//TRIM(my_size_default)//" characters long.", &
                          usage=name_kw_parameters//" A B", type_of_var=char_t, &
                          n_var=-1, repeats=.TRUE.)
      CALL section_add_keyword(gen_section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name=name_kw_units, &
                          description="Optionally, defines CP2K unit string for each of the given "// &
                          "parameter value(s). A valid CP2K unit string is one of the strings that "// &
                          "is listed under the References/Units section of the CP2K manual, or a "// &
                          "combination of multiple of these strings with `*`, `/` and `^`. It is "// &
                          "assumed that the corresponding parameter value is specified in this unit "// &
                          "which needs to be converted first to the internal unit during input parsing. "// &
                          "Note that the meaningfulness of the unit is not checked. If omitted, all "// &
                          "values will be treated as already in the internal unit and used literally.", &
                          usage=name_kw_units//" eV*angstrom^-1 angstrom", type_of_var=char_t, &
                          n_var=-1, repeats=.TRUE.)
      CALL section_add_keyword(gen_section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name=name_kw_values, &
                          description="Defines the value(s) of the parameter in the functional form.", &
                          usage=name_kw_values//" 50 1.0E+2", type_of_var=real_t, &
                          n_var=-1, repeats=.TRUE., unit_str=default_values_unit_str)
      CALL section_add_keyword(gen_section, keyword)
      CALL keyword_release(keyword)

      IF (.NOT. has_kw_variables) THEN
         CALL keyword_create(keyword, __LOCATION__, name=name_kw_variables, &
                             description="Defines the symbols for the variable(s) in the functional "// &
                             "form, delimited by whitespaces. Symbols are case sensitive and each is "// &
                             "at most "//TRIM(my_size_default)//" characters long.", &
                             usage=name_kw_variables//" X Y", type_of_var=char_t, &
                             n_var=-1)
         CALL section_add_keyword(gen_section, keyword)
         CALL keyword_release(keyword)
      END IF

   END SUBROUTINE create_function_section

! **************************************************************************************************
!> \brief Reads from the input structure all information for generic functions
!> \param gen_section ...
!> \param func_name ...
!> \param xfunction ...
!> \param parameters ...
!> \param values ...
!> \param var_values ...
!> \param i_rep_sec ...
!> \param input_variables ...
!> \par History
!>      10.2026 created by migrating the get_generic_info in force_fields_util.F to here
! **************************************************************************************************
   SUBROUTINE parse_function_section(gen_section, func_name, xfunction, parameters, values, &
                                     var_values, i_rep_sec, input_variables)
      TYPE(section_vals_type), POINTER                   :: gen_section
      CHARACTER(LEN=*), INTENT(IN)                       :: func_name
      CHARACTER(LEN=default_path_length), INTENT(OUT)    :: xfunction
      CHARACTER(LEN=default_string_length), &
         DIMENSION(:), POINTER                           :: parameters
      REAL(KIND=dp), DIMENSION(:), POINTER               :: values
      REAL(KIND=dp), DIMENSION(:), OPTIONAL, POINTER     :: var_values
      INTEGER, INTENT(IN), OPTIONAL                      :: i_rep_sec
      CHARACTER(LEN=*), DIMENSION(:), OPTIONAL           :: input_variables

      CHARACTER(LEN=default_string_length), &
         DIMENSION(:), POINTER                           :: my_par, my_par_tmp, my_units, &
                                                            my_units_tmp, my_var
      INTEGER                                            :: i, ind, irep, isize, j, mydim, n_par, &
                                                            n_units, n_val, nblank
      LOGICAL                                            :: check
      REAL(KIND=dp), DIMENSION(:), POINTER               :: my_val, my_val_tmp

      NULLIFY (my_var, my_par, my_val, my_par_tmp, my_val_tmp)
      NULLIFY (my_units)
      NULLIFY (my_units_tmp)
      IF (ASSOCIATED(parameters)) THEN
         DEALLOCATE (parameters)
      END IF
      IF (ASSOCIATED(values)) THEN
         DEALLOCATE (values)
      END IF
      irep = 1
      IF (PRESENT(i_rep_sec)) irep = i_rep_sec
      mydim = 0
      CALL section_vals_val_get(gen_section, TRIM(func_name), i_rep_section=irep, c_val=xfunction)
      CALL compress(xfunction, full=.TRUE.)
      IF (PRESENT(input_variables)) THEN
         ALLOCATE (my_var(SIZE(input_variables)))
         my_var = input_variables
      ELSE
         CALL section_vals_val_get(gen_section, name_kw_variables, i_rep_section=irep, c_vals=my_var)
      END IF
      IF (ASSOCIATED(my_var)) THEN
         mydim = SIZE(my_var)
      END IF
      ! Check for presence of Parameters
      CALL section_vals_val_get(gen_section, name_kw_parameters, i_rep_section=irep, n_rep_val=n_par)
      CALL section_vals_val_get(gen_section, name_kw_values, i_rep_section=irep, n_rep_val=n_val)
      check = (n_par > 0) .EQV. (n_val > 0)
      CPASSERT(check)
      CALL section_vals_val_get(gen_section, name_kw_units, i_rep_section=irep, n_rep_val=n_units)
      IF (n_par > 0) THEN
         ! Parameters
         ALLOCATE (my_par(0))
         ALLOCATE (my_val(0))
         DO i = 1, n_par
            isize = SIZE(my_par)
            CALL section_vals_val_get(gen_section, name_kw_parameters, i_rep_section=irep, i_rep_val=i, c_vals=my_par_tmp)
            nblank = COUNT(my_par_tmp == "")
            CALL reallocate(my_par, 1, isize + SIZE(my_par_tmp) - nblank)
            ind = 0
            DO j = 1, SIZE(my_par_tmp)
               IF (my_par_tmp(j) == "") CYCLE
               ind = ind + 1
               my_par(isize + ind) = my_par_tmp(j)
            END DO
         END DO
         DO i = 1, n_val
            isize = SIZE(my_val)
            CALL section_vals_val_get(gen_section, name_kw_values, i_rep_section=irep, i_rep_val=i, r_vals=my_val_tmp)
            CALL reallocate(my_val, 1, isize + SIZE(my_val_tmp))
            my_val(isize + 1:isize + SIZE(my_val_tmp)) = my_val_tmp
         END DO
         CPASSERT(SIZE(my_par) == SIZE(my_val))
         ! Optionally read the units for each parameter value
         ALLOCATE (my_units(0))
         IF (n_units > 0) THEN
            DO i = 1, n_units
               isize = SIZE(my_units)
               CALL section_vals_val_get(gen_section, name_kw_units, i_rep_section=irep, i_rep_val=i, c_vals=my_units_tmp)
               nblank = COUNT(my_units_tmp == "")
               CALL reallocate(my_units, 1, isize + SIZE(my_units_tmp) - nblank)
               ind = 0
               DO j = 1, SIZE(my_units_tmp)
                  IF (my_units_tmp(j) == "") CYCLE
                  ind = ind + 1
                  my_units(isize + ind) = my_units_tmp(j)
               END DO
            END DO
            CPASSERT(SIZE(my_units) == SIZE(my_val))
         END IF
         mydim = mydim + SIZE(my_val)
         IF (SIZE(my_val) == 0) THEN
            DEALLOCATE (my_par)
            DEALLOCATE (my_val)
            DEALLOCATE (my_units)
         END IF
      END IF
      ! Handle trivial case of a null function defined
      ALLOCATE (parameters(mydim))
      ALLOCATE (values(mydim))
      IF (mydim > 0) THEN
         parameters(1:SIZE(my_var)) = my_var
         values(1:SIZE(my_var)) = 0.0_dp
         IF (PRESENT(var_values)) THEN
            CPASSERT(SIZE(var_values) == SIZE(my_var))
            values(1:SIZE(my_var)) = var_values
         END IF
         IF (ASSOCIATED(my_val)) THEN
            DO i = 1, SIZE(my_val)
               parameters(SIZE(my_var) + i) = my_par(i)
               IF (n_units > 0) THEN
                  values(SIZE(my_var) + i) = cp_unit_to_cp2k(my_val(i), TRIM(ADJUSTL(my_units(i))))
               ELSE
                  values(SIZE(my_var) + i) = my_val(i)
               END IF
            END DO
         END IF
      END IF
      IF (ASSOCIATED(my_par)) THEN
         DEALLOCATE (my_par)
      END IF
      IF (ASSOCIATED(my_val)) THEN
         DEALLOCATE (my_val)
      END IF
      IF (ASSOCIATED(my_units)) THEN
         DEALLOCATE (my_units)
      END IF
      IF (PRESENT(input_variables)) THEN
         DEALLOCATE (my_var)
      END IF
   END SUBROUTINE parse_function_section

END MODULE input_cp2k_function
