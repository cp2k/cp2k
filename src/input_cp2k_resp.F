!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief function that builds the resp section of the input
!> \par History
!>      02.2007 created
!> \author Joost VandeVondele
! *****************************************************************************
MODULE input_cp2k_resp
  USE cp_output_handling,              ONLY: cp_print_key_section_create
  USE cp_units,                        ONLY: cp_unit_to_cp2k
  USE f77_blas
  USE input_constants,                 ONLY: &
       add_last_numeric, do_resp_minus_x_dir, do_resp_minus_y_dir, &
       do_resp_minus_z_dir, do_resp_x_dir, do_resp_y_dir, do_resp_z_dir, &
       gaussian, high_print_level, low_print_level
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
  USE string_utilities,                ONLY: s2a
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_resp'

  PUBLIC :: create_resp_section

CONTAINS

! *****************************************************************************
!> \brief Creates the RESP section
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE create_resp_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_resp_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: subsection

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="RESP",&
            description="Requests a RESP fit of charges. When using a periodic "//&
            "Poisson solver and a periodic cell, the periodic RESP routines are "//&
            "used. If the Hartree potential matches with the one of an isolated "//&
            "system (i.e. isolated Poisson solver and big, nonperiodic cells), "//&
            "the nonperiodic RESP routines are automatically used. The subsections "//&
            "NONPERIODIC_SYS and PERIODIC_SYS do not request a nonperiodic or periodic "//&
            "solution, they only determine the sampling of the fitting points. All "//&
            "restraints are harmonic!",&
            n_keywords=2, n_subsections=2, repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword, subsection)

       CALL keyword_create(keyword, name="stride",&
            description="The stride (X,Y,Z) used to write the cube file "//&
            "(larger values result in smaller cube files). You can provide "//&
            "3 numbers (for X,Y,Z) or 1 number valid for all components.",&
            usage="STRIDE 2 2 2",n_var=-1,default_i_vals=(/2,2,2/), type_of_var=integer_t,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="INTEGER_TOTAL_CHARGE",&
            description="Forces the total charge to be integer",&
            usage="INTEGER_TOTAL_CHARGE TRUE",&
            default_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RESTRAIN_HEAVIES_TO_ZERO",&
            description="Restrain non-hydrogen atoms to zero.",&
            usage="RESTRAIN_HEAVIES_TO_ZERO FALSE",&
            default_l_val=.TRUE., error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RESTRAIN_HEAVIES_STRENGTH",&
            description="If defined, enforce the restraint of non-hydrogen "//&
            "atoms to zero. Its value is the strength of the restraint on "//&
            "the heavy atoms.",&
            usage="RESTRAIN_HEAVIES_STRENGTH 0.0001 ",&
            default_r_val=1.0E-6_dp ,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="WIDTH",&
            description="Specifies the value of the width of the Gaussian "//&
                        "charge distribution carried by each atom. Needs only "//&
                        "to be specified when using a periodic Poisson solver.",&
            usage="WIDTH <real> ",n_var=1,type_of_var=real_t,&
            default_r_val=cp_unit_to_cp2k(value=11.249_dp,unit_str="angstrom^-2",error=error), &
            unit_str="angstrom^-2",error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL create_constraint_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_restraint_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_nonperiodic_sys_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_periodic_sys_section(subsection,error)
       CALL section_add_subsection(section, subsection, error=error)
       CALL section_release(subsection,error=error)

       CALL create_print_resp_section(subsection, error=error)
       CALL section_add_subsection(section,subsection,error=error)
       CALL section_release(subsection,error=error)

    END IF
  END SUBROUTINE create_resp_section

! *****************************************************************************
!> \brief specifies constraints to be satisfied in a resp fit
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE create_constraint_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_constraint_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="CONSTRAINT",&
            description="specifies a linear constraint on the fitted charges."//&
                        "This can be used to give equal values to equivalent atoms." //&
                        "sum over atom_list c_i * q_i = t",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="TARGET",&
            description="the target value for the constraint",&
            usage="TARGET 0.0",&
            n_var=1,default_r_val=0.0_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="EQUAL_CHARGES",&
            description="All atoms in ATOM_LIST are constrained to have the "//&
            "same charges. When using this keyword, TARGET and ATOM_COEF do "//& 
            "not need to be set and will be ignored. Instead of using this "//&
            "keyword, the constraint section could be repeated.",&
            usage="EQUAL_CHARGES",&
            default_l_val=.FALSE.,lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="ATOM_LIST",&
            description="Defines the list of atoms involved in this constraint",&
            usage="ATOM_LIST 3 4",&
            type_of_var=integer_t, n_var=-1,required=.TRUE.,repeats=.TRUE.,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="ATOM_COEF",&
            description="Defines the coefficient of the atom in this "//&
                        "linear constraint",&
            usage="ATOM_COEF 1.0 -1.0",&
            type_of_var=real_t, n_var=-1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_constraint_section

! *****************************************************************************
!> \brief specifies restraints to be added to a resp fit
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Joost VandeVondele
! *****************************************************************************
  SUBROUTINE create_restraint_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_restraint_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="RESTRAINT",&
            description="specifies a restraint on the fitted charges."//&
                        "This can be used to restrain values to zero." //&
                        "s*(sum over atom_list q_i - t)**2",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="TARGET",&
            description="the target value for the restraint",&
            usage="TARGET 0.0",&
            n_var=1,default_r_val=0.0_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="STRENGTH",&
            description="the target value for the constraint",&
            usage="STRENGTH 0.001",&
            n_var=1,default_r_val=0.001_dp,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="ATOM_LIST",&
            description="Defines the list of atoms involved in this restraint",&
            usage="ATOM_LIST 3 4",&
            type_of_var=integer_t, n_var=-1, required=.TRUE.,repeats=.TRUE.,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword,name="ATOM_COEF",&
            description="Defines the coefficient of the atom in this "//&
                        "linear restraint. If given, the restraint will be: "//&
                        "s*(sum over atom_list c_i * q_i - t)**2 ",&
            usage="ATOM_COEF 1.0 -1.0",&
            type_of_var=real_t, n_var=-1, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_restraint_section

! *****************************************************************************
!> \brief specifies the parameter for sampling the resp fitting points for 
!>        non-periodic systems 
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Dorothea Golze
! *****************************************************************************
  SUBROUTINE create_nonperiodic_sys_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'create_nonperiodic_sys_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="NONPERIODIC_SYS",&
            description="Specifies the parameter for sampling the RESP fitting points "//&
                        "for non-periodic systems, i.e. systems that do not involve " //&
                        "surfaces. This section can be used with periodic and nonperiodic "//&
                        "Poisson solvers, it only affects the sampling of grid points. "//&
                        "All grid points in the shell defined by rmin and rmax are accepted "//&
                        "for fitting.",&
            n_keywords=1, n_subsections=0, repeats=.FALSE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="X_LOW",&
            description="Specifies the lower boundary of the box along X used to sample the potential.",&
            usage="X_LOW  -15.", type_of_var=real_t, n_var=1, unit_str='angstrom', error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="X_HI",&
            description="Specifies the upper boundary of the box along X used to sample the potential.",&
            usage="X_HI  5.", type_of_var=real_t, n_var=1, unit_str='angstrom', error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="Y_LOW",&
            description="Specifies the lower boundary of the box along Y used to sample the potential.",&
            usage="Y_LOW  -15.", type_of_var=real_t, n_var=1, unit_str='angstrom', error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="Y_HI",&
            description="Specifies the upper boundary of the box along Y used to sample the potential.",&
            usage="Y_HI  5.", type_of_var=real_t, n_var=1, unit_str='angstrom', error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="Z_LOW",&
            description="Specifies the lower boundary of the box along Z used to sample the potential.",&
            usage="Z_LOW  -15.", type_of_var=real_t, n_var=1, unit_str='angstrom', error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="Z_HI",&
            description="Specifies the upper boundary of the box along Z used to sample the potential.",&
            usage="Z_HI  5.", type_of_var=real_t, n_var=1, unit_str='angstrom', error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RMAX",&
            description="Specifies the maximum distance a fit point is away from an atom. "//&
                        "Valid for all atomic kinds for which no RMAX_KIND are specified.",&
            usage="RMAX 2.5",&
            default_r_val=cp_unit_to_cp2k(value=2.5_dp, unit_str="angstrom", error=error),&
            unit_str='angstrom', error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RMIN",&
            description="Specifies the minimum distance a fit point is away from an atom. "//&
                        "Valid for all atomic kinds for which no RMIN_KIND are specified.",&
            usage="RMIN 2.1",&
            default_r_val=cp_unit_to_cp2k(value=2.1_dp, unit_str="angstrom", error=error),&
            unit_str='angstrom', error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RMAX_KIND",&
            description="Specifies the maximum distance a fit point is away from an atom "//&
                        "of a given kind",&
            usage="RMAX 2.5 Br", required=.FALSE., repeats=.TRUE.,&
            n_var=-1,type_of_var=char_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RMIN_KIND",&
            description="Specifies the minimum distance a fit point is away from an atom "//&
                        "of a given kind",&
            usage="RMIN 2.1 Br", required=.FALSE., repeats=.TRUE.,&
            n_var=-1,type_of_var=char_t, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_nonperiodic_sys_section
! *****************************************************************************
!> \brief specifies the parameter for sampling the resp fitting points for 
!>        periodic systems, i.e. systems that involve surfaces 
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Dorothea Golze
! *****************************************************************************
  SUBROUTINE create_periodic_sys_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_periodic_sys_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="PERIODIC_SYS",&
            description="Specifies the parameter for sampling the RESP fitting "//&
                        "points for periodic systems, i.e. systems that involve " //&
                       "surfaces. This section can only be used with periodic  "//&
                        "Poisson solver and cell. To see, which grid points were "//&
                        "used, switch on COORD_FIT_POINTS in the PRINT section.",&
            n_keywords=1, n_subsections=0, repeats=.TRUE., required=.TRUE.,&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword,name="ATOM_LIST",&
            description="Specifies the list of indexes of atoms used to define "//&
                        "the region for the RESP fitting. The list should "//&
                        "contain indexes of atoms of the first surface layer." ,&
            usage="ATOM_LIST 1 2 3 or 1..3",&
            type_of_var=integer_t, n_var=-1,required=.TRUE.,repeats=.TRUE.,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="RANGE",&
            description="Range where the fitting points are sampled. A range of "//&
            "3 to 5 Angstroms means that the fitting points are sampled in the region "//&
            "of 3 to 5 Angstroms above the surface which is defined by atom indexes given "//&
            "in ATOM_LIST.",&
            usage="RANGE <real> <real>",unit_str="angstrom",n_var=2,type_of_var=real_t,&
            required=.TRUE.,error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="LENGTH",&
            description="Length of the sampling box, i.e. a box of this length and "//&
            "the height specified by RANGE is defined above each surface atom given "//&
            "in ATOM_LIST. The grid points in the boxes are accepted as fitting point. "//&
            "Should be in the range of the nearest neighbour distance (a bit larger to be "//&
            "on the safe side). Allows for a refined sampling of grid points in case of "//&
            "corrugated surfaces.",&
            usage="LENGTH <real> ",unit_str="angstrom",n_var=1,type_of_var=real_t,&
            default_r_val=cp_unit_to_cp2k(value=3.0_dp,unit_str="angstrom", error=error),error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       CALL keyword_create(keyword, name="SURF_DIRECTION",&
            description="Specifies what above the surface means. Defines the direction.",&
            usage="SURF_DIRECTION Z",&
            enum_c_vals=s2a( "X", "Y", "Z", "-X", "-Y", "-Z"),&
            enum_i_vals=(/ do_resp_x_dir, do_resp_y_dir, do_resp_z_dir,&
                           do_resp_minus_x_dir, do_resp_minus_y_dir, do_resp_minus_z_dir /),&
            enum_desc=s2a("surface layers are piled up in x-direction", &
                          "surface layers are piled up in y-direction",&
                          "surface layers are piled up in z-direction",&
                          "surface layers are piled up in -x-direction",&
                          "surface layers are piled up in -y-direction",&
                          "surface layers are piled up in -z-direction"),&
            default_i_val=do_resp_z_dir, error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

    END IF

  END SUBROUTINE create_periodic_sys_section

! *****************************************************************************
!> \brief create the resp print section 
!> \param section the section to create
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Dorothea Golze
! *****************************************************************************
  SUBROUTINE create_print_resp_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_print_resp_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key

    failure=.FALSE.
    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       NULLIFY(print_key, keyword)
       CALL section_create(section,name="print",&
            description="Section of possible print options specific for the RESP code.",&
            n_keywords=0, n_subsections=1, repeats=.FALSE., required=.TRUE.,&
            error=error)

       CALL cp_print_key_section_create(print_key,"PROGRAM_RUN_INFO",&
            description="Controls the printing of information regarding the run.",&
            print_level=low_print_level,filename="__STD_OUT__",&
            error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"COORD_FIT_POINTS",&
            description="Controls the printing of the coordinates of the "// &
            "grid points used for periodic RESP fitting. This section "//&
            "is intended to be only used for testing (you can get large files).",&
            print_level=high_print_level,add_last=add_last_numeric,&
            filename="RESP_FIT_POINTS",&
            common_iter_levels=3, error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"RESP_CHARGES_TO_FILE",&
            description="Controls the printing of the RESP charges "// &
            "to a file.",&
            print_level=high_print_level,add_last=add_last_numeric,&
            filename="RESP_CHARGES",&
            common_iter_levels=3, error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)

       CALL cp_print_key_section_create(print_key,"V_RESP_CUBE",&
            description="Controls the printing of the potential generated "//&
            "by the RESP CHARGES to a cube file. Prints the relative "//&
            "root-mean-square (RRMS) and root-mean-square (RMS) errors.",&
            print_level=high_print_level,add_last=add_last_numeric,&
            filename="RESP_POTENTIAL",&
            common_iter_levels=3, error=error)
       CALL keyword_create(keyword, name="stride",&
            description="The stride (X,Y,Z) used to write the cube file "//&
            "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"//&
            " 1 number valid for all components.",&
            usage="STRIDE 2 2 2",n_var=-1,default_i_vals=(/2,2,2/), type_of_var=integer_t,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL keyword_create(keyword, name="APPEND",&
            description="append the cube files when they already exist",&
            default_l_val=.FALSE., lone_keyword_l_val=.TRUE.,error=error)
       CALL section_add_keyword(print_key,keyword,error=error)
       CALL keyword_release(keyword,error=error)
       CALL section_add_subsection(section,print_key,error=error)
       CALL section_release(print_key,error=error)
   ENDIF
  END SUBROUTINE create_print_resp_section

END MODULE input_cp2k_resp
