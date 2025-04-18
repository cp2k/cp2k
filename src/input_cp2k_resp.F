!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief function that builds the resp section of the input
!> \par History
!>      02.2007 created
!> \author Joost VandeVondele
! **************************************************************************************************
MODULE input_cp2k_resp
   USE bibliography,                    ONLY: Campana2009,&
                                              Golze2015,&
                                              Rappe1992
   USE cp_output_handling,              ONLY: add_last_numeric,&
                                              cp_print_key_section_create,&
                                              high_print_level,&
                                              low_print_level
   USE cp_units,                        ONLY: cp_unit_to_cp2k
   USE input_constants,                 ONLY: &
        do_resp_minus_x_dir, do_resp_minus_y_dir, do_resp_minus_z_dir, do_resp_x_dir, &
        do_resp_y_dir, do_resp_z_dir, gaussian, use_cambridge_vdw_radii, use_uff_vdw_radii
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
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_resp'

   PUBLIC :: create_resp_section

CONTAINS

! **************************************************************************************************
!> \brief Creates the RESP section
!> \param section the section to create
!> \author Joost VandeVondele
! **************************************************************************************************
   SUBROUTINE create_resp_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="RESP", &
                          description="Requests a RESP fit of charges. When using a periodic "// &
                          "Poisson solver and a periodic cell, the periodic RESP routines are "// &
                          "used. If the Hartree potential matches with the one of an isolated "// &
                          "system (i.e. isolated Poisson solver and big, nonperiodic cells), "// &
                          "the nonperiodic RESP routines are automatically used. All restraints "// &
                          "are harmonic!", &
                          n_keywords=2, n_subsections=2, repeats=.FALSE., citations=(/Golze2015/))

      NULLIFY (keyword, subsection)

      CALL keyword_create(keyword, __LOCATION__, name="stride", &
                          description="The stride (X,Y,Z) used to write the cube file "// &
                          "(larger values result in smaller cube files). You can provide "// &
                          "3 numbers (for X,Y,Z) or 1 number valid for all components.", &
                          usage="STRIDE 2 2 2", n_var=-1, default_i_vals=(/2, 2, 2/), type_of_var=integer_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="INTEGER_TOTAL_CHARGE", &
                          description="Forces the total charge to be integer", &
                          usage="INTEGER_TOTAL_CHARGE TRUE", &
                          default_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RESTRAIN_HEAVIES_TO_ZERO", &
                          description="Restrain non-hydrogen atoms to zero.", &
                          usage="RESTRAIN_HEAVIES_TO_ZERO FALSE", &
                          default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RESTRAIN_HEAVIES_STRENGTH", &
                          description="If defined, enforce the restraint of non-hydrogen "// &
                          "atoms to zero. Its value is the strength of the restraint on "// &
                          "the heavy atoms.", &
                          usage="RESTRAIN_HEAVIES_STRENGTH 0.0001 ", &
                          default_r_val=1.0E-6_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="WIDTH", &
                          description="Specifies the value of the width of the Gaussian "// &
                          "charge distribution carried by each atom. Needs only "// &
                          "to be specified when using a periodic Poisson solver.", &
                          usage="WIDTH <real> ", n_var=1, type_of_var=real_t, &
                          default_r_val=cp_unit_to_cp2k(value=11.249_dp, unit_str="angstrom^-2"), &
                          unit_str="angstrom^-2")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="USE_REPEAT_METHOD", &
                          description="Fits the variance of the potential, i.e. the deviation "// &
                          "from the mean value of the potential within the selected "// &
                          "range. The evaluation of the potentials is still treated "// &
                          "within the GPW approach as described in [Golze2015]. "// &
                          "When used in conjunction with INTEGER_TOTAL_CHARGE = T "// &
                          "and SPHERE_SAMPLING, the results will be very similar to "// &
                          "the REPEAT charges given in [Campana2009]. In most "// &
                          "cases switching on this option gives reasonable "// &
                          "atomic charges without the need to define any "// &
                          "restraints. Note that by switching on this option, "// &
                          "RESTRAIN_HEAVIES_TO_ZERO will be switched off. ", &
                          usage="USE_REPEAT_METHOD", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE., &
                          citations=(/Campana2009/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL create_constraint_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_restraint_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_sphere_sampling_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_slab_sampling_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_print_resp_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_resp_section

! **************************************************************************************************
!> \brief specifies constraints to be satisfied in a resp fit
!> \param section the section to create
!> \author Joost VandeVondele
! **************************************************************************************************
   SUBROUTINE create_constraint_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="CONSTRAINT", &
                          description="specifies a linear constraint on the fitted charges. "// &
                          "This can be used to give equal values to equivalent atoms. "// &
                          "sum over atom_list c_i * q_i = t", &
                          n_keywords=1, n_subsections=0, repeats=.TRUE.)

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="TARGET", &
                          description="the target value for the constraint", &
                          usage="TARGET 0.0", &
                          n_var=1, default_r_val=0.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EQUAL_CHARGES", &
                          description="All atoms in ATOM_LIST are constrained to have the "// &
                          "same charges. When using this keyword, TARGET and ATOM_COEF do "// &
                          "not need to be set and will be ignored. Instead of using this "// &
                          "keyword, the constraint section could be repeated.", &
                          usage="EQUAL_CHARGES", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ATOM_LIST", &
                          description="Defines the list of atoms involved in this constraint", &
                          usage="ATOM_LIST 3 4", &
                          type_of_var=integer_t, n_var=-1, repeats=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ATOM_COEF", &
                          description="Defines the coefficient of the atom in this "// &
                          "linear constraint", &
                          usage="ATOM_COEF 1.0 -1.0", &
                          type_of_var=real_t, n_var=-1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_constraint_section

! **************************************************************************************************
!> \brief specifies restraints to be added to a resp fit
!> \param section the section to create
!> \author Joost VandeVondele
! **************************************************************************************************
   SUBROUTINE create_restraint_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="RESTRAINT", &
                          description="specifies a restraint on the fitted charges. "// &
                          "This can be used to restrain values to zero. "// &
                          "s*(sum over atom_list q_i - t)**2", &
                          n_keywords=1, n_subsections=0, repeats=.TRUE.)

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="TARGET", &
                          description="the target value for the restraint", &
                          usage="TARGET 0.0", &
                          n_var=1, default_r_val=0.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="STRENGTH", &
                          description="the target value for the constraint", &
                          usage="STRENGTH 0.001", &
                          n_var=1, default_r_val=0.001_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ATOM_LIST", &
                          description="Defines the list of atoms involved in this restraint", &
                          usage="ATOM_LIST 3 4", &
                          type_of_var=integer_t, n_var=-1, repeats=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ATOM_COEF", &
                          description="Defines the coefficient of the atom in this "// &
                          "linear restraint. If given, the restraint will be: "// &
                          "s*(sum over atom_list c_i * q_i - t)**2 ", &
                          usage="ATOM_COEF 1.0 -1.0", &
                          type_of_var=real_t, n_var=-1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_restraint_section

! **************************************************************************************************
!> \brief specifies the parameter for sampling the resp fitting points for
!>        molecular structures; sampling in spheres around the atoms
!> \param section the section to create
!> \author Dorothea Golze
! **************************************************************************************************
   SUBROUTINE create_sphere_sampling_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="SPHERE_SAMPLING", &
                          description="Specifies the parameter for sampling the RESP fitting points "// &
                          "for molecular structures, i.e. systems that do not involve "// &
                          "surfaces. Fitting points are sampled in spheres around the "// &
                          "atom. All grid points in the shell defined by rmin and rmax "// &
                          "are accepted for fitting. Default is that rmin is the vdW "// &
                          "radius and rmax=100.0*vdW_radius, which can be overwritten "// &
                          "by the keywords below.", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="X_LOW", &
                          description="Specifies the lower boundary of the box along X used to "// &
                          "sample the potential. Only valid for nonperiodic RESP fitting.", &
                          usage="X_LOW  -15.", type_of_var=real_t, n_var=1, unit_str='angstrom')
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="X_HI", &
                          description="Specifies the upper boundary of the box along X used to "// &
                          "sample the potential. Only valid for nonperiodic RESP fitting.", &
                          usage="X_HI  5.", type_of_var=real_t, n_var=1, unit_str='angstrom')
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="Y_LOW", &
                          description="Specifies the lower boundary of the box along Y used to "// &
                          "sample the potential. Only valid for nonperiodic RESP fitting.", &
                          usage="Y_LOW  -15.", type_of_var=real_t, n_var=1, unit_str='angstrom')
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="Y_HI", &
                          description="Specifies the upper boundary of the box along Y used to "// &
                          "sample the potential. Only valid for nonperiodic RESP fitting.", &
                          usage="Y_HI  5.", type_of_var=real_t, n_var=1, unit_str='angstrom')
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="Z_LOW", &
                          description="Specifies the lower boundary of the box along Z used to "// &
                          "sample the potential. Only valid for nonperiodic RESP fitting.", &
                          usage="Z_LOW  -15.", type_of_var=real_t, n_var=1, unit_str='angstrom')
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="Z_HI", &
                          description="Specifies the upper boundary of the box along Z used to "// &
                          "sample the potential. Only valid for nonperiodic RESP fitting.", &
                          usage="Z_HI  5.", type_of_var=real_t, n_var=1, unit_str='angstrom')
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="AUTO_VDW_RADII_TABLE", &
                          description="Select which vdW radii table to use for automatic "// &
                          "determination of RMIN_KIND and RMAX_KIND if those "// &
                          "are not declared explicitly", &
                          usage="AUTO_VDW_RADII_TABLE UFF", &
                          default_i_val=use_cambridge_vdw_radii, &
                          enum_c_vals=s2a("CAMBRIDGE", &
                                          "UFF"), &
                          enum_desc=s2a("Cambridge Structural Database", &
                                        "Universal Force Field: "// &
                                        "Rappe et al. J. Am. Chem. Soc. 114, 10024 (1992)"), &
                          enum_i_vals=(/use_cambridge_vdw_radii, &
                                        use_uff_vdw_radii/), &
                          citations=(/Rappe1992/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="AUTO_RMAX_SCALE", &
                          description="IF RMAX or RMAX_KIND keywords are not present, defines the "// &
                          "maximumn distance a fit point is away from an atom based on "// &
                          "the formula: rmax(kind) = AUTO_RMAX_SCALE * vdW_radius(kind). "// &
                          "The van der Waals radiii of the elements are based on data from "// &
                          "table chosen by AUTO_VDW_RADII_TABLE.", &
                          usage="AUTO_RMAX_SCALE 60.0", &
                          default_r_val=100.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="AUTO_RMIN_SCALE", &
                          description="IF RMIN or RMIN_KIND keywords are not present, defines the "// &
                          "minimum distance a fit point is away from an atom based on "// &
                          "the formula: rmin(kind) = AUTO_RMIN_SCALE * vdW_radius(kind). "// &
                          "The van der Waals radii of the elements are based on data from "// &
                          "table chosen by AUTO_VDW_RADII_TABLE.", &
                          usage="AUTO_RMIN_SCALE 1.5", &
                          default_r_val=1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RMAX", &
                          description="Specifies the maximum distance a fit point is away from an atom. "// &
                          "Valid for all atomic kinds for which no RMAX_KIND are specified.", &
                          usage="RMAX 2.5", &
                          default_r_val=cp_unit_to_cp2k(value=2.5_dp, unit_str="angstrom"), &
                          unit_str='angstrom')
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RMIN", &
                          description="Specifies the minimum distance a fit point is away from an atom. "// &
                          "Valid for all atomic kinds for which no RMIN_KIND are specified.", &
                          usage="RMIN 2.1", &
                          default_r_val=cp_unit_to_cp2k(value=2.1_dp, unit_str="angstrom"), &
                          unit_str='angstrom')
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RMAX_KIND", &
                          description="Specifies the maximum distance a fit point is away from an atom "// &
                          "of a given kind", &
                          usage="RMAX_KIND 2.5 Br", repeats=.TRUE., &
                          n_var=-1, type_of_var=char_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RMIN_KIND", &
                          description="Specifies the minimum distance a fit point is away from an atom "// &
                          "of a given kind", &
                          usage="RMIN_KIND 2.1 Br", repeats=.TRUE., &
                          n_var=-1, type_of_var=char_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_sphere_sampling_section

! **************************************************************************************************
!> \brief specifies the parameter for sampling the resp fitting points for
!>        slab-like periodic systems, i.e. systems that involve surfaces
!> \param section the section to create
!> \author Dorothea Golze
! **************************************************************************************************
   SUBROUTINE create_slab_sampling_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="SLAB_SAMPLING", &
                          description="Specifies the parameter for sampling the RESP fitting "// &
                          "points for slab-like periodic systems, i.e. systems that "// &
                          "involve surfaces. This section can only be used with periodic  "// &
                          "Poisson solver and cell. To see, which grid points were "// &
                          "used, switch on COORD_FIT_POINTS in the PRINT section.", &
                          n_keywords=1, n_subsections=0, repeats=.TRUE.)

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ATOM_LIST", &
                          description="Specifies the list of indexes of atoms used to define "// &
                          "the region for the RESP fitting. The list should "// &
                          "contain indexes of atoms of the first surface layer.", &
                          usage="ATOM_LIST 1 2 3 or 1..3", &
                          type_of_var=integer_t, n_var=-1, repeats=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="RANGE", &
                          description="Range where the fitting points are sampled. A range of "// &
                          "3 to 5 Angstroms means that the fitting points are sampled in the region "// &
                          "of 3 to 5 Angstroms above the surface which is defined by atom indexes given "// &
                          "in ATOM_LIST.", &
                          usage="RANGE <real> <real>", unit_str="angstrom", n_var=2, type_of_var=real_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="LENGTH", &
                          description="Length of the sampling box, i.e. a box of this length and "// &
                          "the height specified by RANGE is defined above each surface atom given "// &
                          "in ATOM_LIST. The grid points in the boxes are accepted as fitting point. "// &
                          "Should be in the range of the nearest neighbour distance (a bit larger to be "// &
                          "on the safe side). Allows for a refined sampling of grid points in case of "// &
                          "corrugated surfaces.", &
                          usage="LENGTH <real> ", unit_str="angstrom", n_var=1, type_of_var=real_t, &
                          default_r_val=cp_unit_to_cp2k(value=3.0_dp, unit_str="angstrom"))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SURF_DIRECTION", &
                          description="Specifies what above the surface means. Defines the direction.", &
                          usage="SURF_DIRECTION Z", &
                          enum_c_vals=s2a("X", "Y", "Z", "-X", "-Y", "-Z"), &
                          enum_i_vals=(/do_resp_x_dir, do_resp_y_dir, do_resp_z_dir, &
                                        do_resp_minus_x_dir, do_resp_minus_y_dir, do_resp_minus_z_dir/), &
                          enum_desc=s2a("surface layers are piled up in x-direction", &
                                        "surface layers are piled up in y-direction", &
                                        "surface layers are piled up in z-direction", &
                                        "surface layers are piled up in -x-direction", &
                                        "surface layers are piled up in -y-direction", &
                                        "surface layers are piled up in -z-direction"), &
                          default_i_val=do_resp_z_dir)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_slab_sampling_section

! **************************************************************************************************
!> \brief create the resp print section
!> \param section the section to create
!> \author Dorothea Golze
! **************************************************************************************************
   SUBROUTINE create_print_resp_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: print_key

      CPASSERT(.NOT. ASSOCIATED(section))
      NULLIFY (print_key, keyword)
      CALL section_create(section, __LOCATION__, name="print", &
                          description="Section of possible print options specific for the RESP code.", &
                          n_keywords=0, n_subsections=1, repeats=.FALSE.)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "PROGRAM_RUN_INFO", &
                                       description="Controls the printing of information regarding the run.", &
                                       print_level=low_print_level, filename="__STD_OUT__")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "COORD_FIT_POINTS", &
                                       description="Controls the printing of the coordinates of the "// &
                                       "grid points used for periodic RESP fitting. This section "// &
                                       "is intended to be only used for testing (you can get large files).", &
                                       print_level=high_print_level, add_last=add_last_numeric, &
                                       filename="RESP_FIT_POINTS", &
                                       common_iter_levels=3)
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "RESP_CHARGES_TO_FILE", &
                                       description="Controls the printing of the RESP charges "// &
                                       "to a file.", &
                                       print_level=high_print_level, add_last=add_last_numeric, &
                                       filename="RESP_CHARGES", &
                                       common_iter_levels=3)
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "V_RESP_CUBE", &
                                       description="Controls the printing of the potential generated "// &
                                       "by the RESP CHARGES to a cube file. Prints the relative "// &
                                       "root-mean-square (RRMS) and root-mean-square (RMS) errors.", &
                                       print_level=high_print_level, add_last=add_last_numeric, &
                                       filename="RESP_POTENTIAL", &
                                       common_iter_levels=3)
      CALL keyword_create(keyword, __LOCATION__, name="stride", &
                          description="The stride (X,Y,Z) used to write the cube file "// &
                          "(larger values result in smaller cube files). You can provide 3 numbers (for X,Y,Z) or"// &
                          " 1 number valid for all components.", &
                          usage="STRIDE 2 2 2", n_var=-1, default_i_vals=(/2, 2, 2/), type_of_var=integer_t)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="APPEND", &
                          description="append the cube files when they already exist", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(print_key, keyword)
      CALL keyword_release(keyword)
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)
   END SUBROUTINE create_print_resp_section

END MODULE input_cp2k_resp
