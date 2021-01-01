!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief function that build the input sections for external [potential, density VXC]
!> \par History
!>      10.2005 moved out of input_cp2k [fawzi]
!>      10.2020 moved out of input_cp2k_dft [JGH]
!> \author fawzi
! **************************************************************************************************
MODULE input_cp2k_external
   USE bibliography,                    ONLY: Tozer1996,&
                                              Zhao1994
   USE input_constants,                 ONLY: use_coulomb,&
                                              use_diff,&
                                              use_no
   USE input_keyword_types,             ONLY: keyword_create,&
                                              keyword_release,&
                                              keyword_type
   USE input_section_types,             ONLY: section_add_keyword,&
                                              section_add_subsection,&
                                              section_create,&
                                              section_release,&
                                              section_type
   USE input_val_types,                 ONLY: char_t,&
                                              lchar_t,&
                                              real_t
   USE kinds,                           ONLY: dp
   USE string_utilities,                ONLY: s2a
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_external'

   PUBLIC :: create_ext_pot_section, create_ext_den_section, create_ext_vxc_section

CONTAINS

! **************************************************************************************************
!> \brief Creates the section for applying an electrostatic external potential
!> \param section ...
!> \date 12.2009
!> \author teo
! **************************************************************************************************
   SUBROUTINE create_ext_pot_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="EXTERNAL_POTENTIAL", &
                          description="Section controlling the presence of an electrostatic "// &
                          "external potential dependent on the atomic positions (X,Y,Z). "// &
                          "As the external potential is currently applied via a grid, "// &
                          "it only works with DFT based methods (GPW/GAPW) that already use "// &
                          "a grid based approach to solve the Poisson equation.", &
                          n_keywords=7, n_subsections=0, repeats=.FALSE.)
      NULLIFY (keyword, subsection)

      CALL keyword_create(keyword, __LOCATION__, name="FUNCTION", &
                          description="Specifies the functional form in mathematical notation. Variables must be the atomic "// &
                          "coordinates (X,Y,Z) of the grid.", usage="FUNCTION  X^2+Y^2+Z^2+LOG(ABS(X+Y))", &
                          type_of_var=lchar_t, n_var=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PARAMETERS", &
                          description="Defines the parameters of the functional form", &
                          usage="PARAMETERS a b D", type_of_var=char_t, &
                          n_var=-1, repeats=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="VALUES", &
                          description="Defines the values of parameter of the functional form", &
                          usage="VALUES ", type_of_var=real_t, &
                          n_var=-1, repeats=.TRUE., unit_str="internal_cp2k")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="UNITS", &
                          description="Optionally, allows to define valid CP2K unit strings for each parameter value. "// &
                          "It is assumed that the corresponding parameter value is specified in this unit.", &
                          usage="UNITS angstrom eV*angstrom^-1 angstrom^1 K", type_of_var=char_t, &
                          n_var=-1, repeats=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="STATIC", &
                          description="Specifies the external potential as STATIC or time dependent. At the moment "// &
                          "only static potentials are implemented.", &
                          usage="STATIC T", default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DX", &
                          description="Parameter used for computing the derivative with the Ridders' method.", &
                          usage="DX <REAL>", default_r_val=0.1_dp, unit_str="bohr")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ERROR_LIMIT", &
                          description="Checks that the error in computing the derivative is not larger than "// &
                          "the value set; in case error is larger a warning message is printed.", &
                          usage="ERROR_LIMIT <REAL>", default_r_val=1.0E-12_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      !keyword for reading the external potential from cube file
      CALL keyword_create(keyword, __LOCATION__, name="READ_FROM_CUBE", &
                          description="Switch for reading the external potential from file pot.cube. The values "// &
                          "of the potential must be on the grid points of the realspace grid.", &
                          usage="READ_FROM_CUBE T", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      !keyword for scaling the external potential that is read from file by a constant factor
      CALL keyword_create(keyword, __LOCATION__, name="SCALING_FACTOR", &
                          description="A factor for scaling the the external potential that is read from file."// &
                          "The value of the potential at each grid point is multiplied by this factor.", &
                          usage="SCALING_FACTOR <REAL>", default_r_val=1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL create_maxwell_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_ext_pot_section

! **************************************************************************************************
!> \brief Creates the section for applying an electrostatic external potential
!> \param section ...
!> \date 12.2009
!> \author teo
! **************************************************************************************************
   SUBROUTINE create_maxwell_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="MAXWELL", &
                          description="Section controlling the calculation of an electrostatic "// &
                          "external potential calculated from Maxwell equations. ", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)
      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="TEST_LOGICAL", &
                          description="Test for logical value", &
                          usage="TEST_LOGICAL T", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="TEST_REAL", &
                          description="TEST for Real", &
                          usage="TEST_REAL <REAL>", default_r_val=1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="TEST_INTEGER", &
                          description="TEST for Int", &
                          usage="TEST_INTEGER <INT>", default_i_val=0)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_maxwell_section

! **************************************************************************************************
!> \brief ZMP Creates the section for reading user supplied external density
!> \param section ...
!> \date 03.2011
!> \author D. Varsano [daniele.varsano@nano.cnr.it]
! **************************************************************************************************
   SUBROUTINE create_ext_den_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="EXTERNAL_DENSITY", &
                          description="Section for the use of the ZMP technique on external densities.", &
                          n_keywords=4, n_subsections=0, repeats=.FALSE., &
                          citations=(/Zhao1994, Tozer1996/))
      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FILE_DENSITY", &
                          description="Specifies the filename containing the target density in *.cube format."// &
                          "In the MGRID section it must be imposed NGRID 1, as it works with only"// &
                          "one grid. The number of points in each direction, and the spacing must"// &
                          "be previously defined choosing the plane waves cut-off in section MGRID"// &
                          "keyword CUTOFF, and the cube dimension in section SUBSYS / CELL / keyword ABC", &
                          usage="DENSITY_FILE_NAME <FILENAME>", &
                          type_of_var=char_t, default_c_val="RHO_O.dat", n_var=-1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="LAMBDA", &
                          description="Lagrange multiplier defined in the constraint ZMP method. When starting, use"// &
                          "small values when starting from scratch (around 5,10). Then gradually increase"// &
                          "the values depending, restarting from the previous calculation with the smaller"// &
                          "value. To choose the progressive values of LAMBDA look at the convergence of the"// &
                          " eigenvalues.", &
                          usage="DX <REAL>", default_r_val=10.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ZMP_CONSTRAINT", &
                          description="Specify which kind of constraint to solve the ZMP equation. The COULOMB default"// &
                          "option is more stable.", &
                          usage="ZMP_CONSTRAINT <CHAR>", &
                          enum_c_vals=s2a("COULOMB", "DIFF", "NONE"), &
                          enum_i_vals=(/use_coulomb, use_diff, use_no/), &
                          enum_desc=s2a("Coulomb constraint, integral of [rho_0(r)-rho(r)]/|r-r'|", &
                                        "Simple constraint, [rho_0(r)-rho(r)]", &
                                        "No constrain imposed"), &
                          default_i_val=use_coulomb)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FERMI_AMALDI", &
                          description="Add the Fermi-Amaldi contribution to the Hartree potential."// &
                          "It leads to a more stable convergence.", &
                          usage="FERMI_AMALDI <LOGICAL>", &
                          repeats=.FALSE., &
                          n_var=1, &
                          default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_ext_den_section

! **************************************************************************************************
!> \brief ZMP Creates the section for creating the external v_xc
!> \param section ...
!> \date 03.2011
!> \author D. Varsano [daniele.varsano@nano.cnr.it]
! **************************************************************************************************
   SUBROUTINE create_ext_vxc_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="EXTERNAL_VXC", &
                          description="SCF convergence with external v_xc calculated through previous ZMP"// &
                          "calculation", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)
      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="FILE_VXC", &
                          description="The *.cube filename containing the v_xc potential. This works only"// &
                          "with NGRID 1 imposed in the MGRID section. The number of points in each"// &
                          "direction, and the spacing must equal to those previously used in the ZMP"// &
                          "calculation and defined through the plane wave cut-off and the cube dimension"// &
                          "respectively set in section MGRID / keyword CUTOFF, and in section SUBSYS /"// &
                          "CELL / keyword ABC", &
                          usage="FILE_VXC <FILENAME>", &
                          type_of_var=char_t, default_c_val="VXC_O.dat", n_var=-1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)
   END SUBROUTINE create_ext_vxc_section

END MODULE input_cp2k_external
