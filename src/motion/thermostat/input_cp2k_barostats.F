!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \par History
!>      10.2005 split input_cp2k into smaller modules [fawzi]
!> \author teo & fawzi
! **************************************************************************************************
MODULE input_cp2k_barostats
   USE barostat_types,                  ONLY: do_clv_x,&
                                              do_clv_xy,&
                                              do_clv_xyz,&
                                              do_clv_xz,&
                                              do_clv_y,&
                                              do_clv_yz,&
                                              do_clv_z
   USE cp_output_handling,              ONLY: cp_print_key_section_create,&
                                              high_print_level
   USE cp_units,                        ONLY: cp_unit_to_cp2k
   USE input_cp2k_thermostats,          ONLY: create_mass_section,&
                                              create_thermostat_section,&
                                              create_velocity_section
   USE input_keyword_types,             ONLY: keyword_create,&
                                              keyword_release,&
                                              keyword_type
   USE input_section_types,             ONLY: section_add_keyword,&
                                              section_add_subsection,&
                                              section_create,&
                                              section_release,&
                                              section_type
   USE input_val_types,                 ONLY: real_t
   USE kinds,                           ONLY: dp
   USE string_utilities,                ONLY: s2a
#include "../../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_barostats'

   PUBLIC :: create_barostat_section

!***
CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param section will contain the coeff section
!> \author teo
! **************************************************************************************************
   SUBROUTINE create_barostat_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection, thermo_section

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="barostat", &
                          description="Parameters of barostat.", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword, subsection, thermo_section)
      CALL keyword_create(keyword, __LOCATION__, name="PRESSURE", &
                          description="Initial pressure", &
                          usage="PRESSURE real", &
                          default_r_val=0._dp, unit_str='bar')
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="TIMECON", &
                          description="Barostat time constant", &
                          usage="TIMECON real", &
                          default_r_val=cp_unit_to_cp2k(1000.0_dp, "fs"), &
                          unit_str='fs')
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="TEMPERATURE", &
                          description="Barostat initial temperature. If not set, the ensemble temperature is used instead.", &
                          usage="TEMPERATURE real", type_of_var=real_t, &
                          unit_str='K')
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="TEMP_TOL", &
                          description="Maximum oscillation of the Barostat temperature imposed by rescaling.", &
                          usage="TEMP_TOL real", default_r_val=0._dp, &
                          unit_str='K')
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="VIRIAL", &
                          description="For NPT_F only: allows the screening of one or more components of the virial in order"// &
                          " to relax the cell only along specific cartesian axis", &
                          usage="VIRIAL (XYZ | X | Y | Z | XY| XZ | YZ)", &
                          enum_c_vals=s2a("XYZ", "X", "Y", "Z", "XY", "XZ", "YZ"), &
                          enum_i_vals=(/do_clv_xyz, do_clv_x, do_clv_y, do_clv_z, do_clv_xy, do_clv_xz, do_clv_yz/), &
                          default_i_val=do_clv_xyz)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL create_velocity_section(subsection, "BAROSTAT")
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_mass_section(subsection, "BAROSTAT")
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_thermostat_section(thermo_section, coupled_thermostat=.TRUE.)
      CALL section_add_subsection(section, thermo_section)
      CALL section_release(thermo_section)

      CALL create_print_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_barostat_section

! **************************************************************************************************
!> \brief Creates print section for barostat section
!> \param section ...
!> \author teo [tlaino] - University of Zurich - 02.2008
! **************************************************************************************************
   SUBROUTINE create_print_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(section_type), POINTER                        :: print_key

      CPASSERT(.NOT. ASSOCIATED(section))
      NULLIFY (print_key)
      CALL section_create(section, __LOCATION__, name="PRINT", &
                          description="Collects all print_keys for barostat", &
                          n_keywords=1, n_subsections=0, repeats=.FALSE.)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "ENERGY", &
                                       description="Controls the output of kinetic energy, and potential energy "// &
                                       "of the defined barostat.", print_level=high_print_level, common_iter_levels=1, &
                                       filename="")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)
   END SUBROUTINE create_print_section

END MODULE input_cp2k_barostats
