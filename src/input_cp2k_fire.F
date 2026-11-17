!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Input shared by the geometry, cell and band FIRE optimizers.
! **************************************************************************************************
MODULE input_cp2k_fire
   USE bibliography,                    ONLY: Bitzek2006,&
                                              Guenole2020
   USE fire_minimizer,                  ONLY: fire_parameters_type
   USE input_keyword_types,             ONLY: keyword_create,&
                                              keyword_release,&
                                              keyword_type
   USE input_section_types,             ONLY: section_add_keyword,&
                                              section_create,&
                                              section_type,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: create_fire_section, read_fire_section

CONTAINS

! **************************************************************************************************
!> \brief Create the FIRE subsection.
!> \param section subsection to create
!> \param band_opt add a band iteration limit
!> \param cell_opt add a cell mass override
! **************************************************************************************************
   SUBROUTINE create_fire_section(section, band_opt, cell_opt)
      TYPE(section_type), POINTER                        :: section
      LOGICAL, INTENT(IN)                                :: band_opt, cell_opt

      TYPE(fire_parameters_type), PARAMETER :: defaults = fire_parameters_type()

      TYPE(keyword_type), POINTER                        :: keyword

      NULLIFY (keyword)
      CALL section_create(section, __LOCATION__, name="FIRE", &
                          description="Parameters shared by the FIRE and FIRE2 optimizers.", &
                          citations=[Bitzek2006, Guenole2020], repeats=.FALSE.)
      CALL keyword_create(keyword, __LOCATION__, name="TIMESTEP", unit_str="fs", &
                          description="Initial timestep of the fictitious dynamics, using unit atomic masses (1 amu).", &
                          default_r_val=defaults%timestep)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="MAX_TIMESTEP_FACTOR", &
                          description="Maximum timestep divided by the initial timestep.", &
                          default_r_val=defaults%max_timestep_factor)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="MAX_STEP", unit_str="angstrom", &
                          description="Maximum length of the forward step for any particle in optimizer coordinates. "// &
                          "In cell optimization the norm of the reduced cell step is also capped by this value.", &
                          default_r_val=defaults%max_step)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)
      IF (cell_opt) THEN
         CALL keyword_create(keyword, __LOCATION__, name="CELL_MASS_FACTOR", &
                             description="Fictitious cell mass relative to the atomic mass used by FIRE. "// &
                             "Zero selects the cube root of the number of particles.", default_r_val=0.0_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
      END IF
      IF (band_opt) THEN
         CALL keyword_create(keyword, __LOCATION__, name="MAX_STEPS", &
                             description="Maximum number of band optimization steps.", default_i_val=100)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
      END IF
   END SUBROUTINE create_fire_section

! **************************************************************************************************
!> \brief Read and validate options common to FIRE and FIRE2.
!> \param section values of the FIRE subsection
!> \param params integrator parameters
! **************************************************************************************************
   SUBROUTINE read_fire_section(section, params)
      TYPE(section_vals_type), POINTER                   :: section
      TYPE(fire_parameters_type), INTENT(OUT)            :: params

      CALL section_vals_val_get(section, "TIMESTEP", r_val=params%timestep)
      CALL section_vals_val_get(section, "MAX_TIMESTEP_FACTOR", r_val=params%max_timestep_factor)
      CALL section_vals_val_get(section, "MAX_STEP", r_val=params%max_step)
      IF (.NOT. (params%timestep > 0.0_dp .AND. params%max_step > 0.0_dp)) THEN
         CPABORT("FIRE requires positive TIMESTEP and MAX_STEP.")
      END IF
      IF (.NOT. (params%max_timestep_factor >= 1.0_dp)) THEN
         CPABORT("FIRE requires MAX_TIMESTEP_FACTOR >= 1.")
      END IF
   END SUBROUTINE read_fire_section

END MODULE input_cp2k_fire
