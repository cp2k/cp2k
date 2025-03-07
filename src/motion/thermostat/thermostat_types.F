!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Thermostat structure: module containing thermostat available for MD
!> \author teo [tlaino] - University of Zurich - 09.2007
! **************************************************************************************************
MODULE thermostat_types
   USE al_system_types,                 ONLY: al_dealloc,&
                                              al_init,&
                                              al_system_type
   USE csvr_system_types,               ONLY: csvr_dealloc,&
                                              csvr_init,&
                                              csvr_system_type
   USE extended_system_types,           ONLY: lnhc_dealloc,&
                                              lnhc_init,&
                                              lnhc_parameters_type
   USE gle_system_types,                ONLY: gle_dealloc,&
                                              gle_init,&
                                              gle_type
   USE input_constants,                 ONLY: do_region_global,&
                                              do_region_massive,&
                                              do_thermo_al,&
                                              do_thermo_csvr,&
                                              do_thermo_gle,&
                                              do_thermo_no_communication,&
                                              do_thermo_nose
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE simpar_types,                    ONLY: simpar_type
#include "../../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: thermostats_type, &
             thermostat_type, &
             allocate_thermostats, &
             release_thermostats, &
             create_thermostat_type, &
             release_thermostat_type, &
             thermostat_info_type, &
             release_thermostat_info, &
             set_thermostats

! **************************************************************************************************
!> \brief Define thermostat types
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      10.2007 created [tlaino] - Teodoro Laino - University of Zurich
!> \author Teodoro Laino
! **************************************************************************************************
   TYPE thermostats_type
      TYPE(thermostat_info_type), POINTER :: thermostat_info_part => NULL()
      TYPE(thermostat_info_type), POINTER :: thermostat_info_shell => NULL()
! cjm
      TYPE(thermostat_info_type), POINTER :: thermostat_info_fast => NULL()
      TYPE(thermostat_type), POINTER       :: thermostat_fast => NULL()
      TYPE(thermostat_info_type), POINTER :: thermostat_info_slow => NULL()
      TYPE(thermostat_type), POINTER       :: thermostat_slow => NULL()
! cjm
      TYPE(thermostat_type), POINTER       :: thermostat_part => NULL()
      TYPE(thermostat_type), POINTER       :: thermostat_coef => NULL()
      TYPE(thermostat_type), POINTER       :: thermostat_shell => NULL()
      TYPE(thermostat_type), POINTER       :: thermostat_baro => NULL()
   END TYPE thermostats_type

   ! Single thermostat_type
! **************************************************************************************************
   TYPE thermostat_type
      INTEGER :: type_of_thermostat = do_thermo_nose
      CHARACTER(LEN=default_string_length) :: label = ""
      TYPE(lnhc_parameters_type), POINTER :: nhc => NULL()
      TYPE(csvr_system_type), POINTER     :: csvr => NULL()
      TYPE(al_system_type), POINTER     :: al => NULL()
      TYPE(gle_type), POINTER      ::  gle => NULL()
      TYPE(section_vals_type), POINTER    :: section => NULL()
   END TYPE thermostat_type

   ! Global info type
! **************************************************************************************************
   TYPE thermostat_info_type
      INTEGER                              :: sum_of_thermostats = 0
      INTEGER                              :: number_of_thermostats = 0
      INTEGER                              :: dis_type = do_thermo_no_communication
      INTEGER, POINTER, DIMENSION(:)       :: map_loc_thermo_gen => NULL()
   END TYPE thermostat_info_type

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'thermostat_types'

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param thermostats ...
!> \par History
!>      09.2007 created [tlaino]
!> \author Teodoro Laino
! **************************************************************************************************
   SUBROUTINE allocate_thermostats(thermostats)
      TYPE(thermostats_type), INTENT(OUT)                :: thermostats

      ! Thermostats Info
      ALLOCATE (thermostats%thermostat_info_part)
      ALLOCATE (thermostats%thermostat_info_shell)
!cjm
      ALLOCATE (thermostats%thermostat_info_fast)
      ALLOCATE (thermostats%thermostat_info_slow)
!cjm

   END SUBROUTINE allocate_thermostats

! **************************************************************************************************
!> \brief ...
!> \param thermostats ...
!> \par History
!>      09.2007 created [tlaino]
!> \author Teodoro Laino
! **************************************************************************************************
   SUBROUTINE release_thermostats(thermostats)
      TYPE(thermostats_type), INTENT(INOUT)              :: thermostats

      IF (ASSOCIATED(thermostats%thermostat_info_part)) THEN
         CALL release_thermostat_info(thermostats%thermostat_info_part)
         DEALLOCATE (thermostats%thermostat_info_part)
      END IF
      IF (ASSOCIATED(thermostats%thermostat_info_shell)) THEN
         CALL release_thermostat_info(thermostats%thermostat_info_shell)
         DEALLOCATE (thermostats%thermostat_info_shell)
      END IF
      IF (ASSOCIATED(thermostats%thermostat_info_fast)) THEN
         CALL release_thermostat_info(thermostats%thermostat_info_fast)
         DEALLOCATE (thermostats%thermostat_info_fast)
      END IF
      IF (ASSOCIATED(thermostats%thermostat_info_slow)) THEN
         CALL release_thermostat_info(thermostats%thermostat_info_slow)
         DEALLOCATE (thermostats%thermostat_info_slow)
      END IF
      IF (ASSOCIATED(thermostats%thermostat_fast)) THEN
         CALL release_thermostat_type(thermostats%thermostat_fast)
         DEALLOCATE (thermostats%thermostat_fast)
      END IF
      IF (ASSOCIATED(thermostats%thermostat_slow)) THEN
         CALL release_thermostat_type(thermostats%thermostat_slow)
         DEALLOCATE (thermostats%thermostat_slow)
      END IF
      IF (ASSOCIATED(thermostats%thermostat_part)) THEN
         CALL release_thermostat_type(thermostats%thermostat_part)
         DEALLOCATE (thermostats%thermostat_part)
      END IF
      IF (ASSOCIATED(thermostats%thermostat_shell)) THEN
         CALL release_thermostat_type(thermostats%thermostat_shell)
         DEALLOCATE (thermostats%thermostat_shell)
      END IF
      IF (ASSOCIATED(thermostats%thermostat_baro)) THEN
         CALL release_thermostat_type(thermostats%thermostat_baro)
         DEALLOCATE (thermostats%thermostat_baro)
      END IF
      IF (ASSOCIATED(thermostats%thermostat_coef)) THEN
         CALL release_thermostat_type(thermostats%thermostat_coef)
         DEALLOCATE (thermostats%thermostat_coef)
      END IF

   END SUBROUTINE release_thermostats

! **************************************************************************************************
!> \brief Create a thermostat type
!> \param thermostat ...
!> \param simpar ...
!> \param section ...
!> \param skip_region ...
!> \param label ...
!> \par History
!>      09.2007 created [tlaino]
!> \author Teodoro Laino
! **************************************************************************************************
   SUBROUTINE create_thermostat_type(thermostat, simpar, section, skip_region, label)
      TYPE(thermostat_type), INTENT(OUT)                 :: thermostat
      TYPE(simpar_type), POINTER                         :: simpar
      TYPE(section_vals_type), POINTER                   :: section
      LOGICAL, INTENT(IN), OPTIONAL                      :: skip_region
      CHARACTER(LEN=*), INTENT(IN)                       :: label

      INTEGER                                            :: region
      LOGICAL                                            :: skip_region_loc
      TYPE(section_vals_type), POINTER                   :: al_section, csvr_section, gle_section, &
                                                            nose_section

      skip_region_loc = .FALSE.
      IF (PRESENT(skip_region)) skip_region_loc = skip_region
      thermostat%section => section
      thermostat%label = label
      region = do_region_global

      CALL section_vals_val_get(section, "TYPE", i_val=thermostat%type_of_thermostat)
      IF (.NOT. skip_region_loc) CALL section_vals_val_get(section, "REGION", i_val=region)
      IF (thermostat%type_of_thermostat == do_thermo_nose) THEN
         nose_section => section_vals_get_subs_vals(section, "NOSE")
         ALLOCATE (thermostat%nhc)
         CALL lnhc_init(thermostat%nhc, nose_section)
         thermostat%nhc%region = region
      ELSE IF (thermostat%type_of_thermostat == do_thermo_csvr) THEN
         csvr_section => section_vals_get_subs_vals(section, "CSVR")
         ALLOCATE (thermostat%csvr)
         CALL csvr_init(thermostat%csvr, simpar, csvr_section)
         thermostat%csvr%region = region
      ELSE IF (thermostat%type_of_thermostat == do_thermo_al) THEN
         al_section => section_vals_get_subs_vals(section, "AD_LANGEVIN")
         ALLOCATE (thermostat%al)
         CALL al_init(thermostat%al, simpar, al_section)
         thermostat%al%region = region
      ELSE IF (thermostat%type_of_thermostat == do_thermo_gle) THEN
         gle_section => section_vals_get_subs_vals(section, "GLE")
         ALLOCATE (thermostat%gle)
         CALL gle_init(thermostat%gle, dt=simpar%dt, temp=simpar%temp_ext, &
                       section=gle_section)
         thermostat%gle%region = region
         CPASSERT(region == do_region_massive)
      END IF

   END SUBROUTINE create_thermostat_type

! **************************************************************************************************
!> \brief ...
!> \param thermostat_info ...
!> \par History
!>      10.2007 created [tlaino]
!> \author Teodoro Laino
! **************************************************************************************************
   SUBROUTINE release_thermostat_info(thermostat_info)
      TYPE(thermostat_info_type), INTENT(INOUT)          :: thermostat_info

      IF (ASSOCIATED(thermostat_info%map_loc_thermo_gen)) THEN
         DEALLOCATE (thermostat_info%map_loc_thermo_gen)
      END IF

   END SUBROUTINE release_thermostat_info

! **************************************************************************************************
!> \brief ...
!> \param thermostat ...
!> \par History
!>      09.2007 created [tlaino]
!> \author Teodoro Laino
! **************************************************************************************************
   SUBROUTINE release_thermostat_type(thermostat)
      TYPE(thermostat_type), INTENT(INOUT)               :: thermostat

      NULLIFY (thermostat%section)
      IF (ASSOCIATED(thermostat%nhc)) THEN
         CALL lnhc_dealloc(thermostat%nhc)
      END IF
      IF (ASSOCIATED(thermostat%csvr)) THEN
         CALL csvr_dealloc(thermostat%csvr)
      END IF
      IF (ASSOCIATED(thermostat%al)) THEN
         CALL al_dealloc(thermostat%al)
      END IF
      IF (ASSOCIATED(thermostat%gle)) THEN
         CALL gle_dealloc(thermostat%gle)
      END IF

   END SUBROUTINE release_thermostat_type

! **************************************************************************************************
!> \brief access internal structures of thermostats
!> \param thermostats ...
!> \param dt_fact ...
!> \par History
!>      10.2008 created [tlaino]
!> \author Teodoro Laino [tlaino] - University of Zurich
! **************************************************************************************************
   SUBROUTINE set_thermostats(thermostats, dt_fact)
      TYPE(thermostats_type), POINTER                    :: thermostats
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: dt_fact

      IF (ASSOCIATED(thermostats)) THEN
         IF (PRESENT(dt_fact)) THEN
            ! Particles
!cjm
            IF (ASSOCIATED(thermostats%thermostat_fast)) THEN
               SELECT CASE (thermostats%thermostat_fast%type_of_thermostat)
               CASE (do_thermo_nose)
                  thermostats%thermostat_fast%nhc%dt_fact = dt_fact
               END SELECT
            END IF
            IF (ASSOCIATED(thermostats%thermostat_slow)) THEN
               SELECT CASE (thermostats%thermostat_slow%type_of_thermostat)
               CASE (do_thermo_nose)
                  thermostats%thermostat_slow%nhc%dt_fact = dt_fact
               END SELECT
            END IF
!cjm
            IF (ASSOCIATED(thermostats%thermostat_part)) THEN
               SELECT CASE (thermostats%thermostat_part%type_of_thermostat)
               CASE (do_thermo_nose)
                  thermostats%thermostat_part%nhc%dt_fact = dt_fact
               CASE (do_thermo_csvr)
                  thermostats%thermostat_part%csvr%dt_fact = dt_fact
               CASE (do_thermo_al)
                  thermostats%thermostat_part%al%dt_fact = dt_fact
               CASE (do_thermo_gle)
                  thermostats%thermostat_part%gle%dt_fact = dt_fact
               END SELECT
            END IF
            ! Coefficients
            IF (ASSOCIATED(thermostats%thermostat_coef)) THEN
               SELECT CASE (thermostats%thermostat_coef%type_of_thermostat)
               CASE (do_thermo_nose)
                  thermostats%thermostat_coef%nhc%dt_fact = dt_fact
               CASE (do_thermo_csvr)
                  thermostats%thermostat_coef%csvr%dt_fact = dt_fact
               END SELECT
            END IF
            ! Shell
            IF (ASSOCIATED(thermostats%thermostat_shell)) THEN
               SELECT CASE (thermostats%thermostat_shell%type_of_thermostat)
               CASE (do_thermo_nose)
                  thermostats%thermostat_shell%nhc%dt_fact = dt_fact
               CASE (do_thermo_csvr)
                  thermostats%thermostat_shell%csvr%dt_fact = dt_fact
               END SELECT
            END IF
            ! Baro
            IF (ASSOCIATED(thermostats%thermostat_baro)) THEN
               SELECT CASE (thermostats%thermostat_baro%type_of_thermostat)
               CASE (do_thermo_nose)
                  thermostats%thermostat_baro%nhc%dt_fact = dt_fact
               CASE (do_thermo_csvr)
                  thermostats%thermostat_baro%csvr%dt_fact = dt_fact
               END SELECT
            END IF
         END IF
      END IF
   END SUBROUTINE set_thermostats

END MODULE thermostat_types
