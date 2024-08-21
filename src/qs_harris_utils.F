!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Harris method environment setup and handling
!> \par History
!>       2024.07 created
!> \author JGH
! **************************************************************************************************
MODULE qs_harris_utils
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE input_constants,                 ONLY: hden_atomic,&
                                              hfun_harris,&
                                              horb_default
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get
   USE qs_environment_types,            ONLY: qs_environment_type
   USE qs_harris_types,                 ONLY: harris_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_harris_utils'

   PUBLIC :: harris_env_create
   PUBLIC :: harris_write_input

CONTAINS

! **************************************************************************************************
!> \brief Allocates and intitializes harris_env
!> \param qs_env The QS environment
!> \param harris_env The Harris method environment (the object to create)
!> \param harris_section The Harris method input section
!> \par History
!>       2024.07 created
!> \author JGH
! **************************************************************************************************
   SUBROUTINE harris_env_create(qs_env, harris_env, harris_section)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(harris_type), POINTER                         :: harris_env
      TYPE(section_vals_type), OPTIONAL, POINTER         :: harris_section

      CPASSERT(.NOT. ASSOCIATED(harris_env))
      ALLOCATE (harris_env)
      CALL init_harris_env(qs_env, harris_env, harris_section)

   END SUBROUTINE harris_env_create

! **************************************************************************************************
!> \brief Initializes Harris method environment
!> \param qs_env The QS environment
!> \param harris_env The Harris method environment
!> \param harris_section The Harris method input section
!> \par History
!>       2024.07 created
!> \author JGH
! **************************************************************************************************
   SUBROUTINE init_harris_env(qs_env, harris_env, harris_section)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(harris_type), POINTER                         :: harris_env
      TYPE(section_vals_type), OPTIONAL, POINTER         :: harris_section

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'init_harris_env'

      INTEGER                                            :: handle, unit_nr
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)

      IF (qs_env%harris_method) THEN

         CPASSERT(PRESENT(harris_section))
         ! get a useful output_unit
         logger => cp_get_default_logger()
         IF (logger%para_env%is_source()) THEN
            unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
         ELSE
            unit_nr = -1
         END IF

         CALL section_vals_val_get(harris_section, "ENERGY_FUNCTIONAL", &
                                   i_val=harris_env%energy_functional)
         CALL section_vals_val_get(harris_section, "DENSITY_SOURCE", &
                                   i_val=harris_env%density_source)
         CALL section_vals_val_get(harris_section, "ORBITAL_BASIS", &
                                   i_val=harris_env%orbital_basis)
         !
         CALL section_vals_val_get(harris_section, "DEBUG_FORCES", &
                                   l_val=harris_env%debug_forces)
         CALL section_vals_val_get(harris_section, "DEBUG_STRESS", &
                                   l_val=harris_env%debug_stress)

      END IF

      CALL timestop(handle)

   END SUBROUTINE init_harris_env

! **************************************************************************************************
!> \brief Print out the Harris method input section
!>
!> \param harris_env ...
!> \par History
!>       2024.07 created [JGH]
!> \author JGH
! **************************************************************************************************
   SUBROUTINE harris_write_input(harris_env)
      TYPE(harris_type), POINTER                         :: harris_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'harris_write_input'

      INTEGER                                            :: handle, unit_nr
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      IF (logger%para_env%is_source()) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      END IF

      IF (unit_nr > 0) THEN

         WRITE (unit_nr, '(/,T2,A)') &
            "!"//REPEAT("-", 29)//"   Harris Model    "//REPEAT("-", 29)//"!"

         ! Type of energy functional
         SELECT CASE (harris_env%energy_functional)
         CASE (hfun_harris)
            WRITE (unit_nr, '(T2,A,T61,A20)') "Energy Functional: ", "Harris"
         END SELECT
         ! density source
         SELECT CASE (harris_env%density_source)
         CASE (hden_atomic)
            WRITE (unit_nr, '(T2,A,T61,A20)') "Harris model density: ", " Atomic kind density"
         END SELECT
         ! orbital basis
         SELECT CASE (harris_env%orbital_basis)
         CASE (horb_default)
            WRITE (unit_nr, '(T2,A,T61,A20)') "Harris model basis: ", "Atomic kind orbitals"
         END SELECT

         WRITE (unit_nr, '(T2,A)') REPEAT("-", 79)
         WRITE (unit_nr, '()')

      END IF ! unit_nr

      CALL timestop(handle)

   END SUBROUTINE harris_write_input

END MODULE qs_harris_utils
