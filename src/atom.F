!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
MODULE atom
   USE atom_basis,                      ONLY: atom_basis_opt
   USE atom_energy,                     ONLY: atom_energy_opt
   USE atom_pseudo,                     ONLY: atom_pseudo_opt
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE header,                          ONLY: atom_footer,&
                                              atom_header
   USE input_constants,                 ONLY: atom_basis_run,&
                                              atom_energy_run,&
                                              atom_no_run,&
                                              atom_pseudo_run
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length
   USE periodic_table,                  ONLY: nelem,&
                                              ptable
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC  :: atom_code

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'atom'

CONTAINS

! **************************************************************************************************
!> \brief Driver routine to perform atomic calculations.
!> \param root_section  root input section
!> \par History
!>    * 08.2008 created [Juerg Hutter]
! **************************************************************************************************
   SUBROUTINE atom_code(root_section)
      TYPE(section_vals_type), POINTER                   :: root_section

      CHARACTER(len=*), PARAMETER                        :: routineN = 'atom_code'

      INTEGER                                            :: handle, iw, run_type_id
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: atom_section

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      NULLIFY (atom_section)
      atom_section => section_vals_get_subs_vals(root_section, "ATOM")

      iw = cp_print_key_unit_nr(logger, atom_section, "PRINT%PROGRAM_BANNER", extension=".log")
      CALL atom_header(iw)
      CALL cp_print_key_finished_output(iw, logger, atom_section, "PRINT%PROGRAM_BANNER")

      CALL atom_test(atom_section)

      CALL section_vals_val_get(atom_section, "RUN_TYPE", i_val=run_type_id)
      SELECT CASE (run_type_id)
      CASE (atom_no_run)
         ! do (almost) nothing
      CASE (atom_energy_run)
         CALL atom_energy_opt(atom_section)
      CASE (atom_basis_run)
         CALL atom_basis_opt(atom_section)
      CASE (atom_pseudo_run)
         CALL atom_pseudo_opt(atom_section)
      CASE default
         CPABORT("")
      END SELECT

      iw = cp_print_key_unit_nr(logger, atom_section, "PRINT%PROGRAM_BANNER", extension=".log")
      CALL atom_footer(iw)
      CALL cp_print_key_finished_output(iw, logger, atom_section, "PRINT%PROGRAM_BANNER")

      CALL timestop(handle)

   END SUBROUTINE atom_code

! **************************************************************************************************
!> \brief Check consistency between the element symbol and its atomic number.
!> \param atom_section  ATOM input section
!> \par History
!>    * 08.2008 created [Juerg Hutter]
! **************************************************************************************************
   SUBROUTINE atom_test(atom_section)
      TYPE(section_vals_type), POINTER                   :: atom_section

      CHARACTER(len=*), PARAMETER                        :: routineN = 'atom_test'

      CHARACTER(len=2)                                   :: elem
      CHARACTER(len=default_string_length)               :: z_string
      INTEGER                                            :: handle, i, z
      LOGICAL                                            :: explicit_elem, explicit_z

      CALL timeset(routineN, handle)

      CALL section_vals_val_get(atom_section, "ATOMIC_NUMBER", i_val=z, explicit=explicit_z)
      CALL section_vals_val_get(atom_section, "ELEMENT", c_val=elem, explicit=explicit_elem)

      IF (explicit_z .AND. (z <= 0 .AND. z > nelem)) THEN
         ! an explicit atomic number is not found in the periodic table
         WRITE (z_string, '(I0)') z
         CALL cp_abort(__LOCATION__, &
                       "The element with the atomic number "//TRIM(z_string)//" is not found in the periodic table.")
      END IF

      IF (explicit_elem) THEN
         ! check that the element symbol is part of the periodic table
         DO i = 1, nelem
            IF (ptable(i)%symbol == elem) EXIT
         END DO

         IF (i > nelem) THEN
            CALL cp_abort(__LOCATION__, &
                          "The element symbol ("//TRIM(elem)//") is not found in the periodic table.")
         END IF
      END IF

      IF (explicit_z .AND. explicit_elem) THEN
         ! check that the element symbol read from the input file
         ! matches for the explicitly given atomic number
         IF (ptable(z)%symbol /= elem) THEN
            WRITE (z_string, '(I0)') z
            CALL cp_abort(__LOCATION__, &
                          "The element symbol ("//TRIM(elem)// &
                          ") contradicts with the explicitly given atomic number ("// &
                          TRIM(z_string)//").")
         END IF
      ELSE IF (.NOT. (explicit_z .OR. explicit_elem)) THEN
         ! default (implicit) element symbol and atomic number are usually consistent,
         ! but check them just in case
         CPASSERT(ptable(z)%symbol == elem)
      END IF

      CALL timestop(handle)
   END SUBROUTINE atom_test

END MODULE atom
