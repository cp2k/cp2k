!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Atomic Polarization Tensor calculation by dF/d(E-field) finite differences
!> \author Leo Decking, Hossam Elgabarty
! **************************************************************************************************

MODULE qs_apt_fdiff_methods
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_log_handling,                 ONLY: cp_add_default_logger,&
                                              cp_get_default_logger,&
                                              cp_logger_create,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_release,&
                                              cp_logger_set,&
                                              cp_logger_type,&
                                              cp_rm_default_logger
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                              cp_subsys_type
   USE force_env_types,                 ONLY: force_env_get,&
                                              force_env_type
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE particle_list_types,             ONLY: particle_list_type
   USE qs_apt_fdiff_types,              ONLY: apt_fdiff_point_type,&
                                              apt_fdiff_points_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force,                        ONLY: qs_calc_energy_force
   USE qs_subsys_types,                 ONLY: qs_subsys_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_apt_fdiff_methods'
   LOGICAL, PARAMETER, PRIVATE :: debug_this_module = .FALSE.
   PUBLIC :: apt_fdiff

CONTAINS

! **************************************************************************************************
!> \brief Perform the 2PNT symmetric difference
!> \param apt_fdiff_points ...
!> \param ap_tensor ...
!> \param natoms ...
!> \author Leo Decking
! **************************************************************************************************
   SUBROUTINE fdiff_2pnt(apt_fdiff_points, ap_tensor, natoms)
      TYPE(apt_fdiff_points_type)                        :: apt_fdiff_points
      REAL(kind=dp), DIMENSION(:, :, :)                  :: ap_tensor
      INTEGER, INTENT(IN)                                :: natoms

      INTEGER                                            :: i, j, n

      DO j = 1, 3 ! axis force
         DO i = 1, 3 ! axis field
            DO n = 1, natoms
               ap_tensor(n, i, j) = (apt_fdiff_points%point_field(i, 1)%forces(n, j) - &
                                     apt_fdiff_points%point_field(i, 2)%forces(n, j)) &
                                    /(2*apt_fdiff_points%field_strength)
            END DO
         END DO
      END DO

   END SUBROUTINE fdiff_2pnt

! **************************************************************************************************
!> \brief ...
!> \param apt_fdiff_point ...
!> \param particles ...
!> \author Leo Decking
! **************************************************************************************************
   SUBROUTINE get_forces(apt_fdiff_point, particles)
      TYPE(apt_fdiff_point_type)                         :: apt_fdiff_point
      TYPE(particle_list_type), POINTER                  :: particles

      INTEGER                                            :: i

      CPASSERT(ASSOCIATED(particles))

      DO i = 1, particles%n_els
         apt_fdiff_point%forces(i, 1:3) = particles%els(i)%f(1:3)
      END DO

   END SUBROUTINE get_forces

! **************************************************************************************************
!> \brief Calculate Atomic Polarization Tensors by dF/d(E-field) finite differences
!> \param force_env ...
!> \author Leo Decking, Hossam Elgabarty
! **************************************************************************************************
   SUBROUTINE apt_fdiff(force_env)

      TYPE(force_env_type), POINTER                      :: force_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'apt_fdiff'

      INTEGER                                            :: apt_log, fd_method, handle, i, j, &
                                                            log_unit, n, natoms, output_fdiff_scf
      REAL(kind=dp)                                      :: born_sum
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: ap_tensor
      TYPE(apt_fdiff_points_type)                        :: apt_fdiff_points
      TYPE(cp_logger_type), POINTER                      :: logger, logger_apt
      TYPE(cp_subsys_type), POINTER                      :: cp_subsys
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_subsys_type), POINTER                      :: subsys
      TYPE(section_vals_type), POINTER                   :: dcdr_section

      CALL timeset(routineN, handle)

      NULLIFY (qs_env, logger, dcdr_section, logger_apt, para_env, dft_control, subsys)
      NULLIFY (particles)

      CPASSERT(ASSOCIATED(force_env))

      CALL force_env_get(force_env, subsys=cp_subsys, qs_env=qs_env, para_env=para_env)

      CALL cp_subsys_get(cp_subsys, particles=particles)
      CPASSERT(ASSOCIATED(particles))
      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control)!, para_env=para_env)

      natoms = particles%n_els
      IF (dft_control%apply_period_efield .OR. dft_control%apply_efield) THEN
         CPABORT("APT calculation not available in the presence of an external electric field")
      END IF

      logger => cp_get_default_logger()

      log_unit = cp_logger_get_default_io_unit(logger)

      dcdr_section => section_vals_get_subs_vals(qs_env%input, "PROPERTIES%LINRES%DCDR")

      apt_log = cp_print_key_unit_nr(logger, dcdr_section, "PRINT%APT", &
                                     extension=".data", middle_name="apt", log_filename=.FALSE., &
                                     file_position="APPEND", file_status="UNKNOWN")

      output_fdiff_scf = cp_print_key_unit_nr(logger, dcdr_section, "PRINT%APT", &
                                              extension=".scfLog", middle_name="apt", log_filename=.FALSE., &
                                              file_position="APPEND", file_status="UNKNOWN")

      CALL cp_logger_create(logger_apt, para_env=para_env, print_level=0, &
                            default_global_unit_nr=output_fdiff_scf, &
                            close_global_unit_on_dealloc=.TRUE.)

      CALL cp_logger_set(logger_apt, global_filename="APT_localLog")

      CALL cp_add_default_logger(logger_apt)

      IF (output_fdiff_scf > 0) THEN
         WRITE (output_fdiff_scf, '(T2,A)') &
            '!----------------------------------------------------------------------------!'
         WRITE (output_fdiff_scf, '(/,T2,A)') "SCF log for finite difference steps"
         WRITE (output_fdiff_scf, '(/,T2,A)') &
            '!----------------------------------------------------------------------------!'
      END IF

      IF (log_unit > 0) THEN
         WRITE (log_unit, '(/,T2,A)') &
            '!----------------------------------------------------------------------------!'
         WRITE (log_unit, '(/,T10,A)') "Computing Atomic polarization tensors using finite differences"
         WRITE (log_unit, '(T2,A)') "  "
      END IF

      CALL section_vals_val_get(dcdr_section, "APT_FD_DE", r_val=apt_fdiff_points%field_strength)
      CALL section_vals_val_get(dcdr_section, "APT_FD_METHOD", i_val=fd_method)

      dft_control%apply_period_efield = .TRUE.
      ALLOCATE (dft_control%period_efield)
      dft_control%period_efield%displacement_field = .FALSE.

      DO i = 1, 3
         dft_control%period_efield%polarisation(1:3) = (/0.0_dp, 0.0_dp, 0.0_dp/)
         dft_control%period_efield%polarisation(i) = 1.0_dp

         IF (log_unit > 0) THEN
            WRITE (log_unit, '(T2,A)') "  "
            WRITE (log_unit, "(T2,A)") "Computing forces under efield in direction +/-"//ACHAR(i + 119)
         END IF

         DO j = 1, 2 ! 1 -> positive, 2 -> negative
            IF (j == 1) THEN
               dft_control%period_efield%strength = apt_fdiff_points%field_strength
            ELSE
               dft_control%period_efield%strength = -1.0*apt_fdiff_points%field_strength
            END IF

            CALL qs_calc_energy_force(qs_env, calc_force=.TRUE., consistent_energies=.TRUE., linres=.FALSE.)

            ALLOCATE (apt_fdiff_points%point_field(i, j)%forces(natoms, 1:3))
            CALL get_forces(apt_fdiff_points%point_field(i, j), particles)

         END DO
      END DO

      IF (output_fdiff_scf > 0) THEN
         WRITE (output_fdiff_scf, '(/,T2,A)') &
            '!----------------------------------------------------------------------------!'
         WRITE (output_fdiff_scf, '(/,T2,A)') "Finite differences done!"
         WRITE (output_fdiff_scf, '(/,T2,A)') &
            '!----------------------------------------------------------------------------!'
      END IF

      CALL cp_print_key_finished_output(output_fdiff_scf, logger_apt, dcdr_section, "PRINT%APT")
      CALL cp_logger_release(logger_apt)
      CALL cp_rm_default_logger()

      ALLOCATE (ap_tensor(natoms, 3, 3))
      CALL fdiff_2pnt(apt_fdiff_points, ap_tensor, natoms)

      !Print
      born_sum = 0.0_dp
      IF (apt_log > 0) THEN
         DO n = 1, natoms
            born_sum = born_sum + (ap_tensor(n, 1, 1) + ap_tensor(n, 1, 1) + ap_tensor(n, 1, 1))/3.0
            WRITE (apt_log, "(I6, A6, F20.10)") n, particles%els(n)%atomic_kind%element_symbol, &
               (ap_tensor(n, 1, 1) + ap_tensor(n, 2, 2) + ap_tensor(n, 3, 3))/3.0
            WRITE (apt_log, '(F20.10,F20.10,F20.10)') &
               ap_tensor(n, 1, 1), ap_tensor(n, 1, 2), ap_tensor(n, 1, 3)
            WRITE (apt_log, '(F20.10,F20.10,F20.10)') &
               ap_tensor(n, 2, 1), ap_tensor(n, 2, 2), ap_tensor(n, 2, 3)
            WRITE (apt_log, '(F20.10,F20.10,F20.10)') &
               ap_tensor(n, 3, 1), ap_tensor(n, 3, 2), ap_tensor(n, 3, 3)
         END DO
         WRITE (apt_log, '(/,A20, F20.10)') "Sum of Born charges:", born_sum
         WRITE (log_unit, '(/,A30, F20.10)') "Checksum (Acoustic Sum Rule):", born_sum
      END IF
      DEALLOCATE (ap_tensor)
      DEALLOCATE (dft_control%period_efield)

      CALL cp_print_key_finished_output(apt_log, logger, dcdr_section, "PRINT%APT")

      dft_control%apply_period_efield = .FALSE.
      qs_env%linres_run = .FALSE.

      IF (log_unit > 0) THEN
         WRITE (log_unit, '(T2,A)') "  "
         WRITE (log_unit, '(T2,A)') "  "
         WRITE (log_unit, '(T22,A)') "APT calculation Done!"
         WRITE (log_unit, '(/,T2,A)') &
            '!----------------------------------------------------------------------------!'
      END IF

      CALL timestop(handle)

   END SUBROUTINE apt_fdiff

END MODULE qs_apt_fdiff_methods
