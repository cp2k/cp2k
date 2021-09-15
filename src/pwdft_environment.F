!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Methods and functions on the PWDFT environment
!> \par History
!>      07.2018 initial create
!> \author JHU
! **************************************************************************************************
MODULE pwdft_environment
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cell_types,                      ONLY: cell_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE cp_symmetry,                     ONLY: write_symmetry
   USE distribution_1d_types,           ONLY: distribution_1d_release,&
                                              distribution_1d_type
   USE distribution_methods,            ONLY: distribute_molecules_1d
   USE header,                          ONLY: sirius_header
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: m_flush
   USE molecule_kind_types,             ONLY: molecule_kind_type
   USE molecule_types,                  ONLY: molecule_type
   USE particle_methods,                ONLY: write_particle_distances,&
                                              write_qs_particle_coordinates,&
                                              write_structure_data
   USE particle_types,                  ONLY: particle_type
   USE pwdft_environment_types,         ONLY: pwdft_env_get,&
                                              pwdft_env_set,&
                                              pwdft_environment_type
   USE qs_energy_types,                 ONLY: allocate_qs_energy,&
                                              qs_energy_type
   USE qs_kind_types,                   ONLY: qs_kind_type,&
                                              write_qs_kind_set
   USE qs_subsys_methods,               ONLY: qs_subsys_create
   USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                              qs_subsys_release,&
                                              qs_subsys_set,&
                                              qs_subsys_type
   USE sirius_interface,                ONLY: cp_sirius_create_env,&
                                              cp_sirius_energy_force,&
                                              cp_sirius_update_context
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pwdft_environment'

! *** Public subroutines ***

   PUBLIC :: pwdft_init, pwdft_calc_energy_force

CONTAINS

! **************************************************************************************************
!> \brief Initialize the pwdft environment
!> \param pwdft_env The pwdft environment to retain
!> \param root_section ...
!> \param para_env ...
!> \param force_env_section ...
!> \param subsys_section ...
!> \param use_motion_section ...
!> \par History
!>      03.2018 initial create
!> \author JHU
! **************************************************************************************************
   SUBROUTINE pwdft_init(pwdft_env, root_section, para_env, force_env_section, subsys_section, &
                         use_motion_section)
      TYPE(pwdft_environment_type), POINTER              :: pwdft_env
      TYPE(section_vals_type), POINTER                   :: root_section
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(section_vals_type), POINTER                   :: force_env_section, subsys_section
      LOGICAL, INTENT(IN)                                :: use_motion_section

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pwdft_init'

      INTEGER                                            :: handle, iw, natom
      LOGICAL                                            :: use_ref_cell
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: my_cell, my_cell_ref
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(distribution_1d_type), POINTER                :: local_molecules, local_particles
      TYPE(molecule_kind_type), DIMENSION(:), POINTER    :: molecule_kind_set
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_subsys_type), POINTER                      :: qs_subsys
      TYPE(section_vals_type), POINTER                   :: pwdft_section, xc_section

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(pwdft_env))
      CPASSERT(ASSOCIATED(force_env_section))

      IF (.NOT. ASSOCIATED(subsys_section)) THEN
         subsys_section => section_vals_get_subs_vals(force_env_section, "SUBSYS")
      END IF
      pwdft_section => section_vals_get_subs_vals(force_env_section, "PW_DFT")

      ! retrieve the functionals parameters
      xc_section => section_vals_get_subs_vals(force_env_section, "DFT%XC%XC_FUNCTIONAL")

      CALL pwdft_env_set(pwdft_env=pwdft_env, pwdft_input=pwdft_section, &
                         force_env_input=force_env_section, xc_input=xc_section)

      ! parallel environment
      CALL pwdft_env_set(pwdft_env=pwdft_env, para_env=para_env)

      NULLIFY (qs_subsys)
      CALL qs_subsys_create(qs_subsys, para_env, &
                            force_env_section=force_env_section, &
                            subsys_section=subsys_section, &
                            use_motion_section=use_motion_section, &
                            root_section=root_section)

      ! Distribute molecules and atoms
      NULLIFY (local_molecules)
      NULLIFY (local_particles)
      CALL qs_subsys_get(qs_subsys, particle_set=particle_set, &
                         atomic_kind_set=atomic_kind_set, &
                         molecule_set=molecule_set, &
                         molecule_kind_set=molecule_kind_set)

      CALL distribute_molecules_1d(atomic_kind_set=atomic_kind_set, &
                                   particle_set=particle_set, &
                                   local_particles=local_particles, &
                                   molecule_kind_set=molecule_kind_set, &
                                   molecule_set=molecule_set, &
                                   local_molecules=local_molecules, &
                                   force_env_section=force_env_section)

      CALL qs_subsys_set(qs_subsys, local_molecules=local_molecules, &
                         local_particles=local_particles)
      CALL distribution_1d_release(local_particles)
      CALL distribution_1d_release(local_molecules)

      CALL pwdft_env_set(pwdft_env=pwdft_env, qs_subsys=qs_subsys)
      CALL qs_subsys_release(qs_subsys)

      CALL pwdft_env_get(pwdft_env=pwdft_env, qs_subsys=qs_subsys)
      CALL qs_subsys_get(qs_subsys, cell=my_cell, cell_ref=my_cell_ref, use_ref_cell=use_ref_cell)

      NULLIFY (logger)
      logger => cp_get_default_logger()
      iw = cp_logger_get_default_io_unit(logger)
      CALL sirius_header(iw)

      NULLIFY (energy)
      CALL allocate_qs_energy(energy)
      CALL qs_subsys_set(qs_subsys, energy=energy)

      CALL qs_subsys_get(qs_subsys, particle_set=particle_set, &
                         qs_kind_set=qs_kind_set, &
                         atomic_kind_set=atomic_kind_set, &
                         molecule_set=molecule_set, &
                         molecule_kind_set=molecule_kind_set)

      ! init energy, force, stress arrays
      CALL qs_subsys_get(qs_subsys, natom=natom)
      ALLOCATE (pwdft_env%energy)
      ALLOCATE (pwdft_env%forces(natom, 3))
      pwdft_env%forces(1:natom, 1:3) = 0.0_dp
      pwdft_env%stress(1:3, 1:3) = 0.0_dp

      ! Print the atomic kind set
      CALL write_qs_kind_set(qs_kind_set, subsys_section)
      ! Print the atomic coordinates
      CALL write_qs_particle_coordinates(particle_set, qs_kind_set, subsys_section, label="PWDFT/SIRIUS")
      ! Print the interatomic distances
      CALL write_particle_distances(particle_set, my_cell, subsys_section)
      ! Print the requested structure data
      CALL write_structure_data(particle_set, my_cell, subsys_section)
      ! Print symmetry information
      CALL write_symmetry(particle_set, my_cell, subsys_section)

      ! Sirius initialization
      CALL cp_sirius_create_env(pwdft_env)

      IF (iw > 0) THEN
         WRITE (iw, '(A)') "SIRIUS INIT: FINISHED"
      END IF
      IF (iw > 0) CALL m_flush(iw)

      CALL timestop(handle)

   END SUBROUTINE pwdft_init

! **************************************************************************************************
!> \brief Calculate energy and forces within the PWDFT/SIRIUS code
!> \param pwdft_env The pwdft environment to retain
!> \param calculate_forces ...
!> \param calculate_stress ...
!> \par History
!>      03.2018 initial create
!> \author JHU
! **************************************************************************************************
   SUBROUTINE pwdft_calc_energy_force(pwdft_env, calculate_forces, calculate_stress)
      TYPE(pwdft_environment_type), POINTER              :: pwdft_env
      LOGICAL, INTENT(IN)                                :: calculate_forces, calculate_stress

      CHARACTER(len=*), PARAMETER :: routineN = 'pwdft_calc_energy_force'

      INTEGER                                            :: handle, iatom, iw, natom
      REAL(KIND=dp), DIMENSION(1:3, 1:3)                 :: stress
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: force
      TYPE(cell_type), POINTER                           :: my_cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_subsys_type), POINTER                      :: qs_subsys
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(pwdft_env))

      NULLIFY (logger)
      logger => cp_get_default_logger()
      iw = cp_logger_get_default_io_unit(logger)

      ! update context for new positions/cell
      ! Sirius initialization
      CALL cp_sirius_update_context(pwdft_env)
      IF (iw > 0) THEN
         WRITE (iw, '(A)') "SIRIUS UPDATE CONTEXT : FINISHED"
      END IF
      IF (iw > 0) CALL m_flush(iw)

      ! calculate energy and forces/stress
      CALL cp_sirius_energy_force(pwdft_env, calculate_forces, calculate_stress)

      IF (calculate_forces) THEN
         CALL pwdft_env_get(pwdft_env=pwdft_env, qs_subsys=qs_subsys)
         CALL pwdft_env_get(pwdft_env=pwdft_env, forces=force)
         CALL qs_subsys_get(qs_subsys, particle_set=particle_set, natom=natom)
         DO iatom = 1, natom
            particle_set(iatom)%f(1:3) = -force(iatom, 1:3)
         END DO
      END IF

      IF (calculate_stress) THEN
         ! i need to retrieve the volume of the unit cell for the stress tensor
         CALL qs_subsys_get(qs_subsys, cell=my_cell, virial=virial)
         CALL pwdft_env_get(pwdft_env=pwdft_env, stress=stress)
         virial%pv_virial(1:3, 1:3) = -stress(1:3, 1:3)*my_cell%deth
      END IF

      CALL timestop(handle)

   END SUBROUTINE pwdft_calc_energy_force
END MODULE pwdft_environment
