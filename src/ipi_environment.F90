!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Methods and functions on the i–PI environment
!> \par History
!>      03.2024 initial create
!> \author Sebastian Seidenath (sebastian.seidenath@uni-jena.de)
! **************************************************************************************************
MODULE ipi_environment
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cell_types,                      ONLY: cell_type
   USE cp_subsys_methods,               ONLY: cp_subsys_create
   USE cp_subsys_types,                 ONLY: cp_subsys_type
   USE distribution_1d_types,           ONLY: distribution_1d_release,&
                                              distribution_1d_type
   USE distribution_methods,            ONLY: distribute_molecules_1d
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE ipi_environment_types,           ONLY: ipi_env_set,&
                                              ipi_environment_type
   USE ipi_server,                      ONLY: start_server
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE molecule_kind_types,             ONLY: molecule_kind_type,&
                                              write_molecule_kind_set
   USE molecule_types,                  ONLY: molecule_type
   USE particle_methods,                ONLY: write_fist_particle_coordinates,&
                                              write_particle_distances,&
                                              write_structure_data
   USE particle_types,                  ONLY: particle_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ipi_environment'

! *** Public subroutines ***

   PUBLIC :: ipi_init

CONTAINS

! **************************************************************************************************
!> \brief Initialize the ipi environment
!> \param ipi_env The ipi environment to retain
!> \param root_section ...
!> \param para_env ...
!> \param force_env_section ...
!> \param subsys_section ...
!> \par History
!>      03.2006 initial create
!> \author Sebastian Seidenath (sebastian.seidenath@uni-jena.de)
! **************************************************************************************************
   SUBROUTINE ipi_init(ipi_env, root_section, para_env, force_env_section, &
                       subsys_section)
      TYPE(ipi_environment_type), POINTER                :: ipi_env
      TYPE(section_vals_type), POINTER                   :: root_section
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(section_vals_type), POINTER                   :: force_env_section, subsys_section

      CHARACTER(len=*), PARAMETER                        :: routineN = 'ipi_init'

      INTEGER                                            :: handle
      TYPE(cp_subsys_type), POINTER                      :: subsys
      TYPE(section_vals_type), POINTER                   :: driver_section, motion_section

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(ipi_env))

      ! nullifying  pointers
      NULLIFY (subsys)

      IF (.NOT. ASSOCIATED(subsys_section)) THEN
         subsys_section => section_vals_get_subs_vals(force_env_section, "SUBSYS")
      END IF

      CALL ipi_env_set(ipi_env=ipi_env, force_env_input=force_env_section)

      CALL cp_subsys_create(subsys, para_env, root_section)

      CALL ipi_init_subsys(ipi_env=ipi_env, subsys=subsys, subsys_section=subsys_section)

      motion_section => section_vals_get_subs_vals(root_section, "MOTION")
      driver_section => section_vals_get_subs_vals(motion_section, "DRIVER")
      CALL start_server(para_env=para_env, driver_section=driver_section, ipi_env=ipi_env)

      CALL timestop(handle)

   END SUBROUTINE ipi_init

! **************************************************************************************************
!> \brief Initialize the ipi environment
!> \param ipi_env The ipi environment
!> \param subsys the subsys
!> \param subsys_section ...
!> \par History
!>      03.2024 initial create
!> \author Sebastian Seidenath (sebastian.seidenath@uni-jena.de)
! **************************************************************************************************
   SUBROUTINE ipi_init_subsys(ipi_env, subsys, subsys_section)
      TYPE(ipi_environment_type), POINTER                :: ipi_env
      TYPE(cp_subsys_type), POINTER                      :: subsys
      TYPE(section_vals_type), POINTER                   :: subsys_section

      CHARACTER(len=*), PARAMETER                        :: routineN = 'ipi_init_subsys'

      INTEGER                                            :: handle, natom
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: my_cell, my_cell_ref
      TYPE(distribution_1d_type), POINTER                :: local_molecules, local_particles
      TYPE(molecule_kind_type), DIMENSION(:), POINTER    :: molecule_kind_set
      TYPE(molecule_type), DIMENSION(:), POINTER         :: molecule_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)

      NULLIFY (atomic_kind_set, molecule_kind_set, particle_set, molecule_set, &
               local_molecules, local_particles, my_cell, my_cell_ref)

      particle_set => subsys%particles%els
      atomic_kind_set => subsys%atomic_kinds%els
      molecule_kind_set => subsys%molecule_kinds%els
      molecule_set => subsys%molecules%els
      my_cell => subsys%cell
      my_cell_ref => subsys%cell_ref

!   *** Print the molecule kind set ***
      CALL write_molecule_kind_set(molecule_kind_set, subsys_section)

!   *** Print the atomic coordinates ***
      CALL write_fist_particle_coordinates(particle_set, subsys_section)
      CALL write_particle_distances(particle_set, cell=my_cell, &
                                    subsys_section=subsys_section)
      CALL write_structure_data(particle_set, cell=my_cell, &
                                input_section=subsys_section)

!   *** Distribute molecules and atoms using the new data structures ***
      CALL distribute_molecules_1d(atomic_kind_set=atomic_kind_set, &
                                   particle_set=particle_set, &
                                   local_particles=local_particles, &
                                   molecule_kind_set=molecule_kind_set, &
                                   molecule_set=molecule_set, &
                                   local_molecules=local_molecules, &
                                   force_env_section=ipi_env%force_env_input)

      natom = SIZE(particle_set)

      ALLOCATE (ipi_env%ipi_forces(3, natom))
      ipi_env%ipi_forces(:, :) = 0.0_dp

      CALL ipi_env_set(ipi_env=ipi_env, subsys=subsys, &
                       cell_ref=my_cell_ref, &
                       local_molecules=local_molecules, &
                       local_particles=local_particles)

      CALL distribution_1d_release(local_particles)
      CALL distribution_1d_release(local_molecules)

      CALL timestop(handle)

   END SUBROUTINE ipi_init_subsys

END MODULE ipi_environment
