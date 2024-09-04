!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief The environment for the empirical interatomic potential methods.
!> \par History
!>      03.2024 initial create
!> \author Sebastian Seidenath (sebastian.seidenath@uni-jena.de)
! **************************************************************************************************
MODULE ipi_environment_types
   USE atomic_kind_list_types,          ONLY: atomic_kind_list_create,&
                                              atomic_kind_list_release,&
                                              atomic_kind_list_type
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cell_types,                      ONLY: cell_release,&
                                              cell_retain,&
                                              cell_type
   USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                              cp_subsys_release,&
                                              cp_subsys_set,&
                                              cp_subsys_type
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE input_section_types,             ONLY: section_vals_release,&
                                              section_vals_retain,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE molecule_kind_list_types,        ONLY: molecule_kind_list_create,&
                                              molecule_kind_list_release,&
                                              molecule_kind_list_type
   USE molecule_kind_types,             ONLY: molecule_kind_type
   USE molecule_list_types,             ONLY: molecule_list_create,&
                                              molecule_list_release,&
                                              molecule_list_type
   USE molecule_types,                  ONLY: molecule_type
   USE particle_list_types,             ONLY: particle_list_create,&
                                              particle_list_release,&
                                              particle_list_type
   USE particle_types,                  ONLY: particle_type
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ipi_environment_types'

   ! *** Public data types ***
   PUBLIC :: ipi_environment_type

   ! *** Public subroutines ***
   PUBLIC :: ipi_env_release, &
             ipi_env_set, &
             ipi_env_get, &
             ipi_env_create

! **************************************************************************************************
!> \brief The i–PI environment
!> \param ipi_energy The total ipi energy
!> \param ipi_forces The final ipi forces [eV/A]
!> \param subsystem The particles, molecules,... of this environment
!> \param force_env_input Pointer to the force_env input section
!> \param cell_ref The reference simulation cell
!> \par History
!>      03.2024 initial create
!> \author Sebastian Seidenath (sebastian.seidenath@uni-jena.de)
! **************************************************************************************************
   TYPE ipi_environment_type
      REAL(KIND=dp)                                    :: ipi_energy = 0.0_dp
      REAL(KIND=dp), DIMENSION(:, :), POINTER          :: ipi_forces => Null()
      TYPE(cp_subsys_type), POINTER                    :: subsys => Null()
      TYPE(section_vals_type), POINTER                 :: force_env_input => Null()
      TYPE(cell_type), POINTER                         :: cell_ref => Null()
      INTEGER                                          :: sockfd = -1
   END TYPE ipi_environment_type

CONTAINS

! **************************************************************************************************
!> \brief Releases the given ipi environment (see doc/ReferenceCounting.html)
!> \param ipi_env The ipi environment to release
!> \par History
!>      03.2024 initial create
!> \author Sebastian Seidenath (sebastian.seidenath@uni-jena.de)
! **************************************************************************************************
   SUBROUTINE ipi_env_release(ipi_env)
      TYPE(ipi_environment_type), INTENT(INOUT)          :: ipi_env

      IF (ASSOCIATED(ipi_env%ipi_forces)) THEN
         DEALLOCATE (ipi_env%ipi_forces)
      END IF
      IF (ASSOCIATED(ipi_env%subsys)) THEN
         CALL cp_subsys_release(ipi_env%subsys)
      END IF
      IF (ASSOCIATED(ipi_env%force_env_input)) THEN
         CALL section_vals_release(ipi_env%force_env_input)
      END IF
      IF (ASSOCIATED(ipi_env%cell_ref)) THEN
         CALL cell_release(ipi_env%cell_ref)
      END IF
   END SUBROUTINE ipi_env_release

! **************************************************************************************************
!> \brief Returns various attributes of the ipi environment
!> \param ipi_env The enquired ipi environment
!> \param ipi_energy The total ipi energy
!> \param ipi_forces The final ipi forces [eV/A]
!> \param subsys the particles, molecules,... of this environment
!> \param atomic_kind_set The set of all atomic kinds involved
!> \param particle_set The set of all particles
!> \param local_particles All particles on this particular node
!> \param molecule_kind_set The set of all different molecule kinds involved
!> \param molecule_set The set of all molecules
!> \param local_molecules All molecules on this particular node
!> \param force_env_input Pointer to the force_env input section
!> \param cell The simulation cell
!> \param cell_ref The reference simulation cell
!> \param virial Dummy virial pointer
!> \param sockfd File descriptor of the communications socket
!> \par History
!>      03.2024 initial create
!> \author Sebastian Seidenath (sebastian.seidenath@uni-jena.de)
! **************************************************************************************************
   SUBROUTINE ipi_env_get(ipi_env, ipi_energy, ipi_forces, subsys, &
                          atomic_kind_set, particle_set, local_particles, &
                          molecule_kind_set, molecule_set, local_molecules, &
                          force_env_input, cell, cell_ref, virial, sockfd)

      TYPE(ipi_environment_type), INTENT(IN)             :: ipi_env
      REAL(kind=dp), OPTIONAL                            :: ipi_energy
      REAL(KIND=dp), DIMENSION(:, :), OPTIONAL, POINTER  :: ipi_forces
      TYPE(cp_subsys_type), OPTIONAL, POINTER            :: subsys
      TYPE(atomic_kind_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: atomic_kind_set
      TYPE(particle_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: particle_set
      TYPE(distribution_1d_type), OPTIONAL, POINTER      :: local_particles
      TYPE(molecule_kind_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: molecule_kind_set
      TYPE(molecule_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: molecule_set
      TYPE(distribution_1d_type), OPTIONAL, POINTER      :: local_molecules
      TYPE(section_vals_type), OPTIONAL, POINTER         :: force_env_input
      TYPE(cell_type), OPTIONAL, POINTER                 :: cell, cell_ref
      TYPE(virial_type), OPTIONAL, POINTER               :: virial
      INTEGER, OPTIONAL                                  :: sockfd

      TYPE(atomic_kind_list_type), POINTER               :: atomic_kinds
      TYPE(molecule_kind_list_type), POINTER             :: molecule_kinds
      TYPE(molecule_list_type), POINTER                  :: molecules
      TYPE(particle_list_type), POINTER                  :: particles

      NULLIFY (atomic_kinds, particles, molecules, molecule_kinds)

      IF (PRESENT(ipi_energy)) ipi_energy = ipi_env%ipi_energy
      IF (PRESENT(ipi_forces)) ipi_forces = ipi_env%ipi_forces
      IF (PRESENT(subsys)) subsys => ipi_env%subsys
      CALL cp_subsys_get(ipi_env%subsys, &
                         atomic_kinds=atomic_kinds, &
                         particles=particles, &
                         molecule_kinds=molecule_kinds, &
                         molecules=molecules, &
                         local_molecules=local_molecules, &
                         local_particles=local_particles, &
                         virial=virial, &
                         cell=cell)
      IF (PRESENT(atomic_kind_set)) atomic_kind_set => atomic_kinds%els
      IF (PRESENT(particle_set)) particle_set => particles%els
      IF (PRESENT(molecule_kind_set)) molecule_kind_set => molecule_kinds%els
      IF (PRESENT(molecule_set)) molecule_set => molecules%els

      IF (PRESENT(force_env_input)) force_env_input => ipi_env%force_env_input
      IF (PRESENT(cell_ref)) cell_ref => ipi_env%cell_ref
      IF (PRESENT(sockfd)) sockfd = ipi_env%sockfd

   END SUBROUTINE ipi_env_get

! **************************************************************************************************
!> \brief Sets various attributes of the ipi environment
!> \param ipi_env The enquired ipi environment
!> \param ipi_energy The total ipi energy
!> \param ipi_forces The final ipi forces [eV/A]
!> \param subsys the particles, molecules,... of this environment
!> \param atomic_kind_set The set of all atomic kinds involved
!> \param particle_set The set of all particles
!> \param local_particles All particles on this particular node
!> \param molecule_kind_set The set of all different molecule kinds involved
!> \param molecule_set The set of all molecules
!> \param local_molecules All molecules on this particular node
!> \param force_env_input Pointer to the force_env input section
!> \param cell_ref The reference simulation cell
!> \param sockfd File descriptor of the communications socket
!> \par History
!>      03.2024 initial create
!> \author Sebastian Seidenath (sebastian.seidenath@uni-jena.de)
! **************************************************************************************************
   SUBROUTINE ipi_env_set(ipi_env, ipi_energy, ipi_forces, subsys, &
                          atomic_kind_set, particle_set, local_particles, &
                          molecule_kind_set, molecule_set, local_molecules, &
                          force_env_input, cell_ref, sockfd)

      TYPE(ipi_environment_type), INTENT(INOUT)          :: ipi_env
      REAL(KIND=dp), OPTIONAL                            :: ipi_energy
      REAL(KIND=dp), DIMENSION(:, :), OPTIONAL, POINTER  :: ipi_forces
      TYPE(cp_subsys_type), OPTIONAL, POINTER            :: subsys
      TYPE(atomic_kind_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: atomic_kind_set
      TYPE(particle_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: particle_set
      TYPE(distribution_1d_type), OPTIONAL, POINTER      :: local_particles
      TYPE(molecule_kind_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: molecule_kind_set
      TYPE(molecule_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: molecule_set
      TYPE(distribution_1d_type), OPTIONAL, POINTER      :: local_molecules
      TYPE(section_vals_type), OPTIONAL, POINTER         :: force_env_input
      TYPE(cell_type), OPTIONAL, POINTER                 :: cell_ref
      INTEGER, OPTIONAL                                  :: sockfd

      TYPE(atomic_kind_list_type), POINTER               :: atomic_kinds
      TYPE(molecule_kind_list_type), POINTER             :: molecule_kinds
      TYPE(molecule_list_type), POINTER                  :: molecules
      TYPE(particle_list_type), POINTER                  :: particles

      IF (PRESENT(ipi_energy)) ipi_env%ipi_energy = ipi_energy
      IF (PRESENT(ipi_forces)) ipi_env%ipi_forces = ipi_forces
      IF (PRESENT(subsys)) THEN
         IF (ASSOCIATED(ipi_env%subsys)) THEN
            IF (.NOT. ASSOCIATED(ipi_env%subsys, subsys)) THEN
               CALL cp_subsys_release(ipi_env%subsys)
            END IF
         END IF
         ipi_env%subsys => subsys
      END IF
      IF (PRESENT(atomic_kind_set)) THEN
         CALL atomic_kind_list_create(atomic_kinds, els_ptr=atomic_kind_set)
         CALL cp_subsys_set(ipi_env%subsys, atomic_kinds=atomic_kinds)
         CALL atomic_kind_list_release(atomic_kinds)
      END IF
      IF (PRESENT(particle_set)) THEN
         CALL particle_list_create(particles, els_ptr=particle_set)
         CALL cp_subsys_set(ipi_env%subsys, particles=particles)
         CALL particle_list_release(particles)
      END IF
      IF (PRESENT(molecule_kind_set)) THEN
         CALL molecule_kind_list_create(molecule_kinds, els_ptr=molecule_kind_set)
         CALL cp_subsys_set(ipi_env%subsys, molecule_kinds=molecule_kinds)
         CALL molecule_kind_list_release(molecule_kinds)
      END IF
      IF (PRESENT(molecule_set)) THEN
         CALL molecule_list_create(molecules, els_ptr=molecule_set)
         CALL cp_subsys_set(ipi_env%subsys, molecules=molecules)
         CALL molecule_list_release(molecules)
      END IF
      IF (PRESENT(local_particles)) THEN
         CALL cp_subsys_set(ipi_env%subsys, local_particles=local_particles)
      END IF
      IF (PRESENT(local_molecules)) THEN
         CALL cp_subsys_set(ipi_env%subsys, local_molecules=local_molecules)
      END IF

      IF (PRESENT(force_env_input)) THEN
         CALL section_vals_retain(force_env_input)
         CALL section_vals_release(ipi_env%force_env_input)
         ipi_env%force_env_input => force_env_input
      END IF
      IF (PRESENT(cell_ref)) THEN
         CALL cell_retain(cell_ref)
         CALL cell_release(ipi_env%cell_ref)
         ipi_env%cell_ref => cell_ref
      END IF
      IF (PRESENT(sockfd)) ipi_env%sockfd = sockfd
   END SUBROUTINE ipi_env_set

! **************************************************************************************************
!> \brief Reinitializes the ipi environment
!> \param ipi_env The ipi environment to be reinitialized
!> \par History
!>      03.2024 initial create
!> \author Sebastian Seidenath (sebastian.seidenath@uni-jena.de)
! **************************************************************************************************
   SUBROUTINE ipi_env_clear(ipi_env)

      TYPE(ipi_environment_type), INTENT(INOUT)          :: ipi_env

      IF (ASSOCIATED(ipi_env%ipi_forces)) THEN
         ipi_env%ipi_forces(:, :) = 0.0_dp
      END IF
      IF (ASSOCIATED(ipi_env%subsys)) THEN
         CALL cp_subsys_release(ipi_env%subsys)
      END IF
      IF (ASSOCIATED(ipi_env%force_env_input)) THEN
         CALL section_vals_release(ipi_env%force_env_input)
      END IF
      IF (ASSOCIATED(ipi_env%cell_ref)) THEN
         CALL cell_release(ipi_env%cell_ref)
      END IF
   END SUBROUTINE ipi_env_clear

! **************************************************************************************************
!> \brief Creates the ipi environment
!> \param ipi_env The ipi environment to be created
!> \par History
!>      03.2024 initial create
!> \author Sebastian Seidenath (sebastian.seidenath@uni-jena.de)
! **************************************************************************************************
   SUBROUTINE ipi_env_create(ipi_env)

      TYPE(ipi_environment_type), INTENT(OUT)            :: ipi_env

      NULLIFY (ipi_env%ipi_forces)
      NULLIFY (ipi_env%subsys)
      NULLIFY (ipi_env%force_env_input)
      NULLIFY (ipi_env%cell_ref)

      ipi_env%ipi_energy = 0_dp
      ipi_env%sockfd = 0 ! stdinp
      CALL ipi_env_clear(ipi_env)
   END SUBROUTINE ipi_env_create

END MODULE ipi_environment_types
