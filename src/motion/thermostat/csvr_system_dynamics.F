!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \author 10.2007 [tlaino] - Teodoro Laino - University of Zurich
! **************************************************************************************************
MODULE csvr_system_dynamics

   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE csvr_system_types,               ONLY: csvr_system_type
   USE csvr_system_utils,               ONLY: rescaling_factor
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE extended_system_types,           ONLY: map_info_type,&
                                              npt_info_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_comm_type
   USE molecule_kind_types,             ONLY: molecule_kind_type
   USE molecule_types,                  ONLY: molecule_type
   USE parallel_rng_types,              ONLY: rng_stream_type
   USE particle_types,                  ONLY: particle_type
   USE thermostat_utils,                ONLY: ke_region_baro,&
                                              ke_region_particles,&
                                              ke_region_shells,&
                                              vel_rescale_baro,&
                                              vel_rescale_particles,&
                                              vel_rescale_shells
#include "../../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   LOGICAL, PARAMETER :: debug_this_module = .FALSE.
   PUBLIC :: csvr_particles, &
             csvr_barostat, &
             csvr_shells

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'csvr_system_dynamics'

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param csvr ...
!> \param npt ...
!> \param group ...
!> \author 10.2007 [tlaino] - Teodoro Laino - University of Zurich
! **************************************************************************************************
   SUBROUTINE csvr_barostat(csvr, npt, group)

      TYPE(csvr_system_type), POINTER                    :: csvr
      TYPE(npt_info_type), DIMENSION(:, :), &
         INTENT(INOUT)                                   :: npt
      TYPE(mp_comm_type), INTENT(IN)                     :: group

      CHARACTER(len=*), PARAMETER                        :: routineN = 'csvr_barostat'

      INTEGER                                            :: handle
      TYPE(map_info_type), POINTER                       :: map_info

      CALL timeset(routineN, handle)
      map_info => csvr%map_info

      ! Compute the kinetic energy of the barostat
      CALL ke_region_baro(map_info, npt, group)

      ! Apply the Canonical Sampling through Velocity Rescaling
      CALL do_csvr(csvr, map_info)

      ! Now scale the particle velocities
      CALL vel_rescale_baro(map_info, npt)

      ! Re-Compute the kinetic energy of the barostat
      CALL ke_region_baro(map_info, npt, group)

      ! Compute thermostat energy
      CALL do_csvr_eval_energy(csvr, map_info)

      CALL timestop(handle)
   END SUBROUTINE csvr_barostat

! **************************************************************************************************
!> \brief ...
!> \param csvr ...
!> \param molecule_kind_set ...
!> \param molecule_set ...
!> \param particle_set ...
!> \param local_molecules ...
!> \param group ...
!> \param shell_adiabatic ...
!> \param shell_particle_set ...
!> \param core_particle_set ...
!> \param vel ...
!> \param shell_vel ...
!> \param core_vel ...
!> \author 10.2007 [tlaino] - Teodoro Laino - University of Zurich
! **************************************************************************************************
   SUBROUTINE csvr_particles(csvr, molecule_kind_set, molecule_set, &
                             particle_set, local_molecules, group, shell_adiabatic, &
                             shell_particle_set, core_particle_set, vel, shell_vel, core_vel)

      TYPE(csvr_system_type), POINTER                    :: csvr
      TYPE(molecule_kind_type), POINTER                  :: molecule_kind_set(:)
      TYPE(molecule_type), POINTER                       :: molecule_set(:)
      TYPE(particle_type), POINTER                       :: particle_set(:)
      TYPE(distribution_1d_type), POINTER                :: local_molecules
      TYPE(mp_comm_type), INTENT(IN)                     :: group
      LOGICAL, INTENT(IN), OPTIONAL                      :: shell_adiabatic
      TYPE(particle_type), OPTIONAL, POINTER             :: shell_particle_set(:), &
                                                            core_particle_set(:)
      REAL(KIND=dp), INTENT(INOUT), OPTIONAL             :: vel(:, :), shell_vel(:, :), &
                                                            core_vel(:, :)

      CHARACTER(len=*), PARAMETER                        :: routineN = 'csvr_particles'

      INTEGER                                            :: handle
      LOGICAL                                            :: my_shell_adiabatic
      TYPE(map_info_type), POINTER                       :: map_info

      CALL timeset(routineN, handle)
      my_shell_adiabatic = .FALSE.
      IF (PRESENT(shell_adiabatic)) my_shell_adiabatic = shell_adiabatic
      map_info => csvr%map_info

      ! Compute the kinetic energy for the region to thermostat
      CALL ke_region_particles(map_info, particle_set, molecule_kind_set, &
                               local_molecules, molecule_set, group, vel)

      ! Apply the Canonical Sampling through Velocity Rescaling
      CALL do_csvr(csvr, map_info)

      ! Now scale the particle velocities
      CALL vel_rescale_particles(map_info, molecule_kind_set, molecule_set, particle_set, &
                                 local_molecules, my_shell_adiabatic, shell_particle_set, core_particle_set, &
                                 vel, shell_vel, core_vel)

      ! Re-Compute the kinetic energy for the region to thermostat
      CALL ke_region_particles(map_info, particle_set, molecule_kind_set, &
                               local_molecules, molecule_set, group, vel)

      ! Compute thermostat energy
      CALL do_csvr_eval_energy(csvr, map_info)

      CALL timestop(handle)
   END SUBROUTINE csvr_particles

! **************************************************************************************************
!> \brief ...
!> \param csvr ...
!> \param atomic_kind_set ...
!> \param particle_set ...
!> \param local_particles ...
!> \param group ...
!> \param shell_particle_set ...
!> \param core_particle_set ...
!> \param vel ...
!> \param shell_vel ...
!> \param core_vel ...
!> \author 10.2007 [tlaino] - Teodoro Laino - University of Zurich
! **************************************************************************************************
   SUBROUTINE csvr_shells(csvr, atomic_kind_set, particle_set, local_particles, &
                          group, shell_particle_set, core_particle_set, vel, shell_vel, core_vel)

      TYPE(csvr_system_type), POINTER                    :: csvr
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind_set(:)
      TYPE(particle_type), POINTER                       :: particle_set(:)
      TYPE(distribution_1d_type), POINTER                :: local_particles
      TYPE(mp_comm_type), INTENT(IN)                     :: group
      TYPE(particle_type), OPTIONAL, POINTER             :: shell_particle_set(:), &
                                                            core_particle_set(:)
      REAL(KIND=dp), INTENT(INOUT), OPTIONAL             :: vel(:, :), shell_vel(:, :), &
                                                            core_vel(:, :)

      CHARACTER(len=*), PARAMETER                        :: routineN = 'csvr_shells'

      INTEGER                                            :: handle
      TYPE(map_info_type), POINTER                       :: map_info

      CALL timeset(routineN, handle)
      map_info => csvr%map_info

      ! Compute the kinetic energy of the region to thermostat
      CALL ke_region_shells(map_info, particle_set, atomic_kind_set, local_particles, &
                            group, core_particle_set, shell_particle_set, core_vel, shell_vel)

      ! Apply the Canonical Sampling through Velocity Rescaling
      CALL do_csvr(csvr, map_info)

      ! Now scale the particle velocities
      CALL vel_rescale_shells(map_info, atomic_kind_set, particle_set, local_particles, &
                              shell_particle_set, core_particle_set, shell_vel, core_vel, vel)

      ! Re-Compute the kinetic energy of the region to thermostat
      CALL ke_region_shells(map_info, particle_set, atomic_kind_set, local_particles, &
                            group, core_particle_set, shell_particle_set, core_vel, shell_vel)

      ! Compute thermostat energy
      CALL do_csvr_eval_energy(csvr, map_info)

      CALL timestop(handle)
   END SUBROUTINE csvr_shells

! **************************************************************************************************
!> \brief ...
!> \param csvr ...
!> \param map_info ...
!> \author 10.2007 [tlaino] - Teodoro Laino - University of Zurich
! **************************************************************************************************
   SUBROUTINE do_csvr(csvr, map_info)
      TYPE(csvr_system_type), POINTER                    :: csvr
      TYPE(map_info_type), POINTER                       :: map_info

      INTEGER                                            :: i, imap, ndeg
      REAL(KIND=dp)                                      :: kin_energy, kin_target, taut
      TYPE(rng_stream_type), POINTER                     :: rng_stream

      DO i = 1, csvr%loc_num_csvr
         imap = map_info%map_index(i)
         kin_energy = map_info%s_kin(imap)
         csvr%nvt(i)%region_kin_energy = kin_energy
         kin_energy = kin_energy*0.5_dp
         kin_target = csvr%nvt(i)%nkt*0.5_dp
         ndeg = csvr%nvt(i)%degrees_of_freedom
         taut = csvr%tau_csvr/csvr%dt_fact
         rng_stream => csvr%nvt(i)%gaussian_rng_stream
         map_info%v_scale(imap) = rescaling_factor(kin_energy, kin_target, ndeg, taut, &
                                                   rng_stream)
      END DO

   END SUBROUTINE do_csvr

! **************************************************************************************************
!> \brief ...
!> \param csvr ...
!> \param map_info ...
!> \author 10.2007 [tlaino] - Teodoro Laino - University of Zurich
! **************************************************************************************************
   SUBROUTINE do_csvr_eval_energy(csvr, map_info)
      TYPE(csvr_system_type), POINTER                    :: csvr
      TYPE(map_info_type), POINTER                       :: map_info

      INTEGER                                            :: i, imap
      REAL(KIND=dp)                                      :: kin_energy_ar, kin_energy_br

      DO i = 1, csvr%loc_num_csvr
         imap = map_info%map_index(i)
         kin_energy_br = csvr%nvt(i)%region_kin_energy
         kin_energy_ar = map_info%s_kin(imap)
         csvr%nvt(i)%thermostat_energy = csvr%nvt(i)%thermostat_energy + &
                                         0.5_dp*(kin_energy_br - kin_energy_ar)
      END DO

   END SUBROUTINE do_csvr_eval_energy

END MODULE csvr_system_dynamics
