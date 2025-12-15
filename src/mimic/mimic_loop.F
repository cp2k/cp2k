!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Module containing the main loop for simulations with the MiMiC framework
!> \par History
!>      05.2025 Created [AA]
!> \author Andrej Antalik
! **************************************************************************************************

MODULE mimic_loop

   USE force_env_methods,               ONLY: force_env_calc_energy_force
   USE force_env_types,                 ONLY: force_env_type
   USE mcl_requests,                    ONLY: &
        MCL_COMPUTE_FORCES, MCL_EXIT, MCL_RECV_PARTICLE_POSITIONS, &
        MCL_RECV_POTENTIAL_ON_GRIDPOINTS, MCL_SEND_BOX_GRIDPOINT_COORDS, &
        MCL_SEND_BOX_NUM_GRIDPOINTS, MCL_SEND_BOX_ORIGIN, MCL_SEND_BOX_VECTORS, &
        MCL_SEND_CLIENT_APIVER, MCL_SEND_CLIENT_ID, MCL_SEND_CLIENT_NAME, MCL_SEND_CLIENT_RUNTYPE, &
        MCL_SEND_DENSITY, MCL_SEND_ENERGY, MCL_SEND_NUCLEAR_CHARGES, MCL_SEND_NUM_CONSTR_ANGLES, &
        MCL_SEND_NUM_CONSTR_BONDS, MCL_SEND_NUM_FRAGMENTS, MCL_SEND_NUM_PARTICLES, &
        MCL_SEND_NUM_PARTICLES_IN_FRAGMENTS, MCL_SEND_NUM_PARTICLE_SPECIES, &
        MCL_SEND_PARTICLE_FORCES, MCL_SEND_PARTICLE_IDS_IN_FRAGMENTS, MCL_SEND_PARTICLE_POSITIONS, &
        MCL_SEND_PARTICLE_SPECIES_IDS, MCL_SEND_REFERENCE_CHARGES, MCL_SEND_SPECIES_ELEMENTS, &
        MCL_SEND_SPECIES_LABELS, MCL_SEND_SPECIES_MASSES
   USE mimic_communicator,              ONLY: mimic_communicator_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: do_mimic_loop

CONTAINS

! **************************************************************************************************
!> \brief The main loop for a MiMiC run
!> \param force_env ...
!> \par History
!>      5.2025 created [AA]
! **************************************************************************************************
   SUBROUTINE do_mimic_loop(force_env)
      TYPE(force_env_type), POINTER                      :: force_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'do_mimic_loop'

      CHARACTER(LEN=32)                                  :: request_str
      INTEGER                                            :: handle, request
      LOGICAL                                            :: is_last_step
      TYPE(mimic_communicator_type)                      :: mimic_comm

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(force_env))

      CALL mimic_comm%initialize(force_env)

      is_last_step = .FALSE.
      DO WHILE (.NOT. is_last_step)
         request = mimic_comm%receive_request()

         SELECT CASE (request)
         CASE (MCL_SEND_CLIENT_ID)
            CALL mimic_comm%send_client_info("id")
         CASE (MCL_SEND_CLIENT_NAME)
            CALL mimic_comm%send_client_info("name")
         CASE (MCL_SEND_CLIENT_RUNTYPE)
            CALL mimic_comm%send_client_info("run_type")
         CASE (MCL_SEND_CLIENT_APIVER)
            CALL mimic_comm%send_client_info("api_version")
         CASE (MCL_SEND_NUM_PARTICLES)
            CALL mimic_comm%send_value("num_atoms")
         CASE (MCL_SEND_NUM_PARTICLE_SPECIES)
            CALL mimic_comm%send_value("num_kinds")
         CASE (MCL_SEND_PARTICLE_SPECIES_IDS)
            CALL mimic_comm%send_atom_info("kinds")
         CASE (MCL_SEND_NUM_FRAGMENTS)
            CALL mimic_comm%send_value("num_fragments")
         CASE (MCL_SEND_NUM_PARTICLES_IN_FRAGMENTS)
            CALL mimic_comm%send_value("num_atoms_in_fragments")
         CASE (MCL_SEND_NUM_CONSTR_BONDS)
            CALL mimic_comm%send_value("num_bonds")
         CASE (MCL_SEND_NUM_CONSTR_ANGLES)
            CALL mimic_comm%send_value("num_angles")
         CASE (MCL_SEND_SPECIES_MASSES)
            CALL mimic_comm%send_kind_info("masses")
         CASE (MCL_SEND_PARTICLE_IDS_IN_FRAGMENTS)
            CALL mimic_comm%send_atom_info("ids")
         CASE (MCL_SEND_SPECIES_ELEMENTS)
            CALL mimic_comm%send_kind_info("elements")
         CASE (MCL_SEND_SPECIES_LABELS)
            CALL mimic_comm%send_kind_info("labels")
         CASE (MCL_SEND_PARTICLE_POSITIONS)
            CALL mimic_comm%send_positions()
         CASE (MCL_RECV_PARTICLE_POSITIONS)
            CALL mimic_comm%receive_positions()
         CASE (MCL_SEND_PARTICLE_FORCES)
            CALL mimic_comm%send_forces()
         CASE (MCL_SEND_ENERGY)
            CALL mimic_comm%send_value("energy")
         CASE (MCL_SEND_NUCLEAR_CHARGES)
            CALL mimic_comm%send_kind_info("nuclear_charges")
         CASE (MCL_SEND_BOX_NUM_GRIDPOINTS)
            CALL mimic_comm%send_box_info("num_gridpoints")
         CASE (MCL_SEND_BOX_ORIGIN)
            CALL mimic_comm%send_box_info("origin")
         CASE (MCL_SEND_BOX_VECTORS)
            CALL mimic_comm%send_box_info("box_vectors")
         CASE (MCL_SEND_BOX_GRIDPOINT_COORDS)
            CALL mimic_comm%send_grid_coordinates()
         CASE (MCL_SEND_DENSITY)
            CALL mimic_comm%send_density()
         CASE (MCL_RECV_POTENTIAL_ON_GRIDPOINTS)
            CALL mimic_comm%receive_potential()
         CASE (MCL_SEND_REFERENCE_CHARGES)
            CALL mimic_comm%send_result("hirshfeld_charges")
         CASE (MCL_COMPUTE_FORCES)
            CALL force_env_calc_energy_force(force_env, calc_force=.TRUE.)
         CASE (MCL_EXIT)
            is_last_step = .TRUE.
         CASE DEFAULT
            WRITE (request_str, "(I0)") request
            CPABORT("Unrecognized MiMiC request: "//TRIM(request_str))
         END SELECT
      END DO

      CALL mimic_comm%finalize()

      CALL timestop(handle)

   END SUBROUTINE do_mimic_loop

END MODULE mimic_loop

