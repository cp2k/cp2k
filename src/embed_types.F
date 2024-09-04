!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE embed_types
   USE atomic_kind_list_types,          ONLY: atomic_kind_list_create,&
                                              atomic_kind_list_release,&
                                              atomic_kind_list_type
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cell_types,                      ONLY: cell_release,&
                                              cell_retain,&
                                              cell_type
   USE cp_fm_types,                     ONLY: cp_fm_type
   USE cp_log_handling,                 ONLY: cp_logger_p_type,&
                                              cp_logger_release
   USE cp_result_types,                 ONLY: cp_result_type
   USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                              cp_subsys_release,&
                                              cp_subsys_set,&
                                              cp_subsys_type
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE input_section_types,             ONLY: section_vals_release,&
                                              section_vals_retain,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE lri_environment_types,           ONLY: lri_kind_type
   USE message_passing,                 ONLY: mp_para_env_p_type,&
                                              mp_para_env_release,&
                                              mp_para_env_type
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
   USE pw_types,                        ONLY: pw_r3d_rs_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: opt_embed_pot_type, opt_dmfet_pot_type, embed_env_type

! *** Public subroutines ***

   PUBLIC :: get_embed_env, &
             set_embed_env, &
             embed_env_create, &
             embed_env_release

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'embed_types'

! **************************************************************************************************
!> \brief Type containing main data for embedding potential optimization
!> \author Vladimir Rybkin
! **************************************************************************************************

   TYPE opt_embed_pot_type
      TYPE(cp_fm_type), POINTER                 :: embed_pot_coef => NULL(), embed_pot_grad => NULL(), &
                                                   prev_step => NULL(), step => NULL(), embed_pot_hess => NULL(), &
                                                   prev_embed_pot_coef => NULL(), prev_embed_pot_grad => NULL(), &
                                                   prev_embed_pot_hess => NULL(), kinetic_mat => NULL()
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)  :: w_func, max_diff, int_diff, int_diff_square, &
                                                   max_grid_step, max_subsys_dens_diff
      INTEGER                                   :: n_iter = -1, i_iter = -1, dimen_aux = -1, &
                                                   last_accepted = -1, dimen_var_aux = -1
      REAL(KIND=dp)                             :: lambda = 0.0_dp, allowed_decrease = 0.0_dp, &
                                                   max_trad = 0.0_dp, min_trad = 0.0_dp, &
                                                   grad_norm = 0.0_dp, vw_cutoff = 0.0_dp, &
                                                   vw_smooth_cutoff_range = 0.0_dp, eta = 0.0_dp
      TYPE(pw_r3d_rs_type), POINTER                  :: const_pot => NULL(), prev_embed_pot => NULL(), &
                                                        prev_spin_embed_pot => NULL(), pot_diff => NULL()
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER    :: prev_grid_grad => NULL(), prev_grid_pot => NULL(), &
                                                        prev_subsys_dens => NULL(), v_w => NULL()
      REAL(KIND=dp)                             :: reg_term = 0.0_dp, trust_rad = 0.0_dp, &
                                                   conv_max = 0.0_dp, conv_int = 0.0_dp, &
                                                   conv_max_spin = 0.0_dp, conv_int_spin = 0.0_dp, step_len = 0.0_dp
      LOGICAL                                   :: accept_step = .FALSE., newton_step = .FALSE., &
                                                   level_shift = .FALSE., steep_desc = .FALSE., &
                                                   add_const_pot = .FALSE., converged = .FALSE., read_embed_pot = .FALSE., &
                                                   read_embed_pot_cube = .FALSE., change_spin = .FALSE., &
                                                   open_shell_embed = .FALSE., &
                                                   grid_opt = .FALSE., leeuwen = .FALSE., fab = .FALSE., &
                                                   coulomb_guess = .FALSE., fermi_amaldi = .FALSE., &
                                                   diff_guess = .FALSE.
      INTEGER, ALLOCATABLE, DIMENSION(:)        :: all_nspins
      TYPE(lri_kind_type), DIMENSION(:), POINTER :: lri => NULL() ! Needed to store integrals

   END TYPE opt_embed_pot_type

! **************************************************************************************************
!> \brief Type containing main data for matrix embedding potential optimization
!> \author Vladimir Rybkin
! **************************************************************************************************

   TYPE opt_dmfet_pot_type
      TYPE(cp_fm_type), POINTER                 :: dmfet_pot => NULL(), dm_1 => NULL(), dm_2 => NULL(), &
                                                   dm_total => NULL(), dm_diff => NULL(), &
                                                   dmfet_pot_beta => NULL(), dm_total_beta => NULL(), &
                                                   dm_diff_beta => NULL(), dm_subsys => NULL(), dm_subsys_beta => NULL()
      REAL(KIND=dp)                             :: trust_rad = 0.0_dp, conv_max = 0.0_dp, &
                                                   conv_int = 0.0_dp, conv_max_beta = 0.0_dp, &
                                                   conv_int_beta = 0.0_dp
      LOGICAL                                   :: accept_step = .FALSE., converged = .FALSE., &
                                                   read_dmfet_pot = .FALSE., &
                                                   change_spin = .FALSE., open_shell_embed = .FALSE.
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)  :: w_func, max_diff, int_diff
      INTEGER                                   :: n_iter = -1, i_iter = -1, last_accepted = -1
      INTEGER, ALLOCATABLE, DIMENSION(:)        :: all_nspins

   END TYPE opt_dmfet_pot_type

! **************************************************************************************************
!> \brief Embedding environment type
!> \author Vladimir Rybkin
! **************************************************************************************************

   TYPE embed_env_type
      TYPE(cell_type), POINTER                         :: cell_ref => NULL()
      TYPE(mp_para_env_type), POINTER                  :: para_env => NULL()
      TYPE(cp_subsys_type), POINTER                    :: subsys => NULL()
      TYPE(section_vals_type), POINTER                 :: input => NULL()
      REAL(KIND=dp), DIMENSION(:), POINTER             :: energies => NULL()
      ! Parallelization of multiple force_eval
      INTEGER                                          :: ngroups = -1
      INTEGER, DIMENSION(:), POINTER                   :: group_distribution => NULL()
      TYPE(mp_para_env_p_type), DIMENSION(:), POINTER  :: sub_para_env => NULL()
      TYPE(cp_logger_p_type), DIMENSION(:), POINTER    :: sub_logger => NULL()
      ! Densities from sunbsystem
      REAL(KIND=dp)                                    :: pot_energy = 0.0_dp
   END TYPE embed_env_type

CONTAINS

! **************************************************************************************************
!> \brief Get the embed environment.
!> \param embed_env ...
!> \param atomic_kind_set ...
!> \param particle_set ...
!> \param local_particles ...
!> \param local_molecules ...
!> \param molecule_kind_set ...
!> \param molecule_set ...
!> \param cell ...
!> \param cell_ref ...
!> \param para_env ...
!> \param sub_para_env ...
!> \param subsys ...
!> \param input ...
!> \param results ...
!> \param pot_energy ...
! **************************************************************************************************
   SUBROUTINE get_embed_env(embed_env, atomic_kind_set, particle_set, &
                            local_particles, local_molecules, molecule_kind_set, &
                            molecule_set, cell, cell_ref, &
                            para_env, sub_para_env, subsys, &
                            input, results, pot_energy)

      TYPE(embed_env_type), INTENT(IN)                   :: embed_env
      TYPE(atomic_kind_type), OPTIONAL, POINTER          :: atomic_kind_set(:)
      TYPE(particle_type), OPTIONAL, POINTER             :: particle_set(:)
      TYPE(distribution_1d_type), OPTIONAL, POINTER      :: local_particles, local_molecules
      TYPE(molecule_kind_type), OPTIONAL, POINTER        :: molecule_kind_set(:)
      TYPE(molecule_type), OPTIONAL, POINTER             :: molecule_set(:)
      TYPE(cell_type), OPTIONAL, POINTER                 :: cell, cell_ref
      TYPE(mp_para_env_type), OPTIONAL, POINTER          :: para_env
      TYPE(mp_para_env_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: sub_para_env
      TYPE(cp_subsys_type), OPTIONAL, POINTER            :: subsys
      TYPE(section_vals_type), OPTIONAL, POINTER         :: input
      TYPE(cp_result_type), OPTIONAL, POINTER            :: results
      REAL(KIND=dp), OPTIONAL                            :: pot_energy

      TYPE(atomic_kind_list_type), POINTER               :: atomic_kinds
      TYPE(molecule_kind_list_type), POINTER             :: molecule_kinds
      TYPE(molecule_list_type), POINTER                  :: molecules
      TYPE(particle_list_type), POINTER                  :: particles

      NULLIFY (atomic_kinds, particles, molecules, molecule_kinds)
      CPASSERT(ASSOCIATED(embed_env%subsys))

      IF (PRESENT(input)) input => embed_env%input
      IF (PRESENT(cell_ref)) cell_ref => embed_env%cell_ref
      IF (PRESENT(para_env)) para_env => embed_env%para_env
      IF (PRESENT(sub_para_env)) sub_para_env => embed_env%sub_para_env
      IF (PRESENT(subsys)) subsys => embed_env%subsys
      CALL cp_subsys_get(embed_env%subsys, &
                         atomic_kinds=atomic_kinds, &
                         local_molecules=local_molecules, &
                         local_particles=local_particles, &
                         particles=particles, &
                         molecule_kinds=molecule_kinds, &
                         molecules=molecules, &
                         results=results, &
                         cell=cell)
      IF (PRESENT(atomic_kind_set)) atomic_kind_set => atomic_kinds%els
      IF (PRESENT(particle_set)) particle_set => particles%els
      IF (PRESENT(molecule_kind_set)) molecule_kind_set => molecule_kinds%els
      IF (PRESENT(molecule_set)) molecule_set => molecules%els
      IF (PRESENT(pot_energy)) pot_energy = embed_env%pot_energy

   END SUBROUTINE get_embed_env

! **************************************************************************************************
!> \brief ...
!> \param embed_env ...
!> \param atomic_kind_set ...
!> \param particle_set ...
!> \param local_particles ...
!> \param local_molecules ...
!> \param molecule_kind_set ...
!> \param molecule_set ...
!> \param cell_ref ...
!> \param subsys ...
!> \param input ...
!> \param sub_para_env ...
! **************************************************************************************************
   SUBROUTINE set_embed_env(embed_env, atomic_kind_set, particle_set, &
                            local_particles, local_molecules, molecule_kind_set, &
                            molecule_set, cell_ref, subsys, &
                            input, sub_para_env)

      TYPE(embed_env_type), INTENT(INOUT)                :: embed_env
      TYPE(atomic_kind_type), OPTIONAL, POINTER          :: atomic_kind_set(:)
      TYPE(particle_type), OPTIONAL, POINTER             :: particle_set(:)
      TYPE(distribution_1d_type), OPTIONAL, POINTER      :: local_particles, local_molecules
      TYPE(molecule_kind_type), OPTIONAL, POINTER        :: molecule_kind_set(:)
      TYPE(molecule_type), OPTIONAL, POINTER             :: molecule_set(:)
      TYPE(cell_type), OPTIONAL, POINTER                 :: cell_ref
      TYPE(cp_subsys_type), OPTIONAL, POINTER            :: subsys
      TYPE(section_vals_type), OPTIONAL, POINTER         :: input
      TYPE(mp_para_env_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: sub_para_env

      TYPE(atomic_kind_list_type), POINTER               :: atomic_kinds
      TYPE(molecule_kind_list_type), POINTER             :: molecule_kinds
      TYPE(molecule_list_type), POINTER                  :: molecules
      TYPE(particle_list_type), POINTER                  :: particles

      IF (PRESENT(cell_ref)) THEN
         CALL cell_retain(cell_ref)
         CALL cell_release(embed_env%cell_ref)
         embed_env%cell_ref => cell_ref
      END IF
      IF (PRESENT(input)) THEN
         CALL section_vals_retain(input)
         CALL section_vals_release(embed_env%input)
         embed_env%input => input
      END IF
      IF (PRESENT(subsys)) THEN
         IF (ASSOCIATED(embed_env%subsys)) THEN
            IF (.NOT. ASSOCIATED(embed_env%subsys, subsys)) THEN
               CALL cp_subsys_release(embed_env%subsys)
            END IF
         END IF
         embed_env%subsys => subsys
      END IF
      IF (PRESENT(sub_para_env)) THEN
         embed_env%sub_para_env => sub_para_env
      END IF
      IF (PRESENT(atomic_kind_set)) THEN
         CALL atomic_kind_list_create(atomic_kinds, &
                                      els_ptr=atomic_kind_set)
         CALL cp_subsys_set(embed_env%subsys, &
                            atomic_kinds=atomic_kinds)
         CALL atomic_kind_list_release(atomic_kinds)
      END IF
      IF (PRESENT(particle_set)) THEN
         CALL particle_list_create(particles, &
                                   els_ptr=particle_set)
         CALL cp_subsys_set(embed_env%subsys, &
                            particles=particles)
         CALL particle_list_release(particles)
      END IF
      IF (PRESENT(local_particles)) THEN
         CALL cp_subsys_set(embed_env%subsys, &
                            local_particles=local_particles)
      END IF
      IF (PRESENT(local_molecules)) THEN
         CALL cp_subsys_set(embed_env%subsys, &
                            local_molecules=local_molecules)
      END IF
      IF (PRESENT(molecule_kind_set)) THEN
         CALL molecule_kind_list_create(molecule_kinds, els_ptr=molecule_kind_set)
         CALL cp_subsys_set(embed_env%subsys, molecule_kinds=molecule_kinds)
         CALL molecule_kind_list_release(molecule_kinds)
      END IF
      IF (PRESENT(molecule_set)) THEN
         CALL molecule_list_create(molecules, els_ptr=molecule_set)
         CALL cp_subsys_set(embed_env%subsys, molecules=molecules)
         CALL molecule_list_release(molecules)
      END IF

   END SUBROUTINE set_embed_env

! **************************************************************************************************
!> \brief ...
!> \param embed_env ...
!> \param para_env the parallel environment for the qs_env
!> \author Vladimir Rybkin 02.2018
! **************************************************************************************************
   SUBROUTINE embed_env_create(embed_env, para_env)
      TYPE(embed_env_type), INTENT(OUT)                  :: embed_env
      TYPE(mp_para_env_type), INTENT(IN), TARGET         :: para_env

      embed_env%para_env => para_env
      CALL embed_env%para_env%retain()

   END SUBROUTINE embed_env_create

! **************************************************************************************************
!> \brief ...
!> \param embed_env ...
!> \author Vladimir Rybkin 02.2018
! **************************************************************************************************
   SUBROUTINE embed_env_release(embed_env)
      TYPE(embed_env_type), INTENT(INOUT)                :: embed_env

      INTEGER                                            :: i, ngroups

      ngroups = embed_env%ngroups
      DO i = 1, ngroups
         IF (ASSOCIATED(embed_env%sub_para_env(i)%para_env)) THEN
            CALL cp_logger_release(embed_env%sub_logger(i)%p)
            CALL mp_para_env_release(embed_env%sub_para_env(i)%para_env)
         END IF
      END DO
      DEALLOCATE (embed_env%sub_para_env)
      DEALLOCATE (embed_env%sub_logger)
      DEALLOCATE (embed_env%energies)
      CALL cell_release(embed_env%cell_ref)
      CALL mp_para_env_release(embed_env%para_env)
      CALL cp_subsys_release(embed_env%subsys)
      CALL section_vals_release(embed_env%input)
      IF (ASSOCIATED(embed_env%group_distribution)) THEN
         DEALLOCATE (embed_env%group_distribution)
      END IF

   END SUBROUTINE embed_env_release

END MODULE embed_types

