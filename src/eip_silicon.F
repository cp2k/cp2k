!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Empirical interatomic potentials for Silicon
!> \note
!>      Stefan Goedecker's OpenMP implementation of Bazant's EDIP & Lenosky's
!>      empirical interatomic potentials for Silicon.
!> \par History
!>      03.2006 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! **************************************************************************************************
MODULE eip_silicon
   USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE cell_types,                      ONLY: cell_type,&
                                              get_cell
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                              cp_subsys_type
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE eip_environment_types,           ONLY: eip_env_get,&
                                              eip_environment_type
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: angstrom,&
                                              evolt

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads
#include "./base/base_uses.f90"

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'eip_silicon'

   ! *** Public subroutines ***
   PUBLIC :: eip_bazant, eip_lenosky, eip_stillinger_weber, eip_tersoff

!***

CONTAINS

! **************************************************************************************************
!> \brief Interface routine of Goedecker's Bazant EDIP to CP2K
!> \param eip_env ...
!> \par Literature
!>      http://www-math.mit.edu/~bazant/EDIP
!>      M.Z. Bazant & E. Kaxiras: Modeling of Covalent Bonding in Solids by
!>                                Inversion of Cohesive Energy Curves;
!>                                Phys. Rev. Lett. 77, 4370 (1996)
!>      M.Z. Bazant, E. Kaxiras and J.F. Justo: Environment-dependent interatomic
!>                                              potential for bulk silicon;
!>                                              Phys. Rev. B 56, 8542-8552 (1997)
!>      S. Goedecker: Optimization and parallelization of a force field for silicon
!>                    using OpenMP; CPC 148, 1 (2002)
!> \par History
!>      03.2006 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! **************************************************************************************************
   SUBROUTINE eip_bazant(eip_env)
      TYPE(eip_environment_type), POINTER                :: eip_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'eip_bazant'

      INTEGER                                            :: handle, i, iparticle, iparticle_kind, &
                                                            iparticle_local, iw, natom, &
                                                            nparticle_kind, nparticle_local
      REAL(KIND=dp)                                      :: ekin, ener, ener_var, mass
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: rxyz
      REAL(KIND=dp), DIMENSION(3)                        :: abc
      TYPE(atomic_kind_list_type), POINTER               :: atomic_kinds
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_subsys_type), POINTER                      :: subsys
      TYPE(distribution_1d_type), POINTER                :: local_particles
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(section_vals_type), POINTER                   :: eip_section

!   ------------------------------------------------------------------------

      CALL timeset(routineN, handle)

      NULLIFY (cell, particle_set, eip_section, logger, atomic_kinds, &
               atomic_kind, local_particles, subsys, atomic_kind_set, para_env)

      ekin = 0.0_dp

      logger => cp_get_default_logger()

      CPASSERT(ASSOCIATED(eip_env))

      CALL eip_env_get(eip_env=eip_env, cell=cell, particle_set=particle_set, &
                       subsys=subsys, local_particles=local_particles, &
                       atomic_kind_set=atomic_kind_set)
      CALL get_cell(cell=cell, abc=abc)

      eip_section => section_vals_get_subs_vals(eip_env%force_env_input, "EIP")
      natom = SIZE(particle_set)
      !natom = local_particles%n_el(1)

      ALLOCATE (rxyz(3, natom))

      DO i = 1, natom
         !iparticle = local_particles%list(1)%array(i)
         rxyz(:, i) = particle_set(i)%r(:)*angstrom
      END DO

      CALL eip_bazant_silicon(nat=natom, alat=abc*angstrom, rxyz0=rxyz, &
                              fxyz=eip_env%eip_forces, ener=ener, &
                              coord=eip_env%coord_avg, ener_var=ener_var, &
                              coord_var=eip_env%coord_var, count=eip_env%count)

      !CALL get_part_ke(md_env, tbmd_energy%E_kinetic, int_grp=globalenv%para_env)
      CALL cp_subsys_get(subsys=subsys, atomic_kinds=atomic_kinds)

      nparticle_kind = atomic_kinds%n_els

      DO iparticle_kind = 1, nparticle_kind
         atomic_kind => atomic_kind_set(iparticle_kind)
         CALL get_atomic_kind(atomic_kind=atomic_kind, mass=mass)
         nparticle_local = local_particles%n_el(iparticle_kind)
         DO iparticle_local = 1, nparticle_local
            iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
            ekin = ekin + 0.5_dp*mass* &
                   (particle_set(iparticle)%v(1)*particle_set(iparticle)%v(1) &
                    + particle_set(iparticle)%v(2)*particle_set(iparticle)%v(2) &
                    + particle_set(iparticle)%v(3)*particle_set(iparticle)%v(3))
         END DO
      END DO

      ! sum all contributions to energy over calculated parts on all processors
      CALL cp_subsys_get(subsys=subsys, para_env=para_env)
      CALL para_env%sum(ekin)
      eip_env%eip_kinetic_energy = ekin

      eip_env%eip_potential_energy = ener/evolt
      eip_env%eip_energy = eip_env%eip_kinetic_energy + eip_env%eip_potential_energy
      eip_env%eip_energy_var = ener_var/evolt

      DO i = 1, natom
         particle_set(i)%f(:) = eip_env%eip_forces(:, i)/evolt*angstrom
      END DO

      DEALLOCATE (rxyz)

      ! Print
      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                           eip_section, "PRINT%ENERGIES"), cp_p_file)) THEN
         iw = cp_print_key_unit_nr(logger, eip_section, "PRINT%ENERGIES", &
                                   extension=".mmLog")

         CALL eip_print_energies(eip_env=eip_env, output_unit=iw)
         CALL cp_print_key_finished_output(iw, logger, eip_section, &
                                           "PRINT%ENERGIES")
      END IF

      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                           eip_section, "PRINT%ENERGIES_VAR"), cp_p_file)) THEN
         iw = cp_print_key_unit_nr(logger, eip_section, "PRINT%ENERGIES_VAR", &
                                   extension=".mmLog")

         CALL eip_print_energy_var(eip_env=eip_env, output_unit=iw)
         CALL cp_print_key_finished_output(iw, logger, eip_section, &
                                           "PRINT%ENERGIES_VAR")
      END IF

      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                           eip_section, "PRINT%FORCES"), cp_p_file)) THEN
         iw = cp_print_key_unit_nr(logger, eip_section, "PRINT%FORCES", &
                                   extension=".mmLog")

         CALL eip_print_forces(eip_env=eip_env, output_unit=iw)
         CALL cp_print_key_finished_output(iw, logger, eip_section, &
                                           "PRINT%FORCES")
      END IF

      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                           eip_section, "PRINT%COORD_AVG"), cp_p_file)) THEN
         iw = cp_print_key_unit_nr(logger, eip_section, "PRINT%COORD_AVG", &
                                   extension=".mmLog")

         CALL eip_print_coord_avg(eip_env=eip_env, output_unit=iw)
         CALL cp_print_key_finished_output(iw, logger, eip_section, &
                                           "PRINT%COORD_AVG")
      END IF

      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                           eip_section, "PRINT%COORD_VAR"), cp_p_file)) THEN
         iw = cp_print_key_unit_nr(logger, eip_section, "PRINT%COORD_VAR", &
                                   extension=".mmLog")

         CALL eip_print_coord_var(eip_env=eip_env, output_unit=iw)
         CALL cp_print_key_finished_output(iw, logger, eip_section, &
                                           "PRINT%COORD_VAR")
      END IF

      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                           eip_section, "PRINT%COUNT"), cp_p_file)) THEN
         iw = cp_print_key_unit_nr(logger, eip_section, "PRINT%COUNT", &
                                   extension=".mmLog")

         CALL eip_print_count(eip_env=eip_env, output_unit=iw)
         CALL cp_print_key_finished_output(iw, logger, eip_section, &
                                           "PRINT%COUNT")
      END IF

      CALL timestop(handle)

   END SUBROUTINE eip_bazant

! **************************************************************************************************
!> \brief Interface routine of Goedecker's Lenosky force field to CP2K
!> \param eip_env ...
!> \par Literature
!>      T. Lenosky, et. al.: Highly optimized empirical potential model of silicon;
!>                           Modelling Simul. Sci. Eng., 8 (2000)
!>      S. Goedecker: Optimization and parallelization of a force field for silicon
!>                    using OpenMP; CPC 148, 1 (2002)
!> \par History
!>      03.2006 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! **************************************************************************************************
   SUBROUTINE eip_lenosky(eip_env)
      TYPE(eip_environment_type), POINTER                :: eip_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'eip_lenosky'

      INTEGER                                            :: handle, i, iparticle, iparticle_kind, &
                                                            iparticle_local, iw, natom, &
                                                            nparticle_kind, nparticle_local
      REAL(KIND=dp)                                      :: ekin, ener, ener_var, mass
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: rxyz
      REAL(KIND=dp), DIMENSION(3)                        :: abc
      TYPE(atomic_kind_list_type), POINTER               :: atomic_kinds
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_subsys_type), POINTER                      :: subsys
      TYPE(distribution_1d_type), POINTER                :: local_particles
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(section_vals_type), POINTER                   :: eip_section

!   ------------------------------------------------------------------------

      CALL timeset(routineN, handle)

      NULLIFY (cell, particle_set, eip_section, logger, atomic_kinds, &
               atomic_kind, local_particles, subsys, atomic_kind_set, para_env)

      ekin = 0.0_dp

      logger => cp_get_default_logger()

      CPASSERT(ASSOCIATED(eip_env))

      CALL eip_env_get(eip_env=eip_env, cell=cell, particle_set=particle_set, &
                       subsys=subsys, local_particles=local_particles, &
                       atomic_kind_set=atomic_kind_set)
      CALL get_cell(cell=cell, abc=abc)

      eip_section => section_vals_get_subs_vals(eip_env%force_env_input, "EIP")
      natom = SIZE(particle_set)
      !natom = local_particles%n_el(1)

      ALLOCATE (rxyz(3, natom))

      DO i = 1, natom
         !iparticle = local_particles%list(1)%array(i)
         rxyz(:, i) = particle_set(i)%r(:)*angstrom
      END DO

      CALL eip_lenosky_silicon(nat=natom, alat=abc*angstrom, rxyz0=rxyz, &
                               fxyz=eip_env%eip_forces, ener=ener, &
                               coord=eip_env%coord_avg, ener_var=ener_var, &
                               coord_var=eip_env%coord_var, count=eip_env%count)

      !CALL get_part_ke(md_env, tbmd_energy%E_kinetic, int_grp=globalenv%para_env)
      CALL cp_subsys_get(subsys=subsys, atomic_kinds=atomic_kinds)

      nparticle_kind = atomic_kinds%n_els

      DO iparticle_kind = 1, nparticle_kind
         atomic_kind => atomic_kind_set(iparticle_kind)
         CALL get_atomic_kind(atomic_kind=atomic_kind, mass=mass)
         nparticle_local = local_particles%n_el(iparticle_kind)
         DO iparticle_local = 1, nparticle_local
            iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
            ekin = ekin + 0.5_dp*mass* &
                   (particle_set(iparticle)%v(1)*particle_set(iparticle)%v(1) &
                    + particle_set(iparticle)%v(2)*particle_set(iparticle)%v(2) &
                    + particle_set(iparticle)%v(3)*particle_set(iparticle)%v(3))
         END DO
      END DO

      ! sum all contributions to energy over calculated parts on all processors
      CALL cp_subsys_get(subsys=subsys, para_env=para_env)
      CALL para_env%sum(ekin)
      eip_env%eip_kinetic_energy = ekin

      eip_env%eip_potential_energy = ener/evolt
      eip_env%eip_energy = eip_env%eip_kinetic_energy + eip_env%eip_potential_energy
      eip_env%eip_energy_var = ener_var/evolt

      DO i = 1, natom
         particle_set(i)%f(:) = eip_env%eip_forces(:, i)/evolt*angstrom
      END DO

      DEALLOCATE (rxyz)

      ! Print
      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                           eip_section, "PRINT%ENERGIES"), cp_p_file)) THEN
         iw = cp_print_key_unit_nr(logger, eip_section, "PRINT%ENERGIES", &
                                   extension=".mmLog")

         CALL eip_print_energies(eip_env=eip_env, output_unit=iw)
         CALL cp_print_key_finished_output(iw, logger, eip_section, &
                                           "PRINT%ENERGIES")
      END IF

      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                           eip_section, "PRINT%ENERGIES_VAR"), cp_p_file)) THEN
         iw = cp_print_key_unit_nr(logger, eip_section, "PRINT%ENERGIES_VAR", &
                                   extension=".mmLog")

         CALL eip_print_energy_var(eip_env=eip_env, output_unit=iw)
         CALL cp_print_key_finished_output(iw, logger, eip_section, &
                                           "PRINT%ENERGIES_VAR")
      END IF

      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                           eip_section, "PRINT%FORCES"), cp_p_file)) THEN
         iw = cp_print_key_unit_nr(logger, eip_section, "PRINT%FORCES", &
                                   extension=".mmLog")

         CALL eip_print_forces(eip_env=eip_env, output_unit=iw)
         CALL cp_print_key_finished_output(iw, logger, eip_section, &
                                           "PRINT%FORCES")
      END IF

      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                           eip_section, "PRINT%COORD_AVG"), cp_p_file)) THEN
         iw = cp_print_key_unit_nr(logger, eip_section, "PRINT%COORD_AVG", &
                                   extension=".mmLog")

         CALL eip_print_coord_avg(eip_env=eip_env, output_unit=iw)
         CALL cp_print_key_finished_output(iw, logger, eip_section, &
                                           "PRINT%COORD_AVG")
      END IF

      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                           eip_section, "PRINT%COORD_VAR"), cp_p_file)) THEN
         iw = cp_print_key_unit_nr(logger, eip_section, "PRINT%COORD_VAR", &
                                   extension=".mmLog")

         CALL eip_print_coord_var(eip_env=eip_env, output_unit=iw)
         CALL cp_print_key_finished_output(iw, logger, eip_section, &
                                           "PRINT%COORD_VAR")
      END IF

      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                           eip_section, "PRINT%COUNT"), cp_p_file)) THEN
         iw = cp_print_key_unit_nr(logger, eip_section, "PRINT%COUNT", &
                                   extension=".mmLog")

         CALL eip_print_count(eip_env=eip_env, output_unit=iw)
         CALL cp_print_key_finished_output(iw, logger, eip_section, &
                                           "PRINT%COUNT")
      END IF

      CALL timestop(handle)

   END SUBROUTINE eip_lenosky

! **************************************************************************************************
!> \brief Interface routine of the Stillinger-Weber force field to CP2K
!> \param eip_env ...
!> \par Literature
!>      F.H. Stillinger and T.A. Weber:
!>      Computer simulation of local order in condensed phases of silicon;
!>      Phys. Rev. B 31, 5262 (1985)
!> \par History
!>      04.2026 added [codex]
! **************************************************************************************************
   SUBROUTINE eip_stillinger_weber(eip_env)
      TYPE(eip_environment_type), POINTER                :: eip_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'eip_stillinger_weber'

      INTEGER                                            :: handle, i, iparticle, iparticle_kind, &
                                                            iparticle_local, iw, natom, &
                                                            nparticle_kind, nparticle_local
      REAL(KIND=dp)                                      :: ekin, ener, mass
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: rxyz
      REAL(KIND=dp), DIMENSION(3)                        :: abc
      TYPE(atomic_kind_list_type), POINTER               :: atomic_kinds
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_subsys_type), POINTER                      :: subsys
      TYPE(distribution_1d_type), POINTER                :: local_particles
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(section_vals_type), POINTER                   :: eip_section

      CALL timeset(routineN, handle)

      NULLIFY (cell, particle_set, eip_section, logger, atomic_kinds, &
               atomic_kind, local_particles, subsys, atomic_kind_set, para_env)

      ekin = 0.0_dp

      logger => cp_get_default_logger()

      CPASSERT(ASSOCIATED(eip_env))

      CALL eip_env_get(eip_env=eip_env, cell=cell, particle_set=particle_set, &
                       subsys=subsys, local_particles=local_particles, &
                       atomic_kind_set=atomic_kind_set)
      CALL get_cell(cell=cell, abc=abc)

      eip_section => section_vals_get_subs_vals(eip_env%force_env_input, "EIP")
      natom = SIZE(particle_set)

      ALLOCATE (rxyz(3, natom))

      DO i = 1, natom
         rxyz(:, i) = particle_set(i)%r(:)*angstrom
      END DO

      CALL eip_stillinger_weber_silicon(nat=natom, alat=abc*angstrom, &
                                             rxyz0=rxyz, fxyz=eip_env%eip_forces, &
                                             etot=ener, count=eip_env%count)

      eip_env%coord_avg = 0.0_dp
      eip_env%coord_var = 0.0_dp

      CALL cp_subsys_get(subsys=subsys, atomic_kinds=atomic_kinds)

      nparticle_kind = atomic_kinds%n_els

      DO iparticle_kind = 1, nparticle_kind
         atomic_kind => atomic_kind_set(iparticle_kind)
         CALL get_atomic_kind(atomic_kind=atomic_kind, mass=mass)
         nparticle_local = local_particles%n_el(iparticle_kind)
         DO iparticle_local = 1, nparticle_local
            iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
            ekin = ekin + 0.5_dp*mass* &
                   (particle_set(iparticle)%v(1)*particle_set(iparticle)%v(1) &
                    + particle_set(iparticle)%v(2)*particle_set(iparticle)%v(2) &
                    + particle_set(iparticle)%v(3)*particle_set(iparticle)%v(3))
         END DO
      END DO

      CALL cp_subsys_get(subsys=subsys, para_env=para_env)
      CALL para_env%sum(ekin)
      eip_env%eip_kinetic_energy = ekin

      eip_env%eip_potential_energy = ener/evolt
      eip_env%eip_energy = eip_env%eip_kinetic_energy + eip_env%eip_potential_energy
      eip_env%eip_energy_var = 0.0_dp

      DO i = 1, natom
         particle_set(i)%f(:) = eip_env%eip_forces(:, i)/evolt*angstrom
      END DO

      DEALLOCATE (rxyz)

      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                           eip_section, "PRINT%ENERGIES"), cp_p_file)) THEN
         iw = cp_print_key_unit_nr(logger, eip_section, "PRINT%ENERGIES", &
                                   extension=".mmLog")

         CALL eip_print_energies(eip_env=eip_env, output_unit=iw)
         CALL cp_print_key_finished_output(iw, logger, eip_section, &
                                           "PRINT%ENERGIES")
      END IF

      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                           eip_section, "PRINT%ENERGIES_VAR"), cp_p_file)) THEN
         iw = cp_print_key_unit_nr(logger, eip_section, "PRINT%ENERGIES_VAR", &
                                   extension=".mmLog")

         CALL eip_print_energy_var(eip_env=eip_env, output_unit=iw)
         CALL cp_print_key_finished_output(iw, logger, eip_section, &
                                           "PRINT%ENERGIES_VAR")
      END IF

      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                           eip_section, "PRINT%FORCES"), cp_p_file)) THEN
         iw = cp_print_key_unit_nr(logger, eip_section, "PRINT%FORCES", &
                                   extension=".mmLog")

         CALL eip_print_forces(eip_env=eip_env, output_unit=iw)
         CALL cp_print_key_finished_output(iw, logger, eip_section, &
                                           "PRINT%FORCES")
      END IF

      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                           eip_section, "PRINT%COORD_AVG"), cp_p_file)) THEN
         iw = cp_print_key_unit_nr(logger, eip_section, "PRINT%COORD_AVG", &
                                   extension=".mmLog")

         CALL eip_print_coord_avg(eip_env=eip_env, output_unit=iw)
         CALL cp_print_key_finished_output(iw, logger, eip_section, &
                                           "PRINT%COORD_AVG")
      END IF

      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                           eip_section, "PRINT%COORD_VAR"), cp_p_file)) THEN
         iw = cp_print_key_unit_nr(logger, eip_section, "PRINT%COORD_VAR", &
                                   extension=".mmLog")

         CALL eip_print_coord_var(eip_env=eip_env, output_unit=iw)
         CALL cp_print_key_finished_output(iw, logger, eip_section, &
                                           "PRINT%COORD_VAR")
      END IF

      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                           eip_section, "PRINT%COUNT"), cp_p_file)) THEN
         iw = cp_print_key_unit_nr(logger, eip_section, "PRINT%COUNT", &
                                   extension=".mmLog")

         CALL eip_print_count(eip_env=eip_env, output_unit=iw)
         CALL cp_print_key_finished_output(iw, logger, eip_section, &
                                           "PRINT%COUNT")
      END IF

      CALL timestop(handle)

   END SUBROUTINE eip_stillinger_weber

! **************************************************************************************************
!> \brief Interface routine of the Tersoff force field to CP2K
!> \param eip_env ...
!> \par Literature
!>      J. Tersoff:
!>      New empirical approach for the structure and energy of covalent systems;
!>      Phys. Rev. Lett. 61, 2879 (1988)
!>      J. Tersoff:
!>      Modeling solid-state chemistry: Interatomic potentials for multicomponent systems;
!>      Phys. Rev. B 39, 5566 (1989)
!> \par History
!>      04.2026 added [codex]
! **************************************************************************************************
   SUBROUTINE eip_tersoff(eip_env)
      TYPE(eip_environment_type), POINTER                :: eip_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'eip_tersoff'

      INTEGER                                            :: handle, i, iparticle, iparticle_kind, &
                                                            iparticle_local, iw, natom, &
                                                            nparticle_kind, nparticle_local
      REAL(KIND=dp)                                      :: ekin, ener, mass
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: rxyz
      REAL(KIND=dp), DIMENSION(3)                        :: abc
      TYPE(atomic_kind_list_type), POINTER               :: atomic_kinds
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_subsys_type), POINTER                      :: subsys
      TYPE(distribution_1d_type), POINTER                :: local_particles
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(section_vals_type), POINTER                   :: eip_section

      CALL timeset(routineN, handle)

      NULLIFY (cell, particle_set, eip_section, logger, atomic_kinds, &
               atomic_kind, local_particles, subsys, atomic_kind_set, para_env)

      ekin = 0.0_dp

      logger => cp_get_default_logger()

      CPASSERT(ASSOCIATED(eip_env))

      CALL eip_env_get(eip_env=eip_env, cell=cell, particle_set=particle_set, &
                       subsys=subsys, local_particles=local_particles, &
                       atomic_kind_set=atomic_kind_set)
      CALL get_cell(cell=cell, abc=abc)

      eip_section => section_vals_get_subs_vals(eip_env%force_env_input, "EIP")
      natom = SIZE(particle_set)

      ALLOCATE (rxyz(3, natom))

      DO i = 1, natom
         rxyz(:, i) = particle_set(i)%r(:)*angstrom
      END DO

      CALL eip_tersoff_silicon(nat=natom, alat=abc*angstrom, rxyz=rxyz, &
                                    fxyz=eip_env%eip_forces, etot=ener, &
                                    count=eip_env%count)

      eip_env%coord_avg = 0.0_dp
      eip_env%coord_var = 0.0_dp

      CALL cp_subsys_get(subsys=subsys, atomic_kinds=atomic_kinds)

      nparticle_kind = atomic_kinds%n_els

      DO iparticle_kind = 1, nparticle_kind
         atomic_kind => atomic_kind_set(iparticle_kind)
         CALL get_atomic_kind(atomic_kind=atomic_kind, mass=mass)
         nparticle_local = local_particles%n_el(iparticle_kind)
         DO iparticle_local = 1, nparticle_local
            iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
            ekin = ekin + 0.5_dp*mass* &
                   (particle_set(iparticle)%v(1)*particle_set(iparticle)%v(1) &
                    + particle_set(iparticle)%v(2)*particle_set(iparticle)%v(2) &
                    + particle_set(iparticle)%v(3)*particle_set(iparticle)%v(3))
         END DO
      END DO

      CALL cp_subsys_get(subsys=subsys, para_env=para_env)
      CALL para_env%sum(ekin)
      eip_env%eip_kinetic_energy = ekin

      eip_env%eip_potential_energy = ener/evolt
      eip_env%eip_energy = eip_env%eip_kinetic_energy + eip_env%eip_potential_energy
      eip_env%eip_energy_var = 0.0_dp

      DO i = 1, natom
         particle_set(i)%f(:) = eip_env%eip_forces(:, i)/evolt*angstrom
      END DO

      DEALLOCATE (rxyz)

      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                           eip_section, "PRINT%ENERGIES"), cp_p_file)) THEN
         iw = cp_print_key_unit_nr(logger, eip_section, "PRINT%ENERGIES", &
                                   extension=".mmLog")

         CALL eip_print_energies(eip_env=eip_env, output_unit=iw)
         CALL cp_print_key_finished_output(iw, logger, eip_section, &
                                           "PRINT%ENERGIES")
      END IF

      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                           eip_section, "PRINT%ENERGIES_VAR"), cp_p_file)) THEN
         iw = cp_print_key_unit_nr(logger, eip_section, "PRINT%ENERGIES_VAR", &
                                   extension=".mmLog")

         CALL eip_print_energy_var(eip_env=eip_env, output_unit=iw)
         CALL cp_print_key_finished_output(iw, logger, eip_section, &
                                           "PRINT%ENERGIES_VAR")
      END IF

      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                           eip_section, "PRINT%FORCES"), cp_p_file)) THEN
         iw = cp_print_key_unit_nr(logger, eip_section, "PRINT%FORCES", &
                                   extension=".mmLog")

         CALL eip_print_forces(eip_env=eip_env, output_unit=iw)
         CALL cp_print_key_finished_output(iw, logger, eip_section, &
                                           "PRINT%FORCES")
      END IF

      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                           eip_section, "PRINT%COORD_AVG"), cp_p_file)) THEN
         iw = cp_print_key_unit_nr(logger, eip_section, "PRINT%COORD_AVG", &
                                   extension=".mmLog")

         CALL eip_print_coord_avg(eip_env=eip_env, output_unit=iw)
         CALL cp_print_key_finished_output(iw, logger, eip_section, &
                                           "PRINT%COORD_AVG")
      END IF

      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                           eip_section, "PRINT%COORD_VAR"), cp_p_file)) THEN
         iw = cp_print_key_unit_nr(logger, eip_section, "PRINT%COORD_VAR", &
                                   extension=".mmLog")

         CALL eip_print_coord_var(eip_env=eip_env, output_unit=iw)
         CALL cp_print_key_finished_output(iw, logger, eip_section, &
                                           "PRINT%COORD_VAR")
      END IF

      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                           eip_section, "PRINT%COUNT"), cp_p_file)) THEN
         iw = cp_print_key_unit_nr(logger, eip_section, "PRINT%COUNT", &
                                   extension=".mmLog")

         CALL eip_print_count(eip_env=eip_env, output_unit=iw)
         CALL cp_print_key_finished_output(iw, logger, eip_section, &
                                           "PRINT%COUNT")
      END IF

      CALL timestop(handle)

   END SUBROUTINE eip_tersoff

! **************************************************************************************************
!> \brief Print routine for the EIP energies
!> \param eip_env The eip environment of matter
!> \param output_unit The output unit
!> \par History
!>      03.2006 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!> \note
!>      As usual the EIP energies differ from the DFT energies!
!>      Only the relative energy differences are correctly reproduced.
! **************************************************************************************************
   SUBROUTINE eip_print_energies(eip_env, output_unit)
      TYPE(eip_environment_type), POINTER                :: eip_env
      INTEGER, INTENT(IN)                                :: output_unit

!   ------------------------------------------------------------------------

      IF (output_unit > 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,(T3,A,T55,F25.14))") &
            "Kinetic energy [Hartree]:        ", eip_env%eip_kinetic_energy, &
            "Potential energy [Hartree]:      ", eip_env%eip_potential_energy, &
            "Total EIP energy [Hartree]:      ", eip_env%eip_energy
      END IF

   END SUBROUTINE eip_print_energies

! **************************************************************************************************
!> \brief Print routine for the variance of the energy/atom
!> \param eip_env The eip environment of matter
!> \param output_unit The output unit
!> \par History
!>      03.2006 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! **************************************************************************************************
   SUBROUTINE eip_print_energy_var(eip_env, output_unit)
      TYPE(eip_environment_type), POINTER                :: eip_env
      INTEGER, INTENT(IN)                                :: output_unit

      INTEGER                                            :: unit_nr

!   ------------------------------------------------------------------------

      unit_nr = output_unit

      IF (unit_nr > 0) THEN

         WRITE (unit_nr, *) ""
         WRITE (unit_nr, *) "The variance of the EIP energy/atom!"
         WRITE (unit_nr, *) ""
         WRITE (unit_nr, *) eip_env%eip_energy_var

      END IF

   END SUBROUTINE eip_print_energy_var

! **************************************************************************************************
!> \brief Print routine for the forces
!> \param eip_env The eip environment of matter
!> \param output_unit The output unit
!> \par History
!>      03.2006 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! **************************************************************************************************
   SUBROUTINE eip_print_forces(eip_env, output_unit)
      TYPE(eip_environment_type), POINTER                :: eip_env
      INTEGER, INTENT(IN)                                :: output_unit

      INTEGER                                            :: iatom, natom, unit_nr
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

!   ------------------------------------------------------------------------

      NULLIFY (particle_set)

      unit_nr = output_unit

      IF (unit_nr > 0) THEN

         CALL eip_env_get(eip_env=eip_env, particle_set=particle_set)

         natom = SIZE(particle_set)

         WRITE (unit_nr, *) ""
         WRITE (unit_nr, *) "The EIP forces!"
         WRITE (unit_nr, *) ""
         WRITE (unit_nr, *) "Total EIP forces [Hartree/Bohr]"
         DO iatom = 1, natom
            WRITE (unit_nr, *) eip_env%eip_forces(1:3, iatom)
         END DO

      END IF

   END SUBROUTINE eip_print_forces

! **************************************************************************************************
!> \brief Print routine for the average coordination number
!> \param eip_env The eip environment of matter
!> \param output_unit The output unit
!> \par History
!>      03.2006 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! **************************************************************************************************
   SUBROUTINE eip_print_coord_avg(eip_env, output_unit)
      TYPE(eip_environment_type), POINTER                :: eip_env
      INTEGER, INTENT(IN)                                :: output_unit

      INTEGER                                            :: unit_nr

!   ------------------------------------------------------------------------

      unit_nr = output_unit

      IF (unit_nr > 0) THEN

         WRITE (unit_nr, *) ""
         WRITE (unit_nr, *) "The average coordination number!"
         WRITE (unit_nr, *) ""
         WRITE (unit_nr, *) eip_env%coord_avg

      END IF

   END SUBROUTINE eip_print_coord_avg

! **************************************************************************************************
!> \brief Print routine for the variance of the coordination number
!> \param eip_env The eip environment of matter
!> \param output_unit The output unit
!> \par History
!>      03.2006 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! **************************************************************************************************
   SUBROUTINE eip_print_coord_var(eip_env, output_unit)
      TYPE(eip_environment_type), POINTER                :: eip_env
      INTEGER, INTENT(IN)                                :: output_unit

      INTEGER                                            :: unit_nr

!   ------------------------------------------------------------------------

      unit_nr = output_unit

      IF (unit_nr > 0) THEN

         WRITE (unit_nr, *) ""
         WRITE (unit_nr, *) "The variance of the coordination number!"
         WRITE (unit_nr, *) ""
         WRITE (unit_nr, *) eip_env%coord_var

      END IF

   END SUBROUTINE eip_print_coord_var

! **************************************************************************************************
!> \brief Print routine for the function call counter
!> \param eip_env The eip environment of matter
!> \param output_unit The output unit
!> \par History
!>      03.2006 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! **************************************************************************************************
   SUBROUTINE eip_print_count(eip_env, output_unit)
      TYPE(eip_environment_type), POINTER                :: eip_env
      INTEGER, INTENT(IN)                                :: output_unit

      INTEGER                                            :: unit_nr

!   ------------------------------------------------------------------------

      unit_nr = output_unit

      IF (unit_nr > 0) THEN

         WRITE (unit_nr, *) ""
         WRITE (unit_nr, *) "The function call counter!"
         WRITE (unit_nr, *) ""
         WRITE (unit_nr, *) eip_env%count

      END IF

   END SUBROUTINE eip_print_count

! **************************************************************************************************
!> \brief Bazant's EDIP (environment-dependent interatomic potential) for Silicon
!>      by Stefan Goedecker
!> \param nat number of atoms
!> \param alat lattice constants of the orthorombic box containing the particles
!> \param rxyz0 atomic positions in Angstrom, may be modified on output.
!>               If an atom is outside the box the program will bring it back
!>               into the box by translations through alat
!> \param fxyz forces in eV/A
!> \param ener total energy in eV
!> \param coord average coordination number
!> \param ener_var variance of the energy/atom
!> \param coord_var variance of the coordination number
!> \param count count is increased by one per call, has to be initialized
!>                to 0.e0_dp before first call of eip_bazant
!> \par Literature
!>      http://www-math.mit.edu/~bazant/EDIP
!>      M.Z. Bazant & E. Kaxiras: Modeling of Covalent Bonding in Solids by
!>                                Inversion of Cohesive Energy Curves;
!>                                Phys. Rev. Lett. 77, 4370 (1996)
!>      M.Z. Bazant, E. Kaxiras and J.F. Justo: Environment-dependent interatomic
!>                                              potential for bulk silicon;
!>                                              Phys. Rev. B 56, 8542-8552 (1997)
!>      S. Goedecker: Optimization and parallelization of a force field for silicon
!>                    using OpenMP; CPC 148, 1 (2002)
!> \par History
!>      03.2006 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! **************************************************************************************************
   SUBROUTINE eip_bazant_silicon(nat, alat, rxyz0, fxyz, ener, coord, ener_var, &
                                 coord_var, count)

      INTEGER                                            :: nat
      REAL(KIND=dp)                                      :: alat, rxyz0, fxyz, ener, coord, &
                                                            ener_var, coord_var, count

      DIMENSION rxyz0(3, nat), fxyz(3, nat), alat(3)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :) :: rxyz
      INTEGER, ALLOCATABLE, DIMENSION(:, :)       :: lsta
      INTEGER, ALLOCATABLE, DIMENSION(:)         :: lstb
      INTEGER, ALLOCATABLE, DIMENSION(:)         :: lay
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :, :)   :: icell
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :) :: rel
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :) :: txyz
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :) :: s2, s3, sz
      INTEGER, ALLOCATABLE, DIMENSION(:)         :: num2, num3, numz

      REAL(KIND=dp) :: coord2, cut, cut2, ener2, rlc1i, rlc2i, rlc3i, tcoord, &
                       tcoord2, tener, tener2
      INTEGER       :: iam, iat, iat1, iat2, ii, i, il, in, indlst, indlstx, istop, &
                       istopg, l2, l3, laymx, ll1, ll2, ll3, lot, max_nbrs, myspace, &
                       l1, myspaceout, ncx, nn, nnbrx, npr

!        cut=par_a
      cut = 3.1213820e0_dp + 1.e-14_dp

      IF (count == 0) OPEN (unit=10, file='bazant.mon', status='unknown')
      count = count + 1.e0_dp

! linear scaling calculation of verlet list
      ll1 = INT(alat(1)/cut)
      IF (ll1 < 1) CPABORT("alat(1) too small")
      ll2 = INT(alat(2)/cut)
      IF (ll2 < 1) CPABORT("alat(2) too small")
      ll3 = INT(alat(3)/cut)
      IF (ll3 < 1) CPABORT("alat(3) too small")

! determine number of threads
      npr = 1
!$OMP PARALLEL PRIVATE(iam)  SHARED (npr) DEFAULT(NONE)
!$    iam = omp_get_thread_num()
!$    if (iam .eq. 0) npr = omp_get_num_threads()
!$OMP END PARALLEL

! linear scaling calculation of verlet list

      IF (npr <= 1) THEN !serial if too few processors to gain by parallelizing

! set ncx for serial case, ncx for parallel case set below
         ncx = 16
         loop_ncx_s: DO
            ALLOCATE (icell(0:ncx, -1:ll1, -1:ll2, -1:ll3))
            icell(0, -1:ll1, -1:ll2, -1:ll3) = 0
            rlc1i = ll1/alat(1)
            rlc2i = ll2/alat(2)
            rlc3i = ll3/alat(3)

            loop_iat_s: DO iat = 1, nat
               rxyz0(1, iat) = MODULO(MODULO(rxyz0(1, iat), alat(1)), alat(1))
               rxyz0(2, iat) = MODULO(MODULO(rxyz0(2, iat), alat(2)), alat(2))
               rxyz0(3, iat) = MODULO(MODULO(rxyz0(3, iat), alat(3)), alat(3))
               l1 = INT(rxyz0(1, iat)*rlc1i)
               l2 = INT(rxyz0(2, iat)*rlc2i)
               l3 = INT(rxyz0(3, iat)*rlc3i)

               ii = icell(0, l1, l2, l3)
               ii = ii + 1
               icell(0, l1, l2, l3) = ii
               IF (ii > ncx) THEN
                  WRITE (10, *) count, 'NCX too small', ncx
                  DEALLOCATE (icell)
                  ncx = ncx*2
                  CYCLE loop_ncx_s
               END IF
               icell(ii, l1, l2, l3) = iat
            END DO loop_iat_s
            EXIT loop_ncx_s
         END DO loop_ncx_s

      ELSE ! parallel case

! periodization of particles can be done in parallel
!$OMP PARALLEL DO SHARED (alat,nat,rxyz0) PRIVATE(iat) DEFAULT(NONE)
         DO iat = 1, nat
            rxyz0(1, iat) = MODULO(MODULO(rxyz0(1, iat), alat(1)), alat(1))
            rxyz0(2, iat) = MODULO(MODULO(rxyz0(2, iat), alat(2)), alat(2))
            rxyz0(3, iat) = MODULO(MODULO(rxyz0(3, iat), alat(3)), alat(3))
         END DO
!$OMP END PARALLEL DO

! assignment to cell is done serially
! set ncx for parallel case, ncx for serial case set above
         ncx = 16
         loop_ncx_p: DO
            ALLOCATE (icell(0:ncx, -1:ll1, -1:ll2, -1:ll3))
            icell(0, -1:ll1, -1:ll2, -1:ll3) = 0

            rlc1i = ll1/alat(1)
            rlc2i = ll2/alat(2)
            rlc3i = ll3/alat(3)

            loop_iat_p: DO iat = 1, nat
               l1 = INT(rxyz0(1, iat)*rlc1i)
               l2 = INT(rxyz0(2, iat)*rlc2i)
               l3 = INT(rxyz0(3, iat)*rlc3i)
               ii = icell(0, l1, l2, l3)
               ii = ii + 1
               icell(0, l1, l2, l3) = ii
               IF (ii > ncx) THEN
                  WRITE (10, *) count, 'NCX too small', ncx
                  DEALLOCATE (icell)
                  ncx = ncx*2
                  CYCLE loop_ncx_p
               END IF
               icell(ii, l1, l2, l3) = iat
            END DO loop_iat_p
            EXIT loop_ncx_p
         END DO loop_ncx_p

      END IF

! duplicate all atoms within boundary layer
      laymx = ncx*(2*ll1*ll2 + 2*ll1*ll3 + 2*ll2*ll3 + 4*ll1 + 4*ll2 + 4*ll3 + 8)
      nn = nat + laymx
      ALLOCATE (rxyz(3, nn), lay(nn))
      DO iat = 1, nat
         lay(iat) = iat
         rxyz(1, iat) = rxyz0(1, iat)
         rxyz(2, iat) = rxyz0(2, iat)
         rxyz(3, iat) = rxyz0(3, iat)
      END DO
      il = nat
! xy plane
      DO l2 = 0, ll2 - 1
      DO l1 = 0, ll1 - 1

         in = icell(0, l1, l2, 0)
         icell(0, l1, l2, ll3) = in
         DO ii = 1, in
            i = icell(ii, l1, l2, 0)
            il = il + 1
            IF (il > nn) CPABORT("enlarge laymx")
            lay(il) = i
            icell(ii, l1, l2, ll3) = il
            rxyz(1, il) = rxyz(1, i)
            rxyz(2, il) = rxyz(2, i)
            rxyz(3, il) = rxyz(3, i) + alat(3)
         END DO

         in = icell(0, l1, l2, ll3 - 1)
         icell(0, l1, l2, -1) = in
         DO ii = 1, in
            i = icell(ii, l1, l2, ll3 - 1)
            il = il + 1
            IF (il > nn) CPABORT("enlarge laymx")
            lay(il) = i
            icell(ii, l1, l2, -1) = il
            rxyz(1, il) = rxyz(1, i)
            rxyz(2, il) = rxyz(2, i)
            rxyz(3, il) = rxyz(3, i) - alat(3)
         END DO

      END DO
      END DO

! yz plane
      DO l3 = 0, ll3 - 1
      DO l2 = 0, ll2 - 1

         in = icell(0, 0, l2, l3)
         icell(0, ll1, l2, l3) = in
         DO ii = 1, in
            i = icell(ii, 0, l2, l3)
            il = il + 1
            IF (il > nn) CPABORT("enlarge laymx")
            lay(il) = i
            icell(ii, ll1, l2, l3) = il
            rxyz(1, il) = rxyz(1, i) + alat(1)
            rxyz(2, il) = rxyz(2, i)
            rxyz(3, il) = rxyz(3, i)
         END DO

         in = icell(0, ll1 - 1, l2, l3)
         icell(0, -1, l2, l3) = in
         DO ii = 1, in
            i = icell(ii, ll1 - 1, l2, l3)
            il = il + 1
            IF (il > nn) CPABORT("enlarge laymx")
            lay(il) = i
            icell(ii, -1, l2, l3) = il
            rxyz(1, il) = rxyz(1, i) - alat(1)
            rxyz(2, il) = rxyz(2, i)
            rxyz(3, il) = rxyz(3, i)
         END DO

      END DO
      END DO

! xz plane
      DO l3 = 0, ll3 - 1
      DO l1 = 0, ll1 - 1

         in = icell(0, l1, 0, l3)
         icell(0, l1, ll2, l3) = in
         DO ii = 1, in
            i = icell(ii, l1, 0, l3)
            il = il + 1
            IF (il > nn) CPABORT("enlarge laymx")
            lay(il) = i
            icell(ii, l1, ll2, l3) = il
            rxyz(1, il) = rxyz(1, i)
            rxyz(2, il) = rxyz(2, i) + alat(2)
            rxyz(3, il) = rxyz(3, i)
         END DO

         in = icell(0, l1, ll2 - 1, l3)
         icell(0, l1, -1, l3) = in
         DO ii = 1, in
            i = icell(ii, l1, ll2 - 1, l3)
            il = il + 1
            IF (il > nn) CPABORT("enlarge laymx")
            lay(il) = i
            icell(ii, l1, -1, l3) = il
            rxyz(1, il) = rxyz(1, i)
            rxyz(2, il) = rxyz(2, i) - alat(2)
            rxyz(3, il) = rxyz(3, i)
         END DO

      END DO
      END DO

! x axis
      DO l1 = 0, ll1 - 1

         in = icell(0, l1, 0, 0)
         icell(0, l1, ll2, ll3) = in
         DO ii = 1, in
            i = icell(ii, l1, 0, 0)
            il = il + 1
            IF (il > nn) CPABORT("enlarge laymx")
            lay(il) = i
            icell(ii, l1, ll2, ll3) = il
            rxyz(1, il) = rxyz(1, i)
            rxyz(2, il) = rxyz(2, i) + alat(2)
            rxyz(3, il) = rxyz(3, i) + alat(3)
         END DO

         in = icell(0, l1, 0, ll3 - 1)
         icell(0, l1, ll2, -1) = in
         DO ii = 1, in
            i = icell(ii, l1, 0, ll3 - 1)
            il = il + 1
            IF (il > nn) CPABORT("enlarge laymx")
            lay(il) = i
            icell(ii, l1, ll2, -1) = il
            rxyz(1, il) = rxyz(1, i)
            rxyz(2, il) = rxyz(2, i) + alat(2)
            rxyz(3, il) = rxyz(3, i) - alat(3)
         END DO

         in = icell(0, l1, ll2 - 1, 0)
         icell(0, l1, -1, ll3) = in
         DO ii = 1, in
            i = icell(ii, l1, ll2 - 1, 0)
            il = il + 1
            IF (il > nn) CPABORT("enlarge laymx")
            lay(il) = i
            icell(ii, l1, -1, ll3) = il
            rxyz(1, il) = rxyz(1, i)
            rxyz(2, il) = rxyz(2, i) - alat(2)
            rxyz(3, il) = rxyz(3, i) + alat(3)
         END DO

         in = icell(0, l1, ll2 - 1, ll3 - 1)
         icell(0, l1, -1, -1) = in
         DO ii = 1, in
            i = icell(ii, l1, ll2 - 1, ll3 - 1)
            il = il + 1
            IF (il > nn) CPABORT("enlarge laymx")
            lay(il) = i
            icell(ii, l1, -1, -1) = il
            rxyz(1, il) = rxyz(1, i)
            rxyz(2, il) = rxyz(2, i) - alat(2)
            rxyz(3, il) = rxyz(3, i) - alat(3)
         END DO

      END DO

! y axis
      DO l2 = 0, ll2 - 1

         in = icell(0, 0, l2, 0)
         icell(0, ll1, l2, ll3) = in
         DO ii = 1, in
            i = icell(ii, 0, l2, 0)
            il = il + 1
            IF (il > nn) CPABORT("enlarge laymx")
            lay(il) = i
            icell(ii, ll1, l2, ll3) = il
            rxyz(1, il) = rxyz(1, i) + alat(1)
            rxyz(2, il) = rxyz(2, i)
            rxyz(3, il) = rxyz(3, i) + alat(3)
         END DO

         in = icell(0, 0, l2, ll3 - 1)
         icell(0, ll1, l2, -1) = in
         DO ii = 1, in
            i = icell(ii, 0, l2, ll3 - 1)
            il = il + 1
            IF (il > nn) CPABORT("enlarge laymx")
            lay(il) = i
            icell(ii, ll1, l2, -1) = il
            rxyz(1, il) = rxyz(1, i) + alat(1)
            rxyz(2, il) = rxyz(2, i)
            rxyz(3, il) = rxyz(3, i) - alat(3)
         END DO

         in = icell(0, ll1 - 1, l2, 0)
         icell(0, -1, l2, ll3) = in
         DO ii = 1, in
            i = icell(ii, ll1 - 1, l2, 0)
            il = il + 1
            IF (il > nn) CPABORT("enlarge laymx")
            lay(il) = i
            icell(ii, -1, l2, ll3) = il
            rxyz(1, il) = rxyz(1, i) - alat(1)
            rxyz(2, il) = rxyz(2, i)
            rxyz(3, il) = rxyz(3, i) + alat(3)
         END DO

         in = icell(0, ll1 - 1, l2, ll3 - 1)
         icell(0, -1, l2, -1) = in
         DO ii = 1, in
            i = icell(ii, ll1 - 1, l2, ll3 - 1)
            il = il + 1
            IF (il > nn) CPABORT("enlarge laymx")
            lay(il) = i
            icell(ii, -1, l2, -1) = il
            rxyz(1, il) = rxyz(1, i) - alat(1)
            rxyz(2, il) = rxyz(2, i)
            rxyz(3, il) = rxyz(3, i) - alat(3)
         END DO

      END DO

! z axis
      DO l3 = 0, ll3 - 1

         in = icell(0, 0, 0, l3)
         icell(0, ll1, ll2, l3) = in
         DO ii = 1, in
            i = icell(ii, 0, 0, l3)
            il = il + 1
            IF (il > nn) CPABORT("enlarge laymx")
            lay(il) = i
            icell(ii, ll1, ll2, l3) = il
            rxyz(1, il) = rxyz(1, i) + alat(1)
            rxyz(2, il) = rxyz(2, i) + alat(2)
            rxyz(3, il) = rxyz(3, i)
         END DO

         in = icell(0, ll1 - 1, 0, l3)
         icell(0, -1, ll2, l3) = in
         DO ii = 1, in
            i = icell(ii, ll1 - 1, 0, l3)
            il = il + 1
            IF (il > nn) CPABORT("enlarge laymx")
            lay(il) = i
            icell(ii, -1, ll2, l3) = il
            rxyz(1, il) = rxyz(1, i) - alat(1)
            rxyz(2, il) = rxyz(2, i) + alat(2)
            rxyz(3, il) = rxyz(3, i)
         END DO

         in = icell(0, 0, ll2 - 1, l3)
         icell(0, ll1, -1, l3) = in
         DO ii = 1, in
            i = icell(ii, 0, ll2 - 1, l3)
            il = il + 1
            IF (il > nn) CPABORT("enlarge laymx")
            lay(il) = i
            icell(ii, ll1, -1, l3) = il
            rxyz(1, il) = rxyz(1, i) + alat(1)
            rxyz(2, il) = rxyz(2, i) - alat(2)
            rxyz(3, il) = rxyz(3, i)
         END DO

         in = icell(0, ll1 - 1, ll2 - 1, l3)
         icell(0, -1, -1, l3) = in
         DO ii = 1, in
            i = icell(ii, ll1 - 1, ll2 - 1, l3)
            il = il + 1
            IF (il > nn) CPABORT("enlarge laymx")
            lay(il) = i
            icell(ii, -1, -1, l3) = il
            rxyz(1, il) = rxyz(1, i) - alat(1)
            rxyz(2, il) = rxyz(2, i) - alat(2)
            rxyz(3, il) = rxyz(3, i)
         END DO

      END DO

! corners
      in = icell(0, 0, 0, 0)
      icell(0, ll1, ll2, ll3) = in
      DO ii = 1, in
         i = icell(ii, 0, 0, 0)
         il = il + 1
         IF (il > nn) CPABORT("enlarge laymx")
         lay(il) = i
         icell(ii, ll1, ll2, ll3) = il
         rxyz(1, il) = rxyz(1, i) + alat(1)
         rxyz(2, il) = rxyz(2, i) + alat(2)
         rxyz(3, il) = rxyz(3, i) + alat(3)
      END DO

      in = icell(0, ll1 - 1, 0, 0)
      icell(0, -1, ll2, ll3) = in
      DO ii = 1, in
         i = icell(ii, ll1 - 1, 0, 0)
         il = il + 1
         IF (il > nn) CPABORT("enlarge laymx")
         lay(il) = i
         icell(ii, -1, ll2, ll3) = il
         rxyz(1, il) = rxyz(1, i) - alat(1)
         rxyz(2, il) = rxyz(2, i) + alat(2)
         rxyz(3, il) = rxyz(3, i) + alat(3)
      END DO

      in = icell(0, 0, ll2 - 1, 0)
      icell(0, ll1, -1, ll3) = in
      DO ii = 1, in
         i = icell(ii, 0, ll2 - 1, 0)
         il = il + 1
         IF (il > nn) CPABORT("enlarge laymx")
         lay(il) = i
         icell(ii, ll1, -1, ll3) = il
         rxyz(1, il) = rxyz(1, i) + alat(1)
         rxyz(2, il) = rxyz(2, i) - alat(2)
         rxyz(3, il) = rxyz(3, i) + alat(3)
      END DO

      in = icell(0, ll1 - 1, ll2 - 1, 0)
      icell(0, -1, -1, ll3) = in
      DO ii = 1, in
         i = icell(ii, ll1 - 1, ll2 - 1, 0)
         il = il + 1
         IF (il > nn) CPABORT("enlarge laymx")
         lay(il) = i
         icell(ii, -1, -1, ll3) = il
         rxyz(1, il) = rxyz(1, i) - alat(1)
         rxyz(2, il) = rxyz(2, i) - alat(2)
         rxyz(3, il) = rxyz(3, i) + alat(3)
      END DO

      in = icell(0, 0, 0, ll3 - 1)
      icell(0, ll1, ll2, -1) = in
      DO ii = 1, in
         i = icell(ii, 0, 0, ll3 - 1)
         il = il + 1
         IF (il > nn) CPABORT("enlarge laymx")
         lay(il) = i
         icell(ii, ll1, ll2, -1) = il
         rxyz(1, il) = rxyz(1, i) + alat(1)
         rxyz(2, il) = rxyz(2, i) + alat(2)
         rxyz(3, il) = rxyz(3, i) - alat(3)
      END DO

      in = icell(0, ll1 - 1, 0, ll3 - 1)
      icell(0, -1, ll2, -1) = in
      DO ii = 1, in
         i = icell(ii, ll1 - 1, 0, ll3 - 1)
         il = il + 1
         IF (il > nn) CPABORT("enlarge laymx")
         lay(il) = i
         icell(ii, -1, ll2, -1) = il
         rxyz(1, il) = rxyz(1, i) - alat(1)
         rxyz(2, il) = rxyz(2, i) + alat(2)
         rxyz(3, il) = rxyz(3, i) - alat(3)
      END DO

      in = icell(0, 0, ll2 - 1, ll3 - 1)
      icell(0, ll1, -1, -1) = in
      DO ii = 1, in
         i = icell(ii, 0, ll2 - 1, ll3 - 1)
         il = il + 1
         IF (il > nn) CPABORT("enlarge laymx")
         lay(il) = i
         icell(ii, ll1, -1, -1) = il
         rxyz(1, il) = rxyz(1, i) + alat(1)
         rxyz(2, il) = rxyz(2, i) - alat(2)
         rxyz(3, il) = rxyz(3, i) - alat(3)
      END DO

      in = icell(0, ll1 - 1, ll2 - 1, ll3 - 1)
      icell(0, -1, -1, -1) = in
      DO ii = 1, in
         i = icell(ii, ll1 - 1, ll2 - 1, ll3 - 1)
         il = il + 1
         IF (il > nn) CPABORT("enlarge laymx")
         lay(il) = i
         icell(ii, -1, -1, -1) = il
         rxyz(1, il) = rxyz(1, i) - alat(1)
         rxyz(2, il) = rxyz(2, i) - alat(2)
         rxyz(3, il) = rxyz(3, i) - alat(3)
      END DO

      ALLOCATE (lsta(2, nat))
      nnbrx = 12
      loop_nnbrx: DO
         ALLOCATE (lstb(nnbrx*nat), rel(5, nnbrx*nat))

         indlstx = 0

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP PRIVATE(iat,cut2,iam,ii,indlst,l1,l2,l3,myspace,npr) &
!$OMP SHARED (indlstx,nat,nn,nnbrx,ncx,ll1,ll2,ll3,icell,lsta,lstb,lay, &
!$OMP rel,rxyz,cut,myspaceout)

         npr = 1
!$       npr = omp_get_num_threads()
         iam = 0
!$       iam = omp_get_thread_num()

         cut2 = cut**2
! assign contiguous portions of the arrays lstb and rel to the threads
         myspace = (nat*nnbrx)/npr
         IF (iam == 0) myspaceout = myspace
! Verlet list, relative positions
         indlst = 0
         loop_l3: DO l3 = 0, ll3 - 1
            loop_l2: DO l2 = 0, ll2 - 1
               loop_l1: DO l1 = 0, ll1 - 1
                  loop_ii: DO ii = 1, icell(0, l1, l2, l3)
                     iat = icell(ii, l1, l2, l3)
                     IF (((iat - 1)*npr)/nat == iam) THEN
!                          write(*,*) 'sublstiat:iam,iat',iam,iat
                        lsta(1, iat) = iam*myspace + indlst + 1
                        CALL sublstiat_b(iat, nn, ncx, ll1, ll2, ll3, l1, l2, l3, myspace, &
                                         rxyz, icell, lstb(iam*myspace + 1), lay, &
                                         rel(1, iam*myspace + 1), cut2, indlst)
                        lsta(2, iat) = iam*myspace + indlst
!                          write(*,'(a,4(x,i3),100(x,i2))') &
!                           'iam,iat,lsta',iam,iat,lsta(1,iat),lsta(2,iat), &
!                           (lstb(j),j=lsta(1,iat),lsta(2,iat))
                     END IF
                  END DO loop_ii
               END DO loop_l1
            END DO loop_l2
         END DO loop_l3
!$OMP ATOMIC UPDATE
         indlstx = MAX(indlstx, indlst)
!$OMP END ATOMIC
!$OMP END PARALLEL

         IF (indlstx >= myspaceout) THEN
            WRITE (10, *) count, 'NNBRX too small', nnbrx
            DEALLOCATE (lstb, rel)
            nnbrx = 3*nnbrx/2
            CYCLE loop_nnbrx
         END IF
         EXIT loop_nnbrx
      END DO loop_nnbrx

      istopg = 0

!$OMP PARALLEL DEFAULT(NONE)  &
!$OMP PRIVATE(iam,npr,iat,iat1,iat2,lot,istop,tcoord,tcoord2, &
!$OMP tener,tener2,txyz,s2,s3,sz,num2,num3,numz,max_nbrs) &
!$OMP SHARED (nat,nnbrx,lsta,lstb,rel,ener,ener2,fxyz,coord,coord2,istopg)

      npr = 1
!$    npr = omp_get_num_threads()
      iam = 0
!$    iam = omp_get_thread_num()

      max_nbrs = 30

      IF (npr /= 1) THEN
! PARALLEL CASE
! create temporary private scalars for reduction sum on energies and
!        temporary private array for reduction sum on forces
!$OMP CRITICAL(omp_eip_bazant_silicon)
         ALLOCATE (txyz(3, nat), s2(max_nbrs, 8), s3(max_nbrs, 7), sz(max_nbrs, 6), &
                   num2(max_nbrs), num3(max_nbrs), numz(max_nbrs))
!$OMP END CRITICAL(omp_eip_bazant_silicon)
         IF (iam == 0) THEN
            ener = 0.e0_dp
            ener2 = 0.e0_dp
            coord = 0.e0_dp
            coord2 = 0.e0_dp
         END IF
!$OMP DO
         DO iat = 1, nat
            fxyz(1, iat) = 0.e0_dp
            fxyz(2, iat) = 0.e0_dp
            fxyz(3, iat) = 0.e0_dp
         END DO
!$OMP BARRIER

! Each thread treats at most lot atoms
         lot = INT(REAL(nat, KIND=dp)/REAL(npr, KIND=dp) + .999999999999e0_dp)
         iat1 = iam*lot + 1
         iat2 = MIN((iam + 1)*lot, nat)
!       write(*,*) 'subfeniat:iat1,iat2,iam',iat1,iat2,iam
         CALL subfeniat_b(iat1, iat2, nat, lsta, lstb, rel, tener, tener2, &
                          tcoord, tcoord2, nnbrx, txyz, max_nbrs, istop, &
                          s2(1, 1), s2(1, 2), s2(1, 3), s2(1, 4), s2(1, 5), s2(1, 6), s2(1, 7), s2(1, 8), &
                          num2, s3(1, 1), s3(1, 2), s3(1, 3), s3(1, 4), s3(1, 5), s3(1, 6), s3(1, 7), &
                          num3, sz(1, 1), sz(1, 2), sz(1, 3), sz(1, 4), sz(1, 5), sz(1, 6), numz)

!$OMP CRITICAL(omp_eip_bazant_silicon)
         ener = ener + tener
         ener2 = ener2 + tener2
         coord = coord + tcoord
         coord2 = coord2 + tcoord2
         istopg = istopg + istop
         DO iat = 1, nat
            fxyz(1, iat) = fxyz(1, iat) + txyz(1, iat)
            fxyz(2, iat) = fxyz(2, iat) + txyz(2, iat)
            fxyz(3, iat) = fxyz(3, iat) + txyz(3, iat)
         END DO
         DEALLOCATE (txyz, s2, s3, sz, num2, num3, numz)
!$OMP END CRITICAL(omp_eip_bazant_silicon)

      ELSE
! SERIAL CASE
         iat1 = 1
         iat2 = nat
         ALLOCATE (s2(max_nbrs, 8), s3(max_nbrs, 7), sz(max_nbrs, 6), &
                   num2(max_nbrs), num3(max_nbrs), numz(max_nbrs))
         CALL subfeniat_b(iat1, iat2, nat, lsta, lstb, rel, ener, ener2, &
                          coord, coord2, nnbrx, fxyz, max_nbrs, istopg, &
                          s2(1, 1), s2(1, 2), s2(1, 3), s2(1, 4), s2(1, 5), s2(1, 6), s2(1, 7), s2(1, 8), &
                          num2, s3(1, 1), s3(1, 2), s3(1, 3), s3(1, 4), s3(1, 5), s3(1, 6), s3(1, 7), &
                          num3, sz(1, 1), sz(1, 2), sz(1, 3), sz(1, 4), sz(1, 5), sz(1, 6), numz)
         DEALLOCATE (s2, s3, sz, num2, num3, numz)

      END IF
!$OMP END PARALLEL

!         write(*,*) 'ener,norm force', &
!                    ener,DNRM2(3*nat,fxyz,1)
      IF (istopg > 0) CPABORT("DIMENSION ERROR (see WARNING above)")
      ener_var = ener2/nat - (ener/nat)**2
      coord = coord/nat
      coord_var = coord2/nat - coord**2

      DEALLOCATE (rxyz, icell, lay, lsta, lstb, rel)

   END SUBROUTINE eip_bazant_silicon

! **************************************************************************************************
!> \brief ...
!> \param iat1 ...
!> \param iat2 ...
!> \param nat ...
!> \param lsta ...
!> \param lstb ...
!> \param rel ...
!> \param ener ...
!> \param ener2 ...
!> \param coord ...
!> \param coord2 ...
!> \param nnbrx ...
!> \param ff ...
!> \param max_nbrs ...
!> \param istop ...
!> \param s2_t0 ...
!> \param s2_t1 ...
!> \param s2_t2 ...
!> \param s2_t3 ...
!> \param s2_dx ...
!> \param s2_dy ...
!> \param s2_dz ...
!> \param s2_r ...
!> \param num2 ...
!> \param s3_g ...
!> \param s3_dg ...
!> \param s3_rinv ...
!> \param s3_dx ...
!> \param s3_dy ...
!> \param s3_dz ...
!> \param s3_r ...
!> \param num3 ...
!> \param sz_df ...
!> \param sz_sum ...
!> \param sz_dx ...
!> \param sz_dy ...
!> \param sz_dz ...
!> \param sz_r ...
!> \param numz ...
! **************************************************************************************************
   SUBROUTINE subfeniat_b(iat1, iat2, nat, lsta, lstb, rel, ener, ener2, &
                          coord, coord2, nnbrx, ff, max_nbrs, istop, &
                          s2_t0, s2_t1, s2_t2, s2_t3, s2_dx, s2_dy, s2_dz, s2_r, &
                          num2, s3_g, s3_dg, s3_rinv, s3_dx, s3_dy, s3_dz, s3_r, &
                          num3, sz_df, sz_sum, sz_dx, sz_dy, sz_dz, sz_r, numz)
! This subroutine is a modification of a subroutine that is available at
! http://www-math.mit.edu/~bazant/EDIP/ and for which Martin Z. Bazant
! and Harvard University have a 1997 copyright.
! The modifications were done by S. Goedecker on April 10, 2002.
! The routines are included with the permission of M. Bazant into this package.

!  ------------------------- VARIABLE DECLARATIONS -------------------------
      INTEGER                                            :: iat1, iat2, nat, lsta(2, nat)
      REAL(KIND=dp)                                      :: ener, ener2, coord, coord2
      INTEGER                                            :: nnbrx
      REAL(KIND=dp)                                      :: rel(5, nnbrx*nat)
      INTEGER                                            :: lstb(nnbrx*nat)
      REAL(KIND=dp)                                      :: ff(3, nat)
      INTEGER                                            :: max_nbrs, istop
      REAL(KIND=dp) :: s2_t0(max_nbrs), s2_t1(max_nbrs), s2_t2(max_nbrs), s2_t3(max_nbrs), &
         s2_dx(max_nbrs), s2_dy(max_nbrs), s2_dz(max_nbrs), s2_r(max_nbrs)
      INTEGER                                            :: num2(max_nbrs)
      REAL(KIND=dp) :: s3_g(max_nbrs), s3_dg(max_nbrs), s3_rinv(max_nbrs), s3_dx(max_nbrs), &
         s3_dy(max_nbrs), s3_dz(max_nbrs), s3_r(max_nbrs)
      INTEGER                                            :: num3(max_nbrs)
      REAL(KIND=dp)                                      :: sz_df(max_nbrs), sz_sum(max_nbrs), &
                                                            sz_dx(max_nbrs), sz_dy(max_nbrs), &
                                                            sz_dz(max_nbrs), sz_r(max_nbrs)
      INTEGER                                            :: numz(max_nbrs)

      INTEGER                                            :: i, j, k, l, n, n2, n3, nj, nk, nl, nz
      REAL(KIND=dp) :: bmc, cmbinv, coord_iat, dEdrl, dEdrlx, dEdrly, dEdrlz, den, dhdl, dHdx, &
         dp1, dtau, dV2dZ, dV2ijx, dV2ijy, dV2ijz, dV2j, dV3dZ, dV3l, dV3ljx, dV3ljy, dV3ljz, &
         dV3lkx, dV3lky, dV3lkz, dV3rij, dV3rijx, dV3rijy, dV3rijz, dV3rik, dV3rikx, dV3riky, &
         dV3rikz, dwinv, dx, dxdZ, dy, dz, ener_iat, fjx, fjy, fjz, fkx, fky, fkz, fZ, H, lcos, &
         muhalf, par_a, par_alp, par_b, par_bet, par_bg, par_c, par_cap_A, par_cap_B, par_delta, &
         par_eta, par_gam, par_lam, par_mu, par_palp, par_Qo, par_rh, par_sig, pZ, Qort, r, rinv, &
         rmainv, rmbinv, tau, temp0, temp1, u1, u2, u3, u4, u5, winv, x, xarg
      REAL(KIND=dp) :: xinv, xinv3, Z

!   size of s2[]
!   atom ID numbers for s2[]
!   size of s3[]
!   atom ID numbers for s3[]
!   size of sz[]
!   atom ID numbers for sz[]
!   indices for the store arrays
!   EDIP parameters

      par_cap_A = 5.6714030e0_dp
      par_cap_B = 2.0002804e0_dp
      par_rh = 1.2085196e0_dp
      par_a = 3.1213820e0_dp
      par_sig = 0.5774108e0_dp
      par_lam = 1.4533108e0_dp
      par_gam = 1.1247945e0_dp
      par_b = 3.1213820e0_dp
      par_c = 2.5609104e0_dp
      par_delta = 78.7590539e0_dp
      par_mu = 0.6966326e0_dp
      par_Qo = 312.1341346e0_dp
      par_palp = 1.4074424e0_dp
      par_bet = 0.0070975e0_dp
      par_alp = 3.1083847e0_dp

      u1 = -0.165799e0_dp
      u2 = 32.557e0_dp
      u3 = 0.286198e0_dp
      u4 = 0.66e0_dp

      par_bg = par_a
      par_eta = par_delta/par_Qo

      DO i = 1, nat
         ff(1, i) = 0.0e0_dp
         ff(2, i) = 0.0e0_dp
         ff(3, i) = 0.0e0_dp
      END DO

      coord = 0.e0_dp
      coord2 = 0.e0_dp
      ener = 0.e0_dp
      ener2 = 0.e0_dp
      istop = 0

!   COMBINE COEFFICIENTS

      Qort = SQRT(par_Qo)
      muhalf = par_mu*0.5e0_dp
      u5 = u2*u4
      bmc = par_b - par_c
      cmbinv = 1.0e0_dp/(par_c - par_b)

!  --- LEVEL 1: OUTER LOOP OVER ATOMS ---

      atoms: DO i = iat1, iat2

!   RESET COORDINATION AND NEIGHBOR NUMBERS

         coord_iat = 0.e0_dp
         ener_iat = 0.e0_dp
         Z = 0.0e0_dp
         n2 = 1
         n3 = 1
         nz = 1

!  --- LEVEL 2: LOOP PREPASS OVER PAIRS ---

         DO n = lsta(1, i), lsta(2, i)
            j = lstb(n)

!   PARTS OF TWO-BODY INTERACTION r<par_a

            num2(n2) = j
            dx = -rel(1, n)
            dy = -rel(2, n)
            dz = -rel(3, n)
            r = rel(4, n)
            rinv = rel(5, n)
            rmainv = 1.e0_dp/(r - par_a)
            s2_t0(n2) = par_cap_A*EXP(par_sig*rmainv)
            s2_t1(n2) = (par_cap_B*rinv)**par_rh
            s2_t2(n2) = par_rh*rinv
            s2_t3(n2) = par_sig*rmainv*rmainv
            s2_dx(n2) = dx
            s2_dy(n2) = dy
            s2_dz(n2) = dz
            s2_r(n2) = r
            n2 = n2 + 1
            IF (n2 > max_nbrs) THEN
               WRITE (*, *) 'WARNING enlarge max_nbrs'
               istop = 1
               RETURN
            END IF

! coordination number calculated with soft cutoff between first
! nearest neighbor and midpoint of first and second nearest neighbor
            IF (r <= 2.36e0_dp) THEN
               coord_iat = coord_iat + 1.e0_dp
            ELSE IF (r >= 3.12e0_dp) THEN
            ELSE
               xarg = (r - 2.36e0_dp)*(1.e0_dp/(3.12e0_dp - 2.36e0_dp))
               coord_iat = coord_iat + (2*xarg + 1.e0_dp)*(xarg - 1.e0_dp)**2
            END IF

!   RADIAL PARTS OF THREE-BODY INTERACTION r<par_b

            IF (r < par_bg) THEN

               num3(n3) = j
               rmbinv = 1.e0_dp/(r - par_bg)
               temp1 = par_gam*rmbinv
               temp0 = EXP(temp1)
               s3_g(n3) = temp0
               s3_dg(n3) = -rmbinv*temp1*temp0
               s3_dx(n3) = dx
               s3_dy(n3) = dy
               s3_dz(n3) = dz
               s3_rinv(n3) = rinv
               s3_r(n3) = r
               n3 = n3 + 1
               IF (n3 > max_nbrs) THEN
                  WRITE (*, *) 'WARNING enlarge max_nbrs'
                  istop = 1
                  RETURN
               END IF

!   COORDINATION AND NEIGHBOR FUNCTION par_c<r<par_b

               IF (r < par_b) THEN
                  IF (r < par_c) THEN
                     Z = Z + 1.e0_dp
                  ELSE
                     xinv = bmc/(r - par_c)
                     xinv3 = xinv*xinv*xinv
                     den = 1.e0_dp/(1 - xinv3)
                     temp1 = par_alp*den
                     fZ = EXP(temp1)
                     Z = Z + fZ
                     numz(nz) = j
                     sz_df(nz) = fZ*temp1*den*3.e0_dp*xinv3*xinv*cmbinv
!   df/dr
                     sz_dx(nz) = dx
                     sz_dy(nz) = dy
                     sz_dz(nz) = dz
                     sz_r(nz) = r
                     nz = nz + 1
                     IF (nz > max_nbrs) THEN
                        WRITE (*, *) 'WARNING enlarge max_nbrs'
                        istop = 1
                        RETURN
                     END IF
                  END IF
!  r < par_C
               END IF
!  r < par_b
            END IF
!  r < par_bg
         END DO

!   ZERO ACCUMULATION ARRAY FOR ENVIRONMENT FORCES

         DO nl = 1, nz - 1
            sz_sum(nl) = 0.e0_dp
         END DO

!   ENVIRONMENT-DEPENDENCE OF PAIR INTERACTION

         temp0 = par_bet*Z
         pZ = par_palp*EXP(-temp0*Z)
!   bond order
         dp1 = -2.e0_dp*temp0*pZ
!   derivative of bond order

!  --- LEVEL 2: LOOP FOR PAIR INTERACTIONS ---

         DO nj = 1, n2 - 1

            temp0 = s2_t1(nj) - pZ

!   two-body energy V2(rij,Z)

            ener_iat = ener_iat + temp0*s2_t0(nj)

!   two-body forces

            dV2j = -s2_t0(nj)*(s2_t1(nj)*s2_t2(nj) + temp0*s2_t3(nj))
!   dV2/dr
            dV2ijx = dV2j*s2_dx(nj)
            dV2ijy = dV2j*s2_dy(nj)
            dV2ijz = dV2j*s2_dz(nj)
            ff(1, i) = ff(1, i) + dV2ijx
            ff(2, i) = ff(2, i) + dV2ijy
            ff(3, i) = ff(3, i) + dV2ijz
            j = num2(nj)
            ff(1, j) = ff(1, j) - dV2ijx
            ff(2, j) = ff(2, j) - dV2ijy
            ff(3, j) = ff(3, j) - dV2ijz

!  --- LEVEL 3: LOOP FOR PAIR COORDINATION FORCES ---

            dV2dZ = -dp1*s2_t0(nj)
            DO nl = 1, nz - 1
               sz_sum(nl) = sz_sum(nl) + dV2dZ
            END DO

         END DO

!   COORDINATION-DEPENDENCE OF THREE-BODY INTERACTION

         winv = Qort*EXP(-muhalf*Z)
!   inverse width of angular function
         dwinv = -muhalf*winv
!   its derivative
         temp0 = EXP(-u4*Z)
         tau = u1 + u2*temp0*(u3 - temp0)
!   -cosine of angular minimum
         dtau = u5*temp0*(2*temp0 - u3)
!   its derivative

!  --- LEVEL 2: FIRST LOOP FOR THREE-BODY INTERACTIONS ---

         DO nj = 1, n3 - 2

            j = num3(nj)

!  --- LEVEL 3: SECOND LOOP FOR THREE-BODY INTERACTIONS ---

            DO nk = nj + 1, n3 - 1

               k = num3(nk)

!   angular function h(l,Z)

               lcos = s3_dx(nj)*s3_dx(nk) + s3_dy(nj)*s3_dy(nk) + s3_dz(nj)*s3_dz(nk)
               x = (lcos + tau)*winv
               temp0 = EXP(-x*x)

               H = par_lam*(1 - temp0 + par_eta*x*x)
               dHdx = 2*par_lam*x*(temp0 + par_eta)

               dhdl = dHdx*winv

!   three-body energy

               temp1 = s3_g(nj)*s3_g(nk)
               ener_iat = ener_iat + temp1*H

!   (-) radial force on atom j

               dV3rij = s3_dg(nj)*s3_g(nk)*H
               dV3rijx = dV3rij*s3_dx(nj)
               dV3rijy = dV3rij*s3_dy(nj)
               dV3rijz = dV3rij*s3_dz(nj)
               fjx = dV3rijx
               fjy = dV3rijy
               fjz = dV3rijz

!   (-) radial force on atom k

               dV3rik = s3_g(nj)*s3_dg(nk)*H
               dV3rikx = dV3rik*s3_dx(nk)
               dV3riky = dV3rik*s3_dy(nk)
               dV3rikz = dV3rik*s3_dz(nk)
               fkx = dV3rikx
               fky = dV3riky
               fkz = dV3rikz

!   (-) angular force on j

               dV3l = temp1*dhdl
               dV3ljx = dV3l*(s3_dx(nk) - lcos*s3_dx(nj))*s3_rinv(nj)
               dV3ljy = dV3l*(s3_dy(nk) - lcos*s3_dy(nj))*s3_rinv(nj)
               dV3ljz = dV3l*(s3_dz(nk) - lcos*s3_dz(nj))*s3_rinv(nj)
               fjx = fjx + dV3ljx
               fjy = fjy + dV3ljy
               fjz = fjz + dV3ljz

!   (-) angular force on k

               dV3lkx = dV3l*(s3_dx(nj) - lcos*s3_dx(nk))*s3_rinv(nk)
               dV3lky = dV3l*(s3_dy(nj) - lcos*s3_dy(nk))*s3_rinv(nk)
               dV3lkz = dV3l*(s3_dz(nj) - lcos*s3_dz(nk))*s3_rinv(nk)
               fkx = fkx + dV3lkx
               fky = fky + dV3lky
               fkz = fkz + dV3lkz

!   apply radial + angular forces to i, j, k

               ff(1, j) = ff(1, j) - fjx
               ff(2, j) = ff(2, j) - fjy
               ff(3, j) = ff(3, j) - fjz
               ff(1, k) = ff(1, k) - fkx
               ff(2, k) = ff(2, k) - fky
               ff(3, k) = ff(3, k) - fkz
               ff(1, i) = ff(1, i) + fjx + fkx
               ff(2, i) = ff(2, i) + fjy + fky
               ff(3, i) = ff(3, i) + fjz + fkz

!   prefactor for 4-body forces from coordination
               dxdZ = dwinv*(lcos + tau) + winv*dtau
               dV3dZ = temp1*dHdx*dxdZ

!  --- LEVEL 4: LOOP FOR THREE-BODY COORDINATION FORCES ---

               DO nl = 1, nz - 1
                  sz_sum(nl) = sz_sum(nl) + dV3dZ
               END DO
            END DO
         END DO

!  --- LEVEL 2: LOOP TO APPLY COORDINATION FORCES ---

         DO nl = 1, nz - 1

            dEdrl = sz_sum(nl)*sz_df(nl)
            dEdrlx = dEdrl*sz_dx(nl)
            dEdrly = dEdrl*sz_dy(nl)
            dEdrlz = dEdrl*sz_dz(nl)
            ff(1, i) = ff(1, i) + dEdrlx
            ff(2, i) = ff(2, i) + dEdrly
            ff(3, i) = ff(3, i) + dEdrlz
            l = numz(nl)
            ff(1, l) = ff(1, l) - dEdrlx
            ff(2, l) = ff(2, l) - dEdrly
            ff(3, l) = ff(3, l) - dEdrlz

         END DO

         coord = coord + coord_iat
         coord2 = coord2 + coord_iat**2
         ener = ener + ener_iat
         ener2 = ener2 + ener_iat**2

      END DO atoms

      RETURN
   END SUBROUTINE subfeniat_b

! **************************************************************************************************
!> \brief ...
!> \param iat ...
!> \param nn ...
!> \param ncx ...
!> \param ll1 ...
!> \param ll2 ...
!> \param ll3 ...
!> \param l1 ...
!> \param l2 ...
!> \param l3 ...
!> \param myspace ...
!> \param rxyz ...
!> \param icell ...
!> \param lstb ...
!> \param lay ...
!> \param rel ...
!> \param cut2 ...
!> \param indlst ...
! **************************************************************************************************
   SUBROUTINE sublstiat_b(iat, nn, ncx, ll1, ll2, ll3, l1, l2, l3, myspace, &
                          rxyz, icell, lstb, lay, rel, cut2, indlst)
! finds the neighbours of atom iat (specified by lsta and lstb) and and
! the relative position rel of iat with respect to these neighbours
      INTEGER                                            :: iat, nn, ncx, ll1, ll2, ll3, l1, l2, l3, &
                                                            myspace
      REAL(KIND=dp)                                      :: rxyz
      INTEGER                                            :: icell, lstb, lay
      REAL(KIND=dp)                                      :: rel, cut2
      INTEGER                                            :: indlst

      DIMENSION rxyz(3, nn), lay(nn), icell(0:ncx, -1:ll1, -1:ll2, -1:ll3), &
         lstb(0:myspace - 1), rel(5, 0:myspace - 1)

      INTEGER       :: jat, k1, k2, k3, jj
      REAL(KIND=dp) :: rr2, tt, tti, xrel, yrel, zrel

      DO k3 = l3 - 1, l3 + 1
      DO k2 = l2 - 1, l2 + 1
      DO k1 = l1 - 1, l1 + 1
      DO jj = 1, icell(0, k1, k2, k3)
         jat = icell(jj, k1, k2, k3)
         IF (jat == iat) CYCLE
         xrel = rxyz(1, iat) - rxyz(1, jat)
         yrel = rxyz(2, iat) - rxyz(2, jat)
         zrel = rxyz(3, iat) - rxyz(3, jat)
         rr2 = xrel**2 + yrel**2 + zrel**2
         IF (rr2 <= cut2) THEN
            indlst = MIN(indlst, myspace - 1)
            lstb(indlst) = lay(jat)
!        write(*,*) 'iat,indlst,lay(jat)',iat,indlst,lay(jat)
            tt = SQRT(rr2)
            tti = 1.e0_dp/tt
            rel(1, indlst) = xrel*tti
            rel(2, indlst) = yrel*tti
            rel(3, indlst) = zrel*tti
            rel(4, indlst) = tt
            rel(5, indlst) = tti
            indlst = indlst + 1
         END IF
      END DO
      END DO
      END DO
      END DO

      RETURN
   END SUBROUTINE sublstiat_b

! **************************************************************************************************
!> \brief Lenosky's "highly optimized empirical potential model of silicon"
!>      by Stefan Goedecker
!> \param nat number of atoms
!> \param alat lattice constants of the orthorombic box containing the particles
!> \param rxyz0 atomic positions in Angstrom, may be modified on output.
!>               If an atom is outside the box the program will bring it back
!>               into the box by translations through alat
!> \param fxyz forces in eV/A
!> \param ener total energy in eV
!> \param coord average coordination number
!> \param ener_var variance of the energy/atom
!> \param coord_var variance of the coordination number
!> \param count count is increased by one per call, has to be initialized
!>                to 0.e0_dp before first call of eip_bazant
!> \par Literature
!>      T. Lenosky, et. al.: Highly optimized empirical potential model of silicon;
!>                           Modeling Simul. Sci. Eng., 8 (2000)
!>      S. Goedecker: Optimization and parallelization of a force field for silicon
!>                    using OpenMP; CPC 148, 1 (2002)
!> \par History
!>      03.2006 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! **************************************************************************************************
   SUBROUTINE eip_lenosky_silicon(nat, alat, rxyz0, fxyz, ener, coord, ener_var, &
                                  coord_var, count)

      INTEGER                                            :: nat
      REAL(KIND=dp)                                      :: alat, rxyz0, fxyz, ener, coord, &
                                                            ener_var, coord_var, count

      DIMENSION rxyz0(3, nat), fxyz(3, nat), alat(3)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :) :: rxyz
      INTEGER, ALLOCATABLE, DIMENSION(:, :)       :: lsta
      INTEGER, ALLOCATABLE, DIMENSION(:)         :: lstb
      INTEGER, ALLOCATABLE, DIMENSION(:)         :: lay
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :, :)   :: icell
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :) :: rel
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :) :: txyz
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :) :: f2ij, f3ij, f3ik

      REAL(KIND=dp) :: coord2, cut, cut2, ener2, tcoord, &
                       tcoord2, tener, tener2
      INTEGER       :: i, iam, iat, iat1, iat2, ii, il, in, indlst, indlstx, &
                       istop, istopg, l1, l2, l3, ll1, ll2, ll3, lot, ncx, nn, &
                       nnbrx, npjkx, npjx, laymx, npr, rlc1i, rlc2i, rlc3i, &
                       myspace, myspaceout

!        tmax_phi= 0.4500000e+01_dp
!        cut=tmax_phi
      cut = 0.4500000e+01_dp

      IF (count == 0) OPEN (unit=10, file='lenosky.mon', status='unknown')
      count = count + 1.e0_dp

! linear scaling calculation of verlet list
      ll1 = INT(alat(1)/cut)
      IF (ll1 < 1) CPABORT("alat(1) too small")
      ll2 = INT(alat(2)/cut)
      IF (ll2 < 1) CPABORT("alat(2) too small")
      ll3 = INT(alat(3)/cut)
      IF (ll3 < 1) CPABORT("alat(3) too small")

! determine number of threads
      npr = 1
!$OMP PARALLEL PRIVATE(iam)  SHARED (npr) DEFAULT(NONE)
!$    iam = omp_get_thread_num()
!$    if (iam .eq. 0) npr = omp_get_num_threads()
!$OMP END PARALLEL

! linear scaling calculation of verlet list

      IF (npr <= 1) THEN !serial if too few processors to gain by parallelizing

! set ncx for serial case, ncx for parallel case set below
         ncx = 16
         loop_ncx_s: DO
            ALLOCATE (icell(0:ncx, -1:ll1, -1:ll2, -1:ll3))
            icell(0, -1:ll1, -1:ll2, -1:ll3) = 0
            rlc1i = INT(ll1/alat(1))
            rlc2i = INT(ll2/alat(2))
            rlc3i = INT(ll3/alat(3))

            loop_iat_s: DO iat = 1, nat
               rxyz0(1, iat) = MODULO(MODULO(rxyz0(1, iat), alat(1)), alat(1))
               rxyz0(2, iat) = MODULO(MODULO(rxyz0(2, iat), alat(2)), alat(2))
               rxyz0(3, iat) = MODULO(MODULO(rxyz0(3, iat), alat(3)), alat(3))
               l1 = INT(rxyz0(1, iat)*rlc1i)
               l2 = INT(rxyz0(2, iat)*rlc2i)
               l3 = INT(rxyz0(3, iat)*rlc3i)

               ii = icell(0, l1, l2, l3)
               ii = ii + 1
               icell(0, l1, l2, l3) = ii
               IF (ii > ncx) THEN
                  WRITE (10, *) count, 'NCX too small', ncx
                  DEALLOCATE (icell)
                  ncx = ncx*2
                  CYCLE loop_ncx_s
               END IF
               icell(ii, l1, l2, l3) = iat
            END DO loop_iat_s
            EXIT loop_ncx_s
         END DO loop_ncx_s

      ELSE ! parallel case

! periodization of particles can be done in parallel
!$OMP PARALLEL DO SHARED (alat,nat,rxyz0) PRIVATE(iat) DEFAULT(NONE)
         DO iat = 1, nat
            rxyz0(1, iat) = MODULO(MODULO(rxyz0(1, iat), alat(1)), alat(1))
            rxyz0(2, iat) = MODULO(MODULO(rxyz0(2, iat), alat(2)), alat(2))
            rxyz0(3, iat) = MODULO(MODULO(rxyz0(3, iat), alat(3)), alat(3))
         END DO
!$OMP END PARALLEL DO

! assignment to cell is done serially
! set ncx for parallel case, ncx for serial case set above
         ncx = 16
         loop_ncx_p: DO
            ALLOCATE (icell(0:ncx, -1:ll1, -1:ll2, -1:ll3))
            icell(0, -1:ll1, -1:ll2, -1:ll3) = 0
            rlc1i = INT(ll1/alat(1))
            rlc2i = INT(ll2/alat(2))
            rlc3i = INT(ll3/alat(3))

            loop_iat_p: DO iat = 1, nat
               l1 = INT(rxyz0(1, iat)*rlc1i)
               l2 = INT(rxyz0(2, iat)*rlc2i)
               l3 = INT(rxyz0(3, iat)*rlc3i)
               ii = icell(0, l1, l2, l3)
               ii = ii + 1
               icell(0, l1, l2, l3) = ii
               IF (ii > ncx) THEN
                  WRITE (10, *) count, 'NCX too small', ncx
                  DEALLOCATE (icell)
                  ncx = ncx*2
                  CYCLE loop_ncx_p
               END IF
               icell(ii, l1, l2, l3) = iat
            END DO loop_iat_p
            EXIT loop_ncx_p
         END DO loop_ncx_p

      END IF

! duplicate all atoms within boundary layer
      laymx = ncx*(2*ll1*ll2 + 2*ll1*ll3 + 2*ll2*ll3 + 4*ll1 + 4*ll2 + 4*ll3 + 8)
      nn = nat + laymx
      ALLOCATE (rxyz(3, nn), lay(nn))
      DO iat = 1, nat
         lay(iat) = iat
         rxyz(1, iat) = rxyz0(1, iat)
         rxyz(2, iat) = rxyz0(2, iat)
         rxyz(3, iat) = rxyz0(3, iat)
      END DO
      il = nat
! xy plane
      DO l2 = 0, ll2 - 1
      DO l1 = 0, ll1 - 1

         in = icell(0, l1, l2, 0)
         icell(0, l1, l2, ll3) = in
         DO ii = 1, in
            i = icell(ii, l1, l2, 0)
            il = il + 1
            IF (il > nn) CPABORT("enlarge laymx")
            lay(il) = i
            icell(ii, l1, l2, ll3) = il
            rxyz(1, il) = rxyz(1, i)
            rxyz(2, il) = rxyz(2, i)
            rxyz(3, il) = rxyz(3, i) + alat(3)
         END DO

         in = icell(0, l1, l2, ll3 - 1)
         icell(0, l1, l2, -1) = in
         DO ii = 1, in
            i = icell(ii, l1, l2, ll3 - 1)
            il = il + 1
            IF (il > nn) CPABORT("enlarge laymx")
            lay(il) = i
            icell(ii, l1, l2, -1) = il
            rxyz(1, il) = rxyz(1, i)
            rxyz(2, il) = rxyz(2, i)
            rxyz(3, il) = rxyz(3, i) - alat(3)
         END DO

      END DO
      END DO

! yz plane
      DO l3 = 0, ll3 - 1
      DO l2 = 0, ll2 - 1

         in = icell(0, 0, l2, l3)
         icell(0, ll1, l2, l3) = in
         DO ii = 1, in
            i = icell(ii, 0, l2, l3)
            il = il + 1
            IF (il > nn) CPABORT("enlarge laymx")
            lay(il) = i
            icell(ii, ll1, l2, l3) = il
            rxyz(1, il) = rxyz(1, i) + alat(1)
            rxyz(2, il) = rxyz(2, i)
            rxyz(3, il) = rxyz(3, i)
         END DO

         in = icell(0, ll1 - 1, l2, l3)
         icell(0, -1, l2, l3) = in
         DO ii = 1, in
            i = icell(ii, ll1 - 1, l2, l3)
            il = il + 1
            IF (il > nn) CPABORT("enlarge laymx")
            lay(il) = i
            icell(ii, -1, l2, l3) = il
            rxyz(1, il) = rxyz(1, i) - alat(1)
            rxyz(2, il) = rxyz(2, i)
            rxyz(3, il) = rxyz(3, i)
         END DO

      END DO
      END DO

! xz plane
      DO l3 = 0, ll3 - 1
      DO l1 = 0, ll1 - 1

         in = icell(0, l1, 0, l3)
         icell(0, l1, ll2, l3) = in
         DO ii = 1, in
            i = icell(ii, l1, 0, l3)
            il = il + 1
            IF (il > nn) CPABORT("enlarge laymx")
            lay(il) = i
            icell(ii, l1, ll2, l3) = il
            rxyz(1, il) = rxyz(1, i)
            rxyz(2, il) = rxyz(2, i) + alat(2)
            rxyz(3, il) = rxyz(3, i)
         END DO

         in = icell(0, l1, ll2 - 1, l3)
         icell(0, l1, -1, l3) = in
         DO ii = 1, in
            i = icell(ii, l1, ll2 - 1, l3)
            il = il + 1
            IF (il > nn) CPABORT("enlarge laymx")
            lay(il) = i
            icell(ii, l1, -1, l3) = il
            rxyz(1, il) = rxyz(1, i)
            rxyz(2, il) = rxyz(2, i) - alat(2)
            rxyz(3, il) = rxyz(3, i)
         END DO

      END DO
      END DO

! x axis
      DO l1 = 0, ll1 - 1

         in = icell(0, l1, 0, 0)
         icell(0, l1, ll2, ll3) = in
         DO ii = 1, in
            i = icell(ii, l1, 0, 0)
            il = il + 1
            IF (il > nn) CPABORT("enlarge laymx")
            lay(il) = i
            icell(ii, l1, ll2, ll3) = il
            rxyz(1, il) = rxyz(1, i)
            rxyz(2, il) = rxyz(2, i) + alat(2)
            rxyz(3, il) = rxyz(3, i) + alat(3)
         END DO

         in = icell(0, l1, 0, ll3 - 1)
         icell(0, l1, ll2, -1) = in
         DO ii = 1, in
            i = icell(ii, l1, 0, ll3 - 1)
            il = il + 1
            IF (il > nn) CPABORT("enlarge laymx")
            lay(il) = i
            icell(ii, l1, ll2, -1) = il
            rxyz(1, il) = rxyz(1, i)
            rxyz(2, il) = rxyz(2, i) + alat(2)
            rxyz(3, il) = rxyz(3, i) - alat(3)
         END DO

         in = icell(0, l1, ll2 - 1, 0)
         icell(0, l1, -1, ll3) = in
         DO ii = 1, in
            i = icell(ii, l1, ll2 - 1, 0)
            il = il + 1
            IF (il > nn) CPABORT("enlarge laymx")
            lay(il) = i
            icell(ii, l1, -1, ll3) = il
            rxyz(1, il) = rxyz(1, i)
            rxyz(2, il) = rxyz(2, i) - alat(2)
            rxyz(3, il) = rxyz(3, i) + alat(3)
         END DO

         in = icell(0, l1, ll2 - 1, ll3 - 1)
         icell(0, l1, -1, -1) = in
         DO ii = 1, in
            i = icell(ii, l1, ll2 - 1, ll3 - 1)
            il = il + 1
            IF (il > nn) CPABORT("enlarge laymx")
            lay(il) = i
            icell(ii, l1, -1, -1) = il
            rxyz(1, il) = rxyz(1, i)
            rxyz(2, il) = rxyz(2, i) - alat(2)
            rxyz(3, il) = rxyz(3, i) - alat(3)
         END DO

      END DO

! y axis
      DO l2 = 0, ll2 - 1

         in = icell(0, 0, l2, 0)
         icell(0, ll1, l2, ll3) = in
         DO ii = 1, in
            i = icell(ii, 0, l2, 0)
            il = il + 1
            IF (il > nn) CPABORT("enlarge laymx")
            lay(il) = i
            icell(ii, ll1, l2, ll3) = il
            rxyz(1, il) = rxyz(1, i) + alat(1)
            rxyz(2, il) = rxyz(2, i)
            rxyz(3, il) = rxyz(3, i) + alat(3)
         END DO

         in = icell(0, 0, l2, ll3 - 1)
         icell(0, ll1, l2, -1) = in
         DO ii = 1, in
            i = icell(ii, 0, l2, ll3 - 1)
            il = il + 1
            IF (il > nn) CPABORT("enlarge laymx")
            lay(il) = i
            icell(ii, ll1, l2, -1) = il
            rxyz(1, il) = rxyz(1, i) + alat(1)
            rxyz(2, il) = rxyz(2, i)
            rxyz(3, il) = rxyz(3, i) - alat(3)
         END DO

         in = icell(0, ll1 - 1, l2, 0)
         icell(0, -1, l2, ll3) = in
         DO ii = 1, in
            i = icell(ii, ll1 - 1, l2, 0)
            il = il + 1
            IF (il > nn) CPABORT("enlarge laymx")
            lay(il) = i
            icell(ii, -1, l2, ll3) = il
            rxyz(1, il) = rxyz(1, i) - alat(1)
            rxyz(2, il) = rxyz(2, i)
            rxyz(3, il) = rxyz(3, i) + alat(3)
         END DO

         in = icell(0, ll1 - 1, l2, ll3 - 1)
         icell(0, -1, l2, -1) = in
         DO ii = 1, in
            i = icell(ii, ll1 - 1, l2, ll3 - 1)
            il = il + 1
            IF (il > nn) CPABORT("enlarge laymx")
            lay(il) = i
            icell(ii, -1, l2, -1) = il
            rxyz(1, il) = rxyz(1, i) - alat(1)
            rxyz(2, il) = rxyz(2, i)
            rxyz(3, il) = rxyz(3, i) - alat(3)
         END DO

      END DO

! z axis
      DO l3 = 0, ll3 - 1

         in = icell(0, 0, 0, l3)
         icell(0, ll1, ll2, l3) = in
         DO ii = 1, in
            i = icell(ii, 0, 0, l3)
            il = il + 1
            IF (il > nn) CPABORT("enlarge laymx")
            lay(il) = i
            icell(ii, ll1, ll2, l3) = il
            rxyz(1, il) = rxyz(1, i) + alat(1)
            rxyz(2, il) = rxyz(2, i) + alat(2)
            rxyz(3, il) = rxyz(3, i)
         END DO

         in = icell(0, ll1 - 1, 0, l3)
         icell(0, -1, ll2, l3) = in
         DO ii = 1, in
            i = icell(ii, ll1 - 1, 0, l3)
            il = il + 1
            IF (il > nn) CPABORT("enlarge laymx")
            lay(il) = i
            icell(ii, -1, ll2, l3) = il
            rxyz(1, il) = rxyz(1, i) - alat(1)
            rxyz(2, il) = rxyz(2, i) + alat(2)
            rxyz(3, il) = rxyz(3, i)
         END DO

         in = icell(0, 0, ll2 - 1, l3)
         icell(0, ll1, -1, l3) = in
         DO ii = 1, in
            i = icell(ii, 0, ll2 - 1, l3)
            il = il + 1
            IF (il > nn) CPABORT("enlarge laymx")
            lay(il) = i
            icell(ii, ll1, -1, l3) = il
            rxyz(1, il) = rxyz(1, i) + alat(1)
            rxyz(2, il) = rxyz(2, i) - alat(2)
            rxyz(3, il) = rxyz(3, i)
         END DO

         in = icell(0, ll1 - 1, ll2 - 1, l3)
         icell(0, -1, -1, l3) = in
         DO ii = 1, in
            i = icell(ii, ll1 - 1, ll2 - 1, l3)
            il = il + 1
            IF (il > nn) CPABORT("enlarge laymx")
            lay(il) = i
            icell(ii, -1, -1, l3) = il
            rxyz(1, il) = rxyz(1, i) - alat(1)
            rxyz(2, il) = rxyz(2, i) - alat(2)
            rxyz(3, il) = rxyz(3, i)
         END DO

      END DO

! corners
      in = icell(0, 0, 0, 0)
      icell(0, ll1, ll2, ll3) = in
      DO ii = 1, in
         i = icell(ii, 0, 0, 0)
         il = il + 1
         IF (il > nn) CPABORT("enlarge laymx")
         lay(il) = i
         icell(ii, ll1, ll2, ll3) = il
         rxyz(1, il) = rxyz(1, i) + alat(1)
         rxyz(2, il) = rxyz(2, i) + alat(2)
         rxyz(3, il) = rxyz(3, i) + alat(3)
      END DO

      in = icell(0, ll1 - 1, 0, 0)
      icell(0, -1, ll2, ll3) = in
      DO ii = 1, in
         i = icell(ii, ll1 - 1, 0, 0)
         il = il + 1
         IF (il > nn) CPABORT("enlarge laymx")
         lay(il) = i
         icell(ii, -1, ll2, ll3) = il
         rxyz(1, il) = rxyz(1, i) - alat(1)
         rxyz(2, il) = rxyz(2, i) + alat(2)
         rxyz(3, il) = rxyz(3, i) + alat(3)
      END DO

      in = icell(0, 0, ll2 - 1, 0)
      icell(0, ll1, -1, ll3) = in
      DO ii = 1, in
         i = icell(ii, 0, ll2 - 1, 0)
         il = il + 1
         IF (il > nn) CPABORT("enlarge laymx")
         lay(il) = i
         icell(ii, ll1, -1, ll3) = il
         rxyz(1, il) = rxyz(1, i) + alat(1)
         rxyz(2, il) = rxyz(2, i) - alat(2)
         rxyz(3, il) = rxyz(3, i) + alat(3)
      END DO

      in = icell(0, ll1 - 1, ll2 - 1, 0)
      icell(0, -1, -1, ll3) = in
      DO ii = 1, in
         i = icell(ii, ll1 - 1, ll2 - 1, 0)
         il = il + 1
         IF (il > nn) CPABORT("enlarge laymx")
         lay(il) = i
         icell(ii, -1, -1, ll3) = il
         rxyz(1, il) = rxyz(1, i) - alat(1)
         rxyz(2, il) = rxyz(2, i) - alat(2)
         rxyz(3, il) = rxyz(3, i) + alat(3)
      END DO

      in = icell(0, 0, 0, ll3 - 1)
      icell(0, ll1, ll2, -1) = in
      DO ii = 1, in
         i = icell(ii, 0, 0, ll3 - 1)
         il = il + 1
         IF (il > nn) CPABORT("enlarge laymx")
         lay(il) = i
         icell(ii, ll1, ll2, -1) = il
         rxyz(1, il) = rxyz(1, i) + alat(1)
         rxyz(2, il) = rxyz(2, i) + alat(2)
         rxyz(3, il) = rxyz(3, i) - alat(3)
      END DO

      in = icell(0, ll1 - 1, 0, ll3 - 1)
      icell(0, -1, ll2, -1) = in
      DO ii = 1, in
         i = icell(ii, ll1 - 1, 0, ll3 - 1)
         il = il + 1
         IF (il > nn) CPABORT("enlarge laymx")
         lay(il) = i
         icell(ii, -1, ll2, -1) = il
         rxyz(1, il) = rxyz(1, i) - alat(1)
         rxyz(2, il) = rxyz(2, i) + alat(2)
         rxyz(3, il) = rxyz(3, i) - alat(3)
      END DO

      in = icell(0, 0, ll2 - 1, ll3 - 1)
      icell(0, ll1, -1, -1) = in
      DO ii = 1, in
         i = icell(ii, 0, ll2 - 1, ll3 - 1)
         il = il + 1
         IF (il > nn) CPABORT("enlarge laymx")
         lay(il) = i
         icell(ii, ll1, -1, -1) = il
         rxyz(1, il) = rxyz(1, i) + alat(1)
         rxyz(2, il) = rxyz(2, i) - alat(2)
         rxyz(3, il) = rxyz(3, i) - alat(3)
      END DO

      in = icell(0, ll1 - 1, ll2 - 1, ll3 - 1)
      icell(0, -1, -1, -1) = in
      DO ii = 1, in
         i = icell(ii, ll1 - 1, ll2 - 1, ll3 - 1)
         il = il + 1
         IF (il > nn) CPABORT("enlarge laymx")
         lay(il) = i
         icell(ii, -1, -1, -1) = il
         rxyz(1, il) = rxyz(1, i) - alat(1)
         rxyz(2, il) = rxyz(2, i) - alat(2)
         rxyz(3, il) = rxyz(3, i) - alat(3)
      END DO

      ALLOCATE (lsta(2, nat))
      nnbrx = 36
      loop_nnbrx: DO
         ALLOCATE (lstb(nnbrx*nat), rel(5, nnbrx*nat))

         indlstx = 0

!$OMP PARALLEL DEFAULT(NONE)  &
!$OMP PRIVATE(iat,cut2,iam,ii,indlst,l1,l2,l3,myspace,npr) &
!$OMP SHARED (indlstx,nat,nn,nnbrx,ncx,ll1,ll2,ll3,icell,lsta,lstb,lay, &
!$OMP rel,rxyz,cut,myspaceout)

         npr = 1
!$       npr = omp_get_num_threads()
         iam = 0
!$       iam = omp_get_thread_num()

         cut2 = cut**2
! assign contiguous portions of the arrays lstb and rel to the threads
         myspace = (nat*nnbrx)/npr
         IF (iam == 0) myspaceout = myspace
! Verlet list, relative positions
         indlst = 0
         loop_l3: DO l3 = 0, ll3 - 1
            loop_l2: DO l2 = 0, ll2 - 1
               loop_l1: DO l1 = 0, ll1 - 1
                  loop_ii: DO ii = 1, icell(0, l1, l2, l3)
                     iat = icell(ii, l1, l2, l3)
                     IF (((iat - 1)*npr)/nat == iam) THEN
!                          write(*,*) 'sublstiat:iam,iat',iam,iat
                        lsta(1, iat) = iam*myspace + indlst + 1
                        CALL sublstiat_l(iat, nn, ncx, ll1, ll2, ll3, l1, l2, l3, myspace, &
                                         rxyz, icell, lstb(iam*myspace + 1), lay, &
                                         rel(1, iam*myspace + 1), cut2, indlst)
                        lsta(2, iat) = iam*myspace + indlst
!                          write(*,'(a,4(x,i3),100(x,i2))') &
!                                'iam,iat,lsta',iam,iat,lsta(1,iat),lsta(2,iat), &
!                                (lstb(j),j=lsta(1,iat),lsta(2,iat))
                     END IF
                  END DO loop_ii
               END DO loop_l1
            END DO loop_l2
         END DO loop_l3

!$OMP ATOMIC UPDATE
         indlstx = MAX(indlstx, indlst)
!$OMP END ATOMIC
!$OMP END PARALLEL

         IF (indlstx >= myspaceout) THEN
            WRITE (10, *) count, 'NNBRX too  small', nnbrx
            DEALLOCATE (lstb, rel)
            nnbrx = 3*nnbrx/2
            CYCLE loop_nnbrx
         END IF
         EXIT loop_nnbrx
      END DO loop_nnbrx

      istopg = 0
!$OMP PARALLEL DEFAULT(NONE)  &
!$OMP PRIVATE(iam,npr,iat,iat1,iat2,lot,istop,tcoord,tcoord2, &
!$OMP tener,tener2,txyz,f2ij,f3ij,f3ik,npjx,npjkx) &
!$OMP SHARED (nat,nnbrx,lsta,lstb,rel,ener,ener2,fxyz,coord,coord2,istopg)

      npr = 1
!$    npr = omp_get_num_threads()
      iam = 0
!$    iam = omp_get_thread_num()

      npjx = 300; npjkx = 6000

      IF (npr /= 1) THEN
! PARALLEL CASE
! create temporary private scalars for reduction sum on energies and
!        temporary private array for reduction sum on forces
!$OMP CRITICAL(omp_eip_lenosky_silicon)
         ALLOCATE (txyz(3, nat), f2ij(3, npjx), f3ij(3, npjkx), f3ik(3, npjkx))
!$OMP END CRITICAL(omp_eip_lenosky_silicon)
         IF (iam == 0) THEN
            ener = 0.e0_dp
            ener2 = 0.e0_dp
            coord = 0.e0_dp
            coord2 = 0.e0_dp
         END IF
!$OMP DO
         DO iat = 1, nat
            fxyz(1, iat) = 0.e0_dp
            fxyz(2, iat) = 0.e0_dp
            fxyz(3, iat) = 0.e0_dp
         END DO
!$OMP BARRIER

! Each thread treats at most lot atoms
         lot = INT(REAL(nat, KIND=dp)/REAL(npr, KIND=dp) + .999999999999e0_dp)
         iat1 = iam*lot + 1
         iat2 = MIN((iam + 1)*lot, nat)
!       write(*,*) 'subfeniat:iat1,iat2,iam',iat1,iat2,iam
         CALL subfeniat_l(iat1, iat2, nat, lsta, lstb, rel, tener, tener2, &
                          tcoord, tcoord2, nnbrx, txyz, f2ij, npjx, f3ij, npjkx, f3ik, istop)
!$OMP CRITICAL(omp_eip_lenosky_silicon)
         ener = ener + tener
         ener2 = ener2 + tener2
         coord = coord + tcoord
         coord2 = coord2 + tcoord2
         istopg = istopg + istop
         DO iat = 1, nat
            fxyz(1, iat) = fxyz(1, iat) + txyz(1, iat)
            fxyz(2, iat) = fxyz(2, iat) + txyz(2, iat)
            fxyz(3, iat) = fxyz(3, iat) + txyz(3, iat)
         END DO
         DEALLOCATE (txyz, f2ij, f3ij, f3ik)
!$OMP END CRITICAL(omp_eip_lenosky_silicon)

      ELSE
! SERIAL CASE
         iat1 = 1
         iat2 = nat
         ALLOCATE (f2ij(3, npjx), f3ij(3, npjkx), f3ik(3, npjkx))
         CALL subfeniat_l(iat1, iat2, nat, lsta, lstb, rel, ener, ener2, &
                          coord, coord2, nnbrx, fxyz, f2ij, npjx, f3ij, npjkx, f3ik, istopg)
         DEALLOCATE (f2ij, f3ij, f3ik)

      END IF
!$OMP END PARALLEL

!         write(*,*) 'ener,norm force', &
!                    ener,DNRM2(3*nat,fxyz,1)
      IF (istopg > 0) CPABORT("DIMENSION ERROR (see WARNING above)")
      ener_var = ener2/nat - (ener/nat)**2
      coord = coord/nat
      coord_var = coord2/nat - coord**2

      DEALLOCATE (rxyz, icell, lay, lsta, lstb, rel)

   END SUBROUTINE eip_lenosky_silicon

! **************************************************************************************************
!> \brief ...
!> \param iat1 ...
!> \param iat2 ...
!> \param nat ...
!> \param lsta ...
!> \param lstb ...
!> \param rel ...
!> \param tener ...
!> \param tener2 ...
!> \param tcoord ...
!> \param tcoord2 ...
!> \param nnbrx ...
!> \param txyz ...
!> \param f2ij ...
!> \param npjx ...
!> \param f3ij ...
!> \param npjkx ...
!> \param f3ik ...
!> \param istop ...
! **************************************************************************************************
   SUBROUTINE subfeniat_l(iat1, iat2, nat, lsta, lstb, rel, tener, tener2, &
                          tcoord, tcoord2, nnbrx, txyz, f2ij, npjx, f3ij, npjkx, f3ik, istop)
! for a subset of atoms iat1 to iat2 the routine calculates the (partial) forces
! txyz acting on these atoms as well as on the atoms (jat, kat) interacting
! with them and their contribution to the energy (tener).
! In addition the coordination number tcoord and the second moment of the
! local energy tener2 and coordination number tcoord2 are returned
      INTEGER                                            :: iat1, iat2, nat, lsta, lstb
      REAL(KIND=dp)                                      :: rel, tener, tener2, tcoord, tcoord2
      INTEGER                                            :: nnbrx
      REAL(KIND=dp)                                      :: txyz, f2ij
      INTEGER                                            :: npjx
      REAL(KIND=dp)                                      :: f3ij
      INTEGER                                            :: npjkx
      REAL(KIND=dp)                                      :: f3ik
      INTEGER                                            :: istop

      DIMENSION lsta(2, nat), lstb(nnbrx*nat), rel(5, nnbrx*nat), txyz(3, nat)
      DIMENSION f2ij(3, npjx), f3ij(3, npjkx), f3ik(3, npjkx)
      REAL(KIND=dp), PARAMETER :: tmin_phi = 0.1500000e+01_dp
      REAL(KIND=dp), PARAMETER :: tmax_phi = 0.4500000e+01_dp
      REAL(KIND=dp), PARAMETER :: hi_phi = 3.00000000000e0_dp
      REAL(KIND=dp), PARAMETER :: hsixth_phi = 5.55555555555556e-002_dp
      REAL(KIND=dp), PARAMETER :: h2sixth_phi = 1.85185185185185e-002_dp
      REAL(KIND=dp), PARAMETER, DIMENSION(0:9) :: cof_phi = &
                                                  [0.69299400000000e+01_dp, -0.43995000000000e+00_dp, &
                                                   -0.17012300000000e+01_dp, -0.16247300000000e+01_dp, &
                                                   -0.99696000000000e+00_dp, -0.27391000000000e+00_dp, &
                                                   -0.24990000000000e-01_dp, -0.17840000000000e-01_dp, &
                                                   -0.96100000000000e-02_dp, 0.00000000000000e+00_dp]
      REAL(KIND=dp), PARAMETER, DIMENSION(0:9) :: dof_phi = &
                                                  [0.16533229480429e+03_dp, 0.39415410391417e+02_dp, &
                                                   0.68710036300407e+01_dp, 0.53406950884203e+01_dp, &
                                                   0.15347960162782e+01_dp, -0.63347591535331e+01_dp, &
                                                   -0.17987794021458e+01_dp, 0.47429676211617e+00_dp, &
                                                   -0.40087646318907e-01_dp, -0.23942617684055e+00_dp]
      REAL(KIND=dp), PARAMETER :: tmin_rho = 0.1500000e+01_dp
      REAL(KIND=dp), PARAMETER :: tmax_rho = 0.3500000e+01_dp
      REAL(KIND=dp), PARAMETER :: hi_rho = 5.00000000000e0_dp
      REAL(KIND=dp), PARAMETER :: hsixth_rho = 3.33333333333333e-002_dp
      REAL(KIND=dp), PARAMETER :: h2sixth_rho = 6.66666666666667e-003_dp
      REAL(KIND=dp), PARAMETER, DIMENSION(0:10) :: cof_rho = &
                                                   [0.13747000000000e+00_dp, -0.14831000000000e+00_dp, &
                                                    -0.55972000000000e+00_dp, -0.73110000000000e+00_dp, &
                                                    -0.76283000000000e+00_dp, -0.72918000000000e+00_dp, &
                                                    -0.66620000000000e+00_dp, -0.57328000000000e+00_dp, &
                                                    -0.40690000000000e+00_dp, -0.16662000000000e+00_dp, &
                                                    0.00000000000000e+00_dp]
      REAL(KIND=dp), PARAMETER, DIMENSION(0:10) :: dof_rho = &
                                                   [-0.32275496741918e+01_dp, -0.64119006516165e+01_dp, &
                                                    0.10030652280658e+02_dp, 0.22937915289857e+01_dp, &
                                                    0.17416816033995e+01_dp, 0.54648205741626e+00_dp, &
                                                    0.47189016693543e+00_dp, 0.20569572748420e+01_dp, &
                                                    0.23192807336964e+01_dp, -0.24908020962757e+00_dp, &
                                                    -0.12371959895186e+02_dp]
      REAL(KIND=dp), PARAMETER :: tmin_fff = 0.1500000e+01_dp
      REAL(KIND=dp), PARAMETER :: tmax_fff = 0.3500000e+01_dp
      REAL(KIND=dp), PARAMETER :: hi_fff = 4.50000000000e0_dp
      REAL(KIND=dp), PARAMETER :: hsixth_fff = 3.70370370370370e-002_dp
      REAL(KIND=dp), PARAMETER :: h2sixth_fff = 8.23045267489712e-003_dp
      REAL(KIND=dp), PARAMETER, DIMENSION(0:9) :: cof_fff = &
                                                  [0.12503100000000e+01_dp, 0.86821000000000e+00_dp, &
                                                   0.60846000000000e+00_dp, 0.48756000000000e+00_dp, &
                                                   0.44163000000000e+00_dp, 0.37610000000000e+00_dp, &
                                                   0.27145000000000e+00_dp, 0.14814000000000e+00_dp, &
                                                   0.48550000000000e-01_dp, 0.00000000000000e+00_dp]
      REAL(KIND=dp), PARAMETER, DIMENSION(0:9) :: dof_fff = &
                                                  [0.27904652711432e+02_dp, -0.45230754228635e+01_dp, &
                                                   0.50531739800222e+01_dp, 0.11806545027747e+01_dp, &
                                                   -0.66693699112098e+00_dp, -0.89430653829079e+00_dp, &
                                                   -0.50891685571587e+00_dp, 0.66278396115427e+00_dp, &
                                                   0.73976101109878e+00_dp, 0.25795319944506e+01_dp]
      REAL(KIND=dp), PARAMETER :: tmin_uuu = -0.1770930e+01_dp
      REAL(KIND=dp), PARAMETER :: tmax_uuu = 0.7908520e+01_dp
      REAL(KIND=dp), PARAMETER :: hi_uuu = 0.723181585730594e0_dp
      REAL(KIND=dp), PARAMETER :: hsixth_uuu = 0.230463095238095e0_dp
      REAL(KIND=dp), PARAMETER :: h2sixth_uuu = 0.318679429600340e0_dp
      REAL(KIND=dp), PARAMETER, DIMENSION(0:7) :: cof_uuu = &
                                                  [-0.10749300000000e+01_dp, -0.20045000000000e+00_dp, &
                                                   0.41422000000000e+00_dp, 0.87939000000000e+00_dp, &
                                                   0.12668900000000e+01_dp, 0.16299800000000e+01_dp, &
                                                   0.19773800000000e+01_dp, 0.23961800000000e+01_dp]
      REAL(KIND=dp), PARAMETER, DIMENSION(0:7) :: dof_uuu = &
                                                  [-0.14827125747284e+00_dp, -0.14922155328475e+00_dp, &
                                                   -0.70113224223509e-01_dp, -0.39449020349230e-01_dp, &
                                                   -0.15815242579643e-01_dp, 0.26112640061855e-01_dp, &
                                                   -0.13786974745095e+00_dp, 0.74941595372657e+00_dp]
      REAL(KIND=dp), PARAMETER :: tmin_ggg = -0.1000000e+01_dp
      REAL(KIND=dp), PARAMETER :: tmax_ggg = 0.8001400e+00_dp
      REAL(KIND=dp), PARAMETER :: hi_ggg = 3.88858644327663e0_dp
      REAL(KIND=dp), PARAMETER :: hsixth_ggg = 4.28604761904762e-002_dp
      REAL(KIND=dp), PARAMETER :: h2sixth_ggg = 1.10221225156463e-002_dp
      REAL(KIND=dp), PARAMETER, DIMENSION(0:7) :: cof_ggg = &
                                                  [0.52541600000000e+01_dp, 0.23591500000000e+01_dp, &
                                                   0.11959500000000e+01_dp, 0.12299500000000e+01_dp, &
                                                   0.20356500000000e+01_dp, 0.34247400000000e+01_dp, &
                                                   0.49485900000000e+01_dp, 0.56179900000000e+01_dp]
      REAL(KIND=dp), PARAMETER, DIMENSION(0:7) :: dof_ggg = &
                                                  [0.15826876132396e+02_dp, 0.31176239377907e+02_dp, &
                                                   0.16589446539683e+02_dp, 0.11083892500520e+02_dp, &
                                                   0.90887216383860e+01_dp, 0.54902279653967e+01_dp, &
                                                   -0.18823313223755e+02_dp, -0.77183416481005e+01_dp]

      REAL(KIND=dp) :: a2_fff, a2_ggg, a_fff, a_ggg, b2_fff, b2_ggg, b_fff, &
                       b_ggg, cof1_fff, cof1_ggg, cof2_fff, cof2_ggg, cof3_fff, &
                       cof3_ggg, cof4_fff, cof4_ggg, cof_fff_khi, cof_fff_klo, &
                       cof_ggg_khi, cof_ggg_klo, coord_iat, costheta, dens, &
                       dens2, dens3, dof_fff_khi, dof_fff_klo, dof_ggg_khi, &
                       dof_ggg_klo, e_phi, e_uuu, ener_iat, ep_phi, ep_uuu, &
                       fij, fijp, fik, fikp, fxij, fxik, fyij, fyik, fzij, fzik, &
                       gjik, gjikp, rho, rhop, rij, rik, sij, sik, t1, t2, t3, t4, &
                       tt, tt_fff, tt_ggg, xarg, ypt1_fff, ypt1_ggg, ypt2_fff, &
                       ypt2_ggg, yt1_fff, yt1_ggg, yt2_fff, yt2_ggg

      INTEGER       :: iat, jat, jbr, jcnt, jkcnt, kat, kbr, khi_fff, khi_ggg, &
                       klo_fff, klo_ggg

! initialize temporary private scalars for reduction sum on energies and
! private workarray txyz for forces forces
      tener = 0.e0_dp
      tener2 = 0.e0_dp
      tcoord = 0.e0_dp
      tcoord2 = 0.e0_dp
      istop = 0
      DO iat = 1, nat
         txyz(1, iat) = 0.e0_dp
         txyz(2, iat) = 0.e0_dp
         txyz(3, iat) = 0.e0_dp
      END DO

! calculation of forces, energy

      forces_and_energy: DO iat = iat1, iat2

         dens2 = 0.e0_dp
         dens3 = 0.e0_dp
         jcnt = 0
         jkcnt = 0
         coord_iat = 0.e0_dp
         ener_iat = 0.e0_dp
         calculate: DO jbr = lsta(1, iat), lsta(2, iat)
            jat = lstb(jbr)
            jcnt = jcnt + 1
            IF (jcnt > npjx) THEN
               WRITE (*, *) 'WARNING: enlarge npjx'
               istop = 1
               RETURN
            END IF

            fxij = rel(1, jbr)
            fyij = rel(2, jbr)
            fzij = rel(3, jbr)
            rij = rel(4, jbr)
            sij = rel(5, jbr)

! coordination number calculated with soft cutoff between first
! nearest neighbor and midpoint of first and second nearest neighbor
            IF (rij <= 2.36e0_dp) THEN
               coord_iat = coord_iat + 1.e0_dp
            ELSE IF (rij >= 3.12e0_dp) THEN
            ELSE
               xarg = (rij - 2.36e0_dp)*(1.e0_dp/(3.12e0_dp - 2.36e0_dp))
               coord_iat = coord_iat + (2*xarg + 1.e0_dp)*(xarg - 1.e0_dp)**2
            END IF

! pairpotential term
            CALL splint(cof_phi, dof_phi, tmin_phi, tmax_phi, &
                        hsixth_phi, h2sixth_phi, hi_phi, 10, rij, e_phi, ep_phi)
            ener_iat = ener_iat + (e_phi*.5e0_dp)
            txyz(1, iat) = txyz(1, iat) - fxij*(ep_phi*.5e0_dp)
            txyz(2, iat) = txyz(2, iat) - fyij*(ep_phi*.5e0_dp)
            txyz(3, iat) = txyz(3, iat) - fzij*(ep_phi*.5e0_dp)
            txyz(1, jat) = txyz(1, jat) + fxij*(ep_phi*.5e0_dp)
            txyz(2, jat) = txyz(2, jat) + fyij*(ep_phi*.5e0_dp)
            txyz(3, jat) = txyz(3, jat) + fzij*(ep_phi*.5e0_dp)

! 2 body embedding term
            CALL splint(cof_rho, dof_rho, tmin_rho, tmax_rho, &
                        hsixth_rho, h2sixth_rho, hi_rho, 11, rij, rho, rhop)
            dens2 = dens2 + rho
            f2ij(1, jcnt) = fxij*rhop
            f2ij(2, jcnt) = fyij*rhop
            f2ij(3, jcnt) = fzij*rhop

! 3 body embedding term
            CALL splint(cof_fff, dof_fff, tmin_fff, tmax_fff, &
                        hsixth_fff, h2sixth_fff, hi_fff, 10, rij, fij, fijp)

            embed_3body: DO kbr = lsta(1, iat), lsta(2, iat)
               kat = lstb(kbr)
               IF (kat < jat) THEN
                  jkcnt = jkcnt + 1
                  IF (jkcnt > npjkx) THEN
                     WRITE (*, *) 'WARNING: enlarge npjkx', npjkx
                     istop = 1
                     RETURN
                  END IF

! begin unoptimized original version:
!        fxik=rel(1,kbr)
!        fyik=rel(2,kbr)
!        fzik=rel(3,kbr)
!        rik=rel(4,kbr)
!        sik=rel(5,kbr)
!
!        call splint(cof_fff,dof_fff,tmin_fff,tmax_fff, &
!             hsixth_fff,h2sixth_fff,hi_fff,10,rik,fik,fikp)
!        costheta=fxij*fxik+fyij*fyik+fzij*fzik
!        call splint(cof_ggg,dof_ggg,tmin_ggg,tmax_ggg, &
!             hsixth_ggg,h2sixth_ggg,hi_ggg,8,costheta,gjik,gjikp)
! end unoptimized original version:

! begin optimized version
                  rik = rel(4, kbr)
                  IF (rik > tmax_fff) THEN
                     fikp = 0.e0_dp; fik = 0.e0_dp
                     gjik = 0.e0_dp; gjikp = 0.e0_dp; sik = 0.e0_dp
                     costheta = 0.e0_dp; fxik = 0.e0_dp; fyik = 0.e0_dp; fzik = 0.e0_dp
                  ELSE IF (rik < tmin_fff) THEN
                     fxik = rel(1, kbr)
                     fyik = rel(2, kbr)
                     fzik = rel(3, kbr)
                     costheta = fxij*fxik + fyij*fyik + fzij*fzik
                     sik = rel(5, kbr)
                     fikp = hi_fff*(cof_fff(1) - cof_fff(0)) - &
                            (dof_fff(1) + 2.e0_dp*dof_fff(0))*hsixth_fff
                     fik = cof_fff(0) + (rik - tmin_fff)*fikp
                     tt_ggg = (costheta - tmin_ggg)*hi_ggg
                     IF (costheta > tmax_ggg) THEN
                        gjikp = hi_ggg*(cof_ggg(8 - 1) - cof_ggg(8 - 2)) + &
                                (2.e0_dp*dof_ggg(8 - 1) + dof_ggg(8 - 2))*hsixth_ggg
                        gjik = cof_ggg(8 - 1) + (costheta - tmax_ggg)*gjikp
                     ELSE
                        klo_ggg = INT(tt_ggg)
                        khi_ggg = klo_ggg + 1
                        cof_ggg_klo = cof_ggg(klo_ggg)
                        dof_ggg_klo = dof_ggg(klo_ggg)
                        b_ggg = tt_ggg - klo_ggg
                        a_ggg = 1.e0_dp - b_ggg
                        cof_ggg_khi = cof_ggg(khi_ggg)
                        dof_ggg_khi = dof_ggg(khi_ggg)
                        b2_ggg = b_ggg*b_ggg
                        gjik = a_ggg*cof_ggg_klo
                        gjikp = cof_ggg_khi - cof_ggg_klo
                        a2_ggg = a_ggg*a_ggg
                        cof1_ggg = a2_ggg - 1.e0_dp
                        cof2_ggg = b2_ggg - 1.e0_dp
                        gjik = gjik + b_ggg*cof_ggg_khi
                        gjikp = hi_ggg*gjikp
                        cof3_ggg = 3.e0_dp*b2_ggg
                        cof4_ggg = 3.e0_dp*a2_ggg
                        cof1_ggg = a_ggg*cof1_ggg
                        cof2_ggg = b_ggg*cof2_ggg
                        cof3_ggg = cof3_ggg - 1.e0_dp
                        cof4_ggg = cof4_ggg - 1.e0_dp
                        yt1_ggg = cof1_ggg*dof_ggg_klo
                        yt2_ggg = cof2_ggg*dof_ggg_khi
                        ypt1_ggg = cof3_ggg*dof_ggg_khi
                        ypt2_ggg = cof4_ggg*dof_ggg_klo
                        gjik = gjik + (yt1_ggg + yt2_ggg)*h2sixth_ggg
                        gjikp = gjikp + (ypt1_ggg - ypt2_ggg)*hsixth_ggg
                     END IF
                  ELSE
                     fxik = rel(1, kbr)
                     tt_fff = rik - tmin_fff
                     costheta = fxij*fxik
                     fyik = rel(2, kbr)
                     tt_fff = tt_fff*hi_fff
                     costheta = costheta + fyij*fyik
                     fzik = rel(3, kbr)
                     klo_fff = INT(tt_fff)
                     costheta = costheta + fzij*fzik
                     sik = rel(5, kbr)
                     tt_ggg = (costheta - tmin_ggg)*hi_ggg
                     IF (costheta > tmax_ggg) THEN
                        gjikp = hi_ggg*(cof_ggg(8 - 1) - cof_ggg(8 - 2)) + &
                                (2.e0_dp*dof_ggg(8 - 1) + dof_ggg(8 - 2))*hsixth_ggg
                        gjik = cof_ggg(8 - 1) + (costheta - tmax_ggg)*gjikp
                        khi_fff = klo_fff + 1
                        cof_fff_klo = cof_fff(klo_fff)
                        dof_fff_klo = dof_fff(klo_fff)
                        b_fff = tt_fff - klo_fff
                        a_fff = 1.e0_dp - b_fff
                        cof_fff_khi = cof_fff(khi_fff)
                        dof_fff_khi = dof_fff(khi_fff)
                        b2_fff = b_fff*b_fff
                        fik = a_fff*cof_fff_klo
                        fikp = cof_fff_khi - cof_fff_klo
                        a2_fff = a_fff*a_fff
                        cof1_fff = a2_fff - 1.e0_dp
                        cof2_fff = b2_fff - 1.e0_dp
                        fik = fik + b_fff*cof_fff_khi
                        fikp = hi_fff*fikp
                        cof3_fff = 3.e0_dp*b2_fff
                        cof4_fff = 3.e0_dp*a2_fff
                        cof1_fff = a_fff*cof1_fff
                        cof2_fff = b_fff*cof2_fff
                        cof3_fff = cof3_fff - 1.e0_dp
                        cof4_fff = cof4_fff - 1.e0_dp
                        yt1_fff = cof1_fff*dof_fff_klo
                        yt2_fff = cof2_fff*dof_fff_khi
                        ypt1_fff = cof3_fff*dof_fff_khi
                        ypt2_fff = cof4_fff*dof_fff_klo
                        fik = fik + (yt1_fff + yt2_fff)*h2sixth_fff
                        fikp = fikp + (ypt1_fff - ypt2_fff)*hsixth_fff
                     ELSE
                        klo_ggg = INT(tt_ggg)
                        khi_ggg = klo_ggg + 1
                        khi_fff = klo_fff + 1
                        cof_ggg_klo = cof_ggg(klo_ggg)
                        cof_fff_klo = cof_fff(klo_fff)
                        dof_ggg_klo = dof_ggg(klo_ggg)
                        dof_fff_klo = dof_fff(klo_fff)
                        b_ggg = tt_ggg - klo_ggg
                        b_fff = tt_fff - klo_fff
                        a_ggg = 1.e0_dp - b_ggg
                        a_fff = 1.e0_dp - b_fff
                        cof_ggg_khi = cof_ggg(khi_ggg)
                        cof_fff_khi = cof_fff(khi_fff)
                        dof_ggg_khi = dof_ggg(khi_ggg)
                        dof_fff_khi = dof_fff(khi_fff)
                        b2_ggg = b_ggg*b_ggg
                        b2_fff = b_fff*b_fff
                        gjik = a_ggg*cof_ggg_klo
                        fik = a_fff*cof_fff_klo
                        gjikp = cof_ggg_khi - cof_ggg_klo
                        fikp = cof_fff_khi - cof_fff_klo
                        a2_ggg = a_ggg*a_ggg
                        a2_fff = a_fff*a_fff
                        cof1_ggg = a2_ggg - 1.e0_dp
                        cof1_fff = a2_fff - 1.e0_dp
                        cof2_ggg = b2_ggg - 1.e0_dp
                        cof2_fff = b2_fff - 1.e0_dp
                        gjik = gjik + b_ggg*cof_ggg_khi
                        fik = fik + b_fff*cof_fff_khi
                        gjikp = hi_ggg*gjikp
                        fikp = hi_fff*fikp
                        cof3_ggg = 3.e0_dp*b2_ggg
                        cof3_fff = 3.e0_dp*b2_fff
                        cof4_ggg = 3.e0_dp*a2_ggg
                        cof4_fff = 3.e0_dp*a2_fff
                        cof1_ggg = a_ggg*cof1_ggg
                        cof1_fff = a_fff*cof1_fff
                        cof2_ggg = b_ggg*cof2_ggg
                        cof2_fff = b_fff*cof2_fff
                        cof3_ggg = cof3_ggg - 1.e0_dp
                        cof3_fff = cof3_fff - 1.e0_dp
                        cof4_ggg = cof4_ggg - 1.e0_dp
                        cof4_fff = cof4_fff - 1.e0_dp
                        yt1_ggg = cof1_ggg*dof_ggg_klo
                        yt1_fff = cof1_fff*dof_fff_klo
                        yt2_ggg = cof2_ggg*dof_ggg_khi
                        yt2_fff = cof2_fff*dof_fff_khi
                        ypt1_ggg = cof3_ggg*dof_ggg_khi
                        ypt1_fff = cof3_fff*dof_fff_khi
                        ypt2_ggg = cof4_ggg*dof_ggg_klo
                        ypt2_fff = cof4_fff*dof_fff_klo
                        gjik = gjik + (yt1_ggg + yt2_ggg)*h2sixth_ggg
                        fik = fik + (yt1_fff + yt2_fff)*h2sixth_fff
                        gjikp = gjikp + (ypt1_ggg - ypt2_ggg)*hsixth_ggg
                        fikp = fikp + (ypt1_fff - ypt2_fff)*hsixth_fff
                     END IF
                  END IF
! end optimized version

                  tt = fij*fik
                  dens3 = dens3 + tt*gjik

                  t1 = fijp*fik*gjik
                  t2 = sij*(tt*gjikp)
                  f3ij(1, jkcnt) = fxij*t1 + (fxik - fxij*costheta)*t2
                  f3ij(2, jkcnt) = fyij*t1 + (fyik - fyij*costheta)*t2
                  f3ij(3, jkcnt) = fzij*t1 + (fzik - fzij*costheta)*t2

                  t3 = fikp*fij*gjik
                  t4 = sik*(tt*gjikp)
                  f3ik(1, jkcnt) = fxik*t3 + (fxij - fxik*costheta)*t4
                  f3ik(2, jkcnt) = fyik*t3 + (fyij - fyik*costheta)*t4
                  f3ik(3, jkcnt) = fzik*t3 + (fzij - fzik*costheta)*t4
               END IF

            END DO embed_3body
         END DO calculate

         dens = dens2 + dens3
         CALL splint(cof_uuu, dof_uuu, tmin_uuu, tmax_uuu, &
                     hsixth_uuu, h2sixth_uuu, hi_uuu, 8, dens, e_uuu, ep_uuu)
         ener_iat = ener_iat + e_uuu

! Only now ep_uu is known and the forces can be calculated, lets loop again
         jcnt = 0
         jkcnt = 0
         loop_again: DO jbr = lsta(1, iat), lsta(2, iat)
            jat = lstb(jbr)
            jcnt = jcnt + 1
            txyz(1, iat) = txyz(1, iat) - ep_uuu*f2ij(1, jcnt)
            txyz(2, iat) = txyz(2, iat) - ep_uuu*f2ij(2, jcnt)
            txyz(3, iat) = txyz(3, iat) - ep_uuu*f2ij(3, jcnt)
            txyz(1, jat) = txyz(1, jat) + ep_uuu*f2ij(1, jcnt)
            txyz(2, jat) = txyz(2, jat) + ep_uuu*f2ij(2, jcnt)
            txyz(3, jat) = txyz(3, jat) + ep_uuu*f2ij(3, jcnt)

! 3 body embedding term
            DO kbr = lsta(1, iat), lsta(2, iat)
               kat = lstb(kbr)
               IF (kat < jat) THEN
                  jkcnt = jkcnt + 1

                  txyz(1, iat) = txyz(1, iat) - ep_uuu*(f3ij(1, jkcnt) + f3ik(1, jkcnt))
                  txyz(2, iat) = txyz(2, iat) - ep_uuu*(f3ij(2, jkcnt) + f3ik(2, jkcnt))
                  txyz(3, iat) = txyz(3, iat) - ep_uuu*(f3ij(3, jkcnt) + f3ik(3, jkcnt))
                  txyz(1, jat) = txyz(1, jat) + ep_uuu*f3ij(1, jkcnt)
                  txyz(2, jat) = txyz(2, jat) + ep_uuu*f3ij(2, jkcnt)
                  txyz(3, jat) = txyz(3, jat) + ep_uuu*f3ij(3, jkcnt)
                  txyz(1, kat) = txyz(1, kat) + ep_uuu*f3ik(1, jkcnt)
                  txyz(2, kat) = txyz(2, kat) + ep_uuu*f3ik(2, jkcnt)
                  txyz(3, kat) = txyz(3, kat) + ep_uuu*f3ik(3, jkcnt)
               END IF
            END DO

         END DO loop_again

!        write(*,'(a,i4,x,e19.12,x,e10.3)') 'iat,ener_iat,coord_iat', &
!                                       iat,ener_iat,coord_iat
         tener = tener + ener_iat
         tener2 = tener2 + ener_iat**2
         tcoord = tcoord + coord_iat
         tcoord2 = tcoord2 + coord_iat**2

      END DO forces_and_energy

   END SUBROUTINE subfeniat_l

! **************************************************************************************************
!> \brief ...
!> \param iat ...
!> \param nn ...
!> \param ncx ...
!> \param ll1 ...
!> \param ll2 ...
!> \param ll3 ...
!> \param l1 ...
!> \param l2 ...
!> \param l3 ...
!> \param myspace ...
!> \param rxyz ...
!> \param icell ...
!> \param lstb ...
!> \param lay ...
!> \param rel ...
!> \param cut2 ...
!> \param indlst ...
! **************************************************************************************************
   SUBROUTINE sublstiat_l(iat, nn, ncx, ll1, ll2, ll3, l1, l2, l3, myspace, &
                          rxyz, icell, lstb, lay, rel, cut2, indlst)
! finds the neighbours of atom iat (specified by lsta and lstb) and and
! the relative position rel of iat with respect to these neighbours
      INTEGER                                            :: iat, nn, ncx, ll1, ll2, ll3, l1, l2, l3, &
                                                            myspace
      REAL(KIND=dp)                                      :: rxyz
      INTEGER                                            :: icell, lstb, lay
      REAL(KIND=dp)                                      :: rel, cut2
      INTEGER                                            :: indlst

      DIMENSION rxyz(3, nn), lay(nn), icell(0:ncx, -1:ll1, -1:ll2, -1:ll3), &
         lstb(0:myspace - 1), rel(5, 0:myspace - 1)

      INTEGER       :: jat, jj, k1, k2, k3
      REAL(KIND=dp) :: rr2, tt, xrel, yrel, zrel, tti

      loop_k3: DO k3 = l3 - 1, l3 + 1
         loop_k2: DO k2 = l2 - 1, l2 + 1
            loop_k1: DO k1 = l1 - 1, l1 + 1
               loop_jj: DO jj = 1, icell(0, k1, k2, k3)
                  jat = icell(jj, k1, k2, k3)
                  IF (jat == iat) CYCLE loop_k3
                  xrel = rxyz(1, iat) - rxyz(1, jat)
                  yrel = rxyz(2, iat) - rxyz(2, jat)
                  zrel = rxyz(3, iat) - rxyz(3, jat)
                  rr2 = xrel**2 + yrel**2 + zrel**2
                  IF (rr2 <= cut2) THEN
                     indlst = MIN(indlst, myspace - 1)
                     lstb(indlst) = lay(jat)
!                       write(*,*) 'iat,indlst,lay(jat)',iat,indlst,lay(jat)
                     tt = SQRT(rr2)
                     tti = 1.e0_dp/tt
                     rel(1, indlst) = xrel*tti
                     rel(2, indlst) = yrel*tti
                     rel(3, indlst) = zrel*tti
                     rel(4, indlst) = tt
                     rel(5, indlst) = tti
                     indlst = indlst + 1
                  END IF
               END DO loop_jj
            END DO loop_k1
         END DO loop_k2
      END DO loop_k3

      RETURN
   END SUBROUTINE sublstiat_l

! **************************************************************************************************
!> \brief ...
!> \param ya ...
!> \param y2a ...
!> \param tmin ...
!> \param tmax ...
!> \param hsixth ...
!> \param h2sixth ...
!> \param hi ...
!> \param n ...
!> \param x ...
!> \param y ...
!> \param yp ...
! **************************************************************************************************
   SUBROUTINE splint(ya, y2a, tmin, tmax, hsixth, h2sixth, hi, n, x, y, yp)
      REAL(KIND=dp)                                      :: ya, y2a, tmin, tmax, hsixth, h2sixth, hi
      INTEGER                                            :: n
      REAL(KIND=dp)                                      :: x, y, yp

      DIMENSION y2a(0:n - 1), ya(0:n - 1)
      REAL(KIND=dp) :: a, a2, b, b2, cof1, cof2, cof3, cof4, tt, &
                       y2a_khi, ya_klo, y2a_klo, ya_khi, ypt1, ypt2, yt1, yt2
      INTEGER :: klo, khi

! interpolate if the argument is outside the cubic spline interval [tmin,tmax]
      tt = (x - tmin)*hi
      IF (x < tmin) THEN
         yp = hi*(ya(1) - ya(0)) - &
              (y2a(1) + 2.e0_dp*y2a(0))*hsixth
         y = ya(0) + (x - tmin)*yp
      ELSE IF (x > tmax) THEN
         yp = hi*(ya(n - 1) - ya(n - 2)) + &
              (2.e0_dp*y2a(n - 1) + y2a(n - 2))*hsixth
         y = ya(n - 1) + (x - tmax)*yp
! otherwise evaluate cubic spline
      ELSE
         klo = INT(tt)
         khi = klo + 1
         ya_klo = ya(klo)
         y2a_klo = y2a(klo)
         b = tt - klo
         a = 1.e0_dp - b
         ya_khi = ya(khi)
         y2a_khi = y2a(khi)
         b2 = b*b
         y = a*ya_klo
         yp = ya_khi - ya_klo
         a2 = a*a
         cof1 = a2 - 1.e0_dp
         cof2 = b2 - 1.e0_dp
         y = y + b*ya_khi
         yp = hi*yp
         cof3 = 3.e0_dp*b2
         cof4 = 3.e0_dp*a2
         cof1 = a*cof1
         cof2 = b*cof2
         cof3 = cof3 - 1.e0_dp
         cof4 = cof4 - 1.e0_dp
         yt1 = cof1*y2a_klo
         yt2 = cof2*y2a_khi
         ypt1 = cof3*y2a_khi
         ypt2 = cof4*y2a_klo
         y = y + (yt1 + yt2)*h2sixth
         yp = yp + (ypt1 - ypt2)*hsixth
      END IF
      RETURN
   END SUBROUTINE splint



! **************************************************************************************************
!> Additional EIP kernels consolidated here to match the original Bazant/Lenosky layout.
! **************************************************************************************************

subroutine eip_stillinger_weber_silicon(nat,alat, rxyz0, fxyz, etot, count)
!*****************************************************************************************
! This subroutine evaluates the Stillinger Weber Silicon potential with linear scaling
! COPYRIGHT
!    Copyright (C) 2009 AIST, UNIBAS
!    This file is distributed under the terms of the
!    GNU General Public License, see 
!    http://www.gnu.org/copyleft/gpl.txt .
!
! Implementation: Original version was written by Tetsuya Morishita, AIST Tsukuba (JP) 
!                 Improved by M. Amsler, S. Goedecker, Basel University (CH), 2009
!
! Note:
!
!     aa is the parameter A  given on page 5263 of PRB 31, 5262 (1985).
!     bb is the parameter B  given on page 5263 of PRB 31, 5262 (1985).
!     ra is the parameter a  given on page 5263 of PRB 31, 5262 (1985).
!     gam and ramda are the parameters gamma and lambda 
!     for the 3-body term, respectively (see Eq. (2.5) in the paper).
!    
! Input:
!     nat, integer: the number of atoms
!     alat,  real(8), dim(3)   : the three edges of the orthoromic simulation cell, periodic boundaries are applied
!                                and atoms outside the cell will be brought back into the box
!     rxyz, real(8), dim(3,nat) : the xyz cartesian components of the atomic positions in Angstroem    
!                                
! Output:
!     fxyz,  real(8), dim(3,nat): the xyz cartesian forces in on the corresponding atomic components n eV/A  
!     etot, real(8)            : total potential energy, 2-body and 3-body, in eV
!     count, real(8)           : increased by 1.d0 at each call of this subroutine,
!                                needs to be initialized to 0.d0 before calling this routine for the first time
!
! Other variables:
!     p:  the 2-body potential energy
!     p3: the 3-body potential energy 
!     fx(i) , fy(i) , fz(i)  are the 2-body forces on atom i.
!     fx3(i), fy3(i), fz3(i) are the 3-body forces on atom i.
!     fxyz(3,nat) contains both 2-body and 3-body forces  
!
!     All units follow the description in PRB 31, 5262 (1985).
!*****************************************************************************************

        implicit real(8) (a-h,o-z)
        integer                                  :: nat, nnbrx
        integer                                  :: i, ii, j, nij, k
        integer                                  :: iam, iat, il, in, indlst, indlstx
        integer                                  :: ipb, l1, l2, l3, laymx
        integer                                  :: ll1, ll2, ll3, myspace, myspaceout
        integer                                  :: ncx, ndat, nn, npjkx, npjx, npr
        real(8)               :: rxyz0(3,nat)
        real(8)               :: fxyz(3,nat), etot
        real(8)               :: x(nat), y(nat), z(nat), fx(nat), fy(nat), fz(nat), fx3(nat), fy3(nat), fz3(nat)
        real(8)               :: lx, ly, lz, gam, ramda, aa, ra, raa, bb, xij, yij, zij, vl, f, pe, p, p3
        real(8)               :: xjj, yjj, zjj, rij, xik, yik, zik, xkk, ykk, zkk, force
        real(8)               :: rik, xjk, yjk, zjk, rjk, cosjik, cosijk, cosikj, cosjik3, cosijk3, cosikj3
        real(8)               :: rija, rika, rjka, refi, refj, refk
        real(8)               :: hi, hj, hk, hixij0, hixij1, hixij, hixik0, hixik, HKXIK1
        real(8)               :: hjxij0, hjxij1, hjxij, hkxik0, hixik1, hkxik
        real(8)               :: hjxjk0, hjxjk1, hjxjk, hkxkj0, hkxkj1, hkxkj
        real(8)               :: hiyij0, hiyij1
        real(8)               :: hiyij, hiyik0, hiyik1, hiyik, hjyij0, hjyij1, hjyij
        real(8)               :: hkyik0, hkyik1, hkyik, hjyjk0, hjyjk1, hjyjk, hkykj0, hkykj1, hkykj
        real(8)               :: HIZIJ0, HIZIJ1, HIZIJ, HIZIK0, HIZIK1, HIZIK, HJZIJ0, HJZIJ1, HJZIJ, HKZIK0, HKZIK1, HKZIK
        real(8)               :: HJZJK0, HJZJK1, HJZJK, HKZKJ0, HKZKJ1, HKZKJ, PV3
        real(8)               :: c, d, uc
        real(8)               :: xbar, ybar, zbar, xmax, ymax, zmax, xmin, ymin, zmin, xx, yy, zz
        real(8)               :: alat(3), count, sigma, eps, esigma, isigma, c4, crainv
        real(8)               :: invrij, invrik, invrjk, invrija, invrika, invrjka
        real(8), ALLOCATABLE, DIMENSION(:,:)     :: rxyz
        real(8), ALLOCATABLE, DIMENSION(:,:)     :: rel
        real(8), ALLOCATABLE, DIMENSION(:,:)     :: txyz
        real(8), ALLOCATABLE, DIMENSION(:,:)     :: f2ij,f3ij,f3ik
        integer, ALLOCATABLE, DIMENSION(:,:)     :: lsta
        integer, ALLOCATABLE, DIMENSION(:)       :: lstb
        integer, ALLOCATABLE, DIMENSION(:)       :: lay
        integer, ALLOCATABLE, DIMENSION(:,:,:,:) :: icell
        parameter(gam=1.2d0,ramda=21.0d0,aa=7.049556277d0)
        parameter(ra=1.8d0,bb=0.6022245584d0)
        parameter(uc=4.5d0)
        parameter(sigma=2.0951d0,eps=2.167239428587d0)

        count=count+1.d0
        cut=sigma*ra*2.d0
        isigma=1.d0/sigma
        esigma=eps*isigma

! linear scaling calculation of verlet list, only serial 
        ll1=int(alat(1)/cut)
        if (ll1.lt.1) stop 'alat(1) too small'
        ll2=int(alat(2)/cut)
        if (ll2.lt.1) stop 'alat(2) too small'
        ll3=int(alat(3)/cut)
        if (ll3.lt.1) stop 'alat(3) too small'

        npr=1
        ncx=29
1234    ncx=ncx*2
        allocate(icell(0:ncx,-1:ll1,-1:ll2,-1:ll3))
        do 984,l3=-1,ll3
        do 984,l2=-1,ll2
        do 984,l1=-1,ll1
984     icell(0,l1,l2,l3)=0
        rlc1i=ll1/alat(1)
        rlc2i=ll2/alat(2)
        rlc3i=ll3/alat(3)

        do 983,iat=1,nat
        rxyz0(1,iat)=modulo(modulo(rxyz0(1,iat),alat(1)),alat(1))
        rxyz0(2,iat)=modulo(modulo(rxyz0(2,iat),alat(2)),alat(2))
        rxyz0(3,iat)=modulo(modulo(rxyz0(3,iat),alat(3)),alat(3))
        l1=int(rxyz0(1,iat)*rlc1i)
        l2=int(rxyz0(2,iat)*rlc2i)
        l3=int(rxyz0(3,iat)*rlc3i)

        ii=icell(0,l1,l2,l3)
        ii=ii+1
        icell(0,l1,l2,l3)=ii
        if (ii.gt.ncx) then
        write(10,*) count,'NCX too small',ncx
        deallocate(icell)
        goto 1234
        endif
        icell(ii,l1,l2,l3)=iat
983     continue


! duplicate all atoms within boundary layer
        laymx=ncx*(2*ll1*ll2+2*ll1*ll3+2*ll2*ll3+4*ll1+4*ll2+4*ll3+8)
        nn=nat+laymx
        allocate(rxyz(3,nn),lay(nn))
        do  iat=1,nat
        lay(iat)=iat
        rxyz(1,iat)=rxyz0(1,iat)
        rxyz(2,iat)=rxyz0(2,iat)
        rxyz(3,iat)=rxyz0(3,iat)
        enddo
        il=nat
! xy plane
        do l2=0,ll2-1
        do l1=0,ll1-1

        in=icell(0,l1,l2,0)
        icell(0,l1,l2,ll3)=in
        do ii=1,in
        i=icell(ii,l1,l2,0)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,l1,l2,ll3)=il
        rxyz(1,il)=rxyz(1,i)
        rxyz(2,il)=rxyz(2,i)
        rxyz(3,il)=rxyz(3,i)+alat(3)
        enddo

        in=icell(0,l1,l2,ll3-1)
        icell(0,l1,l2,-1)=in
        do ii=1,in
        i=icell(ii,l1,l2,ll3-1)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,l1,l2,-1)=il
        rxyz(1,il)=rxyz(1,i)
        rxyz(2,il)=rxyz(2,i)
        rxyz(3,il)=rxyz(3,i)-alat(3)
        enddo

        enddo
        enddo


! yz plane
        do l3=0,ll3-1
        do l2=0,ll2-1

        in=icell(0,0,l2,l3)
        icell(0,ll1,l2,l3)=in
        do ii=1,in
        i=icell(ii,0,l2,l3)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,ll1,l2,l3)=il
        rxyz(1,il)=rxyz(1,i)+alat(1)
        rxyz(2,il)=rxyz(2,i)
        rxyz(3,il)=rxyz(3,i)
        enddo

        in=icell(0,ll1-1,l2,l3)
        icell(0,-1,l2,l3)=in
        do ii=1,in
        i=icell(ii,ll1-1,l2,l3)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,-1,l2,l3)=il
        rxyz(1,il)=rxyz(1,i)-alat(1)
        rxyz(2,il)=rxyz(2,i)
        rxyz(3,il)=rxyz(3,i)
        enddo

        enddo
        enddo


! xz plane
        do l3=0,ll3-1
        do l1=0,ll1-1

        in=icell(0,l1,0,l3)
        icell(0,l1,ll2,l3)=in
        do ii=1,in
        i=icell(ii,l1,0,l3)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,l1,ll2,l3)=il
        rxyz(1,il)=rxyz(1,i)
        rxyz(2,il)=rxyz(2,i)+alat(2)
        rxyz(3,il)=rxyz(3,i)
        enddo

        in=icell(0,l1,ll2-1,l3)
        icell(0,l1,-1,l3)=in
        do ii=1,in
        i=icell(ii,l1,ll2-1,l3)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,l1,-1,l3)=il
        rxyz(1,il)=rxyz(1,i)
        rxyz(2,il)=rxyz(2,i)-alat(2)
        rxyz(3,il)=rxyz(3,i)
        enddo

        enddo
        enddo


! x axis
        do l1=0,ll1-1

        in=icell(0,l1,0,0)
        icell(0,l1,ll2,ll3)=in
        do ii=1,in
        i=icell(ii,l1,0,0)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,l1,ll2,ll3)=il
        rxyz(1,il)=rxyz(1,i)
        rxyz(2,il)=rxyz(2,i)+alat(2)
        rxyz(3,il)=rxyz(3,i)+alat(3)
        enddo

        in=icell(0,l1,0,ll3-1)
        icell(0,l1,ll2,-1)=in
        do ii=1,in
        i=icell(ii,l1,0,ll3-1)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,l1,ll2,-1)=il
        rxyz(1,il)=rxyz(1,i)
        rxyz(2,il)=rxyz(2,i)+alat(2)
        rxyz(3,il)=rxyz(3,i)-alat(3)
        enddo

        in=icell(0,l1,ll2-1,0)
        icell(0,l1,-1,ll3)=in
        do ii=1,in
        i=icell(ii,l1,ll2-1,0)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,l1,-1,ll3)=il
        rxyz(1,il)=rxyz(1,i)
        rxyz(2,il)=rxyz(2,i)-alat(2)
        rxyz(3,il)=rxyz(3,i)+alat(3)
        enddo

        in=icell(0,l1,ll2-1,ll3-1)
        icell(0,l1,-1,-1)=in
        do ii=1,in
        i=icell(ii,l1,ll2-1,ll3-1)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,l1,-1,-1)=il
        rxyz(1,il)=rxyz(1,i)
        rxyz(2,il)=rxyz(2,i)-alat(2)
        rxyz(3,il)=rxyz(3,i)-alat(3)
        enddo

        enddo


! y axis
        do l2=0,ll2-1

        in=icell(0,0,l2,0)
        icell(0,ll1,l2,ll3)=in
        do ii=1,in
        i=icell(ii,0,l2,0)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,ll1,l2,ll3)=il
        rxyz(1,il)=rxyz(1,i)+alat(1)
        rxyz(2,il)=rxyz(2,i)
        rxyz(3,il)=rxyz(3,i)+alat(3)
        enddo

        in=icell(0,0,l2,ll3-1)
        icell(0,ll1,l2,-1)=in
        do ii=1,in
        i=icell(ii,0,l2,ll3-1)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,ll1,l2,-1)=il
        rxyz(1,il)=rxyz(1,i)+alat(1)
        rxyz(2,il)=rxyz(2,i)
        rxyz(3,il)=rxyz(3,i)-alat(3)
        enddo

        in=icell(0,ll1-1,l2,0)
        icell(0,-1,l2,ll3)=in
        do ii=1,in
        i=icell(ii,ll1-1,l2,0)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,-1,l2,ll3)=il
        rxyz(1,il)=rxyz(1,i)-alat(1)
        rxyz(2,il)=rxyz(2,i)
        rxyz(3,il)=rxyz(3,i)+alat(3)
        enddo

        in=icell(0,ll1-1,l2,ll3-1)
        icell(0,-1,l2,-1)=in
        do ii=1,in
        i=icell(ii,ll1-1,l2,ll3-1)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,-1,l2,-1)=il
        rxyz(1,il)=rxyz(1,i)-alat(1)
        rxyz(2,il)=rxyz(2,i)
        rxyz(3,il)=rxyz(3,i)-alat(3)
        enddo

        enddo


! z axis
        do l3=0,ll3-1

        in=icell(0,0,0,l3)
        icell(0,ll1,ll2,l3)=in
        do ii=1,in
        i=icell(ii,0,0,l3)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,ll1,ll2,l3)=il
        rxyz(1,il)=rxyz(1,i)+alat(1)
        rxyz(2,il)=rxyz(2,i)+alat(2)
        rxyz(3,il)=rxyz(3,i)
        enddo

        in=icell(0,ll1-1,0,l3)
        icell(0,-1,ll2,l3)=in
        do ii=1,in
        i=icell(ii,ll1-1,0,l3)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,-1,ll2,l3)=il
        rxyz(1,il)=rxyz(1,i)-alat(1)
        rxyz(2,il)=rxyz(2,i)+alat(2)
        rxyz(3,il)=rxyz(3,i)
        enddo

        in=icell(0,0,ll2-1,l3)
        icell(0,ll1,-1,l3)=in
        do ii=1,in
        i=icell(ii,0,ll2-1,l3)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,ll1,-1,l3)=il
        rxyz(1,il)=rxyz(1,i)+alat(1)
        rxyz(2,il)=rxyz(2,i)-alat(2)
        rxyz(3,il)=rxyz(3,i)
        enddo

        in=icell(0,ll1-1,ll2-1,l3)
        icell(0,-1,-1,l3)=in
        do ii=1,in
        i=icell(ii,ll1-1,ll2-1,l3)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,-1,-1,l3)=il
        rxyz(1,il)=rxyz(1,i)-alat(1)
        rxyz(2,il)=rxyz(2,i)-alat(2)
        rxyz(3,il)=rxyz(3,i)
        enddo

        enddo


! corners
        in=icell(0,0,0,0)
        icell(0,ll1,ll2,ll3)=in
        do ii=1,in
        i=icell(ii,0,0,0)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,ll1,ll2,ll3)=il
        rxyz(1,il)=rxyz(1,i)+alat(1)
        rxyz(2,il)=rxyz(2,i)+alat(2)
        rxyz(3,il)=rxyz(3,i)+alat(3)
        enddo

        in=icell(0,ll1-1,0,0)
        icell(0,-1,ll2,ll3)=in
        do ii=1,in
        i=icell(ii,ll1-1,0,0)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,-1,ll2,ll3)=il
        rxyz(1,il)=rxyz(1,i)-alat(1)
        rxyz(2,il)=rxyz(2,i)+alat(2)
        rxyz(3,il)=rxyz(3,i)+alat(3)
        enddo

        in=icell(0,0,ll2-1,0)
        icell(0,ll1,-1,ll3)=in
        do ii=1,in
        i=icell(ii,0,ll2-1,0)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,ll1,-1,ll3)=il
        rxyz(1,il)=rxyz(1,i)+alat(1)
        rxyz(2,il)=rxyz(2,i)-alat(2)
        rxyz(3,il)=rxyz(3,i)+alat(3)
        enddo

        in=icell(0,ll1-1,ll2-1,0)
        icell(0,-1,-1,ll3)=in
        do ii=1,in
        i=icell(ii,ll1-1,ll2-1,0)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,-1,-1,ll3)=il
        rxyz(1,il)=rxyz(1,i)-alat(1)
        rxyz(2,il)=rxyz(2,i)-alat(2)
        rxyz(3,il)=rxyz(3,i)+alat(3)
        enddo

        in=icell(0,0,0,ll3-1)
        icell(0,ll1,ll2,-1)=in
        do ii=1,in
        i=icell(ii,0,0,ll3-1)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,ll1,ll2,-1)=il
        rxyz(1,il)=rxyz(1,i)+alat(1)
        rxyz(2,il)=rxyz(2,i)+alat(2)
        rxyz(3,il)=rxyz(3,i)-alat(3)
        enddo
        
        in=icell(0,ll1-1,0,ll3-1)
        icell(0,-1,ll2,-1)=in
        do ii=1,in
        i=icell(ii,ll1-1,0,ll3-1)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,-1,ll2,-1)=il
        rxyz(1,il)=rxyz(1,i)-alat(1)
        rxyz(2,il)=rxyz(2,i)+alat(2)
        rxyz(3,il)=rxyz(3,i)-alat(3)
        enddo

        in=icell(0,0,ll2-1,ll3-1)
        icell(0,ll1,-1,-1)=in
        do ii=1,in
        i=icell(ii,0,ll2-1,ll3-1)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,ll1,-1,-1)=il
        rxyz(1,il)=rxyz(1,i)+alat(1)
        rxyz(2,il)=rxyz(2,i)-alat(2)
        rxyz(3,il)=rxyz(3,i)-alat(3)
        enddo

        in=icell(0,ll1-1,ll2-1,ll3-1)
        icell(0,-1,-1,-1)=in
        do ii=1,in
        i=icell(ii,ll1-1,ll2-1,ll3-1)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,-1,-1,-1)=il
        rxyz(1,il)=rxyz(1,i)-alat(1)
        rxyz(2,il)=rxyz(2,i)-alat(2)
        rxyz(3,il)=rxyz(3,i)-alat(3)
        enddo



        allocate(lsta(2,nat))
        nnbrx=300
2345    nnbrx=3*nnbrx/2
        allocate(lstb(nnbrx*nat),rel(5,nnbrx*nat))

        indlstx=0


        npr=1
        iam=0

        cut2=cut**2
! assign contiguous portions of the arrays lstb and rel to the threads (this version only contains one thread)
        myspace=(nat*nnbrx)/npr  
        if (iam.eq.0) myspaceout=myspace
! Verlet list, relative positions
        indlst=0
      do 6000,l3=0,ll3-1
      do 6000,l2=0,ll2-1
      do 6000,l1=0,ll1-1
      do 6600,ii=1,icell(0,l1,l2,l3)
        iat=icell(ii,l1,l2,l3)
        if ( ((iat-1)*npr)/nat .eq. iam) then
        lsta(1,iat)=iam*myspace+indlst+1
        call sw_sublstiat_l(iat,nn,ncx,ll1,ll2,ll3,l1,l2,l3,myspace, & 
             rxyz,icell,lstb(iam*myspace+1),lay,rel(1,iam*myspace+1),cut2,indlst)
        lsta(2,iat)=iam*myspace+indlst
        ipb=lsta(1,iat)
        ndat=lsta(2,iat)-lsta(1,iat)+1
        endif

6600    continue
6000    continue
        indlstx=max(indlstx,indlst)

           if (indlstx.ge.myspaceout) then
               write(10,*) count,'NNBRX too  small', nnbrx
               deallocate(lstb,rel)
               goto 2345
           endif

        npr=1
        iam=0
        npjx=300 ; npjkx=6000
!end of creating pairlist part------------------------------------------------------------

!start energy and force calculation-------------------------------------------------------
!set all variables to zero
        p=0.0d0
        p3=0.0d0
        pv3=0.0d0
        fx(:)=0.0d0
        fy(:)=0.0d0
        fz(:)=0.0d0
        fx3(:)=0.0d0
        fy3(:)=0.0d0
        fz3(:)=0.0d0
!-----------------------------------------------------------------------------------------
!     triple loop for the 2 and 3-body forces
!     do 20 i
!     do 30 j
!     do 40 k
!     the pairlists lsta and lstb are used for the perodic boundary conditions
!-----------------------------------------------------------------------------------------

      do 20 i=1,nat       !DO 20
      call sw_subfeniat_l(i,nat,nnbrx,rel,cut2,p,p3,fx,fy,fz,fx3,fy3,fz3,lstb,lsta,isigma,sigma)
   20 continue
!-----------------------------------------------------------------------------------------
!*****if necessary,********
      do i=1,nat
        fx(i)=fx(i)+fx3(i)
        fy(i)=fy(i)+fy3(i)
        fz(i)=fz(i)+fz3(i)
      end do
!-----------------------------------------------------------------------------------------

        do i=1,nat
           fxyz(1,i)=fx(i)*esigma
           fxyz(2,i)=fy(i)*esigma
           fxyz(3,i)=fz(i)*esigma
         enddo
         etot=(p+p3)*eps
         deallocate(rxyz,icell,lay,lsta,lstb,rel)
end
!End of the force calculation-------------------------------------------------------------

    REAL(8) FUNCTION f(c)
        real(8) :: c, fm, ra, bb, aa, c4, crainv
        parameter(aa=7.049556277d0)
        parameter(ra=1.8d0,bb=0.6022245584d0)

        if((c-ra).lt.0.d0) then
        crainv=1.0d0/(c-ra)
        c4=c*c*c*c
        f=aa*bb*4.0d0/(c4*c)*dexp(crainv)+aa*(bb/(c4)-1.0d0)*dexp(crainv)*crainv*crainv
        else
        f=0.d0
        endif

        return
        end function f

     REAL(8) function pe(d)
        real(8), parameter :: aa=7.049556277d0
        real(8), parameter :: bb=0.6022245584d0
        real(8), parameter :: ra=1.8d0
        real(8) :: d

        if((d-ra).lt.0.d0) then  
          pe=aa*(bb/(d*d*d*d)-1.0d0)*dexp(1.0d0/(d-ra))
        else
          pe=0.d0
        endif
        return
      end function pe


!------------------------------------------------------------------------------------------
        subroutine sw_sublstiat_l(iat,nn,ncx,ll1,ll2,ll3,l1,l2,l3,myspace, &
                   rxyz,icell,lstb,lay,rel,cut2,indlst)
! finds the neighbours of atom iat (specified by lsta and lstb) and and
! the relative position rel of iat with respect to these neighbours
        implicit real(8) (a-h,o-z)
        integer :: iat, nn, ncx, ll1, ll2, ll3, l1, l2, l3, myspace, indlst
        integer :: jat, jj, k1, k2, k3
        integer :: lay(nn), icell(0:ncx,-1:ll1,-1:ll2,-1:ll3), lstb(0:myspace-1)
        real(8) :: rxyz(3,nn), rel(5,0:myspace-1), cut2

        do 6363,k3=l3-1,l3+1
        do 6363,k2=l2-1,l2+1
        do 6363,k1=l1-1,l1+1
        do 6363,jj=1,icell(0,k1,k2,k3)
          jat=icell(jj,k1,k2,k3)
          if (jat.eq.iat) goto 6363
          xrel= rxyz(1,iat)-rxyz(1,jat)
          yrel= rxyz(2,iat)-rxyz(2,jat)
          zrel= rxyz(3,iat)-rxyz(3,jat)
          rr2=xrel**2 + yrel**2 + zrel**2
          if ( rr2 .le. cut2 ) then
           indlst=min(indlst,myspace-1)
           lstb(indlst)=lay(jat)
!        write(6,*) 'iat,indlst,lay(jat)',iat,indlst,lay(jat)
           tt=sqrt(rr2)
           tti=1.d0/tt
           rel(1,indlst)=xrel*tti
           rel(2,indlst)=yrel*tti
           rel(3,indlst)=zrel*tti
           rel(4,indlst)=tt
           rel(5,indlst)=tti
           indlst= indlst+1
          endif
6363    continue

        return
        end


subroutine sw_subfeniat_l(i,nat,nnbrx,rel,cut2,p,p3,fx,fy,fz,fx3,fy3,fz3,lstb,lsta,isigma,sigma)
implicit none
        integer, intent(IN)   :: i,nat,nnbrx,lsta(2,nat),lstb(nnbrx*nat)
        integer               :: Ipb,Ipe,k,m,l,j,nij
        real(8), intent(IN)    :: rel(5,nnbrx*nat), cut2, isigma, sigma
        real(8), intent(INOUT) :: fx(nat), fy(nat), fz(nat), fx3(nat), fy3(nat), fz3(nat), p3, p
        real(8)                :: hi, hj, hk, hixij0, hixij1, hixij, hixik0, hixik, HKXIK1
        real(8)                :: hjxij0, hjxij1, hjxij, hkxik0, hixik1, hkxik
        real(8)                :: hjxjk0, hjxjk1, hjxjk, hkxkj0, hkxkj1, hkxkj
        real(8)                :: hiyij0, hiyij1
        real(8)                :: hiyij, hiyik0, hiyik1, hiyik, hjyij0, hjyij1, hjyij
        real(8)                :: hkyik0, hkyik1, hkyik, hjyjk0, hjyjk1, hjyjk, hkykj0, hkykj1, hkykj
        real(8)                :: HIZIJ0, HIZIJ1, HIZIJ, HIZIK0, HIZIK1, HIZIK, HJZIJ0, HJZIJ1, HJZIJ, HKZIK0, HKZIK1, HKZIK
        real(8)                :: HJZJK0, HJZJK1, HJZJK, HKZKJ0, HKZKJ1, HKZKJ, PV3
        real(8)                :: rij, xik, yik, zik, xij, yij, zij, f, pe
        real(8)                :: rik, xjk, yjk, zjk, rjk, cosjik, cosijk, cosikj, cosjik3, cosijk3, cosikj3
        real(8)                :: rija, rika, rjka, refi, refj, refk
        real(8)                :: invrik, invrjk, invrij, gam, ramda, aa, ra, bb, uc, invrika, invrjka, invrija, force, c4, crainv

        parameter(gam=1.2d0,ramda=21.0d0,aa=7.049556277d0)
        parameter(ra=1.8d0,bb=0.6022245584d0)
        parameter(uc=4.5d0)
        
         Ipb = lsta(1,i)
         Ipe = lsta(2,i)

      do 30 l=Ipb,Ipe,+1 !DO 30
            j=lstb(l)
         if (j.le.i) go to 30
         nij=0
         rij=rel(4,l)*isigma
         invrij=rel(5,l)*sigma
         xij=rel(1,l)*rij
         yij=rel(2,l)*rij
         zij=rel(3,l)*rij


        if(rij.ge.2.d0*ra) go to 30
        if(rij.ge.ra) go to 31
        
        crainv=1.0d0/(rij-ra)
        c4=rij*rij*rij*rij
        force=aa*bb*4.0d0/(c4*rij)*dexp(crainv)+aa*(bb/(c4)-1.0d0)*dexp(crainv)*crainv*crainv

        fx(i)=force*xij*invrij+fx(i)
        fy(i)=force*yij*invrij+fy(i)
        fz(i)=force*zij*invrij+fz(i)

        fx(j)=-force*xij*invrij+fx(j)
        fy(j)=-force*yij*invrij+fy(j)
        fz(j)=-force*zij*invrij+fz(j)

        p=p+aa*(bb/(rij*rij*rij*rij)-1.0d0)*dexp(1.0d0/(rij-ra))

        nij=1

   31   continue

         do 40 m=Ipb,Ipe,+1 !DO 40
         k=lstb(m)
         if (k.le.j) goto 40
         invrik=rel(5,m)*sigma
         rik=rel(4,m)*isigma
         xik=rel(1,m)*rik
         yik=rel(2,m)*rik
         zik=rel(3,m)*rik

        if((rik.ge.ra).and.(nij.eq.0)) go to 40

        xjk=xik-xij
        yjk=yik-yij
        zjk=zik-zij

        rjk=sqrt(xjk*xjk+yjk*yjk+zjk*zjk)
        invrjk=1.d0/rjk

        if((rjk.ge.ra).and.(nij.eq.0)) go to 40
        cosjik=(xij*xik+yij*yik+zij*zik)*(invrij*invrik)
        cosijk=(-xij*xjk-yij*yjk-zij*zjk)*(invrij*invrjk)
        cosikj=(xik*xjk+yik*yjk+zik*zjk)*(invrik*invrjk)
        cosjik3=cosjik+1.0d0/3.0d0
        cosijk3=cosijk+1.0d0/3.0d0
        cosikj3=cosikj+1.0d0/3.0d0

        rija=rij-ra
        rika=rik-ra
        rjka=rjk-ra

        invrija=1.d0/rija
        invrika=1.d0/rika
        invrjka=1.d0/rjka

      if(rija.ge.0.0d0) then
      refi=0.0d0
      refj=0.0d0
      refk=ramda*exp(gam*invrika+gam*invrjka)
      else if((rija.lt.0.0d0).and.(rika.lt.0.0d0)) then
       if(rjka.lt.0.0d0) then
       refi=ramda*exp(gam*invrija+gam*invrika)
       refj=ramda*exp(gam*invrija+gam*invrjka)
       refk=ramda*exp(gam*invrika+gam*invrjka)
       else
       refi=ramda*exp(gam*invrija+gam*invrika)
       refj=0.0d0
       refk=0.0d0
       end if
      else if((rija.lt.0.0d0).and.(rjka.lt.0.0d0)) then
      refi=0.0d0
      refj=ramda*exp(gam*invrija+gam*invrjka)
      refk=0.0d0
else
goto 40
      refi=0.0d0
      refj=0.0d0
      refk=0.0d0
      end if

      hi=refi*cosjik3*cosjik3
      hj=refj*cosijk3*cosijk3
      hk=refk*cosikj3*cosikj3
      p3=p3+hi+hj+hk

      hixij0=2.0d0*(xik*invrik-xij*cosjik*invrij)
      hixij1=gam*xij*cosjik3*(invrija*invrija)
      hixij=refi*cosjik3*(hixij0-hixij1)*invrij
      hixik0=2.0d0*(xij*invrij-xik*cosjik*invrik)
      hixik1=gam*xik*cosjik3*(invrika*invrika)
      hixik=refi*cosjik3*(hixik0-hixik1)*invrik
      hjxij0=2.0d0*(-xjk*invrjk-xij*cosijk*invrij)
      hjxij1=gam*xij*cosijk3*(invrija*invrija)
      hjxij=refj*cosijk3*(hjxij0-hjxij1)*invrij
      hkxik0=2.0d0*(xjk*invrjk-xik*cosikj*invrik)
      hkxik1=gam*xik*cosikj3*(invrika*invrika)
      hkxik=refk*cosikj3*(hkxik0-hkxik1)*invrik
      hjxjk0=2.0d0*(-xij*invrij-xjk*cosijk*invrjk)
      hjxjk1=gam*xjk*cosijk3*(invrjka*invrjka)
      hjxjk=refj*cosijk3*(hjxjk0-hjxjk1)*invrjk
      hkxkj0=2.0d0*(-xik*invrik+xjk*cosikj*invrjk)
      hkxkj1=gam*xjk*cosikj3*(invrjka*invrjka)
      hkxkj=refk*cosikj3*(hkxkj0+hkxkj1)*invrjk


      hiyij0=2.0d0*(yik*invrik-yij*cosjik*invrij)
      hiyij1=gam*yij*cosjik3*(invrija*invrija)
      hiyij=refi*cosjik3*(hiyij0-hiyij1)*invrij
      hiyik0=2.0d0*(yij*invrij-yik*cosjik*invrik)
      hiyik1=gam*yik*cosjik3*(invrika*invrika)
      hiyik=refi*cosjik3*(hiyik0-hiyik1)*invrik
      hjyij0=2.0d0*(-yjk*invrjk-yij*cosijk*invrij)
      hjyij1=gam*yij*cosijk3*(invrija*invrija)
      hjyij=refj*cosijk3*(hjyij0-hjyij1)*invrij
      hkyik0=2.0d0*(yjk*invrjk-yik*cosikj*invrik)
      hkyik1=gam*yik*cosikj3*(invrika*invrika)
      hkyik=refk*cosikj3*(hkyik0-hkyik1)*invrik
      hjyjk0=2.0d0*(-yij*invrij-yjk*cosijk*invrjk)
      hjyjk1=gam*yjk*cosijk3*(invrjka*invrjka)
      hjyjk=refj*cosijk3*(hjyjk0-hjyjk1)*invrjk
      hkykj0=2.0d0*(-yik*invrik+yjk*cosikj*invrjk)
      hkykj1=gam*yjk*cosikj3*(invrjka*invrjka)
      hkykj=refk*cosikj3*(hkykj0+hkykj1)*invrjk

      hizij0=2.0d0*(zik*invrik-zij*cosjik*invrij)
      hizij1=gam*zij*cosjik3*(invrija*invrija)
      hizij=refi*cosjik3*(hizij0-hizij1)*invrij
      hizik0=2.0d0*(zij*invrij-zik*cosjik*invrik)
      hizik1=gam*zik*cosjik3*(invrika*invrika)
      hizik=refi*cosjik3*(hizik0-hizik1)*invrik
      hjzij0=2.0d0*(-zjk*invrjk-zij*cosijk*invrij)
      hjzij1=gam*zij*cosijk3*(invrija*invrija)
      hjzij=refj*cosijk3*(hjzij0-hjzij1)*invrij
      hkzik0=2.0d0*(zjk*invrjk-zik*cosikj*invrik)
      hkzik1=gam*zik*cosikj3*(invrika*invrika)
      hkzik=refk*cosikj3*(hkzik0-hkzik1)*invrik
      hjzjk0=2.0d0*(-zij*invrij-zjk*cosijk*invrjk)
      hjzjk1=gam*zjk*cosijk3*(invrjka*invrjka)
      hjzjk=refj*cosijk3*(hjzjk0-hjzjk1)*invrjk
      hkzkj0=2.0d0*(-zik*invrik+zjk*cosikj*invrjk)
      hkzkj1=gam*zjk*cosikj3*(invrjka*invrjka)
      hkzkj=refk*cosikj3*(hkzkj0+hkzkj1)*invrjk

        fx3(i)=fx3(i)-hixij-hixik-hjxij-hkxik
        fy3(i)=fy3(i)-hiyij-hiyik-hjyij-hkyik
        fz3(i)=fz3(i)-hizij-hizik-hjzij-hkzik

        fx3(j)=fx3(j)+hixij+hjxij-hjxjk+hkxkj
        fy3(j)=fy3(j)+hiyij+hjyij-hjyjk+hkykj
        fz3(j)=fz3(j)+hizij+hjzij-hjzjk+hkzkj

        fx3(k)=fx3(k)+hixik+hkxik-hkxkj+hjxjk
        fy3(k)=fy3(k)+hiyik+hkyik-hkykj+hjyjk
        fz3(k)=fz3(k)+hizik+hkzik-hkzkj+hjzjk
   40   continue
   30   continue
end subroutine

subroutine eip_tersoff_silicon(nat,alat,rxyz,fxyz,etot,count)
!*****************************************************************************************
! This subroutine evaluates the Tersoff Silicon potential with linear scaling
! COPYRIGHT
!    Copyright (C) 2009 AIST, UNIBAS
!    This file is distributed under the terms of the
!    GNU General Public License, see 
!    http://www.gnu.org/copyleft/gpl.txt .
!
! Implementation: Original version was written by Kengo Nishio, AIST Tsukuba (JP) 
!                 Improved by M. Amsler, S. Goedecker, Basel University (CH), 2009
!
! Note:
!     Parameters and functional form from PRL 61, 2879 (1988) and PRB 39, 5566 (1989)
!     
! Input:
!     nat, integer: the number of atoms
!     alat,  real(8), dim(3)   : the three edges of the orthoromic simulation cell, periodic boundaries are applied
!                                and atoms outside the cell will be brought back into the box
!     rxyz, real(8), dim(3,nat) : the xyz cartesian components of the atomic positions in Angstroem    
!                                
! Output:
!     fxyz,  real(8), dim(3,nat): the xyz cartesian forces in on the corresponding atomic components n eV/A  
!     etot, real(8)            : total potential energy, 2-body and 3-body, in eV
!     count, real(8)           : increased by 1.d0 at each call of this subroutine,
!                                needs to be initialized to 0.d0 before calling this routine for the first time
!*****************************************************************************************
      implicit none
      integer   :: nat, iat
      real(8)   :: etot,count
      real(8)                                   :: rxyz(3,nat),fxyz(3,nat), alat(3)
      real(8),              dimension(1:2,1:2)  :: R1, R2
      real(8),              dimension(1:2,1:2)  :: Cr, Ca, alr, ala, X
      real(8),              dimension(1:2)      :: Pmass
      integer                                   :: Nmol, Npmax, NNmax
      real(8)                                   :: xbox, ybox, zbox
      integer,allocatable,  dimension(:)        :: Kinds
      real(8),allocatable,  dimension(:)        :: XYZRrefdf
      real(8),allocatable,  dimension(:)        :: UadUrdf
      real(8)                                   :: Urtot
      real(8),              dimension(1:2)      :: Pn, Co_bcd, bcsq, dsq, h
      real(8)                                   :: Uatot
      real(8),allocatable,  dimension(:)        :: dkEij
      integer,allocatable,  dimension(:,:)      :: lsta
      integer,allocatable,  dimension(:)        :: lstb
      integer:: nnbrx,nnbrxt
      integer i,j
      integer Nptot, NNtot
      real(8) pi
      real(8) Xi, Yi, Zi
      real(8) Xij, Yij, Zij, Rij, Rreij
      real(8) R1ij, R2ij
      real(8) PL1, PL2

      real(8) Crij, Caij, alrij, alaij
      real(8) fij, dfij
      real(8) xhalf, yhalf, zhalf

      real(8) COSijk, Gi, dGi, fdG, fdGcos
      real(8) Co_pa, Co_mb1, Co_mb2
      real(8) Co_cdi, Co_hcosi, Co_dhcosi
      real(8) Pni, bcsqi, dsqi, hi

      real(8) Bij, Eij
      real(8) djEij
      real(8) dXjEij2, dYjEij2, dZjEij2
      real(8) Co1_dkEij, Co2_dkEij

      real(8), parameter :: C_Cr = 1.3936d3    !Parameter for carbon, not used in this version       
      real(8), parameter :: C_Ca = 3.4674d2    !Parameter for carbon, not used in this version
      real(8), parameter :: C_alr= 3.4879d0    !Parameter for carbon, not used in this version
      real(8), parameter :: C_ala= 2.2119d0    !Parameter for carbon, not used in this version
      real(8), parameter :: C_b  = 1.5724d-7   !Parameter for carbon, not used in this version
      real(8), parameter :: C_n  = 7.2751d-1   !Parameter for carbon, not used in this version
      real(8), parameter :: C_c  = 3.8049d4    !Parameter for carbon, not used in this version
      real(8), parameter :: C_d  = 4.3484d0    !Parameter for carbon, not used in this version
      real(8), parameter :: C_h  =-5.7058d-1   !Parameter for carbon, not used in this version
      real(8), parameter :: C_R1 = 1.8d0       !Parameter for carbon, not used in this version 
      real(8), parameter :: C_R2 = 2.1d0       !Parameter for carbon, not used in this version 
      real(8), parameter :: C_mass = 12.0d0    !Parameter for carbon, not used in this version 

      real(8), parameter :: Si_Cr = 1.8308d3
      real(8), parameter :: Si_Ca = 4.7118d2
      real(8), parameter :: Si_alr= 2.4799d0
      real(8), parameter :: Si_ala= 1.7322d0
      real(8), parameter :: Si_b  = 1.1000d-6
      real(8), parameter :: Si_n  = 7.8734d-1
      real(8), parameter :: Si_c  = 1.0039d5
      real(8), parameter :: Si_d  = 1.6217d1
      real(8), parameter :: Si_h  =-5.9825d-1
      real(8), parameter :: Si_R1 = 2.7d0
      real(8), parameter :: Si_R2 = 3.3d0       !Increased Cutoff, originally 3.0d0
      real(8), parameter :: Si_mass = 28.0855d0

        count=count+1.d0

        do iat=1,nat
        rxyz(1,iat)=modulo(modulo(rxyz(1,iat),alat(1)),alat(1))
        rxyz(2,iat)=modulo(modulo(rxyz(2,iat),alat(2)),alat(2))
        rxyz(3,iat)=modulo(modulo(rxyz(3,iat),alat(3)),alat(3))
        enddo

        allocate(Kinds(1:nat))

        nnbrx=24
        nnbrxt=3*nnbrx/2
        nnmax=nnbrxt*nat
        npmax=nnbrxt*nat
        allocate(lsta(2,nat),lstb(nnbrxt*nat))
        allocate(XYZRrefdf(1:6*Npmax), UadUrdf(1:3*Npmax), dkEij(1:3*NNmax))

        do i=1,nat
          kinds(i)=2                             !Since all atoms are Si, all of kind 2
        enddo
        fxyz=0.0d0
        xbox=alat(1); ybox=alat(2); zbox=alat(3)
        call tersoff_parameters(R1,R2,Cr,Ca,alr,ala,X,Pn,Co_bcd,bcsq,dsq,h,Pmass)
        call tersoff_pairlist_energy_forces(nat,Npmax,NNmax,xbox,ybox,zbox,Kinds,rxyz,R1,R2,Cr,Ca,alr,ala,X, XYZRrefdf,UadUrdf,Urtot,lsta,lstb,nnbrx,&
                           Pn,Co_bcd,bcsq,dsq,h,fxyz,Uatot,dkEij,count)
        etot=Urtot+Uatot
        deallocate(Kinds,XYZRrefdf,UadUrdf,dkEij,lsta,lstb)
        end
!-----------------------------------------------------------------------------------------
      Subroutine tersoff_parameters(R1,R2,Cr,Ca,alr,ala,X,Pn,Co_bcd,bcsq,dsq,h,Pmass)
      implicit none

      real(8), intent(out), dimension(1:2,1:2) :: R1, R2
      real(8), intent(out), dimension(1:2,1:2) :: Cr, Ca, alr, ala, X
      real(8), intent(out), dimension(1:2) :: Pn, Co_bcd, bcsq, dsq, h
      real(8), intent(out), dimension(1:2) :: Pmass

      real(8), parameter :: C_Cr = 1.3936d3    !Parameter for carbon, not used in this version       
      real(8), parameter :: C_Ca = 3.4674d2    !Parameter for carbon, not used in this version
      real(8), parameter :: C_alr= 3.4879d0    !Parameter for carbon, not used in this version
      real(8), parameter :: C_ala= 2.2119d0    !Parameter for carbon, not used in this version
      real(8), parameter :: C_b  = 1.5724d-7   !Parameter for carbon, not used in this version
      real(8), parameter :: C_n  = 7.2751d-1   !Parameter for carbon, not used in this version
      real(8), parameter :: C_c  = 3.8049d4    !Parameter for carbon, not used in this version
      real(8), parameter :: C_d  = 4.3484d0    !Parameter for carbon, not used in this version
      real(8), parameter :: C_h  =-5.7058d-1   !Parameter for carbon, not used in this version
      real(8), parameter :: C_R1 = 1.8d0       !Parameter for carbon, not used in this version 
      real(8), parameter :: C_R2 = 2.1d0       !Parameter for carbon, not used in this version 
      real(8), parameter :: C_mass = 12.0d0    !Parameter for carbon, not used in this version 

      real(8), parameter :: Si_Cr = 1.8308d3
      real(8), parameter :: Si_Ca = 4.7118d2
      real(8), parameter :: Si_alr= 2.4799d0
      real(8), parameter :: Si_ala= 1.7322d0
      real(8), parameter :: Si_b  = 1.1000d-6
      real(8), parameter :: Si_n  = 7.8734d-1
      real(8), parameter :: Si_c  = 1.0039d5
      real(8), parameter :: Si_d  = 1.6217d1
      real(8), parameter :: Si_h  =-5.9825d-1
      real(8), parameter :: Si_R1 = 2.7d0
      real(8), parameter :: Si_R2 = 3.3d0      !Increased Cutoff, originally 3.0d0
      real(8), parameter :: Si_mass = 28.0855d0


      Cr(1,1) =  C_Cr
      Cr(2,2) = Si_Cr
      Cr(1,2) = dsqrt( Cr(1,1)*Cr(2,2) )
      Cr(2,1) = Cr(1,2)

      Ca(1,1) =  C_Ca
      Ca(2,2) = Si_Ca
      Ca(1,2) = dsqrt( Ca(1,1)*Ca(2,2) )
      Ca(2,1) = Ca(1,2)

      R1(1,1) =  C_R1
      R1(2,2) = Si_R1
      R1(1,2) = dsqrt( R1(1,1)*R1(2,2) )
      R1(2,1) = R1(1,2)

      R2(1,1) =  C_R2
      R2(2,2) = Si_R2
      R2(1,2) = dsqrt( R2(1,1)*R2(2,2) )
      R2(2,1) = R2(1,2)

      X(1,1)  = 1.0d0
      X(2,2)  = 1.0d0
      X(1,2)  = 0.9776d0
      X(2,1)  = 0.9776d0

      alr(1,1) =  C_alr
      alr(2,2) = Si_alr
      alr(1,2) = 0.5d0*( alr(1,1)+alr(2,2) )
      alr(2,1) = alr(1,2)

      ala(1,1) =  C_ala
      ala(2,2) = Si_ala
      ala(1,2) = 0.5d0*( ala(1,1)+ala(2,2) )
      ala(2,1) = ala(1,2)

      Pn(1) =  C_n
      Pn(2) = Si_n

      Co_bcd(1) =  C_b*( 1.0d0 + C_c*C_c  /(C_d*C_d)   )
      Co_bcd(2) = Si_b*( 1.0d0 + Si_c*Si_c/(Si_d*Si_d) )

      bcsq(1) = C_b  * C_c  * C_c
      bcsq(2) = Si_b * Si_c * Si_c

      dsq(1) =  C_d *  C_d
      dsq(2) = Si_d * Si_d

      h(1) =  C_h
      h(2) = Si_h

      Pmass(1)=C_mass
      Pmass(2)=Si_mass


      RETURN

      STOP
      END
!-----------------------------------------------------------------------------------------
      Subroutine tersoff_pairlist_energy_forces(Nmol,Npmax,NNmax,xbox,ybox,zbox,Kinds,R,R1,R2,Cr,Ca,alr,ala,X,XYZRrefdf,UadUrdf,Urtot,lsta,lstb,nnbrx,&
                           Pn,Co_bcd,bcsq,dsq,h,F,Uatot,dkEij,count)
        implicit real(8) (a-h,o-z)
        integer, intent(in)                       :: Nmol, Npmax, NNmax
        integer, intent(inout)                    :: nnbrx
        real(8), intent(in)                       :: xbox, ybox, zbox,count
        integer, intent(in),  dimension(1:Nmol)   :: Kinds
        real(8), intent(in),  dimension(1:3*Nmol) :: R
        real(8), intent(in),  dimension(1:2,1:2)  :: R1, R2 
        real(8), intent(in),  dimension(1:2,1:2)  :: Cr, Ca, alr, ala
        real(8), intent(in),  dimension(1:2,1:2)  :: X
        real(8), intent(in),  dimension(1:2)      ::Pn, Co_bcd, bcsq, dsq, h
        real(8), intent(out), dimension(1:6*Npmax):: XYZRrefdf
        real(8), intent(out), dimension(1:3*Npmax):: UadUrdf
        real(8), intent(out)                      :: Urtot
        real(8), intent(out), dimension(1:3*Nmol) :: F
        real(8), intent(out)                      :: Uatot
        
        integer i,j
        integer :: iam, iat, ii, il, in, indlst, indlstx, istopg
        integer :: jat, l1, l2, l3, laymx, ll1, ll2, ll3
        integer :: myspace, myspaceout, nat, ncx, ndat, nn, npjkx, npjx, npr
        integer Ipt3,Jpt3,Nppt3,Nppt6
        integer Ki,Kj
        integer Nptot, NNtot
        integer :: ij, ijpt3, ijpt6, ik, ikpt6
        integer Ipb, Ipe
        integer Kpt3
        integer Nkpt3
        integer                                  :: lsta(2,Nmol),lstb(nnbrx*Nmol)
        integer, ALLOCATABLE, DIMENSION(:)       :: lay
        integer, ALLOCATABLE, DIMENSION(:,:,:,:) :: icell
  
        real(8) Ur, Ua
        real(8) pi
        real(8) Xi, Yi, Zi
        real(8) Xij, Yij, Zij, Rij
        real(8) R1ij, R2ij
        real(8) PL1, PL2
        real(8) Crij, Caij, alrij, alaij
        real(8) fij, dfij
        real(8) xhalf, yhalf, zhalf
        real(8) XRreik, YRreik,  ZRreik, Rreik, fik, dfik
        real(8) COSijk, Gi, dGi, fdG, fdGcos
        real(8) Co_pa, Co_mb1, Co_mb2
        real(8) Co_cdi, Co_hcosi, Co_dhcosi
        real(8) Pni, bcsqi, dsqi, hi
        real(8) Bij, Eij
        real(8) djEij
        real(8) dXjEij2, dYjEij2, dZjEij2
        real(8) Co1_dkEij, Co2_dkEij
        real(8) dFxi, dFyi, dFzi
        real(8) dFxj, dFyj, dFzj
        real(8) dFxk, dFyk, dFzk
        real(8), dimension(1:3*NNmax)::dkEij
        real(8), ALLOCATABLE, DIMENSION(:,:)     :: rxyz
        real(8), ALLOCATABLE, DIMENSION(:,:)     :: rel
        real(8), ALLOCATABLE, DIMENSION(:,:)     :: txyz
        real(8), ALLOCATABLE, DIMENSION(:,:)     :: f2ij,f3ij,f3ik
        dimension rxyz0(3,Nmol),fxyz(3,Nmol),alat(3)


      pi=dacos(-1.0d0)

      xhalf=0.5d0*xbox
      yhalf=0.5d0*ybox
      zhalf=0.5d0*zbox
      
      nat=Nmol
      
      alat(1)=xbox
      alat(2)=ybox
      alat(3)=zbox
      
      
      do iat=1,nat
      jat=3*(iat-1)
      rxyz0(1,iat)=R(jat+1)
      rxyz0(2,iat)=R(jat+2)
      rxyz0(3,iat)=R(jat+3)        
      enddo

      cut=R2(2,2)-1.d-9


! linear scaling calculation of verlet list
        ll1=int(alat(1)/cut)
        if (ll1.lt.1) stop 'alat(1) too small'
        ll2=int(alat(2)/cut)
        if (ll2.lt.1) stop 'alat(2) too small'
        ll3=int(alat(3)/cut)
        if (ll3.lt.1) stop 'alat(3) too small'

! determine number of threadsi (this version is only singlethreaded)
        npr=1
! linear scaling calculation of verlet list

        ncx=8
1234    ncx=ncx*2
        allocate(icell(0:ncx,-1:ll1,-1:ll2,-1:ll3))
        do 984,l3=-1,ll3
        do 984,l2=-1,ll2
        do 984,l1=-1,ll1
984     icell(0,l1,l2,l3)=0
        rlc1i=ll1/alat(1)
        rlc2i=ll2/alat(2)
        rlc3i=ll3/alat(3)

        do 983,iat=1,nat
        l1=int(rxyz0(1,iat)*rlc1i)
        l2=int(rxyz0(2,iat)*rlc2i)
        l3=int(rxyz0(3,iat)*rlc3i)

        ii=icell(0,l1,l2,l3)
        ii=ii+1
        icell(0,l1,l2,l3)=ii
        if (ii.gt.ncx) then
        write(10,*) count,'NCX too small',ncx
        deallocate(icell)
        goto 1234
        endif
        icell(ii,l1,l2,l3)=iat
983     continue

! duplicate all atoms within boundary layer
        laymx=ncx*(2*ll1*ll2+2*ll1*ll3+2*ll2*ll3+4*ll1+4*ll2+4*ll3+8)
        nn=nat+laymx
        allocate(rxyz(3,nn),lay(nn))
        do  iat=1,nat
        lay(iat)=iat
        rxyz(1,iat)=rxyz0(1,iat)
        rxyz(2,iat)=rxyz0(2,iat)
        rxyz(3,iat)=rxyz0(3,iat)
        enddo
        il=nat
! xy plane
        do l2=0,ll2-1
        do l1=0,ll1-1

        in=icell(0,l1,l2,0)
        icell(0,l1,l2,ll3)=in
        do ii=1,in
        i=icell(ii,l1,l2,0)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,l1,l2,ll3)=il
        rxyz(1,il)=rxyz(1,i)
        rxyz(2,il)=rxyz(2,i)
        rxyz(3,il)=rxyz(3,i)+alat(3)
        enddo

        in=icell(0,l1,l2,ll3-1)
        icell(0,l1,l2,-1)=in
        do ii=1,in
        i=icell(ii,l1,l2,ll3-1)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,l1,l2,-1)=il
        rxyz(1,il)=rxyz(1,i)
        rxyz(2,il)=rxyz(2,i)
        rxyz(3,il)=rxyz(3,i)-alat(3)
        enddo

        enddo
        enddo


! yz plane
        do l3=0,ll3-1
        do l2=0,ll2-1

        in=icell(0,0,l2,l3)
        icell(0,ll1,l2,l3)=in
        do ii=1,in
        i=icell(ii,0,l2,l3)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,ll1,l2,l3)=il
        rxyz(1,il)=rxyz(1,i)+alat(1)
        rxyz(2,il)=rxyz(2,i)
        rxyz(3,il)=rxyz(3,i)
        enddo

        in=icell(0,ll1-1,l2,l3)
        icell(0,-1,l2,l3)=in
        do ii=1,in
        i=icell(ii,ll1-1,l2,l3)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,-1,l2,l3)=il
        rxyz(1,il)=rxyz(1,i)-alat(1)
        rxyz(2,il)=rxyz(2,i)
        rxyz(3,il)=rxyz(3,i)
        enddo

        enddo
        enddo


! xz plane
        do l3=0,ll3-1
        do l1=0,ll1-1

        in=icell(0,l1,0,l3)
        icell(0,l1,ll2,l3)=in
        do ii=1,in
        i=icell(ii,l1,0,l3)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,l1,ll2,l3)=il
        rxyz(1,il)=rxyz(1,i)
        rxyz(2,il)=rxyz(2,i)+alat(2)
        rxyz(3,il)=rxyz(3,i)
        enddo

        in=icell(0,l1,ll2-1,l3)
        icell(0,l1,-1,l3)=in
        do ii=1,in
        i=icell(ii,l1,ll2-1,l3)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,l1,-1,l3)=il
        rxyz(1,il)=rxyz(1,i)
        rxyz(2,il)=rxyz(2,i)-alat(2)
        rxyz(3,il)=rxyz(3,i)
        enddo

        enddo
        enddo


! x axis
        do l1=0,ll1-1

        in=icell(0,l1,0,0)
        icell(0,l1,ll2,ll3)=in
        do ii=1,in
        i=icell(ii,l1,0,0)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,l1,ll2,ll3)=il
        rxyz(1,il)=rxyz(1,i)
        rxyz(2,il)=rxyz(2,i)+alat(2)
        rxyz(3,il)=rxyz(3,i)+alat(3)
        enddo

        in=icell(0,l1,0,ll3-1)
        icell(0,l1,ll2,-1)=in
        do ii=1,in
        i=icell(ii,l1,0,ll3-1)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,l1,ll2,-1)=il
        rxyz(1,il)=rxyz(1,i)
        rxyz(2,il)=rxyz(2,i)+alat(2)
        rxyz(3,il)=rxyz(3,i)-alat(3)
        enddo

        in=icell(0,l1,ll2-1,0)
        icell(0,l1,-1,ll3)=in
        do ii=1,in
        i=icell(ii,l1,ll2-1,0)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,l1,-1,ll3)=il
        rxyz(1,il)=rxyz(1,i)
        rxyz(2,il)=rxyz(2,i)-alat(2)
        rxyz(3,il)=rxyz(3,i)+alat(3)
        enddo

        in=icell(0,l1,ll2-1,ll3-1)
        icell(0,l1,-1,-1)=in
        do ii=1,in
        i=icell(ii,l1,ll2-1,ll3-1)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,l1,-1,-1)=il
        rxyz(1,il)=rxyz(1,i)
        rxyz(2,il)=rxyz(2,i)-alat(2)
        rxyz(3,il)=rxyz(3,i)-alat(3)
        enddo

        enddo


! y axis
        do l2=0,ll2-1

        in=icell(0,0,l2,0)
        icell(0,ll1,l2,ll3)=in
        do ii=1,in
        i=icell(ii,0,l2,0)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,ll1,l2,ll3)=il
        rxyz(1,il)=rxyz(1,i)+alat(1)
        rxyz(2,il)=rxyz(2,i)
        rxyz(3,il)=rxyz(3,i)+alat(3)
        enddo

        in=icell(0,0,l2,ll3-1)
        icell(0,ll1,l2,-1)=in
        do ii=1,in
        i=icell(ii,0,l2,ll3-1)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,ll1,l2,-1)=il
        rxyz(1,il)=rxyz(1,i)+alat(1)
        rxyz(2,il)=rxyz(2,i)
        rxyz(3,il)=rxyz(3,i)-alat(3)
        enddo

        in=icell(0,ll1-1,l2,0)
        icell(0,-1,l2,ll3)=in
        do ii=1,in
        i=icell(ii,ll1-1,l2,0)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,-1,l2,ll3)=il
        rxyz(1,il)=rxyz(1,i)-alat(1)
        rxyz(2,il)=rxyz(2,i)
        rxyz(3,il)=rxyz(3,i)+alat(3)
        enddo

        in=icell(0,ll1-1,l2,ll3-1)
        icell(0,-1,l2,-1)=in
        do ii=1,in
        i=icell(ii,ll1-1,l2,ll3-1)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,-1,l2,-1)=il
        rxyz(1,il)=rxyz(1,i)-alat(1)
        rxyz(2,il)=rxyz(2,i)
        rxyz(3,il)=rxyz(3,i)-alat(3)
        enddo

        enddo


! z axis
        do l3=0,ll3-1

        in=icell(0,0,0,l3)
        icell(0,ll1,ll2,l3)=in
        do ii=1,in
        i=icell(ii,0,0,l3)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,ll1,ll2,l3)=il
        rxyz(1,il)=rxyz(1,i)+alat(1)
        rxyz(2,il)=rxyz(2,i)+alat(2)
        rxyz(3,il)=rxyz(3,i)
        enddo

        in=icell(0,ll1-1,0,l3)
        icell(0,-1,ll2,l3)=in
        do ii=1,in
        i=icell(ii,ll1-1,0,l3)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,-1,ll2,l3)=il
        rxyz(1,il)=rxyz(1,i)-alat(1)
        rxyz(2,il)=rxyz(2,i)+alat(2)
        rxyz(3,il)=rxyz(3,i)
        enddo

        in=icell(0,0,ll2-1,l3)
        icell(0,ll1,-1,l3)=in
        do ii=1,in
        i=icell(ii,0,ll2-1,l3)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,ll1,-1,l3)=il
        rxyz(1,il)=rxyz(1,i)+alat(1)
        rxyz(2,il)=rxyz(2,i)-alat(2)
        rxyz(3,il)=rxyz(3,i)
        enddo

        in=icell(0,ll1-1,ll2-1,l3)
        icell(0,-1,-1,l3)=in
        do ii=1,in
        i=icell(ii,ll1-1,ll2-1,l3)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,-1,-1,l3)=il
        rxyz(1,il)=rxyz(1,i)-alat(1)
        rxyz(2,il)=rxyz(2,i)-alat(2)
        rxyz(3,il)=rxyz(3,i)
        enddo

        enddo


! corners
        in=icell(0,0,0,0)
        icell(0,ll1,ll2,ll3)=in
        do ii=1,in
        i=icell(ii,0,0,0)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,ll1,ll2,ll3)=il
        rxyz(1,il)=rxyz(1,i)+alat(1)
        rxyz(2,il)=rxyz(2,i)+alat(2)
        rxyz(3,il)=rxyz(3,i)+alat(3)
        enddo

        in=icell(0,ll1-1,0,0)
        icell(0,-1,ll2,ll3)=in
        do ii=1,in
        i=icell(ii,ll1-1,0,0)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,-1,ll2,ll3)=il
        rxyz(1,il)=rxyz(1,i)-alat(1)
        rxyz(2,il)=rxyz(2,i)+alat(2)
        rxyz(3,il)=rxyz(3,i)+alat(3)
        enddo

        in=icell(0,0,ll2-1,0)
        icell(0,ll1,-1,ll3)=in
        do ii=1,in
        i=icell(ii,0,ll2-1,0)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,ll1,-1,ll3)=il
        rxyz(1,il)=rxyz(1,i)+alat(1)
        rxyz(2,il)=rxyz(2,i)-alat(2)
        rxyz(3,il)=rxyz(3,i)+alat(3)
        enddo

        in=icell(0,ll1-1,ll2-1,0)
        icell(0,-1,-1,ll3)=in
        do ii=1,in
        i=icell(ii,ll1-1,ll2-1,0)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,-1,-1,ll3)=il
        rxyz(1,il)=rxyz(1,i)-alat(1)
        rxyz(2,il)=rxyz(2,i)-alat(2)
        rxyz(3,il)=rxyz(3,i)+alat(3)
        enddo

        in=icell(0,0,0,ll3-1)
        icell(0,ll1,ll2,-1)=in
        do ii=1,in
        i=icell(ii,0,0,ll3-1)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,ll1,ll2,-1)=il
        rxyz(1,il)=rxyz(1,i)+alat(1)
        rxyz(2,il)=rxyz(2,i)+alat(2)
        rxyz(3,il)=rxyz(3,i)-alat(3)
        enddo
        
        in=icell(0,ll1-1,0,ll3-1)
        icell(0,-1,ll2,-1)=in
        do ii=1,in
        i=icell(ii,ll1-1,0,ll3-1)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,-1,ll2,-1)=il
        rxyz(1,il)=rxyz(1,i)-alat(1)
        rxyz(2,il)=rxyz(2,i)+alat(2)
        rxyz(3,il)=rxyz(3,i)-alat(3)
        enddo

        in=icell(0,0,ll2-1,ll3-1)
        icell(0,ll1,-1,-1)=in
        do ii=1,in
        i=icell(ii,0,ll2-1,ll3-1)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,ll1,-1,-1)=il
        rxyz(1,il)=rxyz(1,i)+alat(1)
        rxyz(2,il)=rxyz(2,i)-alat(2)
        rxyz(3,il)=rxyz(3,i)-alat(3)
        enddo

        in=icell(0,ll1-1,ll2-1,ll3-1)
        icell(0,-1,-1,-1)=in
        do ii=1,in
        i=icell(ii,ll1-1,ll2-1,ll3-1)
        il=il+1
        if (il.gt.nn) stop 'enlarge laymx'
        lay(il)=i
        icell(ii,-1,-1,-1)=il
        rxyz(1,il)=rxyz(1,i)-alat(1)
        rxyz(2,il)=rxyz(2,i)-alat(2)
        rxyz(3,il)=rxyz(3,i)-alat(3)
        enddo

        nnbrx=3*nnbrx/2
        allocate(rel(5,nnbrx*nat))
        indlstx=0

        npr=1
        iam=0

        cut2=cut**2
! assign contiguous portions of the arrays lstb and rel to the threads
        myspace=(nat*nnbrx)/npr  
        if (iam.eq.0) myspaceout=myspace
! Verlet list, relative positions
        indlst=0
      do 6000,l3=0,ll3-1
      do 6000,l2=0,ll2-1
      do 6000,l1=0,ll1-1
      do 6600,ii=1,icell(0,l1,l2,l3)
        iat=icell(ii,l1,l2,l3)
        if ( ((iat-1)*npr)/nat .eq. iam) then
!       write(6,*) 'sublstiat:iam,iat',iam,iat
        lsta(1,iat)=iam*myspace+indlst+1
        call tersoff_sublstiat_l(iat,nn,ncx,ll1,ll2,ll3,l1,l2,l3,myspace, & 
             rxyz,icell,lstb(iam*myspace+1),lay,rel(1,iam*myspace+1),cut2,indlst)
        lsta(2,iat)=iam*myspace+indlst
        ipb=lsta(1,iat)
        ndat=lsta(2,iat)-lsta(1,iat)+1
        endif

6600    continue
6000    continue
        indlstx=max(indlstx,indlst)

           if (indlstx.ge.myspaceout) stop 'NNBRX too small'
        npr=1
        iam=0

        npjx=300 ; npjkx=6000
        istopg=0
!end of creating pairlist part------------------------------------------------------------
!Energy-----------------------------------------------------------------------------------
      Urtot=0.0d0
      Nptot=0
      
      F=0.0d0
      Uatot=0.0d0
      DO_I:DO i=1, Nmol
         call tersoff_subeniat_l(i,Nmol,Npmax,NNmax,Kinds,X,R,R1,R2,Cr,Ca,alr,ala,XYZRrefdf,UadUrdf,Urtot,lsta,lstb,nnbrx,rel,pi)
      END DO DO_I
      
      Urtot = 0.5d0*Urtot
!Force------------------------------------------------------------------------------------
      F=0.0d0
      Uatot=0.0d0

      DO_If:DO i=1, Nmol
         call tersoff_subfiat_l(i,Nmol,Npmax,NNmax,Kinds,Pn,Co_bcd,bcsq,dsq,h,XYZRrefdf,UadUrdf,F,Uatot,dkEij,lsta,lstb,nnbrx)
      END DO DO_If

      F     = 0.5d0*F
      Uatot = 0.5d0*Uatot
!-----------------------------------------------------------------------------------------
      deallocate(rxyz,icell,lay,rel)
      RETURN
      STOP 
      END


!-----------------------------------------------------------------------------------------
        subroutine tersoff_sublstiat_l(iat,nn,ncx,ll1,ll2,ll3,l1,l2,l3,myspace, & 
                   rxyz,icell,lstb,lay,rel,cut2,indlst)
! finds the neighbours of atom iat (specified by lsta and lstb) and and
! the relative position rel of iat with respect to these neighbours
        implicit real(8) (a-h,o-z)
        integer :: iat, nn, ncx, ll1, ll2, ll3, l1, l2, l3, myspace, indlst
        integer :: jat, jj, k1, k2, k3
        integer :: lay(nn), icell(0:ncx,-1:ll1,-1:ll2,-1:ll3), lstb(0:myspace-1)
        real(8) :: rxyz(3,nn), rel(5,0:myspace-1), cut2

        do 6363,k3=l3-1,l3+1
        do 6363,k2=l2-1,l2+1
        do 6363,k1=l1-1,l1+1
        do 6363,jj=1,icell(0,k1,k2,k3)
          jat=icell(jj,k1,k2,k3)
          if (jat.eq.iat) goto 6363
          xrel= rxyz(1,iat)-rxyz(1,jat)
          yrel= rxyz(2,iat)-rxyz(2,jat)
          zrel= rxyz(3,iat)-rxyz(3,jat)  
          rr2=xrel**2 + yrel**2 + zrel**2
          if ( rr2 .le. cut2 ) then
           indlst=min(indlst,myspace-1)
           lstb(indlst)=lay(jat)
!        write(6,*) 'iat,indlst,lay(jat)',iat,indlst,lay(jat)
           tt=sqrt(rr2)
           tti=1.d0/tt
           rel(1,indlst)=xrel*tti
           rel(2,indlst)=yrel*tti
           rel(3,indlst)=zrel*tti
           rel(4,indlst)=tt
           rel(5,indlst)=tti
           indlst= indlst+1
          endif
6363        continue

        return
        end

subroutine tersoff_subeniat_l(i,Nmol,Npmax,NNmax,Kinds,X,R,R1,R2,Cr,Ca,alr,ala,XYZRrefdf,UadUrdf,Urtot,lsta,lstb,nnbrx,rel,pi)
      implicit none
      integer, intent(in)                       :: Nmol, Npmax, NNmax, nnbrx
      integer, intent(in)                       :: lsta(2,Nmol), lstb(nnbrx*Nmol)
      real(8), intent(in)                       :: rel(5,nnbrx*Nmol)
      integer, intent(in),  dimension(1:Nmol)   :: Kinds
      real(8), intent(in),  dimension(1:3*Nmol) :: R
      real(8), intent(in),  dimension(1:2,1:2)  :: R1, R2 
      real(8), intent(in),  dimension(1:2,1:2)  :: Cr, Ca, alr, ala
      real(8), intent(in),  dimension(1:2,1:2)  :: X
      real(8), intent(inout), dimension(1:6*Npmax):: XYZRrefdf
      real(8), intent(inout), dimension(1:3*Npmax):: UadUrdf
      real(8), intent(inout):: Urtot
      real(8) Ur, Ua
      integer i, j, l
      integer Ipt3,Jpt3,Nppt3,Nppt6
      integer Ki,Kj
      integer Nptot, NNtot, nat
      real(8) pi
      real(8) Xi, Yi, Zi
      real(8) Xij, Yij, Zij, Rij, Rreij
      real(8) R1ij, R2ij
      real(8) PL1, PL2
      real(8) Crij, Caij, alrij, alaij
      real(8) fij, dfij
      real(8) xhalf, yhalf, zhalf
      real(8) XRreik, YRreik,  ZRreik, Rreik, fik, dfik

!     #######################################
!     # Calculate XYZRrefdf, UadUrdf, Urtot #
!     #######################################
         nat=Nmol
         Ki=Kinds(i)
         Ipt3=3*(i-1)
         Xi=R(Ipt3+1)
         Yi=R(Ipt3+2)
         Zi=R(Ipt3+3)

      DO_J:  do l=lsta(1,i),lsta(2,i)
         j=lstb(l)

         Kj=Kinds(j)
         R2ij=R2(Ki,Kj)
         Jpt3=3*(j-1)
         Rij=rel(4,l)
         Xij=rel(1,l)
         Yij=rel(2,l)
         Zij=rel(3,l)
         Nptot=l
 
         Nppt3=3*(Nptot-1)
         Nppt6=6*(Nptot-1)
         Rreij=rel(5,l)

         XYZRrefdf(Nppt6+1)=Xij
         XYZRrefdf(Nppt6+2)=Yij
         XYZRrefdf(Nppt6+3)=Zij
         XYZRrefdf(Nppt6+4)=Rreij

         alrij = alr(Ki,Kj)
         alaij = ala(Ki,Kj)

         Ur= Cr(Ki,Kj)*dexp(-alrij*Rij)
         Ua=-Ca(Ki,Kj)*dexp(-alaij*Rij) * X(Ki,Kj)
         R1ij=R1(Ki,Kj)

      IF (Rij <= R1ij) THEN
          XYZRrefdf(Nppt6+5)=1.0d0
          XYZRrefdf(Nppt6+6)=0.0d0
          Urtot=Urtot+Ur
          UadUrdf(Nppt3+1)=Ua
          UadUrdf(Nppt3+2)=-alrij*Ur
          UadUrdf(Nppt3+3)=-alaij*Ua
      ELSE
          PL1=pi/(R2ij-R1ij)
          PL2=PL1*(Rij-R1ij)
          fij =0.5d0+0.5d0    *dcos(PL2)
          dfij=     -0.5d0*PL1*dsin(PL2)
          XYZRrefdf(Nppt6+5)=fij
          XYZRrefdf(Nppt6+6)=dfij
          Urtot=Urtot+fij*Ur
          UadUrdf(Nppt3+1)=fij*Ua
          UadUrdf(Nppt3+2)=(dfij-alrij*fij)*Ur
          UadUrdf(Nppt3+3)=(dfij-alaij*fij)*Ua
      END IF
      END DO DO_J
end subroutine

subroutine tersoff_subfiat_l(i,Nmol,Npmax,NNmax,Kinds,Pn,Co_bcd,bcsq,dsq,h,XYZRrefdf,UadUrdf,F,Uatot,dkEij,lsta,lstb,nnbrx)
      implicit none
      integer,intent(in)                     ::Nmol, Npmax, NNmax
      integer,intent(in),dimension(1:Nmol)   ::Kinds
      real(8),intent(in),dimension(1:2)      ::Pn, Co_bcd, bcsq, dsq, h
      real(8),intent(in),dimension(1:6*Npmax)::XYZRrefdf
      real(8),intent(in),dimension(1:3*Npmax)::UadUrdf
      real(8),intent(inout),dimension(1:3*Nmol)::F
      real(8),intent(inout)                    ::Uatot
      real(8),dimension(1:3*NNmax)           ::dkEij
      integer i
      integer Ki
      integer :: ij, ijpt3, ijpt6, ik, ikpt6
      integer Ipb, Ipe
      integer Ipt3, Jpt3, Kpt3
      integer Nkpt3

      real(8) XRreij, YRreij,  ZRreij, Rreij, fij, dfij 
      real(8) XRreik, YRreik,  ZRreik, Rreik, fik, dfik
      real(8) COSijk, Gi, dGi, fdG, fdGcos
      real(8) Co_pa, Co_mb1, Co_mb2
      real(8) Co_cdi, Co_hcosi, Co_dhcosi
      real(8) Pni, bcsqi, dsqi, hi
      real(8) Bij, Eij
      real(8) djEij
      real(8) dXjEij2, dYjEij2, dZjEij2
      real(8) Co1_dkEij, Co2_dkEij
      real(8) dFxi, dFyi, dFzi
      real(8) dFxj, dFyj, dFzj
      real(8) dFxk, dFyk, dFzk
      real(8) Ua
      integer, intent(in)                      :: nnbrx
      integer, intent(in)                      :: lsta(2,Nmol), lstb(nnbrx*Nmol)

         Ipb = lsta(1,i) 
         Ipe = lsta(2,i)

         Ki=Kinds(i)
         bcsqi=bcsq(Ki)
         dsqi=dsq(Ki)
         hi  =h(Ki)
         Pni =Pn(Ki)

         Co_cdi=Co_bcd(Ki) 

         dFxi=0.0d0
         dFyi=0.0d0
         dFzi=0.0d0

      DO_J:DO ij=Ipb,Ipe,+1

         IJpt3=3*(ij-1)
         IJpt6=6*(ij-1)

         XRreij = XYZRrefdf(IJpt6+1)
         YRreij = XYZRrefdf(IJpt6+2)
         ZRreij = XYZRrefdf(IJpt6+3)
         Rreij  = XYZRrefdf(IJpt6+4)
         fij    = XYZRrefdf(IJpt6+5)
         dfij   = XYZRrefdf(IJpt6+6)

         Eij     = 0.0d0
         djEij   = 0.0d0
         dXjEij2 = 0.0d0
         dYjEij2 = 0.0d0
         dZjEij2 = 0.0d0

         Nkpt3=-3
      DO_K:DO ik=Ipb,Ipe,+1

         Nkpt3=Nkpt3+3

      IKIJ:IF (ik /= ij ) THEN

         IKpt6=6*(ik-1)

         XRreik = XYZRrefdf(IKpt6+1)
         YRreik = XYZRrefdf(IKpt6+2)
         ZRreik = XYZRrefdf(IKpt6+3)
         Rreik  = XYZRrefdf(IKpt6+4)
         fik    = XYZRrefdf(IKpt6+5)
         dfik   = XYZRrefdf(IKpt6+6)

         COSijk = XRreij*XRreik + YRreij*YRreik + ZRreij*ZRreik

         Co_hcosi  = hi - COSijk
         Co_dhcosi = 1.0d0 / ( dsqi + Co_hcosi*Co_hcosi )
         Gi        = -bcsqi*Co_dhcosi
         dGi       = 2.0d0*Co_hcosi*Co_dhcosi*Gi
         Gi        = Gi + Co_cdi

         Eij = Eij + fik*Gi

         fdG    = fik * dGi
         fdGcos = fdG * COSijk

         djEij    = djEij   + fdGcos

         dXjEij2 = dXjEij2 + fdG*XRreik
         dYjEij2 = dYjEij2 + fdG*YRreik
         dZjEij2 = dZjEij2 + fdG*ZRreik


         Co1_dkEij = -dfik*Gi + fdGcos*Rreik
         Co2_dkEij = -fdG*Rreik

         dkEij(Nkpt3+1) = Co1_dkEij * XRreik + Co2_dkEij * XRreij
         dkEij(Nkpt3+2) = Co1_dkEij * YRreik + Co2_dkEij * YRreij
         dkEij(Nkpt3+3) = Co1_dkEij * ZRreik + Co2_dkEij * ZRreij

      ELSE 
         dkEij(Nkpt3+1) = 0.0d0
         dkEij(Nkpt3+2) = 0.0d0
         dkEij(Nkpt3+3) = 0.0d0
      END IF IKIJ

      END DO DO_K


      Bij  = 1.0d0 + Eij**Pni
      Ua   = UadUrdf(IJpt3+1) * Bij**(-0.5d0/Pni)
      Uatot = Uatot + Ua

      Co_pa  = UadUrdf(IJpt3+2) + UadUrdf(IJpt3+3)*Bij**(-0.5d0/Pni)


      CEij:IF ( Nkpt3 > 0 ) THEN 

         Co_mb1 = Ua * 0.5d0*Eij**(Pni-1.0d0) / Bij
         Co_mb2 = Co_mb1*Rreij

            Nkpt3=-3
         DO ik=Ipb,Ipe,+1

            Nkpt3=Nkpt3+3
            dFxk = Co_mb1 * dkEij(Nkpt3+1)
            dFyk = Co_mb1 * dkEij(Nkpt3+2)
            dFzk = Co_mb1 * dkEij(Nkpt3+3)

            Kpt3=3*(lstb(ik)-1)
            F(Kpt3+1) = F(Kpt3+1) + dFxk
            F(Kpt3+2) = F(Kpt3+2) + dFyk
            F(Kpt3+3) = F(Kpt3+3) + dFzk

            dFxi = dFxi + dFxk
            dFyi = dFyi + dFyk
            dFzi = dFzi + dFzk

         END DO

         dFxj = Co_pa*XRreij + Co_mb2*( XRreij*djEij - dXjEij2 )
         dFyj = Co_pa*YRreij + Co_mb2*( YRreij*djEij - dYjEij2 )
         dFzj = Co_pa*ZRreij + Co_mb2*( ZRreij*djEij - dZjEij2 )

      ELSE 

         dFxj = Co_pa*XRreij 
         dFyj = Co_pa*YRreij 
         dFzj = Co_pa*ZRreij 

      END IF CEij



      Jpt3=3*(lstb(ij)-1)

      F(Jpt3+1) = F(Jpt3+1) + dFxj
      F(Jpt3+2) = F(Jpt3+2) + dFyj
      F(Jpt3+3) = F(Jpt3+3) + dFzj

      dFxi = dFxi + dFxj
      dFyi = dFyi + dFyj
      dFzi = dFzi + dFzj

      END DO DO_J

         Ipt3 = 3*(i-1)
         F(Ipt3+1) = F(Ipt3+1) - dFxi
         F(Ipt3+2) = F(Ipt3+2) - dFyi
         F(Ipt3+3) = F(Ipt3+3) - dFzi
end subroutine

END MODULE eip_silicon
