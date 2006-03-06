!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2006  CP2K developers group                            !
!-----------------------------------------------------------------------------!

!!****h* cp2k/eip_silicon [1.0] *
!!
!!   NAME
!!     eip_silicon
!!
!!   FUNCTION
!!     Empirical interatomic potentials for Silicon
!!
!!   NOTES
!!     Stefan Goedecker's OpenMP implementation of Bazant's EDIP & Lenosky's
!!     empirical interatomic potentials for Silicon. 
!!     
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     03.2006 initial create [tdk]
!!
!!   SOURCE
!****************************************************************************
MODULE eip_silicon
  USE atomic_kind_list_types,          ONLY: atomic_kind_list_type
  USE atomic_kind_types,               ONLY: atomic_kind_type, &
                                             get_atomic_kind
  USE cell_types,                      ONLY: get_cell, &
                                             cell_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output, &
                                             cp_print_key_should_output, &
                                             cp_print_key_unit_nr, &
                                             cp_p_file
  USE cp_subsystem_types,              ONLY: cp_subsys_get, &
                                             cp_subsystem_type
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE eip_environment_types,           ONLY: eip_environment_type, &
                                             eip_env_get
  USE force_env_types,                 ONLY: force_env_type, &
                                             force_env_get
  USE input_section_types,             ONLY: section_vals_get_subs_vals, &
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: angstrom, &
                                             evolt
  USE timings,                         ONLY: timeset, &
                                             timestop

#include "cp_common_uses.h"

  !IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'eip_silicon'

  ! *** Public subroutines ***
  PUBLIC :: eip_bazant, eip_lenosky

!***
!****************************************************************************

CONTAINS

!!****f* eip_silicon/eip_bazant [1.0] *
!!
!!   NAME
!!     eip_bazant
!!
!!   FUNCTION
!!     Interface routine of Goedecker's Bazant EDIP to CP2K
!!
!!   REFERENCES
!!     http://www-math.mit.edu/~bazant/EDIP
!!     M.Z. Bazant & E. Kaxiras: Modeling of Covalent Bonding in Solids by 
!!                               Inversion of Cohesive Energy Curves;
!!                               Phys. Rev. Lett. 77, 4370 (1996)
!!     M.Z. Bazant, E. Kaxiras and J.F. Justo: Environment-dependent interatomic 
!!                                             potential for bulk silicon;
!!                                             Phys. Rev. B 56, 8542-8552 (1997)
!!     S. Goedecker: Optimization and parallelization of a force field for silicon 
!!                   using OpenMP; CPC 148, 1 (2002)
!!
!!   INPUTS
!!     - force_env: The force environment
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     03.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE eip_bazant(force_env, error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'eip_bazant', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    INTEGER                                  :: handle, &
                                                natom, &
                                                stat, &
                                                i, iw, &
                                                iparticle_kind, &
                                                nparticle_kind, &
                                                iparticle_local, &
                                                nparticle_local, &
                                                iparticle
    TYPE(eip_environment_type), POINTER      :: eip_env
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), DIMENSION(3)              :: abc
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    REAL(KIND=dp), DIMENSION(:,:), &
      ALLOCATABLE                            :: rxyz
    REAL(KIND=dp)                            :: ener, &
                                                ener_var, &
                                                ekin, &
                                                mass
    TYPE(section_vals_type), POINTER         :: eip_section
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(cp_subsystem_type), POINTER         :: subsys
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set

!   ------------------------------------------------------------------------

    CALL timeset(routineN, "I", "", handle)

    NULLIFY(eip_env, cell, particle_set, eip_section, logger, atomic_kinds, &
            atomic_kind, local_particles, subsys, atomic_kind_set)

    ekin = 0.0_dp

    failure = .FALSE.
    logger => cp_error_get_logger(error)

    CPPrecondition(ASSOCIATED(force_env), cp_failure_level, routineP, error, failure)
 
    IF (.NOT. failure) THEN

      CALL force_env_get(force_env=force_env, eip_env=eip_env, error=error)
      CALL eip_env_get(eip_env=eip_env, cell=cell, particle_set=particle_set, &
                       subsystem=subsys, local_particles=local_particles, &
                       atomic_kind_set=atomic_kind_set, error=error)
      CALL get_cell(cell=cell, abc=abc)

      eip_section => section_vals_get_subs_vals(eip_env%force_env_input, "EIP", &
                                                      error=error)
      natom = SIZE(particle_set)
      !natom = local_particles%n_el(1)

      ALLOCATE(rxyz(3,natom), stat=stat)
      CPPostcondition(stat == 0, cp_failure_level, routineP, error, failure)

      DO i = 1,natom
        !iparticle = local_particles%list(1)%array(i)
        rxyz(:,i) = particle_set(i)%r(:)*angstrom
      END DO

      CALL eip_bazant_silicon(nat=natom, alat=abc*angstrom, rxyz0=rxyz, &
                              fxyz=eip_env%eip_forces, ener=ener, &
                              coord=eip_env%coord_avg, ener_var=ener_var, &
                              coord_var=eip_env%coord_var, count=eip_env%count)

      !CALL get_part_ke(md_env, tbmd_energy%E_kinetic, int_grp=globalenv%para_env%group)
      CALL cp_subsys_get(subsys=subsys, atomic_kinds=atomic_kinds)

      nparticle_kind = atomic_kinds%n_els

      DO iparticle_kind = 1,nparticle_kind
        atomic_kind => atomic_kind_set(iparticle_kind)
        CALL get_atomic_kind(atomic_kind=atomic_kind, mass=mass)
        nparticle_local = local_particles%n_el(iparticle_kind)
        DO iparticle_local = 1,nparticle_local
          iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
          ekin = ekin +  0.5_dp * mass * &
                     (particle_set(iparticle)%v(1)* particle_set(iparticle)%v(1) &
                    + particle_set(iparticle)%v(2)* particle_set(iparticle)%v(2) &
                    + particle_set(iparticle)%v(3)* particle_set(iparticle)%v(3))
        END DO
      END DO

      ! sum all contributions to energy over calculated parts on all processors
      CALL mp_sum(ekin, force_env%globenv%para_env%group)
      eip_env%eip_kinetic_energy = ekin

      eip_env%eip_potential_energy = ener/evolt
      eip_env%eip_energy = eip_env%eip_kinetic_energy + eip_env%eip_potential_energy
      eip_env%eip_energy_var = ener_var/evolt

      DO i = 1,natom
        particle_set(i)%f(:) = eip_env%eip_forces(:,i)/evolt*angstrom
      END DO

      DEALLOCATE(rxyz, stat=stat)
      CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)

      ! Print
      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                eip_section, "PRINT%ENERGIES", error=error), cp_p_file)) THEN
        iw = cp_print_key_unit_nr(logger, eip_section, "PRINT%ENERGIES", &
                                  extension=".mmLog", error=error)

        CALL eip_print_energies(eip_env=eip_env, output_unit=iw, error=error)
        CALL cp_print_key_finished_output(iw, logger, eip_section, &
                                          "PRINT%ENERGIES", error=error)
      END IF

      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                eip_section, "PRINT%ENERGIES_VAR", error=error), cp_p_file)) THEN
        iw = cp_print_key_unit_nr(logger, eip_section, "PRINT%ENERGIES_VAR", &
                                  extension=".mmLog", error=error)

        CALL eip_print_energy_var(eip_env=eip_env, output_unit=iw, error=error)
        CALL cp_print_key_finished_output(iw, logger, eip_section, &
                                          "PRINT%ENERGIES_VAR", error=error)
      END IF

      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                eip_section, "PRINT%FORCES", error=error), cp_p_file)) THEN
        iw = cp_print_key_unit_nr(logger, eip_section, "PRINT%FORCES", &
                                  extension=".mmLog", error=error)

        CALL eip_print_forces(eip_env=eip_env, output_unit=iw, error=error)
        CALL cp_print_key_finished_output(iw, logger, eip_section, &
                                          "PRINT%FORCES", error=error)
      END IF

      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                eip_section, "PRINT%COORD_AVG", error=error), cp_p_file)) THEN
        iw = cp_print_key_unit_nr(logger, eip_section, "PRINT%COORD_AVG", &
                                  extension=".mmLog", error=error)

        CALL eip_print_coord_avg(eip_env=eip_env, output_unit=iw, error=error)
        CALL cp_print_key_finished_output(iw, logger, eip_section, &
                                          "PRINT%COORD_AVG", error=error)
      END IF

      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                eip_section, "PRINT%COORD_VAR", error=error), cp_p_file)) THEN
        iw = cp_print_key_unit_nr(logger, eip_section, "PRINT%COORD_VAR", &
                                  extension=".mmLog", error=error)

        CALL eip_print_coord_var(eip_env=eip_env, output_unit=iw, error=error)
        CALL cp_print_key_finished_output(iw, logger, eip_section, &
                                          "PRINT%COORD_VAR", error=error)
      END IF

      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                eip_section, "PRINT%COUNT", error=error), cp_p_file)) THEN
        iw = cp_print_key_unit_nr(logger, eip_section, "PRINT%COUNT", &
                                  extension=".mmLog", error=error)

        CALL eip_print_count(eip_env=eip_env, output_unit=iw, error=error)
        CALL cp_print_key_finished_output(iw, logger, eip_section, &
                                          "PRINT%COUNT", error=error)
      END IF

    END IF

    CALL timestop(0.0_dp, handle)

  END SUBROUTINE eip_bazant
!****************************************************************************

!!****f* eip_silicon/eip_lenosky [1.0] *
!!
!!   NAME
!!     eip_lenosky
!!
!!   FUNCTION
!!     Interface routine of Goedecker's Lenosky force field to CP2K
!!
!!   REFERENCES
!!     T. Lenosky, et. al.: Highly optimized empirical potential model of silicon;
!!                          Modelling Simul. Sci. Eng., 8 (2000)
!!     S. Goedecker: Optimization and parallelization of a force field for silicon 
!!                   using OpenMP; CPC 148, 1 (2002)
!!
!!   INPUTS
!!     - force_env: The force environment
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     03.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE eip_lenosky(force_env, error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'eip_lenosky', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    INTEGER                                  :: handle, &
                                                natom, &
                                                stat, &
                                                i, iw, &
                                                iparticle_kind, &
                                                nparticle_kind, &
                                                iparticle_local, &
                                                nparticle_local, &
                                                iparticle
    TYPE(eip_environment_type), POINTER      :: eip_env
    TYPE(cell_type), POINTER                 :: cell
    REAL(KIND=dp), DIMENSION(3)              :: abc
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    REAL(KIND=dp), DIMENSION(:,:), &
      ALLOCATABLE                            :: rxyz
    REAL(KIND=dp)                            :: ener, &
                                                ener_var, &
                                                ekin, &
                                                mass
    TYPE(section_vals_type), POINTER         :: eip_section
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(atomic_kind_list_type), POINTER     :: atomic_kinds
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(cp_subsystem_type), POINTER         :: subsys
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set

!   ------------------------------------------------------------------------

    CALL timeset(routineN, "I", "", handle)

    NULLIFY(eip_env, cell, particle_set, eip_section, logger, atomic_kinds, &
            atomic_kind, local_particles, subsys, atomic_kind_set)

    ekin = 0.0_dp

    failure = .FALSE.
    logger => cp_error_get_logger(error)  
  
    CPPrecondition(ASSOCIATED(force_env), cp_failure_level, routineP, error, failure)
 
    IF (.NOT. failure) THEN

      CALL force_env_get(force_env=force_env, eip_env=eip_env, error=error)
      CALL eip_env_get(eip_env=eip_env, cell=cell, particle_set=particle_set, &
                       subsystem=subsys, local_particles=local_particles, &
                       atomic_kind_set=atomic_kind_set, error=error)
      CALL get_cell(cell=cell, abc=abc)

      eip_section => section_vals_get_subs_vals(eip_env%force_env_input, "EIP", &
                                                      error=error)
      natom = SIZE(particle_set)
      !natom = local_particles%n_el(1)

      ALLOCATE(rxyz(3,natom), stat=stat)
      CPPostcondition(stat == 0, cp_failure_level, routineP, error, failure)

      DO i = 1,natom
        !iparticle = local_particles%list(1)%array(i)
        rxyz(:,i) = particle_set(i)%r(:)*angstrom
      END DO

      CALL eip_lenosky_silicon(nat=natom, alat=abc*angstrom, rxyz0=rxyz, &
                               fxyz=eip_env%eip_forces, ener=ener, &
                               coord=eip_env%coord_avg, ener_var=ener_var, &
                               coord_var=eip_env%coord_var, count=eip_env%count)

      !CALL get_part_ke(md_env, tbmd_energy%E_kinetic, int_grp=globalenv%para_env%group)
      CALL cp_subsys_get(subsys=subsys, atomic_kinds=atomic_kinds)

      nparticle_kind = atomic_kinds%n_els

      DO iparticle_kind = 1,nparticle_kind
        atomic_kind => atomic_kind_set(iparticle_kind)
        CALL get_atomic_kind(atomic_kind=atomic_kind, mass=mass)
        nparticle_local = local_particles%n_el(iparticle_kind)
        DO iparticle_local = 1,nparticle_local
          iparticle = local_particles%list(iparticle_kind)%array(iparticle_local)
          ekin = ekin +  0.5_dp * mass * &
                     (particle_set(iparticle)%v(1)* particle_set(iparticle)%v(1) &
                    + particle_set(iparticle)%v(2)* particle_set(iparticle)%v(2) &
                    + particle_set(iparticle)%v(3)* particle_set(iparticle)%v(3))
        END DO
      END DO

      ! sum all contributions to energy over calculated parts on all processors
      CALL mp_sum(ekin, force_env%globenv%para_env%group)
      eip_env%eip_kinetic_energy = ekin

      eip_env%eip_potential_energy = ener/evolt
      eip_env%eip_energy = eip_env%eip_kinetic_energy + eip_env%eip_potential_energy
      eip_env%eip_energy_var = ener_var/evolt

      DO i = 1,natom
        particle_set(i)%f(:) = eip_env%eip_forces(:,i)/evolt*angstrom
      END DO

      DEALLOCATE(rxyz, stat=stat)
      CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)

      ! Print
      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                eip_section, "PRINT%ENERGIES", error=error), cp_p_file)) THEN
        iw = cp_print_key_unit_nr(logger, eip_section, "PRINT%ENERGIES", &
                                  extension=".mmLog", error=error)

        CALL eip_print_energies(eip_env=eip_env, output_unit=iw, error=error)
        CALL cp_print_key_finished_output(iw, logger, eip_section, &
                                          "PRINT%ENERGIES", error=error)
      END IF

      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                eip_section, "PRINT%ENERGIES_VAR", error=error), cp_p_file)) THEN
        iw = cp_print_key_unit_nr(logger, eip_section, "PRINT%ENERGIES_VAR", &
                                  extension=".mmLog", error=error)

        CALL eip_print_energy_var(eip_env=eip_env, output_unit=iw, error=error)
        CALL cp_print_key_finished_output(iw, logger, eip_section, &
                                          "PRINT%ENERGIES_VAR", error=error)
      END IF

      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                eip_section, "PRINT%FORCES", error=error), cp_p_file)) THEN
        iw = cp_print_key_unit_nr(logger, eip_section, "PRINT%FORCES", &
                                  extension=".mmLog", error=error)

        CALL eip_print_forces(eip_env=eip_env, output_unit=iw, error=error)
        CALL cp_print_key_finished_output(iw, logger, eip_section, &
                                          "PRINT%FORCES", error=error)
      END IF

      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                eip_section, "PRINT%COORD_AVG", error=error), cp_p_file)) THEN
        iw = cp_print_key_unit_nr(logger, eip_section, "PRINT%COORD_AVG", &
                                  extension=".mmLog", error=error)

        CALL eip_print_coord_avg(eip_env=eip_env, output_unit=iw, error=error)
        CALL cp_print_key_finished_output(iw, logger, eip_section, &
                                          "PRINT%COORD_AVG", error=error)
      END IF

      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                eip_section, "PRINT%COORD_VAR", error=error), cp_p_file)) THEN
        iw = cp_print_key_unit_nr(logger, eip_section, "PRINT%COORD_VAR", &
                                  extension=".mmLog", error=error)

        CALL eip_print_coord_var(eip_env=eip_env, output_unit=iw, error=error)
        CALL cp_print_key_finished_output(iw, logger, eip_section, &
                                          "PRINT%COORD_VAR", error=error)
      END IF

      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                eip_section, "PRINT%COUNT", error=error), cp_p_file)) THEN
        iw = cp_print_key_unit_nr(logger, eip_section, "PRINT%COUNT", &
                                  extension=".mmLog", error=error)

        CALL eip_print_count(eip_env=eip_env, output_unit=iw, error=error)
        CALL cp_print_key_finished_output(iw, logger, eip_section, &
                                          "PRINT%COUNT", error=error)
      END IF

    END IF

    CALL timestop(0.0_dp, handle)

  END SUBROUTINE eip_lenosky
!****************************************************************************

!!****f* eip_silicon/eip_print_energies [1.0] *
!!
!!   NAME
!!     eip_print_energies
!!
!!   FUNCTION
!!     Print routine for the EIP energies
!!
!!   NOTES
!!     As usual the EIP energies differ from the DFT energies!
!!     Only the relative energy differneces are correctly reproduced.
!!
!!   INPUTS
!!     - eip_env: The eip environment of matter
!!     - output_unit: The output unit
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     03.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE eip_print_energies(eip_env, output_unit, error)
    TYPE(eip_environment_type), POINTER      :: eip_env
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'eip_print_energies', &
      routineP = moduleN//':'//routineN

!   ------------------------------------------------------------------------

    IF (output_unit > 0) THEN
      WRITE (UNIT=output_unit,FMT="(/,(T3,A,T55,F25.14))")&
             "Kinetic energy [Hartree]:        ",eip_env%eip_kinetic_energy, &
             "Potential energy [Hartree]:      ",eip_env%eip_potential_energy, &
             "Total EIP energy [Hartree]:      ",eip_env%eip_energy
    END IF

  END SUBROUTINE eip_print_energies
!***************************************************************************

!!****f* eip_silicon/eip_print_energy_var [1.0] *
!!
!!   NAME
!!     eip_print_energy_var
!!
!!   FUNCTION
!!     Print routine for the variance of the energy/atom
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - eip_env: The eip environment of matter
!!     - output_unit: The output unit
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     03.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE eip_print_energy_var(eip_env, output_unit, error)
    TYPE(eip_environment_type), POINTER      :: eip_env
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'eip_print_energy_var', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: unit_nr

!   ------------------------------------------------------------------------

    unit_nr = output_unit

    IF (unit_nr > 0) THEN

      WRITE (unit_nr,*) ""
      WRITE (unit_nr,*) "The variance of the EIP energy/atom!"
      WRITE (unit_nr,*) ""
      WRITE (unit_nr,*) eip_env%eip_energy_var

    END IF

  END SUBROUTINE eip_print_energy_var
!***************************************************************************

!!****f* eip_silicon/eip_print_forces [1.0] *
!!
!!   NAME
!!     eip_print_forces
!!
!!   FUNCTION
!!     Print routine for the forces
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - eip_env: The eip environment of matter
!!     - output_unit: The output unit
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     03.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE eip_print_forces(eip_env, output_unit, error)
    TYPE(eip_environment_type), POINTER      :: eip_env
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'eip_print_forces', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: iatom, natom, unit_nr
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set

!   ------------------------------------------------------------------------

    NULLIFY(particle_set)

    unit_nr = output_unit

    IF (unit_nr > 0) THEN

      CALL eip_env_get(eip_env=eip_env, particle_set=particle_set, error=error)

      natom = SIZE(particle_set)

      WRITE (unit_nr,*) ""
      WRITE (unit_nr,*) "The EIP forces!"
      WRITE (unit_nr,*) ""
      WRITE(unit_nr, *) "Total EIP forces [Hartree/Bohr]"
      DO iatom = 1, natom
        WRITE (unit_nr,*) eip_env%eip_forces(1:3,iatom)
      END DO

    END IF

  END SUBROUTINE eip_print_forces
!***************************************************************************

!!****f* eip_silicon/eip_print_coord_avg [1.0] *
!!
!!   NAME
!!     eip_print_coord_avg
!!
!!   FUNCTION
!!     Print routine for the average coordination number
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - eip_env: The eip environment of matter
!!     - output_unit: The output unit
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     03.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE eip_print_coord_avg(eip_env, output_unit, error)
    TYPE(eip_environment_type), POINTER      :: eip_env
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'eip_print_coord_avg', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: unit_nr

!   ------------------------------------------------------------------------

    unit_nr = output_unit

    IF (unit_nr > 0) THEN

      WRITE (unit_nr,*) ""
      WRITE (unit_nr,*) "The average coordination number!"
      WRITE (unit_nr,*) ""
      WRITE (unit_nr,*) eip_env%coord_avg

    END IF

  END SUBROUTINE eip_print_coord_avg
!***************************************************************************

!!****f* eip_silicon/eip_print_coord_var [1.0] *
!!
!!   NAME
!!     eip_print_coord_var
!!
!!   FUNCTION
!!     Print routine for the variance of the coordination number
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - eip_env: The eip environment of matter
!!     - output_unit: The output unit
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     03.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE eip_print_coord_var(eip_env, output_unit, error)
    TYPE(eip_environment_type), POINTER      :: eip_env
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'eip_print_coord_var', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: unit_nr

!   ------------------------------------------------------------------------

    unit_nr = output_unit

    IF (unit_nr > 0) THEN

      WRITE (unit_nr,*) ""
      WRITE (unit_nr,*) "The variance of the coordination number!"
      WRITE (unit_nr,*) ""
      WRITE (unit_nr,*) eip_env%coord_var

    END IF

  END SUBROUTINE eip_print_coord_var
!***************************************************************************

!!****f* eip_silicon/eip_print_count [1.0] *
!!
!!   NAME
!!     eip_print_count
!!
!!   FUNCTION
!!     Print routine for the function call counter
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - eip_env: The eip environment of matter
!!     - output_unit: The output unit
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     03.2006 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE eip_print_count(eip_env, output_unit, error)
    TYPE(eip_environment_type), POINTER      :: eip_env
    INTEGER, INTENT(IN)                      :: output_unit
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'eip_print_count', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: unit_nr

!   ------------------------------------------------------------------------

    unit_nr = output_unit

    IF (unit_nr > 0) THEN

      WRITE (unit_nr,*) ""
      WRITE (unit_nr,*) "The function call counter!"
      WRITE (unit_nr,*) ""
      WRITE (unit_nr,*) eip_env%count

    END IF

  END SUBROUTINE eip_print_count
!***************************************************************************

!!****f* eip_silicon/eip_bazant_silicon [1.0] *
!!
!!   NAME
!!     eip_bazant_silicon
!!
!!   FUNCTION
!!     Bazant's EDIP (environment-dependent interatomic potential) for Silicon 
!!     by Stefan Goedecker
!!
!!   REFERENCES
!!     http://www-math.mit.edu/~bazant/EDIP
!!     M.Z. Bazant & E. Kaxiras: Modeling of Covalent Bonding in Solids by 
!!                               Inversion of Cohesive Energy Curves;
!!                               Phys. Rev. Lett. 77, 4370 (1996)
!!     M.Z. Bazant, E. Kaxiras and J.F. Justo: Environment-dependent interatomic 
!!                                             potential for bulk silicon;
!!                                             Phys. Rev. B 56, 8542-8552 (1997)
!!     S. Goedecker: Optimization and parallelization of a force field for silicon 
!!                   using OpenMP; CPC 148, 1 (2002)
!!
!!   INPUTS
!!     - nat: number of atoms
!!     - alat: lattice constants of the orthorombic box containing the particles
!!     - rxyz0: atomic positions in Angstrom
!!
!!   OUTPUTS
!!     - rxyz0: rxyz0 may be modified on output.
!!              If an atom is outside the box the program will bring it back 
!!              into the box by translations through alat
!!     - fxyz: forces in eV/A
!!     - ener: total energy in eV
!!     - coord: average coordination number
!!     - ener_var: variance of the energy/atom
!!     - coord_var: variance of the coordination number
!!     - count: count is increased by one per call, has to be initialized 
!               to 0.d0 before first call of eip_bazant
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     03.2006 initial create [tdk]
!!
!!*** **********************************************************************
    subroutine eip_bazant_silicon(nat,alat,rxyz0,fxyz,ener,coord,ener_var,&
                                  coord_var,count)

        implicit REAL(KIND=dp) (a-h,o-z)
!$      interface
!$        integer ( kind=4 ) function omp_get_num_threads ( )
!$        end function omp_get_num_threads
!$      end interface
!$      interface
!$        integer ( kind=4 ) function omp_get_thread_num ( )
!$        end function omp_get_thread_num
!$      end interface

        dimension rxyz0(3,nat),fxyz(3,nat),alat(3)
        REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:) :: rxyz
        integer, ALLOCATABLE, DIMENSION(:,:)       :: lsta
        integer, ALLOCATABLE, DIMENSION(:)         :: lstb
        integer, ALLOCATABLE, DIMENSION(:)         :: lay
        integer, ALLOCATABLE, DIMENSION(:,:,:,:)   :: icell
        REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:) :: rel
        REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:) :: txyz
        REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:) :: s2,s3,sz
        integer, ALLOCATABLE, DIMENSION(:)         :: num2,num3,numz


!        cut=par_a
        cut= 3.1213820d0 + 1.d-14

        if (count.eq.0)  open(unit=10,file='bazant.mon',status='unknown')
        count=count+1.d0

! linear scaling calculation of verlet list
        ll1=int(alat(1)/cut)
        if (ll1.lt.1) stop 'alat(1) too small'
        ll2=int(alat(2)/cut)
        if (ll2.lt.1) stop 'alat(2) too small'
        ll3=int(alat(3)/cut)
        if (ll3.lt.1) stop 'alat(3) too small'

! determine number of threads
        npr=1
!$omp parallel private(iam)  shared (npr)
!$       iam=omp_get_thread_num()
!$       if (iam.eq.0) npr=omp_get_num_threads()
!$omp end parallel 

! linear scaling calculation of verlet list

     if (npr.le.1) then !serial if too few processors to gain by parallelizing

! set ncx for serial case, ncx for parallel case set below
        ncx=8
1234    ncx=ncx*2
        allocate(icell(0:ncx,-1:ll1,-1:ll2,-1:ll3))
        do l3=-1,ll3
        do l2=-1,ll2
        do l1=-1,ll1
        icell(0,l1,l2,l3)=0
        enddo
        enddo
        enddo
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

     else  ! parallel case

! periodization of particles can be done in parallel
!$omp parallel do shared (alat,nat,rxyz0) private(iat)
        do 5983,iat=1,nat
        rxyz0(1,iat)=modulo(modulo(rxyz0(1,iat),alat(1)),alat(1))
        rxyz0(2,iat)=modulo(modulo(rxyz0(2,iat),alat(2)),alat(2))
        rxyz0(3,iat)=modulo(modulo(rxyz0(3,iat),alat(3)),alat(3))
5983    continue
!$omp end parallel do

! assignment to cell is done serially
! set ncx for parallel case, ncx for serial case set above
        ncx=8
4321    ncx=ncx*2
        allocate(icell(0:ncx,-1:ll1,-1:ll2,-1:ll3))
        do l3=-1,ll3
        do l2=-1,ll2
        do l1=-1,ll1
        icell(0,l1,l2,l3)=0
        enddo
        enddo
        enddo

        rlc1i=ll1/alat(1)
        rlc2i=ll2/alat(2)
        rlc3i=ll3/alat(3)

        do 6983,iat=1,nat
        l1=int(rxyz0(1,iat)*rlc1i)
        l2=int(rxyz0(2,iat)*rlc2i)
        l3=int(rxyz0(3,iat)*rlc3i)
        ii=icell(0,l1,l2,l3)
        ii=ii+1
        icell(0,l1,l2,l3)=ii
        if (ii.gt.ncx) then
        write(10,*) count,'NCX too small',ncx
        deallocate(icell)
        goto 4321
        endif
        icell(ii,l1,l2,l3)=iat
6983    continue

    endif


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
        nnbrx=8
2345    nnbrx=3*nnbrx/2
        allocate(lstb(nnbrx*nat),rel(5,nnbrx*nat))

        indlstx=0

!$omp parallel  &
!$omp private(iat,cut2,iam,ii,indlst,l1,l2,l3,myspace,npr) &
!$omp shared (indlstx,nat,nn,nnbrx,ncx,ll1,ll2,ll3,icell,lsta,lstb,lay, & 
!$omp rel,rxyz,cut,myspaceout)


        npr=1
!$       npr=omp_get_num_threads()
        iam=0
!$       iam=omp_get_thread_num()

        cut2=cut**2
! assign contiguous portions of the arrays lstb and rel to the threads
        myspace=(nat*nnbrx)/npr  
        if (iam.eq.0) myspaceout=myspace
! Verlet list, relative positions
        indlst=0
      do 6003,l3=0,ll3-1
      do 6002,l2=0,ll2-1
      do 6001,l1=0,ll1-1
      do 6000,ii=1,icell(0,l1,l2,l3)
        iat=icell(ii,l1,l2,l3)
        if ( ((iat-1)*npr)/nat .eq. iam) then
!       write(6,*) 'sublstiat:iam,iat',iam,iat
        lsta(1,iat)=iam*myspace+indlst+1
        call sublstiat_b(iat,nn,ncx,ll1,ll2,ll3,l1,l2,l3,myspace, & 
             rxyz,icell,lstb(iam*myspace+1),lay,rel(1,iam*myspace+1),cut2,indlst)
        lsta(2,iat)=iam*myspace+indlst
!        write(6,'(a,4(x,i3),100(x,i2))') & 
!               'iam,iat,lsta',iam,iat,lsta(1,iat),lsta(2,iat), & 
!                    (lstb(j),j=lsta(1,iat),lsta(2,iat))
        endif

6000    continue
6001    continue
6002    continue
6003    continue
!$omp critical
        indlstx=max(indlstx,indlst)
!$omp end critical
!$omp end parallel  

           if (indlstx.ge.myspaceout) then
               write(10,*) count,'NNBRX too  small', nnbrx
               deallocate(lstb,rel)
               goto 2345
           endif

!$omp parallel  &
!$omp private(iam,npr,iat,iat1,iat2,lot,istop,tcoord,tcoord2, & 
!$omp tener,tener2,txyz,s2,s3,sz,num2,num3,numz,max_nbrs) &
!$omp shared (nat,nnbrx,lsta,lstb,rel,ener,ener2,fxyz,coord,coord2,istopg)

        npr=1
!$       npr=omp_get_num_threads()
        iam=0
!$       iam=omp_get_thread_num()

         max_nbrs=30
         istopg=0

        if (npr.ne.1) then 
! PARALLEL CASE
! create temporary private scalars for reduction sum on energies and 
!        temporary private array for reduction sum on forces
!$omp critical
        allocate(txyz(3,nat),s2(max_nbrs,8),s3(max_nbrs,7),sz(max_nbrs,6),  & 
                 num2(max_nbrs),num3(max_nbrs),numz(max_nbrs))
!$omp end critical
        if (iam.eq.0) then
        ener=0.d0
        ener2=0.d0
        coord=0.d0
        coord2=0.d0
        endif
!$omp do
        do iat=1,nat
        fxyz(1,iat)=0.d0
        fxyz(2,iat)=0.d0
        fxyz(3,iat)=0.d0
        enddo
!$omp barrier

! Each thread treats at most lot atoms
        lot=int(float(nat)/float(npr)+.999999999999d0)
        iat1=iam*lot+1
        iat2=min((iam+1)*lot,nat)
!       write(6,*) 'subfeniat:iat1,iat2,iam',iat1,iat2,iam
        call subfeniat_b(iat1,iat2,nat,lsta,lstb,rel,tener,tener2,  &
          tcoord,tcoord2,nnbrx,txyz,max_nbrs,istop,  &
          s2(1,1),s2(1,2),s2(1,3),s2(1,4),s2(1,5),s2(1,6),s2(1,7),s2(1,8),  &
          num2,s3(1,1),s3(1,2),s3(1,3),s3(1,4),s3(1,5),s3(1,6),s3(1,7),  &
          num3,sz(1,1),sz(1,2),sz(1,3),sz(1,4),sz(1,5),sz(1,6),numz)

!$omp critical
        ener=ener+tener
        ener2=ener2+tener2
        coord=coord+tcoord
        coord2=coord2+tcoord2
        istopg=istopg+istop
        do 8093,iat=1,nat
        fxyz(1,iat)=fxyz(1,iat)+txyz(1,iat)
        fxyz(2,iat)=fxyz(2,iat)+txyz(2,iat)
        fxyz(3,iat)=fxyz(3,iat)+txyz(3,iat)
8093    continue
        deallocate(txyz,s2,s3,sz,num2,num3,numz)
!$omp end critical

        else
! SERIAL CASE
        iat1=1
        iat2=nat
        allocate(s2(max_nbrs,8),s3(max_nbrs,7),sz(max_nbrs,6),  & 
                 num2(max_nbrs),num3(max_nbrs),numz(max_nbrs))
        call subfeniat_b(iat1,iat2,nat,lsta,lstb,rel,ener,ener2,  &
          coord,coord2,nnbrx,fxyz,max_nbrs,istopg,  &
          s2(1,1),s2(1,2),s2(1,3),s2(1,4),s2(1,5),s2(1,6),s2(1,7),s2(1,8),  &
          num2,s3(1,1),s3(1,2),s3(1,3),s3(1,4),s3(1,5),s3(1,6),s3(1,7),  &
          num3,sz(1,1),sz(1,2),sz(1,3),sz(1,4),sz(1,5),sz(1,6),numz)
        deallocate(s2,s3,sz,num2,num3,numz)

        endif
!$omp end parallel

!         write(6,*) 'ener,norm force', & 
!                    ener,DNRM2(3*nat,fxyz,1)
        if (istopg.gt.0) stop 'DIMENSION ERROR (see WARNING above)' 
        ener_var=ener2/nat-(ener/nat)**2 
        coord=coord/nat
        coord_var=coord2/nat-coord**2 

        deallocate(rxyz,icell,lay,lsta,lstb,rel)

        end subroutine eip_bazant_silicon



        subroutine subfeniat_b(iat1,iat2,nat,lsta,lstb,rel,ener,ener2,  &
          coord,coord2,nnbrx,ff,max_nbrs,istop,  &
          s2_t0,s2_t1,s2_t2,s2_t3,s2_dx,s2_dy,s2_dz,s2_r,  &
          num2,s3_g,s3_dg,s3_rinv,s3_dx,s3_dy,s3_dz,s3_r,  &
          num3,sz_df,sz_sum,sz_dx,sz_dy,sz_dz,sz_r,numz)
! This subroutine is a modification of a subroutine that is available at 
! http://www-math.mit.edu/~bazant/EDIP/ and for which Martin Z. Bazant 
! and Harvard University have a 1997 copyright.
! The modifications were done by S. Goedecker on April 10, 2002. 
! The routines are included with the permission of M. Bazant into this package.
        
        implicit none
!  ------------------------- VARIABLE DECLARATIONS -------------------------
          integer       iat1,iat2,nat
          REAL(KIND=dp) ener,ener2,coord,coord2
          REAL(KIND=dp) xarg,coord_iat,ener_iat
          REAL(KIND=dp) ff(3,nat)

          REAL(KIND=dp) par_cap_A,par_cap_B,par_rh,par_a,par_sig,par_lam,par_gam, &
                        par_b,par_c,par_delta,par_mu,par_Qo,par_palp, &
                        par_bet,par_alp,par_bg,par_eta,u1,u2,u3,u4,u5

          integer       nnbrx,max_nbrs,istop
          integer       lsta(2,nat),lstb(nnbrx*nat)
          REAL(KIND=dp) rel(5,nnbrx*nat)

          integer       i,j,k,l,n
          REAL(KIND=dp) dx,dy,dz,r
          REAL(KIND=dp) rinv,rmainv,xinv,xinv3,den,Z,fZ
          REAL(KIND=dp) dV2j,dV2ijx,dV2ijy,dV2ijz,pZ,dp1
          REAL(KIND=dp) temp0,temp1
          REAL(KIND=dp) Qort,muhalf
          REAL(KIND=dp) rmbinv,winv,dwinv,tau,dtau,lcos,x,H,dHdx,dhdl
          REAL(KIND=dp) dV3rij,dV3rijx,dV3rijy,dV3rijz
          REAL(KIND=dp) dV3rik,dV3rikx,dV3riky,dV3rikz
          REAL(KIND=dp) dV3l,dV3ljx,dV3ljy,dV3ljz,dV3lkx,dV3lky,dV3lkz
          REAL(KIND=dp) dV2dZ,dxdZ,dV3dZ
          REAL(KIND=dp) dEdrl,dEdrlx,dEdrly,dEdrlz
          REAL(KIND=dp) bmc,cmbinv
          REAL(KIND=dp) fjx,fjy,fjz,fkx,fky,fkz
        
          REAL(KIND=dp) s2_t0(max_nbrs)
          REAL(KIND=dp) s2_t1(max_nbrs)
          REAL(KIND=dp) s2_t2(max_nbrs)
          REAL(KIND=dp) s2_t3(max_nbrs)
          REAL(KIND=dp) s2_dx(max_nbrs)
          REAL(KIND=dp) s2_dy(max_nbrs)
          REAL(KIND=dp) s2_dz(max_nbrs)
          REAL(KIND=dp) s2_r(max_nbrs)
          integer       n2                
!   size of s2[]
          integer       num2(max_nbrs)  
!   atom ID numbers for s2[]
        
          REAL(KIND=dp) s3_g(max_nbrs)
          REAL(KIND=dp) s3_dg(max_nbrs)
          REAL(KIND=dp) s3_rinv(max_nbrs)
          REAL(KIND=dp) s3_dx(max_nbrs)
          REAL(KIND=dp) s3_dy(max_nbrs)
          REAL(KIND=dp) s3_dz(max_nbrs)
          REAL(KIND=dp) s3_r(max_nbrs)
        
          integer       n3                
!   size of s3[]
          integer       num3(max_nbrs)  
!   atom ID numbers for s3[]
        
          REAL(KIND=dp) sz_df(max_nbrs)
          REAL(KIND=dp) sz_sum(max_nbrs)
          REAL(KIND=dp) sz_dx(max_nbrs)
          REAL(KIND=dp) sz_dy(max_nbrs)
          REAL(KIND=dp) sz_dz(max_nbrs)
          REAL(KIND=dp) sz_r(max_nbrs)
          integer       nz
!   size of sz[]
          integer       numz(max_nbrs)  
!   atom ID numbers for sz[]
        
          integer       nj,nk,nl         
!   indices for the store arrays
        
!   EDIP parameters
          par_cap_A = 5.6714030d0
          par_cap_B = 2.0002804d0
          par_rh = 1.2085196d0
          par_a = 3.1213820d0
          par_sig = 0.5774108d0
          par_lam = 1.4533108d0
          par_gam = 1.1247945d0
          par_b = 3.1213820d0
          par_c = 2.5609104d0
          par_delta = 78.7590539d0
          par_mu = 0.6966326d0
          par_Qo = 312.1341346d0
          par_palp = 1.4074424d0
          par_bet = 0.0070975d0
          par_alp = 3.1083847d0

          u1 = -0.165799d0
          u2 = 32.557d0
          u3 = 0.286198d0
          u4 = 0.66d0

          par_bg=par_a
          par_eta = par_delta/par_Qo

          do i=1, nat
            ff(1,i) = 0.0d0
            ff(2,i) = 0.0d0
            ff(3,i) = 0.0d0
          end do
        
          coord=0.d0
          coord2=0.d0
          ener=0.d0
          ener2=0.d0
          istop=0
        
          
!   COMBINE COEFFICIENTS
        
          Qort = sqrt(par_Qo)
          muhalf = par_mu*0.5D0
          u5 = u2*u4
          bmc = par_b-par_c
          cmbinv = 1.0D0/(par_c-par_b)
        
        
          
!  --- LEVEL 1: OUTER LOOP OVER ATOMS ---
        
          do 1000, i= iat1, iat2
            
!   RESET COORDINATION AND NEIGHBOR NUMBERS
        
            coord_iat=0.d0
            ener_iat=0.d0
            Z = 0.0d0
            n2 = 1
            n3 = 1
            nz = 1
        
            
!  --- LEVEL 2: LOOP PREPASS OVER PAIRS ---
        
            do n=lsta(1,i),lsta(2,i)
              j=lstb(n)
        
                
!   PARTS OF TWO-BODY INTERACTION r<par_a
        
                num2(n2) = j
                dx = -rel(1,n)
                dy = -rel(2,n)
                dz = -rel(3,n)
                r=rel(4,n)
                rinv=rel(5,n)
                rmainv = 1.d0/(r-par_a)
                s2_t0(n2) = par_cap_A*dexp(par_sig*rmainv)
                s2_t1(n2) = (par_cap_B*rinv)**par_rh
                s2_t2(n2) = par_rh*rinv
                s2_t3(n2) = par_sig*rmainv*rmainv
                s2_dx(n2) = dx
                s2_dy(n2) = dy
                s2_dz(n2) = dz
                 s2_r(n2) = r
                n2 = n2 + 1
                if (n2.gt.max_nbrs) then
                write(6,*) 'WARNING enlarge max_nbrs'
                istop=1
                return
                endif

! coordination number calculated with soft cutoff between first 
! nearest neighbor and midpoint of first and second nearest neighbor
        if (r.le.2.36d0) then
        coord_iat=coord_iat+1.d0
        else if (r.ge.3.12d0) then
        else
        xarg=(r-2.36d0)*(1.d0/(3.12d0-2.36d0))
        coord_iat=coord_iat+(2*xarg+1.d0)*(xarg-1.d0)**2
        endif

                
!   RADIAL PARTS OF THREE-BODY INTERACTION r<par_b
        
                if(r .lt. par_bg)  then
        
                  num3(n3) = j
                  rmbinv = 1.d0/(r-par_bg)
                  temp1 = par_gam*rmbinv
                  temp0 = dexp(temp1)
                  s3_g(n3) = temp0
                  s3_dg(n3) = -rmbinv*temp1*temp0
                  s3_dx(n3) = dx
                  s3_dy(n3) = dy
                  s3_dz(n3) = dz
                  s3_rinv(n3) = rinv
                  s3_r(n3) = r
                  n3 = n3 + 1
                  if (n3.gt.max_nbrs) then
                  write(6,*) 'WARNING enlarge max_nbrs'
                  istop=1
                  return
                  endif
        
                  
!   COORDINATION AND NEIGHBOR FUNCTION par_c<r<par_b
        
                  if(r .lt. par_b) then
                    if(r .lt. par_c) then
                    Z = Z + 1.d0
                   else
                    xinv = bmc/(r-par_c)
                    xinv3 = xinv*xinv*xinv
                    den = 1.d0/(1 - xinv3)
                    temp1 = par_alp*den
                    fZ = dexp(temp1)
                    Z = Z + fZ
                    numz(nz) = j
                    sz_df(nz) = fZ*temp1*den*3.d0*xinv3*xinv*cmbinv   
!   df/dr
                    sz_dx(nz) = dx
                    sz_dy(nz) = dy
                    sz_dz(nz) = dz
                    sz_r(nz) = r
                    nz = nz + 1
                    if (nz.gt.max_nbrs) then
                    write(6,*) 'WARNING enlarge max_nbrs'
                    istop=1
                    return
                    endif      
                   end if 
!  r < par_C
                  end if 
!  r < par_b
                  end if 
!  r < par_bg
              end do
        
              
!   ZERO ACCUMULATION ARRAY FOR ENVIRONMENT FORCES
        
              do nl=1, nz-1
                sz_sum(nl)=0.d0
              end do
        
              
!   ENVIRONMENT-DEPENDENCE OF PAIR INTERACTION
        
              temp0 = par_bet*Z
              pZ = par_palp*dexp(-temp0*Z)         
!   bond order
              dp1 = -2.d0*temp0*pZ         
!   derivative of bond order
        
        
            
!  --- LEVEL 2: LOOP FOR PAIR INTERACTIONS ---
        
            do nj=1, n2-1
        
              temp0 = s2_t1(nj) - pZ
        
              
!   two-body energy V2(rij,Z)
        
              ener_iat = ener_iat + temp0*s2_t0(nj)
              
!   two-body forces
        
              dV2j = - s2_t0(nj) * (s2_t1(nj)*s2_t2(nj) + temp0 * s2_t3(nj))   
!   dV2/dr
              dV2ijx = dV2j * s2_dx(nj)
              dV2ijy = dV2j * s2_dy(nj)
              dV2ijz = dV2j * s2_dz(nj)
              ff(1,i) = ff(1,i) + dV2ijx
              ff(2,i) = ff(2,i) + dV2ijy
              ff(3,i) = ff(3,i) + dV2ijz
              j = num2(nj)
              ff(1,j) = ff(1,j) - dV2ijx
              ff(2,j) = ff(2,j) - dV2ijy
              ff(3,j) = ff(3,j) - dV2ijz
        
              
              
!  --- LEVEL 3: LOOP FOR PAIR COORDINATION FORCES ---
        
              dV2dZ = - dp1 * s2_t0(nj)
              do nl=1, nz-1
                 sz_sum(nl) =  sz_sum(nl) + dV2dZ
              end do
        
            end do
        
              
!   COORDINATION-DEPENDENCE OF THREE-BODY INTERACTION
        
              winv = Qort*exp(-muhalf*Z) 
!   inverse width of angular function
              dwinv = -muhalf*winv       
!   its derivative
              temp0 = exp(-u4*Z)
              tau = u1+u2*temp0*(u3-temp0) 
!   -cosine of angular minimum
              dtau = u5*temp0*(2*temp0-u3) 
!   its derivative
        
            
!  --- LEVEL 2: FIRST LOOP FOR THREE-BODY INTERACTIONS ---
        
            do nj=1, n3-2
        
              j=num3(nj)
        
              
!  --- LEVEL 3: SECOND LOOP FOR THREE-BODY INTERACTIONS ---
        
              do nk=nj+1, n3-1
        
                k=num3(nk)
        
                
!   angular function h(l,Z)
        
                lcos=s3_dx(nj)*s3_dx(nk)+s3_dy(nj)*s3_dy(nk)+s3_dz(nj)*s3_dz(nk)
                x = (lcos + tau)*winv
                temp0 = exp(-x*x)
        
                H = par_lam*(1 - temp0 + par_eta*x*x)
                dHdx = 2*par_lam*x*(temp0 + par_eta)
        
                dhdl = dHdx*winv
        
                
!   three-body energy
        
                temp1 = s3_g(nj) * s3_g(nk)
                ener_iat = ener_iat + temp1*H
        
                
!   (-) radial force on atom j
        
                dV3rij = s3_dg(nj) * s3_g(nk) * H
                dV3rijx = dV3rij * s3_dx(nj)
                dV3rijy = dV3rij * s3_dy(nj)
                dV3rijz = dV3rij * s3_dz(nj)
                fjx = dV3rijx
                fjy = dV3rijy
                fjz = dV3rijz
        
                
!   (-) radial force on atom k
        
                dV3rik = s3_g(nj) * s3_dg(nk) * H
                dV3rikx = dV3rik * s3_dx(nk)
                dV3riky = dV3rik * s3_dy(nk)
                dV3rikz = dV3rik * s3_dz(nk)
                fkx = dV3rikx
                fky = dV3riky
                fkz = dV3rikz
        
                
!   (-) angular force on j
        
                dV3l = temp1*dhdl
                dV3ljx = dV3l * (s3_dx(nk) - lcos * s3_dx(nj)) * s3_rinv(nj)
                dV3ljy = dV3l * (s3_dy(nk) - lcos * s3_dy(nj)) * s3_rinv(nj)
                dV3ljz = dV3l * (s3_dz(nk) - lcos * s3_dz(nj)) * s3_rinv(nj)
                fjx = fjx + dV3ljx
                fjy = fjy + dV3ljy
                fjz = fjz + dV3ljz
        
                
!   (-) angular force on k
        
                dV3lkx = dV3l * (s3_dx(nj) - lcos * s3_dx(nk)) * s3_rinv(nk)
                dV3lky = dV3l * (s3_dy(nj) - lcos * s3_dy(nk)) * s3_rinv(nk)
                dV3lkz = dV3l * (s3_dz(nj) - lcos * s3_dz(nk)) * s3_rinv(nk)
                fkx = fkx + dV3lkx
                fky = fky + dV3lky
                fkz = fkz + dV3lkz
        
                
!   apply radial + angular forces to i, j, k
        
                ff(1,j) = ff(1,j) - fjx
                ff(2,j) = ff(2,j) - fjy
                ff(3,j) = ff(3,j) - fjz
                ff(1,k) = ff(1,k) - fkx
                ff(2,k) = ff(2,k) - fky
                ff(3,k) = ff(3,k) - fkz
                ff(1,i) = ff(1,i) + fjx + fkx
                ff(2,i) = ff(2,i) + fjy + fky
                ff(3,i) = ff(3,i) + fjz + fkz
                
        
                
!   prefactor for 4-body forces from coordination
                  dxdZ = dwinv*(lcos + tau) + winv*dtau
                  dV3dZ = temp1*dHdx*dxdZ
        
                
!  --- LEVEL 4: LOOP FOR THREE-BODY COORDINATION FORCES ---
        
                  do nl=1, nz-1
                    sz_sum(nl) = sz_sum(nl) + dV3dZ
                  end do
              end do
            end do
        
            
!  --- LEVEL 2: LOOP TO APPLY COORDINATION FORCES ---
        
            do nl=1, nz-1
        
                dEdrl = sz_sum(nl) * sz_df(nl)
                dEdrlx = dEdrl * sz_dx(nl)
                dEdrly = dEdrl * sz_dy(nl)
                dEdrlz = dEdrl * sz_dz(nl)
                ff(1,i) = ff(1,i) + dEdrlx
                ff(2,i) = ff(2,i) + dEdrly
                ff(3,i) = ff(3,i) + dEdrlz
                l = numz(nl)
                ff(1,l) = ff(1,l) - dEdrlx
                ff(2,l) = ff(2,l) - dEdrly
                ff(3,l) = ff(3,l) - dEdrlz
        
                
            end do

        coord=coord+coord_iat
        coord2=coord2+coord_iat**2
        ener = ener + ener_iat
        ener2 = ener2 + ener_iat**2
        
1000      continue
        

        return
        end subroutine subfeniat_b



        subroutine sublstiat_b(iat,nn,ncx,ll1,ll2,ll3,l1,l2,l3,myspace, & 
                   rxyz,icell,lstb,lay,rel,cut2,indlst)
! finds the neighbours of atom iat (specified by lsta and lstb) and and
! the relative position rel of iat with respect to these neighbours
        implicit REAL(KIND=dp) (a-h,o-z)
        dimension rxyz(3,nn),lay(nn),icell(0:ncx,-1:ll1,-1:ll2,-1:ll3), &
                  lstb(0:myspace-1),rel(5,0:myspace-1)

        do 6003,k3=l3-1,l3+1
        do 6002,k2=l2-1,l2+1
        do 6001,k1=l1-1,l1+1
        do 6000,jj=1,icell(0,k1,k2,k3)
          jat=icell(jj,k1,k2,k3)
          if (jat.eq.iat) goto 6000
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
6000        continue
6001        continue
6002        continue
6003        continue

        return
        end subroutine sublstiat_b
!****************************************************************************

!!****f* eip_silicon/eip_lenosky_silicon [1.0] *
!!
!!   NAME
!!     eip_lenosky_silicon
!!
!!   FUNCTION
!!     Lenosky's "highly optimized empirical potential model of silicon" 
!!     by Stefan Goedecker
!!
!!   REFERENCES
!!     T. Lenosky, et. al.: Highly optimized empirical potential model of silicon;
!!                          Modelling Simul. Sci. Eng., 8 (2000)
!!     S. Goedecker: Optimization and parallelization of a force field for silicon 
!!                   using OpenMP; CPC 148, 1 (2002)
!!
!!   INPUTS
!!     - nat: number of atoms
!!     - alat: lattice constants of the orthorombic box containing the particles
!!     - rxyz0: atomic positions in Angstrom
!!
!!   OUTPUTS
!!     - rxyz0: rxyz0 may be modified on output.
!!              If an atom is outside the box the program will bring it back 
!!              into the box by translations through alat
!!     - fxyz: forces in eV/A
!!     - ener: total energy in eV
!!     - coord: average coordination number
!!     - ener_var: variance of the energy/atom
!!     - coord_var: variance of the coordination number
!!     - count: count is increased by one per call, has to be initialized 
!               to 0.d0 before first call of eip_bazant
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     03.2006 initial create [tdk]
!!
!!*** **********************************************************************
    subroutine eip_lenosky_silicon(nat,alat,rxyz0,fxyz,ener,coord,ener_var,&
                                   coord_var,count)

        implicit REAL(KIND=dp) (a-h,o-z)
!$      interface
!$        integer ( kind=4 ) function omp_get_num_threads ( )
!$        end function omp_get_num_threads
!$      end interface
!$      interface
!$        integer ( kind=4 ) function omp_get_thread_num ( )
!$        end function omp_get_thread_num
!$      end interface

        dimension rxyz0(3,nat),fxyz(3,nat),alat(3)
        REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:) :: rxyz
        integer, ALLOCATABLE, DIMENSION(:,:)       :: lsta
        integer, ALLOCATABLE, DIMENSION(:)         :: lstb
        integer, ALLOCATABLE, DIMENSION(:)         :: lay
        integer, ALLOCATABLE, DIMENSION(:,:,:,:)   :: icell
        REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:) :: rel
        REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:) :: txyz
        REAL(KIND=dp), ALLOCATABLE, DIMENSION(:,:) :: f2ij,f3ij,f3ik

!        tmax_phi= 0.4500000d+01
!        cut=tmax_phi
        cut= 0.4500000d+01

        if (count.eq.0)  open(unit=10,file='lenosky.mon',status='unknown')
        count=count+1.d0

! linear scaling calculation of verlet list
        ll1=int(alat(1)/cut)
        if (ll1.lt.1) stop 'alat(1) too small'
        ll2=int(alat(2)/cut)
        if (ll2.lt.1) stop 'alat(2) too small'
        ll3=int(alat(3)/cut)
        if (ll3.lt.1) stop 'alat(3) too small'

! determine number of threads
        npr=1
!$omp parallel private(iam)  shared (npr)
!$       iam=omp_get_thread_num()
!$       if (iam.eq.0) npr=omp_get_num_threads()
!$omp end parallel 

! linear scaling calculation of verlet list

     if (npr.le.1) then !serial if too few processors to gain by parallelizing

! set ncx for serial case, ncx for parallel case set below
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

     else  ! parallel case

! periodization of particles can be done in parallel
!$omp parallel do shared (alat,nat,rxyz0) private(iat)
        do 5983,iat=1,nat
        rxyz0(1,iat)=modulo(modulo(rxyz0(1,iat),alat(1)),alat(1))
        rxyz0(2,iat)=modulo(modulo(rxyz0(2,iat),alat(2)),alat(2))
        rxyz0(3,iat)=modulo(modulo(rxyz0(3,iat),alat(3)),alat(3))
5983    continue
!$omp end parallel do

! assignment to cell is done serially
! set ncx for parallel case, ncx for serial case set above
        ncx=8
4321    ncx=ncx*2
        allocate(icell(0:ncx,-1:ll1,-1:ll2,-1:ll3))
        do 3984,l3=-1,ll3
        do 3984,l2=-1,ll2
        do 3984,l1=-1,ll1
3984    icell(0,l1,l2,l3)=0

        rlc1i=ll1/alat(1)
        rlc2i=ll2/alat(2)
        rlc3i=ll3/alat(3)

        do 6983,iat=1,nat
        l1=int(rxyz0(1,iat)*rlc1i)
        l2=int(rxyz0(2,iat)*rlc2i)
        l3=int(rxyz0(3,iat)*rlc3i)
        ii=icell(0,l1,l2,l3)
        ii=ii+1
        icell(0,l1,l2,l3)=ii
        if (ii.gt.ncx) then
        write(10,*) count,'NCX too small',ncx
        deallocate(icell)
        goto 4321
        endif
        icell(ii,l1,l2,l3)=iat
6983    continue

    endif


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
        nnbrx=24
2345    nnbrx=3*nnbrx/2
        allocate(lstb(nnbrx*nat),rel(5,nnbrx*nat))

        indlstx=0

!$omp parallel  &
!$omp private(iat,cut2,iam,ii,indlst,l1,l2,l3,myspace,npr) &
!$omp shared (indlstx,nat,nn,nnbrx,ncx,ll1,ll2,ll3,icell,lsta,lstb,lay, & 
!$omp rel,rxyz,cut,myspaceout)


        npr=1
!$       npr=omp_get_num_threads()
        iam=0
!$       iam=omp_get_thread_num()

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
        call sublstiat_l(iat,nn,ncx,ll1,ll2,ll3,l1,l2,l3,myspace, & 
             rxyz,icell,lstb(iam*myspace+1),lay,rel(1,iam*myspace+1),cut2,indlst)
        lsta(2,iat)=iam*myspace+indlst
!        write(6,'(a,4(x,i3),100(x,i2))') & 
!               'iam,iat,lsta',iam,iat,lsta(1,iat),lsta(2,iat), & 
!                    (lstb(j),j=lsta(1,iat),lsta(2,iat))
        endif

6600    continue
6000    continue
!$omp critical
        indlstx=max(indlstx,indlst)
!$omp end critical
!$omp end parallel  

           if (indlstx.ge.myspaceout) then
               write(10,*) count,'NNBRX too  small', nnbrx
               deallocate(lstb,rel)
               goto 2345
           endif

!$omp parallel  &
!$omp private(iam,npr,iat,iat1,iat2,lot,istop,tcoord,tcoord2, & 
!$omp tener,tener2,txyz,f2ij,f3ij,f3ik,npjx,npjkx) &
!$omp shared (nat,nnbrx,lsta,lstb,rel,ener,ener2,fxyz,coord,coord2,istopg)

        npr=1
!$       npr=omp_get_num_threads()
        iam=0
!$       iam=omp_get_thread_num()

        npjx=300 ; npjkx=6000
        istopg=0

        if (npr.ne.1) then 
! PARALLEL CASE
! create temporary private scalars for reduction sum on energies and 
!        temporary private array for reduction sum on forces
!$omp critical
        allocate(txyz(3,nat),f2ij(3,npjx),f3ij(3,npjkx),f3ik(3,npjkx))
!$omp end critical
        if (iam.eq.0) then
        ener=0.d0
        ener2=0.d0
        coord=0.d0
        coord2=0.d0
        endif
!$omp do
        do 121,iat=1,nat
        fxyz(1,iat)=0.d0
        fxyz(2,iat)=0.d0
121     fxyz(3,iat)=0.d0
!$omp barrier    

! Each thread treats at most lot atoms
        lot=int(float(nat)/float(npr)+.999999999999d0)
        iat1=iam*lot+1
        iat2=min((iam+1)*lot,nat)
!       write(6,*) 'subfeniat:iat1,iat2,iam',iat1,iat2,iam
        call subfeniat_l(iat1,iat2,nat,lsta,lstb,rel,tener,tener2, &
               tcoord,tcoord2,nnbrx,txyz,f2ij,npjx,f3ij,npjkx,f3ik,istop)
!$omp critical
        ener=ener+tener
        ener2=ener2+tener2
        coord=coord+tcoord
        coord2=coord2+tcoord2
        istopg=istopg+istop
        do 8093,iat=1,nat
        fxyz(1,iat)=fxyz(1,iat)+txyz(1,iat)
        fxyz(2,iat)=fxyz(2,iat)+txyz(2,iat)
        fxyz(3,iat)=fxyz(3,iat)+txyz(3,iat)
8093    continue
        deallocate(txyz,f2ij,f3ij,f3ik)
!$omp end critical

        else
! SERIAL CASE
        iat1=1
        iat2=nat
        allocate(f2ij(3,npjx),f3ij(3,npjkx),f3ik(3,npjkx))
        call subfeniat_l(iat1,iat2,nat,lsta,lstb,rel,ener,ener2, &
               coord,coord2,nnbrx,fxyz,f2ij,npjx,f3ij,npjkx,f3ik,istopg)
        deallocate(f2ij,f3ij,f3ik)

        endif
!$omp end parallel

!         write(6,*) 'ener,norm force', & 
!                    ener,DNRM2(3*nat,fxyz,1)
        if (istopg.gt.0) stop 'DIMENSION ERROR (see WARNING above)'
        ener_var=ener2/nat-(ener/nat)**2 
        coord=coord/nat
        coord_var=coord2/nat-coord**2 

        deallocate(rxyz,icell,lay,lsta,lstb,rel)

        end subroutine eip_lenosky_silicon


        subroutine subfeniat_l(iat1,iat2,nat,lsta,lstb,rel,tener,tener2, &
               tcoord,tcoord2,nnbrx,txyz,f2ij,npjx,f3ij,npjkx,f3ik,istop)
! for a subset of atoms iat1 to iat2 the routine calculates the (partial) forces
! txyz acting on these atoms as well as on the atoms (jat, kat) interacting
! with them and their contribution to the energy (tener).
! In addition the coordination number tcoord and the second moment of the
! local energy tener2 and coordination number tcoord2 are returned
        implicit REAL(KIND=dp) (a-h,o-z)
        dimension lsta(2,nat),lstb(nnbrx*nat),rel(5,nnbrx*nat),txyz(3,nat)
        dimension f2ij(3,npjx),f3ij(3,npjkx),f3ik(3,npjkx)
        REAL(KIND=dp) :: tmin_phi= 0.1500000d+01  
        REAL(KIND=dp) :: tmax_phi= 0.4500000d+01 
        REAL(KIND=dp) :: hi_phi= 3.00000000000d0 
        REAL(KIND=dp) :: hsixth_phi=5.55555555555556D-002 
        REAL(KIND=dp) :: h2sixth_phi=1.85185185185185D-002
        REAL(KIND=dp), dimension(0:9) :: cof_phi  = &
                     (/ 0.69299400000000d+01, -0.43995000000000d+00, &
                       -0.17012300000000d+01, -0.16247300000000d+01, &
                       -0.99696000000000d+00, -0.27391000000000d+00, &
                       -0.24990000000000d-01, -0.17840000000000d-01, &
                       -0.96100000000000d-02,  0.00000000000000d+00 /)
        REAL(KIND=dp), dimension(0:9) :: dof_phi  = &
                     (/ 0.16533229480429d+03,  0.39415410391417d+02, &
                        0.68710036300407d+01,  0.53406950884203d+01, &
                        0.15347960162782d+01, -0.63347591535331d+01, &
                       -0.17987794021458d+01,  0.47429676211617d+00, &
                       -0.40087646318907d-01, -0.23942617684055d+00 /)
        REAL(KIND=dp) :: tmin_rho= 0.1500000d+01  
        REAL(KIND=dp) :: tmax_rho= 0.3500000d+01 
        REAL(KIND=dp) :: hi_rho= 5.00000000000d0 
        REAL(KIND=dp) :: hsixth_rho=3.33333333333333D-002 
        REAL(KIND=dp) :: h2sixth_rho=6.66666666666667D-003
        REAL(KIND=dp), dimension(0:10) :: cof_rho =  &
                     (/ 0.13747000000000d+00, -0.14831000000000d+00, &
                       -0.55972000000000d+00, -0.73110000000000d+00, &
                       -0.76283000000000d+00, -0.72918000000000d+00, &
                       -0.66620000000000d+00, -0.57328000000000d+00, &
                       -0.40690000000000d+00, -0.16662000000000d+00, &
                        0.00000000000000d+00 /)
        REAL(KIND=dp), dimension(0:10) :: dof_rho =  &
                    (/ -0.32275496741918d+01, -0.64119006516165d+01, &
                        0.10030652280658d+02,  0.22937915289857d+01, &
                        0.17416816033995d+01,  0.54648205741626d+00, &
                        0.47189016693543d+00,  0.20569572748420d+01, &
                        0.23192807336964d+01, -0.24908020962757d+00, &
                       -0.12371959895186d+02 /)
        REAL(KIND=dp) :: tmin_fff= 0.1500000d+01  
        REAL(KIND=dp) :: tmax_fff= 0.3500000d+01 
        REAL(KIND=dp) :: hi_fff= 4.50000000000d0 
        REAL(KIND=dp) :: hsixth_fff=3.70370370370370D-002 
        REAL(KIND=dp) :: h2sixth_fff=8.23045267489712D-003
        REAL(KIND=dp), dimension(0:9) :: cof_fff  = &
                     (/ 0.12503100000000d+01,  0.86821000000000d+00, &
                        0.60846000000000d+00,  0.48756000000000d+00, &
                        0.44163000000000d+00,  0.37610000000000d+00, &
                        0.27145000000000d+00,  0.14814000000000d+00, &
                        0.48550000000000d-01,  0.00000000000000d+00 /)
        REAL(KIND=dp), dimension(0:9) :: dof_fff  = &
                     (/ 0.27904652711432d+02, -0.45230754228635d+01, &
                        0.50531739800222d+01,  0.11806545027747d+01, &
                       -0.66693699112098d+00, -0.89430653829079d+00, &
                       -0.50891685571587d+00,  0.66278396115427d+00, &
                        0.73976101109878d+00,  0.25795319944506d+01 /)
        REAL(KIND=dp) :: tmin_uuu= -0.1770930d+01  
        REAL(KIND=dp) :: tmax_uuu= 0.7908520d+01 
        REAL(KIND=dp) :: hi_uuu= 0.723181585730594d0 
        REAL(KIND=dp) :: hsixth_uuu=0.230463095238095d0 
        REAL(KIND=dp) :: h2sixth_uuu=0.318679429600340d0
        REAL(KIND=dp), dimension(0:7) :: cof_uuu  = &
                    (/ -0.10749300000000d+01, -0.20045000000000d+00, &
                        0.41422000000000d+00,  0.87939000000000d+00, &
                        0.12668900000000d+01,  0.16299800000000d+01, &
                        0.19773800000000d+01,  0.23961800000000d+01 /)
        REAL(KIND=dp), dimension(0:7) :: dof_uuu  = &
                    (/ -0.14827125747284d+00, -0.14922155328475d+00, &
                       -0.70113224223509d-01, -0.39449020349230d-01, &
                       -0.15815242579643d-01,  0.26112640061855d-01, &
                       -0.13786974745095d+00,  0.74941595372657d+00 /)
        REAL(KIND=dp) :: tmin_ggg= -0.1000000d+01  
        REAL(KIND=dp) :: tmax_ggg= 0.8001400d+00 
        REAL(KIND=dp) :: hi_ggg= 3.88858644327663d0 
        REAL(KIND=dp) :: hsixth_ggg=4.28604761904762D-002 
        REAL(KIND=dp) :: h2sixth_ggg=1.10221225156463D-002
        REAL(KIND=dp), dimension(0:7) :: cof_ggg  = &
                     (/ 0.52541600000000d+01,  0.23591500000000d+01, &
                        0.11959500000000d+01,  0.12299500000000d+01, &
                        0.20356500000000d+01,  0.34247400000000d+01, &
                        0.49485900000000d+01,  0.56179900000000d+01 /)
        REAL(KIND=dp), dimension(0:7) :: dof_ggg  = &
                     (/ 0.15826876132396d+02,  0.31176239377907d+02, &
                        0.16589446539683d+02,  0.11083892500520d+02, &
                        0.90887216383860d+01,  0.54902279653967d+01, &
                       -0.18823313223755d+02, -0.77183416481005d+01 /)

! initialize temporary private scalars for reduction sum on energies and 
! private workarray txyz for forces forces
        tener=0.d0
        tener2=0.d0
        tcoord=0.d0
        tcoord2=0.d0
        istop=0
        do 121,iat=1,nat
        txyz(1,iat)=0.d0
        txyz(2,iat)=0.d0
121     txyz(3,iat)=0.d0


! calculation of forces, energy

        do 1000,iat=iat1,iat2

        dens2=0.d0
        dens3=0.d0
        jcnt=0
        jkcnt=0
        coord_iat=0.d0
        ener_iat=0.d0
        do 2000,jbr=lsta(1,iat),lsta(2,iat)
        jat=lstb(jbr)
        jcnt=jcnt+1
        if (jcnt.gt.npjx) then 
            write(6,*) 'WARNING: enlarge npjx'  
            istop=1
            return
        endif

        fxij=rel(1,jbr)
        fyij=rel(2,jbr)
        fzij=rel(3,jbr)
        rij=rel(4,jbr)
        sij=rel(5,jbr)

! coordination number calculated with soft cutoff between first
! nearest neighbor and midpoint of first and second nearest neighbor
        if (rij.le.2.36d0) then
        coord_iat=coord_iat+1.d0
        else if (rij.ge.3.12d0) then
        else
        xarg=(rij-2.36d0)*(1.d0/(3.12d0-2.36d0))
        coord_iat=coord_iat+(2*xarg+1.d0)*(xarg-1.d0)**2
        endif

! pairpotential term        
        call splint(cof_phi,dof_phi,tmin_phi,tmax_phi, &
           hsixth_phi,h2sixth_phi,hi_phi,10,rij,e_phi,ep_phi)
        ener_iat=ener_iat+(e_phi*.5d0)
        txyz(1,iat)=txyz(1,iat)-fxij*(ep_phi*.5d0)
        txyz(2,iat)=txyz(2,iat)-fyij*(ep_phi*.5d0)
        txyz(3,iat)=txyz(3,iat)-fzij*(ep_phi*.5d0)
        txyz(1,jat)=txyz(1,jat)+fxij*(ep_phi*.5d0)
        txyz(2,jat)=txyz(2,jat)+fyij*(ep_phi*.5d0)
        txyz(3,jat)=txyz(3,jat)+fzij*(ep_phi*.5d0)

! 2 body embedding term        
        call splint(cof_rho,dof_rho,tmin_rho,tmax_rho, &
             hsixth_rho,h2sixth_rho,hi_rho,11,rij,rho,rhop)
        dens2=dens2+rho
        f2ij(1,jcnt)=fxij*rhop
        f2ij(2,jcnt)=fyij*rhop
        f2ij(3,jcnt)=fzij*rhop

! 3 body embedding term        
        call splint(cof_fff,dof_fff,tmin_fff,tmax_fff, & 
             hsixth_fff,h2sixth_fff,hi_fff,10,rij,fij,fijp)

        do 3000,kbr=lsta(1,iat),lsta(2,iat)
        kat=lstb(kbr)
        if (kat.lt.jat) then
        jkcnt=jkcnt+1
        if (jkcnt.gt.npjkx) then 
            write(6,*) 'WARNING: enlarge npjkx',npjkx 
            istop=1
            return
        endif

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
        rik=rel(4,kbr)
      if (rik.gt.tmax_fff) then
        fikp=0.d0 ; fik=0.d0
        gjik=0.d0 ;  gjikp=0.d0 ; sik=0.d0
        costheta=0.d0 ; fxik=0.d0 ; fyik=0.d0 ; fzik=0.d0
      else if (rik.lt.tmin_fff) then
        fxik=rel(1,kbr)
        fyik=rel(2,kbr)
        fzik=rel(3,kbr)
        costheta=fxij*fxik+fyij*fyik+fzij*fzik
        sik=rel(5,kbr)
        fikp=hi_fff*(cof_fff(1)-cof_fff(0)) -  & 
             ( dof_fff(1)+2.d0*dof_fff(0) )*hsixth_fff
        fik=cof_fff(0) + (rik-tmin_fff)*fikp
        tt_ggg=(costheta-tmin_ggg)*hi_ggg
        if (costheta.gt.tmax_ggg) then
           gjikp=hi_ggg*(cof_ggg(8-1)-cof_ggg(8-2)) + & 
                ( 2.d0*dof_ggg(8-1)+dof_ggg(8-2) )*hsixth_ggg
                 gjik=cof_ggg(8-1) + (costheta-tmax_ggg)*gjikp
        else
           klo_ggg=tt_ggg
           khi_ggg=klo_ggg+1
           cof_ggg_klo=cof_ggg(klo_ggg)
           dof_ggg_klo=dof_ggg(klo_ggg)
           b_ggg=tt_ggg-klo_ggg
           a_ggg=1.d0-b_ggg
           cof_ggg_khi=cof_ggg(khi_ggg)
           dof_ggg_khi=dof_ggg(khi_ggg)
           b2_ggg=b_ggg*b_ggg
           gjik=a_ggg*cof_ggg_klo
           gjikp=cof_ggg_khi-cof_ggg_klo 
           a2_ggg=a_ggg*a_ggg
           cof1_ggg=a2_ggg-1.d0
           cof2_ggg=b2_ggg-1.d0
           gjik=gjik+b_ggg*cof_ggg_khi
           gjikp=hi_ggg*gjikp
           cof3_ggg=3.d0*b2_ggg
           cof4_ggg=3.d0*a2_ggg
           cof1_ggg=a_ggg*cof1_ggg
           cof2_ggg=b_ggg*cof2_ggg
           cof3_ggg=cof3_ggg-1.d0
           cof4_ggg=cof4_ggg-1.d0
           yt1_ggg=cof1_ggg*dof_ggg_klo
           yt2_ggg=cof2_ggg*dof_ggg_khi
           ypt1_ggg=cof3_ggg*dof_ggg_khi
           ypt2_ggg=cof4_ggg*dof_ggg_klo
           gjik=gjik + (yt1_ggg+yt2_ggg)*h2sixth_ggg
           gjikp=gjikp + ( ypt1_ggg - ypt2_ggg )*hsixth_ggg
        endif
      else
        fxik=rel(1,kbr)
        tt_fff=rik-tmin_fff
        costheta=fxij*fxik
        fyik=rel(2,kbr)
        tt_fff=tt_fff*hi_fff
        costheta=costheta+fyij*fyik
        fzik=rel(3,kbr)
        klo_fff=tt_fff
        costheta=costheta+fzij*fzik
        sik=rel(5,kbr)
        tt_ggg=(costheta-tmin_ggg)*hi_ggg
        if (costheta.gt.tmax_ggg) then
          gjikp=hi_ggg*(cof_ggg(8-1)-cof_ggg(8-2)) + & 
                ( 2.d0*dof_ggg(8-1)+dof_ggg(8-2) )*hsixth_ggg
                gjik=cof_ggg(8-1) + (costheta-tmax_ggg)*gjikp
          khi_fff=klo_fff+1
          cof_fff_klo=cof_fff(klo_fff)
          dof_fff_klo=dof_fff(klo_fff)
          b_fff=tt_fff-klo_fff
          a_fff=1.d0-b_fff
          cof_fff_khi=cof_fff(khi_fff)
          dof_fff_khi=dof_fff(khi_fff)
          b2_fff=b_fff*b_fff
          fik=a_fff*cof_fff_klo
          fikp=cof_fff_khi-cof_fff_klo 
          a2_fff=a_fff*a_fff
          cof1_fff=a2_fff-1.d0
          cof2_fff=b2_fff-1.d0
          fik=fik+b_fff*cof_fff_khi
          fikp=hi_fff*fikp
          cof3_fff=3.d0*b2_fff
          cof4_fff=3.d0*a2_fff
          cof1_fff=a_fff*cof1_fff
          cof2_fff=b_fff*cof2_fff
          cof3_fff=cof3_fff-1.d0
          cof4_fff=cof4_fff-1.d0
          yt1_fff=cof1_fff*dof_fff_klo
          yt2_fff=cof2_fff*dof_fff_khi
          ypt1_fff=cof3_fff*dof_fff_khi
          ypt2_fff=cof4_fff*dof_fff_klo
          fik=fik + (yt1_fff+yt2_fff)*h2sixth_fff
          fikp=fikp + ( ypt1_fff - ypt2_fff )*hsixth_fff
         else
              klo_ggg=tt_ggg
              khi_ggg=klo_ggg+1
           khi_fff=klo_fff+1
              cof_ggg_klo=cof_ggg(klo_ggg)
           cof_fff_klo=cof_fff(klo_fff)
              dof_ggg_klo=dof_ggg(klo_ggg)
           dof_fff_klo=dof_fff(klo_fff)
              b_ggg=tt_ggg-klo_ggg
           b_fff=tt_fff-klo_fff
              a_ggg=1.d0-b_ggg
           a_fff=1.d0-b_fff
              cof_ggg_khi=cof_ggg(khi_ggg)
           cof_fff_khi=cof_fff(khi_fff)
              dof_ggg_khi=dof_ggg(khi_ggg)
           dof_fff_khi=dof_fff(khi_fff)
              b2_ggg=b_ggg*b_ggg
           b2_fff=b_fff*b_fff
              gjik=a_ggg*cof_ggg_klo
           fik=a_fff*cof_fff_klo
              gjikp=cof_ggg_khi-cof_ggg_klo 
           fikp=cof_fff_khi-cof_fff_klo 
              a2_ggg=a_ggg*a_ggg
           a2_fff=a_fff*a_fff
              cof1_ggg=a2_ggg-1.d0
           cof1_fff=a2_fff-1.d0
              cof2_ggg=b2_ggg-1.d0
           cof2_fff=b2_fff-1.d0
              gjik=gjik+b_ggg*cof_ggg_khi
           fik=fik+b_fff*cof_fff_khi
              gjikp=hi_ggg*gjikp
           fikp=hi_fff*fikp
              cof3_ggg=3.d0*b2_ggg
           cof3_fff=3.d0*b2_fff
              cof4_ggg=3.d0*a2_ggg
           cof4_fff=3.d0*a2_fff
              cof1_ggg=a_ggg*cof1_ggg
           cof1_fff=a_fff*cof1_fff
              cof2_ggg=b_ggg*cof2_ggg
           cof2_fff=b_fff*cof2_fff
              cof3_ggg=cof3_ggg-1.d0
           cof3_fff=cof3_fff-1.d0
              cof4_ggg=cof4_ggg-1.d0
           cof4_fff=cof4_fff-1.d0
              yt1_ggg=cof1_ggg*dof_ggg_klo
           yt1_fff=cof1_fff*dof_fff_klo
              yt2_ggg=cof2_ggg*dof_ggg_khi
           yt2_fff=cof2_fff*dof_fff_khi
              ypt1_ggg=cof3_ggg*dof_ggg_khi
           ypt1_fff=cof3_fff*dof_fff_khi
              ypt2_ggg=cof4_ggg*dof_ggg_klo
           ypt2_fff=cof4_fff*dof_fff_klo
              gjik=gjik + (yt1_ggg+yt2_ggg)*h2sixth_ggg
           fik=fik + (yt1_fff+yt2_fff)*h2sixth_fff
              gjikp=gjikp + ( ypt1_ggg - ypt2_ggg )*hsixth_ggg
           fikp=fikp + ( ypt1_fff - ypt2_fff )*hsixth_fff
         endif
      endif
! end optimized version

        tt=fij*fik
        dens3=dens3+tt*gjik

        t1=fijp*fik*gjik
        t2=sij*(tt*gjikp)
        f3ij(1,jkcnt)=fxij*t1 + (fxik-fxij*costheta)*t2
        f3ij(2,jkcnt)=fyij*t1 + (fyik-fyij*costheta)*t2
        f3ij(3,jkcnt)=fzij*t1 + (fzik-fzij*costheta)*t2

        t3=fikp*fij*gjik
        t4=sik*(tt*gjikp)
        f3ik(1,jkcnt)=fxik*t3 + (fxij-fxik*costheta)*t4
        f3ik(2,jkcnt)=fyik*t3 + (fyij-fyik*costheta)*t4
        f3ik(3,jkcnt)=fzik*t3 + (fzij-fzik*costheta)*t4
        endif
3000        continue
2000        continue
        dens=dens2+dens3
        call splint(cof_uuu,dof_uuu,tmin_uuu,tmax_uuu, & 
             hsixth_uuu,h2sixth_uuu,hi_uuu,8,dens,e_uuu,ep_uuu)
        ener_iat=ener_iat+e_uuu

! Only now ep_uu is known and the forces can be calculated, lets loop again
        jcnt=0
        jkcnt=0
        do 2200,jbr=lsta(1,iat),lsta(2,iat)
        jat=lstb(jbr)
        jcnt=jcnt+1
        txyz(1,iat)=txyz(1,iat)-ep_uuu*f2ij(1,jcnt)
        txyz(2,iat)=txyz(2,iat)-ep_uuu*f2ij(2,jcnt)
        txyz(3,iat)=txyz(3,iat)-ep_uuu*f2ij(3,jcnt)
        txyz(1,jat)=txyz(1,jat)+ep_uuu*f2ij(1,jcnt)
        txyz(2,jat)=txyz(2,jat)+ep_uuu*f2ij(2,jcnt)
        txyz(3,jat)=txyz(3,jat)+ep_uuu*f2ij(3,jcnt)
        
! 3 body embedding term        
        do 3300,kbr=lsta(1,iat),lsta(2,iat)
        kat=lstb(kbr)
        if (kat.lt.jat) then
        jkcnt=jkcnt+1

        txyz(1,iat)=txyz(1,iat)-ep_uuu*(f3ij(1,jkcnt)+f3ik(1,jkcnt))
        txyz(2,iat)=txyz(2,iat)-ep_uuu*(f3ij(2,jkcnt)+f3ik(2,jkcnt))
        txyz(3,iat)=txyz(3,iat)-ep_uuu*(f3ij(3,jkcnt)+f3ik(3,jkcnt))
        txyz(1,jat)=txyz(1,jat)+ep_uuu*f3ij(1,jkcnt)
        txyz(2,jat)=txyz(2,jat)+ep_uuu*f3ij(2,jkcnt)
        txyz(3,jat)=txyz(3,jat)+ep_uuu*f3ij(3,jkcnt)
        txyz(1,kat)=txyz(1,kat)+ep_uuu*f3ik(1,jkcnt)
        txyz(2,kat)=txyz(2,kat)+ep_uuu*f3ik(2,jkcnt)
        txyz(3,kat)=txyz(3,kat)+ep_uuu*f3ik(3,jkcnt)


        endif
3300        continue
2200        continue

!        write(6,'(a,i4,x,e19.12,x,e10.3)') 'iat,ener_iat,coord_iat', & 
!                                       iat,ener_iat,coord_iat
        tener=tener+ener_iat
        tener2=tener2+ener_iat**2
        tcoord=tcoord+coord_iat
        tcoord2=tcoord2+coord_iat**2

1000        continue

        return
        end subroutine subfeniat_l


        subroutine sublstiat_l(iat,nn,ncx,ll1,ll2,ll3,l1,l2,l3,myspace, & 
                   rxyz,icell,lstb,lay,rel,cut2,indlst)
! finds the neighbours of atom iat (specified by lsta and lstb) and and
! the relative position rel of iat with respect to these neighbours
        implicit REAL(KIND=dp) (a-h,o-z)
        dimension rxyz(3,nn),lay(nn),icell(0:ncx,-1:ll1,-1:ll2,-1:ll3), &
                  lstb(0:myspace-1),rel(5,0:myspace-1)

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
        end subroutine sublstiat_l

        subroutine splint(ya,y2a,tmin,tmax,hsixth,h2sixth,hi,n,x,y,yp)
        implicit REAL(KIND=dp) (a-h,o-z)
        dimension y2a(0:n-1),ya(0:n-1)
      
! interpolate if the argument is outside the cubic spline interval [tmin,tmax]
        tt=(x-tmin)*hi
        if (x.lt.tmin) then
          yp=hi*(ya(1)-ya(0)) -  &
          ( y2a(1)+2.d0*y2a(0) )*hsixth
          y=ya(0) + (x-tmin)*yp
        else if (x.gt.tmax) then
          yp=hi*(ya(n-1)-ya(n-2)) +  &
          ( 2.d0*y2a(n-1)+y2a(n-2) )*hsixth
          y=ya(n-1) + (x-tmax)*yp
! otherwise evaluate cubic spline
        else
          klo=tt
          khi=klo+1
          ya_klo=ya(klo)
          y2a_klo=y2a(klo)
          b=tt-klo
          a=1.d0-b
          ya_khi=ya(khi)
          y2a_khi=y2a(khi)
          b2=b*b
          y=a*ya_klo
          yp=ya_khi-ya_klo 
          a2=a*a
          cof1=a2-1.d0
          cof2=b2-1.d0
          y=y+b*ya_khi
          yp=hi*yp
          cof3=3.d0*b2
          cof4=3.d0*a2
          cof1=a*cof1
          cof2=b*cof2
          cof3=cof3-1.d0
          cof4=cof4-1.d0
          yt1=cof1*y2a_klo
          yt2=cof2*y2a_khi
          ypt1=cof3*y2a_khi
          ypt2=cof4*y2a_klo
          y=y + (yt1+yt2)*h2sixth
          yp=yp + ( ypt1 - ypt2 )*hsixth
        endif
      return
      end subroutine splint
!****************************************************************************

END MODULE eip_silicon
