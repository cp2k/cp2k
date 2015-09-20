!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      Torsions added(DG)05-Dec-2000
!>      Variable names changed(DG)05-Dec-2000
!>      CJM SEPT-12-2002: int_env is now passed
!>      CJM NOV-30-2003: only uses fist_env
!> \author CJM & JGH
! *****************************************************************************
MODULE fist_force
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE atprop_types,                    ONLY: atprop_type
  USE cell_types,                      ONLY: cell_type,&
                                             init_cell,&
                                             pbc
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_result_methods,               ONLY: cp_results_erase,&
                                             put_results
  USE cp_result_types,                 ONLY: cp_result_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE cp_units,                        ONLY: cp_unit_from_cp2k
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE ewald_environment_types,         ONLY: ewald_env_get,&
                                             ewald_environment_type
  USE ewald_pw_methods,                ONLY: ewald_pw_grid_update
  USE ewald_pw_types,                  ONLY: ewald_pw_type
  USE ewalds,                          ONLY: ewald_evaluate,&
                                             ewald_print,&
                                             ewald_self,&
                                             ewald_self_atom
  USE ewalds_multipole,                ONLY: ewald_multipole_evaluate
  USE exclusion_types,                 ONLY: exclusion_type
  USE fist_energy_types,               ONLY: fist_energy_type
  USE fist_environment_types,          ONLY: fist_env_get,&
                                             fist_environment_type
  USE fist_intra_force,                ONLY: force_intra_control
  USE fist_neighbor_list_control,      ONLY: list_control
  USE fist_nonbond_env_types,          ONLY: fist_nonbond_env_type
  USE fist_nonbond_force,              ONLY: bonded_correct_gaussian,&
                                             force_nonbond
  USE fist_pol_scf,                    ONLY: fist_pol_evaluate
  USE input_constants,                 ONLY: do_fist_pol_none
  USE input_section_types,             ONLY: section_get_ival,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE manybody_eam,                    ONLY: density_nonbond
  USE manybody_potential,              ONLY: energy_manybody,&
                                             force_nonbond_manybody
  USE mathconstants,                   ONLY: twopi
  USE message_passing,                 ONLY: mp_sum
  USE molecule_kind_types,             ONLY: molecule_kind_type
  USE molecule_types_new,              ONLY: molecule_type
  USE moments_utils,                   ONLY: get_reference_point
  USE multipole_types,                 ONLY: multipole_type
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: debye
  USE pme,                             ONLY: pme_evaluate
  USE pw_poisson_types,                ONLY: do_ewald_ewald,&
                                             do_ewald_none,&
                                             do_ewald_pme,&
                                             do_ewald_spme
  USE shell_potential_types,           ONLY: shell_kind_type
  USE spme,                            ONLY: spme_evaluate
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE virial_types,                    ONLY: virial_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: fist_calc_energy_force
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'fist_force'

! *****************************************************************************
  TYPE debug_variables_type
     REAL(KIND=dp)                          :: pot_manybody, pot_bend, pot_torsion
     REAL(KIND=dp)                          :: pot_nonbond, pot_g, pot_bond
     REAL(KIND=dp)                          :: pot_imptors, pot_urey_bradley
     REAL(KIND=dp)                          :: pot_opbend
     REAL(KIND=dp), DIMENSION(:,:), POINTER :: f_nonbond, f_g, f_bond, f_bend,&
                                                 f_torsion, f_imptors, f_ub,&
                                                 f_opbend
     REAL(KIND=dp), DIMENSION(3,3)         :: pv_nonbond, pv_g, pv_bond, pv_ub,&
                                                 pv_bend, pv_torsion, pv_imptors,&
                                                 pv_opbend
  END TYPE debug_variables_type

CONTAINS

! *****************************************************************************
!> \brief Calculates the total potential energy, total force, and the
!>      total pressure tensor from the potentials
!> \param fist_env ...
!> \param debug ...
!> \par History
!>      Harald Forbert(Dec-2000): Changes for multiple linked lists
!>      cjm, 20-Feb-2001: box_ref used to initialize ewald.  Now
!>                        have consistent restarts with npt and ewald
!>      JGH(15-Mar-2001): box_change replaces ensemble parameter
!>                          Call ewald_setup if box has changed
!>                          Consistent setup for PME and SPME
!>      cjm, 28-Feb-2006: box_change is gone
!> \author CJM & JGH
! *****************************************************************************
  SUBROUTINE fist_calc_energy_force(fist_env, debug)
    TYPE(fist_environment_type), POINTER     :: fist_env
    TYPE(debug_variables_type), &
      INTENT(INOUT), OPTIONAL                :: debug

    CHARACTER(len=*), PARAMETER :: routineN = 'fist_calc_energy_force', &
      routineP = moduleN//':'//routineN

    INTEGER :: do_ipol, ewald_type, fg_coulomb_size, handle, i, ii, ikind, &
      iparticle_kind, iparticle_local, iw, iw2, j, natoms, nlocal_particles, &
      node, nparticle_kind, nparticle_local, nshell, shell_index, stat
    LOGICAL                                  :: do_multipoles, failure, &
                                                shell_model_ad, &
                                                shell_present, use_virial
    REAL(KIND=dp) :: mass, massc, masss, pot_bend, pot_bond, pot_imptors, &
      pot_manybody, pot_nonbond, pot_opbend, pot_shell, pot_torsion, &
      pot_urey_bradley, umass, vg_coulomb, xdum1, xdum2, xdum3
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :) :: f_nonbond, f_total, &
      fcore_nonbond, fcore_total, fg_coulomb, fgcore_coulomb, &
      fgshell_coulomb, fshell_nonbond, fshell_total
    REAL(KIND=dp), DIMENSION(3, 3) :: ident, pv_bc, pv_bend, pv_bond, pv_g, &
      pv_imptors, pv_nonbond, pv_opbend, pv_torsion, pv_urey_bradley
    REAL(KIND=dp), DIMENSION(:), POINTER     :: e_coulomb
    REAL(KIND=dp), DIMENSION(:, :, :), &
      POINTER                                :: pv_coulomb
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(atprop_type), POINTER               :: atprop_env
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(distribution_1d_type), POINTER      :: local_molecules, &
                                                local_particles
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(exclusion_type), DIMENSION(:), &
      POINTER                                :: exclusions
    TYPE(fist_energy_type), POINTER          :: thermo
    TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
    TYPE(molecule_kind_type), DIMENSION(:), &
      POINTER                                :: molecule_kind_set
    TYPE(molecule_type), DIMENSION(:), &
      POINTER                                :: molecule_set
    TYPE(multipole_type), POINTER            :: multipoles
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: core_particle_set, &
                                                particle_set, &
                                                shell_particle_set
    TYPE(section_vals_type), POINTER         :: force_env_section, &
                                                mm_section, print_section
    TYPE(shell_kind_type), POINTER           :: shell
    TYPE(virial_type), POINTER               :: virial

    CALL timeset(routineN,handle)
    failure = .FALSE.
    NULLIFY(logger)
    NULLIFY(subsys,virial,atprop_env,para_env,force_env_section)
    logger => cp_get_default_logger()

    CALL fist_env_get(fist_env,&
                      subsys=subsys,&
                      para_env=para_env,&
                      input=force_env_section)

    CALL cp_subsys_get(subsys,&
                       virial=virial,&
                       atprop=atprop_env)

    use_virial = virial%pv_availability.AND.(.NOT.virial%pv_numer)
    NULLIFY(atomic_kind,atomic_kind_set,cell,ewald_pw,ewald_env,&
         fist_nonbond_env,mm_section,local_molecules,local_particles,&
         molecule_kind_set,molecule_set,particle_set, print_section,&
         shell, shell_particle_set, core_particle_set, thermo, multipoles,&
         e_coulomb, pv_coulomb)

    mm_section => section_vals_get_subs_vals(force_env_section,"MM")
    iw = cp_print_key_unit_nr(logger,mm_section,"PRINT%DERIVATIVES",&
         extension=".mmLog")
    iw2= cp_print_key_unit_nr(logger,mm_section,"PRINT%EWALD_INFO",&
         extension=".mmLog")

    CALL fist_env_get(fist_env, ewald_pw=ewald_pw, ewald_env=ewald_env, &
         local_particles=local_particles, particle_set=particle_set, &
         atomic_kind_set=atomic_kind_set, molecule_set=molecule_set, &
         local_molecules=local_molecules, thermo=thermo, &
         molecule_kind_set=molecule_kind_set, fist_nonbond_env=fist_nonbond_env,&
         cell=cell, shell_model=shell_present, shell_model_ad=shell_model_ad, &
         multipoles=multipoles, exclusions=exclusions)

    CALL ewald_env_get(ewald_env, ewald_type=ewald_type, do_multipoles=do_multipoles,&
         do_ipol=do_ipol)

    ! Initialize ewald grids
    CALL init_cell(cell)
    CALL ewald_pw_grid_update(ewald_pw, ewald_env, cell%hmat)

    natoms = SIZE(particle_set)
    nlocal_particles = 0
    nparticle_kind = SIZE(atomic_kind_set)
    DO ikind = 1, nparticle_kind
       nlocal_particles = nlocal_particles + local_particles%n_el(ikind)
    ENDDO

    ALLOCATE(f_nonbond(3,natoms), STAT=stat)
    CPASSERT(stat==0)

    nshell = 0
    IF(shell_present)THEN
       CALL fist_env_get(fist_env, shell_particle_set=shell_particle_set,  &
            core_particle_set=core_particle_set)
       CPASSERT(ASSOCIATED(shell_particle_set))
       CPASSERT(ASSOCIATED(core_particle_set))
       nshell = SIZE(shell_particle_set)
       ALLOCATE(fshell_nonbond(3, nshell), STAT = stat)
       CPASSERT(stat==0)
       ALLOCATE(fcore_nonbond(3, nshell), STAT = stat)
       CPASSERT(stat==0)
    ELSE
       NULLIFY(shell_particle_set,core_particle_set)
    END IF

    IF(fist_nonbond_env%do_nonbonded)THEN
       ! First check with list_control to update neighbor lists
       IF (ASSOCIATED(exclusions)) THEN
          CALL list_control(atomic_kind_set, particle_set, local_particles, &
               cell, fist_nonbond_env,para_env, mm_section, shell_particle_set,&
               core_particle_set, exclusions=exclusions)
       ELSE
          CALL list_control(atomic_kind_set, particle_set, local_particles, &
               cell, fist_nonbond_env,para_env, mm_section, shell_particle_set,&
               core_particle_set)
       END IF
    END IF

    ! Initialize force, energy and pressure tensor arrays
    DO i = 1, natoms
       particle_set(i)% f(1)= 0.0_dp
       particle_set(i)% f(2)= 0.0_dp
       particle_set(i)% f(3)= 0.0_dp
    ENDDO
    IF(nshell>0)THEN
       DO i=1,nshell
          shell_particle_set(i)% f(1)= 0.0_dp
          shell_particle_set(i)% f(2)= 0.0_dp
          shell_particle_set(i)% f(3)= 0.0_dp
          core_particle_set(i)% f(1)= 0.0_dp
          core_particle_set(i)% f(2)= 0.0_dp
          core_particle_set(i)% f(3)= 0.0_dp
       END DO
    ENDIF

    IF(use_virial)THEN
       pv_bc            = 0.0_dp
       pv_bond          = 0.0_dp
       pv_bend          = 0.0_dp
       pv_torsion       = 0.0_dp
       pv_imptors       = 0.0_dp
       pv_opbend        = 0.0_dp
       pv_urey_bradley  = 0.0_dp
       pv_nonbond       = 0.0_dp
       pv_g             = 0.0_dp
       virial%pv_virial = 0.0_dp
    END IF

    pot_nonbond       = 0.0_dp
    pot_manybody      = 0.0_dp
    pot_bond          = 0.0_dp
    pot_bend          = 0.0_dp
    pot_torsion       = 0.0_dp
    pot_opbend        = 0.0_dp
    pot_imptors       = 0.0_dp
    pot_urey_bradley  = 0.0_dp
    pot_shell         = 0.0_dp
    vg_coulomb        = 0.0_dp
    thermo%pot        = 0.0_dp
    thermo%harm_shell = 0.0_dp

    ! Get real-space non-bonded forces:
    IF(iw>0)THEN
       WRITE(iw,'(A)')" FIST:: FORCES IN INPUT..."
       WRITE(iw,'(3f15.9)')((particle_set(i)% f(j),j=1,3), i = 1, SIZE(particle_set))
    END IF

    IF(fist_nonbond_env%do_nonbonded)THEN
       ! Compute density for EAM
       CALL density_nonbond(fist_nonbond_env, particle_set, cell, para_env)

       ! Compute embedding function and manybody energy
       CALL energy_manybody(fist_nonbond_env, atomic_kind_set, local_particles, particle_set,&
            cell, pot_manybody, para_env , mm_section)

       ! Nonbond contribution + Manybody Forces
       IF(shell_present)THEN
          CALL force_nonbond(fist_nonbond_env, ewald_env, particle_set, cell, &
               pot_nonbond, f_nonbond, pv_nonbond, &
               fshell_nonbond=fshell_nonbond, fcore_nonbond=fcore_nonbond, &
               atprop_env=atprop_env,&
               atomic_kind_set=atomic_kind_set, use_virial=use_virial)
       ELSE
          CALL force_nonbond(fist_nonbond_env, ewald_env, particle_set, cell, &
               pot_nonbond, f_nonbond, pv_nonbond, atprop_env=atprop_env,&
               atomic_kind_set=atomic_kind_set, use_virial=use_virial)
          CALL force_nonbond_manybody(fist_nonbond_env, particle_set, cell, f_nonbond, pv_nonbond,&
               use_virial=use_virial)
       END IF
    ELSE
       f_nonbond  = 0.0_dp
       IF(shell_present)THEN
          fshell_nonbond = 0.0_dp
          fcore_nonbond  = 0.0_dp
       END IF
       IF(use_virial)THEN
          pv_nonbond = 0.0_dp
       END IF
    END IF

    IF(iw>0)THEN
       WRITE(iw,'(A)')" FIST:: NONBOND + R-SPACE ELECTROSTATIC FORCES ..."
       WRITE(iw,'(3f15.9)')f_nonbond
       IF(shell_present .AND. shell_model_ad)THEN
          WRITE(iw,'(A)')" FIST:: SHELL NONBOND + R-SPACE ELECTROSTATIC FORCES ..."
          WRITE(iw,'(3f15.9)')fshell_nonbond
          WRITE(iw,'(A)')" FIST:: CORE NONBOND + R-SPACE ELECTROSTATIC FORCES ..."
          WRITE(iw,'(3f15.9)')fcore_nonbond
       END IF
    END IF

    ! Get g-space non-bonded forces:
    IF(ewald_type /= do_ewald_none)THEN
       ! Determine size of the forces array
       SELECT CASE(ewald_type)
       CASE(do_ewald_ewald)
          fg_coulomb_size = nlocal_particles
       CASE DEFAULT
          fg_coulomb_size = natoms
       END SELECT
       ! Allocate and zeroing arrays
       ALLOCATE(fg_coulomb(3, fg_coulomb_size), STAT=stat)
       CPASSERT(stat==0)
       fg_coulomb = 0.0_dp
       IF(shell_present)THEN
          ALLOCATE(fgshell_coulomb(3, nshell), STAT=stat)
          CPASSERT(stat==0)
          ALLOCATE(fgcore_coulomb(3, nshell), STAT=stat)
          CPASSERT(stat==0)
          fgshell_coulomb = 0.0_dp
          fgcore_coulomb  = 0.0_dp
       END IF
       IF(shell_present.AND.do_multipoles)THEN
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="Multipoles and Core-Shell model not implemented.",&
               error_level=cp_failure_level)
       END IF
       ! If not multipole: Compute self-interaction and neutralizing background
       ! for multipoles is handled separately..
       IF(.NOT.do_multipoles)THEN
          CALL ewald_self(ewald_env, cell, atomic_kind_set, local_particles, &
               thermo%e_self, thermo%e_neut, fist_nonbond_env%charges)
          IF(atprop_env%energy) THEN
            CALL ewald_self_atom(ewald_env, atomic_kind_set, local_particles, &
               atprop_env%atener, fist_nonbond_env%charges)
            atprop_env%atener = atprop_env%atener + thermo%e_neut/SIZE(atprop_env%atener)
          END IF
       END IF

       ! Polarizable force-field
       IF(do_ipol /= do_fist_pol_none)THEN
          ! Check if an array of chagres was provided and in case abort due to lack of implementation
          IF (ASSOCIATED(fist_nonbond_env%charges)) THEN
             CALL cp_unimplemented_error(fromWhere=routineP, &
                  message="Polarizable force field and array charges not implemented!", &
                  error_level=cp_fatal_level)
          END IF
          ! Converge the dipoles self-consistently
          CALL fist_pol_evaluate(atomic_kind_set, multipoles, ewald_env, &
               ewald_pw, fist_nonbond_env, cell, particle_set, &
               local_particles, thermo, vg_coulomb, pot_nonbond, f_nonbond, &
               fg_coulomb, use_virial, pv_g, pv_nonbond, mm_section, do_ipol)
       ELSE
          ! Non-Polarizable force-field
          SELECT CASE(ewald_type)
          CASE(do_ewald_ewald)
             ! Parallelisation over atoms --> allocate local atoms
             IF(shell_present)THEN
                ! Check if an array of chagres was provided and in case abort due to lack of implementation
                IF (ASSOCIATED(fist_nonbond_env%charges)) THEN
                   CALL cp_unimplemented_error(fromWhere=routineP, &
                        message="Core-Shell and array charges not implemented!", &
                        error_level=cp_fatal_level)
                END IF
                IF(do_multipoles)THEN
                   CALL cp_unimplemented_error(fromWhere=routineP, &
                        message="Multipole Ewald and CORE-SHELL not yet implemented within Ewald sum!", &
                        error_level=cp_fatal_level)
                ELSE
                   CALL cp_unimplemented_error(fromWhere=routineP, &
                        message="Core-Shell model not yet implemented within Ewald sums.", &
                        error_level=cp_failure_level)
                END IF
             ELSE
                IF(do_multipoles)THEN
                   ! Check if an array of chagres was provided and in case abort due to lack of implementation
                   IF (ASSOCIATED(fist_nonbond_env%charges)) THEN
                      CALL cp_unimplemented_error(fromWhere=routineP, &
                           message="Multipole Ewald and array charges not implemented!", &
                           error_level=cp_fatal_level)
                   END IF
                   CALL ewald_multipole_evaluate(ewald_env, ewald_pw, fist_nonbond_env, cell, &
                        particle_set, local_particles, vg_coulomb, pot_nonbond, thermo%e_neut,&
                        thermo%e_self, multipoles%task, do_correction_bonded=.TRUE., do_forces=.TRUE.,&
                        do_stress=use_virial, do_efield=.FALSE., radii=multipoles%radii,&
                        charges=multipoles%charges, dipoles=multipoles%dipoles,&
                        quadrupoles=multipoles%quadrupoles, forces_local=fg_coulomb,&
                        forces_glob=f_nonbond, pv_local=pv_g, pv_glob=pv_nonbond,iw=iw2,&
                        do_debug=.TRUE.,atomic_kind_set=atomic_kind_set, mm_section=mm_section)
                ELSE
                   IF(atprop_env%energy) THEN
                      ALLOCATE(e_coulomb(fg_coulomb_size), STAT=stat)
                      CPASSERT(stat==0)
                   END IF
                   IF(atprop_env%stress) THEN
                      ALLOCATE(pv_coulomb(3,3,fg_coulomb_size), STAT=stat)
                      CPASSERT(stat==0)
                   END IF
                   CALL ewald_evaluate(ewald_env, ewald_pw, cell, atomic_kind_set, particle_set,&
                        local_particles, fg_coulomb, vg_coulomb, pv_g ,use_virial=use_virial,&
                        charges=fist_nonbond_env%charges, e_coulomb=e_coulomb, &
                        pv_coulomb=pv_coulomb)
                END IF
             END IF
          CASE(do_ewald_pme)
             ! Parallelisation over grids --> allocate all atoms
             IF(shell_present)THEN
                ! Check if an array of chagres was provided and in case abort due to lack of implementation
                IF (ASSOCIATED(fist_nonbond_env%charges)) THEN
                   CALL cp_unimplemented_error(fromWhere=routineP, &
                        message="Core-Shell and array charges not implemented!", &
                        error_level=cp_fatal_level)
                END IF
                IF(do_multipoles)THEN
                   CALL cp_unimplemented_error(fromWhere=routineP, &
                        message="Multipole Ewald and CORE-SHELL not yet implemented within a PME scheme!", &
                        error_level=cp_fatal_level)
                ELSE
                   CALL pme_evaluate(ewald_env, ewald_pw, cell, particle_set, vg_coulomb, fg_coulomb,&
                        pv_g,shell_particle_set=shell_particle_set, core_particle_set=core_particle_set,&
                        fgshell_coulomb=fgshell_coulomb, fgcore_coulomb=fgcore_coulomb, use_virial=use_virial,&
                        atprop=atprop_env)
                   CALL mp_sum(fgshell_coulomb, para_env%group)
                   CALL mp_sum(fgcore_coulomb, para_env%group)
                END IF
             ELSE
                IF(do_multipoles)THEN
                   CALL cp_unimplemented_error(fromWhere=routineP, &
                        message="Multipole Ewald not yet implemented within a PME scheme!", &
                        error_level=cp_fatal_level)
                ELSE
                   CALL pme_evaluate(ewald_env, ewald_pw, cell, particle_set, vg_coulomb, fg_coulomb,&
                        pv_g, use_virial=use_virial, charges=fist_nonbond_env%charges, &
                        atprop=atprop_env)
                END IF
             END IF
             CALL mp_sum(fg_coulomb, para_env%group)
          CASE(do_ewald_spme)
             ! Parallelisation over grids --> allocate all atoms
             IF(shell_present)THEN
                ! Check if an array of chagres was provided and in case abort due to lack of implementation
                IF (ASSOCIATED(fist_nonbond_env%charges)) THEN
                   CALL cp_unimplemented_error(fromWhere=routineP, &
                        message="Core-Shell and array charges not implemented!", &
                        error_level=cp_fatal_level)
                END IF
                IF(do_multipoles)THEN
                   CALL cp_unimplemented_error(fromWhere=routineP, &
                        message="Multipole Ewald and CORE-SHELL not yet implemented within a SPME scheme!", &
                        error_level=cp_fatal_level)
                ELSE
                   CALL spme_evaluate(ewald_env, ewald_pw, cell, particle_set, fg_coulomb, vg_coulomb,&
                        pv_g,shell_particle_set=shell_particle_set, core_particle_set=core_particle_set,&
                        fgshell_coulomb=fgshell_coulomb, fgcore_coulomb=fgcore_coulomb,use_virial=use_virial,&
                        atprop=atprop_env)
                   CALL mp_sum(fgshell_coulomb, para_env%group)
                   CALL mp_sum(fgcore_coulomb, para_env%group)
                END IF
             ELSE
                IF(do_multipoles)THEN
                   CALL cp_unimplemented_error(fromWhere=routineP, &
                        message="Multipole Ewald not yet implemented within a SPME scheme!", &
                        error_level=cp_fatal_level)
                ELSE
                   CALL spme_evaluate(ewald_env, ewald_pw, cell, particle_set, fg_coulomb, vg_coulomb,&
                        pv_g, use_virial=use_virial, charges=fist_nonbond_env%charges, &
                        atprop=atprop_env)
                END IF
             END IF
             CALL mp_sum(fg_coulomb, para_env%group)
          END SELECT
       END IF

       ! Subtract the interaction between screening charges. This is a
       ! correction in real-space and depends on the neighborlists. Therefore
       ! it is only executed if fist_nonbond_env%do_nonbonded is set.
       IF(fist_nonbond_env%do_nonbonded)THEN
          ! Correction for core-shell model
          IF(shell_present)THEN
             CALL bonded_correct_gaussian(fist_nonbond_env, atomic_kind_set, &
                  local_particles, particle_set, ewald_env, thermo%e_bonded, &
                  pv_bc, shell_particle_set=shell_particle_set, &
                  core_particle_set=core_particle_set, atprop_env=atprop_env, cell=cell, &
                  use_virial=use_virial)
          ELSE
             IF(.NOT. do_multipoles)THEN
                CALL bonded_correct_gaussian(fist_nonbond_env, &
                     atomic_kind_set, local_particles, particle_set, &
                     ewald_env, thermo%e_bonded, pv_bc=pv_bc, atprop_env=atprop_env, cell=cell, &
                     use_virial=use_virial)
             END IF
          END IF
       END IF

       IF(.NOT. do_multipoles)THEN
          ! Multipole code has its own printing routine.
          CALL ewald_print(iw2, pot_nonbond, vg_coulomb, thermo%e_self, thermo%e_neut,&
               thermo%e_bonded)
       END IF
    ELSE
       IF(use_virial)THEN
          pv_g  = 0.0_dp
          pv_bc = 0.0_dp
       END IF
       thermo%e_neut = 0.0_dp
    END IF

    IF(iw>0)THEN
       IF(ALLOCATED(fg_coulomb))THEN
          WRITE(iw,'(A)')" FIST:: NONBONDED ELECTROSTATIC FORCES IN G-SPACE..."
          WRITE(iw,'(3f15.9)')((fg_coulomb(j,i),j=1,3),i = 1, SIZE(fg_coulomb,2))
       END IF
       IF(shell_present .AND. shell_model_ad .AND. ALLOCATED(fgshell_coulomb))THEN
          WRITE(iw,'(A)')" FIST:: SHELL NONBONDED ELECTROSTATIC FORCES IN G-SPACE..."
          WRITE(iw,'(3f15.9)')((fgshell_coulomb(j,i),j=1,3),i = 1, SIZE(fgshell_coulomb,2))
          WRITE(iw,'(A)')" FIST:: CORE NONBONDED ELECTROSTATIC FORCES IN G-SPACE..."
          WRITE(iw,'(3f15.9)')((fgcore_coulomb(j,i),j=1,3),i = 1, SIZE(fgcore_coulomb,2))
       END IF
    END IF

    ! Get intramolecular forces
    IF(PRESENT(debug))THEN
       CALL force_intra_control(molecule_set, molecule_kind_set, &
            local_molecules, particle_set, shell_particle_set, &
            core_particle_set, pot_bond, pot_bend, pot_urey_bradley, &
            pot_torsion, pot_imptors, pot_opbend, pot_shell, pv_bond, pv_bend, &
            pv_urey_bradley, pv_torsion, pv_imptors, pv_opbend, &
            debug%f_bond, debug%f_bend, debug%f_torsion, debug%f_ub, &
            debug%f_imptors, debug%f_opbend, cell, use_virial, atprop_env)

    ELSE
       CALL force_intra_control(molecule_set, molecule_kind_set, &
            local_molecules, particle_set, shell_particle_set, &
            core_particle_set, pot_bond, pot_bend, pot_urey_bradley, &
            pot_torsion, pot_imptors, pot_opbend, pot_shell, pv_bond, pv_bend, &
            pv_urey_bradley, pv_torsion, pv_imptors, pv_opbend, &
            cell=cell, use_virial=use_virial, atprop_env=atprop_env)
    ENDIF

    IF(iw>0)THEN
       xdum1 = cp_unit_from_cp2k(pot_bond,"kcalmol")
       xdum2 = cp_unit_from_cp2k(pot_bend,"kcalmol")
       xdum3 = cp_unit_from_cp2k(pot_urey_bradley,"kcalmol")
       WRITE(iw,'(A)')" FIST energy contributions in kcal/mol:"
       WRITE(iw,'(1x,"BOND    = ",f13.4,'//&
                 '2x,"ANGLE   = ",f13.4,'//&
                 '2x,"UBRAD   = ",f13.4)')xdum1, xdum2, xdum3
       xdum1 = cp_unit_from_cp2k(pot_torsion,"kcalmol")
       xdum2 = cp_unit_from_cp2k(pot_imptors,"kcalmol")
       xdum3 = cp_unit_from_cp2k(pot_opbend,"kcalmol")
       WRITE(iw,'(1x,"TORSION = ",f13.4,'//&
                 '2x,"IMPTORS = ",f13.4,'//&
                 '2x,"OPBEND  = ",f13.4)')xdum1, xdum2, xdum3

       WRITE(iw,'(A)')" FIST:: CORRECTED BONDED ELECTROSTATIC FORCES + INTERNAL FORCES..."
       WRITE(iw,'(3f15.9)')((particle_set(i)%f(j),j=1,3), i = 1, SIZE(particle_set))
       IF(shell_present .AND. shell_model_ad)THEN
          WRITE(iw,'(A)')" FIST:: SHELL CORRECTED BONDED ELECTROSTATIC FORCES + INTERNAL FORCES..."
          WRITE(iw,'(3f15.9)')((shell_particle_set(i)%f(j),j=1,3), i = 1, SIZE(shell_particle_set))
          WRITE(iw,'(A)')" FIST:: CORE CORRECTED BONDED ELECTROSTATIC FORCES + INTERNAL FORCES..."
          WRITE(iw,'(3f15.9)')((core_particle_set(i)%f(j),j=1,3), i = 1, SIZE(core_particle_set))
       END IF
    END IF

    ! add up all the potential energies
    thermo%pot = pot_nonbond + pot_bond + pot_bend + pot_torsion + pot_opbend + &
         pot_imptors + pot_urey_bradley + pot_manybody + pot_shell

    CALL mp_sum(thermo%pot, para_env%group)
    IF(shell_present)THEN
       thermo%harm_shell = pot_shell
       CALL mp_sum(thermo%harm_shell, para_env%group)
    END IF
    ! add g-space contributions if needed
    IF(ewald_type /= do_ewald_none)THEN
       ! e_self, e_neut, and ebonded are already summed over all processors
       ! vg_coulomb is not calculated in parallel
       thermo%e_gspace = vg_coulomb
       thermo%pot = thermo%pot + thermo%e_self + thermo%e_neut
       thermo%pot = thermo%pot + vg_coulomb + thermo%e_bonded
       ! add the induction energy of the dipoles for polarizable models
       IF(do_ipol /= do_fist_pol_none)thermo%pot = thermo%pot + thermo%e_induction
    END IF

    ! add up all the forces
    !
    ! nonbonded forces might be calculated for atoms not on this node
    ! ewald forces are strictly local -> sum only over pnode
    ! We first sum the forces in f_nonbond, this allows for a more efficient
    ! global sum in the parallel code and in the end copy them back to part
    ALLOCATE(f_total(3,natoms), STAT=stat)
    CPASSERT(stat==0)
    f_total = 0.0_dp
    DO i = 1, natoms
       f_total(1, i)= particle_set(i)%f(1)+ f_nonbond(1, i)
       f_total(2, i)= particle_set(i)%f(2)+ f_nonbond(2, i)
       f_total(3, i)= particle_set(i)%f(3)+ f_nonbond(3, i)
    END DO
    IF(shell_present)THEN
       ALLOCATE(fshell_total(3,nshell), STAT=stat)
       CPASSERT(stat==0)
       ALLOCATE(fcore_total(3,nshell), STAT=stat)
       CPASSERT(stat==0)
       DO i = 1, nshell
          fshell_total(1, i)= shell_particle_set(i)%f(1)+ fshell_nonbond(1, i)
          fshell_total(2, i)= shell_particle_set(i)%f(2)+ fshell_nonbond(2, i)
          fshell_total(3, i)= shell_particle_set(i)%f(3)+ fshell_nonbond(3, i)
          fcore_total(1, i)= core_particle_set(i)%f(1)+ fcore_nonbond(1, i)
          fcore_total(2, i)= core_particle_set(i)%f(2)+ fcore_nonbond(2, i)
          fcore_total(3, i)= core_particle_set(i)%f(3)+ fcore_nonbond(3, i)
       END DO
    END IF

    IF(iw>0)THEN
       WRITE(iw,'(A)')" FIST::(1)INTERNAL + ELECTROSTATIC BONDED + NONBONDED"
       WRITE(iw,'(3f15.9)')((f_total(j, i),j=1,3), i = 1, natoms)
       IF(shell_present .AND. shell_model_ad)THEN
          WRITE(iw,'(A)')" FIST::(1)SHELL INTERNAL + ELECTROSTATIC BONDED + NONBONDED"
          WRITE(iw,'(3f15.9)')((fshell_total(j, i),j=1,3), i = 1, nshell)
          WRITE(iw,'(A)')" FIST::(1)CORE INTERNAL + ELECTROSTATIC BONDED + NONBONDED"
          WRITE(iw,'(3f15.9)')((fcore_total(j, i),j=1,3), i = 1, nshell)
       END IF
    END IF

    ! Adding in the reciprocal forces: EWALD is a special case because of distrubted data
    IF(ewald_type == do_ewald_ewald)THEN
       node = 0
       DO iparticle_kind = 1, nparticle_kind
          nparticle_local = local_particles%n_el(iparticle_kind)
          DO iparticle_local = 1, nparticle_local
             ii = local_particles%list(iparticle_kind)%array(iparticle_local)
             node = node + 1

             f_total(1, ii)= f_total(1, ii)+ fg_coulomb(1, node)
             f_total(2, ii)= f_total(2, ii)+ fg_coulomb(2, node)
             f_total(3, ii)= f_total(3, ii)+ fg_coulomb(3, node)
             IF(PRESENT(debug))THEN
                debug%f_g(1, ii)=  debug%f_g(1, ii)+ fg_coulomb(1, node)
                debug%f_g(2, ii)=  debug%f_g(2, ii)+ fg_coulomb(2, node)
                debug%f_g(3, ii)=  debug%f_g(3, ii)+ fg_coulomb(3, node)
             ENDIF
             IF(atprop_env%energy) THEN
                atprop_env%atener(ii) = atprop_env%atener(ii) + e_coulomb(node)
             END IF
             IF(atprop_env%stress) THEN
                atprop_env%atstress(:,:,ii) = atprop_env%atstress(:,:,ii) + pv_coulomb(:,:,node)
             END IF
          END DO
       END DO
       IF(atprop_env%energy) THEN
         DEALLOCATE(e_coulomb, STAT = stat)
         CPASSERT(stat==0)
       END IF
       IF(atprop_env%stress) THEN
         DEALLOCATE(pv_coulomb, STAT = stat)
         CPASSERT(stat==0)
       END IF
    END IF

    IF(iw>0)THEN
       WRITE(iw,'(A)')" FIST::(2)TOTAL FORCES(1)+ ELECTROSTATIC FORCES"
       WRITE(iw,'(3f15.9)')((f_total(j, i),j=1,3), i = 1, natoms)
       IF(shell_present .AND. shell_model_ad)THEN
          WRITE(iw,'(A)')" FIST::(2)SHELL TOTAL FORCES(1)+ ELECTROSTATIC FORCES "
          WRITE(iw,'(3f15.9)')((fshell_total(j, i),j=1,3), i = 1, nshell)
          WRITE(iw,'(A)')" FIST::(2)CORE TOTAL FORCES(1)+ ELECTROSTATIC FORCES"
          WRITE(iw,'(3f15.9)')((fcore_total(j, i),j=1,3), i = 1, nshell)
       END IF
    END IF

    IF(use_virial)THEN
       ! Add up all the pressure tensors
       IF(ewald_type == do_ewald_none)THEN
          virial%pv_virial = pv_nonbond + pv_bond + pv_bend + pv_torsion + pv_imptors + pv_urey_bradley
          CALL mp_sum(virial%pv_virial, para_env%group)
       ELSE
          ident = 0.0_dp
          DO i = 1, 3
             ident(i, i)= 1.0_dp
          END DO

          virial%pv_virial = pv_nonbond + pv_bond + pv_bend + pv_torsion + pv_imptors + pv_urey_bradley + pv_bc
          CALL mp_sum(virial%pv_virial, para_env%group)

          virial%pv_virial = virial%pv_virial + ident * thermo%e_neut
          virial%pv_virial = virial%pv_virial + pv_g
       END IF
    END IF

    ! Sum total forces
    CALL mp_sum(f_total, para_env%group)
    IF(shell_present .AND. shell_model_ad)THEN
       CALL mp_sum(fshell_total, para_env%group)
       CALL mp_sum(fcore_total, para_env%group)
    END IF

    ! Assign to particle_set
    SELECT CASE(ewald_type)
    CASE(do_ewald_spme, do_ewald_pme)
       IF(shell_present .AND. shell_model_ad)THEN
          DO i = 1, natoms
             IF(particle_set(i)%shell_index/=0)THEN
                atomic_kind => particle_set(i)%atomic_kind
                CALL get_atomic_kind(atomic_kind=atomic_kind, shell=shell, mass=mass)
                umass = 1.0_dp/mass
                massc = shell%mass_core*umass
                masss = shell%mass_shell*umass
                shell_index = particle_set(i)%shell_index
                fshell_total(1,shell_index)= fshell_total(1,shell_index)+ particle_set(i)%f(1)*masss
                fshell_total(2,shell_index)= fshell_total(2,shell_index)+ particle_set(i)%f(2)*masss
                fshell_total(3,shell_index)= fshell_total(3,shell_index)+ particle_set(i)%f(3)*masss
                fcore_total(1,shell_index) = fcore_total(1,shell_index)+ particle_set(i)%f(1)*massc
                fcore_total(2,shell_index) = fcore_total(2,shell_index)+ particle_set(i)%f(2)*massc
                fcore_total(3,shell_index) = fcore_total(3,shell_index)+ particle_set(i)%f(3)*massc
             ELSE
                particle_set(i)%f(1)= f_total(1, i)+ fg_coulomb(1, i)
                particle_set(i)%f(2)= f_total(2, i)+ fg_coulomb(2, i)
                particle_set(i)%f(3)= f_total(3, i)+ fg_coulomb(3, i)
             END IF
          END DO

          DO i = 1,nshell
             shell_particle_set(i)%f(1)= fshell_total(1,i)+ fgshell_coulomb(1,i)
             shell_particle_set(i)%f(2)= fshell_total(2,i)+ fgshell_coulomb(2,i)
             shell_particle_set(i)%f(3)= fshell_total(3,i)+ fgshell_coulomb(3,i)
             core_particle_set(i)%f(1) = fcore_total(1,i) + fgcore_coulomb(1,i)
             core_particle_set(i)%f(2) = fcore_total(2,i) + fgcore_coulomb(2,i)
             core_particle_set(i)%f(3) = fcore_total(3,i) + fgcore_coulomb(3,i)
          END DO

       ELSEIF(shell_present .AND. .NOT. shell_model_ad)THEN
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="Non adiabatic shell-model not implemented.",&
               error_level=cp_failure_level)
       ELSE
          DO i = 1, natoms
             particle_set(i)%f(1)= f_total(1, i)+ fg_coulomb(1, i)
             particle_set(i)%f(2)= f_total(2, i)+ fg_coulomb(2, i)
             particle_set(i)%f(3)= f_total(3, i)+ fg_coulomb(3, i)
          END DO
       END IF
    CASE(do_ewald_ewald, do_ewald_none)
       IF(shell_present .AND. shell_model_ad)THEN
          DO i = 1, natoms
             IF(particle_set(i)%shell_index/=0)THEN
                atomic_kind => particle_set(i)%atomic_kind
                CALL get_atomic_kind(atomic_kind=atomic_kind, shell=shell, mass=mass)
                umass = 1.0_dp/mass
                massc = shell%mass_core*umass
                masss = shell%mass_shell*umass
                shell_index = particle_set(i)%shell_index
                fshell_total(1,shell_index)= fshell_total(1,shell_index)+ particle_set(i)%f(1)*masss
                fshell_total(2,shell_index)= fshell_total(2,shell_index)+ particle_set(i)%f(2)*masss
                fshell_total(3,shell_index)= fshell_total(3,shell_index)+ particle_set(i)%f(3)*masss
                fcore_total(1,shell_index) = fcore_total(1,shell_index) + particle_set(i)%f(1)*massc
                fcore_total(2,shell_index) = fcore_total(2,shell_index) + particle_set(i)%f(2)*massc
                fcore_total(3,shell_index) = fcore_total(3,shell_index) + particle_set(i)%f(3)*massc
             ELSE
                particle_set(i)%f(1)= f_total(1, i)
                particle_set(i)%f(2)= f_total(2, i)
                particle_set(i)%f(3)= f_total(3, i)
             END IF
          END DO
          DO i = 1,nshell
             shell_particle_set(i)%f(1)= fshell_total(1,i)
             shell_particle_set(i)%f(2)= fshell_total(2,i)
             shell_particle_set(i)%f(3)= fshell_total(3,i)
             core_particle_set(i)%f(1)= fcore_total(1,i)
             core_particle_set(i)%f(2)= fcore_total(2,i)
             core_particle_set(i)%f(3)= fcore_total(3,i)
          END DO
       ELSEIF(shell_present .AND. .NOT. shell_model_ad)THEN
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="Non adiabatic shell-model not implemented.",&
               error_level=cp_failure_level)
       ELSE
          DO i = 1, natoms
             particle_set(i)%f(1)= f_total(1, i)
             particle_set(i)%f(2)= f_total(2, i)
             particle_set(i)%f(3)= f_total(3, i)
          END DO
       END IF
    END SELECT

    IF(iw>0)THEN
       WRITE(iw,'(A)')" FIST::(3)TOTAL FORCES - THE END..."
       WRITE(iw,'(3f15.9)')((particle_set(i)%f(j),j=1,3), i = 1, natoms)
       IF(shell_present .AND. shell_model_ad)THEN
          WRITE(iw,'(A)')" FIST::(3)SHELL TOTAL FORCES - THE END..."
          WRITE(iw,'(3f15.9)')((shell_particle_set(i)% f(j),j=1,3), i = 1, nshell)
          WRITE(iw,'(A)')" FIST::(3)CORE TOTAL FORCES - THE END..."
          WRITE(iw,'(3f15.9)')((core_particle_set(i)% f(j),j=1,3), i = 1, nshell)
       END IF
       WRITE(iw,'(A,f15.9)')"Energy after FIST calculation.. exiting now ::",thermo%pot
    END IF
    !
    ! If we are doing debugging, check if variables are present and assign
    !
    IF(PRESENT(debug))THEN
       CALL mp_sum(pot_manybody, para_env%group)
       debug%pot_manybody = pot_manybody
       CALL mp_sum(pot_nonbond, para_env%group)
       debug%pot_nonbond = pot_nonbond
       CALL mp_sum(pot_bond, para_env%group)
       debug%pot_bond = pot_bond
       CALL mp_sum(pot_bend, para_env%group)
       debug%pot_bend = pot_bend
       CALL mp_sum(pot_torsion, para_env%group)
       debug%pot_torsion = pot_torsion
       CALL mp_sum(pot_imptors, para_env%group)
       debug%pot_imptors = pot_imptors
       CALL mp_sum(pot_opbend, para_env%group)
       debug%pot_opbend = pot_opbend
       CALL mp_sum(pot_urey_bradley, para_env%group)
       debug%pot_urey_bradley = pot_urey_bradley
       CALL mp_sum(f_nonbond, para_env%group)
       debug%f_nonbond = f_nonbond
       IF(use_virial)THEN
          CALL mp_sum(pv_nonbond, para_env%group)
          debug%pv_nonbond = pv_nonbond
          CALL mp_sum(pv_bond, para_env%group)
          debug%pv_bond = pv_bond
          CALL mp_sum(pv_bend, para_env%group)
          debug%pv_bend = pv_bend
          CALL mp_sum(pv_torsion, para_env%group)
          debug%pv_torsion = pv_torsion
          CALL mp_sum(pv_imptors, para_env%group)
          debug%pv_imptors = pv_imptors
          CALL mp_sum(pv_urey_bradley, para_env%group)
          debug%pv_ub = pv_urey_bradley
       END IF
       SELECT CASE(ewald_type)
       CASE(do_ewald_spme, do_ewald_pme)
          debug%pot_g = vg_coulomb
          debug%pv_g = pv_g
          debug%f_g = fg_coulomb
       CASE(do_ewald_ewald)
          debug%pot_g = vg_coulomb
          IF(use_virial)debug%pv_g = pv_g
       CASE default
          debug%pot_g = 0.0_dp
          debug%f_g = 0.0_dp
          IF(use_virial)debug%pv_g = 0.0_dp
       END SELECT
    END IF

    ! print properties if requested
    print_section => section_vals_get_subs_vals(mm_section,"PRINT")
    CALL print_fist(fist_env, print_section, atomic_kind_set, particle_set, cell)

    ! deallocating all local variables
    IF(ALLOCATED(fg_coulomb))THEN
       DEALLOCATE(fg_coulomb, STAT = stat)
       CPASSERT(stat==0)
    END IF
    IF(ALLOCATED(f_total))THEN
       DEALLOCATE(f_total, STAT = stat)
       CPASSERT(stat==0)
    END IF
    DEALLOCATE(f_nonbond, STAT = stat)
    CPASSERT(stat==0)
    IF(shell_present)THEN
       DEALLOCATE(fshell_total, STAT = stat)
       CPASSERT(stat==0)
       IF(ewald_type /= do_ewald_none)THEN
          DEALLOCATE(fgshell_coulomb, STAT = stat)
          CPASSERT(stat==0)
       END IF
       DEALLOCATE(fshell_nonbond, STAT = stat)
       CPASSERT(stat==0)
    END IF
    CALL cp_print_key_finished_output(iw,logger,mm_section,"PRINT%DERIVATIVES")
    CALL cp_print_key_finished_output(iw2,logger,mm_section,"PRINT%EWALD_INFO")
    CALL timestop(handle)

END SUBROUTINE fist_calc_energy_force

! *****************************************************************************
!> \brief Print properties number according the requests in input file
!> \param fist_env ...
!> \param print_section ...
!> \param atomic_kind_set ...
!> \param particle_set ...
!> \param cell ...
!> \par History
!>      [01.2006] created
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE print_fist(fist_env, print_section, atomic_kind_set, particle_set, cell)
    TYPE(fist_environment_type), POINTER     :: fist_env
    TYPE(section_vals_type), POINTER         :: print_section
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cell_type), POINTER                 :: cell

    INTEGER                                  :: unit_nr
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
    TYPE(section_vals_type), POINTER         :: print_key

    NULLIFY(logger,print_key,fist_nonbond_env)
    logger => cp_get_default_logger()
    print_key => section_vals_get_subs_vals(print_section,"dipole")
    CALL fist_env_get(fist_env, fist_nonbond_env=fist_nonbond_env)
    IF(BTEST(cp_print_key_should_output(logger%iter_info,print_key),&
         cp_p_file))THEN
       unit_nr=cp_print_key_unit_nr(logger,print_section,"dipole",&
            extension=".data",middle_name="MM_DIPOLE",log_filename=.FALSE.)
       CALL fist_dipole(fist_env, print_section, atomic_kind_set, particle_set,&
            cell, unit_nr, fist_nonbond_env%charges)
       CALL cp_print_key_finished_output(unit_nr,logger,print_key)
    END IF

  END SUBROUTINE print_fist

! *****************************************************************************
!> \brief Evaluates the Dipole of a classical charge distribution(point-like)
!>      possibly using the berry phase formalism
!> \param fist_env ...
!> \param print_section ...
!> \param atomic_kind_set ...
!> \param particle_set ...
!> \param cell ...
!> \param unit_nr ...
!> \param charges ...
!> \par History
!>      [01.2006] created
!>      [12.2007] tlaino - University of Zurich - debug and extended
!> \author Teodoro Laino
! *****************************************************************************
  SUBROUTINE fist_dipole(fist_env, print_section, atomic_kind_set, particle_set,&
             cell, unit_nr, charges)
    TYPE(fist_environment_type), POINTER     :: fist_env
    TYPE(section_vals_type), POINTER         :: print_section
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(cell_type), POINTER                 :: cell
    INTEGER, INTENT(IN)                      :: unit_nr
    REAL(KIND=dp), DIMENSION(:), OPTIONAL, &
      POINTER                                :: charges

    CHARACTER(LEN=default_string_length)     :: description, dipole_type
    COMPLEX(KIND=dp)                         :: dzeta, dzphase(3), zeta, &
                                                zphase(3)
    COMPLEX(KIND=dp), DIMENSION(3)           :: dggamma, ggamma
    INTEGER                                  :: i, iparticle_kind, j, &
                                                reference
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: do_berry
    REAL(KIND=dp) :: charge_tot, ci(3), dci(3), dipole(3), dipole_deriv(3), &
      drcc(3), dria(3), dtheta, gvec(3), q, rcc(3), ria(3), theta, tmp(3), &
      via(3)
    REAL(KIND=dp), DIMENSION(:), POINTER     :: ref_point
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_result_type), POINTER            :: results

    NULLIFY(atomic_kind)
    ! Reference point
    reference = section_get_ival(print_section,keyword_name="DIPOLE%REFERENCE")
    NULLIFY(ref_point)
    description='[DIPOLE]'
    CALL section_vals_val_get(print_section,"DIPOLE%REF_POINT",r_vals=ref_point)
    CALL section_vals_val_get(print_section,"DIPOLE%PERIODIC",l_val=do_berry)

    CALL get_reference_point(rcc,drcc,fist_env=fist_env,reference=reference,ref_point=ref_point)

    ! Dipole deriv will be the derivative of the Dipole(dM/dt=\sum e_j v_j)
    dipole_deriv = 0.0_dp
    dipole       = 0.0_dp
    IF(do_berry)THEN
       dipole_type="[BERRY PHASE]"
       rcc = pbc(rcc,cell)
       charge_tot = 0._dp
       IF (ASSOCIATED(charges)) THEN
          charge_tot=SUM(charges)
       ELSE
          DO i = 1,SIZE(particle_set)
             atomic_kind => particle_set(i)%atomic_kind
             CALL get_atomic_kind(atomic_kind=atomic_kind,qeff=q)
             charge_tot = charge_tot + q
          END DO
       END IF
       ria  = twopi * MATMUL(cell%h_inv, rcc)
       zphase  = CMPLX(COS(ria), SIN(ria), dp)**charge_tot

       dria = twopi * MATMUL(cell%h_inv, drcc)
       dzphase = charge_tot * CMPLX(-SIN(ria), COS(ria), dp)**(charge_tot-1.0_dp)* dria

       ggamma  = CMPLX(1.0_dp,0.0_dp,KIND=dp)
       dggamma = CMPLX(0.0_dp,0.0_dp,KIND=dp)
       DO iparticle_kind = 1, SIZE(atomic_kind_set)
          atomic_kind => atomic_kind_set(iparticle_kind)
          CALL get_atomic_kind(atomic_kind=atomic_kind, qeff=q, atom_list=atom_list)

          DO i = 1, SIZE(atom_list)
             ria = particle_set(atom_list(i))%r(:)
             ria = pbc(ria,cell)
             via = particle_set(atom_list(i))%v(:)
             IF (ASSOCIATED(charges)) q=charges(atom_list(i))
             DO j = 1, 3
                gvec  = twopi*cell%h_inv(j,:)
                theta = SUM(ria(:)*gvec(:))
                dtheta= SUM(via(:)*gvec(:))
                zeta  = CMPLX(COS(theta),SIN(theta),KIND=dp)**(-q)
                dzeta = -q * CMPLX(-SIN(theta),COS(theta),KIND=dp)**(-q-1.0_dp)* dtheta
                dggamma(j)= dggamma(j)* zeta + ggamma(j)* dzeta
                ggamma(j) = ggamma(j)* zeta
             END DO
          ENDDO
       END DO
       dggamma = dggamma * zphase + ggamma * dzphase
       ggamma = ggamma * zphase
       IF(ALL(REAL(ggamma,KIND=dp)/= 0.0_dp))THEN
          tmp = AIMAG(ggamma)/REAL(ggamma,KIND=dp)
          ci  = ATAN(tmp)
          dci =(1.0_dp /(1.0_dp + tmp**2))* &
              (AIMAG(dggamma)*REAL(ggamma,KIND=dp)-AIMAG(ggamma)*REAL(dggamma,KIND=dp))/(REAL(ggamma,KIND=dp))**2

          dipole       = MATMUL(cell%hmat, ci)/twopi
          dipole_deriv = MATMUL(cell%hmat, dci)/twopi
       END IF
       CALL fist_env_get(fist_env=fist_env,results=results)
       CALL cp_results_erase(results,description)
       CALL put_results(results,description,dipole)
    ELSE
       dipole_type="[Non Periodic]"
       DO i = 1, SIZE(particle_set)
          atomic_kind => particle_set(i)%atomic_kind
          ria = particle_set(i)%r(:) ! no pbc(particle_set(i)%r(:),cell) so that the total dipole is the sum of the molecular dipoles
          CALL get_atomic_kind(atomic_kind=atomic_kind,qeff=q)
          IF (ASSOCIATED(charges)) q=charges(i)
          dipole = dipole - q *(ria-rcc)
          dipole_deriv(:)= dipole_deriv(:)- q *(particle_set(i)%v(:)- drcc)
       END DO
       CALL fist_env_get(fist_env=fist_env,results=results)
       CALL cp_results_erase(results,description)
       CALL put_results(results,description,dipole)
    END IF
    IF(unit_nr>0)THEN
       WRITE(unit_nr,'(1X,A,T48,3F11.6)')"MM DIPOLE "//TRIM(dipole_type)//"(A.U.)|",dipole
       WRITE(unit_nr,'(1X,A,T48,3F11.6)')"MM DIPOLE "//TRIM(dipole_type)//"(Debye)|",dipole*debye
       WRITE(unit_nr,'(1X,A,T48,3F11.6)')"MM DIPOLE "//TRIM(dipole_type)//" DERIVATIVE(A.U.)|",dipole_deriv
    END IF

  END SUBROUTINE fist_dipole

END MODULE fist_force
