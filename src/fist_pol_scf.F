!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2009  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      CJM APRIL-30-2009: only uses fist_env
!>      Teodoro Laino [tlaino] - 05.2009 : Generalization to different Ewald
!>                                         methods (initial framework)
!> \author CJM 
! *****************************************************************************

MODULE fist_pol_scf
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind
  USE cell_types,                      ONLY: cell_type
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE ewald_environment_types,         ONLY: ewald_env_get,&
                                             ewald_environment_type
  USE ewald_pw_types,                  ONLY: ewald_pw_type
  USE ewalds_multipole,                ONLY: ewald_multipole_evaluate
  USE fist_energy_types,               ONLY: fist_energy_type
  USE fist_nonbond_env_types,          ONLY: fist_nonbond_env_type
  USE input_constants,                 ONLY: do_ewald_ewald,&
                                             do_ewald_pme,&
                                             do_ewald_spme
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp
  USE multipole_types,                 ONLY: multipole_type
  USE particle_types,                  ONLY: particle_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'fist_pol_scf'

  PUBLIC :: fist_pol_evaluate

CONTAINS

! *****************************************************************************
!> \brief Main driver for evaluating energy/forces in a polarizable forcefield
! *****************************************************************************
  SUBROUTINE fist_pol_evaluate (atomic_kind_set, multipoles, ewald_env, ewald_pw,&
       fist_nonbond_env, cell, particle_set, local_particles, thermo, vg_coulomb,&
       pot_nonbond, f_nonbond, fg_coulomb,use_virial, pv_g, pv_nonbond, mm_section,&
       error)

    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(multipole_type), POINTER            :: multipoles
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
    TYPE(cell_type), POINTER                 :: cell
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(fist_energy_type), POINTER          :: thermo
    REAL(KIND=dp)                            :: vg_coulomb, pot_nonbond
    REAL(KIND=dp), DIMENSION(:, :)           :: f_nonbond, fg_coulomb
    LOGICAL, INTENT(IN)                      :: use_virial
    REAL(KIND=dp), DIMENSION(3, 3)           :: pv_g, pv_nonbond
    TYPE(section_vals_type), POINTER         :: mm_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fist_pol_evaluate', &
      routineP = moduleN//':'//routineN

    INTEGER :: ewald_type, handle, i, iatom, ii, ikind, iter, iw, iw2, &
      max_ipol_iter, natom_of_kind, natoms, nkind, ntot, stat
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: failure, iwarn
    REAL(KIND=dp)                            :: a_pol, diff, eps_pol, &
                                                pot_nonbond_local
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: efield1, local_dipoles
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)
    failure = .FALSE.
    NULLIFY(logger,atomic_kind)
    logger => cp_error_get_logger(error)

    iw = cp_print_key_unit_nr(logger,mm_section,"PRINT%ITER_INFO",&
         extension=".mmLog",error=error)
    iw2= cp_print_key_unit_nr(logger,mm_section,"PRINT%EWALD_INFO",&
         extension=".mmLog",error=error)

    CALL ewald_env_get (ewald_env, max_ipol_iter=max_ipol_iter, eps_pol=eps_pol,&
         ewald_type=ewald_type, error=error)

    natoms = SIZE ( particle_set )
    ALLOCATE (efield1(3,natoms), STAT=stat )
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (local_dipoles(3,natoms), STAT=stat )
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    pot_nonbond_local=0.0_dp
    local_dipoles    =0.0_dp

    ! Evaluate the electrostatic with Ewald schemes
    CALL eval_pol_ewald(ewald_type, ewald_env, ewald_pw, fist_nonbond_env, cell, &
         particle_set, local_particles, vg_coulomb, pot_nonbond_local, thermo,&
         multipoles, do_correction_bonded=.TRUE., do_forces=.FALSE., &
         do_stress=.FALSE., do_efield=.TRUE., iw2=iw2, do_debug=.FALSE.,&
         atomic_kind_set=atomic_kind_set, mm_section=mm_section, efield1=efield1,&
         error=error)

    nkind = SIZE(atomic_kind_set)
    ntot = 0
    IF (iw>0) WRITE(iw,FMT='(/,T14,"Iteration",9X,"Conv.",T61,"Electrostatic Energy")')
    pol_scf: DO iter=1, max_ipol_iter
       diff=0.0_dp
       pot_nonbond_local=0.0_dp      ! dummmy variable to avoid duplication over iterations
       thermo % e_induction = 0.0_dp ! always zero (otherwise compiled over iterations)
       DO ikind = 1, nkind
          atomic_kind => atomic_kind_set ( ikind )
          CALL get_atomic_kind ( atomic_kind, a_pol = a_pol, natom=natom_of_kind, atom_list=atom_list )
          IF ( a_pol == 0 ) CYCLE
          ! finding the number of atoms with polarizable sites
          ntot = ntot + natom_of_kind
          DO iatom = 1, natom_of_kind
             ii = atom_list(iatom)
             local_dipoles(:,ii)=a_pol*efield1(:,ii)
             DO i=1, 3
                diff=diff+(local_dipoles(i,ii)-multipoles%dipoles(i,ii))**2
             END DO
             multipoles%dipoles(:,ii)=local_dipoles(:,ii)
             ! Compute induction term while we are here
             thermo%e_induction = thermo%e_induction +  &
                  DOT_PRODUCT(multipoles%dipoles( :, ii ), multipoles%dipoles( :, ii ))/a_pol/2.0_dp 

          END DO
       END DO
       diff = SQRT (diff/REAL(ntot,KIND=dp))
       IF (iw>0) THEN
          WRITE(iw,FMT='(T5,"POL_SCF|",5X,I5,5X,E12.6,T61,F20.10)')iter,diff,vg_coulomb+pot_nonbond_local
       END IF
       IF ( diff <= eps_pol ) THEN 
          IF (iw>0) WRITE(iw,FMT='(T5,"POL_SCF|",1X,"Self-consistent Polarization achieved.",/)')
          EXIT pol_scf
       END IF

       iwarn = ((diff >eps_pol) .AND. (iter==max_ipol_iter))
       IF (iwarn.AND.iw>0) WRITE(iw,FMT='(T5,"POL_SCF|",1X,"Self-consistent Polarization not converged!",/)')
       CALL cp_assert(.NOT.iwarn,cp_warning_level,cp_assertion_failed,routineP,&
            "Self-consistent Polarization not converged! "//&
CPSourceFileRef,&
            only_ionode=.TRUE.)

       pot_nonbond_local=0.0_dp
       CALL eval_pol_ewald(ewald_type, ewald_env, ewald_pw, fist_nonbond_env, cell, &
            particle_set, local_particles, vg_coulomb, pot_nonbond_local, thermo,&
            multipoles, do_correction_bonded=.TRUE., do_forces=.FALSE., &
            do_stress=.FALSE., do_efield=.TRUE., iw2=iw2, do_debug=.FALSE.,&
            atomic_kind_set=atomic_kind_set, mm_section=mm_section,&
            efield1=efield1, error=error)
    END DO pol_scf

    ! Now evaluate after convergence to obtain forces and converged energies
    CALL eval_pol_ewald(ewald_type, ewald_env, ewald_pw, fist_nonbond_env, cell,&
         particle_set, local_particles, vg_coulomb, pot_nonbond, thermo,&
         multipoles, do_correction_bonded=.TRUE., do_forces=.TRUE., &
         do_stress=use_virial, do_efield=.FALSE., iw2=iw2, do_debug=.FALSE.,&
         atomic_kind_set=atomic_kind_set, mm_section=mm_section,&
         forces_local=fg_coulomb, forces_glob=f_nonbond, &
         pv_local=pv_g, pv_glob=pv_nonbond, error=error)

    ! Deallocate working arrays
    DEALLOCATE(efield1,stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(local_dipoles,stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_print_key_finished_output(iw2,logger,mm_section,&
         "PRINT%EWALD_INFO",error=error)
    CALL cp_print_key_finished_output(iw,logger,mm_section,&
         "PRINT%ITER_INFO",error=error)

    CALL timestop(handle)
  END SUBROUTINE fist_pol_evaluate

! *****************************************************************************
!> \brief Main driver for evaluating electrostatic in polarible forcefields
!>        All the dependence on the Ewald method should go here!
!> \author Teodoro Laino [tlaino] 05.2009
! *****************************************************************************
  SUBROUTINE eval_pol_ewald(ewald_type, ewald_env, ewald_pw, fist_nonbond_env,&
       cell, particle_set, local_particles, vg_coulomb, pot_nonbond, thermo,&
       multipoles, do_correction_bonded, do_forces, do_stress, do_efield, iw2,&
       do_debug, atomic_kind_set, mm_section, efield0, efield1, forces_local,&
       forces_glob, pv_local, pv_glob, error)

    INTEGER, INTENT(IN)                      :: ewald_type
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
    TYPE(cell_type), POINTER                 :: cell
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(distribution_1d_type), POINTER      :: local_particles
    REAL(KIND=dp)                            :: vg_coulomb, pot_nonbond
    TYPE(fist_energy_type), POINTER          :: thermo
    TYPE(multipole_type), POINTER            :: multipoles
    LOGICAL, INTENT(IN)                      :: do_correction_bonded, &
                                                do_forces, do_stress, &
                                                do_efield
    INTEGER, INTENT(IN)                      :: iw2
    LOGICAL, INTENT(IN)                      :: do_debug
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(section_vals_type), POINTER         :: mm_section
    REAL(KIND=dp), DIMENSION(:), OPTIONAL    :: efield0
    REAL(KIND=dp), DIMENSION(:, :), OPTIONAL :: efield1
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(INOUT), OPTIONAL                :: forces_local, forces_glob, &
                                                pv_local, pv_glob
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'eval_pol_ewald', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure

    CALL timeset(routineN,handle)
    failure = .FALSE.

    SELECT CASE(ewald_type)
    CASE(do_ewald_ewald)
       CALL ewald_multipole_evaluate(ewald_env, ewald_pw, fist_nonbond_env, cell, &
            particle_set, local_particles, vg_coulomb, pot_nonbond, thermo%e_neut,&
            thermo%e_self, multipoles%task, do_correction_bonded=do_correction_bonded,&
            do_forces=do_forces,do_stress=do_stress, do_efield=do_efield,&
            charges=multipoles%charges, dipoles=multipoles%dipoles,&
            quadrupoles=multipoles%quadrupoles, forces_local=forces_local,&
            forces_glob=forces_glob, pv_local=pv_local, pv_glob=pv_glob, iw=iw2,&
            do_debug=do_debug, atomic_kind_set=atomic_kind_set, mm_section=mm_section,&
            efield0=efield0, efield1=efield1, error=error)
    CASE(do_ewald_pme)
       CALL cp_unimplemented_error(fromWhere=routineP, &
            message="Multipole Ewald not yet implemented within a PME scheme!", &
            error=error, error_level=cp_fatal_level)
    CASE(do_ewald_spme)
       CALL cp_unimplemented_error(fromWhere=routineP, &
            message="Multipole Ewald not yet implemented within a SPME scheme!", &
            error=error, error_level=cp_fatal_level)
    END SELECT
    CALL timestop(handle)
  END SUBROUTINE eval_pol_ewald

END MODULE fist_pol_scf
