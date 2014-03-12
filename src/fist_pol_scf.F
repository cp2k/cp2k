!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
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
  USE input_constants,                 ONLY: do_fist_pol_cg,&
                                             do_fist_pol_sc
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE multipole_types,                 ONLY: multipole_type
  USE particle_types,                  ONLY: particle_type
  USE pw_poisson_types,                ONLY: do_ewald_ewald,&
                                             do_ewald_pme,&
                                             do_ewald_spme
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE
  LOGICAL, PRIVATE                     :: debug_this_module=.FALSE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'fist_pol_scf'

  PUBLIC :: fist_pol_evaluate

CONTAINS

! *****************************************************************************
!> \brief Main driver for evaluating energy/forces in a polarizable forcefield
!> \author Toon.Verstraelen@gmail.com (2010-03-01)
! *****************************************************************************
  SUBROUTINE fist_pol_evaluate (atomic_kind_set, multipoles, ewald_env, &
       ewald_pw, fist_nonbond_env, cell, particle_set, local_particles, &
       thermo, vg_coulomb, pot_nonbond, f_nonbond, fg_coulomb,use_virial, &
       pv_g, pv_nonbond, mm_section, do_ipol, error)

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
    INTEGER                                  :: do_ipol
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'fist_pol_evaluate', &
      routineP = moduleN//':'//routineN

    SELECT CASE (do_ipol)
    CASE (do_fist_pol_sc)
      CALL fist_pol_evaluate_sc(atomic_kind_set, multipoles, ewald_env, &
        ewald_pw, fist_nonbond_env, cell, particle_set, local_particles, &
        thermo, vg_coulomb, pot_nonbond, f_nonbond, fg_coulomb,use_virial, &
        pv_g, pv_nonbond, mm_section, error)
    CASE (do_fist_pol_cg)
      CALL fist_pol_evaluate_cg(atomic_kind_set, multipoles, ewald_env, &
        ewald_pw, fist_nonbond_env, cell, particle_set, local_particles, &
        thermo, vg_coulomb, pot_nonbond, f_nonbond, fg_coulomb,use_virial, &
        pv_g, pv_nonbond, mm_section, error)
    END SELECT

  END SUBROUTINE fist_pol_evaluate

! *****************************************************************************
!> \brief Self-consistent solver for a polarizable force-field
!> \author Toon.Verstraelen@gmail.com (2010-03-01)
!> \note
!>    Method: Given an initial guess of the induced dipoles, the electrostatic
!>    field is computed at each dipole. Then new induced dipoles are computed
!>    following p = alpha x E. This is repeated until a convergence criterion is
!>    met. The convergence is measured as the RSMD of the derivatives of the
!>    electrostatic energy (including dipole self-energy) towards the components
!>    of the dipoles.
! *****************************************************************************
  SUBROUTINE fist_pol_evaluate_sc (atomic_kind_set, multipoles, ewald_env, ewald_pw,&
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

    CHARACTER(len=*), PARAMETER :: routineN = 'fist_pol_evaluate_sc', &
      routineP = moduleN//':'//routineN

    INTEGER :: ewald_type, handle, i, iatom, ii, ikind, iter, iw, iw2, j, &
      max_ipol_iter, natom_of_kind, natoms, nkind, ntot, stat
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: failure, iwarn
    REAL(KIND=dp)                            :: apol, cpol, eps_pol, &
                                                pot_nonbond_local, rmsd, &
                                                tmp_trace
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: efield1, efield2
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
    ALLOCATE (efield2(9,natoms), STAT=stat )
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    nkind = SIZE(atomic_kind_set)
    IF (iw>0) WRITE(iw,FMT='(/,T5,"POL_SCF|","Method: self-consistent")')
    IF (iw>0) WRITE(iw,FMT='(T5,"POL_SCF|","  Iteration",7X,"Conv.",T49,"Electrostatic & Induction Energy")')
    pol_scf: DO iter=1, max_ipol_iter
       ! Evaluate the electrostatic with Ewald schemes
       CALL eval_pol_ewald(ewald_type, ewald_env, ewald_pw, fist_nonbond_env, cell, &
            particle_set, local_particles, vg_coulomb, pot_nonbond_local, thermo,&
            multipoles, do_correction_bonded=.TRUE., do_forces=.FALSE., &
            do_stress=.FALSE., do_efield=.TRUE., iw2=iw2, do_debug=.FALSE.,&
            atomic_kind_set=atomic_kind_set, mm_section=mm_section,&
            efield1=efield1, efield2=efield2, error=error)
       CALL mp_sum(pot_nonbond_local,logger%para_env%group)

       ! compute the new dipoles, qudrupoles, and check for convergence
       ntot = 0
       rmsd = 0.0_dp
       thermo%e_induction = 0.0_dp
       DO ikind = 1, nkind
         atomic_kind => atomic_kind_set(ikind)
         CALL get_atomic_kind(atomic_kind, apol=apol, cpol=cpol, natom=natom_of_kind, atom_list=atom_list)
          ! ignore atoms with dipole and quadrupole polarizability zero
          IF ( apol == 0 .AND. cpol == 0 ) CYCLE
          ! increment counter correctly
          IF ( apol /= 0 ) ntot = ntot + natom_of_kind
          IF ( cpol /= 0 ) ntot = ntot + natom_of_kind

          DO iatom = 1, natom_of_kind
             ii = atom_list(iatom)
             IF ( apol /= 0 ) THEN
               DO i=1, 3
                ! the rmsd of the derivatives of the energy towards the
                ! components of the atomic dipole moments
                 rmsd = rmsd + (multipoles%dipoles(i,ii)/apol - efield1(i,ii))**2
               END DO
             END IF
             IF ( cpol /= 0 ) THEN
                rmsd = rmsd + (multipoles%quadrupoles(1,1,ii)/cpol - efield2(1,ii))**2
                rmsd = rmsd + (multipoles%quadrupoles(2,1,ii)/cpol - efield2(2,ii))**2
                rmsd = rmsd + (multipoles%quadrupoles(3,1,ii)/cpol - efield2(3,ii))**2
                rmsd = rmsd + (multipoles%quadrupoles(1,2,ii)/cpol - efield2(4,ii))**2
                rmsd = rmsd + (multipoles%quadrupoles(2,2,ii)/cpol - efield2(5,ii))**2
                rmsd = rmsd + (multipoles%quadrupoles(3,2,ii)/cpol - efield2(6,ii))**2
                rmsd = rmsd + (multipoles%quadrupoles(1,3,ii)/cpol - efield2(7,ii))**2
                rmsd = rmsd + (multipoles%quadrupoles(2,3,ii)/cpol - efield2(8,ii))**2
                rmsd = rmsd + (multipoles%quadrupoles(3,3,ii)/cpol - efield2(9,ii))**2
             END IF
! compute dipole
             multipoles%dipoles(:,ii) = apol*efield1(:,ii)
! compute quadrupole
             IF ( cpol /= 0 ) THEN
               multipoles%quadrupoles(1,1,ii) = cpol*efield2(1,ii)
               multipoles%quadrupoles(2,1,ii) = cpol*efield2(2,ii)
               multipoles%quadrupoles(3,1,ii) = cpol*efield2(3,ii)
               multipoles%quadrupoles(1,2,ii) = cpol*efield2(4,ii)
               multipoles%quadrupoles(2,2,ii) = cpol*efield2(5,ii)
               multipoles%quadrupoles(3,2,ii) = cpol*efield2(6,ii)
               multipoles%quadrupoles(1,3,ii) = cpol*efield2(7,ii)
               multipoles%quadrupoles(2,3,ii) = cpol*efield2(8,ii)
               multipoles%quadrupoles(3,3,ii) = cpol*efield2(9,ii)
             END IF
        ! Compute the new induction term while we are here
             IF ( apol /= 0 ) THEN
               thermo%e_induction = thermo%e_induction + &
                    DOT_PRODUCT(multipoles%dipoles(:,ii), &
                                multipoles%dipoles(:,ii))/apol/2.0_dp
             END IF
             IF ( cpol /= 0 ) THEN
               tmp_trace = 0._dp
               DO i = 1, 3
                 DO j = 1, 3
                    tmp_trace = tmp_trace  + &
                    multipoles%quadrupoles ( i, j, ii ) * multipoles%quadrupoles ( i, j, ii )
                 END DO
               END DO
               thermo%e_induction = thermo%e_induction + tmp_trace/cpol/6.0_dp
             END IF
          END DO
       END DO
       rmsd = SQRT (rmsd/REAL(ntot,KIND=dp))
       IF (iw>0) THEN
          ! print the energy that is minimized (this is electrostatic + induction)
          WRITE(iw,FMT='(T5,"POL_SCF|",5X,I5,5X,E12.6,T61,F20.10)') iter, &
            rmsd, vg_coulomb + pot_nonbond_local + thermo%e_induction
       END IF
       IF ( rmsd <= eps_pol ) THEN
          IF (iw>0) WRITE(iw,FMT='(T5,"POL_SCF|",1X,"Self-consistent Polarization achieved.")')
          EXIT pol_scf
       END IF

       iwarn = ((rmsd > eps_pol) .AND. (iter == max_ipol_iter))
       IF (iwarn.AND.iw>0) WRITE(iw,FMT='(T5,"POL_SCF|",1X,"Self-consistent Polarization not converged!")')
       CALL cp_assert(.NOT.iwarn,cp_warning_level,cp_assertion_failed,routineP,&
            "Self-consistent Polarization not converged! "//&
            CPSourceFileRef,&
            only_ionode=.TRUE.)
    END DO pol_scf

    ! Now evaluate after convergence to obtain forces and converged energies
    CALL eval_pol_ewald(ewald_type, ewald_env, ewald_pw, fist_nonbond_env, cell,&
         particle_set, local_particles, vg_coulomb, pot_nonbond_local, thermo,&
         multipoles, do_correction_bonded=.TRUE., do_forces=.TRUE., &
         do_stress=use_virial, do_efield=.FALSE., iw2=iw2, do_debug=.FALSE.,&
         atomic_kind_set=atomic_kind_set, mm_section=mm_section,&
         forces_local=fg_coulomb, forces_glob=f_nonbond, &
         pv_local=pv_g, pv_glob=pv_nonbond, error=error)
    pot_nonbond=pot_nonbond+pot_nonbond_local
    CALL mp_sum(pot_nonbond_local,logger%para_env%group)

    IF (iw>0) THEN
      ! print the energy that is minimized (this is electrostatic + induction)
      WRITE(iw,FMT='(T5,"POL_SCF|",5X,"Final",T61,F20.10,/)') &
        vg_coulomb + pot_nonbond_local + thermo%e_induction
    END IF

    ! Deallocate working arrays
    DEALLOCATE(efield1,stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(efield2,stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_print_key_finished_output(iw2,logger,mm_section,&
         "PRINT%EWALD_INFO",error=error)
    CALL cp_print_key_finished_output(iw,logger,mm_section,&
         "PRINT%ITER_INFO",error=error)

    CALL timestop(handle)
  END SUBROUTINE fist_pol_evaluate_sc


! *****************************************************************************
!> \brief Conjugate-gradient solver for a polarizable force-field
!> \author Toon.Verstraelen@gmail.com (2010-03-01)
!> \note
!>     Method: The dipoles are found by minimizing the sum of the electrostatic
!>     and the induction energy directly using a conjugate gradient method. This
!>     routine assumes that the energy is a quadratic function of the dipoles.
!>     Finding the minimum is then done by solving a linear system. This will
!>     not work for polarizable force fields that include hyperpolarizability.
!>
!>     The implementation of the conjugate gradient solver for linear systems
!>     is described in chapter 2.7 Sparse Linear Systems, under the section
!>     "Conjugate Gradient Method for a Sparse System". Although the inducible
!>     dipoles are the solution of a dense linear system, the same algorithm is
!>     still recommended for this situation. One does not have access to the
!>     entire hardness kernel to compute the solution with conventional linear
!>     algebra routines, but one only has a function that computes the dot
!>     product of the hardness kernel and a vector. (This is the routine that
!>     computes the electrostatic field at the atoms for a given vector of
!>     inducible dipoles.) Given such function, the conjugate gradient method
!>     is an efficient way to compute the solution of a linear system, and it
!>     scales well towards many degrees of freedom in terms of efficiency and
!>     memory usage.
! *****************************************************************************
  SUBROUTINE fist_pol_evaluate_cg (atomic_kind_set, multipoles, ewald_env, ewald_pw,&
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

    CHARACTER(len=*), PARAMETER :: routineN = 'fist_pol_evaluate_cg', &
      routineP = moduleN//':'//routineN

    INTEGER :: ewald_type, handle, i, iatom, ii, ikind, iter, iw, iw2, &
      max_ipol_iter, natom_of_kind, natoms, nkind, ntot, stat
    INTEGER, DIMENSION(:), POINTER           :: atom_list
    LOGICAL                                  :: failure, iwarn
    REAL(KIND=dp)                            :: alpha, apol, beta, denom, &
                                                eps_pol, pot_nonbond_local, &
                                                rmsd
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: conjugate, conjugate_applied, &
                                                efield1, efield1_ext, &
                                                residual, tmp_dipoles
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

    ! allocate work arrays
    natoms = SIZE ( particle_set )
    ALLOCATE (efield1(3,natoms), STAT=stat )
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (tmp_dipoles(3,natoms), STAT=stat )
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (residual(3,natoms), STAT=stat )
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (conjugate(3,natoms), STAT=stat )
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (conjugate_applied(3,natoms), STAT=stat )
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (efield1_ext(3,natoms), STAT=stat )
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    ! Compute the 'external' electrostatic field (all inducible dipoles
    ! equal to zero). this is required for the conjugate gradient solver.
    ! We assume that all dipoles are inducible dipoles.
    tmp_dipoles = multipoles%dipoles ! backup of dipoles
    multipoles%dipoles = 0.0_dp
    CALL eval_pol_ewald(ewald_type, ewald_env, ewald_pw, fist_nonbond_env, cell, &
         particle_set, local_particles, vg_coulomb, pot_nonbond_local, thermo, &
         multipoles, do_correction_bonded=.TRUE., do_forces=.FALSE., &
         do_stress=.FALSE., do_efield=.TRUE., iw2=iw2, do_debug=.FALSE., &
         atomic_kind_set=atomic_kind_set, mm_section=mm_section, &
         efield1=efield1_ext, error=error)
    multipoles%dipoles = tmp_dipoles ! restore backup

    ! Compute the electric field with the initial guess of the dipoles.
    CALL eval_pol_ewald(ewald_type, ewald_env, ewald_pw, fist_nonbond_env, cell, &
        particle_set, local_particles, vg_coulomb, pot_nonbond_local, thermo, &
        multipoles, do_correction_bonded=.TRUE., do_forces=.FALSE., &
        do_stress=.FALSE., do_efield=.TRUE., iw2=iw2, do_debug=.FALSE., &
        atomic_kind_set=atomic_kind_set, mm_section=mm_section, &
        efield1=efield1, error=error)

    ! Compute the first residual explicitly.
    nkind = SIZE(atomic_kind_set)
    ntot = 0
    residual = 0.0_dp
    DO ikind = 1, nkind
       atomic_kind => atomic_kind_set(ikind)
       CALL get_atomic_kind(atomic_kind, apol=apol, natom=natom_of_kind, atom_list=atom_list)
       ! ignore atoms with polarizability zero
       IF ( apol == 0 ) CYCLE
       ntot = ntot + natom_of_kind
       DO iatom = 1, natom_of_kind
          ii = atom_list(iatom)
          DO i=1, 3
             ! residual = b - A x
             residual(i,ii) = efield1(i,ii) - multipoles%dipoles(i,ii)/apol
          END DO
       END DO
    END DO
    ! The first conjugate residual is equal to the residual.
    conjugate = residual

    IF (iw>0) WRITE(iw,FMT='(/,T5,"POL_SCF|","Method: conjugate-gradient")')
    IF (iw>0) WRITE(iw,FMT='(T5,"POL_SCF|","  Iteration",7X,"Conv.",T49,"Electrostatic & Induction Energy")')
    pol_scf: DO iter=1, max_ipol_iter
       IF (debug_this_module) THEN
          ! In principle the residual must not be computed explicitly any more. It
          ! is obtained in an indirect way below. When the DEBUG flag is set, the
          ! explicit residual is computed and compared with the implicitly derived
          ! residual as a double check.
          CALL eval_pol_ewald(ewald_type, ewald_env, ewald_pw, fist_nonbond_env, cell, &
               particle_set, local_particles, vg_coulomb, pot_nonbond_local, thermo, &
               multipoles, do_correction_bonded=.TRUE., do_forces=.FALSE., &
               do_stress=.FALSE., do_efield=.TRUE., iw2=iw2, do_debug=.FALSE., &
               atomic_kind_set=atomic_kind_set, mm_section=mm_section, &
               efield1=efield1, error=error)
          ! inapropriate use of denom to check the error on the residual
          denom = 0.0_dp
       END IF
       rmsd = 0.0_dp
       ! Compute the rmsd of the residual.
       DO ikind = 1, nkind
          atomic_kind => atomic_kind_set(ikind)
          CALL get_atomic_kind(atomic_kind, apol=apol, natom=natom_of_kind, atom_list=atom_list)
          ! ignore atoms with polarizability zero
          IF ( apol == 0 ) CYCLE
          DO iatom = 1, natom_of_kind
             ii = atom_list(iatom)
             DO i=1, 3
               ! residual = b - A x
               rmsd = rmsd + residual(i,ii)**2
               IF (debug_this_module) THEN
                  denom = denom + (residual(i,ii) - (efield1(i,ii) - &
                          multipoles%dipoles(i,ii)/apol))**2
               END IF
             END DO
          END DO
       END DO
       rmsd = SQRT(rmsd/ntot)
       IF (iw>0) THEN
          WRITE(iw,FMT='(T5,"POL_SCF|",5X,I5,5X,E12.6,T67,"(not computed)")') iter, rmsd
          IF (debug_this_module) THEN
             denom = SQRT(denom/ntot)
             WRITE(iw,FMT='(T5,"POL_SCF|",5X,"Error on implicit residual:",5X,E12.6)') denom
          END IF
       END IF

       ! Apply the hardness kernel to the conjugate residual.
       ! We assume that all non-zero dipoles are inducible dipoles.
       tmp_dipoles = multipoles%dipoles ! backup of dipoles
       multipoles%dipoles = conjugate
       CALL eval_pol_ewald(ewald_type, ewald_env, ewald_pw, fist_nonbond_env, cell, &
            particle_set, local_particles, vg_coulomb, pot_nonbond_local, thermo, &
            multipoles, do_correction_bonded=.TRUE., do_forces=.FALSE., &
            do_stress=.FALSE., do_efield=.TRUE., iw2=iw2, do_debug=.FALSE., &
            atomic_kind_set=atomic_kind_set, mm_section=mm_section, &
            efield1=conjugate_applied, error=error)
       multipoles%dipoles = tmp_dipoles ! restore backup
       conjugate_applied = efield1_ext - conjugate_applied

       ! Finish conjugate_applied and compute alpha from the conjugate gradient algorithm.
       alpha = 0.0_dp
       denom = 0.0_dp
       DO ikind = 1, nkind
          atomic_kind => atomic_kind_set(ikind)
          CALL get_atomic_kind(atomic_kind, apol=apol, natom=natom_of_kind, atom_list=atom_list)
          ! ignore atoms with polarizability zero
          IF ( apol == 0 ) CYCLE
          DO iatom = 1, natom_of_kind
             ii = atom_list(iatom)
             DO i=1, 3
               conjugate_applied(i,ii) = conjugate_applied(i,ii) + conjugate(i,ii)/apol
             END DO
             alpha = alpha + DOT_PRODUCT(residual(:,ii), residual(:,ii))
             denom = denom + DOT_PRODUCT(conjugate(:,ii), conjugate_applied(:,ii))
          END DO
       END DO
       alpha = alpha/denom

       ! Compute the new residual and beta from the conjugate gradient method.
       beta = 0.0_dp
       denom = 0.0_dp
       DO ikind = 1, nkind
          atomic_kind => atomic_kind_set(ikind)
          CALL get_atomic_kind(atomic_kind, apol=apol, natom=natom_of_kind, atom_list=atom_list)
          IF ( apol == 0 ) CYCLE
          DO iatom = 1, natom_of_kind
             ii = atom_list(iatom)
             denom = denom + DOT_PRODUCT(residual(:,ii), residual(:,ii))
             DO i=1, 3
               residual(i,ii) = residual(i,ii) - alpha*conjugate_applied(i,ii)
             END DO
             beta = beta + DOT_PRODUCT(residual(:,ii), residual(:,ii))
          END DO
       END DO
       beta = beta/denom

       ! Compute the new dipoles, the new conjugate residual, and the induction
       ! energy.
       thermo%e_induction = 0.0_dp
       DO ikind = 1, nkind
          atomic_kind => atomic_kind_set(ikind)
          CALL get_atomic_kind(atomic_kind, apol=apol, natom=natom_of_kind, atom_list=atom_list)
          ! ignore atoms with polarizability zero
          IF ( apol == 0 ) CYCLE
          DO iatom = 1, natom_of_kind
             ii = atom_list(iatom)
             DO i=1, 3
               multipoles%dipoles(i,ii) = multipoles%dipoles(i,ii) + alpha*conjugate(i,ii)
               conjugate(i,ii) = residual(i,ii) + beta*conjugate(i,ii)
               thermo%e_induction = thermo%e_induction + multipoles%dipoles(i,ii)**2/apol/2.0_dp
             END DO
          END DO
       END DO

       ! Quit if rmsd is low enough.
       IF ( rmsd <= eps_pol ) THEN
          IF (iw>0) WRITE(iw,FMT='(T5,"POL_SCF|",1X,"Self-consistent Polarization achieved.")')
          EXIT pol_scf
       END IF

       ! Print warning when not converged.
       iwarn = ((rmsd > eps_pol) .AND. (iter == max_ipol_iter))
       IF (iwarn.AND.iw>0) WRITE(iw,FMT='(T5,"POL_SCF|",1X,"Self-consistent Polarization not converged!")')
       CALL cp_assert(.NOT.iwarn,cp_warning_level,cp_assertion_failed,routineP,&
            "Self-consistent Polarization not converged! "//&
            CPSourceFileRef,&
            only_ionode=.TRUE.)
    END DO pol_scf

    IF (debug_this_module) THEN
       ! Now evaluate after convergence to obtain forces and converged energies
       CALL eval_pol_ewald(ewald_type, ewald_env, ewald_pw, fist_nonbond_env, cell,&
            particle_set, local_particles, vg_coulomb, pot_nonbond_local, thermo,&
            multipoles, do_correction_bonded=.TRUE., do_forces=.TRUE., &
            do_stress=use_virial, do_efield=.TRUE., iw2=iw2, do_debug=.FALSE.,&
            atomic_kind_set=atomic_kind_set, mm_section=mm_section,&
            forces_local=fg_coulomb, forces_glob=f_nonbond, &
            pv_local=pv_g, pv_glob=pv_nonbond, efield1=efield1, error=error)

       ! Do a final check on the convergence: compute the residual explicitely
       rmsd = 0.0_dp
       DO ikind = 1, nkind
          atomic_kind => atomic_kind_set(ikind)
          CALL get_atomic_kind(atomic_kind, apol=apol, natom=natom_of_kind, atom_list=atom_list)
          ! ignore atoms with polarizability zero
          IF ( apol == 0 ) CYCLE
          DO iatom = 1, natom_of_kind
             ii = atom_list(iatom)
             DO i=1, 3
                ! residual = b - A x
                rmsd = rmsd + (efield1(i,ii) - multipoles%dipoles(i,ii)/apol)**2
             END DO
          END DO
       END DO
       rmsd = SQRT(rmsd/ntot)
       IF (iw>0) WRITE(iw,FMT='(T5,"POL_SCF|",1X,"Final RMSD of residual:",5X,E12.6)') rmsd
       ! Stop program when congergence is not reached after all
       IF ( rmsd > eps_pol ) THEN
         CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
              "Error in the conjugate gradient method for self-consistent polarization! "//&
              CPSourceFileRef,&
              only_ionode=.TRUE.)
       END IF
    ELSE
       ! Now evaluate after convergence to obtain forces and converged energies
       CALL eval_pol_ewald(ewald_type, ewald_env, ewald_pw, fist_nonbond_env, cell,&
            particle_set, local_particles, vg_coulomb, pot_nonbond_local, thermo,&
            multipoles, do_correction_bonded=.TRUE., do_forces=.TRUE., &
            do_stress=use_virial, do_efield=.FALSE., iw2=iw2, do_debug=.FALSE.,&
            atomic_kind_set=atomic_kind_set, mm_section=mm_section,&
            forces_local=fg_coulomb, forces_glob=f_nonbond, &
            pv_local=pv_g, pv_glob=pv_nonbond, error=error)
    ENDIF
    pot_nonbond=pot_nonbond+pot_nonbond_local
    CALL mp_sum(pot_nonbond_local, logger%para_env%group)

    IF (iw>0) THEN
      WRITE(iw,FMT='(T5,"POL_SCF|",5X,"Final",T61,F20.10,/)') &
        vg_coulomb + pot_nonbond_local + thermo%e_induction
    END IF

    ! Deallocate working arrays
    DEALLOCATE(efield1,stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(tmp_dipoles,stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(residual,stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(conjugate,stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(conjugate_applied,stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE(efield1_ext,stat=stat)
    CPPrecondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_print_key_finished_output(iw2,logger,mm_section,&
         "PRINT%EWALD_INFO",error=error)
    CALL cp_print_key_finished_output(iw,logger,mm_section,&
         "PRINT%ITER_INFO",error=error)

    CALL timestop(handle)
  END SUBROUTINE fist_pol_evaluate_cg

! *****************************************************************************
!> \brief Main driver for evaluating electrostatic in polarible forcefields
!>        All the dependence on the Ewald method should go here!
!> \author Teodoro Laino [tlaino] 05.2009
! *****************************************************************************
  SUBROUTINE eval_pol_ewald(ewald_type, ewald_env, ewald_pw, fist_nonbond_env,&
       cell, particle_set, local_particles, vg_coulomb, pot_nonbond, thermo,&
       multipoles, do_correction_bonded, do_forces, do_stress, do_efield, iw2,&
       do_debug, atomic_kind_set, mm_section, efield0, efield1, efield2, forces_local,&
       forces_glob, pv_local, pv_glob, error)

    INTEGER, INTENT(IN)                      :: ewald_type
    TYPE(ewald_environment_type), POINTER    :: ewald_env
    TYPE(ewald_pw_type), POINTER             :: ewald_pw
    TYPE(fist_nonbond_env_type), POINTER     :: fist_nonbond_env
    TYPE(cell_type), POINTER                 :: cell
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(distribution_1d_type), POINTER      :: local_particles
    REAL(KIND=dp), INTENT(OUT)               :: vg_coulomb, pot_nonbond
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
    REAL(KIND=dp), DIMENSION(:, :), OPTIONAL :: efield1, efield2
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

    pot_nonbond = 0.0_dp ! Initialization..
    vg_coulomb  = 0.0_dp ! Initialization..
    SELECT CASE(ewald_type)
    CASE(do_ewald_ewald)
       CALL ewald_multipole_evaluate(ewald_env, ewald_pw, fist_nonbond_env, cell, &
            particle_set, local_particles, vg_coulomb, pot_nonbond, thermo%e_neut,&
            thermo%e_self, multipoles%task, do_correction_bonded=do_correction_bonded,&
            do_forces=do_forces,do_stress=do_stress, do_efield=do_efield,&
            radii=multipoles%radii, charges=multipoles%charges,&
            dipoles=multipoles%dipoles, quadrupoles=multipoles%quadrupoles,&
            forces_local=forces_local, forces_glob=forces_glob, pv_local=pv_local,&
            pv_glob=pv_glob, iw=iw2, do_debug=do_debug, atomic_kind_set=atomic_kind_set,&
            mm_section=mm_section, efield0=efield0, efield1=efield1, efield2=efield2, &
           error=error)
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
