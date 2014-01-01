!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Performs the Harris functional force correction
!> \par History
!>      06.2005 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
MODULE harris_force

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE cell_types,                      ONLY: cell_type
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_p_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE global_types,                    ONLY: global_environment_type
  USE harris_env_types,                ONLY: harris_env_get,&
                                             harris_env_type
  USE harris_force_types,              ONLY: harris_force_type
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_methods,                      ONLY: pw_axpy,&
                                             pw_copy,&
                                             pw_scale,&
                                             pw_transfer,&
                                             pw_zero
  USE pw_poisson_methods,              ONLY: pw_poisson_solve
  USE pw_poisson_types,                ONLY: pw_poisson_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_p_type,&
                                             pw_release
  USE qs_core_energies,                ONLY: calculate_ecore_overlap,&
                                             calculate_ecore_self
  USE qs_core_hamiltonian,             ONLY: build_core_hamiltonian_matrix
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_force_types,                  ONLY: deallocate_qs_force,&
                                             duplicate_qs_force,&
                                             qs_force_type,&
                                             zero_qs_force
  USE qs_integrate_potential,          ONLY: integrate_v_core_rspace,&
                                             integrate_v_rspace
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE qs_rho_methods,                  ONLY: diff_rho_type
  USE qs_rho_types,                    ONLY: qs_rho_get,&
                                             qs_rho_type
  USE qs_vxc,                          ONLY: qs_vxc_create
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc,                              ONLY: xc_calc_2nd_deriv,&
                                             xc_prep_2nd_deriv
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                             xc_dset_release
  USE xc_rho_set_types,                ONLY: xc_rho_set_release,&
                                             xc_rho_set_type
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'harris_force'

  ! *** Public subroutines ***
  PUBLIC :: harris_force_correction, &
            harris_force_EVal, &
            harris_calc_nsc_force

!***

CONTAINS

! *****************************************************************************
!> \brief Performs the harris functional force correction
!> \param qs_env The QS environment of matter
!> \param harris_force The harris force type in which all necessary informations
!>                      are stored
!> \param globenv The global environment from which the para_env is extracted
!> \param error variable to control error logging, stopping,...
!>               see module cp_error_handling
!> \note
!>      Performs the force calculation according to the Harris energy functional
!> \par History
!>      06.2005 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
  SUBROUTINE harris_force_correction(qs_env, harris_env, globenv, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(harris_env_type), POINTER           :: harris_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'harris_force_correction', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iatom, ikind, &
                                                ispin, natom, nkind, nspins, &
                                                stat, unit_nr
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind, kind_of, &
                                                natom_of_kind
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: Ehartree, Exc
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, rho_ao
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(harris_force_type), POINTER         :: harris_force
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type)                          :: rho_tot_gspace, &
                                                v_hartree_gspace, &
                                                v_hartree_rspace
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: my_rho_r, rho_r, &
                                                v_rspace_new, v_tau_rspace, &
                                                v_xc
    TYPE(pw_p_type), POINTER                 :: rho_core
    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force, force_tmp
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(section_vals_type), POINTER         :: input, xc_section
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    TYPE(xc_rho_set_type), POINTER           :: rho_set

!E_overlap_core
!   ------------------------------------------------------------------------

    CALL timeset(routineN,handle)

    failure = .FALSE.
    NULLIFY(force, force_tmp, v_rspace_new, v_tau_rspace, rho, matrix_ks, &
            pw_env, auxbas_pw_pool, rho_core, cell, dft_control, xc_section, &
            atomic_kind_set, particle_set, deriv_set, rho_set, my_rho_r, &
            v_xc, harris_force, rho_r, ks_env, poisson_env, input, rho_ao)
    para_env=>qs_env%para_env
    logger => cp_error_get_logger(error)

    CPPrecondition(ASSOCIATED(qs_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(qs_env%ref_count>0, cp_failure_level, routineP, error, failure)
    CPPrecondition(ASSOCIATED(harris_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(harris_env%ref_count>0, cp_failure_level, routineP, error, failure)
    CPPrecondition(ASSOCIATED(globenv), cp_failure_level, routineP, error, failure)
    CPPrecondition(globenv%ref_count>0, cp_failure_level, routineP, error, failure)

    IF (.NOT. failure) THEN
      CALL harris_env_get(harris_env=harris_env, harris_force=harris_force, error=error)
      CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set, force=force, &
                      dft_control=dft_control, particle_set=particle_set, &
                      cell=cell, input=input, error=error)

      xc_section => section_vals_get_subs_vals(input, "DFT%XC", error=error)

      natom = SIZE(particle_set)
      nkind = SIZE(force)
      nspins = dft_control%nspins

      ALLOCATE (atom_of_kind(natom),STAT=stat)
      CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

      ALLOCATE (kind_of(natom),STAT=stat)
      CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

      ALLOCATE (natom_of_kind(nkind),STAT=stat)
      CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                               atom_of_kind=atom_of_kind, kind_of=kind_of, &
                               natom_of_kind=natom_of_kind)

      CALL duplicate_qs_force(qs_force_input=force, qs_force_output=force_tmp, &
                              natom_of_kind=natom_of_kind)
      CALL zero_qs_force(force)

      ! *** d/dR[Sum of eigenvalues] = d/dR[trace(H_in * rho_out)] *** !
      CALL build_core_hamiltonian_matrix(qs_env=qs_env,calculate_forces=.TRUE.,error=error)
      CALL calculate_ecore_self(qs_env,error=error)
      CALL calculate_ecore_overlap(qs_env, para_env, &
                                           calculate_forces=.TRUE.,&
                                           error=error)

      DO ikind = 1,SIZE(force)
        CALL mp_sum(force(ikind)%kinetic, para_env%group)
        CALL mp_sum(force(ikind)%overlap, para_env%group)
        CALL mp_sum(force(ikind)%gth_ppl, para_env%group)
        CALL mp_sum(force(ikind)%gth_nlcc, para_env%group)
        CALL mp_sum(force(ikind)%gth_ppnl, para_env%group)
        CALL mp_sum(force(ikind)%core_overlap, para_env%group)
      END DO

      DO iatom = 1, natom
        ikind = kind_of(iatom)
        i = atom_of_kind(iatom)
        harris_force%f_kinetic(iatom, 1:3) = force(ikind)%kinetic(1:3,i)
        harris_force%f_overlap(iatom, 1:3) = force(ikind)%overlap(1:3,i)
        harris_force%f_gth_pp(iatom, 1:3)  = force(ikind)%gth_ppl(1:3,i) &
                                           + force(ikind)%gth_ppnl(1:3,i) + force(ikind)%gth_nlcc(1:3,i)
        harris_force%f_ovrl(iatom, 1:3) = force(ikind)%core_overlap(1:3,i)
        harris_force%f_self(iatom, 1:3) = 0.0_dp
      END DO

      CALL zero_qs_force(force)

      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env, rho=rho, rho_core=rho_core, &
                      ks_env=ks_env, matrix_ks=matrix_ks, error=error)

      CALL pw_env_get(pw_env=pw_env, auxbas_pw_pool=auxbas_pw_pool, &
                      poisson_env=poisson_env, error=error)
      CALL pw_pool_create_pw(auxbas_pw_pool, v_hartree_gspace%pw, &
                              use_data = COMPLEXDATA1D, &
                              in_space = RECIPROCALSPACE, error=error)
      CALL pw_pool_create_pw(auxbas_pw_pool, rho_tot_gspace%pw, &
                              use_data=COMPLEXDATA1D, &
                              in_space=RECIPROCALSPACE, error=error)

      CALL pw_pool_create_pw(auxbas_pw_pool, v_hartree_rspace%pw, &
                              use_data=REALDATA3D, in_space=REALSPACE,error=error)

      ! *** Calculation of the forces due to the Hartree energy *** !
      CALL pw_copy(rho_core%pw,rho_tot_gspace%pw,error=error)
      DO ispin = 1,nspins
        CALL pw_axpy(harris_env%rho%rho_g(ispin)%pw,rho_tot_gspace%pw,error=error)
        ! old code seems to be wrong for spin polarized case
        ! CALL pw_add(rho_core%pw, harris_env%rho%rho_g(ispin)%pw, rho_tot_gspace%pw)
      END DO

      CALL pw_poisson_solve(poisson_env, rho_tot_gspace%pw, ehartree=Ehartree, &
                            vhartree=v_hartree_gspace%pw,error=error)

      CALL pw_transfer(v_hartree_gspace%pw, v_hartree_rspace%pw, error=error)

      CALL qs_vxc_create(qs_env=qs_env, rho_struct=qs_env%rho, xc_section=xc_section, &
                         vxc_rho=v_rspace_new, vxc_tau=v_tau_rspace, exc=Exc, &
                         just_energy=.FALSE., error=error)

      rho_ao => qs_env%rho%rho_ao

      DO ispin = 1,nspins
        CALL pw_copy(v_rspace_new(ispin)%pw,v_rspace_new(ispin)%pw,error=error)
        CALL pw_axpy(v_hartree_rspace%pw,v_rspace_new(ispin)%pw,error=error)
        CALL pw_scale(v_rspace_new(ispin)%pw, v_rspace_new(ispin)%pw%pw_grid%dvol, error=error)

        CALL integrate_v_rspace(v_rspace=v_rspace_new(ispin), p=rho_ao(ispin), &
                                h=matrix_ks(ispin), qs_env=qs_env, calculate_forces=.TRUE., &
                                gapw=.FALSE., error=error)
      END DO

      DO ikind = 1,SIZE(force)
        CALL mp_sum(force(ikind)%rho_elec, para_env%group)
      END DO

      DO iatom = 1, natom
        ikind = kind_of(iatom)
        i = atom_of_kind(iatom)
        harris_force%f_hartree(iatom, 1:3) = force(ikind)%rho_elec(1:3,i)
      END DO

      CALL zero_qs_force(force)

      ! *** Calculation of the forces of the xc-integral *** !
      CALL qs_rho_get(rho_struct=rho, rho_r=rho_r, error=error)
      ALLOCATE(my_rho_r(SIZE(rho_r)), stat=stat)
      CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
      DO ispin=1, SIZE(rho_r)
        my_rho_r(ispin)%pw => harris_env%rho%rho_r(ispin)%pw
      END DO
      CALL xc_prep_2nd_deriv(deriv_set=deriv_set, rho_set=rho_set, &
                             rho_r=my_rho_r, pw_pool=auxbas_pw_pool, &
                             xc_section=xc_section, cell=cell, error=error)
      DEALLOCATE(my_rho_r, stat=stat)
      CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)

      ! ** v_xc ** !
      ALLOCATE(v_xc(nspins), stat=stat)
      CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
      DO ispin=1, nspins
        NULLIFY(v_xc(ispin)%pw)
        CALL pw_pool_create_pw(auxbas_pw_pool, v_xc(ispin)%pw, use_data=REALDATA3D, &
                               in_space=REALSPACE, error=error)
        CALL pw_zero(v_xc(ispin)%pw, error=error)
      END DO

      CALL xc_calc_2nd_deriv(v_xc, deriv_set=deriv_set, rho_set=rho_set, &
                             rho1_set=rho_set, pw_pool=auxbas_pw_pool, &
                             xc_section=xc_section, gapw=.FALSE., error=error)

      CALL xc_dset_release(deriv_set, error=error)

      DO ispin = 1,nspins
        v_xc(ispin)%pw%cr3d = v_xc(ispin)%pw%cr3d * v_xc(ispin)%pw%pw_grid%dvol
        v_rspace_new(ispin)%pw%cr3d = (qs_env%rho%rho_r(ispin)%pw%cr3d &
                                    - harris_env%rho%rho_r(ispin)%pw%cr3d) &
                                    * v_xc(ispin)%pw%cr3d
        CALL pw_release(v_xc(ispin)%pw, error=error)
      END DO

      DEALLOCATE(v_xc,stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      CALL xc_rho_set_release(rho_set,error=error)

      ! + v_hartree[rho_out-rho_in]
      CALL pw_zero(rho_tot_gspace%pw, error=error)
      DO ispin = 1,nspins
        CALL pw_axpy(qs_env%rho%rho_g(ispin)%pw,rho_tot_gspace%pw,error=error)
        CALL pw_axpy(harris_env%rho%rho_g(ispin)%pw,rho_tot_gspace%pw,alpha=-1._dp,error=error)
        ! this original code seems to be wrong for spin polarized cases (overwrites rho_tot_gspace)
        ! CALL pw_subtract(qs_env%rho%rho_g(ispin)%pw, harris_env%rho%rho_g(ispin)%pw, rho_tot_gspace%pw)
      END DO

      CALL pw_poisson_solve(poisson_env, rho_tot_gspace%pw, ehartree=Ehartree, &
                             vhartree=v_hartree_gspace%pw,error=error)

      CALL pw_transfer(v_hartree_gspace%pw, v_hartree_rspace%pw, error=error)


      rho_ao => harris_env%rho%rho_ao

      DO ispin = 1,nspins
        CALL pw_copy(v_rspace_new(ispin)%pw,v_rspace_new(ispin)%pw,error=error)
        CALL pw_axpy(v_hartree_rspace%pw,v_rspace_new(ispin)%pw,error=error)
        CALL pw_scale(v_rspace_new(ispin)%pw, v_rspace_new(ispin)%pw%pw_grid%dvol, error=error)

        CALL integrate_v_rspace(v_rspace=v_rspace_new(ispin), p=rho_ao(ispin), &
                                h=matrix_ks(ispin), qs_env=qs_env, calculate_forces=.TRUE., &
                                gapw=.FALSE., error=error)
      END DO

      IF (ASSOCIATED(v_rspace_new)) THEN
        DO ispin = 1,nspins
          CALL pw_pool_give_back_pw(auxbas_pw_pool, v_rspace_new(ispin)%pw, error=error)
        END DO
        DEALLOCATE(v_rspace_new, stat=stat)
        CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)
      END IF
      IF (ASSOCIATED(v_tau_rspace)) THEN
        DO ispin = 1,nspins
          CALL pw_pool_give_back_pw(auxbas_pw_pool, v_tau_rspace(ispin)%pw, error=error)
        END DO
        DEALLOCATE(v_tau_rspace, stat=stat)
        CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)
      END IF

      DO ikind = 1,SIZE(force)
        CALL mp_sum(force(ikind)%rho_elec, para_env%group)
      END DO

      DO iatom = 1, natom
        ikind = kind_of(iatom)
        i = atom_of_kind(iatom)
        harris_force%f_integral_vxc(iatom, 1:3) = force(ikind)%rho_elec(1:3,i)
      END DO

      DO iatom = 1, natom
        ikind = kind_of(iatom)
        i = atom_of_kind(iatom)
        harris_force%f_trace(iatom, 1:3) = harris_force%f_kinetic(iatom, 1:3) &
                                         + harris_force%f_gth_pp(iatom, 1:3) &
                                         + harris_force%f_overlap(iatom, 1:3) &
                                         + harris_force%f_hartree(iatom, 1:3) &
                                         + harris_force%f_integral_vxc(iatom, 1:3)
      END DO

      DO iatom = 1,natom
        harris_force%f_hartree(iatom, 1:3) = 0.0_dp
      END DO

      CALL zero_qs_force(force)

      ! *** The forces due to n_core *** !
      CALL pw_copy(rho_core%pw,rho_tot_gspace%pw,error=error)
      DO ispin = 1,nspins
        CALL pw_axpy(qs_env%rho%rho_g(ispin)%pw,rho_tot_gspace%pw,error=error)
        ! old code seems to be wrong for spin polarized case
        ! CALL pw_add(rho_core%pw, qs_env%rho%rho_g(ispin)%pw, rho_tot_gspace%pw)
      END DO

      CALL pw_poisson_solve(poisson_env, rho_tot_gspace%pw, ehartree=Ehartree, &
                            vhartree=v_hartree_gspace%pw,error=error)

      CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_tot_gspace%pw, error=error)

      CALL pw_transfer(v_hartree_gspace%pw, v_hartree_rspace%pw, error=error)
      CALL pw_scale(v_hartree_rspace%pw, v_hartree_rspace%pw%pw_grid%dvol, error=error)

      CALL integrate_v_core_rspace(v_hartree_rspace, qs_env,error=error)

      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_gspace%pw, error=error)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_rspace%pw, error=error)

      DO ikind = 1,SIZE(force)
        CALL mp_sum(force(ikind)%rho_core, para_env%group)
      END DO

      DO iatom = 1, natom
        ikind = kind_of(iatom)
        i = atom_of_kind(iatom)
        harris_force%f_rho_core(iatom, 1:3) = force(ikind)%rho_core(1:3,i)
      END DO

      CALL deallocate_qs_force(force)
      CALL set_qs_env(qs_env=qs_env, force=force_tmp, error=error)

      CALL get_qs_env(qs_env=qs_env, force=force, error=error)
      harris_force%f_total(:) = 0.0_dp

      DO iatom = 1, natom
        ikind = kind_of(iatom)
        i = atom_of_kind(iatom)

        harris_force%f_EII(iatom, 1:3) = harris_force%f_ovrl(iatom, 1:3) &
                                       + harris_force%f_self(iatom, 1:3) &
                                       + harris_force%f_rho_core(iatom, 1:3)
        harris_force%f_harris(iatom, 1:3) = harris_force%f_trace(iatom, 1:3) &
                                          + harris_force%f_EII(iatom, 1:3)
        harris_force%f_total(1:3) = harris_force%f_total(1:3) &
                                  + harris_force%f_harris(iatom, 1:3)
      END DO

      ! Output
      unit_nr=cp_logger_get_default_io_unit(logger)
      IF (unit_nr>0) THEN
        WRITE (unit_nr,*) ""; WRITE (unit_nr, *) ""
        WRITE (unit_nr,*) "The Harris functional force correction is performed!"
        WRITE (unit_nr,*) ""

        WRITE (unit_nr, *) "F_Harris"
        DO iatom = 1, natom
          WRITE (unit_nr,*) harris_force%f_harris(iatom, 1:3)
        END DO
        WRITE (unit_nr, *) "F_Total"
        WRITE (unit_nr, *) harris_force%f_total(1:3)
      END IF

      DEALLOCATE (atom_of_kind,STAT=stat)
      CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)

      DEALLOCATE (kind_of,STAT=stat)
      CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)

      DEALLOCATE (natom_of_kind,STAT=stat)
      CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)

    END IF

    CALL timestop(handle)

  END SUBROUTINE harris_force_correction

! *****************************************************************************
!> \brief Performs the harris functional force correction
!> \param qs_env The QS environment of matter
!> \param harris_force The harris force type in which all necessary informations
!>                      are stored
!> \param globenv The global environment from which the para_env is extracted
!> \param error variable to control error logging, stopping,...
!>               see module cp_error_handling
!> \note
!>      Performs the force calculation according to the Harris energy functional
!> \par History
!>      11.2005 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
  SUBROUTINE harris_force_EVal(qs_env, harris_env, globenv, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(harris_env_type), POINTER           :: harris_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'harris_force_EVal', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iatom, ikind, &
                                                ispin, natom, nkind, nspins, &
                                                stat, unit_nr
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind, kind_of, &
                                                natom_of_kind
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: E_overlap_core, Ehartree, Exc
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, rho_ao
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(harris_force_type), POINTER         :: harris_force
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type)                          :: rho_tot_gspace, &
                                                v_hartree_gspace, &
                                                v_hartree_rspace
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: my_rho_r, rho_r, &
                                                v_rspace_new, v_tau_rspace, &
                                                v_xc
    TYPE(pw_p_type), POINTER                 :: rho_core
    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force, force_tmp
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(section_vals_type), POINTER         :: input, xc_section
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    TYPE(xc_rho_set_type), POINTER           :: rho_set

!   ------------------------------------------------------------------------

    CALL timeset(routineN,handle)
    para_env=>qs_env%para_env

    failure = .FALSE.
    NULLIFY(force, force_tmp, v_rspace_new, v_tau_rspace, rho, matrix_ks, &
            pw_env, auxbas_pw_pool, rho_core, cell, dft_control, xc_section, &
            atomic_kind_set, particle_set, deriv_set, rho_set, my_rho_r, &
            v_xc, harris_force, rho_r, ks_env, poisson_env, input)
    logger => cp_error_get_logger(error)

    CPPrecondition(ASSOCIATED(qs_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(qs_env%ref_count>0, cp_failure_level, routineP, error, failure)
    CPPrecondition(ASSOCIATED(harris_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(harris_env%ref_count>0, cp_failure_level, routineP, error, failure)
    CPPrecondition(ASSOCIATED(globenv), cp_failure_level, routineP, error, failure)
    CPPrecondition(globenv%ref_count>0, cp_failure_level, routineP, error, failure)

    IF (.NOT. failure) THEN
      CALL harris_env_get(harris_env=harris_env, harris_force=harris_force, error=error)
      CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set, force=force, &
                      dft_control=dft_control, cell=cell, particle_set=particle_set, &
                      input=input, error=error)

      xc_section => section_vals_get_subs_vals(input, "DFT%XC", error=error)

      natom = SIZE(particle_set)
      nkind = SIZE(force)
      nspins = dft_control%nspins

      ALLOCATE (atom_of_kind(natom),STAT=stat)
      CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

      ALLOCATE (kind_of(natom),STAT=stat)
      CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

      ALLOCATE (natom_of_kind(nkind),STAT=stat)
      CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                               atom_of_kind=atom_of_kind, kind_of=kind_of, &
                               natom_of_kind=natom_of_kind)

      CALL duplicate_qs_force(qs_force_input=force, qs_force_output=force_tmp, &
                              natom_of_kind=natom_of_kind)
      CALL zero_qs_force(force)

      ! *** d/dR[Sum of eigenvalues] = d/dR[trace(H_in * rho_out)] *** !
      CALL build_core_hamiltonian_matrix(qs_env=qs_env,calculate_forces=.TRUE.,error=error)

      DO ikind = 1,SIZE(force)
        CALL mp_sum(force(ikind)%kinetic, para_env%group)
        CALL mp_sum(force(ikind)%overlap, para_env%group)
        CALL mp_sum(force(ikind)%gth_ppl, para_env%group)
        CALL mp_sum(force(ikind)%gth_nlcc, para_env%group)
        CALL mp_sum(force(ikind)%gth_ppnl, para_env%group)
      END DO

      DO iatom = 1, natom
        ikind = kind_of(iatom)
        i = atom_of_kind(iatom)
        harris_force%f_kinetic(iatom, 1:3) = force(ikind)%kinetic(1:3,i)
        harris_force%f_overlap(iatom, 1:3) = force(ikind)%overlap(1:3,i)
        harris_force%f_gth_pp(iatom, 1:3)  = force(ikind)%gth_ppl(1:3,i) &
                                           + force(ikind)%gth_ppnl(1:3,i) + force(ikind)%gth_nlcc(1:3,i)
      END DO

      CALL zero_qs_force(force)

      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env, rho=rho, rho_core=rho_core, &
                      ks_env=ks_env, matrix_ks=matrix_ks, &
                      error=error)

      CALL pw_env_get(pw_env=pw_env, auxbas_pw_pool=auxbas_pw_pool, &
                      poisson_env=poisson_env, error=error)
      CALL pw_pool_create_pw(auxbas_pw_pool, v_hartree_gspace%pw, &
                              use_data = COMPLEXDATA1D, &
                              in_space = RECIPROCALSPACE, error=error)
      CALL pw_pool_create_pw(auxbas_pw_pool, rho_tot_gspace%pw, &
                              use_data=COMPLEXDATA1D, &
                              in_space=RECIPROCALSPACE, error=error)

      CALL pw_pool_create_pw(auxbas_pw_pool, v_hartree_rspace%pw, &
                              use_data=REALDATA3D, in_space=REALSPACE,error=error)

      ! *** Calculation of the forces due to the Hartree energy *** !
      CALL pw_copy(rho_core%pw,rho_tot_gspace%pw,error=error)
      DO ispin = 1,nspins
        CALL pw_axpy(harris_env%rho%rho_g(ispin)%pw,rho_tot_gspace%pw,error=error)
        ! old code seems to be wrong for spin polarized case
        ! CALL pw_add(rho_core%pw, harris_env%rho%rho_g(ispin)%pw, rho_tot_gspace%pw)
      END DO

      CALL pw_poisson_solve(poisson_env, rho_tot_gspace%pw, ehartree=Ehartree, &
                            vhartree=v_hartree_gspace%pw,error=error)

      CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_tot_gspace%pw, error=error)

      CALL pw_transfer(v_hartree_gspace%pw, v_hartree_rspace%pw, error=error)

      CALL qs_vxc_create(qs_env=qs_env, rho_struct=qs_env%rho, xc_section=xc_section, &
                         vxc_rho=v_rspace_new, vxc_tau=v_tau_rspace, exc=Exc, &
                         just_energy=.FALSE., error=error)

      rho_ao => qs_env%rho%rho_ao

      DO ispin = 1,nspins
        CALL pw_copy(v_rspace_new(ispin)%pw,v_rspace_new(ispin)%pw,error=error)
        CALL pw_axpy(v_hartree_rspace%pw,v_rspace_new(ispin)%pw,error=error)
        CALL pw_scale(v_rspace_new(ispin)%pw, v_rspace_new(ispin)%pw%pw_grid%dvol, error=error)

        CALL integrate_v_rspace(v_rspace=v_rspace_new(ispin), p=rho_ao(ispin), &
                                h=matrix_ks(ispin), qs_env=qs_env, calculate_forces=.TRUE., &
                                gapw=.FALSE., error=error)
      END DO

      DO ikind = 1,SIZE(force)
        CALL mp_sum(force(ikind)%rho_elec, para_env%group)
      END DO

      DO iatom = 1, natom
        ikind = kind_of(iatom)
        i = atom_of_kind(iatom)
        harris_force%f_hartree(iatom, 1:3) = force(ikind)%rho_elec(1:3,i)
      END DO

      CALL zero_qs_force(force)

      ! ** Force computation of f_V  ** !
      DO ispin = 1,nspins
        CALL pw_poisson_solve(poisson_env, qs_env%rho%rho_g(ispin)%pw, ehartree=Ehartree, &
                              vhartree=v_hartree_gspace%pw,error=error)
      END DO

      CALL pw_transfer(v_hartree_gspace%pw, v_hartree_rspace%pw, error=error)
      CALL pw_scale(v_hartree_rspace%pw, v_hartree_rspace%pw%pw_grid%dvol, error=error)

      CALL integrate_v_core_rspace(v_hartree_rspace, qs_env,error=error)

      DO ikind = 1,SIZE(force)
        CALL mp_sum(force(ikind)%rho_core, para_env%group)
      END DO

      DO iatom = 1, natom
        ikind = kind_of(iatom)
        i = atom_of_kind(iatom)
        harris_force%f_V(iatom, 1:3) = force(ikind)%rho_core(1:3,i)
      END DO

      CALL zero_qs_force(force)

      ! *** Calculation of the forces of the xc-integral *** !
      CALL qs_rho_get(rho_struct=rho, rho_r=rho_r, error=error)
      ALLOCATE(my_rho_r(SIZE(rho_r)), stat=stat)
      CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
      DO ispin=1, SIZE(rho_r)
        my_rho_r(ispin)%pw => harris_env%rho%rho_r(ispin)%pw
      END DO
      CALL xc_prep_2nd_deriv(deriv_set=deriv_set, rho_set=rho_set, &
                             rho_r=my_rho_r, pw_pool=auxbas_pw_pool, &
                             xc_section=xc_section, cell=cell, error=error)
      DEALLOCATE(my_rho_r, stat=stat)
      CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)

      ! ** v_xc ** !
      ALLOCATE(v_xc(nspins), stat=stat)
      CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
      DO ispin=1, nspins
        NULLIFY(v_xc(ispin)%pw)
        CALL pw_pool_create_pw(auxbas_pw_pool, v_xc(ispin)%pw, use_data=REALDATA3D, &
                               in_space=REALSPACE, error=error)
        CALL pw_zero(v_xc(ispin)%pw, error=error)
      END DO

      CALL xc_calc_2nd_deriv(v_xc, deriv_set=deriv_set, rho_set=rho_set, &
                             rho1_set=rho_set, pw_pool=auxbas_pw_pool, &
                             xc_section=xc_section, gapw=.FALSE., error=error)

      CALL xc_dset_release(deriv_set, error=error)

      DO ispin = 1,nspins
        v_xc(ispin)%pw%cr3d = v_xc(ispin)%pw%cr3d * v_xc(ispin)%pw%pw_grid%dvol
        v_rspace_new(ispin)%pw%cr3d = rho_r(ispin)%pw%cr3d * v_xc(ispin)%pw%cr3d

        CALL pw_release(v_xc(ispin)%pw, error=error)
      END DO

      DEALLOCATE(v_xc,stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      CALL xc_rho_set_release(rho_set,error=error)

      ! + v_hartree[rho_out]
      DO ispin = 1,nspins
        CALL pw_poisson_solve(poisson_env, qs_env%rho%rho_g(ispin)%pw, ehartree=Ehartree, &
                               vhartree=v_hartree_gspace%pw,error=error)
      END DO

      CALL pw_transfer(v_hartree_gspace%pw, v_hartree_rspace%pw, error=error)


    rho_ao => harris_env%rho%rho_ao

      DO ispin = 1,nspins
        CALL pw_copy(v_rspace_new(ispin)%pw,v_rspace_new(ispin)%pw,error=error)
        CALL pw_axpy(v_hartree_rspace%pw,v_rspace_new(ispin)%pw,error=error)
        CALL pw_scale(v_rspace_new(ispin)%pw, v_rspace_new(ispin)%pw%pw_grid%dvol, error=error)

        CALL integrate_v_rspace(v_rspace=v_rspace_new(ispin), p=rho_ao(ispin), &
                                h=matrix_ks(ispin), qs_env=qs_env, calculate_forces=.TRUE., &
                                gapw=.FALSE., error=error)
      END DO

      IF (ASSOCIATED(v_rspace_new)) THEN
        DO ispin = 1,nspins
          CALL pw_pool_give_back_pw(auxbas_pw_pool,v_rspace_new(ispin)%pw,error=error)
        END DO
        DEALLOCATE(v_rspace_new,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      END IF
      IF (ASSOCIATED(v_tau_rspace)) THEN
        DO ispin = 1,nspins
          CALL pw_pool_give_back_pw(auxbas_pw_pool,v_tau_rspace(ispin)%pw,error=error)
        END DO
        DEALLOCATE(v_tau_rspace,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      END IF

      DO ikind = 1,SIZE(force)
        CALL mp_sum(force(ikind)%rho_elec, para_env%group)
      END DO

      DO iatom = 1, natom
        ikind = kind_of(iatom)
        i = atom_of_kind(iatom)
        harris_force%f_integral_vxc(iatom, 1:3) = force(ikind)%rho_elec(1:3,i)
      END DO

      DO iatom = 1, natom
        ikind = kind_of(iatom)
        i = atom_of_kind(iatom)
        harris_force%f_trace(iatom, 1:3) = harris_force%f_kinetic(iatom, 1:3) &
                                         + harris_force%f_gth_pp(iatom, 1:3) &
                                         + harris_force%f_overlap(iatom, 1:3) &
                                         + harris_force%f_V(iatom, 1:3) &
                                         + harris_force%f_hartree(iatom, 1:3) &
                                         + harris_force%f_integral_vxc(iatom, 1:3)
      END DO

      DO iatom = 1,natom
        harris_force%f_integral_vxc(iatom, 1:3) = 0.0_dp
      END DO

      CALL zero_qs_force(force)
      ! *** End of force calculation due to the sum_of_eigenvalues *** !

      CALL harris_env_get(harris_env=harris_env, rho=rho, error=error)

      ! *** The forces due to the core-core repulsion *** !
      CALL pw_poisson_solve(poisson_env, rho_core%pw, ehartree=Ehartree, &
                            vhartree=v_hartree_gspace%pw,error=error)

      CALL pw_transfer(v_hartree_gspace%pw, v_hartree_rspace%pw, error=error)
      CALL pw_scale(v_hartree_rspace%pw, v_hartree_rspace%pw%pw_grid%dvol, error=error)

      CALL integrate_v_core_rspace(v_hartree_rspace, qs_env,error=error)

      DO ikind = 1,SIZE(force)
        CALL mp_sum(force(ikind)%rho_core, para_env%group)
      END DO

      DO iatom = 1, natom
        ikind = kind_of(iatom)
        i = atom_of_kind(iatom)
        harris_force%f_rho_core(iatom, 1:3) = force(ikind)%rho_core(1:3,i)
      END DO

      CALL zero_qs_force(force)

      ! *** The forces due to exchange and correlation *** !
      CALL qs_vxc_create(qs_env=qs_env, rho_struct=qs_env%rho, xc_section=xc_section, &
                         vxc_rho=v_rspace_new, vxc_tau=v_tau_rspace, exc=Exc, &
                         just_energy=.FALSE., error=error)

      ! *** Calculation of the forces of the xc-integral *** !
      CALL qs_rho_get(rho, rho_r=rho_r, error=error)
      ALLOCATE(my_rho_r(SIZE(rho_r)), stat=stat)
      CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
      DO ispin=1, SIZE(rho_r)
        my_rho_r(ispin)%pw => harris_env%rho%rho_r(ispin)%pw
      END DO
      CALL xc_prep_2nd_deriv(deriv_set=deriv_set, rho_set=rho_set, &
                             rho_r=my_rho_r, pw_pool=auxbas_pw_pool, &
                             xc_section=xc_section, cell=cell, error=error)
      DEALLOCATE(my_rho_r, stat=stat)
      CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)

      ! ** v_xc ** !
      nspins = dft_control%nspins
      ALLOCATE(v_xc(nspins), stat=stat)
      CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
      DO ispin=1, nspins
        NULLIFY(v_xc(ispin)%pw)
        CALL pw_pool_create_pw(auxbas_pw_pool, v_xc(ispin)%pw, use_data=REALDATA3D, &
                               in_space=REALSPACE, error=error)
        CALL pw_zero(v_xc(ispin)%pw, error=error)
      END DO

      CALL xc_calc_2nd_deriv(v_xc, deriv_set=deriv_set, rho_set=rho_set, &
                             rho1_set=rho_set, pw_pool=auxbas_pw_pool, &
                             xc_section=xc_section, gapw=.FALSE., error=error)

      CALL xc_dset_release(deriv_set, error=error)

      DO ispin = 1,nspins
        v_xc(ispin)%pw%cr3d = v_xc(ispin)%pw%cr3d * v_xc(ispin)%pw%pw_grid%dvol
        v_rspace_new(ispin)%pw%cr3d = rho_r(ispin)%pw%cr3d * v_xc(ispin)%pw%cr3d
        CALL pw_scale(v_rspace_new(ispin)%pw, v_rspace_new(ispin)%pw%pw_grid%dvol, error=error)

        CALL pw_release(v_xc(ispin)%pw, error=error)
      END DO

      DEALLOCATE(v_xc,stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      CALL xc_rho_set_release(rho_set,error=error)

      DO ispin = 1,nspins
        CALL pw_poisson_solve(poisson_env, rho%rho_g(ispin)%pw, ehartree=Ehartree, &
                              vhartree=v_hartree_gspace%pw,error=error)
      END DO

      CALL pw_transfer(v_hartree_gspace%pw, v_hartree_rspace%pw, error=error)
      CALL pw_scale(v_hartree_rspace%pw, v_hartree_rspace%pw%pw_grid%dvol, error=error)

    rho_ao => harris_env%rho%rho_ao
      DO ispin = 1,nspins
        CALL pw_copy(v_rspace_new(ispin)%pw,v_rspace_new(ispin)%pw,error=error)
        CALL pw_axpy(v_hartree_rspace%pw,v_rspace_new(ispin)%pw,error=error)

        CALL integrate_v_rspace(v_rspace=v_rspace_new(ispin), p=rho_ao(ispin), &
                                h=matrix_ks(ispin), qs_env=qs_env, calculate_forces=.TRUE., &
                                gapw=.FALSE., error=error)
      END DO

      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_gspace%pw, error=error)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_rspace%pw, error=error)

      IF (ASSOCIATED(v_rspace_new)) THEN
        DO ispin = 1,nspins
          CALL pw_pool_give_back_pw(auxbas_pw_pool,v_rspace_new(ispin)%pw,error=error)
        END DO
        DEALLOCATE(v_rspace_new,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      END IF
      IF (ASSOCIATED(v_tau_rspace)) THEN
        DO ispin = 1,nspins
          CALL pw_pool_give_back_pw(auxbas_pw_pool,v_tau_rspace(ispin)%pw,error=error)
        END DO
        DEALLOCATE(v_tau_rspace,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      END IF

      DO ikind = 1,SIZE(force)
        CALL mp_sum(force(ikind)%rho_elec, para_env%group)
      END DO

      DO iatom = 1, natom
        ikind = kind_of(iatom)
        i = atom_of_kind(iatom)
        harris_force%f_integral_vxc(iatom, 1:3) = force(ikind)%rho_elec(1:3,i)
      END DO

      CALL zero_qs_force(force)

      CALL calculate_ecore_overlap(qs_env=qs_env, para_env=para_env, &
                                   calculate_forces=.TRUE., E_overlap_core=E_overlap_core,&
                                   error=error)

      DO iatom = 1, natom
        ikind = kind_of(iatom)
        i = atom_of_kind(iatom)
        harris_force%f_ovrl(iatom, 1:3) = force(ikind)%core_overlap(1:3,i)
      END DO

      CALL deallocate_qs_force(force)
      CALL set_qs_env(qs_env=qs_env, force=force_tmp, error=error)

      CALL get_qs_env(qs_env=qs_env, force=force, error=error)
      harris_force%f_total(:) = 0.0_dp

      DO ikind = 1,SIZE(force)
        CALL mp_sum(force(ikind)%core_overlap, para_env%group)
      END DO

      DO iatom = 1, natom
        ikind = kind_of(iatom)
        i = atom_of_kind(iatom)

        !harris_force%f_ovrl(iatom, 1:3) = force(ikind)%core_overlap(1:3,i)
        harris_force%f_self(iatom, 1:3) = 0.0_dp
        harris_force%f_EII(iatom, 1:3) = harris_force%f_ovrl(iatom, 1:3) &
                                       + harris_force%f_self(iatom, 1:3) &
                                       + harris_force%f_rho_core(iatom, 1:3)
        harris_force%f_harris(iatom, 1:3) = harris_force%f_trace(iatom, 1:3) &
                                          - harris_force%f_integral_vxc(iatom, 1:3) &
                                          + harris_force%f_EII(iatom, 1:3)
        harris_force%f_total(1:3) = harris_force%f_total(1:3) &
                                  + harris_force%f_harris(iatom, 1:3)
      END DO

      ! Output
      unit_nr=cp_logger_get_default_io_unit(logger)
      IF (unit_nr>0) THEN
        WRITE (unit_nr,*) ""; WRITE (unit_nr, *) ""
        WRITE (unit_nr,*) "The Harris functional force correction is performed!"
        WRITE (unit_nr,*) ""

        WRITE (unit_nr, *) "F_Harris"
        DO iatom = 1, natom
          WRITE (unit_nr,*) harris_force%f_harris(iatom, 1:3)
        END DO
        WRITE (unit_nr, *) "F_Total"
        WRITE (unit_nr, *) harris_force%f_total(1:3)
      END IF

      DEALLOCATE (atom_of_kind,STAT=stat)
      CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)

      DEALLOCATE (kind_of,STAT=stat)
      CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)

      DEALLOCATE (natom_of_kind,STAT=stat)
      CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)

    END IF

    CALL timestop(handle)

  END SUBROUTINE harris_force_EVal

! *****************************************************************************
!> \brief Calculates the non-self-consistent force correction, in order to keep
!>      the forces consistent with the according energies.
!> \param error variable to control error logging, stopping,...
!>               see module cp_error_handling
!> \par History
!>      11.2005 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
  SUBROUTINE harris_calc_nsc_force(qs_env, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'harris_calc_nsc_force', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, iatom, ikind, &
                                                ispin, natom, nkind, nspins, &
                                                stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind, kind_of, &
                                                natom_of_kind
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: Ehartree, Exc
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, rho_ao
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(harris_env_type), POINTER           :: harris_env
    TYPE(harris_force_type), POINTER         :: harris_force
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type)                          :: v_hartree_gspace, &
                                                v_hartree_rspace
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: my_rho_r, rho_r, &
                                                v_rspace_new, v_tau_rspace, &
                                                v_xc
    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force, force_tmp
    TYPE(section_vals_type), POINTER         :: input, xc_section
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    TYPE(xc_rho_set_type), POINTER           :: rho_set

!   ------------------------------------------------------------------------

    CALL timeset(routineN,handle)

    failure = .FALSE.
    para_env=>qs_env%para_env
    NULLIFY(harris_env, harris_force, force, force_tmp, particle_set, &
            v_rspace_new, v_tau_rspace, rho_r, my_rho_r, v_xc, dft_control, &
            cell, deriv_set, rho_set, auxbas_pw_pool, matrix_ks, pw_env, &
            poisson_env, xc_section, input)

    CPPrecondition(ASSOCIATED(qs_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(qs_env%ref_count>0, cp_failure_level, routineP, error, failure)

    IF (.NOT. failure) THEN
      CALL get_qs_env(qs_env=qs_env, harris_env=harris_env, dft_control=dft_control, &
                      atomic_kind_set=atomic_kind_set, force=force, pw_env=pw_env, &
                      particle_set=particle_set, matrix_ks=matrix_ks, cell=cell, &
                      error=error)

      xc_section => section_vals_get_subs_vals(input, "DFT%XC", error=error)

      CALL harris_env_get(harris_env=harris_env, harris_force=harris_force, error=error)

      natom = SIZE(particle_set)
      nkind = SIZE(force)
      nspins = dft_control%nspins

      ALLOCATE (atom_of_kind(natom),STAT=stat)
      CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

      ALLOCATE (kind_of(natom),STAT=stat)
      CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

      ALLOCATE (natom_of_kind(nkind),STAT=stat)
      CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                               atom_of_kind=atom_of_kind, kind_of=kind_of, &
                               natom_of_kind=natom_of_kind)

      CALL diff_rho_type(rho_input1=qs_env%rho, rho_input2=harris_env%rho, &
                         rho_output=harris_env%rho_diff, qs_env=qs_env, error=error)

      CALL duplicate_qs_force(qs_force_input=force, qs_force_output=force_tmp, &
                              natom_of_kind=natom_of_kind)
      CALL zero_qs_force(force)

      CALL pw_env_get(pw_env=pw_env, auxbas_pw_pool=auxbas_pw_pool, &
                      poisson_env=poisson_env, error=error)
      CALL pw_pool_create_pw(auxbas_pw_pool, v_hartree_gspace%pw, &
                              use_data = COMPLEXDATA1D, &
                              in_space = RECIPROCALSPACE, error=error)
      CALL pw_pool_create_pw(auxbas_pw_pool, v_hartree_rspace%pw, &
                              use_data=REALDATA3D, in_space=REALSPACE,error=error)

      ! *** Calculation of the forces due to the Hartree energy *** !
      DO ispin = 1,nspins
        CALL pw_poisson_solve(poisson_env, harris_env%rho_diff%rho_g(ispin)%pw, &
                              ehartree=Ehartree, vhartree=v_hartree_gspace%pw,error=error)
      END DO

      CALL pw_transfer(v_hartree_gspace%pw, v_hartree_rspace%pw, error=error)
      CALL pw_scale(v_hartree_rspace%pw, v_hartree_rspace%pw%pw_grid%dvol, error=error)

    rho_ao => harris_env%rho_diff%rho_ao
      DO ispin = 1,nspins
        CALL integrate_v_rspace(v_rspace=v_hartree_rspace, &
                                p=rho_ao(ispin), &
                                h=matrix_ks(ispin), qs_env=qs_env, &
                                calculate_forces=.TRUE., gapw=.FALSE., &
                                error=error)
      END DO

      DO ikind = 1,SIZE(force)
        CALL mp_sum(force(ikind)%rho_elec, para_env%group)
      END DO

      DO iatom = 1, natom
        ikind = kind_of(iatom)
        i = atom_of_kind(iatom)
        harris_force%f_nsc(iatom, 1:3) = force(ikind)%rho_elec(1:3,i)
      END DO

      CALL zero_qs_force(force)

      DO ispin = 1,nspins
        CALL pw_poisson_solve(poisson_env, harris_env%rho_diff%rho_g(ispin)%pw, &
                              ehartree=Ehartree, vhartree=v_hartree_gspace%pw,error=error)
      END DO

      CALL pw_transfer(v_hartree_gspace%pw, v_hartree_rspace%pw, error=error)
      CALL pw_scale(v_hartree_rspace%pw, v_hartree_rspace%pw%pw_grid%dvol, error=error)


      rho_ao => harris_env%rho_diff%rho_ao
      DO ispin = 1,nspins
        CALL integrate_v_rspace(v_rspace=v_hartree_rspace, &
                                p=rho_ao(ispin), &
                                h=matrix_ks(ispin), qs_env=qs_env, &
                                calculate_forces=.TRUE., gapw=.FALSE., &
                                error=error)
      END DO

      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_gspace%pw, error=error)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_rspace%pw, error=error)

      DO ikind = 1,SIZE(force)
        CALL mp_sum(force(ikind)%rho_elec, para_env%group)
      END DO

      DO iatom = 1, natom
        ikind = kind_of(iatom)
        i = atom_of_kind(iatom)
        harris_force%f_nsc(iatom, 1:3) = harris_force%f_nsc(iatom, 1:3) &
                                       + force(ikind)%rho_elec(1:3,i)
      END DO

      CALL zero_qs_force(force)

      ! *** Calculation of the forces of the xc-integral *** !

      CALL qs_vxc_create(qs_env=qs_env, rho_struct=harris_env%rho_diff, xc_section=xc_section, &
                         vxc_rho=v_rspace_new, vxc_tau=v_tau_rspace, exc=Exc, &
                         just_energy=.FALSE., error=error)

      rho_ao => harris_env%rho_diff%rho_ao

      DO ispin = 1,nspins
        CALL pw_scale(v_rspace_new(ispin)%pw, v_rspace_new(ispin)%pw%pw_grid%dvol, error=error)
        CALL integrate_v_rspace(v_rspace=v_rspace_new(ispin), &
                                p=rho_ao(ispin), &
                                h=matrix_ks(ispin), qs_env=qs_env, &
                                calculate_forces=.TRUE., gapw=.FALSE., &
                                error=error)
      END DO

      DO ikind = 1,SIZE(force)
        CALL mp_sum(force(ikind)%rho_elec, para_env%group)
      END DO

      DO iatom = 1, natom
        ikind = kind_of(iatom)
        i = atom_of_kind(iatom)
        !harris_force%f_nsc(iatom, 1:3) = force(ikind)%rho_elec(1:3,i)
      END DO

      CALL zero_qs_force(force)

      CALL qs_rho_get(rho_struct=harris_env%rho_diff, rho_r=rho_r, error=error)
      !CALL qs_rho_get(rho_struct=harris_env%rho_diff, rho_r=rho_r, error=error)
      ALLOCATE(my_rho_r(SIZE(rho_r)), stat=stat)
      CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
      DO ispin=1, SIZE(rho_r)
        my_rho_r(ispin)%pw => rho_r(ispin)%pw
      END DO
      CALL xc_prep_2nd_deriv(deriv_set=deriv_set, rho_set=rho_set, &
                             rho_r=my_rho_r, pw_pool=auxbas_pw_pool, &
                             xc_section=xc_section, cell=cell, error=error)
      DEALLOCATE(my_rho_r, stat=stat)
      CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)

      ! ** v_xc ** !
      ALLOCATE(v_xc(nspins), stat=stat)
      CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
      DO ispin=1, nspins
        NULLIFY(v_xc(ispin)%pw)
        CALL pw_pool_create_pw(auxbas_pw_pool, v_xc(ispin)%pw, use_data=REALDATA3D, &
                               in_space=REALSPACE, error=error)
        CALL pw_zero(v_xc(ispin)%pw, error=error)
      END DO

      CALL xc_calc_2nd_deriv(v_xc, deriv_set=deriv_set, rho_set=rho_set, &
                             rho1_set=rho_set, pw_pool=auxbas_pw_pool, &
                             xc_section=xc_section, gapw=.FALSE., error=error)

      CALL xc_dset_release(deriv_set, error=error)

      DO ispin = 1,nspins
        v_xc(ispin)%pw%cr3d = v_xc(ispin)%pw%cr3d * v_xc(ispin)%pw%pw_grid%dvol
        v_rspace_new(ispin)%pw%cr3d = rho_r(ispin)%pw%cr3d * v_xc(ispin)%pw%cr3d
        CALL pw_scale(v_rspace_new(ispin)%pw, v_rspace_new(ispin)%pw%pw_grid%dvol, error=error)

        CALL pw_release(v_xc(ispin)%pw, error=error)
      END DO

      DEALLOCATE(v_xc,stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      CALL xc_rho_set_release(rho_set,error=error)

    rho_ao => harris_env%rho_diff%rho_ao

      DO ispin = 1,nspins
        CALL integrate_v_rspace(v_rspace=v_rspace_new(ispin), &
                                p=rho_ao(ispin), &
                                h=matrix_ks(ispin), qs_env=qs_env, &
                                calculate_forces=.TRUE., gapw=.FALSE., &
                                error=error)
      END DO

      IF (ASSOCIATED(v_rspace_new)) THEN
        DO ispin = 1,nspins
          CALL pw_pool_give_back_pw(auxbas_pw_pool,v_rspace_new(ispin)%pw,error=error)
        END DO
        DEALLOCATE(v_rspace_new,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      END IF
      IF (ASSOCIATED(v_tau_rspace)) THEN
        DO ispin = 1,nspins
          CALL pw_pool_give_back_pw(auxbas_pw_pool,v_tau_rspace(ispin)%pw,error=error)
        END DO
        DEALLOCATE(v_tau_rspace,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      END IF

      DO ikind = 1,SIZE(force)
        CALL mp_sum(force(ikind)%rho_elec, para_env%group)
      END DO

      DO iatom = 1, natom
        ikind = kind_of(iatom)
        i = atom_of_kind(iatom)
        !harris_force%f_nsc(iatom, 1:3) = harris_force%f_nsc(iatom,1:3) + force(ikind)%rho_elec(1:3,i)
      END DO

    END IF

    CALL timestop(handle)

  END SUBROUTINE harris_calc_nsc_force

! *****************************************************************************
!> \brief Performs a test if the forces in the harris functional due to
!>      rho_core sums up to force%rho_core.
!> \param qs_env The QS environment of matter
!> \param harris_force The harris force type in which all necessary informations
!>                      are stored
!> \param error variable to control error logging, stopping,...
!>               see module cp_error_handling
!> \note
!>      Only for internal debugging purposes
!> \par History
!>      06.2005 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
  SUBROUTINE harris_force_test_rho_core(qs_env, harris_force, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(harris_force_type), POINTER         :: harris_force
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'harris_force_test_rho_core', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iatom, ikind, natom, &
                                                nkind, stat, unit_nr
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind, kind_of, &
                                                natom_of_kind
    LOGICAL                                  :: failure
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force

!   ------------------------------------------------------------------------

    failure = .FALSE.
    logger => cp_error_get_logger(error)

    CPPrecondition(ASSOCIATED(qs_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(qs_env%ref_count>0, cp_failure_level, routineP, error, failure)
    CPPrecondition(ASSOCIATED(harris_force), cp_failure_level, routineP, error, failure)
    CPPrecondition(harris_force%ref_count>0, cp_failure_level, routineP, error, failure)

    unit_nr=cp_logger_get_default_io_unit(logger)
    IF (unit_nr>0) THEN
      WRITE (unit_nr,*) ""; WRITE (unit_nr, *) ""
      WRITE (unit_nr,*) "The Harris force correction rho_core test is performed!"
      WRITE (unit_nr,*) ""
    END IF

    IF (.NOT. failure) THEN
      CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set, force=force, &
                      particle_set=particle_set, error=error)

      natom = SIZE(particle_set)
      nkind = SIZE(force)

      ALLOCATE (atom_of_kind(natom),STAT=stat)
      CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

      ALLOCATE (kind_of(natom),STAT=stat)
      CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

      ALLOCATE (natom_of_kind(nkind),STAT=stat)
      CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                               atom_of_kind=atom_of_kind, kind_of=kind_of, &
                               natom_of_kind=natom_of_kind)

      IF (unit_nr>0) THEN
        DO  iatom = 1, natom
           ikind = kind_of(iatom)
           i = atom_of_kind(iatom)
           WRITE (unit_nr,*) force(ikind)%rho_core(1:3,i)
           WRITE (unit_nr,*) harris_force%f_rho_core(iatom,1:3) + &
           harris_force%f_cross_integrate_v_core(iatom,1:3)
        END DO
      END IF

      DEALLOCATE (atom_of_kind,STAT=stat)
      CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)

      DEALLOCATE (kind_of,STAT=stat)
      CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)

      DEALLOCATE (natom_of_kind,STAT=stat)
      CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)

    END IF

  END SUBROUTINE harris_force_test_rho_core

! *****************************************************************************
!> \brief Performs a test if the forces in the harris functional due to
!>      rho_elec sums up to force%rho_elec.
!> \param qs_env The QS environment of matter
!> \param harris_force The harris force type in which all necessary informations
!>                      are stored
!> \param error variable to control error logging, stopping,...
!>               see module cp_error_handling
!> \note
!>      Only for internal debugging purposes
!> \par History
!>      06.2005 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
  SUBROUTINE harris_force_test_rho_elec(qs_env, harris_force, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(harris_force_type), POINTER         :: harris_force
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'harris_force_test_rho_elec', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iatom, ikind, natom, &
                                                nkind, stat, unit_nr
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind, kind_of, &
                                                natom_of_kind
    LOGICAL                                  :: failure
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force

!   ------------------------------------------------------------------------

    failure = .FALSE.
    logger => cp_error_get_logger(error)

    CPPrecondition(ASSOCIATED(qs_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(qs_env%ref_count>0, cp_failure_level, routineP, error, failure)
    CPPrecondition(ASSOCIATED(harris_force), cp_failure_level, routineP, error, failure)
    CPPrecondition(harris_force%ref_count>0, cp_failure_level, routineP, error, failure)

    unit_nr=cp_logger_get_default_io_unit(logger)
    IF (unit_nr>0) THEN
          WRITE (unit_nr,*) ""; WRITE (unit_nr, *) ""
          WRITE (unit_nr,*) "The Harris force correction rho_elec test is performed!"
          WRITE (unit_nr,*) ""
    END IF

    IF (.NOT. failure) THEN
      CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set, force=force, &
                      particle_set=particle_set, error=error)

      natom = SIZE(particle_set)
      nkind = SIZE(force)

      ALLOCATE (atom_of_kind(natom),STAT=stat)
      CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

      ALLOCATE (kind_of(natom),STAT=stat)
      CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

      ALLOCATE (natom_of_kind(nkind),STAT=stat)
      CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                               atom_of_kind=atom_of_kind, kind_of=kind_of, &
                               natom_of_kind=natom_of_kind)

      IF (unit_nr>0) THEN
        DO  iatom = 1, natom
          ikind = kind_of(iatom)
          i = atom_of_kind(iatom)
          WRITE (unit_nr,*) force(ikind)%rho_elec(1:3,i)
          WRITE (unit_nr,*) harris_force%f_hartree(iatom,1:3) + &
             harris_force%f_xc(iatom,1:3) + &
             harris_force%f_cross_integrate_v(iatom,1:3)
        END DO
      END IF

      DEALLOCATE (atom_of_kind,STAT=stat)
      CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)

      DEALLOCATE (kind_of,STAT=stat)
      CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)

      DEALLOCATE (natom_of_kind,STAT=stat)
      CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)

    END IF

  END SUBROUTINE harris_force_test_rho_elec

! *****************************************************************************
!> \param qs_env The QS environment of matter
!> \param harris_force The harris force type in which all necessary informations
!>                      are stored
!> \param error variable to control error logging, stopping,...
!>               see module cp_error_handling
!> \note
!>      Only for internal debugging purposes
!> \par History
!>      06.2005 initial create [tdk]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! *****************************************************************************
  SUBROUTINE harris_force_test_integral_vxc(qs_env, harris_force, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(harris_force_type), POINTER         :: harris_force
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: &
      routineN = 'harris_force_test_integral_vxc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, iatom, ikind, ispin, &
                                                natom, nkind, nspins, stat, &
                                                unit_nr
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: atom_of_kind, kind_of, &
                                                natom_of_kind
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: Ehartree, Exc
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, rho_ao
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type)                          :: rho_tot_gspace, &
                                                v_hartree_gspace, &
                                                v_hartree_rspace
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: v_rspace_new, v_tau_rspace
    TYPE(pw_p_type), POINTER                 :: rho_core
    TYPE(pw_poisson_type), POINTER           :: poisson_env
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force, force_tmp
    TYPE(qs_rho_type), POINTER               :: rho
    TYPE(section_vals_type), POINTER         :: input, xc_section

!   ------------------------------------------------------------------------

    failure = .FALSE.
    NULLIFY(rho, v_rspace_new, matrix_ks, auxbas_pw_pool, pw_env, &
            rho_core, force_tmp, v_tau_rspace, poisson_env, dft_control, input,&
            xc_section)
    logger => cp_error_get_logger(error)

    CPPrecondition(ASSOCIATED(qs_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(qs_env%ref_count>0, cp_failure_level, routineP, error, failure)
    CPPrecondition(ASSOCIATED(harris_force), cp_failure_level, routineP, error, failure)
    CPPrecondition(harris_force%ref_count>0, cp_failure_level, routineP, error, failure)

    unit_nr=cp_logger_get_default_io_unit(logger)
    IF (unit_nr>0) THEN
          WRITE (unit_nr,*) ""; WRITE (unit_nr, *) ""
          WRITE (unit_nr,*) "The Harris force correction integral_vxc test is performed!"
          WRITE (unit_nr,*) ""
    END IF

    IF (.NOT. failure) THEN
      CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set, force=force, &
                      dft_control=dft_control, input=input, particle_set=particle_set, error=error)

      xc_section => section_vals_get_subs_vals(input, "DFT%XC", error=error)

      natom = SIZE(particle_set)
      nkind = SIZE(force)
      nspins = dft_control%nspins

      ALLOCATE (atom_of_kind(natom),STAT=stat)
      CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

      ALLOCATE (kind_of(natom),STAT=stat)
      CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

      ALLOCATE (natom_of_kind(nkind),STAT=stat)
      CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                               atom_of_kind=atom_of_kind, kind_of=kind_of, &
                               natom_of_kind=natom_of_kind)

      ! *** Debug Code: Calculate the forces of the sum of eigenvalues *** !
      CALL duplicate_qs_force(qs_force_input=force, qs_force_output=force_tmp, &
                              natom_of_kind=natom_of_kind)
      CALL zero_qs_force(force)

      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env, rho=rho, rho_core=rho_core, &
                      matrix_ks=matrix_ks, error=error)

      CALL pw_env_get(pw_env=pw_env, auxbas_pw_pool=auxbas_pw_pool, &
                      poisson_env=poisson_env, error=error)
      CALL pw_pool_create_pw(auxbas_pw_pool, v_hartree_gspace%pw, &
                              use_data = COMPLEXDATA1D,&
                              in_space = RECIPROCALSPACE, error=error)
      CALL pw_pool_create_pw(auxbas_pw_pool, rho_tot_gspace%pw, &
                              use_data=COMPLEXDATA1D, &
                              in_space=RECIPROCALSPACE, error=error)
      CALL pw_copy(rho_core%pw,rho_tot_gspace%pw, error=error)
      DO ispin = 1,nspins
        CALL pw_axpy(rho%rho_g(ispin)%pw,rho_tot_gspace%pw,error=error)
        ! old code seems to be wrong for spin polarized case
        ! CALL pw_add(rho_core%pw, rho%rho_g(ispin)%pw, rho_tot_gspace%pw)
      END DO

      CALL pw_poisson_solve(poisson_env, rho_tot_gspace%pw, ehartree=Ehartree, &
                            vhartree=v_hartree_gspace%pw,error=error)

      CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_tot_gspace%pw, &
                                   error=error)

      CALL pw_pool_create_pw(auxbas_pw_pool, v_hartree_rspace%pw, &
                              use_data=REALDATA3D, in_space=REALSPACE,error=error)
      CALL pw_transfer(v_hartree_gspace%pw, v_hartree_rspace%pw, error=error)

      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_gspace%pw, &
                                   error=error)

      NULLIFY(v_rspace_new)
      CALL qs_vxc_create(qs_env=qs_env, rho_struct=qs_env%rho, xc_section=xc_section, &
                         vxc_rho=v_rspace_new, vxc_tau=v_tau_rspace, exc=Exc, &
                         just_energy = .FALSE., error=error)

      ! sum-up: v_hartree_rspace = v_hartree_rspace + v_rspace_new(ispin)
      DO ispin = 1,nspins
        CALL pw_axpy(v_rspace_new(ispin)%pw, v_hartree_rspace%pw,error=error)
        CALL pw_scale(v_hartree_rspace%pw, v_hartree_rspace%pw%pw_grid%dvol, error=error)
      END DO

      IF (ASSOCIATED(v_rspace_new)) THEN
        DO ispin = 1,nspins
          CALL pw_pool_give_back_pw(auxbas_pw_pool, v_rspace_new(ispin)%pw, error=error)
        END DO
        DEALLOCATE(v_rspace_new,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      END IF
      IF (ASSOCIATED(v_tau_rspace)) THEN
        DO ispin = 1,nspins
          CALL pw_pool_give_back_pw(auxbas_pw_pool, v_tau_rspace(ispin)%pw, error=error)
        END DO
        DEALLOCATE(v_tau_rspace,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      END IF

    rho_ao => rho%rho_ao
      DO ispin = 1,nspins
        CALL integrate_v_rspace(v_rspace=v_hartree_rspace, p=rho_ao(ispin), &
                                h=matrix_ks(ispin), qs_env=qs_env, &
                                calculate_forces=.TRUE., gapw=.FALSE.,error=error)
      END DO

      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_rspace%pw, error=error)

      DO iatom = 1, natom
        ikind = kind_of(iatom)
        i = atom_of_kind(iatom)
        WRITE (unit_nr,*) force(ikind)%rho_elec(1:3,i)
      END DO

      CALL deallocate_qs_force(force)
      CALL set_qs_env(qs_env=qs_env, force=force_tmp,error=error)
      CALL get_qs_env(qs_env=qs_env, force=force, error=error)

      IF (unit_nr>0) THEN
        WRITE (unit_nr, *) ""
        DO iatom = 1, natom
          ikind = kind_of(iatom)
          i = atom_of_kind(iatom)
          WRITE (unit_nr, *) force(ikind)%rho_elec(1:3,i)
          WRITE (unit_nr, *) harris_force%f_cross_integrate_v(iatom,1:3) &
                           + harris_force%f_hartree(iatom,1:3) &
                           + harris_force%f_xc(iatom,1:3)
        END DO
        WRITE (unit_nr, *) ""
        DO iatom = 1, natom
          ikind = kind_of(iatom)
          i = atom_of_kind(iatom)
          WRITE (unit_nr, *) harris_force%f_cross_integrate_v(iatom,1:3)
          WRITE (unit_nr, *) harris_force%f_cross_integrate_v_core(iatom,1:3)
          WRITE (unit_nr, *) harris_force%f_V(iatom,1:3)
        END DO
        WRITE (unit_nr, *) ""
        DO  iatom = 1, natom
          ikind = kind_of(iatom)
          i = atom_of_kind(iatom)
          WRITE (unit_nr,*) harris_force%f_trace(iatom,1:3)
          WRITE (unit_nr,*) force(ikind)%kinetic(1:3,i) &
               + force(ikind)%rho_elec(1:3,i) &
               + harris_force%f_delta_integral_vxc(iatom,1:3) &
               + harris_force%f_hartree(iatom,1:3) &
               + harris_force%f_cross_integrate_v_core(iatom,1:3) &
               + force(ikind)%overlap(1:3,i) &
               + force(ikind)%gth_ppl(1:3,i) &
               + force(ikind)%gth_nlcc(1:3,i) &
               + force(ikind)%gth_ppnl(1:3,i)
        END DO

      END IF

      DEALLOCATE (atom_of_kind,STAT=stat)
      CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)

      DEALLOCATE (kind_of,STAT=stat)
      CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)

      DEALLOCATE (natom_of_kind,STAT=stat)
      CPPostconditionNoFail(stat==0, cp_warning_level, routineP, error)

    END IF

  END SUBROUTINE harris_force_test_integral_vxc

END MODULE harris_force
