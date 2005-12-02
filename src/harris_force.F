!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000-2005  CP2K developers group                            !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/harris_functional [1.0] *
!!
!!   NAME
!!     harris_functional
!!
!!   FUNCTION
!!     Performs the Harris functional force correction
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     06.2005 initial create [tdk]
!!
!!   SOURCE
!****************************************************************************
MODULE harris_force

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind_set
  USE coefficient_types,               ONLY: coeff_type, &
                                             coeff_transform_space,&
                                             coeff_scale, &
                                             coeff_add, &
                                             coeff_sumup
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_error_handling,               ONLY: cp_a_l, &
                                             cp_assert, &
                                             cp_assertion_failed, &
                                             cp_debug, &
                                             cp_error_get_logger, &
                                             cp_error_message, &
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level, &
                                             cp_note_level, &
                                             cp_to_string, &
                                             cp_warning_level,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type
  USE global_types,                    ONLY: global_environment_type
  USE green_types,                     ONLY: greens_function_type
  USE harris_env_types,                ONLY: harris_env_type, &
                                             harris_env_get
  USE harris_force_types,              ONLY: harris_force_type
  USE hartree,                         ONLY: calculate_hartree
  USE kinds,                           ONLY: dp
  USE message_passing,                 ONLY: mp_sum
  USE particle_types,                  ONLY: particle_type
  USE pw_env_types,                    ONLY: pw_env_type, &
                                             pw_env_get
  USE pw_pool_types,                   ONLY: pw_pool_type, &
                                             pw_pool_init_coeff, &
                                             pw_pool_give_back_coeff, &
                                             pw_pool_create_pw
  USE pws,                             ONLY: cleanup_pw_poisson_solver, &
                                             init_pw_poisson_solver
  USE pw_types,                        ONLY: pw_p_type, &
                                             COMPLEXDATA1D, &
                                             RECIPROCALSPACE, &
                                             REALDATA3D, &
                                             REALSPACE, &
                                             pw_zero, &
                                             pw_release
  !USE qs_core_energies,                ONLY: calculate_ecore_overlap
  USE qs_core_hamiltonian,             ONLY: build_core_hamiltonian_matrix
  USE qs_environment_types,            ONLY: get_qs_env, &
                                             set_qs_env, &
                                             qs_environment_type
  USE qs_force_types,                  ONLY: qs_force_type, &
                                             deallocate_qs_force, &
                                             zero_qs_force, &
                                             duplicate_qs_force
  USE qs_integrate_potential,          ONLY: integrate_v_rspace, &
                                             integrate_v_core_rspace
  USE qs_ks_methods,                   ONLY: qs_vxc_create
  USE qs_rho_methods,                  ONLY: diff_rho_type
  USE qs_rho_types,                    ONLY: qs_rho_type, &
                                             qs_rho_get
  USE simulation_cell,                 ONLY: cell_type
  USE sparse_matrix_types,             ONLY: real_matrix_p_type, &
                                             transfer_matrix
  USE timings,                         ONLY: timeset, &
                                             timestop
  USE xc,                              ONLY: xc_prep_2nd_deriv_new, &
                                             xc_calc_2nd_deriv_new
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type, &
                                             xc_dset_release
  USE xc_rho_set_types,                ONLY: xc_rho_set_type, &
                                             xc_rho_set_release

  IMPLICIT NONE
  PRIVATE
  
  CHARACTER(len=*), PARAMETER :: moduleN = 'harris_force'

  ! *** Public subroutines ***
  PUBLIC :: harris_force_correction, &
            harris_force_EVal, &
            harris_calc_nsc_force

!***
!****************************************************************************

CONTAINS

!!****f* harris_force/harris_force_correction [1.0] *
!!
!!   NAME
!!     harris_force_correction
!!
!!   FUNCTION
!!     Performs the harris functional force correction
!!
!!   NOTES
!!     Performs the force calculation according to the Harris energy functional
!!
!!   INPUTS
!!     - qs_env: The QS environment of matter
!!     - harris_force: The harris force type in which all necessary informations 
!!                     are stored
!!     - globenv: The global environment from which the para_env is extracted
!!     - fast: Suppress some double calculations among the Harris force terms.
!!     - error: variable to control error logging, stopping,... 
!!              see module cp_error_handling 
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     06.2005 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE harris_force_correction(qs_env, harris_env, globenv, fast, error)

    TYPE(qs_environment_type), POINTER              :: qs_env
    TYPE(harris_force_type),  POINTER               :: harris_force
    TYPE(global_environment_type), POINTER          :: globenv
    LOGICAL, INTENT(IN), OPTIONAL                   :: fast
    TYPE(cp_error_type), INTENT(INOUT), OPTIONAL    :: error
    TYPE(harris_env_type), POINTER                  :: harris_env
    TYPE(qs_force_type), DIMENSION(:), POINTER      :: force, &
                                                       force_tmp
    INTEGER                                         :: natom, &
                                                       handle, &
                                                       nkind, &
                                                       i, &
                                                       iatom, &
                                                       ikind, &
                                                       ispin, &
                                                       stat, &
                                                       nspins
    INTEGER, ALLOCATABLE, DIMENSION(:)              :: atom_of_kind, &
                                                       natom_of_kind, &
                                                       kind_of
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                       :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                       :: particle_set
    TYPE(coeff_type), DIMENSION(:), POINTER         :: v_rspace_new, &
                                                       v_tau_rspace
    REAL(KIND=dp)                                   :: Exc, &
                                                       Ehartree
                                                       !E_overlap_core
    TYPE(qs_rho_type), POINTER                      :: rho
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: matrix_ks, &
                                                       matrix_h
    TYPE(greens_function_type), POINTER             :: green
    TYPE(pw_env_type), POINTER                      :: pw_env
    TYPE(pw_pool_type), POINTER                     :: auxbas_pw_pool
    TYPE(coeff_type)                                :: v_hartree_gspace, &
                                                       v_hartree_rspace
    TYPE(coeff_type), POINTER                       :: rho_core
    TYPE(cell_type), POINTER                        :: cell
    TYPE(dft_control_type), POINTER                 :: dft_control
    TYPE(xc_derivative_set_type), POINTER           :: deriv_set
    TYPE(xc_rho_set_type), POINTER                  :: rho_set
    TYPE(pw_p_type), DIMENSION(:), POINTER          :: my_rho_r, &
                                                       v_xc
    TYPE(coeff_type), DIMENSION(:), POINTER         :: rho_r
    LOGICAL                                         :: failure, &
                                                       fast_flag

    CHARACTER(len=*), PARAMETER :: routineN = 'harris_force_correction', &
      routineP = moduleN//':'//routineN
    TYPE(cp_logger_type), POINTER :: logger
    INTEGER :: unit_nr

!   ------------------------------------------------------------------------

    CALL timeset(routineN, "I", "", handle)

    failure = .FALSE.
    NULLIFY(force, force_tmp, v_rspace_new, v_tau_rspace, rho, matrix_ks, &
            green, pw_env, auxbas_pw_pool, rho_core, cell, dft_control, &
            atomic_kind_set, particle_set, deriv_set, rho_set, my_rho_r, &
            v_xc, matrix_h, harris_force, rho_r)
    logger => cp_error_get_logger(error)

    IF (PRESENT(fast)) THEN
      fast_flag = fast
    ELSE
      fast_flag = .FALSE.
    END IF

    CPPrecondition(ASSOCIATED(qs_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(qs_env%ref_count>0, cp_failure_level, routineP, error, failure)
    CPPrecondition(ASSOCIATED(harris_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(harris_env%ref_count>0, cp_failure_level, routineP, error, failure)
    CPPrecondition(ASSOCIATED(globenv), cp_failure_level, routineP, error, failure)
    CPPrecondition(globenv%ref_count>0, cp_failure_level, routineP, error, failure)

    IF (.NOT. failure) THEN
      CALL harris_env_get(harris_env=harris_env, harris_force=harris_force, error=error)
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

      CALL duplicate_qs_force(qs_force_input=force, qs_force_output=force_tmp, &
                              natom_of_kind=natom_of_kind)
      CALL zero_qs_force(force)
      

      ! *** d/dR[Sum of eigenvalues] = d/dR[trace(H_in * rho_out)] *** !
      CALL build_core_hamiltonian_matrix(qs_env=qs_env, globenv=globenv, &
                                         calculate_forces=.TRUE.)

      DO ikind = 1,SIZE(force)
        CALL mp_sum(force(ikind)%kinetic, globenv%group)
        CALL mp_sum(force(ikind)%overlap, globenv%group)
        CALL mp_sum(force(ikind)%gth_ppl, globenv%group)
        CALL mp_sum(force(ikind)%gth_ppnl, globenv%group)
      END DO

      DO iatom = 1, natom
        ikind = kind_of(iatom)
        i = atom_of_kind(iatom)
        harris_force%f_kinetic(iatom, 1:3) = force(ikind)%kinetic(1:3,i)
        harris_force%f_overlap(iatom, 1:3) = force(ikind)%overlap(1:3,i)
        harris_force%f_gth_pp(iatom, 1:3)  = force(ikind)%gth_ppl(1:3,i) &
                                           + force(ikind)%gth_ppnl(1:3,i)
      END DO

      CALL zero_qs_force(force)

      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env, rho=rho, rho_core=rho_core, &
                      green=green, matrix_h=matrix_h, matrix_ks=matrix_ks, error=error)

      CALL pw_env_get(pw_env=pw_env, auxbas_pw_pool=auxbas_pw_pool, error=error)
      CALL pw_pool_init_coeff(auxbas_pw_pool, v_hartree_gspace, &
                              use_data = COMPLEXDATA1D, &
                              in_space = RECIPROCALSPACE, error=error)
      CALL pw_pool_init_coeff(auxbas_pw_pool, v_hartree_rspace, &
                              use_data=REALDATA3D, in_space=REALSPACE)

      ! *** Calculation of the forces due to the Hartree energy *** !
      CALL init_pw_poisson_solver(green,error=error)
      CALL calculate_hartree(qs_env%rho%rho_g(1), ehartree=Ehartree, vhartree=v_hartree_gspace)
      CALL cleanup_pw_poisson_solver(green,error=error)

      CALL coeff_transform_space(v_hartree_gspace, v_hartree_rspace)
      CALL coeff_scale(v_hartree_rspace, v_hartree_rspace%pw%pw_grid%dvol)

      CALL integrate_v_rspace(v_rspace=v_hartree_rspace, p=harris_env%rho%rho_ao(1), &
                              h=matrix_ks(1), qs_env=qs_env, calculate_forces=.TRUE., &
                              gapw=.FALSE., error=error)

      DO ikind = 1,SIZE(force)
        CALL mp_sum(force(ikind)%rho_elec, globenv%group)
      END DO

      DO iatom = 1, natom
        ikind = kind_of(iatom)
        i = atom_of_kind(iatom)
        harris_force%f_hartree(iatom, 1:3) = force(ikind)%rho_elec(1:3,i)
      END DO

      CALL zero_qs_force(force)

      CALL init_pw_poisson_solver(green,error=error)
      CALL calculate_hartree(harris_env%rho%rho_g(1), ehartree=Ehartree, &
                             vhartree=v_hartree_gspace)
      CALL cleanup_pw_poisson_solver(green,error=error)

      CALL coeff_transform_space(v_hartree_gspace, v_hartree_rspace)
      CALL coeff_scale(v_hartree_rspace, v_hartree_rspace%pw%pw_grid%dvol)

      CALL integrate_v_rspace(v_rspace=v_hartree_rspace, p=qs_env%rho%rho_ao(1), &
                              h=matrix_ks(1), qs_env=qs_env, calculate_forces=.TRUE., &
                              gapw=.FALSE.,error=error)

      DO ikind = 1,SIZE(force)
        CALL mp_sum(force(ikind)%rho_elec, globenv%group)
      END DO

      DO iatom = 1, natom
        ikind = kind_of(iatom)
        i = atom_of_kind(iatom)
        harris_force%f_hartree(iatom, 1:3) = harris_force%f_hartree(iatom, 1:3) + &
             force(ikind)%rho_elec(1:3,i)
      END DO

      CALL zero_qs_force(force)

      ! ** Force computation of f_V  ** !
      CALL init_pw_poisson_solver(green, error=error)
      CALL calculate_hartree(qs_env%rho%rho_g(1), ehartree=Ehartree, vhartree=v_hartree_gspace)
      CALL cleanup_pw_poisson_solver(green, error=error)

      CALL coeff_transform_space(v_hartree_gspace, v_hartree_rspace)
      CALL coeff_scale(v_hartree_rspace, v_hartree_rspace%pw%pw_grid%dvol)

      CALL integrate_v_core_rspace(v_hartree_rspace, qs_env)

      DO ikind = 1,SIZE(force)
        CALL mp_sum(force(ikind)%rho_core, globenv%group)
      END DO
      
      DO iatom = 1, natom
        ikind = kind_of(iatom)
        i = atom_of_kind(iatom)
        harris_force%f_cross_integrate_v_core(iatom, 1:3) = force(ikind)%rho_core(1:3,i)
      END DO

      CALL zero_qs_force(force)

      CALL transfer_matrix(matrix_h(1)%matrix, matrix_ks(1)%matrix)
            
      CALL init_pw_poisson_solver(green,error=error)
      CALL calculate_hartree(rho_core, ehartree=Ehartree, &
                             vhartree=v_hartree_gspace)
      CALL cleanup_pw_poisson_solver(green,error=error)

      CALL coeff_transform_space(v_hartree_gspace, v_hartree_rspace)
      CALL coeff_scale(v_hartree_rspace, v_hartree_rspace%pw%pw_grid%dvol)

      CALL integrate_v_rspace(v_rspace=v_hartree_rspace, p=qs_env%rho%rho_ao(1), &
                              h=matrix_ks(1), qs_env=qs_env, calculate_forces=.TRUE., &
                              gapw=.FALSE.,error=error)

      CALL pw_pool_give_back_coeff(auxbas_pw_pool, v_hartree_gspace, &
                                   error=error)
      CALL pw_pool_give_back_coeff(auxbas_pw_pool, v_hartree_rspace, error=error)

      DO ikind = 1,SIZE(force)
        CALL mp_sum(force(ikind)%rho_elec, globenv%group)
      END DO

      DO iatom = 1, natom
        ikind = kind_of(iatom)
        i = atom_of_kind(iatom)
        harris_force%f_cross_integrate_v(iatom, 1:3) = force(ikind)%rho_elec(1:3,i)
        harris_force%f_V(iatom,1:3) = harris_force%f_cross_integrate_v(iatom,1:3) &
                                    + harris_force%f_cross_integrate_v_core(iatom,1:3)
      END DO

      CALL zero_qs_force(force)

      ! *** Calculation of the forces of the xc-integral *** !
      CALL qs_vxc_create(qs_env=qs_env, vxc_rho=v_rspace_new, vxc_tau=v_tau_rspace, &
                         exc=Exc, just_energy=.FALSE., error=error)

      CALL coeff_scale(v_rspace_new(1), v_rspace_new(1)%pw%pw_grid%dvol)
      CALL integrate_v_rspace(v_rspace=v_rspace_new(1), p=qs_env%rho%rho_ao(1), &
                              h=matrix_ks(1), qs_env=qs_env, calculate_forces=.TRUE., &
                              gapw=.FALSE., error=error)

      DO ikind = 1,SIZE(force)
        CALL mp_sum(force(ikind)%rho_elec, globenv%group)
      END DO

      DO iatom = 1, natom
        ikind = kind_of(iatom)
        i = atom_of_kind(iatom)
        harris_force%f_xc(iatom, 1:3) = force(ikind)%rho_elec(1:3,i)
      END DO

      CALL zero_qs_force(force)

      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, cell=cell, error=error)

      CALL qs_rho_get(rho_struct=rho, rho_r=rho_r, error=error)
      ALLOCATE(my_rho_r(SIZE(rho_r)), stat=stat)
      CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
      DO ispin=1, SIZE(rho_r)
        my_rho_r(ispin)%pw => harris_env%rho%rho_r(ispin)%pw
      END DO
      CALL xc_prep_2nd_deriv_new(deriv_set=deriv_set, rho_set=rho_set, &
                                 rho_r=my_rho_r, pw_pool=auxbas_pw_pool, &
                                 dft_control=dft_control, cell=cell, error=error)
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
        CALL pw_zero(v_xc(ispin)%pw)
      END DO

      CALL xc_calc_2nd_deriv_new(v_xc, deriv_set=deriv_set, rho_set=rho_set, &
                                rho1_set=rho_set, pw_pool=auxbas_pw_pool, &
                                dft_control=dft_control, gapw=.FALSE., error=error)

      CALL xc_dset_release(deriv_set, error=error)

      v_xc(1)%pw%cr3d = v_xc(1)%pw%cr3d * v_xc(1)%pw%pw_grid%dvol
      v_rspace_new(1)%pw%cr3d = rho_r(1)%pw%cr3d * v_xc(1)%pw%cr3d
      CALL coeff_scale(v_rspace_new(1), v_rspace_new(1)%pw%pw_grid%dvol)

      CALL pw_release(v_xc(1)%pw, error=error)
      DEALLOCATE(v_xc,stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      
      CALL xc_rho_set_release(rho_set)

      CALL integrate_v_rspace(v_rspace=v_rspace_new(1), p=harris_env%rho%rho_ao(1), &
                              h=matrix_ks(1), qs_env=qs_env, calculate_forces=.TRUE., &
                              gapw=.FALSE., error=error)

      IF (ASSOCIATED(v_rspace_new)) THEN
        CALL pw_pool_give_back_coeff(auxbas_pw_pool,v_rspace_new(1),error=error)
        DEALLOCATE(v_rspace_new,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      END IF
      IF (ASSOCIATED(v_tau_rspace)) THEN
        CALL pw_pool_give_back_coeff(auxbas_pw_pool,v_tau_rspace(1),error=error)
        DEALLOCATE(v_tau_rspace,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      END IF

      DO ikind = 1,SIZE(force)
        CALL mp_sum(force(ikind)%rho_elec, globenv%group)
      END DO

      DO iatom = 1, natom
        ikind = kind_of(iatom)
        i = atom_of_kind(iatom)
        harris_force%f_delta_integral_vxc(iatom, 1:3) = force(ikind)%rho_elec(1:3,i)
        harris_force%f_integral_vxc(iatom, 1:3) = harris_force%f_xc(iatom,1:3) &
             + harris_force%f_delta_integral_vxc(iatom,1:3)
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
        harris_force%f_kinetic(iatom, 1:3) = 0.0_dp
        harris_force%f_gth_pp(iatom, 1:3) = 0.0_dp
        harris_force%f_overlap(iatom, 1:3) = 0.0_dp
        harris_force%f_V(iatom, 1:3) = 0.0_dp
        harris_force%f_integral_vxc(iatom, 1:3) = 0.0_dp
        harris_force%f_xc(iatom, 1:3) = 0.0_dp
        harris_force%f_delta_integral_vxc(iatom, 1:3) = 0.0_dp
        harris_force%f_cross_integrate_v(iatom, 1:3) = 0.0_dp
        harris_force%f_cross_integrate_v_core(iatom, 1:3) = 0.0_dp
      END DO

      !CALL zero_qs_force(force)
      !
      !CALL calculate_ecore_overlap(qs_env=qs_env, globenv=globenv, &
      !                             calculate_forces=.TRUE., E_overlap_core=E_overlap_core)
      !
      !DO iatom = 1, natom
      !  ikind = kind_of(iatom)
      !  i = atom_of_kind(iatom)
      !  harris_force%f_ovrl(iatom, 1:3) = force(ikind)%core_overlap(1:3,i)
      !END DO

      CALL zero_qs_force(force)

      NULLIFY(rho, rho_core, pw_env, green, auxbas_pw_pool, cell, dft_control, &
              rho_r, my_rho_r, deriv_set, rho_set, v_xc, matrix_ks, matrix_h, &
              v_rspace_new, v_tau_rspace)
      
      ! *** End of force calculation due to the sum_of_eigenvalues *** !


      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env, rho_core=rho_core, &
                      green=green, error=error)
      CALL harris_env_get(harris_env=harris_env, rho=rho, error=error)

      CALL pw_env_get(pw_env=pw_env, auxbas_pw_pool=auxbas_pw_pool, error=error)


      ! *** The forces due to hartree core-core repulsion *** !
      CALL pw_pool_init_coeff(auxbas_pw_pool, v_hartree_gspace, &
                              use_data = COMPLEXDATA1D,&
                              in_space = RECIPROCALSPACE, error=error)
      !CALL pw_pool_init_coeff(auxbas_pw_pool, rho_tot_gspace, &
      !                        use_data=COMPLEXDATA1D, &
      !                        in_space=RECIPROCALSPACE, error=error)
      
      !CALL coeff_add(rho_core, rho%rho_g(1), rho_tot_gspace)

      CALL init_pw_poisson_solver(green,error=error)
      CALL calculate_hartree(rho_core, ehartree=Ehartree, vhartree=v_hartree_gspace)
      !CALL calculate_hartree(rho_tot_gspace, ehartree=Ehartree, vhartree=v_hartree_gspace)
      CALL cleanup_pw_poisson_solver(green,error=error)

      !CALL pw_pool_give_back_coeff(auxbas_pw_pool, rho_tot_gspace, &
      !                             error=error)

      CALL pw_pool_init_coeff(auxbas_pw_pool, v_hartree_rspace, &
                              use_data=REALDATA3D, in_space=REALSPACE)
      CALL coeff_transform_space(v_hartree_gspace, v_hartree_rspace)
      CALL coeff_scale(v_hartree_rspace, v_hartree_rspace%pw%pw_grid%dvol)

      !CALL pw_pool_give_back_coeff(auxbas_pw_pool, v_hartree_gspace, &
      !                             error=error)
      !CALL integrate_v_rspace(v_rspace=v_hartree_rspace, p=rho%rho_ao(1), &
      !                        h=matrix_ks(1), qs_env=qs_env, calculate_forces=.TRUE., &
      !                        gapw=.FALSE.,error=error)
      CALL integrate_v_core_rspace(v_hartree_rspace, qs_env)
      !CALL pw_pool_give_back_coeff(auxbas_pw_pool, v_hartree_rspace, error=error)

      DO ikind = 1,SIZE(force)
        CALL mp_sum(force(ikind)%rho_core, globenv%group)
      END DO
      
      DO iatom = 1, natom
        ikind = kind_of(iatom)
        i = atom_of_kind(iatom)
        harris_force%f_rho_core(iatom, 1:3) = force(ikind)%rho_core(1:3,i)
      END DO

      CALL zero_qs_force(force)

      
      !CALL calculate_ecore_overlap(qs_env=qs_env, globenv=globenv, &
      !                             calculate_forces=.TRUE., E_overlap_core=E_overlap_core)
      !
      !DO iatom = 1, natom
      !  ikind = kind_of(iatom)
      !  i = atom_of_kind(iatom)
      !  harris_force%f_ovrl(iatom, 1:3) = force(ikind)%core_overlap(1:3,i)
      !END DO
      !
      !CALL zero_qs_force(force)

      CALL get_qs_env(qs_env=qs_env, matrix_h=matrix_h, matrix_ks=matrix_ks, error=error)

      ! *** The forces due to the hartree electron-electron repulsion *** !
      !CALL pw_pool_init_coeff(auxbas_pw_pool, v_hartree_gspace, &
      !                        use_data = COMPLEXDATA1D,&
      !                        in_space = RECIPROCALSPACE, error=error)
      !CALL pw_pool_init_coeff(auxbas_pw_pool, rho_tot_gspace, &
      !                        use_data=COMPLEXDATA1D, &
      !                        in_space=RECIPROCALSPACE, error=error)

      !CALL coeff_add(rho_core, rho%rho_g(1), rho_tot_gspace)

      CALL init_pw_poisson_solver(green,error=error)
      !CALL calculate_hartree(rho_tot_gspace, ehartree=Ehartree, &
      !                       vhartree=v_hartree_gspace)
      CALL calculate_hartree(rho%rho_g(1), ehartree=Ehartree, vhartree=v_hartree_gspace)
      CALL cleanup_pw_poisson_solver(green,error=error)

      !CALL pw_pool_give_back_coeff(auxbas_pw_pool, rho_tot_gspace, &
      !                             error=error)

      !CALL pw_pool_init_coeff(auxbas_pw_pool, v_hartree_rspace, &
      !                        use_data=REALDATA3D, in_space=REALSPACE)
      CALL coeff_transform_space(v_hartree_gspace, v_hartree_rspace)
      CALL coeff_scale(v_hartree_rspace, v_hartree_rspace%pw%pw_grid%dvol)

      CALL pw_pool_give_back_coeff(auxbas_pw_pool, v_hartree_gspace, &
                                   error=error)
      CALL integrate_v_rspace(v_rspace=v_hartree_rspace, p=rho%rho_ao(1), &
                              h=matrix_ks(1), qs_env=qs_env, calculate_forces=.TRUE., &
                              gapw=.FALSE.,error=error)
      CALL pw_pool_give_back_coeff(auxbas_pw_pool, v_hartree_rspace, error=error)

      DO ikind = 1,SIZE(force)
        CALL mp_sum(force(ikind)%rho_elec, globenv%group)
      END DO

      DO iatom = 1, natom
        ikind = kind_of(iatom)
        i = atom_of_kind(iatom)
        harris_force%f_hartree(iatom, 1:3) = force(ikind)%rho_elec(1:3,i)
      END DO

      CALL zero_qs_force(force)


      ! *** The forces due to exchange and correlation *** !
      CALL qs_vxc_create(qs_env=qs_env, vxc_rho=v_rspace_new, vxc_tau=v_tau_rspace, &
                         exc=Exc, just_energy=.FALSE., harris=.FALSE., error=error)
      !CALL qs_vxc_create(qs_env=qs_env, vxc_rho=v_rspace_new, vxc_tau=v_tau_rspace, &
      !                   exc=Exc, just_energy=.FALSE., harris=.TRUE., error=error)
      IF (.NOT. fast_flag) THEN
        CALL coeff_scale(v_rspace_new(1), v_rspace_new(1)%pw%pw_grid%dvol)
        CALL integrate_v_rspace(v_rspace=v_rspace_new(1), p=rho%rho_ao(1), &
                                h=matrix_ks(1), qs_env=qs_env, calculate_forces=.TRUE., &
                                gapw=.FALSE., error=error)

        DO ikind = 1,SIZE(force)
          CALL mp_sum(force(ikind)%rho_elec, globenv%group)
        END DO

        DO iatom = 1, natom
          ikind = kind_of(iatom)
          i = atom_of_kind(iatom)
          harris_force%f_xc(iatom, 1:3) = force(ikind)%rho_elec(1:3,i)
        END DO

        CALL zero_qs_force(force)
      END IF

      ! *** Calculation of the forces of the xc-integral *** !
      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, cell=cell, error=error)

      CALL qs_rho_get(rho, rho_r=rho_r, error=error)
      ALLOCATE(my_rho_r(SIZE(rho_r)), stat=stat)
      CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
      DO ispin=1, SIZE(rho_r)
        my_rho_r(ispin)%pw => harris_env%rho%rho_r(ispin)%pw
      END DO
      CALL xc_prep_2nd_deriv_new(deriv_set=deriv_set, rho_set=rho_set, &
                                 rho_r=my_rho_r, pw_pool=auxbas_pw_pool, &
                                 dft_control=dft_control, cell=cell, error=error)
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
        CALL pw_zero(v_xc(ispin)%pw)
      END DO

      CALL xc_calc_2nd_deriv_new(v_xc, deriv_set=deriv_set, rho_set=rho_set, &
                                rho1_set=rho_set, pw_pool=auxbas_pw_pool, &
                                dft_control=dft_control, gapw=.FALSE., error=error)

      CALL xc_dset_release(deriv_set, error=error)

      !CALL coeff_init(v_rspace_new(1), pw=v_xc(1)%pw, error=error)
      !v_rspace_new(1)%pw%cr3d = v_rspace_new(1)%pw%cr3d * &
      !                          v_rspace_new(1)%pw%pw_grid%dvol
      !v_rspace_new(1)%pw%cr3d = v_rsapce_new(1)%pw%cr3d * rho_r(1)%pw%cr3d

      v_xc(1)%pw%cr3d = v_xc(1)%pw%cr3d * v_xc(1)%pw%pw_grid%dvol
      v_rspace_new(1)%pw%cr3d = rho_r(1)%pw%cr3d * v_xc(1)%pw%cr3d
      CALL coeff_scale(v_rspace_new(1), v_rspace_new(1)%pw%pw_grid%dvol)

      CALL pw_release(v_xc(1)%pw, error=error)
      DEALLOCATE(v_xc,stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      
      CALL xc_rho_set_release(rho_set)

      CALL integrate_v_rspace(v_rspace=v_rspace_new(1), p=harris_env%rho%rho_ao(1), &
                              h=matrix_ks(1), qs_env=qs_env, calculate_forces=.TRUE., &
                              gapw=.FALSE., error=error)

      IF (ASSOCIATED(v_rspace_new)) THEN
        CALL pw_pool_give_back_coeff(auxbas_pw_pool,v_rspace_new(1),error=error)
        DEALLOCATE(v_rspace_new,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      END IF
      IF (ASSOCIATED(v_tau_rspace)) THEN
        CALL pw_pool_give_back_coeff(auxbas_pw_pool,v_tau_rspace(1),error=error)
        DEALLOCATE(v_tau_rspace,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      END IF

      DO ikind = 1,SIZE(force)
        CALL mp_sum(force(ikind)%rho_elec, globenv%group)
      END DO

      DO iatom = 1, natom
        ikind = kind_of(iatom)
        i = atom_of_kind(iatom)
        harris_force%f_delta_integral_vxc(iatom, 1:3) = force(ikind)%rho_elec(1:3,i)
        IF (fast_flag) THEN
          harris_force%f_integral_vxc(iatom, 1:3) = harris_force%f_delta_integral_vxc(iatom,1:3)
        ELSE
          harris_force%f_integral_vxc(iatom, 1:3) = harris_force%f_xc(iatom,1:3) &
               + harris_force%f_delta_integral_vxc(iatom,1:3)
        END IF
      END DO

      !CALL zero_qs_force(force)
      !
      !CALL build_core_hamiltonian_matrix(qs_env=qs_env, globenv=globenv, &
      !                                   calculate_forces=.TRUE., harris=.TRUE.)
      !
      !DO iatom = 1, natom
      !  ikind = kind_of(iatom)
      !  i = atom_of_kind(iatom)
      !  harris_force%f_overlap(iatom, 1:3) = force(ikind)%overlap(1:3,i)
      !  harris_force%f_trace(iatom, 1:3) = harris_force%f_trace(iatom, 1:3) &
      !                                   + force(ikind)%overlap(1:3, i)
      !END DO

      !CALL zero_qs_force(force)
      !
      !CALL calculate_ecore_overlap(qs_env=qs_env, globenv=globenv, &
      !                             calculate_forces=.TRUE., E_overlap_core=E_overlap_core)
      !
      !DO iatom = 1, natom
      !  ikind = kind_of(iatom)
      !  i = atom_of_kind(iatom)
      !  harris_force%f_ovrl(iatom, 1:3) = force(ikind)%core_overlap(1:3,i)
      !END DO
      
      CALL deallocate_qs_force(force)
      CALL set_qs_env(qs_env=qs_env, force=force_tmp, error=error)
      
      CALL get_qs_env(qs_env=qs_env, force=force, error=error)
      harris_force%f_total(:) = 0.0_dp

      DO ikind = 1,SIZE(force)
        CALL mp_sum(force(ikind)%core_overlap, globenv%group)
      END DO

      DO iatom = 1, natom
        ikind = kind_of(iatom)
        i = atom_of_kind(iatom)

        harris_force%f_ovrl(iatom, 1:3) = force(ikind)%core_overlap(1:3,i)
        harris_force%f_self(iatom, 1:3) = 0.0_dp
        harris_force%f_EII(iatom, 1:3) = harris_force%f_ovrl(iatom, 1:3) &
                                       + harris_force%f_self(iatom, 1:3) &
                                       + harris_force%f_rho_core(iatom, 1:3)
        IF (fast_flag) THEN
          harris_force%f_harris(iatom, 1:3) = harris_force%f_trace(iatom, 1:3) &
                                            - harris_force%f_hartree(iatom, 1:3) &
                                            - harris_force%f_integral_vxc(iatom, 1:3) &
                                            + harris_force%f_EII(iatom, 1:3)
        ELSE
          harris_force%f_harris(iatom, 1:3) = harris_force%f_trace(iatom, 1:3) &
                                            - harris_force%f_hartree(iatom, 1:3) &
                                            - harris_force%f_integral_vxc(iatom, 1:3) &
                                            + harris_force%f_xc(iatom, 1:3) &
                                            + harris_force%f_EII(iatom, 1:3)
        END IF
        harris_force%f_total(1:3) = harris_force%f_total(1:3) &
                                  + harris_force%f_harris(iatom, 1:3)
      END DO


      ! Output
      IF (logger%para_env%source==logger%para_env%mepos) THEN
        unit_nr=cp_logger_get_default_unit_nr(logger,local=.false.)
        WRITE (unit_nr,*) ""; WRITE (unit_nr, *) ""
        WRITE (unit_nr,*) "The Harris functional force correction is performed!"
        WRITE (unit_nr,*) ""
      
        !WRITE (unit_nr, *) "F_Trace"
        !DO iatom = 1, natom
        !  WRITE (unit_nr,*) harris_force%f_trace(iatom, 1:3)
        !END DO
        !WRITE (unit_nr, *) "F_Hartree"
        !DO iatom = 1, natom
        !  WRITE (unit_nr,*) harris_force%f_hartree(iatom, 1:3)
        !END DO
        !WRITE (unit_nr, *) "F_Integral_vxc"
        !DO iatom = 1, natom
        !  WRITE (unit_nr,*) harris_force%f_integral_vxc(iatom, 1:3)
        !END DO
        !WRITE (unit_nr, *) "F_Exc"
        !DO iatom = 1, natom
        !  WRITE (unit_nr,*) harris_force%f_xc(iatom, 1:3)
        !END DO
        !WRITE (unit_nr, *) "F_ovrl"
        !DO iatom = 1, natom
        !  WRITE (unit_nr,*) harris_force%f_ovrl(iatom, 1:3)
        !END DO
        !WRITE (unit_nr, *) "F_Self"
        !DO iatom = 1, natom
        !  WRITE (unit_nr,*) harris_force%f_self(iatom, 1:3)
        !END DO
        !WRITE (unit_nr, *) "F_Rho_Core"
        !DO iatom = 1, natom
        !  WRITE (unit_nr,*) harris_force%f_rho_core(iatom, 1:3)
        !END DO
        !WRITE (unit_nr, *) "F_EII"
        !DO iatom = 1, natom
        !  WRITE (unit_nr,*) harris_force%f_EII(iatom, 1:3)
        !END DO
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

    CALL timestop(0.0_dp, handle)

  END SUBROUTINE harris_force_correction
!***************************************************************************

!!****f* harris_force/harris_force_EVal [1.0] *
!!
!!   NAME
!!     harris_force_EVal
!!
!!   FUNCTION
!!     Performs the harris functional force correction
!!
!!   NOTES
!!     Performs the force calculation according to the Harris energy functional
!!
!!   INPUTS
!!     - qs_env: The QS environment of matter
!!     - harris_force: The harris force type in which all necessary informations 
!!                     are stored
!!     - globenv: The global environment from which the para_env is extracted
!!     - fast: Suppress some double calculations among the Harris force terms.
!!     - error: variable to control error logging, stopping,... 
!!              see module cp_error_handling 
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     11.2005 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE harris_force_EVal(qs_env, harris_env, globenv, fast, error)

    TYPE(qs_environment_type), POINTER              :: qs_env
    TYPE(harris_force_type),  POINTER               :: harris_force
    TYPE(global_environment_type), POINTER          :: globenv
    LOGICAL, INTENT(IN), OPTIONAL                   :: fast
    TYPE(cp_error_type), INTENT(INOUT), OPTIONAL    :: error
    TYPE(harris_env_type), POINTER                  :: harris_env
    TYPE(qs_force_type), DIMENSION(:), POINTER      :: force, &
                                                       force_tmp
    INTEGER                                         :: natom, &
                                                       handle, &
                                                       nkind, &
                                                       i, &
                                                       iatom, &
                                                       ikind, &
                                                       ispin, &
                                                       stat, &
                                                       nspins
    INTEGER, ALLOCATABLE, DIMENSION(:)              :: atom_of_kind, &
                                                       natom_of_kind, &
                                                       kind_of
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                       :: atomic_kind_set
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                       :: particle_set
    TYPE(coeff_type), DIMENSION(:), POINTER         :: v_rspace_new, &
                                                       v_tau_rspace
    REAL(KIND=dp)                                   :: Exc, &
                                                       Ehartree
                                                       !E_overlap_core
    TYPE(qs_rho_type), POINTER                      :: rho
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: matrix_ks, &
                                                       matrix_h
    TYPE(greens_function_type), POINTER             :: green
    TYPE(pw_env_type), POINTER                      :: pw_env
    TYPE(pw_pool_type), POINTER                     :: auxbas_pw_pool
    TYPE(coeff_type)                                :: v_hartree_gspace, &
                                                       v_hartree_rspace, &
                                                       rho_tot_gspace
    TYPE(coeff_type), POINTER                       :: rho_core
    TYPE(cell_type), POINTER                        :: cell
    TYPE(dft_control_type), POINTER                 :: dft_control
    TYPE(xc_derivative_set_type), POINTER           :: deriv_set
    TYPE(xc_rho_set_type), POINTER                  :: rho_set
    TYPE(pw_p_type), DIMENSION(:), POINTER          :: my_rho_r, &
                                                       v_xc
    TYPE(coeff_type), DIMENSION(:), POINTER         :: rho_r
    LOGICAL                                         :: failure, &
                                                       fast_flag

    CHARACTER(len=*), PARAMETER :: routineN = 'harris_force_EVal', &
      routineP = moduleN//':'//routineN
    TYPE(cp_logger_type), POINTER :: logger
    INTEGER :: unit_nr

!   ------------------------------------------------------------------------

    CALL timeset(routineN, "I", "", handle)

    failure = .FALSE.
    NULLIFY(force, force_tmp, v_rspace_new, v_tau_rspace, rho, matrix_ks, &
            green, pw_env, auxbas_pw_pool, rho_core, cell, dft_control, &
            atomic_kind_set, particle_set, deriv_set, rho_set, my_rho_r, &
            v_xc, matrix_h, harris_force, rho_r)
    logger => cp_error_get_logger(error)

    IF (PRESENT(fast)) THEN
      fast_flag = fast
    ELSE
      fast_flag = .FALSE.
    END IF

    CPPrecondition(ASSOCIATED(qs_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(qs_env%ref_count>0, cp_failure_level, routineP, error, failure)
    CPPrecondition(ASSOCIATED(harris_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(harris_env%ref_count>0, cp_failure_level, routineP, error, failure)
    CPPrecondition(ASSOCIATED(globenv), cp_failure_level, routineP, error, failure)
    CPPrecondition(globenv%ref_count>0, cp_failure_level, routineP, error, failure)

    IF (.NOT. failure) THEN
      CALL harris_env_get(harris_env=harris_env, harris_force=harris_force, error=error)
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

      CALL duplicate_qs_force(qs_force_input=force, qs_force_output=force_tmp, &
                              natom_of_kind=natom_of_kind)
      CALL zero_qs_force(force)
      

      ! *** d/dR[Sum of eigenvalues] = d/dR[trace(H_in * rho_out)] *** !
      CALL build_core_hamiltonian_matrix(qs_env=qs_env, globenv=globenv, &
                                         calculate_forces=.TRUE.)

      DO ikind = 1,SIZE(force)
        CALL mp_sum(force(ikind)%kinetic, globenv%group)
        CALL mp_sum(force(ikind)%overlap, globenv%group)
        CALL mp_sum(force(ikind)%gth_ppl, globenv%group)
        CALL mp_sum(force(ikind)%gth_ppnl, globenv%group)
      END DO

      DO iatom = 1, natom
        ikind = kind_of(iatom)
        i = atom_of_kind(iatom)
        harris_force%f_kinetic(iatom, 1:3) = force(ikind)%kinetic(1:3,i)
        harris_force%f_overlap(iatom, 1:3) = force(ikind)%overlap(1:3,i)
        harris_force%f_gth_pp(iatom, 1:3)  = force(ikind)%gth_ppl(1:3,i) &
                                           + force(ikind)%gth_ppnl(1:3,i)
      END DO

      CALL zero_qs_force(force)

      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env, rho=rho, rho_core=rho_core, &
                      green=green, matrix_h=matrix_h, matrix_ks=matrix_ks, error=error)

      CALL pw_env_get(pw_env=pw_env, auxbas_pw_pool=auxbas_pw_pool, error=error)
      CALL pw_pool_init_coeff(auxbas_pw_pool, v_hartree_gspace, &
                              use_data = COMPLEXDATA1D, &
                              in_space = RECIPROCALSPACE, error=error)
      CALL pw_pool_init_coeff(auxbas_pw_pool, rho_tot_gspace, &
                              use_data=COMPLEXDATA1D, &
                              in_space=RECIPROCALSPACE, error=error)

      CALL pw_pool_init_coeff(auxbas_pw_pool, v_hartree_rspace, &
                              use_data=REALDATA3D, in_space=REALSPACE)

      ! *** Calculation of the forces due to the Hartree energy *** !
      CALL coeff_add(rho_core, harris_env%rho%rho_g(1), rho_tot_gspace)
      
      CALL init_pw_poisson_solver(green,error=error)
      CALL calculate_hartree(rho_tot_gspace, ehartree=Ehartree, vhartree=v_hartree_gspace)
      CALL cleanup_pw_poisson_solver(green,error=error)

      CALL pw_pool_give_back_coeff(auxbas_pw_pool, rho_tot_gspace, error=error)

      CALL coeff_transform_space(v_hartree_gspace, v_hartree_rspace)

      CALL qs_vxc_create(qs_env=qs_env, vxc_rho=v_rspace_new, vxc_tau=v_tau_rspace, &
                         exc=Exc, just_energy=.FALSE., error=error)

      CALL coeff_add(v_rspace_new(1), v_hartree_rspace, v_rspace_new(1))
      CALL coeff_scale(v_rspace_new(1), v_rspace_new(1)%pw%pw_grid%dvol)

      CALL integrate_v_rspace(v_rspace=v_rspace_new(1), p=qs_env%rho%rho_ao(1), &
                              h=matrix_ks(1), qs_env=qs_env, calculate_forces=.TRUE., &
                              gapw=.FALSE., error=error)

      DO ikind = 1,SIZE(force)
        CALL mp_sum(force(ikind)%rho_elec, globenv%group)
      END DO

      DO iatom = 1, natom
        ikind = kind_of(iatom)
        i = atom_of_kind(iatom)
        harris_force%f_hartree(iatom, 1:3) = force(ikind)%rho_elec(1:3,i)
      END DO

      CALL zero_qs_force(force)

      ! ** Force computation of f_V  ** !
      CALL init_pw_poisson_solver(green, error=error)
      CALL calculate_hartree(qs_env%rho%rho_g(1), ehartree=Ehartree, vhartree=v_hartree_gspace)
      CALL cleanup_pw_poisson_solver(green, error=error)

      CALL coeff_transform_space(v_hartree_gspace, v_hartree_rspace)
      CALL coeff_scale(v_hartree_rspace, v_hartree_rspace%pw%pw_grid%dvol)

      CALL integrate_v_core_rspace(v_hartree_rspace, qs_env)

      DO ikind = 1,SIZE(force)
        CALL mp_sum(force(ikind)%rho_core, globenv%group)
      END DO
      
      DO iatom = 1, natom
        ikind = kind_of(iatom)
        i = atom_of_kind(iatom)
        harris_force%f_V(iatom, 1:3) = force(ikind)%rho_core(1:3,i)
      END DO

      CALL zero_qs_force(force)

      ! *** Calculation of the forces of the xc-integral *** !
      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, cell=cell, error=error)

      CALL qs_rho_get(rho_struct=rho, rho_r=rho_r, error=error)
      ALLOCATE(my_rho_r(SIZE(rho_r)), stat=stat)
      CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
      DO ispin=1, SIZE(rho_r)
        my_rho_r(ispin)%pw => harris_env%rho%rho_r(ispin)%pw
      END DO
      CALL xc_prep_2nd_deriv_new(deriv_set=deriv_set, rho_set=rho_set, &
                                 rho_r=my_rho_r, pw_pool=auxbas_pw_pool, &
                                 dft_control=dft_control, cell=cell, error=error)
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
        CALL pw_zero(v_xc(ispin)%pw)
      END DO

      CALL xc_calc_2nd_deriv_new(v_xc, deriv_set=deriv_set, rho_set=rho_set, &
                                rho1_set=rho_set, pw_pool=auxbas_pw_pool, &
                                dft_control=dft_control, gapw=.FALSE., error=error)

      CALL xc_dset_release(deriv_set, error=error)

      v_xc(1)%pw%cr3d = v_xc(1)%pw%cr3d * v_xc(1)%pw%pw_grid%dvol
      v_rspace_new(1)%pw%cr3d = rho_r(1)%pw%cr3d * v_xc(1)%pw%cr3d
      CALL coeff_scale(v_rspace_new(1), v_rspace_new(1)%pw%pw_grid%dvol)

      CALL pw_release(v_xc(1)%pw, error=error)
      DEALLOCATE(v_xc,stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      
      CALL xc_rho_set_release(rho_set)

      ! + rho_out
      CALL init_pw_poisson_solver(green,error=error)
      CALL calculate_hartree(qs_env%rho%rho_g(1), ehartree=Ehartree, &
                             vhartree=v_hartree_gspace)
      CALL cleanup_pw_poisson_solver(green,error=error)

      CALL coeff_transform_space(v_hartree_gspace, v_hartree_rspace)
      CALL coeff_scale(v_hartree_rspace, v_hartree_rspace%pw%pw_grid%dvol)

      CALL coeff_add(v_rspace_new(1), v_hartree_rspace, v_rspace_new(1))
      
      CALL integrate_v_rspace(v_rspace=v_rspace_new(1), p=harris_env%rho%rho_ao(1), &
                              h=matrix_ks(1), qs_env=qs_env, calculate_forces=.TRUE., &
                              gapw=.FALSE., error=error)

      IF (ASSOCIATED(v_rspace_new)) THEN
        CALL pw_pool_give_back_coeff(auxbas_pw_pool,v_rspace_new(1),error=error)
        DEALLOCATE(v_rspace_new,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      END IF
      IF (ASSOCIATED(v_tau_rspace)) THEN
        CALL pw_pool_give_back_coeff(auxbas_pw_pool,v_tau_rspace(1),error=error)
        DEALLOCATE(v_tau_rspace,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      END IF

      DO ikind = 1,SIZE(force)
        CALL mp_sum(force(ikind)%rho_elec, globenv%group)
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
        harris_force%f_hartree(iatom, 1:3) = 0.0_dp
      END DO

      CALL zero_qs_force(force)     
      ! *** End of force calculation due to the sum_of_eigenvalues *** !


      CALL harris_env_get(harris_env=harris_env, rho=rho, error=error)

      ! *** The forces due to hartree core-core repulsion *** !
      CALL init_pw_poisson_solver(green,error=error)
      CALL calculate_hartree(rho_core, ehartree=Ehartree, vhartree=v_hartree_gspace)
      CALL cleanup_pw_poisson_solver(green,error=error)

      CALL coeff_transform_space(v_hartree_gspace, v_hartree_rspace)
      CALL coeff_scale(v_hartree_rspace, v_hartree_rspace%pw%pw_grid%dvol)

      CALL integrate_v_core_rspace(v_hartree_rspace, qs_env)

      DO ikind = 1,SIZE(force)
        CALL mp_sum(force(ikind)%rho_core, globenv%group)
      END DO
      
      DO iatom = 1, natom
        ikind = kind_of(iatom)
        i = atom_of_kind(iatom)
        harris_force%f_rho_core(iatom, 1:3) = force(ikind)%rho_core(1:3,i)
      END DO

      CALL zero_qs_force(force)
      

      ! *** The forces due to exchange and correlation *** !
      CALL qs_vxc_create(qs_env=qs_env, vxc_rho=v_rspace_new, vxc_tau=v_tau_rspace, &
                         exc=Exc, just_energy=.FALSE., error=error)

      ! *** Calculation of the forces of the xc-integral *** !
      CALL qs_rho_get(rho, rho_r=rho_r, error=error)
      ALLOCATE(my_rho_r(SIZE(rho_r)), stat=stat)
      CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
      DO ispin=1, SIZE(rho_r)
        my_rho_r(ispin)%pw => harris_env%rho%rho_r(ispin)%pw
      END DO
      CALL xc_prep_2nd_deriv_new(deriv_set=deriv_set, rho_set=rho_set, &
                                 rho_r=my_rho_r, pw_pool=auxbas_pw_pool, &
                                 dft_control=dft_control, cell=cell, error=error)
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
        CALL pw_zero(v_xc(ispin)%pw)
      END DO

      CALL xc_calc_2nd_deriv_new(v_xc, deriv_set=deriv_set, rho_set=rho_set, &
                                rho1_set=rho_set, pw_pool=auxbas_pw_pool, &
                                dft_control=dft_control, gapw=.FALSE., error=error)

      CALL xc_dset_release(deriv_set, error=error)

      v_xc(1)%pw%cr3d = v_xc(1)%pw%cr3d * v_xc(1)%pw%pw_grid%dvol
      v_rspace_new(1)%pw%cr3d = rho_r(1)%pw%cr3d * v_xc(1)%pw%cr3d
      CALL coeff_scale(v_rspace_new(1), v_rspace_new(1)%pw%pw_grid%dvol)

      CALL pw_release(v_xc(1)%pw, error=error)
      DEALLOCATE(v_xc,stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      
      CALL xc_rho_set_release(rho_set)

      CALL init_pw_poisson_solver(green,error=error)
      CALL calculate_hartree(rho%rho_g(1), ehartree=Ehartree, vhartree=v_hartree_gspace)
      CALL cleanup_pw_poisson_solver(green,error=error)

      CALL coeff_transform_space(v_hartree_gspace, v_hartree_rspace)
      CALL coeff_scale(v_hartree_rspace, v_hartree_rspace%pw%pw_grid%dvol)

      CALL coeff_add(v_rspace_new(1), v_hartree_rspace, v_rspace_new(1))

      CALL integrate_v_rspace(v_rspace=v_rspace_new(1), p=harris_env%rho%rho_ao(1), &
                              h=matrix_ks(1), qs_env=qs_env, calculate_forces=.TRUE., &
                              gapw=.FALSE., error=error)

      IF (ASSOCIATED(v_rspace_new)) THEN
        CALL pw_pool_give_back_coeff(auxbas_pw_pool,v_rspace_new(1),error=error)
        DEALLOCATE(v_rspace_new,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      END IF
      IF (ASSOCIATED(v_tau_rspace)) THEN
        CALL pw_pool_give_back_coeff(auxbas_pw_pool,v_tau_rspace(1),error=error)
        DEALLOCATE(v_tau_rspace,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      END IF

      DO ikind = 1,SIZE(force)
        CALL mp_sum(force(ikind)%rho_elec, globenv%group)
      END DO

      DO iatom = 1, natom
        ikind = kind_of(iatom)
        i = atom_of_kind(iatom)
        harris_force%f_integral_vxc(iatom, 1:3) = force(ikind)%rho_elec(1:3,i)
      END DO

      !CALL zero_qs_force(force)
      !
      !CALL calculate_ecore_overlap(qs_env=qs_env, globenv=globenv, &
      !                             calculate_forces=.TRUE., E_overlap_core=E_overlap_core)
      !
      !DO iatom = 1, natom
      !  ikind = kind_of(iatom)
      !  i = atom_of_kind(iatom)
      !  harris_force%f_ovrl(iatom, 1:3) = force(ikind)%core_overlap(1:3,i)
      !END DO
      
      CALL deallocate_qs_force(force)
      CALL set_qs_env(qs_env=qs_env, force=force_tmp, error=error)
      
      CALL get_qs_env(qs_env=qs_env, force=force, error=error)
      harris_force%f_total(:) = 0.0_dp

      DO ikind = 1,SIZE(force)
        CALL mp_sum(force(ikind)%core_overlap, globenv%group)
      END DO

      DO iatom = 1, natom
        ikind = kind_of(iatom)
        i = atom_of_kind(iatom)

        harris_force%f_ovrl(iatom, 1:3) = force(ikind)%core_overlap(1:3,i)
        harris_force%f_self(iatom, 1:3) = 0.0_dp
        harris_force%f_EII(iatom, 1:3) = harris_force%f_ovrl(iatom, 1:3) &
                                       + harris_force%f_self(iatom, 1:3) &
                                       + harris_force%f_rho_core(iatom, 1:3)
        harris_force%f_harris(iatom, 1:3) = harris_force%f_trace(iatom, 1:3) &
                                          - harris_force%f_hartree(iatom, 1:3) &
                                          - harris_force%f_integral_vxc(iatom, 1:3) &
                                          + harris_force%f_EII(iatom, 1:3)
        harris_force%f_total(1:3) = harris_force%f_total(1:3) &
                                  + harris_force%f_harris(iatom, 1:3)
      END DO

      ! Output
      IF (logger%para_env%source==logger%para_env%mepos) THEN
        unit_nr=cp_logger_get_default_unit_nr(logger,local=.false.)
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

    CALL timestop(0.0_dp, handle)

  END SUBROUTINE harris_force_EVal
!***************************************************************************

!!****f* harris_force/harris_calc_nsc_force [1.0] *
!!
!!   NAME
!!     harris_calc_nsc_force
!!
!!   FUNCTION
!!     Calculates the non-self-consistent force correction, in order to keep
!!     the forces consistent with the according energies.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - error: variable to control error logging, stopping,... 
!!              see module cp_error_handling 
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     11.2005 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE harris_calc_nsc_force(qs_env, globenv, error)
    TYPE(qs_environment_type), POINTER              :: qs_env
    TYPE(global_environment_type), POINTER          :: globenv
    TYPE(cp_error_type), INTENT(INOUT), OPTIONAL    :: error

    INTEGER                                         :: handle, &
                                                       natom, &
                                                       nkind, &
                                                       stat, &
                                                       ikind, &
                                                       iatom, &
                                                       i, &
                                                       ispin, &
                                                       nspins
    LOGICAL                                         :: failure
    INTEGER, ALLOCATABLE, DIMENSION(:)              :: atom_of_kind, &
                                                       natom_of_kind, &
                                                       kind_of
    TYPE(harris_env_type), POINTER                  :: harris_env
    TYPE(harris_force_type),  POINTER               :: harris_force
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                       :: atomic_kind_set
    TYPE(qs_force_type), DIMENSION(:), POINTER      :: force, &
                                                       force_tmp
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                       :: particle_set
    REAL(KIND=dp)                                   :: Ehartree, &
                                                       Exc
    TYPE(coeff_type)                                :: v_hartree_gspace, &
                                                       v_hartree_rspace
    TYPE(greens_function_type), POINTER             :: green
    TYPE(coeff_type), DIMENSION(:), POINTER         :: v_rspace_new, &
                                                       v_tau_rspace, &
                                                       rho_r
    TYPE(pw_p_type), DIMENSION(:), POINTER          :: my_rho_r, &
                                                       v_xc
    TYPE(dft_control_type), POINTER                 :: dft_control
    TYPE(cell_type), POINTER                        :: cell
    TYPE(xc_derivative_set_type), POINTER           :: deriv_set
    TYPE(xc_rho_set_type), POINTER                  :: rho_set
    TYPE(pw_pool_type), POINTER                     :: auxbas_pw_pool
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: matrix_ks
    TYPE(pw_env_type), POINTER                      :: pw_env

    CHARACTER(len=*), PARAMETER :: routineN = 'harris_calc_nsc_force', &
      routineP = moduleN//':'//routineN

!   ------------------------------------------------------------------------

    CALL timeset(routineN, "I", "", handle)

    failure = .FALSE.
    NULLIFY(harris_env, harris_force, force, force_tmp, particle_set, green, &
            v_rspace_new, v_tau_rspace, rho_r, my_rho_r, v_xc, dft_control, &
            cell, deriv_set, rho_set, auxbas_pw_pool, matrix_ks, pw_env)

    CPPrecondition(ASSOCIATED(qs_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(qs_env%ref_count>0, cp_failure_level, routineP, error, failure)

    IF (.NOT. failure) THEN
      CALL get_qs_env(qs_env=qs_env, harris_env=harris_env, green=green, &
                      atomic_kind_set=atomic_kind_set, force=force, pw_env=pw_env, &
                      particle_set=particle_set, matrix_ks=matrix_ks, error=error)
      CALL harris_env_get(harris_env=harris_env, harris_force=harris_force, error=error)

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

      CALL diff_rho_type(rho_input1=qs_env%rho, rho_input2=harris_env%rho, &
                         rho_output=harris_env%rho_diff, qs_env=qs_env, error=error)

      CALL duplicate_qs_force(qs_force_input=force, qs_force_output=force_tmp, &
                              natom_of_kind=natom_of_kind)
      CALL zero_qs_force(force)

      CALL pw_env_get(pw_env=pw_env, auxbas_pw_pool=auxbas_pw_pool, error=error)
      CALL pw_pool_init_coeff(auxbas_pw_pool, v_hartree_gspace, &
                              use_data = COMPLEXDATA1D, &
                              in_space = RECIPROCALSPACE, error=error)
      CALL pw_pool_init_coeff(auxbas_pw_pool, v_hartree_rspace, &
                              use_data=REALDATA3D, in_space=REALSPACE)

      ! *** Calculation of the forces due to the Hartree energy *** !
      CALL init_pw_poisson_solver(green,error=error)
      CALL calculate_hartree(harris_env%rho_diff%rho_g(1), ehartree=Ehartree, &
                             vhartree=v_hartree_gspace)
      CALL cleanup_pw_poisson_solver(green,error=error)

      CALL coeff_transform_space(v_hartree_gspace, v_hartree_rspace)
      CALL coeff_scale(v_hartree_rspace, v_hartree_rspace%pw%pw_grid%dvol)

      CALL integrate_v_rspace(v_rspace=v_hartree_rspace, p=harris_env%rho_diff%rho_ao(1), &
                              h=matrix_ks(1), qs_env=qs_env, calculate_forces=.TRUE., &
                              gapw=.FALSE., error=error)

      DO ikind = 1,SIZE(force)
        CALL mp_sum(force(ikind)%rho_elec, globenv%group)
      END DO

      DO iatom = 1, natom
        ikind = kind_of(iatom)
        i = atom_of_kind(iatom)
        harris_force%f_nsc(iatom, 1:3) = force(ikind)%rho_elec(1:3,i)
      END DO

      CALL zero_qs_force(force)

      CALL init_pw_poisson_solver(green,error=error)
      CALL calculate_hartree(harris_env%rho_diff%rho_g(1), ehartree=Ehartree, &
                             vhartree=v_hartree_gspace)
      CALL cleanup_pw_poisson_solver(green,error=error)

      CALL coeff_transform_space(v_hartree_gspace, v_hartree_rspace)
      CALL coeff_scale(v_hartree_rspace, v_hartree_rspace%pw%pw_grid%dvol)

      CALL integrate_v_rspace(v_rspace=v_hartree_rspace, p=harris_env%rho_diff%rho_ao(1), &
                              h=matrix_ks(1), qs_env=qs_env, calculate_forces=.TRUE., &
                              gapw=.FALSE.,error=error)

      DO ikind = 1,SIZE(force)
        CALL mp_sum(force(ikind)%rho_elec, globenv%group)
      END DO

      DO iatom = 1, natom
        ikind = kind_of(iatom)
        i = atom_of_kind(iatom)
        harris_force%f_nsc(iatom, 1:3) = harris_force%f_nsc(iatom, 1:3) &
                                       + force(ikind)%rho_elec(1:3,i)
      END DO

      CALL zero_qs_force(force)

      ! *** Calculation of the forces of the xc-integral *** !
      CALL qs_vxc_create(qs_env=qs_env, vxc_rho=v_rspace_new, vxc_tau=v_tau_rspace, &
                         exc=Exc, just_energy=.FALSE., harris=.TRUE., nsc_force=.TRUE., &
                         error=error)
      !CALL qs_vxc_create(qs_env=qs_env, vxc_rho=v_rspace_new, vxc_tau=v_tau_rspace, &
      !                   exc=Exc, just_energy=.FALSE., harris=.TRUE., nsc_force=.FALSE., &
      !                   error=error)
      !CALL qs_vxc_create(qs_env=qs_env, vxc_rho=v_rspace_new, vxc_tau=v_tau_rspace, &
      !                   exc=Exc, just_energy=.FALSE., harris=.FALSE., error=error)
      CALL coeff_scale(v_rspace_new(1), v_rspace_new(1)%pw%pw_grid%dvol)
      CALL integrate_v_rspace(v_rspace=v_rspace_new(1), p=harris_env%rho_diff%rho_ao(1), &
                              h=matrix_ks(1), qs_env=qs_env, calculate_forces=.TRUE., &
                              gapw=.FALSE., error=error)

      DO ikind = 1,SIZE(force)
        CALL mp_sum(force(ikind)%rho_elec, globenv%group)
      END DO

      DO iatom = 1, natom
        ikind = kind_of(iatom)
        i = atom_of_kind(iatom)
        !harris_force%f_nsc(iatom, 1:3) = force(ikind)%rho_elec(1:3,i)
      END DO

      CALL zero_qs_force(force)

      CALL get_qs_env(qs_env=qs_env, dft_control=dft_control, cell=cell, error=error)

      CALL qs_rho_get(rho_struct=harris_env%rho_diff, rho_r=rho_r, error=error)
      !CALL qs_rho_get(rho_struct=harris_env%rho_diff, rho_r=rho_r, error=error)
      ALLOCATE(my_rho_r(SIZE(rho_r)), stat=stat)
      CPPostcondition(stat==0, cp_failure_level, routineP, error, failure)
      DO ispin=1, SIZE(rho_r)
        my_rho_r(ispin)%pw => rho_r(ispin)%pw
      END DO
      CALL xc_prep_2nd_deriv_new(deriv_set=deriv_set, rho_set=rho_set, &
                                 rho_r=my_rho_r, pw_pool=auxbas_pw_pool, &
                                 dft_control=dft_control, cell=cell, error=error)
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
        CALL pw_zero(v_xc(ispin)%pw)
      END DO

      CALL xc_calc_2nd_deriv_new(v_xc, deriv_set=deriv_set, rho_set=rho_set, &
                                rho1_set=rho_set, pw_pool=auxbas_pw_pool, &
                                dft_control=dft_control, gapw=.FALSE., error=error)

      CALL xc_dset_release(deriv_set, error=error)

      v_xc(1)%pw%cr3d = v_xc(1)%pw%cr3d * v_xc(1)%pw%pw_grid%dvol
      v_rspace_new(1)%pw%cr3d = rho_r(1)%pw%cr3d * v_xc(1)%pw%cr3d
      CALL coeff_scale(v_rspace_new(1), v_rspace_new(1)%pw%pw_grid%dvol)

      CALL pw_release(v_xc(1)%pw, error=error)
      DEALLOCATE(v_xc,stat=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      
      CALL xc_rho_set_release(rho_set)

      CALL integrate_v_rspace(v_rspace=v_rspace_new(1), p=harris_env%rho_diff%rho_ao(1), &
                              h=matrix_ks(1), qs_env=qs_env, calculate_forces=.TRUE., &
                              gapw=.FALSE., error=error)


      IF (ASSOCIATED(v_rspace_new)) THEN
        CALL pw_pool_give_back_coeff(auxbas_pw_pool,v_rspace_new(1),error=error)
        DEALLOCATE(v_rspace_new,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      END IF
      IF (ASSOCIATED(v_tau_rspace)) THEN
        CALL pw_pool_give_back_coeff(auxbas_pw_pool,v_tau_rspace(1),error=error)
        DEALLOCATE(v_tau_rspace,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      END IF

      DO ikind = 1,SIZE(force)
        CALL mp_sum(force(ikind)%rho_elec, globenv%group)
      END DO

      DO iatom = 1, natom
        ikind = kind_of(iatom)
        i = atom_of_kind(iatom)
        !harris_force%f_nsc(iatom, 1:3) = harris_force%f_nsc(iatom,1:3) + force(ikind)%rho_elec(1:3,i)
      END DO

    END IF

    CALL timestop(0.0_dp, handle)

  END SUBROUTINE harris_calc_nsc_force

!!****f* harris_force/harris_force_test_rho_core [1.0] *
!!
!!   NAME
!!     harris_force_test_rho_core
!!
!!   FUNCTION
!!     Performs a test if the forces in the harris functional due to 
!!     rho_core sums up to force%rho_core.
!!
!!   NOTES
!!     Only for internal debugging purposes
!!
!!   INPUTS
!!     - qs_env: The QS environment of matter
!!     - harris_force: The harris force type in which all necessary informations 
!!                     are stored
!!     - error: variable to control error logging, stopping,... 
!!              see module cp_error_handling 
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     06.2005 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE harris_force_test_rho_core(qs_env, harris_force, error)
    TYPE(qs_environment_type), POINTER              :: qs_env
    TYPE(harris_force_type),  POINTER               :: harris_force
    TYPE(cp_error_type), INTENT(INOUT), OPTIONAL    :: error
    INTEGER                                         :: natom, &
                                                       nkind, &
                                                       stat, &
                                                       iatom, &
                                                       ikind, &
                                                       i
    INTEGER, ALLOCATABLE, DIMENSION(:)              :: atom_of_kind, &
                                                       natom_of_kind, &
                                                       kind_of
    TYPE(particle_type), DIMENSION(:), POINTER      :: particle_set
    TYPE(atomic_kind_type), DIMENSION(:), POINTER   :: atomic_kind_set
    TYPE(qs_force_type), DIMENSION(:), POINTER      :: force
    LOGICAL                                         :: failure

    CHARACTER(len=*), PARAMETER :: routineN = 'harris_force_test_rho_core', &
      routineP = moduleN//':'//routineN
    TYPE(cp_logger_type), POINTER :: logger
    INTEGER :: unit_nr

!   ------------------------------------------------------------------------

    failure = .FALSE.
    logger => cp_error_get_logger(error)

    CPPrecondition(ASSOCIATED(qs_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(qs_env%ref_count>0, cp_failure_level, routineP, error, failure)
    CPPrecondition(ASSOCIATED(harris_force), cp_failure_level, routineP, error, failure)
    CPPrecondition(harris_force%ref_count>0, cp_failure_level, routineP, error, failure)
    
    IF (logger%para_env%source==logger%para_env%mepos) THEN
          unit_nr=cp_logger_get_default_unit_nr(logger,local=.false.)
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

      IF (logger%para_env%source==logger%para_env%mepos) THEN
        unit_nr=cp_logger_get_default_unit_nr(logger,local=.false.)
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
!***************************************************************************

!!****f* harris_force/harris_force_test_rho_elec [1.0] *
!!
!!   NAME
!!     harris_force_test_rho_core
!!
!!   FUNCTION
!!     Performs a test if the forces in the harris functional due to 
!!     rho_elec sums up to force%rho_elec.
!!
!!   NOTES
!!     Only for internal debugging purposes
!!
!!   INPUTS
!!     - qs_env: The QS environment of matter
!!     - harris_force: The harris force type in which all necessary informations 
!!                     are stored
!!     - error: variable to control error logging, stopping,... 
!!              see module cp_error_handling 
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     06.2005 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE harris_force_test_rho_elec(qs_env, harris_force, error)
    TYPE(qs_environment_type), POINTER              :: qs_env
    TYPE(harris_force_type),  POINTER               :: harris_force
    TYPE(cp_error_type), INTENT(INOUT), OPTIONAL    :: error
    INTEGER                                         :: natom, &
                                                       nkind, &
                                                       stat, &
                                                       iatom, &
                                                       ikind, &
                                                       i
    INTEGER, ALLOCATABLE, DIMENSION(:)              :: atom_of_kind, &
                                                       natom_of_kind, &
                                                       kind_of
    TYPE(particle_type), DIMENSION(:), POINTER      :: particle_set
    TYPE(atomic_kind_type), DIMENSION(:), POINTER   :: atomic_kind_set
    TYPE(qs_force_type), DIMENSION(:), POINTER      :: force
    LOGICAL                                         :: failure

    CHARACTER(len=*), PARAMETER :: routineN = 'harris_force_test_rho_elec', &
      routineP = moduleN//':'//routineN
    TYPE(cp_logger_type), POINTER :: logger
    INTEGER :: unit_nr
    
!   ------------------------------------------------------------------------

    failure = .FALSE.
    logger => cp_error_get_logger(error)

    CPPrecondition(ASSOCIATED(qs_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(qs_env%ref_count>0, cp_failure_level, routineP, error, failure)
    CPPrecondition(ASSOCIATED(harris_force), cp_failure_level, routineP, error, failure)
    CPPrecondition(harris_force%ref_count>0, cp_failure_level, routineP, error, failure)
    
    IF (logger%para_env%source==logger%para_env%mepos) THEN
          unit_nr=cp_logger_get_default_unit_nr(logger,local=.false.)
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

      IF (logger%para_env%source==logger%para_env%mepos) THEN
        unit_nr=cp_logger_get_default_unit_nr(logger,local=.false.)
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
!***************************************************************************

!!****f* harris_force/harris_force_test_integral_vxc [1.0] *
!!
!!   NAME
!!     harris_force_test_rho_core
!!
!!   FUNCTION
!!     
!!   NOTES
!!     Only for internal debugging purposes
!!
!!   INPUTS
!!     - qs_env: The QS environment of matter
!!     - harris_force: The harris force type in which all necessary informations 
!!                     are stored
!!     - error: variable to control error logging, stopping,... 
!!              see module cp_error_handling 
!!
!!   AUTHOR
!!     Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
!!
!!   MODIFICATION HISTORY
!!     06.2005 initial create [tdk]
!!
!!*** **********************************************************************
  SUBROUTINE harris_force_test_integral_vxc(qs_env, harris_force, error)
    TYPE(qs_environment_type), POINTER              :: qs_env
    TYPE(harris_force_type),  POINTER               :: harris_force
    TYPE(cp_error_type), INTENT(INOUT), OPTIONAL    :: error
    INTEGER                                         :: natom, &
                                                       nkind, &
                                                       stat, &
                                                       iatom, &
                                                       ikind, &
                                                       i
    INTEGER, ALLOCATABLE, DIMENSION(:)              :: atom_of_kind, &
                                                       natom_of_kind, &
                                                       kind_of
    TYPE(particle_type), DIMENSION(:), POINTER      :: particle_set
    TYPE(atomic_kind_type), DIMENSION(:), POINTER   :: atomic_kind_set
    TYPE(qs_force_type), DIMENSION(:), POINTER      :: force, &
                                                       force_tmp
    TYPE(qs_rho_type), POINTER                      :: rho
    TYPE(coeff_type), DIMENSION(:), POINTER         :: v_rspace_new, &
                                                       v_tau_rspace
    TYPE(coeff_type)                                :: v_hartree_rspace, &
                                                       v_hartree_gspace, &
                                                       rho_tot_gspace
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: matrix_ks
    TYPE(greens_function_type), POINTER             :: green
    REAL(KIND=dp)                                   :: Ehartree, &
                                                       Exc
    TYPE(pw_pool_type), POINTER                     :: auxbas_pw_pool
    TYPE(pw_env_type), POINTER                      :: pw_env
    TYPE(coeff_type), POINTER                       :: rho_core
    LOGICAL                                         :: failure

    CHARACTER(len=*), PARAMETER :: routineN = 'harris_force_test_integral_vxc', &
      routineP = moduleN//':'//routineN
    TYPE(cp_logger_type), POINTER :: logger
    INTEGER :: unit_nr
    
!   ------------------------------------------------------------------------

    failure = .FALSE.
    NULLIFY(rho, v_rspace_new, matrix_ks, green, auxbas_pw_pool, pw_env, &
            rho_core, force_tmp, v_tau_rspace)
    logger => cp_error_get_logger(error)

    CPPrecondition(ASSOCIATED(qs_env), cp_failure_level, routineP, error, failure)
    CPPrecondition(qs_env%ref_count>0, cp_failure_level, routineP, error, failure)
    CPPrecondition(ASSOCIATED(harris_force), cp_failure_level, routineP, error, failure)
    CPPrecondition(harris_force%ref_count>0, cp_failure_level, routineP, error, failure)
    
    IF (logger%para_env%source==logger%para_env%mepos) THEN
          unit_nr=cp_logger_get_default_unit_nr(logger,local=.false.)
          WRITE (unit_nr,*) ""; WRITE (unit_nr, *) ""
          WRITE (unit_nr,*) "The Harris force correction integral_vxc test is performed!"
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
      
      ! *** Debug Code: Calculate the forces of the sum of eigenvalues *** !
      CALL duplicate_qs_force(qs_force_input=force, qs_force_output=force_tmp, &
                              natom_of_kind=natom_of_kind)
      CALL zero_qs_force(force)

      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env, rho=rho, rho_core=rho_core, &
                      matrix_ks=matrix_ks, green=green, error=error)
      CALL pw_env_get(pw_env=pw_env, auxbas_pw_pool=auxbas_pw_pool, error=error)
      CALL pw_pool_init_coeff(auxbas_pw_pool, v_hartree_gspace, &
                              use_data = COMPLEXDATA1D,&
                              in_space = RECIPROCALSPACE, error=error)
      CALL pw_pool_init_coeff(auxbas_pw_pool, rho_tot_gspace, &
                              use_data=COMPLEXDATA1D, &
                              in_space=RECIPROCALSPACE, error=error)
            
      CALL coeff_add(rho_core, rho%rho_g(1), rho_tot_gspace)
      
      CALL init_pw_poisson_solver(green,error=error)
      CALL calculate_hartree(rho_tot_gspace, ehartree=Ehartree, vhartree=v_hartree_gspace)
      !CALL calculate_hartree(rho%rho_g(1), ehartree=Ehartree, vhartree=v_hartree_gspace)
      CALL cleanup_pw_poisson_solver(green,error=error)

      CALL pw_pool_give_back_coeff(auxbas_pw_pool, rho_tot_gspace, &
                                   error=error)
      
      CALL pw_pool_init_coeff(auxbas_pw_pool, v_hartree_rspace, &
                              use_data=REALDATA3D, in_space=REALSPACE)
      CALL coeff_transform_space(v_hartree_gspace, v_hartree_rspace)
      !CALL coeff_scale(v_hartree_rspace, v_hartree_rspace%pw%pw_grid%dvol)
      
      CALL pw_pool_give_back_coeff(auxbas_pw_pool, v_hartree_gspace, &
                                   error=error)
      
      NULLIFY(v_rspace_new)
      CALL qs_vxc_create(qs_env=qs_env, vxc_rho = v_rspace_new, vxc_tau=v_tau_rspace, &
                         exc=Exc, just_energy = .FALSE., error=error)
      
      !CALL coeff_scale(v_rspace_new(1), v_rspace_new(1)%pw%pw_grid%dvol)
      
      ! sum-up: v_hartree_rspace = v_hartree_rspace + v_rspace_new(1)
      CALL coeff_sumup(v_rspace_new(1), v_hartree_rspace)
      CALL coeff_scale(v_hartree_rspace, v_hartree_rspace%pw%pw_grid%dvol)

      IF (ASSOCIATED(v_rspace_new)) THEN
        CALL pw_pool_give_back_coeff(auxbas_pw_pool,v_rspace_new(1),error=error)
      END IF
      IF (ASSOCIATED(v_tau_rspace)) THEN
        CALL pw_pool_give_back_coeff(auxbas_pw_pool,v_tau_rspace(1),error=error)
      END IF
      
      CALL integrate_v_rspace(v_rspace=v_hartree_rspace, p=rho%rho_ao(1), &
                              h=matrix_ks(1), qs_env=qs_env, calculate_forces=.TRUE., &
                              gapw=.FALSE.,error=error)
      CALL pw_pool_give_back_coeff(auxbas_pw_pool, v_hartree_rspace, error=error)
      
      DO iatom = 1, natom
        ikind = kind_of(iatom)
        i = atom_of_kind(iatom)
        WRITE (unit_nr,*) force(ikind)%rho_elec(1:3,i)
      END DO
      
      CALL deallocate_qs_force(force)
      CALL set_qs_env(qs_env=qs_env, force=force_tmp)
      CALL get_qs_env(qs_env=qs_env, force=force, error=error)

      IF (logger%para_env%source==logger%para_env%mepos) THEN
        unit_nr=cp_logger_get_default_unit_nr(logger,local=.false.)
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
               + force(ikind)%gth_ppnl(1:3,i)
        END DO
  
      END IF

      IF (ASSOCIATED(v_rspace_new)) THEN
        DEALLOCATE(v_rspace_new,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
      END IF
      
      IF (ASSOCIATED(v_tau_rspace)) THEN
        DEALLOCATE(v_tau_rspace,stat=stat)
        CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
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
