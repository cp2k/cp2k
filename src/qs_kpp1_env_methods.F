!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/qs_kpp1_env_methods [1.0] *
!!
!!   NAME
!!     qs_kpp1_env_methods
!!
!!   FUNCTION
!!     module that builds the second order perturbation kernel
!!     kpp1 = delta_rho|_P delta_rho|_P E drho(P1) drho
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     07.2002 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE qs_kpp1_env_methods
  USE coefficient_types,               ONLY: coeff_add,&
                                             coeff_allocate,&
                                             coeff_copy,&
                                             coeff_deallocate,&
                                             coeff_sumup,&
                                             coeff_transform_space,&
                                             coeff_type,&
                                             coeff_zero
  USE cp_block_matrix,                 ONLY: cp_block_matrix_p_type,&
                                             cp_block_matrix_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_error_get_logger,&
                                             cp_error_type,&
                                             cp_unimplemented_error,&
                                             cp_unimplemented_error_nr
  USE cp_fm_types,                     ONLY: cp_full_matrix_p_type,&
                                             cp_full_matrix_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_log,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_output_handling,              ONLY: cp_finished_output,&
                                             cp_unitnr_for_output,&
                                             cp_would_output
  USE cube_utils,                      ONLY: cube_info_type
  USE dft_types,                       ONLY: dft_control_type
  USE gaussian_gridlevels,             ONLY: gridlevel_info_type
  USE global_types,                    ONLY: DEBUG,&
                                             global_environment_type
  USE greens_fn,                       ONLY: pw_green_fn_init
  USE hartree,                         ONLY: calculate_hartree
  USE kinds,                           ONLY: wp=>dp
  USE l_utils,                         ONLY: l_info_type
  USE message_passing,                 ONLY: mp_allgather,&
                                             mp_max,&
                                             mp_range,&
                                             mp_sum,&
                                             mp_sync
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_grid_types,                   ONLY: HALFSPACE,&
                                             pw_grid_type
  USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                             pw_pool_give_back_coeff,&
                                             pw_pool_give_back_pw,&
                                             pw_pool_init_coeff,&
                                             pw_pool_p_type,&
                                             pw_pool_type,&
                                             pw_pools_give_back_coeffs,&
                                             pw_pools_init_coeffs
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             COMPLEXDATA3D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_copy,&
                                             pw_derive,&
                                             pw_p_type,&
                                             pw_sumup,&
                                             pw_transfer,&
                                             pw_type,&
                                             pw_zero
  USE pws,                             ONLY: cleanup_pw_poisson_solver,&
                                             init_pw_poisson_solver,&
                                             pw_poisson_solver
  USE qs_charges_types,                ONLY: qs_charges_create,&
                                             qs_charges_release,&
                                             qs_charges_type
  USE qs_collocate_density,            ONLY: calculate_rho_core,&
                                             calculate_rho_elec,&
                                             calculate_total_rho
  USE qs_energy_types,                 ONLY: allocate_qs_energy,&
                                             deallocate_qs_energy,&
                                             qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_integrate_potential,          ONLY: integrate_v_core_rspace,&
                                             integrate_v_rspace
  USE qs_overlap,                      ONLY: write_blacs_matrix,&
                                             write_sparse_matrix
  USE qs_kpp1_env_types,                      ONLY: kpp1_release,&
                                             kpp1_retain,&
                                             qs_kpp1_env_p_type,&
                                             qs_kpp1_env_type
  USE qs_rho_methods,                  ONLY: qs_rho_update_rho
  USE qs_rho_types,                    ONLY: qs_rho_did_change,&
                                             qs_rho_get,&
                                             qs_rho_type
  USE simulation_cell,                 ONLY: cell_type,&
                                             get_cell
  USE sparse_matrix_types,             ONLY: add_matrices,&
                                             allocate_matrix,&
                                             allocate_matrix_set,&
                                             deallocate_matrix,&
                                             deallocate_matrix_set,&
                                             first_block_node,&
                                             get_block_node,&
                                             get_matrix_info,&
                                             next_block_node,&
                                             real_block_node_type,&
                                             real_matrix_p_type,&
                                             real_matrix_type,&
                                             replicate_matrix_structure,&
                                             set_matrix,&
                                             symmetrise_diagonal_blocks,&
                                             transfer_matrix
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc_functionals,                  ONLY: xc_create_2nd_deriv_info,&
                                             xc_get_pot_size,&
                                             xc_lda,&
                                             xc_set
  IMPLICIT NONE

  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qs_kpp1_env_methods'
  INTEGER, PRIVATE, SAVE :: last_kpp1_id_nr=0

  PUBLIC :: kpp1_create, &
       kpp1_calc_k_p_p1, kpp1_did_change
!!***
!****************************************************************************

CONTAINS

!!****f* qs_kpp1_env_methods/kpp1_create [1.0] *
!!
!!   NAME
!!     kpp1_create
!!
!!   SYNOPSIS
!!     Subroutine kpp1_create(kpp1_env, qs_env, rebuild_each, error)
!!       Type(qs_kpp1_env_type), Pointer:: kpp1_env
!!       Type(qs_environment_type), Intent (IN), Target:: qs_env
!!       Integer, Intent (IN), Optional:: rebuild_each
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine kpp1_create
!!
!!   FUNCTION
!!     allocates and initializes a kpp1_env
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - kpp1_env: the environement to initialize
!!     - qs_env: the qs environment in which this environment lives
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     07.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE kpp1_create(kpp1_env,qs_env, rebuild_each, error)
    TYPE(qs_kpp1_env_type), POINTER :: kpp1_env
    TYPE(qs_environment_type), INTENT(in), TARGET :: qs_env
    INTEGER, INTENT(in), OPTIONAL :: rebuild_each
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    INTEGER :: stat
    CHARACTER(len=*), PARAMETER :: routineN='kpp1_create',&
         routineP=moduleN//':'//routineN

    failure=.FALSE.

    ALLOCATE(kpp1_env,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       NULLIFY(kpp1_env%v_rspace, kpp1_env%v_ao, kpp1_env%drho_r,&
            kpp1_env%deriv_xc, kpp1_env%spin_pot, kpp1_env%grad_pot,&
            kpp1_env%ndiag_term)
       kpp1_env%ref_count=1
       kpp1_env%rebuild_each=5
       IF (PRESENT(rebuild_each)) kpp1_env%rebuild_each=rebuild_each
       kpp1_env%updates_since_rebuild=-1
       last_kpp1_id_nr=last_kpp1_id_nr+1
       kpp1_env%id_nr=last_kpp1_id_nr
       kpp1_env%iter=0
    END IF
  END SUBROUTINE kpp1_create
!***************************************************************************

!!****f* qs_p_utils/p_calc_k_p_p1 [1.0] *
!!
!!   NAME
!!     p_calc_k_p_p1
!!
!!   FUNCTION
!!     calcualtes the k_p_p1 kernel of the pertubation theory
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - kpp1_env: the kpp1 environement in which perform the calculation
!!     - qs_env: kpp1's qs_env
!!     - k_p_p1: the sparse matrix that will contain the kernel k_p_p1
!!     - rho: the density where to evaluate the derivatives (i.e. p along
!!       with with its grid representations, that must be valid)
!!     - rho1: the density that represent the first direction along which
!!       you should evaluate the derivatives
!!     - rebuild: if the potential should be reevaluated from scratch
!!       looks at rebuild_each
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     10.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE kpp1_calc_k_p_p1(kpp1_env,qs_env,k_p_p1,rho,rho1,&
       rebuild,error)
    TYPE(qs_kpp1_env_type), POINTER                  :: kpp1_env
    TYPE(qs_environment_type), INTENT(in), TARGET      :: qs_env
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER    :: k_p_p1
    TYPE(qs_rho_type), POINTER                         :: rho, rho1
    LOGICAL, OPTIONAL                                  :: rebuild
    TYPE(cp_error_type), OPTIONAL, INTENT(inout)       :: error

    LOGICAL :: failure, my_rebuild
    INTEGER :: handle, nspins, ispin, stat, output_unit=6
    CHARACTER(len=*), PARAMETER :: routineN='kpp1_calc_k_p_p1',&
         routineP=moduleN//':'//routineN
    REAL(kind=wp) :: total_rho_gspace, energy_hartree
    TYPE(dft_control_type), POINTER               :: dft_control
    TYPE(pw_pool_p_type), DIMENSION(:), POINTER   :: pw_pools
    TYPE(pw_pool_type), POINTER                   :: auxbas_pw_pool
    TYPE(pw_env_type), POINTER                    :: pw_env
    TYPE(coeff_type), DIMENSION(:), POINTER       :: v_rspace_new,&
         rho_r,rho_g,rho1_r, rho1_g
    TYPE(coeff_type), DIMENSION(:), POINTER       :: mgrid_gspace, &
         mgrid_rspace
    TYPE(coeff_type)                              :: v_hartree_gspace,&
         v_hartree_rspace,&
         rho1_tot_gspace
    TYPE(cp_logger_type), POINTER                 :: logger
    TYPE(cell_type), POINTER :: cell

    CALL timeset(routineN,"I","",handle)
    failure=.FALSE.
    NULLIFY(dft_control,auxbas_pw_pool, pw_pools,&
         pw_env,v_rspace_new,rho_r,rho_g,rho1_r, rho1_g,&
         mgrid_rspace, mgrid_gspace,cell)
    logger => cp_error_get_logger(error)

    CPPrecondition(ASSOCIATED(kpp1_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(k_p_p1),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(rho),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(rho1),cp_failure_level,routineP,error,failure)

    IF (.NOT.failure) THEN
       CPPrecondition(kpp1_env%ref_count>0,cp_failure_level,routineP,error,failure)
       CALL kpp1_check_i_alloc(kpp1_env,qs_env=qs_env,error=error)

       CALL get_qs_env(qs_env=qs_env,&
            dft_control=dft_control,&
            pw_env=pw_env,&
            cell=cell,&
            error=error)

       nspins=dft_control%nspins
       CPPrecondition(SIZE(k_p_p1)>=nspins,cp_failure_level,routineP,error,failure)
       CPPrecondition(ASSOCIATED(pw_env),cp_failure_level,routineP,error,failure)
    END IF

    IF (.NOT.failure) THEN
       kpp1_env%iter=kpp1_env%iter+1

! decides if to rebuild or just update
       my_rebuild=kpp1_env%updates_since_rebuild > kpp1_env%rebuild_each
       IF (PRESENT(rebuild)) my_rebuild=rebuild
       IF (kpp1_env%updates_since_rebuild==-1) THEN
          my_rebuild=.TRUE.
       END IF
       IF (my_rebuild) THEN
          kpp1_env%updates_since_rebuild=0
       ELSE
          kpp1_env%updates_since_rebuild=kpp1_env%updates_since_rebuild+1
       END IF
    END IF

! gets the tmp grids
    IF (.NOT. failure) THEN
       CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,&
            pw_pools=pw_pools)
       ALLOCATE(v_rspace_new(nspins), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       ALLOCATE(mgrid_rspace(SIZE(pw_pools)), mgrid_gspace(SIZE(pw_pools)),&
            stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF (.NOT.failure) THEN
       CALL pw_pool_init_coeff(auxbas_pw_pool,v_hartree_gspace,&
            use_data = COMPLEXDATA1D,&
            in_space = RECIPROCALSPACE, error=error)
       CALL pw_pool_init_coeff(auxbas_pw_pool,v_hartree_rspace,&
            use_data = REALDATA3D,&
            in_space = REALSPACE, error=error)
       CALL pw_pools_init_coeffs(pw_pools,mgrid_gspace,&
            use_data = COMPLEXDATA1D,&
            in_space = RECIPROCALSPACE, error=error)
       CALL pw_pools_init_coeffs(pw_pools,mgrid_rspace,&
            use_data = REALDATA3D,&
            in_space = REALSPACE, error=error)
    END IF

! *** calculate the hartree potential on the total density ***
    IF (.NOT. failure) THEN
       CALL pw_pool_init_coeff(auxbas_pw_pool,rho1_tot_gspace,&
            use_data = COMPLEXDATA1D,&
            in_space = RECIPROCALSPACE, error=error)

       CALL coeff_copy(rho1%rho_g(1),rho1_tot_gspace)
       DO ispin=2,nspins
          CALL coeff_sumup(rho1%rho_g(ispin),rho1_tot_gspace)
       END DO

       IF (logger%print_keys%total_densities) THEN
          total_rho_gspace = calculate_total_rho(rho1_tot_gspace)
!FM          if (cp_would_output(logger,outputName="kpp1_charges",&
!FM               fromWhere=routineP,iter=iter,&
!FM               local=.false.,error=error)) THEN
          IF (logger%para_env%mepos==logger%para_env%source) THEN
             output_unit=cp_unitnr_for_output(logger,&
                  outputName="kpp1_charges",&
                  fromWhere=routineP, &
                  iter=200*(kpp1_env%id_nr-1)+kpp1_env%iter,&
                  local=.FALSE.,error=error)
             WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
                  "KPP1 total charge density (r-space):",&
                  SUM(rho1%tot_rho_r),&
                  "KPP1 total charge density (g-space):",&
                  total_rho_gspace
             CALL cp_finished_output(logger, outputName="kpp1_charges",&
                  fromWhere=routineP, &
                  iter=200*(kpp1_env%id_nr-1)+kpp1_env%iter, &
                  unit_nr=output_unit,local=.FALSE.,error=error)
          END IF
       END IF

       CALL init_pw_poisson_solver(cell%green,error=error)
       CALL calculate_hartree(rho1_tot_gspace,energy_hartree,&
            v_hartree_gspace)
       CALL cleanup_pw_poisson_solver(cell%green,error=error)

       CALL pw_pool_give_back_coeff(auxbas_pw_pool,rho1_tot_gspace,&
            error=error)

! *** calculate the xc potential ***       

       DO ispin=1,nspins
          CALL pw_pool_init_coeff(auxbas_pw_pool,v_rspace_new(ispin),&
               use_data = REALDATA3D,&
               in_space = REALSPACE, error=error)
          CALL coeff_zero(v_rspace_new(ispin))
       END DO

       IF (cp_would_output(logger,outputName="kpp1_energies",&
            fromWhere=routineP,iter=200*(kpp1_env%id_nr-1)+kpp1_env%iter,&
            local=.FALSE.,error=error)) THEN
          IF (logger%para_env%mepos==logger%para_env%source) THEN
             output_unit=cp_unitnr_for_output(logger,&
                  outputName="kpp1_energies",&
                  fromWhere=routineP, &
                  iter=200*(kpp1_env%id_nr-1)+kpp1_env%iter,&
                  local=.FALSE.,error=error)
             WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
                  "KPP1 Hartree energy:                           ",&
                  energy_hartree
             CALL cp_finished_output(logger, outputName="kpp1_energies",&
                  fromWhere=routineP, &
                  iter=200*(kpp1_env%id_nr-1)+kpp1_env%iter, &
                  unit_nr=output_unit,local=.FALSE.,error=error)
          END IF
       END IF

       CALL qs_rho_get(rho,rho_r=rho_r,rho_g=rho_g,error=error)
       CALL qs_rho_get(rho1,rho_r=rho1_r,rho_g=rho1_g,error=error)
       CALL xc_calc_2nd_deriv(rho_r=rho_r,drho_r=kpp1_env%drho_r,&
            rho1_r=rho1_r,&
            rho1_g=rho1_g,&
            dft_control=dft_control,&
            pot=kpp1_env%deriv_xc,&
            spin_pot=kpp1_env%spin_pot,grad_pot=kpp1_env%grad_pot,&
            ndiag_term=kpp1_env%ndiag_term,&
            v_rspace_new=v_rspace_new,&
            pw_pool=auxbas_pw_pool,&
            error=error)

       CALL coeff_transform_space(v_hartree_gspace,&
            v_hartree_rspace)

       v_hartree_rspace%pw%cr3d(:,:,:) =&
            v_hartree_rspace%pw%pw_grid%dvol*&
            v_hartree_rspace%pw%cr3d(:,:,:)

       DO ispin=1,nspins
          v_rspace_new(ispin)%pw%cr3d(:,:,:) =&
               v_rspace_new(ispin)%pw%pw_grid%dvol*&
               v_rspace_new(ispin)%pw%cr3d(:,:,:) +&
               v_hartree_rspace%pw%cr3d(:,:,:)

          IF (my_rebuild) CALL set_matrix(kpp1_env%v_ao(ispin)%matrix,0.0_wp)
          CALL integrate_v_rspace(v_rspace=v_rspace_new(ispin),&
               mgrid_rspace=mgrid_rspace,&
               mgrid_gspace=mgrid_gspace,&
               auxbas_grid=pw_env%auxbas_grid,&
               gridlevel_info=pw_env%gridlevel_info,&
               v_rspace_old=kpp1_env%v_rspace(ispin),&
               cube_info=pw_env%cube_info,&
               l_info=pw_env%l_info,&
               rebuild=my_rebuild,p=rho%rho_ao(ispin),&
               h=kpp1_env%v_ao(ispin),&
               qs_env=qs_env,&
               calculate_forces=.FALSE.,error=error)
       END DO

       CALL pw_pool_give_back_coeff(auxbas_pw_pool,v_hartree_gspace,&
            error=error)
       CALL pw_pool_give_back_coeff(auxbas_pw_pool,v_hartree_rspace,&
            error=error)
       DO ispin=1,nspins
          CALL pw_pool_give_back_coeff(auxbas_pw_pool,v_rspace_new(ispin),&
               error=error)
       END DO
       CALL pw_pools_give_back_coeffs(pw_pools,mgrid_gspace,&
            error=error)
       CALL pw_pools_give_back_coeffs(pw_pools,mgrid_rspace,&
            error=error)
       DEALLOCATE(v_rspace_new, stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       DEALLOCATE(mgrid_rspace, mgrid_gspace, stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF

    CALL timestop(0.0_wp,handle)
  END SUBROUTINE kpp1_calc_k_p_p1
!***************************************************************************

!!****f* qs_kpp1_env_methods/kpp1_check_i_alloc [1.0] *
!!
!!   NAME
!!     kpp1_check_i_alloc
!!
!!   FUNCTION
!!     checks that the intenal storage is allocated, and allocs it if needed
!!
!!   NOTES
!!     private routine
!!
!!   INPUTS
!!     - kpp1_env: the environment to check
!!     - qs_env: the qs environment this kpp1_env lives in
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     10.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE kpp1_check_i_alloc(kpp1_env, qs_env, error)
    TYPE(qs_kpp1_env_type), POINTER :: kpp1_env
    TYPE(qs_environment_type), INTENT(in), TARGET :: qs_env
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    INTEGER :: stat, ispin, nspins
    CHARACTER(len=*), PARAMETER :: routineN='kpp1_check_i_alloc',&
         routineP=moduleN//':'//routineN
    TYPE(dft_control_type), POINTER :: dft_control
    TYPE(pw_env_type), POINTER :: pw_env
    TYPE(pw_pool_type), POINTER :: auxbas_pw_pool
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: s
    type(qs_rho_type), pointer :: rho
    type(coeff_type), dimension(:), pointer :: rho_r, rho_g

    failure=.FALSE.
    NULLIFY(dft_control,pw_env,auxbas_pw_pool,s,rho,rho_r,rho_g)

    CPPrecondition(ASSOCIATED(kpp1_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(kpp1_env%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       IF (.NOT.ASSOCIATED(kpp1_env%v_rspace).OR.&
            .NOT.ASSOCIATED(kpp1_env%v_ao)) THEN
! force rebuild if allocs something
          kpp1_env%updates_since_rebuild=-1
       END IF

       CALL get_qs_env(qs_env, pw_env=pw_env, dft_control=dft_control,&
            s=s,error=error, rho=rho)
       nspins=dft_control%nspins
       CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, error=error)

       IF (.NOT.ASSOCIATED(kpp1_env%v_rspace)) THEN
          ALLOCATE(kpp1_env%v_rspace(nspins),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          IF (.NOT.failure) THEN
             DO ispin=1,nspins
                CALL pw_pool_init_coeff(auxbas_pw_pool, &
                     kpp1_env%v_rspace(ispin),&
                     use_data=REALDATA3D, in_space=REALSPACE,error=error)
             END DO
          END IF
       END IF

       IF (.NOT.ASSOCIATED(kpp1_env%v_ao)) THEN
          ALLOCATE(kpp1_env%v_ao(nspins),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          IF (.NOT. failure) THEN
             DO ispin=1,nspins
                NULLIFY(kpp1_env%v_ao(ispin)%matrix)
                CALL replicate_matrix_structure(s(1)%matrix,&
                     kpp1_env%v_ao(ispin)%matrix,&
                     "kpp1%v_ao-"//ADJUSTL(cp_to_string(ispin)))
             END DO
          END IF
       END IF
       
       if (.not.(associated(kpp1_env%deriv_xc).and.&
            associated(kpp1_env%spin_pot).and.&
            associated(kpp1_env%grad_pot).and.&
            associated(kpp1_env%ndiag_term))) then
          call cp_assert(.not.associated(kpp1_env%deriv_xc).and.&
               .not.associated(kpp1_env%spin_pot).and.&
               .not.associated(kpp1_env%grad_pot).and.&
               .not.associated(kpp1_env%ndiag_term),&
               cp_failure_level,cp_assertion_failed,routineP,&
               " internal error xc_deriv related thing should be "//&
               "either all valid or all invalid "//&
               CPSourceFileRef,&
               error,failure)
          
          call qs_rho_get(rho, rho_r=rho_r, rho_g=rho_g, error=error)
          call xc_prep_2nd_deriv(rho_r=rho_r,rho_g=rho_g,&
               dft_control=dft_control,pw_pool=auxbas_pw_pool,&
               drho_r=kpp1_env%drho_r,pot=kpp1_env%deriv_xc,&
               spin_pot=kpp1_env%spin_pot,grad_pot=kpp1_env%grad_pot,&
               ndiag_term=kpp1_env%ndiag_term,error=error)
       end if
    END IF
  END SUBROUTINE kpp1_check_i_alloc
!***************************************************************************

!!****f* qs_kpp1_env_methods/kpp1_did_change [1.0] *
!!
!!   NAME
!!     kpp1_did_change
!!
!!   SYNOPSIS
!!     Subroutine kpp1_did_change(kpp1_env, s_struct_changed, grid_changed,&
!!         error)
!!       Type(qs_kpp1_env_type), Pointer:: kpp1_env
!!       Logical, Intent (IN), Optional:: s_struct_changed, grid_changed
!!       Type(cp_error_type), Optional, Intent (INOUT):: error
!!     End Subroutine kpp1_did_change
!!
!!   FUNCTION
!!     function to advise of changes either in the grids or in the the
!!     matrix structure (to update kpp1_env internal cache)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - kpp1_env: the kpp1_env
!!     - grid_changed: true if the grids have changed and have to be rebuilt
!!     - s_struct_changed: true if the structure of the s matrix has changed
!!     - psi0_changed: true if the value of psi0 (qs_env%rho) has changed
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     11.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE kpp1_did_change(kpp1_env, s_struct_changed, grid_changed,&
       psi0_changed, error)
    TYPE(qs_kpp1_env_type), POINTER :: kpp1_env
    LOGICAL, INTENT(in), OPTIONAL :: s_struct_changed, grid_changed, psi0_changed
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='kpp1_did_change',&
         routineP=moduleN//':'//routineN
    INTEGER :: i, stat
    logical :: my_s_struct_changed, my_psi0_changed

    failure=.FALSE.
    my_s_struct_changed=.false.
    my_psi0_changed=.false.

    if (present(s_struct_changed)) my_s_struct_changed=s_struct_changed
    if (present(psi0_changed)) my_psi0_changed=psi0_changed

    CPPrecondition(ASSOCIATED(kpp1_env),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPrecondition(kpp1_env%ref_count>0,cp_failure_level,routineP,error,failure)
    END IF
    IF (.not.failure) THEN
       IF (my_s_struct_changed) THEN
          kpp1_env%updates_since_rebuild=-1
          IF (ASSOCIATED(kpp1_env%v_ao)) THEN
             CALL deallocate_matrix_set(kpp1_env%v_ao)
          END IF
       END IF
       if (my_s_struct_changed.or.my_psi0_changed) then
          if (associated(kpp1_env%drho_r)) then
             deallocate(kpp1_env%drho_r, stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          end if
          if (associated(kpp1_env%deriv_xc)) then
             deallocate(kpp1_env%deriv_xc,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          end if
          if (associated(kpp1_env%spin_pot)) then
             deallocate(kpp1_env%spin_pot, stat=stat)
          end if
          if (associated(kpp1_env%grad_pot)) then
             deallocate(kpp1_env%grad_pot, stat=stat)
          end if
          if (associated(kpp1_env%ndiag_term)) then
             deallocate(kpp1_env%ndiag_term, stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          end if
       end if
       IF (PRESENT(grid_changed)) THEN
          IF (grid_changed) THEN
             kpp1_env%updates_since_rebuild=-1
             IF (ASSOCIATED(kpp1_env%v_rspace)) THEN
                DO i=1,SIZE(kpp1_env%v_rspace)
                   CALL coeff_deallocate(kpp1_env%v_rspace(i))
                END DO
                DEALLOCATE(kpp1_env%v_rspace,stat=stat)
                CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
             END IF
          END IF
       END IF
    END IF
  END SUBROUTINE kpp1_did_change
!***************************************************************************


!!****f* qs_kpp1_env_methods/xc_calc_2nd_deriv [1.0] *
!!
!!   NAME
!!     xc_calc_2nd_deriv
!!
!!   FUNCTION
!!     Does the first prepartions for the calculation of the 2nd deriv
!!     Calculates the second derivative of E_xc at rho in the direction
!!     rho1  (if you see the second derivative as bilinear form)
!!     partial_rho|_(rho=rho) partial_rho|_(rho=rho) E_xc drho(rho1)drho
!!     The other direction is still indetermined, thus it returns
!!     a potential (partial integration is performed to reduce it to
!!     function of rho, removing the dependence from its partial derivs)
!!
!!   NOTES
!!     drho_r, pot, spin_pot, grad_pot, ndiag_term, should be as returned by
!!     xc_prep_2nd_deriv.
!!     Private, I am not sure I want make this ugly interface public...
!!
!!   INPUTS
!!     - rho_r: point where to evaluate the (functional) derivative
!!       of exc in r space
!!     - drho_r(idir,ispin): derivatives of rho_r in in the x,y,z
!!       directions in r-space (should be allocated only if you use a
!!       gradient corrected xc potential)
!!     - rho1_r, rho1_g: direction of the first derivative in r and g space
!!     - pot: the derivative of the potential, already divised by the gradient
!!       if nedeed (grad_pot true).
!!     - spin_pot, grad_pot, ndiag_term: info about the content of pot
!!       (see xc_create_2nd_deriv_info)
!!     - dft_control: informations about the functional that should be
!!       derived
!!     - v_rspace_new: will contain the new potential. The array iterates
!!       on the spins. The grids must be already present.
!!     - pw_pool: the pool for the grids
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     11.2002 created [fawzi]
!!     12.2002 many bugfixes, split in prep and calc [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE xc_calc_2nd_deriv(rho_r,drho_r,rho1_r,rho1_g,&
       pot,spin_pot,grad_pot,ndiag_term,dft_control,&
       v_rspace_new,pw_pool,error)
    TYPE(coeff_type), DIMENSION(:), POINTER :: rho_r, rho1_r, rho1_g
    TYPE(pw_p_type), DIMENSION(:,:),pointer :: drho_r
    REAL ( kind=wp ) , DIMENSION (:,:), pointer :: pot
    INTEGER, DIMENSION(:,:), POINTER :: spin_pot
    LOGICAL, DIMENSION(:,:), POINTER :: grad_pot
    LOGICAL, DIMENSION(:), POINTER :: ndiag_term
    TYPE(dft_control_type), POINTER :: dft_control
    TYPE(coeff_type), DIMENSION(:), POINTER :: v_rspace_new
    TYPE(pw_pool_type), POINTER :: pw_pool
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure, gradient_f, crossterms
    CHARACTER(len=*), PARAMETER :: routineN='xc_calc_2nd_deriv',&
         routineP=moduleN//':'//routineN
    INTEGER :: handle, nspins, spin, ii, i,j,k, n(3), n_tot, idir, ipot,ntot,&
         stat
    TYPE(pw_type), POINTER :: tmp_g, tmp_r, v_deriv_g
    TYPE(pw_p_type), DIMENSION(:), ALLOCATABLE :: v_r_to_deriv, drho1_r
    INTEGER, DIMENSION (3,3) :: nd 
    INTEGER, DIMENSION (2,3) :: bo
    TYPE(cp_logger_type), POINTER :: logger

    CALL timeset(routineN,"I","",handle)

    failure=.FALSE.
    logger => cp_error_get_logger(error)
    NULLIFY(tmp_g, tmp_r, v_deriv_g)
    nd = RESHAPE ((/1,0,0,0,1,0,0,0,1/),(/3,3/))

    CPPrecondition(ASSOCIATED(rho_r),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(rho1_r),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(v_rspace_new),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(dft_control),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(pw_pool),cp_failure_level,routineP,error,failure)
    CPPrecondition(associated(spin_pot),cp_failure_level,routineP,error,failure)
    CPPrecondition(associated(grad_pot),cp_failure_level,routineP,error,failure)
    CPPrecondition(associated(ndiag_term),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       gradient_f=dft_control%gradient_functional
       nspins=dft_control%nspins
       crossterms=ANY(dft_control%ftags(2,:)).and.nspins>1

       IF (gradient_f) THEN
          CPPrecondition(ASSOCIATED(drho_r),cp_failure_level,routineP,error,failure)
          CPPrecondition(ASSOCIATED(rho1_g),cp_failure_level,routineP,error,failure)
       END IF       
       CPPrecondition(SIZE(v_rspace_new)==nspins,cp_failure_level,routineP,error,failure)
    END IF
    CALL cp_assert(dft_control%integration_grid%pw_default,&
         cp_failure_level,cp_unimplemented_error_nr,routineP,&
         "calculation on a special grid unimplemented "//&
         CPSourceFileRef,&
         error,failure)

    IF (.NOT.failure) THEN
       n(1:3) = rho_r(1)%pw%pw_grid%npts_local (1:3)
       bo = rho_r(1)%pw%pw_grid%bounds_local  
       ntot = PRODUCT ( n )
    end IF
    
    if (gradient_f) then

! alloc & calc drho1_r=<grad rho,grad rho1>
       IF (.NOT.failure) THEN
          ALLOCATE(drho1_r(nspins),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF (.NOT.failure) THEN
          CALL pw_pool_create_pw(pw_pool,tmp_g, &
               use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,&
               error=error)
          CALL pw_pool_create_pw(pw_pool,tmp_r, &
               use_data=REALDATA3D,in_space=REALSPACE,&
               error=error)
          DO spin=1,nspins
             NULLIFY(drho1_r(spin)%pw)
             CALL pw_pool_create_pw(pw_pool,drho1_r(spin)%pw, &
                  use_data=REALDATA3D, in_space=REALSPACE, &
                  error=error)
             DO idir=1,3
                CALL pw_copy ( rho1_g(spin)%pw, tmp_g )
                CALL pw_derive ( tmp_g, nd(:,idir) )
                CALL pw_transfer ( tmp_g, tmp_r )

                IF (idir==1) THEN
                   DO k = bo(1,3), bo(2,3)
                      DO j = bo(1,2), bo(2,2)
                         DO i = bo(1,1), bo(2,1)
                            drho1_r(spin)%pw%cr3d(i,j,k) = &
                                 tmp_r%cr3d(i,j,k)*&
                                 drho_r(idir,spin)%pw%cr3d(i,j,k)
                         END DO
                      END DO
                   END DO
                ELSE
                   DO k = bo(1,3), bo(2,3)
                      DO j = bo(1,2), bo(2,2)
                         DO i = bo(1,1), bo(2,1)
                            drho1_r(spin)%pw%cr3d(i,j,k) = &
                                 drho1_r(spin)%pw%cr3d(i,j,k)+&
                                 tmp_r%cr3d(i,j,k)*&
                                 drho_r(idir,spin)%pw%cr3d(i,j,k)
                         END DO
                      END DO
                   END DO
                END IF
             END DO
          END DO
          CALL pw_pool_give_back_pw(pw_pool,tmp_g,&
               error=error)
          CALL pw_pool_give_back_pw(pw_pool,tmp_r,&
               error=error)
       END IF

! alloc v_r_to_deriv
       IF (.NOT.failure) THEN
          ALLOCATE(v_r_to_deriv(nspins),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF (.NOT.failure) THEN
          DO spin=1,nspins
             NULLIFY(v_r_to_deriv(spin)%pw)
             CALL pw_pool_create_pw(pw_pool,v_r_to_deriv(spin)%pw,&
                  use_data=REALDATA3D, in_space=REALSPACE,error=error)
             CALL pw_zero(v_r_to_deriv(spin)%pw)
          END DO
       END IF
    END IF

! calc v_rspace_new and v_r_to_deriv
    IF (.NOT.failure) THEN
       DO ipot=1,size(pot,2)
          IF (grad_pot(1,ipot).AND.grad_pot(2,ipot)) THEN
             CPPrecondition(gradient_f,cp_failure_level,routineP,error,failure)
             IF (ndiag_term(ipot)) THEN
                CALL cp_assert(spin_pot(1,ipot)/=spin_pot(2,ipot),&
                     cp_failure_level,cp_assertion_failed,routineP,&
                     "error inconsistent spin_pot info "//&
                     CPSourceFileRef,&
                     error,failure)

                ii=0
                DO k = bo(1,3), bo(2,3)
                   DO j = bo(1,2), bo(2,2)
                      DO i = bo(1,1), bo(2,1)
                         ii=ii+1
                         v_r_to_deriv(spin_pot(2,ipot))%pw%cr3d(i,j,k) = &
                              v_r_to_deriv(spin_pot(2,ipot))%pw%cr3d(i,j,k)+&
                              drho1_r(spin_pot(1,ipot))%pw%cr3d(i,j,k)*&
                              pot(ii,ipot)
                      END DO
                   END DO
                END DO

                ii=0
                DO k = bo(1,3), bo(2,3)
                   DO j = bo(1,2), bo(2,2)
                      DO i = bo(1,1), bo(2,1)
                         ii=ii+1
                         v_r_to_deriv(spin_pot(1,ipot))%pw%cr3d(i,j,k) = &
                              v_r_to_deriv(spin_pot(1,ipot))%pw%cr3d(i,j,k)+&
                              drho1_r(spin_pot(2,ipot))%pw%cr3d(i,j,k)*&
                              pot(ii,ipot)
                      END DO
                   END DO
                END DO
             ELSE
                CALL cp_assert(spin_pot(1,ipot)==spin_pot(2,ipot),&
                     cp_failure_level,cp_assertion_failed,routineP,&
                     "error inconsistent spin_pot info "//&
                     CPSourceFileRef,&
                     error,failure)

                ii=0
                DO k = bo(1,3), bo(2,3)
                   DO j = bo(1,2), bo(2,2)
                      DO i = bo(1,1), bo(2,1)
                         ii=ii+1
                         v_r_to_deriv(spin_pot(1,ipot))%pw%cr3d(i,j,k) = &
                              v_r_to_deriv(spin_pot(1,ipot))%pw%cr3d(i,j,k)+&
                              drho1_r(spin_pot(2,ipot))%pw%cr3d(i,j,k)*&
                              pot(ii,ipot)
                      END DO
                   END DO
                END DO
             END IF
          ELSE IF (grad_pot(1,ipot)) THEN
             CPPrecondition(ndiag_term(ipot),cp_failure_level,routineP,error,failure)

             ii=0
             DO k = bo(1,3), bo(2,3)
                DO j = bo(1,2), bo(2,2)
                   DO i = bo(1,1), bo(2,1)
                      ii=ii+1
                      v_rspace_new(spin_pot(2,ipot))%pw%cr3d(i,j,k) = &
                           v_rspace_new(spin_pot(2,ipot))%pw%cr3d(i,j,k)+&
                           drho1_r(spin_pot(1,ipot))%pw%cr3d(i,j,k)*&
                           pot(ii,ipot)
                   END DO
                END DO
             END DO

             ii=0
             DO k = bo(1,3), bo(2,3)
                DO j = bo(1,2), bo(2,2)
                   DO i = bo(1,1), bo(2,1)
                      ii=ii+1
                      v_r_to_deriv(spin_pot(1,ipot))%pw%cr3d(i,j,k) = &
                           v_r_to_deriv(spin_pot(1,ipot))%pw%cr3d(i,j,k)+&
                           rho1_r(spin_pot(2,ipot))%pw%cr3d(i,j,k)*&
                           pot(ii,ipot)
                   END DO
                END DO
             END DO
          ELSE IF (grad_pot(2,ipot)) THEN
             CPPrecondition(ndiag_term(ipot),cp_failure_level,routineP,error,failure)

             ii=0
             DO k = bo(1,3), bo(2,3)
                DO j = bo(1,2), bo(2,2)
                   DO i = bo(1,1), bo(2,1)
                      ii=ii+1
                      v_r_to_deriv(spin_pot(2,ipot))%pw%cr3d(i,j,k) = &
                           v_r_to_deriv(spin_pot(2,ipot))%pw%cr3d(i,j,k)+&
                           rho1_r(spin_pot(1,ipot))%pw%cr3d(i,j,k)*&
                           pot(ii,ipot)
                   END DO
                END DO
             END DO

             ii=0
             DO k = bo(1,3), bo(2,3)
                DO j = bo(1,2), bo(2,2)
                   DO i = bo(1,1), bo(2,1)
                      ii=ii+1
                      v_rspace_new(spin_pot(1,ipot))%pw%cr3d(i,j,k) = &
                           v_rspace_new(spin_pot(1,ipot))%pw%cr3d(i,j,k)+&
                           drho1_r(spin_pot(2,ipot))%pw%cr3d(i,j,k)*&
                           pot(ii,ipot)
                   END DO
                END DO
             END DO
          ELSE
! no grad
             IF (ndiag_term(ipot)) THEN
                CALL cp_assert(spin_pot(1,ipot)/=spin_pot(2,ipot),&
                     cp_failure_level,cp_assertion_failed,routineP,&
                     "error inconsistent spin_pot info "//&
                     CPSourceFileRef,&
                     error,failure)

                ii=0
                DO k = bo(1,3), bo(2,3)
                   DO j = bo(1,2), bo(2,2)
                      DO i = bo(1,1), bo(2,1)
                         ii=ii+1
                         v_rspace_new(spin_pot(2,ipot))%pw%cr3d(i,j,k) = &
                              v_rspace_new(spin_pot(2,ipot))%pw%cr3d(i,j,k)+&
                              rho1_r(spin_pot(1,ipot))%pw%cr3d(i,j,k)*&
                              pot(ii,ipot)
                      END DO
                   END DO
                END DO

                ii=0
                DO k = bo(1,3), bo(2,3)
                   DO j = bo(1,2), bo(2,2)
                      DO i = bo(1,1), bo(2,1)
                         ii=ii+1
                         v_rspace_new(spin_pot(1,ipot))%pw%cr3d(i,j,k) = &
                              v_rspace_new(spin_pot(1,ipot))%pw%cr3d(i,j,k)+&
                              rho1_r(spin_pot(2,ipot))%pw%cr3d(i,j,k)*&
                              pot(ii,ipot)
                      END DO
                   END DO
                END DO
             ELSE
                CALL cp_assert(spin_pot(1,ipot)==spin_pot(2,ipot),&
                     cp_failure_level,cp_assertion_failed,routineP,&
                     "error inconsistent spin_pot info "//&
                     CPSourceFileRef,&
                     error,failure)
                ii=0
                DO k = bo(1,3), bo(2,3)
                   DO j = bo(1,2), bo(2,2)
                      DO i = bo(1,1), bo(2,1)
                         ii=ii+1
                         v_rspace_new(spin_pot(1,ipot))%pw%cr3d(i,j,k) = &
                              v_rspace_new(spin_pot(1,ipot))%pw%cr3d(i,j,k)+&
                              rho1_r(spin_pot(2,ipot))%pw%cr3d(i,j,k)*&
                              pot(ii,ipot)
                      END DO
                   END DO
                END DO
             END IF
          END IF
       END DO
    END IF

    IF (gradient_f) THEN
       IF (.NOT.failure) THEN
! give back drho1_r
          DO spin=1,nspins
             CALL pw_pool_give_back_pw(pw_pool,drho1_r(spin)%pw,&
                  error=error)
          END DO
          DEALLOCATE(drho1_r,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
    END IF

! perform partial integration and filter v_rspace_new
    IF (.NOT.failure) THEN
       IF (gradient_f) THEN
          CALL pw_pool_create_pw(pw_pool,tmp_g, &
               use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,&
               error=error)
       END IF
       CALL pw_pool_create_pw(pw_pool,v_deriv_g, &
            use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,&
            error=error)

       DO spin=1,nspins
          CALL pw_transfer(v_rspace_new(spin)%pw,v_deriv_g)
          
          IF (gradient_f) THEN
             DO idir=1,3
                
                DO k = bo(1,3), bo(2,3)
                   DO j = bo(1,2), bo(2,2)
                      DO i = bo(1,1), bo(2,1)
                         v_rspace_new(spin)%pw%cr3d(i,j,k)= &
                              -v_r_to_deriv(spin)%pw%cr3d(i,j,k)* &
                              drho_r(idir,spin)%pw%cr3d(i,j,k)
                      END DO
                   END DO
                END DO
                
                CALL pw_transfer(v_rspace_new(spin)%pw,tmp_g)
                CALL pw_derive(tmp_g,nd(:,idir))
                CALL pw_sumup(tmp_g,v_deriv_g)
             END DO
          END IF

          CALL pw_transfer(v_deriv_g,v_rspace_new(spin)%pw)
       END DO

       IF (gradient_f) THEN
          CALL pw_pool_give_back_pw(pw_pool,tmp_g,error=error)
       END IF
       CALL pw_pool_give_back_pw(pw_pool,v_deriv_g,error=error)

       IF (gradient_f) THEN
! dealloc v_r_to_deriv
          DO spin=1,nspins
             CALL pw_pool_give_back_pw(pw_pool,v_r_to_deriv(spin)%pw,&
                  error=error)
          END DO
          DEALLOCATE(v_r_to_deriv, stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
       
    END IF

    CALL timestop(0.0_wp,handle)
  END SUBROUTINE xc_calc_2nd_deriv
!***************************************************************************

!!****f* qs_kpp1_env_methods/xc_prep_2nd_deriv [1.0] *
!!
!!   NAME
!!     xc_prep_2nd_deriv
!!
!!   FUNCTION
!!     Does the first prepartions for the calculation of the 2nd deriv
!!     done by xc_calc_2nd_deriv
!!
!!   NOTES
!!     private, I am not sure I want make this ugly interface public...
!!
!!   INPUTS
!!     - rho_r,rho_g: point where to evaluate the (functional) derivative
!!       of exc in r and g space
!!     - dft_control: informations about the functional that should be
!!       derived
!!     - drho_r(idir,ispin) :: derivative of rho_r in the direction x,y,z
!!       in the r space (allocated if the xc potential is gradient corrected)
!!     - pot: the potential (allocated)
!!     - spin_pot, grad_pot, ndiag_term: info about the content of pot
!!       (allocated, see xc_create_2nd_deriv_info)
!!     - pw_pool: the pool for the grids
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     11.2002 created [fawzi]
!!     12.2002 many bugfixes, split in prep and calc [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE xc_prep_2nd_deriv(rho_r,rho_g,dft_control,pw_pool,&
       drho_r,pot,spin_pot,grad_pot,ndiag_term,error)
    TYPE(coeff_type), DIMENSION(:), POINTER :: rho_r, rho_g
    TYPE(dft_control_type), POINTER :: dft_control
    TYPE(pw_p_type), DIMENSION(:,:),pointer :: drho_r
    REAL ( kind=wp ) , DIMENSION (:,:), pointer :: pot
    INTEGER, DIMENSION(:,:), POINTER :: spin_pot
    LOGICAL, DIMENSION(:,:), POINTER :: grad_pot
    LOGICAL, DIMENSION(:), POINTER :: ndiag_term
    TYPE(pw_pool_type), POINTER :: pw_pool
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure, gradient_f, crossterms
    CHARACTER(len=*), PARAMETER :: routineN='xc_prep_2nd_deriv',&
         routineP=moduleN//':'//routineN
    INTEGER :: handle, nspins, spin, ii, i,j,k, n(3), n_tot, idir, ipot,ntot,&
         stat, order,npot
    TYPE(pw_type), POINTER :: tmp_g
    REAL ( kind=wp ) , DIMENSION (:,:), ALLOCATABLE :: rgrid, drgrid
    INTEGER, DIMENSION (3,3) :: nd 
    INTEGER, DIMENSION (2,3) :: bo
    TYPE(cp_logger_type), POINTER :: logger

    CALL timeset(routineN,"I","",handle)

    failure=.FALSE.
    logger => cp_error_get_logger(error)
    NULLIFY(tmp_g)
    nullify(drho_r, pot, spin_pot, grad_pot, ndiag_term)
    nd = RESHAPE ((/1,0,0,0,1,0,0,0,1/),(/3,3/))

    CPPrecondition(ASSOCIATED(rho_r),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(dft_control),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(pw_pool),cp_failure_level,routineP,error,failure)

    IF (.NOT. failure) THEN
       gradient_f=dft_control%gradient_functional
       nspins=dft_control%nspins
       crossterms=ANY(dft_control%ftags(2,:)).and.nspins>1

       IF (gradient_f) THEN
          CPPrecondition(ASSOCIATED(rho_g),cp_failure_level,routineP,error,failure)
       END IF       
    END IF
    CALL cp_assert(dft_control%integration_grid%pw_default,&
         cp_failure_level,cp_unimplemented_error_nr,routineP,&
         "calculation on a special grid unimplemented "//&
         CPSourceFileRef,&
         error,failure)

    IF (.NOT.failure) THEN
       n(1:3) = rho_r(1)%pw%pw_grid%npts_local (1:3)
       bo = rho_r(1)%pw%pw_grid%bounds_local  
       ntot = PRODUCT ( n )

       ALLOCATE(rgrid(ntot,nspins),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       IF (gradient_f) THEN
          ALLOCATE(drgrid(ntot,nspins),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          ALLOCATE(drho_r(3,nspins),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
    END IF
    IF (.NOT.failure) THEN

! set rgrid
       ii=0
       DO spin=1,nspins
          DO k = bo(1,3), bo(2,3)
             DO j = bo(1,2), bo(2,2)
                DO i = bo(1,1), bo(2,1)
                   ii = ii + 1
                   rgrid(ii,spin) = rho_r(spin)%pw%cr3d(i,j,k)
                END DO
             END DO
          END DO
       END DO

       IF (gradient_f) THEN
! calc drho_r=grad rho
          CALL pw_pool_create_pw(pw_pool,tmp_g, &
               use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,&
               error=error)
          DO spin=1,nspins
             DO idir=1,3
                NULLIFY(drho_r(idir,spin)%pw)
                CALL pw_pool_create_pw(pw_pool,drho_r(i,spin)%pw, &
                     use_data=REALDATA3D, in_space=REALSPACE, &
                     error=error)
                CALL pw_copy ( rho_g(spin)%pw, tmp_g )
                CALL pw_derive ( tmp_g, nd(:,idir) )
                CALL pw_transfer ( tmp_g, drho_r(i,spin)%pw )
             END DO
          END DO
          CALL pw_pool_give_back_pw(pw_pool,tmp_g,&
               error=error)

! calc drgrid=||grad rho||
          DO spin=1,nspins
             ii = 0
             DO k = bo(1,3), bo(2,3)
                DO j = bo(1,2), bo(2,2)
                   DO i = bo(1,1), bo(2,1)
                      ii = ii + 1
                      drgrid(ii,spin) = SQRT(drho_r(1,spin)%pw%cr3d(i,j,k)**2+&
                           drho_r(2,spin)%pw%cr3d(i,j,k)**2+&
                           drho_r(3,spin)%pw%cr3d(i,j,k)**2)
                   END DO
                END DO
             END DO
          END DO
       END IF
    END IF

    IF (.NOT.failure) THEN
       order = -2
! allocate potential
       npot = xc_get_pot_size ( nspins>1, gradient_f, crossterms, order )
       ALLOCATE ( pot ( ntot, npot ), STAT = stat )
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF

    IF (.NOT.failure) THEN
       pot = 0._wp
! calculate functionals
       CALL xc_set ( dft_control%density_cut, dft_control%gradient_cut )
       SELECT CASE ( nspins )
       CASE (1)
          IF (gradient_f) THEN
             CALL xc_lda ( functional=dft_control%functional,&
                  ftags=dft_control%ftags, &
                  rho=rgrid(:,1), drho=drgrid(:,1), pot=pot, &
                  order=order, error=error )
          ELSE
             CALL xc_lda ( functional=dft_control%functional,&
                  ftags=dft_control%ftags, &
                  rho=rgrid(:,1), pot=pot, &
                  order=order, error=error )
          END IF
       CASE (2)
          CALL cp_unimplemented_error(fromWhere=routineP,&
               message='lsd not yet implemented '//&
               CPSourceFileRef,&
               error=error)
       CASE default
          CALL cp_unimplemented_error(fromWhere=routineP,&
               message='ONLY 1 or 2 spin supported '//&
               CPSourceFileRef,&
               error=error)
       END SELECT

       CALL xc_create_2nd_deriv_info(lsd=nspins>1,gradient=gradient_f,&
            crossterms=crossterms,spin_pot=spin_pot,grad_pot=grad_pot,&
            ndiag_term=ndiag_term,error=error)
    END IF

! dealloc rgrid
    IF (.not.failure) THEN
       DEALLOCATE(rgrid,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF

! divide the potential by the gradient where necessary
    IF (.NOT.failure) THEN
       DO ipot=1,npot
          IF (grad_pot(1,ipot).AND.grad_pot(2,ipot)) THEN
             CPPrecondition(gradient_f,cp_failure_level,routineP,error,failure)
             IF (ndiag_term(ipot)) THEN
                CALL cp_assert(spin_pot(1,ipot)/=spin_pot(2,ipot),&
                     cp_failure_level,cp_assertion_failed,routineP,&
                     "error inconsistent spin_pot info "//&
                     CPSourceFileRef,&
                     error,failure)

                DO ii=1,ntot
                   if (drgrid(ii,spin_pot(1,ipot))<dft_control%gradient_cut.or.&
                        drgrid(ii,spin_pot(2,ipot))<dft_control%gradient_cut) then
                      pot(ii,ipot)=0.0_wp
                   else
                      pot(ii,ipot)=pot(ii,ipot)/(drgrid(ii,spin_pot(1,ipot))*&
                           drgrid(ii,spin_pot(2,ipot)))
                   end if
                END DO
             ELSE
                CALL cp_assert(spin_pot(1,ipot)==spin_pot(2,ipot),&
                     cp_failure_level,cp_assertion_failed,routineP,&
                     "error inconsistent spin_pot info "//&
                     CPSourceFileRef,&
                     error,failure)

                DO ii=1,ntot
                   if (drgrid(ii,spin_pot(1,ipot))<dft_control%gradient_cut) then
                      pot(ii,ipot)=0.0_wp
                   else
                      pot(ii,ipot)=pot(ii,ipot)/drgrid(ii,spin_pot(1,ipot))**2
                   end if
                END DO
              END IF
          ELSE IF (grad_pot(1,ipot)) THEN
             CPPrecondition(gradient_f,cp_failure_level,routineP,error,failure)
             CPPrecondition(ndiag_term(ipot),cp_failure_level,routineP,error,failure)
             DO ii=1,ntot
                if (drgrid(ii,spin_pot(1,ipot))<dft_control%gradient_cut.or.&
                     drgrid(ii,spin_pot(2,ipot))<dft_control%gradient_cut) then
                   pot(ii,ipot)=0.0_wp
                else
                   pot(ii,ipot)=pot(ii,ipot)/(drgrid(ii,spin_pot(1,ipot))*&
                        drgrid(ii,spin_pot(2,ipot)))
                end if
             END DO
          ELSE IF (grad_pot(2,ipot)) THEN
             CPPrecondition(gradient_f,cp_failure_level,routineP,error,failure)
             CPPrecondition(ndiag_term(ipot),cp_failure_level,routineP,error,failure)

             DO ii=1,ntot
                if (drgrid(ii,spin_pot(1,ipot))<dft_control%gradient_cut) then
                   pot(ii,ipot)=0.0_wp
                else
                   pot(ii,ipot)=pot(ii,ipot)/drgrid(ii,spin_pot(2,ipot))**2
                end if
             END DO
          ELSE
! no grad
             IF (ndiag_term(ipot)) THEN
                CALL cp_assert(spin_pot(1,ipot)/=spin_pot(2,ipot),&
                     cp_failure_level,cp_assertion_failed,routineP,&
                     "error inconsistent spin_pot info "//&
                     CPSourceFileRef,&
                     error,failure)
             ELSE
                CALL cp_assert(spin_pot(1,ipot)==spin_pot(2,ipot),&
                     cp_failure_level,cp_assertion_failed,routineP,&
                     "error inconsistent spin_pot info "//&
                     CPSourceFileRef,&
                     error,failure)
             END IF
          END IF
       END DO
    END IF

    IF (gradient_f) THEN
       IF (.NOT.failure) THEN
! dealloc drgrid
          DEALLOCATE(drgrid,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
    END IF

    CALL timestop(0.0_wp,handle)
  END SUBROUTINE xc_prep_2nd_deriv
!***************************************************************************

END MODULE qs_kpp1_env_methods
