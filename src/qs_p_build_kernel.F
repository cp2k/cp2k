
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/qs_p_build_kernel [1.0] *
!!
!!   NAME
!!     qs_p_build_kernel
!!
!!   FUNCTION
!!     module that builds the second order perturbation kernel
!!     K_P_P1 = delta_rho|_P delta_rho|_P E drho(P1) drho
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
MODULE qs_p_build_kernel
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
  USE cp_error_handling,               ONLY: cp_assert,&
       cp_assertion_failed,&
       cp_error_type,&
       cp_unimplemented_error,&
       cp_unimplemented_error_nr,&
       cp_error_get_logger
  USE cp_fm_types,                     ONLY: cp_full_matrix_p_type,&
       cp_full_matrix_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
       cp_logger_type,&
       cp_note_level,&
       cp_to_string,&
       cp_warning_level
  use cp_output_handling, only: cp_would_output, cp_finished_output,&
       cp_unitnr_for_output
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
  USE pws,                             ONLY: init_pw_poisson_solver,&
       pw_poisson_solver
  USE qs_build_ks_matrix,              ONLY: qs_ks_did_change,&
       qs_ks_env_type,&
       qs_ks_update_qs_env
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
  USE qs_p_types,                      ONLY: qs_K_P_P1_env_type
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
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qs_p_build_kernel'
  INTEGER, PRIVATE, SAVE :: last_kpp1_id_nr=0

  PUBLIC :: kpp1_create, kpp1_retain, kpp1_release, &
       kpp1_calc_k_p_p1, calculate_2nd_xc_deriv
!!***
!****************************************************************************

CONTAINS

!!****f* qs_p_build_kernel/kpp1_create [1.0] *
!!
!!   NAME
!!     kpp1_create
!!
!!   FUNCTION
!!     allocates and initializes a K_P_P1_env
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
    TYPE(qs_K_P_P1_env_type), POINTER :: kpp1_env
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
       NULLIFY(kpp1_env%v_rspace, kpp1_env%v_ao)
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

!!****f* qs_p_build_kernel/kpp1_release [1.0] *
!!
!!   NAME
!!     kpp1_release
!!
!!   FUNCTION
!!     releases a K_P_P1_env (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - kpp1_env: the environement to release
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
  SUBROUTINE kpp1_release(kpp1_env,error)
    TYPE(qs_K_P_P1_env_type), POINTER :: kpp1_env
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='kpp1_release',&
         routineP=moduleN//':'//routineN
    INTEGER :: ispin, stat
    failure=.FALSE.

    IF (ASSOCIATED(kpp1_env)) THEN
       CPPrecondition(kpp1_env%ref_count>0,cp_failure_level,routineP,error,failure)
       kpp1_env%ref_count=kpp1_env%ref_count-1
       IF (kpp1_env%ref_count<1) THEN
          IF (ASSOCIATED(kpp1_env%v_rspace)) THEN
             DO ispin=1,SIZE(kpp1_env%v_rspace)
                CALL coeff_deallocate(kpp1_env%v_rspace(ispin))
             END DO
             DEALLOCATE(kpp1_env%v_rspace,stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
          IF (ASSOCIATED(kpp1_env%v_ao)) THEN
             DO ispin=1,SIZE(kpp1_env%v_ao)
                IF (ASSOCIATED(kpp1_env%v_ao(ispin)%matrix)) THEN
                   CALL deallocate_matrix(kpp1_env%v_ao(ispin)%matrix)
                END IF
             END DO
             DEALLOCATE(kpp1_env%v_ao, stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          END IF
       END IF
    END IF
    NULLIFY(kpp1_env)
  END SUBROUTINE kpp1_release
!***************************************************************************

!!****f* qs_p_build_kernel/kpp1_retain [1.0] *
!!
!!   NAME
!!     kpp1_retain
!!
!!   FUNCTION
!!     retains a K_P_P1_env (see doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - kpp1_env: the environement to retain
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
  SUBROUTINE kpp1_retain(kpp1_env,error)
    TYPE(qs_K_P_P1_env_type), POINTER :: kpp1_env
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='kpp1_retain',&
         routineP=moduleN//':'//routineN
    failure=.FALSE.

    CPPrecondition(ASSOCIATED(kpp1_env),cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       CPPrecondition(kpp1_env%ref_count>0,cp_failure_level,routineP,error,failure)
       kpp1_env%ref_count=kpp1_env%ref_count+1
    END IF
  END SUBROUTINE kpp1_retain
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
    TYPE(qs_K_P_P1_env_type), POINTER                  :: kpp1_env
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
    TYPE(coeff_type), DIMENSION(:), POINTER       :: v_rspace_new
    TYPE(coeff_type), DIMENSION(:), POINTER       :: mgrid_gspace, &
         mgrid_rspace
    TYPE(coeff_type)                              :: v_hartree_gspace,&
         v_hartree_rspace,&
         rho1_tot_gspace
    TYPE(cp_logger_type), POINTER                 :: logger

    CALL timeset(moduleN,"I","",handle)
    failure=.FALSE.
    logger => cp_error_get_logger(error)
    NULLIFY(dft_control,auxbas_pw_pool, pw_pools,&
         pw_env, mgrid_rspace, mgrid_gspace)

    CPPrecondition(ASSOCIATED(kpp1_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(k_p_p1),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(rho),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(rho1),cp_failure_level,routineP,error,failure)

    IF (.NOT.failure) THEN
       CPPrecondition(kpp1_env%ref_count>0,cp_failure_level,routineP,error,failure)
       CPPrecondition(SIZE(k_p_p1)>=nspins,cp_failure_level,routineP,error,failure)

       CALL get_qs_env(qs_env=qs_env,&
            dft_control=dft_control,&
            pw_env=pw_env,&
            error=error)

       nspins=dft_control%nspins
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

       CALL calculate_hartree(rho1_tot_gspace,energy_hartree,&
            v_hartree_gspace)

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

       CALL calculate_2nd_xc_deriv(rho=rho,rho1=rho1,&
            qs_env=qs_env,&
            v_rspace_new=v_rspace_new,&
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

!!****f* qs_p_build_kernel/kpp1_check_i_alloc [1.0] *
!!
!!   NAME
!!     kpp1_check_i_alloc
!!
!!   FUNCTION
!!     checks that the intenal storage is allocated, and allocs it if needed
!!
!!   NOTES
!!     -
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
    TYPE(qs_K_P_P1_env_type), POINTER :: kpp1_env
    TYPE(qs_environment_type), INTENT(in), TARGET :: qs_env
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    INTEGER :: stat, ispin, nspins
    CHARACTER(len=*), PARAMETER :: routineN='kpp1_check_i_alloc',&
         routineP=moduleN//':'//routineN
    TYPE(dft_control_type), POINTER :: dft_control
    TYPE(pw_env_type), POINTER :: pw_env
    TYPE(pw_pool_type), POINTER :: auxbas_pw_pool
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: h

    failure=.FALSE.
    NULLIFY(dft_control,pw_env,auxbas_pw_pool,h)

    IF (.NOT. failure) THEN
       IF (.NOT.ASSOCIATED(kpp1_env%v_rspace).OR.&
            .NOT.ASSOCIATED(kpp1_env%v_ao)) THEN
! force rebuild if allocs something
          kpp1_env%updates_since_rebuild=-1
       END IF

       CALL get_qs_env(qs_env, pw_env=pw_env, dft_control=dft_control,&
            h=h,error=error)
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
                CALL replicate_matrix_structure(h(1)%matrix,&
                     kpp1_env%v_ao(ispin)%matrix,&
                     "kpp1%v_ao-"//ADJUSTL(cp_to_string(ispin)))
             END DO
          END IF
       END IF
    END IF
  END SUBROUTINE kpp1_check_i_alloc
!***************************************************************************

!!****f* qs_p_build_kernel/kpp1_calculate_xc_deriv [1.0] *
!!
!!   NAME
!!     kpp1_calculate_xc_deriv
!!
!!   FUNCTION
!!     Calculates the second derivative of E_xc at rho in the direction
!!     rho1.
!!      partial_rho|_rho=rho partial_rho|rho=rho E_xc drho(rho1)drho
!!     The other direction is still indetermined, thus it returns
!!     a map (array of pws), one for each spin and (if necessary) also 
!!     gradient)
!!
!!   NOTES
!!     I tried to wrote a general code, but to be as efficent as possible,
!!     sometime readability suffered a little
!!
!!   INPUTS
!!     - rho: point where to evaluate the (functional) derivative of exc
!!     - rho1: direction of the ,qs_env,v_rspace_new,error)
!!     - qs_env: qs environment
!!     - v_rspace_new: will contain the new potential. The array iterates
!!       on the spins. The grids must be already present.
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
  SUBROUTINE calculate_2nd_xc_deriv(rho,rho1,qs_env,v_rspace_new,error)
    TYPE(qs_rho_type), POINTER :: rho
    TYPE(qs_rho_type), POINTER :: rho1
    TYPE(qs_environment_type), INTENT(in), TARGET :: qs_env
    TYPE(coeff_type), DIMENSION(:), POINTER :: v_rspace_new
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure, gradient_f, crossterms
    CHARACTER(len=*), PARAMETER :: routineN='calculate_2nd_xc_deriv',&
         routineP=moduleN//':'//routineN
    INTEGER :: handle, nspins, spin, ii, i,j,k, n(3), n_tot, idir, ipot,ntot,&
         stat, order,npot
    TYPE(pw_pool_type), POINTER :: auxbas_pw_pool
    TYPE(pw_pool_p_type), DIMENSION(:), POINTER :: pw_pools
    TYPE(pw_env_type), POINTER :: pw_env
    TYPE(pw_p_type), DIMENSION(:,:),ALLOCATABLE :: drho_r
    TYPE(pw_type), POINTER :: tmp_g, tmp_r, v_deriv_g
    TYPE(pw_p_type), DIMENSION(:), ALLOCATABLE :: v_r_to_deriv, drho1_r
    TYPE(coeff_type), DIMENSION(:), POINTER :: rho_r, rho_g, rho1_r, rho1_g
    TYPE(dft_control_type), POINTER :: dft_control
    REAL ( kind=wp ) , DIMENSION (:,:), ALLOCATABLE :: rgrid, drgrid
    INTEGER, DIMENSION (3,3) :: nd 
    INTEGER, DIMENSION (2,3) :: bo
    INTEGER, DIMENSION(:,:), POINTER :: spin_pot
    LOGICAL, DIMENSION(:,:), POINTER :: grad_pot
    LOGICAL, DIMENSION(:), POINTER :: ndiag_term
    REAL ( kind=wp ) , DIMENSION (:,:), ALLOCATABLE :: pot

    CALL timeset(routineN,"I","",handle)

    failure=.FALSE.
    NULLIFY(auxbas_pw_pool, pw_pools, pw_env, tmp_g, rho_r, rho_g, rho1_r,&
         rho1_g, spin_pot,grad_pot,ndiag_term,dft_control)
    nd = RESHAPE ((/1,0,0,0,1,0,0,0,1/),(/3,3/))
    CALL get_qs_env(qs_env,dft_control=dft_control,pw_env=pw_env,error=error)
    gradient_f=dft_control%gradient_functional
    crossterms=ANY(dft_control%ftags(2,:))
    nspins=dft_control%nspins
    CALL pw_env_get(pw_env,auxbas_pw_pool=auxbas_pw_pool,pw_pools=pw_pools,&
         error=error)

    CPPrecondition(ASSOCIATED(rho),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(rho1),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(v_rspace_new),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPrecondition(SIZE(v_rspace_new)==nspins,cp_failure_level,routineP,error,failure)
    END IF
    CALL cp_assert(dft_control%integration_grid%pw_default,&
         cp_failure_level,cp_unimplemented_error_nr,routineP,&
         "calculation on a special grid unimplemented "//&
         CPSourceFileRef,&
         error,failure)

    IF (.NOT.failure) THEN
       CALL qs_rho_get(rho,rho_r=rho_r, rho_g=rho_g, error=error)
       CALL qs_rho_get(rho1,rho_r=rho1_r, rho_g=rho1_g, error=error)
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
          CALL pw_pool_create_pw(auxbas_pw_pool,tmp_g, &
               use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,&
               error=error)
          DO spin=1,nspins
             DO idir=1,3
                NULLIFY(drho_r(idir,spin)%pw)
                CALL pw_pool_create_pw(auxbas_pw_pool,drho_r(i,spin)%pw, &
                     use_data=REALDATA3D, in_space=REALSPACE, &
                     error=error)
                CALL pw_copy ( rho_g(spin)%pw, tmp_g )
                CALL pw_derive ( tmp_g, nd(:,idir) )
                CALL pw_transfer ( tmp_g, drho_r(i,spin)%pw )
             END DO
          END DO
          CALL pw_pool_give_back_pw(auxbas_pw_pool,tmp_g,&
               error=error)

! calc drgrid=||grad rho||
          DO spin=1,nspins
             ii = 0
             DO k = bo(1,3), bo(2,3)
                DO j = bo(1,2), bo(2,2)
                   DO i = bo(1,1), bo(2,1)
                      ii = ii + 1
                      drgrid(ii,spin) = SQRT(drho_r(1,spin)%pw%cr3d(i,j,k)**2 + &
                           drho_r(2,spin)%pw%cr3d(i,j,k)**2 + &
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
          CALL xc_lda ( dft_control%functional, dft_control%ftags, &
               rgrid(:,1), drgrid(:,1), pot, order )
       CASE (2)
          CALL cp_unimplemented_error(fromWhere=routineP,&
               message='lsd not yet implemented '//&
               CPSourceFileRef,&
               error=error)
       CASE default
          CALL cp_unimplemented_error(fromWhere=routineP,&
               message='only 1 or 2 spin supported '//&
               CPSourceFileRef,&
               error=error)
       END SELECT

       CALL xc_create_2nd_deriv_info(lsd=nspins>1,gradient=gradient_f,&
            crossterms=crossterms,spin_pot=spin_pot,grad_pot=grad_pot,&
            ndiag_term=ndiag_term,error=error)
    END IF

    IF (gradient_f) THEN

! alloc & calc drho1_r=<grad rho,grad rho1>
       IF (.NOT.failure) THEN
          ALLOCATE(drho1_r(nspins),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF (.NOT.failure) THEN
          CALL pw_pool_create_pw(auxbas_pw_pool,tmp_g, &
               use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,&
               error=error)
          CALL pw_pool_create_pw(auxbas_pw_pool,tmp_r, &
               use_data=REALDATA3D,in_space=REALSPACE,&
               error=error)
          DO spin=1,nspins
             NULLIFY(drho1_r(spin)%pw)
             IF (nspins==1) THEN ! lda can use tmp_r
                drho1_r(spin)%pw => tmp_r
             ELSE
                CALL pw_pool_create_pw(auxbas_pw_pool,drho1_r(spin)%pw, &
                     use_data=REALDATA3D, in_space=REALSPACE, &
                     error=error)
             END IF
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
          CALL pw_pool_give_back_pw(auxbas_pw_pool,tmp_g,&
               error=error)
          IF (nspins==1) THEN
             NULLIFY(tmp_r)
          ELSE
             CALL pw_pool_give_back_pw(auxbas_pw_pool,tmp_r,&
                  error=error)
          END IF
       END IF

! alloc v_r_to_deriv
       IF (.NOT.failure) THEN
          ALLOCATE(v_r_to_deriv(nspins),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       END IF
       IF (.NOT.failure) THEN
          DO spin=1,nspins
             NULLIFY(v_r_to_deriv(spin)%pw)
             CALL pw_pool_create_pw(auxbas_pw_pool,v_r_to_deriv(spin)%pw,&
                  use_data=REALDATA3D, in_space=REALSPACE,error=error)
             CALL pw_zero(v_r_to_deriv(spin)%pw)
          END DO
       END IF
    END IF

! calc v_rspace_new and v_r_to_deriv
    IF (.NOT.failure) THEN
       DO ipot=1,npot
          IF (grad_pot(1,ipot).AND.grad_pot(2,ipot)) THEN
             IF (ndiag_term(ipot)) THEN
                DO ii=0,ntot
                   pot(ii,ipot)=pot(ii,ipot)/(drgrid(ii,spin_pot(1,ipot))*&
                        drgrid(ii,spin_pot(2,ipot)))
                END DO

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
                              pot(ii,ipot)/drgrid(ii,spin_pot(1,ipot))**2
                      END DO
                   END DO
                END DO
             END IF
          ELSE IF (grad_pot(1,ipot)) THEN
             CPPrecondition(ndiag_term(ipot),cp_failure_level,routineP,error,failure)
             DO ii=0,ntot
                pot(ii,ipot)=pot(ii,ipot)/drgrid(ii,spin_pot(1,ipot))
             END DO

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

             DO ii=0,ntot
                pot(ii,ipot)=pot(ii,ipot)/drgrid(ii,spin_pot(2,ipot))
             END DO

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
                              pot(ii,ipot)/drgrid(ii,spin_pot(1,ipot))**2
                      END DO
                   END DO
                END DO
             END IF
          END IF
       END DO
    END IF

! give back drho1_r
    IF (gradient_f) THEN
       IF (.NOT.failure) THEN
          DO spin=1,nspins
             CALL pw_pool_give_back_pw(auxbas_pw_pool,drho1_r(spin)%pw,&
                  error=error)
          END DO
          DEALLOCATE(drho1_r,stat=stat)
          CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
       END IF
    END IF


! dealloc pot
    IF (.NOT.failure) THEN
       DEALLOCATE(pot,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
    END IF

! perform partial integration
    IF (.NOT.failure) THEN
       CALL pw_pool_create_pw(auxbas_pw_pool,tmp_g, &
            use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,&
            error=error)
       CALL pw_pool_create_pw(auxbas_pw_pool,v_deriv_g, &
            use_data=COMPLEXDATA1D,in_space=RECIPROCALSPACE,&
            error=error)
       CALL pw_pool_create_pw(auxbas_pw_pool,tmp_r, &
            use_data=REALDATA3D,in_space=REALSPACE,&
            error=error)

       DO spin=1,nspins
          DO idir=1,3

             DO k = bo(1,3), bo(2,3)
                DO j = bo(1,2), bo(2,2)
                   DO i = bo(1,1), bo(2,1)
                      tmp_r%cr3d(i,j,k)= &
                           -v_r_to_deriv(spin)%pw%cr3d(i,j,k)* &
                           drho_r(idir,spin)%pw%cr3d(i,j,k)
                   END DO
                END DO
             END DO

             CALL pw_transfer(tmp_r,tmp_g)
             CALL pw_derive(tmp_g,nd(:,idir))
             IF (idir==1) THEN
                CALL pw_copy(tmp_g,v_deriv_g)
             ELSE
                CALL pw_sumup(v_deriv_g,tmp_g)
             END IF
          END DO
          CALL pw_transfer(v_deriv_g,tmp_r)
          CALL pw_sumup(v_rspace_new(spin)%pw,tmp_r)
       END DO

       CALL pw_pool_give_back_pw(auxbas_pw_pool,tmp_g,error=error)
       CALL pw_pool_give_back_pw(auxbas_pw_pool,v_deriv_g,error=error)
       CALL pw_pool_give_back_pw(auxbas_pw_pool,tmp_r,error=error)

! dealloc v_rspace_new and drho_r
       DO spin=1,nspins
          CALL pw_pool_give_back_pw(auxbas_pw_pool,v_rspace_new(spin)%pw,&
               error=error)
          DO idir=1,3
             CALL pw_pool_give_back_pw(auxbas_pw_pool,drho_r(idir,spin)%pw,&
                  error=error)
          END DO
       END DO
       DEALLOCATE(v_rspace_new,drho_r,stat=stat)
       CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)

    END IF

    CALL timestop(0.0_wp,handle)
  END SUBROUTINE calculate_2nd_xc_deriv
!***************************************************************************

END MODULE qs_p_build_kernel
