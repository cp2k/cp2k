
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
  USE qs_p_types,                      ONLY: kpp1_release,&
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
  USE xc,                              ONLY: xc_calc_2nd_deriv
  IMPLICIT NONE

  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qs_p_build_kernel'
  INTEGER, PRIVATE, SAVE :: last_kpp1_id_nr=0

  PUBLIC :: kpp1_create, &
       kpp1_calc_k_p_p1, kpp1_did_change
!!***
!****************************************************************************

CONTAINS

!!****f* qs_p_build_kernel/kpp1_create [1.0] *
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
    logger => cp_error_get_logger(error)
    NULLIFY(dft_control,auxbas_pw_pool, pw_pools,&
         pw_env,v_rspace_new,rho_r,rho_g,rho1_r, rho1_g,&
         mgrid_rspace, mgrid_gspace,logger,cell)

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
       CALL xc_calc_2nd_deriv(rho_r=rho_r,rho1_r=rho1_r,&
            rho_g=rho_g, rho1_g=rho1_g,&
            dft_control=dft_control,&
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

!!****f* qs_p_build_kernel/kpp1_check_i_alloc [1.0] *
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

    failure=.FALSE.
    NULLIFY(dft_control,pw_env,auxbas_pw_pool,s)

    CPPrecondition(ASSOCIATED(kpp1_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(kpp1_env%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       IF (.NOT.ASSOCIATED(kpp1_env%v_rspace).OR.&
            .NOT.ASSOCIATED(kpp1_env%v_ao)) THEN
! force rebuild if allocs something
          kpp1_env%updates_since_rebuild=-1
       END IF

       CALL get_qs_env(qs_env, pw_env=pw_env, dft_control=dft_control,&
            s=s,error=error)
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
    END IF
  END SUBROUTINE kpp1_check_i_alloc
!***************************************************************************

!!****f* qs_p_build_kernel/kpp1_did_change [1.0] *
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
SUBROUTINE kpp1_did_change(kpp1_env, s_struct_changed, grid_changed,error)
  TYPE(qs_kpp1_env_type), POINTER :: kpp1_env
  LOGICAL, INTENT(in), OPTIONAL :: s_struct_changed, grid_changed
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='kpp1_did_change',&
        routineP=moduleN//':'//routineN
  INTEGER :: i, stat

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(kpp1_env),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     CPPrecondition(kpp1_env%ref_count>0,cp_failure_level,routineP,error,failure)
  END IF
  IF (.not.failure) THEN
     IF (PRESENT(s_struct_changed)) THEN
        IF (s_struct_changed) THEN
           kpp1_env%updates_since_rebuild=-1
           IF (ASSOCIATED(kpp1_env%v_ao)) THEN
              CALL deallocate_matrix_set(kpp1_env%v_ao)
           END IF
        END IF
     END IF
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

END MODULE qs_p_build_kernel
