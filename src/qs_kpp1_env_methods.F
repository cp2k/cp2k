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
  USE coefficient_types,               ONLY: coeff_allocate,&
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
  USE kinds,                           ONLY: dp
  USE l_utils,                         ONLY: l_info_type
  USE qs_ks_methods,                   ONLY: qs_ks_build_kohn_sham_matrix
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
  USE qs_collocate_density,            ONLY: calculate_total_rho
  USE qs_energy_types,                 ONLY: allocate_qs_energy,&
                                             deallocate_qs_energy,&
                                             qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_integrate_potential,          ONLY: integrate_v_core_rspace,&
                                             integrate_v_rspace
  USE qs_kpp1_env_types,               ONLY: kpp1_release,&
                                             kpp1_retain,&
                                             qs_kpp1_env_p_type,&
                                             qs_kpp1_env_type
  USE qs_rho_methods,                  ONLY: qs_rho_update_rho
  USE qs_rho_types,                    ONLY: qs_rho_did_change,&
                                             qs_rho_get,&
                                             qs_rho_type
  USE qs_tddfpt_types,                 ONLY: tddfpt_control_type
  USE simulation_cell,                 ONLY: cell_type,&
                                             get_cell
  USE sparse_matrix_types,             ONLY: add_matrices,&
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
                                             transfer_matrix,&
                                             cp_sm_scale_and_add
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc,                              ONLY: xc_prep_2nd_deriv1, xc_calc_2nd_deriv1
  IMPLICIT NONE

  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='qs_kpp1_env_methods'
  INTEGER, PRIVATE, SAVE :: last_kpp1_id_nr=0

  PUBLIC :: kpp1_create, &
            kpp1_calc_k_p_p1, &
            kpp1_calc_k_p_p1_fdiff, &
            kpp1_did_change

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
       kpp1_env%print_count=0
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
    
    IMPLICIT NONE

    TYPE(qs_kpp1_env_type), POINTER                  :: kpp1_env
    TYPE(qs_environment_type), INTENT(in), TARGET      :: qs_env
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER    :: k_p_p1
    TYPE(qs_rho_type), POINTER                         :: rho, rho1
    LOGICAL, intent(in), OPTIONAL                      :: rebuild
    TYPE(cp_error_type), OPTIONAL, INTENT(inout)       :: error

    LOGICAL :: failure, my_rebuild
    INTEGER :: handle, nspins, ispin, stat, output_unit=6, ntot
    CHARACTER(len=*), PARAMETER :: routineN='kpp1_calc_k_p_p1',&
         routineP=moduleN//':'//routineN
    REAL(KIND = dp) :: total_rho_gspace, energy_hartree, scale_fac
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

       IF (.NOT.(dft_control%do_tddfpt_calculation .AND. &
           dft_control%tddfpt_control%restricted_excitations_type == "TRIPLET" )) THEN
          CALL init_pw_poisson_solver(cell%green,error=error)
          CALL calculate_hartree(rho1_tot_gspace, &
                                 energy_hartree, & 
                                 v_hartree_gspace)
          CALL cleanup_pw_poisson_solver(cell%green,error=error)
          CALL coeff_transform_space(v_hartree_gspace,&
                                     v_hartree_rspace)
       END IF

       CALL pw_pool_give_back_coeff(auxbas_pw_pool,rho1_tot_gspace,&
            error=error)

! *** calculate the xc potential ***       

       DO ispin=1,nspins
          CALL pw_pool_init_coeff(auxbas_pw_pool,v_rspace_new(ispin),&
               use_data = REALDATA3D,&
               in_space = REALSPACE, error=error)
          CALL coeff_zero(v_rspace_new(ispin))
       END DO


!!!       IF (cp_would_output(logger,outputName="kpp1_energies",&
!!!            fromWhere=routineP,iter=200*(kpp1_env%id_nr-1)+kpp1_env%iter,&
!!!            local=.FALSE.,error=error)) THEN
!!!          IF (logger%para_env%mepos==logger%para_env%source) THEN
!!!             output_unit=cp_unitnr_for_output(logger,&
!!!                  outputName="kpp1_energies",&
!!!                  fromWhere=routineP, &
!!!                  iter=200*(kpp1_env%id_nr-1)+kpp1_env%iter,&
!!!                  local=.FALSE.,error=error)
!!!             WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
!!!                  "KPP1 Hartree energy:                           ",&
!!!                  energy_hartree
!!!             CALL cp_finished_output(logger, outputName="kpp1_energies",&
!!!                  fromWhere=routineP, &
!!!                  iter=200*(kpp1_env%id_nr-1)+kpp1_env%iter, &
!!!                  unit_nr=output_unit,local=.FALSE.,error=error)
!!!          END IF
!!!       END IF

       CALL qs_rho_get(rho,rho_r=rho_r,rho_g=rho_g,error=error)
       CALL qs_rho_get(rho1,rho_r=rho1_r,rho_g=rho1_g,error=error)
       CALL xc_calc_2nd_deriv1(rho_r=rho_r,drho_r=kpp1_env%drho_r,&
            rho1_r=rho1_r,&
            rho1_g=rho1_g,&
            dft_control=dft_control,&
            pot=kpp1_env%deriv_xc,&
            spin_pot=kpp1_env%spin_pot,grad_pot=kpp1_env%grad_pot,&
            ndiag_term=kpp1_env%ndiag_term,&
            v_rspace_new=v_rspace_new,&
            pw_pool=auxbas_pw_pool,&
            error=error)


       !-------------------------------!
       ! Add both hartree and xc terms !
       !-------------------------------!
       DO ispin=1,nspins

          ntot = SIZE(v_rspace_new(ispin)%pw%cr3d,1) * &
                 SIZE(v_rspace_new(ispin)%pw%cr3d,2) * &
                 SIZE(v_rspace_new(ispin)%pw%cr3d,3)


          CALL dscal(ntot, v_rspace_new(ispin)%pw%pw_grid%dvol, &
                           v_rspace_new(ispin)%pw%cr3d,1)

          IF (dft_control%do_tddfpt_calculation .AND. nspins==1) THEN
             IF (dft_control%tddfpt_control%restricted_excitations_type == "SINGLET") THEN
                ! add 2 times the hartree term in the singlet case
                CALL daxpy(ntot, 2.0_dp*v_hartree_rspace%pw%pw_grid%dvol, &
                     v_hartree_rspace%pw%cr3d, 1, v_rspace_new(ispin)%pw%cr3d, 1)
             END IF
             ! do not add the hartree term for triplets
          ELSE
             ! just add it if not tddfpt
             CALL daxpy(ntot, v_hartree_rspace%pw%pw_grid%dvol, &
                  v_hartree_rspace%pw%cr3d, 1, v_rspace_new(ispin)%pw%cr3d, 1)
          END IF

          ! remove kpp1_env%v_ao and work directly on k_p_p1 ?
          IF (my_rebuild) CALL set_matrix(kpp1_env%v_ao(ispin)%matrix,0.0_dp)
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
          call transfer_matrix(kpp1_env%v_ao(ispin)%matrix,&
               k_p_p1(ispin)%matrix)
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

    CALL timestop(0.0_dp,handle)
  END SUBROUTINE kpp1_calc_k_p_p1
!***************************************************************************


!!****f* qs_p_utils/p_calc_k_p_p1_fdiff [1.0] *
!!
!!   NAME
!!     p_calc_k_p_p1_fdiff
!!
!!   FUNCTION
!!     calcualtes the k_p_p1 kernel of the pertubation theory with finite
!!     differences
!!
!!   NOTES
!!     useful for testing purposes.
!!     rescale my_diff depending on the norm of rho1?
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
!!     - diff: the amount of the finite difference step
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     01.2003 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE kpp1_calc_k_p_p1_fdiff(kpp1_env,qs_env,k_p_p1,rho,rho1,&
       rebuild, diff, error)
    TYPE(qs_kpp1_env_type), POINTER                    :: kpp1_env
    TYPE(qs_environment_type), INTENT(inout), TARGET   :: qs_env
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER    :: k_p_p1
    TYPE(qs_rho_type), POINTER                         :: rho, rho1
    LOGICAL, intent(in), OPTIONAL                      :: rebuild
    real(KIND = dp), intent(in), optional                 :: diff
    TYPE(cp_error_type), OPTIONAL, INTENT(inout)       :: error

    logical :: failure
    CHARACTER(len=*), PARAMETER :: routineN='kpp1_calc_k_p_p1_fdiff',&
         routineP=moduleN//':'//routineN
    type(real_matrix_p_type), dimension(:), pointer :: ks_2, s, rho_ao
    type(dft_control_type), pointer :: dft_control
    type(qs_energy_type), pointer :: qs_energy
    integer :: stat, ispin
    real(KIND = dp):: scaling, my_diff

    nullify(ks_2,s,dft_control,qs_energy)
    my_diff=1.0e-6
    if (present(diff)) my_diff=diff
    call allocate_qs_energy(qs_energy)

    call qs_rho_get(rho, rho_ao=rho_ao,error=error)
    call get_qs_env(qs_env, dft_control=dft_control,&
         s=s,error=error)

    ! rho = rho0+h/2*rho1
    my_diff=my_diff/2.0_dp
    do ispin=1,dft_control%nspins
       call cp_sm_scale_and_add(rho%rho_ao(ispin)%matrix,&
            matrix_b=rho1%rho_ao(ispin)%matrix,beta=my_diff)
       rho%rho_r(ispin)%pw%cr3d=rho%rho_r(ispin)%pw%cr3d+&
            my_diff*rho1%rho_r(ispin)%pw%cr3d
       rho%rho_g(ispin)%pw%cc=rho%rho_g(ispin)%pw%cc+&
            my_diff*rho1%rho_g(ispin)%pw%cc
    end do

    call qs_ks_build_kohn_sham_matrix(ks_env=qs_env%ks_env,&
         qs_env=qs_env,ks=k_p_p1,rho=rho,energy=qs_energy,&
         rebuild=.true.,calculate_forces=.false.,&
         just_energy=.false.,error=error)

    allocate(ks_2(dft_control%nspins),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    do ispin=1,dft_control%nspins
       nullify(ks_2(ispin)%matrix)
       call replicate_matrix_structure(s(1)%matrix,ks_2(ispin)%matrix,&
            target_name="tmp_ks2-"//adjustl(cp_to_string(ispin)))
    end do

    ! rho = rho0-h/2*rho1
    my_diff=-2.0_dp*my_diff
    do ispin=1,dft_control%nspins
       call cp_sm_scale_and_add(rho%rho_ao(ispin)%matrix,&
            matrix_b=rho1%rho_ao(ispin)%matrix,beta=my_diff)
       rho%rho_r(ispin)%pw%cr3d=rho%rho_r(ispin)%pw%cr3d+&
            my_diff*rho1%rho_r(ispin)%pw%cr3d
       rho%rho_g(ispin)%pw%cc=rho%rho_g(ispin)%pw%cc+&
            my_diff*rho1%rho_g(ispin)%pw%cc
    end do

    call qs_ks_build_kohn_sham_matrix(ks_env=qs_env%ks_env,&
         qs_env=qs_env,ks=ks_2,rho=rho,energy=qs_energy,&
         rebuild=.true.,calculate_forces=.false.,&
         just_energy=.false.,error=error)

    ! rho = rho0
    my_diff=-0.5_dp*my_diff
    do ispin=1,dft_control%nspins
       call cp_sm_scale_and_add(rho%rho_ao(ispin)%matrix,&
            matrix_b=rho1%rho_ao(ispin)%matrix,beta=my_diff)
       rho%rho_r(ispin)%pw%cr3d=rho%rho_r(ispin)%pw%cr3d+&
            my_diff*rho1%rho_r(ispin)%pw%cr3d
       rho%rho_g(ispin)%pw%cc=rho%rho_g(ispin)%pw%cc+&
            my_diff*rho1%rho_g(ispin)%pw%cc
    end do

    ! k_p_p1=(H(rho0+h/2 rho1)-H(rho0-h/2 rho1))/h
    do ispin=1,dft_control%nspins
       call cp_sm_scale_and_add(matrix_a=k_p_p1(ispin)%matrix, &
            matrix_b=ks_2(ispin)%matrix,&
            beta=-1.0_dp)
       call cp_sm_scale_and_add(alpha=0.5_dp/my_diff,&
            matrix_a=k_p_p1(ispin)%matrix)
    end do

    call deallocate_matrix_set(ks_2)
    call deallocate_qs_energy(qs_energy)
  end SUBROUTINE kpp1_calc_k_p_p1_fdiff

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
!!     04.2003 add optional dft_control [tchassai]
!!
!!*** **********************************************************************
  SUBROUTINE kpp1_check_i_alloc(kpp1_env, qs_env, error)

    IMPLICIT NONE

    TYPE(qs_kpp1_env_type), POINTER                 :: kpp1_env
    TYPE(qs_environment_type), INTENT(in), TARGET   :: qs_env
    TYPE(cp_error_type), OPTIONAL, INTENT(inout)    :: error
                                                    
    LOGICAL                                         :: failure
    INTEGER                                         :: stat, ispin, nspins
    TYPE(dft_control_type), POINTER                 :: dft_control
    TYPE(pw_env_type), POINTER                      :: pw_env
    TYPE(pw_pool_type), POINTER                     :: auxbas_pw_pool
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: s
    type(qs_rho_type), pointer                      :: rho
    type(coeff_type), dimension(:), pointer         :: rho_r, rho_g

    CHARACTER(len=*), PARAMETER :: routineN='kpp1_check_i_alloc',&
                                   routineP=moduleN//':'//routineN

    ! ------------------------------------------------------------------

    failure=.FALSE.
    NULLIFY(pw_env,auxbas_pw_pool,s,rho,rho_r,rho_g)

    CPPrecondition(ASSOCIATED(kpp1_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(kpp1_env%ref_count>0,cp_failure_level,routineP,error,failure)

    IF (.NOT. failure) THEN


       IF (.NOT.ASSOCIATED(kpp1_env%v_rspace).OR.&
            .NOT.ASSOCIATED(kpp1_env%v_ao)) THEN
          ! force rebuild if allocs something
          kpp1_env%updates_since_rebuild=-1
       END IF

       CALL get_qs_env(qs_env, pw_env=pw_env, dft_control=dft_control,&
                       s=s, error=error, rho=rho)

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
          
          !-----------------------------------------!
          ! now the derivative is actually computed !
          !-----------------------------------------!
          call qs_rho_get(rho, rho_r=rho_r, rho_g=rho_g, error=error)
          call xc_prep_2nd_deriv1(rho_r=rho_r, rho_g=rho_g,&
                                  dft_control=dft_control, pw_pool=auxbas_pw_pool,&
                                  drho_r=kpp1_env%drho_r, pot=kpp1_env%deriv_xc,&
                                  spin_pot=kpp1_env%spin_pot, grad_pot=kpp1_env%grad_pot,&
                                  ndiag_term=kpp1_env%ndiag_term, error=error)
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
! ***************************************************************************

END MODULE qs_kpp1_env_methods
