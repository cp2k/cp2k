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
  USE coefficient_types,               ONLY: coeff_copy,&
                                             coeff_deallocate,&
                                             coeff_init,&
                                             coeff_sumup,&
                                             coeff_transform_space,&
                                             coeff_type
  USE cp_control_types,                ONLY: dft_control_type,&
                                             xc_control_type
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_error_get_logger,&
                                             cp_error_type
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cp_output_handling,              ONLY: cp_finished_output,&
                                             cp_unitnr_for_output
  USE cp_subsystem_types,              ONLY: cp_subsystem_p_type
  USE hartree,                         ONLY: calculate_hartree
  USE hartree_local_methods,           ONLY: Vh_1c_gg_integrals
  USE hartree_local_2centers,          ONLY: integrate_vhg0_g0
  USE hartree_local_3centers,          ONLY: integrate_a_vhg0_b
  USE kinds,                           ONLY: dp
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_pool_types,                   ONLY: pw_pool_give_back_coeff,&
                                             pw_pool_init_coeff,&
                                             pw_pool_p_type,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_create,&
                                             pw_p_type,&
                                             pw_release,&
                                             pw_transfer
  USE pws,                             ONLY: cleanup_pw_poisson_solver,&
                                             init_pw_poisson_solver
  USE qs_collocate_density,            ONLY: calculate_total_rho
  USE qs_energy_types,                 ONLY: allocate_qs_energy,&
                                             deallocate_qs_energy,&
                                             qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_gapw_densities,               ONLY: prepare_gapw_den
  USE qs_integrate_potential,          ONLY: integrate_v_rspace
  USE qs_kpp1_env_types,               ONLY: qs_kpp1_env_type
  USE qs_ks_atom,                      ONLY: update_ks_atom
  USE qs_ks_methods,                   ONLY: qs_ks_build_kohn_sham_matrix
  USE qs_p_env_types,                  ONLY: qs_p_env_type
  USE qs_rho_types,                    ONLY: qs_rho_get,&
                                             qs_rho_type
  USE qs_rho0_ggrid,                   ONLY: integrate_vhg0_rspace
  USE qs_vxc_atom,                     ONLY: calculate_xc_2nd_deriv_atom
  USE simulation_cell,                 ONLY: cell_type
  USE sparse_matrix_types,             ONLY: cp_sm_scale_and_add,&
                                             deallocate_matrix_set,&
                                             real_matrix_p_type,&
                                             replicate_matrix_structure,&
                                             set_matrix,&
                                             transfer_matrix
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc,                              ONLY: xc_calc_2nd_deriv_new,&
                                             xc_prep_2nd_deriv_new
  USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                             xc_dset_release
  USE xc_derivatives,                  ONLY: xc_functional_get_info
  USE xc_rho_set_types,                ONLY: xc_rho_cflags_setall,&
                                             xc_rho_cflags_type,&
                                             xc_rho_set_create,&
                                             xc_rho_set_release,&
                                             xc_rho_set_type,&
                                             xc_rho_set_update

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
!!     Subroutine kpp1_create(kpp1_env, qs_env,  error)
!!       Type(qs_kpp1_env_type), Pointer:: kpp1_env
!!       Type(qs_environment_type), Intent (IN), Target:: qs_env
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
  SUBROUTINE kpp1_create(kpp1_env,qs_env, error)
    TYPE(qs_kpp1_env_type), POINTER          :: kpp1_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kpp1_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    ALLOCATE(kpp1_env,stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       NULLIFY(kpp1_env%v_rspace, kpp1_env%v_ao, kpp1_env%drho_r,&
            kpp1_env%rho_set, &
            kpp1_env%deriv_set, kpp1_env%spin_pot, kpp1_env%grad_pot,&
            kpp1_env%ndiag_term)
       kpp1_env%ref_count=1
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
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!*** **********************************************************************
  SUBROUTINE kpp1_calc_k_p_p1(kpp1_env, p_env, qs_env, k_p_p1, rho, rho1, error)
    

    TYPE(qs_kpp1_env_type), POINTER          :: kpp1_env
    TYPE(qs_p_env_type), POINTER             :: p_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: k_p_p1
    TYPE(qs_rho_type), POINTER               :: rho, rho1
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kpp1_calc_k_p_p1', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ispin, nspins, &
                                                ntot, output_unit, stat
    INTEGER, DIMENSION(2, 3)                 :: bo
    LOGICAL                                  :: failure, lsd, ionode, gapw
    REAL(KIND=dp)                            :: energy_hartree, &
                                                total_rho_gspace, &
                                                energy_hartree_1c,&
                                                energy_hartree_2c,&
                                                energy_hartree_3c
    TYPE(cell_type), POINTER                 :: cell
    TYPE(coeff_type)                         :: rho1_tot_gspace, &
                                                v_hartree_gspace, &
                                                v_hartree_rspace
    TYPE(coeff_type), DIMENSION(:), POINTER  :: rho1_r, &
                                                v_rspace_new
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_subsystem_p_type), &
      DIMENSION(:), POINTER                  :: atoms
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(pw_env_type), POINTER               :: pw_env
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho1_g_pw, rho1_r_pw, tau_pw, &
                                                v_xc
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pw_pools
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    TYPE(xc_control_type)                    :: xc_control
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set
    TYPE(xc_rho_cflags_type)                 :: needs
    TYPE(xc_rho_set_type), POINTER           :: rho1_set, rho_set

    CALL timeset(routineN,"I","",handle)
    failure=.FALSE.

    NULLIFY(dft_control, auxbas_pw_pool, pw_pools, pw_env, v_rspace_new, &
            rho1_r, rho1_g_pw, tau_pw, cell, v_xc, deriv_set, rho_set, rho1_set)
    logger => cp_error_get_logger(error)
    ionode = logger%para_env%mepos==logger%para_env%source

    gapw = qs_env%dft_control%qs_control%gapw

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

       nspins = dft_control%nspins
       CPPrecondition(SIZE(k_p_p1)>=nspins,cp_failure_level,routineP,error,failure)

       lsd = (dft_control%nspins == 2)
       xc_control = qs_env%dft_control%xc_control
       IF (dft_control%do_tddfpt_calculation) &
            xc_control = qs_env%dft_control%tddfpt_control%xc_control
    END IF

    IF (.NOT.failure) THEN
       kpp1_env%iter=kpp1_env%iter+1
    END IF

! gets the tmp grids
    IF (.NOT. failure) THEN
       CPPrecondition(ASSOCIATED(pw_env),cp_failure_level,routineP,error,failure)
       CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool,&
            pw_pools=pw_pools)
       ALLOCATE(v_rspace_new(nspins), stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    END IF
    IF (.NOT.failure) THEN
       CALL pw_pool_init_coeff(auxbas_pw_pool,v_hartree_gspace,&
            use_data = COMPLEXDATA1D,&
            in_space = RECIPROCALSPACE, error=error)
       CALL pw_pool_init_coeff(auxbas_pw_pool,v_hartree_rspace,&
            use_data = REALDATA3D,&
            in_space = REALSPACE, error=error)
    END IF

    IF (gapw) &
       CALL prepare_gapw_den(qs_env, .TRUE., p_env, error)


! *** calculate the hartree potential on the total density ***
    IF (.NOT. failure) THEN

       CALL pw_pool_init_coeff(auxbas_pw_pool, rho1_tot_gspace,&
            use_data = COMPLEXDATA1D,&
            in_space = RECIPROCALSPACE, error=error)

       CALL coeff_copy(rho1%rho_g(1),rho1_tot_gspace)
       DO ispin=2,nspins
          CALL coeff_sumup(rho1%rho_g(ispin), rho1_tot_gspace)
       END DO
       IF (qs_env%dft_control%qs_control%gapw) &
            CALL coeff_sumup(p_env%local_rho_set%rho0_mpole%rho0_s_gs, rho1_tot_gspace)

       CALL print_densities(kpp1_env, rho1, rho1_tot_gspace, error=error)

       IF (.NOT.(dft_control%do_tddfpt_calculation .AND. &
           dft_control%tddfpt_control%res_etype == "TRIPLET" )) THEN
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

       CALL qs_rho_get(rho1, rho_r=rho1_r, error=error)
       
       IF (dft_control%do_tddfpt_calculation .AND. nspins == 1) THEN
          ! In this case we need to do a unrestricted calculation to
          ! differenciate beween singlets and triplets
          lsd = .TRUE.
          ALLOCATE(rho1_r_pw(2))
          DO ispin=1, 2
             CALL pw_create(rho1_r_pw(ispin)%pw, rho1_r(1)%pw%pw_grid, &
                  rho1_r(1)%pw%in_use, rho1_r(1)%pw%in_space)
             CALL pw_transfer(rho1_r(1)%pw, rho1_r_pw(ispin)%pw)
             CALL dscal(PRODUCT(rho1_r_pw(ispin)%pw%pw_grid%npts_local), &
                  0.5_dp, rho1_r_pw(ispin)%pw%cr3d, 1)
          END DO
       ELSE
          ALLOCATE(rho1_r_pw(nspins))
          DO ispin=1, nspins
             rho1_r_pw(ispin)%pw => rho1_r(ispin)%pw
          END DO
       END IF
       NULLIFY(tau_pw)
       ALLOCATE(atoms(1))
       atoms(1)%subsys => qs_env%subsys

       !------!
       ! rho1 !
       !------!
       bo = rho1_r(1)%pw%pw_grid%bounds_local
       ! create the place where to store the argument for the functionals
       CALL xc_rho_set_create(rho1_set,bo,rho_cutoff=xc_control%density_cut,&
                              drho_cutoff=xc_control%gradient_cut,tau_cutoff=xc_control%tau_cut,&
                              error=error)
       ! assesses the argument needed by the functionals
       CALL xc_rho_cflags_setall(needs,.FALSE.,error=error)
       DO i=1,SIZE(xc_control%functionals)
          CALL xc_functional_get_info(xc_control%functionals(i),lsd=lsd,&
               needs=needs,&
               error=error)
       END DO
       ! needed basic components
       IF (lsd) THEN
          needs%rho_spin=.TRUE.
          needs%tau_spin=needs%tau_spin.or.needs%tau
       ELSE
          needs%rho=.TRUE.
       END IF
       IF (needs%drho.OR.needs%drho_spin.OR.needs%norm_drho.OR.&
            needs%norm_drho_spin) THEN
          IF (lsd) THEN
             needs%drho_spin=.TRUE.
          ELSE
             needs%drho=.TRUE.
          END IF
       END IF
       ! calculate the agruments needed by the functionals
       CALL xc_rho_set_update(rho1_set, rho1_r_pw, rho1_g_pw, tau_pw, needs,& 
                              dft_control, cell, auxbas_pw_pool, atoms, error)

       CALL  xc_calc_2nd_deriv_new(v_xc, kpp1_env%deriv_set, kpp1_env%rho_set, &
            rho1_set, auxbas_pw_pool, dft_control, error=error)
       do ispin=1,nspins
          call coeff_init(v_rspace_new(ispin), pw=v_xc(ispin)%pw, error=error)
       end do
       deallocate(v_xc,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       
       IF (gapw) CALL calculate_xc_2nd_deriv_atom(p_env, qs_env, error)

       CALL xc_rho_set_release(rho1_set)
       IF (dft_control%do_tddfpt_calculation .AND. nspins == 1) THEN
          DO ispin=1, 2
             CALL pw_release(rho1_r_pw(ispin)%pw)
          END DO
       END IF
       DEALLOCATE(rho1_r_pw, atoms,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       ntot = SIZE(v_hartree_rspace%pw%cr3d,1) * &
              SIZE(v_hartree_rspace%pw%cr3d,2) * &
              SIZE(v_hartree_rspace%pw%cr3d,3)
       CALL dscal(ntot, v_hartree_rspace%pw%pw_grid%dvol, &
                        v_hartree_rspace%pw%cr3d,1)

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
             IF (dft_control%tddfpt_control%res_etype == "SINGLET") THEN
                ! add 2 times the hartree term in the singlet case
                CALL daxpy(ntot, 2.0_dp, &
                     v_hartree_rspace%pw%cr3d, 1, v_rspace_new(ispin)%pw%cr3d, 1)
             END IF
             ! do not add the hartree term for triplets
          ELSE
             ! just add it if not tddfpt
             CALL daxpy(ntot, 1.0_dp, &
                  v_hartree_rspace%pw%cr3d, 1, v_rspace_new(ispin)%pw%cr3d, 1)
          END IF

          ! remove kpp1_env%v_ao and work directly on k_p_p1 ?
          CALL set_matrix(kpp1_env%v_ao(ispin)%matrix,0.0_dp)
          CALL integrate_v_rspace(v_rspace=v_rspace_new(ispin),&
               p=rho%rho_ao(ispin),&
               h=kpp1_env%v_ao(ispin),&
               qs_env=qs_env,&
               calculate_forces=.FALSE.,error=error)
          call transfer_matrix(kpp1_env%v_ao(ispin)%matrix,&
               k_p_p1(ispin)%matrix)
       END DO

       IF(gapw) THEN
          CALL Vh_1c_gg_integrals(qs_env,energy_hartree_1c, .TRUE., p_env)
          CALL integrate_vhg0_rspace(qs_env, v_hartree_rspace, &
                                     .FALSE., .TRUE., p_env, error=error)
          
          CALL integrate_vhg0_g0(qs_env,energy_hartree_2c,.FALSE.,.FALSE., &
                                 .TRUE., p_env)
          DO ispin = 1,nspins
             CALL integrate_a_vhg0_b(qs_env,k_p_p1(ispin),rho%rho_ao(ispin),&
                                     ispin,energy_hartree_3c,.FALSE.,.FALSE.,&
                                     .TRUE.,p_env,error)
          END DO
!         ***  Add single atom contributions to the KS matrix ***
             DO ispin = 1,nspins
                CALL update_ks_atom(qs_env,k_p_p1(ispin),rho%rho_ao(ispin),&
                     ispin,.FALSE.,.TRUE.,p_env,error)
             ENDDO
       END IF
       
       CALL pw_pool_give_back_coeff(auxbas_pw_pool,v_hartree_gspace,&
            error=error)
       CALL pw_pool_give_back_coeff(auxbas_pw_pool,v_hartree_rspace,&
            error=error)
       DO ispin=1,nspins
          CALL pw_pool_give_back_coeff(auxbas_pw_pool,v_rspace_new(ispin),&
               error=error)
       END DO
       DEALLOCATE(v_rspace_new, stat=stat)
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
       diff, error)
    TYPE(qs_kpp1_env_type), POINTER          :: kpp1_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: k_p_p1
    TYPE(qs_rho_type), POINTER               :: rho, rho1
    real(KIND=dp), intent(in), optional      :: diff
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kpp1_calc_k_p_p1_fdiff', &
      routineP = moduleN//':'//routineN

    integer                                  :: ispin, stat
    logical                                  :: failure
    real(KIND=dp)                            :: my_diff
    type(dft_control_type), pointer          :: dft_control
    type(qs_energy_type), pointer            :: qs_energy
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: ks_2, rho_ao, matrix_s

    nullify(ks_2,matrix_s,dft_control,qs_energy)
    my_diff=1.0e-6_dp
    if (present(diff)) my_diff=diff
    call allocate_qs_energy(qs_energy)

    call qs_rho_get(rho, rho_ao=rho_ao,error=error)
    call get_qs_env(qs_env, dft_control=dft_control,&
         matrix_s=matrix_s,error=error)

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
         qs_env=qs_env,ks_matrix=k_p_p1,rho=rho,energy=qs_energy,&
         calculate_forces=.false.,&
         just_energy=.false.,error=error)

    allocate(ks_2(dft_control%nspins),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    do ispin=1,dft_control%nspins
       nullify(ks_2(ispin)%matrix)
       call replicate_matrix_structure(matrix_s(1)%matrix,ks_2(ispin)%matrix,&
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
         qs_env=qs_env,ks_matrix=ks_2,rho=rho,energy=qs_energy,&
         calculate_forces=.false.,&
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
!!*** **********************************************************************
  SUBROUTINE kpp1_check_i_alloc(kpp1_env, qs_env, error)


    TYPE(qs_kpp1_env_type), POINTER          :: kpp1_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kpp1_check_i_alloc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ispin, nspins, stat
    LOGICAL                                  :: failure
    type(cell_type), pointer                 :: cell
    type(coeff_type), dimension(:), pointer  :: rho_r
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(pw_env_type), POINTER               :: pw_env
    type(pw_p_type), dimension(:), pointer   :: my_rho_r
    TYPE(pw_pool_type), POINTER              :: auxbas_pw_pool
    type(qs_rho_type), pointer               :: rho
    TYPE(real_matrix_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(xc_derivative_set_type), POINTER    :: deriv_set

! ------------------------------------------------------------------

    failure=.FALSE.
    NULLIFY(pw_env,auxbas_pw_pool,matrix_s,rho,rho_r)

    CPPrecondition(ASSOCIATED(kpp1_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(kpp1_env%ref_count>0,cp_failure_level,routineP,error,failure)

    IF (.NOT. failure) THEN


       CALL get_qs_env(qs_env, pw_env=pw_env, dft_control=dft_control,&
                       matrix_s=matrix_s, cell=cell, error=error, rho=rho)

       nspins=dft_control%nspins

       CALL pw_env_get(pw_env, auxbas_pw_pool = auxbas_pw_pool, error=error)

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
                CALL replicate_matrix_structure(matrix_s(1)%matrix,&
                     kpp1_env%v_ao(ispin)%matrix,&
                     "kpp1%v_ao-"//ADJUSTL(cp_to_string(ispin)))
             END DO
          END IF
       END IF
       
       if (.not.associated(kpp1_env%deriv_set)) then

          call qs_rho_get(rho, rho_r=rho_r, error=error)

          NULLIFY(deriv_set)
          allocate(my_rho_r(size(rho_r)),stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
          do ispin=1,size(rho_r)
             my_rho_r(ispin)%pw => rho_r(ispin)%pw
          end do
          call xc_prep_2nd_deriv_new(kpp1_env%deriv_set, kpp1_env%rho_set, &
               my_rho_r, auxbas_pw_pool, &
               dft_control, cell=cell, error=error)
          deallocate(my_rho_r,stat=stat)
          CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
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
    TYPE(qs_kpp1_env_type), POINTER          :: kpp1_env
    LOGICAL, INTENT(in), OPTIONAL            :: s_struct_changed, &
                                                grid_changed, psi0_changed
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'kpp1_did_change', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure, my_psi0_changed, &
                                                my_s_struct_changed

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
          IF (ASSOCIATED(kpp1_env%v_ao)) THEN
             CALL deallocate_matrix_set(kpp1_env%v_ao)
          END IF
       END IF
       if (my_s_struct_changed.or.my_psi0_changed) then
          if (associated(kpp1_env%drho_r)) then
             deallocate(kpp1_env%drho_r, stat=stat)
             CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
          end if
          if (associated(kpp1_env%deriv_set)) then
             call xc_dset_release(kpp1_env%deriv_set, error=error)
             nullify(kpp1_env%deriv_set)
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

  SUBROUTINE print_densities(kpp1_env, rho1, rho1_tot_gspace, error)
    
    ! arguments
    TYPE(qs_kpp1_env_type), POINTER          :: kpp1_env
    TYPE(qs_rho_type), POINTER               :: rho1
    TYPE(coeff_type), INTENT(IN)             :: rho1_tot_gspace
    TYPE(cp_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error

    ! locals
    INTEGER                                  :: output_unit
    LOGICAL                                  :: ionode
    REAL(KIND=dp)                            :: total_rho_gspace
    TYPE(cp_logger_type), POINTER            :: logger

    ! parameters
    CHARACTER(len=*), PARAMETER :: routineN = 'print_densities', &
                                   routineP = moduleN//':'//routineN

    NULLIFY(logger)
    logger => cp_error_get_logger(error)
    ionode = (logger%para_env%mepos==logger%para_env%source)

    IF (.NOT.ionode) RETURN

    IF (logger%print_keys%total_densities) THEN
       total_rho_gspace = calculate_total_rho(rho1_tot_gspace)
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

  END SUBROUTINE print_densities

END MODULE qs_kpp1_env_methods
