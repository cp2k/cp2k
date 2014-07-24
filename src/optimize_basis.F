!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!
MODULE optimize_basis
  USE admm_methods,                    ONLY: fit_mo_coeffs_diag
  USE cp_fm_types,                     ONLY: cp_fm_p_type,&
                                             cp_fm_release
  USE cp_para_env,                     ONLY: cp_para_env_create,&
                                             cp_para_env_release
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE f77_blas
  USE f77_interface,                   ONLY: create_force_env,&
                                             destroy_force_env,&
                                             f_env_add_defaults,&
                                             f_env_get_from_id,&
                                             f_env_rm_defaults,&
                                             f_env_type
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE input_cp2k_read,                 ONLY: empty_initial_variables,&
                                             read_input
  USE input_section_types,             ONLY: section_type,&
                                             section_vals_release,&
                                             section_vals_type
  USE kinds,                           ONLY: default_path_length,&
                                             dp
  USE machine,                         ONLY: m_chdir,&
                                             m_getcwd,&
                                             m_walltime
  USE message_passing,                 ONLY: mp_bcast,&
                                             mp_comm_free,&
                                             mp_comm_split,&
                                             mp_sum,&
                                             mp_sync
  USE optbas_fenv_manipulation,        ONLY: allocate_mo_sets,&
                                             calculate_ks_matrix,&
                                             calculate_overlap,&
                                             calculate_overlap_inverse,&
                                             create_opt_admm_env,&
                                             modify_input_settings,&
                                             update_basis_set
  USE optbas_opt_utils,                ONLY: evaluate_energy,&
                                             evaluate_fval
  USE optimize_basis_types,            ONLY: basis_optimization_type,&
                                             deallocate_basis_optimization_type,&
                                             subset_type
  USE optimize_basis_utils,            ONLY: get_set_and_basis_id,&
                                             optimize_basis_init_read_input,&
                                             update_derived_basis_sets
  USE powell,                          ONLY: powell_optimize
  USE qs_environment_types,            ONLY: qs_environment_type,&
                                             set_qs_env
  USE qs_neighbor_lists,               ONLY: build_qs_neighbor_lists
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: run_optimize_basis

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'optimize_basis'

CONTAINS

! *****************************************************************************
!> \brief main entry point for methods aimed at optimizing basis sets
!> \param input_declaration ...
!> \param root_section ...
!> \param para_env ...
!> \param error ...
!> \author Florian Schiffmann
! *****************************************************************************
  SUBROUTINE run_optimize_basis(input_declaration,root_section,para_env,error)
    TYPE(section_type), POINTER              :: input_declaration
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'run_optimize_basis', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(basis_optimization_type)            :: opt_bas
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)
    logger => cp_get_default_logger()
    CALL cp_error_init(error,stop_level=cp_fatal_level, logger=logger)

    CALL optimize_basis_init_read_input(opt_bas,root_section,para_env,error)

    CALL driver_para_opt_basis(opt_bas,input_declaration,para_env,error)

    CALL deallocate_basis_optimization_type(opt_bas)
    CALL cp_error_dealloc_ref(error)
    CALL timestop(handle)


  END SUBROUTINE run_optimize_basis

! *****************************************************************************
!> \brief driver routine for the parallel part of the method
!> \param opt_bas ...
!> \param input_declaration ...
!> \param para_env ...
!> \param error ...
!> \author Florian Schiffmann
! *****************************************************************************

  SUBROUTINE driver_para_opt_basis(opt_bas,input_declaration,para_env,error)
    TYPE(basis_optimization_type)            :: opt_bas
    TYPE(section_type), POINTER              :: input_declaration
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'driver_para_opt_basis', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, n_groups_created, &
                                                opt_group
    INTEGER, DIMENSION(:), POINTER           :: group_distribution_p
    INTEGER, &
      DIMENSION(0:para_env%num_pe-1), TARGET :: group_distribution

    CALL timeset(routineN,handle)
    group_distribution_p=>group_distribution
    CALL mp_comm_split(para_env%group, opt_group, n_groups_created, group_distribution_p,&
         n_subgroups=SIZE(opt_bas%group_partition),group_partition=opt_bas%group_partition)
    opt_bas%opt_id = group_distribution(para_env%mepos)+1

    CALL driver_optimization_para_low(opt_bas,input_declaration,para_env,opt_group,error)

    CALL mp_comm_free(opt_group)
    CALL timestop(handle)

  END SUBROUTINE driver_para_opt_basis

! *****************************************************************************
!> \brief low level optimization routine includes initialization of the subsytems
!>        powell optimizer and deallocation of the various force envs
!> \param opt_bas ...
!> \param input_declaration ...
!> \param para_env_top ...
!> \param mpi_comm_opt ...
!> \param top_error ...
!> \author Florian Schiffmann
! *****************************************************************************

  SUBROUTINE driver_optimization_para_low(opt_bas,input_declaration,para_env_top,mpi_comm_opt,top_error)
    TYPE(basis_optimization_type)            :: opt_bas
    TYPE(section_type), POINTER              :: input_declaration
    TYPE(cp_para_env_type), POINTER          :: para_env_top
    INTEGER                                  :: mpi_comm_opt
    TYPE(cp_error_type), INTENT(INOUT)       :: top_error

    CHARACTER(len=*), PARAMETER :: routineN = 'driver_optimization_para_low', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, icalc, iopt, mp_id, &
                                                stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: f_env_id
    LOGICAL                                  :: failure, write_basis
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: tot_time
    TYPE(cp_error_type)                      :: error
    TYPE(cp_fm_p_type), ALLOCATABLE, &
      DIMENSION(:)                           :: matrix_S_inv
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(f_env_type), POINTER                :: f_env

    NULLIFY(f_env)

    CALL timeset(routineN,handle)
    failure=.FALSE.

    ! ======  initialize the f_env and precompute some matrices ===== 
    mp_id=opt_bas%opt_id
    NULLIFY(para_env, f_env)
    ALLOCATE(f_env_id(SIZE(opt_bas%comp_group(mp_id)%member_list)))
    ALLOCATE(tot_time(opt_bas%ncombinations*opt_bas%ntraining_sets))
    ALLOCATE(matrix_s_inv(SIZE(opt_bas%comp_group(mp_id)%member_list)))

    CALL cp_error_init(error, stop_level=cp_fatal_level)
    CALL cp_para_env_create(para_env, group=mpi_comm_opt,&
         owns_group=.FALSE.,error=error)


    CALL init_training_force_envs(opt_bas,f_env_id,input_declaration,matrix_s_inv,para_env,mpi_comm_opt,error)
 
    CALL init_free_vars(opt_bas,error)
    tot_time=0.0_dp
    
    ! ======= The real optimization loop  =======
    DO iopt=0,opt_bas%powell_param%maxfun
      CALL compute_residuum_vectors(opt_bas,f_env_id,matrix_S_inv,tot_time,&
                                    para_env_top,para_env,iopt,error,top_error)
      IF(para_env_top%ionode)&
      CALL powell_optimize (opt_bas%powell_param%nvar, opt_bas%x_opt , opt_bas%powell_param)
      IF(.NOT.para_env_top%ionode)opt_bas%x_opt=0.0_dp
      CALL mp_bcast(opt_bas%powell_param%state,para_env_top%source,para_env_top%group)
      CALL mp_sum(opt_bas%x_opt,para_env_top%group)
      CALL update_free_vars(opt_bas,error)
      write_basis=MOD(iopt,opt_bas%write_frequency)==0
      CALL update_derived_basis_sets(opt_bas,write_basis,opt_bas%output_basis_file,&
                                       para_env_top,top_error)
      IF(opt_bas%powell_param%state==-1)EXIT
    END DO
 
    ! ======= Update the basis set and print the final basis  =======
    IF(para_env_top%ionode)THEN
       opt_bas%powell_param%state = 8
       CALL powell_optimize (opt_bas%powell_param%nvar, opt_bas%x_opt , opt_bas%powell_param)
    END IF
 
    IF(.NOT.para_env_top%ionode)opt_bas%x_opt=0.0_dp
    CALL mp_sum(opt_bas%x_opt,para_env_top%group)
    CALL update_free_vars(opt_bas,error)
    CALL update_derived_basis_sets(opt_bas,.TRUE.,opt_bas%output_basis_file,&
                                         para_env_top,top_error)
 
     ! ======  get rid of the f_env again =====          
    
     DO icalc=SIZE(opt_bas%comp_group(mp_id)%member_list),1,-1
        CALL f_env_get_from_id(f_env_id(icalc),f_env,error)
        CALL destroy_force_env(f_env_id(icalc), stat)
        CALL cp_fm_release(matrix_s_inv(icalc)%matrix,error)
     END DO
     DEALLOCATE(f_env_id); DEALLOCATE(tot_time); DEALLOCATE(matrix_S_inv) 
     CALL cp_para_env_release(para_env,error=error)
     CALL cp_error_dealloc_ref(error)   
     CALL timestop(handle)
 
  END SUBROUTINE driver_optimization_para_low

! *****************************************************************************
!> \brief compute all ingredients for powell optimizer. Rho_diff,
!>        condition number, energy,... for all ttraining sets in 
!>        the computational group
!> \param opt_bas ...
!> \param f_env_id ...
!> \param matrix_S_inv ...
!> \param tot_time ...
!> \param para_env_top ...
!> \param para_env ...
!> \param iopt ...
!> \param error ...
!> \param top_error ...
!>\author Florian Schiffmann
! *****************************************************************************

  SUBROUTINE compute_residuum_vectors(opt_bas,f_env_id,matrix_S_inv,tot_time,&
                                      para_env_top,para_env,iopt,error,top_error)
    TYPE(basis_optimization_type)            :: opt_bas
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: f_env_id
    TYPE(cp_fm_p_type), DIMENSION(:)         :: matrix_S_inv
    REAL(KIND=dp), DIMENSION(:)              :: tot_time
    TYPE(cp_para_env_type), POINTER          :: para_env_top, para_env
    INTEGER                                  :: iopt
    TYPE(cp_error_type)                      :: error, top_error

    CHARACTER(len=*), PARAMETER :: routineN = 'compute_residuum_vectors', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: bas_id, handle, icalc, icomb, &
                                                mp_id, my_id, ncalc, set_id
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: cond_vec, energy, f_vec, &
                                                my_time, start_time
    TYPE(f_env_type), POINTER                :: f_env
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(qs_environment_type), POINTER       :: qs_env

    CALL timeset(routineN,handle)
 
    ncalc=opt_bas%ncombinations*opt_bas%ntraining_sets
    ALLOCATE(f_vec(ncalc)); ALLOCATE(my_time(ncalc)); ALLOCATE(cond_vec(ncalc)); ALLOCATE(energy(ncalc))
    f_vec=0.0_dp ; cond_vec=0.0_dp; my_time=0.0_dp; energy=0.0_dp
    mp_id=opt_bas%opt_id
    ALLOCATE(start_time(SIZE(opt_bas%comp_group(mp_id)%member_list)))

    DO icalc=1,SIZE(opt_bas%comp_group(mp_id)%member_list)
        my_id=opt_bas%comp_group(mp_id)%member_list(icalc)+1
        ! setup timings
        start_time(icalc)=m_walltime()
 
        CALL get_set_and_basis_id(opt_bas%comp_group(mp_id)%member_list(icalc),opt_bas,set_id,bas_id)
        CALL f_env_get_from_id(f_env_id(icalc),f_env,error)
        force_env => f_env%force_env
        CALL force_env_get(force_env,qs_env=qs_env,error=error)
 
        CALL update_basis_set(opt_bas,set_id,bas_id,qs_env,error)
        CALL build_qs_neighbor_lists(qs_env,para_env,molecular=.FALSE.,force_env_section=qs_env%input,error=error)
        CALL calculate_overlap(qs_env,"S_AB_AUX",error)
        CALL fit_mo_coeffs_diag(qs_env, qs_env%admm_env, qs_env%para_env, qs_env%matrix_s_aux_fit,&
                             qs_env%matrix_s_aux_fit_vs_orb, qs_env%mos, qs_env%mos_aux_fit, .TRUE. , error)
        CALL evaluate_fval(qs_env%mos,qs_env%mos_aux_fit,qs_env%matrix_s_aux_fit_vs_orb(1)%matrix,&
                         qs_env%matrix_s_aux_fit(1)%matrix,qs_env%admm_env,f_vec(my_id),cond_vec(my_id),error)
        CALL evaluate_energy(qs_env%mos_aux_fit,qs_env%matrix_ks,matrix_s_inv(icalc)%matrix,qs_env%admm_env%Q,&
                             qs_env%admm_env%work_aux_aux,energy(my_id),error)
 
        my_time(my_id)=m_walltime()-start_time(icalc)
        IF(.NOT.para_env%ionode)THEN
            f_vec=0.0_dp; cond_vec=0.0_dp; my_time=0.0_dp; energy=0.0_dp
        END IF
     END DO

     DEALLOCATE(start_time)
 
     CALL mp_sum(f_vec,para_env_top%group); CALL mp_sum(cond_vec,para_env_top%group); CALL mp_sum(my_time,para_env_top%group)
     opt_bas%powell_param%f=0.0_dp
     DO icalc=1,SIZE(f_vec)
        icomb=MOD(icalc-1,opt_bas%ncombinations)
        opt_bas%powell_param%f=opt_bas%powell_param%f+&
              (f_vec(icalc)+energy(icalc))*opt_bas%fval_weight(icomb)
        IF(opt_bas%use_condition_number)&
          opt_bas%powell_param%f=opt_bas%powell_param%f+&
              LOG(cond_vec(icalc))*opt_bas%condition_weight(icomb)
     END DO

     CALL mp_sync(para_env_top%group)

     ! output info if required 
     CALL output_opt_info(f_vec,cond_vec,my_time,tot_time,opt_bas,iopt,para_env_top,top_error)
     DEALLOCATE(f_vec); DEALLOCATE(my_time); DEALLOCATE(cond_vec); DEALLOCATE(energy)

     CALL timestop(handle)

  END SUBROUTINE compute_residuum_vectors

! *****************************************************************************
!> \brief create the force_envs for every input in the computational group
!> \param opt_bas ...
!> \param f_env_id ...
!> \param input_declaration ...
!> \param matrix_s_inv ...
!> \param para_env ...
!> \param mpi_comm_opt ...
!> \param error ...
!>\author Florian Schiffmann
! *****************************************************************************

  SUBROUTINE init_training_force_envs(opt_bas,f_env_id,input_declaration,matrix_s_inv,para_env,mpi_comm_opt,error)

    TYPE(basis_optimization_type)            :: opt_bas
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: f_env_id
    TYPE(section_type), POINTER              :: input_declaration
    TYPE(cp_fm_p_type), DIMENSION(:)         :: matrix_S_inv
    TYPE(cp_para_env_type), POINTER          :: para_env
    INTEGER                                  :: mpi_comm_opt
    TYPE(cp_error_type)                      :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'init_training_force_envs', &
      routineP = moduleN//':'//routineN

    CHARACTER(len=default_path_length)       :: main_dir
    INTEGER                                  :: bas_id, handle, icalc, ierr, &
                                                mp_id, set_id, stat
    LOGICAL                                  :: failure
    TYPE(cp_error_type), ALLOCATABLE, &
      DIMENSION(:)                           :: suberror
    TYPE(f_env_type), POINTER                :: f_env
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(section_vals_type), POINTER         :: input_file

    CALL timeset(routineN,handle)

    mp_id=opt_bas%opt_id
    CALL m_getcwd(main_dir)
    ALLOCATE(suberror(SIZE(opt_bas%comp_group(mp_id)%member_list)))

    ! ======= Create f_env for all calculations in MPI group =======
    DO icalc=1,SIZE(opt_bas%comp_group(mp_id)%member_list)
       NULLIFY(input_file)
       ! parse the input of the training sets
       CALL get_set_and_basis_id(opt_bas%comp_group(mp_id)%member_list(icalc),opt_bas,set_id,bas_id)
       CALL m_chdir(TRIM(opt_bas%training_dir(set_id)),ierr)
       CPPostcondition(ierr==0,cp_failure_level,routineP,error,failure)
       input_file => read_input(input_declaration,&
                                opt_bas%training_input(set_id),&
                                initial_variables=empty_initial_variables, &
                                para_env=para_env,error=error)
       CALL cp_error_check(error,failure)

       CALL modify_input_settings(opt_bas,bas_id,input_file,error)
       CALL create_force_env(f_env_id(icalc), &
                             input_declaration=input_declaration,&
                             input_path=opt_bas%training_input(set_id),&
                             input=input_file, &
                             output_path="scrap_information", &
                             mpi_comm=mpi_comm_opt,&
                             ierr=stat)

       ! some weirdness with the default stacks defaults have to addded to get the 
       ! correct default program name this causes trouble with the timer stack if kept
       CALL f_env_add_defaults(f_env_id(icalc), f_env, suberror(icalc), failure)
       force_env => f_env%force_env
       CALL force_env_get(force_env,qs_env=qs_env,error=error)
       CALL set_qs_env(qs_env=qs_env,atprop=force_env%atprop_env,error=error)
       CALL allocate_mo_sets(qs_env,error)
       CALL f_env_rm_defaults(f_env, suberror(icalc),stat)

       CALL build_qs_neighbor_lists(qs_env,para_env,molecular=.FALSE.,&
                                    force_env_section=qs_env%input,error=error)
       CALL calculate_overlap(qs_env,"S_AB",error)
       CALL calculate_overlap_inverse(qs_env%matrix_s(1)%matrix,matrix_s_inv(icalc)%matrix,&
                                      para_env,qs_env%blacs_env,error)
       CALL calculate_ks_matrix(qs_env,suberror(icalc))

       CALL create_opt_admm_env(qs_env,error)
       CALL section_vals_release(input_file,error=error)
       CALL m_chdir(TRIM(ADJUSTL(main_dir)),ierr)
    END DO

    DEALLOCATE(suberror)

    CALL timestop(handle)

  END SUBROUTINE init_training_force_envs

! *****************************************************************************
!> \brief variable update from the powell vector for all sets
!> \param opt_bas ...
!> \param error ...
!> \author Florian Schiffmann
! *****************************************************************************

  SUBROUTINE update_free_vars(opt_bas,error)
    TYPE(basis_optimization_type)            :: opt_bas
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'update_free_vars', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ikind, iset, ix

    CALL timeset(routineN,handle)
    ix=0
    DO ikind=1,opt_bas%nkind
       DO iset=1,opt_bas%kind_basis(ikind)%flex_basis(0)%nsets
          CALL update_subset_freevars(opt_bas%kind_basis(ikind)%flex_basis(0)%subset(iset),ix,opt_bas%x_opt,error)
       END DO
    END DO
    CALL timestop(handle)

  END SUBROUTINE update_free_vars

! *****************************************************************************
!> \brief low level update for the basis sets. Exponents are transformed according to constraint
!> \param subset ...
!> \param ix ...
!> \param x ...
!> \param error ...
!> \author Florian Schiffmann
! *****************************************************************************

  SUBROUTINE update_subset_freevars(subset,ix,x,error)
    TYPE(subset_type)                        :: subset
    INTEGER                                  :: ix
    REAL(KIND=dp), DIMENSION(:)              :: x
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'update_subset_freevars', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, icon1, icon2, icont, &
                                                iexp, il, istart
    REAL(KIND=dp)                            :: fermi_f, gs_scale

    CALL timeset(routineN,handle)
       DO iexp=1,subset%nexp
          IF(subset%opt_exps(iexp))THEN
             ix=ix+1
             subset%exps(iexp)=ABS(x(ix))
             IF(subset%exp_has_const(iexp))THEN
                 !use a fermi function to keep expoenents in a given range around their initial value
                 fermi_f=1.0_dp/ (EXP( (x(ix)-1.0_dp)/0.5_dp ) +1.0_dp)
                 subset%exps(iexp)=(2.0_dp*fermi_f-1.0_dp)*subset%exp_const(iexp)%var_fac*subset%exp_const(iexp)%init+&
                                   subset%exp_const(iexp)%init
             ELSE
                 
             END IF
          END IF
          DO icont=1,subset%ncon_tot
             IF(subset%opt_coeff(iexp,icont))THEN
                ix=ix+1
                subset%coeff(iexp,icont)=x(ix)
             END IF
          END DO
       END DO

      ! orthonormalize contraction coefficients using gram schmidt
      istart=1
      DO il=1,subset%nl
         DO icon1=istart,istart+subset%l(il)-2
            DO icon2=icon1+1,istart+subset%l(il)-1
               gs_scale=DOT_PRODUCT(subset%coeff(:,icon2),subset%coeff(:,icon1))/&
                        DOT_PRODUCT(subset%coeff(:,icon1),subset%coeff(:,icon1))
               subset%coeff(:,icon2)=subset%coeff(:,icon2)-gs_scale*subset%coeff(:,icon1)
            END DO
         END DO
         istart=istart+subset%l(il)
      END DO

       DO icon1=1,subset%ncon_tot
          subset%coeff(:,icon1)=subset%coeff(:,icon1)/&
                SQRT(DOT_PRODUCT(subset%coeff(:,icon1),subset%coeff(:,icon1)))
       END DO
    CALL timestop(handle)

  END SUBROUTINE update_subset_freevars

! *****************************************************************************
!> \brief variable initialization for the powell vector for all sets
!> \param opt_bas ...
!> \param error ...
!> \author Florian Schiffmann
! *****************************************************************************

  SUBROUTINE init_free_vars(opt_bas,error)
    TYPE(basis_optimization_type)            :: opt_bas
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'init_free_vars', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ikind, iset, ix

    CALL timeset(routineN,handle)
    ix=0
    DO ikind=1,opt_bas%nkind
       DO iset=1,opt_bas%kind_basis(ikind)%flex_basis(0)%nsets
          CALL init_subset_freevars(opt_bas%kind_basis(ikind)%flex_basis(0)%subset(iset),ix,opt_bas%x_opt,error)
       END DO
    END DO
    CALL timestop(handle)

    END SUBROUTINE init_free_vars

! *****************************************************************************
!> \brief variable initialization for the powell vector from low level informations
!>        constraint exponents will be mapped on a fermi function 
!> \param subset ...
!> \param ix ...
!> \param x ...
!> \param error ...
!> \author Florian Schiffmann
! *****************************************************************************

    SUBROUTINE init_subset_freevars(subset,ix,x,error)
    TYPE(subset_type)                        :: subset
    INTEGER                                  :: ix
    REAL(KIND=dp), DIMENSION(:)              :: x
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'init_subset_freevars', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, icont, iexp
    REAL(KIND=dp)                            :: fract

    CALL timeset(routineN,handle)

       DO iexp=1,subset%nexp
          IF(subset%opt_exps(iexp))THEN
             ix=ix+1
             x(ix)=subset%exps(iexp)
             IF(subset%exp_has_const(iexp))THEN
                 IF(subset%exp_const(iexp)%const_type==0)THEN
                    fract=1.0_dp+(subset%exps(iexp)-subset%exp_const(iexp)%init)/&
                          (subset%exp_const(iexp)%init*subset%exp_const(iexp)%var_fac)
                    x(ix)=0.5_dp*LOG((2.0_dp/fract-1.0_dp))+1.0_dp
                 END IF
                 IF(subset%exp_const(iexp)%const_type==1)THEN
                    x(ix)=1.0_dp
                 END IF
             END IF
          END IF
          DO icont=1,subset%ncon_tot
             IF(subset%opt_coeff(iexp,icont))THEN
                ix=ix+1
                x(ix)=subset%coeff(iexp,icont)
             END IF
          END DO
       END DO
    CALL timestop(handle)

    END SUBROUTINE init_subset_freevars

! *****************************************************************************
!> \brief commuticates all info to the master and assembles the output
!> \param f_vec ...
!> \param cond_vec ...
!> \param my_time ...
!> \param tot_time ...
!> \param opt_bas ...
!> \param iopt ...
!> \param para_env_top ...
!> \param error ...
!> \author Florian Schiffmann
! *****************************************************************************

   SUBROUTINE  output_opt_info(f_vec,cond_vec,my_time,tot_time,opt_bas,iopt,para_env_top,error)
    REAL(KIND=dp), DIMENSION(:)              :: f_vec, cond_vec, my_time, &
                                                tot_time
    TYPE(basis_optimization_type)            :: opt_bas
    INTEGER                                  :: iopt
    TYPE(cp_para_env_type), POINTER          :: para_env_top
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'output_opt_info', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ibasis, icalc, iset, &
                                                unit_nr
    TYPE(cp_logger_type), POINTER            :: logger

    CALL timeset(routineN,handle)
    logger => cp_error_get_logger(error)
    
    tot_time=tot_time+my_time

    unit_nr=-1
    IF(para_env_top%ionode.AND.(MOD(iopt,opt_bas%write_frequency)==0.OR.iopt==opt_bas%powell_param%maxfun))&
       unit_nr=cp_logger_get_default_unit_nr(logger)
     
    IF(unit_nr.gt.0)THEN
       WRITE(unit_nr,'(1X,A,I8)')"BASOPT| Information at iteration number:",iopt
       WRITE(unit_nr,'(1X,A)')"BASOPT| Training set | Combination | Rho difference | Condition num. | Time"
       WRITE(unit_nr,'(1X,A)')"BASOPT| -----------------------------------------------------------------------"
       icalc=0
       DO iset=1,opt_bas%ntraining_sets
          DO ibasis=1,opt_bas%ncombinations
             icalc=icalc+1
             WRITE(unit_nr,'(1X,A,2(5X,I3,5X,A),2(1X,E14.8,1X,A),1X,F8.3)')&
                'BASOPT| ',iset,"|",ibasis,"|",f_vec(icalc),"|",cond_vec(icalc),"|",tot_time(icalc)
          END DO
       END DO
       WRITE(unit_nr,'(1X,A)')"BASOPT| -----------------------------------------------------------------------"
       WRITE(unit_nr,'(1X,A,E14.8)')"BASOPT| Total residuum value: ", opt_bas%powell_param%f
       WRITE(unit_nr,'(A)')""
    END IF
    CALL timestop(handle)
   END SUBROUTINE output_opt_info

END MODULE optimize_basis
     
