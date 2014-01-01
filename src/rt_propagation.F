!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Routines for the real time propagation.
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

MODULE rt_propagation
  USE cell_types,                      ONLY: cell_type
  USE cp_cfm_basic_linalg,             ONLY: cp_cfm_column_scale,&
                                             cp_cfm_gemm
  USE cp_cfm_diag,                     ONLY: cp_cfm_heevd
  USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                             cp_cfm_release,&
                                             cp_cfm_to_cfm,&
                                             cp_cfm_type
  USE cp_control_types,                ONLY: rtp_control_type
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_copy,&
                                             cp_dbcsr_init
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                             cp_dbcsr_deallocate_matrix,&
                                             cp_dbcsr_sm_fm_multiply
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_p_type,&
                                             cp_dbcsr_type
  USE cp_external_control,             ONLY: external_control
  USE cp_fm_basic_linalg,              ONLY: cp_fm_gemm,&
                                             cp_fm_scale_and_add,&
                                             cp_fm_upper_to_full
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                             cp_fm_cholesky_invert,&
                                             cp_fm_cholesky_reduce,&
                                             cp_fm_cholesky_restore
  USE cp_fm_diag,                      ONLY: cp_fm_syevd
  USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type
  USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                             cp_fm_struct_release,&
                                             cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_p_type,&
                                             cp_fm_release,&
                                             cp_fm_set_all,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_output_handling,              ONLY: cp_add_iter_level,&
                                             cp_iterate,&
                                             cp_rm_iter_level
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE efield_utils,                    ONLY: calculate_ecore_efield
  USE force_env_methods,               ONLY: force_env_calc_energy_force
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE input_constants,                 ONLY: extr_s_ks,&
                                             real_time_propagation,&
                                             use_restart_wfn,&
                                             use_rt_restart,&
                                             use_scf_wfn
  USE input_cp2k_restarts,             ONLY: write_restart
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_val_set
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_walltime
  USE mathconstants,                   ONLY: twopi
  USE md_environment_types,            ONLY: md_environment_type
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_energy_utils,                 ONLY: qs_energies_init
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_external_potential,           ONLY: external_c_potential,&
                                             external_e_potential
  USE qs_ks_methods,                   ONLY: qs_ks_did_change,&
                                             qs_ks_update_qs_env
  USE qs_matrix_pools,                 ONLY: mpools_get
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             init_mo_set,&
                                             mo_set_p_type
  USE qs_moments,                      ONLY: build_berry_moment_matrix
  USE rt_propagation_methods,          ONLY: calc_update_rho,&
                                             propagation_step,&
                                             put_data_to_history
  USE rt_propagation_output,           ONLY: rt_prop_output
  USE rt_propagation_types,            ONLY: get_rtp,&
                                             rt_prop_create,&
                                             rt_prop_type,&
                                             rtp_history_create
  USE rt_propagation_utils,            ONLY: get_restart_wfn
  USE rt_propagator_init,              ONLY: init_propagators
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rt_propagation'

  PUBLIC :: rt_prop_setup,&
       init_ehrenfest_md,&
       run_propagation, &
       rt_write_input_restart


CONTAINS

! *****************************************************************************
!> \brief creates rtp_type, gets the initial state, either by reading MO's
!>        from file or calling SCF run
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE rt_prop_setup(force_env,error)
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rt_prop_setup', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    INTEGER                                  :: aspc_order, homo, ispin, stat
    LOGICAL                                  :: extrapolate_H, failure
    TYPE(cp_fm_pool_p_type), DIMENSION(:), &
      POINTER                                :: ao_mo_fm_pools, &
                                                ao_mo_fm_pools_aux_fit
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_coeff_aux_fit
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(rt_prop_type), POINTER              :: rtp
    TYPE(rtp_control_type), POINTER          :: rtp_control
    TYPE(section_vals_type), POINTER         :: md_section, motion_section

    failure=.FALSE.
    NULLIFY(qs_env,rtp_control)

    CALL force_env_get(force_env=force_env,qs_env=qs_env,globenv=globenv,error=error)
    CALL set_qs_env(qs_env=qs_env,atprop=force_env%atprop_env,error=error)
    rtp_control=>qs_env%dft_control%rtp_control

    SELECT CASE(rtp_control%initial_wfn)
    CASE(use_scf_wfn)
       qs_env%sim_time=0.0_dp
       qs_env%sim_step=0
       CALL force_env_calc_energy_force(force_env,calc_force=.FALSE.,&
            consistent_energies=.TRUE.,error=error)
       qs_env%run_rtp=.TRUE.
       ALLOCATE(qs_env%rtp, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL rt_prop_create(qs_env%rtp,qs_env%mos,qs_env%mpools,qs_env%dft_control,error)

    CASE(use_restart_wfn,use_rt_restart)
       CALL qs_energies_init(qs_env, .FALSE. , error)
       DO ispin=1,SIZE(qs_env%mos)
          CALL get_mo_set(qs_env%mos(ispin)%mo_set,mo_coeff=mo_coeff,homo=homo)
          IF (.NOT.ASSOCIATED(mo_coeff)) THEN
             CALL init_mo_set(qs_env%mos(ispin)%mo_set,&
                  qs_env%mpools%ao_mo_fm_pools(ispin)%pool,&
                  name="qs_env"//TRIM(ADJUSTL(cp_to_string(qs_env%id_nr)))//&
                  "%mo"//TRIM(ADJUSTL(cp_to_string(ispin))),&
                  error=error)
          END IF
       END DO

       IF(qs_env%dft_control%do_admm) THEN
          CALL mpools_get(qs_env%mpools_aux_fit, ao_mo_fm_pools=ao_mo_fm_pools_aux_fit,&
               error=error)
          CPPrecondition(ASSOCIATED(qs_env%mos_aux_fit),cp_failure_level,routineP,error,failure)
          IF (.NOT.failure) THEN
             DO ispin=1,SIZE(qs_env%mos_aux_fit)
                CALL get_mo_set(qs_env%mos_aux_fit(ispin)%mo_set,mo_coeff=mo_coeff_aux_fit,homo=homo)
                IF (.NOT.ASSOCIATED(mo_coeff_aux_fit)) THEN
                   CALL init_mo_set(qs_env%mos_aux_fit(ispin)%mo_set,&
                        ao_mo_fm_pools_aux_fit(ispin)%pool,&
                        name="qs_env"//TRIM(ADJUSTL(cp_to_string(qs_env%id_nr)))//&
                        "%mo_aux_fit"//TRIM(ADJUSTL(cp_to_string(ispin))),&
                        error=error)
                END IF
             END DO
          END IF
       END IF
       ALLOCATE(qs_env%rtp, stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

       CALL rt_prop_create(qs_env%rtp,qs_env%mos,qs_env%mpools,qs_env%dft_control,error)

       CALL get_restart_wfn(qs_env,error)
       qs_env%run_rtp=.TRUE.
    END SELECT

    CALL get_qs_env(qs_env=qs_env,rtp=rtp,error=error)
    aspc_order=rtp_control%aspc_order
    extrapolate_H=(rtp_control%extrapolation==extr_s_ks)
    CALL rtp_history_create(qs_env%rtp,extrapolate_H,aspc_order,error)

    motion_section => section_vals_get_subs_vals(force_env%root_section,"MOTION",error=error)
    md_section     => section_vals_get_subs_vals(motion_section,"MD",error=error)
    CALL section_vals_val_get(md_section,"TIMESTEP",r_val=qs_env%rtp%dt,error=error)
    CALL section_vals_val_get(md_section,"STEP_START_VAL",i_val=qs_env%rtp%i_start,error=error)
    CALL section_vals_val_get(md_section,"STEPS",i_val=rtp%nsteps,error=error)


    IF(rtp_control%fixed_ions)THEN
       CALL init_propagation_run(qs_env,error)
       CALL run_propagation(qs_env,force_env,globenv,error)
    ELSE
       CALL init_ehrenfest_md(force_env,qs_env,error)
    END IF


  END SUBROUTINE rt_prop_setup

! *****************************************************************************
!> \brief calculates the matrices needed in the first step of RTP
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE init_propagation_run(qs_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    INTEGER                                  :: i
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mos_new, mos_old
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(rt_prop_type), POINTER              :: rtp
    TYPE(rtp_control_type), POINTER          :: rtp_control

    CALL get_qs_env(qs_env,&
         mos=mos,&
         rtp=rtp,&
         error=error)
    rtp_control=>qs_env%dft_control%rtp_control
    CALL get_rtp(rtp=rtp,mos_old=mos_old,mos_new=mos_new,error=error)

    IF(rtp_control%initial_wfn==use_scf_wfn)THEN
      IF (rtp_control%apply_delta_pulse) THEN
         IF (rtp_control%periodic) THEN
           CALL apply_delta_pulse_periodic(qs_env,mos_old,mos_new,error)
         ELSE
           CALL apply_delta_pulse(qs_env,mos_old,mos_new,error)
         ENDIF
      ELSE
         DO i=1,SIZE(mos)
             CALL cp_fm_to_fm(mos(i)%mo_set%mo_coeff,mos_old(2*i-1)%matrix,error)
             CALL cp_fm_set_all(mos_old(2*i)%matrix,zero,zero,error)
         END DO
      ENDIF
    END IF

    DO i=1,SIZE(mos_old)
       CALL cp_fm_to_fm(mos_old(i)%matrix,mos_new(i)%matrix,error)
    END DO

!    IF(rtp_control%initial_wfn==use_restart_wfn.OR.rtp_control%initial_wfn==use_rt_restart)THEN
       CALL calc_update_rho(qs_env,error)
       CALL qs_ks_update_qs_env(qs_env%ks_env,qs_env=qs_env,calculate_forces=.FALSE.,&
            error=error)
!    END IF

    CALL init_propagators(qs_env,error)

  END SUBROUTINE init_propagation_run

! *****************************************************************************
!> \brief uses perturbation theory to get the proper initial conditions
!> \author Joost & Martin (2011)
! *****************************************************************************

  SUBROUTINE apply_delta_pulse_periodic(qs_env,mos_old,mos_new,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mos_old, mos_new
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'apply_delta_pulse_periodic', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), DIMENSION(:), POINTER  :: eigenvalues_sqrt
    INTEGER                                  :: handle, icol, idir, irow, &
                                                ispin, nao, ncol_local, nmo, &
                                                nrow_global, nrow_local, nvirt
    INTEGER, DIMENSION(:), POINTER           :: col_indices, row_indices
    REAL(KIND=dp)                            :: factor
    REAL(KIND=dp), DIMENSION(3)              :: kvec
    REAL(kind=dp), DIMENSION(:), POINTER     :: eigenvalues
    REAL(KIND=dp), DIMENSION(:, :), POINTER  :: local_data
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_cfm_type), POINTER               :: oo_c, oo_v, oo_vt
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_s
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: eigenvectors, mat_ks, &
                                                mat_tmp, momentum, oo_1, &
                                                oo_2, S_chol, S_inv, virtuals
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(rt_prop_type), POINTER              :: rtp

    CALL timeset(routineN,handle)

    ! we need the overlap and ks matrix for a full diagionalization
    CALL get_qs_env(qs_env, cell=cell, mos=mos, rtp=rtp, matrix_s=matrix_s, matrix_ks=matrix_ks, error=error)
    CALL get_rtp(rtp=rtp,S_inv=S_inv,error=error)
    CALL cp_fm_get_info(S_inv,nrow_global=nrow_global,error=error)

    CALL cp_fm_create(S_chol, matrix_struct=S_inv%matrix_struct, name="S_chol", error=error)
    CALL copy_dbcsr_to_fm(matrix_s(1)%matrix,S_chol,error=error)
    CALL cp_fm_cholesky_decompose(S_chol,error=error)
    NULLIFY(mat_ks,eigenvectors,mat_tmp)
    CALL cp_fm_create(mat_ks, matrix_struct=S_inv%matrix_struct, name="mat_ks", error=error)
    CALL cp_fm_create(eigenvectors, matrix_struct=S_inv%matrix_struct, name="eigenvectors", error=error)

    DO ispin=1,SIZE(matrix_ks)
       ALLOCATE(eigenvalues(nrow_global))
       CALL cp_fm_create(mat_tmp, matrix_struct=S_inv%matrix_struct, name="mat_tmp", error=error)

       CALL copy_dbcsr_to_fm(matrix_ks(ispin)%matrix,mat_ks,error=error)
       CALL cp_fm_cholesky_reduce(mat_ks,S_chol,error=error)
       CALL cp_fm_syevd(mat_ks,mat_tmp,eigenvalues,error=error)
       CALL cp_fm_cholesky_restore(mat_tmp,nrow_global,S_chol,eigenvectors,"SOLVE",error=error)

       ! virtuals
       CALL get_mo_set(mo_set=mos(ispin)%mo_set, nao=nao, nmo=nmo)
       nvirt=nao-nmo
       CALL cp_fm_struct_create(fm_struct_tmp, para_env=S_inv%matrix_struct%para_env, context=S_inv%matrix_struct%context,&
                                nrow_global=nrow_global, ncol_global=nvirt, error=error)
       CALL cp_fm_create(virtuals, matrix_struct=fm_struct_tmp, name="virtuals", error=error)
       CALL cp_fm_struct_release(fm_struct_tmp,error=error)
       CALL cp_fm_to_fm(eigenvectors,virtuals,nvirt,nmo+1,1)

       ! occupied
       CALL cp_fm_to_fm(eigenvectors,mos_old(2*ispin-1)%matrix,nmo,1,1)

       CALL cp_fm_struct_create(fm_struct_tmp, para_env=S_inv%matrix_struct%para_env, context=S_inv%matrix_struct%context,&
                                nrow_global=nvirt, ncol_global=nmo, error=error)
       CALL cp_fm_create(momentum, matrix_struct=fm_struct_tmp, name="momentum", error=error)
       CALL cp_fm_struct_release(fm_struct_tmp,error=error)

       ! the momentum operator (in a given direction)
       CALL cp_fm_set_all( mos_new(2*ispin-1)%matrix,0.0_dp,error=error)

       ! the prefactor (strength of the electric field)
       kvec(:) = cell%h_inv(1,:)*qs_env%dft_control%rtp_control%delta_pulse_direction(1)+&
                 cell%h_inv(2,:)*qs_env%dft_control%rtp_control%delta_pulse_direction(2)+&
                 cell%h_inv(3,:)*qs_env%dft_control%rtp_control%delta_pulse_direction(3)
       kvec = - kvec * twopi * qs_env%dft_control%rtp_control%delta_pulse_scale

       DO idir=1,3
          factor= kvec(idir)
          IF (factor.NE.0.0_dp) THEN
             CALL cp_dbcsr_sm_fm_multiply(matrix_s(idir+1)%matrix,mos_old(2*ispin-1)%matrix,&
                                          mos_old(2*ispin)%matrix,ncol=nmo,error=error)
             CALL cp_fm_scale_and_add(1.0_dp, mos_new(2*ispin-1)%matrix,factor,mos_old(2*ispin)%matrix,error)
          ENDIF
       ENDDO

       CALL cp_fm_gemm('T','N',nvirt,nmo,nao,1.0_dp, virtuals, mos_new(2*ispin-1)%matrix,0.0_dp,momentum,error=error)


       ! the tricky bit ... rescale by the eigenvalue difference
       CALL cp_fm_get_info(momentum, nrow_local=nrow_local, ncol_local=ncol_local,&
                           row_indices=row_indices,col_indices=col_indices,local_data=local_data,error=error)
       DO icol=1,ncol_local
          DO irow=1,nrow_local
             factor= 1 / (eigenvalues(col_indices(icol))-eigenvalues(nmo+row_indices(irow)))
             local_data(irow,icol)=factor*local_data(irow,icol)
          ENDDO
       ENDDO
       CALL cp_fm_release(mat_tmp,error=error)
       DEALLOCATE(eigenvalues)


       ! now obtain the initial condition in mos_old
       CALL cp_fm_to_fm(eigenvectors,mos_old(2*ispin-1)%matrix,nmo,1,1)
       CALL cp_fm_gemm("N","N",nao,nmo,nvirt,1.0_dp,virtuals,momentum,0.0_dp,mos_old(2*ispin)%matrix,error=error)

       CALL cp_fm_release(virtuals, error=error)
       CALL cp_fm_release(momentum, error=error)

       ! orthonormalize afterwards
       CALL cp_fm_struct_create(fm_struct_tmp, para_env=S_inv%matrix_struct%para_env, context=S_inv%matrix_struct%context,&
                                nrow_global=nmo, ncol_global=nmo, error=error)
       CALL cp_fm_create(oo_1, matrix_struct=fm_struct_tmp, name="oo_1", error=error)
       CALL cp_fm_create(oo_2, matrix_struct=fm_struct_tmp, name="oo_2", error=error)
       CALL cp_fm_struct_release(fm_struct_tmp, error=error)

       CALL cp_fm_create(mat_tmp, matrix_struct=mos_old(2*ispin-1)%matrix%matrix_struct, name="tmp_mat", error=error)
       ! get the complex overlap matrix
       ! x^T S x + y^T S y + i (-y^TS x+x^T S y)
       CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix,mos_old(2*ispin-1)%matrix,&
                                          mat_tmp,ncol=nmo,error=error)

       CALL cp_fm_gemm("T","N",nmo,nmo,nao,1.0_dp,mos_old(2*ispin-1)%matrix,mat_tmp,0.0_dp,oo_1,error=error)
       CALL cp_fm_gemm("T","N",nmo,nmo,nao,-1.0_dp,mos_old(2*ispin)%matrix,mat_tmp,0.0_dp,oo_2,error=error)

       CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix,mos_old(2*ispin)%matrix,&
                                          mat_tmp,ncol=nmo,error=error)
       CALL cp_fm_gemm("T","N",nmo,nmo,nao,1.0_dp,mos_old(2*ispin)%matrix,mat_tmp,1.0_dp,oo_1,error=error)
       CALL cp_fm_gemm("T","N",nmo,nmo,nao,1.0_dp,mos_old(2*ispin-1)%matrix,mat_tmp,1.0_dp,oo_2,error=error)
       CALL cp_fm_release(mat_tmp, error=error)

       CALL cp_cfm_create ( oo_c, oo_1 % matrix_struct,error=error)
       CALL cp_cfm_create ( oo_v, oo_1 % matrix_struct,error=error)
       CALL cp_cfm_create ( oo_vt, oo_1 % matrix_struct,error=error)
       oo_c % local_data = CMPLX(oo_1%local_data,oo_2%local_data,KIND=dp)

       ! compute inv(sqrt(overlap))
       ALLOCATE(eigenvalues(nmo))
       ALLOCATE(eigenvalues_sqrt(nmo))
       CALL cp_cfm_heevd(oo_c,oo_v,eigenvalues,error)
       eigenvalues_sqrt=1/SQRT(eigenvalues)
       CALL cp_cfm_to_cfm(oo_v,oo_vt,error=error)
       CALL cp_cfm_column_scale(oo_v,eigenvalues_sqrt)
       DEALLOCATE(eigenvalues)
       DEALLOCATE(eigenvalues_sqrt)
       CALL cp_cfm_gemm('N','C',nmo,nmo,nmo,(1.0_dp,0.0_dp),&
            oo_v,oo_vt,(0.0_dp,0.0_dp),oo_c,error=error)
       oo_1%local_data=REAL(oo_c%local_data,KIND=dp)
       oo_2%local_data=AIMAG(oo_c%local_data)
       CALL cp_cfm_release(oo_c,error=error)
       CALL cp_cfm_release(oo_v,error=error)
       CALL cp_cfm_release(oo_vt,error=error)

       ! use this to compute the orthonormal vectors
       CALL cp_fm_gemm("N","N",nao,nmo,nmo,1.0_dp,mos_old(2*ispin-1)%matrix,oo_1,0.0_dp,mos_new(2*ispin-1)%matrix,error=error)
       CALL cp_fm_gemm("N","N",nao,nmo,nmo,1.0_dp,mos_old(2*ispin-1)%matrix,oo_2,0.0_dp,mos_new(2*ispin)%matrix,error=error)

       CALL cp_fm_gemm("N","N",nao,nmo,nmo,-1.0_dp,mos_old(2*ispin)%matrix,oo_2,0.0_dp,mos_old(2*ispin-1)%matrix,error=error)
       CALL cp_fm_scale_and_add(1.0_dp,mos_old(2*ispin-1)%matrix,1.0_dp,mos_new(2*ispin-1)%matrix,error)

       CALL cp_fm_gemm("N","N",nao,nmo,nmo,1.0_dp,mos_old(2*ispin)%matrix,oo_1,1.0_dp,mos_new(2*ispin)%matrix,error=error)
       CALL cp_fm_to_fm(mos_new(2*ispin)%matrix,mos_old(2*ispin)%matrix,error)

       CALL cp_fm_release(oo_1,error=error)
       CALL cp_fm_release(oo_2,error=error)
    END DO

    CALL cp_fm_release(S_chol, error=error)
    CALL cp_fm_release(mat_ks, error=error)
    CALL cp_fm_release(eigenvectors, error=error)

    CALL timestop(handle)

  END SUBROUTINE apply_delta_pulse_periodic

! *****************************************************************************
!> \brief applies exp(ikr) to the wavefunction.... stored in mos_old...
!> \author Joost & Martin (2011)
! *****************************************************************************

  SUBROUTINE apply_delta_pulse(qs_env,mos_old,mos_new,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mos_old, mos_new
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'apply_delta_pulse', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp), DIMENSION(:), POINTER  :: eigenvalues_sqrt
    INTEGER                                  :: handle, i, nao, nmo
    REAL(KIND=dp), DIMENSION(3)              :: kvec
    REAL(kind=dp), DIMENSION(:), POINTER     :: eigenvalues
    TYPE(cell_type), POINTER                 :: cell
    TYPE(cp_cfm_type), POINTER               :: oo_c, oo_v, oo_vt
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s
    TYPE(cp_dbcsr_type), POINTER             :: cosmat, sinmat
    TYPE(cp_fm_struct_type), POINTER         :: fm_struct_tmp
    TYPE(cp_fm_type), POINTER                :: mat_S, oo_1, oo_2, S_inv, tmp
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(rt_prop_type), POINTER              :: rtp

    CALL timeset(routineN,handle)

    ! we need the inverse overlap

    CALL get_qs_env(qs_env, mos=mos, rtp=rtp, matrix_s=matrix_s, error=error)
    CALL get_rtp(rtp=rtp,S_inv=S_inv,error=error)

    CALL cp_fm_create(tmp, matrix_struct=S_inv%matrix_struct, name="tmp_mat", error=error)

    CALL copy_dbcsr_to_fm(matrix_s(1)%matrix,S_inv,error=error)
    CALL cp_fm_cholesky_decompose(S_inv,error=error)
    CALL cp_fm_cholesky_invert(S_inv,error=error)
    CALL cp_fm_upper_to_full(S_inv,tmp,error=error)

    CALL cp_fm_create(mat_S, matrix_struct=S_inv%matrix_struct, name="mat_S", error=error)
    CALL copy_dbcsr_to_fm(matrix_s(1)%matrix,mat_S,error=error)
    CALL cp_fm_upper_to_full(mat_S,tmp,error=error)

    CALL cp_fm_release(tmp,error=error)

    ! we need the berry matrix
    CALL get_qs_env(qs_env, cell=cell, error=error)

    ! direction ... unscaled, this will yield a exp(ikr) that is periodic with the cell
    kvec(:) = cell%h_inv(1,:)*qs_env%dft_control%rtp_control%delta_pulse_direction(1)+&
              cell%h_inv(2,:)*qs_env%dft_control%rtp_control%delta_pulse_direction(2)+&
              cell%h_inv(3,:)*qs_env%dft_control%rtp_control%delta_pulse_direction(3)
    kvec = - kvec * twopi
    ! scaling will make the things not periodic with the cell, which would only be good for gas phase systems ?
    kvec(:) = qs_env%dft_control%rtp_control%delta_pulse_scale * kvec

    ALLOCATE(cosmat, sinmat)
    CALL cp_dbcsr_init(cosmat, error=error)
    CALL cp_dbcsr_init(sinmat, error=error)
    CALL cp_dbcsr_copy(cosmat,matrix_s(1)%matrix,'COS MOM',error=error)
    CALL cp_dbcsr_copy(sinmat,matrix_s(1)%matrix,'SIN MOM',error=error)
    CALL build_berry_moment_matrix(qs_env,cosmat,sinmat,kvec,error)

    ! apply inv(S)*operator to C
    DO i=1,SIZE(mos)
       CALL get_mo_set(mos(i)%mo_set, nao=nao,nmo=nmo)
       CALL cp_dbcsr_sm_fm_multiply(cosmat, mos(i)%mo_set%mo_coeff, mos_new(2*i-1)%matrix, ncol=nmo, error=error)
       CALL cp_dbcsr_sm_fm_multiply(sinmat, mos(i)%mo_set%mo_coeff, mos_new(2*i)%matrix, ncol=nmo, error=error)

       CALL cp_fm_gemm("N","N",nao,nmo,nao,1.0_dp,S_inv,mos_new(2*i-1)%matrix,0.0_dp,mos_old(2*i-1)%matrix,error=error)
       CALL cp_fm_gemm("N","N",nao,nmo,nao,1.0_dp,S_inv,mos_new(2*i)%matrix,0.0_dp,mos_old(2*i)%matrix,error=error)

       ! in a finite basis, unfortunately, inv(S)*operator is not unitary, so orthonormalize afterwards
       CALL cp_fm_struct_create(fm_struct_tmp, para_env=S_inv%matrix_struct%para_env, context=S_inv%matrix_struct%context,&
                                nrow_global=nmo, ncol_global=nmo, error=error)
       CALL cp_fm_create(oo_1, matrix_struct=fm_struct_tmp, name="oo_1", error=error)
       CALL cp_fm_create(oo_2, matrix_struct=fm_struct_tmp, name="oo_2", error=error)
       CALL cp_fm_struct_release(fm_struct_tmp, error=error)

       CALL cp_fm_create(tmp, matrix_struct=mos_old(2*i-1)%matrix%matrix_struct, name="tmp_mat", error=error)
       ! get the complex overlap matrix
       ! x^T S x + y^T S y + i (-y^TS x+x^T S y)
       CALL cp_fm_gemm("N","N",nao,nmo,nao,1.0_dp,mat_S,mos_old(2*i-1)%matrix,0.0_dp,tmp,error=error)
       CALL cp_fm_gemm("T","N",nmo,nmo,nao,1.0_dp,mos_old(2*i-1)%matrix,tmp,0.0_dp,oo_1,error=error)
       CALL cp_fm_gemm("T","N",nmo,nmo,nao,-1.0_dp,mos_old(2*i)%matrix,tmp,0.0_dp,oo_2,error=error)

       CALL cp_fm_gemm("N","N",nao,nmo,nao,1.0_dp,mat_S,mos_old(2*i)%matrix,0.0_dp,tmp,error=error)
       CALL cp_fm_gemm("T","N",nmo,nmo,nao,1.0_dp,mos_old(2*i)%matrix,tmp,1.0_dp,oo_1,error=error)
       CALL cp_fm_gemm("T","N",nmo,nmo,nao,1.0_dp,mos_old(2*i-1)%matrix,tmp,1.0_dp,oo_2,error=error)
       CALL cp_fm_release(tmp, error=error)

       CALL cp_cfm_create ( oo_c, oo_1 % matrix_struct,error=error)
       CALL cp_cfm_create ( oo_v, oo_1 % matrix_struct,error=error)
       CALL cp_cfm_create ( oo_vt, oo_1 % matrix_struct,error=error)
       oo_c % local_data = CMPLX(oo_1%local_data,oo_2%local_data,KIND=dp)

       ! compute inv(sqrt(overlap))
       ALLOCATE(eigenvalues(nmo))
       ALLOCATE(eigenvalues_sqrt(nmo))
       CALL cp_cfm_heevd(oo_c,oo_v,eigenvalues,error)
       eigenvalues_sqrt=1/SQRT(eigenvalues)
       CALL cp_cfm_to_cfm(oo_v,oo_vt,error=error)
       CALL cp_cfm_column_scale(oo_v,eigenvalues_sqrt)
       DEALLOCATE(eigenvalues)
       DEALLOCATE(eigenvalues_sqrt)
       CALL cp_cfm_gemm('N','C',nmo,nmo,nmo,(1.0_dp,0.0_dp),&
            oo_v,oo_vt,(0.0_dp,0.0_dp),oo_c,error=error)
       oo_1%local_data=REAL(oo_c%local_data,KIND=dp)
       oo_2%local_data=AIMAG(oo_c%local_data)
       CALL cp_cfm_release(oo_c,error=error)
       CALL cp_cfm_release(oo_v,error=error)
       CALL cp_cfm_release(oo_vt,error=error)

       ! use this to compute the orthonormal vectors
       CALL cp_fm_gemm("N","N",nao,nmo,nmo,1.0_dp,mos_old(2*i-1)%matrix,oo_1,0.0_dp,mos_new(2*i-1)%matrix,error=error)
       CALL cp_fm_gemm("N","N",nao,nmo,nmo,1.0_dp,mos_old(2*i-1)%matrix,oo_2,0.0_dp,mos_new(2*i)%matrix,error=error)

       CALL cp_fm_gemm("N","N",nao,nmo,nmo,-1.0_dp,mos_old(2*i)%matrix,oo_2,0.0_dp,mos_old(2*i-1)%matrix,error=error)
       CALL cp_fm_scale_and_add(1.0_dp,mos_old(2*i-1)%matrix,1.0_dp,mos_new(2*i-1)%matrix,error)

       CALL cp_fm_gemm("N","N",nao,nmo,nmo,1.0_dp,mos_old(2*i)%matrix,oo_1,1.0_dp,mos_new(2*i)%matrix,error=error)
       CALL cp_fm_to_fm(mos_new(2*i)%matrix,mos_old(2*i)%matrix,error)

       CALL cp_fm_release(oo_1,error=error)
       CALL cp_fm_release(oo_2,error=error)
    END DO

    CALL cp_fm_release(mat_S, error=error)

    CALL cp_dbcsr_deallocate_matrix(cosmat, error)
    CALL cp_dbcsr_deallocate_matrix(sinmat, error)


    CALL timestop(handle)

  END SUBROUTINE apply_delta_pulse

! *****************************************************************************
!> \brief calculates the matrices needed in the first step of EMD
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE init_ehrenfest_md(force_env,qs_env,error)

    TYPE(force_env_type), POINTER            :: force_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    INTEGER                                  :: i
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mos_new, mos_old
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(rt_prop_type), POINTER              :: rtp
    TYPE(rtp_control_type), POINTER          :: rtp_control

     CALL get_qs_env(qs_env,&
         mos=mos,&
         rtp=rtp,&
         energy=energy,&
         error=error)
     rtp_control=>qs_env%dft_control%rtp_control
     CALL get_rtp(rtp=rtp,mos_old=mos_old,mos_new=mos_new,error=error)

     IF(rtp_control%initial_wfn==use_scf_wfn)THEN
        DO i=1,SIZE(mos)
            CALL cp_fm_to_fm(mos(i)%mo_set%mo_coeff,mos_old(2*i-1)%matrix,error)
            CALL cp_fm_set_all(mos_old(2*i)%matrix,zero,zero,error)
        END DO
     END IF

     DO i=1,SIZE(mos_old)
        CALL cp_fm_to_fm(mos_old(i)%matrix,mos_new(i)%matrix,error)
     END DO

     qs_env%dft_control%rtp_control%initial_step=.TRUE.

     CALL force_env_calc_energy_force(force_env,calc_force=.TRUE.,&
          error=error)

     qs_env%dft_control%rtp_control%initial_step=.FALSE.
     rtp%energy_old=energy%total


  END SUBROUTINE init_ehrenfest_md

! *****************************************************************************
!> \brief performes the real RTP run, gets information from MD section
!>        uses MD as iteration level
!> \author Florian Schiffmann (02.09)
! *****************************************************************************

  SUBROUTINE run_propagation(qs_env,force_env,globenv,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'run_propagation', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: aspc_order, handle, i, &
                                                i_iter, i_step, ihist, &
                                                max_iter, max_steps
    LOGICAL                                  :: failure, should_stop
    REAL(Kind=dp)                            :: delta_iter, eps_ener, &
                                                time_iter_start, &
                                                time_iter_stop, used_time
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: exp_H_new, exp_H_old, &
                                                mos_new, mos_old
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(rt_prop_type), POINTER              :: rtp

    failure=.FALSE.
    should_stop=.FALSE.
    CALL timeset(routineN,handle)
    logger   => cp_error_get_logger(error)
    max_steps=qs_env%rtp%nsteps
    max_iter=qs_env%dft_control%rtp_control%max_iter
    eps_ener=qs_env%dft_control%rtp_control%eps_ener


    CALL get_qs_env(qs_env=qs_env,rtp=rtp,energy=energy,error=error)
    aspc_order=qs_env%dft_control%rtp_control%aspc_order

    CALL get_rtp(rtp=rtp,&
                 exp_H_old=exp_H_old,&
                 exp_H_new=exp_H_new,&
                 mos_old=mos_old,&
                 mos_new=mos_new,&
                 error=error)

    rtp%energy_old=energy%total
    time_iter_start=m_walltime()
    CALL cp_add_iter_level(logger%iter_info,"MD",error=error)
    CALL cp_iterate(logger%iter_info,iter_nr=0,error=error)
    DO i_step=rtp%i_start+1,max_steps
       qs_env%energy%efield_core=0.0_dp

       qs_env%sim_time=REAL(i_step,dp)*rtp%dt
       qs_env%sim_step=i_step
       rtp%istep=i_step-rtp%i_start
       CALL calculate_ecore_efield(qs_env,.FALSE.,error)
       CALL external_c_potential(qs_env,calculate_forces=.FALSE.,error=error)
       CALL external_e_potential(qs_env,error=error)

       CALL cp_iterate(logger%iter_info,last=(i_step==max_steps),iter_nr=i_step,error=error)
       rtp%converged=.FALSE.
       DO i_iter=1,max_iter
          IF(i_step==rtp%i_start+1.AND.i_iter==2.AND.qs_env%dft_control%rtp_control%hfx_redistribute)&
              CALL qs_ks_did_change(qs_env%ks_env,s_mstruct_changed=.TRUE., error=error)
          rtp%iter=i_iter
          CALL propagation_step(qs_env,delta_iter, error)
          rtp%energy_new=energy%total
          rtp%converged=(delta_iter.LT.eps_ener)

          IF(rtp%converged)EXIT
          CALL rt_prop_output(qs_env,real_time_propagation,delta_iter,error=error)

       END DO

       IF(rtp%converged)THEN
          ihist=MOD(rtp%istep,aspc_order)+1
          DO i=1,SIZE(mos_new)
             CALL cp_fm_to_fm(mos_new(i)%matrix,mos_old(i)%matrix,error)
             CALL cp_fm_to_fm(exp_H_new(i)%matrix,exp_H_old(i)%matrix,error)
          END DO
          IF(rtp%history%mos_or_H==2)&
             CALL  put_data_to_history(qs_env,exp_H=exp_H_new,ihist=ihist,error=error)
          IF(rtp%history%mos_or_H==1)&
             CALL put_data_to_history(qs_env,mos=mos_new,ihist=ihist,error=error)
          CALL external_control(should_stop,"MD",globenv=globenv,error=error)
          IF (should_stop)CALL cp_iterate(logger%iter_info,last=.TRUE.,iter_nr=i_step,error=error)
          time_iter_stop=m_walltime()
          used_time= time_iter_stop - time_iter_start
          time_iter_start=time_iter_stop
          CALL rt_prop_output(qs_env,real_time_propagation,delta_iter=delta_iter,used_time=used_time,error=error)
          CALL rt_write_input_restart(force_env=force_env,error=error)
          IF (should_stop)EXIT
       ELSE
          EXIT
       END IF


    END DO
    CALL cp_rm_iter_level(logger%iter_info,"MD",error=error)
    IF(.NOT.rtp%converged)&
         CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
         routineP,"propagation did not converge, either increase MAX_ITER or use a smaller TIMESTEP",&
         error,failure)
    CALL timestop(handle)

  END SUBROUTINE run_propagation

  SUBROUTINE rt_write_input_restart(md_env,force_env,error)
    TYPE(md_environment_type), OPTIONAL, &
      POINTER                                :: md_env
    TYPE(force_env_type), POINTER            :: force_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rt_write_input_restart', &
      routineP = moduleN//':'//routineN

    TYPE(section_vals_type), POINTER         :: motion_section, root_section, &
                                                rt_section

    root_section => force_env%root_section
    motion_section => section_vals_get_subs_vals(root_section,"MOTION",error=error)
    rt_section => section_vals_get_subs_vals(root_section,"FORCE_EVAL%DFT%REAL_TIME_PROPAGATION",error=error)
    CALL section_vals_val_set(rt_section,"INITIAL_WFN",i_val=use_rt_restart,error=error)
    ! coming from RTP
    IF (.NOT. PRESENT(md_env)) THEN
     CALL section_vals_val_set(motion_section,"MD%STEP_START_VAL",i_val=force_env%qs_env%sim_step,error=error)
    ENDIF

    CALL write_restart(md_env=md_env,root_section=root_section,error=error)

  END SUBROUTINE rt_write_input_restart

END MODULE rt_propagation
