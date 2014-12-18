
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Utilities for rtp in combination with admm methods
!>        adapted routines from admm_method (author Manuel Guidon)
!>
!> \par History    Use new "force only" overlap routine [07.2014,JGH]
!> \author Florian Schiffmann
! *****************************************************************************
MODULE rtp_admm_methods
  USE admm_types,                      ONLY: admm_env_create,&
                                             admm_type
  USE cp_control_types,                ONLY: admm_control_type,&
                                             dft_control_type
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_add, cp_dbcsr_col_block_sizes, cp_dbcsr_copy, &
       cp_dbcsr_create, cp_dbcsr_deallocate_matrix, cp_dbcsr_desymmetrize, &
       cp_dbcsr_distribution, cp_dbcsr_get_data_size, cp_dbcsr_get_data_type, &
       cp_dbcsr_get_matrix_type, cp_dbcsr_get_num_blocks, cp_dbcsr_init, &
       cp_dbcsr_p_type, cp_dbcsr_release, cp_dbcsr_row_block_sizes, &
       cp_dbcsr_set, cp_dbcsr_type, dbcsr_type_no_symmetry
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                             copy_fm_to_dbcsr,&
                                             cp_dbcsr_plus_fm_fm_t
  USE cp_fm_basic_linalg,              ONLY: cp_fm_upper_to_full
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                             cp_fm_cholesky_invert
  USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                             cp_fm_p_type,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_gemm_interface,               ONLY: cp_gemm
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE hfx_admm_utils,                  ONLY: create_admm_xc_section
  USE input_constants,                 ONLY: do_admm_basis_projection,&
                                             do_admm_purify_none,&
                                             use_aux_fit_basis_set
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: one,&
                                             zero
  USE particle_types,                  ONLY: particle_type
  USE pw_types,                        ONLY: pw_p_type
  USE qs_collocate_density,            ONLY: calculate_rho_elec
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type
  USE qs_rho_types,                    ONLY: qs_rho_get,&
                                             qs_rho_set,&
                                             qs_rho_type
  USE rt_propagation_types,            ONLY: get_rtp,&
                                             rt_prop_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  ! *** Public subroutines ***
  PUBLIC :: rtp_admm_calc_rho_aux,rtp_admm_merge_ks_matrix

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rtp_admm_methods'

CONTAINS

! *****************************************************************************
!> \brief  Compute the ADMM density matrix in case of rtp (complex MO's)
!>
!> \param qs_env ...
!> \param error ...
!> \par History
! *****************************************************************************
  SUBROUTINE  rtp_admm_calc_rho_aux(qs_env,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rtp_admm_calc_rho_aux', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, nspins
    LOGICAL                                  :: failure, s_mstruct_changed
    REAL(KIND=dp), DIMENSION(:), POINTER     :: tot_rho_r_aux
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_p_aux, &
                                                matrix_p_aux_im, &
                                                matrix_s_aux_fit, &
                                                matrix_s_aux_fit_vs_orb
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: rtp_coeff_aux_fit
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos, mos_aux_fit
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho_g_aux, rho_r_aux
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_rho_type), POINTER               :: rho, rho_aux_fit
    TYPE(rt_prop_type), POINTER              :: rtp

    failure=.FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(admm_env, matrix_p_aux, matrix_p_aux_im, mos, &
         mos_aux_fit, para_env,matrix_s_aux_fit,matrix_s_aux_fit_vs_orb, rho,&
         ks_env,dft_control,tot_rho_r_aux,rho_r_aux, rho_g_aux)

    CALL get_qs_env(qs_env,&
                    admm_env=admm_env,&
                    ks_env=ks_env,&
                    dft_control=dft_control,&
                    para_env=para_env,&
                    matrix_s_aux_fit=matrix_s_aux_fit,&
                    matrix_s_aux_fit_vs_orb=matrix_s_aux_fit_vs_orb,&
                    mos=mos,&
                    mos_aux_fit=mos_aux_fit,&
                    rtp=rtp,&
                    rho=rho,&
                    rho_aux_fit=rho_aux_fit,&
                    s_mstruct_changed=s_mstruct_changed,&
                    error=error)

    nspins = dft_control%nspins

    CALL get_rtp(rtp=rtp,admm_mos=rtp_coeff_aux_fit,error=error)
    CALL rtp_admm_fit_mo_coeffs(qs_env, admm_env, dft_control%admm_control, para_env, &
                                matrix_s_aux_fit, matrix_s_aux_fit_vs_orb,&
                                mos, mos_aux_fit,rtp,rtp_coeff_aux_fit,&
                                s_mstruct_changed, error=error)

    DO ispin=1,nspins
       CALL qs_rho_get(rho_aux_fit,&
                       rho_ao=matrix_p_aux,&
                       rho_ao_im=matrix_p_aux_im,&
                       rho_r=rho_r_aux,&
                       rho_g=rho_g_aux,&
                       tot_rho_r=tot_rho_r_aux,&
                       error=error)

       CALL rtp_admm_calculate_dm(admm_env,rtp_coeff_aux_fit,&
                                  matrix_p_aux(ispin)%matrix,&
                                  matrix_p_aux_im(ispin)%matrix,&
                                  ispin,error)

       CALL calculate_rho_elec(matrix_p=matrix_p_aux(ispin)%matrix,&
                               rho=rho_r_aux(ispin),&
                               rho_gspace=rho_g_aux(ispin),&
                               total_rho=tot_rho_r_aux(ispin),&
                               ks_env=ks_env,soft_valid=.FALSE.,&
                               basis_set_id=use_aux_fit_basis_set, error=error)
    END DO
    CALL set_qs_env(qs_env,admm_env=admm_env,error=error)
    CALL qs_rho_set(rho_aux_fit, rho_r_valid=.TRUE., rho_g_valid=.TRUE., error=error)

    CALL timestop(handle)

  END SUBROUTINE rtp_admm_calc_rho_aux

! *****************************************************************************
!> \brief ...
!> \param admm_env ...
!> \param rtp_coeff_aux_fit ...
!> \param density_matrix_aux ...
!> \param density_matrix_aux_im ...
!> \param ispin ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE rtp_admm_calculate_dm(admm_env,rtp_coeff_aux_fit, density_matrix_aux,&
                                   density_matrix_aux_im, ispin,error)
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: rtp_coeff_aux_fit
    TYPE(cp_dbcsr_type), POINTER             :: density_matrix_aux, &
                                                density_matrix_aux_im
    INTEGER, INTENT(in)                      :: ispin
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rtp_admm_calculate_dm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure

    failure = .FALSE.

    CALL timeset(routineN,handle)

    SELECT CASE(admm_env%purification_method)
    CASE(do_admm_purify_none)
       CALL calculate_rtp_admm_density(density_matrix_aux,density_matrix_aux_im,&
            rtp_coeff_aux_fit,ispin,error)
    CASE DEFAULT
       CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
                      "only purification NONE possible with RTP/EMD at the moment"//&
                       CPSourceFileRef,&
                       only_ionode=.TRUE.)
    END SELECT

    CALL timestop(handle)

  END SUBROUTINE rtp_admm_calculate_dm

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param admm_env ...
!> \param admm_control ...
!> \param para_env ...
!> \param matrix_s_aux_fit ...
!> \param matrix_s_mixed ...
!> \param mos ...
!> \param mos_aux_fit ...
!> \param rtp ...
!> \param rtp_coeff_aux_fit ...
!> \param geometry_did_change ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE rtp_admm_fit_mo_coeffs(qs_env, admm_env, admm_control, para_env, matrix_s_aux_fit, matrix_s_mixed, &
                                mos, mos_aux_fit,rtp,rtp_coeff_aux_fit, geometry_did_change, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(admm_control_type), POINTER         :: admm_control
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s_aux_fit, &
                                                matrix_s_mixed
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos, mos_aux_fit
    TYPE(rt_prop_type), POINTER              :: rtp
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: rtp_coeff_aux_fit
    LOGICAL, INTENT(IN)                      :: geometry_did_change
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rtp_admm_fit_mo_coeffs', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, natoms
    LOGICAL                                  :: recalc_S
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER         :: input, xc_section

    CALL timeset(routineN,handle)

    NULLIFY(xc_section, particle_set)

    IF (.NOT.(ASSOCIATED(admm_env) )) THEN
      ! setup admm environment
      CALL get_qs_env(qs_env, input=input, particle_set=particle_set, error=error)
      natoms = SIZE(particle_set,1)
      CALL admm_env_create(admm_env, admm_control, mos, mos_aux_fit, &
                           para_env, natoms, error)
      xc_section => section_vals_get_subs_vals(input,"DFT%XC",error=error)
      CALL create_admm_xc_section(qs_env=qs_env, xc_section=xc_section, &
                                  admm_env=admm_env, error=error)

      IF(admm_control%method /= do_admm_basis_projection) THEN
           CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
                "RTP requires BASIS_PROJECTION."//&
                CPSourceFileRef, only_ionode=.TRUE.)
      END IF
    END IF

    recalc_S=geometry_did_change.OR.(rtp%iter==0.AND.(rtp%istep==rtp%i_start))

    SELECT CASE(admm_env%purification_method)
    CASE(do_admm_purify_none)
       CALL rtp_fit_mo_coeffs_none(qs_env, admm_env, para_env, matrix_s_aux_fit, matrix_s_mixed, &
            mos, mos_aux_fit,rtp,rtp_coeff_aux_fit,recalc_S , error)
    CASE DEFAULT
       CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
            "Purification method not implemented in combination with RTP"//&
            CPSourceFileRef,&
            only_ionode=.TRUE.)
     END SELECT

    CALL timestop(handle)

  END SUBROUTINE rtp_admm_fit_mo_coeffs
! *****************************************************************************
!> \brief Calculates the MO coefficients for the auxiliary fitting basis set
!>        by minimizing int (psi_i - psi_aux_i)^2 using Lagrangian Multipliers
!>
!> \param qs_env ...
!> \param admm_env The ADMM env
!> \param para_env The parallel env
!> \param matrix_s_aux_fit the overlap matrix of the auxiliary fitting basis set
!> \param matrix_s_mixed the mixed overlap matrix of the auxiliary fitting basis
!>        set and the orbital basis set
!> \param mos the MO's of the orbital basis set
!> \param mos_aux_fit the MO's of the auxiliary fitting basis set
!> \param rtp ...
!> \param rtp_coeff_aux_fit ...
!> \param geometry_did_change flag to indicate if the geomtry changed
!> \param error ...
!> \par History
!>      05.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE rtp_fit_mo_coeffs_none(qs_env, admm_env, para_env, matrix_s_aux_fit, matrix_s_mixed, &
                                mos, mos_aux_fit,rtp,rtp_coeff_aux_fit, geometry_did_change, error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s_aux_fit, &
                                                matrix_s_mixed
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos, mos_aux_fit
    TYPE(rt_prop_type), POINTER              :: rtp
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: rtp_coeff_aux_fit
    LOGICAL, INTENT(IN)                      :: geometry_did_change
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rtp_fit_mo_coeffs_none', &
      routineP = moduleN//':'//routineN
    COMPLEX(KIND=dp), PARAMETER              :: one = (1.0_dp,0.0_dp) , &
                                                zero = (0.0_dp,0.0_dp)

    INTEGER                                  :: handle, ispin, nao_aux_fit, &
                                                nao_orb, natoms, nmo, &
                                                nmo_mos, nspins
    REAL(KIND=dp), DIMENSION(:), POINTER     :: occ_num, occ_num_aux
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mos_new
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_coeff_aux_fit
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(section_vals_type), POINTER         :: input, xc_section

    CALL timeset(routineN,handle)

    NULLIFY(dft_control, particle_set)

    IF (.NOT.(ASSOCIATED(admm_env) )) THEN
      CALL get_qs_env(qs_env,input=input,particle_set=particle_set,dft_control=dft_control,error=error)
      natoms = SIZE(particle_set,1)
      CALL admm_env_create(admm_env, dft_control%admm_control, mos, mos_aux_fit, para_env, natoms, error)
      xc_section => section_vals_get_subs_vals(input,"DFT%XC",error=error)
      CALL create_admm_xc_section(qs_env=qs_env, xc_section=xc_section, & 
                                  admm_env=admm_env,error=error)
    END IF

    nao_aux_fit = admm_env%nao_aux_fit
    nao_orb = admm_env%nao_orb
    nspins = SIZE(mos)


    ! *** This part only depends on overlap matrices ==> needs only to be calculated if the geometry changed

    IF( geometry_did_change ) THEN
      CALL copy_dbcsr_to_fm(matrix_s_aux_fit(1)%matrix,admm_env%S_inv,error)
      CALL cp_fm_upper_to_full(admm_env%S_inv,admm_env%work_aux_aux,error=error)
      CALL cp_fm_to_fm(admm_env%S_inv, admm_env%S, error=error)

      CALL copy_dbcsr_to_fm(matrix_s_mixed(1)%matrix,admm_env%Q,error)

      !! Calculate S'_inverse
      CALL cp_fm_cholesky_decompose(admm_env%S_inv,error=error)
      CALL cp_fm_cholesky_invert(admm_env%S_inv,error=error)
      !! Symmetrize the guy
      CALL cp_fm_upper_to_full(admm_env%S_inv,admm_env%work_aux_aux,error=error)
      !! Calculate A=S'^(-1)*P
      CALL cp_gemm('N','N',nao_aux_fit,nao_orb,nao_aux_fit,&
                    1.0_dp,admm_env%S_inv,admm_env%Q,0.0_dp,&
                    admm_env%A,error)
   END IF

    ! *** Calculate the mo_coeffs for the fitting basis
    DO ispin=1,nspins
      nmo = admm_env%nmo(ispin)
      IF( nmo == 0 ) CYCLE
      !! Lambda = C^(T)*B*C
      CALL get_rtp (rtp=rtp,mos_new=mos_new,error=error)
      CALL get_mo_set(mos(ispin)%mo_set,mo_coeff=mo_coeff, occupation_numbers=occ_num,nmo=nmo_mos)
      CALL get_mo_set(mos_aux_fit(ispin)%mo_set,mo_coeff=mo_coeff_aux_fit,&
           occupation_numbers=occ_num_aux)

      CALL cp_gemm('N','N',nao_aux_fit,nmo,nao_orb,&
                      1.0_dp,admm_env%A,mos_new(2*ispin-1)%matrix,0.0_dp,&
                      rtp_coeff_aux_fit(2*ispin-1)%matrix,error)
      CALL cp_gemm('N','N',nao_aux_fit,nmo,nao_orb,&
                      1.0_dp,admm_env%A,mos_new(2*ispin)%matrix,0.0_dp,&
                      rtp_coeff_aux_fit(2*ispin)%matrix,error)

      CALL cp_fm_to_fm(rtp_coeff_aux_fit(2*ispin-1)%matrix,mo_coeff_aux_fit,error)
    END DO

    CALL timestop(handle)

  END SUBROUTINE rtp_fit_mo_coeffs_none


! *****************************************************************************
!> \brief ...
!> \param density_matrix_aux ...
!> \param density_matrix_aux_im ...
!> \param rtp_coeff_aux_fit ...
!> \param ispin ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE calculate_rtp_admm_density(density_matrix_aux,density_matrix_aux_im,&
                                   rtp_coeff_aux_fit,ispin, error)

    TYPE(cp_dbcsr_type), POINTER             :: density_matrix_aux, &
                                                density_matrix_aux_im
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: rtp_coeff_aux_fit
    INTEGER, INTENT(in)                      :: ispin
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_rtp_admm_density', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    INTEGER                                  :: handle, im, ncol, re
    REAL(KIND=dp)                            :: alpha

    CALL timeset(routineN,handle)

    re=2*ispin-1 ; im =2*ispin
    alpha=3.0_dp-REAL(SIZE(rtp_coeff_aux_fit)/2,dp)
    CALL cp_dbcsr_set(density_matrix_aux,0.0_dp,error=error)
    CALL cp_fm_get_info(rtp_coeff_aux_fit(re)%matrix,ncol_global=ncol,error=error)
    CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=density_matrix_aux,&
                               matrix_v=rtp_coeff_aux_fit(re)%matrix,&
                               matrix_g=rtp_coeff_aux_fit(re)%matrix,&
                               ncol=ncol,&
                               alpha=alpha,error=error)

    ! It is actually complex conjugate but i*i=-1 therfore it must be added
    CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=density_matrix_aux,&
                               matrix_v=rtp_coeff_aux_fit(im)%matrix,&
                               matrix_g=rtp_coeff_aux_fit(im)%matrix,&
                               ncol=ncol,&
                               alpha=alpha,error=error)

!   compute the imaginary part of the dm
    CALL cp_dbcsr_set(density_matrix_aux_im,0.0_dp,error=error)
    CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=density_matrix_aux_im,&
                               matrix_v=rtp_coeff_aux_fit(im)%matrix,&
                               matrix_g=rtp_coeff_aux_fit(re)%matrix,&
                               ncol=ncol,&
                               alpha=alpha,error=error)
    alpha=-alpha
    CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=density_matrix_aux_im,&
                               matrix_v=rtp_coeff_aux_fit(re)%matrix,&
                               matrix_g=rtp_coeff_aux_fit(im)%matrix,&
                               ncol=ncol,&
                               alpha=alpha,error=error)

    CALL timestop(handle)

  END SUBROUTINE calculate_rtp_admm_density

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE rtp_admm_merge_ks_matrix(qs_env, error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rtp_admm_merge_ks_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin
    LOGICAL                                  :: failure
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_ks_aux_fit, &
                                                matrix_ks_aux_fit_im, &
                                                matrix_ks_im
    TYPE(dft_control_type), POINTER          :: dft_control

    failure = .FALSE.
    NULLIFY(admm_env, dft_control,&
            matrix_ks, matrix_ks_im, matrix_ks_aux_fit, matrix_ks_aux_fit_im)
    CALL timeset(routineN,handle)

    CALL get_qs_env(qs_env,&
                    admm_env=admm_env,&
                    dft_control=dft_control,&
                    matrix_ks=matrix_ks,&
                    matrix_ks_im=matrix_ks_im,&
                    matrix_ks_aux_fit=matrix_ks_aux_fit,&
                    matrix_ks_aux_fit_im=matrix_ks_aux_fit_im,&
                    error=error)

    DO ispin=1, dft_control%nspins

       SELECT CASE(admm_env%purification_method)
       CASE(do_admm_purify_none)
          CALL rt_merge_ks_matrix_none(ispin, admm_env, &
                                       matrix_ks, matrix_ks_aux_fit, error)
          CALL rt_merge_ks_matrix_none(ispin, admm_env,&
                                       matrix_ks_im, matrix_ks_aux_fit_im, error)
       CASE DEFAULT
          CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
                         "only purification NONE possible with RTP/EMD at the moment"//&
                          CPSourceFileRef, only_ionode=.TRUE.)
       END SELECT

    ENDDO !spin loop
    CALL timestop(handle)

  END SUBROUTINE rtp_admm_merge_ks_matrix


! *****************************************************************************
!> \brief ...
!> \param ispin ...
!> \param admm_env ...
!> \param matrix_ks ...
!> \param matrix_ks_aux_fit ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE rt_merge_ks_matrix_none(ispin, admm_env,  &
                                  matrix_ks, matrix_ks_aux_fit,  error)
    INTEGER, INTENT(IN)                      :: ispin
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_ks_aux_fit
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rt_merge_ks_matrix_none', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nao_aux_fit, nao_orb, &
                                                nmo
    INTEGER, SAVE                            :: counter = 0
    LOGICAL                                  :: failure
    TYPE(cp_dbcsr_type)                      :: matrix_ks_nosym
    TYPE(cp_dbcsr_type), POINTER             :: matrix_k_tilde

    failure = .FALSE.

    CALL timeset(routineN,handle)

    counter = counter + 1
    nao_aux_fit = admm_env%nao_aux_fit
    nao_orb = admm_env%nao_orb
    nmo = admm_env%nmo(ispin)
    CALL cp_dbcsr_init (matrix_ks_nosym, error)
    CALL cp_dbcsr_create (matrix_ks_nosym, template=matrix_ks_aux_fit(ispin)%matrix,&
                 matrix_type=dbcsr_type_no_symmetry, error=error)
    CALL cp_dbcsr_set(matrix_ks_nosym,0.0_dp,error)
    CALL cp_dbcsr_desymmetrize(matrix_ks_aux_fit(ispin)%matrix,matrix_ks_nosym,error)

    CALL copy_dbcsr_to_fm(matrix_ks_nosym,admm_env%K(ispin)%matrix,error)

    !! K*A
    CALL cp_gemm('N', 'N',  nao_aux_fit, nao_orb, nao_aux_fit,&
                    1.0_dp,admm_env%K(ispin)%matrix,admm_env%A,0.0_dp,&
                    admm_env%work_aux_orb,error)
    !! A^T*K*A
    CALL cp_gemm('T', 'N',  nao_orb, nao_orb, nao_aux_fit,&
                    1.0_dp,admm_env%A,admm_env%work_aux_orb,0.0_dp,&
                    admm_env%work_orb_orb,error)


    NULLIFY(matrix_k_tilde)
    ALLOCATE(matrix_k_tilde)
    CALL cp_dbcsr_init (matrix_k_tilde, error)
    CALL cp_dbcsr_create(matrix_k_tilde, 'MATRIX K_tilde', &
         cp_dbcsr_distribution(matrix_ks(ispin)%matrix), cp_dbcsr_get_matrix_type(matrix_ks_aux_fit(ispin)%matrix),&
         cp_dbcsr_row_block_sizes(matrix_ks(ispin)%matrix),&
         cp_dbcsr_col_block_sizes(matrix_ks(ispin)%matrix), &
         cp_dbcsr_get_num_blocks(matrix_ks(ispin)%matrix), cp_dbcsr_get_data_size(matrix_ks(ispin)%matrix),&
         cp_dbcsr_get_data_type(matrix_ks(ispin)%matrix), error=error)

    CALL cp_dbcsr_copy(matrix_k_tilde, matrix_ks(ispin)%matrix, error=error)
    CALL cp_dbcsr_set(matrix_k_tilde, 0.0_dp, error)
    CALL copy_fm_to_dbcsr(admm_env%work_orb_orb, matrix_k_tilde, keep_sparsity=.TRUE.,&
         error=error)

    CALL cp_dbcsr_add(matrix_ks(ispin)%matrix, matrix_k_tilde, 1.0_dp, 1.0_dp, error)

    CALL cp_dbcsr_deallocate_matrix(matrix_k_tilde,error)
    CALL cp_dbcsr_release(matrix_ks_nosym,error)

    CALL timestop(handle)

  END SUBROUTINE rt_merge_ks_matrix_none

END MODULE rtp_admm_methods
