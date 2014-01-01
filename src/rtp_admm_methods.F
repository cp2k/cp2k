
!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Utilities for rtp in combination with admm methods
!>        adapted routines from admm_method (author Manuel Guidon)
!>
!> \par History
!> \author Florian Schiffmann
! *****************************************************************************
MODULE rtp_admm_methods
  USE admm_methods,                    ONLY: create_admm_xc_section
  USE admm_types,                      ONLY: admm_env_create,&
                                             admm_type
  USE cp_control_types,                ONLY: admm_control_type,&
                                             dft_control_type
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_add, cp_dbcsr_col_block_sizes, cp_dbcsr_copy, &
       cp_dbcsr_create, cp_dbcsr_distribution, cp_dbcsr_get_data_size, &
       cp_dbcsr_get_data_type, cp_dbcsr_get_num_blocks, cp_dbcsr_init, &
       cp_dbcsr_row_block_sizes, cp_dbcsr_scale, cp_dbcsr_set
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                             copy_fm_to_dbcsr,&
                                             cp_dbcsr_deallocate_matrix,&
                                             cp_dbcsr_plus_fm_fm_t
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_p_type,&
                                             cp_dbcsr_type
  USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                             cp_fm_gemm,&
                                             cp_fm_schur_product,&
                                             cp_fm_upper_to_full
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                             cp_fm_cholesky_invert,&
                                             cp_fm_cholesky_reduce,&
                                             cp_fm_cholesky_restore
  USE cp_fm_diag,                      ONLY: cp_fm_syevd
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_p_type,&
                                             cp_fm_set_all,&
                                             cp_fm_set_element,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_fm_vect,                      ONLY: cp_fm_vect_dealloc
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE dbcsr_types,                     ONLY: dbcsr_type_symmetric
  USE input_constants,                 ONLY: do_admm_block_density_matrix,&
                                             do_admm_purify_cauchy,&
                                             do_admm_purify_none,&
                                             use_aux_fit_basis_set
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE mathconstants
  USE qs_collocate_density,            ONLY: calculate_rho_elec
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE qs_mo_methods,                   ONLY: calculate_density_matrix
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type,&
                                             mo_set_type
  USE rt_propagation_types,            ONLY: get_rtp,&
                                             rt_prop_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  ! *** Public subroutines ***
  PUBLIC :: rtp_admm_density_matrix,rtp_admm_merge_ks_matrix

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rtp_admm_methods'

CONTAINS

! *****************************************************************************
!> \brief  Compute the ADMM density matrix in case of rtp (complex MO's)
!>
!> \par History
! *****************************************************************************
  SUBROUTINE  rtp_admm_density_matrix(qs_env,ks_env,dft_control,input,error)

    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(section_vals_type), POINTER         :: input
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rtp_admm_density_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, n_rep_hf, &
                                                nspins
    LOGICAL                                  :: do_hfx, failure
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_p, matrix_p_aux_fit
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: rtp_coeff_aux_fit
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos, mos_aux_fit
    TYPE(rt_prop_type), POINTER              :: rtp
    TYPE(section_vals_type), POINTER         :: hfx_sections

    failure=.FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(admm_env, hfx_sections, matrix_p, matrix_p_aux_fit, mos, &
         mos_aux_fit, para_env)

    hfx_sections => section_vals_get_subs_vals(input,"DFT%XC%HF",error=error)
    CALL section_vals_get(hfx_sections,explicit=do_hfx,error=error)
    !! ** ADMM can only be used with HFX
    IF ( .NOT. do_hfx ) THEN
       CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
                      "Wavefunction fitting requested without Hartree-Fock."//&
                       CPSourceFileRef,&
                       only_ionode=.TRUE.)
    END IF

    ! ** Method only implemented for GPW
    IF( dft_control%qs_control%gapw ) THEN
      CALL cp_unimplemented_error(fromWhere=routineP, &
                               message="ADMM only implemented for GPW", &
                               error=error, error_level=cp_failure_level)
    END IF

    nspins=dft_control%nspins
    CALL get_qs_env(qs_env,mos_aux_fit=mos_aux_fit,&
                    rtp=rtp,mos=mos,admm_env=admm_env,&
                    para_env=para_env,error=error)
    CALL section_vals_get(hfx_sections,n_repetition=n_rep_hf,error=error)
    IF ( n_rep_hf > 1 ) THEN
       CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
                      "ADMM can handle only one HF section."//&
                       CPSourceFileRef,&
                       only_ionode=.TRUE.)
    END IF
    CALL get_rtp(rtp=rtp,admm_mos=rtp_coeff_aux_fit,error=error)
    CALL rtp_admm_fit_mo_coeffs(qs_env, admm_env, dft_control%admm_control, para_env, &
                            qs_env%matrix_s_aux_fit, qs_env%matrix_s_aux_fit_vs_orb,&
                            mos, mos_aux_fit,rtp,rtp_coeff_aux_fit,&
                            ks_env%s_mstruct_changed, error=error)

    DO ispin=1,nspins
       matrix_p_aux_fit =>qs_env%rho_aux_fit%rho_ao
       matrix_p => qs_env%rho%rho_ao
       CALL rtp_admm_calculate_dm(admm_env,rtp_coeff_aux_fit,mos_aux_fit(ispin)%mo_set,&
                                       matrix_p_aux_fit(ispin)%matrix,&
                                       ispin,nspins,error)
       CALL calculate_rho_elec(matrix_p=qs_env%rho_aux_fit%rho_ao(ispin)%matrix,&
            rho=qs_env%rho_aux_fit%rho_r(ispin),&
            rho_gspace=qs_env%rho_aux_fit%rho_g(ispin),&
            total_rho=qs_env%rho_aux_fit%tot_rho_r(ispin),&
            qs_env=qs_env,soft_valid=.FALSE.,&
            basis_set_id=use_aux_fit_basis_set, error=error)
    END DO
    CALL cp_fm_vect_dealloc(rtp_coeff_aux_fit,error)
    CALL set_qs_env(qs_env,admm_env=admm_env,error=error)

    qs_env%rho_aux_fit%rho_r_valid=.TRUE.
    qs_env%rho_aux_fit%rho_g_valid=.TRUE.

    CALL timestop(handle)

  END SUBROUTINE rtp_admm_density_matrix

  SUBROUTINE rtp_admm_calculate_dm(admm_env,rtp_coeff_aux_fit,mo_set, density_matrix_aux,&
                                           ispin, nspins, error)
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: rtp_coeff_aux_fit
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_dbcsr_type), POINTER             :: density_matrix_aux
    INTEGER                                  :: ispin, nspins
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rtp_admm_calculate_dm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure

    failure = .FALSE.

    CALL timeset(routineN,handle)

    SELECT CASE(admm_env%purification_method)
    CASE(do_admm_purify_none)
       CALL calculate_rtp_admm_density(mo_set,&
            density_matrix_aux,&
            rtp_coeff_aux_fit,ispin,error)
    CASE(do_admm_purify_cauchy)
       CALL rtp_purify_dm_cauchy(admm_env,rtp_coeff_aux_fit,mo_set,density_matrix_aux , ispin, error)
    END SELECT

    CALL timestop(handle)

  END SUBROUTINE rtp_admm_calculate_dm

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

    INTEGER                                  :: handle
    LOGICAL                                  :: recalc_S
    TYPE(section_vals_type), POINTER         :: input, xc_section

    CALL timeset(routineN,handle)

    NULLIFY(xc_section)

    IF (.NOT.(ASSOCIATED(admm_env) )) THEN
      CALL admm_env_create(mos, mos_aux_fit, &
                           para_env, admm_env,&
                           error)
      CALL get_qs_env(qs_env, input=input, error=error)
      xc_section => section_vals_get_subs_vals(input,"DFT%XC",error=error)
      CALL create_admm_xc_section(qs_env, xc_section, admm_env, error)
      admm_env%method_id = admm_control%method_id
      admm_env%purification_method = admm_control%purification_method
      IF( admm_env%method_id == do_admm_block_density_matrix) THEN
         CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
              "Blocking not implemented with RTP."//&
              CPSourceFileRef,&
              only_ionode=.TRUE.)
      END IF
    END IF
    recalc_S=geometry_did_change.OR.(rtp%iter==0.AND.(rtp%istep==rtp%i_start))
    SELECT CASE(admm_env%purification_method)
    CASE(do_admm_purify_none)
       CALL rtp_fit_mo_coeffs_none(qs_env, admm_env, para_env, matrix_s_aux_fit, matrix_s_mixed, &
            mos, mos_aux_fit,rtp,rtp_coeff_aux_fit,recalc_S , error)
    CASE(do_admm_purify_cauchy)
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
!> \param admm_env The ADMM env
!> \param para_env The parallel env
!> \param matrix_s_aux_fit the overlap matrix of the auxiliary fitting basis set
!> \param matrix_s_mixed the mixed overlap matrix of the auxiliary fitting basis
!>        set and the orbital basis set
!> \param mos the MO's of the orbital basis set
!> \param mos_aux_fit the MO's of the auxiliary fitting basis set
!> \param geometry_did_change flag to indicate if the geomtry changed
!> \param error
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

    INTEGER                                  :: handle, ispin, istat, &
                                                nao_aux_fit, nao_orb, nmo, &
                                                nmo_mos, nspins
    LOGICAL                                  :: failure
    REAL(KIND=dp), DIMENSION(:), POINTER     :: occ_num, occ_num_aux
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mos_new
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_coeff_aux_fit
    TYPE(section_vals_type), POINTER         :: input, xc_section

    CALL timeset(routineN,handle)

    IF (.NOT.(ASSOCIATED(admm_env) )) THEN
      CALL admm_env_create(mos, mos_aux_fit, &
                           para_env, admm_env,&
                           error)
      CALL get_qs_env(qs_env, input=input, error=error)
      xc_section => section_vals_get_subs_vals(input,"DFT%XC",error=error)
      CALL create_admm_xc_section(qs_env, xc_section, admm_env, error)
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
      CALL cp_fm_gemm('N','N',nao_aux_fit,nao_orb,nao_aux_fit,&
                    1.0_dp,admm_env%S_inv,admm_env%Q,0.0_dp,&
                    admm_env%A,error)
   END IF

    ! *** Calculate the mo_coeffs for the fitting basis
   ALLOCATE(rtp_coeff_aux_fit(2*nspins), STAT=istat)
   CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
    DO ispin=1,nspins
      nmo = admm_env%nmo(ispin)
      IF( nmo == 0 ) CYCLE
      !! Lambda = C^(T)*B*C
      mos_new=>qs_env%rtp%mos_new
      CALL get_mo_set(mos(ispin)%mo_set,mo_coeff=mo_coeff, occupation_numbers=occ_num,nmo=nmo_mos)
      CALL get_mo_set(mos_aux_fit(ispin)%mo_set,mo_coeff=mo_coeff_aux_fit,&
           occupation_numbers=occ_num_aux)

      CALL cp_fm_create(matrix=rtp_coeff_aux_fit(2*ispin-1)%matrix,&
           matrix_struct=mo_coeff_aux_fit%matrix_struct,error=error)
      CALL cp_fm_create(matrix=rtp_coeff_aux_fit(2*ispin)%matrix,&
           matrix_struct=mo_coeff_aux_fit%matrix_struct,error=error)

      CALL cp_fm_gemm('N','N',nao_aux_fit,nmo,nao_orb,&
                      1.0_dp,admm_env%A,mos_new(2*ispin-1)%matrix,0.0_dp,&
                      rtp_coeff_aux_fit(2*ispin-1)%matrix,error)
      CALL cp_fm_gemm('N','N',nao_aux_fit,nmo,nao_orb,&
                      1.0_dp,admm_env%A,mos_new(2*ispin)%matrix,0.0_dp,&
                      rtp_coeff_aux_fit(2*ispin)%matrix,error)


      CALL cp_fm_to_fm(rtp_coeff_aux_fit(2*ispin-1)%matrix,mo_coeff_aux_fit,error)
    END DO

    CALL timestop(handle)

  END SUBROUTINE rtp_fit_mo_coeffs_none


  SUBROUTINE calculate_rtp_admm_density(mo_set, density_matrix_aux,&
                                   rtp_coeff_aux_fit,ispin, error)

    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_dbcsr_type), POINTER             :: density_matrix_aux
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: rtp_coeff_aux_fit
    INTEGER                                  :: ispin
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_rtp_admm_density', &
      routineP = moduleN//':'//routineN
    REAL(KIND=dp), PARAMETER                 :: one = 1.0_dp , zero = 0.0_dp

    INTEGER                                  :: handle, im_part, re_part
    TYPE(cp_dbcsr_type), POINTER             :: tmp_rho_im, tmp_rho_re

    CALL timeset(routineN,handle)

    ALLOCATE(tmp_rho_im, tmp_rho_re)
    CALL cp_dbcsr_init(tmp_rho_re,error=error)
    CALL cp_dbcsr_init(tmp_rho_im,error=error)
    CALL cp_dbcsr_copy(tmp_rho_re,density_matrix_aux,"tmp_tho_rts_re",error=error)
    CALL cp_dbcsr_copy(tmp_rho_im,density_matrix_aux,"tmp_tho_rts_im",error=error)

    re_part=2*ispin-1
    im_part=2*ispin
    CALL cp_fm_to_fm(rtp_coeff_aux_fit(im_part)%matrix,mo_set%mo_coeff,error)
    CALL calculate_density_matrix(mo_set,tmp_rho_im,error=error)
    CALL cp_fm_to_fm(rtp_coeff_aux_fit(re_part)%matrix,mo_set%mo_coeff,error)
    CALL calculate_density_matrix(mo_set,tmp_rho_re,error=error)

    CALL cp_dbcsr_add(density_matrix_aux,tmp_rho_re,&
         alpha_scalar=zero,beta_scalar=one,error=error)
    CALL cp_dbcsr_add(density_matrix_aux,tmp_rho_im,&
         alpha_scalar=one,beta_scalar=one,error=error)

    CALL cp_dbcsr_deallocate_matrix(tmp_rho_re,error)
    CALL cp_dbcsr_deallocate_matrix(tmp_rho_im,error)

    CALL timestop(handle)

  END SUBROUTINE calculate_rtp_admm_density

  SUBROUTINE rtp_purify_dm_cauchy(admm_env,rtp_coeff_aux_fit,mo_set,density_matrix,ispin,error)

    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: rtp_coeff_aux_fit
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_dbcsr_type), POINTER             :: density_matrix
    INTEGER                                  :: ispin
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'rtp_purify_dm_cauchy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, nao_aux_fit, &
                                                nao_orb, nmo, nspins
    REAL(KIND=dp)                            :: pole
    TYPE(cp_fm_type), POINTER                :: mo_coeff_aux_fit

    CALL timeset(routineN,handle)

    nao_aux_fit = admm_env%nao_aux_fit
    nao_orb = admm_env%nao_orb
    nmo = admm_env%nmo(ispin)

    nspins = SIZE(admm_env%P_to_be_purified)

    CALL get_mo_set(mo_set=mo_set, mo_coeff=mo_coeff_aux_fit)

    !! * For the time beeing, get the P to be purified from the mo_coeffs
    !! * This needs to be replaced with the a block modified P

    CALL cp_fm_gemm('N','T',nao_aux_fit,nao_aux_fit,nmo,&
                     1.0_dp,rtp_coeff_aux_fit(2*ispin-1)%matrix,&
                     rtp_coeff_aux_fit(2*ispin-1)%matrix,0.0_dp,&
                     admm_env%P_to_be_purified(ispin)%matrix,error)
    CALL cp_fm_gemm('N','T',nao_aux_fit,nao_aux_fit,nmo,&
                     1.0_dp,rtp_coeff_aux_fit(2*ispin)%matrix,&
                     rtp_coeff_aux_fit(2*ispin)%matrix,1.0_dp,&
                     admm_env%P_to_be_purified(ispin)%matrix,error)

    CALL cp_fm_to_fm(admm_env%S, admm_env%work_aux_aux, error=error)
    CALL cp_fm_to_fm(admm_env%P_to_be_purified(ispin)%matrix, admm_env%work_aux_aux2, error=error)

    CALL cp_fm_cholesky_decompose(admm_env%work_aux_aux,error=error)

    CALL cp_fm_cholesky_reduce(admm_env%work_aux_aux2, admm_env%work_aux_aux, itype=3, error=error)

    CALL cp_fm_syevd(admm_env%work_aux_aux2,admm_env%R_purify(ispin)%matrix,&
                     admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data,error=error)

    CALL cp_fm_cholesky_restore(admm_env%R_purify(ispin)%matrix, nao_aux_fit,admm_env%work_aux_aux, &
                                admm_env%work_aux_aux3,op="MULTIPLY",pos="LEFT", transa="T", error=error)

    CALL cp_fm_to_fm(admm_env%work_aux_aux3, admm_env%R_purify(ispin)%matrix, error=error)

    ! *** Construct Matrix M for Hadamard Product
    CALL cp_fm_set_all(admm_env%M_purify(ispin)%matrix,0.0_dp,error=error)
    pole = 0.0_dp
    DO i=1,nao_aux_fit
      pole = Heaviside(admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data(i)-0.5_dp)
      CALL cp_fm_set_element(admm_env%M_purify(ispin)%matrix,i,i,pole,error)
    END DO
    CALL cp_fm_upper_to_full(admm_env%M_purify(ispin)%matrix,admm_env%work_aux_aux,error=error)

    CALL copy_dbcsr_to_fm(density_matrix,admm_env%work_aux_aux3,error)
    CALL cp_fm_upper_to_full(admm_env%work_aux_aux3,admm_env%work_aux_aux,error=error)

    ! ** S^(-1)*R
    CALL cp_fm_gemm('N','N',nao_aux_fit,nao_aux_fit,nao_aux_fit,&
                    1.0_dp,admm_env%S_inv,admm_env%R_purify(ispin)%matrix,0.0_dp,&
                    admm_env%work_aux_aux,error)
    ! ** S^(-1)*R*M
    CALL cp_fm_gemm('N','N',nao_aux_fit,nao_aux_fit,nao_aux_fit,&
                    1.0_dp,admm_env%work_aux_aux,admm_env%M_purify(ispin)%matrix,0.0_dp,&
                    admm_env%work_aux_aux2,error)
    ! ** S^(-1)*R*M*R^T*S^(-1)
    CALL cp_fm_gemm('N','T',nao_aux_fit,nao_aux_fit,nao_aux_fit,&
                    1.0_dp,admm_env%work_aux_aux2,admm_env%work_aux_aux,0.0_dp,&
                    admm_env%work_aux_aux3,error)

    CALL copy_fm_to_dbcsr(admm_env%work_aux_aux3, density_matrix,keep_sparsity=.TRUE.,error=error)

    IF( nspins == 1 ) THEN
      CALL cp_dbcsr_scale(density_matrix, 2.0_dp, error=error)
    END IF

    CALL timestop(handle)

  END SUBROUTINE rtp_purify_dm_cauchy
  SUBROUTINE calculate_dm_mo_no_diag(admm_env,mo_set,density_matrix,ispin,error)
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_dbcsr_type), POINTER             :: density_matrix
    INTEGER                                  :: ispin
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_dm_mo_no_diag', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nao_aux_fit, nmo
    REAL(KIND=dp)                            :: alpha

    CALL timeset(routineN,handle)


    CALL cp_dbcsr_set(density_matrix,0.0_dp,error=error)
    nao_aux_fit = admm_env%nao_aux_fit
    nmo = admm_env%nmo(ispin)
    CALL cp_fm_to_fm(admm_env%C_hat(ispin)%matrix, admm_env%work_aux_nmo(ispin)%matrix, error=error)
    CALL cp_fm_column_scale(admm_env%work_aux_nmo(ispin)%matrix,mo_set%occupation_numbers(1:mo_set%homo))

    CALL cp_fm_gemm('N','N',nao_aux_fit,nmo,nmo,&
                    1.0_dp,admm_env%work_aux_nmo(ispin)%matrix,admm_env%lambda_inv(ispin)%matrix,0.0_dp,&
                    admm_env%work_aux_nmo2(ispin)%matrix,error)

    IF ( .NOT. mo_set%uniform_occupation ) THEN ! not all orbitals 1..homo are equally occupied
      alpha=1.0_dp
      CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=density_matrix,&
                              matrix_v=admm_env%C_hat(ispin)%matrix,&
                              matrix_g=admm_env%work_aux_nmo2(ispin)%matrix,&
                              ncol=mo_set%homo,&
                              alpha=alpha,error=error)
    ELSE
      alpha=1.0_dp
      CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=density_matrix,&
                              matrix_v=admm_env%C_hat(ispin)%matrix,&
                              matrix_g=admm_env%work_aux_nmo2(ispin)%matrix,&
                              ncol=mo_set%homo,&
                              alpha=alpha,error=error)
    ENDIF

    CALL timestop(handle)

  END SUBROUTINE calculate_dm_mo_no_diag
  FUNCTION Heaviside(x)
    REAL(KIND=dp), INTENT(IN)                :: x
    REAL(KIND=dp)                            :: Heaviside

    IF( x < 0.0_dp ) THEN
      Heaviside = 0.0_dp
    ELSE
      Heaviside = 1.0_dp
    END IF
  END FUNCTION Heaviside



  SUBROUTINE rtp_admm_merge_ks_matrix(ispin,rtp, admm_env, mo_set, mo_coeff, mo_coeff_aux_fit, &
                                             matrix_ks, matrix_ks_aux_fit, matrix_s, &
                                             matrix_p_aux_fit, matrix_p, error)
    INTEGER, INTENT(IN)                      :: ispin
    TYPE(rt_prop_type), POINTER              :: rtp
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_coeff_aux_fit
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_ks_aux_fit, &
                                                matrix_s, matrix_p_aux_fit, &
                                                matrix_p
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rtp_admm_merge_ks_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mos_new

    failure = .FALSE.

    CALL timeset(routineN,handle)

    SELECT CASE(admm_env%purification_method)
    CASE(do_admm_purify_none)
       CALL rt_merge_ks_matrix_none(ispin, admm_env, mo_set, mo_coeff,&
            mo_coeff_aux_fit, matrix_ks, matrix_ks_aux_fit, &
            matrix_s, matrix_p_aux_fit, error)
    CASE(do_admm_purify_cauchy)
       CALL get_rtp (rtp=rtp,mos_new=mos_new,error=error)
        CALL rt_merge_ks_matrix_cauchy(ispin,mos_new, admm_env, mo_set, mo_coeff,&
             mo_coeff_aux_fit, matrix_ks, matrix_ks_aux_fit, &
             matrix_s, matrix_p_aux_fit, matrix_p, error)
    END SELECT

    CALL timestop(handle)

  END SUBROUTINE rtp_admm_merge_ks_matrix


  SUBROUTINE rt_merge_ks_matrix_cauchy(ispin,mos_new, admm_env, mo_set, mo_coeff, mo_coeff_aux_fit, &
                                    matrix_ks, matrix_ks_aux_fit, matrix_s, matrix_p_aux_fit,&
                                    matrix_p, error)
    INTEGER, INTENT(IN)                      :: ispin
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mos_new
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_coeff_aux_fit
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_ks_aux_fit, &
                                                matrix_s, matrix_p_aux_fit, &
                                                matrix_p
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rt_merge_ks_matrix_cauchy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, nao_aux_fit, &
                                                nao_orb, nmo
    INTEGER, SAVE                            :: counter = 0
    LOGICAL                                  :: failure
    REAL(dp)                                 :: eig_diff, pole, tmp
    TYPE(cp_dbcsr_type), POINTER             :: matrix_k_tilde

    failure = .FALSE.

    CALL timeset(routineN,handle)

    counter = counter + 1
    nao_aux_fit = admm_env%nao_aux_fit
    nao_orb = admm_env%nao_orb
    nmo = admm_env%nmo(ispin)

    !** Get P from mo_coeffs, otherwise we have troubles with occupation numbers ...
    CALL cp_fm_gemm('N', 'T',  nao_orb, nao_orb, nmo,&
                    1.0_dp, mos_new(2*ispin-1)%matrix,&
                    mos_new(2*ispin-1)%matrix, 0.0_dp,&
                    admm_env%work_orb_orb,error)

    CALL cp_fm_gemm('N', 'T',  nao_orb, nao_orb, nmo,&
                    1.0_dp, mos_new(2*ispin)%matrix,&
                    mos_new(2*ispin)%matrix, 1.0_dp,&
                    admm_env%work_orb_orb,error)
    !! A*P
    CALL cp_fm_gemm('N', 'N',  nao_aux_fit, nao_orb, nao_orb,&
                    1.0_dp, admm_env%A, admm_env%work_orb_orb, 0.0_dp,&
                    admm_env%work_aux_orb2,error)
    !! A*P*A^T
    CALL cp_fm_gemm('N', 'T',  nao_aux_fit, nao_aux_fit, nao_orb,&
                    1.0_dp, admm_env%work_aux_orb2, admm_env%A, 0.0_dp,&
                    admm_env%P_to_be_purified(ispin)%matrix,error)


    CALL cp_fm_to_fm(admm_env%S, admm_env%work_aux_aux, error=error)
    CALL cp_fm_to_fm(admm_env%P_to_be_purified(ispin)%matrix, admm_env%work_aux_aux2, error=error)

    CALL cp_fm_cholesky_decompose(admm_env%work_aux_aux,error=error)

    CALL cp_fm_cholesky_reduce(admm_env%work_aux_aux2, admm_env%work_aux_aux, itype=3, error=error)

    CALL cp_fm_syevd(admm_env%work_aux_aux2,admm_env%R_purify(ispin)%matrix,&
                     admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data,error=error)

    CALL cp_fm_cholesky_restore(admm_env%R_purify(ispin)%matrix, nao_aux_fit,admm_env%work_aux_aux, &
                                admm_env%work_aux_aux3,op="MULTIPLY",pos="LEFT", transa="T", error=error)

    CALL cp_fm_to_fm(admm_env%work_aux_aux3, admm_env%R_purify(ispin)%matrix, error=error)

    ! *** Construct Matrix M for Hadamard Product
    pole = 0.0_dp
    DO i=1,nao_aux_fit
      DO j=i,nao_aux_fit
        eig_diff = ( admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data(i) -&
                     admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data(j) )
        ! *** two eigenvalues could be the degenerated. In that case use 2nd order formula for the poles
        IF( ABS(eig_diff) == 0.0_dp ) THEN
          pole = delta(admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data(i)-0.5_dp)
          CALL cp_fm_set_element(admm_env%M_purify(ispin)%matrix,i,j,pole,error)
        ELSE
          pole = 1.0_dp/(admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data(i)-&
                         admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data(j))
          tmp = Heaviside(admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data(i)-0.5_dp)
          tmp = tmp - Heaviside(admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data(j)-0.5_dp)
          pole = tmp*pole
          CALL cp_fm_set_element(admm_env%M_purify(ispin)%matrix,i,j,pole,error)
        END IF
      END DO
    END DO
    CALL cp_fm_upper_to_full(admm_env%M_purify(ispin)%matrix,admm_env%work_aux_aux,error=error)

    CALL copy_dbcsr_to_fm(matrix_ks_aux_fit(ispin)%matrix,admm_env%K(ispin)%matrix,error)
    CALL cp_fm_upper_to_full(admm_env%K(ispin)%matrix,admm_env%work_aux_aux,error=error)


    !! S^(-1)*R
    CALL cp_fm_gemm('N', 'N',  nao_aux_fit, nao_aux_fit, nao_aux_fit,&
                    1.0_dp,admm_env%S_inv,admm_env%R_purify(ispin)%matrix,0.0_dp,&
                    admm_env%work_aux_aux,error)
    !! K*S^(-1)*R
    CALL cp_fm_gemm('N', 'N',  nao_aux_fit, nao_aux_fit, nao_aux_fit,&
                    1.0_dp,admm_env%K(ispin)%matrix,admm_env%work_aux_aux,0.0_dp,&
                    admm_env%work_aux_aux2,error)
    !! R^T*S^(-1)*K*S^(-1)*R
    CALL cp_fm_gemm('T', 'N',  nao_aux_fit, nao_aux_fit, nao_aux_fit,&
                    1.0_dp,admm_env%work_aux_aux,admm_env%work_aux_aux2,0.0_dp,&
                    admm_env%work_aux_aux3,error)
    !! R^T*S^(-1)*K*S^(-1)*R x M
    CALL cp_fm_schur_product(admm_env%work_aux_aux3, admm_env%M_purify(ispin)%matrix,&
                             admm_env%work_aux_aux,error)

    !! R^T*A
    CALL cp_fm_gemm('T', 'N',  nao_aux_fit, nao_orb, nao_aux_fit,&
                    1.0_dp, admm_env%R_purify(ispin)%matrix, admm_env%A, 0.0_dp,&
                    admm_env%work_aux_orb,error)

    !! (R^T*S^(-1)*K*S^(-1)*R x M) * R^T*A
    CALL cp_fm_gemm('N', 'N',  nao_aux_fit, nao_orb, nao_aux_fit,&
                    1.0_dp, admm_env%work_aux_aux, admm_env%work_aux_orb, 0.0_dp,&
                    admm_env%work_aux_orb2,error)
    !! A^T*R*(R^T*S^(-1)*K*S^(-1)*R x M) * R^T*A
    CALL cp_fm_gemm('T', 'N',  nao_orb, nao_orb, nao_aux_fit,&
                    1.0_dp, admm_env%work_aux_orb, admm_env%work_aux_orb2, 0.0_dp,&
                    admm_env%work_orb_orb,error)


    NULLIFY(matrix_k_tilde)
    ALLOCATE(matrix_k_tilde)
    CALL cp_dbcsr_init (matrix_k_tilde, error)
    CALL cp_dbcsr_create(matrix_k_tilde, 'MATRIX K_tilde', &
         cp_dbcsr_distribution(matrix_ks(ispin)%matrix), dbcsr_type_symmetric, cp_dbcsr_row_block_sizes(matrix_ks(ispin)%matrix),&
         cp_dbcsr_col_block_sizes(matrix_ks(ispin)%matrix), cp_dbcsr_get_num_blocks(matrix_ks(ispin)%matrix),&
         cp_dbcsr_get_data_size( matrix_ks(ispin)%matrix),&
         cp_dbcsr_get_data_type(matrix_ks(ispin)%matrix), error=error)

    CALL cp_fm_to_fm(admm_env%work_orb_orb, admm_env%ks_to_be_merged(ispin)%matrix, error=error)

    CALL cp_dbcsr_copy(matrix_k_tilde, matrix_ks(ispin)%matrix, error=error)
    CALL cp_dbcsr_set(matrix_k_tilde, 0.0_dp, error)
    CALL copy_fm_to_dbcsr(admm_env%work_orb_orb, matrix_k_tilde, keep_sparsity=.TRUE., error=error)

    CALL cp_dbcsr_add(matrix_ks(ispin)%matrix, matrix_k_tilde, 1.0_dp, 1.0_dp, error)

    CALL cp_dbcsr_deallocate_matrix(matrix_k_tilde,error)

    CALL timestop(handle)

  END SUBROUTINE rt_merge_ks_matrix_cauchy





  SUBROUTINE rt_merge_ks_matrix_none(ispin, admm_env, mo_set, mo_coeff, mo_coeff_aux_fit, &
                                  matrix_ks, matrix_ks_aux_fit, matrix_s, matrix_p_aux_fit, error)
    INTEGER, INTENT(IN)                      :: ispin
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_coeff_aux_fit
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_ks_aux_fit, &
                                                matrix_s, matrix_p_aux_fit
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rt_merge_ks_matrix_none', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nao_aux_fit, nao_orb, &
                                                nmo
    INTEGER, SAVE                            :: counter = 0
    LOGICAL                                  :: failure
    TYPE(cp_dbcsr_type), POINTER             :: matrix_k_tilde

    failure = .FALSE.

    CALL timeset(routineN,handle)

    counter = counter + 1
    nao_aux_fit = admm_env%nao_aux_fit
    nao_orb = admm_env%nao_orb
    nmo = admm_env%nmo(ispin)

    CALL copy_dbcsr_to_fm(matrix_ks_aux_fit(ispin)%matrix,admm_env%K(ispin)%matrix,error)
    CALL cp_fm_upper_to_full(admm_env%K(ispin)%matrix,admm_env%work_aux_aux,error=error)

    !! K*A
    CALL cp_fm_gemm('N', 'N',  nao_aux_fit, nao_orb, nao_aux_fit,&
                    1.0_dp,admm_env%K(ispin)%matrix,admm_env%A,0.0_dp,&
                    admm_env%work_aux_orb,error)
    !! A^T*K*A
    CALL cp_fm_gemm('T', 'N',  nao_orb, nao_orb, nao_aux_fit,&
                    1.0_dp,admm_env%A,admm_env%work_aux_orb,0.0_dp,&
                    admm_env%work_orb_orb,error)


    NULLIFY(matrix_k_tilde)
    ALLOCATE(matrix_k_tilde)
    CALL cp_dbcsr_init (matrix_k_tilde, error)
    CALL cp_dbcsr_create(matrix_k_tilde, 'MATRIX K_tilde', &
         cp_dbcsr_distribution(matrix_ks(ispin)%matrix), dbcsr_type_symmetric, cp_dbcsr_row_block_sizes(matrix_ks(ispin)%matrix),&
         cp_dbcsr_col_block_sizes(matrix_ks(ispin)%matrix), &
         cp_dbcsr_get_num_blocks(matrix_ks(ispin)%matrix), cp_dbcsr_get_data_size(matrix_ks(ispin)%matrix),&
         cp_dbcsr_get_data_type(matrix_ks(ispin)%matrix), error=error)

    CALL cp_dbcsr_copy(matrix_k_tilde, matrix_ks(ispin)%matrix, error=error)
    CALL cp_dbcsr_set(matrix_k_tilde, 0.0_dp, error)
    CALL copy_fm_to_dbcsr(admm_env%work_orb_orb, matrix_k_tilde, keep_sparsity=.TRUE.,&
         error=error)

    CALL cp_dbcsr_add(matrix_ks(ispin)%matrix, matrix_k_tilde, 1.0_dp, 1.0_dp, error)

    CALL cp_dbcsr_deallocate_matrix(matrix_k_tilde,error)

    CALL timestop(handle)

  END SUBROUTINE rt_merge_ks_matrix_none

  FUNCTION delta(x)
    REAL(KIND=dp), INTENT(IN)                :: x
    REAL(KIND=dp)                            :: delta

    IF( x == 0.0_dp) THEN
      delta = 1.0_dp
    ELSE
      delta = 0.0_dp
    END IF

  END FUNCTION delta
END MODULE rtp_admm_methods
