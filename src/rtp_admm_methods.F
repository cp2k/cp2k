
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
  USE admm_methods,                    ONLY: create_admm_xc_section
  USE admm_types,                      ONLY: admm_env_create,&
                                             admm_type
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE cp_control_types,                ONLY: admm_control_type,&
                                             dft_control_type
  USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_add, cp_dbcsr_col_block_sizes, cp_dbcsr_copy, &
       cp_dbcsr_create, cp_dbcsr_deallocate_matrix, cp_dbcsr_desymmetrize, &
       cp_dbcsr_distribution, cp_dbcsr_get_data_size, cp_dbcsr_get_data_type, &
       cp_dbcsr_get_matrix_type, cp_dbcsr_get_num_blocks, cp_dbcsr_init, &
       cp_dbcsr_p_type, cp_dbcsr_release, cp_dbcsr_row_block_sizes, &
       cp_dbcsr_set, cp_dbcsr_type, dbcsr_type_no_symmetry
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                             copy_fm_to_dbcsr,&
                                             cp_dbcsr_plus_fm_fm_t,&
                                             cp_dbcsr_sm_fm_multiply
  USE cp_fm_basic_linalg,              ONLY: cp_fm_upper_to_full
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                             cp_fm_cholesky_invert
  USE cp_fm_struct,                    ONLY: cp_fm_struct_type
  USE cp_fm_types,                     ONLY: cp_fm_create,&
                                             cp_fm_get_info,&
                                             cp_fm_p_type,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_fm_vect,                      ONLY: cp_fm_vect_dealloc
  USE cp_gemm_interface,               ONLY: cp_gemm
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_constants,                 ONLY: do_admm_basis_projection,&
                                             do_admm_purify_none,&
                                             use_aux_fit_basis_set,&
                                             use_orb_basis_set
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE mathconstants,                   ONLY: one,&
                                             zero
  USE particle_types,                  ONLY: particle_type
  USE qs_collocate_density,            ONLY: calculate_rho_elec
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_force_types,                  ONLY: add_qs_force,&
                                             qs_force_type
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type,&
                                             mo_set_type
  USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
  USE qs_overlap,                      ONLY: build_overlap_force
  USE qs_rho_types,                    ONLY: qs_rho_type
  USE rt_propagation_types,            ONLY: get_rtp,&
                                             rt_prop_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  ! *** Public subroutines ***
  PUBLIC :: rtp_admm_density_matrix,rtp_admm_merge_ks_matrix, rt_admm_force

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rtp_admm_methods'

CONTAINS

! *****************************************************************************
!> \brief  Compute the ADMM density matrix in case of rtp (complex MO's)
!>
!> \param qs_env ...
!> \param ks_env ...
!> \param dft_control ...
!> \param input ...
!> \param error ...
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
    LOGICAL                                  :: do_hfx, failure, &
                                                s_mstruct_changed
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_p, matrix_p_aux_fit, &
                                                matrix_s_aux_fit, &
                                                matrix_s_aux_fit_vs_orb
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: rtp_coeff_aux_fit
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos, mos_aux_fit
    TYPE(qs_rho_type), POINTER               :: rho, rho_aux_fit
    TYPE(rt_prop_type), POINTER              :: rtp
    TYPE(section_vals_type), POINTER         :: hfx_sections

    failure=.FALSE.
    CALL timeset(routineN,handle)
    NULLIFY(admm_env, hfx_sections, matrix_p, matrix_p_aux_fit, mos, &
         mos_aux_fit, para_env,matrix_s_aux_fit,matrix_s_aux_fit_vs_orb, rho)

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
    CALL get_qs_env(qs_env,&
                    mos_aux_fit=mos_aux_fit,&
                    rtp=rtp,&
                    mos=mos,&
                    admm_env=admm_env,&
                    para_env=para_env,&
                    matrix_s_aux_fit=matrix_s_aux_fit,&
                    matrix_s_aux_fit_vs_orb=matrix_s_aux_fit_vs_orb,&
                    s_mstruct_changed=s_mstruct_changed,&
                    rho=rho,&
                    rho_aux_fit=rho_aux_fit,&
                    error=error)

    CALL section_vals_get(hfx_sections,n_repetition=n_rep_hf,error=error)
    IF ( n_rep_hf > 1 ) THEN
       CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
                      "ADMM can handle only one HF section."//&
                       CPSourceFileRef,&
                       only_ionode=.TRUE.)
    END IF
    CALL get_rtp(rtp=rtp,admm_mos=rtp_coeff_aux_fit,error=error)
    CALL rtp_admm_fit_mo_coeffs(qs_env, admm_env, dft_control%admm_control, para_env, &
                            matrix_s_aux_fit, matrix_s_aux_fit_vs_orb,&
                            mos, mos_aux_fit,rtp,rtp_coeff_aux_fit,&
                            s_mstruct_changed, error=error)

    DO ispin=1,nspins
       matrix_p_aux_fit => rho_aux_fit%rho_ao
       matrix_p => rho%rho_ao
       CALL rtp_admm_calculate_dm(admm_env,rtp_coeff_aux_fit,mos_aux_fit(ispin)%mo_set,&
                                       matrix_p_aux_fit(ispin)%matrix,&
                                       rtp%rtp_hfx%matrix_p_admm_im(ispin)%matrix,&
                                       ispin,nspins,error)
       CALL calculate_rho_elec(matrix_p=rho_aux_fit%rho_ao(ispin)%matrix,&
            rho=rho_aux_fit%rho_r(ispin),&
            rho_gspace=rho_aux_fit%rho_g(ispin),&
            total_rho=rho_aux_fit%tot_rho_r(ispin),&
            qs_env=qs_env,soft_valid=.FALSE.,&
            basis_set_id=use_aux_fit_basis_set, error=error)
    END DO
    CALL set_qs_env(qs_env,admm_env=admm_env,error=error)
    rho_aux_fit%rho_r_valid=.TRUE.
    rho_aux_fit%rho_g_valid=.TRUE.

    CALL timestop(handle)

  END SUBROUTINE rtp_admm_density_matrix

! *****************************************************************************
!> \brief ...
!> \param admm_env ...
!> \param rtp_coeff_aux_fit ...
!> \param mo_set ...
!> \param density_matrix_aux ...
!> \param density_matrix_aux_im ...
!> \param ispin ...
!> \param nspins ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE rtp_admm_calculate_dm(admm_env,rtp_coeff_aux_fit,mo_set, density_matrix_aux,&
                                   density_matrix_aux_im, ispin, nspins, error)
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: rtp_coeff_aux_fit
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_dbcsr_type), POINTER             :: density_matrix_aux, &
                                                density_matrix_aux_im
    INTEGER, INTENT(in)                      :: ispin, nspins
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
      CALL create_admm_xc_section(qs_env, xc_section, admm_env, error)

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
!> \param x ...
! *****************************************************************************
  FUNCTION Heaviside(x)
    REAL(KIND=dp), INTENT(IN)                :: x
    REAL(KIND=dp)                            :: Heaviside

    IF( x < 0.0_dp ) THEN
      Heaviside = 0.0_dp
    ELSE
      Heaviside = 1.0_dp
    END IF
  END FUNCTION Heaviside



! *****************************************************************************
!> \brief ...
!> \param ispin ...
!> \param rtp ...
!> \param admm_env ...
!> \param mo_set ...
!> \param mo_coeff ...
!> \param mo_coeff_aux_fit ...
!> \param matrix_ks ...
!> \param matrix_ks_aux_fit ...
!> \param matrix_s ...
!> \param matrix_p_aux_fit ...
!> \param matrix_p ...
!> \param error ...
! *****************************************************************************
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

    failure = .FALSE.

    CALL timeset(routineN,handle)

    SELECT CASE(admm_env%purification_method)
    CASE(do_admm_purify_none)
       CALL rt_merge_ks_matrix_none(ispin, admm_env, &
            matrix_ks, matrix_ks_aux_fit, &
            error)
       CALL rt_merge_ks_matrix_none(ispin, admm_env, &
            rtp%rtp_hfx%matrix_ks_im, rtp%rtp_hfx%matrix_ks_admm_im, &
            error)
    CASE DEFAULT
       CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,routineP,&
                      "only purification NONE possible with RTP/EMD at the moment"//&
                       CPSourceFileRef,&
                       only_ionode=.TRUE.)
    END SELECT

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

! *****************************************************************************
!> \brief ...
!> \param x ...
! *****************************************************************************
  FUNCTION delta(x)
    REAL(KIND=dp), INTENT(IN)                :: x
    REAL(KIND=dp)                            :: delta

    IF( x == 0.0_dp) THEN
      delta = 1.0_dp
    ELSE
      delta = 0.0_dp
    END IF

  END FUNCTION delta
 
! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE rt_admm_force(qs_env,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rt_admm_force', &
      routineP = moduleN//':'//routineN

    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: KS_aux_re, matrix_s_aux_fit, &
                                                matrix_s_aux_fit_vs_orb
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mos, mos_admm
    TYPE(rt_prop_type), POINTER              :: rtp

    CALL get_qs_env(qs_env=qs_env, admm_env=admm_env, matrix_ks_aux_fit=KS_aux_re, matrix_s_aux_fit=matrix_s_aux_fit,&
                     matrix_s_aux_fit_vs_orb= matrix_s_aux_fit_vs_orb , rtp=rtp, error=error)

    CALL get_rtp(rtp=rtp,mos_new=mos,admm_mos=mos_admm,error=error)

    ! currently only none option
    CALL rt_admm_forces_none(qs_env,admm_env,KS_aux_re,rtp%rtp_hfx%matrix_ks_admm_im, &
                             matrix_s_aux_fit, matrix_s_aux_fit_vs_orb , mos_admm,mos,error)

  END SUBROUTINE rt_admm_force

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param admm_env ...
!> \param KS_aux_re ...
!> \param KS_aux_im ...
!> \param matrix_s_aux_fit ...
!> \param matrix_s_aux_fit_vs_orb ...
!> \param mos_admm ...
!> \param mos ...
!> \param error ...
! *****************************************************************************
  SUBROUTINE rt_admm_forces_none(qs_env,admm_env,KS_aux_re,KS_aux_im,matrix_s_aux_fit, matrix_s_aux_fit_vs_orb,mos_admm,mos,error)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: KS_aux_re, KS_aux_im, &
                                                matrix_s_aux_fit, &
                                                matrix_s_aux_fit_vs_orb
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mos_admm, mos
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'rt_admm_forces_none', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: im, ispin, istat, nao, natom, &
                                                naux, nmo, re
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: admm_force
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_type), POINTER             :: matrix_w_q, matrix_w_s
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: tmp_aux_aux, tmp_aux_mo, &
                                                tmp_aux_mo1, tmp_aux_nao
    TYPE(cp_fm_struct_type), POINTER         :: mstruct
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_aux_fit_asymm, &
                                                sab_aux_fit_vs_orb
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force

    NULLIFY(sab_aux_fit_asymm, sab_aux_fit_vs_orb)
    failure = .FALSE.
!   CALL cp_fm_create(tmp_aux_aux,admm_env%fm_struct_tmp,name="fm matrix",error=error)

    CALL get_qs_env(qs_env,&
                    sab_aux_fit_asymm=sab_aux_fit_asymm,&
                    sab_aux_fit_vs_orb=sab_aux_fit_vs_orb,&
                    error=error)

    ALLOCATE(matrix_w_s)
    CALL cp_dbcsr_init (matrix_w_s, error)
    CALL cp_dbcsr_create(matrix_w_s, 'W MATRIX AUX S', &
         cp_dbcsr_distribution(matrix_s_aux_fit(1)%matrix), dbcsr_type_no_symmetry, &
         cp_dbcsr_row_block_sizes(matrix_s_aux_fit(1)%matrix),&
         cp_dbcsr_col_block_sizes(matrix_s_aux_fit(1)%matrix), &
         cp_dbcsr_get_num_blocks(matrix_s_aux_fit(1)%matrix), &
         cp_dbcsr_get_data_size(matrix_s_aux_fit(1)%matrix),&
         cp_dbcsr_get_data_type(matrix_s_aux_fit(1)%matrix), &
         error=error)
    CALL cp_dbcsr_alloc_block_from_nbl(matrix_w_s,sab_aux_fit_asymm,error=error)

    ALLOCATE(matrix_w_q)
    CALL cp_dbcsr_init(matrix_w_q, error=error)
    CALL cp_dbcsr_copy(matrix_w_q,matrix_s_aux_fit_vs_orb(1)%matrix,&
                    "W MATRIX AUX Q",error=error)

    DO ispin=1,SIZE(KS_aux_re)
       re=2*ispin-1; im=2*ispin
       naux=admm_env%nao_aux_fit; nmo=admm_env%nmo(ispin); nao=admm_env%nao_orb
       
       ALLOCATE(tmp_aux_aux(2),tmp_aux_nao(2),tmp_aux_mo(2),tmp_aux_mo1(2))
       CALL cp_fm_create(tmp_aux_aux(1)%matrix,admm_env%work_aux_aux%matrix_struct,name="taa",error=error)
       CALL cp_fm_create(tmp_aux_aux(2)%matrix,admm_env%work_aux_aux%matrix_struct,name="taa",error=error)
       CALL cp_fm_create(tmp_aux_nao(1)%matrix,admm_env%work_aux_orb%matrix_struct,name="tao",error=error)
       CALL cp_fm_create(tmp_aux_nao(2)%matrix,admm_env%work_aux_orb%matrix_struct,name="tao",error=error)
       mstruct => admm_env%work_aux_nmo(ispin)%matrix%matrix_struct
       CALL cp_fm_create(tmp_aux_mo(1)%matrix,mstruct,name="tam",error=error)
       CALL cp_fm_create(tmp_aux_mo(2)%matrix,mstruct,name="tam",error=error)
       CALL cp_fm_create(tmp_aux_mo1(1)%matrix,mstruct,name="tam",error=error)
       CALL cp_fm_create(tmp_aux_mo1(2)%matrix,mstruct,name="tam",error=error)

! First calculate H=KS_aux*C~, real part ends on work_aux_aux2, imaginary part ends at work_aux_aux3
       CALL cp_dbcsr_sm_fm_multiply(KS_aux_re(ispin)%matrix,mos_admm(re)%matrix,tmp_aux_mo(re)%matrix,nmo,4.0_dp,0.0_dp,error)
       CALL cp_dbcsr_sm_fm_multiply(KS_aux_re(ispin)%matrix,mos_admm(im)%matrix,tmp_aux_mo(im)%matrix,nmo,4.0_dp,0.0_dp,error)
       CALL cp_dbcsr_sm_fm_multiply(KS_aux_im(ispin)%matrix,mos_admm(im)%matrix,tmp_aux_mo(re)%matrix,nmo,-4.0_dp,1.0_dp,error)
       CALL cp_dbcsr_sm_fm_multiply(KS_aux_im(ispin)%matrix,mos_admm(re)%matrix,tmp_aux_mo(im)%matrix,nmo,4.0_dp,1.0_dp,error)
 
! Next step compute S-1*H
       CALL cp_gemm('N','N',naux,nmo,naux,1.0_dp,admm_env%S_inv,tmp_aux_mo(re)%matrix,0.0_dp,tmp_aux_mo1(re)%matrix,error)
       CALL cp_gemm('N','N',naux,nmo,naux,1.0_dp,admm_env%S_inv,tmp_aux_mo(im)%matrix,0.0_dp,tmp_aux_mo1(im)%matrix,error)

! Here we go on with Ws=S-1*H * C^H (take care of sign of the imaginary part!!!)

       CALL cp_gemm("N","T",naux,nao,nmo,-1.0_dp,tmp_aux_mo1(re)%matrix, mos(re)%matrix, 0.0_dp,&
                       tmp_aux_nao(re)%matrix, error)
       CALL cp_gemm("N","T",naux,nao,nmo,-1.0_dp,tmp_aux_mo1(im)%matrix, mos(im)%matrix, 1.0_dp,&
                       tmp_aux_nao(re)%matrix, error)
       CALL cp_gemm("N","T",naux,nao,nmo,1.0_dp,tmp_aux_mo1(re)%matrix, mos(im)%matrix, 0.0_dp,&
                       tmp_aux_nao(im)%matrix, error)
       CALL cp_gemm("N","T",naux,nao,nmo,-1.0_dp,tmp_aux_mo1(im)%matrix, mos(re)%matrix, 1.0_dp,&
                       tmp_aux_nao(im)%matrix, error)

! Let's do the final bit  Wq=S-1*H * C^H * A^T
       CALL cp_gemm('N','T',naux,naux,nao,-1.0_dp,tmp_aux_nao(re)%matrix,admm_env%A,0.0_dp,tmp_aux_aux(re)%matrix,error)     
       CALL cp_gemm('N','T',naux,naux,nao,-1.0_dp,tmp_aux_nao(im)%matrix,admm_env%A,0.0_dp,tmp_aux_aux(im)%matrix,error)     

       ! *** copy to sparse matrix
       CALL copy_fm_to_dbcsr(tmp_aux_nao(re)%matrix, matrix_w_q,keep_sparsity=.TRUE.,&
            error=error)


       ! *** copy to sparse matrix
       CALL copy_fm_to_dbcsr(tmp_aux_aux(re)%matrix, matrix_w_s,keep_sparsity=.TRUE.,&
            error=error)

! *** This can be done in one call w_total = w_alpha + w_beta
       ! allocate force vector
       CALL get_qs_env(qs_env=qs_env,natom=natom,error=error)
       ALLOCATE(admm_force(3,natom),STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       admm_force = 0.0_dp
       CALL build_overlap_force(qs_env, admm_force,&
            basis_set_id_a=use_aux_fit_basis_set, basis_set_id_b=use_aux_fit_basis_set, &
            sab_nl=sab_aux_fit_asymm, matrix_p=matrix_w_s, error=error)
       CALL build_overlap_force(qs_env, admm_force,&
            basis_set_id_a=use_aux_fit_basis_set, basis_set_id_b=use_orb_basis_set, &
            sab_nl=sab_aux_fit_vs_orb, matrix_p=matrix_w_q, error=error)
       ! add forces
       CALL get_qs_env(qs_env=qs_env,atomic_kind_set=atomic_kind_set,&
                       force=force,error=error)
       CALL add_qs_force(admm_force, force, "overlap_admm", atomic_kind_set, error)
       DEALLOCATE(admm_force,STAT=istat)
       CPPostcondition(istat==0,cp_failure_level,routineP,error,failure)
       
       ! *** Deallocated weighted density matrices
       CALL cp_dbcsr_deallocate_matrix(matrix_w_s,error)
       CALL cp_dbcsr_deallocate_matrix(matrix_w_q,error)
       CALL cp_fm_vect_dealloc(tmp_aux_aux,error)
       CALL cp_fm_vect_dealloc(tmp_aux_nao,error)
       CALL cp_fm_vect_dealloc(tmp_aux_mo,error)
       CALL cp_fm_vect_dealloc(tmp_aux_mo1,error)
     END DO
    
  END SUBROUTINE rt_admm_forces_none


END MODULE rtp_admm_methods
