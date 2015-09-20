!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Contains ADMM methods which require molecular orbitals
!> \par History
!>      04.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
MODULE admm_methods
  USE admm_types,                      ONLY: admm_type
  USE atomic_kind_types,               ONLY: atomic_kind_type
  USE bibliography,                    ONLY: Merlot2014,&
                                             cite_reference
  USE cp_control_types,                ONLY: dft_control_type
  USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_add, cp_dbcsr_col_block_sizes, cp_dbcsr_copy, &
       cp_dbcsr_create, cp_dbcsr_deallocate_matrix, cp_dbcsr_desymmetrize, &
       cp_dbcsr_distribution, cp_dbcsr_get_block_p, cp_dbcsr_get_data_size, &
       cp_dbcsr_get_data_type, cp_dbcsr_init, cp_dbcsr_iterator, &
       cp_dbcsr_iterator_blocks_left, cp_dbcsr_iterator_next_block, &
       cp_dbcsr_iterator_start, cp_dbcsr_iterator_stop, cp_dbcsr_p_type, &
       cp_dbcsr_row_block_sizes, cp_dbcsr_scale, cp_dbcsr_set, &
       cp_dbcsr_trace, cp_dbcsr_type, dbcsr_type_no_symmetry, &
       dbcsr_type_symmetric
  USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                             copy_fm_to_dbcsr,&
                                             cp_dbcsr_plus_fm_fm_t
  USE cp_dbcsr_output,                 ONLY: cp_dbcsr_write_sparse_matrix
  USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                             cp_fm_scale,&
                                             cp_fm_scale_and_add,&
                                             cp_fm_schur_product,&
                                             cp_fm_upper_to_full
  USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                             cp_fm_cholesky_invert,&
                                             cp_fm_cholesky_reduce,&
                                             cp_fm_cholesky_restore
  USE cp_fm_diag,                      ONLY: cp_fm_syevd
  USE cp_fm_types,                     ONLY: cp_fm_p_type,&
                                             cp_fm_set_all,&
                                             cp_fm_set_element,&
                                             cp_fm_to_fm,&
                                             cp_fm_type
  USE cp_gemm_interface,               ONLY: cp_gemm
  USE cp_output_handling,              ONLY: cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE input_constants,                 ONLY: do_admm_exch_scaling_merlot,&
                                             do_admm_exch_scaling_none,&
                                             do_admm_purify_cauchy,&
                                             do_admm_purify_cauchy_subspace,&
                                             do_admm_purify_mo_diag,&
                                             do_admm_purify_mo_no_diag,&
                                             do_admm_purify_none
  USE input_section_types,             ONLY: section_vals_type
  USE kinds,                           ONLY: dp
  USE pw_types,                        ONLY: pw_p_type
  USE qs_collocate_density,            ONLY: calculate_rho_elec
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE qs_force_types,                  ONLY: add_qs_force,&
                                             qs_force_type
  USE qs_ks_types,                     ONLY: qs_ks_env_type
  USE qs_mo_types,                     ONLY: get_mo_set,&
                                             mo_set_p_type,&
                                             mo_set_type
  USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
  USE qs_overlap,                      ONLY: build_overlap_force
  USE qs_rho_types,                    ONLY: qs_rho_get,&
                                             qs_rho_set,&
                                             qs_rho_type
  USE qs_vxc,                          ONLY: qs_vxc_create
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: admm_mo_calc_rho_aux,&
            admm_mo_merge_ks_matrix,&
            admm_mo_merge_derivs,&
            calc_mixed_overlap_force,&
            calc_aux_mo_derivs_none,&
            scale_dm

  PUBLIC :: admm_fit_mo_coeffs

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'admm_methods'

  CONTAINS

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
! *****************************************************************************
  SUBROUTINE admm_mo_calc_rho_aux(qs_env)
    TYPE(qs_environment_type), POINTER       :: qs_env

    CHARACTER(len=*), PARAMETER :: routineN = 'admm_mo_calc_rho_aux', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin
    LOGICAL                                  :: s_mstruct_changed
    REAL(KIND=dp), DIMENSION(:), POINTER     :: tot_rho_r_aux
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s, matrix_s_aux_fit, &
                                                matrix_s_aux_fit_vs_orb, &
                                                rho_ao, rho_ao_aux
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos, mos_aux_fit
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho_g_aux, rho_r_aux
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_rho_type), POINTER               :: rho, rho_aux_fit

    CALL timeset(routineN,handle)

    NULLIFY(ks_env, admm_env, mos, mos_aux_fit, matrix_s_aux_fit,&
            matrix_s_aux_fit_vs_orb,matrix_s, rho, rho_aux_fit, para_env)
    NULLIFY(rho_g_aux, rho_r_aux, rho_ao, rho_ao_aux, tot_rho_r_aux)

    CALL get_qs_env(qs_env,&
                    ks_env=ks_env,&
                    admm_env=admm_env,&
                    dft_control=dft_control,&
                    mos_aux_fit=mos_aux_fit,&
                    mos=mos,&
                    matrix_s_aux_fit=matrix_s_aux_fit,&
                    matrix_s_aux_fit_vs_orb=matrix_s_aux_fit_vs_orb,&
                    matrix_s=matrix_s,&
                    para_env=para_env,&
                    s_mstruct_changed=s_mstruct_changed,&
                    rho=rho,&
                    rho_aux_fit=rho_aux_fit)

    CALL qs_rho_get(rho, rho_ao=rho_ao)
    CALL qs_rho_get(rho_aux_fit,&
                    rho_ao=rho_ao_aux,&
                    rho_g=rho_g_aux,&
                    rho_r=rho_r_aux,&
                    tot_rho_r=tot_rho_r_aux)

    ! convert mos from full to dbcsr matrices
    DO ispin=1, dft_control%nspins
       IF(mos(ispin)%mo_set%use_mo_coeff_b) THEN
          CALL copy_dbcsr_to_fm(mos(ispin)%mo_set%mo_coeff_b,mos(ispin)%mo_set%mo_coeff)
       ENDIF
    ENDDO

    ! fit mo coeffcients
    CALL admm_fit_mo_coeffs(admm_env, matrix_s_aux_fit, matrix_s_aux_fit_vs_orb,&
                            mos, mos_aux_fit, s_mstruct_changed)

    DO ispin=1,dft_control%nspins
       IF(admm_env%block_dm) THEN
          CALL blockify_density_matrix(admm_env,&
                                       density_matrix=rho_ao(ispin)%matrix,&
                                       density_matrix_aux=rho_ao_aux(ispin)%matrix,&
                                       ispin=ispin,&
                                       nspins=dft_control%nspins)
    
       ELSE

         ! Here, the auxiliary DM gets calculated and is written into rho_aux_fit%... 
          CALL calculate_dm_mo_no_diag(admm_env,&
                                       mo_set=mos(ispin)%mo_set,&
                                       overlap_matrix=matrix_s_aux_fit(1)%matrix,&
                                       density_matrix=rho_ao_aux(ispin)%matrix,&
                                       overlap_matrix_large=matrix_s(1)%matrix,&
                                       density_matrix_large=rho_ao(ispin)%matrix,&
                                       ispin=ispin)


       ENDIF


       IF(admm_env%purification_method == do_admm_purify_cauchy)&
         CALL purify_dm_cauchy(admm_env,&
                               mo_set=mos_aux_fit(ispin)%mo_set,&
                               density_matrix=rho_ao_aux(ispin)%matrix,&
                               ispin=ispin,&
                               blocked=admm_env%block_dm)

       CALL calculate_rho_elec(ks_env=ks_env,&
                               matrix_p=rho_ao_aux(ispin)%matrix,&
                               rho=rho_r_aux(ispin),&
                               rho_gspace=rho_g_aux(ispin),&
                               total_rho=tot_rho_r_aux(ispin),&
                               soft_valid=.FALSE.,&
                               basis_type="AUX_FIT")
    END DO

    IF(dft_control%nspins==1 ) THEN
      admm_env%gsi(3) = admm_env%gsi(1)
    ELSE
      admm_env%gsi(3) = (admm_env%gsi(1)+admm_env%gsi(2))/2.0_dp
    END IF
 
    CALL qs_rho_set(rho_aux_fit, rho_r_valid=.TRUE., rho_g_valid=.TRUE.)

    CALL timestop(handle)

  END SUBROUTINE admm_mo_calc_rho_aux

! *****************************************************************************
!> \brief ...
!> \param qs_env ...
! *****************************************************************************
  SUBROUTINE admm_mo_merge_ks_matrix(qs_env)
    TYPE(qs_environment_type), POINTER       :: qs_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'admm_mo_merge_ks_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    TYPE(admm_type), POINTER                 :: admm_env

    CALL timeset(routineN,handle)
    NULLIFY(admm_env)

    CALL get_qs_env(qs_env, admm_env=admm_env)

    SELECT CASE(admm_env%purification_method)
      CASE(do_admm_purify_cauchy)
        CALL merge_ks_matrix_cauchy(qs_env)

      CASE(do_admm_purify_cauchy_subspace)
        CALL merge_ks_matrix_cauchy_subspace(qs_env)

      CASE(do_admm_purify_none)
         CALL merge_ks_matrix_none(qs_env)

      CASE(do_admm_purify_mo_diag, do_admm_purify_mo_no_diag)
         !do nothing
      CASE DEFAULT
         STOP "admm_mo_merge_ks_matrix: unknown purification method"
    END SELECT

    CALL timestop(handle)

  END SUBROUTINE admm_mo_merge_ks_matrix

! *****************************************************************************
!> \brief ...
!> \param ispin ...
!> \param admm_env ...
!> \param mo_set ...
!> \param mo_coeff ...
!> \param mo_coeff_aux_fit ...
!> \param mo_derivs ...
!> \param mo_derivs_aux_fit ...
!> \param matrix_ks_aux_fit ...
! *****************************************************************************
  SUBROUTINE admm_mo_merge_derivs(ispin, admm_env, mo_set, mo_coeff, mo_coeff_aux_fit, mo_derivs, &
                                  mo_derivs_aux_fit, matrix_ks_aux_fit)
    INTEGER, INTENT(IN)                      :: ispin
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_coeff_aux_fit
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mo_derivs, mo_derivs_aux_fit
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks_aux_fit

    CHARACTER(LEN=*), PARAMETER :: routineN = 'admm_mo_merge_derivs', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)

    SELECT CASE(admm_env%purification_method)
      CASE(do_admm_purify_mo_diag)
        CALL merge_mo_derivs_diag(ispin, admm_env, mo_set, mo_coeff, mo_coeff_aux_fit,&
                                  mo_derivs,mo_derivs_aux_fit, matrix_ks_aux_fit)

      CASE(do_admm_purify_mo_no_diag)
        CALL merge_mo_derivs_no_diag(ispin, admm_env, mo_set, mo_derivs, matrix_ks_aux_fit)

      CASE(do_admm_purify_none, do_admm_purify_cauchy,do_admm_purify_cauchy_subspace)
          !do nothing
      CASE DEFAULT
          STOP "admm_mo_merge_derivs: unknown purification method"
    END SELECT

    CALL timestop(handle)

  END SUBROUTINE admm_mo_merge_derivs


! *****************************************************************************
!> \brief ...
!> \param admm_env ...
!> \param matrix_s_aux_fit ...
!> \param matrix_s_mixed ...
!> \param mos ...
!> \param mos_aux_fit ...
!> \param geometry_did_change ...
! *****************************************************************************
  SUBROUTINE admm_fit_mo_coeffs(admm_env, matrix_s_aux_fit, matrix_s_mixed, &
                                mos, mos_aux_fit, geometry_did_change)

    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s_aux_fit, &
                                                matrix_s_mixed
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos, mos_aux_fit
    LOGICAL, INTENT(IN)                      :: geometry_did_change

    CHARACTER(LEN=*), PARAMETER :: routineN = 'admm_fit_mo_coeffs', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle

    CALL timeset(routineN,handle)

    CALL fit_mo_coeffs(admm_env, matrix_s_aux_fit, matrix_s_mixed, &
                                 mos, geometry_did_change,&
                                 blocked=admm_env%block_fit)

    SELECT CASE(admm_env%purification_method)
       CASE(do_admm_purify_mo_no_diag, do_admm_purify_cauchy_subspace)
         CALL purify_mo_cholesky(admm_env, mos, mos_aux_fit)

       CASE(do_admm_purify_mo_diag)
         CALL purify_mo_diag(admm_env, mos, mos_aux_fit)

       CASE DEFAULT
         CALL purify_mo_none(admm_env, mos, mos_aux_fit)
    END SELECT

    CALL timestop(handle)

  END SUBROUTINE admm_fit_mo_coeffs


! *****************************************************************************
!> \brief ...
!> \param admm_env ...
!> \param matrix_s_aux_fit ...
!> \param matrix_s_mixed ...
!> \param mos ...
!> \param geometry_did_change ...
!> \param blocked ...
! *****************************************************************************
  SUBROUTINE fit_mo_coeffs(admm_env, matrix_s_aux_fit, matrix_s_mixed, &
                                   mos,  geometry_did_change, blocked)
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s_aux_fit, &
                                                matrix_s_mixed
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    LOGICAL, INTENT(IN)                      :: geometry_did_change, blocked

    CHARACTER(LEN=*), PARAMETER :: routineN = 'fit_mo_coeffs', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, handle, iatom, jatom, &
                                                nao_aux_fit, nao_orb, nspins
    REAL(dp), DIMENSION(:, :), POINTER       :: sparse_block
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_dbcsr_type), POINTER             :: matrix_s_tilde

    CALL timeset(routineN,handle)

    nao_aux_fit = admm_env%nao_aux_fit
    nao_orb = admm_env%nao_orb
    nspins = SIZE(mos)


    ! *** This part only depends on overlap matrices ==> needs only to be calculated if the geometry changed

    IF( geometry_did_change ) THEN
      IF(.NOT. blocked) THEN
         CALL copy_dbcsr_to_fm(matrix_s_aux_fit(1)%matrix,admm_env%S_inv)
      ELSE
         NULLIFY(matrix_s_tilde)
         ALLOCATE(matrix_s_tilde)
         CALL cp_dbcsr_init (matrix_s_tilde)
         CALL cp_dbcsr_create(matrix_s_tilde, 'MATRIX s_tilde', &
              cp_dbcsr_distribution(matrix_s_aux_fit(1)%matrix), dbcsr_type_symmetric, &
              cp_dbcsr_row_block_sizes(matrix_s_aux_fit(1)%matrix),&
              cp_dbcsr_col_block_sizes(matrix_s_aux_fit(1)%matrix), &
              cp_dbcsr_get_data_size(matrix_s_aux_fit(1)%matrix),&
              cp_dbcsr_get_data_type(matrix_s_aux_fit(1)%matrix))

         CALL cp_dbcsr_copy(matrix_s_tilde, matrix_s_aux_fit(1)%matrix)

         CALL cp_dbcsr_iterator_start(iter, matrix_s_tilde)
         DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
           CALL cp_dbcsr_iterator_next_block(iter, iatom, jatom, sparse_block, blk)
           IF( admm_env%block_map(iatom,jatom) == 0 ) THEN
               sparse_block = 0.0_dp
           END IF
         END DO
         CALL cp_dbcsr_iterator_stop(iter)
         CALL copy_dbcsr_to_fm(matrix_s_tilde,admm_env%S_inv)
         CALL cp_dbcsr_deallocate_matrix(matrix_s_tilde)
      ENDIF

      CALL cp_fm_upper_to_full(admm_env%S_inv,admm_env%work_aux_aux)
      CALL cp_fm_to_fm(admm_env%S_inv, admm_env%S)

      CALL copy_dbcsr_to_fm(matrix_s_mixed(1)%matrix,admm_env%Q)

      !! Calculate S'_inverse
      CALL cp_fm_cholesky_decompose(admm_env%S_inv)
      CALL cp_fm_cholesky_invert(admm_env%S_inv)
      !! Symmetrize the guy
      CALL cp_fm_upper_to_full(admm_env%S_inv,admm_env%work_aux_aux)

      !! Calculate A=S'^(-1)*Q
      IF(blocked) THEN
        CALL cp_fm_set_all(admm_env%A, 0.0_dp, 1.0_dp)
      ELSE
        CALL cp_gemm('N','N',nao_aux_fit,nao_orb,nao_aux_fit,&
                      1.0_dp,admm_env%S_inv,admm_env%Q,0.0_dp,&
                      admm_env%A)

         ! this multiplication is apparent not need for purify_none
         !! B=Q^(T)*A
         CALL cp_gemm('T','N',nao_orb,nao_orb,nao_aux_fit,&
                         1.0_dp,admm_env%Q,admm_env%A,0.0_dp,&
                         admm_env%B)
      ENDIF
    END IF

    CALL timestop(handle)

  END SUBROUTINE fit_mo_coeffs


! *****************************************************************************
!> \brief Calculates the MO coefficients for the auxiliary fitting basis set
!>        by minimizing int (psi_i - psi_aux_i)^2 using Lagrangian Multipliers
!>
!> \param admm_env The ADMM env
!> \param mos the MO's of the orbital basis set
!> \param mos_aux_fit the MO's of the auxiliary fitting basis set
!> \par History
!>      05.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE purify_mo_cholesky(admm_env, mos, mos_aux_fit)

    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos, mos_aux_fit

    CHARACTER(LEN=*), PARAMETER :: routineN = 'purify_mo_cholesky', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, nao_aux_fit, &
                                                nao_orb, nmo, nspins
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_coeff_aux_fit

    CALL timeset(routineN,handle)

    nao_aux_fit = admm_env%nao_aux_fit
    nao_orb = admm_env%nao_orb
    nspins = SIZE(mos)

    ! *** Calculate the mo_coeffs for the fitting basis
    DO ispin=1,nspins
      nmo = admm_env%nmo(ispin)
      IF( nmo == 0 ) CYCLE
      !! Lambda = C^(T)*B*C
      CALL get_mo_set(mos(ispin)%mo_set,mo_coeff=mo_coeff)
      CALL get_mo_set(mos_aux_fit(ispin)%mo_set,mo_coeff=mo_coeff_aux_fit)
      CALL cp_gemm('N','N',nao_orb,nmo,nao_orb,&
                      1.0_dp,admm_env%B,mo_coeff,0.0_dp,&
                      admm_env%work_orb_nmo(ispin)%matrix)
      CALL cp_gemm('T','N',nmo,nmo,nao_orb,&
                      1.0_dp,mo_coeff,admm_env%work_orb_nmo(ispin)%matrix,0.0_dp,&
                      admm_env%lambda(ispin)%matrix)
      CALL cp_fm_to_fm(admm_env%lambda(ispin)%matrix, admm_env%work_nmo_nmo1(ispin)%matrix)



      CALL cp_fm_cholesky_decompose(admm_env%work_nmo_nmo1(ispin)%matrix)
      CALL cp_fm_cholesky_invert(admm_env%work_nmo_nmo1(ispin)%matrix)
      !! Symmetrize the guy
      CALL cp_fm_upper_to_full(admm_env%work_nmo_nmo1(ispin)%matrix,admm_env%lambda_inv(ispin)%matrix)
      CALL cp_fm_to_fm(admm_env%work_nmo_nmo1(ispin)%matrix,admm_env%lambda_inv(ispin)%matrix)


      !! ** C_hat = AC
      CALL cp_gemm('N', 'N',  nao_aux_fit, nmo, nao_orb,&
                      1.0_dp,admm_env%A,mo_coeff,0.0_dp,&
                      admm_env%C_hat(ispin)%matrix)
      CALL cp_fm_to_fm(admm_env%C_hat(ispin)%matrix, mo_coeff_aux_fit)

    END DO

    CALL timestop(handle)

  END SUBROUTINE purify_mo_cholesky


! *****************************************************************************
!> \brief Calculates the MO coefficients for the auxiliary fitting basis set
!>        by minimizing int (psi_i - psi_aux_i)^2 using Lagrangian Multipliers
!>
!> \param admm_env The ADMM env
!> \param mos the MO's of the orbital basis set
!> \param mos_aux_fit the MO's of the auxiliary fitting basis set
!> \par History
!>      05.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE purify_mo_diag(admm_env, mos, mos_aux_fit)

    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos, mos_aux_fit

    CHARACTER(LEN=*), PARAMETER :: routineN = 'purify_mo_diag', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, ispin, &
                                                nao_aux_fit, nao_orb, nmo, &
                                                nspins
    REAL(dp), ALLOCATABLE, DIMENSION(:)      :: eig_work
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_coeff_aux_fit

    CALL timeset(routineN,handle)

    nao_aux_fit = admm_env%nao_aux_fit
    nao_orb = admm_env%nao_orb
    nspins = SIZE(mos)

    ! *** Calculate the mo_coeffs for the fitting basis
    DO ispin=1,nspins
      nmo = admm_env%nmo(ispin)
      IF( nmo == 0 ) CYCLE
      !! Lambda = C^(T)*B*C
      CALL get_mo_set(mos(ispin)%mo_set,mo_coeff=mo_coeff)
      CALL get_mo_set(mos_aux_fit(ispin)%mo_set,mo_coeff=mo_coeff_aux_fit)
      CALL cp_gemm('N','N',nao_orb,nmo,nao_orb,&
                      1.0_dp,admm_env%B,mo_coeff,0.0_dp,&
                      admm_env%work_orb_nmo(ispin)%matrix)
      CALL cp_gemm('T','N',nmo,nmo,nao_orb,&
                      1.0_dp,mo_coeff,admm_env%work_orb_nmo(ispin)%matrix,0.0_dp,&
                      admm_env%lambda(ispin)%matrix)
      CALL cp_fm_to_fm(admm_env%lambda(ispin)%matrix, admm_env%work_nmo_nmo1(ispin)%matrix)


      CALL cp_fm_syevd(admm_env%work_nmo_nmo1(ispin)%matrix,admm_env%R(ispin)%matrix,&
                       admm_env%eigvals_lambda(ispin)%eigvals%data)
      ALLOCATE(eig_work(nmo))
      DO i=1,nmo
        eig_work(i) = 1.0_dp/SQRT(admm_env%eigvals_lambda(ispin)%eigvals%data(i))
      END DO
      CALL cp_fm_to_fm(admm_env%R(ispin)%matrix, admm_env%work_nmo_nmo1(ispin)%matrix)
      CALL cp_fm_column_scale(admm_env%work_nmo_nmo1(ispin)%matrix,eig_work)
      CALL cp_gemm('N','T',nmo,nmo,nmo,&
                      1.0_dp,admm_env%work_nmo_nmo1(ispin)%matrix,admm_env%R(ispin)%matrix,0.0_dp,&
                      admm_env%lambda_inv_sqrt(ispin)%matrix)
      CALL cp_gemm('N','N',nao_orb,nmo,nmo,&
                      1.0_dp,mo_coeff,admm_env%lambda_inv_sqrt(ispin)%matrix,0.0_dp,&
                      admm_env%work_orb_nmo(ispin)%matrix)
      CALL cp_gemm('N','N',nao_aux_fit,nmo,nao_orb,&
                      1.0_dp,admm_env%A,admm_env%work_orb_nmo(ispin)%matrix, 0.0_dp,&
                      mo_coeff_aux_fit)

      CALL cp_fm_to_fm(mo_coeff_aux_fit, admm_env%C_hat(ispin)%matrix)
      CALL cp_fm_set_all(admm_env%lambda_inv(ispin)%matrix,0.0_dp,1.0_dp)
      DEALLOCATE(eig_work)
    END DO

    CALL timestop(handle)

  END SUBROUTINE purify_mo_diag


! *****************************************************************************
!> \brief ...
!> \param admm_env ...
!> \param mos ...
!> \param mos_aux_fit ...
! *****************************************************************************
  SUBROUTINE purify_mo_none(admm_env, mos, mos_aux_fit)
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos, mos_aux_fit

    CHARACTER(LEN=*), PARAMETER :: routineN = 'purify_mo_none', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, nao_aux_fit, &
                                                nao_orb, nmo, nmo_mos, nspins
    REAL(KIND=dp), DIMENSION(:), POINTER     :: occ_num, occ_num_aux
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_coeff_aux_fit

    CALL timeset(routineN,handle)

    nao_aux_fit = admm_env%nao_aux_fit
    nao_orb = admm_env%nao_orb
    nspins = SIZE(mos)

    DO ispin = 1,nspins
      nmo = admm_env%nmo(ispin)
      CALL get_mo_set(mos(ispin)%mo_set,mo_coeff=mo_coeff, occupation_numbers=occ_num,nmo=nmo_mos)
      CALL get_mo_set(mos_aux_fit(ispin)%mo_set,mo_coeff=mo_coeff_aux_fit,&
           occupation_numbers=occ_num_aux)

      CALL cp_gemm('N','N',nao_aux_fit,nmo,nao_orb,&
                      1.0_dp,admm_env%A,mo_coeff,0.0_dp,&
                      mo_coeff_aux_fit)
      CALL cp_fm_to_fm(mo_coeff_aux_fit, admm_env%C_hat(ispin)%matrix)

      occ_num_aux(1:nmo) = occ_num(1:nmo)
      ! XXXX should only be done first time XXXX
      CALL cp_fm_set_all(admm_env%lambda(ispin)%matrix,0.0_dp,1.0_dp)
      CALL cp_fm_set_all(admm_env%lambda_inv(ispin)%matrix,0.0_dp,1.0_dp)
      CALL cp_fm_set_all(admm_env%lambda_inv_sqrt(ispin)%matrix,0.0_dp,1.0_dp)
    END DO

    CALL timestop(handle)

  END SUBROUTINE purify_mo_none


! *****************************************************************************
!> \brief ...
!> \param admm_env ...
!> \param mo_set ...
!> \param density_matrix ...
!> \param ispin ...
!> \param blocked ...
! *****************************************************************************
  SUBROUTINE purify_dm_cauchy(admm_env,mo_set,density_matrix,ispin,blocked)

    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_dbcsr_type), POINTER             :: density_matrix
    INTEGER                                  :: ispin
    LOGICAL, INTENT(IN)                      :: blocked

    CHARACTER(len=*), PARAMETER :: routineN = 'purify_dm_cauchy', &
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

    IF(.NOT. blocked) THEN
       CALL cp_gemm('N','T',nao_aux_fit,nao_aux_fit,nmo,&
                     1.0_dp,mo_coeff_aux_fit,mo_coeff_aux_fit,0.0_dp,&
                     admm_env%P_to_be_purified(ispin)%matrix)
    ENDIF

    CALL cp_fm_to_fm(admm_env%S, admm_env%work_aux_aux)
    CALL cp_fm_to_fm(admm_env%P_to_be_purified(ispin)%matrix, admm_env%work_aux_aux2)

    CALL cp_fm_cholesky_decompose(admm_env%work_aux_aux)

    CALL cp_fm_cholesky_reduce(admm_env%work_aux_aux2, admm_env%work_aux_aux, itype=3)

    CALL cp_fm_syevd(admm_env%work_aux_aux2,admm_env%R_purify(ispin)%matrix,&
                     admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data)

    CALL cp_fm_cholesky_restore(admm_env%R_purify(ispin)%matrix, nao_aux_fit,admm_env%work_aux_aux, &
                                admm_env%work_aux_aux3,op="MULTIPLY",pos="LEFT", transa="T")

    CALL cp_fm_to_fm(admm_env%work_aux_aux3, admm_env%R_purify(ispin)%matrix)

    ! *** Construct Matrix M for Hadamard Product
    CALL cp_fm_set_all(admm_env%M_purify(ispin)%matrix,0.0_dp)
    pole = 0.0_dp
    DO i=1,nao_aux_fit
      pole = Heaviside(admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data(i)-0.5_dp)
      CALL cp_fm_set_element(admm_env%M_purify(ispin)%matrix,i,i,pole)
    END DO
    CALL cp_fm_upper_to_full(admm_env%M_purify(ispin)%matrix,admm_env%work_aux_aux)

    CALL copy_dbcsr_to_fm(density_matrix,admm_env%work_aux_aux3)
    CALL cp_fm_upper_to_full(admm_env%work_aux_aux3,admm_env%work_aux_aux)

    ! ** S^(-1)*R
    CALL cp_gemm('N','N',nao_aux_fit,nao_aux_fit,nao_aux_fit,&
                    1.0_dp,admm_env%S_inv,admm_env%R_purify(ispin)%matrix,0.0_dp,&
                    admm_env%work_aux_aux)
    ! ** S^(-1)*R*M
    CALL cp_gemm('N','N',nao_aux_fit,nao_aux_fit,nao_aux_fit,&
                    1.0_dp,admm_env%work_aux_aux,admm_env%M_purify(ispin)%matrix,0.0_dp,&
                    admm_env%work_aux_aux2)
    ! ** S^(-1)*R*M*R^T*S^(-1)
    CALL cp_gemm('N','T',nao_aux_fit,nao_aux_fit,nao_aux_fit,&
                    1.0_dp,admm_env%work_aux_aux2,admm_env%work_aux_aux,0.0_dp,&
                    admm_env%work_aux_aux3)

    CALL copy_fm_to_dbcsr(admm_env%work_aux_aux3, density_matrix,keep_sparsity=.TRUE.)

    IF( nspins == 1 ) THEN
      CALL cp_dbcsr_scale(density_matrix, 2.0_dp)
    END IF

    CALL timestop(handle)

  END SUBROUTINE purify_dm_cauchy


! *****************************************************************************
!> \brief ...
!> \param qs_env ...
! *****************************************************************************
  SUBROUTINE merge_ks_matrix_cauchy(qs_env)
    TYPE(qs_environment_type), POINTER       :: qs_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'merge_ks_matrix_cauchy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, handle, i, iatom, ispin, &
                                                j, jatom, nao_aux_fit, &
                                                nao_orb, nmo
    REAL(dp)                                 :: eig_diff, pole, tmp
    REAL(dp), DIMENSION(:, :), POINTER       :: sparse_block
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_ks_aux_fit
    TYPE(cp_dbcsr_type), POINTER             :: matrix_k_tilde
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos

    CALL timeset(routineN,handle)
    NULLIFY(admm_env, dft_control, matrix_ks, matrix_ks_aux_fit, mos, mo_coeff)

    CALL get_qs_env(qs_env,&
                    admm_env=admm_env,&
                    dft_control=dft_control,&
                    matrix_ks=matrix_ks,&
                    matrix_ks_aux_fit=matrix_ks_aux_fit,&
                    mos=mos)

    DO ispin=1, dft_control%nspins
       nao_aux_fit = admm_env%nao_aux_fit
       nao_orb = admm_env%nao_orb
       nmo = admm_env%nmo(ispin)
       CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff)

       IF(.NOT. admm_env%block_dm) THEN
          !** Get P from mo_coeffs, otherwise we have troubles with occupation numbers ...
          CALL cp_gemm('N', 'T',  nao_orb, nao_orb, nmo,&
                          1.0_dp, mo_coeff, mo_coeff, 0.0_dp,&
                          admm_env%work_orb_orb)

          !! A*P
          CALL cp_gemm('N', 'N',  nao_aux_fit, nao_orb, nao_orb,&
                          1.0_dp, admm_env%A, admm_env%work_orb_orb, 0.0_dp,&
                          admm_env%work_aux_orb2)
          !! A*P*A^T
          CALL cp_gemm('N', 'T',  nao_aux_fit, nao_aux_fit, nao_orb,&
                          1.0_dp, admm_env%work_aux_orb2, admm_env%A, 0.0_dp,&
                          admm_env%P_to_be_purified(ispin)%matrix)

       ENDIF

       CALL cp_fm_to_fm(admm_env%S, admm_env%work_aux_aux)
       CALL cp_fm_to_fm(admm_env%P_to_be_purified(ispin)%matrix, admm_env%work_aux_aux2)

       CALL cp_fm_cholesky_decompose(admm_env%work_aux_aux)

       CALL cp_fm_cholesky_reduce(admm_env%work_aux_aux2, admm_env%work_aux_aux, itype=3)

       CALL cp_fm_syevd(admm_env%work_aux_aux2,admm_env%R_purify(ispin)%matrix,&
                        admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data)

       CALL cp_fm_cholesky_restore(admm_env%R_purify(ispin)%matrix, nao_aux_fit,admm_env%work_aux_aux, &
                                   admm_env%work_aux_aux3,op="MULTIPLY",pos="LEFT", transa="T")

       CALL cp_fm_to_fm(admm_env%work_aux_aux3, admm_env%R_purify(ispin)%matrix)

       ! *** Construct Matrix M for Hadamard Product
       pole = 0.0_dp
       DO i=1,nao_aux_fit
         DO j=i,nao_aux_fit
           eig_diff = ( admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data(i) -&
                        admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data(j) )
           ! *** two eigenvalues could be the degenerated. In that case use 2nd order formula for the poles
           IF( ABS(eig_diff) == 0.0_dp ) THEN
             pole = delta(admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data(i)-0.5_dp)
             CALL cp_fm_set_element(admm_env%M_purify(ispin)%matrix,i,j,pole)
           ELSE
             pole = 1.0_dp/(admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data(i)-&
                            admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data(j))
             tmp = Heaviside(admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data(i)-0.5_dp)
             tmp = tmp - Heaviside(admm_env%eigvals_P_to_be_purified(ispin)%eigvals%data(j)-0.5_dp)
             pole = tmp*pole
             CALL cp_fm_set_element(admm_env%M_purify(ispin)%matrix,i,j,pole)
           END IF
         END DO
       END DO
       CALL cp_fm_upper_to_full(admm_env%M_purify(ispin)%matrix,admm_env%work_aux_aux)

       CALL copy_dbcsr_to_fm(matrix_ks_aux_fit(ispin)%matrix,admm_env%K(ispin)%matrix)
       CALL cp_fm_upper_to_full(admm_env%K(ispin)%matrix,admm_env%work_aux_aux)

       !! S^(-1)*R
       CALL cp_gemm('N', 'N',  nao_aux_fit, nao_aux_fit, nao_aux_fit,&
                       1.0_dp,admm_env%S_inv,admm_env%R_purify(ispin)%matrix,0.0_dp,&
                       admm_env%work_aux_aux)
       !! K*S^(-1)*R
       CALL cp_gemm('N', 'N',  nao_aux_fit, nao_aux_fit, nao_aux_fit,&
                       1.0_dp,admm_env%K(ispin)%matrix,admm_env%work_aux_aux,0.0_dp,&
                       admm_env%work_aux_aux2)
       !! R^T*S^(-1)*K*S^(-1)*R
       CALL cp_gemm('T', 'N',  nao_aux_fit, nao_aux_fit, nao_aux_fit,&
                       1.0_dp,admm_env%work_aux_aux,admm_env%work_aux_aux2,0.0_dp,&
                       admm_env%work_aux_aux3)
       !! R^T*S^(-1)*K*S^(-1)*R x M
       CALL cp_fm_schur_product(admm_env%work_aux_aux3, admm_env%M_purify(ispin)%matrix,&
                                admm_env%work_aux_aux)

       !! R^T*A
       CALL cp_gemm('T', 'N',  nao_aux_fit, nao_orb, nao_aux_fit,&
                       1.0_dp, admm_env%R_purify(ispin)%matrix, admm_env%A, 0.0_dp,&
                       admm_env%work_aux_orb)

       !! (R^T*S^(-1)*K*S^(-1)*R x M) * R^T*A
       CALL cp_gemm('N', 'N',  nao_aux_fit, nao_orb, nao_aux_fit,&
                       1.0_dp, admm_env%work_aux_aux, admm_env%work_aux_orb, 0.0_dp,&
                       admm_env%work_aux_orb2)
       !! A^T*R*(R^T*S^(-1)*K*S^(-1)*R x M) * R^T*A
       CALL cp_gemm('T', 'N',  nao_orb, nao_orb, nao_aux_fit,&
                       1.0_dp, admm_env%work_aux_orb, admm_env%work_aux_orb2, 0.0_dp,&
                       admm_env%work_orb_orb)

       NULLIFY(matrix_k_tilde)
       ALLOCATE(matrix_k_tilde)
       CALL cp_dbcsr_init (matrix_k_tilde)
       CALL cp_dbcsr_create(matrix_k_tilde, 'MATRIX K_tilde', &
            cp_dbcsr_distribution(matrix_ks(ispin)%matrix), &
            dbcsr_type_symmetric,&
            cp_dbcsr_row_block_sizes(matrix_ks(ispin)%matrix),&
            cp_dbcsr_col_block_sizes(matrix_ks(ispin)%matrix),&
            cp_dbcsr_get_data_size( matrix_ks(ispin)%matrix),&
            cp_dbcsr_get_data_type(matrix_ks(ispin)%matrix))

       CALL cp_fm_to_fm(admm_env%work_orb_orb, admm_env%ks_to_be_merged(ispin)%matrix)

       CALL cp_dbcsr_copy(matrix_k_tilde, matrix_ks(ispin)%matrix)
       CALL cp_dbcsr_set(matrix_k_tilde, 0.0_dp)
       CALL copy_fm_to_dbcsr(admm_env%work_orb_orb, matrix_k_tilde, keep_sparsity=.TRUE.)

       IF(admm_env%block_dm) THEN
          ! ** now loop through the list and nullify blocks
          CALL cp_dbcsr_iterator_start(iter, matrix_k_tilde)
          DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
            CALL cp_dbcsr_iterator_next_block(iter, iatom, jatom, sparse_block, blk)
            IF( admm_env%block_map(iatom,jatom) == 0 ) THEN
               sparse_block = 0.0_dp
            END IF
          END DO
          CALL cp_dbcsr_iterator_stop(iter)
       ENDIF

       CALL cp_dbcsr_add(matrix_ks(ispin)%matrix, matrix_k_tilde, 1.0_dp, 1.0_dp)

       CALL cp_dbcsr_deallocate_matrix(matrix_k_tilde)

    END DO !spin-loop

    CALL timestop(handle)

  END SUBROUTINE merge_ks_matrix_cauchy


! *****************************************************************************
!> \brief ...
!> \param qs_env ...
! *****************************************************************************
  SUBROUTINE merge_ks_matrix_cauchy_subspace(qs_env)
    TYPE(qs_environment_type), POINTER       :: qs_env

    CHARACTER(LEN=*), PARAMETER :: &
      routineN = 'merge_ks_matrix_cauchy_subspace', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, nao_aux_fit, &
                                                nao_orb, nmo
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks, matrix_ks_aux_fit
    TYPE(cp_dbcsr_type), POINTER             :: matrix_k_tilde
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_coeff_aux_fit
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos, mos_aux_fit

    CALL timeset(routineN,handle)
    NULLIFY(admm_env, dft_control, matrix_ks, matrix_ks_aux_fit, mos, mos_aux_fit,&
            mo_coeff, mo_coeff_aux_fit)

    CALL get_qs_env(qs_env,&
                    admm_env=admm_env,&
                    dft_control=dft_control,&
                    matrix_ks=matrix_ks,&
                    matrix_ks_aux_fit=matrix_ks_aux_fit,&
                    mos=mos,&
                    mos_aux_fit=mos_aux_fit)

    DO ispin=1, dft_control%nspins
       nao_aux_fit = admm_env%nao_aux_fit
       nao_orb = admm_env%nao_orb
       nmo = admm_env%nmo(ispin)
       CALL get_mo_set(mo_set=mos(ispin)%mo_set,mo_coeff=mo_coeff)
       CALL get_mo_set(mo_set=mos_aux_fit(ispin)%mo_set,mo_coeff=mo_coeff_aux_fit)

       !! Calculate Lambda^{-2}
       CALL cp_fm_to_fm(admm_env%lambda(ispin)%matrix, admm_env%work_nmo_nmo1(ispin)%matrix)
       CALL cp_fm_cholesky_decompose(admm_env%work_nmo_nmo1(ispin)%matrix)
       CALL cp_fm_cholesky_invert(admm_env%work_nmo_nmo1(ispin)%matrix)
       !! Symmetrize the guy
       CALL cp_fm_upper_to_full(admm_env%work_nmo_nmo1(ispin)%matrix,admm_env%lambda_inv2(ispin)%matrix)
       !! Take square
       CALL cp_gemm('N', 'T',  nmo, nmo, nmo,&
                       1.0_dp,admm_env%work_nmo_nmo1(ispin)%matrix,admm_env%work_nmo_nmo1(ispin)%matrix,0.0_dp,&
                       admm_env%lambda_inv2(ispin)%matrix)

       !! ** C_hat = AC
       CALL cp_gemm('N', 'N',  nao_aux_fit, nmo, nao_orb,&
                       1.0_dp,admm_env%A,mo_coeff,0.0_dp,&
                       admm_env%C_hat(ispin)%matrix)

       !! calc P_tilde from C_hat
       CALL cp_gemm('N', 'N',  nao_aux_fit, nmo, nmo,&
                       1.0_dp,admm_env%C_hat(ispin)%matrix,admm_env%lambda_inv(ispin)%matrix,0.0_dp,&
                       admm_env%work_aux_nmo(ispin)%matrix)

       CALL cp_gemm('N', 'T',  nao_aux_fit, nao_aux_fit, nmo,&
                       1.0_dp,admm_env%C_hat(ispin)%matrix,  admm_env%work_aux_nmo(ispin)%matrix,0.0_dp,&
                       admm_env%P_tilde(ispin)%matrix)

       !! ** C_hat*Lambda^{-2}
       CALL cp_gemm('N', 'N',  nao_aux_fit, nmo, nmo,&
                       1.0_dp,admm_env%C_hat(ispin)%matrix,admm_env%lambda_inv2(ispin)%matrix,0.0_dp,&
                       admm_env%work_aux_nmo(ispin)%matrix)

       !! ** C_hat*Lambda^{-2}*C_hat^T
       CALL cp_gemm('N', 'T',  nao_aux_fit, nao_aux_fit, nmo,&
                       1.0_dp,admm_env%work_aux_nmo(ispin)%matrix,admm_env%C_hat(ispin)%matrix,0.0_dp,&
                       admm_env%work_aux_aux)

       !! ** S*C_hat*Lambda^{-2}*C_hat^T
       CALL cp_gemm('N', 'N',  nao_aux_fit, nao_aux_fit, nao_aux_fit,&
                       1.0_dp,admm_env%S,admm_env%work_aux_aux,0.0_dp,&
                       admm_env%work_aux_aux2)

       CALL copy_dbcsr_to_fm(matrix_ks_aux_fit(ispin)%matrix,admm_env%K(ispin)%matrix)
       CALL cp_fm_upper_to_full(admm_env%K(ispin)%matrix,admm_env%work_aux_aux)

       !! ** S*C_hat*Lambda^{-2}*C_hat^T*H_tilde
       CALL cp_gemm('N', 'N',  nao_aux_fit, nao_aux_fit, nao_aux_fit,&
                       1.0_dp,admm_env%work_aux_aux2,admm_env%K(ispin)%matrix,0.0_dp,&
                       admm_env%work_aux_aux)

       !! ** P_tilde*S
       CALL cp_gemm('N', 'N',  nao_aux_fit, nao_aux_fit, nao_aux_fit,&
                       1.0_dp,admm_env%P_tilde(ispin)%matrix,admm_env%S,0.0_dp,&
                       admm_env%work_aux_aux2)

       !! ** -S*C_hat*Lambda^{-2}*C_hat^T*H_tilde*P_tilde*S
       CALL cp_gemm('N', 'N',  nao_aux_fit, nao_aux_fit, nao_aux_fit,&
                       -1.0_dp,admm_env%work_aux_aux,admm_env%work_aux_aux2,0.0_dp,&
                       admm_env%work_aux_aux3)

       !! ** -S*C_hat*Lambda^{-2}*C_hat^T*H_tilde*P_tilde*S+S*C_hat*Lambda^{-2}*C_hat^T*H_tilde
       CALL cp_fm_scale_and_add(1.0_dp,admm_env%work_aux_aux3,1.0_dp,admm_env%work_aux_aux)

       !! first_part*A
       CALL cp_gemm('N', 'N',  nao_aux_fit, nao_orb, nao_aux_fit,&
                       1.0_dp,admm_env%work_aux_aux3,admm_env%A,0.0_dp,&
                       admm_env%work_aux_orb)

       !! + first_part^T*A
       CALL cp_gemm('T', 'N',  nao_aux_fit, nao_orb, nao_aux_fit,&
                       1.0_dp,admm_env%work_aux_aux3,admm_env%A,1.0_dp,&
                       admm_env%work_aux_orb)

       !! A^T*(first+seccond)=H
       CALL cp_gemm('T', 'N',  nao_orb, nao_orb, nao_aux_fit,&
                       1.0_dp,admm_env%A,admm_env%work_aux_orb,0.0_dp,&
                       admm_env%work_orb_orb)

       NULLIFY(matrix_k_tilde)
       ALLOCATE(matrix_k_tilde)
       CALL cp_dbcsr_init (matrix_k_tilde)
       CALL cp_dbcsr_create(matrix_k_tilde, 'MATRIX K_tilde', &
            cp_dbcsr_distribution(matrix_ks(ispin)%matrix),&
            dbcsr_type_symmetric,&
            cp_dbcsr_row_block_sizes(matrix_ks(ispin)%matrix),&
            cp_dbcsr_col_block_sizes(matrix_ks(ispin)%matrix), &
            cp_dbcsr_get_data_size(matrix_ks(ispin)%matrix),&
            cp_dbcsr_get_data_type(matrix_ks(ispin)%matrix))

       CALL cp_fm_to_fm(admm_env%work_orb_orb, admm_env%ks_to_be_merged(ispin)%matrix)

       CALL cp_dbcsr_copy(matrix_k_tilde, matrix_ks(ispin)%matrix)
       CALL cp_dbcsr_set(matrix_k_tilde, 0.0_dp)
       CALL copy_fm_to_dbcsr(admm_env%work_orb_orb, matrix_k_tilde, keep_sparsity=.TRUE.)

       CALL cp_gemm('N', 'N',  nao_orb, nmo, nao_orb,&
                       1.0_dp,admm_env%work_orb_orb,mo_coeff,0.0_dp,&
                       admm_env%mo_derivs_tmp(ispin)%matrix)

       CALL cp_dbcsr_add(matrix_ks(ispin)%matrix, matrix_k_tilde, 1.0_dp, 1.0_dp)

       CALL cp_dbcsr_deallocate_matrix(matrix_k_tilde)

    ENDDO !spin loop
    CALL timestop(handle)

  END SUBROUTINE merge_ks_matrix_cauchy_subspace


! *****************************************************************************
!> \brief Calculates the product Kohn-Sham-Matrix x mo_coeff for the auxiliary
!>        basis set and transforms it into the orbital basis. This is needed
!>        in order to use OT
!>
!> \param ispin which spin to transform
!> \param admm_env The ADMM env
!> \param mo_set ...
!> \param mo_coeff the MO coefficients from the orbital basis set
!> \param mo_coeff_aux_fit the MO coefficients from the auxiliary fitting basis set
!> \param mo_derivs KS x mo_coeff from the orbital basis set to which we add the
!>        auxiliary basis set part
!> \param mo_derivs_aux_fit ...
!> \param matrix_ks_aux_fit the Kohn-Sham matrix from the auxiliary fitting basis set
!> \par History
!>      05.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE merge_mo_derivs_diag(ispin, admm_env, mo_set, mo_coeff, mo_coeff_aux_fit, mo_derivs, &
                                  mo_derivs_aux_fit, matrix_ks_aux_fit)
    INTEGER, INTENT(IN)                      :: ispin
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_fm_type), POINTER                :: mo_coeff, mo_coeff_aux_fit
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mo_derivs, mo_derivs_aux_fit
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks_aux_fit

    CHARACTER(LEN=*), PARAMETER :: routineN = 'merge_mo_derivs_diag', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, i, j, nao_aux_fit, &
                                                nao_orb, nmo
    REAL(dp)                                 :: eig_diff, pole, tmp32, tmp52, &
                                                tmp72, tmp92
    REAL(dp), DIMENSION(:), POINTER          :: occupation_numbers, &
                                                scaling_factor

    CALL timeset(routineN,handle)

    nao_aux_fit = admm_env%nao_aux_fit
    nao_orb = admm_env%nao_orb
    nmo = admm_env%nmo(ispin)

    CALL copy_dbcsr_to_fm(matrix_ks_aux_fit(ispin)%matrix,admm_env%K(ispin)%matrix)
    CALL cp_fm_upper_to_full(admm_env%K(ispin)%matrix,admm_env%work_aux_aux)

    CALL cp_gemm('N','N', nao_aux_fit, nmo, nao_aux_fit,&
                    1.0_dp,admm_env%K(ispin)%matrix,mo_coeff_aux_fit,0.0_dp,&
                    admm_env%H(ispin)%matrix)

    CALL get_mo_set(mo_set=mo_set, occupation_numbers=occupation_numbers)
    ALLOCATE(scaling_factor(SIZE(occupation_numbers)))
    scaling_factor = 2.0_dp*occupation_numbers

    CALL cp_fm_column_scale(admm_env%H(ispin)%matrix,scaling_factor)

    CALL cp_fm_to_fm(admm_env%H(ispin)%matrix, mo_derivs_aux_fit(ispin)%matrix)

    ! *** Add first term
    CALL cp_gemm('N','T', nao_aux_fit, nmo, nmo,&
                    1.0_dp,admm_env%H(ispin)%matrix,admm_env%lambda_inv_sqrt(ispin)%matrix,0.0_dp,&
                    admm_env%work_aux_nmo(ispin)%matrix)
    CALL cp_gemm('T','N', nao_orb, nmo, nao_aux_fit,&
                    1.0_dp,admm_env%A,admm_env%work_aux_nmo(ispin)%matrix,0.0_dp,&
                    admm_env%mo_derivs_tmp(ispin)%matrix)


    ! *** Construct Matrix M for Hadamard Product
    pole = 0.0_dp
    DO i=1,nmo
      DO j=i,nmo
        eig_diff = ( admm_env%eigvals_lambda(ispin)%eigvals%data(i) -&
                     admm_env%eigvals_lambda(ispin)%eigvals%data(j) )
        ! *** two eigenvalues could be the degenerated. In that case use 2nd order formula for the poles
        IF( ABS(eig_diff) < 0.0001_dp ) THEN
          tmp32 = 1.0_dp/SQRT(admm_env%eigvals_lambda(ispin)%eigvals%data(j))**3
          tmp52 = tmp32/admm_env%eigvals_lambda(ispin)%eigvals%data(j)*eig_diff
          tmp72 = tmp52/admm_env%eigvals_lambda(ispin)%eigvals%data(j)*eig_diff
          tmp92 = tmp72/admm_env%eigvals_lambda(ispin)%eigvals%data(j)*eig_diff

          pole = -0.5_dp*tmp32 + 3.0_dp/8.0_dp*tmp52 - 5.0_dp/16.0_dp*tmp72 + 35.0_dp/128.0_dp*tmp92
          CALL cp_fm_set_element(admm_env%M(ispin)%matrix,i,j,pole)
        ELSE
          pole = 1.0_dp/SQRT(admm_env%eigvals_lambda(ispin)%eigvals%data(i))
          pole = pole - 1.0_dp/SQRT(admm_env%eigvals_lambda(ispin)%eigvals%data(j))
          pole = pole/(admm_env%eigvals_lambda(ispin)%eigvals%data(i)-&
                       admm_env%eigvals_lambda(ispin)%eigvals%data(j))
          CALL cp_fm_set_element(admm_env%M(ispin)%matrix,i,j,pole)
        END IF
      END DO
    END DO
    CALL cp_fm_upper_to_full(admm_env%M(ispin)%matrix,admm_env%work_nmo_nmo1(ispin)%matrix)

    ! *** 2nd term to be added to fm_H

    !! Part 1: B^(T)*C* R*[R^(T)*c^(T)*A^(T)*H_aux_fit*R x M]*R^(T)
    !! Part 2: B*C*(R*[R^(T)*c^(T)*A^(T)*H_aux_fit*R x M]*R^(T))^(T)

    ! *** H'*R
    CALL cp_gemm('N','N', nao_aux_fit, nmo, nmo,&
                    1.0_dp,admm_env%H(ispin)%matrix,admm_env%R(ispin)%matrix,0.0_dp,&
                    admm_env%work_aux_nmo(ispin)%matrix)
    ! *** A^(T)*H'*R
    CALL cp_gemm('T','N', nao_orb, nmo, nao_aux_fit,&
                    1.0_dp,admm_env%A,admm_env%work_aux_nmo(ispin)%matrix,0.0_dp,&
                    admm_env%work_orb_nmo(ispin)%matrix)
    ! *** c^(T)*A^(T)*H'*R
    CALL cp_gemm('T','N', nmo, nmo, nao_orb,&
                    1.0_dp,mo_coeff,admm_env%work_orb_nmo(ispin)%matrix,0.0_dp,&
                    admm_env%work_nmo_nmo1(ispin)%matrix)
    ! *** R^(T)*c^(T)*A^(T)*H'*R
    CALL cp_gemm('T','N', nmo, nmo, nmo,&
                    1.0_dp,admm_env%R(ispin)%matrix,admm_env%work_nmo_nmo1(ispin)%matrix,0.0_dp,&
                    admm_env%work_nmo_nmo2(ispin)%matrix)
    ! *** R^(T)*c^(T)*A^(T)*H'*R x M
    CALL cp_fm_schur_product(admm_env%work_nmo_nmo2(ispin)%matrix,&
                             admm_env%M(ispin)%matrix,admm_env%work_nmo_nmo1(ispin)%matrix)
    ! *** R* (R^(T)*c^(T)*A^(T)*H'*R x M)
    CALL cp_gemm('N','N', nmo, nmo, nmo,&
                    1.0_dp,admm_env%R(ispin)%matrix,admm_env%work_nmo_nmo1(ispin)%matrix,0.0_dp,&
                    admm_env%work_nmo_nmo2(ispin)%matrix)

    ! *** R* (R^(T)*c^(T)*A^(T)*H'*R x M) *R^(T)
    CALL cp_gemm('N','T', nmo, nmo, nmo,&
                    1.0_dp,admm_env%work_nmo_nmo2(ispin)%matrix,admm_env%R(ispin)%matrix,0.0_dp,&
                    admm_env%R_schur_R_t(ispin)%matrix)

    ! *** B^(T)*c
    CALL cp_gemm('T','N', nao_orb, nmo, nao_orb,&
                    1.0_dp,admm_env%B,mo_coeff,0.0_dp,&
                    admm_env%work_orb_nmo(ispin)%matrix)

    ! *** Add first term to fm_H
    ! *** B^(T)*c* R* (R^(T)*c^(T)*A^(T)*H'*R x M) *R^(T)
    CALL cp_gemm('N','N', nao_orb, nmo, nmo,&
                    1.0_dp,admm_env%work_orb_nmo(ispin)%matrix,admm_env%R_schur_R_t(ispin)%matrix,1.0_dp,&
                    admm_env%mo_derivs_tmp(ispin)%matrix)

    ! *** Add second term to fm_H
    ! *** B*C *[ R* (R^(T)*c^(T)*A^(T)*H'*R x M) *R^(T)]^(T)
    CALL cp_gemm('N','T', nao_orb, nmo, nmo,&
                    1.0_dp,admm_env%work_orb_nmo(ispin)%matrix,admm_env%R_schur_R_t(ispin)%matrix,1.0_dp,&
                    admm_env%mo_derivs_tmp(ispin)%matrix)

    DO i = 1,SIZE(scaling_factor)
      scaling_factor(i) = 1.0_dp/scaling_factor(i)
    END DO

    CALL cp_fm_column_scale(admm_env%mo_derivs_tmp(ispin)%matrix,scaling_factor)

    CALL cp_fm_scale_and_add(1.0_dp,mo_derivs(ispin)%matrix,1.0_dp,admm_env%mo_derivs_tmp(ispin)%matrix)

    DEALLOCATE(scaling_factor)

    CALL timestop(handle)

  END SUBROUTINE merge_mo_derivs_diag


! *****************************************************************************
!> \brief ...
!> \param qs_env ...
! *****************************************************************************
  SUBROUTINE merge_ks_matrix_none(qs_env)
    TYPE(qs_environment_type), POINTER       :: qs_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'merge_ks_matrix_none', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, handle, iatom, ispin, &
                                                jatom, nao_aux_fit, nao_orb, &
                                                nmo
    REAL(dp), DIMENSION(:, :), POINTER       :: sparse_block
    REAL(KIND=dp)                            :: ener_k(2), ener_x(2), &
                                                gsi_square, trace_tmp, &
                                                trace_tmp_two
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER :: matrix_ks, &
      matrix_ks_aux_fit, matrix_ks_aux_fit_dft, matrix_ks_aux_fit_hfx, &
      matrix_s, matrix_s_aux_fit, rho_ao, rho_ao_aux
    TYPE(cp_dbcsr_type), POINTER             :: matrix_k_tilde, &
                                                matrix_ks_aux_fit_admms_tmp, &
                                                matrix_TtsT
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_rho_type), POINTER               :: rho, rho_aux_fit

    CALL timeset(routineN,handle)
    NULLIFY(admm_env, dft_control, matrix_ks, matrix_ks_aux_fit, matrix_ks_aux_fit_dft, &
            matrix_ks_aux_fit_hfx, matrix_s, matrix_s_aux_fit, rho_ao, rho_ao_aux, matrix_k_tilde,&
            matrix_TtsT, matrix_ks_aux_fit_admms_tmp, rho, rho_aux_fit, sparse_block, para_env, energy)

    CALL get_qs_env(qs_env,&
                    admm_env=admm_env,&
                    dft_control=dft_control,&
                    matrix_ks=matrix_ks,&
                    matrix_ks_aux_fit=matrix_ks_aux_fit,&
                    matrix_ks_aux_fit_dft=matrix_ks_aux_fit_dft,&
                    matrix_ks_aux_fit_hfx=matrix_ks_aux_fit_hfx,&
                    !mos=mos,&
                    !mos_aux_fit=mos_aux_fit,&
                    rho=rho,& 
                    rho_aux_fit=rho_aux_fit,&
                    matrix_s=matrix_s,&
                    matrix_s_aux_fit=matrix_s_aux_fit,&
                    energy=energy,&
                    para_env=para_env)

    CALL qs_rho_get(rho, rho_ao=rho_ao)
    CALL qs_rho_get(rho_aux_fit,&
                    rho_ao=rho_ao_aux)


    DO ispin=1, dft_control%nspins
       IF(admm_env%block_dm) THEN
          CALL cp_dbcsr_iterator_start(iter, matrix_ks_aux_fit(ispin)%matrix)
          DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
            CALL cp_dbcsr_iterator_next_block(iter, iatom, jatom, sparse_block, blk)
            IF( admm_env%block_map(iatom,jatom) == 0 ) THEN
              sparse_block = 0.0_dp
            END IF
          END DO
          CALL cp_dbcsr_iterator_stop(iter)
          CALL cp_dbcsr_add(matrix_ks(ispin)%matrix, matrix_ks_aux_fit(ispin)%matrix, 1.0_dp, 1.0_dp)

       ELSE

          nao_aux_fit = admm_env%nao_aux_fit
          nao_orb = admm_env%nao_orb
          nmo = admm_env%nmo(ispin)

          ! ADMMS: different matrix for calculating A^(T)*K*A, see Eq. (37) Merlot
          IF( admm_env%charge_constrain .AND. &
              (admm_env%scaling_model == do_admm_exch_scaling_merlot ) ) THEN
              NULLIFY(matrix_ks_aux_fit_admms_tmp)
              ALLOCATE(matrix_ks_aux_fit_admms_tmp)
              CALL cp_dbcsr_init (matrix_ks_aux_fit_admms_tmp)
              CALL cp_dbcsr_create(matrix_ks_aux_fit_admms_tmp,template=matrix_ks_aux_fit(ispin)%matrix,&
                                   name='matrix_ks_aux_fit_admms_tmp', matrix_type='s')
              ! matrix_ks_aux_fit_admms_tmp = k(d_Q)
              CALL cp_dbcsr_copy(matrix_ks_aux_fit_admms_tmp, matrix_ks_aux_fit_hfx(ispin)%matrix)

              ! matrix_ks_aux_fit_admms_tmp = k(d_Q) - gsi^2/3 x(d_Q)
              CALL cp_dbcsr_add(matrix_ks_aux_fit_admms_tmp, matrix_ks_aux_fit_dft(ispin)%matrix, &
                                1.0_dp, -(admm_env%gsi(ispin))**(2.0_dp/3.0_dp))    
              CALL copy_dbcsr_to_fm(matrix_ks_aux_fit_admms_tmp,admm_env%K(ispin)%matrix)
              CALL cp_dbcsr_deallocate_matrix(matrix_ks_aux_fit_admms_tmp)
          ELSE 
              CALL copy_dbcsr_to_fm(matrix_ks_aux_fit(ispin)%matrix,admm_env%K(ispin)%matrix)
          END IF

           
          CALL cp_fm_upper_to_full(admm_env%K(ispin)%matrix,admm_env%work_aux_aux)

          !! K*A
          CALL cp_gemm('N', 'N',  nao_aux_fit, nao_orb, nao_aux_fit,&
                          1.0_dp,admm_env%K(ispin)%matrix,admm_env%A,0.0_dp,&
                          admm_env%work_aux_orb)
          !! A^T*K*A
          CALL cp_gemm('T', 'N',  nao_orb, nao_orb, nao_aux_fit,&
                          1.0_dp,admm_env%A,admm_env%work_aux_orb,0.0_dp,&
                          admm_env%work_orb_orb)

          NULLIFY(matrix_k_tilde)
          ALLOCATE(matrix_k_tilde)
          CALL cp_dbcsr_init (matrix_k_tilde)
          CALL cp_dbcsr_create(matrix_k_tilde,template=matrix_ks(ispin)%matrix,&
                               name='MATRIX K_tilde', matrix_type='S')
          CALL cp_dbcsr_copy(matrix_k_tilde, matrix_ks(ispin)%matrix)
          CALL cp_dbcsr_set(matrix_k_tilde, 0.0_dp)
          CALL copy_fm_to_dbcsr(admm_env%work_orb_orb, matrix_k_tilde, keep_sparsity=.TRUE.)

          ! Scale matrix_K_tilde here. Then, the scaling has to be done for forces seperately
          ! Scale matrix_K_tilde by gsi for ADMMQ and ADMMS (Eqs. (27), (37) in Merlot, 2014)
          IF( admm_env%charge_constrain ) THEN
              CALL cp_dbcsr_scale(matrix_k_tilde, admm_env%gsi(ispin))
          END IF

          ! Scale matrix_K_tilde by gsi^2 for ADMMP (Eq. (35) in Merlot, 2014)
          IF( (.NOT. admm_env%charge_constrain) .AND. &
              (admm_env%scaling_model==do_admm_exch_scaling_merlot) ) THEN
              gsi_square = (admm_env%gsi(ispin))*(admm_env%gsi(ispin))
              CALL cp_dbcsr_scale(matrix_k_tilde, gsi_square)
          END IF

          admm_env%lambda_merlot(ispin) = 0

          ! Calculate LAMBDA according to Merlot, 1. IF: ADMMQ, 2. IF: ADMMP, 3. IF: ADMMS,
          IF( admm_env%charge_constrain .AND. &
               (admm_env%scaling_model == do_admm_exch_scaling_none ) ) THEN
               CALL cp_dbcsr_trace(matrix_ks_aux_fit(ispin)%matrix, rho_ao_aux(ispin)%matrix, &
                    trace_tmp, 'T', 'N') 

               ! Factor of 2 is missing compared to Eq. 28 in Merlot due to
               ! Tr(ds) = N in the code \neq 2N in Merlot
               admm_env%lambda_merlot(ispin) = trace_tmp/(admm_env%n_large_basis(ispin))

          ELSE IF( ( .NOT. admm_env%charge_constrain) .AND. &
                    (admm_env%scaling_model == do_admm_exch_scaling_merlot ) ) THEN
              IF(dft_control%nspins==2) THEN
                 CALL calc_spin_dep_aux_exch_ener(qs_env=qs_env, admm_env=admm_env, ener_k_ispin=ener_k(ispin), &
                                                  ener_x_ispin=ener_x(ispin), ispin=ispin)
                 admm_env%lambda_merlot(ispin) = 2.0_dp*(admm_env%gsi(ispin))**2* &
                                                 (ener_k(ispin)+ener_x(ispin))/(admm_env%n_large_basis(ispin))

              ELSE
                  admm_env%lambda_merlot(ispin) = 2.0_dp*(admm_env%gsi(ispin))**2* &
                                                  (energy%ex+energy%exc_aux_fit)/(admm_env%n_large_basis(ispin))
              ENDIF

          ELSE IF( admm_env%charge_constrain .AND. &
              (admm_env%scaling_model == do_admm_exch_scaling_merlot ) ) THEN
               CALL cp_dbcsr_trace(matrix_ks_aux_fit_hfx(ispin)%matrix, rho_ao_aux(ispin)%matrix, &
                    trace_tmp, 'T', 'N')
               CALL cp_dbcsr_trace(matrix_ks_aux_fit_dft(ispin)%matrix, rho_ao_aux(ispin)%matrix, &
                    trace_tmp_two, 'T', 'N')
               ! For ADMMS open-shell case we need k and x (Merlot) separately since gsi(a)\=gsi(b) 
               IF(dft_control%nspins==2) THEN
                  CALL calc_spin_dep_aux_exch_ener(qs_env=qs_env, admm_env=admm_env, ener_k_ispin=ener_k(ispin), &
                                                   ener_x_ispin=ener_x(ispin), ispin=ispin)
                  admm_env%lambda_merlot(ispin) = (trace_tmp+2.0_dp/3.0_dp*((admm_env%gsi(ispin))**(2.0_dp/3.0_dp))*ener_x(ispin)-&
                                           ((admm_env%gsi(ispin))**(2.0_dp/3.0_dp))*trace_tmp_two)/(admm_env%n_large_basis(ispin))

               ELSE
                  admm_env%lambda_merlot(ispin) = (trace_tmp+(admm_env%gsi(ispin))**(2.0_dp/3.0_dp)*&
                                                  (2.0_dp/3.0_dp*energy%exc_aux_fit-trace_tmp_two))/&
                                                  (admm_env%n_large_basis(ispin))
               END IF
          END IF   


          ! Calculate variational distribution to KS matrix according
          ! to Eqs. (27), (35) and (37) in Merlot, 2014

          IF( admm_env%charge_constrain .OR. &
              (admm_env%scaling_model == do_admm_exch_scaling_merlot ) ) THEN
             
              !! T^T*s_aux*T in (27) Merlot (T=A), as calculating A^T*K*A few lines above
              CALL copy_dbcsr_to_fm(matrix_s_aux_fit(1)%matrix,admm_env%work_aux_aux4)
              CALL cp_fm_upper_to_full(admm_env%work_aux_aux4,admm_env%work_aux_aux5)

              ! s_aux*T
              CALL cp_gemm('N', 'N',  nao_aux_fit, nao_orb, nao_aux_fit,&
                          1.0_dp,admm_env%work_aux_aux4,admm_env%A,0.0_dp,&
                          admm_env%work_aux_orb3)
              ! T^T*s_aux*T
              CALL cp_gemm('T', 'N',  nao_orb, nao_orb, nao_aux_fit,&
                          1.0_dp,admm_env%A,admm_env%work_aux_orb3,0.0_dp,&
                          admm_env%work_orb_orb3)


              NULLIFY(matrix_TtsT)
              ALLOCATE(matrix_TtsT)
              CALL cp_dbcsr_init (matrix_TtsT)
              CALL cp_dbcsr_create(matrix_TtsT,template=matrix_ks(ispin)%matrix,&
                               name='MATRIX TtsT', matrix_type='S')
              CALL cp_dbcsr_copy(matrix_TtsT, matrix_ks(ispin)%matrix)
              CALL cp_dbcsr_set(matrix_TtsT, 0.0_dp)
              CALL copy_fm_to_dbcsr(admm_env%work_orb_orb3, matrix_TtsT, keep_sparsity=.TRUE.)
 

              !Add -(gsi)*Lambda*TtsT and Lambda*S to the KS matrix according to Merlot2014

              IF( admm_env%scaling_model == do_admm_exch_scaling_merlot .OR. &
                       admm_env%charge_constrain) THEN
                  CALL cp_dbcsr_add(matrix_ks(ispin)%matrix, matrix_TtsT, 1.0_dp, &
                                    (-admm_env%lambda_merlot(ispin))*admm_env%gsi(ispin))
              END IF

              CALL cp_dbcsr_add(matrix_ks(ispin)%matrix, matrix_s(1)%matrix, 1.0_dp, admm_env%lambda_merlot(ispin))

              CALL cp_dbcsr_deallocate_matrix(matrix_TtsT)

 
          END IF    

          CALL cp_dbcsr_add(matrix_ks(ispin)%matrix, matrix_k_tilde, 1.0_dp, 1.0_dp)

          CALL cp_dbcsr_deallocate_matrix(matrix_k_tilde)

       ENDIF
    ENDDO !spin loop
     

    ! Scale energy for ADMMP and ADMMS
    IF(admm_env%scaling_model == do_admm_exch_scaling_merlot) THEN
        IF( .NOT. admm_env%charge_constrain ) THEN
     !       ener_k = ener_k*(admm_env%gsi(1))*(admm_env%gsi(1))
     !       ener_x = ener_x*(admm_env%gsi(1))*(admm_env%gsi(1))
     !        PRINT *, 'energy%ex = ', energy%ex 
            IF(dft_control%nspins==2) THEN
              energy%exc_aux_fit = 0.0_dp
              energy%ex = 0.0_dp
              DO ispin=1, dft_control%nspins
                energy%exc_aux_fit = energy%exc_aux_fit + (admm_env%gsi(ispin))**2.0_dp*ener_x(ispin)
                energy%ex = energy%ex + (admm_env%gsi(ispin))**2.0_dp*ener_k(ispin)
              END DO
            ELSE
              energy%exc_aux_fit = (admm_env%gsi(1))**2.0_dp*energy%exc_aux_fit
              energy%ex = (admm_env%gsi(1))**2.0_dp*energy%ex
            END IF

        ELSE IF( admm_env%charge_constrain ) THEN
            IF(dft_control%nspins==2) THEN
              energy%exc_aux_fit = 0.0_dp
              DO ispin=1, dft_control%nspins
                energy%exc_aux_fit = energy%exc_aux_fit + (admm_env%gsi(ispin))**(2.0_dp/3.0_dp)*ener_x(ispin)
              END DO
            ELSE 
              energy%exc_aux_fit = (admm_env%gsi(1))**(2.0_dp/3.0_dp)*energy%exc_aux_fit 
            END IF
        END IF
    END IF


    CALL timestop(handle)

  END SUBROUTINE merge_ks_matrix_none


! *****************************************************************************
!> \brief Calculate exchange correction energy (Merlot2014 Eqs. 32, 33) for every spin 
!> \param qs_env ...
!> \param admm_env ...
!> \param ener_k_ispin exact ispin (Fock) exchange in auxiliary basis
!> \param ener_x_ispin ispin DFT exchange in auxiliary basis
!> \param ispin ...
!> \author Jan Wilhelm, 12/2014
! *****************************************************************************
  SUBROUTINE calc_spin_dep_aux_exch_ener(qs_env, admm_env, ener_k_ispin, ener_x_ispin, &
                          ispin)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(admm_type), POINTER                 :: admm_env
    REAL(dp), INTENT(INOUT)                  :: ener_k_ispin, ener_x_ispin
    INTEGER, INTENT(IN)                      :: ispin

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_spin_dep_aux_exch_ener', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, myspin
    REAL(KIND=dp), DIMENSION(:), POINTER     :: tot_rho_r
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks_aux_fit_hfx, &
                                                rho_ao_aux, rho_ao_aux_buffer
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: rho_g, rho_r, v_rspace_dummy, &
                                                v_tau_rspace_dummy
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_rho_type), POINTER               :: rho_aux_fit, &
                                                rho_aux_fit_buffer
    TYPE(section_vals_type), POINTER         :: xc_section_aux

    CALL timeset(routineN,handle)


    NULLIFY(ks_env, rho_aux_fit, rho_aux_fit_buffer, &
           xc_section_aux, v_rspace_dummy, v_tau_rspace_dummy,&
           rho_ao_aux, rho_ao_aux_buffer, dft_control, &
           matrix_ks_aux_fit_hfx)

    NULLIFY(rho_g, rho_r, tot_rho_r)


    CALL get_qs_env(qs_env,&
                    ks_env=ks_env,&
                    rho_aux_fit=rho_aux_fit,&
                    rho_aux_fit_buffer=rho_aux_fit_buffer,&
                    dft_control=dft_control,&
                    matrix_ks_aux_fit_hfx=matrix_ks_aux_fit_hfx)

    CALL qs_rho_get(rho_aux_fit,&
                    rho_ao=rho_ao_aux)


    CALL qs_rho_get(rho_aux_fit_buffer,&
                    rho_ao=rho_ao_aux_buffer,&
                    rho_g=rho_g,&
                    rho_r=rho_r,&
                    tot_rho_r=tot_rho_r)


!   Calculate rho_buffer = rho_aux(ispin) to get exchange of ispin electrons
    CALL cp_dbcsr_set(rho_ao_aux_buffer(1)%matrix,0.0_dp)
    CALL cp_dbcsr_set(rho_ao_aux_buffer(2)%matrix,0.0_dp)
    CALL cp_dbcsr_add(rho_ao_aux_buffer(ispin)%matrix, &
                      rho_ao_aux(ispin)%matrix, 0.0_dp, 1.0_dp)


    ! integration for getting the spin dependent density has to done for both spins!
    DO myspin=1, dft_control%nspins
      CALL calculate_rho_elec(ks_env=ks_env,&
                              matrix_p=rho_ao_aux_buffer(myspin)%matrix,&
                              rho=rho_r(myspin),&
                              rho_gspace=rho_g(myspin),&
                              total_rho=tot_rho_r(myspin),&
                              soft_valid=.FALSE.,&
                              basis_type="AUX_FIT")

    END DO

    ! Write changes in buffer density matrix
    CALL qs_rho_set(rho_aux_fit_buffer, rho_r_valid=.TRUE., rho_g_valid=.TRUE.)


    xc_section_aux => admm_env%xc_section_aux

    ener_x_ispin = 0.0_dp

    CALL qs_vxc_create(ks_env=ks_env, rho_struct=rho_aux_fit_buffer, xc_section=xc_section_aux, &
         vxc_rho=v_rspace_dummy, vxc_tau=v_tau_rspace_dummy, exc=ener_x_ispin, &
         just_energy=.TRUE.)

    ener_k_ispin = 0.0_dp

    !! ** Calculate the exchange energy
    CALL cp_dbcsr_trace(matrix_ks_aux_fit_hfx(ispin)%matrix, rho_ao_aux_buffer(ispin)%matrix, &
                        ener_k_ispin)

    ! Divide exchange for indivivual spin by two, since the ener_k_ispin originally is total
    ! exchange of alpha and beta
    ener_k_ispin = ener_k_ispin/2.0_dp

    CALL timestop(handle)

  END SUBROUTINE calc_spin_dep_aux_exch_ener 

! *****************************************************************************
!> \brief Scale density matrix by gsi(ispin), is neede for force scaling in ADMMP 
!> \param qs_env ...
!> \param rho_ao_orb ...
!> \param scale_back ...
!> \author Jan Wilhelm, 12/2014
! *****************************************************************************
  SUBROUTINE scale_dm(qs_env, rho_ao_orb, scale_back)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_p_type), DIMENSION(:, :), &
      POINTER                                :: rho_ao_orb
    LOGICAL, INTENT(IN)                      :: scale_back

    CHARACTER(LEN=*), PARAMETER :: routineN = 'scale_dm', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, img, ispin
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(dft_control_type), POINTER          :: dft_control

    CALL timeset(routineN,handle)

    NULLIFY(admm_env, dft_control)

    CALL get_qs_env(qs_env,&
                    admm_env=admm_env,&
                    dft_control=dft_control)

    ! only for ADMMP 
    IF (admm_env%scaling_model == do_admm_exch_scaling_merlot .AND. &
                 .NOT. admm_env%charge_constrain ) THEN
      DO ispin=1, dft_control%nspins
        DO img=1, dft_control%nimages
          IF(scale_back) THEN
            CALL cp_dbcsr_scale(rho_ao_orb(ispin,img)%matrix, 1.0_dp/admm_env%gsi(ispin))     
          ELSE
            CALL cp_dbcsr_scale(rho_ao_orb(ispin,img)%matrix, admm_env%gsi(ispin))
          END IF
        END DO
      END DO
    END IF

    CALL timestop(handle)

  END SUBROUTINE scale_dm


! *****************************************************************************
!> \brief ...
!> \param ispin ...
!> \param admm_env ...
!> \param mo_set ...
!> \param mo_coeff_aux_fit ...
!> \param matrix_ks_aux_fit ...
! *****************************************************************************
  SUBROUTINE calc_aux_mo_derivs_none(ispin, admm_env, mo_set, mo_coeff_aux_fit,  &
                                  matrix_ks_aux_fit)
    INTEGER, INTENT(IN)                      :: ispin
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_fm_type), POINTER                :: mo_coeff_aux_fit
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks_aux_fit

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_aux_mo_derivs_none', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nao_aux_fit, nao_orb, &
                                                nmo
    REAL(dp), DIMENSION(:), POINTER          :: occupation_numbers, &
                                                scaling_factor

    CALL timeset(routineN,handle)

    nao_aux_fit = admm_env%nao_aux_fit
    nao_orb = admm_env%nao_orb
    nmo = admm_env%nmo(ispin)

    ! just calculate the mo derivs in the aux basis  
    ! only needs to be done on the converged ks matrix for the force calc
    ! Note with OT and purification NONE, the merging of the derivs
    ! happens implicitly because the KS matrices have been already been merged
    ! and adding them here would be double counting.

    CALL copy_dbcsr_to_fm(matrix_ks_aux_fit(ispin)%matrix,admm_env%K(ispin)%matrix)
    CALL cp_fm_upper_to_full(admm_env%K(ispin)%matrix,admm_env%work_aux_aux)

    CALL cp_gemm('N','N', nao_aux_fit, nmo, nao_aux_fit,&
                    1.0_dp,admm_env%K(ispin)%matrix,mo_coeff_aux_fit,0.0_dp,&
                    admm_env%H(ispin)%matrix)

    CALL get_mo_set(mo_set=mo_set, occupation_numbers=occupation_numbers)
    ALLOCATE(scaling_factor(SIZE(occupation_numbers)))

    scaling_factor = 2.0_dp*occupation_numbers

    CALL cp_fm_column_scale(admm_env%H(ispin)%matrix,scaling_factor)

    DEALLOCATE(scaling_factor)

    CALL timestop(handle)

  END SUBROUTINE calc_aux_mo_derivs_none


! *****************************************************************************
!> \brief ...
!> \param ispin ...
!> \param admm_env ...
!> \param mo_set ...
!> \param mo_derivs ...
!> \param matrix_ks_aux_fit ...
! *****************************************************************************
  SUBROUTINE merge_mo_derivs_no_diag(ispin, admm_env, mo_set, mo_derivs, matrix_ks_aux_fit)
    INTEGER, INTENT(IN)                      :: ispin
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_fm_p_type), DIMENSION(:), &
      POINTER                                :: mo_derivs
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_ks_aux_fit

    CHARACTER(LEN=*), PARAMETER :: routineN = 'merge_mo_derivs_no_diag', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nao_aux_fit, nao_orb, &
                                                nmo
    REAL(dp), DIMENSION(:), POINTER          :: occupation_numbers, &
                                                scaling_factor

    CALL timeset(routineN,handle)

    nao_aux_fit = admm_env%nao_aux_fit
    nao_orb = admm_env%nao_orb
    nmo = admm_env%nmo(ispin)

    CALL copy_dbcsr_to_fm(matrix_ks_aux_fit(ispin)%matrix,admm_env%K(ispin)%matrix)
    CALL cp_fm_upper_to_full(admm_env%K(ispin)%matrix,admm_env%work_aux_aux)

    CALL get_mo_set(mo_set=mo_set, occupation_numbers=occupation_numbers)
    ALLOCATE(scaling_factor(SIZE(occupation_numbers)))
    scaling_factor = 0.5_dp


    !! ** calculate first part
    CALL cp_gemm('N', 'N',  nao_aux_fit, nmo, nmo,&
                    1.0_dp,admm_env%C_hat(ispin)%matrix,admm_env%lambda_inv(ispin)%matrix,0.0_dp,&
                    admm_env%work_aux_nmo(ispin)%matrix)
    CALL cp_gemm('N', 'N',  nao_aux_fit, nmo, nao_aux_fit,&
                    1.0_dp,admm_env%K(ispin)%matrix,admm_env%work_aux_nmo(ispin)%matrix,0.0_dp,&
                    admm_env%work_aux_nmo2(ispin)%matrix)
    CALL cp_gemm('T', 'N',  nao_orb, nmo, nao_aux_fit,&
                    2.0_dp,admm_env%A,admm_env%work_aux_nmo2(ispin)%matrix,0.0_dp,&
                    admm_env%mo_derivs_tmp(ispin)%matrix)
    !! ** calculate second part
    CALL cp_gemm('T', 'N',  nmo, nmo, nao_aux_fit,&
                    1.0_dp,admm_env%work_aux_nmo(ispin)%matrix,admm_env%work_aux_nmo2(ispin)%matrix,0.0_dp,&
                    admm_env%work_orb_orb)
    CALL cp_gemm('N', 'N',  nao_aux_fit, nmo, nmo,&
                    1.0_dp,admm_env%C_hat(ispin)%matrix,admm_env%work_orb_orb,0.0_dp,&
                    admm_env%work_aux_orb)
    CALL cp_gemm('N', 'N',  nao_aux_fit, nmo, nao_aux_fit,&
                    1.0_dp,admm_env%S,admm_env%work_aux_orb,0.0_dp,&
                    admm_env%work_aux_nmo(ispin)%matrix)
    CALL cp_gemm('T', 'N',  nao_orb, nmo, nao_aux_fit,&
                    -2.0_dp,admm_env%A,admm_env%work_aux_nmo(ispin)%matrix,1.0_dp,&
                    admm_env%mo_derivs_tmp(ispin)%matrix)

    CALL cp_fm_column_scale(admm_env%mo_derivs_tmp(ispin)%matrix,scaling_factor)

    CALL cp_fm_scale_and_add(1.0_dp,mo_derivs(ispin)%matrix,1.0_dp,admm_env%mo_derivs_tmp(ispin)%matrix)

    DEALLOCATE(scaling_factor)

    CALL timestop(handle)

  END SUBROUTINE merge_mo_derivs_no_diag


! *****************************************************************************
!> \brief Calculates contribution of forces due to basis transformation
!>
!>        dE/dR = dE/dC'*dC'/dR
!>        dE/dC = Ks'*c'*occ = H'
!>
!>        dC'/dR = - tr(A*lambda^(-1/2)*H'^(T)*S^(-1) * dS'/dR)
!>                 - tr(A*C*Y^(T)*C^(T)*Q^(T)*A^(T) * dS'/dR)
!>                 + tr(C*lambda^(-1/2)*H'^(T)*S^(-1) * dQ/dR)
!>                 + tr(A*C*Y^(T)*c^(T) * dQ/dR)
!>                 + tr(C*Y^(T)*C^(T)*A^(T) * dQ/dR)
!>
!>        where
!>
!>        A = S'^(-1)*Q
!>        lambda = C^(T)*B*C
!>        B = Q^(T)*A
!>        Y = R*[ (R^(T)*C^(T)*A^(T)*H'*R) xx M ]*R^(T)
!>        lambda = R*D*R^(T)
!>        Mij = Poles-Matrix (see above)
!>        xx = schur product
!>
!> \param qs_env the QS environment
!> \par History
!>      05.2008 created [Manuel Guidon]
!> \author Manuel Guidon
! *****************************************************************************
  SUBROUTINE calc_mixed_overlap_force(qs_env)

    TYPE(qs_environment_type), POINTER       :: qs_env

    CHARACTER(LEN=*), PARAMETER :: routineN = 'calc_mixed_overlap_force', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, ispin, iw, &
                                                nao_aux_fit, nao_orb, natom, &
                                                neighbor_list_id, nmo
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :)                        :: admm_force
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(cp_dbcsr_p_type), DIMENSION(:), &
      POINTER                                :: matrix_s, matrix_s_aux_fit, &
                                                matrix_s_aux_fit_vs_orb, &
                                                rho_ao, rho_ao_aux
    TYPE(cp_dbcsr_type), POINTER             :: matrix_rho_aux_desymm_tmp, &
                                                matrix_w_q, matrix_w_s
    TYPE(cp_fm_type), POINTER                :: mo_coeff
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dft_control_type), POINTER          :: dft_control
    TYPE(mo_set_p_type), DIMENSION(:), &
      POINTER                                :: mos
    TYPE(neighbor_list_set_p_type), &
      DIMENSION(:), POINTER                  :: sab_aux_fit, &
                                                sab_aux_fit_asymm, &
                                                sab_aux_fit_vs_orb, sab_orb
    TYPE(qs_energy_type), POINTER            :: energy
    TYPE(qs_force_type), DIMENSION(:), &
      POINTER                                :: force
    TYPE(qs_ks_env_type), POINTER            :: ks_env
    TYPE(qs_rho_type), POINTER               :: rho, rho_aux_fit

    CALL timeset(routineN,handle)
    failure=.FALSE.

    NULLIFY(admm_env, logger, dft_control, para_env, mos, mo_coeff, matrix_w_q, matrix_w_s,&
            rho, rho_aux_fit, energy, sab_aux_fit, sab_aux_fit_asymm, &
            sab_aux_fit_vs_orb, sab_orb, ks_env, matrix_s_aux_fit, matrix_s_aux_fit_vs_orb, matrix_s)

    CALL get_qs_env(qs_env,&
                    admm_env=admm_env,&
                    ks_env=ks_env,&
                    dft_control=dft_control,&
                    matrix_s_aux_fit=matrix_s_aux_fit,&
                    matrix_s_aux_fit_vs_orb=matrix_s_aux_fit_vs_orb,&
                    matrix_s=matrix_s,&
                    neighbor_list_id=neighbor_list_id,&
                    rho=rho,&
                    rho_aux_fit=rho_aux_fit,&
                    energy=energy,&
                    sab_orb=sab_orb,&
                    sab_aux_fit=sab_aux_fit,&
                    sab_aux_fit_asymm=sab_aux_fit_asymm,&
                    sab_aux_fit_vs_orb=sab_aux_fit_vs_orb,&
                    mos=mos,&
                    para_env=para_env)

    CALL qs_rho_get(rho, rho_ao=rho_ao)
    CALL qs_rho_get(rho_aux_fit,&
                    rho_ao=rho_ao_aux)


    nao_aux_fit = admm_env%nao_aux_fit
    nao_orb = admm_env%nao_orb

    logger => cp_get_default_logger()

    ! *** forces are only implemented for mo_diag or none and basis_projection ***
    IF (admm_env%block_dm) THEN
      CPPostcondition(.FALSE.,cp_failure_level,routineP,failure)
    ENDIF

    IF (.NOT.(admm_env%purification_method==do_admm_purify_mo_diag .OR. &
               admm_env%purification_method==do_admm_purify_none)) THEN
      CPPostcondition(.FALSE.,cp_failure_level,routineP,failure)
    ENDIF

    ! *** Create sparse work matrices

    ALLOCATE(matrix_w_s)
    CALL cp_dbcsr_init (matrix_w_s)
    CALL cp_dbcsr_create(matrix_w_s, 'W MATRIX AUX S', &
         cp_dbcsr_distribution(matrix_s_aux_fit(1)%matrix), dbcsr_type_no_symmetry, &
         cp_dbcsr_row_block_sizes(matrix_s_aux_fit(1)%matrix),&
         cp_dbcsr_col_block_sizes(matrix_s_aux_fit(1)%matrix), &
         cp_dbcsr_get_data_size(matrix_s_aux_fit(1)%matrix),&
         cp_dbcsr_get_data_type(matrix_s_aux_fit(1)%matrix))
    CALL cp_dbcsr_alloc_block_from_nbl(matrix_w_s,sab_aux_fit_asymm)

    ALLOCATE(matrix_w_q)
    CALL cp_dbcsr_init(matrix_w_q)
    CALL cp_dbcsr_copy(matrix_w_q,matrix_s_aux_fit_vs_orb(1)%matrix,&
                    "W MATRIX AUX Q")

    DO ispin=1, dft_control%nspins
        nmo = admm_env%nmo(ispin)
        CALL get_mo_set(mo_set=mos(ispin)%mo_set, mo_coeff=mo_coeff)

        ! *** S'^(-T)*H'
        IF (.NOT. admm_env%purification_method==do_admm_purify_none) THEN
        CALL cp_gemm('T', 'N', nao_aux_fit, nmo, nao_aux_fit,&
                        1.0_dp,admm_env%S_inv,qs_env%mo_derivs_aux_fit(ispin)%matrix,0.0_dp,&
                        admm_env%work_aux_nmo(ispin)%matrix)
        ELSE

        CALL cp_gemm('T', 'N', nao_aux_fit, nmo, nao_aux_fit,&
                        1.0_dp,admm_env%S_inv,admm_env%H(ispin)%matrix,0.0_dp,&
                        admm_env%work_aux_nmo(ispin)%matrix)
        END IF

        ! *** S'^(-T)*H'*Lambda^(-T/2)
        CALL cp_gemm('N', 'T',  nao_aux_fit,nmo, nmo,&
                        1.0_dp,admm_env%work_aux_nmo(ispin)%matrix,admm_env%lambda_inv_sqrt(ispin)%matrix,0.0_dp,&
                        admm_env%work_aux_nmo2(ispin)%matrix)

        ! *** C*Lambda^(-1/2)*H'^(T)*S'^(-1) minus sign due to force = -dE/dR
        CALL cp_gemm('N', 'T',  nao_aux_fit, nao_orb, nmo,&
                        -1.0_dp,admm_env%work_aux_nmo2(ispin)%matrix,mo_coeff,0.0_dp,&
                        admm_env%work_aux_orb)

        ! *** A*C*Lambda^(-1/2)*H'^(T)*S'^(-1), minus sign to recover from above
        CALL cp_gemm('N', 'T',  nao_aux_fit, nao_aux_fit, nao_orb,&
                        -1.0_dp,admm_env%work_aux_orb,admm_env%A,0.0_dp,&
                        admm_env%work_aux_aux)


        IF (.NOT. (admm_env%purification_method==do_admm_purify_none)) THEN
           ! *** C*Y
           CALL cp_gemm('N', 'N',  nao_orb, nmo, nmo,&
                        1.0_dp,mo_coeff,admm_env%R_schur_R_t(ispin)%matrix,0.0_dp,&
                        admm_env%work_orb_nmo(ispin)%matrix)
           ! *** C*Y^(T)*C^(T)
           CALL cp_gemm('N', 'T',  nao_orb, nao_orb, nmo,&
                        1.0_dp,mo_coeff,admm_env%work_orb_nmo(ispin)%matrix,0.0_dp,&
                        admm_env%work_orb_orb)
           ! *** A*C*Y^(T)*C^(T) Add to work aux_orb, minus sign due to force = -dE/dR
           CALL cp_gemm('N', 'N',  nao_aux_fit, nao_orb, nao_orb,&
                        -1.0_dp,admm_env%A,admm_env%work_orb_orb,1.0_dp,&
                        admm_env%work_aux_orb)

           ! *** C*Y^(T)
           CALL cp_gemm('N', 'T',  nao_orb, nmo, nmo,&
                        1.0_dp,mo_coeff,admm_env%R_schur_R_t(ispin)%matrix,0.0_dp,&
                        admm_env%work_orb_nmo(ispin)%matrix)
           ! *** C*Y*C^(T)
           CALL cp_gemm('N', 'T',  nao_orb, nao_orb, nmo,&
                        1.0_dp,mo_coeff,admm_env%work_orb_nmo(ispin)%matrix,0.0_dp,&
                        admm_env%work_orb_orb)
           ! *** A*C*Y*C^(T) Add to work aux_orb, minus sign due to -dE/dR
           CALL cp_gemm('N', 'N',  nao_aux_fit, nao_orb, nao_orb,&
                        -1.0_dp,admm_env%A,admm_env%work_orb_orb,1.0_dp,&
                        admm_env%work_aux_orb)
        END IF

        ! Add derivative contribution matrix*dQ/dR in additional last term in
        ! Eq. (26,32, 33) in Merlot2014 to the force
        ! ADMMS
        IF ( admm_env%scaling_model == do_admm_exch_scaling_merlot .AND. &
             admm_env%charge_constrain  )  THEN
           ! *** scale admm_env%work_aux_orb by gsi due to inner derivative
           CALL cp_fm_scale(admm_env%gsi(ispin), admm_env%work_aux_orb)
           ! ***  as in ADMMP only with different sign
           CALL cp_gemm('N', 'T',  nao_orb, nao_orb, nmo, &
                   4.0_dp*( admm_env%gsi(ispin) )*admm_env%lambda_merlot(ispin)/dft_control%nspins, &
                   mo_coeff,mo_coeff,0.0_dp, admm_env%work_orb_orb2)

           ! *** prefactor*A*C*C^(T) Add to work aux_orb
           CALL cp_gemm('N', 'N', nao_aux_fit, nao_orb, nao_orb,&
                        1.0_dp, admm_env%A, admm_env%work_orb_orb2,1.0_dp,&
                        admm_env%work_aux_orb)


        ! ADMMP  
        ELSE IF (admm_env%scaling_model == do_admm_exch_scaling_merlot .AND. &
                 .NOT. admm_env%charge_constrain ) THEN
           ! *** prefactor*C*C^(T), nspins since 2/n_spin*C*C^(T)=P 
           CALL cp_gemm('N', 'T',  nao_orb, nao_orb, nmo, &
                   -4.0_dp*( admm_env%gsi(ispin) )*admm_env%lambda_merlot(ispin)/dft_control%nspins, &
                   mo_coeff,mo_coeff,0.0_dp, admm_env%work_orb_orb2)

           ! *** prefactor*A*C*C^(T) Add to work aux_orb
           CALL cp_gemm('N', 'N', nao_aux_fit, nao_orb, nao_orb,&
                        1.0_dp, admm_env%A, admm_env%work_orb_orb2,1.0_dp,&
                        admm_env%work_aux_orb)          
  

        ! ADMMQ 
        ELSE IF (admm_env%scaling_model == do_admm_exch_scaling_none .AND. &
                 admm_env%charge_constrain ) THEN
           ! *** scale admm_env%work_aux_orb by gsi due to inner derivative
           CALL cp_fm_scale(admm_env%gsi(ispin), admm_env%work_aux_orb) 
           ! ***  as in ADMMP only with different sign
           CALL cp_gemm('N', 'T',  nao_orb, nao_orb, nmo, &
                   4.0_dp*( admm_env%gsi(ispin) )*admm_env%lambda_merlot(ispin)/dft_control%nspins, &
                   mo_coeff,mo_coeff,0.0_dp, admm_env%work_orb_orb2)

           ! *** prefactor*A*C*C^(T) Add to work aux_orb
           CALL cp_gemm('N', 'N', nao_aux_fit, nao_orb, nao_orb,&
                        1.0_dp, admm_env%A, admm_env%work_orb_orb2,1.0_dp,&
                        admm_env%work_aux_orb)
        END IF



        ! *** copy to sparse matrix
        CALL copy_fm_to_dbcsr(admm_env%work_aux_orb, matrix_w_q,keep_sparsity=.TRUE.)

        IF (.NOT. (admm_env%purification_method==do_admm_purify_none)) THEN
           ! *** A*C*Y^(T)*C^(T)
           CALL cp_gemm('N', 'N',  nao_aux_fit, nao_orb, nao_orb,&
                        1.0_dp,admm_env%A,admm_env%work_orb_orb,0.0_dp,&
                        admm_env%work_aux_orb)
           ! *** A*C*Y^(T)*C^(T)*A^(T) add to aux_aux, minus sign cancels
           CALL cp_gemm('N', 'T',  nao_aux_fit, nao_aux_fit, nao_orb,&
                        1.0_dp,admm_env%work_aux_orb,admm_env%A,1.0_dp,&
                        admm_env%work_aux_aux)
        END IF

        ! *** copy to sparse matrix
        CALL copy_fm_to_dbcsr(admm_env%work_aux_aux, matrix_w_s,keep_sparsity=.TRUE.)

        ! Add derivative of Eq. (33) with respect to s_aux Merlot2014 to the force
        IF ( admm_env%scaling_model == do_admm_exch_scaling_merlot .OR. &
             admm_env%charge_constrain )  THEN

           !Create desymmetrized auxiliary density matrix           
           NULLIFY(matrix_rho_aux_desymm_tmp) 
           ALLOCATE(matrix_rho_aux_desymm_tmp)
           CALL cp_dbcsr_init (matrix_rho_aux_desymm_tmp)
           CALL cp_dbcsr_create(matrix_rho_aux_desymm_tmp, 'Rho_aux non-symm', &
                cp_dbcsr_distribution(matrix_s_aux_fit(1)%matrix), dbcsr_type_no_symmetry, &
                cp_dbcsr_row_block_sizes(matrix_s_aux_fit(1)%matrix),&
                cp_dbcsr_col_block_sizes(matrix_s_aux_fit(1)%matrix), &
                cp_dbcsr_get_data_size(matrix_s_aux_fit(1)%matrix),&
                cp_dbcsr_get_data_type(matrix_s_aux_fit(1)%matrix))

           CALL cp_dbcsr_desymmetrize(rho_ao_aux(ispin)%matrix, matrix_rho_aux_desymm_tmp)
          
 
           ! ADMMS 1. scale original matrix_w_s by gsi due to inner deriv.
           !       2. add derivative of variational term with resp. to s
           IF ( admm_env%scaling_model == do_admm_exch_scaling_merlot .AND. &
                admm_env%charge_constrain ) THEN
             CALL cp_dbcsr_scale(matrix_w_s, admm_env%gsi(ispin))
             CALL cp_dbcsr_add(matrix_w_s, matrix_rho_aux_desymm_tmp, 1.0_dp, &
                                -admm_env%lambda_merlot(ispin))
           
          
           ! ADMMP add derivative of variational term with resp. to s
           ELSE IF ( admm_env%scaling_model == do_admm_exch_scaling_merlot .AND. & 
                     .NOT. admm_env%charge_constrain )  THEN

             CALL cp_dbcsr_add(matrix_w_s, matrix_rho_aux_desymm_tmp, 1.0_dp, &
                                (admm_env%gsi(ispin))*admm_env%lambda_merlot(ispin))


           ! ADMMQ 1. scale original matrix_w_s by gsi due to inner deriv.
           !       2. add derivative of variational term with resp. to s
           ELSE IF ( admm_env%scaling_model == do_admm_exch_scaling_none .AND.  &
                     admm_env%charge_constrain )  THEN
             CALL cp_dbcsr_scale(matrix_w_s, admm_env%gsi(ispin))
             CALL cp_dbcsr_add(matrix_w_s, matrix_rho_aux_desymm_tmp, 1.0_dp, &
                                -admm_env%lambda_merlot(ispin))
         
           END IF

           CALL cp_dbcsr_deallocate_matrix(matrix_rho_aux_desymm_tmp)


        END IF



        ! allocate force vector
        CALL get_qs_env(qs_env=qs_env,natom=natom)
        ALLOCATE(admm_force(3,natom))
        admm_force = 0.0_dp
        CALL build_overlap_force(ks_env, admm_force,&
             basis_type_a="AUX_FIT", basis_type_b="AUX_FIT", &
             sab_nl=sab_aux_fit_asymm, matrix_p=matrix_w_s)
        CALL build_overlap_force(ks_env, admm_force,&
             basis_type_a="AUX_FIT", basis_type_b="ORB", &
             sab_nl=sab_aux_fit_vs_orb, matrix_p=matrix_w_q)

        ! Add contribution of original basis set for ADMMQ
        IF ( .NOT. admm_env%scaling_model == do_admm_exch_scaling_merlot .AND. admm_env%charge_constrain )  THEN
           CALL cp_dbcsr_scale(rho_ao(ispin)%matrix, -admm_env%lambda_merlot(ispin))
           CALL build_overlap_force(ks_env, admm_force,&
                basis_type_a="ORB", basis_type_b="ORB", &
                sab_nl=sab_orb, matrix_p=rho_ao(ispin)%matrix)
           CALL cp_dbcsr_scale(rho_ao(ispin)%matrix, -1.0_dp/admm_env%lambda_merlot(ispin))
       END IF

        ! Add contribution of original basis set for ADMMP
        IF ( admm_env%scaling_model == do_admm_exch_scaling_merlot .AND. .NOT. admm_env%charge_constrain )  THEN
           CALL cp_dbcsr_scale(rho_ao(ispin)%matrix, admm_env%lambda_merlot(ispin))
           CALL build_overlap_force(ks_env, admm_force,&
                basis_type_a="ORB", basis_type_b="ORB", &
                sab_nl=sab_orb, matrix_p=rho_ao(ispin)%matrix)
           CALL cp_dbcsr_scale(rho_ao(ispin)%matrix, 1.0_dp/admm_env%lambda_merlot(ispin))
       END IF

        ! Add contribution of original basis set for ADMMS
        IF ( admm_env%scaling_model == do_admm_exch_scaling_merlot .AND. admm_env%charge_constrain )  THEN
           CALL cp_dbcsr_scale(rho_ao(ispin)%matrix, -admm_env%lambda_merlot(ispin))
           CALL build_overlap_force(ks_env, admm_force,&
                basis_type_a="ORB", basis_type_b="ORB", &
                sab_nl=sab_orb, matrix_p=rho_ao(ispin)%matrix)
           CALL cp_dbcsr_scale(rho_ao(ispin)%matrix, -1.0_dp/admm_env%lambda_merlot(ispin))
       END IF





        ! add forces
        CALL get_qs_env(qs_env=qs_env,atomic_kind_set=atomic_kind_set,&
                        force=force)
        CALL add_qs_force(admm_force, force, "overlap_admm", atomic_kind_set)
        DEALLOCATE(admm_force)

        IF (BTEST(cp_print_key_should_output(logger%iter_info,&
             qs_env%input,"DFT%PRINT%AO_MATRICES/W_MATRIX_AUX_FIT"),cp_p_file)) THEN
           iw = cp_print_key_unit_nr(logger,qs_env%input,"DFT%PRINT%AO_MATRICES/W_MATRIX_AUX_FIT",&
                extension=".Log")
           CALL cp_dbcsr_write_sparse_matrix(matrix_w_s,4,6,qs_env,para_env,output_unit=iw)
           CALL cp_print_key_finished_output(iw,logger,qs_env%input,&
                "DFT%PRINT%AO_MATRICES/W_MATRIX_AUX_FIT")
        END IF
        IF (BTEST(cp_print_key_should_output(logger%iter_info,&
             qs_env%input,"DFT%PRINT%AO_MATRICES/W_MATRIX_AUX_FIT"),cp_p_file)) THEN
           iw = cp_print_key_unit_nr(logger,qs_env%input,"DFT%PRINT%AO_MATRICES/W_MATRIX_AUX_FIT",&
                extension=".Log")
           CALL cp_dbcsr_write_sparse_matrix(matrix_w_q,4,6,qs_env,para_env,output_unit=iw)
           CALL cp_print_key_finished_output(iw,logger,qs_env%input,&
                "DFT%PRINT%AO_MATRICES/W_MATRIX_AUX_FIT")
        END IF

    ENDDO !spin loop

    ! *** Deallocated weighted density matrices
    CALL cp_dbcsr_deallocate_matrix(matrix_w_s)
    CALL cp_dbcsr_deallocate_matrix(matrix_w_q)

    CALL timestop(handle)

  END SUBROUTINE calc_mixed_overlap_force


! *****************************************************************************
!> \brief ...
!> \param admm_env environment of auxiliary DM 
!> \param mo_set ...
!> \param density_matrix auxiliary DM
!> \param overlap_matrix auxiliary OM
!> \param density_matrix_large DM of the original basis
!> \param overlap_matrix_large overlap matrix of original basis
!> \param ispin ...
! *****************************************************************************
  SUBROUTINE calculate_dm_mo_no_diag(admm_env,mo_set,density_matrix,overlap_matrix,&
             density_matrix_large,overlap_matrix_large,ispin)
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(mo_set_type), POINTER               :: mo_set
    TYPE(cp_dbcsr_type), POINTER             :: density_matrix, &
                                                overlap_matrix, &
                                                density_matrix_large, &
                                                overlap_matrix_large
    INTEGER                                  :: ispin

    CHARACTER(len=*), PARAMETER :: routineN = 'calculate_dm_mo_no_diag', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nao_aux_fit, nmo
    REAL(KIND=dp)                            :: alpha, nel_tmp_aux

! Number of electrons in the aux. DM

    CALL timeset(routineN,handle)

    CALL cp_dbcsr_set(density_matrix,0.0_dp)
    nao_aux_fit = admm_env%nao_aux_fit
    nmo = admm_env%nmo(ispin)
    CALL cp_fm_to_fm(admm_env%C_hat(ispin)%matrix, admm_env%work_aux_nmo(ispin)%matrix)
    CALL cp_fm_column_scale(admm_env%work_aux_nmo(ispin)%matrix,mo_set%occupation_numbers(1:mo_set%homo))

    CALL cp_gemm('N','N',nao_aux_fit,nmo,nmo,&
                    1.0_dp,admm_env%work_aux_nmo(ispin)%matrix,admm_env%lambda_inv(ispin)%matrix,0.0_dp,&
                    admm_env%work_aux_nmo2(ispin)%matrix)

    ! The following IF doesn't do anything unless !alpha=mo_set%maxocc is uncommented.
    IF ( .NOT. mo_set%uniform_occupation ) THEN ! not all orbitals 1..homo are equally occupied
      alpha=1.0_dp
      CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=density_matrix,&
                              matrix_v=admm_env%C_hat(ispin)%matrix,&
                              matrix_g=admm_env%work_aux_nmo2(ispin)%matrix,&
                              ncol=mo_set%homo,&
                              alpha=alpha)
    ELSE
      alpha=1.0_dp
      !alpha=mo_set%maxocc
      CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=density_matrix,&
                              matrix_v=admm_env%C_hat(ispin)%matrix,&
                              matrix_g=admm_env%work_aux_nmo2(ispin)%matrix,&
                              ncol=mo_set%homo,&
                              alpha=alpha)
    ENDIF


    !  The following IF checks whether gsi needs to be calculated. This is the case if
    !   the auxiliary density matrix gets scaled
    !   according to Eq. 22 (Merlot) or a scaling of exchange_correction is employed, Eq. 35 (Merlot).  
    IF( admm_env%charge_constrain .OR. (admm_env%scaling_model == do_admm_exch_scaling_merlot ) ) THEN

        CALL cite_reference(Merlot2014)

        admm_env%n_large_basis(3) = 0.0_dp

        ! Calculate number of electrons in the original density matrix, transposing doesn't matter
        ! since both matrices are symmetric
        CALL cp_dbcsr_trace(density_matrix_large, overlap_matrix_large, admm_env%n_large_basis(ispin), 'T', 'N')
        admm_env%n_large_basis(3) = admm_env%n_large_basis(3) + admm_env%n_large_basis(ispin)
        ! Calculate number of electrons in the auxiliary density matrix
        CALL cp_dbcsr_trace(density_matrix, overlap_matrix, nel_tmp_aux, 'T', 'N')
        admm_env%gsi(ispin) = admm_env%n_large_basis(ispin)/nel_tmp_aux

        IF( admm_env%charge_constrain ) THEN
          ! multiply aux. DM with gsi to get the scaled DM (Merlot, Eq. 21)
          CALL cp_dbcsr_scale(density_matrix, admm_env%gsi(ispin))
        ENDIF

    ENDIF


   
    CALL timestop(handle)

  END SUBROUTINE calculate_dm_mo_no_diag

! *****************************************************************************
!> \brief ...
!> \param admm_env ...
!> \param density_matrix ...
!> \param density_matrix_aux ...
!> \param ispin ...
!> \param nspins ...
! *****************************************************************************
  SUBROUTINE blockify_density_matrix(admm_env,density_matrix, density_matrix_aux,&
                                     ispin, nspins)
    TYPE(admm_type), POINTER                 :: admm_env
    TYPE(cp_dbcsr_type), POINTER             :: density_matrix, &
                                                density_matrix_aux
    INTEGER                                  :: ispin, nspins

    CHARACTER(len=*), PARAMETER :: routineN = 'blockify_density_matrix', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: blk, handle, iatom, jatom
    LOGICAL                                  :: found
    REAL(dp), DIMENSION(:, :), POINTER       :: sparse_block, sparse_block_aux
    TYPE(cp_dbcsr_iterator)                  :: iter

    CALL timeset(routineN,handle)

   ! ** set blocked density matrix to 0
   CALL cp_dbcsr_set(density_matrix_aux, 0.0_dp)

   ! ** now loop through the list and copy corresponding blocks
   CALL cp_dbcsr_iterator_start(iter, density_matrix)
   DO WHILE (cp_dbcsr_iterator_blocks_left (iter))
     CALL cp_dbcsr_iterator_next_block(iter, iatom, jatom, sparse_block, blk)
     IF( admm_env%block_map(iatom,jatom) == 1 ) THEN
       CALL cp_dbcsr_get_block_p(density_matrix_aux,&
                              row=iatom,col=jatom,BLOCK=sparse_block_aux,found=found)
       IF( found ) THEN
         sparse_block_aux = sparse_block
       END IF

     END IF
   END DO
   CALL cp_dbcsr_iterator_stop(iter)

   CALL copy_dbcsr_to_fm(density_matrix_aux,admm_env%P_to_be_purified(ispin)%matrix)
   CALL cp_fm_upper_to_full(admm_env%P_to_be_purified(ispin)%matrix, admm_env%work_orb_orb2)

   IF( nspins == 1 ) THEN
     CALL cp_fm_scale(0.5_dp, admm_env%P_to_be_purified(ispin)%matrix)
   END IF

   CALL timestop(handle)
  END SUBROUTINE blockify_density_matrix

! *****************************************************************************
!> \brief ...
!> \param x ...
!> \retval delta ...
! *****************************************************************************
  FUNCTION delta(x)
    REAL(KIND=dp), INTENT(IN)                :: x
    REAL(KIND=dp)                            :: delta

    IF( x == 0.0_dp) THEN !TODO: exact comparison of reals?
      delta = 1.0_dp
    ELSE
      delta = 0.0_dp
    END IF

  END FUNCTION delta

! *****************************************************************************
!> \brief ...
!> \param x ...
!> \retval Heaviside ...
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

END MODULE admm_methods
