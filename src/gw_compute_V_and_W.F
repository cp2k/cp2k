!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Coulomb matrices and screened interactions for imaginary-time GW.
! **************************************************************************************************
MODULE gw_compute_V_and_W
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cell_types,                      ONLY: cell_type
   USE constants_operator,              ONLY: operator_coulomb
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_add_on_diag,&
                                              cp_cfm_uplo_to_full
   USE cp_cfm_cholesky,                 ONLY: cp_cfm_cholesky_decompose,&
                                              cp_cfm_cholesky_invert
   USE cp_cfm_diag,                     ONLY: cp_cfm_power
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_info,&
                                              cp_cfm_release,&
                                              cp_cfm_to_cfm,&
                                              cp_cfm_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_add, dbcsr_create, dbcsr_deallocate_matrix, dbcsr_distribution_release, &
        dbcsr_distribution_type, dbcsr_p_type, dbcsr_release, dbcsr_set, dbcsr_type, &
        dbcsr_type_no_symmetry
   USE cp_dbcsr_contrib,                ONLY: dbcsr_reserve_all_blocks
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              cp_dbcsr_dist2d_to_dist,&
                                              dbcsr_deallocate_matrix_set
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_fm_basic_linalg,              ONLY: cp_fm_add_on_diag,&
                                              cp_fm_scale_and_add,&
                                              cp_fm_uplo_to_full
   USE cp_fm_redistribution,            ONLY: cp_fm_to_replicated_array
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_read_unformatted,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_type,&
                                              cp_fm_write_unformatted
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE distribution_2d_types,           ONLY: distribution_2d_type
   USE gw_utils_fm,                     ONLY: cfm_contract_ABA,&
                                              fm_contract_ABA,&
                                              fm_invert,&
                                              fm_sqrt
   USE input_constants,                 ONLY: rtp_method_bse,&
                                              rtp_method_bse_linearized
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE kpoint_coulomb_2c,               ONLY: build_2c_coulomb_matrix_kp
   USE libint_2c_3c,                    ONLY: libint_potential_type
   USE machine,                         ONLY: m_walltime
   USE mathconstants,                   ONLY: z_one,&
                                              z_zero
   USE mp2_ri_2c,                       ONLY: RI_2c_integral_mat_kp
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE post_scf_bandstructure_utils,    ONLY: MIC_contribution_from_ikp,&
                                              cfm_ikp_from_fm_Gamma
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type,&
                                              release_neighbor_list_sets
   USE qs_tensors,                      ONLY: build_2c_integrals,&
                                              build_2c_neighbor_lists
   USE rpa_gw_kpoints_util,             ONLY: cfm_screened_interaction_from_eps_cholesky
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_compute_V_and_W'
   PUBLIC :: compute_W, &
             get_W_MIC, &
             compute_V_k_by_lattice_sum, &
             compute_MinvVsqrt_Vsqrt, &
             compute_fm_W_MIC_freq_j_from_chi, &
             apply_extrapol_factor, &
             fm_read, &
             fm_write, &
             get_V_tr_R

CONTAINS

! **************************************************************************************************
!> \brief Computes the screened Coulomb interaction on the imaginary-time grid, entirely in the
!>        RI auxiliary (PQ) basis:
!>          χ_PQ(iω) = Σ_τ w(ω,τ) cos(ωτ) χ_PQ(iτ)                  (cosine transform)
!>          ε(iω)    = Id - V^0.5 M^-1 χ(iω) M^-1 V^0.5             (dielectric function)
!>          W(iω)    = V^0.5 ( ε^-1(iω) - Id ) V^0.5                (correlation part only)
!>          W(iτ)    = Σ_ω w̃(τ,ω) cos(ωτ) W(iω)                    (back transform)
!>          W(iτ)   <- M^-1 W(iτ) M^-1                              (fold in the RI metric)
!>        where V is the bare Coulomb matrix and M the RI metric.
!> \param bs_env ...
!> \param qs_env ...
!> \param mat_chi_Gamma_tau ...
!> \param fm_W_time ...
! **************************************************************************************************
   SUBROUTINE compute_W(bs_env, qs_env, mat_chi_Gamma_tau, fm_W_time)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_chi_Gamma_tau
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_W_time

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_W'

      INTEGER                                            :: handle, i_t, j_w
      REAL(KIND=dp)                                      :: t1
      TYPE(cp_fm_type)                                   :: fm_M_inv_V_sqrt, fm_V, fm_V_sqrt

      CALL timeset(routineN, handle)

      t1 = m_walltime()

      CALL create_fm_W_MIC_time(bs_env, fm_W_time)

      ! 1. Allocate V and M matrices
      CALL cp_fm_create(fm_V, bs_env%fm_RI_RI%matrix_struct)
      CALL cp_fm_create(fm_V_sqrt, bs_env%fm_RI_RI%matrix_struct)
      CALL cp_fm_create(fm_M_inv_V_sqrt, bs_env%fm_RI_RI%matrix_struct)

      ! Compute V and M^-1 * V^0.5
      CALL compute_V_MinvVsqrt(bs_env, qs_env, fm_V, fm_V_sqrt, fm_M_inv_V_sqrt)

      ! 2. Loop over frequencies
      DO j_w = 1, bs_env%num_time_freq_points
         ! Fourier transformation of χ_PQ(iτ) to χ_PQ(iω_j)
         CALL compute_fm_chi_Gamma_freq(bs_env, bs_env%fm_chi_Gamma_freq, j_w, mat_chi_Gamma_tau)

         ! ε(iω_j) = Id - V^0.5*M^-1*χ(iω_j)*M^-1*V^0.5
         ! W(iω_j) = V^0.5*(ε^-1(iω_j)-Id)*V^0.5
         CALL compute_fm_W_freq(bs_env, bs_env%fm_chi_Gamma_freq, fm_V_sqrt, &
                                fm_M_inv_V_sqrt, bs_env%fm_W_MIC_freq)

         ! Fourier transform from W_PQ^MIC(iω_j) to W_PQ^MIC(iτ)
         CALL Fourier_transform_w_to_t(bs_env, fm_W_time, bs_env%fm_W_MIC_freq, j_w)
      END DO

      ! M^-1(k=0) W^MIC(iτ) M^-1(k=0) -> fm_W_time
      CALL fm_contract_ABA(bs_env%fm_Minv_Gamma, fm_W_time)

      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,T58,A,F7.1,A)') &
            'Computed W(iτ),', ' Execution time', m_walltime() - t1, ' s'
      END IF

      CALL dbcsr_deallocate_matrix_set(mat_chi_Gamma_tau)

      ! Cleanup
      CALL cp_fm_release(fm_V)
      CALL cp_fm_release(fm_V_sqrt)
      CALL cp_fm_release(fm_M_inv_V_sqrt)

      ! Marek : Fourier transform W^MIC(itau) back to get it at a specific im.frequency point - iomega = 0
      IF (bs_env%rtp_method == rtp_method_bse) THEN
         t1 = m_walltime()
         CALL cp_fm_create(bs_env%fm_W_MIC_freq_zero, bs_env%fm_W_MIC_freq%matrix_struct)
         ! Set to zero
         CALL cp_fm_set_all(bs_env%fm_W_MIC_freq_zero, 0.0_dp)
         ! Sum over all times
         DO i_t = 1, bs_env%num_time_freq_points
            ! Add the relevant structure with correct weight
            CALL cp_fm_scale_and_add(1.0_dp, bs_env%fm_W_MIC_freq_zero, &
                                     bs_env%time_frequency_grid%time_weights_at_zero_frequency(i_t), fm_W_time(i_t))
         END DO
         ! Done, save to file
         CALL fm_write(bs_env%fm_W_MIC_freq_zero, 0, "W_freq_rtp", qs_env)
         ! Report calculation
         IF (bs_env%unit_nr > 0) THEN
            WRITE (bs_env%unit_nr, '(T2,A,T57,A,F7.1,A)') &
               'Computed W(0),', ' Execution time', m_walltime() - t1, ' s'
         END IF
      END IF

      IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(A)') ' '

      CALL timestop(handle)

   END SUBROUTINE compute_W

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param qs_env ...
!> \param mat_chi_Gamma_tau ...
!> \param fm_W_MIC_time ...
! **************************************************************************************************
   SUBROUTINE get_W_MIC(bs_env, qs_env, mat_chi_Gamma_tau, fm_W_MIC_time)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_chi_Gamma_tau
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_W_MIC_time

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'get_W_MIC'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (bs_env%all_W_exist) THEN
         CALL read_W_MIC_time(bs_env, mat_chi_Gamma_tau, fm_W_MIC_time)
      ELSE
         CALL compute_W_MIC(bs_env, qs_env, mat_chi_Gamma_tau, fm_W_MIC_time)
      END IF

      CALL timestop(handle)

   END SUBROUTINE get_W_MIC

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param qs_env ...
!> \param mat_chi_Gamma_tau ...
!> \param fm_W_MIC_time ...
! **************************************************************************************************
   SUBROUTINE compute_W_MIC(bs_env, qs_env, mat_chi_Gamma_tau, fm_W_MIC_time)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_chi_Gamma_tau
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_W_MIC_time

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_W_MIC'

      INTEGER                                            :: handle, i_t, ikp, ikp_batch, &
                                                            ikp_in_batch, j_w
      REAL(KIND=dp)                                      :: t1
      TYPE(cp_cfm_type)                                  :: cfm_M_inv_V_sqrt_ikp, cfm_V_sqrt_ikp
      TYPE(cp_cfm_type), ALLOCATABLE, DIMENSION(:)       :: cfm_V_kp

      CALL timeset(routineN, handle)

      CALL create_fm_W_MIC_time(bs_env, fm_W_MIC_time)

      DO ikp_batch = 1, bs_env%num_chi_eps_W_batches

         t1 = m_walltime()

         ! Compute V_PQ(k) = sum_R e^(ikR) <phi_P, cell 0 | 1/r | phi_Q, cell R>
         CALL compute_V_k_by_lattice_sum(bs_env, qs_env, cfm_V_kp, ikp_batch)

         DO ikp_in_batch = 1, bs_env%nkp_chi_eps_W_batch

            ikp = (ikp_batch - 1)*bs_env%nkp_chi_eps_W_batch + ikp_in_batch

            IF (ikp > bs_env%nkp_chi_eps_W_orig_plus_extra) THEN
               CYCLE
            END IF

            CALL compute_MinvVsqrt_Vsqrt(bs_env, qs_env, cfm_V_kp, &
                                         cfm_V_sqrt_ikp, cfm_M_inv_V_sqrt_ikp, ikp)

            CALL bs_env%para_env%sync()
            CALL cp_cfm_release(cfm_V_kp(ikp))

            DO j_w = 1, bs_env%num_time_freq_points

               ! check if we need this (ikp, ω_j) combination for approximate k-point extrapolation
               IF (bs_env%approx_kp_extrapol .AND. j_w > 1 .AND. &
                   ikp > bs_env%nkp_chi_eps_W_orig) CYCLE

               CALL compute_fm_W_MIC_freq_j(bs_env, qs_env, bs_env%fm_W_MIC_freq, j_w, ikp, &
                                            mat_chi_Gamma_tau, cfm_M_inv_V_sqrt_ikp, &
                                            cfm_V_sqrt_ikp)

               ! Fourier trafo from W_PQ^MIC(iω_j) to W_PQ^MIC(iτ)
               CALL Fourier_transform_w_to_t(bs_env, fm_W_MIC_time, bs_env%fm_W_MIC_freq, j_w)

            END DO ! ω_j

         END DO ! ikp_in_batch

         DEALLOCATE (cfm_V_kp)

         IF (bs_env%unit_nr > 0) THEN
            WRITE (bs_env%unit_nr, '(T2,A,I12,A,I3,A,F10.1,A)') &
               'Computed W(iτ,k) for k-point batch', &
               ikp_batch, ' /', bs_env%num_chi_eps_W_batches, &
               ', Execution time', m_walltime() - t1, ' s'
         END IF

      END DO ! ikp_batch

      IF (bs_env%approx_kp_extrapol) THEN
         CALL apply_extrapol_factor(bs_env, fm_W_MIC_time)
      END IF

      ! M^-1(k=0) W^MIC(iτ) M^-1(k=0) -> fm_W_MIC_time
      CALL fm_contract_ABA(bs_env%fm_Minv_Gamma, fm_W_MIC_time)

      DO i_t = 1, bs_env%num_time_freq_points
         CALL fm_write(fm_W_MIC_time(i_t), i_t, bs_env%W_time_name, qs_env)
      END DO

      CALL cp_cfm_release(cfm_M_inv_V_sqrt_ikp)
      CALL cp_cfm_release(cfm_V_sqrt_ikp)
      CALL dbcsr_deallocate_matrix_set(mat_chi_Gamma_tau)

      ! W^MIC(iω=0) from W^MIC(iτ) for the post-GW consumers; the file is the RT-BSE restart's
      IF (bs_env%rtp_method == rtp_method_bse .OR. &
          bs_env%rtp_method == rtp_method_bse_linearized) THEN
         CALL sum_W_MIC_freq_zero(bs_env, fm_W_MIC_time)
         CALL fm_write(bs_env%fm_W_MIC_freq_zero, 0, "W_freq_rtp", qs_env)
      ELSE IF (bs_env%keep_W_freq_zero) THEN
         CALL sum_W_MIC_freq_zero(bs_env, fm_W_MIC_time)
      END IF

      IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(A)') ' '

      CALL timestop(handle)

   END SUBROUTINE compute_W_MIC

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param qs_env ...
!> \param cfm_V_kp ...
!> \param ikp_batch ...
! **************************************************************************************************
   SUBROUTINE compute_V_k_by_lattice_sum(bs_env, qs_env, cfm_V_kp, ikp_batch)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_cfm_type), ALLOCATABLE, DIMENSION(:)       :: cfm_V_kp
      INTEGER                                            :: ikp_batch

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_V_k_by_lattice_sum'

      INTEGER                                            :: handle, ikp, ikp_end, ikp_start, &
                                                            nkp_chi_eps_W_batch, re_im
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_type)                                   :: fm_V_kp_work
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_V_kp
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      nkp_chi_eps_W_batch = bs_env%nkp_chi_eps_W_batch

      ikp_start = (ikp_batch - 1)*bs_env%nkp_chi_eps_W_batch + 1
      ikp_end = MIN(ikp_batch*bs_env%nkp_chi_eps_W_batch, bs_env%kpoints_chi_eps_W%nkp)

      NULLIFY (mat_V_kp)
      ALLOCATE (mat_V_kp(ikp_start:ikp_end, 2))

      DO re_im = 1, 2
         DO ikp = ikp_start, ikp_end
            NULLIFY (mat_V_kp(ikp, re_im)%matrix)
            ALLOCATE (mat_V_kp(ikp, re_im)%matrix)
            CALL dbcsr_create(mat_V_kp(ikp, re_im)%matrix, template=bs_env%mat_RI_RI%matrix)
            CALL dbcsr_reserve_all_blocks(mat_V_kp(ikp, re_im)%matrix)
            CALL dbcsr_set(mat_V_kp(ikp, re_im)%matrix, 0.0_dp)
         END DO ! ikp
      END DO ! re_im

      CALL get_qs_env(qs_env=qs_env, &
                      particle_set=particle_set, &
                      cell=cell, &
                      qs_kind_set=qs_kind_set, &
                      atomic_kind_set=atomic_kind_set)

      IF (ikp_end <= bs_env%nkp_chi_eps_W_orig) THEN

         ! 1. 2c Coulomb integrals for the first "original" k-point grid
         bs_env%kpoints_chi_eps_W%nkp_grid = bs_env%nkp_grid_chi_eps_W_orig

      ELSE IF (ikp_start > bs_env%nkp_chi_eps_W_orig .AND. &
               ikp_end <= bs_env%nkp_chi_eps_W_orig_plus_extra) THEN

         ! 2. 2c Coulomb integrals for the second "extrapolation" k-point grid
         bs_env%kpoints_chi_eps_W%nkp_grid = bs_env%nkp_grid_chi_eps_W_extra

      ELSE

         CPABORT("Error with k-point parallelization.")

      END IF

      CALL build_2c_coulomb_matrix_kp(mat_V_kp, &
                                      bs_env%kpoints_chi_eps_W, &
                                      basis_type="RI_AUX", &
                                      cell=cell, &
                                      particle_set=particle_set, &
                                      qs_kind_set=qs_kind_set, &
                                      atomic_kind_set=atomic_kind_set, &
                                      size_lattice_sum=bs_env%size_lattice_sum_V, &
                                      operator_type=operator_coulomb, &
                                      ikp_start=ikp_start, &
                                      ikp_end=ikp_end)

      bs_env%kpoints_chi_eps_W%nkp_grid = bs_env%nkp_grid_chi_eps_W_orig

      ALLOCATE (cfm_V_kp(ikp_start:ikp_end))
      CALL cp_fm_create(fm_V_kp_work, bs_env%fm_RI_RI%matrix_struct)
      DO ikp = ikp_start, ikp_end
         CALL cp_cfm_create(cfm_V_kp(ikp), bs_env%fm_RI_RI%matrix_struct)
         CALL cp_fm_set_all(fm_V_kp_work, 0.0_dp)
         CALL copy_dbcsr_to_fm(mat_V_kp(ikp, 1)%matrix, fm_V_kp_work)
         cfm_V_kp(ikp)%local_data = CMPLX(fm_V_kp_work%local_data, 0.0_dp, KIND=dp)

         CALL cp_fm_set_all(fm_V_kp_work, 0.0_dp)
         CALL copy_dbcsr_to_fm(mat_V_kp(ikp, 2)%matrix, fm_V_kp_work)
         cfm_V_kp(ikp)%local_data = CMPLX(REAL(cfm_V_kp(ikp)%local_data, KIND=dp), &
                                          fm_V_kp_work%local_data, KIND=dp)
      END DO
      CALL cp_fm_release(fm_V_kp_work)
      CALL dbcsr_deallocate_matrix_set(mat_V_kp)

      CALL timestop(handle)

   END SUBROUTINE compute_V_k_by_lattice_sum

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param qs_env ...
!> \param cfm_V_kp ...
!> \param cfm_V_sqrt_ikp ...
!> \param cfm_M_inv_V_sqrt_ikp ...
!> \param ikp ...
! **************************************************************************************************
   SUBROUTINE compute_MinvVsqrt_Vsqrt(bs_env, qs_env, cfm_V_kp, cfm_V_sqrt_ikp, &
                                      cfm_M_inv_V_sqrt_ikp, ikp)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_cfm_type), ALLOCATABLE, DIMENSION(:)       :: cfm_V_kp
      TYPE(cp_cfm_type)                                  :: cfm_V_sqrt_ikp, cfm_M_inv_V_sqrt_ikp
      INTEGER                                            :: ikp

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_MinvVsqrt_Vsqrt'

      INTEGER                                            :: handle, info, n_RI
      TYPE(cp_cfm_type)                                  :: cfm_M_inv_ikp, cfm_work
      TYPE(cp_cfm_type), ALLOCATABLE, DIMENSION(:)       :: cfm_M_ikp

      CALL timeset(routineN, handle)

      n_RI = bs_env%n_RI

      ! get here M(k) and write it to cfm_M_ikp
      CALL RI_2c_integral_mat_kp(qs_env, cfm_M_ikp, cfm_V_kp(ikp)%matrix_struct, n_RI, bs_env%ri_metric, &
                                 bs_env%kpoints_chi_eps_W, regularization_RI=bs_env%regularization_RI, &
                                 ikp_ext=ikp, do_build_cell_index=(ikp == 1))

      IF (ikp == 1) THEN
         CALL cp_cfm_create(cfm_V_sqrt_ikp, cfm_V_kp(ikp)%matrix_struct)
         CALL cp_cfm_create(cfm_M_inv_V_sqrt_ikp, cfm_V_kp(ikp)%matrix_struct)
      END IF
      CALL cp_cfm_create(cfm_M_inv_ikp, cfm_V_kp(ikp)%matrix_struct)

      CALL cp_cfm_to_cfm(cfm_M_ikp(1), cfm_M_inv_ikp)
      CALL cp_cfm_to_cfm(cfm_V_kp(ikp), cfm_V_sqrt_ikp)

      CALL cp_cfm_release(cfm_M_ikp(1))
      DEALLOCATE (cfm_M_ikp)

      CALL cp_cfm_create(cfm_work, cfm_V_kp(ikp)%matrix_struct)

      ! M(k) -> M^-1(k)
      CALL cp_cfm_to_cfm(cfm_M_inv_ikp, cfm_work)
      CALL cp_cfm_cholesky_decompose(matrix=cfm_M_inv_ikp, n=n_RI, info_out=info)
      IF (info == 0) THEN
         ! successful Cholesky decomposition
         CALL cp_cfm_cholesky_invert(cfm_M_inv_ikp)
         ! symmetrize the result
         CALL cp_cfm_uplo_to_full(cfm_M_inv_ikp)
      ELSE
         ! Cholesky decomposition not successful: use expensive diagonalization
         CALL cp_cfm_power(cfm_work, threshold=bs_env%eps_eigval_mat_RI, exponent=-1.0_dp)
         CALL cp_cfm_to_cfm(cfm_work, cfm_M_inv_ikp)
      END IF

      ! V(k) -> L(k) with L^H(k)*L(k) = V(k) [L(k) can be just considered to be V^0.5(k)]
      CALL cp_cfm_to_cfm(cfm_V_sqrt_ikp, cfm_work)
      CALL cp_cfm_cholesky_decompose(matrix=cfm_V_sqrt_ikp, n=n_RI, info_out=info)
      IF (info == 0) THEN
         ! successful Cholesky decomposition
         CALL clean_lower_part(cfm_V_sqrt_ikp)
      ELSE
         ! Cholesky decomposition not successful: use expensive diagonalization
         CALL cp_cfm_power(cfm_work, threshold=0.0_dp, exponent=0.5_dp)
         CALL cp_cfm_to_cfm(cfm_work, cfm_V_sqrt_ikp)
      END IF
      CALL cp_cfm_release(cfm_work)

      ! get M^-1(k)*V^0.5(k)
      CALL parallel_gemm("N", "C", n_RI, n_RI, n_RI, z_one, cfm_M_inv_ikp, cfm_V_sqrt_ikp, &
                         z_zero, cfm_M_inv_V_sqrt_ikp)

      CALL cp_cfm_release(cfm_M_inv_ikp)

      CALL timestop(handle)

   END SUBROUTINE compute_MinvVsqrt_Vsqrt

! **************************************************************************************************
!> \brief Compute the screened interaction for one frequency and k-point.
!> \param bs_env band-structure environment
!> \param cfm_chi_eps_W_ikp_freq_j input chi; overwritten by W
!> \param cfm_V_sqrt_ikp square root of the Coulomb matrix
!> \param cfm_M_inv_V_sqrt_ikp metric-weighted Coulomb factor
! **************************************************************************************************
   SUBROUTINE compute_cfm_W_ikp_freq_j(bs_env, cfm_chi_eps_W_ikp_freq_j, cfm_V_sqrt_ikp, &
                                       cfm_M_inv_V_sqrt_ikp)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_cfm_type), INTENT(INOUT)                   :: cfm_chi_eps_W_ikp_freq_j
      TYPE(cp_cfm_type), INTENT(IN)                      :: cfm_V_sqrt_ikp, cfm_M_inv_V_sqrt_ikp

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_cfm_W_ikp_freq_j'

      INTEGER                                            :: handle, info

      CALL timeset(routineN, handle)

      ! 1. epsilon(iw,k) = Id - V^0.5(k)*M^-1(k)*chi(iw,k)*M^-1(k)*V^0.5(k)

      ! 1. a) Form V^0.5(k)*M^-1(k)*chi(iw,k)*M^-1(k)*V^0.5(k).
      CALL cfm_contract_ABA(cfm_M_inv_V_sqrt_ikp, cfm_chi_eps_W_ikp_freq_j)

      ! 1. b) Add the identity matrix.
      CALL cp_cfm_add_on_diag(cfm_chi_eps_W_ikp_freq_j, z_one)

      ! 2. W(iw,k) = V^0.5(k)*(epsilon^-1(iw,k)-Id)*V^0.5(k)

      ! 2. a) Cholesky-factorize epsilon before inversion.
      CALL cp_cfm_cholesky_decompose(matrix=cfm_chi_eps_W_ikp_freq_j, n=bs_env%n_RI, info_out=info)
      CPASSERT(info == 0)

      ! 2. b) Invert epsilon, subtract Id, then contract with V^0.5 on both sides.
      CALL cfm_screened_interaction_from_eps_cholesky(cfm_chi_eps_W_ikp_freq_j, cfm_V_sqrt_ikp)

      CALL timestop(handle)

   END SUBROUTINE compute_cfm_W_ikp_freq_j

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param mat_chi_Gamma_tau ...
!> \param fm_W_MIC_time ...
! **************************************************************************************************
   SUBROUTINE read_W_MIC_time(bs_env, mat_chi_Gamma_tau, fm_W_MIC_time)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_chi_Gamma_tau
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_W_MIC_time

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'read_W_MIC_time'

      INTEGER                                            :: handle, i_t
      REAL(KIND=dp)                                      :: t1

      CALL timeset(routineN, handle)

      CALL dbcsr_deallocate_matrix_set(mat_chi_Gamma_tau)
      CALL create_fm_W_MIC_time(bs_env, fm_W_MIC_time)

      DO i_t = 1, bs_env%num_time_freq_points

         t1 = m_walltime()

         CALL fm_read(fm_W_MIC_time(i_t), bs_env, bs_env%W_time_name, i_t)

         IF (bs_env%unit_nr > 0) THEN
            WRITE (bs_env%unit_nr, '(T2,A,I5,A,I3,A,F10.1,A)') &
               'Read W^MIC(iτ) from file for time point  ', i_t, ' /', bs_env%num_time_freq_points, &
               ', Execution time', m_walltime() - t1, ' s'
         END IF

      END DO

      IF (bs_env%unit_nr > 0) WRITE (bs_env%unit_nr, '(A)') ' '

      ! W^MIC(iω=0) for the post-GW consumers, built here because fm_W_MIC_time does not outlive
      ! the W stage: RT-BSE reads its own restart file, PROPERTIES%BSE sums the W(iτ) just read.
      IF (bs_env%rtp_method == rtp_method_bse .OR. &
          bs_env%rtp_method == rtp_method_bse_linearized) THEN
         CALL cp_fm_create(bs_env%fm_W_MIC_freq_zero, bs_env%fm_W_MIC_freq%matrix_struct)
         t1 = m_walltime()
         CALL fm_read(bs_env%fm_W_MIC_freq_zero, bs_env, "W_freq_rtp", 0)
         IF (bs_env%unit_nr > 0) THEN
            WRITE (bs_env%unit_nr, '(T2,A,I3,A,I3,A,F10.1,A)') &
               'Read W^MIC(f=0) from file for freq. point  ', 1, ' /', 1, &
               ', Execution time', m_walltime() - t1, ' s'
         END IF
      ELSE IF (bs_env%keep_W_freq_zero) THEN
         ! no RTP file in a LRBSE-only run: W(iτ) was just read from the W_time files
         CALL sum_W_MIC_freq_zero(bs_env, fm_W_MIC_time)
      END IF

      CALL timestop(handle)

   END SUBROUTINE read_W_MIC_time

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param qs_env ...
!> \param fm_W_MIC_freq_j ...
!> \param j_w ...
!> \param ikp ...
!> \param mat_chi_Gamma_tau ...
!> \param cfm_M_inv_V_sqrt_ikp ...
!> \param cfm_V_sqrt_ikp ...
! **************************************************************************************************
   SUBROUTINE compute_fm_W_MIC_freq_j(bs_env, qs_env, fm_W_MIC_freq_j, j_w, ikp, mat_chi_Gamma_tau, &
                                      cfm_M_inv_V_sqrt_ikp, cfm_V_sqrt_ikp)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type)                                   :: fm_W_MIC_freq_j
      INTEGER                                            :: j_w, ikp
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_chi_Gamma_tau
      TYPE(cp_cfm_type)                                  :: cfm_M_inv_V_sqrt_ikp, cfm_V_sqrt_ikp

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_fm_W_MIC_freq_j'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! 1. Fourier transformation of χ_PQ(iτ,k=0) to χ_PQ(iω_j,k=0)
      CALL compute_fm_chi_Gamma_freq(bs_env, bs_env%fm_chi_Gamma_freq, j_w, mat_chi_Gamma_tau)

      ! 2. χ_PQ(iω_j,k=0) -> W_PQ^MIC(iω_j)
      CALL compute_fm_W_MIC_freq_j_from_chi(bs_env, qs_env, fm_W_MIC_freq_j, j_w, ikp, &
                                            cfm_M_inv_V_sqrt_ikp, cfm_V_sqrt_ikp)

      CALL timestop(handle)

   END SUBROUTINE compute_fm_W_MIC_freq_j

! **************************************************************************************************
!> \brief W_PQ^MIC(iω_j) from χ_PQ(iω_j,k=0) in bs_env%fm_chi_Gamma_freq
!> \param bs_env ...
!> \param qs_env ...
!> \param fm_W_MIC_freq_j ...
!> \param j_w ...
!> \param ikp ...
!> \param cfm_M_inv_V_sqrt_ikp ...
!> \param cfm_V_sqrt_ikp ...
! **************************************************************************************************
   SUBROUTINE compute_fm_W_MIC_freq_j_from_chi(bs_env, qs_env, fm_W_MIC_freq_j, j_w, ikp, &
                                               cfm_M_inv_V_sqrt_ikp, cfm_V_sqrt_ikp)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type)                                   :: fm_W_MIC_freq_j
      INTEGER                                            :: j_w, ikp
      TYPE(cp_cfm_type)                                  :: cfm_M_inv_V_sqrt_ikp, cfm_V_sqrt_ikp

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_fm_W_MIC_freq_j_from_chi'

      INTEGER                                            :: handle
      TYPE(cp_cfm_type)                                  :: cfm_chi_eps_W_ikp_freq_j

      CALL timeset(routineN, handle)

      CALL cp_fm_set_all(fm_W_MIC_freq_j, 0.0_dp)

      ! 1. Get χ_PQ(iω_j,k_i) from χ_PQ(iω_j,k=0) using the minimum image convention
      CALL cfm_ikp_from_fm_Gamma(cfm_chi_eps_W_ikp_freq_j, bs_env%fm_chi_Gamma_freq, &
                                 ikp, qs_env, bs_env%kpoints_chi_eps_W, "RI_AUX")

      ! 2. Remove all negative eigenvalues from χ_PQ(iω_j,k_i)
      CALL cp_cfm_power(cfm_chi_eps_W_ikp_freq_j, threshold=0.0_dp, exponent=1.0_dp)

      ! 3. ε(iω_j,k_i) = Id - V^0.5(k_i)*M^-1(k_i)*χ(iω_j,k_i)*M^-1(k_i)*V^0.5(k_i)
      CALL compute_cfm_W_ikp_freq_j(bs_env, cfm_chi_eps_W_ikp_freq_j, cfm_V_sqrt_ikp, &
                                    cfm_M_inv_V_sqrt_ikp)

      ! 4. k-point integration W_PQ(iω_j, k_i) to W_PQ^MIC(iω_j)
      SELECT CASE (bs_env%approx_kp_extrapol)
      CASE (.FALSE.)
         ! default: standard k-point extrapolation
         CALL MIC_contribution_from_ikp(bs_env, qs_env, fm_W_MIC_freq_j, &
                                        cfm_chi_eps_W_ikp_freq_j, ikp, &
                                        bs_env%kpoints_chi_eps_W, "RI_AUX")
      CASE (.TRUE.)
         ! for approximate kpoint extrapolation: get W_PQ^MIC(iω_1) with and without k-point
         ! extrapolation to compute the extrapolation factor f_PQ for every PQ-matrix element,
         ! f_PQ = (W_PQ^MIC(iω_1) with extrapolation) / (W_PQ^MIC(iω_1) without extrapolation)

         ! for ω_1, we compute the k-point extrapolated result using all k-points
         IF (j_w == 1) THEN

            ! k-point extrapolated
            CALL MIC_contribution_from_ikp(bs_env, qs_env, bs_env%fm_W_MIC_freq_1_extra, &
                                           cfm_chi_eps_W_ikp_freq_j, ikp, &
                                           bs_env%kpoints_chi_eps_W, &
                                           "RI_AUX")
            ! non-kpoint extrapolated
            IF (ikp <= bs_env%nkp_chi_eps_W_orig) THEN
               CALL MIC_contribution_from_ikp(bs_env, qs_env, bs_env%fm_W_MIC_freq_1_no_extra, &
                                              cfm_chi_eps_W_ikp_freq_j, ikp, &
                                              bs_env%kpoints_chi_eps_W, &
                                              "RI_AUX", wkp_ext=bs_env%wkp_orig)
            END IF

         END IF

         ! for all ω_j, we need to compute W^MIC without k-point extrpolation
         IF (ikp <= bs_env%nkp_chi_eps_W_orig) THEN
            CALL MIC_contribution_from_ikp(bs_env, qs_env, fm_W_MIC_freq_j, &
                                           cfm_chi_eps_W_ikp_freq_j, &
                                           ikp, bs_env%kpoints_chi_eps_W, "RI_AUX", &
                                           wkp_ext=bs_env%wkp_orig)
         END IF
      END SELECT

      CALL cp_cfm_release(cfm_chi_eps_W_ikp_freq_j)

      CALL timestop(handle)

   END SUBROUTINE compute_fm_W_MIC_freq_j_from_chi

! **************************************************************************************************
!> \brief ...
!> \param cfm_mat ...
! **************************************************************************************************
   SUBROUTINE clean_lower_part(cfm_mat)
      TYPE(cp_cfm_type)                                  :: cfm_mat

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'clean_lower_part'

      INTEGER                                            :: handle, i_row, j_col, j_global, &
                                                            ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices

      CALL timeset(routineN, handle)

      CALL cp_cfm_get_info(matrix=cfm_mat, &
                           nrow_local=nrow_local, ncol_local=ncol_local, &
                           row_indices=row_indices, col_indices=col_indices)

      DO j_col = 1, ncol_local
         j_global = col_indices(j_col)
         DO i_row = 1, nrow_local
            IF (j_global < row_indices(i_row)) cfm_mat%local_data(i_row, j_col) = z_zero
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE clean_lower_part

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param fm_W_MIC_time ...
! **************************************************************************************************
   SUBROUTINE apply_extrapol_factor(bs_env, fm_W_MIC_time)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_W_MIC_time

      CHARACTER(LEN=*), PARAMETER :: routineN = 'apply_extrapol_factor'
      REAL(KIND=dp), PARAMETER                           :: eps_W_no_extra_1 = 1.0E-13_dp

      INTEGER                                            :: handle, i, i_t, j, ncol_local, nrow_local
      REAL(KIND=dp)                                      :: extrapol_factor, W_extra_1, W_no_extra_1

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(matrix=fm_W_MIC_time(1), nrow_local=nrow_local, ncol_local=ncol_local)

      DO i_t = 1, SIZE(fm_W_MIC_time)
         DO j = 1, ncol_local
            DO i = 1, nrow_local

               W_extra_1 = bs_env%fm_W_MIC_freq_1_extra%local_data(i, j)
               W_no_extra_1 = bs_env%fm_W_MIC_freq_1_no_extra%local_data(i, j)

               IF (ABS(W_no_extra_1) > eps_W_no_extra_1) THEN
                  extrapol_factor = ABS(W_extra_1/W_no_extra_1)
               ELSE
                  extrapol_factor = 1.0_dp
               END IF

               ! reset extrapolation factor if it is very large
               IF (extrapol_factor > 10.0_dp) extrapol_factor = 1.0_dp

               fm_W_MIC_time(i_t)%local_data(i, j) = fm_W_MIC_time(i_t)%local_data(i, j) &
                                                     *extrapol_factor
            END DO
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE apply_extrapol_factor

! **************************************************************************************************
!> \brief Computes V, V^0.5, and M^-1 V^0.5 for the RI-RS dielectric function.
!>        The Coulomb matrix V is constructed by the RI-RS k-point path. The inverse metric
!>        M^-1(k=0) is precomputed once in gw_utils and stored in bs_env.
!> \param bs_env ...
!> \param qs_env ...
!> \param fm_V Coulomb matrix V(k=0)
!> \param fm_V_sqrt symmetric factor V^0.5
!> \param fm_Minv_Vsqrt product M^-1(k=0) V^0.5
! **************************************************************************************************
   SUBROUTINE compute_V_MinvVsqrt(bs_env, qs_env, fm_V, fm_V_sqrt, fm_Minv_Vsqrt)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_V, fm_V_sqrt, fm_Minv_Vsqrt

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_V_MinvVsqrt'

      INTEGER                                            :: handle
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_V_kp
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      IF (bs_env%auto_ri%enabled) THEN
         ! -------------------------------------------------------------------
         ! 1a. The optimized AB functions span two atoms, so their previously transformed
         !     Coulomb matrix cannot be rebuilt by the atom-local k-point integral routine.
         ! -------------------------------------------------------------------
         CALL cp_fm_to_fm(bs_env%auto_ri%V_pq, fm_V)
      ELSE
         ! -------------------------------------------------------------------
         ! 1b. Build Coulomb Matrix V(k=0) using the kp-routine but only for ikp=1.
         ! -------------------------------------------------------------------
         CALL get_qs_env(qs_env=qs_env, particle_set=particle_set, cell=cell, &
                         qs_kind_set=qs_kind_set)
         atomic_kind_set => bs_env%ri_rs%atomic_kind_set

         ALLOCATE (mat_V_kp(1:1, 1:2))
         NULLIFY (mat_V_kp(1, 1)%matrix, mat_V_kp(1, 2)%matrix)
         ALLOCATE (mat_V_kp(1, 1)%matrix, mat_V_kp(1, 2)%matrix)
         CALL dbcsr_create(mat_V_kp(1, 1)%matrix, template=bs_env%mat_RI_RI%matrix)
         CALL dbcsr_reserve_all_blocks(mat_V_kp(1, 1)%matrix)
         CALL dbcsr_set(mat_V_kp(1, 1)%matrix, 0.0_dp)
         CALL dbcsr_create(mat_V_kp(1, 2)%matrix, template=bs_env%mat_RI_RI%matrix)
         CALL dbcsr_reserve_all_blocks(mat_V_kp(1, 2)%matrix)
         ! The dummy imaginary part is required only by the k-point routine interface.
         CALL dbcsr_set(mat_V_kp(1, 2)%matrix, 0.0_dp)

         bs_env%kpoints_chi_eps_W%nkp_grid = bs_env%nkp_grid_chi_eps_W_orig
         CALL build_2c_coulomb_matrix_kp(mat_V_kp, bs_env%kpoints_chi_eps_W, "RI_AUX", cell, &
                                         particle_set, qs_kind_set, atomic_kind_set, &
                                         bs_env%size_lattice_sum_V, operator_coulomb, 1, 1)
         CALL copy_dbcsr_to_fm(mat_V_kp(1, 1)%matrix, fm_V)

         CALL dbcsr_deallocate_matrix(mat_V_kp(1, 1)%matrix)
         CALL dbcsr_deallocate_matrix(mat_V_kp(1, 2)%matrix)
         DEALLOCATE (mat_V_kp)
      END IF

      ! -----------------------------------------------------------------------
      ! 2. V -> V^0.5.
      ! -----------------------------------------------------------------------
      CALL fm_sqrt(fm_V, fm_V_sqrt, bs_env%eps_eigval_mat_RI, bs_env%unit_nr)

      ! -----------------------------------------------------------------------
      ! 3. M^-1(k=0) V^0.5.
      ! -----------------------------------------------------------------------
      CALL parallel_gemm("N", "T", bs_env%n_RI, bs_env%n_RI, bs_env%n_RI, 1.0_dp, &
                         bs_env%fm_Minv_Gamma, fm_V_sqrt, 0.0_dp, fm_Minv_Vsqrt)

      CALL timestop(handle)

   END SUBROUTINE compute_V_MinvVsqrt

! **************************************************************************************************
!> \brief Computes the screened interaction at one imaginary frequency:
!>          ε(iω_j)  = Id - (M^-1 V^0.5)^T χ(iω_j) (M^-1 V^0.5)
!>          W(iω_j)  = V^0.5^T ( ε^-1(iω_j) - Id ) V^0.5
!>        ε is inverted via Cholesky; if that fails due to conditioning, via
!>        eigendecomposition (cp_fm_power) with eigenvalue filtering.
!> \param bs_env ...
!> \param fm_chi_freq_j ...
!> \param fm_V_sqrt ...
!> \param fm_Minv_Vsqrt ...
!> \param fm_W_freq_j ...
! **************************************************************************************************
   SUBROUTINE compute_fm_W_freq(bs_env, fm_chi_freq_j, fm_V_sqrt, fm_Minv_Vsqrt, fm_W_freq_j)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_chi_freq_j, fm_V_sqrt, fm_Minv_Vsqrt
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_W_freq_j

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'compute_fm_W_freq'

      INTEGER                                            :: handle, n_RI
      TYPE(cp_fm_type)                                   :: fm_eps_freq_j, fm_work

      CALL timeset(routineN, handle)

      n_RI = bs_env%n_RI

      CALL cp_fm_create(fm_eps_freq_j, fm_chi_freq_j%matrix_struct)
      CALL cp_fm_create(fm_work, fm_chi_freq_j%matrix_struct)

      ! -----------------------------------------------------------------------
      ! 1. ε(iω_j) = Id - (M^-1 * V^0.5)^T * χ(iω_j) * (M^-1 * V^0.5)
      ! -----------------------------------------------------------------------
      ! work = χ(iω_j) * (M^-1 * V^0.5)
      CALL parallel_gemm('N', 'N', n_RI, n_RI, n_RI, 1.0_dp, &
                         fm_chi_freq_j, fm_Minv_Vsqrt, 0.0_dp, fm_work)

      ! eps_work = (M^-1 * V^0.5)^T * work
      CALL parallel_gemm('T', 'N', n_RI, n_RI, n_RI, 1.0_dp, &
                         fm_Minv_Vsqrt, fm_work, 0.0_dp, fm_eps_freq_j)

      ! ε(iω_j) = Id - eps_work  -->  -eps_work + Id
      CALL cp_fm_add_on_diag(fm_eps_freq_j, 1.0_dp)

      ! Force perfect symmetry before Cholesky to avoid info != 0 due to GEMM noise
      CALL cp_fm_uplo_to_full(fm_eps_freq_j, fm_work)

      ! -----------------------------------------------------------------------
      ! 2. W(iω_j) = V^0.5^T * (ε^-1(iω_j) - Id) * V^0.5
      ! -----------------------------------------------------------------------

      ! a) Invert ε by Cholesky decomposition or, if that fails, by diagonalization.
      CALL fm_invert(fm_eps_freq_j, bs_env%eps_eigval_mat_RI, bs_env%unit_nr)

      ! b) ε^-1(iω_j) - Id
      CALL cp_fm_add_on_diag(fm_eps_freq_j, -1.0_dp)

      ! c) W(iω_j) = V^0.5^T * (ε^-1(iω_j) - Id) * V^0.5
      CALL fm_contract_ABA(fm_V_sqrt, fm_eps_freq_j, fm_W_freq_j)

      ! Cleanup
      CALL cp_fm_release(fm_work)
      CALL cp_fm_release(fm_eps_freq_j)

      CALL timestop(handle)

   END SUBROUTINE compute_fm_W_freq

! **************************************************************************************************
!> \brief ...
!> \param fm ...
!> \param bs_env ...
!> \param mat_name ...
!> \param idx ...
! **************************************************************************************************
   SUBROUTINE fm_read(fm, bs_env, mat_name, idx)
      TYPE(cp_fm_type)                                   :: fm
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      CHARACTER(LEN=*)                                   :: mat_name
      INTEGER                                            :: idx

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'fm_read'

      CHARACTER(LEN=default_path_length)                 :: f_chi
      INTEGER                                            :: handle, unit_nr

      CALL timeset(routineN, handle)

      unit_nr = -1
      IF (bs_env%para_env%is_source()) THEN

         IF (idx < 10) THEN
            WRITE (f_chi, '(3A,I1,A)') TRIM(bs_env%prefix), TRIM(mat_name), "_0", idx, ".matrix"
         ELSE IF (idx < 100) THEN
            WRITE (f_chi, '(3A,I2,A)') TRIM(bs_env%prefix), TRIM(mat_name), "_", idx, ".matrix"
         ELSE
            CPABORT('Please implement more than 99 time/frequency points.')
         END IF

         CALL open_file(file_name=TRIM(f_chi), file_action="READ", file_form="UNFORMATTED", &
                        file_position="REWIND", file_status="OLD", unit_number=unit_nr)

      END IF

      CALL cp_fm_read_unformatted(fm, unit_nr)

      IF (bs_env%para_env%is_source()) CALL close_file(unit_number=unit_nr)

      CALL timestop(handle)

   END SUBROUTINE fm_read

! **************************************************************************************************
!> \brief Write a GW restart matrix when requested by the print key.
!> \param fm Matrix to write, or reusable conversion buffer when matrix is supplied.
!> \param matrix_index ...
!> \param matrix_name ...
!> \param qs_env ...
!> \param matrix Optional DBCSR source, converted to fm only when output is enabled.
! **************************************************************************************************
   SUBROUTINE fm_write(fm, matrix_index, matrix_name, qs_env, matrix)
      TYPE(cp_fm_type)                                   :: fm
      INTEGER                                            :: matrix_index
      CHARACTER(LEN=*)                                   :: matrix_name
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_type), INTENT(IN), OPTIONAL             :: matrix

      CHARACTER(LEN=*), PARAMETER :: key = 'PROPERTIES%BANDSTRUCTURE%GW%PRINT%RESTART', &
         routineN = 'fm_write'

      CHARACTER(LEN=default_path_length)                 :: filename
      INTEGER                                            :: handle, unit_nr
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: input

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, input=input)

      logger => cp_get_default_logger()

      IF (BTEST(cp_print_key_should_output(logger%iter_info, input, key), cp_p_file)) THEN

         IF (PRESENT(matrix)) THEN
            CALL cp_fm_set_all(fm, 0.0_dp)
            CALL copy_dbcsr_to_fm(matrix, fm)
         END IF

         IF (matrix_index < 10) THEN
            WRITE (filename, '(3A,I1)') "RESTART_", matrix_name, "_0", matrix_index
         ELSE IF (matrix_index < 100) THEN
            WRITE (filename, '(3A,I2)') "RESTART_", matrix_name, "_", matrix_index
         ELSE
            CPABORT('Please implement more than 99 time/frequency points.')
         END IF

         unit_nr = cp_print_key_unit_nr(logger, input, key, extension=".matrix", &
                                        file_form="UNFORMATTED", middle_name=TRIM(filename), &
                                        file_position="REWIND", file_action="WRITE")

         CALL cp_fm_write_unformatted(fm, unit_nr)
         IF (unit_nr > 0) THEN
            CALL close_file(unit_nr)
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE fm_write

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param fm_chi_Gamma_freq ...
!> \param j_w ...
!> \param mat_chi_Gamma_tau ...
! **************************************************************************************************
   SUBROUTINE compute_fm_chi_Gamma_freq(bs_env, fm_chi_Gamma_freq, j_w, mat_chi_Gamma_tau)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_fm_type)                                   :: fm_chi_Gamma_freq
      INTEGER                                            :: j_w
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_chi_Gamma_tau

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_fm_chi_Gamma_freq'

      INTEGER                                            :: handle, i_t
      REAL(KIND=dp)                                      :: weight_ij

      CALL timeset(routineN, handle)

      CALL dbcsr_set(bs_env%mat_RI_RI%matrix, 0.0_dp)

      DO i_t = 1, bs_env%num_time_freq_points

         weight_ij = bs_env%time_frequency_grid%cosine_time_to_frequency(j_w, i_t)

         ! Fourier transform
         CALL dbcsr_add(bs_env%mat_RI_RI%matrix, mat_chi_Gamma_tau(i_t)%matrix, &
                        1.0_dp, weight_ij)

      END DO

      CALL copy_dbcsr_to_fm(bs_env%mat_RI_RI%matrix, fm_chi_Gamma_freq)

      CALL timestop(handle)

   END SUBROUTINE compute_fm_chi_Gamma_freq

! **************************************************************************************************
!> \brief W^MIC_PQ(iω=0) = sum_j w_j W^MIC_PQ(iτ_j) with the time weights w_j of the grid at zero
!>        frequency, into bs_env%fm_W_MIC_freq_zero
!> \param bs_env time grid and its weights; receives fm_W_MIC_freq_zero
!> \param fm_W_MIC_time W^MIC_PQ(iτ_j) for every time point j
! **************************************************************************************************
   SUBROUTINE sum_W_MIC_freq_zero(bs_env, fm_W_MIC_time)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_fm_type), DIMENSION(:)                     :: fm_W_MIC_time

      CHARACTER(LEN=*), PARAMETER :: routineN = 'sum_W_MIC_freq_zero'

      INTEGER                                            :: handle, i_t
      REAL(KIND=dp)                                      :: t_start

      CALL timeset(routineN, handle)

      t_start = m_walltime()
      CALL cp_fm_create(bs_env%fm_W_MIC_freq_zero, bs_env%fm_W_MIC_freq%matrix_struct)
      CALL cp_fm_set_all(bs_env%fm_W_MIC_freq_zero, 0.0_dp)
      DO i_t = 1, bs_env%num_time_freq_points
         CALL cp_fm_scale_and_add(1.0_dp, bs_env%fm_W_MIC_freq_zero, &
                                  bs_env%time_frequency_grid%time_weights_at_zero_frequency(i_t), &
                                  fm_W_MIC_time(i_t))
      END DO
      IF (bs_env%unit_nr > 0) THEN
         WRITE (bs_env%unit_nr, '(T2,A,I11,A,I3,A,F10.1,A)') &
            'Computed W(ω=0,k) for k-point batch', &
            1, ' /', 1, &
            ', Execution time', m_walltime() - t_start, ' s'
      END IF

      CALL timestop(handle)

   END SUBROUTINE sum_W_MIC_freq_zero

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param fm_W_MIC_time ...
! **************************************************************************************************
   SUBROUTINE create_fm_W_MIC_time(bs_env, fm_W_MIC_time)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_W_MIC_time

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_fm_W_MIC_time'

      INTEGER                                            :: handle, i_t

      CALL timeset(routineN, handle)

      ALLOCATE (fm_W_MIC_time(bs_env%num_time_freq_points))
      DO i_t = 1, bs_env%num_time_freq_points
         CALL cp_fm_create(fm_W_MIC_time(i_t), bs_env%fm_RI_RI%matrix_struct, set_zero=.TRUE.)
      END DO

      CALL timestop(handle)

   END SUBROUTINE create_fm_W_MIC_time

! **************************************************************************************************
!> \brief ...
!> \param bs_env ...
!> \param fm_W_MIC_time ...
!> \param fm_W_MIC_freq_j ...
!> \param j_w ...
! **************************************************************************************************
   SUBROUTINE Fourier_transform_w_to_t(bs_env, fm_W_MIC_time, fm_W_MIC_freq_j, j_w)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_W_MIC_time
      TYPE(cp_fm_type)                                   :: fm_W_MIC_freq_j
      INTEGER                                            :: j_w

      CHARACTER(LEN=*), PARAMETER :: routineN = 'Fourier_transform_w_to_t'

      INTEGER                                            :: handle, i_t
      REAL(KIND=dp)                                      :: weight_ij

      CALL timeset(routineN, handle)

      DO i_t = 1, bs_env%num_time_freq_points

         weight_ij = bs_env%time_frequency_grid%cosine_frequency_to_time(i_t, j_w)

         ! actual Fourier transform
         CALL cp_fm_scale_and_add(alpha=1.0_dp, matrix_a=fm_W_MIC_time(i_t), &
                                  beta=weight_ij, matrix_b=fm_W_MIC_freq_j)

      END DO

      CALL timestop(handle)

   END SUBROUTINE Fourier_transform_w_to_t

! **************************************************************************************************
!> \brief ...
!> \param V_tr_R ...
!> \param pot_type ...
!> \param regularization_RI ...
!> \param bs_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE get_V_tr_R(V_tr_R, pot_type, regularization_RI, bs_env, qs_env)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: V_tr_R
      TYPE(libint_potential_type)                        :: pot_type
      REAL(KIND=dp)                                      :: regularization_RI
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'get_V_tr_R'

      INTEGER                                            :: handle, img, nimages_scf_desymm
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: sizes_RI
      INTEGER, DIMENSION(:), POINTER                     :: col_bsize, row_bsize
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_V_tr_R
      TYPE(dbcsr_distribution_type)                      :: dbcsr_dist
      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:)        :: mat_V_tr_R
      TYPE(distribution_2d_type), POINTER                :: dist_2d
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_RI
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      NULLIFY (sab_RI, dist_2d)

      CALL get_qs_env(qs_env=qs_env, &
                      blacs_env=blacs_env, &
                      distribution_2d=dist_2d, &
                      qs_kind_set=qs_kind_set, &
                      particle_set=particle_set)

      ALLOCATE (sizes_RI(bs_env%n_atom))
      CALL get_particle_set(particle_set, qs_kind_set, nsgf=sizes_RI, basis=bs_env%basis_set_RI)
      CALL build_2c_neighbor_lists(sab_RI, bs_env%basis_set_RI, bs_env%basis_set_RI, &
                                   pot_type, "2c_nl_RI", qs_env, sym_ij=.FALSE., &
                                   dist_2d=dist_2d)
      CALL cp_dbcsr_dist2d_to_dist(dist_2d, dbcsr_dist)
      ALLOCATE (row_bsize(SIZE(sizes_RI)))
      ALLOCATE (col_bsize(SIZE(sizes_RI)))
      row_bsize(:) = sizes_RI
      col_bsize(:) = sizes_RI

      nimages_scf_desymm = bs_env%nimages_scf_desymm
      ALLOCATE (mat_V_tr_R(nimages_scf_desymm))
      CALL dbcsr_create(mat_V_tr_R(1), "(RI|RI)", dbcsr_dist, dbcsr_type_no_symmetry, &
                        row_bsize, col_bsize)
      DEALLOCATE (row_bsize, col_bsize)

      DO img = 2, nimages_scf_desymm
         CALL dbcsr_create(mat_V_tr_R(img), template=mat_V_tr_R(1))
      END DO

      CALL build_2c_integrals(mat_V_tr_R, 0.0_dp, qs_env, sab_RI, bs_env%basis_set_RI, &
                              bs_env%basis_set_RI, pot_type, do_kpoints=.TRUE., &
                              ext_kpoints=bs_env%kpoints_scf_desymm, &
                              regularization_RI=regularization_RI)

      ALLOCATE (fm_V_tr_R(nimages_scf_desymm))
      DO img = 1, nimages_scf_desymm
         CALL cp_fm_create(fm_V_tr_R(img), bs_env%fm_RI_RI%matrix_struct)
         CALL copy_dbcsr_to_fm(mat_V_tr_R(img), fm_V_tr_R(img))
         CALL dbcsr_release(mat_V_tr_R(img))
      END DO

      IF (.NOT. ALLOCATED(V_tr_R)) THEN
         ALLOCATE (V_tr_R(bs_env%n_RI, bs_env%n_RI, nimages_scf_desymm))
      END IF

      CALL cp_fm_to_replicated_array(fm_V_tr_R, V_tr_R)

      CALL cp_fm_release(fm_V_tr_R)
      CALL dbcsr_distribution_release(dbcsr_dist)
      CALL release_neighbor_list_sets(sab_RI)

      CALL timestop(handle)

   END SUBROUTINE get_V_tr_R

END MODULE gw_compute_V_and_W
