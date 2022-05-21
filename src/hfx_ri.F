!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief RI-methods for HFX
! **************************************************************************************************

MODULE hfx_ri

   USE arnoldi_api,                     ONLY: arnoldi_extremal
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              real_to_scaled
   USE cp_array_utils,                  ONLY: cp_1d_r_p_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_cholesky,               ONLY: cp_dbcsr_cholesky_decompose,&
                                              cp_dbcsr_cholesky_invert
   USE cp_dbcsr_diag,                   ONLY: cp_dbcsr_power
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              cp_dbcsr_dist2d_to_dist,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_add_on_diag, dbcsr_complete_redistribute, dbcsr_copy, dbcsr_create, &
        dbcsr_desymmetrize, dbcsr_distribution_get, dbcsr_distribution_release, &
        dbcsr_distribution_type, dbcsr_dot, dbcsr_filter, dbcsr_frobenius_norm, dbcsr_get_info, &
        dbcsr_get_num_blocks, dbcsr_multiply, dbcsr_p_type, dbcsr_release, dbcsr_scale, &
        dbcsr_type, dbcsr_type_antisymmetric, dbcsr_type_no_symmetry, dbcsr_type_symmetric
   USE dbt_api,                         ONLY: &
        dbt_batched_contract_finalize, dbt_batched_contract_init, dbt_clear, dbt_contract, &
        dbt_copy, dbt_copy_matrix_to_tensor, dbt_copy_tensor_to_matrix, dbt_create, dbt_destroy, &
        dbt_filter, dbt_get_block, dbt_get_info, dbt_get_num_blocks_total, &
        dbt_iterator_blocks_left, dbt_iterator_next_block, dbt_iterator_start, dbt_iterator_stop, &
        dbt_iterator_type, dbt_mp_environ_pgrid, dbt_nd_mp_comm, dbt_pgrid_create, &
        dbt_pgrid_destroy, dbt_pgrid_type, dbt_type
   USE distribution_2d_types,           ONLY: distribution_2d_type
   USE hfx_types,                       ONLY: hfx_ri_init,&
                                              hfx_ri_release,&
                                              hfx_ri_type
   USE input_constants,                 ONLY: hfx_ri_do_2c_cholesky,&
                                              hfx_ri_do_2c_diag,&
                                              hfx_ri_do_2c_iter
   USE input_cp2k_hfx,                  ONLY: ri_mo,&
                                              ri_pmat
   USE iterate_matrix,                  ONLY: invert_hotelling,&
                                              matrix_sqrt_newton_schulz
   USE kinds,                           ONLY: default_string_length,&
                                              dp,&
                                              int_8
   USE machine,                         ONLY: m_walltime
   USE message_passing,                 ONLY: mp_cart_create,&
                                              mp_environ,&
                                              mp_sum,&
                                              mp_sync
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_integral_utils,               ONLY: basis_set_list_setup
   USE qs_interactions,                 ONLY: init_interaction_radii_orb_basis
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_loc_methods,                  ONLY: qs_loc_driver
   USE qs_loc_types,                    ONLY: get_qs_loc_env,&
                                              qs_loc_env_create,&
                                              qs_loc_env_release
   USE qs_loc_utils,                    ONLY: qs_loc_control_init,&
                                              qs_loc_init
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type,&
                                              mo_set_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type,&
                                              release_neighbor_list_sets
   USE qs_tensors,                      ONLY: &
        build_2c_derivatives, build_2c_integrals, build_2c_neighbor_lists, build_3c_derivatives, &
        build_3c_integrals, build_3c_neighbor_lists, compress_tensor, decompress_tensor, &
        get_tensor_occupancy, neighbor_list_3c_destroy
   USE qs_tensors_types,                ONLY: create_2c_tensor,&
                                              create_3c_tensor,&
                                              create_tensor_batches,&
                                              distribution_3d_create,&
                                              distribution_3d_type,&
                                              neighbor_list_3c_type,&
                                              split_block_sizes
   USE util,                            ONLY: sort
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: hfx_ri_update_ks, hfx_ri_update_forces, get_force_from_3c_trace, get_2c_der_force, get_idx_to_atom

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_ri'
CONTAINS

! **************************************************************************************************
!> \brief Switches the RI_FLAVOR from MO to RHO or vice-versa
!> \param ri_data ...
!> \param qs_env ...
!> \note As a side product, the ri_data is mostly reinitialized and the integrals recomputed
! **************************************************************************************************
   SUBROUTINE switch_ri_flavor(ri_data, qs_env)
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'switch_ri_flavor'

      INTEGER                                            :: handle, n_mem, new_flavor
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (qs_kind_set, particle_set, atomic_kind_set, para_env, dft_control)

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, para_env=para_env, dft_control=dft_control, atomic_kind_set=atomic_kind_set, &
                      particle_set=particle_set, qs_kind_set=qs_kind_set)

      CALL hfx_ri_release(ri_data, write_stats=.FALSE.)

      IF (ri_data%flavor == ri_pmat) THEN
         new_flavor = ri_mo
      ELSE
         new_flavor = ri_pmat
      END IF
      ri_data%flavor = new_flavor

      n_mem = ri_data%n_mem_input
      ri_data%n_mem_input = ri_data%n_mem_flavor_switch
      ri_data%n_mem_flavor_switch = n_mem

      CALL hfx_ri_init(ri_data, qs_kind_set, particle_set, atomic_kind_set, para_env)

      !Need to recalculate the integrals in the new flavor
      !TODO: should we backup the integrals and symmetrize/desymmetrize them instead of recomputing ?!?
      !      that only makes sense if actual integral calculation is not negligible
      IF (ri_data%flavor == ri_pmat) THEN
         CALL hfx_ri_pre_scf_Pmat(qs_env, ri_data)
      ELSE
         CALL hfx_ri_pre_scf_mo(qs_env, ri_data, dft_control%nspins)
      END IF

      IF (ri_data%unit_nr > 0) THEN
         IF (ri_data%flavor == ri_pmat) THEN
            WRITE (ri_data%unit_nr, '(T2,A)') "HFX_RI_INFO| temporarily switched to RI_FLAVOR RHO"
         ELSE
            WRITE (ri_data%unit_nr, '(T2,A)') "HFX_RI_INFO| temporarily switched to RI_FLAVOR MO"
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE switch_ri_flavor

! **************************************************************************************************
!> \brief Pre-SCF steps in MO flavor of RI HFX
!>
!> Calculate 2-center & 3-center integrals (see hfx_ri_pre_scf_calc_tensors) and contract
!> K(P, S) = sum_R K_2(P, R)^{-1} K_1(R, S)^{1/2}
!> B(mu, lambda, R) = sum_P int_3c(mu, lambda, P) K(P, R)
!> \param qs_env ...
!> \param ri_data ...
!> \param nspins ...
! **************************************************************************************************
   SUBROUTINE hfx_ri_pre_scf_mo(qs_env, ri_data, nspins)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      INTEGER, INTENT(IN)                                :: nspins

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'hfx_ri_pre_scf_mo'

      INTEGER                                            :: handle, handle2, ispin, n_dependent, &
                                                            unit_nr, unit_nr_dbcsr
      REAL(KIND=dp)                                      :: threshold
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_type), DIMENSION(1) :: dbcsr_work_1, dbcsr_work_2, t_2c_int_mat, t_2c_op_pot, &
         t_2c_op_pot_sqrt, t_2c_op_pot_sqrt_inv, t_2c_op_RI, t_2c_op_RI_inv
      TYPE(dbt_type), DIMENSION(1)                       :: t_2c_int, t_2c_work
      TYPE(dbt_type), DIMENSION(1, 1)                    :: t_3c_int

      CALL timeset(routineN, handle)

      unit_nr_dbcsr = ri_data%unit_nr_dbcsr
      unit_nr = ri_data%unit_nr

      CALL get_qs_env(qs_env, para_env=para_env, blacs_env=blacs_env)

      CALL timeset(routineN//"_int", handle2)

      CALL hfx_ri_pre_scf_calc_tensors(qs_env, ri_data, t_2c_op_RI, t_2c_op_pot, t_3c_int)

      CALL timestop(handle2)

      CALL timeset(routineN//"_2c", handle2)
      IF (.NOT. ri_data%same_op) THEN
         SELECT CASE (ri_data%t2c_method)
         CASE (hfx_ri_do_2c_iter)
            CALL dbcsr_create(t_2c_op_RI_inv(1), template=t_2c_op_RI(1), matrix_type=dbcsr_type_no_symmetry)
            threshold = MAX(ri_data%filter_eps, 1.0e-12_dp)
            CALL invert_hotelling(t_2c_op_RI_inv(1), t_2c_op_RI(1), threshold=threshold, silent=.FALSE.)
         CASE (hfx_ri_do_2c_cholesky)
            CALL dbcsr_copy(t_2c_op_RI_inv(1), t_2c_op_RI(1))
            CALL cp_dbcsr_cholesky_decompose(t_2c_op_RI_inv(1), para_env=para_env, blacs_env=blacs_env)
            CALL cp_dbcsr_cholesky_invert(t_2c_op_RI_inv(1), para_env=para_env, blacs_env=blacs_env, upper_to_full=.TRUE.)
         CASE (hfx_ri_do_2c_diag)
            CALL dbcsr_copy(t_2c_op_RI_inv(1), t_2c_op_RI(1))
            CALL cp_dbcsr_power(t_2c_op_RI_inv(1), -1.0_dp, ri_data%eps_eigval, n_dependent, &
                                para_env, blacs_env, verbose=ri_data%unit_nr_dbcsr > 0)
         END SELECT

         IF (ri_data%check_2c_inv) THEN
            CALL check_inverse(t_2c_op_RI_inv(1), t_2c_op_RI(1), unit_nr=unit_nr)
         END IF

         CALL dbcsr_release(t_2c_op_RI(1))

         SELECT CASE (ri_data%t2c_method)
         CASE (hfx_ri_do_2c_iter)
            CALL dbcsr_create(t_2c_op_pot_sqrt(1), template=t_2c_op_pot(1), matrix_type=dbcsr_type_symmetric)
            CALL dbcsr_create(t_2c_op_pot_sqrt_inv(1), template=t_2c_op_pot(1), matrix_type=dbcsr_type_symmetric)
            CALL matrix_sqrt_newton_schulz(t_2c_op_pot_sqrt(1), t_2c_op_pot_sqrt_inv(1), t_2c_op_pot(1), &
                                           ri_data%filter_eps, ri_data%t2c_sqrt_order, ri_data%eps_lanczos, &
                                           ri_data%max_iter_lanczos)

            CALL dbcsr_release(t_2c_op_pot_sqrt_inv(1))
         CASE (hfx_ri_do_2c_diag, hfx_ri_do_2c_cholesky)
            CALL dbcsr_copy(t_2c_op_pot_sqrt(1), t_2c_op_pot(1))
            CALL cp_dbcsr_power(t_2c_op_pot_sqrt(1), 0.5_dp, ri_data%eps_eigval, n_dependent, &
                                para_env, blacs_env, verbose=ri_data%unit_nr_dbcsr > 0)
         END SELECT

         !We need S^-1 and (P|Q) for the forces.
         CALL dbt_create(t_2c_op_RI_inv(1), t_2c_work(1))
         CALL dbt_copy_matrix_to_tensor(t_2c_op_RI_inv(1), t_2c_work(1))
         CALL dbt_copy(t_2c_work(1), ri_data%t_2c_inv(1, 1), move_data=.TRUE.)
         CALL dbt_destroy(t_2c_work(1))
         CALL dbt_filter(ri_data%t_2c_inv(1, 1), ri_data%filter_eps)

         CALL dbt_create(t_2c_op_pot(1), t_2c_work(1))
         CALL dbt_copy_matrix_to_tensor(t_2c_op_pot(1), t_2c_work(1))
         CALL dbt_copy(t_2c_work(1), ri_data%t_2c_pot(1, 1), move_data=.TRUE.)
         CALL dbt_destroy(t_2c_work(1))
         CALL dbt_filter(ri_data%t_2c_pot(1, 1), ri_data%filter_eps)

         IF (ri_data%check_2c_inv) THEN
            CALL check_sqrt(t_2c_op_pot(1), matrix_sqrt=t_2c_op_pot_sqrt(1), unit_nr=unit_nr)
         END IF
         CALL dbcsr_create(t_2c_int_mat(1), template=t_2c_op_pot(1), matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_multiply("N", "N", 1.0_dp, t_2c_op_RI_inv(1), t_2c_op_pot_sqrt(1), &
                             0.0_dp, t_2c_int_mat(1), filter_eps=ri_data%filter_eps)
         CALL dbcsr_release(t_2c_op_RI_inv(1))
         CALL dbcsr_release(t_2c_op_pot_sqrt(1))
      ELSE
         SELECT CASE (ri_data%t2c_method)
         CASE (hfx_ri_do_2c_iter)
            CALL dbcsr_create(t_2c_int_mat(1), template=t_2c_op_pot(1), matrix_type=dbcsr_type_symmetric)
            CALL dbcsr_create(t_2c_op_pot_sqrt(1), template=t_2c_op_pot(1), matrix_type=dbcsr_type_symmetric)
            CALL matrix_sqrt_newton_schulz(t_2c_op_pot_sqrt(1), t_2c_int_mat(1), t_2c_op_pot(1), &
                                           ri_data%filter_eps, ri_data%t2c_sqrt_order, ri_data%eps_lanczos, &
                                           ri_data%max_iter_lanczos)
            CALL dbcsr_release(t_2c_op_pot_sqrt(1))
         CASE (hfx_ri_do_2c_diag, hfx_ri_do_2c_cholesky)
            CALL dbcsr_copy(t_2c_int_mat(1), t_2c_op_pot(1))
            CALL cp_dbcsr_power(t_2c_int_mat(1), -0.5_dp, ri_data%eps_eigval, n_dependent, &
                                para_env, blacs_env, verbose=ri_data%unit_nr_dbcsr > 0)
         END SELECT
         IF (ri_data%check_2c_inv) THEN
            CALL check_sqrt(t_2c_op_pot(1), matrix_sqrt_inv=t_2c_int_mat(1), unit_nr=unit_nr)
         END IF

         !We need (P|Q)^-1 for the forces
         CALL dbcsr_copy(dbcsr_work_1(1), t_2c_int_mat(1))
         CALL dbcsr_create(dbcsr_work_2(1), template=t_2c_int_mat(1))
         CALL dbcsr_multiply("N", "N", 1.0_dp, dbcsr_work_1(1), t_2c_int_mat(1), 0.0_dp, dbcsr_work_2(1))
         CALL dbcsr_release(dbcsr_work_1(1))
         CALL dbt_create(dbcsr_work_2(1), t_2c_work(1))
         CALL dbt_copy_matrix_to_tensor(dbcsr_work_2(1), t_2c_work(1))
         CALL dbcsr_release(dbcsr_work_2(1))
         CALL dbt_copy(t_2c_work(1), ri_data%t_2c_inv(1, 1), move_data=.TRUE.)
         CALL dbt_destroy(t_2c_work(1))
         CALL dbt_filter(ri_data%t_2c_inv(1, 1), ri_data%filter_eps)
      END IF

      CALL dbcsr_release(t_2c_op_pot(1))

      CALL dbt_create(t_2c_int_mat(1), t_2c_int(1), name="(RI|RI)")
      CALL dbt_copy_matrix_to_tensor(t_2c_int_mat(1), t_2c_int(1))
      CALL dbcsr_release(t_2c_int_mat(1))
      DO ispin = 1, nspins
         CALL dbt_copy(t_2c_int(1), ri_data%t_2c_int(ispin, 1))
      END DO
      CALL dbt_destroy(t_2c_int(1))
      CALL timestop(handle2)

      CALL timeset(routineN//"_3c", handle2)
      CALL dbt_copy(t_3c_int(1, 1), ri_data%t_3c_int_ctr_1(1, 1), order=[2, 1, 3], move_data=.TRUE.)
      CALL dbt_filter(ri_data%t_3c_int_ctr_1(1, 1), ri_data%filter_eps)
      CALL dbt_copy(ri_data%t_3c_int_ctr_1(1, 1), ri_data%t_3c_int_ctr_2(1, 1))
      CALL dbt_destroy(t_3c_int(1, 1))
      CALL timestop(handle2)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param matrix_1 ...
!> \param matrix_2 ...
!> \param name ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE check_inverse(matrix_1, matrix_2, name, unit_nr)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_1, matrix_2
      CHARACTER(len=*), INTENT(IN), OPTIONAL             :: name
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(len=default_string_length)               :: name_prv
      REAL(KIND=dp)                                      :: error, frob_matrix, frob_matrix_base
      TYPE(dbcsr_type)                                   :: matrix_tmp

      IF (PRESENT(name)) THEN
         name_prv = name
      ELSE
         CALL dbcsr_get_info(matrix_1, name=name_prv)
      END IF

      CALL dbcsr_create(matrix_tmp, template=matrix_1)
      CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_1, matrix_2, &
                          0.0_dp, matrix_tmp)
      frob_matrix_base = dbcsr_frobenius_norm(matrix_tmp)
      CALL dbcsr_add_on_diag(matrix_tmp, -1.0_dp)
      frob_matrix = dbcsr_frobenius_norm(matrix_tmp)
      error = frob_matrix/frob_matrix_base
      IF (unit_nr > 0) THEN
         WRITE (UNIT=unit_nr, FMT="(T3,A,A,A,T73,ES8.1)") &
            "HFX_RI_INFO| Error for INV(", TRIM(name_prv), "):", error
      END IF

      CALL dbcsr_release(matrix_tmp)
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param matrix_sqrt ...
!> \param matrix_sqrt_inv ...
!> \param name ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE check_sqrt(matrix, matrix_sqrt, matrix_sqrt_inv, name, unit_nr)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      TYPE(dbcsr_type), INTENT(IN), OPTIONAL             :: matrix_sqrt, matrix_sqrt_inv
      CHARACTER(len=*), INTENT(IN), OPTIONAL             :: name
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(len=default_string_length)               :: name_prv
      REAL(KIND=dp)                                      :: frob_matrix
      TYPE(dbcsr_type)                                   :: matrix_copy, matrix_tmp

      IF (PRESENT(name)) THEN
         name_prv = name
      ELSE
         CALL dbcsr_get_info(matrix, name=name_prv)
      END IF
      IF (PRESENT(matrix_sqrt)) THEN
         CALL dbcsr_create(matrix_tmp, template=matrix)
         CALL dbcsr_copy(matrix_copy, matrix_sqrt)
         CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_sqrt, matrix_copy, &
                             0.0_dp, matrix_tmp)
         CALL dbcsr_add(matrix_tmp, matrix, 1.0_dp, -1.0_dp)
         frob_matrix = dbcsr_frobenius_norm(matrix_tmp)
         IF (unit_nr > 0) THEN
            WRITE (UNIT=unit_nr, FMT="(T3,A,A,A,T73,ES8.1)") &
               "HFX_RI_INFO| Error for SQRT(", TRIM(name_prv), "):", frob_matrix
         END IF
         CALL dbcsr_release(matrix_tmp)
         CALL dbcsr_release(matrix_copy)
      END IF

      IF (PRESENT(matrix_sqrt_inv)) THEN
         CALL dbcsr_create(matrix_tmp, template=matrix)
         CALL dbcsr_copy(matrix_copy, matrix_sqrt_inv)
         CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_sqrt_inv, matrix_copy, &
                             0.0_dp, matrix_tmp)
         CALL check_inverse(matrix_tmp, matrix, name="SQRT("//TRIM(name_prv)//")", unit_nr=unit_nr)
         CALL dbcsr_release(matrix_tmp)
         CALL dbcsr_release(matrix_copy)
      END IF

   END SUBROUTINE

! **************************************************************************************************
!> \brief Calculate 2-center and 3-center integrals
!>
!> 2c: K_1(P, R) = (P|v1|R) and K_2(P, R) = (P|v2|R)
!> 3c: int_3c(mu, lambda, P) = (mu lambda |v2| P)
!> v_1 is HF operator, v_2 is RI metric
!> \param qs_env ...
!> \param ri_data ...
!> \param t_2c_int_RI K_2(P, R)
!> \param t_2c_int_pot K_1(P, R)
!> \param t_3c_int int_3c(mu, lambda, P)
! **************************************************************************************************
   SUBROUTINE hfx_ri_pre_scf_calc_tensors(qs_env, ri_data, t_2c_int_RI, t_2c_int_pot, t_3c_int)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(dbcsr_type), DIMENSION(1), INTENT(OUT)        :: t_2c_int_RI, t_2c_int_pot
      TYPE(dbt_type), DIMENSION(1, 1), INTENT(OUT)       :: t_3c_int

      CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_ri_pre_scf_calc_tensors'

      INTEGER                                            :: handle, i_mem, ibasis, mp_comm_t3c, &
                                                            n_mem, natom, nkind
      INTEGER, ALLOCATABLE, DIMENSION(:) :: dist_AO_1, dist_AO_2, dist_RI, &
         ends_array_mc_block_int, ends_array_mc_int, sizes_AO, sizes_RI, &
         starts_array_mc_block_int, starts_array_mc_int
      INTEGER, DIMENSION(3)                              :: pcoord, pdims
      INTEGER, DIMENSION(:), POINTER                     :: col_bsize, row_bsize
      LOGICAL                                            :: converged
      REAL(dp)                                           :: max_ev, min_ev
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dbcsr_distribution_type)                      :: dbcsr_dist
      TYPE(dbt_type)                                     :: t_3c_tmp
      TYPE(dbt_type), DIMENSION(1, 1)                    :: t_3c_int_batched
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(distribution_2d_type), POINTER                :: dist_2d
      TYPE(distribution_3d_type)                         :: dist_3d
      TYPE(gto_basis_set_p_type), ALLOCATABLE, &
         DIMENSION(:), TARGET                            :: basis_set_AO, basis_set_RI
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis, ri_basis
      TYPE(neighbor_list_3c_type)                        :: nl_3c
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: nl_2c_pot, nl_2c_RI
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env

      CALL timeset(routineN, handle)
      NULLIFY (col_bsize, row_bsize, dist_2d, nl_2c_pot, nl_2c_RI, &
               particle_set, qs_kind_set, ks_env)

      CALL get_qs_env(qs_env, natom=natom, nkind=nkind, qs_kind_set=qs_kind_set, particle_set=particle_set, &
                      distribution_2d=dist_2d, ks_env=ks_env, dft_control=dft_control)

      ALLOCATE (sizes_RI(natom), sizes_AO(natom))
      ALLOCATE (basis_set_RI(nkind), basis_set_AO(nkind))
      CALL basis_set_list_setup(basis_set_RI, ri_data%ri_basis_type, qs_kind_set)
      CALL get_particle_set(particle_set, qs_kind_set, nsgf=sizes_RI, basis=basis_set_RI)
      CALL basis_set_list_setup(basis_set_AO, ri_data%orb_basis_type, qs_kind_set)

      CALL get_particle_set(particle_set, qs_kind_set, nsgf=sizes_AO, basis=basis_set_AO)

      DO ibasis = 1, SIZE(basis_set_AO)
         orb_basis => basis_set_AO(ibasis)%gto_basis_set
         ri_basis => basis_set_RI(ibasis)%gto_basis_set
         ! interaction radii should be based on eps_pgf_orb controlled in RI section
         ! (since hartree-fock needs very tight eps_pgf_orb for Kohn-Sham/Fock matrix but eps_pgf_orb
         ! can be much looser in RI HFX since no systematic error is introduced with tensor sparsity)
         CALL init_interaction_radii_orb_basis(orb_basis, ri_data%eps_pgf_orb)
         CALL init_interaction_radii_orb_basis(ri_basis, ri_data%eps_pgf_orb)
      END DO

      n_mem = ri_data%n_mem
      CALL create_tensor_batches(sizes_AO, n_mem, starts_array_mc_int, ends_array_mc_int, &
                                 starts_array_mc_block_int, ends_array_mc_block_int)

      DEALLOCATE (starts_array_mc_int, ends_array_mc_int)

      CALL create_3c_tensor(t_3c_int_batched(1, 1), dist_RI, dist_AO_1, dist_AO_2, ri_data%pgrid, &
                            sizes_RI, sizes_AO, sizes_AO, map1=[1], map2=[2, 3], &
                            name="(RI | AO AO)")

      CALL get_qs_env(qs_env, nkind=nkind, particle_set=particle_set, atomic_kind_set=atomic_kind_set)
      CALL dbt_mp_environ_pgrid(ri_data%pgrid, pdims, pcoord)
      CALL mp_cart_create(ri_data%pgrid%mp_comm_2d, 3, pdims, pcoord, mp_comm_t3c)
      CALL distribution_3d_create(dist_3d, dist_RI, dist_AO_1, dist_AO_2, &
                                  nkind, particle_set, mp_comm_t3c, own_comm=.TRUE.)
      DEALLOCATE (dist_RI, dist_AO_1, dist_AO_2)

      CALL create_3c_tensor(t_3c_int(1, 1), dist_RI, dist_AO_1, dist_AO_2, ri_data%pgrid, &
                            ri_data%bsizes_RI_split, ri_data%bsizes_AO_split, ri_data%bsizes_AO_split, &
                            map1=[1], map2=[2, 3], &
                            name="O (RI AO | AO)")

      ! create 3c tensor for storage of ints

      CALL build_3c_neighbor_lists(nl_3c, basis_set_RI, basis_set_AO, basis_set_AO, dist_3d, ri_data%ri_metric, &
                                   "HFX_3c_nl", qs_env, op_pos=1, sym_jk=.TRUE., own_dist=.TRUE.)

      DO i_mem = 1, n_mem
         CALL build_3c_integrals(t_3c_int_batched, ri_data%filter_eps/2, qs_env, nl_3c, &
                                 basis_set_RI, basis_set_AO, basis_set_AO, &
                                 ri_data%ri_metric, int_eps=ri_data%eps_schwarz, op_pos=1, &
                                 desymmetrize=.FALSE., &
                                 bounds_j=[starts_array_mc_block_int(i_mem), ends_array_mc_block_int(i_mem)])
         CALL dbt_copy(t_3c_int_batched(1, 1), t_3c_int(1, 1), summation=.TRUE., move_data=.TRUE.)
         CALL dbt_filter(t_3c_int(1, 1), ri_data%filter_eps/2)
      END DO

      CALL dbt_destroy(t_3c_int_batched(1, 1))

      CALL neighbor_list_3c_destroy(nl_3c)

      CALL dbt_create(t_3c_int(1, 1), t_3c_tmp)

      IF (ri_data%flavor == ri_pmat) THEN ! desymmetrize
         ! desymmetrize
         CALL dbt_copy(t_3c_int(1, 1), t_3c_tmp)
         CALL dbt_copy(t_3c_tmp, t_3c_int(1, 1), order=[1, 3, 2], summation=.TRUE., move_data=.TRUE.)

         ! For RI-RHO filter_eps_storage is reserved for screening tensor contracted with RI-metric
         ! with RI metric but not to bare integral tensor
         CALL dbt_filter(t_3c_int(1, 1), ri_data%filter_eps)
      ELSE
         CALL dbt_filter(t_3c_int(1, 1), ri_data%filter_eps_storage/2)
      END IF

      CALL dbt_destroy(t_3c_tmp)

      CALL build_2c_neighbor_lists(nl_2c_pot, basis_set_RI, basis_set_RI, ri_data%hfx_pot, &
                                   "HFX_2c_nl_pot", &
                                   qs_env, sym_ij=.TRUE., &
                                   dist_2d=dist_2d)

      CALL cp_dbcsr_dist2d_to_dist(dist_2d, dbcsr_dist)
      ALLOCATE (row_bsize(SIZE(sizes_RI)))
      ALLOCATE (col_bsize(SIZE(sizes_RI)))
      row_bsize(:) = sizes_RI
      col_bsize(:) = sizes_RI

      CALL dbcsr_create(t_2c_int_pot(1), "(R|P) HFX", dbcsr_dist, dbcsr_type_symmetric, &
                        row_bsize, col_bsize, reuse_arrays=.TRUE.)

      CALL dbcsr_distribution_release(dbcsr_dist)

      CALL build_2c_integrals(t_2c_int_pot, ri_data%filter_eps_2c, qs_env, nl_2c_pot, basis_set_RI, basis_set_RI, &
                              ri_data%hfx_pot)
      CALL release_neighbor_list_sets(nl_2c_pot)

      IF (.NOT. ri_data%same_op) THEN
         CALL build_2c_neighbor_lists(nl_2c_RI, basis_set_RI, basis_set_RI, ri_data%ri_metric, &
                                      "HFX_2c_nl_RI", &
                                      qs_env, sym_ij=.TRUE., &
                                      dist_2d=dist_2d)

         CALL dbcsr_create(t_2c_int_RI(1), template=t_2c_int_pot(1), matrix_type=dbcsr_type_symmetric, name="(R|P) RI")
         CALL build_2c_integrals(t_2c_int_RI, ri_data%filter_eps_2c, qs_env, nl_2c_RI, basis_set_RI, basis_set_RI, &
                                 ri_data%ri_metric)

         CALL release_neighbor_list_sets(nl_2c_RI)
      END IF

      DO ibasis = 1, SIZE(basis_set_AO)
         orb_basis => basis_set_AO(ibasis)%gto_basis_set
         ri_basis => basis_set_RI(ibasis)%gto_basis_set
         ! reset interaction radii of orb basis
         CALL init_interaction_radii_orb_basis(orb_basis, dft_control%qs_control%eps_pgf_orb)
         CALL init_interaction_radii_orb_basis(ri_basis, dft_control%qs_control%eps_pgf_orb)
      END DO

      IF (ri_data%calc_condnum) THEN
         CALL arnoldi_extremal(t_2c_int_pot(1), max_ev, min_ev, threshold=ri_data%eps_lanczos, &
                               max_iter=ri_data%max_iter_lanczos, converged=converged)

         IF (.NOT. converged) THEN
            CPWARN("Condition number estimate of (P|Q) (HFX potential) is not reliable (not converged).")
         END IF

         IF (ri_data%unit_nr > 0) THEN
            WRITE (ri_data%unit_nr, '(T2,A)') "2-Norm Condition Number of (P|Q) integrals (HFX potential)"
            IF (min_ev > 0) THEN
               WRITE (ri_data%unit_nr, '(T4,A,ES11.3E3,T32,A,ES11.3E3,A4,ES11.3E3,T63,A,F8.4)') &
                  "CN : max/min ev: ", max_ev, " / ", min_ev, "=", max_ev/min_ev, "Log(2-CN):", LOG10(max_ev/min_ev)
            ELSE
               WRITE (ri_data%unit_nr, '(T4,A,ES11.3E3,T32,A,ES11.3E3,T63,A)') &
                  "CN : max/min ev: ", max_ev, " / ", min_ev, "Log(CN): infinity"
            END IF
         END IF

         IF (.NOT. ri_data%same_op) THEN
            CALL arnoldi_extremal(t_2c_int_RI(1), max_ev, min_ev, threshold=ri_data%eps_lanczos, &
                                  max_iter=ri_data%max_iter_lanczos, converged=converged)

            IF (.NOT. converged) THEN
               CPWARN("Condition number estimate of (P|Q) matrix (RI metric) is not reliable (not converged).")
            END IF

            IF (ri_data%unit_nr > 0) THEN
               WRITE (ri_data%unit_nr, '(T2,A)') "2-Norm Condition Number of (P|Q) integrals (RI metric)"
               IF (min_ev > 0) THEN
                  WRITE (ri_data%unit_nr, '(T4,A,ES11.3E3,T32,A,ES11.3E3,A4,ES11.3E3,T63,A,F8.4)') &
                     "CN : max/min ev: ", max_ev, " / ", min_ev, "=", max_ev/min_ev, "Log(2-CN):", LOG10(max_ev/min_ev)
               ELSE
                  WRITE (ri_data%unit_nr, '(T4,A,ES11.3E3,T32,A,ES11.3E3,T63,A)') &
                     "CN : max/min ev: ", max_ev, " / ", min_ev, "Log(CN): infinity"
               END IF
            END IF
         END IF
      END IF

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Pre-SCF steps in rho flavor of RI HFX
!>
!> K(P, S) = sum_{R,Q} K_2(P, R)^{-1} K_1(R, Q) K_2(Q, S)^{-1}
!> Calculate B(mu, lambda, R) = sum_P int_3c(mu, lambda, P) K(P, R)
!> \param qs_env ...
!> \param ri_data ...
! **************************************************************************************************
   SUBROUTINE hfx_ri_pre_scf_Pmat(qs_env, ri_data)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data

      CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_ri_pre_scf_Pmat'

      INTEGER                                            :: handle, handle2, i_mem, j_mem, &
                                                            n_dependent, unit_nr, unit_nr_dbcsr
      INTEGER(int_8)                                     :: nflop, nze, nze_O
      INTEGER, DIMENSION(2, 1)                           :: bounds_i
      INTEGER, DIMENSION(2, 2)                           :: bounds_j
      INTEGER, DIMENSION(3)                              :: dims_3c
      REAL(KIND=dp)                                      :: compression_factor, memory_3c, occ, &
                                                            threshold
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_type), DIMENSION(1)                     :: t_2c_int_mat, t_2c_op_pot, t_2c_op_RI, &
                                                            t_2c_tmp, t_2c_tmp_2
      TYPE(dbt_type)                                     :: t_3c_2
      TYPE(dbt_type), DIMENSION(1)                       :: t_2c_int, t_2c_work
      TYPE(dbt_type), DIMENSION(1, 1)                    :: t_3c_int_1

      CALL timeset(routineN, handle)

      unit_nr_dbcsr = ri_data%unit_nr_dbcsr
      unit_nr = ri_data%unit_nr

      CALL get_qs_env(qs_env, para_env=para_env, blacs_env=blacs_env)

      CALL timeset(routineN//"_int", handle2)

      CALL hfx_ri_pre_scf_calc_tensors(qs_env, ri_data, t_2c_op_RI, t_2c_op_pot, t_3c_int_1)

      CALL dbt_copy(t_3c_int_1(1, 1), ri_data%t_3c_int_ctr_3(1, 1), order=[1, 2, 3], move_data=.TRUE.)

      CALL dbt_destroy(t_3c_int_1(1, 1))

      CALL timestop(handle2)

      CALL timeset(routineN//"_2c", handle2)

      IF (ri_data%same_op) t_2c_op_RI(1) = t_2c_op_pot(1)
      CALL dbcsr_create(t_2c_int_mat(1), template=t_2c_op_RI(1), matrix_type=dbcsr_type_no_symmetry)
      threshold = MAX(ri_data%filter_eps, 1.0e-12_dp)

      SELECT CASE (ri_data%t2c_method)
      CASE (hfx_ri_do_2c_iter)
         CALL invert_hotelling(t_2c_int_mat(1), t_2c_op_RI(1), &
                               threshold=threshold, silent=.FALSE.)
      CASE (hfx_ri_do_2c_cholesky)
         CALL dbcsr_copy(t_2c_int_mat(1), t_2c_op_RI(1))
         CALL cp_dbcsr_cholesky_decompose(t_2c_int_mat(1), para_env=para_env, blacs_env=blacs_env)
         CALL cp_dbcsr_cholesky_invert(t_2c_int_mat(1), para_env=para_env, blacs_env=blacs_env, upper_to_full=.TRUE.)
      CASE (hfx_ri_do_2c_diag)
         CALL dbcsr_copy(t_2c_int_mat(1), t_2c_op_RI(1))
         CALL cp_dbcsr_power(t_2c_int_mat(1), -1.0_dp, ri_data%eps_eigval, n_dependent, &
                             para_env, blacs_env, verbose=ri_data%unit_nr_dbcsr > 0)
      END SELECT

      IF (ri_data%check_2c_inv) THEN
         CALL check_inverse(t_2c_int_mat(1), t_2c_op_RI(1), unit_nr=unit_nr)
      END IF

      !Need to save the (P|Q)^-1 tensor for forces (inverse metric if not same_op)
      CALL dbt_create(t_2c_int_mat(1), t_2c_work(1))
      CALL dbt_copy_matrix_to_tensor(t_2c_int_mat(1), t_2c_work(1))
      CALL dbt_copy(t_2c_work(1), ri_data%t_2c_inv(1, 1), move_data=.TRUE.)
      CALL dbt_destroy(t_2c_work(1))
      CALL dbt_filter(ri_data%t_2c_inv(1, 1), ri_data%filter_eps)
      IF (.NOT. ri_data%same_op) THEN
         !Also save the RI (P|Q) integral
         CALL dbt_create(t_2c_op_pot(1), t_2c_work(1))
         CALL dbt_copy_matrix_to_tensor(t_2c_op_pot(1), t_2c_work(1))
         CALL dbt_copy(t_2c_work(1), ri_data%t_2c_pot(1, 1), move_data=.TRUE.)
         CALL dbt_destroy(t_2c_work(1))
         CALL dbt_filter(ri_data%t_2c_pot(1, 1), ri_data%filter_eps)
      END IF

      IF (ri_data%same_op) THEN
         CALL dbcsr_release(t_2c_op_pot(1))
      ELSE
         CALL dbcsr_create(t_2c_tmp(1), template=t_2c_op_RI(1), matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(t_2c_tmp_2(1), template=t_2c_op_RI(1), matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_release(t_2c_op_RI(1))
         CALL dbcsr_multiply('N', 'N', 1.0_dp, t_2c_int_mat(1), t_2c_op_pot(1), 0.0_dp, t_2c_tmp(1), &
                             filter_eps=ri_data%filter_eps)

         CALL dbcsr_release(t_2c_op_pot(1))
         CALL dbcsr_multiply('N', 'N', 1.0_dp, t_2c_tmp(1), t_2c_int_mat(1), 0.0_dp, t_2c_tmp_2(1), &
                             filter_eps=ri_data%filter_eps)
         CALL dbcsr_release(t_2c_tmp(1))
         CALL dbcsr_release(t_2c_int_mat(1))
         t_2c_int_mat(1) = t_2c_tmp_2(1)
      END IF

      CALL dbt_create(t_2c_int_mat(1), t_2c_int(1), name="(RI|RI)")
      CALL dbt_copy_matrix_to_tensor(t_2c_int_mat(1), t_2c_int(1))
      CALL dbcsr_release(t_2c_int_mat(1))
      CALL dbt_copy(t_2c_int(1), ri_data%t_2c_int(1, 1), move_data=.TRUE.)
      CALL dbt_destroy(t_2c_int(1))
      CALL dbt_filter(ri_data%t_2c_int(1, 1), ri_data%filter_eps)

      CALL timestop(handle2)

      CALL dbt_create(ri_data%t_3c_int_ctr_3(1, 1), t_3c_2)

      CALL dbt_get_info(ri_data%t_3c_int_ctr_3(1, 1), nfull_total=dims_3c)

      memory_3c = 0.0_dp
      nze_O = 0

      DO i_mem = 1, ri_data%n_mem_RI
         bounds_i(:, 1) = [ri_data%starts_array_RI_mem(i_mem), ri_data%ends_array_RI_mem(i_mem)]
         CALL dbt_batched_contract_init(ri_data%t_2c_int(1, 1))
         CALL dbt_batched_contract_init(ri_data%t_3c_int_ctr_3(1, 1))
         CALL dbt_batched_contract_init(t_3c_2)
         DO j_mem = 1, ri_data%n_mem
            bounds_j(:, 1) = [ri_data%starts_array_mem(j_mem), ri_data%ends_array_mem(j_mem)]
            bounds_j(:, 2) = [1, dims_3c(3)]
            CALL timeset(routineN//"_RIx3C", handle2)
            CALL dbt_contract(1.0_dp, ri_data%t_2c_int(1, 1), ri_data%t_3c_int_ctr_3(1, 1), &
                              0.0_dp, t_3c_2, &
                              contract_1=[2], notcontract_1=[1], &
                              contract_2=[1], notcontract_2=[2, 3], &
                              map_1=[1], map_2=[2, 3], filter_eps=ri_data%filter_eps_storage, &
                              bounds_2=bounds_i, &
                              bounds_3=bounds_j, &
                              unit_nr=unit_nr_dbcsr, &
                              flop=nflop)

            ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop
            CALL timestop(handle2)

            CALL timeset(routineN//"_copy_2", handle2)
            CALL dbt_copy(t_3c_2, ri_data%t_3c_int_ctr_1(1, 1), order=[2, 1, 3], move_data=.TRUE.)

            CALL get_tensor_occupancy(ri_data%t_3c_int_ctr_1(1, 1), nze, occ)
            nze_O = nze_O + nze

            CALL compress_tensor(ri_data%t_3c_int_ctr_1(1, 1), ri_data%blk_indices(j_mem, i_mem)%ind, &
                                 ri_data%store_3c(j_mem, i_mem), ri_data%filter_eps_storage, memory_3c)

            CALL timestop(handle2)
         END DO
         CALL dbt_batched_contract_finalize(ri_data%t_2c_int(1, 1))
         CALL dbt_batched_contract_finalize(ri_data%t_3c_int_ctr_3(1, 1))
         CALL dbt_batched_contract_finalize(t_3c_2)
      END DO

      CALL mp_sum(memory_3c, para_env%group)
      compression_factor = REAL(nze_O, dp)*1.0E-06*8.0_dp/memory_3c

      IF (unit_nr > 0) THEN
         WRITE (UNIT=unit_nr, FMT="((T3,A,T66,F11.2,A4))") &
            "MEMORY_INFO| Memory for 3-center integrals (compressed):", memory_3c, ' MiB'

         WRITE (UNIT=unit_nr, FMT="((T3,A,T60,F21.2))") &
            "MEMORY_INFO| Compression factor:                  ", compression_factor
      END IF

      CALL dbt_clear(ri_data%t_2c_int(1, 1))
      CALL dbt_destroy(t_3c_2)

      CALL dbt_copy(ri_data%t_3c_int_ctr_3(1, 1), ri_data%t_3c_int_ctr_2(1, 1), order=[2, 1, 3], move_data=.TRUE.)

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Sorts 2d indices w.r.t. rows and columns
!> \param blk_ind ...
! **************************************************************************************************
   SUBROUTINE sort_unique_blkind_2d(blk_ind)
      INTEGER, ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: blk_ind

      INTEGER                                            :: end_ind, iblk, iblk_all, irow, nblk, &
                                                            ncols, start_ind
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ind_1, ind_2, sort_1, sort_2
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: blk_ind_tmp

      nblk = SIZE(blk_ind, 1)

      ALLOCATE (sort_1(nblk))
      ALLOCATE (ind_1(nblk))

      sort_1(:) = blk_ind(:, 1)
      CALL sort(sort_1, nblk, ind_1)

      blk_ind(:, :) = blk_ind(ind_1, :)

      start_ind = 1

      DO WHILE (start_ind <= nblk)
         irow = blk_ind(start_ind, 1)
         end_ind = start_ind

         IF (end_ind + 1 <= nblk) THEN
         DO WHILE (blk_ind(end_ind + 1, 1) == irow)
            end_ind = end_ind + 1
            IF (end_ind + 1 > nblk) EXIT
         END DO
         END IF

         ncols = end_ind - start_ind + 1
         ALLOCATE (sort_2(ncols))
         ALLOCATE (ind_2(ncols))
         sort_2(:) = blk_ind(start_ind:end_ind, 2)
         CALL sort(sort_2, ncols, ind_2)
         ind_2 = ind_2 + start_ind - 1

         blk_ind(start_ind:end_ind, :) = blk_ind(ind_2, :)
         start_ind = end_ind + 1

         DEALLOCATE (sort_2, ind_2)
      END DO

      ALLOCATE (blk_ind_tmp(nblk, 2))
      blk_ind_tmp = 0

      iblk = 0
      DO iblk_all = 1, nblk
         IF (iblk >= 1) THEN
            IF (ALL(blk_ind_tmp(iblk, :) == blk_ind(iblk_all, :))) THEN
               CYCLE
            END IF
         END IF
         iblk = iblk + 1
         blk_ind_tmp(iblk, :) = blk_ind(iblk_all, :)
      END DO
      nblk = iblk

      DEALLOCATE (blk_ind)
      ALLOCATE (blk_ind(nblk, 2))

      blk_ind(:, :) = blk_ind_tmp(:nblk, :)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param ri_data ...
!> \param ks_matrix ...
!> \param ehfx ...
!> \param mos ...
!> \param rho_ao ...
!> \param geometry_did_change ...
!> \param nspins ...
!> \param hf_fraction ...
! **************************************************************************************************
   SUBROUTINE hfx_ri_update_ks(qs_env, ri_data, ks_matrix, ehfx, mos, rho_ao, &
                               geometry_did_change, nspins, hf_fraction)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ks_matrix
      REAL(KIND=dp), INTENT(OUT)                         :: ehfx
      TYPE(mo_set_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: mos
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao
      LOGICAL, INTENT(IN)                                :: geometry_did_change
      INTEGER, INTENT(IN)                                :: nspins
      REAL(KIND=dp), INTENT(IN)                          :: hf_fraction

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'hfx_ri_update_ks'

      CHARACTER(1)                                       :: mtype
      INTEGER                                            :: handle, handle2, i, ispin, j
      INTEGER(int_8)                                     :: nblks
      INTEGER, DIMENSION(2)                              :: homo
      LOGICAL                                            :: is_antisymmetric
      REAL(dp)                                           :: etmp, fac
      REAL(KIND=dp), DIMENSION(:), POINTER               :: mo_eigenvalues
      TYPE(cp_1d_r_p_type), DIMENSION(:), POINTER        :: occupied_evals
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: homo_localized, moloc_coeff, &
                                                            occupied_orbs
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: my_ks_matrix, my_rho_ao
      TYPE(dbcsr_type), DIMENSION(2)                     :: mo_coeff_b
      TYPE(dbcsr_type), POINTER                          :: mo_coeff_b_tmp
      TYPE(mo_set_type), POINTER                         :: mo_set

      CALL timeset(routineN, handle)

      IF (nspins == 1) THEN
         fac = 0.5_dp*hf_fraction
      ELSE
         fac = 1.0_dp*hf_fraction
      END IF

      !If incoming assymetric matrices, need to convert to normal
      NULLIFY (my_ks_matrix, my_rho_ao)
      CALL dbcsr_get_info(ks_matrix(1, 1)%matrix, matrix_type=mtype)
      is_antisymmetric = mtype == dbcsr_type_antisymmetric
      IF (is_antisymmetric) THEN
         ALLOCATE (my_ks_matrix(SIZE(ks_matrix, 1), SIZE(ks_matrix, 2)))
         ALLOCATE (my_rho_ao(SIZE(rho_ao, 1), SIZE(rho_ao, 2)))

         DO i = 1, SIZE(ks_matrix, 1)
            DO j = 1, SIZE(ks_matrix, 2)
               ALLOCATE (my_ks_matrix(i, j)%matrix, my_rho_ao(i, j)%matrix)
               CALL dbcsr_create(my_ks_matrix(i, j)%matrix, template=ks_matrix(i, j)%matrix, &
                                 matrix_type=dbcsr_type_no_symmetry)
               CALL dbcsr_desymmetrize(ks_matrix(i, j)%matrix, my_ks_matrix(i, j)%matrix)
               CALL dbcsr_create(my_rho_ao(i, j)%matrix, template=rho_ao(i, j)%matrix, &
                                 matrix_type=dbcsr_type_no_symmetry)
               CALL dbcsr_desymmetrize(rho_ao(i, j)%matrix, my_rho_ao(i, j)%matrix)
            END DO
         END DO
      ELSE
         my_ks_matrix => ks_matrix
         my_rho_ao => rho_ao
      END IF

      !Case analysis on RI_FLAVOR: we switch if the input flavor is MO, there is no provided MO, and
      !                            the current flavor is not yet RHO. We switch back to MO if there are
      !                            MOs available and the current flavor is actually RHO
      IF (ri_data%input_flavor == ri_mo .AND. (.NOT. PRESENT(mos)) .AND. ri_data%flavor == ri_mo) THEN
         CALL switch_ri_flavor(ri_data, qs_env)
      ELSE IF (ri_data%input_flavor == ri_mo .AND. PRESENT(mos) .AND. ri_data%flavor == ri_pmat) THEN
         CALL switch_ri_flavor(ri_data, qs_env)
      END IF

      SELECT CASE (ri_data%flavor)
      CASE (ri_mo)
         CPASSERT(PRESENT(mos))
         CALL timeset(routineN//"_MO", handle2)

         IF (ri_data%do_loc) THEN
            ALLOCATE (occupied_orbs(nspins))
            ALLOCATE (occupied_evals(nspins))
            ALLOCATE (homo_localized(nspins))
         END IF
         DO ispin = 1, nspins
            NULLIFY (mo_coeff_b_tmp)
            mo_set => mos(ispin)%mo_set
            CPASSERT(mo_set%uniform_occupation)
            CALL get_mo_set(mo_set=mo_set, mo_coeff=mo_coeff, eigenvalues=mo_eigenvalues, mo_coeff_b=mo_coeff_b_tmp)

            IF (.NOT. ri_data%do_loc) THEN
               IF (.NOT. mo_set%use_mo_coeff_b) CALL copy_fm_to_dbcsr(mo_coeff, mo_coeff_b_tmp)
               CALL dbcsr_copy(mo_coeff_b(ispin), mo_coeff_b_tmp)
            ELSE
               IF (mo_set%use_mo_coeff_b) CALL copy_dbcsr_to_fm(mo_coeff_b_tmp, mo_coeff)
               CALL dbcsr_create(mo_coeff_b(ispin), template=mo_coeff_b_tmp)
            END IF

            IF (ri_data%do_loc) THEN
               occupied_orbs(ispin)%matrix => mo_coeff
               occupied_evals(ispin)%array => mo_eigenvalues
               CALL cp_fm_create(homo_localized(ispin)%matrix, occupied_orbs(ispin)%matrix%matrix_struct)
               CALL cp_fm_to_fm(occupied_orbs(ispin)%matrix, homo_localized(ispin)%matrix)
            END IF
         END DO

         IF (ri_data%do_loc) THEN
            CALL qs_loc_env_create(ri_data%qs_loc_env)
            CALL qs_loc_control_init(ri_data%qs_loc_env, ri_data%loc_subsection, do_homo=.TRUE.)
            CALL qs_loc_init(qs_env, ri_data%qs_loc_env, ri_data%loc_subsection, homo_localized)
            DO ispin = 1, nspins
               CALL qs_loc_driver(qs_env, ri_data%qs_loc_env, ri_data%print_loc_subsection, ispin, &
                                  ext_mo_coeff=homo_localized(ispin)%matrix)
            END DO
            CALL get_qs_loc_env(qs_loc_env=ri_data%qs_loc_env, moloc_coeff=moloc_coeff)

            DO ispin = 1, nspins
               CALL cp_fm_release(homo_localized(ispin)%matrix)
            END DO

            DEALLOCATE (occupied_orbs, occupied_evals, homo_localized)
         END IF

         DO ispin = 1, nspins
            mo_set => mos(ispin)%mo_set
            IF (ri_data%do_loc) THEN
               CALL copy_fm_to_dbcsr(moloc_coeff(ispin)%matrix, mo_coeff_b(ispin))
            END IF
            CALL dbcsr_scale(mo_coeff_b(ispin), SQRT(mo_set%maxocc))
            homo(ispin) = mo_set%homo
         END DO

         IF (ri_data%do_loc) CALL qs_loc_env_release(ri_data%qs_loc_env)
         CALL timestop(handle2)

         CALL hfx_ri_update_ks_mo(qs_env, ri_data, my_ks_matrix, mo_coeff_b, homo, &
                                  geometry_did_change, nspins, fac)
      CASE (ri_pmat)

         NULLIFY (para_env)
         CALL get_qs_env(qs_env, para_env=para_env)
         DO ispin = 1, SIZE(my_rho_ao, 1)
            nblks = dbcsr_get_num_blocks(my_rho_ao(ispin, 1)%matrix)
            CALL mp_sum(nblks, para_env%group)
            IF (nblks == 0) THEN
               CPABORT("received empty density matrix")
            END IF
         END DO

         CALL hfx_ri_update_ks_pmat(qs_env, ri_data, my_ks_matrix, my_rho_ao, &
                                    geometry_did_change, nspins, fac)

      END SELECT

      DO ispin = 1, nspins
         CALL dbcsr_release(mo_coeff_b(ispin))
      END DO

      DO ispin = 1, nspins
         CALL dbcsr_filter(my_ks_matrix(ispin, 1)%matrix, ri_data%filter_eps)
      END DO

      CALL timeset(routineN//"_energy", handle2)
      ! Calculate the exchange energy
      ehfx = 0.0_dp
      DO ispin = 1, nspins
         CALL dbcsr_dot(my_ks_matrix(ispin, 1)%matrix, my_rho_ao(ispin, 1)%matrix, &
                        etmp)
         ehfx = ehfx + 0.5_dp*etmp

      END DO
      CALL timestop(handle2)

      !Anti-symmetric case
      IF (is_antisymmetric) THEN
         DO i = 1, SIZE(ks_matrix, 1)
            DO j = 1, SIZE(ks_matrix, 2)
               CALL dbcsr_complete_redistribute(my_ks_matrix(i, j)%matrix, ks_matrix(i, j)%matrix)
               CALL dbcsr_complete_redistribute(my_rho_ao(i, j)%matrix, rho_ao(i, j)%matrix)
            END DO
         END DO
         CALL dbcsr_deallocate_matrix_set(my_ks_matrix)
         CALL dbcsr_deallocate_matrix_set(my_rho_ao)
      END IF

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Calculate Fock (AKA Kohn-Sham) matrix in MO flavor
!>
!> C(mu, i) (MO coefficients)
!> M(mu, i, R) = sum_nu B(mu, nu, R) C(nu, i)
!> KS(mu, lambda) = sum_{i,R} M(mu, i, R) M(lambda, i, R)
!> \param qs_env ...
!> \param ri_data ...
!> \param ks_matrix ...
!> \param mo_coeff C(mu, i)
!> \param homo ...
!> \param geometry_did_change ...
!> \param nspins ...
!> \param fac ...
! **************************************************************************************************
   SUBROUTINE hfx_ri_update_ks_mo(qs_env, ri_data, ks_matrix, mo_coeff, &
                                  homo, geometry_did_change, nspins, fac)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(dbcsr_p_type), DIMENSION(:, :)                :: ks_matrix
      TYPE(dbcsr_type), DIMENSION(:), INTENT(IN)         :: mo_coeff
      INTEGER, DIMENSION(:)                              :: homo
      LOGICAL, INTENT(IN)                                :: geometry_did_change
      INTEGER, INTENT(IN)                                :: nspins
      REAL(dp), INTENT(IN)                               :: fac

      CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_ri_update_ks_mo'

      INTEGER                                            :: bsize, bsum, comm_2d, handle, handle2, &
                                                            i_mem, iblock, iproc, ispin, n_mem, &
                                                            n_mos, nblock, nproc, unit_nr_dbcsr
      INTEGER(int_8)                                     :: nblks, nflop
      INTEGER, ALLOCATABLE, DIMENSION(:) :: batch_ranges_1, batch_ranges_2, dist1, dist2, dist3, &
         mem_end, mem_end_block_1, mem_end_block_2, mem_size, mem_start, mem_start_block_1, &
         mem_start_block_2, mo_bsizes_1, mo_bsizes_2
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: bounds
      INTEGER, DIMENSION(2)                              :: pdims_2d
      INTEGER, DIMENSION(3)                              :: pdims
      LOGICAL                                            :: do_initialize
      REAL(dp)                                           :: t1, t2
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_distribution_type)                      :: ks_dist
      TYPE(dbt_pgrid_type)                               :: pgrid, pgrid_2d
      TYPE(dbt_type)                                     :: ks_t, ks_t_mat, mo_coeff_t, &
                                                            mo_coeff_t_split
      TYPE(dbt_type), DIMENSION(1, 1)                    :: t_3c_int_mo_1, t_3c_int_mo_2

      CALL timeset(routineN, handle)

      CPASSERT(SIZE(ks_matrix, 2) == 1)

      unit_nr_dbcsr = ri_data%unit_nr_dbcsr

      IF (geometry_did_change) THEN
         CALL hfx_ri_pre_scf_mo(qs_env, ri_data, nspins)
      END IF

      nblks = dbt_get_num_blocks_total(ri_data%t_3c_int_ctr_1(1, 1))
      IF (nblks == 0) THEN
         CPABORT("3-center integrals are not available (first call requires geometry_did_change=.TRUE.)")
      END IF

      DO ispin = 1, nspins
         nblks = dbt_get_num_blocks_total(ri_data%t_2c_int(ispin, 1))
         IF (nblks == 0) THEN
            CPABORT("2-center integrals are not available (first call requires geometry_did_change=.TRUE.)")
         END IF
      END DO

      IF (.NOT. ALLOCATED(ri_data%t_3c_int_mo)) THEN
         do_initialize = .TRUE.
         CPASSERT(.NOT. ALLOCATED(ri_data%t_3c_ctr_RI))
         CPASSERT(.NOT. ALLOCATED(ri_data%t_3c_ctr_KS))
         CPASSERT(.NOT. ALLOCATED(ri_data%t_3c_ctr_KS_copy))
         ALLOCATE (ri_data%t_3c_int_mo(nspins, 1, 1))
         ALLOCATE (ri_data%t_3c_ctr_RI(nspins, 1, 1))
         ALLOCATE (ri_data%t_3c_ctr_KS(nspins, 1, 1))
         ALLOCATE (ri_data%t_3c_ctr_KS_copy(nspins, 1, 1))
      ELSE
         do_initialize = .FALSE.
      END IF

      CALL get_qs_env(qs_env, para_env=para_env)

      ALLOCATE (bounds(2, 1))

      CALL dbcsr_get_info(ks_matrix(1, 1)%matrix, distribution=ks_dist)
      CALL dbcsr_distribution_get(ks_dist, group=comm_2d, nprows=pdims_2d(1), npcols=pdims_2d(2))

      pgrid_2d = dbt_nd_mp_comm(comm_2d, [1], [2], pdims_2d=pdims_2d)

      CALL create_2c_tensor(ks_t, dist1, dist2, pgrid_2d, ri_data%bsizes_AO_fit, ri_data%bsizes_AO_fit, &
                            name="(AO | AO)")

      DEALLOCATE (dist1, dist2)

      CALL mp_sync(para_env%group)
      t1 = m_walltime()

      DO ispin = 1, nspins

         CALL dbcsr_get_info(mo_coeff(ispin), nfullcols_total=n_mos)
         ALLOCATE (mo_bsizes_2(n_mos))
         mo_bsizes_2 = 1

         CALL create_tensor_batches(mo_bsizes_2, ri_data%n_mem, mem_start, mem_end, &
                                    mem_start_block_2, mem_end_block_2)
         n_mem = ri_data%n_mem
         ALLOCATE (mem_size(n_mem))

         DO i_mem = 1, n_mem
            bsize = SUM(mo_bsizes_2(mem_start_block_2(i_mem):mem_end_block_2(i_mem)))
            mem_size(i_mem) = bsize
         END DO

         CALL split_block_sizes(mem_size, mo_bsizes_1, ri_data%max_bsize_MO)
         ALLOCATE (mem_start_block_1(n_mem))
         ALLOCATE (mem_end_block_1(n_mem))
         nblock = SIZE(mo_bsizes_1)
         iblock = 0
         DO i_mem = 1, n_mem
            bsum = 0
            DO
               iblock = iblock + 1
               CPASSERT(iblock <= nblock)
               bsum = bsum + mo_bsizes_1(iblock)
               IF (bsum == mem_size(i_mem)) THEN
                  IF (i_mem == 1) THEN
                     mem_start_block_1(i_mem) = 1
                  ELSE
                     mem_start_block_1(i_mem) = mem_end_block_1(i_mem - 1) + 1
                  END IF
                  mem_end_block_1(i_mem) = iblock
                  EXIT
               END IF
            END DO
         END DO

         ALLOCATE (batch_ranges_1(ri_data%n_mem + 1))
         batch_ranges_1(:ri_data%n_mem) = mem_start_block_1(:)
         batch_ranges_1(ri_data%n_mem + 1) = mem_end_block_1(ri_data%n_mem) + 1

         ALLOCATE (batch_ranges_2(ri_data%n_mem + 1))
         batch_ranges_2(:ri_data%n_mem) = mem_start_block_2(:)
         batch_ranges_2(ri_data%n_mem + 1) = mem_end_block_2(ri_data%n_mem) + 1

         CALL mp_environ(nproc, iproc, para_env%group)

         CALL create_3c_tensor(t_3c_int_mo_1(1, 1), dist1, dist2, dist3, ri_data%pgrid_1, &
                               ri_data%bsizes_AO_split, ri_data%bsizes_RI_split, mo_bsizes_1, &
                               [1, 2], [3], &
                               name="(AO RI | MO)")

         DEALLOCATE (dist1, dist2, dist3)

         CALL create_3c_tensor(t_3c_int_mo_2(1, 1), dist1, dist2, dist3, ri_data%pgrid_2, &
                               mo_bsizes_1, ri_data%bsizes_RI_split, ri_data%bsizes_AO_split, &
                               [1], [2, 3], &
                               name="(MO | RI AO)")

         DEALLOCATE (dist1, dist2, dist3)

         CALL create_2c_tensor(mo_coeff_t_split, dist1, dist2, pgrid_2d, ri_data%bsizes_AO_split, mo_bsizes_1, &
                               name="(AO | MO)")

         DEALLOCATE (dist1, dist2)

         CPASSERT(homo(ispin)/ri_data%n_mem > 0)

         IF (do_initialize) THEN
            pdims(:) = 0

            CALL dbt_pgrid_create(para_env%group, pdims, pgrid, &
                                  tensor_dims=[SIZE(ri_data%bsizes_RI_fit), &
                                               (homo(ispin) - 1)/ri_data%n_mem + 1, &
                                               SIZE(ri_data%bsizes_AO_fit)])
            CALL create_3c_tensor(ri_data%t_3c_int_mo(ispin, 1, 1), dist1, dist2, dist3, pgrid, &
                                  ri_data%bsizes_RI_fit, mo_bsizes_2, ri_data%bsizes_AO_fit, &
                                  [1], [2, 3], &
                                  name="(RI | MO AO)")

            DEALLOCATE (dist1, dist2, dist3)

            CALL create_3c_tensor(ri_data%t_3c_ctr_KS(ispin, 1, 1), dist1, dist2, dist3, pgrid, &
                                  ri_data%bsizes_RI_fit, mo_bsizes_2, ri_data%bsizes_AO_fit, &
                                  [1, 2], [3], &
                                  name="(RI MO | AO)")
            DEALLOCATE (dist1, dist2, dist3)
            CALL dbt_pgrid_destroy(pgrid)

            CALL dbt_create(ri_data%t_3c_int_mo(ispin, 1, 1), ri_data%t_3c_ctr_RI(ispin, 1, 1), name="(RI | MO AO)")
            CALL dbt_create(ri_data%t_3c_ctr_KS(ispin, 1, 1), ri_data%t_3c_ctr_KS_copy(ispin, 1, 1))
         END IF

         CALL dbt_create(mo_coeff(ispin), mo_coeff_t, name="MO coeffs")
         CALL dbt_copy_matrix_to_tensor(mo_coeff(ispin), mo_coeff_t)
         CALL dbt_copy(mo_coeff_t, mo_coeff_t_split, move_data=.TRUE.)
         CALL dbt_filter(mo_coeff_t_split, ri_data%filter_eps_mo)
         CALL dbt_destroy(mo_coeff_t)

         CALL dbt_batched_contract_init(ks_t)
         CALL dbt_batched_contract_init(ri_data%t_3c_ctr_KS(ispin, 1, 1), batch_range_2=batch_ranges_2)
         CALL dbt_batched_contract_init(ri_data%t_3c_ctr_KS_copy(ispin, 1, 1), batch_range_2=batch_ranges_2)

         CALL dbt_batched_contract_init(ri_data%t_2c_int(ispin, 1))
         CALL dbt_batched_contract_init(ri_data%t_3c_int_mo(ispin, 1, 1), batch_range_2=batch_ranges_2)
         CALL dbt_batched_contract_init(ri_data%t_3c_ctr_RI(ispin, 1, 1), batch_range_2=batch_ranges_2)

         DO i_mem = 1, n_mem

            bounds(:, 1) = [mem_start(i_mem), mem_end(i_mem)]

            CALL dbt_batched_contract_init(mo_coeff_t_split)
            CALL dbt_batched_contract_init(ri_data%t_3c_int_ctr_1(1, 1))
            CALL dbt_batched_contract_init(t_3c_int_mo_1(1, 1), &
                                           batch_range_3=batch_ranges_1)
            CALL timeset(routineN//"_MOx3C_R", handle2)
            CALL dbt_contract(1.0_dp, mo_coeff_t_split, ri_data%t_3c_int_ctr_1(1, 1), &
                              0.0_dp, t_3c_int_mo_1(1, 1), &
                              contract_1=[1], notcontract_1=[2], &
                              contract_2=[3], notcontract_2=[1, 2], &
                              map_1=[3], map_2=[1, 2], &
                              bounds_2=bounds, &
                              filter_eps=ri_data%filter_eps_mo/2, &
                              unit_nr=unit_nr_dbcsr, &
                              move_data=.FALSE., &
                              flop=nflop)

            ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop

            CALL timestop(handle2)
            CALL dbt_batched_contract_finalize(mo_coeff_t_split)
            CALL dbt_batched_contract_finalize(ri_data%t_3c_int_ctr_1(1, 1))
            CALL dbt_batched_contract_finalize(t_3c_int_mo_1(1, 1))

            CALL timeset(routineN//"_copy_1", handle2)
            CALL dbt_copy(t_3c_int_mo_1(1, 1), ri_data%t_3c_int_mo(ispin, 1, 1), order=[3, 1, 2], move_data=.TRUE.)
            CALL timestop(handle2)

            CALL dbt_batched_contract_init(mo_coeff_t_split)
            CALL dbt_batched_contract_init(ri_data%t_3c_int_ctr_2(1, 1))
            CALL dbt_batched_contract_init(t_3c_int_mo_2(1, 1), &
                                           batch_range_1=batch_ranges_1)

            CALL timeset(routineN//"_MOx3C_L", handle2)
            CALL dbt_contract(1.0_dp, mo_coeff_t_split, ri_data%t_3c_int_ctr_2(1, 1), &
                              0.0_dp, t_3c_int_mo_2(1, 1), &
                              contract_1=[1], notcontract_1=[2], &
                              contract_2=[1], notcontract_2=[2, 3], &
                              map_1=[1], map_2=[2, 3], &
                              bounds_2=bounds, &
                              filter_eps=ri_data%filter_eps_mo/2, &
                              unit_nr=unit_nr_dbcsr, &
                              move_data=.FALSE., &
                              flop=nflop)

            ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop

            CALL timestop(handle2)

            CALL dbt_batched_contract_finalize(mo_coeff_t_split)
            CALL dbt_batched_contract_finalize(ri_data%t_3c_int_ctr_2(1, 1))
            CALL dbt_batched_contract_finalize(t_3c_int_mo_2(1, 1))

            CALL timeset(routineN//"_copy_1", handle2)
            CALL dbt_copy(t_3c_int_mo_2(1, 1), ri_data%t_3c_int_mo(ispin, 1, 1), order=[2, 1, 3], &
                          summation=.TRUE., move_data=.TRUE.)

            CALL dbt_filter(ri_data%t_3c_int_mo(ispin, 1, 1), ri_data%filter_eps_mo)
            CALL timestop(handle2)

            CALL timeset(routineN//"_RIx3C", handle2)

            CALL dbt_contract(1.0_dp, ri_data%t_2c_int(ispin, 1), ri_data%t_3c_int_mo(ispin, 1, 1), &
                              0.0_dp, ri_data%t_3c_ctr_RI(ispin, 1, 1), &
                              contract_1=[1], notcontract_1=[2], &
                              contract_2=[1], notcontract_2=[2, 3], &
                              map_1=[1], map_2=[2, 3], filter_eps=ri_data%filter_eps, &
                              unit_nr=unit_nr_dbcsr, &
                              flop=nflop)

            ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop

            CALL timestop(handle2)

            CALL timeset(routineN//"_copy_2", handle2)

            ! note: this copy should not involve communication (same block sizes, same 3d distribution on same process grid)
            CALL dbt_copy(ri_data%t_3c_ctr_RI(ispin, 1, 1), ri_data%t_3c_ctr_KS(ispin, 1, 1), move_data=.TRUE.)
            CALL dbt_copy(ri_data%t_3c_ctr_KS(ispin, 1, 1), ri_data%t_3c_ctr_KS_copy(ispin, 1, 1))
            CALL timestop(handle2)

            CALL timeset(routineN//"_3Cx3C", handle2)
            CALL dbt_contract(-fac, ri_data%t_3c_ctr_KS(ispin, 1, 1), ri_data%t_3c_ctr_KS_copy(ispin, 1, 1), &
                              1.0_dp, ks_t, &
                              contract_1=[1, 2], notcontract_1=[3], &
                              contract_2=[1, 2], notcontract_2=[3], &
                              map_1=[1], map_2=[2], filter_eps=ri_data%filter_eps/n_mem, &
                              unit_nr=unit_nr_dbcsr, move_data=.TRUE., &
                              flop=nflop)

            ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop

            CALL timestop(handle2)
         END DO

         CALL dbt_batched_contract_finalize(ks_t)
         CALL dbt_batched_contract_finalize(ri_data%t_3c_ctr_KS(ispin, 1, 1))
         CALL dbt_batched_contract_finalize(ri_data%t_3c_ctr_KS_copy(ispin, 1, 1))

         CALL dbt_batched_contract_finalize(ri_data%t_2c_int(ispin, 1))
         CALL dbt_batched_contract_finalize(ri_data%t_3c_int_mo(ispin, 1, 1))
         CALL dbt_batched_contract_finalize(ri_data%t_3c_ctr_RI(ispin, 1, 1))

         CALL dbt_destroy(t_3c_int_mo_1(1, 1))
         CALL dbt_destroy(t_3c_int_mo_2(1, 1))
         CALL dbt_clear(ri_data%t_3c_int_mo(ispin, 1, 1))

         CALL dbt_destroy(mo_coeff_t_split)

         CALL dbt_filter(ks_t, ri_data%filter_eps)

         CALL dbt_create(ks_matrix(ispin, 1)%matrix, ks_t_mat)
         CALL dbt_copy(ks_t, ks_t_mat, move_data=.TRUE.)
         CALL dbt_copy_tensor_to_matrix(ks_t_mat, ks_matrix(ispin, 1)%matrix, summation=.TRUE.)
         CALL dbt_destroy(ks_t_mat)

         DEALLOCATE (mem_end, mem_start, mo_bsizes_2, mem_size, mem_start_block_1, mem_end_block_1, &
                     mem_start_block_2, mem_end_block_2, batch_ranges_1, batch_ranges_2)

      END DO

      CALL dbt_pgrid_destroy(pgrid_2d)
      CALL dbt_destroy(ks_t)

      CALL mp_sync(para_env%group)
      t2 = m_walltime()

      ri_data%dbcsr_time = ri_data%dbcsr_time + t2 - t1

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief Calculate Fock (AKA Kohn-Sham) matrix in rho flavor
!>
!> M(mu, lambda, R) = sum_{nu} int_3c(mu, nu, R) P(nu, lambda)
!> KS(mu, lambda) = sum_{nu,R} B(mu, nu, R) M(lambda, nu, R)
!> \param qs_env ...
!> \param ri_data ...
!> \param ks_matrix ...
!> \param rho_ao ...
!> \param geometry_did_change ...
!> \param nspins ...
!> \param fac ...
! **************************************************************************************************
   SUBROUTINE hfx_ri_update_ks_Pmat(qs_env, ri_data, ks_matrix, rho_ao, &
                                    geometry_did_change, nspins, fac)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(dbcsr_p_type), DIMENSION(:, :)                :: ks_matrix, rho_ao
      LOGICAL, INTENT(IN)                                :: geometry_did_change
      INTEGER, INTENT(IN)                                :: nspins
      REAL(dp), INTENT(IN)                               :: fac

      CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_ri_update_ks_Pmat'

      INTEGER                                            :: handle, handle2, i_mem, ispin, j_mem, &
                                                            n_mem, n_mem_RI, unit_nr, unit_nr_dbcsr
      INTEGER(int_8)                                     :: flops_ks_max, flops_p_max, nblks, nflop, &
                                                            nze, nze_3c, nze_3c_1, nze_3c_2, &
                                                            nze_ks, nze_rho
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: batch_ranges_AO, batch_ranges_RI, dist1, &
                                                            dist2
      INTEGER, DIMENSION(2, 1)                           :: bounds_i
      INTEGER, DIMENSION(2, 2)                           :: bounds_ij, bounds_j
      INTEGER, DIMENSION(3)                              :: dims_3c
      REAL(dp)                                           :: memory_3c, occ, occ_3c, occ_3c_1, &
                                                            occ_3c_2, occ_ks, occ_rho, t1, t2, &
                                                            unused
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbt_type)                                     :: ks_t, ks_tmp, rho_ao_t, rho_ao_tmp, &
                                                            t_3c_1, t_3c_3, tensor_old

      CALL timeset(routineN, handle)

      NULLIFY (para_env)

      ! get a useful output_unit
      unit_nr_dbcsr = ri_data%unit_nr_dbcsr
      unit_nr = ri_data%unit_nr

      CALL get_qs_env(qs_env, para_env=para_env)

      CPASSERT(SIZE(ks_matrix, 2) == 1)

      IF (geometry_did_change) THEN
         CALL hfx_ri_pre_scf_Pmat(qs_env, ri_data)
      END IF

      nblks = dbt_get_num_blocks_total(ri_data%t_3c_int_ctr_2(1, 1))
      IF (nblks == 0) THEN
         CPABORT("3-center integrals are not available (first call requires geometry_did_change=.TRUE.)")
      END IF

      n_mem = ri_data%n_mem
      n_mem_RI = ri_data%n_mem_RI

      CALL dbt_create(ks_matrix(1, 1)%matrix, ks_tmp)
      CALL dbt_create(rho_ao(1, 1)%matrix, rho_ao_tmp)

      CALL create_2c_tensor(rho_ao_t, dist1, dist2, ri_data%pgrid_2d, &
                            ri_data%bsizes_AO_split, ri_data%bsizes_AO_split, &
                            name="(AO | AO)")
      DEALLOCATE (dist1, dist2)

      CALL create_2c_tensor(ks_t, dist1, dist2, ri_data%pgrid_2d, &
                            ri_data%bsizes_AO_split, ri_data%bsizes_AO_split, &
                            name="(AO | AO)")
      DEALLOCATE (dist1, dist2)

      CALL dbt_create(ri_data%t_3c_int_ctr_2(1, 1), t_3c_1)
      CALL dbt_create(ri_data%t_3c_int_ctr_1(1, 1), t_3c_3)

      CALL mp_sync(para_env%group)
      t1 = m_walltime()

      flops_ks_max = 0; flops_p_max = 0

      ALLOCATE (batch_ranges_RI(ri_data%n_mem_RI + 1))
      ALLOCATE (batch_ranges_AO(ri_data%n_mem + 1))
      batch_ranges_RI(:ri_data%n_mem_RI) = ri_data%starts_array_RI_mem_block(:)
      batch_ranges_RI(ri_data%n_mem_RI + 1) = ri_data%ends_array_RI_mem_block(ri_data%n_mem_RI) + 1
      batch_ranges_AO(:ri_data%n_mem) = ri_data%starts_array_mem_block(:)
      batch_ranges_AO(ri_data%n_mem + 1) = ri_data%ends_array_mem_block(ri_data%n_mem) + 1

      memory_3c = 0.0_dp
      DO ispin = 1, nspins

         CALL get_tensor_occupancy(ri_data%t_3c_int_ctr_2(1, 1), nze_3c, occ_3c)

         nze_rho = 0
         occ_rho = 0.0_dp
         nze_3c_1 = 0
         occ_3c_1 = 0.0_dp
         nze_3c_2 = 0
         occ_3c_2 = 0.0_dp

         CALL dbt_copy_matrix_to_tensor(rho_ao(ispin, 1)%matrix, rho_ao_tmp)
         CALL dbt_copy(rho_ao_tmp, rho_ao_t, move_data=.TRUE.)

         CALL get_tensor_occupancy(rho_ao_t, nze_rho, occ_rho)

         CALL dbt_batched_contract_init(ri_data%t_3c_int_ctr_1(1, 1), batch_range_1=batch_ranges_AO, &
                                        batch_range_2=batch_ranges_RI)
         CALL dbt_batched_contract_init(t_3c_3, batch_range_1=batch_ranges_AO, batch_range_2=batch_ranges_RI)

         CALL dbt_create(ri_data%t_3c_int_ctr_1(1, 1), tensor_old)

         DO i_mem = 1, n_mem

            CALL dbt_batched_contract_init(rho_ao_t)
            CALL dbt_batched_contract_init(ri_data%t_3c_int_ctr_2(1, 1), batch_range_2=batch_ranges_RI, &
                                           batch_range_3=batch_ranges_AO)
            CALL dbt_batched_contract_init(t_3c_1, batch_range_2=batch_ranges_RI, batch_range_3=batch_ranges_AO)
            DO j_mem = 1, n_mem_RI

               CALL timeset(routineN//"_Px3C", handle2)

               CALL dbt_get_info(t_3c_1, nfull_total=dims_3c)
               bounds_i(:, 1) = [ri_data%starts_array_mem(i_mem), ri_data%ends_array_mem(i_mem)]
               bounds_j(:, 1) = [1, dims_3c(1)]
               bounds_j(:, 2) = [ri_data%starts_array_RI_mem(j_mem), ri_data%ends_array_RI_mem(j_mem)]

               CALL dbt_contract(1.0_dp, rho_ao_t, ri_data%t_3c_int_ctr_2(1, 1), &
                                 0.0_dp, t_3c_1, &
                                 contract_1=[2], notcontract_1=[1], &
                                 contract_2=[3], notcontract_2=[1, 2], &
                                 map_1=[3], map_2=[1, 2], filter_eps=ri_data%filter_eps, &
                                 bounds_2=bounds_i, &
                                 bounds_3=bounds_j, &
                                 unit_nr=unit_nr_dbcsr, &
                                 flop=nflop)

               CALL timestop(handle2)

               ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop

               CALL get_tensor_occupancy(t_3c_1, nze, occ)
               nze_3c_1 = nze_3c_1 + nze
               occ_3c_1 = occ_3c_1 + occ

               CALL timeset(routineN//"_copy_2", handle2)
               CALL dbt_copy(t_3c_1, t_3c_3, order=[3, 2, 1], move_data=.TRUE.)
               CALL timestop(handle2)

               bounds_ij(:, 1) = [ri_data%starts_array_mem(i_mem), ri_data%ends_array_mem(i_mem)]
               bounds_ij(:, 2) = [ri_data%starts_array_RI_mem(j_mem), ri_data%ends_array_RI_mem(j_mem)]

               CALL decompress_tensor(tensor_old, ri_data%blk_indices(i_mem, j_mem)%ind, &
                                      ri_data%store_3c(i_mem, j_mem), ri_data%filter_eps_storage)

               CALL dbt_copy(tensor_old, ri_data%t_3c_int_ctr_1(1, 1), move_data=.TRUE.)

               CALL get_tensor_occupancy(ri_data%t_3c_int_ctr_1(1, 1), nze, occ)
               nze_3c_2 = nze_3c_2 + nze
               occ_3c_2 = occ_3c_2 + occ
               CALL timeset(routineN//"_KS", handle2)
               CALL dbt_batched_contract_init(ks_t)
               CALL dbt_contract(-fac, ri_data%t_3c_int_ctr_1(1, 1), t_3c_3, &
                                 1.0_dp, ks_t, &
                                 contract_1=[1, 2], notcontract_1=[3], &
                                 contract_2=[1, 2], notcontract_2=[3], &
                                 map_1=[1], map_2=[2], filter_eps=ri_data%filter_eps/n_mem, &
                                 bounds_1=bounds_ij, &
                                 unit_nr=unit_nr_dbcsr, &
                                 flop=nflop, move_data=.TRUE.)

               CALL dbt_batched_contract_finalize(ks_t, unit_nr=unit_nr_dbcsr)
               CALL timestop(handle2)

               ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop

            END DO
            CALL dbt_batched_contract_finalize(rho_ao_t, unit_nr=unit_nr_dbcsr)
            CALL dbt_batched_contract_finalize(ri_data%t_3c_int_ctr_2(1, 1))
            CALL dbt_batched_contract_finalize(t_3c_1)
         END DO
         CALL dbt_batched_contract_finalize(ri_data%t_3c_int_ctr_1(1, 1))
         CALL dbt_batched_contract_finalize(t_3c_3)

         DO i_mem = 1, n_mem
            DO j_mem = 1, n_mem_RI
               ASSOCIATE (blk_indices => ri_data%blk_indices(i_mem, j_mem), t_3c => ri_data%t_3c_int_ctr_1(1, 1))
                  CALL decompress_tensor(tensor_old, blk_indices%ind, &
                                         ri_data%store_3c(i_mem, j_mem), ri_data%filter_eps_storage)
                  CALL dbt_copy(tensor_old, t_3c, move_data=.TRUE.)

                  unused = 0
                  CALL compress_tensor(t_3c, blk_indices%ind, ri_data%store_3c(i_mem, j_mem), &
                                       ri_data%filter_eps_storage, unused)
               END ASSOCIATE
            END DO
         END DO

         CALL dbt_destroy(tensor_old)

         CALL dbt_clear(rho_ao_t)
         CALL get_tensor_occupancy(ks_t, nze_ks, occ_ks)

         CALL dbt_copy(ks_t, ks_tmp)
         CALL dbt_clear(ks_t)
         CALL dbt_copy_tensor_to_matrix(ks_tmp, ks_matrix(ispin, 1)%matrix, summation=.TRUE.)
         CALL dbt_clear(ks_tmp)

         IF (unit_nr > 0 .AND. geometry_did_change) THEN
            WRITE (unit_nr, '(T6,A,T63,ES7.1,1X,A1,1X,F7.3,A1)') &
               'Occupancy of density matrix P:', REAL(nze_rho, dp), '/', occ_rho*100, '%'
            WRITE (unit_nr, '(T6,A,T63,ES7.1,1X,A1,1X,F7.3,A1)') &
               'Occupancy of 3c ints:', REAL(nze_3c, dp), '/', occ_3c*100, '%'
            WRITE (unit_nr, '(T6,A,T63,ES7.1,1X,A1,1X,F7.3,A1)') &
               'Occupancy after contraction with K:', REAL(nze_3c_2, dp), '/', occ_3c_2*100, '%'
            WRITE (unit_nr, '(T6,A,T63,ES7.1,1X,A1,1X,F7.3,A1)') &
               'Occupancy after contraction with P:', REAL(nze_3c_1, dp), '/', occ_3c_1*100, '%'
            WRITE (unit_nr, '(T6,A,T63,ES7.1,1X,A1,1X,F7.3,A1)') &
               'Occupancy of Kohn-Sham matrix:', REAL(nze_ks, dp), '/', occ_ks*100, '%'
         END IF

      END DO

      CALL mp_sync(para_env%group)
      t2 = m_walltime()

      ri_data%dbcsr_time = ri_data%dbcsr_time + t2 - t1

      CALL dbt_destroy(t_3c_1)
      CALL dbt_destroy(t_3c_3)

      CALL dbt_destroy(rho_ao_t)
      CALL dbt_destroy(rho_ao_tmp)
      CALL dbt_destroy(ks_t)
      CALL dbt_destroy(ks_tmp)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief Implementation based on the MO flavor
!> \param qs_env ...
!> \param ri_data ...
!> \param nspins ...
!> \param hf_fraction ...
!> \param mo_coeff ...
!> \param use_virial ...
!> \note There is no response code for forces with the MO flavor
! **************************************************************************************************
   SUBROUTINE hfx_ri_forces_mo(qs_env, ri_data, nspins, hf_fraction, mo_coeff, use_virial)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      INTEGER, INTENT(IN)                                :: nspins
      REAL(dp), INTENT(IN)                               :: hf_fraction
      TYPE(dbcsr_type), DIMENSION(:), INTENT(IN)         :: mo_coeff
      LOGICAL, INTENT(IN), OPTIONAL                      :: use_virial

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'hfx_ri_forces_mo'

      INTEGER :: handle, i_mem, i_xyz, ispin, j_mem, j_xyz, k_mem, k_xyz, n_mem, n_mem_input, &
         n_mem_input_RI, n_mem_RI, n_mem_RI_fit, n_mos, natom, unit_nr_dbcsr
      INTEGER(int_8)                                     :: nflop
      INTEGER, ALLOCATABLE, DIMENSION(:) :: atom_of_kind, batch_blk_end, batch_blk_start, &
         batch_end, batch_end_RI, batch_end_RI_fit, batch_ranges, batch_ranges_RI, &
         batch_ranges_RI_fit, batch_start, batch_start_RI, batch_start_RI_fit, bsizes_MO, dist1, &
         dist2, dist3, idx_to_at_AO, idx_to_at_RI, kind_of
      INTEGER, DIMENSION(2, 1)                           :: bounds_ctr_1d
      INTEGER, DIMENSION(2, 2)                           :: bounds_ctr_2d
      INTEGER, DIMENSION(3)                              :: pdims
      LOGICAL                                            :: use_virial_prv
      REAL(dp)                                           :: pref, spin_fac, t1, t2
      REAL(dp), DIMENSION(3, 3)                          :: work_virial
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dbt_pgrid_type)                               :: pgrid_1, pgrid_2
      TYPE(dbt_type) :: t_2c_RI, t_2c_RI_inv, t_2c_RI_met, t_2c_RI_PQ, t_2c_tmp, t_3c_0, t_3c_1, &
         t_3c_2, t_3c_3, t_3c_4, t_3c_5, t_3c_6, t_3c_ao_ri_ao, t_3c_ao_ri_mo, t_3c_desymm, &
         t_3c_mo_ri_ao, t_3c_mo_ri_mo, t_3c_ri_ao_ao, t_3c_RI_ctr, t_3c_ri_mo_mo, &
         t_3c_ri_mo_mo_fit, t_3c_work, t_mo_coeff, t_mo_cpy
      TYPE(dbt_type), DIMENSION(3) :: t_2c_der_metric, t_2c_der_RI, t_2c_MO_AO, t_2c_MO_AO_ctr, &
         t_3c_der_AO, t_3c_der_AO_ctr_1, t_3c_der_RI, t_3c_der_RI_ctr_1, t_3c_der_RI_ctr_2, &
         t_3c_tmp_RI
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(virial_type), POINTER                         :: virial

      ! 1) Precompute the derivatives that are needed (3c, 3c RI and metric)
      ! 2) Go over batched of occupied MOs so as to save memory and optimize contractions
      ! 3) Contract all 3c integrals and derivatives with MO coeffs
      ! 4) Contract relevant quantities with the inverse 2c RI (metric or pot)
      ! 5) First force contribution with the 2c RI derivative d/dx (Q|R)
      ! 6) If metric, do the additional contraction  with S_pq^-1 (Q|R)
      ! 7) Do the force contribution due to 3c integrals (a'b|P) and (ab|P')
      ! 8) If metric, do the last force contribution due to d/dx S^-1 (First contract (ab|P), then S^-1)

      use_virial_prv = .FALSE.
      IF (PRESENT(use_virial)) use_virial_prv = use_virial

      unit_nr_dbcsr = ri_data%unit_nr_dbcsr

      CALL get_qs_env(qs_env, natom=natom, particle_set=particle_set, &
                      atomic_kind_set=atomic_kind_set, virial=virial, &
                      cell=cell, force=force, matrix_s=matrix_s, &
                      para_env=para_env)

      pdims(:) = 0
      CALL dbt_pgrid_create(para_env%group, pdims, pgrid_1, tensor_dims=[SIZE(ri_data%bsizes_AO_split), &
                                                                         SIZE(ri_data%bsizes_RI_split), &
                                                                         SIZE(ri_data%bsizes_AO_split)])
      pdims(:) = 0
      CALL dbt_pgrid_create(para_env%group, pdims, pgrid_2, tensor_dims=[SIZE(ri_data%bsizes_RI_split), &
                                                                         SIZE(ri_data%bsizes_AO_split), &
                                                                         SIZE(ri_data%bsizes_AO_split)])

      CALL create_3c_tensor(t_3c_ao_ri_ao, dist1, dist2, dist3, pgrid_1, &
                            ri_data%bsizes_AO_split, ri_data%bsizes_RI_split, ri_data%bsizes_AO_split, &
                            [1, 2], [3], name="(AO RI | AO)")
      DEALLOCATE (dist1, dist2, dist3)
      CALL create_3c_tensor(t_3c_ri_ao_ao, dist1, dist2, dist3, pgrid_2, &
                            ri_data%bsizes_RI_split, ri_data%bsizes_AO_split, ri_data%bsizes_AO_split, &
                            [1], [2, 3], name="(RI | AO AO)")
      DEALLOCATE (dist1, dist2, dist3)

      ! 1) Precompute the derivatives
      CALL precalc_derivatives(t_3c_tmp_RI, t_3c_der_AO, t_2c_der_RI, t_2c_der_metric, &
                               t_3c_ri_ao_ao, t_3c_ao_ri_ao, ri_data, qs_env)
      DO i_xyz = 1, 3
         CALL dbt_create(t_3c_ao_ri_ao, t_3c_der_RI(i_xyz))
         CALL dbt_copy(t_3c_tmp_RI(i_xyz), t_3c_der_RI(i_xyz), order=[2, 1, 3], move_data=.TRUE.)
         CALL dbt_destroy(t_3c_tmp_RI(i_xyz))
      END DO

      ! Get the 3c integrals (desymmetrized)
      CALL dbt_create(t_3c_ao_ri_ao, t_3c_desymm)
      CALL dbt_copy(ri_data%t_3c_int_ctr_1(1, 1), t_3c_desymm)
      CALL dbt_copy(ri_data%t_3c_int_ctr_1(1, 1), t_3c_desymm, order=[3, 2, 1], &
                    summation=.TRUE., move_data=.TRUE.)

      CALL dbt_destroy(t_3c_ao_ri_ao)
      CALL dbt_destroy(t_3c_ri_ao_ao)

      ! Some utilities
      spin_fac = 0.5_dp
      IF (nspins == 2) spin_fac = 1.0_dp

      ALLOCATE (idx_to_at_RI(SIZE(ri_data%bsizes_RI_split)))
      CALL get_idx_to_atom(idx_to_at_RI, ri_data%bsizes_RI_split, ri_data%bsizes_RI)

      ALLOCATE (idx_to_at_AO(SIZE(ri_data%bsizes_AO_split)))
      CALL get_idx_to_atom(idx_to_at_AO, ri_data%bsizes_AO_split, ri_data%bsizes_AO)

      ALLOCATE (atom_of_kind(natom), kind_of(natom))
      CALL get_atomic_kind_set(atomic_kind_set, kind_of=kind_of, atom_of_kind=atom_of_kind)

      ! 2-center RI tensors
      CALL create_2c_tensor(t_2c_RI, dist1, dist2, ri_data%pgrid_2d, &
                            ri_data%bsizes_RI_split, ri_data%bsizes_RI_split, name="(RI | RI)")
      DEALLOCATE (dist1, dist2)

      CALL create_2c_tensor(t_2c_RI_PQ, dist1, dist2, ri_data%pgrid_2d, &
                            ri_data%bsizes_RI_fit, ri_data%bsizes_RI_fit, name="(RI | RI)")
      DEALLOCATE (dist1, dist2)

      IF (.NOT. ri_data%same_op) THEN
         !precompute the (P|Q)*S^-1 product
         CALL dbt_create(t_2c_RI_PQ, t_2c_RI_inv)
         CALL dbt_create(t_2c_RI_PQ, t_2c_RI_met)
         CALL dbt_create(ri_data%t_2c_inv(1, 1), t_2c_tmp)

         CALL dbt_contract(1.0_dp, ri_data%t_2c_inv(1, 1), ri_data%t_2c_pot(1, 1), &
                           0.0_dp, t_2c_tmp, &
                           contract_1=[2], notcontract_1=[1], &
                           contract_2=[1], notcontract_2=[2], &
                           map_1=[1], map_2=[2], filter_eps=ri_data%filter_eps, &
                           unit_nr=unit_nr_dbcsr, flop=nflop)
         ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop

         CALL dbt_copy(t_2c_tmp, t_2c_RI_inv, move_data=.TRUE.)
         CALL dbt_destroy(t_2c_tmp)
      END IF

      !3 loops in MO force evaluations. To be consistent with input MEMORY_CUT, need to take the cubic root
      !No need to cut memory further because SCF tensors alrady dense
      n_mem_input = FLOOR((ri_data%n_mem_input - 0.1_dp)**(1._dp/3._dp)) + 1
      n_mem_input_RI = FLOOR((ri_data%n_mem_input - 0.1_dp)/n_mem_input**2) + 1

      !batches on RI_split and RI_fit blocks
      n_mem_RI = n_mem_input_RI
      CALL create_tensor_batches(ri_data%bsizes_RI_split, n_mem_RI, batch_start_RI, batch_end_RI, &
                                 batch_blk_start, batch_blk_end)
      ALLOCATE (batch_ranges_RI(n_mem_RI + 1))
      batch_ranges_RI(1:n_mem_RI) = batch_blk_start(1:n_mem_RI)
      batch_ranges_RI(n_mem_RI + 1) = batch_blk_end(n_mem_RI) + 1
      DEALLOCATE (batch_blk_start, batch_blk_end)

      n_mem_RI_fit = n_mem_input_RI
      CALL create_tensor_batches(ri_data%bsizes_RI_fit, n_mem_RI_fit, batch_start_RI_fit, batch_end_RI_fit, &
                                 batch_blk_start, batch_blk_end)
      ALLOCATE (batch_ranges_RI_fit(n_mem_RI_fit + 1))
      batch_ranges_RI_fit(1:n_mem_RI_fit) = batch_blk_start(1:n_mem_RI_fit)
      batch_ranges_RI_fit(n_mem_RI_fit + 1) = batch_blk_end(n_mem_RI_fit) + 1
      DEALLOCATE (batch_blk_start, batch_blk_end)

      DO ispin = 1, nspins
         work_virial = 0.0_dp

         ! 2 )Prepare the batches for this spin
         CALL dbcsr_get_info(mo_coeff(ispin), nfullcols_total=n_mos)
         !note: optimized GPU block size for SCF is 64x1x64. Here we do 8x8x64
         CALL split_block_sizes([n_mos], bsizes_MO, max_size=FLOOR(SQRT(ri_data%max_bsize_MO - 0.1)) + 1)

         !batching on MO blocks
         n_mem = n_mem_input
         CALL create_tensor_batches(bsizes_MO, n_mem, batch_start, batch_end, &
                                    batch_blk_start, batch_blk_end)
         ALLOCATE (batch_ranges(n_mem + 1))
         batch_ranges(1:n_mem) = batch_blk_start(1:n_mem)
         batch_ranges(n_mem + 1) = batch_blk_end(n_mem) + 1
         DEALLOCATE (batch_blk_start, batch_blk_end)

         ! Initialize the different tensors needed (Note: keep MO coeffs as (MO | AO) for less transpose)
         CALL create_2c_tensor(t_mo_coeff, dist1, dist2, ri_data%pgrid_2d, bsizes_MO, &
                               ri_data%bsizes_AO_split, name="MO coeffs")
         DEALLOCATE (dist1, dist2)
         CALL dbt_create(mo_coeff(ispin), t_2c_tmp, name="MO coeffs")
         CALL dbt_copy_matrix_to_tensor(mo_coeff(ispin), t_2c_tmp)
         CALL dbt_copy(t_2c_tmp, t_mo_coeff, order=[2, 1], move_data=.TRUE.)
         CALL dbt_destroy(t_2c_tmp)

         CALL dbt_create(t_mo_coeff, t_mo_cpy)
         CALL dbt_copy(t_mo_coeff, t_mo_cpy)
         DO i_xyz = 1, 3
            CALL dbt_create(t_mo_coeff, t_2c_MO_AO_ctr(i_xyz))
            CALL dbt_create(t_mo_coeff, t_2c_MO_AO(i_xyz))
         END DO

         CALL create_3c_tensor(t_3c_ao_ri_mo, dist1, dist2, dist3, pgrid_1, ri_data%bsizes_AO_split, &
                               ri_data%bsizes_RI_split, bsizes_MO, [1, 2], [3], name="(AO RI| MO)")
         DEALLOCATE (dist1, dist2, dist3)

         CALL dbt_create(t_3c_ao_ri_mo, t_3c_0)
         CALL dbt_destroy(t_3c_ao_ri_mo)

         CALL create_3c_tensor(t_3c_mo_ri_ao, dist1, dist2, dist3, pgrid_1, bsizes_MO, ri_data%bsizes_RI_split, &
                               ri_data%bsizes_AO_split, [1, 2], [3], name="(MO RI | AO)")
         DEALLOCATE (dist1, dist2, dist3)
         CALL dbt_create(t_3c_mo_ri_ao, t_3c_1)

         DO i_xyz = 1, 3
            CALL dbt_create(t_3c_mo_ri_ao, t_3c_der_RI_ctr_1(i_xyz))
            CALL dbt_create(t_3c_mo_ri_ao, t_3c_der_AO_ctr_1(i_xyz))
         END DO

         CALL create_3c_tensor(t_3c_mo_ri_mo, dist1, dist2, dist3, pgrid_1, bsizes_MO, &
                               ri_data%bsizes_RI_split, bsizes_MO, [1, 2], [3], name="(MO RI | MO)")
         DEALLOCATE (dist1, dist2, dist3)
         CALL dbt_create(t_3c_mo_ri_mo, t_3c_work)

         CALL create_3c_tensor(t_3c_ri_mo_mo, dist1, dist2, dist3, pgrid_2, ri_data%bsizes_RI_split, &
                               bsizes_MO, bsizes_MO, [1], [2, 3], name="(RI| MO MO)")
         DEALLOCATE (dist1, dist2, dist3)

         CALL dbt_create(t_3c_ri_mo_mo, t_3c_2)
         CALL dbt_create(t_3c_ri_mo_mo, t_3c_3)
         CALL dbt_create(t_3c_ri_mo_mo, t_3c_RI_ctr)
         DO i_xyz = 1, 3
            CALL dbt_create(t_3c_ri_mo_mo, t_3c_der_RI_ctr_2(i_xyz))
         END DO

         !Very large RI_fit blocks => new pgrid to make sure distribution is ideal
         pdims(:) = 0
         CALL create_3c_tensor(t_3c_ri_mo_mo_fit, dist1, dist2, dist3, pgrid_2, ri_data%bsizes_RI_fit, &
                               bsizes_MO, bsizes_MO, [1], [2, 3], name="(RI| MO MO)")
         DEALLOCATE (dist1, dist2, dist3)

         CALL dbt_create(t_3c_ri_mo_mo_fit, t_3c_4)
         CALL dbt_create(t_3c_ri_mo_mo_fit, t_3c_5)
         CALL dbt_create(t_3c_ri_mo_mo_fit, t_3c_6)

         CALL dbt_batched_contract_init(t_3c_desymm, batch_range_2=batch_ranges_RI)
         CALL dbt_batched_contract_init(t_3c_0, batch_range_2=batch_ranges_RI, batch_range_3=batch_ranges)

         DO i_xyz = 1, 3
            CALL dbt_batched_contract_init(t_3c_der_AO(i_xyz), batch_range_2=batch_ranges_RI)
            CALL dbt_batched_contract_init(t_3c_der_RI(i_xyz), batch_range_2=batch_ranges_RI)
         END DO

         CALL mp_sync(para_env%group)
         t1 = m_walltime()

         ! 2) Loop over batches
         DO i_mem = 1, n_mem

            bounds_ctr_1d(1, 1) = batch_start(i_mem)
            bounds_ctr_1d(2, 1) = batch_end(i_mem)

            bounds_ctr_2d(1, 1) = 1
            bounds_ctr_2d(2, 1) = SUM(ri_data%bsizes_AO)

            ! 3) Do the first AO to MO contraction here
            CALL timeset(routineN//"_AO2MO_1", handle)
            CALL dbt_batched_contract_init(t_mo_coeff)
            DO k_mem = 1, n_mem_RI
               bounds_ctr_2d(1, 2) = batch_start_RI(k_mem)
               bounds_ctr_2d(2, 2) = batch_end_RI(k_mem)

               CALL dbt_contract(1.0_dp, t_mo_coeff, t_3c_desymm, &
                                 1.0_dp, t_3c_0, &
                                 contract_1=[2], notcontract_1=[1], &
                                 contract_2=[3], notcontract_2=[1, 2], &
                                 map_1=[3], map_2=[1, 2], filter_eps=ri_data%filter_eps, &
                                 bounds_2=bounds_ctr_1d, &
                                 bounds_3=bounds_ctr_2d, &
                                 unit_nr=unit_nr_dbcsr, flop=nflop)
               ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop
            END DO
            CALL dbt_copy(t_3c_0, t_3c_1, order=[3, 2, 1], move_data=.TRUE.)

            DO i_xyz = 1, 3
               DO k_mem = 1, n_mem_RI
                  bounds_ctr_2d(1, 2) = batch_start_RI(k_mem)
                  bounds_ctr_2d(2, 2) = batch_end_RI(k_mem)

                  CALL dbt_contract(1.0_dp, t_mo_coeff, t_3c_der_AO(i_xyz), &
                                    1.0_dp, t_3c_0, &
                                    contract_1=[2], notcontract_1=[1], &
                                    contract_2=[3], notcontract_2=[1, 2], &
                                    map_1=[3], map_2=[1, 2], filter_eps=ri_data%filter_eps, &
                                    bounds_2=bounds_ctr_1d, &
                                    bounds_3=bounds_ctr_2d, &
                                    unit_nr=unit_nr_dbcsr, flop=nflop)
                  ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop
               END DO
               CALL dbt_copy(t_3c_0, t_3c_der_AO_ctr_1(i_xyz), order=[3, 2, 1], move_data=.TRUE.)

               DO k_mem = 1, n_mem_RI
                  bounds_ctr_2d(1, 2) = batch_start_RI(k_mem)
                  bounds_ctr_2d(2, 2) = batch_end_RI(k_mem)

                  CALL dbt_contract(1.0_dp, t_mo_coeff, t_3c_der_RI(i_xyz), &
                                    1.0_dp, t_3c_0, &
                                    contract_1=[2], notcontract_1=[1], &
                                    contract_2=[3], notcontract_2=[1, 2], &
                                    map_1=[3], map_2=[1, 2], filter_eps=ri_data%filter_eps, &
                                    bounds_2=bounds_ctr_1d, &
                                    bounds_3=bounds_ctr_2d, &
                                    unit_nr=unit_nr_dbcsr, flop=nflop)
                  ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop
               END DO
               CALL dbt_copy(t_3c_0, t_3c_der_RI_ctr_1(i_xyz), order=[3, 2, 1], move_data=.TRUE.)
            END DO
            CALL dbt_batched_contract_finalize(t_mo_coeff)
            CALL timestop(handle)

            CALL dbt_batched_contract_init(t_3c_1, batch_range_1=batch_ranges, batch_range_2=batch_ranges_RI)
            CALL dbt_batched_contract_init(t_3c_work, batch_range_1=batch_ranges, batch_range_2=batch_ranges_RI, &
                                           batch_range_3=batch_ranges)
            CALL dbt_batched_contract_init(t_3c_2, batch_range_2=batch_ranges, batch_range_3=batch_ranges)
            CALL dbt_batched_contract_init(t_3c_3, batch_range_1=batch_ranges_RI, &
                                           batch_range_2=batch_ranges, batch_range_3=batch_ranges)

            CALL dbt_batched_contract_init(t_3c_4, batch_range_1=batch_ranges_RI_fit, &
                                           batch_range_2=batch_ranges, batch_range_3=batch_ranges)
            CALL dbt_batched_contract_init(t_3c_5, batch_range_2=batch_ranges, batch_range_3=batch_ranges)

            DO i_xyz = 1, 3
               CALL dbt_batched_contract_init(t_3c_der_RI_ctr_1(i_xyz), batch_range_1=batch_ranges, &
                                              batch_range_2=batch_ranges_RI)
               CALL dbt_batched_contract_init(t_3c_der_AO_ctr_1(i_xyz), batch_range_1=batch_ranges, &
                                              batch_range_2=batch_ranges_RI)

            END DO

            IF (.NOT. ri_data%same_op) THEN
               CALL dbt_batched_contract_init(t_3c_6, batch_range_2=batch_ranges, batch_range_3=batch_ranges)
            END IF

            DO j_mem = 1, n_mem

               bounds_ctr_1d(1, 1) = batch_start(j_mem)
               bounds_ctr_1d(2, 1) = batch_end(j_mem)

               bounds_ctr_2d(1, 1) = batch_start(i_mem)
               bounds_ctr_2d(2, 1) = batch_end(i_mem)

               ! 3) Do the second AO to MO contraction here, followed by the S^-1 contraction
               CALL timeset(routineN//"_AO2MO_2", handle)
               CALL dbt_batched_contract_init(t_mo_coeff)
               DO k_mem = 1, n_mem_RI
                  bounds_ctr_2d(1, 2) = batch_start_RI(k_mem)
                  bounds_ctr_2d(2, 2) = batch_end_RI(k_mem)

                  CALL dbt_contract(1.0_dp, t_mo_coeff, t_3c_1, &
                                    1.0_dp, t_3c_work, &
                                    contract_1=[2], notcontract_1=[1], &
                                    contract_2=[3], notcontract_2=[1, 2], &
                                    map_1=[3], map_2=[1, 2], filter_eps=ri_data%filter_eps, &
                                    bounds_2=bounds_ctr_1d, &
                                    bounds_3=bounds_ctr_2d, &
                                    unit_nr=unit_nr_dbcsr, flop=nflop)
                  ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop
               END DO
               CALL dbt_batched_contract_finalize(t_mo_coeff)
               CALL timestop(handle)

               bounds_ctr_2d(1, 1) = batch_start(i_mem)
               bounds_ctr_2d(2, 1) = batch_end(i_mem)
               bounds_ctr_2d(1, 2) = batch_start(j_mem)
               bounds_ctr_2d(2, 2) = batch_end(j_mem)

               ! 4) Contract 3c MO integrals with S^-1 as well
               CALL timeset(routineN//"_2c_inv", handle)
               CALL dbt_copy(t_3c_work, t_3c_3, order=[2, 1, 3], move_data=.TRUE.)
               DO k_mem = 1, n_mem_RI
                  bounds_ctr_1d(1, 1) = batch_start_RI(k_mem)
                  bounds_ctr_1d(2, 1) = batch_end_RI(k_mem)

                  CALL dbt_batched_contract_init(ri_data%t_2c_inv(1, 1))
                  CALL dbt_contract(1.0_dp, ri_data%t_2c_inv(1, 1), t_3c_3, &
                                    1.0_dp, t_3c_2, &
                                    contract_1=[2], notcontract_1=[1], &
                                    contract_2=[1], notcontract_2=[2, 3], &
                                    map_1=[1], map_2=[2, 3], filter_eps=ri_data%filter_eps, &
                                    bounds_1=bounds_ctr_1d, &
                                    bounds_3=bounds_ctr_2d, &
                                    unit_nr=unit_nr_dbcsr, flop=nflop)
                  ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop
                  CALL dbt_batched_contract_finalize(ri_data%t_2c_inv(1, 1))
               END DO
               CALL dbt_copy(t_3c_ri_mo_mo, t_3c_3)
               CALL timestop(handle)

               !Only contract (ab|P') with MO coeffs since need AO rep for the force of (a'b|P)
               bounds_ctr_1d(1, 1) = batch_start(j_mem)
               bounds_ctr_1d(2, 1) = batch_end(j_mem)

               bounds_ctr_2d(1, 1) = batch_start(i_mem)
               bounds_ctr_2d(2, 1) = batch_end(i_mem)

               CALL timeset(routineN//"_AO2MO_2", handle)
               CALL dbt_batched_contract_init(t_mo_coeff)
               DO i_xyz = 1, 3
                  DO k_mem = 1, n_mem_RI
                     bounds_ctr_2d(1, 2) = batch_start_RI(k_mem)
                     bounds_ctr_2d(2, 2) = batch_end_RI(k_mem)

                     CALL dbt_contract(1.0_dp, t_mo_coeff, t_3c_der_RI_ctr_1(i_xyz), &
                                       1.0_dp, t_3c_work, &
                                       contract_1=[2], notcontract_1=[1], &
                                       contract_2=[3], notcontract_2=[1, 2], &
                                       map_1=[3], map_2=[1, 2], filter_eps=ri_data%filter_eps, &
                                       bounds_2=bounds_ctr_1d, &
                                       bounds_3=bounds_ctr_2d, &
                                       unit_nr=unit_nr_dbcsr, flop=nflop)
                     ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop
                  END DO
                  CALL dbt_copy(t_3c_work, t_3c_der_RI_ctr_2(i_xyz), order=[2, 1, 3], move_data=.TRUE.)
               END DO
               CALL dbt_batched_contract_finalize(t_mo_coeff)
               CALL timestop(handle)

               bounds_ctr_2d(1, 1) = batch_start(i_mem)
               bounds_ctr_2d(2, 1) = batch_end(i_mem)
               bounds_ctr_2d(1, 2) = batch_start(j_mem)
               bounds_ctr_2d(2, 2) = batch_end(j_mem)

               ! 5) Force due to d/dx (P|Q)
               CALL timeset(routineN//"_PQ_der", handle)
               CALL dbt_copy(t_3c_2, t_3c_4, move_data=.TRUE.)
               CALL dbt_copy(t_3c_4, t_3c_5)
               DO k_mem = 1, n_mem_RI_fit
                  bounds_ctr_1d(1, 1) = batch_start_RI_fit(k_mem)
                  bounds_ctr_1d(2, 1) = batch_end_RI_fit(k_mem)

                  CALL dbt_batched_contract_init(t_2c_RI_PQ)
                  CALL dbt_contract(1.0_dp, t_3c_4, t_3c_5, &
                                    1.0_dp, t_2c_RI_PQ, &
                                    contract_1=[2, 3], notcontract_1=[1], &
                                    contract_2=[2, 3], notcontract_2=[1], &
                                    bounds_1=bounds_ctr_2d, &
                                    bounds_2=bounds_ctr_1d, &
                                    map_1=[1], map_2=[2], filter_eps=ri_data%filter_eps, &
                                    unit_nr=unit_nr_dbcsr, flop=nflop)
                  ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop
                  CALL dbt_batched_contract_finalize(t_2c_RI_PQ)
               END DO
               CALL timestop(handle)

               ! 6) If metric, do the additional contraction  with S_pq^-1 (Q|R) (not on the derivatives)
               IF (.NOT. ri_data%same_op) THEN
                  CALL timeset(routineN//"_metric", handle)
                  DO k_mem = 1, n_mem_RI_fit
                     bounds_ctr_1d(1, 1) = batch_start_RI_fit(k_mem)
                     bounds_ctr_1d(2, 1) = batch_end_RI_fit(k_mem)

                     CALL dbt_batched_contract_init(t_2c_RI_inv)
                     CALL dbt_contract(1.0_dp, t_2c_RI_inv, t_3c_4, &
                                       1.0_dp, t_3c_6, &
                                       contract_1=[2], notcontract_1=[1], &
                                       contract_2=[1], notcontract_2=[2, 3], &
                                       bounds_1=bounds_ctr_1d, &
                                       bounds_3=bounds_ctr_2d, &
                                       map_1=[1], map_2=[2, 3], filter_eps=ri_data%filter_eps, &
                                       unit_nr=unit_nr_dbcsr, flop=nflop)
                     ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop
                     CALL dbt_batched_contract_finalize(t_2c_RI_inv)
                  END DO
                  CALL dbt_copy(t_3c_6, t_3c_4, move_data=.TRUE.)

                  ! 8) and get the force due to d/dx S^-1
                  DO k_mem = 1, n_mem_RI_fit
                     bounds_ctr_1d(1, 1) = batch_start_RI_fit(k_mem)
                     bounds_ctr_1d(2, 1) = batch_end_RI_fit(k_mem)

                     CALL dbt_batched_contract_init(t_2c_RI_met)
                     CALL dbt_contract(1.0_dp, t_3c_4, t_3c_5, &
                                       1.0_dp, t_2c_RI_met, &
                                       contract_1=[2, 3], notcontract_1=[1], &
                                       contract_2=[2, 3], notcontract_2=[1], &
                                       bounds_1=bounds_ctr_2d, &
                                       bounds_2=bounds_ctr_1d, &
                                       map_1=[1], map_2=[2], filter_eps=ri_data%filter_eps, &
                                       unit_nr=unit_nr_dbcsr, flop=nflop)
                     ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop
                     CALL dbt_batched_contract_finalize(t_2c_RI_met)
                  END DO
                  CALL timestop(handle)
               END IF
               CALL dbt_copy(t_3c_ri_mo_mo_fit, t_3c_5)

               ! 7) Do the force contribution due to 3c integrals (a'b|P) and (ab|P')

               ! (ab|P')
               CALL timeset(routineN//"_3c_RI", handle)
               pref = -0.5_dp*2.0_dp*hf_fraction*spin_fac
               CALL dbt_copy(t_3c_4, t_3c_RI_ctr, move_data=.TRUE.)
               IF (use_virial_prv) THEN
                  CALL get_force_from_3c_trace(force, t_3c_RI_ctr, t_3c_der_RI_ctr_2, atom_of_kind, kind_of, &
                                               idx_to_at_RI, pref, work_virial, cell, particle_set)
               ELSE
                  CALL get_force_from_3c_trace(force, t_3c_RI_ctr, t_3c_der_RI_ctr_2, atom_of_kind, kind_of, &
                                               idx_to_at_RI, pref)
               END IF
               CALL timestop(handle)

               ! (a'b|P) Note that derivative remains in AO rep until the actual force evaluation,
               ! which also prevents doing a direct 3-center trace
               bounds_ctr_2d(1, 1) = batch_start(i_mem)
               bounds_ctr_2d(2, 1) = batch_end(i_mem)

               bounds_ctr_1d(1, 1) = batch_start(j_mem)
               bounds_ctr_1d(2, 1) = batch_end(j_mem)

               CALL timeset(routineN//"_3c_AO", handle)
               CALL dbt_copy(t_3c_RI_ctr, t_3c_work, order=[2, 1, 3], move_data=.TRUE.)
               DO i_xyz = 1, 3

                  CALL dbt_batched_contract_init(t_2c_MO_AO_ctr(i_xyz))
                  DO k_mem = 1, n_mem_RI
                     bounds_ctr_2d(1, 2) = batch_start_RI(k_mem)
                     bounds_ctr_2d(2, 2) = batch_end_RI(k_mem)

                     CALL dbt_contract(1.0_dp, t_3c_work, t_3c_der_AO_ctr_1(i_xyz), &
                                       1.0_dp, t_2c_MO_AO_ctr(i_xyz), &
                                       contract_1=[1, 2], notcontract_1=[3], &
                                       contract_2=[1, 2], notcontract_2=[3], &
                                       map_1=[1], map_2=[2], filter_eps=ri_data%filter_eps, &
                                       bounds_1=bounds_ctr_2d, &
                                       bounds_2=bounds_ctr_1d, &
                                       unit_nr=unit_nr_dbcsr, flop=nflop)
                     ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop
                  END DO
                  CALL dbt_batched_contract_finalize(t_2c_MO_AO_ctr(i_xyz))
               END DO
               CALL timestop(handle)

            END DO !j_mem
            CALL dbt_batched_contract_finalize(t_3c_1)
            CALL dbt_batched_contract_finalize(t_3c_work)
            CALL dbt_batched_contract_finalize(t_3c_2)
            CALL dbt_batched_contract_finalize(t_3c_3)
            CALL dbt_batched_contract_finalize(t_3c_4)
            CALL dbt_batched_contract_finalize(t_3c_5)

            DO i_xyz = 1, 3
               CALL dbt_batched_contract_finalize(t_3c_der_RI_ctr_1(i_xyz))
               CALL dbt_batched_contract_finalize(t_3c_der_AO_ctr_1(i_xyz))
            END DO

            IF (.NOT. ri_data%same_op) THEN
               CALL dbt_batched_contract_finalize(t_3c_6)
            END IF

         END DO !i_mem
         CALL dbt_batched_contract_finalize(t_3c_desymm)
         CALL dbt_batched_contract_finalize(t_3c_0)

         DO i_xyz = 1, 3
            CALL dbt_batched_contract_finalize(t_3c_der_AO(i_xyz))
            CALL dbt_batched_contract_finalize(t_3c_der_RI(i_xyz))
         END DO

         !Force contribution due to 3-center AO derivatives (a'b|P)
         pref = -0.5_dp*4.0_dp*hf_fraction*spin_fac
         DO i_xyz = 1, 3
            CALL dbt_copy(t_2c_MO_AO_ctr(i_xyz), t_2c_MO_AO(i_xyz), move_data=.TRUE.) !ensures matching distributions
            IF (use_virial_prv) THEN
               CALL get_mo_ao_force(force, t_mo_cpy, t_2c_MO_AO(i_xyz), atom_of_kind, kind_of, idx_to_at_AO, pref, &
                                    i_xyz, work_virial, cell, particle_set)
            ELSE
               CALL get_mo_ao_force(force, t_mo_cpy, t_2c_MO_AO(i_xyz), atom_of_kind, kind_of, idx_to_at_AO, pref, i_xyz)
            END IF
            CALL dbt_clear(t_2c_MO_AO(i_xyz))
         END DO

         !Force contribution of d/dx (P|Q)
         pref = 0.5_dp*hf_fraction*spin_fac
         IF (.NOT. ri_data%same_op) pref = -pref

         !Making sure dists of the t_2c_RI tensors match
         CALL dbt_copy(t_2c_RI_PQ, t_2c_RI, move_data=.TRUE.)
         IF (use_virial_prv) THEN
            CALL get_2c_der_force(force, t_2c_RI, t_2c_der_RI, atom_of_kind, &
                                  kind_of, idx_to_at_RI, pref, work_virial, cell, particle_set)
         ELSE
            CALL get_2c_der_force(force, t_2c_RI, t_2c_der_RI, atom_of_kind, &
                                  kind_of, idx_to_at_RI, pref)

         END IF
         CALL dbt_clear(t_2c_RI)

         !Force contribution due to the inverse metric
         IF (.NOT. ri_data%same_op) THEN
            pref = 0.5_dp*2.0_dp*hf_fraction*spin_fac

            CALL dbt_copy(t_2c_RI_met, t_2c_RI, move_data=.TRUE.)
            IF (use_virial_prv) THEN
               CALL get_2c_der_force(force, t_2c_RI, t_2c_der_metric, atom_of_kind, &
                                     kind_of, idx_to_at_RI, pref, work_virial, cell, particle_set)
            ELSE
               CALL get_2c_der_force(force, t_2c_RI, t_2c_der_metric, atom_of_kind, &
                                     kind_of, idx_to_at_RI, pref)
            END IF
            CALL dbt_clear(t_2c_RI)
         END IF

         IF (use_virial_prv) THEN
            DO k_xyz = 1, 3
               DO j_xyz = 1, 3
                  DO i_xyz = 1, 3
                     virial%pv_fock_4c(i_xyz, j_xyz) = virial%pv_fock_4c(i_xyz, j_xyz) &
                                                       + work_virial(i_xyz, k_xyz)*cell%hmat(j_xyz, k_xyz)
                  END DO
               END DO
            END DO
         END IF

         CALL dbt_destroy(t_3c_0)
         CALL dbt_destroy(t_3c_1)
         CALL dbt_destroy(t_3c_2)
         CALL dbt_destroy(t_3c_3)
         CALL dbt_destroy(t_3c_4)
         CALL dbt_destroy(t_3c_5)
         CALL dbt_destroy(t_3c_6)
         CALL dbt_destroy(t_3c_work)
         CALL dbt_destroy(t_3c_RI_ctr)
         CALL dbt_destroy(t_3c_mo_ri_ao)
         CALL dbt_destroy(t_3c_mo_ri_mo)
         CALL dbt_destroy(t_3c_ri_mo_mo)
         CALL dbt_destroy(t_3c_ri_mo_mo_fit)
         CALL dbt_destroy(t_mo_coeff)
         CALL dbt_destroy(t_mo_cpy)
         DO i_xyz = 1, 3
            CALL dbt_destroy(t_2c_MO_AO(i_xyz))
            CALL dbt_destroy(t_2c_MO_AO_ctr(i_xyz))
            CALL dbt_destroy(t_3c_der_RI_ctr_1(i_xyz))
            CALL dbt_destroy(t_3c_der_AO_ctr_1(i_xyz))
            CALL dbt_destroy(t_3c_der_RI_ctr_2(i_xyz))
         END DO
         DEALLOCATE (batch_ranges, batch_start, batch_end)
      END DO !ispin

      ! Clean-up
      CALL dbt_pgrid_destroy(pgrid_1)
      CALL dbt_pgrid_destroy(pgrid_2)
      CALL dbt_destroy(t_3c_desymm)
      CALL dbt_destroy(t_2c_RI)
      CALL dbt_destroy(t_2c_RI_PQ)
      IF (.NOT. ri_data%same_op) THEN
         CALL dbt_destroy(t_2c_RI_met)
         CALL dbt_destroy(t_2c_RI_inv)
      END IF
      DO i_xyz = 1, 3
         CALL dbt_destroy(t_3c_der_AO(i_xyz))
         CALL dbt_destroy(t_3c_der_RI(i_xyz))
         CALL dbt_destroy(t_2c_der_RI(i_xyz))
         IF (.NOT. ri_data%same_op) CALL dbt_destroy(t_2c_der_metric(i_xyz))
      END DO
      CALL dbt_copy(ri_data%t_3c_int_ctr_2(1, 1), ri_data%t_3c_int_ctr_1(1, 1))

      CALL mp_sync(para_env%group)
      t2 = m_walltime()

      ri_data%dbcsr_time = ri_data%dbcsr_time + t2 - t1

   END SUBROUTINE hfx_ri_forces_mo

! **************************************************************************************************
!> \brief New sparser implementation
!> \param qs_env ...
!> \param ri_data ...
!> \param nspins ...
!> \param hf_fraction ...
!> \param rho_ao ...
!> \param rho_ao_resp ...
!> \param use_virial ...
!> \param resp_only ...
!> \param rescale_factor ...
! **************************************************************************************************
   SUBROUTINE hfx_ri_forces_Pmat(qs_env, ri_data, nspins, hf_fraction, rho_ao, rho_ao_resp, &
                                 use_virial, resp_only, rescale_factor)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      INTEGER, INTENT(IN)                                :: nspins
      REAL(dp), INTENT(IN)                               :: hf_fraction
      TYPE(dbcsr_p_type), DIMENSION(:, :)                :: rho_ao
      TYPE(dbcsr_p_type), DIMENSION(:), OPTIONAL         :: rho_ao_resp
      LOGICAL, INTENT(IN), OPTIONAL                      :: use_virial, resp_only
      REAL(dp), INTENT(IN), OPTIONAL                     :: rescale_factor

      CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_ri_forces_Pmat'

      INTEGER                                            :: handle, i_mem, i_spin, i_xyz, j_mem, &
                                                            j_xyz, k_xyz, n_mem, natom, &
                                                            unit_nr_dbcsr
      INTEGER(int_8)                                     :: nflop
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind, batch_end, batch_ranges, &
                                                            batch_start, dist1, dist2, dist3, &
                                                            idx_to_at_AO, idx_to_at_RI, kind_of
      INTEGER, DIMENSION(2, 1)                           :: ibounds, jbounds
      INTEGER, DIMENSION(2, 2)                           :: ijbounds
      INTEGER, DIMENSION(2, 3)                           :: bounds_cpy
      LOGICAL                                            :: do_resp, resp_only_prv, use_virial_prv
      REAL(dp)                                           :: pref, spin_fac, t1, t2
      REAL(dp), DIMENSION(3, 3)                          :: work_virial
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_type)                                   :: dbcsr_tmp
      TYPE(dbt_type) :: rho_ao_1, rho_ao_2, t_2c_RI, t_2c_RI_tmp, t_2c_tmp, t_3c_1, t_3c_2, &
         t_3c_3, t_3c_4, t_3c_ao_ri_ao, t_3c_help_1, t_3c_help_2, t_3c_int, t_3c_int_2, &
         t_3c_ri_ao_ao, t_3c_sparse, t_R, t_SVS
      TYPE(dbt_type), DIMENSION(3)                       :: t_2c_der_metric, t_2c_der_RI, &
                                                            t_3c_der_AO, t_3c_der_RI, t_3c_tmp
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(virial_type), POINTER                         :: virial

      !The idea is the following: we need to compute the gradients
      ! d/dx [P_ab P_cd (acP) S^-1_PQ (Q|R) S^-1_RS (Sbd)]
      ! Which we do in a few steps:
      ! 1) Contract the density matrices with the 3c integrals: M_acS = P_ab P_cd (Sbd)
      ! 2) Calculate the 3c contributions: d/dx (acP) [S^-1_PQ (Q|R) S^-1_RS M_acS]
      !    For maximum perf, we first multiply all 2c matrices together, than contract with retain_sparsity
      ! 3) Contract the 3c integrals and the M tensor together in order to only work with 2c quantities:
      !    R_PS = (acP) M_acS
      ! 4) From there, we can easily calculate the 2c contributions to the force:
      !    Potential: [S^-1*R*S^-1]_QR d/dx (Q|R)
      !    Metric:    [S^-1*R*S^-1*(Q|R)*S^-1]_UV d/dx S_UV

      NULLIFY (particle_set, virial, cell, force, atomic_kind_set)

      use_virial_prv = .FALSE.
      IF (PRESENT(use_virial)) use_virial_prv = use_virial

      do_resp = .FALSE.
      IF (PRESENT(rho_ao_resp)) THEN
         IF (ASSOCIATED(rho_ao_resp(1)%matrix)) do_resp = .TRUE.
      END IF

      resp_only_prv = .FALSE.
      IF (PRESENT(resp_only)) resp_only_prv = resp_only

      unit_nr_dbcsr = ri_data%unit_nr_dbcsr

      CALL get_qs_env(qs_env, natom=natom, particle_set=particle_set, &
                      atomic_kind_set=atomic_kind_set, virial=virial, &
                      cell=cell, force=force, para_env=para_env)

      CALL create_3c_tensor(t_3c_ao_ri_ao, dist1, dist2, dist3, ri_data%pgrid_1, &
                            ri_data%bsizes_AO_split, ri_data%bsizes_RI_split, ri_data%bsizes_AO_split, &
                            [1, 2], [3], name="(AO RI | AO)")
      DEALLOCATE (dist1, dist2, dist3)

      CALL create_3c_tensor(t_3c_ri_ao_ao, dist1, dist2, dist3, ri_data%pgrid_2, &
                            ri_data%bsizes_RI_split, ri_data%bsizes_AO_split, ri_data%bsizes_AO_split, &
                            [1], [2, 3], name="(RI | AO AO)")
      DEALLOCATE (dist1, dist2, dist3)

      ! Precompute the derivatives
      CALL precalc_derivatives(t_3c_der_RI, t_3c_tmp, t_2c_der_RI, t_2c_der_metric, &
                               t_3c_ri_ao_ao, t_3c_ao_ri_ao, ri_data, qs_env)

      CALL dbt_create(t_3c_ri_ao_ao, t_3c_sparse)
      DO i_xyz = 1, 3
         CALL dbt_create(t_3c_ri_ao_ao, t_3c_der_AO(i_xyz))
         CALL dbt_copy(t_3c_tmp(i_xyz), t_3c_der_AO(i_xyz), order=[2, 1, 3], move_data=.TRUE.)
         CALL dbt_destroy(t_3c_tmp(i_xyz))

         CALL dbt_copy(t_3c_der_RI(i_xyz), t_3c_sparse, summation=.TRUE.)
         CALL dbt_copy(t_3c_der_AO(i_xyz), t_3c_sparse, summation=.TRUE.)
         CALL dbt_copy(t_3c_der_AO(i_xyz), t_3c_sparse, order=[1, 3, 2], summation=.TRUE.)
      END DO

      ! Some utilities
      spin_fac = 0.5_dp
      IF (nspins == 2) spin_fac = 1.0_dp
      IF (PRESENT(rescale_factor)) spin_fac = spin_fac*rescale_factor

      ALLOCATE (idx_to_at_RI(SIZE(ri_data%bsizes_RI_split)))
      CALL get_idx_to_atom(idx_to_at_RI, ri_data%bsizes_RI_split, ri_data%bsizes_RI)

      ALLOCATE (idx_to_at_AO(SIZE(ri_data%bsizes_AO_split)))
      CALL get_idx_to_atom(idx_to_at_AO, ri_data%bsizes_AO_split, ri_data%bsizes_AO)

      ALLOCATE (atom_of_kind(natom), kind_of(natom))
      CALL get_atomic_kind_set(atomic_kind_set, kind_of=kind_of, atom_of_kind=atom_of_kind)

      ! Go over batches of the 2 AO indices to save memory
      n_mem = ri_data%n_mem
      ALLOCATE (batch_start(n_mem), batch_end(n_mem))
      batch_start(:) = ri_data%starts_array_mem(:)
      batch_end(:) = ri_data%ends_array_mem(:)

      ALLOCATE (batch_ranges(n_mem + 1))
      batch_ranges(:n_mem) = ri_data%starts_array_mem_block(:)
      batch_ranges(n_mem + 1) = ri_data%ends_array_mem_block(n_mem) + 1

      ! Pre-create all the needed tensors
      CALL create_2c_tensor(rho_ao_1, dist1, dist2, ri_data%pgrid_2d, &
                            ri_data%bsizes_AO_split, ri_data%bsizes_AO_split, &
                            name="(AO | AO)")
      DEALLOCATE (dist1, dist2)
      CALL dbt_create(rho_ao_1, rho_ao_2)

      CALL create_2c_tensor(t_2c_RI, dist1, dist2, ri_data%pgrid_2d, &
                            ri_data%bsizes_RI_split, ri_data%bsizes_RI_split, name="(RI | RI)")
      DEALLOCATE (dist1, dist2)
      CALL dbt_create(t_2c_RI, t_SVS)
      CALL dbt_create(t_2c_RI, t_R)
      CALL dbt_create(t_2c_RI, t_2c_RI_tmp)

      CALL dbt_create(t_3c_ao_ri_ao, t_3c_1)
      CALL dbt_create(t_3c_ao_ri_ao, t_3c_2)
      CALL dbt_create(t_3c_ri_ao_ao, t_3c_3)
      CALL dbt_create(t_3c_ri_ao_ao, t_3c_4)
      CALL dbt_create(t_3c_ri_ao_ao, t_3c_help_1)
      CALL dbt_create(t_3c_ri_ao_ao, t_3c_help_2)

      CALL dbt_create(t_3c_ao_ri_ao, t_3c_int)
      CALL dbt_copy(ri_data%t_3c_int_ctr_2(1, 1), t_3c_int)

      CALL dbt_create(t_3c_ri_ao_ao, t_3c_int_2)
      CALL dbt_copy(t_3c_int, t_3c_int_2, order=[2, 1, 3])

      CALL mp_sync(para_env%group)
      t1 = m_walltime()

      !Pre-calculate the necessary 2-center quantities
      IF (.NOT. ri_data%same_op) THEN
         !S^-1 * V * S^-1
         CALL dbt_contract(1.0_dp, ri_data%t_2c_inv(1, 1), ri_data%t_2c_pot(1, 1), 0.0_dp, t_2c_RI, &
                           contract_1=[2], notcontract_1=[1], &
                           contract_2=[1], notcontract_2=[2], &
                           map_1=[1], map_2=[2], filter_eps=ri_data%filter_eps, &
                           unit_nr=unit_nr_dbcsr, flop=nflop)
         ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop

         CALL dbt_contract(1.0_dp, t_2c_RI, ri_data%t_2c_inv(1, 1), 0.0_dp, t_SVS, &
                           contract_1=[2], notcontract_1=[1], &
                           contract_2=[1], notcontract_2=[2], &
                           map_1=[1], map_2=[2], filter_eps=ri_data%filter_eps, &
                           unit_nr=unit_nr_dbcsr, flop=nflop)
         ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop
      ELSE
         ! Simply V^-1
         CALL dbt_copy(ri_data%t_2c_inv(1, 1), t_SVS)
      END IF

      CALL dbt_batched_contract_init(t_3c_int, batch_range_1=batch_ranges, batch_range_3=batch_ranges)
      CALL dbt_batched_contract_init(t_3c_int_2, batch_range_2=batch_ranges, batch_range_3=batch_ranges)
      CALL dbt_batched_contract_init(t_3c_1, batch_range_1=batch_ranges, batch_range_3=batch_ranges)
      CALL dbt_batched_contract_init(t_3c_2, batch_range_1=batch_ranges, batch_range_3=batch_ranges)
      CALL dbt_batched_contract_init(t_3c_3, batch_range_2=batch_ranges, batch_range_3=batch_ranges)
      CALL dbt_batched_contract_init(t_3c_4, batch_range_2=batch_ranges, batch_range_3=batch_ranges)

      DO i_spin = 1, nspins

         !Prepare Pmat in tensor format
         CALL dbt_create(rho_ao(i_spin, 1)%matrix, t_2c_tmp)
         CALL dbt_copy_matrix_to_tensor(rho_ao(i_spin, 1)%matrix, t_2c_tmp)
         CALL dbt_copy(t_2c_tmp, rho_ao_1, move_data=.TRUE.)
         CALL dbt_destroy(t_2c_tmp)

         IF (.NOT. do_resp) THEN
            CALL dbt_copy(rho_ao_1, rho_ao_2)
         ELSE IF (do_resp .AND. resp_only_prv) THEN

            CALL dbt_create(rho_ao_resp(i_spin)%matrix, t_2c_tmp)
            CALL dbt_copy_matrix_to_tensor(rho_ao_resp(i_spin)%matrix, t_2c_tmp)
            CALL dbt_copy(t_2c_tmp, rho_ao_2)
            !symmetry allows to take 2*P_resp rasther than explicitely take all cross products
            CALL dbt_copy(t_2c_tmp, rho_ao_2, summation=.TRUE., move_data=.TRUE.)
            CALL dbt_destroy(t_2c_tmp)
         ELSE

            !if not resp_only, need P-P_resp and P+P_resp
            CALL dbt_copy(rho_ao_1, rho_ao_2)
            CALL dbcsr_create(dbcsr_tmp, template=rho_ao_resp(i_spin)%matrix)
            CALL dbcsr_add(dbcsr_tmp, rho_ao_resp(i_spin)%matrix, 0.0_dp, -1.0_dp)
            CALL dbt_create(dbcsr_tmp, t_2c_tmp)
            CALL dbt_copy_matrix_to_tensor(dbcsr_tmp, t_2c_tmp)
            CALL dbt_copy(t_2c_tmp, rho_ao_1, summation=.TRUE., move_data=.TRUE.)
            CALL dbcsr_release(dbcsr_tmp)

            CALL dbt_copy_matrix_to_tensor(rho_ao_resp(i_spin)%matrix, t_2c_tmp)
            CALL dbt_copy(t_2c_tmp, rho_ao_2, summation=.TRUE., move_data=.TRUE.)
            CALL dbt_destroy(t_2c_tmp)

         END IF
         work_virial = 0.0_dp

         CALL timeset(routineN//"_3c", handle)
         !Start looping of the batches
         DO i_mem = 1, n_mem
            ibounds(:, 1) = [batch_start(i_mem), batch_end(i_mem)]

            CALL dbt_batched_contract_init(rho_ao_1)
            CALL dbt_contract(1.0_dp, rho_ao_1, t_3c_int, 0.0_dp, t_3c_1, &
                              contract_1=[1], notcontract_1=[2], &
                              contract_2=[3], notcontract_2=[1, 2], &
                              map_1=[3], map_2=[1, 2], filter_eps=ri_data%filter_eps, &
                              bounds_2=ibounds, unit_nr=unit_nr_dbcsr, flop=nflop)
            ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop
            CALL dbt_batched_contract_finalize(rho_ao_1)

            CALL dbt_copy(t_3c_1, t_3c_2, order=[3, 2, 1], move_data=.TRUE.)

            DO j_mem = 1, n_mem
               jbounds(:, 1) = [batch_start(j_mem), batch_end(j_mem)]

               CALL dbt_batched_contract_init(rho_ao_2)
               CALL dbt_contract(1.0_dp, rho_ao_2, t_3c_2, 0.0_dp, t_3c_1, &
                                 contract_1=[1], notcontract_1=[2], &
                                 contract_2=[3], notcontract_2=[1, 2], &
                                 map_1=[3], map_2=[1, 2], filter_eps=ri_data%filter_eps, &
                                 bounds_2=jbounds, unit_nr=unit_nr_dbcsr, flop=nflop)
               ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop
               CALL dbt_batched_contract_finalize(rho_ao_2)

               bounds_cpy(:, 1) = [1, SUM(ri_data%bsizes_RI)]
               bounds_cpy(:, 2) = [batch_start(i_mem), batch_end(i_mem)]
               bounds_cpy(:, 3) = [batch_start(j_mem), batch_end(j_mem)]
               CALL dbt_copy(t_3c_1, t_3c_3, order=[2, 1, 3], move_data=.TRUE.)
               CALL dbt_copy(t_3c_sparse, t_3c_4, bounds=bounds_cpy)

               !Contract with the 2-center product S^-1 * V * S^-1 while keeping sparsity of derivatives
               CALL dbt_batched_contract_init(t_SVS)
               CALL dbt_contract(1.0_dp, t_SVS, t_3c_3, 0.0_dp, t_3c_4, &
                                 contract_1=[2], notcontract_1=[1], &
                                 contract_2=[1], notcontract_2=[2, 3], &
                                 map_1=[1], map_2=[2, 3], filter_eps=ri_data%filter_eps, &
                                 retain_sparsity=.TRUE., unit_nr=unit_nr_dbcsr, flop=nflop)
               ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop
               CALL dbt_batched_contract_finalize(t_SVS)

               CALL dbt_copy(t_3c_4, t_3c_help_1, move_data=.TRUE., summation=.TRUE.)

               ijbounds(:, 1) = ibounds(:, 1)
               ijbounds(:, 2) = jbounds(:, 1)

               !Contract R_PS = (acP) M_acS
               CALL dbt_batched_contract_init(t_R)
               CALL dbt_contract(1.0_dp, t_3c_int_2, t_3c_3, 1.0_dp, t_R, &
                                 contract_1=[2, 3], notcontract_1=[1], &
                                 contract_2=[2, 3], notcontract_2=[1], &
                                 map_1=[1], map_2=[2], filter_eps=ri_data%filter_eps, &
                                 bounds_1=ijbounds, unit_nr=unit_nr_dbcsr, flop=nflop)
               ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop
               CALL dbt_batched_contract_finalize(t_R)

               CALL dbt_copy(t_3c_4, t_3c_3)
            END DO !j_mem
         END DO !i_mem

         !The force from the 3c derivatives
         pref = -0.5_dp*2.0_dp*hf_fraction*spin_fac
         IF (use_virial_prv) THEN
            CALL get_force_from_3c_trace(force, t_3c_help_1, t_3c_der_RI, atom_of_kind, kind_of, &
                                         idx_to_at_RI, pref, work_virial, cell, particle_set)
         ELSE
            CALL get_force_from_3c_trace(force, t_3c_help_1, t_3c_der_RI, atom_of_kind, kind_of, &
                                         idx_to_at_RI, pref)
         END IF

         pref = -0.5_dp*4.0_dp*hf_fraction*spin_fac
         IF (do_resp) pref = 0.5_dp*pref
         IF (use_virial_prv) THEN
            CALL get_force_from_3c_trace(force, t_3c_help_1, t_3c_der_AO, atom_of_kind, kind_of, &
                                         idx_to_at_AO, pref, work_virial, cell, particle_set, deriv_dim=2)
         ELSE
            CALL get_force_from_3c_trace(force, t_3c_help_1, t_3c_der_AO, atom_of_kind, kind_of, &
                                         idx_to_at_AO, pref, deriv_dim=2)
         END IF

         IF (do_resp) THEN
            CALL dbt_copy(t_3c_help_1, t_3c_help_2, order=[1, 3, 2], move_data=.TRUE.)
            IF (use_virial_prv) THEN
               CALL get_force_from_3c_trace(force, t_3c_help_2, t_3c_der_AO, atom_of_kind, kind_of, &
                                            idx_to_at_AO, pref, work_virial, cell, particle_set, deriv_dim=2)
            ELSE
               CALL get_force_from_3c_trace(force, t_3c_help_2, t_3c_der_AO, atom_of_kind, kind_of, &
                                            idx_to_at_AO, pref, deriv_dim=2)
            END IF
         END IF
         CALL timestop(handle)

         CALL timeset(routineN//"_2c", handle)
         !Now we deal with all the 2-center quantities
         !Calculate S^-1 * R * S^-1
         CALL dbt_contract(1.0_dp, ri_data%t_2c_inv(1, 1), t_R, 0.0_dp, t_2c_RI_tmp, &
                           contract_1=[2], notcontract_1=[1], &
                           contract_2=[1], notcontract_2=[2], &
                           map_1=[1], map_2=[2], filter_eps=ri_data%filter_eps, &
                           unit_nr=unit_nr_dbcsr, flop=nflop)
         ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop

         CALL dbt_contract(1.0_dp, t_2c_RI_tmp, ri_data%t_2c_inv(1, 1), 0.0_dp, t_2c_RI, &
                           contract_1=[2], notcontract_1=[1], &
                           contract_2=[1], notcontract_2=[2], &
                           map_1=[1], map_2=[2], filter_eps=ri_data%filter_eps, &
                           unit_nr=unit_nr_dbcsr, flop=nflop)
         ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop

         !Calculate the potential contribution to the force: [S^-1*R*S^-1]_QR d/dx (Q|R)
         pref = 0.5_dp*hf_fraction*spin_fac
         IF (.NOT. ri_data%same_op) pref = -pref
         IF (use_virial_prv) THEN
            CALL get_2c_der_force(force, t_2c_RI, t_2c_der_RI, atom_of_kind, &
                                  kind_of, idx_to_at_RI, pref, work_virial, cell, particle_set)
         ELSE
            CALL get_2c_der_force(force, t_2c_RI, t_2c_der_RI, atom_of_kind, &
                                  kind_of, idx_to_at_RI, pref)

         END IF

         !And that from the metric: [S^-1*R*S^-1*(Q|R)*S^-1]_UV d/dx S_UV
         IF (.NOT. ri_data%same_op) THEN
            CALL dbt_contract(1.0_dp, t_2c_RI, ri_data%t_2c_pot(1, 1), 0.0_dp, t_2c_RI_tmp, &
                              contract_1=[2], notcontract_1=[1], &
                              contract_2=[1], notcontract_2=[2], &
                              map_1=[1], map_2=[2], filter_eps=ri_data%filter_eps, &
                              unit_nr=unit_nr_dbcsr, flop=nflop)
            ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop

            CALL dbt_contract(1.0_dp, t_2c_RI_tmp, ri_data%t_2c_inv(1, 1), 0.0_dp, t_2c_RI, &
                              contract_1=[2], notcontract_1=[1], &
                              contract_2=[1], notcontract_2=[2], &
                              map_1=[1], map_2=[2], filter_eps=ri_data%filter_eps, &
                              unit_nr=unit_nr_dbcsr, flop=nflop)
            ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop

            pref = 0.5_dp*2.0_dp*hf_fraction*spin_fac
            IF (use_virial_prv) THEN
               CALL get_2c_der_force(force, t_2c_RI, t_2c_der_metric, atom_of_kind, &
                                     kind_of, idx_to_at_RI, pref, work_virial, cell, particle_set)
            ELSE
               CALL get_2c_der_force(force, t_2c_RI, t_2c_der_metric, atom_of_kind, &
                                     kind_of, idx_to_at_RI, pref)
            END IF
         END IF
         CALL dbt_clear(t_2c_RI)
         CALL dbt_clear(t_2c_RI_tmp)
         CALL dbt_clear(t_R)
         CALL dbt_clear(t_3c_help_1)
         CALL dbt_clear(t_3c_help_2)
         CALL timestop(handle)

         IF (use_virial_prv) THEN
            DO k_xyz = 1, 3
               DO j_xyz = 1, 3
                  DO i_xyz = 1, 3
                     virial%pv_fock_4c(i_xyz, j_xyz) = virial%pv_fock_4c(i_xyz, j_xyz) &
                                                       + work_virial(i_xyz, k_xyz)*cell%hmat(j_xyz, k_xyz)
                  END DO
               END DO
            END DO
         END IF
      END DO !i_spin

      CALL dbt_batched_contract_finalize(t_3c_int)
      CALL dbt_batched_contract_finalize(t_3c_int_2)
      CALL dbt_batched_contract_finalize(t_3c_1)
      CALL dbt_batched_contract_finalize(t_3c_2)
      CALL dbt_batched_contract_finalize(t_3c_3)
      CALL dbt_batched_contract_finalize(t_3c_4)

      CALL mp_sync(para_env%group)
      t2 = m_walltime()

      CALL dbt_copy(t_3c_int, ri_data%t_3c_int_ctr_2(1, 1), move_data=.TRUE.)

      !clean-up
      CALL dbt_destroy(rho_ao_1)
      CALL dbt_destroy(rho_ao_2)
      CALL dbt_destroy(t_3c_ao_ri_ao)
      CALL dbt_destroy(t_3c_ri_ao_ao)
      CALL dbt_destroy(t_3c_int)
      CALL dbt_destroy(t_3c_int_2)
      CALL dbt_destroy(t_3c_1)
      CALL dbt_destroy(t_3c_2)
      CALL dbt_destroy(t_3c_3)
      CALL dbt_destroy(t_3c_4)
      CALL dbt_destroy(t_3c_help_1)
      CALL dbt_destroy(t_3c_help_2)
      CALL dbt_destroy(t_3c_sparse)
      CALL dbt_destroy(t_SVS)
      CALL dbt_destroy(t_R)
      CALL dbt_destroy(t_2c_RI)
      CALL dbt_destroy(t_2c_RI_tmp)

      DO i_xyz = 1, 3
         CALL dbt_destroy(t_3c_der_RI(i_xyz))
         CALL dbt_destroy(t_3c_der_AO(i_xyz))
         CALL dbt_destroy(t_2c_der_RI(i_xyz))
         IF (.NOT. ri_data%same_op) CALL dbt_destroy(t_2c_der_metric(i_xyz))
      END DO

   END SUBROUTINE hfx_ri_forces_Pmat

! **************************************************************************************************
!> \brief the general routine that calls the relevant force code
!> \param qs_env ...
!> \param ri_data ...
!> \param nspins ...
!> \param hf_fraction ...
!> \param rho_ao ...
!> \param rho_ao_resp ...
!> \param mos ...
!> \param use_virial ...
!> \param resp_only ...
!> \param rescale_factor ...
! **************************************************************************************************
   SUBROUTINE hfx_ri_update_forces(qs_env, ri_data, nspins, hf_fraction, rho_ao, rho_ao_resp, &
                                   mos, use_virial, resp_only, rescale_factor)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      INTEGER, INTENT(IN)                                :: nspins
      REAL(KIND=dp), INTENT(IN)                          :: hf_fraction
      TYPE(dbcsr_p_type), DIMENSION(:, :), OPTIONAL      :: rho_ao
      TYPE(dbcsr_p_type), DIMENSION(:), OPTIONAL         :: rho_ao_resp
      TYPE(mo_set_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: mos
      LOGICAL, INTENT(IN), OPTIONAL                      :: use_virial, resp_only
      REAL(dp), INTENT(IN), OPTIONAL                     :: rescale_factor

      CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_ri_update_forces'

      INTEGER                                            :: handle, ispin
      INTEGER, DIMENSION(2)                              :: homo
      REAL(KIND=dp), DIMENSION(:), POINTER               :: mo_eigenvalues
      TYPE(cp_1d_r_p_type), DIMENSION(:), POINTER        :: occupied_evals
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: homo_localized, moloc_coeff, &
                                                            occupied_orbs
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(dbcsr_type), DIMENSION(2)                     :: mo_coeff_b
      TYPE(dbcsr_type), POINTER                          :: mo_coeff_b_tmp
      TYPE(mo_set_type), POINTER                         :: mo_set

      CALL timeset(routineN, handle)

      SELECT CASE (ri_data%flavor)
      CASE (ri_mo)

         IF (ri_data%do_loc) THEN
            ALLOCATE (occupied_orbs(nspins))
            ALLOCATE (occupied_evals(nspins))
            ALLOCATE (homo_localized(nspins))
         END IF
         DO ispin = 1, nspins
            NULLIFY (mo_coeff_b_tmp)
            mo_set => mos(ispin)%mo_set
            CPASSERT(mo_set%uniform_occupation)
            CALL get_mo_set(mo_set=mo_set, mo_coeff=mo_coeff, eigenvalues=mo_eigenvalues, mo_coeff_b=mo_coeff_b_tmp)

            IF (.NOT. ri_data%do_loc) THEN
               IF (.NOT. mo_set%use_mo_coeff_b) CALL copy_fm_to_dbcsr(mo_coeff, mo_coeff_b_tmp)
               CALL dbcsr_copy(mo_coeff_b(ispin), mo_coeff_b_tmp)
            ELSE
               IF (mo_set%use_mo_coeff_b) CALL copy_dbcsr_to_fm(mo_coeff_b_tmp, mo_coeff)
               CALL dbcsr_create(mo_coeff_b(ispin), template=mo_coeff_b_tmp)
            END IF

            IF (ri_data%do_loc) THEN
               occupied_orbs(ispin)%matrix => mo_coeff
               occupied_evals(ispin)%array => mo_eigenvalues
               CALL cp_fm_create(homo_localized(ispin)%matrix, occupied_orbs(ispin)%matrix%matrix_struct)
               CALL cp_fm_to_fm(occupied_orbs(ispin)%matrix, homo_localized(ispin)%matrix)
            END IF
         END DO

         IF (ri_data%do_loc) THEN
            CALL qs_loc_env_create(ri_data%qs_loc_env)
            CALL qs_loc_control_init(ri_data%qs_loc_env, ri_data%loc_subsection, do_homo=.TRUE.)
            CALL qs_loc_init(qs_env, ri_data%qs_loc_env, ri_data%loc_subsection, homo_localized)
            DO ispin = 1, nspins
               CALL qs_loc_driver(qs_env, ri_data%qs_loc_env, ri_data%print_loc_subsection, ispin, &
                                  ext_mo_coeff=homo_localized(ispin)%matrix)
            END DO
            CALL get_qs_loc_env(qs_loc_env=ri_data%qs_loc_env, moloc_coeff=moloc_coeff)

            DO ispin = 1, nspins
               CALL cp_fm_release(homo_localized(ispin)%matrix)
            END DO

            DEALLOCATE (occupied_orbs, occupied_evals, homo_localized)

         END IF

         DO ispin = 1, nspins
            mo_set => mos(ispin)%mo_set
            IF (ri_data%do_loc) THEN
               CALL copy_fm_to_dbcsr(moloc_coeff(ispin)%matrix, mo_coeff_b(ispin))
            END IF
            CALL dbcsr_scale(mo_coeff_b(ispin), SQRT(mo_set%maxocc))
            homo(ispin) = mo_set%homo
         END DO

         IF (ri_data%do_loc) CALL qs_loc_env_release(ri_data%qs_loc_env)

         CALL hfx_ri_forces_mo(qs_env, ri_data, nspins, hf_fraction, mo_coeff_b, use_virial)

      CASE (ri_pmat)

         CALL hfx_ri_forces_Pmat(qs_env, ri_data, nspins, hf_fraction, rho_ao, rho_ao_resp, use_virial, &
                                 resp_only, rescale_factor)
      END SELECT

      DO ispin = 1, nspins
         CALL dbcsr_release(mo_coeff_b(ispin))
      END DO

      CALL timestop(handle)

   END SUBROUTINE hfx_ri_update_forces

! **************************************************************************************************
!> \brief Calculate the derivatives tensors for the force, in a format fit for contractions
!> \param t_3c_der_RI format based on template
!> \param t_3c_der_AO format based on template
!> \param t_2c_der_RI format based on standard atomic block sizes
!> \param t_2c_der_metric format based on standard atomic block sizes
!> \param ri_ao_ao_template ...
!> \param ao_ri_ao_template ...
!> \param ri_data ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE precalc_derivatives(t_3c_der_RI, t_3c_der_AO, t_2c_der_RI, t_2c_der_metric, &
                                  ri_ao_ao_template, ao_ri_ao_template, ri_data, qs_env)

      TYPE(dbt_type), DIMENSION(3), INTENT(OUT)          :: t_3c_der_RI, t_3c_der_AO, t_2c_der_RI, &
                                                            t_2c_der_metric
      TYPE(dbt_type), INTENT(INOUT)                      :: ri_ao_ao_template, ao_ri_ao_template
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'precalc_derivatives'

      INTEGER                                            :: handle, i_mem, i_xyz, ibasis, &
                                                            mp_comm_t3c, n_mem, nkind
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: dist1, dist2, dist_AO_1, dist_AO_2, &
                                                            dist_RI, dummy_end, dummy_start, &
                                                            end_blocks, start_blocks
      INTEGER, DIMENSION(3)                              :: pcoord, pdims
      INTEGER, DIMENSION(:), POINTER                     :: col_bsize, row_bsize
      TYPE(dbcsr_distribution_type)                      :: dbcsr_dist
      TYPE(dbcsr_type), DIMENSION(1, 3)                  :: t_2c_der_metric_prv, t_2c_der_RI_prv
      TYPE(dbt_type)                                     :: t_2c_template, t_2c_tmp, t_3c_template
      TYPE(dbt_type), DIMENSION(1, 1, 3)                 :: t_3c_der_AO_prv, t_3c_der_RI_prv
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(distribution_2d_type), POINTER                :: dist_2d
      TYPE(distribution_3d_type)                         :: dist_3d
      TYPE(gto_basis_set_p_type), ALLOCATABLE, &
         DIMENSION(:), TARGET                            :: basis_set_AO, basis_set_RI
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis, ri_basis
      TYPE(neighbor_list_3c_type)                        :: nl_3c
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: nl_2c
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (qs_kind_set, orb_basis, dist_2d, nl_2c, particle_set, dft_control)

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, nkind=nkind, qs_kind_set=qs_kind_set, distribution_2d=dist_2d, &
                      particle_set=particle_set, dft_control=dft_control)

      ALLOCATE (basis_set_RI(nkind), basis_set_AO(nkind))
      CALL basis_set_list_setup(basis_set_RI, ri_data%ri_basis_type, qs_kind_set)
      CALL get_particle_set(particle_set, qs_kind_set, basis=basis_set_RI)
      CALL basis_set_list_setup(basis_set_AO, ri_data%orb_basis_type, qs_kind_set)
      CALL get_particle_set(particle_set, qs_kind_set, basis=basis_set_AO)

      DO ibasis = 1, SIZE(basis_set_AO)
         orb_basis => basis_set_AO(ibasis)%gto_basis_set
         CALL init_interaction_radii_orb_basis(orb_basis, ri_data%eps_pgf_orb)
         ri_basis => basis_set_RI(ibasis)%gto_basis_set
         CALL init_interaction_radii_orb_basis(ri_basis, ri_data%eps_pgf_orb)
      END DO

      !Dealing with the 3c derivatives
      CALL create_3c_tensor(t_3c_template, dist_RI, dist_AO_1, dist_AO_2, ri_data%pgrid, &
                            ri_data%bsizes_RI, ri_data%bsizes_AO, ri_data%bsizes_AO, &
                            map1=[1], map2=[2, 3], &
                            name="der (RI AO | AO)")

      DO i_xyz = 1, 3
         CALL dbt_create(t_3c_template, t_3c_der_RI_prv(1, 1, i_xyz))
         CALL dbt_create(t_3c_template, t_3c_der_AO_prv(1, 1, i_xyz))
      END DO
      CALL dbt_destroy(t_3c_template)

      CALL dbt_mp_environ_pgrid(ri_data%pgrid, pdims, pcoord)
      CALL mp_cart_create(ri_data%pgrid%mp_comm_2d, 3, pdims, pcoord, mp_comm_t3c)
      CALL distribution_3d_create(dist_3d, dist_RI, dist_AO_1, dist_AO_2, &
                                  nkind, particle_set, mp_comm_t3c, own_comm=.TRUE.)
      DEALLOCATE (dist_RI, dist_AO_1, dist_AO_2)

      CALL build_3c_neighbor_lists(nl_3c, basis_set_RI, basis_set_AO, basis_set_AO, dist_3d, ri_data%ri_metric, &
                                   "HFX_3c_nl", qs_env, op_pos=1, sym_jk=.TRUE., own_dist=.TRUE.)

      !Output tensor must be in a format fit for contraction, with splitted blocks
      DO i_xyz = 1, 3
         CALL dbt_create(ri_ao_ao_template, t_3c_der_RI(i_xyz)) ! (RI | AO AO) format
         CALL dbt_create(ao_ri_ao_template, t_3c_der_AO(i_xyz)) !(AO RI | AO) format
      END DO

      n_mem = ri_data%n_mem
      CALL create_tensor_batches(ri_data%bsizes_AO, n_mem, dummy_start, dummy_end, &
                                 start_blocks, end_blocks)
      DEALLOCATE (dummy_start, dummy_end)

      DO i_mem = 1, n_mem
         CALL build_3c_derivatives(t_3c_der_RI_prv, t_3c_der_AO_prv, ri_data%filter_eps, qs_env, &
                                   nl_3c, basis_set_RI, basis_set_AO, basis_set_AO, &
                                   ri_data%ri_metric, der_eps=ri_data%eps_schwarz_forces, op_pos=1, &
                                   bounds_j=[start_blocks(i_mem), end_blocks(i_mem)])

         DO i_xyz = 1, 3
            CALL dbt_copy(t_3c_der_RI_prv(1, 1, i_xyz), t_3c_der_RI(i_xyz), &
                          move_data=.TRUE., summation=.TRUE.)
            CALL dbt_filter(t_3c_der_RI(i_xyz), ri_data%filter_eps)

            !put AO derivative as first index (intermediat step needed)
            CALL dbt_copy(t_3c_der_AO_prv(1, 1, i_xyz), ao_ri_ao_template, order=[2, 1, 3], move_data=.TRUE.)
            CALL dbt_copy(ao_ri_ao_template, t_3c_der_AO(i_xyz), order=[3, 2, 1], move_data=.TRUE., summation=.TRUE.)
            CALL dbt_filter(t_3c_der_AO(i_xyz), ri_data%filter_eps)
         END DO
      END DO
      CALL dbt_clear(ao_ri_ao_template)

      CALL neighbor_list_3c_destroy(nl_3c)

      DO i_xyz = 1, 3
         CALL dbt_destroy(t_3c_der_RI_prv(1, 1, i_xyz))
         CALL dbt_destroy(t_3c_der_AO_prv(1, 1, i_xyz))
      END DO

      !Deal with the 2-center derivatives
      CALL cp_dbcsr_dist2d_to_dist(dist_2d, dbcsr_dist)
      ALLOCATE (row_bsize(SIZE(ri_data%bsizes_RI)))
      ALLOCATE (col_bsize(SIZE(ri_data%bsizes_RI)))
      row_bsize(:) = ri_data%bsizes_RI
      col_bsize(:) = ri_data%bsizes_RI

      CALL build_2c_neighbor_lists(nl_2c, basis_set_RI, basis_set_RI, ri_data%hfx_pot, &
                                   "HFX_2c_nl_pot", qs_env, sym_ij=.TRUE., dist_2d=dist_2d)

      DO i_xyz = 1, 3
         CALL dbcsr_create(t_2c_der_RI_prv(1, i_xyz), "(R|P) HFX der", dbcsr_dist, &
                           dbcsr_type_antisymmetric, row_bsize, col_bsize)
      END DO

      CALL build_2c_derivatives(t_2c_der_RI_prv, ri_data%filter_eps_2c, qs_env, nl_2c, basis_set_RI, &
                                basis_set_RI, ri_data%hfx_pot)
      CALL release_neighbor_list_sets(nl_2c)

      !copy 2c derivative tensor into the standard format
      CALL create_2c_tensor(t_2c_template, dist1, dist2, ri_data%pgrid_2d, ri_data%bsizes_RI_split, &
                            ri_data%bsizes_RI_split, name='(RI| RI)')
      DEALLOCATE (dist1, dist2)

      DO i_xyz = 1, 3
         CALL dbt_create(t_2c_der_RI_prv(1, i_xyz), t_2c_tmp)
         CALL dbt_copy_matrix_to_tensor(t_2c_der_RI_prv(1, i_xyz), t_2c_tmp)

         CALL dbt_create(t_2c_template, t_2c_der_RI(i_xyz))
         CALL dbt_copy(t_2c_tmp, t_2c_der_RI(i_xyz), move_data=.TRUE.)

         CALL dbt_destroy(t_2c_tmp)
         CALL dbcsr_release(t_2c_der_RI_prv(1, i_xyz))
      END DO

      !Repeat with the metric, if required
      IF (.NOT. ri_data%same_op) THEN

         CALL build_2c_neighbor_lists(nl_2c, basis_set_RI, basis_set_RI, ri_data%ri_metric, &
                                      "HFX_2c_nl_RI", qs_env, sym_ij=.TRUE., dist_2d=dist_2d)

         DO i_xyz = 1, 3
            CALL dbcsr_create(t_2c_der_metric_prv(1, i_xyz), "(R|P) HFX der", dbcsr_dist, &
                              dbcsr_type_antisymmetric, row_bsize, col_bsize)
         END DO

         CALL build_2c_derivatives(t_2c_der_metric_prv, ri_data%filter_eps_2c, qs_env, nl_2c, &
                                   basis_set_RI, basis_set_RI, ri_data%ri_metric)
         CALL release_neighbor_list_sets(nl_2c)

         DO i_xyz = 1, 3
            CALL dbt_create(t_2c_der_metric_prv(1, i_xyz), t_2c_tmp)
            CALL dbt_copy_matrix_to_tensor(t_2c_der_metric_prv(1, i_xyz), t_2c_tmp)

            CALL dbt_create(t_2c_template, t_2c_der_metric(i_xyz))
            CALL dbt_copy(t_2c_tmp, t_2c_der_metric(i_xyz), move_data=.TRUE.)

            CALL dbt_destroy(t_2c_tmp)
            CALL dbcsr_release(t_2c_der_metric_prv(1, i_xyz))
         END DO

      END IF

      CALL dbt_destroy(t_2c_template)
      CALL dbcsr_distribution_release(dbcsr_dist)
      DEALLOCATE (row_bsize, col_bsize)

      DO ibasis = 1, SIZE(basis_set_AO)
         orb_basis => basis_set_AO(ibasis)%gto_basis_set
         ri_basis => basis_set_RI(ibasis)%gto_basis_set
         CALL init_interaction_radii_orb_basis(orb_basis, dft_control%qs_control%eps_pgf_orb)
         CALL init_interaction_radii_orb_basis(ri_basis, dft_control%qs_control%eps_pgf_orb)
      END DO

      CALL timestop(handle)

   END SUBROUTINE precalc_derivatives

! **************************************************************************************************
!> \brief This routines calculates the force contribution from a trace over 3D tensors, i.e.
!>        force = sum_ijk A_ijk B_ijk. An iteration over the blocks is made, which index determin
!>        the atom on which the force acts
!> \param force ...
!> \param t_3c_contr ...
!> \param t_3c_der ...
!> \param atom_of_kind ...
!> \param kind_of ...
!> \param idx_to_at ...
!> \param pref ...
!> \param work_virial ...
!> \param cell ...
!> \param particle_set ...
!> \param do_mp2 ...
!> \param deriv_dim the dimension of the tensor corresponding to the derivative (by default 1)
! **************************************************************************************************
   SUBROUTINE get_force_from_3c_trace(force, t_3c_contr, t_3c_der, atom_of_kind, kind_of, idx_to_at, &
                                      pref, work_virial, cell, particle_set, do_mp2, deriv_dim)

      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(dbt_type), INTENT(INOUT)                      :: t_3c_contr
      TYPE(dbt_type), DIMENSION(3), INTENT(INOUT)        :: t_3c_der
      INTEGER, DIMENSION(:), INTENT(IN)                  :: atom_of_kind, kind_of, idx_to_at
      REAL(dp), INTENT(IN)                               :: pref
      REAL(dp), DIMENSION(3, 3), INTENT(INOUT), OPTIONAL :: work_virial
      TYPE(cell_type), OPTIONAL, POINTER                 :: cell
      TYPE(particle_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: particle_set
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_mp2
      INTEGER, INTENT(IN), OPTIONAL                      :: deriv_dim

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_force_from_3c_trace'

      INTEGER                                            :: deriv_dim_prv, handle, i_xyz, iat, &
                                                            iat_of_kind, ikind, j_xyz
      INTEGER, DIMENSION(3)                              :: ind
      LOGICAL                                            :: do_mp2_prv, found, use_virial
      REAL(dp)                                           :: new_force
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :), TARGET  :: contr_blk, der_blk
      REAL(dp), DIMENSION(3)                             :: scoord
      TYPE(dbt_iterator_type)                            :: iter

      CALL timeset(routineN, handle)

      use_virial = .FALSE.
      IF (PRESENT(work_virial) .AND. PRESENT(cell) .AND. PRESENT(particle_set)) use_virial = .TRUE.

      do_mp2_prv = .FALSE.
      IF (PRESENT(do_mp2)) do_mp2_prv = do_mp2

      deriv_dim_prv = 1
      IF (PRESENT(deriv_dim)) deriv_dim_prv = deriv_dim

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(t_3c_der,t_3c_contr,work_virial,force,use_virial,do_mp2_prv,deriv_dim_prv) &
!$OMP SHARED(pref,idx_to_at,atom_of_kind,kind_of,particle_set,cell) &
!$OMP PRIVATE(i_xyz,j_xyz,iter,ind,der_blk,contr_blk,found,new_force,iat,iat_of_kind,ikind,scoord)
      DO i_xyz = 1, 3
         CALL dbt_iterator_start(iter, t_3c_der(i_xyz))
         DO WHILE (dbt_iterator_blocks_left(iter))
            CALL dbt_iterator_next_block(iter, ind)

            CALL dbt_get_block(t_3c_der(i_xyz), ind, der_blk, found)
            CPASSERT(found)
            CALL dbt_get_block(t_3c_contr, ind, contr_blk, found)

            IF (found) THEN

               !take the trace of the blocks
               new_force = pref*SUM(der_blk(:, :, :)*contr_blk(:, :, :))

               !the first index of the derivative tensor defines the atom
               iat = idx_to_at(ind(deriv_dim_prv))
               iat_of_kind = atom_of_kind(iat)
               ikind = kind_of(iat)

               IF (.NOT. do_mp2_prv) THEN
!$OMP ATOMIC
                  force(ikind)%fock_4c(i_xyz, iat_of_kind) = force(ikind)%fock_4c(i_xyz, iat_of_kind) &
                                                             + new_force
               ELSE
!$OMP ATOMIC
                  force(ikind)%mp2_non_sep(i_xyz, iat_of_kind) = force(ikind)%mp2_non_sep(i_xyz, iat_of_kind) &
                                                                 + new_force
               END IF

               IF (use_virial) THEN
                  CALL real_to_scaled(scoord, particle_set(iat)%r, cell)

                  DO j_xyz = 1, 3
!$OMP ATOMIC
                     work_virial(i_xyz, j_xyz) = work_virial(i_xyz, j_xyz) + new_force*scoord(j_xyz)
                  END DO
               END IF
               DEALLOCATE (contr_blk)
            END IF
            DEALLOCATE (der_blk)
         END DO !iter
         CALL dbt_iterator_stop(iter)
      END DO
!$OMP END PARALLEL
      CALL timestop(handle)

   END SUBROUTINE get_force_from_3c_trace

! **************************************************************************************************
!> \brief Get the force from a contraction of type SUM_a,beta (a|beta') C_a,beta, where beta is an AO
!>        and a is a MO
!> \param force ...
!> \param t_mo_coeff ...
!> \param t_2c_MO_AO ...
!> \param atom_of_kind ...
!> \param kind_of ...
!> \param idx_to_at ...
!> \param pref ...
!> \param i_xyz ...
!> \param work_virial ...
!> \param cell ...
!> \param particle_set ...
! **************************************************************************************************
   SUBROUTINE get_MO_AO_force(force, t_mo_coeff, t_2c_MO_AO, atom_of_kind, kind_of, idx_to_at, &
                              pref, i_xyz, work_virial, cell, particle_set)

      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(dbt_type), INTENT(INOUT)                      :: t_mo_coeff, t_2c_MO_AO
      INTEGER, DIMENSION(:), INTENT(IN)                  :: atom_of_kind, kind_of, idx_to_at
      REAL(dp), INTENT(IN)                               :: pref
      INTEGER, INTENT(IN)                                :: i_xyz
      REAL(dp), DIMENSION(3, 3), INTENT(INOUT), OPTIONAL :: work_virial
      TYPE(cell_type), OPTIONAL, POINTER                 :: cell
      TYPE(particle_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: particle_set

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'get_MO_AO_force'

      INTEGER                                            :: handle, iat, iat_of_kind, ikind, j_xyz
      INTEGER, DIMENSION(2)                              :: ind
      LOGICAL                                            :: found, use_virial
      REAL(dp)                                           :: new_force
      REAL(dp), ALLOCATABLE, DIMENSION(:, :), TARGET     :: mo_ao_blk, mo_coeff_blk
      REAL(dp), DIMENSION(3)                             :: scoord
      TYPE(dbt_iterator_type)                            :: iter

      CALL timeset(routineN, handle)

      use_virial = .FALSE.
      IF (PRESENT(work_virial) .AND. PRESENT(cell) .AND. PRESENT(particle_set)) use_virial = .TRUE.

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(t_2c_MO_AO,t_mo_coeff,pref,use_virial,particle_set,cell,force,work_virial) &
!$OMP SHARED(idx_to_at,atom_of_kind,kind_of,i_xyz) &
!$OMP PRIVATE(iter,ind,mo_ao_blk,mo_coeff_blk,found,new_force,iat,iat_of_kind,ikind,scoord,j_xyz)
      CALL dbt_iterator_start(iter, t_2c_MO_AO)
      DO WHILE (dbt_iterator_blocks_left(iter))
         CALL dbt_iterator_next_block(iter, ind)

         CALL dbt_get_block(t_2c_MO_AO, ind, mo_ao_blk, found)
         CPASSERT(found)
         CALL dbt_get_block(t_mo_coeff, ind, mo_coeff_blk, found)

         IF (found) THEN

            new_force = pref*SUM(mo_ao_blk(:, :)*mo_coeff_blk(:, :))

            iat = idx_to_at(ind(2)) !AO index is column index
            iat_of_kind = atom_of_kind(iat)
            ikind = kind_of(iat)

!$OMP ATOMIC
            force(ikind)%fock_4c(i_xyz, iat_of_kind) = force(ikind)%fock_4c(i_xyz, iat_of_kind) &
                                                       + new_force

            IF (use_virial) THEN

               CALL real_to_scaled(scoord, particle_set(iat)%r, cell)

               DO j_xyz = 1, 3
!$OMP ATOMIC
                  work_virial(i_xyz, j_xyz) = work_virial(i_xyz, j_xyz) + new_force*scoord(j_xyz)
               END DO
            END IF

            DEALLOCATE (mo_coeff_blk)
         END IF

         DEALLOCATE (mo_ao_blk)
      END DO !iter
      CALL dbt_iterator_stop(iter)
!$OMP END PARALLEL

      CALL timestop(handle)

   END SUBROUTINE get_MO_AO_force

! **************************************************************************************************
!> \brief Update the forces due to the derivative of the a 2-center product d/dR (Q|R)
!> \param force ...
!> \param t_2c_contr A precontracted tensor containing sum_abcdPS (ab|P)(P|Q)^-1 (R|S)^-1 (S|cd) P_ac P_bd
!> \param t_2c_der the d/dR (Q|R) tensor, in all 3 cartesian directions
!> \param atom_of_kind ...
!> \param kind_of ...
!> \param idx_to_at ...
!> \param pref ...
!> \param work_virial ...
!> \param cell ...
!> \param particle_set ...
!> \param do_mp2 ...
!> \param do_ovlp ...
!> \note IMPORTANT: t_tc_contr and t_2c_der need to have the same distribution
! **************************************************************************************************
   SUBROUTINE get_2c_der_force(force, t_2c_contr, t_2c_der, atom_of_kind, kind_of, idx_to_at, &
                               pref, work_virial, cell, particle_set, do_mp2, do_ovlp)

      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(dbt_type), INTENT(INOUT)                      :: t_2c_contr
      TYPE(dbt_type), DIMENSION(3), INTENT(INOUT)        :: t_2c_der
      INTEGER, DIMENSION(:), INTENT(IN)                  :: atom_of_kind, kind_of, idx_to_at
      REAL(dp), INTENT(IN)                               :: pref
      REAL(dp), DIMENSION(3, 3), INTENT(INOUT), OPTIONAL :: work_virial
      TYPE(cell_type), OPTIONAL, POINTER                 :: cell
      TYPE(particle_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: particle_set
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_mp2, do_ovlp

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'get_2c_der_force'

      INTEGER                                            :: handle, i_xyz, iat, iat_of_kind, ikind, &
                                                            j_xyz, jat, jat_of_kind, jkind
      INTEGER, DIMENSION(2)                              :: ind
      LOGICAL                                            :: do_mp2_prv, do_ovlp_prv, found, &
                                                            use_virial
      REAL(dp)                                           :: new_force
      REAL(dp), ALLOCATABLE, DIMENSION(:, :), TARGET     :: contr_blk, der_blk
      REAL(dp), DIMENSION(3)                             :: scoord
      TYPE(dbt_iterator_type)                            :: iter

      !Loop over the blocks of d/dR (Q|R), contract with the corresponding block of t_2c_contr and
      !update the relevant force

      CALL timeset(routineN, handle)

      use_virial = .FALSE.
      IF (PRESENT(work_virial) .AND. PRESENT(cell) .AND. PRESENT(particle_set)) use_virial = .TRUE.

      do_mp2_prv = .FALSE.
      IF (PRESENT(do_mp2)) do_mp2_prv = do_mp2

      do_ovlp_prv = .FALSE.
      IF (PRESENT(do_ovlp)) do_ovlp_prv = do_ovlp

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(t_2c_der,t_2c_contr,work_virial,force,use_virial,do_mp2_prv,do_ovlp_prv) &
!$OMP SHARED(pref,idx_to_at,atom_of_kind,kind_of,particle_set,cell) &
!$OMP PRIVATE(i_xyz,j_xyz,iter,ind,der_blk,contr_blk,found,new_force) &
!$OMP PRIVATE(iat,jat,iat_of_kind,jat_of_kind,ikind,jkind,scoord)
      DO i_xyz = 1, 3
         CALL dbt_iterator_start(iter, t_2c_der(i_xyz))
         DO WHILE (dbt_iterator_blocks_left(iter))
            CALL dbt_iterator_next_block(iter, ind)

            IF (ind(1) == ind(2)) CYCLE

            CALL dbt_get_block(t_2c_der(i_xyz), ind, der_blk, found)
            CPASSERT(found)
            CALL dbt_get_block(t_2c_contr, ind, contr_blk, found)

            IF (found) THEN

               !an element of d/dR (Q|R) corresponds to 2 things because of translational invariance
               !(Q'| R) = - (Q| R'), once wrt the center on Q, and once on R
               new_force = pref*SUM(der_blk(:, :)*contr_blk(:, :))

               iat = idx_to_at(ind(1))
               iat_of_kind = atom_of_kind(iat)
               ikind = kind_of(iat)

               IF (do_mp2_prv) THEN
!$OMP ATOMIC
                  force(ikind)%mp2_non_sep(i_xyz, iat_of_kind) = force(ikind)%mp2_non_sep(i_xyz, iat_of_kind) &
                                                                 + new_force
               ELSE IF (do_ovlp_prv) THEN
!$OMP ATOMIC
                  force(ikind)%overlap(i_xyz, iat_of_kind) = force(ikind)%overlap(i_xyz, iat_of_kind) &
                                                             + new_force
               ELSE
!$OMP ATOMIC
                  force(ikind)%fock_4c(i_xyz, iat_of_kind) = force(ikind)%fock_4c(i_xyz, iat_of_kind) &
                                                             + new_force
               END IF

               IF (use_virial) THEN

                  CALL real_to_scaled(scoord, particle_set(iat)%r, cell)

                  DO j_xyz = 1, 3
!$OMP ATOMIC
                     work_virial(i_xyz, j_xyz) = work_virial(i_xyz, j_xyz) + new_force*scoord(j_xyz)
                  END DO
               END IF

               jat = idx_to_at(ind(2))
               jat_of_kind = atom_of_kind(jat)
               jkind = kind_of(jat)

               IF (do_mp2_prv) THEN
!$OMP ATOMIC
                  force(jkind)%mp2_non_sep(i_xyz, jat_of_kind) = force(jkind)%mp2_non_sep(i_xyz, jat_of_kind) &
                                                                 - new_force
               ELSE IF (do_ovlp_prv) THEN
!$OMP ATOMIC
                  force(jkind)%overlap(i_xyz, jat_of_kind) = force(jkind)%overlap(i_xyz, jat_of_kind) &
                                                             - new_force
               ELSE
!$OMP ATOMIC
                  force(jkind)%fock_4c(i_xyz, jat_of_kind) = force(jkind)%fock_4c(i_xyz, jat_of_kind) &
                                                             - new_force
               END IF

               IF (use_virial) THEN

                  CALL real_to_scaled(scoord, particle_set(jat)%r, cell)

                  DO j_xyz = 1, 3
!$OMP ATOMIC
                     work_virial(i_xyz, j_xyz) = work_virial(i_xyz, j_xyz) - new_force*scoord(j_xyz)
                  END DO
               END IF

               DEALLOCATE (contr_blk)
            END IF

            DEALLOCATE (der_blk)
         END DO !iter
         CALL dbt_iterator_stop(iter)

      END DO !i_xyz
!$OMP END PARALLEL
      CALL timestop(handle)

   END SUBROUTINE get_2c_der_force

! **************************************************************************************************
!> \brief a small utility function that returns the atom corresponding to a block of a split tensor
!> \param idx_to_at ...
!> \param bsizes_split ...
!> \param bsizes_orig ...
!> \return ...
! **************************************************************************************************
   SUBROUTINE get_idx_to_atom(idx_to_at, bsizes_split, bsizes_orig)
      INTEGER, DIMENSION(:), INTENT(INOUT)               :: idx_to_at
      INTEGER, DIMENSION(:), INTENT(IN)                  :: bsizes_split, bsizes_orig

      INTEGER                                            :: full_sum, iat, iblk, split_sum

      iat = 1
      full_sum = bsizes_orig(iat)
      split_sum = 0
      DO iblk = 1, SIZE(bsizes_split)
         split_sum = split_sum + bsizes_split(iblk)

         IF (split_sum .GT. full_sum) THEN
            iat = iat + 1
            full_sum = full_sum + bsizes_orig(iat)
         END IF

         idx_to_at(iblk) = iat
      END DO

   END SUBROUTINE get_idx_to_atom

! **************************************************************************************************
!> \brief Function for calculating sqrt of a matrix
!> \param values ...
!> \return ...
! **************************************************************************************************
   FUNCTION my_sqrt(values)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: values
      REAL(KIND=dp), DIMENSION(SIZE(values))             :: my_sqrt

      my_sqrt = SQRT(values)
   END FUNCTION

! **************************************************************************************************
!> \brief Function for calculation inverse sqrt of a matrix
!> \param values ...
!> \return ...
! **************************************************************************************************
   FUNCTION my_invsqrt(values)
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: values
      REAL(KIND=dp), DIMENSION(SIZE(values))             :: my_invsqrt

      my_invsqrt = SQRT(1.0_dp/values)
   END FUNCTION

END MODULE
