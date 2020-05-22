!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief RI-methods for HFX
! **************************************************************************************************

MODULE hfx_ri

   USE arnoldi_api,                     ONLY: arnoldi_extremal
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE cp_array_utils,                  ONLY: cp_1d_r_p_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_cholesky,               ONLY: cp_dbcsr_cholesky_decompose,&
                                              cp_dbcsr_cholesky_invert
   USE cp_dbcsr_diag,                   ONLY: cp_dbcsr_power
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              cp_dbcsr_dist2d_to_dist
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_add_on_diag, dbcsr_copy, dbcsr_create, dbcsr_distribution_get, &
        dbcsr_distribution_release, dbcsr_distribution_type, dbcsr_dot, dbcsr_filter, &
        dbcsr_frobenius_norm, dbcsr_get_info, dbcsr_get_num_blocks, dbcsr_multiply, dbcsr_p_type, &
        dbcsr_release, dbcsr_scalar, dbcsr_scale, dbcsr_type, dbcsr_type_no_symmetry, &
        dbcsr_type_symmetric
   USE dbcsr_tensor_api,                ONLY: &
        dbcsr_t_batched_contract_finalize, dbcsr_t_batched_contract_init, dbcsr_t_clear, &
        dbcsr_t_contract, dbcsr_t_copy, dbcsr_t_copy_matrix_to_tensor, &
        dbcsr_t_copy_tensor_to_matrix, dbcsr_t_create, dbcsr_t_destroy, dbcsr_t_filter, &
        dbcsr_t_get_info, dbcsr_t_get_num_blocks, dbcsr_t_get_num_blocks_total, &
        dbcsr_t_mp_environ_pgrid, dbcsr_t_nd_mp_comm, dbcsr_t_pgrid_change_dims, &
        dbcsr_t_pgrid_create, dbcsr_t_pgrid_destroy, dbcsr_t_pgrid_type, &
        dbcsr_t_reserved_block_indices, dbcsr_t_type
   USE distribution_2d_types,           ONLY: distribution_2d_type
   USE hfx_types,                       ONLY: hfx_ri_type
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
        build_2c_integrals, build_2c_neighbor_lists, build_3c_integrals, build_3c_neighbor_lists, &
        compress_tensor, decompress_tensor, get_tensor_occupancy, neighbor_list_3c_destroy, &
        tensor_change_pgrid
   USE qs_tensors_types,                ONLY: create_2c_tensor,&
                                              create_3c_tensor,&
                                              create_tensor_batches,&
                                              distribution_3d_create,&
                                              distribution_3d_type,&
                                              neighbor_list_3c_type,&
                                              split_block_sizes
   USE util,                            ONLY: sort
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: hfx_ri_update_ks

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'hfx_ri'
CONTAINS

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

      CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_ri_pre_scf_mo', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, handle2, ispin, n_dependent, &
                                                            unit_nr, unit_nr_dbcsr
      REAL(KIND=dp)                                      :: threshold
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_t_type), DIMENSION(1)                   :: t_2c_int
      TYPE(dbcsr_type), DIMENSION(1)                     :: t_2c_int_mat, t_2c_op_pot, &
                                                            t_2c_op_pot_sqrt, &
                                                            t_2c_op_pot_sqrt_inv, t_2c_op_RI, &
                                                            t_2c_op_RI_inv

      CALL timeset(routineN, handle)

      unit_nr_dbcsr = ri_data%unit_nr_dbcsr
      unit_nr = ri_data%unit_nr

      CALL get_qs_env(qs_env, para_env=para_env, blacs_env=blacs_env)

      CALL timeset(routineN//"_int", handle2)

      CALL hfx_ri_pre_scf_calc_tensors(qs_env, ri_data, t_2c_op_RI, t_2c_op_pot, ri_data%t_3c_int)

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
         ENDIF

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
         IF (ri_data%check_2c_inv) THEN
            CALL check_sqrt(t_2c_op_pot(1), matrix_sqrt=t_2c_op_pot_sqrt(1), unit_nr=unit_nr)
         ENDIF
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
         ENDIF
      ENDIF

      CALL dbcsr_release(t_2c_op_pot(1))

      CALL dbcsr_t_create(t_2c_int_mat(1), t_2c_int(1), name="(RI|RI)")
      CALL dbcsr_t_copy_matrix_to_tensor(t_2c_int_mat(1), t_2c_int(1))
      CALL dbcsr_release(t_2c_int_mat(1))
      DO ispin = 1, nspins
         IF (ispin > 1) CALL dbcsr_t_create(ri_data%t_2c_int(1, 1), ri_data%t_2c_int(ispin, 1))
         CALL dbcsr_t_copy(t_2c_int(1), ri_data%t_2c_int(ispin, 1))
      ENDDO
      CALL dbcsr_t_destroy(t_2c_int(1))
      CALL timestop(handle2)

      CALL timeset(routineN//"_3c", handle2)
      CALL dbcsr_t_copy(ri_data%t_3c_int(1, 1), ri_data%t_3c_int_ctr_1(1, 1), order=[2, 1, 3], move_data=.TRUE.)
      CALL dbcsr_t_destroy(ri_data%t_3c_int(1, 1))
      CALL dbcsr_t_filter(ri_data%t_3c_int_ctr_1(1, 1), ri_data%filter_eps)
      CALL dbcsr_t_copy(ri_data%t_3c_int_ctr_1(1, 1), ri_data%t_3c_int_ctr_2(1, 1))
      DEALLOCATE (ri_data%t_3c_int)
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
      ENDIF

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
      ENDIF

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
      ENDIF
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
         ENDIF
         CALL dbcsr_release(matrix_tmp)
         CALL dbcsr_release(matrix_copy)
      ENDIF

      IF (PRESENT(matrix_sqrt_inv)) THEN
         CALL dbcsr_create(matrix_tmp, template=matrix)
         CALL dbcsr_copy(matrix_copy, matrix_sqrt_inv)
         CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_sqrt_inv, matrix_copy, &
                             0.0_dp, matrix_tmp)
         CALL check_inverse(matrix_tmp, matrix, name="SQRT("//TRIM(name_prv)//")", unit_nr=unit_nr)
         CALL dbcsr_release(matrix_tmp)
         CALL dbcsr_release(matrix_copy)
      ENDIF

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
      TYPE(dbcsr_t_type), DIMENSION(1, 1), INTENT(OUT)   :: t_3c_int

      CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_ri_pre_scf_calc_tensors', &
         routineP = moduleN//':'//routineN

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
      TYPE(dbcsr_t_type)                                 :: t_3c_tmp
      TYPE(dbcsr_t_type), DIMENSION(1, 1)                :: t_3c_int_batched
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(distribution_2d_type), POINTER                :: dist_2d
      TYPE(distribution_3d_type)                         :: dist_3d
      TYPE(gto_basis_set_p_type), ALLOCATABLE, &
         DIMENSION(:), TARGET                            :: basis_set_AO, basis_set_RI
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis
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
         ! interaction radii should be based on eps_pgf_orb controlled in RI section
         ! (since hartree-fock needs very tight eps_pgf_orb for Kohn-Sham/Fock matrix but eps_pgf_orb
         ! can be much looser in RI HFX since no systematic error is introduced with tensor sparsity)
         CALL init_interaction_radii_orb_basis(orb_basis, ri_data%eps_pgf_orb)
      ENDDO

      n_mem = ri_data%n_mem
      CALL create_tensor_batches(sizes_AO, n_mem, starts_array_mc_int, ends_array_mc_int, &
                                 starts_array_mc_block_int, ends_array_mc_block_int)

      DEALLOCATE (starts_array_mc_int, ends_array_mc_int)

      CALL create_3c_tensor(t_3c_int_batched(1, 1), dist_RI, dist_AO_1, dist_AO_2, ri_data%pgrid, &
                            sizes_RI, sizes_AO, sizes_AO, map1=[1], map2=[2, 3], &
                            starts_array_block_2=starts_array_mc_block_int, ends_array_block_2=ends_array_mc_block_int, &
                            name="(RI | AO AO)")

      CALL get_qs_env(qs_env, nkind=nkind, particle_set=particle_set, atomic_kind_set=atomic_kind_set)
      CALL dbcsr_t_mp_environ_pgrid(ri_data%pgrid, pdims, pcoord)
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
                                   "HFX_3c_nl", &
                                   qs_env, op_pos=1, sym_jk=.TRUE., own_dist=.TRUE.)

      DO i_mem = 1, n_mem
         CALL build_3c_integrals(t_3c_int_batched, ri_data%filter_eps/2, qs_env, nl_3c, &
                                 basis_set_RI, basis_set_AO, basis_set_AO, &
                                 ri_data%ri_metric, int_eps=ri_data%eps_schwarz, op_pos=1, &
                                 desymmetrize=.FALSE., &
                                 bounds_j=[starts_array_mc_block_int(i_mem), ends_array_mc_block_int(i_mem)])
         CALL dbcsr_t_copy(t_3c_int_batched(1, 1), t_3c_int(1, 1), summation=.TRUE., move_data=.TRUE.)
         CALL dbcsr_t_filter(t_3c_int(1, 1), ri_data%filter_eps/2)
      ENDDO

      CALL dbcsr_t_destroy(t_3c_int_batched(1, 1))

      CALL neighbor_list_3c_destroy(nl_3c)

      CALL dbcsr_t_create(t_3c_int(1, 1), t_3c_tmp)

      IF (ri_data%flavor == ri_pmat) THEN ! desymmetrize
         ! desymmetrize
         CALL dbcsr_t_copy(t_3c_int(1, 1), t_3c_tmp)
         CALL dbcsr_t_copy(t_3c_tmp, t_3c_int(1, 1), order=[1, 3, 2], summation=.TRUE., move_data=.TRUE.)

         ! For RI-RHO filter_eps_storage is reserved for screening tensor contracted with RI-metric
         ! with RI metric but not to bare integral tensor
         CALL dbcsr_t_filter(t_3c_int(1, 1), ri_data%filter_eps)
      ELSE
         CALL dbcsr_t_filter(t_3c_int(1, 1), ri_data%filter_eps_storage/2)
      ENDIF

      CALL dbcsr_t_destroy(t_3c_tmp)

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
      ENDIF

      DO ibasis = 1, SIZE(basis_set_AO)
         orb_basis => basis_set_AO(ibasis)%gto_basis_set
         ! reset interaction radii of orb basis
         CALL init_interaction_radii_orb_basis(orb_basis, dft_control%qs_control%eps_pgf_orb)
      ENDDO

      IF (ri_data%calc_condnum) THEN
         CALL arnoldi_extremal(t_2c_int_pot(1), max_ev, min_ev, threshold=ri_data%eps_lanczos, &
                               max_iter=ri_data%max_iter_lanczos, converged=converged)

         IF (.NOT. converged) THEN
            CPWARN("Condition number estimate of (P|Q) (HFX potential) is not reliable (not converged).")
         ENDIF

         IF (ri_data%unit_nr > 0) THEN
            WRITE (ri_data%unit_nr, '(T2,A)') "2-Norm Condition Number of (P|Q) integrals (HFX potential)"
            IF (min_ev > 0) THEN
               WRITE (ri_data%unit_nr, '(T4,A,ES11.3E3,T32,A,ES11.3E3,A4,ES11.3E3,T63,A,F8.4)') &
                  "CN : max/min ev: ", max_ev, " / ", min_ev, "=", max_ev/min_ev, "Log(2-CN):", LOG10(max_ev/min_ev)
            ELSE
               WRITE (ri_data%unit_nr, '(T4,A,ES11.3E3,T32,A,ES11.3E3,T63,A)') &
                  "CN : max/min ev: ", max_ev, " / ", min_ev, "Log(CN): infinity"
            ENDIF
         ENDIF

         IF (.NOT. ri_data%same_op) THEN
            CALL arnoldi_extremal(t_2c_int_RI(1), max_ev, min_ev, threshold=ri_data%eps_lanczos, &
                                  max_iter=ri_data%max_iter_lanczos, converged=converged)

            IF (.NOT. converged) THEN
               CPWARN("Condition number estimate of (P|Q) matrix (RI metric) is not reliable (not converged).")
            ENDIF

            IF (ri_data%unit_nr > 0) THEN
               WRITE (ri_data%unit_nr, '(T2,A)') "2-Norm Condition Number of (P|Q) integrals (RI metric)"
               IF (min_ev > 0) THEN
                  WRITE (ri_data%unit_nr, '(T4,A,ES11.3E3,T32,A,ES11.3E3,A4,ES11.3E3,T63,A,F8.4)') &
                     "CN : max/min ev: ", max_ev, " / ", min_ev, "=", max_ev/min_ev, "Log(2-CN):", LOG10(max_ev/min_ev)
               ELSE
                  WRITE (ri_data%unit_nr, '(T4,A,ES11.3E3,T32,A,ES11.3E3,T63,A)') &
                     "CN : max/min ev: ", max_ev, " / ", min_ev, "Log(CN): infinity"
               ENDIF
            ENDIF
         ENDIF
      ENDIF

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

      CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_ri_pre_scf_Pmat', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, handle2, i_mem, j_mem, &
                                                            n_dependent, unit_nr, unit_nr_dbcsr
      INTEGER(int_8)                                     :: nflop
      INTEGER, DIMENSION(2, 1)                           :: bounds_i
      INTEGER, DIMENSION(2, 2)                           :: bounds_j
      INTEGER, DIMENSION(3)                              :: dims_3c
      REAL(KIND=dp)                                      :: threshold
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_t_type)                                 :: t_3c_2
      TYPE(dbcsr_t_type), DIMENSION(1)                   :: t_2c_int
      TYPE(dbcsr_t_type), DIMENSION(1, 1)                :: t_3c_int_1
      TYPE(dbcsr_type), DIMENSION(1)                     :: t_2c_int_mat, t_2c_op_pot, t_2c_op_RI, &
                                                            t_2c_tmp, t_2c_tmp_2

      CALL timeset(routineN, handle)

      unit_nr_dbcsr = ri_data%unit_nr_dbcsr
      unit_nr = ri_data%unit_nr

      CALL get_qs_env(qs_env, para_env=para_env, blacs_env=blacs_env)

      CALL timeset(routineN//"_int", handle2)
      CALL hfx_ri_pre_scf_calc_tensors(qs_env, ri_data, t_2c_op_RI, t_2c_op_pot, t_3c_int_1)

      CALL dbcsr_t_copy(t_3c_int_1(1, 1), ri_data%t_3c_int_ctr_3(1, 1), order=[1, 2, 3], move_data=.TRUE.)

      CALL dbcsr_t_destroy(t_3c_int_1(1, 1))

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
      ENDIF

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
      ENDIF

      CALL dbcsr_t_create(t_2c_int_mat(1), t_2c_int(1), name="(RI|RI)")
      CALL dbcsr_t_copy_matrix_to_tensor(t_2c_int_mat(1), t_2c_int(1))
      CALL dbcsr_release(t_2c_int_mat(1))
      CALL dbcsr_t_copy(t_2c_int(1), ri_data%t_2c_int(1, 1))
      CALL dbcsr_t_destroy(t_2c_int(1))
      CALL dbcsr_t_filter(ri_data%t_2c_int(1, 1), ri_data%filter_eps)

      CALL timestop(handle2)

      CALL dbcsr_t_create(ri_data%t_3c_int_ctr_3(1, 1), t_3c_2)

      CALL dbcsr_t_get_info(ri_data%t_3c_int_ctr_3(1, 1), nfull_total=dims_3c)

      DO i_mem = 1, ri_data%n_mem
         bounds_i(:, 1) = [ri_data%starts_array_RI_mem(i_mem), ri_data%ends_array_RI_mem(i_mem)]
         CALL dbcsr_t_batched_contract_init(ri_data%t_2c_int(1, 1))
         DO j_mem = 1, ri_data%n_mem
            bounds_j(:, 1) = [ri_data%starts_array_mem(j_mem), ri_data%ends_array_mem(j_mem)]
            bounds_j(:, 2) = [1, dims_3c(3)]
            CALL timeset(routineN//"_RIx3C", handle2)
            CALL dbcsr_t_contract(dbcsr_scalar(1.0_dp), ri_data%t_2c_int(1, 1), ri_data%t_3c_int_ctr_3(1, 1), &
                                  dbcsr_scalar(0.0_dp), t_3c_2, &
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
            CALL dbcsr_t_copy(t_3c_2, ri_data%t_3c_int_ctr_1(j_mem, i_mem), order=[2, 1, 3], move_data=.TRUE.)

            IF (ALLOCATED(ri_data%blk_indices(j_mem, i_mem)%ind)) DEALLOCATE (ri_data%blk_indices(j_mem, i_mem)%ind)
            ALLOCATE (ri_data%blk_indices(j_mem, i_mem)%ind(dbcsr_t_get_num_blocks(ri_data%t_3c_int_ctr_1(j_mem, i_mem)), 3))
            CALL dbcsr_t_reserved_block_indices(ri_data%t_3c_int_ctr_1(j_mem, i_mem), ri_data%blk_indices(j_mem, i_mem)%ind)
            CALL compress_tensor(ri_data%t_3c_int_ctr_1(j_mem, i_mem), ri_data%store_3c(j_mem, i_mem), ri_data%filter_eps_storage)

            CALL timestop(handle2)
         ENDDO
         CALL dbcsr_t_batched_contract_finalize(ri_data%t_2c_int(1, 1))
      ENDDO

      CALL dbcsr_t_clear(ri_data%t_2c_int(1, 1))
      CALL dbcsr_t_destroy(t_3c_2)

      CALL dbcsr_t_copy(ri_data%t_3c_int_ctr_3(1, 1), ri_data%t_3c_int_ctr_2(1, 1), order=[2, 1, 3], move_data=.TRUE.)

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
         ENDDO
         ENDIF

         ncols = end_ind - start_ind + 1
         ALLOCATE (sort_2(ncols))
         ALLOCATE (ind_2(ncols))
         sort_2(:) = blk_ind(start_ind:end_ind, 2)
         CALL sort(sort_2, ncols, ind_2)
         ind_2 = ind_2 + start_ind - 1

         blk_ind(start_ind:end_ind, :) = blk_ind(ind_2, :)
         start_ind = end_ind + 1

         DEALLOCATE (sort_2, ind_2)
      ENDDO

      ALLOCATE (blk_ind_tmp(nblk, 2))
      blk_ind_tmp = 0

      iblk = 0
      DO iblk_all = 1, nblk
         IF (iblk >= 1) THEN
            IF (ALL(blk_ind_tmp(iblk, :) == blk_ind(iblk_all, :))) THEN
               CYCLE
            ENDIF
         ENDIF
         iblk = iblk + 1
         blk_ind_tmp(iblk, :) = blk_ind(iblk_all, :)
      ENDDO
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
      TYPE(dbcsr_p_type), DIMENSION(:, :), INTENT(INOUT) :: ks_matrix
      REAL(KIND=dp), INTENT(OUT)                         :: ehfx
      TYPE(mo_set_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: mos
      TYPE(dbcsr_p_type), DIMENSION(:, :)                :: rho_ao
      LOGICAL, INTENT(IN)                                :: geometry_did_change
      INTEGER, INTENT(IN)                                :: nspins
      REAL(KIND=dp), INTENT(IN)                          :: hf_fraction

      CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_ri_update_ks', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, handle2, ispin
      INTEGER(int_8)                                     :: nblks
      INTEGER, DIMENSION(2)                              :: homo
      REAL(dp)                                           :: etmp, fac
      REAL(KIND=dp), DIMENSION(:), POINTER               :: mo_eigenvalues
      TYPE(cp_1d_r_p_type), DIMENSION(:), POINTER        :: occupied_evals
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: homo_localized, moloc_coeff, &
                                                            occupied_orbs
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_type), DIMENSION(2)                     :: mo_coeff_b
      TYPE(dbcsr_type), POINTER                          :: mo_coeff_b_tmp
      TYPE(mo_set_type), POINTER                         :: mo_set

      CALL timeset(routineN, handle)

      IF (nspins == 1) THEN
         fac = 0.5_dp*hf_fraction
      ELSE
         fac = 1.0_dp*hf_fraction
      END IF

      SELECT CASE (ri_data%flavor)
      CASE (ri_mo)
         CPASSERT(PRESENT(mos))
         CALL timeset(routineN//"_MO", handle2)

         IF (ri_data%do_loc) THEN
            ALLOCATE (occupied_orbs(nspins))
            ALLOCATE (occupied_evals(nspins))
            ALLOCATE (homo_localized(nspins))
         ENDIF
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
            ENDIF

            IF (ri_data%do_loc) THEN
               occupied_orbs(ispin)%matrix => mo_coeff
               occupied_evals(ispin)%array => mo_eigenvalues
               CALL cp_fm_create(homo_localized(ispin)%matrix, occupied_orbs(ispin)%matrix%matrix_struct)
               CALL cp_fm_to_fm(occupied_orbs(ispin)%matrix, homo_localized(ispin)%matrix)
            ENDIF
         ENDDO

         IF (ri_data%do_loc) THEN
            CALL qs_loc_env_create(ri_data%qs_loc_env)
            CALL qs_loc_control_init(ri_data%qs_loc_env, ri_data%loc_subsection, do_homo=.TRUE.)
            CALL qs_loc_init(qs_env, ri_data%qs_loc_env, ri_data%loc_subsection, homo_localized)
            DO ispin = 1, nspins
               CALL qs_loc_driver(qs_env, ri_data%qs_loc_env, ri_data%print_loc_subsection, ispin, &
                                  ext_mo_coeff=homo_localized(ispin)%matrix)
            ENDDO
            CALL get_qs_loc_env(qs_loc_env=ri_data%qs_loc_env, moloc_coeff=moloc_coeff)

            DO ispin = 1, nspins
               CALL cp_fm_release(homo_localized(ispin)%matrix)
            ENDDO

            DEALLOCATE (occupied_orbs, occupied_evals, homo_localized)

         ENDIF

         DO ispin = 1, nspins
            mo_set => mos(ispin)%mo_set
            IF (ri_data%do_loc) THEN
               CALL copy_fm_to_dbcsr(moloc_coeff(ispin)%matrix, mo_coeff_b(ispin))
            ENDIF
            CALL dbcsr_scale(mo_coeff_b(ispin), SQRT(mo_set%maxocc))
            homo(ispin) = mo_set%homo
         ENDDO

         IF (ri_data%do_loc) CALL qs_loc_env_release(ri_data%qs_loc_env)
         CALL timestop(handle2)

         CALL hfx_ri_update_ks_mo(qs_env, ri_data, ks_matrix, mo_coeff_b, homo, &
                                  geometry_did_change, nspins)
      CASE (ri_pmat)

         NULLIFY (para_env)
         CALL get_qs_env(qs_env, para_env=para_env)
         DO ispin = 1, SIZE(rho_ao, 1)
            nblks = dbcsr_get_num_blocks(rho_ao(ispin, 1)%matrix)
            CALL mp_sum(nblks, para_env%group)
            IF (nblks == 0) THEN
               CPABORT("received empty density matrix")
            ENDIF
         ENDDO

         CALL hfx_ri_update_ks_pmat(qs_env, ri_data, ks_matrix, rho_ao, &
                                    geometry_did_change, nspins)

      END SELECT

      DO ispin = 1, nspins
         CALL dbcsr_release(mo_coeff_b(ispin))
      ENDDO

      DO ispin = 1, nspins
         CALL dbcsr_scale(ks_matrix(ispin, 1)%matrix, -fac)
         CALL dbcsr_filter(ks_matrix(ispin, 1)%matrix, ri_data%filter_eps)
      ENDDO

      CALL timeset(routineN//"_energy", handle2)
      ! Calculate the exchange energy
      ehfx = 0.0_dp
      DO ispin = 1, nspins
         CALL dbcsr_dot(ks_matrix(ispin, 1)%matrix, rho_ao(ispin, 1)%matrix, &
                        etmp)
         ehfx = ehfx + 0.5_dp*etmp

      ENDDO
      CALL timestop(handle2)

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
! **************************************************************************************************
   SUBROUTINE hfx_ri_update_ks_mo(qs_env, ri_data, ks_matrix, mo_coeff, &
                                  homo, geometry_did_change, nspins)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(dbcsr_p_type), DIMENSION(:, :)                :: ks_matrix
      TYPE(dbcsr_type), DIMENSION(:), INTENT(IN)         :: mo_coeff
      INTEGER, DIMENSION(:)                              :: homo
      LOGICAL, INTENT(IN)                                :: geometry_did_change
      INTEGER, INTENT(IN)                                :: nspins

      CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_ri_update_ks_mo', &
         routineP = moduleN//':'//routineN
      REAL(dp), PARAMETER                                :: lb_ratio = 0.1_dp

      INTEGER :: bsize, bsum, comm_2d, count, handle, handle2, i_mem, iblock, iproc, ispin, &
         max_pdim_mo, n_mem, n_mos, nblock, nproc, nproc_rem, pdim_AO, pdim_AO_, pdim_AO__, &
         pdim_mo, tdim_mo, unit_nr_dbcsr
      INTEGER(int_8)                                     :: nblks, nflop
      INTEGER, ALLOCATABLE, DIMENSION(:) :: dist1, dist2, dist3, mem_end, mem_end_block_1, &
         mem_end_block_2, mem_size, mem_start, mem_start_block_1, mem_start_block_2, mo_bsizes_1, &
         mo_bsizes_2
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: bounds
      INTEGER, DIMENSION(2)                              :: pdims_2d
      INTEGER, DIMENSION(3)                              :: pcoord, pdims, pdims_AO, pdims_RI
      LOGICAL                                            :: do_initialize
      REAL(dp)                                           :: t1, t2
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_distribution_type)                      :: ks_dist
      TYPE(dbcsr_t_pgrid_type)                           :: pgrid, pgrid_2d
      TYPE(dbcsr_t_pgrid_type), POINTER                  :: pgrid_opt_KS, pgrid_opt_mo_L, &
                                                            pgrid_opt_mo_R, pgrid_opt_RI
      TYPE(dbcsr_t_type)                                 :: ks_t, ks_t_mat, mo_coeff_t, &
                                                            mo_coeff_t_split
      TYPE(dbcsr_t_type), DIMENSION(1, 1)                :: t_3c_int_mo_1, t_3c_int_mo_2

      CALL timeset(routineN, handle)

      CPASSERT(SIZE(ks_matrix, 2) == 1)

      NULLIFY (pgrid_opt_RI, pgrid_opt_KS, pgrid_opt_mo_L, pgrid_opt_mo_R)

      unit_nr_dbcsr = ri_data%unit_nr_dbcsr

      IF (geometry_did_change) THEN
         CALL hfx_ri_pre_scf_mo(qs_env, ri_data, nspins)
      ENDIF

      nblks = dbcsr_t_get_num_blocks_total(ri_data%t_3c_int_ctr_1(1, 1))
      IF (nblks == 0) THEN
         CPABORT("3-center integrals are not available (first call requires geometry_did_change=.TRUE.)")
      ENDIF

      DO ispin = 1, nspins
         nblks = dbcsr_t_get_num_blocks_total(ri_data%t_2c_int(ispin, 1))
         IF (nblks == 0) THEN
            CPABORT("2-center integrals are not available (first call requires geometry_did_change=.TRUE.)")
         ENDIF
      ENDDO

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
      ENDIF

      CALL get_qs_env(qs_env, para_env=para_env)

      ALLOCATE (bounds(2, 1))

      CALL dbcsr_get_info(ks_matrix(1, 1)%matrix, distribution=ks_dist)
      CALL dbcsr_distribution_get(ks_dist, group=comm_2d, nprows=pdims_2d(1), npcols=pdims_2d(2))

      pgrid_2d = dbcsr_t_nd_mp_comm(comm_2d, [1], [2], pdims_2d=pdims_2d)

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
         ENDDO

         CALL split_block_sizes(mem_size, mo_bsizes_1, ri_data%min_bsize_MO)
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
                  ENDIF
                  mem_end_block_1(i_mem) = iblock
                  EXIT
               ENDIF
            ENDDO
         ENDDO

         CALL mp_environ(nproc, iproc, para_env%group)

         CALL create_3c_tensor(t_3c_int_mo_1(1, 1), dist1, dist2, dist3, ri_data%pgrid_1, &
                               ri_data%bsizes_AO_split, ri_data%bsizes_RI_split, mo_bsizes_1, &
                               [1, 2], [3], &
                               starts_array_block_3=mem_start_block_1, ends_array_block_3=mem_end_block_1, &
                               name="(AO RI | MO)")

         DEALLOCATE (dist1, dist2, dist3)

         CALL create_3c_tensor(t_3c_int_mo_2(1, 1), dist1, dist2, dist3, ri_data%pgrid_2, &
                               mo_bsizes_1, ri_data%bsizes_RI_split, ri_data%bsizes_AO_split, &
                               [1], [2, 3], &
                               starts_array_block_1=mem_start_block_1, ends_array_block_1=mem_end_block_1, &
                               name="(MO | RI AO)")

         DEALLOCATE (dist1, dist2, dist3)

         CALL create_2c_tensor(mo_coeff_t_split, dist1, dist2, pgrid_2d, ri_data%bsizes_AO_split, mo_bsizes_1, &
                               name="(AO | MO)")

         DEALLOCATE (dist1, dist2)

         CPASSERT(homo(ispin)/ri_data%n_mem > 0)

         IF (do_initialize) THEN
            pdims(:) = 0

            CALL dbcsr_t_pgrid_create(para_env%group, pdims, pgrid, &
                                      tensor_dims=[SIZE(ri_data%bsizes_RI_fit), &
                                                   (homo(ispin) - 1)/ri_data%n_mem + 1, &
                                                   SIZE(ri_data%bsizes_AO_fit)])
            CALL create_3c_tensor(ri_data%t_3c_int_mo(ispin, 1, 1), dist1, dist2, dist3, pgrid, &
                                  ri_data%bsizes_RI_fit, mo_bsizes_2, ri_data%bsizes_AO_fit, &
                                  [1], [2, 3], &
                                  starts_array_block_2=mem_start_block_2, ends_array_block_2=mem_end_block_2, &
                                  name="(RI | MO AO)")

            DEALLOCATE (dist1, dist2, dist3)

            CALL create_3c_tensor(ri_data%t_3c_ctr_KS(ispin, 1, 1), dist1, dist2, dist3, pgrid, &
                                  ri_data%bsizes_RI_fit, mo_bsizes_2, ri_data%bsizes_AO_fit, &
                                  [1, 2], [3], &
                                  starts_array_block_2=mem_start_block_2, ends_array_block_2=mem_end_block_2, &
                                  name="(RI MO | AO)")
            DEALLOCATE (dist1, dist2, dist3)
            CALL dbcsr_t_pgrid_destroy(pgrid)

            CALL dbcsr_t_create(ri_data%t_3c_int_mo(ispin, 1, 1), ri_data%t_3c_ctr_RI(ispin, 1, 1), name="(RI | MO AO)")
            CALL dbcsr_t_create(ri_data%t_3c_ctr_KS(ispin, 1, 1), ri_data%t_3c_ctr_KS_copy(ispin, 1, 1))
         ENDIF

         CALL dbcsr_t_create(mo_coeff(ispin), mo_coeff_t, name="MO coeffs")
         CALL dbcsr_t_copy_matrix_to_tensor(mo_coeff(ispin), mo_coeff_t)
         CALL dbcsr_t_copy(mo_coeff_t, mo_coeff_t_split, move_data=.TRUE.)
         CALL dbcsr_t_filter(mo_coeff_t_split, ri_data%filter_eps_mo)
         CALL dbcsr_t_destroy(mo_coeff_t)

         DO i_mem = 1, n_mem

            IF (i_mem == 2 .AND. do_initialize) THEN
               CPASSERT(ASSOCIATED(pgrid_opt_RI))
               CPASSERT(ASSOCIATED(pgrid_opt_KS))
               CALL dbcsr_t_mp_environ_pgrid(pgrid_opt_RI, pdims_RI, pcoord)
               CALL dbcsr_t_mp_environ_pgrid(pgrid_opt_KS, pdims_AO, pcoord)

               nproc_rem = pdims_RI(2)*pdims_RI(3)
               pdim_AO = pdims_AO(3)
               pdim_MO = MAX(nproc_rem/pdim_AO, 1)

               ! modify pdim_AO until it divides remaining number of processes
               ! and make sure that process dimension corresponding to MOs is much smaller than tensor dimension (for good load balancing)
               pdim_AO_ = pdim_AO
               tdim_mo = (homo(ispin) - 1)/n_mem + 1
               max_pdim_mo = MAX(INT(REAL(tdim_mo, dp)*lb_ratio), 1)
               count = 0

               DO WHILE (MOD(nproc_rem, pdim_AO_) /= 0 .OR. &
                         (pdim_MO >= max_pdim_mo &
                          .AND. MOD(tdim_mo, pdim_MO) /= 0))
                  pdim_AO__ = pdim_AO - (count/2 + 1)*(-1)**count
                  IF (pdim_AO__ > 0) THEN
                     pdim_AO_ = pdim_AO__
                     pdim_MO = MAX(nproc_rem/pdim_AO_, 1)
                  ENDIF
                  count = count + 1
               ENDDO
               pdim_AO = pdim_AO_

               pdims = [pdims_RI(1), nproc_rem/pdim_AO, pdim_AO]

               CALL dbcsr_t_pgrid_change_dims(pgrid_opt_RI, pdims)
               CALL tensor_change_pgrid(ri_data%t_3c_int_mo(ispin, 1, 1), pgrid_opt_RI, nodata=.TRUE., &
                                        starts_array_mc_block_2=mem_start_block_2, ends_array_mc_block_2=mem_end_block_2, &
                                        unit_nr=unit_nr_dbcsr)
               CALL tensor_change_pgrid(ri_data%t_3c_ctr_RI(ispin, 1, 1), pgrid_opt_RI, nodata=.TRUE., &
                                        starts_array_mc_block_2=mem_start_block_2, ends_array_mc_block_2=mem_end_block_2)
               CALL tensor_change_pgrid(ri_data%t_3c_ctr_KS(ispin, 1, 1), pgrid_opt_RI, nodata=.TRUE., &
                                        starts_array_mc_block_2=mem_start_block_2, ends_array_mc_block_2=mem_end_block_2)
               CALL tensor_change_pgrid(ri_data%t_3c_ctr_KS_copy(ispin, 1, 1), pgrid_opt_RI, nodata=.TRUE., &
                                        starts_array_mc_block_2=mem_start_block_2, ends_array_mc_block_2=mem_end_block_2)
               CALL tensor_change_pgrid(t_3c_int_mo_1(1, 1), pgrid_opt_mo_R, nodata=.TRUE., &
                                        starts_array_mc_block_3=mem_start_block_1, ends_array_mc_block_3=mem_end_block_1, &
                                        unit_nr=unit_nr_dbcsr)
               CALL tensor_change_pgrid(ri_data%t_3c_int_ctr_1(1, 1), pgrid_opt_mo_R)
               CALL tensor_change_pgrid(t_3c_int_mo_2(1, 1), pgrid_opt_mo_L, nodata=.TRUE., &
                                        starts_array_mc_block_1=mem_start_block_1, ends_array_mc_block_1=mem_end_block_1, &
                                        unit_nr=unit_nr_dbcsr)
               CALL tensor_change_pgrid(ri_data%t_3c_int_ctr_2(1, 1), pgrid_opt_mo_L)

               ! initialize batched contraction here since process grid have changed
               CALL dbcsr_t_batched_contract_init(ri_data%t_2c_int(ispin, 1))
               CALL dbcsr_t_batched_contract_init(ks_t)

            ELSEIF (i_mem == 2) THEN
               CALL dbcsr_t_batched_contract_init(ks_t)
               IF (geometry_did_change) THEN
                  CALL dbcsr_t_batched_contract_finalize(ri_data%t_2c_int(ispin, 1))
                  CALL dbcsr_t_batched_contract_init(ri_data%t_2c_int(ispin, 1))
               ENDIF
            ENDIF

            bounds(:, 1) = [mem_start(i_mem), mem_end(i_mem)]

            IF (ASSOCIATED(pgrid_opt_mo_R)) THEN
               CALL dbcsr_t_pgrid_destroy(pgrid_opt_mo_R)
               DEALLOCATE (pgrid_opt_mo_R)
            ENDIF

            CALL timeset(routineN//"_MOx3C_R", handle2)
            CALL dbcsr_t_contract(dbcsr_scalar(1.0_dp), mo_coeff_t_split, ri_data%t_3c_int_ctr_1(1, 1), &
                                  dbcsr_scalar(0.0_dp), t_3c_int_mo_1(1, 1), &
                                  contract_1=[1], notcontract_1=[2], &
                                  contract_2=[3], notcontract_2=[1, 2], &
                                  map_1=[3], map_2=[1, 2], &
                                  bounds_2=bounds, &
                                  pgrid_opt_3=pgrid_opt_mo_R, &
                                  filter_eps=ri_data%filter_eps_mo/2, &
                                  unit_nr=unit_nr_dbcsr, &
                                  move_data=.FALSE., &
                                  flop=nflop)

            ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop

            CALL timestop(handle2)
            CALL timeset(routineN//"_copy_1", handle2)
            CALL dbcsr_t_copy(t_3c_int_mo_1(1, 1), ri_data%t_3c_int_mo(ispin, 1, 1), order=[3, 1, 2], move_data=.TRUE.)
            CALL timestop(handle2)

            IF (ASSOCIATED(pgrid_opt_mo_L)) THEN
               CALL dbcsr_t_pgrid_destroy(pgrid_opt_mo_L)
               DEALLOCATE (pgrid_opt_mo_L)
            ENDIF

            CALL timeset(routineN//"_MOx3C_L", handle2)
            CALL dbcsr_t_contract(dbcsr_scalar(1.0_dp), mo_coeff_t_split, ri_data%t_3c_int_ctr_2(1, 1), &
                                  dbcsr_scalar(0.0_dp), t_3c_int_mo_2(1, 1), &
                                  contract_1=[1], notcontract_1=[2], &
                                  contract_2=[1], notcontract_2=[2, 3], &
                                  map_1=[1], map_2=[2, 3], &
                                  bounds_2=bounds, &
                                  pgrid_opt_3=pgrid_opt_mo_L, &
                                  filter_eps=ri_data%filter_eps_mo/2, &
                                  unit_nr=unit_nr_dbcsr, &
                                  move_data=.FALSE., &
                                  flop=nflop)

            ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop

            CALL timestop(handle2)
            CALL timeset(routineN//"_copy_1", handle2)
            CALL dbcsr_t_copy(t_3c_int_mo_2(1, 1), ri_data%t_3c_int_mo(ispin, 1, 1), order=[2, 1, 3], &
                              summation=.TRUE., move_data=.TRUE.)

            CALL dbcsr_t_filter(ri_data%t_3c_int_mo(ispin, 1, 1), ri_data%filter_eps_mo)
            CALL timestop(handle2)

            IF (ASSOCIATED(pgrid_opt_RI)) THEN
               CALL dbcsr_t_pgrid_destroy(pgrid_opt_RI)
               DEALLOCATE (pgrid_opt_RI)
            ENDIF

            CALL timeset(routineN//"_RIx3C", handle2)

            CALL dbcsr_t_contract(dbcsr_scalar(1.0_dp), ri_data%t_2c_int(ispin, 1), ri_data%t_3c_int_mo(ispin, 1, 1), &
                                  dbcsr_scalar(0.0_dp), ri_data%t_3c_ctr_RI(ispin, 1, 1), &
                                  contract_1=[1], notcontract_1=[2], &
                                  contract_2=[1], notcontract_2=[2, 3], &
                                  map_1=[1], map_2=[2, 3], filter_eps=ri_data%filter_eps, &
                                  pgrid_opt_3=pgrid_opt_RI, &
                                  unit_nr=unit_nr_dbcsr, &
                                  flop=nflop)

            ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop

            CALL timestop(handle2)

            CALL timeset(routineN//"_copy_2", handle2)

            ! note: this copy should not involve communication (same block sizes, same 3d distribution on same process grid)
            CALL dbcsr_t_copy(ri_data%t_3c_ctr_RI(ispin, 1, 1), ri_data%t_3c_ctr_KS(ispin, 1, 1), move_data=.TRUE.)
            CALL dbcsr_t_copy(ri_data%t_3c_ctr_KS(ispin, 1, 1), ri_data%t_3c_ctr_KS_copy(ispin, 1, 1))
            CALL timestop(handle2)

            IF (ASSOCIATED(pgrid_opt_KS)) THEN
               CALL dbcsr_t_pgrid_destroy(pgrid_opt_KS)
               DEALLOCATE (pgrid_opt_KS)
            ENDIF

            CALL timeset(routineN//"_3Cx3C", handle2)
            CALL dbcsr_t_contract(dbcsr_scalar(1.0_dp), ri_data%t_3c_ctr_KS(ispin, 1, 1), ri_data%t_3c_ctr_KS_copy(ispin, 1, 1), &
                                  dbcsr_scalar(1.0_dp), ks_t, &
                                  contract_1=[1, 2], notcontract_1=[3], &
                                  contract_2=[1, 2], notcontract_2=[3], &
                                  map_1=[1], map_2=[2], filter_eps=ri_data%filter_eps/n_mem, &
                                  pgrid_opt_1=pgrid_opt_KS, &
                                  unit_nr=unit_nr_dbcsr, move_data=.TRUE., &
                                  flop=nflop)

            ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop

            CALL timestop(handle2)
         ENDDO

         IF (ASSOCIATED(pgrid_opt_mo_R)) THEN
            CALL dbcsr_t_pgrid_destroy(pgrid_opt_mo_R)
            DEALLOCATE (pgrid_opt_mo_R)
         ENDIF

         IF (ASSOCIATED(pgrid_opt_mo_L)) THEN
            CALL dbcsr_t_pgrid_destroy(pgrid_opt_mo_L)
            DEALLOCATE (pgrid_opt_mo_L)
         ENDIF

         IF (ASSOCIATED(pgrid_opt_RI)) THEN
            CALL dbcsr_t_pgrid_destroy(pgrid_opt_RI)
            DEALLOCATE (pgrid_opt_RI)
         ENDIF

         IF (ASSOCIATED(pgrid_opt_KS)) THEN
            CALL dbcsr_t_pgrid_destroy(pgrid_opt_KS)
            DEALLOCATE (pgrid_opt_KS)
         ENDIF

         !CALL dbcsr_t_batched_contract_finalize(ri_data%t_2c_int(1))
         CALL dbcsr_t_batched_contract_finalize(ks_t, unit_nr=unit_nr_dbcsr)

         CALL dbcsr_t_destroy(t_3c_int_mo_1(1, 1))
         CALL dbcsr_t_destroy(t_3c_int_mo_2(1, 1))
         CALL dbcsr_t_clear(ri_data%t_3c_int_mo(ispin, 1, 1))

         CALL dbcsr_t_destroy(mo_coeff_t_split)

         CALL dbcsr_t_filter(ks_t, ri_data%filter_eps)

         CALL dbcsr_t_create(ks_matrix(ispin, 1)%matrix, ks_t_mat)
         CALL dbcsr_t_copy(ks_t, ks_t_mat, move_data=.TRUE.)
         CALL dbcsr_t_copy_tensor_to_matrix(ks_t_mat, ks_matrix(ispin, 1)%matrix, summation=.TRUE.)
         CALL dbcsr_t_destroy(ks_t_mat)

         DEALLOCATE (mem_end, mem_start, mo_bsizes_2, mem_size, mem_start_block_1, mem_end_block_1, &
                     mem_start_block_2, mem_end_block_2)

      ENDDO

      CALL dbcsr_t_pgrid_destroy(pgrid_2d)
      CALL dbcsr_t_destroy(ks_t)

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
! **************************************************************************************************
   SUBROUTINE hfx_ri_update_ks_Pmat(qs_env, ri_data, ks_matrix, rho_ao, &
                                    geometry_did_change, nspins)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data
      TYPE(dbcsr_p_type), DIMENSION(:, :)                :: ks_matrix, rho_ao
      LOGICAL, INTENT(IN)                                :: geometry_did_change
      INTEGER, INTENT(IN)                                :: nspins

      CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_ri_update_ks_Pmat', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, handle2, i_mem, ispin, j_mem, &
                                                            n_mem, unit_nr, unit_nr_dbcsr
      INTEGER(int_8)                                     :: flops_ks_max, flops_p_max, nblks, nflop, &
                                                            nze, nze_3c, nze_3c_1, nze_3c_2, &
                                                            nze_ks, nze_rho
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: dist1, dist2
      INTEGER, DIMENSION(2, 1)                           :: bounds_i
      INTEGER, DIMENSION(2, 2)                           :: bounds_ij, bounds_j
      INTEGER, DIMENSION(2, 3)                           :: bounds_3c
      INTEGER, DIMENSION(3)                              :: dims_3c
      REAL(dp)                                           :: occ, occ_3c, occ_3c_1, occ_3c_2, occ_ks, &
                                                            occ_rho, t1, t2
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_t_pgrid_type), POINTER                  :: pgrid_opt
      TYPE(dbcsr_t_type)                                 :: ks_t, ks_tmp, rho_ao_t, rho_ao_tmp, &
                                                            t_3c_1, t_3c_3

      CALL timeset(routineN, handle)

      NULLIFY (pgrid_opt, para_env)

      ! get a useful output_unit
      unit_nr_dbcsr = ri_data%unit_nr_dbcsr
      unit_nr = ri_data%unit_nr

      CALL get_qs_env(qs_env, para_env=para_env)

      CPASSERT(SIZE(ks_matrix, 2) == 1)

      IF (geometry_did_change) THEN
         CALL hfx_ri_pre_scf_Pmat(qs_env, ri_data)
      ENDIF

      nblks = dbcsr_t_get_num_blocks_total(ri_data%t_3c_int_ctr_2(1, 1))
      IF (nblks == 0) THEN
         CPABORT("3-center integrals are not available (first call requires geometry_did_change=.TRUE.)")
      ENDIF

      n_mem = ri_data%n_mem

      DO i_mem = 1, n_mem
         DO j_mem = 1, n_mem
            IF (ASSOCIATED(ri_data%pgrid_1)) CALL tensor_change_pgrid(ri_data%t_3c_int_ctr_1(i_mem, j_mem), &
                                                                      ri_data%pgrid_1, &
                                                                      starts_array_mc_block_1=ri_data%starts_array_mem_block, &
                                                                      ends_array_mc_block_1=ri_data%ends_array_mem_block, &
                                                                      starts_array_mc_block_2=ri_data%starts_array_RI_mem_block, &
                                                                      ends_array_mc_block_2=ri_data%ends_array_RI_mem_block, &
                                                                      nodata=.FALSE., &
                                                                      blk_indices=ri_data%blk_indices(i_mem, j_mem)%ind, &
                                                                      compressed=ri_data%store_3c(i_mem, j_mem), &
                                                                      eps=ri_data%filter_eps_storage, &
                                                                      unit_nr=unit_nr_dbcsr)
         ENDDO
      ENDDO

      IF (ASSOCIATED(ri_data%pgrid_2)) CALL tensor_change_pgrid(ri_data%t_3c_int_ctr_2(1, 1), &
                                                                ri_data%pgrid_2, &
                                                                starts_array_mc_block_2=ri_data%starts_array_RI_mem_block, &
                                                                ends_array_mc_block_2=ri_data%ends_array_RI_mem_block, &
                                                                starts_array_mc_block_3=ri_data%starts_array_mem_block, &
                                                                ends_array_mc_block_3=ri_data%ends_array_mem_block, &
                                                                nodata=.FALSE., &
                                                                unit_nr=unit_nr_dbcsr)

      CALL dbcsr_t_create(ks_matrix(1, 1)%matrix, ks_tmp)
      CALL dbcsr_t_create(rho_ao(1, 1)%matrix, rho_ao_tmp)

      CALL create_2c_tensor(rho_ao_t, dist1, dist2, ri_data%pgrid_2d, &
                            ri_data%bsizes_AO_split, ri_data%bsizes_AO_split, &
                            name="(AO | AO)")
      DEALLOCATE (dist1, dist2)

      CALL create_2c_tensor(ks_t, dist1, dist2, ri_data%pgrid_2d, &
                            ri_data%bsizes_AO_split, ri_data%bsizes_AO_split, &
                            name="(AO | AO)")
      DEALLOCATE (dist1, dist2)

      CALL dbcsr_t_create(ri_data%t_3c_int_ctr_2(1, 1), t_3c_1)
      CALL dbcsr_t_create(ri_data%t_3c_int_ctr_1(1, 1), t_3c_3)

      CALL mp_sync(para_env%group)
      t1 = m_walltime()

      flops_ks_max = 0; flops_p_max = 0

      DO ispin = 1, nspins

         CALL get_tensor_occupancy(ri_data%t_3c_int_ctr_2(1, 1), nze_3c, occ_3c)

         nze_rho = 0
         occ_rho = 0.0_dp
         nze_3c_1 = 0
         occ_3c_1 = 0.0_dp
         nze_3c_2 = 0
         occ_3c_2 = 0.0_dp

         CALL dbcsr_t_copy_matrix_to_tensor(rho_ao(ispin, 1)%matrix, rho_ao_tmp)
         CALL dbcsr_t_copy(rho_ao_tmp, rho_ao_t, move_data=.TRUE.)

         CALL get_tensor_occupancy(rho_ao_t, nze_rho, occ_rho)

         CALL dbcsr_t_batched_contract_init(ks_t)

         DO i_mem = 1, n_mem

            CALL dbcsr_t_batched_contract_init(rho_ao_t)
            DO j_mem = 1, n_mem

               CALL timeset(routineN//"_Px3C", handle2)

               CALL dbcsr_t_get_info(t_3c_1, nfull_total=dims_3c)
               bounds_i(:, 1) = [ri_data%starts_array_mem(i_mem), ri_data%ends_array_mem(i_mem)]
               bounds_j(:, 1) = [1, dims_3c(1)]
               bounds_j(:, 2) = [ri_data%starts_array_RI_mem(j_mem), ri_data%ends_array_RI_mem(j_mem)]

               CALL dbcsr_t_contract(dbcsr_scalar(1.0_dp), rho_ao_t, ri_data%t_3c_int_ctr_2(1, 1), &
                                     dbcsr_scalar(0.0_dp), t_3c_1, &
                                     contract_1=[2], notcontract_1=[1], &
                                     contract_2=[3], notcontract_2=[1, 2], &
                                     map_1=[3], map_2=[1, 2], filter_eps=ri_data%filter_eps, &
                                     bounds_2=bounds_i, &
                                     bounds_3=bounds_j, &
                                     unit_nr=unit_nr_dbcsr, &
                                     pgrid_opt_2=pgrid_opt, &
                                     flop=nflop)

               CALL timestop(handle2)

               ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop

               IF (nflop > flops_p_max) THEN
                  CPASSERT(ASSOCIATED(pgrid_opt))
                  IF (ASSOCIATED(ri_data%pgrid_2)) THEN
                     CALL dbcsr_t_pgrid_destroy(ri_data%pgrid_2)
                     DEALLOCATE (ri_data%pgrid_2)
                  ENDIF
                  ri_data%pgrid_2 => pgrid_opt
                  NULLIFY (pgrid_opt)
                  flops_p_max = nflop
               ELSEIF (ASSOCIATED(pgrid_opt)) THEN
                  CALL dbcsr_t_pgrid_destroy(pgrid_opt)
                  DEALLOCATE (pgrid_opt)
               ENDIF

               CALL get_tensor_occupancy(t_3c_1, nze, occ)
               nze_3c_1 = nze_3c_1 + nze
               occ_3c_1 = occ_3c_1 + occ

               CALL timeset(routineN//"_copy_2", handle2)
               CALL dbcsr_t_copy(t_3c_1, t_3c_3, order=[3, 2, 1], move_data=.TRUE.)
               CALL timestop(handle2)

               bounds_3c(:, 1) = [ri_data%starts_array_mem(i_mem), ri_data%ends_array_mem(i_mem)]
               bounds_3c(:, 2) = [ri_data%starts_array_RI_mem(j_mem), ri_data%ends_array_RI_mem(j_mem)]
               bounds_3c(:, 3) = [1, dims_3c(3)]

               bounds_ij(:, 1) = [ri_data%starts_array_mem(i_mem), ri_data%ends_array_mem(i_mem)]
               bounds_ij(:, 2) = [ri_data%starts_array_RI_mem(j_mem), ri_data%ends_array_RI_mem(j_mem)]

               CALL decompress_tensor(ri_data%t_3c_int_ctr_1(i_mem, j_mem), ri_data%blk_indices(i_mem, j_mem)%ind, &
                                      ri_data%store_3c(i_mem, j_mem), ri_data%filter_eps_storage)

               CALL get_tensor_occupancy(ri_data%t_3c_int_ctr_1(i_mem, j_mem), nze, occ)
               nze_3c_2 = nze_3c_2 + nze
               occ_3c_2 = occ_3c_2 + occ
               CALL timeset(routineN//"_KS", handle2)
               CALL dbcsr_t_contract(dbcsr_scalar(1.0_dp), ri_data%t_3c_int_ctr_1(i_mem, j_mem), t_3c_3, &
                                     dbcsr_scalar(1.0_dp), ks_t, &
                                     contract_1=[1, 2], notcontract_1=[3], &
                                     contract_2=[1, 2], notcontract_2=[3], &
                                     map_1=[1], map_2=[2], filter_eps=ri_data%filter_eps/n_mem, &
                                     bounds_1=bounds_ij, &
                                     unit_nr=unit_nr_dbcsr, &
                                     pgrid_opt_1=pgrid_opt, &
                                     flop=nflop, move_data=.TRUE.)

               CALL timestop(handle2)

               ri_data%dbcsr_nflop = ri_data%dbcsr_nflop + nflop

               IF (nflop > flops_ks_max) THEN
                  CPASSERT(ASSOCIATED(pgrid_opt))
                  IF (ASSOCIATED(ri_data%pgrid_1)) THEN
                     CALL dbcsr_t_pgrid_destroy(ri_data%pgrid_1)
                     DEALLOCATE (ri_data%pgrid_1)
                  ENDIF
                  ri_data%pgrid_1 => pgrid_opt
                  NULLIFY (pgrid_opt)
                  flops_ks_max = nflop
               ELSEIF (ASSOCIATED(pgrid_opt)) THEN
                  CALL dbcsr_t_pgrid_destroy(pgrid_opt)
                  DEALLOCATE (pgrid_opt)
               ENDIF

            ENDDO
            CALL dbcsr_t_batched_contract_finalize(rho_ao_t, unit_nr=unit_nr_dbcsr)
         ENDDO
         CALL dbcsr_t_batched_contract_finalize(ks_t, unit_nr=unit_nr_dbcsr)
         CALL dbcsr_t_clear(rho_ao_t)
         CALL get_tensor_occupancy(ks_t, nze_ks, occ_ks)

         CALL dbcsr_t_copy(ks_t, ks_tmp)
         CALL dbcsr_t_clear(ks_t)
         CALL dbcsr_t_copy_tensor_to_matrix(ks_tmp, ks_matrix(ispin, 1)%matrix)
         CALL dbcsr_t_clear(ks_tmp)

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
         ENDIF

      ENDDO

      CALL mp_sync(para_env%group)
      t2 = m_walltime()

      ri_data%dbcsr_time = ri_data%dbcsr_time + t2 - t1

      CALL dbcsr_t_destroy(t_3c_1)
      CALL dbcsr_t_destroy(t_3c_3)

      CALL dbcsr_t_destroy(rho_ao_t)
      CALL dbcsr_t_destroy(rho_ao_tmp)
      CALL dbcsr_t_destroy(ks_t)
      CALL dbcsr_t_destroy(ks_tmp)

      CALL timestop(handle)

   END SUBROUTINE

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
