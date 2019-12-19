!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief RI-methods for HFX
! **************************************************************************************************

MODULE hfx_ri

   USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE cp_array_utils,                  ONLY: cp_1d_r_p_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_cholesky,               ONLY: cp_dbcsr_cholesky_decompose,&
                                              cp_dbcsr_cholesky_invert
   USE cp_dbcsr_func,                   ONLY: cp_dbcsr_matrix_function
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
        dbcsr_frobenius_norm, dbcsr_get_info, dbcsr_multiply, dbcsr_p_type, dbcsr_release, &
        dbcsr_scalar, dbcsr_scale, dbcsr_type, dbcsr_type_no_symmetry, dbcsr_type_real_8, &
        dbcsr_type_symmetric
   USE dbcsr_tensor_api,                ONLY: &
        dbcsr_t_clear, dbcsr_t_contract, dbcsr_t_copy, dbcsr_t_copy_matrix_to_tensor, &
        dbcsr_t_copy_tensor_to_matrix, dbcsr_t_create, dbcsr_t_destroy, dbcsr_t_filter, &
        dbcsr_t_nd_mp_comm, dbcsr_t_pgrid_create, dbcsr_t_pgrid_destroy, dbcsr_t_pgrid_type, &
        dbcsr_t_type
   USE distribution_2d_types,           ONLY: distribution_2d_type
   USE hfx_types,                       ONLY: hfx_ri_type
   USE input_constants,                 ONLY: hfx_ri_do_2c_diag,&
                                              hfx_ri_do_2c_iter
   USE input_cp2k_hfx,                  ONLY: ri_mo,&
                                              ri_pmat
   USE iterate_matrix,                  ONLY: invert_hotelling,&
                                              matrix_sqrt_newton_schulz
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE message_passing,                 ONLY: mp_environ
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
   USE qs_tensors,                      ONLY: build_2c_integrals,&
                                              build_2c_neighbor_lists,&
                                              build_3c_integrals,&
                                              build_3c_neighbor_lists,&
                                              neighbor_list_3c_destroy,&
                                              tensor_change_pgrid
   USE qs_tensors_types,                ONLY: create_2c_tensor,&
                                              create_3c_tensor,&
                                              neighbor_list_3c_type,&
                                              split_block_sizes
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
! **************************************************************************************************
   SUBROUTINE hfx_ri_pre_scf_mo(qs_env, ri_data)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(hfx_ri_type), INTENT(INOUT)                   :: ri_data

      CHARACTER(LEN=*), PARAMETER :: routineN = 'hfx_ri_pre_scf_mo', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, handle2, unit_nr, unit_nr_dbcsr
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
         CASE (hfx_ri_do_2c_diag)
            CALL dbcsr_copy(t_2c_op_RI_inv(1), t_2c_op_RI(1))
            CALL cp_dbcsr_cholesky_decompose(t_2c_op_RI_inv(1), para_env=para_env, blacs_env=blacs_env)
            CALL cp_dbcsr_cholesky_invert(t_2c_op_RI_inv(1), para_env=para_env, blacs_env=blacs_env, upper_to_full=.TRUE.)
            CALL dbcsr_filter(t_2c_op_RI_inv(1), ri_data%filter_eps)
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
         CASE (hfx_ri_do_2c_diag)
            CALL dbcsr_create(t_2c_op_pot_sqrt(1), template=t_2c_op_pot(1), matrix_type=dbcsr_type_symmetric)

            CALL cp_dbcsr_matrix_function(t_2c_op_pot(1), my_sqrt, t_2c_op_pot_sqrt(1), para_env=para_env, blacs_env=blacs_env)
            CALL dbcsr_filter(t_2c_op_pot(1), ri_data%filter_eps)
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
         CASE (hfx_ri_do_2c_diag)
            CALL dbcsr_create(t_2c_int_mat(1), template=t_2c_op_pot(1), matrix_type=dbcsr_type_symmetric)
            CALL cp_dbcsr_matrix_function(t_2c_op_pot(1), my_invsqrt, t_2c_int_mat(1), para_env=para_env, blacs_env=blacs_env)
            CALL dbcsr_filter(t_2c_op_pot(1), ri_data%filter_eps)
         END SELECT
         IF (ri_data%check_2c_inv) THEN
            CALL check_sqrt(t_2c_op_pot(1), matrix_sqrt_inv=t_2c_int_mat(1), unit_nr=unit_nr)
         ENDIF
      ENDIF

      CALL dbcsr_release(t_2c_op_pot(1))

      CALL dbcsr_t_create(t_2c_int_mat(1), t_2c_int(1), name="(RI|RI)")
      CALL dbcsr_t_copy_matrix_to_tensor(t_2c_int_mat(1), t_2c_int(1))
      CALL dbcsr_release(t_2c_int_mat(1))
      CALL dbcsr_t_copy(t_2c_int(1), ri_data%t_2c_int(1), move_data=.TRUE.)
      CALL dbcsr_t_destroy(t_2c_int(1))
      CALL timestop(handle2)

      CALL timeset(routineN//"_3c", handle2)
      CALL dbcsr_t_copy(ri_data%t_3c_int(1, 1), ri_data%t_3c_int_ctr_1(1, 1), order=[2, 1, 3], move_data=.TRUE.)
      CALL dbcsr_t_filter(ri_data%t_3c_int_ctr_1(1, 1), ri_data%filter_eps)
      CALL dbcsr_t_destroy(ri_data%t_3c_int(1, 1))
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

      INTEGER                                            :: handle, ibasis, natom, nkind
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: sizes_AO, sizes_RI
      INTEGER, DIMENSION(:), POINTER                     :: col_bsize, row_bsize
      TYPE(dbcsr_distribution_type)                      :: dbcsr_dist
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(distribution_2d_type), POINTER                :: dist_2d
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

      CALL dbcsr_t_create(t_3c_int(1, 1), "(RI | AO AO)", ri_data%dist, [1], [2, 3], &
                          dbcsr_type_real_8, sizes_RI, sizes_AO, sizes_AO)

      CALL build_3c_neighbor_lists(nl_3c, basis_set_RI, basis_set_AO, basis_set_AO, ri_data%dist_3d, ri_data%ri_metric, &
                                   "HFX_3c_nl", &
                                   qs_env, op_pos=1, sym_jk=.TRUE.)

      CALL build_3c_integrals(t_3c_int, ri_data%filter_eps_3c, qs_env, nl_3c, basis_set_RI, basis_set_AO, basis_set_AO, &
                              ri_data%ri_metric, int_eps=ri_data%eps_schwarz, op_pos=1, &
                              desymmetrize=ri_data%flavor == ri_pmat)

      CALL neighbor_list_3c_destroy(nl_3c)

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

      INTEGER                                            :: handle, handle2, unit_nr, unit_nr_dbcsr
      REAL(KIND=dp)                                      :: threshold
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_t_pgrid_type), POINTER                  :: pgrid_opt
      TYPE(dbcsr_t_type), DIMENSION(1)                   :: t_2c_int
      TYPE(dbcsr_t_type), DIMENSION(1, 1)                :: t_3c_int_1, t_3c_int_2
      TYPE(dbcsr_type), DIMENSION(1)                     :: t_2c_int_mat, t_2c_op_pot, t_2c_op_RI, &
                                                            t_2c_tmp, t_2c_tmp_2

      CALL timeset(routineN, handle)

      NULLIFY (pgrid_opt)

      unit_nr_dbcsr = ri_data%unit_nr_dbcsr
      unit_nr = ri_data%unit_nr

      CALL get_qs_env(qs_env, para_env=para_env, blacs_env=blacs_env)

      CALL timeset(routineN//"_int", handle2)
      CALL hfx_ri_pre_scf_calc_tensors(qs_env, ri_data, t_2c_op_RI, t_2c_op_pot, t_3c_int_1)

      IF (ri_data%same_op) t_2c_op_RI(1) = t_2c_op_pot(1)

      IF (ASSOCIATED(ri_data%pgrid_3)) CALL tensor_change_pgrid(ri_data%t_3c_int_ctr_3(1, 1), &
                                                                ri_data%pgrid_3, nodata=.TRUE., unit_nr=unit_nr_dbcsr)
      IF (ASSOCIATED(ri_data%pgrid_2)) CALL tensor_change_pgrid(ri_data%t_3c_int_ctr_2(1, 1), &
                                                                ri_data%pgrid_2, nodata=.TRUE., unit_nr=unit_nr_dbcsr)

      CALL dbcsr_t_copy(t_3c_int_1(1, 1), ri_data%t_3c_int_ctr_3(1, 1))
      CALL dbcsr_t_filter(ri_data%t_3c_int_ctr_3(1, 1), ri_data%filter_eps)
      CALL dbcsr_t_copy(t_3c_int_1(1, 1), ri_data%t_3c_int_ctr_2(1, 1), order=[2, 1, 3], move_data=.TRUE.)
      CALL dbcsr_t_destroy(t_3c_int_1(1, 1))
      CALL dbcsr_t_filter(ri_data%t_3c_int_ctr_2(1, 1), ri_data%filter_eps)

      CALL timestop(handle2)

      CALL timeset(routineN//"_2c", handle2)
      CALL dbcsr_create(t_2c_int_mat(1), template=t_2c_op_RI(1), matrix_type=dbcsr_type_no_symmetry)
      threshold = MAX(ri_data%filter_eps, 1.0e-12_dp)

      SELECT CASE (ri_data%t2c_method)
      CASE (hfx_ri_do_2c_iter)
         CALL invert_hotelling(t_2c_int_mat(1), t_2c_op_RI(1), &
                               threshold=threshold, silent=.FALSE.)
      CASE (hfx_ri_do_2c_diag)
         CALL dbcsr_copy(t_2c_int_mat(1), t_2c_op_RI(1))
         CALL cp_dbcsr_cholesky_decompose(t_2c_int_mat(1), para_env=para_env, blacs_env=blacs_env)
         CALL cp_dbcsr_cholesky_invert(t_2c_int_mat(1), para_env=para_env, blacs_env=blacs_env, upper_to_full=.TRUE.)
         CALL dbcsr_filter(t_2c_int_mat(1), ri_data%filter_eps)
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

      CALL timestop(handle2)
      CALL timeset(routineN//"_3c", handle2)

      CALL dbcsr_t_create(ri_data%t_3c_int_ctr_3(1, 1), t_3c_int_2(1, 1), name="(RI | AO AO)")
      CALL dbcsr_t_contract(dbcsr_scalar(1.0_dp), t_2c_int(1), ri_data%t_3c_int_ctr_3(1, 1), &
                            dbcsr_scalar(0.0_dp), t_3c_int_2(1, 1), &
                            contract_1=[2], notcontract_1=[1], &
                            contract_2=[1], notcontract_2=[2, 3], &
                            map_1=[1], map_2=[2, 3], filter_eps=ri_data%filter_eps, &
                            unit_nr=unit_nr_dbcsr, move_data=.TRUE., &
                            pgrid_opt_2=pgrid_opt)
      CALL dbcsr_t_clear(ri_data%t_3c_int_ctr_3(1, 1))
      CALL dbcsr_t_destroy(t_2c_int(1))

      CPASSERT(ASSOCIATED(pgrid_opt))
      IF (ASSOCIATED(ri_data%pgrid_3)) THEN
         CALL dbcsr_t_pgrid_destroy(ri_data%pgrid_3)
         DEALLOCATE (ri_data%pgrid_3)
      ENDIF
      ri_data%pgrid_3 => pgrid_opt

      IF (ASSOCIATED(ri_data%pgrid_1)) CALL tensor_change_pgrid(ri_data%t_3c_int_ctr_1(1, 1), ri_data%pgrid_1, &
                                                                nodata=.TRUE., unit_nr=unit_nr_dbcsr)
      CALL dbcsr_t_copy(t_3c_int_2(1, 1), ri_data%t_3c_int_ctr_1(1, 1), order=[2, 1, 3], move_data=.TRUE.)
      CALL dbcsr_t_destroy(t_3c_int_2(1, 1))
      CALL dbcsr_t_filter(ri_data%t_3c_int_ctr_1(1, 1), ri_data%filter_eps)

      CALL timestop(handle2)

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Calculate Fock (AKA Kohn-Sham) matrix and calculate Hartree-Fock exchange energy based on RI expansion.
!> \param qs_env ...
!> \param ri_data RI parameters and intermediate tensor data
!> \param ks_matrix Fock matrix
!> \param ehfx exchange energy
!> \param mos MO coefficients
!> \param rho_ao Density matrix
!> \param geometry_did_change flag that indicates we have to recalc integrals
!> \param nspins Number of spins
!> \param hf_fraction Fraction of exact exchange
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
      INTEGER, DIMENSION(2)                              :: homo
      REAL(dp)                                           :: etmp, fac
      REAL(KIND=dp), DIMENSION(:), POINTER               :: mo_eigenvalues
      TYPE(cp_1d_r_p_type), DIMENSION(:), POINTER        :: occupied_evals
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: homo_localized, moloc_coeff, &
                                                            occupied_orbs
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
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

      INTEGER                                            :: comm_2d, handle, handle2, iproc, ispin, &
                                                            n_mos, nproc, pdim_mo, unit_nr_dbcsr
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: dist1, dist2, dist3, mo_bsizes
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: bounds
      INTEGER, DIMENSION(2)                              :: pdims_2d
      INTEGER, DIMENSION(3)                              :: pdims
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_distribution_type)                      :: ks_dist
      TYPE(dbcsr_t_pgrid_type)                           :: pgrid_1, pgrid_2, pgrid_2d
      TYPE(dbcsr_t_type)                                 :: ks_t, mo_coeff_t, mo_coeff_t_split
      TYPE(dbcsr_t_type), DIMENSION(1, 1)                :: t_3c_ctr_KS, t_3c_ctr_KS_copy, &
                                                            t_3c_ctr_RI, t_3c_int_mo_1, &
                                                            t_3c_int_mo_2, t_3c_int_mo_3

      CALL timeset(routineN, handle)

      CPASSERT(SIZE(ks_matrix, 2) == 1)

      unit_nr_dbcsr = ri_data%unit_nr_dbcsr

      IF (geometry_did_change) THEN
         CALL hfx_ri_pre_scf_mo(qs_env, ri_data)
      ENDIF

      CALL get_qs_env(qs_env, para_env=para_env)

      ALLOCATE (bounds(2, 1))
      DO ispin = 1, nspins
         CALL dbcsr_get_info(mo_coeff(ispin), nfullcols_total=n_mos)
         CALL split_block_sizes([n_mos], mo_bsizes, ri_data%max_bsize_mo)

         CALL mp_environ(nproc, iproc, para_env%group)

         ! restrict process grid dimension corresponding to MO index to be smaller than the number of occupied MO blocks
         pdim_mo = MIN((homo(ispin) - 1)/ri_data%max_bsize_mo + 1, nproc)
         DO WHILE (MOD(nproc, pdim_mo) .NE. 0)
            pdim_mo = pdim_mo - 1
         ENDDO

         pdims = [nproc/pdim_mo, 1, pdim_mo]
         CALL dbcsr_t_pgrid_create(para_env%group, pdims, pgrid_1)


         pdims = [nproc/pdim_mo, pdim_mo, 1]
         CALL dbcsr_t_pgrid_create(para_env%group, pdims, pgrid_2)

         CALL create_3c_tensor(t_3c_int_mo_1(1, 1), dist1, dist2, dist3, ri_data%pgrid_1, &
                               ri_data%bsizes_AO, ri_data%bsizes_RI_split, mo_bsizes, &
                               [1, 2], [3], name="(AO RI | MO)")

         DEALLOCATE (dist1, dist2, dist3)

         CALL create_3c_tensor(t_3c_int_mo_2(1, 1), dist1, dist2, dist3, ri_data%pgrid_2, &
                               mo_bsizes, ri_data%bsizes_RI_split, ri_data%bsizes_AO, &
                               [1], [2, 3], "(MO | RI AO)")

         DEALLOCATE (dist1, dist2, dist3)

         CALL create_3c_tensor(t_3c_int_mo_3(1, 1), dist1, dist2, dist3, pgrid_2, &
                               ri_data%bsizes_RI_split, mo_bsizes, ri_data%bsizes_AO, &
                               [1], [2, 3], name="(RI | MO AO)")

         DEALLOCATE (dist1, dist2, dist3)
         CALL create_3c_tensor(t_3c_ctr_KS(1, 1), dist1, dist2, dist3, pgrid_2, &
                               ri_data%bsizes_RI_split, mo_bsizes, ri_data%bsizes_AO, &
                               [1, 2], [3], name="(RI MO | AO)")
         DEALLOCATE (dist1, dist2, dist3)
         CALL dbcsr_t_pgrid_destroy(pgrid_2)

         CALL dbcsr_t_create(t_3c_ctr_KS(1, 1), t_3c_ctr_KS_copy(1, 1))

         CALL dbcsr_t_create(ks_matrix(1, 1)%matrix, ks_t)

         ! AO block sizes are split in 3c tensor, we need to split blocks of Fock matrix as well
         ! also split MO coefficients
         CALL dbcsr_get_info(ks_matrix(1, 1)%matrix, distribution=ks_dist)
         CALL dbcsr_distribution_get(ks_dist, group=comm_2d, nprows=pdims_2d(1), npcols=pdims_2d(2))
         pgrid_2d = dbcsr_t_nd_mp_comm(comm_2d, [1], [2], pdims_2d=pdims_2d)

         CALL create_2c_tensor(mo_coeff_t_split, dist1, dist2, pgrid_2d, ri_data%bsizes_AO, mo_bsizes, &
                               name="(AO | MO)")

         DEALLOCATE (dist1, dist2)
         CALL dbcsr_t_pgrid_destroy(pgrid_2d)

         CALL dbcsr_t_create(mo_coeff(ispin), mo_coeff_t, name="MO coeffs")
         CALL dbcsr_t_copy_matrix_to_tensor(mo_coeff(ispin), mo_coeff_t)
         CALL dbcsr_t_copy(mo_coeff_t, mo_coeff_t_split, move_data=.TRUE.)
         CALL dbcsr_t_filter(mo_coeff_t_split, ri_data%filter_eps_mo)
         CALL dbcsr_t_destroy(mo_coeff_t)

         bounds(1, 1) = 1
         bounds(2, 1) = homo(ispin)

         CALL timeset(routineN//"_MOx3C_R", handle2)
         CALL dbcsr_t_contract(dbcsr_scalar(1.0_dp), mo_coeff_t_split, ri_data%t_3c_int_ctr_1(1, 1), &
                               dbcsr_scalar(0.0_dp), t_3c_int_mo_1(1, 1), &
                               contract_1=[1], notcontract_1=[2], &
                               contract_2=[3], notcontract_2=[1, 2], &
                               map_1=[3], map_2=[1, 2], &
                               bounds_2=bounds, &
                               filter_eps=ri_data%filter_eps_mo/2, &
                               unit_nr=unit_nr_dbcsr, &
                               move_data=.FALSE.)

         CALL timestop(handle2)
         CALL timeset(routineN//"_copy_1", handle2)
         CALL dbcsr_t_copy(t_3c_int_mo_1(1, 1), t_3c_int_mo_3(1, 1), order=[3, 1, 2], move_data=.TRUE.)
         CALL dbcsr_t_destroy(t_3c_int_mo_1(1, 1))
         CALL timestop(handle2)

         CALL timeset(routineN//"_copy_symm", handle2)
         CALL dbcsr_t_copy(ri_data%t_3c_int_ctr_1(1, 1), ri_data%t_3c_int_ctr_2(1, 1), move_data=.TRUE.)
         CALL dbcsr_t_clear(ri_data%t_3c_int_ctr_1(1, 1))
         CALL timestop(handle2)

         CALL timeset(routineN//"_MOx3C_L", handle2)
         CALL dbcsr_t_contract(dbcsr_scalar(1.0_dp), mo_coeff_t_split, ri_data%t_3c_int_ctr_2(1, 1), &
                               dbcsr_scalar(0.0_dp), t_3c_int_mo_2(1, 1), &
                               contract_1=[1], notcontract_1=[2], &
                               contract_2=[1], notcontract_2=[2, 3], &
                               map_1=[1], map_2=[2, 3], &
                               bounds_2=bounds, &
                               filter_eps=ri_data%filter_eps_mo/2, &
                               unit_nr=unit_nr_dbcsr, &
                               move_data=.FALSE.)
         CALL timestop(handle2)
         CALL timeset(routineN//"_copy_1", handle2)
         CALL dbcsr_t_copy(t_3c_int_mo_2(1, 1), t_3c_int_mo_3(1, 1), order=[2, 1, 3], summation=.TRUE., move_data=.TRUE.)
         CALL dbcsr_t_filter(t_3c_int_mo_3(1, 1), ri_data%filter_eps_mo)
         CALL dbcsr_t_destroy(t_3c_int_mo_2(1, 1))
         CALL timestop(handle2)

         CALL timeset(routineN//"_copy_symm", handle2)
         CALL dbcsr_t_copy(ri_data%t_3c_int_ctr_2(1, 1), ri_data%t_3c_int_ctr_1(1, 1), move_data=.TRUE.)
         CALL dbcsr_t_clear(ri_data%t_3c_int_ctr_2(1, 1))
         CALL timestop(handle2)

         CALL dbcsr_t_destroy(mo_coeff_t_split)

         CALL timeset(routineN//"_RIx3C", handle2)
         CALL dbcsr_t_create(t_3c_int_mo_3(1, 1), t_3c_ctr_RI(1, 1), name="(RI | MO AO)")
         CALL dbcsr_t_contract(dbcsr_scalar(1.0_dp), ri_data%t_2c_int(1), t_3c_int_mo_3(1, 1), &
                               dbcsr_scalar(0.0_dp), t_3c_ctr_RI(1, 1), &
                               contract_1=[1], notcontract_1=[2], &
                               contract_2=[1], notcontract_2=[2, 3], &
                               map_1=[1], map_2=[2, 3], filter_eps=ri_data%filter_eps, &
                               unit_nr=unit_nr_dbcsr)
         CALL dbcsr_t_destroy(t_3c_int_mo_3(1, 1))
         CALL timestop(handle2)

         CALL timeset(routineN//"_copy_2", handle2)

         ! note: this copy should not involve communication (same block sizes, same 3d distribution on same process grid)
         CALL dbcsr_t_copy(t_3c_ctr_RI(1, 1), t_3c_ctr_KS(1, 1), move_data=.TRUE.)
         CALL dbcsr_t_destroy(t_3c_ctr_RI(1, 1))
         CALL dbcsr_t_copy(t_3c_ctr_KS(1, 1), t_3c_ctr_KS_copy(1, 1))
         CALL timestop(handle2)

         CALL timeset(routineN//"_3Cx3C", handle2)
         CALL dbcsr_t_contract(dbcsr_scalar(1.0_dp), t_3c_ctr_KS(1, 1), t_3c_ctr_KS_copy(1, 1), &
                               dbcsr_scalar(0.0_dp), ks_t, &
                               contract_1=[1, 2], notcontract_1=[3], &
                               contract_2=[1, 2], notcontract_2=[3], &
                               map_1=[1], map_2=[2], filter_eps=ri_data%filter_eps, &
                               unit_nr=unit_nr_dbcsr, move_data=.TRUE.)
         CALL timestop(handle2)

         CALL dbcsr_t_destroy(t_3c_ctr_KS(1, 1))
         CALL dbcsr_t_destroy(t_3c_ctr_KS_copy(1, 1))

         CALL dbcsr_t_filter(ks_t, ri_data%filter_eps)
         CALL dbcsr_t_copy_tensor_to_matrix(ks_t, ks_matrix(ispin, 1)%matrix, summation=.TRUE.)
         CALL dbcsr_t_destroy(ks_t)

      ENDDO

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

      INTEGER                                            :: handle, handle2, ispin, unit_nr_dbcsr
      TYPE(dbcsr_t_pgrid_type), POINTER                  :: pgrid_opt
      TYPE(dbcsr_t_type)                                 :: ks_t, rho_ao_t
      TYPE(dbcsr_t_type), DIMENSION(1, 1)                :: t_3c_1, t_3c_2

      CALL timeset(routineN, handle)

      NULLIFY (pgrid_opt)

      ! get a useful output_unit

      unit_nr_dbcsr = ri_data%unit_nr_dbcsr

      CPASSERT(SIZE(ks_matrix, 2) == 1)

      IF (geometry_did_change) THEN
         CALL hfx_ri_pre_scf_Pmat(qs_env, ri_data)
      ELSE
         IF (ASSOCIATED(ri_data%pgrid_2)) CALL tensor_change_pgrid(ri_data%t_3c_int_ctr_2(1, 1), &
                                                                   ri_data%pgrid_2, unit_nr=unit_nr_dbcsr)
         IF (ASSOCIATED(ri_data%pgrid_1)) CALL tensor_change_pgrid(ri_data%t_3c_int_ctr_1(1, 1), &
                                                                   ri_data%pgrid_1, unit_nr=unit_nr_dbcsr)
      ENDIF

      CALL dbcsr_t_create(ri_data%t_3c_int_ctr_2(1, 1), t_3c_1(1, 1))

      CALL dbcsr_t_create(ks_matrix(1, 1)%matrix, ks_t)
      CALL dbcsr_t_create(rho_ao(1, 1)%matrix, rho_ao_t)

      DO ispin = 1, nspins
         CALL dbcsr_t_copy_matrix_to_tensor(rho_ao(ispin, 1)%matrix, rho_ao_t)

         CALL timeset(routineN//"_Px3C", handle2)
         CALL dbcsr_t_contract(dbcsr_scalar(1.0_dp), rho_ao_t, ri_data%t_3c_int_ctr_2(1, 1), &
                               dbcsr_scalar(0.0_dp), t_3c_1(1, 1), &
                               contract_1=[2], notcontract_1=[1], &
                               contract_2=[1], notcontract_2=[2, 3], &
                               map_1=[1], map_2=[2, 3], filter_eps=ri_data%filter_eps, &
                               unit_nr=unit_nr_dbcsr, &
                               pgrid_opt_2=pgrid_opt)

         CALL timestop(handle2)

         CALL timeset(routineN//"_copy", handle2)
         CALL dbcsr_t_create(ri_data%t_3c_int_ctr_1(1, 1), t_3c_2(1, 1))
         CALL dbcsr_t_copy(t_3c_1(1, 1), t_3c_2(1, 1), move_data=.TRUE.)
         CALL dbcsr_t_clear(t_3c_1(1, 1))
         CALL dbcsr_t_filter(t_3c_2(1, 1), ri_data%filter_eps)
         CALL timestop(handle2)

         CPASSERT(ASSOCIATED(pgrid_opt))
         IF (ASSOCIATED(ri_data%pgrid_2)) THEN
            CALL dbcsr_t_pgrid_destroy(ri_data%pgrid_2)
            DEALLOCATE (ri_data%pgrid_2)
         ENDIF
         ri_data%pgrid_2 => pgrid_opt

         NULLIFY (pgrid_opt)

         CALL timeset(routineN//"_KS", handle2)
         CALL dbcsr_t_contract(dbcsr_scalar(1.0_dp), ri_data%t_3c_int_ctr_1(1, 1), t_3c_2(1, 1), &
                               dbcsr_scalar(0.0_dp), ks_t, &
                               contract_1=[1, 2], notcontract_1=[3], &
                               contract_2=[1, 2], notcontract_2=[3], &
                               map_1=[1], map_2=[2], filter_eps=ri_data%filter_eps, &
                               unit_nr=unit_nr_dbcsr, &
                               pgrid_opt_1=pgrid_opt)
         CALL timestop(handle2)

         CALL dbcsr_t_destroy(t_3c_2(1, 1))

         CALL dbcsr_t_copy_tensor_to_matrix(ks_t, ks_matrix(ispin, 1)%matrix, summation=.TRUE.)
         CALL dbcsr_t_clear(ks_t)

         CPASSERT(ASSOCIATED(pgrid_opt))
         IF (ASSOCIATED(ri_data%pgrid_1)) THEN
            CALL dbcsr_t_pgrid_destroy(ri_data%pgrid_1)
            DEALLOCATE (ri_data%pgrid_1)
         ENDIF
         ri_data%pgrid_1 => pgrid_opt

      ENDDO

      CALL dbcsr_t_destroy(t_3c_1(1, 1))
      CALL dbcsr_t_destroy(rho_ao_t)
      CALL dbcsr_t_destroy(ks_t)

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
