!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for kpoint treatment in GW
!> \par History
!>      04.2019 created [Jan Wilhelm]
! **************************************************************************************************
MODULE rpa_gw_kpoints
   USE basis_set_types,                 ONLY: gto_basis_set_p_type
   USE cell_types,                      ONLY: cell_type,&
                                              get_cell,&
                                              pbc
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_cholesky_invert,&
                                              cp_cfm_gemm,&
                                              cp_cfm_scale_and_add,&
                                              cp_cfm_scale_and_add_fm,&
                                              cp_cfm_transpose
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_get_info,&
                                              cp_cfm_p_type,&
                                              cp_cfm_release,&
                                              cp_cfm_set_all,&
                                              cp_cfm_type
   USE cp_dbcsr_operations,             ONLY: copy_fm_to_dbcsr
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add
   USE cp_fm_types,                     ONLY: cp_fm_copy_general,&
                                              cp_fm_create,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_create, dbcsr_dot, dbcsr_get_block_p, dbcsr_get_num_blocks, &
        dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, dbcsr_iterator_start, &
        dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_p_type, dbcsr_release, dbcsr_set, &
        dbcsr_type
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE mathconstants,                   ONLY: gaussi,&
                                              twopi,&
                                              z_one,&
                                              z_zero
   USE mathlib,                         ONLY: invmat
   USE message_passing,                 ONLY: mp_sum
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_integral_utils,               ONLY: basis_set_list_setup
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE rpa_im_time,                     ONLY: compute_transl_dm
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rpa_gw_kpoints'

   PUBLIC :: compute_Wc_real_space_tau_GW, compute_Wc_kp_tau_GW, &
             compute_wkp_W

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param fm_mat_W_tau ...
!> \param cfm_mat_Q ...
!> \param fm_mat_L_re ...
!> \param fm_mat_L_im ...
!> \param dimen_RI ...
!> \param num_integ_points ...
!> \param jquad ...
!> \param ikp ...
!> \param tj ...
!> \param tau_tj ...
!> \param weights_cos_tf_w_to_t ...
!> \param ikp_local ...
!> \param para_env ...
!> \param kpoints ...
!> \param qs_env ...
!> \param wkp_W ...
!> \param mat_SinvVSinv ...
!> \param do_W_and_not_V ...
! **************************************************************************************************
   SUBROUTINE compute_Wc_real_space_tau_GW(fm_mat_W_tau, cfm_mat_Q, fm_mat_L_re, fm_mat_L_im, &
                                           dimen_RI, num_integ_points, jquad, &
                                           ikp, tj, tau_tj, weights_cos_tf_w_to_t, ikp_local, &
                                           para_env, kpoints, qs_env, wkp_W, mat_SinvVSinv, do_W_and_not_V)

      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: fm_mat_W_tau
      TYPE(cp_cfm_type), POINTER                         :: cfm_mat_Q
      TYPE(cp_fm_type), POINTER                          :: fm_mat_L_re, fm_mat_L_im
      INTEGER                                            :: dimen_RI, num_integ_points, jquad, ikp
      REAL(KIND=dp), DIMENSION(:)                        :: tj
      REAL(KIND=dp), DIMENSION(0:num_integ_points)       :: tau_tj
      REAL(KIND=dp), DIMENSION(:, :)                     :: weights_cos_tf_w_to_t
      INTEGER, DIMENSION(:)                              :: ikp_local
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), DIMENSION(:), POINTER               :: wkp_W
      TYPE(dbcsr_p_type)                                 :: mat_SinvVSinv
      LOGICAL                                            :: do_W_and_not_V

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_Wc_real_space_tau_GW', &
         routineP = moduleN//':'//routineN

      INTEGER :: handle, handle2, i_global, iatom, iatom_old, icell, iiB, iquad, irow, j_global, &
         jatom, jatom_old, jcol, jjB, jkp, LLL, natom, ncol_local, nkind, nkp, nrow_local, &
         num_cells, xcell, ycell, zcell
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_from_RI_index
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_blk_end, row_blk_start, &
                                                            row_indices
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell
      LOGICAL                                            :: do_V_and_not_W
      REAL(KIND=dp) :: abs_rab_cell, arg, contribution, coskl, cutoff_exp, d_0, omega, sinkl, &
         sum_exp, sum_exp_k_im, sum_exp_k_re, tau, weight, weight_im, weight_re
      REAL(KIND=dp), DIMENSION(3)                        :: cell_vector, rab_cell_i
      REAL(KIND=dp), DIMENSION(3, 3)                     :: hmat
      REAL(KIND=dp), DIMENSION(:), POINTER               :: wkp
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: xkp
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_cfm_type), POINTER                         :: cfm_mat_L, cfm_mat_work, cfm_mat_work_2
      TYPE(cp_fm_type), POINTER                          :: fm_dummy, fm_mat_work_global, &
                                                            fm_mat_work_local
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set_RI_tmp
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      CALL timeset(routineN//"_1", handle2)

      NULLIFY (cfm_mat_work)
      CALL cp_cfm_create(cfm_mat_work, cfm_mat_Q%matrix_struct)
      CALL cp_cfm_set_all(cfm_mat_work, z_zero)

      NULLIFY (cfm_mat_work_2)
      CALL cp_cfm_create(cfm_mat_work_2, cfm_mat_Q%matrix_struct)
      CALL cp_cfm_set_all(cfm_mat_work_2, z_zero)

      NULLIFY (cfm_mat_L)
      CALL cp_cfm_create(cfm_mat_L, cfm_mat_Q%matrix_struct)
      CALL cp_cfm_set_all(cfm_mat_L, z_zero)

      ! Copy fm_mat_L_re and fm_mat_L_re to cfm_mat_L
      CALL cp_cfm_scale_and_add_fm(z_zero, cfm_mat_L, z_one, fm_mat_L_re)
      CALL cp_cfm_scale_and_add_fm(z_one, cfm_mat_L, gaussi, fm_mat_L_im)

      NULLIFY (fm_mat_work_global)
      CALL cp_fm_create(fm_mat_work_global, fm_mat_W_tau(1)%matrix%matrix_struct)
      CALL cp_fm_set_all(fm_mat_work_global, 0.0_dp)

      NULLIFY (fm_mat_work_local)
      CALL cp_fm_create(fm_mat_work_local, cfm_mat_Q%matrix_struct)
      CALL cp_fm_set_all(fm_mat_work_local, 0.0_dp)

      CALL timestop(handle2)

      IF (do_W_and_not_V) THEN

         CALL timeset(routineN//"_2", handle2)

         ! calculate [1+Q(iw')]^-1
         CALL cp_cfm_cholesky_invert(cfm_mat_Q)

         ! symmetrize the result
         CALL own_cfm_upper_to_full(cfm_mat_Q, cfm_mat_work)

         ! subtract exchange part by subtracing identity matrix from epsilon
         CALL cp_cfm_get_info(matrix=cfm_mat_Q, &
                              nrow_local=nrow_local, &
                              ncol_local=ncol_local, &
                              row_indices=row_indices, &
                              col_indices=col_indices)

         DO jjB = 1, ncol_local
            j_global = col_indices(jjB)
            DO iiB = 1, nrow_local
               i_global = row_indices(iiB)
               IF (j_global == i_global .AND. i_global <= dimen_RI) THEN
                  cfm_mat_Q%local_data(iiB, jjB) = cfm_mat_Q%local_data(iiB, jjB) - z_one
               END IF
            END DO
         END DO

         CALL timestop(handle2)

         CALL timeset(routineN//"_3.1", handle2)

         ! work = epsilon(iw,k)*L^H(k)
         CALL cp_cfm_gemm('N', 'C', dimen_RI, dimen_RI, dimen_RI, z_one, cfm_mat_Q, cfm_mat_L, &
                          z_zero, cfm_mat_work)

         ! W(iw,k) = L(k)*work
         CALL cp_cfm_gemm('N', 'N', dimen_RI, dimen_RI, dimen_RI, z_one, cfm_mat_L, cfm_mat_work, &
                          z_zero, cfm_mat_work_2)

         CALL timestop(handle2)

      ELSE

         ! S^-1(k)V(k)S^-1(k) = L(k)*L(k)^H
         CALL cp_cfm_gemm('N', 'C', dimen_RI, dimen_RI, dimen_RI, z_one, cfm_mat_L, cfm_mat_L, &
                          z_zero, cfm_mat_work_2)

      END IF

      CALL timeset(routineN//"_4", handle2)

      CALL get_kpoint_info(kpoints, xkp=xkp, wkp=wkp, nkp=nkp)
      index_to_cell => kpoints%index_to_cell
      num_cells = SIZE(index_to_cell, 2)
      d_0 = qs_env%mp2_env%ri_rpa_im_time%cutoff
      cutoff_exp = 10000.0_dp
      CALL cp_cfm_set_all(cfm_mat_work, z_zero)

      NULLIFY (qs_kind_set, cell, particle_set)
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, cell=cell, natom=natom, nkind=nkind, &
                      particle_set=particle_set)

      ALLOCATE (row_blk_start(natom))
      ALLOCATE (row_blk_end(natom))
      ALLOCATE (basis_set_RI_tmp(nkind))
      CALL basis_set_list_setup(basis_set_RI_tmp, "RI_AUX", qs_kind_set)
      CALL get_particle_set(particle_set, qs_kind_set, first_sgf=row_blk_start, last_sgf=row_blk_end, &
                            basis=basis_set_RI_tmp)
      DEALLOCATE (basis_set_RI_tmp)
      ALLOCATE (atom_from_RI_index(dimen_RI))
      DO LLL = 1, dimen_RI
         DO iatom = 1, natom
            IF (LLL >= row_blk_start(iatom) .AND. LLL <= row_blk_end(iatom)) THEN
               atom_from_RI_index(LLL) = iatom
            END IF
         END DO
      END DO
      CALL get_cell(cell=cell, h=hmat)
      iatom_old = 0
      jatom_old = 0

      CALL cp_cfm_get_info(matrix=cfm_mat_Q, &
                           nrow_local=nrow_local, &
                           ncol_local=ncol_local, &
                           row_indices=row_indices, &
                           col_indices=col_indices)

      DO irow = 1, nrow_local
         DO jcol = 1, ncol_local

            iatom = atom_from_RI_index(row_indices(irow))
            jatom = atom_from_RI_index(col_indices(jcol))

            IF (iatom .NE. iatom_old .OR. jatom .NE. jatom_old) THEN

               sum_exp = 0.0_dp
               sum_exp_k_re = 0.0_dp
               sum_exp_k_im = 0.0_dp

               DO icell = 1, num_cells

                  xcell = index_to_cell(1, icell)
                  ycell = index_to_cell(2, icell)
                  zcell = index_to_cell(3, icell)

                  arg = REAL(xcell, dp)*xkp(1, ikp) + REAL(ycell, dp)*xkp(2, ikp) + REAL(zcell, dp)*xkp(3, ikp)

                  coskl = wkp_W(ikp)*COS(twopi*arg)
                  sinkl = wkp_W(ikp)*SIN(twopi*arg)

                  cell_vector(1:3) = MATMUL(hmat, REAL(index_to_cell(1:3, icell), dp))

                  rab_cell_i(1:3) = pbc(particle_set(iatom)%r(1:3), cell) - &
                                    (pbc(particle_set(jatom)%r(1:3), cell) + cell_vector(1:3))

                  abs_rab_cell = SQRT(rab_cell_i(1)**2 + rab_cell_i(2)**2 + rab_cell_i(3)**2)

                  IF (abs_rab_cell/d_0 < cutoff_exp) THEN
                     sum_exp = sum_exp + EXP(-abs_rab_cell/d_0)
                     sum_exp_k_re = sum_exp_k_re + EXP(-abs_rab_cell/d_0)*coskl
                     sum_exp_k_im = sum_exp_k_im + EXP(-abs_rab_cell/d_0)*sinkl
                  END IF

               END DO

               weight_re = sum_exp_k_re/sum_exp
               weight_im = sum_exp_k_im/sum_exp

               iatom_old = iatom
               jatom_old = jatom

            END IF

            contribution = weight_re*REAL(cfm_mat_work_2%local_data(irow, jcol)) + &
                           weight_im*AIMAG(cfm_mat_work_2%local_data(irow, jcol))

            fm_mat_work_local%local_data(irow, jcol) = fm_mat_work_local%local_data(irow, jcol) + contribution

         END DO
      END DO

      CALL timestop(handle2)

      CALL timeset(routineN//"_5", handle2)

      IF (do_W_and_not_V) THEN

         IF (SUM(ikp_local) > nkp) THEN

            CALL cp_fm_copy_general(fm_mat_work_local, fm_mat_work_global, para_env)

            DO iquad = 1, num_integ_points

               omega = tj(jquad)
               tau = tau_tj(iquad)
               weight = weights_cos_tf_w_to_t(iquad, jquad)*COS(tau*omega)

               IF (jquad == 1 .AND. ikp == 1) THEN
                  CALL cp_fm_set_all(matrix=fm_mat_W_tau(iquad)%matrix, alpha=0.0_dp)
               END IF

               CALL cp_fm_scale_and_add(alpha=1.0_dp, matrix_a=fm_mat_W_tau(iquad)%matrix, beta=weight, matrix_b=fm_mat_work_global)

            END DO

         ELSE

            DO jkp = 1, nkp

               IF (ANY(ikp_local(:) == jkp)) THEN
                  CALL cp_fm_copy_general(fm_mat_work_local, fm_mat_work_global, para_env)
               ELSE
                  NULLIFY (fm_dummy)
                  CALL cp_fm_copy_general(fm_dummy, fm_mat_work_global, para_env)
               END IF

               DO iquad = 1, num_integ_points

                  omega = tj(jquad)
                  tau = tau_tj(iquad)
                  weight = weights_cos_tf_w_to_t(iquad, jquad)*COS(tau*omega)

                  IF (jquad == 1 .AND. jkp == 1) THEN
                     CALL cp_fm_set_all(matrix=fm_mat_W_tau(iquad)%matrix, alpha=0.0_dp)
                  END IF

                  CALL cp_fm_scale_and_add(alpha=1.0_dp, matrix_a=fm_mat_W_tau(iquad)%matrix, beta=weight, &
                                           matrix_b=fm_mat_work_global)

               END DO

            END DO

         END IF

      END IF

      do_V_and_not_W = .NOT. do_W_and_not_V
      IF (do_V_and_not_W) THEN

         IF (SUM(ikp_local) > nkp) THEN
            CALL cp_fm_copy_general(fm_mat_work_local, fm_mat_work_global, para_env)
            CALL fm_mat_work_global_to_mat_SinvVSinv(mat_SinvVSinv, fm_mat_work_global)
         ELSE
            DO jkp = 1, nkp
               IF (ANY(ikp_local(:) == jkp)) THEN
                  CALL cp_fm_copy_general(fm_mat_work_local, fm_mat_work_global, para_env)
               ELSE
                  NULLIFY (fm_dummy)
                  CALL cp_fm_copy_general(fm_dummy, fm_mat_work_global, para_env)
               END IF
               CALL fm_mat_work_global_to_mat_SinvVSinv(mat_SinvVSinv, fm_mat_work_global)
            END DO
         END IF
      END IF

      CALL cp_cfm_release(cfm_mat_work)
      CALL cp_cfm_release(cfm_mat_work_2)
      CALL cp_cfm_release(cfm_mat_L)
      CALL cp_fm_release(fm_mat_work_global)
      CALL cp_fm_release(fm_mat_work_local)
      DEALLOCATE (atom_from_RI_index)
      DEALLOCATE (row_blk_start)
      DEALLOCATE (row_blk_end)

      CALL timestop(handle2)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param mat_SinvVSinv ...
!> \param fm_mat_work_global ...
! **************************************************************************************************
   SUBROUTINE fm_mat_work_global_to_mat_SinvVSinv(mat_SinvVSinv, fm_mat_work_global)

      TYPE(dbcsr_p_type)                                 :: mat_SinvVSinv
      TYPE(cp_fm_type), POINTER                          :: fm_mat_work_global

      CHARACTER(LEN=*), PARAMETER :: routineN = 'fm_mat_work_global_to_mat_SinvVSinv', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      TYPE(dbcsr_p_type)                                 :: mat_work

      CALL timeset(routineN, handle)

      NULLIFY (mat_work%matrix)
      ALLOCATE (mat_work%matrix)
      CALL dbcsr_create(mat_work%matrix, template=mat_SinvVSinv%matrix)

      CALL copy_fm_to_dbcsr(fm_mat_work_global, mat_work%matrix, keep_sparsity=.FALSE.)

      CALL dbcsr_add(mat_SinvVSinv%matrix, mat_work%matrix, 1.0_dp, 1.0_dp)

      CALL dbcsr_release(mat_work%matrix)
      DEALLOCATE (mat_work%matrix)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param cfm_mat_W_kp_tau ...
!> \param cfm_mat_Q ...
!> \param fm_mat_L_re ...
!> \param fm_mat_L_im ...
!> \param dimen_RI ...
!> \param num_integ_points ...
!> \param jquad ...
!> \param ikp ...
!> \param tj ...
!> \param tau_tj ...
!> \param weights_cos_tf_w_to_t ...
! **************************************************************************************************
   SUBROUTINE compute_Wc_kp_tau_GW(cfm_mat_W_kp_tau, cfm_mat_Q, fm_mat_L_re, fm_mat_L_im, &
                                   dimen_RI, num_integ_points, jquad, &
                                   ikp, tj, tau_tj, weights_cos_tf_w_to_t)

      TYPE(cp_cfm_p_type), DIMENSION(:, :), POINTER      :: cfm_mat_W_kp_tau
      TYPE(cp_cfm_type), POINTER                         :: cfm_mat_Q
      TYPE(cp_fm_type), POINTER                          :: fm_mat_L_re, fm_mat_L_im
      INTEGER                                            :: dimen_RI, num_integ_points, jquad, ikp
      REAL(KIND=dp), DIMENSION(:)                        :: tj
      REAL(KIND=dp), DIMENSION(0:num_integ_points)       :: tau_tj
      REAL(KIND=dp), DIMENSION(:, :)                     :: weights_cos_tf_w_to_t

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_Wc_kp_tau_GW', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, handle2, i_global, iiB, iquad, &
                                                            j_global, jjB, ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: omega, tau, weight
      TYPE(cp_cfm_type), POINTER                         :: cfm_mat_L, cfm_mat_work

      CALL timeset(routineN, handle)

      NULLIFY (cfm_mat_work)
      CALL cp_cfm_create(cfm_mat_work, fm_mat_L_re%matrix_struct)
      CALL cp_cfm_set_all(cfm_mat_work, z_zero)

      NULLIFY (cfm_mat_L)
      CALL cp_cfm_create(cfm_mat_L, fm_mat_L_re%matrix_struct)
      CALL cp_cfm_set_all(cfm_mat_L, z_zero)

      CALL timeset(routineN//"_cholesky_inv", handle2)

      ! calculate [1+Q(iw')]^-1
      CALL cp_cfm_cholesky_invert(cfm_mat_Q)

      ! symmetrize the result
      CALL own_cfm_upper_to_full(cfm_mat_Q, cfm_mat_work)

      ! subtract exchange part by subtracing identity matrix from epsilon
      CALL cp_cfm_get_info(matrix=cfm_mat_Q, &
                           nrow_local=nrow_local, &
                           ncol_local=ncol_local, &
                           row_indices=row_indices, &
                           col_indices=col_indices)

      DO jjB = 1, ncol_local
         j_global = col_indices(jjB)
         DO iiB = 1, nrow_local
            i_global = row_indices(iiB)
            IF (j_global == i_global .AND. i_global <= dimen_RI) THEN
               cfm_mat_Q%local_data(iiB, jjB) = cfm_mat_Q%local_data(iiB, jjB) - z_one
            END IF
         END DO
      END DO

      CALL timestop(handle2)

      ! Copy fm_mat_L_re and fm_mat_L_re to cfm_mat_L
      CALL cp_cfm_scale_and_add_fm(z_zero, cfm_mat_L, z_one, fm_mat_L_re)
      CALL cp_cfm_scale_and_add_fm(z_one, cfm_mat_L, gaussi, fm_mat_L_im)

      ! work = epsilon(iw,k)*L^H(k)
      CALL cp_cfm_gemm('N', 'C', dimen_RI, dimen_RI, dimen_RI, z_one, cfm_mat_Q, cfm_mat_L, &
                       z_zero, cfm_mat_work)

      ! W(iw,k) = L(k)*work
      CALL cp_cfm_gemm('N', 'N', dimen_RI, dimen_RI, dimen_RI, z_one, cfm_mat_L, cfm_mat_work, &
                       z_zero, cfm_mat_Q)

      DO iquad = 1, num_integ_points
         omega = tj(jquad)
         tau = tau_tj(iquad)
         weight = weights_cos_tf_w_to_t(iquad, jquad)*COS(tau*omega)
         CALL cp_cfm_scale_and_add(alpha=z_one, matrix_a=cfm_mat_W_kp_tau(ikp, iquad)%matrix, &
                                   beta=CMPLX(weight, KIND=dp), matrix_b=cfm_mat_Q)
      END DO

      CALL cp_cfm_release(cfm_mat_work)
      CALL cp_cfm_release(cfm_mat_L)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param wkp_W ...
!> \param kpoints ...
!> \param h_mat ...
!> \param h_inv ...
!> \param exp_kpoints ...
!> \param periodic ...
! **************************************************************************************************
   SUBROUTINE compute_wkp_W(wkp_W, kpoints, h_mat, h_inv, exp_kpoints, periodic)
      REAL(KIND=dp), DIMENSION(:), POINTER               :: wkp_W
      TYPE(kpoint_type), POINTER                         :: kpoints
      REAL(KIND=dp), DIMENSION(3, 3)                     :: h_mat, h_inv
      REAL(KIND=dp)                                      :: exp_kpoints
      INTEGER, DIMENSION(3)                              :: periodic

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_wkp_W'

      INTEGER                                            :: handle, i_dim, i_x, ikp, info, j_y, k_z, &
                                                            n_x, n_y, n_z, nkp, nsuperfine, &
                                                            num_lin_eqs
      REAL(KIND=dp)                                      :: a_vec_dot_k_vec, integral, k_sq, weight
      REAL(KIND=dp), DIMENSION(3)                        :: a_vec, k_vec, x_vec
      REAL(KIND=dp), DIMENSION(:), POINTER               :: right_side, wkp
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: matrix_lin_eqs, xkp

      CALL timeset(routineN, handle)

      CALL get_kpoint_info(kpoints, xkp=xkp, wkp=wkp, nkp=nkp)

      ! we determine the kpoint weights of the Monkhors Pack mesh new
      ! such that the functions 1/k^2, 1/k and const are integrated exactly
      ! in the Brillouin zone
      ! this is done by minimizing sum_i |w_i|^2 where w_i are the weights of
      ! the i-th kpoint under the following constraints:
      ! 1) 1/k^2, 1/k and const are integrated exactly
      ! 2) the kpoint weights of kpoints with identical absolute value are
      !    the same, of e.g. (1/8,3/8,3/8) same weight as for (3/8,1/8,3/8)
      ! for 1d and 2d materials: we use normal Monkhorst-Pack mesh, checked
      ! by SUM(periodic) == 3

      IF (exp_kpoints < 2.0_dp .AND. SUM(periodic) == 3) THEN

         ! first, compute the integral of f(k)=1/k^2 and 1/k on super fine grid
         nsuperfine = 500
         integral = 0.0_dp
         IF (exp_kpoints > 0.0_dp) exp_kpoints = -2.0_dp

         ! actually, there is the factor *det_3x3(h_inv) missing to account for the
         ! integration volume but for wkp det_3x3(h_inv) is needed
         weight = 2.0_dp/(REAL(nsuperfine, dp))**3
         DO i_x = 1, nsuperfine
            DO j_y = 1, nsuperfine
               DO k_z = 1, nsuperfine/2

                  x_vec = (/REAL(i_x - nsuperfine/2, dp) - 0.5_dp, &
                            REAL(j_y - nsuperfine/2, dp) - 0.5_dp, &
                            REAL(k_z - nsuperfine/2, dp) - 0.5_dp/)/ &
                          REAL(nsuperfine, dp)
                  k_vec = MATMUL(h_inv(1:3, 1:3), x_vec)
                  k_sq = k_vec(1)**2 + k_vec(2)**2 + k_vec(3)**2
                  integral = integral + weight*k_sq**(exp_kpoints*0.5_dp)
               END DO
            END DO
         END DO

         num_lin_eqs = nkp + 2

         ALLOCATE (matrix_lin_eqs(num_lin_eqs, num_lin_eqs))
         matrix_lin_eqs(:, :) = 0.0_dp

         DO ikp = 1, nkp

            k_vec = MATMUL(h_inv(1:3, 1:3), xkp(1:3, ikp))
            k_sq = k_vec(1)**2 + k_vec(2)**2 + k_vec(3)**2

            matrix_lin_eqs(ikp, ikp) = 2.0_dp
            matrix_lin_eqs(ikp, nkp + 1) = 1.0_dp
            matrix_lin_eqs(ikp, nkp + 2) = 1.0_dp*k_sq**(exp_kpoints*0.5_dp)

            matrix_lin_eqs(nkp + 1, ikp) = 1.0_dp
            matrix_lin_eqs(nkp + 2, ikp) = 1.0_dp*k_sq**(exp_kpoints*0.5_dp)

         END DO

         CALL invmat(matrix_lin_eqs, info)
         ! check whether inversion was successful
         CPASSERT(info == 0)

         ALLOCATE (wkp_W(num_lin_eqs))

         ALLOCATE (right_side(num_lin_eqs))
         right_side = 0.0_dp
         right_side(nkp + 1) = 1.0_dp
         right_side(nkp + 2) = integral

         wkp_W(1:num_lin_eqs) = MATMUL(matrix_lin_eqs, right_side)

         DEALLOCATE (matrix_lin_eqs, right_side)

      ELSE IF (exp_kpoints < 2.0_dp .AND. SUM(periodic) == 1) THEN

         ! first, compute the integral of f(k)=1/k^2 and 1/k on super fine grid
         nsuperfine = 5000
         integral = 0.0_dp

         ! actually, there is the factor *det_3x3(h_inv) missing to account for the
         ! integration volume but for wkp det_3x3(h_inv) is needed
         weight = 1.0_dp/REAL(nsuperfine, dp)
         IF (periodic(1) == 1) THEN
            n_x = nsuperfine
         ELSE
            n_x = 1
         END IF
         IF (periodic(2) == 1) THEN
            n_y = nsuperfine
         ELSE
            n_y = 1
         END IF
         IF (periodic(3) == 1) THEN
            n_z = nsuperfine
         ELSE
            n_z = 1
         END IF

         a_vec = MATMUL(h_mat(1:3, 1:3), &
                        (/REAL(periodic(1), dp), REAL(periodic(2), dp), REAL(periodic(3), dp)/))

         DO i_x = 1, n_x
            DO j_y = 1, n_y
               DO k_z = 1, n_z

                  x_vec = (/REAL(i_x - nsuperfine/2, dp) - 0.5_dp, &
                            REAL(j_y - nsuperfine/2, dp) - 0.5_dp, &
                            REAL(k_z - nsuperfine/2, dp) - 0.5_dp/)/ &
                          REAL(nsuperfine, dp)

                  DO i_dim = 1, 3
                     IF (periodic(i_dim) == 0) THEN
                        x_vec(i_dim) = 0.0_dp
                     END IF
                  END DO

                  k_vec = MATMUL(h_inv(1:3, 1:3), x_vec)
                  a_vec_dot_k_vec = a_vec(1)*k_vec(1) + a_vec(2)*k_vec(2) + a_vec(3)*k_vec(3)
                  integral = integral + weight*LOG(2.0_dp - 2.0_dp*COS(a_vec_dot_k_vec))
               END DO
            END DO
         END DO

         num_lin_eqs = nkp + 2

         ALLOCATE (matrix_lin_eqs(num_lin_eqs, num_lin_eqs))
         matrix_lin_eqs(:, :) = 0.0_dp

         DO ikp = 1, nkp

            k_vec = MATMUL(h_inv(1:3, 1:3), xkp(1:3, ikp))
            k_sq = k_vec(1)**2 + k_vec(2)**2 + k_vec(3)**2

            matrix_lin_eqs(ikp, ikp) = 2.0_dp
            matrix_lin_eqs(ikp, nkp + 1) = 1.0_dp

            a_vec_dot_k_vec = a_vec(1)*k_vec(1) + a_vec(2)*k_vec(2) + a_vec(3)*k_vec(3)
            matrix_lin_eqs(ikp, nkp + 2) = LOG(2.0_dp - 2.0_dp*COS(a_vec_dot_k_vec))

            matrix_lin_eqs(nkp + 1, ikp) = 1.0_dp
            matrix_lin_eqs(nkp + 2, ikp) = LOG(2.0_dp - 2.0_dp*COS(a_vec_dot_k_vec))

         END DO

         CALL invmat(matrix_lin_eqs, info)
         ! check whether inversion was successful
         CPASSERT(info == 0)

         ALLOCATE (wkp_W(num_lin_eqs))

         ALLOCATE (right_side(num_lin_eqs))
         right_side = 0.0_dp
         right_side(nkp + 1) = 1.0_dp
         right_side(nkp + 2) = integral

         wkp_W(1:num_lin_eqs) = MATMUL(matrix_lin_eqs, right_side)

         DEALLOCATE (matrix_lin_eqs, right_side)

      ELSE

         ALLOCATE (wkp_W(nkp))
         wkp_W(:) = wkp(:)

      END IF

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param cfm_mat_Q ...
!> \param cfm_mat_work ...
! **************************************************************************************************
   SUBROUTINE own_cfm_upper_to_full(cfm_mat_Q, cfm_mat_work)

      TYPE(cp_cfm_type), POINTER                         :: cfm_mat_Q, cfm_mat_work

      CHARACTER(LEN=*), PARAMETER :: routineN = 'own_cfm_upper_to_full', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i_global, iiB, j_global, jjB, &
                                                            ncol_local, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices

      CALL timeset(routineN, handle)

      ! get info of fm_mat_Q
      CALL cp_cfm_get_info(matrix=cfm_mat_Q, &
                           nrow_local=nrow_local, &
                           ncol_local=ncol_local, &
                           row_indices=row_indices, &
                           col_indices=col_indices)

      DO jjB = 1, ncol_local
         j_global = col_indices(jjB)
         DO iiB = 1, nrow_local
            i_global = row_indices(iiB)
            IF (j_global < i_global) THEN
               cfm_mat_Q%local_data(iiB, jjB) = z_zero
            END IF
            IF (j_global == i_global) THEN
               cfm_mat_Q%local_data(iiB, jjB) = cfm_mat_Q%local_data(iiB, jjB)/(2.0_dp, 0.0_dp)
            END IF
         END DO
      END DO

      CALL cp_cfm_transpose(cfm_mat_Q, 'C', cfm_mat_work)

      CALL cp_cfm_scale_and_add(z_one, cfm_mat_Q, z_one, cfm_mat_work)

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param mat_contr_W_re ...
!> \param mat_contr_W_im ...
!> \param n_level_gw ...
!> \param jquad ...
!> \param gw_corr_lev_occ ...
!> \param homo ...
!> \param mat_contr_gf_occ_re ...
!> \param mat_contr_gf_occ_im ...
!> \param mat_contr_gf_virt_re ...
!> \param mat_contr_gf_virt_im ...
!> \param vec_Sigma_c_gw_pos_tau ...
!> \param vec_Sigma_c_gw_neg_tau ...
!> \param vec_Sigma_x_gw ...
!> \param ikp ...
!> \param cycle_R1_S2_n_level ...
!> \param cycle_R1_plus_S2_n_level ...
!> \param eps_filter ...
!> \param i_cell_R1_plus_S2 ...
!> \param i_cell_S2 ...
!> \param start_jquad ...
! **************************************************************************************************
   SUBROUTINE trace_for_self_ener(mat_contr_W_re, mat_contr_W_im, n_level_gw, jquad, &
                                  gw_corr_lev_occ, homo, &
                                  mat_contr_gf_occ_re, mat_contr_gf_occ_im, &
                                  mat_contr_gf_virt_re, mat_contr_gf_virt_im, &
                                  vec_Sigma_c_gw_pos_tau, vec_Sigma_c_gw_neg_tau, &
                                  vec_Sigma_x_gw, ikp, &
                                  cycle_R1_S2_n_level, cycle_R1_plus_S2_n_level, &
                                  eps_filter, i_cell_R1_plus_S2, i_cell_S2, &
                                  start_jquad)

      TYPE(dbcsr_type), POINTER                          :: mat_contr_W_re, mat_contr_W_im
      INTEGER                                            :: n_level_gw, jquad, gw_corr_lev_occ, homo
      TYPE(dbcsr_type), POINTER                          :: mat_contr_gf_occ_re, &
                                                            mat_contr_gf_occ_im, &
                                                            mat_contr_gf_virt_re, &
                                                            mat_contr_gf_virt_im
      COMPLEX(KIND=dp), DIMENSION(:, :)                  :: vec_Sigma_c_gw_pos_tau, &
                                                            vec_Sigma_c_gw_neg_tau
      REAL(KIND=dp), DIMENSION(:)                        :: vec_Sigma_x_gw
      INTEGER                                            :: ikp
      LOGICAL, ALLOCATABLE, DIMENSION(:, :, :, :)        :: cycle_R1_S2_n_level
      LOGICAL, ALLOCATABLE, DIMENSION(:, :, :)           :: cycle_R1_plus_S2_n_level
      REAL(KIND=dp)                                      :: eps_filter
      INTEGER                                            :: i_cell_R1_plus_S2, i_cell_S2, start_jquad

      CHARACTER(LEN=*), PARAMETER :: routineN = 'trace_for_self_ener', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, n_level_gw_ref
      REAL(KIND=dp) :: trace_neg_tau_im_1, trace_neg_tau_im_2, trace_neg_tau_re_1, &
         trace_neg_tau_re_2, trace_pos_tau_im_1, trace_pos_tau_im_2, trace_pos_tau_re_1, &
         trace_pos_tau_re_2

      CALL timeset(routineN, handle)

      CALL dbcsr_dot(mat_contr_gf_occ_re, &
                     mat_contr_W_re, &
                     trace_neg_tau_re_1)

      CALL dbcsr_dot(mat_contr_gf_occ_im, &
                     mat_contr_W_im, &
                     trace_neg_tau_re_2)

      CALL dbcsr_dot(mat_contr_gf_occ_re, &
                     mat_contr_W_im, &
                     trace_neg_tau_im_1)

      CALL dbcsr_dot(mat_contr_gf_occ_im, &
                     mat_contr_W_re, &
                     trace_neg_tau_im_2)

      IF (ABS(trace_neg_tau_re_1) + ABS(trace_neg_tau_re_2) + ABS(trace_neg_tau_im_1) + &
          ABS(trace_neg_tau_im_2) < eps_filter) THEN
         cycle_R1_S2_n_level(i_cell_R1_plus_S2, i_cell_S2, n_level_gw, jquad) = .TRUE.
      END IF

      IF (ikp == 1 .AND. jquad == start_jquad .AND. i_cell_S2 == 1) THEN
         cycle_R1_plus_S2_n_level(i_cell_R1_plus_S2, n_level_gw, jquad) = .TRUE.
      END IF

      IF (ABS(trace_neg_tau_re_1) + ABS(trace_neg_tau_re_2) + ABS(trace_neg_tau_im_1) + &
          ABS(trace_neg_tau_im_2) > eps_filter) THEN
         cycle_R1_plus_S2_n_level(i_cell_R1_plus_S2, n_level_gw, jquad) = .FALSE.
      END IF

      IF (jquad == 0) THEN

         n_level_gw_ref = n_level_gw + homo - gw_corr_lev_occ

         vec_Sigma_x_gw(n_level_gw_ref) = vec_Sigma_x_gw(n_level_gw_ref) + trace_neg_tau_re_1 - trace_neg_tau_re_2

      ELSE

         vec_Sigma_c_gw_neg_tau(n_level_gw, jquad) = vec_Sigma_c_gw_neg_tau(n_level_gw, jquad) + &
                                                     CMPLX(trace_neg_tau_re_1 - trace_neg_tau_re_2, &
                                                           trace_neg_tau_im_1 + trace_neg_tau_im_2, dp)

      END IF

      CALL dbcsr_dot(mat_contr_gf_virt_re, &
                     mat_contr_W_re, &
                     trace_pos_tau_re_1)

      CALL dbcsr_dot(mat_contr_gf_virt_im, &
                     mat_contr_W_im, &
                     trace_pos_tau_re_2)

      CALL dbcsr_dot(mat_contr_gf_virt_re, &
                     mat_contr_W_im, &
                     trace_pos_tau_im_1)

      CALL dbcsr_dot(mat_contr_gf_virt_im, &
                     mat_contr_W_re, &
                     trace_pos_tau_im_2)

      IF (jquad > 0) THEN

         vec_Sigma_c_gw_pos_tau(n_level_gw, jquad) = vec_Sigma_c_gw_pos_tau(n_level_gw, jquad) + &
                                                     CMPLX(trace_pos_tau_re_1 - trace_pos_tau_re_2, &
                                                           trace_pos_tau_im_1 + trace_pos_tau_im_2, dp)

      END IF

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param mat_I_muP_occ_re ...
!> \param mat_I_muP_virt_re ...
!> \param mat_I_muP_occ_im ...
!> \param mat_I_muP_virt_im ...
!> \param mat_contr_gf_occ_re ...
!> \param mat_contr_gf_occ_im ...
!> \param mat_contr_gf_virt_re ...
!> \param mat_contr_gf_virt_im ...
!> \param index_to_cell_3c ...
!> \param kpoints ...
!> \param ikp ...
!> \param x_cell_R1 ...
!> \param y_cell_R1 ...
!> \param z_cell_R1 ...
! **************************************************************************************************
   SUBROUTINE trafo_I_T_R1_plus_S2_to_M_R1_S2(mat_I_muP_occ_re, mat_I_muP_virt_re, &
                                              mat_I_muP_occ_im, mat_I_muP_virt_im, &
                                              mat_contr_gf_occ_re, mat_contr_gf_occ_im, &
                                              mat_contr_gf_virt_re, mat_contr_gf_virt_im, &
                                              index_to_cell_3c, kpoints, ikp, &
                                              x_cell_R1, y_cell_R1, z_cell_R1)

      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mat_I_muP_occ_re, mat_I_muP_virt_re, &
                                                            mat_I_muP_occ_im, mat_I_muP_virt_im
      TYPE(dbcsr_type), POINTER                          :: mat_contr_gf_occ_re, &
                                                            mat_contr_gf_occ_im, &
                                                            mat_contr_gf_virt_re, &
                                                            mat_contr_gf_virt_im
      INTEGER, DIMENSION(:, :)                           :: index_to_cell_3c
      TYPE(kpoint_type), POINTER                         :: kpoints
      INTEGER                                            :: ikp, x_cell_R1, y_cell_R1, z_cell_R1

      CHARACTER(LEN=*), PARAMETER :: routineN = 'trafo_I_T_R1_plus_S2_to_M_R1_S2', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i_cell_T, num_cells_3c, xcell, &
                                                            ycell, zcell
      REAL(KIND=dp)                                      :: arg, coskl, sinkl
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: xkp

      CALL timeset(routineN, handle)

      num_cells_3c = SIZE(mat_I_muP_occ_re)
      CALL get_kpoint_info(kpoints, xkp=xkp)

      CALL dbcsr_set(mat_contr_gf_occ_re, 0.0_dp)
      CALL dbcsr_set(mat_contr_gf_occ_im, 0.0_dp)
      CALL dbcsr_set(mat_contr_gf_virt_re, 0.0_dp)
      CALL dbcsr_set(mat_contr_gf_virt_im, 0.0_dp)

      DO i_cell_T = 1, num_cells_3c

         xcell = index_to_cell_3c(1, i_cell_T) - x_cell_R1
         ycell = index_to_cell_3c(2, i_cell_T) - y_cell_R1
         zcell = index_to_cell_3c(3, i_cell_T) - z_cell_R1
         arg = REAL(xcell, dp)*xkp(1, ikp) + REAL(ycell, dp)*xkp(2, ikp) + REAL(zcell, dp)*xkp(3, ikp)
         coskl = COS(twopi*arg)
         sinkl = SIN(twopi*arg)

         CALL dbcsr_scale_and_add_local(mat_contr_gf_occ_re, mat_I_muP_occ_re(i_cell_T)%matrix, coskl)
         CALL dbcsr_scale_and_add_local(mat_contr_gf_occ_re, mat_I_muP_occ_im(i_cell_T)%matrix, -sinkl)
         CALL dbcsr_scale_and_add_local(mat_contr_gf_occ_im, mat_I_muP_occ_re(i_cell_T)%matrix, sinkl)
         CALL dbcsr_scale_and_add_local(mat_contr_gf_occ_im, mat_I_muP_occ_im(i_cell_T)%matrix, coskl)

         CALL dbcsr_scale_and_add_local(mat_contr_gf_virt_re, mat_I_muP_virt_re(i_cell_T)%matrix, coskl)
         CALL dbcsr_scale_and_add_local(mat_contr_gf_virt_re, mat_I_muP_virt_im(i_cell_T)%matrix, -sinkl)
         CALL dbcsr_scale_and_add_local(mat_contr_gf_virt_im, mat_I_muP_virt_re(i_cell_T)%matrix, sinkl)
         CALL dbcsr_scale_and_add_local(mat_contr_gf_virt_im, mat_I_muP_virt_im(i_cell_T)%matrix, coskl)

      END DO

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param mat_A ...
!> \param mat_B ...
!> \param beta ...
! **************************************************************************************************
   SUBROUTINE dbcsr_scale_and_add_local(mat_A, mat_B, beta)
      TYPE(dbcsr_type), POINTER                          :: mat_A, mat_B
      REAL(KIND=dp)                                      :: beta

      INTEGER                                            :: col, row
      LOGICAL                                            :: found
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: block_to_compute, data_block
      TYPE(dbcsr_iterator_type)                          :: iter

      CALL dbcsr_iterator_start(iter, mat_B)
      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, row, col, data_block)

         NULLIFY (block_to_compute)
         CALL dbcsr_get_block_p(matrix=mat_A, &
                                row=row, col=col, block=block_to_compute, found=found)

         CPASSERT(found)

         block_to_compute(:, :) = block_to_compute(:, :) + beta*data_block(:, :)

      END DO

      CALL dbcsr_iterator_stop(iter)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param mat_dm_occ_S ...
!> \param mat_dm_virt_S ...
!> \param qs_env ...
!> \param num_integ_points ...
!> \param e_fermi ...
!> \param eps_filter ...
!> \param tau_tj ...
!> \param num_cells_dm ...
!> \param index_to_cell_dm ...
!> \param has_blocks_mat_dm_occ_S ...
!> \param has_blocks_mat_dm_virt_S ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE compute_G_real_space(mat_dm_occ_S, mat_dm_virt_S, qs_env, num_integ_points, &
                                   e_fermi, eps_filter, tau_tj, num_cells_dm, index_to_cell_dm, &
                                   has_blocks_mat_dm_occ_S, has_blocks_mat_dm_virt_S, para_env)

      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: mat_dm_occ_S, mat_dm_virt_S
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER                                            :: num_integ_points
      REAL(KIND=dp)                                      :: e_fermi, eps_filter
      REAL(KIND=dp), DIMENSION(0:num_integ_points)       :: tau_tj
      INTEGER                                            :: num_cells_dm
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell_dm
      LOGICAL, ALLOCATABLE, DIMENSION(:, :)              :: has_blocks_mat_dm_occ_S, &
                                                            has_blocks_mat_dm_virt_S
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_G_real_space', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i_cell, ispin, jquad, nblks
      REAL(KIND=dp)                                      :: tau

      CALL timeset(routineN, handle)

      ispin = 1

      ! get denity matrix for exchange self-energy
      tau = 0.0_dp
      CALL compute_transl_dm(mat_dm_occ_S, qs_env, ispin, num_integ_points, 0, e_fermi, tau, &
                             eps_filter, num_cells_dm, index_to_cell_dm, &
                             remove_occ=.FALSE., remove_virt=.TRUE., first_jquad=0)

      DO i_cell = 1, num_cells_dm

         nblks = dbcsr_get_num_blocks(mat_dm_occ_S(0, i_cell)%matrix)
         CALL mp_sum(nblks, para_env%group)
         IF (nblks == 0) has_blocks_mat_dm_occ_S(0, i_cell) = .FALSE.
         IF (nblks > 0) has_blocks_mat_dm_occ_S(0, i_cell) = .TRUE.

      END DO

      DO jquad = 1, num_integ_points

         tau = tau_tj(jquad)

         CALL compute_transl_dm(mat_dm_occ_S, qs_env, ispin, num_integ_points, jquad, e_fermi, tau, &
                                eps_filter, num_cells_dm, index_to_cell_dm, &
                                remove_occ=.FALSE., remove_virt=.TRUE., first_jquad=0)

         CALL compute_transl_dm(mat_dm_virt_S, qs_env, ispin, num_integ_points, jquad, e_fermi, tau, &
                                eps_filter, num_cells_dm, index_to_cell_dm, &
                                remove_occ=.TRUE., remove_virt=.FALSE., first_jquad=1)

         DO i_cell = 1, num_cells_dm

            nblks = dbcsr_get_num_blocks(mat_dm_occ_S(jquad, i_cell)%matrix)
            CALL mp_sum(nblks, para_env%group)
            IF (nblks == 0) has_blocks_mat_dm_occ_S(jquad, i_cell) = .FALSE.
            IF (nblks > 0) has_blocks_mat_dm_occ_S(jquad, i_cell) = .TRUE.

            nblks = dbcsr_get_num_blocks(mat_dm_virt_S(jquad, i_cell)%matrix)
            CALL mp_sum(nblks, para_env%group)
            IF (nblks == 0) has_blocks_mat_dm_virt_S(jquad, i_cell) = .FALSE.
            IF (nblks > 0) has_blocks_mat_dm_virt_S(jquad, i_cell) = .TRUE.

         END DO

      END DO

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param index_to_cell_R1_plus_S2 ...
!> \param cell_to_index_3c ...
!> \param num_cells_R1_plus_S2 ...
!> \param kpoints ...
!> \param unit_nr ...
!> \param maxcell ...
!> \param num_cells_dm ...
! **************************************************************************************************
   SUBROUTINE compute_cell_vec_for_R1_plus_S2_or_R1(index_to_cell_R1_plus_S2, cell_to_index_3c, &
                                                    num_cells_R1_plus_S2, kpoints, &
                                                    unit_nr, maxcell, num_cells_dm)

      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: index_to_cell_R1_plus_S2
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: cell_to_index_3c
      INTEGER                                            :: num_cells_R1_plus_S2
      TYPE(kpoint_type), POINTER                         :: kpoints
      INTEGER                                            :: unit_nr, maxcell, num_cells_dm

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_cell_vec_for_R1_plus_S2_or_R1', &
         routineP = moduleN//':'//routineN

      INTEGER :: bound_1, bound_2, bound_3, bound_4, bound_5, bound_6, handle, i_cell_S_1, &
         size_set, x_cell_R1_plus_S2, x_cell_R1_plus_S2_minus_S1, y_cell_R1_plus_S2, &
         y_cell_R1_plus_S2_minus_S1, z_cell_R1_plus_S2, z_cell_R1_plus_S2_minus_S1
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: index_to_cell_tmp
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell_dm
      LOGICAL                                            :: already_there

      CALL timeset(routineN, handle)

      index_to_cell_dm => kpoints%index_to_cell

      ALLOCATE (index_to_cell_R1_plus_S2(3, 1))
      index_to_cell_R1_plus_S2 = 0

      bound_1 = LBOUND(cell_to_index_3c, 1)
      bound_2 = UBOUND(cell_to_index_3c, 1)
      bound_3 = LBOUND(cell_to_index_3c, 2)
      bound_4 = UBOUND(cell_to_index_3c, 2)
      bound_5 = LBOUND(cell_to_index_3c, 3)
      bound_6 = UBOUND(cell_to_index_3c, 3)

      maxcell = 8

      DO x_cell_R1_plus_S2 = -maxcell, maxcell
         DO y_cell_R1_plus_S2 = -maxcell, maxcell
            DO z_cell_R1_plus_S2 = -maxcell, maxcell

               already_there = .FALSE.

               DO i_cell_S_1 = 1, num_cells_dm

                  IF (already_there) CYCLE

                  x_cell_R1_plus_S2_minus_S1 = x_cell_R1_plus_S2 - index_to_cell_dm(1, i_cell_S_1)
                  y_cell_R1_plus_S2_minus_S1 = y_cell_R1_plus_S2 - index_to_cell_dm(2, i_cell_S_1)
                  z_cell_R1_plus_S2_minus_S1 = z_cell_R1_plus_S2 - index_to_cell_dm(3, i_cell_S_1)

                  IF (x_cell_R1_plus_S2_minus_S1 .GE. bound_1 .AND. &
                      x_cell_R1_plus_S2_minus_S1 .LE. bound_2 .AND. &
                      y_cell_R1_plus_S2_minus_S1 .GE. bound_3 .AND. &
                      y_cell_R1_plus_S2_minus_S1 .LE. bound_4 .AND. &
                      z_cell_R1_plus_S2_minus_S1 .GE. bound_5 .AND. &
                      z_cell_R1_plus_S2_minus_S1 .LE. bound_6) THEN

                     size_set = SIZE(index_to_cell_R1_plus_S2, 2)

                     IF (size_set == 1 .AND. &
                         index_to_cell_R1_plus_S2(1, size_set) == 0 .AND. &
                         index_to_cell_R1_plus_S2(2, size_set) == 0 .AND. &
                         index_to_cell_R1_plus_S2(3, size_set) == 0) THEN

                        index_to_cell_R1_plus_S2(1, size_set) = x_cell_R1_plus_S2
                        index_to_cell_R1_plus_S2(2, size_set) = y_cell_R1_plus_S2
                        index_to_cell_R1_plus_S2(3, size_set) = z_cell_R1_plus_S2

                     ELSE

                        ALLOCATE (index_to_cell_tmp(3, size_set))
                        index_to_cell_tmp(1:3, 1:size_set) = index_to_cell_R1_plus_S2(1:3, 1:size_set)
                        DEALLOCATE (index_to_cell_R1_plus_S2)
                        ALLOCATE (index_to_cell_R1_plus_S2(3, size_set + 1))
                        index_to_cell_R1_plus_S2(1:3, 1:size_set) = index_to_cell_tmp(1:3, 1:size_set)
                        index_to_cell_R1_plus_S2(1, size_set + 1) = x_cell_R1_plus_S2
                        index_to_cell_R1_plus_S2(2, size_set + 1) = y_cell_R1_plus_S2
                        index_to_cell_R1_plus_S2(3, size_set + 1) = z_cell_R1_plus_S2
                        DEALLOCATE (index_to_cell_tmp)
                        already_there = .TRUE.

                     END IF

                  END IF

               END DO

            END DO
         END DO
      END DO

      num_cells_R1_plus_S2 = SIZE(index_to_cell_R1_plus_S2, 2)

      IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
         "GW_INFO| Number of periodic images for R_1+S_2 and R_2 sum in self-energy:", num_cells_R1_plus_S2

      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param index_to_cell_R2 ...
!> \param cell_to_index_R2 ...
!> \param num_cells_R2 ...
!> \param unit_nr ...
!> \param index_to_cell_dm ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE compute_cell_vec_for_R2(index_to_cell_R2, cell_to_index_R2, num_cells_R2, unit_nr, &
                                      index_to_cell_dm, qs_env)

      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: index_to_cell_R2
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: cell_to_index_R2
      INTEGER                                            :: num_cells_R2, unit_nr
      INTEGER, DIMENSION(:, :), POINTER                  :: index_to_cell_dm
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compute_cell_vec_for_R2', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: bound_1, bound_2, bound_3, bound_4, &
                                                            bound_5, bound_6, handle, icell, &
                                                            num_cells_dm
      TYPE(kpoint_type), POINTER                         :: kpoints

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, kpoints=kpoints)

      index_to_cell_dm => kpoints%index_to_cell

      num_cells_dm = SIZE(index_to_cell_dm, 2)

      ALLOCATE (index_to_cell_R2(3, num_cells_dm))
      index_to_cell_R2(1:3, 1:num_cells_dm) = index_to_cell_dm(1:3, 1:num_cells_dm)

      num_cells_R2 = SIZE(index_to_cell_R2, 2)

      IF (unit_nr > 0) WRITE (UNIT=unit_nr, FMT="(T3,A,T75,i6)") &
         "GW_INFO| Number of periodic images for R_2 W-sum in self-energy:", num_cells_R2

      bound_1 = MINVAL(index_to_cell_R2(1, :))
      bound_2 = MAXVAL(index_to_cell_R2(1, :))
      bound_3 = MINVAL(index_to_cell_R2(2, :))
      bound_4 = MAXVAL(index_to_cell_R2(2, :))
      bound_5 = MINVAL(index_to_cell_R2(3, :))
      bound_6 = MAXVAL(index_to_cell_R2(3, :))

      ALLOCATE (cell_to_index_R2(bound_1:bound_2, bound_3:bound_4, bound_5:bound_6))
      cell_to_index_R2(:, :, :) = 0

      DO icell = 1, SIZE(index_to_cell_R2, 2)

         cell_to_index_R2(index_to_cell_R2(1, icell), &
                          index_to_cell_R2(2, icell), &
                          index_to_cell_R2(3, icell)) = icell

      END DO

      CALL timestop(handle)

   END SUBROUTINE

END MODULE rpa_gw_kpoints
