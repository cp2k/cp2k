!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for calculating RI-MP2 gradients
!> \par History
!>      10.2013 created [Mauro Del Ben]
! **************************************************************************************************
MODULE mp2_ri_grad_util
!
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type,&
                                              dbcsr_type,&
                                              dbcsr_type_no_symmetry
   USE cp_dbcsr_operations,             ONLY: copy_fm_to_dbcsr,&
                                              cp_dbcsr_m_by_n_from_template
   USE cp_fm_basic_linalg,              ONLY: cp_fm_geadd
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_get,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE group_dist_types,                ONLY: create_group_dist,&
                                              get_group_dist,&
                                              group_dist_d1_type,&
                                              release_group_dist
   USE kinds,                           ONLY: dp
   USE libint_2c_3c,                    ONLY: compare_potential_types
   USE message_passing,                 ONLY: mp_para_env_type,&
                                              mp_request_null,&
                                              mp_request_type,&
                                              mp_waitall
   USE mp2_types,                       ONLY: integ_mat_buffer_type,&
                                              mp2_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE util,                            ONLY: get_limit
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mp2_ri_grad_util'

   PUBLIC :: complete_gamma, array2fm, fm2array, prepare_redistribution, create_dbcsr_gamma

   TYPE index_map
      INTEGER, DIMENSION(:, :), ALLOCATABLE :: map
   END TYPE

CONTAINS

! **************************************************************************************************
!> \brief complete the calculation of the Gamma matrices
!> \param mp2_env ...
!> \param B_ia_Q ...
!> \param dimen_RI ...
!> \param homo ...
!> \param virtual ...
!> \param para_env ...
!> \param para_env_sub ...
!> \param ngroup ...
!> \param my_group_L_size ...
!> \param my_group_L_start ...
!> \param my_group_L_end ...
!> \param my_B_size ...
!> \param my_B_virtual_start ...
!> \param gd_array ...
!> \param gd_B_virtual ...
!> \param kspin ...
!> \author Mauro Del Ben
! **************************************************************************************************
   SUBROUTINE complete_gamma(mp2_env, B_ia_Q, dimen_RI, homo, virtual, para_env, para_env_sub, ngroup, &
                             my_group_L_size, my_group_L_start, my_group_L_end, &
                             my_B_size, my_B_virtual_start, gd_array, gd_B_virtual, kspin)

      TYPE(mp2_type)                                     :: mp2_env
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(INOUT)                                   :: B_ia_Q
      INTEGER, INTENT(IN)                                :: dimen_RI, homo, virtual
      TYPE(mp_para_env_type), INTENT(IN), POINTER        :: para_env, para_env_sub
      INTEGER, INTENT(IN)                                :: ngroup, my_group_L_size, &
                                                            my_group_L_start, my_group_L_end, &
                                                            my_B_size, my_B_virtual_start
      TYPE(group_dist_d1_type), INTENT(IN)               :: gd_array, gd_B_virtual
      INTEGER, INTENT(IN)                                :: kspin

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'complete_gamma'

      INTEGER :: dimen_ia, handle, i, ispin, kkB, my_ia_end, my_ia_size, my_ia_start, my_P_end, &
         my_P_size, my_P_start, nspins, pos_group, pos_sub
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: pos_info
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: group_grid_2_mepos, mepos_2_grid_group
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: BIb_C_2D, Gamma_2D, Gamma_PQ
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_ia, fm_struct_RI
      TYPE(cp_fm_type) :: fm_Gamma, fm_Gamma_PQ, fm_Gamma_PQ_2, fm_Gamma_PQ_temp, &
         fm_Gamma_PQ_temp_2, fm_ia_P, fm_Y, operator_half, PQ_half
      TYPE(group_dist_d1_type)                           :: gd_array_new, gd_ia, gd_ia_new, gd_P, &
                                                            gd_P_new

      CALL timeset(routineN, handle)

      ! reshape the data into a 2D array
      ! reorder the local data in such a way to help the next stage of matrix creation
      ! now the data inside the group are divided into a ia x K matrix
      dimen_ia = homo*virtual
      CALL create_group_dist(gd_ia, para_env_sub%num_pe, dimen_ia)
      CALL get_group_dist(gd_ia, para_env_sub%mepos, my_ia_start, my_ia_end, my_ia_size)

      CALL mat_3d_to_2d(B_ia_Q, BIb_C_2D, homo, my_B_size, virtual, my_B_virtual_start, &
                        my_ia_start, my_ia_end, my_ia_size, my_group_L_size, para_env_sub, gd_B_virtual)

      CALL mat_3d_to_2d_gamma(mp2_env%ri_grad%Gamma_P_ia(kspin)%array, Gamma_2D, homo, my_B_size, virtual, my_B_virtual_start, &
                              my_ia_start, my_ia_end, my_ia_size, my_group_L_size, para_env_sub, gd_B_virtual)

      ! create the processor map and size arrays
      CALL create_group_dist(gd_ia_new, para_env%num_pe)
      CALL create_group_dist(gd_array_new, para_env%num_pe)
      CALL create_group_dist(gd_P, para_env_sub%num_pe, dimen_RI)
      CALL create_group_dist(gd_P_new, para_env%num_pe)

      CALL get_group_dist(gd_P, para_env_sub%mepos, my_P_start, my_P_end, my_P_size)

      CALL prepare_redistribution(para_env, para_env_sub, ngroup, &
                                  group_grid_2_mepos, mepos_2_grid_group, pos_info=pos_info)

      DO i = 0, para_env%num_pe - 1
         ! calculate position of the group
         pos_group = i/para_env_sub%num_pe
         ! calculate position in the subgroup
         pos_sub = pos_info(i)
         ! 1 -> rows, 2 -> cols
         CALL get_group_dist(gd_ia, pos_sub, gd_ia_new, i)
         CALL get_group_dist(gd_array, pos_group, gd_array_new, i)
         CALL get_group_dist(gd_P, pos_sub, gd_P_new, i)
      END DO

      ! create the blacs env
      NULLIFY (blacs_env)
      CALL cp_blacs_env_create(blacs_env=blacs_env, para_env=para_env)

      NULLIFY (fm_struct_ia)
      CALL cp_fm_struct_create(fm_struct_ia, context=blacs_env, nrow_global=dimen_ia, &
                               ncol_global=dimen_RI, para_env=para_env)

      ! create the fm matrix Gamma
      CALL array2fm(Gamma_2D, fm_struct_ia, my_ia_start, my_ia_end, &
                    my_group_L_start, my_group_L_end, &
                    gd_ia_new, gd_array_new, &
                    group_grid_2_mepos, para_env_sub%num_pe, ngroup, &
                    fm_Y)
      ! create the fm matrix B_ia_P
      CALL array2fm(BIb_C_2D, fm_struct_ia, my_ia_start, my_ia_end, &
                    my_group_L_start, my_group_L_end, &
                    gd_ia_new, gd_array_new, &
                    group_grid_2_mepos, para_env_sub%num_pe, ngroup, &
                    fm_ia_P)

      NULLIFY (fm_struct_RI)
      CALL cp_fm_struct_create(fm_struct_RI, context=blacs_env, nrow_global=dimen_RI, &
                               ncol_global=dimen_RI, para_env=para_env)

      ! since we will need (P|Q)^(-1/2) in the future, make a copy
      CALL array2fm(mp2_env%ri_grad%PQ_half, fm_struct_RI, my_P_start, my_P_end, &
                    my_group_L_start, my_group_L_end, &
                    gd_P_new, gd_array_new, &
                    group_grid_2_mepos, para_env_sub%num_pe, ngroup, &
                    PQ_half, do_release_mat=.FALSE.)

      IF (.NOT. compare_potential_types(mp2_env%ri_metric, mp2_env%potential_parameter)) THEN
         CALL array2fm(mp2_env%ri_grad%operator_half, fm_struct_RI, my_P_start, my_P_end, &
                       my_group_L_start, my_group_L_end, &
                       gd_P_new, gd_array_new, &
                       group_grid_2_mepos, para_env_sub%num_pe, ngroup, &
                       operator_half, do_release_mat=.FALSE.)
      END IF

      CALL release_group_dist(gd_P_new)
      CALL release_group_dist(gd_ia_new)
      CALL release_group_dist(gd_array_new)

      ! complete the gamma matrix Gamma_ia^P = sum_Q (Y_ia^Q * (Q|P)^(-1/2))
      IF (compare_potential_types(mp2_env%ri_metric, mp2_env%potential_parameter)) THEN
         CALL cp_fm_create(fm_Gamma, fm_struct_ia, name="fm_Gamma")
         CALL cp_fm_struct_release(fm_struct_ia)
         ! perform the matrix multiplication
         CALL parallel_gemm(transa="N", transb="T", m=dimen_ia, n=dimen_RI, k=dimen_RI, alpha=1.0_dp, &
                            matrix_a=fm_Y, matrix_b=PQ_half, beta=0.0_dp, &
                            matrix_c=fm_Gamma)
         ! release the Y matrix
         CALL cp_fm_release(fm_Y)

         ! complete gamma small (fm_Gamma_PQ)
         ! create temp matrix
         CALL cp_fm_create(fm_Gamma_PQ_temp, fm_struct_RI, name="fm_Gamma_PQ_temp")
         CALL parallel_gemm(transa="T", transb="N", m=dimen_RI, n=dimen_RI, k=dimen_ia, alpha=1.0_dp, &
                            matrix_a=fm_Gamma, matrix_b=fm_ia_P, beta=0.0_dp, &
                            matrix_c=fm_Gamma_PQ_temp)
         CALL cp_fm_release(fm_ia_P)
         ! create fm_Gamma_PQ matrix
         CALL cp_fm_create(fm_Gamma_PQ, fm_struct_RI, name="fm_Gamma_PQ")
         CALL cp_fm_struct_release(fm_struct_RI)
         ! perform matrix multiplication
         CALL parallel_gemm(transa="N", transb="T", m=dimen_RI, n=dimen_RI, k=dimen_RI, alpha=1.0_dp, &
                            matrix_a=fm_Gamma_PQ_temp, matrix_b=PQ_half, beta=0.0_dp, &
                            matrix_c=fm_Gamma_PQ)
         CALL cp_fm_release(fm_Gamma_PQ_temp)
         CALL cp_fm_release(PQ_half)

      ELSE
         ! create temp matrix
         CALL cp_fm_create(fm_Gamma_PQ_temp, fm_struct_RI, name="fm_Gamma_PQ_temp")
         CALL parallel_gemm(transa="T", transb="N", m=dimen_RI, n=dimen_RI, k=dimen_ia, alpha=1.0_dp, &
                            matrix_a=fm_Y, matrix_b=fm_ia_P, beta=0.0_dp, &
                            matrix_c=fm_Gamma_PQ_temp)
         CALL cp_fm_release(fm_ia_P)

         CALL cp_fm_create(fm_Gamma, fm_struct_ia, name="fm_Gamma")
         CALL cp_fm_struct_release(fm_struct_ia)
         ! perform the matrix multiplication
         CALL parallel_gemm(transa="N", transb="T", m=dimen_ia, n=dimen_RI, k=dimen_RI, alpha=1.0_dp, &
                            matrix_a=fm_Y, matrix_b=PQ_half, beta=0.0_dp, &
                            matrix_c=fm_Gamma)
         CALL cp_fm_release(fm_Y)

         CALL cp_fm_create(fm_Gamma_PQ_temp_2, fm_struct_RI, name="fm_Gamma_PQ_temp_2")
         CALL parallel_gemm(transa="N", transb="T", m=dimen_RI, n=dimen_RI, k=dimen_RI, alpha=1.0_dp, &
                            matrix_a=fm_Gamma_PQ_temp, matrix_b=operator_half, beta=0.0_dp, &
                            matrix_c=fm_Gamma_PQ_temp_2)

         CALL cp_fm_create(fm_Gamma_PQ_2, fm_struct_RI, name="fm_Gamma_PQ_2")
         CALL parallel_gemm(transa="N", transb="N", m=dimen_RI, n=dimen_RI, k=dimen_RI, alpha=1.0_dp, &
                            matrix_a=PQ_half, matrix_b=fm_Gamma_PQ_temp_2, beta=0.0_dp, &
                            matrix_c=fm_Gamma_PQ_temp)
         CALL cp_fm_to_fm(fm_Gamma_PQ_temp, fm_Gamma_PQ_2)
         CALL cp_fm_geadd(1.0_dp, "T", fm_Gamma_PQ_temp, 1.0_dp, fm_Gamma_PQ_2)
         CALL cp_fm_release(fm_Gamma_PQ_temp)
         CALL cp_fm_release(PQ_half)

         CALL cp_fm_create(fm_Gamma_PQ, fm_struct_RI)
         CALL cp_fm_struct_release(fm_struct_RI)
         CALL parallel_gemm(transa="N", transb="N", m=dimen_RI, n=dimen_RI, k=dimen_RI, alpha=-1.0_dp, &
                            matrix_a=operator_half, matrix_b=fm_Gamma_PQ_temp_2, beta=0.0_dp, &
                            matrix_c=fm_Gamma_PQ)
         CALL cp_fm_release(operator_half)
         CALL cp_fm_release(fm_Gamma_PQ_temp_2)
      END IF

      ! now redistribute the data, in this case we go back
      ! to the original way the integrals were distributed
      CALL fm2array(Gamma_2D, my_ia_size, my_ia_start, my_ia_end, &
                    my_group_L_size, my_group_L_start, my_group_L_end, &
                    group_grid_2_mepos, mepos_2_grid_group, &
                    para_env_sub%num_pe, ngroup, &
                    fm_Gamma)

      ALLOCATE (Gamma_PQ(my_P_size, my_group_L_size))
      CALL fm2array(Gamma_PQ, my_P_size, my_P_start, my_P_end, &
                    my_group_L_size, my_group_L_start, my_group_L_end, &
                    group_grid_2_mepos, mepos_2_grid_group, &
                    para_env_sub%num_pe, ngroup, &
                    fm_Gamma_PQ)
      IF (.NOT. ALLOCATED(mp2_env%ri_grad%Gamma_PQ)) THEN
         CALL MOVE_ALLOC(Gamma_PQ, mp2_env%ri_grad%Gamma_PQ)
      ELSE
         mp2_env%ri_grad%Gamma_PQ(:, :) = mp2_env%ri_grad%Gamma_PQ + Gamma_PQ
         DEALLOCATE (Gamma_PQ)
      END IF

      IF (.NOT. compare_potential_types(mp2_env%ri_metric, mp2_env%potential_parameter)) THEN
         ALLOCATE (Gamma_PQ(my_P_size, my_group_L_size))
         CALL fm2array(Gamma_PQ, my_P_size, my_P_start, my_P_end, &
                       my_group_L_size, my_group_L_start, my_group_L_end, &
                       group_grid_2_mepos, mepos_2_grid_group, &
                       para_env_sub%num_pe, ngroup, &
                       fm_Gamma_PQ_2)
         IF (.NOT. ALLOCATED(mp2_env%ri_grad%Gamma_PQ_2)) THEN
            CALL MOVE_ALLOC(Gamma_PQ, mp2_env%ri_grad%Gamma_PQ_2)
         ELSE
            mp2_env%ri_grad%Gamma_PQ_2(:, :) = mp2_env%ri_grad%Gamma_PQ_2 + Gamma_PQ
            DEALLOCATE (Gamma_PQ)
         END IF
      END IF

      ! allocate G_P_ia (DBCSR)
      IF (.NOT. ALLOCATED(mp2_env%ri_grad%G_P_ia)) THEN
         nspins = SIZE(mp2_env%ri_grad%mo_coeff_o)
         ALLOCATE (mp2_env%ri_grad%G_P_ia(my_group_L_size, nspins))
         DO ispin = 1, nspins
         DO kkB = 1, my_group_L_size
            NULLIFY (mp2_env%ri_grad%G_P_ia(kkB, ispin)%matrix)
         END DO
         END DO
      END IF

      ! create the Gamma_ia_P in DBCSR style
      CALL create_dbcsr_gamma(Gamma_2D, homo, virtual, dimen_ia, para_env_sub, &
                              my_ia_start, my_ia_end, my_group_L_size, gd_ia, &
                              mp2_env%ri_grad%G_P_ia(:, kspin), mp2_env%ri_grad%mo_coeff_o(1)%matrix)

      DEALLOCATE (pos_info)
      DEALLOCATE (group_grid_2_mepos, mepos_2_grid_group)
      CALL release_group_dist(gd_ia)
      CALL release_group_dist(gd_P)

      ! release blacs_env
      CALL cp_blacs_env_release(blacs_env)

      CALL timestop(handle)

   END SUBROUTINE complete_gamma

! **************************************************************************************************
!> \brief Redistribute a 3d matrix to a 2d matrix
!> \param B_ia_Q ...
!> \param BIb_C_2D ...
!> \param homo ...
!> \param my_B_size ...
!> \param virtual ...
!> \param my_B_virtual_start ...
!> \param my_ia_start ...
!> \param my_ia_end ...
!> \param my_ia_size ...
!> \param my_group_L_size ...
!> \param para_env_sub ...
!> \param gd_B_virtual ...
! **************************************************************************************************
   SUBROUTINE mat_3d_to_2d(B_ia_Q, BIb_C_2D, homo, my_B_size, virtual, my_B_virtual_start, &
                           my_ia_start, my_ia_end, my_ia_size, my_group_L_size, para_env_sub, gd_B_virtual)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(INOUT)                                   :: B_ia_Q
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: BIb_C_2D
      INTEGER, INTENT(IN)                                :: homo, my_B_size, virtual, &
                                                            my_B_virtual_start, my_ia_start, &
                                                            my_ia_end, my_ia_size, my_group_L_size
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env_sub
      TYPE(group_dist_d1_type), INTENT(IN)               :: gd_B_virtual

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'mat_3d_to_2d'

      INTEGER :: handle, ia_global, iiB, jjB, proc_receive, proc_send, proc_shift, rec_B_size, &
         rec_B_virtual_end, rec_B_virtual_start
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: BIb_C_rec

      CALL timeset(routineN, handle)

      ALLOCATE (BIb_C_2D(my_ia_size, my_group_L_size))
      BIb_C_2D = 0.0_dp

      DO iiB = 1, homo
         DO jjB = 1, my_B_size
            ia_global = (iiB - 1)*virtual + my_B_virtual_start + jjB - 1
            IF (ia_global >= my_ia_start .AND. ia_global <= my_ia_end) THEN
               BIb_C_2D(ia_global - my_ia_start + 1, :) = B_ia_Q(iiB, jjB, 1:my_group_L_size)
            END IF
         END DO
      END DO

      DO proc_shift = 1, para_env_sub%num_pe - 1
         proc_send = MODULO(para_env_sub%mepos + proc_shift, para_env_sub%num_pe)
         proc_receive = MODULO(para_env_sub%mepos - proc_shift, para_env_sub%num_pe)

         CALL get_group_dist(gd_B_virtual, proc_receive, rec_B_virtual_start, rec_B_virtual_end, rec_B_size)

         ALLOCATE (BIb_C_rec(homo, rec_B_size, my_group_L_size))
         BIb_C_rec = 0.0_dp

         CALL para_env_sub%sendrecv(B_ia_Q, proc_send, &
                                    BIb_C_rec, proc_receive)

         DO iiB = 1, homo
            DO jjB = 1, rec_B_size
               ia_global = (iiB - 1)*virtual + rec_B_virtual_start + jjB - 1
               IF (ia_global >= my_ia_start .AND. ia_global <= my_ia_end) THEN
                  BIb_C_2D(ia_global - my_ia_start + 1, :) = BIb_C_rec(iiB, jjB, 1:my_group_L_size)
               END IF
            END DO
         END DO

         DEALLOCATE (BIb_C_rec)
      END DO
      DEALLOCATE (B_ia_Q)

      CALL timestop(handle)
   END SUBROUTINE mat_3d_to_2d

! **************************************************************************************************
!> \brief Redistribute a 3d matrix to a 2d matrix, specialized for Gamma_P_ia to account for a different order of indices
!> \param B_ia_Q ...
!> \param BIb_C_2D ...
!> \param homo ...
!> \param my_B_size ...
!> \param virtual ...
!> \param my_B_virtual_start ...
!> \param my_ia_start ...
!> \param my_ia_end ...
!> \param my_ia_size ...
!> \param my_group_L_size ...
!> \param para_env_sub ...
!> \param gd_B_virtual ...
! **************************************************************************************************
   SUBROUTINE mat_3d_to_2d_gamma(B_ia_Q, BIb_C_2D, homo, my_B_size, virtual, my_B_virtual_start, &
                                 my_ia_start, my_ia_end, my_ia_size, my_group_L_size, para_env_sub, gd_B_virtual)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(INOUT)                                   :: B_ia_Q
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: BIb_C_2D
      INTEGER, INTENT(IN)                                :: homo, my_B_size, virtual, &
                                                            my_B_virtual_start, my_ia_start, &
                                                            my_ia_end, my_ia_size, my_group_L_size
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env_sub
      TYPE(group_dist_d1_type), INTENT(IN)               :: gd_B_virtual

      CHARACTER(LEN=*), PARAMETER :: routineN = 'mat_3d_to_2d_gamma'

      INTEGER :: handle, ia_global, iiB, jjB, proc_receive, proc_send, proc_shift, rec_B_size, &
         rec_B_virtual_end, rec_B_virtual_start
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: BIb_C_rec

      CALL timeset(routineN, handle)

      ALLOCATE (BIb_C_2D(my_ia_size, my_group_L_size))
      BIb_C_2D = 0.0_dp

      DO iiB = 1, homo
         DO jjB = 1, my_B_size
            ia_global = (iiB - 1)*virtual + my_B_virtual_start + jjB - 1
            IF (ia_global >= my_ia_start .AND. ia_global <= my_ia_end) THEN
               BIb_C_2D(ia_global - my_ia_start + 1, :) = B_ia_Q(jjB, iiB, 1:my_group_L_size)
            END IF
         END DO
      END DO

      DO proc_shift = 1, para_env_sub%num_pe - 1
         proc_send = MODULO(para_env_sub%mepos + proc_shift, para_env_sub%num_pe)
         proc_receive = MODULO(para_env_sub%mepos - proc_shift, para_env_sub%num_pe)

         CALL get_group_dist(gd_B_virtual, proc_receive, rec_B_virtual_start, rec_B_virtual_end, rec_B_size)

         ALLOCATE (BIb_C_rec(rec_B_size, homo, my_group_L_size))
         BIb_C_rec = 0.0_dp

         CALL para_env_sub%sendrecv(B_ia_Q, proc_send, &
                                    BIb_C_rec, proc_receive)

         DO iiB = 1, homo
            DO jjB = 1, rec_B_size
               ia_global = (iiB - 1)*virtual + rec_B_virtual_start + jjB - 1
               IF (ia_global >= my_ia_start .AND. ia_global <= my_ia_end) THEN
                  BIb_C_2D(ia_global - my_ia_start + 1, :) = BIb_C_rec(jjB, iiB, 1:my_group_L_size)
               END IF
            END DO
         END DO

         DEALLOCATE (BIb_C_rec)
      END DO
      DEALLOCATE (B_ia_Q)

      CALL timestop(handle)
   END SUBROUTINE mat_3d_to_2d_gamma

! **************************************************************************************************
!> \brief redistribute local part of array to fm
!> \param mat2D ...
!> \param fm_struct ...
!> \param my_start_row ...
!> \param my_end_row ...
!> \param my_start_col ...
!> \param my_end_col ...
!> \param gd_row ...
!> \param gd_col ...
!> \param group_grid_2_mepos ...
!> \param ngroup_row ...
!> \param ngroup_col ...
!> \param fm_mat ...
!> \param integ_group_size ...
!> \param color_group ...
!> \param do_release_mat whether to release the array (default: yes)
! **************************************************************************************************
   SUBROUTINE array2fm(mat2D, fm_struct, my_start_row, my_end_row, &
                       my_start_col, my_end_col, &
                       gd_row, gd_col, &
                       group_grid_2_mepos, ngroup_row, ngroup_col, &
                       fm_mat, integ_group_size, color_group, do_release_mat)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: mat2D
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      INTEGER, INTENT(IN)                                :: my_start_row, my_end_row, my_start_col, &
                                                            my_end_col
      TYPE(group_dist_d1_type), INTENT(IN)               :: gd_row, gd_col
      INTEGER, ALLOCATABLE, DIMENSION(:, :), INTENT(IN)  :: group_grid_2_mepos
      INTEGER, INTENT(IN)                                :: ngroup_row, ngroup_col
      TYPE(cp_fm_type), INTENT(OUT)                      :: fm_mat
      INTEGER, INTENT(IN), OPTIONAL                      :: integ_group_size, color_group
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_release_mat

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'array2fm'

      INTEGER :: dummy_proc, end_col_block, end_row_block, handle, handle2, i_global, i_local, &
         i_sub, iiB, iii, itmp(2), j_global, j_local, j_sub, jjB, my_num_col_blocks, &
         my_num_row_blocks, mypcol, myprow, ncol_local, npcol, nprow, nrow_local, num_cols, &
         num_rec_cols, num_rows, number_of_rec, number_of_send, proc_receive, proc_send, &
         proc_shift, rec_col_end, rec_col_size, rec_col_start, rec_counter, rec_pcol, rec_prow, &
         rec_row_end, rec_row_size, rec_row_start, ref_send_pcol, ref_send_prow, send_counter, &
         send_pcol, send_prow, size_rec_buffer, size_send_buffer, start_col_block, start_row_block
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: iii_vet, index_col_rec, map_rec_size, &
                                                            map_send_size
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: blocks_ranges_col, blocks_ranges_row, &
                                                            grid_2_mepos, grid_ref_2_send_pos, &
                                                            mepos_2_grid
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: convert_pos, my_do_release_mat
      REAL(KIND=dp)                                      :: part_col, part_row
      TYPE(integ_mat_buffer_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: buffer_rec, buffer_send
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(mp_request_type), ALLOCATABLE, DIMENSION(:)   :: req_send

      CALL timeset(routineN, handle)

      my_do_release_mat = .TRUE.
      IF (PRESENT(do_release_mat)) my_do_release_mat = do_release_mat

      CALL cp_fm_struct_get(fm_struct, para_env=para_env, nrow_global=num_rows, ncol_global=num_cols)

      ! create the full matrix, (num_rows x num_cols)
      CALL cp_fm_create(fm_mat, fm_struct, name="fm_mat")
      CALL cp_fm_set_all(matrix=fm_mat, alpha=0.0_dp)

      ! start filling procedure
      ! fill the matrix
      CALL cp_fm_get_info(matrix=fm_mat, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)
      myprow = fm_mat%matrix_struct%context%mepos(1)
      mypcol = fm_mat%matrix_struct%context%mepos(2)
      nprow = fm_mat%matrix_struct%context%num_pe(1)
      npcol = fm_mat%matrix_struct%context%num_pe(2)

      ! start the communication part
      ! 0) create array containing the processes position
      !    and supporting infos
      CALL timeset(routineN//"_info", handle2)
      ALLOCATE (grid_2_mepos(0:nprow - 1, 0:npcol - 1))
      grid_2_mepos = 0
      ALLOCATE (mepos_2_grid(2, 0:para_env%num_pe - 1))
      grid_2_mepos(myprow, mypcol) = para_env%mepos
      ! sum infos
      CALL para_env%sum(grid_2_mepos)
      CALL para_env%allgather([myprow, mypcol], mepos_2_grid)

      ! 1) loop over my local data and define a map for the proc to send data
      ALLOCATE (map_send_size(0:para_env%num_pe - 1))
      map_send_size = 0
      dummy_proc = 0
      DO jjB = my_start_col, my_end_col
         send_pcol = fm_mat%matrix_struct%g2p_col(jjB)
         DO iiB = my_start_row, my_end_row
            send_prow = fm_mat%matrix_struct%g2p_row(iiB)
            proc_send = grid_2_mepos(send_prow, send_pcol)
            map_send_size(proc_send) = map_send_size(proc_send) + 1
         END DO
      END DO

      ! 2) loop over my local data of fm_mat and define a map for the proc from which rec data
      ALLOCATE (map_rec_size(0:para_env%num_pe - 1))
      map_rec_size = 0
      part_row = REAL(num_rows, KIND=dp)/REAL(ngroup_row, KIND=dp)
      part_col = REAL(num_cols, KIND=dp)/REAL(ngroup_col, KIND=dp)

      ! In some cases we have to convert global positions to positions in a subgroup (RPA)
      convert_pos = .FALSE.
      IF (PRESENT(integ_group_size) .AND. PRESENT(color_group)) convert_pos = .TRUE.

      DO jjB = 1, nrow_local
         j_global = row_indices(jjB)
         ! check the group holding this element
         ! dirty way, if someone has a better idea ...
         rec_prow = INT(REAL(j_global - 1, KIND=dp)/part_row)
         rec_prow = MAX(0, rec_prow)
         rec_prow = MIN(rec_prow, ngroup_row - 1)
         DO
            itmp = get_limit(num_rows, ngroup_row, rec_prow)
            IF (j_global >= itmp(1) .AND. j_global <= itmp(2)) EXIT
            IF (j_global < itmp(1)) rec_prow = rec_prow - 1
            IF (j_global > itmp(2)) rec_prow = rec_prow + 1
         END DO

         IF (convert_pos) THEN
            ! if the group is not in the same RPA group cycle
            IF ((rec_prow/integ_group_size) .NE. color_group) CYCLE
            ! convert global position to position into the RPA group
            rec_prow = MOD(rec_prow, integ_group_size)
         END IF

         DO iiB = 1, ncol_local
            i_global = col_indices(iiB)
            ! check the process in the group holding this element
            rec_pcol = INT(REAL(i_global - 1, KIND=dp)/part_col)
            rec_pcol = MAX(0, rec_pcol)
            rec_pcol = MIN(rec_pcol, ngroup_col - 1)
            DO
               itmp = get_limit(num_cols, ngroup_col, rec_pcol)
               IF (i_global >= itmp(1) .AND. i_global <= itmp(2)) EXIT
               IF (i_global < itmp(1)) rec_pcol = rec_pcol - 1
               IF (i_global > itmp(2)) rec_pcol = rec_pcol + 1
            END DO

            proc_receive = group_grid_2_mepos(rec_prow, rec_pcol)

            map_rec_size(proc_receive) = map_rec_size(proc_receive) + 1

         END DO ! i_global
      END DO ! j_global

      ! 3) check if the local data has to be stored in the new fm_mat
      IF (map_rec_size(para_env%mepos) > 0) THEN
         DO jjB = 1, ncol_local
            j_global = col_indices(jjB)
            IF (j_global >= my_start_col .AND. j_global <= my_end_col) THEN
               DO iiB = 1, nrow_local
                  i_global = row_indices(iiB)
                  IF (i_global >= my_start_row .AND. i_global <= my_end_row) THEN
                     fm_mat%local_data(iiB, jjB) = mat2D(i_global - my_start_row + 1, j_global - my_start_col + 1)
                  END IF
               END DO
            END IF
         END DO
      END IF
      CALL timestop(handle2)

      ! 4) reorder data in the send_buffer
      CALL timeset(routineN//"_buffer_send", handle2)
      ! count the number of messages to send
      number_of_send = 0
      DO proc_shift = 1, para_env%num_pe - 1
         proc_send = MODULO(para_env%mepos + proc_shift, para_env%num_pe)
         IF (map_send_size(proc_send) > 0) THEN
            number_of_send = number_of_send + 1
         END IF
      END DO
      ! allocate the structure that will hold the messages to be sent
      ALLOCATE (buffer_send(number_of_send))
      ! grid_ref_2_send_pos is an array (map) that given a pair
      ! (ref_send_prow,ref_send_pcol) returns
      ! the position in the buffer_send associated to that process
      ALLOCATE (grid_ref_2_send_pos(0:nprow - 1, 0:npcol - 1))
      grid_ref_2_send_pos = 0
      ! finalize the allocation of buffer_send with the actual size
      ! of each message (actual size is size_send_buffer)
      send_counter = 0
      DO proc_shift = 1, para_env%num_pe - 1
         proc_send = MODULO(para_env%mepos + proc_shift, para_env%num_pe)
         size_send_buffer = map_send_size(proc_send)
         IF (map_send_size(proc_send) > 0) THEN
            send_counter = send_counter + 1
            ! allocate the sending buffer (msg)
            ALLOCATE (buffer_send(send_counter)%msg(size_send_buffer))
            buffer_send(send_counter)%msg = 0.0_dp
            buffer_send(send_counter)%proc = proc_send
            ! get the pointer to prow, pcol of the process that has
            ! to receive this message
            ref_send_prow = mepos_2_grid(1, proc_send)
            ref_send_pcol = mepos_2_grid(2, proc_send)
            ! save the rank of the process that has to receive this message
            grid_ref_2_send_pos(ref_send_prow, ref_send_pcol) = send_counter
         END IF
      END DO

      ! loop over the locally held data and fill the buffer_send
      ! for doing that we need an array that keep track if the
      ! sequential increase of the index for each message
      ALLOCATE (iii_vet(number_of_send))
      iii_vet = 0
      DO iiB = my_start_row, my_end_row
         send_prow = fm_mat%matrix_struct%g2p_row(iiB)
         DO jjB = my_start_col, my_end_col
            send_pcol = fm_mat%matrix_struct%g2p_col(jjB)
            ! we don't need to send to ourselves
            IF (grid_2_mepos(send_prow, send_pcol) == para_env%mepos) CYCLE

            send_counter = grid_ref_2_send_pos(send_prow, send_pcol)
            iii_vet(send_counter) = iii_vet(send_counter) + 1
            iii = iii_vet(send_counter)
            buffer_send(send_counter)%msg(iii) = mat2D(iiB - my_start_row + 1, jjB - my_start_col + 1)
         END DO
      END DO

      DEALLOCATE (iii_vet)
      DEALLOCATE (grid_ref_2_send_pos)
      IF (my_do_release_mat) DEALLOCATE (mat2D)
      CALL timestop(handle2)

      ! 5) similarly to what done for the buffer_send
      !    create the buffer for receive, post the message with irecv
      !    and send the messages non-blocking
      CALL timeset(routineN//"_isendrecv", handle2)
      ! count the number of messages to be received
      number_of_rec = 0
      DO proc_shift = 1, para_env%num_pe - 1
         proc_receive = MODULO(para_env%mepos - proc_shift, para_env%num_pe)
         IF (map_rec_size(proc_receive) > 0) THEN
            number_of_rec = number_of_rec + 1
         END IF
      END DO

      ALLOCATE (buffer_rec(number_of_rec))

      rec_counter = 0
      DO proc_shift = 1, para_env%num_pe - 1
         proc_receive = MODULO(para_env%mepos - proc_shift, para_env%num_pe)
         size_rec_buffer = map_rec_size(proc_receive)
         IF (map_rec_size(proc_receive) > 0) THEN
            rec_counter = rec_counter + 1
            ! prepare the buffer for receive
            ALLOCATE (buffer_rec(rec_counter)%msg(size_rec_buffer))
            buffer_rec(rec_counter)%msg = 0.0_dp
            buffer_rec(rec_counter)%proc = proc_receive
            ! post the message to be received
            CALL para_env%irecv(buffer_rec(rec_counter)%msg, proc_receive, &
                                buffer_rec(rec_counter)%msg_req)
         END IF
      END DO

      ! send messages
      ALLOCATE (req_send(number_of_send))
      send_counter = 0
      DO proc_shift = 1, para_env%num_pe - 1
         proc_send = MODULO(para_env%mepos + proc_shift, para_env%num_pe)
         IF (map_send_size(proc_send) > 0) THEN
            send_counter = send_counter + 1
            CALL para_env%isend(buffer_send(send_counter)%msg, proc_send, &
                                buffer_send(send_counter)%msg_req)
            req_send(send_counter) = buffer_send(send_counter)%msg_req
         END IF
      END DO
      CALL timestop(handle2)

      ! 6) fill the fm_mat matrix with the messages received from the
      !    other processes
      CALL timeset(routineN//"_fill", handle2)
      ! In order to perform this step efficiently first we have to know the
      ! ranges of the blocks over which a given process hold its local data.
      ! Start with the rows ...
      my_num_row_blocks = 1
      DO iiB = 1, nrow_local - 1
         IF (ABS(row_indices(iiB + 1) - row_indices(iiB)) == 1) CYCLE
         my_num_row_blocks = my_num_row_blocks + 1
      END DO
      ALLOCATE (blocks_ranges_row(2, my_num_row_blocks))
      blocks_ranges_row = 0
      blocks_ranges_row(1, 1) = row_indices(1)
      iii = 1
      DO iiB = 1, nrow_local - 1
         IF (ABS(row_indices(iiB + 1) - row_indices(iiB)) == 1) CYCLE
         iii = iii + 1
         blocks_ranges_row(2, iii - 1) = row_indices(iiB)
         blocks_ranges_row(1, iii) = row_indices(iiB + 1)
      END DO
      blocks_ranges_row(2, my_num_row_blocks) = row_indices(MAX(nrow_local, 1))
      ! ... and columns
      my_num_col_blocks = 1
      DO jjB = 1, ncol_local - 1
         IF (ABS(col_indices(jjB + 1) - col_indices(jjB)) == 1) CYCLE
         my_num_col_blocks = my_num_col_blocks + 1
      END DO
      ALLOCATE (blocks_ranges_col(2, my_num_col_blocks))
      blocks_ranges_col = 0
      blocks_ranges_col(1, 1) = col_indices(1)
      iii = 1
      DO jjB = 1, ncol_local - 1
         IF (ABS(col_indices(jjB + 1) - col_indices(jjB)) == 1) CYCLE
         iii = iii + 1
         blocks_ranges_col(2, iii - 1) = col_indices(jjB)
         blocks_ranges_col(1, iii) = col_indices(jjB + 1)
      END DO
      blocks_ranges_col(2, my_num_col_blocks) = col_indices(MAX(ncol_local, 1))

      rec_counter = 0
      DO proc_shift = 1, para_env%num_pe - 1
         proc_receive = MODULO(para_env%mepos - proc_shift, para_env%num_pe)
         size_rec_buffer = map_rec_size(proc_receive)

         IF (map_rec_size(proc_receive) > 0) THEN
            rec_counter = rec_counter + 1

            CALL get_group_dist(gd_col, proc_receive, rec_col_start, rec_col_end, rec_col_size)

            ! precompute the number of received columns and relative index
            num_rec_cols = 0
            DO jjB = 1, my_num_col_blocks
               start_col_block = MAX(blocks_ranges_col(1, jjB), rec_col_start)
               end_col_block = MIN(blocks_ranges_col(2, jjB), rec_col_end)
               DO j_sub = start_col_block, end_col_block
                  num_rec_cols = num_rec_cols + 1
               END DO
            END DO
            ALLOCATE (index_col_rec(num_rec_cols))
            index_col_rec = 0
            iii = 0
            DO jjB = 1, my_num_col_blocks
               start_col_block = MAX(blocks_ranges_col(1, jjB), rec_col_start)
               end_col_block = MIN(blocks_ranges_col(2, jjB), rec_col_end)
               DO j_sub = start_col_block, end_col_block
                  iii = iii + 1
                  j_local = fm_mat%matrix_struct%g2l_col(j_sub)
                  index_col_rec(iii) = j_local
               END DO
            END DO

            CALL get_group_dist(gd_row, proc_receive, rec_row_start, rec_row_end, rec_row_size)

            ! wait for the message
            CALL buffer_rec(rec_counter)%msg_req%wait()

            ! fill the local data in fm_mat
            iii = 0
            DO iiB = 1, my_num_row_blocks
               start_row_block = MAX(blocks_ranges_row(1, iiB), rec_row_start)
               end_row_block = MIN(blocks_ranges_row(2, iiB), rec_row_end)
               DO i_sub = start_row_block, end_row_block
                  i_local = fm_mat%matrix_struct%g2l_row(i_sub)
                  DO jjB = 1, num_rec_cols
                     iii = iii + 1
                     j_local = index_col_rec(jjB)
                     fm_mat%local_data(i_local, j_local) = buffer_rec(rec_counter)%msg(iii)
                  END DO
               END DO
            END DO

            DEALLOCATE (buffer_rec(rec_counter)%msg)
            DEALLOCATE (index_col_rec)
         END IF
      END DO
      DEALLOCATE (buffer_rec)

      DEALLOCATE (blocks_ranges_row)
      DEALLOCATE (blocks_ranges_col)

      CALL timestop(handle2)

      ! 7) Finally wait for all messeges to be sent
      CALL timeset(routineN//"_waitall", handle2)
      CALL mp_waitall(req_send(:))
      DO send_counter = 1, number_of_send
         DEALLOCATE (buffer_send(send_counter)%msg)
      END DO
      DEALLOCATE (buffer_send)
      CALL timestop(handle2)

      DEALLOCATE (map_send_size)
      DEALLOCATE (map_rec_size)
      DEALLOCATE (grid_2_mepos)
      DEALLOCATE (mepos_2_grid)

      CALL timestop(handle)

   END SUBROUTINE array2fm

! **************************************************************************************************
!> \brief redistribute fm to local part of array
!> \param mat2D ...
!> \param my_rows ...
!> \param my_start_row ...
!> \param my_end_row ...
!> \param my_cols ...
!> \param my_start_col ...
!> \param my_end_col ...
!> \param group_grid_2_mepos ...
!> \param mepos_2_grid_group ...
!> \param ngroup_row ...
!> \param ngroup_col ...
!> \param fm_mat ...
! **************************************************************************************************
   SUBROUTINE fm2array(mat2D, my_rows, my_start_row, my_end_row, &
                       my_cols, my_start_col, my_end_col, &
                       group_grid_2_mepos, mepos_2_grid_group, &
                       ngroup_row, ngroup_col, &
                       fm_mat)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: mat2D
      INTEGER, INTENT(IN)                                :: my_rows, my_start_row, my_end_row, &
                                                            my_cols, my_start_col, my_end_col
      INTEGER, ALLOCATABLE, DIMENSION(:, :), INTENT(IN)  :: group_grid_2_mepos, mepos_2_grid_group
      INTEGER, INTENT(IN)                                :: ngroup_row, ngroup_col
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_mat

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'fm2array'

      INTEGER :: dummy_proc, handle, handle2, i_global, iiB, iii, itmp(2), j_global, jjB, mypcol, &
         myprow, ncol_local, npcol, nprow, nrow_local, num_cols, num_rec_rows, num_rows, &
         number_of_rec, number_of_send, proc_receive, proc_send, proc_shift, rec_col_size, &
         rec_counter, rec_pcol, rec_prow, rec_row_size, ref_send_pcol, ref_send_prow, &
         send_counter, send_pcol, send_prow, size_rec_buffer, size_send_buffer
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: iii_vet, index_row_rec, map_rec_size, &
                                                            map_send_size
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: grid_2_mepos, grid_ref_2_send_pos, &
                                                            mepos_2_grid, sizes
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: part_col, part_row
      TYPE(integ_mat_buffer_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: buffer_rec, buffer_send
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(mp_request_type), ALLOCATABLE, DIMENSION(:)   :: req_send

      CALL timeset(routineN, handle)

      ! allocate the array
      ALLOCATE (mat2D(my_rows, my_cols))
      mat2D = 0.0_dp

      ! start procedure
      ! get information of from the full matrix
      CALL cp_fm_get_info(matrix=fm_mat, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices, &
                          nrow_global=num_rows, &
                          ncol_global=num_cols, &
                          para_env=para_env)
      myprow = fm_mat%matrix_struct%context%mepos(1)
      mypcol = fm_mat%matrix_struct%context%mepos(2)
      nprow = fm_mat%matrix_struct%context%num_pe(1)
      npcol = fm_mat%matrix_struct%context%num_pe(2)

      ! start the communication part
      ! 0) create array containing the processes position
      !    and supporting infos
      CALL timeset(routineN//"_info", handle2)
      ALLOCATE (grid_2_mepos(0:nprow - 1, 0:npcol - 1))
      grid_2_mepos = 0
      ALLOCATE (mepos_2_grid(2, 0:para_env%num_pe - 1))

      ! fill the info array
      grid_2_mepos(myprow, mypcol) = para_env%mepos

      ! sum infos
      CALL para_env%sum(grid_2_mepos)
      CALL para_env%allgather([myprow, mypcol], mepos_2_grid)

      ! 1) loop over my local data and define a map for the proc to send data
      ALLOCATE (map_send_size(0:para_env%num_pe - 1))
      map_send_size = 0
      part_row = REAL(num_rows, KIND=dp)/REAL(ngroup_row, KIND=dp)
      part_col = REAL(num_cols, KIND=dp)/REAL(ngroup_col, KIND=dp)

      DO jjB = 1, ncol_local
         j_global = col_indices(jjB)
         ! check the group holding this element
         ! dirty way, if someone has a better idea ...
         send_pcol = INT(REAL(j_global - 1, KIND=dp)/part_col)
         send_pcol = MAX(0, send_pcol)
         send_pcol = MIN(send_pcol, ngroup_col - 1)
         DO
            itmp = get_limit(num_cols, ngroup_col, send_pcol)
            IF (j_global >= itmp(1) .AND. j_global <= itmp(2)) EXIT
            IF (j_global < itmp(1)) send_pcol = send_pcol - 1
            IF (j_global > itmp(2)) send_pcol = send_pcol + 1
         END DO

         DO iiB = 1, nrow_local
            i_global = row_indices(iiB)
            ! check the process in the group holding this element
            send_prow = INT(REAL(i_global - 1, KIND=dp)/part_row)
            send_prow = MAX(0, send_prow)
            send_prow = MIN(send_prow, ngroup_row - 1)
            DO
               itmp = get_limit(num_rows, ngroup_row, send_prow)
               IF (i_global >= itmp(1) .AND. i_global <= itmp(2)) EXIT
               IF (i_global < itmp(1)) send_prow = send_prow - 1
               IF (i_global > itmp(2)) send_prow = send_prow + 1
            END DO

            proc_send = group_grid_2_mepos(send_prow, send_pcol)

            map_send_size(proc_send) = map_send_size(proc_send) + 1

         END DO ! i_global
      END DO ! j_global

      ! 2) loop over my local data of the array and define a map for the proc from which rec data
      ALLOCATE (map_rec_size(0:para_env%num_pe - 1))
      map_rec_size = 0
      dummy_proc = 0
      DO jjB = my_start_col, my_end_col
         rec_pcol = fm_mat%matrix_struct%g2p_col(jjB)
         DO iiB = my_start_row, my_end_row
            rec_prow = fm_mat%matrix_struct%g2p_row(iiB)
            proc_receive = grid_2_mepos(rec_prow, rec_pcol)
            map_rec_size(proc_receive) = map_rec_size(proc_receive) + 1
         END DO
      END DO

      ! 3) check if the local data in the fm_mat has to be stored in the new array
      IF (map_rec_size(para_env%mepos) > 0) THEN
         DO jjB = 1, ncol_local
            j_global = col_indices(jjB)
            IF (j_global >= my_start_col .AND. j_global <= my_end_col) THEN
               DO iiB = 1, nrow_local
                  i_global = row_indices(iiB)
                  IF (i_global >= my_start_row .AND. i_global <= my_end_row) THEN
                     mat2D(i_global - my_start_row + 1, j_global - my_start_col + 1) = fm_mat%local_data(iiB, jjB)
                  END IF
               END DO
            END IF
         END DO
      END IF
      CALL timestop(handle2)

      ! 4) reorder data in the send_buffer
      CALL timeset(routineN//"_buffer_send", handle2)
      ! count the number of messages to send
      number_of_send = 0
      DO proc_shift = 1, para_env%num_pe - 1
         proc_send = MODULO(para_env%mepos + proc_shift, para_env%num_pe)
         IF (map_send_size(proc_send) > 0) THEN
            number_of_send = number_of_send + 1
         END IF
      END DO
      ! allocate the structure that will hold the messages to be sent
      ALLOCATE (buffer_send(number_of_send))
      ! grid_ref_2_send_pos is an array (map) that given a pair
      ! (ref_send_prow,ref_send_pcol) returns
      ! the position in the buffer_send associated to that process

      ALLOCATE (grid_ref_2_send_pos(0:ngroup_row - 1, 0:ngroup_col - 1))
      grid_ref_2_send_pos = 0

      ! finalize the allocation of buffer_send with the actual size
      ! of each message (actual size is size_send_buffer)
      send_counter = 0
      DO proc_shift = 1, para_env%num_pe - 1
         proc_send = MODULO(para_env%mepos + proc_shift, para_env%num_pe)
         size_send_buffer = map_send_size(proc_send)
         IF (map_send_size(proc_send) > 0) THEN
            send_counter = send_counter + 1
            ! allocate the sending buffer (msg)
            ALLOCATE (buffer_send(send_counter)%msg(size_send_buffer))
            buffer_send(send_counter)%msg = 0.0_dp
            buffer_send(send_counter)%proc = proc_send
            ! get the pointer to prow, pcol of the process that has
            ! to receive this message
            ref_send_prow = mepos_2_grid_group(1, proc_send)
            ref_send_pcol = mepos_2_grid_group(2, proc_send)
            ! save the rank of the process that has to receive this message
            grid_ref_2_send_pos(ref_send_prow, ref_send_pcol) = send_counter
         END IF
      END DO

      ! loop over the locally held data and fill the buffer_send
      ! for doing that we need an array that keep track if the
      ! sequential increase of the index for each message
      ALLOCATE (iii_vet(number_of_send))
      iii_vet = 0
      DO jjB = 1, ncol_local
         j_global = col_indices(jjB)
         ! check the group holding this element
         ! dirty way, if someone has a better idea ...
         send_pcol = INT(REAL(j_global - 1, KIND=dp)/part_col)
         send_pcol = MAX(0, send_pcol)
         send_pcol = MIN(send_pcol, ngroup_col - 1)
         DO
            itmp = get_limit(num_cols, ngroup_col, send_pcol)
            IF (j_global >= itmp(1) .AND. j_global <= itmp(2)) EXIT
            IF (j_global < itmp(1)) send_pcol = send_pcol - 1
            IF (j_global > itmp(2)) send_pcol = send_pcol + 1
         END DO

         DO iiB = 1, nrow_local
            i_global = row_indices(iiB)
            ! check the process in the group holding this element
            send_prow = INT(REAL(i_global - 1, KIND=dp)/part_row)
            send_prow = MAX(0, send_prow)
            send_prow = MIN(send_prow, ngroup_row - 1)
            DO
               itmp = get_limit(num_rows, ngroup_row, send_prow)
               IF (i_global >= itmp(1) .AND. i_global <= itmp(2)) EXIT
               IF (i_global < itmp(1)) send_prow = send_prow - 1
               IF (i_global > itmp(2)) send_prow = send_prow + 1
            END DO
            ! we don't need to send to ourselves
            IF (group_grid_2_mepos(send_prow, send_pcol) == para_env%mepos) CYCLE

            send_counter = grid_ref_2_send_pos(send_prow, send_pcol)
            iii_vet(send_counter) = iii_vet(send_counter) + 1
            iii = iii_vet(send_counter)
            buffer_send(send_counter)%msg(iii) = fm_mat%local_data(iiB, jjB)
         END DO
      END DO

      DEALLOCATE (iii_vet)
      DEALLOCATE (grid_ref_2_send_pos)
      CALL timestop(handle2)

      ! 5) similarly to what done for the buffer_send
      !    create the buffer for receive, post the message with irecv
      !    and send the messages with non-blocking
      CALL timeset(routineN//"_isendrecv", handle2)
      ! count the number of messages to be received
      number_of_rec = 0
      DO proc_shift = 1, para_env%num_pe - 1
         proc_receive = MODULO(para_env%mepos - proc_shift, para_env%num_pe)
         IF (map_rec_size(proc_receive) > 0) THEN
            number_of_rec = number_of_rec + 1
         END IF
      END DO

      ALLOCATE (buffer_rec(number_of_rec))

      rec_counter = 0
      DO proc_shift = 1, para_env%num_pe - 1
         proc_receive = MODULO(para_env%mepos - proc_shift, para_env%num_pe)
         size_rec_buffer = map_rec_size(proc_receive)
         IF (map_rec_size(proc_receive) > 0) THEN
            rec_counter = rec_counter + 1
            ! prepare the buffer for receive
            ALLOCATE (buffer_rec(rec_counter)%msg(size_rec_buffer))
            buffer_rec(rec_counter)%msg = 0.0_dp
            buffer_rec(rec_counter)%proc = proc_receive
            ! post the message to be received
            CALL para_env%irecv(buffer_rec(rec_counter)%msg, proc_receive, &
                                buffer_rec(rec_counter)%msg_req)
         END IF
      END DO

      ! send messages
      ALLOCATE (req_send(number_of_send))
      send_counter = 0
      DO proc_shift = 1, para_env%num_pe - 1
         proc_send = MODULO(para_env%mepos + proc_shift, para_env%num_pe)
         IF (map_send_size(proc_send) > 0) THEN
            send_counter = send_counter + 1
            CALL para_env%isend(buffer_send(send_counter)%msg, proc_send, &
                                buffer_send(send_counter)%msg_req)
            req_send(send_counter) = buffer_send(send_counter)%msg_req
         END IF
      END DO
      CALL timestop(handle2)

      ! 6) fill the fm_mat matrix with the messages received from the
      !    other processes
      CALL timeset(routineN//"_fill", handle2)
      CALL cp_fm_get_info(matrix=fm_mat, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local)
      ALLOCATE (sizes(2, 0:para_env%num_pe - 1))
      CALL para_env%allgather([nrow_local, ncol_local], sizes)
      iiB = MAXVAL(sizes(1, :))
      ALLOCATE (index_row_rec(iiB))
      index_row_rec = 0
      rec_counter = 0
      DO proc_shift = 1, para_env%num_pe - 1
         proc_receive = MODULO(para_env%mepos - proc_shift, para_env%num_pe)
         size_rec_buffer = map_rec_size(proc_receive)

         IF (map_rec_size(proc_receive) > 0) THEN
            rec_counter = rec_counter + 1

            rec_col_size = sizes(2, proc_receive)
            rec_row_size = sizes(1, proc_receive)

            ! precompute the indeces of the rows
            num_rec_rows = 0
            DO iiB = 1, rec_row_size
               i_global = fm_mat%matrix_struct%l2g_row(iiB, mepos_2_grid(1, proc_receive))
               IF (i_global >= my_start_row .AND. i_global <= my_end_row) THEN
                  num_rec_rows = num_rec_rows + 1
                  index_row_rec(num_rec_rows) = i_global
               END IF
            END DO

            CALL buffer_rec(rec_counter)%msg_req%wait()

            iii = 0
            DO jjB = 1, rec_col_size
               j_global = fm_mat%matrix_struct%l2g_col(jjB, mepos_2_grid(2, proc_receive))
               IF (j_global >= my_start_col .AND. j_global <= my_end_col) THEN
                  DO iiB = 1, num_rec_rows
                     i_global = index_row_rec(iiB)
                     iii = iii + 1
                     mat2D(i_global - my_start_row + 1, j_global - my_start_col + 1) = buffer_rec(rec_counter)%msg(iii)
                  END DO
               END IF
            END DO

            DEALLOCATE (buffer_rec(rec_counter)%msg)
         END IF
      END DO
      DEALLOCATE (buffer_rec)
      DEALLOCATE (index_row_rec)
      CALL cp_fm_release(fm_mat)
      CALL timestop(handle2)

      ! 7) Finally wait for all messeges to be sent
      CALL timeset(routineN//"_waitall", handle2)
      CALL mp_waitall(req_send(:))
      DO send_counter = 1, number_of_send
         DEALLOCATE (buffer_send(send_counter)%msg)
      END DO
      DEALLOCATE (buffer_send)
      CALL timestop(handle2)

      CALL timestop(handle)

   END SUBROUTINE fm2array

! **************************************************************************************************
!> \brief redistribute 2D representation of 3d tensor to a set of dbcsr matrices
!> \param Gamma_2D ...
!> \param homo ...
!> \param virtual ...
!> \param dimen_ia ...
!> \param para_env_sub ...
!> \param my_ia_start ...
!> \param my_ia_end ...
!> \param my_group_L_size ...
!> \param gd_ia ...
!> \param dbcsr_Gamma ...
!> \param mo_coeff_o ...
! **************************************************************************************************
   SUBROUTINE create_dbcsr_gamma(Gamma_2D, homo, virtual, dimen_ia, para_env_sub, &
                                 my_ia_start, my_ia_end, my_group_L_size, &
                                 gd_ia, dbcsr_Gamma, mo_coeff_o)
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: Gamma_2D
      INTEGER, INTENT(IN)                                :: homo, virtual, dimen_ia
      TYPE(mp_para_env_type), INTENT(IN), POINTER        :: para_env_sub
      INTEGER, INTENT(IN)                                :: my_ia_start, my_ia_end, my_group_L_size
      TYPE(group_dist_d1_type), INTENT(IN)               :: gd_ia
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(INOUT)    :: dbcsr_Gamma
      TYPE(dbcsr_type), INTENT(INOUT)                    :: mo_coeff_o

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_dbcsr_gamma'

      INTEGER :: dummy_proc, handle, i_global, i_local, iaia, iiB, iii, itmp(2), j_global, &
         j_local, jjB, jjj, kkB, mypcol, myprow, ncol_local, npcol, nprow, nrow_local, &
         number_of_rec, number_of_send, proc_receive, proc_send, proc_shift, rec_counter, &
         rec_iaia_end, rec_iaia_size, rec_iaia_start, rec_pcol, rec_prow, ref_send_pcol, &
         ref_send_prow, send_counter, send_pcol, send_prow, size_rec_buffer, size_send_buffer
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: iii_vet, map_rec_size, map_send_size
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: grid_2_mepos, grid_ref_2_send_pos, &
                                                            indeces_map_my, mepos_2_grid
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      REAL(KIND=dp)                                      :: part_ia
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: fm_ia
      TYPE(index_map), ALLOCATABLE, DIMENSION(:)         :: indeces_rec
      TYPE(integ_mat_buffer_type), ALLOCATABLE, &
         DIMENSION(:)                                    :: buffer_rec, buffer_send
      TYPE(mp_request_type), ALLOCATABLE, DIMENSION(:)   :: req_send

      CALL timeset(routineN, handle)

      ! create sub blacs env
      NULLIFY (blacs_env)
      CALL cp_blacs_env_create(blacs_env=blacs_env, para_env=para_env_sub)

      ! create the fm_ia buffer matrix
      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, context=blacs_env, nrow_global=homo, &
                               ncol_global=virtual, para_env=para_env_sub)
      CALL cp_fm_create(fm_ia, fm_struct, name="fm_ia")
      ! release structure
      CALL cp_fm_struct_release(fm_struct)
      ! release blacs_env
      CALL cp_blacs_env_release(blacs_env)

      ! get array information
      CALL cp_fm_get_info(matrix=fm_ia, &
                          nrow_local=nrow_local, &
                          ncol_local=ncol_local, &
                          row_indices=row_indices, &
                          col_indices=col_indices)
      myprow = fm_ia%matrix_struct%context%mepos(1)
      mypcol = fm_ia%matrix_struct%context%mepos(2)
      nprow = fm_ia%matrix_struct%context%num_pe(1)
      npcol = fm_ia%matrix_struct%context%num_pe(2)

      ! 0) create array containing the processes position and supporting infos
      ALLOCATE (grid_2_mepos(0:nprow - 1, 0:npcol - 1))
      grid_2_mepos = 0
      ALLOCATE (mepos_2_grid(2, 0:para_env_sub%num_pe - 1))
      ! fill the info array
      grid_2_mepos(myprow, mypcol) = para_env_sub%mepos
      ! sum infos
      CALL para_env_sub%sum(grid_2_mepos)
      CALL para_env_sub%allgather([myprow, mypcol], mepos_2_grid)

      ! loop over local index range and define the sending map
      ALLOCATE (map_send_size(0:para_env_sub%num_pe - 1))
      map_send_size = 0
      dummy_proc = 0
      DO iaia = my_ia_start, my_ia_end
         i_global = (iaia - 1)/virtual + 1
         j_global = MOD(iaia - 1, virtual) + 1
         send_prow = fm_ia%matrix_struct%g2p_row(i_global)
         send_pcol = fm_ia%matrix_struct%g2p_col(j_global)
         proc_send = grid_2_mepos(send_prow, send_pcol)
         map_send_size(proc_send) = map_send_size(proc_send) + 1
      END DO

      ! loop over local data of fm_ia and define the receiving map
      ALLOCATE (map_rec_size(0:para_env_sub%num_pe - 1))
      map_rec_size = 0
      part_ia = REAL(dimen_ia, KIND=dp)/REAL(para_env_sub%num_pe, KIND=dp)

      DO iiB = 1, nrow_local
         i_global = row_indices(iiB)
         DO jjB = 1, ncol_local
            j_global = col_indices(jjB)
            iaia = (i_global - 1)*virtual + j_global
            proc_receive = INT(REAL(iaia - 1, KIND=dp)/part_ia)
            proc_receive = MAX(0, proc_receive)
            proc_receive = MIN(proc_receive, para_env_sub%num_pe - 1)
            DO
               itmp = get_limit(dimen_ia, para_env_sub%num_pe, proc_receive)
               IF (iaia >= itmp(1) .AND. iaia <= itmp(2)) EXIT
               IF (iaia < itmp(1)) proc_receive = proc_receive - 1
               IF (iaia > itmp(2)) proc_receive = proc_receive + 1
            END DO
            map_rec_size(proc_receive) = map_rec_size(proc_receive) + 1
         END DO
      END DO

      ! allocate the buffer for sending data
      number_of_send = 0
      DO proc_shift = 1, para_env_sub%num_pe - 1
         proc_send = MODULO(para_env_sub%mepos + proc_shift, para_env_sub%num_pe)
         IF (map_send_size(proc_send) > 0) THEN
            number_of_send = number_of_send + 1
         END IF
      END DO
      ! allocate the structure that will hold the messages to be sent
      ALLOCATE (buffer_send(number_of_send))
      ! and the map from the grid of processess to the message position
      ALLOCATE (grid_ref_2_send_pos(0:nprow - 1, 0:npcol - 1))
      grid_ref_2_send_pos = 0
      ! finally allocate each message
      send_counter = 0
      DO proc_shift = 1, para_env_sub%num_pe - 1
         proc_send = MODULO(para_env_sub%mepos + proc_shift, para_env_sub%num_pe)
         size_send_buffer = map_send_size(proc_send)
         IF (map_send_size(proc_send) > 0) THEN
            send_counter = send_counter + 1
            ! allocate the sending buffer (msg)
            ALLOCATE (buffer_send(send_counter)%msg(size_send_buffer))
            buffer_send(send_counter)%proc = proc_send
            ! get the pointer to prow, pcol of the process that has
            ! to receive this message
            ref_send_prow = mepos_2_grid(1, proc_send)
            ref_send_pcol = mepos_2_grid(2, proc_send)
            ! save the rank of the process that has to receive this message
            grid_ref_2_send_pos(ref_send_prow, ref_send_pcol) = send_counter
         END IF
      END DO

      ! allocate the buffer for receiving data
      number_of_rec = 0
      DO proc_shift = 1, para_env_sub%num_pe - 1
         proc_receive = MODULO(para_env_sub%mepos - proc_shift, para_env_sub%num_pe)
         IF (map_rec_size(proc_receive) > 0) THEN
            number_of_rec = number_of_rec + 1
         END IF
      END DO
      ! allocate the structure that will hold the messages to be received
      ! and relative indeces
      ALLOCATE (buffer_rec(number_of_rec))
      ALLOCATE (indeces_rec(number_of_rec))
      ! finally allocate each message and fill the array of indeces
      rec_counter = 0
      DO proc_shift = 1, para_env_sub%num_pe - 1
         proc_receive = MODULO(para_env_sub%mepos - proc_shift, para_env_sub%num_pe)
         size_rec_buffer = map_rec_size(proc_receive)
         IF (map_rec_size(proc_receive) > 0) THEN
            rec_counter = rec_counter + 1
            ! prepare the buffer for receive
            ALLOCATE (buffer_rec(rec_counter)%msg(size_rec_buffer))
            buffer_rec(rec_counter)%proc = proc_receive
            ! create the indeces array
            ALLOCATE (indeces_rec(rec_counter)%map(2, size_rec_buffer))
            indeces_rec(rec_counter)%map = 0
            CALL get_group_dist(gd_ia, proc_receive, rec_iaia_start, rec_iaia_end, rec_iaia_size)
            iii = 0
            DO iaia = rec_iaia_start, rec_iaia_end
               i_global = (iaia - 1)/virtual + 1
               j_global = MOD(iaia - 1, virtual) + 1
               rec_prow = fm_ia%matrix_struct%g2p_row(i_global)
               rec_pcol = fm_ia%matrix_struct%g2p_col(j_global)
               IF (grid_2_mepos(rec_prow, rec_pcol) /= para_env_sub%mepos) CYCLE
               iii = iii + 1
               i_local = fm_ia%matrix_struct%g2l_row(i_global)
               j_local = fm_ia%matrix_struct%g2l_col(j_global)
               indeces_rec(rec_counter)%map(1, iii) = i_local
               indeces_rec(rec_counter)%map(2, iii) = j_local
            END DO
         END IF
      END DO
      ! and create the index map for my local data
      IF (map_rec_size(para_env_sub%mepos) > 0) THEN
         size_rec_buffer = map_rec_size(para_env_sub%mepos)
         ALLOCATE (indeces_map_my(2, size_rec_buffer))
         indeces_map_my = 0
         iii = 0
         DO iaia = my_ia_start, my_ia_end
            i_global = (iaia - 1)/virtual + 1
            j_global = MOD(iaia - 1, virtual) + 1
            rec_prow = fm_ia%matrix_struct%g2p_row(i_global)
            rec_pcol = fm_ia%matrix_struct%g2p_col(j_global)
            IF (grid_2_mepos(rec_prow, rec_pcol) /= para_env_sub%mepos) CYCLE
            iii = iii + 1
            i_local = fm_ia%matrix_struct%g2l_row(i_global)
            j_local = fm_ia%matrix_struct%g2l_col(j_global)
            indeces_map_my(1, iii) = i_local
            indeces_map_my(2, iii) = j_local
         END DO
      END IF

      ! auxiliary vector of indeces for the send buffer
      ALLOCATE (iii_vet(number_of_send))
      ! vector for the send requests
      ALLOCATE (req_send(number_of_send))
      ! loop over auxiliary basis function and redistribute into a fm
      ! and then compy the fm into a dbcsr matrix
      DO kkB = 1, my_group_L_size
         ! zero the matries of the buffers and post the messages to be received
         CALL cp_fm_set_all(matrix=fm_ia, alpha=0.0_dp)
         rec_counter = 0
         DO proc_shift = 1, para_env_sub%num_pe - 1
            proc_receive = MODULO(para_env_sub%mepos - proc_shift, para_env_sub%num_pe)
            IF (map_rec_size(proc_receive) > 0) THEN
               rec_counter = rec_counter + 1
               buffer_rec(rec_counter)%msg = 0.0_dp
               CALL para_env_sub%irecv(buffer_rec(rec_counter)%msg, proc_receive, &
                                       buffer_rec(rec_counter)%msg_req)
            END IF
         END DO
         ! fill the sending buffer and send the messages
         DO send_counter = 1, number_of_send
            buffer_send(send_counter)%msg = 0.0_dp
         END DO
         iii_vet = 0
         jjj = 0
         DO iaia = my_ia_start, my_ia_end
            i_global = (iaia - 1)/virtual + 1
            j_global = MOD(iaia - 1, virtual) + 1
            send_prow = fm_ia%matrix_struct%g2p_row(i_global)
            send_pcol = fm_ia%matrix_struct%g2p_col(j_global)
            proc_send = grid_2_mepos(send_prow, send_pcol)
            ! we don't need to send to ourselves
            IF (grid_2_mepos(send_prow, send_pcol) == para_env_sub%mepos) THEN
               ! filling fm_ia with local data
               jjj = jjj + 1
               i_local = indeces_map_my(1, jjj)
               j_local = indeces_map_my(2, jjj)
               fm_ia%local_data(i_local, j_local) = Gamma_2D(iaia - my_ia_start + 1, kkB)
            ELSE
               send_counter = grid_ref_2_send_pos(send_prow, send_pcol)
               iii_vet(send_counter) = iii_vet(send_counter) + 1
               iii = iii_vet(send_counter)
               buffer_send(send_counter)%msg(iii) = Gamma_2D(iaia - my_ia_start + 1, kkB)
            END IF
         END DO
         req_send = mp_request_null
         send_counter = 0
         DO proc_shift = 1, para_env_sub%num_pe - 1
            proc_send = MODULO(para_env_sub%mepos + proc_shift, para_env_sub%num_pe)
            IF (map_send_size(proc_send) > 0) THEN
               send_counter = send_counter + 1
               CALL para_env_sub%isend(buffer_send(send_counter)%msg, proc_send, &
                                       buffer_send(send_counter)%msg_req)
               req_send(send_counter) = buffer_send(send_counter)%msg_req
            END IF
         END DO

         ! receive the massages and fill the fm_ia
         rec_counter = 0
         DO proc_shift = 1, para_env_sub%num_pe - 1
            proc_receive = MODULO(para_env_sub%mepos - proc_shift, para_env_sub%num_pe)
            size_rec_buffer = map_rec_size(proc_receive)
            IF (map_rec_size(proc_receive) > 0) THEN
               rec_counter = rec_counter + 1
               ! wait for the message
               CALL buffer_rec(rec_counter)%msg_req%wait()
               DO iii = 1, size_rec_buffer
                  i_local = indeces_rec(rec_counter)%map(1, iii)
                  j_local = indeces_rec(rec_counter)%map(2, iii)
                  fm_ia%local_data(i_local, j_local) = buffer_rec(rec_counter)%msg(iii)
               END DO
            END IF
         END DO

         ! wait all
         CALL mp_waitall(req_send(:))

         ! now create the DBCSR matrix and copy fm_ia into it
         ALLOCATE (dbcsr_Gamma(kkB)%matrix)
         CALL cp_dbcsr_m_by_n_from_template(dbcsr_Gamma(kkB)%matrix, &
                                            template=mo_coeff_o, &
                                            m=homo, n=virtual, sym=dbcsr_type_no_symmetry)
         CALL copy_fm_to_dbcsr(fm_ia, dbcsr_Gamma(kkB)%matrix, keep_sparsity=.FALSE.)

      END DO

      ! deallocate stuff
      DEALLOCATE (Gamma_2D)
      DEALLOCATE (iii_vet)
      DEALLOCATE (req_send)
      IF (map_rec_size(para_env_sub%mepos) > 0) THEN
         DEALLOCATE (indeces_map_my)
      END IF
      DO rec_counter = 1, number_of_rec
         DEALLOCATE (indeces_rec(rec_counter)%map)
         DEALLOCATE (buffer_rec(rec_counter)%msg)
      END DO
      DEALLOCATE (indeces_rec)
      DEALLOCATE (buffer_rec)
      DO send_counter = 1, number_of_send
         DEALLOCATE (buffer_send(send_counter)%msg)
      END DO
      DEALLOCATE (buffer_send)
      DEALLOCATE (map_send_size)
      DEALLOCATE (map_rec_size)
      DEALLOCATE (grid_2_mepos)
      DEALLOCATE (mepos_2_grid)

      ! release buffer matrix
      CALL cp_fm_release(fm_ia)

      CALL timestop(handle)

   END SUBROUTINE create_dbcsr_gamma

! **************************************************************************************************
!> \brief prepare array for redistribution
!> \param para_env ...
!> \param para_env_sub ...
!> \param ngroup ...
!> \param group_grid_2_mepos ...
!> \param mepos_2_grid_group ...
!> \param pos_info ...
! **************************************************************************************************
   SUBROUTINE prepare_redistribution(para_env, para_env_sub, ngroup, &
                                     group_grid_2_mepos, mepos_2_grid_group, &
                                     pos_info)
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env, para_env_sub
      INTEGER, INTENT(IN)                                :: ngroup
      INTEGER, ALLOCATABLE, DIMENSION(:, :), INTENT(OUT) :: group_grid_2_mepos, mepos_2_grid_group
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT), &
         OPTIONAL                                        :: pos_info

      INTEGER                                            :: i, pos_group, pos_sub
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: my_pos_info

      ALLOCATE (my_pos_info(0:para_env%num_pe - 1))
      CALL para_env%allgather(para_env_sub%mepos, my_pos_info)

      ALLOCATE (group_grid_2_mepos(0:para_env_sub%num_pe - 1, 0:ngroup - 1))
      group_grid_2_mepos = 0
      ALLOCATE (mepos_2_grid_group(2, 0:para_env%num_pe - 1))
      mepos_2_grid_group = 0

      DO i = 0, para_env%num_pe - 1
         ! calculate position of the group
         pos_group = i/para_env_sub%num_pe
         ! calculate position in the subgroup
         pos_sub = my_pos_info(i)
         ! fill the map from the grid of groups to process
         group_grid_2_mepos(pos_sub, pos_group) = i
         ! and the opposite, from the global pos to the grid pos
         mepos_2_grid_group(1, i) = pos_sub
         mepos_2_grid_group(2, i) = pos_group
      END DO

      IF (PRESENT(pos_info)) CALL move_alloc(my_pos_info, pos_info)

   END SUBROUTINE prepare_redistribution

END MODULE mp2_ri_grad_util
