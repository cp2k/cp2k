!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
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
   USE cp_dbcsr_redistribution,         ONLY: cp_array_to_dbcsr_slices
   USE cp_fm_basic_linalg,              ONLY: cp_fm_geadd
   USE cp_fm_redistribution,            ONLY: cp_array_to_fm,&
                                              cp_fm_to_array
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_release,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE group_dist_types,                ONLY: create_group_dist,&
                                              get_group_dist,&
                                              group_dist_d1_type,&
                                              release_group_dist
   USE kinds,                           ONLY: dp
   USE libint_2c_3c,                    ONLY: compare_potential_types
   USE message_passing,                 ONLY: mp_para_env_type
   USE mp2_types,                       ONLY: mp2_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mp2_ri_grad_util'

   PUBLIC :: complete_gamma

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
!> \param my_group_L_size ...
!> \param my_group_L_start ...
!> \param my_group_L_end ...
!> \param my_B_size ...
!> \param my_B_virtual_start ...
!> \param gd_B_virtual ...
!> \param kspin ...
!> \author Mauro Del Ben
! **************************************************************************************************
   SUBROUTINE complete_gamma(mp2_env, B_ia_Q, dimen_RI, homo, virtual, para_env, para_env_sub, my_group_L_size, &
                             my_group_L_start, my_group_L_end, my_B_size, my_B_virtual_start, gd_B_virtual, kspin)

      TYPE(mp2_type)                                     :: mp2_env
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(INOUT)                                   :: B_ia_Q
      INTEGER, INTENT(IN)                                :: dimen_RI, homo, virtual
      TYPE(mp_para_env_type), INTENT(IN), POINTER        :: para_env, para_env_sub
      INTEGER, INTENT(IN)                                :: my_group_L_size, my_group_L_start, &
                                                            my_group_L_end, my_B_size, &
                                                            my_B_virtual_start
      TYPE(group_dist_d1_type), INTENT(IN)               :: gd_B_virtual
      INTEGER, INTENT(IN)                                :: kspin

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'complete_gamma'

      INTEGER                                            :: dimen_ia, handle, ispin, kkB, my_ia_end, &
                                                            my_ia_size, my_ia_start, my_P_end, &
                                                            my_P_size, my_P_start, nspins
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: BIb_C_2D, Gamma_2D, Gamma_PQ
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_ia, fm_struct_RI
      TYPE(cp_fm_type) :: fm_Gamma, fm_Gamma_PQ, fm_Gamma_PQ_2, fm_Gamma_PQ_temp, &
         fm_Gamma_PQ_temp_2, fm_ia_P, fm_Y, operator_half, PQ_half
      TYPE(group_dist_d1_type)                           :: gd_ia, gd_P

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

      CALL create_group_dist(gd_P, para_env_sub%num_pe, dimen_RI)
      CALL get_group_dist(gd_P, para_env_sub%mepos, my_P_start, my_P_end, my_P_size)

      ! create the blacs env
      NULLIFY (blacs_env)
      CALL cp_blacs_env_create(blacs_env=blacs_env, para_env=para_env)

      NULLIFY (fm_struct_ia)
      CALL cp_fm_struct_create(fm_struct_ia, context=blacs_env, nrow_global=dimen_ia, &
                               ncol_global=dimen_RI, para_env=para_env)

      ! create the fm matrix Gamma
      CALL cp_array_to_fm(Gamma_2D, fm_struct_ia, my_ia_start, my_ia_end, my_group_L_start, my_group_L_end, &
                          fm_Y)
      ! create the fm matrix B_ia_P
      CALL cp_array_to_fm(BIb_C_2D, fm_struct_ia, my_ia_start, my_ia_end, my_group_L_start, my_group_L_end, &
                          fm_ia_P)

      NULLIFY (fm_struct_RI)
      CALL cp_fm_struct_create(fm_struct_RI, context=blacs_env, nrow_global=dimen_RI, &
                               ncol_global=dimen_RI, para_env=para_env)

      ! since we will need (P|Q)^(-1/2) in the future, make a copy
      CALL cp_array_to_fm(mp2_env%ri_grad%PQ_half, fm_struct_RI, my_P_start, my_P_end, my_group_L_start, &
                          my_group_L_end, PQ_half, do_release_mat=.FALSE.)

      IF (.NOT. compare_potential_types(mp2_env%ri_metric, mp2_env%potential_parameter)) THEN
         CALL cp_array_to_fm(mp2_env%ri_grad%operator_half, fm_struct_RI, my_P_start, my_P_end, &
                             my_group_L_start, my_group_L_end, operator_half, do_release_mat=.FALSE.)
      END IF

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
      CALL cp_fm_to_array(Gamma_2D, my_ia_start, my_ia_end, my_group_L_start, my_group_L_end, fm_Gamma)

      CALL cp_fm_to_array(Gamma_PQ, my_P_start, my_P_end, my_group_L_start, my_group_L_end, fm_Gamma_PQ)
      IF (.NOT. ALLOCATED(mp2_env%ri_grad%Gamma_PQ)) THEN
         CALL MOVE_ALLOC(Gamma_PQ, mp2_env%ri_grad%Gamma_PQ)
      ELSE
         mp2_env%ri_grad%Gamma_PQ(:, :) = mp2_env%ri_grad%Gamma_PQ + Gamma_PQ
         DEALLOCATE (Gamma_PQ)
      END IF

      IF (.NOT. compare_potential_types(mp2_env%ri_metric, mp2_env%potential_parameter)) THEN
         CALL cp_fm_to_array(Gamma_PQ, my_P_start, my_P_end, my_group_L_start, my_group_L_end, fm_Gamma_PQ_2)
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
      CALL cp_array_to_dbcsr_slices(Gamma_2D, homo, virtual, para_env_sub, my_ia_start, my_ia_end, &
                                    mp2_env%ri_grad%mo_coeff_o(1)%matrix, &
                                    matrix_ptrs=mp2_env%ri_grad%G_P_ia(:, kspin))

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

END MODULE mp2_ri_grad_util
