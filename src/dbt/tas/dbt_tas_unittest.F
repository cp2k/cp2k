!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2022 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Unit testing for tall-and-skinny matrices
!> \author Patrick Seewald
! **************************************************************************************************
PROGRAM dbt_tas_unittest
   USE dbcsr_api,                       ONLY: dbcsr_finalize_lib,&
                                              dbcsr_init_lib
   USE dbm_api,                         ONLY: dbm_get_name,&
                                              dbm_library_finalize,&
                                              dbm_library_init
   USE dbt_tas_base,                    ONLY: dbt_tas_create,&
                                              dbt_tas_destroy,&
                                              dbt_tas_info,&
                                              dbt_tas_nblkcols_total,&
                                              dbt_tas_nblkrows_total
   USE dbt_tas_io,                      ONLY: dbt_tas_write_split_info
   USE dbt_tas_test,                    ONLY: dbt_tas_random_bsizes,&
                                              dbt_tas_reset_randmat_seed,&
                                              dbt_tas_setup_test_matrix,&
                                              dbt_tas_test_mm
   USE dbt_tas_types,                   ONLY: dbt_tas_type
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE machine,                         ONLY: default_output_unit
   USE message_passing,                 ONLY: mp_comm_free,&
                                              mp_environ,&
                                              mp_world_finalize,&
                                              mp_world_init
   USE offload_api,                     ONLY: offload_get_device_count,&
                                              offload_set_chosen_device
#include "../../base/base_uses.f90"

   IMPLICIT NONE

   INTEGER(KIND=int_8), PARAMETER :: m = 100, k = 20, n = 10
   TYPE(dbt_tas_type)             :: A, B, C, At, Bt, Ct, A_out, B_out, C_out, At_out, Bt_out, Ct_out
   INTEGER, DIMENSION(m)          :: bsize_m
   INTEGER, DIMENSION(n)          :: bsize_n
   INTEGER, DIMENSION(k)          :: bsize_k
   REAL(KIND=dp), PARAMETER   :: sparsity = 0.1
   INTEGER                        :: mp_comm, numnodes, mynode, io_unit
   INTEGER                        :: mp_comm_A, mp_comm_At, mp_comm_B, mp_comm_Bt, mp_comm_C, mp_comm_Ct
   REAL(KIND=dp), PARAMETER   :: filter_eps = 1.0E-08

   CALL mp_world_init(mp_comm)

   CALL mp_environ(numnodes, mynode, mp_comm)

   ! Select active offload device when available.
   IF (offload_get_device_count() > 0) THEN
      CALL offload_set_chosen_device(MOD(mynode, offload_get_device_count()))
   END IF

   io_unit = -1
   IF (mynode .EQ. 0) io_unit = default_output_unit

   CALL dbcsr_init_lib(mp_comm, io_unit) ! Needed for DBM_VALIDATE_AGAINST_DBCSR.
   CALL dbm_library_init()

   CALL dbt_tas_reset_randmat_seed()

   CALL dbt_tas_random_bsizes([13, 8, 5, 25, 12], 2, bsize_m)
   CALL dbt_tas_random_bsizes([3, 78, 33, 12, 3, 15], 1, bsize_n)
   CALL dbt_tas_random_bsizes([9, 64, 23, 2], 3, bsize_k)

   CALL mp_environ(numnodes, mynode, mp_comm)

   CALL dbt_tas_setup_test_matrix(A, mp_comm_A, mp_comm, m, k, bsize_m, bsize_k, [5, 1], "A", sparsity)
   CALL dbt_tas_setup_test_matrix(At, mp_comm_At, mp_comm, k, m, bsize_k, bsize_m, [3, 8], "A^t", sparsity)
   CALL dbt_tas_setup_test_matrix(B, mp_comm_B, mp_comm, n, m, bsize_n, bsize_m, [3, 2], "B", sparsity)
   CALL dbt_tas_setup_test_matrix(Bt, mp_comm_Bt, mp_comm, m, n, bsize_m, bsize_n, [1, 3], "B^t", sparsity)
   CALL dbt_tas_setup_test_matrix(C, mp_comm_C, mp_comm, k, n, bsize_k, bsize_n, [5, 7], "C", sparsity)
   CALL dbt_tas_setup_test_matrix(Ct, mp_comm_Ct, mp_comm, n, k, bsize_n, bsize_k, [1, 1], "C^t", sparsity)

   CALL dbt_tas_create(A, A_out)
   CALL dbt_tas_create(At, At_out)
   CALL dbt_tas_create(B, B_out)
   CALL dbt_tas_create(Bt, Bt_out)
   CALL dbt_tas_create(C, C_out)
   CALL dbt_tas_create(Ct, Ct_out)

   IF (mynode == 0) WRITE (io_unit, '(A)') "DBM TALL-AND-SKINNY MATRICES"
   IF (mynode == 0) WRITE (io_unit, '(1X, A, 1X, A, I10, 1X, A, 1X, I10)') "Split info for matrix", &
      TRIM(dbm_get_name(A%matrix)), &
      dbt_tas_nblkrows_total(A), 'X', dbt_tas_nblkcols_total(A)
   CALL dbt_tas_write_split_info(dbt_tas_info(A), io_unit, name="A")
   IF (mynode == 0) WRITE (io_unit, '(1X, A, 1X, A, I10, 1X, A, 1X, I10)') "Split info for matrix", &
      TRIM(dbm_get_name(At%matrix)), &
      dbt_tas_nblkrows_total(At), 'X', dbt_tas_nblkcols_total(At)
   CALL dbt_tas_write_split_info(dbt_tas_info(At), io_unit, name="At")
   IF (mynode == 0) WRITE (io_unit, '(1X, A, 1X, A, I10, 1X, A, 1X, I10)') "Split info for matrix", &
      TRIM(dbm_get_name(B%matrix)), &
      dbt_tas_nblkrows_total(B), 'X', dbt_tas_nblkcols_total(B)
   CALL dbt_tas_write_split_info(dbt_tas_info(B), io_unit, name="B")
   IF (mynode == 0) WRITE (io_unit, '(1X, A, 1X, A, I10, 1X, A, 1X, I10)') "Split info for matrix", &
      TRIM(dbm_get_name(Bt%matrix)), &
      dbt_tas_nblkrows_total(Bt), 'X', dbt_tas_nblkcols_total(Bt)
   CALL dbt_tas_write_split_info(dbt_tas_info(Bt), io_unit, name="Bt")
   IF (mynode == 0) WRITE (io_unit, '(1X, A, 1X, A, I10, 1X, A, 1X, I10)') "Split info for matrix", &
      TRIM(dbm_get_name(C%matrix)), &
      dbt_tas_nblkrows_total(C), 'X', dbt_tas_nblkcols_total(C)
   CALL dbt_tas_write_split_info(dbt_tas_info(C), io_unit, name="C")
   IF (mynode == 0) WRITE (io_unit, '(1X, A, 1X, A, I10, 1X, A, 1X, I10)') "Split info for matrix", &
      TRIM(dbm_get_name(Ct%matrix)), &
      dbt_tas_nblkrows_total(Ct), 'X', dbt_tas_nblkcols_total(Ct)
   CALL dbt_tas_write_split_info(dbt_tas_info(Ct), io_unit, name="Ct")

   CALL dbt_tas_test_mm(.FALSE., .FALSE., .FALSE., B, A, Ct_out, unit_nr=io_unit, filter_eps=filter_eps)
   CALL dbt_tas_test_mm(.TRUE., .FALSE., .FALSE., Bt, A, Ct_out, unit_nr=io_unit, filter_eps=filter_eps)
   CALL dbt_tas_test_mm(.FALSE., .TRUE., .FALSE., B, At, Ct_out, unit_nr=io_unit, filter_eps=filter_eps)
   CALL dbt_tas_test_mm(.TRUE., .TRUE., .FALSE., Bt, At, Ct_out, unit_nr=io_unit, filter_eps=filter_eps)
   CALL dbt_tas_test_mm(.FALSE., .FALSE., .TRUE., B, A, C_out, unit_nr=io_unit, filter_eps=filter_eps)
   CALL dbt_tas_test_mm(.TRUE., .FALSE., .TRUE., Bt, A, C_out, unit_nr=io_unit, filter_eps=filter_eps)
   CALL dbt_tas_test_mm(.FALSE., .TRUE., .TRUE., B, At, C_out, unit_nr=io_unit, filter_eps=filter_eps)
   CALL dbt_tas_test_mm(.TRUE., .TRUE., .TRUE., Bt, At, C_out, unit_nr=io_unit, filter_eps=filter_eps)

   CALL dbt_tas_test_mm(.FALSE., .FALSE., .FALSE., A, C, Bt_out, unit_nr=io_unit, filter_eps=filter_eps)
   CALL dbt_tas_test_mm(.TRUE., .FALSE., .FALSE., At, C, Bt_out, unit_nr=io_unit, filter_eps=filter_eps)
   CALL dbt_tas_test_mm(.FALSE., .TRUE., .FALSE., A, Ct, Bt_out, unit_nr=io_unit, filter_eps=filter_eps)
   CALL dbt_tas_test_mm(.TRUE., .TRUE., .FALSE., At, Ct, Bt_out, unit_nr=io_unit, filter_eps=filter_eps)

   CALL dbt_tas_test_mm(.FALSE., .FALSE., .TRUE., A, C, B_out, unit_nr=io_unit, filter_eps=filter_eps)
   CALL dbt_tas_test_mm(.TRUE., .FALSE., .TRUE., At, C, B_out, unit_nr=io_unit, filter_eps=filter_eps)
   CALL dbt_tas_test_mm(.FALSE., .TRUE., .TRUE., A, Ct, B_out, unit_nr=io_unit, filter_eps=filter_eps)
   CALL dbt_tas_test_mm(.TRUE., .TRUE., .TRUE., At, Ct, B_out, unit_nr=io_unit, filter_eps=filter_eps)

   CALL dbt_tas_test_mm(.FALSE., .FALSE., .FALSE., C, B, At_out, unit_nr=io_unit, filter_eps=filter_eps)
   CALL dbt_tas_test_mm(.TRUE., .FALSE., .FALSE., Ct, B, At_out, unit_nr=io_unit, filter_eps=filter_eps)
   CALL dbt_tas_test_mm(.FALSE., .TRUE., .FALSE., C, Bt, At_out, unit_nr=io_unit, filter_eps=filter_eps)
   CALL dbt_tas_test_mm(.TRUE., .TRUE., .FALSE., Ct, Bt, At_out, unit_nr=io_unit, filter_eps=filter_eps)

   CALL dbt_tas_test_mm(.FALSE., .FALSE., .TRUE., C, B, A_out, unit_nr=io_unit, filter_eps=filter_eps)
   CALL dbt_tas_test_mm(.TRUE., .FALSE., .TRUE., Ct, B, A_out, unit_nr=io_unit, filter_eps=filter_eps)
   CALL dbt_tas_test_mm(.FALSE., .TRUE., .TRUE., C, Bt, A_out, unit_nr=io_unit, filter_eps=filter_eps)
   CALL dbt_tas_test_mm(.TRUE., .TRUE., .TRUE., Ct, Bt, A_out, unit_nr=io_unit, filter_eps=filter_eps)

   CALL dbt_tas_destroy(A)
   CALL dbt_tas_destroy(At)
   CALL dbt_tas_destroy(B)
   CALL dbt_tas_destroy(Bt)
   CALL dbt_tas_destroy(C)
   CALL dbt_tas_destroy(Ct)
   CALL dbt_tas_destroy(A_out)
   CALL dbt_tas_destroy(At_out)
   CALL dbt_tas_destroy(B_out)
   CALL dbt_tas_destroy(Bt_out)
   CALL dbt_tas_destroy(C_out)
   CALL dbt_tas_destroy(Ct_out)

   CALL mp_comm_free(mp_comm_A)
   CALL mp_comm_free(mp_comm_At)
   CALL mp_comm_free(mp_comm_B)
   CALL mp_comm_free(mp_comm_Bt)
   CALL mp_comm_free(mp_comm_C)
   CALL mp_comm_free(mp_comm_Ct)

   CALL dbm_library_finalize()
   CALL dbcsr_finalize_lib() ! Needed for DBM_VALIDATE_AGAINST_DBCSR.
   CALL mp_world_finalize()

END PROGRAM
