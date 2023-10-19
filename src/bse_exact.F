!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for GW + Bethe-Salpeter for computing electronic excitations
!> \par History
!>      04.2017 created [Jan Wilhelm]
! **************************************************************************************************
MODULE bse_exact
    !MG check if all are necessary; copied from bse.F
    USE cp_fm_basic_linalg,              ONLY: cp_fm_upper_to_full
    USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                               cp_fm_cholesky_invert
    USE cp_fm_types,                     ONLY: cp_fm_create,&
                                               cp_fm_get_info,&
                                               cp_fm_release,&
                                               cp_fm_set_all,&
                                               cp_fm_to_fm,&
                                               cp_fm_type
    USE group_dist_types,                ONLY: get_group_dist,&
                                               group_dist_d1_type
    USE kinds,                           ONLY: dp
    USE message_passing,                 ONLY: mp_para_env_type,&
                                               mp_request_type
    USE mp2_types,                       ONLY: integ_mat_buffer_type
    USE parallel_gemm_api,               ONLY: parallel_gemm
    USE physcon,                         ONLY: evolt
    USE rpa_communication,               ONLY: communicate_buffer
#include "./base/base_uses.f90"
 
    IMPLICIT NONE
 
    PRIVATE
 
    CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'bse_exact'
 
    PUBLIC :: create_and_diagonalize_A
 
 CONTAINS

    SUBROUTINE create_and_diagonalize_A(B_bar_ijQ_bse_local, B_abQ_bse_local, B_iaQ_bse_local,&
                                        homo,virtual)

        REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
        INTENT(OUT)                                     :: B_bar_ijQ_bse_local, B_abQ_bse_local, &
                                                           B_iaQ_bse_local
        INTEGER,INTENT(IN)                              :: homo, virtual

        CHARACTER(LEN=*), PARAMETER                         :: routineN = 'do_subspace_iterations'
        INTEGER                                             :: handle
        TYPE(cp_fm_type)                                    :: fm_mat_A
        TYPE(cp_fm_struct_type), POINTER                    :: fm_struct_A

        CALL timeset(routineN, handle)

        CALL cp_fm_struct_create(fm_struct_A, context=blacs_env, nrow_global=homo*virtual, &
                                ncol_global=homo*virtual, para_env=para_env_RPA, &
                                nrow_block=nrow_block_mat, ncol_block=ncol_block_mat,)
        CALL cp_fm_create(fm_mat_A,fm_struct_A)

        !Hier dann parallel_gemm (a la mult_B_with_W_and_fill_local_3c_arrays)
        !dafür braucht man evtl fm_mat als Input (also fm_mat_S statt lokale B Matrizen)
        
        CALL timestop(handle)

    END SUBROUTINE

END MODULE bse_exact
