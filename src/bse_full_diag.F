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
MODULE bse_full_diag
    !MG check if all are necessary; copied from bse.F

    USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                               cp_blacs_env_release,&
                                               cp_blacs_env_type
    USE cp_fm_basic_linalg,              ONLY: cp_fm_upper_to_full
    USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                               cp_fm_cholesky_invert
    USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                               cp_fm_struct_get, &
                                               cp_fm_struct_release,&
                                               cp_fm_struct_type
    USE cp_fm_types,                     ONLY: cp_fm_create,&
                                               cp_fm_get_info,&
                                               cp_fm_release,&
                                               cp_fm_set_all,&
                                               cp_fm_to_fm,&
                                               cp_fm_type, &
                                               cp_fm_write_formatted !MG
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
 
    CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'bse_full_diag'
 
    PUBLIC :: create_and_diagonalize_A
 
 CONTAINS

    SUBROUTINE create_and_diagonalize_A(fm_mat_S, fm_mat_S_bar_ij_bse, fm_mat_S_ab_bse,&
                                        Eigenval, &
                                        homo,virtual, dimen_RI, bse_spin_config, unit_nr, &
                                        para_env, para_env_RPA, blacs_env_BSE)

        TYPE(cp_fm_type), INTENT(IN)                        :: fm_mat_S, fm_mat_S_bar_ij_bse, &
                                                                fm_mat_S_ab_bse
        REAL(KIND=dp), DIMENSION(:)                         :: Eigenval
        INTEGER,INTENT(IN)                                  :: homo, virtual, dimen_RI, bse_spin_config
        !MG for debugging  
        INTEGER                                             :: unit_nr
        TYPE(mp_para_env_type), INTENT(INOUT)               :: para_env
        TYPE(mp_para_env_type), INTENT(IN), POINTER         :: para_env_RPA
        TYPE(cp_blacs_env_type), POINTER                    :: blacs_env_BSE

        CHARACTER(LEN=*), PARAMETER                         :: routineN = 'create_and_diagonalize_A'
        INTEGER                                             :: handle, alpha, ia, jb, ii, jj, & 
                                                               nrow_local_A, ncol_local_A, &
                                                               nrow_local_W, ncol_local_W, &
                                                               i_occ_row, a_virt_row, &
                                                               j_occ_col, b_virt_col
        TYPE(cp_blacs_env_type), POINTER                    :: blacs_env, blacs_env_RPA
        TYPE(cp_fm_struct_type), POINTER                    :: fm_struct_A,fm_struct_W,fm_struct_V
        TYPE(cp_fm_type)                                    :: fm_mat_A, fm_mat_W, fm_mat_V
        INTEGER, DIMENSION(:), POINTER                      :: col_indices_A, row_indices_A, &
                                                               col_indices_W, row_indices_W
        REAL(KIND=dp)                                       :: eigen_diff

        CALL timeset(routineN, handle)

        !Determines factor of exchange term, depending on requested spin configuration (cf. input_constants.F)
        IF (bse_spin_config == 0) THEN
            !Case of singlet
            alpha = 2
        ELSE 
            !Case of triplet
            alpha = 0
        END IF


        ! create the blacs env for BSE matrices (NOT fm_mat_S related parallel_gemms!)
        NULLIFY (blacs_env)
        CALL cp_blacs_env_create(blacs_env=blacs_env, para_env=para_env)

        !We have to use the same blacs_env for A as for the matrices fm_mat_S from RPA
        CALL cp_fm_struct_create(fm_struct_A, context=blacs_env_BSE, nrow_global=homo*virtual, &
                                ncol_global=homo*virtual, para_env=para_env_RPA)
        CALL cp_fm_create(fm_mat_A,fm_struct_A,name="fm_mat_A_iajb_BSE")
        CALL cp_fm_set_all(fm_mat_A, 0.0_dp)

        CALL cp_fm_struct_create(fm_struct_v, context=blacs_env_BSE, nrow_global=homo*virtual, &
                                ncol_global=homo*virtual, para_env=para_env_RPA)
        CALL cp_fm_create(fm_mat_V,fm_struct_v,name="fm_mat_V_iajb_BSE")
        CALL cp_fm_set_all(fm_mat_V, 0.0_dp)

        CALL cp_fm_struct_create(fm_struct_W, context=blacs_env, nrow_global=homo**2, &
                                ncol_global=virtual**2, para_env=para_env)
        CALL cp_fm_create(fm_mat_W,fm_struct_W,name="fm_mat_W_ijab_BSE")
        CALL cp_fm_set_all(fm_mat_W, 0.0_dp)

        CALL cp_fm_get_info(matrix=fm_mat_A, &
                            nrow_local=nrow_local_A, &
                            ncol_local=ncol_local_A, &
                            row_indices=row_indices_A, &
                            col_indices=col_indices_A)
        
        CALL cp_fm_get_info(matrix=fm_mat_W, &
                            nrow_local=nrow_local_W, &
                            ncol_local=ncol_local_W, &
                            row_indices=row_indices_W, &
                            col_indices=col_indices_W)

        
        !print *,"MG: on proc",para_env%mepos,"indices row",row_indices_A
        !print *,"MG: on proc",para_env%mepos,"indices col",col_indices_A
        

        !Create A matrix from GW Energies, v_iajb and W_ijab (different blacs_env!)
        !v_iajb
        CALL parallel_gemm(transa="T", transb="N", m=homo*virtual, n=homo*virtual, k=dimen_RI, alpha=1.0_dp, &
                         matrix_a=fm_mat_S, matrix_b=fm_mat_S, beta=0.0_dp, &
                         matrix_c=fm_mat_V)
                         
        CALL cp_fm_write_formatted(fm_mat_V,unit_nr,"MG: print full V in fulldiag")

        !W_ijab
        CALL parallel_gemm(transa="T", transb="N", m=homo**2, n=virtual**2, k=dimen_RI, alpha=1.0_dp, &
                         matrix_a=fm_mat_S_bar_ij_bse, matrix_b=fm_mat_S_ab_bse, beta=0.0_dp, &
                         matrix_c=fm_mat_W)
                         
        CALL cp_fm_write_formatted(fm_mat_W,unit_nr,"MG: print full W in fulldiag")

        !A_iajb
        DO ii=1,nrow_local_A
            DO jj=1,ncol_local_A
                !Take indices for (un-)occ. orbitals from combined index ia=row_indices_A(ii)
                i_occ_row = (row_indices_A(ii)-1)/virtual+1
                a_virt_row = MOD(row_indices_A(ii)-1,virtual)+1
                j_occ_col = (col_indices_A(jj)-1)/virtual+1
                b_virt_col = MOD(col_indices_A(jj)-1,virtual)+1
                !Get contributions from fm_mat_V and fm_mat_W
                fm_mat_A%local_data(ii,jj) = fm_mat_A%local_data(ii,jj) + alpha*fm_mat_V%local_data(ii,jj)
                !Take care of different indices in A and W
                fm_mat_A%local_data(ii,jj) = fm_mat_A%local_data(ii,jj) &
                                    - fm_mat_W%local_data(i_occ_row*j_occ_col,a_virt_row*b_virt_col)
                

                !fm_mat_A%local_data(ii,jj) = fm_mat_A%local_data(ii,jj) + fm_mat_W

                !Add energy differences on diagonal
                IF (row_indices_A(ii) == col_indices_A(jj)) THEN 
                    !print *, "MG: ia, i_occ_row, a_virt_row", row_indices_A(ii), i_occ_row, a_virt_row
                    eigen_diff = Eigenval(a_virt_row + homo) - Eigenval(i_occ_row)
                    fm_mat_A%local_data(ii,jj) = fm_mat_A%local_data(ii,jj) + eigen_diff
                END IF

            END DO
        END DO


        CALL cp_fm_write_formatted(fm_mat_A,unit_nr,"MG: print full A in fulldiag")



        CALL cp_fm_struct_release(fm_struct_A)
        CALL cp_fm_release(fm_mat_A)
        CALL cp_fm_struct_release(fm_struct_V)
        CALL cp_fm_release(fm_mat_V)
        CALL cp_fm_struct_release(fm_struct_W)
        CALL cp_fm_release(fm_mat_W)
        CALL timestop(handle)

    END SUBROUTINE

END MODULE bse_full_diag
