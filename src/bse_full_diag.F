!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines for the full diagonalization of GW + Bethe-Salpeter for computing
!> \brief electronic excitations
!> \par History
!>      10.2023 created [Maximilian Graml]
! **************************************************************************************************
MODULE bse_full_diag

   USE bse_util,                        ONLY: redist_fm
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_type
   USE cp_cfm_diag,                     ONLY: cp_cfm_heevd
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_release,&
                                              cp_cfm_type,&
                                              cp_fm_to_cfm
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose
   USE cp_fm_diag,                      ONLY: choose_eigv_solver,&
                                              cp_fm_power
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_to_fm_submat,&
                                              cp_fm_type
   USE input_constants,                 ONLY: bse_singlet,&
                                              bse_triplet
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE physcon,                         ONLY: evolt
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'bse_full_diag'

   PUBLIC :: create_A, diagonalize_A, create_B, create_ABBA, diagonalize_ABBA_via_skew, create_skew_of_ABBA, &
             create_hermitian_form_of_ABBA, diagonalize_C, reduce_BSE_matrices

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param fm_mat_S ...
!> \param fm_mat_S_bar_ij_bse ...
!> \param fm_mat_S_ab_bse ...
!> \param fm_mat_A ...
!> \param Eigenval ...
!> \param homo ...
!> \param virtual ...
!> \param dimen_RI ...
!> \param bse_spin_config ...
!> \param unit_nr ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE create_A(fm_mat_S, fm_mat_S_bar_ij_bse, fm_mat_S_ab_bse, &
                       fm_mat_A, Eigenval, &
                       homo, virtual, dimen_RI, bse_spin_config, unit_nr, &
                       para_env)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat_S, fm_mat_S_bar_ij_bse, &
                                                            fm_mat_S_ab_bse
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_mat_A
      REAL(KIND=dp), DIMENSION(:)                        :: Eigenval
      INTEGER, INTENT(IN)                                :: homo, virtual, dimen_RI, bse_spin_config
      INTEGER                                            :: unit_nr
      TYPE(mp_para_env_type), INTENT(INOUT)              :: para_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'create_A'

      INTEGER                                            :: a_virt_row, handle, i_occ_row, &
                                                            i_row_global, ii, j_col_global, jj, &
                                                            ncol_local_A, nrow_local_A
      INTEGER, DIMENSION(4)                              :: reordering
      INTEGER, DIMENSION(:), POINTER                     :: col_indices_A, row_indices_A
      REAL(KIND=dp)                                      :: alpha, eigen_diff
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_A, fm_struct_W
      TYPE(cp_fm_type)                                   :: fm_mat_W

      CALL timeset(routineN, handle)

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T3,A1,T4,A11)') '|', '-Creating A'
      END IF

      !Determines factor of exchange term, depending on requested spin configuration (cf. input_constants.F)
      SELECT CASE (bse_spin_config)
      CASE (bse_singlet)
         alpha = 2.0_dp
      CASE (bse_triplet)
         alpha = 0.0_dp
      END SELECT

      ! create the blacs env for ij matrices (NOT fm_mat_S%matrix_struct related parallel_gemms!)
      NULLIFY (blacs_env)
      CALL cp_blacs_env_create(blacs_env=blacs_env, para_env=para_env)

      !We have to use the same blacs_env for A as for the matrices fm_mat_S from RPA
      !Logic: A_iajb = (eps_a-eps_i) delta_ij delta _ab + alpha * v_iajb - W_ijab
      ! We create v_iajb and W_ijab, then we communicate entries from local W_ijab
      ! to the full matrix v_iajb. By adding these and the energy diffenences: v_iajb -> A_iajb
      ! We use the A matrix already from the start instead of v
      CALL cp_fm_struct_create(fm_struct_a, context=fm_mat_S%matrix_struct%context, nrow_global=homo*virtual, &
                               ncol_global=homo*virtual, para_env=fm_mat_S%matrix_struct%para_env)
      CALL cp_fm_create(fm_mat_A, fm_struct_a, name="fm_A_iajb")
      CALL cp_fm_set_all(fm_mat_A, 0.0_dp)

      CALL cp_fm_struct_create(fm_struct_W, context=blacs_env, nrow_global=homo**2, &
                               ncol_global=virtual**2, para_env=para_env)
      CALL cp_fm_create(fm_mat_W, fm_struct_W, name="fm_W_ijab")
      CALL cp_fm_set_all(fm_mat_W, 0.0_dp)

      !Create A matrix from GW Energies, v_iajb and W_ijab (different blacs_env!)
      !v_iajb, which is directly initialized in A (with a factor of alpha)
      !MG to do: speed up by IF(bse_spin_config) to check for necessity of gemm?
      CALL parallel_gemm(transa="T", transb="N", m=homo*virtual, n=homo*virtual, k=dimen_RI, alpha=alpha, &
                         matrix_a=fm_mat_S, matrix_b=fm_mat_S, beta=0.0_dp, &
                         matrix_c=fm_mat_A)

      !CALL cp_fm_write_formatted(fm_mat_A, unit_nr, "MG: alpha* Full V")
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T3,A1,T8,A16)') '|', 'Allocated A_iajb'
      END IF

      !W_ijab
      CALL parallel_gemm(transa="T", transb="N", m=homo**2, n=virtual**2, k=dimen_RI, alpha=1.0_dp, &
                         matrix_a=fm_mat_S_bar_ij_bse, matrix_b=fm_mat_S_ab_bse, beta=0.0_dp, &
                         matrix_c=fm_mat_W)
      !CALL cp_fm_write_formatted(fm_mat_W, unit_nr, "MG: Full W")
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T3,A1,T8,A16)') '|', 'Allocated W_ijab'
      END IF

      ! We start by moving data from local parts of W_ijab to the full matrix A_iajb using buffers
      CALL cp_fm_get_info(matrix=fm_mat_A, &
                          nrow_local=nrow_local_A, &
                          ncol_local=ncol_local_A, &
                          row_indices=row_indices_A, &
                          col_indices=col_indices_A)
      ! Writing -1.0_dp * W_ijab to A_iajb, i.e. beta = -1.0_dp,
      ! nrow_secidx_in = homo, ncol_secidx_in = virtual !from W_ijab
      ! nrow_secidx_out = virtual, ncol_secidx_out = virtual !from A_iajb
      reordering = (/1, 3, 2, 4/)
      CALL redist_fm(fm_mat_A, fm_mat_W, -1.0_dp, homo, virtual, &
                     virtual, virtual, unit_nr, reordering)
      !full matrix W is not needed anymore, release it to save memory
      CALL cp_fm_struct_release(fm_struct_W)
      CALL cp_fm_release(fm_mat_W)
      CALL cp_fm_struct_release(fm_struct_A)

      !Now add the energy differences to A
      DO ii = 1, nrow_local_A

         i_row_global = row_indices_A(ii)

         DO jj = 1, ncol_local_A

            j_col_global = col_indices_A(jj)

            IF (i_row_global == j_col_global) THEN
               i_occ_row = (i_row_global - 1)/virtual + 1
               a_virt_row = MOD(i_row_global - 1, virtual) + 1
               eigen_diff = Eigenval(a_virt_row + homo) - Eigenval(i_occ_row)
               fm_mat_A%local_data(ii, jj) = fm_mat_A%local_data(ii, jj) + eigen_diff

            END IF
         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE create_A

! **************************************************************************************************
!> \brief ...
!> \param fm_mat_S ...
!> \param fm_mat_S_bar_ia_bse ...
!> \param fm_mat_B ...
!> \param homo ...
!> \param virtual ...
!> \param dimen_RI ...
!> \param bse_spin_config ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE create_B(fm_mat_S, fm_mat_S_bar_ia_bse, fm_mat_B, &
                       homo, virtual, dimen_RI, bse_spin_config, unit_nr)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat_S, fm_mat_S_bar_ia_bse
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_mat_B
      INTEGER, INTENT(IN)                                :: homo, virtual, dimen_RI, bse_spin_config
      INTEGER                                            :: unit_nr

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'create_B'

      INTEGER                                            :: handle
      INTEGER, DIMENSION(4)                              :: reordering
      REAL(KIND=dp)                                      :: alpha
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_v
      TYPE(cp_fm_type)                                   :: fm_mat_W

      CALL timeset(routineN, handle)

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T3,A1,T4,A11)') '|', '-Creating B'
      END IF

      !Determines factor of exchange term, depending on requested spin configuration (cf. input_constants.F)
      SELECT CASE (bse_spin_config)
      CASE (bse_singlet)
         alpha = 2.0_dp
      CASE (bse_triplet)
         alpha = 0.0_dp
      END SELECT

      CALL cp_fm_struct_create(fm_struct_v, context=fm_mat_S%matrix_struct%context, nrow_global=homo*virtual, &
                               ncol_global=homo*virtual, para_env=fm_mat_S%matrix_struct%para_env)
      CALL cp_fm_create(fm_mat_B, fm_struct_v, name="fm_B_iajb")
      CALL cp_fm_set_all(fm_mat_B, 0.0_dp)

      ! CALL cp_fm_struct_create(fm_struct_W, context=blacs_env_BSE, nrow_global=homo*virtual, &
      !                          ncol_global=homo*virtual, para_env=para_env_RPA)
      CALL cp_fm_create(fm_mat_W, fm_struct_v, name="fm_W_ibaj")
      CALL cp_fm_set_all(fm_mat_W, 0.0_dp)

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T3,A1,T8,A16)') '|', 'Allocated B_iajb'
      END IF

      !MG to do: speed up by IF(bse_spin_config) to check for necessity of gemm?
      CALL parallel_gemm(transa="T", transb="N", m=homo*virtual, n=homo*virtual, k=dimen_RI, alpha=alpha, &
                         matrix_a=fm_mat_S, matrix_b=fm_mat_S, beta=0.0_dp, &
                         matrix_c=fm_mat_B)

      !W_ibaj
      CALL parallel_gemm(transa="T", transb="N", m=homo*virtual, n=homo*virtual, k=dimen_RI, alpha=1.0_dp, &
                         matrix_a=fm_mat_S_bar_ia_bse, matrix_b=fm_mat_S, beta=0.0_dp, &
                         matrix_c=fm_mat_W)

      !from W_iajb to W_ibja(formally: W_ibaj, but our internal indexorder is different)
      reordering = (/1, 4, 3, 2/)
      CALL redist_fm(fm_mat_B, fm_mat_W, -1.0_dp, virtual, virtual, &
                     virtual, virtual, unit_nr, reordering)

      CALL cp_fm_release(fm_mat_W)
      CALL timestop(handle)

   END SUBROUTINE create_B

! **************************************************************************************************
!> \brief ...
!> \param fm_mat_A ...
!> \param fm_mat_B ...
!> \param fm_mat_ABBA ...
!> \param para_env ...
!> \param blacs_env_BSE ...
!> \param homo ...
!> \param virtual ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE create_ABBA(fm_mat_A, fm_mat_B, fm_mat_ABBA, para_env, &
                          blacs_env_BSE, homo, virtual, unit_nr)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat_A, fm_mat_B
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_mat_ABBA
      TYPE(mp_para_env_type), INTENT(IN), POINTER        :: para_env
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_BSE
      INTEGER, INTENT(IN)                                :: homo, virtual, unit_nr

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'create_ABBA'

      INTEGER                                            :: handle
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_abba

      CALL timeset(routineN, handle)

      !Create ABBA
      CALL cp_fm_struct_create(fm_struct_abba, context=blacs_env_BSE, nrow_global=2*homo*virtual, &
                               ncol_global=2*homo*virtual, para_env=para_env)
      CALL cp_fm_create(fm_mat_ABBA, fm_struct_abba, name="fm_ABBA_BSE")
      CALL cp_fm_set_all(fm_mat_ABBA, 0.0_dp)

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T3,A1,T8,A14)') '|', 'Allocated ABBA'
      END IF

      !Copy each part
      CALL cp_fm_to_fm_submat(fm_mat_A, fm_mat_ABBA, homo*virtual, homo*virtual, &
                              1, 1, 1, 1)
      CALL cp_fm_to_fm_submat(fm_mat_B, fm_mat_ABBA, homo*virtual, homo*virtual, &
                              1, 1, 1, homo*virtual + 1)
      !Change sign for second row
      fm_mat_A%local_data = -fm_mat_A%local_data
      fm_mat_B%local_data = -fm_mat_B%local_data
      !Now fill the second row
      CALL cp_fm_to_fm_submat(fm_mat_B, fm_mat_ABBA, homo*virtual, homo*virtual, &
                              1, 1, homo*virtual + 1, 1)
      CALL cp_fm_to_fm_submat(fm_mat_A, fm_mat_ABBA, homo*virtual, homo*virtual, &
                              1, 1, homo*virtual + 1, homo*virtual + 1)

      CALL timestop(handle)
   END SUBROUTINE

   ! **************************************************************************************************
!> \brief Prepare ABBA for Diag. following Penke et al., Parallel Computing 2020
!> \brief https://doi.org/10.1016/j.parco.2020.102639
!> \param fm_mat_A ...
!> \param fm_mat_B ...
!> \param fm_mat_C ...
!> \param homo ...
!> \param virtual ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE create_hermitian_form_of_ABBA(fm_mat_A, fm_mat_B, fm_mat_C, &
                                            homo, virtual, unit_nr)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat_A, fm_mat_B
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_mat_C
      INTEGER, INTENT(IN)                                :: homo, virtual, unit_nr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_hermitian_form_of_ABBA'

      INTEGER                                            :: dim_mat, handle, n_dependent
      TYPE(cp_fm_type)                                   :: fm_mat_work_diff, &
                                                            fm_mat_work_diff_eigvec, &
                                                            fm_mat_work_product, fm_mat_work_sum

      CALL timeset(routineN, handle)

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T3,A1,T4,A46)') '|', '-Creating C for hermitian solution of full BSE'
      END IF
      !Create C as A+B
      CALL cp_fm_create(fm_mat_C, fm_mat_A%matrix_struct)
      CALL cp_fm_to_fm(fm_mat_A, fm_mat_C)

      !Create work matrix, which will hold A+b and A-B
      CALL cp_fm_create(fm_mat_work_sum, fm_mat_A%matrix_struct)
      CALL cp_fm_to_fm(fm_mat_A, fm_mat_work_sum)
      CALL cp_fm_create(fm_mat_work_diff, fm_mat_A%matrix_struct)
      CALL cp_fm_to_fm(fm_mat_A, fm_mat_work_diff)

      CALL cp_fm_create(fm_mat_work_diff_eigvec, fm_mat_A%matrix_struct)
      CALL cp_fm_create(fm_mat_work_product, fm_mat_A%matrix_struct)

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T3,A1,T8,A19)') '|', 'Created work arrays'
      END IF
      !Add/Substract B
      CALL cp_fm_scale_and_add(1.0_dp, fm_mat_work_sum, 1.0_dp, fm_mat_B)
      CALL cp_fm_scale_and_add(1.0_dp, fm_mat_work_diff, -1.0_dp, fm_mat_B)

      !Now raise A-B to 1/2, no quenching of eigenvectors, hence threshold=0.0_dp
      CALL cp_fm_power(fm_mat_work_diff, fm_mat_work_diff_eigvec, 0.5_dp, 0.0_dp, n_dependent)
      !(A-B)**-0.5 is written to the eigvec - work matrix
      CALL cp_fm_release(fm_mat_work_diff_eigvec)

      !Create C = (A-B)^0.5 (A+B) (A-B)^0.5
      dim_mat = homo*virtual
      CALL parallel_gemm("N", "N", dim_mat, dim_mat, dim_mat, 1.0_dp, fm_mat_work_diff, fm_mat_work_sum, 0.0_dp, &
                         fm_mat_work_product)
      CALL parallel_gemm("N", "N", dim_mat, dim_mat, dim_mat, 1.0_dp, fm_mat_work_product, fm_mat_work_diff, 0.0_dp, &
                         fm_mat_C)
      CALL cp_fm_release(fm_mat_work_product)
      CALL cp_fm_release(fm_mat_work_sum)
      CALL cp_fm_release(fm_mat_work_diff)
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T3,A1,T8,A36)') '|', 'Filled C=(A-B)^-0.5 (A+B) (A-B)^-0.5'
      END IF

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param fm_mat_C ...
!> \param homo ...
!> \param virtual ...
!> \param bse_spin_config ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE diagonalize_C(fm_mat_C, homo, virtual, bse_spin_config, unit_nr)

      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_mat_C
      INTEGER, INTENT(IN)                                :: homo, virtual, bse_spin_config
      INTEGER                                            :: unit_nr

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'diagonalize_C'

      INTEGER                                            :: diag_info, handle
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Exc_ens
      TYPE(cp_fm_type)                                   :: fm_mat_eigvec

      CALL timeset(routineN, handle)

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T3,A1,T8,A15)') '|', 'Diagonalizing C'
      END IF

      !We have now the full matrix C**2, distributed over all ranks
      !Now: Diagonalize it
      CALL cp_fm_create(fm_mat_eigvec, fm_mat_C%matrix_struct)

      ALLOCATE (Exc_ens(homo*virtual))

      CALL choose_eigv_solver(fm_mat_C, fm_mat_eigvec, Exc_ens, diag_info)
      CPASSERT(diag_info == 0)
      Exc_ens = SQRT(Exc_ens)
      CALL success_message(Exc_ens, bse_spin_config, unit_nr, homo*virtual, .FALSE.)

      CALL cp_fm_release(fm_mat_eigvec)
      DEALLOCATE (Exc_ens)

      CALL timestop(handle)

   END SUBROUTINE diagonalize_C

! **************************************************************************************************
!> \brief Prepare ABBA for Diag. following Penke et al., Parallel Computing 2020
!> \brief https://doi.org/10.1016/j.parco.2020.102639
!> \param fm_mat_A ...
!> \param fm_mat_B ...
!> \param fm_mat_skew ...
!> \param para_env ...
!> \param blacs_env_BSE ...
!> \param homo ...
!> \param virtual ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE create_skew_of_ABBA(fm_mat_A, fm_mat_B, fm_mat_skew, para_env, &
                                  blacs_env_BSE, homo, virtual, unit_nr)

      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat_A, fm_mat_B
      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_mat_skew
      TYPE(mp_para_env_type), INTENT(IN), POINTER        :: para_env
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_BSE
      INTEGER, INTENT(IN)                                :: homo, virtual, unit_nr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'create_skew_of_ABBA'

      INTEGER                                            :: handle, info_cholesky
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_skew
      TYPE(cp_fm_type)                                   :: fm_mat_work

      CALL timeset(routineN, handle)

      !Create ABBA
      CALL cp_fm_struct_create(fm_struct_skew, context=blacs_env_BSE, nrow_global=2*homo*virtual, &
                               ncol_global=2*homo*virtual, para_env=para_env)
      CALL cp_fm_create(fm_mat_skew, fm_struct_skew, name="fm_ABBA_skew")
      CALL cp_fm_set_all(fm_mat_skew, 0.0_dp)

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T3,A1,T8,A14)') '|', 'Allocated skew of ABBA'
      END IF

      !Create work matrix, which will hold A+B and A-B
      CALL cp_fm_create(fm_mat_work, fm_mat_A%matrix_struct)
      CALL cp_fm_to_fm(fm_mat_A, fm_mat_work)

      fm_mat_work%local_data = fm_mat_work%local_data + fm_mat_B%local_data

      !Copy A+B to 11-block
      CALL cp_fm_to_fm_submat(fm_mat_work, fm_mat_skew, homo*virtual, homo*virtual, &
                              1, 1, 1, 1)

      !Now A-B
      CALL cp_fm_to_fm(fm_mat_A, fm_mat_work)

      fm_mat_work%local_data = fm_mat_work%local_data - fm_mat_B%local_data

      !Copy A-B to 22-block
      CALL cp_fm_to_fm_submat(fm_mat_work, fm_mat_skew, homo*virtual, homo*virtual, &
                              1, 1, homo*virtual + 1, homo*virtual + 1)

      !Now we decompose it via Cholesky -> Matrix is replaced by its decomposition!
      CALL cp_fm_cholesky_decompose(fm_mat_skew, 2*homo*virtual, info_cholesky)
      CPASSERT(info_cholesky == 0)
      !.... Paper uses lower triangular?
      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Diagonalize ABBA following Penke et al., Parallel Computing 2020
!> \brief https://doi.org/10.1016/j.parco.2020.102639
!> \param fm_mat_skew ...
!> \param homo ...
!> \param virtual ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE diagonalize_ABBA_via_skew(fm_mat_skew, homo, virtual, unit_nr)

      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_mat_skew
      INTEGER, INTENT(IN)                                :: homo, virtual, unit_nr

      CHARACTER(LEN=*), PARAMETER :: routineN = 'diagonalize_ABBA_via_skew'

      INTEGER                                            :: diag_info, handle
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues, Exc_ens
      TYPE(cp_cfm_type)                                  :: cfm_mat_eigvec, cfm_mat_skew
      TYPE(cp_fm_type)                                   :: fm_mat_eigvec

      CALL timeset(routineN, handle)
      ALLOCATE (Exc_ens(2*homo*virtual))
      CALL cp_fm_create(fm_mat_eigvec, fm_mat_skew%matrix_struct)

      CALL choose_eigv_solver(fm_mat_skew, fm_mat_eigvec, Exc_ens, diag_info)
      CPASSERT(diag_info == 0)
      PRINT *, "MG - Full BSE eigenvalues from fm (no B)", Exc_ens*evolt
      DEALLOCATE (Exc_ens)
      !Use existing cfm routines for benchmark
      CALL cp_cfm_create(cfm_mat_skew, fm_mat_skew%matrix_struct, "cfm_ABBA_skew")
      CALL cp_cfm_create(cfm_mat_eigvec, fm_mat_skew%matrix_struct, "cfm_ABBA_eigvec")

      !Create hermitian matrix by multiplication with imaginary unit
      CALL cp_fm_to_cfm(msourcei=fm_mat_skew, mtarget=cfm_mat_skew)
      ALLOCATE (eigenvalues(2*homo*virtual))
      CALL cp_cfm_heevd(cfm_mat_skew, cfm_mat_eigvec, eigenvalues)
      PRINT *, "MG evolt", evolt
      PRINT *, "MG - Full BSE eigenvalues from cfm (no B)", eigenvalues*evolt
      DEALLOCATE (eigenvalues)
      CALL cp_cfm_release(cfm_mat_eigvec)
      CALL cp_cfm_release(cfm_mat_skew)
      CALL timestop(handle)

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param fm_mat_A ...
!> \param homo ...
!> \param virtual ...
!> \param bse_spin_config ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE diagonalize_A(fm_mat_A, homo, virtual, bse_spin_config, unit_nr)

      TYPE(cp_fm_type), INTENT(INOUT)                    :: fm_mat_A
      INTEGER, INTENT(IN)                                :: homo, virtual, bse_spin_config
      INTEGER                                            :: unit_nr

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'diagonalize_A'

      INTEGER                                            :: diag_info, handle
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Exc_ens
      TYPE(cp_fm_type)                                   :: fm_mat_eigvec

      CALL timeset(routineN, handle)
      !Continue with formatting of subroutine create_A
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T3,A1,T8,A15)') '|', 'Diagonalizing A'
      END IF

      !We have now the full matrix A_iajb, distributed over all ranks
      !Now: Diagonalize it
      CALL cp_fm_create(fm_mat_eigvec, fm_mat_A%matrix_struct)

      ALLOCATE (Exc_ens(homo*virtual))

      CALL choose_eigv_solver(fm_mat_A, fm_mat_eigvec, Exc_ens, diag_info)
      CPASSERT(diag_info == 0)
      CALL success_message(Exc_ens, bse_spin_config, unit_nr, homo*virtual, .TRUE.)

      CALL cp_fm_release(fm_mat_eigvec)
      DEALLOCATE (Exc_ens)

      CALL timestop(handle)

   END SUBROUTINE diagonalize_A

! **************************************************************************************************
!> \brief ...
!> \param Exc_ens ...
!> \param bse_spin_config ...
!> \param unit_nr ...
!> \param num_exc_en ...
!> \param flag_TDA ...
! **************************************************************************************************
   SUBROUTINE success_message(Exc_ens, bse_spin_config, unit_nr, num_exc_en, flag_TDA)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: Exc_ens
      INTEGER                                            :: bse_spin_config, unit_nr, num_exc_en
      LOGICAL, OPTIONAL                                  :: flag_TDA

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'success_message'

      CHARACTER(LEN=10)                                  :: info_approximation, multiplet
      INTEGER                                            :: handle, i
      REAL(KIND=dp)                                      :: alpha

      CALL timeset(routineN, handle)

      !Prepare variables for printing
      IF (bse_spin_config == 0) THEN
         multiplet = "Singlet"
         alpha = 2.0_dp
      ELSE
         multiplet = "Triplet"
         alpha = 0.0_dp
      END IF
      IF (.NOT. PRESENT(flag_TDA)) THEN
         flag_TDA = .FALSE.
      END IF
      IF (flag_TDA) THEN
         info_approximation = "TDA "
      ELSE
         info_approximation = "full"
      END IF

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, *) ' '
         WRITE (unit_nr, '(T3,A)') '******************************************************************************'
         WRITE (unit_nr, '(T3,A)') '**                                                                          **'
         WRITE (unit_nr, '(T3,A30,A4,A44)') '**                        BSE-', info_approximation, &
            ' EXCITONIC ENERGIES                       **'
         WRITE (unit_nr, '(T3,A)') '**                                                                          **'
         WRITE (unit_nr, '(T3,A)') '******************************************************************************'
         WRITE (unit_nr, '(T3,A)') ' '
         WRITE (unit_nr, '(T3,A)') ' '
         WRITE (unit_nr, '(T3,A)') ' The excitation energies are calculated by fully diagonalizing: '
         WRITE (unit_nr, '(T3,A)') ' '
         WRITE (unit_nr, '(T3,A)') '    A_iajb   =  (E_a-E_i) delta_ij delta_ab   +  alpha * v_iajb   -  W_ijab   '
         IF (.NOT. flag_TDA) THEN
            WRITE (unit_nr, '(T3,A)') '    B_iajb   =  alpha * v_iajb   -  W_ibaj   '
         END IF
         WRITE (unit_nr, '(T3,A)') ' '
         WRITE (unit_nr, '(T3,A48,A7,A12,F3.1)') &
            ' The spin-dependent factor is for the requested ', multiplet, " is alpha = ", alpha
         WRITE (unit_nr, '(T3,A)') ' '
         WRITE (unit_nr, '(T3,A16,T50,A22)') &
            ' Excitonic level', 'Excitation energy (eV)'
         !prints actual energies values
         DO i = 1, MIN(15, num_exc_en)
            WRITE (unit_nr, '(T3,I16,T50,F22.3)') i, Exc_ens(i)*evolt
         END DO
         WRITE (unit_nr, '(T3,A)') ' '
      END IF

      CALL timestop(handle)
   END SUBROUTINE success_message


   SUBROUTINE reduce_BSE_matrices(fm_mat_S, fm_mat_S_ij_bse, fm_mat_S_ab_bse,&
                                  fm_mat_S_reduced, Eigenval, homo, virtual, dimen_RI, &
                                  homo_red, virtual_red)

   TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat_S, fm_mat_S_ij_bse, fm_mat_S_ab_bse, &
                                                         fm_mat_S_reduced
      REAL(KIND=dp), DIMENSION(:)                     :: Eigenval
   INTEGER, INTENT(IN)                                :: homo, virtual, dimen_RI
   INTEGER, INTENT(OUT)                               :: homo_red, virtual_red

   CHARACTER(LEN=*), PARAMETER                        :: routineN = 'reduce_BSE_matrices'

   INTEGER :: handle


   CALL timeset(routineN, handle)



   CALL timestop(handle)

   END SUBROUTINE reduce_BSE_matrices

END MODULE bse_full_diag
