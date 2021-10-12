!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief collects routines that calculate density matrices
!> \note
!>      first version : most routines imported
!> \author JGH (2020-01)
! **************************************************************************************************
MODULE qs_density_matrices
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              cp_dbcsr_plus_fm_fm_t,&
                                              cp_dbcsr_sm_fm_multiply
   USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale,&
                                              cp_fm_scale_and_add,&
                                              cp_fm_symm,&
                                              cp_fm_transpose,&
                                              cp_fm_upper_to_full
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_to_fm,&
                                              cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_copy,&
                                              dbcsr_multiply,&
                                              dbcsr_release,&
                                              dbcsr_scale_by_vector,&
                                              dbcsr_set,&
                                              dbcsr_type
   USE kinds,                           ONLY: dp
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_density_matrices'

   PUBLIC :: calculate_density_matrix
   PUBLIC :: calculate_w_matrix, calculate_w_matrix_ot
   PUBLIC :: calculate_wz_matrix, calculate_whz_matrix
   PUBLIC :: calculate_wx_matrix, calculate_xwx_matrix

   INTERFACE calculate_density_matrix
      MODULE PROCEDURE calculate_dm_sparse
   END INTERFACE

   INTERFACE calculate_w_matrix
      MODULE PROCEDURE calculate_w_matrix_1, calculate_w_matrix_roks
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief   Calculate the density matrix
!> \param mo_set ...
!> \param density_matrix ...
!> \param use_dbcsr ...
!> \param retain_sparsity ...
!> \date    06.2002
!> \par History
!>       - Fractional occupied orbitals (MK)
!> \author  Joost VandeVondele
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE calculate_dm_sparse(mo_set, density_matrix, use_dbcsr, retain_sparsity)

      TYPE(mo_set_type), POINTER                         :: mo_set
      TYPE(dbcsr_type), POINTER                          :: density_matrix
      LOGICAL, INTENT(IN), OPTIONAL                      :: use_dbcsr, retain_sparsity

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_dm_sparse', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      LOGICAL                                            :: my_retain_sparsity, my_use_dbcsr
      REAL(KIND=dp)                                      :: alpha
      TYPE(cp_fm_type), POINTER                          :: fm_tmp
      TYPE(dbcsr_type)                                   :: dbcsr_tmp

      CALL timeset(routineN, handle)

      my_use_dbcsr = .FALSE.
      IF (PRESENT(use_dbcsr)) my_use_dbcsr = use_dbcsr
      my_retain_sparsity = .TRUE.
      IF (PRESENT(retain_sparsity)) my_retain_sparsity = retain_sparsity
      IF (my_use_dbcsr) THEN
         IF (.NOT. ASSOCIATED(mo_set%mo_coeff_b)) THEN
            CPABORT("mo_coeff_b NOT ASSOCIATED")
         END IF
      END IF

      CALL dbcsr_set(density_matrix, 0.0_dp)

      IF (.NOT. mo_set%uniform_occupation) THEN ! not all orbitals 1..homo are equally occupied

         IF (my_use_dbcsr) THEN
            CALL dbcsr_copy(dbcsr_tmp, mo_set%mo_coeff_b)
            CALL dbcsr_scale_by_vector(dbcsr_tmp, mo_set%occupation_numbers(1:mo_set%homo), &
                                       side='right')
            CALL dbcsr_multiply("N", "T", 1.0_dp, mo_set%mo_coeff_b, dbcsr_tmp, &
                                1.0_dp, density_matrix, retain_sparsity=my_retain_sparsity, &
                                last_k=mo_set%homo)
            CALL dbcsr_release(dbcsr_tmp)
         ELSE
            NULLIFY (fm_tmp)
            CALL cp_fm_create(fm_tmp, mo_set%mo_coeff%matrix_struct)
            CALL cp_fm_to_fm(mo_set%mo_coeff, fm_tmp)
            CALL cp_fm_column_scale(fm_tmp, mo_set%occupation_numbers(1:mo_set%homo))
            alpha = 1.0_dp
            CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=density_matrix, &
                                       matrix_v=mo_set%mo_coeff, &
                                       matrix_g=fm_tmp, &
                                       ncol=mo_set%homo, &
                                       alpha=alpha)
            CALL cp_fm_release(fm_tmp)
         ENDIF
      ELSE
         IF (my_use_dbcsr) THEN
            CALL dbcsr_multiply("N", "T", mo_set%maxocc, mo_set%mo_coeff_b, mo_set%mo_coeff_b, &
                                1.0_dp, density_matrix, retain_sparsity=my_retain_sparsity, &
                                last_k=mo_set%homo)
         ELSE
            alpha = mo_set%maxocc
            CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=density_matrix, &
                                       matrix_v=mo_set%mo_coeff, &
                                       ncol=mo_set%homo, &
                                       alpha=alpha)
         ENDIF
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE calculate_dm_sparse

! **************************************************************************************************
!> \brief Calculate the W matrix from the MO eigenvectors, MO eigenvalues,
!>       and the MO occupation numbers. Only works if they are eigenstates
!> \param mo_set type containing the full matrix of the MO and the eigenvalues
!> \param w_matrix sparse matrix
!>        error
!> \par History
!>         Creation (03.03.03,MK)
!>         Modification that computes it as a full block, several times (e.g. 20)
!>               faster at the cost of some additional memory
!> \author MK
! **************************************************************************************************
   SUBROUTINE calculate_w_matrix_1(mo_set, w_matrix)

      TYPE(mo_set_type), POINTER                         :: mo_set
      TYPE(dbcsr_type), POINTER                          :: w_matrix

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_w_matrix_1', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, imo
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigocc
      TYPE(cp_fm_type), POINTER                          :: weighted_vectors

      CALL timeset(routineN, handle)
      NULLIFY (weighted_vectors)

      CALL dbcsr_set(w_matrix, 0.0_dp)
      CALL cp_fm_create(weighted_vectors, mo_set%mo_coeff%matrix_struct, "weighted_vectors")
      CALL cp_fm_to_fm(mo_set%mo_coeff, weighted_vectors)

      ! scale every column with the occupation
      ALLOCATE (eigocc(mo_set%homo))

      DO imo = 1, mo_set%homo
         eigocc(imo) = mo_set%eigenvalues(imo)*mo_set%occupation_numbers(imo)
      ENDDO
      CALL cp_fm_column_scale(weighted_vectors, eigocc)
      DEALLOCATE (eigocc)

      CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=w_matrix, &
                                 matrix_v=mo_set%mo_coeff, &
                                 matrix_g=weighted_vectors, &
                                 ncol=mo_set%homo)

      CALL cp_fm_release(weighted_vectors)

      CALL timestop(handle)

   END SUBROUTINE calculate_w_matrix_1

! **************************************************************************************************
!> \brief Calculate the W matrix from the MO coefs, MO derivs
!>        could overwrite the mo_derivs for increased memory efficiency
!> \param mo_set type containing the full matrix of the MO coefs
!>        mo_deriv:
!> \param mo_deriv ...
!> \param w_matrix sparse matrix
!> \param s_matrix sparse matrix for the overlap
!>        error
!> \par History
!>         Creation (JV)
!> \author MK
! **************************************************************************************************
   SUBROUTINE calculate_w_matrix_ot(mo_set, mo_deriv, w_matrix, s_matrix)

      TYPE(mo_set_type), POINTER                         :: mo_set
      TYPE(dbcsr_type), POINTER                          :: mo_deriv, w_matrix, s_matrix

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_w_matrix_ot', &
         routineP = moduleN//':'//routineN
      LOGICAL, PARAMETER                                 :: check_gradient = .FALSE., &
                                                            do_symm = .FALSE.

      INTEGER                                            :: handle, iounit, ncol_block, ncol_global, &
                                                            nrow_block, nrow_global
      REAL(KIND=dp), DIMENSION(:), POINTER               :: occupation_numbers, scaling_factor
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_tmp
      TYPE(cp_fm_type), POINTER                          :: gradient, h_block, h_block_t, &
                                                            weighted_vectors
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)
      NULLIFY (weighted_vectors, h_block, fm_struct_tmp)

      CALL cp_fm_get_info(matrix=mo_set%mo_coeff, &
                          ncol_global=ncol_global, &
                          nrow_global=nrow_global, &
                          nrow_block=nrow_block, &
                          ncol_block=ncol_block)

      CALL cp_fm_create(weighted_vectors, mo_set%mo_coeff%matrix_struct, "weighted_vectors")
      CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=ncol_global, ncol_global=ncol_global, &
                               para_env=mo_set%mo_coeff%matrix_struct%para_env, &
                               context=mo_set%mo_coeff%matrix_struct%context)
      CALL cp_fm_create(h_block, fm_struct_tmp, name="h block")
      IF (do_symm) CALL cp_fm_create(h_block_t, fm_struct_tmp, name="h block t")
      CALL cp_fm_struct_release(fm_struct_tmp)

      CALL get_mo_set(mo_set=mo_set, occupation_numbers=occupation_numbers)
      ALLOCATE (scaling_factor(SIZE(occupation_numbers)))
      scaling_factor = 2.0_dp*occupation_numbers
      CALL copy_dbcsr_to_fm(mo_deriv, weighted_vectors)
      CALL cp_fm_column_scale(weighted_vectors, scaling_factor)
      DEALLOCATE (scaling_factor)

      ! the convention seems to require the half here, the factor of two is presumably taken care of
      ! internally in qs_core_hamiltonian
      CALL cp_gemm('T', 'N', ncol_global, ncol_global, nrow_global, 0.5_dp, &
                   mo_set%mo_coeff, weighted_vectors, 0.0_dp, h_block)

      IF (do_symm) THEN
         ! at the minimum things are anyway symmetric, but numerically it might not be the case
         ! needs some investigation to find out if using this is better
         CALL cp_fm_transpose(h_block, h_block_t)
         CALL cp_fm_scale_and_add(0.5_dp, h_block, 0.5_dp, h_block_t)
      ENDIF

      ! this could overwrite the mo_derivs to save the weighted_vectors
      CALL cp_gemm('N', 'N', nrow_global, ncol_global, ncol_global, 1.0_dp, &
                   mo_set%mo_coeff, h_block, 0.0_dp, weighted_vectors)

      CALL dbcsr_set(w_matrix, 0.0_dp)
      CALL cp_dbcsr_plus_fm_fm_t(sparse_matrix=w_matrix, &
                                 matrix_v=mo_set%mo_coeff, &
                                 matrix_g=weighted_vectors, &
                                 ncol=mo_set%homo)

      IF (check_gradient) THEN
         CALL cp_fm_create(gradient, mo_set%mo_coeff%matrix_struct, "gradient")
         CALL cp_dbcsr_sm_fm_multiply(s_matrix, weighted_vectors, &
                                      gradient, ncol_global)

         ALLOCATE (scaling_factor(SIZE(occupation_numbers)))
         scaling_factor = 2.0_dp*occupation_numbers
         CALL copy_dbcsr_to_fm(mo_deriv, weighted_vectors)
         CALL cp_fm_column_scale(weighted_vectors, scaling_factor)
         DEALLOCATE (scaling_factor)

         logger => cp_get_default_logger()
         IF (logger%para_env%ionode) THEN
            iounit = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
            WRITE (iounit, *) " maxabs difference ", &
               MAXVAL(ABS(weighted_vectors%local_data - 2.0_dp*gradient%local_data))
         ENDIF
         CALL cp_fm_release(gradient)
      ENDIF

      IF (do_symm) CALL cp_fm_release(h_block_t)
      CALL cp_fm_release(weighted_vectors)
      CALL cp_fm_release(h_block)

      CALL timestop(handle)

   END SUBROUTINE calculate_w_matrix_ot

! **************************************************************************************************
!> \brief Calculate the energy-weighted density matrix W if ROKS is active.
!>        The W matrix is returned in matrix_w.
!> \param mo_set ...
!> \param matrix_ks ...
!> \param matrix_p ...
!> \param matrix_w ...
!> \author 04.05.06,MK
! **************************************************************************************************
   SUBROUTINE calculate_w_matrix_roks(mo_set, matrix_ks, matrix_p, matrix_w)
      TYPE(mo_set_type), POINTER                         :: mo_set
      TYPE(dbcsr_type), POINTER                          :: matrix_ks, matrix_p, matrix_w

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_w_matrix_roks', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, nao
      TYPE(cp_blacs_env_type), POINTER                   :: context
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: c, ks, p, work
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      NULLIFY (c)
      NULLIFY (context)
      NULLIFY (fm_struct)
      NULLIFY (ks)
      NULLIFY (p)
      NULLIFY (para_env)
      NULLIFY (work)

      CALL get_mo_set(mo_set=mo_set, mo_coeff=c)
      CALL cp_fm_get_info(c, context=context, nrow_global=nao, para_env=para_env)
      CALL cp_fm_struct_create(fm_struct, context=context, nrow_global=nao, &
                               ncol_global=nao, para_env=para_env)
      CALL cp_fm_create(ks, fm_struct, name="Kohn-Sham matrix")
      CALL cp_fm_create(p, fm_struct, name="Density matrix")
      CALL cp_fm_create(work, fm_struct, name="Work matrix")
      CALL cp_fm_struct_release(fm_struct)
      CALL copy_dbcsr_to_fm(matrix_ks, ks)
      CALL copy_dbcsr_to_fm(matrix_p, p)
      CALL cp_fm_upper_to_full(p, work)
      CALL cp_fm_symm("L", "U", nao, nao, 1.0_dp, ks, p, 0.0_dp, work)
      CALL cp_gemm("T", "N", nao, nao, nao, 1.0_dp, p, work, 0.0_dp, ks)
      CALL dbcsr_set(matrix_w, 0.0_dp)
      CALL copy_fm_to_dbcsr(ks, matrix_w, keep_sparsity=.TRUE.)

      CALL cp_fm_release(work)
      CALL cp_fm_release(p)
      CALL cp_fm_release(ks)

      CALL timestop(handle)

   END SUBROUTINE calculate_w_matrix_roks

! **************************************************************************************************
!> \brief Calculate the response W matrix from the MO eigenvectors, MO eigenvalues,
!>       and the MO occupation numbers. Only works if they are eigenstates
!> \param mo_set type containing the full matrix of the MO and the eigenvalues
!> \param psi1 response orbitals
!> \param ks_matrix Kohn-Sham sparse matrix
!> \param w_matrix sparse matrix
!> \par History
!>               adapted from calculate_w_matrix_1
!> \author JGH
! **************************************************************************************************
   SUBROUTINE calculate_wz_matrix(mo_set, psi1, ks_matrix, w_matrix)

      TYPE(mo_set_type), POINTER                         :: mo_set
      TYPE(cp_fm_type), POINTER                          :: psi1
      TYPE(dbcsr_type), POINTER                          :: ks_matrix, w_matrix

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_wz_matrix', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ncol, nrow
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_tmp
      TYPE(cp_fm_type), POINTER                          :: ksmat, scrv

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(matrix=mo_set%mo_coeff, ncol_global=ncol, nrow_global=nrow)
      CALL cp_fm_create(scrv, mo_set%mo_coeff%matrix_struct, "scr vectors")
      CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=ncol, ncol_global=ncol, &
                               para_env=mo_set%mo_coeff%matrix_struct%para_env, &
                               context=mo_set%mo_coeff%matrix_struct%context)
      CALL cp_fm_create(ksmat, fm_struct_tmp, name="KS")
      CALL cp_fm_struct_release(fm_struct_tmp)
      CALL cp_dbcsr_sm_fm_multiply(ks_matrix, mo_set%mo_coeff, scrv, ncol)
      CALL cp_gemm("T", "N", ncol, ncol, nrow, 1.0_dp, mo_set%mo_coeff, scrv, 0.0_dp, ksmat)
      CALL cp_gemm("N", "N", nrow, ncol, ncol, 1.0_dp, mo_set%mo_coeff, ksmat, 0.0_dp, scrv)
      CALL dbcsr_set(w_matrix, 0.0_dp)
      CALL cp_dbcsr_plus_fm_fm_t(w_matrix, matrix_v=scrv, matrix_g=psi1, &
                                 ncol=mo_set%homo, alpha=0.5_dp)
      CALL cp_dbcsr_plus_fm_fm_t(w_matrix, matrix_v=psi1, matrix_g=scrv, &
                                 ncol=mo_set%homo, alpha=0.5_dp)
      CALL cp_fm_release(scrv)
      CALL cp_fm_release(ksmat)

      CALL timestop(handle)

   END SUBROUTINE calculate_wz_matrix

! **************************************************************************************************
!> \brief Calculate the Wz matrix from the MO eigenvectors, MO eigenvalues,
!>       and the MO occupation numbers. Only works if they are eigenstates
!> \param c0vec ...
!> \param hzm ...
!> \param w_matrix sparse matrix
!> \param focc ...
!> \par History
!>               adapted from calculate_w_matrix_1
!> \author JGH
! **************************************************************************************************
   SUBROUTINE calculate_whz_matrix(c0vec, hzm, w_matrix, focc)

      TYPE(cp_fm_type), POINTER                          :: c0vec
      TYPE(dbcsr_type), POINTER                          :: hzm, w_matrix
      REAL(KIND=dp), INTENT(IN)                          :: focc

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_whz_matrix', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, nao, norb
      REAL(KIND=dp)                                      :: falpha
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct, fm_struct_mat
      TYPE(cp_fm_type), POINTER                          :: chcmat, hcvec

      CALL timeset(routineN, handle)

      falpha = focc

      CALL cp_fm_create(hcvec, c0vec%matrix_struct, "hcvec")
      CALL cp_fm_get_info(hcvec, matrix_struct=fm_struct, nrow_global=nao, ncol_global=norb)
      CALL cp_fm_struct_create(fm_struct_mat, context=fm_struct%context, nrow_global=norb, &
                               ncol_global=norb, para_env=fm_struct%para_env)
      CALL cp_fm_create(chcmat, fm_struct_mat)
      CALL cp_fm_struct_release(fm_struct_mat)

      CALL cp_dbcsr_sm_fm_multiply(hzm, c0vec, hcvec, norb)
      CALL cp_gemm("T", "N", norb, norb, nao, 1.0_dp, c0vec, hcvec, 0.0_dp, chcmat)
      CALL cp_gemm("N", "N", nao, norb, norb, 1.0_dp, c0vec, chcmat, 0.0_dp, hcvec)

      CALL cp_dbcsr_plus_fm_fm_t(w_matrix, matrix_v=hcvec, matrix_g=c0vec, ncol=norb, alpha=falpha)

      CALL cp_fm_release(hcvec)
      CALL cp_fm_release(chcmat)

      CALL timestop(handle)

   END SUBROUTINE calculate_whz_matrix

! **************************************************************************************************
!> \brief Calculate the excited state W matrix from the MO eigenvectors, KS matrix
!> \param mos_occ ...
!> \param xvec ...
!> \param ks_matrix ...
!> \param w_matrix ...
!> \par History
!>               adapted from calculate_wz_matrix
!> \author JGH
! **************************************************************************************************
   SUBROUTINE calculate_wx_matrix(mos_occ, xvec, ks_matrix, w_matrix)

      TYPE(cp_fm_type), POINTER                          :: mos_occ, xvec
      TYPE(dbcsr_type), POINTER                          :: ks_matrix, w_matrix

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_wx_matrix', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ncol, nrow
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_tmp
      TYPE(cp_fm_type), POINTER                          :: ksmat, scrv

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(matrix=mos_occ, ncol_global=ncol, nrow_global=nrow)
      CALL cp_fm_create(scrv, mos_occ%matrix_struct, "scr vectors")
      CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=ncol, ncol_global=ncol, &
                               para_env=mos_occ%matrix_struct%para_env, &
                               context=mos_occ%matrix_struct%context)
      CALL cp_fm_create(ksmat, fm_struct_tmp, name="KS")
      CALL cp_fm_struct_release(fm_struct_tmp)
      CALL cp_dbcsr_sm_fm_multiply(ks_matrix, mos_occ, scrv, ncol)
      CALL cp_gemm("T", "N", ncol, ncol, nrow, 1.0_dp, mos_occ, scrv, 0.0_dp, ksmat)
      CALL cp_gemm("N", "N", nrow, ncol, ncol, 1.0_dp, xvec, ksmat, 0.0_dp, scrv)
      CALL cp_dbcsr_plus_fm_fm_t(w_matrix, matrix_v=scrv, matrix_g=xvec, ncol=ncol, alpha=0.5_dp)
      CALL cp_dbcsr_plus_fm_fm_t(w_matrix, matrix_v=xvec, matrix_g=scrv, ncol=ncol, alpha=0.5_dp)
      CALL cp_fm_release(scrv)
      CALL cp_fm_release(ksmat)

      CALL timestop(handle)

   END SUBROUTINE calculate_wx_matrix

! **************************************************************************************************
!> \brief Calculate the excited state W matrix from the MO eigenvectors, KS matrix
!> \param mos_occ ...
!> \param xvec ...
!> \param s_matrix ...
!> \param ks_matrix ...
!> \param w_matrix ...
!> \param eval ...
!> \par History
!>               adapted from calculate_wz_matrix
!> \author JGH
! **************************************************************************************************
   SUBROUTINE calculate_xwx_matrix(mos_occ, xvec, s_matrix, ks_matrix, w_matrix, eval)

      TYPE(cp_fm_type), POINTER                          :: mos_occ, xvec
      TYPE(dbcsr_type), POINTER                          :: s_matrix, ks_matrix, w_matrix
      REAL(KIND=dp), INTENT(IN)                          :: eval

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_xwx_matrix', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ncol, nrow
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_tmp
      TYPE(cp_fm_type), POINTER                          :: scrv, xsxmat

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(matrix=mos_occ, ncol_global=ncol, nrow_global=nrow)
      CALL cp_fm_create(scrv, mos_occ%matrix_struct, "scr vectors")
      CALL cp_fm_struct_create(fm_struct_tmp, nrow_global=ncol, ncol_global=ncol, &
                               para_env=mos_occ%matrix_struct%para_env, &
                               context=mos_occ%matrix_struct%context)
      CALL cp_fm_create(xsxmat, fm_struct_tmp, name="XSX")
      CALL cp_fm_struct_release(fm_struct_tmp)

      CALL cp_dbcsr_sm_fm_multiply(ks_matrix, xvec, scrv, ncol, 1.0_dp, 0.0_dp)
      CALL cp_dbcsr_sm_fm_multiply(s_matrix, xvec, scrv, ncol, eval, -1.0_dp)
      CALL cp_gemm("T", "N", ncol, ncol, nrow, 1.0_dp, xvec, scrv, 0.0_dp, xsxmat)

      CALL cp_gemm("N", "N", nrow, ncol, ncol, 1.0_dp, mos_occ, xsxmat, 0.0_dp, scrv)
      CALL cp_dbcsr_plus_fm_fm_t(w_matrix, matrix_v=scrv, matrix_g=mos_occ, ncol=ncol, alpha=0.5_dp)
      CALL cp_dbcsr_plus_fm_fm_t(w_matrix, matrix_v=mos_occ, matrix_g=scrv, ncol=ncol, alpha=0.5_dp)

      CALL cp_fm_release(scrv)
      CALL cp_fm_release(xsxmat)

      CALL timestop(handle)

   END SUBROUTINE calculate_xwx_matrix

END MODULE qs_density_matrices
