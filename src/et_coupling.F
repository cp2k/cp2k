!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief calculates the electron transfer coupling elements
!>      Wu, Van Voorhis, JCP 125, 164105 (2006)
!> \author fschiff (01.2007)
! **************************************************************************************************
MODULE et_coupling
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_basic_linalg,              ONLY: cp_fm_invert,&
                                              cp_fm_transpose
   USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                              fm_pool_get_el_struct
   USE cp_fm_struct,                    ONLY: cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type,&
                                              cp_to_string
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_p_type
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: diamat_all
   USE message_passing,                 ONLY: mp_sum
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_matrix_pools,                 ONLY: mpools_get
   USE qs_mo_types,                     ONLY: get_mo_set
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'et_coupling'
   LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .FALSE.

! *** Public subroutines ***

   PUBLIC :: calc_et_coupling

CONTAINS
! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE calc_et_coupling(qs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'calc_et_coupling'

      INTEGER                                            :: handle, i, iw, j, k, nao, ncol_local, &
                                                            nmo, nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: col_indices, row_indices
      LOGICAL                                            :: is_spin_constraint
      REAL(KIND=dp)                                      :: Sda, strength, Waa, Wbb, Wda
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: a, b, S_det
      REAL(KIND=dp), DIMENSION(2)                        :: eigenv
      REAL(KIND=dp), DIMENSION(2, 2)                     :: S_mat, tmp_mat, U, W_mat
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: rest_MO
      TYPE(cp_fm_pool_p_type), DIMENSION(:), POINTER     :: mo_mo_fm_pools
      TYPE(cp_fm_struct_type), POINTER                   :: mo_mo_fmstruct
      TYPE(cp_fm_type), POINTER                          :: inverse_mat, SMO, Tinverse, tmp2
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(section_vals_type), POINTER                   :: et_coupling_section

      NULLIFY (tmp2, SMO, rest_MO, mo_mo_fmstruct, rest_MO, energy, matrix_s, dft_control, para_env)

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      et_coupling_section => section_vals_get_subs_vals(qs_env%input, &
                                                        "PROPERTIES%ET_COUPLING")

      CALL get_qs_env(qs_env, dft_control=dft_control, para_env=para_env)

      iw = cp_print_key_unit_nr(logger, et_coupling_section, "PROGRAM_RUN_INFO", &
                                extension=".log")

      is_spin_constraint = .FALSE.
      ALLOCATE (rest_MO(2))
      ALLOCATE (a(dft_control%nspins))
      ALLOCATE (b(dft_control%nspins))
      ALLOCATE (S_det(dft_control%nspins))

      CALL mpools_get(qs_env%mpools, mo_mo_fm_pools=mo_mo_fm_pools)
      mo_mo_fmstruct => fm_pool_get_el_struct(mo_mo_fm_pools(1)%pool)
      DO i = 1, dft_control%nspins
         mo_mo_fmstruct => fm_pool_get_el_struct(mo_mo_fm_pools(i)%pool)

         CALL get_mo_set(mo_set=qs_env%mos(i)%mo_set, &
                         nao=nao, &
                         nmo=nmo)

         CALL cp_fm_create(matrix=tmp2, &
                           matrix_struct=qs_env%mos(i)%mo_set%mo_coeff%matrix_struct, &
                           name="ET_TMP"//TRIM(ADJUSTL(cp_to_string(2)))//"MATRIX")
         CALL cp_fm_create(matrix=inverse_mat, &
                           matrix_struct=mo_mo_fmstruct, &
                           name="INVERSE"//TRIM(ADJUSTL(cp_to_string(2)))//"MATRIX")
         CALL cp_fm_create(matrix=Tinverse, &
                           matrix_struct=mo_mo_fmstruct, &
                           name="T_INVERSE"//TRIM(ADJUSTL(cp_to_string(2)))//"MATRIX")
         CALL cp_fm_create(matrix=SMO, &
                           matrix_struct=mo_mo_fmstruct, &
                           name="ET_SMO"//TRIM(ADJUSTL(cp_to_string(1)))//"MATRIX")
         DO j = 1, 2
            NULLIFY (rest_MO(j)%matrix)
            CALL cp_fm_create(matrix=rest_MO(j)%matrix, &
                              matrix_struct=mo_mo_fmstruct, &
                              name="ET_rest_MO"//TRIM(ADJUSTL(cp_to_string(j)))//"MATRIX")
         END DO

!   calculate MO-overlap

         CALL get_qs_env(qs_env, matrix_s=matrix_s)
         CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, qs_env%et_coupling%et_mo_coeff(i)%matrix, &
                                      tmp2, nmo, 1.0_dp, 0.0_dp)
         CALL cp_gemm('T', 'N', nmo, nmo, nao, 1.0_dp, &
                      qs_env%mos(i)%mo_set%mo_coeff, &
                      tmp2, 0.0_dp, SMO)

!    calculate the MO-representation of the restraint matrix A

         CALL cp_dbcsr_sm_fm_multiply(qs_env%et_coupling%rest_mat(1)%matrix, &
                                      qs_env%et_coupling%et_mo_coeff(i)%matrix, &
                                      tmp2, nmo, 1.0_dp, 0.0_dp)

         CALL cp_gemm('T', 'N', nmo, nmo, nao, 1.0_dp, &
                      qs_env%mos(i)%mo_set%mo_coeff, &
                      tmp2, 0.0_dp, rest_MO(1)%matrix)

!    calculate the MO-representation of the restraint matrix D

         CALL cp_dbcsr_sm_fm_multiply(qs_env%et_coupling%rest_mat(2)%matrix, &
                                      qs_env%mos(i)%mo_set%mo_coeff, &
                                      tmp2, nmo, 1.0_dp, 0.0_dp)

         CALL cp_gemm('T', 'N', nmo, nmo, nao, 1.0_dp, &
                      qs_env%et_coupling%et_mo_coeff(i)%matrix, &
                      tmp2, 0.0_dp, rest_MO(2)%matrix)

         CALL cp_fm_invert(SMO, inverse_mat, S_det(i))

         CALL cp_fm_get_info(inverse_mat, nrow_local=nrow_local, ncol_local=ncol_local, &
                             row_indices=row_indices, col_indices=col_indices)
         b(i) = 0.0_dp

         DO j = 1, ncol_local
            DO k = 1, nrow_local
               b(i) = b(i) + rest_MO(2)%matrix%local_data(k, j)*inverse_mat%local_data(k, j)
            END DO
         END DO

         CALL cp_fm_transpose(inverse_mat, Tinverse)
         a(i) = 0.0_dp
         DO j = 1, ncol_local
            DO k = 1, nrow_local
               a(i) = a(i) + rest_MO(1)%matrix%local_data(k, j)*Tinverse%local_data(k, j)
            END DO
         END DO
         IF (is_spin_constraint) THEN
            a(i) = -a(i)
            b(i) = -b(i)
         END IF
         CALL mp_sum(a(i), para_env%group)

         CALL mp_sum(b(i), para_env%group)

         CALL cp_fm_release(tmp2)
         CALL cp_fm_release(rest_MO(1)%matrix)
         CALL cp_fm_release(rest_MO(2)%matrix)
         CALL cp_fm_release(SMO)
         CALL cp_fm_release(Tinverse)
         CALL cp_fm_release(inverse_mat)
      END DO
      DEALLOCATE (rest_MO)

!    solve eigenstates for the projector matrix

      IF (dft_control%nspins == 2) THEN
         Sda = S_det(1)*S_det(2)
         Wda = ((a(1) + a(2)) + (b(1) + b(2)))*0.5_dp*Sda
      ELSE
         Sda = S_det(1)**2
         Wda = (a(1) + b(1))*Sda
      END IF

      IF (dft_control%qs_control%ddapc_restraint) THEN
         Waa = qs_env%et_coupling%order_p
         Wbb = dft_control%qs_control%ddapc_restraint_control(1)%ddapc_restraint_control%ddapc_order_p
         strength = dft_control%qs_control%ddapc_restraint_control(1)%ddapc_restraint_control%strength
      END IF

!!   construct S and W   !!!
      S_mat(1, 1) = 1.0_dp
      S_mat(2, 2) = 1.0_dp
      S_mat(2, 1) = Sda
      S_mat(1, 2) = Sda

      W_mat(1, 1) = Wbb
      W_mat(2, 2) = Waa
      W_mat(2, 1) = Wda
      W_mat(1, 2) = Wda

!!  solve WC=SCN
      CALL diamat_all(S_mat, eigenv, .TRUE.)
      ! U = S**(-1/2)
      U = 0.0_dp
      U(1, 1) = 1.0_dp/SQRT(eigenv(1))
      U(2, 2) = 1.0_dp/SQRT(eigenv(2))
      tmp_mat = MATMUL(U, TRANSPOSE(S_mat))
      U = MATMUL(S_mat, tmp_mat)
      tmp_mat = MATMUL(W_mat, U)
      W_mat = MATMUL(U, tmp_mat)
      CALL diamat_all(W_mat, eigenv, .TRUE.)
      tmp_mat = MATMUL(U, W_mat)

      CALL get_qs_env(qs_env, energy=energy)
      W_mat(1, 1) = energy%total
      W_mat(2, 2) = qs_env%et_coupling%energy
      a(1) = (energy%total + strength*Wbb)*Sda - strength*Wda
      a(2) = (qs_env%et_coupling%energy + qs_env%et_coupling%e1*Waa)*Sda - qs_env%et_coupling%e1*Wda
      W_mat(1, 2) = (a(1) + a(2))*0.5_dp
      W_mat(2, 1) = W_mat(1, 2)

      S_mat = MATMUL(W_mat, (tmp_mat))
      W_mat = MATMUL(TRANSPOSE(tmp_mat), S_mat)

      IF (iw > 0) THEN
         WRITE (iw, *)
         WRITE (iw, '(T3,A,T60,(3X,F12.6))') 'Strength of constraint A          :', qs_env%et_coupling%e1
         WRITE (iw, '(T3,A,T60,(3X,F12.6))') 'Strength of constraint B          :', strength
         WRITE (iw, '(T3,A,T60,(3X,F12.6))') 'Final target value of constraint A:', Waa
         WRITE (iw, '(T3,A,T60,(3X,F12.6))') 'Final target value of constraint B:', Wbb
         WRITE (iw, *)
         WRITE (iw, '(T3,A,T60,(3X,F12.6))') &
            'Diabatic electronic coupling matrix element(mHartree):', ABS(W_mat(1, 2)*1000.0_dp)

      ENDIF

      CALL dbcsr_deallocate_matrix_set(qs_env%et_coupling%rest_mat)

      CALL cp_print_key_finished_output(iw, logger, et_coupling_section, &
                                        "PROGRAM_RUN_INFO")
      CALL timestop(handle)
   END SUBROUTINE calc_et_coupling

END MODULE et_coupling

