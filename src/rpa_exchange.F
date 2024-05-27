!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Auxiliary routines needed for RPA-exchange
!>        given blacs_env to another
!> \par History
!>      09.2016 created [Vladimir Rybkin]
!>      03.2019 Renamed [Frederick Stein]
!>      03.2019 Moved Functions from rpa_ri_gpw.F [Frederick Stein]
!>      04.2024 Added open-shell calculations, SOSEX [Frederick Stein]
!> \author Vladimir Rybkin
! **************************************************************************************************
MODULE rpa_exchange
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_fm_basic_linalg,              ONLY: cp_fm_column_scale
   USE cp_fm_diag,                      ONLY: choose_eigv_solver
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_to_fm_submat_general,&
                                              cp_fm_type
   USE cp_dbcsr_api,                       ONLY: &
        dbcsr_copy, dbcsr_create, dbcsr_dot, dbcsr_frobenius_norm, dbcsr_get_info, dbcsr_init_p, &
        dbcsr_multiply, dbcsr_p_type, dbcsr_release, dbcsr_set, dbcsr_trace, dbcsr_type, &
        dbcsr_type_no_symmetry
   USE group_dist_types,                ONLY: create_group_dist,&
                                              get_group_dist,&
                                              group_dist_d1_type,&
                                              group_dist_proc,&
                                              release_group_dist
   USE input_constants,                 ONLY: rpa_exchange_axk,&
                                              rpa_exchange_none,&
                                              rpa_exchange_sosex
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE local_gemm_api,                  ONLY: LOCAL_GEMM_PU_GPU
   USE mathconstants,                   ONLY: sqrthalf
   USE message_passing,                 ONLY: mp_para_env_release,&
                                              mp_para_env_type,&
                                              mp_proc_null
   USE mp2_types,                       ONLY: mp2_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE rpa_util,                        ONLY: calc_fm_mat_S_rpa,&
                                              remove_scaling_factor_rpa
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rpa_exchange'

   PUBLIC :: rpa_exchange_work_type, rpa_exchange_needed_mem

   TYPE rpa_exchange_work_type
      PRIVATE
      INTEGER :: exchange_correction = rpa_exchange_none
      INTEGER, DIMENSION(:), ALLOCATABLE :: homo, virtual, dimen_ia
      TYPE(group_dist_d1_type) :: aux_func_dist = group_dist_d1_type()
      INTEGER, DIMENSION(:), ALLOCATABLE :: aux2send
      INTEGER :: dimen_RI = 0
      INTEGER :: block_size = 0
      INTEGER :: group_size = 0
      INTEGER :: color_sub = 0
      INTEGER :: ngroup = 0
      REAL(KIND=dp) :: eps_filter = 0.0_dp
      TYPE(cp_fm_type) :: fm_mat_Q_tmp = cp_fm_type()
      TYPE(cp_fm_type) :: fm_mat_R_half_gemm = cp_fm_type()
      TYPE(cp_fm_type) :: fm_mat_U = cp_fm_type()
      TYPE(mp_para_env_type), POINTER :: para_env_sub => NULL()
   CONTAINS
      PROCEDURE, PASS(exchange_work), NON_OVERRIDABLE :: create => rpa_exchange_work_create
      PROCEDURE, PASS(exchange_work), NON_OVERRIDABLE :: compute => rpa_exchange_work_compute
      PROCEDURE, PASS(exchange_work), NON_OVERRIDABLE :: release => rpa_exchange_work_release
      PROCEDURE, PASS(exchange_work), NON_OVERRIDABLE :: redistribute_into_subgroups
      PROCEDURE, PASS(exchange_work), NON_OVERRIDABLE :: redistribute_Gamma
   END TYPE

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param mp2_env ...
!> \param homo ...
!> \param virtual ...
!> \param dimen_RI ...
!> \param para_env ...
!> \param mem_per_rank ...
!> \param mem_per_repl ...
! **************************************************************************************************
   SUBROUTINE rpa_exchange_needed_mem(mp2_env, homo, virtual, dimen_RI, para_env, mem_per_rank, mem_per_repl)
      TYPE(mp2_type), INTENT(IN)                         :: mp2_env
      INTEGER, DIMENSION(:), INTENT(IN)                  :: homo, virtual
      INTEGER, INTENT(IN)                                :: dimen_RI
      TYPE(mp_para_env_type), INTENT(IN)                 :: para_env
      REAL(KIND=dp), INTENT(INOUT)                       :: mem_per_rank, mem_per_repl

      INTEGER                                            :: block_size

      ! We need the block size and if it is unknown, an upper bound
      block_size = mp2_env%ri_rpa%exchange_block_size
      IF (block_size <= 0) block_size = MAX(1, (dimen_RI + para_env%num_pe - 1)/para_env%num_pe)

      ! storage of product matrix (upper bound only as it depends on the square of the potential still unknown block size)
      mem_per_rank = mem_per_rank + REAL(MAXVAL(homo), KIND=dp)**2*block_size**2*8.0_dp/(1024_dp**2)

      ! work arrays R (2x) and U, copies of Gamma (2x), communication buffer (as expensive as Gamma)
      mem_per_repl = mem_per_repl + 3.0_dp*dimen_RI*dimen_RI*8.0_dp/(1024_dp**2) &
                     + 3.0_dp*MAXVAL(homo*virtual)*dimen_RI*8.0_dp/(1024_dp**2)
   END SUBROUTINE rpa_exchange_needed_mem

! **************************************************************************************************
!> \brief ...
!> \param exchange_work ...
!> \param mp2_env ...
!> \param dimen_RI ...
!> \param fm_mat_S ...
!> \param fm_mat_Q ...
!> \param fm_mat_Q_gemm ...
!> \param homo ...
!> \param virtual ...
! **************************************************************************************************
   SUBROUTINE rpa_exchange_work_create(exchange_work, mp2_env, dimen_RI, &
                                       fm_mat_S, fm_mat_Q, fm_mat_Q_gemm, homo, virtual)
      CLASS(rpa_exchange_work_type), INTENT(INOUT) :: exchange_work
      TYPE(mp2_type), INTENT(INOUT) :: mp2_env
      INTEGER, INTENT(IN) :: dimen_RI
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN) :: fm_mat_S
      TYPE(cp_fm_type), INTENT(IN) :: fm_mat_Q, fm_mat_Q_gemm
      INTEGER, DIMENSION(SIZE(fm_mat_S)), INTENT(IN) :: homo, virtual

      INTEGER :: nspins, group_size, aux_global, aux_local, my_process_row, proc
      INTEGER, DIMENSION(:), POINTER :: row_indices, aux_distribution_fm
      TYPE(cp_blacs_env_type), POINTER :: context

      exchange_work%exchange_correction = mp2_env%ri_rpa%exchange_correction

      IF (exchange_work%exchange_correction == rpa_exchange_none) RETURN

      ASSOCIATE (para_env => fm_mat_S(1)%matrix_struct%para_env)
         group_size = mp2_env%ri_rpa%exchange_group_size
! If automatic choice is requested, choose a suitable group size
         IF (mp2_env%ri_rpa%use_hfx_implementation) THEN
            IF (group_size <= 0) group_size = para_env%num_pe
         ELSE
            IF (group_size <= 0 .OR. dimen_RI*group_size < para_env%num_pe) &
               group_size = (para_env%num_pe - 1)/dimen_RI + 1
         END IF
         exchange_work%ngroup = para_env%num_pe/group_size
         exchange_work%color_sub = para_env%mepos/group_size
         ALLOCATE (exchange_work%para_env_sub)
         IF (group_size < para_env%num_pe) THEN
            CALL exchange_work%para_env_sub%from_split(para_env, exchange_work%color_sub)
         ELSE
            ! We use all processes of the subgroup, so let us just duplicate the communicator to save ressources
            ! Pointer assignment is disliked by the compiler
            CALL exchange_work%para_env_sub%from_dup(para_env)
         END IF
      END ASSOCIATE

      CALL cp_fm_get_info(fm_mat_S(1), row_indices=row_indices, nrow_locals=aux_distribution_fm, context=context)
      CALL context%get(my_process_row=my_process_row)

      CALL create_group_dist(exchange_work%aux_func_dist, exchange_work%ngroup, dimen_RI)
      ALLOCATE (exchange_work%aux2send(0:exchange_work%ngroup - 1))
      exchange_work%aux2send = 0
      DO aux_local = 1, aux_distribution_fm(my_process_row)
         aux_global = row_indices(aux_local)
         proc = group_dist_proc(exchange_work%aux_func_dist, aux_global)
         exchange_work%aux2send(proc) = exchange_work%aux2send(proc) + 1
      END DO

      exchange_work%eps_filter = mp2_env%mp2_gpw%eps_filter

      nspins = SIZE(fm_mat_S)

      ALLOCATE (exchange_work%homo(nspins), exchange_work%virtual(nspins), exchange_work%dimen_ia(nspins))
      exchange_work%homo(:) = homo
      exchange_work%virtual(:) = virtual
      exchange_work%dimen_ia(:) = homo*virtual
      exchange_work%dimen_RI = dimen_RI

      exchange_work%block_size = mp2_env%ri_rpa%exchange_block_size
      IF (exchange_work%block_size <= 0) exchange_work%block_size = dimen_RI

      CALL cp_fm_create(exchange_work%fm_mat_U, fm_mat_Q%matrix_struct, name="fm_mat_U")
      CALL cp_fm_create(exchange_work%fm_mat_Q_tmp, fm_mat_Q%matrix_struct, name="fm_mat_Q_tmp")
      CALL cp_fm_create(exchange_work%fm_mat_R_half_gemm, fm_mat_Q_gemm%matrix_struct)
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param exchange_work ...
! **************************************************************************************************
   SUBROUTINE rpa_exchange_work_release(exchange_work)
      CLASS(rpa_exchange_work_type), INTENT(INOUT) :: exchange_work

      IF (ALLOCATED(exchange_work%homo)) DEALLOCATE (exchange_work%homo)
      IF (ALLOCATED(exchange_work%virtual)) DEALLOCATE (exchange_work%virtual)
      IF (ALLOCATED(exchange_work%dimen_ia)) DEALLOCATE (exchange_work%dimen_ia)
      IF (ASSOCIATED(exchange_work%para_env_sub)) THEN
         CALL mp_para_env_release(exchange_work%para_env_sub)
      END IF
      CALL release_group_dist(exchange_work%aux_func_dist)
      IF (ALLOCATED(exchange_work%aux2send)) DEALLOCATE(exchange_work%aux2send)
      CALL cp_fm_release(exchange_work%fm_mat_Q_tmp)
      CALL cp_fm_release(exchange_work%fm_mat_U)
      CALL cp_fm_release(exchange_work%fm_mat_R_half_gemm)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Main driver for RPA-exchange energies
!> \param exchange_work ...
!> \param fm_mat_Q ...
!> \param eig ...
!> \param fm_mat_S ...
!> \param omega ...
!> \param e_exchange_corr ...          exchange energy correctrion for a quadrature point
!> \param mp2_env ...
!> \author Vladimir Rybkin, 07/2016
! **************************************************************************************************
   SUBROUTINE rpa_exchange_work_compute(exchange_work, fm_mat_Q, eig, fm_mat_S, omega, &
                                        e_exchange_corr, mp2_env)
      CLASS(rpa_exchange_work_type), INTENT(INOUT) :: exchange_work
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat_Q
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: eig
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: fm_mat_S
      REAL(KIND=dp), INTENT(IN)                          :: omega
      REAL(KIND=dp), INTENT(INOUT)                       :: e_exchange_corr
      TYPE(mp2_type), INTENT(INOUT) :: mp2_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'rpa_exchange_work_compute'
      REAL(KIND=dp), PARAMETER                           :: thresh = 0.0000001_dp

      INTEGER :: handle, iib, ispin, nspins, P, Q, L_size_Gamma, hom, virt, i, &
                 send_proc, recv_proc, recv_size, max_aux_size, proc_shift, dimen_ia, dimen_RI, &
                 block_size, P_start, P_end, P_size, Q_start, Q_size, Q_end, handle2, my_aux_size, my_virt
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenval
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), TARGET :: mat_Gamma_3_3D
      REAL(KIND=dp), POINTER, DIMENSION(:), CONTIGUOUS :: mat_Gamma_3_1D
      REAL(KIND=dp), POINTER, DIMENSION(:, :), CONTIGUOUS :: mat_Gamma_3_2D
      REAL(KIND=dp), ALLOCATABLE, TARGET, DIMENSION(:) :: recv_buffer_1D
      REAL(KIND=dp), POINTER, DIMENSION(:, :), CONTIGUOUS :: recv_buffer_2D
      REAL(KIND=dp), POINTER, DIMENSION(:, :, :), CONTIGUOUS :: recv_buffer_3D
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :) :: mat_B_iaP
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), TARGET :: product_matrix_1D
      REAL(KIND=dp), POINTER, DIMENSION(:, :), CONTIGUOUS :: product_matrix_2D
      REAL(KIND=dp), POINTER, DIMENSION(:, :, :, :), CONTIGUOUS :: product_matrix_4D
      TYPE(cp_fm_type)        :: fm_mat_Gamma_3
      TYPE(mp_para_env_type), POINTER :: para_env
      TYPE(group_dist_d1_type)                           :: virt_dist

      IF (exchange_work%exchange_correction == rpa_exchange_none) RETURN

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(fm_mat_Q, para_env=para_env, ncol_global=dimen_RI)

      nspins = SIZE(fm_mat_S)

      ! Eigenvalues
      ALLOCATE (eigenval(dimen_RI))
      eigenval = 0.0_dp

      CALL cp_fm_set_all(matrix=exchange_work%fm_mat_Q_tmp, alpha=0.0_dp)
      CALL cp_fm_set_all(matrix=exchange_work%fm_mat_U, alpha=0.0_dp)

      ! Copy Q to Q_tmp
      CALL cp_fm_to_fm(fm_mat_Q, exchange_work%fm_mat_Q_tmp)
      ! Diagonalize Q
      CALL choose_eigv_solver(exchange_work%fm_mat_Q_tmp, exchange_work%fm_mat_U, eigenval)

      ! Calculate diagonal matrix for R_half

      ! Manipulate eigenvalues to get diagonal matrix
      IF (exchange_work%exchange_correction == rpa_exchange_axk) THEN
         DO iib = 1, dimen_RI
            IF (ABS(eigenval(iib)) .GE. thresh) THEN
               eigenval(iib) = &
                  SQRT((1.0_dp/(eigenval(iib)**2))*LOG(1.0_dp + eigenval(iib)) &
                       - 1.0_dp/(eigenval(iib)*(eigenval(iib) + 1.0_dp)))
            ELSE
               eigenval(iib) = sqrthalf
            END IF
         END DO
      ELSE IF (exchange_work%exchange_correction == rpa_exchange_sosex) THEN
         DO iib = 1, dimen_RI
            IF (ABS(eigenval(iib)) .GE. thresh) THEN
               eigenval(iib) = &
                  SQRT(-(1.0_dp/(eigenval(iib)**2))*LOG(1.0_dp + eigenval(iib)) &
                       + 1.0_dp/eigenval(iib))
            ELSE
               eigenval(iib) = sqrthalf
            END IF
         END DO
      ELSE
         CPABORT("Unknown RPA exchange correction")
      END IF

      ! fm_mat_U now contains some sqrt of the required matrix-valued function
      CALL cp_fm_column_scale(exchange_work%fm_mat_U, eigenval)

      ! Release memory
      DEALLOCATE (eigenval)

      ! Redistribute fm_mat_U for "rectangular" multiplication: ia*P P*P
      CALL cp_fm_set_all(matrix=exchange_work%fm_mat_R_half_gemm, alpha=0.0_dp)

      CALL cp_fm_to_fm_submat_general(exchange_work%fm_mat_U, exchange_work%fm_mat_R_half_gemm, dimen_RI, &
                                      dimen_RI, 1, 1, 1, 1, exchange_work%fm_mat_U%matrix_struct%context)

      CALL get_group_dist(exchange_work%aux_func_dist, exchange_work%color_sub, sizes=my_aux_size)

      e_exchange_corr = 0.0_dp
      max_aux_size = (dimen_RI + para_env%num_pe - 1)/para_env%num_pe

      IF (.NOT. mp2_env%ri_rpa%use_hfx_implementation) THEN
         ! local_gemm_ctx has a very footprint the first time this routine is
         ! called.
         CALL mp2_env%local_gemm_ctx%create(LOCAL_GEMM_PU_GPU)
         CALL mp2_env%local_gemm_ctx%set_op_threshold_gpu(128*128*128*2)
      END IF

      DO ispin = 1, nspins
         hom = exchange_work%homo(ispin)
         virt = exchange_work%virtual(ispin)
         dimen_ia = hom*virt
         IF (hom < 1 .OR. virt < 1) CYCLE

         CALL cp_fm_create(fm_mat_Gamma_3, fm_mat_S(ispin)%matrix_struct)
         CALL cp_fm_set_all(matrix=fm_mat_Gamma_3, alpha=0.0_dp)

         ! Update G with a new value of Omega: in practice, it is G*S

         ! Scale fm_work_iaP
         CALL calc_fm_mat_S_rpa(fm_mat_S(ispin), .TRUE., virt, eig(:, ispin), &
                                hom, omega, 0.0_dp)

         ! Calculate Gamma_3: Gamma_3 = G*S*R^(1/2) = G*S*R^(1/2)
         CALL parallel_gemm(transa="T", transb="N", m=dimen_RI, n=dimen_ia, k=dimen_RI, alpha=1.0_dp, &
                            matrix_a=exchange_work%fm_mat_R_half_gemm, matrix_b=fm_mat_S(ispin), beta=0.0_dp, &
                            matrix_c=fm_mat_Gamma_3)

         IF (mp2_env%ri_rpa%use_hfx_implementation) THEN
! Redistribute fm_mat_Gamma_3

! Transform to AO basis

! Contract with exchange kernel
         ELSE
            CALL create_group_dist(virt_dist, exchange_work%para_env_sub%num_pe, virt)

            ! Remove extra factor from S after the multiplication (to return to the original matrix)
            CALL remove_scaling_factor_rpa(fm_mat_S(ispin), virt, eig(:, ispin), hom, omega)

            CALL exchange_work%redistribute_into_subgroups(fm_mat_Gamma_3, mat_Gamma_3_3D, ispin, virt_dist)
            CALL cp_fm_release(fm_mat_Gamma_3)

            mat_Gamma_3_1D(1:INT(virt, KIND=int_8)*hom*my_aux_size) => mat_Gamma_3_3D(:, :, 1:my_aux_size)
            mat_Gamma_3_2D(1:virt, 1:hom*my_aux_size) => mat_Gamma_3_1D(1:INT(virt, KIND=int_8)*hom*my_aux_size)

            ! We need only the pure matrix
            CALL remove_scaling_factor_rpa(fm_mat_S(ispin), virt, eig(:, ispin), hom, omega)

            ! Reorder matrix from (P, i*a) -> (a, i, P) with P being distributed within subgroups
            CALL exchange_work%redistribute_into_subgroups(fm_mat_S(ispin), mat_B_iaP, ispin, virt_dist)

            ! Return to the original tensor
            CALL calc_fm_mat_S_rpa(fm_mat_S(ispin), .TRUE., virt, eig(:, ispin), hom, omega, 0.0_dp)

            L_size_Gamma = SIZE(mat_Gamma_3_3D, 3)
            my_virt = SIZE(mat_Gamma_3_3D, 1)
            block_size = exchange_work%block_size

            ALLOCATE (product_matrix_1D(INT(hom*MIN(block_size, L_size_gamma), KIND=int_8)* &
                                        INT(hom*MIN(block_size, max_aux_size), KIND=int_8)))
            ALLOCATE (recv_buffer_1D(INT(virt, KIND=int_8)*hom*max_aux_size))
            recv_buffer_2D(1:my_virt, 1:hom*max_aux_size) => recv_buffer_1D(1:INT(virt, KIND=int_8)*hom*max_aux_size)
            recv_buffer_3D(1:my_virt, 1:hom, 1:max_aux_size) => recv_buffer_1D(1:INT(virt, KIND=int_8)*hom*max_aux_size)
            DO proc_shift = 0, para_env%num_pe - 1, exchange_work%para_env_sub%num_pe
               send_proc = MODULO(para_env%mepos + proc_shift, para_env%num_pe)
               recv_proc = MODULO(para_env%mepos - proc_shift, para_env%num_pe)

               CALL get_group_dist(exchange_work%aux_func_dist, recv_proc, sizes=recv_size)

               IF (recv_size == 0) recv_proc = mp_proc_null

               CALL para_env%sendrecv(mat_B_iaP, send_proc, recv_buffer_3D(:, :, 1:recv_size), recv_proc)

               IF (recv_size == 0) CYCLE

               DO P_start = 1, L_size_Gamma, block_size
                  P_end = MIN(L_size_Gamma, P_start + block_size - 1)
                  P_size = P_end - P_start + 1
                  DO Q_start = 1, recv_size, block_size
                     Q_end = MIN(recv_size, Q_start + block_size - 1)
                     Q_size = Q_end - Q_start + 1

                     ! Reassign product_matrix pointers to enforce contiguity of target array
                     product_matrix_2D(1:hom*P_size, 1:hom*Q_size) => &
                        product_matrix_1D(1:INT(hom*P_size, KIND=int_8)*INT(hom*Q_size, KIND=int_8))
                     product_matrix_4D(1:hom, 1:P_size, 1:hom, 1:Q_size) => &
                        product_matrix_1D(1:INT(hom*P_size, KIND=int_8)*INT(hom*Q_size, KIND=int_8))

                     CALL timeset(routineN//"_gemm", handle2)
                     CALL mp2_env%local_gemm_ctx%gemm("T", "N", hom*P_size, hom*Q_size, my_virt, 1.0_dp, &
                                                      mat_Gamma_3_2D(:, hom*(P_start - 1) + 1:hom*P_end), my_virt, &
                                                      recv_buffer_2D(:, hom*(Q_start - 1) + 1:hom*Q_end), my_virt, &
                                                      0.0_dp, product_matrix_2D, hom*P_size)
                     CALL timestop(handle2)

                     CALL timeset(routineN//"_energy", handle2)
!$OMP PARALLEL DO DEFAULT(NONE) SHARED(P_size, Q_size, hom, product_matrix_4D) &
!$OMP             COLLAPSE(3) REDUCTION(+: e_exchange_corr) PRIVATE(P, Q, i)
                     DO P = 1, P_size
                     DO Q = 1, Q_size
                     DO i = 1, hom
                       e_exchange_corr = e_exchange_corr + DOT_PRODUCT(product_matrix_4D(i, P, :, Q), product_matrix_4D(:, P, i, Q))
                     END DO
                     END DO
                     END DO
                     CALL timestop(handle2)
                  END DO
               END DO
            END DO

            CALL release_group_dist(virt_dist)
            DEALLOCATE (mat_B_iaP, mat_Gamma_3_3D, product_matrix_1D, recv_buffer_1D)
         END IF
      END DO

      IF (.NOT. mp2_env%ri_rpa%use_hfx_implementation) THEN
         CALL mp2_env%local_gemm_ctx%destroy()
      END IF
      IF (nspins == 2) e_exchange_corr = e_exchange_corr*2.0_dp
      IF (nspins == 1) e_exchange_corr = e_exchange_corr*4.0_dp

      CALL timestop(handle)

   END SUBROUTINE rpa_exchange_work_compute

! **************************************************************************************************
!> \brief ...
!> \param fm_mat ...
!> \param mat ...
!> \param ncols ...
!> \param nrows ...
!> \param aux_func_dist ...
!> \param virt_dist ...
!> \param para_env_sub ...
! **************************************************************************************************
   SUBROUTINE redistribute_into_subgroups(exchange_work, fm_mat, mat, ispin, virt_dist)
   CLASS(rpa_exchange_work_type), INTENT(IN) :: exchange_work
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT)                                     :: mat
      INTEGER, INTENT(IN)                                :: ispin
      TYPE(group_dist_d1_type), INTENT(IN)               :: virt_dist

      CHARACTER(LEN=*), PARAMETER :: routineN = 'redistribute_into_subgroups'

      INTEGER :: aux_counter, aux_global, aux_local, aux_proc, avirt, dimen_RI, handle, handle2, &
         ia_global, ia_local, iocc, max_number_recv, max_number_send, my_aux_end, my_aux_size, &
         my_aux_start, my_process_column, my_process_row, my_virt_end, my_virt_size, &
         my_virt_start, proc, proc_shift, recv_proc, send_proc, virt_counter, virt_proc, group_size
      INTEGER, ALLOCATABLE, DIMENSION(:) :: data2send, recv_col_indices, &
         recv_row_indices, send_aux_indices, send_indices, send_virt_indices, virt2send
      INTEGER, DIMENSION(2)                              :: recv_shape
      INTEGER, DIMENSION(:), POINTER                     :: aux_distribution_fm, col_indices, &
                                                            ia_distribution_fm, row_indices
      INTEGER, DIMENSION(:, :), POINTER                  :: mpi2blacs
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), TARGET   :: recv_buffer, send_buffer
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         POINTER                                         :: recv_ptr, send_ptr
      TYPE(cp_blacs_env_type), POINTER                   :: context
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      CALL cp_fm_get_info(matrix=fm_mat, &
                          nrow_locals=aux_distribution_fm, &
                          col_indices=col_indices, &
                          row_indices=row_indices, &
                          ncol_locals=ia_distribution_fm, &
                          context=context, &
                          nrow_global=dimen_RI, &
                          para_env=para_env)

      IF (exchange_work%homo(ispin) <= 0 .OR. exchange_work%virtual(ispin) <= 0) THEN
         ALLOCATE (mat(exchange_work%homo(ispin), exchange_work%virtual(ispin), dimen_RI))
         CALL timestop(handle)
         RETURN
      END IF

      group_size = exchange_work%para_env_sub%num_pe

      CALL timeset(routineN//"_prep", handle2)
      CALL get_group_dist(exchange_work%aux_func_dist, exchange_work%color_sub, my_aux_start, my_aux_end, my_aux_size)
      CALL get_group_dist(virt_dist, exchange_work%para_env_sub%mepos, my_virt_start, my_virt_end, my_virt_size)
      CALL context%get(my_process_column=my_process_column, my_process_row=my_process_row, mpi2blacs=mpi2blacs)

      ! Determine the number of columns to send
      ALLOCATE (send_aux_indices(MAXVAL(exchange_work%aux2send)))
      ALLOCATE (virt2send(0:group_size - 1))
      virt2send = 0
      DO ia_local = 1, ia_distribution_fm(my_process_column)
         ia_global = col_indices(ia_local)
         avirt = MOD(ia_global - 1, exchange_work%virtual(ispin)) + 1
         proc = group_dist_proc(virt_dist, avirt)
         virt2send(proc) = virt2send(proc) + 1
      END DO
      ALLOCATE (send_virt_indices(MAXVAL(virt2send)))
      max_number_send = MAXVAL(exchange_work%aux2send)*MAXVAL(virt2send)

      ALLOCATE (send_buffer(INT(max_number_send, KIND=int_8)*exchange_work%homo(ispin)), send_indices(max_number_send))

      ALLOCATE (data2send(0:para_env%num_pe - 1))
      DO aux_proc = 0, exchange_work%ngroup - 1
      DO virt_proc = 0, group_size - 1
         data2send(aux_proc*group_size + virt_proc) = exchange_work%aux2send(aux_proc)*virt2send(virt_proc)
      END DO
      END DO

      max_number_recv = exchange_work%homo(ispin)*MAXVAL(data2send)
      CALL para_env%max(max_number_recv)
      ALLOCATE (recv_buffer(max_number_recv))

      ALLOCATE (mat(my_virt_size, exchange_work%homo(ispin), my_aux_size))

      CALL timestop(handle2)

      CALL timeset(routineN//"_own", handle2)
      ! Start with own data
      DO aux_local = 1, aux_distribution_fm(my_process_row)
         aux_global = row_indices(aux_local)
         IF (aux_global < my_aux_start .OR. aux_global > my_aux_end) CYCLE
         DO ia_local = 1, ia_distribution_fm(my_process_column)
            ia_global = fm_mat%matrix_struct%col_indices(ia_local)

            iocc = (ia_global - 1)/exchange_work%virtual(ispin) + 1
            avirt = MOD(ia_global - 1, exchange_work%virtual(ispin)) + 1

            IF (my_virt_start > avirt .OR. my_virt_end < avirt) CYCLE

            mat(avirt - my_virt_start + 1, iocc, aux_global - my_aux_start + 1) = fm_mat%local_data(aux_local, ia_local)
         END DO
      END DO
      CALL timestop(handle2)

      DO proc_shift = 1, para_env%num_pe - 1
         send_proc = MODULO(para_env%mepos + proc_shift, para_env%num_pe)
         recv_proc = MODULO(para_env%mepos - proc_shift, para_env%num_pe)

         CALL timeset(routineN//"_pack_buffer", handle2)
         send_ptr(1:virt2send(MOD(send_proc, group_size)), &
         1:exchange_work%aux2send(send_proc/group_size)) => &
            send_buffer(1:INT(virt2send(MOD(send_proc, group_size)), KIND=int_8)*&
            exchange_work%aux2send(send_proc/group_size))
! Pack send buffer
         aux_counter = 0
         DO aux_local = 1, aux_distribution_fm(my_process_row)
            aux_global = row_indices(aux_local)
            proc = group_dist_proc(exchange_work%aux_func_dist, aux_global)
            IF (proc /= send_proc/group_size) CYCLE
            aux_counter = aux_counter + 1
            virt_counter = 0
            DO ia_local = 1, ia_distribution_fm(my_process_column)
               ia_global = col_indices(ia_local)
               avirt = MOD(ia_global - 1, exchange_work%virtual(ispin)) + 1

               proc = group_dist_proc(virt_dist, avirt)
               IF (proc /= MOD(send_proc, group_size)) CYCLE
               virt_counter = virt_counter + 1
               send_ptr(virt_counter, aux_counter) = fm_mat%local_data(aux_local, ia_local)
               send_virt_indices(virt_counter) = ia_global
            END DO
            send_aux_indices(aux_counter) = aux_global
         END DO
         CALL timestop(handle2)

         CALL timeset(routineN//"_ex_size", handle2)
         recv_shape = [1, 1]
         CALL para_env%sendrecv(SHAPE(send_ptr), send_proc, recv_shape, recv_proc)
         CALL timestop(handle2)

         IF (SIZE(send_ptr) == 0) send_proc = mp_proc_null
         IF (PRODUCT(recv_shape) == 0) recv_proc = mp_proc_null

         CALL timeset(routineN//"_ex_idx", handle2)
         ALLOCATE (recv_row_indices(recv_shape(1)), recv_col_indices(recv_shape(2)))
         CALL para_env%sendrecv(send_virt_indices(1:virt_counter), send_proc, recv_row_indices, recv_proc)
         CALL para_env%sendrecv(send_aux_indices(1:aux_counter), send_proc, recv_col_indices, recv_proc)
         CALL timestop(handle2)

         ! Prepare pointer to recv buffer (consider transposition while packing the send buffer)
         recv_ptr(1:recv_shape(1), 1:MAX(1, recv_shape(2))) => recv_buffer(1:recv_shape(1)*MAX(1, recv_shape(2)))

         CALL timeset(routineN//"_sendrecv", handle2)
! Perform communication
         CALL para_env%sendrecv(send_ptr, send_proc, recv_ptr, recv_proc)
         CALL timestop(handle2)

         IF (recv_proc == mp_proc_null) THEN
            DEALLOCATE (recv_row_indices, recv_col_indices)
            CYCLE
         END IF

         CALL timeset(routineN//"_unpack", handle2)
! Unpack receive buffer
         DO aux_local = 1, SIZE(recv_col_indices)
            aux_global = recv_col_indices(aux_local)

            DO ia_local = 1, SIZE(recv_row_indices)
               ia_global = recv_row_indices(ia_local)

               iocc = (ia_global - 1)/exchange_work%virtual(ispin) + 1
               avirt = MOD(ia_global - 1, exchange_work%virtual(ispin)) + 1

               mat(avirt, iocc, aux_global - my_aux_start + 1) = recv_ptr(ia_local, aux_local)
            END DO
         END DO
         CALL timestop(handle2)

         IF (ALLOCATED(recv_row_indices)) DEALLOCATE (recv_row_indices)
         IF (ALLOCATED(recv_col_indices)) DEALLOCATE (recv_col_indices)
      END DO

      DEALLOCATE (send_aux_indices, send_virt_indices)

      CALL timestop(handle)

   END SUBROUTINE redistribute_into_subgroups

! **************************************************************************************************
!> \brief ...
!> \param fm_mat_Gamma_3 ...
!> \param mat_Gamma_3 ...
!> \param para_env_sub ...
!> \param template_dbcsr ...
! **************************************************************************************************
   SUBROUTINE redistribute_Gamma(exchange_work, fm_mat_Gamma_3, mat_Gamma_3, template_dbcsr)
   CLASS(rpa_exchange_work_type), INTENT(IN) :: exchange_work
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat_Gamma_3
      TYPE(dbcsr_type), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: mat_Gamma_3
      TYPE(dbcsr_type), INTENT(IN)                       :: template_dbcsr

      !TYPE(cp_blacs_env), POINTER :: blacs_env_sub

MARK_USED(exchange_work)
      MARK_USED(fm_mat_Gamma_3)
      MARK_USED(mat_Gamma_3)
      MARK_USED(template_dbcsr)

      !CALL cp_blacs_env_create(blacs_env_sub, para_env_sub)

      ! Redistribute fm_mat_Gamma_3 into subgroups (P, ia) -> (a, i, P) (both fm)
      ! Assign P to different subgroups
      ! Distribute a, i within subgroups according to block-cyclic layout

      ! Copy fm to dbcsr
   END SUBROUTINE redistribute_Gamma

END MODULE rpa_exchange
