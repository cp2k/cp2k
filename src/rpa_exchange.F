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
                                              group_dist_proc
   USE input_constants,                 ONLY: rpa_exchange_axk,&
                                              rpa_exchange_none,&
                                              rpa_exchange_sosex
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE local_gemm_api,                  ONLY: LOCAL_GEMM_PU_GPU
   USE mathconstants,                   ONLY: sqrthalf
   USE message_passing,                 ONLY: mp_para_env_type
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
      INTEGER :: dimen_RI = 0
      INTEGER :: block_size = 0
      REAL(KIND=dp) :: eps_filter = 0.0_dp
      TYPE(cp_fm_type) :: fm_mat_Q_tmp = cp_fm_type()
      TYPE(cp_fm_type) :: fm_mat_R_half_gemm = cp_fm_type()
      TYPE(cp_fm_type) :: fm_mat_U = cp_fm_type()
   CONTAINS
      PROCEDURE, PASS(exchange_work), NON_OVERRIDABLE :: create => rpa_exchange_work_create
      PROCEDURE, PASS(exchange_work), NON_OVERRIDABLE :: compute => rpa_exchange_work_compute
      PROCEDURE, PASS(exchange_work), NON_OVERRIDABLE :: release => rpa_exchange_work_release
   END TYPE

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param mp2_env ...
!> \param homo ...
!> \param virtual ...
!> \param dimen_RI ...
!> \param global_comm ...
!> \param mem_per_rank ...
!> \param mem_per_repl ...
! **************************************************************************************************
   SUBROUTINE rpa_exchange_needed_mem(mp2_env, homo, virtual, dimen_RI, global_comm, mem_per_rank, mem_per_repl)
      TYPE(mp2_type), INTENT(IN)                         :: mp2_env
      INTEGER, DIMENSION(:), INTENT(IN)                  :: homo, virtual
      INTEGER, INTENT(IN)                                :: dimen_RI
      TYPE(mp_para_env_type), INTENT(IN)                 :: global_comm
      REAL(KIND=dp), INTENT(INOUT)                       :: mem_per_rank, mem_per_repl

      INTEGER                                            :: block_size

      ! We need the block size and if it is unknown, an upper bound
      block_size = mp2_env%ri_rpa%exchange_block_size
      IF (block_size <= 0) block_size = MAX(1, dimen_RI/global_comm%num_pe*mp2_env%ri_rpa%rpa_num_quad_points)

      ! storage of product matrix (upper bound only as it depends on the square of the potential still unknown block size)
      mem_per_rank = mem_per_rank + REAL(MAXVAL(homo), KIND=dp)**2*block_size**2

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

      INTEGER :: nspins

      exchange_work%exchange_correction = mp2_env%ri_rpa%exchange_correction

      IF (exchange_work%exchange_correction == rpa_exchange_none) RETURN

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
                 block_size, P_start, P_end, P_size, Q_start, Q_size, Q_end
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

      e_exchange_corr = 0.0_dp
      max_aux_size = (dimen_RI + para_env%num_pe - 1)/para_env%num_pe

      ! local_gemm_ctx has a very footprint the first time this routine is
      ! called.
      CALL mp2_env%local_gemm_ctx%create(LOCAL_GEMM_PU_GPU)
      CALL mp2_env%local_gemm_ctx%set_op_threshold_gpu(128*128*128*2)

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

         ! Remove extra factor from S after the multiplication (to return to the original matrix)
         CALL remove_scaling_factor_rpa(fm_mat_S(ispin), virt, eig(:, ispin), hom, omega)

         CALL redistribute_into_subgroups(fm_mat_Gamma_3, mat_Gamma_3_3D, hom, virt)
         CALL cp_fm_release(fm_mat_Gamma_3)

         mat_Gamma_3_1D(1:INT(virt, KIND=int_8)*hom*max_aux_size) => mat_Gamma_3_3D(:, :, 1:max_aux_size)
         mat_Gamma_3_2D(1:virt, 1:hom*max_aux_size) => mat_Gamma_3_1D(1:INT(virt, KIND=int_8)*hom*max_aux_size)

         ! We need only the pure matrix
         CALL remove_scaling_factor_rpa(fm_mat_S(ispin), virt, eig(:, ispin), hom, omega)

         ! Reorder matrix from (P, i*a) -> (a, i, P) with P being distributed within subgroups
         CALL redistribute_into_subgroups(fm_mat_S(ispin), mat_B_iaP, hom, virt)

         ! Return to the original tensor
         CALL calc_fm_mat_S_rpa(fm_mat_S(ispin), .TRUE., virt, eig(:, ispin), hom, omega, 0.0_dp)

         L_size_Gamma = SIZE(mat_Gamma_3_3D, 3)
         block_size = exchange_work%block_size

         ALLOCATE (product_matrix_1D(INT(hom*block_size, KIND=int_8)*INT(hom*block_size, KIND=int_8)))
         ALLOCATE (recv_buffer_1D(INT(virt, KIND=int_8)*hom*max_aux_size))
         recv_buffer_2D(1:virt, 1:hom*max_aux_size) => recv_buffer_1D(1:INT(virt, KIND=int_8)*hom*max_aux_size)
         recv_buffer_3D(1:virt, 1:hom, 1:max_aux_size) => recv_buffer_1D(1:INT(virt, KIND=int_8)*hom*max_aux_size)
         DO proc_shift = 0, para_env%num_pe - 1
            send_proc = MODULO(para_env%mepos + proc_shift, para_env%num_pe)
            recv_proc = MODULO(para_env%mepos - proc_shift, para_env%num_pe)

            IF (recv_proc < MOD(dimen_RI, para_env%num_pe)) THEN
               recv_size = (dimen_RI + para_env%num_pe - 1)/para_env%num_pe
            ELSE
               recv_size = dimen_RI/para_env%num_pe
            END IF

            CALL para_env%sendrecv(mat_B_iaP, send_proc, recv_buffer_3D(:, :, 1:recv_size), recv_proc)

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

                  CALL mp2_env%local_gemm_ctx%gemm("T", "N", hom*P_size, hom*Q_size, virt, 1.0_dp, &
                                                   mat_Gamma_3_2D(:, hom*(P_start - 1) + 1:hom*P_end), virt, &
                                                   recv_buffer_2D(:, hom*(Q_start - 1) + 1:hom*Q_end), virt, &
                                                   0.0_dp, product_matrix_2D, hom*P_size)
                  DO P = 1, P_size
                  DO Q = 1, Q_size
                  DO i = 1, hom
                     e_exchange_corr = e_exchange_corr + DOT_PRODUCT(product_matrix_4D(i, P, :, Q), product_matrix_4D(:, P, i, Q))
                  END DO
                  END DO
                  END DO
               END DO
            END DO
         END DO

         DEALLOCATE (mat_B_iaP, mat_Gamma_3_3D, product_matrix_1D, recv_buffer_1D)
      END DO
      IF (nspins == 2) e_exchange_corr = e_exchange_corr*2.0_dp
      IF (nspins == 1) e_exchange_corr = e_exchange_corr*4.0_dp

      CALL mp2_env%local_gemm_ctx%destroy()

      CALL timestop(handle)

   END SUBROUTINE rpa_exchange_work_compute

! **************************************************************************************************
!> \brief ...
!> \param fm_mat ...
!> \param mat ...
!> \param ncols ...
!> \param nrows ...
! **************************************************************************************************
   SUBROUTINE redistribute_into_subgroups(fm_mat, mat, ncols, nrows)
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_mat
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT)                                     :: mat
      INTEGER, INTENT(IN)                                :: ncols, nrows

      CHARACTER(LEN=*), PARAMETER :: routineN = 'redistribute_into_subgroups'

      INTEGER :: aux_global, aux_local, avirt, dimen_RI, handle, ia_global, ia_local, iocc, &
         max_number_recv, max_number_send, my_aux_end, my_aux_size, my_aux_start, &
         my_process_column, my_process_row, proc, proc_shift, recv_cols, recv_counter, recv_proc, &
         recv_proc_col, recv_proc_row, recv_rows, send_counter, send_proc
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: data2recv, data2send
      INTEGER, DIMENSION(:), POINTER                     :: aux_distribution_fm, col_indices, &
                                                            ia_distribution_fm, row_indices
      INTEGER, DIMENSION(:, :), POINTER                  :: mpi2blacs
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), TARGET   :: recv_buffer
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: send_buffer
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         POINTER                                         :: recv_ptr
      TYPE(cp_blacs_env_type), POINTER                   :: context
      TYPE(group_dist_d1_type)                           :: aux_func_dist
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

      IF (ncols <= 0 .OR. nrows <= 0) THEN
         ALLOCATE (mat(ncols, nrows, dimen_RI))
         CALL timestop(handle)
         RETURN
      END IF

      CALL create_group_dist(aux_func_dist, para_env%num_pe, dimen_RI)
      CALL get_group_dist(aux_func_dist, para_env%mepos, my_aux_start, my_aux_end, my_aux_size)

      CALL context%get(my_process_column=my_process_column, my_process_row=my_process_row, &
                       mpi2blacs=mpi2blacs)

      ! Determine the number of columns to send
      ALLOCATE (data2send(0:para_env%num_pe - 1))
      data2send = 0
      DO aux_local = 1, aux_distribution_fm(my_process_row)
         aux_global = row_indices(aux_local)
         proc = group_dist_proc(aux_func_dist, aux_global)
         data2send(proc) = data2send(proc) + 1
      END DO
      max_number_send = MAXVAL(data2send)

      ! Determine the number of Auxiliary functions and occ/virt pairs to recv
      ALLOCATE (data2recv(0:para_env%num_pe - 1))
      data2recv = 0
      DO aux_global = my_aux_start, my_aux_end
         proc = fm_mat%matrix_struct%g2p_row(aux_global)
         data2recv(proc) = data2recv(proc) + 1
      END DO
      max_number_recv = MAXVAL(data2recv*ia_distribution_fm(mpi2blacs(1, :)))

      ALLOCATE (send_buffer(ia_distribution_fm(my_process_column), max_number_send))
      ALLOCATE (recv_buffer(max_number_recv))

      ALLOCATE (mat(nrows, ncols, my_aux_size))

      ! Start with own data
      recv_counter = 0
      DO aux_local = 1, my_aux_size
         aux_global = row_indices(aux_local)
         IF (aux_global < my_aux_start .OR. aux_global > my_aux_end) CYCLE
         DO ia_local = 1, ia_distribution_fm(my_process_column)
            ia_global = fm_mat%matrix_struct%col_indices(ia_local)

            iocc = (ia_global - 1)/nrows + 1
            avirt = MOD(ia_global - 1, nrows) + 1

            mat(avirt, iocc, aux_global - my_aux_start + 1) = fm_mat%local_data(aux_local, ia_local)
         END DO

      END DO

      DO proc_shift = 0, para_env%num_pe - 1
         send_proc = MODULO(para_env%mepos + proc_shift, para_env%num_pe)
         recv_proc = MODULO(para_env%mepos - proc_shift, para_env%num_pe)

! Pack send buffer
         send_counter = 0
         DO aux_local = 1, aux_distribution_fm(my_process_row)
            aux_global = row_indices(aux_local)
            proc = group_dist_proc(aux_func_dist, aux_global)
            IF (proc /= send_proc) CYCLE
            send_counter = send_counter + 1
            send_buffer(:, send_counter) = fm_mat%local_data(aux_local, :)
         END DO

         ! Prepare pointer to recv buffer
         recv_proc_col = mpi2blacs(2, recv_proc)
         recv_rows = ia_distribution_fm(recv_proc_col)
         recv_cols = data2recv(recv_proc)
         recv_ptr(1:recv_rows, 1:recv_cols) => recv_buffer(1:recv_rows*recv_cols)

! Perform communication
         CALL para_env%sendrecv(send_buffer(:, 1:data2send(send_proc)), send_proc, &
                                recv_ptr, recv_proc)

! Unpack receive buffer
         recv_counter = 0
         DO aux_global = my_aux_start, my_aux_end
            proc = fm_mat%matrix_struct%g2p_row(aux_global)
            IF (proc /= recv_proc) CYCLE

            recv_counter = recv_counter + 1

            DO ia_local = 1, recv_rows
               ia_global = fm_mat%matrix_struct%l2g_col(ia_local, recv_proc_row)

               iocc = (ia_global - 1)/nrows + 1
               avirt = MOD(ia_global - 1, nrows) + 1

               mat(avirt, iocc, aux_global - my_aux_start + 1) = recv_ptr(ia_local, recv_counter)
            END DO

         END DO
      END DO

      CALL timestop(handle)

   END SUBROUTINE redistribute_into_subgroups

END MODULE rpa_exchange
