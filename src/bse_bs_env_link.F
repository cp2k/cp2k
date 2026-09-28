!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Bethe-Salpeter equation on the RI-AO GW of PROPERTIES%BANDSTRUCTURE%GW: reads
!>        PROPERTIES%BSE, builds the window-sized RI slabs and the effective inverse dielectric
!>        matrix from bs_env, and hands them to the solver core through its windowed seam.
!>        Non-periodic, Gamma point only. Reads no real-time propagation input.
!> \author Max Graml
! **************************************************************************************************
MODULE bse_bs_env_link
   USE bse_input,                       ONLY: read_bse_section
   USE bse_main,                        ONLY: bse_solve,&
                                              prepare_bse_env_windowed
   USE bse_matvec,                      ONLY: mem_fraction
   USE bse_types,                       ONLY: bse_env_release,&
                                              bse_env_type
   USE bse_util,                        ONLY: bse_window_from_dft
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_create, dbcsr_distribution_get, dbcsr_distribution_new, dbcsr_distribution_release, &
        dbcsr_distribution_type, dbcsr_finalize, dbcsr_get_stored_coordinates, dbcsr_put_block, &
        dbcsr_release, dbcsr_type, dbcsr_type_no_symmetry, dbcsr_work_create
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add
   USE cp_fm_diag,                      ONLY: cp_fm_power
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_submatrix,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm,&
                                              cp_fm_to_fm_submat,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE dbt_api,                         ONLY: &
        dbt_batched_contract_finalize, dbt_batched_contract_init, dbt_clear, dbt_contract, &
        dbt_create, dbt_destroy, dbt_filter, dbt_get_block, dbt_get_stored_coordinates, &
        dbt_iterator_blocks_left, dbt_iterator_next_block, dbt_iterator_start, dbt_iterator_stop, &
        dbt_iterator_type, dbt_pgrid_create, dbt_pgrid_destroy, dbt_pgrid_type, dbt_put_block, &
        dbt_reserve_blocks, dbt_type
   USE gw_tensor_large_cell_Gamma,      ONLY: compute_3c_integrals
   USE gw_utils_communication,          ONLY: local_matrix_to_global_matrix
   USE gw_utils_fm,                     ONLY: fm_contract_ABA
   USE input_constants,                 ONLY: bse_memcheck_off,&
                                              bse_memcheck_warn
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE libint_2c_3c,                    ONLY: libint_potential_type
   USE message_passing,                 ONLY: mp_mem_avail_per_rank_GB
   USE mo_window,                       ONLY: mo_window_type
   USE mp2_ri_2c,                       ONLY: RI_2c_integral_mat
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE post_scf_bandstructure_types,    ONLY: post_scf_bandstructure_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_tensors,                      ONLY: neighbor_list_3c_destroy
   USE qs_tensors_types,                ONLY: create_2c_tensor,&
                                              create_3c_tensor

!$ USE OMP_LIB, ONLY: omp_get_num_threads, omp_get_thread_num
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'bse_bs_env_link'

   ! MO columns per tensor block of the window coefficients
   INTEGER, PARAMETER, PRIVATE :: mo_chunk = 64
   ! eigenvalue floor of the RI Coulomb matrix for its symmetric roots; the GW inverts its metric
   ! with no floor on the RI-AO path (eps_eigval_mat_RI = 0), so the roots need their own
   REAL(KIND=dp), PARAMETER, PRIVATE :: eps_ri_coulomb = 1.0E-10_dp

   PUBLIC :: bse_from_bs_env

CONTAINS

! **************************************************************************************************
!> \brief Entry from the properties driver after the bandstructure block: checks, section, window,
!>        screening, per-spin slabs, seam, solve, release
!> \param qs_env ...
!> \param bse_section PROPERTIES%BSE
! **************************************************************************************************
   SUBROUTINE bse_from_bs_env(qs_env, bse_section)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: bse_section

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'bse_from_bs_env'

      INTEGER                                            :: handle, ispin, n_mo_real, n_RI, n_spin, &
                                                            n_window, unit_nr
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: homo_full, homo_red, virt_red
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: eigenval_dft, eigenval_qp
      TYPE(bse_env_type)                                 :: bse_env
      TYPE(cp_fm_type)                                   :: fm_eps_inv, fm_T
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_S_ab, fm_S_ia, fm_S_ij
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(mo_window_type)                               :: window
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      unit_nr = cp_logger_get_default_io_unit(logger)
      NULLIFY (bs_env)
      CALL get_qs_env(qs_env, bs_env=bs_env)
      CPASSERT(ASSOCIATED(bs_env))
      ! the requirement reader at GW init refused RI_RS and periodic runs already
      CPASSERT(bs_env%do_gw .AND. .NOT. bs_env%do_periodic .AND. .NOT. bs_env%do_gw_ri_rs)
      CPASSERT(bs_env%keep_W_freq_zero .AND. bs_env%keep_3c_nl)
      CPASSERT(ASSOCIATED(bs_env%fm_W_MIC_freq_zero%matrix_struct))
      CPASSERT(bs_env%nkp_bs_and_DOS == 1)
      ! every RI batch of the slabs goes to one tensor group; unset groups would skip them all silently
      CPASSERT(bs_env%num_tensor_groups >= 1 .AND. bs_env%tensor_group_color >= 0)

      n_spin = bs_env%n_spin
      n_RI = bs_env%n_RI
      ! the canonical orthogonalisation parks the linearly dependent directions above n_mo_retained
      n_mo_real = MIN(bs_env%n_mo_retained, bs_env%n_ao)

      CALL read_bse_section(bse_section, bse_env)

      CALL link_window(bs_env, bse_env, n_mo_real, window, homo_red, virt_red, unit_nr)
      n_window = window%last_mo - window%first_mo + 1

      CALL link_screening(qs_env, bs_env, unit_nr, fm_T, fm_eps_inv)
      ! the GW kept W(0) for the screening and nl_3c for the slabs; each goes once read
      CALL cp_fm_release(bs_env%fm_W_MIC_freq_zero)

      ALLOCATE (fm_S_ia(n_spin), fm_S_ij(n_spin), fm_S_ab(n_spin))
      ALLOCATE (eigenval_qp(n_window, n_spin), eigenval_dft(n_window, n_spin), homo_full(n_spin))
      DO ispin = 1, n_spin
         CALL link_slabs(qs_env, bs_env, bse_env, ispin, window, homo_red(ispin), virt_red(ispin), fm_T, &
                         fm_S_ij(ispin), fm_S_ia(ispin), fm_S_ab(ispin), unit_nr)
         eigenval_qp(:, ispin) = bs_env%eigenval_GW(window%first_mo:window%last_mo, 1, ispin)
         eigenval_dft(:, ispin) = bs_env%eigenval_scf_Gamma(window%first_mo:window%last_mo, ispin)
         homo_full(ispin) = bs_env%n_occ(ispin)
      END DO
      CALL cp_fm_release(fm_T)
      CALL neighbor_list_3c_destroy(bs_env%nl_3c)

      CALL prepare_bse_env_windowed(bse_env, fm_S_ia, fm_S_ij, fm_S_ab, fm_eps_inv, &
                                    eigenval_qp, eigenval_dft, homo_full, homo_red, virt_red, &
                                    n_RI, qs_env, unit_nr)
      CALL cp_fm_release(fm_eps_inv)
      CALL bse_solve(bse_env, bs_env%fm_mo_coeff_Gamma(1:n_spin), unit_nr)
      CALL bse_env_release(bse_env)

      CALL timestop(handle)

   END SUBROUTINE bse_from_bs_env

! **************************************************************************************************
!> \brief The active window on the DFT axis over the real (retained) MOs, one absolute range for
!>        every spin, and the occupied and virtual counts inside it per spin
!> \param bs_env ...
!> \param bse_env its cutoffs, read from the section already
!> \param n_mo_real MOs that are states, not linear-dependency placeholders
!> \param window first and last MO of the window
!> \param homo_red occupied levels in the window per spin, allocated here
!> \param virt_red virtual levels in the window per spin, allocated here
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE link_window(bs_env, bse_env, n_mo_real, window, homo_red, virt_red, unit_nr)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(bse_env_type), INTENT(IN)                     :: bse_env
      INTEGER, INTENT(IN)                                :: n_mo_real
      TYPE(mo_window_type), INTENT(OUT)                  :: window
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: homo_red, virt_red
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'link_window'

      INTEGER                                            :: handle, ispin, n_spin
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: homo, n_mo

      CALL timeset(routineN, handle)

      n_spin = bs_env%n_spin
      ALLOCATE (n_mo(n_spin), homo(n_spin), homo_red(n_spin), virt_red(n_spin))
      n_mo(:) = n_mo_real
      homo(:) = bs_env%n_occ(1:n_spin)
      CALL bse_window_from_dft(bs_env%eigenval_scf_Gamma(1:n_mo_real, 1:n_spin), n_mo, homo, &
                               bse_env%bse_cutoff_occ, bse_env%bse_cutoff_empty, window, warn=.TRUE.)
      DO ispin = 1, n_spin
         homo_red(ispin) = homo(ispin) - window%first_mo + 1
         virt_red(ispin) = window%last_mo - homo(ispin)
         IF (homo_red(ispin) < 1) CPABORT("PROPERTIES%BSE: ENERGY_CUTOFF_OCC leaves no occupied level")
         IF (virt_red(ispin) < 1) CPABORT("PROPERTIES%BSE: ENERGY_CUTOFF_EMPTY leaves no empty level")
      END DO
      IF (unit_nr > 0) WRITE (unit_nr, '(T2,A4,T7,A,T61,2I10)') 'BSE|', &
         'Active window on the DFT axis (first, last MO)', window%first_mo, window%last_mo
      DEALLOCATE (n_mo, homo)

      CALL timestop(handle)

   END SUBROUTINE link_window

! **************************************************************************************************
!> \brief The RI transform and the effective inverse dielectric matrix from the GW's static W,
!>        with V_PQ the bare Coulomb and M_PQ the RI metric in the RI basis:
!>        T_QP = sum_R [M^-1]_QR [V^1/2]_RP, so that the slabs S_P,pq = sum_Q (pq|Q) T_QP give
!>        sum_P S_P,pq S_P,rs = sum_QR (pq|Q) [M^-1 V M^-1]_QR (R|rs);
!>        with the stored W~_PQ = [M^-1 W^c(0) M^-1]_PQ,
!>        ε^-1_eff,PQ = δ_PQ + sum_RS [V^-1/2]_PR [M W~ M]_RS [V^-1/2]_SQ, so that
!>        sum_PQ S_P,pq ε^-1_eff,PQ S_Q,rs = sum_RS (pq|R) [M^-1 (V + W^c(0)) M^-1]_RS (S|rs),
!>        the full static W of the GW. Symmetric roots V^±1/2 drop RI directions under the threshold
!> \param qs_env ...
!> \param bs_env ...
!> \param unit_nr ...
!> \param fm_T T = M^-1 V^1/2, RI x RI, created here
!> \param fm_eps_inv ε^-1_eff, RI x RI, created here
! **************************************************************************************************
   SUBROUTINE link_screening(qs_env, bs_env, unit_nr, fm_T, fm_eps_inv)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: unit_nr
      TYPE(cp_fm_type), INTENT(OUT)                      :: fm_T, fm_eps_inv

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'link_screening'

      INTEGER                                            :: handle, n_dependent, n_RI
      REAL(KIND=dp)                                      :: threshold
      TYPE(cp_fm_type)                                   :: fm_V_msqrt, fm_V_sqrt, fm_Wc, fm_work
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: fm_M, fm_V
      TYPE(libint_potential_type)                        :: coulomb_op

      CALL timeset(routineN, handle)

      n_RI = bs_env%n_RI
      ! bare Coulomb V (coulomb_op keeps its default, the full 1/r) and the RI metric M
      CALL RI_2c_integral_mat(qs_env, fm_V, bs_env%fm_RI_RI, n_RI, coulomb_op, do_kpoints=.FALSE.)
      CALL RI_2c_integral_mat(qs_env, fm_M, bs_env%fm_RI_RI, n_RI, bs_env%ri_metric, do_kpoints=.FALSE.)

      ! V^1/2 and V^-1/2 by symmetric roots, near-null RI directions dropped by the threshold
      threshold = bs_env%eps_eigval_mat_RI
      IF (threshold <= 0.0_dp) threshold = eps_ri_coulomb
      CALL cp_fm_create(fm_V_sqrt, bs_env%fm_RI_RI%matrix_struct)
      CALL cp_fm_create(fm_V_msqrt, bs_env%fm_RI_RI%matrix_struct)
      CALL cp_fm_create(fm_work, bs_env%fm_RI_RI%matrix_struct)
      CALL cp_fm_to_fm(fm_V(1, 1), fm_V_sqrt)
      CALL cp_fm_power(fm_V_sqrt, fm_work, 0.5_dp, threshold, n_dependent)
      CALL cp_fm_to_fm(fm_V(1, 1), fm_V_msqrt)
      CALL cp_fm_power(fm_V_msqrt, fm_work, -0.5_dp, threshold, n_dependent)
      IF (unit_nr > 0 .AND. n_dependent > 0) WRITE (unit_nr, '(T2,A4,T7,A,T71,I10)') 'BSE|', &
         'RI Coulomb directions dropped by the eigenvalue threshold', n_dependent

      ! T = M^-1 V^1/2
      CALL cp_fm_create(fm_T, bs_env%fm_RI_RI%matrix_struct)
      CALL parallel_gemm("N", "N", n_RI, n_RI, n_RI, 1.0_dp, bs_env%fm_Minv_Gamma, fm_V_sqrt, &
                         0.0_dp, fm_T)

      ! ε^-1_eff = 1 + V^-1/2 (M W~ M) V^-1/2
      CALL cp_fm_create(fm_Wc, bs_env%fm_RI_RI%matrix_struct)
      CALL cp_fm_to_fm(bs_env%fm_W_MIC_freq_zero, fm_Wc)
      CALL fm_contract_ABA(fm_M(1, 1), fm_Wc)
      CALL fm_contract_ABA(fm_V_msqrt, fm_Wc)
      CALL cp_fm_create(fm_eps_inv, bs_env%fm_RI_RI%matrix_struct)
      CALL cp_fm_set_all(fm_eps_inv, 0.0_dp, 1.0_dp)
      CALL cp_fm_scale_and_add(1.0_dp, fm_eps_inv, 1.0_dp, fm_Wc)

      CALL cp_fm_release(fm_V)
      CALL cp_fm_release(fm_M)
      CALL cp_fm_release(fm_V_sqrt)
      CALL cp_fm_release(fm_V_msqrt)
      CALL cp_fm_release(fm_Wc)
      CALL cp_fm_release(fm_work)

      CALL timestop(handle)

   END SUBROUTINE link_screening

! **************************************************************************************************
!> \brief The window columns of the MO coefficients, C(μ, p) for p in the window, as a 2-index
!>        tensor (AO atom | MO block) on the tensor group's grid, the MO blocking the caller's
!> \param bs_env ...
!> \param ispin ...
!> \param window ...
!> \param sizes_mo MO block sizes of the tensor, summing to the window
!> \param t_C the tensor, created here
! **************************************************************************************************
   SUBROUTINE link_mo_tensor(bs_env, ispin, window, sizes_mo, t_C)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: ispin
      TYPE(mo_window_type), INTENT(IN)                   :: window
      INTEGER, DIMENSION(:), INTENT(IN)                  :: sizes_mo
      TYPE(dbt_type), INTENT(OUT)                        :: t_C

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'link_mo_tensor'

      INTEGER :: ao0, handle, i_atom, i_blk, i_chunk, i_first, i_last, iproc, ithread, mo0, n_blk, &
         n_c, n_chunk, n_mu, n_window, nblk_per_thread, nthreads
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: blk_ind_1, blk_ind_2, dist_1, dist_2
      INTEGER, DIMENSION(2)                              :: pdims_2d
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: blk, c_win
      TYPE(dbt_pgrid_type)                               :: pgrid_2d

      CALL timeset(routineN, handle)

      n_window = window%last_mo - window%first_mo + 1
      n_chunk = SIZE(sizes_mo)
      CPASSERT(SUM(sizes_mo) == n_window)
      ! the window columns replicated, the source of the local blocks below
      ALLOCATE (c_win(bs_env%n_ao, n_window))
      CALL cp_fm_get_submatrix(bs_env%fm_mo_coeff_Gamma(ispin), c_win, start_row=1, &
                               start_col=window%first_mo, n_rows=bs_env%n_ao, n_cols=n_window)

      pdims_2d = 0
      CALL dbt_pgrid_create(bs_env%para_env_tensor, pdims_2d, pgrid_2d)
      CALL create_2c_tensor(t_C, dist_1, dist_2, pgrid_2d, bs_env%sizes_AO, sizes_mo, name="(AO | MO)")
      DEALLOCATE (dist_1, dist_2)
      CALL dbt_pgrid_destroy(pgrid_2d)

      ! reserve and fill the blocks this rank stores
      n_blk = 0
      DO i_atom = 1, bs_env%n_atom
         DO i_chunk = 1, n_chunk
            CALL dbt_get_stored_coordinates(t_C, [i_atom, i_chunk], iproc)
            IF (iproc == bs_env%para_env_tensor%mepos) n_blk = n_blk + 1
         END DO
      END DO
      ALLOCATE (blk_ind_1(n_blk), blk_ind_2(n_blk))
      n_blk = 0
      DO i_atom = 1, bs_env%n_atom
         DO i_chunk = 1, n_chunk
            CALL dbt_get_stored_coordinates(t_C, [i_atom, i_chunk], iproc)
            IF (iproc /= bs_env%para_env_tensor%mepos) CYCLE
            n_blk = n_blk + 1
            blk_ind_1(n_blk) = i_atom
            blk_ind_2(n_blk) = i_chunk
         END DO
      END DO
      ! DBM reserves inside one parallel region, every thread over its own slice of the list
!$OMP PARALLEL DEFAULT(NONE) SHARED(t_C, blk_ind_1, blk_ind_2, n_blk) &
!$OMP PRIVATE(nthreads, ithread, nblk_per_thread, i_first, i_last)
      nthreads = 1
      ithread = 0
!$    nthreads = omp_get_num_threads()
!$    ithread = omp_get_thread_num()
      nblk_per_thread = n_blk/nthreads + 1
      i_first = ithread*nblk_per_thread + 1
      i_last = MIN(i_first + nblk_per_thread - 1, n_blk)
      CALL dbt_reserve_blocks(t_C, blk_ind_1(i_first:i_last), blk_ind_2(i_first:i_last))
!$OMP END PARALLEL
      DO i_blk = 1, n_blk
         i_atom = blk_ind_1(i_blk)
         i_chunk = blk_ind_2(i_blk)
         ao0 = bs_env%i_ao_start_from_atom(i_atom)
         n_mu = bs_env%sizes_AO(i_atom)
         mo0 = SUM(sizes_mo(1:i_chunk - 1)) + 1
         n_c = sizes_mo(i_chunk)
         ALLOCATE (blk(n_mu, n_c))
         blk(:, :) = c_win(ao0:ao0 + n_mu - 1, mo0:mo0 + n_c - 1)
         CALL dbt_put_block(t_C, [i_atom, i_chunk], [n_mu, n_c], blk)
         DEALLOCATE (blk)
      END DO
      CALL dbt_filter(t_C, bs_env%eps_filter)
      DEALLOCATE (blk_ind_1, blk_ind_2, c_win)

      CALL timestop(handle)

   END SUBROUTINE link_mo_tensor

! **************************************************************************************************
!> \brief Block sizes of the window index p: chunks of at most mo_chunk levels inside the occupied
!>        and inside the virtual sector, so that no block straddles the two
!> \param n_o ...
!> \param n_v ...
!> \param sizes the chunk sizes, occupied chunks first, allocated here
!> \param n_chunk_o occupied chunks
! **************************************************************************************************
   SUBROUTINE link_sector_chunks(n_o, n_v, sizes, n_chunk_o)
      INTEGER, INTENT(IN)                                :: n_o, n_v
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: sizes
      INTEGER, INTENT(OUT)                               :: n_chunk_o

      INTEGER                                            :: n_chunk_v

      n_chunk_o = (n_o - 1)/mo_chunk + 1
      n_chunk_v = (n_v - 1)/mo_chunk + 1
      ALLOCATE (sizes(n_chunk_o + n_chunk_v))
      sizes(:) = mo_chunk
      sizes(n_chunk_o) = n_o - (n_chunk_o - 1)*mo_chunk
      sizes(n_chunk_o + n_chunk_v) = n_v - (n_chunk_v - 1)*mo_chunk

   END SUBROUTINE link_sector_chunks

! **************************************************************************************************
!> \brief Contiguous RI atom ranges of at most max_fn RI functions; an atom larger than that forms
!>        its own batch
!> \param sizes_RI RI functions per atom
!> \param max_fn RI functions per batch the memory budget allows
!> \param batch_first first atom per batch, allocated here
!> \param batch_last last atom per batch, allocated here
! **************************************************************************************************
   SUBROUTINE link_ri_batches(sizes_RI, max_fn, batch_first, batch_last)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: sizes_RI
      INTEGER, INTENT(IN)                                :: max_fn
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: batch_first, batch_last

      INTEGER                                            :: i_atom, n_atom, n_batch, n_in
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: first, last

      n_atom = SIZE(sizes_RI)
      ALLOCATE (first(n_atom), last(n_atom))
      n_batch = 0
      n_in = 0
      DO i_atom = 1, n_atom
         IF (n_in > 0 .AND. n_in + sizes_RI(i_atom) > max_fn) THEN
            last(n_batch) = i_atom - 1
            n_in = 0
         END IF
         IF (n_in == 0) THEN
            n_batch = n_batch + 1
            first(n_batch) = i_atom
         END IF
         n_in = n_in + sizes_RI(i_atom)
      END DO
      last(n_batch) = n_atom
      ALLOCATE (batch_first(n_batch), batch_last(n_batch))
      batch_first(:) = first(1:n_batch)
      batch_last(:) = last(1:n_batch)
      DEALLOCATE (first, last)

   END SUBROUTINE link_ri_batches

! **************************************************************************************************
!> \brief The memory budget of the slab build per rank and the largest RI batch it allows:
!>        BSE_ITERAT%MEMORY_BUDGET_GB as given, else mem_fraction of the free memory per rank measured
!>        now (its minimum over the ranks, so every rank cuts the same batches); MEMORY_CHECK decides
!>        what happens when one RI atom alone exceeds the budget
!> \param bs_env ...
!> \param bse_env MEMORY_BUDGET_GB and MEMORY_CHECK
!> \param n_w MOs in the window
!> \param unit_nr ...
!> \param budget_bytes the budget; the RI transform chunks its work matrix by it as well
!> \param max_fn RI functions per batch within the budget, at least one
! **************************************************************************************************
   SUBROUTINE link_batch_budget(bs_env, bse_env, n_w, unit_nr, budget_bytes, max_fn)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(bse_env_type), INTENT(IN)                     :: bse_env
      INTEGER, INTENT(IN)                                :: n_w, unit_nr
      INTEGER(KIND=int_8), INTENT(OUT)                   :: budget_bytes
      INTEGER, INTENT(OUT)                               :: max_fn

      CHARACTER(LEN=16)                                  :: bud_str, need_str
      INTEGER(KIND=int_8)                                :: bytes_per_fn, largest_bytes, &
                                                            one_atom_bytes
      REAL(KIND=dp)                                      :: budget_GB, mem_avail_GB

      budget_GB = bse_env%memory_budget_gb
      IF (budget_GB < 0.0_dp) THEN
         CALL mp_mem_avail_per_rank_GB(bs_env%para_env, mem_avail_GB)
         budget_GB = mem_fraction*mem_avail_GB
      END IF
      budget_bytes = INT(MAX(budget_GB, 0.0_dp)*1.0E9_dp, KIND=int_8)

      ! per RI function of a batch, all distributed: on the ranks of its tensor group the 3c
      ! integrals (P|μν), the half-transform (P μ|q) and the product (P|qp); on all ranks the block
      ! matrices of the batch, their global copies (two inside the group-to-global hop) and the
      ! batch's full matrix, four (P|qp); the contractions' own scratch is not counted
      bytes_per_fn = 8_int_8*((INT(bs_env%n_ao, int_8)**2 + INT(bs_env%n_ao, int_8)*INT(n_w, int_8) + &
                               INT(n_w, int_8)**2)/INT(bs_env%para_env_tensor%num_pe, int_8) + &
                              4_int_8*INT(n_w, int_8)**2/INT(bs_env%para_env%num_pe, int_8))
      max_fn = INT(MAX(1_int_8, MIN(INT(SUM(bs_env%sizes_RI), int_8), budget_bytes/bytes_per_fn)))
      one_atom_bytes = INT(MAXVAL(bs_env%sizes_RI), int_8)*bytes_per_fn
      largest_bytes = MAX(INT(max_fn, int_8)*bytes_per_fn, one_atom_bytes)

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A4,T7,A,T67,F14.3)') 'BSE|', 'Memory budget of the slab build per MPI rank (GB)', &
            budget_GB
         WRITE (unit_nr, '(T2,A4,T7,A,T67,F14.3)') 'BSE|', 'Largest RI batch per MPI rank (GB)', &
            REAL(largest_bytes, dp)*1.0E-9_dp
      END IF

      ! nothing smaller than one atom per batch exists, so CLAMP has nothing left to clamp
      IF (one_atom_bytes > budget_bytes .AND. bse_env%memory_check /= bse_memcheck_off) THEN
         WRITE (need_str, '(F16.3)') REAL(one_atom_bytes, dp)*1.0E-9_dp
         WRITE (bud_str, '(F16.3)') budget_GB
         IF (bse_env%memory_check == bse_memcheck_warn) THEN
            CALL cp_warn(__LOCATION__, &
                         "PROPERTIES%BSE: the slab batch of one RI atom needs "//TRIM(ADJUSTL(need_str))// &
                         " GB per MPI rank, the budget is "//TRIM(ADJUSTL(bud_str))//" GB; one atom per batch")
         ELSE
            CALL cp_abort(__LOCATION__, &
                          "PROPERTIES%BSE: not enough memory for the slab batch of one RI atom: "// &
                          TRIM(ADJUSTL(need_str))//" GB per MPI rank, the budget is "//TRIM(ADJUSTL(bud_str))// &
                          " GB. Raise BSE_ITERAT%MEMORY_BUDGET_GB or use more MPI ranks.")
         END IF
      END IF

   END SUBROUTINE link_batch_budget

! **************************************************************************************************
!> \brief The three window-sized slabs of one spin, S_P,pq = sum_Q (pq|Q) T_QP with
!>        (pq|Q) = sum_μν C_μp (Q|μν) C_νq, batched over RI atoms: the 3c integrals of the batch
!>        are rebuilt, half-transformed by a tensor contraction over ν, the second half runs on the
!>        blocks this rank holds, the batch is summed over all ranks and its rows scattered into
!>        the slabs; the RI transform by T follows in column chunks. Column conventions of the
!>        solver core: ia = (i-1) n_v + a, ij = (i-1) n_o + j, ab = (a-1) n_v + b
!> \param qs_env ...
!> \param bs_env ...
!> \param bse_env its memory budget and memory check
!> \param ispin ...
!> \param window ...
!> \param n_o occupied levels in the window
!> \param n_v virtual levels in the window
!> \param fm_T the RI transform, RI x RI
!> \param fm_S_ij created here, RI x n_o*n_o
!> \param fm_S_ia created here, RI x n_o*n_v
!> \param fm_S_ab created here, RI x n_v*n_v
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE link_slabs(qs_env, bs_env, bse_env, ispin, window, n_o, n_v, fm_T, fm_S_ij, fm_S_ia, &
                         fm_S_ab, unit_nr)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(bse_env_type), INTENT(IN)                     :: bse_env
      INTEGER, INTENT(IN)                                :: ispin
      TYPE(mo_window_type), INTENT(IN)                   :: window
      INTEGER, INTENT(IN)                                :: n_o, n_v
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_T
      TYPE(cp_fm_type), INTENT(OUT)                      :: fm_S_ij, fm_S_ia, fm_S_ab
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'link_slabs'

      INTEGER                                            :: a1, a2, handle, i_batch, max_fn, &
                                                            n_batch, n_chunk_o, n_RI, n_w
      INTEGER(KIND=int_8)                                :: batch_bytes
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: batch_first, batch_last, dist_1, dist_2, &
                                                            dist_3, sizes_p, sizes_q
      INTEGER, DIMENSION(3)                              :: pdims_3d
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: struct_ab, struct_ia, struct_ij
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist
      TYPE(dbt_pgrid_type)                               :: pgrid_3d
      TYPE(dbt_type)                                     :: t_3c, t_C_p, t_C_q, t_S, t_T

      CALL timeset(routineN, handle)

      NULLIFY (blacs_env, dbcsr_dist, struct_ab, struct_ia, struct_ij)
      CALL get_qs_env(qs_env, blacs_env=blacs_env, dbcsr_dist=dbcsr_dist)
      n_RI = bs_env%n_RI
      n_w = n_o + n_v
      CALL cp_fm_struct_create(struct_ij, para_env=bs_env%para_env, context=blacs_env, &
                               nrow_global=n_RI, ncol_global=n_o*n_o)
      CALL cp_fm_struct_create(struct_ia, para_env=bs_env%para_env, context=blacs_env, &
                               nrow_global=n_RI, ncol_global=n_o*n_v)
      CALL cp_fm_struct_create(struct_ab, para_env=bs_env%para_env, context=blacs_env, &
                               nrow_global=n_RI, ncol_global=n_v*n_v)
      CALL cp_fm_create(fm_S_ij, struct_ij, name="fm_S_ij_newgw", set_zero=.TRUE.)
      CALL cp_fm_create(fm_S_ia, struct_ia, name="fm_S_ia_newgw", set_zero=.TRUE.)
      CALL cp_fm_create(fm_S_ab, struct_ab, name="fm_S_ab_newgw", set_zero=.TRUE.)
      CALL cp_fm_struct_release(struct_ij)
      CALL cp_fm_struct_release(struct_ia)
      CALL cp_fm_struct_release(struct_ab)

      ! the window index q in two blocks, occupied and virtual, so that every (P, q, p) block of the
      ! product covers the whole q range of one pair type; p in chunks inside each sector
      ALLOCATE (sizes_q(2))
      sizes_q(:) = [n_o, n_v]
      CALL link_sector_chunks(n_o, n_v, sizes_p, n_chunk_o)
      CALL link_mo_tensor(bs_env, ispin, window, sizes_q, t_C_q)
      CALL link_mo_tensor(bs_env, ispin, window, sizes_p, t_C_p)

      CALL link_batch_budget(bs_env, bse_env, n_w, unit_nr, batch_bytes, max_fn)
      CALL link_ri_batches(bs_env%sizes_RI, max_fn, batch_first, batch_last)
      n_batch = SIZE(batch_first)
      IF (unit_nr > 0) WRITE (unit_nr, '(T2,A4,T7,A,T71,I10)') 'BSE|', 'RI atom batches of the slabs', n_batch

      ! (RI AO | MO) of the half-transform and (RI | MO MO) of the product on the tensor group's grid
      pdims_3d = 0
      CALL dbt_pgrid_create(bs_env%para_env_tensor, pdims_3d, pgrid_3d)
      CALL create_3c_tensor(t_T, dist_1, dist_2, dist_3, pgrid_3d, bs_env%sizes_RI, bs_env%sizes_AO, &
                            sizes_q, map1=[1, 2], map2=[3], name="(RI AO | MO)")
      DEALLOCATE (dist_1, dist_2, dist_3)
      CALL dbt_pgrid_destroy(pgrid_3d)
      pdims_3d = 0
      CALL dbt_pgrid_create(bs_env%para_env_tensor, pdims_3d, pgrid_3d, &
                            tensor_dims=[bs_env%n_atom, SIZE(sizes_q), SIZE(sizes_p)])
      CALL create_3c_tensor(t_S, dist_1, dist_2, dist_3, pgrid_3d, bs_env%sizes_RI, sizes_q, sizes_p, &
                            map1=[1], map2=[2, 3], name="(RI | MO MO)")
      DEALLOCATE (dist_1, dist_2, dist_3)
      CALL dbt_pgrid_destroy(pgrid_3d)
      CALL dbt_create(bs_env%t_RI_AO__AO, t_3c)

      ! C_q contracts once per batch in a layout-compatible contraction: the batched context keeps its replicated
      ! copy across the batches; on the small operand only, so that no grid changes under the neighbour list
      ! of the 3c integrals. No context on C_p: no layout of T fits both contractions, so the second remaps T anyway
      CALL dbt_batched_contract_init(t_C_q)

      DO i_batch = 1, n_batch
         a1 = batch_first(i_batch)
         a2 = batch_last(i_batch)
         ! one tensor group per batch, round robin; every rank joins the hop into the slabs
         IF (MODULO(i_batch - 1, bs_env%num_tensor_groups) == bs_env%tensor_group_color) THEN
            CALL compute_3c_integrals(qs_env, bs_env, t_3c, atoms_RI=[a1, a2])
            ! T_Pμq = sum_ν (P|μν) C_νq
            CALL dbt_contract(alpha=1.0_dp, tensor_1=t_3c, tensor_2=t_C_q, beta=0.0_dp, tensor_3=t_T, &
                              contract_1=[3], notcontract_1=[1, 2], map_1=[1, 2], &
                              contract_2=[1], notcontract_2=[2], map_2=[3], &
                              filter_eps=bs_env%eps_filter)
            CALL dbt_clear(t_3c)
            ! S_Pqp = sum_μ T_Pμq C_μp, distributed over the group
            CALL dbt_contract(alpha=1.0_dp, tensor_1=t_T, tensor_2=t_C_p, beta=0.0_dp, tensor_3=t_S, &
                              contract_1=[2], notcontract_1=[1, 3], map_1=[1, 2], &
                              contract_2=[1], notcontract_2=[2], map_2=[3], &
                              filter_eps=bs_env%eps_filter)
            CALL dbt_clear(t_T)
         END IF
         CALL link_slab_rows(bs_env, blacs_env, dbcsr_dist, t_S, n_o, n_v, sizes_p, n_chunk_o, a1, a2, &
                             fm_S_ij, fm_S_ia, fm_S_ab)
         CALL dbt_clear(t_S)
      END DO

      CALL dbt_batched_contract_finalize(t_C_q)
      CALL dbt_destroy(t_3c)
      CALL dbt_destroy(t_T)
      CALL dbt_destroy(t_S)
      CALL dbt_destroy(t_C_q)
      CALL dbt_destroy(t_C_p)
      DEALLOCATE (sizes_q, sizes_p, batch_first, batch_last)

      ! S_P,pq = sum_Q T_QP (Q|pq), in column chunks, in place
      CALL apply_ri_transform(fm_T, fm_S_ij, n_RI, n_o*n_o, batch_bytes)
      CALL apply_ri_transform(fm_T, fm_S_ia, n_RI, n_o*n_v, batch_bytes)
      CALL apply_ri_transform(fm_T, fm_S_ab, n_RI, n_v*n_v, batch_bytes)

      CALL timestop(handle)

   END SUBROUTINE link_slabs

! **************************************************************************************************
!> \brief The batch product S(P, q, p) into rows P0.. of the three slabs: the local (P, q, p) blocks
!>        into one block matrix per pair type on the tensor group's grid, the group-to-global hop of
!>        the GW, the block-cyclic copy, the row copy into the slab; every rank takes part, a rank
!>        outside the batch's group with no blocks
!> \param bs_env ...
!> \param blacs_env context of the slabs
!> \param dbcsr_dist a distribution over all ranks, the template of the global block matrices
!> \param t_S S(P, q, p) of the RI atoms a1..a2, q in two blocks (occupied, virtual), p in sector chunks
!> \param n_o ...
!> \param n_v ...
!> \param sizes_p the p chunks, occupied chunks first
!> \param n_chunk_o occupied chunks of p
!> \param a1 first RI atom of the batch
!> \param a2 last RI atom of the batch
!> \param fm_S_ij ...
!> \param fm_S_ia ...
!> \param fm_S_ab ...
! **************************************************************************************************
   SUBROUTINE link_slab_rows(bs_env, blacs_env, dbcsr_dist, t_S, n_o, n_v, sizes_p, n_chunk_o, a1, a2, &
                             fm_S_ij, fm_S_ia, fm_S_ab)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist
      TYPE(dbt_type), INTENT(INOUT)                      :: t_S
      INTEGER, INTENT(IN)                                :: n_o, n_v
      INTEGER, DIMENSION(:), INTENT(IN)                  :: sizes_p
      INTEGER, INTENT(IN)                                :: n_chunk_o, a1, a2
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_S_ij, fm_S_ia, fm_S_ab

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'link_slab_rows'
      INTEGER, PARAMETER                                 :: kind_ab = 3, kind_ia = 2, kind_ij = 1

      INTEGER                                            :: bp, bq, group_handle, handle, ib, ikind, &
                                                            iproc, jb, my_rank, n_chunk, nb_P, nP, &
                                                            npcol, npcol_g, nprow, nprow_g, P0, &
                                                            rank
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: coords
      INTEGER, DIMENSION(3)                              :: blk_size, bp0_of, bq_of, ind, n_cols, &
                                                            nb_col
      INTEGER, DIMENSION(:), POINTER                     :: col_dist, col_dist_g, col_sizes, &
                                                            row_dist, row_dist_g, row_sizes
      INTEGER, DIMENSION(:, :), POINTER                  :: pgrid
      LOGICAL                                            :: found
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: blk
      TYPE(cp_fm_struct_type), POINTER                   :: struct
      TYPE(cp_fm_type)                                   :: fm_batch
      TYPE(dbcsr_distribution_type), DIMENSION(3)        :: dist_glob, dist_loc
      TYPE(dbcsr_type), DIMENSION(3)                     :: mat_glob, mat_loc
      TYPE(dbt_iterator_type)                            :: iter

      CALL timeset(routineN, handle)

      NULLIFY (col_dist, col_dist_g, col_sizes, pgrid, row_dist, row_dist_g, row_sizes, struct)
      P0 = bs_env%i_RI_start_from_atom(a1)
      nP = bs_env%i_RI_end_from_atom(a2) - P0 + 1
      nb_P = a2 - a1 + 1
      n_chunk = SIZE(sizes_p)
      n_cols = [n_o*n_o, n_o*n_v, n_v*n_v]
      nb_col = [n_chunk_o, n_chunk_o, n_chunk - n_chunk_o]
      bq_of = [1, 2, 2]
      ! column block jb of a pair type is the p chunk bp = bp0_of + jb
      bp0_of = [0, 0, n_chunk_o]

      ! the group grid of the tensor as the process grid of the local matrices, its ranks placed by
      ! their coordinates, so that every block the tensor holds is local to its matrix
      ASSOCIATE (cart => t_S%pgrid%mp_comm_2d)
         my_rank = cart%mepos
         group_handle = cart%get_handle()
         nprow = cart%num_pe_cart(1)
         npcol = cart%num_pe_cart(2)
         ALLOCATE (coords(2, 0:cart%num_pe - 1))
         CALL cart%allgather(cart%mepos_cart, coords)
         ALLOCATE (pgrid(0:nprow - 1, 0:npcol - 1))
         DO rank = 0, cart%num_pe - 1
            pgrid(coords(1, rank), coords(2, rank)) = rank
         END DO
      END ASSOCIATE
      ALLOCATE (row_dist(nb_P), row_sizes(nb_P), row_dist_g(nb_P))
      CALL dbcsr_distribution_get(dbcsr_dist, nprows=nprow_g, npcols=npcol_g)
      DO ib = 1, nb_P
         CALL dbt_get_stored_coordinates(t_S, [a1 + ib - 1, 1, 1], iproc)
         row_dist(ib) = coords(1, iproc)
         row_sizes(ib) = bs_env%sizes_RI(a1 + ib - 1)
         row_dist_g(ib) = MOD(ib - 1, nprow_g)
      END DO
      DO ikind = 1, 3
         ALLOCATE (col_dist(nb_col(ikind)), col_sizes(nb_col(ikind)), col_dist_g(nb_col(ikind)))
         DO jb = 1, nb_col(ikind)
            bp = bp0_of(ikind) + jb
            CALL dbt_get_stored_coordinates(t_S, [a1, bq_of(ikind), bp], iproc)
            col_dist(jb) = coords(2, iproc)
            ! the q range of the pair type times the p chunk
            col_sizes(jb) = MERGE(n_o, n_v, ikind == kind_ij)*sizes_p(bp)
            col_dist_g(jb) = MOD(jb - 1, npcol_g)
         END DO
         CALL dbcsr_distribution_new(dist_loc(ikind), group=group_handle, pgrid=pgrid, &
                                     row_dist=row_dist, col_dist=col_dist)
         CALL dbcsr_create(mat_loc(ikind), name="S batch, group", dist=dist_loc(ikind), &
                           matrix_type=dbcsr_type_no_symmetry, row_blk_size=row_sizes, col_blk_size=col_sizes)
         ! the threaded put_block below cannot create the work matrix itself
         CALL dbcsr_work_create(mat_loc(ikind), work_mutable=.TRUE.)
         ! the global matrix of the batch, round robin over the process grid of the slabs
         CALL dbcsr_distribution_new(dist_glob(ikind), template=dbcsr_dist, row_dist=row_dist_g, &
                                     col_dist=col_dist_g)
         CALL dbcsr_create(mat_glob(ikind), name="S batch", dist=dist_glob(ikind), &
                           matrix_type=dbcsr_type_no_symmetry, row_blk_size=row_sizes, col_blk_size=col_sizes)
         DEALLOCATE (col_dist, col_sizes, col_dist_g)
      END DO
      DEALLOCATE (pgrid, coords, row_dist, row_sizes, row_dist_g)

      ! the local (P, q, p) blocks as (P, qp) blocks, q fastest, the column order of the slabs; the
      ! iterator needs the full thread team
!$OMP PARALLEL DEFAULT(NONE) SHARED(t_S, mat_loc, n_chunk_o, a1, bp0_of, my_rank) &
!$OMP PRIVATE(iter, ind, blk_size, blk, found, bq, bp, ikind, ib, jb, iproc)
      CALL dbt_iterator_start(iter, t_S)
      DO WHILE (dbt_iterator_blocks_left(iter))
         CALL dbt_iterator_next_block(iter, ind, blk_size=blk_size)
         bq = ind(2)
         bp = ind(3)
         ! pair types of the (q block, p chunk) blocks: (occ, occ) ij, (virt, occ) ia, (virt, virt) ab;
         ! (occ, virt) is the ai block, which no slab holds
         IF (bq == 1 .AND. bp > n_chunk_o) CYCLE
         IF (bq == 1) THEN
            ikind = kind_ij
         ELSE IF (bp <= n_chunk_o) THEN
            ikind = kind_ia
         ELSE
            ikind = kind_ab
         END IF
         ALLOCATE (blk(blk_size(1), blk_size(2), blk_size(3)))
         CALL dbt_get_block(t_S, ind, blk_size, blk, found)
         IF (found) THEN
            ib = ind(1) - a1 + 1
            jb = bp - bp0_of(ikind)
            CALL dbcsr_get_stored_coordinates(mat_loc(ikind), ib, jb, iproc)
            CPASSERT(iproc == my_rank)
            CALL dbcsr_put_block(mat_loc(ikind), ib, jb, RESHAPE(blk, [blk_size(1), blk_size(2)*blk_size(3)]))
         END IF
         DEALLOCATE (blk)
      END DO
      CALL dbt_iterator_stop(iter)
!$OMP END PARALLEL

      DO ikind = 1, 3
         CALL dbcsr_finalize(mat_loc(ikind))
         CALL local_matrix_to_global_matrix(mat_loc(ikind), mat_glob(ikind), bs_env%para_env)
         CALL cp_fm_struct_create(struct, para_env=bs_env%para_env, context=blacs_env, &
                                  nrow_global=nP, ncol_global=n_cols(ikind))
         CALL cp_fm_create(fm_batch, struct)
         CALL cp_fm_struct_release(struct)
         CALL copy_dbcsr_to_fm(mat_glob(ikind), fm_batch)
         SELECT CASE (ikind)
         CASE (kind_ij)
            CALL cp_fm_to_fm_submat(fm_batch, fm_S_ij, nP, n_cols(ikind), 1, 1, P0, 1)
         CASE (kind_ia)
            CALL cp_fm_to_fm_submat(fm_batch, fm_S_ia, nP, n_cols(ikind), 1, 1, P0, 1)
         CASE (kind_ab)
            CALL cp_fm_to_fm_submat(fm_batch, fm_S_ab, nP, n_cols(ikind), 1, 1, P0, 1)
         END SELECT
         CALL cp_fm_release(fm_batch)
         CALL dbcsr_release(mat_loc(ikind))
         CALL dbcsr_release(mat_glob(ikind))
         CALL dbcsr_distribution_release(dist_loc(ikind))
         CALL dbcsr_distribution_release(dist_glob(ikind))
      END DO

      CALL timestop(handle)

   END SUBROUTINE link_slab_rows

! **************************************************************************************************
!> \brief slab_P,c <- sum_Q T_QP slab_Q,c over the columns c of one chunk, with one RI x n_chunk
!>        work matrix
!> \param fm_T ...
!> \param fm_slab ...
!> \param n_RI ...
!> \param n_cols ...
!> \param batch_bytes per-rank budget of the work matrix, which is distributed over the ranks
! **************************************************************************************************
   SUBROUTINE apply_ri_transform(fm_T, fm_slab, n_RI, n_cols, batch_bytes)
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_T, fm_slab
      INTEGER, INTENT(IN)                                :: n_RI, n_cols
      INTEGER(KIND=int_8), INTENT(IN)                    :: batch_bytes

      CHARACTER(LEN=*), PARAMETER :: routineN = 'apply_ri_transform'

      INTEGER                                            :: c0, handle, n_chunk, nc
      TYPE(cp_fm_struct_type), POINTER                   :: struct_work
      TYPE(cp_fm_type)                                   :: fm_work

      CALL timeset(routineN, handle)

      n_chunk = INT(MAX(1_int_8, MIN(INT(n_cols, int_8), &
                                     batch_bytes*INT(fm_slab%matrix_struct%para_env%num_pe, int_8) &
                                     /(8_int_8*INT(n_RI, int_8)))))
      NULLIFY (struct_work)
      CALL cp_fm_struct_create(struct_work, para_env=fm_slab%matrix_struct%para_env, &
                               context=fm_slab%matrix_struct%context, nrow_global=n_RI, ncol_global=n_chunk)
      CALL cp_fm_create(fm_work, struct_work)
      CALL cp_fm_struct_release(struct_work)
      DO c0 = 1, n_cols, n_chunk
         nc = MIN(n_chunk, n_cols - c0 + 1)
         CALL parallel_gemm("T", "N", n_RI, nc, n_RI, 1.0_dp, fm_T, fm_slab, 0.0_dp, fm_work, &
                            b_first_col=c0, c_first_col=1)
         CALL cp_fm_to_fm_submat(fm_work, fm_slab, nrow=n_RI, ncol=nc, s_firstrow=1, s_firstcol=1, &
                                 t_firstrow=1, t_firstcol=c0)
      END DO
      CALL cp_fm_release(fm_work)

      CALL timestop(handle)

   END SUBROUTINE apply_ri_transform

END MODULE bse_bs_env_link
