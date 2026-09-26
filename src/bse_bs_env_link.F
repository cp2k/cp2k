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
   USE bse_types,                       ONLY: bse_env_release,&
                                              bse_env_type
   USE bse_util,                        ONLY: bse_window_from_dft
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add
   USE cp_fm_diag,                      ONLY: cp_fm_power
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_submatrix,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_set_submatrix,&
                                              cp_fm_to_fm,&
                                              cp_fm_to_fm_submat,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE dbt_api,                         ONLY: &
        dbt_clear, dbt_contract, dbt_create, dbt_destroy, dbt_filter, dbt_get_block, &
        dbt_get_stored_coordinates, dbt_iterator_blocks_left, dbt_iterator_next_block, &
        dbt_iterator_start, dbt_iterator_stop, dbt_iterator_type, dbt_pgrid_create, &
        dbt_pgrid_destroy, dbt_pgrid_type, dbt_put_block, dbt_reserve_blocks, dbt_type
   USE gw_tensor_large_cell_Gamma,      ONLY: compute_3c_integrals
   USE gw_utils_fm,                     ONLY: fm_contract_ABA
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE libint_2c_3c,                    ONLY: libint_potential_type
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
   ! bytes of one RI batch of (P, p, q) transformed integrals and of one column chunk of a slab
   INTEGER(KIND=int_8), PARAMETER, PRIVATE :: slab_batch_bytes = 2_int_8**28
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
         CALL link_slabs(qs_env, bs_env, ispin, window, homo_red(ispin), virt_red(ispin), fm_T, &
                         fm_S_ij(ispin), fm_S_ia(ispin), fm_S_ab(ispin))
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
!> \brief The window columns of the MO coefficients, C(μ, p) for p in the window, as a replicated
!>        array and as a 2-index tensor (AO atom | MO chunk) on the tensor group's grid
!> \param bs_env ...
!> \param ispin ...
!> \param window ...
!> \param sizes_mo MO block sizes of the tensor, allocated here
!> \param c_win C(n_ao, n_window), allocated here
!> \param t_C the tensor, created here
! **************************************************************************************************
   SUBROUTINE link_mo_tensor(bs_env, ispin, window, sizes_mo, c_win, t_C)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: ispin
      TYPE(mo_window_type), INTENT(IN)                   :: window
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: sizes_mo
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: c_win
      TYPE(dbt_type), INTENT(OUT)                        :: t_C

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'link_mo_tensor'

      INTEGER :: ao0, handle, i_atom, i_blk, i_chunk, i_first, i_last, iproc, ithread, mo0, n_blk, &
         n_c, n_chunk, n_mu, n_window, nblk_per_thread, nthreads
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: blk_ind_1, blk_ind_2, dist_1, dist_2
      INTEGER, DIMENSION(2)                              :: pdims_2d
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: blk
      TYPE(dbt_pgrid_type)                               :: pgrid_2d

      CALL timeset(routineN, handle)

      n_window = window%last_mo - window%first_mo + 1
      ALLOCATE (c_win(bs_env%n_ao, n_window))
      CALL cp_fm_get_submatrix(bs_env%fm_mo_coeff_Gamma(ispin), c_win, start_row=1, &
                               start_col=window%first_mo, n_rows=bs_env%n_ao, n_cols=n_window)

      ! MO chunks of at most mo_chunk columns
      n_chunk = (n_window - 1)/mo_chunk + 1
      ALLOCATE (sizes_mo(n_chunk))
      sizes_mo(:) = mo_chunk
      sizes_mo(n_chunk) = n_window - (n_chunk - 1)*mo_chunk

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
         mo0 = (i_chunk - 1)*mo_chunk + 1
         n_c = sizes_mo(i_chunk)
         ALLOCATE (blk(n_mu, n_c))
         blk(:, :) = c_win(ao0:ao0 + n_mu - 1, mo0:mo0 + n_c - 1)
         CALL dbt_put_block(t_C, [i_atom, i_chunk], [n_mu, n_c], blk)
         DEALLOCATE (blk)
      END DO
      CALL dbt_filter(t_C, bs_env%eps_filter)
      DEALLOCATE (blk_ind_1, blk_ind_2)

      CALL timestop(handle)

   END SUBROUTINE link_mo_tensor

! **************************************************************************************************
!> \brief Contiguous RI atom ranges whose (P, p, q) batch stays under slab_batch_bytes; an atom
!>        larger than the budget forms its own batch
!> \param sizes_RI RI functions per atom
!> \param n_w MOs in the window
!> \param batch_first first atom per batch, allocated here
!> \param batch_last last atom per batch, allocated here
! **************************************************************************************************
   SUBROUTINE link_ri_batches(sizes_RI, n_w, batch_first, batch_last)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: sizes_RI
      INTEGER, INTENT(IN)                                :: n_w
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: batch_first, batch_last

      INTEGER                                            :: i_atom, max_fn, n_atom, n_batch, n_in
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: first, last

      n_atom = SIZE(sizes_RI)
      max_fn = MAX(1, INT(slab_batch_bytes/(8_int_8*INT(n_w, int_8)*INT(n_w, int_8))))
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
!> \brief The three window-sized slabs of one spin, S_P,pq = sum_Q (pq|Q) T_QP with
!>        (pq|Q) = sum_μν C_μp (Q|μν) C_νq, batched over RI atoms: the 3c integrals of the batch
!>        are rebuilt, half-transformed by a tensor contraction over ν, the second half runs on the
!>        blocks this rank holds, the batch is summed over all ranks and its rows scattered into
!>        the slabs; the RI transform by T follows in column chunks. Column conventions of the
!>        solver core: ia = (i-1) n_v + a, ij = (i-1) n_o + j, ab = (a-1) n_v + b
!> \param qs_env ...
!> \param bs_env ...
!> \param ispin ...
!> \param window ...
!> \param n_o occupied levels in the window
!> \param n_v virtual levels in the window
!> \param fm_T the RI transform, RI x RI
!> \param fm_S_ij created here, RI x n_o*n_o
!> \param fm_S_ia created here, RI x n_o*n_v
!> \param fm_S_ab created here, RI x n_v*n_v
! **************************************************************************************************
   SUBROUTINE link_slabs(qs_env, bs_env, ispin, window, n_o, n_v, fm_T, fm_S_ij, fm_S_ia, fm_S_ab)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: ispin
      TYPE(mo_window_type), INTENT(IN)                   :: window
      INTEGER, INTENT(IN)                                :: n_o, n_v
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_T
      TYPE(cp_fm_type), INTENT(OUT)                      :: fm_S_ij, fm_S_ia, fm_S_ab

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'link_slabs'

      INTEGER                                            :: a, a1, a2, b, handle, i, i_batch, j, &
                                                            mu0, n_batch, n_RI, n_w, nP, p, P0, &
                                                            p_off, q0
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: batch_first, batch_last, dist_1, dist_2, &
                                                            dist_3, sizes_mo
      INTEGER, DIMENSION(3)                              :: blk_size, ind, pdims_3d
      LOGICAL                                            :: found
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: c_win, cols, tmp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: blk, s_loc
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: struct_ab, struct_ia, struct_ij
      TYPE(dbt_iterator_type)                            :: iter
      TYPE(dbt_pgrid_type)                               :: pgrid_3d
      TYPE(dbt_type)                                     :: t_3c, t_C, t_T

      CALL timeset(routineN, handle)

      NULLIFY (blacs_env, struct_ab, struct_ia, struct_ij)
      CALL get_qs_env(qs_env, blacs_env=blacs_env)
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

      CALL link_mo_tensor(bs_env, ispin, window, sizes_mo, c_win, t_C)

      CALL link_ri_batches(bs_env%sizes_RI, n_w, batch_first, batch_last)
      n_batch = SIZE(batch_first)

      ! (RI AO | MO) result of the half-transform on the tensor group's grid
      pdims_3d = 0
      CALL dbt_pgrid_create(bs_env%para_env_tensor, pdims_3d, pgrid_3d)
      CALL create_3c_tensor(t_T, dist_1, dist_2, dist_3, pgrid_3d, bs_env%sizes_RI, bs_env%sizes_AO, &
                            sizes_mo, map1=[1, 2], map2=[3], name="(RI AO | MO)")
      DEALLOCATE (dist_1, dist_2, dist_3)
      CALL dbt_pgrid_destroy(pgrid_3d)
      CALL dbt_create(bs_env%t_RI_AO__AO, t_3c)

      DO i_batch = 1, n_batch
         a1 = batch_first(i_batch)
         a2 = batch_last(i_batch)
         P0 = bs_env%i_RI_start_from_atom(a1)
         nP = bs_env%i_RI_end_from_atom(a2) - P0 + 1
         ! s_loc(p, q, P): the (p, q) slice of one RI function is contiguous, the DGEMM operand below
         ALLOCATE (s_loc(n_w, n_w, nP))
         s_loc(:, :, :) = 0.0_dp

         ! one tensor group per batch, round robin; every rank joins the sum below
         IF (MODULO(i_batch - 1, bs_env%num_tensor_groups) == bs_env%tensor_group_color) THEN
            CALL compute_3c_integrals(qs_env, bs_env, t_3c, atoms_RI=[a1, a2])
            ! T_Pμq = sum_ν (P|μν) C_νq
            CALL dbt_contract(alpha=1.0_dp, tensor_1=t_3c, tensor_2=t_C, beta=0.0_dp, tensor_3=t_T, &
                              contract_1=[3], notcontract_1=[1, 2], map_1=[1, 2], &
                              contract_2=[1], notcontract_2=[2], map_2=[3], &
                              filter_eps=bs_env%eps_filter)
            CALL dbt_clear(t_3c)
            ! s_loc(p, q, P) += sum_μ C_μp T_Pμq over the blocks this rank holds; the iterator
            ! needs the full thread team, the accumulation into one (p, q, P) slice is serialised
!$OMP PARALLEL DEFAULT(NONE) SHARED(t_T, bs_env, c_win, s_loc, P0, n_w) &
!$OMP PRIVATE(iter, ind, blk_size, blk, found, p_off, mu0, q0, tmp, p)
            CALL dbt_iterator_start(iter, t_T)
            DO WHILE (dbt_iterator_blocks_left(iter))
               CALL dbt_iterator_next_block(iter, ind, blk_size=blk_size)
               ALLOCATE (blk(blk_size(1), blk_size(2), blk_size(3)))
               CALL dbt_get_block(t_T, ind, blk_size, blk, found)
               IF (found) THEN
                  p_off = bs_env%i_RI_start_from_atom(ind(1)) - P0
                  mu0 = bs_env%i_ao_start_from_atom(ind(2))
                  q0 = (ind(3) - 1)*mo_chunk + 1
                  ALLOCATE (tmp(blk_size(2), blk_size(3)))
                  DO p = 1, blk_size(1)
                     tmp(:, :) = blk(p, :, :)
                     ! A = C(μ0.., 1:n_w) with leading dimension n_ao, transposed: n_w x n_μ
                     ! B = tmp, n_μ x n_c; C = s_loc(1:n_w, q0.., P) with leading dimension n_w
!$OMP CRITICAL(bse_link_slabs_sum)
                     CALL dgemm('T', 'N', n_w, blk_size(3), blk_size(2), 1.0_dp, &
                                c_win(mu0, 1), bs_env%n_ao, tmp, blk_size(2), 1.0_dp, &
                                s_loc(1, q0, p_off + p), n_w)
!$OMP END CRITICAL(bse_link_slabs_sum)
                  END DO
                  DEALLOCATE (tmp)
               END IF
               DEALLOCATE (blk)
            END DO
            CALL dbt_iterator_stop(iter)
!$OMP END PARALLEL
            CALL dbt_clear(t_T)
         END IF

         CALL bs_env%para_env%sum(s_loc)

         ! the batch rows into the three slabs (every rank writes the elements it owns)
         ALLOCATE (cols(nP, n_o*n_o))
         DO i = 1, n_o
            DO j = 1, n_o
               cols(:, (i - 1)*n_o + j) = s_loc(i, j, :)
            END DO
         END DO
         CALL cp_fm_set_submatrix(fm_S_ij, cols, start_row=P0, start_col=1, n_rows=nP, n_cols=n_o*n_o)
         DEALLOCATE (cols)
         ALLOCATE (cols(nP, n_o*n_v))
         DO i = 1, n_o
            DO a = 1, n_v
               cols(:, (i - 1)*n_v + a) = s_loc(i, n_o + a, :)
            END DO
         END DO
         CALL cp_fm_set_submatrix(fm_S_ia, cols, start_row=P0, start_col=1, n_rows=nP, n_cols=n_o*n_v)
         DEALLOCATE (cols)
         ALLOCATE (cols(nP, n_v*n_v))
         DO a = 1, n_v
            DO b = 1, n_v
               cols(:, (a - 1)*n_v + b) = s_loc(n_o + a, n_o + b, :)
            END DO
         END DO
         CALL cp_fm_set_submatrix(fm_S_ab, cols, start_row=P0, start_col=1, n_rows=nP, n_cols=n_v*n_v)
         DEALLOCATE (cols, s_loc)
      END DO

      CALL dbt_destroy(t_3c)
      CALL dbt_destroy(t_T)
      CALL dbt_destroy(t_C)
      DEALLOCATE (c_win, sizes_mo, batch_first, batch_last)

      ! S_P,pq = sum_Q T_QP (Q|pq), in column chunks, in place
      CALL apply_ri_transform(fm_T, fm_S_ij, n_RI, n_o*n_o)
      CALL apply_ri_transform(fm_T, fm_S_ia, n_RI, n_o*n_v)
      CALL apply_ri_transform(fm_T, fm_S_ab, n_RI, n_v*n_v)

      CALL timestop(handle)

   END SUBROUTINE link_slabs

! **************************************************************************************************
!> \brief slab_P,c <- sum_Q T_QP slab_Q,c over the columns c of one chunk, with one RI x n_chunk
!>        work matrix
!> \param fm_T ...
!> \param fm_slab ...
!> \param n_RI ...
!> \param n_cols ...
! **************************************************************************************************
   SUBROUTINE apply_ri_transform(fm_T, fm_slab, n_RI, n_cols)
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_T, fm_slab
      INTEGER, INTENT(IN)                                :: n_RI, n_cols

      CHARACTER(LEN=*), PARAMETER :: routineN = 'apply_ri_transform'

      INTEGER                                            :: c0, handle, n_chunk, nc
      TYPE(cp_fm_struct_type), POINTER                   :: struct_work
      TYPE(cp_fm_type)                                   :: fm_work

      CALL timeset(routineN, handle)

      n_chunk = MIN(n_cols, MAX(1, INT(slab_batch_bytes/(8_int_8*INT(n_RI, int_8)))))
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
