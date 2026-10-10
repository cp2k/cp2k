!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Bethe-Salpeter equation on the RI-AO GW, the GW of PROPERTIES%BANDSTRUCTURE%GW with explicit
!>        three-center RI integrals over atomic orbitals (no RI_RS): reads PROPERTIES%BSE, builds the
!>        window-sized RI slabs and the effective inverse dielectric matrix from bs_env, and hands
!>        them to the solver core through prepare_bse_env_from_bandstructure_gw.
!>        Non-periodic systems only. Reads no real-time propagation input.
!> \par History
!>      09.2026 created [Maximilian Graml]
!> \author Maximilian Graml
! **************************************************************************************************
MODULE bse_bs_env_link
   USE bse_input,                       ONLY: read_bse_section
   USE bse_main,                        ONLY: bse_solve,&
                                              prepare_bse_env_from_bandstructure_gw
   USE bse_matvec,                      ONLY: bse_memory_left_GB,&
                                              bse_resolve_memory_per_proc
   USE bse_types,                       ONLY: bse_env_release,&
                                              bse_env_type
   USE bse_util,                        ONLY: bse_window_from_dft
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_create, dbcsr_distribution_get, dbcsr_distribution_new, dbcsr_distribution_release, &
        dbcsr_distribution_type, dbcsr_finalize, dbcsr_get_info, dbcsr_get_stored_coordinates, &
        dbcsr_put_block, dbcsr_release, dbcsr_type, dbcsr_type_no_symmetry, dbcsr_work_create
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
   USE cp_dbcsr_redistribution,         ONLY: cp_dbcsr_from_subgroups
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale_and_add
   USE cp_fm_diag,                      ONLY: cp_fm_power
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
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
        dbt_create, dbt_destroy, dbt_get_block, dbt_get_stored_coordinates, &
        dbt_iterator_blocks_left, dbt_iterator_next_block, dbt_iterator_start, dbt_iterator_stop, &
        dbt_iterator_type, dbt_pgrid_create, dbt_pgrid_destroy, dbt_pgrid_type, dbt_type
   USE dbt_matrix_redistribution,       ONLY: cp_fm_to_dbt
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
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'bse_bs_env_link'

   ! block size of the MO index p in the tensors: the AO and RI indices come in atomic blocks, the MO
   ! index has none. A block is the unit the tensor grid distributes and filters, so the size bounds
   ! the largest block of the MO product, (RI atom) x (q sector) x mo_block_size
   INTEGER, PARAMETER, PRIVATE :: mo_block_size = 64

   PUBLIC :: bse_from_bs_env

CONTAINS

! **************************************************************************************************
!> \brief Entry from the properties driver after the bandstructure block. Hands the solver core what
!>        its kernels are built from, for p, q in the active window:
!>        B^P_pq = sum_QR (pq|Q) [M^-1]_QR [V^1/2]_RP,
!>        ε^-1_eff,PQ = δ_PQ + sum_RS [V^-1/2]_PR W^c_RS(iω=0) [V^-1/2]_SQ,
!>        and the quasiparticle energies ε_p of the GW, so that, with the core's
!>        \bar{B}^P_ij = sum_Q ε^-1_eff,PQ B^Q_ij,
!>        A_ia,jb = δ_ij δ_ab (ε_a-ε_i) + α sum_P B^P_ia B^P_jb - w sum_P \bar{B}^P_ij B^P_ab
!>        B_ia,jb =                       α sum_P B^P_ia B^P_jb - w sum_P \bar{B}^P_ib B^P_ja
!>        V_PQ: bare Coulomb and M_PQ: RI metric in the RI basis, W^c = W - V: the correlation part
!>        of the GW's static screened interaction. Order: checks, section, window, screening,
!>        per-spin slabs, seam, solve, release
!> \param qs_env source of bs_env; the slab build reads its blacs context and matrix distribution
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
      TYPE(cp_fm_type)                                   :: fm_eps_inv, fm_M_inv_V_sqrt
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)        :: fm_S_ab, fm_S_ia, fm_S_ij
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(mo_window_type)                               :: window
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      unit_nr = cp_logger_get_default_io_unit(logger)
      NULLIFY (bs_env)
      CALL get_qs_env(qs_env, bs_env=bs_env)
      ! the driver calls this after the bandstructure block, which creates bs_env
      IF (.NOT. ASSOCIATED(bs_env)) CPABORT("PROPERTIES%BSE found no GW environment (bs_env)")
      ! check_cp2k_input refused each of these at input parsing with a message that says what to change;
      ! Gamma is the one k-point that carries the MOs and the energies read below
      IF (.NOT. bs_env%do_gw .OR. bs_env%do_periodic .OR. bs_env%do_gw_ri_rs .OR. bs_env%nkp_bs_and_DOS /= 1) THEN
         CPABORT("PROPERTIES%BSE needs a non-periodic, Gamma-only GW without RI_RS")
      END IF
      ! read_bse_input told the GW at its init to keep its static W and the 3c neighbour list
      IF (.NOT. (bs_env%keep_W_freq_zero .AND. bs_env%keep_3c_nl) .OR. &
          .NOT. ASSOCIATED(bs_env%fm_W_MIC_freq_zero%matrix_struct)) THEN
         CPABORT("PROPERTIES%BSE: the GW kept no static W(0) or no 3c neighbour list")
      END IF
      ! every RI batch of the slabs goes to one tensor group; unset groups would skip them all silently
      IF (bs_env%num_tensor_groups < 1 .OR. bs_env%tensor_group_color < 0) THEN
         CPABORT("PROPERTIES%BSE: the tensor groups of the GW are not set")
      END IF

      n_spin = bs_env%n_spin
      n_RI = bs_env%n_RI
      ! the canonical orthogonalisation parks the linearly dependent directions above n_mo_retained
      n_mo_real = MIN(bs_env%n_mo_retained, bs_env%n_ao)

      CALL read_bse_section(bse_section, bse_env)
      CALL bse_resolve_memory_per_proc(bse_env%memory_per_proc_gb, bs_env%para_env, unit_nr)

      CALL link_window(bs_env, bse_env, n_mo_real, window, homo_red, virt_red, unit_nr)
      n_window = window%last_mo - window%first_mo + 1

      CALL link_screening(qs_env, bs_env, unit_nr, fm_M_inv_V_sqrt, fm_eps_inv)
      ! the GW kept W(0) for the screening and nl_3c for the slabs; each goes once read
      CALL cp_fm_release(bs_env%fm_W_MIC_freq_zero)

      ALLOCATE (fm_S_ia(n_spin), fm_S_ij(n_spin), fm_S_ab(n_spin))
      ALLOCATE (eigenval_qp(n_window, n_spin), eigenval_dft(n_window, n_spin), homo_full(n_spin))
      DO ispin = 1, n_spin
         CALL link_slabs(qs_env, bs_env, bse_env, ispin, window, homo_red(ispin), virt_red(ispin), fm_M_inv_V_sqrt, &
                         fm_S_ij(ispin), fm_S_ia(ispin), fm_S_ab(ispin), unit_nr)
         eigenval_qp(:, ispin) = bs_env%eigenval_GW(window%first_mo:window%last_mo, 1, ispin)
         eigenval_dft(:, ispin) = bs_env%eigenval_scf_Gamma(window%first_mo:window%last_mo, ispin)
         homo_full(ispin) = bs_env%n_occ(ispin)
      END DO
      CALL cp_fm_release(fm_M_inv_V_sqrt)
      CALL neighbor_list_3c_destroy(bs_env%nl_3c)

      CALL prepare_bse_env_from_bandstructure_gw(bse_env, fm_S_ia, fm_S_ij, fm_S_ab, fm_eps_inv, &
                                                 eigenval_qp, eigenval_dft, homo_full, homo_red, virt_red, &
                                                 n_RI, qs_env, unit_nr)
      CALL cp_fm_release(fm_eps_inv)
      CALL bse_solve(bse_env, bs_env%fm_mo_coeff_Gamma(1:n_spin), unit_nr)
      CALL bse_env_release(bse_env)

      CALL timestop(handle)

   END SUBROUTINE bse_from_bs_env

! **************************************************************************************************
!> \brief The active window on the DFT axis over the real (retained) MOs, one absolute range for
!>        every spin, and the levels inside it per spin: homo_red = n_occ - first_mo + 1 occupied,
!>        virt_red = last_mo - n_occ virtual, with n_occ the occupied MOs of the spin
!> \param bs_env DFT eigenvalues and occupations at Gamma
!> \param bse_env its cutoffs, read from the section already
!> \param n_mo_real MOs that are states, not linear-dependency placeholders
!> \param window first and last MO of the window
!> \param homo_red occupied levels in the window per spin, allocated here
!> \param virt_red virtual levels in the window per spin, allocated here
!> \param unit_nr output unit, positive on the writing rank only
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
!>        [M^-1 V^1/2]_QP = sum_R [M^-1]_QR [V^1/2]_RP, so that the slabs
!>        B^P_pq = sum_Q (pq|Q) [M^-1 V^1/2]_QP give
!>        sum_P B^P_pq B^P_rs = sum_QR (pq|Q) [M^-1 V M^-1]_QR (R|rs);
!>        with the stored W~_PQ = sum_RS [M^-1]_PR W^c_RS(iω=0) [M^-1]_SQ,
!>        ε^-1_eff,PQ = δ_PQ + sum_RS [V^-1/2]_PR W^c_RS(iω=0) [V^-1/2]_SQ, so that
!>        sum_PQ B^P_pq ε^-1_eff,PQ B^Q_rs = sum_RS (pq|R) [M^-1 (V + W^c(iω=0)) M^-1]_RS (S|rs),
!>        the full static W of the GW. The GW builds the same factor M^-1 V^1/2 per frequency from a
!>        Cholesky factor of V (compute_MinvVsqrt_Vsqrt); the symmetric roots here give V^-1/2 as
!>        well, which the stored W~ needs, with the eigenvalue threshold of the GW (eps_eigval_mat_RI).
!>        Two diagonalisations of an RI x RI matrix, once per run
!> \param qs_env for the 2c RI integrals
!> \param bs_env RI metric and its inverse, the stored static W~, the RI x RI matrix structure
!> \param unit_nr output unit, positive on the writing rank only
!> \param fm_M_inv_V_sqrt [M^-1 V^1/2]_QP, RI x RI, created here
!> \param fm_eps_inv ε^-1_eff,PQ, RI x RI, created here
! **************************************************************************************************
   SUBROUTINE link_screening(qs_env, bs_env, unit_nr, fm_M_inv_V_sqrt, fm_eps_inv)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      INTEGER, INTENT(IN)                                :: unit_nr
      TYPE(cp_fm_type), INTENT(OUT)                      :: fm_M_inv_V_sqrt, fm_eps_inv

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'link_screening'

      INTEGER                                            :: handle, n_dependent, n_RI
      TYPE(cp_fm_type)                                   :: fm_V_inv_sqrt, fm_V_sqrt, fm_Wc, fm_work
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:, :)     :: fm_M, fm_V
      TYPE(libint_potential_type)                        :: coulomb_op

      CALL timeset(routineN, handle)

      n_RI = bs_env%n_RI
      ! bare Coulomb V_PQ (coulomb_op keeps its default, the full 1/r) and the RI metric M_PQ
      CALL RI_2c_integral_mat(qs_env, fm_V, bs_env%fm_RI_RI%matrix_struct, n_RI, coulomb_op)
      CALL RI_2c_integral_mat(qs_env, fm_M, bs_env%fm_RI_RI%matrix_struct, n_RI, bs_env%ri_metric)

      ! [V^1/2]_PQ and [V^-1/2]_PQ by symmetric roots; eigenvalues of V under the GW's threshold are dropped
      CALL cp_fm_create(fm_V_sqrt, bs_env%fm_RI_RI%matrix_struct)
      CALL cp_fm_create(fm_V_inv_sqrt, bs_env%fm_RI_RI%matrix_struct)
      CALL cp_fm_create(fm_work, bs_env%fm_RI_RI%matrix_struct)
      CALL cp_fm_to_fm(fm_V(1, 1), fm_V_sqrt)
      CALL cp_fm_power(fm_V_sqrt, fm_work, 0.5_dp, bs_env%eps_eigval_mat_RI, n_dependent)
      CALL cp_fm_to_fm(fm_V(1, 1), fm_V_inv_sqrt)
      CALL cp_fm_power(fm_V_inv_sqrt, fm_work, -0.5_dp, bs_env%eps_eigval_mat_RI, n_dependent)
      IF (unit_nr > 0 .AND. n_dependent > 0) WRITE (unit_nr, '(T2,A4,T7,A,T71,I10)') 'BSE|', &
         'RI Coulomb directions dropped by the eigenvalue threshold', n_dependent

      ! [M^-1 V^1/2]_QP = sum_R [M^-1]_QR [V^1/2]_RP
      CALL cp_fm_create(fm_M_inv_V_sqrt, bs_env%fm_RI_RI%matrix_struct)
      CALL parallel_gemm("N", "N", n_RI, n_RI, n_RI, 1.0_dp, bs_env%fm_Minv_Gamma, fm_V_sqrt, &
                         0.0_dp, fm_M_inv_V_sqrt)

      ! W^c_PQ(iω=0) = sum_RS M_PR W~_RS M_SQ, in place on a copy of the stored W~
      CALL cp_fm_create(fm_Wc, bs_env%fm_RI_RI%matrix_struct)
      CALL cp_fm_to_fm(bs_env%fm_W_MIC_freq_zero, fm_Wc)
      CALL fm_contract_ABA(fm_M(1, 1), fm_Wc)
      ! ε^-1_eff,PQ = δ_PQ + sum_RS [V^-1/2]_PR W^c_RS(iω=0) [V^-1/2]_SQ, the sum in the work matrix
      CALL fm_contract_ABA(fm_V_inv_sqrt, fm_Wc, fm_work)
      CALL cp_fm_create(fm_eps_inv, bs_env%fm_RI_RI%matrix_struct)
      CALL cp_fm_set_all(fm_eps_inv, 0.0_dp, 1.0_dp)
      CALL cp_fm_scale_and_add(1.0_dp, fm_eps_inv, 1.0_dp, fm_work)

      CALL cp_fm_release(fm_V)
      CALL cp_fm_release(fm_M)
      CALL cp_fm_release(fm_V_sqrt)
      CALL cp_fm_release(fm_V_inv_sqrt)
      CALL cp_fm_release(fm_Wc)
      CALL cp_fm_release(fm_work)

      CALL timestop(handle)

   END SUBROUTINE link_screening

! **************************************************************************************************
!> \brief The window columns of the MO coefficients, C_μp for p in the window, as a 2-index
!>        tensor (AO atom | MO block) on the tensor group's grid, the MO blocking the caller's.
!>        The coefficients are a full matrix on all ranks, the tensor lives on the ranks of one
!>        tensor group, and no copy between the two exists. The routine takes the route of the GW's
!>        own matrices (cp_fm_to_dbt: full matrix, DBCSR matrix on all ranks, DBCSR matrix on
!>        the group, tensor) and creates the two DBCSR matrices of that route in the
!>        (AO atom | MO block) blocking, for which the GW holds no template. A small object,
!>        n_ao x n_window: the steps are there for the change of communicator, not for memory
!> \param bs_env MO coefficients at Gamma, AO block sizes, the tensor group
!> \param dbcsr_dist distribution of the SCF matrices, the template of the global matrix
!> \param ispin spin channel
!> \param window first and last MO of the active window
!> \param sizes_mo MO block sizes of the tensor, summing to the window
!> \param t_C the tensor, created here
! **************************************************************************************************
   SUBROUTINE link_mo_tensor(bs_env, dbcsr_dist, ispin, window, sizes_mo, t_C)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist
      INTEGER, INTENT(IN)                                :: ispin
      TYPE(mo_window_type), INTENT(IN)                   :: window
      INTEGER, DIMENSION(:), INTENT(IN)                  :: sizes_mo
      TYPE(dbt_type), INTENT(OUT)                        :: t_C

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'link_mo_tensor'

      INTEGER                                            :: handle, n_window
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: dist_1, dist_2
      INTEGER, DIMENSION(2)                              :: pdims_2d
      TYPE(cp_fm_struct_type), POINTER                   :: struct_C_window
      TYPE(cp_fm_type)                                   :: fm_C_window
      TYPE(dbcsr_distribution_type)                      :: dist_global, dist_group, dist_tmpl
      TYPE(dbcsr_type)                                   :: mat_global, mat_group
      TYPE(dbt_pgrid_type)                               :: pgrid_2d

      CALL timeset(routineN, handle)

      n_window = window%last_mo - window%first_mo + 1
      ! the MO blocks cover the window exactly
      IF (SUM(sizes_mo) /= n_window) CPABORT("MO block sizes do not sum to the active window")

      pdims_2d = 0
      CALL dbt_pgrid_create(bs_env%para_env_tensor, pdims_2d, pgrid_2d)
      CALL create_2c_tensor(t_C, dist_1, dist_2, pgrid_2d, bs_env%sizes_AO, sizes_mo, name="(AO | MO)")
      DEALLOCATE (dist_1, dist_2)
      CALL dbt_pgrid_destroy(pgrid_2d)

      ! the window columns as a full matrix, then the route of the GW's own matrices (cp_fm_to_dbt):
      ! global DBCSR, group DBCSR, tensor; the coefficients are never replicated
      NULLIFY (struct_C_window)
      CALL cp_fm_struct_create(struct_C_window, template_fmstruct=bs_env%fm_mo_coeff_Gamma(ispin)%matrix_struct, &
                               ncol_global=n_window)
      CALL cp_fm_create(fm_C_window, struct_C_window)
      CALL cp_fm_struct_release(struct_C_window)
      CALL cp_fm_to_fm_submat(bs_env%fm_mo_coeff_Gamma(ispin), fm_C_window, bs_env%n_ao, n_window, 1, window%first_mo, 1, 1)
      CALL link_matrix_round_robin(dbcsr_dist, bs_env%sizes_AO, sizes_mo, "C window", dist_global, mat_global)
      ! the group matrix on the grid of the GW's AO matrices of the tensor group
      CALL dbcsr_get_info(bs_env%mat_ao_ao_tensor%matrix, distribution=dist_tmpl)
      CALL link_matrix_round_robin(dist_tmpl, bs_env%sizes_AO, sizes_mo, "C window, group", dist_group, mat_group)
      CALL cp_fm_to_dbt(fm_C_window, mat_global, mat_group, t_C, bs_env%para_env, &
                        bs_env%para_env_tensor%num_pe, bs_env%eps_filter)

      CALL cp_fm_release(fm_C_window)
      CALL dbcsr_release(mat_global)
      CALL dbcsr_release(mat_group)
      CALL dbcsr_distribution_release(dist_global)
      CALL dbcsr_distribution_release(dist_group)

      CALL timestop(handle)

   END SUBROUTINE link_mo_tensor

! **************************************************************************************************
!> \brief A DBCSR matrix with the given blocking on the communicator and process grid of a template
!>        distribution, the blocks round robin over the grid rows and columns
!> \param dist_tmpl the template
!> \param row_sizes row block sizes
!> \param col_sizes column block sizes
!> \param matrix_name name of the matrix
!> \param matrix_dist the distribution of the matrix, created here, released after the matrix
!> \param new_matrix the matrix, created here
! **************************************************************************************************
   SUBROUTINE link_matrix_round_robin(dist_tmpl, row_sizes, col_sizes, matrix_name, matrix_dist, new_matrix)
      TYPE(dbcsr_distribution_type), INTENT(IN)          :: dist_tmpl
      INTEGER, DIMENSION(:), INTENT(IN)                  :: row_sizes, col_sizes
      CHARACTER(LEN=*), INTENT(IN)                       :: matrix_name
      TYPE(dbcsr_distribution_type), INTENT(OUT)         :: matrix_dist
      TYPE(dbcsr_type), INTENT(OUT)                      :: new_matrix

      INTEGER                                            :: i_block, npcols, nprows
      INTEGER, DIMENSION(:), POINTER                     :: col_blk, col_dist, row_blk, row_dist

      NULLIFY (col_blk, col_dist, row_blk, row_dist)
      CALL dbcsr_distribution_get(dist_tmpl, nprows=nprows, npcols=npcols)
      ALLOCATE (row_dist(SIZE(row_sizes)), col_dist(SIZE(col_sizes)), row_blk(SIZE(row_sizes)), col_blk(SIZE(col_sizes)))
      row_dist(:) = [(MOD(i_block - 1, nprows), i_block=1, SIZE(row_sizes))]
      col_dist(:) = [(MOD(i_block - 1, npcols), i_block=1, SIZE(col_sizes))]
      row_blk(:) = row_sizes(:)
      col_blk(:) = col_sizes(:)
      CALL dbcsr_distribution_new(matrix_dist, template=dist_tmpl, row_dist=row_dist, col_dist=col_dist)
      CALL dbcsr_create(new_matrix, name=matrix_name, dist=matrix_dist, matrix_type=dbcsr_type_no_symmetry, &
                        row_blk_size=row_blk, col_blk_size=col_blk)
      DEALLOCATE (row_dist, col_dist, row_blk, col_blk)

   END SUBROUTINE link_matrix_round_robin

! **************************************************************************************************
!> \brief Block sizes of the window index p: blocks of at most mo_block_size levels inside the
!>        occupied and inside the virtual sector, so that no block straddles the two
!> \param n_occ occupied levels in the window
!> \param n_virt virtual levels in the window
!> \param block_sizes levels per block, occupied blocks first, allocated here
!> \param n_block_occ occupied blocks
! **************************************************************************************************
   SUBROUTINE link_sector_blocks(n_occ, n_virt, block_sizes, n_block_occ)
      INTEGER, INTENT(IN)                                :: n_occ, n_virt
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: block_sizes
      INTEGER, INTENT(OUT)                               :: n_block_occ

      INTEGER                                            :: n_block_virt

      n_block_occ = (n_occ - 1)/mo_block_size + 1
      n_block_virt = (n_virt - 1)/mo_block_size + 1
      ALLOCATE (block_sizes(n_block_occ + n_block_virt))
      block_sizes(:) = mo_block_size
      block_sizes(n_block_occ) = n_occ - (n_block_occ - 1)*mo_block_size
      block_sizes(n_block_occ + n_block_virt) = n_virt - (n_block_virt - 1)*mo_block_size

   END SUBROUTINE link_sector_blocks

! **************************************************************************************************
!> \brief Contiguous RI atom ranges of at most max_RI_per_batch RI functions; an atom larger than
!>        that forms its own batch
!> \param sizes_RI RI functions per atom
!> \param max_RI_per_batch RI functions per batch the memory left allows
!> \param batch_first first atom per batch, allocated here
!> \param batch_last last atom per batch, allocated here
! **************************************************************************************************
   SUBROUTINE link_ri_batches(sizes_RI, max_RI_per_batch, batch_first, batch_last)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: sizes_RI
      INTEGER, INTENT(IN)                                :: max_RI_per_batch
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: batch_first, batch_last

      INTEGER                                            :: i_atom, n_atom, n_batch, n_RI_in_batch
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: first_atom, last_atom

      n_atom = SIZE(sizes_RI)
      ALLOCATE (first_atom(n_atom), last_atom(n_atom))
      n_batch = 0
      n_RI_in_batch = 0
      DO i_atom = 1, n_atom
         IF (n_RI_in_batch > 0 .AND. n_RI_in_batch + sizes_RI(i_atom) > max_RI_per_batch) THEN
            last_atom(n_batch) = i_atom - 1
            n_RI_in_batch = 0
         END IF
         IF (n_RI_in_batch == 0) THEN
            n_batch = n_batch + 1
            first_atom(n_batch) = i_atom
         END IF
         n_RI_in_batch = n_RI_in_batch + sizes_RI(i_atom)
      END DO
      last_atom(n_batch) = n_atom
      ALLOCATE (batch_first(n_batch), batch_last(n_batch))
      batch_first(:) = first_atom(1:n_batch)
      batch_last(:) = last_atom(1:n_batch)
      DEALLOCATE (first_atom, last_atom)

   END SUBROUTINE link_ri_batches

! **************************************************************************************************
!> \brief What MEMORY_PER_PROC leaves to the slab build per rank and the largest RI batch it allows:
!>        per RI function of a batch, on the ranks of its tensor group the 3c integrals (P|μν), the
!>        half-transform (P μ|q) and the product (P|qp), and on all ranks the block matrices of the
!>        batch, their global copies and the batch's full matrix; one atom per batch with a warning
!>        when one RI atom alone exceeds what is left
!> \param bs_env sizes and communicators
!> \param bse_env the resolved MEMORY_PER_PROC
!> \param n_window levels in the active window
!> \param unit_nr output unit, positive on the writing rank only
!> \param budget_bytes what is left, in bytes; the RI transform chunks its work matrix by it as well
!> \param max_RI_per_batch RI functions per batch within what is left, at least one
! **************************************************************************************************
   SUBROUTINE link_batch_budget(bs_env, bse_env, n_window, unit_nr, budget_bytes, max_RI_per_batch)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(bse_env_type), INTENT(IN)                     :: bse_env
      INTEGER, INTENT(IN)                                :: n_window, unit_nr
      INTEGER(KIND=int_8), INTENT(OUT)                   :: budget_bytes
      INTEGER, INTENT(OUT)                               :: max_RI_per_batch

      CHARACTER(LEN=16)                                  :: left_str, need_str
      INTEGER(KIND=int_8)                                :: bytes_per_RI, largest_batch_bytes, &
                                                            one_atom_bytes
      REAL(KIND=dp)                                      :: mem_left_GB

      mem_left_GB = bse_memory_left_GB(bse_env%memory_per_proc_gb, bs_env%para_env)
      budget_bytes = INT(mem_left_GB*1.0E9_dp, KIND=int_8)
      ! per RI function of a batch, all distributed: on the ranks of its tensor group the 3c
      ! integrals (P|μν), the half-transform (P μ|q) and the product (P|qp); on all ranks the block
      ! matrices of the batch, their global copies (two inside the group-to-global hop) and the
      ! batch's full matrix, four (P|qp); the contractions' own scratch is not counted
      bytes_per_RI = 8_int_8*((INT(bs_env%n_ao, int_8)**2 + INT(bs_env%n_ao, int_8)*INT(n_window, int_8) + &
                               INT(n_window, int_8)**2)/INT(bs_env%para_env_tensor%num_pe, int_8) + &
                              4_int_8*INT(n_window, int_8)**2/INT(bs_env%para_env%num_pe, int_8))
      max_RI_per_batch = INT(MAX(1_int_8, MIN(INT(SUM(bs_env%sizes_RI), int_8), budget_bytes/bytes_per_RI)))
      one_atom_bytes = INT(MAXVAL(bs_env%sizes_RI), int_8)*bytes_per_RI
      largest_batch_bytes = MAX(INT(max_RI_per_batch, int_8)*bytes_per_RI, one_atom_bytes)
      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A4,T7,A,T67,F14.3)') 'BSE|', 'Memory per MPI rank for the BSE, MEMORY_PER_PROC (GB)', &
            bse_env%memory_per_proc_gb
         WRITE (unit_nr, '(T2,A4,T7,A,T67,F14.3)') 'BSE|', 'Left for the 3c integrals after the memory held (GB)', &
            mem_left_GB
         WRITE (unit_nr, '(T2,A4,T7,A,T67,F14.3)') 'BSE|', 'Largest RI batch per MPI rank (GB)', &
            REAL(largest_batch_bytes, dp)*1.0E-9_dp
      END IF
      ! nothing smaller than one atom per batch exists: run with it and say so
      IF (one_atom_bytes > budget_bytes) THEN
         WRITE (need_str, '(F16.3)') REAL(one_atom_bytes, dp)*1.0E-9_dp
         WRITE (left_str, '(F16.3)') mem_left_GB
         CALL cp_warn(__LOCATION__, &
                      "PROPERTIES%BSE: the 3c integral batch of one RI atom needs "//TRIM(ADJUSTL(need_str))// &
                      " GB per MPI rank, "//TRIM(ADJUSTL(left_str))//" GB are left of MEMORY_PER_PROC; one atom "// &
                      "per batch. Raise MEMORY_PER_PROC or use more MPI ranks.")
      END IF

   END SUBROUTINE link_batch_budget

! **************************************************************************************************
!> \brief The three window-sized slabs of one spin, B^P_pq = sum_Q (pq|Q) [M^-1 V^1/2]_QP with
!>        (Q|pq) = sum_μν C_μp (Q|μν) C_νq, batched over RI atoms: the 3c integrals of the batch
!>        are rebuilt and transformed over ν, then over μ, by two tensor contractions on the batch's
!>        tensor group; the product goes to all ranks and its rows are copied into the slabs; the
!>        RI transform by M^-1 V^1/2 follows in column chunks. Column conventions of the solver core:
!>        ia = (i-1) n_virt + a, ij = (i-1) n_occ + j, ab = (a-1) n_virt + b
!> \param qs_env blacs context, distribution of the SCF matrices, environment of the 3c integrals
!> \param bs_env RI and AO sizes, the 3c neighbour list and tensor template, the tensor groups
!> \param bse_env its memory budget and memory check
!> \param ispin spin channel
!> \param window first and last MO of the active window
!> \param n_occ occupied levels in the window
!> \param n_virt virtual levels in the window
!> \param fm_M_inv_V_sqrt the RI transform [M^-1 V^1/2]_QP, RI x RI
!> \param fm_S_ij B^P_ij, RI x n_occ*n_occ, created here
!> \param fm_S_ia B^P_ia, RI x n_occ*n_virt, created here
!> \param fm_S_ab B^P_ab, RI x n_virt*n_virt, created here
!> \param unit_nr output unit, positive on the writing rank only
! **************************************************************************************************
   SUBROUTINE link_slabs(qs_env, bs_env, bse_env, ispin, window, n_occ, n_virt, fm_M_inv_V_sqrt, fm_S_ij, fm_S_ia, &
                         fm_S_ab, unit_nr)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(bse_env_type), INTENT(IN)                     :: bse_env
      INTEGER, INTENT(IN)                                :: ispin
      TYPE(mo_window_type), INTENT(IN)                   :: window
      INTEGER, INTENT(IN)                                :: n_occ, n_virt
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_M_inv_V_sqrt
      TYPE(cp_fm_type), INTENT(OUT)                      :: fm_S_ij, fm_S_ia, fm_S_ab
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'link_slabs'

      INTEGER                                            :: atom_first, atom_last, handle, i_batch, &
                                                            max_RI_per_batch, n_batch, &
                                                            n_block_occ, n_RI, n_window
      INTEGER(KIND=int_8)                                :: budget_bytes
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: batch_first, batch_last, dist_1, dist_2, &
                                                            dist_3, sizes_p, sizes_q
      INTEGER, DIMENSION(3)                              :: pdims_3d
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: struct_ab, struct_ia, struct_ij
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist
      TYPE(dbt_pgrid_type)                               :: pgrid_3d
      TYPE(dbt_type)                                     :: t_3c, t_3c_half, t_3c_mo, t_C_p, t_C_q

      CALL timeset(routineN, handle)

      NULLIFY (blacs_env, dbcsr_dist, struct_ab, struct_ia, struct_ij)
      CALL get_qs_env(qs_env, blacs_env=blacs_env, dbcsr_dist=dbcsr_dist)
      n_RI = bs_env%n_RI
      n_window = n_occ + n_virt
      CALL cp_fm_struct_create(struct_ij, para_env=bs_env%para_env, context=blacs_env, &
                               nrow_global=n_RI, ncol_global=n_occ*n_occ)
      CALL cp_fm_struct_create(struct_ia, para_env=bs_env%para_env, context=blacs_env, &
                               nrow_global=n_RI, ncol_global=n_occ*n_virt)
      CALL cp_fm_struct_create(struct_ab, para_env=bs_env%para_env, context=blacs_env, &
                               nrow_global=n_RI, ncol_global=n_virt*n_virt)
      ! the slabs collect (Q|pq) batch by batch; the RI transform at the end makes them B^P_pq in place
      CALL cp_fm_create(fm_S_ij, struct_ij, name="fm_S_ij", set_zero=.TRUE.)
      CALL cp_fm_create(fm_S_ia, struct_ia, name="fm_S_ia", set_zero=.TRUE.)
      CALL cp_fm_create(fm_S_ab, struct_ab, name="fm_S_ab", set_zero=.TRUE.)
      CALL cp_fm_struct_release(struct_ij)
      CALL cp_fm_struct_release(struct_ia)
      CALL cp_fm_struct_release(struct_ab)

      ! the window index q in two blocks, occupied and virtual, so that every (P, q, p) block of the
      ! product covers the whole q range of one pair type; p in blocks inside each sector
      ALLOCATE (sizes_q(2))
      sizes_q(:) = [n_occ, n_virt]
      CALL link_sector_blocks(n_occ, n_virt, sizes_p, n_block_occ)
      CALL link_mo_tensor(bs_env, dbcsr_dist, ispin, window, sizes_q, t_C_q)
      CALL link_mo_tensor(bs_env, dbcsr_dist, ispin, window, sizes_p, t_C_p)

      CALL link_batch_budget(bs_env, bse_env, n_window, unit_nr, budget_bytes, max_RI_per_batch)
      CALL link_ri_batches(bs_env%sizes_RI, max_RI_per_batch, batch_first, batch_last)
      n_batch = SIZE(batch_first)
      IF (unit_nr > 0) WRITE (unit_nr, '(T2,A4,T7,A,T71,I10)') 'BSE|', 'RI atom batches of the 3c integral transformation', n_batch

      ! (RI AO | MO) of the half-transform and (RI | MO MO) of the product on the tensor group's grid
      pdims_3d = 0
      CALL dbt_pgrid_create(bs_env%para_env_tensor, pdims_3d, pgrid_3d)
      CALL create_3c_tensor(t_3c_half, dist_1, dist_2, dist_3, pgrid_3d, bs_env%sizes_RI, bs_env%sizes_AO, &
                            sizes_q, map1=[1, 2], map2=[3], name="(RI AO | MO)")
      DEALLOCATE (dist_1, dist_2, dist_3)
      CALL dbt_pgrid_destroy(pgrid_3d)
      pdims_3d = 0
      CALL dbt_pgrid_create(bs_env%para_env_tensor, pdims_3d, pgrid_3d, &
                            tensor_dims=[bs_env%n_atom, SIZE(sizes_q), SIZE(sizes_p)])
      CALL create_3c_tensor(t_3c_mo, dist_1, dist_2, dist_3, pgrid_3d, bs_env%sizes_RI, sizes_q, sizes_p, &
                            map1=[1], map2=[2, 3], name="(RI | MO MO)")
      DEALLOCATE (dist_1, dist_2, dist_3)
      CALL dbt_pgrid_destroy(pgrid_3d)
      CALL dbt_create(bs_env%t_RI_AO__AO, t_3c)

      ! C_q contracts once per batch in a layout-compatible contraction: the batched context keeps its replicated
      ! copy across the batches; on the small operand only, so that no grid changes under the neighbour list
      ! of the 3c integrals. No context on C_p: no layout of (Pμ|q) fits both contractions, so the second remaps it
      CALL dbt_batched_contract_init(t_C_q)

      DO i_batch = 1, n_batch
         atom_first = batch_first(i_batch)
         atom_last = batch_last(i_batch)
         ! one tensor group per batch, round robin; every rank joins the hop into the slabs
         IF (MODULO(i_batch - 1, bs_env%num_tensor_groups) == bs_env%tensor_group_color) THEN
            CALL compute_3c_integrals(qs_env, bs_env, t_3c, atoms_RI=[atom_first, atom_last])
            ! (Pμ|q) = sum_ν (P|μν) C_νq
            CALL dbt_contract(alpha=1.0_dp, tensor_1=t_3c, tensor_2=t_C_q, beta=0.0_dp, tensor_3=t_3c_half, &
                              contract_1=[3], notcontract_1=[1, 2], map_1=[1, 2], &
                              contract_2=[1], notcontract_2=[2], map_2=[3], &
                              filter_eps=bs_env%eps_filter)
            CALL dbt_clear(t_3c)
            ! (P|qp) = sum_μ (Pμ|q) C_μp, distributed over the group
            CALL dbt_contract(alpha=1.0_dp, tensor_1=t_3c_half, tensor_2=t_C_p, beta=0.0_dp, tensor_3=t_3c_mo, &
                              contract_1=[2], notcontract_1=[1, 3], map_1=[1, 2], &
                              contract_2=[1], notcontract_2=[2], map_2=[3], &
                              filter_eps=bs_env%eps_filter)
            CALL dbt_clear(t_3c_half)
         END IF
         CALL link_slab_rows(bs_env, blacs_env, dbcsr_dist, t_3c_mo, n_occ, n_virt, sizes_p, n_block_occ, atom_first, atom_last, &
                             fm_S_ij, fm_S_ia, fm_S_ab)
         CALL dbt_clear(t_3c_mo)
      END DO

      CALL dbt_batched_contract_finalize(t_C_q)
      CALL dbt_destroy(t_3c)
      CALL dbt_destroy(t_3c_half)
      CALL dbt_destroy(t_3c_mo)
      CALL dbt_destroy(t_C_q)
      CALL dbt_destroy(t_C_p)
      DEALLOCATE (sizes_q, sizes_p, batch_first, batch_last)

      ! B^P_pq = sum_Q [M^-1 V^1/2]_QP (Q|pq), in column chunks, in place
      CALL apply_ri_transform(fm_M_inv_V_sqrt, fm_S_ij, n_RI, n_occ*n_occ, budget_bytes)
      CALL apply_ri_transform(fm_M_inv_V_sqrt, fm_S_ia, n_RI, n_occ*n_virt, budget_bytes)
      CALL apply_ri_transform(fm_M_inv_V_sqrt, fm_S_ab, n_RI, n_virt*n_virt, budget_bytes)

      CALL timestop(handle)

   END SUBROUTINE link_slabs

! **************************************************************************************************
!> \brief The batch product (P|qp) into the rows of its RI functions in the three slabs: the local
!>        (P, q, p) blocks into one block matrix per pair type on the tensor group's grid, the
!>        group-to-global hop of the GW, the block-cyclic copy, the row copy into the slab; every
!>        rank takes part, a rank outside the batch's group with no blocks. The pair (q, p) lands
!>        in column (p-1) n_q + q of its slab, n_q = n_occ for ij and n_virt for ia and ab:
!>        ij = (i-1) n_occ + j at (q, p) = (j, i), ia = (i-1) n_virt + a at (a, i),
!>        ab = (a-1) n_virt + b at (b, a)
!> \param bs_env RI offsets and sizes per atom, the communicator of all ranks
!> \param blacs_env context of the slabs
!> \param dbcsr_dist a distribution over all ranks, the template of the global block matrices
!> \param t_3c_mo (P|qp) of the RI atoms atom_first..atom_last, q in two blocks (occupied, virtual),
!>        p in sector blocks
!> \param n_occ occupied levels in the window
!> \param n_virt virtual levels in the window
!> \param sizes_p the p blocks, occupied blocks first
!> \param n_block_occ occupied blocks of p
!> \param atom_first first RI atom of the batch
!> \param atom_last last RI atom of the batch
!> \param fm_S_ij slab of the ij pairs; the rows of the batch are written
!> \param fm_S_ia slab of the ia pairs; the rows of the batch are written
!> \param fm_S_ab slab of the ab pairs; the rows of the batch are written
! **************************************************************************************************
   SUBROUTINE link_slab_rows(bs_env, blacs_env, dbcsr_dist, t_3c_mo, n_occ, n_virt, sizes_p, n_block_occ, atom_first, atom_last, &
                             fm_S_ij, fm_S_ia, fm_S_ab)
      TYPE(post_scf_bandstructure_type), POINTER         :: bs_env
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist
      TYPE(dbt_type), INTENT(INOUT)                      :: t_3c_mo
      INTEGER, INTENT(IN)                                :: n_occ, n_virt
      INTEGER, DIMENSION(:), INTENT(IN)                  :: sizes_p
      INTEGER, INTENT(IN)                                :: n_block_occ, atom_first, atom_last
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_S_ij, fm_S_ia, fm_S_ab

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'link_slab_rows'
      INTEGER, PARAMETER                                 :: pair_ab = 3, pair_ia = 2, pair_ij = 1

      INTEGER :: col_block, group_handle, handle, i_pair, my_rank, n_atom_batch, n_block_p, &
         n_RI_batch, npcol, nprow, owner_rank, p_block, q_block, rank, RI_first, row_block
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: cart_coords
      INTEGER, DIMENSION(3)                              :: blk_index, blk_size, n_col_blocks, &
                                                            n_cols_pair, p_block_offset, &
                                                            q_block_of_pair
      INTEGER, DIMENSION(:), POINTER                     :: col_dist, col_sizes, row_dist, row_sizes
      INTEGER, DIMENSION(:, :), POINTER                  :: pgrid
      LOGICAL                                            :: found
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: block_Pqp
      TYPE(cp_fm_struct_type), POINTER                   :: struct_batch
      TYPE(cp_fm_type)                                   :: fm_batch
      TYPE(dbcsr_distribution_type), DIMENSION(3)        :: dist_global, dist_group
      TYPE(dbcsr_type), DIMENSION(3)                     :: mat_global, mat_group
      TYPE(dbt_iterator_type)                            :: iter

      CALL timeset(routineN, handle)

      NULLIFY (col_dist, col_sizes, pgrid, row_dist, row_sizes, struct_batch)
      RI_first = bs_env%i_RI_start_from_atom(atom_first)
      n_RI_batch = bs_env%i_RI_end_from_atom(atom_last) - RI_first + 1
      n_atom_batch = atom_last - atom_first + 1
      n_block_p = SIZE(sizes_p)
      n_cols_pair = [n_occ*n_occ, n_occ*n_virt, n_virt*n_virt]
      n_col_blocks = [n_block_occ, n_block_occ, n_block_p - n_block_occ]
      q_block_of_pair = [1, 2, 2]
      ! column block col_block of a pair type is the p block p_block = p_block_offset + col_block
      p_block_offset = [0, 0, n_block_occ]

      ! the group grid of the tensor as the process grid of the local matrices, its ranks placed by
      ! their coordinates, so that every block the tensor holds is local to its matrix
      ASSOCIATE (cart => t_3c_mo%pgrid%mp_comm_2d)
         my_rank = cart%mepos
         group_handle = cart%get_handle()
         nprow = cart%num_pe_cart(1)
         npcol = cart%num_pe_cart(2)
         ALLOCATE (cart_coords(2, 0:cart%num_pe - 1), pgrid(0:nprow - 1, 0:npcol - 1))
         DO rank = 0, cart%num_pe - 1
            CALL cart%coords(rank, cart_coords(:, rank))
            pgrid(cart_coords(1, rank), cart_coords(2, rank)) = rank
         END DO
      END ASSOCIATE
      ALLOCATE (row_dist(n_atom_batch), row_sizes(n_atom_batch))
      DO row_block = 1, n_atom_batch
         CALL dbt_get_stored_coordinates(t_3c_mo, [atom_first + row_block - 1, 1, 1], owner_rank)
         row_dist(row_block) = cart_coords(1, owner_rank)
         row_sizes(row_block) = bs_env%sizes_RI(atom_first + row_block - 1)
      END DO
      DO i_pair = 1, 3
         ALLOCATE (col_dist(n_col_blocks(i_pair)), col_sizes(n_col_blocks(i_pair)))
         DO col_block = 1, n_col_blocks(i_pair)
            p_block = p_block_offset(i_pair) + col_block
            CALL dbt_get_stored_coordinates(t_3c_mo, [atom_first, q_block_of_pair(i_pair), p_block], owner_rank)
            col_dist(col_block) = cart_coords(2, owner_rank)
            ! the q range of the pair type times the p block
            col_sizes(col_block) = MERGE(n_occ, n_virt, i_pair == pair_ij)*sizes_p(p_block)
         END DO
         CALL dbcsr_distribution_new(dist_group(i_pair), group=group_handle, pgrid=pgrid, &
                                     row_dist=row_dist, col_dist=col_dist)
         CALL dbcsr_create(mat_group(i_pair), name="(P|qp) batch, group", dist=dist_group(i_pair), &
                           matrix_type=dbcsr_type_no_symmetry, row_blk_size=row_sizes, col_blk_size=col_sizes)
         ! the threaded put_block below cannot create the work matrix itself
         CALL dbcsr_work_create(mat_group(i_pair), work_mutable=.TRUE.)
         ! the global matrix of the batch, round robin over the process grid of the slabs
         CALL link_matrix_round_robin(dbcsr_dist, row_sizes, col_sizes, "(P|qp) batch", dist_global(i_pair), mat_global(i_pair))
         DEALLOCATE (col_dist, col_sizes)
      END DO
      DEALLOCATE (pgrid, cart_coords, row_dist, row_sizes)

      ! the local (P, q, p) blocks as (P, qp) blocks, q fastest, the column order of the slabs; the
      ! iterator needs the full thread team
!$OMP PARALLEL DEFAULT(NONE) SHARED(t_3c_mo, mat_group, n_block_occ, atom_first, p_block_offset, my_rank) &
!$OMP PRIVATE(iter, blk_index, blk_size, block_Pqp, found, q_block, p_block, i_pair, row_block, col_block, owner_rank)
      CALL dbt_iterator_start(iter, t_3c_mo)
      DO WHILE (dbt_iterator_blocks_left(iter))
         CALL dbt_iterator_next_block(iter, blk_index, blk_size=blk_size)
         q_block = blk_index(2)
         p_block = blk_index(3)
         ! pair types of the (q block, p block) blocks: (occ, occ) ij, (virt, occ) ia, (virt, virt) ab;
         ! (occ, virt) is the ai block, which no slab holds
         IF (q_block == 1 .AND. p_block > n_block_occ) CYCLE
         IF (q_block == 1) THEN
            i_pair = pair_ij
         ELSE IF (p_block <= n_block_occ) THEN
            i_pair = pair_ia
         ELSE
            i_pair = pair_ab
         END IF
         ALLOCATE (block_Pqp(blk_size(1), blk_size(2), blk_size(3)))
         CALL dbt_get_block(t_3c_mo, blk_index, blk_size, block_Pqp, found)
         IF (found) THEN
            row_block = blk_index(1) - atom_first + 1
            col_block = p_block - p_block_offset(i_pair)
            CALL dbcsr_get_stored_coordinates(mat_group(i_pair), row_block, col_block, owner_rank)
            ! the group matrix was laid out by the tensor's coordinates: the rank holding the block owns it
            IF (owner_rank /= my_rank) CPABORT("(P|qp) block is not local to its group matrix")
            CALL dbcsr_put_block(mat_group(i_pair), row_block, col_block, &
                                 RESHAPE(block_Pqp, [blk_size(1), blk_size(2)*blk_size(3)]))
         END IF
         DEALLOCATE (block_Pqp)
      END DO
      CALL dbt_iterator_stop(iter)
!$OMP END PARALLEL

      DO i_pair = 1, 3
         CALL dbcsr_finalize(mat_group(i_pair))
         CALL cp_dbcsr_from_subgroups(mat_group(i_pair), mat_global(i_pair), bs_env%para_env)
         CALL cp_fm_struct_create(struct_batch, para_env=bs_env%para_env, context=blacs_env, &
                                  nrow_global=n_RI_batch, ncol_global=n_cols_pair(i_pair))
         CALL cp_fm_create(fm_batch, struct_batch)
         CALL cp_fm_struct_release(struct_batch)
         CALL copy_dbcsr_to_fm(mat_global(i_pair), fm_batch)
         ! known bottleneck: cp_fm_to_fm_submat copies column by column, one pdcopy per pair column; the
         ! memory it takes on single ranks can stay until the end of the run and set their peak. A block
         ! redistribution of fm_batch into the rows of the slab would replace it
         SELECT CASE (i_pair)
         CASE (pair_ij)
            CALL cp_fm_to_fm_submat(fm_batch, fm_S_ij, n_RI_batch, n_cols_pair(i_pair), 1, 1, RI_first, 1)
         CASE (pair_ia)
            CALL cp_fm_to_fm_submat(fm_batch, fm_S_ia, n_RI_batch, n_cols_pair(i_pair), 1, 1, RI_first, 1)
         CASE (pair_ab)
            CALL cp_fm_to_fm_submat(fm_batch, fm_S_ab, n_RI_batch, n_cols_pair(i_pair), 1, 1, RI_first, 1)
         END SELECT
         CALL cp_fm_release(fm_batch)
         CALL dbcsr_release(mat_group(i_pair))
         CALL dbcsr_release(mat_global(i_pair))
         CALL dbcsr_distribution_release(dist_group(i_pair))
         CALL dbcsr_distribution_release(dist_global(i_pair))
      END DO

      CALL timestop(handle)

   END SUBROUTINE link_slab_rows

! **************************************************************************************************
!> \brief B^P_c = sum_Q [M^-1 V^1/2]_QP (Q|c) for every column c of a slab, in place: the columns pass in
!>        chunks through one RI x chunk work matrix of at most budget_bytes per rank,
!>        max_cols_chunk = budget_bytes N_all/(8 n_RI) with N_all the ranks of the slab
!> \param fm_M_inv_V_sqrt the RI transform [M^-1 V^1/2]_QP, RI x RI
!> \param fm_slab (Q|c) on entry, B^P_c on exit
!> \param n_RI RI functions, the rows of the slab
!> \param n_cols columns of the slab
!> \param budget_bytes per-rank budget of the work matrix, which is distributed over the ranks
! **************************************************************************************************
   SUBROUTINE apply_ri_transform(fm_M_inv_V_sqrt, fm_slab, n_RI, n_cols, budget_bytes)
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_M_inv_V_sqrt, fm_slab
      INTEGER, INTENT(IN)                                :: n_RI, n_cols
      INTEGER(KIND=int_8), INTENT(IN)                    :: budget_bytes

      CHARACTER(LEN=*), PARAMETER :: routineN = 'apply_ri_transform'

      INTEGER                                            :: col_first, handle, max_cols_chunk, &
                                                            n_cols_chunk
      TYPE(cp_fm_struct_type), POINTER                   :: struct_work
      TYPE(cp_fm_type)                                   :: fm_work

      CALL timeset(routineN, handle)

      max_cols_chunk = INT(MAX(1_int_8, MIN(INT(n_cols, int_8), &
                                            budget_bytes*INT(fm_slab%matrix_struct%para_env%num_pe, int_8) &
                                            /(8_int_8*INT(n_RI, int_8)))))
      NULLIFY (struct_work)
      CALL cp_fm_struct_create(struct_work, para_env=fm_slab%matrix_struct%para_env, &
                               context=fm_slab%matrix_struct%context, nrow_global=n_RI, ncol_global=max_cols_chunk)
      CALL cp_fm_create(fm_work, struct_work)
      CALL cp_fm_struct_release(struct_work)
      DO col_first = 1, n_cols, max_cols_chunk
         n_cols_chunk = MIN(max_cols_chunk, n_cols - col_first + 1)
         ! work_P,c = sum_Q [M^-1 V^1/2]_QP slab_Q,c for the columns c of this chunk, then back into the
         ! slab; the copy back is column by column (one pdcopy each), see link_slab_rows
         CALL parallel_gemm("T", "N", n_RI, n_cols_chunk, n_RI, 1.0_dp, fm_M_inv_V_sqrt, fm_slab, 0.0_dp, fm_work, &
                            b_first_col=col_first, c_first_col=1)
         CALL cp_fm_to_fm_submat(fm_work, fm_slab, nrow=n_RI, ncol=n_cols_chunk, s_firstrow=1, s_firstcol=1, &
                                 t_firstrow=1, t_firstcol=col_first)
      END DO
      CALL cp_fm_release(fm_work)

      CALL timestop(handle)

   END SUBROUTINE apply_ri_transform

END MODULE bse_bs_env_link
