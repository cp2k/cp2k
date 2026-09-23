!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Matrix-free application of the BSE matrices A and B to trial vectors from RI slabs that
!>        are sliced along the RI index over all MPI ranks
!> \par History
!>      09.2026 created [Maximilian Graml]
! **************************************************************************************************
MODULE bse_matvec

   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_get_submatrix,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_to_fm_submat_general,&
                                              cp_fm_type
   USE input_constants,                 ONLY: bse_precond_full_diag
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_mem_avail_per_rank_GB,&
                                              mp_para_env_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'bse_matvec'

   PUBLIC :: bse_matvec_env_type, bse_matvec_create, bse_matvec_release, bse_matvec_apply, &
             bse_matvec_diagonal, bse_matvec_subblock, bse_matvec_vector_struct, bse_matvec_selfcheck

! **************************************************************************************************
!> \brief RI slabs of one spin channel, every rank holding n_ri_loc whole RI slices
!> \param alpha prefactor of the exchange term (2 singlet, 0 triplet)
!> \param w_fac prefactor of the screened term (0 switches it off)
!> \param eps_diff ε_a - ε_i at compound index ia = (i-1)*virt + a
!> \param B_ia B^P_ia at (a, i, P)
!> \param B_bar_ia \bar{B}^P_ia at (a, i, P), only for do_abba
!> \param B_bar_ij \bar{B}^P_ij at (j, i, P)
!> \param B_ab B^P_ab at (b, a, P)
!> \param row_count rows of a trial vector owned by each rank, in rank order
!> \param row_displ first row of each rank minus one, so that row_count and row_displ describe
!>        the contiguous ia range per rank that allgatherv and sum_scatter need
!> \param blacs_env npe x 1 process grid of the slabs and of all trial vectors
!> \param block_cols trial vectors per pass in bse_matvec_apply, -1 to size it from free memory
!> \param mem_fraction share of the free memory per rank the iterative solver may use, for the
!>        pass width at block_cols -1 and for the subspace ceiling's default budget
! **************************************************************************************************
   TYPE bse_matvec_env_type
      INTEGER                                            :: homo = 0, virt = 0, n_ov = 0, &
                                                            n_ri_loc = 0, block_cols = -1
      ! mem_fraction: half of the free memory, since the probe counts reclaimable cache and is optimistic
      REAL(KIND=dp)                                      :: alpha = 0.0_dp, w_fac = 0.0_dp, &
                                                            mem_fraction = 0.5_dp
      LOGICAL                                            :: do_abba = .FALSE.
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: row_count, row_displ
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eps_diff
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: B_ia, B_bar_ia, B_bar_ij, B_ab
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env => NULL()
      TYPE(mp_para_env_type), POINTER                    :: para_env => NULL()
   END TYPE bse_matvec_env_type

CONTAINS

! **************************************************************************************************
!> \brief Moves the RI slabs onto an npe x 1 process grid, such that every rank owns whole RI
!>        slices, and stores them as local 3-index arrays: the contraction over each P in
!>        bse_matvec_apply is then a local DGEMM and only the final sum over P crosses the ranks
!> \param mv_env the environment created here
!> \param fm_S_ia B^P_ia, N_RI x homo*virt
!> \param fm_S_bar_ij \bar{B}^P_ij, the slab that enters W_ij,ab, N_RI x homo*homo
!> \param fm_S_ab B^P_ab, N_RI x virt*virt
!> \param fm_S_bar_ia \bar{B}^P_ia, the slab that enters W_ib,aj, N_RI x homo*virt
!> \param eps_reduced quasiparticle energies of the homo+virt active levels
!> \param homo occupied levels of the active window
!> \param virt virtual levels of the active window
!> \param alpha prefactor of the exchange term, 2 singlet, 0 triplet
!> \param w_fac prefactor of the screened term, 0 switches it off
!> \param do_abba slice \bar{B}^P_ia as well, for the application of B
!> \param unit_nr output unit, positive on the writing rank only
!> \param block_cols trial vectors per pass of bse_matvec_apply; absent or -1 sizes it from the free memory
! **************************************************************************************************
   SUBROUTINE bse_matvec_create(mv_env, fm_S_ia, fm_S_bar_ij, fm_S_ab, fm_S_bar_ia, eps_reduced, &
                                homo, virt, alpha, w_fac, do_abba, unit_nr, block_cols)

      TYPE(bse_matvec_env_type), INTENT(OUT)             :: mv_env
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_S_ia, fm_S_bar_ij, fm_S_ab, &
                                                            fm_S_bar_ia
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: eps_reduced
      INTEGER, INTENT(IN)                                :: homo, virt
      REAL(KIND=dp), INTENT(IN)                          :: alpha, w_fac
      LOGICAL, INTENT(IN)                                :: do_abba
      INTEGER, INTENT(IN)                                :: unit_nr
      INTEGER, INTENT(IN), OPTIONAL                      :: block_cols

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'bse_matvec_create'

      CHARACTER(LEN=12)                                  :: n_idle_str, n_ri_str, npe_str
      INTEGER                                            :: a, handle, i, n_ri, n_row_block
      REAL(KIND=dp)                                      :: mem_slabs_GB

      CALL timeset(routineN, handle)

      mv_env%homo = homo
      mv_env%virt = virt
      IF (PRESENT(block_cols)) mv_env%block_cols = block_cols
      mv_env%n_ov = homo*virt
      mv_env%alpha = alpha
      mv_env%w_fac = w_fac
      mv_env%do_abba = do_abba
      mv_env%para_env => fm_S_ia%matrix_struct%para_env

      CALL cp_blacs_env_create(mv_env%blacs_env, mv_env%para_env, &
                               grid_2d=[mv_env%para_env%num_pe, 1])

      ! one contiguous ia range per rank, so that a trial vector is gathered with allgatherv and
      ! its image scattered back with sum_scatter, rather than replicated by two allreduces
      n_row_block = contiguous_row_block(mv_env%n_ov, mv_env%para_env%num_pe)
      ALLOCATE (mv_env%row_count(mv_env%para_env%num_pe), mv_env%row_displ(mv_env%para_env%num_pe))
      DO i = 1, mv_env%para_env%num_pe
         mv_env%row_displ(i) = MIN(mv_env%n_ov, (i - 1)*n_row_block)
         mv_env%row_count(i) = MIN(mv_env%n_ov, i*n_row_block) - mv_env%row_displ(i)
      END DO

      ! eps_diff_ia = ε_a - ε_i at the compound index ia = (i-1)*virt + a, the row index of every trial vector
      ALLOCATE (mv_env%eps_diff(mv_env%n_ov))
      DO i = 1, homo
         DO a = 1, virt
            mv_env%eps_diff((i - 1)*virt + a) = eps_reduced(homo + a) - eps_reduced(i)
         END DO
      END DO

      CALL slice_slab(mv_env, fm_S_ia, virt, homo, mv_env%B_ia)
      IF (w_fac /= 0.0_dp) THEN
         CALL slice_slab(mv_env, fm_S_bar_ij, homo, homo, mv_env%B_bar_ij)
         CALL slice_slab(mv_env, fm_S_ab, virt, virt, mv_env%B_ab)
         IF (do_abba) CALL slice_slab(mv_env, fm_S_bar_ia, virt, homo, mv_env%B_bar_ia)
      END IF
      mv_env%n_ri_loc = SIZE(mv_env%B_ia, 3)

      CALL cp_fm_get_info(fm_S_ia, nrow_global=n_ri)
      IF (mv_env%para_env%num_pe > n_ri .AND. unit_nr > 0) THEN
         WRITE (npe_str, '(I0)') mv_env%para_env%num_pe
         WRITE (n_ri_str, '(I0)') n_ri
         WRITE (n_idle_str, '(I0)') mv_env%para_env%num_pe - n_ri
         CALL cp_warn(__LOCATION__, &
                      "BSE iterative solver: more MPI ranks ("//TRIM(npe_str)// &
                      ") than RI basis functions ("//TRIM(n_ri_str)//"); "// &
                      TRIM(n_idle_str)//" ranks hold no RI slice and idle in the kernel application.")
      END IF

      IF (unit_nr > 0) THEN
         mem_slabs_GB = REAL(homo*virt, dp)
         IF (w_fac /= 0.0_dp) THEN
            mem_slabs_GB = mem_slabs_GB + REAL(homo, dp)**2 + REAL(virt, dp)**2
            IF (do_abba) mem_slabs_GB = mem_slabs_GB + REAL(homo*virt, dp)
         END IF
         mem_slabs_GB = mem_slabs_GB*REAL(mv_env%n_ri_loc, dp)*8.0E-9_dp
         WRITE (unit_nr, '(T2,A4,T7,A,T71,I10)') 'BSE|', &
            'Max. number of RI functions per MPI rank', mv_env%n_ri_loc
         WRITE (unit_nr, '(T2,A4,T7,A,T67,F14.3)') 'BSE|', &
            'Memory of the RI slabs per MPI rank (GB)', mem_slabs_GB
      END IF

      CALL timestop(handle)

   END SUBROUTINE bse_matvec_create

! **************************************************************************************************
!> \brief Rows of a trial vector per rank, the smallest block that leaves no rank beyond the last
!> \param n_ov number of transitions, the global row count
!> \param npe ranks of the npe x 1 grid
!> \return rows of every full block; the last rank takes what remains, possibly none
! **************************************************************************************************
   PURE FUNCTION contiguous_row_block(n_ov, npe) RESULT(n_row_block)

      INTEGER, INTENT(IN)                                :: n_ov, npe
      INTEGER                                            :: n_row_block

      n_row_block = (n_ov + npe - 1)/npe

   END FUNCTION contiguous_row_block

! **************************************************************************************************
!> \brief Redistributes one N_RI x n_fast*n_slow slab to whole RI rows per rank and copies the owned
!>        rows into a contiguous (n_fast, n_slow, n_ri_loc) array, slab_loc(x, y, p) = slab(P_p, (y-1)*n_fast + x)
!>        for the p-th owned RI index P_p, so that every slice slab_loc(:, :, p) is a n_fast x n_slow DGEMM operand
!> \param mv_env process grid and communicator of the slices
!> \param fm_slab the slab on the grid of the caller, N_RI x n_fast*n_slow
!> \param n_fast fast index of the pair index of the slab
!> \param n_slow slow index of the pair index of the slab
!> \param slab_loc the owned slices, allocated here
! **************************************************************************************************
   SUBROUTINE slice_slab(mv_env, fm_slab, n_fast, n_slow, slab_loc)

      TYPE(bse_matvec_env_type), INTENT(IN)              :: mv_env
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_slab
      INTEGER, INTENT(IN)                                :: n_fast, n_slow
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT)                                     :: slab_loc

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'slice_slab'

      INTEGER                                            :: handle, n_ri, n_ri_loc, npair, p
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: fm_sliced

      CALL timeset(routineN, handle)

      ! row count from the struct: the slabs may carry fewer rows than the full RI basis
      CALL cp_fm_get_info(fm_slab, nrow_global=n_ri, ncol_global=npair)
      CPASSERT(npair == n_fast*n_slow)

      ! one-row blocks on the npe x 1 grid spread the RI rows cyclically over the ranks; the copy
      ! below turns each strided local_data row into a contiguous (n_fast, n_slow) slice
      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, para_env=mv_env%para_env, context=mv_env%blacs_env, &
                               nrow_global=n_ri, ncol_global=npair, nrow_block=1, &
                               force_block=.TRUE.)
      CALL cp_fm_create(fm_sliced, fm_struct, name="fm_slab_ri_sliced")
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_to_fm_submat_general(fm_slab, fm_sliced, n_ri, npair, 1, 1, 1, 1, &
                                      fm_slab%matrix_struct%context)

      CALL cp_fm_get_info(fm_sliced, nrow_local=n_ri_loc)
      ALLOCATE (slab_loc(n_fast, n_slow, n_ri_loc))
      DO p = 1, n_ri_loc
         slab_loc(:, :, p) = RESHAPE(fm_sliced%local_data(p, 1:npair), [n_fast, n_slow])
      END DO
      CALL cp_fm_release(fm_sliced)

      CALL timestop(handle)

   END SUBROUTINE slice_slab

! **************************************************************************************************
!> \brief Frees the sliced slabs, the transition energies and the process grid of mv_env
!> \param mv_env the environment released
! **************************************************************************************************
   SUBROUTINE bse_matvec_release(mv_env)

      TYPE(bse_matvec_env_type), INTENT(INOUT)           :: mv_env

      IF (ALLOCATED(mv_env%row_count)) DEALLOCATE (mv_env%row_count)
      IF (ALLOCATED(mv_env%row_displ)) DEALLOCATE (mv_env%row_displ)
      IF (ALLOCATED(mv_env%eps_diff)) DEALLOCATE (mv_env%eps_diff)
      IF (ALLOCATED(mv_env%B_ia)) DEALLOCATE (mv_env%B_ia)
      IF (ALLOCATED(mv_env%B_bar_ia)) DEALLOCATE (mv_env%B_bar_ia)
      IF (ALLOCATED(mv_env%B_bar_ij)) DEALLOCATE (mv_env%B_bar_ij)
      IF (ALLOCATED(mv_env%B_ab)) DEALLOCATE (mv_env%B_ab)
      IF (ASSOCIATED(mv_env%blacs_env)) CALL cp_blacs_env_release(mv_env%blacs_env)
      NULLIFY (mv_env%para_env)

   END SUBROUTINE bse_matvec_release

! **************************************************************************************************
!> \brief Matrix structure of a block of trial vectors: rows ia distributed, all columns local
!> \param mv_env process grid and communicator of the trial vectors
!> \param ncol_global number of trial vectors
!> \param fm_struct created here, released by the caller
! **************************************************************************************************
   SUBROUTINE bse_matvec_vector_struct(mv_env, ncol_global, fm_struct)

      TYPE(bse_matvec_env_type), INTENT(IN)              :: mv_env
      INTEGER, INTENT(IN)                                :: ncol_global
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct

      NULLIFY (fm_struct)
      ! the row block matches mv_env%row_count, which bse_matvec_apply's collectives rely on
      CALL cp_fm_struct_create(fm_struct, para_env=mv_env%para_env, context=mv_env%blacs_env, &
                               nrow_global=mv_env%n_ov, ncol_global=ncol_global, &
                               nrow_block=contiguous_row_block(mv_env%n_ov, mv_env%para_env%num_pe), &
                               force_block=.TRUE.)

   END SUBROUTINE bse_matvec_vector_struct

! **************************************************************************************************
!> \brief Applies A (and B) to ncol trial vectors without forming an N_ov x N_ov object,
!>        (A Z)_ia = (ε_a-ε_i) Z_ia + α sum_P B^P_ia sum_jb B^P_jb Z_jb - w sum_Pjb \bar{B}^P_ij B^P_ab Z_jb
!>        (B Z)_ia =                  α sum_P B^P_ia sum_jb B^P_jb Z_jb - w sum_Pjb \bar{B}^P_ib B^P_ja Z_jb
!>        α: exchange prefactor (2 singlet, 0 triplet), w: prefactor of the screened term,
!>        \bar{B}: slab bound by the caller of bse_matvec_create (B itself for TDHF and ALPHA
!>        screening, sum_Q [1+Q(0)]^-1_PQ B^Q otherwise); the sums over P run over the local RI
!>        slices and are completed by a sum over the ranks
!> \param mv_env slabs, prefactors and transition energies
!> \param fm_Z trial vectors, columns first_col .. first_col+ncol-1 are read
!> \param first_col first trial vector read, and first column written
!> \param ncol number of trial vectors
!> \param fm_AZ receives A Z in the same columns
!> \param fm_BZ receives B Z in the same columns, or from first_col_BZ on
!> \param first_col_BZ first column of fm_BZ written, first_col by default
! **************************************************************************************************
   SUBROUTINE bse_matvec_apply(mv_env, fm_Z, first_col, ncol, fm_AZ, fm_BZ, first_col_BZ)

      TYPE(bse_matvec_env_type), INTENT(IN), TARGET      :: mv_env
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_Z
      INTEGER, INTENT(IN)                                :: first_col, ncol
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_AZ
      TYPE(cp_fm_type), INTENT(IN), OPTIONAL             :: fm_BZ
      INTEGER, INTENT(IN), OPTIONAL                      :: first_col_BZ

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'bse_matvec_apply'

      INTEGER                                            :: col, col_shift_B, handle, handle2, homo, &
                                                            ia, iloc, k, me, n_ov, n_ri_loc, &
                                                            n_work_bufs, nb, nb_max, nrow_local, &
                                                            p, virt
      INTEGER, DIMENSION(:), POINTER                     :: row_indices
      LOGICAL                                            :: do_B
      REAL(KIND=dp)                                      :: mem_avail_GB
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: R_loc
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), TARGET   :: BabZ_buf, RA_buf, RB_buf, Z_buf
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: BjaZ_ab, T_Pk
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         POINTER                                         :: B_ia_2d, BabZ_wide, RA, RB, Z, Z_wide
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :, :), &
         POINTER                                         :: BabZ_3d, RA_3d, RB_3d, Z_3d

      CALL timeset(routineN, handle)

      homo = mv_env%homo
      virt = mv_env%virt
      n_ov = mv_env%n_ov
      n_ri_loc = mv_env%n_ri_loc
      do_B = PRESENT(fm_BZ)
      IF (do_B) THEN
         ! B Z needs the \bar{B}^P_ia slab, which is sliced for do_abba only
         CPASSERT(mv_env%do_abba)
      END IF
      col_shift_B = 0
      IF (PRESENT(first_col_BZ)) col_shift_B = first_col_BZ - first_col

      CALL cp_fm_get_info(fm_Z, nrow_local=nrow_local, row_indices=row_indices)
      ! the collectives below read the rows of one rank as one range, as bse_matvec_vector_struct
      ! lays them out; a vector built elsewhere would scatter into the wrong rows
      me = mv_env%para_env%mepos + 1
      CPASSERT(nrow_local == mv_env%row_count(me))
      CPASSERT(nrow_local == 0 .OR. row_indices(1) == mv_env%row_displ(me) + 1)
      ALLOCATE (R_loc(nrow_local))

      ! columns per pass from the memory of the replicated work arrays Z, B_ab Z, R^A (and R^B)
      n_work_bufs = 3
      IF (do_B) n_work_bufs = 4
      IF (mv_env%block_cols > 0) THEN
         ! pinned by input: the replicated work arrays are then a known n_work_bufs*n_ov*nb doubles
         nb_max = MIN(ncol, mv_env%block_cols)
      ELSE
         CALL mp_mem_avail_per_rank_GB(mv_env%para_env, mem_avail_GB)
         nb_max = ncol
         IF (mem_avail_GB > 0.0_dp) THEN
            nb_max = INT(MIN(REAL(ncol, dp), &
                             mv_env%mem_fraction*mem_avail_GB*1.0E9_dp/(8.0_dp*REAL(n_work_bufs, dp)*REAL(n_ov, dp))))
         END IF
      END IF
      nb_max = MAX(nb_max, 1)

      ALLOCATE (Z_buf(n_ov*nb_max), RA_buf(n_ov*nb_max), BabZ_buf(n_ov*nb_max))
      IF (do_B) THEN
         ALLOCATE (RB_buf(n_ov*nb_max), BjaZ_ab(virt, virt))
      END IF
      ALLOCATE (T_Pk(n_ri_loc, nb_max))
      B_ia_2d(1:n_ov, 1:n_ri_loc) => mv_env%B_ia

      DO col = first_col, first_col + ncol - 1, nb_max
         nb = MIN(nb_max, first_col + ncol - col)
         Z(1:n_ov, 1:nb) => Z_buf(1:n_ov*nb)
         Z_3d(1:virt, 1:homo, 1:nb) => Z_buf(1:n_ov*nb)
         Z_wide(1:virt, 1:homo*nb) => Z_buf(1:n_ov*nb)
         RA(1:n_ov, 1:nb) => RA_buf(1:n_ov*nb)
         RA_3d(1:virt, 1:homo, 1:nb) => RA_buf(1:n_ov*nb)
         BabZ_3d(1:virt, 1:homo, 1:nb) => BabZ_buf(1:n_ov*nb)
         BabZ_wide(1:virt, 1:homo*nb) => BabZ_buf(1:n_ov*nb)
         IF (do_B) THEN
            RB(1:n_ov, 1:nb) => RB_buf(1:n_ov*nb)
            RB_3d(1:virt, 1:homo, 1:nb) => RB_buf(1:n_ov*nb)
         END IF

         ! Z_ia,k on every rank: every rank owns a contiguous ia range, so the column arrives sorted
         CALL timeset(routineN//"_gather", handle2)
         DO k = 1, nb
            CALL mv_env%para_env%allgatherv(fm_Z%local_data(1:nrow_local, col + k - 1), Z(:, k), &
                                            mv_env%row_count, mv_env%row_displ)
         END DO
         CALL timestop(handle2)

         ! exchange, identical in A Z and B Z: R_ia,k = α sum_P B^P_ia t_Pk, t_Pk = sum_jb B^P_jb Z_jb,k
         CALL timeset(routineN//"_exchange", handle2)
         RA(:, :) = 0.0_dp
         IF (mv_env%alpha /= 0.0_dp .AND. n_ri_loc > 0) THEN
            CALL DGEMM('T', 'N', n_ri_loc, nb, n_ov, 1.0_dp, B_ia_2d, n_ov, Z, n_ov, &
                       0.0_dp, T_Pk, n_ri_loc)
            CALL DGEMM('N', 'N', n_ov, nb, n_ri_loc, mv_env%alpha, B_ia_2d, n_ov, T_Pk, n_ri_loc, &
                       0.0_dp, RA, n_ov)
         END IF
         IF (do_B) RB(:, :) = RA(:, :)
         CALL timestop(handle2)

         CALL timeset(routineN//"_W", handle2)
         IF (mv_env%w_fac /= 0.0_dp) THEN
            DO p = 1, n_ri_loc
               ! A: R_ia,k -= w sum_jb \bar{B}^P_ij B^P_ab Z_jb,k, BabZ_ib,k = sum_b B^P_ab Z_jb,k first
               CALL DGEMM('T', 'N', virt, homo*nb, virt, 1.0_dp, mv_env%B_ab(:, :, p), virt, &
                          Z_wide, virt, 0.0_dp, BabZ_wide, virt)
               DO k = 1, nb
                  CALL DGEMM('N', 'N', virt, homo, homo, -mv_env%w_fac, BabZ_3d(:, :, k), virt, &
                             mv_env%B_bar_ij(:, :, p), homo, 1.0_dp, RA_3d(:, :, k), virt)
               END DO
               ! B: R_ia,k -= w sum_jb \bar{B}^P_ib B^P_ja Z_jb,k, BjaZ_ab = sum_j B^P_ja Z_jb,k first
               IF (do_B) THEN
                  DO k = 1, nb
                     CALL DGEMM('N', 'T', virt, virt, homo, 1.0_dp, mv_env%B_ia(:, :, p), virt, &
                                Z_3d(:, :, k), virt, 0.0_dp, BjaZ_ab, virt)
                     CALL DGEMM('N', 'N', virt, homo, virt, -mv_env%w_fac, BjaZ_ab, virt, &
                                mv_env%B_bar_ia(:, :, p), virt, 1.0_dp, RB_3d(:, :, k), virt)
                  END DO
               END IF
            END DO
         END IF
         CALL timestop(handle2)

         ! the sum over the ranks lands on the owner of each row, never replicated
         CALL timeset(routineN//"_reduce", handle2)
         ! (A Z)_ia,k = (ε_a-ε_i) Z_ia,k + R^A_ia,k and (B Z)_ia,k = R^B_ia,k on the local rows
         DO k = 1, nb
            CALL mv_env%para_env%sum_scatter(RA(:, k:k), R_loc, mv_env%row_count)
            DO iloc = 1, nrow_local
               ia = row_indices(iloc)
               fm_AZ%local_data(iloc, col + k - 1) = R_loc(iloc) + mv_env%eps_diff(ia)*Z(ia, k)
            END DO
            IF (do_B) THEN
               CALL mv_env%para_env%sum_scatter(RB(:, k:k), R_loc, mv_env%row_count)
               fm_BZ%local_data(1:nrow_local, col + col_shift_B + k - 1) = R_loc(1:nrow_local)
            END IF
         END DO
         CALL timestop(handle2)
      END DO

      DEALLOCATE (Z_buf, RA_buf, BabZ_buf, T_Pk)
      IF (do_B) THEN
         DEALLOCATE (RB_buf, BjaZ_ab)
      END IF

      CALL timestop(handle)

   END SUBROUTINE bse_matvec_apply

! **************************************************************************************************
!> \brief Diagonal used by the Davidson correction, either ε_a-ε_i or the full diagonal
!>        A_ia,ia = ε_a-ε_i + α sum_P (B^P_ia)^2 - w sum_P \bar{B}^P_ii B^P_aa
!>        with α, w and \bar{B} as in bse_matvec_apply
!> \param mv_env slabs, prefactors and transition energies
!> \param precond_kind bse_precond_full_diag or bse_precond_qp_diff
!> \param diag replicated, size N_ov
! **************************************************************************************************
   SUBROUTINE bse_matvec_diagonal(mv_env, precond_kind, diag)

      TYPE(bse_matvec_env_type), INTENT(IN)              :: mv_env
      INTEGER, INTENT(IN)                                :: precond_kind
      REAL(KIND=dp), DIMENSION(:), INTENT(OUT)           :: diag

      CHARACTER(LEN=*), PARAMETER :: routineN = 'bse_matvec_diagonal'

      INTEGER                                            :: a, handle, i, ia, p

      CALL timeset(routineN, handle)

      diag(:) = 0.0_dp
      IF (precond_kind == bse_precond_full_diag) THEN
         DO p = 1, mv_env%n_ri_loc
            DO i = 1, mv_env%homo
               DO a = 1, mv_env%virt
                  ia = (i - 1)*mv_env%virt + a
                  diag(ia) = diag(ia) + mv_env%alpha*mv_env%B_ia(a, i, p)**2
                  IF (mv_env%w_fac /= 0.0_dp) THEN
                     diag(ia) = diag(ia) - mv_env%w_fac*mv_env%B_bar_ij(i, i, p)*mv_env%B_ab(a, a, p)
                  END IF
               END DO
            END DO
         END DO
         CALL mv_env%para_env%sum(diag)
      END IF
      diag(:) = diag(:) + mv_env%eps_diff(:)

      CALL timestop(handle)

   END SUBROUTINE bse_matvec_diagonal

! **************************************************************************************************
!> \brief Exact A (and B) on a list of transitions, replicated on every rank,
!>        A_kl = δ_kl (ε_a-ε_i) + α sum_P B^P_ia B^P_jb - w sum_P \bar{B}^P_ij B^P_ab
!>        B_kl =                  α sum_P B^P_ia B^P_jb - w sum_P \bar{B}^P_ib B^P_ja
!>        with k = (i,a), l = (j,b) and α, w, \bar{B} as in bse_matvec_apply
!> \param mv_env slabs, prefactors and transition energies
!> \param ia_list global transition indices ia = (i-1)*virt + a of the block
!> \param A_sub SIZE(ia_list) x SIZE(ia_list)
!> \param B_sub same, only formed when present
! **************************************************************************************************
   SUBROUTINE bse_matvec_subblock(mv_env, ia_list, A_sub, B_sub)

      TYPE(bse_matvec_env_type), INTENT(IN)              :: mv_env
      INTEGER, DIMENSION(:), INTENT(IN)                  :: ia_list
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT)        :: A_sub
      REAL(KIND=dp), DIMENSION(:, :), INTENT(OUT), &
         OPTIONAL                                        :: B_sub

      CHARACTER(LEN=*), PARAMETER :: routineN = 'bse_matvec_subblock'

      INTEGER                                            :: handle, k, l, n_sub, p
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: a_virt, i_occ
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: B_Pk

      CALL timeset(routineN, handle)

      n_sub = SIZE(ia_list)
      ALLOCATE (i_occ(n_sub), a_virt(n_sub), B_Pk(MAX(mv_env%n_ri_loc, 1), n_sub))
      DO k = 1, n_sub
         i_occ(k) = (ia_list(k) - 1)/mv_env%virt + 1
         a_virt(k) = ia_list(k) - (i_occ(k) - 1)*mv_env%virt
      END DO

      ! exchange, identical in A and B: α sum_P B_Pk B_Pl with B_Pk = B^P_ia at ia = ia_list(k)
      A_sub(:, :) = 0.0_dp
      IF (mv_env%alpha /= 0.0_dp .AND. mv_env%n_ri_loc > 0) THEN
         DO k = 1, n_sub
            B_Pk(:, k) = mv_env%B_ia(a_virt(k), i_occ(k), :)
         END DO
         CALL DGEMM('T', 'N', n_sub, n_sub, mv_env%n_ri_loc, mv_env%alpha, B_Pk, SIZE(B_Pk, 1), B_Pk, SIZE(B_Pk, 1), &
                    0.0_dp, A_sub, n_sub)
      END IF
      IF (PRESENT(B_sub)) B_sub(:, :) = A_sub(:, :)

      IF (mv_env%w_fac /= 0.0_dp) THEN
         DO p = 1, mv_env%n_ri_loc
            DO l = 1, n_sub
               DO k = 1, n_sub
                A_sub(k, l) = A_sub(k, l) - mv_env%w_fac*mv_env%B_bar_ij(i_occ(l), i_occ(k), p)*mv_env%B_ab(a_virt(l), a_virt(k), p)
               END DO
            END DO
            IF (PRESENT(B_sub)) THEN
               DO l = 1, n_sub
                  DO k = 1, n_sub
                B_sub(k, l) = B_sub(k, l) - mv_env%w_fac*mv_env%B_ia(a_virt(k), i_occ(l), p)*mv_env%B_bar_ia(a_virt(l), i_occ(k), p)
                  END DO
               END DO
            END IF
         END DO
      END IF

      CALL mv_env%para_env%sum(A_sub)
      IF (PRESENT(B_sub)) CALL mv_env%para_env%sum(B_sub)
      DO k = 1, n_sub
         A_sub(k, k) = A_sub(k, k) + mv_env%eps_diff(ia_list(k))
      END DO

      DEALLOCATE (i_occ, a_virt, B_Pk)

      CALL timestop(handle)

   END SUBROUTINE bse_matvec_subblock

! **************************************************************************************************
!> \brief Debug check of the matrix-free application against the explicit matrices A (and B),
!>        dev = max_ia,k |(A Z)_ia,k - sum_jb A_ia,jb Z_jb,k| over up to eight unit vectors Z_ia,k = δ_ia,k
!>        and one dense vector Z_ia = sin(ia), the same for B, and max_ia |d_ia - A_ia,ia| for the diagonal
!> \param mv_env the environment under test
!> \param fm_A_explicit A from create_A_and_B, N_ov x N_ov on the grid of the caller
!> \param unit_nr output unit, positive on the writing rank only
!> \param fm_B_explicit B, present for an ABBA run
! **************************************************************************************************
   SUBROUTINE bse_matvec_selfcheck(mv_env, fm_A_explicit, unit_nr, fm_B_explicit)

      TYPE(bse_matvec_env_type), INTENT(IN)              :: mv_env
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_A_explicit
      INTEGER, INTENT(IN)                                :: unit_nr
      TYPE(cp_fm_type), INTENT(IN), OPTIONAL             :: fm_B_explicit

      CHARACTER(LEN=*), PARAMETER :: routineN = 'bse_matvec_selfcheck'

      INTEGER                                            :: handle, ia, iloc, k, n_ov, n_unit, ncol, &
                                                            nrow_local
      INTEGER, DIMENSION(:), POINTER                     :: row_indices
      REAL(KIND=dp)                                      :: dev_A, dev_B, dev_diag
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: diag
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: mv_result, ref_matrix, ref_result, Z
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: fm_AZ, fm_BZ, fm_Z

      CALL timeset(routineN, handle)

      n_ov = mv_env%n_ov
      ! eight unit vectors read eight columns of the explicit matrix exactly; the dense column covers the whole contraction
      n_unit = MIN(n_ov, 8)
      dev_B = 0.0_dp
      ncol = n_unit + 1

      ! trial vectors: e_k for k <= n_unit, then one dense column Z_ia = sin(ia)
      ALLOCATE (Z(n_ov, ncol))
      Z(:, :) = 0.0_dp
      DO k = 1, n_unit
         Z(k, k) = 1.0_dp
      END DO
      DO ia = 1, n_ov
         Z(ia, ncol) = SIN(REAL(ia, dp))
      END DO

      CALL bse_matvec_vector_struct(mv_env, ncol, fm_struct)
      CALL cp_fm_create(fm_Z, fm_struct, name="fm_Z_selfcheck")
      CALL cp_fm_create(fm_AZ, fm_struct, name="fm_AZ_selfcheck")
      CALL cp_fm_set_all(fm_AZ, 0.0_dp)
      IF (PRESENT(fm_B_explicit)) THEN
         CALL cp_fm_create(fm_BZ, fm_struct, name="fm_BZ_selfcheck")
         CALL cp_fm_set_all(fm_BZ, 0.0_dp)
      END IF
      CALL cp_fm_struct_release(fm_struct)

      CALL cp_fm_get_info(fm_Z, nrow_local=nrow_local, row_indices=row_indices)
      DO k = 1, ncol
         DO iloc = 1, nrow_local
            fm_Z%local_data(iloc, k) = Z(row_indices(iloc), k)
         END DO
      END DO

      IF (PRESENT(fm_B_explicit)) THEN
         CALL bse_matvec_apply(mv_env, fm_Z, 1, ncol, fm_AZ, fm_BZ)
      ELSE
         CALL bse_matvec_apply(mv_env, fm_Z, 1, ncol, fm_AZ)
      END IF

      ALLOCATE (ref_matrix(n_ov, n_ov), mv_result(n_ov, ncol), ref_result(n_ov, ncol), diag(n_ov))

      CALL cp_fm_get_submatrix(fm_A_explicit, ref_matrix)
      CALL cp_fm_get_submatrix(fm_AZ, mv_result)
      ref_result(:, :) = MATMUL(ref_matrix, Z)
      dev_A = MAXVAL(ABS(mv_result - ref_result))
      CALL bse_matvec_diagonal(mv_env, bse_precond_full_diag, diag)
      dev_diag = 0.0_dp
      DO ia = 1, n_ov
         dev_diag = MAX(dev_diag, ABS(diag(ia) - ref_matrix(ia, ia)))
      END DO

      IF (PRESENT(fm_B_explicit)) THEN
         CALL cp_fm_get_submatrix(fm_B_explicit, ref_matrix)
         CALL cp_fm_get_submatrix(fm_BZ, mv_result)
         ref_result(:, :) = MATMUL(ref_matrix, Z)
         dev_B = MAXVAL(ABS(mv_result - ref_result))
      END IF

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, '(T2,A10,T13,A,T59,ES22.6)') 'BSE|DEBUG|', &
            'Max deviation matvec vs explicit A', dev_A
         IF (PRESENT(fm_B_explicit)) THEN
            WRITE (unit_nr, '(T2,A10,T13,A,T59,ES22.6)') 'BSE|DEBUG|', &
               'Max deviation matvec vs explicit B', dev_B
         END IF
         WRITE (unit_nr, '(T2,A10,T13,A,T59,ES22.6)') 'BSE|DEBUG|', &
            'Max deviation diagonal vs explicit A', dev_diag
      END IF

      DEALLOCATE (Z, ref_matrix, mv_result, ref_result, diag)
      CALL cp_fm_release(fm_Z)
      CALL cp_fm_release(fm_AZ)
      IF (PRESENT(fm_B_explicit)) CALL cp_fm_release(fm_BZ)

      CALL timestop(handle)

   END SUBROUTINE bse_matvec_selfcheck

END MODULE bse_matvec
