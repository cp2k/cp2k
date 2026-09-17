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
             bse_matvec_diagonal, bse_matvec_vector_struct, bse_matvec_selfcheck

! **************************************************************************************************
!> \brief RI slabs of one spin channel, every rank holding n_ri_loc whole RI slices
!> \param alpha prefactor of the exchange term (2 singlet, 0 triplet)
!> \param w_fac prefactor of the screened term (0 switches it off)
!> \param eps_diff ε_a - ε_i at compound index ia = (i-1)*virt + a
!> \param B_ia B^P_ia at (a, i, P)
!> \param Bbar_ia \bar{B}^P_ia at (a, i, P), only for do_abba
!> \param Bbar_ij \bar{B}^P_ij at (j, i, P)
!> \param B_ab B^P_ab at (b, a, P)
!> \param blacs_env npe x 1 process grid of the slabs and of all trial vectors
! **************************************************************************************************
   TYPE bse_matvec_env_type
      INTEGER                                            :: homo = 0, virt = 0, n_ov = 0, &
                                                            n_ri_loc = 0
      REAL(KIND=dp)                                      :: alpha = 0.0_dp, w_fac = 0.0_dp
      LOGICAL                                            :: do_abba = .FALSE.
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eps_diff
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: B_ia, Bbar_ia, Bbar_ij, B_ab
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env => NULL()
      TYPE(mp_para_env_type), POINTER                    :: para_env => NULL()
   END TYPE bse_matvec_env_type

CONTAINS

! **************************************************************************************************
!> \brief Moves the RI slabs onto an npe x 1 process grid, such that every rank owns whole RI
!>        slices, and stores them as local 3-index arrays
!> \param mv_env ...
!> \param fm_S_ia B^P_ia, N_RI x homo*virt
!> \param fm_Sw_ij slab that enters W_ij,ab as \bar{B}^P_ij, N_RI x homo*homo
!> \param fm_S_ab B^P_ab, N_RI x virt*virt
!> \param fm_Sw_ia slab that enters W_ib,aj as \bar{B}^P_ia, N_RI x homo*virt
!> \param eps_reduced quasiparticle energies of the homo+virt active levels
!> \param homo ...
!> \param virt ...
!> \param alpha ...
!> \param w_fac ...
!> \param do_abba ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE bse_matvec_create(mv_env, fm_S_ia, fm_Sw_ij, fm_S_ab, fm_Sw_ia, eps_reduced, &
                                homo, virt, alpha, w_fac, do_abba, unit_nr)

      TYPE(bse_matvec_env_type), INTENT(OUT)             :: mv_env
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_S_ia, fm_Sw_ij, fm_S_ab, fm_Sw_ia
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: eps_reduced
      INTEGER, INTENT(IN)                                :: homo, virt
      REAL(KIND=dp), INTENT(IN)                          :: alpha, w_fac
      LOGICAL, INTENT(IN)                                :: do_abba
      INTEGER, INTENT(IN)                                :: unit_nr

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'bse_matvec_create'

      CHARACTER(LEN=12)                                  :: n_idle_str, n_ri_str, npe_str
      INTEGER                                            :: a, handle, i, n_ri
      REAL(KIND=dp)                                      :: mem_slabs_GB

      CALL timeset(routineN, handle)

      mv_env%homo = homo
      mv_env%virt = virt
      mv_env%n_ov = homo*virt
      mv_env%alpha = alpha
      mv_env%w_fac = w_fac
      mv_env%do_abba = do_abba
      mv_env%para_env => fm_S_ia%matrix_struct%para_env

      CALL cp_blacs_env_create(mv_env%blacs_env, mv_env%para_env, &
                               grid_2d=[mv_env%para_env%num_pe, 1])

      ALLOCATE (mv_env%eps_diff(mv_env%n_ov))
      DO i = 1, homo
         DO a = 1, virt
            mv_env%eps_diff((i - 1)*virt + a) = eps_reduced(homo + a) - eps_reduced(i)
         END DO
      END DO

      CALL slice_slab(mv_env, fm_S_ia, virt, homo, mv_env%B_ia)
      IF (w_fac /= 0.0_dp) THEN
         CALL slice_slab(mv_env, fm_Sw_ij, homo, homo, mv_env%Bbar_ij)
         CALL slice_slab(mv_env, fm_S_ab, virt, virt, mv_env%B_ab)
         IF (do_abba) CALL slice_slab(mv_env, fm_Sw_ia, virt, homo, mv_env%Bbar_ia)
      END IF
      mv_env%n_ri_loc = SIZE(mv_env%B_ia, 3)

      CALL cp_fm_get_info(fm_S_ia, nrow_global=n_ri)
      IF (mv_env%para_env%num_pe > n_ri) THEN
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
!> \brief Redistributes one N_RI x d1*d2 slab to whole RI rows per rank and copies the owned
!>        rows into a contiguous (d1, d2, n_ri_loc) array
!> \param mv_env ...
!> \param fm_slab ...
!> \param d1 fast index of the pair index of the slab
!> \param d2 slow index of the pair index of the slab
!> \param slab_loc ...
! **************************************************************************************************
   SUBROUTINE slice_slab(mv_env, fm_slab, d1, d2, slab_loc)

      TYPE(bse_matvec_env_type), INTENT(IN)              :: mv_env
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_slab
      INTEGER, INTENT(IN)                                :: d1, d2
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :), &
         INTENT(OUT)                                     :: slab_loc

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'slice_slab'

      INTEGER                                            :: handle, n_ri, n_ri_loc, npair, p
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: fm_sliced

      CALL timeset(routineN, handle)

      ! row count from the struct: the slabs may carry fewer rows than the full RI basis
      CALL cp_fm_get_info(fm_slab, nrow_global=n_ri, ncol_global=npair)
      CPASSERT(npair == d1*d2)

      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, para_env=mv_env%para_env, context=mv_env%blacs_env, &
                               nrow_global=n_ri, ncol_global=npair, nrow_block=1, &
                               force_block=.TRUE.)
      CALL cp_fm_create(fm_sliced, fm_struct, name="fm_slab_ri_sliced")
      CALL cp_fm_struct_release(fm_struct)
      CALL cp_fm_to_fm_submat_general(fm_slab, fm_sliced, n_ri, npair, 1, 1, 1, 1, &
                                      fm_slab%matrix_struct%context)

      CALL cp_fm_get_info(fm_sliced, nrow_local=n_ri_loc)
      ALLOCATE (slab_loc(d1, d2, n_ri_loc))
      DO p = 1, n_ri_loc
         slab_loc(:, :, p) = RESHAPE(fm_sliced%local_data(p, 1:npair), [d1, d2])
      END DO
      CALL cp_fm_release(fm_sliced)

      CALL timestop(handle)

   END SUBROUTINE slice_slab

! **************************************************************************************************
!> \brief ...
!> \param mv_env ...
! **************************************************************************************************
   SUBROUTINE bse_matvec_release(mv_env)

      TYPE(bse_matvec_env_type), INTENT(INOUT)           :: mv_env

      IF (ALLOCATED(mv_env%eps_diff)) DEALLOCATE (mv_env%eps_diff)
      IF (ALLOCATED(mv_env%B_ia)) DEALLOCATE (mv_env%B_ia)
      IF (ALLOCATED(mv_env%Bbar_ia)) DEALLOCATE (mv_env%Bbar_ia)
      IF (ALLOCATED(mv_env%Bbar_ij)) DEALLOCATE (mv_env%Bbar_ij)
      IF (ALLOCATED(mv_env%B_ab)) DEALLOCATE (mv_env%B_ab)
      IF (ASSOCIATED(mv_env%blacs_env)) CALL cp_blacs_env_release(mv_env%blacs_env)
      NULLIFY (mv_env%para_env)

   END SUBROUTINE bse_matvec_release

! **************************************************************************************************
!> \brief Matrix structure of a block of trial vectors: rows ia distributed, all columns local
!> \param mv_env ...
!> \param ncol_global ...
!> \param fm_struct ...
! **************************************************************************************************
   SUBROUTINE bse_matvec_vector_struct(mv_env, ncol_global, fm_struct)

      TYPE(bse_matvec_env_type), INTENT(IN)              :: mv_env
      INTEGER, INTENT(IN)                                :: ncol_global
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct

      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, para_env=mv_env%para_env, context=mv_env%blacs_env, &
                               nrow_global=mv_env%n_ov, ncol_global=ncol_global)

   END SUBROUTINE bse_matvec_vector_struct

! **************************************************************************************************
!> \brief Applies A (and B) to ncol trial vectors without forming an N_ov x N_ov object,
!>        (A Z)_ia = (ε_a-ε_i) Z_ia + α sum_P B^P_ia sum_jb B^P_jb Z_jb - w sum_Pjb \bar{B}^P_ij B^P_ab Z_jb
!>        (B Z)_ia =                  α sum_P B^P_ia sum_jb B^P_jb Z_jb - w sum_Pjb \bar{B}^P_ib B^P_ja Z_jb
!>        α: exchange prefactor (2 singlet, 0 triplet), w: prefactor of the screened term,
!>        \bar{B}: slab bound by the caller of bse_matvec_create (B itself for TDHF and ALPHA
!>        screening, sum_Q [1+Q(0)]^-1_PQ B^Q otherwise); the sums over P run over the local RI
!>        slices and are completed by a sum over the ranks
!> \param mv_env ...
!> \param fm_Z trial vectors, columns first_col .. first_col+ncol-1 are read
!> \param first_col ...
!> \param ncol ...
!> \param fm_AZ receives A Z in the same columns
!> \param fm_BZ receives B Z in the same columns
! **************************************************************************************************
   SUBROUTINE bse_matvec_apply(mv_env, fm_Z, first_col, ncol, fm_AZ, fm_BZ)

      TYPE(bse_matvec_env_type), INTENT(IN), TARGET      :: mv_env
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_Z
      INTEGER, INTENT(IN)                                :: first_col, ncol
      TYPE(cp_fm_type), INTENT(IN)                       :: fm_AZ
      TYPE(cp_fm_type), INTENT(IN), OPTIONAL             :: fm_BZ

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'bse_matvec_apply'
      REAL(KIND=dp), PARAMETER                           :: mem_fraction = 0.25_dp

      INTEGER                                            :: col, handle, handle2, homo, ia, iloc, k, &
                                                            n_arrays, n_ov, n_ri_loc, nb, nb_max, &
                                                            nrow_local, p, virt
      INTEGER, DIMENSION(:), POINTER                     :: row_indices
      LOGICAL                                            :: do_B
      REAL(KIND=dp)                                      :: mem_avail_GB
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:), TARGET   :: RA_buf, RB_buf, tmp_buf, Z_buf
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: T_Pk, tmp_ab
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :), &
         POINTER                                         :: B_ia_2d, RA, RB, tmp_wide, Z, Z_wide
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :, :), &
         POINTER                                         :: RA_3d, RB_3d, tmp_3d, Z_3d

      CALL timeset(routineN, handle)

      homo = mv_env%homo
      virt = mv_env%virt
      n_ov = mv_env%n_ov
      n_ri_loc = mv_env%n_ri_loc
      do_B = PRESENT(fm_BZ)
      IF (do_B) THEN
         CPASSERT(mv_env%do_abba)
      END IF

      CALL cp_fm_get_info(fm_Z, nrow_local=nrow_local, row_indices=row_indices)

      ! columns per pass from the memory of the replicated work arrays Z, tmp, R^A (and R^B)
      n_arrays = 3
      IF (do_B) n_arrays = 4
      CALL mp_mem_avail_per_rank_GB(mv_env%para_env, mem_avail_GB)
      nb_max = ncol
      IF (mem_avail_GB > 0.0_dp) THEN
         nb_max = INT(MIN(REAL(ncol, dp), &
                          mem_fraction*mem_avail_GB*1.0E9_dp/(8.0_dp*REAL(n_arrays, dp)*REAL(n_ov, dp))))
         nb_max = MAX(nb_max, 1)
      END IF

      ALLOCATE (Z_buf(n_ov*nb_max), RA_buf(n_ov*nb_max), tmp_buf(n_ov*nb_max))
      IF (do_B) THEN
         ALLOCATE (RB_buf(n_ov*nb_max), tmp_ab(virt, virt))
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
         tmp_3d(1:virt, 1:homo, 1:nb) => tmp_buf(1:n_ov*nb)
         tmp_wide(1:virt, 1:homo*nb) => tmp_buf(1:n_ov*nb)
         IF (do_B) THEN
            RB(1:n_ov, 1:nb) => RB_buf(1:n_ov*nb)
            RB_3d(1:virt, 1:homo, 1:nb) => RB_buf(1:n_ov*nb)
         END IF

         ! Z_ia,k on every rank
         CALL timeset(routineN//"_gather", handle2)
         Z(:, :) = 0.0_dp
         DO k = 1, nb
            DO iloc = 1, nrow_local
               Z(row_indices(iloc), k) = fm_Z%local_data(iloc, col + k - 1)
            END DO
         END DO
         CALL mv_env%para_env%sum(Z)
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
               ! A: R_ia,k -= w sum_jb \bar{B}^P_ij B^P_ab Z_jb,k
               CALL DGEMM('T', 'N', virt, homo*nb, virt, 1.0_dp, mv_env%B_ab(:, :, p), virt, &
                          Z_wide, virt, 0.0_dp, tmp_wide, virt)
               DO k = 1, nb
                  CALL DGEMM('N', 'N', virt, homo, homo, -mv_env%w_fac, tmp_3d(:, :, k), virt, &
                             mv_env%Bbar_ij(:, :, p), homo, 1.0_dp, RA_3d(:, :, k), virt)
               END DO
               ! B: R_ia,k -= w sum_jb \bar{B}^P_ib B^P_ja Z_jb,k
               IF (do_B) THEN
                  DO k = 1, nb
                     CALL DGEMM('N', 'T', virt, virt, homo, 1.0_dp, mv_env%B_ia(:, :, p), virt, &
                                Z_3d(:, :, k), virt, 0.0_dp, tmp_ab, virt)
                     CALL DGEMM('N', 'N', virt, homo, virt, -mv_env%w_fac, tmp_ab, virt, &
                                mv_env%Bbar_ia(:, :, p), virt, 1.0_dp, RB_3d(:, :, k), virt)
                  END DO
               END IF
            END DO
         END IF
         CALL timestop(handle2)

         CALL timeset(routineN//"_reduce", handle2)
         CALL mv_env%para_env%sum(RA)
         IF (do_B) CALL mv_env%para_env%sum(RB)
         DO k = 1, nb
            DO iloc = 1, nrow_local
               ia = row_indices(iloc)
               fm_AZ%local_data(iloc, col + k - 1) = RA(ia, k) + mv_env%eps_diff(ia)*Z(ia, k)
               IF (do_B) fm_BZ%local_data(iloc, col + k - 1) = RB(ia, k)
            END DO
         END DO
         CALL timestop(handle2)
      END DO

      DEALLOCATE (Z_buf, RA_buf, tmp_buf, T_Pk)
      IF (do_B) THEN
         DEALLOCATE (RB_buf, tmp_ab)
      END IF

      CALL timestop(handle)

   END SUBROUTINE bse_matvec_apply

! **************************************************************************************************
!> \brief Diagonal used by the Davidson correction, either ε_a-ε_i or the full diagonal
!>        A_ia,ia = ε_a-ε_i + α sum_P (B^P_ia)^2 - w sum_P \bar{B}^P_ii B^P_aa
!>        with α, w and \bar{B} as in bse_matvec_apply
!> \param mv_env ...
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
                     diag(ia) = diag(ia) - mv_env%w_fac*mv_env%Bbar_ij(i, i, p)*mv_env%B_ab(a, a, p)
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
!> \brief Debug check of the matrix-free application against the explicit matrices A (and B):
!>        up to eight unit vectors, one dense vector and the diagonal
!> \param mv_env ...
!> \param fm_A_explicit ...
!> \param unit_nr ...
!> \param fm_B_explicit ...
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
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: M_explicit, MZ, MZ_ref, Z
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type)                                   :: fm_AZ, fm_BZ, fm_Z

      CALL timeset(routineN, handle)

      n_ov = mv_env%n_ov
      n_unit = MIN(n_ov, 8)
      dev_B = 0.0_dp
      ncol = n_unit + 1

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

      ALLOCATE (M_explicit(n_ov, n_ov), MZ(n_ov, ncol), MZ_ref(n_ov, ncol), diag(n_ov))

      CALL cp_fm_get_submatrix(fm_A_explicit, M_explicit)
      CALL cp_fm_get_submatrix(fm_AZ, MZ)
      MZ_ref(:, :) = MATMUL(M_explicit, Z)
      dev_A = MAXVAL(ABS(MZ - MZ_ref))
      CALL bse_matvec_diagonal(mv_env, bse_precond_full_diag, diag)
      dev_diag = 0.0_dp
      DO ia = 1, n_ov
         dev_diag = MAX(dev_diag, ABS(diag(ia) - M_explicit(ia, ia)))
      END DO

      IF (PRESENT(fm_B_explicit)) THEN
         CALL cp_fm_get_submatrix(fm_B_explicit, M_explicit)
         CALL cp_fm_get_submatrix(fm_BZ, MZ)
         MZ_ref(:, :) = MATMUL(M_explicit, Z)
         dev_B = MAXVAL(ABS(MZ - MZ_ref))
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

      DEALLOCATE (Z, M_explicit, MZ, MZ_ref, diag)
      CALL cp_fm_release(fm_Z)
      CALL cp_fm_release(fm_AZ)
      IF (PRESENT(fm_B_explicit)) CALL cp_fm_release(fm_BZ)

      CALL timestop(handle)

   END SUBROUTINE bse_matvec_selfcheck

END MODULE bse_matvec
