!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routine to return block diagonal density matrix. Blocks correspond to the atomic densities
!> \par History
!>       2006.03 Moved here from qs_scf.F [Joost VandeVondele]
!>       2022.05 split from qs_initial_guess.F to break circular dependency [Harald Forbert]
! **************************************************************************************************
MODULE qs_atomic_block
   USE atom_kind_orbitals,              ONLY: calculate_atomic_orbitals
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_add_on_diag, dbcsr_dot, dbcsr_get_info, dbcsr_iterator_blocks_left, &
        dbcsr_iterator_next_block, dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, &
        dbcsr_p_type, dbcsr_scale, dbcsr_set, dbcsr_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE qs_kind_types,                   ONLY: qs_kind_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_atomic_block'

   PUBLIC ::  calculate_atomic_block_dm

   TYPE atom_matrix_type
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER   :: mat
   END TYPE atom_matrix_type

CONTAINS

! **************************************************************************************************
!> \brief returns a block diagonal density matrix. Blocks correspond to the atomic densities.
!> \param pmatrix ...
!> \param matrix_s ...
!> \param atomic_kind_set ...
!> \param qs_kind_set ...
!> \param nspin ...
!> \param nelectron_spin ...
!> \param ounit ...
!> \param para_env ...
! **************************************************************************************************
   SUBROUTINE calculate_atomic_block_dm(pmatrix, matrix_s, atomic_kind_set, &
                                        qs_kind_set, nspin, nelectron_spin, ounit, para_env)
      TYPE(dbcsr_p_type), DIMENSION(:), INTENT(INOUT)    :: pmatrix
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_s
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      INTEGER, INTENT(IN)                                :: nspin
      INTEGER, DIMENSION(:), INTENT(IN)                  :: nelectron_spin
      INTEGER, INTENT(IN)                                :: ounit
      TYPE(mp_para_env_type)                             :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'calculate_atomic_block_dm'

      INTEGER                                            :: blk, handle, icol, ikind, irow, ispin, &
                                                            nc, nkind, nocc(2)
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: nok
      REAL(dp), DIMENSION(:, :), POINTER                 :: pdata
      REAL(KIND=dp)                                      :: rds, rscale, trps1
      TYPE(atom_matrix_type), ALLOCATABLE, DIMENSION(:)  :: pmat
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_type), POINTER                          :: matrix_p
      TYPE(qs_kind_type), POINTER                        :: qs_kind

      CALL timeset(routineN, handle)

      nkind = SIZE(atomic_kind_set)
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, kind_of=kind_of)
      ALLOCATE (pmat(nkind))
      ALLOCATE (nok(2, nkind))

      ! precompute the atomic blocks corresponding to spherical atoms
      DO ikind = 1, nkind
         atomic_kind => atomic_kind_set(ikind)
         qs_kind => qs_kind_set(ikind)
         NULLIFY (pmat(ikind)%mat)
         IF (ounit > 0) THEN
            WRITE (UNIT=ounit, FMT="(/,T2,A)") &
               "Guess for atomic kind: "//TRIM(atomic_kind%name)
         END IF
         CALL calculate_atomic_orbitals(atomic_kind, qs_kind, iunit=ounit, &
                                        pmat=pmat(ikind)%mat, nocc=nocc)
         nok(1:2, ikind) = nocc(1:2)
      END DO

      rscale = 1.0_dp
      IF (nspin == 2) rscale = 0.5_dp

      DO ispin = 1, nspin
         IF ((ounit > 0) .AND. (nspin > 1)) THEN
            WRITE (UNIT=ounit, FMT="(/,T2,A,I0)") "Spin ", ispin
         END IF

         matrix_p => pmatrix(ispin)%matrix
         CALL dbcsr_set(matrix_p, 0.0_dp)

         nocc(ispin) = 0
         CALL dbcsr_iterator_start(iter, matrix_p)
         DO WHILE (dbcsr_iterator_blocks_left(iter))
            CALL dbcsr_iterator_next_block(iter, irow, icol, pdata, blk)
            ikind = kind_of(irow)
            IF (icol .EQ. irow) THEN
               IF (ispin == 1) THEN
                  pdata(:, :) = pmat(ikind)%mat(:, :, 1)*rscale + &
                                pmat(ikind)%mat(:, :, 2)*rscale
               ELSE
                  pdata(:, :) = pmat(ikind)%mat(:, :, 1)*rscale - &
                                pmat(ikind)%mat(:, :, 2)*rscale
               END IF
               nocc(ispin) = nocc(ispin) + nok(ispin, ikind)
            END IF
         END DO
         CALL dbcsr_iterator_stop(iter)

         CALL dbcsr_dot(matrix_p, matrix_s, trps1)
         rds = 0.0_dp
         ! could be a ghost-atoms-only simulation
         IF (nelectron_spin(ispin) > 0) THEN
            rds = REAL(nelectron_spin(ispin), dp)/trps1
         END IF
         CALL dbcsr_scale(matrix_p, rds)

         IF (ounit > 0) THEN
            IF (nspin > 1) THEN
               WRITE (UNIT=ounit, FMT="(T2,A,I1)") &
                  "Re-scaling the density matrix to get the right number of electrons for spin ", ispin
            ELSE
               WRITE (UNIT=ounit, FMT="(T2,A)") &
                  "Re-scaling the density matrix to get the right number of electrons"
            END IF
            WRITE (ounit, '(T19,A,T44,A,T67,A)') "# Electrons", "Trace(P)", "Scaling factor"
            WRITE (ounit, '(T20,I10,T40,F12.3,T67,F14.3)') nelectron_spin(ispin), trps1, rds
         END IF

         IF (nspin > 1) THEN
            CALL para_env%sum(nocc)
            IF (nelectron_spin(ispin) > nocc(ispin)) THEN
               rds = 0.99_dp
               CALL dbcsr_scale(matrix_p, rds)
               rds = (1.0_dp - rds)*nelectron_spin(ispin)
               CALL dbcsr_get_info(matrix_p, nfullcols_total=nc)
               rds = rds/REAL(nc, KIND=dp)
               CALL dbcsr_add_on_diag(matrix_p, rds)
               IF (ounit > 0) THEN
                  WRITE (UNIT=ounit, FMT="(T4,A,/,T4,A,T59,F20.12)") &
                     "More MOs than initial guess orbitals detected", &
                     "Add constant to diagonal elements ", rds
               END IF
            END IF
         END IF

      END DO

      DO ikind = 1, nkind
         IF (ASSOCIATED(pmat(ikind)%mat)) THEN
            DEALLOCATE (pmat(ikind)%mat)
         END IF
      END DO
      DEALLOCATE (pmat)

      DEALLOCATE (kind_of, nok)

      CALL timestop(handle)

   END SUBROUTINE calculate_atomic_block_dm

END MODULE qs_atomic_block
