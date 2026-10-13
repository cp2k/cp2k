!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Block-addressable, immutable compressed storage for RI-HFX tensors.
!> Each block starts on a word boundary. This avoids shared decoding cursors and the
!> per-bit-width caches of the four-center integral stream. Scratch space is bounded
!> by one tensor block and a small integer buffer per worker, not by the tensor size.
! **************************************************************************************************
MODULE hfx_ri_compression
   USE dbt_api,                         ONLY: &
        dbt_blk_sizes, dbt_clear, dbt_get_block, dbt_get_num_blocks, dbt_iterator_next_block, &
        dbt_iterator_num_blocks, dbt_iterator_start, dbt_iterator_stop, dbt_iterator_type, &
        dbt_put_block, dbt_reserve_blocks, dbt_type
   USE hfx_compression_core_methods,    ONLY: bits2ints_specific,&
                                              ints2bits_specific
   USE kinds,                           ONLY: dp,&
                                              int_8
#include "./base/base_uses.f90"

!$ USE OMP_LIB, ONLY: omp_get_num_threads, omp_get_thread_num

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: hfx_ri_compression_type, compress_ri_tensor, decompress_ri_tensor, clear_ri_compressed_tensor

   TYPE hfx_ri_compression_type
      INTEGER, ALLOCATABLE                              :: nbits(:)
      INTEGER(int_8), ALLOCATABLE                       :: offsets(:), DATA(:)
   END TYPE hfx_ri_compression_type

   INTEGER, PARAMETER                                   :: chunk_size = 1024

CONTAINS

! **************************************************************************************************
!> \brief Release compressed payload and metadata.
!> \param compressed ...
! **************************************************************************************************
   SUBROUTINE clear_ri_compressed_tensor(compressed)
      TYPE(hfx_ri_compression_type), INTENT(INOUT)       :: compressed

      IF (ALLOCATED(compressed%nbits)) DEALLOCATE (compressed%nbits)
      IF (ALLOCATED(compressed%offsets)) DEALLOCATE (compressed%offsets)
      IF (ALLOCATED(compressed%data)) DEALLOCATE (compressed%data)
   END SUBROUTINE clear_ri_compressed_tensor

! **************************************************************************************************
!> \brief Compress local blocks in parallel, then clear the input tensor.
!> \param tensor ...
!> \param blk_indices Block indices, in the same order as the compressed metadata.
!> \param compressed ...
!> \param eps Storage threshold; the existing RI-HFX quantization rule is preserved.
!> \param memory Accumulated payload and array storage in MiB (excluding allocator overhead).
! **************************************************************************************************
   SUBROUTINE compress_ri_tensor(tensor, blk_indices, compressed, eps, memory)
      TYPE(dbt_type), INTENT(INOUT)                      :: tensor
      INTEGER, ALLOCATABLE, INTENT(INOUT)                :: blk_indices(:, :)
      TYPE(hfx_ri_compression_type), INTENT(INOUT)       :: compressed
      REAL(dp), INTENT(IN)                               :: eps
      REAL(dp), INTENT(INOUT)                            :: memory

      CHARACTER(LEN=*), PARAMETER :: routineN = 'compress_ri_tensor'

      INTEGER                                            :: handle, i, iblk, nblk, nlocal, offset, &
                                                            shared_offset
      INTEGER(int_8)                                     :: nwords
      INTEGER, DIMENSION(3)                              :: ind
      LOGICAL                                            :: found
      REAL(dp), ALLOCATABLE, TARGET                      :: block(:, :, :)
      REAL(dp), POINTER                                  :: values(:)
      TYPE(dbt_iterator_type)                            :: iter

      CALL timeset(routineN, handle)
      CPASSERT(eps > 0.0_dp)
      CALL clear_ri_compressed_tensor(compressed)
      IF (ALLOCATED(blk_indices)) DEALLOCATE (blk_indices)
      nblk = dbt_get_num_blocks(tensor)
      ALLOCATE (blk_indices(nblk, 3), compressed%nbits(nblk), compressed%offsets(nblk + 1))
      shared_offset = 0

      ! The DBT iterator partitions shards between the actual OpenMP workers.
!$OMP PARALLEL DEFAULT(NONE) SHARED(tensor,blk_indices,compressed,eps,shared_offset) &
!$OMP PRIVATE(iter,nlocal,offset,iblk,i,ind,block,values,found)
      CALL dbt_iterator_start(iter, tensor)
      nlocal = dbt_iterator_num_blocks(iter)
!$OMP CRITICAL(hfx_ri_compression_offset)
      offset = shared_offset
      shared_offset = shared_offset + nlocal
!$OMP END CRITICAL(hfx_ri_compression_offset)
      DO iblk = 1, nlocal
         i = offset + iblk
         CALL dbt_iterator_next_block(iter, ind)
         blk_indices(i, :) = ind
         CALL dbt_get_block(tensor, ind, block, found)
         CPASSERT(found)
         values(1:SIZE(block)) => block
         compressed%nbits(i) = block_nbits(values, eps)
         compressed%offsets(i) = (INT(SIZE(block), int_8)*compressed%nbits(i) + 63_int_8)/64_int_8
         NULLIFY (values)
         DEALLOCATE (block)
      END DO
      CALL dbt_iterator_stop(iter)
!$OMP END PARALLEL
      CPASSERT(shared_offset == nblk)

      ! Only the small metadata prefix sum is serial; no tensor-sized integer copy.
      nwords = 0_int_8
      DO i = 1, nblk
         nwords = nwords + compressed%offsets(i)
         compressed%offsets(i) = nwords - compressed%offsets(i) + 1_int_8
      END DO
      compressed%offsets(nblk + 1) = nwords + 1_int_8
      ALLOCATE (compressed%data(nwords))

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(DYNAMIC,1) &
!$OMP SHARED(tensor,blk_indices,compressed,eps,nblk) PRIVATE(i,ind,block,values,found)
      DO i = 1, nblk
         IF (compressed%nbits(i) == 0) CYCLE
         ind = blk_indices(i, :)
         CALL dbt_get_block(tensor, ind, block, found)
         CPASSERT(found)
         values(1:SIZE(block)) => block
         CALL pack_block(values, eps, compressed%nbits(i), &
                         compressed%data(compressed%offsets(i):compressed%offsets(i + 1) - 1))
         NULLIFY (values)
         DEALLOCATE (block)
      END DO
!$OMP END PARALLEL DO

      CALL dbt_clear(tensor)
      memory = memory + REAL(nwords*8_int_8 + INT(nblk, int_8)*24_int_8 + 8_int_8, dp)/1024.0_dp**2
      CALL timestop(handle)
   END SUBROUTINE compress_ri_tensor

! **************************************************************************************************
!> \brief Decode blocks independently; the compressed input has no mutable cursor.
!> \param tensor Empty tensor receiving the decoded blocks.
!> \param blk_indices ...
!> \param compressed ...
!> \param eps ...
!> \param sparsity_only Reserve zero blocks without decoding their values.
! **************************************************************************************************
   SUBROUTINE decompress_ri_tensor(tensor, blk_indices, compressed, eps, sparsity_only)
      TYPE(dbt_type), INTENT(INOUT)                      :: tensor
      INTEGER, INTENT(IN)                                :: blk_indices(:, :)
      TYPE(hfx_ri_compression_type), INTENT(IN)          :: compressed
      REAL(dp), INTENT(IN)                               :: eps
      LOGICAL, INTENT(IN), OPTIONAL                      :: sparsity_only

      CHARACTER(LEN=*), PARAMETER :: routineN = 'decompress_ri_tensor'

      INTEGER                                            :: a, b, handle, i, me, nblk, nt
      INTEGER, DIMENSION(3)                              :: ind, shape_block
      LOGICAL                                            :: pattern_only
      REAL(dp), ALLOCATABLE, TARGET                      :: block(:, :, :)
      REAL(dp), POINTER                                  :: values(:)

      CALL timeset(routineN, handle)
      pattern_only = .FALSE.
      IF (PRESENT(sparsity_only)) pattern_only = sparsity_only
      nblk = SIZE(blk_indices, 1)
      CPASSERT(nblk == SIZE(compressed%nbits))

!$OMP PARALLEL DEFAULT(NONE) SHARED(tensor,blk_indices,compressed,eps,nblk,pattern_only) &
!$OMP PRIVATE(A,b,i,me,nt,ind,shape_block,block,values)
      me = 0
      nt = 1
!$    me = omp_get_thread_num()
!$    nt = omp_get_num_threads()
      ! Disjoint ranges, including the nblk < nt and zero-block cases.
      a = INT(INT(nblk, int_8)*me/nt) + 1
      b = INT(INT(nblk, int_8)*(me + 1)/nt)
      CALL dbt_reserve_blocks(tensor, blk_indices(a:b, :))
!$OMP BARRIER
      IF (.NOT. pattern_only) THEN
!$OMP DO SCHEDULE(DYNAMIC,1)
         DO i = 1, nblk
            ind = blk_indices(i, :)
            CALL dbt_blk_sizes(tensor, ind, shape_block)
            ALLOCATE (block(shape_block(1), shape_block(2), shape_block(3)))
            values(1:SIZE(block)) => block
            CALL unpack_block(compressed%data(compressed%offsets(i):compressed%offsets(i + 1) - 1), &
                              eps, compressed%nbits(i), values)
            CALL dbt_put_block(tensor, ind, shape_block, block)
            NULLIFY (values)
            DEALLOCATE (block)
         END DO
!$OMP END DO
      END IF
!$OMP END PARALLEL
      CALL timestop(handle)
   END SUBROUTINE decompress_ri_tensor

! **************************************************************************************************
!> \brief Signed width from the largest quantized magnitude; zero needs no payload.
!> \param values ...
!> \param eps ...
!> \return ...
! **************************************************************************************************
   FUNCTION block_nbits(values, eps) RESULT(nbits)
      REAL(dp), INTENT(IN)                               :: values(:), eps
      INTEGER                                            :: nbits

      INTEGER(int_8)                                     :: largest_int
      REAL(dp)                                           :: largest, scaled

      largest = MAXVAL(ABS(values))
      nbits = 0
      IF (largest <= eps) RETURN
      scaled = largest*(1.0_dp/eps)
      IF (scaled >= SCALE(1.0_dp, 63)) THEN
         CPABORT("Overflow during RI tensor compression. Increase EPS_FILTER or EPS_STORAGE_SCALING")
      END IF
      largest_int = NINT(scaled, KIND=int_8)
      nbits = INT(BIT_SIZE(largest_int) - LEADZ(largest_int)) + 1
   END FUNCTION block_nbits

! **************************************************************************************************
!> \brief Pack a block using the existing HFX bit-packing kernels.
!> \param values ...
!> \param eps ...
!> \param nbits ...
!> \param packed ...
! **************************************************************************************************
   SUBROUTINE pack_block(values, eps, nbits, packed)
      REAL(dp), INTENT(IN)                               :: values(:), eps
      INTEGER, INTENT(IN)                                :: nbits
      INTEGER(int_8), INTENT(OUT)                        :: packed(:)

      INTEGER                                            :: first, i, n
      INTEGER(int_8)                                     :: offset, shift, work(chunk_size)
      REAL(dp)                                           :: eps_inv

      eps_inv = 1.0_dp/eps
      shift = 0_int_8
      IF (nbits < 64) shift = ISHFT(1_int_8, nbits - 1)
      offset = 1_int_8
      DO first = 1, SIZE(values), chunk_size
         n = MIN(chunk_size, SIZE(values) - first + 1)
         DO i = 1, n
            work(i) = 0_int_8
            IF (ABS(values(first + i - 1)) > eps) work(i) = NINT(values(first + i - 1)*eps_inv, KIND=int_8)
            work(i) = work(i) + shift
         END DO
         IF (nbits == 64) THEN
            packed(offset:offset + n - 1) = work(1:n)
         ELSE
            CALL ints2bits_specific(nbits, n, packed(offset:), work)
         END IF
         ! chunk_size is a multiple of 64, so only the final chunk needs padding.
         offset = offset + (INT(n, int_8)*nbits + 63_int_8)/64_int_8
      END DO
   END SUBROUTINE pack_block

! **************************************************************************************************
!> \brief Unpack a block with fixed-size, thread-private integer scratch.
!> \param packed ...
!> \param eps ...
!> \param nbits ...
!> \param values ...
! **************************************************************************************************
   SUBROUTINE unpack_block(packed, eps, nbits, values)
      INTEGER(int_8), INTENT(IN)                         :: packed(:)
      REAL(dp), INTENT(IN)                               :: eps
      INTEGER, INTENT(IN)                                :: nbits
      REAL(dp), INTENT(OUT)                              :: values(:)

      INTEGER                                            :: first, n
      INTEGER(int_8)                                     :: offset, shift, work(chunk_size)

      IF (nbits == 0) THEN
         values(:) = 0.0_dp
         RETURN
      END IF
      shift = 0_int_8
      IF (nbits < 64) shift = ISHFT(1_int_8, nbits - 1)
      offset = 1_int_8
      DO first = 1, SIZE(values), chunk_size
         n = MIN(chunk_size, SIZE(values) - first + 1)
         IF (nbits == 64) THEN
            work(1:n) = packed(offset:offset + n - 1)
         ELSE
            CALL bits2ints_specific(nbits, n, packed(offset:), work)
         END IF
         values(first:first + n - 1) = REAL(work(1:n) - shift, dp)*eps
         offset = offset + (INT(n, int_8)*nbits + 63_int_8)/64_int_8
      END DO
   END SUBROUTINE unpack_block

END MODULE hfx_ri_compression
