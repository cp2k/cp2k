!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief General methods for testing DBCSR tensors.
!> \author Patrick Seewald
! **************************************************************************************************
MODULE dbcsr_tensor_test

#:include "dbcsr_tensor.fypp"
#:setvar maxdim maxrank
#:setvar ndims range(2,maxdim+1)

   USE dbcsr_api,                       ONLY: ${uselist(dtype_float_param)}$
   USE dbcsr_tensor,                    ONLY: &
        dbcsr_t_copy, dbcsr_t_get_block, dbcsr_t_iterator_type, dbcsr_t_iterator_blocks_left, &
        dbcsr_t_iterator_next_block, dbcsr_t_iterator_start, dbcsr_t_iterator_stop, dbcsr_t_ndims, &
        dbcsr_t_reserve_blocks, dbcsr_t_get_stored_coordinates, dbcsr_t_put_block
   USE dbcsr_tensor_block,              ONLY: block_nd
   USE dbcsr_tensor_types,              ONLY: dbcsr_t_create,&
                                              dbcsr_t_destroy,&
                                              dbcsr_t_type,&
                                              dbcsr_t_distribution_type,&
                                              dbcsr_t_distribution_destroy,&
                                              dims_tensor,&
                                              dbcsr_t_distribution_new,&
                                              dbcsr_t_nd_mp_comm
   USE kinds,                           ONLY: ${uselist(dtype_float_prec)}$,&
                                              default_string_length
   USE message_passing,                 ONLY: mp_environ, &
                                              mp_bcast, &
                                              mp_comm_free
   USE allocate_wrap,                   ONLY: allocate_any
   USE dbcsr_tensor_index,              ONLY: combine_index
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_tensor_test'

   PUBLIC :: &
      test_tensor_formats

CONTAINS

! **************************************************************************************************
!> \brief check if two (arbitrarily mapped and distributed) tensors are equal.
!> \param tensor1 ...
!> \param tensor2 ...
!> \retval dbcsr_t_equal ...
! **************************************************************************************************
   FUNCTION dbcsr_t_equal(tensor1, tensor2)
      TYPE(dbcsr_t_type), INTENT(INOUT)          :: tensor1, tensor2
      LOGICAL                                    :: dbcsr_t_equal

      INTEGER                                    :: blk
      TYPE(dbcsr_t_type)                         :: tensor2_tmp
      TYPE(dbcsr_t_iterator_type)                :: iter
      TYPE(block_nd)                             :: blk_data1, blk_data2
      INTEGER, DIMENSION(dbcsr_t_ndims(tensor1)) :: blk_size, ind_nd

      ! create a copy of tensor2 that has exact same data format as tensor1
      CALL dbcsr_t_create(tensor1, tensor2_tmp)
      CALL dbcsr_t_reserve_blocks(tensor1, tensor2_tmp)
      CALL dbcsr_t_copy(tensor2, tensor2_tmp)

      dbcsr_t_equal = .TRUE.

      CALL dbcsr_t_iterator_start(iter, tensor1)

      DO WHILE (dbcsr_t_iterator_blocks_left(iter))
         CALL dbcsr_t_iterator_next_block(iter, ind_nd, blk, blk_size=blk_size)
         CALL dbcsr_t_get_block(tensor1, ind_nd, blk_data1)
         CALL dbcsr_t_get_block(tensor2_tmp, ind_nd, blk_data2)

         IF (.NOT. blocks_equal(blk_data1, blk_data2)) THEN
            dbcsr_t_equal = .FALSE.
         ENDIF
      ENDDO

      CALL dbcsr_t_iterator_stop(iter)
      CALL dbcsr_t_destroy(tensor2_tmp)
   END FUNCTION

! **************************************************************************************************
!> \brief check if two blocks are equal
!> \param block1 ...
!> \param block2 ...
!> \retval blocks_equal ...
! **************************************************************************************************
   PURE FUNCTION blocks_equal(block1, block2)
      TYPE(block_nd), INTENT(IN) :: block1, block2
      LOGICAL                    :: blocks_equal

      SELECT CASE (block1%data_type)
#:for dprec, dparam, dsuffix in list(zip(dtype_float_prec, dtype_float_param, dtype_float_suffix))
      CASE (${dparam}$)
         blocks_equal = MAXVAL(ABS(block1%${dsuffix}$%blk-block2%${dsuffix}$%blk)) .LT. 1.0E-12_${dprec}$
#:endfor
      END SELECT

   END FUNCTION

! **************************************************************************************************
!> \brief Compute factorial
!> \param n ...
!> \retval factorial ...
! **************************************************************************************************
   PURE FUNCTION factorial(n)
      INTEGER, INTENT(IN) :: n
      INTEGER             :: k
      INTEGER             :: factorial
      factorial = PRODUCT((/(k, k=1, n)/))
   END FUNCTION

! **************************************************************************************************
!> \brief Compute all permutations p of (1, 2, ..., n)
!> \param n ...
!> \param p ...
! **************************************************************************************************
   SUBROUTINE permute(n, p)
      INTEGER, INTENT(IN)                              :: n
      INTEGER                                          :: i, c
      INTEGER, DIMENSION(n)                            :: pp
      INTEGER, DIMENSION(n, factorial(n)), INTENT(OUT) :: p

      pp = [(i, i=1, n)]
      c = 1
      CALL perm(1)
   CONTAINS
      RECURSIVE SUBROUTINE perm(i)
         INTEGER, INTENT(IN) :: i
         INTEGER :: j, t
         IF (i == n) THEN
            p(:, c) = pp(:)
            c = c+1
         ELSE
            DO j = i, n
               t = pp(i)
               pp(i) = pp(j)
               pp(j) = t
               call perm(i+1)
               t = pp(i)
               pp(i) = pp(j)
               pp(j) = t
            END DO
         END IF
      END SUBROUTINE
   END SUBROUTINE

! **************************************************************************************************
!> \brief Test equivalence of all tensor formats, using a random distribution.
!> \param ndims tensor rank
!> \param mp_comm_2d communicator
!> \param unit_nr output unit, needs to be a valid unit number on all mpi ranks
!> \param verbose if .TRUE., print all tensor blocks
!> \param blk_size_1 block sizes along dimension 1
!> \param blk_ind_1 index along dimension 1 of non-zero blocks
! **************************************************************************************************
   SUBROUTINE test_tensor_formats(ndims, mp_comm_2d, unit_nr, verbose, &
                                  ${varlist("blk_size")}$, &
                                  ${varlist("blk_ind")}$)
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: ${varlist("blk_size")}$
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: ${varlist("blk_ind")}$
      INTEGER, INTENT(IN)                         :: ndims
      INTEGER, INTENT(IN)                         :: unit_nr
      LOGICAL, INTENT(IN)                         :: verbose
      INTEGER, INTENT(IN)                         :: mp_comm_2d
      INTEGER                                     :: mp_comm_ref, mp_comm_test
      TYPE(dbcsr_t_distribution_type)             :: dist1, dist2
      TYPE(dbcsr_t_type)                          :: tensor1, tensor2
      INTEGER                                     :: isep, iblk
      INTEGER, DIMENSION(:), ALLOCATABLE          :: ${varlist("dist1")}$, &
                                                     ${varlist("dist2")}$
      INTEGER                                     :: nblks, imap
      INTEGER, DIMENSION(ndims)                   :: pdims_ref, pdims_test, myploc
      LOGICAL                                     :: eql
      INTEGER                                     :: iperm, idist, icount
      INTEGER, DIMENSION(:), ALLOCATABLE          :: map1, map2, map1_ref, map2_ref
      INTEGER, DIMENSION(ndims, factorial(ndims)) :: perm
      LOGICAL, DIMENSION(ndims)                   :: periods
      INTEGER                                     :: io_unit
      INTEGER                                     :: mynode, numnodes
      CHARACTER(LEN=default_string_length)        :: tensor_name

      CALL mp_environ(numnodes, mynode, mp_comm_2d)

      io_unit = 0
      IF (mynode .EQ. 0) io_unit = unit_nr

      CALL permute(ndims, perm)
      ALLOCATE(map1_ref(ndims/2), map2_ref(ndims-ndims/2))
      map1_ref(:) = perm(1:ndims/2, 1)
      map2_ref(:) = perm(ndims/2+1:ndims, 1)

      mp_comm_ref = dbcsr_t_nd_mp_comm(mp_comm_2d, map1_ref, map2_ref)
      CALL mp_environ(numnodes, mynode, mp_comm_ref)
      CALL mp_environ(mp_comm_ref, ndims, pdims_ref, myploc, periods)

      IF(io_unit > 0) THEN
         WRITE(io_unit, *)
         WRITE(io_unit, '(A)') repeat("-",80)
         WRITE(io_unit, '(A,1X,I1)') "Testing matrix representations of tensor rank", ndims
         WRITE(io_unit, '(A)') repeat("-",80)
         WRITE(io_unit,'(A)') "Block sizes:"

#:for dim in range(1, maxdim+1)
         IF (ndims>=${dim}$) THEN
            WRITE(io_unit,'(T4,A,1X,I1,A,1X)', advance='no') 'Dim', ${dim}$,':'
            DO iblk=1, SIZE(blk_size_${dim}$)
               WRITE(io_unit,'(I2,1X)', advance='no') blk_size_${dim}$(iblk)
            ENDDO
            WRITE(io_unit,*)
         ENDIF
#:endfor

         WRITE(io_unit,'(A)') "Non-zero blocks:"
         DO iblk=1,SIZE(blk_ind_1)
#:for ndim in ndims
            IF (ndims==${ndim}$) THEN
               WRITE(io_unit,'(T4,A, I3, A, ${ndim}$I3, 1X, A)') &
                  'Block',iblk, ': (', ${varlist("blk_ind", nmax=ndim, suffix='(iblk)')}$, ')'
            ENDIF
#:endfor
         ENDDO

         WRITE(io_unit, *)
         WRITE(io_unit,'(A,1X)', advance='no') "Reference map:"
         WRITE (io_unit, '(A1,1X)', advance='no') "("
         DO imap = 1, SIZE(map1_ref)
            WRITE (io_unit, '(I1,1X)', advance='no') map1_ref(imap)
         ENDDO
         WRITE (io_unit, '(A1,1X)', advance='no') "|"
         DO imap = 1, SIZE(map2_ref)
            WRITE (io_unit, '(I1,1X)', advance='no') map2_ref(imap)
         ENDDO
         WRITE (io_unit, '(A1)') ")"

      ENDIF

      icount = 0
      DO iperm = 1, factorial(ndims)
         DO isep = 1, ndims-1
            icount = icount+1

            ALLOCATE(map1(isep), map2(ndims-isep))
            map1(:) = perm(1:isep,iperm)
            map2(:) = perm(isep+1:ndims,iperm)

            mp_comm_test = dbcsr_t_nd_mp_comm(mp_comm_2d, map1, map2)
            CALL mp_environ(numnodes, mynode, mp_comm_test)
            CALL mp_environ(mp_comm_test, ndims, pdims_test, myploc, periods)

#:for dim in range(1, maxdim+1)
            IF (${dim}$ <= ndims) THEN
               nblks = SIZE(blk_size_${dim}$)
               CALL random_dist(dist1_${dim}$, nblks, pdims_test(${dim}$), mp_comm_test)
               CALL random_dist(dist2_${dim}$, nblks, pdims_ref(${dim}$), mp_comm_ref)
            ENDIF
#:endfor

            IF(io_unit > 0) THEN
               WRITE(io_unit, *)
               WRITE(tensor_name, '(A,1X,I3,1X)') "Test", icount
               WRITE(io_unit, '(A,A,1X)', advance='no') TRIM(tensor_name),':'
               WRITE (io_unit, '(A1,1X)', advance='no') "("
               DO imap = 1, SIZE(map1)
                  WRITE (io_unit, '(I1,1X)', advance='no') map1(imap)
               ENDDO
               WRITE (io_unit, '(A1,1X)', advance='no') "|"
               DO imap = 1, SIZE(map2)
                  WRITE (io_unit, '(I1,1X)', advance='no') map2(imap)
               ENDDO
               WRITE (io_unit, '(A1)') ")"

               WRITE(io_unit,'(T4,A)') "Reference distribution:"
#:for dim in range(1, maxdim+1)
               IF (${dim}$ <= ndims) THEN
                  WRITE(io_unit, '(T7,A,1X)', advance='no') "Dist vec ${dim}$:"
                  DO idist = 1, SIZE(dist2_${dim}$)
                     WRITE(io_unit, '(I2,1X)', advance='no') dist2_${dim}$(idist)
                  ENDDO
                  WRITE(io_unit, *)
               ENDIF
#:endfor

               WRITE(io_unit,'(T4,A)') "Test distribution:"
#:for dim in range(1, maxdim+1)
               IF (${dim}$ <= ndims) THEN
                  WRITE(io_unit, '(T7,A,1X)', advance='no') "Dist vec ${dim}$:"
                  DO idist = 1, SIZE(dist2_${dim}$)
                     WRITE(io_unit, '(I2,1X)', advance='no') dist1_${dim}$(idist)
                  ENDDO
                  WRITE(io_unit, *)
               ENDIF
#:endfor
            ENDIF

#:for ndim in ndims
            IF (ndims==${ndim}$) THEN
               CALL dbcsr_t_distribution_new(dist2, mp_comm_ref, mp_comm_2d, map1_ref, map2_ref, ${varlist("dist2", nmax=ndim)}$)
               CALL dbcsr_t_create(tensor2, "Ref", dist2, map1_ref, map2_ref, &
                                   dbcsr_type_real_8, ${varlist("blk_size", nmax=ndim)}$)
               CALL allocate_and_fill_test_tensor(tensor2, mp_comm_ref, ${varlist("blk_ind", nmax=ndim)}$)
            ENDIF
#:endfor

            IF(verbose) CALL write_tensor_blocks(tensor2, io_unit, unit_nr)

#:for ndim in ndims
            IF (ndims==${ndim}$) THEN
               CALL dbcsr_t_distribution_new(dist1, mp_comm_test, mp_comm_2d, map1, map2, ${varlist("dist1", nmax=ndim)}$)
               CALL dbcsr_t_create(tensor1, tensor_name, dist1, map1, map2, &
                                   dbcsr_type_real_8, ${varlist("blk_size", nmax=ndim)}$)
               CALL allocate_and_fill_test_tensor(tensor1, mp_comm_test, ${varlist("blk_ind", nmax=ndim)}$)
            ENDIF
#:endfor

            IF(verbose) CALL write_tensor_blocks(tensor1, io_unit, unit_nr)

            eql = dbcsr_t_equal(tensor1, tensor2)


            IF (.NOT. eql) THEN
               IF (io_unit > 0) WRITE(io_unit,'(A,1X,A)') TRIM(tensor_name), 'failed!'
               CPABORT('')
            ELSE
               IF (io_unit > 0) WRITE(io_unit,'(A,1X,A)') TRIM(tensor_name), 'passed!'
            ENDIF
            DEALLOCATE(map1, map2)

            CALL dbcsr_t_destroy(tensor1)
            CALL dbcsr_t_distribution_destroy(dist1)

            CALL mp_comm_free(mp_comm_test)
            CALL dbcsr_t_destroy(tensor2)
            CALL dbcsr_t_distribution_destroy(dist2)

         ENDDO
      ENDDO
      CALL mp_comm_free(mp_comm_ref)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Create test distribution
!> \param dist_array ...
!> \param dist_size ...
!> \param seed ...
!> \param nbins ...
! **************************************************************************************************
   SUBROUTINE random_dist(dist_array, dist_size, nbins, mp_comm)
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(out) :: dist_array
      INTEGER, INTENT(IN)                             :: dist_size, nbins, mp_comm
      REAL, DIMENSION(dist_size)                      :: rn
      INTEGER                                         :: numnodes, mynode

      CALL mp_environ(numnodes, mynode, mp_comm)

      IF (mynode .EQ. 0) THEN
         CALL RANDOM_NUMBER(rn)
      ENDIF
      CALL mp_bcast(rn, 0, mp_comm)

      ALLOCATE (dist_array(dist_size))

      dist_array(:) = FLOOR(rn*nbins)

   END SUBROUTINE random_dist

! **************************************************************************************************
!> \brief Allocate and fill test tensor - entries are enumerated by their index s.t. they only depend
!>        on global properties of the tensor but not on distribution, matrix representation, etc.
!> \param tensor ...
!> \param mp_comm communicator
!> \param blk_ind_1 index along dimension 1 of non-zero blocks
! **************************************************************************************************
   SUBROUTINE allocate_and_fill_test_tensor(tensor, mp_comm, ${varlist("blk_ind")}$)
      TYPE(dbcsr_t_type), INTENT(INOUT)                  :: tensor
      INTEGER, INTENT(IN)                                :: mp_comm
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL        :: ${varlist("blk_ind")}$
      INTEGER                                            :: blk, numnodes, mynode

      INTEGER                                            :: i, ib, my_nblks_alloc, nblks_alloc, proc
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ${varlist("my_blk_ind")}$
      INTEGER, DIMENSION(dbcsr_t_ndims(tensor))          :: blk_index, blk_offset, blk_size, &
                                                            tensor_dims
      INTEGER, DIMENSION(:,:), ALLOCATABLE               :: ind_nd
#:for ndim in ndims
      REAL(KIND=real_8), ALLOCATABLE, &
         DIMENSION(${shape_colon(ndim)}$)                :: blk_values_${ndim}$
#:endfor
      TYPE(dbcsr_t_iterator_type)                        :: iterator

      nblks_alloc = SIZE(blk_ind_1)
      CALL mp_environ(numnodes, mynode, mp_comm)

      ALLOCATE(ind_nd(nblks_alloc, dbcsr_t_ndims(tensor)))
      my_nblks_alloc = 0
      DO ib = 1, nblks_alloc
#:for ndim in ndims
         IF(dbcsr_t_ndims(tensor) == ${ndim}$) THEN
            ind_nd(ib,:) = [${varlist("blk_ind", nmax=ndim, suffix="(ib)")}$]
         ENDIF
#:endfor
         CALL dbcsr_t_get_stored_coordinates(tensor, ind_nd(ib,:), proc)
         IF (proc == mynode) THEN
            my_nblks_alloc = my_nblks_alloc+1
         ENDIF
      ENDDO

#:for dim in range(1, maxdim+1)
      IF(dbcsr_t_ndims(tensor) >= ${dim}$) THEN
         ALLOCATE (my_blk_ind_${dim}$(my_nblks_alloc))
      ENDIF
#:endfor

      i = 0
      DO ib = 1, nblks_alloc
         CALL dbcsr_t_get_stored_coordinates(tensor, ind_nd(ib,:), proc)
         IF (proc == mynode) THEN
            i = i+1
#:for dim in range(1, maxdim+1)
            IF(dbcsr_t_ndims(tensor) >= ${dim}$) THEN
               my_blk_ind_${dim}$(i) = blk_ind_${dim}$(ib)
            ENDIF
#:endfor
         ENDIF
      ENDDO

#:for ndim in ndims
      IF(dbcsr_t_ndims(tensor) == ${ndim}$) THEN
         CALL dbcsr_t_reserve_blocks(tensor, ${varlist("my_blk_ind", nmax=ndim)}$)
      ENDIF
#:endfor

      CALL dbcsr_t_iterator_start(iterator, tensor)
      DO WHILE (dbcsr_t_iterator_blocks_left(iterator))
         CALL dbcsr_t_iterator_next_block(iterator, blk_index, blk, blk_size=blk_size, blk_offset=blk_offset)

#:for ndim in ndims
         IF(dbcsr_t_ndims(tensor) == ${ndim}$) THEN
            CALL allocate_any(blk_size, blk_values_${ndim}$)
            CALL dims_tensor(tensor, tensor_dims)
            CALL enumerate_block_elements(blk_size, blk_offset, tensor_dims, blk_${ndim}$=blk_values_${ndim}$)
            CALL dbcsr_t_put_block(tensor, blk_index, blk_size, blk_values_${ndim}$)
            DEALLOCATE (blk_values_${ndim}$)
         ENDIF
#:endfor
      ENDDO
      CALL dbcsr_t_iterator_stop(iterator)

   END SUBROUTINE

! **************************************************************************************************
!> \brief Enumerate tensor entries in block
!> \param blk_size size of block
!> \param blk_offset block offset (indices of first element)
!> \param tensor_size global tensor sizes
! **************************************************************************************************
   SUBROUTINE enumerate_block_elements(blk_size, blk_offset, tensor_size, ${varlist("blk", nmin=2)}$)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: blk_size, blk_offset, tensor_size
#:for ndim in ndims
      REAL(KIND=real_8), DIMENSION(${shape_colon(ndim)}$), &
         OPTIONAL, INTENT(OUT)                           :: blk_${ndim}$
#:endfor
      INTEGER                                            :: ndim
      INTEGER, DIMENSION(SIZE(blk_size))                 :: arr_ind, tens_ind
      INTEGER                                            :: ${varlist("i")}$

      ndim = SIZE(tensor_size)

#:for ndim in ndims
      IF (ndim == ${ndim}$) THEN
#:for idim in range(ndim,0,-1)
            DO i_${idim}$ = 1, blk_size(${idim}$)
#:endfor
               arr_ind(:) = [${varlist("i", nmax=ndim)}$]
               tens_ind(:) = arr_ind(:) + blk_offset(:) - 1
               blk_${ndim}$(${arrlist("arr_ind", nmax=ndim)}$) = combine_index(tens_ind, tensor_size)
#:for idim in range(ndim,0,-1)
            ENDDO
#:endfor
      ENDIF
#:endfor

   END SUBROUTINE

! **************************************************************************************************
!> \brief Write all tensor blocks
!> \param tensor ...
!> \param io_unit_master for global output
!> \param io_unit_all for local output
! **************************************************************************************************
   SUBROUTINE write_tensor_blocks(tensor, io_unit_master, io_unit_all)
      TYPE(dbcsr_t_type), INTENT(INOUT)                  :: tensor
      INTEGER, INTENT(IN)                                :: io_unit_master, io_unit_all

      INTEGER                                            :: blk
      INTEGER, DIMENSION(dbcsr_t_ndims(tensor))          :: blk_index, blk_size
#:for ndim in ndims
      REAL(KIND=real_8), ALLOCATABLE, &
         DIMENSION(${shape_colon(ndim)}$)                :: blk_values_${ndim}$
#:endfor
      TYPE(dbcsr_t_iterator_type)                        :: iterator
      INTEGER                                            :: proc

      IF(io_unit_master > 0) THEN
         WRITE (io_unit_master, '(T7,A)') "(block index) @ process: (array index) value"
      ENDIF
      CALL dbcsr_t_iterator_start(iterator, tensor)
      DO WHILE (dbcsr_t_iterator_blocks_left(iterator))
         CALL dbcsr_t_iterator_next_block(iterator, blk_index, blk, blk_size=blk_size)
         CALL dbcsr_t_get_stored_coordinates(tensor, blk_index, proc)
#:for ndim in ndims
         IF(dbcsr_t_ndims(tensor) == ${ndim}$) THEN
            CALL allocate_any(blk_size, blk_values_${ndim}$)
            CALL dbcsr_t_get_block(tensor, blk_index, blk_size, blk_values_${ndim}$)
            CALL write_block(tensor%name, blk_size, blk_index, proc, io_unit_all, blk_values_${ndim}$=blk_values_${ndim}$)
            DEALLOCATE (blk_values_${ndim}$)
         ENDIF
#:endfor
      ENDDO
      CALL dbcsr_t_iterator_stop(iterator)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Write a tensor block
!> \param name tensor name
!> \param blk_size block size
!> \param blk_index block index
!> \param proc which process am I
!> \param blk_values_2 block values
!> \param unit_nr unit number
! **************************************************************************************************
   SUBROUTINE write_block(name, blk_size, blk_index, proc, unit_nr, &
                          ${varlist("blk_values",nmin=2)}$)

      CHARACTER(LEN=*), INTENT(IN)                       :: name
      INTEGER, DIMENSION(:), INTENT(IN)                  :: blk_size, blk_index
#:for ndim in ndims
      REAL(KIND=real_8), &
         DIMENSION(${arrlist("blk_size", nmax=ndim)}$), &
         INTENT(IN), OPTIONAL                            :: blk_values_${ndim}$
#:endfor
      INTEGER, INTENT(IN)                                :: unit_nr, proc
      INTEGER                                            :: ${varlist("i")}$
      INTEGER                                            :: ndim

      ndim = SIZE(blk_size)

      IF (unit_nr > 0) THEN
#:for ndim in ndims
         IF (ndim == ${ndim}$) THEN
#:for idim in range(ndim,0,-1)
            DO i_${idim}$ = 1, blk_size(${idim}$)
#:endfor
               WRITE (unit_nr, '(T7,A,T16,A,${ndim}$I3,1X,A,1X,I3,A,1X,A,${ndim}$I3,1X,A,1X,I7)') &
                  TRIM(name),"(", blk_index,") @", proc, ':', &
                  "(",${varlist("i", nmax=ndim)}$,")", &
                  INT(blk_values_${ndim}$(${varlist("i", nmax=ndim)}$))
#:for idim in range(ndim,0,-1)
            ENDDO
#:endfor
         ENDIF
#:endfor
      ENDIF
   END SUBROUTINE

END MODULE
