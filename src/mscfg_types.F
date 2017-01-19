!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Types used to generate the molecular SCF guess
!> \par History
!>       10.2014 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
MODULE mscfg_types
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_complete_redistribute, dbcsr_create, dbcsr_distribution_get, &
        dbcsr_distribution_new, dbcsr_distribution_obj, dbcsr_distribution_release, &
        dbcsr_finalize, dbcsr_get_info, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_mp_obj, &
        dbcsr_nfullcols_total, dbcsr_nfullrows_total, dbcsr_release, dbcsr_reserve_block2d, &
        dbcsr_set, dbcsr_type, dbcsr_type_no_symmetry, dbcsr_work_create
   USE kinds,                           ONLY: dp
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mscfg_types'

   INTEGER, PARAMETER, PUBLIC               :: mscfg_max_moset_size = 2

   ! Public types
   PUBLIC :: molecular_scf_guess_env_type

   ! Public subroutines
   PUBLIC :: molecular_scf_guess_env_init, &
             molecular_scf_guess_env_destroy, &
             get_matrix_from_submatrices

   ! Contains data pertaining to molecular_scf_guess calculations
   TYPE molecular_scf_guess_env_type

      ! Useful flags to pass around
      LOGICAL                                           :: is_fast_dirty, &
                                                           is_crystal

      ! Real data
      INTEGER                                           :: nfrags
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE          :: energy_of_frag
      INTEGER, DIMENSION(:), ALLOCATABLE                :: nmosets_of_frag
      TYPE(dbcsr_type), DIMENSION(:, :), ALLOCATABLE  :: mos_of_frag

   END TYPE

CONTAINS

! **************************************************************************************************
!> \brief Allocates data
!> \param env ...
!> \param nfrags   number of entries
!> \par History
!>       2014.10 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE molecular_scf_guess_env_init(env, nfrags)

      TYPE(molecular_scf_guess_env_type)                 :: env
      INTEGER, INTENT(IN)                                :: nfrags

      CHARACTER(len=*), PARAMETER :: routineN = 'molecular_scf_guess_env_init', &
         routineP = moduleN//':'//routineN

! check if the number of fragments is already set
!IF (env%nfrags.ne.0) THEN
!   ! do not allow re-initialization
!   ! to prevent recursive calls
!   CPPostcondition(.FALSE.,cp_failure_level,routineP,failure)
!ENDIF

      env%nfrags = nfrags
      IF (nfrags .GT. 0) THEN
         ALLOCATE (env%energy_of_frag(nfrags))
         ALLOCATE (env%nmosets_of_frag(nfrags))
         ALLOCATE (env%mos_of_frag(nfrags, mscfg_max_moset_size))
      ENDIF

   END SUBROUTINE molecular_scf_guess_env_init

! **************************************************************************************************
!> \brief Destroyes both data and environment
!> \param env ...
!> \par History
!>       2014.10 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE molecular_scf_guess_env_destroy(env)

      TYPE(molecular_scf_guess_env_type)                 :: env

      CHARACTER(len=*), PARAMETER :: routineN = 'molecular_scf_guess_env_destroy', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ifrag, jfrag

      IF (ALLOCATED(env%mos_of_frag)) THEN
         DO ifrag = 1, SIZE(env%mos_of_frag, 1)
            DO jfrag = 1, env%nmosets_of_frag(ifrag)
               CALL dbcsr_release(env%mos_of_frag(ifrag, jfrag))
            ENDDO
         ENDDO
         DEALLOCATE (env%mos_of_frag)
      ENDIF
      IF (ALLOCATED(env%energy_of_frag)) DEALLOCATE (env%energy_of_frag)
      IF (ALLOCATED(env%nmosets_of_frag)) DEALLOCATE (env%nmosets_of_frag)

      env%nfrags = 0

   END SUBROUTINE molecular_scf_guess_env_destroy

! **************************************************************************************************
!> \brief Creates a distributed matrix from MOs on fragments
!> \param mscfg_env   env containing MOs of fragments
!> \param matrix_out   all existing blocks will be deleted!
!> \param iset   which set of MOs in mscfg_env has to be converted (e.g. spin)
!> \par History
!>       10.2014 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE get_matrix_from_submatrices(mscfg_env, matrix_out, iset)

      TYPE(molecular_scf_guess_env_type), INTENT(IN)     :: mscfg_env
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_out
      INTEGER, INTENT(IN)                                :: iset

      CHARACTER(len=*), PARAMETER :: routineN = 'get_matrix_from_submatrices', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ifrag
      INTEGER, DIMENSION(2)                              :: matrix_size, offset, submatrix_size
      TYPE(dbcsr_type)                                   :: matrix_temp

      CALL timeset(routineN, handle)

      CPASSERT(iset .LE. mscfg_max_moset_size)

      CALL dbcsr_create(matrix_temp, &
                        template=matrix_out, &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_set(matrix_out, 0.0_dp)

      matrix_size(1) = dbcsr_nfullrows_total(matrix_out)
      matrix_size(2) = dbcsr_nfullcols_total(matrix_out)

      ! assume that the initial offset is zero
      offset(1) = 0
      offset(2) = 0

      DO ifrag = 1, mscfg_env%nfrags

         CPASSERT(iset .LE. mscfg_env%nmosets_of_frag(ifrag))

         submatrix_size(1) = dbcsr_nfullrows_total(mscfg_env%mos_of_frag(ifrag, iset))
         submatrix_size(2) = dbcsr_nfullcols_total(mscfg_env%mos_of_frag(ifrag, iset))

         CALL copy_submatrix_into_matrix(mscfg_env%mos_of_frag(ifrag, iset), &
                                         matrix_temp, offset, submatrix_size, matrix_size)

         CALL dbcsr_add(matrix_out, matrix_temp, 1.0_dp, 1.0_dp)

         offset(1) = offset(1)+submatrix_size(1)
         offset(2) = offset(2)+submatrix_size(2)

      ENDDO

      ! Check that the accumulated size of submatrices
      ! is exactly the same as the size of the big matrix
      ! This is to prevent unexpected conversion errors
      ! If however such conversion is intended - remove these safeguards
      CPASSERT(offset(1) .EQ. matrix_size(1))
      CPASSERT(offset(2) .EQ. matrix_size(2))

      CALL dbcsr_release(matrix_temp)

      CALL timestop(handle)

   END SUBROUTINE get_matrix_from_submatrices

! **************************************************************************************************
!> \brief Copies a distributed dbcsr submatrix into a distributed dbcsr matrix
!> \param submatrix_in ...
!> \param matrix_out   all existing blocks will be deleted!
!> \param offset ...
!> \param submatrix_size ...
!> \param matrix_size ...
!> \par History
!>       10.2014 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE copy_submatrix_into_matrix(submatrix_in, matrix_out, &
                                         offset, submatrix_size, matrix_size)

      TYPE(dbcsr_type), INTENT(IN)                       :: submatrix_in
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix_out
      INTEGER, DIMENSION(2), INTENT(IN)                  :: offset, submatrix_size, matrix_size

      CHARACTER(len=*), PARAMETER :: routineN = 'copy_submatrix_into_matrix', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: add_blocks_after, dimen, iblock_col, &
                                                            iblock_row, iblock_size, nblocks, &
                                                            nblocks_new, start_index, trailing_size
      INTEGER, DIMENSION(2)                              :: add_blocks_before
      INTEGER, DIMENSION(:), POINTER :: blk_distr, blk_sizes, block_sizes_new, cluster_distr, &
         cluster_distr_new, col_cluster_new, col_distr_new, col_sizes_new, distr_new_array, &
         row_cluster_new, row_distr_new, row_sizes_new
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: data_p, p_new_block
      TYPE(dbcsr_distribution_obj)                       :: dist_new, dist_qs
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_mp_obj)                                 :: mp_env
      TYPE(dbcsr_type)                                   :: matrix_new

! obtain distribution of the submatrix

      CALL dbcsr_get_info(submatrix_in, distribution=dist_qs)

      DO dimen = 1, 2 ! 1 - row, 2 - column dimension

         add_blocks_before(dimen) = 0
         add_blocks_after = 0
         start_index = 1
         trailing_size = matrix_size(dimen)-offset(dimen)-submatrix_size(dimen)
         IF (offset(dimen) .GT. 0) THEN
            add_blocks_before(dimen) = add_blocks_before(dimen)+1
            start_index = 2
         ENDIF
         IF (trailing_size .GT. 0) THEN
            add_blocks_after = add_blocks_after+1
         ENDIF

         NULLIFY (cluster_distr_new)

         IF (dimen == 1) THEN !rows
            CALL dbcsr_distribution_get(dist_qs, row_dist=blk_distr, row_cluster_dist=cluster_distr)
            CALL dbcsr_get_info(submatrix_in, row_blk_size=blk_sizes)
         ELSE !columns
            CALL dbcsr_distribution_get(dist_qs, col_dist=blk_distr, col_cluster_dist=cluster_distr)
            CALL dbcsr_get_info(submatrix_in, col_blk_size=blk_sizes)
         ENDIF
         nblocks = SIZE(blk_sizes) ! number of blocks in the small matrix

         nblocks_new = nblocks+add_blocks_before(dimen)+add_blocks_after
         ALLOCATE (block_sizes_new(nblocks_new))
         ALLOCATE (distr_new_array(nblocks_new))
         IF (ASSOCIATED(cluster_distr)) THEN
            ALLOCATE (cluster_distr_new(nblocks_new))
         END IF
         IF (add_blocks_before(dimen) .GT. 0) THEN
            block_sizes_new(1) = offset(dimen)
            distr_new_array(1) = 0
            IF (ASSOCIATED(cluster_distr)) THEN
               cluster_distr_new(1) = 0
            END IF
         ENDIF
         block_sizes_new(start_index:nblocks+start_index-1) = blk_sizes(1:nblocks)
         distr_new_array(start_index:nblocks+start_index-1) = blk_distr(1:nblocks)
         IF (ASSOCIATED(cluster_distr)) THEN
            cluster_distr_new(start_index:nblocks+start_index-1) = cluster_distr(1:nblocks)
         END IF
         IF (add_blocks_after .GT. 0) THEN
            block_sizes_new(nblocks_new) = trailing_size
            distr_new_array(nblocks_new) = 0
            IF (ASSOCIATED(cluster_distr)) THEN
               cluster_distr_new(nblocks_new) = 0
            END IF
         ENDIF

         ! create final arrays
         IF (dimen == 1) THEN !rows
            row_sizes_new => block_sizes_new
            row_distr_new => distr_new_array
            row_cluster_new => cluster_distr_new
         ELSE !columns
            col_sizes_new => block_sizes_new
            col_distr_new => distr_new_array
            col_cluster_new => cluster_distr_new
         ENDIF
      ENDDO ! both rows and columns are done

      ! Create the distribution
      CALL dbcsr_distribution_get(dist_qs, mp_env=mp_env)
      CALL dbcsr_distribution_new(dist_new, &
                                  mp_env, row_distr_new, col_distr_new, &
                                  row_cluster_new, col_cluster_new, &
                                  reuse_arrays=.TRUE.)

      ! Create big the matrix
      CALL dbcsr_create(matrix_new, "BIG_AND_FAKE", &
                        dist_new, dbcsr_type_no_symmetry, &
                        row_sizes_new, col_sizes_new, &
                        reuse_arrays=.TRUE.)
      CALL dbcsr_distribution_release(dist_new)

      !CALL dbcsr_finalize(matrix_new)

      ! copy blocks of the small matrix to the big matrix
      !mynode = dbcsr_mp_mynode(dbcsr_distribution_mp(dbcsr_distribution(matrix_new)))
      CALL dbcsr_work_create(matrix_new, work_mutable=.TRUE.)

      ! iterate over local blocks of the small matrix
      CALL dbcsr_iterator_start(iter, submatrix_in)

      DO WHILE (dbcsr_iterator_blocks_left(iter))

         CALL dbcsr_iterator_next_block(iter, iblock_row, iblock_col, data_p, row_size=iblock_size)

         ! it is important that distribution of the big matrix is the same as
         ! that of the small matrix but has the same number of columns and rows
         ! as the super-system matrix. this is necessary for complete redistribute
         ! to work
         NULLIFY (p_new_block)
         CALL dbcsr_reserve_block2d(matrix_new, &
                                    iblock_row+add_blocks_before(1), &
                                    iblock_col+add_blocks_before(2), &
                                    p_new_block)

         CPASSERT(ASSOCIATED(p_new_block))
         CPASSERT(SIZE(p_new_block, 1) .EQ. SIZE(data_p, 1))
         CPASSERT(SIZE(p_new_block, 2) .EQ. SIZE(data_p, 2))

         p_new_block(:, :) = data_p(:, :)

      ENDDO
      CALL dbcsr_iterator_stop(iter)

      CALL dbcsr_finalize(matrix_new)

      ! finally call complete redistribute to get the matrix of the entire system
      CALL dbcsr_set(matrix_out, 0.0_dp)
      CALL dbcsr_complete_redistribute(matrix_new, matrix_out)
      CALL dbcsr_release(matrix_new)

   END SUBROUTINE copy_submatrix_into_matrix

END MODULE mscfg_types

