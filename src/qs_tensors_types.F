!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Utility methods to build 3-center integral tensors of various types.
! **************************************************************************************************

MODULE qs_tensors_types
   USE cp_array_utils,                  ONLY: cp_1d_i_p_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_create,&
                                              cp_blacs_env_release,&
                                              cp_blacs_env_type
   USE cp_para_env,                     ONLY: cp_para_env_create,&
                                              cp_para_env_release
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_type_real_8
   USE dbcsr_tensor_api,                ONLY: dbcsr_t_create,&
                                              dbcsr_t_distribution_destroy,&
                                              dbcsr_t_distribution_new,&
                                              dbcsr_t_distribution_type,&
                                              dbcsr_t_mp_environ_pgrid,&
                                              dbcsr_t_pgrid_type,&
                                              dbcsr_t_type
   USE distribution_2d_types,           ONLY: distribution_2d_create_prv => distribution_2d_create,&
                                              distribution_2d_release,&
                                              distribution_2d_type
   USE message_passing,                 ONLY: mp_cart_sub,&
                                              mp_comm_free,&
                                              mp_environ
   USE particle_types,                  ONLY: particle_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_iterator_p_type,&
                                              neighbor_list_set_p_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tensors_types'

   PUBLIC :: distribution_3d_type, neighbor_list_3c_type, neighbor_list_3c_iterator_type, &
             distribution_2d_create, distribution_3d_create, distribution_3d_destroy, cyclic_tensor_dist, &
             split_block_sizes, create_3c_tensor, create_2c_tensor, contiguous_tensor_dist

   INTEGER, PARAMETER, PUBLIC :: symmetric_none = 0, symmetric_ij = 1, symmetric_jk = 2, symmetrik_ik = 3, symmetric_ijk = 4

   INTEGER, PARAMETER, PUBLIC :: libcusmm_maxblocksize = 45
   !! upper limit for block sizes still covered by lib cusmm

   TYPE distribution_3d_type
      TYPE(distribution_2d_type), POINTER :: dist_2d_1 => NULL(), dist_2d_2 => NULL()
      INTEGER :: comm_3d, comm_2d_1, comm_2d_2
      LOGICAL :: owns_comm
   END TYPE distribution_3d_type

   TYPE neighbor_list_3c_type
      TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: ij_list => NULL(), jk_list => NULL()
      INTEGER :: sym
      TYPE(distribution_3d_type) :: dist_3d
      LOGICAL :: owns_dist
   END TYPE

   TYPE neighbor_list_3c_iterator_type
      TYPE(neighbor_list_iterator_p_type), DIMENSION(:), POINTER :: iter_ij => NULL()
      TYPE(neighbor_list_iterator_p_type), DIMENSION(:), POINTER :: iter_jk => NULL()
      INTEGER                                                    :: iter_level
      TYPE(neighbor_list_3c_type)                                :: ijk_nl
      INTEGER, DIMENSION(2)                                      :: bounds_i, bounds_j, bounds_k
   END TYPE

CONTAINS
! **************************************************************************************************
!> \brief Create a 3d distribution
!> \param dist_3d 3d distribution object
!> \param dist1 distribution vector along 1st process grid dimension
!> \param dist2 distribution vector along 2nd process grid dimension
!> \param dist3 distribution vector along 3rd process grid dimension
!> \param nkind ...
!> \param particle_set ...
!> \param mp_comm_3d MPI communicator with a 3d cartesian topology
!> \param own_comm Whether mp_comm_3d should be owned by dist_3d (default false)
! **************************************************************************************************
   SUBROUTINE distribution_3d_create(dist_3d, dist1, dist2, dist3, nkind, particle_set, mp_comm_3d, own_comm)
      TYPE(distribution_3d_type)                         :: dist_3d
      INTEGER, DIMENSION(:), INTENT(IN)                  :: dist1, dist2, dist3
      INTEGER, INTENT(IN)                                :: nkind
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      INTEGER, INTENT(IN)                                :: mp_comm_3d
      LOGICAL, INTENT(IN), OPTIONAL                      :: own_comm

      CHARACTER(len=*), PARAMETER :: routineN = 'distribution_3d_create', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: comm_2d_1, comm_2d_2, handle, nproc_1, &
                                                            nproc_2
      INTEGER, DIMENSION(2)                              :: mp_coor_1, mp_coor_2, mp_dims_1, &
                                                            mp_dims_2

      CALL timeset(routineN, handle)

      IF (PRESENT(own_comm)) THEN
         IF (own_comm) dist_3d%comm_3d = mp_comm_3d
         dist_3d%owns_comm = own_comm
      ELSE
         dist_3d%owns_comm = .FALSE.
      ENDIF

      CALL mp_cart_sub(mp_comm_3d, [.TRUE., .TRUE., .FALSE.], comm_2d_1)
      CALL mp_cart_sub(mp_comm_3d, [.FALSE., .TRUE., .TRUE.], comm_2d_2)

      CALL mp_environ(nproc_1, mp_dims_1, mp_coor_1, comm_2d_1)
      CALL mp_environ(nproc_2, mp_dims_2, mp_coor_2, comm_2d_2)

      CPASSERT(mp_coor_1(2) == mp_coor_2(1))

      CALL distribution_2d_create(dist_3d%dist_2d_1, dist1, dist2, nkind, particle_set, comm_2d_1)
      CALL distribution_2d_create(dist_3d%dist_2d_2, dist2, dist3, nkind, particle_set, comm_2d_2)

      dist_3d%comm_2d_1 = comm_2d_1
      dist_3d%comm_2d_2 = comm_2d_2

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Destroy a 3d distribution
!> \param dist ...
! **************************************************************************************************
   SUBROUTINE distribution_3d_destroy(dist)
      TYPE(distribution_3d_type)                         :: dist

      CHARACTER(len=*), PARAMETER :: routineN = 'distribution_3d_destroy', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)
      CALL distribution_2d_release(dist%dist_2d_1)
      CALL distribution_2d_release(dist%dist_2d_2)
      CALL mp_comm_free(dist%comm_2d_1)
      CALL mp_comm_free(dist%comm_2d_2)
      IF (dist%owns_comm) CALL mp_comm_free(dist%comm_3d)

      NULLIFY (dist%dist_2d_1, dist%dist_2d_2)

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief Create a 2d distribution. This mainly wraps distribution_2d_create
!>        for consistency with distribution_3d_create.
!> \param dist_2d 2d distribution object
!> \param dist1 distribution vector along 1st process grid dimension
!> \param dist2 distribution vector along 2nd process grid dimension
!> \param nkind ...
!> \param particle_set ...
!> \param mp_comm_2d MPI communicator with a 3d cartesian topology
! **************************************************************************************************
   SUBROUTINE distribution_2d_create(dist_2d, dist1, dist2, nkind, particle_set, mp_comm_2d)
      TYPE(distribution_2d_type), POINTER                :: dist_2d
      INTEGER, DIMENSION(:), INTENT(IN)                  :: dist1, dist2
      INTEGER, INTENT(IN)                                :: nkind
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      INTEGER, INTENT(IN)                                :: mp_comm_2d

      CHARACTER(len=*), PARAMETER :: routineN = 'distribution_2d_create', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iatom, ikind, n, natom, nproc
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: nparticle_local_col, nparticle_local_row
      INTEGER, DIMENSION(2)                              :: mp_coor, mp_dims
      INTEGER, DIMENSION(:, :), POINTER                  :: dist1_prv, dist2_prv
      TYPE(cp_1d_i_p_type), DIMENSION(:), POINTER        :: local_particle_col, local_particle_row
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_para_env_type), POINTER                    :: para_env

      NULLIFY (blacs_env, local_particle_col, local_particle_row, para_env)

      CALL timeset(routineN, handle)

      CALL mp_environ(nproc, mp_dims, mp_coor, mp_comm_2d)
      CALL cp_para_env_create(para_env, mp_comm_2d, owns_group=.FALSE.)
      CALL cp_blacs_env_create(blacs_env, para_env, &
                               grid_2d=mp_dims)

      CPASSERT(blacs_env%mepos(1) == mp_coor(1))
      CPASSERT(blacs_env%mepos(2) == mp_coor(2))

      CALL cp_para_env_release(para_env)

      natom = SIZE(particle_set)
      ALLOCATE (dist1_prv(natom, 2), dist2_prv(natom, 2))
      dist1_prv(:, 1) = dist1
      dist2_prv(:, 1) = dist2

      ALLOCATE (local_particle_col(nkind), local_particle_row(nkind))
      ALLOCATE (nparticle_local_row(nkind), nparticle_local_col(nkind))
      nparticle_local_row = 0; nparticle_local_col = 0

      DO iatom = 1, natom
         ikind = particle_set(iatom)%atomic_kind%kind_number

         IF (dist1_prv(iatom, 1) == mp_coor(1)) nparticle_local_row(ikind) = nparticle_local_row(ikind) + 1
         IF (dist2_prv(iatom, 1) == mp_coor(2)) nparticle_local_col(ikind) = nparticle_local_col(ikind) + 1
      END DO

      DO ikind = 1, nkind
         n = nparticle_local_row(ikind)
         ALLOCATE (local_particle_row(ikind)%array(n))

         n = nparticle_local_col(ikind)
         ALLOCATE (local_particle_col(ikind)%array(n))
      END DO

      nparticle_local_row = 0; nparticle_local_col = 0
      DO iatom = 1, natom
         ikind = particle_set(iatom)%atomic_kind%kind_number

         IF (dist1_prv(iatom, 1) == mp_coor(1)) THEN
            nparticle_local_row(ikind) = nparticle_local_row(ikind) + 1
            local_particle_row(ikind)%array(nparticle_local_row(ikind)) = iatom
         END IF
         IF (dist2_prv(iatom, 1) == mp_coor(2)) THEN
            nparticle_local_col(ikind) = nparticle_local_col(ikind) + 1
            local_particle_col(ikind)%array(nparticle_local_col(ikind)) = iatom
         END IF
      END DO

      CALL distribution_2d_create_prv(dist_2d, row_distribution_ptr=dist1_prv, &
                                      col_distribution_ptr=dist2_prv, local_rows_ptr=local_particle_row, &
                                      local_cols_ptr=local_particle_col, blacs_env=blacs_env)

      CALL cp_blacs_env_release(blacs_env)

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief cyclic distribution of weighted elements
!> \param nel ...
!> \param nbin ...
!> \param weights ...
!> \param dist ...
! **************************************************************************************************
   SUBROUTINE cyclic_tensor_dist(nel, nbin, weights, dist)
      INTEGER, INTENT(IN)                                :: nel, nbin
      INTEGER, DIMENSION(nel), INTENT(IN)                :: weights
      INTEGER, DIMENSION(nel), INTENT(OUT)               :: dist

      INTEGER                                            :: ibin, iel, niter
      INTEGER, DIMENSION(nbin)                           :: occup

      occup(:) = 0
      ibin = 0
      DO iel = 1, nel
         niter = 0
         ibin = MOD(ibin + 1, nbin)
         DO WHILE (occup(ibin + 1) + weights(iel) .GE. MAXVAL(occup))
            IF (MINLOC(occup, DIM=1) == ibin + 1) EXIT
            ibin = MOD(ibin + 1, nbin)
            niter = niter + 1
         ENDDO
         dist(iel) = ibin
         occup(ibin + 1) = occup(ibin + 1) + weights(iel)
      ENDDO

   END SUBROUTINE cyclic_tensor_dist

! **************************************************************************************************
!> \brief contiguous distribution of weighted elements
!> \param nel ...
!> \param nbin ...
!> \param weights ...
!> \param limits_start ...
!> \param limits_end ...
!> \param dist ...
! **************************************************************************************************
   SUBROUTINE contiguous_tensor_dist(nel, nbin, weights, limits_start, limits_end, dist)
      INTEGER, INTENT(IN)                                :: nel
      INTEGER, INTENT(INOUT)                             :: nbin
      INTEGER, DIMENSION(nel), INTENT(IN)                :: weights
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT), &
         OPTIONAL                                        :: limits_start, limits_end
      INTEGER, DIMENSION(nel), INTENT(OUT), OPTIONAL     :: dist

      INTEGER                                            :: el_end, el_start, end_weight, ibin, &
                                                            nel_div, nel_rem, nel_split, nel_w, &
                                                            w_partialsum
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: lim_e, lim_s

      ALLOCATE (lim_s(nbin), lim_e(nbin))
      lim_s = 0; lim_e = 0

      nel_w = SUM(weights)
      nel_div = nel_w/nbin
      nel_rem = MOD(nel_w, nbin)

      w_partialsum = 0
      el_end = 0
      end_weight = 0
      DO ibin = 1, nbin
         nel_split = nel_div
         IF (ibin <= nel_rem) THEN
            nel_split = nel_split + 1
         ENDIF
         el_start = el_end + 1
         el_end = el_start
         w_partialsum = w_partialsum + weights(el_end)
         end_weight = end_weight + nel_split
         DO WHILE (w_partialsum < end_weight)
            !IF (ABS(w_partialsum + weights(el_end) - end_weight) > ABS(w_partialsum - end_weight)) EXIT
            el_end = el_end + 1
            w_partialsum = w_partialsum + weights(el_end)
            IF (el_end == nel) EXIT
         ENDDO

         IF (PRESENT(dist)) dist(el_start:el_end) = ibin - 1
         lim_s(ibin) = el_start
         lim_e(ibin) = el_end

         IF (el_end == nel) EXIT
      ENDDO

      IF (PRESENT(limits_start) .AND. PRESENT(limits_end)) THEN
         ALLOCATE (limits_start(ibin)); limits_start(:ibin) = lim_s(:ibin)
         ALLOCATE (limits_end(ibin)); limits_end(:ibin) = lim_e(:ibin)
      ENDIF

      nbin = ibin

   END SUBROUTINE contiguous_tensor_dist

! **************************************************************************************************
!> \brief ...
!> \param t3c Create 3-center tensor with load balanced default distribution.
!> \param dist_1 ...
!> \param dist_2 ...
!> \param dist_3 ...
!> \param pgrid ...
!> \param sizes_1 ...
!> \param sizes_2 ...
!> \param sizes_3 ...
!> \param map1 ...
!> \param map2 ...
!> \param starts_array_block_1 ...
!> \param ends_array_block_1 ...
!> \param starts_array_block_2 ...
!> \param ends_array_block_2 ...
!> \param starts_array_block_3 ...
!> \param ends_array_block_3 ...
!> \param name ...
! **************************************************************************************************
   SUBROUTINE create_3c_tensor(t3c, dist_1, dist_2, dist_3, pgrid, sizes_1, sizes_2, sizes_3, map1, map2, &
                               starts_array_block_1, ends_array_block_1, &
                               starts_array_block_2, ends_array_block_2, &
                               starts_array_block_3, ends_array_block_3, &
                               name)
      TYPE(dbcsr_t_type), INTENT(OUT)                    :: t3c
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: dist_1, dist_2, dist_3
      TYPE(dbcsr_t_pgrid_type), INTENT(IN)               :: pgrid
      INTEGER, DIMENSION(:), INTENT(IN)                  :: sizes_1, sizes_2, sizes_3, map1, map2
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: starts_array_block_1, ends_array_block_1, &
         starts_array_block_2, ends_array_block_2, starts_array_block_3, ends_array_block_3
      CHARACTER(len=*), INTENT(IN)                       :: name

      CHARACTER(len=*), PARAMETER :: routineN = 'create_3c_tensor', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, imem, nbatches_1, nbatches_2, &
                                                            nbatches_3, size_1, size_2, size_3, &
                                                            size_cut
      INTEGER, DIMENSION(3)                              :: pcoord, pdims
      LOGICAL                                            :: mem_aware_1, mem_aware_2, mem_aware_3
      TYPE(dbcsr_t_distribution_type)                    :: dist

      CALL timeset(routineN, handle)

      mem_aware_1 = PRESENT(starts_array_block_1) .AND. PRESENT(ends_array_block_1)
      mem_aware_2 = PRESENT(starts_array_block_2) .AND. PRESENT(ends_array_block_2)
      mem_aware_3 = PRESENT(starts_array_block_3) .AND. PRESENT(ends_array_block_3)

      IF (mem_aware_1) nbatches_1 = SIZE(starts_array_block_1)
      IF (mem_aware_2) nbatches_2 = SIZE(starts_array_block_2)
      IF (mem_aware_3) nbatches_3 = SIZE(starts_array_block_3)

      CALL dbcsr_t_mp_environ_pgrid(pgrid, pdims, pcoord)

      size_1 = SIZE(sizes_1)
      size_2 = SIZE(sizes_2)
      size_3 = SIZE(sizes_3)

      ALLOCATE (dist_1(size_1))
      ALLOCATE (dist_2(size_2))
      ALLOCATE (dist_3(size_3))

      IF (mem_aware_1) THEN
         DO imem = 1, nbatches_1
            size_cut = ends_array_block_1(imem) - starts_array_block_1(imem) + 1
            CALL cyclic_tensor_dist( &
               size_cut, pdims(1), &
               sizes_1(starts_array_block_1(imem):ends_array_block_1(imem)), &
               dist_1(starts_array_block_1(imem):ends_array_block_1(imem)))
         ENDDO
      ELSE
         CALL cyclic_tensor_dist(size_1, pdims(1), sizes_1, dist_1)
      ENDIF

      IF (mem_aware_2) THEN
         DO imem = 1, nbatches_2
            size_cut = ends_array_block_2(imem) - starts_array_block_2(imem) + 1
            CALL cyclic_tensor_dist( &
               size_cut, pdims(2), &
               sizes_2(starts_array_block_2(imem):ends_array_block_2(imem)), &
               dist_2(starts_array_block_2(imem):ends_array_block_2(imem)))
         ENDDO
      ELSE
         CALL cyclic_tensor_dist(size_2, pdims(2), sizes_2, dist_2)
      ENDIF

      IF (mem_aware_3) THEN
         DO imem = 1, nbatches_3
            size_cut = ends_array_block_3(imem) - starts_array_block_3(imem) + 1
            CALL cyclic_tensor_dist( &
               size_cut, pdims(3), &
               sizes_3(starts_array_block_3(imem):ends_array_block_3(imem)), &
               dist_3(starts_array_block_3(imem):ends_array_block_3(imem)))
         ENDDO
      ELSE
         CALL cyclic_tensor_dist(size_3, pdims(3), sizes_3, dist_3)
      ENDIF

      CALL dbcsr_t_distribution_new(dist, pgrid, dist_1, dist_2, dist_3)
      CALL dbcsr_t_create(t3c, name, dist, map1, map2, dbcsr_type_real_8, sizes_1, &
                          sizes_2, sizes_3)
      CALL dbcsr_t_distribution_destroy(dist)

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param t2c ...
!> \param dist_1 ...
!> \param dist_2 ...
!> \param pgrid ...
!> \param sizes_1 ...
!> \param sizes_2 ...
!> \param order ...
!> \param name ...
! **************************************************************************************************
   SUBROUTINE create_2c_tensor(t2c, dist_1, dist_2, pgrid, sizes_1, sizes_2, order, name)
      TYPE(dbcsr_t_type), INTENT(OUT)                    :: t2c
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: dist_1, dist_2
      TYPE(dbcsr_t_pgrid_type), INTENT(IN)               :: pgrid
      INTEGER, DIMENSION(:), INTENT(IN)                  :: sizes_1, sizes_2
      INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL        :: order
      CHARACTER(len=*), INTENT(IN)                       :: name

      CHARACTER(len=*), PARAMETER :: routineN = 'create_2c_tensor', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, size_1, size_2
      INTEGER, DIMENSION(2)                              :: order_in, pcoord, pdims
      TYPE(dbcsr_t_distribution_type)                    :: dist

      CALL timeset(routineN, handle)

      IF (PRESENT(order)) THEN
         order_in = order
      ELSE
         order_in = [1, 2]
      ENDIF

      CALL dbcsr_t_mp_environ_pgrid(pgrid, pdims, pcoord)

      size_1 = SIZE(sizes_1)
      size_2 = SIZE(sizes_2)

      ALLOCATE (dist_1(size_1))
      ALLOCATE (dist_2(size_2))

      CALL cyclic_tensor_dist(size_1, pdims(1), sizes_1, dist_1)
      CALL cyclic_tensor_dist(size_2, pdims(2), sizes_2, dist_2)

      CALL dbcsr_t_distribution_new(dist, pgrid, dist_1, dist_2)
      CALL dbcsr_t_create(t2c, name, dist, [order_in(1)], [order_in(2)], dbcsr_type_real_8, sizes_1, sizes_2)
      CALL dbcsr_t_distribution_destroy(dist)

      CALL timestop(handle)
   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param blk_sizes ...
!> \param blk_sizes_split ...
!> \param max_size ...
! **************************************************************************************************
   SUBROUTINE split_block_sizes(blk_sizes, blk_sizes_split, max_size)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: blk_sizes
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(OUT)    :: blk_sizes_split
      INTEGER, INTENT(IN)                                :: max_size

      INTEGER                                            :: blk_remainder, i, isplit, isplit_sum, &
                                                            nsplit

      isplit_sum = 0
      DO i = 1, SIZE(blk_sizes)
         nsplit = (blk_sizes(i) + max_size - 1)/max_size
         isplit_sum = isplit_sum + nsplit
      ENDDO

      ALLOCATE (blk_sizes_split(isplit_sum))

      isplit_sum = 0
      DO i = 1, SIZE(blk_sizes)
         nsplit = (blk_sizes(i) + max_size - 1)/max_size
         blk_remainder = blk_sizes(i)
         DO isplit = 1, nsplit
            isplit_sum = isplit_sum + 1
            blk_sizes_split(isplit_sum) = MIN(max_size, blk_remainder)
            blk_remainder = blk_remainder - max_size
         ENDDO
      ENDDO

   END SUBROUTINE split_block_sizes

END MODULE
