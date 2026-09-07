!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Assembly of complex k-point operators from real-space DBCSR matrices.
!>        The output is a complex full matrix on the distribution of the k-point group.
!>        The module splits this work into a start step and a finish step.
!>        One process can start the MPI transfers of many operators before it finishes
!>        any of them. The STANDARD k-point driver uses this order.
!>        The routines start, finish and get are collective over the
!>        environment-wide communicator. Every rank must call the same
!>        sequence of service calls.
!>        A context references its host kpoint environment without owning it,
!>        so it must not outlive that environment. The calls are not thread
!>        safe: the caller must serialize them, because their collective
!>        communication tolerates no concurrent service call. A separate
!>        context per thread only keeps the scratch from aliasing.
! **************************************************************************************************
MODULE qs_kpoint_operators
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_scale_and_add_fm
   USE cp_cfm_types,                    ONLY: cp_cfm_to_fm,&
                                              cp_cfm_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_create, dbcsr_deallocate_matrix, dbcsr_desymmetrize, dbcsr_p_type, dbcsr_set, &
        dbcsr_type, dbcsr_type_antisymmetric, dbcsr_type_no_symmetry, dbcsr_type_symmetric
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr_bc
   USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                              fm_pool_create_fm,&
                                              fm_pool_give_back_fm
   USE cp_fm_struct,                    ONLY: cp_fm_struct_equivalent,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: copy_info_type,&
                                              cp_fm_cleanup_copy_general,&
                                              cp_fm_create,&
                                              cp_fm_finish_copy_general,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_start_copy_general,&
                                              cp_fm_type
   USE kinds,                           ONLY: dp
   USE kpoint_methods,                  ONLY: rskp_grid_type,&
                                              rskp_transform,&
                                              rskp_transform_grid_extract
   USE kpoint_types,                    ONLY: kpoint_type
   USE mathconstants,                   ONLY: gaussi,&
                                              z_one,&
                                              z_zero
   USE message_passing,                 ONLY: mp_para_env_type
   USE qs_matrix_pools,                 ONLY: mpools_get,&
                                              qs_matrix_pools_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_kpoint_operators'

   PUBLIC :: kpoint_operator_cfm_to_mo
   PUBLIC :: kpoint_operator_context_create, kpoint_operator_context_release
   PUBLIC :: kpoint_operator_start, kpoint_operator_finish
   PUBLIC :: kpoint_operator_get

   INTEGER, PARAMETER, PUBLIC :: kpoint_slot_ks = 1, &
                                 kpoint_slot_s = 2, &
                                 kpoint_slot_t = 3

   !> Row of rsmat that carries the spin-free operator image. S and T arrive
   !> with a single spin row. Callers launch their one instance under this
   !> spin key and read this row.
   INTEGER, PARAMETER, PUBLIC :: kpoint_spin_free = 1

   !> Number of defined slots; grows only inside this module.
   INTEGER, PARAMETER, PRIVATE :: kpoint_num_slots = kpoint_slot_t - kpoint_slot_ks + 1

   !> Part indices of one transfer: its real and its imaginary half.
   INTEGER, PARAMETER, PRIVATE :: kpoint_part_re = 1, &
                                  kpoint_part_im = 2

   !> Instance states of one (k point, spin, slot) transfer.
   !> State transitions happen only inside the service routines.
   INTEGER, PARAMETER, PRIVATE :: kpoint_slot_idle = 0, &
                                  kpoint_slot_in_flight = 1

! **************************************************************************************************
!> \brief Caller-owned assembly context for one driver call.
!>        Created and released in the same driver routine and passed to helpers
!>        as an argument; never stored on the host kpoint_type.
!>        Components are private and touched only by this module.
! **************************************************************************************************
   TYPE, PUBLIC :: kpoint_operator_context_type
      PRIVATE
      ! Instance spin space of the owning driver call: the bound of the spin
      ! loop that start and finish key. Declared once at creation. It never
      ! comes from the operators: the first matrix dimension of a spin-free
      ! operator such as S counts real-space blocks, not spins, so callers
      ! address operator rows through matrix_row instead.
      INTEGER                                            :: nspin = 0
      ! Host references, read fresh at every use: sab_nl and cell_to_index
      ! are reallocated at every energy evaluation.
      TYPE(kpoint_type), POINTER                          :: kpoints => NULL()
      TYPE(cp_fm_struct_type), POINTER                    :: ao_ao_fmstruct => NULL()
      TYPE(qs_matrix_pools_type), POINTER                 :: group_mpools => NULL()
      TYPE(cp_fm_pool_p_type), DIMENSION(:), POINTER      :: group_fm_pools => NULL()
      ! Identity of the neighbor lists the work trio was built from; every
      ! buffer use asserts it.
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                          :: trio_sab => NULL()
      ! Owned scratch, created on first use and released with the context.
      TYPE(dbcsr_type), POINTER                           :: rmatrix => NULL()
      TYPE(dbcsr_type), POINTER                           :: cmatrix => NULL()
      TYPE(dbcsr_type), POINTER                           :: tmpmat => NULL()
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)         :: fmwork
      TYPE(cp_fm_type), ALLOCATABLE                       :: fmlocal
      TYPE(cp_fm_type)                                    :: fmdummy
      ! Bookkeeping. info is keyed (local k point x group, spin, slot, part),
      ! part kpoint_part_re holding the real half and kpoint_part_im holding
      ! the imaginary half. Status is keyed (local k point, spin, slot)
      ! because start and finish treat all groups of one instance
      ! atomically. Allocated on first start.
      TYPE(copy_info_type), ALLOCATABLE, DIMENSION(:, :, :, :) :: info
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)            :: status
   END TYPE kpoint_operator_context_type

CONTAINS

! **************************************************************************************************
!> \brief Transform one real-space operator to one k point and densify the
!>        result: the real part into fm_re and the imaginary part into fm_im.
!>        Densify means: copy the blocks of a DBCSR matrix into a full matrix.
!>        Stateless: fills caller-provided scratch.
!> \param rsmat real-space image matrices of the operator
!> \param matrix_row row of rsmat to transform (kpoint_spin_free for
!>        spin-free matrices such as S)
!> \param ik global k-point index
!> \param xkp coordinates of this k point
!> \param cell_to_index mapping of cell indices to real-space index
!> \param sab_nl neighbor lists defining the real-space sparsity
!> \param grid prepared reciprocal-grid cache of the caller. Required when use_grid is set
!> \param use_grid extract from grid instead of a direct phase sum
!> \param rmatrix symmetric DBCSR workspace holding the real part of the k-point matrix
!> \param cmatrix antisymmetric DBCSR workspace holding the imaginary part of the k-point matrix
!> \param tmpmat unsymmetric DBCSR workspace
!> \param fm_re full-matrix workspace holding the real part
!> \param fm_im full-matrix workspace holding the imaginary part
! **************************************************************************************************
   SUBROUTINE kpoint_operator_densify(rsmat, matrix_row, ik, xkp, cell_to_index, sab_nl, grid, use_grid, &
                                      rmatrix, cmatrix, tmpmat, fm_re, fm_im)

      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rsmat
      INTEGER, INTENT(IN)                                :: matrix_row, ik
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: xkp
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(rskp_grid_type), INTENT(IN), OPTIONAL         :: grid
      LOGICAL, INTENT(IN)                                :: use_grid
      TYPE(dbcsr_type)                                   :: rmatrix, cmatrix, tmpmat
      TYPE(cp_fm_type)                                   :: fm_re, fm_im

      CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_operator_densify'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (use_grid) THEN
         CPASSERT(PRESENT(grid))
         CALL rskp_transform_grid_extract(grid, ik, rmatrix, cmatrix)
      ELSE
         CALL dbcsr_set(rmatrix, 0.0_dp)
         CALL dbcsr_set(cmatrix, 0.0_dp)
         CALL rskp_transform(rmatrix=rmatrix, cmatrix=cmatrix, rsmat=rsmat, ispin=matrix_row, &
                             xkp=xkp, cell_to_index=cell_to_index, sab_nl=sab_nl)
      END IF
      CALL dbcsr_desymmetrize(rmatrix, tmpmat)
      CALL copy_dbcsr_to_fm(tmpmat, fm_re)
      CALL dbcsr_desymmetrize(cmatrix, tmpmat)
      CALL copy_dbcsr_to_fm(tmpmat, fm_im)

      CALL timestop(handle)

   END SUBROUTINE kpoint_operator_densify

! **************************************************************************************************
!> \brief Transform one real-space operator row to one k point, densify the
!>        result, and start the transfer of both parts: into target_re and
!>        target_im on the group that owns this k point, into fmdummy on
!>        every other group so the collective calls stay balanced.
!>        Densify means: copy the blocks of a DBCSR matrix into a full matrix.
!>        Each transfer receives its data in a private buffer. For this reason,
!>        several launches can use the same target matrix.
!>        The caller owns the transfer states: finish or clean up every started
!>        copy through info_re and info_im.
!> \param rsmat real-space image matrices of the operator
!> \param matrix_row row of rsmat to transform (kpoint_spin_free for
!>        spin-free matrices such as S)
!> \param ik global k-point index
!> \param xkp coordinates of this k point
!> \param cell_to_index mapping of cell indices to real-space index
!> \param sab_nl neighbor lists defining the real-space sparsity
!> \param grid prepared reciprocal-grid cache of the caller. Required when use_grid is set
!> \param use_grid extract from grid instead of a direct phase sum
!> \param rmatrix symmetric DBCSR workspace holding the real part of the k-point matrix
!> \param cmatrix antisymmetric DBCSR workspace holding the imaginary part of the k-point matrix
!> \param tmpmat unsymmetric DBCSR workspace
!> \param fm_re full-matrix workspace of the global communicator holding the real part
!> \param fm_im full-matrix workspace of the global communicator holding the imaginary part
!> \param target_re transfer target of the real part on the group that owns this k point
!> \param target_im transfer target of the imaginary part on the owning group
!> \param fmdummy transfer target on the other groups. This keeps the collective
!>        calls balanced
!> \param my_kpgrp true when the calling group owns this k point
!> \param para_env communicator that spans the source and the target distributions
!> \param info_re transfer state of the real part
!> \param info_im transfer state of the imaginary part
! **************************************************************************************************
   SUBROUTINE kpoint_operator_launch(rsmat, matrix_row, ik, xkp, cell_to_index, sab_nl, &
                                     grid, use_grid, rmatrix, cmatrix, tmpmat, fm_re, fm_im, &
                                     target_re, target_im, fmdummy, my_kpgrp, para_env, &
                                     info_re, info_im)

      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rsmat
      INTEGER, INTENT(IN)                                :: matrix_row, ik
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: xkp
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(rskp_grid_type), INTENT(IN), OPTIONAL         :: grid
      LOGICAL, INTENT(IN)                                :: use_grid
      TYPE(dbcsr_type)                                   :: rmatrix, cmatrix, tmpmat
      TYPE(cp_fm_type)                                   :: fm_re, fm_im, target_re, target_im, &
                                                            fmdummy
      LOGICAL, INTENT(IN)                                :: my_kpgrp
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(copy_info_type)                               :: info_re, info_im

      CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_operator_launch'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL kpoint_operator_densify(rsmat, matrix_row, ik, xkp, cell_to_index, sab_nl, grid, use_grid, &
                                   rmatrix, cmatrix, tmpmat, fm_re, fm_im)

      IF (my_kpgrp) THEN
         CALL cp_fm_start_copy_general(fm_re, target_re, para_env, info_re)
         CALL cp_fm_start_copy_general(fm_im, target_im, para_env, info_im)
      ELSE
         CALL cp_fm_start_copy_general(fm_re, fmdummy, para_env, info_re)
         CALL cp_fm_start_copy_general(fm_im, fmdummy, para_env, info_im)
      END IF

      CALL timestop(handle)

   END SUBROUTINE kpoint_operator_launch

! **************************************************************************************************
!> \brief Finish one transfer that kpoint_operator_launch started and add the received
!>        matrix to a complex matrix.
!>        The finish of the real part overwrites cmat.
!>        The finish of the imaginary part adds to cmat.
!>        Finish the real part first. A real-part finish after an imaginary-part
!>        finish erases the imaginary contribution.
!>        The routine consumes fmlocal before it returns, so the next finish can
!>        reuse the same matrix.
!> \param fmlocal transfer target filled by the finished copy
!> \param info transfer state of one part
!> \param cmat complex matrix accumulating the operator
!> \param imaginary_part add fmlocal as the imaginary part instead of the real part
! **************************************************************************************************
   SUBROUTINE kpoint_operator_raw_finish_cfm(fmlocal, info, cmat, imaginary_part)

      TYPE(cp_fm_type)                                   :: fmlocal
      TYPE(copy_info_type)                               :: info
      TYPE(cp_cfm_type)                                  :: cmat
      LOGICAL, INTENT(IN)                                :: imaginary_part

      CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_operator_raw_finish_cfm'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CALL cp_fm_finish_copy_general(fmlocal, info)
      IF (imaginary_part) THEN
         CALL cp_cfm_scale_and_add_fm(z_one, cmat, gaussi, fmlocal)
      ELSE
         CALL cp_cfm_scale_and_add_fm(z_zero, cmat, z_one, fmlocal)
      END IF

      CALL timestop(handle)

   END SUBROUTINE kpoint_operator_raw_finish_cfm

! **************************************************************************************************
!> \brief Split a complex matrix into the real and imaginary MO sets and copy
!>        the eigenvalues to the imaginary set.
!> \param cmat complex MO coefficients
!> \param mo_re real part MO set
!> \param mo_im imaginary part MO set
!> \param eigenvalues eigenvalues of mo_re, copied to mo_im
! **************************************************************************************************
   SUBROUTINE kpoint_operator_cfm_to_mo(cmat, mo_re, mo_im, eigenvalues)

      TYPE(cp_cfm_type)                                  :: cmat
      TYPE(mo_set_type)                                  :: mo_re, mo_im
      REAL(KIND=dp), DIMENSION(:)                        :: eigenvalues

      CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_operator_cfm_to_mo'

      INTEGER                                            :: handle
      TYPE(cp_fm_type), POINTER                          :: imos, rmos

      CALL timeset(routineN, handle)

      CALL get_mo_set(mo_re, mo_coeff=rmos)
      CALL get_mo_set(mo_im, mo_coeff=imos)
      mo_im%eigenvalues = eigenvalues
      CALL cp_cfm_to_fm(cmat, rmos, imos)

      CALL timestop(handle)

   END SUBROUTINE kpoint_operator_cfm_to_mo

! **************************************************************************************************
!> \brief Create an empty operator context.
!>        The context owns no buffers after this call; they are created on first
!>        use and released with the context.
!> \param ctx context to initialize
!> \param kpoints host k-point environment, referenced, not copied
!> \param ao_ao_fm any square full matrix on the environment-wide communicator;
!>        its matrix struct is the template for the context work buffers,
!>        which densify on the source side of transfers that cross into the
!>        k-point groups
!> \param nspin instance spin space of the owning driver call: the bound of
!>        the spin loop its start and finish calls iterate. Fixed for the
!>        lifetime of the context
! **************************************************************************************************
   SUBROUTINE kpoint_operator_context_create(ctx, kpoints, ao_ao_fm, nspin)

      TYPE(kpoint_operator_context_type), INTENT(OUT)    :: ctx
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(cp_fm_type), INTENT(IN)                       :: ao_ao_fm
      INTEGER, INTENT(IN)                                :: nspin

      CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_operator_context_create'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CPASSERT(nspin >= 1)
      ctx%nspin = nspin
      CPASSERT(ASSOCIATED(kpoints))
      ctx%kpoints => kpoints
      CALL cp_fm_get_info(ao_ao_fm, matrix_struct=ctx%ao_ao_fmstruct)
      CPASSERT(ASSOCIATED(ctx%ao_ao_fmstruct))
      ! The template must live on the transfer world. A template from
      ! another communicator fails during the transfers, far from this call.
      CPASSERT(ASSOCIATED(ctx%ao_ao_fmstruct%para_env, kpoints%blacs_env_all%para_env))

      CALL timestop(handle)

   END SUBROUTINE kpoint_operator_context_create

! **************************************************************************************************
!> \brief Release the context together with every buffer it created.
!>        Call once, at the exit of the routine that created the context.
!>        Every started instance must be idle at this point. The call aborts
!>        when a transfer is still in flight, because nothing would ever wait
!>        for its MPI requests.
!> \param ctx context to release
! **************************************************************************************************
   SUBROUTINE kpoint_operator_context_release(ctx)

      TYPE(kpoint_operator_context_type), INTENT(INOUT)  :: ctx

      CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_operator_context_release'

      INTEGER                                            :: handle, iwork

      CALL timeset(routineN, handle)

      ! Every started instance must have been finished: an unfinished start
      ! holds MPI requests that nothing would ever wait for.
      IF (ALLOCATED(ctx%status)) THEN
         CPASSERT(ALL(ctx%status == kpoint_slot_idle))
      END IF

      IF (ASSOCIATED(ctx%rmatrix)) CALL dbcsr_deallocate_matrix(ctx%rmatrix)
      IF (ASSOCIATED(ctx%cmatrix)) CALL dbcsr_deallocate_matrix(ctx%cmatrix)
      IF (ASSOCIATED(ctx%tmpmat)) CALL dbcsr_deallocate_matrix(ctx%tmpmat)
      NULLIFY (ctx%rmatrix, ctx%cmatrix, ctx%tmpmat, ctx%trio_sab)
      IF (ALLOCATED(ctx%fmwork)) THEN
         DO iwork = 1, SIZE(ctx%fmwork)
            CALL cp_fm_release(ctx%fmwork(iwork))
         END DO
         DEALLOCATE (ctx%fmwork)
      END IF
      IF (ALLOCATED(ctx%fmlocal)) THEN
         IF (ASSOCIATED(ctx%group_fm_pools)) THEN
            CALL fm_pool_give_back_fm(ctx%group_fm_pools(1)%pool, ctx%fmlocal)
         ELSE
            CALL cp_fm_release(ctx%fmlocal)
         END IF
         DEALLOCATE (ctx%fmlocal)
      END IF
      IF (ALLOCATED(ctx%info)) DEALLOCATE (ctx%info)
      IF (ALLOCATED(ctx%status)) DEALLOCATE (ctx%status)
      NULLIFY (ctx%kpoints, ctx%ao_ao_fmstruct, ctx%group_mpools, ctx%group_fm_pools)

      CALL timestop(handle)

   END SUBROUTINE kpoint_operator_context_release

! **************************************************************************************************
!> \brief Create the context buffers that the first start needs: the DBCSR work
!>        trio with blocks from the current neighbor lists, two densify buffers
!>        on the environment-wide struct, the group transfer target from the
!>        group pools, and the per-instance bookkeeping.
!>        The buffers assume the context lifetime of one driver call: the
!>        host neighbor lists do not change within it, and every buffer use
!>        asserts this identity. Buffers from an earlier call stay in place.
!>        Called from the service routines only.
!> \param ctx context that receives the buffers
!> \param rsmat real-space image matrices of one operator. The call uses the
!>        first matrix as the DBCSR template of the work trio. It reads no
!>        matrix data. The bookkeeping spin space comes from the value that
!>        context creation stored, never from the row count of rsmat
! **************************************************************************************************
   SUBROUTINE kpoint_operator_ensure(ctx, rsmat)

      TYPE(kpoint_operator_context_type), INTENT(INOUT)  :: ctx
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rsmat

      CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_operator_ensure'

      INTEGER                                            :: handle, kplocal, ngrp

      CALL timeset(routineN, handle)

      kplocal = ctx%kpoints%kp_range(2) - ctx%kpoints%kp_range(1) + 1
      ngrp = ctx%kpoints%nkp_groups

      ! The trio must still match the host neighbor lists: they are
      ! reallocated at every energy evaluation, while the context spans one
      ! driver call within it.
      IF (ASSOCIATED(ctx%rmatrix)) THEN
         CPASSERT(ASSOCIATED(ctx%trio_sab, ctx%kpoints%sab_nl))
      END IF
      IF (.NOT. ASSOCIATED(ctx%rmatrix)) THEN
         ALLOCATE (ctx%rmatrix, ctx%cmatrix, ctx%tmpmat)
         CALL dbcsr_create(ctx%rmatrix, template=rsmat(1, 1)%matrix, matrix_type=dbcsr_type_symmetric)
         CALL dbcsr_create(ctx%cmatrix, template=rsmat(1, 1)%matrix, matrix_type=dbcsr_type_antisymmetric)
         CALL dbcsr_create(ctx%tmpmat, template=rsmat(1, 1)%matrix, matrix_type=dbcsr_type_no_symmetry)
         CALL cp_dbcsr_alloc_block_from_nbl(ctx%rmatrix, ctx%kpoints%sab_nl)
         CALL cp_dbcsr_alloc_block_from_nbl(ctx%cmatrix, ctx%kpoints%sab_nl)
         ctx%trio_sab => ctx%kpoints%sab_nl
      END IF

      IF (.NOT. ALLOCATED(ctx%fmwork)) THEN
         ALLOCATE (ctx%fmwork(kpoint_part_im))
         CALL cp_fm_create(ctx%fmwork(kpoint_part_re), ctx%ao_ao_fmstruct, &
                           name="KPOINT-OPERATOR-WORK-RE")
         CALL cp_fm_create(ctx%fmwork(kpoint_part_im), ctx%ao_ao_fmstruct, &
                           name="KPOINT-OPERATOR-WORK-IM")
      END IF

      IF (.NOT. ALLOCATED(ctx%fmlocal)) THEN
         ! The group pool owns the storage; the release returns the matrix
         ! to the pool instead of destroying it.
         ctx%group_mpools => ctx%kpoints%mpools
         CALL mpools_get(ctx%group_mpools, ao_ao_fm_pools=ctx%group_fm_pools)
         ALLOCATE (ctx%fmlocal)
         CALL fm_pool_create_fm(ctx%group_fm_pools(1)%pool, ctx%fmlocal)
      END IF

      ! The bookkeeping spin space is the instance spin space declared at
      ! context creation. The operator row count never enters the sizing:
      ! spin-free operators such as S carry their real-space blocks in the
      ! first matrix dimension, and the row choice is made through
      ! matrix_row, so a fixed spin space serves every operator of the call.
      IF (.NOT. ALLOCATED(ctx%status)) THEN
         ALLOCATE (ctx%status(kplocal, ctx%nspin, kpoint_num_slots))
         ctx%status = kpoint_slot_idle
         ! One copy_info per in-flight part: start takes it INTENT(OUT),
         ! so two concurrent starts cannot share one instance.
         ALLOCATE (ctx%info(kplocal*ngrp, ctx%nspin, kpoint_num_slots, kpoint_part_im))
      END IF

      CALL timestop(handle)

   END SUBROUTINE kpoint_operator_ensure

! **************************************************************************************************
!> \brief Start the assembly of one operator for one local k point and spin and
!>        launch its transfers to every k-point group. The caller index ikp is
!>        local; the walk over groups happens here, so no caller can express a
!>        group-local loop. Finish every started instance with
!>        kpoint_operator_finish.
!>        Collective over the environment-wide communicator: every rank must
!>        call the same sequence of service calls.
!> \param ctx context created by kpoint_operator_context_create
!> \param ikp local k-point index within kp_range
!> \param ispin instance spin: keys the bookkeeping
!> \param rsmat real-space image matrices of the operator, read at start
!>        time only. The matching finish does not read them
!> \param slot destination slot of the operator (one of kpoint_slot_ks/s/t);
!>        the slot must be idle
!> \param grid prepared reciprocal-grid cache of the caller
!> \param use_grid extract from grid instead of a direct phase sum
!> \param matrix_row row of rsmat to transform. Spin-free operators such
!>        as S and T must pass kpoint_spin_free: their first matrix
!>        dimension counts real-space blocks, not spins, so passing the
!>        instance spin of a spin loop reads a different block and no size
!>        check can catch it
! **************************************************************************************************
   SUBROUTINE kpoint_operator_start(ctx, ikp, ispin, rsmat, slot, grid, use_grid, matrix_row)

      TYPE(kpoint_operator_context_type), INTENT(INOUT)  :: ctx
      INTEGER, INTENT(IN)                                :: ikp, ispin
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rsmat
      INTEGER, INTENT(IN)                                :: slot
      TYPE(rskp_grid_type), INTENT(IN), OPTIONAL         :: grid
      LOGICAL, INTENT(IN)                                :: use_grid
      INTEGER, INTENT(IN)                                :: matrix_row

      CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_operator_start'

      INTEGER                                            :: handle, igroup, ik, indx, kplocal
      INTEGER, DIMENSION(2)                              :: kp_range
      INTEGER, DIMENSION(:, :), POINTER                  :: kp_dist
      LOGICAL                                            :: my_kpgrp
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(ctx%kpoints))
      CPASSERT(ASSOCIATED(ctx%kpoints%sab_nl))
      ! The instance spin is bounded by the spin space declared at creation;
      ! the matrix row by this operator's own rows. The two bounds differ for
      ! spin-free operators such as S, whose first dimension counts
      ! real-space blocks.
      CPASSERT(ispin >= 1 .AND. ispin <= ctx%nspin)
      CPASSERT(matrix_row >= 1 .AND. matrix_row <= SIZE(rsmat, 1))
      kplocal = ctx%kpoints%kp_range(2) - ctx%kpoints%kp_range(1) + 1
      CPASSERT(ikp >= 1 .AND. ikp <= kplocal)
      CPASSERT(slot >= kpoint_slot_ks .AND. slot <= kpoint_slot_t)
      IF (ALLOCATED(ctx%status)) THEN
         CPASSERT(ctx%status(ikp, ispin, slot) == kpoint_slot_idle)
      END IF

      CALL kpoint_operator_ensure(ctx, rsmat)

      kp_range = ctx%kpoints%kp_range
      kp_dist => ctx%kpoints%kp_dist
      para_env => ctx%kpoints%blacs_env_all%para_env
      DO igroup = 1, ctx%kpoints%nkp_groups
         ik = kp_dist(1, igroup) + ikp - 1
         my_kpgrp = (ik >= kp_range(1) .AND. ik <= kp_range(2))
         indx = ikp + (igroup - 1)*kplocal
         CALL kpoint_operator_launch(rsmat, matrix_row, ik, ctx%kpoints%xkp(1:3, ik), &
                                     ctx%kpoints%cell_to_index, ctx%kpoints%sab_nl, grid, use_grid, &
                                     ctx%rmatrix, ctx%cmatrix, ctx%tmpmat, &
                                     ctx%fmwork(kpoint_part_re), ctx%fmwork(kpoint_part_im), &
                                     target_re=ctx%fmlocal, target_im=ctx%fmlocal, fmdummy=ctx%fmdummy, &
                                     my_kpgrp=my_kpgrp, para_env=para_env, &
                                     info_re=ctx%info(indx, ispin, slot, kpoint_part_re), &
                                     info_im=ctx%info(indx, ispin, slot, kpoint_part_im))
      END DO
      ctx%status(ikp, ispin, slot) = kpoint_slot_in_flight

      CALL timestop(handle)

   END SUBROUTINE kpoint_operator_start

! **************************************************************************************************
!> \brief Finish the transfers of one started instance: on the group that owns
!>        the k point, merge the received parts into cmat (real part first,
!>        overwrite then add); on every group, reclaim the send side.
!>        The slot returns to idle.
!>        Collective over the environment-wide communicator: every rank must
!>        call the same sequence of service calls.
!> \param ctx context that started the instance
!> \param ikp local k-point index of the start call
!> \param ispin spin component of the start call
!> \param slot slot of the start call; the slot must be in flight
!> \param cmat caller-allocated complex matrix that receives the operator.
!>        It must share the distribution of the k-point group
! **************************************************************************************************
   SUBROUTINE kpoint_operator_finish(ctx, ikp, ispin, slot, cmat)

      TYPE(kpoint_operator_context_type), INTENT(INOUT)  :: ctx
      INTEGER, INTENT(IN)                                :: ikp, ispin, slot
      TYPE(cp_cfm_type)                                  :: cmat

      CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_operator_finish'

      INTEGER                                            :: handle, igroup, ik, indx, kplocal
      INTEGER, DIMENSION(2)                              :: kp_range
      INTEGER, DIMENSION(:, :), POINTER                  :: kp_dist
      LOGICAL                                            :: my_kpgrp

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(ctx%kpoints))
      CPASSERT(ALLOCATED(ctx%status))
      CPASSERT(ctx%status(ikp, ispin, slot) == kpoint_slot_in_flight)
      ! The merge writes cmat in place. A foreign distribution would fail
      ! inside the BLACS calls, far from this argument.
      CPASSERT(cp_fm_struct_equivalent(cmat%matrix_struct, ctx%fmlocal%matrix_struct))

      kplocal = ctx%kpoints%kp_range(2) - ctx%kpoints%kp_range(1) + 1
      kp_range = ctx%kpoints%kp_range
      kp_dist => ctx%kpoints%kp_dist
      DO igroup = 1, ctx%kpoints%nkp_groups
         ik = kp_dist(1, igroup) + ikp - 1
         my_kpgrp = (ik >= kp_range(1) .AND. ik <= kp_range(2))
         indx = ikp + (igroup - 1)*kplocal
         IF (my_kpgrp) THEN
            CALL kpoint_operator_raw_finish_cfm(ctx%fmlocal, ctx%info(indx, ispin, slot, kpoint_part_re), &
                                                cmat, .FALSE.)
            CALL kpoint_operator_raw_finish_cfm(ctx%fmlocal, ctx%info(indx, ispin, slot, kpoint_part_im), &
                                                cmat, .TRUE.)
         END IF
         CALL cp_fm_cleanup_copy_general(ctx%info(indx, ispin, slot, kpoint_part_re))
         CALL cp_fm_cleanup_copy_general(ctx%info(indx, ispin, slot, kpoint_part_im))
      END DO
      ctx%status(ikp, ispin, slot) = kpoint_slot_idle

      CALL timestop(handle)

   END SUBROUTINE kpoint_operator_finish

! **************************************************************************************************
!> \brief Assemble one operator for one local k point and spin in one call:
!>        run the single-group transfer of the service for every group and
!>        finish it immediately, through transfer states local to this call.
!>        The persistent slot machinery stays with start and finish.
!>        Exactly one output flavor must be present: cmat receives the complex
!>        operator (real part overwrites, imaginary part adds) through the
!>        context buffer; the fm pair receives the two parts separately,
!>        transferred directly into fm_re and fm_im on the owning group (any
!>        distribution: the general copy redistributes); or the dbcsr pair
!>        receives the two parts as block-cyclic matrices built by the service,
!>        replacing the caller-side copy_fm_to_dbcsr_bc detour. The dbcsr
!>        matrices must be freshly initialized (dbcsr_init_p) and the cmat
!>        flavor must share the distribution of the k-point group.
!>        Collective over the environment-wide communicator: every rank must
!>        call the same sequence of service calls.
!> \param ctx context created by kpoint_operator_context_create
!> \param ikp local k-point index within kp_range
!> \param ispin instance spin: keys the bookkeeping
!> \param rsmat real-space image matrices of the operator
!> \param cmat complex matrix output
!> \param fm_re real-part matrix output
!> \param fm_im imaginary-part matrix output
!> \param mat_re real-part block-cyclic dbcsr output
!> \param mat_im imaginary-part block-cyclic dbcsr output
!> \param matrix_row row of rsmat to transform. Spin-free operators such
!>        as S and T must pass kpoint_spin_free: their first matrix
!>        dimension counts real-space blocks, not spins, so passing the
!>        instance spin of a spin loop reads a different block and no size
!>        check can catch it
! **************************************************************************************************
   SUBROUTINE kpoint_operator_get(ctx, ikp, ispin, rsmat, cmat, fm_re, fm_im, &
                                  mat_re, mat_im, matrix_row)

      TYPE(kpoint_operator_context_type), INTENT(INOUT)  :: ctx
      INTEGER, INTENT(IN)                                :: ikp, ispin
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rsmat
      TYPE(cp_cfm_type), OPTIONAL                        :: cmat
      TYPE(cp_fm_type), OPTIONAL                         :: fm_re, fm_im
      TYPE(dbcsr_type), INTENT(INOUT), OPTIONAL          :: mat_re, mat_im
      INTEGER, INTENT(IN)                                :: matrix_row

      CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_operator_get'

      INTEGER                                            :: handle, igroup, ik, kplocal
      INTEGER, DIMENSION(2)                              :: kp_range
      INTEGER, DIMENSION(:, :), POINTER                  :: kp_dist
      LOGICAL                                            :: my_kpgrp
      TYPE(copy_info_type)                               :: info_im, info_re
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      CPASSERT(COUNT([PRESENT(cmat), PRESENT(fm_re), PRESENT(mat_re)]) == 1)
      IF (PRESENT(fm_re) .OR. PRESENT(fm_im)) THEN
         CPASSERT(PRESENT(fm_re) .AND. PRESENT(fm_im))
      END IF
      IF (PRESENT(mat_re) .OR. PRESENT(mat_im)) THEN
         CPASSERT(PRESENT(mat_re) .AND. PRESENT(mat_im))
      END IF

      CPASSERT(ASSOCIATED(ctx%kpoints))
      CPASSERT(ASSOCIATED(ctx%kpoints%sab_nl))
      CPASSERT(ispin >= 1 .AND. ispin <= ctx%nspin)
      CPASSERT(matrix_row >= 1 .AND. matrix_row <= SIZE(rsmat, 1))
      CALL kpoint_operator_ensure(ctx, rsmat)
      kplocal = ctx%kpoints%kp_range(2) - ctx%kpoints%kp_range(1) + 1
      CPASSERT(ikp >= 1 .AND. ikp <= kplocal)
      IF (PRESENT(cmat)) THEN
         ! The merge writes cmat in place. A foreign distribution would fail
         ! inside the BLACS calls, far from this argument.
         CPASSERT(cp_fm_struct_equivalent(cmat%matrix_struct, ctx%fmlocal%matrix_struct))
      END IF

      kp_range = ctx%kpoints%kp_range
      kp_dist => ctx%kpoints%kp_dist
      para_env => ctx%kpoints%blacs_env_all%para_env
      DO igroup = 1, ctx%kpoints%nkp_groups
         ik = kp_dist(1, igroup) + ikp - 1
         my_kpgrp = (ik >= kp_range(1) .AND. ik <= kp_range(2))
         IF (PRESENT(cmat)) THEN
            ! cmat flavor: both parts pass through the context buffer; the
            ! real part overwrites cmat, the imaginary part adds to it
            CALL kpoint_operator_launch(rsmat, matrix_row, ik, ctx%kpoints%xkp(1:3, ik), &
                                        ctx%kpoints%cell_to_index, ctx%kpoints%sab_nl, &
                                        use_grid=.FALSE., rmatrix=ctx%rmatrix, cmatrix=ctx%cmatrix, &
                                        tmpmat=ctx%tmpmat, fm_re=ctx%fmwork(kpoint_part_re), &
                                        fm_im=ctx%fmwork(kpoint_part_im), target_re=ctx%fmlocal, &
                                        target_im=ctx%fmlocal, fmdummy=ctx%fmdummy, my_kpgrp=my_kpgrp, &
                                        para_env=para_env, info_re=info_re, info_im=info_im)
            IF (my_kpgrp) THEN
               CALL kpoint_operator_raw_finish_cfm(ctx%fmlocal, info_re, cmat, .FALSE.)
               CALL kpoint_operator_raw_finish_cfm(ctx%fmlocal, info_im, cmat, .TRUE.)
            END IF
         ELSE IF (PRESENT(fm_re)) THEN
            ! fm flavor: transfer each part straight into the caller matrices
            CALL kpoint_operator_launch(rsmat, matrix_row, ik, ctx%kpoints%xkp(1:3, ik), &
                                        ctx%kpoints%cell_to_index, ctx%kpoints%sab_nl, &
                                        use_grid=.FALSE., rmatrix=ctx%rmatrix, cmatrix=ctx%cmatrix, &
                                        tmpmat=ctx%tmpmat, fm_re=ctx%fmwork(kpoint_part_re), &
                                        fm_im=ctx%fmwork(kpoint_part_im), target_re=fm_re, &
                                        target_im=fm_im, fmdummy=ctx%fmdummy, my_kpgrp=my_kpgrp, &
                                        para_env=para_env, info_re=info_re, info_im=info_im)
            IF (my_kpgrp) THEN
               CALL cp_fm_finish_copy_general(fm_re, info_re)
               CALL cp_fm_finish_copy_general(fm_im, info_im)
            END IF
         ELSE
            ! dbcsr flavor: receive both parts through the context buffer,
            ! then emit the block-cyclic pair
            CALL kpoint_operator_launch(rsmat, matrix_row, ik, ctx%kpoints%xkp(1:3, ik), &
                                        ctx%kpoints%cell_to_index, ctx%kpoints%sab_nl, &
                                        use_grid=.FALSE., rmatrix=ctx%rmatrix, cmatrix=ctx%cmatrix, &
                                        tmpmat=ctx%tmpmat, fm_re=ctx%fmwork(kpoint_part_re), &
                                        fm_im=ctx%fmwork(kpoint_part_im), target_re=ctx%fmlocal, &
                                        target_im=ctx%fmlocal, fmdummy=ctx%fmdummy, my_kpgrp=my_kpgrp, &
                                        para_env=para_env, info_re=info_re, info_im=info_im)
            IF (my_kpgrp) THEN
               CALL cp_fm_finish_copy_general(ctx%fmlocal, info_re)
               CALL copy_fm_to_dbcsr_bc(ctx%fmlocal, mat_re)
               CALL cp_fm_finish_copy_general(ctx%fmlocal, info_im)
               CALL copy_fm_to_dbcsr_bc(ctx%fmlocal, mat_im)
            END IF
         END IF
         CALL cp_fm_cleanup_copy_general(info_re)
         CALL cp_fm_cleanup_copy_general(info_im)
      END DO

      CALL timestop(handle)

   END SUBROUTINE kpoint_operator_get

END MODULE qs_kpoint_operators
