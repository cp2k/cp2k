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
!>        One process can start many operators before it finishes any of them.
!>        The STANDARD k-point driver uses this order.
!>        The routines start, finish and get are collective over the
!>        environment-wide communicator. Every rank must call the same
!>        sequence of service calls.
!>        A context references its host kpoint environment without owning it,
!>        so it must not outlive that environment.
!>        An explicitly prepared, immutable operator batch can instead be built
!>        from bounded context-owned real-space snapshots inside each group.
!>        The public round and output-ownership contracts stay the same.
!>        The calls are not thread safe: the caller must serialize them, because their collective
!>        communication tolerates no concurrent service call. A separate
!>        context per thread only keeps the scratch from aliasing.
! **************************************************************************************************
MODULE qs_kpoint_operators
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_scale_and_add_fm
   USE cp_cfm_types,                    ONLY: cp_cfm_to_fm,&
                                              cp_cfm_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_create, dbcsr_deallocate_matrix, dbcsr_get_info, dbcsr_get_readonly_block_p, &
        dbcsr_p_type, dbcsr_set, dbcsr_type, dbcsr_type_antisymmetric, dbcsr_type_symmetric
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr_bc,&
                                              dbcsr_fm_copy_plan_type
   USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                              fm_pool_create_fm,&
                                              fm_pool_give_back_fm
   USE cp_fm_struct,                    ONLY: cp_fm_struct_equivalent,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_retain,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: copy_info_type,&
                                              cp_fm_cleanup_copy_general,&
                                              cp_fm_create,&
                                              cp_fm_finish_copy_general,&
                                              cp_fm_get_info,&
                                              cp_fm_release,&
                                              cp_fm_start_copy_general,&
                                              cp_fm_type
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE kpoint_methods,                  ONLY: kp_transform_plan_create,&
                                              kp_transform_plan_type,&
                                              rskp_grid_type,&
                                              rskp_transform,&
                                              rskp_transform_grid_extract
   USE kpoint_types,                    ONLY: kpoint_type
   USE mathconstants,                   ONLY: gaussi,&
                                              twopi,&
                                              z_one,&
                                              z_zero
   USE message_passing,                 ONLY: mp_para_env_type
   USE qs_matrix_pools,                 ONLY: mpools_get
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
   PUBLIC :: kpoint_operator_get, kpoint_operator_prepare

   INTEGER, PARAMETER, PUBLIC :: kpoint_slot_ks = 1, &
                                 kpoint_slot_s = 2, &
                                 kpoint_slot_t = 3

   !> Row of rsmat that carries the spin-free operator image. S and T carry
   !> no spin copies: drivers launch one instance per spin of their loop and
   !> read this same row in every instance.
   INTEGER, PARAMETER, PUBLIC :: kpoint_spin_free = 1

   !> Number of defined slots; grows only inside this module.
   INTEGER, PARAMETER, PRIVATE :: kpoint_num_slots = kpoint_slot_t - kpoint_slot_ks + 1

   !> Part indices of one transfer: its real and its imaginary half.
   INTEGER, PARAMETER, PRIVATE :: kpoint_part_re = 1, &
                                  kpoint_part_im = 2

   !> Instance states of one (k point, spin, slot) transfer.
   !> State transitions happen only inside the service routines.
   INTEGER, PARAMETER, PRIVATE :: kpoint_slot_idle = 0, &
                                  kpoint_slot_in_flight = 1, &
                                  kpoint_slot_local = 2

   ! Structural route shared by all captured rows. Whole blocks go to one
   ! representative group, then only equivalent tiles are replicated.
   TYPE kpoint_operator_rs_route_type
      TYPE(kp_transform_plan_type) :: traversal = kp_transform_plan_type()
      TYPE(cp_fm_struct_type), POINTER :: matrix_struct => NULL()
      INTEGER :: nimage = 0, nvalue = 0
      INTEGER(KIND=int_8) :: storage = 0
      INTEGER, ALLOCATABLE :: row_size(:), col_size(:), ENTRY(:, :), first(:), row_local(:)
      INTEGER, ALLOCATABLE :: send_entry(:), send_first(:), scount(:, :), rcount(:, :), &
                              sdisp(:, :), rdisp(:, :)
      REAL(KIND=dp), ALLOCATABLE :: send_values(:)
   END TYPE kpoint_operator_rs_route_type

   ! Immutable values, not per-k matrices. Metadata belong to the shared route;
   ! source identifies starts belonging to the explicitly prepared batch only.
   TYPE kpoint_operator_rs_type
      INTEGER :: matrix_row = 0
      TYPE(dbcsr_p_type), POINTER :: source(:, :) => NULL()
      REAL(KIND=dp), ALLOCATABLE :: values(:)
      COMPLEX(KIND=dp), ALLOCATABLE :: phase(:)
      LOGICAL, ALLOCATABLE :: nonzero(:)
   END TYPE kpoint_operator_rs_type

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
      ! operator such as S counts derivatives, and the images live in the
      ! second, so callers address operator rows through matrix_row instead.
      INTEGER                                            :: nspin = 0
      ! Borrowed host/distribution references; the context cannot span an energy evaluation.
      TYPE(kpoint_type), POINTER                          :: kpoints => NULL()
      TYPE(cp_fm_struct_type), POINTER                    :: ao_ao_fmstruct => NULL()
      TYPE(cp_fm_pool_p_type), DIMENSION(:), POINTER      :: group_fm_pools => NULL()
      ! Identity of the neighbor lists used by scratch and prepared snapshots.
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                          :: work_sab => NULL()
      ! Owned scratch, created on first use and released with the context.
      TYPE(dbcsr_type), POINTER                           :: rmatrix => NULL()
      TYPE(dbcsr_type), POINTER                           :: cmatrix => NULL()
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)         :: fmwork
      TYPE(cp_fm_type), ALLOCATABLE                       :: fmlocal, fmimag
      ! Plans share the lifetime and fixed layouts of the two densify workspaces.
      TYPE(dbcsr_fm_copy_plan_type), ALLOCATABLE          :: copy_plan(:)
      ! Optional immutable operator snapshots. No per-k matrices or outstanding
      ! MPI requests are needed for starts served by this backend.
      TYPE(kpoint_operator_rs_type), ALLOCATABLE          :: rs(:, :)
      TYPE(kpoint_operator_rs_route_type), ALLOCATABLE    :: rs_route
      ! Bookkeeping. info is keyed (global k point, spin, slot, part),
      ! part kpoint_part_re holding the real half and kpoint_part_im holding
      ! the imaginary half. Status is keyed (global round, spin, slot)
      ! because start and finish treat all groups of one instance
      ! atomically. Allocated on first start.
      TYPE(copy_info_type), ALLOCATABLE, DIMENSION(:, :, :, :) :: info
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)            :: status
   END TYPE kpoint_operator_context_type

CONTAINS

! **************************************************************************************************
!> \brief Transform and densify one operator using the context-owned workspaces and routing plans.
!> \param ctx initialized transfer context
!> \param rsmat real-space operator images
!> \param matrix_row operator row, independent of the instance spin
!> \param ik global k-point index
!> \param grid caller-owned reciprocal-grid cache
!> \param use_grid extract from grid instead of a direct phase sum
! **************************************************************************************************
   SUBROUTINE kpoint_operator_densify(ctx, rsmat, matrix_row, ik, grid, use_grid)

      TYPE(kpoint_operator_context_type), INTENT(INOUT)  :: ctx
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rsmat
      INTEGER, INTENT(IN)                                :: matrix_row, ik
      TYPE(rskp_grid_type), INTENT(IN), OPTIONAL         :: grid
      LOGICAL, INTENT(IN)                                :: use_grid

      CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_operator_densify'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (use_grid) THEN
         CPASSERT(PRESENT(grid))
         CALL rskp_transform_grid_extract(grid, ik, ctx%rmatrix, ctx%cmatrix)
      ELSE
         CALL dbcsr_set(ctx%rmatrix, 0.0_dp)
         CALL dbcsr_set(ctx%cmatrix, 0.0_dp)
         CALL rskp_transform(rmatrix=ctx%rmatrix, cmatrix=ctx%cmatrix, rsmat=rsmat, ispin=matrix_row, &
                             xkp=ctx%kpoints%xkp(:, ik), cell_to_index=ctx%kpoints%cell_to_index, &
                             sab_nl=ctx%kpoints%sab_nl)
      END IF
      CALL copy_dbcsr_to_fm(ctx%rmatrix, ctx%fmwork(kpoint_part_re), ctx%copy_plan(kpoint_part_re))
      CALL copy_dbcsr_to_fm(ctx%cmatrix, ctx%fmwork(kpoint_part_im), ctx%copy_plan(kpoint_part_im))

      CALL timestop(handle)

   END SUBROUTINE kpoint_operator_densify

! **************************************************************************************************
!> \brief Launch a densified operator's two parts to its owning group.
!>        Absent explicit targets, receive both parts through fmlocal. Each copy
!>        owns its receive buffer until finish, so the targets can alias.
!>        Non-owning groups participate with an empty dummy target.
!> \param kpoints host environment, including group ownership and transfer communicator
!> \param ik global k-point index
!> \param fmwork real/imaginary source matrices on the global distribution
!> \param fmlocal default receiver on the k-point group distribution
!> \param info_re real-part transfer state, consumed by finish and cleanup
!> \param info_im imaginary-part transfer state, consumed by finish and cleanup
!> \param target_re optional caller-owned real-part target
!> \param target_im optional caller-owned imaginary-part target; must accompany target_re
! **************************************************************************************************
   SUBROUTINE kpoint_operator_launch(kpoints, ik, fmwork, fmlocal, info_re, info_im, target_re, target_im)

      TYPE(kpoint_type), POINTER                         :: kpoints
      INTEGER, INTENT(IN)                                :: ik
      TYPE(cp_fm_type), DIMENSION(:), INTENT(IN)         :: fmwork
      TYPE(cp_fm_type), INTENT(IN)                       :: fmlocal
      TYPE(copy_info_type)                               :: info_re, info_im
      TYPE(cp_fm_type), INTENT(IN), OPTIONAL             :: target_re, target_im

      CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_operator_launch'

      INTEGER                                            :: handle
      TYPE(cp_fm_type)                                   :: fmdummy
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)
      CPASSERT(PRESENT(target_re) .EQV. PRESENT(target_im))

      para_env => kpoints%blacs_env_all%para_env
      IF (ik >= kpoints%kp_range(1) .AND. ik <= kpoints%kp_range(2)) THEN
         IF (PRESENT(target_re)) THEN
            CALL cp_fm_start_copy_general(fmwork(kpoint_part_re), target_re, para_env, info_re)
            CALL cp_fm_start_copy_general(fmwork(kpoint_part_im), target_im, para_env, info_im)
         ELSE
            CALL cp_fm_start_copy_general(fmwork(kpoint_part_re), fmlocal, para_env, info_re)
            CALL cp_fm_start_copy_general(fmwork(kpoint_part_im), fmlocal, para_env, info_im)
         END IF
      ELSE
         CALL cp_fm_start_copy_general(fmwork(kpoint_part_re), fmdummy, para_env, info_re)
         CALL cp_fm_start_copy_general(fmwork(kpoint_part_im), fmdummy, para_env, info_im)
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
      mo_im%eigenvalues(:) = eigenvalues
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
      NULLIFY (ctx%rmatrix, ctx%cmatrix, ctx%work_sab)
      IF (ALLOCATED(ctx%copy_plan)) DEALLOCATE (ctx%copy_plan)
      IF (ALLOCATED(ctx%rs)) DEALLOCATE (ctx%rs)
      IF (ALLOCATED(ctx%rs_route)) THEN
         CALL cp_fm_struct_release(ctx%rs_route%matrix_struct)
         DEALLOCATE (ctx%rs_route)
      END IF
      IF (ALLOCATED(ctx%fmwork)) THEN
         DO iwork = 1, SIZE(ctx%fmwork)
            CALL cp_fm_release(ctx%fmwork(iwork))
         END DO
         DEALLOCATE (ctx%fmwork)
      END IF
      IF (ALLOCATED(ctx%fmimag)) THEN
         CALL cp_fm_release(ctx%fmimag)
         DEALLOCATE (ctx%fmimag)
      END IF
      IF (ALLOCATED(ctx%fmlocal)) THEN
         CALL fm_pool_give_back_fm(ctx%group_fm_pools(1)%pool, ctx%fmlocal)
         DEALLOCATE (ctx%fmlocal)
      END IF
      IF (ALLOCATED(ctx%info)) DEALLOCATE (ctx%info)
      IF (ALLOCATED(ctx%status)) DEALLOCATE (ctx%status)
      NULLIFY (ctx%kpoints, ctx%ao_ao_fmstruct, ctx%group_fm_pools)

      CALL timestop(handle)

   END SUBROUTINE kpoint_operator_context_release

! **************************************************************************************************
!> \brief Create the context scratch that the first service call needs: the
!>        DBCSR work pair with blocks from the current neighbor lists, two
!>        densify buffers on the environment-wide struct, and the group
!>        transfer target from the group pools.
!>        The per-instance slot bookkeeping is not scratch: start creates it
!>        on its first call, because the synchronous get call never touches it.
!>        The buffers assume the context lifetime of one driver call: the
!>        host neighbor lists do not change within it, and every buffer use
!>        asserts this identity. Buffers from an earlier call stay in place.
!>        Called from the service routines only.
!> \param ctx context that receives the buffers
!> \param rsmat real-space image matrices of one operator. The call uses the
!>        first matrix as the DBCSR template of the work pair. It reads no
!>        matrix data.
! **************************************************************************************************
   SUBROUTINE kpoint_operator_ensure(ctx, rsmat)

      TYPE(kpoint_operator_context_type), INTENT(INOUT)  :: ctx
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rsmat

      CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_operator_ensure'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      ! The pair must still match the host neighbor lists: they are
      ! reallocated at every energy evaluation, while the context spans one
      ! driver call within it.
      IF (ASSOCIATED(ctx%work_sab)) THEN
         CPASSERT(ASSOCIATED(ctx%work_sab, ctx%kpoints%sab_nl))
      END IF
      IF (.NOT. ASSOCIATED(ctx%rmatrix)) THEN
         ALLOCATE (ctx%rmatrix, ctx%cmatrix)
         CALL dbcsr_create(ctx%rmatrix, template=rsmat(1, 1)%matrix, matrix_type=dbcsr_type_symmetric)
         CALL dbcsr_create(ctx%cmatrix, template=rsmat(1, 1)%matrix, matrix_type=dbcsr_type_antisymmetric)
         CALL cp_dbcsr_alloc_block_from_nbl(ctx%rmatrix, ctx%kpoints%sab_nl)
         CALL cp_dbcsr_alloc_block_from_nbl(ctx%cmatrix, ctx%kpoints%sab_nl)
         ctx%work_sab => ctx%kpoints%sab_nl
      END IF

      IF (.NOT. ALLOCATED(ctx%fmwork)) THEN
         ALLOCATE (ctx%fmwork(kpoint_part_im), ctx%copy_plan(kpoint_part_im))
         CALL cp_fm_create(ctx%fmwork(kpoint_part_re), ctx%ao_ao_fmstruct, &
                           name="KPOINT-OPERATOR-WORK-RE")
         CALL cp_fm_create(ctx%fmwork(kpoint_part_im), ctx%ao_ao_fmstruct, &
                           name="KPOINT-OPERATOR-WORK-IM")
      END IF

      CALL kpoint_operator_ensure_group(ctx)

      CALL timestop(handle)

   END SUBROUTINE kpoint_operator_ensure

! **************************************************************************************************
!> \brief Acquire only the group-local output workspace; both backends share it.
!> \param ctx operator context
! **************************************************************************************************
   SUBROUTINE kpoint_operator_ensure_group(ctx)

      TYPE(kpoint_operator_context_type), INTENT(INOUT)  :: ctx

      IF (.NOT. ALLOCATED(ctx%fmlocal)) THEN
         CALL mpools_get(ctx%kpoints%mpools, ao_ao_fm_pools=ctx%group_fm_pools)
         ALLOCATE (ctx%fmlocal)
         CALL fm_pool_create_fm(ctx%group_fm_pools(1)%pool, ctx%fmlocal)
      END IF

   END SUBROUTINE kpoint_operator_ensure_group

! **************************************************************************************************
!> \brief Prepare an immutable operator slot for group-local assembly.
!>        Refresh after the source values change, with every instance idle.
!>        Other prepared slots and structural routing remain valid. Without
!>        matrix_row, capture one row per instance spin; with matrix_row, store
!>        one row shared by all instance spins. Failure clears this slot on
!>        all ranks and retains the existing transfer fallback.
!> \param ctx context spanning the current initialized k-point environment
!> \param rsmat current real-space operator images
!> \param slot slot to replace (KS, S or T)
!> \param prepared optional collective success flag
!> \param matrix_row optional explicit spin-free operator row
!> \param TARGET optional group-local FM layout for prepared get outputs;
!>        otherwise use the host AO/AO pool layout. One layout per context
! **************************************************************************************************
   SUBROUTINE kpoint_operator_prepare(ctx, rsmat, slot, prepared, matrix_row, TARGET)

      TYPE(kpoint_operator_context_type), INTENT(INOUT)  :: ctx
      TYPE(dbcsr_p_type), POINTER                        :: rsmat(:, :)
      INTEGER, INTENT(IN)                                :: slot
      LOGICAL, INTENT(OUT), OPTIONAL                     :: prepared
      INTEGER, INTENT(IN), OPTIONAL                      :: matrix_row
      TYPE(cp_fm_type), INTENT(IN), OPTIONAL             :: target

      CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_operator_prepare'
      INTEGER(KIND=int_8), PARAMETER :: max_bytes = 64_int_8*1024_int_8**2

      INTEGER                                            :: compatible, handle, ispin, jslot, nrow, &
                                                            row
      INTEGER(KIND=int_8)                                :: retained, storage
      LOGICAL                                            :: ready
      TYPE(cp_fm_struct_type), POINTER                   :: fms
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)
      IF (PRESENT(prepared)) prepared = .FALSE.
      CPASSERT(ASSOCIATED(ctx%kpoints))
      CPASSERT(slot >= kpoint_slot_ks .AND. slot <= kpoint_slot_t)
      CPASSERT(ASSOCIATED(ctx%kpoints%sab_nl))
      IF (ASSOCIATED(ctx%work_sab)) THEN
         CPASSERT(ASSOCIATED(ctx%work_sab, ctx%kpoints%sab_nl))
      END IF
      IF (ALLOCATED(ctx%status)) THEN
         CPASSERT(ALL(ctx%status == kpoint_slot_idle))
      END IF
      IF (.NOT. ALLOCATED(ctx%rs)) ALLOCATE (ctx%rs(ctx%nspin, kpoint_num_slots))
      DO ispin = 1, ctx%nspin
         CALL kpoint_operator_rs_clear(ctx%rs(ispin, slot))
      END DO
      IF (ctx%kpoints%nkp_groups <= 1) THEN
         CALL timestop(handle)
         RETURN
      END IF
      nrow = ctx%nspin
      IF (PRESENT(matrix_row)) THEN
         CPASSERT(matrix_row >= 1 .AND. matrix_row <= SIZE(rsmat, 1))
         nrow = 1
      ELSE
         CPASSERT(SIZE(rsmat, 1) == ctx%nspin)
      END IF
      para_env => ctx%kpoints%blacs_env_all%para_env
      CALL kpoint_operator_ensure_group(ctx)
      fms => ctx%fmlocal%matrix_struct
      IF (PRESENT(TARGET)) fms => target%matrix_struct
      CPASSERT(ASSOCIATED(fms))
      ! Routing uses the host group-rank numbering. Different BLACS layouts
      ! on that communicator are supported, foreign communicators are not.
      compatible = MERGE(1, 0, ASSOCIATED(fms%para_env, ctx%kpoints%para_env_kp))
      CALL para_env%min(compatible)
      IF (compatible == 0) THEN
         CALL timestop(handle)
         RETURN
      END IF
      IF (.NOT. ALLOCATED(ctx%rs_route)) THEN
         ALLOCATE (ctx%rs_route)
         CALL kpoint_operator_rs_plan_create(ctx%kpoints, rsmat, fms, max_bytes, ctx%rs_route, ready)
         IF (.NOT. ready) THEN
            DEALLOCATE (ctx%rs_route)
            CALL timestop(handle)
            RETURN
         END IF
         ctx%work_sab => ctx%kpoints%sab_nl
      END IF
      compatible = MERGE(1, 0, cp_fm_struct_equivalent(fms, ctx%rs_route%matrix_struct))
      DO ispin = 1, nrow
         row = ispin
         IF (PRESENT(matrix_row)) row = matrix_row
         IF (.NOT. kpoint_operator_rs_compatible(ctx%rs_route, rsmat, row)) compatible = 0
      END DO
      CALL para_env%min(compatible)
      retained = 0_int_8
      DO jslot = kpoint_slot_ks, kpoint_slot_t
         DO ispin = 1, ctx%nspin
            IF (.NOT. ALLOCATED(ctx%rs(ispin, jslot)%values)) CYCLE
            retained = retained + 8_int_8*SIZE(ctx%rs(ispin, jslot)%values, KIND=int_8) + &
                       20_int_8*SIZE(ctx%rs(ispin, jslot)%phase, KIND=int_8)
         END DO
      END DO
      storage = ctx%rs_route%storage + retained + INT(nrow, int_8)* &
                (8_int_8*ctx%rs_route%nvalue + 20_int_8*SIZE(ctx%rs_route%entry, 2))
      CALL para_env%max(storage)
      ready = compatible == 1 .AND. storage <= max_bytes
      IF (ready) THEN
         DO ispin = 1, nrow
            row = ispin
            IF (PRESENT(matrix_row)) row = matrix_row
            CALL kpoint_operator_rs_capture(rsmat, row, ctx%rs_route, ctx%kpoints, ctx%rs(ispin, slot))
         END DO
      END IF
      IF (PRESENT(prepared)) prepared = ready
      CALL timestop(handle)

   END SUBROUTINE kpoint_operator_prepare

! **************************************************************************************************
!> \brief Release one captured row without invalidating other slots or routing.
!> \param snapshot row to clear
! **************************************************************************************************
   SUBROUTINE kpoint_operator_rs_clear(snapshot)

      TYPE(kpoint_operator_rs_type), INTENT(INOUT)       :: snapshot

      IF (ALLOCATED(snapshot%values)) DEALLOCATE (snapshot%values, snapshot%phase, snapshot%nonzero)
      NULLIFY (snapshot%source)
      snapshot%matrix_row = 0

   END SUBROUTINE kpoint_operator_rs_clear

! **************************************************************************************************
!> \brief Match an explicitly prepared row, sharing spin-free rows across spins.
!> \param ctx operator context
!> \param source live array identity identifying the prepared operation
!> \param ispin instance spin
!> \param slot prepared slot
!> \param matrix_row operator row
!> \return snapshot spin index, or zero for the transfer backend
! **************************************************************************************************
   FUNCTION kpoint_operator_rs_match(ctx, source, ispin, slot, matrix_row) RESULT(irs)

      TYPE(kpoint_operator_context_type), INTENT(IN)     :: ctx
      TYPE(dbcsr_p_type), POINTER                        :: source(:, :)
      INTEGER, INTENT(IN)                                :: ispin, slot, matrix_row
      INTEGER                                            :: irs

      INTEGER                                            :: i

      irs = 0
      IF (.NOT. ALLOCATED(ctx%rs)) RETURN
      i = ispin
      IF (.NOT. ASSOCIATED(ctx%rs(i, slot)%source)) i = 1
      IF (ASSOCIATED(ctx%rs(i, slot)%source, source) .AND. ctx%rs(i, slot)%matrix_row == matrix_row) irs = i

   END FUNCTION kpoint_operator_rs_match

! **************************************************************************************************
!> \brief Cache whole-block routing to the FM tiles of one representative group.
!>        Equivalent ranks in the other groups receive only their tile's blocks.
!>        Boundary blocks may overlap several tiles. Missing blocks are packed
!>        as zero, permitting different operator sparsities on a fixed AO layout.
!> \param kpoints initialized parallel environment and neighbor lists
!> \param source operator layout template
!> \param matrix_struct target group layout
!> \param max_bytes per-rank budget for routing and capture buffers
!> \param route context-owned structural route
!> \param ready false collectively if unsupported or over budget
! **************************************************************************************************
   SUBROUTINE kpoint_operator_rs_plan_create(kpoints, source, matrix_struct, max_bytes, route, ready)

      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(dbcsr_p_type), POINTER                        :: source(:, :)
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct
      INTEGER(KIND=int_8), INTENT(IN)                    :: max_bytes
      TYPE(kpoint_operator_rs_route_type), INTENT(OUT)   :: route
      LOGICAL, INTENT(OUT)                               :: ready

      CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_operator_rs_plan_create'

      INTEGER                                            :: compatible, dest, handle, i, iblock, &
                                                            nblock, ncol, nrow, p, pos, rank
      INTEGER(KIND=int_8)                                :: nr(2), ns(2), storage
      INTEGER(KIND=int_8), ALLOCATABLE                   :: recv_size(:, :), send_size(:, :)
      INTEGER, ALLOCATABLE                               :: identities(:, :), owner(:), &
                                                            position(:, :), send_metadata(:, :)
      INTEGER, POINTER                                   :: col_size(:), row_size(:)
      LOGICAL, ALLOCATABLE                               :: needed(:)
      TYPE(mp_para_env_type), POINTER                    :: inter, world

      CALL timeset(routineN, handle)
      ready = .FALSE.
      world => kpoints%blacs_env_all%para_env
      inter => kpoints%para_env_inter_kp
      CALL dbcsr_get_info(source(1, 1)%matrix, row_blk_size=row_size, col_blk_size=col_size)
      CALL kp_transform_plan_create(route%traversal, kpoints%sab_nl, kpoints%cell_to_index, &
                                    SIZE(source, 2), block_template=source(1, 1)%matrix)
      compatible = MERGE(1, 0, route%traversal%symmetric)
      IF (SUM(row_size) /= matrix_struct%nrow_global .OR. &
          SUM(col_size) /= matrix_struct%ncol_global) compatible = 0
      IF (SIZE(row_size) == SIZE(col_size)) THEN
         IF (ANY(row_size /= col_size)) compatible = 0
      ELSE
         compatible = 0
      END IF
      CALL world%min(compatible)
      IF (compatible == 0) THEN
         CALL timestop(handle)
         RETURN
      END IF
      route%nimage = SIZE(source, 2)
      ALLOCATE (route%row_size(SIZE(row_size)), route%col_size(SIZE(col_size)))
      route%row_size(:) = row_size
      route%col_size(:) = col_size
      ALLOCATE (identities(2, 0:world%num_pe - 1), owner(0:kpoints%para_env_kp%num_pe - 1), &
                needed(0:kpoints%para_env_kp%num_pe - 1))
      CALL world%allgather([kpoints%para_env_kp%mepos, inter%mepos], identities)
      owner(:) = -1
      DO rank = 0, world%num_pe - 1
         IF (identities(2, rank) == 0) owner(identities(1, rank)) = rank
      END DO
      CPASSERT(ALL(owner >= 0))
      ALLOCATE (send_size(2, 0:world%num_pe - 1), recv_size(2, 0:world%num_pe - 1))
      send_size(:, :) = 0_int_8
      DO i = 1, route%traversal%nentry
         nrow = row_size(route%traversal%row(i))
         ncol = col_size(route%traversal%col(i))
         CALL kpoint_operator_rs_owners(route, matrix_struct, i, needed)
         DO p = 0, SIZE(owner) - 1
            IF (.NOT. needed(p)) CYCLE
            dest = owner(p)
            send_size(1, dest) = send_size(1, dest) + 1_int_8
            send_size(2, dest) = send_size(2, dest) + INT(nrow, int_8)*ncol
         END DO
      END DO
      CALL world%alltoall(send_size, recv_size, 2)
      ! Only inter-rank zero receives from the world. The other groups hold
      ! equivalent FM tiles, as required by kpoint_env_initialize.
      ns = SUM(send_size, DIM=2)
      nr = SUM(recv_size, DIM=2)
      CALL inter%bcast(nr, 0)
      ! Conservative peak: traversal, retained route, temporary metadata and
      ! one snapshot. The full batch budget is checked again by prepare.
      storage = 40_int_8*route%traversal%nentry + &
                4_int_8*(SIZE(row_size) + SIZE(col_size) + matrix_struct%nrow_global) + &
                48_int_8*ns(1) + 8_int_8*ns(2) + 36_int_8*nr(1) + &
                128_int_8*world%num_pe + 8_int_8*SIZE(owner) + 8_int_8
      route%storage = storage
      storage = storage + 8_int_8*nr(2) + 20_int_8*nr(1)
      CALL world%max(storage)
      IF (storage > max_bytes) THEN
         CALL timestop(handle)
         RETURN
      END IF
      ! The byte bound also keeps all counts/displacements within default INTEGER.
      nblock = INT(nr(1))
      route%nvalue = INT(nr(2))
      ALLOCATE (route%scount(2, 0:world%num_pe - 1), route%rcount(2, 0:world%num_pe - 1), &
                route%sdisp(2, 0:world%num_pe - 1), route%rdisp(2, 0:world%num_pe - 1))
      route%scount(:, :) = INT(send_size)
      route%rcount(:, :) = INT(recv_size)
      IF (SIZE(owner) == 1) CALL inter%bcast(route%rcount, 0)
      route%sdisp(:, 0) = 0
      route%rdisp(:, 0) = 0
      DO rank = 1, world%num_pe - 1
         route%sdisp(:, rank) = route%sdisp(:, rank - 1) + route%scount(:, rank - 1)
         route%rdisp(:, rank) = route%rdisp(:, rank - 1) + route%rcount(:, rank - 1)
      END DO
      ALLOCATE (route%send_entry(INT(ns(1))), route%send_first(INT(ns(1)) + 1), &
                route%send_values(INT(ns(2))), send_metadata(8, INT(ns(1))), &
                position(2, 0:world%num_pe - 1), route%entry(8, nblock), &
                route%first(nblock + 1), route%row_local(matrix_struct%nrow_global))
      position(:, :) = route%sdisp
      DO i = 1, route%traversal%nentry
         nrow = row_size(route%traversal%row(i))
         ncol = col_size(route%traversal%col(i))
         CALL kpoint_operator_rs_owners(route, matrix_struct, i, needed)
         DO p = 0, SIZE(owner) - 1
            IF (.NOT. needed(p)) CYCLE
            dest = owner(p)
            iblock = position(1, dest) + 1
            pos = position(2, dest) + 1
            route%send_entry(iblock) = i
            route%send_first(iblock) = pos
            send_metadata(1:4, iblock) = [route%traversal%row_offset(i), route%traversal%col_offset(i), nrow, ncol]
            send_metadata(5:7, iblock) = NINT(route%traversal%symmetry_sign(i))*route%traversal%cell(:, i)
            send_metadata(8, iblock) = MERGE(1, 0, route%traversal%row(i) /= route%traversal%col(i))
            position(:, dest) = position(:, dest) + [1, nrow*ncol]
         END DO
      END DO
      route%send_first(INT(ns(1)) + 1) = INT(ns(2)) + 1
      IF (SIZE(owner) == 1) THEN
         CALL world%allgatherv(send_metadata, route%entry, 8*route%rcount(1, :), 8*route%rdisp(1, :))
      ELSE
         CALL world%alltoall(send_metadata, 8*route%scount(1, :), 8*route%sdisp(1, :), &
                             route%entry, 8*route%rcount(1, :), 8*route%rdisp(1, :))
         CALL inter%bcast(route%entry, 0)
      END IF
      route%first(1) = 1
      DO i = 1, nblock
         route%first(i + 1) = route%first(i) + route%entry(3, i)*route%entry(4, i)
      END DO
      CPASSERT(route%first(nblock + 1) == route%nvalue + 1)
      route%row_local(:) = 0
      DO i = 1, SIZE(matrix_struct%row_indices)
         route%row_local(matrix_struct%row_indices(i)) = i
      END DO
      route%matrix_struct => matrix_struct
      CALL cp_fm_struct_retain(route%matrix_struct)
      ready = .TRUE.
      CALL timestop(handle)

   END SUBROUTINE kpoint_operator_rs_plan_create

! **************************************************************************************************
!> \brief Find group-local FM owners intersecting a stored block or its mirror.
!> \param route source AO block layout and neighbor traversal
!> \param fms target FM layout within each group
!> \param i traversal entry
!> \param needed recipient group-local ranks, without duplicate mirror sends
! **************************************************************************************************
   SUBROUTINE kpoint_operator_rs_owners(route, fms, i, needed)

      TYPE(kpoint_operator_rs_route_type), INTENT(IN)    :: route
      TYPE(cp_fm_struct_type), INTENT(IN)                :: fms
      INTEGER, INTENT(IN)                                :: i
      LOGICAL, INTENT(OUT)                               :: needed(0:)

      INTEGER                                            :: c0, j, nc, nr, pc, pr, r0
      LOGICAL :: rows(2, 0:fms%context%num_pe(1) - 1), cols(2, 0:fms%context%num_pe(2) - 1)

      r0 = route%traversal%row_offset(i)
      c0 = route%traversal%col_offset(i)
      nr = route%row_size(route%traversal%row(i))
      nc = route%col_size(route%traversal%col(i))
      rows(:, :) = .FALSE.
      cols(:, :) = .FALSE.
      DO j = r0, r0 + nr - 1
         rows(1, fms%g2p_row(j)) = .TRUE.
         cols(2, fms%g2p_col(j)) = .TRUE.
      END DO
      DO j = c0, c0 + nc - 1
         rows(2, fms%g2p_row(j)) = .TRUE.
         cols(1, fms%g2p_col(j)) = .TRUE.
      END DO
      IF (route%traversal%row(i) == route%traversal%col(i)) rows(2, :) = .FALSE.
      needed(:) = .FALSE.
      DO pc = 0, fms%context%num_pe(2) - 1
         DO pr = 0, fms%context%num_pe(1) - 1
            needed(fms%context%blacs2mpi(pr, pc)) = &
               (rows(1, pr) .AND. cols(1, pc)) .OR. (rows(2, pr) .AND. cols(2, pc))
         END DO
      END DO

   END SUBROUTINE kpoint_operator_rs_owners

! **************************************************************************************************
!> \brief Check the source AO layout, not its values or sparse block presence.
!> \param route fixed context route
!> \param source current operator images
!> \param matrix_row row to capture
!> \return whether the row has the route's image count and AO layout
! **************************************************************************************************
   FUNCTION kpoint_operator_rs_compatible(route, source, matrix_row) RESULT(compatible)

      TYPE(kpoint_operator_rs_route_type), INTENT(IN)    :: route
      TYPE(dbcsr_p_type), POINTER                        :: source(:, :)
      INTEGER, INTENT(IN)                                :: matrix_row
      LOGICAL                                            :: compatible

      INTEGER                                            :: image
      INTEGER, POINTER                                   :: col_size(:), row_size(:)

      compatible = .FALSE.
      IF (SIZE(source, 2) /= route%nimage) RETURN
      IF (matrix_row < 1 .OR. matrix_row > SIZE(source, 1)) RETURN
      DO image = 1, SIZE(source, 2)
         CALL dbcsr_get_info(source(matrix_row, image)%matrix, row_blk_size=row_size, col_blk_size=col_size)
         IF (SIZE(row_size) /= SIZE(route%row_size) .OR. SIZE(col_size) /= SIZE(route%col_size)) RETURN
         IF (ANY(row_size /= route%row_size) .OR. ANY(col_size /= route%col_size)) RETURN
      END DO
      compatible = .TRUE.

   END FUNCTION kpoint_operator_rs_compatible

! **************************************************************************************************
!> \brief Capture current values using cached routing, with no per-k communication.
!> \param source current real-space operator images
!> \param matrix_row row to snapshot
!> \param route context-owned metadata, counts and reusable packing buffer
!> \param kpoints initialized host parallel environment
!> \param snapshot immutable tile-local values
! **************************************************************************************************
   SUBROUTINE kpoint_operator_rs_capture(source, matrix_row, route, kpoints, snapshot)

      TYPE(dbcsr_p_type), POINTER                        :: source(:, :)
      INTEGER, INTENT(IN)                                :: matrix_row
      TYPE(kpoint_operator_rs_route_type), INTENT(INOUT) :: route
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(kpoint_operator_rs_type), INTENT(OUT)         :: snapshot

      CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_operator_rs_capture'

      INTEGER                                            :: col, handle, i, iblock, ncol, nrow, pos, &
                                                            row
      LOGICAL                                            :: found
      REAL(KIND=dp), POINTER                             :: block(:, :)
      TYPE(mp_para_env_type), POINTER                    :: inter, world

      CALL timeset(routineN, handle)
      world => kpoints%blacs_env_all%para_env
      inter => kpoints%para_env_inter_kp
      ALLOCATE (snapshot%values(route%nvalue), snapshot%phase(SIZE(route%entry, 2)), &
                snapshot%nonzero(SIZE(route%entry, 2)))
      route%send_values(:) = 0.0_dp
!$OMP PARALLEL DO DEFAULT(NONE) SHARED(source, matrix_row, route) &
!$OMP PRIVATE(iblock, i, block, found, nrow, ncol, pos, row, col) SCHEDULE(STATIC)
      DO iblock = 1, SIZE(route%send_entry)
         i = route%send_entry(iblock)
         CALL dbcsr_get_readonly_block_p(source(matrix_row, route%traversal%image(i))%matrix, &
                                         route%traversal%row(i), route%traversal%col(i), block, found)
         IF (.NOT. found) CYCLE
         nrow = route%row_size(route%traversal%row(i))
         ncol = route%col_size(route%traversal%col(i))
         CPASSERT(SIZE(block, 1) == nrow .AND. SIZE(block, 2) == ncol)
         pos = route%send_first(iblock)
         DO col = 1, ncol
            DO row = 1, nrow
               route%send_values(pos) = block(row, col)
               pos = pos + 1
            END DO
         END DO
      END DO
!$OMP END PARALLEL DO
      IF (kpoints%para_env_kp%num_pe == 1) THEN
         CALL world%allgatherv(route%send_values, snapshot%values, route%rcount(2, :), route%rdisp(2, :))
      ELSE
         CALL world%alltoall(route%send_values, route%scount(2, :), route%sdisp(2, :), &
                             snapshot%values, route%rcount(2, :), route%rdisp(2, :))
         CALL inter%bcast(snapshot%values, 0)
      END IF
      DO iblock = 1, SIZE(route%entry, 2)
         snapshot%nonzero(iblock) = ANY(snapshot%values(route%first(iblock):route%first(iblock + 1) - 1) /= 0.0_dp)
      END DO
      snapshot%source => source
      snapshot%matrix_row = matrix_row
      CALL timestop(handle)

   END SUBROUTINE kpoint_operator_rs_capture

! **************************************************************************************************
!> \brief Build this group's k point into its CFM tiles without MPI or BLACS.
!>        Columns are independent OpenMP work items. The snapshot's phase
!>        workspace is reused, so builds through the same context are serialized.
!> \param snapshot prepared real-space contributions
!> \param route structural block routing and target tile map
!> \param xkp the owning group's actual k-point coordinates
!> \param TARGET optional complex matrix on the prepared group distribution
!> \param fm_re optional real output, with fm_im instead of target
!> \param fm_im optional imaginary output on the same distribution
! **************************************************************************************************
   SUBROUTINE kpoint_operator_rs_build(snapshot, route, xkp, TARGET, fm_re, fm_im)

      TYPE(kpoint_operator_rs_type), INTENT(INOUT)       :: snapshot
      TYPE(kpoint_operator_rs_route_type), INTENT(IN)    :: route
      REAL(KIND=dp), INTENT(IN)                          :: xkp(3)
      TYPE(cp_cfm_type), INTENT(INOUT), OPTIONAL         :: target
      TYPE(cp_fm_type), INTENT(INOUT), OPTIONAL          :: fm_re, fm_im

      CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_operator_rs_build'

      COMPLEX(KIND=dp)                                   :: value
      INTEGER                                            :: c0, col, handle, iblock, j, local_row, &
                                                            nc, nr, pos, r0, row
      REAL(KIND=dp)                                      :: arg

      CALL timeset(routineN, handle)
      DO iblock = 1, SIZE(route%entry, 2)
         IF (.NOT. snapshot%nonzero(iblock)) CYCLE
         arg = twopi*DOT_PRODUCT(REAL(route%entry(5:7, iblock), dp), xkp)
         snapshot%phase(iblock) = CMPLX(COS(arg), SIN(arg), KIND=dp)
      END DO
      IF (PRESENT(TARGET)) THEN
         target%local_data(:, :) = z_zero
      ELSE
         CPASSERT(PRESENT(fm_re) .AND. PRESENT(fm_im))
         fm_re%local_data(:, :) = 0.0_dp
         fm_im%local_data(:, :) = 0.0_dp
      END IF
!$OMP PARALLEL DO DEFAULT(NONE) SHARED(snapshot, route, target, fm_re, fm_im) &
!$OMP PRIVATE(j, col, iblock, row, local_row, r0, c0, nr, nc, pos, value) SCHEDULE(STATIC)
      DO j = 1, SIZE(route%matrix_struct%col_indices)
         col = route%matrix_struct%col_indices(j)
         DO iblock = 1, SIZE(route%entry, 2)
            IF (.NOT. snapshot%nonzero(iblock)) CYCLE
            r0 = route%entry(1, iblock)
            c0 = route%entry(2, iblock)
            nr = route%entry(3, iblock)
            nc = route%entry(4, iblock)
            IF (col >= c0 .AND. col < c0 + nc) THEN
               pos = route%first(iblock) + (col - c0)*nr
               DO row = r0, r0 + nr - 1
                  local_row = route%row_local(row)
                  IF (local_row == 0) CYCLE
                  value = snapshot%phase(iblock)*snapshot%values(pos + row - r0)
                  IF (PRESENT(TARGET)) THEN
                     target%local_data(local_row, j) = target%local_data(local_row, j) + value
                  ELSE
                     fm_re%local_data(local_row, j) = fm_re%local_data(local_row, j) + REAL(value, dp)
                     fm_im%local_data(local_row, j) = fm_im%local_data(local_row, j) + AIMAG(value)
                  END IF
               END DO
            END IF
            IF (route%entry(8, iblock) == 0) CYCLE
            IF (col >= r0 .AND. col < r0 + nr) THEN
               pos = route%first(iblock) + col - r0
               DO row = c0, c0 + nc - 1
                  local_row = route%row_local(row)
                  IF (local_row == 0) CYCLE
                  value = CONJG(snapshot%phase(iblock))*snapshot%values(pos + (row - c0)*nr)
                  IF (PRESENT(TARGET)) THEN
                     target%local_data(local_row, j) = target%local_data(local_row, j) + value
                  ELSE
                     fm_re%local_data(local_row, j) = fm_re%local_data(local_row, j) + REAL(value, dp)
                     fm_im%local_data(local_row, j) = fm_im%local_data(local_row, j) + AIMAG(value)
                  END IF
               END DO
            END IF
         END DO
      END DO
!$OMP END PARALLEL DO
      CALL timestop(handle)

   END SUBROUTINE kpoint_operator_rs_build

! **************************************************************************************************
!> \brief Start one operator for one global round and spin. The service walks
!>        all groups with a k point in this round; groups without one still
!>        participate in the transfers. Every rank must iterate through the
!>        maximum local k-point count. Finish every started instance with
!>        kpoint_operator_finish. The first start also creates the slot
!>        bookkeeping of the context.
!>        Collective over the environment-wide communicator: every rank must
!>        call the same sequence of service calls.
!> \param ctx context created by kpoint_operator_context_create
!> \param ikp global round, up to the maximum k-point count across groups
!> \param ispin instance spin: keys the bookkeeping
!> \param rsmat real-space image matrices of the operator. A prepared operator
!>        batch uses the values captured by prepare; otherwise start reads
!>        the current values. The matching finish never reads rsmat
!> \param slot optional prepared slot: opt in to its immutable snapshot.
!>        Without slot, get always reads current source values. A missing or
!>        mismatched snapshot retains the original transfer path. Prepared FM
!>        targets must match the layout selected at prepare; inactive groups
!>        do not inspect or modify their output arguments.
!> \param grid prepared reciprocal-grid cache of the caller
!> \param use_grid extract from grid instead of a direct phase sum
!> \param matrix_row row of rsmat to transform. Spin-free operators such
!>        as S and T must pass kpoint_spin_free: their first matrix
!>        dimension counts derivatives, not spins. A spin index is out of
!>        bounds in plain SCF, but with derivatives stored it reads another
!>        valid row silently
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

      INTEGER                                            :: handle, igroup, ik, irs, kpmax
      INTEGER, DIMENSION(:, :), POINTER                  :: kp_dist

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(ctx%kpoints))
      CPASSERT(ASSOCIATED(ctx%kpoints%sab_nl))
      ! The instance spin is bounded by the spin space declared at creation;
      ! the matrix row by this operator's own rows. The two bounds differ for
      ! spin-free operators such as S, whose first dimension counts
      ! real-space blocks.
      CPASSERT(ispin >= 1 .AND. ispin <= ctx%nspin)
      CPASSERT(matrix_row >= 1 .AND. matrix_row <= SIZE(rsmat, 1))
      kpmax = MAXVAL(ctx%kpoints%kp_dist(2, :) - ctx%kpoints%kp_dist(1, :) + 1)
      CPASSERT(ikp >= 1 .AND. ikp <= kpmax)
      CPASSERT(slot >= kpoint_slot_ks .AND. slot <= kpoint_slot_t)
      ! State is per round; allocate expensive copy descriptors only when
      ! an actual distributed transfer is launched.
      IF (.NOT. ALLOCATED(ctx%status)) THEN
         ALLOCATE (ctx%status(kpmax, ctx%nspin, kpoint_num_slots))
         ctx%status(:, :, :) = kpoint_slot_idle
      END IF
      CPASSERT(ctx%status(ikp, ispin, slot) == kpoint_slot_idle)
      IF (.NOT. use_grid) THEN
         irs = kpoint_operator_rs_match(ctx, rsmat, ispin, slot, matrix_row)
         IF (irs > 0) THEN
            CPASSERT(ASSOCIATED(ctx%work_sab, ctx%kpoints%sab_nl))
            ctx%status(ikp, ispin, slot) = kpoint_slot_local
            CALL timestop(handle)
            RETURN
         END IF
      END IF
      CALL kpoint_operator_ensure(ctx, rsmat)

      ! Only the transfer backend needs copy descriptors. A prepared-only
      ! or get-only context never carries them. Size by the declared spin
      ! space, not by the possibly spin-free operator's row count.
      IF (.NOT. ALLOCATED(ctx%info)) THEN
         ! One copy_info per in-flight part: start takes it INTENT(OUT),
         ! so two concurrent starts cannot share one instance.
         ALLOCATE (ctx%info(SIZE(ctx%kpoints%xkp, 2), ctx%nspin, &
                            kpoint_num_slots, kpoint_part_im))
      END IF

      kp_dist => ctx%kpoints%kp_dist
      DO igroup = 1, ctx%kpoints%nkp_groups
         ik = kp_dist(1, igroup) + ikp - 1
         IF (ik > kp_dist(2, igroup)) CYCLE
         CALL kpoint_operator_densify(ctx, rsmat, matrix_row, ik, grid, use_grid)
         CALL kpoint_operator_launch(ctx%kpoints, ik, ctx%fmwork, ctx%fmlocal, &
                                     ctx%info(ik, ispin, slot, kpoint_part_re), &
                                     ctx%info(ik, ispin, slot, kpoint_part_im))
      END DO
      ctx%status(ikp, ispin, slot) = kpoint_slot_in_flight

      CALL timestop(handle)

   END SUBROUTINE kpoint_operator_start

! **************************************************************************************************
!> \brief Finish one started instance: build its prepared group-local operator,
!>        or merge received parts into cmat (real part first, overwrite then add)
!>        and reclaim the send side on every group.
!>        The slot returns to idle.
!>        Collective over the environment-wide communicator: every rank must
!>        call the same sequence of service calls.
!> \param ctx context that started the instance
!> \param ikp global round of the start call
!> \param ispin spin component of the start call
!> \param slot slot of the start call; the slot must be in flight
!> \param cmat caller-allocated complex matrix that receives the operator.
!>        It must share the distribution of the k-point group. Unchanged on
!>        groups without a k point in this round
! **************************************************************************************************
   SUBROUTINE kpoint_operator_finish(ctx, ikp, ispin, slot, cmat)

      TYPE(kpoint_operator_context_type), INTENT(INOUT)  :: ctx
      INTEGER, INTENT(IN)                                :: ikp, ispin, slot
      TYPE(cp_cfm_type)                                  :: cmat

      CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_operator_finish'

      INTEGER                                            :: handle, igroup, ik, irs, kpmax
      INTEGER, DIMENSION(2)                              :: kp_range
      INTEGER, DIMENSION(:, :), POINTER                  :: kp_dist
      LOGICAL                                            :: layout_ok, my_kpgrp

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(ctx%kpoints))
      CPASSERT(ALLOCATED(ctx%status))
      kpmax = SIZE(ctx%status, 1)
      CPASSERT(ikp >= 1 .AND. ikp <= kpmax)
      CPASSERT(ispin >= 1 .AND. ispin <= ctx%nspin)
      CPASSERT(slot >= kpoint_slot_ks .AND. slot <= kpoint_slot_t)
      CPASSERT(ctx%status(ikp, ispin, slot) /= kpoint_slot_idle)
      ! The merge writes cmat in place. A foreign distribution would fail
      ! inside the BLACS calls, far from this argument.

      kp_range = ctx%kpoints%kp_range
      IF (ctx%status(ikp, ispin, slot) == kpoint_slot_local) THEN
         layout_ok = cp_fm_struct_equivalent(cmat%matrix_struct, ctx%rs_route%matrix_struct)
         CPASSERT(layout_ok)
         irs = ispin
         IF (.NOT. ASSOCIATED(ctx%rs(irs, slot)%source)) irs = 1
         ik = kp_range(1) + ikp - 1
         IF (ik <= kp_range(2)) THEN
            CALL kpoint_operator_rs_build(ctx%rs(irs, slot), ctx%rs_route, ctx%kpoints%xkp(:, ik), cmat)
         END IF
         ctx%status(ikp, ispin, slot) = kpoint_slot_idle
         CALL timestop(handle)
         RETURN
      END IF
      CPASSERT(ctx%status(ikp, ispin, slot) == kpoint_slot_in_flight)
      layout_ok = cp_fm_struct_equivalent(cmat%matrix_struct, ctx%fmlocal%matrix_struct)
      CPASSERT(layout_ok)
      kp_dist => ctx%kpoints%kp_dist
      DO igroup = 1, ctx%kpoints%nkp_groups
         ik = kp_dist(1, igroup) + ikp - 1
         IF (ik > kp_dist(2, igroup)) CYCLE
         my_kpgrp = (ik >= kp_range(1) .AND. ik <= kp_range(2))
         IF (my_kpgrp) THEN
            CALL kpoint_operator_raw_finish_cfm(ctx%fmlocal, ctx%info(ik, ispin, slot, kpoint_part_re), &
                                                cmat, .FALSE.)
            CALL kpoint_operator_raw_finish_cfm(ctx%fmlocal, ctx%info(ik, ispin, slot, kpoint_part_im), &
                                                cmat, .TRUE.)
         END IF
         CALL cp_fm_cleanup_copy_general(ctx%info(ik, ispin, slot, kpoint_part_re))
         CALL cp_fm_cleanup_copy_general(ctx%info(ik, ispin, slot, kpoint_part_im))
      END DO
      ctx%status(ikp, ispin, slot) = kpoint_slot_idle

      CALL timestop(handle)

   END SUBROUTINE kpoint_operator_finish

! **************************************************************************************************
!> \brief Assemble one operator for one global round and spin in one call:
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
!> \param ikp global round, up to the maximum k-point count across groups
!> \param ispin instance spin: bounds-checked against the context spin space
!> \param rsmat real-space image matrices of the operator
!> \param cmat complex matrix output
!> \param fm_re real-part matrix output
!> \param fm_im imaginary-part matrix output
!> \param mat_re real-part block-cyclic dbcsr output
!> \param mat_im imaginary-part block-cyclic dbcsr output
!> \param matrix_row row of rsmat to transform. Spin-free operators such
!>        as S and T must pass kpoint_spin_free: their first matrix
!>        dimension counts derivatives, not spins. A spin index is out of
!>        bounds in plain SCF, but with derivatives stored it reads another
!>        valid row silently
!> \param slot optional prepared slot: opt in to its immutable snapshot.
!>        Without slot, get always reads current source values. A missing or
!>        mismatched snapshot retains the original transfer path. Prepared FM
!>        targets must match the layout selected at prepare; inactive groups
!>        do not inspect or modify their output arguments.
! **************************************************************************************************
   SUBROUTINE kpoint_operator_get(ctx, ikp, ispin, rsmat, cmat, fm_re, fm_im, &
                                  mat_re, mat_im, matrix_row, slot)

      TYPE(kpoint_operator_context_type), INTENT(INOUT)  :: ctx
      INTEGER, INTENT(IN)                                :: ikp, ispin
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rsmat
      TYPE(cp_cfm_type), OPTIONAL                        :: cmat
      TYPE(cp_fm_type), OPTIONAL                         :: fm_re, fm_im
      TYPE(dbcsr_type), INTENT(INOUT), OPTIONAL          :: mat_re, mat_im
      INTEGER, INTENT(IN)                                :: matrix_row
      INTEGER, INTENT(IN), OPTIONAL                      :: slot

      CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_operator_get'

      INTEGER                                            :: handle, igroup, ik, irs, kpmax
      INTEGER, DIMENSION(2)                              :: kp_range
      INTEGER, DIMENSION(:, :), POINTER                  :: kp_dist
      LOGICAL                                            :: layout_ok, my_kpgrp
      TYPE(copy_info_type)                               :: info_im, info_re

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
      kpmax = MAXVAL(ctx%kpoints%kp_dist(2, :) - ctx%kpoints%kp_dist(1, :) + 1)
      CPASSERT(ikp >= 1 .AND. ikp <= kpmax)
      IF (PRESENT(slot)) THEN
         CPASSERT(slot >= kpoint_slot_ks .AND. slot <= kpoint_slot_t)
         irs = kpoint_operator_rs_match(ctx, rsmat, ispin, slot, matrix_row)
         IF (irs > 0) THEN
            ik = ctx%kpoints%kp_range(1) + ikp - 1
            IF (ik <= ctx%kpoints%kp_range(2)) THEN
               IF (PRESENT(cmat)) THEN
                  layout_ok = cp_fm_struct_equivalent(cmat%matrix_struct, ctx%rs_route%matrix_struct)
                  CPASSERT(layout_ok)
                  CALL kpoint_operator_rs_build(ctx%rs(irs, slot), ctx%rs_route, ctx%kpoints%xkp(:, ik), cmat)
               ELSE IF (PRESENT(fm_re)) THEN
                  layout_ok = cp_fm_struct_equivalent(fm_re%matrix_struct, ctx%rs_route%matrix_struct)
                  CPASSERT(layout_ok)
                  layout_ok = cp_fm_struct_equivalent(fm_im%matrix_struct, ctx%rs_route%matrix_struct)
                  CPASSERT(layout_ok)
                  CALL kpoint_operator_rs_build(ctx%rs(irs, slot), ctx%rs_route, ctx%kpoints%xkp(:, ik), &
                                                fm_re=fm_re, fm_im=fm_im)
               ELSE
                  layout_ok = cp_fm_struct_equivalent(ctx%fmlocal%matrix_struct, ctx%rs_route%matrix_struct)
                  CPASSERT(layout_ok)
                  IF (.NOT. ALLOCATED(ctx%fmimag)) THEN
                     ALLOCATE (ctx%fmimag)
                     CALL cp_fm_create(ctx%fmimag, ctx%fmlocal%matrix_struct)
                  END IF
                  CALL kpoint_operator_rs_build(ctx%rs(irs, slot), ctx%rs_route, ctx%kpoints%xkp(:, ik), &
                                                fm_re=ctx%fmlocal, fm_im=ctx%fmimag)
                  CALL copy_fm_to_dbcsr_bc(ctx%fmlocal, mat_re)
                  CALL copy_fm_to_dbcsr_bc(ctx%fmimag, mat_im)
               END IF
            END IF
            CALL timestop(handle)
            RETURN
         END IF
      END IF
      CALL kpoint_operator_ensure(ctx, rsmat)
      IF (PRESENT(cmat)) THEN
         ! The merge writes cmat in place. A foreign distribution would fail
         ! inside the BLACS calls, far from this argument.
         layout_ok = cp_fm_struct_equivalent(cmat%matrix_struct, ctx%fmlocal%matrix_struct)
         CPASSERT(layout_ok)
      END IF

      kp_range = ctx%kpoints%kp_range
      kp_dist => ctx%kpoints%kp_dist
      DO igroup = 1, ctx%kpoints%nkp_groups
         ik = kp_dist(1, igroup) + ikp - 1
         IF (ik > kp_dist(2, igroup)) CYCLE
         my_kpgrp = (ik >= kp_range(1) .AND. ik <= kp_range(2))
         ! Absent fm targets are forwarded as absent optional arguments. Both
         ! CFM and DBCSR outputs use the same context-owned receive buffer.
         CALL kpoint_operator_densify(ctx, rsmat, matrix_row, ik, use_grid=.FALSE.)
         CALL kpoint_operator_launch(ctx%kpoints, ik, ctx%fmwork, ctx%fmlocal, info_re, info_im, &
                                     target_re=fm_re, target_im=fm_im)
         IF (my_kpgrp) THEN
            IF (PRESENT(cmat)) THEN
               CALL kpoint_operator_raw_finish_cfm(ctx%fmlocal, info_re, cmat, .FALSE.)
               CALL kpoint_operator_raw_finish_cfm(ctx%fmlocal, info_im, cmat, .TRUE.)
            ELSE IF (PRESENT(fm_re)) THEN
               CALL cp_fm_finish_copy_general(fm_re, info_re)
               CALL cp_fm_finish_copy_general(fm_im, info_im)
            ELSE
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
