!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Assembly of real and complex k-point operators from real-space DBCSR matrices.
!>        SCF slots and spin rows adapt the generic kpoint_transform machinery.
!>        Prepare captures values collectively; evaluation reads immutable snapshots
!>        into group-local or explicitly prepared world-distributed destinations.
!>        Budget-limited GEMM slots use collective single-k transforms and transfers.
!>        Start/finish retain the round bookkeeping used by SCF drivers.
!>        Serialize context mutation and collective preparation. Parallel evaluation
!>        requires private destinations and scratch. Release before the host environment.
! **************************************************************************************************
MODULE qs_kpoint_operators
   USE cp_cfm_types,                    ONLY: cp_cfm_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_init_p,&
                                              dbcsr_p_type,&
                                              dbcsr_type
   USE cp_dbcsr_operations,             ONLY: copy_fm_to_dbcsr_bc
   USE cp_fm_pool_types,                ONLY: cp_fm_pool_type,&
                                              fm_pool_create_fm,&
                                              fm_pool_give_back_fm
   USE cp_fm_struct,                    ONLY: cp_fm_struct_equivalent,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE kpoint_transform,                ONLY: &
        cell_to_k_backend_evaluate_batch, cell_to_k_backend_prepare, cell_to_k_sparse_fft, &
        rskp_backend_fft, rskp_backend_gemm, rskp_batch_size, rskp_capture, rskp_compatible, &
        rskp_prepared_transform_type, rskp_route_create, rskp_route_release, rskp_route_type, &
        rskp_transform_evaluate, rskp_transform_release, rskp_transform_stream_prepare, &
        rskp_values_release, rskp_values_type
   USE kpoint_types,                    ONLY: kpoint_env_type,&
                                              kpoint_type
   USE message_passing,                 ONLY: mp_para_env_type
   USE qs_matrix_pools,                 ONLY: mpools_get
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type

!$ USE omp_lib, ONLY: omp_get_max_threads, omp_get_num_threads
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_kpoint_operators'

   PUBLIC :: kpoint_operator_context_create, kpoint_operator_context_release
   PUBLIC :: kpoint_operator_start, kpoint_operator_finish
   PUBLIC :: kpoint_operator_get, kpoint_operator_prepare
   PUBLIC :: kpoint_operator_evaluate
   PUBLIC :: kpoint_operator_batch_size, kpoint_operator_evaluate_batch
   PUBLIC :: kpoint_operator_store, kpoint_operator_get_local, kpoint_operator_get_all

   INTEGER, PARAMETER, PUBLIC :: kpoint_slot_ks = 1, &
                                 kpoint_slot_s = 2, &
                                 kpoint_slot_t = 3

   !> Row of rsmat that carries the spin-free operator image. S and T carry
   !> no spin copies: drivers launch one instance per spin of their loop and
   !> read this same row in every instance.
   INTEGER, PARAMETER, PUBLIC :: kpoint_spin_free = 1

   !> Number of defined slots; grows only inside this module.
   INTEGER, PARAMETER, PRIVATE :: kpoint_num_slots = kpoint_slot_t - kpoint_slot_ks + 1

   !> Instance states of one (k point, spin, slot) round.
   !> State transitions happen only inside the service routines.
   INTEGER, PARAMETER, PRIVATE :: kpoint_slot_idle = 0, &
                                  kpoint_slot_local = 1

   ! Structural element/image route shared by captured GEMM/FFT rows. The
   ! borrowed source neighbor lists must remain alive and unchanged with it.
   ! Immutable captured values; FFT slots also own the local reciprocal matrices.
   ! Shared metadata belong to the route; source identifies explicitly prepared starts.
   INTEGER(KIND=int_8), PARAMETER :: rs_snapshot_bytes = 64_int_8*1024_int_8**2
   ! Temporary benchmark allowance: shared phases, worker scratch and output tiles.
   INTEGER(KIND=int_8), PARAMETER :: rs_batch_bytes = 64_int_8*1024_int_8**2

! **************************************************************************************************
!> \brief Caller-owned assembly context for a fixed geometry and parallel layout.
!>        OT retains it for one inner SCF loop; other drivers use one call.
!>        Release before changing the neighbor lists, matrix pools or k-point groups.
!>        Components are private and touched only by this module.
! **************************************************************************************************
   TYPE, PUBLIC :: kpoint_operator_context_type
      PRIVATE
      INTEGER :: nspin = 0
      LOGICAL :: overlap_stored = .FALSE.
      TYPE(kpoint_type), POINTER :: kpoints => NULL()
      TYPE(cp_fm_pool_type), POINTER :: group_fm_pool => NULL()
      TYPE(cp_fm_type), ALLOCATABLE :: fmlocal, fmimag
      TYPE(rskp_values_type), ALLOCATABLE :: rs(:, :)
      TYPE(rskp_route_type), ALLOCATABLE :: rs_route
      TYPE(rskp_prepared_transform_type), POINTER :: stream(:) => NULL()
      LOGICAL :: streaming(kpoint_num_slots) = .FALSE.
      INTEGER, ALLOCATABLE :: status(:, :, :)
   END TYPE kpoint_operator_context_type

CONTAINS

! **************************************************************************************************
!> \brief Create an empty operator context.
!>        The context owns no buffers after this call; they are created on first
!>        use and released with the context.
!> \param ctx context to initialize
!> \param kpoints host k-point environment, referenced, not copied
!> \param ao_ao_fm optional legacy context template, checked for communicator compatibility only
!> \param nspin instance spin space of the owning driver call: the bound of
!>        the spin loop its start and finish calls iterate. Fixed for the
!>        lifetime of the context
! **************************************************************************************************
   SUBROUTINE kpoint_operator_context_create(ctx, kpoints, ao_ao_fm, nspin)

      TYPE(kpoint_operator_context_type), INTENT(OUT)    :: ctx
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(cp_fm_type), INTENT(IN), OPTIONAL             :: ao_ao_fm
      INTEGER, INTENT(IN)                                :: nspin

      CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_operator_context_create'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CPASSERT(nspin >= 1)
      ctx%nspin = nspin
      CPASSERT(ASSOCIATED(kpoints))
      ctx%kpoints => kpoints
      IF (PRESENT(ao_ao_fm)) THEN
         IF (.NOT. ASSOCIATED(ao_ao_fm%matrix_struct%para_env, kpoints%blacs_env_all%para_env)) THEN
            CPABORT("R-to-K template must use the world communicator")
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE kpoint_operator_context_create

! **************************************************************************************************
!> \brief Release the context together with every buffer it created.
!>        Call once, at the end of the fixed-layout scope that created the context.
!>        Every started instance must have been finished.
!> \param ctx context to release
! **************************************************************************************************
   SUBROUTINE kpoint_operator_context_release(ctx)

      TYPE(kpoint_operator_context_type), INTENT(INOUT)  :: ctx

      CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_operator_context_release'

      INTEGER                                            :: handle, ispin, slot

      CALL timeset(routineN, handle)

      IF (ALLOCATED(ctx%status)) THEN
         CPASSERT(ALL(ctx%status == kpoint_slot_idle))
      END IF

      IF (ALLOCATED(ctx%rs)) THEN
         DO slot = 1, kpoint_num_slots
            DO ispin = 1, ctx%nspin
               CALL rskp_values_release(ctx%rs(ispin, slot))
            END DO
         END DO
         DEALLOCATE (ctx%rs)
      END IF
      IF (ASSOCIATED(ctx%stream)) THEN
         DO slot = 1, kpoint_num_slots
            CALL rskp_transform_release(ctx%stream(slot))
         END DO
         DEALLOCATE (ctx%stream)
      END IF
      ctx%streaming = .FALSE.
      IF (ALLOCATED(ctx%rs_route)) THEN
         CALL rskp_route_release(ctx%rs_route)
         DEALLOCATE (ctx%rs_route)
      END IF
      IF (ALLOCATED(ctx%fmimag)) THEN
         CALL cp_fm_release(ctx%fmimag)
         DEALLOCATE (ctx%fmimag)
      END IF
      IF (ALLOCATED(ctx%fmlocal)) THEN
         CALL fm_pool_give_back_fm(ctx%group_fm_pool, ctx%fmlocal)
         DEALLOCATE (ctx%fmlocal)
      END IF
      IF (ALLOCATED(ctx%status)) DEALLOCATE (ctx%status)
      NULLIFY (ctx%kpoints, ctx%group_fm_pool)
      ctx%overlap_stored = .FALSE.

      CALL timestop(handle)

   END SUBROUTINE kpoint_operator_context_release

! **************************************************************************************************
!> \brief Acquire only the group-local output workspace; both backends share it.
!> \param ctx operator context
! **************************************************************************************************
   SUBROUTINE kpoint_operator_ensure_group(ctx)

      TYPE(kpoint_operator_context_type), INTENT(INOUT)  :: ctx

      IF (.NOT. ALLOCATED(ctx%fmlocal)) THEN
         CALL mpools_get(ctx%kpoints%mpools, ao_ao_fm_pool=ctx%group_fm_pool)
         ALLOCATE (ctx%fmlocal)
         CALL fm_pool_create_fm(ctx%group_fm_pool, ctx%fmlocal)
      END IF

   END SUBROUTINE kpoint_operator_ensure_group

! **************************************************************************************************
!> \brief Prepare an immutable operator slot for group-local assembly.
!>        Refresh after the source values change, with every instance idle.
!>        Other prepared slots and structural routing remain valid. Without
!>        matrix_row, capture one row per instance spin; with matrix_row, store
!>        one row shared by all instance spins. mirror_sign controls the
!>        transpose relation of the transformed operator.
!> \param ctx context spanning the current initialized k-point environment
!> \param rsmat current real-space operator images
!> \param slot slot to replace (KS, S or T)
!> \param prepared true for independent cached evaluation; false for collective streaming
!> \param matrix_row optional explicit spin-free operator row
!> \param TARGET optional FM object defining the prepared-output layout
!> \param target_struct optional FM structure defining the prepared-output
!>        layout, as an alternative to TARGET
!> \param sab_nl optional source neighbor lists; defaults to the host AO lists
!> \param mirror_sign optional transpose sign (+1 Hermitian, -1 anti-Hermitian)
! **************************************************************************************************
   SUBROUTINE kpoint_operator_prepare(ctx, rsmat, slot, prepared, matrix_row, TARGET, target_struct, sab_nl, &
                                      mirror_sign)

      TYPE(kpoint_operator_context_type), INTENT(INOUT)  :: ctx
      TYPE(dbcsr_p_type), POINTER                        :: rsmat(:, :)
      INTEGER, INTENT(IN)                                :: slot
      LOGICAL, INTENT(OUT), OPTIONAL                     :: prepared
      INTEGER, INTENT(IN), OPTIONAL                      :: matrix_row
      TYPE(cp_fm_type), INTENT(IN), OPTIONAL             :: target
      TYPE(cp_fm_struct_type), OPTIONAL, POINTER         :: target_struct
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         OPTIONAL, POINTER                               :: sab_nl
      INTEGER, INTENT(IN), OPTIONAL                      :: mirror_sign

      CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_operator_prepare'

      INTEGER                                            :: compatible, handle, ispin, jslot, nk, &
                                                            nrow, output_mirror_sign, row
      INTEGER(KIND=int_8)                                :: fft_storage, retained, storage, &
                                                            workspace_bytes
      INTEGER, ALLOCATABLE                               :: source_rows(:)
      LOGICAL                                            :: over_budget, ready, world_target
      TYPE(cp_fm_struct_type), POINTER                   :: fms
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: source_sab

      CALL timeset(routineN, handle)
      IF (PRESENT(prepared)) prepared = .FALSE.
      CPASSERT(ASSOCIATED(ctx%kpoints))
      CPASSERT(.NOT. (PRESENT(TARGET) .AND. PRESENT(target_struct)))
      output_mirror_sign = 1
      IF (PRESENT(mirror_sign)) output_mirror_sign = mirror_sign
      CPASSERT(ABS(output_mirror_sign) == 1)
      IF (PRESENT(sab_nl)) THEN
         CPASSERT(ASSOCIATED(sab_nl))
         source_sab => sab_nl
      ELSE
         source_sab => ctx%kpoints%sab_nl
      END IF
      CPASSERT(ASSOCIATED(source_sab))
      CPASSERT(slot >= kpoint_slot_ks .AND. slot <= kpoint_slot_t)
      IF (ALLOCATED(ctx%status)) THEN
         CPASSERT(ALL(ctx%status == kpoint_slot_idle))
      END IF
      IF (.NOT. ALLOCATED(ctx%rs)) THEN
         ALLOCATE (ctx%rs(ctx%nspin, kpoint_num_slots))
      END IF
      IF (.NOT. ASSOCIATED(ctx%stream)) ALLOCATE (ctx%stream(kpoint_num_slots))
      CALL rskp_transform_release(ctx%stream(slot))
      ctx%streaming(slot) = .FALSE.
      DO ispin = 1, ctx%nspin
         CALL rskp_values_release(ctx%rs(ispin, slot))
      END DO
      nrow = ctx%nspin
      IF (PRESENT(matrix_row)) THEN
         CPASSERT(matrix_row >= 1 .AND. matrix_row <= SIZE(rsmat, 1))
         nrow = 1
      ELSE
         CPASSERT(SIZE(rsmat, 1) == ctx%nspin)
      END IF
      para_env => ctx%kpoints%blacs_env_all%para_env
      IF (PRESENT(target_struct)) THEN
         CPASSERT(ASSOCIATED(target_struct))
         fms => target_struct
      ELSE IF (PRESENT(TARGET)) THEN
         fms => target%matrix_struct
      ELSE
         CALL kpoint_operator_ensure_group(ctx)
         fms => ctx%fmlocal%matrix_struct
      END IF
      CPASSERT(ASSOCIATED(fms))
      ! Group-local targets are routed through group leaders and broadcast
      ! within each group; world targets are routed directly to world ranks.
      world_target = ASSOCIATED(fms%para_env, para_env)
      compatible = MERGE(1, 0, world_target .OR. &
                         ASSOCIATED(fms%para_env, ctx%kpoints%para_env_kp))
      CALL para_env%min(compatible)
      IF (compatible == 0) THEN
         CPABORT("R-to-K output must use the host group or world communicator")
      END IF
      snapshot_preparation: BLOCK
         ready = .FALSE.
         IF (.NOT. ALLOCATED(ctx%rs_route)) THEN
            ALLOCATE (ctx%rs_route)
            CALL rskp_route_create(rsmat, source_sab, ctx%kpoints%cell_to_index, fms, &
                                   ctx%kpoints%blacs_env_all%para_env, ctx%kpoints%para_env_kp, &
                                   ctx%kpoints%para_env_inter_kp, ctx%kpoints%rskp_transform, &
                                   rs_snapshot_bytes, ctx%rs_route, ready, over_budget)
            IF (.NOT. ready) THEN
               IF (.NOT. over_budget .OR. ctx%kpoints%rskp_transform /= rskp_backend_gemm) THEN
                  CPABORT("R-to-K route is incompatible or forced FFT exceeds the memory budget")
               END IF
               CALL rskp_route_release(ctx%rs_route)
               DEALLOCATE (ctx%rs_route)
               EXIT snapshot_preparation
            END IF
            IF (ctx%rs_route%transform%backend == rskp_backend_fft .OR. ctx%rs_route%world_target) THEN
               CALL cell_to_k_backend_prepare(ctx%rs_route, ctx%kpoints%xkp, ctx%kpoints%nkp_grid, rs_batch_bytes)
            ELSE
               CALL cell_to_k_backend_prepare(ctx%rs_route, &
                                              ctx%kpoints%xkp(:, ctx%kpoints%kp_range(1):ctx%kpoints%kp_range(2)), &
                                              ctx%kpoints%nkp_grid, rs_batch_bytes)
            END IF
         END IF
         compatible = MERGE(1, 0, cp_fm_struct_equivalent(fms, ctx%rs_route%matrix_struct))
         IF (world_target .NEQV. ctx%rs_route%world_target) compatible = 0
         DO ispin = 1, nrow
            row = ispin
            IF (PRESENT(matrix_row)) row = matrix_row
            IF (.NOT. rskp_compatible(ctx%rs_route, rsmat, row, source_sab)) compatible = 0
         END DO
         CALL para_env%min(compatible)
         ! Fixed element/scatter metadata belong to the route; packed snapshots own only values.
         retained = 0_int_8
         DO jslot = kpoint_slot_ks, kpoint_slot_t
            DO ispin = 1, ctx%nspin
               IF (ALLOCATED(ctx%rs(ispin, jslot)%packed)) THEN
                  retained = retained + 8_int_8*SIZE(ctx%rs(ispin, jslot)%packed, KIND=int_8)
               END IF
               IF (ALLOCATED(ctx%rs(ispin, jslot)%active)) THEN
                  retained = retained + INT(STORAGE_SIZE(.FALSE.)/8, int_8)*SIZE(ctx%rs(ispin, jslot)%active, KIND=int_8)
               END IF
            END DO
         END DO
         storage = ctx%rs_route%storage + retained
         storage = storage + 8_int_8*nrow*( &
                   ctx%rs_route%transform%packed_first(SIZE(ctx%rs_route%transform%packed_first)) - 1)
         IF (ctx%rs_route%transform%backend == rskp_backend_gemm) THEN
            storage = storage + INT(STORAGE_SIZE(.FALSE.)/8, int_8)*4_int_8*nrow*SIZE(ctx%rs_route%transform%phase_count, 2)
         END IF
         CALL para_env%max(storage)
         IF (compatible == 0) THEN
            CPABORT("R-to-K snapshot has an incompatible layout")
         END IF
         ready = storage <= rs_snapshot_bytes
         IF (.NOT. ready) THEN
            IF (ctx%kpoints%rskp_transform /= rskp_backend_gemm) THEN
               CPABORT("FFT snapshot exceeds budget")
            END IF
            EXIT snapshot_preparation
         END IF
         IF (ready) THEN
            fft_storage = 0_int_8
            IF (ctx%rs_route%transform%backend == rskp_backend_fft) THEN
               DO jslot = kpoint_slot_ks, kpoint_slot_t
                  DO ispin = 1, ctx%nspin
                     IF (ALLOCATED(ctx%rs(ispin, jslot)%fft_values)) THEN
                        fft_storage = fft_storage + 16_int_8*SIZE(ctx%rs(ispin, jslot)%fft_values, KIND=int_8)
                     END IF
                  END DO
               END DO
               nk = SIZE(ctx%kpoints%xkp, 2)
               IF (.NOT. world_target) nk = SIZE(ctx%kpoints%kp_env)
               fft_storage = fft_storage + 16_int_8*nrow*nk*fms%local_leading_dimension*MAX(1, ctx%rs_route%ncol)
               workspace_bytes = 512_int_8*1024_int_8**2 - storage - fft_storage
               IF (workspace_bytes <= 0_int_8) THEN
                  CPABORT("Forced lattice FFT cache exceeds the memory budget")
               END IF
            END IF
            DO ispin = 1, nrow
               row = ispin
               IF (PRESENT(matrix_row)) row = matrix_row
               CALL rskp_capture(rsmat, row, ctx%rs_route, ctx%rs(ispin, slot), mirror_sign=output_mirror_sign)
               IF (ctx%rs_route%transform%backend == rskp_backend_fft) THEN
                  ALLOCATE (ctx%rs(ispin, slot)%fft_values(fms%local_leading_dimension, MAX(1, ctx%rs_route%ncol), nk))
                  CALL cell_to_k_sparse_fft(ctx%rs(ispin, slot), ctx%rs_route, ctx%kpoints%xkp, &
                                            ctx%kpoints%nkp_grid, ctx%kpoints%kp_dist, ctx%kpoints%para_env_inter_kp, &
                                            ctx%rs(ispin, slot)%fft_values, workspace_bytes)
                  ! The reciprocal cache no longer needs the replicated real-space values.
                  DEALLOCATE (ctx%rs(ispin, slot)%packed)
               END IF
            END DO
         END IF
      END BLOCK snapshot_preparation
      IF (.NOT. ready) THEN
         ALLOCATE (source_rows(nrow))
         DO ispin = 1, nrow
            row = ispin
            IF (PRESENT(matrix_row)) row = matrix_row
            source_rows(ispin) = row
            ctx%rs(ispin, slot)%source => rsmat
            ctx%rs(ispin, slot)%matrix_row = row
         END DO
         CALL rskp_transform_stream_prepare(ctx%stream(slot), rsmat, source_rows, source_sab, &
                                            ctx%kpoints%cell_to_index, fms, para_env, ctx%kpoints%para_env_kp, &
                                            ctx%kpoints%para_env_inter_kp, ctx%kpoints%kp_dist, &
                                            mirror_sign=output_mirror_sign)
         ctx%streaming(slot) = .TRUE.
      END IF
      IF (PRESENT(prepared)) prepared = ready
      CALL timestop(handle)

   END SUBROUTINE kpoint_operator_prepare

! **************************************************************************************************
!> \brief Match an explicitly prepared row, sharing spin-free rows across spins.
!> \param ctx operator context
!> \param source live array identity identifying the prepared operation
!> \param ispin instance spin
!> \param slot prepared slot
!> \param matrix_row operator row
!> \return snapshot spin index, or zero when the prepared source does not match
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
!> \brief Evaluate an immutable prepared operator with exclusively caller-owned output and scratch.
!>        Cached slots support independent worker evaluation. Streaming slots require
!>        collective rounds outside worker regions. Source blocks must stay unchanged.
!> \param ctx prepared read-only context
!> \param ispin instance spin
!> \param slot operator slot
!> \param ik global k index owned by this group
!> \param TARGET optional private complex singleton output
!> \param fm_re optional private real output instead of TARGET
!> \param fm_im optional private imaginary output, accompanying fm_re for a complex operator
!> \param round collective streaming round; k is a local round shared by all groups
! **************************************************************************************************
   SUBROUTINE kpoint_operator_evaluate(ctx, ispin, slot, ik, TARGET, fm_re, fm_im, round)
      TYPE(kpoint_operator_context_type), INTENT(IN)     :: ctx
      INTEGER, INTENT(IN)                                :: ispin, slot, ik
      TYPE(cp_cfm_type), INTENT(INOUT), OPTIONAL         :: target
      TYPE(cp_fm_type), INTENT(INOUT), OPTIONAL, TARGET  :: fm_re, fm_im
      LOGICAL, INTENT(IN), OPTIONAL                      :: round

      COMPLEX(KIND=dp), CONTIGUOUS, POINTER              :: batch(:, :, :)
      INTEGER                                            :: first, nc, nk, nworkers, requested, row, &
                                                            width
      INTEGER(KIND=int_8)                                :: work_size
      REAL(KIND=dp), ALLOCATABLE                         :: work(:)
      TYPE(cp_fm_p_type)                                 :: output(1, 2)

      CPASSERT(PRESENT(TARGET) .NEQV. PRESENT(fm_re))
      CPASSERT(.NOT. PRESENT(fm_im) .OR. PRESENT(fm_re))
      IF (ctx%streaming(slot)) THEN
         row = ispin
         IF (.NOT. ASSOCIATED(ctx%rs(row, slot)%source)) row = 1
         CALL rskp_transform_evaluate(ctx%stream(slot), row, ik, ctx%kpoints%xkp, &
                                      fm_re=fm_re, fm_im=fm_im, cfm=TARGET, round=round)
         IF (PRESENT(round)) THEN
            IF (round .AND. ik > SIZE(ctx%kpoints%kp_env)) RETURN
         END IF
         IF (ctx%kpoints%use_real_wfn) THEN
            IF (PRESENT(TARGET)) target%local_data(:, :) = CMPLX(REAL(target%local_data, dp), 0.0_dp, dp)
            IF (PRESENT(fm_im)) fm_im%local_data(:, :) = 0.0_dp
         END IF
         RETURN
      END IF
      IF (ctx%rs_route%world_target) THEN
         CPASSERT(ik >= 1 .AND. ik <= ctx%kpoints%nkp)
      ELSE
         CPASSERT(ik >= ctx%kpoints%kp_range(1))
         CPASSERT(ik <= ctx%kpoints%kp_range(2))
      END IF
      row = ispin
      IF (.NOT. ASSOCIATED(ctx%rs(row, slot)%source)) row = 1
      CPASSERT(ASSOCIATED(ctx%rs(row, slot)%source))
      requested = 1
!$    requested = omp_get_num_threads()
      nworkers = requested
      first = ik
      IF (.NOT. ctx%rs_route%world_target) first = ik - ctx%kpoints%kp_range(1) + 1
      IF (PRESENT(TARGET)) THEN
         CPASSERT(cp_fm_struct_equivalent(target%matrix_struct, ctx%rs_route%matrix_struct))
         CALL kpoint_operator_batch_size(ctx, nworkers, 1, width, work_size)
      ELSE
         nc = 1
         IF (PRESENT(fm_im) .OR. AIMAG(ctx%rs(row, slot)%scale) /= 0.0_dp) nc = 2
         nk = SIZE(ctx%kpoints%kp_env)
         IF (ctx%rs_route%world_target) nk = ctx%kpoints%nkp
         CALL rskp_batch_size(ctx%rs_route, nk, nc, nworkers, 0, rs_batch_bytes, width, work_size)
      END IF
      IF (nworkers < requested) THEN
         CPABORT("R-to-K singleton workers exceed the budget; use batched assembly")
      END IF
      ALLOCATE (work(work_size/width))
      IF (PRESENT(TARGET)) THEN
         batch(1:SIZE(target%local_data, 1), 1:SIZE(target%local_data, 2), 1:1) => target%local_data
         CALL kpoint_operator_evaluate_batch(ctx, row, slot, ik, batch, work)
      ELSE
         output(1, 1)%matrix => fm_re
         nc = 1
         IF (PRESENT(fm_im)) THEN
            output(1, 2)%matrix => fm_im
            nc = 2
         END IF
         CALL cell_to_k_backend_evaluate_batch(ctx%rs(row, slot), ctx%rs_route, &
                                               ctx%kpoints%xkp(:, ik:ik), first, work, output=output(:, :nc))
      END IF
      IF (ctx%kpoints%use_real_wfn .AND. AIMAG(ctx%rs(row, slot)%scale) == 0.0_dp) THEN
         IF (PRESENT(TARGET)) target%local_data(:, :) = CMPLX(REAL(target%local_data, dp), 0.0_dp, dp)
         IF (PRESENT(fm_im)) fm_im%local_data(:, :) = 0.0_dp
      END IF
   END SUBROUTINE kpoint_operator_evaluate

! **************************************************************************************************
!> \brief Mark a prepared operator for a round, including inactive groups.
!>        Finish every started instance before preparing another snapshot.
!> \param ctx context created by kpoint_operator_context_create
!> \param ikp global round, up to the maximum k-point count across groups
!> \param ispin instance spin: keys the bookkeeping
!> \param rsmat source identity of the immutable snapshot captured by prepare
!> \param slot prepared instance slot
!> \param matrix_row row of rsmat to transform. Spin-free operators such
!>        as S and T must pass kpoint_spin_free: their first matrix
!>        dimension counts derivatives, not spins. A spin index is out of
!>        bounds in plain SCF, but with derivatives stored it reads another
!>        valid row silently
! **************************************************************************************************
   SUBROUTINE kpoint_operator_start(ctx, ikp, ispin, rsmat, slot, matrix_row)

      TYPE(kpoint_operator_context_type), INTENT(INOUT)  :: ctx
      INTEGER, INTENT(IN)                                :: ikp, ispin
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rsmat
      INTEGER, INTENT(IN)                                :: slot, matrix_row

      CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_operator_start'

      INTEGER                                            :: handle, irs, kpmax

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
      ! State is per round; no communication is deferred to finish.
      IF (.NOT. ALLOCATED(ctx%status)) THEN
         ALLOCATE (ctx%status(kpmax, ctx%nspin, kpoint_num_slots))
         ctx%status(:, :, :) = kpoint_slot_idle
      END IF
      CPASSERT(ctx%status(ikp, ispin, slot) == kpoint_slot_idle)
      irs = kpoint_operator_rs_match(ctx, rsmat, ispin, slot, matrix_row)
      IF (irs == 0) THEN
         CPABORT("Operator start requires a matching prepared source row")
      END IF
      ctx%status(ikp, ispin, slot) = kpoint_slot_local

      CALL timestop(handle)

   END SUBROUTINE kpoint_operator_start

! **************************************************************************************************
!> \brief Evaluate a started instance on its owning group and return the slot to idle.
!> \param ctx context that started the instance
!> \param ikp global round of the start call
!> \param ispin spin component of the start call
!> \param slot slot of the start call
!> \param cmat caller-allocated complex matrix that receives the operator.
!>        It must share the distribution of the k-point group. Unchanged on
!>        groups without a k point in this round
!> \param fm_re optional real FM output instead of cmat for real wavefunctions, on the group layout
! **************************************************************************************************
   SUBROUTINE kpoint_operator_finish(ctx, ikp, ispin, slot, cmat, fm_re)

      TYPE(kpoint_operator_context_type), INTENT(INOUT)  :: ctx
      INTEGER, INTENT(IN)                                :: ikp, ispin, slot
      TYPE(cp_cfm_type), OPTIONAL                        :: cmat
      TYPE(cp_fm_type), OPTIONAL                         :: fm_re

      CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_operator_finish'

      INTEGER                                            :: handle, ik, irs, kpmax
      INTEGER, DIMENSION(2)                              :: kp_range

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(ctx%kpoints))
      CPASSERT(ALLOCATED(ctx%status))
      kpmax = SIZE(ctx%status, 1)
      CPASSERT(ikp >= 1 .AND. ikp <= kpmax)
      CPASSERT(ispin >= 1 .AND. ispin <= ctx%nspin)
      CPASSERT(slot >= kpoint_slot_ks .AND. slot <= kpoint_slot_t)
      CPASSERT(ctx%status(ikp, ispin, slot) /= kpoint_slot_idle)
      CPASSERT(PRESENT(cmat) .NEQV. PRESENT(fm_re))
      CPASSERT(.NOT. PRESENT(fm_re) .OR. ctx%kpoints%use_real_wfn)
      kp_range = ctx%kpoints%kp_range
      CPASSERT(ctx%status(ikp, ispin, slot) == kpoint_slot_local)
      irs = ispin
      IF (.NOT. ASSOCIATED(ctx%rs(irs, slot)%source)) irs = 1
      ik = kp_range(1) + ikp - 1
      IF (ctx%streaming(slot)) THEN
         CALL kpoint_operator_evaluate(ctx, irs, slot, ikp, cmat, fm_re=fm_re, round=.TRUE.)
      ELSE IF (ik <= kp_range(2)) THEN
         CALL kpoint_operator_evaluate(ctx, irs, slot, ik, cmat, fm_re=fm_re)
      END IF
      ctx%status(ikp, ispin, slot) = kpoint_slot_idle

      CALL timestop(handle)

   END SUBROUTINE kpoint_operator_finish

! **************************************************************************************************
!> \brief Evaluate one operator round into CFM, split FM or block-cyclic DBCSR.
!>        Destinations must match the prepared layout. DBCSR matrices must be
!>        freshly initialized with dbcsr_init_p. Inactive groups leave outputs unchanged.
!>        Without a slot, get collectively refreshes the source before evaluation.
!> \param ctx context created by kpoint_operator_context_create
!> \param ikp global round, up to the maximum k-point count across groups
!> \param ispin instance spin: bounds-checked against the context spin space
!> \param rsmat real-space image matrices of the operator
!> \param cmat complex matrix output
!> \param fm_re real-part matrix output
!> \param fm_im imaginary-part output; may be omitted for real wavefunctions
!> \param mat_re real-part block-cyclic dbcsr output
!> \param mat_im imaginary-part block-cyclic dbcsr output
!> \param matrix_row row of rsmat to transform. Spin-free operators such
!>        as S and T must pass kpoint_spin_free: their first matrix
!>        dimension counts derivatives, not spins. A spin index is out of
!>        bounds in plain SCF, but with derivatives stored it reads another
!>        valid row silently
!> \param slot optional prepared slot; a missing or mismatched snapshot aborts
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

      INTEGER                                            :: handle, ik, irs, kpmax, selected_slot

      CALL timeset(routineN, handle)

      CPASSERT(COUNT([PRESENT(cmat), PRESENT(fm_re), PRESENT(mat_re)]) == 1)
      IF (PRESENT(fm_re) .OR. PRESENT(fm_im)) THEN
         CPASSERT(PRESENT(fm_re))
         CPASSERT(PRESENT(fm_im) .OR. ctx%kpoints%use_real_wfn)
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
      selected_slot = kpoint_slot_ks
      IF (PRESENT(slot)) THEN
         selected_slot = slot
      ELSE
         CALL kpoint_operator_prepare(ctx, rsmat, selected_slot, matrix_row=matrix_row)
      END IF
      CPASSERT(selected_slot >= kpoint_slot_ks .AND. selected_slot <= kpoint_slot_t)
      irs = kpoint_operator_rs_match(ctx, rsmat, ispin, selected_slot, matrix_row)
      IF (irs == 0) THEN
         CPABORT("Operator get requires a matching prepared source row")
      END IF
      ik = ctx%kpoints%kp_range(1) + ikp - 1
      IF (ctx%streaming(selected_slot)) THEN
         IF (PRESENT(mat_re)) THEN
            CALL kpoint_operator_ensure_group(ctx)
            IF (.NOT. ALLOCATED(ctx%fmimag)) THEN
               ALLOCATE (ctx%fmimag)
               CALL cp_fm_create(ctx%fmimag, ctx%fmlocal%matrix_struct)
            END IF
            CALL kpoint_operator_evaluate(ctx, irs, selected_slot, ikp, &
                                          fm_re=ctx%fmlocal, fm_im=ctx%fmimag, round=.TRUE.)
            IF (ik <= ctx%kpoints%kp_range(2)) THEN
               CALL copy_fm_to_dbcsr_bc(ctx%fmlocal, mat_re)
               CALL copy_fm_to_dbcsr_bc(ctx%fmimag, mat_im)
            END IF
         ELSE
            CALL kpoint_operator_evaluate(ctx, irs, selected_slot, ikp, &
                                          TARGET=cmat, fm_re=fm_re, fm_im=fm_im, round=.TRUE.)
         END IF
      ELSE IF (ik <= ctx%kpoints%kp_range(2)) THEN
         IF (PRESENT(mat_re)) THEN
            CALL kpoint_operator_ensure_group(ctx)
            IF (.NOT. ALLOCATED(ctx%fmimag)) THEN
               ALLOCATE (ctx%fmimag)
               CALL cp_fm_create(ctx%fmimag, ctx%fmlocal%matrix_struct)
            END IF
            CALL kpoint_operator_evaluate(ctx, irs, selected_slot, ik, fm_re=ctx%fmlocal, fm_im=ctx%fmimag)
            CALL copy_fm_to_dbcsr_bc(ctx%fmlocal, mat_re)
            CALL copy_fm_to_dbcsr_bc(ctx%fmimag, mat_im)
         ELSE
            CALL kpoint_operator_evaluate(ctx, irs, selected_slot, ik, TARGET=cmat, fm_re=fm_re, fm_im=fm_im)
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE kpoint_operator_get

! **************************************************************************************************
!> \brief Read a cached operator inside a group-local OT loop.
!>        Every group layout requires collective cache preparation before entering the loop.
!>        This routine copies the prepared cache, without recomputing the real-space transform.
!> \param matrix_rs current real-space operator
!> \param kpoints host environment
!> \param kp local k-point state
!> \param ispin source row (one for spin-free S/T)
!> \param cache_re prepared real part on the k-point group
!> \param cache_im prepared imaginary part on the k-point group
!> \param matrix_re newly initialized sparse real output
!> \param matrix_im newly initialized sparse imaginary output
! **************************************************************************************************
   SUBROUTINE kpoint_operator_get_local(matrix_rs, kpoints, kp, ispin, cache_re, cache_im, &
                                        matrix_re, matrix_im)
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_rs
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(kpoint_env_type), POINTER                     :: kp
      INTEGER, INTENT(IN)                                :: ispin
      TYPE(cp_fm_type), INTENT(IN)                       :: cache_re, cache_im
      TYPE(dbcsr_type), POINTER                          :: matrix_re, matrix_im

      CPASSERT(ASSOCIATED(matrix_rs))
      CPASSERT(ispin >= 1 .AND. ispin <= SIZE(matrix_rs, 1))
      CPASSERT(ASSOCIATED(kpoints))
      CPASSERT(ASSOCIATED(kp))

      CALL dbcsr_init_p(matrix_re)
      CALL dbcsr_init_p(matrix_im)
      CALL copy_fm_to_dbcsr_bc(cache_re, matrix_re)
      CALL copy_fm_to_dbcsr_bc(cache_im, matrix_im)
   END SUBROUTINE kpoint_operator_get_local

! **************************************************************************************************
!> \brief Fill caller-owned matrices for one operator row, including independent history snapshots.
!>        Preparation is collective. One-rank groups evaluate k points
!>        concurrently, with no MPI or matrix allocation/refcounts inside the worker loop.
!> \param kpoints host environment
!> \param ao_ao_fm global-layout template
!> \param rsmat current real-space images
!> \param slot operator slot
!> \param cmat allocated group-local outputs, one per locally owned k point
!> \param matrix_row source row; S/T use kpoint_spin_free
! **************************************************************************************************
   SUBROUTINE kpoint_operator_get_all(kpoints, ao_ao_fm, rsmat, slot, cmat, matrix_row)
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(cp_fm_type), INTENT(IN)                       :: ao_ao_fm
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rsmat
      INTEGER, INTENT(IN)                                :: slot
      TYPE(cp_cfm_type), DIMENSION(:), INTENT(INOUT)     :: cmat
      INTEGER, INTENT(IN)                                :: matrix_row

      CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_operator_get_all'

      INTEGER                                            :: handle, ik, ikp, kplocal, kpmax, &
                                                            nworkers, width
      INTEGER(KIND=int_8)                                :: work_size
      LOGICAL                                            :: prepared
      TYPE(kpoint_operator_context_type)                 :: ctx

      CALL timeset(routineN, handle)
      kplocal = SIZE(kpoints%kp_env)
      CPASSERT(SIZE(cmat) == kplocal .AND. kplocal > 0)
      CALL kpoint_operator_context_create(ctx, kpoints, ao_ao_fm, 1)
      CALL kpoint_operator_prepare(ctx, rsmat, slot, prepared=prepared, matrix_row=matrix_row)
      IF (prepared .AND. kpoints%para_env_kp%num_pe == 1) THEN
         nworkers = 1
!$       nworkers = MIN(omp_get_max_threads(), kplocal)
         CALL kpoint_operator_batch_size(ctx, nworkers, 1, width, work_size)
!$OMP PARALLEL DEFAULT(NONE) NUM_THREADS(nworkers) SHARED(ctx,kpoints,cmat,slot,kplocal) &
!$OMP PRIVATE(ikp,ik)
!$OMP DO SCHEDULE(DYNAMIC, 1)
         DO ikp = 1, kplocal
            ik = kpoints%kp_range(1) + ikp - 1
            CALL kpoint_operator_evaluate(ctx, 1, slot, ik, cmat(ikp))
         END DO
!$OMP END DO
!$OMP END PARALLEL
      ELSE
         kpmax = MAXVAL(kpoints%kp_dist(2, :) - kpoints%kp_dist(1, :) + 1)
         DO ikp = 1, kpmax
            CALL kpoint_operator_get(ctx, ikp, 1, rsmat, cmat=cmat(MIN(ikp, kplocal)), &
                                     matrix_row=matrix_row, slot=slot)
         END DO
      END IF
      CALL kpoint_operator_context_release(ctx)
      CALL timestop(handle)
   END SUBROUTINE kpoint_operator_get_all

! **************************************************************************************************
!> \brief Borrow the group-local OT cache for one operator slot.
!> \param kp owning k-point environment
!> \param slot KS, overlap or kinetic operator slot
!> \param ispin spin index, used only for the KS slot
!> \return non-owning real/imaginary FM array
! **************************************************************************************************
   FUNCTION kpoint_operator_cache(kp, slot, ispin) RESULT(cache)
      TYPE(kpoint_env_type), INTENT(IN)                  :: kp
      INTEGER, INTENT(IN)                                :: slot, ispin
      TYPE(cp_fm_type), POINTER                          :: cache(:)

      SELECT CASE (slot)
      CASE (kpoint_slot_ks)
         cache => kp%ot_hmat(:, ispin)
      CASE (kpoint_slot_s)
         cache => kp%ot_smat
      CASE (kpoint_slot_t)
         cache => kp%ot_tmat
      END SELECT
   END FUNCTION kpoint_operator_cache

! **************************************************************************************************
!> \brief Fill existing split-FM caches through the common prepared batch evaluator.
!>        Each worker owns a k tile and spin. Direct/mirrored destinations
!>        share contractions; packed input is immutable and reused across k tiles.
!> \param ctx prepared group-local context and owned output caches
!> \param slot operator slot
!> \param nrows number of captured spin rows
! **************************************************************************************************
   SUBROUTINE kpoint_operator_store_prepared(ctx, slot, nrows)
      TYPE(kpoint_operator_context_type), INTENT(INOUT)  :: ctx
      INTEGER, INTENT(IN)                                :: slot, nrows

      CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_operator_store_prepared'

      INTEGER                                            :: first, handle, ic, ikp, ispin, last, nc, &
                                                            nk, nworkers, width
      INTEGER(KIND=int_8)                                :: work_size
      REAL(KIND=dp), ALLOCATABLE                         :: work(:)
      TYPE(cp_fm_p_type), ALLOCATABLE                    :: output(:, :, :)
      TYPE(cp_fm_type), POINTER                          :: cache(:)

      CALL timeset(routineN, handle)
      nk = SIZE(ctx%kpoints%kp_env)
      IF (ctx%streaming(slot)) THEN
         DO first = 1, MAXVAL(ctx%kpoints%kp_dist(2, :) - ctx%kpoints%kp_dist(1, :) + 1)
            DO ispin = 1, nrows
               IF (first <= nk) THEN
                  cache => kpoint_operator_cache(ctx%kpoints%kp_env(first)%kpoint_env, slot, ispin)
                  IF (ctx%kpoints%use_real_wfn) THEN
                     CALL kpoint_operator_evaluate(ctx, ispin, slot, first, fm_re=cache(1), round=.TRUE.)
                  ELSE
                     CALL kpoint_operator_evaluate(ctx, ispin, slot, first, &
                                                   fm_re=cache(1), fm_im=cache(2), round=.TRUE.)
                  END IF
               ELSE
                  CALL kpoint_operator_evaluate(ctx, ispin, slot, first, fm_re=ctx%fmlocal, round=.TRUE.)
               END IF
            END DO
         END DO
         CALL timestop(handle)
         RETURN
      END IF
      IF (nk == 0) THEN
         CALL timestop(handle)
         RETURN
      END IF
      nc = MERGE(1, 2, ctx%kpoints%use_real_wfn)
      ALLOCATE (output(nk, nrows, nc))
      DO ikp = 1, nk
         DO ispin = 1, nrows
            cache => kpoint_operator_cache(ctx%kpoints%kp_env(ikp)%kpoint_env, slot, ispin)
            DO ic = 1, nc
               output(ikp, ispin, ic)%matrix => cache(ic)
            END DO
         END DO
      END DO
      nworkers = 1
!$    nworkers = omp_get_max_threads()
      CALL rskp_batch_size(ctx%rs_route, nk, nc, nworkers, 0, rs_batch_bytes, width, work_size)
!$OMP PARALLEL DEFAULT(NONE) NUM_THREADS(nworkers) &
!$OMP SHARED(ctx,slot,output,nk,nrows,width,work_size) PRIVATE(work,first,last,ispin)
      ALLOCATE (work(work_size))
!$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC, 1)
      DO first = 1, nk, width
         DO ispin = 1, nrows
            last = MIN(nk, first + width - 1)
            CALL cell_to_k_backend_evaluate_batch(ctx%rs(ispin, slot), ctx%rs_route, &
                                                  ctx%kpoints%xkp(:, ctx%kpoints%kp_range(1) + first - 1: &
                                                                  ctx%kpoints%kp_range(1) + last - 1), &
                                                  first, work, output=output(first:last, ispin, :))
         END DO
      END DO
!$OMP END DO
      DEALLOCATE (work)
!$OMP END PARALLEL
      CALL timestop(handle)
   END SUBROUTINE kpoint_operator_store_prepared

! **************************************************************************************************
!> \brief Size local k tiles, budgeting operator outputs and block scratch over all workers.
!>        A block exceeding the budget is processed by one worker with one k point.
!>        OT writes existing caches and distributes block/spin jobs; dense callers distribute k tiles.
!> \param ctx prepared group-local context
!> \param nworkers requested workers, reduced to fit the minimum working set
!> \param nmatrices simultaneously retained complex operator tiles; zero writes existing split-FM caches
!> \param width number of k points per worker
!> \param work_size real block scratch elements per worker
! **************************************************************************************************
   SUBROUTINE kpoint_operator_batch_size(ctx, nworkers, nmatrices, width, work_size)
      TYPE(kpoint_operator_context_type), INTENT(IN)     :: ctx
      INTEGER, INTENT(INOUT)                             :: nworkers
      INTEGER, INTENT(IN)                                :: nmatrices
      INTEGER, INTENT(OUT)                               :: width
      INTEGER(KIND=int_8), INTENT(OUT)                   :: work_size

      INTEGER                                            :: nc, nk

      CPASSERT(nmatrices >= 0 .AND. nworkers > 0)
      nk = SIZE(ctx%kpoints%kp_env)
      IF (ctx%rs_route%world_target) nk = ctx%kpoints%nkp
      nc = MERGE(1, 2, ctx%kpoints%use_real_wfn)
      CALL rskp_batch_size(ctx%rs_route, nk, nc, nworkers, nmatrices, rs_batch_bytes, width, work_size)
   END SUBROUTINE kpoint_operator_batch_size

! **************************************************************************************************
!> \brief Deliver a complex k tile using real block/image contractions or the prepared FFT cache.
!>        The prepared route and snapshots are immutable; workers share no output or scratch.
!> \param ctx prepared group-local context
!> \param ispin operator spin
!> \param slot prepared operator slot
!> \param first_k first global k index
!> \param batch complex output (local row, local column, k point), including padding
!> \param work private real scratch sized by kpoint_operator_batch_size
! **************************************************************************************************
   SUBROUTINE kpoint_operator_evaluate_batch(ctx, ispin, slot, first_k, batch, work)
      TYPE(kpoint_operator_context_type), INTENT(IN)     :: ctx
      INTEGER, INTENT(IN)                                :: ispin, slot, first_k
      COMPLEX(KIND=dp), CONTIGUOUS, INTENT(OUT)          :: batch(:, :, :)
      REAL(KIND=dp), CONTIGUOUS, INTENT(INOUT), TARGET   :: work(:)

      CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_operator_evaluate_batch'

      INTEGER                                            :: first, handle, spin, width

      CALL timeset(routineN, handle)
      CPASSERT(SIZE(batch, 1) >= ctx%rs_route%nrow)
      CPASSERT(SIZE(batch, 2) >= ctx%rs_route%ncol)
      width = SIZE(batch, 3)
      IF (ctx%rs_route%world_target) THEN
         CPASSERT(first_k >= 1 .AND. first_k + width - 1 <= ctx%kpoints%nkp)
      ELSE
         CPASSERT(first_k >= ctx%kpoints%kp_range(1))
         CPASSERT(first_k + width - 1 <= ctx%kpoints%kp_range(2))
      END IF
      spin = ispin
      IF (.NOT. ASSOCIATED(ctx%rs(spin, slot)%source)) spin = 1
      CPASSERT(ASSOCIATED(ctx%rs(spin, slot)%source))
      first = first_k
      IF (.NOT. ctx%rs_route%world_target) first = first_k - ctx%kpoints%kp_range(1) + 1
      CALL cell_to_k_backend_evaluate_batch(ctx%rs(spin, slot), ctx%rs_route, &
                                            ctx%kpoints%xkp(:, first_k:first_k + width - 1), &
                                            first, work, batch)
      CALL timestop(handle)
   END SUBROUTINE kpoint_operator_evaluate_batch

! **************************************************************************************************
!> \brief Prepare owned split-FM operators for OT without entering an SCF eigensolver.
!>        H (or the initialization density) retains all source spins and is refreshed on each call.
!>        S is stored once per context; the caller controls when to prepare T.
!>        Outputs survive the context. Only routing and scratch persist between calls.
!> \param ctx context with fixed geometry, AO layout and spin capacity
!> \param matrix_ks Hamiltonian or initialization density
!> \param matrix_s spin-free overlap
!> \param matrix_t optional spin-free kinetic operator; absent leaves the existing kinetic cache alone
! **************************************************************************************************
   SUBROUTINE kpoint_operator_store(ctx, matrix_ks, matrix_s, matrix_t)
      TYPE(kpoint_operator_context_type), INTENT(INOUT)  :: ctx
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks, matrix_s
      TYPE(dbcsr_p_type), DIMENSION(:, :), OPTIONAL, &
         POINTER                                         :: matrix_t

      CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_operator_store'
      INTEGER                                            :: handle, ic, ikp, ispin, ncomp, nspin
      LOGICAL                                            :: kinetic, reuse
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(kpoint_type), POINTER                         :: kpoints

      CALL timeset(routineN, handle)
      kpoints => ctx%kpoints
      nspin = SIZE(matrix_ks, 1)
      CPASSERT(nspin == 1 .OR. nspin == ctx%nspin)
      ncomp = MERGE(1, 2, kpoints%use_real_wfn)
      kinetic = PRESENT(matrix_t)
      IF (kinetic) kinetic = ASSOCIATED(matrix_t)
      CALL kpoint_operator_ensure_group(ctx)
      matrix_struct => ctx%fmlocal%matrix_struct
      DO ikp = 1, SIZE(kpoints%kp_env)
         kp => kpoints%kp_env(ikp)%kpoint_env
         reuse = .FALSE.
         IF (ASSOCIATED(kp%ot_hmat)) THEN
            IF (ALL(SHAPE(kp%ot_hmat) == [ncomp, nspin])) THEN
               reuse = cp_fm_struct_equivalent(kp%ot_hmat(1, 1)%matrix_struct, matrix_struct)
            END IF
         END IF
         IF (.NOT. reuse) THEN
            CALL cp_fm_release(kp%ot_hmat)
            ALLOCATE (kp%ot_hmat(ncomp, nspin))
            DO ispin = 1, nspin
               DO ic = 1, ncomp
                  CALL cp_fm_create(kp%ot_hmat(ic, ispin), matrix_struct)
               END DO
            END DO
         END IF
         CALL allocate_spin_free(kp%ot_smat)
         IF (kinetic) CALL allocate_spin_free(kp%ot_tmat)
      END DO
      CALL store_slot(matrix_ks, kpoint_slot_ks, nspin)
      IF (.NOT. ctx%overlap_stored) THEN
         CALL store_slot(matrix_s, kpoint_slot_s, 1)
         ctx%overlap_stored = .TRUE.
      END IF
      IF (kinetic) CALL store_slot(matrix_t, kpoint_slot_t, 1)
      CALL timestop(handle)

   CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param cache ...
! **************************************************************************************************
      SUBROUTINE allocate_spin_free(cache)
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: cache

      INTEGER                                            :: ic
      LOGICAL                                            :: reuse

         reuse = .FALSE.
         IF (ASSOCIATED(cache)) THEN
            IF (SIZE(cache) == ncomp) reuse = cp_fm_struct_equivalent(cache(1)%matrix_struct, matrix_struct)
         END IF
         IF (reuse) RETURN
         CALL cp_fm_release(cache)
         ALLOCATE (cache(ncomp))
         DO ic = 1, ncomp
            CALL cp_fm_create(cache(ic), matrix_struct)
         END DO
      END SUBROUTINE allocate_spin_free

! **************************************************************************************************
!> \brief ...
!> \param rsmat ...
!> \param slot ...
!> \param nrows ...
! **************************************************************************************************
      SUBROUTINE store_slot(rsmat, slot, nrows)
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rsmat
      INTEGER, INTENT(IN)                                :: slot, nrows

      INTEGER                                            :: ispin

         IF (slot == kpoint_slot_ks .AND. nrows == ctx%nspin) THEN
            CALL kpoint_operator_prepare(ctx, rsmat, slot)
         ELSE
            CALL kpoint_operator_prepare(ctx, rsmat, slot, matrix_row=kpoint_spin_free)
         END IF
         CALL kpoint_operator_store_prepared(ctx, slot, nrows)
         ! The split-FM caches now own the results. Do not retain source block
         ! pointers or duplicate snapshots across Hamiltonian rebuilds.
         IF (ALLOCATED(ctx%rs)) THEN
            DO ispin = 1, ctx%nspin
               CALL rskp_values_release(ctx%rs(ispin, slot))
            END DO
         END IF
         CALL rskp_transform_release(ctx%stream(slot))
         ctx%streaming(slot) = .FALSE.
      END SUBROUTINE store_slot

   END SUBROUTINE kpoint_operator_store

END MODULE qs_kpoint_operators
