!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Prepared real-space/reciprocal-space transforms and distributed AO adapters.
!>        Routes fix geometry and ownership. Captured source rows own numerical
!>        values; evaluation reads them without communication or SCF state.
! **************************************************************************************************
MODULE kpoint_transform
   USE cp_cfm_types,                    ONLY: cp_cfm_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_get_info,&
                                              dbcsr_get_readonly_block_p,&
                                              dbcsr_p_type,&
                                              dbcsr_type
   USE cp_fm_struct,                    ONLY: cp_fm_struct_equivalent,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_retain,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_p_type,&
                                              cp_fm_type
   USE input_constants,                 ONLY: rskp_backend_fft,&
                                              rskp_backend_gemm
   USE kinds,                           ONLY: dp,&
                                              int_8
   USE kpoint_lattice_fft,              ONLY: cell_to_k_grid_fft,&
                                              cell_to_k_grid_prepare,&
                                              cell_to_k_grid_release,&
                                              cell_to_k_grid_work_type
   USE local_gemm_api,                  ONLY: local_gemm_ctxt_type
   USE mathconstants,                   ONLY: twopi,&
                                              z_one,&
                                              z_zero
   USE message_passing,                 ONLY: mp_para_env_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              get_neighbor_list_set_p,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
   USE util,                            ONLY: sort
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'kpoint_transform'
   INTEGER(KIND=int_8), PARAMETER :: fft_batch_bytes = 32_int_8*1024_int_8**2
   INTEGER(KIND=int_8), PARAMETER, PRIVATE :: rskp_route_bytes = 64_int_8*1024_int_8**2, &
                                              rskp_fft_bytes = 512_int_8*1024_int_8**2, &
                                              rskp_gemm_bytes = 1000_int_8*1024_int_8**2, &
                                              dbcsr_batch_bytes = 1000_int_8*1024_int_8**2
   PUBLIC :: kp_transform_plan_type, kp_transform_plan_create, rskp_route_type, rskp_values_type
   PUBLIC :: rskp_route_create, rskp_route_release, rskp_capture, rskp_compatible, rskp_values_release
   PUBLIC :: rskp_backend_gemm, rskp_backend_fft, cell_to_k_backend_prepare
   PUBLIC :: rskp_batch_size, cell_to_k_backend_evaluate_batch
   PUBLIC :: rskp_dbcsr_transform_batch
   PUBLIC :: rskp_prepared_transform_type, rskp_transform_prepare
   PUBLIC :: rskp_transform_evaluate, rskp_transform_release
   PUBLIC :: cell_to_k_sparse_fft, k_grid_to_cell_contract, rskp_phase_factors, rskp_work_matrix

   TYPE, PUBLIC :: kp_transform_plan_type
      INTEGER :: nentry = 0, ngroup = 0
      LOGICAL :: symmetric = .FALSE.
      INTEGER, ALLOCATABLE, DIMENSION(:) :: col, col_offset, group_start, image, row, row_offset
      INTEGER, ALLOCATABLE, DIMENSION(:, :) :: cell
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: symmetry_sign
   END TYPE kp_transform_plan_type

   ! Numerical plans borrow no FM, DBCSR or SCF objects.
   TYPE rskp_backend_plan_type
      INTEGER :: backend = rskp_backend_gemm
      INTEGER, ALLOCATABLE :: packed_first(:), fold_map(:, :)
      INTEGER, ALLOCATABLE :: phase_index(:), phase_cells(:, :), phase_count(:, :), phase_first(:)
      INTEGER(KIND=int_8) :: max_per_k = 0_int_8
      REAL(KIND=dp), ALLOCATABLE :: phases(:, :, :)
      TYPE(cell_to_k_grid_work_type) :: fft_work
   END TYPE rskp_backend_plan_type

   ! Stable canonical AO route and destination mapping, shared by both backends.
   TYPE rskp_route_type
      TYPE(cp_fm_struct_type), POINTER :: matrix_struct => NULL()
      TYPE(mp_para_env_type), POINTER :: world => NULL(), group => NULL(), inter => NULL()
      TYPE(kp_transform_plan_type) :: traversal = kp_transform_plan_type()
      TYPE(neighbor_list_set_p_type), POINTER :: source_sab(:) => NULL()
      TYPE(rskp_backend_plan_type) :: transform
      INTEGER :: nrow = 0, ncol = 0, nimage = 0, nvalue = 0
      LOGICAL :: world_target = .FALSE.
      INTEGER(KIND=int_8) :: storage = 0_int_8
      INTEGER, ALLOCATABLE :: entry(:, :), block_first(:), block_order(:), element_first(:)
      INTEGER(KIND=int_8), ALLOCATABLE :: element_offset(:, :)
      INTEGER, ALLOCATABLE :: row_size(:), col_size(:), first(:)
      INTEGER, ALLOCATABLE :: send_entry(:), send_first(:), scount(:, :), rcount(:, :), &
                              sdisp(:, :), rdisp(:, :), send_element(:)
      REAL(KIND=dp), ALLOCATABLE :: send_values(:), recv_values(:)
   END TYPE rskp_route_type

   ! One source row: backend payload plus the source identity used by adapters.
   TYPE rskp_values_type
      REAL(KIND=dp), ALLOCATABLE :: packed(:)
      LOGICAL, ALLOCATABLE :: active(:, :)
      COMPLEX(KIND=dp), ALLOCATABLE :: fft_values(:, :, :)
      TYPE(dbcsr_p_type), POINTER :: source(:, :) => NULL()
      INTEGER :: matrix_row = 0, mirror_sign = 1
      COMPLEX(KIND=dp) :: scale = z_one
   END TYPE rskp_values_type

   ! Slot-free owner for captured source rows sharing one target-local R-to-K map.
   TYPE, PUBLIC :: rskp_prepared_transform_type
      PRIVATE
      TYPE(rskp_route_type) :: route
      TYPE(rskp_values_type), ALLOCATABLE :: values(:)
      INTEGER, ALLOCATABLE :: kp_dist(:, :)
      REAL(KIND=dp), ALLOCATABLE :: work(:)
   END TYPE rskp_prepared_transform_type

CONTAINS

! **************************************************************************************************
!> \brief Prepare source rows for a generic group- or world-distributed target.
!>        Own route, snapshots, singleton scratch and any FFT cache without SCF slots.
!> \param prepared transform handle to initialize
!> \param source real-space source rows and images
!> \param source_rows source rows sharing the route/backend plan
!> \param source_sab neighbor-list geometry of source
!> \param cell_to_index translation-to-image map
!> \param target_struct output matrix layout
!> \param world source-matrix communicator
!> \param group optional target-group communicator; omit for a world target
!> \param inter optional communicator joining equivalent group ranks; omit for a world target
!> \param backend selected GEMM or FFT backend
!> \param xkp full global k-point coordinates
!> \param nkp_grid regular-grid dimensions
!> \param kp_dist optional global k ranges by inter-group rank; omit for a world target
!> \param mirror_sign transpose parity (+1 Hermitian, -1 anti-Hermitian)
!> \param scale complex prefactor applied to direct and mirrored contributions
! **************************************************************************************************
   SUBROUTINE rskp_transform_prepare(prepared, source, source_rows, source_sab, cell_to_index, target_struct, &
                                     world, group, inter, backend, xkp, nkp_grid, kp_dist, mirror_sign, scale)

      TYPE(rskp_prepared_transform_type), INTENT(INOUT)  :: prepared
      TYPE(dbcsr_p_type), POINTER                        :: source(:, :)
      INTEGER, INTENT(IN)                                :: source_rows(:)
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: source_sab
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      TYPE(cp_fm_struct_type), POINTER                   :: target_struct
      TYPE(mp_para_env_type), POINTER                    :: world
      TYPE(mp_para_env_type), OPTIONAL, POINTER          :: group, inter
      INTEGER, INTENT(IN)                                :: backend
      REAL(KIND=dp), INTENT(IN)                          :: xkp(:, :)
      INTEGER, INTENT(IN)                                :: nkp_grid(3)
      INTEGER, INTENT(IN), OPTIONAL                      :: kp_dist(:, :), mirror_sign
      COMPLEX(KIND=dp), INTENT(IN), OPTIONAL             :: scale

      CHARACTER(LEN=*), PARAMETER :: routineN = 'rskp_transform_prepare'

      INTEGER                                            :: first, handle, irow, nk_local
      INTEGER(KIND=int_8)                                :: max_backend_bytes, phase_bytes, &
                                                            retained_bytes, workspace_bytes
      INTEGER, ALLOCATABLE                               :: eval_kp_dist(:, :)
      LOGICAL                                            :: ready
      REAL(KIND=dp), ALLOCATABLE                         :: phase_xkp(:, :)
      TYPE(mp_para_env_type), POINTER                    :: eval_group, eval_inter

      CALL timeset(routineN, handle)
      CALL rskp_transform_release(prepared)
      CPASSERT(SIZE(source_rows) > 0)
      CPASSERT(ALL(source_rows >= 1) .AND. ALL(source_rows <= SIZE(source, 1)))
      CPASSERT(SIZE(xkp, 1) == 3 .AND. SIZE(xkp, 2) > 0)
      CPASSERT(ASSOCIATED(target_struct) .AND. ASSOCIATED(world))
      IF (PRESENT(group) .EQV. PRESENT(inter)) THEN
         IF (PRESENT(group)) THEN
            CPASSERT(PRESENT(kp_dist))
            CPASSERT(ASSOCIATED(group) .AND. ASSOCIATED(inter))
            eval_group => group
            eval_inter => inter
         ELSE
            CPASSERT(.NOT. PRESENT(kp_dist))
            CPASSERT(ASSOCIATED(target_struct%para_env, world))
            eval_group => world
            eval_inter => world
         END IF
      ELSE
         CPABORT("R-to-K group and inter communicators must be supplied together")
      END IF
      IF (PRESENT(kp_dist)) THEN
         CPASSERT(SIZE(kp_dist, 1) == 2)
         ALLOCATE (eval_kp_dist, source=kp_dist)
      ELSE
         ALLOCATE (eval_kp_dist(2, 1))
         eval_kp_dist(:, 1) = [1, SIZE(xkp, 2)]
      END IF
      max_backend_bytes = rskp_gemm_bytes
      IF (backend == rskp_backend_fft) max_backend_bytes = rskp_fft_bytes
      ALLOCATE (prepared%kp_dist, source=eval_kp_dist)
      ALLOCATE (prepared%values(SIZE(source_rows)))
      CALL rskp_route_create(source, source_sab, cell_to_index, target_struct, world, eval_group, eval_inter, backend, &
                             rskp_route_bytes, prepared%route, ready)
      IF (.NOT. ready) THEN
         CPABORT("R-to-K route is unsupported, incompatible with the target, or over budget")
      END IF

      phase_bytes = max_backend_bytes
      IF (backend == rskp_backend_fft .OR. prepared%route%world_target) THEN
         workspace_bytes = max_backend_bytes
         IF (backend == rskp_backend_gemm) workspace_bytes = phase_bytes
         CALL cell_to_k_backend_prepare(prepared%route, xkp, nkp_grid, workspace_bytes)
      ELSE
         first = eval_kp_dist(1, eval_inter%mepos + 1)
         ALLOCATE (phase_xkp, source=xkp(:, first:eval_kp_dist(2, eval_inter%mepos + 1)))
         CALL cell_to_k_backend_prepare(prepared%route, phase_xkp, nkp_grid, phase_bytes)
      END IF
      retained_bytes = prepared%route%storage
      DO irow = 1, SIZE(source_rows)
         CALL rskp_capture(source, source_rows(irow), prepared%route, prepared%values(irow), mirror_sign, scale)
         retained_bytes = retained_bytes + 8_int_8*SIZE(prepared%values(irow)%packed, KIND=int_8)
         IF (ALLOCATED(prepared%values(irow)%active)) THEN
            retained_bytes = retained_bytes + INT(STORAGE_SIZE(.FALSE.)/8, int_8)* &
                             SIZE(prepared%values(irow)%active, KIND=int_8)
         END IF
      END DO
      IF (retained_bytes > rskp_route_bytes) CPABORT("Prepared R-to-K source exceeds the memory budget")
      IF (backend == rskp_backend_fft) THEN
         nk_local = SIZE(xkp, 2)
         IF (.NOT. prepared%route%world_target) THEN
            nk_local = eval_kp_dist(2, eval_inter%mepos + 1) - &
                       eval_kp_dist(1, eval_inter%mepos + 1) + 1
         END IF
         retained_bytes = retained_bytes + 16_int_8*SIZE(source_rows, KIND=int_8)*nk_local* &
                          target_struct%local_leading_dimension*MAX(1, prepared%route%ncol)
         workspace_bytes = max_backend_bytes - retained_bytes
         IF (workspace_bytes <= 0_int_8) CPABORT("Prepared lattice FFT exceeds the memory budget")
         DO irow = 1, SIZE(source_rows)
            ALLOCATE (prepared%values(irow)%fft_values(target_struct%local_leading_dimension, &
                                                       MAX(1, prepared%route%ncol), nk_local))
            CALL cell_to_k_sparse_fft(prepared%values(irow), prepared%route, xkp, nkp_grid, &
                                      eval_kp_dist, eval_inter, &
                                      prepared%values(irow)%fft_values, workspace_bytes)
            DEALLOCATE (prepared%values(irow)%packed)
         END DO
         ALLOCATE (prepared%work(0))
      ELSE
         ALLOCATE (prepared%work(2_int_8*prepared%route%transform%max_per_k))
      END IF
      CALL timestop(handle)
   END SUBROUTINE rskp_transform_prepare

! **************************************************************************************************
!> \brief Evaluate one k point directly into an FM, split-FM or CFM target.
!>        The handle owns mutable scratch; concurrent evaluations need separate handles.
!> \param prepared initialized transform handle
!> \param source_index one-based source row index supplied at preparation
!> \param k global k-point index
!> \param xkp full global k-point coordinates
!> \param fm_re real output matrix for split-FM storage
!> \param fm_im optional imaginary output matrix for split-FM storage
!> \param cfm optional complex full-matrix output
! **************************************************************************************************
   SUBROUTINE rskp_transform_evaluate(prepared, source_index, k, xkp, fm_re, fm_im, cfm)
      TYPE(rskp_prepared_transform_type), INTENT(INOUT)  :: prepared
      INTEGER, INTENT(IN)                                :: source_index, k
      REAL(KIND=dp), INTENT(IN)                          :: xkp(:, :)
      TYPE(cp_fm_type), INTENT(INOUT), OPTIONAL, TARGET  :: fm_re, fm_im
      TYPE(cp_cfm_type), INTENT(INOUT), OPTIONAL, TARGET :: cfm

      COMPLEX(KIND=dp), CONTIGUOUS, POINTER              :: batch(:, :, :)
      INTEGER                                            :: first_local, ncomp
      TYPE(cp_fm_p_type)                                 :: output(1, 2)

      CPASSERT(ALLOCATED(prepared%values) .AND. ALLOCATED(prepared%work))
      CPASSERT(source_index >= 1 .AND. source_index <= SIZE(prepared%values))
      CPASSERT(k >= 1 .AND. k <= SIZE(xkp, 2))
      first_local = k
      IF (.NOT. prepared%route%world_target) THEN
         CPASSERT(k >= prepared%kp_dist(1, prepared%route%inter%mepos + 1))
         CPASSERT(k <= prepared%kp_dist(2, prepared%route%inter%mepos + 1))
         first_local = k - prepared%kp_dist(1, prepared%route%inter%mepos + 1) + 1
      END IF
      CPASSERT(PRESENT(fm_re) .NEQV. PRESENT(cfm))
      CPASSERT(.NOT. PRESENT(fm_im) .OR. PRESENT(fm_re))
      IF (PRESENT(cfm)) THEN
         CPASSERT(cp_fm_struct_equivalent(cfm%matrix_struct, prepared%route%matrix_struct))
         batch(1:SIZE(cfm%local_data, 1), 1:SIZE(cfm%local_data, 2), 1:1) => cfm%local_data
         CALL cell_to_k_backend_evaluate_batch(prepared%values(source_index), prepared%route, &
                                               xkp(:, k:k), first_local, prepared%work, batch=batch)
      ELSE
         output(1, 1)%matrix => fm_re
         ncomp = 1
         IF (PRESENT(fm_im)) THEN
            output(1, 2)%matrix => fm_im
            ncomp = 2
         END IF
         CALL cell_to_k_backend_evaluate_batch(prepared%values(source_index), prepared%route, &
                                               xkp(:, k:k), first_local, prepared%work, output=output(:, 1:ncomp))
      END IF
   END SUBROUTINE rskp_transform_evaluate

! **************************************************************************************************
!> \brief Release the route, snapshot and cached backend state of a prepared transform.
!> \param prepared transform handle to release
! **************************************************************************************************
   SUBROUTINE rskp_transform_release(prepared)
      TYPE(rskp_prepared_transform_type), INTENT(INOUT)  :: prepared

      INTEGER                                            :: irow

      IF (ALLOCATED(prepared%values)) THEN
         DO irow = 1, SIZE(prepared%values)
            CALL rskp_values_release(prepared%values(irow))
         END DO
         DEALLOCATE (prepared%values)
      END IF
      CALL rskp_route_release(prepared%route)
      IF (ALLOCATED(prepared%kp_dist)) DEALLOCATE (prepared%kp_dist)
      IF (ALLOCATED(prepared%work)) DEALLOCATE (prepared%work)
   END SUBROUTINE rskp_transform_release

! **************************************************************************************************
!> \brief Transform all real-space spin rows to a full k-grid in distributed DBCSR matrices.
!>        The real and imaginary output blocks retain their caller-provided sparsity and
!>        distribution. GEMM contracts all images of an AO block together; FFT transforms
!>        element tiles while reusing one prepared lattice geometry.
!> \param source real-space image matrices, by spin and image
!> \param real_k symmetric real output matrices, by spin and k point
!> \param imag_k antisymmetric imaginary output matrices, by spin and k point
!> \param xkp full reciprocal grid
!> \param cell_to_index translation-to-image map
!> \param sab_nl neighbor-list geometry of source
!> \param nkp_grid dimensions of the regular reciprocal grid
!> \param backend selected bulk R-to-K backend
! **************************************************************************************************
   SUBROUTINE rskp_dbcsr_transform_batch(source, real_k, imag_k, xkp, cell_to_index, sab_nl, nkp_grid, backend)

      TYPE(dbcsr_p_type), POINTER                        :: source(:, :), real_k(:, :), imag_k(:, :)
      REAL(KIND=dp), INTENT(IN)                          :: xkp(:, :)
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      INTEGER, INTENT(IN)                                :: nkp_grid(3), backend

      CHARACTER(LEN=*), PARAMETER :: routineN = 'rskp_dbcsr_transform_batch'

      COMPLEX(KIND=dp), ALLOCATABLE                      :: fft_output(:, :, :)
      INTEGER :: col, first, first_k, g, handle, i, ik, ispin, last, max_images, max_nelem, ncol, &
         nelem, nentry, nfold, ngroup, nimages, nk, nphase, nrow, nspin, odd, offset, panels(4), &
         pos, row_first, row_last, tile_first, tile_size, tile_width, width
      INTEGER(KIND=int_8)                                :: available_bytes, fft_tile_budget, &
                                                            nfactor, panel_bytes, per_k_bytes, &
                                                            per_tile_bytes
      INTEGER, ALLOCATABLE                               :: element_first(:), group_first(:), &
                                                            order(:), signed_cells(:, :)
      INTEGER, DIMENSION(:), POINTER                     :: col_size, row_size
      LOGICAL                                            :: active(4), found, found_im, found_re
      REAL(KIND=dp), ALLOCATABLE, TARGET                 :: packed(:), work(:)
      REAL(KIND=dp), CONTIGUOUS, POINTER                 :: factors(:, :), fft_input(:, :, :), &
                                                            image_re(:, :), RESULT(:, :)
      REAL(KIND=dp), POINTER                             :: iblock(:, :), rblock(:, :), &
                                                            result_block(:, :), rsblock(:, :)
      TYPE(cell_to_k_grid_work_type)                     :: fft_work
      TYPE(kp_transform_plan_type)                       :: plan
      TYPE(rskp_backend_plan_type)                       :: transform

      CALL timeset(routineN, handle)
      nspin = SIZE(source, 1)
      nk = SIZE(xkp, 2)
      CPASSERT(nk > 0)
      IF (backend /= rskp_backend_gemm .AND. backend /= rskp_backend_fft) THEN
         CPABORT("Unknown R-to-K backend")
      END IF
      CPASSERT(nspin > 0 .AND. SIZE(real_k, 1) == nspin .AND. SIZE(imag_k, 1) == nspin)
      CPASSERT(SIZE(real_k, 2) == nk .AND. SIZE(imag_k, 2) == nk)
      CPASSERT(SIZE(source, 2) > 0 .AND. ASSOCIATED(sab_nl) .AND. ASSOCIATED(cell_to_index))
      DO ik = 1, nk
         DO ispin = 1, nspin
            CALL dbcsr_set(real_k(ispin, ik)%matrix, 0.0_dp)
            CALL dbcsr_set(imag_k(ispin, ik)%matrix, 0.0_dp)
         END DO
      END DO

      CALL kp_transform_plan_create(plan, sab_nl, cell_to_index, SIZE(source, 2), group_entries=.TRUE.)
      nentry = plan%nentry
      IF (nentry == 0) THEN
         CALL timestop(handle)
         RETURN
      END IF
      CALL dbcsr_get_info(source(1, 1)%matrix, row_blk_size=row_size, col_blk_size=col_size)
      ALLOCATE (group_first(nentry + 1), element_first(nentry + 1), order(nentry), signed_cells(3, nentry))
      DO i = 1, nentry
         order(i) = i
         signed_cells(:, i) = NINT(plan%symmetry_sign(i))*plan%cell(:, i)
      END DO
      max_images = 0
      max_nelem = 0
      ngroup = 0
      element_first(1) = 1
      first = 1
      DO WHILE (first <= nentry)
         last = first
         DO WHILE (last < nentry)
            IF (plan%row(last + 1) /= plan%row(first) .OR. plan%col(last + 1) /= plan%col(first)) EXIT
            last = last + 1
         END DO
         ngroup = ngroup + 1
         group_first(ngroup) = first
         panel_bytes = INT(row_size(plan%row(first)), int_8)*col_size(plan%col(first))*nspin
         IF (panel_bytes + element_first(ngroup) > HUGE(nelem)) THEN
            CPABORT("R-to-K DBCSR panel offsets exceed the integer range")
         END IF
         element_first(ngroup + 1) = element_first(ngroup) + INT(panel_bytes)
         max_nelem = MAX(max_nelem, INT(panel_bytes))
         max_images = MAX(max_images, last - first + 1)
         first = last + 1
      END DO
      group_first(ngroup + 1) = nentry + 1
      panel_bytes = 8_int_8*max_nelem*max_images
      IF (panel_bytes >= dbcsr_batch_bytes) THEN
         CPABORT("R-to-K DBCSR input panel exceeds the memory budget")
      END IF
      available_bytes = dbcsr_batch_bytes - panel_bytes
      IF (backend == rskp_backend_gemm) THEN
         CALL cell_to_k_phase_prepare(transform, signed_cells, group_first(1:ngroup + 1), order, &
                                      element_first(1:ngroup + 1), xkp, available_bytes, 0_int_8)
         IF (ALLOCATED(transform%phases)) THEN
            available_bytes = available_bytes - 8_int_8*SIZE(transform%phases, KIND=int_8)
         END IF
         per_k_bytes = 16_int_8*transform%max_per_k
         tile_width = INT(MIN(INT(nk, int_8), available_bytes/MAX(1_int_8, per_k_bytes)))
         IF (tile_width < 1) CPABORT("R-to-K DBCSR GEMM workspace exceeds the memory budget")
         ALLOCATE (work(per_k_bytes/8*tile_width))
      ELSE
         available_bytes = available_bytes - fft_batch_bytes
         IF (available_bytes <= 0_int_8) CPABORT("R-to-K DBCSR FFT workspace exceeds the memory budget")
         fft_tile_budget = MIN(fft_batch_bytes, available_bytes)
         CALL cell_to_k_grid_prepare(fft_work, signed_cells, xkp, nkp_grid)
         IF (.NOT. fft_work%regular) THEN
            CPABORT("Forced lattice FFT needs a supported complete regular grid")
         END IF
      END IF
      ALLOCATE (packed(panel_bytes/8))

      DO g = 1, ngroup
         first = group_first(g)
         last = group_first(g + 1) - 1
         nimages = last - first + 1
         nrow = row_size(plan%row(first))
         ncol = col_size(plan%col(first))
         nelem = nrow*ncol
         IF (backend == rskp_backend_gemm) THEN
            panels = transform%phase_count([1, 1, 2, 3], g)
            nfold = SUM(panels)
         ELSE
            nfold = nimages
         END IF
         image_re(1:nelem*nspin, 1:nfold) => packed(1:nelem*nspin*nfold)
         image_re = 0.0_dp
         DO i = first, last
            pos = i - first + 1
            odd = 0
            IF (backend == rskp_backend_gemm) THEN
               pos = (transform%fold_map(1, i) - transform%packed_first(g))/(nelem*nspin) + 1
               IF (transform%fold_map(2, i) /= 0) THEN
                  odd = (ABS(transform%fold_map(2, i)) - transform%packed_first(g))/(nelem*nspin) + 1
               END IF
            END IF
            DO ispin = 1, nspin
               CALL dbcsr_get_readonly_block_p(source(ispin, plan%image(i))%matrix, plan%row(first), plan%col(first), &
                                               rsblock, found)
               IF (.NOT. found) CYCLE
               CPASSERT(SIZE(rsblock, 1) == nrow .AND. SIZE(rsblock, 2) == ncol)
               image_re((ispin - 1)*nelem + 1:ispin*nelem, pos) = &
                  image_re((ispin - 1)*nelem + 1:ispin*nelem, pos) + RESHAPE(rsblock, [nelem])
               IF (odd == 0) CYCLE
               image_re((ispin - 1)*nelem + 1:ispin*nelem, odd) = &
                  image_re((ispin - 1)*nelem + 1:ispin*nelem, odd) + &
                  REAL(SIGN(1, transform%fold_map(2, i)), dp)*RESHAPE(rsblock, [nelem])
            END DO
         END DO
         IF (.NOT. ANY(image_re /= 0.0_dp)) CYCLE

         IF (backend == rskp_backend_gemm) THEN
            pos = 1
            DO i = 1, 4
               active(i) = ANY(image_re(:, pos:pos + panels(i) - 1) /= 0.0_dp)
               pos = pos + panels(i)
            END DO
            nphase = MAXVAL(transform%phase_count(1:2, g))
            DO first_k = 1, nk, tile_width
               width = MIN(tile_width, nk - first_k + 1)
               nfactor = 2_int_8*nphase*width
               factors(1:nphase, 1:2*width) => work(1:nfactor)
               RESULT(1:nelem*nspin, 1:2*width) => work(nfactor + 1:nfactor + 2_int_8*nelem*nspin*width)
               CALL cell_to_k_block_gemm(packed(1:nelem*nspin*nfold), transform, g, &
                                         xkp(:, first_k:first_k + width - 1), first_k, factors, RESULT, active)
               DO ik = 1, width
                  DO ispin = 1, nspin
                     CALL dbcsr_get_block_p(real_k(ispin, first_k + ik - 1)%matrix, plan%row(first), plan%col(first), &
                                            rblock, found)
                     IF (found) THEN
                        CPASSERT(SIZE(rblock, 1) == nrow .AND. SIZE(rblock, 2) == ncol)
                        result_block(1:nrow, 1:ncol) => RESULT((ispin - 1)*nelem + 1:ispin*nelem, ik)
                        rblock = result_block
                     END IF
                     CALL dbcsr_get_block_p(imag_k(ispin, first_k + ik - 1)%matrix, plan%row(first), plan%col(first), &
                                            iblock, found)
                     IF (found) THEN
                        CPASSERT(SIZE(iblock, 1) == nrow .AND. SIZE(iblock, 2) == ncol)
                        result_block(1:nrow, 1:ncol) => RESULT((ispin - 1)*nelem + 1:ispin*nelem, width + ik)
                        iblock = result_block
                     END IF
                  END DO
               END DO
            END DO
         ELSE
            per_tile_bytes = 16_int_8*nspin*nk
            IF (per_tile_bytes > fft_tile_budget) THEN
               CPABORT("R-to-K DBCSR FFT element tile exceeds the memory budget")
            END IF
            tile_width = MIN(nelem, INT(fft_tile_budget/per_tile_bytes))
            fft_input(1:nelem, 1:nspin, 1:nimages) => packed(1:nelem*nspin*nimages)
            ALLOCATE (fft_output(tile_width, nspin, nk))
            DO tile_first = 1, nelem, tile_width
               tile_size = MIN(tile_width, nelem - tile_first + 1)
               CALL cell_to_k_grid_fft(fft_input(tile_first:tile_first + tile_size - 1, :, :), &
                                       signed_cells(:, first:last), xkp, &
                                       nkp_grid, fft_output(1:tile_size, :, :), require_fft=.TRUE., work=fft_work, &
                                       cell_first=first)
               DO ik = 1, nk
                  DO ispin = 1, nspin
                     CALL dbcsr_get_block_p(real_k(ispin, ik)%matrix, plan%row(first), plan%col(first), &
                                            rblock, found_re)
                     IF (found_re) THEN
                        CPASSERT(SIZE(rblock, 1) == nrow .AND. SIZE(rblock, 2) == ncol)
                     END IF
                     CALL dbcsr_get_block_p(imag_k(ispin, ik)%matrix, plan%row(first), plan%col(first), &
                                            iblock, found_im)
                     IF (found_im) THEN
                        CPASSERT(SIZE(iblock, 1) == nrow .AND. SIZE(iblock, 2) == ncol)
                     END IF
                     IF (.NOT. (found_re .OR. found_im)) CYCLE
                     ! Write each column segment directly; FFT tiles may split a column.
                     DO col = (tile_first - 1)/nrow + 1, (tile_first + tile_size - 2)/nrow + 1
                        offset = (col - 1)*nrow
                        row_first = MAX(1, tile_first - offset)
                        row_last = MIN(nrow, tile_first + tile_size - 1 - offset)
                        pos = offset + row_first - tile_first + 1
                        IF (found_re) rblock(row_first:row_last, col) = &
                           REAL(fft_output(pos:pos + row_last - row_first, ispin, ik), dp)
                        IF (found_im) iblock(row_first:row_last, col) = &
                           AIMAG(fft_output(pos:pos + row_last - row_first, ispin, ik))
                     END DO
                  END DO
               END DO
            END DO
            DEALLOCATE (fft_output)
         END IF
      END DO
      CALL cell_to_k_grid_release(fft_work)
      CALL timestop(handle)
   END SUBROUTINE rskp_dbcsr_transform_batch

! **************************************************************************************************
!> \brief Build the immutable neighbor-list traversal shared by R-to-K and K-to-R transforms.
!> \param plan ...
!> \param sab_nl ...
!> \param cell_to_index ...
!> \param nimg ...
!> \param block_template ...
!> \param group_entries ...
! **************************************************************************************************
   SUBROUTINE kp_transform_plan_create(plan, sab_nl, cell_to_index, nimg, block_template, group_entries)

      TYPE(kp_transform_plan_type), INTENT(OUT)          :: plan
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      INTEGER, INTENT(IN)                                :: nimg
      TYPE(dbcsr_type), INTENT(IN), OPTIONAL             :: block_template
      LOGICAL, INTENT(IN), OPTIONAL                      :: group_entries

      INTEGER                                            :: i, iatom, icell, icol, igroup, irow, &
                                                            jatom, nblock
      INTEGER(KIND=int_8), ALLOCATABLE, DIMENSION(:)     :: key
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: order
      INTEGER, DIMENSION(3)                              :: cell
      INTEGER, DIMENSION(:), POINTER                     :: col_offsets, row_offsets
      LOGICAL                                            :: do_grouping, do_symmetric, store_offsets
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: iterator

      CALL get_neighbor_list_set_p(neighbor_list_sets=sab_nl, symmetric=do_symmetric)
      plan%symmetric = do_symmetric
      do_grouping = .FALSE.
      IF (PRESENT(group_entries)) do_grouping = group_entries
      store_offsets = PRESENT(block_template)
      IF (store_offsets) THEN
         CALL dbcsr_get_info(block_template, row_blk_offset=row_offsets, col_blk_offset=col_offsets)
      END IF
      CALL neighbor_list_iterator_create(iterator, sab_nl)
      DO WHILE (neighbor_list_iterate(iterator) == 0)
         CALL get_iterator_info(iterator, cell=cell)
         icell = cell_to_index(cell(1), cell(2), cell(3))
         IF (icell >= 1 .AND. icell <= nimg) plan%nentry = plan%nentry + 1
      END DO
      CALL neighbor_list_iterator_release(iterator)

      ALLOCATE (plan%row(plan%nentry), plan%col(plan%nentry), plan%image(plan%nentry), &
                plan%cell(3, plan%nentry), plan%symmetry_sign(plan%nentry))
      IF (store_offsets) THEN
         ALLOCATE (plan%row_offset(plan%nentry), plan%col_offset(plan%nentry))
      END IF

      i = 0
      CALL neighbor_list_iterator_create(iterator, sab_nl)
      DO WHILE (neighbor_list_iterate(iterator) == 0)
         CALL get_iterator_info(iterator, iatom=iatom, jatom=jatom, cell=cell)
         icell = cell_to_index(cell(1), cell(2), cell(3))
         IF (icell < 1 .OR. icell > nimg) CYCLE
         i = i + 1
         irow = iatom
         icol = jatom
         plan%symmetry_sign(i) = 1.0_dp
         IF (do_symmetric .AND. iatom > jatom) THEN
            irow = jatom
            icol = iatom
            plan%symmetry_sign(i) = -1.0_dp
         END IF
         plan%row(i) = irow
         plan%col(i) = icol
         IF (store_offsets) THEN
            plan%row_offset(i) = row_offsets(irow)
            plan%col_offset(i) = col_offsets(icol)
         END IF
         plan%image(i) = icell
         plan%cell(:, i) = cell
      END DO
      CALL neighbor_list_iterator_release(iterator)
      CPASSERT(i == plan%nentry)

      IF (do_grouping .AND. plan%nentry > 0) THEN
         nblock = MAX(MAXVAL(plan%row), MAXVAL(plan%col))
         ALLOCATE (key(plan%nentry), order(plan%nentry))
         DO i = 1, plan%nentry
            key(i) = INT(plan%image(i), KIND=int_8) + INT(nimg, KIND=int_8)* &
                     (INT(plan%row(i) - 1, KIND=int_8) + INT(nblock, KIND=int_8)* &
                      INT(plan%col(i) - 1, KIND=int_8))
         END DO
         CALL sort(key, plan%nentry, order)
         plan%row(:) = plan%row(order)
         plan%col(:) = plan%col(order)
         IF (store_offsets) THEN
            plan%row_offset(:) = plan%row_offset(order)
            plan%col_offset(:) = plan%col_offset(order)
         END IF
         plan%image(:) = plan%image(order)
         plan%cell(:, :) = plan%cell(:, order)
         plan%symmetry_sign(:) = plan%symmetry_sign(order)

         plan%ngroup = 1
         DO i = 2, plan%nentry
            IF (key(i) /= key(i - 1)) plan%ngroup = plan%ngroup + 1
         END DO
         ALLOCATE (plan%group_start(plan%ngroup + 1))
         igroup = 1
         plan%group_start(igroup) = 1
         DO i = 2, plan%nentry
            IF (key(i) /= key(i - 1)) THEN
               igroup = igroup + 1
               plan%group_start(igroup) = i
            END IF
         END DO
         plan%group_start(plan%ngroup + 1) = plan%nentry + 1
         DEALLOCATE (key, order)
      END IF

   END SUBROUTINE kp_transform_plan_create

! **************************************************************************************************
!> \brief Contract one real-space image panel with real and optional imaginary phase panels.
!> \param nrow number of AO elements in the panel
!> \param nimages number of real-space images in the panel
!> \param nk number of k-point columns
!> \param values_re image values contributing to the real output
!> \param phase_re real phase factors
!> \param result_re real output panel
!> \param beta output accumulation factor
!> \param values_im optional image values contributing to the imaginary output
!> \param phase_im optional imaginary phase factors
!> \param result_im optional imaginary output panel
! **************************************************************************************************
   SUBROUTINE rskp_gemm_contract(nrow, nimages, nk, values_re, phase_re, result_re, beta, &
                                 values_im, phase_im, result_im)
      INTEGER, INTENT(IN)                                :: nrow, nimages, nk
      REAL(KIND=dp), CONTIGUOUS, INTENT(IN)              :: values_re(:, :), phase_re(:, :)
      REAL(KIND=dp), CONTIGUOUS, INTENT(INOUT)           :: result_re(:, :)
      REAL(KIND=dp), INTENT(IN)                          :: beta
      REAL(KIND=dp), CONTIGUOUS, INTENT(IN), OPTIONAL    :: values_im(:, :), phase_im(:, :)
      REAL(KIND=dp), CONTIGUOUS, INTENT(INOUT), OPTIONAL :: result_im(:, :)

      CPASSERT(nrow > 0 .AND. nimages > 0 .AND. nk > 0)
      CPASSERT(SIZE(values_re, 1) >= nrow .AND. SIZE(values_re, 2) >= nimages)
      CPASSERT(SIZE(phase_re, 1) >= nimages .AND. SIZE(phase_re, 2) >= nk)
      CPASSERT(SIZE(result_re, 1) >= nrow .AND. SIZE(result_re, 2) >= nk)
      CPASSERT(PRESENT(values_im) .EQV. PRESENT(phase_im))
      CPASSERT(PRESENT(values_im) .EQV. PRESENT(result_im))

      CALL dgemm('N', 'N', nrow, nk, nimages, 1.0_dp, values_re, SIZE(values_re, 1), &
                 phase_re, SIZE(phase_re, 1), beta, result_re, SIZE(result_re, 1))
      IF (PRESENT(values_im)) THEN
         CPASSERT(SIZE(values_im, 1) >= nrow .AND. SIZE(values_im, 2) >= nimages)
         CPASSERT(SIZE(phase_im, 1) >= nimages .AND. SIZE(phase_im, 2) >= nk)
         CPASSERT(SIZE(result_im, 1) >= nrow .AND. SIZE(result_im, 2) >= nk)
         CALL dgemm('N', 'N', nrow, nk, nimages, 1.0_dp, values_im, SIZE(values_im, 1), &
                    phase_im, SIZE(phase_im, 1), beta, result_im, SIZE(result_im, 1))
      END IF
   END SUBROUTINE rskp_gemm_contract

! **************************************************************************************************
!> \brief Cache routing to the FM tiles of one representative group.
!>        Equivalent ranks in the other groups receive only their tile's blocks.
!>        Both backends receive the union of locally owned direct/mirrored AO elements.
!>        Missing source blocks are packed as zero on a fixed geometry.
!> \param source operator layout template
!> \param source_sab neighbor lists belonging to source
!> \param cell_to_index map from lattice translations to source image indices
!> \param matrix_struct target group or world FM layout
!> \param world communicator distributing the source DBCSR matrices
!> \param group communicator of group-local FM targets
!> \param inter communicator joining equivalent FM ranks across groups
!> \param backend GEMM or FFT
!> \param max_bytes per-rank budget for routing and capture buffers
!> \param route context-owned structural route
!> \param ready false collectively if unsupported or over budget
! **************************************************************************************************
   SUBROUTINE rskp_route_create(source, source_sab, cell_to_index, matrix_struct, world, group, inter, &
                                backend, max_bytes, route, ready)

      TYPE(dbcsr_p_type), POINTER                        :: source(:, :)
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: source_sab
      INTEGER, POINTER                                   :: cell_to_index(:, :, :)
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct
      TYPE(mp_para_env_type), POINTER                    :: world, group, inter
      INTEGER, INTENT(IN)                                :: backend
      INTEGER(KIND=int_8), INTENT(IN)                    :: max_bytes
      TYPE(rskp_route_type), INTENT(OUT)                 :: route
      LOGICAL, INTENT(OUT)                               :: ready

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'rskp_route_create'

      INTEGER                                            :: compatible, dest, handle, i, iblock, &
                                                            layout(11), layout_max(11), &
                                                            layout_min(11), nblock, ncol, nowner, &
                                                            nrow, p, pos, rank
      INTEGER(KIND=int_8)                                :: nr(2), ns(2), storage
      INTEGER(KIND=int_8), ALLOCATABLE                   :: recv_size(:, :), send_size(:, :)
      INTEGER, ALLOCATABLE                               :: identities(:, :), owner(:), &
                                                            position(:, :), send_metadata(:, :), &
                                                            value_count(:), value_start(:)
      INTEGER, POINTER                                   :: col_size(:), row_size(:)
      LOGICAL, ALLOCATABLE                               :: needed(:)

      CALL timeset(routineN, handle)
      ready = .FALSE.
      route%world => world
      route%group => group
      route%inter => inter
      CALL dbcsr_get_info(source(1, 1)%matrix, row_blk_size=row_size, col_blk_size=col_size)
      CALL kp_transform_plan_create(route%traversal, source_sab, cell_to_index, &
                                    SIZE(source, 2), block_template=source(1, 1)%matrix)
      route%source_sab => source_sab
      route%transform%backend = backend
      route%world_target = ASSOCIATED(matrix_struct%para_env, world)
      nowner = group%num_pe
      IF (route%world_target) nowner = world%num_pe
      compatible = MERGE(1, 0, route%world_target .OR. &
                         ASSOCIATED(matrix_struct%para_env, group))
      IF (.NOT. route%world_target) THEN
         layout = [matrix_struct%nrow_global, matrix_struct%ncol_global, &
                   matrix_struct%nrow_block, matrix_struct%ncol_block, matrix_struct%first_p_pos, &
                   matrix_struct%context%num_pe, matrix_struct%context%mepos, matrix_struct%local_leading_dimension]
         layout_min = layout
         layout_max = layout
         CALL inter%min(layout_min)
         CALL inter%max(layout_max)
         IF (ANY(layout_min /= layout_max)) compatible = 0
      END IF
      IF (SUM(row_size) /= matrix_struct%nrow_global .OR. &
          SUM(col_size) /= matrix_struct%ncol_global) compatible = 0
      IF (route%traversal%symmetric) THEN
         IF (matrix_struct%nrow_global /= matrix_struct%ncol_global .OR. SIZE(row_size) /= SIZE(col_size)) THEN
            compatible = 0
         ELSE IF (ANY(row_size /= col_size)) THEN
            compatible = 0
         END IF
      END IF
      CALL world%min(compatible)
      IF (compatible == 0) THEN
         CALL cp_warn(__LOCATION__, &
                      "Prepared R-to-K transform has an incompatible source or target layout")
         CALL timestop(handle)
         RETURN
      END IF
      route%nimage = SIZE(source, 2)
      ALLOCATE (route%row_size(SIZE(row_size)), route%col_size(SIZE(col_size)))
      route%row_size(:) = row_size
      route%col_size(:) = col_size
      ALLOCATE (identities(2, 0:world%num_pe - 1), owner(0:nowner - 1), &
                needed(0:nowner - 1), value_count(0:nowner - 1), value_start(0:nowner - 1))
      CALL world%allgather([group%mepos, inter%mepos], identities)
      owner(:) = -1
      IF (route%world_target) THEN
         DO rank = 0, world%num_pe - 1
            owner(rank) = rank
         END DO
      ELSE
         DO rank = 0, world%num_pe - 1
            IF (identities(2, rank) == 0) owner(identities(1, rank)) = rank
         END DO
      END IF
      CPASSERT(ALL(owner >= 0))
      ALLOCATE (send_size(2, 0:world%num_pe - 1), recv_size(2, 0:world%num_pe - 1))
      send_size(:, :) = 0_int_8
      DO i = 1, route%traversal%nentry
         nrow = row_size(route%traversal%row(i))
         ncol = col_size(route%traversal%col(i))
         CALL rskp_owners(route, matrix_struct, i, needed, value_count)
         DO p = 0, SIZE(owner) - 1
            IF (.NOT. needed(p)) CYCLE
            dest = owner(p)
            send_size(1, dest) = send_size(1, dest) + 1_int_8
            send_size(2, dest) = send_size(2, dest) + INT(value_count(p), int_8)
         END DO
      END DO
      CALL world%alltoall(send_size, recv_size, 2)
      ! Only inter-rank zero receives from the world. The other groups hold
      ! equivalent FM tiles, as required by kpoint_env_initialize.
      ns = SUM(send_size, DIM=2)
      nr = SUM(recv_size, DIM=2)
      IF (.NOT. route%world_target) CALL inter%bcast(nr, 0)
      ! Conservative peak: traversal, retained route, temporary metadata and
      ! one snapshot. The full batch budget is checked again by prepare.
      storage = 40_int_8*route%traversal%nentry + &
                4_int_8*(SIZE(row_size) + SIZE(col_size) + matrix_struct%nrow_global + matrix_struct%ncol_global) + &
                48_int_8*ns(1) + MERGE(0_int_8, 8_int_8*ns(2), world%num_pe == 1) + 68_int_8*nr(1) + &
                128_int_8*world%num_pe + 16_int_8*SIZE(owner) + 8_int_8
      ! Singleton groups also retain the AO-block group boundaries.
      IF (SIZE(owner) == 1) storage = storage + 4_int_8*(nr(1) + 1_int_8)
      ! Group permutations, packed offsets, direct/mirror element maps and
      ! sorting scratch. The cropped receive buffer is shared by all slots.
      storage = storage + 52_int_8*(nr(1) + 1_int_8) + 16_int_8*nr(2)
      IF (world%num_pe > 1) storage = storage + 8_int_8*nr(2)
      IF (SIZE(owner) > 1) storage = storage + 4_int_8*ns(2)
      route%storage = storage
      storage = storage + 128_int_8*nr(1)
      ! Temporary canonical-image classification and one packed snapshot.
      storage = storage + 128_int_8*nr(1) + 8_int_8
      ! Borrowed values do not consume the snapshot budget; still bound integer route offsets.
      IF (MAXVAL(ns) >= INT(HUGE(0), int_8) .OR. MAXVAL(nr) >= INT(HUGE(0), int_8)) THEN
         storage = MAX(storage, max_bytes + 1_int_8)
      END IF
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
      IF (SIZE(owner) == 1 .AND. .NOT. route%world_target) CALL inter%bcast(route%rcount, 0)
      route%sdisp(:, 0) = 0
      route%rdisp(:, 0) = 0
      DO rank = 1, world%num_pe - 1
         route%sdisp(:, rank) = route%sdisp(:, rank - 1) + route%scount(:, rank - 1)
         route%rdisp(:, rank) = route%rdisp(:, rank - 1) + route%rcount(:, rank - 1)
      END DO
      ALLOCATE (route%send_entry(INT(ns(1))), route%send_first(INT(ns(1)) + 1), &
                route%send_values(MERGE(0, INT(ns(2)), world%num_pe == 1)), send_metadata(8, INT(ns(1))), &
                position(2, 0:world%num_pe - 1), route%entry(8, nblock), &
                route%first(nblock + 1))
      IF (SIZE(owner) > 1) ALLOCATE (route%send_element(INT(ns(2))))
      IF (world%num_pe > 1) ALLOCATE (route%recv_values(route%nvalue))
      position(:, :) = route%sdisp
      DO i = 1, route%traversal%nentry
         nrow = row_size(route%traversal%row(i))
         ncol = col_size(route%traversal%col(i))
         IF (ALLOCATED(route%send_element)) THEN
            value_start(:) = position(2, owner)
            CALL rskp_owners(route, matrix_struct, i, needed, value_count, &
                             value_start, route%send_element)
         ELSE
            CALL rskp_owners(route, matrix_struct, i, needed, value_count)
         END IF
         DO p = 0, SIZE(owner) - 1
            IF (.NOT. needed(p)) CYCLE
            dest = owner(p)
            iblock = position(1, dest) + 1
            pos = position(2, dest) + 1
            route%send_entry(iblock) = i
            route%send_first(iblock) = pos
            send_metadata(1:4, iblock) = [route%traversal%row_offset(i), route%traversal%col_offset(i), nrow, ncol]
            send_metadata(5:7, iblock) = NINT(route%traversal%symmetry_sign(i))*route%traversal%cell(:, i)
            send_metadata(8, iblock) = MERGE(1, 0, route%traversal%symmetric .AND. &
                                             route%traversal%row(i) /= route%traversal%col(i))
            position(:, dest) = position(:, dest) + [1, value_count(p)]
         END DO
      END DO
      route%send_first(INT(ns(1)) + 1) = INT(ns(2)) + 1
      IF (SIZE(owner) == 1 .AND. world%num_pe == 1) THEN
         CALL world%allgatherv(send_metadata, route%entry, 8*route%rcount(1, :), 8*route%rdisp(1, :))
      ELSE IF (.NOT. route%world_target .AND. SIZE(owner) == 1) THEN
         CALL world%allgatherv(send_metadata, route%entry, 8*route%rcount(1, :), 8*route%rdisp(1, :))
      ELSE
         CALL world%alltoall(send_metadata, 8*route%scount(1, :), 8*route%sdisp(1, :), &
                             route%entry, 8*route%rcount(1, :), 8*route%rdisp(1, :))
         IF (.NOT. route%world_target) CALL inter%bcast(route%entry, 0)
      END IF
      route%matrix_struct => matrix_struct
      CALL cp_fm_struct_retain(route%matrix_struct)
      CALL rskp_tile_map(route)
      ready = .TRUE.
      CALL timestop(handle)

   END SUBROUTINE rskp_route_create

! **************************************************************************************************
!> \brief Map routed blocks to local FM ownership and stable atom-pair groups.
!>        FM indices are ordered globally; dummy indices on empty tiles lie outside the matrix.
!> \param route fixed routing and target distribution
! **************************************************************************************************
   SUBROUTINE rskp_tile_map(route)
      TYPE(rskp_route_type), INTENT(INOUT)               :: route

      INTEGER                                            :: c0, col, element, first, g, i, image, m, &
                                                            nc, ngroup, nr, overlap, pos, r0, row
      INTEGER(KIND=int_8), ALLOCATABLE                   :: key(:)
      INTEGER, ALLOCATABLE                               :: element_count(:)
      LOGICAL                                            :: direct, mirror
      INTEGER :: rows(0:route%matrix_struct%nrow_global), cols(0:route%matrix_struct%ncol_global)

      rows(:) = 0
      cols(:) = 0
      DO i = 1, SIZE(route%matrix_struct%row_indices)
         r0 = route%matrix_struct%row_indices(i)
         IF (r0 >= 1 .AND. r0 <= UBOUND(rows, 1)) rows(r0) = 1
      END DO
      DO i = 1, SIZE(route%matrix_struct%col_indices)
         c0 = route%matrix_struct%col_indices(i)
         IF (c0 >= 1 .AND. c0 <= UBOUND(cols, 1)) cols(c0) = 1
      END DO
      DO i = 1, UBOUND(rows, 1)
         rows(i) = rows(i) + rows(i - 1)
      END DO
      DO i = 1, UBOUND(cols, 1)
         cols(i) = cols(i) + cols(i - 1)
      END DO
      route%nrow = rows(UBOUND(rows, 1))
      route%ncol = cols(UBOUND(cols, 1))
      ALLOCATE (route%block_order(SIZE(route%entry, 2)), element_count(SIZE(route%entry, 2)))
      BLOCK
         ! Source-rank reception need not keep a pair contiguous on distributed FM layouts.
         ALLOCATE (key(SIZE(route%entry, 2)))
         DO i = 1, SIZE(key)
            key(i) = INT(route%entry(1, i), int_8) + &
                     INT(route%matrix_struct%nrow_global, int_8)*(route%entry(2, i) - 1)
         END DO
         IF (SIZE(key) > 0) CALL sort(key, SIZE(key), route%block_order)
         ngroup = 0
         IF (SIZE(key) > 0) ngroup = 1 + COUNT(key(2:) /= key(:SIZE(key) - 1))
         ALLOCATE (route%block_first(ngroup + 1))
         route%block_first(1) = 1
         g = 1
         DO i = 2, SIZE(key)
            IF (key(i) == key(i - 1)) CYCLE
            g = g + 1
            route%block_first(g) = i
         END DO
         route%block_first(ngroup + 1) = SIZE(key) + 1
      END BLOCK
      ALLOCATE (route%element_first(ngroup + 1))
      route%element_first(1) = 1
      pos = 1
      DO g = 1, ngroup
         first = route%block_order(route%block_first(g))
         r0 = route%entry(1, first)
         c0 = route%entry(2, first)
         nr = route%entry(3, first)
         nc = route%entry(4, first)
         m = (rows(r0 + nr - 1) - rows(r0 - 1))*(cols(c0 + nc - 1) - cols(c0 - 1))
         IF (route%entry(8, first) /= 0) THEN
            ! Subtract the intersection: a rank owning both orientations
            ! needs one input value and one GEMM row for that element.
            overlap = COUNT(rows(r0:r0 + nr - 1) > rows(r0 - 1:r0 + nr - 2) .AND. &
                            cols(r0:r0 + nr - 1) > cols(r0 - 1:r0 + nr - 2))* &
                      COUNT(rows(c0:c0 + nc - 1) > rows(c0 - 1:c0 + nc - 2) .AND. &
                            cols(c0:c0 + nc - 1) > cols(c0 - 1:c0 + nc - 2))
            m = m + (rows(c0 + nc - 1) - rows(c0 - 1))*(cols(r0 + nr - 1) - cols(r0 - 1)) - overlap
         END IF
         route%element_first(g + 1) = route%element_first(g) + m
         DO i = route%block_first(g), route%block_first(g + 1) - 1
            image = route%block_order(i)
            element_count(image) = m
            pos = pos + m
         END DO
      END DO
      route%first(1) = 1
      DO i = 1, SIZE(element_count)
         route%first(i + 1) = route%first(i) + element_count(i)
      END DO
      CPASSERT(pos == route%nvalue + 1)
      CPASSERT(route%first(SIZE(element_count) + 1) == pos)
      ALLOCATE (route%element_offset(2, route%element_first(ngroup + 1) - 1))
      route%element_offset(:, :) = 0_int_8
      DO g = 1, ngroup
         first = route%block_order(route%block_first(g))
         r0 = route%entry(1, first)
         c0 = route%entry(2, first)
         nr = route%entry(3, first)
         nc = route%entry(4, first)
         element = route%element_first(g)
         DO col = c0, c0 + nc - 1
            DO row = r0, r0 + nr - 1
               direct = rows(row) > rows(row - 1) .AND. cols(col) > cols(col - 1)
               mirror = .FALSE.
               IF (route%entry(8, first) /= 0) THEN
                  mirror = rows(col) > rows(col - 1) .AND. cols(row) > cols(row - 1)
               END IF
               IF (.NOT. (direct .OR. mirror)) CYCLE
               IF (direct) route%element_offset(1, element) = INT(rows(row), int_8) + &
                                                              INT(cols(col) - 1, int_8)*route%matrix_struct%local_leading_dimension
               IF (mirror) route%element_offset(2, element) = INT(rows(col), int_8) + &
                                                              INT(cols(row) - 1, int_8)*route%matrix_struct%local_leading_dimension
               element = element + 1
            END DO
         END DO
         CPASSERT(element == route%element_first(g + 1))
      END DO
      IF (route%transform%backend == rskp_backend_fft) THEN
         ALLOCATE (route%transform%packed_first(ngroup + 1), route%transform%fold_map(2, SIZE(route%entry, 2)), source=0)
         pos = 1
         DO g = 1, ngroup
            route%transform%packed_first(g) = pos
            m = route%element_first(g + 1) - route%element_first(g)
            DO i = route%block_first(g), route%block_first(g + 1) - 1
               route%transform%fold_map(1, route%block_order(i)) = pos
               pos = pos + m
            END DO
         END DO
         route%transform%packed_first(ngroup + 1) = pos
      END IF
   END SUBROUTINE rskp_tile_map

! **************************************************************************************************
!> \brief Find group-local FM owners and optionally pack their cropped source indices.
!> \param route source AO block layout and neighbor traversal
!> \param fms target FM layout within each group
!> \param i traversal entry
!> \param needed recipient group-local ranks, without duplicate mirror sends
!> \param value_count cropped value counts by target rank
!> \param value_start optional offsets into the value-index buffer
!> \param value_index optional column-major source indices, packed by recipient
! **************************************************************************************************
   SUBROUTINE rskp_owners(route, fms, i, needed, value_count, value_start, value_index)

      TYPE(rskp_route_type), INTENT(IN)                  :: route
      TYPE(cp_fm_struct_type), INTENT(IN)                :: fms
      INTEGER, INTENT(IN)                                :: i
      LOGICAL, INTENT(OUT)                               :: needed(0:)
      INTEGER, INTENT(OUT)                               :: value_count(0:)
      INTEGER, INTENT(IN), OPTIONAL                      :: value_start(0:)
      INTEGER, INTENT(INOUT), OPTIONAL                   :: value_index(:)

      INTEGER                                            :: c0, col, direct, mirror, nc, nr, r0, row
      LOGICAL                                            :: mirror_block

      r0 = route%traversal%row_offset(i)
      c0 = route%traversal%col_offset(i)
      nr = route%row_size(route%traversal%row(i))
      nc = route%col_size(route%traversal%col(i))
      mirror_block = route%traversal%symmetric .AND. &
                     route%traversal%row(i) /= route%traversal%col(i)
      CPASSERT(PRESENT(value_start) .EQV. PRESENT(value_index))
      value_count(:) = 0
      IF (SIZE(needed) == 1) THEN
         needed(:) = .TRUE.
         value_count(0) = nr*nc
         RETURN
      END IF
      ! Send each canonical source element at most once to a recipient,
      ! even when that rank owns both its direct and mirrored destinations.
      DO col = 0, nc - 1
         DO row = 0, nr - 1
            direct = fms%context%blacs2mpi(fms%g2p_row(r0 + row), fms%g2p_col(c0 + col))
            mirror = direct
            IF (mirror_block) THEN
               mirror = fms%context%blacs2mpi(fms%g2p_row(c0 + col), fms%g2p_col(r0 + row))
            END IF
            value_count(direct) = value_count(direct) + 1
            IF (PRESENT(value_index)) THEN
               value_index(value_start(direct) + value_count(direct)) = row + col*nr + 1
            END IF
            IF (mirror == direct) CYCLE
            value_count(mirror) = value_count(mirror) + 1
            IF (PRESENT(value_index)) THEN
               value_index(value_start(mirror) + value_count(mirror)) = row + col*nr + 1
            END IF
         END DO
      END DO
      needed(:) = value_count > 0
   END SUBROUTINE rskp_owners

! **************************************************************************************************
!> \brief Check the source AO layout, not its values or sparse block presence.
!> \param route fixed context route
!> \param source current operator images
!> \param matrix_row row to capture
!> \param source_sab source neighbor list, required to match the route's geometry
!> \return whether the row has the route's image count and AO layout
! **************************************************************************************************
   FUNCTION rskp_compatible(route, source, matrix_row, source_sab) RESULT(compatible)

      TYPE(rskp_route_type), INTENT(IN)                  :: route
      TYPE(dbcsr_p_type), POINTER                        :: source(:, :)
      INTEGER, INTENT(IN)                                :: matrix_row
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: source_sab
      LOGICAL                                            :: compatible

      INTEGER                                            :: image
      INTEGER, POINTER                                   :: col_size(:), row_size(:)

      compatible = .FALSE.
      IF (.NOT. ASSOCIATED(route%source_sab, source_sab)) RETURN
      IF (SIZE(source, 2) /= route%nimage) RETURN
      IF (matrix_row < 1 .OR. matrix_row > SIZE(source, 1)) RETURN
      DO image = 1, SIZE(source, 2)
         CALL dbcsr_get_info(source(matrix_row, image)%matrix, row_blk_size=row_size, col_blk_size=col_size)
         IF (SIZE(row_size) /= SIZE(route%row_size) .OR. SIZE(col_size) /= SIZE(route%col_size)) RETURN
         IF (ANY(row_size /= route%row_size) .OR. ANY(col_size /= route%col_size)) RETURN
      END DO
      compatible = .TRUE.

   END FUNCTION rskp_compatible

! **************************************************************************************************
!> \brief Capture current values using cached routing, with no per-k communication.
!>        GEMM/FFT snapshots own packed local AO-element/image matrices. A shared
!>        cropped receive buffer is reused across slots and operator refreshes.
!> \param source current real-space operator images
!> \param matrix_row row to snapshot
!> \param route context-owned metadata, counts and reusable packing buffer
!> \param snapshot immutable tile-local values
!> \param mirror_sign transpose parity of the unscaled operator (+1 or -1)
!> \param scale complex prefactor, applied after forming direct and mirrored values
! **************************************************************************************************
   SUBROUTINE rskp_capture(source, matrix_row, route, snapshot, mirror_sign, scale)

      TYPE(dbcsr_p_type), POINTER                        :: source(:, :)
      INTEGER, INTENT(IN)                                :: matrix_row
      TYPE(rskp_route_type), INTENT(INOUT), TARGET       :: route
      TYPE(rskp_values_type), INTENT(OUT), TARGET        :: snapshot
      INTEGER, INTENT(IN), OPTIONAL                      :: mirror_sign
      COMPLEX(KIND=dp), INTENT(IN), OPTIONAL             :: scale

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'rskp_capture'

      INTEGER                                            :: col, element, g, handle, i, iblock, m, &
                                                            ncol, nrow, odd, panels(4), pos, row
      LOGICAL                                            :: found
      REAL(KIND=dp)                                      :: sense
      REAL(KIND=dp), CONTIGUOUS, POINTER                 :: received(:)
      REAL(KIND=dp), POINTER                             :: block(:, :)
      TYPE(mp_para_env_type), POINTER                    :: inter, world

      CALL timeset(routineN, handle)
      world => route%world
      inter => route%inter
      m = route%transform%packed_first(SIZE(route%transform%packed_first)) - 1
      CPASSERT(m <= route%nvalue)
      ALLOCATE (snapshot%packed(m), source=0.0_dp)
      snapshot%source => source
      snapshot%matrix_row = matrix_row
      IF (PRESENT(mirror_sign)) snapshot%mirror_sign = mirror_sign
      CPASSERT(ABS(snapshot%mirror_sign) == 1)
      IF (PRESENT(scale)) snapshot%scale = scale
      IF (world%num_pe == 1) THEN
         ! Add each image to its fixed raw-FFT or parity-GEMM panels.
         DO iblock = 1, SIZE(route%entry, 2)
            i = route%send_entry(iblock)
            CALL dbcsr_get_readonly_block_p(source(matrix_row, route%traversal%image(i))%matrix, &
                                            route%traversal%row(i), route%traversal%col(i), block, found)
            IF (.NOT. found) CYCLE
            pos = route%transform%fold_map(1, iblock)
            odd = ABS(route%transform%fold_map(2, iblock))
            sense = REAL(SIGN(1, route%transform%fold_map(2, iblock)), dp)
            DO col = 1, SIZE(block, 2)
               DO row = 1, SIZE(block, 1)
                  snapshot%packed(pos) = snapshot%packed(pos) + block(row, col)
                  IF (odd > 0) THEN
                     snapshot%packed(odd) = snapshot%packed(odd) + sense*block(row, col)
                     odd = odd + 1
                  END IF
                  pos = pos + 1
               END DO
            END DO
         END DO
      ELSE
         received => route%recv_values
         route%send_values(:) = 0.0_dp
!$OMP PARALLEL DO DEFAULT(NONE) SHARED(source, matrix_row, route) &
!$OMP PRIVATE(iblock, i, block, found, nrow, ncol, pos, row, col, element) SCHEDULE(STATIC)
         DO iblock = 1, SIZE(route%send_entry)
            i = route%send_entry(iblock)
            CALL dbcsr_get_readonly_block_p(source(matrix_row, route%traversal%image(i))%matrix, &
                                            route%traversal%row(i), route%traversal%col(i), block, found)
            IF (.NOT. found) CYCLE
            nrow = route%row_size(route%traversal%row(i))
            ncol = route%col_size(route%traversal%col(i))
            CPASSERT(SIZE(block, 1) == nrow .AND. SIZE(block, 2) == ncol)
            pos = route%send_first(iblock)
            IF (ALLOCATED(route%send_element)) THEN
               DO element = pos, route%send_first(iblock + 1) - 1
                  row = MODULO(route%send_element(element) - 1, nrow) + 1
                  col = (route%send_element(element) - 1)/nrow + 1
                  route%send_values(element) = block(row, col)
               END DO
               CYCLE
            END IF
            DO col = 1, ncol
               DO row = 1, nrow
                  route%send_values(pos) = block(row, col)
                  pos = pos + 1
               END DO
            END DO
         END DO
!$OMP END PARALLEL DO
         IF (world%num_pe == 1 .OR. &
             (.NOT. route%world_target .AND. route%group%num_pe == 1)) THEN
            CALL world%allgatherv(route%send_values, received, route%rcount(2, :), route%rdisp(2, :))
         ELSE
            CALL world%alltoall(route%send_values, route%scount(2, :), route%sdisp(2, :), &
                                received, route%rcount(2, :), route%rdisp(2, :))
            IF (.NOT. route%world_target) CALL inter%bcast(received, 0)
         END IF
         BLOCK
            DO iblock = 1, SIZE(route%entry, 2)
               pos = route%transform%fold_map(1, iblock)
               odd = ABS(route%transform%fold_map(2, iblock))
               m = route%first(iblock + 1) - route%first(iblock)
               snapshot%packed(pos:pos + m - 1) = snapshot%packed(pos:pos + m - 1) + &
                                                  received(route%first(iblock):route%first(iblock + 1) - 1)
               IF (odd == 0) CYCLE
               sense = REAL(SIGN(1, route%transform%fold_map(2, iblock)), dp)
               snapshot%packed(odd:odd + m - 1) = snapshot%packed(odd:odd + m - 1) + &
                                                  sense*received(route%first(iblock):route%first(iblock + 1) - 1)
            END DO
         END BLOCK
      END IF
      IF (route%transform%backend == rskp_backend_gemm) THEN
         ALLOCATE (snapshot%active(4, SIZE(route%transform%phase_count, 2)))
         DO g = 1, SIZE(snapshot%active, 2)
            m = route%element_first(g + 1) - route%element_first(g)
            panels = route%transform%phase_count([1, 1, 2, 3], g)
            pos = route%transform%packed_first(g)
            DO i = 1, 4
               snapshot%active(i, g) = ANY(snapshot%packed(pos:pos + m*panels(i) - 1) /= 0.0_dp)
               pos = pos + m*panels(i)
            END DO
         END DO
      END IF
      CALL timestop(handle)

   END SUBROUTINE rskp_capture

! **************************************************************************************************
!> \brief Release one captured row without invalidating other slots or routing.
!> \param snapshot row to clear
! **************************************************************************************************
   SUBROUTINE rskp_values_release(snapshot)

      TYPE(rskp_values_type), INTENT(INOUT)              :: snapshot

      IF (ALLOCATED(snapshot%packed)) DEALLOCATE (snapshot%packed)
      IF (ALLOCATED(snapshot%active)) DEALLOCATE (snapshot%active)
      IF (ALLOCATED(snapshot%fft_values)) DEALLOCATE (snapshot%fft_values)
      NULLIFY (snapshot%source)
      snapshot%matrix_row = 0
      snapshot%mirror_sign = 1
      snapshot%scale = z_one

   END SUBROUTINE rskp_values_release

! **************************************************************************************************
!> \brief Release backend-owned geometry and FFT work buffers from a route.
!> \param route route whose backend plan is no longer needed
! **************************************************************************************************
   SUBROUTINE rskp_route_release(route)
      TYPE(rskp_route_type), INTENT(INOUT)               :: route

      CALL cell_to_k_grid_release(route%transform%fft_work)
      IF (ASSOCIATED(route%matrix_struct)) CALL cp_fm_struct_release(route%matrix_struct)
      route = rskp_route_type()
   END SUBROUTINE rskp_route_release

! **************************************************************************************************
!> \brief FFT captured raw-element panels and retain the destination group's k points.
!>        Reuse the sparse route and the dense layout consumed by the batched solvers.
!> \param snapshot immutable packed raw-image panels
!> \param route fixed element/image layout and equivalent target tiles across groups
!> \param xkp full reciprocal grid
!> \param nkp_grid reciprocal-grid dimensions
!> \param kp_dist global k-point range of each group
!> \param inter communicator joining equivalent FM ranks across groups
!> \param batch complex local operators (row, column, local k point), including padding
!> \param max_bytes available workspace, excluding the caller-owned destination
! **************************************************************************************************
   SUBROUTINE cell_to_k_sparse_fft(snapshot, route, xkp, nkp_grid, kp_dist, inter, batch, max_bytes)
      TYPE(rskp_values_type), INTENT(IN), TARGET         :: snapshot
      TYPE(rskp_route_type), INTENT(INOUT)               :: route
      REAL(KIND=dp), INTENT(IN)                          :: xkp(:, :)
      INTEGER, INTENT(IN)                                :: nkp_grid(3), kp_dist(:, :)
      TYPE(mp_para_env_type), INTENT(IN)                 :: inter
      COMPLEX(KIND=dp), CONTIGUOUS, INTENT(OUT), TARGET  :: batch(:, :, :)
      INTEGER(KIND=int_8), INTENT(IN)                    :: max_bytes

      COMPLEX(KIND=dp), ALLOCATABLE                      :: received(:)
      COMPLEX(KIND=dp), ALLOCATABLE, TARGET              :: transformed(:, :, :)
      COMPLEX(KIND=dp), CONTIGUOUS, POINTER              :: output_values(:, :), send(:)
      INTEGER                                            :: element, first, group, i, igroup, ik, &
                                                            last, my_peer, nelem, nimages, nk, &
                                                            nlocal, npeers, ntile, nwave, offset, &
                                                            start, width
      INTEGER(KIND=int_8)                                :: available, dest, fixed_bytes, ngrid, &
                                                            per_element, workspace_bytes
      INTEGER, ALLOCATABLE                               :: ranges(:, :), rcount(:), rdisp(:), &
                                                            scount(:), sdisp(:)
      REAL(KIND=dp), CONTIGUOUS, POINTER                 :: values(:, :, :)

      CPASSERT(route%transform%backend == rskp_backend_fft)
      CPASSERT(SIZE(batch, 1) == route%matrix_struct%local_leading_dimension)
      output_values(1:SIZE(batch, 1, KIND=int_8)*SIZE(batch, 2), 1:SIZE(batch, 3)) => batch
      CPASSERT(route%transform%fft_work%regular)
      ngrid = PRODUCT(INT(route%transform%fft_work%nfft, int_8))
      nk = SIZE(xkp, 2)
      nlocal = SIZE(batch, 3)
      npeers = 1
      my_peer = 0
      IF (.NOT. route%world_target) THEN
         npeers = inter%num_pe
         my_peer = inter%mepos
      END IF
      ALLOCATE (ranges(2, npeers))
      ranges(:, 1) = [1, nk]
      IF (.NOT. route%world_target) ranges(:, :) = kp_dist
      CPASSERT(nlocal == ranges(2, my_peer + 1) - ranges(1, my_peer + 1) + 1)
      available = max_bytes
      IF (.NOT. route%world_target) CALL inter%min(available)
      ALLOCATE (scount(npeers), sdisp(npeers), rcount(npeers), rdisp(npeers))
      batch = z_zero
      DO igroup = 1, SIZE(route%block_first) - 1
         first = route%block_first(igroup)
         last = route%block_first(igroup + 1) - 1
         nelem = route%element_first(igroup + 1) - route%element_first(igroup)
         nimages = last - first + 1
         values(1:nelem, 1:1, 1:nimages) => &
            snapshot%packed(route%transform%packed_first(igroup):route%transform%packed_first(igroup + 1) - 1)
         IF (.NOT. ANY(values /= 0.0_dp)) CYCLE
         ! Each group transforms different AO elements; exchange only the owned k ranges.
         fixed_bytes = 16_int_8*nk + 28_int_8*SIZE(route%entry, 2) + 16_int_8*npeers
         per_element = 16_int_8*nk + 48_int_8*ngrid + &
                       16_int_8*npeers*MAXVAL(ranges(2, :) - ranges(1, :) + 1)
         width = INT(MAX(0_int_8, MIN(INT((nelem + npeers - 1)/npeers, int_8), &
                                      MAX(1_int_8, fft_batch_bytes/(48_int_8*ngrid)), &
                                      (available - fixed_bytes)/per_element, INT(HUGE(0), int_8)/nk, &
                                      INT(HUGE(0), int_8)/(INT(npeers, int_8)* &
                                                           MAXVAL(ranges(2, :) - ranges(1, :) + 1)))))
         IF (width == 0) THEN
            CPABORT("Forced lattice FFT block exceeds the memory budget")
         END IF
         DO element = 1, nelem, width*npeers
            start = element + my_peer*width
            ntile = MAX(0, MIN(width, nelem - start + 1))
            nwave = MIN(width*npeers, nelem - element + 1)
            ALLOCATE (transformed(ntile, 1, nk), received(nwave*nlocal))
            IF (ntile > 0) THEN
               workspace_bytes = available - fixed_bytes - 16_int_8*(SIZE(transformed, KIND=int_8) + &
                                                                     SIZE(received, KIND=int_8))
               CALL cell_to_k_grid_fft(values(start:start + ntile - 1, :, :), &
                                       route%entry(5:7, route%block_order(first:last)), &
                                       xkp, nkp_grid, transformed, &
                                       require_fft=.TRUE., work=route%transform%fft_work, cell_first=first, &
                                       max_workspace_bytes=workspace_bytes)
            END IF
            send(1:SIZE(transformed)) => transformed
            scount(:) = ntile*(ranges(2, :) - ranges(1, :) + 1)
            sdisp(:) = ntile*(ranges(1, :) - 1)
            offset = 0
            DO group = 1, npeers
               rcount(group) = MAX(0, MIN(width, nelem - element - (group - 1)*width + 1))*nlocal
               rdisp(group) = offset
               offset = offset + rcount(group)
            END DO
            IF (npeers > 1) THEN
               CALL inter%alltoall(send, scount, sdisp, received, rcount, rdisp)
            ELSE
               received(:) = send
            END IF
            DO group = 1, npeers
               ntile = rcount(group)/MAX(1, nlocal)
               DO ik = 1, nlocal
                  DO i = 1, ntile
                     start = route%element_first(igroup) + element + (group - 1)*width + i - 2
                     dest = route%element_offset(1, start)
                     IF (dest > 0) output_values(dest, ik) = snapshot%scale*received(rdisp(group) + (ik - 1)*ntile + i)
                     dest = route%element_offset(2, start)
                     IF (dest > 0) output_values(dest, ik) = snapshot%scale*REAL(snapshot%mirror_sign, dp)* &
                                                             CONJG(received(rdisp(group) + (ik - 1)*ntile + i))
                  END DO
               END DO
            END DO
            DEALLOCATE (transformed, received)
         END DO
      END DO
   END SUBROUTINE cell_to_k_sparse_fft

! **************************************************************************************************
!> \brief Fix parity panels per canonical atom pair and cache their shared phases.
!>        Reserve one complex worker and caller-specified retained output storage.
!>        Oversized tables use the same contraction kernel with tile-local phases.
!> \param plan geometry-only dictionary and optional immutable phase table
!> \param translations signed image translations
!> \param group_first boundaries of canonical AO groups in order
!> \param order image permutation grouped by canonical AO pair
!> \param element_first cumulative row counts of the packed group panels
!> \param xkp fixed group-local reciprocal coordinates
!> \param max_bytes shared-phase and worker/output-workspace budget
!> \param output_bytes workspace reserved for retained output tiles
! **************************************************************************************************
   SUBROUTINE cell_to_k_phase_prepare(plan, translations, group_first, order, element_first, xkp, max_bytes, output_bytes)
      TYPE(rskp_backend_plan_type), INTENT(INOUT)        :: plan
      INTEGER, INTENT(IN)                                :: translations(:, :), group_first(:), &
                                                            order(:), element_first(:)
      REAL(KIND=dp), INTENT(IN)                          :: xkp(:, :)
      INTEGER(KIND=int_8), INTENT(IN)                    :: max_bytes, output_bytes

      INTEGER                                            :: c, first, g, i, iblock, ik, j, m, n, &
                                                            nfold, nk, npair, nr, nsingle, pos
      INTEGER(KIND=int_8)                                :: cache_bytes, minimum_bytes
      INTEGER, ALLOCATABLE                               :: cells(:, :), classes(:), &
                                                            global_class(:), local_class(:), &
                                                            slots(:)
      LOGICAL, ALLOCATABLE                               :: members(:, :)
      REAL(KIND=dp)                                      :: arg

      CALL rskp_image_classes(translations, classes, plan%phase_cells)
      n = SIZE(group_first) - 1
      ALLOCATE (plan%phase_count(3, n), plan%phase_first(n + 1), plan%packed_first(n + 1), &
                plan%phase_index(SIZE(classes)), plan%fold_map(2, SIZE(classes)))
      plan%fold_map = 0
      plan%phase_first(1) = 1
      plan%packed_first(1) = 1
      plan%max_per_k = 0_int_8
      DO g = 1, n
         first = group_first(g)
         CALL rskp_image_classes(translations(:, order(first:group_first(g + 1) - 1)), &
                                 local_class, cells)
         ALLOCATE (members(2, SIZE(cells, 2)), slots(SIZE(cells, 2)), global_class(SIZE(cells, 2)))
         members = .FALSE.
         DO i = 1, SIZE(local_class)
            c = local_class(i)
            IF (c == 0) CYCLE
            members(MERGE(1, 2, c > 0), ABS(c)) = .TRUE.
            iblock = order(first + i - 1)
            global_class(ABS(c)) = ABS(classes(iblock))
         END DO
         npair = COUNT(ALL(members, DIM=1))
         nsingle = SIZE(cells, 2) - npair
         plan%phase_count(:, g) = [npair, nsingle, MERGE(1, 0, ANY(local_class == 0))]
         i = 0
         j = 0
         DO c = 1, SIZE(cells, 2)
            IF (ALL(members(:, c))) THEN
               i = i + 1
               slots(c) = i
               plan%phase_index(plan%phase_first(g) + i - 1) = global_class(c)
            ELSE
               j = j + 1
               slots(c) = 2*npair + j
               plan%phase_index(plan%phase_first(g) + npair + j - 1) = MERGE(1, -1, members(1, c))*global_class(c)
            END IF
         END DO
         m = element_first(g + 1) - element_first(g)
         nfold = 2*npair + nsingle + plan%phase_count(3, g)
         plan%phase_first(g + 1) = plan%phase_first(g) + npair + nsingle
         IF (INT(plan%packed_first(g), int_8) + INT(m, int_8)*nfold > HUGE(m)) THEN
            CPABORT("R-to-K folded panel offsets exceed the integer range")
         END IF
         plan%packed_first(g + 1) = plan%packed_first(g) + m*nfold
         plan%max_per_k = MAX(plan%max_per_k, INT(m, int_8) + MAX(npair, nsingle))
         DO i = 1, SIZE(local_class)
            iblock = order(first + i - 1)
            c = local_class(i)
            pos = nfold
            IF (c /= 0) pos = slots(ABS(c))
            plan%fold_map(1, iblock) = plan%packed_first(g) + (pos - 1)*m
            IF (c == 0) CYCLE
            IF (ALL(members(:, ABS(c)))) THEN
               plan%fold_map(2, iblock) = SIGN(plan%packed_first(g) + (npair + pos - 1)*m, c)
            END IF
         END DO
         DEALLOCATE (members, slots, global_class)
      END DO

      nr = SIZE(plan%phase_cells, 2)
      nk = SIZE(xkp, 2)
      cache_bytes = 16_int_8*nr*nk
      minimum_bytes = 16_int_8*plan%max_per_k + output_bytes
      IF (cache_bytes + minimum_bytes > max_bytes) RETURN
      ALLOCATE (plan%phases(nr, nk, 2))
      DO ik = 1, nk
!$OMP SIMD PRIVATE(arg)
         DO i = 1, nr
            arg = twopi*DOT_PRODUCT(REAL(plan%phase_cells(:, i), dp), xkp(:, ik))
            plan%phases(i, ik, 1) = COS(arg)
            plan%phases(i, ik, 2) = SIN(arg)
         END DO
      END DO
   END SUBROUTINE cell_to_k_phase_prepare

! **************************************************************************************************
!> \brief Prepare the selected R-to-K backend plan once for a fixed route.
!>        Both backends retain their geometry with the route; the snapshots then
!>        provide only values and use the common evaluate interface.
!> \param route fixed sparse AO route and backend workspace
!> \param xkp reciprocal coordinates for the route's k-point domain
!> \param nkp_grid regular-grid dimensions required by the FFT backend
!> \param max_bytes phase-cache budget for the GEMM backend
! **************************************************************************************************
   SUBROUTINE cell_to_k_backend_prepare(route, xkp, nkp_grid, max_bytes)
      TYPE(rskp_route_type), INTENT(INOUT)               :: route
      REAL(KIND=dp), INTENT(IN)                          :: xkp(:, :)
      INTEGER, INTENT(IN)                                :: nkp_grid(3)
      INTEGER(KIND=int_8), INTENT(IN)                    :: max_bytes

      SELECT CASE (route%transform%backend)
      CASE (rskp_backend_gemm)
         CALL cell_to_k_phase_prepare(route%transform, route%entry(5:7, :), route%block_first, route%block_order, &
                                      route%element_first, xkp, max_bytes, &
                                      32_int_8*route%matrix_struct%local_leading_dimension*MAX(1, route%ncol))
      CASE (rskp_backend_fft)
         CALL cell_to_k_grid_prepare(route%transform%fft_work, route%entry(5:7, route%block_order), xkp, nkp_grid)
         IF (.NOT. route%transform%fft_work%regular) THEN
            CPABORT("Forced lattice FFT needs a supported complete regular grid")
         END IF
      CASE DEFAULT
         CPABORT("Unknown R-to-K backend")
      END SELECT
   END SUBROUTINE cell_to_k_backend_prepare

! **************************************************************************************************
!> \brief Contract parity-folded pairs and concatenate phases for unpaired translations.
!> \param packed immutable parity panels of one AO group
!> \param plan phase dictionary shared by the FM and DBCSR adapters
!> \param igroup canonical atom-pair group
!> \param xkp k coordinates for this tile
!> \param first_k first index in the route-local phase table
!> \param factors cosine columns followed by optional sine columns
!> \param RESULT real columns followed by optional imaginary columns
!> \param active nonzero even, odd, single and origin panels, determined at capture
! **************************************************************************************************
   SUBROUTINE cell_to_k_block_gemm(packed, plan, igroup, xkp, first_k, factors, RESULT, active)
      REAL(KIND=dp), CONTIGUOUS, INTENT(IN), TARGET      :: packed(:)
      TYPE(rskp_backend_plan_type), INTENT(IN)           :: plan
      INTEGER, INTENT(IN)                                :: igroup
      REAL(KIND=dp), INTENT(IN)                          :: xkp(:, :)
      INTEGER, INTENT(IN)                                :: first_k
      REAL(KIND=dp), CONTIGUOUS, INTENT(OUT)             :: factors(:, :), RESULT(:, :)
      LOGICAL, INTENT(IN)                                :: active(4)

      INTEGER                                            :: first, i, ik, nactive, nelem, nk, npair, &
                                                            nsingle, pass, phase, start
      LOGICAL                                            :: use_cosine, use_sine
      REAL(KIND=dp)                                      :: arg, beta(2)
      REAL(KIND=dp), CONTIGUOUS, POINTER                 :: values(:, :), values_im(:, :)

      nelem = SIZE(RESULT, 1)
      nk = SIZE(xkp, 2)
      CPASSERT(SIZE(RESULT, 2) == nk .OR. SIZE(RESULT, 2) == 2*nk)
      npair = plan%phase_count(1, igroup)
      nsingle = plan%phase_count(2, igroup)
      first = plan%phase_first(igroup)
      start = 1
      IF (ALLOCATED(plan%phases)) THEN
         CPASSERT(first_k >= 1 .AND. first_k + nk - 1 <= SIZE(plan%phases, 2))
      END IF
      beta = 0.0_dp
      DO pass = 1, 2
         nactive = MERGE(npair, nsingle, pass == 1)
         IF (nactive == 0) CYCLE
         use_cosine = MERGE(active(1), active(3), pass == 1)
         use_sine = SIZE(RESULT, 2) == 2*nk .AND. MERGE(active(2), active(3), pass == 1)
         IF (.NOT. (use_cosine .OR. use_sine)) THEN
            start = start + nelem*nactive*MERGE(2, 1, pass == 1)
            first = first + nactive
            CYCLE
         END IF
         DO ik = 1, nk
            DO i = 1, nactive
               phase = plan%phase_index(first + i - 1)
               IF (ALLOCATED(plan%phases)) THEN
                  IF (use_cosine) factors(i, ik) = plan%phases(ABS(phase), first_k + ik - 1, 1)
                  IF (use_sine) THEN
                     factors(i, nk + ik) = REAL(SIGN(1, phase), dp)*plan%phases(ABS(phase), first_k + ik - 1, 2)
                  END IF
               ELSE
                  arg = twopi*DOT_PRODUCT(REAL(plan%phase_cells(:, ABS(phase)), dp), xkp(:, ik))
                  IF (use_cosine) factors(i, ik) = COS(arg)
                  IF (use_sine) factors(i, nk + ik) = REAL(SIGN(1, phase), dp)*SIN(arg)
               END IF
            END DO
         END DO
         values(1:nelem, 1:nactive) => packed(start:start + nelem*nactive - 1)
         IF (pass == 1) THEN
            start = start + nelem*npair
            IF (active(1)) THEN
               CALL rskp_gemm_contract(nelem, npair, nk, values, factors(:, 1:nk), RESULT(:, 1:nk), 0.0_dp)
               beta(1) = 1.0_dp
            END IF
            IF (SIZE(RESULT, 2) == 2*nk .AND. active(2)) THEN
               values_im(1:nelem, 1:npair) => packed(start:start + nelem*npair - 1)
               CALL rskp_gemm_contract(nelem, npair, nk, values_im, factors(:, nk + 1:2*nk), &
                                       RESULT(:, nk + 1:2*nk), 0.0_dp)
               beta(2) = 1.0_dp
            END IF
            start = start + nelem*npair
         ELSE
            IF (SIZE(RESULT, 2) == 2*nk .AND. beta(1) == beta(2)) THEN
               CALL rskp_gemm_contract(nelem, nsingle, 2*nk, values, factors, RESULT, beta(1))
            ELSE
               CALL rskp_gemm_contract(nelem, nsingle, nk, values, factors(:, 1:nk), RESULT(:, 1:nk), beta(1))
               IF (SIZE(RESULT, 2) == 2*nk) THEN
                  CALL rskp_gemm_contract(nelem, nsingle, nk, values, factors(:, nk + 1:2*nk), &
                                          RESULT(:, nk + 1:2*nk), beta(2))
               END IF
            END IF
            beta = 1.0_dp
            start = start + nelem*nsingle
         END IF
         first = first + nactive
      END DO
      IF (beta(1) == 0.0_dp) RESULT(:, 1:nk) = 0.0_dp
      IF (SIZE(RESULT, 2) == 2*nk .AND. beta(2) == 0.0_dp) RESULT(:, nk + 1:2*nk) = 0.0_dp
      IF (active(4)) THEN
         DO ik = 1, nk
            RESULT(:, ik) = RESULT(:, ik) + packed(start:start + nelem - 1)
         END DO
      END IF
   END SUBROUTINE cell_to_k_block_gemm

! **************************************************************************************************
!> \brief Bound worker scratch and output tiles for either prepared backend.
!>        Packed input belongs to the snapshot and is shared by every k tile.
!> \param route borrowed local AO-element layout
!> \param nk local k-point count
!> \param ncomp contraction components for split-FM outputs (two for a complex source scale)
!> \param nworkers requested workers, reduced to fit the workspace
!> \param nmatrices retained complex operator tiles; zero writes existing split-FM outputs
!> \param max_bytes shared-phase, worker-scratch and output-tile budget
!> \param width k points per worker
!> \param work_size real scratch elements per worker
! **************************************************************************************************
   SUBROUTINE rskp_batch_size(route, nk, ncomp, nworkers, nmatrices, max_bytes, width, work_size)
      TYPE(rskp_route_type), INTENT(IN)                  :: route
      INTEGER, INTENT(IN)                                :: nk, ncomp
      INTEGER, INTENT(INOUT)                             :: nworkers
      INTEGER, INTENT(IN)                                :: nmatrices
      INTEGER(KIND=int_8), INTENT(IN)                    :: max_bytes
      INTEGER, INTENT(OUT)                               :: width
      INTEGER(KIND=int_8), INTENT(OUT)                   :: work_size

      INTEGER                                            :: nc
      INTEGER(KIND=int_8)                                :: available_bytes, output_bytes, &
                                                            per_k_bytes

      CPASSERT(nworkers > 0 .AND. nk > 0 .AND. nmatrices >= 0)
      available_bytes = max_bytes
      IF (ALLOCATED(route%transform%phases)) available_bytes = available_bytes - 8_int_8*SIZE(route%transform%phases, KIND=int_8)
      output_bytes = 16_int_8*nmatrices*route%matrix_struct%local_leading_dimension*MAX(1, route%ncol)
      nc = 2
      IF (nmatrices == 0) nc = ncomp
      per_k_bytes = MAX(1_int_8, 8_int_8*nc*route%transform%max_per_k + output_bytes)
      nworkers = MIN(nworkers, nk)
      nworkers = MIN(nworkers, INT(MAX(1_int_8, available_bytes/per_k_bytes)))
      width = INT(MIN(INT(nk, int_8), MAX(1_int_8, available_bytes/(nworkers*per_k_bytes))))
      width = MIN(width, (nk + nworkers - 1)/nworkers)
      work_size = INT(nc, int_8)*width*route%transform%max_per_k
   END SUBROUTINE rskp_batch_size

! **************************************************************************************************
!> \brief Assemble a local complex operator tile directly from packed sparse inputs.
!> \param snapshot immutable BLAS-layout input values
!> \param route borrowed local AO-element layout
!> \param xkp coordinates for this k tile
!> \param first_k first index in the route-local phase table
!> \param work caller-owned real phase and result scratch
!> \param batch optional complex output (row, column, k point), including padding
!> \param output optional existing split-FM destinations (k point, component)
! **************************************************************************************************
   SUBROUTINE cell_to_k_gemm_batch(snapshot, route, xkp, first_k, work, batch, output)
      TYPE(rskp_values_type), INTENT(IN), TARGET         :: snapshot
      TYPE(rskp_route_type), INTENT(IN)                  :: route
      REAL(KIND=dp), INTENT(IN)                          :: xkp(:, :)
      INTEGER, INTENT(IN)                                :: first_k
      REAL(KIND=dp), CONTIGUOUS, INTENT(INOUT), TARGET   :: work(:)
      COMPLEX(KIND=dp), CONTIGUOUS, INTENT(OUT), &
         OPTIONAL, TARGET                                :: batch(:, :, :)
      TYPE(cp_fm_p_type), INTENT(IN), OPTIONAL           :: output(:, :)

      COMPLEX(KIND=dp)                                   :: value
      COMPLEX(KIND=dp), CONTIGUOUS, POINTER              :: complex_values(:, :)
      INTEGER                                            :: element, first, ic, igroup, j, ld, nc, &
                                                            nelem, nphase, width
      INTEGER(KIND=int_8)                                :: dest, nfactor
      REAL(KIND=dp), CONTIGUOUS, POINTER                 :: factors(:, :), real_values(:), &
                                                            RESULT(:, :)

      CPASSERT(route%transform%backend == rskp_backend_gemm)
      CPASSERT(PRESENT(batch) .NEQV. PRESENT(output))
      width = SIZE(xkp, 2)
      ld = route%matrix_struct%local_leading_dimension
      nc = 2
      IF (PRESENT(batch)) THEN
         CPASSERT(SIZE(batch, 1) == route%matrix_struct%local_leading_dimension)
         CPASSERT(SIZE(batch, 2) >= route%ncol .AND. SIZE(batch, 3) == width)
         batch = z_zero
         complex_values(1:SIZE(batch, 1, KIND=int_8)*SIZE(batch, 2), 1:width) => batch
      ELSE
         nc = SIZE(output, 2)
         CPASSERT(SIZE(output, 1) == width .AND. (nc == 1 .OR. nc == 2))
         DO j = 1, width
            DO ic = 1, nc
               CPASSERT(SIZE(output(j, ic)%matrix%local_data, 1) == ld)
               output(j, ic)%matrix%local_data = 0.0_dp
            END DO
         END DO
         IF (AIMAG(snapshot%scale) /= 0.0_dp) nc = 2
      END IF
      DO igroup = 1, SIZE(route%transform%phase_count, 2)
         first = route%element_first(igroup)
         nelem = route%element_first(igroup + 1) - first
         IF (nelem == 0) CYCLE
         IF (.NOT. ANY(snapshot%active(:, igroup))) CYCLE
         nphase = MAXVAL(route%transform%phase_count(1:2, igroup))
         nfactor = INT(nc, int_8)*nphase*width
         CPASSERT(SIZE(work, KIND=int_8) >= nfactor + INT(nc, int_8)*nelem*width)
         factors(1:nphase, 1:nc*width) => work(1:nfactor)
         RESULT(1:nelem, 1:nc*width) => work(nfactor + 1:nfactor + INT(nc, int_8)*nelem*width)
         CALL cell_to_k_block_gemm(snapshot%packed(route%transform%packed_first(igroup): &
                                                   route%transform%packed_first(igroup + 1) - 1), &
                                   route%transform, igroup, xkp, first_k, factors, RESULT, snapshot%active(:, igroup))
         DO j = 1, width
            IF (PRESENT(batch)) THEN
!$OMP SIMD PRIVATE(dest,value)
               DO element = first, route%element_first(igroup + 1) - 1
                  value = CMPLX(RESULT(element - first + 1, j), RESULT(element - first + 1, width + j), KIND=dp)
                  dest = route%element_offset(1, element)
                  IF (dest > 0) complex_values(dest, j) = snapshot%scale*value
                  dest = route%element_offset(2, element)
                  IF (dest > 0) complex_values(dest, j) = snapshot%scale*REAL(snapshot%mirror_sign, dp)*CONJG(value)
               END DO
            ELSE
               DO ic = 1, SIZE(output, 2)
                  real_values(1:SIZE(output(j, ic)%matrix%local_data, KIND=int_8)) => output(j, ic)%matrix%local_data
!$OMP SIMD PRIVATE(dest,value)
                  DO element = first, route%element_first(igroup + 1) - 1
                     dest = route%element_offset(1, element)
                     IF (AIMAG(snapshot%scale) == 0.0_dp) THEN
                        IF (dest > 0) real_values(dest) = REAL(snapshot%scale, dp)* &
                                                          RESULT(element - first + 1, j + (ic - 1)*width)
                     ELSE
                        value = snapshot%scale*CMPLX(RESULT(element - first + 1, j), &
                                                     RESULT(element - first + 1, width + j), dp)
                        IF (dest > 0) real_values(dest) = MERGE(REAL(value, dp), AIMAG(value), ic == 1)
                     END IF
                     dest = route%element_offset(2, element)
                     IF (AIMAG(snapshot%scale) == 0.0_dp) THEN
                        IF (dest > 0) real_values(dest) = REAL(snapshot%scale, dp)*REAL(snapshot%mirror_sign, dp)* &
                                                          MERGE(1.0_dp, -1.0_dp, ic == 1)* &
                                                          RESULT(element - first + 1, j + (ic - 1)*width)
                     ELSE
                        value = snapshot%scale*REAL(snapshot%mirror_sign, dp)* &
                                CMPLX(RESULT(element - first + 1, j), -RESULT(element - first + 1, width + j), dp)
                        IF (dest > 0) real_values(dest) = MERGE(REAL(value, dp), AIMAG(value), ic == 1)
                     END IF
                  END DO
               END DO
            END IF
         END DO
      END DO
   END SUBROUTINE cell_to_k_gemm_batch

! **************************************************************************************************
!> \brief Evaluate a k-point tile through the selected prepared R-to-K backend.
!>        GEMM contracts the retained real-space panel; FFT returns a slice of
!>        the reciprocal-space snapshot materialized during capture.
!> \param snapshot backend input or prepared reciprocal-space values
!> \param route fixed sparse layout and backend plan
!> \param xkp reciprocal coordinates for this k tile
!> \param first_local first k index in the route-local domain
!> \param work caller-owned GEMM scratch; unused by the FFT backend
!> \param batch complex output tile (row, column, k point)
!> \param output optional split-FM output tile (k point, component), instead of batch
! **************************************************************************************************
   SUBROUTINE cell_to_k_backend_evaluate_batch(snapshot, route, xkp, first_local, work, batch, output)
      TYPE(rskp_values_type), INTENT(IN), TARGET         :: snapshot
      TYPE(rskp_route_type), INTENT(IN)                  :: route
      REAL(KIND=dp), INTENT(IN)                          :: xkp(:, :)
      INTEGER, INTENT(IN)                                :: first_local
      REAL(KIND=dp), CONTIGUOUS, INTENT(INOUT), TARGET   :: work(:)
      COMPLEX(KIND=dp), CONTIGUOUS, INTENT(OUT), &
         OPTIONAL                                        :: batch(:, :, :)
      TYPE(cp_fm_p_type), INTENT(IN), OPTIONAL           :: output(:, :)

      INTEGER                                            :: ic, j, width
      LOGICAL                                            :: compatible

      CPASSERT(PRESENT(batch) .NEQV. PRESENT(output))
      width = SIZE(xkp, 2)
      IF (PRESENT(output)) THEN
         CPASSERT(SIZE(output, 1) == width)
         CPASSERT(SIZE(output, 2) >= 1 .AND. SIZE(output, 2) <= 2)
         DO j = 1, width
            DO ic = 1, SIZE(output, 2)
               compatible = cp_fm_struct_equivalent(output(j, ic)%matrix%matrix_struct, route%matrix_struct)
               CPASSERT(compatible)
            END DO
         END DO
      END IF
      SELECT CASE (route%transform%backend)
      CASE (rskp_backend_gemm)
         CALL cell_to_k_gemm_batch(snapshot, route, xkp, first_local, work, batch=batch, output=output)
      CASE (rskp_backend_fft)
         CPASSERT(ALLOCATED(snapshot%fft_values))
         CPASSERT(first_local >= 1)
         CPASSERT(first_local + width - 1 <= SIZE(snapshot%fft_values, 3))
         IF (PRESENT(batch)) THEN
            CPASSERT(SIZE(batch, 1) == SIZE(snapshot%fft_values, 1))
            CPASSERT(SIZE(batch, 2) == SIZE(snapshot%fft_values, 2))
            CPASSERT(SIZE(batch, 3) == width)
            batch = snapshot%fft_values(:, :, first_local:first_local + width - 1)
         ELSE
            DO j = 1, width
               output(j, 1)%matrix%local_data = REAL(snapshot%fft_values(:, :, first_local + j - 1), dp)
               IF (SIZE(output, 2) == 2) THEN
                  output(j, 2)%matrix%local_data = AIMAG(snapshot%fft_values(:, :, first_local + j - 1))
               END IF
            END DO
         END IF
      CASE DEFAULT
         CPABORT("Unknown R-to-K backend")
      END SELECT
   END SUBROUTINE cell_to_k_backend_evaluate_batch

! **************************************************************************************************
!> \brief Weighted K-to-R contraction of one AO block, before spatial rotations and MPI.
!>        The caller supplies the atom mapping and image multiplicities. Phases include
!>        the gauge shift, time reversal and conjugation of a transposed source block.
!> \param source borrowed local matrices by k point, real/imaginary component and spin
!> \param tiles local FM intervals (first/last row and column) by atom, k, component and spin
!> \param ikp_index local source indices for this symmetry contribution
!> \param xkp reciprocal coordinates for these contributions
!> \param weight original k-point weights, including symmetry multiplicity
!> \param images translations, real multiplicities and imaginary sign sums
!> \param source_pair source AO atom pair
!> \param block_shape source AO row and column sizes
!> \param block_offset global first row and column of the source block
!> \param shift column-minus-row gauge shift
!> \param time_reversal conjugate the source before applying its gauge phase
!> \param trans conjugate the imaginary contribution on transposition
!> \param packed reusable thread-local source workspace
!> \param phases reusable thread-local phase workspace
!> \param sums allocated flattened source blocks by image and spin
!> \param gemm_ctx thread-owned host GEMM context
! **************************************************************************************************
   SUBROUTINE k_grid_to_cell_contract(source, tiles, ikp_index, xkp, weight, images, source_pair, &
                                      block_shape, block_offset, shift, time_reversal, trans, &
                                      packed, phases, sums, gemm_ctx)
      TYPE(cp_fm_p_type), INTENT(IN)                     :: source(:, :, :)
      INTEGER, INTENT(IN)                                :: tiles(:, :, :, :, :), ikp_index(:)
      REAL(KIND=dp), INTENT(IN)                          :: xkp(:, :), weight(:)
      INTEGER, INTENT(IN)                                :: images(:, :), source_pair(2), &
                                                            block_shape(2), block_offset(2), &
                                                            shift(3)
      LOGICAL, INTENT(IN)                                :: time_reversal, trans
      REAL(KIND=dp), ALLOCATABLE, INTENT(INOUT)          :: packed(:, :), phases(:, :), sums(:, :, :)
      TYPE(local_gemm_ctxt_type), INTENT(INOUT)          :: gemm_ctx

      INTEGER(KIND=int_8), PARAMETER                     :: max_elements = 65536_int_8

      INTEGER                                            :: a, b, batch, column, first, i, ic, &
                                                            ifirst, ii, ikp, ilast, image, ispin, &
                                                            j, jj, k, last, nc, nclass, nelement, &
                                                            nimage, nrow, nsine, nspin
      INTEGER, ALLOCATABLE                               :: cells(:, :), image_class(:), &
                                                            sine_index(:)
      LOGICAL                                            :: contributes, has_zero, identity_gauge, &
                                                            real_only
      REAL(KIND=dp)                                      :: arg, sine
      REAL(KIND=dp), ALLOCATABLE                         :: beta(:), gauge(:), &
                                                            gauge_factors(:, :, :), &
                                                            reduced(:, :, :), zero_sum(:, :)
      TYPE(cp_fm_type), POINTER                          :: matrix

      nc = SIZE(source, 2)
      nspin = SIZE(source, 3)
      real_only = nc == 1
      CPASSERT(nc == 1 .OR. nc == 2)
      CPASSERT(SIZE(xkp, 2) == SIZE(ikp_index) .AND. SIZE(weight) == SIZE(ikp_index))
      a = source_pair(1)
      b = source_pair(2)
      nelement = block_shape(1)*block_shape(2)
      nimage = SIZE(images, 2)
      IF (ALLOCATED(sums)) THEN
         IF (ANY(SHAPE(sums) /= [nelement, nimage, nspin])) DEALLOCATE (sums)
      END IF
      IF (.NOT. ALLOCATED(sums)) ALLOCATE (sums(nelement, nimage, nspin))
      IF (nelement == 0 .OR. nimage == 0 .OR. SIZE(ikp_index) == 0) THEN
         sums = 0.0_dp
         RETURN
      END IF
      CALL rskp_image_classes(images(1:3, :), image_class, cells)
      nclass = SIZE(cells, 2)
      has_zero = ANY(image_class == 0)
      ! Only classes whose output multiplicities use the imaginary source need sine.
      ALLOCATE (sine_index(nclass), source=0)
      IF (.NOT. real_only) THEN
         DO image = 1, nimage
            i = ABS(image_class(image))
            IF (i > 0 .AND. images(5, image) /= 0) sine_index(i) = 1
         END DO
      END IF
      nsine = 0
      DO i = 1, nclass
         IF (sine_index(i) == 0) CYCLE
         nsine = nsine + 1
         sine_index(i) = nsine
      END DO
      identity_gauge = ALL(shift == 0) .AND. .NOT. time_reversal
      IF (.NOT. identity_gauge .AND. nc == 2) ALLOCATE (gauge(nelement))
      ALLOCATE (reduced(nelement, nclass + nsine, nspin))
      ALLOCATE (beta(nspin), source=0.0_dp)
      IF (has_zero) ALLOCATE (zero_sum(nelement, nspin), source=0.0_dp)
      ! Bound packed input and the shared phase tile; all spins reuse the same phases.
      batch = MAX(1, MIN(128, SIZE(ikp_index), &
                         INT(max_elements/(INT(nc, int_8)*nelement + nclass + nsine))))
      CALL rskp_work_matrix(packed, nelement, nc*batch)
      CALL rskp_work_matrix(phases, batch, nclass + nsine)
      IF (.NOT. identity_gauge) ALLOCATE (gauge_factors(2, 2, batch))
      DO first = 1, SIZE(ikp_index), batch
         last = MIN(first + batch - 1, SIZE(ikp_index))
         DO image = 1, nclass
            DO k = first, last
               arg = twopi*DOT_PRODUCT(REAL(cells(:, image), dp), xkp(:, k))
               phases(k - first + 1, image) = weight(k)*COS(arg)
               IF (sine_index(image) > 0) THEN
                  phases(k - first + 1, nclass + sine_index(image)) = weight(k)*SIN(arg)
               END IF
            END DO
         END DO
         IF (.NOT. identity_gauge) THEN
            DO k = first, last
               CALL rskp_phase_factors(shift, xkp(:, k), time_reversal, real_only, &
                                       gauge_factors(:, :, k - first + 1))
            END DO
         END IF
         DO ispin = 1, nspin
            packed = 0.0_dp
            contributes = .FALSE.
            DO k = first, last
               ikp = ikp_index(k)
               DO ic = 1, nc
                  ifirst = tiles(1, a, ikp, ic, ispin)
                  ilast = tiles(2, a, ikp, ic, ispin)
                  IF (ifirst > ilast .OR. tiles(3, b, ikp, ic, ispin) > tiles(4, b, ikp, ic, ispin)) CYCLE
                  contributes = .TRUE.
                  column = k - first + 1 + (ic - 1)*batch
                  matrix => source(ikp, ic, ispin)%matrix
                  nrow = ilast - ifirst + 1
                  DO j = tiles(3, b, ikp, ic, ispin), tiles(4, b, ikp, ic, ispin)
                     jj = matrix%matrix_struct%col_indices(j) - block_offset(2)
                     ii = matrix%matrix_struct%row_indices(ifirst) - block_offset(1) + 1 + jj*block_shape(1)
                     IF (matrix%matrix_struct%row_indices(ilast) - &
                         matrix%matrix_struct%row_indices(ifirst) == nrow - 1) THEN
                        packed(ii:ii + nrow - 1, column) = matrix%local_data(ifirst:ilast, j)
                     ELSE
                        DO i = ifirst, ilast
                           ii = matrix%matrix_struct%row_indices(i) - block_offset(1) + 1 + jj*block_shape(1)
                           packed(ii, column) = matrix%local_data(i, j)
                        END DO
                     END IF
                  END DO
               END DO
               column = k - first + 1
               IF (.NOT. identity_gauge) THEN
                  IF (nc == 1) THEN
                     packed(:, column) = gauge_factors(1, 1, column)*packed(:, column)
                  ELSE
                     gauge(:) = gauge_factors(1, 1, column)*packed(:, column) + &
                                gauge_factors(1, 2, column)*packed(:, batch + column)
                     packed(:, batch + column) = gauge_factors(2, 1, column)*packed(:, column) + &
                                                 gauge_factors(2, 2, column)*packed(:, batch + column)
                     packed(:, column) = gauge
                  END IF
               END IF
               IF (has_zero) zero_sum(:, ispin) = zero_sum(:, ispin) + weight(k)*packed(:, column)
            END DO
            IF (.NOT. contributes .OR. nclass == 0) CYCLE
            CALL gemm_ctx%gemm('N', 'N', nelement, nclass, last - first + 1, 1.0_dp, &
                               packed, nelement, phases, batch, beta(ispin), reduced(:, :, ispin), nelement)
            IF (nsine > 0) THEN
               CALL gemm_ctx%gemm('N', 'N', nelement, nsine, last - first + 1, 1.0_dp, &
                                  packed(:, batch + 1:), nelement, phases(:, nclass + 1:), batch, &
                                  beta(ispin), reduced(:, nclass + 1:, ispin), nelement)
            END IF
            beta(ispin) = 1.0_dp
         END DO
      END DO
      ! Expand the shared parity sums with each image's original multiplicities.
      DO ispin = 1, nspin
         IF (beta(ispin) == 0.0_dp) reduced(:, :, ispin) = 0.0_dp
         DO image = 1, nimage
            i = ABS(image_class(image))
            IF (i == 0) THEN
               sums(:, image, ispin) = REAL(images(4, image), dp)*zero_sum(:, ispin)
            ELSE
               sums(:, image, ispin) = REAL(images(4, image), dp)*reduced(:, i, ispin)
               IF (sine_index(i) == 0) CYCLE
               sine = REAL(images(5, image)*SIGN(1, image_class(image)), dp)*MERGE(-1.0_dp, 1.0_dp, trans)
               sums(:, image, ispin) = sums(:, image, ispin) + sine*reduced(:, nclass + sine_index(i), ispin)
            END IF
         END DO
      END DO
   END SUBROUTINE k_grid_to_cell_contract

! **************************************************************************************************
!> \brief Classify translations by +/-R, without assuming operator or k-point symmetry.
!> \param cells integer translations
!> \param classes signed class for each input; zero denotes the origin
!> \param representatives canonical nonzero translations, one per class
! **************************************************************************************************
   SUBROUTINE rskp_image_classes(cells, classes, representatives)
      INTEGER, INTENT(IN)                                :: cells(:, :)
      INTEGER, ALLOCATABLE, INTENT(OUT)                  :: classes(:), representatives(:, :)

      INTEGER                                            :: d, i, n, nclass
      INTEGER, ALLOCATABLE                               :: first(:), keys(:, :), order(:), sense(:)

      n = SIZE(cells, 2)
      ALLOCATE (keys(3, n), order(n), sense(n), first(n), classes(n))
      DO i = 1, n
         sense(i) = 1
         DO d = 1, 3
            IF (cells(d, i) == 0) CYCLE
            sense(i) = SIGN(1, cells(d, i))
            EXIT
         END DO
         keys(:, i) = sense(i)*cells(:, i)
      END DO
      IF (n > 0) CALL sort(keys, 1, n, 1, 3, order)
      nclass = 0
      DO i = 1, n
         IF (ALL(keys(:, i) == 0)) THEN
            classes(order(i)) = 0
            CYCLE
         END IF
         IF (nclass == 0) THEN
            nclass = 1
            first(nclass) = i
         ELSE IF (ANY(keys(:, i) /= keys(:, first(nclass)))) THEN
            nclass = nclass + 1
            first(nclass) = i
         END IF
         classes(order(i)) = sense(order(i))*nclass
      END DO
      ALLOCATE (representatives(3, nclass))
      representatives(:, :) = keys(:, first(1:nclass))
   END SUBROUTINE rskp_image_classes

! **************************************************************************************************
!> \brief Real two-component representation of the Bloch phase and time reversal.
!> \param shift column-minus-row cell shift with the selected phase direction
!> \param xkp target k-point coordinates
!> \param time_reversal conjugate the source before applying the phase
!> \param real_only reject phases not representable by real wavefunctions
!> \param phase coefficients mapping source [Re,Im] to transformed [Re,Im]
! **************************************************************************************************
   SUBROUTINE rskp_phase_factors(shift, xkp, time_reversal, real_only, phase)
      INTEGER, INTENT(IN)                                :: shift(3)
      REAL(KIND=dp), INTENT(IN)                          :: xkp(3)
      LOGICAL, INTENT(IN)                                :: time_reversal, real_only
      REAL(KIND=dp), INTENT(OUT)                         :: phase(2, 2)

      REAL(KIND=dp)                                      :: arg, cosine, sine

      arg = twopi*(REAL(shift(1), dp)*xkp(1) + REAL(shift(2), dp)*xkp(2) + REAL(shift(3), dp)*xkp(3))
      cosine = COS(arg)
      sine = SIN(arg)
      IF (real_only .AND. ABS(sine) > 1.e-12_dp) THEN
         CALL cp_abort(__LOCATION__, "Real k-point wavefunctions cannot represent symmetry phases")
      END IF
      phase(:, 1) = [cosine, -sine]
      phase(:, 2) = MERGE(-1.0_dp, 1.0_dp, time_reversal)*[sine, cosine]
   END SUBROUTINE rskp_phase_factors

! **************************************************************************************************
!> \brief Allocate a dense work matrix with the requested shape
!> \param work dense work matrix
!> \param nrow number of rows
!> \param ncol number of columns
! **************************************************************************************************
   SUBROUTINE rskp_work_matrix(work, nrow, ncol)

      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(INOUT)                                   :: work
      INTEGER, INTENT(IN)                                :: nrow, ncol

      IF (ALLOCATED(work)) THEN
         IF (SIZE(work, 1) == nrow .AND. SIZE(work, 2) == ncol) RETURN
         DEALLOCATE (work)
      END IF
      ALLOCATE (work(nrow, ncol))

   END SUBROUTINE rskp_work_matrix

END MODULE kpoint_transform
