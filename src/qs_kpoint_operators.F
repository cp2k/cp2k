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
!> \author MI
! **************************************************************************************************
MODULE qs_kpoint_operators
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_scale_and_add_fm
   USE cp_cfm_types,                    ONLY: cp_cfm_to_fm,&
                                              cp_cfm_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_deallocate_matrix,&
                                              dbcsr_desymmetrize,&
                                              dbcsr_p_type,&
                                              dbcsr_set,&
                                              dbcsr_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
   USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                              fm_pool_give_back_fm
   USE cp_fm_struct,                    ONLY: cp_fm_struct_type
   USE cp_fm_types,                     ONLY: copy_info_type,&
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
   USE qs_matrix_pools,                 ONLY: qs_matrix_pools_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_kpoint_operators'

   PUBLIC :: kpoint_operator_raw_start, kpoint_operator_raw_finish_cfm, kpoint_operator_cfm_to_mo
   PUBLIC :: kpoint_operator_context_create, kpoint_operator_context_release

   INTEGER, PARAMETER, PUBLIC :: kpoint_slot_ks = 1, &
                                 kpoint_slot_s = 2, &
                                 kpoint_slot_t = 3

   !> Number of defined slots; grows only inside this module.
   INTEGER, PARAMETER, PRIVATE :: kpoint_num_slots = 3

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
      ! Host references, read fresh at every use: sab_nl and cell_to_index
      ! are reallocated at every energy evaluation.
      TYPE(kpoint_type), POINTER                          :: kpoints => NULL()
      TYPE(cp_fm_struct_type), POINTER                    :: ao_ao_fmstruct => NULL()
      TYPE(qs_matrix_pools_type), POINTER                 :: group_mpools => NULL()
      TYPE(cp_fm_pool_p_type), DIMENSION(:), POINTER      :: group_fm_pools => NULL()
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                          :: trio_sab => NULL()
      ! Owned scratch, created on first use and released with the context.
      TYPE(dbcsr_type), POINTER                           :: rmatrix => NULL()
      TYPE(dbcsr_type), POINTER                           :: cmatrix => NULL()
      TYPE(dbcsr_type), POINTER                           :: tmpmat => NULL()
      TYPE(cp_fm_type), ALLOCATABLE, DIMENSION(:)         :: fmwork
      TYPE(cp_fm_type), ALLOCATABLE                       :: fmlocal
      TYPE(cp_fm_type)                                    :: fmdummy
      TYPE(copy_info_type), ALLOCATABLE, DIMENSION(:, :, :) :: info
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)            :: status
   END TYPE kpoint_operator_context_type

CONTAINS

! **************************************************************************************************
!> \brief Transform one real-space operator to one k point and start the transfer of
!>        the result to the k-point group.
!>        The routine densifies the real part into fm_re and the imaginary part into
!>        fm_im. Densify means: copy the blocks of a DBCSR matrix into a full matrix.
!>        It then starts the transfer of both parts to the target matrix.
!>        Each transfer receives its data in a private buffer. For this reason, several
!>        starts can use the same target matrix.
!>        Finish a start with kpoint_operator_raw_finish_cfm when its target received data.
!>        On the groups that do not own this k point, my_kpgrp routes the transfer
!>        to fmdummy; clean such a start up with cp_fm_cleanup_copy_general and
!>        do not finish it.
!> \param rsmat real-space image matrices of the operator
!> \param ispin spin component of rsmat (1 for spin-free matrices such as S)
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
!> \param fmlocal transfer target on the group that owns this k point
!> \param fmdummy transfer target on the other groups. This keeps the collective
!>        calls balanced.
!> \param my_kpgrp true when the calling group owns this k point
!> \param para_env communicator that spans the source and the target distributions
!> \param info_re transfer state of the real part
!> \param info_im transfer state of the imaginary part
! **************************************************************************************************
   SUBROUTINE kpoint_operator_raw_start(rsmat, ispin, ik, xkp, cell_to_index, sab_nl, grid, use_grid, &
                                        rmatrix, cmatrix, tmpmat, fm_re, fm_im, fmlocal, fmdummy, my_kpgrp, &
                                        para_env, info_re, info_im)

      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rsmat
      INTEGER, INTENT(IN)                                :: ispin, ik
      REAL(KIND=dp), DIMENSION(3), INTENT(IN)            :: xkp
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(rskp_grid_type), INTENT(IN), OPTIONAL         :: grid
      LOGICAL, INTENT(IN)                                :: use_grid
      TYPE(dbcsr_type)                                   :: rmatrix, cmatrix, tmpmat
      TYPE(cp_fm_type)                                   :: fm_re, fm_im, fmlocal, fmdummy
      LOGICAL, INTENT(IN)                                :: my_kpgrp
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(copy_info_type)                               :: info_re, info_im

      CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_operator_raw_start'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (use_grid) THEN
         CPASSERT(PRESENT(grid))
         CALL rskp_transform_grid_extract(grid, ik, rmatrix, cmatrix)
      ELSE
         CALL dbcsr_set(rmatrix, 0.0_dp)
         CALL dbcsr_set(cmatrix, 0.0_dp)
         CALL rskp_transform(rmatrix=rmatrix, cmatrix=cmatrix, rsmat=rsmat, ispin=ispin, &
                             xkp=xkp, cell_to_index=cell_to_index, sab_nl=sab_nl)
      END IF
      CALL dbcsr_desymmetrize(rmatrix, tmpmat)
      CALL copy_dbcsr_to_fm(tmpmat, fm_re)
      CALL dbcsr_desymmetrize(cmatrix, tmpmat)
      CALL copy_dbcsr_to_fm(tmpmat, fm_im)

      IF (my_kpgrp) THEN
         CALL cp_fm_start_copy_general(fm_re, fmlocal, para_env, info_re)
         CALL cp_fm_start_copy_general(fm_im, fmlocal, para_env, info_im)
      ELSE
         CALL cp_fm_start_copy_general(fm_re, fmdummy, para_env, info_re)
         CALL cp_fm_start_copy_general(fm_im, fmdummy, para_env, info_im)
      END IF

      CALL timestop(handle)

   END SUBROUTINE kpoint_operator_raw_start

! **************************************************************************************************
!> \brief Finish one transfer that kpoint_operator_raw_start started and add the received
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
!>        its matrix struct is the template for the context work buffers
! **************************************************************************************************
   SUBROUTINE kpoint_operator_context_create(ctx, kpoints, ao_ao_fm)

      TYPE(kpoint_operator_context_type), INTENT(OUT)    :: ctx
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(cp_fm_type), INTENT(IN)                       :: ao_ao_fm

      CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_operator_context_create'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(kpoints))
      ctx%kpoints => kpoints
      CALL cp_fm_get_info(ao_ao_fm, matrix_struct=ctx%ao_ao_fmstruct)
      CPASSERT(ASSOCIATED(ctx%ao_ao_fmstruct))

      CALL timestop(handle)

   END SUBROUTINE kpoint_operator_context_create

! **************************************************************************************************
!> \brief Release the context together with every buffer it created.
!>        Call once, at the exit of the routine that created the context.
!> \param ctx context to release
! **************************************************************************************************
   SUBROUTINE kpoint_operator_context_release(ctx)

      TYPE(kpoint_operator_context_type), INTENT(INOUT)  :: ctx

      CHARACTER(LEN=*), PARAMETER :: routineN = 'kpoint_operator_context_release'

      INTEGER                                            :: handle, iwork

      CALL timeset(routineN, handle)

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

END MODULE qs_kpoint_operators
