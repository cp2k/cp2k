!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Storage of past states of the qs_env.
!>      Methods to interpolate (or actually normally extrapolate) the
!>      new guess for density and wavefunctions.
!> \note
!>      Most of the last snapshot should actually be in qs_env, but taking
!>      advantage of it would make the programming much convoluted
!> \par History
!>      02.2003 created [fawzi]
!>      11.2003 Joost VandeVondele : Implemented Nth order PS extrapolation
!>      02.2005 modified for KG_GPW [MI]
!> \author fawzi
! **************************************************************************************************
MODULE qs_wf_history_methods
   USE bibliography,                    ONLY: Kolafa2004,&
                                              Kuhne2007,&
                                              VandeVondele2005a,&
                                              cite_reference
   USE cp_cfm_basic_linalg,             ONLY: cp_cfm_column_scale,&
                                              cp_cfm_gemm,&
                                              cp_cfm_scale_and_add,&
                                              cp_cfm_scale_and_add_fm,&
                                              cp_cfm_triangular_multiply
   USE cp_cfm_cholesky,                 ONLY: cp_cfm_cholesky_decompose
   USE cp_cfm_diag,                     ONLY: cp_cfm_heevd
   USE cp_cfm_types,                    ONLY: cp_cfm_create,&
                                              cp_cfm_release,&
                                              cp_cfm_set_all,&
                                              cp_cfm_to_cfm,&
                                              cp_cfm_to_fm,&
                                              cp_cfm_type,&
                                              cp_fm_to_cfm
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_create, dbcsr_deallocate_matrix, dbcsr_desymmetrize, &
        dbcsr_p_type, dbcsr_set, dbcsr_type, dbcsr_type_antisymmetric, dbcsr_type_no_symmetry, &
        dbcsr_type_symmetric
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale,&
                                              cp_fm_scale_and_add
   USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                              cp_fm_pool_type,&
                                              fm_pool_create_fm,&
                                              fm_pool_give_back_fm,&
                                              fm_pools_create_fm_vect,&
                                              fm_pools_give_back_fm_vect
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: &
        copy_info_type, cp_fm_cleanup_copy_general, cp_fm_create, cp_fm_finish_copy_general, &
        cp_fm_get_info, cp_fm_release, cp_fm_start_copy_general, cp_fm_to_fm, cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type,&
                                              cp_to_string
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr,&
                                              low_print_level
   USE input_constants,                 ONLY: &
        wfi_aspc_nr, wfi_frozen_method_nr, wfi_linear_p_method_nr, wfi_linear_ps_method_nr, &
        wfi_linear_wf_method_nr, wfi_ps_method_nr, wfi_use_guess_method_nr, &
        wfi_use_prev_p_method_nr, wfi_use_prev_rho_r_method_nr, wfi_use_prev_wf_method_nr
   USE kinds,                           ONLY: dp
   USE kpoint_methods,                  ONLY: kpoint_density_matrices,&
                                              kpoint_density_transform,&
                                              kpoint_set_mo_occupation,&
                                              rskp_transform
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_env_type,&
                                              kpoint_type
   USE mathconstants,                   ONLY: gaussi,&
                                              z_one,&
                                              z_zero
   USE mathlib,                         ONLY: binomial
   USE message_passing,                 ONLY: mp_para_env_type
   USE parallel_gemm_api,               ONLY: parallel_gemm
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_copy
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_density_matrices,             ONLY: calculate_density_matrix
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_ks_types,                     ONLY: qs_ks_did_change
   USE qs_matrix_pools,                 ONLY: mpools_get,&
                                              qs_matrix_pools_type
   USE qs_mo_methods,                   ONLY: make_basis_cholesky,&
                                              make_basis_lowdin,&
                                              make_basis_simple,&
                                              make_basis_sm
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_rho_methods,                  ONLY: qs_rho_update_rho
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_set,&
                                              qs_rho_type
   USE qs_scf_types,                    ONLY: ot_method_nr,&
                                              qs_scf_env_type
   USE qs_wf_history_types,             ONLY: qs_wf_history_type,&
                                              qs_wf_snapshot_type,&
                                              wfi_get_snapshot,&
                                              wfi_release
   USE scf_control_types,               ONLY: scf_control_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   LOGICAL, PRIVATE, PARAMETER :: debug_this_module = .TRUE.
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_wf_history_methods'

   PUBLIC :: wfi_create, wfi_update, wfi_create_for_kp, &
             wfi_extrapolate, wfi_get_method_label, &
             reorthogonalize_vectors, wfi_purge_history

CONTAINS

! **************************************************************************************************
!> \brief allocates and initialize a wavefunction snapshot
!> \param snapshot the snapshot to create
!> \par History
!>      02.2003 created [fawzi]
!>      02.2005 added wf_mol [MI]
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE wfs_create(snapshot)
      TYPE(qs_wf_snapshot_type), INTENT(OUT)             :: snapshot

      NULLIFY (snapshot%wf, snapshot%rho_r, &
               snapshot%rho_g, snapshot%rho_ao, snapshot%rho_ao_kp, &
               snapshot%overlap, snapshot%wf_kp, snapshot%overlap_cfm_kp, &
               snapshot%rho_frozen)
      snapshot%dt = 1.0_dp
   END SUBROUTINE wfs_create

! **************************************************************************************************
!> \brief updates the given snapshot
!> \param snapshot the snapshot to be updated
!> \param wf_history the history
!> \param qs_env the qs_env that should be snapshotted
!> \param dt the time of the snapshot (wrt. to the previous snapshot)
!> \par History
!>      02.2003 created [fawzi]
!>      02.2005 added kg_fm_mol_set for KG_GPW [MI]
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE wfs_update(snapshot, wf_history, qs_env, dt)
      TYPE(qs_wf_snapshot_type), POINTER                 :: snapshot
      TYPE(qs_wf_history_type), POINTER                  :: wf_history
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), INTENT(in), OPTIONAL                :: dt

      CHARACTER(len=*), PARAMETER                        :: routineN = 'wfs_update'

      INTEGER                                            :: handle, ic, igroup, ik, ikp, img, &
                                                            indx_ft, ispin, kplocal, nc, nimg, &
                                                            nkp_all, nkp_grps, nspin_kp, nspins
      INTEGER, DIMENSION(2)                              :: kp_range
      INTEGER, DIMENSION(:, :), POINTER                  :: kp_dist
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: my_kpgrp
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: xkp
      TYPE(copy_info_type), ALLOCATABLE, DIMENSION(:, :) :: info_ft
      TYPE(cp_fm_pool_p_type), DIMENSION(:), POINTER     :: ao_ao_fm_pools, ao_mo_pools
      TYPE(cp_fm_struct_type), POINTER                   :: ao_ao_struct_ft
      TYPE(cp_fm_type)                                   :: fmdummy_ft, fmlocal_ft
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s, rho_ao
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s_kp, rho_ao_kp
      TYPE(dbcsr_type), POINTER                          :: cmat_ft, rmat_ft, tmpmat_ft
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(mp_para_env_type), POINTER                    :: para_env_ft
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(pw_c1d_gs_type), DIMENSION(:), POINTER        :: rho_g
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_r3d_rs_type), DIMENSION(:), POINTER        :: rho_r
      TYPE(qs_matrix_pools_type), POINTER                :: mpools_kp
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_scf_env_type), POINTER                     :: scf_env

      CALL timeset(routineN, handle)

      NULLIFY (pw_env, auxbas_pw_pool, ao_mo_pools, ao_ao_fm_pools, dft_control, mos, mo_coeff, &
               rho, rho_r, rho_g, rho_ao, matrix_s, matrix_s_kp, kpoints, kp, &
               kp_dist, cell_to_index, xkp, sab_nl, scf_env, mpools_kp, para_env_ft, &
               rmat_ft, cmat_ft, tmpmat_ft, ao_ao_struct_ft)
      CALL get_qs_env(qs_env, pw_env=pw_env, &
                      dft_control=dft_control, rho=rho)
      CALL mpools_get(qs_env%mpools, ao_mo_fm_pools=ao_mo_pools)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)

      CPASSERT(ASSOCIATED(wf_history))
      CPASSERT(ASSOCIATED(dft_control))
      IF (.NOT. ASSOCIATED(snapshot)) THEN
         ALLOCATE (snapshot)
         CALL wfs_create(snapshot)
      END IF
      CPASSERT(wf_history%ref_count > 0)

      nspins = dft_control%nspins
      snapshot%dt = 1.0_dp
      IF (PRESENT(dt)) snapshot%dt = dt
      IF (wf_history%store_wf) THEN
         CALL get_qs_env(qs_env, mos=mos)
         IF (.NOT. ASSOCIATED(snapshot%wf)) THEN
            CALL fm_pools_create_fm_vect(ao_mo_pools, snapshot%wf, &
                                         name="ws_snap-ws")
            CPASSERT(nspins == SIZE(snapshot%wf))
         END IF
         DO ispin = 1, nspins
            CALL get_mo_set(mos(ispin), mo_coeff=mo_coeff)
            CALL cp_fm_to_fm(mo_coeff, snapshot%wf(ispin))
         END DO
      ELSE
         CALL fm_pools_give_back_fm_vect(ao_mo_pools, snapshot%wf)
      END IF

      IF (wf_history%store_rho_r) THEN
         CALL qs_rho_get(rho, rho_r=rho_r)
         CPASSERT(ASSOCIATED(rho_r))
         IF (.NOT. ASSOCIATED(snapshot%rho_r)) THEN
            ALLOCATE (snapshot%rho_r(nspins))
            DO ispin = 1, nspins
               CALL auxbas_pw_pool%create_pw(snapshot%rho_r(ispin))
            END DO
         END IF
         DO ispin = 1, nspins
            CALL pw_copy(rho_r(ispin), snapshot%rho_r(ispin))
         END DO
      ELSE IF (ASSOCIATED(snapshot%rho_r)) THEN
         DO ispin = 1, SIZE(snapshot%rho_r)
            CALL auxbas_pw_pool%give_back_pw(snapshot%rho_r(ispin))
         END DO
         DEALLOCATE (snapshot%rho_r)
      END IF

      IF (wf_history%store_rho_g) THEN
         CALL qs_rho_get(rho, rho_g=rho_g)
         CPASSERT(ASSOCIATED(rho_g))
         IF (.NOT. ASSOCIATED(snapshot%rho_g)) THEN
            ALLOCATE (snapshot%rho_g(nspins))
            DO ispin = 1, nspins
               CALL auxbas_pw_pool%create_pw(snapshot%rho_g(ispin))
            END DO
         END IF
         DO ispin = 1, nspins
            CALL pw_copy(rho_g(ispin), snapshot%rho_g(ispin))
         END DO
      ELSE IF (ASSOCIATED(snapshot%rho_g)) THEN
         DO ispin = 1, SIZE(snapshot%rho_g)
            CALL auxbas_pw_pool%give_back_pw(snapshot%rho_g(ispin))
         END DO
         DEALLOCATE (snapshot%rho_g)
      END IF

      IF (ASSOCIATED(snapshot%rho_ao)) THEN ! the sparsity might be different
         ! (future struct:check)
         CALL dbcsr_deallocate_matrix_set(snapshot%rho_ao)
      END IF
      IF (wf_history%store_rho_ao) THEN
         CALL qs_rho_get(rho, rho_ao=rho_ao)
         CPASSERT(ASSOCIATED(rho_ao))

         CALL dbcsr_allocate_matrix_set(snapshot%rho_ao, nspins)
         DO ispin = 1, nspins
            ALLOCATE (snapshot%rho_ao(ispin)%matrix)
            CALL dbcsr_copy(snapshot%rho_ao(ispin)%matrix, rho_ao(ispin)%matrix)
         END DO
      END IF

      IF (ASSOCIATED(snapshot%rho_ao_kp)) THEN ! the sparsity might be different
         ! (future struct:check)
         CALL dbcsr_deallocate_matrix_set(snapshot%rho_ao_kp)
      END IF
      IF (wf_history%store_rho_ao_kp) THEN
         CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp)
         CPASSERT(ASSOCIATED(rho_ao_kp))

         nimg = dft_control%nimages
         CALL dbcsr_allocate_matrix_set(snapshot%rho_ao_kp, nspins, nimg)
         DO ispin = 1, nspins
            DO img = 1, nimg
               ALLOCATE (snapshot%rho_ao_kp(ispin, img)%matrix)
               CALL dbcsr_copy(snapshot%rho_ao_kp(ispin, img)%matrix, &
                               rho_ao_kp(ispin, img)%matrix)
            END DO
         END DO
      END IF

      IF (ASSOCIATED(snapshot%overlap)) THEN ! the sparsity might be different
         ! (future struct:check)
         CALL dbcsr_deallocate_matrix(snapshot%overlap)
      END IF
      IF (wf_history%store_overlap) THEN
         CALL get_qs_env(qs_env, matrix_s=matrix_s)
         CPASSERT(ASSOCIATED(matrix_s))
         CPASSERT(ASSOCIATED(matrix_s(1)%matrix))
         ALLOCATE (snapshot%overlap)
         CALL dbcsr_copy(snapshot%overlap, matrix_s(1)%matrix)
      END IF

      CALL get_qs_env(qs_env, kpoints=kpoints)
      IF (ASSOCIATED(kpoints)) THEN
         IF (ASSOCIATED(kpoints%kp_env)) THEN
            ! --- k-point WFN snapshot: store complex MO coefficients per local k-point ---
            IF (wf_history%store_wf_kp) THEN
               CALL get_kpoint_info(kpoints, kp_range=kp_range)
               kplocal = kp_range(2) - kp_range(1) + 1
               nspin_kp = SIZE(kpoints%kp_env(1)%kpoint_env%mos, 2)
               nc = SIZE(kpoints%kp_env(1)%kpoint_env%mos, 1) ! 2=complex, 1=real

               IF (ASSOCIATED(snapshot%wf_kp)) THEN
                  DO ikp = 1, SIZE(snapshot%wf_kp, 1)
                     DO ic = 1, SIZE(snapshot%wf_kp, 2)
                        DO ispin = 1, SIZE(snapshot%wf_kp, 3)
                           CALL cp_fm_release(snapshot%wf_kp(ikp, ic, ispin))
                        END DO
                     END DO
                  END DO
                  DEALLOCATE (snapshot%wf_kp)
               END IF

               ALLOCATE (snapshot%wf_kp(kplocal, nc, nspin_kp))
               DO ikp = 1, kplocal
                  kp => kpoints%kp_env(ikp)%kpoint_env
                  DO ispin = 1, nspin_kp
                     DO ic = 1, nc
                        CALL get_mo_set(kp%mos(ic, ispin), mo_coeff=mo_coeff)
                        CALL cp_fm_create(snapshot%wf_kp(ikp, ic, ispin), &
                                          mo_coeff%matrix_struct, &
                                          name="wfkp_snap")
                        CALL cp_fm_to_fm(mo_coeff, snapshot%wf_kp(ikp, ic, ispin))
                     END DO
                  END DO
               END DO
            END IF

            ! --- k-point overlap snapshot: Fourier-transform S(R)→S(k) NOW and store as cfm ---
            ! This is critical: we MUST transform at snapshot time using the CURRENT neighbor
            ! list. Storing S(R) and re-transforming later would use a stale neighbor list,
            ! producing wrong S(k) when the neighbor list changes during MD.
            IF (wf_history%store_overlap_kp) THEN
               CALL get_qs_env(qs_env, matrix_s_kp=matrix_s_kp, scf_env=scf_env)
               CALL get_kpoint_info(kpoints, nkp=nkp_all, xkp=xkp, kp_range=kp_range, &
                                    nkp_groups=nkp_grps, kp_dist=kp_dist, &
                                    sab_nl=sab_nl, cell_to_index=cell_to_index)
               kplocal = kp_range(2) - kp_range(1) + 1
               para_env_ft => kpoints%blacs_env_all%para_env

               ! Allocate dbcsr work matrices for FT (same pattern as do_general_diag_kp)
               ALLOCATE (rmat_ft, cmat_ft, tmpmat_ft)
               CALL dbcsr_create(rmat_ft, template=matrix_s_kp(1, 1)%matrix, &
                                 matrix_type=dbcsr_type_symmetric)
               CALL dbcsr_create(cmat_ft, template=matrix_s_kp(1, 1)%matrix, &
                                 matrix_type=dbcsr_type_antisymmetric)
               CALL dbcsr_create(tmpmat_ft, template=matrix_s_kp(1, 1)%matrix, &
                                 matrix_type=dbcsr_type_no_symmetry)
               CALL cp_dbcsr_alloc_block_from_nbl(rmat_ft, sab_nl)
               CALL cp_dbcsr_alloc_block_from_nbl(cmat_ft, sab_nl)

               ! Get kp-subgroup FM from pool
               CALL get_kpoint_info(kpoints, mpools=mpools_kp)
               CALL mpools_get(mpools_kp, ao_ao_fm_pools=ao_ao_fm_pools)
               CALL fm_pool_create_fm(ao_ao_fm_pools(1)%pool, fmlocal_ft)

               ! Release old snapshot if present
               IF (ASSOCIATED(snapshot%overlap_cfm_kp)) THEN
                  DO ikp = 1, SIZE(snapshot%overlap_cfm_kp)
                     CALL cp_cfm_release(snapshot%overlap_cfm_kp(ikp))
                  END DO
                  DEALLOCATE (snapshot%overlap_cfm_kp)
               END IF
               ALLOCATE (snapshot%overlap_cfm_kp(kplocal))

               CALL cp_fm_get_info(fmlocal_ft, matrix_struct=ao_ao_struct_ft)

               ! Communication info array
               ALLOCATE (info_ft(kplocal*nkp_grps, 2))

               ! Phase A: Start async FT + redistribute for each k-point
               indx_ft = 0
               DO ikp = 1, kplocal
                  DO igroup = 1, nkp_grps
                     ik = kp_dist(1, igroup) + ikp - 1
                     my_kpgrp = (ik >= kp_range(1) .AND. ik <= kp_range(2))
                     indx_ft = indx_ft + 1

                     CALL dbcsr_set(rmat_ft, 0.0_dp)
                     CALL dbcsr_set(cmat_ft, 0.0_dp)
                     CALL rskp_transform(rmatrix=rmat_ft, cmatrix=cmat_ft, rsmat=matrix_s_kp, &
                                         ispin=1, xkp=xkp(1:3, ik), &
                                         cell_to_index=cell_to_index, sab_nl=sab_nl)
                     CALL dbcsr_desymmetrize(rmat_ft, tmpmat_ft)
                     CALL copy_dbcsr_to_fm(tmpmat_ft, scf_env%scf_work1(1))
                     CALL dbcsr_desymmetrize(cmat_ft, tmpmat_ft)
                     CALL copy_dbcsr_to_fm(tmpmat_ft, scf_env%scf_work1(2))

                     IF (my_kpgrp) THEN
                        CALL cp_fm_start_copy_general(scf_env%scf_work1(1), fmlocal_ft, &
                                                      para_env_ft, info_ft(indx_ft, 1))
                        CALL cp_fm_start_copy_general(scf_env%scf_work1(2), fmlocal_ft, &
                                                      para_env_ft, info_ft(indx_ft, 2))
                     ELSE
                        CALL cp_fm_start_copy_general(scf_env%scf_work1(1), fmdummy_ft, &
                                                      para_env_ft, info_ft(indx_ft, 1))
                        CALL cp_fm_start_copy_general(scf_env%scf_work1(2), fmdummy_ft, &
                                                      para_env_ft, info_ft(indx_ft, 2))
                     END IF
                  END DO
               END DO

               ! Phase B: Finish communication and assemble S(k) as cfm
               indx_ft = 0
               DO ikp = 1, kplocal
                  CALL cp_cfm_create(snapshot%overlap_cfm_kp(ikp), ao_ao_struct_ft)
                  CALL cp_cfm_set_all(snapshot%overlap_cfm_kp(ikp), z_zero)
                  DO igroup = 1, nkp_grps
                     ik = kp_dist(1, igroup) + ikp - 1
                     my_kpgrp = (ik >= kp_range(1) .AND. ik <= kp_range(2))
                     indx_ft = indx_ft + 1
                     IF (my_kpgrp) THEN
                        CALL cp_fm_finish_copy_general(fmlocal_ft, info_ft(indx_ft, 1))
                        CALL cp_cfm_scale_and_add_fm(z_zero, snapshot%overlap_cfm_kp(ikp), &
                                                     z_one, fmlocal_ft)
                        CALL cp_fm_finish_copy_general(fmlocal_ft, info_ft(indx_ft, 2))
                        CALL cp_cfm_scale_and_add_fm(z_one, snapshot%overlap_cfm_kp(ikp), &
                                                     gaussi, fmlocal_ft)
                     END IF
                  END DO
               END DO

               ! Cleanup
               DO indx_ft = 1, kplocal*nkp_grps
                  CALL cp_fm_cleanup_copy_general(info_ft(indx_ft, 1))
                  CALL cp_fm_cleanup_copy_general(info_ft(indx_ft, 2))
               END DO
               DEALLOCATE (info_ft)
               CALL fm_pool_give_back_fm(ao_ao_fm_pools(1)%pool, fmlocal_ft)
               CALL dbcsr_deallocate_matrix(rmat_ft)
               CALL dbcsr_deallocate_matrix(cmat_ft)
               CALL dbcsr_deallocate_matrix(tmpmat_ft)
            END IF
         END IF
      END IF

      IF (wf_history%store_frozen_density) THEN
         ! do nothing
         ! CALL deallocate_matrix_set(snapshot%rho_frozen%rho_ao)
      END IF

      CALL timestop(handle)

   END SUBROUTINE wfs_update

! **************************************************************************************************
!> \brief ...
!> \param wf_history ...
!> \param interpolation_method_nr the tag of the method used for
!>        the extrapolation of the initial density for the next md step
!>        (see qs_wf_history_types:wfi_*_method_nr)
!> \param extrapolation_order ...
!> \param has_unit_metric ...
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE wfi_create(wf_history, interpolation_method_nr, extrapolation_order, &
                         has_unit_metric)
      TYPE(qs_wf_history_type), POINTER                  :: wf_history
      INTEGER, INTENT(in)                                :: interpolation_method_nr, &
                                                            extrapolation_order
      LOGICAL, INTENT(IN)                                :: has_unit_metric

      CHARACTER(len=*), PARAMETER                        :: routineN = 'wfi_create'

      INTEGER                                            :: handle, i

      CALL timeset(routineN, handle)

      ALLOCATE (wf_history)
      wf_history%ref_count = 1
      wf_history%memory_depth = 0
      wf_history%snapshot_count = 0
      wf_history%last_state_index = 1
      wf_history%store_wf = .FALSE.
      wf_history%store_rho_r = .FALSE.
      wf_history%store_rho_g = .FALSE.
      wf_history%store_rho_ao = .FALSE.
      wf_history%store_rho_ao_kp = .FALSE.
      wf_history%store_overlap = .FALSE.
      wf_history%store_wf_kp = .FALSE.
      wf_history%store_overlap_kp = .FALSE.
      wf_history%store_frozen_density = .FALSE.
      NULLIFY (wf_history%past_states)

      wf_history%interpolation_method_nr = interpolation_method_nr

      SELECT CASE (wf_history%interpolation_method_nr)
      CASE (wfi_use_guess_method_nr)
         wf_history%memory_depth = 0
      CASE (wfi_use_prev_wf_method_nr)
         wf_history%memory_depth = 0
      CASE (wfi_use_prev_p_method_nr)
         wf_history%memory_depth = 1
         wf_history%store_rho_ao = .TRUE.
      CASE (wfi_use_prev_rho_r_method_nr)
         wf_history%memory_depth = 1
         wf_history%store_rho_ao = .TRUE.
      CASE (wfi_linear_wf_method_nr)
         wf_history%memory_depth = 2
         wf_history%store_wf = .TRUE.
      CASE (wfi_linear_p_method_nr)
         wf_history%memory_depth = 2
         wf_history%store_rho_ao = .TRUE.
      CASE (wfi_linear_ps_method_nr)
         wf_history%memory_depth = 2
         wf_history%store_wf = .TRUE.
         IF (.NOT. has_unit_metric) wf_history%store_overlap = .TRUE.
      CASE (wfi_ps_method_nr)
         CALL cite_reference(VandeVondele2005a)
         wf_history%memory_depth = extrapolation_order + 1
         wf_history%store_wf = .TRUE.
         IF (.NOT. has_unit_metric) wf_history%store_overlap = .TRUE.
      CASE (wfi_frozen_method_nr)
         wf_history%memory_depth = 1
         wf_history%store_frozen_density = .TRUE.
      CASE (wfi_aspc_nr)
         wf_history%memory_depth = extrapolation_order + 2
         wf_history%store_wf = .TRUE.
         wf_history%store_wf_kp = .TRUE.
         IF (.NOT. has_unit_metric) THEN
            wf_history%store_overlap = .TRUE.
            wf_history%store_overlap_kp = .TRUE.
         END IF
      CASE default
         CALL cp_abort(__LOCATION__, &
                       "Unknown interpolation method: "// &
                       TRIM(ADJUSTL(cp_to_string(interpolation_method_nr))))
         wf_history%interpolation_method_nr = wfi_use_prev_rho_r_method_nr
      END SELECT
      ALLOCATE (wf_history%past_states(wf_history%memory_depth))

      DO i = 1, SIZE(wf_history%past_states)
         NULLIFY (wf_history%past_states(i)%snapshot)
      END DO

      CALL timestop(handle)
   END SUBROUTINE wfi_create

! **************************************************************************************************
!> \brief Adapts wf_history storage flags for k-point calculations.
!>        For ASPC, switches from Gamma WFN storage to k-point WFN storage.
!>        Other WFN-based methods remain blocked.
!> \param wf_history ...
!> \par History
!>      06.2015 created [jhu]
!> \author jhu
! **************************************************************************************************
   SUBROUTINE wfi_create_for_kp(wf_history)
      TYPE(qs_wf_history_type), POINTER                  :: wf_history

      CPASSERT(ASSOCIATED(wf_history))
      IF (wf_history%store_rho_ao) THEN
         wf_history%store_rho_ao_kp = .TRUE.
         wf_history%store_rho_ao = .FALSE.
      END IF
      ! KP-compatible ASPC: switch from Gamma WFN storage to KP WFN storage
      IF (wf_history%store_wf_kp) THEN
         wf_history%store_wf = .FALSE.
         wf_history%store_overlap = .FALSE.
         ! store_wf_kp and store_overlap_kp remain TRUE
      ELSE
         ! Non-ASPC methods are still temporarily blocked
         IF (wf_history%store_wf) THEN
            CPABORT("WFN based interpolation method not implemented for kpoints.")
         END IF
         IF (wf_history%store_overlap) THEN
            CPABORT("Inconsistent interpolation method for kpoints.")
         END IF
      END IF
      IF (wf_history%store_frozen_density) THEN
         CPABORT("Frozen density initialization method not possible for kpoints.")
      END IF

   END SUBROUTINE wfi_create_for_kp

! **************************************************************************************************
!> \brief returns a string describing the interpolation method
!> \param method_nr ...
!> \return ...
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
   FUNCTION wfi_get_method_label(method_nr) RESULT(res)
      INTEGER, INTENT(in)                                :: method_nr
      CHARACTER(len=30)                                  :: res

      res = "unknown"
      SELECT CASE (method_nr)
      CASE (wfi_use_prev_p_method_nr)
         res = "previous_p"
      CASE (wfi_use_prev_wf_method_nr)
         res = "previous_wf"
      CASE (wfi_use_prev_rho_r_method_nr)
         res = "previous_rho_r"
      CASE (wfi_use_guess_method_nr)
         res = "initial_guess"
      CASE (wfi_linear_wf_method_nr)
         res = "mo linear"
      CASE (wfi_linear_p_method_nr)
         res = "P linear"
      CASE (wfi_linear_ps_method_nr)
         res = "PS linear"
      CASE (wfi_ps_method_nr)
         res = "PS Nth order"
      CASE (wfi_frozen_method_nr)
         res = "frozen density approximation"
      CASE (wfi_aspc_nr)
         res = "ASPC"
      CASE default
         CALL cp_abort(__LOCATION__, &
                       "Unknown interpolation method: "// &
                       TRIM(ADJUSTL(cp_to_string(method_nr))))
      END SELECT
   END FUNCTION wfi_get_method_label

! **************************************************************************************************
!> \brief calculates the new starting state for the scf for the next
!>      wf optimization
!> \param wf_history the previous history needed to extrapolate
!> \param qs_env the qs env with the latest result, and that will contain
!>        the new starting state
!> \param dt the time at which to extrapolate (wrt. to the last snapshot)
!> \param extrapolation_method_nr returns the extrapolation method used
!> \param orthogonal_wf ...
!> \par History
!>      02.2003 created [fawzi]
!>      11.2003 Joost VandeVondele : Implemented Nth order PS extrapolation
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE wfi_extrapolate(wf_history, qs_env, dt, extrapolation_method_nr, &
                              orthogonal_wf)
      TYPE(qs_wf_history_type), POINTER                  :: wf_history
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), INTENT(IN)                          :: dt
      INTEGER, INTENT(OUT), OPTIONAL                     :: extrapolation_method_nr
      LOGICAL, INTENT(OUT), OPTIONAL                     :: orthogonal_wf

      CHARACTER(len=*), PARAMETER                        :: routineN = 'wfi_extrapolate'

      INTEGER                                            :: actual_extrapolation_method_nr, handle, &
                                                            i, img, io_unit, ispin, k, n, nmo, &
                                                            nvec, print_level
      LOGICAL                                            :: do_kpoints, my_orthogonal_wf, use_overlap
      REAL(KIND=dp)                                      :: alpha, t0, t1, t2
      TYPE(cp_fm_pool_p_type), DIMENSION(:), POINTER     :: ao_mo_fm_pools
      TYPE(cp_fm_struct_type), POINTER                   :: matrix_struct, matrix_struct_new
      TYPE(cp_fm_type)                                   :: csc, fm_tmp
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: rho_ao, rho_frozen_ao
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao_kp
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mos
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_wf_snapshot_type), POINTER                 :: t0_state, t1_state

      NULLIFY (mos, ao_mo_fm_pools, t0_state, t1_state, mo_coeff, &
               rho, rho_ao, rho_frozen_ao)

      use_overlap = wf_history%store_overlap

      CALL timeset(routineN, handle)
      logger => cp_get_default_logger()
      print_level = logger%iter_info%print_level
      io_unit = cp_print_key_unit_nr(logger, qs_env%input, "DFT%SCF%PRINT%PROGRAM_RUN_INFO", &
                                     extension=".scfLog")

      CPASSERT(ASSOCIATED(wf_history))
      CPASSERT(wf_history%ref_count > 0)
      CPASSERT(ASSOCIATED(qs_env))
      CALL get_qs_env(qs_env, mos=mos, rho=rho, do_kpoints=do_kpoints)
      CALL mpools_get(qs_env%mpools, ao_mo_fm_pools=ao_mo_fm_pools)
      ! chooses the method for this extrapolation
      IF (wf_history%snapshot_count < 1) THEN
         actual_extrapolation_method_nr = wfi_use_guess_method_nr
      ELSE
         actual_extrapolation_method_nr = wf_history%interpolation_method_nr
      END IF

      SELECT CASE (actual_extrapolation_method_nr)
      CASE (wfi_linear_wf_method_nr)
         IF (wf_history%snapshot_count < 2) THEN
            actual_extrapolation_method_nr = wfi_use_prev_wf_method_nr
         END IF
      CASE (wfi_linear_p_method_nr)
         IF (wf_history%snapshot_count < 2) THEN
            actual_extrapolation_method_nr = wfi_use_prev_wf_method_nr
         END IF
      CASE (wfi_linear_ps_method_nr)
         IF (wf_history%snapshot_count < 2) THEN
            actual_extrapolation_method_nr = wfi_use_prev_wf_method_nr
         END IF
      END SELECT

      IF (PRESENT(extrapolation_method_nr)) &
         extrapolation_method_nr = actual_extrapolation_method_nr
      my_orthogonal_wf = .FALSE.

      SELECT CASE (actual_extrapolation_method_nr)
      CASE (wfi_frozen_method_nr)
         CPASSERT(.NOT. do_kpoints)
         t0_state => wfi_get_snapshot(wf_history, wf_index=1)
         CPASSERT(ASSOCIATED(t0_state%rho_frozen))

         nvec = MIN(wf_history%memory_depth, wf_history%snapshot_count)
         CALL wfi_set_history_variables(qs_env=qs_env, nvec=nvec)

         CALL qs_rho_get(t0_state%rho_frozen, rho_ao=rho_frozen_ao)
         CALL qs_rho_get(rho, rho_ao=rho_ao)
         DO ispin = 1, SIZE(rho_frozen_ao)
            CALL dbcsr_copy(rho_ao(ispin)%matrix, &
                            rho_frozen_ao(ispin)%matrix, &
                            keep_sparsity=.TRUE.)
         END DO
         !FM updating rho_ao directly with t0_state%rho_ao would have the
         !FM wrong matrix structure
         CALL qs_rho_update_rho(rho, qs_env=qs_env)
         CALL qs_ks_did_change(qs_env%ks_env, rho_changed=.TRUE.)

         my_orthogonal_wf = .FALSE.
      CASE (wfi_use_prev_rho_r_method_nr)
         t0_state => wfi_get_snapshot(wf_history, wf_index=1)
         nvec = MIN(wf_history%memory_depth, wf_history%snapshot_count)
         CALL wfi_set_history_variables(qs_env=qs_env, nvec=nvec)
         IF (do_kpoints) THEN
            CPASSERT(ASSOCIATED(t0_state%rho_ao_kp))
            CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp)
            DO ispin = 1, SIZE(t0_state%rho_ao_kp, 1)
               DO img = 1, SIZE(t0_state%rho_ao_kp, 2)
                  IF (img > SIZE(rho_ao_kp, 2)) THEN
                     CPWARN("Change in cell neighborlist: might affect quality of initial guess")
                  ELSE
                     CALL dbcsr_copy(rho_ao_kp(ispin, img)%matrix, &
                                     t0_state%rho_ao_kp(ispin, img)%matrix, &
                                     keep_sparsity=.TRUE.)
                  END IF
               END DO
            END DO
         ELSE
            CPASSERT(ASSOCIATED(t0_state%rho_ao))
            CALL qs_rho_get(rho, rho_ao=rho_ao)
            DO ispin = 1, SIZE(t0_state%rho_ao)
               CALL dbcsr_copy(rho_ao(ispin)%matrix, &
                               t0_state%rho_ao(ispin)%matrix, &
                               keep_sparsity=.TRUE.)
            END DO
         END IF
         ! Why is rho_g valid at this point ?
         CALL qs_rho_set(rho, rho_g_valid=.TRUE.)

         ! does nothing
      CASE (wfi_use_prev_wf_method_nr)
         my_orthogonal_wf = .TRUE.
         nvec = MIN(wf_history%memory_depth, wf_history%snapshot_count)
         CALL wfi_set_history_variables(qs_env=qs_env, nvec=nvec)

         IF (do_kpoints) THEN
            CALL wfi_use_prev_wf_kp(qs_env, io_unit, print_level)
         ELSE
            CALL qs_rho_get(rho, rho_ao=rho_ao)
            DO ispin = 1, SIZE(mos)
               CALL get_mo_set(mos(ispin), mo_coeff=mo_coeff, nmo=nmo)
               CALL reorthogonalize_vectors(qs_env, v_matrix=mo_coeff, n_col=nmo)
               CALL calculate_density_matrix(mo_set=mos(ispin), density_matrix=rho_ao(ispin)%matrix)
            END DO
            CALL qs_rho_update_rho(rho, qs_env=qs_env)
            CALL qs_ks_did_change(qs_env%ks_env, rho_changed=.TRUE.)
         END IF

      CASE (wfi_use_prev_p_method_nr)
         t0_state => wfi_get_snapshot(wf_history, wf_index=1)
         nvec = MIN(wf_history%memory_depth, wf_history%snapshot_count)
         CALL wfi_set_history_variables(qs_env=qs_env, nvec=nvec)
         IF (do_kpoints) THEN
            CPASSERT(ASSOCIATED(t0_state%rho_ao_kp))
            CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp)
            DO ispin = 1, SIZE(t0_state%rho_ao_kp, 1)
               DO img = 1, SIZE(t0_state%rho_ao_kp, 2)
                  IF (img > SIZE(rho_ao_kp, 2)) THEN
                     CPWARN("Change in cell neighborlist: might affect quality of initial guess")
                  ELSE
                     CALL dbcsr_copy(rho_ao_kp(ispin, img)%matrix, &
                                     t0_state%rho_ao_kp(ispin, img)%matrix, &
                                     keep_sparsity=.TRUE.)
                  END IF
               END DO
            END DO
         ELSE
            CPASSERT(ASSOCIATED(t0_state%rho_ao))
            CALL qs_rho_get(rho, rho_ao=rho_ao)
            DO ispin = 1, SIZE(t0_state%rho_ao)
               CALL dbcsr_copy(rho_ao(ispin)%matrix, &
                               t0_state%rho_ao(ispin)%matrix, &
                               keep_sparsity=.TRUE.)
            END DO
         END IF
         !FM updating rho_ao directly with t0_state%rho_ao would have the
         !FM wrong matrix structure
         CALL qs_rho_update_rho(rho, qs_env=qs_env)
         CALL qs_ks_did_change(qs_env%ks_env, rho_changed=.TRUE.)

      CASE (wfi_use_guess_method_nr)
         !FM more clean to do it here, but it
         !FM might need to read a file (restart) and thus globenv
         !FM I do not want globenv here, thus done by the caller
         !FM (btw. it also needs the eigensolver, and unless you relocate it
         !FM gives circular dependencies)
         nvec = MIN(wf_history%memory_depth, wf_history%snapshot_count)
         CALL wfi_set_history_variables(qs_env=qs_env, nvec=nvec)
      CASE (wfi_linear_wf_method_nr)
         CPASSERT(.NOT. do_kpoints)
         t0_state => wfi_get_snapshot(wf_history, wf_index=2)
         t1_state => wfi_get_snapshot(wf_history, wf_index=1)
         CPASSERT(ASSOCIATED(t0_state))
         CPASSERT(ASSOCIATED(t1_state))
         CPASSERT(ASSOCIATED(t0_state%wf))
         CPASSERT(ASSOCIATED(t1_state%wf))
         nvec = MIN(wf_history%memory_depth, wf_history%snapshot_count)
         CALL wfi_set_history_variables(qs_env=qs_env, nvec=nvec)

         my_orthogonal_wf = .TRUE.
         t0 = 0.0_dp
         t1 = t1_state%dt
         t2 = t1 + dt
         CALL qs_rho_get(rho, rho_ao=rho_ao)
         DO ispin = 1, SIZE(mos)
            CALL get_mo_set(mos(ispin), mo_coeff=mo_coeff, &
                            nmo=nmo)
            CALL cp_fm_scale_and_add(alpha=0.0_dp, &
                                     matrix_a=mo_coeff, &
                                     matrix_b=t1_state%wf(ispin), &
                                     beta=(t2 - t0)/(t1 - t0))
            ! this copy should be unnecessary
            CALL cp_fm_scale_and_add(alpha=1.0_dp, &
                                     matrix_a=mo_coeff, &
                                     beta=(t1 - t2)/(t1 - t0), matrix_b=t0_state%wf(ispin))
            CALL reorthogonalize_vectors(qs_env, &
                                         v_matrix=mo_coeff, &
                                         n_col=nmo)
            CALL calculate_density_matrix(mo_set=mos(ispin), &
                                          density_matrix=rho_ao(ispin)%matrix)
         END DO
         CALL qs_rho_update_rho(rho, qs_env=qs_env)

         CALL qs_ks_did_change(qs_env%ks_env, &
                               rho_changed=.TRUE.)
      CASE (wfi_linear_p_method_nr)
         t0_state => wfi_get_snapshot(wf_history, wf_index=2)
         t1_state => wfi_get_snapshot(wf_history, wf_index=1)
         CPASSERT(ASSOCIATED(t0_state))
         CPASSERT(ASSOCIATED(t1_state))
         IF (do_kpoints) THEN
            CPASSERT(ASSOCIATED(t0_state%rho_ao_kp))
            CPASSERT(ASSOCIATED(t1_state%rho_ao_kp))
         ELSE
            CPASSERT(ASSOCIATED(t0_state%rho_ao))
            CPASSERT(ASSOCIATED(t1_state%rho_ao))
         END IF
         nvec = MIN(wf_history%memory_depth, wf_history%snapshot_count)
         CALL wfi_set_history_variables(qs_env=qs_env, nvec=nvec)

         t0 = 0.0_dp
         t1 = t1_state%dt
         t2 = t1 + dt
         IF (do_kpoints) THEN
            CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp)
            DO ispin = 1, SIZE(rho_ao_kp, 1)
               DO img = 1, SIZE(rho_ao_kp, 2)
                  IF (img > SIZE(t0_state%rho_ao_kp, 2) .OR. &
                      img > SIZE(t1_state%rho_ao_kp, 2)) THEN
                     CPWARN("Change in cell neighborlist: might affect quality of initial guess")
                  ELSE
                     CALL dbcsr_add(rho_ao_kp(ispin, img)%matrix, t1_state%rho_ao_kp(ispin, img)%matrix, &
                                    alpha_scalar=0.0_dp, beta_scalar=(t2 - t0)/(t1 - t0)) ! this copy should be unnecessary
                     CALL dbcsr_add(rho_ao_kp(ispin, img)%matrix, t0_state%rho_ao_kp(ispin, img)%matrix, &
                                    alpha_scalar=1.0_dp, beta_scalar=(t1 - t2)/(t1 - t0))
                  END IF
               END DO
            END DO
         ELSE
            CALL qs_rho_get(rho, rho_ao=rho_ao)
            DO ispin = 1, SIZE(rho_ao)
               CALL dbcsr_add(rho_ao(ispin)%matrix, t1_state%rho_ao(ispin)%matrix, &
                              alpha_scalar=0.0_dp, beta_scalar=(t2 - t0)/(t1 - t0)) ! this copy should be unnecessary
               CALL dbcsr_add(rho_ao(ispin)%matrix, t0_state%rho_ao(ispin)%matrix, &
                              alpha_scalar=1.0_dp, beta_scalar=(t1 - t2)/(t1 - t0))
            END DO
         END IF
         CALL qs_rho_update_rho(rho, qs_env=qs_env)
         CALL qs_ks_did_change(qs_env%ks_env, rho_changed=.TRUE.)
      CASE (wfi_linear_ps_method_nr)
         ! wf not calculated, extract with PSC renormalized?
         ! use wf_linear?
         CPASSERT(.NOT. do_kpoints)
         t0_state => wfi_get_snapshot(wf_history, wf_index=2)
         t1_state => wfi_get_snapshot(wf_history, wf_index=1)
         CPASSERT(ASSOCIATED(t0_state))
         CPASSERT(ASSOCIATED(t1_state))
         CPASSERT(ASSOCIATED(t0_state%wf))
         CPASSERT(ASSOCIATED(t1_state%wf))
         IF (wf_history%store_overlap) THEN
            CPASSERT(ASSOCIATED(t0_state%overlap))
            CPASSERT(ASSOCIATED(t1_state%overlap))
         END IF
         nvec = MIN(wf_history%memory_depth, wf_history%snapshot_count)
         IF (nvec >= wf_history%memory_depth) THEN
            IF ((qs_env%scf_control%max_scf_hist /= 0) .AND. (qs_env%scf_control%eps_scf_hist /= 0)) THEN
               qs_env%scf_control%max_scf = qs_env%scf_control%max_scf_hist
               qs_env%scf_control%eps_scf = qs_env%scf_control%eps_scf_hist
               qs_env%scf_control%outer_scf%have_scf = .FALSE.
            ELSE IF (qs_env%scf_control%max_scf_hist /= 0) THEN
               qs_env%scf_control%max_scf = qs_env%scf_control%max_scf_hist
               qs_env%scf_control%outer_scf%have_scf = .FALSE.
            ELSE IF (qs_env%scf_control%eps_scf_hist /= 0) THEN
               qs_env%scf_control%eps_scf = qs_env%scf_control%eps_scf_hist
            END IF
         END IF

         my_orthogonal_wf = .TRUE.
         ! use PS_2=2 PS_1-PS_0
         ! C_2 comes from using PS_2 as a projector acting on C_1
         CALL qs_rho_get(rho, rho_ao=rho_ao)
         DO ispin = 1, SIZE(mos)
            NULLIFY (mo_coeff, matrix_struct, matrix_struct_new)
            CALL get_mo_set(mo_set=mos(ispin), mo_coeff=mo_coeff)
            CALL cp_fm_get_info(mo_coeff, nrow_global=n, ncol_global=k, &
                                matrix_struct=matrix_struct)
            CALL cp_fm_struct_create(matrix_struct_new, template_fmstruct=matrix_struct, &
                                     nrow_global=k, ncol_global=k)
            CALL cp_fm_create(csc, matrix_struct_new)
            CALL cp_fm_struct_release(matrix_struct_new)

            IF (use_overlap) THEN
               CALL cp_dbcsr_sm_fm_multiply(t0_state%overlap, t1_state%wf(ispin), mo_coeff, k)
               CALL parallel_gemm('T', 'N', k, k, n, 1.0_dp, t0_state%wf(ispin), mo_coeff, 0.0_dp, csc)
            ELSE
               CALL parallel_gemm('T', 'N', k, k, n, 1.0_dp, t0_state%wf(ispin), &
                                  t1_state%wf(ispin), 0.0_dp, csc)
            END IF
            CALL parallel_gemm('N', 'N', n, k, k, 1.0_dp, t0_state%wf(ispin), csc, 0.0_dp, mo_coeff)
            CALL cp_fm_release(csc)
            CALL cp_fm_scale_and_add(-1.0_dp, mo_coeff, 2.0_dp, t1_state%wf(ispin))
            CALL reorthogonalize_vectors(qs_env, &
                                         v_matrix=mo_coeff, &
                                         n_col=k)
            CALL calculate_density_matrix(mo_set=mos(ispin), &
                                          density_matrix=rho_ao(ispin)%matrix)
         END DO
         CALL qs_rho_update_rho(rho, qs_env=qs_env)
         CALL qs_ks_did_change(qs_env%ks_env, rho_changed=.TRUE.)

      CASE (wfi_ps_method_nr)
         CPASSERT(.NOT. do_kpoints)
         ! figure out the actual number of vectors to use in the extrapolation:
         nvec = MIN(wf_history%memory_depth, wf_history%snapshot_count)
         CPASSERT(nvec > 0)
         IF (nvec >= wf_history%memory_depth) THEN
            IF ((qs_env%scf_control%max_scf_hist /= 0) .AND. (qs_env%scf_control%eps_scf_hist /= 0)) THEN
               qs_env%scf_control%max_scf = qs_env%scf_control%max_scf_hist
               qs_env%scf_control%eps_scf = qs_env%scf_control%eps_scf_hist
               qs_env%scf_control%outer_scf%have_scf = .FALSE.
            ELSE IF (qs_env%scf_control%max_scf_hist /= 0) THEN
               qs_env%scf_control%max_scf = qs_env%scf_control%max_scf_hist
               qs_env%scf_control%outer_scf%have_scf = .FALSE.
            ELSE IF (qs_env%scf_control%eps_scf_hist /= 0) THEN
               qs_env%scf_control%eps_scf = qs_env%scf_control%eps_scf_hist
            END IF
         END IF

         my_orthogonal_wf = .TRUE.
         DO ispin = 1, SIZE(mos)
            NULLIFY (mo_coeff, matrix_struct, matrix_struct_new)
            CALL get_mo_set(mo_set=mos(ispin), mo_coeff=mo_coeff)
            CALL cp_fm_get_info(mo_coeff, nrow_global=n, ncol_global=k, &
                                matrix_struct=matrix_struct)
            CALL cp_fm_create(fm_tmp, matrix_struct)
            CALL cp_fm_struct_create(matrix_struct_new, template_fmstruct=matrix_struct, &
                                     nrow_global=k, ncol_global=k)
            CALL cp_fm_create(csc, matrix_struct_new)
            CALL cp_fm_struct_release(matrix_struct_new)
            ! first the most recent
            t1_state => wfi_get_snapshot(wf_history, wf_index=1)
            CALL cp_fm_to_fm(t1_state%wf(ispin), mo_coeff)
            alpha = nvec
            CALL cp_fm_scale(alpha, mo_coeff)
            CALL qs_rho_get(rho, rho_ao=rho_ao)
            DO i = 2, nvec
               t0_state => wfi_get_snapshot(wf_history, wf_index=i)
               IF (use_overlap) THEN
                  CALL cp_dbcsr_sm_fm_multiply(t0_state%overlap, t1_state%wf(ispin), fm_tmp, k)
                  CALL parallel_gemm('T', 'N', k, k, n, 1.0_dp, t0_state%wf(ispin), fm_tmp, 0.0_dp, csc)
               ELSE
                  CALL parallel_gemm('T', 'N', k, k, n, 1.0_dp, t0_state%wf(ispin), &
                                     t1_state%wf(ispin), 0.0_dp, csc)
               END IF
               CALL parallel_gemm('N', 'N', n, k, k, 1.0_dp, t0_state%wf(ispin), csc, 0.0_dp, fm_tmp)
               alpha = -1.0_dp*alpha*REAL(nvec - i + 1, dp)/REAL(i, dp)
               CALL cp_fm_scale_and_add(1.0_dp, mo_coeff, alpha, fm_tmp)
            END DO

            CALL cp_fm_release(csc)
            CALL cp_fm_release(fm_tmp)
            CALL reorthogonalize_vectors(qs_env, &
                                         v_matrix=mo_coeff, &
                                         n_col=k)
            CALL calculate_density_matrix(mo_set=mos(ispin), &
                                          density_matrix=rho_ao(ispin)%matrix)
         END DO
         CALL qs_rho_update_rho(rho, qs_env=qs_env)
         CALL qs_ks_did_change(qs_env%ks_env, rho_changed=.TRUE.)

      CASE (wfi_aspc_nr)
         CALL cite_reference(Kolafa2004)
         CALL cite_reference(Kuhne2007)
         ! figure out the actual number of vectors to use in the extrapolation:
         nvec = MIN(wf_history%memory_depth, wf_history%snapshot_count)
         CPASSERT(nvec > 0)
         IF (nvec >= wf_history%memory_depth) THEN
            IF ((qs_env%scf_control%max_scf_hist /= 0) .AND. &
                (qs_env%scf_control%eps_scf_hist /= 0)) THEN
               qs_env%scf_control%max_scf = qs_env%scf_control%max_scf_hist
               qs_env%scf_control%eps_scf = qs_env%scf_control%eps_scf_hist
               qs_env%scf_control%outer_scf%have_scf = .FALSE.
            ELSE IF (qs_env%scf_control%max_scf_hist /= 0) THEN
               qs_env%scf_control%max_scf = qs_env%scf_control%max_scf_hist
               qs_env%scf_control%outer_scf%have_scf = .FALSE.
            ELSE IF (qs_env%scf_control%eps_scf_hist /= 0) THEN
               qs_env%scf_control%eps_scf = qs_env%scf_control%eps_scf_hist
            END IF
         END IF

         IF (do_kpoints) THEN
            ! ---------- k-point ASPC: complex WFN extrapolation ----------
            CALL wfi_extrapolate_aspc_kp(wf_history, qs_env, nvec, io_unit, print_level)
            my_orthogonal_wf = .TRUE.
         ELSE
            ! ---------- Gamma-point ASPC: original real WFN extrapolation ----------
            my_orthogonal_wf = .TRUE.
            CALL qs_rho_get(rho, rho_ao=rho_ao)
            DO ispin = 1, SIZE(mos)
               NULLIFY (mo_coeff, matrix_struct, matrix_struct_new)
               CALL get_mo_set(mo_set=mos(ispin), mo_coeff=mo_coeff)
               CALL cp_fm_get_info(mo_coeff, &
                                   nrow_global=n, &
                                   ncol_global=k, &
                                   matrix_struct=matrix_struct)
               CALL cp_fm_create(fm_tmp, matrix_struct, set_zero=.TRUE.)
               CALL cp_fm_struct_create(matrix_struct_new, &
                                        template_fmstruct=matrix_struct, &
                                        nrow_global=k, &
                                        ncol_global=k)
               CALL cp_fm_create(csc, matrix_struct_new, set_zero=.TRUE.)
               CALL cp_fm_struct_release(matrix_struct_new)
               ! first the most recent
               t1_state => wfi_get_snapshot(wf_history, &
                                            wf_index=1)
               CALL cp_fm_to_fm(t1_state%wf(ispin), mo_coeff)
               alpha = REAL(4*nvec - 2, KIND=dp)/REAL(nvec + 1, KIND=dp)
               IF ((io_unit > 0) .AND. (print_level > low_print_level)) THEN
                  WRITE (UNIT=io_unit, FMT="(/,T2,A,/,/,T3,A,I0,/,/,T3,A2,I0,A4,F10.6)") &
                     "Parameters for the always stable predictor-corrector (ASPC) method:", &
                     "ASPC order: ", MAX(nvec - 2, 0), &
                     "B(", 1, ") = ", alpha
               END IF
               CALL cp_fm_scale(alpha, mo_coeff)

               DO i = 2, nvec
                  t0_state => wfi_get_snapshot(wf_history, wf_index=i)
                  IF (use_overlap) THEN
                     CALL cp_dbcsr_sm_fm_multiply(t0_state%overlap, t1_state%wf(ispin), fm_tmp, k)
                     CALL parallel_gemm('T', 'N', k, k, n, 1.0_dp, t0_state%wf(ispin), fm_tmp, 0.0_dp, csc)
                  ELSE
                     CALL parallel_gemm('T', 'N', k, k, n, 1.0_dp, t0_state%wf(ispin), &
                                        t1_state%wf(ispin), 0.0_dp, csc)
                  END IF
                  CALL parallel_gemm('N', 'N', n, k, k, 1.0_dp, t0_state%wf(ispin), csc, 0.0_dp, fm_tmp)
                  alpha = (-1.0_dp)**(i + 1)*REAL(i, KIND=dp)* &
                          binomial(2*nvec, nvec - i)/binomial(2*nvec - 2, nvec - 1)
                  IF ((io_unit > 0) .AND. (print_level > low_print_level)) THEN
                     WRITE (UNIT=io_unit, FMT="(T3,A2,I0,A4,F10.6)") &
                        "B(", i, ") = ", alpha
                  END IF
                  CALL cp_fm_scale_and_add(1.0_dp, mo_coeff, alpha, fm_tmp)
               END DO
               CALL cp_fm_release(csc)
               CALL cp_fm_release(fm_tmp)
               CALL reorthogonalize_vectors(qs_env, &
                                            v_matrix=mo_coeff, &
                                            n_col=k)
               CALL calculate_density_matrix(mo_set=mos(ispin), &
                                             density_matrix=rho_ao(ispin)%matrix)
            END DO
            CALL qs_rho_update_rho(rho, qs_env=qs_env)
            CALL qs_ks_did_change(qs_env%ks_env, rho_changed=.TRUE.)
         END IF ! do_kpoints

      CASE default
         CALL cp_abort(__LOCATION__, &
                       "Unknown interpolation method: "// &
                       TRIM(ADJUSTL(cp_to_string(wf_history%interpolation_method_nr))))
      END SELECT
      IF (PRESENT(orthogonal_wf)) orthogonal_wf = my_orthogonal_wf
      CALL cp_print_key_finished_output(io_unit, logger, qs_env%input, &
                                        "DFT%SCF%PRINT%PROGRAM_RUN_INFO")
      CALL timestop(handle)
   END SUBROUTINE wfi_extrapolate

! **************************************************************************************************
!> \brief Reorthogonalizes the wavefunctions from the previous step for k-points
!>        using the current S(k) metric and rebuilds the density matrix.
!> \param qs_env The QS environment
!> \param io_unit output unit
!> \param print_level print level
! **************************************************************************************************
   SUBROUTINE wfi_use_prev_wf_kp(qs_env, io_unit, print_level)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: io_unit, print_level

      CHARACTER(len=*), PARAMETER :: routineN = 'wfi_use_prev_wf_kp'

      COMPLEX(KIND=dp), ALLOCATABLE, DIMENSION(:)        :: col_scaling
      INTEGER                                            :: chol_info, handle, igroup, ik, ikp, &
                                                            indx, ispin, j, kplocal, nao, nkp, &
                                                            nkp_groups, nmo, nspin
      INTEGER, DIMENSION(2)                              :: kp_range
      INTEGER, DIMENSION(:, :), POINTER                  :: kp_dist
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: my_kpgrp, use_real_wfn
      REAL(KIND=dp)                                      :: eval_thresh
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigenvalues
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: xkp
      TYPE(copy_info_type), ALLOCATABLE, DIMENSION(:, :) :: info
      TYPE(cp_cfm_type)                                  :: cfm_evecs, cfm_mhalf, cfm_nao_nmo_work, &
                                                            cmos_new, csc_cfm
      TYPE(cp_cfm_type), ALLOCATABLE, DIMENSION(:)       :: csmat_cur
      TYPE(cp_fm_pool_p_type), DIMENSION(:), POINTER     :: ao_ao_fm_pools_kp
      TYPE(cp_fm_struct_type), POINTER                   :: ao_ao_struct, nmo_nmo_struct
      TYPE(cp_fm_type)                                   :: fmdummy, fmlocal
      TYPE(cp_fm_type), POINTER                          :: imos, mo_coeff, rmos
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s_kp, rho_ao_kp
      TYPE(dbcsr_type), POINTER                          :: cmatrix_db, rmatrix, tmpmat
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(qs_matrix_pools_type), POINTER                :: mpools_kp
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control

      CALL timeset(routineN, handle)

      NULLIFY (kpoints, dft_control, matrix_s_kp, scf_env, scf_control, rho, sab_nl, kp, mo_coeff, rmos, imos)

      CALL get_qs_env(qs_env, kpoints=kpoints, dft_control=dft_control, &
                      matrix_s_kp=matrix_s_kp, scf_env=scf_env, &
                      scf_control=scf_control, rho=rho)
      CALL get_kpoint_info(kpoints, nkp=nkp, xkp=xkp, use_real_wfn=use_real_wfn, &
                           kp_range=kp_range, nkp_groups=nkp_groups, kp_dist=kp_dist, &
                           sab_nl=sab_nl, cell_to_index=cell_to_index)
      kplocal = kp_range(2) - kp_range(1) + 1

      IF (use_real_wfn) THEN
         CALL timestop(handle)
         RETURN
      END IF

      kp => kpoints%kp_env(1)%kpoint_env
      nspin = SIZE(kp%mos, 2)
      CALL get_mo_set(kp%mos(1, 1), nao=nao, nmo=nmo, mo_coeff=mo_coeff)

      IF ((io_unit > 0) .AND. (print_level > low_print_level)) THEN
         WRITE (UNIT=io_unit, FMT="(/,T2,A)") &
            "Using previous wavefunctions as initial guess for k-points (with reorthogonalization)"
      END IF

      ! Allocate dbcsr work matrices
      ALLOCATE (rmatrix, cmatrix_db, tmpmat)
      CALL dbcsr_create(rmatrix, template=matrix_s_kp(1, 1)%matrix, matrix_type=dbcsr_type_symmetric)
      CALL dbcsr_create(cmatrix_db, template=matrix_s_kp(1, 1)%matrix, matrix_type=dbcsr_type_antisymmetric)
      CALL dbcsr_create(tmpmat, template=matrix_s_kp(1, 1)%matrix, matrix_type=dbcsr_type_no_symmetry)
      CALL cp_dbcsr_alloc_block_from_nbl(rmatrix, sab_nl)
      CALL cp_dbcsr_alloc_block_from_nbl(cmatrix_db, sab_nl)

      CALL get_kpoint_info(kpoints, mpools=mpools_kp)
      CALL mpools_get(mpools_kp, ao_ao_fm_pools=ao_ao_fm_pools_kp)
      CALL fm_pool_create_fm(ao_ao_fm_pools_kp(1)%pool, fmlocal)
      CALL cp_fm_get_info(fmlocal, matrix_struct=ao_ao_struct)

      CALL cp_cfm_create(cmos_new, mo_coeff%matrix_struct)
      CALL cp_cfm_create(cfm_nao_nmo_work, mo_coeff%matrix_struct)

      NULLIFY (nmo_nmo_struct)
      CALL cp_fm_struct_create(nmo_nmo_struct, template_fmstruct=mo_coeff%matrix_struct, &
                               nrow_global=nmo, ncol_global=nmo)
      CALL cp_cfm_create(csc_cfm, nmo_nmo_struct)
      CALL cp_fm_struct_release(nmo_nmo_struct)

      para_env => kpoints%blacs_env_all%para_env
      ALLOCATE (info(kplocal*nkp_groups, 2))

      ALLOCATE (csmat_cur(kplocal))
      DO ikp = 1, kplocal
         CALL cp_cfm_create(csmat_cur(ikp), ao_ao_struct)
      END DO

      ! Phase A: Fourier Transform S(R) -> S(k)
      indx = 0
      DO ikp = 1, kplocal
         DO igroup = 1, nkp_groups
            ik = kp_dist(1, igroup) + ikp - 1
            my_kpgrp = (ik >= kp_range(1) .AND. ik <= kp_range(2))
            indx = indx + 1

            CALL dbcsr_set(rmatrix, 0.0_dp)
            CALL dbcsr_set(cmatrix_db, 0.0_dp)
            CALL rskp_transform(rmatrix=rmatrix, cmatrix=cmatrix_db, rsmat=matrix_s_kp, &
                                ispin=1, xkp=xkp(1:3, ik), cell_to_index=cell_to_index, sab_nl=sab_nl)
            CALL dbcsr_desymmetrize(rmatrix, tmpmat)
            CALL copy_dbcsr_to_fm(tmpmat, scf_env%scf_work1(1))
            CALL dbcsr_desymmetrize(cmatrix_db, tmpmat)
            CALL copy_dbcsr_to_fm(tmpmat, scf_env%scf_work1(2))

            IF (my_kpgrp) THEN
               CALL cp_fm_start_copy_general(scf_env%scf_work1(1), fmlocal, para_env, info(indx, 1))
               CALL cp_fm_start_copy_general(scf_env%scf_work1(2), fmlocal, para_env, info(indx, 2))
            ELSE
               CALL cp_fm_start_copy_general(scf_env%scf_work1(1), fmdummy, para_env, info(indx, 1))
               CALL cp_fm_start_copy_general(scf_env%scf_work1(2), fmdummy, para_env, info(indx, 2))
            END IF
         END DO
      END DO

      ! Finish Communication
      indx = 0
      DO ikp = 1, kplocal
         DO igroup = 1, nkp_groups
            ik = kp_dist(1, igroup) + ikp - 1
            my_kpgrp = (ik >= kp_range(1) .AND. ik <= kp_range(2))
            indx = indx + 1
            IF (my_kpgrp) THEN
               CALL cp_fm_finish_copy_general(fmlocal, info(indx, 1))
               CALL cp_cfm_scale_and_add_fm(z_zero, csmat_cur(ikp), z_one, fmlocal)
               CALL cp_fm_finish_copy_general(fmlocal, info(indx, 2))
               CALL cp_cfm_scale_and_add_fm(z_one, csmat_cur(ikp), gaussi, fmlocal)
            END IF
         END DO
      END DO

      DO indx = 1, kplocal*nkp_groups
         CALL cp_fm_cleanup_copy_general(info(indx, 1))
         CALL cp_fm_cleanup_copy_general(info(indx, 2))
      END DO

      ! Phase B: Orthogonalize current MOs w.r.t. new S(k)
      ALLOCATE (eigenvalues(nmo))
      eval_thresh = 1.0E-12_dp

      DO ikp = 1, kplocal
         kp => kpoints%kp_env(ikp)%kpoint_env
         DO ispin = 1, nspin
            CALL get_mo_set(kp%mos(1, ispin), mo_coeff=rmos)
            CALL get_mo_set(kp%mos(2, ispin), mo_coeff=imos)
            CALL cp_fm_to_cfm(rmos, imos, cmos_new)

            CALL cp_cfm_gemm('N', 'N', nao, nmo, nao, z_one, &
                             csmat_cur(ikp), cmos_new, z_zero, cfm_nao_nmo_work)
            CALL cp_cfm_gemm('C', 'N', nmo, nmo, nao, z_one, &
                             cmos_new, cfm_nao_nmo_work, z_zero, csc_cfm)

            CALL cp_cfm_cholesky_decompose(csc_cfm, info_out=chol_info)
            IF (chol_info == 0) THEN
               CALL cp_cfm_triangular_multiply(csc_cfm, cmos_new, side='R', invert_tr=.TRUE., uplo_tr='U')
            ELSE
               CALL cp_cfm_gemm('C', 'N', nmo, nmo, nao, z_one, cmos_new, cfm_nao_nmo_work, z_zero, csc_cfm)
               CALL cp_cfm_create(cfm_evecs, csc_cfm%matrix_struct)
               CALL cp_cfm_create(cfm_mhalf, csc_cfm%matrix_struct)
               CALL cp_cfm_heevd(csc_cfm, cfm_evecs, eigenvalues)
               CALL cp_cfm_to_cfm(cfm_evecs, cfm_mhalf)
               ALLOCATE (col_scaling(nmo))
               DO j = 1, nmo
                  IF (eigenvalues(j) > eval_thresh) THEN
                     col_scaling(j) = CMPLX(1.0_dp/SQRT(eigenvalues(j)), 0.0_dp, KIND=dp)
                  ELSE
                     col_scaling(j) = z_zero
                  END IF
               END DO
               CALL cp_cfm_column_scale(cfm_mhalf, col_scaling)
               DEALLOCATE (col_scaling)
               CALL cp_cfm_gemm('N', 'C', nmo, nmo, nmo, z_one, cfm_mhalf, cfm_evecs, z_zero, csc_cfm)
               CALL cp_cfm_gemm('N', 'N', nao, nmo, nmo, z_one, cmos_new, csc_cfm, z_zero, cfm_nao_nmo_work)
               CALL cp_cfm_to_cfm(cfm_nao_nmo_work, cmos_new)
               CALL cp_cfm_release(cfm_evecs)
               CALL cp_cfm_release(cfm_mhalf)
            END IF
            CALL cp_cfm_to_fm(cmos_new, rmos, imos)
         END DO
      END DO
      DEALLOCATE (eigenvalues)

      ! Phase C: Rebuild Density Matrix P(R)
      CALL kpoint_set_mo_occupation(kpoints, scf_control%smear)
      CALL kpoint_density_matrices(kpoints)
      CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp)
      CALL kpoint_density_transform(kpoints, rho_ao_kp, .FALSE., &
                                    matrix_s_kp(1, 1)%matrix, sab_nl, scf_env%scf_work1)
      CALL qs_rho_update_rho(rho, qs_env=qs_env)
      CALL qs_ks_did_change(qs_env%ks_env, rho_changed=.TRUE.)

      ! Cleanup
      DO ikp = 1, kplocal
         CALL cp_cfm_release(csmat_cur(ikp))
      END DO
      DEALLOCATE (csmat_cur)
      DEALLOCATE (info)
      CALL cp_cfm_release(cmos_new)
      CALL cp_cfm_release(cfm_nao_nmo_work)
      CALL cp_cfm_release(csc_cfm)
      CALL fm_pool_give_back_fm(ao_ao_fm_pools_kp(1)%pool, fmlocal)
      CALL dbcsr_deallocate_matrix(rmatrix)
      CALL dbcsr_deallocate_matrix(cmatrix_db)
      CALL dbcsr_deallocate_matrix(tmpmat)

      CALL timestop(handle)
   END SUBROUTINE wfi_use_prev_wf_kp

! **************************************************************************************************
!> \brief Performs ASPC wavefunction extrapolation for k-point calculations.
!>        Applies ASPC coefficients to complex MO coefficients at each k-point,
!>        with subspace alignment via historical overlap matrices.
!>        Delegates final orthogonalization and density building to wfi_use_prev_wf_kp.
!> \param wf_history  wavefunction history buffer
!> \param qs_env      QS environment
!> \param nvec        number of history snapshots to use
!> \param io_unit     output unit for logging
!> \param print_level current print level
! **************************************************************************************************
   SUBROUTINE wfi_extrapolate_aspc_kp(wf_history, qs_env, nvec, io_unit, print_level)
      TYPE(qs_wf_history_type), POINTER                  :: wf_history
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: nvec, io_unit, print_level

      CHARACTER(len=*), PARAMETER :: routineN = 'wfi_extrapolate_aspc_kp'

      INTEGER                                            :: handle, i, ikp, ispin, kplocal, nao, &
                                                            nmo, nspin
      INTEGER, DIMENSION(2)                              :: kp_range
      LOGICAL                                            :: use_real_wfn
      REAL(KIND=dp)                                      :: alpha_coeff
      TYPE(cp_cfm_type)                                  :: cfm_nao_nmo_work, cmos_1, cmos_i, &
                                                            cmos_new, csc_cfm
      TYPE(cp_fm_struct_type), POINTER                   :: nmo_nmo_struct
      TYPE(cp_fm_type), POINTER                          :: imos, mo_coeff, rmos
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(qs_wf_snapshot_type), POINTER                 :: t0_state, t1_state

      CALL timeset(routineN, handle)
      NULLIFY (kpoints, kp, mo_coeff, rmos, imos, t0_state, t1_state, nmo_nmo_struct)

      CALL get_qs_env(qs_env, kpoints=kpoints)
      CALL get_kpoint_info(kpoints, use_real_wfn=use_real_wfn, kp_range=kp_range)
      kplocal = kp_range(2) - kp_range(1) + 1

      IF (use_real_wfn) THEN
         CPWARN("ASPC with k-points requires complex wavefunctions; falling back.")
         CALL timestop(handle)
         RETURN
      END IF

      kp => kpoints%kp_env(1)%kpoint_env
      nspin = SIZE(kp%mos, 2)
      CALL get_mo_set(kp%mos(1, 1), nao=nao, nmo=nmo, mo_coeff=mo_coeff)

      IF ((io_unit > 0) .AND. (print_level > low_print_level)) THEN
         WRITE (UNIT=io_unit, FMT="(/,T2,A,/,T3,A,I0)") &
            "ASPC wavefunction extrapolation for k-points", &
            "ASPC order: ", MAX(nvec - 2, 0)
      END IF

      CALL cp_cfm_create(cmos_new, mo_coeff%matrix_struct)
      CALL cp_cfm_create(cmos_1, mo_coeff%matrix_struct)
      CALL cp_cfm_create(cmos_i, mo_coeff%matrix_struct)
      CALL cp_cfm_create(cfm_nao_nmo_work, mo_coeff%matrix_struct)

      CALL cp_fm_struct_create(nmo_nmo_struct, template_fmstruct=mo_coeff%matrix_struct, &
                               nrow_global=nmo, ncol_global=nmo)
      CALL cp_cfm_create(csc_cfm, nmo_nmo_struct)
      CALL cp_fm_struct_release(nmo_nmo_struct)

      ! =======================================================================
      ! Phase 1: Initialize C_new(k) = B(1) * C_1(k)
      ! =======================================================================
      t1_state => wfi_get_snapshot(wf_history, wf_index=1)
      alpha_coeff = REAL(4*nvec - 2, KIND=dp)/REAL(nvec + 1, KIND=dp)

      IF ((io_unit > 0) .AND. (print_level > low_print_level)) THEN
         WRITE (UNIT=io_unit, FMT="(T3,A2,I0,A4,F10.6)") "B(", 1, ") = ", alpha_coeff
      END IF

      DO ikp = 1, kplocal
         kp => kpoints%kp_env(ikp)%kpoint_env
         DO ispin = 1, nspin
            CALL get_mo_set(kp%mos(1, ispin), mo_coeff=rmos)
            CALL get_mo_set(kp%mos(2, ispin), mo_coeff=imos)
            CALL cp_fm_to_fm(t1_state%wf_kp(ikp, 1, ispin), rmos)
            CALL cp_fm_to_fm(t1_state%wf_kp(ikp, 2, ispin), imos)
            CALL cp_fm_scale(alpha_coeff, rmos)
            CALL cp_fm_scale(alpha_coeff, imos)
         END DO
      END DO

      ! =======================================================================
      ! Phase 2: Accumulate historical snapshots C_new += B(i) * C_proj(k)
      ! =======================================================================
      DO i = 2, nvec
         t0_state => wfi_get_snapshot(wf_history, wf_index=i)
         alpha_coeff = (-1.0_dp)**(i + 1)*REAL(i, KIND=dp)* &
                       binomial(2*nvec, nvec - i)/binomial(2*nvec - 2, nvec - 1)

         IF ((io_unit > 0) .AND. (print_level > low_print_level)) THEN
            WRITE (UNIT=io_unit, FMT="(T3,A2,I0,A4,F10.6)") "B(", i, ") = ", alpha_coeff
         END IF

         DO ikp = 1, kplocal
            kp => kpoints%kp_env(ikp)%kpoint_env
            DO ispin = 1, nspin
               CALL cp_fm_to_cfm(t1_state%wf_kp(ikp, 1, ispin), t1_state%wf_kp(ikp, 2, ispin), cmos_1)
               CALL cp_fm_to_cfm(t0_state%wf_kp(ikp, 1, ispin), t0_state%wf_kp(ikp, 2, ispin), cmos_i)

               ! Subspace projection: C_proj = C_i * (C_i^dag S_i C_1)
               CALL cp_cfm_gemm('N', 'N', nao, nmo, nao, z_one, &
                                t0_state%overlap_cfm_kp(ikp), cmos_1, z_zero, cfm_nao_nmo_work)
               CALL cp_cfm_gemm('C', 'N', nmo, nmo, nao, z_one, &
                                cmos_i, cfm_nao_nmo_work, z_zero, csc_cfm)
               CALL cp_cfm_gemm('N', 'N', nao, nmo, nmo, z_one, &
                                cmos_i, csc_cfm, z_zero, cfm_nao_nmo_work)

               CALL get_mo_set(kp%mos(1, ispin), mo_coeff=rmos)
               CALL get_mo_set(kp%mos(2, ispin), mo_coeff=imos)
               CALL cp_fm_to_cfm(rmos, imos, cmos_new)
               CALL cp_cfm_scale_and_add(z_one, cmos_new, CMPLX(alpha_coeff, 0.0_dp, KIND=dp), cfm_nao_nmo_work)
               CALL cp_cfm_to_fm(cmos_new, rmos, imos)
            END DO
         END DO
      END DO

      CALL cp_cfm_release(cmos_new)
      CALL cp_cfm_release(cmos_1)
      CALL cp_cfm_release(cmos_i)
      CALL cp_cfm_release(cfm_nao_nmo_work)
      CALL cp_cfm_release(csc_cfm)

      ! =======================================================================
      ! Phase 3: Delegate Orthogonalization and Density Building
      ! =======================================================================
      ! We pass io_unit = 0 to suppress the redundant "Using previous..." print
      ! inside wfi_use_prev_wf_kp, keeping the ASPC log output perfectly clean.
      CALL wfi_use_prev_wf_kp(qs_env, 0, print_level)

      CALL timestop(handle)

   END SUBROUTINE wfi_extrapolate_aspc_kp

! **************************************************************************************************
!> \brief Decides if scf control variables has to changed due
!>      to using a WF extrapolation.
!> \param qs_env The QS environment
!> \param nvec ...
!> \par History
!>      11.2006 created [TdK]
!> \author Thomas D. Kuehne (tkuehne@phys.chem.ethz.ch)
! **************************************************************************************************
   ELEMENTAL SUBROUTINE wfi_set_history_variables(qs_env, nvec)
      TYPE(qs_environment_type), INTENT(INOUT)           :: qs_env
      INTEGER, INTENT(IN)                                :: nvec

      IF (nvec >= qs_env%wf_history%memory_depth) THEN
         IF ((qs_env%scf_control%max_scf_hist /= 0) .AND. (qs_env%scf_control%eps_scf_hist /= 0)) THEN
            qs_env%scf_control%max_scf = qs_env%scf_control%max_scf_hist
            qs_env%scf_control%eps_scf = qs_env%scf_control%eps_scf_hist
            qs_env%scf_control%outer_scf%have_scf = .FALSE.
         ELSE IF (qs_env%scf_control%max_scf_hist /= 0) THEN
            qs_env%scf_control%max_scf = qs_env%scf_control%max_scf_hist
            qs_env%scf_control%outer_scf%have_scf = .FALSE.
         ELSE IF (qs_env%scf_control%eps_scf_hist /= 0) THEN
            qs_env%scf_control%eps_scf = qs_env%scf_control%eps_scf_hist
            qs_env%scf_control%outer_scf%eps_scf = qs_env%scf_control%eps_scf_hist
         END IF
      END IF

   END SUBROUTINE wfi_set_history_variables

! **************************************************************************************************
!> \brief updates the snapshot buffer, taking a new snapshot
!> \param wf_history the history buffer to update
!> \param qs_env the qs_env we get the info from
!> \param dt ...
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE wfi_update(wf_history, qs_env, dt)
      TYPE(qs_wf_history_type), POINTER                  :: wf_history
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(KIND=dp), INTENT(in)                          :: dt

      CPASSERT(ASSOCIATED(wf_history))
      CPASSERT(wf_history%ref_count > 0)
      CPASSERT(ASSOCIATED(qs_env))

      wf_history%snapshot_count = wf_history%snapshot_count + 1
      IF (wf_history%memory_depth > 0) THEN
         wf_history%last_state_index = MODULO(wf_history%snapshot_count, &
                                              wf_history%memory_depth) + 1
         CALL wfs_update(snapshot=wf_history%past_states &
                         (wf_history%last_state_index)%snapshot, wf_history=wf_history, &
                         qs_env=qs_env, dt=dt)
      END IF
   END SUBROUTINE wfi_update

! **************************************************************************************************
!> \brief reorthogonalizes the mos
!> \param qs_env the qs_env in which to orthogonalize
!> \param v_matrix the vectors to orthogonalize
!> \param n_col number of column of v to orthogonalize
!> \par History
!>      04.2003 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE reorthogonalize_vectors(qs_env, v_matrix, n_col)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), INTENT(IN)                       :: v_matrix
      INTEGER, INTENT(in), OPTIONAL                      :: n_col

      CHARACTER(len=*), PARAMETER :: routineN = 'reorthogonalize_vectors'

      INTEGER                                            :: handle, my_n_col
      LOGICAL                                            :: has_unit_metric, &
                                                            ortho_contains_cholesky, &
                                                            smearing_is_used
      TYPE(cp_fm_pool_type), POINTER                     :: maxao_maxmo_fm_pool
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_matrix_pools_type), POINTER                :: mpools
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control

      NULLIFY (scf_env, scf_control, maxao_maxmo_fm_pool, matrix_s, mpools, dft_control)
      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(qs_env))

      CALL cp_fm_get_info(v_matrix, ncol_global=my_n_col)
      IF (PRESENT(n_col)) my_n_col = n_col
      CALL get_qs_env(qs_env, mpools=mpools, &
                      scf_env=scf_env, &
                      scf_control=scf_control, &
                      matrix_s=matrix_s, &
                      dft_control=dft_control)
      CALL mpools_get(mpools, maxao_maxmo_fm_pool=maxao_maxmo_fm_pool)
      IF (ASSOCIATED(scf_env)) THEN
         ortho_contains_cholesky = (scf_env%method /= ot_method_nr) .AND. &
                                   (scf_env%cholesky_method > 0) .AND. &
                                   ASSOCIATED(scf_env%ortho)
      ELSE
         ortho_contains_cholesky = .FALSE.
      END IF

      CALL get_qs_env(qs_env, has_unit_metric=has_unit_metric)
      smearing_is_used = .FALSE.
      IF (dft_control%smear) THEN
         smearing_is_used = .TRUE.
      END IF

      IF (has_unit_metric) THEN
         CALL make_basis_simple(v_matrix, my_n_col)
      ELSE IF (smearing_is_used) THEN
         CALL make_basis_lowdin(vmatrix=v_matrix, ncol=my_n_col, &
                                matrix_s=matrix_s(1)%matrix)
      ELSE IF (ortho_contains_cholesky) THEN
         CALL make_basis_cholesky(vmatrix=v_matrix, ncol=my_n_col, &
                                  ortho=scf_env%ortho)
      ELSE
         CALL make_basis_sm(v_matrix, my_n_col, matrix_s(1)%matrix)
      END IF
      CALL timestop(handle)
   END SUBROUTINE reorthogonalize_vectors

! **************************************************************************************************
!> \brief purges wf_history retaining only the latest snapshot
!> \param qs_env the qs env with the latest result, and that will contain
!>        the purged wf_history
!> \par History
!>      05.2016 created [Nico Holmberg]
!> \author Nico Holmberg
! **************************************************************************************************
   SUBROUTINE wfi_purge_history(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'wfi_purge_history'

      INTEGER                                            :: handle, io_unit, print_level
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_wf_history_type), POINTER                  :: wf_history

      NULLIFY (dft_control, wf_history)

      CALL timeset(routineN, handle)
      logger => cp_get_default_logger()
      print_level = logger%iter_info%print_level
      io_unit = cp_print_key_unit_nr(logger, qs_env%input, "DFT%SCF%PRINT%PROGRAM_RUN_INFO", &
                                     extension=".scfLog")

      CPASSERT(ASSOCIATED(qs_env))
      CPASSERT(ASSOCIATED(qs_env%wf_history))
      CPASSERT(qs_env%wf_history%ref_count > 0)
      CALL get_qs_env(qs_env, dft_control=dft_control)

      SELECT CASE (qs_env%wf_history%interpolation_method_nr)
      CASE (wfi_use_guess_method_nr, wfi_use_prev_wf_method_nr, &
            wfi_use_prev_p_method_nr, wfi_use_prev_rho_r_method_nr, &
            wfi_frozen_method_nr)
         ! do nothing
      CASE (wfi_linear_wf_method_nr, wfi_linear_p_method_nr, &
            wfi_linear_ps_method_nr, wfi_ps_method_nr, &
            wfi_aspc_nr)
         IF (qs_env%wf_history%snapshot_count >= 2) THEN
            IF (debug_this_module .AND. io_unit > 0) &
               WRITE (io_unit, FMT="(T2,A)") "QS| Purging WFN history"
            CALL wfi_create(wf_history, interpolation_method_nr= &
                            dft_control%qs_control%wf_interpolation_method_nr, &
                            extrapolation_order=dft_control%qs_control%wf_extrapolation_order, &
                            has_unit_metric=qs_env%has_unit_metric)
            CALL set_qs_env(qs_env=qs_env, &
                            wf_history=wf_history)
            CALL wfi_release(wf_history)
            CALL wfi_update(qs_env%wf_history, qs_env=qs_env, dt=1.0_dp)
         END IF
      CASE DEFAULT
         CPABORT("Unknown extrapolation method.")
      END SELECT
      CALL timestop(handle)

   END SUBROUTINE wfi_purge_history

END MODULE qs_wf_history_methods
