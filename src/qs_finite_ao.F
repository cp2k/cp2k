!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Shared finite-system AO position and SOC operators for post-SCF analyses.
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
MODULE qs_finite_ao
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_create,&
                                              dbcsr_p_type,&
                                              dbcsr_set
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE input_constants,                 ONLY: do_method_gapw,&
                                              do_method_gpw
   USE kinds,                           ONLY: dp
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_moments,                      ONLY: build_local_moment_matrix
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE soc_pseudopotential_methods,     ONLY: V_SOC_xyz_from_pseudopotential
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_finite_ao'
   PUBLIC :: finite_ao_positions, finite_ao_soc
CONTAINS

! **************************************************************************************************
!> \brief Analytic covariant first moments in the coordinate frame of SUBSYS/COORD.
!> \param qs_env Finite Quickstep environment
!> \param positions Allocated X, Y, Z matrices, released by caller
! **************************************************************************************************
   SUBROUTINE finite_ao_positions(qs_env, positions)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), POINTER                        :: positions(:)

      INTEGER                                            :: axis
      LOGICAL                                            :: do_kpoints
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), POINTER                        :: matrix_s(:)
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(neighbor_list_set_p_type), POINTER            :: sab_orb(:)

      CALL get_qs_env(qs_env, cell=cell, dft_control=dft_control, do_kpoints=do_kpoints, &
                      matrix_s=matrix_s, sab_orb=sab_orb)
      IF (ANY(cell%perd /= 0) .OR. do_kpoints) THEN
         CPABORT('Finite AO positions require PERIODIC NONE without KPOINTS.')
      END IF
      SELECT CASE (dft_control%qs_control%method_id)
      CASE (do_method_gpw, do_method_gapw)
      CASE DEFAULT
         CPABORT('Finite AO positions currently support GPW/GAPW only.')
      END SELECT
      ALLOCATE (positions(3))
      DO axis = 1, 3
         ALLOCATE (positions(axis)%matrix)
         CALL dbcsr_create(positions(axis)%matrix, template=matrix_s(1)%matrix, name='AO POSITION')
         CALL cp_dbcsr_alloc_block_from_nbl(positions(axis)%matrix, sab_orb)
         CALL dbcsr_set(positions(axis)%matrix, 0.0_dp)
      END DO
      CALL build_local_moment_matrix(qs_env, positions, 1, ref_point=[0.0_dp, 0.0_dp, 0.0_dp])
   END SUBROUTINE finite_ao_positions
! **************************************************************************************************
!> \brief Validate and construct post-SCF GTH SOC components shared by finite AO analyses.
!> \param qs_env Restricted SCF environment with GTH SOC potentials
!> \param soc_matrices Newly allocated SOC matrices, released by caller
! **************************************************************************************************
   SUBROUTINE finite_ao_soc(qs_env, soc_matrices)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), POINTER                        :: soc_matrices(:, :)

      INTEGER                                            :: ik
      TYPE(dbcsr_p_type), POINTER                        :: matrix_ks(:, :)
      TYPE(neighbor_list_set_p_type), POINTER            :: sap_ppnl(:)
      TYPE(qs_kind_type), POINTER                        :: qs_kind_set(:)

      CALL get_qs_env(qs_env, matrix_ks_kp=matrix_ks, qs_kind_set=qs_kind_set, sap_ppnl=sap_ppnl)
      IF (SIZE(matrix_ks, 1) /= 1) CPABORT('AO SOC requires a restricted SCF calculation.')
      IF (.NOT. ASSOCIATED(sap_ppnl)) CPABORT('Finite AO SOC requires nonlocal GTH SOC projectors.')
      DO ik = 1, SIZE(qs_kind_set)
         IF (.NOT. ASSOCIATED(qs_kind_set(ik)%gth_potential)) THEN
            CPABORT('Finite AO SOC requires GTH potentials.')
         END IF
         IF (.NOT. qs_kind_set(ik)%gth_potential%soc) CPABORT('GTH potential lacks SOC parameters.')
      END DO
      CALL V_SOC_xyz_from_pseudopotential(qs_env, soc_matrices)
   END SUBROUTINE finite_ao_soc
END MODULE qs_finite_ao
