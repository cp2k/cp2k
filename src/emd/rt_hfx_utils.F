!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Utility functions that are needed for RTP/EMD in combination with
!>        HF or hybrid functionals (needs to deal with imaginary KS and P
!> \par History
!>      2014 created [fschiff]
!> \author Florina Schiffmann
! **************************************************************************************************
MODULE rt_hfx_utils
   USE admm_types,                      ONLY: get_admm_env,&
                                              set_admm_env
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE dbcsr_api,                       ONLY: dbcsr_create,&
                                              dbcsr_p_type,&
                                              dbcsr_set,&
                                              dbcsr_type,&
                                              dbcsr_type_antisymmetric
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_ks_types,                     ONLY: qs_ks_env_type,&
                                              set_ks_env
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_set,&
                                              qs_rho_type
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rt_hfx_utils'

   PUBLIC :: rtp_hfx_rebuild

!***
CONTAINS
! **************************************************************************************************
!> \brief rebuilds the structures of P and KS (imaginary) in case S changed
!> \param qs_env ...
!> \author Florian Schiffmann
! **************************************************************************************************
   SUBROUTINE rtp_hfx_rebuild(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks_aux_im, matrix_s_aux, &
                                                            rho_aux_ao_im
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_aux
      TYPE(qs_rho_type), POINTER                         :: rho_aux

      NULLIFY (dft_control)
      NULLIFY (sab_aux, rho_aux, rho_aux_ao_im, matrix_ks_aux_im, matrix_s_aux)

      CALL get_qs_env(qs_env, &
                      dft_control=dft_control)

      IF (dft_control%do_admm) THEN
         CALL get_admm_env(qs_env%admm_env, &
                           matrix_s_aux_fit=matrix_s_aux, &
                           sab_aux_fit=sab_aux, &
                           rho_aux_fit=rho_aux, &
                           matrix_ks_aux_fit_im=matrix_ks_aux_im)
         CALL qs_rho_get(rho_aux, rho_ao_im=rho_aux_ao_im)
         CALL rebuild_matrices(rho_aux_ao_im, matrix_ks_aux_im, sab_aux, matrix_s_aux, &
                               dft_control%nspins)
         CALL set_admm_env(qs_env%admm_env, matrix_ks_aux_fit_im=matrix_ks_aux_im)
         CALL qs_rho_set(rho_aux, rho_ao_im=rho_aux_ao_im)
      END IF

   END SUBROUTINE rtp_hfx_rebuild

! **************************************************************************************************
!> \brief does the actual rebuilding of P and KS (imaginary) in case S changed
!> \param matrix_p ...
!> \param matrix_ks ...
!> \param sab_orb ...
!> \param matrix_s ...
!> \param nspins ...
!> \author Florian Schiffmann
! **************************************************************************************************

   SUBROUTINE rebuild_matrices(matrix_p, matrix_ks, sab_orb, matrix_s, nspins)
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_p, matrix_ks
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      INTEGER, INTENT(in)                                :: nspins

      INTEGER                                            :: i

      IF (ASSOCIATED(matrix_p)) THEN
         CALL dbcsr_deallocate_matrix_set(matrix_p)
      END IF
      ! Create a new density matrix set
      CALL dbcsr_allocate_matrix_set(matrix_p, nspins)
      DO i = 1, nspins
         ALLOCATE (matrix_p(i)%matrix)
         CALL dbcsr_create(matrix=matrix_p(i)%matrix, template=matrix_s(1)%matrix, &
                           name="Imaginary density matrix", matrix_type=dbcsr_type_antisymmetric, nze=0)
         CALL cp_dbcsr_alloc_block_from_nbl(matrix_p(i)%matrix, sab_orb)
         CALL dbcsr_set(matrix_p(i)%matrix, 0.0_dp)
      END DO

      IF (ASSOCIATED(matrix_ks)) THEN
         CALL dbcsr_deallocate_matrix_set(matrix_ks)
      END IF
      ! Create a new density matrix set
      CALL dbcsr_allocate_matrix_set(matrix_ks, nspins)
      DO i = 1, nspins
         ALLOCATE (matrix_ks(i)%matrix)
         CALL dbcsr_create(matrix=matrix_ks(i)%matrix, template=matrix_s(1)%matrix, &
                           name="Imaginary Kohn-Sham matrix", matrix_type=dbcsr_type_antisymmetric, nze=0)
         CALL cp_dbcsr_alloc_block_from_nbl(matrix_ks(i)%matrix, sab_orb)
         CALL dbcsr_set(matrix_ks(i)%matrix, 0.0_dp)
      END DO

   END SUBROUTINE rebuild_matrices

! **************************************************************************************************
!> \brief (Re)-allocate imaginary parts of density matrix and KS-matrix as kpoint_transitional_types
!> \param qs_env ...
!> \param rho ...
!> \param ks_env ...
!> \param nspins ...
!> \param nimages ...
!> \param rebuild_matrix_p ...
!> \param rebuild_matrix_h ...
! **************************************************************************************************
   SUBROUTINE rebuild_matrices_im_kp(qs_env, rho, ks_env, nspins, nimages, rebuild_matrix_p, rebuild_matrix_h)
      TYPE(qs_environment_type), INTENT(inout), POINTER  :: qs_env
      TYPE(qs_rho_type), INTENT(INOUT), POINTER          :: rho
      TYPE(qs_ks_env_type), INTENT(INOUT), POINTER       :: ks_env
      INTEGER, INTENT(in), OPTIONAL                      :: nspins, nimages
      LOGICAL, INTENT(in), OPTIONAL                      :: rebuild_matrix_p, rebuild_matrix_h

      CHARACTER(LEN=default_string_length)               :: headline
      INTEGER                                            :: image, ispin
      LOGICAL                                            :: my_rebuild_h, my_rebuild_p
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_h_im, matrix_h_kp, matrix_ks_kp, &
                                                            matrix_p_im, rho_ao_kp
      TYPE(dbcsr_type), POINTER                          :: refmatrix
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb

      my_rebuild_p = .FALSE.
      IF (PRESENT(rebuild_matrix_p)) my_rebuild_p = rebuild_matrix_p
      my_rebuild_h = .FALSE.
      IF (PRESENT(rebuild_matrix_h)) my_rebuild_h = rebuild_matrix_h

      NULLIFY (matrix_ks_kp, rho_ao_kp, sab_orb, matrix_h_kp)
      CALL get_qs_env(qs_env, matrix_ks_kp=matrix_ks_kp, matrix_h_kp=matrix_h_kp, sab_orb=sab_orb)
      CALL qs_rho_get(rho, rho_ao_kp=rho_ao_kp)

      IF (my_rebuild_p) THEN
         NULLIFY (matrix_p_im)
         CALL qs_rho_get(rho, rho_ao_im_kp=matrix_p_im)

         CALL dbcsr_allocate_matrix_set(matrix_p_im, nspins, nimages)
         refmatrix => rho_ao_kp(1, 1)%matrix
         DO ispin = 1, nspins
            DO image = 1, nimages
               IF (nspins > 1) THEN
                  IF (ispin == 1) THEN
                     headline = "IMAGINARY PART OF DENSITY MATRIX FOR ALPHA SPIN"
                  ELSE
                     headline = "IMAGINARY PART OF DENSITY MATRIX FOR BETA SPIN"
                  END IF
               ELSE
                  headline = "IMAGINARY PART OF DENSITY MATRIX"
               END IF
               ALLOCATE (matrix_p_im(ispin, image)%matrix)
               CALL dbcsr_create(matrix=matrix_p_im(ispin, image)%matrix, template=refmatrix, &
                                 name=TRIM(headline), matrix_type=dbcsr_type_antisymmetric, nze=0)
               CALL cp_dbcsr_alloc_block_from_nbl(matrix_p_im(ispin, image)%matrix, sab_orb)
               CALL dbcsr_set(matrix_p_im(ispin, image)%matrix, 0.0_dp)
            END DO
         END DO
         CALL qs_rho_set(rho, rho_ao_im_kp=matrix_p_im)
      END IF

      IF (my_rebuild_h) THEN
         NULLIFY (matrix_h_im)
         CALL get_qs_env(qs_env, matrix_h_im_kp=matrix_h_im)

         CALL dbcsr_allocate_matrix_set(matrix_h_im, 1, nimages)
         refmatrix => matrix_h_kp(1, 1)%matrix

         DO image = 1, nimages
            headline = "IMAGINARY CORE HAMILTONIAN MATRIX"
            ALLOCATE (matrix_h_im(1, image)%matrix)
            CALL dbcsr_create(matrix=matrix_h_im(1, image)%matrix, template=refmatrix, &
                              name=TRIM(headline), matrix_type=dbcsr_type_antisymmetric, nze=0)
            CALL cp_dbcsr_alloc_block_from_nbl(matrix_h_im(1, image)%matrix, sab_orb)
            CALL dbcsr_set(matrix_h_im(1, image)%matrix, 0.0_dp)
         END DO
         CALL set_ks_env(ks_env, matrix_h_im_kp=matrix_h_im)
      END IF

   END SUBROUTINE

END MODULE rt_hfx_utils
