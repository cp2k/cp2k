!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Types needed for cubic-scaling RPA and SOS-Laplace-MP2 forces
!> \author Augustin Bussy
! **************************************************************************************************
MODULE rpa_im_time_force_types
   USE cp_dbcsr_operations,             ONLY: dbcsr_deallocate_matrix_set
   USE dbcsr_api,                       ONLY: dbcsr_p_type,&
                                              dbcsr_release,&
                                              dbcsr_type
   USE dbt_api,                         ONLY: dbt_destroy,&
                                              dbt_type
   USE hfx_types,                       ONLY: block_ind_type,&
                                              dealloc_containers,&
                                              hfx_compression_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type,&
                                              release_neighbor_list_sets
   USE qs_tensors,                      ONLY: neighbor_list_3c_destroy
   USE qs_tensors_types,                ONLY: neighbor_list_3c_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'rpa_im_time_force_types'

   PUBLIC :: im_time_force_type, im_time_force_release

   !Type definitions
   TYPE im_time_force_type

      !The various 2-center integral derivatives
      TYPE(dbt_type), DIMENSION(3)  :: t_2c_der_metric, t_2c_der_pot, t_2c_der_ovlp

      !The 3-center integral derivatives (with the RI metric operator)
      TYPE(dbt_type), DIMENSION(3)  :: t_3c_der_AO, & ! (RI| AO deriv_AO)
                                       t_3c_der_RI    ! (deriv_RI| AO AO)

      !The compressed 3-center derivatives
      TYPE(hfx_compression_type), ALLOCATABLE, DIMENSION(:, :)   :: t_3c_der_AO_comp, t_3c_der_RI_comp
      TYPE(block_ind_type), ALLOCATABLE, DIMENSION(:, :)         :: t_3c_der_AO_ind, t_3c_der_RI_ind

      !The RI related 2-center quantities
      TYPE(dbt_type) :: t_2c_pot_psqrt, t_2c_inv_metric, t_2c_K, t_2c_pot_msqrt
      TYPE(dbcsr_type) :: inv_ovlp, G_PQ

      !The occupied and virtual density matrices (standard block size, one for each spin)
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER :: P_occ => NULL(), P_virt => NULL()

      !The weighted sum of the O(tau) matrices for thre response
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER :: sum_O_tau => NULL()

      !The weigthed sum of the YP matrices for the trace with the Fockian derivative
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER :: sum_YP_tau => NULL()

      !Split block size info
      INTEGER, DIMENSION(:), ALLOCATABLE  :: bsizes_RI_split, bsizes_AO_split

      !Keep track of atom index for splitted blocks
      INTEGER, DIMENSION(:), ALLOCATABLE  :: idx_to_at_AO, idx_to_at_RI

      !Is it a periodic calculation
      LOGICAL :: do_periodic = .FALSE.

      !Necessary stuff for the on-the fly calculation of the virial
      TYPE(neighbor_list_set_p_type), DIMENSION(:), POINTER :: nl_2c_ovlp => Null(), &
                                                               nl_2c_met => Null(), &
                                                               nl_2c_pot => Null()
      TYPE(neighbor_list_3c_type), POINTER                  :: nl_3c => Null()
      TYPE(dbcsr_type), POINTER                             :: RI_virial_pot => Null(), &
                                                               RI_virial_met => Null()
      TYPE(dbt_type), POINTER                               :: t_3c_virial => Null(), &
                                                               t_3c_virial_split => Null()

   END TYPE im_time_force_type

CONTAINS

! **************************************************************************************************
!> \brief Cleans everything up
!> \param force_data ...
! **************************************************************************************************
   SUBROUTINE im_time_force_release(force_data)

      TYPE(im_time_force_type)                           :: force_data

      INTEGER                                            :: dummy_int, i, i_xyz, j

      CALL dbt_destroy(force_data%t_2c_pot_psqrt)
      CALL dbt_destroy(force_data%t_2c_pot_msqrt)
      CALL dbt_destroy(force_data%t_2c_K)

      CALL dbcsr_release(force_data%inv_ovlp)
      CALL dbcsr_release(force_data%G_PQ)
      CALL dbcsr_deallocate_matrix_set(force_data%P_virt)
      CALL dbcsr_deallocate_matrix_set(force_data%P_occ)
      CALL dbcsr_deallocate_matrix_set(force_data%sum_O_tau)
      CALL dbcsr_deallocate_matrix_set(force_data%sum_YP_tau)

      DO i_xyz = 1, 3
         IF (.NOT. force_data%do_periodic) CALL dbt_destroy(force_data%t_2c_der_pot(i_xyz))
         CALL dbt_destroy(force_data%t_2c_der_ovlp(i_xyz))

         CALL dbt_destroy(force_data%t_3c_der_AO(i_xyz))
         CALL dbt_destroy(force_data%t_3c_der_RI(i_xyz))
      END DO

      CALL dbt_destroy(force_data%t_2c_inv_metric)
      DO i_xyz = 1, 3
         CALL dbt_destroy(force_data%t_2c_der_metric(i_xyz))
      END DO

      DO i = 1, SIZE(force_data%t_3c_der_AO_comp, 1)
         DO j = 1, SIZE(force_data%t_3c_der_AO_comp, 2)
            CALL dealloc_containers(force_data%t_3c_der_AO_comp(i, j), dummy_int)
         END DO
      END DO
      DO i = 1, SIZE(force_data%t_3c_der_RI_comp, 1)
         DO j = 1, SIZE(force_data%t_3c_der_RI_comp, 2)
            CALL dealloc_containers(force_data%t_3c_der_RI_comp(i, j), dummy_int)
         END DO
      END DO
      DEALLOCATE (force_data%t_3c_der_AO_ind, force_data%t_3c_der_RI_ind)

      IF (ASSOCIATED(force_data%nl_2c_ovlp)) THEN
         CALL release_neighbor_list_sets(force_data%nl_2c_ovlp)
      END IF

      IF (ASSOCIATED(force_data%nl_2c_pot)) THEN
         CALL release_neighbor_list_sets(force_data%nl_2c_pot)
      END IF

      IF (ASSOCIATED(force_data%nl_2c_met)) THEN
         CALL release_neighbor_list_sets(force_data%nl_2c_met)
      END IF

      IF (ASSOCIATED(force_data%nl_3c)) THEN
         CALL neighbor_list_3c_destroy(force_data%nl_3c)
         DEALLOCATE (force_data%nl_3c)
      END IF

      IF (ASSOCIATED(force_data%RI_virial_pot)) THEN
         CALL dbcsr_release(force_data%RI_virial_pot)
         DEALLOCATE (force_data%RI_virial_pot)
      END IF

      IF (ASSOCIATED(force_data%RI_virial_met)) THEN
         CALL dbcsr_release(force_data%RI_virial_met)
         DEALLOCATE (force_data%RI_virial_met)
      END IF

      IF (ASSOCIATED(force_data%t_3c_virial)) THEN
         CALL dbt_destroy(force_data%t_3c_virial)
         DEALLOCATE (force_data%t_3c_virial)
      END IF

      IF (ASSOCIATED(force_data%t_3c_virial_split)) THEN
         CALL dbt_destroy(force_data%t_3c_virial_split)
         DEALLOCATE (force_data%t_3c_virial_split)
      END IF

   END SUBROUTINE im_time_force_release

END MODULE rpa_im_time_force_types
