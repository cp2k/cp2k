!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief types for tblite
!> \author JVP
!> \history creation 09.2024
! **************************************************************************************************

MODULE tblite_types

#if defined(__TBLITE)
   USE mctc_io, ONLY: structure_type, new
   USE tblite_xtb_calculator, ONLY: xtb_calculator
   USE tblite_container, ONLY: container_cache
   USE tblite_wavefunction_type, ONLY: wavefunction_type
   USE tblite_scf_potential, ONLY: potential_type
#endif

   USE cp_dbcsr_api, ONLY: dbcsr_p_type
   USE cp_dbcsr_operations, ONLY: dbcsr_deallocate_matrix_set
   USE kinds, ONLY: dp

#include "./base/base_uses.f90"
   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tblite_types'

   PUBLIC :: tblite_type, deallocate_tblite_type, allocate_tblite_type

   TYPE tblite_type

      LOGICAL                                            :: use_virial = .FALSE.
      INTEGER, ALLOCATABLE                               :: el_num(:)

      REAL(KIND=dp), DIMENSION(3, 3)                     :: sigma = -1.0_dp
      REAL(KIND=dp), ALLOCATABLE                         :: e_hal(:)
      REAL(KIND=dp), ALLOCATABLE                         :: e_rep(:)
      REAL(KIND=dp), ALLOCATABLE                         :: e_disp(:)
      REAL(KIND=dp), ALLOCATABLE                         :: e_scd(:)
      REAL(KIND=dp), ALLOCATABLE                         :: e_es(:)
      REAL(KIND=dp), ALLOCATABLE                         :: selfenergy(:)
      REAL(KIND=dp), ALLOCATABLE                         :: cn(:)

      REAL(KIND=dp), ALLOCATABLE                         :: grad(:, :)
      !derivatives w. r. t. coordination number
      REAL(KIND=dp), ALLOCATABLE                         :: dsedcn(:)
      REAL(KIND=dp), ALLOCATABLE                         :: dcndr(:, :, :)
      REAL(KIND=dp), ALLOCATABLE                         :: dcndL(:, :, :)

      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: dipbra => Null(), &
                                                            dipket => Null()
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: quadbra => Null(), &
                                                            quadket => Null()

#if defined(__TBLITE)
      TYPE(structure_type)                               :: mol = structure_type()
      TYPE(xtb_calculator)                               :: calc = xtb_calculator()
      TYPE(container_cache)                              :: cache = container_cache()
      TYPE(container_cache)                              :: dcache = container_cache()
      TYPE(wavefunction_type)                            :: wfn = wavefunction_type()
      TYPE(potential_type)                               :: pot = potential_type()
#endif

   END TYPE tblite_type

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param tb_tblite ...
! **************************************************************************************************
   SUBROUTINE allocate_tblite_type(tb_tblite)

      TYPE(tblite_type), POINTER                         :: tb_tblite

      CALL deallocate_tblite_type(tb_tblite)

      ALLOCATE (tb_tblite)

   END SUBROUTINE allocate_tblite_type

! **************************************************************************************************
!> \brief ...
!> \param tb_tblite ...
! **************************************************************************************************
   SUBROUTINE deallocate_tblite_type(tb_tblite)

      TYPE(tblite_type), POINTER                         :: tb_tblite

      IF (ASSOCIATED(tb_tblite)) THEN
         IF (ALLOCATED(tb_tblite%el_num)) DEALLOCATE (tb_tblite%el_num)

         IF (ALLOCATED(tb_tblite%e_hal)) DEALLOCATE (tb_tblite%e_hal)
         IF (ALLOCATED(tb_tblite%e_rep)) DEALLOCATE (tb_tblite%e_rep)
         IF (ALLOCATED(tb_tblite%e_disp)) DEALLOCATE (tb_tblite%e_disp)
         IF (ALLOCATED(tb_tblite%e_scd)) DEALLOCATE (tb_tblite%e_scd)
         IF (ALLOCATED(tb_tblite%e_es)) DEALLOCATE (tb_tblite%e_es)
         IF (ALLOCATED(tb_tblite%selfenergy)) DEALLOCATE (tb_tblite%selfenergy)
         IF (ALLOCATED(tb_tblite%cn)) DEALLOCATE (tb_tblite%cn)

         IF (ALLOCATED(tb_tblite%grad)) DEALLOCATE (tb_tblite%grad)
         IF (ALLOCATED(tb_tblite%dsedcn)) DEALLOCATE (tb_tblite%dsedcn)
         IF (ALLOCATED(tb_tblite%dcndr)) DEALLOCATE (tb_tblite%dcndr)
         IF (ALLOCATED(tb_tblite%dcndL)) DEALLOCATE (tb_tblite%dcndL)

         IF (ASSOCIATED(tb_tblite%dipbra)) &
            CALL dbcsr_deallocate_matrix_set(tb_tblite%dipbra)
         IF (ASSOCIATED(tb_tblite%dipket)) &
            CALL dbcsr_deallocate_matrix_set(tb_tblite%dipket)
         IF (ASSOCIATED(tb_tblite%quadbra)) &
            CALL dbcsr_deallocate_matrix_set(tb_tblite%quadbra)
         IF (ASSOCIATED(tb_tblite%quadket)) &
            CALL dbcsr_deallocate_matrix_set(tb_tblite%quadket)

         DEALLOCATE (tb_tblite)
      END IF

   END SUBROUTINE deallocate_tblite_type

END MODULE tblite_types
