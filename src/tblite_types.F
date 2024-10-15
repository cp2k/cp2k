!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief interface to tblite [Johann Pototschnig]
!> \History creation 09.2024
! **************************************************************************************************
MODULE tblite_types

#if defined(__TBLITE)
   USE mctc_io, ONLY: structure_type, new
   USE tblite_xtb_calculator, ONLY: xtb_calculator
   USE tblite_container, ONLY: container_cache
   USE tblite_wavefunction_type, ONLY: wavefunction_type
   USE tblite_scf_potential, ONLY: potential_type
#endif

   USE kinds, ONLY: dp

#include "./base/base_uses.f90"
   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tblite_types'

   PUBLIC :: tblite_type, deallocate_tblite_type, allocate_tblite_type

   TYPE tblite_type

      INTEGER, ALLOCATABLE                               :: el_num(:)

      REAL(KIND=dp), DIMENSION(3, 3)                     :: sigma
      REAL(KIND=dp), ALLOCATABLE                         :: e_hal(:)
      REAL(KIND=dp), ALLOCATABLE                         :: e_rep(:)
      REAL(KIND=dp), ALLOCATABLE                         :: e_disp(:)
      REAL(KIND=dp), ALLOCATABLE                         :: e_scd(:)
      REAL(KIND=dp), ALLOCATABLE                         :: e_es(:)
      REAL(KIND=dp), ALLOCATABLE                         :: selfenergy(:)
      REAL(KIND=dp), ALLOCATABLE                         :: cn(:)

      REAL(KIND=dp), ALLOCATABLE                         :: grad(:, :)
      REAL(KIND=dp), ALLOCATABLE                         :: dsedcn(:)
      REAL(KIND=dp), ALLOCATABLE                         :: dcndr(:, :, :)
      REAL(KIND=dp), ALLOCATABLE                         :: dcndL(:, :, :)

#if defined(__TBLITE)
      TYPE(structure_type)                               :: mol
      TYPE(xtb_calculator)                               :: calc
      TYPE(container_cache)                              :: dcache, cache
      TYPE(wavefunction_type)                            :: wfn
      TYPE(potential_type)                               :: pot
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

         DEALLOCATE (tb_tblite)
      END IF

   END SUBROUTINE deallocate_tblite_type

END MODULE tblite_types
