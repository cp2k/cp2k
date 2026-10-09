!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
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
   USE mctc_io, ONLY: structure_type
   USE tblite_xtb_calculator, ONLY: xtb_calculator
   USE tblite_container_cache, ONLY: container_cache
   USE tblite_scf_mixer, ONLY: mixer_type
   USE tblite_wavefunction_type, ONLY: wavefunction_type
   USE tblite_scf_potential, ONLY: potential_type
   USE tblite_param, ONLY: param_record
#endif

   USE cp_dbcsr_api, ONLY: dbcsr_p_type
   USE cp_dbcsr_operations, ONLY: dbcsr_deallocate_matrix_set
   USE input_constants, ONLY: tblite_mixer_damping_default, &
                              tblite_mixer_iterations_default, &
                              tblite_mixer_max_weight_default, &
                              tblite_mixer_min_weight_default, &
                              tblite_mixer_omega0_default, &
                              tblite_mixer_weight_factor_default
   USE kinds, ONLY: dp

#include "../base/base_uses.f90"
   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tblite_types'

   PUBLIC :: tblite_type, deallocate_tblite_type, allocate_tblite_type

   TYPE tblite_type

      LOGICAL                                            :: use_virial = .FALSE.
      INTEGER, ALLOCATABLE                               :: el_num(:)

      REAL(KIND=dp)                                      :: accuracy = 1.0_dp
      REAL(KIND=dp), DIMENSION(3, 3)                     :: sigma = -1.0_dp
      ! tblite contributions are currently replicated on every MPI rank.
      ! Force/virial assembly accounts for CP2K's subsequent MPI reductions.
      REAL(KIND=dp), ALLOCATABLE                         :: e_hal(:)
      REAL(KIND=dp), ALLOCATABLE                         :: e_rep(:)
      REAL(KIND=dp), ALLOCATABLE                         :: e_disp(:)
      REAL(KIND=dp), ALLOCATABLE                         :: e_scd(:)
      REAL(KIND=dp), ALLOCATABLE                         :: e_es(:)
      REAL(KIND=dp), ALLOCATABLE                         :: e_int(:)
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
      ! Spin-resolved density snapshot used by the UKS/LSD force path.
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: rho_ao_kp_ref => Null()
      INTEGER                                            :: mixer_memory = tblite_mixer_iterations_default
      REAL(KIND=dp)                                      :: mixer_damping = tblite_mixer_damping_default
      REAL(KIND=dp)                                      :: mixer_max_weight = tblite_mixer_max_weight_default
      REAL(KIND=dp)                                      :: mixer_min_weight = tblite_mixer_min_weight_default
      REAL(KIND=dp)                                      :: mixer_omega0 = tblite_mixer_omega0_default
      REAL(KIND=dp)                                      :: mixer_weight_factor = tblite_mixer_weight_factor_default

#if defined(__TBLITE)
      TYPE(structure_type)                               :: mol = structure_type()
      TYPE(xtb_calculator)                               :: calc = xtb_calculator()
      ! Geometry-dependent caches: Coulomb, dispersion and additional interactions.
      ! Refreshed by tb_init_ham, reused for SCC potentials, energies and gradients.
      TYPE(container_cache)                              :: cache = container_cache()
      TYPE(container_cache)                              :: dcache = container_cache()
      TYPE(container_cache)                              :: icache = container_cache()
      TYPE(wavefunction_type)                            :: wfn = wavefunction_type()
      TYPE(potential_type)                               :: pot = potential_type()
      CLASS(mixer_type), ALLOCATABLE                     :: mixer
      TYPE(param_record), ALLOCATABLE                    :: param
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
         ! Release pointer-owned DBCSR objects; allocatable components are released
         ! automatically when the parent object is deallocated.
         IF (ASSOCIATED(tb_tblite%dipbra)) THEN
            CALL dbcsr_deallocate_matrix_set(tb_tblite%dipbra)
         END IF
         IF (ASSOCIATED(tb_tblite%dipket)) THEN
            CALL dbcsr_deallocate_matrix_set(tb_tblite%dipket)
         END IF
         IF (ASSOCIATED(tb_tblite%quadbra)) THEN
            CALL dbcsr_deallocate_matrix_set(tb_tblite%quadbra)
         END IF
         IF (ASSOCIATED(tb_tblite%quadket)) THEN
            CALL dbcsr_deallocate_matrix_set(tb_tblite%quadket)
         END IF
         IF (ASSOCIATED(tb_tblite%rho_ao_kp_ref)) THEN
            CALL dbcsr_deallocate_matrix_set(tb_tblite%rho_ao_kp_ref)
         END IF
         DEALLOCATE (tb_tblite)
         NULLIFY (tb_tblite)
      END IF

   END SUBROUTINE deallocate_tblite_type

END MODULE tblite_types
