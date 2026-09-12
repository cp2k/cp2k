!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Input and persistent data for automatic RI basis optimization.
!> \par History
!>      09.2026 created [Jan Wilhelm]
! **************************************************************************************************
MODULE gw_auto_ri_types
   USE base_hooks,                      ONLY: timeset,&
                                              timestop
   USE cp_fm_types,                     ONLY: cp_fm_release,&
                                              cp_fm_type
   USE kinds,                           ONLY: dp

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_auto_ri_types'

   PUBLIC :: auto_ri_release, auto_ri_type

   TYPE auto_ri_type

      ! Input parameters.
      REAL(KIND=dp)                   :: occ_energy_window = -1.0_dp
      REAL(KIND=dp)                   :: neighbor_radius = -1.0_dp
      LOGICAL                         :: enabled = .FALSE.
      REAL(KIND=dp)                   :: ri_ao_ratio = -1.0_dp

      ! Two-center matrices M_pq and V_pq of optimized RI basis set
      TYPE(cp_fm_type)                :: V_pq = cp_fm_type()
      TYPE(cp_fm_type)                :: M_pq_inv = cp_fm_type()
      LOGICAL                         :: ready = .FALSE.
      INTEGER, ALLOCATABLE            :: sizes_opt_RI(:)

      ! U_Pp_AB contains the RI-basis contraction coefficients:
      ! for A=B: φ_p^A(r)  = Σ_P U_Pp^AA φ_P^A(r)    (φ_P^A(r): large ref. RI basis)
      ! for A≠B: φ_p^AB(r) = Σ_P U_Pp^AB φ_P^AB^⊥(r) (see forthcoming paper)
      ! RI basis contains contractions of Gaussians on neighboring atoms!!
      REAL(KIND=dp), ALLOCATABLE      :: U_Pp_AB(:)
      INTEGER                         :: AB_block_count = 0
      INTEGER, ALLOCATABLE            :: AB_atom_A(:), AB_atom_B(:)
      INTEGER, ALLOCATABLE            :: AB_first_p_A(:), AB_first_p_B(:)
      INTEGER, ALLOCATABLE            :: AB_size_opt_RI_to_A(:)
      INTEGER, ALLOCATABLE            :: AB_size_ref_RI(:)
      INTEGER, ALLOCATABLE            :: AB_size_opt_RI(:)
      INTEGER, ALLOCATABLE            :: U_Pp_AB_offset(:)

   END TYPE auto_ri_type

CONTAINS

! **************************************************************************************************
!> \brief Releases all data owned by an optimized AUTO_RI basis.
!> \param auto_ri AUTO_RI input and optimized-basis data
! **************************************************************************************************
   SUBROUTINE auto_ri_release(auto_ri)
      TYPE(auto_ri_type), INTENT(INOUT)                  :: auto_ri

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'auto_ri_release'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      IF (ALLOCATED(auto_ri%sizes_opt_RI)) DEALLOCATE (auto_ri%sizes_opt_RI)
      IF (ALLOCATED(auto_ri%AB_atom_A)) DEALLOCATE (auto_ri%AB_atom_A)
      IF (ALLOCATED(auto_ri%AB_atom_B)) DEALLOCATE (auto_ri%AB_atom_B)
      IF (ALLOCATED(auto_ri%AB_first_p_A)) DEALLOCATE (auto_ri%AB_first_p_A)
      IF (ALLOCATED(auto_ri%AB_first_p_B)) DEALLOCATE (auto_ri%AB_first_p_B)
      IF (ALLOCATED(auto_ri%AB_size_opt_RI_to_A)) DEALLOCATE (auto_ri%AB_size_opt_RI_to_A)
      IF (ALLOCATED(auto_ri%AB_size_ref_RI)) DEALLOCATE (auto_ri%AB_size_ref_RI)
      IF (ALLOCATED(auto_ri%AB_size_opt_RI)) DEALLOCATE (auto_ri%AB_size_opt_RI)
      IF (ALLOCATED(auto_ri%U_Pp_AB_offset)) DEALLOCATE (auto_ri%U_Pp_AB_offset)
      IF (ALLOCATED(auto_ri%U_Pp_AB)) DEALLOCATE (auto_ri%U_Pp_AB)
      CALL cp_fm_release(auto_ri%V_pq)
      CALL cp_fm_release(auto_ri%M_pq_inv)
      auto_ri%AB_block_count = 0
      auto_ri%ready = .FALSE.

      CALL timestop(handle)

   END SUBROUTINE auto_ri_release

END MODULE gw_auto_ri_types
