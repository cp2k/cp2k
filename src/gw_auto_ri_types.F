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
      LOGICAL                         :: enabled = .FALSE.
      REAL(KIND=dp)                   :: ri_ao_ratio = -1.0_dp

      ! Two-center matrices M_pq and V_pq of optimized RI basis set
      TYPE(cp_fm_type)                :: V_pq = cp_fm_type()
      TYPE(cp_fm_type)                :: M_pq_inv = cp_fm_type()
      LOGICAL                         :: ready = .FALSE.
      INTEGER, ALLOCATABLE            :: sizes_opt_RI(:)
      INTEGER, ALLOCATABLE            :: sizes_ref_RI(:)
      REAL(KIND=dp), ALLOCATABLE      :: U_Pp_atom(:)
      INTEGER, ALLOCATABLE            :: U_Pp_atom_offset(:)

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
      IF (ALLOCATED(auto_ri%sizes_ref_RI)) DEALLOCATE (auto_ri%sizes_ref_RI)
      IF (ALLOCATED(auto_ri%U_Pp_atom)) DEALLOCATE (auto_ri%U_Pp_atom)
      IF (ALLOCATED(auto_ri%U_Pp_atom_offset)) DEALLOCATE (auto_ri%U_Pp_atom_offset)
      CALL cp_fm_release(auto_ri%V_pq)
      CALL cp_fm_release(auto_ri%M_pq_inv)
      auto_ri%ready = .FALSE.

      CALL timestop(handle)

   END SUBROUTINE auto_ri_release

END MODULE gw_auto_ri_types
