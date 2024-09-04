!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Definition of gCP types for DFT calculations
!> \author JGH (20.10.2018)
! **************************************************************************************************
MODULE qs_gcp_types

   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type,&
                                              release_neighbor_list_sets
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_gcp_types'

! **************************************************************************************************
   TYPE qs_gcp_kind_type
      INTEGER                                :: za = -1
      REAL(KIND=dp)                          :: asto = -1.0_dp
      REAL(KIND=dp)                          :: rcsto = -1.0_dp
      INTEGER                                :: nq = -1
      REAL(KIND=dp)                          :: nbvirt = -1.0_dp
      REAL(KIND=dp)                          :: eamiss = -1.0_dp
      REAL(KIND=dp), DIMENSION(6)            :: al = -1.0_dp
      REAL(KIND=dp), DIMENSION(6)            :: cl = -1.0_dp
   END TYPE qs_gcp_kind_type
! **************************************************************************************************
   TYPE qs_gcp_type
      LOGICAL                                :: do_gcp = .FALSE.
      LOGICAL                                :: verbose = .FALSE. !extended output
      CHARACTER(LEN=default_string_length)   :: parameter_file_name = ""
      ! parameter input
      CHARACTER(LEN=default_string_length), &
         DIMENSION(:), POINTER               :: kind_type => NULL()
      REAL(KIND=dp), DIMENSION(:), POINTER   :: ea => NULL()
      !global parameters
      REAL(KIND=dp)                          :: alpha = -1.0_dp, beta = -1.0_dp, sigma = -1.0_dp, eta = -1.0_dp
      !neighborlist
      TYPE(neighbor_list_set_p_type), &
         DIMENSION(:), POINTER               :: sab_gcp => Null() ! neighborlists for pair interactions
      !kind information
      TYPE(qs_gcp_kind_type), DIMENSION(:), &
         POINTER                             :: gcp_kind => NULL() ! atomic kind parameters
   END TYPE qs_gcp_type
! **************************************************************************************************

   PUBLIC :: qs_gcp_type
   PUBLIC :: qs_gcp_release

! **************************************************************************************************
CONTAINS
! **************************************************************************************************
!> \brief ...
!> \param gcp_env ...
! **************************************************************************************************
   SUBROUTINE qs_gcp_release(gcp_env)

      TYPE(qs_gcp_type), POINTER                         :: gcp_env

      IF (ASSOCIATED(gcp_env)) THEN

         CALL release_neighbor_list_sets(gcp_env%sab_gcp)
         IF (ASSOCIATED(gcp_env%kind_type)) THEN
            DEALLOCATE (gcp_env%kind_type)
         END IF
         IF (ASSOCIATED(gcp_env%ea)) THEN
            DEALLOCATE (gcp_env%ea)
         END IF

         IF (ASSOCIATED(gcp_env%gcp_kind)) THEN
            DEALLOCATE (gcp_env%gcp_kind)
         END IF

         DEALLOCATE (gcp_env)

      END IF

   END SUBROUTINE qs_gcp_release

! **************************************************************************************************

END MODULE qs_gcp_types

