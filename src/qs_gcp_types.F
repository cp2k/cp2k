!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
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
      INTEGER                                :: za
      REAL(KIND=dp)                          :: asto
      REAL(KIND=dp)                          :: rcsto
      INTEGER                                :: nq
      REAL(KIND=dp)                          :: nbvirt
      REAL(KIND=dp)                          :: eamiss
      REAL(KIND=dp), DIMENSION(6)            :: al
      REAL(KIND=dp), DIMENSION(6)            :: cl
   END TYPE qs_gcp_kind_type
! **************************************************************************************************
   TYPE qs_gcp_type
      LOGICAL                                :: do_gcp = .FALSE.
      LOGICAL                                :: verbose !extended output
      CHARACTER(LEN=default_string_length)   :: parameter_file_name
      ! parameter input
      CHARACTER(LEN=default_string_length), &
         DIMENSION(:), POINTER               :: kind_type => NULL()
      REAL(KIND=dp), DIMENSION(:), POINTER   :: ea => NULL()
      !global parameters
      REAL(KIND=dp)                          :: alpha, beta, sigma, eta
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

      CHARACTER(len=*), PARAMETER                        :: routineN = 'qs_gcp_release'

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

