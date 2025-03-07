!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Storage to keep precomputed surface Green's functions
! **************************************************************************************************
MODULE negf_green_cache
   USE cp_cfm_types,                    ONLY: cp_cfm_release,&
                                              cp_cfm_type
   USE kinds,                           ONLY: dp
   USE util,                            ONLY: sort
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'negf_green_cache'
   LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .TRUE.

   PUBLIC :: green_functions_cache_type

   PUBLIC :: green_functions_cache_expand, &
             green_functions_cache_reorder, &
             green_functions_cache_release

! **************************************************************************************************
!> \brief Storage to keep surface Green's functions.
!> \author Sergey Chulkov
! **************************************************************************************************
   TYPE green_functions_cache_type
      !> retarded surface Green's functions [ncontacts, nnodes]
      TYPE(cp_cfm_type), ALLOCATABLE, DIMENSION(:, :)   :: g_surf_contacts
      !> list of points over the normalised interval [-1 .. 1].
      !> Coordinates of actual point where Green's functions were evaluated
      !> can be obtained by using an appropriate rescale_nodes_*() subroutine
      !> from the module 'negf_integr_utils'.
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: tnodes
   END TYPE green_functions_cache_type

CONTAINS
! **************************************************************************************************
!> \brief Reallocate storage so it can handle extra 'nnodes_extra' items for each contact.
!> \param cache        storage to expand
!> \param ncontacts    number of contacts
!> \param nnodes_extra number of items to add
!> \author Sergey Chulkov
! **************************************************************************************************
   SUBROUTINE green_functions_cache_expand(cache, ncontacts, nnodes_extra)
      TYPE(green_functions_cache_type), INTENT(inout)    :: cache
      INTEGER, INTENT(in)                                :: ncontacts, nnodes_extra

      INTEGER                                            :: nentries_exist
      LOGICAL                                            :: is_alloc
      TYPE(cp_cfm_type), ALLOCATABLE, DIMENSION(:, :)    :: g_surf_contacts

      is_alloc = ALLOCATED(cache%g_surf_contacts)

      IF (is_alloc) THEN
         CPASSERT(SIZE(cache%g_surf_contacts, 1) == ncontacts)
         nentries_exist = SIZE(cache%g_surf_contacts, 2)

      ELSE
         nentries_exist = 0
      END IF

      ALLOCATE (g_surf_contacts(ncontacts, nentries_exist + nnodes_extra))

      IF (is_alloc) THEN
         g_surf_contacts(1:ncontacts, 1:nentries_exist) = cache%g_surf_contacts(1:ncontacts, 1:nentries_exist)
         DEALLOCATE (cache%g_surf_contacts)
      END IF

      CALL MOVE_ALLOC(g_surf_contacts, cache%g_surf_contacts)
   END SUBROUTINE green_functions_cache_expand

! **************************************************************************************************
!> \brief Sort cached items in ascending order.
!> \param cache        storage to reorder
!> \param tnodes       coordinate of items in storage
!> \author Sergey Chulkov
! **************************************************************************************************
   SUBROUTINE green_functions_cache_reorder(cache, tnodes)
      TYPE(green_functions_cache_type), INTENT(inout)    :: cache
      REAL(kind=dp), DIMENSION(:), INTENT(in)            :: tnodes

      INTEGER                                            :: ind_new, ind_old, ncontacts, nnodes
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: indices
      TYPE(cp_cfm_type), ALLOCATABLE, DIMENSION(:, :)    :: g_surf_contacts

      nnodes = SIZE(tnodes)

      CPASSERT(ALLOCATED(cache%g_surf_contacts))
      CPASSERT(SIZE(cache%g_surf_contacts, 2) == nnodes)

      ncontacts = SIZE(cache%g_surf_contacts, 1)

      IF (ALLOCATED(cache%tnodes)) DEALLOCATE (cache%tnodes)

      ALLOCATE (g_surf_contacts(ncontacts, nnodes))
      ALLOCATE (cache%tnodes(nnodes))
      ALLOCATE (indices(nnodes))

      cache%tnodes(:) = tnodes(:)
      CALL sort(cache%tnodes, nnodes, indices)

      DO ind_new = 1, nnodes
         ind_old = indices(ind_new)
         g_surf_contacts(1:ncontacts, ind_new) = cache%g_surf_contacts(1:ncontacts, ind_old)
      END DO

      CALL MOVE_ALLOC(g_surf_contacts, cache%g_surf_contacts)
   END SUBROUTINE green_functions_cache_reorder

! **************************************************************************************************
!> \brief Release storage.
!> \param cache        storage to release
!> \author Sergey Chulkov
! **************************************************************************************************
   SUBROUTINE green_functions_cache_release(cache)
      TYPE(green_functions_cache_type), INTENT(inout)    :: cache

      INTEGER                                            :: icontact, ipoint, ncontacts

      IF (ALLOCATED(cache%tnodes)) DEALLOCATE (cache%tnodes)

      IF (ALLOCATED(cache%g_surf_contacts)) THEN
         ncontacts = SIZE(cache%g_surf_contacts, 1)
         DO ipoint = SIZE(cache%g_surf_contacts, 2), 1, -1
            DO icontact = ncontacts, 1, -1
               CALL cp_cfm_release(cache%g_surf_contacts(icontact, ipoint))
            END DO
         END DO

         DEALLOCATE (cache%g_surf_contacts)
      END IF
   END SUBROUTINE green_functions_cache_release
END MODULE negf_green_cache

