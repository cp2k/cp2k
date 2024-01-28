!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \note
!>      If parallel mode is distributed certain combination of
!>      "in_use" and "in_space" can not be used.
!>      For performance reasons it would be better to have the loops
!>      over g-vectros in the gather/scatter routines in new subprograms
!>      with the actual arrays (also the addressing) in the parameter list
!> \par History
!>      JGH (29-Dec-2000) : Changes for parallel use
!>      JGH (13-Mar-2001) : added timing calls
!>      JGH (26-Feb-2003) : OpenMP enabled
!>      JGH (17-Nov-2007) : Removed mass arrays
!>      JGH (01-Dec-2007) : Removed and renamed routines
!>      03.2008 [tlaino] : Splitting pw_types into pw_types and pw_methods
!> \author apsi
! **************************************************************************************************
MODULE pw_types

   USE cp_log_handling, ONLY: cp_to_string
   USE kinds, ONLY: dp
   USE pw_grid_types, ONLY: pw_grid_type
   USE pw_grids, ONLY: pw_grid_release, &
                       pw_grid_retain
#include "../base/base_uses.f90"

   IMPLICIT NONE

   #:include 'pw_types.fypp'

   PRIVATE
   #:for kind in pw_kinds
      PUBLIC :: pw_${kind}$_type, pw_${kind}$_p_type
   #:endfor

   ! Flags for the structure member 'in_space'
   ! NOSPACE only for internal usage, we enforce the choice of a proper space
   INTEGER, PARAMETER, PUBLIC :: NOSPACE = 371
   INTEGER, PARAMETER, PUBLIC :: REALSPACE = 372, RECIPROCALSPACE = 373

   #:for kind, type in zip(pw_kinds, pw_types)
! **************************************************************************************************
      TYPE pw_${kind}$_type
         ${type}$, CONTIGUOUS, POINTER :: array => NULL()
         TYPE(pw_grid_type), POINTER :: pw_grid => NULL()
         INTEGER :: in_space = NOSPACE
      CONTAINS
         PROCEDURE, PUBLIC, NON_OVERRIDABLE :: create => pw_create_${kind}$
         PROCEDURE, PUBLIC, NON_OVERRIDABLE :: release => pw_release_${kind}$
      END TYPE pw_${kind}$_type

! **************************************************************************************************
      TYPE pw_${kind}$_p_type
         TYPE(pw_${kind}$_type), POINTER :: pw => NULL()
      END TYPE pw_${kind}$_p_type
   #:endfor

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pw_types'
   LOGICAL, PARAMETER, PRIVATE :: debug_this_module = .FALSE.

CONTAINS

   #:for kind, type in zip(pw_kinds, pw_types)
! **************************************************************************************************
!> \brief releases the given pw
!> \param pw the pw to release
!> \par History
!>      04.2003 created [fawzi]
!> \author fawzi
!> \note
!>      see doc/ReferenceCounting.html
! **************************************************************************************************
      SUBROUTINE pw_release_${kind}$ (pw)
         CLASS(pw_${kind}$_type), INTENT(INOUT)                       :: pw

         IF (ASSOCIATED(pw%array)) DEALLOCATE (pw%array)
         CALL pw_grid_release(pw%pw_grid)
      END SUBROUTINE pw_release_${kind}$

! **************************************************************************************************
!> \brief allocates and initializes pw_r3d_type
!> \param pw the type that will bw allocated and initialized
!> \param pw_grid ...
!> \param use_data which kind of data will be used
!> \param in_space in which space the pw is (real or reciprocal)
!> \param array_ptr pointer with to data
!> \par History
!>      11.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
      SUBROUTINE pw_create_${kind}$ (pw, pw_grid, in_space, array_ptr)
         CLASS(pw_${kind}$_type), INTENT(INOUT)                         :: pw
         TYPE(pw_grid_type), INTENT(IN), POINTER            :: pw_grid
         INTEGER, INTENT(in)                                :: in_space
         ${type}$, CONTIGUOUS, INTENT(IN), OPTIONAL, POINTER                   :: array_ptr

         CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_create_${kind}$'

         INTEGER                                            :: handle

         CALL timeset(routineN, handle)

         IF (in_space /= REALSPACE .AND. in_space /= RECIPROCALSPACE) &
            CPABORT("Choose a proper space in which the grid lives!")

         ! Ensure a clean grid to prevent memory leaks
         CALL pw%release()

         pw%pw_grid => pw_grid
         CALL pw_grid_retain(pw%pw_grid)

         #:if kind[1]=="1"
         IF (PRESENT(array_ptr)) THEN
         IF (ASSOCIATED(array_ptr)) THEN
            CPASSERT(SIZE(array_ptr) == pw%pw_grid%ngpts_cut_local)
            pw%array(1:pw%pw_grid%ngpts_cut_local) => array_ptr
         END IF
         END IF
         IF (.NOT. ASSOCIATED(pw%array)) ALLOCATE (pw%array(pw%pw_grid%ngpts_cut_local))
         #:elif kind[1]=="3"
         ASSOCIATE (bounds => pw%pw_grid%bounds_local)
            IF (PRESENT(array_ptr)) THEN
               IF (ASSOCIATED(array_ptr)) THEN
                  IF (ALL(bounds(1, :) <= bounds(2, :))) THEN
                     CPASSERT(ALL(LBOUND(array_ptr) == bounds(1, :)))
                     CPASSERT(ALL(UBOUND(array_ptr) == bounds(2, :)))
                  END IF
                  pw%array => array_ptr
               END IF
            END IF
            IF (.NOT. ASSOCIATED(pw%array)) THEN
               ALLOCATE (pw%array( &
                         bounds(1, 1):bounds(2, 1), &
                         bounds(1, 2):bounds(2, 2), &
                         bounds(1, 3):bounds(2, 3)))
            END IF
         END ASSOCIATE
         #:endif
         pw%in_space = in_space
         CALL timestop(handle)
      END SUBROUTINE pw_create_${kind}$
   #:endfor

END MODULE pw_types
