!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Manages a pool of grids (to be used for example as tmp objects),
!>      but can also be used to instantiate grids that are never given back.
!>
!>      Multigrid pools are just an array of pw_pools
!> \note
!>      The pool could also work without pointers (doing = each time),
!>      but I find it *very* ugly.
!>
!>      The pool could be integrated into pw_grid_type, I don't know if
!>      it would be a good or bad idea (but would add a circular dependence
!>      between pw and pw_grid types).
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
MODULE pw_pool_types
   #:include 'pw_types.fypp'
   #:for kind in pw_kinds
      USE cp_linked_list_pw, ONLY: &
cp_sll_${kind[1:]}$_${kind[0]}$_dealloc, cp_sll_${kind[1:]}$_${kind[0]}$_get_first_el, cp_sll_${kind[1:]}$_${kind[0]}$_get_length, &
         cp_sll_${kind[1:]}$_${kind[0]}$_insert_el, cp_sll_${kind[1:]}$_${kind[0]}$_next, cp_sll_${kind[1:]}$_${kind[0]}$_rm_first_el, cp_sll_${kind[1:]}$_${kind[0]}$_type, &
         cp_sll_pw_${kind}$_dealloc, cp_sll_pw_${kind}$_get_first_el, cp_sll_pw_${kind}$_get_length, cp_sll_pw_${kind}$_insert_el, &
         cp_sll_pw_${kind}$_next, cp_sll_pw_${kind}$_rm_first_el, cp_sll_pw_${kind}$_type
   #:endfor
   USE kinds, ONLY: dp
   USE pw_grid_types, ONLY: pw_grid_type
   USE pw_grids, ONLY: pw_grid_compare, &
                       pw_grid_release, &
                       pw_grid_retain
   USE pw_types, ONLY: pw_create, &
#:for kind in pw_kinds
                       pw_${kind}$_type, &
      #:endfor
                       pw_release
#include "../base/base_uses.f90"

   #:include 'pw_types.fypp'

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pw_pool_types'
   INTEGER, PARAMETER :: default_max_cache = 75, max_max_cache = 150

   PUBLIC :: pw_pool_type, pw_pool_p_type
   PUBLIC :: pw_pool_create, pw_pool_retain, pw_pool_release, &
             pw_pool_create_pw, pw_pool_give_back_pw, &
             pw_pool_create_array, pw_pool_give_back_array
   PUBLIC :: pw_pools_copy, pw_pools_dealloc, &
             pw_pools_create_pws, pw_pools_give_back_pws

! **************************************************************************************************
!> \brief Manages a pool of grids (to be used for example as tmp objects),
!>      but can also be used to instantiate grids that are never given back.
!> \param ref_count reference count (see /cp2k/doc/ReferenceCounting.html)
!> \param r1d_pw, r3d_pw, c1d_pw, c3d_pw: liked list with
!>        the cached grids of the corresponding type
!> \note
!>      As of now I would like replace the linked lists by arrays
!>      (no annoying list elements that are allocated would show up when
!>      tracking leaks) [fawzi]
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   TYPE pw_pool_type
      INTEGER :: ref_count = 0, max_cache = 0
      TYPE(pw_grid_type), POINTER :: pw_grid => NULL()
      #:for kind in pw_kinds
         TYPE(cp_sll_pw_${kind}$_type), POINTER :: ${kind}$_pw => NULL()
         TYPE(cp_sll_${kind[1:]}$_${kind[0]}$_type), POINTER :: ${kind}$_array => NULL()
      #:endfor
   END TYPE pw_pool_type

! **************************************************************************************************
!> \brief to create arrays of pools
!> \param pool the pool
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   TYPE pw_pool_p_type
      TYPE(pw_pool_type), POINTER :: pool => NULL()
   END TYPE pw_pool_p_type

   INTERFACE try_pop
      #:for kind in pw_kinds
         MODULE PROCEDURE try_pop_${kind}$
      #:endfor
   END INTERFACE

   INTERFACE pw_pool_create_pw
      #:for kind in pw_kinds
         MODULE PROCEDURE pw_pool_create_pw_${kind}$
      #:endfor
   END INTERFACE

   INTERFACE pw_pool_give_back_pw
      #:for kind in pw_kinds
         MODULE PROCEDURE pw_pool_give_back_pw_${kind}$
      #:endfor
   END INTERFACE

   INTERFACE pw_pools_create_pws
      #:for kind in pw_kinds
         MODULE PROCEDURE pw_pools_create_pws_${kind}$
      #:endfor
   END INTERFACE

   INTERFACE pw_pools_give_back_pws
      #:for kind in pw_kinds
         MODULE PROCEDURE pw_pools_give_back_pws_${kind}$
      #:endfor
   END INTERFACE

   INTERFACE pw_pool_create_array
      #:for kind in pw_kinds
         MODULE PROCEDURE pw_pool_create_array_${kind}$
      #:endfor
   END INTERFACE

   INTERFACE pw_pool_give_back_array
      #:for kind in pw_kinds
         MODULE PROCEDURE pw_pool_give_back_array_${kind}$
      #:endfor
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief creates a pool for pw
!> \param pool the pool to create
!> \param pw_grid the grid that is used to create the pw
!> \param max_cache ...
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE pw_pool_create(pool, pw_grid, max_cache)
      TYPE(pw_pool_type), POINTER                        :: pool
      TYPE(pw_grid_type), POINTER                        :: pw_grid
      INTEGER, OPTIONAL                                  :: max_cache

      ALLOCATE (pool)
      pool%pw_grid => pw_grid
      CALL pw_grid_retain(pw_grid)
      pool%ref_count = 1
      pool%max_cache = default_max_cache
      IF (PRESENT(max_cache)) pool%max_cache = max_cache
      pool%max_cache = MIN(max_max_cache, pool%max_cache)
   END SUBROUTINE pw_pool_create

! **************************************************************************************************
!> \brief retains the pool (see cp2k/doc/ReferenceCounting.html)
!> \param pool the pool to retain
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE pw_pool_retain(pool)
      TYPE(pw_pool_type), POINTER                        :: pool

      CPASSERT(ASSOCIATED(pool))
      CPASSERT(pool%ref_count > 0)

      pool%ref_count = pool%ref_count + 1
   END SUBROUTINE pw_pool_retain

! **************************************************************************************************
!> \brief deallocates all the cached grids
!> \param pool the pool to flush
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE pw_pool_flush_cache(pool)
      TYPE(pw_pool_type), INTENT(INOUT)                  :: pool

      #:for kind, type in zip(pw_kinds, pw_types)
         ${type}$, CONTIGUOUS, POINTER                      :: array_att_${kind}$
         TYPE(cp_sll_${kind[1:]}$_${kind[0]}$_type), POINTER :: array_iterator_${kind}$
         TYPE(cp_sll_pw_${kind}$_type), POINTER             :: iterator_${kind}$
         TYPE(pw_${kind}$_type), POINTER                    :: pw_el_${kind}$
      #:endfor

      #:for kind in pw_kinds
         NULLIFY (array_att_${kind}$, array_iterator_${kind}$, iterator_${kind}$, pw_el_${kind}$)
      #:endfor

      #:for kind in pw_kinds
         iterator_${kind}$ => pool%${kind}$_pw
         DO
            IF (.NOT. cp_sll_pw_${kind}$_next(iterator_${kind}$, el_att=pw_el_${kind}$)) EXIT
            CALL pw_release(pw_el_${kind}$)
            DEALLOCATE (pw_el_${kind}$)
         END DO
         CALL cp_sll_pw_${kind}$_dealloc(pool%${kind}$_pw)

         array_iterator_${kind}$ => pool%${kind}$_array
         DO
            IF (.NOT. cp_sll_${kind[1:]}$_${kind[0]}$_next(array_iterator_${kind}$, el_att=array_att_${kind}$) &
                ) EXIT
            DEALLOCATE (array_att_${kind}$)
         END DO
         CALL cp_sll_${kind[1:]}$_${kind[0]}$_dealloc(pool%${kind}$_array)
      #:endfor

   END SUBROUTINE pw_pool_flush_cache

! **************************************************************************************************
!> \brief releases the given pool (see cp2k/doc/ReferenceCounting.html)
!> \param pool the pool to release
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE pw_pool_release(pool)
      TYPE(pw_pool_type), POINTER                        :: pool

      IF (ASSOCIATED(pool)) THEN
         CPASSERT(pool%ref_count > 0)
         pool%ref_count = pool%ref_count - 1
         IF (pool%ref_count == 0) THEN
            CALL pw_pool_flush_cache(pool)
            CALL pw_grid_release(pool%pw_grid)

            DEALLOCATE (pool)
         END IF
      END IF
      NULLIFY (pool)
   END SUBROUTINE pw_pool_release

   #:for kind, type in zip(pw_kinds, pw_types)

! **************************************************************************************************
!> \brief tries to pop an element from the given list (no error on failure)
!> \param list the list to pop
!> \return ...
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      private function
! **************************************************************************************************
      FUNCTION try_pop_${kind}$ (list) RESULT(res)
         TYPE(cp_sll_pw_${kind}$_type), POINTER             :: list
         TYPE(pw_${kind}$_type), POINTER                    :: res

         IF (ASSOCIATED(list)) THEN
            res => cp_sll_pw_${kind}$_get_first_el(list)
            CALL cp_sll_pw_${kind}$_rm_first_el(list)
         ELSE
            NULLIFY (res)
         END IF
      END FUNCTION try_pop_${kind}$

! **************************************************************************************************
!> \brief returns a pw, allocating it if none is in the pool
!> \param pool the pool from where you get the pw
!> \param pw will contain the new pw
!> \param in_space in which space it is: REALSPACE, RECIPROCALSPACE
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
      SUBROUTINE pw_pool_create_pw_${kind}$ (pool, pw, in_space)
         TYPE(pw_pool_type), INTENT(IN)                     :: pool
         TYPE(pw_${kind}$_type), INTENT(OUT)                :: pw
         INTEGER, INTENT(in)                                :: in_space

         CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_pool_create_pw_${kind}$'

         INTEGER                                            :: handle
         ${type}$, CONTIGUOUS, &
            POINTER                                         :: array_ptr
         TYPE(pw_${kind}$_type), POINTER                    :: el

         CALL timeset(routineN, handle)
         NULLIFY (array_ptr)

         el => try_pop(pool%${kind}$_pw)
         IF (.NOT. ASSOCIATED(el)) THEN
            IF (ASSOCIATED(pool%${kind}$_array)) THEN
               array_ptr => cp_sll_${kind[1:]}$_${kind[0]}$_get_first_el(pool%${kind}$_array)
               CALL cp_sll_${kind[1:]}$_${kind[0]}$_rm_first_el(pool%${kind}$_array)
            END IF
         END IF

         IF (.NOT. ASSOCIATED(el)) THEN
            CALL pw_create(pw, pool%pw_grid, in_space=in_space, &
                           array_ptr=array_ptr)
         ELSE
            pw = el
            DEALLOCATE (el)
            pw%in_space = in_space
         END IF

         CALL timestop(handle)

      END SUBROUTINE pw_pool_create_pw_${kind}$

! **************************************************************************************************
!> \brief returns the pw to the pool
!> \param pool the pool where to reintegrate the pw
!> \param pw the pw to give back
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
      SUBROUTINE pw_pool_give_back_pw_${kind}$ (pool, pw)
         TYPE(pw_pool_type), INTENT(IN)                     :: pool
         TYPE(pw_${kind}$_type), INTENT(INOUT)              :: pw

         CHARACTER(len=*), PARAMETER :: routineN = 'pw_pool_give_back_pw_${kind}$'

         INTEGER                                            :: handle
         TYPE(pw_${kind}$_type), POINTER                    :: el

         CALL timeset(routineN, handle)
         IF (ASSOCIATED(pw%pw_grid)) THEN
            IF (.NOT. pw_grid_compare(pw%pw_grid, pool%pw_grid)) THEN
               CALL pw_release(pw)
            ELSE
               NULLIFY (el)
               ALLOCATE (el)
               el = pw

               IF (ASSOCIATED(pw%array)) THEN
                  IF (cp_sll_pw_${kind}$_get_length(pool%${kind}$_pw) < pool%max_cache) THEN
                     CALL cp_sll_pw_${kind}$_insert_el(pool%${kind}$_pw, el=el)
                  ELSE
                     IF (max_max_cache >= 0) &
                        CPWARN("hit max_cache")
                     CALL pw_release(el)
                     DEALLOCATE (el)
                  END IF
               ELSE
                  CALL pw_release(el)
                  DEALLOCATE (el)
               END IF
            END IF
         ELSE
            IF (ASSOCIATED(pw%array)) DEALLOCATE (pw%array)
         END IF
         pw = pw_${kind}$_type()
         CALL timestop(handle)
      END SUBROUTINE pw_pool_give_back_pw_${kind}$

! **************************************************************************************************
!> \brief creates a multigrid structure
!> \param pools the multigrid pool (i.e. an array of pw_pool)
!> \param pws the multigrid of coefficent you want to initialize
!> \param in_space ...
!> \par History
!>      07.2004 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
      SUBROUTINE pw_pools_create_pws_${kind}$ (pools, pws, in_space)
         TYPE(pw_pool_p_type), DIMENSION(:), INTENT(IN)     :: pools
         TYPE(pw_${kind}$_type), ALLOCATABLE, DIMENSION(:), &
            INTENT(OUT)                                     :: pws
         INTEGER, INTENT(in)                                :: in_space

         INTEGER                                            :: i

         ALLOCATE (pws(SIZE(pools)))
         DO i = 1, SIZE(pools)
            CALL pw_pool_create_pw_${kind}$ (pools(i)%pool, pws(i), &
                                             in_space=in_space)
         END DO
      END SUBROUTINE pw_pools_create_pws_${kind}$

! **************************************************************************************************
!> \brief returns the pw part of the coefficients into the pools
!> \param pools the pools that will cache the pws %pw
!> \param pws the coefficients to give back
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
      SUBROUTINE pw_pools_give_back_pws_${kind}$ (pools, pws)
         TYPE(pw_pool_p_type), DIMENSION(:), INTENT(IN)     :: pools
         TYPE(pw_${kind}$_type), ALLOCATABLE, DIMENSION(:), &
            INTENT(INOUT)                                   :: pws

         INTEGER                                            :: i

         CPASSERT(SIZE(pws) == SIZE(pools))
         DO i = 1, SIZE(pools)
            CALL pw_pool_give_back_pw_${kind}$ (pools(i)%pool, pws(i))
         END DO
         DEALLOCATE (pws)
      END SUBROUTINE pw_pools_give_back_pws_${kind}$

! **************************************************************************************************
!> \brief returns an array of coefficients as the one used by pw
!>      allocating it if none is present in the pool
!> \param pw_pool the pool that caches the cr3d
!> \param array the pointer that will contain the array
!> \par History
!>      11.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
      SUBROUTINE pw_pool_create_array_${kind}$ (pw_pool, array)
         TYPE(pw_pool_type), INTENT(IN)                     :: pw_pool
         ${type}$, POINTER                                  :: array

         TYPE(pw_${kind}$_type), POINTER                    :: pw

         NULLIFY (pw)

         IF (ASSOCIATED(pw_pool%${kind}$_array)) THEN
            array => cp_sll_${kind[1:]}$_${kind[0]}$_get_first_el(pw_pool%${kind}$_array)
            CALL cp_sll_${kind[1:]}$_${kind[0]}$_rm_first_el(pw_pool%${kind}$_array)
         ELSE
            pw => try_pop(pw_pool%${kind}$_pw)
            IF (ASSOCIATED(pw)) THEN
               array => pw%array
               NULLIFY (pw%array)
               CALL pw_release(pw)
               DEALLOCATE (pw)
            END IF
         END IF
         IF (.NOT. ASSOCIATED(array)) THEN
            #:if kind.endswith('1d')
               ALLOCATE (array(pw_pool%pw_grid%ngpts_cut_local))
            #:else
               ALLOCATE (array(pw_pool%pw_grid%bounds_local(1, 1):pw_pool%pw_grid%bounds_local(2, 1), &
                               pw_pool%pw_grid%bounds_local(1, 2):pw_pool%pw_grid%bounds_local(2, 2), &
                               pw_pool%pw_grid%bounds_local(1, 3):pw_pool%pw_grid%bounds_local(2, 3)))
            #:endif
         END IF
      END SUBROUTINE pw_pool_create_array_${kind}$

! **************************************************************************************************
!> \brief returns an array of coefficients as the one used by pw,
!>      allocating it if none is present in the pool
!> \param pw_pool the pool that caches the cr3d
!> \param array the pointer that will contain the array
!> \param accept_non_compatible if true deallocates the non compatible
!>        arrays passed in, if false (the default) stops with an error
!> \par History
!>      11.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
      SUBROUTINE pw_pool_give_back_array_${kind}$ (pw_pool, array, accept_non_compatible)
         TYPE(pw_pool_type), INTENT(IN)                     :: pw_pool
         ${type}$, CONTIGUOUS, POINTER                      :: array
         LOGICAL, INTENT(in), OPTIONAL                      :: accept_non_compatible

         LOGICAL                                            :: compatible, my_accept_non_compatible

         my_accept_non_compatible = .FALSE.
         IF (PRESENT(accept_non_compatible)) my_accept_non_compatible = accept_non_compatible

         IF (ASSOCIATED(array)) THEN
            #:if kind.endswith('3d')
               compatible = ALL(MERGE(pw_pool%pw_grid%bounds_local(1, :) == LBOUND(array) .AND. &
                                      pw_pool%pw_grid%bounds_local(2, :) == UBOUND(array), &
                                      pw_pool%pw_grid%bounds_local(2, :) < pw_pool%pw_grid%bounds_local(1, :), &
                                      UBOUND(array) >= LBOUND(array)))
            #:else
               compatible = UBOUND(array, 1) >= LBOUND(array, 1) .AND. &
                            (1 == LBOUND(array, 1) .AND. pw_pool%pw_grid%ngpts_cut_local == UBOUND(array, 1)) .OR. &
                            (UBOUND(array, 1) < LBOUND(array, 1) .AND. pw_pool%pw_grid%ngpts_cut_local < 0)
            #:endif
            CPASSERT(compatible .OR. my_accept_non_compatible)
            IF (compatible) THEN
               IF (cp_sll_${kind[1:]}$_${kind[0]}$_get_length(pw_pool%${kind}$_array) < pw_pool%max_cache) THEN
                  CALL cp_sll_${kind[1:]}$_${kind[0]}$_insert_el(pw_pool%${kind}$_array, el=array)
               ELSE
                  IF (max_max_cache >= 0) &
                     CPWARN("hit max_cache")
                  DEALLOCATE (array)
               END IF
            ELSE
               DEALLOCATE (array)
            END IF
         ELSE
            CPASSERT(my_accept_non_compatible)
         END IF
         NULLIFY (array)
      END SUBROUTINE pw_pool_give_back_array_${kind}$

   #:endfor

! **************************************************************************************************
!> \brief copies a multigrid pool, the underlying pools are shared
!> \param source_pools the pools to copy
!> \param target_pools will hold the copy of the pools
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE pw_pools_copy(source_pools, target_pools)
      TYPE(pw_pool_p_type), DIMENSION(:), INTENT(IN)     :: source_pools
      TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: target_pools

      INTEGER                                            :: i

      ALLOCATE (target_pools(SIZE(source_pools)))
      DO i = 1, SIZE(source_pools)
         target_pools(i)%pool => source_pools(i)%pool
         CALL pw_pool_retain(source_pools(i)%pool)
      END DO
   END SUBROUTINE pw_pools_copy

! **************************************************************************************************
!> \brief deallocates the given pools (releasing each of the underlying
!>      pools)
!> \param pools the pols to deallocate
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE pw_pools_dealloc(pools)
      TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pools

      INTEGER                                            :: i

      IF (ASSOCIATED(pools)) THEN
         DO i = 1, SIZE(pools)
            CALL pw_pool_release(pools(i)%pool)
         END DO
         DEALLOCATE (pools)
      END IF
      NULLIFY (pools)
   END SUBROUTINE pw_pools_dealloc

END MODULE pw_pool_types
