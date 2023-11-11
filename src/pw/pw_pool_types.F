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
   USE cp_linked_list_pw,               ONLY: &
        cp_sll_3d_r_dealloc, cp_sll_3d_r_get_first_el, cp_sll_3d_r_get_length, &
        cp_sll_3d_r_insert_el, cp_sll_3d_r_next, cp_sll_3d_r_rm_first_el, cp_sll_3d_r_type, &
        cp_sll_pw_dealloc, cp_sll_pw_get_first_el, cp_sll_pw_get_length, cp_sll_pw_insert_el, &
        cp_sll_pw_next, cp_sll_pw_rm_first_el, cp_sll_pw_type
   USE kinds,                           ONLY: dp
   USE pw_grid_types,                   ONLY: pw_grid_type
   USE pw_grids,                        ONLY: pw_grid_compare,&
                                              pw_grid_release,&
                                              pw_grid_retain
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              COMPLEXDATA3D,&
                                              REALDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_type
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pw_pool_types'
   INTEGER, PARAMETER :: default_max_cache = 75, max_max_cache = 150

   PUBLIC :: pw_pool_type, pw_pool_p_type
   PUBLIC :: pw_pool_create, pw_pool_release
   PUBLIC :: pw_pools_copy, pw_pools_dealloc, &
             pw_pools_create_pws, pw_pools_give_back_pws

! **************************************************************************************************
!> \brief Manages a pool of grids (to be used for example as tmp objects),
!>      but can also be used to instantiate grids that are never given back.
!> \param ref_count reference count (see /cp2k/doc/ReferenceCounting.html)
!> \param real 1d_pw, real3d_pw, complex1d_pw, complex3d_pw: liked list with
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
      TYPE(cp_sll_pw_type), POINTER :: real1d_pw => NULL(), real3d_pw => NULL(), &
                                       complex1d_pw => NULL(), complex3d_pw => NULL()
      TYPE(cp_sll_3d_r_type), POINTER :: real3d_array => NULL()
   CONTAINS
      PROCEDURE, PUBLIC, NON_OVERRIDABLE :: retain => pw_pool_retain
      PROCEDURE, PUBLIC, NON_OVERRIDABLE :: create_pw => pw_pool_create_pw
      PROCEDURE, PUBLIC, NON_OVERRIDABLE :: give_back_pw => pw_pool_give_back_pw
      PROCEDURE, PUBLIC, NON_OVERRIDABLE :: create_cr3d => pw_pool_create_cr3d
      PROCEDURE, PUBLIC, NON_OVERRIDABLE :: give_back_cr3d => pw_pool_give_back_cr3d
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
      NULLIFY (pool%real1d_pw, pool%real3d_pw, &
               pool%complex1d_pw, pool%complex3d_pw, pool%real3d_array)
   END SUBROUTINE pw_pool_create

! **************************************************************************************************
!> \brief retains the pool (see cp2k/doc/ReferenceCounting.html)
!> \param pool the pool to retain
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE pw_pool_retain(pool)
      CLASS(pw_pool_type), INTENT(INOUT)                  :: pool

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

      REAL(kind=dp), CONTIGUOUS, DIMENSION(:, :, :), &
         POINTER                                         :: array_att
      TYPE(cp_sll_3d_r_type), POINTER                    :: array_iterator
      TYPE(cp_sll_pw_type), POINTER                      :: iterator
      TYPE(pw_type), POINTER                             :: pw_el

      NULLIFY (iterator, array_iterator, pw_el, array_att)
      iterator => pool%real1d_pw
      DO
         IF (.NOT. cp_sll_pw_next(iterator, el_att=pw_el)) EXIT
         CALL pw_el%release()
         DEALLOCATE (pw_el)
      END DO
      CALL cp_sll_pw_dealloc(pool%real1d_pw)

      iterator => pool%real3d_pw
      DO
         IF (.NOT. cp_sll_pw_next(iterator, el_att=pw_el)) EXIT
         CALL pw_el%release()
         DEALLOCATE (pw_el)
      END DO
      CALL cp_sll_pw_dealloc(pool%real3d_pw)

      iterator => pool%complex1d_pw
      DO
         IF (.NOT. cp_sll_pw_next(iterator, el_att=pw_el)) EXIT
         CALL pw_el%release()
         DEALLOCATE (pw_el)
      END DO
      CALL cp_sll_pw_dealloc(pool%complex1d_pw)

      iterator => pool%complex3d_pw
      DO
         IF (.NOT. cp_sll_pw_next(iterator, el_att=pw_el)) EXIT
         CALL pw_el%release()
         DEALLOCATE (pw_el)
      END DO
      CALL cp_sll_pw_dealloc(pool%complex3d_pw)

      array_iterator => pool%real3d_array
      DO
         IF (.NOT. cp_sll_3d_r_next(array_iterator, el_att=array_att) &
             ) EXIT
         DEALLOCATE (array_att)
      END DO
      CALL cp_sll_3d_r_dealloc(pool%real3d_array)

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
   FUNCTION try_pop(list) RESULT(res)
      TYPE(cp_sll_pw_type), POINTER                      :: list
      TYPE(pw_type), POINTER                             :: res

      IF (ASSOCIATED(list)) THEN
         res => cp_sll_pw_get_first_el(list)
         CALL cp_sll_pw_rm_first_el(list)
      ELSE
         NULLIFY (res)
      END IF
   END FUNCTION try_pop

! **************************************************************************************************
!> \brief returns a pw, allocating it if none is in the pool
!> \param pool the pool from where you get the pw
!> \param pw will contain the new pw
!> \param use_data which data it uses: REALDATA1D, COMPLEXDATA1D,
!>        REALDATA3D, COMPLEXDATA3D
!> \param in_space in which space it is: REALSPACE, RECIPROCALSPACE
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE pw_pool_create_pw(pool, pw, use_data, in_space)
      CLASS(pw_pool_type), INTENT(IN)                     :: pool
      TYPE(pw_type), INTENT(OUT)                         :: pw
      INTEGER, INTENT(in)                                :: use_data, in_space

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_pool_create_pw'

      INTEGER                                            :: handle
      REAL(kind=dp), CONTIGUOUS, DIMENSION(:, :, :), &
         POINTER                                         :: cr3d_ptr
      TYPE(pw_type), POINTER                             :: el

      CALL timeset(routineN, handle)
      NULLIFY (cr3d_ptr)

      SELECT CASE (use_data)
      CASE (REALDATA1D)
         el => try_pop(pool%real1d_pw)
      CASE (REALDATA3D)
         el => try_pop(pool%real3d_pw)
         IF (.NOT. ASSOCIATED(el)) THEN
            IF (ASSOCIATED(pool%real3d_array)) THEN
               cr3d_ptr => cp_sll_3d_r_get_first_el(pool%real3d_array)
               CALL cp_sll_3d_r_rm_first_el(pool%real3d_array)
            END IF
         END IF
      CASE (COMPLEXDATA1D)
         el => try_pop(pool%complex1d_pw)
      CASE (COMPLEXDATA3D)
         el => try_pop(pool%complex3d_pw)
      CASE default
! unknown use_data
         CPABORT("")
      END SELECT

      IF (.NOT. ASSOCIATED(el)) THEN
         CALL pw%create(pool%pw_grid, use_data=use_data, &
                        in_space=in_space, cr3d_ptr=cr3d_ptr)
      ELSE
         pw = el
         DEALLOCATE (el)
         IF (in_space /= REALSPACE .AND. in_space /= RECIPROCALSPACE) &
            CPABORT("A PW grid either lives in real space or in reciprocal space!")
         pw%in_space = in_space
      END IF

      CALL timestop(handle)

   END SUBROUTINE pw_pool_create_pw

! **************************************************************************************************
!> \brief returns the pw to the pool
!> \param pool the pool where to reintegrate the pw
!> \param pw the pw to give back
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE pw_pool_give_back_pw(pool, pw)
      CLASS(pw_pool_type), INTENT(IN)                     :: pool
      TYPE(pw_type), INTENT(INOUT)                       :: pw

      CHARACTER(len=*), PARAMETER :: routineN = 'pw_pool_give_back_pw'

      INTEGER                                            :: handle
      TYPE(pw_type), POINTER                             :: el

      CALL timeset(routineN, handle)
      IF (ASSOCIATED(pw%pw_grid)) THEN
         IF (.NOT. pw_grid_compare(pw%pw_grid, pool%pw_grid)) THEN
            CALL pw%release()
         ELSE
            NULLIFY (el)
            ALLOCATE (el)
            el = pw

            SELECT CASE (el%in_use)
            CASE (REALDATA1D)
               IF (cp_sll_pw_get_length(pool%real1d_pw) < pool%max_cache) THEN
                  CALL cp_sll_pw_insert_el(pool%real1d_pw, el=el)
               ELSE
                  IF (max_max_cache >= 0) &
                     CPWARN("hit max_cache")
                  CALL el%release()
                  DEALLOCATE (el)
               END IF
            CASE (REALDATA3D)
               IF (ASSOCIATED(el%cr3d)) THEN
                  IF (cp_sll_pw_get_length(pool%real3d_pw) < pool%max_cache) THEN
                     CALL cp_sll_pw_insert_el(pool%real3d_pw, el=el)
                  ELSE
                     IF (max_max_cache >= 0) &
                        CPWARN("hit max_cache")
                     CALL el%release()
                     DEALLOCATE (el)
                  END IF
               ELSE
                  CALL el%release()
                  DEALLOCATE (el)
               END IF
            CASE (COMPLEXDATA1D)
               IF (cp_sll_pw_get_length(pool%complex1d_pw) < pool%max_cache) THEN
                  CALL cp_sll_pw_insert_el(pool%complex1d_pw, el=el)
               ELSE
                  IF (max_max_cache >= 0) &
                     CPWARN("hit max_cache")
                  CALL el%release()
                  DEALLOCATE (el)
               END IF
            CASE (COMPLEXDATA3D)
               IF (cp_sll_pw_get_length(pool%complex3d_pw) < pool%max_cache) THEN
                  CALL cp_sll_pw_insert_el(pool%complex3d_pw, el=el)
               ELSE
                  IF (max_max_cache >= 0) &
                     CPWARN("hit max_cache")
                  CALL el%release()
                  DEALLOCATE (el)
               END IF
            CASE default
               ! unknown in_use
               CPABORT("")
            END SELECT
            !FM so that if someone tries to use a pw that is in the pool
            !FM (s)he gets problems
         END IF
      ELSE
         IF (ASSOCIATED(pw%cc)) DEALLOCATE (pw%cc)
         IF (ASSOCIATED(pw%cc3d)) DEALLOCATE (pw%cc3d)
         IF (ASSOCIATED(pw%cr)) DEALLOCATE (pw%cr)
         IF (ASSOCIATED(pw%cr3d)) DEALLOCATE (pw%cr3d)
      END IF
      NULLIFY (pw%cc, pw%cc3d, pw%cr, pw%cr3d, pw%pw_grid)
      CALL timestop(handle)
   END SUBROUTINE pw_pool_give_back_pw

! **************************************************************************************************
!> \brief returns a 3d real array of coefficients as the one used by pw with
!>      REALDATA3D, allocating it if none is present in the pool
!> \param pw_pool the pool that caches the cr3d
!> \param cr3d the pointer that will contain the array
!> \par History
!>      11.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE pw_pool_create_cr3d(pw_pool, cr3d)
      CLASS(pw_pool_type), INTENT(IN)                     :: pw_pool
      REAL(kind=dp), DIMENSION(:, :, :), POINTER         :: cr3d

      TYPE(pw_type), POINTER                             :: pw

      NULLIFY (pw)
      IF (ASSOCIATED(pw_pool%real3d_array)) THEN
         cr3d => cp_sll_3d_r_get_first_el(pw_pool%real3d_array)
         CALL cp_sll_3d_r_rm_first_el(pw_pool%real3d_array)
      ELSE
         pw => try_pop(pw_pool%real3d_pw)
         IF (ASSOCIATED(pw)) THEN
            cr3d => pw%cr3d
            NULLIFY (pw%cr3d)
            CALL pw%release()
            DEALLOCATE (pw)
         END IF
      END IF
      IF (.NOT. ASSOCIATED(cr3d)) THEN
         ALLOCATE (cr3d(pw_pool%pw_grid%bounds_local(1, 1):pw_pool%pw_grid%bounds_local(2, 1), &
                        pw_pool%pw_grid%bounds_local(1, 2):pw_pool%pw_grid%bounds_local(2, 2), &
                        pw_pool%pw_grid%bounds_local(1, 3):pw_pool%pw_grid%bounds_local(2, 3)))
      END IF
   END SUBROUTINE pw_pool_create_cr3d

! **************************************************************************************************
!> \brief returns a 3d real array of coefficients as the one used by pw with
!>      REALDATA3D, allocating it if none is present in the pool
!> \param pw_pool the pool that caches the cr3d
!> \param cr3d the pointer that will contain the array
!> \par History
!>      11.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE pw_pool_give_back_cr3d(pw_pool, cr3d)
      CLASS(pw_pool_type), INTENT(IN)                     :: pw_pool
      REAL(kind=dp), CONTIGUOUS, DIMENSION(:, :, :), &
         POINTER                                         :: cr3d

      LOGICAL                                            :: compatible

      IF (ASSOCIATED(cr3d)) THEN
         compatible = ALL(MERGE(pw_pool%pw_grid%bounds_local(1, :) == LBOUND(cr3d) .AND. &
                                pw_pool%pw_grid%bounds_local(2, :) == UBOUND(cr3d), &
                                pw_pool%pw_grid%bounds_local(2, :) < pw_pool%pw_grid%bounds_local(1, :), &
                                UBOUND(cr3d) >= LBOUND(cr3d)))
         IF (compatible) THEN
            IF (cp_sll_3d_r_get_length(pw_pool%real3d_array) < pw_pool%max_cache) THEN
               CALL cp_sll_3d_r_insert_el(pw_pool%real3d_array, el=cr3d)
            ELSE
               IF (max_max_cache >= 0) &
                  CPWARN("hit max_cache")
               DEALLOCATE (cr3d)
            END IF
         ELSE
            DEALLOCATE (cr3d)
         END IF
      END IF
      NULLIFY (cr3d)
   END SUBROUTINE pw_pool_give_back_cr3d

! **************************************************************************************************
!> \brief creates a multigrid structure
!> \param pools the multigrid pool (i.e. an array of pw_pool)
!> \param pws the multigrid of coefficent you want to initialize
!> \param use_data which data it uses: REALDATA1D, COMPLEXDATA1D,
!>        REALDATA3D, COMPLEXDATA3D
!> \param in_space ...
!> \par History
!>      07.2004 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE pw_pools_create_pws(pools, pws, use_data, in_space)
      TYPE(pw_pool_p_type), DIMENSION(:), INTENT(IN)     :: pools
      TYPE(pw_type), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: pws
      INTEGER, INTENT(in)                                :: use_data
      INTEGER, INTENT(in), OPTIONAL                      :: in_space

      INTEGER                                            :: i

      ALLOCATE (pws(SIZE(pools)))
      DO i = 1, SIZE(pools)
         CALL pools(i)%pool%create_pw(pws(i), use_data, &
                                      in_space=in_space)
      END DO
   END SUBROUTINE pw_pools_create_pws

! **************************************************************************************************
!> \brief returns the pw part of the coefficients into the pools
!> \param pools the pools that will cache the pws %pw
!> \param pws the coefficients to give back
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE pw_pools_give_back_pws(pools, pws)
      TYPE(pw_pool_p_type), DIMENSION(:), INTENT(IN)     :: pools
      TYPE(pw_type), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: pws

      INTEGER                                            :: i

      CPASSERT(SIZE(pws) == SIZE(pools))
      DO i = 1, SIZE(pools)
         CALL pools(i)%pool%give_back_pw(pws(i))
      END DO
      DEALLOCATE (pws)
   END SUBROUTINE pw_pools_give_back_pws

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
         CALL source_pools(i)%pool%retain()
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
