!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
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
! *****************************************************************************
MODULE pw_pool_types
  USE cp_linked_list_3d_r,             ONLY: cp_sll_3d_r_dealloc,&
                                             cp_sll_3d_r_get_first_el,&
                                             cp_sll_3d_r_get_length,&
                                             cp_sll_3d_r_insert_el,&
                                             cp_sll_3d_r_next,&
                                             cp_sll_3d_r_rm_first_el,&
                                             cp_sll_3d_r_type
  USE cp_linked_list_pw,               ONLY: cp_sll_pw_dealloc,&
                                             cp_sll_pw_get_first_el,&
                                             cp_sll_pw_get_length,&
                                             cp_sll_pw_insert_el,&
                                             cp_sll_pw_insert_ordered,&
                                             cp_sll_pw_next,&
                                             cp_sll_pw_rm_first_el,&
                                             cp_sll_pw_type
  USE kinds,                           ONLY: dp
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_grids,                        ONLY: pw_grid_compare,&
                                             pw_grid_release,&
                                             pw_grid_retain
  USE pw_methods,                      ONLY: pw_write
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             COMPLEXDATA3D,&
                                             REALDATA1D,&
                                             REALDATA3D,&
                                             pw_create,&
                                             pw_p_type,&
                                             pw_release,&
                                             pw_type
  USE timings,                         ONLY: print_stack,&
                                             timeset,&
                                             timestop
#include "../common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.FALSE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pw_pool_types'
  INTEGER, SAVE, PRIVATE :: last_pw_pool_id_nr=0
  INTEGER, PARAMETER :: default_max_cache=75, max_max_cache=150

  PUBLIC :: pw_pool_type, pw_pool_p_type
  PUBLIC :: pw_pool_create, pw_pool_retain, pw_pool_release,&
            pw_pool_create_pw, pw_pool_give_back_pw,&
            pw_pool_create_cr3d, pw_pool_give_back_cr3d
  PUBLIC :: pw_pools_copy, pw_pools_dealloc, pw_pools_flush_cache,&
            pw_pools_create_pws, pw_pools_give_back_pws

! *****************************************************************************
!> \brief Manages a pool of grids (to be used for example as tmp objects),
!>      but can also be used to intantiate grids that are never given back.
!> \param ref_count reference count (see /cp2k/doc/ReferenceCounting.html)
!> \param id_nr number that identifies each pool
!> \param real 1d_pw, real3d_pw, complex1d_pw, complex3d_pw: liked list with
!>        the cached grids of the corresponding type
!> \note
!>      As of now I would like replace the linked lists by arrays
!>      (no annoying list elements that are allocated would show up when
!>      tracking leaks) [fawzi]
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  TYPE pw_pool_type
     INTEGER :: ref_count, id_nr, max_cache
     TYPE(pw_grid_type), POINTER :: pw_grid
     TYPE(cp_sll_pw_type), POINTER :: real1d_pw, real3d_pw,&
          complex1d_pw, complex3d_pw
     TYPE(cp_sll_3d_r_type), POINTER :: real3d_array
  END TYPE pw_pool_type

! *****************************************************************************
!> \brief to create arrays of pools
!> \param pool the pool
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  TYPE pw_pool_p_type
     TYPE(pw_pool_type), POINTER :: pool
  END TYPE pw_pool_p_type

CONTAINS

! *****************************************************************************
!> \brief creates a pool for pw
!> \param pool the pool to create
!> \param pw_grid the grid that is used to create the pw
!> \param max_cache ...
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE pw_pool_create(pool, pw_grid, max_cache)
    TYPE(pw_pool_type), POINTER              :: pool
    TYPE(pw_grid_type), POINTER              :: pw_grid
    INTEGER, OPTIONAL                        :: max_cache

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_pool_create', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger

    failure=.FALSE.
    logger => cp_get_default_logger()

    ALLOCATE(pool)
    pool%pw_grid => pw_grid
    CALL pw_grid_retain(pw_grid)
    last_pw_pool_id_nr=last_pw_pool_id_nr+1
    pool%id_nr=last_pw_pool_id_nr
    pool%ref_count=1
    pool%max_cache=default_max_cache
    IF (PRESENT(max_cache)) pool%max_cache=max_cache
    pool%max_cache=MIN(max_max_cache,pool%max_cache)
    IF (debug_this_module) THEN
       WRITE (unit=cp_logger_get_default_unit_nr(logger,local=.TRUE.),&
            fmt="(' *** pw_pool ',i4,' has been created')") pool%id_nr
       CALL print_stack(cp_logger_get_default_unit_nr(logger,local=.TRUE.))
    END IF
    NULLIFY(pool%real1d_pw, pool%real3d_pw, &
         pool%complex1d_pw, pool%complex3d_pw, pool%real3d_array)
  END SUBROUTINE pw_pool_create

! *****************************************************************************
!> \brief retains the pool (see cp2k/doc/ReferenceCounting.html)
!> \param pool the pool to retain
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE pw_pool_retain(pool)
    TYPE(pw_pool_type), POINTER              :: pool

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_pool_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger

    failure=.FALSE.
    logger => cp_get_default_logger()

    CPASSERT(ASSOCIATED(pool))
    CPASSERT(pool%ref_count>0)

    pool%ref_count=pool%ref_count+1
    IF (debug_this_module) THEN
       WRITE (unit=cp_logger_get_default_unit_nr(logger,local=.TRUE.),&
            fmt="(' *** pw_pool ',i4,' has been retained, ref_count=',i4)")&
            pool%id_nr, pool%ref_count
       CALL print_stack(cp_logger_get_default_unit_nr(logger,local=.TRUE.))
    END IF
  END SUBROUTINE pw_pool_retain

! *****************************************************************************
!> \brief deallocates all the cached grids
!> \param pool the pool to flush
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE pw_pool_flush_cache(pool)
    TYPE(pw_pool_type), POINTER              :: pool

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_pool_flush_cache', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    REAL(kind=dp), DIMENSION(:, :, :), &
      POINTER                                :: array_att
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_sll_3d_r_type), POINTER          :: array_iterator
    TYPE(cp_sll_pw_type), POINTER            :: iterator
    TYPE(pw_type), POINTER                   :: pw_el

    failure=.FALSE.

    CPASSERT(ASSOCIATED(pool))
    CPASSERT(pool%ref_count>0)
    NULLIFY(iterator,array_iterator,pw_el,array_att)
    logger => cp_get_default_logger()
    IF (debug_this_module) THEN
       WRITE (unit=cp_logger_get_default_unit_nr(logger,local=.TRUE.),&
            fmt="(' *** pw_pool ',i4,' is flushing the cache')") pool%id_nr
       CALL print_stack(cp_logger_get_default_unit_nr(logger,local=.TRUE.))
    END IF

    iterator => pool%real1d_pw
    DO
       IF (.NOT.cp_sll_pw_next(iterator,el_att=pw_el)) EXIT
       CPASSERT(pw_el%ref_count==0)
       pw_el%ref_count=1
       CALL pw_release(pw_el)
    END DO
    CALL cp_sll_pw_dealloc(pool%real1d_pw)

    iterator => pool%real3d_pw
    DO
       IF (.NOT.cp_sll_pw_next(iterator,el_att=pw_el)) EXIT
       CPASSERT(pw_el%ref_count==0)
       pw_el%ref_count=1
       CALL pw_release(pw_el)
    END DO
    CALL cp_sll_pw_dealloc(pool%real3d_pw)

    iterator => pool%complex1d_pw
    DO
       IF (.NOT.cp_sll_pw_next(iterator,el_att=pw_el)) EXIT
       CPASSERT(pw_el%ref_count==0)
       pw_el%ref_count=1
       CALL pw_release(pw_el)
    END DO
    CALL cp_sll_pw_dealloc(pool%complex1d_pw)

    iterator => pool%complex3d_pw
    DO
       IF (.NOT.cp_sll_pw_next(iterator,el_att=pw_el)) EXIT
       CPASSERT(pw_el%ref_count==0)
       pw_el%ref_count=1
       CALL pw_release(pw_el)
    END DO
    CALL cp_sll_pw_dealloc(pool%complex3d_pw)

    array_iterator => pool%real3d_array
    DO
       IF (.NOT.cp_sll_3d_r_next(array_iterator,el_att=array_att)&
            ) EXIT
       DEALLOCATE(array_att)
    END DO
    CALL cp_sll_3d_r_dealloc(pool%real3d_array)

  END SUBROUTINE pw_pool_flush_cache

! *****************************************************************************
!> \brief releases the given pool (see cp2k/doc/ReferenceCounting.html)
!> \param pool the pool to release
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE pw_pool_release(pool)
    TYPE(pw_pool_type), POINTER              :: pool

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_pool_release', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger

    failure=.FALSE.
    logger => cp_get_default_logger()

    IF (ASSOCIATED(pool)) THEN
       CPASSERT(pool%ref_count>0)
       pool%ref_count=pool%ref_count-1
       IF (debug_this_module) THEN
          WRITE (unit=cp_logger_get_default_unit_nr(logger,local=.TRUE.),&
               fmt="(' *** pw_pool ',i4,' released ref_count=',i4)") &
               pool%id_nr, pool%ref_count
          CALL print_stack(cp_logger_get_default_unit_nr(logger,local=.TRUE.))
       END IF
       IF (pool%ref_count==0) THEN
          pool%ref_count=1
          CALL pw_pool_flush_cache(pool)
          pool%ref_count=0
          CPASSERT(ASSOCIATED(pool%pw_grid))
          CALL pw_grid_release(pool%pw_grid)

          DEALLOCATE(pool)
       END IF
    END IF
    NULLIFY(pool)
  END SUBROUTINE pw_pool_release

! *****************************************************************************
!> \brief tries to pop an element from the given list (no error on failure)
!> \param list the list to pop
!> \retval res ...
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      private function
! *****************************************************************************
  FUNCTION try_pop(list) RESULT(res)
    TYPE(cp_sll_pw_type), POINTER            :: list
    TYPE(pw_type), POINTER                   :: res

    IF (ASSOCIATED(list)) THEN
       res => cp_sll_pw_get_first_el(list)
       CALL cp_sll_pw_rm_first_el(list)
    ELSE
       NULLIFY(res)
    END IF
  END FUNCTION try_pop

! *****************************************************************************
!> \brief returns a pw, allocating it if none is in the pool
!> \param pool the pool from where you get the pw
!> \param pw will contain the new pw
!> \param use_data which data it uses: REALDATA1D, COMPLEXDATA1D,
!>        REALDATA3D, COMPLEXDATA3D
!> \param in_space in which space it is: REALSPACE, RECIPROCALSPACE
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE pw_pool_create_pw(pool, pw, use_data, in_space)
    TYPE(pw_pool_type), POINTER              :: pool
    TYPE(pw_type), POINTER                   :: pw
    INTEGER, INTENT(in)                      :: use_data
    INTEGER, INTENT(in), OPTIONAL            :: in_space

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_pool_create_pw', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: failure
    REAL(kind=dp), DIMENSION(:, :, :), &
      POINTER                                :: cr3d_ptr
    TYPE(cp_logger_type), POINTER            :: logger

    failure=.FALSE.

    CALL timeset(routineN,handle)
    NULLIFY(pw)
    NULLIFY(cr3d_ptr)
    logger => cp_get_default_logger()
    CPASSERT(ASSOCIATED(pool))
    CPASSERT(pool%ref_count>0)

    SELECT CASE(use_data)
    CASE (REALDATA1D)
       pw => try_pop(pool%real1d_pw)
    CASE (REALDATA3D)
       pw => try_pop(pool%real3d_pw)
       IF (.NOT.ASSOCIATED(pw)) THEN
          IF (ASSOCIATED(pool%real3d_array)) THEN
             cr3d_ptr => cp_sll_3d_r_get_first_el(pool%real3d_array)
             CALL cp_sll_3d_r_rm_first_el(pool%real3d_array)
          END IF
       END IF
    CASE (COMPLEXDATA1D)
       pw => try_pop(pool%complex1d_pw)
    CASE (COMPLEXDATA3D)
       pw => try_pop(pool%complex3d_pw)
    CASE default
! unknown use_data
       CPASSERT(.FALSE.)
    END SELECT

    IF (.NOT.ASSOCIATED(pw)) THEN
       CALL pw_create(pw, pool%pw_grid, use_data=use_data, &
            cr3d_ptr=cr3d_ptr)
       IF (debug_this_module) THEN
          IF (ASSOCIATED(cr3d_ptr)) THEN
             WRITE(unit=cp_logger_get_default_unit_nr(logger,local=.TRUE.),&
                  fmt="(' *** pw_pool ',i4,' created pw ',i4,'reusing 3dr array')")&
                  pool%id_nr,pw%id_nr
          ELSE
             WRITE(unit=cp_logger_get_default_unit_nr(logger,local=.TRUE.),&
                  fmt="(' *** pw_pool ',i4,' created pw ',i4)")&
                  pool%id_nr,pw%id_nr
          END IF
          CALL pw_write(&
               unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.),&
               pw=pw)
          CALL print_stack(cp_logger_get_default_unit_nr(logger,local=.TRUE.))
       END IF
    ELSE
       CPASSERT(pw%ref_count==0)
       pw%ref_count=1
       IF (debug_this_module) THEN
          WRITE(unit=cp_logger_get_default_unit_nr(logger,local=.TRUE.),&
               fmt="(' *** pw_pool ',i4,' created pw reusing old ',i4)")&
               pool%id_nr,pw%id_nr
          CALL pw_write(&
               unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.),&
               pw=pw)
          CALL print_stack(cp_logger_get_default_unit_nr(logger,local=.TRUE.))
       END IF
    END IF

    pw%in_space=0
    IF (PRESENT(in_space)) pw%in_space=in_space

    CALL timestop(handle)

  END SUBROUTINE pw_pool_create_pw

! *****************************************************************************
!> \brief returns the pw to the pool
!> \param pool the pool where to reintegrate the pw
!> \param pw the pw to give back
!> \param accept_non_compatible if non compatible pw should be accepted
!>        (they will be destroied). Defaults to false (and thus stops with
!>        an error)
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE pw_pool_give_back_pw(pool, pw, accept_non_compatible)
    TYPE(pw_pool_type), POINTER              :: pool
    TYPE(pw_type), POINTER                   :: pw
    LOGICAL, INTENT(in), OPTIONAL            :: accept_non_compatible

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_pool_give_back_pw', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle
    LOGICAL                                  :: did_insert, failure, &
                                                my_accept_non_compatible
    TYPE(cp_logger_type), POINTER            :: logger

    failure=.FALSE.
    my_accept_non_compatible=.FALSE.
    logger => cp_get_default_logger()
    IF (PRESENT(accept_non_compatible)) my_accept_non_compatible=accept_non_compatible

    CALL timeset(routineN,handle)
    CPASSERT(ASSOCIATED(pool))
    CPASSERT(pool%ref_count>0)
    IF (.NOT.ASSOCIATED(pw)) THEN
       CPASSERT(my_accept_non_compatible)
       failure=.TRUE.
    END IF
    IF (.NOT. failure) THEN
       CPASSERT(pw%ref_count==1)
       IF (.NOT.pw_grid_compare(pw%pw_grid,pool%pw_grid)) THEN
          IF (debug_this_module) THEN
             WRITE (unit=cp_logger_get_default_unit_nr(logger,local=.TRUE.),&
                  fmt="(' *** pw_pool ',i4,' giving back incompatible pw ',i4)")&
                  pool%id_nr, pw%id_nr
             CALL pw_write(unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.),&
                  pw=pw)
             CALL print_stack(cp_logger_get_default_unit_nr(logger,local=.TRUE.))
          END IF

          CALL cp_assert(my_accept_non_compatible,&
               cp_failure_level, cp_assertion_failed, routineP,&
               "pool cannot reuse pw of another grid "//&
CPSourceFileRef,&
               failure=failure)
          CALL pw_release(pw)
          failure=.TRUE.
       END IF
    END IF

    IF (.NOT. failure) THEN

       IF (debug_this_module) THEN
          WRITE (unit=cp_logger_get_default_unit_nr(logger,local=.TRUE.),&
               fmt="(' *** pw_pool ',i4,' giving back pw ',i4)")&
               pool%id_nr, pw%id_nr
          CALL pw_write(unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.),&
               pw=pw)
          CALL print_stack(cp_logger_get_default_unit_nr(logger,local=.TRUE.))
       END IF

       SELECT CASE(pw%in_use)
       CASE(REALDATA1D)
          IF (cp_sll_pw_get_length(pool%real1d_pw)<pool%max_cache) THEN
             IF (cp_debug) THEN
                CALL cp_sll_pw_insert_ordered(pool%real1d_pw, el=pw,&
                     insert_equals=.FALSE., did_insert=did_insert)
                CPASSERT(did_insert)
             ELSE
                CALL cp_sll_pw_insert_el(pool%real1d_pw, el=pw)
             END IF
          ELSE
             CALL cp_assert(max_max_cache<0,cp_warning_level,cp_assertion_failed,&
                  routineP,"hit max_cache"//&
CPSourceFileRef)
             CALL pw_release(pw)
          END IF
       CASE(REALDATA3D)
          IF (ASSOCIATED(pw%cr3d)) THEN
             IF (cp_sll_pw_get_length(pool%real3d_pw)<pool%max_cache) THEN
                IF (cp_debug) THEN
                   CALL cp_sll_pw_insert_ordered(pool%real3d_pw, el=pw,&
                        insert_equals=.FALSE., did_insert=did_insert)
                   CPASSERT(did_insert)
                ELSE
                   CALL cp_sll_pw_insert_el(pool%real3d_pw, el=pw)
                END IF
             ELSE
                CALL cp_assert(max_max_cache<0,cp_warning_level,cp_assertion_failed,&
                     routineP,"hit max_cache"//&
CPSourceFileRef)
                CALL pw_release(pw)
             END IF
          ELSE
             IF (debug_this_module) THEN
                WRITE (unit=cp_logger_get_default_unit_nr(logger,local=.TRUE.),&
                     fmt="(' *** pw_pool ',i4,' pw ',i4,' cr3d is not associated, discarding')")&
                     pool%id_nr, pw%id_nr
                CALL print_stack(cp_logger_get_default_unit_nr(logger,local=.TRUE.))
             END IF
             CPASSERT(my_accept_non_compatible)
             CALL pw_release(pw)
          END IF
       CASE(COMPLEXDATA1D)
          IF (cp_sll_pw_get_length(pool%complex1d_pw)<pool%max_cache) THEN
             IF (cp_debug) THEN
                CALL cp_sll_pw_insert_ordered(pool%complex1d_pw, el=pw,&
                     insert_equals=.FALSE., did_insert=did_insert)
                CPASSERT(did_insert)
             ELSE
                CALL cp_sll_pw_insert_el(pool%complex1d_pw, el=pw)
             END IF
          ELSE
             CALL cp_assert(max_max_cache<0,cp_warning_level,cp_assertion_failed,&
                  routineP,"hit max_cache"//&
CPSourceFileRef)
             CALL pw_release(pw)
          END IF
       CASE(COMPLEXDATA3D)
          IF (cp_sll_pw_get_length(pool%complex3d_pw)<pool%max_cache) THEN
             IF (cp_debug) THEN
                CALL cp_sll_pw_insert_ordered(pool%complex3d_pw, el=pw,&
                     insert_equals=.FALSE., did_insert=did_insert)
                CPASSERT(did_insert)
             ELSE
                CALL cp_sll_pw_insert_el(pool%complex3d_pw, el=pw)
             END IF
          ELSE
             CALL cp_assert(max_max_cache<0,cp_warning_level,cp_assertion_failed,&
                  routineP,"hit max_cache"//&
CPSourceFileRef)
             CALL pw_release(pw)
          END IF
       CASE default
          ! unknown in_use
          CPASSERT(.FALSE.)
       END SELECT
       IF (ASSOCIATED(pw)) pw%ref_count=0
       !FM so that if someone tries to use a pw that is in the pool
       !FM (s)he gets problems
    END IF
    NULLIFY(pw)
    CALL timestop(handle)
  END SUBROUTINE pw_pool_give_back_pw

! *****************************************************************************
!> \brief returns a 3d real array of coefficents as the one used by pw with
!>      REALDATA3D, allocating it if none is present in the pool
!> \param pw_pool the pool that caches the cr3d
!> \param cr3d the pointer that will contain the array
!> \par History
!>      11.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
SUBROUTINE pw_pool_create_cr3d(pw_pool,cr3d)
    TYPE(pw_pool_type), POINTER              :: pw_pool
    REAL(kind=dp), DIMENSION(:, :, :), &
      POINTER                                :: cr3d

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_pool_create_cr3d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(pw_type), POINTER                   :: pw

  failure=.FALSE.
  NULLIFY(pw)
  logger => cp_get_default_logger()

  CPASSERT(ASSOCIATED(pw_pool))
  CPASSERT(pw_pool%ref_count>0)
  CPASSERT(.NOT.ASSOCIATED(cr3d))
  IF (ASSOCIATED(pw_pool%real3d_array)) THEN
     cr3d => cp_sll_3d_r_get_first_el(pw_pool%real3d_array)
     CALL cp_sll_3d_r_rm_first_el(pw_pool%real3d_array)
  ELSE
     pw => try_pop(pw_pool%real3d_pw)
     IF (ASSOCIATED(pw)) THEN
        CPASSERT(pw%ref_count==0)
        pw%ref_count=1
        cr3d => pw%cr3d
        NULLIFY(pw%cr3d)
        CALL pw_release(pw)
     END IF
  END IF
  IF (.NOT.ASSOCIATED(cr3d)) THEN
     ALLOCATE(cr3d(pw_pool%pw_grid%bounds_local(1,1):pw_pool%pw_grid%bounds_local(2,1),&
          pw_pool%pw_grid%bounds_local(1,2):pw_pool%pw_grid%bounds_local(2,2),&
          pw_pool%pw_grid%bounds_local(1,3):pw_pool%pw_grid%bounds_local(2,3)),&
          stat=stat)
     CPASSERT(stat==0)
     IF (debug_this_module) THEN
        WRITE (unit=cp_logger_get_default_unit_nr(logger,local=.TRUE.),&
             fmt="(' *** pw_pool ',i4,' created cr3d')") pw_pool%id_nr
        CALL print_stack(cp_logger_get_default_unit_nr(logger,local=.TRUE.))
     END IF
  ELSEIF (debug_this_module) THEN
        WRITE (unit=cp_logger_get_default_unit_nr(logger,local=.TRUE.),&
             fmt="(' *** pw_pool ',i4,' created cr3d')") pw_pool%id_nr
       CALL print_stack(cp_logger_get_default_unit_nr(logger,local=.TRUE.))
  END IF
END SUBROUTINE pw_pool_create_cr3d

! *****************************************************************************
!> \brief returns a 3d real array of coefficents as the one used by pw with
!>      REALDATA3D, allocating it if none is present in the pool
!> \param pw_pool the pool that caches the cr3d
!> \param cr3d the pointer that will contain the array
!> \param accept_non_compatible if true deallocates the non compatible
!>        arrays passed in, if false (the default) stops with an error
!> \par History
!>      11.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
SUBROUTINE pw_pool_give_back_cr3d(pw_pool,cr3d,accept_non_compatible)
    TYPE(pw_pool_type), POINTER              :: pw_pool
    REAL(kind=dp), DIMENSION(:, :, :), &
      POINTER                                :: cr3d
    LOGICAL, INTENT(in), OPTIONAL            :: accept_non_compatible

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_pool_give_back_cr3d', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: compatible, failure, &
                                                my_accept_non_compatible
    TYPE(cp_logger_type), POINTER            :: logger

  failure=.FALSE.
  my_accept_non_compatible=.FALSE.
  logger => cp_get_default_logger()
  IF (PRESENT(accept_non_compatible)) my_accept_non_compatible=accept_non_compatible

  CPASSERT(ASSOCIATED(pw_pool))
  CPASSERT(pw_pool%ref_count>0)
  IF (ASSOCIATED(cr3d)) THEN
     IF (debug_this_module) THEN
        WRITE (unit=cp_logger_get_default_unit_nr(logger,local=.TRUE.),&
             fmt="(' *** pw_pool ',i4,' received back a cr3d')") pw_pool%id_nr
        CALL print_stack(cp_logger_get_default_unit_nr(logger,local=.TRUE.))
     END IF
     compatible=ALL( MERGE( pw_pool%pw_grid%bounds_local(1,:)==LBOUND(cr3d).AND. &
                            pw_pool%pw_grid%bounds_local(2,:)==UBOUND(cr3d), &
                     pw_pool%pw_grid%bounds_local(2,:)<pw_pool%pw_grid%bounds_local(1,:), &
                     UBOUND(cr3d)>=LBOUND(cr3d) ) )
     CPASSERT(compatible.OR.my_accept_non_compatible)
     IF (compatible) THEN
        IF (cp_sll_3d_r_get_length(pw_pool%real3d_array)<pw_pool%max_cache) THEN
           CALL cp_sll_3d_r_insert_el(pw_pool%real3d_array, el=cr3d)
        ELSE
           CALL cp_assert(max_max_cache<0,cp_warning_level,cp_assertion_failed,&
                routineP,"hit max_cache"//&
CPSourceFileRef)
           DEALLOCATE(cr3d)
        END IF
     ELSE
        IF (debug_this_module) THEN
           WRITE (unit=cp_logger_get_default_unit_nr(logger,local=.TRUE.),&
                fmt="(' *** pw_pool ',i4,' cr3d not accepted, deallocating')")&
             pw_pool%id_nr
           CALL print_stack(cp_logger_get_default_unit_nr(logger,local=.TRUE.))
        END IF
        DEALLOCATE(cr3d)
     END IF
  ELSE
     IF (debug_this_module) THEN
        WRITE (unit=cp_logger_get_default_unit_nr(logger,local=.TRUE.),&
             fmt="(' *** pw_pool ',i4,' received back a null cr3d')")&
             pw_pool%id_nr
        CALL print_stack(cp_logger_get_default_unit_nr(logger,local=.TRUE.))
     END IF
     CPASSERT(my_accept_non_compatible)
  END IF
  NULLIFY(cr3d)
END SUBROUTINE pw_pool_give_back_cr3d

! *****************************************************************************
!> \brief creates a multigrid structure
!> \param pools the multigrid pool (i.e. an array of pw_pool)
!> \param pws the multigrid of coefficent you want to initialize
!> \param use_data which data it uses: REALDATA1D, COMPLEXDATA1D,
!>        REALDATA3D, COMPLEXDATA3D
!> \param in_space ...
!> \par History
!>      07.2004 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE pw_pools_create_pws(pools, pws, use_data, in_space)
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pools
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: pws
    INTEGER, INTENT(in)                      :: use_data
    INTEGER, INTENT(in), OPTIONAL            :: in_space

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_pools_create_pws', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure

  failure=.FALSE.

  CPASSERT(ASSOCIATED(pools))
  ALLOCATE(pws(SIZE(pools)))
  DO i=1,SIZE(pools)
     NULLIFY(pws(i)%pw)
     CALL pw_pool_create_pw(pools(i)%pool,pws(i)%pw,use_data,&
          in_space=in_space)
  END DO
END SUBROUTINE pw_pools_create_pws

! *****************************************************************************
!> \brief returns the pw part of the coefficents into the pools
!> \param pools the pools that will cache the pws %pw
!> \param pws the coefficents to give back
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE pw_pools_give_back_pws(pools, pws)
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pools
    TYPE(pw_p_type), DIMENSION(:), POINTER   :: pws

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_pools_give_back_pws', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure

  failure=.FALSE.

  CPASSERT(ASSOCIATED(pools))
  CPASSERT(ASSOCIATED(pws))
  CPASSERT(SIZE(pws)==SIZE(pools))
  DO i=1,SIZE(pools)
     CALL pw_pool_give_back_pw(pools(i)%pool,pws(i)%pw)
  END DO
  DEALLOCATE(pws)
END SUBROUTINE pw_pools_give_back_pws

! *****************************************************************************
!> \brief copies a multigrid pool, the underlying pools are shared
!> \param source_pools the pools to copy
!> \param target_pools will hold the copy of the pools
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE pw_pools_copy(source_pools, target_pools)
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: source_pools, target_pools

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_pools_copy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure

  failure=.FALSE.

  CPASSERT(ASSOCIATED(source_pools))
  ALLOCATE(target_pools(SIZE(source_pools)))
  DO i=1,SIZE(source_pools)
     target_pools(i)%pool => source_pools(i)%pool
     CALL pw_pool_retain(source_pools(i)%pool)
  END DO
END SUBROUTINE pw_pools_copy

! *****************************************************************************
!> \brief deallocates the given pools (releasing each of the underlying
!>      pools)
!> \param pools the pols to deallocate
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE pw_pools_dealloc(pools)
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pools

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_pools_dealloc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure

  failure=.FALSE.

  IF (ASSOCIATED(pools)) THEN
     DO i=1,SIZE(pools)
        CALL pw_pool_release(pools(i)%pool)
     END DO
     DEALLOCATE(pools)
  END IF
  NULLIFY(pools)
END SUBROUTINE pw_pools_dealloc

! *****************************************************************************
!> \brief deallocates all the cached grids
!> \param pools the pools to flush
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
SUBROUTINE pw_pools_flush_cache(pools)
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pools

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_pools_flush_cache', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure

  failure=.FALSE.

  CPASSERT(ASSOCIATED(pools))
  DO i=1,SIZE(pools)
     CALL pw_pool_flush_cache(pools(i)%pool)
  END DO
END SUBROUTINE pw_pools_flush_cache

! *****************************************************************************
!> \brief writes the actual contents of the pw_pool
!> \param pw_pool ...
!> \param unit_nr ...
!> \par History
!>      02.2003 created [fawzi]
!> \author fawzi
! *****************************************************************************
SUBROUTINE pw_pool_write(pw_pool,unit_nr)
    TYPE(pw_pool_type), POINTER              :: pw_pool
    INTEGER                                  :: unit_nr

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_pool_write', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(cp_sll_pw_type), POINTER            :: iter
    TYPE(pw_type), POINTER                   :: pw_att

  failure=.FALSE.

  IF (ASSOCIATED(pw_pool)) THEN
     WRITE (unit=unit_nr, &
          fmt="(' <pw_pool>{ id_nr=',i8,', ref_count=',i8,', max_cache=',i8,',')")&
          pw_pool%id_nr,pw_pool%ref_count,pw_pool%max_cache

     WRITE (unit=unit_nr, fmt="('real1d_pw=(')",advance="no")
     iter => pw_pool%real1d_pw
     DO WHILE(cp_sll_pw_next(iter,el_att=pw_att))
        WRITE (unit=unit_nr, fmt="(i8,',')",advance="no") pw_att%id_nr
     END DO
     WRITE (unit=unit_nr, fmt="('),')")

     WRITE (unit=unit_nr, fmt="('real3d_pw=(')",advance="no")
     iter => pw_pool%real3d_pw
     DO WHILE(cp_sll_pw_next(iter,el_att=pw_att))
        WRITE (unit=unit_nr, fmt="(i8,',')",advance="no") pw_att%id_nr
     END DO
     WRITE (unit=unit_nr, fmt="('),')")

     WRITE (unit=unit_nr, fmt="('complex1d_pw=(')",advance="no")
     iter => pw_pool%complex1d_pw
     DO WHILE(cp_sll_pw_next(iter,el_att=pw_att))
        WRITE (unit=unit_nr, fmt="(i8,',')",advance="no") pw_att%id_nr
     END DO
     WRITE (unit=unit_nr, fmt="('),')")

     WRITE (unit=unit_nr, fmt="('complex3d_pw=(')",advance="no")
     iter => pw_pool%complex3d_pw
     DO WHILE(cp_sll_pw_next(iter,el_att=pw_att))
        WRITE (unit=unit_nr, fmt="(i8,',')",advance="no") pw_att%id_nr
     END DO
     WRITE (unit=unit_nr, fmt="(')')")

     WRITE (unit=unit_nr, fmt="('}')")
  ELSE
     WRITE (unit=unit_nr, fmt="('<pw_pool *null*>')")
  END IF
END SUBROUTINE pw_pool_write

END MODULE pw_pool_types
