!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/pw_pool_types [1.0] *
!!
!!   NAME
!!     pw_pool_types
!!
!!   FUNCTION
!!     Manages a pool of grids (to be used for example as tmp objects),
!!     but can also be used to instantiate grids that are never given back.
!!
!!     Multigrid pools are just an array of pw_pools
!!
!!   NOTES
!!     The pool could also work without pointers (doing = each time),
!!     but I find it *very* ugly.
!!
!!     The pool could be integrated into pw_grid_type, I don't know if
!!     it would be a good or bad idea (but would add a circular dependence
!!     between pw and pw_grid types).
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!   SOURCE
!****************************************************************************
MODULE pw_pool_types
  USE coefficient_types,               ONLY: PLANEWAVES,&
                                             coeff_type,&
                                             NOTHING
  USE cp_error_handling,               ONLY: cp_a_l,&
                                             cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_get_logger,&
                                             cp_error_message,&
                                             cp_error_type
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
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_logger_get_default_unit_nr,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE kinds,                           ONLY: dp
  USE machine,                         ONLY: m_loc_r
  USE pw_grid_types,                   ONLY: pw_grid_type
  USE pw_grids,                        ONLY: pw_grid_release,&
                                             pw_grid_retain
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             COMPLEXDATA3D,&
                                             REALDATA1D,&
                                             REALDATA3D,&
                                             pw_create,&
                                             pw_release,&
                                             pw_type,&
                                             pw_write, pw_p_type
  USE timings,                         ONLY: print_stack,&
                                             timeset,&
                                             timestop

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.FALSE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='pw_pool_types'
  INTEGER, SAVE, PRIVATE :: last_pw_pool_id_nr=0
  INTEGER, PARAMETER :: default_max_cache=20, max_max_cache=50

  PUBLIC :: pw_pool_type, pw_pool_p_type
  PUBLIC :: pw_pool_create, pw_pool_retain, pw_pool_release,&
       pw_pool_create_pw, pw_pool_give_back_pw, pw_pool_give_back_coeff,&
       pw_pool_init_coeff, pw_pool_flush_cache, pw_pool_write,&
       pw_pool_create_cr3d, pw_pool_give_back_cr3d
  PUBLIC :: pw_pools_init_coeffs, pw_pools_give_back_coeffs,&
       pw_pools_copy, pw_pools_dealloc, pw_pools_flush_cache,&
       pw_pools_create_pws, pw_pools_give_back_pws
!!***
!****************************************************************************


!!****s* pw_pool_types/pw_pool_type [1.0] *
!!
!!   NAME
!!    pw_pool_type
!!
!!   FUNCTION
!!     Manages a pool of grids (to be used for example as tmp objects),
!!     but can also be used to intantiate grids that are never given back.
!!     
!!
!!   NOTES
!!     As of now I would like replace the linked lists by arrays
!!     (no annoying list elements that are allocated would show up when
!!     tracking leaks) [fawzi]
!!
!!   ATTRIBUTES
!!     - ref_count: reference count (see /cp2k/doc/ReferenceCounting.html)
!!     - id_nr: number that identifies each pool
!!     - real1d_pw, real3d_pw, complex1d_pw, complex3d_pw: liked list with
!!       the cached grids of the corresponding type
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!   SOURCE
!***************************************************************************
  TYPE pw_pool_type
     INTEGER :: ref_count, id_nr, max_cache
     TYPE(pw_grid_type), POINTER :: pw_grid
     TYPE(cp_sll_pw_type), POINTER :: real1d_pw, real3d_pw,&
          complex1d_pw, complex3d_pw
     TYPE(cp_sll_3d_r_type), POINTER :: real3d_array
  END TYPE pw_pool_type
!!***
!****************************************************************************

!!****s* pw_pool_types/pw_pool_p_type [1.0] *
!!
!!   NAME
!!     pw_pool_p_type
!!
!!   FUNCTION
!!     to create arrays of pools
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     - pool: the pool
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!   SOURCE
!***************************************************************************
  TYPE pw_pool_p_type
     TYPE(pw_pool_type), POINTER :: pool
  END TYPE pw_pool_p_type
!!***
!****************************************************************************

CONTAINS

!!****f* pw_pool_types/pw_pool_create [1.0] *
!!
!!   NAME
!!     pw_pool_create
!!
!!   FUNCTION
!!     creates a pool for pw
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pw_pool: the pool to create
!!     - pw_grid: the grid that is used to create the pw
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE pw_pool_create(pool, pw_grid, max_cache, error)
    TYPE(pw_pool_type), POINTER              :: pool
    TYPE(pw_grid_type), POINTER              :: pw_grid
    INTEGER, OPTIONAL                        :: max_cache
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_pool_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger

    failure=.FALSE.
    logger => cp_error_get_logger(error)

    ALLOCATE(pool, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       pool%pw_grid => pw_grid
       CALL pw_grid_retain(pw_grid,error=error)
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
    END IF
  END SUBROUTINE pw_pool_create
!***************************************************************************

!!****f* pw_pool_types/pw_pool_retain [1.0] *
!!
!!   NAME
!!     pw_pool_retain
!!
!!   FUNCTION
!!     retains the pool (see cp2k/doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pool: the pool to retain
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE pw_pool_retain(pool,error)
    TYPE(pw_pool_type), POINTER              :: pool
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_pool_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger

    failure=.FALSE.
    logger => cp_error_get_logger(error)

    CPPrecondition(ASSOCIATED(pool),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPrecondition(pool%ref_count>0,cp_failure_level,routineP,error,failure)

       pool%ref_count=pool%ref_count+1
       IF (debug_this_module) THEN
          WRITE (unit=cp_logger_get_default_unit_nr(logger,local=.TRUE.),&
               fmt="(' *** pw_pool ',i4,' has been retained, ref_count=',i4)")&
               pool%id_nr, pool%ref_count
          CALL print_stack(cp_logger_get_default_unit_nr(logger,local=.TRUE.))
       END IF
    END IF
  END SUBROUTINE pw_pool_retain
!***************************************************************************

!!****f* pw_pool_types/pw_pool_flush_cache [1.0] *
!!
!!   NAME
!!     pw_pool_flush_cache
!!
!!   FUNCTION
!!     deallocates all the cached grids
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pool: the pool to flush
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE pw_pool_flush_cache(pool, error)
    TYPE(pw_pool_type), POINTER              :: pool
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_pool_flush_cache', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure
    REAL(kind=dp), DIMENSION(:, :, :), &
      POINTER                                :: array_att
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(cp_sll_3d_r_type), POINTER          :: array_iterator
    TYPE(cp_sll_pw_type), POINTER            :: iterator
    TYPE(pw_type), POINTER                   :: pw_el

    failure=.FALSE.
 

    CPPrecondition(ASSOCIATED(pool),cp_failure_level,routineP,error,failure)
    CPPrecondition(pool%ref_count>0,cp_failure_level,routineP,error,failure)
    NULLIFY(iterator,array_iterator,pw_el,array_att)
    logger => cp_error_get_logger(error)
    IF (.NOT.failure) THEN
       IF (debug_this_module) THEN
          WRITE (unit=cp_logger_get_default_unit_nr(logger,local=.TRUE.),&
               fmt="(' *** pw_pool ',i4,' is flushing the cache')") pool%id_nr
          CALL print_stack(cp_logger_get_default_unit_nr(logger,local=.TRUE.))
       END IF

       iterator => pool%real1d_pw
       DO
          IF (.NOT.cp_sll_pw_next(iterator,el_att=pw_el,error=error)) EXIT
          CPPrecondition(pw_el%ref_count==0,cp_failure_level,routineP,error,failure)
          pw_el%ref_count=1
          CALL pw_release(pw_el,error=error)
       END DO
       CALL cp_sll_pw_dealloc(pool%real1d_pw,error=error)

       iterator => pool%real3d_pw
       DO
          IF (.NOT.cp_sll_pw_next(iterator,el_att=pw_el,error=error)) EXIT
          CPPrecondition(pw_el%ref_count==0,cp_failure_level,routineP,error,failure)
          pw_el%ref_count=1
          CALL pw_release(pw_el, error=error)
       END DO
       CALL cp_sll_pw_dealloc(pool%real3d_pw,error=error)

       iterator => pool%complex1d_pw
       DO
          IF (.NOT.cp_sll_pw_next(iterator,el_att=pw_el,error=error)) EXIT
          CPPrecondition(pw_el%ref_count==0,cp_failure_level,routineP,error,failure)
          pw_el%ref_count=1
          CALL pw_release(pw_el, error=error)
       END DO
       CALL cp_sll_pw_dealloc(pool%complex1d_pw,error=error)

       iterator => pool%complex3d_pw
       DO
          IF (.NOT.cp_sll_pw_next(iterator,el_att=pw_el,error=error)) EXIT
          CPPrecondition(pw_el%ref_count==0,cp_failure_level,routineP,error,failure)
          pw_el%ref_count=1
          CALL pw_release(pw_el, error=error)
       END DO
       CALL cp_sll_pw_dealloc(pool%complex3d_pw)

       array_iterator => pool%real3d_array
       DO
          IF (.NOT.cp_sll_3d_r_next(array_iterator,el_att=array_att,&
               error=error)) EXIT
          DEALLOCATE(array_att,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END DO
       CALL cp_sll_3d_r_dealloc(pool%real3d_array)

    END IF
  END SUBROUTINE pw_pool_flush_cache
!***************************************************************************

!!****f* pw_pool_types/pw_pool_release [1.0] *
!!
!!   NAME
!!     pw_pool_release
!!
!!   FUNCTION
!!     releases the given pool (see cp2k/doc/ReferenceCounting.html)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pool: the pool to release
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE pw_pool_release(pool,error)
    TYPE(pw_pool_type), POINTER              :: pool
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_pool_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger

    failure=.FALSE.
    logger => cp_error_get_logger(error)

    IF (ASSOCIATED(pool)) THEN
       CPPrecondition(pool%ref_count>0,cp_failure_level,routineP,error,failure)
       pool%ref_count=pool%ref_count-1
       IF (debug_this_module) THEN
          WRITE (unit=cp_logger_get_default_unit_nr(logger,local=.TRUE.),&
               fmt="(' *** pw_pool ',i4,' released ref_count=',i4)") &
               pool%id_nr, pool%ref_count
          CALL print_stack(cp_logger_get_default_unit_nr(logger,local=.TRUE.))
       END IF
       IF (pool%ref_count==0) THEN
          pool%ref_count=1
          CALL pw_pool_flush_cache(pool,error=error)
          pool%ref_count=0
          CPPrecondition(ASSOCIATED(pool%pw_grid),cp_warning_level,routineP,error,failure)
          IF (.NOT.failure) THEN
             CALL pw_grid_release(pool%pw_grid,error=error)
          END IF

          DEALLOCATE(pool,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
    END IF
    NULLIFY(pool)
  END SUBROUTINE pw_pool_release
!***************************************************************************

!!****f* pw_pool_types/try_pop [1.0] *
!!
!!   NAME
!!     try_pop
!!
!!   FUNCTION
!!     tries to pop an element from the given list (no error on failure)
!!
!!   NOTES
!!     private function
!!
!!   INPUTS
!!     - list: the list to pop
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
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
!***************************************************************************

!!****f* pw_pool_types/pw_pool_create_pw [1.0] *
!!
!!   NAME
!!     pw_pool_create_pw
!!
!!   FUNCTION
!!     returns a pw, allocating it if none is in the pool
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pool: the pool from where you get the pw
!!     - pw: will contain the new pw
!!     - use_data: which data it uses: REALDATA1D, COMPLEXDATA1D,
!!       REALDATA3D, COMPLEXDATA3D
!!     - in_space: in which space it is: REALSPACE, RECIPROCALSPACE
!!     - mass_flag: whether to allocate mass coeff (mass_cr or mass_cc)
!!       defaults to false
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE pw_pool_create_pw(pool, pw, use_data, in_space, mass_flag, error)
    TYPE(pw_pool_type), POINTER              :: pool
    TYPE(pw_type), POINTER                   :: pw
    INTEGER, INTENT(in)                      :: use_data
    INTEGER, INTENT(in), OPTIONAL            :: in_space
    LOGICAL, INTENT(in), OPTIONAL            :: mass_flag
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_pool_create_pw', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    LOGICAL                                  :: failure
    REAL(kind=dp), DIMENSION(:, :, :), &
      POINTER                                :: cr3d_ptr
    TYPE(cp_logger_type), POINTER            :: logger

    failure=.FALSE.

    CALL timeset(routineN,"I","",handle)
    NULLIFY(pw)
    NULLIFY(cr3d_ptr)
    logger => cp_error_get_logger(error)
    CPPrecondition(ASSOCIATED(pool),cp_failure_level,routineP,error,failure)
    CPPrecondition(pool%ref_count>0,cp_failure_level,routineP,error,failure)

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
       CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT

    IF (.NOT.ASSOCIATED(pw)) THEN
       CALL pw_create(pw, pool%pw_grid, use_data=use_data, &
            mass_flag=mass_flag, cr3d_ptr=cr3d_ptr)
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
               pw=pw, error=error)
          CALL print_stack(cp_logger_get_default_unit_nr(logger,local=.TRUE.))
       END IF
    ELSE
       CPPrecondition(pw%ref_count==0,cp_failure_level,routineP,error,failure)
       pw%ref_count=1
       IF (debug_this_module) THEN
          WRITE(unit=cp_logger_get_default_unit_nr(logger,local=.TRUE.),&
               fmt="(' *** pw_pool ',i4,' created pw reusing old ',i4)")&
               pool%id_nr,pw%id_nr
          CALL pw_write(&
               unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.),&
               pw=pw, error=error)
          CALL print_stack(cp_logger_get_default_unit_nr(logger,local=.TRUE.))
       END IF
    END IF

    IF (PRESENT(mass_flag)) THEN
       IF (mass_flag) THEN
          CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,&
               routineP,"mass not cached in pool"//&
               CPSourceFileRef,&
               error=error)
          
          IF (use_data==REALDATA1D) THEN
             IF (.NOT.ASSOCIATED(pw%mass_cr)) THEN
                ALLOCATE(pw%mass_cr(pw % pw_grid % ngpts_cut_local),&
                     stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
          ELSE IF (use_data==COMPLEXDATA1D) THEN
             IF (.NOT.ASSOCIATED(pw%mass_cc)) THEN
                ALLOCATE(pw%mass_cc(pw % pw_grid % ngpts_cut_local),&
                     stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
          END IF
       END IF
    END IF
    
    pw%in_space=0
    IF (PRESENT(in_space)) pw%in_space=in_space

    CALL timestop(0.0_dp,handle)

  END SUBROUTINE pw_pool_create_pw
!***************************************************************************

!!****f* pw_pool_types/pw_pool_give_back_pw [1.0] *
!!
!!   NAME
!!     pw_pool_give_back_pw
!!
!!   FUNCTION
!!     returns the pw to the pool
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pool: the pool where to reintegrate the pw
!!     - pw: the pw to give back
!!     - accept_non_compatible: if non compatible pw should be accepted
!!       (they will be destroied). Defaults to false (and thus stops with
!!       an error)
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE pw_pool_give_back_pw(pool, pw, accept_non_compatible, error)
    TYPE(pw_pool_type), POINTER              :: pool
    TYPE(pw_type), POINTER                   :: pw
    LOGICAL, INTENT(in), OPTIONAL            :: accept_non_compatible
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_pool_give_back_pw', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    LOGICAL                                  :: did_insert, failure, &
                                                my_accept_non_compatible
    TYPE(cp_logger_type), POINTER            :: logger

    failure=.FALSE.
    my_accept_non_compatible=.FALSE.
    logger => cp_error_get_logger(error)
    IF (PRESENT(accept_non_compatible)) my_accept_non_compatible=accept_non_compatible

    CALL timeset(routineN,"I","",handle)
    CPPrecondition(ASSOCIATED(pool),cp_failure_level,routineP,error,failure)
    CPPrecondition(pool%ref_count>0,cp_failure_level,routineP,error,failure)
    IF (.NOT.ASSOCIATED(pw)) THEN
       CPPrecondition(my_accept_non_compatible,cp_warning_level,routineP,error,failure)
       failure=.TRUE.
    END IF
    IF (.NOT. failure) THEN
       CPPrecondition(pw%ref_count==1,cp_failure_level,routineP,error,failure)
       IF (pw%pw_grid%identifier/=pool%pw_grid%identifier) THEN
          IF (debug_this_module) THEN
             WRITE (unit=cp_logger_get_default_unit_nr(logger,local=.TRUE.),&
                  fmt="(' *** pw_pool ',i4,' giving back incompatible pw ',i4)")&
                  pool%id_nr, pw%id_nr
             CALL pw_write(unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.),&
                  pw=pw,error=error)
             CALL print_stack(cp_logger_get_default_unit_nr(logger,local=.TRUE.))
          END IF

          CALL cp_assert(my_accept_non_compatible,&
               cp_failure_level, cp_assertion_failed, routineP,&
               "pool cannot reuse pw of another grid "//&
               CPSourceFileRef,&
               error=error,failure=failure)
          CALL pw_release(pw,error=error)
          failure=.TRUE.
       END IF
    END IF

    IF (.NOT. failure) THEN
       IF (ASSOCIATED(pw%mass_cr)) THEN
          DEALLOCATE(pw%mass_cr,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          NULLIFY(pw%mass_cr)
       END IF
       IF (ASSOCIATED(pw%mass_cc)) THEN
          DEALLOCATE(pw%mass_cc,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
          NULLIFY(pw%mass_cc)
       END IF

       IF (debug_this_module) THEN
          WRITE (unit=cp_logger_get_default_unit_nr(logger,local=.TRUE.),&
               fmt="(' *** pw_pool ',i4,' giving back pw ',i4)")&
               pool%id_nr, pw%id_nr
          CALL pw_write(unit_nr=cp_logger_get_default_unit_nr(logger,local=.TRUE.),&
               pw=pw,error=error)
          CALL print_stack(cp_logger_get_default_unit_nr(logger,local=.TRUE.))
       END IF

       SELECT CASE(pw%in_use)
       CASE(REALDATA1D)
          IF (cp_sll_pw_get_length(pool%real1d_pw,error=error)<pool%max_cache) THEN
             IF (cp_debug) THEN
                CALL cp_sll_pw_insert_ordered(pool%real1d_pw, el=pw,&
                     insert_equals=.FALSE., did_insert=did_insert)
                CPAssert(did_insert,cp_failure_level,routineP,error,failure)
             ELSE
                CALL cp_sll_pw_insert_el(pool%real1d_pw, el=pw)
             END IF
          ELSE
             CALL cp_assert(max_max_cache<0,cp_warning_level,cp_assertion_failed,&
                  routineP,"hit max_cache"//&
                  CPSourceFileRef,&
                  error)
             CALL pw_release(pw,error=error)
          END IF
       CASE(REALDATA3D)
          IF (ASSOCIATED(pw%cr3d)) THEN
             IF (cp_sll_pw_get_length(pool%real3d_pw,error=error)<pool%max_cache) THEN
                IF (cp_debug) THEN
                   CALL cp_sll_pw_insert_ordered(pool%real3d_pw, el=pw,&
                        insert_equals=.FALSE., did_insert=did_insert)
                   CPAssert(did_insert,cp_failure_level,routineP,error,failure)
                ELSE
                   CALL cp_sll_pw_insert_el(pool%real3d_pw, el=pw, error=error)
                END IF
             ELSE
                CALL cp_assert(max_max_cache<0,cp_warning_level,cp_assertion_failed,&
                     routineP,"hit max_cache"//&
                     CPSourceFileRef,&
                     error)
                CALL pw_release(pw,error=error)
             END IF
          ELSE
             IF (debug_this_module) THEN
                WRITE (unit=cp_logger_get_default_unit_nr(logger,local=.TRUE.),&
                     fmt="(' *** pw_pool ',i4,' pw ',i4,' cr3d is not associated, discarding')")&
                     pool%id_nr, pw%id_nr
                CALL print_stack(cp_logger_get_default_unit_nr(logger,local=.TRUE.))
             END IF
             CPAssert(my_accept_non_compatible,cp_failure_level,routineP,error,failure)
             CALL pw_release(pw, error=error)
          END IF
       CASE(COMPLEXDATA1D)
          IF (cp_sll_pw_get_length(pool%complex1d_pw,error=error)<pool%max_cache) THEN
             IF (cp_debug) THEN
                CALL cp_sll_pw_insert_ordered(pool%complex1d_pw, el=pw,&
                     insert_equals=.FALSE., did_insert=did_insert)
                CPAssert(did_insert,cp_failure_level,routineP,error,failure)
             ELSE
                CALL cp_sll_pw_insert_el(pool%complex1d_pw, el=pw)
             END IF
          ELSE
             CALL cp_assert(max_max_cache<0,cp_warning_level,cp_assertion_failed,&
                  routineP,"hit max_cache"//&
                  CPSourceFileRef,&
                  error)
             CALL pw_release(pw,error=error)
          END IF
       CASE(COMPLEXDATA3D)
          IF (cp_sll_pw_get_length(pool%complex3d_pw,error=error)<pool%max_cache) THEN
             IF (cp_debug) THEN
                CALL cp_sll_pw_insert_ordered(pool%complex3d_pw, el=pw,&
                     insert_equals=.FALSE., did_insert=did_insert)
                CPAssert(did_insert,cp_failure_level,routineP,error,failure)
             ELSE
                CALL cp_sll_pw_insert_el(pool%complex3d_pw, el=pw)
             END IF
          ELSE
             CALL cp_assert(max_max_cache<0,cp_warning_level,cp_assertion_failed,&
                  routineP,"hit max_cache"//&
                  CPSourceFileRef,&
                  error)
             CALL pw_release(pw,error=error)
          END IF
       CASE default
          ! unknown in_use
          CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
       END SELECT
       IF (ASSOCIATED(pw)) pw%ref_count=0
       !FM so that if someone tries to use a pw that is in the pool
       !FM (s)he gets problems
    END IF
    NULLIFY(pw)
    CALL timestop(0.0_dp, handle)
  END SUBROUTINE pw_pool_give_back_pw
!***************************************************************************

!!****f* pw_pool_types/pw_pool_create_cr3d [1.0] *
!!
!!   NAME
!!     pw_pool_create_cr3d
!!
!!   FUNCTION
!!     returns a 3d real array of coefficents as the one used by pw with
!!     REALDATA3D, allocating it if none is present in the pool
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - pw_pool: the pool that caches the cr3d
!!     - cr3d: the pointer that will contain the array
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     11.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE pw_pool_create_cr3d(pw_pool,cr3d,error)
    TYPE(pw_pool_type), POINTER              :: pw_pool
    REAL(kind=dp), DIMENSION(:, :, :), &
      POINTER                                :: cr3d
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_pool_create_cr3d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: failure
    TYPE(cp_logger_type), POINTER            :: logger
    TYPE(pw_type), POINTER                   :: pw

  failure=.FALSE.
  NULLIFY(pw)
  logger => cp_error_get_logger(error)
  
  CPPrecondition(ASSOCIATED(pw_pool),cp_failure_level,routineP,error,failure)
  CPPrecondition(pw_pool%ref_count>0,cp_failure_level,routineP,error,failure)
  CPPrecondition(.NOT.ASSOCIATED(cr3d),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     IF (ASSOCIATED(pw_pool%real3d_array)) THEN
        cr3d => cp_sll_3d_r_get_first_el(pw_pool%real3d_array)
        CALL cp_sll_3d_r_rm_first_el(pw_pool%real3d_array)
     ELSE
        pw => try_pop(pw_pool%real3d_pw)
        IF (ASSOCIATED(pw)) THEN
           CPPrecondition(pw%ref_count==0,cp_failure_level,routineP,error,failure)
           pw%ref_count=1
           cr3d => pw%cr3d
           NULLIFY(pw%cr3d)
           CALL pw_release(pw, error=error)
        END IF
     END IF
     IF (.NOT.ASSOCIATED(cr3d)) THEN
        ALLOCATE(cr3d(pw_pool%pw_grid%bounds_local(1,1):pw_pool%pw_grid%bounds_local(2,1),&
             pw_pool%pw_grid%bounds_local(1,2):pw_pool%pw_grid%bounds_local(2,2),&
             pw_pool%pw_grid%bounds_local(1,3):pw_pool%pw_grid%bounds_local(2,3)),&
             stat=stat)
        CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
        IF (debug_this_module) THEN
           WRITE (unit=cp_logger_get_default_unit_nr(logger,local=.TRUE.),&
                fmt="(' *** pw_pool ',i4,' created cr3d at 0x',z16.16)")&
                pw_pool%id_nr, m_loc_r(cr3d)
           CALL print_stack(cp_logger_get_default_unit_nr(logger,local=.TRUE.))
        END IF
     ELSEIF (debug_this_module) THEN
           WRITE (unit=cp_logger_get_default_unit_nr(logger,local=.TRUE.),&
                fmt="(' *** pw_pool ',i4,' created cr3d reusing the one at 0x',z16.16)")&
                pw_pool%id_nr, m_loc_r(cr3d)
          CALL print_stack(cp_logger_get_default_unit_nr(logger,local=.TRUE.))
     END IF
  END IF
END SUBROUTINE pw_pool_create_cr3d
!***************************************************************************
!!****f* pw_pool_types/pw_pool_give_back_cr3d [1.0] *
!!
!!   NAME
!!     pw_pool_give_back_cr3d
!!
!!   FUNCTION
!!     returns a 3d real array of coefficents as the one used by pw with
!!     REALDATA3D, allocating it if none is present in the pool
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - pw_pool: the pool that caches the cr3d
!!     - cr3d: the pointer that will contain the array
!!     - accept_non_compatible: if true deallocates the non compatible
!!       arrays passed in, if false (the default) stops with an error
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     11.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE pw_pool_give_back_cr3d(pw_pool,cr3d,accept_non_compatible,error)
    TYPE(pw_pool_type), POINTER              :: pw_pool
    REAL(kind=dp), DIMENSION(:, :, :), &
      POINTER                                :: cr3d
    LOGICAL, INTENT(in), OPTIONAL            :: accept_non_compatible
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_pool_give_back_cr3d', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: stat
    LOGICAL                                  :: compatible, failure, &
                                                my_accept_non_compatible
    TYPE(cp_logger_type), POINTER            :: logger

  failure=.FALSE.
  my_accept_non_compatible=.FALSE.
  logger => cp_error_get_logger(error)
  IF (PRESENT(accept_non_compatible)) my_accept_non_compatible=accept_non_compatible

  CPPrecondition(ASSOCIATED(pw_pool),cp_failure_level,routineP,error,failure)
  CPPrecondition(pw_pool%ref_count>0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     IF (ASSOCIATED(cr3d)) THEN
        IF (debug_this_module) THEN
           WRITE (unit=cp_logger_get_default_unit_nr(logger,local=.TRUE.),&
                fmt="(' *** pw_pool ',i4,' received back a cr3d at 0x',z16.16)")&
                pw_pool%id_nr, m_loc_r(cr3d)
           CALL print_stack(cp_logger_get_default_unit_nr(logger,local=.TRUE.))
        END IF
        compatible=ALL( MERGE( pw_pool%pw_grid%bounds_local(1,:)==LBOUND(cr3d).AND. &
                               pw_pool%pw_grid%bounds_local(2,:)==UBOUND(cr3d), &
                        pw_pool%pw_grid%bounds_local(2,:)<pw_pool%pw_grid%bounds_local(1,:), &
                        UBOUND(cr3d)>=LBOUND(cr3d) ) ) 
        CPPrecondition(compatible.OR.my_accept_non_compatible,cp_failure_level,routineP,error,failure)
        IF (compatible) THEN
           IF (cp_sll_3d_r_get_length(pw_pool%real3d_array,error=error)<pw_pool%max_cache) THEN
              CALL cp_sll_3d_r_insert_el(pw_pool%real3d_array, el=cr3d, error=error)
           ELSE
              CALL cp_assert(max_max_cache<0,cp_warning_level,cp_assertion_failed,&
                   routineP,"hit max_cache"//&
                   CPSourceFileRef,&
                   error)
              DEALLOCATE(cr3d, stat=stat)
              CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
           END IF
        ELSE
           IF (debug_this_module) THEN
              WRITE (unit=cp_logger_get_default_unit_nr(logger,local=.TRUE.),&
                   fmt="(' *** pw_pool ',i4,' cr3d not accepted, deallocating')")&
                pw_pool%id_nr
              CALL print_stack(cp_logger_get_default_unit_nr(logger,local=.TRUE.))
           END IF
           DEALLOCATE(cr3d, stat=stat)
           CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
        END IF
     ELSE
        IF (debug_this_module) THEN
           WRITE (unit=cp_logger_get_default_unit_nr(logger,local=.TRUE.),&
                fmt="(' *** pw_pool ',i4,' received back a null cr3d')")&
                pw_pool%id_nr
           CALL print_stack(cp_logger_get_default_unit_nr(logger,local=.TRUE.))
        END IF
        CPPrecondition(my_accept_non_compatible,cp_failure_level,routineP,error,failure)
     END IF
  END IF
  NULLIFY(cr3d)
END SUBROUTINE pw_pool_give_back_cr3d
!***************************************************************************

!!****f* pw_pool_types/pw_pool_init_coeff [1.0] *
!!
!!   NAME
!!     pw_pool_init_coeff
!!
!!   FUNCTION
!!     initializes a coefficent with the corresponding grid type
!!
!!   NOTES
!!     The coeff will contain random data, call coeff_zero if necessary
!!
!!   INPUTS
!!     - pool: the pool you get the coefficent from
!!     - coeff: the coefficent you want to initialize
!!     - use_data: which data it uses: REALDATA1D, COMPLEXDATA1D,
!!       REALDATA3D, COMPLEXDATA3D
!!     - mass_flag: whether to allocate mass coeff (mass_cr or mass_cc)
!!       defaults to false
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE pw_pool_init_coeff(pool, coeff, use_data, in_space, mass_flag,error)
    TYPE(pw_pool_type), POINTER              :: pool
    TYPE(coeff_type), INTENT(inout)          :: coeff
    INTEGER, INTENT(in)                      :: use_data
    INTEGER, INTENT(in), OPTIONAL            :: in_space
    LOGICAL, INTENT(in), OPTIONAL            :: mass_flag
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_pool_init_coeff', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(pool),cp_failure_level,routineP,error,failure)
  CPPrecondition(pool%ref_count>0,cp_failure_level,routineP,error,failure)
  
  IF (.NOT. failure) THEN
     NULLIFY(coeff%pw)
     CALL pw_pool_create_pw(pool,coeff%pw,use_data,mass_flag=mass_flag,&
          in_space=in_space, error=error)
     coeff%in_use=PLANEWAVES
  END IF
END SUBROUTINE pw_pool_init_coeff
!***************************************************************************

!!****f* pw_pool_types/pw_pool_give_back_coeff [1.0] *
!!
!!   NAME
!!     pw_pool_give_back_coeff
!!
!!   FUNCTION
!!     returns the pw part of the coeff to the pool
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pool: the pool that will cache coeff%pw
!!     - coeff: the coefficent to give back
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE pw_pool_give_back_coeff(pool,coeff,accept_non_compatible,error)
    TYPE(pw_pool_type), POINTER              :: pool
    TYPE(coeff_type), INTENT(INOUT)          :: coeff
    LOGICAL, INTENT(in), OPTIONAL            :: accept_non_compatible
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_pool_give_back_coeff', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

  failure=.FALSE.
  
  IF (coeff%in_use/=NOTHING) THEN
     CALL pw_pool_give_back_pw(pool,coeff%pw,&
          accept_non_compatible=accept_non_compatible,error=error)
     coeff%in_use=NOTHING
  END IF
END SUBROUTINE pw_pool_give_back_coeff
!***************************************************************************

!!****f* pw_pool_types/pw_pools_create_pws [1.0] *
!!
!!   NAME
!!     pw_pools_create_pws
!!
!!   FUNCTION
!!     creates a multigrid structure
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pools: the multigrid pool (i.e. an array of pw_pool)
!!     - pws: the multigrid of coefficent you want to initialize
!!     - use_data: which data it uses: REALDATA1D, COMPLEXDATA1D,
!!       REALDATA3D, COMPLEXDATA3D
!!     - mass_flag: whether to allocate mass coeff (mass_cr or mass_cc)
!!       defaults to false
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     07.2004 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE pw_pools_create_pws(pools, pws, use_data,&
     mass_flag, in_space, error)
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pools
    TYPE(pw_p_type), DIMENSION(:), &
      pointer                                :: pws
    INTEGER, INTENT(in)                      :: use_data
    LOGICAL, INTENT(in), OPTIONAL            :: mass_flag
    INTEGER, INTENT(in), OPTIONAL            :: in_space
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_pools_create_pws', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i,stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(pools),cp_failure_level,routineP,error,failure)
  ALLOCATE(pws(SIZE(pools)),stat=stat)
  CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     DO i=1,SIZE(pools)
        nullify(pws(i)%pw)
        CALL pw_pool_create_pw(pools(i)%pool,pws(i)%pw,use_data,&
             mass_flag=mass_flag,in_space=in_space,error=error)
     END DO
  END IF
END SUBROUTINE pw_pools_create_pws
!***************************************************************************

!!****f* pw_pool_types/pw_pools_give_back_pws [1.0] *
!!
!!   NAME
!!     pw_pools_give_back_pws
!!
!!   FUNCTION
!!     returns the pw part of the coefficents into the pools
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pools: the pools that will cache the pws %pw
!!     - pws: the coefficents to give back
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE pw_pools_give_back_pws(pools, pws, error)
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pools
    TYPE(pw_p_type), DIMENSION(:), &
      pointer                          :: pws
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_pools_give_back_pws', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i,stat
    LOGICAL                                  :: failure

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(pools),cp_failure_level,routineP,error,failure)
  CPPrecondition(associated(pws),cp_failure_level,routineP,error,failure)
  CPPrecondition(SIZE(pws)==SIZE(pools),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     DO i=1,SIZE(pools)
        CALL pw_pool_give_back_pw(pools(i)%pool,pws(i)%pw,error=error)
     END DO
     DEALLOCATE(pws,stat=stat)
     CPPostconditionNoFail(stat==0,cp_warning_level,routineP,error)
  END IF
END SUBROUTINE pw_pools_give_back_pws
!***************************************************************************

!!****f* pw_pool_types/pw_pools_init_coeffs [1.0] *
!!
!!   NAME
!!     pw_pools_init_coeffs
!!
!!   FUNCTION
!!     initializes a multigrid structure of coeff
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pools: the multigrid pool (i.e. an array of pw_pool)
!!     - coeffs: the multigrid of coefficent you want to initialize
!!     - use_data: which data it uses: REALDATA1D, COMPLEXDATA1D,
!!       REALDATA3D, COMPLEXDATA3D
!!     - mass_flag: whether to allocate mass coeff (mass_cr or mass_cc)
!!       defaults to false
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE pw_pools_init_coeffs(pools, coeffs, use_data,&
     mass_flag, in_space, error)
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pools
    TYPE(coeff_type), DIMENSION(:), &
      INTENT(inout)                          :: coeffs
    INTEGER, INTENT(in)                      :: use_data
    LOGICAL, INTENT(in), OPTIONAL            :: mass_flag
    INTEGER, INTENT(in), OPTIONAL            :: in_space
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_pools_init_coeffs', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(pools),cp_failure_level,routineP,error,failure)
  CPPrecondition(SIZE(pools)==SIZE(coeffs),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     DO i=1,SIZE(pools)
        CALL pw_pool_init_coeff(pools(i)%pool,coeffs(i),use_data,&
             mass_flag=mass_flag,in_space=in_space,error=error)
     END DO
  END IF
END SUBROUTINE pw_pools_init_coeffs
!***************************************************************************

!!****f* pw_pool_types/pw_pools_give_back_coeffs [1.0] *
!!
!!   NAME
!!     pw_pools_give_back_coeffs
!!
!!   FUNCTION
!!     returns the pw part of the coefficents into the pools
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pools: the pools that will cache the coeffs %pw
!!     - coeffs: the coefficents to give back
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE pw_pools_give_back_coeffs(pools, coeffs, error)
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pools
    TYPE(coeff_type), DIMENSION(:), &
      POINTER                                :: coeffs
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_pools_give_back_coeffs', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure

  failure=.FALSE.

  CPPrecondition(ASSOCIATED(pools),cp_failure_level,routineP,error,failure)
  CPPrecondition(SIZE(pools)==SIZE(coeffs),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     DO i=1,SIZE(pools)
        CALL pw_pool_give_back_coeff(pools(i)%pool,coeffs(i),error=error)
     END DO
  END IF
END SUBROUTINE pw_pools_give_back_coeffs
!***************************************************************************

!!****f* pw_pools/pw_pools_copy [1.0] *
!!
!!   NAME
!!     pw_pools_copy
!!
!!   FUNCTION
!!     copies a multigrid pool, the underlying pools are shared
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - source_pool: the pool to copy
!!     - target_pool: will hold the copy of the pool
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE pw_pools_copy(source_pools, target_pools, error)
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: source_pools, target_pools
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_pools_copy', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(source_pools),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     ALLOCATE(target_pools(SIZE(source_pools)), stat=stat)
     CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
  END IF
  IF (.NOT.failure) THEN
     DO i=1,SIZE(source_pools)
        target_pools(i)%pool => source_pools(i)%pool
        CALL pw_pool_retain(source_pools(i)%pool, error=error)
     END DO
  END IF
END SUBROUTINE pw_pools_copy
!***************************************************************************

!!****f* pw_pools/pw_pools_dealloc [1.0] *
!!
!!   NAME
!!     pw_pools_dealloc
!!
!!   FUNCTION
!!     deallocates the given pools (releasing each of the underlying
!!     pools)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pools: the pols to deallocate
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE pw_pools_dealloc(pools,error)
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pools
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_pools_dealloc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  IF (ASSOCIATED(pools)) THEN
     DO i=1,SIZE(pools)
        CALL pw_pool_release(pools(i)%pool, error=error)
     END DO
     DEALLOCATE(pools,stat=stat)
     CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
  END IF
  NULLIFY(pools)
END SUBROUTINE pw_pools_dealloc
!***************************************************************************

!!****f* pw_pool_types/pw_pools_flush_cache [1.0] *
!!
!!   NAME
!!     pw_pools_flush_cache
!!
!!   FUNCTION
!!     deallocates all the cached grids
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pools: the pools to flush
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     08.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE pw_pools_flush_cache(pools,error)
    TYPE(pw_pool_p_type), DIMENSION(:), &
      POINTER                                :: pools
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'pw_pools_flush_cache', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i
    LOGICAL                                  :: failure

  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(pools),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     DO i=1,SIZE(pools)
        CALL pw_pool_flush_cache(pools(i)%pool,error=error)
     END DO
  END IF
END SUBROUTINE pw_pools_flush_cache
!***************************************************************************

!!****f* pw_pool_types/pw_pool_write [1.0] *
!!
!!   NAME
!!     pw_pool_write
!!
!!   FUNCTION
!!     writes the actual contents of the pw_pool
!!
!!   NOTES
!!     -
!!
!!   ARGUMENTS
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     fawzi
!!
!!   MODIFICATION HISTORY
!!     02.2003 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE pw_pool_write(pw_pool,unit_nr,error)
    TYPE(pw_pool_type), POINTER              :: pw_pool
    INTEGER                                  :: unit_nr
    TYPE(cp_error_type), INTENT(inout), &
      OPTIONAL                               :: error

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
     DO WHILE(cp_sll_pw_next(iter,el_att=pw_att,error=error))
        WRITE (unit=unit_nr, fmt="(i8,',')",advance="no") pw_att%id_nr
     END DO
     WRITE (unit=unit_nr, fmt="('),')")

     WRITE (unit=unit_nr, fmt="('real3d_pw=(')",advance="no")
     iter => pw_pool%real3d_pw
     DO WHILE(cp_sll_pw_next(iter,el_att=pw_att,error=error))
        WRITE (unit=unit_nr, fmt="(i8,',')",advance="no") pw_att%id_nr
     END DO
     WRITE (unit=unit_nr, fmt="('),')")

     WRITE (unit=unit_nr, fmt="('complex1d_pw=(')",advance="no")
     iter => pw_pool%complex1d_pw
     DO WHILE(cp_sll_pw_next(iter,el_att=pw_att,error=error))
        WRITE (unit=unit_nr, fmt="(i8,',')",advance="no") pw_att%id_nr
     END DO
     WRITE (unit=unit_nr, fmt="('),')")

     WRITE (unit=unit_nr, fmt="('complex3d_pw=(')",advance="no")
     iter => pw_pool%complex3d_pw
     DO WHILE(cp_sll_pw_next(iter,el_att=pw_att,error=error))
        WRITE (unit=unit_nr, fmt="(i8,',')",advance="no") pw_att%id_nr
     END DO
     WRITE (unit=unit_nr, fmt="(')')")
     
     WRITE (unit=unit_nr, fmt="('}')")
  ELSE
     WRITE (unit=unit_nr, fmt="('<pw_pool *null*>')")
  END IF
END SUBROUTINE pw_pool_write
!***************************************************************************

END MODULE pw_pool_types
