!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/pw_pool [1.0] *
!!
!!   NAME
!!     pw_pool
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
MODULE pw_pool
  USE cp_log_handling, ONLY: cp_fatal_level, cp_failure_level,&
       cp_warning_level, cp_note_level, cp_to_string, cp_log
  USE cp_error_handling, ONLY: cp_debug, cp_error_type, cp_error_init,&
       cp_error_dealloc_ref, cp_error_message, cp_assert,&
       cp_assertion_failed, cp_internal_error, cp_unimplemented_error,&
       cp_error_get_logger
  USE kinds, ONLY: wp=>dp
  USE cp_para_types, ONLY: cp_para_env_type
  USE timings, ONLY: timeset, timestop
  USE cp_linked_list_pw, ONLY: cp_sll_pw_type, cp_sll_pw_next,&
       cp_sll_pw_dealloc, cp_sll_pw_get_first_el, cp_sll_pw_rm_first_el,&
       cp_sll_pw_insert_el
  USE pw_grid_types, ONLY: pw_grid_type
  USE pw_grids, ONLY: pw_grid_destruct
  USE pw_types, ONLY: pw_type, REALDATA1D, REALDATA3D, COMPLEXDATA1D,&
       COMPLEXDATA3D, pw_allocate, pw_deallocate
  USE coefficient_types, ONLY: coeff_type
  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PRIVATE, PARAMETER :: moduleN='pw_pool'
  INTEGER, SAVE, PRIVATE :: last_pw_pool_id_nr=0

!***
!****************************************************************************

!!****s* pw_pool/pw_pool_type [1.0] *
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
!!     -
!!
!!   ATTRIBUTES
!!     - to do
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
     INTEGER :: ref_count, id_nr
     TYPE(pw_grid_type), POINTER :: pw_grid
     TYPE(cp_sll_pw_type), POINTER :: real1d_pw, real3d_pw,&
          complex1d_pw, complex3d_pw
  END TYPE pw_pool_type
!!***
!****************************************************************************

  !!****s* pw_pool/pw_pool_p_type [1.0] *
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

!!****f* pw_pool/pw_pool_create [1.0] *
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
  SUBROUTINE pw_pool_create(pool, pw_grid_ptr, error)
    TYPE(pw_pool_type), POINTER :: pool
    TYPE(pw_grid_type), POINTER :: pw_grid_ptr
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='pw_pool_create',&
         routineP=moduleN//':'//routineN
    INTEGER :: stat
    failure=.FALSE.

    ALLOCATE(pool, stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    IF (.NOT.failure) THEN
       pool%pw_grid => pw_grid_ptr
       last_pw_pool_id_nr=last_pw_pool_id_nr+1
       pool%id_nr=last_pw_pool_id_nr
       pool%ref_count=1
       NULLIFY(pool%real1d_pw, pool%real3d_pw, &
            pool%complex1d_pw, pool%complex3d_pw)
    END IF
  END SUBROUTINE pw_pool_create
!***************************************************************************

!!****f* pw_pool/pw_pool_retain [1.0] *
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
    TYPE(pw_pool_type), POINTER :: pool
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='pw_pool_retain',&
         routineP=moduleN//':'//routineN
    failure=.FALSE.

    CPPrecondition(ASSOCIATED(pool),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPrecondition(pool%ref_count>0,cp_failure_level,routineP,error,failure)

       pool%ref_count=pool%ref_count+1
    END IF
  END SUBROUTINE pw_pool_retain
!***************************************************************************

!!****f* pw_pool/pw_pool_release [1.0] *
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
    TYPE(pw_pool_type), POINTER :: pool
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='pw_pool_release',&
         routineP=moduleN//':'//routineN
    TYPE(cp_sll_pw_type), POINTER :: iterator
    TYPE(pw_type), POINTER :: pw_el
    INTEGER :: stat
    failure=.FALSE.

    IF (ASSOCIATED(pool)) THEN
       CPPrecondition(pool%ref_count>0,cp_failure_level,routineP,error,failure)
       pool%ref_count=pool%ref_count-1
       IF (pool%ref_count<1) THEN

          iterator=pool%real1d_pw
          DO
             IF (.NOT.cp_sll_pw_next(iterator,el_att=pw_el,error=error)) EXIT
             CALL pw_deallocate(pw_el)
          END DO
          CALL cp_sll_pw_dealloc(pool%real1d_pw,error=error)

          iterator=pool%real3d_pw
          DO
             IF (.NOT.cp_sll_pw_next(iterator,el_att=pw_el,error=error)) EXIT
             CALL pw_deallocate(pw_el)
          END DO
          CALL cp_sll_pw_dealloc(pool%real3d_pw,error=error)

          iterator=pool%complex1d_pw
          DO
             IF (.NOT.cp_sll_pw_next(iterator,el_att=pw_el,error=error)) EXIT
             CALL pw_deallocate(pw_el)
          END DO
          CALL cp_sll_pw_dealloc(pool%complex1d_pw,error=error)

          iterator=pool%complex3d_pw
          DO
             IF (.NOT.cp_sll_pw_next(iterator,el_att=pw_el,error=error)) EXIT
             CALL pw_deallocate(pw_el)
          END DO
          CALL cp_sll_pw_dealloc(pool%complex3d_pw)

          CPPrecondition(ASSOCIATED(pool%pw_grid),cp_warning_level,routineP,error,failure)
          IF (.NOT.failure) THEN
             CALL pw_grid_destruct(pool%pw_grid)
          END IF

          DEALLOCATE(pool,stat=stat)
          CPPostcondition(stat==0,cp_warning_level,routineP,error,failure)
       END IF
    END IF
    NULLIFY(pool)
  END SUBROUTINE pw_pool_release
!***************************************************************************

!!****f* pw_pool/try_pop [1.0] *
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
    TYPE(pw_type), POINTER ::res
    TYPE(cp_sll_pw_type), POINTER :: list

    IF (ASSOCIATED(list)) THEN
       res => cp_sll_pw_get_first_el(list)
       CALL cp_sll_pw_rm_first_el(list)
    ELSE
       NULLIFY(res)
    END IF
  END FUNCTION try_pop
!***************************************************************************

!!****f* pw_pool/pw_pool_get_pw [1.0] *
!!
!!   NAME
!!     pw_pool_get_pw
!!
!!   FUNCTION
!!     returns a pw, allocating it if none is in the pool
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - pool: the ppol from where you get the pw
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
  FUNCTION pw_pool_get_pw(pool, use_data, mass_flag, error) RESULT(res)
    TYPE(pw_pool_type), POINTER :: pool
    INTEGER, INTENT(in) :: use_data
    LOGICAL, INTENT(in), OPTIONAL :: mass_flag
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
    TYPE(pw_type), POINTER :: res

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='pw_pool_get_pw',&
         routineP=moduleN//':'//routineN
    INTEGER :: stat
    failure=.FALSE.

    CPPrecondition(ASSOCIATED(pool),cp_failure_level,routineP,error,failure)
    CPPrecondition(pool%ref_count>0,cp_failure_level,routineP,error,failure)

    SELECT CASE(use_data)
    CASE (REALDATA1D)
       res => try_pop(pool%real1d_pw)
    CASE (REALDATA3D)
       res => try_pop(pool%real3d_pw)
    CASE (COMPLEXDATA1D)
       res => try_pop(pool%complex1d_pw)
    CASE (COMPLEXDATA3D)
       res => try_pop(pool%complex3d_pw)   
    CASE default
! unknown use_data
       CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
    END SELECT

    IF (.NOT.ASSOCIATED(res)) THEN
       ALLOCATE(res,stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       CALL pw_allocate(res, pool%pw_grid, use_data=use_data, &
            mass_flag=mass_flag)
    END IF

    IF (PRESENT(mass_flag)) THEN
       IF (mass_flag) THEN
          CALL cp_assert(.FALSE.,cp_warning_level,cp_assertion_failed,&
               routineP,"mass not cached in pool"//&
               CPSourceFileRef,&
               error=error)

          IF (use_data==REALDATA1D) THEN
             IF (.NOT.ASSOCIATED(res%mass_cr)) THEN
                ALLOCATE(res%mass_cr(res % pw_grid % ngpts_cut_local),&
                     stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
          ELSE IF (use_data==COMPLEXDATA1D) THEN
             IF (.NOT.ASSOCIATED(res%mass_cc)) THEN
                ALLOCATE(res%mass_cc(res % pw_grid % ngpts_cut_local),&
                     stat=stat)
                CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
             END IF
          END IF
       END IF
    END IF

  END FUNCTION pw_pool_get_pw
!***************************************************************************

!!****f* pw_pool/pw_pool_give_back_pw [1.0] *
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
  SUBROUTINE pw_pool_give_back_pw(pool, pw, error)
    TYPE(pw_pool_type), POINTER :: pool
    TYPE(pw_type), POINTER :: pw

    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='pw_pool_give_back_pw',&
         routineP=moduleN//':'//routineN
    INTEGER :: stat
    failure=.FALSE.

    CPPrecondition(ASSOCIATED(pool),cp_failure_level,routineP,error,failure)
    CPPrecondition(pool%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(pw),cp_warning_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL cp_assert(pw%pw_grid%identifier==pool%pw_grid%identifier,&
            cp_failure_level, cp_assertion_failed, routineP,&
            "pool cannot reuse pw of another grid "//&
            CPSourceFileRef,&
            error=error,failure=failure)
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

       SELECT CASE(pw%in_use)
       CASE(REALDATA1D)
          CALL cp_sll_pw_insert_el(pool%real1d_pw, element=pw)
       CASE(REALDATA3D)
          CALL cp_sll_pw_insert_el(pool%real3d_pw, element=pw)
       CASE(COMPLEXDATA1D)
          CALL cp_sll_pw_insert_el(pool%complex1d_pw, element=pw)
       CASE(COMPLEXDATA3D)
          CALL cp_sll_pw_insert_el(pool%complex3d_pw, element=pw)
       CASE default
! unknown in_use
          CPAssert(.FALSE.,cp_failure_level,routineP,error,failure)
       END SELECT

    END IF
  END SUBROUTINE pw_pool_give_back_pw
!***************************************************************************

  !!****f* pw_pool/pw_pool_init_coeff [1.0] *
!!
!!   NAME
!!     pw_pool_init_coeff
!!
!!   FUNCTION
!!     initializes a coefficent with the corresponding grid type
!!
!!   NOTES
!!     -
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
SUBROUTINE pw_pool_init_coeff(pool, coeff, use_data, mass_flag, error)
    TYPE(pw_pool_type), POINTER :: pool
    TYPE(coeff_type), INTENT(out) :: coeff
    INTEGER, INTENT(in) :: use_data
    LOGICAL, INTENT(in), OPTIONAL :: mass_flag
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='pw_pool_init_coeff',&
        routineP=moduleN//':'//routineN
  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(pool),cp_failure_level,routineP,error,failure)
  CPPrecondition(pool%ref_count>0,cp_failure_level,routineP,error,failure)
  
  IF (.NOT. failure) THEN
     coeff%pw => pw_pool_get_pw(pool,use_data,mass_flag=mass_flag,error=error)
  END IF
END SUBROUTINE pw_pool_init_coeff
!***************************************************************************

!!****f* pw_pool/pw_pool_give_back_coeff [1.0] *
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
SUBROUTINE pw_pool_give_back_coeff(pool,coeff,error)
  TYPE(pw_pool_type), POINTER :: pool
  TYPE(coeff_type), INTENT(in) :: coeff
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='pw_pool_give_back_coeff',&
        routineP=moduleN//':'//routineN
  failure=.FALSE.
  
  CALL pw_pool_give_back_pw(pool,coeff%pw)
END SUBROUTINE pw_pool_give_back_coeff
!***************************************************************************

!!****f* pw_pools_init_coeffs [1.0] *
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
     mass_flag, error)
  TYPE(pw_pool_p_type), DIMENSION(:), POINTER :: pools
  TYPE(coeff_type), DIMENSION(:), INTENT(out)    :: coeffs
  INTEGER, INTENT(in) :: use_data
  LOGICAL, INTENT(in), OPTIONAL :: mass_flag
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='pw_pools_init_coeffs',&
        routineP=moduleN//':'//routineN
  INTEGER :: i
  failure=.FALSE.
  
  CPPrecondition(ASSOCIATED(pools),cp_failure_level,routineP,error,failure)
  CPPrecondition(SIZE(pools)==SIZE(coeffs),cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     DO i=1,SIZE(pools)
        CALL pw_pool_init_coeff(pools(i)%pool,coeffs(i),use_data,&
             mass_flag=mass_flag,error=error)
     END DO
  END IF
END SUBROUTINE pw_pools_init_coeffs
!***************************************************************************

!!****f* pw_pool/pw_pools_give_back_coeffs [1.0] *
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
  TYPE(pw_pool_p_type), DIMENSION(:), POINTER :: pools
  TYPE(coeff_type), DIMENSION(:), INTENT(in) :: coeffs
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='pw_pools_give_back_coeffs',&
        routineP=moduleN//':'//routineN
  INTEGER :: i
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

END MODULE pw_pool
