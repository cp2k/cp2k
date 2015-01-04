!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   DBCSR Memory Pool to avoid slow allocations of accelerator memory
!> \author  Ole Schuett
! *****************************************************************************
MODULE dbcsr_mem_methods

  USE acc_stream,                      ONLY: acc_stream_associated,&
                                             acc_stream_equal,&
                                             acc_stream_type
  USE dbcsr_data_methods_low,          ONLY: dbcsr_data_exists,&
                                             dbcsr_data_get_size,&
                                             internal_data_deallocate
  USE dbcsr_data_types,                ONLY: dbcsr_data_obj,&
                                             dbcsr_mempool_entry_type,&
                                             dbcsr_mempool_type,&
                                             dbcsr_memtype_type
  USE dbcsr_error_handling,            ONLY: dbcsr_error_set,&
                                             dbcsr_error_stop,&
                                             dbcsr_error_type
  USE kinds,                           ONLY: dp

  !$ USE OMP_LIB

  IMPLICIT NONE


  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_mem_methods'

  PUBLIC :: dbcsr_mempool_get, dbcsr_mempool_add, dbcsr_mempool_ensure_capacity
  PUBLIC :: dbcsr_mempool_destruct, dbcsr_mempool_clear
  PUBLIC :: dbcsr_memtype_setup, dbcsr_memtype_equal

  LOGICAL, PARAMETER :: debug = .FALSE.

  CONTAINS

! *****************************************************************************
!> \brief Creates a memory pool.
!> \param pool ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mempool_create(pool)
    TYPE(dbcsr_mempool_type), POINTER        :: pool

    IF(ASSOCIATED(pool)) STOP "dbcsr_mempool_create: pool already allocated"
    ALLOCATE(pool)
    !$ CALL OMP_INIT_LOCK(pool%lock)
    ALLOCATE(pool%root) !root allways allocated, but unsused. Simplifies looping.
  END SUBROUTINE dbcsr_mempool_create


! *****************************************************************************
!> \brief Ensures that mempool has at least the given capacity.
!> \param pool ...
!> \param capacity ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mempool_ensure_capacity(pool, capacity)
    TYPE(dbcsr_mempool_type), POINTER        :: pool
    INTEGER, INTENT(IN)                      :: capacity

    IF(.NOT.ASSOCIATED(pool)) STOP "dbcsr_mempool_ensure_capacity: pool not allocated"
    !$ CALL OMP_SET_LOCK(pool%lock)
    pool%capacity = MAX(pool%capacity, capacity)
    !$ CALL OMP_UNSET_LOCK(pool%lock)

  END SUBROUTINE dbcsr_mempool_ensure_capacity

! *****************************************************************************
!> \brief Picks a suitable data_area from mempool, returns Null() if none found.
!> \param memtype ...
!> \param datatype ...
!> \param datasize ...
!> \param error ...
!> \retval res ...
!> \author Ole Schuett
! *****************************************************************************
  FUNCTION dbcsr_mempool_get(memtype, datatype, datasize, error) RESULT(res)
    TYPE(dbcsr_memtype_type)                 :: memtype
    INTEGER, INTENT(IN)                      :: datatype, datasize
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error
    TYPE(dbcsr_data_obj)                     :: res

    INTEGER                                  :: best_size, s
    TYPE(dbcsr_mempool_entry_type), POINTER  :: best_cur, best_prev, cur, prev
    TYPE(dbcsr_mempool_type), POINTER        :: pool

    pool => memtype%pool
    IF(.NOT.ASSOCIATED(pool)) STOP "dbcsr_mempool_get: pool not allocated"

    !$ CALL OMP_SET_LOCK(pool%lock)
    res%d => Null()
    best_cur => Null()
    best_prev => Null()
    best_size = HUGE(1)
    prev => Null()
    cur => pool%root
    DO WHILE(ASSOCIATED(cur%next))
       prev => cur
       cur => cur%next
       s = dbcsr_data_get_size(cur%area)
       IF(s < datasize) CYCLE
       IF(.NOT. dbcsr_memtype_equal(cur%area%d%memory_type, memtype)) CYCLE
       IF(cur%area%d%data_type /= datatype) CYCLE
       !we found a match
       IF(s < best_size) THEN
          best_cur  => cur
          best_prev => prev
          best_size = s
       ENDIF
    ENDDO

    IF(ASSOCIATED(best_cur)) THEN
       IF(best_cur%area%d%refcount /= 0) STOP "mempool_get: refcount /= 0"
       best_cur%area%d%refcount =  1
       best_prev%next => best_cur%next
       res = best_cur%area
       DEALLOCATE(best_cur)
    END IF
    !$ CALL OMP_UNSET_LOCK(pool%lock)

    IF(.NOT.ASSOCIATED(res%d)) &
       CALL mempool_collect_garbage(pool, error)
  END FUNCTION dbcsr_mempool_get


! *****************************************************************************
!> \brief Adds an unused (refcount==0) data_area to the pool.
!> \param area ...
!> \param error ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE dbcsr_mempool_add(area, error)
    TYPE(dbcsr_data_obj)                     :: area
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    TYPE(dbcsr_mempool_entry_type), POINTER  :: new_entry
    TYPE(dbcsr_mempool_type), POINTER        :: pool

    pool => area%d%memory_type%pool
    IF(.NOT.ASSOCIATED(pool)) STOP "dbcsr_mempool_add: pool not allocated"
    IF(.NOT.dbcsr_data_exists(area, error)) STOP "dbcsr_mempool_add: area not allocated"
    IF(area%d%refcount /= 0) STOP "mempool_add: refcount /= 0"

    CALL mempool_collect_garbage(pool, error)

    !$ CALL OMP_SET_LOCK(pool%lock)
    ALLOCATE(new_entry)
    new_entry%area = area
    new_entry%next => pool%root%next
    pool%root%next => new_entry
    !$ CALL OMP_UNSET_LOCK(pool%lock)
 END SUBROUTINE dbcsr_mempool_add

! *****************************************************************************
!> \brief Ensures that pool_size < max_size, e.g. that there is a free slot.
!> \param pool ...
!> \param error ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE mempool_collect_garbage(pool, error)
    TYPE(dbcsr_mempool_type), POINTER        :: pool
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    INTEGER                                  :: n
    TYPE(dbcsr_mempool_entry_type), POINTER  :: cur, prev

    IF(.NOT.ASSOCIATED(pool)) STOP "mempool_collect_garbage: pool not allocated"

    !$ CALL OMP_SET_LOCK(pool%lock)
    prev => pool%root
    cur => pool%root%next
    n = 0
    DO WHILE(ASSOCIATED(cur))
       n = n + 1  
       IF(n >= pool%capacity) THEN
         CALL internal_data_deallocate(cur%area%d, error)
         DEALLOCATE(cur%area%d)
         prev%next => cur%next
         DEALLOCATE(cur)
         cur => prev%next
       ELSE
         prev => cur
         cur => cur%next
       ENDIF
    ENDDO
    !$ CALL OMP_UNSET_LOCK(pool%lock)
 END SUBROUTINE mempool_collect_garbage


! *****************************************************************************
!> \brief Finalizes mempool, includes deallocation of all contained data_areas.
!> \param pool ...
!> \param error ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE dbcsr_mempool_destruct(pool, error)
    TYPE(dbcsr_mempool_type), POINTER        :: pool
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    IF(.NOT.ASSOCIATED(pool)) STOP "dbcsr_mempool_destruct: pool not allocated"

    CALL dbcsr_mempool_clear(pool, error)

    !$ CALL OMP_DESTROY_LOCK (pool%lock)
    DEALLOCATE(pool%root)
    DEALLOCATE(pool)
    NULLIFY(pool)

 END SUBROUTINE dbcsr_mempool_destruct

! *****************************************************************************
!> \brief Deallocates all data_areas contained in given mempool.
!> \param pool ...
!> \param error ...
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE dbcsr_mempool_clear(pool, error)
    TYPE(dbcsr_mempool_type), POINTER        :: pool
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_mempool_clear'

    INTEGER                                  :: error_handler
    TYPE(dbcsr_mempool_entry_type), POINTER  :: cur, prev

    IF(.NOT.ASSOCIATED(pool)) STOP "dbcsr_mempool_clear: pool not allocated"

    CALL dbcsr_error_set (routineN, error_handler, error)

    !$ CALL OMP_SET_LOCK(pool%lock)
    cur => pool%root%next
    DO WHILE(ASSOCIATED(cur))
       CALL internal_data_deallocate(cur%area%d, error)
       DEALLOCATE(cur%area%d)
       prev => cur
       cur => cur%next
       DEALLOCATE(prev)
    ENDDO
    NULLIFY(pool%root%next)
    !$ CALL OMP_UNSET_LOCK(pool%lock)

    CALL dbcsr_error_stop(error_handler, error)
 END SUBROUTINE dbcsr_mempool_clear

! *****************************************************************************
!> \brief Ensures that given memtype has requested settings.
!> \param memtype ...
!> \param acc_hostalloc ...
!> \param acc_devalloc ...
!> \param mpi ...
!> \param acc_stream ...
!> \param oversize_factor ...
!> \param has_pool ...
!> \param error ...
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_memtype_setup(memtype, acc_hostalloc, acc_devalloc, mpi,&
              acc_stream, oversize_factor, has_pool, error)
    TYPE(dbcsr_memtype_type), INTENT(INOUT)  :: memtype
    LOGICAL, INTENT(IN), OPTIONAL            :: acc_hostalloc, acc_devalloc, &
                                                mpi
    TYPE(acc_stream_type), OPTIONAL          :: acc_stream
    REAL(KIND=dp), OPTIONAL                  :: oversize_factor
    LOGICAL, INTENT(IN), OPTIONAL            :: has_pool
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    LOGICAL                                  :: is_ok, my_has_pool
    TYPE(dbcsr_memtype_type)                 :: aim

! variable aim is initialized with default values from type definition

    my_has_pool = .FALSE.
    IF(PRESENT(has_pool)) my_has_pool = has_pool
    IF(PRESENT(acc_hostalloc))   aim%acc_hostalloc   = acc_hostalloc
    IF(PRESENT(acc_devalloc))    aim%acc_devalloc    = acc_devalloc
    IF(PRESENT(mpi))             aim%mpi             = mpi
    IF(PRESENT(acc_stream))      aim%acc_stream      = acc_stream
    IF(PRESENT(oversize_factor)) aim%oversize_factor = oversize_factor

    IF(.NOT. aim%acc_devalloc.EQV.acc_stream_associated(aim%acc_stream))&
       STOP "dbcsr_memtype_setup: acc_stream missing"

    is_ok = .TRUE.
    is_ok = is_ok .AND. (memtype%acc_hostalloc.EQV.aim%acc_hostalloc)
    is_ok = is_ok .AND. (memtype%acc_devalloc.EQV.aim%acc_devalloc)
    is_ok = is_ok .AND. (memtype%mpi.EQV.aim%mpi)
    is_ok = is_ok .AND. acc_stream_equal(memtype%acc_stream,aim%acc_stream)
    is_ok = is_ok .AND. (memtype%oversize_factor==aim%oversize_factor)
    is_ok = is_ok .AND. (ASSOCIATED(memtype%pool).EQV.my_has_pool)

    IF(.NOT.is_ok) THEN
        IF(ASSOCIATED(memtype%pool)) &
           CALL dbcsr_mempool_destruct(memtype%pool, error)

        memtype%acc_hostalloc = aim%acc_hostalloc
        memtype%acc_devalloc  = aim%acc_devalloc
        memtype%mpi = aim%mpi
        memtype%acc_stream = aim%acc_stream
        memtype%oversize_factor = aim%oversize_factor
        IF(my_has_pool) &
           CALL dbcsr_mempool_create(memtype%pool)
    ENDIF
  END SUBROUTINE dbcsr_memtype_setup


! *****************************************************************************
!> \brief Test if two memtypes are equal
!> \param mt1 ...
!> \param mt2 ...
!> \retval res ...
!> \author Ole Schuett
! *****************************************************************************
  FUNCTION dbcsr_memtype_equal(mt1,mt2) RESULT (res)
    TYPE(dbcsr_memtype_type), INTENT(in)     :: mt1, mt2
    LOGICAL                                  :: res

      res = (mt1%mpi .EQV. mt2%mpi) .AND. &
            (mt1%acc_hostalloc .EQV. mt2%acc_hostalloc) .AND. &
            (mt1%acc_devalloc .EQV. mt2%acc_devalloc) .AND. &
            (ASSOCIATED(mt1%pool) .EQV. ASSOCIATED(mt2%pool)) .AND. &
            (.NOT. ASSOCIATED(mt1%pool) .OR. ASSOCIATED(mt1%pool, mt2%pool))
  END FUNCTION dbcsr_memtype_equal

END MODULE dbcsr_mem_methods
