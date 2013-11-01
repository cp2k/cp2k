!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2013  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   DBCSR Memory Pool to avoid slow allocations e.g. of cuda memory
!> \author  Ole Schuett
! *****************************************************************************
MODULE dbcsr_mem_methods

  USE dbcsr_data_methods_low,          ONLY: dbcsr_data_exists,&
                                             dbcsr_data_get_size,&
                                             internal_data_deallocate
  USE dbcsr_error_handling,            ONLY: dbcsr_error_set,&
                                             dbcsr_error_stop,&
                                             dbcsr_error_type
  USE dbcsr_kinds,                     ONLY: dp
  USE dbcsr_types,                     ONLY: dbcsr_data_obj,&
                                             dbcsr_mempool_entry_type,&
                                             dbcsr_mempool_type,&
                                             dbcsr_memtype_type

  !$ USE OMP_LIB

  IMPLICIT NONE


  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_mem_methods'

  PUBLIC :: dbcsr_mempool_get, dbcsr_mempool_add, dbcsr_mempool_ensure_capacity
  PUBLIC :: dbcsr_mempool_create, dbcsr_mempool_destruct
  PUBLIC :: dbcsr_memtype_setup, dbcsr_memtype_equal

  LOGICAL, PARAMETER :: debug = .FALSE.


  CONTAINS

! *****************************************************************************
! \brief Creates a memory pool.
! \param capacity number of data_areas to hold at most within the pool
! \param name for debug purpose
! \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mempool_create(pool, capacity, name)
    TYPE(dbcsr_mempool_type), POINTER        :: pool
    INTEGER, INTENT(IN)                      :: capacity
    CHARACTER(len=*)                         :: name

    IF(ASSOCIATED(pool)) STOP "dbcsr_mempool_create: pool already allocated"
    ALLOCATE(pool)
    !$ CALL OMP_INIT_LOCK(pool%lock)
    pool%capacity = capacity
    pool%name = name
    ALLOCATE(pool%root) !root allways allocated, but unsused. Simplifies looping.
    NULLIFY(pool%root%next)

    IF(debug) WRITE (*,*) "mempool_create: name=",TRIM(name), " capacity=",capacity
  END SUBROUTINE dbcsr_mempool_create


! *****************************************************************************
! \brief Ensures that mempool has at least the given capacity.
! \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_mempool_ensure_capacity(pool, capacity)
    TYPE(dbcsr_mempool_type), POINTER        :: pool
    INTEGER, INTENT(IN)                      :: capacity

    IF(.NOT.ASSOCIATED(pool)) STOP "dbcsr_mempool_ensure_capacity: pool not allocated"
    !$ CALL OMP_SET_LOCK(pool%lock)
    pool%capacity = MAX(pool%capacity, capacity)
    !$ CALL OMP_UNSET_LOCK(pool%lock)

    IF(debug) WRITE (*,*) "mempool_ensure_capacity: name=",TRIM(pool%name), " capacity=",capacity

  END SUBROUTINE dbcsr_mempool_ensure_capacity

! *****************************************************************************
! \brief Picks a suitable data_area from mempool, returns Null() if none found.
! \author Ole Schuett
! *****************************************************************************
  FUNCTION dbcsr_mempool_get(memtype, datatype, datasize, error) RESULT(res)
    TYPE(dbcsr_memtype_type)                 :: memtype
    INTEGER, INTENT(IN)                      :: datatype, datasize
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error
    TYPE(dbcsr_data_obj)                     :: res

    TYPE(dbcsr_mempool_entry_type), POINTER  :: cur, prev
    TYPE(dbcsr_mempool_type), POINTER        :: pool

    pool => memtype%pool
    IF(.NOT.ASSOCIATED(pool)) STOP "dbcsr_mempool_get: pool not allocated"

    res%d => Null()

    !$ CALL OMP_SET_LOCK(pool%lock)
    prev => Null()
    cur => pool%root
    DO WHILE(ASSOCIATED(cur%next))
       prev => cur
       cur => cur%next
       IF(dbcsr_data_get_size(cur%area) < datasize) CYCLE
       IF(.NOT. dbcsr_memtype_equal(cur%area%d%memory_type, memtype)) CYCLE
       IF(cur%area%d%data_type /= datatype) CYCLE
       !we found a match
       IF(debug) WRITE (*,*) "mempool_get: name=",TRIM(pool%name)," found area_id=",cur%area%d%id
       IF(cur%area%d%refcount /= 0) STOP "mempool_get: refcount /= 0"
       cur%area%d%refcount =  1
       prev%next => cur%next
       res = cur%area
       DEALLOCATE(cur)
       EXIT
    ENDDO
    !$ CALL OMP_UNSET_LOCK(pool%lock)

    IF(.NOT.ASSOCIATED(res%d)) &
       CALL mempool_collect_garbage(pool, error)

    !IF(.NOT.ASSOCIATED(res%d)) WRITE (*,*) "mempool_get missed"
  END FUNCTION dbcsr_mempool_get


! *****************************************************************************
! \brief Adds an unused (refcount==0) data_area to the pool.
! \author Ole Schuett
! *****************************************************************************
 SUBROUTINE dbcsr_mempool_add(area, error)
    TYPE(dbcsr_data_obj)                     :: area
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    TYPE(dbcsr_mempool_entry_type), POINTER  :: new_entry
    TYPE(dbcsr_mempool_type), POINTER        :: pool

    pool => area%d%memory_type%pool
    IF(.NOT.ASSOCIATED(pool)) STOP "dbcsr_mempool_add: pool not allocated"
    IF(.NOT.dbcsr_data_exists(area, error)) STOP "dbcsr_mempool_add: area not allocated"

    IF(debug) WRITE (*,*) "mempool_add: name=",TRIM(pool%name),&
          " area_id=",area%d%id, " size=",dbcsr_data_get_size(area)

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
! \brief Ensures that pool_size < max_size, e.g. that there is a free slot.
! \author Ole Schuett
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
         IF(debug) WRITE (*,*) "mempool_collect_garbage: deallocating, ",&
           "name=", TRIM(pool%name), " size: ", dbcsr_data_get_size(cur%area)
         CALL internal_data_deallocate(cur%area%d, error)
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
! \brief Finalizes mempool, includes deallocation of all contained data_areas.
! \author Ole Schuett
! *****************************************************************************
 SUBROUTINE dbcsr_mempool_destruct(pool, error)
    TYPE(dbcsr_mempool_type), POINTER        :: pool
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_mempool_destruct'

    INTEGER                                  :: error_handler
    TYPE(dbcsr_mempool_entry_type), POINTER  :: cur, prev

    IF(.NOT.ASSOCIATED(pool)) STOP "dbcsr_mempool_destruct: pool not allocated"

    CALL dbcsr_error_set (routineN, error_handler, error)

    IF(debug) WRITE(*,*) "mempool_destruct name=", TRIM(pool%name)

    !$ CALL OMP_SET_LOCK(pool%lock)
    cur => pool%root%next
    DO WHILE(ASSOCIATED(cur))
       CALL internal_data_deallocate(cur%area%d, error)
       prev => cur
       cur => cur%next
       DEALLOCATE(prev)
    ENDDO
    !$ CALL OMP_UNSET_LOCK(pool%lock)

    !$ CALL OMP_DESTROY_LOCK (pool%lock)
    DEALLOCATE(pool%root)
    DEALLOCATE(pool)
    NULLIFY(pool)

    CALL dbcsr_error_stop(error_handler, error)
 END SUBROUTINE dbcsr_mempool_destruct


! *****************************************************************************
! \brief Ensures that given memtype has requested settings.
! \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_memtype_setup(memtype, cuda_hostalloc, cuda_devalloc, mpi, &
             oversize_factor, pool_capacity, pool_name, error)
    TYPE(dbcsr_memtype_type), INTENT(INOUT)  :: memtype
    LOGICAL, INTENT(IN), OPTIONAL            :: cuda_hostalloc, &
                                                cuda_devalloc, mpi
    REAL(KIND=dp), OPTIONAL                  :: oversize_factor
    INTEGER, INTENT(IN), OPTIONAL            :: pool_capacity
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: pool_name
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    INTEGER                                  :: capacity
    LOGICAL                                  :: is_ok
    TYPE(dbcsr_memtype_type)                 :: aim

! variable aim is initialized with default values from type definition

    capacity = 0
    IF(PRESENT(pool_capacity)) THEN
       capacity = pool_capacity
       IF(.NOT. PRESENT(pool_name)) STOP "dbcsr_memtype_setup: pool_name missing"
    ENDIF
    IF(PRESENT(cuda_hostalloc))  aim%cuda_hostalloc  = cuda_hostalloc
    IF(PRESENT(cuda_devalloc))   aim%cuda_devalloc   = cuda_devalloc
    IF(PRESENT(mpi))             aim%mpi             = mpi
    IF(PRESENT(oversize_factor)) aim%oversize_factor = oversize_factor

    is_ok = .TRUE.
    is_ok = is_ok .AND. (memtype%cuda_hostalloc.EQV.aim%cuda_hostalloc)
    is_ok = is_ok .AND. (memtype%cuda_devalloc.EQV.aim%cuda_devalloc)
    is_ok = is_ok .AND. (memtype%mpi.EQV.aim%mpi)
    is_ok = is_ok .AND. (memtype%oversize_factor==aim%oversize_factor)

    is_ok = is_ok .AND. (ASSOCIATED(memtype%pool).EQV.(capacity>0))
    IF(ASSOCIATED(memtype%pool)) THEN
      is_ok = is_ok .AND. (memtype%pool%capacity==capacity)
      is_ok = is_ok .AND. (TRIM(memtype%pool%name)==TRIM(pool_name))
    ENDIF

    IF(.NOT.is_ok) THEN
        IF(ASSOCIATED(memtype%pool)) &
           CALL dbcsr_mempool_destruct(memtype%pool, error)

        memtype%cuda_hostalloc = aim%cuda_hostalloc
        memtype%cuda_devalloc = aim%cuda_devalloc
        memtype%mpi = aim%mpi
        memtype%oversize_factor = aim%oversize_factor
        IF(capacity>0) &
           CALL dbcsr_mempool_create(memtype%pool, capacity=capacity, name=pool_name)
    ENDIF
  END SUBROUTINE dbcsr_memtype_setup


! *****************************************************************************
! \brief Test if two memtypes are equal
! \author Ole Schuett
! *****************************************************************************
  FUNCTION dbcsr_memtype_equal(mt1,mt2) RESULT (res)
    TYPE(dbcsr_memtype_type), INTENT(in)     :: mt1, mt2
    LOGICAL                                  :: res

      res = (mt1%mpi .EQV. mt2%mpi) .AND. &
            (mt1%cuda_hostalloc .EQV. mt2%cuda_hostalloc) .AND. &
            (mt1%cuda_devalloc .EQV. mt2%cuda_devalloc) .AND. &
            (ASSOCIATED(mt1%pool) .EQV. ASSOCIATED(mt2%pool)) .AND. &
            (.NOT. ASSOCIATED(mt1%pool) .OR. ASSOCIATED(mt1%pool, mt2%pool))
  END FUNCTION dbcsr_memtype_equal

END MODULE dbcsr_mem_methods
