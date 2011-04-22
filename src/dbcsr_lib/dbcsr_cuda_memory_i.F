!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2011  Urban Borstnik and the CP2K developers group   !
!-----------------------------------------------------------------------------!

  SUBROUTINE dev_mem_alloc_i (dev_mem, n, error, stat)
    TYPE(dbcsr_cuda_mem_type_i4), &
      INTENT(OUT)                                   :: dev_mem
    INTEGER, INTENT(IN)                             :: n
    TYPE(dbcsr_error_type), INTENT(INOUT), OPTIONAL :: error
    INTEGER, INTENT(OUT), OPTIONAL                  :: stat

    CHARACTER(len=*), PARAMETER :: routineN = 'dev_mem_alloc_i', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle, istat

!   ---------------------------------------------------------------------------

    IF (PRESENT (error)) CALL dbcsr_error_set (routineN, error_handle, error)
    dev_mem%req_len = n
    dev_mem%len = MAX(1, n)
#if defined (__DBCSR_CUDA)
    istat = cuda_dev_mem_alloc_cu(dev_mem%ref,&
         INT(dev_mem%len*int_4_size, KIND=C_SIZE_T))
#else
    istat = -1
#endif
    IF (istat /= 0) THEN
       dev_mem%len = -1
       IF (PRESENT(error)) THEN
          CALL dbcsr_assert (istat, "EQ", 0,&
               dbcsr_fatal_level, dbcsr_internal_error, routineN,&
               "Could not allocate GPU device memory",&
               __LINE__, error=error)
       ENDIF
    ENDIF
    IF (PRESENT (stat)) stat = istat
    IF (PRESENT (error)) CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE dev_mem_alloc_i


  SUBROUTINE dev_mem_dealloc_i (dev_mem, error, stat)
    TYPE(dbcsr_cuda_mem_type_i4), INTENT(INOUT) :: dev_mem
    TYPE(dbcsr_error_type), INTENT(INOUT), OPTIONAL      :: error
    INTEGER, INTENT(OUT), OPTIONAL                       :: stat
    

    CHARACTER(len=*), PARAMETER :: routineN = 'dev_mem_dealloc_i', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle, istat

!   ---------------------------------------------------------------------------

    IF (PRESENT (error)) CALL dbcsr_error_set (routineN, error_handle, error)
    IF (dev_mem%len .GE. 0) THEN
#if defined (__DBCSR_CUDA)
       istat = cuda_dev_mem_dealloc_cu(dev_mem%ref)
#else
       istat = -1
#endif
    ELSE
       istat = 0
    ENDIF
    IF (istat /= 0) THEN
       IF (PRESENT (error)) THEN
          CALL dbcsr_assert (istat, "EQ", 0,&
               dbcsr_fatal_level, dbcsr_internal_error, routineN,&
               "Could not deallocate GPU device memory",&
               __LINE__, error=error)
       ENDIF
    ENDIF
    IF (PRESENT (stat)) stat = istat
    IF (PRESENT (error)) CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE dev_mem_dealloc_i


  SUBROUTINE dev_mem_realloc_i (dev_mem, n, error, stat, memory_crunch)
    TYPE(dbcsr_cuda_mem_type_i4), &
      INTENT(INOUT)                          :: dev_mem
    INTEGER, INTENT(IN)                      :: n
    TYPE(dbcsr_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error
    INTEGER, INTENT(OUT), OPTIONAL           :: stat
    LOGICAL, INTENT(OUT), OPTIONAL           :: memory_crunch

    CHARACTER(len=*), PARAMETER :: routineN = 'dev_mem_realloc_i', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle, istat
    INTEGER(KIND=C_INT)                      :: mem_crunch
    INTEGER(KIND=C_SIZE_T)                   :: old_size

!   ---------------------------------------------------------------------------

    IF (PRESENT (error)) CALL dbcsr_error_set (routineN, error_handle, error)
    old_size = MAX(1, dev_mem%req_len*int_4_size)
    dev_mem%req_len = n
    dev_mem%len = MAX(1, n)
#if defined (__DBCSR_CUDA)
    istat = cuda_dev_mem_realloc_cu(dev_mem%ref,&
         INT(dev_mem%len*int_4_size, KIND=C_SIZE_T),&
         old_size,&
         mem_crunch)
    IF (PRESENT (memory_crunch)) &
       memory_crunch = mem_crunch .NE. 0
#else
    istat = -1
    IF (PRESENT (memory_crunch)) &
         memory_crunch = .FALSE.
#endif
    IF (istat /= 0) THEN
       dev_mem%len = -1
       IF (PRESENT(error)) THEN
          CALL dbcsr_assert (istat, "EQ", 0,&
               dbcsr_fatal_level, dbcsr_internal_error, routineN,&
               "Could not reallocate GPU device memory",&
               __LINE__, error=error)
       ENDIF
    ENDIF
    IF (PRESENT (stat)) stat = istat
    IF (PRESENT (error)) CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE dev_mem_realloc_i


  SUBROUTINE host_mem_alloc_i (host_mem, n, write_combined, all_contexts, error, stat)
    INTEGER(KIND=int_4), DIMENSION(:), &
      POINTER                                :: host_mem
    INTEGER, INTENT(IN)                      :: n
    LOGICAL, INTENT(IN), OPTIONAL            :: write_combined, all_contexts
    TYPE(dbcsr_error_type), INTENT(INOUT), OPTIONAL :: error
    INTEGER, INTENT(OUT), OPTIONAL           :: stat

    CHARACTER(len=*), PARAMETER :: routineN = 'host_mem_alloc_i', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle, istat
    INTEGER(KIND=c_int)                      :: all_flag, wc_flag
    TYPE(C_PTR)                              :: host_mem_c_ptr

!   ---------------------------------------------------------------------------

    IF (PRESENT (error)) CALL dbcsr_error_set (routineN, error_handle, error)
    IF (PRESENT (write_combined)) THEN
       wc_flag = 1
    ELSE
       wc_flag = 0
    ENDIF
    IF (PRESENT (all_contexts)) THEN
       all_flag = 1
    ELSE
       all_flag = 0
    ENDIF
#if defined (__DBCSR_CUDA)
    istat = cuda_host_mem_alloc_cu(host_mem_c_ptr,&
         INT(MAX(1,n)*int_4_size, KIND=C_SIZE_T),&
         wc_flag, all_flag)
#else
    istat = -1
#endif
    IF (istat /= 0 .AND. PRESENT (error)) THEN
       CALL dbcsr_assert (istat, "EQ", 0,&
            dbcsr_fatal_level, dbcsr_internal_error, routineN,&
            "Could not allocate host pinned memory",&
            __LINE__, error=error)
       NULLIFY (host_mem)
    ELSE
#if defined (__DBCSR_CUDA)
       CALL C_F_POINTER (host_mem_c_ptr, host_mem, (/ MAX(1,n) /))
#endif
    ENDIF
    IF (PRESENT (stat)) stat = istat
    IF (PRESENT (error)) CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE host_mem_alloc_i


  SUBROUTINE host_mem_dealloc_i (host_mem, error, stat)
    INTEGER(KIND=int_4), DIMENSION(:), &
      POINTER                                :: host_mem
    TYPE(dbcsr_error_type), INTENT(INOUT), OPTIONAL :: error
    INTEGER, INTENT(OUT), OPTIONAL           :: stat

    CHARACTER(len=*), PARAMETER :: routineN = 'host_mem_dealloc_i', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle, istat
    INTEGER(KIND=int_4), POINTER                         :: host_first_el
    TYPE(C_PTR)                              :: host_mem_c_ptr

!   ---------------------------------------------------------------------------

    IF (PRESENT (error)) CALL dbcsr_error_set (routineN, error_handle, error)
#if defined (__DBCSR_CUDA)
    IF (SIZE (host_mem) .GT. 0) THEN
       host_first_el => host_mem(1)
       host_mem_c_ptr = C_LOC (host_first_el)
    ENDIF
    IF (SIZE (host_mem) .GT. 0) THEN
       istat = cuda_host_mem_dealloc_cu (host_mem_c_ptr)
    ENDIF
#else
    istat = -1
#endif
    IF (istat /= 0 .AND. PRESENT (error)) THEN
       CALL dbcsr_assert (istat, "EQ", 0,&
            dbcsr_fatal_level, dbcsr_internal_error, routineN,&
            "Error deallocating host pinned memory",&
            __LINE__, error=error)
    ENDIF
    IF (PRESENT (stat)) stat = istat
    IF (PRESENT (error)) CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE host_mem_dealloc_i



  SUBROUTINE dev_mem_zero_i (dev_mem, first, last, error)
    TYPE(dbcsr_cuda_mem_type_i4), &
      INTENT(INOUT)                          :: dev_mem
    INTEGER, INTENT(IN)                      :: first, last
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    INTEGER(KIND=C_SIZE_T)                   :: len, off
    INTEGER                                  :: error_handle, istat
    CHARACTER(len=*), PARAMETER :: routineN = 'dev_mem_zero_i', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handle, error)
#if defined (__DBCSR_CUDA)
    IF (last .GE. first) THEN
       off = (first-1) * int_4_size
       len = (last-first+1) * int_4_size
       istat = cuda_dev_zero_cu(dev_mem%ref, off, len)
    ELSE
       istat = 0;
    ENDIF
#else
    istat = -1
#endif
    IF (istat /= 0) THEN
       CALL dbcsr_assert (istat, "EQ", 0,&
            dbcsr_fatal_level, dbcsr_internal_error, routineN,&
            "Error zeroing device memory.",&
            __LINE__, error=error)
    ENDIF
    CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE dev_mem_zero_i
