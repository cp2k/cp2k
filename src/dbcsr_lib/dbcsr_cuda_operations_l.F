!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2013  Urban Borstnik and the CP2K developers group   !
!-----------------------------------------------------------------------------!

  SUBROUTINE do_mm_stack_l (param_stack, stack_size, nparams,&
       a_data, b_data, c_data, c_locks, m_max, n_max, k_max, def_mnk, stream, error)
    TYPE(dbcsr_cuda_mem_type_i4), INTENT(IN) :: param_stack
    INTEGER, INTENT(IN)                      :: stack_size, nparams
    TYPE(dbcsr_cuda_mem_type_i8), INTENT(IN) :: a_data, b_data
    TYPE(dbcsr_cuda_mem_type_i8), &
      INTENT(INOUT)                          :: c_data
    TYPE(dbcsr_cuda_mem_type_i4), INTENT(IN) :: c_locks
    INTEGER, INTENT(IN)                      :: m_max, n_max, k_max
    LOGICAL, INTENT(IN)                      :: def_mnk
    INTEGER, INTENT(IN), OPTIONAL            :: stream
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'do_mm_stack_l', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle, istat
    INTEGER(KIND=C_INT)                      :: mnk, stream_l

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handle, error)
    IF (PRESENT (stream)) THEN
       stream_l = stream
    ELSE
       stream_l = default_stream
    ENDIF
#if defined (__DBCSR_CUDA)
    IF (def_mnk) THEN
       mnk = 1
    ELSE
       mnk = 0
    ENDIF
    istat = cuda_do_stack_cu(param_stack%ref, INT(stack_size, KIND=C_INT),&
         INT(nparams, KIND=C_INT),&
         INT(dbcsr_type_int_8, KIND=C_INT),&
         a_data%ref, b_data%ref, c_data%ref, c_locks%ref,&
         INT(m_max,KIND=C_INT), INT(n_max,KIND=C_INT), INT(k_max,KIND=C_INT),&
         mnk, stream_l)
#else
    istat = -1
#endif
    IF (istat /= 0) THEN
       CALL dbcsr_assert (istat, "EQ", 0,&
            dbcsr_fatal_level, dbcsr_internal_error, routineN,&
            "Error executing CUDA kernel.",&
            __LINE__, error=error)
    ENDIF
    CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE do_mm_stack_l




  SUBROUTINE cp_host_to_dev_l (host_mem, dev_mem, n, async, stream, error)
    INTEGER, INTENT(IN)                         :: n
    TYPE(dbcsr_cuda_mem_type_i8),&
         INTENT(INOUT)                          :: dev_mem
    !INTEGER(KIND=int_8), DIMENSION(1:n), INTENT(IN), TARGET :: host_mem
    INTEGER(KIND=int_8), DIMENSION(:), POINTER              :: host_mem
    LOGICAL, INTENT(IN), OPTIONAL               :: async
    INTEGER, INTENT(IN), OPTIONAL               :: stream
    TYPE(dbcsr_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_host_to_dev_l', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle, istat
    INTEGER(KIND=C_INT)                      :: async_type, stream_l
    INTEGER(KIND=int_8), POINTER                         :: host_first_el
    TYPE(C_PTR)                              :: host_mem_c_ptr

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handle, error)
    IF (careful_mod) &
         CALL dbcsr_assert (n, "LE", dev_mem%len,&
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Device memory too small for copy.",&
         __LINE__, error=error)
    IF (PRESENT (async)) THEN
       IF (async) THEN
          async_type = 1 ! Async
       ELSE
          async_type = 0 ! Sync
       ENDIF
    ELSE
       async_type = 2 ! Sync if Async not possible
    ENDIF
    IF (PRESENT (stream)) THEN
       stream_l = stream
    ELSE
       stream_l = default_stream
    ENDIF
#if defined (__DBCSR_CUDA)
    IF (n .GT. 0) THEN
       host_first_el => host_mem(1)
       host_mem_c_ptr = C_LOC (host_first_el)
    ENDIF
    IF (n .GT. 0) THEN
       istat = cuda_memcpy_h2d_cu(host_mem_c_ptr, dev_mem%ref,&
            INT(n*int_8_size, KIND=C_SIZE_T), async_type, stream_l)
    ELSE
       istat = 0
    ENDIF
#else
    istat = -1
#endif
    IF (istat /= 0) THEN
       CALL dbcsr_assert (istat, "EQ", 0,&
            dbcsr_fatal_level, dbcsr_internal_error, routineN,&
            "Error copying data from host to device memory.",&
            __LINE__, error=error)
    ENDIF
    CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE cp_host_to_dev_l


  SUBROUTINE cp_dev_to_host_l (dev_mem, host_mem, n, async, stream, error)
    TYPE(dbcsr_cuda_mem_type_i8), INTENT(IN) :: dev_mem
    INTEGER, INTENT(IN)                      :: n
    !INTEGER(KIND=int_8), DIMENSION(1:n), INTENT(OUT),&
    !  TARGET                                 :: host_mem
    INTEGER(KIND=int_8), DIMENSION(:), POINTER           :: host_mem
    LOGICAL, INTENT(IN), OPTIONAL            :: async
    INTEGER, INTENT(IN), OPTIONAL            :: stream
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dev_to_host_l', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle, istat
    INTEGER(KIND=C_INT)                      :: async_type, stream_l
    INTEGER(KIND=int_8), POINTER                         :: host_first_el
    TYPE(C_PTR)                              :: host_mem_c_ptr

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handle, error)
    IF (careful_mod) &
         CALL dbcsr_assert (n, "LE", dev_mem%len,&
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Device memory too small for copy.",&
         __LINE__, error=error)
    IF (PRESENT (async)) THEN
       IF (async) THEN
          async_type = 1 ! Async
       ELSE
          async_type = 0 ! Sync
       ENDIF
    ELSE
       async_type = 2 ! Sync if Async not possible
    ENDIF
    IF (PRESENT (stream)) THEN
       stream_l = stream
    ELSE
       stream_l = default_stream
    ENDIF
#if defined (__DBCSR_CUDA)
    IF (n .GT. 0) THEN
       host_first_el => host_mem(1)
       host_mem_c_ptr = C_LOC (host_first_el)
    ENDIF
    IF (n .GT. 0) THEN
       istat = cuda_memcpy_d2h_cu(dev_mem%ref, host_mem_c_ptr,&
            INT(n*int_8_size, KIND=C_SIZE_T), async_type, stream_l)
    ELSE
       istat = 0
    ENDIF
#else
    istat = -1
#endif
    IF (istat /= 0) THEN
       CALL dbcsr_assert (istat, "EQ", 0,&
            dbcsr_fatal_level, dbcsr_internal_error, routineN,&
            "Error copying data from device to host memory.",&
            __LINE__, error=error)
    ENDIF
    CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE cp_dev_to_host_l




  SUBROUTINE cp_host_to_dev_l0 (host_mem, dev_mem, async, stream, error)
    TYPE(dbcsr_cuda_mem_type_i8),&
         INTENT(INOUT)                          :: dev_mem
    INTEGER(KIND=int_8), POINTER                            :: host_mem
    LOGICAL, INTENT(IN), OPTIONAL               :: async
    INTEGER, INTENT(IN), OPTIONAL               :: stream
    TYPE(dbcsr_error_type), INTENT(INOUT)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_host_to_dev_l0', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle, istat
    INTEGER(KIND=C_INT)                      :: async_type, stream_l
    TYPE(C_PTR)                              :: host_mem_c_ptr

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handle, error)
    IF (careful_mod) &
         CALL dbcsr_assert (1, "LE", dev_mem%len,&
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Device memory too small for copy.",&
         __LINE__, error=error)
    IF (PRESENT (async)) THEN
       IF (async) THEN
          async_type = 1 ! Async
       ELSE
          async_type = 0 ! Sync
       ENDIF
    ELSE
       async_type = 2 ! Sync if Async not possible
    ENDIF
    IF (PRESENT (stream)) THEN
       stream_l = stream
    ELSE
       stream_l = default_stream
    ENDIF
#if defined (__DBCSR_CUDA)
    host_mem_c_ptr = C_LOC (host_mem)
    istat = cuda_memcpy_h2d_cu(host_mem_c_ptr, dev_mem%ref,&
         INT(1*int_8_size, KIND=C_SIZE_T), async_type, stream_l)
#else
    istat = -1
#endif
    IF (istat /= 0) THEN
       CALL dbcsr_assert (istat, "EQ", 0,&
            dbcsr_fatal_level, dbcsr_internal_error, routineN,&
            "Error copying data from host to device memory.",&
            __LINE__, error=error)
    ENDIF
    CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE cp_host_to_dev_l0



  SUBROUTINE cp_dev_to_host_l0 (dev_mem, host_mem, async, stream, error)
    TYPE(dbcsr_cuda_mem_type_i8), INTENT(IN) :: dev_mem
    INTEGER(KIND=int_8), POINTER                         :: host_mem
    LOGICAL, INTENT(IN), OPTIONAL            :: async
    INTEGER, INTENT(IN), OPTIONAL            :: stream
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'cp_dev_to_host_l0', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle, istat
    INTEGER(KIND=C_INT)                      :: async_type, stream_l
    TYPE(C_PTR)                              :: host_mem_c_ptr

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handle, error)
    IF (careful_mod) &
         CALL dbcsr_assert (1, "LE", dev_mem%len,&
         dbcsr_fatal_level, dbcsr_internal_error, routineN,&
         "Device memory too small for copy.",&
         __LINE__, error=error)
    IF (PRESENT (async)) THEN
       IF (async) THEN
          async_type = 1 ! Async
       ELSE
          async_type = 0 ! Sync
       ENDIF
    ELSE
       async_type = 2 ! Sync if Async not possible
    ENDIF
    IF (PRESENT (stream)) THEN
       stream_l = stream
    ELSE
       stream_l = default_stream
    ENDIF
#if defined (__DBCSR_CUDA)
    host_mem_c_ptr = C_LOC (host_mem)
    istat = cuda_memcpy_d2h_cu(dev_mem%ref, host_mem_c_ptr,&
         INT(1*int_8_size, KIND=C_SIZE_T), async_type, stream_l)
#else
    istat = -1
#endif
    IF (istat /= 0) THEN
       CALL dbcsr_assert (istat, "EQ", 0,&
            dbcsr_fatal_level, dbcsr_internal_error, routineN,&
            "Error copying data from device to host memory.",&
            __LINE__, error=error)
    ENDIF
    CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE cp_dev_to_host_l0
