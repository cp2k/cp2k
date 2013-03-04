!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2013  Urban Borstnik and the CP2K developers group   !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief
! *****************************************************************************
  SUBROUTINE dev_mem_alloc_z (dev_mem, n, error, stat)
    TYPE(dbcsr_cuda_mem_type_c8), &
      INTENT(OUT)                                   :: dev_mem
    INTEGER, INTENT(IN)                             :: n
    TYPE(dbcsr_error_type), INTENT(INOUT), OPTIONAL :: error
    INTEGER, INTENT(OUT), OPTIONAL                  :: stat

    CHARACTER(len=*), PARAMETER :: routineN = 'dev_mem_alloc_z', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle, istat

!   ---------------------------------------------------------------------------

    IF (PRESENT (error)) CALL dbcsr_error_set (routineN, error_handle, error)
    dev_mem%req_len = n
    dev_mem%len = MAX(1, n)
#if defined (__DBCSR_CUDA)
    istat = cuda_dev_mem_alloc_cu(dev_mem%ref,&
         INT(dev_mem%len*(2*real_8_size), KIND=C_SIZE_T))
#else
    istat = -1
#endif
    IF (istat /= 0) THEN
       dev_mem%len = -1
       STOP "dev_mem_alloc_z: Could not allocate GPU device memory"
    ENDIF
    IF (PRESENT (stat)) stat = istat
    IF (PRESENT (error)) CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE dev_mem_alloc_z


! *****************************************************************************
!> \brief
! *****************************************************************************
  SUBROUTINE dev_mem_dealloc_z (dev_mem, error, stat)
    TYPE(dbcsr_cuda_mem_type_c8), INTENT(INOUT) :: dev_mem
    TYPE(dbcsr_error_type), INTENT(INOUT), OPTIONAL      :: error
    INTEGER, INTENT(OUT), OPTIONAL                       :: stat


    CHARACTER(len=*), PARAMETER :: routineN = 'dev_mem_dealloc_z', &
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
    IF (istat /= 0) &
       STOP "dev_mem_dealloc_z: Could not deallocate GPU device memory"
    IF (PRESENT (stat)) stat = istat
    IF (PRESENT (error)) CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE dev_mem_dealloc_z


! *****************************************************************************
!> \brief
! *****************************************************************************
  SUBROUTINE dev_mem_realloc_z (dev_mem, n, error, stat, memory_crunch)
    TYPE(dbcsr_cuda_mem_type_c8), &
      INTENT(INOUT)                          :: dev_mem
    INTEGER, INTENT(IN)                      :: n
    TYPE(dbcsr_error_type), INTENT(INOUT), &
      OPTIONAL                               :: error
    INTEGER, INTENT(OUT), OPTIONAL           :: stat
    LOGICAL, INTENT(OUT), OPTIONAL           :: memory_crunch

    CHARACTER(len=*), PARAMETER :: routineN = 'dev_mem_realloc_z', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: error_handle, istat
    INTEGER(KIND=C_INT)                      :: mem_crunch
    INTEGER(KIND=C_SIZE_T)                   :: old_size

!   ---------------------------------------------------------------------------

    IF (PRESENT (error)) CALL dbcsr_error_set (routineN, error_handle, error)
    old_size = MAX(1, dev_mem%req_len) * (2*real_8_size)
    dev_mem%req_len = n
    dev_mem%len = MAX(1, n)
#if defined (__DBCSR_CUDA)
    istat = cuda_dev_mem_realloc_cu(dev_mem%ref,&
         INT(dev_mem%len*(2*real_8_size), KIND=C_SIZE_T),&
         INT(old_size, KIND=C_SIZE_T),&
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
       STOP "dev_mem_realloc_z: Could not reallocate GPU device memory"
    ENDIF
    IF (PRESENT (stat)) stat = istat
    IF (PRESENT (error)) CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE dev_mem_realloc_z


! *****************************************************************************
!> \brief
! *****************************************************************************
  SUBROUTINE host_mem_alloc_z (host_mem, n, error)
    COMPLEX(kind=real_8), DIMENSION(:), POINTER           :: host_mem
    INTEGER, INTENT(IN)                      :: n
    TYPE(dbcsr_error_type), INTENT(INOUT), OPTIONAL :: error
    TYPE(C_PTR)                              :: host_mem_c_ptr

    CALL host_mem_alloc_raw(host_mem_c_ptr, MAX(1,n)*(2*real_8_size), error)
#if defined (__DBCSR_CUDA)
    CALL C_F_POINTER (host_mem_c_ptr, host_mem, (/ MAX(1,n) /))
#else
    STOP "host_mem_alloc_z_4D: DBCSR_CUDA not compiled in."
#endif
  END SUBROUTINE host_mem_alloc_z


! *****************************************************************************
!> \brief
! *****************************************************************************
  SUBROUTINE host_mem_alloc_z_2D (host_mem, n1, n2, error)
    COMPLEX(kind=real_8), DIMENSION(:,:), POINTER           :: host_mem
    INTEGER, INTENT(IN)                      :: n1, n2
    TYPE(dbcsr_error_type), INTENT(INOUT), OPTIONAL :: error
    TYPE(C_PTR)                              :: host_mem_c_ptr
    INTEGER                                  :: n_bytes
    n_bytes = MAX(1,n1)*MAX(1,n2)*(2*real_8_size)
    CALL host_mem_alloc_raw(host_mem_c_ptr,n_bytes , error)
#if defined (__DBCSR_CUDA)
    CALL C_F_POINTER (host_mem_c_ptr, host_mem, (/ MAX(1,n1),MAX(1,n2) /))
#else
    STOP "host_mem_alloc_z_4D: DBCSR_CUDA not compiled in."
#endif
  END SUBROUTINE host_mem_alloc_z_2D


! *****************************************************************************
!> \brief
! *****************************************************************************
  SUBROUTINE host_mem_alloc_z_4D (host_mem, n1, n2, n3, n4, error)
    COMPLEX(kind=real_8), DIMENSION(:,:,:,:), POINTER           :: host_mem
    INTEGER, INTENT(IN)                      :: n1, n2, n3, n4
    TYPE(dbcsr_error_type), INTENT(INOUT), OPTIONAL :: error
    TYPE(C_PTR)                              :: host_mem_c_ptr
    INTEGER                                  :: n_bytes
    n_bytes = MAX(1,n1)*MAX(1,n2)*MAX(1,n3)*MAX(1,n4)*(2*real_8_size)
    CALL host_mem_alloc_raw(host_mem_c_ptr,n_bytes , error)
#if defined (__DBCSR_CUDA)
    CALL C_F_POINTER (host_mem_c_ptr, host_mem, &
                               (/ MAX(1,n1),MAX(1,n2),MAX(1,n3),MAX(1,n4) /))
#else
    STOP "host_mem_alloc_z_4D: DBCSR_CUDA not compiled in."
#endif
  END SUBROUTINE host_mem_alloc_z_4D



! *****************************************************************************
!> \brief
! *****************************************************************************
  SUBROUTINE host_mem_dealloc_z (host_mem, error)
    COMPLEX(kind=real_8), DIMENSION(:), &
      POINTER                                :: host_mem
    TYPE(dbcsr_error_type), INTENT(INOUT), OPTIONAL :: error
    CHARACTER(len=*), PARAMETER :: routineN = 'host_mem_dealloc_z', &
      routineP = moduleN//':'//routineN
    INTEGER                                  :: error_handle, istat

    IF (SIZE (host_mem) == 0) RETURN
    IF (PRESENT (error)) CALL dbcsr_error_set (routineN, error_handle, error)
#if defined (__DBCSR_CUDA)
    istat = cuda_host_mem_dealloc_cu(C_LOC(host_mem(1)))
    IF (istat /= 0 ) &
       STOP "host_mem_dealloc_z: Error deallocating host pinned memory"
#else
    STOP "host_mem_dealloc_z: DBCSR_CUDA not compiled in."
#endif
     IF (PRESENT (error)) CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE host_mem_dealloc_z


! *****************************************************************************
!> \brief
! *****************************************************************************
  SUBROUTINE host_mem_dealloc_z_2D (host_mem, error)
    COMPLEX(kind=real_8), DIMENSION(:,:), &
      POINTER                                :: host_mem
    TYPE(dbcsr_error_type), INTENT(INOUT), OPTIONAL :: error
    CHARACTER(len=*), PARAMETER :: routineN = 'host_mem_dealloc_z_2D', &
      routineP = moduleN//':'//routineN
    INTEGER                                  :: error_handle, istat

    IF (SIZE (host_mem) == 0) RETURN
    IF (PRESENT (error)) CALL dbcsr_error_set (routineN, error_handle, error)
#if defined (__DBCSR_CUDA)
    istat = cuda_host_mem_dealloc_cu(C_LOC(host_mem(1,1)))
    IF (istat /= 0 ) &
       STOP "host_mem_dealloc_z_2D: Error deallocating host pinned memory"
#else
    STOP "host_mem_dealloc_z: DBCSR_CUDA not compiled in."
#endif
     IF (PRESENT (error)) CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE host_mem_dealloc_z_2D


! *****************************************************************************
!> \brief
! *****************************************************************************
  SUBROUTINE host_mem_dealloc_z_4D (host_mem, error)
    COMPLEX(kind=real_8), DIMENSION(:,:,:,:), &
      POINTER                                :: host_mem
    TYPE(dbcsr_error_type), INTENT(INOUT), OPTIONAL :: error
    CHARACTER(len=*), PARAMETER :: routineN = 'host_mem_dealloc_z_4D', &
      routineP = moduleN//':'//routineN
    INTEGER                                  :: error_handle, istat

    IF (SIZE (host_mem) == 0) RETURN
    IF (PRESENT (error)) CALL dbcsr_error_set (routineN, error_handle, error)
#if defined (__DBCSR_CUDA)
    istat = cuda_host_mem_dealloc_cu(C_LOC(host_mem(1,1,1,1)))
    IF (istat /= 0 ) &
       STOP "host_mem_dealloc_z_4D: Error deallocating host pinned memory"
#else
    STOP "host_mem_dealloc_z: DBCSR_CUDA not compiled in."
#endif
     IF (PRESENT (error)) CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE host_mem_dealloc_z_4D


! *****************************************************************************
!> \brief
! *****************************************************************************
  SUBROUTINE dev_mem_zero_z (dev_mem, first, last, error)
    TYPE(dbcsr_cuda_mem_type_c8), &
      INTENT(INOUT)                          :: dev_mem
    INTEGER, INTENT(IN)                      :: first, last
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

    INTEGER(KIND=C_SIZE_T)                   :: len, off
    INTEGER                                  :: error_handle, istat
    CHARACTER(len=*), PARAMETER :: routineN = 'dev_mem_zero_z', &
      routineP = moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    CALL dbcsr_error_set (routineN, error_handle, error)
#if defined (__DBCSR_CUDA)
    IF (last .GE. first) THEN
       off = (first-1) * (2*real_8_size)
       len = (last-first+1) * (2*real_8_size)
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
  END SUBROUTINE dev_mem_zero_z
