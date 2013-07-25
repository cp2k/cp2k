!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2013  Urban Borstnik and the CP2K developers group   !
!-----------------------------------------------------------------------------!


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
  SUBROUTINE host_mem_alloc_z_3D (host_mem, n1, n2, n3, error)
    COMPLEX(kind=real_8), DIMENSION(:,:,:), POINTER           :: host_mem
    INTEGER, INTENT(IN)                      :: n1, n2, n3
    TYPE(dbcsr_error_type), INTENT(INOUT), OPTIONAL :: error
    TYPE(C_PTR)                              :: host_mem_c_ptr
    INTEGER                                  :: n_bytes
    n_bytes = MAX(1,n1)*MAX(1,n2)*MAX(1,n3)*(2*real_8_size)
    CALL host_mem_alloc_raw(host_mem_c_ptr,n_bytes , error)
#if defined (__DBCSR_CUDA)
    CALL C_F_POINTER (host_mem_c_ptr, host_mem, &
                               (/ MAX(1,n1),MAX(1,n2),MAX(1,n3) /))
#else
    STOP "host_mem_alloc_z_3D: DBCSR_CUDA not compiled in."
#endif
  END SUBROUTINE host_mem_alloc_z_3D


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
  SUBROUTINE host_mem_dealloc_z_3D (host_mem, error)
    COMPLEX(kind=real_8), DIMENSION(:,:,:), &
      POINTER                                :: host_mem
    TYPE(dbcsr_error_type), INTENT(INOUT), OPTIONAL :: error
    CHARACTER(len=*), PARAMETER :: routineN = 'host_mem_dealloc_z_4D', &
      routineP = moduleN//':'//routineN
    INTEGER                                  :: error_handle, istat

    IF (SIZE (host_mem) == 0) RETURN
    IF (PRESENT (error)) CALL dbcsr_error_set (routineN, error_handle, error)
#if defined (__DBCSR_CUDA)
    istat = cuda_host_mem_dealloc_cu(C_LOC(host_mem(1,1,1)))
    IF (istat /= 0 ) &
       STOP "host_mem_dealloc_z_3D: Error deallocating host pinned memory"
#else
    STOP "host_mem_dealloc_z: DBCSR_CUDA not compiled in."
#endif
     IF (PRESENT (error)) CALL dbcsr_error_stop (error_handle, error)
  END SUBROUTINE host_mem_dealloc_z_3D
  
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

