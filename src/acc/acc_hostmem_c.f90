!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014 the CP2K developers group                       !
!-----------------------------------------------------------------------------!


! *****************************************************************************
!> \brief Allocates 1D fortan-array as cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param n size given in terms of item-count (not bytes!)
!> \param stream ...
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE acc_hostmem_alloc_c (host_mem, n, stream)
    COMPLEX(kind=real_4), DIMENSION(:), POINTER           :: host_mem
    INTEGER, INTENT(IN)                      :: n
    TYPE(acc_stream_type), INTENT(IN)        :: stream
#if defined (__ACC)
    TYPE(C_PTR)                              :: host_mem_c_ptr

    CALL acc_hostmem_alloc_raw(host_mem_c_ptr, MAX(1,n)*(2*real_4_size), stream)
    CALL C_F_POINTER (host_mem_c_ptr, host_mem, (/ MAX(1,n) /))
#else
    STOP "acc_hostmem_alloc_c: ACC not compiled in."
#endif
  END SUBROUTINE acc_hostmem_alloc_c



! *****************************************************************************
!> \brief Allocates 2D fortan-array as cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param n1 sizes given in terms of item-count (not bytes!)
!> \param n2 sizes given in terms of item-count (not bytes!)
!> \param stream ...
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE acc_hostmem_alloc_c_2D (host_mem, n1, n2, stream)
    COMPLEX(kind=real_4), DIMENSION(:,:), POINTER           :: host_mem
    INTEGER, INTENT(IN)                      :: n1, n2
    TYPE(acc_stream_type), INTENT(IN)        :: stream
#if defined (__ACC)
    TYPE(C_PTR)                              :: host_mem_c_ptr
    INTEGER                                  :: n_bytes

    n_bytes = MAX(1,n1)*MAX(1,n2)*(2*real_4_size)
    CALL acc_hostmem_alloc_raw(host_mem_c_ptr, n_bytes, stream)
    CALL C_F_POINTER (host_mem_c_ptr, host_mem, (/ MAX(1,n1),MAX(1,n2) /))
#else
    STOP "acc_hostmem_alloc_c_2D: ACC not compiled in."
#endif
  END SUBROUTINE acc_hostmem_alloc_c_2D


! *****************************************************************************
!> \brief Allocates 3D fortan-array as cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param n1 sizes given in terms of item-count (not bytes!)
!> \param n2 sizes given in terms of item-count (not bytes!)
!> \param n3 sizes given in terms of item-count (not bytes!)
!> \param stream ...
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE acc_hostmem_alloc_c_3D (host_mem, n1, n2, n3, stream)
    COMPLEX(kind=real_4), DIMENSION(:,:,:), POINTER           :: host_mem
    INTEGER, INTENT(IN)                      :: n1, n2, n3
    TYPE(acc_stream_type), INTENT(IN)        :: stream
#if defined (__ACC)
    TYPE(C_PTR)                              :: host_mem_c_ptr
    INTEGER                                  :: n_bytes

    n_bytes = MAX(1,n1)*MAX(1,n2)*MAX(1,n3)*(2*real_4_size)
    CALL acc_hostmem_alloc_raw(host_mem_c_ptr, n_bytes, stream)
    CALL C_F_POINTER (host_mem_c_ptr, host_mem, &
                               (/ MAX(1,n1),MAX(1,n2),MAX(1,n3) /))
#else
    STOP "acc_hostmem_alloc_c_3D: ACC not compiled in."
#endif
  END SUBROUTINE acc_hostmem_alloc_c_3D


! *****************************************************************************
!> \brief Allocates 4D fortan-array as cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param n1 sizes given in terms of item-count (not bytes!)
!> \param n2 sizes given in terms of item-count (not bytes!)
!> \param n3 sizes given in terms of item-count (not bytes!)
!> \param n4 sizes given in terms of item-count (not bytes!)
!> \param stream ...
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE acc_hostmem_alloc_c_4D (host_mem, n1, n2, n3, n4, stream)
    COMPLEX(kind=real_4), DIMENSION(:,:,:,:), POINTER           :: host_mem
    INTEGER, INTENT(IN)                      :: n1, n2, n3, n4
    TYPE(acc_stream_type), INTENT(IN)        :: stream
#if defined (__ACC)
    TYPE(C_PTR)                              :: host_mem_c_ptr
    INTEGER                                  :: n_bytes

    n_bytes = MAX(1,n1)*MAX(1,n2)*MAX(1,n3)*MAX(1,n4)*(2*real_4_size)
    CALL acc_hostmem_alloc_raw(host_mem_c_ptr, n_bytes, stream)
    CALL C_F_POINTER (host_mem_c_ptr, host_mem, &
                               (/ MAX(1,n1),MAX(1,n2),MAX(1,n3),MAX(1,n4) /))
#else
    STOP "acc_hostmem_alloc_c_4D: ACC not compiled in."
#endif
  END SUBROUTINE acc_hostmem_alloc_c_4D



! *****************************************************************************
!> \brief Deallocates a 1D fortan-array, which is cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param stream ...
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE acc_hostmem_dealloc_c (host_mem, stream)
    COMPLEX(kind=real_4), DIMENSION(:), &
      POINTER                                :: host_mem
    TYPE(acc_stream_type), INTENT(IN)        :: stream
    CHARACTER(len=*), PARAMETER :: routineN = 'acc_hostmem_dealloc_c', &
      routineP = moduleN//':'//routineN

    IF (SIZE (host_mem) == 0) RETURN
#if defined (__ACC)
    CALL acc_hostmem_dealloc_raw(C_LOC(host_mem(1)), stream)
#else
    STOP "acc_hostmem_dealloc_c: ACC not compiled in."
#endif
  END SUBROUTINE acc_hostmem_dealloc_c


! *****************************************************************************
!> \brief Deallocates a 2D fortan-array, which is cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param stream ...
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE acc_hostmem_dealloc_c_2D (host_mem, stream)
    COMPLEX(kind=real_4), DIMENSION(:,:), &
      POINTER                                :: host_mem
    TYPE(acc_stream_type), INTENT(IN)        :: stream
    CHARACTER(len=*), PARAMETER :: routineN = 'acc_hostmem_dealloc_c_2D', &
      routineP = moduleN//':'//routineN

    IF (SIZE (host_mem) == 0) RETURN
#if defined (__ACC)
    CALL acc_hostmem_dealloc_raw(C_LOC(host_mem(1,1)), stream)
#else
    STOP "acc_hostmem_dealloc_c: ACC not compiled in."
#endif
  END SUBROUTINE acc_hostmem_dealloc_c_2D


! *****************************************************************************
!> \brief Deallocates a 3D fortan-array, which is cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param stream ...
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE acc_hostmem_dealloc_c_3D (host_mem, stream)
    COMPLEX(kind=real_4), DIMENSION(:,:,:), &
      POINTER                                :: host_mem
    TYPE(acc_stream_type), INTENT(IN)        :: stream
    CHARACTER(len=*), PARAMETER :: routineN = 'acc_hostmem_dealloc_c_3D', &
      routineP = moduleN//':'//routineN

    IF (SIZE (host_mem) == 0) RETURN
#if defined (__ACC)
    CALL acc_hostmem_dealloc_raw(C_LOC(host_mem(1,1,1)), stream)
#else
    STOP "acc_hostmem_dealloc_c: ACC not compiled in."
#endif
  END SUBROUTINE acc_hostmem_dealloc_c_3D


! *****************************************************************************
!> \brief Deallocates a 4D fortan-array, which is cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param stream ...
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE acc_hostmem_dealloc_c_4D (host_mem, stream)
    COMPLEX(kind=real_4), DIMENSION(:,:,:,:), &
      POINTER                                :: host_mem
    TYPE(acc_stream_type), INTENT(IN)        :: stream
    CHARACTER(len=*), PARAMETER :: routineN = 'acc_hostmem_dealloc_c_4D', &
      routineP = moduleN//':'//routineN

    IF (SIZE (host_mem) == 0) RETURN
#if defined (__ACC)
    CALL acc_hostmem_dealloc_raw(C_LOC(host_mem(1,1,1,1)), stream)
#else
    STOP "acc_hostmem_dealloc_c: ACC not compiled in."
#endif
  END SUBROUTINE acc_hostmem_dealloc_c_4D
