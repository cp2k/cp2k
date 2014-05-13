!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014 the CP2K developers group                       !
!-----------------------------------------------------------------------------!


! *****************************************************************************
!> \brief Allocates 1D fortan-array as cuda host-pinned memory.
!> \param n size given in terms of item-count (not bytes!)
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE acc_hostmem_alloc_r (host_mem, n)
    REAL(kind=real_4), DIMENSION(:), POINTER           :: host_mem
    INTEGER, INTENT(IN)                      :: n
    TYPE(C_PTR)                              :: host_mem_c_ptr

    CALL acc_hostmem_alloc_raw(host_mem_c_ptr, MAX(1,n)*real_4_size)
#if defined (__ACC)
    CALL C_F_POINTER (host_mem_c_ptr, host_mem, (/ MAX(1,n) /))
#else
    STOP "acc_hostmem_alloc_r_4D: ACC not compiled in."
#endif
  END SUBROUTINE acc_hostmem_alloc_r



! *****************************************************************************
!> \brief Allocates 2D fortan-array as cuda host-pinned memory.
!> \param n1,n2 sizes given in terms of item-count (not bytes!)
!> \author  Ole Schuett
! *****************************************************************************
SUBROUTINE acc_hostmem_alloc_r_2D (host_mem, n1, n2)
    REAL(kind=real_4), DIMENSION(:,:), POINTER           :: host_mem
    INTEGER, INTENT(IN)                      :: n1, n2
    TYPE(C_PTR)                              :: host_mem_c_ptr
    INTEGER                                  :: n_bytes
    n_bytes = MAX(1,n1)*MAX(1,n2)*real_4_size
    CALL acc_hostmem_alloc_raw(host_mem_c_ptr,n_bytes)
#if defined (__ACC)
    CALL C_F_POINTER (host_mem_c_ptr, host_mem, (/ MAX(1,n1),MAX(1,n2) /))
#else
    STOP "acc_hostmem_alloc_r_4D: ACC not compiled in."
#endif
  END SUBROUTINE acc_hostmem_alloc_r_2D


! *****************************************************************************
!> \brief Allocates 3D fortan-array as cuda host-pinned memory.
!> \param n1,n2,n3 sizes given in terms of item-count (not bytes!)
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE acc_hostmem_alloc_r_3D (host_mem, n1, n2, n3)
    REAL(kind=real_4), DIMENSION(:,:,:), POINTER           :: host_mem
    INTEGER, INTENT(IN)                      :: n1, n2, n3
    TYPE(C_PTR)                              :: host_mem_c_ptr
    INTEGER                                  :: n_bytes
    n_bytes = MAX(1,n1)*MAX(1,n2)*MAX(1,n3)*real_4_size
    CALL acc_hostmem_alloc_raw(host_mem_c_ptr,n_bytes)
#if defined (__ACC)
    CALL C_F_POINTER (host_mem_c_ptr, host_mem, &
                               (/ MAX(1,n1),MAX(1,n2),MAX(1,n3) /))
#else
    STOP "acc_hostmem_alloc_r_3D: ACC not compiled in."
#endif
  END SUBROUTINE acc_hostmem_alloc_r_3D


! *****************************************************************************
!> \brief Allocates 4D fortan-array as cuda host-pinned memory.
!> \param n1,..,n4 sizes given in terms of item-count (not bytes!)
!> \author  Ole Schuett
! *****************************************************************************
SUBROUTINE acc_hostmem_alloc_r_4D (host_mem, n1, n2, n3, n4)
    REAL(kind=real_4), DIMENSION(:,:,:,:), POINTER           :: host_mem
    INTEGER, INTENT(IN)                      :: n1, n2, n3, n4
    TYPE(C_PTR)                              :: host_mem_c_ptr
    INTEGER                                  :: n_bytes
    n_bytes = MAX(1,n1)*MAX(1,n2)*MAX(1,n3)*MAX(1,n4)*real_4_size
    CALL acc_hostmem_alloc_raw(host_mem_c_ptr,n_bytes)
#if defined (__ACC)
    CALL C_F_POINTER (host_mem_c_ptr, host_mem, &
                               (/ MAX(1,n1),MAX(1,n2),MAX(1,n3),MAX(1,n4) /))
#else
    STOP "acc_hostmem_alloc_r_4D: ACC not compiled in."
#endif
  END SUBROUTINE acc_hostmem_alloc_r_4D



! *****************************************************************************
!> \brief Deallocates a 1D fortan-array, which is cuda host-pinned memory.
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE acc_hostmem_dealloc_r (host_mem)
    REAL(kind=real_4), DIMENSION(:), &
      POINTER                                :: host_mem
    CHARACTER(len=*), PARAMETER :: routineN = 'acc_hostmem_dealloc_r', &
      routineP = moduleN//':'//routineN
#if defined (__ACC)
    INTEGER                                  :: istat
#endif

    IF (SIZE (host_mem) == 0) RETURN
#if defined (__ACC)
    istat = cuda_host_mem_dealloc_cu(C_LOC(host_mem(1)))
    IF (istat /= 0 ) &
       STOP "acc_hostmem_dealloc_r: Error deallocating host pinned memory"
#else
    STOP "acc_hostmem_dealloc_r: ACC not compiled in."
#endif
  END SUBROUTINE acc_hostmem_dealloc_r


! *****************************************************************************
!> \brief Deallocates a 2D fortan-array, which is cuda host-pinned memory.
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE acc_hostmem_dealloc_r_2D (host_mem)
    REAL(kind=real_4), DIMENSION(:,:), &
      POINTER                                :: host_mem
    CHARACTER(len=*), PARAMETER :: routineN = 'acc_hostmem_dealloc_r_2D', &
      routineP = moduleN//':'//routineN
#if defined (__ACC)
    INTEGER                                  :: istat
#endif

    IF (SIZE (host_mem) == 0) RETURN
#if defined (__ACC)
    istat = cuda_host_mem_dealloc_cu(C_LOC(host_mem(1,1)))
    IF (istat /= 0 ) &
       STOP "acc_hostmem_dealloc_r_2D: Error deallocating host pinned memory"
#else
    STOP "acc_hostmem_dealloc_r: ACC not compiled in."
#endif
  END SUBROUTINE acc_hostmem_dealloc_r_2D


! *****************************************************************************
!> \brief Deallocates a 3D fortan-array, which is cuda host-pinned memory.
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE acc_hostmem_dealloc_r_3D (host_mem)
    REAL(kind=real_4), DIMENSION(:,:,:), &
      POINTER                                :: host_mem
    CHARACTER(len=*), PARAMETER :: routineN = 'acc_hostmem_dealloc_r_3D', &
      routineP = moduleN//':'//routineN
#if defined (__ACC)
    INTEGER                                  :: istat
#endif

    IF (SIZE (host_mem) == 0) RETURN
#if defined (__ACC)
    istat = cuda_host_mem_dealloc_cu(C_LOC(host_mem(1,1,1)))
    IF (istat /= 0 ) &
       STOP "acc_hostmem_dealloc_r_3D: Error deallocating host pinned memory"
#else
    STOP "acc_hostmem_dealloc_r: ACC not compiled in."
#endif
  END SUBROUTINE acc_hostmem_dealloc_r_3D


! *****************************************************************************
!> \brief Deallocates a 4D fortan-array, which is cuda host-pinned memory.
!> \author  Ole Schuett
! *****************************************************************************
  SUBROUTINE acc_hostmem_dealloc_r_4D (host_mem)
    REAL(kind=real_4), DIMENSION(:,:,:,:), &
      POINTER                                :: host_mem
    CHARACTER(len=*), PARAMETER :: routineN = 'acc_hostmem_dealloc_r_4D', &
      routineP = moduleN//':'//routineN
#if defined (__ACC)
    INTEGER                                  :: istat
#endif

    IF (SIZE (host_mem) == 0) RETURN
#if defined (__ACC)
    istat = cuda_host_mem_dealloc_cu(C_LOC(host_mem(1,1,1,1)))
    IF (istat /= 0 ) &
       STOP "acc_hostmem_dealloc_r_4D: Error deallocating host pinned memory"
#else
    STOP "acc_hostmem_dealloc_r: ACC not compiled in."
#endif
  END SUBROUTINE acc_hostmem_dealloc_r_4D
