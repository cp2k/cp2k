!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2017  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief   Accelerator support
!> \author  Ole Schuett
!> \date    2013-04
! **************************************************************************************************
MODULE acc_hostmem
#if defined (__ACC)
    USE ISO_C_BINDING, ONLY: C_INT, C_SIZE_T, C_PTR, C_LOC, C_F_POINTER
#endif
  USE acc_kinds,                       ONLY: int_4,&
                                             int_4_size,&
                                             int_8,&
                                             int_8_size,&
                                             real_4,&
                                             real_4_size,&
                                             real_8,&
                                             real_8_size
  USE acc_stream,                      ONLY: acc_stream_associated,&
                                             acc_stream_cptr,&
                                             acc_stream_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'acc_hostmem'

   LOGICAL, PARAMETER :: careful_mod = .TRUE.

   PUBLIC :: acc_hostmem_allocate, acc_hostmem_deallocate

   INTERFACE acc_hostmem_allocate
      MODULE PROCEDURE acc_hostmem_alloc_i4, acc_hostmem_alloc_i8
      MODULE PROCEDURE acc_hostmem_alloc_r4, acc_hostmem_alloc_r8
      MODULE PROCEDURE acc_hostmem_alloc_c4, acc_hostmem_alloc_c8
      MODULE PROCEDURE acc_hostmem_alloc_i4_2D, acc_hostmem_alloc_i8_2D
      MODULE PROCEDURE acc_hostmem_alloc_r4_2D, acc_hostmem_alloc_r8_2D
      MODULE PROCEDURE acc_hostmem_alloc_c4_2D, acc_hostmem_alloc_c8_2D
   END INTERFACE

   INTERFACE acc_hostmem_deallocate
      MODULE PROCEDURE acc_hostmem_dealloc_i4, acc_hostmem_dealloc_i8
      MODULE PROCEDURE acc_hostmem_dealloc_r4, acc_hostmem_dealloc_r8
      MODULE PROCEDURE acc_hostmem_dealloc_c4, acc_hostmem_dealloc_c8
      MODULE PROCEDURE acc_hostmem_dealloc_i4_2D, acc_hostmem_dealloc_i8_2D
      MODULE PROCEDURE acc_hostmem_dealloc_r4_2D, acc_hostmem_dealloc_r8_2D
      MODULE PROCEDURE acc_hostmem_dealloc_c4_2D, acc_hostmem_dealloc_c8_2D
   END INTERFACE

#if defined (__ACC)

   INTERFACE
      FUNCTION cuda_host_mem_alloc_cu(mem, n, stream_ptr) RESULT(istat) BIND(C, name="acc_host_mem_allocate")
         IMPORT
         TYPE(C_PTR)                              :: mem
         INTEGER(KIND=C_SIZE_T), INTENT(IN), &
            VALUE                                  :: n
         TYPE(C_PTR), VALUE                       :: stream_ptr
         INTEGER(KIND=C_INT)                      :: istat

      END FUNCTION cuda_host_mem_alloc_cu
   END INTERFACE

   INTERFACE
      FUNCTION cuda_host_mem_dealloc_cu(mem, stream_ptr) RESULT(istat) bind(C, name="acc_host_mem_deallocate")
         IMPORT
         TYPE(C_PTR), VALUE                       :: mem, stream_ptr
         INTEGER(KIND=C_INT)                      :: istat

      END FUNCTION cuda_host_mem_dealloc_cu
   END INTERFACE

#endif

CONTAINS


! **************************************************************************************************
!> \brief Helper-routine performing allocation of host-pinned cuda memory.
!> \param host_mem_c_ptr pointer to allocated memory
!> \param n_bytes number of bytes to allocate
!> \param stream ...
! **************************************************************************************************
#if defined (__ACC)
   SUBROUTINE acc_hostmem_alloc_raw(host_mem_c_ptr, n_bytes, stream)
      TYPE(C_PTR), INTENT(OUT)                           :: host_mem_c_ptr
      INTEGER, INTENT(IN)                                :: n_bytes
      TYPE(acc_stream_type), INTENT(IN)                  :: stream

      CHARACTER(len=*), PARAMETER :: routineN = 'acc_hostmem_alloc_raw', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: istat
      TYPE(C_PTR)                                        :: stream_cptr

      IF (.NOT. acc_stream_associated(stream)) &
         CPABORT("acc_hostmem_alloc_raw: stream not associated")

      stream_cptr = acc_stream_cptr(stream)

      istat = cuda_host_mem_alloc_cu(host_mem_c_ptr, INT(n_bytes, KIND=C_SIZE_T), stream_cptr)
      IF (istat /= 0) &
         CPABORT("acc_hostmem_alloc_raw: Could not allocate host pinned memory")
   END SUBROUTINE acc_hostmem_alloc_raw
#endif

#if defined (__ACC)
! **************************************************************************************************
!> \brief ...
!> \param host_mem_c_ptr ...
!> \param stream ...
! **************************************************************************************************
   SUBROUTINE acc_hostmem_dealloc_raw(host_mem_c_ptr, stream)
      TYPE(C_PTR), INTENT(IN)                            :: host_mem_c_ptr
      TYPE(acc_stream_type), INTENT(IN)                  :: stream

      CHARACTER(len=*), PARAMETER :: routineN = 'acc_hostmem_dealloc_raw', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: istat
      TYPE(C_PTR)                                        :: stream_cptr

      IF (.NOT. acc_stream_associated(stream)) &
         CPABORT("acc_hostmem_dealloc_raw: stream not associated")

      stream_cptr = acc_stream_cptr(stream)

      istat = cuda_host_mem_dealloc_cu(host_mem_c_ptr, stream_cptr)
      IF (istat /= 0) &
         CPABORT("acc_hostmem_dealloc_raw: Could not deallocate host pinned memory")
   END SUBROUTINE acc_hostmem_dealloc_raw
#endif


#:set instances = [ ('i4', 'int_4_size',     'INTEGER(kind=int_4)'), &
                 &  ('i8', 'int_8_size',     'INTEGER(kind=int_8)'), &
                 &  ('r4', 'real_4_size',    'REAL(kind=real_4)'), &
                 &  ('r8', 'real_8_size',    'REAL(kind=real_8)'), &
                 &  ('c4', 'complex_4_size', 'COMPLEX(kind=real_4)'), &
                 &  ('c8', 'complex_8_size', 'COMPLEX(kind=real_8)') ]

#:for nametype, size, type in instances

! **************************************************************************************************
!> \brief Allocates 1D fortan-array as cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param n size given in terms of item-count (not bytes!)
!> \param stream ...
!> \author  Ole Schuett
! **************************************************************************************************
  SUBROUTINE acc_hostmem_alloc_${nametype}$ (host_mem, n, stream)
    ${type}$, DIMENSION(:), POINTER          :: host_mem
    INTEGER, INTENT(IN)                      :: n
    TYPE(acc_stream_type), INTENT(IN)        :: stream
#if defined (__ACC)
    TYPE(C_PTR)                              :: host_mem_c_ptr

    CALL acc_hostmem_alloc_raw(host_mem_c_ptr, MAX(1,n)*${size}$, stream)
    CALL C_F_POINTER (host_mem_c_ptr, host_mem, (/ MAX(1,n) /))
#else
    MARK_USED(host_mem)
    MARK_USED(n)
    MARK_USED(stream)
    CPABORT("acc_hostmem_alloc_${nametype}$: ACC not compiled in.")
#endif
  END SUBROUTINE acc_hostmem_alloc_${nametype}$



! **************************************************************************************************
!> \brief Allocates 2D fortan-array as cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param n1 sizes given in terms of item-count (not bytes!)
!> \param n2 sizes given in terms of item-count (not bytes!)
!> \param stream ...
!> \author  Ole Schuett
! **************************************************************************************************
  SUBROUTINE acc_hostmem_alloc_${nametype}$_2D (host_mem, n1, n2, stream)
    ${type}$, DIMENSION(:,:), POINTER        :: host_mem
    INTEGER, INTENT(IN)                      :: n1, n2
    TYPE(acc_stream_type), INTENT(IN)        :: stream
#if defined (__ACC)
    TYPE(C_PTR)                              :: host_mem_c_ptr
    INTEGER                                  :: n_bytes

    n_bytes = MAX(1,n1)*MAX(1,n2)*${size}$
    CALL acc_hostmem_alloc_raw(host_mem_c_ptr, n_bytes, stream)
    CALL C_F_POINTER (host_mem_c_ptr, host_mem, (/ MAX(1,n1),MAX(1,n2) /))
#else
    MARK_USED(host_mem)
    MARK_USED(n1)
    MARK_USED(n2)
    MARK_USED(stream)
    CPABORT("acc_hostmem_alloc_${nametype}$_2D: ACC not compiled in.")
#endif
  END SUBROUTINE acc_hostmem_alloc_${nametype}$_2D


! **************************************************************************************************
!> \brief Deallocates a 1D fortan-array, which is cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param stream ...
!> \author  Ole Schuett
! **************************************************************************************************
  SUBROUTINE acc_hostmem_dealloc_${nametype}$ (host_mem, stream)
    ${type}$, DIMENSION(:), POINTER          :: host_mem
    TYPE(acc_stream_type), INTENT(IN)        :: stream
    CHARACTER(len=*), PARAMETER :: routineN = 'acc_hostmem_dealloc_${nametype}$', &
      routineP = moduleN//':'//routineN

    IF (SIZE (host_mem) == 0) RETURN
#if defined (__ACC)
    CALL acc_hostmem_dealloc_raw(C_LOC(host_mem(1)), stream)
#else
    MARK_USED(host_mem)
    MARK_USED(stream)
    CPABORT("acc_hostmem_dealloc_${nametype}$: ACC not compiled in.")
#endif
  END SUBROUTINE acc_hostmem_dealloc_${nametype}$


! **************************************************************************************************
!> \brief Deallocates a 2D fortan-array, which is cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param stream ...
!> \author  Ole Schuett
! **************************************************************************************************
  SUBROUTINE acc_hostmem_dealloc_${nametype}$_2D (host_mem, stream)
    ${type}$, DIMENSION(:,:), POINTER        :: host_mem
    TYPE(acc_stream_type), INTENT(IN)        :: stream
    CHARACTER(len=*), PARAMETER :: routineN = 'acc_hostmem_dealloc_${nametype}$_2D', &
      routineP = moduleN//':'//routineN

    IF (SIZE (host_mem) == 0) RETURN
#if defined (__ACC)
    CALL acc_hostmem_dealloc_raw(C_LOC(host_mem(1,1)), stream)
#else
    MARK_USED(host_mem)
    MARK_USED(stream)
    CPABORT("acc_hostmem_dealloc_${nametype}$: ACC not compiled in.")
#endif
  END SUBROUTINE acc_hostmem_dealloc_${nametype}$_2D

#:endfor

END MODULE acc_hostmem
