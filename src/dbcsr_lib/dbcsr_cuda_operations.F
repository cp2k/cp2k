!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group   !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   CUDA support for DBCSR
!> \author  Urban Borstnik
!> \date    2011-04-06
!> \version 1.0
!>
!> <b>Modification history:</b>
!> - Created 2011-04-06
! *****************************************************************************
MODULE dbcsr_cuda_operations
  USE dbcsr_cuda_devmem,               ONLY: dbcsr_cuda_devmem_datatype,&
                                             dbcsr_cuda_devmem_type
  USE dbcsr_cuda_stream,               ONLY: dbcsr_cuda_stream_type
  USE dbcsr_error_handling
  USE dbcsr_mm_types,                  ONLY: dbcsr_ps_width
  USE dbcsr_message_passing,           ONLY: mp_abort

#if defined (__DBCSR_CUDA)
  USE ISO_C_BINDING
  USE dbcsr_cuda_devmem,               ONLY: dbcsr_cuda_devmem_cptr
  USE dbcsr_cuda_stream,               ONLY: dbcsr_cuda_stream_cptr
#endif

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_cuda_operations'

  LOGICAL, PARAMETER :: careful_mod = .FALSE.


  PUBLIC :: dbcsr_cuda_do_mm_stack, dbcsr_cuda_transpose


#if defined (__DBCSR_CUDA)
  INTERFACE
     FUNCTION dc_do_stack_cu(param_stack, stack_size, nparams,&
          data_type,&
          a_data, b_data, c_data, m_max, n_max, k_max, def_mnk,&
          stream_ptr)&
          RESULT (istat) &
          BIND(C, name="dc_do_stack_cu")
       USE ISO_C_BINDING
    TYPE(C_PTR), INTENT(IN), VALUE           :: param_stack
    INTEGER(KIND=C_INT), INTENT(IN), VALUE   :: stack_size, nparams, data_type
    TYPE(C_PTR), INTENT(IN), VALUE           :: a_data, b_data, c_data
    INTEGER(KIND=C_INT), INTENT(IN), VALUE   :: m_max, n_max, k_max, def_mnk
    TYPE(C_PTR), VALUE                       :: stream_ptr
    INTEGER(KIND=C_INT)                      :: istat

     END FUNCTION dc_do_stack_cu
  END INTERFACE
  
  INTERFACE
     FUNCTION dc_do_transpose_stack_cu(trs_stack, offset, nblks, buffer, data_type, m, n, stream_ptr)&
          RESULT (istat) &
          BIND(C, name="dc_do_transpose_stack_cu")
       USE ISO_C_BINDING
    TYPE(C_PTR), INTENT(IN), VALUE           :: trs_stack
    INTEGER(KIND=C_INT), INTENT(IN), VALUE   :: offset, nblks
    TYPE(C_PTR), INTENT(IN), VALUE           :: buffer
    INTEGER(KIND=C_INT), INTENT(IN), VALUE   :: data_type, m, n
    TYPE(C_PTR), VALUE                       :: stream_ptr
    INTEGER(KIND=C_INT)                      :: istat

     END FUNCTION dc_do_transpose_stack_cu
  END INTERFACE
#endif

CONTAINS

! *****************************************************************************
!> \brief Launch a cuda kernel for processing a stack.
! *****************************************************************************
  SUBROUTINE dbcsr_cuda_do_mm_stack(param_stack, stack_size,&
       a_data, b_data, c_data, m_max, n_max, k_max, def_mnk, stream, success, error)
    TYPE(dbcsr_cuda_devmem_type), INTENT(IN) :: param_stack
    INTEGER, INTENT(IN)                      :: stack_size
    TYPE(dbcsr_cuda_devmem_type), INTENT(IN) :: a_data, b_data
    TYPE(dbcsr_cuda_devmem_type), &
      INTENT(INOUT)                          :: c_data
    INTEGER, INTENT(IN)                      :: m_max, n_max, k_max
    LOGICAL, INTENT(IN)                      :: def_mnk
    TYPE(dbcsr_cuda_stream_type), INTENT(IN) :: stream
    LOGICAL, INTENT(INOUT)                   :: success
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

#if ! defined (__DBCSR_CUDA)
    CALL mp_abort("__DBCSR_CUDA not compiled in.")
#else

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_cuda_do_mm_stack', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: datatype, error_handle, istat
    INTEGER(KIND=C_INT)                      :: mnk

   IF (careful_mod) CALL dbcsr_error_set (routineN, error_handle, error)

   datatype = dbcsr_cuda_devmem_datatype(a_data)

   IF (datatype /= dbcsr_cuda_devmem_datatype(b_data) .OR. &
       datatype /= dbcsr_cuda_devmem_datatype(c_data))&
        STOP "dbcsr_cuda_do_mm_stack: datatype missmatch"


    mnk = 0
    IF (def_mnk) mnk = 1

    istat = dc_do_stack_cu(dbcsr_cuda_devmem_cptr(param_stack),&
                             INT(stack_size, KIND=C_INT),&
                             INT(dbcsr_ps_width, KIND=C_INT),&
                             INT(datatype, KIND=C_INT),&
                             dbcsr_cuda_devmem_cptr(a_data),&
                             dbcsr_cuda_devmem_cptr(b_data),&
                             dbcsr_cuda_devmem_cptr(c_data),&
                             INT(m_max,KIND=C_INT),&
                             INT(n_max,KIND=C_INT),&
                             INT(k_max,KIND=C_INT),&
                             mnk, dbcsr_cuda_stream_cptr(stream))

    IF (istat == 100) STOP "dbcsr_cuda_do_mm_stack: failed due to cudaGetLastError()."

    success = (istat == 0) ! false if no suitable kernel was found

    IF (careful_mod) CALL dbcsr_error_stop (error_handle, error)
#endif
  END SUBROUTINE dbcsr_cuda_do_mm_stack


! *****************************************************************************
!> \brief Launch a cuda transpose kernel
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE dbcsr_cuda_transpose(trs_stack, offset, nblks, buffer, m, n, stream, error)
    TYPE(dbcsr_cuda_devmem_type), INTENT(IN) :: trs_stack
    INTEGER, INTENT(IN)                      :: offset
    INTEGER, INTENT(IN)                      :: nblks
    TYPE(dbcsr_cuda_devmem_type), INTENT(IN) :: buffer
    INTEGER, INTENT(IN)                      :: m, n
    TYPE(dbcsr_cuda_stream_type), INTENT(IN) :: stream
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

#if ! defined (__DBCSR_CUDA)
    CALL mp_abort("__DBCSR_CUDA not compiled in.")
#else

    CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_cuda_transpose', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: datatype, error_handle, istat

    IF (careful_mod) CALL dbcsr_error_set (routineN, error_handle, error)

    datatype = dbcsr_cuda_devmem_datatype(buffer)

    istat = dc_do_transpose_stack_cu(dbcsr_cuda_devmem_cptr(trs_stack),&
                             INT(offset, KIND=C_INT),&
                             INT(nblks, KIND=C_INT),&
                             dbcsr_cuda_devmem_cptr(buffer),&
                             INT(datatype, KIND=C_INT),&
                             INT(m,KIND=C_INT),&
                             INT(n,KIND=C_INT),&
                             dbcsr_cuda_stream_cptr(stream))

    IF (istat /= 0) STOP "dbcsr_cuda_transpose: something went wrong."
    IF (careful_mod) CALL dbcsr_error_stop (error_handle, error)
#endif
  END SUBROUTINE dbcsr_cuda_transpose

END MODULE dbcsr_cuda_operations
