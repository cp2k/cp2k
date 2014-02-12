! *****************************************************************************
!> \brief Transfers 1D fortran-array from host to cuda devmem.
!> \author  Ole Schuett
! *****************************************************************************
 SUBROUTINE host2dev_c4_1D(this, hostmem, stream, error)
    TYPE(dbcsr_cuda_devmem_type), INTENT(IN) :: this
    COMPLEX(kind=real_4), DIMENSION(:), POINTER :: hostmem
    TYPE(dbcsr_cuda_stream_type), INTENT(IN) :: stream
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

#if ! defined (__DBCSR_CUDA)
    CALL mp_abort("__DBCSR_CUDA not compiled in.")
#else
    IF(this%datatype /= dbcsr_type_complex_4) CALL mp_abort("dbcsr_cuda_devmem_host2dev: datatype missmatch")
    CALL host2dev_raw(this, C_LOC(hostmem(1)), 2*real_4_size*SIZE(hostmem), stream, error)
#endif
 END SUBROUTINE host2dev_c4_1D


! *****************************************************************************
!> \brief Transfers 2D fortran-array from host to cuda devmem.
!> \author  Ole Schuett
! *****************************************************************************
 SUBROUTINE host2dev_c4_2D(this, hostmem, stream, error)
    TYPE(dbcsr_cuda_devmem_type), INTENT(IN) :: this
    COMPLEX(kind=real_4), DIMENSION(:, :), POINTER         :: hostmem
    TYPE(dbcsr_cuda_stream_type), INTENT(IN) :: stream
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

#if ! defined (__DBCSR_CUDA)
    CALL mp_abort("__DBCSR_CUDA not compiled in.")
#else
    IF(this%datatype /= dbcsr_type_complex_4) CALL mp_abort("dbcsr_cuda_devmem_host2dev: datatype missmatch")
    CALL host2dev_raw(this, C_LOC(hostmem(1,1)), 2*real_4_size*SIZE(hostmem), stream, error)
#endif
 END SUBROUTINE host2dev_c4_2D


! *****************************************************************************
!> \brief Transfers cuda devmem to 1D fortran-array.
!> \author  Ole Schuett
! *****************************************************************************
 SUBROUTINE dev2host_c4_1D(this, hostmem, stream, error)
    TYPE(dbcsr_cuda_devmem_type), INTENT(IN) :: this
    COMPLEX(kind=real_4), DIMENSION(:), POINTER            :: hostmem
    TYPE(dbcsr_cuda_stream_type), INTENT(IN) :: stream
    TYPE(dbcsr_error_type), INTENT(INOUT)    :: error

#if ! defined (__DBCSR_CUDA)
    CALL mp_abort("__DBCSR_CUDA not compiled in.")
#else
    IF(this%datatype /= dbcsr_type_complex_4) CALL mp_abort("dbcsr_cuda_devmem_dev2host: datatype missmatch")
    CALL dev2host_raw(this, C_LOC(hostmem(1)), 2*real_4_size*SIZE(hostmem), stream, error)
#endif
 END SUBROUTINE dev2host_c4_1D

