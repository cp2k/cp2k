! *****************************************************************************
!> \brief Transfers 1D fortran-array from host to cuda devmem.
!> \param[in] this device memory
!> \param hostmem host memory
!> \param[in] stream stream
!> \author  Ole Schuett
! *****************************************************************************
 SUBROUTINE host2dev_c8_1D(this, hostmem, stream)
    TYPE(acc_devmem_type), INTENT(IN) :: this
    COMPLEX(kind=real_8), DIMENSION(:), POINTER :: hostmem
    TYPE(acc_stream_type), INTENT(IN) :: stream

#if ! defined (__ACC)
    STOP "__ACC not compiled in."
#else
    CALL host2dev_raw(this, C_LOC(hostmem(1)), 2*real_8_size*SIZE(hostmem), stream)
#endif
 END SUBROUTINE host2dev_c8_1D


! *****************************************************************************
!> \brief Transfers 2D fortran-array from host to cuda devmem.
!> \param[in] this device memory
!> \param hostmem host memory
!> \param[in] stream stream
!> \author  Ole Schuett
! *****************************************************************************
 SUBROUTINE host2dev_c8_2D(this, hostmem, stream)
    TYPE(acc_devmem_type), INTENT(IN) :: this
    COMPLEX(kind=real_8), DIMENSION(:, :), POINTER         :: hostmem
    TYPE(acc_stream_type), INTENT(IN) :: stream

#if ! defined (__ACC)
    STOP "__ACC not compiled in."
#else
    CALL host2dev_raw(this, C_LOC(hostmem(1,1)), 2*real_8_size*SIZE(hostmem), stream)
#endif
 END SUBROUTINE host2dev_c8_2D


! *****************************************************************************
!> \brief Transfers cuda devmem to 1D fortran-array.
!> \param[in] this device memory
!> \param hostmem host memory
!> \param[in] stream stream
!> \author  Ole Schuett
! *****************************************************************************
 SUBROUTINE dev2host_c8_1D(this, hostmem, stream)
    TYPE(acc_devmem_type), INTENT(IN) :: this
    COMPLEX(kind=real_8), DIMENSION(:), POINTER            :: hostmem
    TYPE(acc_stream_type), INTENT(IN) :: stream

#if ! defined (__ACC)
    STOP "__ACC not compiled in."
#else
    CALL dev2host_raw(this, C_LOC(hostmem(1)), 2*real_8_size*SIZE(hostmem), stream)
#endif
 END SUBROUTINE dev2host_c8_1D

