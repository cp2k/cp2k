! *****************************************************************************
!> \brief Transfers 1D fortran-array from host to cuda devmem.
!> \author  Ole Schuett
! *****************************************************************************
 SUBROUTINE host2dev_r4_1D(this, hostmem, stream)
    TYPE(acc_devmem_type), INTENT(IN) :: this
    REAL(kind=real_4), DIMENSION(:), POINTER :: hostmem
    TYPE(acc_stream_type), INTENT(IN) :: stream

#if ! defined (__ACC)
    STOP "__ACC not compiled in."
#else
    CALL host2dev_raw(this, C_LOC(hostmem(1)), real_4_size*SIZE(hostmem), stream)
#endif
 END SUBROUTINE host2dev_r4_1D


! *****************************************************************************
!> \brief Transfers 2D fortran-array from host to cuda devmem.
!> \author  Ole Schuett
! *****************************************************************************
 SUBROUTINE host2dev_r4_2D(this, hostmem, stream)
    TYPE(acc_devmem_type), INTENT(IN) :: this
    REAL(kind=real_4), DIMENSION(:, :), POINTER         :: hostmem
    TYPE(acc_stream_type), INTENT(IN) :: stream

#if ! defined (__ACC)
    STOP "__ACC not compiled in."
#else
    CALL host2dev_raw(this, C_LOC(hostmem(1,1)), real_4_size*SIZE(hostmem), stream)
#endif
 END SUBROUTINE host2dev_r4_2D


! *****************************************************************************
!> \brief Transfers cuda devmem to 1D fortran-array.
!> \author  Ole Schuett
! *****************************************************************************
 SUBROUTINE dev2host_r4_1D(this, hostmem, stream)
    TYPE(acc_devmem_type), INTENT(IN) :: this
    REAL(kind=real_4), DIMENSION(:), POINTER            :: hostmem
    TYPE(acc_stream_type), INTENT(IN) :: stream

#if ! defined (__ACC)
    STOP "__ACC not compiled in."
#else
    CALL dev2host_raw(this, C_LOC(hostmem(1)), real_4_size*SIZE(hostmem), stream)
#endif
 END SUBROUTINE dev2host_r4_1D

