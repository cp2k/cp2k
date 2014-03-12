!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief  routines involved in management of device memory when cuda support
!> is enabled
!> \par History
!>      02.2008 created [bgl]
!> \author Benjamin G Levine
! *****************************************************************************
! *****************************************************************************
!> \brief allocates GPU device memory to be used later.  Device memory is
!> allocated up front to eliminate the cost of allocating repeatedly later.
!> \author Benjamin G Levine
! *****************************************************************************
  SUBROUTINE cuda_d_gemm(transa, transb, m, n, k, alpha, a, lda, b,&
      ldb, beta, c, ldc)
      CHARACTER transa, transb
      INTEGER m, n, k, lda, ldb, ldc
      REAL*8 alpha, beta
      REAL*8 a(1), b(1), c(1)
#if defined(__CUBLASDP)

!      if ((m.gt.10.and.n.gt.10.and.k.gt.10).or.m.gt.1000.or.n.gt.1000.or.k.gt.1000) then
      IF ((m.gt.10.and.n.gt.10.and.k.gt.10)) THEN
!      write(*,*) 'ben before gpu_d_gemm',m,n,k
!      call flush(6)
         CALL gpu_d_gemm(transa, transb, m, n, k, alpha, a, lda, b,&
         ldb, beta, c, ldc)
      ELSE
!      write(*,*) 'ben before cpu_d_gemm',m,n,k
!      call flush(6)
         CALL cpu_d_gemm(transa, transb, m, n, k, alpha, a, lda, b,&
         ldb, beta, c, ldc)
      ENDIF

#endif

  END SUBROUTINE cuda_d_gemm

