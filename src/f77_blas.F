MODULE f77_blas
!
! Provides interfaces and generic interfaces for BLAS as distributed
! (http://www.netlib.org/blas/) Jan 2006.
!
! provides a central point to switch between generic/poisoned generic/netlib
!
!
! See f77_blas_netlib, f77_blas_poison and f77_blas_generic for further comments
!
!
! this file is the public domain. Joost VandeVondele.
!


! USE at most one of f77_blas_netlib and f77_blas_poison

!USE f77_blas_netlib
USE f77_blas_poison

USE f77_blas_generic
USE f77_blas_extra

END MODULE f77_blas
