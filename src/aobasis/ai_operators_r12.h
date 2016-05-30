!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
!> \brief Interface for the calculation of integrals over s-functions and the s-type auxiliary
!>        integrals using the Obara-Saika (OS) scheme
! **************************************************************************************************
   ABSTRACT INTERFACE
      SUBROUTINE ab_sint_os(v, nmax, zetp, zetq, zetw, rho, rac2, omega)
       USE kinds,                           ONLY: dp
      REAL(KIND=dp), DIMENSION(:, :, :), INTENT(INOUT)   :: v
      INTEGER, INTENT(IN)                                :: nmax
      REAL(KIND=dp), INTENT(IN)                          :: zetp, zetq, zetw, rho, rac2, omega

      END SUBROUTINE ab_sint_os
   END INTERFACE
