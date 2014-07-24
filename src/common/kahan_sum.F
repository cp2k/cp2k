!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2014  CP2K developers group                          !
!                                                                             !
!   This program is free software; you can redistribute it and/or modify      !
!   it under the terms of the GNU General Public License as published by      !
!   the Free Software Foundation; either version 2 of the License, or         !
!   (at your option) any later version.                                       !
!                                                                             !
!   This program is distributed in the hope that it will be useful,           !
!   but WITHOUT ANY WARRANTY; without even the implied warranty of            !
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             !
!   GNU General Public License for more details.                              !
!                                                                             !
!   You should have received a copy of the GNU General Public License         !
!   along with this program; if not, write to the Free Software               !
!   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                 !
!                                                                             !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief sums arrays of real/complex numbers with *much* reduced round-off as compared to
!>      a naive implementation (or the one found in most compiler's SUM intrinsic)
!>      using an implementation of Kahan's algorithm for summing real numbers
!>      that can be used instead of the standard Fortran SUM(array[,mask]).
!>
!>      see also http://en.wikipedia.org/wiki/Kahan_summation_algorithm
!> \note
!>      if the compiler optimises away the 'tricky' bit, no accuracy is gained,
!>      if the compiler uses extended precision inconsistently even worse results might be obtained.
!>      This has not been observed.
!>      This algorithm is not fast, and thus not recommended for cases where round-off is not a
!>      concern but performance is.
!>
!>      the standard intrinsic sum can be 'replaced' using the following use statement
!>
!>      USE kahan_sum, ONLY: sum => kahan_sum
!> \par History
!>      03.2006 [Joost VandeVondele]
!> \author Joost VandeVondele
! *****************************************************************************
MODULE kahan_sum
  USE f77_blas

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: accurate_sum
  INTEGER, PARAMETER :: sp=KIND(0.0),dp=KIND(0.0D0)
  REAL(KIND=sp),    PARAMETER :: szero=0.0_sp
  REAL(KIND=dp),    PARAMETER :: dzero=0.0_dp
  COMPLEX(KIND=sp), PARAMETER :: czero=(0.0_sp,0.0_sp)
  COMPLEX(KIND=dp), PARAMETER :: zzero=(0.0_dp,0.0_dp)
  INTERFACE accurate_sum
    MODULE PROCEDURE &
     kahan_sum_s1,      kahan_sum_d1,      kahan_sum_c1,      kahan_sum_z1,&
     kahan_sum_s2,      kahan_sum_d2,      kahan_sum_c2,      kahan_sum_z2,&
     kahan_sum_s3,      kahan_sum_d3,      kahan_sum_c3,      kahan_sum_z3,&
     kahan_sum_s4,      kahan_sum_d4,      kahan_sum_c4,      kahan_sum_z4,&
     kahan_sum_s5,      kahan_sum_d5,      kahan_sum_c5,      kahan_sum_z5,&
     kahan_sum_s6,      kahan_sum_d6,      kahan_sum_c6,      kahan_sum_z6,&
     kahan_sum_s7,      kahan_sum_d7,      kahan_sum_c7,      kahan_sum_z7
  END INTERFACE accurate_sum
 CONTAINS
! *****************************************************************************
!> \brief ...
!> \param array ...
!> \param mask ...
!> \retval ks ...
! *****************************************************************************
  FUNCTION kahan_sum_s1(array,mask) RESULT(ks)
    REAL(KIND=sp), DIMENSION(:), INTENT(IN)  :: array
    LOGICAL, DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: mask
    REAL(KIND=sp)                            :: ks

    INTEGER                                  :: i1
    REAL(KIND=sp)                            :: c, t, y

    ks = szero ; t = szero ; y = szero ; c = szero

    IF (PRESENT(mask)) THEN
      DO i1=1,SIZE(array,1)
         IF (mask(i1                  )) THEN
           y  = array(i1                  ) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
         END IF
      ENDDO
    ELSE
      DO i1=1,SIZE(array,1)
           y  = array(i1                  ) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
      ENDDO
    ENDIF
  END FUNCTION kahan_sum_s1

! *****************************************************************************
!> \brief ...
!> \param array ...
!> \param mask ...
!> \retval ks ...
! *****************************************************************************
  FUNCTION kahan_sum_d1(array,mask) RESULT(ks)
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: array
    LOGICAL, DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: mask
    REAL(KIND=dp)                            :: ks

    INTEGER                                  :: i1
    REAL(KIND=dp)                            :: c, t, y

    ks = dzero ; t = dzero ; y = dzero ; c = dzero

    IF (PRESENT(mask)) THEN
      DO i1=1,SIZE(array,1)
         IF (mask(i1                  )) THEN
           y  = array(i1                  ) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
         END IF
      ENDDO
    ELSE
      DO i1=1,SIZE(array,1)
           y  = array(i1                  ) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
      ENDDO
    ENDIF
  END FUNCTION kahan_sum_d1

! *****************************************************************************
!> \brief ...
!> \param array ...
!> \param mask ...
!> \retval ks ...
! *****************************************************************************
  FUNCTION kahan_sum_c1(array,mask) RESULT(ks)
    COMPLEX(KIND=sp), DIMENSION(:), &
      INTENT(IN)                             :: array
    LOGICAL, DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: mask
    COMPLEX(KIND=sp)                         :: ks

    COMPLEX(KIND=sp)                         :: c, t, y
    INTEGER                                  :: i1

    ks = czero ; t = czero ; y = czero ; c = czero

    IF (PRESENT(mask)) THEN
      DO i1=1,SIZE(array,1)
         IF (mask(i1                  )) THEN
           y  = array(i1                  ) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
         END IF
      ENDDO
    ELSE
      DO i1=1,SIZE(array,1)
           y  = array(i1                  ) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
      ENDDO
    ENDIF
  END FUNCTION kahan_sum_c1

! *****************************************************************************
!> \brief ...
!> \param array ...
!> \param mask ...
!> \retval ks ...
! *****************************************************************************
  FUNCTION kahan_sum_z1(array,mask) RESULT(ks)
    COMPLEX(KIND=dp), DIMENSION(:), &
      INTENT(IN)                             :: array
    LOGICAL, DIMENSION(:), INTENT(IN), &
      OPTIONAL                               :: mask
    COMPLEX(KIND=dp)                         :: ks

    COMPLEX(KIND=dp)                         :: c, t, y
    INTEGER                                  :: i1

    ks = zzero ; t = zzero ; y = zzero ; c = zzero

    IF (PRESENT(mask)) THEN
      DO i1=1,SIZE(array,1)
         IF (mask(i1                  )) THEN
           y  = array(i1                  ) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
         END IF
      ENDDO
    ELSE
      DO i1=1,SIZE(array,1)
           y  = array(i1                  ) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
      ENDDO
    ENDIF
  END FUNCTION kahan_sum_z1

! *****************************************************************************
!> \brief ...
!> \param array ...
!> \param mask ...
!> \retval ks ...
! *****************************************************************************
  FUNCTION kahan_sum_s2(array,mask) RESULT(ks)
    REAL(KIND=sp), DIMENSION(:, :), &
      INTENT(IN)                             :: array
    LOGICAL, DIMENSION(:, :), INTENT(IN), &
      OPTIONAL                               :: mask
    REAL(KIND=sp)                            :: ks

    INTEGER                                  :: i1, i2
    REAL(KIND=sp)                            :: c, t, y

    ks = szero ; t = szero ; y = szero ; c = szero

    IF (PRESENT(mask)) THEN
      DO i2=1,SIZE(array,2)
      DO i1=1,SIZE(array,1)
         IF (mask(i1,i2               )) THEN
           y  = array(i1,i2               ) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
         END IF
      ENDDO
      ENDDO
    ELSE
      DO i2=1,SIZE(array,2)
      DO i1=1,SIZE(array,1)
           y  = array(i1,i2               ) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
      ENDDO
      ENDDO
    ENDIF
  END FUNCTION kahan_sum_s2

! *****************************************************************************
!> \brief ...
!> \param array ...
!> \param mask ...
!> \retval ks ...
! *****************************************************************************
  FUNCTION kahan_sum_d2(array,mask) RESULT(ks)
    REAL(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: array
    LOGICAL, DIMENSION(:, :), INTENT(IN), &
      OPTIONAL                               :: mask
    REAL(KIND=dp)                            :: ks

    INTEGER                                  :: i1, i2
    REAL(KIND=dp)                            :: c, t, y

    ks = dzero ; t = dzero ; y = dzero ; c = dzero

    IF (PRESENT(mask)) THEN
      DO i2=1,SIZE(array,2)
      DO i1=1,SIZE(array,1)
         IF (mask(i1,i2               )) THEN
           y  = array(i1,i2               ) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
         END IF
      ENDDO
      ENDDO
    ELSE
      DO i2=1,SIZE(array,2)
      DO i1=1,SIZE(array,1)
           y  = array(i1,i2               ) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
      ENDDO
      ENDDO
    ENDIF
  END FUNCTION kahan_sum_d2

! *****************************************************************************
!> \brief ...
!> \param array ...
!> \param mask ...
!> \retval ks ...
! *****************************************************************************
  FUNCTION kahan_sum_c2(array,mask) RESULT(ks)
    COMPLEX(KIND=sp), DIMENSION(:, :), &
      INTENT(IN)                             :: array
    LOGICAL, DIMENSION(:, :), INTENT(IN), &
      OPTIONAL                               :: mask
    COMPLEX(KIND=sp)                         :: ks

    COMPLEX(KIND=sp)                         :: c, t, y
    INTEGER                                  :: i1, i2

    ks = czero ; t = czero ; y = czero ; c = czero

    IF (PRESENT(mask)) THEN
      DO i2=1,SIZE(array,2)
      DO i1=1,SIZE(array,1)
         IF (mask(i1,i2               )) THEN
           y  = array(i1,i2               ) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
         END IF
      ENDDO
      ENDDO
    ELSE
      DO i2=1,SIZE(array,2)
      DO i1=1,SIZE(array,1)
           y  = array(i1,i2               ) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
      ENDDO
      ENDDO
    ENDIF
  END FUNCTION kahan_sum_c2

! *****************************************************************************
!> \brief ...
!> \param array ...
!> \param mask ...
!> \retval ks ...
! *****************************************************************************
  FUNCTION kahan_sum_z2(array,mask) RESULT(ks)
    COMPLEX(KIND=dp), DIMENSION(:, :), &
      INTENT(IN)                             :: array
    LOGICAL, DIMENSION(:, :), INTENT(IN), &
      OPTIONAL                               :: mask
    COMPLEX(KIND=dp)                         :: ks

    COMPLEX(KIND=dp)                         :: c, t, y
    INTEGER                                  :: i1, i2

    ks = zzero ; t = zzero ; y = zzero ; c = zzero

    IF (PRESENT(mask)) THEN
      DO i2=1,SIZE(array,2)
      DO i1=1,SIZE(array,1)
         IF (mask(i1,i2               )) THEN
           y  = array(i1,i2               ) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
         END IF
      ENDDO
      ENDDO
    ELSE
      DO i2=1,SIZE(array,2)
      DO i1=1,SIZE(array,1)
           y  = array(i1,i2               ) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
      ENDDO
      ENDDO
    ENDIF
  END FUNCTION kahan_sum_z2

! *****************************************************************************
!> \brief ...
!> \param array ...
!> \param mask ...
!> \retval ks ...
! *****************************************************************************
  FUNCTION kahan_sum_s3(array,mask) RESULT(ks)
    REAL(KIND=sp), DIMENSION(:, :, :), &
      INTENT(IN)                             :: array
    LOGICAL, DIMENSION(:, :, :), &
      INTENT(IN), OPTIONAL                   :: mask
    REAL(KIND=sp)                            :: ks

    INTEGER                                  :: i1, i2, i3
    REAL(KIND=sp)                            :: c, t, y

    ks = szero ; t = szero ; y = szero ; c = szero

    IF (PRESENT(mask)) THEN
      DO i3=1,SIZE(array,3)
      DO i2=1,SIZE(array,2)
      DO i1=1,SIZE(array,1)
         IF (mask(i1,i2,i3            )) THEN
           y  = array(i1,i2,i3            ) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
         END IF
      ENDDO
      ENDDO
      ENDDO
    ELSE
      DO i3=1,SIZE(array,3)
      DO i2=1,SIZE(array,2)
      DO i1=1,SIZE(array,1)
           y  = array(i1,i2,i3            ) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
      ENDDO
      ENDDO
      ENDDO
    ENDIF
  END FUNCTION kahan_sum_s3

! *****************************************************************************
!> \brief ...
!> \param array ...
!> \param mask ...
!> \retval ks ...
! *****************************************************************************
  FUNCTION kahan_sum_d3(array,mask) RESULT(ks)
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(IN)                             :: array
    LOGICAL, DIMENSION(:, :, :), &
      INTENT(IN), OPTIONAL                   :: mask
    REAL(KIND=dp)                            :: ks

    INTEGER                                  :: i1, i2, i3
    REAL(KIND=dp)                            :: c, t, y

    ks = dzero ; t = dzero ; y = dzero ; c = dzero

    IF (PRESENT(mask)) THEN
      DO i3=1,SIZE(array,3)
      DO i2=1,SIZE(array,2)
      DO i1=1,SIZE(array,1)
         IF (mask(i1,i2,i3            )) THEN
           y  = array(i1,i2,i3            ) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
         END IF
      ENDDO
      ENDDO
      ENDDO
    ELSE
      DO i3=1,SIZE(array,3)
      DO i2=1,SIZE(array,2)
      DO i1=1,SIZE(array,1)
           y  = array(i1,i2,i3            ) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
      ENDDO
      ENDDO
      ENDDO
    ENDIF
  END FUNCTION kahan_sum_d3

! *****************************************************************************
!> \brief ...
!> \param array ...
!> \param mask ...
!> \retval ks ...
! *****************************************************************************
  FUNCTION kahan_sum_c3(array,mask) RESULT(ks)
    COMPLEX(KIND=sp), DIMENSION(:, :, :), &
      INTENT(IN)                             :: array
    LOGICAL, DIMENSION(:, :, :), &
      INTENT(IN), OPTIONAL                   :: mask
    COMPLEX(KIND=sp)                         :: ks

    COMPLEX(KIND=sp)                         :: c, t, y
    INTEGER                                  :: i1, i2, i3

    ks = czero ; t = czero ; y = czero ; c = czero

    IF (PRESENT(mask)) THEN
      DO i3=1,SIZE(array,3)
      DO i2=1,SIZE(array,2)
      DO i1=1,SIZE(array,1)
         IF (mask(i1,i2,i3            )) THEN
           y  = array(i1,i2,i3            ) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
         END IF
      ENDDO
      ENDDO
      ENDDO
    ELSE
      DO i3=1,SIZE(array,3)
      DO i2=1,SIZE(array,2)
      DO i1=1,SIZE(array,1)
           y  = array(i1,i2,i3            ) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
      ENDDO
      ENDDO
      ENDDO
    ENDIF
  END FUNCTION kahan_sum_c3

! *****************************************************************************
!> \brief ...
!> \param array ...
!> \param mask ...
!> \retval ks ...
! *****************************************************************************
  FUNCTION kahan_sum_z3(array,mask) RESULT(ks)
    COMPLEX(KIND=dp), DIMENSION(:, :, :), &
      INTENT(IN)                             :: array
    LOGICAL, DIMENSION(:, :, :), &
      INTENT(IN), OPTIONAL                   :: mask
    COMPLEX(KIND=dp)                         :: ks

    COMPLEX(KIND=dp)                         :: c, t, y
    INTEGER                                  :: i1, i2, i3

    ks = zzero ; t = zzero ; y = zzero ; c = zzero

    IF (PRESENT(mask)) THEN
      DO i3=1,SIZE(array,3)
      DO i2=1,SIZE(array,2)
      DO i1=1,SIZE(array,1)
         IF (mask(i1,i2,i3            )) THEN
           y  = array(i1,i2,i3            ) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
         END IF
      ENDDO
      ENDDO
      ENDDO
    ELSE
      DO i3=1,SIZE(array,3)
      DO i2=1,SIZE(array,2)
      DO i1=1,SIZE(array,1)
           y  = array(i1,i2,i3            ) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
      ENDDO
      ENDDO
      ENDDO
    ENDIF
  END FUNCTION kahan_sum_z3

! *****************************************************************************
!> \brief ...
!> \param array ...
!> \param mask ...
!> \retval ks ...
! *****************************************************************************
  FUNCTION kahan_sum_s4(array,mask) RESULT(ks)
    REAL(KIND=sp), DIMENSION(:, :, :, :), &
      INTENT(IN)                             :: array
    LOGICAL, DIMENSION(:, :, :, :), &
      INTENT(IN), OPTIONAL                   :: mask
    REAL(KIND=sp)                            :: ks

    INTEGER                                  :: i1, i2, i3, i4
    REAL(KIND=sp)                            :: c, t, y

    ks = szero ; t = szero ; y = szero ; c = szero

    IF (PRESENT(mask)) THEN
      DO i4=1,SIZE(array,4)
      DO i3=1,SIZE(array,3)
      DO i2=1,SIZE(array,2)
      DO i1=1,SIZE(array,1)
         IF (mask(i1,i2,i3,i4         )) THEN
           y  = array(i1,i2,i3,i4         ) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
         END IF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
    ELSE
      DO i4=1,SIZE(array,4)
      DO i3=1,SIZE(array,3)
      DO i2=1,SIZE(array,2)
      DO i1=1,SIZE(array,1)
           y  = array(i1,i2,i3,i4         ) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
      ENDDO
      ENDDO
      ENDDO
      ENDDO
    ENDIF
  END FUNCTION kahan_sum_s4

! *****************************************************************************
!> \brief ...
!> \param array ...
!> \param mask ...
!> \retval ks ...
! *****************************************************************************
  FUNCTION kahan_sum_d4(array,mask) RESULT(ks)
    REAL(KIND=dp), DIMENSION(:, :, :, :), &
      INTENT(IN)                             :: array
    LOGICAL, DIMENSION(:, :, :, :), &
      INTENT(IN), OPTIONAL                   :: mask
    REAL(KIND=dp)                            :: ks

    INTEGER                                  :: i1, i2, i3, i4
    REAL(KIND=dp)                            :: c, t, y

    ks = dzero ; t = dzero ; y = dzero ; c = dzero

    IF (PRESENT(mask)) THEN
      DO i4=1,SIZE(array,4)
      DO i3=1,SIZE(array,3)
      DO i2=1,SIZE(array,2)
      DO i1=1,SIZE(array,1)
         IF (mask(i1,i2,i3,i4         )) THEN
           y  = array(i1,i2,i3,i4         ) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
         END IF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
    ELSE
      DO i4=1,SIZE(array,4)
      DO i3=1,SIZE(array,3)
      DO i2=1,SIZE(array,2)
      DO i1=1,SIZE(array,1)
           y  = array(i1,i2,i3,i4         ) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
      ENDDO
      ENDDO
      ENDDO
      ENDDO
    ENDIF
  END FUNCTION kahan_sum_d4

! *****************************************************************************
!> \brief ...
!> \param array ...
!> \param mask ...
!> \retval ks ...
! *****************************************************************************
  FUNCTION kahan_sum_c4(array,mask) RESULT(ks)
    COMPLEX(KIND=sp), &
      DIMENSION(:, :, :, :), INTENT(IN)      :: array
    LOGICAL, DIMENSION(:, :, :, :), &
      INTENT(IN), OPTIONAL                   :: mask
    COMPLEX(KIND=sp)                         :: ks

    COMPLEX(KIND=sp)                         :: c, t, y
    INTEGER                                  :: i1, i2, i3, i4

    ks = czero ; t = czero ; y = czero ; c = czero

    IF (PRESENT(mask)) THEN
      DO i4=1,SIZE(array,4)
      DO i3=1,SIZE(array,3)
      DO i2=1,SIZE(array,2)
      DO i1=1,SIZE(array,1)
         IF (mask(i1,i2,i3,i4         )) THEN
           y  = array(i1,i2,i3,i4         ) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
         END IF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
    ELSE
      DO i4=1,SIZE(array,4)
      DO i3=1,SIZE(array,3)
      DO i2=1,SIZE(array,2)
      DO i1=1,SIZE(array,1)
           y  = array(i1,i2,i3,i4         ) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
      ENDDO
      ENDDO
      ENDDO
      ENDDO
    ENDIF
  END FUNCTION kahan_sum_c4

! *****************************************************************************
!> \brief ...
!> \param array ...
!> \param mask ...
!> \retval ks ...
! *****************************************************************************
  FUNCTION kahan_sum_z4(array,mask) RESULT(ks)
    COMPLEX(KIND=dp), &
      DIMENSION(:, :, :, :), INTENT(IN)      :: array
    LOGICAL, DIMENSION(:, :, :, :), &
      INTENT(IN), OPTIONAL                   :: mask
    COMPLEX(KIND=dp)                         :: ks

    COMPLEX(KIND=dp)                         :: c, t, y
    INTEGER                                  :: i1, i2, i3, i4

    ks = zzero ; t = zzero ; y = zzero ; c = zzero

    IF (PRESENT(mask)) THEN
      DO i4=1,SIZE(array,4)
      DO i3=1,SIZE(array,3)
      DO i2=1,SIZE(array,2)
      DO i1=1,SIZE(array,1)
         IF (mask(i1,i2,i3,i4         )) THEN
           y  = array(i1,i2,i3,i4         ) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
         END IF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
    ELSE
      DO i4=1,SIZE(array,4)
      DO i3=1,SIZE(array,3)
      DO i2=1,SIZE(array,2)
      DO i1=1,SIZE(array,1)
           y  = array(i1,i2,i3,i4         ) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
      ENDDO
      ENDDO
      ENDDO
      ENDDO
    ENDIF
  END FUNCTION kahan_sum_z4

! *****************************************************************************
!> \brief ...
!> \param array ...
!> \param mask ...
!> \retval ks ...
! *****************************************************************************
  FUNCTION kahan_sum_s5(array,mask) RESULT(ks)
    REAL(KIND=sp), &
      DIMENSION(:, :, :, :, :), INTENT(IN)   :: array
    LOGICAL, DIMENSION(:, :, :, :, :), &
      INTENT(IN), OPTIONAL                   :: mask
    REAL(KIND=sp)                            :: ks

    INTEGER                                  :: i1, i2, i3, i4, i5
    REAL(KIND=sp)                            :: c, t, y

    ks = szero ; t = szero ; y = szero ; c = szero

    IF (PRESENT(mask)) THEN
      DO i5=1,SIZE(array,5)
      DO i4=1,SIZE(array,4)
      DO i3=1,SIZE(array,3)
      DO i2=1,SIZE(array,2)
      DO i1=1,SIZE(array,1)
         IF (mask(i1,i2,i3,i4,i5      )) THEN
           y  = array(i1,i2,i3,i4,i5      ) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
         END IF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
    ELSE
      DO i5=1,SIZE(array,5)
      DO i4=1,SIZE(array,4)
      DO i3=1,SIZE(array,3)
      DO i2=1,SIZE(array,2)
      DO i1=1,SIZE(array,1)
           y  = array(i1,i2,i3,i4,i5      ) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
    ENDIF
  END FUNCTION kahan_sum_s5

! *****************************************************************************
!> \brief ...
!> \param array ...
!> \param mask ...
!> \retval ks ...
! *****************************************************************************
  FUNCTION kahan_sum_d5(array,mask) RESULT(ks)
    REAL(KIND=dp), &
      DIMENSION(:, :, :, :, :), INTENT(IN)   :: array
    LOGICAL, DIMENSION(:, :, :, :, :), &
      INTENT(IN), OPTIONAL                   :: mask
    REAL(KIND=dp)                            :: ks

    INTEGER                                  :: i1, i2, i3, i4, i5
    REAL(KIND=dp)                            :: c, t, y

    ks = dzero ; t = dzero ; y = dzero ; c = dzero

    IF (PRESENT(mask)) THEN
      DO i5=1,SIZE(array,5)
      DO i4=1,SIZE(array,4)
      DO i3=1,SIZE(array,3)
      DO i2=1,SIZE(array,2)
      DO i1=1,SIZE(array,1)
         IF (mask(i1,i2,i3,i4,i5      )) THEN
           y  = array(i1,i2,i3,i4,i5      ) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
         END IF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
    ELSE
      DO i5=1,SIZE(array,5)
      DO i4=1,SIZE(array,4)
      DO i3=1,SIZE(array,3)
      DO i2=1,SIZE(array,2)
      DO i1=1,SIZE(array,1)
           y  = array(i1,i2,i3,i4,i5      ) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
    ENDIF
  END FUNCTION kahan_sum_d5

! *****************************************************************************
!> \brief ...
!> \param array ...
!> \param mask ...
!> \retval ks ...
! *****************************************************************************
  FUNCTION kahan_sum_c5(array,mask) RESULT(ks)
    COMPLEX(KIND=sp), &
      DIMENSION(:, :, :, :, :), INTENT(IN)   :: array
    LOGICAL, DIMENSION(:, :, :, :, :), &
      INTENT(IN), OPTIONAL                   :: mask
    COMPLEX(KIND=sp)                         :: ks

    COMPLEX(KIND=sp)                         :: c, t, y
    INTEGER                                  :: i1, i2, i3, i4, i5

    ks = czero ; t = czero ; y = czero ; c = czero

    IF (PRESENT(mask)) THEN
      DO i5=1,SIZE(array,5)
      DO i4=1,SIZE(array,4)
      DO i3=1,SIZE(array,3)
      DO i2=1,SIZE(array,2)
      DO i1=1,SIZE(array,1)
         IF (mask(i1,i2,i3,i4,i5      )) THEN
           y  = array(i1,i2,i3,i4,i5      ) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
         END IF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
    ELSE
      DO i5=1,SIZE(array,5)
      DO i4=1,SIZE(array,4)
      DO i3=1,SIZE(array,3)
      DO i2=1,SIZE(array,2)
      DO i1=1,SIZE(array,1)
           y  = array(i1,i2,i3,i4,i5      ) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
    ENDIF
  END FUNCTION kahan_sum_c5

! *****************************************************************************
!> \brief ...
!> \param array ...
!> \param mask ...
!> \retval ks ...
! *****************************************************************************
  FUNCTION kahan_sum_z5(array,mask) RESULT(ks)
    COMPLEX(KIND=dp), &
      DIMENSION(:, :, :, :, :), INTENT(IN)   :: array
    LOGICAL, DIMENSION(:, :, :, :, :), &
      INTENT(IN), OPTIONAL                   :: mask
    COMPLEX(KIND=dp)                         :: ks

    COMPLEX(KIND=dp)                         :: c, t, y
    INTEGER                                  :: i1, i2, i3, i4, i5

    ks = zzero ; t = zzero ; y = zzero ; c = zzero

    IF (PRESENT(mask)) THEN
      DO i5=1,SIZE(array,5)
      DO i4=1,SIZE(array,4)
      DO i3=1,SIZE(array,3)
      DO i2=1,SIZE(array,2)
      DO i1=1,SIZE(array,1)
         IF (mask(i1,i2,i3,i4,i5      )) THEN
           y  = array(i1,i2,i3,i4,i5      ) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
         END IF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
    ELSE
      DO i5=1,SIZE(array,5)
      DO i4=1,SIZE(array,4)
      DO i3=1,SIZE(array,3)
      DO i2=1,SIZE(array,2)
      DO i1=1,SIZE(array,1)
           y  = array(i1,i2,i3,i4,i5      ) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
    ENDIF
  END FUNCTION kahan_sum_z5

! *****************************************************************************
!> \brief ...
!> \param array ...
!> \param mask ...
!> \retval ks ...
! *****************************************************************************
  FUNCTION kahan_sum_s6(array,mask) RESULT(ks)
    REAL(KIND=sp), &
      DIMENSION(:, :, :, :, :, :), &
      INTENT(IN)                             :: array
    LOGICAL, DIMENSION(:, :, :, :, :, :), &
      INTENT(IN), OPTIONAL                   :: mask
    REAL(KIND=sp)                            :: ks

    INTEGER                                  :: i1, i2, i3, i4, i5, i6
    REAL(KIND=sp)                            :: c, t, y

    ks = szero ; t = szero ; y = szero ; c = szero

    IF (PRESENT(mask)) THEN
      DO i6=1,SIZE(array,6)
      DO i5=1,SIZE(array,5)
      DO i4=1,SIZE(array,4)
      DO i3=1,SIZE(array,3)
      DO i2=1,SIZE(array,2)
      DO i1=1,SIZE(array,1)
         IF (mask(i1,i2,i3,i4,i5,i6   )) THEN
           y  = array(i1,i2,i3,i4,i5,i6   ) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
         END IF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
    ELSE
      DO i6=1,SIZE(array,6)
      DO i5=1,SIZE(array,5)
      DO i4=1,SIZE(array,4)
      DO i3=1,SIZE(array,3)
      DO i2=1,SIZE(array,2)
      DO i1=1,SIZE(array,1)
           y  = array(i1,i2,i3,i4,i5,i6   ) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
    ENDIF
  END FUNCTION kahan_sum_s6

! *****************************************************************************
!> \brief ...
!> \param array ...
!> \param mask ...
!> \retval ks ...
! *****************************************************************************
  FUNCTION kahan_sum_d6(array,mask) RESULT(ks)
    REAL(KIND=dp), &
      DIMENSION(:, :, :, :, :, :), &
      INTENT(IN)                             :: array
    LOGICAL, DIMENSION(:, :, :, :, :, :), &
      INTENT(IN), OPTIONAL                   :: mask
    REAL(KIND=dp)                            :: ks

    INTEGER                                  :: i1, i2, i3, i4, i5, i6
    REAL(KIND=dp)                            :: c, t, y

    ks = dzero ; t = dzero ; y = dzero ; c = dzero

    IF (PRESENT(mask)) THEN
      DO i6=1,SIZE(array,6)
      DO i5=1,SIZE(array,5)
      DO i4=1,SIZE(array,4)
      DO i3=1,SIZE(array,3)
      DO i2=1,SIZE(array,2)
      DO i1=1,SIZE(array,1)
         IF (mask(i1,i2,i3,i4,i5,i6   )) THEN
           y  = array(i1,i2,i3,i4,i5,i6   ) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
         END IF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
    ELSE
      DO i6=1,SIZE(array,6)
      DO i5=1,SIZE(array,5)
      DO i4=1,SIZE(array,4)
      DO i3=1,SIZE(array,3)
      DO i2=1,SIZE(array,2)
      DO i1=1,SIZE(array,1)
           y  = array(i1,i2,i3,i4,i5,i6   ) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
    ENDIF
  END FUNCTION kahan_sum_d6

! *****************************************************************************
!> \brief ...
!> \param array ...
!> \param mask ...
!> \retval ks ...
! *****************************************************************************
  FUNCTION kahan_sum_c6(array,mask) RESULT(ks)
    COMPLEX(KIND=sp), &
      DIMENSION(:, :, :, :, :, :), &
      INTENT(IN)                             :: array
    LOGICAL, DIMENSION(:, :, :, :, :, :), &
      INTENT(IN), OPTIONAL                   :: mask
    COMPLEX(KIND=sp)                         :: ks

    COMPLEX(KIND=sp)                         :: c, t, y
    INTEGER                                  :: i1, i2, i3, i4, i5, i6

    ks = czero ; t = czero ; y = czero ; c = czero

    IF (PRESENT(mask)) THEN
      DO i6=1,SIZE(array,6)
      DO i5=1,SIZE(array,5)
      DO i4=1,SIZE(array,4)
      DO i3=1,SIZE(array,3)
      DO i2=1,SIZE(array,2)
      DO i1=1,SIZE(array,1)
         IF (mask(i1,i2,i3,i4,i5,i6   )) THEN
           y  = array(i1,i2,i3,i4,i5,i6   ) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
         END IF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
    ELSE
      DO i6=1,SIZE(array,6)
      DO i5=1,SIZE(array,5)
      DO i4=1,SIZE(array,4)
      DO i3=1,SIZE(array,3)
      DO i2=1,SIZE(array,2)
      DO i1=1,SIZE(array,1)
           y  = array(i1,i2,i3,i4,i5,i6   ) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
    ENDIF
  END FUNCTION kahan_sum_c6

! *****************************************************************************
!> \brief ...
!> \param array ...
!> \param mask ...
!> \retval ks ...
! *****************************************************************************
  FUNCTION kahan_sum_z6(array,mask) RESULT(ks)
    COMPLEX(KIND=dp), &
      DIMENSION(:, :, :, :, :, :), &
      INTENT(IN)                             :: array
    LOGICAL, DIMENSION(:, :, :, :, :, :), &
      INTENT(IN), OPTIONAL                   :: mask
    COMPLEX(KIND=dp)                         :: ks

    COMPLEX(KIND=dp)                         :: c, t, y
    INTEGER                                  :: i1, i2, i3, i4, i5, i6

    ks = zzero ; t = zzero ; y = zzero ; c = zzero

    IF (PRESENT(mask)) THEN
      DO i6=1,SIZE(array,6)
      DO i5=1,SIZE(array,5)
      DO i4=1,SIZE(array,4)
      DO i3=1,SIZE(array,3)
      DO i2=1,SIZE(array,2)
      DO i1=1,SIZE(array,1)
         IF (mask(i1,i2,i3,i4,i5,i6   )) THEN
           y  = array(i1,i2,i3,i4,i5,i6   ) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
         END IF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
    ELSE
      DO i6=1,SIZE(array,6)
      DO i5=1,SIZE(array,5)
      DO i4=1,SIZE(array,4)
      DO i3=1,SIZE(array,3)
      DO i2=1,SIZE(array,2)
      DO i1=1,SIZE(array,1)
           y  = array(i1,i2,i3,i4,i5,i6   ) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
    ENDIF
  END FUNCTION kahan_sum_z6

! *****************************************************************************
!> \brief ...
!> \param array ...
!> \param mask ...
!> \retval ks ...
! *****************************************************************************
  FUNCTION kahan_sum_s7(array,mask) RESULT(ks)
    REAL(KIND=sp), &
      DIMENSION(:, :, :, :, :, :, :), &
      INTENT(IN)                             :: array
    LOGICAL, &
      DIMENSION(:, :, :, :, :, :, :), &
      INTENT(IN), OPTIONAL                   :: mask
    REAL(KIND=sp)                            :: ks

    INTEGER                                  :: i1, i2, i3, i4, i5, i6, i7
    REAL(KIND=sp)                            :: c, t, y

    ks = szero ; t = szero ; y = szero ; c = szero

    IF (PRESENT(mask)) THEN
      DO i7=1,SIZE(array,7)
      DO i6=1,SIZE(array,6)
      DO i5=1,SIZE(array,5)
      DO i4=1,SIZE(array,4)
      DO i3=1,SIZE(array,3)
      DO i2=1,SIZE(array,2)
      DO i1=1,SIZE(array,1)
         IF (mask(i1,i2,i3,i4,i5,i6,i7)) THEN
           y  = array(i1,i2,i3,i4,i5,i6,i7) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
         END IF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
    ELSE
      DO i7=1,SIZE(array,7)
      DO i6=1,SIZE(array,6)
      DO i5=1,SIZE(array,5)
      DO i4=1,SIZE(array,4)
      DO i3=1,SIZE(array,3)
      DO i2=1,SIZE(array,2)
      DO i1=1,SIZE(array,1)
           y  = array(i1,i2,i3,i4,i5,i6,i7) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
    ENDIF
  END FUNCTION kahan_sum_s7

! *****************************************************************************
!> \brief ...
!> \param array ...
!> \param mask ...
!> \retval ks ...
! *****************************************************************************
  FUNCTION kahan_sum_d7(array,mask) RESULT(ks)
    REAL(KIND=dp), &
      DIMENSION(:, :, :, :, :, :, :), &
      INTENT(IN)                             :: array
    LOGICAL, &
      DIMENSION(:, :, :, :, :, :, :), &
      INTENT(IN), OPTIONAL                   :: mask
    REAL(KIND=dp)                            :: ks

    INTEGER                                  :: i1, i2, i3, i4, i5, i6, i7
    REAL(KIND=dp)                            :: c, t, y

    ks = dzero ; t = dzero ; y = dzero ; c = dzero

    IF (PRESENT(mask)) THEN
      DO i7=1,SIZE(array,7)
      DO i6=1,SIZE(array,6)
      DO i5=1,SIZE(array,5)
      DO i4=1,SIZE(array,4)
      DO i3=1,SIZE(array,3)
      DO i2=1,SIZE(array,2)
      DO i1=1,SIZE(array,1)
         IF (mask(i1,i2,i3,i4,i5,i6,i7)) THEN
           y  = array(i1,i2,i3,i4,i5,i6,i7) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
         END IF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
    ELSE
      DO i7=1,SIZE(array,7)
      DO i6=1,SIZE(array,6)
      DO i5=1,SIZE(array,5)
      DO i4=1,SIZE(array,4)
      DO i3=1,SIZE(array,3)
      DO i2=1,SIZE(array,2)
      DO i1=1,SIZE(array,1)
           y  = array(i1,i2,i3,i4,i5,i6,i7) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
    ENDIF
  END FUNCTION kahan_sum_d7

! *****************************************************************************
!> \brief ...
!> \param array ...
!> \param mask ...
!> \retval ks ...
! *****************************************************************************
  FUNCTION kahan_sum_c7(array,mask) RESULT(ks)
    COMPLEX(KIND=sp), &
      DIMENSION(:, :, :, :, :, :, :), &
      INTENT(IN)                             :: array
    LOGICAL, &
      DIMENSION(:, :, :, :, :, :, :), &
      INTENT(IN), OPTIONAL                   :: mask
    COMPLEX(KIND=sp)                         :: ks

    COMPLEX(KIND=sp)                         :: c, t, y
    INTEGER                                  :: i1, i2, i3, i4, i5, i6, i7

    ks = czero ; t = czero ; y = czero ; c = czero

    IF (PRESENT(mask)) THEN
      DO i7=1,SIZE(array,7)
      DO i6=1,SIZE(array,6)
      DO i5=1,SIZE(array,5)
      DO i4=1,SIZE(array,4)
      DO i3=1,SIZE(array,3)
      DO i2=1,SIZE(array,2)
      DO i1=1,SIZE(array,1)
         IF (mask(i1,i2,i3,i4,i5,i6,i7)) THEN
           y  = array(i1,i2,i3,i4,i5,i6,i7) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
         END IF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
    ELSE
      DO i7=1,SIZE(array,7)
      DO i6=1,SIZE(array,6)
      DO i5=1,SIZE(array,5)
      DO i4=1,SIZE(array,4)
      DO i3=1,SIZE(array,3)
      DO i2=1,SIZE(array,2)
      DO i1=1,SIZE(array,1)
           y  = array(i1,i2,i3,i4,i5,i6,i7) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
    ENDIF
  END FUNCTION kahan_sum_c7

! *****************************************************************************
!> \brief ...
!> \param array ...
!> \param mask ...
!> \retval ks ...
! *****************************************************************************
  FUNCTION kahan_sum_z7(array,mask) RESULT(ks)
    COMPLEX(KIND=dp), &
      DIMENSION(:, :, :, :, :, :, :), &
      INTENT(IN)                             :: array
    LOGICAL, &
      DIMENSION(:, :, :, :, :, :, :), &
      INTENT(IN), OPTIONAL                   :: mask
    COMPLEX(KIND=dp)                         :: ks

    COMPLEX(KIND=dp)                         :: c, t, y
    INTEGER                                  :: i1, i2, i3, i4, i5, i6, i7

    ks = zzero ; t = zzero ; y = zzero ; c = zzero

    IF (PRESENT(mask)) THEN
      DO i7=1,SIZE(array,7)
      DO i6=1,SIZE(array,6)
      DO i5=1,SIZE(array,5)
      DO i4=1,SIZE(array,4)
      DO i3=1,SIZE(array,3)
      DO i2=1,SIZE(array,2)
      DO i1=1,SIZE(array,1)
         IF (mask(i1,i2,i3,i4,i5,i6,i7)) THEN
           y  = array(i1,i2,i3,i4,i5,i6,i7) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
         END IF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
    ELSE
      DO i7=1,SIZE(array,7)
      DO i6=1,SIZE(array,6)
      DO i5=1,SIZE(array,5)
      DO i4=1,SIZE(array,4)
      DO i3=1,SIZE(array,3)
      DO i2=1,SIZE(array,2)
      DO i1=1,SIZE(array,1)
           y  = array(i1,i2,i3,i4,i5,i6,i7) - c
           t  = ks + y
           c  = (t - ks) - y
           ks = t
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
    ENDIF
  END FUNCTION kahan_sum_z7

END MODULE kahan_sum
