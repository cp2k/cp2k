!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

MODULE trajana_fft
   USE trajana_kinds,                   ONLY: dp

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: fft_any_in_place, fft_in_place, next_power_of_two

CONTAINS

   RECURSIVE SUBROUTINE fft_any_in_place(DATA, inverse)
      COMPLEX(dp), INTENT(INOUT)                         :: data(:)
      LOGICAL, INTENT(IN)                                :: inverse

      COMPLEX(dp), ALLOCATABLE                           :: first(:), second(:)
      COMPLEX(dp)                                        :: phase
      INTEGER                                            :: index, n, nfft
      REAL(dp)                                           :: angle, pi

      n = SIZE(DATA)
      IF (n <= 1) RETURN
      IF (inverse) THEN
         DATA = CONJG(DATA)
         CALL fft_any_in_place(DATA, inverse=.FALSE.)
         DATA = CONJG(DATA)/REAL(n, dp)
         RETURN
      END IF
      IF (IAND(n, n - 1) == 0) THEN
         CALL fft_in_place(DATA, inverse=.FALSE.)
         RETURN
      END IF

      nfft = next_power_of_two(2*n - 1)
      ALLOCATE (first(nfft), second(nfft))
      first = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      second = CMPLX(0.0_dp, 0.0_dp, KIND=dp)
      pi = ACOS(-1.0_dp)
      DO index = 0, n - 1
         angle = pi*REAL(index, dp)*REAL(index, dp)/REAL(n, dp)
         phase = CMPLX(COS(angle), -SIN(angle), KIND=dp)
         first(index + 1) = DATA(index + 1)*phase
         phase = CONJG(phase)
         second(index + 1) = phase
         IF (index > 0) second(nfft - index + 1) = phase
      END DO
      CALL fft_in_place(first, inverse=.FALSE.)
      CALL fft_in_place(second, inverse=.FALSE.)
      first = first*second
      CALL fft_in_place(first, inverse=.TRUE.)
      DO index = 0, n - 1
         angle = pi*REAL(index, dp)*REAL(index, dp)/REAL(n, dp)
         phase = CMPLX(COS(angle), -SIN(angle), KIND=dp)
         DATA(index + 1) = first(index + 1)*phase
      END DO
   END SUBROUTINE fft_any_in_place

   PURE INTEGER FUNCTION next_power_of_two(value)
      INTEGER, INTENT(IN)                                :: value

      next_power_of_two = 1
      DO WHILE (next_power_of_two < value)
         next_power_of_two = 2*next_power_of_two
      END DO
   END FUNCTION next_power_of_two

   SUBROUTINE fft_in_place(DATA, inverse)
      COMPLEX(dp), INTENT(INOUT)                         :: data(:)
      LOGICAL, INTENT(IN)                                :: inverse

      COMPLEX(dp)                                        :: even, factor, odd, root, temporary
      INTEGER                                            :: bit, half, i, j, length, n, offset
      REAL(dp)                                           :: angle, pi

      n = SIZE(DATA)
      IF (n <= 1) RETURN
      IF (IAND(n, n - 1) /= 0) ERROR STOP "FFT length is not a power of two"

      j = 1
      DO i = 2, n
         bit = n/2
         DO WHILE (j > bit)
            j = j - bit
            bit = bit/2
         END DO
         j = j + bit
         IF (i < j) THEN
            temporary = DATA(i)
            DATA(i) = DATA(j)
            DATA(j) = temporary
         END IF
      END DO

      pi = ACOS(-1.0_dp)
      length = 2
      DO WHILE (length <= n)
         angle = 2.0_dp*pi/REAL(length, dp)
         IF (.NOT. inverse) angle = -angle
         root = CMPLX(COS(angle), SIN(angle), KIND=dp)
         half = length/2
         DO offset = 1, n, length
            factor = CMPLX(1.0_dp, 0.0_dp, KIND=dp)
            DO j = 0, half - 1
               even = DATA(offset + j)
               odd = factor*DATA(offset + j + half)
               DATA(offset + j) = even + odd
               DATA(offset + j + half) = even - odd
               factor = factor*root
            END DO
         END DO
         length = 2*length
      END DO
      IF (inverse) DATA = DATA/REAL(n, dp)
   END SUBROUTINE fft_in_place

END MODULE trajana_fft
