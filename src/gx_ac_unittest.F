!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief GreenX Analytic continuation unit test
!> \author Stepan Marek
! **************************************************************************************************
PROGRAM gx_ac_unittest
#include "base/base_uses.f90"
#if !defined(__GREENX)
   ! Abort and inform that GreenX was not included in the compilation
   ! Ideally, this will be avoided in testing by the conditional tests
   CPABORT("CP2K not compiled with GreenX link - please recompile to access GXAC.")
#else
   USE kinds, ONLY: dp
   USE gx_ac, ONLY: create_thiele_pade, &
                    evaluate_thiele_pade_at, &
                    free_params, &
                    params

   ! Create the dataset containing the fitting data
   ! Two Lorentzian peaks with some overlap
   COMPLEX(kind=dp)                   :: damp_one = (2, 0), &
                                         damp_two = (4, 0), &
                                         center_one = (2, 0), &
                                         center_two = (8, 0), &
                                         amp_one = (10, 0), &
                                         amp_two = (10, 0), &
                                         min_source = (0, 0), &
                                         min_fit = (0, 0), &
                                         max_source = (10, 0), &
                                         max_fit = (10, 0)
   INTEGER, PARAMETER                 :: n_source = 20, &
                                         n_fit = 100, &
                                         n_param = 10
   INTEGER                            :: i
   COMPLEX(kind=dp)                   :: d_source, &
                                         d_fit
   COMPLEX(kind=dp), DIMENSION(n_source) :: x_source, &
                                            y_source
   COMPLEX(kind=dp), DIMENSION(n_fit)    :: x_fit, &
                                            y_fit
   TYPE(params)                       :: fit_params

   d_source = (max_source - min_source)/CMPLX(n_source - 1, kind=dp)
   d_fit = (max_fit - min_fit)/CMPLX(n_fit - 1, kind=dp)

   PRINT '(A12)', "#Source data"

   DO i = 1, n_source
      x_source(i) = min_source + CMPLX(i - 1, 0.0, kind=dp)*d_source
      y_source(i) = amp_one/(damp_one*damp_one + (x_source(i) - center_one)*(x_source(i) - center_one))
      y_source(i) = y_source(i) + amp_two/(damp_two*damp_two + (x_source(i) - center_two)*(x_source(i) - center_two))
      PRINT '(E20.8E3,E20.8E3)', REAL(x_source(i), kind=dp), REAL(y_source(i), kind=dp)
   END DO

   PRINT '(A9)', "#Fit data"

   ! Fit points created
   ! Now do the actual fitting
   fit_params = create_thiele_pade(n_param, x_source, y_source)

   ! Create the evaluation grid
   DO i = 1, n_fit
      x_fit(i) = min_fit + d_fit*CMPLX(i - 1, 0, kind=dp)
   END DO

   y_fit(1:n_fit) = evaluate_thiele_pade_at(fit_params, x_fit)

   DO i = 1, n_fit
      PRINT '(E20.8E3,E20.8E3)', REAL(x_fit(i), kind=dp), REAL(y_fit(i), kind=dp)
   END DO

   CALL free_params(fit_params)
#endif
END PROGRAM gx_ac_unittest
