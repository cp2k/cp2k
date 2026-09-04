!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Interface to libwignernj for the angular momentum algebra
!> \note
!>      libwignernj carries out all intermediate arithmetic exactly, in a prime
!>      factorization representation, and rounds only once at the end. The
!>      returned coefficients are therefore correct to the last bit, whereas the
!>      factorial-ratio expressions that CP2K used before lost roughly one digit
!>      per two units of l.
!>
!>      The library takes all angular momenta as twice their value so that
!>      half-integers stay exact; the wrappers here take the plain integer l and
!>      m of a spherical harmonic and do the doubling.
!>
!>      The real spherical harmonics S(l,m) follow the Condon-Shortley
!>      convention
!>          S(l, 0)   = Y(l,0)
!>          S(l, m>0) = ( Y(l,-m) + (-1)^m Y(l,m) )/sqrt(2)
!>          S(l, m<0) = ( Y(l,m) - (-1)^|m| Y(l,-m) )*i/sqrt(2)
!>      which is the same convention as spherical_harmonics::y_lm, so no
!>      additional phase is needed on either side.
!> \par History
!>      09.2026 created [Susi Lehtola]
!> \author Susi Lehtola
! **************************************************************************************************
MODULE wignernj_interface

   USE ISO_C_BINDING,                   ONLY: C_DOUBLE,&
                                              C_INT
   USE kinds,                           ONLY: dp
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: real_gaunt, complex_gaunt

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'wignernj_interface'

   INTERFACE
! **************************************************************************************************
!> \brief Gaunt coefficient of three complex spherical harmonics
!> \param tl1 twice the angular momentum of the first harmonic
!> \param tm1 twice the magnetic quantum number of the first harmonic
!> \param tl2 twice the angular momentum of the second harmonic
!> \param tm2 twice the magnetic quantum number of the second harmonic
!> \param tl3 twice the angular momentum of the third harmonic
!> \param tm3 twice the magnetic quantum number of the third harmonic
!> \return integral of Y(l1,m1)*Y(l2,m2)*Y(l3,m3) over the unit sphere
! **************************************************************************************************
      FUNCTION gaunt_c(tl1, tm1, tl2, tm2, tl3, tm3) BIND(C, name="gaunt")
         IMPORT :: C_DOUBLE, C_INT
         INTEGER(KIND=C_INT), VALUE                      :: tl1, tm1, tl2, tm2, tl3, tm3
         REAL(KIND=C_DOUBLE)                             :: gaunt_c
      END FUNCTION gaunt_c

! **************************************************************************************************
!> \brief Gaunt coefficient of three real spherical harmonics
!> \param tl1 twice the angular momentum of the first harmonic
!> \param tm1 twice the magnetic quantum number of the first harmonic
!> \param tl2 twice the angular momentum of the second harmonic
!> \param tm2 twice the magnetic quantum number of the second harmonic
!> \param tl3 twice the angular momentum of the third harmonic
!> \param tm3 twice the magnetic quantum number of the third harmonic
!> \return integral of S(l1,m1)*S(l2,m2)*S(l3,m3) over the unit sphere
! **************************************************************************************************
      FUNCTION gaunt_real_c(tl1, tm1, tl2, tm2, tl3, tm3) BIND(C, name="gaunt_real")
         IMPORT :: C_DOUBLE, C_INT
         INTEGER(KIND=C_INT), VALUE                      :: tl1, tm1, tl2, tm2, tl3, tm3
         REAL(KIND=C_DOUBLE)                             :: gaunt_real_c
      END FUNCTION gaunt_real_c
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief Gaunt coefficient of three real spherical harmonics,
!>        int S(l1,m1) S(l2,m2) S(l3,m3) dOmega
!> \param l1 angular momentum of the first harmonic
!> \param m1 magnetic quantum number of the first harmonic
!> \param l2 angular momentum of the second harmonic
!> \param m2 magnetic quantum number of the second harmonic
!> \param l3 angular momentum of the third harmonic
!> \param m3 magnetic quantum number of the third harmonic
!> \return the Gaunt coefficient, zero if the selection rules are violated
! **************************************************************************************************
   FUNCTION real_gaunt(l1, m1, l2, m2, l3, m3) RESULT(g)
      INTEGER, INTENT(IN)                                :: l1, m1, l2, m2, l3, m3
      REAL(KIND=dp)                                      :: g

      g = REAL(gaunt_real_c(2*l1, 2*m1, 2*l2, 2*m2, 2*l3, 2*m3), KIND=dp)

   END FUNCTION real_gaunt

! **************************************************************************************************
!> \brief Gaunt coefficient of three complex spherical harmonics,
!>        int Y(l1,m1) Y(l2,m2) Y(l3,m3) dOmega
!> \param l1 angular momentum of the first harmonic
!> \param m1 magnetic quantum number of the first harmonic
!> \param l2 angular momentum of the second harmonic
!> \param m2 magnetic quantum number of the second harmonic
!> \param l3 angular momentum of the third harmonic
!> \param m3 magnetic quantum number of the third harmonic
!> \return the Gaunt coefficient, zero if the selection rules are violated
! **************************************************************************************************
   FUNCTION complex_gaunt(l1, m1, l2, m2, l3, m3) RESULT(g)
      INTEGER, INTENT(IN)                                :: l1, m1, l2, m2, l3, m3
      REAL(KIND=dp)                                      :: g

      g = REAL(gaunt_c(2*l1, 2*m1, 2*l2, 2*m2, 2*l3, 2*m3), KIND=dp)

   END FUNCTION complex_gaunt

END MODULE wignernj_interface
