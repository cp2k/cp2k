! This file is part of dftd4.
! SPDX-Identifier: LGPL-3.0-or-later
!
! dftd4 is free software: you can redistribute it and/or modify it under
! the terms of the Lesser GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! dftd4 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! Lesser GNU General Public License for more details.
!
! You should have received a copy of the Lesser GNU General Public License
! along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

!> Utility functions for the dispersion models 
module dftd4_model_utils
   use ieee_arithmetic, only : ieee_is_nan
   use mctc_env, only : wp
   implicit none
   
   public :: is_exceptional, weight_cn, zeta, dzeta, zeta_tanh, dzeta_tanh

contains


!> Check whether we are dealing with an exceptional value, NaN or Inf
elemental function is_exceptional(val)
   real(wp), intent(in) :: val
   logical :: is_exceptional
   is_exceptional = ieee_is_nan(val) .or. abs(val) > huge(val)
end function is_exceptional


elemental function weight_cn(wf,cn,cnref) result(cngw)
   real(wp),intent(in) :: wf, cn, cnref
   real(wp) :: cngw
   intrinsic :: exp
   cngw = exp ( -wf * ( cn - cnref )**2 )
end function weight_cn

!> charge scaling function
elemental function zeta(a, c, qref, qmod)
   real(wp), intent(in) :: a
   real(wp), intent(in) :: c
   real(wp), intent(in) :: qref
   real(wp), intent(in) :: qmod
   real(wp) :: zeta

   intrinsic :: exp

   if (qmod < 0.0_wp) then
      zeta = exp( a )
   else
      zeta = exp( a * ( 1.0_wp - exp( c * ( 1.0_wp - qref/qmod ) ) ) )
   endif

end function zeta

!> derivative of charge scaling function w.r.t. charge
elemental function dzeta(a, c, qref, qmod)
   real(wp), intent(in) :: a
   real(wp), intent(in) :: c
   real(wp), intent(in) :: qref
   real(wp), intent(in) :: qmod
   real(wp) :: dzeta

   intrinsic :: exp

   if (qmod < 0.0_wp) then
      dzeta = 0.0_wp
   else
      dzeta = - a * c * exp( c * ( 1.0_wp - qref/qmod ) ) &
         & * zeta(a,c,qref,qmod) * qref / ( qmod**2 )
   endif

end function dzeta


!> Charge scaling function based on tanh interpolation
elemental function zeta_tanh(a, b, c, d, qref, qmod) result(zeta)
   !> Parameters from the tanh interpolation
   real(wp), intent(in) :: a
   real(wp), intent(in) :: b
   real(wp), intent(in) :: c
   real(wp), intent(in) :: d
   !> Charge in the reference system
   real(wp), intent(in) :: qref
   !> Charge in the actual system
   real(wp), intent(in) :: qmod
   real(wp) :: zeta

   real(wp) :: zeta_mod, zeta_ref

   intrinsic :: dtanh

   zeta_mod = a + b * dtanh( c * (qmod) + d)
   zeta_ref = a + b * dtanh( c * (qref) + d) 
   zeta = zeta_mod/zeta_ref

end function zeta_tanh


!> New charge scaling function derivative 
elemental function dzeta_tanh(a, b, c, d, qref, qmod) result(dzeta)
   !> Parameters from the tanh interpolation
   real(wp), intent(in) :: a
   real(wp), intent(in) :: b
   real(wp), intent(in) :: c
   real(wp), intent(in) :: d
   !> Charge in the reference system
   real(wp), intent(in) :: qref
   !> Charge in the actual system
   real(wp), intent(in) :: qmod
   real(wp) :: dzeta

   real(wp) ::dzeta_mod, zeta_ref

   intrinsic :: dtanh, dcosh

   dzeta_mod = c * b * (1.0_wp/ dcosh( c * (qmod) + d))**2
   zeta_ref = a + b * dtanh( c * (qref) + d) 
   dzeta =  dzeta_mod / zeta_ref 
   
end function dzeta_tanh

end module dftd4_model_utils
