! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later
!
! tblite is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! tblite is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with tblite.  If not, see <https://www.gnu.org/licenses/>.

!> @file tblite/wavefunction/spin.f90
!> Provides conversion routines to change the representation of spin-polarized densities

!> Handling of spin information in the wavefunction
!>
!> Spin is represented as charge and magnetization density in the population based
!> properties, *e.g.* Mulliken partial charges, atomic dipole moments, ..., and
!> in up/down representation for orbials energies, occupation numbers, ...
module tblite_wavefunction_spin
   use mctc_env, only : wp
   implicit none
   private

   public :: magnet_to_updown, updown_to_magnet
   public :: pot_magnet_to_updown, pot_updown_to_magnet

   !> Convert a charge-magnetization representation
   !> into an up-down representation
   !> Applies to extensive (inverse to intensive) properties (density, charge, ...)
   interface magnet_to_updown
      module procedure :: magnet_to_updown_1
      module procedure :: magnet_to_updown_2
      module procedure :: magnet_to_updown_3
      module procedure :: magnet_to_updown_4
   end interface magnet_to_updown

   !> Convert an up-down representation
   !> into a charge-magnetization representation
   !> Applies to extensive (inverse to intensive) properties (density, charge, ...)
   interface updown_to_magnet
      module procedure :: updown_to_magnet_1
      module procedure :: updown_to_magnet_2
      module procedure :: updown_to_magnet_3
      module procedure :: updown_to_magnet_4
   end interface updown_to_magnet


   !> Convert a scalar-magnetic potential representation
   !> into an up-down potential representation
   !> Applies to intensive (inverse to extensive) properties (potential, ...)
   interface pot_magnet_to_updown
      module procedure :: updown_to_magnet_1
      module procedure :: updown_to_magnet_2
      module procedure :: updown_to_magnet_3
      module procedure :: updown_to_magnet_4
   end interface pot_magnet_to_updown

   !> Convert an up-down potential representation
   !> into a scalar-magnetic potential representation
   !> Applies to intensive (inverse to extensive) properties (potential, ...)
   interface pot_updown_to_magnet
      module procedure :: magnet_to_updown_1
      module procedure :: magnet_to_updown_2
      module procedure :: magnet_to_updown_3
      module procedure :: magnet_to_updown_4
   end interface pot_updown_to_magnet

contains

!> Convert a charge-magnetization representation
!> into an up-down representation (for extensive properties)
!> Convert an up-down potential representation
!> into a scalar-magnetic potential representation (for intensive properties)
subroutine magnet_to_updown_1(x)

   !> Data in charge-magnetization representation on entry
   !> transformed to up-down representation on exit
   real(wp), intent(inout) :: x(:)

   if (size(x, 1) == 2) then
      x(1) = 0.5_wp * (x(1) + x(2))
      x(2) = x(1) - x(2)
   end if
end subroutine magnet_to_updown_1

!> Convert a charge-magnetization representation
!> into an up-down representation (for extensive properties)
!> Convert an up-down potential representation
!> into a scalar-magnetic potential representation (for intensive properties)
subroutine magnet_to_updown_2(x)

   !> Data in charge-magnetization representation on entry
   !> transformed to up-down representation on exit
   real(wp), intent(inout) :: x(:, :)

   if (size(x, 2) == 2) then
      x(:, 1) = 0.5_wp * (x(:, 1) + x(:, 2))
      x(:, 2) = x(:, 1) - x(:, 2)
   end if
end subroutine magnet_to_updown_2

!> Convert a charge-magnetization representation
!> into an up-down representation (for extensive properties)
!> Convert an up-down potential representation
!> into a scalar-magnetic potential representation (for intensive properties)
subroutine magnet_to_updown_3(x)

   !> Data in charge-magnetization representation on entry
   !> transformed to up-down representation on exit
   real(wp), intent(inout) :: x(:, :, :)

   if (size(x, 3) == 2) then
      x(:, :, 1) = 0.5_wp * (x(:, :, 1) + x(:, :, 2))
      x(:, :, 2) = x(:, :, 1) - x(:, :, 2)
   end if
end subroutine magnet_to_updown_3

!> Convert a charge-magnetization representation
!> into an up-down representation (for extensive properties)
!> Convert an up-down potential representation
!> into a scalar-magnetic potential representation (for intensive properties)
subroutine magnet_to_updown_4(x)

   !> Data in charge-magnetization representation on entry
   !> transformed to up-down representation on exit
   real(wp), intent(inout) :: x(:, :, :, :)

   if (size(x, 4) == 2) then
      x(:, :, :, 1) = 0.5_wp * (x(:, :, :, 1) + x(:, :, :, 2))
      x(:, :, :, 2) = x(:, :, :, 1) - x(:, :, :, 2)
   end if
end subroutine magnet_to_updown_4

!> Convert an up-down representation
!> into a charge-magnetization representation
!> Convert a scalar-magnetic potential representation
!> into an up-down potential representation
subroutine updown_to_magnet_1(x)

   !> Data in up-down representation on entry
   !> transformed to charge-magnetization representation on exit
   real(wp), intent(inout) :: x(:)

   if (size(x, 1) == 2) then
      x(1) = x(1) + x(2)
      x(2) = x(1) - 2.0_wp * x(2)
   end if
end subroutine updown_to_magnet_1

!> Convert an up-down representation
!> into a charge-magnetization representation
!> Convert a scalar-magnetic potential representation
!> into an up-down potential representation
subroutine updown_to_magnet_2(x)

   !> Data in up-down representation on entry
   !> transformed to charge-magnetization representation on exit
   real(wp), intent(inout) :: x(:, :)

   if (size(x, 2) == 2) then
      x(:, 1) = x(:, 1) + x(:, 2)
      x(:, 2) = x(:, 1) - 2.0_wp * x(:, 2)
   end if
end subroutine updown_to_magnet_2

!> Convert an up-down representation
!> into a charge-magnetization representation
!> Convert a scalar-magnetic potential representation
!> into an up-down potential representation
subroutine updown_to_magnet_3(x)

   !> Data in up-down representation on entry
   !> transformed to charge-magnetization representation on exit
   real(wp), intent(inout) :: x(:, :, :)

   if (size(x, 3) == 2) then
      x(:, :, 1) = x(:, :, 1) + x(:, :, 2)
      x(:, :, 2) = x(:, :, 1) - 2.0_wp * x(:, :, 2)
   end if
end subroutine updown_to_magnet_3

!> Convert an up-down representation
!> into a charge-magnetization representation
!> Convert a scalar-magnetic potential representation
!> into an up-down potential representation
subroutine updown_to_magnet_4(x)

   !> Data in up-down representation on entry
   !> transformed to charge-magnetization representation on exit
   real(wp), intent(inout) :: x(:, :, :, :)

   if (size(x, 4) == 2) then
      x(:, :, :, 1) = x(:, :, :, 1) + x(:, :, :, 2)
      x(:, :, :, 2) = x(:, :, :, 1) - 2.0_wp * x(:, :, :, 2)
   end if
end subroutine updown_to_magnet_4

end module tblite_wavefunction_spin
