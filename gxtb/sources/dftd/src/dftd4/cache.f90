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

!> Cache for precomputed dispersion model properties and coefficients
module dftd4_cache
   use mctc_env, only : wp
   implicit none
   private

   public :: dispersion_cache


   !> Cache for dispersion model data: C6 coefficients, derivatives, and
   !> dynamic polarizabilities
   type :: dispersion_cache

      !> Precomputed C6 coefficients: [nat, nat]
      real(wp), allocatable :: c6(:, :)

      !> Derivative of C6 w.r.t. coordination number: [nat, nat]
      real(wp), allocatable :: dc6dcn(:, :)

      !> Derivative of C6 w.r.t. partial charge: [nat, nat]
      real(wp), allocatable :: dc6dq(:, :)

      !> Dipole-dipole dynamic polarizabilities: [ngrid, nat]
      real(wp), allocatable :: adiw(:, :)

      !> Quadrupole-quadrupole dynamic polarizabilities: [ngrid, nat]
      real(wp), allocatable :: aqiw(:, :)

      !> Derivative of dipole-dipole dynamic polarizabilities
      !> w.r.t. the coordination numbers: [ngrid, nat]
      real(wp), allocatable :: dadiwdcn(:, :)

      !> Derivative of dipole-dipole dynamic polarizabilities
      !> w.r.t. the atomic partial charges: [ngrid, nat]
      real(wp), allocatable :: dadiwdq(:, :)

      !> Derivative of quadrupole-quadrupole dynamic polarizabilities
      !> w.r.t. the coordination numbers: [ngrid, nat]
      real(wp), allocatable :: daqiwdcn(:, :)

      !> Derivative of quadrupole-quadrupole dynamic polarizabilities
      !> w.r.t. the atomic partial charge: [ngrid, nat]
      real(wp), allocatable :: daqiwdq(:, :)

   end type dispersion_cache


end module dftd4_cache
