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

!> @file tblite/scf/mixer/cache.f90
!> Provides a base cache definition for the electronic mixer

!> Base definition of a data container for mutable mixer data
module tblite_scf_mixer_cache
   use mctc_env, only : wp
   implicit none
   private

   !> Container for mixer caches
   type, public :: mixer_cache_container
      !> Actual mixer cache
      class(mixer_cache), allocatable :: raw
   end type mixer_cache_container

   !> Abstract base type for the mixer cache
   type, public, abstract :: mixer_cache
      !> Flag to indicate initial seeding of the mixer cache
      logical :: initialized = .false.

      !> Size of the extrapolated vector
      integer :: ndim = 0
      !> Iteration counter for this mixer
      integer :: iter = 0
      !> Index for setting the input vector
      integer :: iset = 0
      !> Index for setting the residuum vector
      integer :: idif = 0
      !> Index for getting the extrapolated vector
      integer :: iget = 0

      !> Current input vector (last extrapolated vector): [ndim]
      real(wp), allocatable :: x_in(:)
      !> Current residual vector of new iteration minus x_in: [ndim]
      real(wp), allocatable :: r_in(:)
   end type mixer_cache

end module tblite_scf_mixer_cache