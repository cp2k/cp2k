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

!> @file tblite/scf/cache.f90
!> Provides a cache for the SCF iterator

!> Data container for mutable iterator data
module tblite_scf_cache
   use mctc_env, only : wp
   use tblite_scf_info, only : scf_info
   use tblite_scf_mixer_cache, only : mixer_cache_container
   implicit none
   private

   type, public :: iterator_cache
      !> Information about the used SCF contributions
      type(scf_info) :: info
      !> Current SCF iteration
      integer :: iscf = 0
      !> Index of currently used mixer
      integer :: current_mixer = 1
      !> Index of precollecting mixer
      integer :: precollect_mixer = -1
      !> Previous SCF energy
      real(wp) :: elast = 0.0_wp
      !> Container for mixer caches
      type(mixer_cache_container), allocatable :: mcache(:)
   end type iterator_cache

end module tblite_scf_cache