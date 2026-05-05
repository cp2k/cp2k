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

!> @file tblite/scf/mixer/type.f90
!> Base class for electronic mixing

!> Provides base class for electronic mixing routines
module tblite_scf_mixer_type
   use mctc_env, only : error_type, wp
   use tblite_scf_mixer_cache, only : mixer_cache_container
   implicit none
   private

   !> Container for electronic mixers
   type, public :: mixer_container
      !> Actual mixer
      class(mixer_type), allocatable :: raw
   end type mixer_container

   !> Abstract base class for electronic mixing
   type, public, abstract :: mixer_type
      !> Starting iteration for mixer
      integer :: start
      !> Convergence cutoffs to start mixer
      real(wp) :: conv_start
      !> Iterations to precollect before starting
      integer :: precollect
      !> Length of history considered for extrapolation
      integer :: memory
      !> Mode of mixing either density or potential (Fock) based
      integer :: mode
      !> Linear mixing or damping parameter
      real(wp) :: damp
   contains
      !> Update mixer cache      
      procedure(update), deferred :: update
      !> Collect data for mixing
      procedure(collect), deferred :: collect
      !> Apply mixing to the density/potential
      procedure(next), deferred :: next
      !> Set new vector
      generic :: set => set_1d, set_2d, set_3d
      !> Set new vector from 1D array
      procedure :: set_1d
      !> Set new vector from 2D array
      procedure :: set_2d
      !> Set new vector from 3D array
      procedure :: set_3d
      !> Set difference between new and old vector
      generic :: diff => diff_1d, diff_2d, diff_3d
      !> Set difference between new and old vector from 1D array
      procedure :: diff_1d
      !> Set difference between new and old vector from 2D array
      procedure :: diff_2d
      !> Set difference between new and old vector from 3D array
      procedure :: diff_3d
      !> Get extrapolated vector
      generic :: get => get_1d, get_2d, get_3d
      !> Get extrapolated vector as 1D array
      procedure :: get_1d
      !> Get extrapolated vector as 2D array
      procedure :: get_2d
      !> Get extrapolated vector as 3D array
      procedure :: get_3d
      !> Get error metric from mixing
      procedure :: get_error
      !> Remove cached mixer data
      procedure(cleanup), deferred :: cleanup
   end type mixer_type

   abstract interface
      !> Update mixer cache
      subroutine update(self, cache, ndim)
         import :: mixer_type, mixer_cache_container
         !> Instance of the electronic mixer
         class(mixer_type), intent(in) :: self
         !> Cache container for mutable mixer data
         type(mixer_cache_container), intent(inout) :: cache
         !> Size of the extrapolated vector
         integer, intent(in) :: ndim
      end subroutine update

      !> Collect data for mixing
      subroutine collect(self, cache)
         import :: mixer_type, mixer_cache_container
         !> Instance of the electronic mixer
         class(mixer_type), intent(in) :: self
         !> Cache container for mutable mixer data
         type(mixer_cache_container), intent(inout) :: cache
      end subroutine collect

      !> Apply mixing to the density/potential
      subroutine next(self, cache, error)
         import :: mixer_type, mixer_cache_container, error_type
         !> Instance of the electronic mixer
         class(mixer_type), intent(in) :: self
         !> Cache container for mutable mixer data
         type(mixer_cache_container), intent(inout) :: cache
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine next

      !> Remove cached mixer data
      subroutine cleanup(self, cache)
         import :: mixer_type, mixer_cache_container
         !> Instance of the electronic mixer
         class(mixer_type), intent(in) :: self
         !> Cache container for mutable mixer data
         type(mixer_cache_container), intent(inout) :: cache
      end subroutine cleanup
   end interface

contains


!> Set new vector from 1D array
subroutine set_1d(self, cache, x_vec)
   !> Instance of the electronic mixer
   class(mixer_type), intent(in) :: self
   !> Cache container for mutable mixer data
   type(mixer_cache_container), intent(inout) :: cache
   !> Input vector chunk for storage
   real(wp), intent(in) :: x_vec(:)

   ! Pack chunk into current input vector storage
   cache%raw%x_in(cache%raw%iset+1:cache%raw%iset+size(x_vec)) = x_vec
   ! Advance input vector index
   cache%raw%iset = cache%raw%iset + size(x_vec)
end subroutine set_1d

!> Set new vector from 2D array
subroutine set_2d(self, cache, x_vec)
   !> Instance of the electronic mixer
   class(mixer_type), intent(in) :: self
   !> Cache container for mutable mixer data
   type(mixer_cache_container), intent(inout) :: cache
   !> Input 2D array chunk for storage
   real(wp), contiguous, intent(in), target :: x_vec(:, :)
   real(wp), pointer :: x_ptr(:)
   x_ptr(1:size(x_vec)) => x_vec
   call self%set(cache, x_ptr)
end subroutine set_2d

!> Set new vector from 3D array
subroutine set_3d(self, cache, x_vec)
   !> Instance of the electronic mixer
   class(mixer_type), intent(in) :: self
   !> Cache container for mutable mixer data
   type(mixer_cache_container), intent(inout) :: cache
   !> Input 3D array chunk for storage
   real(wp), contiguous, intent(in), target :: x_vec(:, :, :)
   real(wp), pointer :: x_ptr(:)
   x_ptr(1:size(x_vec)) => x_vec
   call self%set(cache, x_ptr)
end subroutine set_3d


!> Set difference between new and old vector from 1D array
subroutine diff_1d(self, cache, x_vec)
   !> Instance of the electronic mixer
   class(mixer_type), intent(in) :: self
   !> Cache container for mutable mixer data
   type(mixer_cache_container), intent(inout) :: cache
   !> New input vector chunk for residual formation
   real(wp), intent(in) :: x_vec(:)

   if (cache%raw%initialized) then
      ! Chunk-wise build of residual dx = x_vec - x_in
      cache%raw%r_in(cache%raw%idif+1:cache%raw%idif+size(x_vec)) = x_vec &
         - cache%raw%x_in(cache%raw%idif+1:cache%raw%idif+size(x_vec))
      ! Advance residual vector index
      cache%raw%idif = cache%raw%idif + size(x_vec)
   end if
end subroutine diff_1d

!> Set difference between new and old vector from 2D array
subroutine diff_2d(self, cache, x_vec)
   !> Instance of the electronic mixer
   class(mixer_type), intent(in) :: self
   !> Cache container for mutable mixer data
   type(mixer_cache_container), intent(inout) :: cache
   !> New input 2D array chunk for residual formation
   real(wp), contiguous, intent(in), target :: x_vec(:, :)
   real(wp), pointer :: x_ptr(:)
   x_ptr(1:size(x_vec)) => x_vec
   call self%diff(cache, x_ptr)
end subroutine diff_2d

!> Set difference between new and old density from 3D array
subroutine diff_3d(self, cache, x_vec)
   !> Instance of the electronic mixer
   class(mixer_type), intent(in) :: self
   !> Cache container for mutable mixer data
   type(mixer_cache_container), intent(inout) :: cache
   !> New input 3D array chunk for residual formation
   real(wp), contiguous, intent(in), target :: x_vec(:, :, :)
   real(wp), pointer :: x_ptr(:)
   x_ptr(1:size(x_vec)) => x_vec
   call self%diff(cache, x_ptr)
end subroutine diff_3d


!> Get extrapolated vector as 1D array
subroutine get_1d(self, cache, x_vec)
   !> Instance of the electronic mixer
   class(mixer_type), intent(in) :: self
   !> Cache container for mutable mixer data
   type(mixer_cache_container), intent(inout) :: cache
   !> Extrapolated vector chunk
   real(wp), intent(out) :: x_vec(:)

   if (cache%raw%initialized) then
      ! Unpack chunk from extrapolated vector storage
      x_vec = cache%raw%x_in(cache%raw%iget+1:cache%raw%iget+size(x_vec))
      ! Advance extrapolated vector index
      cache%raw%iget = cache%raw%iget + size(x_vec)
   end if
end subroutine get_1d

!> Get extrapolated vector as 2D array
subroutine get_2d(self, cache, x_vec)
   !> Instance of the electronic mixer
   class(mixer_type), intent(in) :: self
   !> Cache container for mutable mixer data
   type(mixer_cache_container), intent(inout) :: cache
   !> Extrapolated 2D array chunk
   real(wp), contiguous, intent(out), target :: x_vec(:, :)
   real(wp), pointer :: x_ptr(:)
   x_ptr(1:size(x_vec)) => x_vec
   call self%get(cache, x_ptr)
end subroutine get_2d

!> Get extrapolated vector as 3D array
subroutine get_3d(self, cache, x_vec)
   !> Instance of the electronic mixer
   class(mixer_type), intent(in) :: self
   !> Cache container for mutable mixer data
   type(mixer_cache_container), intent(inout) :: cache
   !> Extrapolated 3D array chunk
   real(wp), contiguous, intent(out), target :: x_vec(:, :, :)
   real(wp), pointer :: x_ptr(:)
   x_ptr(1:size(x_vec)) => x_vec
   call self%get(cache, x_ptr)
end subroutine get_3d


!> Get error metric from mixing
pure function get_error(self, cache) result(error)
   !> Instance of the DIIS mixer
   class(mixer_type), intent(in) :: self
   !> Cache container for mutable mixer data
   type(mixer_cache_container), intent(in) :: cache
   !> Error metric RMS norm of residual
   real(wp) :: error

   error = 0.0_wp 
   ! Compute RMS norm of residual r_in:
   !   error = sqrt( (1/n) * sum_i r_in_i^2 )
   error = sqrt(sum(cache%raw%r_in * cache%raw%r_in) / real(cache%raw%ndim, wp))

end function get_error

end module tblite_scf_mixer_type
