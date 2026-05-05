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

!> @file tblite/scf/mixer/simple.f90
!> Provides an simple mixer implementation

!> Implementation of simple mixing
module tblite_scf_mixer_simple
   use mctc_env, only : wp, error_type, fatal_error
   use tblite_scf_mixer_cache, only : mixer_cache, mixer_cache_container
   use tblite_scf_mixer_input, only : mixer_input
   use tblite_scf_mixer_type, only : mixer_type
   implicit none
   private

   public :: new_simple

   !> Default value for self-consistent iteration mixing
   real(wp), parameter :: mixer_damping_default = 0.4_wp

   !> Default start of the simple mixing
   integer, parameter :: simple_start_default = 1

   !> Configuration for the simple mixer
   type, public, extends(mixer_input) :: simple_input
   end type simple_input

   !> Provide constructor for the simple input
   interface simple_input
      module procedure :: create_simple_input
   end interface simple_input

   !> Electronic mixer using the simple mixing scheme
   type, public, extends(mixer_type) :: simple_mixer
   contains
      !> Update mixer cache
      procedure :: update
      !> Collect data for mixing
      procedure :: collect
      !> Apply mixing to the density/potential
      procedure :: next
      !> Remove cached mixer data
      procedure :: cleanup
   end type simple_mixer

   !> Cache for mutable simple mixer data
   type, public, extends(mixer_cache) :: simple_cache
   end type simple_cache

contains


!> Constructor for simpe input
function create_simple_input(mode, start, conv_start, damp) result(self)
   !> Mode of mixing either density or potential
   integer, intent(in) :: mode
   !> Starting iteration for mixer
   integer, intent(in), optional :: start
   !> Convergence cutoffs to start mixer
   real(wp), intent(in), optional :: conv_start
   !> Damping parameter
   real(wp), intent(in), optional :: damp
   !> Instance of the simple input
   type(simple_input) :: self

   self%mode = mode

   if (present(start)) then
      self%start = start
   else 
      self%start = simple_start_default
   end if

   if (present(conv_start)) then
      self%conv_start = conv_start
   else
      ! Activate only based on iteration count
      self%conv_start = 1.0e-12_wp
   end if

   self%precollect = 0
   self%memory = 1

   if (present(damp)) then
      self%damp = damp
   else
      self%damp = mixer_damping_default
   end if

end function create_simple_input


!> Create new instance of the simple electronic mixer
subroutine new_simple(self, input)
   !> Instance of the mixer
   type(simple_mixer), intent(out) :: self
   !> Configuration of the simple mixer
   type(simple_input), intent(in) :: input

   self%start = input%start
   self%conv_start = input%conv_start
   self%precollect = input%precollect
   self%memory = input%memory
   self%mode = input%mode
   self%damp = input%damp
end subroutine new_simple


!> Update mixer cache for simple mixer
subroutine update(self, cache, ndim)
   !> Instance of the simple mixer
   class(simple_mixer), intent(in) :: self
   !> Cache container for mutable mixer data
   type(mixer_cache_container), intent(inout) :: cache
   !> Size of the extrapolated vector
   integer, intent(in) :: ndim

   type(simple_cache), pointer :: ptr

   call taint(cache, ptr)

   ptr%ndim = ndim

   ! Current input/residual vectors
   if (.not. allocated(ptr%x_in)) then
      allocate(ptr%x_in(ndim), source=0.0_wp)
   end if
   if (.not. allocated(ptr%r_in)) then
      allocate(ptr%r_in(ndim), source=0.0_wp)
   end if

end subroutine update


!> Collect data for mixing
subroutine collect(self, cache)
   !> Instance of the simple mixer
   class(simple_mixer), intent(in) :: self
   !> Cache container for mutable mixer data
   type(mixer_cache_container), intent(inout) :: cache

   type(simple_cache), pointer :: ptr

   call view(cache, ptr)

   ! Update mixer internal iteration counter
   ptr%iter = ptr%iter + 1

end subroutine collect


!> Apply mixing to the density/potential
subroutine next(self, cache, error)
   !> Instance of the simple mixer
   class(simple_mixer), intent(in) :: self
   !> Cache container for mutable mixer data
   type(mixer_cache_container), intent(inout) :: cache
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(simple_cache), pointer :: ptr

   call view(cache, ptr)

   ! Linear simple mixing step:
   !   x_in = x_in + damp * r_in
   ptr%x_in = ptr%x_in + self%damp * ptr%r_in

end subroutine next


!> Remove cached mixer data
subroutine cleanup(self, cache)
   !> Instance of the simple mixer
   class(simple_mixer), intent(in) :: self
   !> Cache container for mutable mixer data
   type(mixer_cache_container), intent(inout) :: cache

   type(simple_cache), pointer :: ptr

   call view(cache, ptr)

   if (allocated(ptr%x_in)) deallocate(ptr%x_in)
   if (allocated(ptr%r_in)) deallocate(ptr%r_in)

end subroutine cleanup


!> Inspect mixer cache and reallocate it in case of type mismatch
subroutine taint(cache, ptr)
   !> Instance of the mixer cache container
   type(mixer_cache_container), target, intent(inout) :: cache
   !> Reference to the mixer cache
   type(simple_cache), pointer, intent(out) :: ptr

   if (allocated(cache%raw)) then
      call view(cache, ptr)
      if (associated(ptr)) return
      deallocate(cache%raw)
   end if

   if (.not.allocated(cache%raw)) then
      block
         type(simple_cache), allocatable :: tmp
         allocate(tmp)
         call move_alloc(tmp, cache%raw)
      end block
   end if

   call view(cache, ptr)
end subroutine taint

!> Return reference to mixer cache after resolving its type
subroutine view(cache, ptr)
   !> Instance of the mixer cache container
   type(mixer_cache_container), target, intent(in) :: cache
   !> Reference to the mixer cache
   type(simple_cache), pointer, intent(out) :: ptr
   nullify(ptr)
   select type(target => cache%raw)
   type is(simple_cache)
      ptr => target
   end select
end subroutine view

end module tblite_scf_mixer_simple
