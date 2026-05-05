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

!> @file tblite/scf/mixer/diis.f90
!> Provides an DIIS mixer implementation

!> Implementation of direct inversion of the iterative subspace mixing
module tblite_scf_mixer_diis
   use mctc_env, only : wp, error_type, fatal_error
   use tblite_lapack, only : getrf, getrs
   use tblite_scf_mixer_cache, only : mixer_cache, mixer_cache_container
   use tblite_scf_mixer_input, only : mixer_input
   use tblite_scf_mixer_type, only : mixer_type
   implicit none
   private

   public :: new_diis

   !> Default number of iterations mixed by DIIS
   integer, parameter :: diis_memory_default = 7

   !> Default start of the DIIS mixing
   integer, parameter :: diis_start_default = 4

   !> Default number of precollected iterations for DIIS
   integer, parameter :: diis_precollect_default = 1

   !> Default value for damping in DIIS (no damping with full DIIS step)
   real(wp), parameter :: diis_damping_default = 1.0_wp

   !> Default output fraction or relaxation parameter (full output mixing)
   real(wp), parameter :: diis_output_fraction_default = 1.0_wp

   !> Configuration for the DIIS mixer
   type, public, extends(mixer_input) :: diis_input
      !> Relaxation or output fraction parameter
      real(wp) :: output_fraction
   end type diis_input

   !> Provide constructor for the DIIS input
   interface diis_input
      module procedure :: create_diis_input
   end interface diis_input

   !> Electronic mixer using the DIIS scheme
   type, public, extends(mixer_type) :: diis_mixer
      !> Relaxation or output fraction parameter
      real(wp) :: output_fraction
   contains
      !> Update mixer cache
      procedure :: update
      !> Collect data for mixing
      procedure :: collect
      !> Apply mixing to the density/potential
      procedure :: next
      !> Remove cached mixer data
      procedure :: cleanup
   end type diis_mixer

   !> Cache for mutable DIIS mixer data
   type, public, extends(mixer_cache) :: diis_cache
      !> History buffer for input vectors (last extrapolated vector): [ndim, memory]
      real(wp), allocatable :: x_hist(:, :)
      !> History buffer for residuum vectors: [ndim, memory]
      real(wp), allocatable :: r_hist(:, :)
      !> Overlap matrix for DIIS: [memory+1, memory+1]
      real(wp), allocatable :: overlap(:, :)
   end type diis_cache

contains


!> Constructor for DIIS input
function create_diis_input(mode, start, conv_start, precollect, memory, damp, &
   & output_fraction) result(self)
   !> Mode of mixing either density or potential
   integer, intent(in) :: mode
   !> Starting iteration for mixer
   integer, intent(in), optional :: start
   !> Convergence cutoffs to start mixer
   real(wp), intent(in), optional :: conv_start
   !> Iterations to precollect before starting
   integer, intent(in), optional :: precollect
   !> Length of history considered for extrapolation
   integer, intent(in), optional :: memory
   !> Damping parameter
   real(wp), intent(in), optional :: damp
   !> Relaxation or output fraction parameter
   real(wp), intent(in), optional :: output_fraction
   !> Instance of the DIIS input
   type(diis_input) :: self

   self%mode = mode

   if (present(start)) then
      self%start = start
   else 
      self%start = diis_start_default
   end if

   if (present(conv_start)) then
      self%conv_start = conv_start
   else
      ! Activate only based on iteration count
      self%conv_start = 1.0e-12_wp
   end if

   if (present(precollect)) then
      self%precollect = precollect
   else
      self%precollect = diis_precollect_default
   end if

   if (present(memory)) then
      self%memory = memory
   else
      self%memory = diis_memory_default
   end if

   if (present(damp)) then
      self%damp = damp
   else
      self%damp = diis_damping_default
   end if

   if (present(output_fraction)) then
      self%output_fraction = output_fraction
   else
      self%output_fraction = diis_output_fraction_default
   end if

end function create_diis_input


!> Create new instance of the DIIS electronic mixer
subroutine new_diis(self, input)
   !> Instance of the mixer
   type(diis_mixer), intent(out) :: self
   !> Configuration of the DIIS mixer
   type(diis_input), intent(in) :: input

   self%start = input%start
   self%conv_start = input%conv_start
   self%precollect = input%precollect
   self%memory = input%memory
   self%mode = input%mode
   self%damp = input%damp
   self%output_fraction = input%output_fraction
end subroutine new_diis


!> Update mixer cache for DIIS mixer
subroutine update(self, cache, ndim)
   class(diis_mixer), intent(in) :: self
   type(mixer_cache_container), intent(inout) :: cache
   integer, intent(in) :: ndim

   type(diis_cache), pointer :: ptr
   logical :: wrong_dim, wrong_mem

   call taint(cache, ptr)

   ptr%ndim = ndim

   ! Current input/residual vectors
   if (.not. allocated(ptr%x_in)) then
      allocate(ptr%x_in(ndim), source=0.0_wp)
   end if
   if (.not. allocated(ptr%r_in)) then
      allocate(ptr%r_in(ndim), source=0.0_wp)
   end if

   ! History arrays for DIIS mixing
   if (.not. allocated(ptr%x_hist)) then
      allocate(ptr%x_hist(ndim, self%memory), source=0.0_wp)
   end if
   if (.not. allocated(ptr%r_hist)) then
      allocate(ptr%r_hist(ndim, self%memory), source=0.0_wp)
   end if

   ! Overlap matrix for DIIS mixing
   if (.not. allocated(ptr%overlap)) then
      allocate(ptr%overlap(self%memory, self%memory), source=0.0_wp)
   end if

end subroutine update


!> Collect data for mixing
subroutine collect(self, cache)
   !> Instance of the DIIS mixer
   class(diis_mixer), intent(in) :: self
   !> Cache container for mutable mixer data
   type(mixer_cache_container), intent(inout) :: cache

   integer :: slot, i, sloti, current_memory

   type(diis_cache), pointer :: ptr

   call view(cache, ptr)

   ! Update mixer internal iteration counter
   ptr%iter = ptr%iter + 1

   ! Slot in the cyclic buffer history for new entry
   slot = mod(ptr%iter - 1, self%memory) + 1

   ! Save the current input vector (last extrapolated vector in the history
   ptr%x_hist(:, slot) = ptr%x_in
   
   ! Save the current residual vector (new - extrapolated) in the history
   ptr%r_hist(:, slot) = ptr%r_in

   ! Update overlap column/row for this slot
   current_memory = min(ptr%iter, self%memory)
   do i = 1, current_memory
      ! Slot j in cyclic buffer history
      sloti = mod((ptr%iter - current_memory) + (i - 1), self%memory) + 1
      ! Calculate overlap entry as the Forbenius norm of the full matrix
      ! to ensure invariance w.r.t. orbital rotations and permutations:
      !   overlap(slot,sloti) = <r(:, slot) | r(:, sloti)>
      ptr%overlap(slot, sloti) = &
         & dot_product(ptr%r_hist(:, slot), ptr%r_hist(:, sloti))
      ptr%overlap(sloti, slot) = ptr%overlap(slot, sloti)
   end do

end subroutine collect


!> Apply mixing to the density/potential
subroutine next(self, cache, error)
   !> Instance of the DIIS mixer
   class(diis_mixer), intent(in) :: self
   !> Cache container for mutable mixer data
   type(mixer_cache_container), intent(inout) :: cache
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: current_memory, start_index, i, j, sloti, slotj, info
   integer, allocatable :: ipiv(:)
   real(wp), allocatable :: B(:, :), rhs(:, :)

   type(diis_cache), pointer :: ptr

   call view(cache, ptr)

   ! Currently available history size
   current_memory = min(ptr%iter, self%memory)

   ! Peform simple mixing if only one history entry is available
   if (current_memory <= 1) then
      ptr%x_in = ptr%x_in + self%damp * ptr%r_in
      return
   end if

   ! Memory for solution of linear system
   allocate(B(current_memory+1, current_memory+1), rhs(current_memory+1, 1), &
      & source=0.0_wp)

   ! Start index in cyclic buffer for the current history
   start_index = ptr%iter - current_memory

   ! Build Pulay DIIS matrix
   do i = 1, current_memory
      ! Slot i in cyclic buffer history
      sloti = mod(start_index + (i - 1), self%memory) + 1

      do j = 1, current_memory
         ! Slot j in cyclic buffer history
         slotj = mod(start_index + (j - 1), self%memory) + 1
         ! Pulay DIIS matrix entry from overlap of residuals
         B(i, j) = ptr%overlap(sloti, slotj)
      end do
      ! Add Lagrange constraint to keep the weights normalized
      B(i, current_memory+1) = -1.0_wp
      B(current_memory+1, i) = -1.0_wp
   end do
   B(current_memory+1, current_memory+1) = 0.0_wp

   ! Right-hand side for the linear system
   rhs(:, 1) = 0.0_wp
   rhs(current_memory+1, 1) = -1.0_wp

   ! LU factorization of Pulay DIIS matrix
   allocate(ipiv(current_memory+1))
   call getrf(B, ipiv, info)
   if (info /= 0) then
      call fatal_error(error, "DIIS error in LU factorization")
      return
   end if

   ! Solve linear system: B * coeff = rhs
   call getrs(B, rhs, ipiv, info)
   if (info /= 0) then
      call fatal_error(error, "DIIS error in linear solve")
      return
   end if

   ! Construct new extrapolated output vector with additional damping
   !   x_in = (1 - damp) * x_in + damp * diis_extrapolation
   ptr%x_in = (1.0_wp - self%damp) * ptr%x_in
   do i = 1, current_memory
      ! Slot i in cyclic buffer history
      sloti = mod(start_index + (i - 1), self%memory) + 1
      ! Form linear combination of either input vectors (last extrapolated vector)
      ! or raw vectors (from last actual diagonalization/Fock matrix formation)
      ! based on the specified relaxation or output fraction parameter:
      !   x_in = x_in + rhs(i, 1) * (x_hist(:, sloti) + output_fraction * r_hist(:, sloti))
      ptr%x_in = ptr%x_in + self%damp * rhs(i, 1) * (ptr%x_hist(:, sloti) &
         & + self%output_fraction * ptr%r_hist(:, sloti))
   end do

end subroutine next


!> Remove cached mixer data
subroutine cleanup(self, cache)
   !> Instance of the Broyden mixer
   class(diis_mixer), intent(in) :: self
   !> Cache container for mutable mixer data
   type(mixer_cache_container), intent(inout) :: cache

   type(diis_cache), pointer :: ptr

   call view(cache, ptr)

   if (allocated(ptr%x_in)) deallocate(ptr%x_in)
   if (allocated(ptr%r_in)) deallocate(ptr%r_in)

   if (allocated(ptr%x_hist)) deallocate(ptr%x_hist)
   if (allocated(ptr%r_hist)) deallocate(ptr%r_hist)
   if (allocated(ptr%overlap)) deallocate(ptr%overlap)

end subroutine cleanup


!> Inspect mixer cache and reallocate it in case of type mismatch
subroutine taint(cache, ptr)
   !> Instance of the mixer cache container
   type(mixer_cache_container), target, intent(inout) :: cache
   !> Reference to the mixer cache
   type(diis_cache), pointer, intent(out) :: ptr

   if (allocated(cache%raw)) then
      call view(cache, ptr)
      if (associated(ptr)) return
      deallocate(cache%raw)
   end if

   if (.not.allocated(cache%raw)) then
      block
         type(diis_cache), allocatable :: tmp
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
   type(diis_cache), pointer, intent(out) :: ptr
   nullify(ptr)
   select type(target => cache%raw)
   type is(diis_cache)
      ptr => target
   end select
end subroutine view

end module tblite_scf_mixer_diis
