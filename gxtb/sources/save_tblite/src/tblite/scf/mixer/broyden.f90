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

!> @file tblite/scf/mixer/broyden.f90
!> Provides a modified Broyden mixer implementation

!> Implementation of a modified Broyden mixing
module tblite_scf_mixer_broyden
   use mctc_env, only : wp, error_type, fatal_error
   use tblite_lapack, only : getrf, getrs
   use tblite_scf_mixer_cache, only : mixer_cache, mixer_cache_container
   use tblite_scf_mixer_input, only : mixer_input
   use tblite_scf_mixer_type, only : mixer_type
   implicit none
   private

   public :: new_broyden

   !> Default value for damping in Broyden mixing
   real(wp), parameter :: broyden_damping_default = 0.4_wp

   !> Default start of the Broyden mixing
   integer, parameter :: broyden_start_default = 2

   !> Configuration for the Broyden mixer
   type, public, extends(mixer_input) :: broyden_input
   end type broyden_input

   !> Provide constructor for the Broyden input
   interface broyden_input
      module procedure :: create_broyden_input
   end interface broyden_input


   !> Electronic mixer using modified Broyden scheme
   type, public, extends(mixer_type) :: broyden_mixer
   contains
      !> Update mixer cache
      procedure :: update
      !> Collect data for mixing
      procedure :: collect
      !> Apply mixing to the density/potential
      procedure :: next
      !> Remove cached mixer data
      procedure :: cleanup
   end type broyden_mixer

   !> Cache for mutable broyden mixer data
   type, public, extends(mixer_cache) :: broyden_cache
      !> Previous vector: [ndim]
      real(wp), allocatable :: x_last(:)
      !> Previous residual vector: [ndim]
      real(wp), allocatable :: r_last(:)

      !> Normalized change in the residual vector between iterations: [ndim, memory]
      real(wp), allocatable :: delta_r(:, :)
      !> Auxiliary vectors for modified Broyden update: [ndim, memory]
      real(wp), allocatable :: aux_vec(:, :)
      !> Overlap matrix between residual change vectors: [memory, memory]
      real(wp), allocatable :: overlap(:, :)
      !> Weight for each history element: [memory]
      real(wp), allocatable :: weights(:)
   end type broyden_cache

   !> Lower bound on the weight
   real(wp), parameter :: min_w = 1.0_wp
   !> Upper bound of the weight
   real(wp), parameter :: max_w = 1.0e5_wp
   !> Scaling factor for the weight computation
   real(wp), parameter :: w_scale = 1.0e-2_wp
   !> Diagonal regularization of the Broyden matrix
   real(wp), parameter :: reg = 1.0e-2_wp

contains


!> Constructor for the Broyden input
function create_broyden_input(mode, start, conv_start, precollect, memory, damp) &
   & result(self)
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
   !> Instance of the broyden input
   type(broyden_input) :: self

   self%mode = mode

   if (present(start)) then
      self%start = start
   else 
      self%start = broyden_start_default
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
      self%precollect = 1
   end if

   if (present(memory)) then
      self%memory = memory
   else
      self%memory = 8
   end if

   if (present(damp)) then
      self%damp = damp
   else
      self%damp = broyden_damping_default
   end if

end function create_broyden_input


!> Create new instance of the Broyden mixer
subroutine new_broyden(self, input)
   !> Instance of the Broyden mixer
   type(broyden_mixer), intent(out) :: self
   !> Configuration of the Broyden mixer
   type(broyden_input), intent(in) :: input

   self%start = input%start
   self%conv_start = input%conv_start
   self%precollect = input%precollect
   self%memory = input%memory
   self%mode = input%mode
   self%damp = input%damp
end subroutine new_broyden


!> Update mixer cache for Broyden mixer
subroutine update(self, cache, ndim)
   !> Instance of the Broyden mixer
   class(broyden_mixer), intent(in) :: self
   !> Cache container for mutable mixer data
   type(mixer_cache_container), intent(inout) :: cache
   !> Size of the extrapolated vector
   integer, intent(in) :: ndim
   
   type(broyden_cache), pointer :: ptr

   call taint(cache, ptr)

   ptr%ndim = ndim

   ! Current and last input/residual vectors
   if (.not. allocated(ptr%x_in)) then
      allocate(ptr%x_in(ndim), source=0.0_wp)
   end if
   if (.not. allocated(ptr%x_last)) then
      allocate(ptr%x_last(ndim), source=0.0_wp)
   end if
   if (.not. allocated(ptr%r_in)) then
      allocate(ptr%r_in(ndim), source=0.0_wp)
   end if
   if (.not. allocated(ptr%r_last)) then
      allocate(ptr%r_last(ndim), source=0.0_wp)
   end if

   ! History arrays for Broyden mixing
   if (.not. allocated(ptr%delta_r)) then
      allocate(ptr%delta_r(ndim, self%memory), source=0.0_wp)
   end if
   if (.not. allocated(ptr%aux_vec)) then
      allocate(ptr%aux_vec(ndim, self%memory), source=0.0_wp)
   end if

   ! Weights and overlap matrix for Broyden mixing
   if (.not. allocated(ptr%weights)) then
      allocate(ptr%weights(self%memory), source=0.0_wp)
   end if
   if (.not. allocated(ptr%overlap)) then
      allocate(ptr%overlap(self%memory, self%memory), source=0.0_wp)
   end if

end subroutine update


!> Collect data for mixing
subroutine collect(self, cache)
   !> Instance of the Broyden mixer
   class(broyden_mixer), intent(in) :: self
   !> Cache container for mutable mixer data
   type(mixer_cache_container), intent(inout) :: cache

   integer :: slot, i, sloti, current_memory
   real(wp) :: norm_r, norm_dr, inv_norm_dr

   type(broyden_cache), pointer :: ptr

   call view(cache, ptr)

   ! Update mixer internal iteration counter
   ptr%iter = ptr%iter + 1

   ! For the first iteration just initialize the previous vectors
   if (ptr%iter == 1) then
      ptr%x_last = ptr%x_in
      ptr%r_last = ptr%r_in
      return
   end if

   ! Slot in the cyclic buffer history for new entry
   slot = mod(ptr%iter - 2, self%memory) + 1

   ! Calculate weight for this slot based on the residual norm
   ! while enforcing lower and upper bounds
   norm_r = sqrt(dot_product(ptr%r_in, ptr%r_in))
   if (norm_r > (w_scale / max_w)) then
      ptr%weights(slot) = w_scale / norm_r
   else
      ptr%weights(slot) = max_w
   end if
   ptr%weights(slot) = min(max(ptr%weights(slot), min_w), max_w)


   ! Calculate change in the residual vector: r_in - r_last
   ptr%delta_r(:, slot) = ptr%r_in - ptr%r_last

   ! Calculate the norm of the residual vector change and inverse: ||r_in - r_last||
   norm_dr = max(sqrt(dot_product(ptr%delta_r(:, slot), ptr%delta_r(:, slot))), &
      & epsilon(1.0_wp))
   inv_norm_dr = 1.0_wp / norm_dr

   ! Calculate the normalized change in the residual vector:
   !   delta_r = (r_in - r_last) / ||r_in - r_last||
   ptr%delta_r(:, slot) = inv_norm_dr * ptr%delta_r(:, slot)

   ! Calculate auxiliary Broyden correction vector
   !   aux_vec = damp * delta_r + (x_in - x_last) / ||r_in - r_last||
   ptr%aux_vec(:, slot) = self%damp * ptr%delta_r(:, slot) &
      & + inv_norm_dr * (ptr%x_in - ptr%x_last)


   ! Update the overlap matrix column/row for the current slot
   current_memory = min(ptr%iter - 1, self%memory)
   do i = 1, current_memory
      ! Slot i in cyclic buffer history
      sloti = mod((ptr%iter - 1 - current_memory) + (i - 1), self%memory) + 1
      ! Calculate overlap entry:
      !   overlap(slot,sloti) = <delta_r(:, slot) | delta_r(:, sloti)>
      ptr%overlap(slot, sloti) = &
         & dot_product(ptr%delta_r(:, slot), ptr%delta_r(:, sloti))
      ptr%overlap(sloti, slot) = ptr%overlap(slot, sloti)
   end do

   ! Save the current vector and residual for the next iteration
   ptr%x_last = ptr%x_in
   ptr%r_last = ptr%r_in

end subroutine collect


!> Apply mixing to the density/potential
subroutine next(self, cache, error)
   !> Instance of the Broyden mixer
   class(broyden_mixer), intent(in) :: self
   !> Cache container for mutable mixer data
   type(mixer_cache_container), intent(inout) :: cache
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: current_memory, start_index, i, sloti, j, slotj, info
   integer, allocatable :: ipiv(:)
   real(wp), allocatable :: beta(:, :), rhs(:, :)

   type(broyden_cache), pointer :: ptr

   call view(cache, ptr)

   ! Currently available history size
   current_memory = min(ptr%iter-1, self%memory)

   ! Peform simple mixing if no history entry is available
   if (current_memory < 1) then
      ptr%x_in = ptr%x_in + self%damp * ptr%r_in
      return
   end if

   ! Memory for solution of linear system
   allocate(beta(current_memory, current_memory), rhs(current_memory, 1), &
      & source=0.0_wp)   

   ! Start index in cyclic buffer for the current history
   start_index = ptr%iter - 1 - current_memory

   ! Build Broyden matrix and right hand side
   do i = 1, current_memory
      ! Slot i in cyclic buffer history
      sloti = mod(start_index + (i - 1), self%memory) + 1

      ! Right hand side from weighted overlap of residual change and current residual
      !   rhs_i   = w_i * <delta_r_i | r_k>
      rhs(i, 1) = ptr%weights(sloti) * dot_product(ptr%delta_r(:, sloti), ptr%r_in)

      do j = 1, current_memory
         slotj = mod(start_index + (j - 1), self%memory) + 1

         ! Broyden matrix from weighted overlap of residual change vectors
         !   beta_ij = w_i*w_j*<delta_r_i|delta_r_j>
         beta(i, j) = ptr%weights(sloti) * ptr%weights(slotj) &
            & * ptr%overlap(sloti, slotj)
      end do

      ! Regularization of the diagonal entries
      !   beta_ii = beta_ii + reg^2
      beta(i, i) = beta(i, i) + reg * reg
   end do

   ! LU factorization of Broyden matrix
   allocate(ipiv(current_memory))
   call getrf(beta, ipiv, info)
   if (info /= 0) then
      call fatal_error(error, "Broyden error in LU factorization")
      return
   end if

   ! Solve linear system: beta^T * coeff = rhs
   call getrs(beta, rhs, ipiv, info, trans="T")
   if (info /= 0) then
      call fatal_error(error, "Broyden error in linear solve")
      return
   end if

   ! Start from linear mixing
   ptr%x_in = ptr%x_in + self%damp * ptr%r_in

   do i = 1, current_memory
      ! Slot i in cyclic buffer history
      sloti = mod(start_index + (i - 1), self%memory) + 1
      ! Apply Broyden correction of iteratrion i
      !   x_in = x_in - w_i * coeff_i * aux_vec_i
      ptr%x_in = ptr%x_in - ptr%weights(sloti) * rhs(i, 1) * ptr%aux_vec(:, sloti)
   end do

end subroutine next


!> Remove cached mixer data
subroutine cleanup(self, cache)
   !> Instance of the Broyden mixer
   class(broyden_mixer), intent(in) :: self
   !> Cache container for mutable mixer data
   type(mixer_cache_container), intent(inout) :: cache

   type(broyden_cache), pointer :: ptr

   call view(cache, ptr)

   if (allocated(ptr%x_in)) deallocate(ptr%x_in)
   if (allocated(ptr%r_in)) deallocate(ptr%r_in)

   if (allocated(ptr%x_last)) deallocate(ptr%x_last)
   if (allocated(ptr%r_last)) deallocate(ptr%r_last)
   if (allocated(ptr%delta_r)) deallocate(ptr%delta_r)
   if (allocated(ptr%aux_vec)) deallocate(ptr%aux_vec)
   if (allocated(ptr%weights)) deallocate(ptr%weights)
   if (allocated(ptr%overlap)) deallocate(ptr%overlap)

end subroutine cleanup


!> Inspect mixer cache and reallocate it in case of type mismatch
subroutine taint(cache, ptr)
   !> Instance of the mixer cache container
   type(mixer_cache_container), target, intent(inout) :: cache
   !> Reference to the mixer cache
   type(broyden_cache), pointer, intent(out) :: ptr

   if (allocated(cache%raw)) then
      call view(cache, ptr)
      if (associated(ptr)) return
      deallocate(cache%raw)
   end if

   if (.not.allocated(cache%raw)) then
      block
         type(broyden_cache), allocatable :: tmp
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
   type(broyden_cache), pointer, intent(out) :: ptr
   nullify(ptr)
   select type(target => cache%raw)
   type is(broyden_cache)
      ptr => target
   end select
end subroutine view

end module tblite_scf_mixer_broyden
