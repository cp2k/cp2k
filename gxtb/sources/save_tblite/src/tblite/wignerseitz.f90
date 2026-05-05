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

!> @file tblite/wignerseitz.f90
!> Provides a Wigner-Seitz cell

!> Implementation of finding the relevant nearest neighbours in a Wigner-Seitz cell
module tblite_wignerseitz
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_cutoff, only : get_lattice_points
   implicit none
   private

   public :: new_wignerseitz_cell

   !> Information on Wigner-Seitz images
   type, public :: wignerseitz_cell
      integer, allocatable :: nimg(:, :)
      integer, allocatable :: tridx(:, :, :)
      real(wp), allocatable :: trans(:, :)
   end type wignerseitz_cell


   !> Small cutoff threshold to create only closest cells
   real(wp), parameter :: thr = sqrt(epsilon(0.0_wp))

   !> Tolerance to consider equivalent images
   real(wp), parameter :: tol = 0.01_wp


contains


subroutine new_wignerseitz_cell(self, mol)

   !> Wigner-Seitz cell instance
   type(wignerseitz_cell), intent(out) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   integer :: iat, jat, ntr, nimg
   integer, allocatable :: tridx(:)
   real(wp) :: vec(3)
   real(wp), allocatable :: trans(:, :)

   call get_lattice_points(mol%periodic, mol%lattice, thr, trans)
   ntr = size(trans, 2)
   allocate(self%nimg(mol%nat, mol%nat), self%tridx(ntr, mol%nat, mol%nat), &
      & tridx(ntr))

   !$omp parallel do default(none) schedule(runtime) collapse(2) &
   !$omp shared(mol, trans, self) private(iat, jat, vec, nimg, tridx)
   do iat = 1, mol%nat
      do jat = 1, mol%nat
         vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat)
         call get_pairs(nimg, trans, vec, tridx)
         self%nimg(jat, iat) = nimg
         self%tridx(:, jat, iat) = tridx
      end do
   end do

   call move_alloc(trans, self%trans)
   
end subroutine new_wignerseitz_cell


subroutine get_pairs(iws, trans, rij, list)
   integer, intent(out) :: iws
   real(wp), intent(in) :: rij(3)
   real(wp), intent(in) :: trans(:, :)
   integer, intent(out) :: list(:)

   logical :: mask(size(list))
   real(wp) :: dist(size(list)), vec(3), r2
   integer :: itr, img, pos

   iws = 0
   img = 0
   list(:) = 0
   mask(:) = .true.

   do itr = 1, size(trans, 2)
      vec(:) = rij - trans(:, itr)
      r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
      if (r2 < thr) cycle
      img = img + 1
      dist(img) = r2
   end do

   if (img == 0) return

   pos = minloc(dist(:img), dim=1)

   r2 = dist(pos)
   mask(pos) = .false.

   iws = 1
   list(iws) = pos
   if (img <= iws) return

   do
      pos = minloc(dist(:img), dim=1, mask=mask(:img))
      if (abs(dist(pos) - r2) > tol) exit
      mask(pos) = .false.
      iws = iws + 1
      list(iws) = pos
   end do

end subroutine get_pairs


end module tblite_wignerseitz
