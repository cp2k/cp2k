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

!> @file tblite/adjlist.f90
!> Sparse neighbour map / adjacency list implementation.

!> Implementation of a sparse neighbour map in compressed sparse row format.
!>
!> A symmetric neighbour map given in dense format like
!>
!>   |   | 1 | 2 | 3 | 4 | 5 | 6 |
!>   |---|---|---|---|---|---|---|
!>   | 1 |   | x |   | x | x |   |
!>   | 2 | x |   | x |   | x | x |
!>   | 3 |   | x |   | x |   | x |
!>   | 4 | x |   | x |   | x | x |
!>   | 5 | x | x |   | x |   |   |
!>   | 6 |   | x | x | x |   |   |
!>
!> Is stored in two compressed array identifying the neighbouring atom `nlat`
!> and its cell index `nltr`. Two index arrays `inl` for the offset
!> and `nnl` for the number of entries map the atomic index to the row index.
!>
!> ```
!> inl   =  0,       3,          7,      10,         14,      17, 20
!> nnl   =  |  2 ->  |  3 ->     |  2 ->  |  3 ->     |  2 ->  |  |
!> nlat  =     2, 4, 5, 1, 3, 5, 6, 2, 4, 6, 1, 3, 5, 6, 1, 2, 4, 2, 3, 4
!> nltr  =     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
!> ```
!>
!> An alternative representation would be to store just the offsets in `inl` and
!> additional beyond the last element the total number of neighbors. However,
!> the indexing is from inl(i) to inl(i+1)-1 could be confusing, therefore
!> two arrays are used for clarity.
module tblite_adjlist
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_resize, only : resize
   implicit none
   private

   public :: adjacency_list, new_adjacency_list

   !> @class adjacency_list
   !> Neighbourlist in CSR format
   type :: adjacency_list
      !> Offset index in the neighbour map
      integer, allocatable :: inl(:)
      !> Number of neighbours for each atom
      integer, allocatable :: nnl(:)
      !> Index of the neighbouring atom
      integer, allocatable :: nlat(:)
      !> Cell index of the neighbouring atom
      integer, allocatable :: nltr(:)
   end type adjacency_list

contains

   !> Create new neighbourlist for a given geometry and cutoff
   subroutine new_adjacency_list(self, mol, trans, cutoff, complete)
      !> Instance of the neighbourlist
      type(adjacency_list), intent(out) :: self
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Translation vectors for all images
      real(wp), intent(in) :: trans(:, :)
      !> Realspace cutoff for neighbourlist generation
      real(wp), intent(in) :: cutoff
      !> Whether a complete or a symmetrical reduced map should be generated
      logical, intent(in), optional :: complete

      logical :: cmplt

      cmplt = .false.
      if (present(complete)) cmplt = complete

      allocate(self%inl(mol%nat), source=0)
      allocate(self%nnl(mol%nat), source=0)
      call generate(mol, trans, cutoff, self%inl, self%nnl, self%nlat, self%nltr, cmplt)
   end subroutine new_adjacency_list

   !> Generator for neighbourlist
   subroutine generate(mol, trans, cutoff, inl, nnl, nlat, nltr, complete)
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Translation vectors for all images
      real(wp), intent(in) :: trans(:, :)
      !> Realspace cutoff for neighbourlist generation
      real(wp), intent(in) :: cutoff
      !> Offset index in the neighbour map
      integer, intent(inout) :: inl(:)
      !> Number of neighbours for each atom
      integer, intent(inout) :: nnl(:)
      !> Index of the neighbouring atom
      integer, allocatable, intent(out) :: nlat(:)
      !> Cell index of the neighbouring atom
      integer, allocatable, intent(out) :: nltr(:)
      !> Whether a complete or a symmetrical reduced map should be generated
      logical, intent(in) :: complete

      integer :: iat, jat, itr, img
      real(wp) :: r2, vec(3), cutoff2

      img = 0
      cutoff2 = cutoff**2

      call resize(nlat, 10*mol%nat)
      call resize(nltr, 10*mol%nat)

      do iat = 1, mol%nat
         inl(iat) = img
         do jat = 1, merge(mol%nat, iat, complete)
            do itr = 1, size(trans, 2)
               vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
               r2 = sum(vec**2)
               if (r2 < epsilon(cutoff2) .or. r2 > cutoff2) cycle
               img = img + 1
               if (size(nlat) < img) call resize(nlat)
               if (size(nltr) < img) call resize(nltr)
               nlat(img) = jat
               nltr(img) = itr
            end do
         end do
         nnl(iat) = img - inl(iat)
      end do

      call resize(nlat, img)
      call resize(nltr, img)

   end subroutine generate

end module tblite_adjlist
