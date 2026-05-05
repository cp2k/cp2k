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

!> @file tblite/post-processing/xtb-ml/convolution.f90
!> Convolution of the density features
module tblite_xtbml_convolution
   use, intrinsic :: iso_fortran_env, only : error_unit
   use mctc_env, only : wp, error_type
   use mctc_data_covrad, only : get_covalent_rad
   use mctc_io, only : structure_type
   use mctc_ncoord, only : ncoord_type, new_ncoord, cn_count

   implicit none
   private
   real(wp) :: k1 = 16.0_wp
   public :: xtbml_convolution_type
   type :: xtbml_convolution_type
      real(wp), allocatable :: rcov(:)
      real(wp), allocatable :: a(:)
      real(wp), allocatable :: cn(:, :)
      integer :: n_a
      real(wp), allocatable :: kernel(:, :, :)
      character(len=:), allocatable :: label
   contains
      procedure :: setup
      procedure :: compute_kernel
      procedure, private :: populate_kernel
      procedure, private :: get_rcov
      procedure, private :: compute_cn
      procedure :: info
   end type xtbml_convolution_type
   character(len=*), parameter :: label = "CN-based convolution"


contains

subroutine setup(self)
   class(xtbml_convolution_type), intent(inout) :: self
   self%label = label
end subroutine setup

!> Compute values of the convolution kernel
subroutine compute_kernel(self, mol)
   !> Convolution container
   class(xtbml_convolution_type), intent(inout) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   call self%get_rcov(mol)
   call self%populate_kernel(mol%id, mol%xyz)
   call self%compute_cn(mol)

end subroutine compute_kernel

!> Compute the CN
subroutine compute_cn(self, mol)
   !> Convolution container 
   class(xtbml_convolution_type) :: self
   !> Molecular structure data
   type(structure_type) :: mol

   class(ncoord_type), allocatable :: ncoord
   type(error_type), allocatable :: error
   integer :: i, n_a

   n_a = size(self%a)
   allocate(self%cn(mol%nat, n_a), source=0.0_wp)
   do i = 1, n_a
      call new_ncoord(ncoord, mol, cn_count%exp, error, rcov=self%rcov*self%a(i))
      if(allocated(error)) then
         write(error_unit, '("[Error]:", 1x, a)') error%message
         error stop
      end if
   
      call ncoord%get_cn(mol, self%cn(:, i))
   end do
end subroutine compute_cn

!> Get the covalent radii
subroutine get_rcov(self, mol)
   !> Convolution container
   class(xtbml_convolution_type) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   if (allocated(self%rcov)) then
      deallocate (self%rcov)
   end if
   allocate (self%rcov(mol%nid), source=0.0_wp)
   self%rcov(:) = get_covalent_rad(mol%num)
end subroutine get_rcov

!> Populate the kernel
subroutine populate_kernel(self, at, xyz)
   !> Convolution container
   class(xtbml_convolution_type) :: self
   !> Atomic numbers
   integer, intent(in) :: at(:)
   !> Atomic coordinates
   real(wp), intent(in) :: xyz(:, :)
   real(wp) :: result
   integer :: i, j, k, n_a, nat
   n_a = size(self%a)
   nat = size(at)
   if (allocated(self%kernel)) then
      deallocate (self%kernel)
   end if
   allocate (self%kernel(nat, nat, n_a), source=0.0_wp)

   !$omp parallel do default(none) collapse(2) &
   !$omp shared(self, nat, at, xyz, n_a) &
   !$omp private(result, i, j, k)
   do k = 1, n_a
      do i = 1, nat
         do j = 1, nat
            if (i /= j) then
               call inv_cn(self, i, j, at, xyz, self%a(k), result)
               self%kernel(i, j, k) = result
            else
               self%kernel(i, j, k) = 1.0_wp
            end if
         end do
      end do
   end do

end subroutine populate_kernel

!> Compute the inverse of the coordination number
subroutine inv_cn(self, a, b, at, xyz, dampening_fact, result)
   !> Convolution container
   type(xtbml_convolution_type) :: self
   !> Atom index a
   integer, intent(in) :: a
   !> Atom index b
   integer, intent(in) :: b
   !> Atomic numbers
   integer, intent(in) :: at(:)
   !> Atomic coordinates
   real(wp), intent(in) :: xyz(:, :)
   !> Dampening factor
   real(wp), intent(in) :: dampening_fact
   !> Result
   real(wp), intent(out) :: result
   real(wp) :: rab(3), r, rco, r2

   result = 0.0_wp

   rab = xyz(:, a) - xyz(:, b)
   r2 = sum(rab**2)
   r = sqrt(r2)

   rco = dampening_fact*(self%rcov(at(a)) + self%rcov(at(b)))

   result = inv_exp_count(k1, r, rco)

end subroutine inv_cn

!> Inversed exponential counting function from D3
pure elemental function inv_exp_count(k, r, r0) result(count)
   real(wp), intent(in) :: k
   real(wp), intent(in) :: r
   real(wp), intent(in) :: r0
   real(wp) :: count
   count = 1.0_wp + exp(-k*(r0/r - 1.0_wp))
end function inv_exp_count

!> Information on the container
pure function info(self, verbosity, indent) result(str)
   !> Instance of the interaction container
   class(xtbml_convolution_type), intent(in) :: self
   !> Verbosity level
   integer, intent(in) :: verbosity
   !> Indentation level
   character(len=*), intent(in) :: indent
   !> Information on the container
   character(len=:), allocatable :: str

   if (allocated(self%label)) then
      str = indent // self%label
   else
      str = "Unknown"
   end if
end function info

end module tblite_xtbml_convolution
