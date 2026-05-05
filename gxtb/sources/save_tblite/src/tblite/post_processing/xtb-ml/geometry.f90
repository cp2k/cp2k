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

!> @file tblite/post-processing/xtb-ml/geometry.f90
!> Geometry based xtbml features
module tblite_xtbml_geometry_based
   use, intrinsic :: iso_fortran_env, only : error_unit
   use mctc_env, only : wp, error_type
   use mctc_io, only : structure_type
   use mctc_ncoord, only : ncoord_type, new_ncoord, cn_count
   use tblite_integral_type, only : integral_type
   use tblite_output_format, only : format_string
   use tblite_container, only : container_cache
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtbml_convolution, only : xtbml_convolution_type
   use tblite_xtbml_feature_type, only : xtbml_feature_type
   implicit none
   private

   integer, parameter ::features = 1
   integer, parameter :: ext_features = 1
   character(len=23), parameter :: label = "geometry-based features"

   type, public, extends(xtbml_feature_type) :: xtbml_geometry_features_type

      real(wp), allocatable ::  cn_atom(:)
      real(wp), allocatable ::  ext_cn(:, :)

   contains
      procedure :: compute_features
      procedure :: compute_extended
      procedure :: setup
   end type xtbml_geometry_features_type

contains

!> Setup the container
subroutine setup(self)
   class(xtbml_geometry_features_type) :: self
   self%label = label
   if (allocated(self%dict)) deallocate(self%dict)
   allocate(self%dict)
   if (allocated(self%dict_ext)) deallocate(self%dict_ext)
   allocate(self%dict_ext)
   if (allocated(self%cn_atom)) deallocate(self%cn_atom)
   if (allocated(self%ext_cn)) deallocate(self%ext_cn)
end subroutine setup

!> Compute geometry-based features
subroutine compute_features(self, mol, wfn, integrals, calc, cache_list)
   !> Instance of feature container
   class(xtbml_geometry_features_type), intent(inout) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Wavefunction strcuture data
   type(wavefunction_type), intent(in) :: wfn
   !> Integral container
   type(integral_type), intent(in) :: integrals
   !> Single-point calculator
   type(xtb_calculator), intent(in) :: calc
   !> Cache list for storing caches of various interactions
   type(container_cache), intent(inout) :: cache_list(:)

   class(ncoord_type), allocatable :: ncoord
   type(error_type), allocatable :: error

   self%n_features = self%n_features + features

   allocate(self%cn_atom(mol%nat))
   call new_ncoord(ncoord, mol, cn_count%exp, error)
   if(allocated(error)) then
      write(error_unit, '("[Error]:", 1x, a)') error%message
      error stop
   end if

   call ncoord%get_cn(mol, self%cn_atom)

   call self%dict%add_entry("CN_A", self%cn_atom)

end subroutine compute_features

subroutine compute_extended(self, mol, wfn, integrals, calc, cache_list, convolution)
   !> Instance of feature container
   class(xtbml_geometry_features_type), intent(inout) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Wavefunction strcuture data
   type(wavefunction_type), intent(in) :: wfn
   !> Integral container
   type(integral_type), intent(in) :: integrals
   !> Single-point calculator
   type(xtb_calculator), intent(in) :: calc
   !> Cache list for storing caches of various interactions
   type(container_cache), intent(inout) :: cache_list(:)
   !> Convolution container
   type(xtbml_convolution_type), intent(in) :: convolution
   character(len=:), allocatable :: tmp_label
   integer :: j

   allocate(self%ext_cn(mol%nat, convolution%n_a))
   call get_ext_cn(self%cn_atom, mol%xyz, self%ext_cn, convolution)
   self%n_features = self%n_features + ext_features

   do j = 1, convolution%n_a
      tmp_label = trim("ext_CN_A"//'_'//adjustl(format_string(convolution%a(j), '(f12.2)')))
      if (tmp_label .eq. "ext_CN_A_1.00") tmp_label = "ext_CN"
      call self%dict_ext%add_entry(tmp_label, self%ext_cn(:, j))
   end do
end subroutine compute_extended

subroutine get_ext_cn(cn, xyz, ext_cn, conv)
   !> Coordinated number
   real(wp), intent(in) :: cn(:)
   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)
   !> Delta CN
   real(wp), intent(out) :: ext_cn(:, :)
   !> Convolution container
   type(xtbml_convolution_type) :: conv
   integer :: i, j, k, nat
   real(wp) :: result
   nat = size(cn)

   ext_cn = 0.0_wp
   !$omp parallel do default(none) collapse(2)&
   !$omp shared(nat, conv, cn, ext_cn)&
   !$omp private(i, j , k, result)
   do k = 1, conv%n_a
      do i = 1, nat
         do j = 1, nat
            if (i == j) cycle
            result = cn(j)/conv%kernel(i, j, k)
            !$omp atomic
            ext_cn(i, k) = ext_cn(i, k) + result
         end do
      end do
   end do
   !$omp end parallel do
end subroutine get_ext_cn

end module tblite_xtbml_geometry_based
