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

!> @file tblite/post-processing/xtb-ml/features.f90
!> Abstract base class for xtb based features
module tblite_xtbml_feature_type
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_container, only : container_cache
   use tblite_context , only : context_type
   use tblite_double_dictionary, only : double_dictionary_type
   use tblite_integral_type, only : integral_type
   use tblite_xtbml_convolution, only : xtbml_convolution_type
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_wavefunction_type, only : wavefunction_type
   implicit none
   private

   type, public, abstract :: xtbml_feature_type
      integer :: n_features
      character(len=:), allocatable :: label
      type(double_dictionary_type), allocatable :: dict, dict_ext
   contains
      procedure(compute_features), deferred :: compute_features
      procedure(compute_extended), deferred :: compute_extended
      procedure :: get_n_features
      procedure :: info
      procedure :: setup
   end type xtbml_feature_type

   character(len=*), parameter :: label = "General feature class"

   abstract interface
      subroutine compute_features(self, mol, wfn, integrals, calc, cache_list)
         import :: wp, wavefunction_type, structure_type, integral_type, xtb_calculator,&
         & container_cache, context_type, xtbml_feature_type
         !> Instance of feature container
         class(xtbml_feature_type), intent(inout) :: self
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
      end subroutine compute_features
      subroutine compute_extended(self, mol, wfn, integrals, calc, cache_list, convolution)
         import :: wp, wavefunction_type, structure_type, integral_type, xtb_calculator,&
         & container_cache, context_type, xtbml_feature_type, xtbml_convolution_type
         !> Instance of feature container
         class(xtbml_feature_type), intent(inout) :: self
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
      end subroutine compute_extended

   end interface

contains

subroutine setup(self)
   class(xtbml_feature_type) :: self
   self%label = label
end subroutine setup

function get_n_features(self) result(n)
   class(xtbml_feature_type) :: self
   integer :: n

   n = self%dict%get_n_entries()
   n = n + self%dict_ext%get_n_entries()
end function get_n_features

pure function info(self, verbosity, indent) result(str)
   !> Instance of the interaction container
   class(xtbml_feature_type), intent(in) :: self
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

end module tblite_xtbml_feature_type
