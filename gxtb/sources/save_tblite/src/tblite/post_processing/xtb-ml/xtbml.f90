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

!> @file tblite/post-processing/xtb-ml/xtbml.f90
!> File to collect all xtbml features
module tblite_post_processing_xtbml_features
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_basis_type, only : basis_type
   use tblite_container, only : container_cache
   use tblite_context, only : context_type
   use tblite_double_dictionary, only : double_dictionary_type
   use tblite_integral_type, only : integral_type
   use tblite_param_xtbml_features, only : xtbml_features_record
   use tblite_post_processing_type, only : post_processing_type
   use tblite_output_format, only : format_string
   use tblite_timer, only : timer_type, format_time
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtbml_convolution, only : xtbml_convolution_type
   use tblite_xtbml_density_based, only : xtbml_density_features_type
   use tblite_xtbml_energy_features, only : xtbml_energy_features_type
   use tblite_xtbml_geometry_based, only : xtbml_geometry_features_type
   use tblite_xtbml_orbital_energy, only : xtbml_orbital_features_type
   implicit none
   private
   public :: xtbml_type, new_xtbml_features
   type, extends(post_processing_type) :: xtbml_type
      type(xtbml_geometry_features_type), allocatable :: geom
      type(xtbml_density_features_type), allocatable :: dens
      type(xtbml_orbital_features_type), allocatable :: orb
      type(xtbml_energy_features_type), allocatable :: energy
      type(xtbml_convolution_type), allocatable :: conv
   contains
      procedure :: compute
      procedure :: info
      procedure :: print_timer
   end type xtbml_type
   character(len=*), parameter :: label = "  xtbml features:"
   type(timer_type) :: timer
contains

subroutine new_xtbml_features(new_xtbml_model, param)
   type(xtbml_features_record) :: param
   type(xtbml_type), intent(inout) :: new_xtbml_model

   new_xtbml_model%label = label

   if (param%xtbml_geometry) then
      allocate(new_xtbml_model%geom)
      call new_xtbml_model%geom%setup()
   end if

   if (param%xtbml_density) then
      allocate(new_xtbml_model%dens)
      call new_xtbml_model%dens%setup()
      new_xtbml_model%dens%return_xyz = param%xtbml_tensor
   end if

   if (param%xtbml_orbital_energy) then
      allocate(new_xtbml_model%orb)
      call new_xtbml_model%orb%setup()
   end if

   if (param%xtbml_energy) then
      allocate(new_xtbml_model%energy)
      call new_xtbml_model%energy%setup()
   end if

   if (param%xtbml_convolution) then
      allocate(new_xtbml_model%conv)
      if (allocated(param%xtbml_a)) then
         new_xtbml_model%conv%a = param%xtbml_a
      else 
         new_xtbml_model%conv%a = (/1.0_wp/)
      end if
      new_xtbml_model%conv%n_a = size(new_xtbml_model%conv%a)
      call new_xtbml_model%conv%setup()
   end if

end subroutine new_xtbml_features

subroutine compute(self, mol, wfn, integrals, calc, cache_list, ctx, prlevel, dict)
   class(xtbml_type),intent(inout) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Wavefunction strcuture data
   type(wavefunction_type), intent(in) :: wfn
   !> integral container for dipole and quadrupole integrals for CAMMs
   type(integral_type), intent(in) :: integrals
   !> Single-point calculator conatiner
   type(xtb_calculator), intent(in) :: calc
   !> Context container for writing to stdout
   type(context_type), intent(inout) :: ctx
   !> Cache list for storing caches of various interactions
   type(container_cache), intent(inout) :: cache_list(:)
   !> Dictionary to store the features
   type(double_dictionary_type), intent(inout) :: dict
   !> Print level
   integer, intent(in) :: prlevel

   call timer%push("total")


   if (allocated(self%geom)) then
      call timer%push("geometry")
      associate(category => self%geom)
         call category%compute_features(mol, wfn, integrals, calc, cache_list)
         dict = dict + category%dict
         deallocate(category%dict)
      end associate
      call timer%pop()
   end if

   if (allocated(self%dens)) then
      call timer%push("density")
      associate(category => self%dens)
         call category%compute_features(mol, wfn, integrals, calc, cache_list)
         dict = dict + category%dict
         deallocate(category%dict)
      end associate
      call timer%pop()
   end if

   if (allocated(self%orb)) then
      call timer%push("orbital energy")
      associate(category => self%orb)
         call category%compute_features(mol, wfn, integrals, calc, cache_list)
         dict = dict + category%dict
         deallocate(category%dict)
      end associate
      call timer%pop()
   end if

   if (allocated(self%energy)) then
      call timer%push("energy")
      associate(category => self%energy)
         call category%compute_features(mol, wfn, integrals, calc, cache_list)
         dict = dict + category%dict
         deallocate(category%dict)
      end associate
      call timer%pop()
   end if

   if (allocated(self%conv)) then

      if (allocated(self%geom)) then
         call timer%push("geometry convolution")
         call self%conv%compute_kernel(mol)
         associate(category => self%geom)
            call category%compute_extended(mol, wfn, integrals, calc, cache_list, self%conv)
            dict = dict + category%dict_ext
            deallocate(category%dict_ext)
         end associate
         call timer%pop()
      else
         call self%conv%compute_kernel(mol)
      end if

      if (allocated(self%dens)) then
         call timer%push("density convolution")
         associate(category => self%dens)
            call category%compute_extended(mol, wfn, integrals, calc, cache_list, self%conv)
            dict = dict + category%dict_ext
            deallocate(category%dict_ext)
         end associate
         call timer%pop()
      end if

      if (allocated(self%orb)) then
         call timer%push("orbital energy convolution")
         associate(category => self%orb)
            call category%compute_extended(mol, wfn, integrals, calc, cache_list, self%conv)
            dict = dict + category%dict_ext
            deallocate(category%dict_ext)
         end associate
         call timer%pop()
      end if
   end if
   call timer%pop()
end subroutine compute

pure function info(self, verbosity, indent) result(str)
   !> Instance of the interaction container
   class(xtbml_type), intent(in) :: self
   !> Verbosity level
   integer, intent(in) :: verbosity
   !> Indentation level
   character(len=*), intent(in) :: indent
   !> Information on the container
   character(len=:), allocatable :: str
   character(len=*), parameter :: nl = new_line('a')
   character(len=:), allocatable :: category_indent
   if (allocated(self%label)) then
      str = self%label
   else
      str = "Unknown"
   end if

   category_indent = indent // " * "

   if (allocated(self%geom)) then
      str = str // nl // self%geom%info(verbosity, category_indent)
   end if

   if (allocated(self%dens)) then
      str = str // nl // self%dens%info(verbosity, category_indent)
   end if

   if (allocated(self%orb)) then
      str = str // nl // self%orb%info(verbosity, category_indent)
   end if

   if (allocated(self%energy)) then
      str = str // nl // self%energy%info(verbosity, category_indent)
   end if

   if (allocated(self%conv)) then
      str = str // nl // self%conv%info(verbosity, category_indent)
   end if
end function info

subroutine print_timer(self, prlevel, ctx)
   !> Instance of the interaction container
   class(xtbml_type), intent(in) :: self
   integer :: prlevel
   type(context_type) :: ctx
   real(wp) :: ttime, stime
   integer :: it
   character(len=*), parameter :: labels(*) = [character(len=20):: &
   & "geometry", "density", "orbital energy", "energy", "geometry convolution", &
      "density convolution", "orbital energy convolution"]

   if (prlevel > 2) then
      call ctx%message("ML features timing details:")
      ttime = timer%get("total")
      call ctx%message(" total:"//repeat(" ", 16)//format_time(ttime))
      do it = 1, size(labels)
         stime = timer%get(labels(it))
         if (stime <= epsilon(0.0_wp)) cycle
         call ctx%message(" - "//labels(it)//format_time(stime) &
         & //" ("//format_string(int(stime/ttime*100), '(i3)')//"%)")
      end do
      call ctx%message("")
   end if
end subroutine print_timer

end module tblite_post_processing_xtbml_features
