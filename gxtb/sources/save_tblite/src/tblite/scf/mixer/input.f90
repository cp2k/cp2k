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

!> @file tblite/scf/mixer/input.f90
!> Provides an electronic mixer implementation

!> Module for processing mixer input parameters
module tblite_scf_mixer_input
   use mctc_env, only : wp
   implicit none
   private

   public :: mixer_kind, mixer_mode

   !> Container for mixer inputs
   type, public :: mixer_input_container
      !> Actual mixer input
      class(mixer_input), allocatable :: raw
   end type mixer_input_container

   !> Input parameters for electronic mixer
   type, public, abstract :: mixer_input
      !> Starting iteration for mixer
      integer :: start
      !> Convergence cutoffs to start mixer
      real(wp) :: conv_start = 1e-12_wp
      !> Iterations to precollect before starting
      integer :: precollect = 0
      !> Length of history considered for extrapolation
      integer :: memory
      !> Mode of mixing either density or potential (Fock)
      integer :: mode
      !> Damping parameter
      real(wp) :: damp
   end type mixer_input

   type enum_mixers
      !> Simple mixer
      integer :: simple = 1
      !> Broyden mixer
      integer :: broyden = 2
      !> DIIS mixer
      integer :: diis = 3
   end type enum_mixers

   type enum_mixer_mode
      !> Mixing density-related properties
      integer :: density = 1
      !> Mixing potential-related properties
      integer :: potential = 2
   end type enum_mixer_mode

   type(enum_mixers), parameter :: mixer_kind = enum_mixers()
   type(enum_mixer_mode), parameter :: mixer_mode = enum_mixer_mode()

end module tblite_scf_mixer_input