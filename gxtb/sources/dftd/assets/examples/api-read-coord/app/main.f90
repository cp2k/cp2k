! This file is part of dftd4.
! SPDX-Identifier: LGPL-3.0-or-later
!
! Simple D4 energy evaluation using the Fortran API.

program d4_energy
   use mctc_env, only : wp, error_type
   use mctc_io, only : structure_type, read_structure
   use dftd4, only : d4_model, new_d4_model, damping_param, &
      & get_rational_damping, realspace_cutoff, get_dispersion
   implicit none

   type(structure_type) :: mol
   type(d4_model) :: model
   class(damping_param), allocatable :: param
   type(realspace_cutoff) :: cutoff
   type(error_type), allocatable :: error
   real(wp) :: energy

   call read_structure(mol, "coord", error)
   if (allocated(error)) then
      write (*, '(a)') trim(error%message)
      stop 1
   end if

   call new_d4_model(error, model, mol)
   if (allocated(error)) then
      write (*, '(a)') trim(error%message)
      stop 1
   end if

   call get_rational_damping("pbe", param)
   call get_dispersion(mol, model, param, cutoff, energy)

   write (*, '(a,f18.12)') "D4 dispersion energy (Hartree): ", energy
end program d4_energy
