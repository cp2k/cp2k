program test_simple_d3
   use, intrinsic :: iso_fortran_env, only : r8 => real64
   use mctc_env, only: error_type
   use mctc_io, only: structure_type, new
   use dftd3, only: d3_model, d3_param, rational_damping_param, get_rational_damping, &
      & new_rational_damping, new_d3_model, get_dispersion, realspace_cutoff
   implicit none

   character(len=:), allocatable :: method
   type(structure_type) :: mol
   type(error_type), allocatable :: error
   integer, allocatable :: num(:)
   real(r8), allocatable :: xyz(:, :)
   real(r8) :: energy

   type(d3_model) :: disp
   type(d3_param) :: inp
   type(rational_damping_param) :: param

   method = 'PBE0'
   num = [6, 1, 1, 1, 1]
   xyz = reshape([ &  ! coordinates in Bohr
     &  0.0000000_r8, -0.0000000_r8,  0.0000000_r8, &
     & -1.1922080_r8,  1.1922080_r8,  1.1922080_r8, &
     &  1.1922080_r8, -1.1922080_r8,  1.1922080_r8, &
     & -1.1922080_r8, -1.1922080_r8, -1.1922080_r8, &
     &  1.1922080_r8,  1.1922080_r8, -1.1922080_r8],&
     & [3, size(num)])
   call new(mol, num, xyz, charge=0.0_r8, uhf=0)

   call get_rational_damping(inp, method, error, s9=1.0_r8)
   if (allocated(error)) then
      print '(2a)', "Error: ", error%message
      return
   end if
   call new_rational_damping(param, inp)
   call new_d3_model(disp, mol)

   call get_dispersion(mol, disp, param, realspace_cutoff(), energy)
   print '(3a, f13.10, a)', 'Dispersion energy for ', method, '-D3(BJ) is ', energy, ' Hartree'

end program test_simple_d3
