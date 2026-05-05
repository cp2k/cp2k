program demo
   use, intrinsic :: iso_fortran_env, only : error_unit
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : new, structure_type
   use dftd4, only : damping_param, get_rational_damping, &
      & get_dispersion, realspace_cutoff, &
      & dispersion_model, new_dispersion_model
   implicit none

   type(error_type), allocatable :: error
   type(structure_type) :: mol
   real(wp) :: energy
   real(wp), allocatable :: xyz(:, :), gradient(:, :), sigma(:, :)
   integer, allocatable :: num(:)
   integer, parameter :: nat = 5

   allocate(gradient(3, nat), sigma(3, 3))
   allocate(num(nat), xyz(3, nat))

   ! Initialize molecular structure
   num(:) = [6, 1, 1, 1, 1]
   xyz(:, :) = reshape([ &
     &  0.00000000000000_wp, -0.00000000000000_wp,  0.00000000000000_wp, &
     & -1.19220800552211_wp,  1.19220800552211_wp,  1.19220800552211_wp, &
     &  1.19220800552211_wp, -1.19220800552211_wp,  1.19220800552211_wp, &
     & -1.19220800552211_wp, -1.19220800552211_wp, -1.19220800552211_wp, &
     &  1.19220800552211_wp,  1.19220800552211_wp, -1.19220800552211_wp], &
     & [3, size(num)])

   call new(mol, num, xyz, charge=0.0_wp, uhf=0)

   ! Run calculation, check for errors, and print results
   call calc_dftd4(error, mol, "pbe", energy, gradient=gradient, &
      & sigma=sigma, model="D4")

   if (allocated(error)) then
      write(error_unit, '("[Error]:", 1x, a)') error%message
      error stop
   end if

   write (*,'(a,f18.12)') "D4 dispersion energy (Hartree): ", energy

contains

subroutine calc_dftd4(error, mol, method, energy, model, gradient, sigma)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Method name for which parameters should be used
   character(len=*), intent(in) :: method
   !> Dispersion energy
   real(wp), intent(out) :: energy
   !> Dispersion model (D4 or D4S)
   character(len=*), intent(in), optional :: model
   !> Dispersion gradient
   real(wp), intent(out), contiguous, optional :: gradient(:, :)
   !> Dispersion virial
   real(wp), intent(out), contiguous, optional :: sigma(:, :)

   class(damping_param), allocatable :: param
   class(dispersion_model), allocatable :: disp

   ! Get D4 damping parameters for given method
   call get_rational_damping(method, param, s9=1.0_wp)
   if (.not.allocated(param)) then
      call fatal_error(error, "No parameters for '"//method//"' available")
      return
   end if

   ! Initialize D4/D4S model
   call new_dispersion_model(error, disp, mol, model=model)
   if (allocated(error)) return

   call get_dispersion(mol, disp, param, realspace_cutoff(), energy, &
      & gradient, sigma)

end subroutine calc_dftd4

end program demo
