module d3_param_scan
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   implicit none
   private

   public :: reaction_type
   public :: scan_param_for_reaction

   type :: reaction_type
      type(structure_type), allocatable :: mol(:)
      real(wp), allocatable :: stochiometry(:)
   end type reaction_type

contains

subroutine get_dispersion_for_reaction(reaction, param, energy)
   use dftd3, only : get_dispersion, damping_param, realspace_cutoff, &
      & d3_model, new_d3_model
   type(reaction_type), intent(in) :: reaction
   class(damping_param), intent(in) :: param
   real(wp), intent(inout) :: energy

   integer :: component
   real(wp) :: energy_component
   type(d3_model) :: disp

   do component = 1, size(reaction%mol)
      call new_d3_model(disp, reaction%mol(component))
      call get_dispersion(reaction%mol(component), disp, param, realspace_cutoff(), &
         & energy_component)
      energy = energy + energy_component * reaction%stochiometry(component)
   end do

end subroutine get_dispersion_for_reaction

subroutine scan_param_for_reaction(error, reaction, method, dft_energy)
   use mctc_env, only : error_type, fatal_error
   use dftd3, only : d3_param
   type(error_type), allocatable :: error
   type(reaction_type), intent(in) :: reaction
   character(*), intent(in) :: method
   real(wp), intent(in), optional :: dft_energy

   logical :: has_param(5)
   real(wp) :: disp_energies(2, 5)
   type(d3_param) :: inp
   integer, parameter :: zero_damping_idx = 1, rational_damping_idx = 2, &
      & mzero_damping_idx = 3, mrational_damping_idx = 4, op_damping_idx = 5
   character(*), parameter :: label(5) = [character(len=20) :: &
      & "D3(0)", "D3(BJ)", "D3M(0)", "D3M(BJ)", "D3(op)"]

   has_param(:) = .false.
   disp_energies(:, :) = 0.0_wp
   if (present(dft_energy)) disp_energies(:, :) = dft_energy

   zero_d: block
      use dftd3, only : zero_damping_param, get_zero_damping, new_zero_damping
      type(zero_damping_param) :: param
      call get_zero_damping(inp, method, error, s9=0.0_wp)
      if (allocated(error)) exit zero_d

      has_param(zero_damping_idx) = .true.
      call new_zero_damping(param, inp)
      call get_dispersion_for_reaction(reaction, param, disp_energies(1, zero_damping_idx))

      inp%s9 = 1.0_wp
      call new_zero_damping(param, inp)
      call get_dispersion_for_reaction(reaction, param, disp_energies(2, zero_damping_idx))
   end block zero_d
   if (allocated(error)) deallocate(error)

   rational_d: block
      use dftd3, only : rational_damping_param, get_rational_damping, new_rational_damping
      type(rational_damping_param) :: param
      call get_rational_damping(inp, method, error, s9=0.0_wp)
      if (allocated(error)) exit rational_d

      has_param(rational_damping_idx) = .true.
      call new_rational_damping(param, inp)
      call get_dispersion_for_reaction(reaction, param, disp_energies(1, rational_damping_idx))

      inp%s9 = 1.0_wp
      call new_rational_damping(param, inp)
      call get_dispersion_for_reaction(reaction, param, disp_energies(2, rational_damping_idx))
   end block rational_d
   if (allocated(error)) deallocate(error)

   mzero_d: block
      use dftd3, only : mzero_damping_param, get_mzero_damping, new_mzero_damping
      type(mzero_damping_param) :: param
      call get_mzero_damping(inp, method, error, s9=0.0_wp)
      if (allocated(error)) exit mzero_d

      has_param(zero_damping_idx) = .true.
      call new_mzero_damping(param, inp)
      call get_dispersion_for_reaction(reaction, param, disp_energies(1, mzero_damping_idx))

      inp%s9 = 1.0_wp
      call new_mzero_damping(param, inp)
      call get_dispersion_for_reaction(reaction, param, disp_energies(2, mzero_damping_idx))
   end block mzero_d
   if (allocated(error)) deallocate(error)

   mrational_d: block
      use dftd3, only : rational_damping_param, get_mrational_damping, new_rational_damping
      type(rational_damping_param) :: param
      call get_mrational_damping(inp, method, error, s9=0.0_wp)
      if (allocated(error)) exit mrational_d

      has_param(mrational_damping_idx) = .true.
      call new_rational_damping(param, inp)
      call get_dispersion_for_reaction(reaction, param, disp_energies(1, mrational_damping_idx))

      inp%s9 = 1.0_wp
      call new_rational_damping(param, inp)
      call get_dispersion_for_reaction(reaction, param, disp_energies(2, mrational_damping_idx))
   end block mrational_d
   if (allocated(error)) deallocate(error)

   optimizedpower_d: block
      use dftd3, only : optimizedpower_damping_param, get_optimizedpower_damping, &
         & new_optimizedpower_damping
      type(optimizedpower_damping_param) :: param
      call get_optimizedpower_damping(inp, method, error, s9=0.0_wp)
      if (allocated(error)) exit optimizedpower_d

      has_param(op_damping_idx) = .true.
      call new_optimizedpower_damping(param, inp)
      call get_dispersion_for_reaction(reaction, param, disp_energies(1, op_damping_idx))

      inp%s9 = 1.0_wp
      call new_optimizedpower_damping(param, inp)
      call get_dispersion_for_reaction(reaction, param, disp_energies(2, op_damping_idx))
   end block optimizedpower_d
   if (allocated(error)) deallocate(error)

   if (.not.any(has_param)) then
      call fatal_error(error, "No parameters found for method '"//method//"'")
      return
   end if

   block
      use mctc_io_convert, only : autokj
      integer :: ipar

      print '(a)', "Energies in kJ/mol"
      print '(66("-"))'
      print '(1x, a, t20, 2a15, a16)', "method", "E(2)", "E(2+3)", "%E(3)"
      print '(66("-"))'
      do ipar = 1, 5
         if (.not.has_param(ipar)) cycle
         print '(1x, a, t20, 3f15.3, "%")', &
            & method // "-" // trim(label(ipar)), &
            & disp_energies(:, ipar) * autokj, &
            & (disp_energies(1, ipar) - disp_energies(2, ipar)) / disp_energies(2, ipar) * 100
      end do
      print '(66("-"))'
   end block

end subroutine scan_param_for_reaction

end module d3_param_scan