program param_scanner
   use mctc_env, only : wp, error_type, get_argument, fatal_error
   use mctc_io, only : structure_type, read_structure
   use d3_param_scan, only : reaction_type, scan_param_for_reaction
   implicit none

   type(reaction_type) :: reaction
   type(error_type), allocatable :: error
   character(:), allocatable :: method
   real(wp), allocatable :: dft_energy
   integer :: n_args, n_mol, iarg, imol

   n_args = command_argument_count()

   if (n_args < 3) then
      print '(a)', "Usage: param-scanner <method> <coeff1> <mol1> ... [dft energy]"
      stop 1
   end if

   n_mol = (n_args - 1) / 2

   call get_argument(1, method)

   allocate(reaction%mol(n_mol))
   allocate(reaction%stochiometry(n_mol))
   imol = 0
   do iarg = 1, 2 * n_mol, 2
      imol = imol + 1
      call read_real(iarg + 1, reaction%stochiometry(imol), error)
      if (allocated(error)) exit
      call read_mol(iarg + 2, reaction%mol(imol), error)
      if (allocated(error)) exit
   end do
   if (.not.allocated(error) .and. 2 * n_mol < n_args - 1) then
      allocate(dft_energy)
      call read_real(n_args, dft_energy, error)
   end if
   if (allocated(error)) then
      print '(a)', error%message
      stop 1
   end if

   call scan_param_for_reaction(error, reaction, method, dft_energy)
   if (allocated(error)) then
      print '(a)', error%message
      stop 1
   end if

contains

subroutine read_mol(idx, mol, error)
   integer, intent(in) :: idx
   type(structure_type), intent(out) :: mol
   type(error_type), allocatable, intent(out) :: error

   character(len=:), allocatable :: tmp

   call get_argument(idx, tmp)
   call read_structure(mol, tmp, error=error)
end subroutine read_mol

subroutine read_real(idx, val, error)
   integer, intent(in) :: idx
   real(wp), intent(out) :: val
   type(error_type), allocatable, intent(out) :: error

   character(len=:), allocatable :: tmp
   integer :: stat

   call get_argument(idx, tmp)
   read(tmp, *, iostat=stat) val
   if (stat /= 0) then
      call fatal_error(error, "Could not read floating point value from '"//tmp//"'")
   end if
end subroutine read_real

end program param_scanner