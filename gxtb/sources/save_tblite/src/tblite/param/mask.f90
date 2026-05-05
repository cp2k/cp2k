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

!> @file tblite/param/mask.f90
!> Provides a general mask for the complete parameters set

!> Collection of the parameter masking
module tblite_param_mask
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io_symbols, only : to_number, symbol_length
   use tblite_param_acp, only : acp_mask
   use tblite_basis, only : basis_set
   use tblite_param_acp, only : acp_mask, count
   use tblite_param_charge, only : charge_mask, count
   use tblite_param_dispersion, only : dispersion_mask, count
   use tblite_param_element, only : element_record, element_mask, count, get_element
   use tblite_param_firstorder, only : firstorder_mask, count
   use tblite_param_fourthorder, only : fourthorder_mask, count
   use tblite_param_halogen, only : halogen_mask, count
   use tblite_param_hamiltonian, only : hamiltonian_mask, count
   use tblite_param_increment, only : increment_mask, count
   use tblite_param_multipole, only : multipole_mask, count
   use tblite_param_exchange, only : exchange_mask, count
   use tblite_param_repulsion, only : repulsion_mask, count
   use tblite_param_serde, only : serde_record
   use tblite_param_shell, only : shell_record, shell_mask, count, get_shell
   use tblite_param_spin, only : spin_mask, count
   use tblite_param_thirdorder, only : thirdorder_mask, count
   use tblite_toml, only : toml_table, toml_array, toml_key, get_value, set_value, &
      & add_table, add_array
   implicit none
   private

   public :: count, allowed_records


   type :: allowed_records
      logical :: hamiltonian = .false.
      logical :: acp = .false.
      logical :: dispersion = .false.
      logical :: increment = .false.
      logical :: repulsion = .false.
      logical :: charge = .false.
      logical :: multipole = .false.
      logical :: halogen = .false.
      logical :: firstorder = .false.
      logical :: thirdorder = .false.
      logical :: fourthorder = .false.
      logical :: exchange = .false.
      logical :: spin = .false.
   end type allowed_records


   !> Definition of the complete parameter mask
   type, public, extends(serde_record) :: param_mask
      !> Definition of the Hamiltonian, always required
      type(hamiltonian_mask), allocatable :: hamiltonian
      !> Definition of the atomic correction potential
      type(acp_mask), allocatable :: acp
      !> Definition of the dispersion correction
      type(dispersion_mask), allocatable :: dispersion
      !> Definition of the core energy increment
      type(increment_mask), allocatable :: increment
      !> Definition of the repulsion contribution
      type(repulsion_mask), allocatable :: repulsion
      !> Definition of the isotropic second-order charge interactions
      type(charge_mask), allocatable :: charge
      !> Definition of the anisotropic second-order multipolar interactions
      type(multipole_mask), allocatable :: multipole
      !> Definition of the exchange interaction
      type(exchange_mask), allocatable :: exchange
      !> Definition of the halogen bonding correction
      type(halogen_mask), allocatable :: halogen
      !> Definition of the isotropic first-order tight-binding interactions
      type(firstorder_mask), allocatable :: firstorder
      !> Definition of the isotropic third-order tight-binding interactions
      type(thirdorder_mask), allocatable :: thirdorder
      !> Definition of the isotropic fourth-order tight-binding interactions
      type(fourthorder_mask), allocatable :: fourthorder
      !> Definition of the spin-polarization interactions
      type(spin_mask), allocatable :: spin
      !> Element specific parameter masks
      type(element_mask), allocatable :: record(:)
      !> Reference to base parametrization
      type(element_record), pointer :: ref(:) => null()
   contains
      !> Read parametrization mask from TOML data structure
      procedure :: load_from_toml
      !> Write parametrization mask to TOML data structure
      procedure :: dump_to_toml
   end type param_mask


   interface count
      module procedure :: count_mask
   end interface count


contains


!> Read parametrization mask from TOML data structure
subroutine load_from_toml(self, table, error)
   !> Instance of the parametrization data
   class(param_mask), intent(inout) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(allowed_records) :: allowed
   type(toml_table), pointer :: child

   call get_value(table, "hamiltonian", child, requested=.false.)
   allowed%hamiltonian = associated(child)
   if (associated(child)) then
      allocate(self%hamiltonian)
   end if

   call get_value(table, "acp", child, requested=.false.)
   allowed%acp = associated(child)
   if (associated(child)) then
      allocate(self%acp)
   end if

   call get_value(table, "dispersion", child, requested=.false.)
   allowed%dispersion = associated(child)
   if (associated(child)) then
      allocate(self%dispersion)
   end if

   call get_value(table, "increment", child, requested=.false.)
   allowed%increment = associated(child)
   if (associated(child)) then
      allocate(self%increment)
   end if   

   call get_value(table, "repulsion", child, requested=.false.)
   allowed%repulsion = associated(child)
   if (associated(child)) then
      allocate(self%repulsion)
   end if

   call get_value(table, "halogen", child, requested=.false.)
   allowed%halogen = associated(child)
   if (associated(child)) then
      allocate(self%halogen)
   end if

   call get_value(table, "charge", child, requested=.false.)
   allowed%charge = associated(child)
   if (associated(child)) then
      allocate(self%charge)
   end if

   call get_value(table, "multipole", child, requested=.false.)
   allowed%multipole = associated(child)
   if (associated(child)) then
      allocate(self%multipole)
   end if

   call get_value(table, "exchange", child, requested=.false.)
   allowed%exchange = associated(child)
   if (associated(child)) then
      allocate(self%exchange)
   end if

   call get_value(table, "spin", child, requested=.false.)
   allowed%spin = associated(child)
   if (associated(child)) then
      allocate(self%spin)
   end if

   call get_value(table, "firstorder", child, requested=.false.)
   allowed%firstorder = associated(child)
   if (associated(child)) then
      allocate(self%firstorder)
   end if

   call get_value(table, "thirdorder", child, requested=.false.)
   allowed%thirdorder = associated(child)
   if (associated(child)) then
      allocate(self%thirdorder)
   end if

   call get_value(table, "fourthorder", child, requested=.false.)
   allowed%fourthorder = associated(child)
   if (associated(child)) then
      allocate(self%fourthorder)
   end if

   call get_value(table, "element", child)
   call masks_from_table(self%record, child, self%ref, allowed, error)
   if (allocated(error)) return
end subroutine load_from_toml


!> Deserialize masks from a table by iterating over all entires
subroutine masks_from_table(mask, table, ref, allowed, error)
   !> List of all element masks
   type(element_mask), allocatable, intent(out) :: mask(:)
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> List of all element records
   type(element_record), intent(in), optional :: ref(:)
   !> Allowed entries in record
   type(allowed_records), intent(in) :: allowed
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: ii, ir
   type(toml_key), allocatable :: list(:)
   type(toml_table), pointer :: child

   call table%get_keys(list)
   allocate(mask(size(list)))
   do ii = 1, size(list)
      call get_value(table, list(ii)%key, child)
      if (present(ref)) then
         call get_element(ref, list(ii)%key, ir)
         if (ir == 0) cycle
      end if
      call mask_from_table(mask(ii), child, ref(ir), allowed, error)
      if (allocated(error)) exit
   end do
   if (allocated(error)) return
end subroutine masks_from_table


!> Read one mask from a table
subroutine mask_from_table(mask, table, ref, allowed, error)
   !> List of all element masks
   type(element_mask), intent(out) :: mask
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Corresponding element record
   type(element_record), intent(in), optional :: ref
   !> Allowed entries in record
   type(allowed_records), intent(in) :: allowed
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table), pointer :: child

   mask%sym = table%key
   mask%num = to_number(table%key)

   ! Increment
   call read_atom_mask(table, "increment", mask%increment, allowed%increment, error)

   ! Basis set
   call read_atom_mask(table, "k0", mask%k0, allowed%hamiltonian, error)
   if (allocated(error)) return
   call read_atom_mask(table, "k1", mask%k1, allowed%hamiltonian, error)
   if (allocated(error)) return
   call read_atom_mask(table, "k2", mask%k2, allowed%hamiltonian, error)
   if (allocated(error)) return
   call read_atom_mask(table, "k3", mask%k3, allowed%hamiltonian, error)
   if (allocated(error)) return
   call read_atom_mask(table, "k0_h0_scale", mask%k0_h0_scale, allowed%hamiltonian, error)
   if (allocated(error)) return
   call read_atom_mask(table, "k1_h0_scale", mask%k1_h0_scale, allowed%hamiltonian, error)
   if (allocated(error)) return
   call read_atom_mask(table, "k2_h0_scale", mask%k2_h0_scale, allowed%hamiltonian, error)
   if (allocated(error)) return
   call read_atom_mask(table, "k3_h0_scale", mask%k3_h0_scale, allowed%hamiltonian, error)
   if (allocated(error)) return

   ! Coordination number
   call read_atom_mask(table, "cn_rcov", mask%cn_rcov, allowed%hamiltonian, error)
   if (allocated(error)) return
   call read_atom_mask(table, "cn_avg", mask%cn_avg, allowed%hamiltonian, error)
   if (allocated(error)) return

   ! H0
   call read_atom_mask(table, "h0_rad", mask%h0_rad, allowed%hamiltonian, error)
   if (allocated(error)) return
   call read_atom_mask(table, "h0_dip_scale", mask%h0_dip_scale, allowed%hamiltonian, error)
   if (allocated(error)) return
   call read_atom_mask(table, "h0_diat_scale_sig", mask%h0_diat_scale_sig, allowed%hamiltonian, error)
   if (allocated(error)) return
   call read_atom_mask(table, "h0_diat_scale_pi", mask%h0_diat_scale_pi, allowed%hamiltonian, error)
   if (allocated(error)) return
   call read_atom_mask(table, "h0_diat_scale_del", mask%h0_diat_scale_del, allowed%hamiltonian, error)
   if (allocated(error)) return
   call read_atom_mask(table, "rvdw_scale", mask%rvdw_scale, allowed%hamiltonian, error)
   if (allocated(error)) return

   ! ACPs
   allocate(mask%acp_expos(mask%n_acp), source=.false.)
   allocate(mask%acp_levels(mask%n_acp), source=.false.)
   if(mask%n_acp > 0) then
      call get_value(table, "acp", child)

      call read_acp_mask(child, "acp_expos", mask%acp_expos, allowed%acp, error)
      if (allocated(error)) return
      call read_acp_mask(child, "acp_levels", mask%acp_levels, allowed%acp, error)
      if (allocated(error)) return
   end if

   ! Repulsion
   call read_atom_mask(table, "zeff", mask%zeff, allowed%repulsion, error)
   if (allocated(error)) return
   call read_atom_mask(table, "arep", mask%alpha, allowed%repulsion, error)
   if (allocated(error)) return
   call read_atom_mask(table, "rep_cn", mask%rep_cn, allowed%repulsion, error)
   if (allocated(error)) return
   call read_atom_mask(table, "rep_q", mask%rep_q, allowed%repulsion, error)
   if (allocated(error)) return
   call read_atom_mask(table, "rep_roffset", mask%rep_roffset, allowed%repulsion, error)
   if (allocated(error)) return
   call read_atom_mask(table, "rep_k1", mask%rep_k1, allowed%repulsion, error)
   if (allocated(error)) return

   ! Halogen bonding correction
   call read_atom_mask(table, "xbond", mask%xbond, allowed%halogen, error)
   if (allocated(error)) return

   ! First-order tight-binding
   call read_atom_mask(table, "ipea_cn", mask%ipea_cn, allowed%firstorder, error)
   if (allocated(error)) return

   ! Second-order tight-binding
   call read_atom_mask(table, "gam", mask%gam, allowed%charge, error)
   if (allocated(error)) return
   call read_atom_mask(table, "gam_cn", mask%gam_cn, allowed%charge, error)
   if (allocated(error)) return

   ! Third-order tight-binding
   call read_atom_mask(table, "gam3", mask%gam3, allowed%thirdorder, error)
   if (allocated(error)) return

   ! Fourth-order tight-binding
   call read_atom_mask(table, "gam4", mask%gam4, allowed%fourthorder, error)
   if (allocated(error)) return

   ! Multipoles
   call read_atom_mask(table, "dkernel", mask%dkernel, allowed%multipole, error)
   if (allocated(error)) return
   call read_atom_mask(table, "qkernel", mask%qkernel, allowed%multipole, error)
   if (allocated(error)) return
   call read_atom_mask(table, "aes_dip_scale", mask%aes_dip_scale, allowed%multipole, error)
   if (allocated(error)) return

   ! Spin polarization
   call read_atom_mask(table, "wll_scale", mask%wll_scale, allowed%spin, error)

   ! Fock exchange
   call read_atom_mask(table, "cscale", mask%cscale, allowed%exchange, error)
   if (allocated(error)) return
   call read_atom_mask(table, "crad", mask%crad, allowed%exchange, error)
   if (allocated(error)) return

   ! Shells
   call get_value(table, "shells", child)
   call shell_masks_from_table(mask%record, child, ref%record, allowed, error)
   if (allocated(error)) return
end subroutine mask_from_table


subroutine read_atom_mask(table, key, mask, default, error)
   type(toml_table), intent(inout) :: table
   character(len=*), intent(in) :: key
   logical, intent(out) :: mask
   logical, intent(in) :: default
   type(error_type), allocatable, intent(out) :: error

   call get_value(table, key, mask, default)
end subroutine read_atom_mask


subroutine read_acp_mask(table, key, mask, default, error)
   type(toml_table), intent(inout) :: table
   character(len=*), intent(in) :: key
   logical, intent(out) :: mask(:)
   logical, intent(in) :: default
   type(error_type), allocatable, intent(out) :: error

   type(toml_array), pointer :: array
   integer :: ii

   call get_value(table, key, array, requested=.false.)
   if (associated(array)) then
      do ii = 1, size(mask)
         call get_value(array, ii, mask(ii))
      end do
   else
      call get_value(table, key, mask(1), default)
      mask(:) = mask(1)
   end if
end subroutine read_acp_mask


!> Deserialize shell masks from a table by iterating over all entires
subroutine shell_masks_from_table(mask, table, ref, allowed, error)
   !> List of all shell masks
   type(shell_mask), allocatable, intent(out) :: mask(:)
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> List of all shell records
   type(shell_record), intent(in), optional :: ref(:)
   !> Allowed entries in record
   type(allowed_records), intent(in) :: allowed
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: ii, ir, ngauss, stat, pqn, lsh
   type(toml_key), allocatable :: list(:)
   type(toml_table), pointer :: child

   call table%get_keys(list)
   allocate(mask(size(list)))
   do ii = 1, size(list)
      call get_value(table, list(ii)%key, child)
      if (present(ref)) then
         pqn = get_pqn(list(ii)%key(1:1))
         lsh = get_lsh(list(ii)%key(2:2))
         call get_shell(ref, pqn, lsh, ir)
         if (ir == 0) cycle
         ngauss = ref(ir)%ngauss
      else
         call get_value(child, "ngauss", ngauss, stat=stat)
         if (stat /= 0) cycle
      end if
      call shell_mask_from_table(mask(ii), child, ngauss, allowed, error)
      if (allocated(error)) exit
   end do
   if (allocated(error)) return
end subroutine shell_masks_from_table


!> Read one shell mask from a table
subroutine shell_mask_from_table(mask, table, ngauss, allowed, error)
   !> List of all element masks
   type(shell_mask), intent(out) :: mask
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Number of Gaussian functions for this shell
   integer, intent(in) :: ngauss
   !> Allowed entries in record
   type(allowed_records), intent(in) :: allowed
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=:), allocatable :: ctmp
   type(toml_table), pointer :: child

   ! Read the shell angular momentum and principal quantum number
   call table%get_key(ctmp)

   mask%pqn = get_pqn(ctmp(1:1))
   if (mask%pqn < 1) then
      call fatal_error(error, "Invalid principal quantum number in "//trim(ctmp))
      return
   end if

   mask%lsh = get_lsh(ctmp(2:2))
   if (mask%lsh < 0) then
      call fatal_error(error, "Invalid angular momentum in "//trim(ctmp))
      return
   end if

   ! Check for basis set input
   mask%basis = basis_set%sto_ng
   call get_value(table, "sto_ng", child, requested=.false.)
   if (.not. associated(child)) then
      mask%basis = basis_set%q_vszp
      call get_value(table, "q_vszp", child, requested=.false.)
   end if
   if (.not. associated(child)) then
      mask%basis = basis_set%custom
      call get_value(table, "basis", child, requested=.false.)
   end if
   if (.not.associated(child)) then
      call fatal_error(error, "No basis set specified in mask.")
      return
   end if

   if (mask%basis == basis_set%sto_ng) then
      ! Stewart STO-NG basis set
      call read_shell_mask(child, "slater", mask%slater, allowed%hamiltonian, error)
      if (allocated(error)) return
   
   else if (mask%basis == basis_set%q_vszp) then
      ! Charge-dependent q-vSZP basis set
      allocate(mask%expos(ngauss), source=.false.)
      allocate(mask%coeffs(ngauss), source=.false.)
      allocate(mask%coeffs_env(ngauss), source=.false.)

      call read_prim_mask(child, "expos", mask%expos, allowed%hamiltonian, error)
      if (allocated(error)) return

      call read_prim_mask(child, "coeffs", mask%coeffs, allowed%hamiltonian, error)
      if (allocated(error)) return

      call read_prim_mask(child, "coeffs_env", mask%coeffs_env, allowed%hamiltonian, error)
      if (allocated(error)) return

   else 
      ! Custom basis set
      allocate(mask%expos(ngauss), source=.false.)
      allocate(mask%coeffs(ngauss), source=.false.)

      call read_prim_mask(child, "expos", mask%expos, allowed%hamiltonian, error)
      if (allocated(error)) return

      call read_prim_mask(child, "coeffs", mask%coeffs, allowed%hamiltonian, error)
      if (allocated(error)) return

   end if 

   ! H0
   call read_shell_mask(table, "level", mask%level, allowed%hamiltonian, error)
   if (allocated(error)) return
   call read_shell_mask(table, "kcn", mask%kcn, allowed%hamiltonian, error)
   if (allocated(error)) return
   call read_shell_mask(table, "shpoly", mask%shpoly, allowed%hamiltonian, error)
   if (allocated(error)) return
   call read_shell_mask(table, "shpoly2", mask%shpoly2, allowed%hamiltonian, error)
   if (allocated(error)) return
   call read_shell_mask(table, "shpoly4", mask%shpoly4, allowed%hamiltonian, error)
   if (allocated(error)) return
   call read_shell_mask(table, "h0_exp_scale", mask%h0_exp_scale, allowed%hamiltonian, error)
   if (allocated(error)) return

   ! First-order tight-binding
   call read_shell_mask(table, "ipea", mask%ipea, allowed%charge, error)
   if (allocated(error)) return

   ! Second-order tight-binding
   call read_shell_mask(table, "lgam", mask%lgam, allowed%charge, error)
   if (allocated(error)) return

   ! Approximated Fock exchange
   call read_shell_mask(table, "lgam_fx", mask%lgam_fx, allowed%exchange, error)
   if (allocated(error)) return
   call read_shell_mask(table, "avg_exp_fx", mask%avg_exp_fx, allowed%exchange, error)
   if (allocated(error)) return

end subroutine shell_mask_from_table


pure function get_pqn(str) result(pqn)
   character, intent(in) :: str
   integer :: pqn
   select case(str)
   case default
      pqn = 0
   case("1")
      pqn = 1
   case("2")
      pqn = 2
   case("3")
      pqn = 3
   case("4")
      pqn = 4
   case("5")
      pqn = 5
   case("6")
      pqn = 6
   case("7")
      pqn = 7
   end select
end function get_pqn

pure function get_lsh(str) result(lsh)
   character, intent(in) :: str
   integer :: lsh
   select case(str)
   case default
      lsh = -1
   case("s", "S")
      lsh = 0
   case("p", "P")
      lsh = 1
   case("d", "D")
      lsh = 2
   case("f", "F")
      lsh = 3
   case("g", "G")
      lsh = 4
   end select
end function get_lsh


subroutine read_shell_mask(table, key, mask, default, error)
   type(toml_table), intent(inout) :: table
   character(len=*), intent(in) :: key
   logical, intent(out) :: mask
   logical, intent(in) :: default
   type(error_type), allocatable, intent(out) :: error

   integer :: stat

   call get_value(table, key, mask, default)
end subroutine read_shell_mask

subroutine read_prim_mask(table, key, mask, default, error)
   type(toml_table), intent(inout) :: table
   character(len=*), intent(in) :: key
   logical, intent(out) :: mask(:)
   logical, intent(in) :: default
   type(error_type), allocatable, intent(out) :: error

   type(toml_array), pointer :: array
   integer :: ii, stat

   call get_value(table, key, array, requested=.false., stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Error while reading "//trim(key)//" mask")
      return
   end if
   if (associated(array)) then
      do ii = 1, size(mask)
         call get_value(array, ii, mask(ii))
      end do
   else
      call get_value(table, key, mask(1), default)
      mask(:) = mask(1)
   end if
end subroutine read_prim_mask


!> Write parametrization mask to TOML datastructure
subroutine dump_to_toml(self, table, error)
   !> Instance of the parametrization data
   class(param_mask), intent(in) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table), pointer :: child

   if (allocated(self%hamiltonian)) then
      call add_table(table, "hamiltonian", child)
   end if

   if (allocated(self%acp)) then
      call add_table(table, "acp", child)
   end if

   if (allocated(self%dispersion)) then
      call add_table(table, "dispersion", child)
   end if

   if (allocated(self%increment)) then
      call add_table(table, "increment", child)
   end if
   
   if (allocated(self%repulsion)) then
      call add_table(table, "repulsion", child)
   end if

   if (allocated(self%halogen)) then
      call add_table(table, "halogen", child)
   end if

   if (allocated(self%charge)) then
      call add_table(table, "charge", child)
   end if

   if (allocated(self%multipole)) then
      call add_table(table, "multipole", child)
   end if

   if (allocated(self%exchange)) then
      call add_table(table, "exchange", child)
   end if

   if (allocated(self%spin)) then
      call add_table(table, "spin", child)
   end if

   if (allocated(self%firstorder)) then
      call add_table(table, "firstorder", child)
   end if

   if (allocated(self%thirdorder)) then
      call add_table(table, "thirdorder", child)
   end if

   if (allocated(self%fourthorder)) then
      call add_table(table, "fourthorder", child)
   end if

   call add_table(table, "element", child)
   call masks_to_table(self%record, child, error)
   if (allocated(error)) return

end subroutine dump_to_toml


!> Serialize masks to a table by iterating over all entries
subroutine masks_to_table(mask, table, error)
   !> List of all element masks
   type(element_mask), intent(in) :: mask(:)
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: ii
   type(toml_table), pointer :: child

   do ii = 1, size(mask)
      call add_table(table, trim(mask(ii)%sym), child)
      call mask_to_table(mask(ii), child, error)
      if (allocated(error)) exit
   end do
   if (allocated(error)) return
end subroutine masks_to_table

!> Read one mask into a table
subroutine mask_to_table(mask, table, error)
   !> List of all element mask
   type(element_mask), intent(in) :: mask
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: i, n_acp
   type(toml_array), pointer :: array
   type(toml_table), pointer :: child

   ! Increment
   call set_value(table, "increment", mask%increment)

   ! Basis set
   call set_value(table, "k0", mask%k0)
   call set_value(table, "k1", mask%k1)
   call set_value(table, "k2", mask%k2)
   call set_value(table, "k3", mask%k3)
   call set_value(table, "k0_h0_scale", mask%k0_h0_scale)
   call set_value(table, "k1_h0_scale", mask%k1_h0_scale)
   call set_value(table, "k2_h0_scale", mask%k2_h0_scale)
   call set_value(table, "k3_h0_scale", mask%k3_h0_scale)
   
   ! Coordination number
   call set_value(table, "cn_rcov", mask%cn_rcov)
   call set_value(table, "cn_avg", mask%cn_avg)

   ! H0
   call set_value(table, "h0_rad", mask%h0_rad)
   call set_value(table, "h0_dip_scale", mask%h0_dip_scale)
   call set_value(table, "h0_diat_scale_sig", mask%h0_diat_scale_sig)
   call set_value(table, "h0_diat_scale_pi", mask%h0_diat_scale_pi)
   call set_value(table, "h0_diat_scale_del", mask%h0_diat_scale_del)
   call set_value(table, "rvdw_scale", mask%rvdw_scale)

   ! ACPs
   n_acp = size(mask%acp_expos)
   if (n_acp > 0) then
      call add_table(table, "acp", child)

      call add_array(table, "acp_expos", array)
      do i = 1, n_acp
         call set_value(array, i, mask%acp_expos(i))
      end do

      call add_array(table, "acp_levels", array)
      do i = 1, n_acp
         call set_value(array, i, mask%acp_levels(i))
      end do
   end if

   ! Repulsion
   call set_value(table, "zeff", mask%zeff)
   call set_value(table, "arep", mask%alpha)
   call set_value(table, "rep_cn", mask%rep_cn)
   call set_value(table, "rep_q", mask%rep_q)
   call set_value(table, "rep_roffset", mask%rep_roffset)
   call set_value(table, "rep_k1", mask%rep_k1)

   ! Halogen bonding correction 
   call set_value(table, "xbond", mask%xbond)

   ! First-order tight-binding
   call set_value(table, "ipea_cn", mask%ipea_cn)

   ! Second-order tight-binding
   call set_value(table, "gam", mask%gam)
   call set_value(table, "gam_cn", mask%gam_cn)

   ! Third-order tight-binding
   call set_value(table, "gam3", mask%gam3)

   ! Fourth-order tight-binding
   call set_value(table, "gam4", mask%gam4)

   ! Multipoles
   call set_value(table, "dkernel", mask%dkernel)
   call set_value(table, "qkernel", mask%qkernel)
   call set_value(table, "aes_dip_scale", mask%aes_dip_scale)

   ! Spin polarization
   call set_value(table, "wll_scale", mask%wll_scale)

   ! Fock exchange
   call set_value(table, "cscale", mask%cscale)
   call set_value(table, "crad", mask%crad)

   ! Shells
   call add_table(table, "shells", child)
   call shell_masks_to_table(mask%record, child, error)
   if (allocated(error)) return

end subroutine mask_to_table

!> Serialize shell masks to a table by iterating over all entries
subroutine shell_masks_to_table(mask, table, error)
   !> List of all shell masks
   type(shell_mask), intent(in) :: mask(:)
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: ii
   type(toml_table), pointer :: child
   character, parameter :: lsh(0:4) = ["s", "p", "d", "f", "g"]
   character, parameter :: pqn(1:7) = ["1", "2", "3", "4", "5", "6", "7"]

   do ii = 1, size(mask)
      call add_table(table, trim(pqn(mask(ii)%pqn)//lsh(mask(ii)%lsh)), child)
      call shell_mask_to_table(mask(ii), child, error)
      if (allocated(error)) exit
   end do
   if (allocated(error)) return
end subroutine shell_masks_to_table


!> Read one shell mask into a table
subroutine shell_mask_to_table(mask, table, error)
   !> List of all shell mask
   type(shell_mask), intent(in) :: mask
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: i, ngauss
   type(toml_array), pointer :: array
   type(toml_table), pointer :: child

   if (mask%basis == basis_set%sto_ng) then
      ! Stewart STO-NG basis set
      call add_table(table, "sto_ng", child)

      call set_value(child, "slater", mask%slater)

   else if (mask%basis == basis_set%q_vszp) then
      ! Charge-dependent q-vSZP basis set
      call add_table(table, "q_vszp", child)

      ngauss = size(mask%expos)   

      call add_array(child, "expos", array)
      do i = 1, ngauss
         call set_value(array, i, mask%expos(i))
      end do   

      call add_array(child, "coeffs", array)
      do i = 1, ngauss
         call set_value(array, i, mask%coeffs(i))
      end do   

      call add_array(child, "coeffs_env", array)
      do i = 1, ngauss
         call set_value(array, i, mask%coeffs_env(i))
      end do

   else
      ! Custom basis set
      call add_table(table, "basis", child)      

      ngauss = size(mask%expos)

      call add_array(child, "expos", array)
      do i = 1, ngauss
         call set_value(array, i, mask%expos(i))
      end do   

      call add_array(child, "coeffs", array)
      do i = 1, ngauss
         call set_value(array, i, mask%coeffs(i))
      end do   

   end if

   ! H0
   call set_value(table, "level", mask%level)
   call set_value(table, "kcn", mask%kcn)   
   call set_value(table, "shpoly", mask%shpoly)
   call set_value(table, "shpoly2", mask%shpoly2)
   call set_value(table, "shpoly4", mask%shpoly4)
   call set_value(table, "h0_exp_scale", mask%h0_exp_scale)

   ! First-order tight-binding
   call set_value(table, "ipea", mask%ipea)

   ! Second-order tight-binding
   call set_value(table, "lgam", mask%lgam)

   ! Approximated Fock exchange
   call set_value(table, "lgam_fx", mask%lgam_fx)
   call set_value(table, "avg_exp_fx", mask%avg_exp_fx)

end subroutine shell_mask_to_table


elemental function count_mask(mask) result(ncount)
   type(param_mask), intent(in) :: mask
   integer :: ncount
   ncount = 0
   if (allocated(mask%hamiltonian)) ncount = ncount + count(mask%hamiltonian)
   if (allocated(mask%acp)) ncount = ncount + count(mask%acp)
   if (allocated(mask%dispersion)) ncount = ncount + count(mask%dispersion)
   if (allocated(mask%increment)) ncount = ncount + count(mask%increment)
   if (allocated(mask%repulsion)) ncount = ncount + count(mask%repulsion)
   if (allocated(mask%halogen)) ncount = ncount + count(mask%halogen)
   if (allocated(mask%charge)) ncount = ncount + count(mask%charge)
   if (allocated(mask%multipole)) ncount = ncount + count(mask%multipole)
   if (allocated(mask%exchange)) ncount = ncount + count(mask%exchange)
   if (allocated(mask%spin)) ncount = ncount + count(mask%spin)
   if (allocated(mask%firstorder)) ncount = ncount + count(mask%firstorder)
   if (allocated(mask%thirdorder)) ncount = ncount + count(mask%thirdorder)
   if (allocated(mask%fourthorder)) ncount = ncount + count(mask%fourthorder)
   if (allocated(mask%record)) ncount = ncount + sum(count(mask%record))
end function count_mask


end module tblite_param_mask
