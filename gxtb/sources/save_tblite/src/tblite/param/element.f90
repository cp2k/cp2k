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

!> @file tblite/param/element.f90
!> Provides records for the element specific parameters

!> Definition of the element specific parameter records
module tblite_param_element
   use mctc_data_atomicrad, only : get_atomic_rad
   use mctc_data_covrad, only : get_covalent_rad
   use mctc_data_paulingen, only : get_pauling_en
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io_symbols, only : to_number, symbol_length
   use tblite_param_serde, only : serde_record
   use tblite_param_shell, only : shell_record, shell_mask, count, get_shell
   use tblite_toml, only : toml_table, toml_key, toml_array, get_value, &
      & set_value, add_array, add_table, len
   implicit none
   private

   public :: count, get_element


   real(wp), parameter :: increment_default = 0.0_wp, xbond_default = 0.0_wp, &
      & mprad_default = 5.0_wp, wll_scale_default= 1.0_wp, ipea_cn_default = 0.0_wp, &
      & gam_cn_default = 0.0_wp, gam4_default = 0.0_wp, cscale_default = 0.0_wp, & 
      & cn_avg_default = 0.0_wp

   character(len=*), parameter :: k_increment = "increment", &
      & k_k0 = "k0", k_k1 = "k1", k_k2 = "k2", k_k3 = "k3", &
      & k_k0_h0_scale = "k0_h0_scale", k_k1_h0_scale = "k1_h0_scale", &
      & k_k2_h0_scale = "k2_h0_scale", k_k3_h0_scale = "k3_h0_scale", &
      & k_cn_rcov = "cn_rcov", k_cn_avg = "cn_avg", &
      & k_h0_dip_scale = "h0_dip_scale", &
      & k_h0_rad = "h0_rad", k_rvdw_scale = "rvdw_scale", &
      & k_h0_diat_scale_sig = "h0_diat_scale_sig", k_h0_diat_scale_pi = "h0_diat_scale_pi", &
      & k_h0_diat_scale_del = "h0_diat_scale_del", &
      & k_acp = "acp", k_l_acp = "l_acp", k_acp_expos = "acp_expos", k_acp_levels = "acp_levels", &  
      & k_zeff = "zeff", k_arep = "arep", k_rep_cn = "rep_cn", k_rep_q = "rep_q", &
      & k_rep_roffset = "rep_roffset", k_rep_k1="rep_k1", &
      & k_xbond = "xbond", &
      & k_ipea_cn = "ipea_cn", &
      & k_gam = "gam", k_gam_cn = "gam_cn", k_gam3 = "gam3", k_gam4 = "gam4", &
      & k_dkernel = "dkernel", k_qkernel = "qkernel", k_mprad = "mprad", k_mpvcn = "mpvcn", &
      & k_aes_dip_scale = "aes_dip_scale", & 
      & k_wll_scale = "wll_scale", &
      & k_cscale = "cscale", k_crad = "crad", &
      & k_en = "en", k_shells = "shells"

   !> Representation of the element specific parameters
   type, public, extends(serde_record) :: element_record
      !> Element symbol of specie represented by this record
      character(len=symbol_length) :: sym = ''
      !> Atomic number of the specie represented by this record
      integer :: num = 0

      !> Atomic increment tot hte tight-binding energy
      real(wp) :: increment = 0.0_wp

      !> Overall charge dependence of q-vSZP
      real(wp) :: k0 = 0.0_wp
      !> Quadratic charge dependence of q-vSZP
      real(wp) :: k1 = 0.0_wp
      !> Square root CN dependence of q-vSZP
      real(wp) :: k2 = 0.0_wp
      !> Mixed charge and CN dependence of q-vSZP
      real(wp) :: k3 = 0.0_wp
      !> H0 scaling of overall charge dependence of q-vSZP
      real(wp) :: k0_h0_scale = 0.0_wp
      !> H0 scaling of quadratic charge dependence of q-vSZP
      real(wp) :: k1_h0_scale = 0.0_wp
      !> H0 scaling of square root CN dependence of q-vSZP
      real(wp) :: k2_h0_scale = 0.0_wp
      !> H0 scaling of mixed charge and CN dependence of q-vSZP
      real(wp) :: k3_h0_scale = 0.0_wp

      !> Covalent radius for the coordination number
      real(wp) :: cn_rcov = 0.0_wp
      !> Average coordination number in reference set
      real(wp) :: cn_avg = 0.0_wp

      !> Atomic radius used in H0
      real(wp) :: h0_rad = 0.0_wp
      !> H0 dipole correction
      real(wp) :: h0_dip_scale = 0.0_wp
      !> H0 diatomic frame scaling for sigma interactions
      real(wp) :: h0_diat_scale_sig = 0.0_wp
      !> H0 diatomic frame scaling for pi interactions
      real(wp) :: h0_diat_scale_pi = 0.0_wp
      !> H0 diatomic frame scaling for delta interactions
      real(wp) :: h0_diat_scale_del = 0.0_wp
      !> Van der Waals radius scaling
      real(wp) :: rvdw_scale = 1.0_wp

      !> Number of atomic correction potential projectors
      integer :: n_acp = 0
      !> Angular momentum for atomic correction potential projector
      integer, allocatable :: l_acp(:)
      !> Atomic correction potential projector exponents
      real(wp), allocatable :: acp_expos(:)
      !> Atomic correction potential projector levels
      real(wp), allocatable :: acp_levels(:)

      !> Effective nuclear charge used in repulsion
      real(wp) :: zeff = 0.0_wp
      !> Repulsion damping exponent
      real(wp) :: alpha = 0.0_wp
      !> Coefficient for linear CN dependence of repulsion
      real(wp) :: rep_cn = 0.0_wp
      !> Coefficient for linear atomic charge dependence of repulsion
      real(wp) :: rep_q = 0.0_wp
      !> Offset atomic radii for repulsion
      real(wp) :: rep_roffset = 0.0_wp
      !> First-order expansion coefficient for the repulsion
      real(wp) :: rep_k1 = 0.0_wp

      !> Halogen bonding strength
      real(wp) :: xbond = 0.0_wp

      !> CN-dependence of IP/EA for first-order tight-binding
      real(wp) :: ipea_cn = 0.0_wp

      !> Chemical hardness / Hubbard parameter
      real(wp) :: gam = 0.0_wp
      !> CN-dependence of Chemical hardness / Hubbard parameter
      real(wp) :: gam_cn = 0.0_wp

      !> Atomic Hubbard derivative
      real(wp) :: gam3 = 0.0_wp

      !> Atomic Hubbard second derivative
      real(wp) :: gam4 = 0.0_wp
      
      !> Dipolar exchange-correlation kernel
      real(wp) :: dkernel = 0.0_wp
      !> Quadrupolar exchange-correlation kernel
      real(wp) :: qkernel = 0.0_wp
      !> Multipole radius
      real(wp) :: mprad = 0.0_wp
      !> Multipole valence CN
      real(wp) :: mpvcn = 0.0_wp
      !> Multipole dipole scaling
      real(wp) :: aes_dip_scale = 0.0_wp

      !> Scaling factor for spin constants
      real(wp) :: wll_scale = 1.0_wp

      !> Bond-order correlation scaling factor
      real(wp) :: cscale = 0.0_wp
      !> Bond-order correlation radius 
      real(wp) :: crad = 0.0_wp

      !> Electronnegativity
      real(wp) :: en = 0.0_wp

      !> Number of valence and polarization shells
      integer :: nsh = 0
      !> Shell specific parameter records
      type(shell_record), allocatable :: record(:)
   contains
      generic :: load => load_from_array
      generic :: dump => dump_to_array
      !> Read parametrization data from TOML data structure
      procedure :: load_from_toml
      !> Write parametrization data to TOML data structure
      procedure :: dump_to_toml
      !> Read parametrization data from parameter array
      procedure, private :: load_from_array
      !> Write parametrization data to parameter array
      procedure, private :: dump_to_array
   end type


   !> Masking for the element record
   type, public :: element_mask
      !> Element symbol of species represented by this record
      character(len=symbol_length) :: sym = ''
      !> Atomic number of the species represented by this record
      integer :: num = 0

      !> Atomic increment tot hte tight-binding energy
      logical :: increment

      !> Overall charge dependence of q-vSZP
      logical :: k0
      !> Quadratic charge dependence of q-vSZP
      logical :: k1
      !> Square root CN dependence of q-vSZP
      logical :: k2
      !> Mixed charge and CN dependence of q-vSZP
      logical :: k3
      !> H0 scaling of overall charge dependence of q-vSZP
      logical :: k0_h0_scale
      !> H0 scaling of quadratic charge dependence of q-vSZP
      logical :: k1_h0_scale
      !> H0 scaling of square root CN dependence of q-vSZP
      logical :: k2_h0_scale
      !> H0 scaling of mixed charge and CN dependence of q-vSZP
      logical :: k3_h0_scale

      !> Covalent radius for the coordination number
      logical :: cn_rcov
      !> Average coordination number in reference set
      logical :: cn_avg

      !> Atomic radius used in H0
      logical :: h0_rad
      !> H0 dipole correction
      logical :: h0_dip_scale
      !> H0 diatomic frame scaling for sigma interactions
      logical :: h0_diat_scale_sig
      !> H0 diatomic frame scaling for pi interactions
      logical :: h0_diat_scale_pi
      !> H0 diatomic frame scaling for delta interactions
      logical :: h0_diat_scale_del
      !> Van der Waals radius scaling
      logical :: rvdw_scale

      !> Number of atomic correction potential projectors
      integer :: n_acp = 0
      !> Atomic correction potential projector exponents
      logical, allocatable :: acp_expos(:)
      !> Atomic correction potential projector levels
      logical, allocatable :: acp_levels(:)

      !> Effective nuclear charge used in repulsion
      logical :: zeff
      !> Repulsion damping exponent
      logical :: alpha
      !> Coefficient for linear CN dependence of repulsion
      logical :: rep_cn
      !> Coefficient for linear atomic charge dependence of repulsion
      logical :: rep_q
      !> Offset atomic radii for repulsion
      logical :: rep_roffset
      !> First-order expansion coefficient for the repulsion
      logical :: rep_k1

      !> Halogen bonding strength
      logical :: xbond

      !> CN-dependence of IP/EA for first-order tight-binding
      logical :: ipea_cn

      !> Chemical hardness / Hubbard parameter
      logical :: gam
      !> CN-dependence of Chemical hardness / Hubbard parameter
      logical :: gam_cn

      !> Atomic Hubbard derivative
      logical :: gam3

      !> Atomic Hubbard second derivative
      logical :: gam4

      !> Dipolar exchange-correlation kernel
      logical :: dkernel
      !> Quadrupolar exchange-correlation kernel
      logical :: qkernel
      !> Multipole dipole scaling
      logical :: aes_dip_scale

      !> Scaling factor for spin constants
      logical :: wll_scale

      !> Bond-order correlation scaling factor
      logical :: cscale
      !> Bond-order correlation radius
      logical :: crad

      !> Shell specific parameter masks
      type(shell_mask), allocatable :: record(:)
   end type element_mask


   interface count
      module procedure :: count_mask
   end interface count


contains


!> Read parametrization data from TOML data structure
subroutine load_from_toml(self, table, error)
   class(element_record), intent(inout) :: self
   type(toml_table), intent(inout) :: table
   type(error_type), allocatable, intent(out) :: error

   character(len=:), allocatable :: sym
   type(toml_table), pointer :: child
   integer :: stat

   call table%get_key(sym)
   self%sym = sym
   self%num = to_number(sym)

   ! Increment
   call get_value(table, k_increment, self%increment, increment_default, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid energy increment for "//trim(self%sym))
      return
   end if

   ! Basis set
   call get_basis(self, table, error)
   if (allocated(error)) return

   ! Coordination number
   call get_value(table, k_cn_rcov, self%cn_rcov, get_covalent_rad(self%sym), stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid covalent radius for CN for "//trim(self%sym))
      return
   end if
   call get_value(table, k_cn_avg, self%cn_avg, cn_avg_default, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid average CN on the reference set for "//trim(self%sym))
      return
   end if

   ! H0
   call get_hamiltonian(self, table, error)
   if (allocated(error)) return

   ! ACPs
   call get_value(table, k_acp, child, requested=.false.)
   if (associated(child)) then
      call get_acp(self, child, error)
      if (allocated(error)) return
   else
      ! No ACPs are defined
      self%n_acp = 0
   end if

   ! Repulsion
   call get_repulsion(self, table, error)
   if (allocated(error)) return

   ! Halogen bonding correction
   call get_value(table, k_xbond, self%xbond, xbond_default, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid halogen bonding strength for "//trim(self%sym))
      return
   end if

   ! First-order tight-binding
   call get_value(table, k_ipea_cn, self%ipea_cn, ipea_cn_default, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid IP/EA CN-dependence for "//trim(self%sym))
      return
   end if

   ! Second-order tight-binding
   call get_value(table, k_gam, self%gam, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid atomic Hubbard parameter for "//trim(self%sym))
      return
   end if

   call get_value(table, k_gam_cn, self%gam_cn, gam_cn_default, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid CN-dependence of atomic Hubbard parameter for "//trim(self%sym))
      return
   end if

   ! Third-order tight-binding
   call get_value(table, k_gam3, self%gam3, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid atomic Hubbard derivative for "//trim(self%sym))
      return
   end if

   ! Fourth-order tight-binding
   call get_value(table, k_gam4, self%gam4, gam4_default, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid atomic Hubbard second derivative for "//trim(self%sym))
      return
   end if

   ! Multipoles
   call get_multipole(self, table, error)
   if (allocated(error)) return

   ! Spin polarization
   call get_value(table, k_wll_scale, self%wll_scale, wll_scale_default, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid spin constant scaling for "//trim(self%sym))
      return
   end if

   ! Static correlation
   call get_value(table, k_cscale, self%cscale, cscale_default, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid bond-order correlation scaling for "//trim(self%sym))
      return
   end if

   call get_value(table, k_crad, self%crad, get_covalent_rad(self%sym), stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid bond-order correlation radius for "//trim(self%sym))
      return
   end if

   ! Electronnegativity
   call get_value(table, k_en, self%en, get_pauling_en(self%sym), stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid electronegativity for "//trim(self%sym))
      return
   end if

   ! Shell specific parameters
   call get_value(table, k_shells, child)
   call records_from_table(self%record, child, error)
   if (allocated(error)) then 
      ! Augment the error message with the element symbol
      error%message = error%message//" for element "//trim(self%sym)
      return
   end if 
   self%nsh = size(self%record)
   if (self%nsh < 1) then
      call fatal_error(error, "No shells defined for "//trim(self%sym))
      return
   end if

end subroutine load_from_toml

subroutine get_basis(self, table, error)
   type(element_record), intent(inout) :: self
   type(toml_table), intent(inout) :: table
   type(error_type), allocatable, intent(out) :: error

   integer :: stat

   ! Atomic charge dependence parameters for q-vSZP
   call get_value(table, k_k0, self%k0, 0.0_wp, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid overall charge dependence for "//trim(self%sym))
      return
   end if
   call get_value(table, k_k1, self%k1, 0.0_wp, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid quadratic charge dependence for "//trim(self%sym))
      return
   end if
   call get_value(table, k_k2, self%k2, 0.0_wp, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid square root CN dependence for "//trim(self%sym))
      return
   end if
   call get_value(table, k_k3, self%k3, 0.0_wp, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid mixed charge and CN dependence for "//trim(self%sym))
      return
   end if

   ! H0 scaling of the charge dependence
   call get_value(table, k_k0_h0_scale, self%k0_h0_scale, 0.0_wp, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid H0 scaling of overall charge dependence for "//trim(self%sym))
      return
   end if
   call get_value(table, k_k1_h0_scale, self%k1_h0_scale, 0.0_wp, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid H0 scaling of quadratic charge dependence for "//trim(self%sym))
      return
   end if  
   call get_value(table, k_k2_h0_scale, self%k2_h0_scale, 0.0_wp, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid H0 scaling of square root CN dependence for "//trim(self%sym))
      return
   end if
   call get_value(table, k_k3_h0_scale, self%k3_h0_scale, 0.0_wp, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid H0 scaling of mixed charge and CN dependence for "//trim(self%sym))
      return
   end if

end subroutine get_basis


subroutine get_hamiltonian(self, table, error)
   type(element_record), intent(inout) :: self
   type(toml_table), intent(inout) :: table
   type(error_type), allocatable, intent(out) :: error

   integer :: stat

   ! H0 covalent radius
   call get_value(table, k_h0_rad, self%h0_rad, get_atomic_rad(self%sym), stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid H0 atomic radius for "//trim(self%sym))
      return
   end if

   ! H0 anisotropic correction
   call get_value(table, k_h0_dip_scale, self%h0_dip_scale, 0.0_wp, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid H0 dipole correction for "//trim(self%sym))
      return
   end if

   ! H0 diatomic frame scaling
   call get_value(table, k_h0_diat_scale_sig, self%h0_diat_scale_sig, 0.0_wp, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid H0 diatomic frame scaling for sigma interactions for "//trim(self%sym))
      return
   end if
   call get_value(table, k_h0_diat_scale_pi, self%h0_diat_scale_pi, 0.0_wp, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid H0 diatomic frame scaling for pi interactions for "//trim(self%sym))
      return
   end if
   call get_value(table, k_h0_diat_scale_del, self%h0_diat_scale_del, 0.0_wp, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid H0 diatomic frame scaling for delta interactions for "//trim(self%sym))
      return
   end if

   ! Van der Waals radius scaling
   call get_value(table, k_rvdw_scale, self%rvdw_scale, 1.0_wp, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid van der Waals radius scaling for "//trim(self%sym))
      return
   end if

end subroutine get_hamiltonian


subroutine get_acp(self, table, error)
   type(element_record), intent(inout) :: self
   type(toml_table), intent(inout) :: table
   type(error_type), allocatable, intent(out) :: error

   integer :: iacp
   type(toml_array), pointer :: array
   character(len=:), allocatable :: tmp

   ! Check the number of ACPs and get angular momentum
   call get_value(table, k_l_acp, array)
   if (.not.associated(array)) then
      call fatal_error(error, "No ACPs defined for "//trim(self%sym))
      return
   end if
   self%n_acp = len(array)
   if (self%n_acp == 0) then
      call fatal_error(error, "No ACPs defined for "//trim(self%sym))
      return
   end if

   allocate(self%l_acp(self%n_acp))
   do iacp = 1, self%n_acp
      call get_value(array, iacp, tmp)
      self%l_acp(iacp) = get_lsh(tmp)
   end do
   if (any(self%l_acp < 0)) then
      call fatal_error(error, "Invalid ACP angular momentum for "//trim(self%sym))
      return
   end if

   ! Get the ACP exponents
   call get_value(table, k_acp_expos, array)
   if (.not.associated(array)) then
      call fatal_error(error, "No ACP exponents for "//trim(self%sym))
      return
   end if
   if (len(array) < self%n_acp) then
      call fatal_error(error, "Insufficient ACP exponents for "//trim(self%sym))
      return
   end if
   allocate(self%acp_expos(self%n_acp))
   do iacp = 1, self%n_acp
      call get_value(array, iacp, self%acp_expos(iacp))
   end do

   ! Get the ACP coefficients
   call get_value(table, k_acp_levels, array)
   if (.not.associated(array)) then
      call fatal_error(error, "No ACP coefficients for "//trim(self%sym))
      return
   end if
   if (len(array) < self%n_acp) then
      call fatal_error(error, "Insufficient ACP coefficients for "//trim(self%sym))
      return
   end if
   allocate(self%acp_levels(self%n_acp))
   do iacp = 1, self%n_acp
      call get_value(array, iacp, self%acp_levels(iacp))
   end do

end subroutine get_acp

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

subroutine get_repulsion(self, table, error)
   type(element_record), intent(inout) :: self
   type(toml_table), intent(inout) :: table
   type(error_type), allocatable, intent(out) :: error

   integer :: stat

   call get_value(table, k_zeff, self%zeff, real(self%num, wp), stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid effective nuclear charge for "//trim(self%sym))
      return
   end if
   call get_value(table, k_arep, self%alpha, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid repulsion exponent for "//trim(self%sym))
      return
   end if
   call get_value(table, k_rep_cn, self%rep_cn, 0.0_wp, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid repulsion linear CN-dependence for "//trim(self%sym))
      return
   end if
   call get_value(table, k_rep_q, self%rep_q, 0.0_wp, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid repulsion linear charge-dependence for "//trim(self%sym))
      return
   end if
   call get_value(table, k_rep_roffset, self%rep_roffset, 0.0_wp, stat=stat)
   if (stat /= 0) then
     call fatal_error(error, "Invalid offset radius for repulsion for "//trim(self%sym))
      return
   end if
   call get_value(table, k_rep_k1, self%rep_k1, 1.0_wp, stat=stat)
   if (stat /= 0) then
     call fatal_error(error, "Invalid first-order expansion coefficient for repulsion for "//trim(self%sym))
      return
   end if
end subroutine get_repulsion

subroutine get_multipole(self, table, error)
   type(element_record), intent(inout) :: self
   type(toml_table), intent(inout) :: table
   type(error_type), allocatable, intent(out) :: error

   integer :: stat

   call get_value(table, k_dkernel, self%dkernel, 0.0_wp, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid dipole kernel for "//trim(self%sym))
      return
   end if
   call get_value(table, k_qkernel, self%qkernel, 0.0_wp, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid quadrupole kernel for "//trim(self%sym))
      return
   end if
   call get_value(table, k_mprad, self%mprad, mprad_default, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid multipole damping radius for "//trim(self%sym))
      return
   end if
   call get_value(table, k_mpvcn, self%mpvcn, 0.0_wp, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid multipole valence CN for "//trim(self%sym))
      return
   end if
   call get_value(table, k_aes_dip_scale, self%aes_dip_scale, 1.0_wp, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Invalid multipole dipole scaling for "//trim(self%sym))
      return
   end if
end subroutine get_multipole


!> Write parametrization data to TOML data structure
subroutine dump_to_toml(self, table, error)
   class(element_record), intent(in) :: self
   type(toml_table), intent(inout) :: table
   type(error_type), allocatable, intent(out) :: error

   integer :: i
   type(toml_array), pointer :: array
   type(toml_table), pointer :: child
   character, parameter :: lsh(0:4) = ["s", "p", "d", "f", "g"]

   ! Increment
   call set_value(table, k_increment, self%increment)

   ! Basis set
   call set_value(table, k_k0, self%k0)
   call set_value(table, k_k1, self%k1)
   call set_value(table, k_k2, self%k2)
   call set_value(table, k_k3, self%k3)
   call set_value(table, k_k0_h0_scale, self%k0_h0_scale)
   call set_value(table, k_k1_h0_scale, self%k1_h0_scale)
   call set_value(table, k_k2_h0_scale, self%k2_h0_scale)
   call set_value(table, k_k3_h0_scale, self%k3_h0_scale)

   ! Coordination number
   call set_value(table, k_cn_rcov, self%cn_rcov)
   call set_value(table, k_cn_avg, self%cn_avg)

   ! H0
   call set_value(table, k_h0_rad, self%h0_rad)
   call set_value(table, k_h0_dip_scale, self%h0_dip_scale)
   call set_value(table, k_h0_diat_scale_sig, self%h0_diat_scale_sig)
   call set_value(table, k_h0_diat_scale_pi, self%h0_diat_scale_pi)
   call set_value(table, k_h0_diat_scale_del, self%h0_diat_scale_del)
   call set_value(table, k_rvdw_scale, self%rvdw_scale)

   ! ACPs
   if (self%n_acp > 0) then
      call add_table(table, k_acp, child)

      call add_array(child, k_l_acp, array)
      do i = 1, self%n_acp
         call set_value(array, i, lsh(self%l_acp(i)))
      end do

      call add_array(child, k_acp_expos, array)
      do i = 1, self%n_acp
         call set_value(array, i, self%acp_expos(i))
      end do

      call add_array(child, k_acp_levels, array)
      do i = 1, self%n_acp
         call set_value(array, i, self%acp_levels(i))
      end do
   end if

   ! Repulsion
   call set_value(table, k_zeff, self%zeff)
   call set_value(table, k_arep, self%alpha)
   call set_value(table, k_rep_cn, self%rep_cn)
   call set_value(table, k_rep_q, self%rep_q)
   call set_value(table, k_rep_roffset, self%rep_roffset)
   call set_value(table, k_rep_k1, self%rep_k1)

   ! Halogen bonding correction
   call set_value(table, k_xbond, self%xbond)

   ! First-order tight-binding
   call set_value(table, k_ipea_cn, self%ipea_cn)

   ! Second-order tight-binding
   call set_value(table, k_gam, self%gam)
   call set_value(table, k_gam_cn, self%gam_cn)

   ! Third-order tight-binding
   call set_value(table, k_gam3, self%gam3)

   ! Fourth-order tight-binding
   call set_value(table, k_gam4, self%gam4)

   ! Multipoles
   call set_value(table, k_dkernel, self%dkernel)
   call set_value(table, k_qkernel, self%qkernel)
   call set_value(table, k_mprad, self%mprad)
   call set_value(table, k_mpvcn, self%mpvcn)
   call set_value(table, k_aes_dip_scale, self%aes_dip_scale)

   ! Spin polarization
   call set_value(table, k_wll_scale, self%wll_scale)

   ! Bond-order correlation
   call set_value(table, k_cscale, self%cscale)
   call set_value(table, k_crad, self%crad)

   ! Electronnegativity
   call set_value(table, k_en, self%en)

   ! Shells
   call add_table(table, k_shells, child)
   call records_to_table(self%record, child, error)
   if (allocated(error)) return
end subroutine dump_to_toml

!> Deserialize records from a table by iterating over all entires
subroutine records_from_table(record, table, error)
   !> List of all shell records
   type(shell_record), allocatable, intent(out) :: record(:)
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: ii
   type(toml_key), allocatable :: list(:)
   type(toml_table), pointer :: child

   call table%get_keys(list)
   allocate(record(size(list)))

   do ii = 1, size(list)
      call get_value(table, list(ii)%key, child)
      call record(ii)%load(child, error)
      if (allocated(error)) exit
   end do
   if (allocated(error)) return
end subroutine records_from_table


!> Serialize records to a table by iterating over all entries
subroutine records_to_table(record, table, error)
   !> List of all element records
   type(shell_record), intent(in) :: record(:)
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: ii
   type(toml_table), pointer :: child
   character, parameter :: lsh(0:4) = ["s", "p", "d", "f", "g"]
   character, parameter :: pqn(1:7) = ["1", "2", "3", "4", "5", "6", "7"]

   do ii = 1, size(record)
      call add_table(table, trim(pqn(record(ii)%pqn)//lsh(record(ii)%lsh)), child)
      call record(ii)%dump(child, error)
      if (allocated(error)) exit
   end do
   if (allocated(error)) return
end subroutine records_to_table


!> Get the position of an element in the parameter records
pure subroutine get_element(record, sym, pos)
   !> Instance of the parametrization records
   type(element_record), intent(in) :: record(:)
   !> Symbol of the element
   character(len=*), intent(in) :: sym
   !> Position in the records
   integer, intent(out) :: pos

   integer :: num
   integer :: ii

   num = to_number(sym)
   pos = 0

   do ii = 1, size(record)
      if (record(ii)%sym == sym) then
         pos = ii
         exit
      end if
   end do
   if (pos /= 0) return

   do ii = 1, size(record)
      if (record(ii)%num == num) then
         pos = ii
         exit
      end if
   end do
   if (pos /= 0) return
end subroutine get_element

!> Read parametrization data from parameter array
subroutine load_from_array(self, array, offset, base, mask, error)
   class(element_record), intent(inout) :: self
   real(wp), intent(in) :: array(:)
   integer, intent(inout) :: offset
   type(element_record), intent(in) :: base
   type(element_mask), intent(in) :: mask
   type(error_type), allocatable, intent(out) :: error

   integer :: ii, ir

   select type(self)
   type is (element_record)
      self = base
   end select

   ! Increment
   call load_atom_par(self%increment, mask%increment, array, offset)

   ! Basis set
   call load_atom_par(self%k0, mask%k0, array, offset)
   call load_atom_par(self%k1, mask%k1, array, offset)
   call load_atom_par(self%k2, mask%k2, array, offset)
   call load_atom_par(self%k3, mask%k3, array, offset)
   call load_atom_par(self%k0_h0_scale, mask%k0_h0_scale, array, offset)
   call load_atom_par(self%k1_h0_scale, mask%k1_h0_scale, array, offset)
   call load_atom_par(self%k2_h0_scale, mask%k2_h0_scale, array, offset)
   call load_atom_par(self%k3_h0_scale, mask%k3_h0_scale, array, offset)

   ! Coordination number
   call load_atom_par(self%cn_rcov, mask%cn_rcov, array, offset)
   call load_atom_par(self%cn_avg, mask%cn_avg, array, offset)

   ! H0
   call load_atom_par(self%h0_rad, mask%h0_rad, array, offset)
   call load_atom_par(self%h0_dip_scale, mask%h0_dip_scale, array, offset)
   call load_atom_par(self%h0_diat_scale_sig, mask%h0_diat_scale_sig, array, offset)
   call load_atom_par(self%h0_diat_scale_pi, mask%h0_diat_scale_pi, array, offset)
   call load_atom_par(self%h0_diat_scale_del, mask%h0_diat_scale_del, array, offset)
   call load_atom_par(self%rvdw_scale, mask%rvdw_scale, array, offset)

   ! ACPs
   if (mask%n_acp > 0) then
      call load_acp_par(self%acp_expos, mask%acp_expos, array, offset)
      call load_acp_par(self%acp_levels, mask%acp_levels, array, offset)
   else 
      ! No ACPs are defined
      self%n_acp = 0
   end if

   ! Repulsion
   call load_atom_par(self%zeff, mask%zeff, array, offset)
   call load_atom_par(self%alpha, mask%alpha, array, offset)
   call load_atom_par(self%rep_cn, mask%rep_cn, array, offset)
   call load_atom_par(self%rep_q, mask%rep_q, array, offset)
   call load_atom_par(self%rep_roffset, mask%rep_roffset, array, offset)
   call load_atom_par(self%rep_k1, mask%rep_k1, array, offset)

   ! Halogen bonding correction
   call load_atom_par(self%xbond, mask%xbond, array, offset)

   ! First-order tight-binding
   call load_atom_par(self%ipea_cn, mask%ipea_cn, array, offset)

   ! Second-order tight-binding
   call load_atom_par(self%gam, mask%gam, array, offset)
   call load_atom_par(self%gam_cn, mask%gam_cn, array, offset)

   ! Third-order tight-bindingy
   call load_atom_par(self%gam3, mask%gam3, array, offset, scale=0.1_wp)

   ! Fourth-order tight-binding
   call load_atom_par(self%gam4, mask%gam4, array, offset, scale=0.1_wp)

   ! Multipoles
   call load_atom_par(self%dkernel, mask%dkernel, array, offset, scale=0.01_wp)
   call load_atom_par(self%qkernel, mask%qkernel, array, offset, scale=0.01_wp)
   call load_atom_par(self%aes_dip_scale, mask%aes_dip_scale, array, offset)

   ! Spin polarization
   call load_atom_par(self%wll_scale, mask%wll_scale, array, offset)

   ! Bond-order correlation
   call load_atom_par(self%cscale, mask%cscale, array, offset)
   call load_atom_par(self%crad, mask%crad, array, offset)

   ! Shells
   do ii = 1, size(mask%record)
      call get_shell(self%record, mask%record(ii)%pqn, mask%record(ii)%lsh, ir)
      call self%record(ir)%load(array, offset, base%record(ir), mask%record(ii), error)
      if (allocated(error)) exit
   end do
   if (allocated(error)) return
end subroutine load_from_array

!> Write parametrization data to parameter array
subroutine dump_to_array(self, array, offset, mask, error)
   class(element_record), intent(in) :: self
   real(wp), intent(inout) :: array(:)
   integer, intent(inout) :: offset
   type(element_mask), intent(in) :: mask
   type(error_type), allocatable, intent(out) :: error

   integer :: ii, ir

   ! Increment
   call dump_atom_par(self%increment, mask%increment, array, offset)
   
   ! Basis set
   call dump_atom_par(self%k0, mask%k0, array, offset)
   call dump_atom_par(self%k1, mask%k1, array, offset)
   call dump_atom_par(self%k2, mask%k2, array, offset)
   call dump_atom_par(self%k3, mask%k3, array, offset)
   call dump_atom_par(self%k0_h0_scale, mask%k0_h0_scale, array, offset)
   call dump_atom_par(self%k1_h0_scale, mask%k1_h0_scale, array, offset)
   call dump_atom_par(self%k2_h0_scale, mask%k2_h0_scale, array, offset)
   call dump_atom_par(self%k3_h0_scale, mask%k3_h0_scale, array, offset)
   
   ! Coordination number
   call dump_atom_par(self%cn_rcov, mask%cn_rcov, array, offset)
   call dump_atom_par(self%cn_avg, mask%cn_avg, array, offset)

   ! H0
   call dump_atom_par(self%h0_rad, mask%h0_rad, array, offset)
   call dump_atom_par(self%h0_dip_scale, mask%h0_dip_scale, array, offset)
   call dump_atom_par(self%h0_diat_scale_sig, mask%h0_diat_scale_sig, array, offset)
   call dump_atom_par(self%h0_diat_scale_pi, mask%h0_diat_scale_pi, array, offset)
   call dump_atom_par(self%h0_diat_scale_del, mask%h0_diat_scale_del, array, offset)
   call dump_atom_par(self%rvdw_scale, mask%rvdw_scale, array, offset)

   ! ACPs
   if (mask%n_acp > 0) then
      call dump_acp_par(self%acp_expos, mask%acp_expos, array, offset)
      call dump_acp_par(self%acp_levels, mask%acp_levels, array, offset)
   end if

   ! Repulsion
   call dump_atom_par(self%zeff, mask%zeff, array, offset)
   call dump_atom_par(self%alpha, mask%alpha, array, offset)
   call dump_atom_par(self%rep_cn, mask%rep_cn, array, offset)
   call dump_atom_par(self%rep_q, mask%rep_q, array, offset)
   call dump_atom_par(self%rep_roffset, mask%rep_roffset, array, offset)
   call dump_atom_par(self%rep_k1, mask%rep_k1, array, offset)

   ! Halogen bonding correction
   call dump_atom_par(self%xbond, mask%xbond, array, offset)

   ! First-order tight-binding
   call dump_atom_par(self%ipea_cn, mask%ipea_cn, array, offset)

   ! Second-order tight-binding
   call dump_atom_par(self%gam, mask%gam, array, offset)
   call dump_atom_par(self%gam_cn, mask%gam_cn, array, offset)

   ! Third-order tight-binding
   call dump_atom_par(self%gam3, mask%gam3, array, offset, scale=0.1_wp)

   ! Fourth-order tight-binding
   call dump_atom_par(self%gam4, mask%gam4, array, offset, scale=0.1_wp)

   ! Multipoles
   call dump_atom_par(self%dkernel, mask%dkernel, array, offset, scale=0.01_wp)
   call dump_atom_par(self%qkernel, mask%qkernel, array, offset, scale=0.01_wp)
   call dump_atom_par(self%aes_dip_scale, mask%aes_dip_scale, array, offset)

   ! Spin polarization
   call dump_atom_par(self%wll_scale, mask%wll_scale, array, offset)

   ! Bond-order correlation
   call dump_atom_par(self%cscale, mask%cscale, array, offset)
   call dump_atom_par(self%crad, mask%crad, array, offset)

   ! Shells
   do ii = 1, size(mask%record)
      call get_shell(self%record, mask%record(ii)%pqn, mask%record(ii)%lsh, ir)
      call self%record(ir)%dump(array, offset, mask%record(ii), error)
      if (allocated(error)) exit
   end do
   if (allocated(error)) return
end subroutine dump_to_array

pure subroutine load_atom_par(par, mask, array, ii, scale)
   real(wp), intent(inout) :: par
   logical, intent(in) :: mask
   real(wp), intent(in) :: array(:)
   integer, intent(inout) :: ii
   real(wp), intent(in), optional :: scale
   real(wp) :: scale_
   scale_ = 1.0_wp
   if (present(scale)) scale_ = scale
   if (mask) then
      ii = ii+1
      par = array(ii) * scale_
   end if
end subroutine load_atom_par

pure subroutine load_acp_par(par, mask, array, ii, scale)
   real(wp), intent(inout) :: par(:)
   logical, intent(in) :: mask(:)
   real(wp), intent(in) :: array(:)
   integer, intent(inout) :: ii
   real(wp), intent(in), optional :: scale
   real(wp) :: scale_
   integer :: iacp
   scale_ = 1.0_wp
   if (present(scale)) scale_ = scale
   do iacp = 1, size(par)
      if (mask(iacp)) then
         ii = ii+1
         par(iacp) = array(ii) * scale_
      end if
   end do
end subroutine load_acp_par

pure subroutine dump_atom_par(par, mask, array, ii, scale)
   real(wp), intent(in) :: par
   logical, intent(in) :: mask
   real(wp), intent(inout) :: array(:)
   integer, intent(inout) :: ii
   real(wp), intent(in), optional :: scale
   real(wp) :: scale_
   scale_ = 1.0_wp
   if (present(scale)) scale_ = scale
   if (mask) then
      ii = ii+1
      array(ii) = par / scale_
   end if
end subroutine dump_atom_par

pure subroutine dump_acp_par(par, mask, array, ii, scale)
   real(wp), intent(in) :: par(:)
   logical, intent(in) :: mask(:)
   real(wp), intent(inout) :: array(:)
   integer, intent(inout) :: ii
   real(wp), intent(in), optional :: scale
   real(wp) :: scale_
   integer :: iacp
   scale_ = 1.0_wp
   if (present(scale)) scale_ = scale
   do iacp = 1, size(par)
      if (mask(iacp)) then
         ii = ii+1
         array(ii) = par(iacp) / scale_
      end if
   end do
end subroutine dump_acp_par

elemental function count_mask(mask) result(ncount)
   type(element_mask), intent(in) :: mask
   integer :: ncount
   ncount = count([ &
      mask%increment, &
      mask%k0, &
      mask%k1, &
      mask%k2, &
      mask%k3, &
      mask%k0_h0_scale, &
      mask%k1_h0_scale, &
      mask%k2_h0_scale, &
      mask%k3_h0_scale, &
      mask%cn_rcov, &
      mask%cn_avg, &
      mask%h0_rad, &
      mask%h0_dip_scale, &
      mask%h0_diat_scale_sig, &
      mask%h0_diat_scale_pi, &
      mask%h0_diat_scale_del, &
      mask%rvdw_scale, &
      mask%zeff, &
      mask%alpha, &
      mask%rep_cn, &
      mask%rep_q, &
      mask%rep_roffset, &
      mask%rep_k1, &
      mask%xbond, &
      mask%ipea_cn, &
      mask%gam, &
      mask%gam_cn, &
      mask%gam3, &
      mask%gam4, &
      mask%dkernel, &
      mask%qkernel, &
      mask%aes_dip_scale, &
      mask%wll_scale, & 
      mask%cscale, &
      mask%crad])
   if (mask%n_acp > 0) then 
      ncount = ncount + count([ &
         mask%acp_expos, &
         mask%acp_levels])
   end if
   if (allocated(mask%record)) ncount = ncount + sum(count(mask%record))
end function count_mask

end module tblite_param_element
