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

!> @file tblite/param/shell.f90
!> Provides records for the shell specific parameters of one element

!> Definition of the shell specific parameter records
module tblite_param_shell
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io_symbols, only : to_number, symbol_length
   use tblite_basis, only : basis_set
   use tblite_param_serde, only : serde_record
   use tblite_toml, only : toml_table, toml_array, get_value, set_value, &
      & add_array, add_table, len
   implicit none
   private

   public :: count, get_shell


   !> The conversion factor from eV to Hartree is used for compatibility with older
   !> versions of xtb
   real(wp), parameter :: evtoau = 1.0_wp / 27.21138505_wp

   integer, parameter :: ngauss_default = 6
   real(wp), parameter :: kcn_default = 0.0_wp, h0_exp_scale_default = 1.0_wp, &
      & shpoly2_default = 0.0_wp, shpoly4_default = 0.0_wp, ipea_default=0.0_wp, &
      & zeffsh_default = 0.0_wp, lgam_default = 0.0_wp, lgam_fx_default=0.0_wp, &
      & avg_exp_fx_default = 1.0_wp

   character, parameter :: lsh(0:4) = ["s", "p", "d", "f", "g"]
   character, parameter :: pqn(1:7) = ["1", "2", "3", "4", "5", "6", "7"]

   character(len=*), parameter :: k_sto_ng = "sto_ng", k_q_vszp = "q_vszp", k_basis = "basis", &
      & k_ngauss = "ngauss", k_slater = "slater", k_expos = "expos", k_coeffs = "coeffs", &
      & k_coeffs_env = "coeffs_env", &
      & k_level = "level", k_kcn = "kcn", k_shpoly = "shpoly", k_shpoly2 = "shpoly2", &
      & k_shpoly4 = "shpoly4", k_h0_exp_scale = "h0_exp_scale", &
      & k_refocc = "refocc", k_zeffsh = "zeffsh", k_ipea = "ipea", k_lgam = "lgam", &
      & k_lgam_fx = "lgam_fx", k_avg_exp_fx = "avg_exp_fx"

   !> Representation of the shell specific parameters
   type, public, extends(serde_record) :: shell_record
      !> Angular momentum for the shell
      integer :: lsh
      !> Principal quantum number for the shell
      integer :: pqn

      !> Basis set type
      integer :: basis
      !> Number of primitive Gaussian functions used in the STO-NG expansion for the shell
      integer :: ngauss
      !> Slater exponents of the STO-NG functions for the shell
      real(wp) :: slater
      !> Exponents for the primitive Gaussian functions
      real(wp), allocatable :: expos(:)
      !> Coefficients for the primitive Gaussian functions
      real(wp), allocatable :: coeffs(:)
      !> Environment-dependence coefficient for the primitive Gaussian functions
      real(wp), allocatable :: coeffs_env(:)

      !> Atomic level energy for the shell
      real(wp) :: level
      !> CN dependent shift of the self energy for the shell
      real(wp) :: kcn
      !> Polynomial enhancement for Hamiltonian elements with square-root R-dependence
      real(wp) :: shpoly
      !> Polynomial enhancement for Hamiltonian elements with linear R-dependence
      real(wp) :: shpoly2
      !> Polynomial enhancement for Hamiltonian elements with square R-dependence
      real(wp) :: shpoly4
      !> Exponent scaling for the q-vSZP basis set in the H0 
      real(wp) :: h0_exp_scale

      !> Reference occupation for the shell
      real(wp) :: refocc

      !> Shell-distributed effective valence nuclear charge
      real(wp) :: zeffsh

      !> IP/EA for first-order tight-binding
      real(wp) :: ipea

      !> Relative chemical hardness for the shell
      real(wp) :: lgam

      !> Shell-scaling of Hubbard parameter for approximated Fock exchange
      real(wp) :: lgam_fx

      !> Averaging exponent for approximated Fock exchange
      real(wp) :: avg_exp_fx

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
   end type shell_record


   !> Masking for the shell record
   type, public :: shell_mask
      !> Angular momentum for the shell
      integer :: lsh
      !> Principal quantum number for the shell
      integer :: pqn

      !> Basis set type
      integer :: basis
      !> Number of primitive Gaussian functions used in the STO-NG expansion for the shell
      integer :: ngauss

      !> Slater exponents of the STO-NG functions for each shell
      logical :: slater
      !> Exponents for the primitive Gaussian functions
      logical, allocatable :: expos(:)
      !> Coefficients for the primitive Gaussian functions
      logical, allocatable :: coeffs(:)
      !> Environment-dependence coefficient for the primitive Gaussian functions
      logical, allocatable :: coeffs_env(:)
      
      !> Atomic level energies for each shell
      logical :: level
      !> CN dependent shift of the self energy for each shell
      logical :: kcn
      !> Polynomial enhancement for Hamiltonian elements with square-root R-dependence
      logical :: shpoly
      !> Polynomial enhancement for Hamiltonian elements with linear R-dependence
      logical :: shpoly2
      !> Polynomial enhancement for Hamiltonian elements with square R-dependence
      logical :: shpoly4
      !> Exponent scaling for the q-vSZP basis set in the H0 
      logical :: h0_exp_scale

      !> IP/EA for first-order tight-binding
      logical :: ipea

      !> Relative chemical hardness for each shell
      logical :: lgam

      !> Shell-scaling of Hubbard parameter for approximated Fock exchange
      logical :: lgam_fx

      !> Averaging exponent for approximated Fock exchange
      logical :: avg_exp_fx

   end type shell_mask


   interface count
      module procedure :: count_mask
   end interface count


contains


!> Read parametrization data from TOML data structure
subroutine load_from_toml(self, table, error)
   class(shell_record), intent(inout) :: self
   type(toml_table), intent(inout) :: table
   type(error_type), allocatable, intent(out) :: error

   character(len=:), allocatable :: ctmp
   integer :: stat

   ! Read the shell angular momentum and principal quantum number
   call table%get_key(ctmp)   

   self%pqn = get_pqn(ctmp(1:1))
   if (self%pqn < 1) then
      call fatal_error(error, "Invalid principal quantum number in "//trim(ctmp))
      return
   end if

   self%lsh = get_lsh(ctmp(2:2))
   if (self%lsh < 0) then
      call fatal_error(error, "Invalid angular momentum in "//trim(ctmp))
      return
   end if

   ! Basis set
   call get_basis(self, table, error)
   if (allocated(error)) return

   ! H0
   call get_h0(self, table, error)
   if (allocated(error)) return

   call get_value(table, k_refocc, self%refocc, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "No reference occupation specified for " &
         & // trim(pqn(self%pqn)//lsh(self%lsh)))
      return
   end if

   ! First-order tight-binding
   call get_value(table, k_zeffsh, self%zeffsh, zeffsh_default, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "No effective valence nuclear charge for " &
         & // trim(pqn(self%pqn)//lsh(self%lsh)))
      return
   end if

   call get_value(table, k_ipea, self%ipea, ipea_default, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Cannot read first order IP/EA for "&
         & // trim(pqn(self%pqn)//lsh(self%lsh)))
      return
   end if

   ! Second-order tight-binding
   call get_value(table, k_lgam, self%lgam, lgam_default, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Cannot read shell hardness for "&
         & // trim(pqn(self%pqn)//lsh(self%lsh)))
      return
   end if

   ! Approximated Fock exchange
   call get_value(table, k_lgam_fx, self%lgam_fx, lgam_fx_default, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Cannot read shell hardness for Fock exchange for " &
         & // trim(pqn(self%pqn)//lsh(self%lsh)))
      return
   end if

   call get_value(table, k_avg_exp_fx, self%avg_exp_fx, avg_exp_fx_default, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Cannot read shell averaging exponent for Fock exchange for " &
         & // trim(pqn(self%pqn)//lsh(self%lsh)))
      return
   end if

end subroutine load_from_toml


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


subroutine get_basis(self, table, error)
   !> Instance of the parametrization data
   type(shell_record), intent(inout) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table), pointer :: child
   type(toml_array), pointer :: array
   integer :: stat
   
   character(len=:), allocatable :: ctmp

   self%basis = basis_set%sto_ng
   call get_value(table, k_sto_ng, child, requested=.false.)
   if (.not. associated(child)) then
      self%basis = basis_set%q_vszp
      call get_value(table, k_q_vszp, child, requested=.false.)
   end if
   if (.not. associated(child)) then
      self%basis = basis_set%custom
      call get_value(table, k_basis, child, requested=.false.)
   end if
   if (.not. associated(child)) then
      call fatal_error(error, "No basis set specified for " & 
         & // trim(pqn(self%pqn)//lsh(self%lsh)))
      return
   end if

   if (self%basis == basis_set%sto_ng) then
      ! Stewart STO-NG basis set
      call get_value(child, k_slater, self%slater, stat=stat)
      if (stat /= 0 .or. self%slater < epsilon(0.0_wp)) then
         call fatal_error(error, "Invalid slater exponents specified for " &
            & // trim(pqn(self%pqn)//lsh(self%lsh))) 
         return
      end if

      call get_value(child, k_ngauss, self%ngauss, ngauss_default, stat=stat)
      if (stat /= 0 .or. self%ngauss < 1 .or. self%ngauss > 6) then
         call fatal_error(error, "Invalid number of primitive gaussian for " &
            & // trim(pqn(self%pqn)//lsh(self%lsh)))
         return
      end if
   
   else if (self%basis == basis_set%q_vszp) then
      ! Charge-dependent q-vSZP basis set
      call get_value(child, k_ngauss, self%ngauss, ngauss_default, stat=stat)
      if (stat /= 0 .or. self%ngauss < 1) then
         call fatal_error(error, "Invalid number of primitive gaussian for " & 
            & // trim(pqn(self%pqn)//lsh(self%lsh)))
         return
      end if

      call get_value(child, k_expos, array)
      if (.not.associated(array)) then
         call fatal_error(error, "No exponents for primitive gaussians specified for " &
            & // trim(pqn(self%pqn)//lsh(self%lsh)))
         return
      end if
      call get_exponents(self, array, error)
      if (allocated(error)) return

      call get_value(child, k_coeffs, array)
      if (.not.associated(array)) then
         call fatal_error(error, "No coefficients for primitive gaussians specified for " &
            & // trim(pqn(self%pqn)//lsh(self%lsh))) 
         return
      end if
      call get_coefficients(self, array, error)
      if (allocated(error)) return
      
      call get_value(child, k_coeffs_env, array)
      if (.not.associated(array)) then
         call fatal_error(error, "No environment coefficients for primitive gaussians specified for " &
            & // trim(pqn(self%pqn)//lsh(self%lsh)))
         return
      end if
      call get_coefficients_env(self, array, error)
      if (allocated(error)) return
   
   else 
      ! Custom basis set
      call get_value(child, k_ngauss, self%ngauss, ngauss_default, stat=stat)
      if (stat /= 0 .or. self%ngauss < 1) then
         call fatal_error(error, "Invalid number of primitive gaussian for "&
            & // trim(pqn(self%pqn)//lsh(self%lsh)))
         return
      end if

      call get_value(child, k_expos, array)
      if (.not.associated(array)) then
         call fatal_error(error, "No exponents for primitive gaussians specified for " &
            & // trim(pqn(self%pqn)//lsh(self%lsh)))
         return
      end if
      call get_exponents(self, array, error)
      if (allocated(error)) return

      call get_value(child, k_coeffs, array)
      if (.not.associated(array)) then
         call fatal_error(error, "No coefficients for primitive gaussians specified for " &
            & // trim(pqn(self%pqn)//lsh(self%lsh)))
         return
      end if
      call get_coefficients(self, array, error)
      if (allocated(error)) return

   end if 

end subroutine get_basis

subroutine get_exponents(self, array, error)
   type(shell_record), intent(inout) :: self
   type(toml_array), intent(inout) :: array
   type(error_type), allocatable, intent(out) :: error

   integer :: i
   real(wp) :: tmp

   if (len(array) < self%ngauss) then
      call fatal_error(error, "Insufficient primitive gaussian exponents for " &
         & // trim(pqn(self%pqn)//lsh(self%lsh)))
      return
   end if
   allocate(self%expos(self%ngauss))

   do i = 1, self%ngauss
      call get_value(array, i, tmp)
      self%expos(i) = tmp
   end do

   if (any(self%expos < epsilon(0.0_wp))) then
      call fatal_error(error, "Invalid primitive gaussian exponents for " &
         & // trim(pqn(self%pqn)//lsh(self%lsh)))
   end if
end subroutine get_exponents

subroutine get_coefficients(self, array, error)
   type(shell_record), intent(inout) :: self
   type(toml_array), intent(inout) :: array
   type(error_type), allocatable, intent(out) :: error

   integer :: i
   real(wp) :: tmp

   if (len(array) < self%ngauss) then
      call fatal_error(error, "Insufficient primitive gaussian coefficients for " &
         & // trim(pqn(self%pqn)//lsh(self%lsh)))
      return
   end if
   allocate(self%coeffs(self%ngauss))

   do i = 1, self%ngauss
      call get_value(array, i, tmp)
      self%coeffs(i) = tmp
   end do
end subroutine get_coefficients

subroutine get_coefficients_env(self, array, error)
   type(shell_record), intent(inout) :: self
   type(toml_array), intent(inout) :: array
   type(error_type), allocatable, intent(out) :: error

   integer :: i
   real(wp) :: tmp

   if (len(array) < self%ngauss) then
      call fatal_error(error, "Insufficient primitive gaussian environment coefficients for " &
         & // trim(pqn(self%pqn)//lsh(self%lsh)))
      return
   end if
   allocate(self%coeffs_env(self%ngauss))

   do i = 1, self%ngauss
      call get_value(array, i, tmp)
      self%coeffs_env(i) = tmp
   end do
end subroutine get_coefficients_env


subroutine get_h0(self, table, error)
   !> Instance of the parametrization data
   type(shell_record), intent(inout) :: self
   !> Data structure
   type(toml_table), intent(inout) :: table
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: stat
   real(wp) :: rtmp

   call get_value(table, k_level, rtmp, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "No atomic levels specified for " &
         & // trim(pqn(self%pqn)//lsh(self%lsh)))
      return
   end if
   self%level = rtmp * evtoau

   call get_value(table, k_kcn, rtmp, kcn_default, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Cannot read CN shift parameter for " &
         & // trim(pqn(self%pqn)//lsh(self%lsh)))
      return
   end if
   self%kcn = rtmp * evtoau

   call get_value(table, k_shpoly, self%shpoly, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "No square-root shell-polynomials specified for " &
         & // trim(pqn(self%pqn)//lsh(self%lsh)))
      return
   end if

   call get_value(table, k_shpoly2, self%shpoly2, shpoly2_default, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "No linear shell-polynomials specified for " &
         & // trim(pqn(self%pqn)//lsh(self%lsh)))
      return
   end if

   call get_value(table, k_shpoly4, self%shpoly4, shpoly4_default, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "No square shell-polynomials specified for " &
         & // trim(pqn(self%pqn)//lsh(self%lsh)))
      return
   end if

   call get_value(table, k_h0_exp_scale, self%h0_exp_scale, h0_exp_scale_default, stat=stat)
   if (stat /= 0) then
      call fatal_error(error, "Cannot read exponent H0 scaling for " &
         & // trim(pqn(self%pqn)//lsh(self%lsh)))
      return
   end if

end subroutine get_h0


!> Write parametrization data to TOML data structure
subroutine dump_to_toml(self, table, error)
   class(shell_record), intent(in) :: self
   type(toml_table), intent(inout) :: table
   type(error_type), allocatable, intent(out) :: error

   integer :: i
   type(toml_array), pointer :: array
   type(toml_table), pointer :: child

   ! Basis set
   if (self%basis == basis_set%sto_ng) then
      ! Stewart STO-NG basis set
      call add_table(table, k_sto_ng, child)      

      call set_value(child, k_slater, self%slater)
      call set_value(child, k_ngauss, self%ngauss)

   else if (self%basis == basis_set%q_vszp) then
      ! Charge-dependent q-vSZP basis set
      call add_table(table, k_q_vszp, child)

      call set_value(child, k_ngauss, self%ngauss)

      call add_array(child, k_expos, array)
      do i = 1, self%ngauss
         call set_value(array, i, self%expos(i))
      end do   

      call add_array(child, k_coeffs, array)
      do i = 1, self%ngauss
         call set_value(array, i, self%coeffs(i))
      end do   

      call add_array(child, k_coeffs_env, array)
      do i = 1, self%ngauss
         call set_value(array, i, self%coeffs_env(i))
      end do

   else
      ! Custom basis set
      call add_table(table, k_basis, child)      
      
      call set_value(child, k_ngauss, self%ngauss)

      call add_array(child, k_expos, array)
      do i = 1, self%ngauss
         call set_value(array, i, self%expos(i))
      end do   

      call add_array(child, k_coeffs, array)
      do i = 1, self%ngauss
         call set_value(array, i, self%coeffs(i))
      end do   

   end if

   ! H0
   call set_value(table, k_level, self%level / evtoau)
   call set_value(table, k_kcn, self%kcn / evtoau)
   call set_value(table, k_shpoly, self%shpoly)
   call set_value(table, k_shpoly2, self%shpoly2)
   call set_value(table, k_shpoly4, self%shpoly4)
   call set_value(table, k_h0_exp_scale, self%h0_exp_scale)

   call set_value(table, k_refocc, self%refocc)
   
   ! First-order tight-binding
   call set_value(table, k_zeffsh, self%zeffsh)
   call set_value(table, k_ipea, self%ipea)
   
   ! Second-order tight-binding
   call set_value(table, k_lgam, self%lgam)
   
   ! Approximated Fock exchange
   call set_value(table, k_lgam_fx, self%lgam_fx)
   call set_value(table, k_avg_exp_fx, self%avg_exp_fx)

end subroutine dump_to_toml


!> Get the position of a shell in the element records
pure subroutine get_shell(record, pqn, lsh, pos)
   !> Instance of the parametrization records
   class(shell_record), intent(in) :: record(:)
   !> Principle quantum number of the shell
   integer, intent(in) :: pqn
   !> Angular momentum of the shell
   integer, intent(in) :: lsh
   !> Position in the records
   integer, intent(out) :: pos

   integer :: ii

   pos = 0
   do ii = 1, size(record)
      if (record(ii)%pqn == pqn .and. record(ii)%lsh == lsh) then
         pos = ii
         exit
      end if
   end do
end subroutine get_shell


!> Read parametrization data from parameter array
subroutine load_from_array(self, array, offset, base, mask, error)
   class(shell_record), intent(inout) :: self
   real(wp), intent(in) :: array(:)
   integer, intent(inout) :: offset
   type(shell_record), intent(in) :: base
   type(shell_mask), intent(in) :: mask
   type(error_type), allocatable, intent(out) :: error

   select type(self)
   type is (shell_record)
      self = base
   end select

   ! Basis set
   if (self%basis == basis_set%sto_ng) then
      ! Stewart STO-NG basis set
      call load_shell_par(self%slater, mask%slater, array, offset)
   else if (self%basis == basis_set%q_vszp) then
      ! Charge-dependent q-vSZP basis set
      call load_prim_par(self%expos, mask%expos, array, offset)
      call load_prim_par(self%coeffs, mask%coeffs, array, offset)
      call load_prim_par(self%coeffs_env, mask%coeffs_env, array, offset)
   else 
      ! Custom basis set
      call load_prim_par(self%expos, mask%expos, array, offset)
      call load_prim_par(self%coeffs, mask%coeffs, array, offset)
   end if

   ! H0
   call load_shell_par(self%level, mask%level, array, offset, scale=evtoau)
   call load_shell_par(self%kcn, mask%kcn, array, offset, scale=evtoau)
   call load_shell_par(self%shpoly, mask%shpoly, array, offset, scale=0.01_wp)
   call load_shell_par(self%shpoly2, mask%shpoly2, array, offset, scale=0.01_wp)
   call load_shell_par(self%shpoly4, mask%shpoly4, array, offset, scale=0.01_wp)
   call load_shell_par(self%h0_exp_scale, mask%h0_exp_scale, array, offset)

   ! First-order tight-binding
   call load_shell_par(self%ipea, mask%ipea, array, offset)

   ! Second-order tight-binding
   call load_shell_par(self%lgam, mask%lgam, array, offset)

   ! Approximated Fock exchange
   call load_shell_par(self%lgam_fx, mask%lgam_fx, array, offset)
   call load_shell_par(self%avg_exp_fx, mask%avg_exp_fx, array, offset)

end subroutine load_from_array

!> Write parametrization data to parameter array
subroutine dump_to_array(self, array, offset, mask, error)
   class(shell_record), intent(in) :: self
   real(wp), intent(inout) :: array(:)
   integer, intent(inout) :: offset
   type(shell_mask), intent(in) :: mask
   type(error_type), allocatable, intent(out) :: error

   ! Basis set
   if (self%basis == basis_set%sto_ng) then
      ! Stewart STO-NG basis set
      call dump_shell_par(self%slater, mask%slater, array, offset)
   else if (self%basis == basis_set%q_vszp) then
      ! Charge-dependent q-vSZP basis set
      call dump_prim_par(self%expos, mask%expos, array, offset)
      call dump_prim_par(self%coeffs, mask%coeffs, array, offset)
      call dump_prim_par(self%coeffs_env, mask%coeffs_env, array, offset)
   else 
      ! Custom basis set
      call dump_prim_par(self%expos, mask%expos, array, offset)
      call dump_prim_par(self%coeffs, mask%coeffs, array, offset)
   end if

   ! H0
   call dump_shell_par(self%level, mask%level, array, offset, scale=evtoau)
   call dump_shell_par(self%kcn, mask%kcn, array, offset, scale=evtoau)
   call dump_shell_par(self%shpoly, mask%shpoly, array, offset, scale=0.01_wp)
   call dump_shell_par(self%shpoly2, mask%shpoly2, array, offset, scale=0.01_wp)
   call dump_shell_par(self%shpoly4, mask%shpoly4, array, offset, scale=0.01_wp)
   call dump_shell_par(self%h0_exp_scale, mask%h0_exp_scale, array, offset)

   ! First-order tight-binding
   call dump_shell_par(self%ipea, mask%ipea, array, offset)

   ! Second-order tight-binding
   call dump_shell_par(self%lgam, mask%lgam, array, offset)

   ! Approximated Fock exchange
   call dump_shell_par(self%lgam_fx, mask%lgam_fx, array, offset)
   call dump_shell_par(self%avg_exp_fx, mask%avg_exp_fx, array, offset)

end subroutine dump_to_array

pure subroutine load_shell_par(par, mask, array, ii, scale)
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
end subroutine load_shell_par

pure subroutine load_prim_par(par, mask, array, ii, scale)
   real(wp), intent(inout) :: par(:)
   logical, intent(in) :: mask(:)
   real(wp), intent(in) :: array(:)
   integer, intent(inout) :: ii
   real(wp), intent(in), optional :: scale
   real(wp) :: scale_
   integer :: iprim
   scale_ = 1.0_wp
   if (present(scale)) scale_ = scale
   do iprim = 1, size(par)
      if (mask(iprim)) then
         ii = ii+1
         par(iprim) = array(ii) * scale_
      end if
   end do
end subroutine load_prim_par

pure subroutine dump_shell_par(par, mask, array, ii, scale)
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
end subroutine dump_shell_par

pure subroutine dump_prim_par(par, mask, array, ii, scale)
   real(wp), intent(in) :: par(:)
   logical, intent(in) :: mask(:)
   real(wp), intent(inout) :: array(:)
   integer, intent(inout) :: ii
   real(wp), intent(in), optional :: scale
   real(wp) :: scale_
   integer :: iprim
   scale_ = 1.0_wp
   if (present(scale)) scale_ = scale
   do iprim = 1, size(par)
      if (mask(iprim)) then
         ii = ii+1
         array(ii) = par(iprim) / scale_
      end if
   end do
end subroutine dump_prim_par

elemental function count_mask(mask) result(ncount)
   type(shell_mask), intent(in) :: mask
   integer :: ncount
   ncount = count([ &
      mask%level, &
      mask%kcn, &
      mask%shpoly, &
      mask%shpoly2, &
      mask%shpoly4, &
      mask%h0_exp_scale, &
      mask%ipea, &
      mask%lgam, &
      mask%lgam_fx, &
      mask%avg_exp_fx &
   ])
   if (mask%basis == basis_set%sto_ng) then
      ncount = ncount + count([ &
         mask%slater &
      ])
   else if (mask%basis == basis_set%q_vszp) then
      ncount = ncount + count([ &
         mask%expos, &
         mask%coeffs, &
         mask%coeffs_env &
      ])
   else
      ncount = ncount + count([ &
         mask%expos, &
         mask%coeffs &
      ])
   end if
end function count_mask

end module tblite_param_shell
