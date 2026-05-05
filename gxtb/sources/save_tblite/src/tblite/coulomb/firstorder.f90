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

!> @file tblite/coulomb/firstorder.f90
!> Provides an onsite first-order tight-binding interaction

!> Isotropic first-order onsite correction
module tblite_coulomb_firstorder
   use mctc_env, only : error_type, wp
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use tblite_container_cache, only : container_cache
   use tblite_coulomb_cache, only : coulomb_cache
   use tblite_coulomb_type, only : coulomb_type
   use tblite_scf_potential, only : potential_type
   use tblite_wavefunction_type, only : wavefunction_type
   implicit none
   private

   public :: new_onsite_firstorder

   !> Onsite correction for first-order charge expansion
   type, public, extends(coulomb_type) :: onsite_firstorder
      !> Whether the first order contribution is shell-dependent
      logical :: shell_resolved
      !> Number of shell for each atom
      integer, allocatable :: nsh_at(:)
      !> Shell offset for each atom
      integer, allocatable :: ish_at(:)
      !> Exponent of IP and EA spliting function
      real(wp) :: split_exp
      !> Slope of IP and EA spliting function
      real(wp) :: split_slope
      !> Offset of IP and EA spliting function
      real(wp) :: split_offset
      !> IP/EA chemical potential for each species and shell
      real(wp), allocatable :: ipea(:, :)
      !> CN-dependence of IP/EA for each species
      real(wp), allocatable :: ipea_cn(:)
   contains
      !> Update container cache
      procedure :: update
      !> Get information about density dependent quantities used in the energy
      procedure :: variable_info
      !> Evaluate selfconsistent energy of the interaction
      procedure :: get_energy
      !> Evaluate charge dependent potential shift from the interaction
      procedure :: get_potential
      !> Evaluate gradient of charge dependent potential shift from the interaction
      procedure :: get_potential_gradient
      !> Evaluate gradient contributions from the selfconsistent interaction
      procedure :: get_gradient
   end type onsite_firstorder

   real(wp), parameter :: tosp = 2.0_wp / sqrt(pi)
   real(wp), parameter :: fosp = 2.0_wp * tosp

   character(len=*), parameter :: label = "onsite first-order tight-binding"

   real(wp), parameter :: cn_reg = 1e-6_wp

contains


!> Create new onsite first-order contribution
subroutine new_onsite_firstorder(self, mol, ipea, ipea_cn, split_exp, &
   & split_slope, split_offset, nshell)
   !> Instance of the first-order tight-binding container
   type(onsite_firstorder), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> IP/EA chemical potential for each species and shell
   real(wp), intent(in) :: ipea(:, :)
   !> CN-dependence of IP/EA for each species and shell
   real(wp), intent(in) :: ipea_cn(:)
   !> Exponent of IP and EA spliting function
   real(wp), intent(in) :: split_exp
   !> Slope of IP and EA spliting function
   real(wp), intent(in) :: split_slope
   !> Offset of IP and EA spliting function
   real(wp), intent(in) :: split_offset
   !> Number of shells for each species
   integer, intent(in), optional :: nshell(:)

   integer :: ind, iat

   self%label = label

   self%ipea = ipea
   self%ipea_cn = ipea_cn

   self%split_exp = split_exp
   self%split_slope = split_slope
   self%split_offset = split_offset

   self%shell_resolved = present(nshell)
   if (present(nshell)) then
      self%nsh_at = nshell(mol%id)

      allocate(self%ish_at(mol%nat))
      ind = 0
      do iat = 1, mol%nat
         self%ish_at(iat) = ind
         ind = ind + self%nsh_at(iat)
      end do
   end if

end subroutine new_onsite_firstorder


!> Update container cache
subroutine update(self, mol, cache)
   !> Instance of the first-order tight-binding container
   class(onsite_firstorder), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache

end subroutine update


!> Evaluate selfconsistent energy of the interaction
subroutine get_energy(self, mol, cache, wfn, energies)
   !> Instance of the first-order tight-binding container
   class(onsite_firstorder), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Electrostatics and tight-binding energy
   real(wp), intent(inout) :: energies(:)

   integer :: iat, izp, ii, ish
   real(wp) :: scale, erf_ip, erf_ea, split, sqrtcni
   type(coulomb_cache), pointer :: ptr

   call view(cache, ptr)

   if (self%shell_resolved) then
      do iat = 1, mol%nat
         izp = mol%id(iat)
         ii = self%ish_at(iat)
         ! Coordination number scaling with regularization
         sqrtcni = sqrt(ptr%cn(iat) + cn_reg**2) - cn_reg
         scale = 1.0_wp + self%ipea_cn(izp) * sqrtcni
         ! IP and EA splitting error function
         erf_ip = erf(self%split_exp * (wfn%qat(iat, 1) - self%split_offset))
         erf_ea = erf(self%split_exp * (wfn%qat(iat, 1) + self%split_offset))
         ! Splitting function 
         split = (2.0_wp + self%split_slope * (erf_ip + erf_ea))

         do ish = 1, self%nsh_at(iat)
            energies(iat) = energies(iat) + 0.5_wp * self%ipea(ish, izp) &
               & * scale * split * wfn%qsh(ii+ish, 1)
         end do
      end do
   else
      do iat = 1, mol%nat
         izp = mol%id(iat)
         ! Coordination number scaling with regularization
         sqrtcni = sqrt(ptr%cn(iat) + cn_reg**2) - cn_reg
         scale = 1.0_wp + self%ipea_cn(izp) * sqrtcni
         ! IP and EA splitting error function
         erf_ip = erf(self%split_exp * (wfn%qat(iat, 1) - self%split_offset))
         erf_ea = erf(self%split_exp * (wfn%qat(iat, 1) + self%split_offset))
         ! Splitting function 
         split = (2.0_wp + self%split_slope * (erf_ip + erf_ea))

         energies(iat) = energies(iat) + 0.5_wp * self%ipea(1, izp) &
            & * scale * split * wfn%qat(iat, 1)
      end do
   end if

end subroutine get_energy


!> Evaluate charge dependent potential shift from the interaction
subroutine get_potential(self, mol, cache, wfn, pot)
   !> Instance of the first-order tight-binding container
   class(onsite_firstorder), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Density dependent potential
   type(potential_type), intent(inout) :: pot

   integer :: iat, izp, ii, ish
   real(wp) :: scale, erf_ip, erf_ea, derf_ip, derf_ea, split, dsplitdq, sqrtcni
   type(coulomb_cache), pointer :: ptr

   call view(cache, ptr)

   if (self%shell_resolved) then
      do iat = 1, mol%nat
         izp = mol%id(iat)
         ii = self%ish_at(iat)
         ! Coordination number scaling with regularization
         sqrtcni = sqrt(ptr%cn(iat) + cn_reg**2) - cn_reg
         scale = 1.0_wp + self%ipea_cn(izp) * sqrtcni
         ! IP and EA splitting error function and derivative
         erf_ip = erf(self%split_exp * (wfn%qat(iat, 1) - self%split_offset))
         derf_ip = tosp * self%split_exp * exp(-(self%split_exp &
            & * (wfn%qat(iat, 1) - self%split_offset))**2)
         erf_ea = erf(self%split_exp * (wfn%qat(iat, 1) + self%split_offset))
         derf_ea = tosp * self%split_exp * exp(-(self%split_exp &
            & * (wfn%qat(iat, 1) + self%split_offset))**2)
         ! Splitting function and derivative
         split = (2.0_wp + self%split_slope * (erf_ip + erf_ea))
         dsplitdq = self%split_slope * (derf_ip + derf_ea)

         do ish = 1, self%nsh_at(iat)
            ! Shell-resolved potential due to first-order charge
            pot%vsh(ii+ish, 1) = pot%vsh(ii+ish, 1) + 0.5_wp &
               & * self%ipea(ish, izp) * scale * split
            ! Atom-resolved potential due to IP/EA splitting
            pot%vat(iat, 1) = pot%vat(iat, 1) + 0.5_wp &
               & * self%ipea(ish, izp) * scale * dsplitdq * wfn%qsh(ii+ish, 1)
         end do
      end do
   else
      do iat = 1, mol%nat
         izp = mol%id(iat)
         ! Coordination number scaling with regularization
         sqrtcni = sqrt(ptr%cn(iat) + cn_reg**2) - cn_reg
         scale = 1.0_wp + self%ipea_cn(izp) * sqrtcni
         ! IP and EA splitting error function and derivative
         erf_ip = erf(self%split_exp * (wfn%qat(iat, 1) - self%split_offset))
         derf_ip = tosp * self%split_exp * exp(-(self%split_exp &
            & * (wfn%qat(iat, 1) - self%split_offset))**2)
         erf_ea = erf(self%split_exp * (wfn%qat(iat, 1) + self%split_offset))
         derf_ea = tosp * self%split_exp * exp(-(self%split_exp &
            & * (wfn%qat(iat, 1) + self%split_offset))**2)
         ! Splitting function and derivative
         split = (2.0_wp + self%split_slope * (erf_ip + erf_ea))
         dsplitdq = self%split_slope * (derf_ip + derf_ea)

         ! Atom-resolved potential due to first-order charge
         pot%vat(iat, 1) = pot%vat(iat, 1) + 0.5_wp * self%ipea(1, izp) &
            & * scale * (split + dsplitdq * wfn%qat(iat, 1) )
      end do
   end if

end subroutine get_potential


!> Evaluate gradient of charge dependent potential shift from the interaction
subroutine get_potential_gradient(self, mol, cache, wfn, pot)
   !> Instance of the first-order tight-binding container
   class(onsite_firstorder), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Density dependent potential
   type(potential_type), intent(inout) :: pot

   integer :: iat, izp, ii, ish
   real(wp) :: scale, dscaledcn, erf_ip, erf_ea, derf_ip, derf_ea, dderf_ip, &
      & dderf_ea, dsplitdq, dsplitdcn, dsplitdqdcn, dsplitdqdq, sqrtcni
   type(coulomb_cache), pointer :: ptr

   call view(cache, ptr)

   if (self%shell_resolved) then
      do iat = 1, mol%nat
         izp = mol%id(iat)
         ii = self%ish_at(iat)
         ! Coordination number scaling and derivative with regularization
         sqrtcni = sqrt(ptr%cn(iat) + cn_reg**2) - cn_reg
         scale = 1.0_wp + self%ipea_cn(izp) * sqrtcni
         dscaledcn = self%ipea_cn(izp) / (2.0_wp * (sqrtcni + cn_reg))
         ! IP and EA splitting error function and derivatives
         erf_ip = erf(self%split_exp * (wfn%qat(iat, 1) - self%split_offset))
         derf_ip = tosp * self%split_exp * exp(-(self%split_exp &
            & * (wfn%qat(iat, 1) - self%split_offset))**2)
         dderf_ip = -fosp * self%split_exp**3 &
            & * (wfn%qat(iat, 1) - self%split_offset) &
            & * exp(-(self%split_exp * (wfn%qat(iat, 1) - self%split_offset))**2)
         erf_ea = erf(self%split_exp * (wfn%qat(iat, 1) + self%split_offset))
         derf_ea = tosp * self%split_exp * exp(-(self%split_exp &
            & * (wfn%qat(iat, 1) + self%split_offset))**2)
         dderf_ea = -fosp * self%split_exp**3 &
            & * (wfn%qat(iat, 1) + self%split_offset) &
            & * exp(-(self%split_exp * (wfn%qat(iat, 1) + self%split_offset))**2)
         ! Splitting function and derivatives
         dsplitdcn = 2.0_wp + self%split_slope * (erf_ip + erf_ea)
         dsplitdq = self%split_slope * (derf_ip + derf_ea)
         dsplitdqdcn = self%split_slope * (derf_ip + derf_ea)
         dsplitdqdq = self%split_slope * (dderf_ip + dderf_ea)

         do ish = 1, self%nsh_at(iat)
            ! Derivatives of the shell-resolved potential (CN and qat)
            pot%dvshdr(:, :, ii+ish, 1) = pot%dvshdr(:, :, ii+ish, 1) + 0.5_wp &
               & * self%ipea(ish, izp) * (dscaledcn * dsplitdcn * ptr%dcndr(:, :, iat) &
               & + scale * dsplitdq * wfn%dqatdr(:, :, iat, 1))
            pot%dvshdL(:, :, ii+ish, 1) = pot%dvshdL(:, :, ii+ish, 1) + 0.5_wp &
               & * self%ipea(ish, izp) * (dscaledcn * dsplitdcn * ptr%dcndL(:, :, iat) &
               & + scale * dsplitdq * wfn%dqatdL(:, :, iat, 1))

            ! Derivatives of the atom-resolved potential 
            ! (CN and qat for the splitting and qsh for the first-order charge)
            pot%dvatdr(:, :, iat, 1) = pot%dvatdr(:, :, iat, 1) + 0.5_wp &
               & * self%ipea(ish, izp) * (wfn%qsh(ii+ish, 1) * dscaledcn * dsplitdqdcn &
               & * ptr%dcndr(:, :, iat) + scale * dsplitdq * wfn%dqshdr(:, :, ii+ish, 1) &
               & + scale * dsplitdqdq * wfn%qsh(ii+ish, 1) * wfn%dqatdr(:, :, iat, 1))
            pot%dvatdL(:, :, iat, 1) = pot%dvatdL(:, :, iat, 1) + 0.5_wp &
               & * self%ipea(ish, izp) * (wfn%qsh(ii+ish, 1) * dscaledcn * dsplitdqdcn &
               & * ptr%dcndL(:, :, iat) + scale * dsplitdq * wfn%dqshdL(:, :, ii+ish, 1) &
               & + scale * dsplitdqdq * wfn%qsh(ii+ish, 1) * wfn%dqatdL(:, :, iat, 1))
         end do
      end do
   else
      do iat = 1, mol%nat
         izp = mol%id(iat)
         ! Coordination number scaling and derivative with regularization
         sqrtcni = sqrt(ptr%cn(iat) + cn_reg**2) - cn_reg
         scale = 1.0_wp + self%ipea_cn(izp) * sqrtcni
         dscaledcn = self%ipea_cn(izp) / (2.0_wp * (sqrtcni + cn_reg))
         ! IP and EA splitting error function and derivatives
         erf_ip = erf(self%split_exp * (wfn%qat(iat, 1) - self%split_offset))
         derf_ip = tosp * self%split_exp * exp(-(self%split_exp &
            & * (wfn%qat(iat, 1) - self%split_offset))**2)
         dderf_ip = -fosp * self%split_exp**3 &
            & * (wfn%qat(iat, 1) - self%split_offset) &
            & * exp(-(self%split_exp * (wfn%qat(iat, 1) - self%split_offset))**2)
         erf_ea = erf(self%split_exp * (wfn%qat(iat, 1) + self%split_offset))
         derf_ea = tosp * self%split_exp * exp(-(self%split_exp &
            & * (wfn%qat(iat, 1) + self%split_offset))**2)
         dderf_ea = -fosp * self%split_exp**3 &
            & * (wfn%qat(iat, 1) + self%split_offset) &
            & * exp(-(self%split_exp * (wfn%qat(iat, 1) + self%split_offset))**2)
         ! Splitting function and derivatives
         dsplitdcn = 2.0_wp + self%split_slope * (erf_ip + erf_ea)
         dsplitdq = self%split_slope * (derf_ip + derf_ea)
         dsplitdqdcn = self%split_slope * (derf_ip + derf_ea)
         dsplitdqdq = self%split_slope * (dderf_ip + dderf_ea)

         ! Derivatives of the atom-resolved potential 
         ! (CN and qat for the splitting and qat for the first-order charge)
         pot%dvatdr(:, :, iat, 1) = pot%dvatdr(:, :, iat, 1) + 0.5_wp &
            & * self%ipea(1, izp) * (dscaledcn * (dsplitdcn + wfn%qat(iat, 1) &
            & * dsplitdqdcn) * ptr%dcndr(:, :, iat) + scale * (2.0_wp * dsplitdq &
            & + dsplitdqdq * wfn%qat(iat, 1)) * wfn%dqatdr(:, :, iat, 1))
         pot%dvatdL(:, :, iat, 1) = pot%dvatdL(:, :, iat, 1) + 0.5_wp &
            & * self%ipea(1, izp) * (dscaledcn *  (dsplitdcn + wfn%qat(iat, 1) &
            & * dsplitdqdcn) * ptr%dcndL(:, :, iat) + scale * (2.0_wp * dsplitdq &
            & + dsplitdqdq * wfn%qat(iat, 1)) * wfn%dqatdL(:, :, iat, 1))
      end do
   end if

end subroutine get_potential_gradient


!> Evaluate gradient contributions from the selfconsistent interaction
subroutine get_gradient(self, mol, cache, wfn, gradient, sigma)
   !> Instance of the first-order tight-binding container
   class(onsite_firstorder), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Molecular gradient of the first-order tight-binding energy
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Strain derivatives of the first-order tight-binding energy
   real(wp), contiguous, intent(inout) :: sigma(:, :)

   integer :: iat, izp, ii, ish
   real(wp) :: dscaledcn, erf_ip, erf_ea, dsplitdcn
   type(coulomb_cache), pointer :: ptr

   call view(cache, ptr)

   if (self%shell_resolved) then
      do iat = 1, mol%nat
         izp = mol%id(iat)
         ii = self%ish_at(iat)
         ! Coordination number scaling derivative with regularization
         dscaledcn = self%ipea_cn(izp) / (2.0_wp * sqrt(ptr%cn(iat) + cn_reg**2))

         ! IP and EA splitting error function
         erf_ip = erf(self%split_exp * (wfn%qat(iat, 1) - self%split_offset))
         erf_ea = erf(self%split_exp * (wfn%qat(iat, 1) + self%split_offset))
         ! Splitting function derivative
         dsplitdcn = 2.0_wp + self%split_slope * (erf_ip + erf_ea)
         
         do ish = 1, self%nsh_at(iat)
            gradient(:, :) = gradient + 0.5_wp * self%ipea(ish, izp) &
               & * wfn%qsh(ii+ish, 1) * dscaledcn * dsplitdcn * ptr%dcndr(:, :, iat)
            sigma(:, :) = sigma + 0.5_wp * self%ipea(ish, izp) &
               & * wfn%qsh(ii+ish, 1) * dscaledcn * dsplitdcn * ptr%dcndL(:, :, iat)
         end do
      end do
   else
      do iat = 1, mol%nat
         izp = mol%id(iat)
         ! Coordination number scaling derivative with regularization
         dscaledcn = self%ipea_cn(izp) / (2.0_wp * sqrt(ptr%cn(iat) + cn_reg**2))
         ! IP and EA splitting error function
         erf_ip = erf(self%split_exp * (wfn%qat(iat, 1) - self%split_offset))
         erf_ea = erf(self%split_exp * (wfn%qat(iat, 1) + self%split_offset))
         ! Splitting function derivative
         dsplitdcn = 2.0_wp + self%split_slope * (erf_ip + erf_ea)
         
         gradient(:, :) = gradient + 0.5_wp * self%ipea(1, izp) &
            & * wfn%qat(iat, 1) * dscaledcn * dsplitdcn * ptr%dcndr(:, :, iat)
         sigma(:, :) = sigma + 0.5_wp * self%ipea(1, izp) &
            & * wfn%qat(iat, 1) * dscaledcn * dsplitdcn * ptr%dcndL(:, :, iat)
      end do
   end if

end subroutine get_gradient


!> Get information about density dependent quantities used in the energy
pure function variable_info(self) result(info)
   use tblite_scf_info, only : scf_info, atom_resolved, shell_resolved
   !> Instance of the first-order tight-binding container
   class(onsite_firstorder), intent(in) :: self
   !> Information on the required potential data
   type(scf_info) :: info

   info = scf_info(charge=merge(shell_resolved, atom_resolved, self%shell_resolved))
end function variable_info


!> Inspect container cache and reallocate it in case of type mismatch
subroutine taint(cache, ptr)
   !> Instance of the container cache
   type(container_cache), target, intent(inout) :: cache
   !> Reference to the coulomb cache
   type(coulomb_cache), pointer, intent(out) :: ptr

   if (allocated(cache%raw)) then
      call view(cache, ptr)
      if (associated(ptr)) return
      deallocate(cache%raw)
   end if

   if (.not.allocated(cache%raw)) then
      block
         type(coulomb_cache), allocatable :: tmp
         allocate(tmp)
         call move_alloc(tmp, cache%raw)
      end block
   end if

   call view(cache, ptr)
end subroutine taint

!> Return reference to container cache after resolving its type
subroutine view(cache, ptr)
   !> Instance of the container cache
   type(container_cache), target, intent(inout) :: cache
   !> Reference to the container cache
   type(coulomb_cache), pointer, intent(out) :: ptr
   nullify(ptr)
   select type(target => cache%raw)
   type is(coulomb_cache)
      ptr => target
   end select
end subroutine view

end module tblite_coulomb_firstorder
