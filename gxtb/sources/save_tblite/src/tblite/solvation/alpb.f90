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

!> @file tblite/solvation/alpb.f90
!> Provides the analytical linearized Poission-Boltzmann model.

!> Analytical linearized Poisson-Boltzmann implicit solvation model.
!>
!> Implements a reaction field model of the Generalized Born type.
module tblite_solvation_alpb
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use mctc_io_math, only : matdet_3x3
   use tblite_blas, only : dot, gemv, symv
   use tblite_container_cache, only : container_cache
   use tblite_mesh_lebedev, only : grid_size, get_angular_grid, list_bisection
   use tblite_scf_info, only : scf_info, atom_resolved
   use tblite_scf_potential, only : potential_type
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_solvation_born, only : born_integrator, new_born_integrator
   use tblite_solvation_data, only : get_vdw_rad_cosmo
   use tblite_solvation_type, only : solvation_type
   use tblite_solvation_cm5, only : get_cm5_charges
   implicit none
   private

   public :: new_alpb
   public :: born_kernel


   !> Possible Born interaction kernels
   type :: enum_born_kernel
      !> Classical Still kernel
      integer :: still = 1
      !> P16 kernel by Lange (JCTC 2012, 8, 1999-2011)
      integer :: p16 = 2
   end type enum_born_kernel

   !> Actual enumerator for the generalized Born kernels
   type(enum_born_kernel), parameter :: born_kernel = enum_born_kernel()


   !> Input for ALPB solvation
   type, public :: alpb_input
      !> Dielectric constant
      real(wp) :: dielectric_const
      !> Scaling factor for Born radii
      real(wp) :: born_scale = 1.0_wp
      !> Offset parameter for Born radii integration
      real(wp) :: born_offset = 0.0_wp
      !> Van-der-Waals radii for all species
      real(wp), allocatable :: rvdw(:)
      !> Dielectric descreening parameter
      real(wp), allocatable :: descreening(:)
      !> Interaction kernel
      integer :: kernel = born_kernel%p16
      !> Use analytical linearized Poisson-Boltzmann model
      logical :: alpb = .true.
      !> Solvent for parameter selection
      character(len=:), allocatable :: solvent
   end type alpb_input

   !> Provide constructor for ALPB input
   interface alpb_input
      module procedure :: create_alpb_input
   end interface alpb_input

   !> Definition of ALPB/GBSA model
   type, public, extends(solvation_type) :: alpb_solvation
      !> Dielectric function
      real(wp) :: keps
      !> Analytical linearized Poisson-Boltzmann constant
      real(wp) :: alpbet
      !> Integrator for Born radii
      type(born_integrator) :: gbobc
      !> Interaction kernel
      integer :: kernel
      !> Use CM5 charges (GFN1-xTB compatibility)
      logical :: useCM5 = .false.
   contains
      !> Update cache from container
      procedure :: update
      !> Return dependency on density
      procedure :: variable_info
      !> Get solvation energy
      procedure :: get_energy
      !> Get solvation potential
      procedure :: get_potential
      !> Get solvation gradient
      procedure :: get_gradient
   end type alpb_solvation

   !> Provide constructor for ALPB solvation
   interface alpb_solvation
      module procedure :: create_alpb
   end interface alpb_solvation


   !> Restart data for ALPB calculation
   type :: alpb_cache
      !> Screening matrix
      real(wp), allocatable :: jmat(:, :)
      !> Scratch array for screening potential intermediates
      real(wp), allocatable :: vat(:)
      !> Born radii
      real(wp), allocatable :: rad(:)
      !> Derivatives of Born radii w.r.t. cartesian displacements
      real(wp), allocatable :: draddr(:, :, :)
      !> Scratch workspace for gradient construction
      real(wp), allocatable :: scratch(:)
      !> Workspace for atomic charges
      real(wp), allocatable :: qscratch(:)
      !> CM5 charges (only required for GFN1 compatibility)
      real(wp), allocatable :: cm5(:)
      !> CM5 charge derivatives
      real(wp), allocatable :: dcm5dr(:,:,:)
   end type alpb_cache


   !> Identifier for container
   character(len=*), parameter :: label = "alpb/gbsa reaction field model"

   real(wp), parameter :: zetaP16 = 1.028_wp
   real(wp), parameter :: zetaP16o16 = zetaP16 / 16.0_wp
   real(wp), parameter :: alpha_alpb = 0.571412_wp


contains


!> Consturctor for ALPB input to properly assign allocatable strings
function create_alpb_input(dielectric_const, solvent, alpb, kernel) result(self)
   !> Dielectric constant
   real(wp), intent(in) :: dielectric_const
   !> Solvent for parameter selection
   character(len=*), intent(in), optional :: solvent
   !> Use analytical linearized Poisson-Boltzmann model
   logical, intent(in), optional :: alpb
   !> Interaction kernel
   integer, intent(in), optional :: kernel

   type(alpb_input) :: self

   self%dielectric_const = dielectric_const

   if (present(solvent)) then 
      self%solvent = solvent
   end if

   if (present(alpb)) then 
      self%alpb = alpb
   end if

   if (present(kernel)) then 
      self%kernel = kernel
   end if

end function create_alpb_input


!> Create new ALPB solvation model
subroutine new_alpb(self, mol, input, method)
   !> Instance of the solvation model
   type(alpb_solvation), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Input for ALPB solvation
   type(alpb_input), intent(in) :: input
   !> Method for parameter selection
   character(len=*), intent(in), optional :: method

   real(wp), allocatable :: rvdw(:)

   self%label = label
   self%alpbet = merge(alpha_alpb / input%dielectric_const, 0.0_wp, input%alpb)
   self%keps = (1.0_wp/input%dielectric_const - 1.0_wp) / (1.0_wp + self%alpbet)
   self%kernel = input%kernel
   if (allocated(input%solvent) .and. present(method)) then
      self%useCM5 = trim(method) == 'gfn1'
   endif

   if (allocated(input%rvdw)) then
      rvdw = input%rvdw
   else
      rvdw = get_vdw_rad_cosmo(mol%num)
   end if

   call new_born_integrator(self%gbobc, mol, rvdw, descreening=input%descreening, &
      & born_scale=input%born_scale, born_offset=input%born_offset)
end subroutine new_alpb


!> Type constructor for ALPB solvation
function create_alpb(mol, input, method) result(self)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Input for ALPB solvation
   type(alpb_input), intent(in) :: input
   !> Method for parameter selection
   character(len=*), intent(in), optional :: method
   !> Instance of the solvation model
   type(alpb_solvation) :: self

   call new_alpb(self, mol, input, method)
end function create_alpb


!> Update cache from container
subroutine update(self, mol, cache)
   !> Instance of the solvation model
   class(alpb_solvation), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache

   type(alpb_cache), pointer :: ptr
   real(wp) :: adet

   call taint(cache, ptr)

   if (.not.allocated(ptr%jmat)) then
      allocate(ptr%jmat(mol%nat, mol%nat))
   end if
   if (.not.allocated(ptr%vat)) then
      allocate(ptr%vat(mol%nat))
   end if
   if (.not.allocated(ptr%rad)) then
      allocate(ptr%rad(mol%nat))
   end if
   if (.not.allocated(ptr%draddr)) then
      allocate(ptr%draddr(3, mol%nat, mol%nat))
   end if
   if (.not.allocated(ptr%qscratch))then
      allocate(ptr%qscratch(mol%nat))
   endif 
   if (self%useCM5)then
      if (.not.allocated(ptr%cm5))then
         allocate(ptr%cm5(mol%nat)) 
      endif
      if (.not.allocated(ptr%dcm5dr))then
         allocate(ptr%dcm5dr(3, mol%nat, mol%nat))
      endif
      call get_cm5_charges(mol, ptr%cm5, ptr%dcm5dr)
   endif
   if (self%useCM5.and..not.allocated(ptr%scratch))then
      allocate(ptr%scratch(mol%nat))
   endif

   call self%gbobc%get_rad(mol, ptr%rad, ptr%draddr)
   ptr%jmat(:, :) = 0.0_wp
   select case(self%kernel)
   case(born_kernel%p16)
      call add_born_mat_p16(mol%nat, mol%xyz, self%keps, ptr%rad, ptr%jmat)
   case(born_kernel%still)
      call add_born_mat_still(mol%nat, mol%xyz, self%keps, ptr%rad, ptr%jmat)
   end select

   if (self%alpbet > 0.0_wp) then
      call get_adet(mol%nat, mol%xyz, self%gbobc%vdwr, adet)

      ptr%jmat(:mol%nat, :mol%nat) = ptr%jmat(:mol%nat, :mol%nat) &
         & + self%keps * self%alpbet / adet
   end if
end subroutine update


!> Get solvation energy
subroutine get_energy(self, mol, cache, wfn, energies)
   !> Instance of the solvation model
   class(alpb_solvation), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Solvation free energy
   real(wp), intent(inout) :: energies(:)

   type(alpb_cache), pointer :: ptr

   call view(cache, ptr)
   
   if(self%useCM5)then
      ptr%qscratch(:) = wfn%qat(:, 1) + ptr%cm5(:)
   else
      ptr%qscratch(:) = wfn%qat(:, 1)
   endif

   call symv(ptr%jmat, ptr%qscratch(:), ptr%vat, alpha=0.5_wp)
   energies(:) = energies + ptr%vat * ptr%qscratch(:)
end subroutine get_energy


!> Get solvation potential
subroutine get_potential(self, mol, cache, wfn, pot)
   !> Instance of the solvation model
   class(alpb_solvation), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Density dependent potential
   type(potential_type), intent(inout) :: pot

   type(alpb_cache), pointer :: ptr

   call view(cache, ptr)

   if(self%useCM5)then
      ptr%qscratch(:) = wfn%qat(:, 1) + ptr%cm5(:)
   else
      ptr%qscratch(:) = wfn%qat(:, 1)
   endif

   call symv(ptr%jmat, ptr%qscratch(:), pot%vat(:, 1), beta=1.0_wp)
end subroutine get_potential


!> Get solvation gradient
subroutine get_gradient(self, mol, cache, wfn, gradient, sigma)
   !> Instance of the solvation model
   class(alpb_solvation), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Molecular gradient of the solvation free energy
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Strain derivatives of the solvation free energy
   real(wp), contiguous, intent(inout) :: sigma(:, :)

   type(alpb_cache), pointer :: ptr
   real(wp) :: energy

   energy = 0.0_wp
   call view(cache, ptr)

   if(self%useCM5)then
      ptr%qscratch(:) = wfn%qat(:, 1) + ptr%cm5(:)
   else
      ptr%qscratch(:) = wfn%qat(:, 1)
   endif

   select case(self%kernel)
   case(born_kernel%p16)
      call add_born_deriv_p16(mol%nat, mol%xyz, &
         & ptr%qscratch(:), self%keps, ptr%rad, ptr%draddr, energy, gradient)
   case(born_kernel%still)
      call add_born_deriv_still(mol%nat, mol%xyz, &
         & ptr%qscratch(:), self%keps, ptr%rad, ptr%draddr, energy, gradient)
   end select

   if (self%alpbet > 0.0_wp) then
      call get_adet_deriv(mol%nat, mol%xyz, self%gbobc%vdwr, self%kEps*self%alpbet, &
         & ptr%qscratch(:), gradient)
   end if

   if(self%useCM5)then
      call gemv(ptr%jmat, ptr%qscratch, ptr%scratch)
      call gemv(ptr%dcm5dr, ptr%scratch, gradient, beta=1.0_wp)
   endif

end subroutine get_gradient


!> Return dependency on density
pure function variable_info(self) result(info)
   !> Instance of the solvation model
   class(alpb_solvation), intent(in) :: self
   !> Information on the required potential data
   type(scf_info) :: info

   info = scf_info(charge=atom_resolved)
end function variable_info


subroutine taint(cache, ptr)
   type(container_cache), target, intent(inout) :: cache
   type(alpb_cache), pointer, intent(out) :: ptr

   if (allocated(cache%raw)) then
      call view(cache, ptr)
      if (associated(ptr)) return
      deallocate(cache%raw)
   end if

   if (.not.allocated(cache%raw)) then
      block
         type(alpb_cache), allocatable :: tmp
         allocate(tmp)
         call move_alloc(tmp, cache%raw)
      end block
   end if

   call view(cache, ptr)
end subroutine taint

subroutine view(cache, ptr)
   type(container_cache), target, intent(inout) :: cache
   type(alpb_cache), pointer, intent(out) :: ptr
   nullify(ptr)
   select type(target => cache%raw)
   type is(alpb_cache)
      ptr => target
   end select
end subroutine view


subroutine add_born_mat_p16(nat, xyz, keps, brad, Amat)
   !> Number of atoms
   integer, intent(in) :: nat
   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)
   !> Dielectric screening
   real(wp), intent(in) :: keps
   !> Born radii
   real(wp), intent(in) :: brad(:)
   !> Interaction matrix
   real(wp), intent(inout) :: Amat(:, :)

   integer :: iat, jat
   real(wp) :: r1, ab, arg, fgb, dfgb, bp, vec(3)

   ! omp parallel do default(none) shared(Amat, ntpair, ppind, ddpair, brad, keps) &
   ! omp private(kk, iat, jat, r1, ab, arg, fgb, dfgb)
   do iat = 1, nat
      do jat = 1, iat - 1
         vec(:) = xyz(:, iat) - xyz(:, jat)
         r1 = norm2(vec)

         ab = sqrt(brad(iat) * brad(jat))
         arg = ab / (ab + zetaP16o16*r1) ! ab / (1 + ζR/(16·ab))
         arg = arg * arg ! ab / (1 + ζR/(16·ab))²
         arg = arg * arg ! ab / (1 + ζR/(16·ab))⁴
         arg = arg * arg ! ab / (1 + ζR/(16·ab))⁸
         arg = arg * arg ! ab / (1 + ζR/(16·ab))¹⁶
         fgb = r1 + ab*arg
         dfgb = 1.0_wp / fgb

         Amat(iat, jat) = keps*dfgb + Amat(iat, jat)
         Amat(jat, iat) = keps*dfgb + Amat(jat, iat)
      enddo
      ! self-energy part
      bp = 1.0_wp/brad(iat)
      Amat(iat, iat) = Amat(iat, iat) + keps*bp
   enddo

end subroutine add_born_mat_p16


subroutine add_born_deriv_p16(nat, xyz, qat, keps, &
      & brad, brdr, energy, gradient)
   !> Number of atoms
   integer, intent(in) :: nat
   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)
   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)
   !> Dielectric screening
   real(wp), intent(in) :: keps
   !> Born radii
   real(wp), intent(in) :: brad(:)
   !> Derivative of Born radii w.r.t. cartesian coordinates
   real(wp), contiguous, intent(in) :: brdr(:, :, :)
   !> Total Born solvation energy
   real(wp), intent(out) :: energy
   !> Deriatives of Born solvation energy
   real(wp), contiguous, intent(inout) :: gradient(:, :)

   integer :: iat, jat
   real(wp) :: vec(3), r2, r1, ab, arg1, arg16, qq, fgb, dfgb, dfgb2, egb
   real(wp) :: dEdbri, dEdbrj, dG(3), ap, bp, dS(3, 3)
   real(wp), allocatable :: dEdbr(:)

   allocate(dEdbr(nat), source = 0.0_wp )

   egb = 0._wp
   dEdbr(:) = 0._wp

   ! GB energy and gradient
   ! omp parallel do default(none) reduction(+:egb, gradient, dEdbr) &
   ! omp private(iat, jat, vec, r1, r2, ab, arg1, arg16, fgb, dfgb, dfgb2, ap, &
   ! omp& bp, qq, dEdbri, dEdbrj, dG, dS) &
   ! omp shared(keps, qat, ntpair, ddpair, ppind, brad)
   do iat = 1, nat
      do jat = 1, iat - 1
         vec(:) = xyz(:, iat) - xyz(:, jat)
         r1 = norm2(vec)
         r2 = r1*r1

         qq = qat(iat)*qat(jat)

         ab = sqrt(brad(iat) * brad(jat))
         arg1 = ab / (ab + zetaP16o16*r1) ! 1 / (1 + ζR/(16·ab))
         arg16 = arg1 * arg1 ! 1 / (1 + ζR/(16·ab))²
         arg16 = arg16 * arg16 ! 1 / (1 + ζR/(16·ab))⁴
         arg16 = arg16 * arg16 ! 1 / (1 + ζR/(16·ab))⁸
         arg16 = arg16 * arg16 ! 1 / (1 + ζR/(16·ab))¹⁶

         fgb = r1 + ab*arg16
         dfgb = 1.0_wp / fgb
         dfgb2 = dfgb * dfgb

         egb = egb + qq*keps*dfgb

         ! (1 - ζ/(1 + Rζ/(16 ab))^17)/(R + ab/(1 + Rζ/(16 ab))¹⁶)²
         ap = (1.0_wp - zetaP16 * arg1 * arg16) * dfgb2
         dG(:) = ap * vec * keps / r1 * qq
         gradient(:, iat) = gradient(:, iat) - dG
         gradient(:, jat) = gradient(:, jat) + dG

         ! -(Rζ/(2·ab²·(1 + Rζ/(16·ab))¹⁷) + 1/(2·ab·(1 + Rζ/(16·ab))¹⁶))/(R + ab/(1 + Rζ/(16·ab))¹⁶)²
         bp = -0.5_wp*(r1 * zetaP16 / ab * arg1 + 1.0_wp) / ab * arg16 * dfgb2
         dEdbri = brad(jat) * bp * keps * qq
         dEdbrj = brad(iat) * bp * keps * qq
         dEdbr(iat) = dEdbr(iat) + dEdbri
         dEdbr(jat) = dEdbr(jat) + dEdbrj

      end do

      ! self-energy part
      bp = 1._wp/brad(iat)
      qq = qat(iat)*bp
      egb = egb + 0.5_wp*qat(iat)*qq*keps
      dEdbri = -0.5_wp*keps*qq*bp
      dEdbr(iat) = dEdbr(iat) + dEdbri*qat(iat)
      !gradient = gradient + brdr(:, :, i) * dEdbri*qat(i)
   enddo

   ! contract with the Born radii derivatives
   call gemv(brdr, dEdbr, gradient, beta=1.0_wp)

   energy = egb

end subroutine add_born_deriv_p16


pure subroutine add_born_mat_still(nat, xyz, keps, brad, Amat)
   !> Number of atoms
   integer, intent(in) :: nat
   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)
   !> Dielectric screening
   real(wp), intent(in) :: keps
   !> Born radii
   real(wp), intent(in) :: brad(:)
   !> Interaction matrix
   real(wp), intent(inout) :: Amat(:, :)

   integer  :: i, j
   real(wp), parameter :: a13=1.0_wp/3.0_wp
   real(wp), parameter :: a4=0.25_wp
   real(wp), parameter :: sqrt2pi = sqrt(2.0_wp/pi)
   real(wp) :: aa, vec(3), r1, r2, bp
   real(wp) :: dd, expd, fgb2, dfgb

   do i = 1, nat
      do j = 1, i - 1
         vec(:) = xyz(:, i) - xyz(:, j)
         r1 = norm2(vec)
         r2 = r1*r1

         aa = brad(i)*brad(j)
         dd = a4*r2/aa
         expd = exp(-dd)
         fgb2 = r2+aa*expd
         dfgb = 1.0_wp/sqrt(fgb2)
         Amat(i, j) = keps*dfgb + Amat(i, j)
         Amat(j, i) = keps*dfgb + Amat(j, i)
      enddo

      ! self-energy part
      bp = 1._wp/brad(i)
      Amat(i, i) = Amat(i, i) + keps*bp
   enddo

end subroutine add_born_mat_still


subroutine add_born_deriv_still(nat, xyz, qat, keps, &
      & brad, brdr, energy, gradient)
   !> Number of atoms
   integer, intent(in) :: nat
   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)
   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)
   !> Dielectric screening
   real(wp), intent(in) :: keps
   !> Born radii
   real(wp), intent(in) :: brad(:)
   !> Derivative of Born radii w.r.t. cartesian coordinates
   real(wp), contiguous, intent(in) :: brdr(:, :, :)
   !> Total Born solvation energy
   real(wp), intent(out) :: energy
   !> Deriatives of Born solvation energy
   real(wp), contiguous, intent(inout) :: gradient(:, :)

   integer :: i, j
   real(wp), parameter :: a13=1._wp/3._wp
   real(wp), parameter :: a4=0.25_wp
   real(wp) :: aa, r2, fgb2
   real(wp) :: qq, dd, expd, dfgb, dfgb2, dfgb3, egb, ap, bp
   real(wp) :: grddbi, grddbj
   real(wp) :: dr(3), r1, vec(3)
   real(wp), allocatable :: grddb(:)

   allocate(grddb(nat), source = 0.0_wp )

   egb = 0._wp
   grddb(:) = 0._wp

   ! GB energy and gradient

   ! compute energy and fgb direct and radii derivatives
   do i = 1, nat
      do j = 1, i - 1
         vec(:) = xyz(:, i) - xyz(:, j)
         r1 = norm2(vec)
         r2 = r1*r1

         ! dielectric scaling of the charges
         qq = qat(i)*qat(j)
         aa = brad(i)*brad(j)
         dd = a4*r2/aa
         expd = exp(-dd)
         fgb2 = r2+aa*expd
         dfgb2 = 1._wp/fgb2
         dfgb = sqrt(dfgb2)
         dfgb3 = dfgb2*dfgb*keps

         egb = egb + qq*keps*dfgb

         ap = (1._wp-a4*expd)*dfgb3
         dr = ap*vec
         gradient(:, i) = gradient(:, i) - dr*qq
         gradient(:, j) = gradient(:, j) + dr*qq

         bp = -0.5_wp*expd*(1._wp+dd)*dfgb3
         grddbi = brad(j)*bp
         grddbj = brad(i)*bp
         grddb(i) = grddb(i) + grddbi*qq
         grddb(j) = grddb(j) + grddbj*qq

      enddo

      ! self-energy part
      bp = 1._wp/brad(i)
      qq = qat(i)*bp
      egb = egb + 0.5_wp*qat(i)*qq*keps
      grddbi = -0.5_wp*keps*qq*bp
      grddb(i) = grddb(i) + grddbi*qat(i)
   enddo

   ! contract with the Born radii derivatives
   call gemv(brdr, grddb, gradient, beta=1.0_wp)

   energy = egb

end subroutine add_born_deriv_still


subroutine get_adet(nat, xyz, rad, aDet)
   !> Number of atoms
   integer, intent(in) :: nat
   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)
   !> Atomic radii
   real(wp), intent(in) :: rad(:)
   !> Shape descriptor of the structure
   real(wp), intent(out) :: aDet

   integer :: iat
   real(wp) :: r2, rad2, rad3, vol, vec(3), center(3), inertia(3, 3)
   real(wp), parameter :: tof = 2.0_wp/5.0_wp, unity(3, 3) = reshape(&
      & [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], &
      & [3, 3])

   vol = 0.0_wp
   center(:) = 0.0_wp
   do iat = 1, nat
      rad2 = rad(iat) * rad(iat)
      rad3 = rad2 * rad(iat)
      vol = vol + rad3
      center(:) = center + xyz(:, iat) * rad3
   end do
   center = center / vol

   inertia(:, :) = 0.0_wp
   do iat = 1, nat
      rad2 = rad(iat) * rad(iat)
      rad3 = rad2 * rad(iat)
      vec(:) = xyz(:, iat) - center
      r2 = sum(vec**2)
      inertia(:, :) = inertia + rad3 * ((r2 + tof*rad2) * unity &
         & - spread(vec, 1, 3) * spread(vec, 2, 3))
   end do

   aDet = sqrt(matdet_3x3(inertia)**(1.0_wp/3.0_wp)/(tof*vol))

end subroutine get_adet


subroutine get_adet_deriv(nAtom, xyz, rad, kEps, qvec, gradient)
   !> Number of atoms
   integer, intent(in) :: nAtom
   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)
   !> Atomic radii
   real(wp), intent(in) :: rad(:)
   real(wp), intent(in) :: kEps
   real(wp), intent(in) :: qvec(:)
   !> Molecular gradient
   real(wp), intent(inout) :: gradient(:, :)

   integer :: iat
   real(wp) :: r2, rad2, rad3, vol, vec(3), center(3), inertia(3, 3), aDet
   real(wp) :: aDeriv(3, 3), qtotal
   real(wp), parameter :: tof = 2.0_wp/5.0_wp, unity(3, 3) = reshape(&
      & [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], &
      & [3, 3])

   qtotal = 0.0_wp
   vol = 0.0_wp
   center(:) = 0.0_wp
   do iat = 1, nAtom
      rad2 = rad(iat) * rad(iat)
      rad3 = rad2 * rad(iat)
      vol = vol + rad3
      center(:) = center + xyz(:, iat) * rad3
      qtotal = qtotal + qvec(iat)
   end do
   center = center / vol

   inertia(:, :) = 0.0_wp
   do iat = 1, nAtom
      rad2 = rad(iat) * rad(iat)
      rad3 = rad2 * rad(iat)
      vec(:) = xyz(:, iat) - center
      r2 = sum(vec**2)
      inertia(:, :) = inertia + rad3 * ((r2 + tof*rad2) * unity &
         & - spread(vec, 1, 3) * spread(vec, 2, 3))
   end do
   aDet = sqrt(matdet_3x3(inertia)**(1.0_wp/3.0_wp)/(tof*vol))

   aDeriv(:, :) = reshape([&
      & inertia(1,1)*(inertia(2,2)+inertia(3,3))-inertia(1,2)**2-inertia(1,3)**2, &
      & inertia(1,2)*inertia(3,3)-inertia(1,3)*inertia(2,3), & ! xy
      & inertia(1,3)*inertia(2,2)-inertia(1,2)*inertia(3,2), & ! xz
      & inertia(1,2)*inertia(3,3)-inertia(1,3)*inertia(2,3), & ! xy
      & inertia(2,2)*(inertia(1,1)+inertia(3,3))-inertia(1,2)**2-inertia(2,3)**2, &
      & inertia(1,1)*inertia(2,3)-inertia(1,2)*inertia(1,3), & ! yz
      & inertia(1,3)*inertia(2,2)-inertia(1,2)*inertia(3,2), & ! xz
      & inertia(1,1)*inertia(2,3)-inertia(1,2)*inertia(1,3), & ! yz
      & inertia(3,3)*(inertia(1,1)+inertia(2,2))-inertia(1,3)**2-inertia(2,3)**2],&
      & shape=[3, 3]) * (250.0_wp / (48.0_wp * vol**3 * aDet**5)) &
      & * (-0.5_wp * kEps * qtotal**2 / aDet**2)

   do iat = 1, nAtom
      rad2 = rad(iat) * rad(iat)
      rad3 = rad2 * rad(iat)
      vec(:) = xyz(:, iat) - center
      gradient(:, iat) = gradient(:, iat) + rad3 * matmul(aderiv, vec)
   end do

end subroutine get_adet_deriv


end module tblite_solvation_alpb
