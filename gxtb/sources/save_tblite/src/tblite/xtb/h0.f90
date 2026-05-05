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

!> @file tblite/xtb/h0.f90
!> Provides the effective core Hamiltonian for xTB.

!> Implementation of the effective core Hamiltonian used in the extended tight binding.
module tblite_xtb_h0
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use tblite_adjlist, only : adjacency_list
   use tblite_basis_cache, only : basis_cache
   use tblite_basis_type, only : basis_type
   use tblite_integral_diat_trafo, only: diat_trafo, diat_trafo_grad
   use tblite_integral_multipole, only : multipole_cgto, multipole_grad_cgto, &
      & maxl, msao, smap
   use tblite_integral_overlap, only : overlap_cgto, overlap_grad_cgto
   use tblite_scf_potential, only : potential_type
   use tblite_xtb_spec, only : tb_h0spec
   implicit none
   private

   public :: new_hamiltonian
   public :: get_selfenergy, get_anisotropy, get_anisotropy_gradient
   public :: get_hamiltonian, get_occupation, get_hamiltonian_gradient, get_number_electrons


   type, public :: tb_hamiltonian
      !> Atomic level information
      real(wp), allocatable :: selfenergy(:, :)
      !> Coordination number dependence of the atomic levels
      real(wp), allocatable :: kcn(:, :)
      !> Electronegativity scaled coordination number dependence of the atomic levels
      real(wp), allocatable :: kcn_en(:, :)
      !> Charge dependence of the atomic levels
      real(wp), allocatable :: kq1(:, :)
      !> Charge dependence of the atomic levels
      real(wp), allocatable :: kq2(:, :)
      !> Enhancement factor to scale the Hamiltonian elements
      real(wp), allocatable :: hscale(:, :, :, :)
      !> Polynomial coefficients for square-root distance dependent enhancement factor
      real(wp), allocatable :: shpoly(:, :)
      !> Polynomial coefficients for linear distance dependent enhancement factor
      real(wp), allocatable :: shpoly2(:, :)
      !> Polynomial coefficients for square distance dependent enhancement factor
      real(wp), allocatable :: shpoly4(:, :)
      !> Atomic radius for polynomial enhancement
      real(wp), allocatable :: rad(:)
      !> Coefficients for dipole Hamiltonian correction
      real(wp), allocatable :: dip_scale(:, :)
      !> Van der Waals radii for anisotropic counting function
      real(wp), allocatable :: rvdw(:, :)
      !> Exponent for the anisotropic counting function
      real(wp) :: aniso_exp
      !> Reference occupation numbers
      real(wp), allocatable :: refocc(:, :)
      !> Diatomic frame scaling of sigma bonding contribution
      real(wp), allocatable :: ksig(:, :)
      !> Diatomic frame scaling of pi bonding contribution
      real(wp), allocatable :: kpi(:, :)
      !> Diatomic frame scaling of delta bonding contribution
      real(wp), allocatable :: kdel(:, :)
      !> Perform diatomic frame scaling for the overlap in the Hamiltonian
      logical :: do_diat_scale = .false.
   end type tb_hamiltonian

   real(wp), parameter :: tosqrtpi = 2.0_wp / sqrt(pi)

contains


!> Constructor for a new Hamiltonian object, consumes a Hamiltonian specification
subroutine new_hamiltonian(self, mol, bas, spec)
   !> Hamiltonian object
   type(tb_hamiltonian), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   class(basis_type), intent(in) :: bas
   !> Hamiltonian specification
   class(tb_h0spec), intent(in) :: spec

   integer :: mshell

   mshell = maxval(bas%nsh_id)
   allocate(self%selfenergy(mshell, mol%nid), &
      & self%kcn(mshell, mol%nid), self%kcn_en(mshell, mol%nid), & 
      & self%kq1(mshell, mol%nid), self%kq2(mshell, mol%nid))
   call spec%get_selfenergy(mol, bas, self%selfenergy)
   call spec%get_cnshift(mol, bas, self%kcn)
   call spec%get_cnenshift(mol, bas, self%kcn_en)
   call spec%get_q1shift(mol, bas, self%kq1)
   call spec%get_q2shift(mol, bas, self%kq2)

   allocate(self%hscale(mshell, mshell, mol%nid, mol%nid))
   call spec%get_hscale(mol, bas, self%hscale)

   allocate(self%rad(mol%nid), self%shpoly(mshell, mol%nid), &
      & self%shpoly2(mshell, mol%nid), self%shpoly4(mshell, mol%nid))
   call spec%get_rad(mol, bas, self%rad)
   call spec%get_shpoly(mol, bas, self%shpoly)
   call spec%get_shpoly2(mol, bas, self%shpoly2)
   call spec%get_shpoly4(mol, bas, self%shpoly4)

   allocate(self%dip_scale(mol%nid, mol%nid), self%rvdw(mol%nid, mol%nid))
   call spec%get_anisotropy(mol, bas, self%dip_scale, self%aniso_exp)
   call spec%get_rvdw(mol, bas, self%rvdw)

   allocate(self%refocc(mshell, mol%nid))
   call spec%get_reference_occ(mol, bas, self%refocc)

   allocate(self%ksig(mol%nid, mol%nid), self%kpi(mol%nid, mol%nid), &
      & self%kdel(mol%nid, mol%nid))
   call spec%get_diat_scale(mol, bas, self%ksig, self%kpi, self%kdel, &
      & self%do_diat_scale)

end subroutine new_hamiltonian


subroutine get_selfenergy(h0, id, ish_at, nshell, cn, cn_en, qat, selfenergy, &
   & dsedcn, dsedcn_en, dsedq)
   !> Hamiltonian object
   type(tb_hamiltonian), intent(in) :: h0
   !> Species identifiers
   integer, intent(in) :: id(:)
   !> Index of first shell per atom
   integer, intent(in) :: ish_at(:)
   !> Number of shells per element
   integer, intent(in) :: nshell(:)
   !> Optional coordination number
   real(wp), intent(in), optional :: cn(:)
   !> Optional electronegativity scaled coordination number
   real(wp), intent(in), optional :: cn_en(:)
   !> Optional atomic partial charge
   real(wp), intent(in), optional :: qat(:)
   !> Reference energy levels 
   real(wp), intent(out) :: selfenergy(:)
   !> Derivative of energy levels w.r.t. coordination number
   real(wp), intent(out), optional :: dsedcn(:)
   !> Derivative of energy levels w.r.t. EN-scaled coordination number
   real(wp), intent(out), optional :: dsedcn_en(:)
   !> Derivative of energy levels w.r.t. atomic charge
   real(wp), intent(out), optional :: dsedq(:)

   integer :: iat, izp, ish, ii

   selfenergy(:) = 0.0_wp
   if (present(dsedcn)) dsedcn(:) = 0.0_wp
   if (present(dsedq)) dsedq(:) = 0.0_wp
   do iat = 1, size(id)
      izp = id(iat)
      ii = ish_at(iat)
      do ish = 1, nshell(izp)
         selfenergy(ii+ish) = h0%selfenergy(ish, izp)
      end do
   end do
   if (present(cn)) then
      if (present(dsedcn)) then
         do iat = 1, size(id)
            izp = id(iat)
            ii = ish_at(iat)
            do ish = 1, nshell(izp)
               selfenergy(ii+ish) = selfenergy(ii+ish) &
                  & - h0%kcn(ish, izp) * cn(iat)
               dsedcn(ii+ish) = -h0%kcn(ish, izp)
            end do
         end do
      else
         do iat = 1, size(id)
            izp = id(iat)
            ii = ish_at(iat)
            do ish = 1, nshell(izp)
               selfenergy(ii+ish) = selfenergy(ii+ish) &
                  & - h0%kcn(ish, izp) * cn(iat)
            end do
         end do
      end if
   end if
   if (present(cn_en)) then
      if (present(dsedcn_en)) then
         do iat = 1, size(id)
            izp = id(iat)
            ii = ish_at(iat)
            do ish = 1, nshell(izp)
               selfenergy(ii+ish) = selfenergy(ii+ish) &
                  & - h0%kcn_en(ish, izp) * cn_en(iat)
               dsedcn_en(ii+ish) = -h0%kcn_en(ish, izp)
            end do
         end do
      else
         do iat = 1, size(id)
            izp = id(iat)
            ii = ish_at(iat)
            do ish = 1, nshell(izp)
               selfenergy(ii+ish) = selfenergy(ii+ish) &
                  & - h0%kcn_en(ish, izp) * cn_en(iat)
            end do
         end do
      end if
   end if
   if (present(qat)) then
      if (present(dsedq)) then
         do iat = 1, size(id)
            izp = id(iat)
            ii = ish_at(iat)
            do ish = 1, nshell(izp)
               selfenergy(ii+ish) = selfenergy(ii+ish) &
                  & - h0%kq1(ish, izp)*qat(iat) - h0%kq2(ish, izp)*qat(iat)**2
               dsedq(ii+ish) = -h0%kq1(ish, izp) - h0%kq2(ish, izp)*2*qat(iat)
            end do
         end do
      else
         do iat = 1, size(id)
            izp = id(iat)
            ii = ish_at(iat)
            do ish = 1, nshell(izp)
               selfenergy(ii+ish) = selfenergy(ii+ish) &
                  & - h0%kq1(ish, izp)*qat(iat) - h0%kq2(ish, izp)*qat(iat)**2
            end do
         end do
      end if
   end if

end subroutine get_selfenergy


subroutine get_anisotropy(h0, mol, aniso_dip)
   !> Hamiltonian object
   type(tb_hamiltonian), intent(in) :: h0
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Anisotropy vectors for dipole Hamiltonian correction
   real(wp), intent(out) :: aniso_dip(:, :)

   integer :: iat, izp, jat, jzp, k
   real(wp) :: vec(3), r2, r, rvdw, count, scale, tmp(3)

   aniso_dip(:, :) = 0.0_wp
   !$omp parallel do schedule(runtime) default(none) &
   !$omp shared(mol, h0, aniso_dip) &
   !$omp private(iat, jat, izp, jzp, k, r2, r, vec, rvdw, count, scale, tmp)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      tmp = 0.0_wp
      do jat = 1, iat-1
         jzp = mol%id(jat)
         vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat)
         r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
         r = sqrt(r2)
         ! Counting function
         rvdw = h0%rvdw(izp, jzp)
         count = 1.0_wp + erf(-h0%aniso_exp * (r - rvdw))
         scale = 0.5_wp * count * h0%dip_scale(izp, jzp)
         ! Anisotropy vector points towards maximum coordination
         tmp(:) = tmp - scale * vec
         do k = 1, 3
            !$omp atomic
            aniso_dip(k, jat) = aniso_dip(k, jat) + scale * vec(k)
         end do
      end do
      do k = 1, 3
         !$omp atomic
         aniso_dip(k, iat) = aniso_dip(k, iat) + tmp(k)
      end do
   end do

end subroutine get_anisotropy


subroutine get_anisotropy_gradient(h0, mol, dEdad, gradient, sigma)
   !> Hamiltonian object
   type(tb_hamiltonian), intent(in) :: h0
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Derivative of the electronic energy w.r.t. the anisotropy vector
   real(wp), intent(in) :: dEdad(:, :)
   !> Derivative of the electronic energy w.r.t. coordinate displacements
   real(wp), intent(inout) :: gradient(:, :)
   !> Derivative of the electronic energy w.r.t. strain deformations
   real(wp), intent(inout) :: sigma(:, :)
   
   integer :: iat, izp, jat, jzp, k, l
   real(wp) :: vec(3), r2, r, rvdw, count, dcount, scale, dscale
   real(wp) :: proj, diff_dEdad(3), dG(3)

   ! Thread-private array for reduction
   ! Set to 0 explicitly as the shared variants are potentially non-zero (inout)
   real(wp), allocatable :: gradient_local(:, :), sigma_local(:, :)

   !$omp parallel default(none) &
   !$omp shared(mol, h0, dEdad, gradient, sigma) &
   !$omp private(iat, jat, izp, jzp, k, l, r2, r, vec, rvdw, count, dcount) &
   !$omp private(scale, dscale, diff_dEdad, proj, dG, gradient_local, sigma_local)
   allocate(gradient_local(size(gradient, 1), size(gradient, 2)), source=0.0_wp)
   allocate(sigma_local(size(sigma, 1), size(sigma, 2)), source=0.0_wp)
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat-1
         jzp = mol%id(jat)
         vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat)
         r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
         r = sqrt(r2)
         ! Counting function and derivative
         rvdw = h0%rvdw(izp, jzp)
         count = 1.0_wp + erf(-h0%aniso_exp * (r - rvdw))
         dcount = - exp(- (h0%aniso_exp * (r - rvdw))**2) &
            & * tosqrtpi * h0%aniso_exp
         scale = 0.5_wp * count * h0%dip_scale(izp, jzp)
         dscale = 0.5_wp * dcount * h0%dip_scale(izp, jzp) / r
         ! Relative energy derivative vector
         diff_dEdad(:) = dEdad(:, jat) - dEdad(:, iat)
         ! Projection onto the interatomic vector
         proj = dot_product(vec, diff_dEdad)
         ! Gradient contribution
         dG(:) = scale * diff_dEdad(:) + dscale * proj * vec(:)
         gradient_local(:, iat) = gradient_local(:, iat) + dG
         gradient_local(:, jat) = gradient_local(:, jat) - dG
      end do
   end do
   !$omp end do
   !$omp critical
   gradient(:, :) = gradient + gradient_local
   sigma(:, :) = sigma + sigma_local
   !$omp end critical
   deallocate(gradient_local, sigma_local)
   !$omp end parallel

end subroutine get_anisotropy_gradient


subroutine get_hamiltonian(mol, trans, list, bas, bcache, h0, selfenergy, &
   & aniso_dip, overlap, dpint, qpint, hamiltonian)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Lattice points within a given realspace cutoff
   real(wp), intent(in) :: trans(:, :)
   !> Neighbour list
   type(adjacency_list), intent(in) :: list
   !> Basis set information
   class(basis_type), intent(in) :: bas
   !> Basis set cache
   type(basis_cache), intent(in) :: bcache
   !> Hamiltonian interaction data
   type(tb_hamiltonian), intent(in) :: h0
   !> Diagonal elememts of the Hamiltonian
   real(wp), intent(in) :: selfenergy(:)
   !> Anisotropy vectors for dipole Hamiltonian correction
   real(wp), intent(in) :: aniso_dip(:, :)
   !> Overlap integral matrix
   real(wp), intent(out) :: overlap(:, :)
   !> Dipole moment integral matrix
   real(wp), intent(out) :: dpint(:, :, :)
   !> Quadrupole moment integral matrix
   real(wp), intent(out) :: qpint(:, :, :)
   !> Effective Hamiltonian
   real(wp), intent(out) :: hamiltonian(:, :)

   integer :: iat, jat, izp, jzp, itr, k, img, inl
   integer :: ish, jsh, is, js, ii, jj, iao, jao, nao, ij, iaosh, jaosh
   real(wp) :: r2, r, rr, rr2, rr4, vec(3), shpolyi, shpolyj, shpoly
   real(wp) :: hij, dtmpj(3), qtmpj(6), hdij, mod_h0_fraction
   real(wp), allocatable :: stmp(:), dtmpi(:, :), qtmpi(:, :), block_overlap(:,:)
   logical :: do_scaled_h0, modify_h0

   overlap(:, :) = 0.0_wp
   dpint(:, :, :) = 0.0_wp
   qpint(:, :, :) = 0.0_wp
   hamiltonian(:, :) = 0.0_wp

   ! Predefine loop-invariant logicals to control H0 scaling
   ! and diatomic frame transformation 
   do_scaled_h0 = allocated(bas%cgto_h0)
   modify_h0 = do_scaled_h0 .or. h0%do_diat_scale

   ! Select if we construct the Hamiltonian with modifications
   mod_h0_fraction = 0.0_wp
   if(modify_h0) then
      mod_h0_fraction = 1.0_wp
   end if

   allocate(stmp(msao(bas%maxl)**2), dtmpi(3, msao(bas%maxl)**2), &
      & qtmpi(6, msao(bas%maxl)**2))
   ! Temporary overlap block between two atoms with maximum size for a single-zeta basis
   allocate(block_overlap(smap(maxval(bas%nsh_id, 1)),smap( maxval(bas%nsh_id, 1))))

   !$omp parallel do schedule(runtime) default(none) &
   !$omp shared(mol, bas, bcache, trans, list, overlap, dpint, qpint, hamiltonian) &
   !$omp shared(h0, selfenergy, aniso_dip) &
   !$omp firstprivate(do_scaled_h0, modify_h0, mod_h0_fraction) &
   !$omp private(iat, jat, izp, jzp, itr, inl, img, is, js, ish, jsh, ii, jj) &
   !$omp private(iao, jao, iaosh, jaosh, nao, ij, k, r2, r, rr, rr2, rr4, vec) &
   !$omp private(shpolyi, shpolyj, shpoly, stmp, dtmpi, qtmpi, dtmpj, qtmpj) &
   !$omp private(hij, hdij, block_overlap)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = bas%ish_at(iat)
      inl = list%inl(iat)
      do img = 1, list%nnl(iat)
         jat = list%nlat(img+inl)
         itr = list%nltr(img+inl)
         jzp = mol%id(jat)
         js = bas%ish_at(jat)
         vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
         r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
         r = sqrt(r2)
         ! Precompute square root, linear, and quadratic distance factors
         rr2 = r / (h0%rad(jzp) + h0%rad(izp))
         rr = sqrt(rr2)
         rr4 = rr2 * r
         ! Get the overlap and dipole integrals for the current diatomic pair
         do ish = 1, bas%nsh_id(izp)
            ii = bas%iao_sh(is+ish)
            iaosh = smap(ish-1)
            do jsh = 1, bas%nsh_id(jzp)
               jj = bas%iao_sh(js+jsh)
               jaosh = smap(jsh-1)

               ! Calculate overlap and multipole integrals
               ! Gauge-origin of the dipole/quadrupole integrals is at atom i
               call multipole_cgto(bas%cgto(jsh, jzp)%raw, bas%cgto(ish, izp)%raw, &
                  & bcache%cgto(jsh, jat), bcache%cgto(ish, iat), &
                  & r2, vec, bas%intcut, stmp, dtmpi, qtmpi)

               ! Calculate shell polynomials
               shpolyi = 1.0_wp + h0%shpoly(ish, izp)*rr &
                  & + h0%shpoly2(ish, izp)*rr2 + h0%shpoly4(ish, izp)*rr4
               shpolyj = 1.0_wp + h0%shpoly(jsh, jzp)*rr &
                  & + h0%shpoly2(jsh, jzp)*rr2 + h0%shpoly4(jsh, jzp)*rr4
               shpoly = shpolyi * shpolyj

               ! Calculate Hamiltonian elements
               hij = 0.5_wp * (1.0_wp - mod_h0_fraction) &
                  & * (selfenergy(is+ish) + selfenergy(js+jsh)) &
                  & * h0%hscale(jsh, ish, jzp, izp) * shpoly

               ! Distribute shell block to Hamiltonian, overlap and multipole matrices
               nao = msao(bas%cgto(jsh, jzp)%raw%ang)
               do iao = 1, msao(bas%cgto(ish, izp)%raw%ang)
                  do jao = 1, nao
                     ij = jao + nao*(iao-1)
                     ! Shift the atom i centered multipole integrals to atom j
                     call shift_operator(vec, stmp(ij), dtmpi(:, ij), qtmpi(:, ij), &
                        & dtmpj, qtmpj)

                     ! Save overlap for possible diatomic frame trafo
                     block_overlap(jaosh+jao, iaosh+iao) = stmp(ij)

                     !$omp atomic
                     overlap(jj+jao, ii+iao) = overlap(jj+jao, ii+iao) &
                        + stmp(ij)

                     hdij = 0.0_wp
                     do k = 1, 3
                        !$omp atomic
                        dpint(k, jj+jao, ii+iao) = dpint(k, jj+jao, ii+iao) &
                           + dtmpi(k, ij)
                        ! Construct dipole Hamiltonian correction
                        hdij = hdij - 0.5_wp * (aniso_dip(k, iat) * dtmpi(k, ij) + &
                           & aniso_dip(k, jat) * dtmpj(k))
                     end do

                     do k = 1, 6
                        !$omp atomic
                        qpint(k, jj+jao, ii+iao) = qpint(k, jj+jao, ii+iao) &
                           + qtmpi(k, ij)
                     end do

                     ! Hamiltonian entry with dipole and quadrupole corrections
                     !$omp atomic
                     hamiltonian(jj+jao, ii+iao) = hamiltonian(jj+jao, ii+iao) &
                        + stmp(ij) * hij + hdij

                     ! Fill lower triangle of the matrices
                     if (iat /= jat) then
                        !$omp atomic
                        overlap(ii+iao, jj+jao) = overlap(ii+iao, jj+jao) &
                           + stmp(ij)

                        do k = 1, 3
                           !$omp atomic
                           dpint(k, ii+iao, jj+jao) = dpint(k, ii+iao, jj+jao) &
                              + dtmpj(k)
                        end do

                        do k = 1, 6
                           !$omp atomic
                           qpint(k, ii+iao, jj+jao) = qpint(k, ii+iao, jj+jao) &
                              + qtmpj(k)
                        end do

                        !$omp atomic
                        hamiltonian(ii+iao, jj+jao) = hamiltonian(ii+iao, jj+jao) &
                           + stmp(ij) * hij + hdij
                     end if
                  end do
               end do

            end do
         end do

         ! Get the overlap integrals for optional scaled H0 basis set
         if (do_scaled_h0) then
            do ish = 1, bas%nsh_id(izp)
               iaosh = smap(ish-1)
               do jsh = 1, bas%nsh_id(jzp)
                  jaosh = smap(jsh-1)

                  ! Use modified H0 basis functions to recalculate the overlap
                  call overlap_cgto(bas%cgto_h0(jsh, jzp)%raw, &
                     & bas%cgto_h0(ish, izp)%raw, bcache%cgto_h0(jsh, jat), &
                     & bcache%cgto_h0(ish, iat), r2, vec, bas%intcut, stmp)
                  
                  ! Distribute shell block to diatomic overlap matrix
                  nao = msao(bas%cgto(jsh, jzp)%raw%ang)
                  do iao = 1, msao(bas%cgto(ish, izp)%raw%ang)
                     do jao = 1, nao
                        ij = jao + nao*(iao-1)
                        block_overlap(jaosh+jao, iaosh+iao) = stmp(ij)
                     end do
                  end do

               end do
            end do
         end if

         ! Optional diatomic frame transformation and scaling of the overlap
         if (h0%do_diat_scale) then
            call diat_trafo(block_overlap, vec, h0%ksig(izp,jzp), h0%kpi(izp,jzp), &
               &  h0%kdel(izp,jzp), bas%nsh_at(jat)-1, bas%nsh_at(iat)-1)
         end if

         ! Optional repeated setup of the Hamiltonian after modifications
         if (modify_h0) then
            do ish = 1, bas%nsh_id(izp)
               ii = bas%iao_sh(is+ish)
               iaosh = smap(ish-1)
               do jsh = 1, bas%nsh_id(jzp)
                  jj = bas%iao_sh(js+jsh)
                  jaosh = smap(jsh-1)

                  ! Recalculate shell polynomials
                  shpolyi = 1.0_wp + h0%shpoly(ish, izp)*rr &
                     & + h0%shpoly2(ish, izp)*rr2 + h0%shpoly4(ish, izp)*rr4
                  shpolyj = 1.0_wp + h0%shpoly(jsh, jzp)*rr &
                     & + h0%shpoly2(jsh, jzp)*rr2 + h0%shpoly4(jsh, jzp)*rr4
                  shpoly = shpolyi * shpolyj

                  ! Recalculate Hamiltonian elements for modified H0
                  hij = 0.5_wp * mod_h0_fraction &
                     & * (selfenergy(is+ish) + selfenergy(js+jsh)) &
                     & * h0%hscale(jsh, ish, jzp, izp) * shpoly

                  ! Distribute shell block to Hamiltonian matrix
                  nao = msao(bas%cgto(jsh, jzp)%raw%ang)
                  do iao = 1, msao(bas%cgto(ish, izp)%raw%ang)
                     do jao = 1, nao
                        ij = jao + nao*(iao-1)

                        ! Add modified Hamiltonian contribution
                        !$omp atomic
                        hamiltonian(jj+jao, ii+iao) = hamiltonian(jj+jao, ii+iao) &
                           + block_overlap(jaosh+jao, iaosh+iao) * hij

                        if (iat /= jat) then
                           !$omp atomic
                           hamiltonian(ii+iao, jj+jao) = hamiltonian(ii+iao, jj+jao) &
                              + block_overlap(jaosh+jao, iaosh+iao) * hij
                        end if
                     end do
                  end do

               end do
            end do
         end if
      end do
   end do

   !$omp parallel do schedule(runtime) default(none) &
   !$omp shared(mol, bas, bcache, overlap, dpint, qpint, hamiltonian) &
   !$omp shared(h0, selfenergy, aniso_dip) &
   !$omp private(iat, izp, is, ish, jsh, ii, jj, iao, jao, nao, ij) &
   !$omp private(vec, stmp, dtmpi, qtmpi, hij, hdij)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = bas%ish_at(iat)
      vec(:) = 0.0_wp
      do ish = 1, bas%nsh_id(izp)
         ii = bas%iao_sh(is+ish)
         do jsh = 1, bas%nsh_id(izp)
            jj = bas%iao_sh(is+jsh)

            ! Calculate onsite overlap and multipole integrals
            call multipole_cgto(bas%cgto(jsh, izp)%raw, bas%cgto(ish, izp)%raw, &
               & bcache%cgto(jsh, iat), bcache%cgto(ish, iat), &
               & 0.0_wp, vec, bas%intcut, stmp, dtmpi, qtmpi)

            ! The shell polynomial is unity, because rr is always zero
            hij = 0.5_wp * (selfenergy(is+ish) + selfenergy(is+jsh))

            nao = msao(bas%cgto(jsh, izp)%raw%ang)
            do iao = 1, msao(bas%cgto(ish, izp)%raw%ang)
               do jao = 1, nao
                  ij = jao + nao*(iao-1)
                  overlap(jj+jao, ii+iao) = overlap(jj+jao, ii+iao) &
                     + stmp(ij)

                  dpint(:, jj+jao, ii+iao) = dpint(:, jj+jao, ii+iao) &
                     + dtmpi(:, ij)
                  ! Construct onsite dipole Hamiltonian correction
                  hdij = - dot_product(aniso_dip(:, iat), dtmpi(:, ij))

                  qpint(:, jj+jao, ii+iao) = qpint(:, jj+jao, ii+iao) &
                     + qtmpi(:, ij)

                  hamiltonian(jj+jao, ii+iao) = hamiltonian(jj+jao, ii+iao) &
                     + stmp(ij) * hij + hdij
               end do
            end do

         end do
      end do
   end do

end subroutine get_hamiltonian


subroutine get_hamiltonian_gradient(mol, trans, list, bas, bcache, h0, &
   & selfenergy, dsedcn, aniso_dip, pot, pmat, xmat, &
   & dEdcn, dEdad, dEdcnbas, dEdqbas, gradient, sigma)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Lattice points within a given realspace cutoff
   real(wp), intent(in) :: trans(:, :)
   !> Neighbour list
   type(adjacency_list), intent(in) :: list
   !> Basis set information
   class(basis_type), intent(in) :: bas
   !> Basis set cache
   type(basis_cache), intent(in) :: bcache
   !> Hamiltonian interaction data
   type(tb_hamiltonian), intent(in) :: h0
   !> Diagonal elememts of the Hamiltonian
   real(wp), intent(in) :: selfenergy(:)
   !> Derivative of the diagonal elements of the Hamiltonian w.r.t. the coordination number
   real(wp), intent(in) :: dsedcn(:)
   !> Anisotropy vectors for dipole Hamiltonian correction
   real(wp), intent(in) :: aniso_dip(:, :)
   !> Density dependent potential shifts on the Hamiltonian
   type(potential_type), intent(in) :: pot
   !> Density matrix
   real(wp), intent(in) :: pmat(:, :, :)
   !> Energy weighted density matrix
   real(wp), intent(in) :: xmat(:, :, :)

   !> Derivative of the electronic energy w.r.t. the coordination number
   real(wp), intent(inout) :: dEdcn(:)
   !> Derivative of the electronic energy w.r.t. the anisotropy vector
   real(wp), intent(inout) :: dEdad(:, :)
   !> Derivative of the electronic energy w.r.t. basis set coordination numbers
   real(wp), intent(inout) :: dEdcnbas(:)
   !> Derivative of the electronic energy w.r.t. basis set atomic charges
   real(wp), intent(inout) :: dEdqbas(:)
   !> Derivative of the electronic energy w.r.t. coordinate displacements
   real(wp), intent(inout) :: gradient(:, :)
   !> Derivative of the electronic energy w.r.t. strain deformations
   real(wp), intent(inout) :: sigma(:, :)

   integer :: iat, jat, izp, jzp, itr, img, inl, spin, nspin, max_nsh
   integer :: ish, jsh, is, js, ii, jj, iao, jao, iaosh, jaosh, nao, ij
   real(wp) :: r, r2, rr, rr2, rr4, vec(3), hij, dtmpi(3), qtmpi(6), dG(3), tmp
   real(wp) :: ddtmpidqeffi(3), ddtmpidqeffj(3), dqtmpidqeffi(6), dqtmpidqeffj(6)
   real(wp) :: shpolyi, shpolyj, shpoly,  dshpolyi, dshpolyj, dshpoly, dsv(3)
   real(wp) :: hscale, hs, sval, dcni, dcnj, dhdcni, dhdcnj, hpij, pij, xij
   real(wp) :: dadi(3), dadj(3), daqi(6), daqj(6), mod_h0_fraction
   real(wp) :: dcnbasi, dcnbasj, dqbasi, dqbasj
   real(wp), allocatable :: stmp(:), dtmpj(:, :), qtmpj(:, :)
   real(wp), allocatable :: dstmp(:, :), ddtmpi(:, :, :), dqtmpi(:, :, :)
   real(wp), allocatable :: ddtmpj(:, :, :), dqtmpj(:, :, :)
   real(wp), allocatable :: block_overlap(:, :), block_doverlap(:, :, :)
   real(wp), allocatable :: block_doverlapdqeffi(:, :), block_doverlapdqeffj(:, :)
   real(wp), allocatable :: dstmpdqeffi(:), dstmpdqeffj(:)
   real(wp), allocatable :: ddtmpjdqeffi(:, :), ddtmpjdqeffj(:, :)
   real(wp), allocatable :: dqtmpjdqeffi(:, :), dqtmpjdqeffj(:, :)
   logical :: do_scaled_h0, modify_h0, compute_qeff_grad

   ! Thread-private array for reduction
   ! Set to 0 explicitly as the shared variants are potentially non-zero (inout)
   real(wp), allocatable :: dEdcn_local(:), dEdad_local(:, :), &
      & dEdcnbas_local(:), dEdqbas_local(:), gradient_local(:, :), sigma_local(:, :)

   ! Determine before the loop if effective charge derivatives are needed
   compute_qeff_grad = bas%charge_dependent

   nspin = size(pmat, 3)

   ! Predefine loop-invariant logicals to control H0 scaling
   ! and diatomic frame transformation 
   do_scaled_h0 = allocated(bas%cgto_h0)
   modify_h0 = do_scaled_h0 .or. h0%do_diat_scale

   ! Select if we construct the Hamiltonian with modifications
   mod_h0_fraction = 0.0_wp
   if(modify_h0) then
      mod_h0_fraction = 1.0_wp 
   end if

   allocate(stmp(msao(bas%maxl)**2), dstmp(3, msao(bas%maxl)**2), &
      & dtmpj(3, msao(bas%maxl)**2), ddtmpi(3, 3, msao(bas%maxl)**2), &
      & qtmpj(6, msao(bas%maxl)**2), dqtmpi(3, 6, msao(bas%maxl)**2), &
      & ddtmpj(3, 3, msao(bas%maxl)**2), dqtmpj(3, 6, msao(bas%maxl)**2))

   if (compute_qeff_grad) then
      allocate(dstmpdqeffi(msao(bas%maxl)**2), dstmpdqeffj(msao(bas%maxl)**2), &
         & ddtmpjdqeffi(3, msao(bas%maxl)**2), ddtmpjdqeffj(3, msao(bas%maxl)**2), &
         & dqtmpjdqeffi(6, msao(bas%maxl)**2), dqtmpjdqeffj(6, msao(bas%maxl)**2))
   end if

   ! Temporary overlap block between two atoms with maximum size for a single-zeta basis
   max_nsh = maxval(bas%nsh_id, 1)
   allocate(block_overlap(smap(max_nsh), smap(max_nsh)), &
      & block_doverlap(3, smap(max_nsh), smap(max_nsh)))
   if (compute_qeff_grad) then
      allocate(block_doverlapdqeffi(smap(max_nsh), smap(max_nsh)), &
         & block_doverlapdqeffj(smap(max_nsh), smap(max_nsh)))
   end if

   !$omp parallel default(none) shared(dEdcn, dEdad, dEdcnbas, dEdqbas) &
   !$omp shared(gradient, sigma, nspin, mol, bas, bcache, trans, h0, selfenergy) &
   !$omp shared(dsedcn, pot, pmat, xmat, list, aniso_dip) &
   !$omp firstprivate(do_scaled_h0, modify_h0, mod_h0_fraction, compute_qeff_grad) &
   !$omp private(iat, jat, izp, jzp, itr, is, js, ish, jsh, ii, jj, iaosh, jaosh) &
   !$omp private(iao, jao, nao, ij, inl, img, spin,  vec, r2, r, rr, rr2, rr4, stmp) &
   !$omp private(dstmp, dtmpi, dtmpj, ddtmpi, ddtmpj, qtmpi, qtmpj, dqtmpi, dqtmpj) &
   !$omp private(dstmpdqeffi, ddtmpjdqeffi, ddtmpidqeffi, dqtmpidqeffi, dqtmpjdqeffi) &
   !$omp private(dstmpdqeffj, ddtmpjdqeffj, ddtmpidqeffj, dqtmpidqeffj, dqtmpjdqeffj) &
   !$omp private(tmp, dcnbasi, dcnbasj, dqbasi, dqbasj, shpolyi, shpolyj, shpoly) &
   !$omp private(dshpolyi, dshpolyj, dshpoly, dsv, hscale, hs, hij, dhdcni, dhdcnj) &
   !$omp private(dG, dcni, dcnj, pij, xij, sval, hpij, dadi, dadj, daqi, daqj) &
   !$omp private(block_overlap, block_doverlap, block_doverlapdqeffi, block_doverlapdqeffj) &
   !$omp private(dEdcn_local, dEdad_local, dEdcnbas_local, dEdqbas_local) & 
   !$omp private(gradient_local, sigma_local)
   allocate(dEdcn_local(size(dEdcn)), source=0.0_wp)
   allocate(dEdad_local(size(dEdad,1), size(dEdad,2)), source=0.0_wp)
   if (compute_qeff_grad) then
      allocate(dEdcnbas_local(size(dEdcnbas)), source=0.0_wp)
      allocate(dEdqbas_local(size(dEdqbas)), source=0.0_wp)
   end if
   allocate(gradient_local(size(gradient,1), size(gradient,2)), source=0.0_wp)
   allocate(sigma_local(size(sigma,1), size(sigma,2)), source=0.0_wp)
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = bas%ish_at(iat)
      inl = list%inl(iat)
      do img = 1, list%nnl(iat)
         jat = list%nlat(img+inl)
         itr = list%nltr(img+inl)
         jzp = mol%id(jat)
         js = bas%ish_at(jat)
         if (iat == jat) cycle
         vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
         r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
         r = sqrt(r2)
         ! Precompute square root, linear, and quadratic distance factors
         rr2 = r / (h0%rad(jzp) + h0%rad(izp))
         rr = sqrt(rr2)
         rr4 = rr2 * r
         ! Reset all temporary accumulating variables
         block_overlap(:, :) = 0.0_wp
         block_doverlap(:, :, :) = 0.0_wp
         if (compute_qeff_grad) then
            block_doverlapdqeffi(:,:) = 0.0_wp
            block_doverlapdqeffj(:,:) = 0.0_wp
         end if
         dG(:) = 0.0_wp
         dcni = 0.0_wp
         dcnj = 0.0_wp
         dadi(:) = 0.0_wp
         dadj(:) = 0.0_wp
         daqi(:) = 0.0_wp
         daqj(:) = 0.0_wp
         dcnbasi = 0.0_wp
         dcnbasj = 0.0_wp
         dqbasi = 0.0_wp
         dqbasj = 0.0_wp
         ! Gradient for current diatomic pair
         do ish = 1, bas%nsh_id(izp)
            ii = bas%iao_sh(is+ish)
            iaosh = smap(ish-1)
            do jsh = 1, bas%nsh_id(jzp)
               jj = bas%iao_sh(js+jsh)
               jaosh = smap(jsh-1)

               ! Calculate overlap and multipole integral derivatives
               ! Guage-origin of the dipole/quadrupole integrals is at atom j
               call multipole_grad_cgto(bas%cgto(jsh, jzp)%raw, bas%cgto(ish, izp)%raw, &
                  & bcache%cgto(jsh, jat), bcache%cgto(ish, iat), r2, vec, bas%intcut, &
                  & stmp, dtmpj, qtmpj, dstmp, ddtmpj, dqtmpj, ddtmpi, dqtmpi, &
                  & dstmpdqeffj, dstmpdqeffi, ddtmpjdqeffj, ddtmpjdqeffi, &
                  & dqtmpjdqeffj, dqtmpjdqeffi)

               ! Calculate shell polynomials and their derivatives
               shpolyi = 1.0_wp + h0%shpoly(ish, izp)*rr &
                  & + h0%shpoly2(ish, izp)*rr2 + h0%shpoly4(ish, izp)*rr4
               dshpolyi = (0.5_wp * h0%shpoly(ish, izp) * rr &
                  & + h0%shpoly2(ish, izp) * rr2 &
                  & + 2.0_wp * h0%shpoly4(ish, izp) * rr4) / r2
               shpolyj = 1.0_wp + h0%shpoly(jsh, jzp)*rr &
                  & + h0%shpoly2(jsh, jzp)*rr2 + h0%shpoly4(jsh, jzp)*rr4
               dshpolyj = (0.5_wp * h0%shpoly(jsh, jzp) * rr &
                  & + h0%shpoly2(jsh, jzp) * rr2 &
                  & + 2.0_wp * h0%shpoly4(jsh, jzp) * rr4) / r2
               shpoly = shpolyi * shpolyj
               dshpoly = dshpolyi * shpolyj + dshpolyj * shpolyi
               dsv(:) = dshpoly / shpoly * vec

               ! Calculate Hamiltonian elements
               hscale = h0%hscale(jsh, ish, jzp, izp)
               hs = hscale * shpoly * (1.0_wp - mod_h0_fraction)
               hij = 0.5_wp * (selfenergy(is+ish) + selfenergy(js+jsh)) * hs 

               ! Calculate Hamiltonian derivatives w.r.t. the CN
               dhdcni = dsedcn(is+ish) * hs
               dhdcnj = dsedcn(js+jsh) * hs

               ! Contract derivatives density matrix for the current shell pair
               nao = msao(bas%cgto(jsh, jzp)%raw%ang)
               do iao = 1, msao(bas%cgto(ish, izp)%raw%ang)
                  do jao = 1, nao
                     ij = jao + nao*(iao-1)
                     ! Shift the atom j centered multipole integrals back to atom i
                     call shift_operator(-vec, stmp(ij), dtmpj(:, ij), qtmpj(:, ij), &
                        & dtmpi, qtmpi)
                     if (compute_qeff_grad) then
                        call shift_operator(-vec, dstmpdqeffj(ij), ddtmpjdqeffj(:, ij), &
                           & dqtmpjdqeffj(:, ij), ddtmpidqeffj, dqtmpidqeffj)
                        call shift_operator(-vec, dstmpdqeffi(ij), ddtmpjdqeffi(:, ij), &
                           & dqtmpjdqeffi(:, ij), ddtmpidqeffi, dqtmpidqeffi)
                     end if

                     ! Overlap and derivatives for possible diatomic frame trafo
                     block_overlap(jaosh+jao, iaosh+iao) = stmp(ij)
                     block_doverlap(:, jaosh+jao, iaosh+iao) = dstmp(:, ij)
                     if (compute_qeff_grad) then
                        block_doverlapdqeffi(jaosh+jao, iaosh+iao) = dstmpdqeffi(ij)
                        block_doverlapdqeffj(jaosh+jao, iaosh+iao) = dstmpdqeffj(ij)
                     end if

                     ! Collect all prefactors for the overlap derivative
                     sval = 0.0_wp
                     do spin = 1, nspin
                        pij = pmat(jj+jao, ii+iao, spin)
                        ! Mulliken partition gradient: (V_i + V_j) * P_ij * dS_ij
                        sval = sval - pij * (pot%vao(jj+jao, spin) &
                           & + pot%vao(ii+iao, spin))
                     end do
                     ! Contributions from the total density matrix
                     pij = pmat(jj+jao, ii+iao, 1)
                     hpij = pij * hij
                     ! H0 and energy weighted density matrix for both spins
                     sval = sval + 2*hpij - 2*xmat(jj+jao, ii+iao, 1)

                     ! Accumulate overlap, shell-polynomial, and multipole-gradient
                     dG(:) = dG + sval * dstmp(:, ij) &
                        & + 2*hpij*stmp(ij) * dsv &
                        & - pij * matmul(ddtmpi(:, :, ij), aniso_dip(:, iat)) &
                        & - pij * matmul(ddtmpj(:, :, ij), aniso_dip(:, jat)) &
                        & - pij * matmul(ddtmpi(:, :, ij), pot%vdp(:, iat, 1)) &
                        & - pij * matmul(ddtmpj(:, :, ij), pot%vdp(:, jat, 1)) &
                        & - pij * matmul(dqtmpi(:, :, ij), pot%vqp(:, iat, 1)) &
                        & - pij * matmul(dqtmpj(:, :, ij), pot%vqp(:, jat, 1))

                     ! Accumulate effective charge derivatives from overlap and multipoles
                     if (compute_qeff_grad) then
                        tmp = (sval * dstmpdqeffi(ij) & 
                           & - pij * dot_product(ddtmpidqeffi(:), aniso_dip(:, iat)) &
                           & - pij * dot_product(ddtmpjdqeffi(:, ij), aniso_dip(:, jat)) &
                           & - pij * dot_product(ddtmpidqeffi(:), pot%vdp(:, iat, 1)) &
                           & - pij * dot_product(ddtmpjdqeffi(:, ij), pot%vdp(:, jat, 1)) &
                           & - pij * dot_product(dqtmpidqeffi(:), pot%vqp(:, iat, 1)) &
                           & - pij * dot_product(dqtmpjdqeffi(:, ij), pot%vqp(:, jat, 1)))
                        dqbasi = dqbasi + bcache%cgto(ish, iat)%dqeffdq * tmp
                        dcnbasi = dcnbasi + bcache%cgto(ish, iat)%dqeffdcn * tmp

                        tmp = (sval * dstmpdqeffj(ij) & 
                           & - pij * dot_product(ddtmpidqeffj(:), aniso_dip(:, iat)) &
                           & - pij * dot_product(ddtmpjdqeffj(:, ij), aniso_dip(:, jat)) &
                           & - pij * dot_product(ddtmpidqeffj(:), pot%vdp(:, iat, 1)) &
                           & - pij * dot_product(ddtmpjdqeffj(:, ij), pot%vdp(:, jat, 1)) &
                           & - pij * dot_product(dqtmpidqeffj(:), pot%vqp(:, iat, 1)) &
                           & - pij * dot_product(dqtmpjdqeffj(:, ij), pot%vqp(:, jat, 1)))
                        dqbasj = dqbasj + bcache%cgto(jsh, jat)%dqeffdq * tmp
                        dcnbasj = dcnbasj + bcache%cgto(jsh, jat)%dqeffdcn * tmp
                     end if

                     ! Accumulate CN derivatives
                     dcni = dcni + dhdcni * pij * stmp(ij)
                     dcnj = dcnj + dhdcnj * pij * stmp(ij)

                     ! Accumulate anisotropy derivatives
                     dadi(:) = dadi - pij * dtmpi(:)
                     dadj(:) = dadj - pij * dtmpj(:, ij)
                     daqi(:) = daqi - pij * qtmpi(:)
                     daqj(:) = daqj - pij * qtmpj(:, ij)
                  end do
               end do

            end do
         end do

         ! Get the overlap integral derivatives for optional scaled H0 basis set 
         if (do_scaled_h0) then
            do ish = 1, bas%nsh_id(izp)
               iaosh = smap(ish-1)
               do jsh = 1, bas%nsh_id(jzp)
                  jaosh = smap(jsh-1)

                  ! Use modified H0 basis functions to recalculate the overlap
                  call overlap_grad_cgto(bas%cgto_h0(jsh, jzp)%raw, &
                     & bas%cgto_h0(ish, izp)%raw, bcache%cgto_h0(jsh, jat), &
                     & bcache%cgto_h0(ish, iat), r2, vec, bas%intcut, stmp, dstmp, &
                     & dstmpdqeffj, dstmpdqeffi)
                  
                  ! Distribute shell block to diatomic overlap matrix
                  nao = msao(bas%cgto(jsh, jzp)%raw%ang)
                  do iao = 1, msao(bas%cgto(ish, izp)%raw%ang)
                     do jao = 1, nao
                        ij = jao + nao*(iao-1)
                        block_overlap(jaosh+jao, iaosh+iao) = stmp(ij)
                        block_doverlap(:, jaosh+jao, iaosh+iao) = dstmp(:, ij)
                        if (compute_qeff_grad) then
                           block_doverlapdqeffj(jaosh+jao, iaosh+iao) = &
                              dstmpdqeffj(ij)
                           block_doverlapdqeffi(jaosh+jao, iaosh+iao) = &
                              dstmpdqeffi(ij)
                        end if
                     end do
                  end do

               end do
            end do
         end if

         ! Optional diatomic frame scaling transformation for overlap and derivatives
         if (h0%do_diat_scale) then
            call diat_trafo_grad(block_overlap, block_doverlap, vec, &
               & h0%ksig(izp, jzp), h0%kpi(izp, jzp), h0%kdel(izp, jzp), &
               & bas%nsh_at(jat)-1, bas%nsh_at(iat)-1)
            if (compute_qeff_grad) then
               call diat_trafo(block_doverlapdqeffi, vec, h0%ksig(izp, jzp), &
                  & h0%kpi(izp, jzp), h0%kdel(izp,jzp), bas%nsh_at(jat)-1, &
                  & bas%nsh_at(iat)-1)
               call diat_trafo(block_doverlapdqeffj, vec, h0%ksig(izp, jzp), &
                  & h0%kpi(izp, jzp), h0%kdel(izp, jzp), bas%nsh_at(jat)-1, &
                  & bas%nsh_at(iat)-1)
            end if
         end if

         ! Optional repeated setup of the Hamiltonian after modifications        
         if (modify_h0) then
            do ish = 1, bas%nsh_id(izp)
               ii = bas%iao_sh(is+ish)
               iaosh = smap(ish-1) 
               do jsh = 1, bas%nsh_id(jzp)
                  jj = bas%iao_sh(js+jsh)
                  jaosh = smap(jsh-1) 

                  ! Recalculate shell polynomials and their derivatives
                  shpolyi = 1.0_wp + h0%shpoly(ish, izp)*rr &
                     & + h0%shpoly2(ish, izp)*rr2 + h0%shpoly4(ish, izp)*rr4
                  dshpolyi = (0.5_wp * h0%shpoly(ish, izp) * rr &
                     & + h0%shpoly2(ish, izp) * rr2 &
                     & + 2.0_wp * h0%shpoly4(ish, izp) * rr4) / r2
                  shpolyj = 1.0_wp + h0%shpoly(jsh, jzp)*rr &
                     & + h0%shpoly2(jsh, jzp)*rr2 + h0%shpoly4(jsh, jzp)*rr4
                  dshpolyj = (0.5_wp * h0%shpoly(jsh, jzp) * rr &
                     & + h0%shpoly2(jsh, jzp) * rr2 &
                     & + 2.0_wp * h0%shpoly4(jsh, jzp) * rr4) / r2
                  shpoly = shpolyi * shpolyj
                  dshpoly = dshpolyi * shpolyj + dshpolyj * shpolyi
                  dsv(:) = dshpoly / shpoly * vec

                  ! Recalculate Hamiltonian elements
                  hscale = h0%hscale(jsh, ish, jzp, izp)
                  hs = hscale * shpoly * mod_h0_fraction
                  hij = 0.5_wp * (selfenergy(is+ish) + selfenergy(js+jsh)) * hs

                  ! Recalculate Hamiltonian derivatives w.r.t. the CN
                  dhdcni = dsedcn(is+ish) * hs
                  dhdcnj = dsedcn(js+jsh) * hs

                  ! Contract derivatives with density matrix for current shell pair
                  nao = msao(bas%cgto(jsh, jzp)%raw%ang)
                  do iao = 1, msao(bas%cgto(ish, izp)%raw%ang)
                     do jao = 1, nao
                        ij = jao + nao*(iao-1)

                        ! Add only H0 gradient with the modified overlap
                        pij = pmat(jj+jao, ii+iao, 1)
                        hpij = 2*pij * hij

                        ! Accumulate gradient from overlap and shell-polynomials
                        dG(:) = dG + hpij * (block_doverlap(:, jaosh+jao, iaosh+iao) &
                           + block_overlap(jaosh+jao, iaosh+iao) * dsv)

                        ! Accumulate effective charge derivatives for the modified H0
                        if (compute_qeff_grad) then
                           tmp = hpij * block_doverlapdqeffi(jaosh+jao, iaosh+iao)
                           dqbasi = dqbasi + bcache%cgto_h0(ish, iat)%dqeffdq * tmp
                           dcnbasi = dcnbasi + bcache%cgto_h0(ish, iat)%dqeffdcn * tmp

                           tmp = hpij * block_doverlapdqeffj(jaosh+jao, iaosh+iao)
                           dqbasj = dqbasj + bcache%cgto_h0(jsh, jat)%dqeffdq * tmp
                           dcnbasj = dcnbasj + bcache%cgto_h0(jsh, jat)%dqeffdcn * tmp
                        end if

                        ! Accumulate CN derivatives
                        dcni = dcni + dhdcni * pij * block_overlap(jaosh+jao, iaosh+iao)
                        dcnj = dcnj + dhdcnj * pij * block_overlap(jaosh+jao, iaosh+iao)
                     end do
                  end do

               end do
            end do
         end if

         ! Collect CN-, anisotropy, effective (H0) charge-derivatives, gradient
         ! and sigma once per atom pair
         dEdcn_local(iat) = dEdcn_local(iat) + dcni
         dEdcn_local(jat) = dEdcn_local(jat) + dcnj
         dEdad_local(:, iat) = dEdad_local(:, iat) + dadi
         dEdad_local(:, jat) = dEdad_local(:, jat) + dadj
         if (compute_qeff_grad) then
            dEdqbas_local(iat) = dEdqbas_local(iat) + dqbasi
            dEdqbas_local(jat) = dEdqbas_local(jat) + dqbasj
            dEdcnbas_local(iat) = dEdcnbas_local(iat) + dcnbasi
            dEdcnbas_local(jat) = dEdcnbas_local(jat) + dcnbasj
         end if
         gradient_local(:, iat) = gradient_local(:, iat) + dG
         gradient_local(:, jat) = gradient_local(:, jat) - dG
         sigma_local(:, :) = sigma_local + 0.5_wp * (spread(vec, 1, 3) &
            & * spread(dG, 2, 3) + spread(dG, 1, 3) * spread(vec, 2, 3))

      end do

      ! Onsite contributions
      vec(:) = 0.0_wp
      dadi(:) = 0.0_wp
      daqi(:) = 0.0_wp
      dcni = 0.0_wp
      dqbasi = 0.0_wp
      dcnbasi = 0.0_wp
      do ish = 1, bas%nsh_id(izp)
         ! Accumulate onsite anisotropy derivatives
         ii = bas%iao_sh(is+ish)
         do jsh = 1, bas%nsh_id(izp)
            jj = bas%iao_sh(is+jsh)

            ! Calculate only the onsite multipole integrals and derivatives
            call multipole_grad_cgto(bas%cgto(jsh, izp)%raw, bas%cgto(ish, izp)%raw, &
               & bcache%cgto(jsh, iat), bcache%cgto(ish, iat), 0.0_wp, vec, bas%intcut, &
               & stmp, dtmpj, qtmpj, dstmp, ddtmpj, dqtmpj, ddtmpi, dqtmpi, &
               & dstmpdqeffj, dstmpdqeffi, ddtmpjdqeffj, ddtmpjdqeffi, &
               & dqtmpjdqeffj, dqtmpjdqeffi)

            nao = msao(bas%cgto(jsh, izp)%raw%ang)
            do iao = 1, msao(bas%cgto(ish, izp)%raw%ang)
               do jao = 1, nao
                  ij = jao + nao*(iao-1)

                  pij = pmat(jj+jao, ii+iao, 1)

                  ! Accumulate onsite anisotropy derivatives
                  dadi(:) = dadi(:) - pij * dtmpj(:, ij)
                  daqi(:) = daqi(:) - pij * qtmpj(:, ij)

                  ! Accumulate onsite dipole/quadrupole effective charge derivatives
                  ! For the overlap this derivative is zero due to normalization
                  if (compute_qeff_grad) then
                     tmp = 2*pij * dot_product(ddtmpjdqeffi(:, ij), aniso_dip(:, iat))
                     dqbasi = dqbasi - bcache%cgto(ish, iat)%dqeffdq * tmp
                     dcnbasi = dcnbasi - bcache%cgto(ish, iat)%dqeffdcn * tmp

                     tmp = 2*pij * dot_product(ddtmpjdqeffi(:, ij), pot%vdp(:, iat, 1))
                     dqbasi = dqbasi - bcache%cgto(ish, iat)%dqeffdq * tmp
                     dcnbasi = dcnbasi - bcache%cgto(ish, iat)%dqeffdcn * tmp

                     tmp = 2*pij * dot_product(dqtmpjdqeffi(:, ij), pot%vqp(:, iat, 1))
                     dqbasi = dqbasi - bcache%cgto(ish, iat)%dqeffdq * tmp
                     dcnbasi = dcnbasi - bcache%cgto(ish, iat)%dqeffdcn * tmp
                  end if
               end do
            end do
         end do

         ! Accumulate onsite CN derivative
         dhdcni = dsedcn(is+ish)
         do iao = 1, msao(bas%cgto(ish, izp)%raw%ang)
            dcni = dcni + dhdcni * pmat(ii+iao, ii+iao, 1)
         end do
      end do
   
      ! Collect onsite CN- and anisotropy-derivatives once per atom
      dEdad_local(:, iat) = dEdad_local(:, iat) + dadi
      dEdcn_local(iat) = dEdcn_local(iat) + dcni
      if (compute_qeff_grad) then
         dEdqbas_local(iat) = dEdqbas_local(iat) + dqbasi
         dEdcnbas_local(iat) = dEdcnbas_local(iat) + dcnbasi
      end if
   end do
   !$omp end do
   !$omp critical (get_hamiltonian_gradient_)
   dEdcn(:) = dEdcn + dEdcn_local
   dEdad(:, :) = dEdad + dEdad_local
   if (compute_qeff_grad) then
      dEdqbas(:) = dEdqbas + dEdqbas_local
      dEdcnbas(:) = dEdcnbas + dEdcnbas_local
   end if
   gradient(:, :) = gradient + gradient_local
   sigma(:, :) = sigma + sigma_local
   !$omp end critical (get_hamiltonian_gradient_)
   if (compute_qeff_grad) then
      deallocate(dEdqbas_local, dEdcnbas_local)
   end if
   deallocate(dEdcn_local, dEdad_local, gradient_local, sigma_local)
   !$omp end parallel

end subroutine get_hamiltonian_gradient


subroutine get_occupation(mol, bas, h0, nocc, n0at, n0sh)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   class(basis_type), intent(in) :: bas
   !> Hamiltonian interaction data
   type(tb_hamiltonian), intent(in) :: h0
   !> Occupation number
   real(wp), intent(out) :: nocc
   !> Reference occupation for each atom
   real(wp), intent(out) :: n0at(:)
   !> Reference occupation for each shell
   real(wp), intent(out) :: n0sh(:)

   integer :: iat, ish, izp, ii

   nocc = -mol%charge
   n0at(:) = 0.0_wp
   n0sh(:) = 0.0_wp
   do iat = 1, mol%nat
      izp = mol%id(iat)
      ii = bas%ish_at(iat)
      do ish = 1, bas%nsh_id(izp)
         nocc = nocc + h0%refocc(ish, izp)
         n0at(iat) = n0at(iat) + h0%refocc(ish, izp)
         n0sh(ii+ish) = n0sh(ii+ish) + h0%refocc(ish, izp)
      end do
   end do

end subroutine get_occupation


subroutine get_number_electrons(mol, bas, h0, nel)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   class(basis_type), intent(in) :: bas
   !> Hamiltonian interaction data
   type(tb_hamiltonian), intent(in) :: h0
   !> Number of electrons
   real(wp), intent(out) :: nel

   integer :: iat, ish, izp

   nel = -mol%charge
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do ish = 1, bas%nsh_id(izp)
         nel = nel + h0%refocc(ish, izp)
      end do
   end do

end subroutine get_number_electrons


!> Shift multipole operator from Ket function (center i) to Bra function (center j),
!> the multipole operator on the Bra function can be assembled from the lower moments
!> on the Ket function and the displacement vector using horizontal shift rules.
pure subroutine shift_operator(vec, s, di, qi, dj, qj)
   !> Displacement vector of center i and j
   real(wp),intent(in) :: vec(:)
   !> Overlap integral between basis functions
   real(wp),intent(in) :: s
   !> Dipole integral with operator on Ket function (center i)
   real(wp),intent(in) :: di(:)
   !> Quadrupole integral with operator on Ket function (center i)
   real(wp),intent(in) :: qi(:)
   !> Dipole integral with operator on Bra function (center j)
   real(wp),intent(out) :: dj(:)
   !> Quadrupole integral with operator on Bra function (center j)
   real(wp),intent(out) :: qj(:)

   real(wp) :: tr

   ! Create dipole operator on Bra function from Ket function and shift contribution
   ! due to monopol displacement
   dj(1) = di(1) + vec(1)*s
   dj(2) = di(2) + vec(2)*s
   dj(3) = di(3) + vec(3)*s

   ! For the quadrupole operator on the Bra function we first construct the shift
   ! contribution from the dipole and monopol displacement, since we have to remove
   ! the trace contribution from the shift and the moment integral on the Ket function
   ! is already traceless
   qj(1) = 2*vec(1)*di(1) + vec(1)**2*s
   qj(3) = 2*vec(2)*di(2) + vec(2)**2*s
   qj(6) = 2*vec(3)*di(3) + vec(3)**2*s
   qj(2) = vec(1)*di(2) + vec(2)*di(1) + vec(1)*vec(2)*s
   qj(4) = vec(1)*di(3) + vec(3)*di(1) + vec(1)*vec(3)*s
   qj(5) = vec(2)*di(3) + vec(3)*di(2) + vec(2)*vec(3)*s
   ! Now collect the trace of the shift contribution
   tr = 0.5_wp * (qj(1) + qj(3) + qj(6))

   ! Finally, assemble the quadrupole operator on the Bra function from the operator
   ! on the Ket function and the traceless shift contribution
   qj(1) = qi(1) + 1.5_wp * qj(1) - tr
   qj(2) = qi(2) + 1.5_wp * qj(2)
   qj(3) = qi(3) + 1.5_wp * qj(3) - tr
   qj(4) = qi(4) + 1.5_wp * qj(4)
   qj(5) = qi(5) + 1.5_wp * qj(5)
   qj(6) = qi(6) + 1.5_wp * qj(6) - tr
end subroutine shift_operator

end module tblite_xtb_h0
