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

!> @file tblite/coulomb/charge/effective.f90
!> Provides an effective Coulomb operator for isotropic electrostatic interactions

!> Isotropic second-order electrostatics using an effective Coulomb operator
module tblite_coulomb_charge_effective
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_math, only : matdet_3x3, matinv_3x3
   use mctc_io_constants, only : pi
   use tblite_blas, only : dot, gemv, symv, gemm
   use tblite_coulomb_cache, only : coulomb_cache
   use tblite_coulomb_ewald, only : get_dir_cutoff, get_rec_cutoff
   use tblite_coulomb_charge_type, only : coulomb_charge_type
   use tblite_cutoff, only : get_lattice_points
   use tblite_scf_potential, only : potential_type
   use tblite_utils_average, only : average_type
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_wignerseitz, only : wignerseitz_cell
   implicit none
   private

   public :: new_effective_coulomb


   !> Effective, Klopman-Ohno-type, second-order electrostatics
   type, public, extends(coulomb_charge_type) :: effective_coulomb
      !> Averager for the hubbard parameters
      type(average_type), allocatable :: average
      !> CN-dependence of Hubbard parameter for each species
      real(wp), allocatable :: hubbard_cn(:)
      !> Exponent of radius dependent hubbard scaling
      real(wp) :: hubbard_exp
      !> Exponent of Coulomb kernel
      real(wp) :: gexp
      !> Long-range cutoff
      real(wp) :: rcut
   contains
      !> Evaluate Coulomb matrix
      procedure :: get_coulomb_matrix
      !> Evaluate uncontracted derivatives of Coulomb matrix
      procedure :: get_coulomb_derivs
   end type effective_coulomb

   real(wp), parameter :: twopi = 2 * pi
   real(wp), parameter :: sqrtpi = sqrt(pi)
   real(wp), parameter :: eps = sqrt(epsilon(0.0_wp))
   real(wp), parameter :: conv = eps
   real(wp), parameter :: cn_reg = 1e-6_wp
   character(len=*), parameter :: label = "isotropic Klopman-Ohno-Mataga-Nishimoto electrostatics"

contains

!> Construct new effective electrostatic interaction container
subroutine new_effective_coulomb(self, mol, gexp, hubbard, &
   & average, hubbard_cn, hubbard_exp, refqsh, nshell)
   !> Instance of the electrostatic container
   type(effective_coulomb), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Exponent of Coulomb kernel
   real(wp), intent(in) :: gexp
   !> Hubbard parameter for all shells and species
   real(wp), intent(in) :: hubbard(:, :)
   !> Averaging function for Hubbard parameter of a shell-pair
   type(average_type), intent(in) :: average
   !> CN-dependence of Hubbard parameters
   real(wp), intent(in), optional :: hubbard_cn(:)
   !> Exponent of radius-dependent hubbard scaling
   real(wp), intent(in), optional :: hubbard_exp
   !> Shell-resolved reference partial charges
   !> for first/pseudo-second-order offsite tight-binding
   real(wp), intent(in), optional :: refqsh(:)
   !> Number of shells for each species
   integer, intent(in), optional :: nshell(:)

   integer :: mshell
   integer :: isp, jsp, ish, jsh, ind, iat

   self%label = label

   self%shell_resolved = present(nshell)
   if (present(nshell)) then
      mshell = maxval(nshell)
      self%nshell = nshell(mol%id)
   else
      mshell = 1
      self%nshell = spread(1, 1, mol%nat)
   end if
   allocate(self%offset(mol%nat))
   ind = 0
   do iat = 1, mol%nat
      self%offset(iat) = ind
      ind = ind + self%nshell(iat)
   end do

   self%gexp = gexp
   self%rcut = 10.0_wp

   if (present(nshell)) then
      allocate(self%hubbard(mshell, mol%nid))
      do isp = 1, mol%nid
         self%hubbard(:, isp) = 0.0_wp
         do ish = 1, nshell(isp)
            self%hubbard(ish, isp) = hubbard(ish, isp)
         end do
      end do
   else
      allocate(self%hubbard(1, mol%nid))
      do isp = 1, mol%nid
         self%hubbard(1, isp) = hubbard(1, isp)
      end do
   end if

   if (present(hubbard_cn)) then
      self%hubbard_cn = hubbard_cn
      self%cn_dep = .true.
   else
      allocate(self%hubbard_cn(mol%nid))
      self%hubbard_cn = 0.0_wp
      self%cn_dep = .false.
   end if

   if (present(hubbard_exp)) then
      self%hubbard_exp = hubbard_exp
   else
      self%hubbard_exp = 0.0_wp
   end if

   self%average = average

   if (present(refqsh)) then
      self%refqsh = refqsh
   end if

end subroutine new_effective_coulomb


!> Evaluate coulomb matrix
subroutine get_coulomb_matrix(self, mol, cache, amat)
   !> Instance of the electrostatic container
   class(effective_coulomb), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(coulomb_cache), intent(inout) :: cache
   !> Coulomb matrix
   real(wp), contiguous, intent(out) :: amat(:, :)

   amat(:, :) = 0.0_wp

   if (any(mol%periodic)) then
      call get_amat_3d(mol, self%nshell, self%offset, self%hubbard, &
         & self%hubbard_cn, self%hubbard_exp, self%average, self%gexp, &
         & cache%cn, self%rcut, cache%wsc, cache%alpha, amat)
   else
      call get_amat_0d(mol, self%nshell, self%offset, self%hubbard, &
         & self%hubbard_cn, self%hubbard_exp, self%average, self%gexp, & 
         & cache%cn, amat)
   end if

end subroutine get_coulomb_matrix


!> Get real lattice vectors
subroutine get_dir_trans(lattice, alpha, conv, trans)
   !> Lattice parameters
   real(wp), intent(in) :: lattice(:, :)
   !> Parameter for Ewald summation
   real(wp), intent(in) :: alpha
   !> Tolerance for Ewald summation
   real(wp), intent(in) :: conv
   !> Translation vectors
   real(wp), allocatable, intent(out) :: trans(:, :)

   call get_lattice_points([.true.], lattice, get_dir_cutoff(alpha, conv), trans)

end subroutine get_dir_trans

!> Get reciprocal lattice translations
subroutine get_rec_trans(lattice, alpha, volume, conv, trans)
   !> Lattice parameters
   real(wp), intent(in) :: lattice(:, :)
   !> Parameter for Ewald summation
   real(wp), intent(in) :: alpha
   !> Cell volume
   real(wp), intent(in) :: volume
   !> Tolerance for Ewald summation
   real(wp), intent(in) :: conv
   !> Translation vectors
   real(wp), allocatable, intent(out) :: trans(:, :)

   real(wp) :: rec_lat(3, 3)

   rec_lat = twopi*transpose(matinv_3x3(lattice))
   call get_lattice_points([.true.], rec_lat, get_rec_cutoff(alpha, volume, conv), trans)
   trans = trans(:, 2:)

end subroutine get_rec_trans


!> Evaluate Coulomb matrix for finite systems
subroutine get_amat_0d(mol, nshell, offset, hubbard, hubbard_cn, hubbard_exp, &
   & average, gexp, cn, amat)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Number of shells for each atom
   integer, intent(in) :: nshell(:)
   !> Index offset for each shell
   integer, intent(in) :: offset(:)
   !> Hubbard parameter for each shell
   real(wp), intent(in) :: hubbard(:, :)
   !> CN-dependence of Hubbard parameters
   real(wp), intent(in) :: hubbard_cn(:)
   !> Exponent of radius-dependent hubbard scaling
   real(wp), intent(in) :: hubbard_exp
   !> Averaging function for Hubbard parameter of a shell-pair
   type(average_type), intent(in) :: average
   !> Exponent of Coulomb kernel
   real(wp), intent(in) :: gexp
   !> Coordination number
   real(wp), intent(in) :: cn(:)
   !> Coulomb matrix
   real(wp), intent(inout) :: amat(:, :)

   integer :: iat, jat, izp, jzp, ii, jj, ish, jsh
   real(wp) :: sqrtcni, sqrtcnj, scalei, scalej
   real(wp) :: vec(3), r1, r1g, gam, hubi, hubj, damp, tmp

   !$omp parallel do default(none) schedule(runtime) shared(amat, mol, cn) &
   !$omp shared(nshell, offset, hubbard, hubbard_cn, hubbard_exp, average, gexp) &
   !$omp private(iat, izp, ii, ish, jat, jzp, jj, jsh, scalei, scalej) &
   !$omp private(sqrtcni, sqrtcnj, vec, r1, r1g, damp, hubi, hubj, gam, tmp)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      ii = offset(iat)
      ! Coordination number scaling witz regularization
      sqrtcni = sqrt(cn(iat) + cn_reg**2) - cn_reg
      scalei = 1.0_wp + hubbard_cn(izp) * sqrtcni
      do jat = 1, iat-1
         jzp = mol%id(jat)
         jj = offset(jat)
         ! Coordination number scaling with regularization
         sqrtcnj = sqrt(cn(jat) + cn_reg**2) - cn_reg
         scalej = 1.0_wp + hubbard_cn(jzp) * sqrtcnj
         vec = mol%xyz(:, jat) - mol%xyz(:, iat)
         r1 = norm2(vec)
         r1g = r1**gexp
         ! Damping factor for distance-dependent hubbard scaling
         damp = exp (-hubbard_exp * r1)
         do ish = 1, nshell(iat)
            ! CN dependent hubbard parameter
            hubi = hubbard(ish, izp) * scalei
            do jsh = 1, nshell(jat)
               ! CN dependent hubbard parameter
               hubj = hubbard(jsh, jzp) * scalej
               ! Damped average hubbard parameter
               gam = average%value(hubi, hubj) / damp
               ! Effective Coulomb interaction
               tmp = 1.0_wp/(r1g + gam**(-gexp))**(1.0_wp/gexp)
               !$omp atomic
               amat(jj+jsh, ii+ish) = amat(jj+jsh, ii+ish) + tmp
               !$omp atomic
               amat(ii+ish, jj+jsh) = amat(ii+ish, jj+jsh) + tmp
            end do
         end do
      end do
      do ish = 1, nshell(iat)
         hubi = hubbard(ish, izp) * scalei
         do jsh = 1, ish-1
            hubj = hubbard(jsh, izp) * scalei
            ! Diagonal average hubbard parameter
            gam = average%value(hubi, hubj)
            !$omp atomic
            amat(ii+jsh, ii+ish) = amat(ii+jsh, ii+ish) + gam
            !$omp atomic
            amat(ii+ish, ii+jsh) = amat(ii+ish, ii+jsh) + gam
         end do
         !$omp atomic
         amat(ii+ish, ii+ish) = amat(ii+ish, ii+ish) + hubi
      end do
   end do

end subroutine get_amat_0d

!> Evaluate the coulomb matrix for 3D systems
subroutine get_amat_3d(mol, nshell, offset, hubbard, hubbard_cn, hubbard_exp, &
   & average, gexp, cn, rcut, wsc, alpha, amat)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Number of shells per atom
   integer, intent(in) :: nshell(:)
   !> Index offset for each atom
   integer, intent(in) :: offset(:)
   !> Hubbard parameter for each shell
   real(wp), intent(in) :: hubbard(:, :)
   !> CN-dependence of Hubbard parameters
   real(wp), intent(in) :: hubbard_cn(:)
   !> Exponent of radius-dependent hubbard scaling
   real(wp), intent(in) :: hubbard_exp
   !> Averaging function for Hubbard parameter of a shell-pair
   type(average_type), intent(in) :: average
   !> Exponent of the interaction kernel
   real(wp), intent(in) :: gexp
   !> Coordination number
   real(wp), intent(in) :: cn(:)
   !> Long-range cutoff
   real(wp), intent(in) :: rcut
   !> Wigner-Seitz cell
   type(wignerseitz_cell), intent(in) :: wsc
   !> Convergence factor
   real(wp), intent(in) :: alpha
   !> Coulomb matrix
   real(wp), intent(inout) :: amat(:, :)

   integer :: iat, jat, izp, jzp, img, ii, jj, ish, jsh
   real(wp) :: sqrtcni, sqrtcnj, scalei, scalej
   real(wp) :: vec(3), gam, hubi, hubj, wsw, dtmp, rtmp, vol, aval
   real(wp), allocatable :: dtrans(:, :), rtrans(:, :)

   vol = abs(matdet_3x3(mol%lattice))
   call get_dir_trans(mol%lattice, alpha, conv, dtrans)
   call get_rec_trans(mol%lattice, alpha, vol, conv, rtrans)

   !$omp parallel do default(none) schedule(runtime) shared(amat) &
   !$omp shared(mol, nshell, offset, hubbard, hubbard_cn, hubbard_exp) &
   !$omp shared(average, gexp, cn, wsc, dtrans, rtrans, alpha, vol, rcut) &
   !$omp private(iat, izp, jat, jzp, ii, jj, ish, jsh, scalei, scalej) &
   !$omp private(sqrtcni, sqrtcnj, hubi, hubj, gam, wsw, vec, dtmp, rtmp, aval)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      ii = offset(iat)
      ! Coordination number scaling with regularization
      sqrtcni = sqrt(cn(iat) + cn_reg**2) - cn_reg
      scalei = 1.0_wp + hubbard_cn(izp) * sqrtcni
      do jat = 1, iat-1
         jzp = mol%id(jat)
         jj = offset(jat)
         ! Coordination number scaling with regularization
         sqrtcnj = sqrt(cn(jat) + cn_reg**2) - cn_reg
         scalej = 1.0_wp + hubbard_cn(jzp) * sqrtcnj
         wsw = 1.0_wp / real(wsc%nimg(jat, iat), wp)
         do img = 1, wsc%nimg(jat, iat)
            vec = mol%xyz(:, iat) - mol%xyz(:, jat) - wsc%trans(:, wsc%tridx(img, jat, iat))
            call get_amat_rec_3d(vec, vol, alpha, 0.0_wp, rtrans, rtmp)
            do ish = 1, nshell(iat)
               hubi = hubbard(ish, izp) * scalei
               do jsh = 1, nshell(jat)
                  hubj = hubbard(jsh, jzp) * scalej
                  gam = average%value(hubi, hubj)
                  call get_amat_dir_3d(vec, gam, gexp, rcut, alpha, dtrans, dtmp)
                  aval = (dtmp + rtmp) * wsw
                  !$omp atomic
                  amat(jj+jsh, ii+ish) = amat(jj+jsh, ii+ish) + aval
                  !$omp atomic
                  amat(ii+ish, jj+jsh) = amat(ii+ish, jj+jsh) + aval
               end do
            end do
         end do
      end do

      wsw = 1.0_wp / real(wsc%nimg(iat, iat), wp)
      do img = 1, wsc%nimg(iat, iat)
         vec = wsc%trans(:, wsc%tridx(img, iat, iat))
         call get_amat_rec_3d(vec, vol, alpha, 0.0_wp, rtrans, rtmp)
         rtmp = rtmp - 2 * alpha / sqrtpi
         do ish = 1, nshell(iat)
            hubi = hubbard(ish, izp) * scalei
            do jsh = 1, ish-1
               hubj = hubbard(jsh, izp) * scalei
               gam = average%value(hubi, hubj)
               call get_amat_dir_3d(vec, gam, gexp, rcut, alpha, dtrans, dtmp)
               aval = (dtmp + rtmp + gam) * wsw
               !$omp atomic
               amat(ii+jsh, ii+ish) = amat(ii+jsh, ii+ish) + aval
               !$omp atomic
               amat(ii+ish, ii+jsh) = amat(ii+ish, ii+jsh) + aval
            end do
            gam = hubi
            call get_amat_dir_3d(vec, gam, gexp, rcut, alpha, dtrans, dtmp)
            aval = (dtmp + rtmp + gam) * wsw
            !$omp atomic
            amat(ii+ish, ii+ish) = amat(ii+ish, ii+ish) + aval
         end do
      end do

   end do

end subroutine get_amat_3d

!> Calculate real space contributions for a pair under 3D periodic boundary conditions
subroutine get_amat_dir_3d(rij, gam, gexp, rcut, alp, trans, amat)
   !> Distance between pair
   real(wp), intent(in) :: rij(3)
   !> Hubbard parameter
   real(wp), intent(in) :: gam
   !> Exponent for interaction kernel
   real(wp), intent(in) :: gexp
   !> Long-range cutoff
   real(wp), intent(in) :: rcut
   !> Convergence factor
   real(wp), intent(in) :: alp
   !> Translation vectors to consider
   real(wp), intent(in) :: trans(:, :)
   !> Interaction matrix element
   real(wp), intent(out) :: amat

   integer :: itr
   real(wp) :: vec(3), r1, tmp, fcut

   amat = 0.0_wp

   do itr = 1, size(trans, 2)
      vec(:) = rij + trans(:, itr)
      r1 = norm2(vec)
      if (r1 < eps) cycle
      fcut = fsmooth(r1, rcut)
      tmp = fcut/(r1**gexp + gam**(-gexp))**(1.0_wp/gexp) + (1.0_wp-fcut)/r1 &
         & - erf(alp*r1)/r1
      amat = amat + tmp
   end do

end subroutine get_amat_dir_3d

!> Calculate reciprocal space contributions for a pair under 3D periodic boundary conditions
subroutine get_amat_rec_3d(rij, vol, alp, qpc, trans, amat)
   !> Distance between pair
   real(wp), intent(in) :: rij(3)
   !> Volume of cell
   real(wp), intent(in) :: vol
   !> Convergence factor
   real(wp), intent(in) :: alp
   !> Pseudo-quadrupole contribution
   real(wp), intent(in) :: qpc
   !> Translation vectors to consider
   real(wp), intent(in) :: trans(:, :)
   !> Interaction matrix element
   real(wp), intent(out) :: amat

   integer :: itr
   real(wp) :: fac, vec(3), g2, gv, expk, cosk

   amat = 0.0_wp
   fac = 4*pi/vol

   do itr = 1, size(trans, 2)
      vec(:) = trans(:, itr)
      g2 = dot_product(vec, vec)
      if (g2 < eps) cycle
      gv = dot_product(rij, vec)
      expk = fac * exp(-0.25_wp*g2/(alp*alp))/g2
      cosk = cos(gv) * expk
      amat = amat + cosk
   end do

end subroutine get_amat_rec_3d

function fsmooth(r1, rcut) result(fcut)
   real(wp), intent(in) :: r1
   real(wp), intent(in) :: rcut
   real(wp) :: fcut

   real(wp), parameter :: offset = 1.0_wp
   real(wp), parameter :: c(*) = [-6.0_wp, 15.0_wp, -10.0_wp, 1.0_wp]
   real(wp) :: xrel

   if (r1 < rcut - offset) then
      fcut = 1.0_wp
   else if (r1 > rcut) then
      fcut = 0.0_wp
   else
      xrel = (r1 - (rcut - offset)) / offset
      fcut = c(1)*xrel**5 + c(2)*xrel**4 + c(3)*xrel**3 + c(4)
   end if
end function fsmooth

function dsmooth(r1, rcut) result(dcut)
   real(wp), intent(in) :: r1
   real(wp), intent(in) :: rcut
   real(wp) :: dcut

   real(wp), parameter :: offset = 1.0_wp
   real(wp), parameter :: c(*) = [-6.0_wp, 15.0_wp, -10.0_wp, 1.0_wp]
   real(wp) :: xrel

   if (r1 < rcut - offset .or. r1 > rcut) then
      dcut = 0.0_wp
   else
      xrel = (r1 - (rcut - offset)) / offset
      dcut = (5*c(1)*xrel**4 + 4*c(2)*xrel**3 + 3*c(3)*xrel**2) / offset
   end if

end function dsmooth


!> Evaluate uncontracted derivatives of Coulomb matrix
subroutine get_coulomb_derivs(self, mol, cache, qat, qsh, dadr, dadL, &
   & atrace, dadcn, datracedcn)
   !> Instance of the electrostatic container
   class(effective_coulomb), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Reusable data container
   type(coulomb_cache), intent(inout) :: cache
   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)
   !> Shell-resolved partial charges
   real(wp), intent(in) :: qsh(:)
   !> Derivative of interactions with respect to cartesian displacements
   real(wp), contiguous, intent(out) :: dadr(:, :, :)
   !> Derivative of interactions with respect to strain deformations
   real(wp), contiguous, intent(out) :: dadL(:, :, :)
   !> On-site derivatives with respect to cartesian displacements
   real(wp), contiguous, intent(out) :: atrace(:, :)
   !> Coulomb matrix derivative w.r.t. the coordination numbers
   real(wp), contiguous, intent(out) :: dadcn(:, :)
   !> On-site derivatives with respect to coordination numbers
   real(wp), contiguous, intent(out) :: datracedcn(:)

   if(self%shell_resolved) then
      if (any(mol%periodic)) then
         call get_damat_3d(mol, self%nshell, self%offset, self%hubbard, & 
            & self%hubbard_cn, self%hubbard_exp, self%average, self%gexp, &
            & cache%cn, self%rcut, cache%wsc, cache%alpha, qsh, &
            & dadr, dadL, atrace, dadcn, datracedcn)
      else
         call get_damat_0d(mol, self%nshell, self%offset, self%hubbard, & 
            & self%hubbard_cn, self%hubbard_exp, self%average, self%gexp, &
            & cache%cn, qsh, dadr, dadL, atrace, dadcn, datracedcn)
      end if
   else
      if (any(mol%periodic)) then
         call get_damat_3d(mol, self%nshell, self%offset, self%hubbard, &
            & self%hubbard_cn, self%hubbard_exp, self%average, self%gexp, &
            & cache%cn, self%rcut, cache%wsc, cache%alpha, qat, &
            & dadr, dadL, atrace, dadcn, datracedcn)
      else
         call get_damat_0d(mol, self%nshell, self%offset, self%hubbard, &
            & self%hubbard_cn, self%hubbard_exp, self%average, self%gexp, &
            & cache%cn, qat, dadr, dadL, atrace, dadcn, datracedcn)
      end if
   end if 

end subroutine get_coulomb_derivs


!> Evaluate uncontracted derivatives of Coulomb matrix for finite system
subroutine get_damat_0d(mol, nshell, offset, hubbard, hubbard_cn, hubbard_exp, &
   & average, gexp, cn, qvec, dadr, dadL, atrace, dadcn, datracedcn)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Number of shells for each atom
   integer, intent(in) :: nshell(:)
   !> Index offset for each shell
   integer, intent(in) :: offset(:)
   !> Hubbard parameter for each shell and species
   real(wp), intent(in) :: hubbard(:, :)
   !> CN-dependence of Hubbard parameters
   real(wp), intent(in) :: hubbard_cn(:)
   !> Exponent of radius-dependent hubbard scaling
   real(wp), intent(in) :: hubbard_exp
   !> Averaging function for Hubbard parameter of a shell-pair
   type(average_type), intent(in) :: average
   !> Exponent of Coulomb kernel
   real(wp), intent(in) :: gexp
   !> Coordination number
   real(wp), intent(in) :: cn(:)
   !> Partial charge vector
   real(wp), intent(in) :: qvec(:)
   !> Derivative of interactions with respect to cartesian displacements
   real(wp), intent(out) :: dadr(:, :, :)
   !> Derivative of interactions with respect to strain deformations
   real(wp), intent(out) :: dadL(:, :, :)
   !> On-site derivatives with respect to cartesian displacements
   real(wp), intent(out) :: atrace(:, :)
   !> Coulomb matrix derivative w.r.t. the coordination numbers
   real(wp), intent(out) :: dadcn(:, :)
   !> On-site derivatives with respect to coordination numbers
   real(wp), intent(out) :: datracedcn(:)

   integer :: iat, jat, izp, jzp, ii, jj, ish, jsh
   real(wp) :: vec(3), r1, r1g, damp, hubi, hubj, dhubidcn, dhubjdcn
   real(wp) :: sqrtcni, sqrtcnj, scalei, scalej, dscaleidcn, dscalejdcn 
   real(wp) :: gam, dgami, dgamj, denom, dtmp, dtmpdgam, dG(3), dS(3, 3)
   real(wp), allocatable :: itrace(:, :), didr(:, :, :), didL(:, :, :)
   real(wp), allocatable :: didcn(:, :), ditracedcn(:)

   atrace(:, :) = 0.0_wp
   dadr(:, :, :) = 0.0_wp
   dadL(:, :, :) = 0.0_wp
   dadcn(:, :) = 0.0_wp
   datracedcn(:) = 0.0_wp

   !$omp parallel default(none) &
   !$omp shared(dadr, dadL, atrace, dadcn, datracedcn, mol, qvec, nshell) &
   !$omp shared(offset, hubbard, hubbard_cn, hubbard_exp, average, gexp, cn) &
   !$omp private(iat, izp, ii, ish, jat, jzp, jj, jsh, r1, r1g, vec, damp) &
   !$omp private(sqrtcni, sqrtcnj, scalei, scalej, dscaleidcn, dscalejdcn) &
   !$omp private(hubi, hubj, dhubidcn, dhubjdcn, gam, dgami, dgamj, denom) &
   !$omp private(dtmp, dtmpdgam, dG, dS, itrace, didr, didL, didcn, ditracedcn)
   allocate(itrace, source=atrace)
   allocate(didr, source=dadr)
   allocate(didL, source=dadL)
   allocate(didcn, source=dadcn)
   allocate(ditracedcn, source=datracedcn)
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      ii = offset(iat)
      ! Coordination number scaling
      sqrtcni = sqrt(cn(iat) + cn_reg**2) - cn_reg
      scalei = 1.0_wp + hubbard_cn(izp) * sqrtcni
      dscaleidcn = 0.5_wp * hubbard_cn(izp) / (sqrtcni + cn_reg)
      do jat = 1, iat-1
         jzp = mol%id(jat)
         jj = offset(jat)
         ! Coordination number scaling
         sqrtcnj = sqrt(cn(jat) + cn_reg**2) - cn_reg
         scalej = 1.0_wp + hubbard_cn(jzp) * sqrtcnj
         dscalejdcn = 0.5_wp * hubbard_cn(jzp) / (sqrtcnj + cn_reg)
         vec = mol%xyz(:, iat) - mol%xyz(:, jat)
         r1 = norm2(vec)
         r1g = r1**gexp
         ! Damping factor for distance-dependent hubbard scaling
         damp = exp(-hubbard_exp * r1)
         do ish = 1, nshell(iat)
            ! CN dependent hubbard parameter and its CN derivative
            hubi = hubbard(ish, izp) * scalei
            dhubidcn = hubbard(ish, izp) * dscaleidcn
            do jsh = 1, nshell(jat)
               ! CN dependent hubbard parameter and its CN derivative
               hubj = hubbard(jsh, jzp) * scalej
               dhubjdcn = hubbard(jsh, jzp) * dscalejdcn
               ! Damped average hubbard parameter and its CN derivative
               gam = average%value(hubi, hubj) / damp
               dgami = (average%deriv(hubi, hubj) * dhubidcn) / damp
               dgamj = (average%deriv(hubj, hubi) * dhubjdcn) / damp
               ! Chain rule derivative of effective Coulomb interaction
               denom = 1.0_wp / (r1g + gam**(-gexp))**(1.0_wp + 1.0_wp/gexp)
               ! Direct position derivative of the effective Coulomb interaction
               dtmp = -(r1**(gexp-1.0_wp) - hubbard_exp * gam**(-gexp)) * denom
               ! Implicit damped averaged hubbard parameter derivative
               dtmpdgam = gam**(-gexp-1.0_wp) * denom
               dG = dtmp*vec/r1
               dS = spread(dG, 1, 3) * spread(vec, 2, 3)
               ! Build once-contracted derivatives ∂γ · q, including 
               ! only other atom contributions as second order is bilinear in q
               itrace(:, ii+ish) = +dG*qvec(jj+jsh) + itrace(:, ii+ish)
               itrace(:, jj+jsh) = -dG*qvec(ii+ish) + itrace(:, jj+jsh)
               didr(:, iat, jj+jsh) = +dG*qvec(ii+ish) + didr(:, iat, jj+jsh)
               didr(:, jat, ii+ish) = -dG*qvec(jj+jsh) + didr(:, jat, ii+ish)
               didL(:, :, jj+jsh) = +dS*qvec(ii+ish) + didL(:, :, jj+jsh)
               didL(:, :, ii+ish) = +dS*qvec(jj+jsh) + didL(:, :, ii+ish)
               didcn(jat, ii+ish) = +dtmpdgam*dgamj*qvec(jj+jsh) + didcn(jat, ii+ish)
               didcn(iat, jj+jsh) = +dtmpdgam*dgami*qvec(ii+ish) + didcn(iat, jj+jsh)
               ditracedcn(ii+ish) = +dtmpdgam*dgami*qvec(jj+jsh) + ditracedcn(ii+ish)
               ditracedcn(jj+jsh) = +dtmpdgam*dgamj*qvec(ii+ish) + ditracedcn(jj+jsh)
            end do
         end do
      end do
      do ish = 1, nshell(iat)
         ! CN dependent hubbard parameter and its CN derivative
         hubi = hubbard(ish, izp) * scalei
         dhubidcn = hubbard(ish, izp) * dscaleidcn
         do jsh = 1, ish-1
            ! CN dependent hubbard parameter and its CN derivative
            hubj = hubbard(jsh, izp) * scalei
            dhubjdcn = hubbard(jsh, izp) * dscaleidcn
            ! CN derivative of diagonal averaged hubbard parameter
            dgami = average%deriv(hubi, hubj) * dhubidcn
            dgamj = average%deriv(hubj, hubi) * dhubjdcn
            didcn(iat, ii+ish) = +dgami*qvec(ii+jsh) + didcn(iat, ii+ish)
            didcn(iat, ii+jsh) = +dgamj*qvec(ii+ish) + didcn(iat, ii+jsh)
            ditracedcn(ii+ish) = -dgami*qvec(ii+jsh) + ditracedcn(ii+ish)
            ditracedcn(ii+jsh) = -dgamj*qvec(ii+ish) + ditracedcn(ii+jsh)
         end do
         didcn(iat, ii+ish) = +0.5_wp*dhubidcn*qvec(ii+ish) + didcn(iat, ii+ish)
         ditracedcn(ii+ish) = -0.5_wp*dhubidcn*qvec(ii+ish) + ditracedcn(ii+ish)
      end do
   end do
   !$omp critical (get_damat_0d_)
   atrace(:, :) = atrace + itrace
   dadr(:, :, :) = dadr + didr
   dadL(:, :, :) = dadL + didL
   dadcn(:, :) = dadcn + didcn
   datracedcn(:) = datracedcn + ditracedcn
   !$omp end critical (get_damat_0d_)
   deallocate(didL, didr, itrace, didcn, ditracedcn)
   !$omp end parallel

end subroutine get_damat_0d

!> Evaluate uncontracted derivatives of Coulomb matrix for 3D periodic system
subroutine get_damat_3d(mol, nshell, offset, hubbard, hubbard_cn, hubbard_exp, &
      & average, gexp, cn, rcut, wsc, alpha, qvec, dadr, dadL, atrace, dadcn, &
      & datracedcn)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Number of shells for each atom
   integer, intent(in) :: nshell(:)
   !> Index offset for each shell
   integer, intent(in) :: offset(:)
   !> Hubbard parameter for each shell and species
   real(wp), intent(in) :: hubbard(:, :)
   !> CN-dependence of Hubbard parameters
   real(wp), intent(in) :: hubbard_cn(:)
   !> Exponent of radius-dependent hubbard scaling
   real(wp), intent(in) :: hubbard_exp
   !> Averaging function for Hubbard parameter of a shell-pair
   type(average_type), intent(in) :: average
   !> Exponent of Coulomb kernel
   real(wp), intent(in) :: gexp
   !> Coordination number
   real(wp), intent(in) :: cn(:)
   !> Long-range cutoff
   real(wp), intent(in) :: rcut
   !> Wigner-Seitz image information
   type(wignerseitz_cell), intent(in) :: wsc
   !> Convergence factor for Ewald sum
   real(wp), intent(in) :: alpha
   !> Partial charge vector
   real(wp), intent(in) :: qvec(:)
   !> Derivative of interactions with respect to cartesian displacements
   real(wp), intent(out) :: dadr(:, :, :)
   !> Derivative of interactions with respect to strain deformations
   real(wp), intent(out) :: dadL(:, :, :)
   !> On-site derivatives with respect to cartesian displacements
   real(wp), intent(out) :: atrace(:, :)
   !> Coulomb matrix derivative w.r.t. the coordination numbers
   real(wp), intent(out) :: dadcn(:, :)
   !> On-site derivatives with respect to coordination numbers
   real(wp), intent(out) :: datracedcn(:)

   integer :: iat, jat, izp, jzp, img, ii, jj, ish, jsh
   real(wp) :: sqrtcni, sqrtcnj, scalei, scalej
   real(wp) :: vol, gam, hubi, hubj, wsw, vec(3), dG(3), dS(3, 3)
   real(wp) :: dGd(3), dSd(3, 3), dGr(3), dSr(3, 3)
   real(wp), allocatable :: itrace(:, :), didr(:, :, :), didL(:, :, :)
   real(wp), allocatable :: didcn(:, :), ditracedcn(:)
   real(wp), allocatable :: dtrans(:, :), rtrans(:, :)

   atrace(:, :) = 0.0_wp
   dadr(:, :, :) = 0.0_wp
   dadL(:, :, :) = 0.0_wp
   dadcn(:, :) = 0.0_wp
   datracedcn(:) = 0.0_wp

   vol = abs(matdet_3x3(mol%lattice))
   call get_dir_trans(mol%lattice, alpha, conv, dtrans)
   call get_rec_trans(mol%lattice, alpha, vol, conv, rtrans)

   !$omp parallel default(none) shared(atrace, dadr, dadL, dadcn, datracedcn) &
   !$omp shared(mol, wsc, alpha, vol, dtrans, rtrans, qvec, nshell, offset) &
   !$omp shared(hubbard, hubbard_cn, hubbard_exp, cn, average, gexp, rcut) &
   !$omp private(iat, izp, jat, jzp, img, ii, jj, ish, jsh, wsw, vec, gam) &
   !$omp private(scalei, scalej, sqrtcni, sqrtcnj, hubi, hubj, dG, dS) &
   !$omp private(dGr, dSr, dGd, dSd, itrace, didr, didL, didcn, ditracedcn)
   allocate(itrace, source=atrace)
   allocate(didr, source=dadr)
   allocate(didL, source=dadL)
   allocate(didcn, source=dadcn)
   allocate(ditracedcn, source=datracedcn)
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      ii = offset(iat)
      ! Coordination number scaling
      sqrtcni = sqrt(cn(iat) + cn_reg**2) - cn_reg
      scalei = 1.0_wp + hubbard_cn(izp) * sqrtcni
      do jat = 1, iat-1
         jzp = mol%id(jat)
         jj = offset(jat)
         ! Coordination number scaling
         sqrtcnj = sqrt(cn(jat) + cn_reg**2) - cn_reg
         scalej = 1.0_wp + hubbard_cn(jzp) * sqrtcnj
         do img = 1, wsc%nimg(jat, iat)
            vec = mol%xyz(:, iat) - mol%xyz(:, jat) - wsc%trans(:, wsc%tridx(img, jat, iat))
            wsw = 1.0_wp / real(wsc%nimg(jat, iat), wp)
            call get_damat_rec_3d(vec, vol, alpha, 0.0_wp, rtrans, dGr, dSr)
            do ish = 1, nshell(iat)
               hubi = hubbard(ish, izp) * scalei
               do jsh = 1, nshell(jat)
                  hubj = hubbard(jsh, jzp) * scalej
                  gam = average%value(hubi, hubj)
                  call get_damat_dir_3d(vec, gam, gexp, rcut, alpha, dtrans, dGd, dSd)
                  dG = (dGd + dGr) * wsw
                  dS = (dSd + dSr) * wsw
                  itrace(:, ii+ish) = +dG*qvec(jj+jsh) + itrace(:, ii+ish)
                  itrace(:, jj+jsh) = -dG*qvec(ii+ish) + itrace(:, jj+jsh)
                  didr(:, iat, jj+jsh) = +dG*qvec(ii+ish) + didr(:, iat, jj+jsh)
                  didr(:, jat, ii+ish) = -dG*qvec(jj+jsh) + didr(:, jat, ii+ish)
                  didL(:, :, jj+jsh) = +dS*qvec(ii+ish) + didL(:, :, jj+jsh)
                  didL(:, :, ii+ish) = +dS*qvec(jj+jsh) + didL(:, :, ii+ish)
               end do
            end do
         end do
      end do

      wsw = 1.0_wp / real(wsc%nimg(iat, iat), wp)
      do img = 1, wsc%nimg(iat, iat)
         vec = wsc%trans(:, wsc%tridx(img, iat, iat))
         call get_damat_rec_3d(vec, vol, alpha, 0.0_wp, rtrans, dGr, dSr)
         do ish = 1, nshell(iat)
            hubi = hubbard(ish, izp) * scalei
            do jsh = 1, ish-1
               hubj = hubbard(jsh, izp) * scalei
               gam = average%value(hubi, hubj)
               call get_damat_dir_3d(vec, gam, gexp, rcut, alpha, dtrans, dGd, dSd)
               dS = (dSd + dSr) * wsw
               didL(:, :, ii+jsh) = +dS*qvec(ii+ish) + didL(:, :, ii+jsh)
               didL(:, :, ii+ish) = +dS*qvec(ii+jsh) + didL(:, :, ii+ish)
            end do
            gam = hubi
            call get_damat_dir_3d(vec, gam, gexp, rcut, alpha, dtrans, dGd, dSd)
            dS = (dSd + dSr) * wsw
            didL(:, :, ii+ish) = +dS*qvec(ii+ish) + didL(:, :, ii+ish)
         end do
      end do
   end do
   !$omp critical (get_damat_3d_)
   atrace(:, :) = atrace + itrace
   dadr(:, :, :) = dadr + didr
   dadL(:, :, :) = dadL + didL
   dadcn(:, :) = dadcn + didcn
   datracedcn(:) = datracedcn + ditracedcn
   !$omp end critical (get_damat_3d_)
   deallocate(didL, didr, itrace, didcn, ditracedcn) 
   !$omp end parallel

end subroutine get_damat_3d

!> Calculate real space contributions for a pair under 3D periodic boundary conditions
subroutine get_damat_dir_3d(rij, gam, gexp, rcut, alp, trans, dg, ds)
   !> Distance between pair
   real(wp), intent(in) :: rij(3)
   !> Hubbard parameter
   real(wp), intent(in) :: gam
   !> Exponent for interaction kernel
   real(wp), intent(in) :: gexp
   !> Long-range cutoff
   real(wp), intent(in) :: rcut
   !> Convergence factor
   real(wp), intent(in) :: alp
   !> Translation vectors to consider
   real(wp), intent(in) :: trans(:, :)
   !> Derivative with respect to cartesian displacements
   real(wp), intent(out) :: dg(3)
   !> Derivative with respect to strain deformations
   real(wp), intent(out) :: ds(3, 3)

   integer :: itr
   real(wp) :: vec(3), r1, r2, gtmp, atmp, alp2, fcut, dcut

   dg(:) = 0.0_wp
   ds(:, :) = 0.0_wp

   alp2 = alp*alp

   do itr = 1, size(trans, 2)
      vec(:) = rij + trans(:, itr)
      r1 = norm2(vec)
      if (r1 < eps) cycle
      r2 = r1*r1
      fcut = fsmooth(r1, rcut)
      dcut = dsmooth(r1, rcut)
      gtmp = 1.0_wp / (r1**gexp + gam**(-gexp))
      gtmp = -r1**(gexp-2.0_wp) * gtmp * gtmp**(1.0_wp/gexp) * fcut - (1.0_wp-fcut)/(r2*r1) &
         & + dcut * (gtmp**(1.0_wp/gexp) * r1 - 1.0_wp) / r2
      atmp = -2*alp*exp(-r2*alp2)/(sqrtpi*r2) + erf(r1*alp)/(r2*r1)
      dg(:) = dg + (gtmp + atmp) * vec
      ds(:, :) = ds + (gtmp + atmp) * spread(vec, 1, 3) * spread(vec, 2, 3)
   end do

end subroutine get_damat_dir_3d

!> Calculate reciprocal space contributions for a pair under 3D periodic boundary conditions
subroutine get_damat_rec_3d(rij, vol, alp, qpc, trans, dg, ds)
   !> Distance between pair
   real(wp), intent(in) :: rij(3)
   !> Cell volume
   real(wp), intent(in) :: vol
   !> Convergence factor
   real(wp), intent(in) :: alp
   !> Pseudo-quadrupole contribution
   real(wp), intent(in) :: qpc
   !> Translation vectors to consider
   real(wp), intent(in) :: trans(:, :)
   !> Derivative with respect to cartesian displacements
   real(wp), intent(out) :: dg(3)
   !> Derivative with respect to strain deformations
   real(wp), intent(out) :: ds(3, 3)

   integer :: itr
   real(wp) :: fac, vec(3), g2, gv, expk, sink, cosk, alp2, fqp
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))

   dg(:) = 0.0_wp
   ds(:, :) = 0.0_wp
   fac = 4*pi/vol
   alp2 = alp*alp
   fqp = 2*qpc*qpc

   do itr = 1, size(trans, 2)
      vec(:) = trans(:, itr)
      g2 = dot_product(vec, vec)
      if (g2 < eps) cycle
      gv = dot_product(rij, vec)
      expk = fac * exp(-0.25_wp*g2/alp2)/g2
      cosk = cos(gv) * expk
      sink = sin(gv) * expk
      dg(:) = dg - sink * vec
      ds(:, :) = ds + cosk &
         & * ((2.0_wp/g2 + 0.5_wp/alp2 + 0.5_wp*fqp) * spread(vec, 1, 3)*spread(vec, 2, 3) &
         &     - unity * (1.0_wp + g2*fqp))
   end do

end subroutine get_damat_rec_3d


end module tblite_coulomb_charge_effective
