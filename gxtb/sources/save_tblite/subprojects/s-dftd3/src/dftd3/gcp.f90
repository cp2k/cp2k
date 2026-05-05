! This file is part of s-dftd3.
! SPDX-Identifier: LGPL-3.0-or-later
!
! s-dftd3 is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! s-dftd3 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with s-dftd3.  If not, see <https://www.gnu.org/licenses/>.

module dftd3_gcp
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use dftd3_gcp_param, only : gcp_param, get_gcp_param
   use dftd3_cutoff, only : realspace_cutoff, get_lattice_points
   implicit none
   private

   public :: gcp_param, get_gcp_param, get_geometric_counterpoise

   !> Geometric counterpoise correction
   interface get_geometric_counterpoise
      module procedure get_geometric_counterpoise
      module procedure get_geometric_counterpoise_atomic
   end interface get_geometric_counterpoise



contains


!> Geometric counterpoise correction
subroutine get_geometric_counterpoise(mol, param, cutoff, energy, gradient, sigma)

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Geometric counterpoise parameters
   type(gcp_param), intent(in) :: param

   !> Realspace cutoffs
   type(realspace_cutoff), intent(in) :: cutoff

   !> Counter-poise energy
   real(wp), intent(inout) :: energy

   !> Counter-poise gradient
   real(wp), intent(inout), optional :: gradient(:, :)

   !> Counter-poise virial
   real(wp), intent(inout), optional :: sigma(:, :)

   real(wp), allocatable :: energies(:)

   allocate(energies(mol%nat), source=0.0_wp)
   call get_geometric_counterpoise_atomic(mol, param, cutoff, energies, gradient, sigma)
   energy = energy + sum(energies)
end subroutine get_geometric_counterpoise


!> Geometric counterpoise correction with atom-resolved energies
subroutine get_geometric_counterpoise_atomic(mol, param, cutoff, energies, gradient, sigma)

   !> Molecular structure data
   class(structure_type), intent(in) :: mol

   !> Geometric counterpoise parameters
   type(gcp_param), intent(in) :: param

   !> Realspace cutoffs
   type(realspace_cutoff), intent(in) :: cutoff

   !> Dispersion energy
   real(wp), intent(inout) :: energies(:)

   !> Dispersion gradient
   real(wp), intent(inout), optional :: gradient(:, :)

   !> Dispersion virial
   real(wp), intent(inout), optional :: sigma(:, :)

   real(wp), allocatable :: lattr(:, :)
   real(wp), parameter :: rexp_base = 0.75_wp, zexp_base = 1.5_wp, rexp_srb = -1.0_wp, zexp_srb = 0.5_wp
   logical :: grad

   grad = present(gradient) .and. present(sigma)

   if (allocated(param%emiss) .and. allocated(param%slater) .and. allocated(param%xv)) then
      call get_lattice_points(mol%periodic, mol%lattice, cutoff%gcp, lattr)
      if (grad) then
         call gcp_deriv(mol, lattr, cutoff%gcp, param%zeff, param%emiss, param%slater, &
            & param%xv, param%rvdw, param%sigma, param%alpha, param%beta, param%damp, &
            & param%dmp_scal, param%dmp_exp, energies, gradient, sigma)
      else
         call gcp_energy(mol, lattr, cutoff%gcp, param%zeff, param%emiss, param%slater, &
            & param%xv, param%rvdw, param%sigma, param%alpha, param%beta, param%damp, &
            & param%dmp_scal, param%dmp_exp, energies)
      end if
   end if

   if (param%srb .or. param%base) then
      call get_lattice_points(mol%periodic, mol%lattice, cutoff%srb, lattr)
   end if

   if (param%srb) then
      if (grad) then
         call srb_deriv(mol, lattr, cutoff%srb, mol%num, param%rvdw_srb, param%rscal, param%qscal, &
            & rexp_srb, zexp_srb, energies, gradient, sigma)
      else
         call srb_energy(mol, lattr, cutoff%srb, mol%num, param%rvdw_srb, param%rscal, param%qscal, &
            & rexp_srb, zexp_srb, energies)
      end if
   end if

   if (param%base) then
      if (grad) then
         call srb_deriv(mol, lattr, cutoff%srb, param%zeff, param%rvdw, param%rscal, param%qscal, &
            & rexp_base, zexp_base, energies, gradient, sigma)
      else
         call srb_energy(mol, lattr, cutoff%srb, param%zeff, param%rvdw, param%rscal, param%qscal, &
            & rexp_base, zexp_base, energies)
      end if
   end if
end subroutine get_geometric_counterpoise_atomic


!> Geometric counterpoise correction
subroutine gcp_energy(mol, trans, cutoff, iz, emiss, slater, xv, rvdw, escal, alpha, beta, &
      & damp, dmp_scal, dmp_exp, energies)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Translation vectors
   real(wp), intent(in) :: trans(:, :)
   !> Distance cutoff
   real(wp), intent(in) :: cutoff
   !> Effective nuclear charges
   integer, intent(in) :: iz(:)
   !> Basis set superposition error per atom
   real(wp), intent(in) :: emiss(:)
   !> Slater exponents
   real(wp), intent(in) :: slater(:)
   !> Number of virtual orbitals
   real(wp), intent(in) :: xv(:)
   !> Van der Waals radii
   real(wp), intent(in) :: rvdw(:, :)
   !> Scaling factor
   real(wp), intent(in) :: escal
   !> Exponential factor
   real(wp), intent(in) :: alpha
   !> Power factor
   real(wp), intent(in) :: beta
   !> Damping flag
   logical, intent(in) :: damp
   !> Damping scaling factor
   real(wp), intent(in) :: dmp_scal
   !> Damping exponent
   real(wp), intent(in) :: dmp_exp
   !> Atom-resolved energy
   real(wp), intent(inout) :: energies(:)

   integer :: iat, jat, jtr, izp, jzp
   real(wp) :: argv
   real(wp) :: xvi, xvj
   real(wp) :: r1, rscal, rscalexp, r0
   real(wp) :: sij, expv, bsse, emi, emj
   real(wp) :: vec(3)
   real(wp) :: dampval, grd_dmp
   real(wp) :: dE

   !$omp parallel do default(none) &
   !$omp reduction(+:energies) &
   !$omp shared(mol, iz, xv, emiss, rvdw, trans, cutoff, alpha, beta, &
   !$omp&       dmp_scal, dmp_exp, escal, slater, damp) &
   !$omp private(izp, jat, jzp, xvi, xvj, emi, emj, r0, jtr, vec, r1, sij, &
   !$omp&        expv, bsse, argv, dE, dampval, grd_dmp, rscal, rscalexp)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      xvi = merge(1.0_wp / sqrt(xv(izp)), 0.0_wp, xv(izp) >= 0.5_wp)

      ! the BSSE due to atom jat, Loop over all j atoms
      do jat = 1, iat
         jzp = mol%id(jat)
         xvj = merge(1.0_wp / sqrt(xv(jzp)), 0.0_wp, xv(jzp) >= 0.5_wp)

         emi = emiss(izp)*xvj*escal
         emj = emiss(jzp)*xvi*escal
         r0 = rvdw(izp, jzp)
         do jtr = 1, size(trans, 2)
            vec(:) = mol%xyz(:, iat) - (mol%xyz(:, jat) + trans(:, jtr))
            r1 = norm2(vec)

            if(r1 > cutoff .or. r1 < epsilon(1.0_wp)) cycle
            ! calulate slater overlap sij
            call ssovl(r1, izp, jzp, iz, slater(izp), slater(jzp), sij)
            ! evaluate gcp central expression
            argv = -alpha*r1**beta
            expv = exp(argv)
            bsse = expv/sqrt(sij)
            if (damp) then
               rscal = r1/r0
               rscalexp = dmp_scal*rscal**dmp_exp
               dampval = (1.0_wp-1.0_wp/(1.0_wp+rscalexp))
            else
               dampval = 1.0_wp
            endif
            dE = bsse*dampval
            energies(iat) = energies(iat) + emi*dE
            if (iat /= jat) then
               energies(jat) = energies(jat) + emj*dE
            end if
         end do
      end do
   end do

end subroutine gcp_energy


!> Geometric counterpoise correction
subroutine gcp_deriv(mol, trans, cutoff, iz, emiss, slater, xv, rvdw, escal, alpha, beta, &
      & damp, dmp_scal, dmp_exp, energies, gradient, sigma)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Translation vectors
   real(wp), intent(in) :: trans(:, :)
   !> Distance cutoff
   real(wp), intent(in) :: cutoff
   !> Effective nuclear charges
   integer, intent(in) :: iz(:)
   !> Basis set superposition error per atom
   real(wp), intent(in) :: emiss(:)
   !> Slater exponents
   real(wp), intent(in) :: slater(:)
   !> Number of virtual orbitals
   real(wp), intent(in) :: xv(:)
   !> Van der Waals radii
   real(wp), intent(in) :: rvdw(:, :)
   !> Scaling factor
   real(wp), intent(in) :: escal
   !> Exponential factor
   real(wp), intent(in) :: alpha
   !> Power factor
   real(wp), intent(in) :: beta
   !> Damping flag
   logical, intent(in) :: damp
   !> Damping scaling factor
   real(wp), intent(in) :: dmp_scal
   !> Damping exponent
   real(wp), intent(in) :: dmp_exp
   !> Atom-resolved energy
   real(wp), intent(inout) :: energies(:)
   !> Molecular gradient
   real(wp), intent(inout) :: gradient(:, :)
   !> Virial
   real(wp), intent(inout) :: sigma(:, :)

   integer :: iat, jat, jtr, izp, jzp
   real(wp) :: argv, argd, ovlpd
   real(wp) :: xvi, xvj
   real(wp) :: r1, rscal, rscalexp, rscalexpm1, r0
   real(wp) :: sij, expv, expd, bsse, emi, emj, emij
   real(wp) :: vec(3), gs(3), gij
   real(wp) :: dampval, grd_dmp
   real(wp) :: dE, dG(3), dS(3, 3)

   !$omp parallel do default(none) &
   !$omp reduction(+:energies, gradient, sigma) &
   !$omp shared(mol, iz, xv, emiss, rvdw, trans, cutoff, alpha, beta, &
   !$omp&       dmp_scal, dmp_exp, escal, slater, damp) &
   !$omp private(izp, jat, jzp, xvi, xvj, emi, emj, r0, jtr, vec, r1, sij, gij, emij, &
   !$omp&        expv, bsse, argv, dE, dampval, grd_dmp, dG, dS, expd, argd, ovlpd, gs, &
   !$omp&        rscal, rscalexp, rscalexpm1)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      xvi = merge(1.0_wp / sqrt(xv(izp)), 0.0_wp, xv(izp) >= 0.5_wp)

      ! the BSSE due to atom jat, Loop over all j atoms
      do jat = 1, iat
         jzp = mol%id(jat)
         xvj = merge(1.0_wp / sqrt(xv(jzp)), 0.0_wp, xv(jzp) >= 0.5_wp)

         emi = emiss(izp)*xvj*escal
         emj = emiss(jzp)*xvi*escal
         emij = emi + emj
         r0 = rvdw(izp, jzp)
         do jtr = 1, size(trans, 2)
            vec(:) = mol%xyz(:, iat) - (mol%xyz(:, jat) + trans(:, jtr))
            r1 = norm2(vec)

            if(r1 > cutoff .or. r1 < epsilon(1.0_wp)) cycle
            ! calulate slater overlap sij
            call ssovl(r1, izp, jzp, iz, slater(izp), slater(jzp), sij)
            ! evaluate gcp central expression
            argv = -alpha*r1**beta
            expv = exp(argv)
            bsse = expv/sqrt(sij)
            if (damp) then
               rscal = r1/r0
               rscalexp = dmp_scal*rscal**dmp_exp
               dampval = (1.0_wp-1.0_wp/(1.0_wp+rscalexp))
            else
               dampval = 1.0_wp
            endif

            call gsovl(r1, izp, jzp, iz, slater(izp), slater(jzp), gij)

            gs(:) = gij*vec
            expd = exp(-alpha*r1**beta)*(-0.5_wp)
            argd = 2d0*alpha*beta*r1**beta*sij/r1
            ovlpd = r1*sij**1.5_wp

            if(damp) then
               rscalexpm1 = rscal**(dmp_exp-1)
               grd_dmp = dmp_scal*dmp_exp*rscalexpm1/r0
               grd_dmp = grd_dmp/(rscalexp+1.0_wp)**2
            endif

            dE = bsse*dampval
            dG = expd*(argd*vec + gs)/ovlpd*emij
            if(damp) then
               dG = dG*dampval+bsse*grd_dmp*(vec/r1)*emij
            endif
            dS = spread(dG, 1, 3) * spread(vec, 2, 3) * 0.5_wp

            energies(iat) = energies(iat) + emi*dE
            if (iat /= jat) then
               energies(jat) = energies(jat) + emj*dE
            end if
            sigma(:, :) = sigma(:, :) + dS
            if (iat /= jat) then
               gradient(:, iat) = gradient(:, iat) + dG
               gradient(:, jat) = gradient(:, jat) - dG
               sigma(:, :) = sigma(:, :) + dS
            end if
         end do
      end do
   end do

end subroutine gcp_deriv


!> Short-range bond length correction for HF-3c
subroutine srb_energy(mol, trans, cutoff, iz, r0ab, rscal, qscal, rexp, zexp, energies)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Translation vectors
   real(wp), intent(in) :: trans(:, :)
   !> Distance cutoff
   real(wp), intent(in) :: cutoff
   !> Effective nuclear charges
   integer, intent(in) :: iz(:)
   !> Van der Waals radii
   real(wp), intent(in) :: r0ab(:, :)
   !> Radii scaling factor
   real(wp), intent(in) :: rscal
   !> Prefactor for the SRB potential
   real(wp), intent(in) :: qscal
   !> Exponent for radii
   real(wp), intent(in) :: rexp
   !> Exponent for charges
   real(wp), intent(in) :: zexp
   !> Atom-resolved energy
   real(wp), intent(inout) :: energies(:)

   real(wp) :: fi, fj, ff, r1, expt
   real(wp) :: r0, vec(3), dE
   integer iat, jat, jtr, izp, jzp

   !$omp parallel do default(none) &
   !$omp reduction(+:energies) &
   !$omp shared(mol, trans, cutoff, iz, r0ab, rexp, zexp, rscal, qscal) &
   !$omp private(izp, jat, jzp, r0, vec, r1, fi, fj, ff, expt, dE)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      fi = real(iz(izp), wp)
      do jat = 1, iat
         jzp = mol%id(jat)
         r0 = rscal*r0ab(izp, jzp)**rexp
         fj = real(iz(jzp), wp)
         ff = -(fi*fj)**zexp

         do jtr = 1, size(trans, 2)
            vec(:) = mol%xyz(:, iat) - (mol%xyz(:, jat) + trans(:, jtr))
            r1 = norm2(vec)
            if(r1 > cutoff .or. r1 < epsilon(1.0_wp)) cycle
            expt = exp(-r0*r1)
            dE = qscal*ff*expt*0.5_wp

            energies(iat) = energies(iat) + dE
            if (iat /= jat) then
               energies(jat) = energies(jat) + dE
            end if
         end do
      end do
   end do

end subroutine srb_energy


!> Short-range bond length correction for HF-3c
subroutine srb_deriv(mol, trans, cutoff, iz, r0ab, rscal, qscal, rexp, zexp, energies, gradient, sigma)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Translation vectors
   real(wp), intent(in) :: trans(:, :)
   !> Distance cutoff
   real(wp), intent(in) :: cutoff
   !> Effective nuclear charges
   integer, intent(in) :: iz(:)
   !> Van der Waals radii
   real(wp), intent(in) :: r0ab(:, :)
   !> Radii scaling factor
   real(wp), intent(in) :: rscal
   !> Prefactor for the SRB potential
   real(wp), intent(in) :: qscal
   !> Exponent for radii
   real(wp), intent(in) :: rexp
   !> Exponent for charges
   real(wp), intent(in) :: zexp
   !> Atom-resolved energy
   real(wp), intent(inout) :: energies(:)
   !> Molecular gradient
   real(wp), intent(inout) :: gradient(:, :)
   !> Molecular virial
   real(wp), intent(inout) :: sigma(:, :)

   real(wp) :: fi, fj, ff, rf, r1, expt
   real(wp) :: r0, vec(3), dE, dG(3), dS(3, 3)
   integer iat, jat, jtr, izp, jzp

   !$omp parallel do default(none) &
   !$omp reduction(+:energies, gradient, sigma) &
   !$omp shared(mol, trans, cutoff, iz, r0ab, rexp, zexp, rscal, qscal) &
   !$omp private(izp, jat, jzp, r0, vec, r1, fi, fj, ff, expt, rf, dE, dG, dS)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      fi = real(iz(izp), wp)
      do jat = 1, iat
         jzp = mol%id(jat)
         r0 = rscal*r0ab(izp, jzp)**rexp
         fj = real(iz(jzp), wp)
         ff = -(fi*fj)**zexp

         do jtr = 1, size(trans, 2)
            vec(:) = mol%xyz(:, iat) - (mol%xyz(:, jat) + trans(:, jtr))
            r1 = norm2(vec)
            if(r1 > cutoff .or. r1 < epsilon(1.0_wp)) cycle
            expt = exp(-r0*r1)
            rf = qscal/r1

            dE = qscal*ff*expt*0.5_wp
            dG(:) = -ff*r0*vec*expt*rf
            dS(:, :) = spread(dG, 1, 3) * spread(vec, 2, 3) * 0.5_wp

            energies(iat) = energies(iat) + dE
            if (iat /= jat) then
               energies(jat) = energies(jat) + dE
            end if
            sigma(:, :) = sigma + dS
            if (iat /= jat) then
               gradient(:, iat) = gradient(:, iat) + dG
               gradient(:, jat) = gradient(:, jat) - dG
               sigma(:, :) = sigma + dS
            end if

         end do
      end do
   end do

end subroutine srb_deriv


!******************************************************************************
!* calculates the s-type overlap integral over 1s, 2s and 3s slater functions
!* added support for 3s functions
!* ovl = overlap integral
!* za  = slater exponent atom A
!* zb  = slater exponent atom B
!* R   = distance between atom A and B
!* Inspired by mopac7.0
!******************************************************************************
subroutine ssovl(r, iat, jat, iz, xza, xzb, ovl)
   integer ii, shell(72)
   logical debug
   real(wp) za, zb, R, ovl, ax, bx, norm, R05
   integer na, nb
   real(wp) xx
   data shell/                 &
   !          h, he
            1, 1               &
   !         li-ne
            , 2, 2, 2, 2, 2, 2, 2, 2, &
   !         na-ar
            3, 3, 3, 3, 3, 3, 3, 3,  &
   ! 4s, 5s will be treated as 3s
   !         k-rn , no f-elements
            54*3/
   ! ...
   real(wp) xza, xzb
   integer iat, jat, iz(*)

   za = xza
   zb = xzb
   na = iz(iat)
   nb = iz(jat)
   debug = .false.
!debug = .true.

! ii selects kind of ovl by multiplying the shell
! kind    <1s|1s>  <2s|1s>  <2s|2s>  <1s|3s>  <2s|3s>  <3s|3s>
! case:      1        2        4       3        6         9
!
   ii = shell(na)*shell(nb)
   if(debug) write(*, *) 'shell', ii

   R05 = R*0.5
   ax = (za+zb)*R05
   bx = (zb-za)*R05

   ! same elements
   if(abs(za-zb) < 0.1) then
      select case (ii)
      case (1)
         ovl = 0.25d0*sqrt((za*zb*R*R)**3)*(A2(ax)*Bint(bx, 0)-Bint(bx, 2)*A0(ax))
      case (2)
         ovl = SQRT(1._wp/3._wp)
         if(shell(na) < shell(nb)) then
         ! <1s|2s>
            norm = SQRT((ZA**3)*(ZB**5))*(R**4)*0.125_wp
            ovl = ovl*norm*(A3(ax)*Bint(bx, 0)-Bint(bx, 3)*A0(ax)+A2(ax)*Bint(bx, 1)-Bint(bx, 2)*A1(ax))
         else
         ! switch za/zb to get <2s|1s>
            xx = za
            za = zb
            zb = xx
            ax = (za+zb)*R05
            bx = (zb-za)*R05
            norm = SQRT((ZA**3)*(ZB**5))*(R**4)*0.125_wp
            ovl = ovl*norm*(A3(ax)*Bint(bx, 0)-Bint(bx, 3)*A0(ax)+A2(ax)*Bint(bx, 1)-Bint(bx, 2)*A1(ax))
         endif
      case (4)
         norm = SQRT((ZA*ZB)**5)*(R**5)*0.0625d0
         ovl = norm* (A4(ax)*Bint(bx, 0)+Bint(bx, 4)*A0(ax)-2.0d0*A2(ax)*Bint(bx, 2))*(1d0/3d0)
      case(3)
         if(shell(na) < shell(nb)) then
            norm = SQRT((ZA**3)*(ZB**7)/7.5_wp)*(R**5)*0.0625_wp
            ovl = norm*(A4(ax)*Bint(bx, 0)-Bint(bx, 4)*A0(ax)+2.d0*(A3(ax)*Bint(bx, 1)-Bint(bx, 3)*A1(ax)))/sqrt(3.d0)
         else
            xx = za
            za = zb
            zb = xx
            ax = (za+zb)*R05
            bx = (zb-za)*R05
            norm = SQRT((ZA**3)*(ZB**7)/7.5_wp)*(R**5)*0.0625_wp
            ovl = norm*(A4(ax)*Bint(bx, 0)-Bint(bx, 4)*A0(ax)+2.d0*(A3(ax)*Bint(bx, 1)-Bint(bx, 3)*A1(ax)))/sqrt(3.d0)
         endif
      case(6)
         if(shell(na) < shell(nb)) then
            norm = SQRT((za**5)*(zb**7)/7.5_wp)*(R**6)*0.03125_wp
            ovl = norm*(A5(ax)*Bint(bx, 0)+A4(ax)*Bint(bx, 1) &
               & -2d0*(A3(ax)*Bint(bx, 2)+A2(ax)*Bint(bx, 3)) &
               & +A1(ax)*Bint(bx, 4)+A0(ax)*Bint(bx, 5))/3.d0
         else
            xx = za
            za = zb
            zb = xx
            ax = (za+zb)*R05
            bx = (zb-za)*R05
            norm = SQRT((za**5)*(zb**7)/7.5_wp)*(R**6)*0.03125_wp
            ovl = norm*(A5(ax)*Bint(bx, 0)+A4(ax)*Bint(bx, 1) &
               & -2d0*(A3(ax)*Bint(bx, 2)+A2(ax)*Bint(bx, 3)) &
               & +A1(ax)*Bint(bx, 4)+A0(ax)*Bint(bx, 5))/3.d0
         endif
      case(9)
         norm = sqrt((ZA*ZB*R*R)**7)/480.d0
         ovl = norm*(A6(ax)*Bint(bx, 0)-3.d0*(A4(ax)*Bint(bx, 2) &
            & -A2(ax)*Bint(bx, 4))-A0(ax)*Bint(bx, 6))/3._wp
      end select
   else ! different elements
      select case (ii)
      case (1)
         norm = 0.25d0*sqrt((za*zb*R*R)**3)
         ovl = (A2(ax)*B0(bx)-B2(bx)*A0(ax))*norm
      case (2)
         ovl = SQRT(1._wp/3._wp)
         if(shell(na) < shell(nb)) then
         ! <1s|2s>
            norm = SQRT((ZA**3)*(ZB**5))*(R**4)*0.125_wp
            ovl = ovl*norm*(A3(ax)*B0(bx)-B3(bx)*A0(ax)+A2(ax)*B1(bx)-B2(bx)*A1(ax))
         else
         ! switch za/zb to get <2s|1s>
            xx = za
            za = zb
            zb = xx
            ax = (za+zb)*R05
            bx = (zb-za)*R05
            norm = SQRT((ZA**3)*(ZB**5))*(R**4)*0.125_wp
            ovl = ovl*norm*(A3(ax)*B0(bx)-B3(bx)*A0(ax)+A2(ax)*B1(bx)-B2(bx)*A1(ax))
         endif
      case (4) ! <2s|2s>
         norm = SQRT((ZA*ZB)**5)*(R**5)*0.0625_wp
         ovl = norm* (A4(ax)*B0(bx)+B4(bx)*A0(ax)-2.0_wp*A2(ax)*B2(bx))*(1d0/3d0)
      case(3)  ! <1s|3s> + <3s|1s>
         if(shell(na) < shell(nb)) then
            norm = SQRT((ZA**3)*(ZB**7)/7.5_wp)*(R**5)*0.0625_wp
            ovl = norm*(A4(ax)*B0(bx)-B4(bx)*A0(ax)+2.d0*(A3(ax)*B1(bx)-B3(bx)*A1(ax)))/sqrt(3.d0)
         else
            xx = za
            za = zb
            zb = xx
            ax = (za+zb)*R05
            bx = (zb-za)*R05
            norm = SQRT((ZA**3)*(ZB**7)/7.5_wp)*(R**5)*0.0625_wp
            ovl = norm*(A4(ax)*B0(bx)-B4(bx)*A0(ax)+2.d0*(A3(ax)*B1(bx)-B3(bx)*A1(ax)))/sqrt(3.d0)
         endif
      case(6)  ! <2s|3s> + <3s|2s>
         if(shell(na) < shell(nb)) then
            norm = SQRT((za**5)*(zb**7)/7.5_wp)*(R**6)*0.03125_wp
            ovl = norm*(A5(ax)*B0(bx)+A4(ax)*B1(bx)-2d0*(A3(ax)*B2(bx)+A2(ax)*B3(bx))+A1(ax)*B4(bx)+A0(ax)*B5(bx))/3.d0
         else
            xx = za
            za = zb
            zb = xx
            ax = (za+zb)*R05
            bx = (zb-za)*R05
            norm = SQRT((za**5)*(zb**7)/7.5_wp)*(R**6)*0.03125_wp
            ovl = norm*(A5(ax)*B0(bx)+A4(ax)*B1(bx)-2.0_wp*(A3(ax)*B2(bx)+A2(ax)*B3(bx))+A1(ax)*B4(bx)+A0(ax)*B5(bx))/3.d0
         endif
      case(9) ! <3s|3>
         norm = sqrt((ZA*ZB*R*R)**7)/1440.d0
         ovl = norm*(A6(ax)*B0(bx)-3.d0*(A4(ax)*B2(bx)-A2(ax)*B4(bx))-A0(ax)*B6(bx))
      end select
   endif
end subroutine ssovl


!****************************************
!* A(x) auxiliary integrals             *
!* Quantenchemie - Ein Lehrgang Vol 5   *
!* p. 570  eq. 11.4.14                  *
!****************************************

real(wp) pure function A0(x)
! Hilfsintegral A_0
implicit none
real(wp), intent(in) :: x
A0 = exp(-x)/x
return
end function

real(wp) pure function A1(x)
! Hilfsintegral A_1
implicit none
real(wp), intent(in) :: x
A1 = ((1+x)*exp(-x))/(x**2)
return
end function


real(wp) pure function A2(x)
! Hilfsintegral A_2
implicit none
real(wp), intent(in) :: x
A2 = ((2d0+2d0*x+x**2)*exp(-x))/x**3
return
end function


real(wp) pure function A3(x)
! Hilfsintegral A_3
implicit none
real(wp), intent(in) :: x
real(wp) xx
real(wp) x2, x3, x4
x2 = x*x
x3 = x2*x
x4 = x3*x
xx = (6d0+6d0*x+3d0*x2+x3)
A3 = (xx*exp(-x))/x4
return
end function


real(wp) pure function A4(x)
! Hilfsintegral A_4
implicit none
real(wp), intent(in) :: x
real(wp) xx
real(wp) x2, x3, x4, x5
x2 = x*x
x3 = x2*x
x4 = x3*x
x5 = x4*x
xx = (24d0+24d0*x+12d0*x2+4d0*x3+x4)
A4 = (xx*exp(-x))/x5
return
end function

real(wp) pure function A5(x)
! Hilfsintegral A_5
implicit none
real(wp), intent(in) :: x
real(wp) xx
real(wp) x2, x3, x4, x5, x6
x2 = x*x
x3 = x2*x
x4 = x3*x
x5 = x4*x
x6 = x5*x
xx = (120d0+120d0*x+60d0*x2+20d0*x3+5d0*x4+x5)
A5 = (xx*exp(-x))/x6
return
end function

real(wp) pure function A6(x)
! Hilfsintegral A_6
implicit none
real(wp), intent(in) :: x
real(wp) xx
real(wp) x2, x3, x4, x5, x6, x7
x2 = x*x
x3 = x2*x
x4 = x3*x
x5 = x4*x
x6 = x5*x
x7 = x6*x
xx = (720d0+720d0*x+360d0*x2+120d0*x3+30d0*x4+6d0*x5+x6)
A6 = (xx*exp(-x))/x7
return
end function



!**************************************
!* B(x) auxiliary integrals           *
!* Quantenchemie - Ein Lehrgang Vol 5 *
!* p. 570  eq. 11.4.14b               *
!**************************************


real(wp) pure function B0(x)
   real(wp), intent(in) :: x
   B0 = (exp(x)-exp(-x))/x
end function

real(wp) pure function B1(x)
   real(wp), intent(in) :: x
   real(wp) x2, x3
   x2 = x*x
   x3 = x2*x
   B1 = ((1.0_wp-x)*exp(x)-(1.0_wp+x)*exp(-x))/x2
end function

real(wp) pure function B2(x)
   real(wp), intent(in) :: x
   real(wp) x2, x3
   x2 = x*x
   x3 = x2*x
   B2 = (((2.0_wp-2*x+x2)*exp(x)) - ((2.0_wp+2.0_wp*x+x2)*exp(-x)))/x3
end function

real(wp) pure function B3(x)
   real(wp), intent(in) :: x
   real(wp) xx, yy
   real(wp) x2, x3, x4
   x2 = x*x
   x3 = x2*x
   x4 = x3*x
   xx = (6.0_wp-6.0_wp*x+3.0_wp*x2-x3)*exp(x)/x4
   yy = (6.0_wp+6.0_wp*x+3.0_wp*x2+x3)*exp(-x)/x4
   B3 = xx-yy
end function


real(wp) pure function B4(x)
   real(wp), intent(in) :: x
   real(wp) xx, yy
   real(wp) x2, x3, x4, x5
   x2 = x*x
   x3 = x2*x
   x4 = x3*x
   x5 = x4*x
   xx = (24.0_wp-24.0_wp*x+12.0_wp*x2-4.0_wp*x3+x4)*exp(x)/x5
   yy = (24.0_wp+24.0_wp*x+12.0_wp*x2+4.0_wp*x3+x4)*exp(-x)/x5
   B4 = xx-yy
end function

real(wp) pure function B5(x)
   real(wp), intent(in) :: x
   real(wp) xx, yy
   real(wp) x2, x3, x4, x5, x6
   x2 = x*x
   x3 = x2*x
   x4 = x3*x
   x5 = x4*x
   x6 = x5*x
   xx = (120.0_wp-120*x+60*x2-20*x3+5*x4-x5)*exp(x)/x6
   yy = (120.0_wp+120*x+60*x2+20*x3+5*x4+x5)*exp(-x)/x6
   B5 = xx-yy
end function

real(wp) function B6(x)
   real(wp), intent(in) :: x
   real(wp) x2, x3, x4, x5, x6, x7, yy, xx
   x2 = x*x
   x3 = x2*x
   x4 = x3*x
   x5 = x4*x
   x6 = x5*x
   x7 = x6*x
   xx = (720.0_wp - 720.0_wp*x+ 360.0_wp*x2 - 120.0_wp*x3 + 30.0_wp*x4 - 6.0_wp*x5 + x6)*exp(x)/x7
   yy = (720.0_wp + 720.0_wp*x + 360.0_wp*x2 + 120.0_wp*x3 + 30.0_wp*x4 + 6.0_wp*x5 + x6)*exp(-x)/x7
   B6 = xx-yy
end function


real(wp) function bint(x, k)
! calculates B_k(x)
! general summation formula
! 'infinite' sum is numerically unstable. 12 terms seem
! accurate enough
implicit none
real(wp), intent(in) :: x
real(wp) xx, yy
integer, intent(in) :: k
integer i
bint = 0

if(abs(x).lt.1e-6) then
do i = 0, k
   bint = (1.d0+(-1d0)**i)/(dble(i)+1.d0)
end do
return
endif

do i = 0, 12
xx = 1d0-((-1d0)**(k+i+1))
yy = dble(fact(i))*dble((k+i+1))
bint = bint+xx/yy*(-x)**i
enddo


end function bint


! faculty function
integer(wp) function fact(N)
implicit none
integer j, n
fact = 1
do j = 2, n
  fact = fact*j
enddo
return

end




subroutine gsovl(r, iat, jat, iz, xza, xzb, g)
   ! GRADIENT
   ! calculates the s-type overlap integral over 1s, 2s, 3s slater functions
   ! ovl = overlap integral
   ! za  = slater exponent atom A
   ! zb  = slater exponent atom B
   ! R   = distance between atom A and B
   integer ii, shell(72)
   logical debug
   real(wp) ax, bx, R05, za, zb, R
   integer na, nb
   data shell/                 &
   !          h, he
            1, 1               &
   !         li-ne
            , 2, 2, 2, 2, 2, 2, 2, 2, &
   !         na-ar
            3, 3, 3, 3, 3, 3, 3, 3,  &
   ! 4s, 5s will be treated as 3s
   !         k-rn , no f-elements
            54*3/
   ! ...
   real(wp) g, Fa, Fb
   !--------------------- set exponents ---------------------------------------
   real(wp) xza, xzb
   real(wp) xx
   integer iat, jat, iz(*)
   logical lsame

   za = xza
   zb = xzb
   na = iz(iat)
   nb = iz(jat)
   !----------------------------------------------------------------------------

   debug = .false.
   !debug = .true.

   ! ii selects kind of ovl by multiplying the shell
   ! kind    <1s|1s>  <2s|1s>  <2s|2s>  <1s|3s>  <2s|3s>  <3s|3s>
   ! case:      1        2        4       3        6         9
   !
   ii = shell(na)*shell(nb)
   if(debug) write(*, *) 'gshell', ii
   R05 = R*0.5_wp
   ax = (za+zb)*R05
   Fa = (za+zb)
   bx = (zb-za)*R05
   Fb = (zb-za)
   lsame = .false.
   !
   ! same elements
   if(abs(za-zb) < 0.1) then
      lsame = .true.
      ! call arguments: gtype(exponents, argumentDeriv., distance, gradient, (Switch shell), sameElement)
      select case (ii)
      case (1)
         call g1s1s(za, zb, Fa, Fb, R, g, lsame)
      case (2)
         if(shell(na) < shell(nb)) then
            call g2s1s(za, zb, Fa, Fb, R, g, .false., lsame)
         else
            xx = za
            za = zb
            zb = xx
            call g2s1s(za, zb, Fa, Fb, R, g, .true., lsame)
         end if
      case (4)
         call g2s2s(za, zb, Fa, Fb, R, g, lsame)
      case(3)
         if(shell(na) < shell(nb)) then
            call g1s3s(za, zb, Fa, Fb, R, g, .false., lsame)
         else
            xx = za
            za = zb
            zb = xx
            call g1s3s(za, zb, Fa, Fb, R, g, .true., lsame)
         end if
      case(6)
         if(shell(na) < shell(nb)) then
            call g2s3s(za, zb, Fa, Fb, R, g, .false., lsame)
         else
            xx = za
            za = zb
            zb = xx
            call g2s3s(za, zb, Fa, Fb, R, g, .true., lsame)
         end if
      case(9)
         call g3s3s(za, zb, Fa, Fb, R, g, lsame)
      end select
   else ! different elements
      lsame = .false.
      select case (ii)
      case (1)
         call g1s1s(za, zb, Fa, Fb, R, g, lsame)
      case (2)  ! <1s|2s>
         if(shell(na) < shell(nb)) then
            call g2s1s(za, zb, Fa, Fb, R, g, .false., lsame)
         else
            xx = za
            za = zb
            zb = xx
            call g2s1s(za, zb, Fa, Fb, R, g, .true., lsame)
         end if
      case (4) ! <2s|2s>
         call g2s2s(za, zb, Fa, Fb, R, g, lsame)
      case(3)  ! <1s|3s> + <3s|1s>
         if(shell(na) < shell(nb)) then
            call g1s3s(za, zb, Fa, Fb, R, g, .false., lsame)
         else
            xx = za
            za = zb
            zb = xx
            call g1s3s(za, zb, Fa, Fb, R, g, .true., lsame)
         end if
      case(6)  ! <2s|3s> + <3s|2s>
         if(shell(na) < shell(nb)) then
            call g2s3s(za, zb, Fa, Fb, R, g, .false., lsame)
         else
            xx = za
            za = zb
            zb = xx
            call g2s3s(za, zb, Fa, Fb, R, g, .true., lsame)
         end if
      case(9) ! <3s|3>
         call g3s3s(za, zb, Fa, Fb, R, g, lsame)
      end select
   end if

end subroutine gsovl


!-------------------------------------------------------------
! Maple was used to find the analy. derivatives of
! the slater integrals (expressions over A, B aux. integrals)
! Optimized fortran code by maple with some human-corrections
!-------------------------------------------------------------
subroutine g1s1s(za, zb, Fa, Fb, R, g, sameElement)
   ! slater overlap derv.
   ! derivative of explicit integral expression
   ! using maple
   implicit real(wp) (t)
   real(wp) za, zb, Fa, Fb
   real(wp) g, R
   logical sameElement

   if(sameElement) then
      t1 = za ** 2
      t3 = zb ** 2
      t5 = t1 * za * t3 * zb
      t6 = R ** 2
      t7 = t6 ** 2
      t10 = Fa * R
      t14 = exp(-0.5_wp * t10)
      t17 = sqrt(t5 * t7 * t6)
      g = -(1.0_wp/3.0_wp) * t5 * t7 / Fa * (0.2D1 + t10) * t14 / t17
   else
      t1 = za ** 2
      t3 = zb ** 2
      t5 = t1 * za * t3 * zb
      t6 = Fb ** 2
      t7 = Fb * R
      t8 = 0.5_wp * t7
      t9 = exp(t8)
      t12 = exp(-t8)
      t15 = t6 * Fa
      t22 = Fa ** 2
      t23 = t22 * t9
      t27 = t22 * t12
      t31 = t6 * Fb
      t32 = R * t31
      t37 = t22 * Fa
      t38 = R * t37
      t43 = R ** 2
      t44 = t43 * t31
      t51 = t43 * t37
      t56 = 0.4D1 * t6 * t9 - 0.4D1 * t6 * t12 + 0.2D1 * t15 * R * t9 -          &
      0.2D1 * t15 * R * t12 - 0.4D1 * t23 + 0.2D1 * t23 * t7 + 0.4D1 * t27       &
      + 0.2D1 * t27 * t7 - 0.2D1 * t32 * t9 - 0.2D1 * t32 * t12 - 0.2D1 * t38 *  &
      t9 + 0.2D1 * t38 * t12 - 0.1D1 * t44 * Fa * t9 - 0.1D1                     &
      * t44 * Fa * t12 + t51 * t9 * Fb + t51 * t12 * Fb
      t61 = exp(-0.5_wp * Fa * R)
      t62 = t43 ** 2
      t65 = sqrt(t5 * t62 * t43)
      g = -0.2D1 * t5 * R * t56 * t61 / t65 / t31 / t37
   end if


end subroutine g1s1s


subroutine g2s1s(za, zb, Fa, Fb, R, g, switch, lsame)
   ! slater overlap derv.
   ! derivative of explicit integral expression
   ! using maple
   implicit real(wp) (t)
   real(wp) za, zb, Fa, Fb
   real(wp) g, R, norm
   logical switch
   logical lsame
   norm = (1d0/24d0)*sqrt(za**3*zb**5*3d0)

   if(switch) then
      Fb = -Fb
   endif

   if(lsame) then

      t1 = Fa * R
      t3 = exp(-0.5000000000_wp * t1)
      t6 = Fa ** 2
      t7 = R ** 2
      g = -0.1000000000D-8 * R * t3 * (0.5333333380D10 + 0.2666666670D10 &
      * t1 + 0.1333333333D10 * t6 * t7) / t6
      g = g*norm

   else

      t3 = exp(-0.5000000000_wp * Fa * R)
      t4 = Fa ** 2
      t5 = t4 * Fa
      t6 = Fb * R
      t7 = 0.5000000000_wp * t6
      t8 = exp(t7)
      t9 = t5 * t8
      t11 = Fb ** 2
      t12 = t11 * Fa
      t15 = exp(-t7)
      t18 = t4 ** 2
      t19 = R * t18
      t22 = t11 ** 2
      t29 = Fb * t4
      t36 = R ** 2
      t37 = t36 * t18
      t44 = t36 * R
      t48 = -0.12D2 * t9 + 0.4D1 * t12 * t8 - 0.4D1 * t12 * t15 &
      - 0.6D1 * t19 * t8 - 0.6D1 * t22 * t8 * R - 0.6D1 * t22 * t15 &
      * R + 0.4D1 * t29 * t15 - 0.4D1 * t29 * t8 + 0.6D1 * t19 * t15 &
      + 0.2D1 * t37 * t8 * Fb + 0.4D1 * t37 * t15 * Fb + t44 * t18 * t15 * t11
      t49 = t5 * t15
      t51 = t11 * Fb
      t58 = t51 * Fa
      t59 = R * t8
      t76 = t36 * t15
      t79 = t22 * Fa
      t87 = 0.12D2 * t49 - 0.12D2 * t51 * t15 - 0.1D1 * t22 * t4 * t15 * t44 &
      + 0.4D1 * t58 * t59 - 0.8D1 * t58 * R * t15 + 0.4D1 * t9 * t6 + 0.8D1 * &
      t49 * t6 + 0.2D1 * t49 * t11 * t36 + 0.4D1 * t11 * t4 * t59 - 0.2D1 * t51 &
      * t4 * t76 - 0.2D1 * t79 * t36 * t8 - 0.4D1 * t79 * t76 + 0.12D2 * t51 * t8
      g = -0.16D2 * t3 * (t48 + t87) / t36 / t22 / t18
      g = g*norm
   end if

end subroutine g2s1s

subroutine g2s2s(za, zb, Fa, Fb, R, g, SameElement)
   ! slater overlap derv.
   ! derivative of explicit integral expression
   ! using maple
   implicit real(wp) (t)
   real(wp) za, zb, Fa, Fb
   real(wp) g, R, norm
   logical SameElement

   norm = 1d0/(16d0*3d0)*SQRT((ZA*ZB)**5)

   if(SameElement) then

      t2 = R ** 2
      t5 = Fa ** 2
      t9 = t5 * Fa
      t10 = t2 ** 2
      t16 = exp(-Fa * R / 0.2D1)
      g = (-0.4266666666D2 * R - 0.2133333333D2 * Fa * t2 - 0.2133333333D1 &
          * t5 * t2 * R - 0.1066666666D1 * t9 * t10) * t16 / t9
      g = g*norm


   else
      t1 = R ** 2
      t3 = 0.3840000000D3 * t1 * Fb
      t4 = t1 * R
      t5 = Fb ** 2
      t7 = 0.6400000000D2 * t4 * t5
      t8 = 0.7680000000D3 * R
      t10 = Fa ** 2
      t11 = t10 ** 2
      t12 = t11 * Fa
      t14 = Fb * R
      t15 = 0.768000000D3 * t14
      t17 = 0.1280000000D3 * t5 * t1
      t21 = 0.256000000D3 * t5 * R
      t22 = t5 * Fb
      t24 = 0.1280000000D3 * t22 * t1
      t26 = t10 * Fa
      t28 = t5 ** 2
      t30 = 0.1280000000D3 * t1 * t28
      t32 = 0.256000000D3 * t22 * R
      t33 = 0.512000000D3 * t5
      t34 = t28 * Fb
      t36 = 0.6400000000D2 * t4 * t34
      t40 = 0.768000000D3 * t28 * R
      t42 = 0.3840000000D3 * t1 * t34
      t45 = 0.1536000000D4 * t28
      t47 = 0.7680000000D3 * t34 * R
      t51 = exp(-0.5_wp * Fa * R)
      t53 = 0.5_wp * t14
      t54 = exp(-t53)
      t68 = exp(t53)
      g = (((t3 + t7 + t8) * t12 + (0.1536000000D4 + t15 + t17) * t11 + &
      (-t21 - t24) * t26 + (t30 - t32 - t33 + t36) * t10 + (t40 + t42) *&
      Fa + t45 + t47) * t51 * t54 + ((t3 - t8 - t7) * t12 + (-0.1536000000D4 &
      + t15 - t17) * t11 + (-t24 + t21) * t26 + (-t30 + t33 - t32 &
      + t36) * t10 + (-t40 + t42) * Fa + t47 - t45) * t51 * t68) / t1 /  &
      t12 / t34


      g = g*norm
   endif

end subroutine g2s2s


subroutine g1s3s(za, zb, Fa, Fb, R, g, switch, lsame)
   ! slater overlap derv.
   ! derivative of explicit integral expression
   ! using maple
   implicit real(wp) (t)
   real(wp) za, zb, Fa, Fb
   real(wp) g, R, norm
   logical switch
   logical lsame

   if(switch) Fb = -Fb

   norm = SQRT((ZA**3)*(ZB**7)/7.5_wp)/(16d0*sqrt(3d0))

   if(lsame) then

      t1 = Fa * R
      t3 = exp(-0.5000000000_wp * t1)
      t4 = R ** 2
      g = -0.1600000000D1 * t3 * t4 * R * (0.2D1 + t1) / Fa
      g = g*norm
   else

      t3 = exp(-0.5000_wp * Fa * R)
      t4 = Fb ** 2
      t5 = t4 ** 2
      t6 = t5 * Fb
      t7 = t6 * Fa
      t8 = R ** 2
      t9 = Fb * R
      t10 = 0.50_wp * t9
      t11 = exp(t10)
      t15 = exp(-t10)
      t16 = t8 * t15
      t19 = Fa ** 2
      t21 = t8 * R
      t22 = t21 * t15
      t25 = t19 * Fa
      t27 = t8 ** 2
      t31 = t19 ** 2
      t32 = t31 * Fa
      t33 = t8 * t32
      t45 = t4 * Fb
      t48 = t31 * t15
      t55 = t4 * t25
      t56 = t11 * R
      t59 = t15 * R
      t62 = t5 * Fa
      t73 = -0.6D1 * t7 * t8 * t11 - 0.18D2 * t7 * t16 - 0.6D1 * t6 * t19 &
      * t22 - 0.1D1 * t6 * t25 * t27 * t15 + 0.6D1 * t33 * t11 * Fb + 0.18D2 &
      * t33 * t15 * Fb + 0.6D1 * t21 * t32 * t15 * t4 + t27 * t32* t15 * t45 &
      + 0.2D1 * t48 * t45 * t21 + 0.12D2 * t48 * t4 * t8 + 0.12D2 * t55 * t56 &
      + 0.12D2 * t55 * t59 + 0.12D2 * t62 * t56 - 0.36D2 * t62 * t59 - 0.12D2 &
      * t5 * t19 * t16 - 0.2D1 * t5 * t25 * t22
      t74 = t31 * t11
      t79 = t45 * t19
      t92 = R * t32
      t95 = t45 * Fa
      t100 = Fb * t25
      t111 = 0.12D2 * t74 * t9 + 0.36D2 * t48 * t9 + 0.12D2 * t79 * t56 - 0.12D2  &
      * t79 * t59 + 0.48D2 * t5 * t11 - 0.24D2 * t6 * t11 * R - 0.24D2 * t6 * t15 &
      * R - 0.24D2 * t92 * t11 + 0.24D2 * t95 * t11 - 0.24D2 * t95 * t15 + 0.24D2 &
      * t100 * t15 + 0.24D2 * t92 * t15 - 0.24D2 * t100 * t11 - 0.48D2 * t5 * t15 &
      - 0.48D2 * t74 + 0.48D2 * t48
      g = -0.32D2 * t3 * (t73 + t111) / t8 / t6 / t32

      g = g*norm
   endif

end subroutine g1s3s


subroutine g2s3s(za, zb, Fa, Fb, R, g, switch, lsame)
   ! slater overlap derv.
   ! derivative of explicit integral expression
   ! using maple
   implicit real(wp) (t)
   real(wp) za, zb, Fa, Fb
   real(wp) g, R, norm
   logical switch
   logical lsame
   norm = sqrt((za**5)*(zb**7)/7.5_wp)/96.d0

   if(switch) Fb = -Fb

   if(lsame) then
      t1 = Fa * R
      t3 = exp(-0.5000000000_wp * t1)
      t6 = Fa ** 2
      t7 = R ** 2
      t14 = t6 ** 2
      t15 = t7 ** 2
      g = -0.2000000000D-8 * R * t3 * (0.1280000000D12 + 0.6400000000D11 &
      * t1 + 0.1280000000D11 * t6 * t7 + 0.1066666670D10 * t6 * Fa * t7 &
      * R + 0.533333333D9 * t14 * t15) / t14
      g = g*norm
   else

      t3 = exp(-0.5_wp * Fa * R)
      t4 = Fb ** 2
      t5 = t4 ** 2
      t6 = Fa ** 2
      t7 = t6 * Fa
      t8 = t5 * t7
      t9 = R ** 2
      t11 = 0.50_wp * Fb * R
      t12 = exp(t11)
      t13 = t9 * t12
      t16 = t6 ** 2
      t17 = t16 * Fa
      t18 = exp(-t11)
      t21 = t5 * Fb
      t28 = t9 * t18
      t32 = t9 * R
      t33 = t32 * t18
      t36 = t5 * t4
      t38 = t9 ** 2
      t39 = t38 * t18
      t41 = t21 * Fa
      t42 = R * t12
      t45 = t16 * t6
      t46 = t4 * Fb
      t49 = t46 * t16
      t52 = -0.6D1 * t8 * t13 + 0.120D3 * t17 * t18 + 0.120D3 * t21 * t18 &
       - 0.120D3 * t17 * t12 - 0.120D3 * t21 * t12 - 0.6D1 * t8 * t28 - 0.2D1 &
      * t5 * t16 * t33 + t36 * t7 * t39 - 0.48D2 * t41 * t42 + t45 * t46 * t39 - 0.6D1 * t49 * t13
      t54 = R * t18
      t60 = t46 * t6
      t63 = Fb * t16
      t66 = t5 * Fa
      t69 = t4 * t7
      t72 = t36 * t6
      t75 = t32 * t12
      t78 = Fb * t9
      t84 = Fb * t17
      t87 = -0.24D2 * t46 * t7 * t54 - 0.24D2 * t5 * t6 * t42 + 0.24D2 *&
       t60 * t12 + 0.24D2 * t63 * t18 - 0.24D2 * t66 * t12 - 0.24D2 * t69 &
      * t18 + 0.9D1 * t72 * t33 + 0.3D1 * t72 * t75 + 0.24D2 * t78 * t45 &
      * t12 - 0.6D1 * t49 * t28 + 0.48D2 * t84 * t42
      t102 = t21 * t6
      t105 = t4 * t17
      t113 = t45 * t4
      t118 = 0.72D2 * t84 * t54 + 0.72D2 * t41 * t54 + 0.36D2 * t78 * t45 &
      * t18 + 0.2D1 * t46 * t17 * t33 + 0.24D2 * t4 * t16 * t42 - 0.6D1   &
      * t102 * t13 - 0.6D1 * t105 * t13 + 0.18D2 * t105 * t28 + 0.2D1 &
      * t21 * t7 * t33 - 0.3D1 * t113 * t75 + 0.9D1 * t113 * t33
      t121 = t36 * Fa
      t130 = R * t45
      t145 = 0.18D2 * t102 * t28 + 0.24D2 * t121 * t13 + 0.36D2 * t121 * &
       t28 - 0.24D2 * t60 * t18 - 0.24D2 * t63 * t12 + 0.60D2 * t130 * t18 &
      + 0.60D2 * t36 * t18 * R + 0.24D2 * t69 * t12 + 0.60D2 * t36 * t12 * &
      R - 0.60D2 * t130 * t12 + 0.24D2 * t66 * t18
      g = 0.128D3 * t3 * (t52 + t87 + t118 + t145) / t9 / t36 / t45

   g = g*norm
   endif

end subroutine g2s3s


subroutine g3s3s(za, zb, Fa, Fb, R, g, SameElement)
   ! slater overlap derv.
   ! derivative of explicit integral expression
   ! using maple
   implicit real(wp) (t)
   real(wp) za, zb, Fa, Fb
   real(wp) g, R, norm
   logical SameElement

   norm = sqrt((ZA*ZB)**7)/1440.d0

   if(SameElement) then

      t1 = Fa * R
      t3 = exp(-0.5000000000_wp * t1)
      t5 = Fa ** 2
      t6 = t5 ** 2
      t7 = t6 * Fa
      t8 = R ** 2
      t9 = t8 ** 2
      g = -0.2000000000D-8 * t3 * R * (0.457142857D9 * t7 * t9 &
      * R + 0.7680000000D12 * t1 + 0.1536000000D12 * t5 * t8 &
      + 0.1280000000D11 * t5 * Fa * t8 * R + 0.914285715D9 * t6 * t9 + 0.1536000000D13) / t7

      g = g*norm
   else

      t3 = exp(-0.5000000000_wp * Fa * R)
      t4 = Fa ** 2
      t5 = t4 ** 2
      t6 = t5 * t4
      t7 = Fb * R
      t8 = 0.5000000000_wp * t7
      t9 = exp(-t8)
      t10 = t6 * t9
      t13 = Fb ** 2
      t14 = t13 * Fb
      t15 = t13 ** 2
      t16 = t15 * t14
      t17 = R ** 2
      t18 = t17 * R
      t19 = t16 * t18
      t23 = exp(t8)
      t24 = t6 * t23
      t27 = t5 * Fa
      t28 = t27 * t13
      t29 = R * t23
      t32 = t6 * t13
      t33 = t17 * t9
      t36 = t15 * Fb
      t37 = t4 * t36
      t38 = t9 * R
      t43 = t17 * t23
      t46 = t4 * Fa
      t47 = t5 * t46
      t48 = t47 * t18
      t52 = t47 * t17
      t65 = 0.120D3 * t10 * t7 - 0.12D2 * t19 * t4 * t9 + 0.120D3 &
      * t24 * t7 + 0.24D2 * t28 * t29 + 0.24D2 * t32 * t33 + 0.24D2 * t37 &
      * t38 - 0.24D2 * t28 * t38 - 0.24D2 * t32 * t43 - 0.12D2 * t48 * t13 &
      * t23 + 0.60D2 * t52 * t23 * Fb + 0.12D2 * t48 * t13 * t9 + 0.60D2 &
      * t52 * t9 * Fb - 0.12D2 * t19 * t4 * t23
      t66 = t17 ** 2
      t67 = t16 * t66
      t74 = t27 * t14
      t77 = t6 * t14
      t78 = t18 * t23
      t81 = t46 * t15
      t86 = t27 * t15
      t89 = t5 * t36
      t90 = t18 * t9
      t97 = t46 * t36
      t104 = -0.1D1 * t67 * t46 * t9 - 0.1D1 * t67 * t46 * t23 - 0.12D2 &
      * t74 * t43 + 0.2D1 * t77 * t78 - 0.24D2 * t81 * t29 + 0.24D2 * t81 &
      * t38 + 0.2D1 * t86 * t78 + 0.2D1 * t89 * t90 - 0.2D1 * t86 * t90 &
      + 0.24D2 * t37 * t29 + 0.12D2 * t97 * t33 + 0.2D1 * t89 * t78 - 0.12D2 * t74 * t33
      t108 = t5 * t14
      t111 = t15 * t13
      t112 = t111 * t4
      t117 = t111 * t46
      t122 = t111 * Fa
      t129 = t4 * t15
      t132 = t47 * R
      t139 = 0.2D1 * t77 * t90 - 0.24D2 * t108 * t38 + 0.24D2 * t112 * t43 &
      - 0.24D2 * t112 * t33 + 0.2D1 * t117 * t78 - 0.2D1 * t117 * t90 + 0.120D3 &
      * t122 * t29 - 0.120D3 * t122 * t38 + 0.12D2 * t97 * t43 - 0.48D2 * t129 &
      * t23 + 0.120D3 * t132 * t9 - 0.120D3 * t132 * t23 + 0.240D3 * t111 * t23
      t140 = t47 * t66
      t145 = t16 * R
      t150 = t16 * t17
      t160 = t5 * t13
      t170 = t140 * t14 * t23 + t140 * t14 * t9 - 0.120D3 * t145 * t9 - 0.24D2 &
      * t108 * t29 - 0.60D2 * t150 * Fa * t23 - 0.240D3 * t111 * t9 - 0.240D3 &
      * t24 + 0.240D3 * t10 + 0.48D2 * t129 * t9 - 0.48D2 * t160 * t9 + 0.48D2 &
      * t160 * t23 - 0.120D3 * t145 * t23 - 0.60D2 * t150 * Fa * t9
      g = -0.768D3 * t3 * (t65 + t104 + t139 + t170) / t17 / t47 / t16

      g = g*norm
   endif

end subroutine g3s3s

end module dftd3_gcp