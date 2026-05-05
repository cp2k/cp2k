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

!> @file tblite/basis/qvszp.f90
!> Provides coefficients and exponents for the q-vSZP basis set

!> Implements the valence-only version of the charge-dependent 
!> valence single-zeta q-vSZP basis set based on: 
!>
!> Marcel Müller, Andreas Hansen, and Stefan Grimme, 
!> "An atom-in-molecule adaptive polarized valence single-ζ 
!> atomic orbital basis for electronic structure calculations"
!> J. Chem. Phys. 159, 164108 (2023). DOI: 10.1063/5.0172373
!> 
!> Final coefficient and exponents taken from: 
!>
!> Marcel Müller, Thomas Froitzheim, Andreas Hansen, and Stefan Grimme, 
!> "Advanced Charge Extended Hückel (CEH) Model and a Consistent 
!> Adaptive Minimal Basis Set for the Elements Z = 1–103"
!> J. Phys. Chem. A 2024, 128, 49, 10723–10736, DOI: 10.1021/acs.jpca.4c06989
!>
!> Thomas Froitzheim, Marcel Müller, Andreas Hansen, and Stefan Grimme, 
!> "The Bond Capacity Electronegativity Equilibration Charge Model (EEQBC) 
!> for the Elements Z=1–103"
!> J. Chem. Phys. 2025, 162, 214109, DOI: 10.1063/5.0268978
module tblite_basis_qvszp
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io_constants, only : pi
   use mctc_io, only : structure_type
   use mctc_ncoord, only : ncoord_type, new_ncoord, cn_count
   use mctc_cutoff, only : get_lattice_points
   use tblite_basis_cache, only : basis_cache, cgto_cache
   use tblite_basis_type, only : basis_type, integral_cutoff, cgto_type, cgto_container
   use tblite_integral_overlap, only : overlap_cgto, overlap_grad_cgto, msao
   use tblite_output_format, only: format_string
   use tblite_wavefunction_type, only : wavefunction_type
   implicit none
   private

   public :: new_qvszp_basis, new_qvszp_cgto, get_qvszp_parameters, maxg

   !> Maximal number of primitives per CGTO
   integer, parameter :: maxg = 12
   !> Maximal number of elements
   integer, parameter :: max_elem = 103
   !> Maximal number of shells
   integer, parameter :: max_shell = 4

   !> Default accuracy for basis set creation
   real(wp), parameter :: default_accuracy = 1.0_wp

   !> Two over pi
   real(wp), parameter :: top = 2.0_wp / pi
   !> Double factorial, see OEIS A001147
   real(wp), parameter :: dfactorial(8) = &
      & [1.0_wp,1.0_wp,3.0_wp,15.0_wp,105.0_wp,945.0_wp,10395.0_wp,135135.0_wp]

   !> Contracted Gaussian type basis function of qvszp type
   type, public, extends(cgto_type) :: qvszp_cgto_type
      !> Primitive-specific strength of the environment specific scaling
      real(wp) :: coeff1(maxg) = 0.0_wp
      !> Overall charge-dependence of the coefficient
      real(wp) :: k0 = 0.0_wp
      !> Squared charge-dependence of the coefficient
      real(wp) :: k1 = 0.0_wp
      !> CN-dependence of the coefficient
      real(wp) :: k2 = 0.0_wp
      !> Mixed charge- and CN-dependence of the coefficient
      real(wp) :: k3 = 0.0_wp
   contains
      !> Update CGTO cache
      procedure :: update => cgto_update
      !> Get normalization factor and its derivatives
      procedure :: get_normalization
      !> Get (scaled) coefficient of the CGTO
      procedure :: get_coeffs
      !> Get coefficient derivatives of the CGTO
      procedure :: get_coeff_derivs
      !> Get effective atomic charge
      procedure :: get_qeff
   end type qvszp_cgto_type

   !> Collection of information regarding the basis set of a system
   type, public, extends(basis_type) :: qvszp_basis_type 
      !> Coordination number for modifying the contraction coefficients
      class(ncoord_type), allocatable :: ncoord
   contains
      !> Update basis cache
      procedure :: update => basis_update
      !> Contract gradients w.r.t. the basis set specific charge and CN
      procedure :: get_basis_gradient
   end type qvszp_basis_type

   !> Number of shells
   integer, parameter :: nshell(max_elem) = [ &
      & 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
      & 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
      & 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, &
      & 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
      & 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, &
      & 4, 4, 4]

   !> Number of primitives per shell
   integer, parameter :: n_prim(max_shell, max_elem) = reshape([&
      & 8, 3, 0, 0, 8, 2, 0, 0, 5, 5, 2, 0, 6, 4, 2, 0, 6, 5, 2, 0, & ! -5
      & 6, 6, 3, 0, 6, 6, 3, 0, 6, 6, 3, 0, 6, 6, 2, 0, 6, 6, 2, 0, & ! -10
      & 4, 4, 2, 0, 4, 3, 2, 0, 5, 4, 2, 0, 5, 4, 2, 0, 5, 5, 2, 0, & ! -15
      & 5, 5, 2, 0, 5, 5, 2, 0, 5, 5, 2, 0, 4, 3, 2, 0, 4, 3, 3, 0, & ! -20
      & 5, 2, 6, 0, 5, 3, 6, 0, 5, 3, 6, 0, 5, 3, 6, 0, 5, 3, 6, 0, & ! -25
      & 5, 3, 6, 0, 5, 3, 6, 0, 5, 3, 6, 0, 5, 3, 6, 0, 5, 3, 2, 0, & ! -30
      & 5, 4, 2, 0, 6, 5, 2, 0, 6, 5, 2, 0, 6, 5, 2, 0, 6, 5, 2, 0, & ! -35
      & 6, 5, 2, 0, 4, 3, 2, 0, 4, 3, 3, 0, 5, 2, 6, 0, 5, 3, 6, 0, & ! -40
      & 5, 3, 6, 0, 5, 3, 6, 0, 5, 2, 6, 0, 5, 3, 6, 0, 5, 3, 6, 0, & ! -45
      & 5, 3, 6, 0, 5, 3, 6, 0, 5, 3, 2, 0, 5, 4, 2, 0, 6, 5, 2, 0, & ! -50
      & 6, 5, 2, 0, 6, 5, 2, 0, 6, 5, 2, 0, 6, 5, 2, 0, 4, 3, 2, 0, & ! -55
      & 4, 3, 3, 0, 4, 2, 5, 0, 7, 4, 7, 7, 7, 4, 6, 7, 7, 4, 6, 7, & ! -60
      & 7, 4, 7, 7, 7, 4, 6, 7, 7, 4, 6, 7, 7, 4, 6, 7, 7, 4, 7, 7, & ! -65
      & 7, 4, 6, 7, 7, 4, 6, 7, 7, 4, 6, 7, 7, 4, 6, 7, 7, 4, 6, 7, & ! -70
      & 7, 4, 6, 7, 4, 3, 5, 0, 5, 3, 5, 0, 5, 3, 5, 0, 5, 3, 5, 0, & ! -75
      & 5, 3, 5, 0, 5, 3, 5, 0, 5, 3, 5, 0, 5, 3, 5, 0, 5, 3, 2, 0, & ! -80
      & 5, 4, 2, 0, 5, 5, 2, 0, 6, 5, 2, 0, 6, 5, 2, 0, 6, 5, 2, 0, & ! -85
      & 6, 5, 2, 0, 9, 9, 4, 0, 9, 9, 4, 0, 7, 4, 6, 7, 7, 4, 6, 7, & ! -90
      & 7, 4, 6, 7, 7, 4, 6, 7, 7, 4, 6, 7, 7, 4, 6, 7, 7, 4, 6, 7, & ! -95
      & 7, 4, 6, 7, 7, 4, 6, 7, 7, 4, 6, 7, 7, 4, 6, 7, 7, 4, 6, 7, & ! -100
      & 7, 4, 6, 7, 7, 4, 6, 7, 7, 4, 6, 7], shape(n_prim))

   !> Angular momentum of each shell
   integer, parameter :: ang_shell(max_shell, max_elem) = reshape([&
      & 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, & ! -5
      & 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, & ! -10
      & 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, & ! -15
      & 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, & ! -20
      & 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, & ! -25
      & 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, & ! -30
      & 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, & ! -35
      & 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, & ! -40
      & 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, & ! -45
      & 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, & ! -50
      & 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, & ! -55
      & 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, & ! -60
      & 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, & ! -65
      & 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, & ! -70
      & 0, 1, 2, 3, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, & ! -75
      & 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, & ! -80
      & 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, & ! -85
      & 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 0, 0, 1, 2, 3, 0, 1, 2, 3, & ! -90
      & 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, & ! -95
      & 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, & ! -100
      & 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3], shape(ang_shell))

   integer, parameter :: principal_quantum_number(max_shell, max_elem) = reshape([&
      & 1, 2, 0, 0, 1, 2, 0, 0, 2, 2, 3, 0, 2, 2, 3, 0, 2, 2, 3, 0, & ! -5
      & 2, 2, 3, 0, 2, 2, 3, 0, 2, 2, 3, 0, 2, 2, 3, 0, 2, 2, 3, 0, & ! -10
      & 3, 3, 3, 0, 3, 3, 3, 0, 3, 3, 3, 0, 3, 3, 3, 0, 3, 3, 3, 0, & ! -15
      & 3, 3, 3, 0, 3, 3, 3, 0, 3, 3, 3, 0, 4, 4, 3, 0, 4, 4, 3, 0, & ! -20
      & 4, 4, 3, 0, 4, 4, 3, 0, 4, 4, 3, 0, 4, 4, 3, 0, 4, 4, 3, 0, & ! -25
      & 4, 4, 3, 0, 4, 4, 3, 0, 4, 4, 3, 0, 4, 4, 3, 0, 4, 4, 3, 0, & ! -30
      & 4, 4, 4, 0, 4, 4, 4, 0, 4, 4, 4, 0, 4, 4, 4, 0, 4, 4, 4, 0, & ! -35
      & 4, 4, 4, 0, 5, 5, 4, 0, 5, 5, 4, 0, 5, 5, 4, 0, 5, 5, 4, 0, & ! -40
      & 5, 5, 4, 0, 5, 5, 4, 0, 5, 5, 4, 0, 5, 5, 4, 0, 5, 5, 4, 0, & ! -45
      & 5, 5, 4, 0, 5, 5, 4, 0, 5, 5, 4, 0, 5, 5, 5, 0, 5, 5, 5, 0, & ! -50
      & 5, 5, 5, 0, 5, 5, 5, 0, 5, 5, 5, 0, 5, 5, 5, 0, 6, 6, 5, 0, & ! -55
      & 6, 6, 5, 0, 6, 6, 5, 0, 6, 6, 5, 4, 6, 6, 5, 4, 6, 6, 5, 4, & ! -60
      & 6, 6, 5, 4, 6, 6, 5, 4, 6, 6, 5, 4, 6, 6, 5, 4, 6, 6, 5, 4, & ! -65
      & 6, 6, 5, 4, 6, 6, 5, 4, 6, 6, 5, 4, 6, 6, 5, 4, 6, 6, 5, 4, & ! -70
      & 6, 6, 5, 4, 6, 6, 5, 0, 6, 6, 5, 0, 6, 6, 5, 0, 6, 6, 5, 0, & ! -75 
      & 6, 6, 5, 0, 6, 6, 5, 0, 6, 6, 5, 0, 6, 6, 5, 0, 6, 6, 5, 0, & ! -80
      & 6, 6, 6, 0, 6, 6, 6, 0, 6, 6, 6, 0, 6, 6, 6, 0, 6, 6, 6, 0, & ! -85
      & 6, 6, 6, 0, 7, 7, 6, 0, 7, 7, 6, 0, 7, 7, 6, 5, 7, 7, 6, 5, & ! -90
      & 7, 7, 6, 5, 7, 7, 6, 5, 7, 7, 6, 5, 7, 7, 6, 5, 7, 7, 6, 5, & ! -95
      & 7, 7, 6, 5, 7, 7, 6, 5, 7, 7, 6, 5, 7, 7, 6, 5, 7, 7, 6, 5, & ! -100
      & 7, 7, 6, 5, 7, 7, 6, 5, 7, 7, 6, 5], shape(principal_quantum_number))

   !> Parameter for overall charge dependence
   real(wp), parameter :: p_k0(max_elem) = [&
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & ! -4
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & ! -8
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & ! -12
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & ! -16
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & ! -20
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & ! -24
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & ! -28
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & ! -32
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & ! -36
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & ! -40
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & ! -44
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & ! -48
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & ! -52
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & ! -56
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & ! -60
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & ! -64
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & ! -68
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & ! -72
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & ! -76
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & ! -80
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & ! -84
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & ! -88
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & ! -92
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & ! -96
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp, & ! -100
      &  1.0000000000_wp,  1.0000000000_wp,  1.0000000000_wp]


   !> Parameter for quadratic charge dependence
   real(wp), parameter :: p_k1(max_elem) = [&
      &  0.2272460655_wp, -0.0596139193_wp, -0.8162219211_wp,  0.2615605768_wp, & ! -4
      &  0.2862843025_wp, -0.0173073504_wp,  0.2392693084_wp,  0.0118176989_wp, & ! -8
      &  0.0124716084_wp, -0.0393471359_wp, -0.3706940822_wp, -0.0865911547_wp, & ! -12
      &  0.2215897181_wp, -0.0315062343_wp,  0.3434101268_wp,  0.2263363232_wp, & ! -16
      & -0.0492370074_wp, -0.0538331339_wp, -0.5261134284_wp, -0.1165840386_wp, & ! -20
      &  0.0844447054_wp, -0.0516881520_wp,  0.0522662650_wp,  0.3083293812_wp, & ! -24
      & -0.0336221573_wp,  0.3426408669_wp,  0.5599083638_wp,  0.0303309809_wp, & ! -28
      & -0.1089879016_wp,  0.2495153629_wp, -0.1199126821_wp,  0.0268359876_wp, & ! -32
      &  0.4647810638_wp,  0.2773352545_wp,  0.0529580888_wp, -0.0772051837_wp, & ! -36
      & -0.5679581136_wp,  1.4183340595_wp,  0.1372990988_wp, -0.0341739234_wp, & ! -40
      & -0.0234715066_wp,  0.7684617366_wp, -0.0873748856_wp,  0.2345703772_wp, & ! -44
      &  0.3916056534_wp,  0.1175717490_wp, -0.0956929785_wp,  0.2552858006_wp, & ! -48
      &  0.0449687823_wp, -0.0474013379_wp,  0.5613791651_wp,  0.2526477340_wp, & ! -52
      & -0.0478303905_wp,  0.1495971686_wp, -0.8369485378_wp, -0.1061959147_wp, & ! -56
      & -0.0176044338_wp,  0.0502888557_wp, -0.0637994570_wp, -0.0635998097_wp, & ! -60
      & -0.0531518530_wp, -0.0968731360_wp, -0.1337199042_wp, -0.0542471443_wp, & ! -64
      & -0.0409157548_wp, -0.0500477946_wp, -0.0903409668_wp, -0.0624137643_wp, & ! -68
      & -0.0964656303_wp, -0.0547844435_wp, -0.0901499353_wp, -0.1044976019_wp, & ! -72
      &  0.0748615545_wp, -0.0409080635_wp, -0.0952820449_wp,  0.0245211400_wp, & ! -76
      & -0.0497361587_wp,  0.3030126210_wp, -0.2156726008_wp, -0.1016401687_wp, & ! -80
      & -0.0512111539_wp, -0.0482926025_wp,  0.3318716965_wp,  0.3413143387_wp, & ! -84
      & -0.0400146567_wp,  0.2200034940_wp,  0.0876924655_wp,  0.1321612468_wp, & ! -88
      &  0.2949305821_wp,  0.4665862315_wp, -0.1050478058_wp,  0.1694553742_wp, & ! -92
      & -0.1161440443_wp,  0.2497086409_wp,  0.0496197475_wp, -0.0424576824_wp, & ! -96
      & -0.0417067664_wp, -0.0367557847_wp, -0.0348215574_wp, -0.0364418652_wp, & ! -100
      & -0.0412106965_wp, -0.0350811514_wp, -0.0609692849_wp]

   !> Parameter for square-root CN dependence
   real(wp), parameter :: p_k2(max_elem) = [&
      &  0.4237989919_wp,  0.9054329237_wp,  0.7273599522_wp,  1.1780944589_wp, & ! -4
      &  1.2432066425_wp,  0.9021610088_wp,  0.1729209683_wp,  0.3137107221_wp, & ! -8
      &  0.1481611609_wp,  0.0315012512_wp,  0.6948564260_wp,  1.5731412842_wp, & ! -12
      &  1.4646004177_wp,  0.7338489851_wp,  0.2672743388_wp,  0.1255489520_wp, & ! -16
      &  0.3795992335_wp,  0.3132264025_wp,  1.2745899845_wp,  1.9855778911_wp, & ! -20
      &  0.3176116613_wp,  1.5288514307_wp,  1.4881571248_wp,  1.3809899002_wp, & ! -24
      &  2.3307533266_wp,  0.3765445274_wp,  0.3755945109_wp, -0.4118221596_wp, & ! -28
      &  0.2692468691_wp,  1.0593568196_wp,  1.8497298969_wp,  1.6964058522_wp, & ! -32
      &  0.1456480496_wp,  0.2315384838_wp,  0.1852541944_wp,  0.3599891154_wp, & ! -36
      &  1.2633917188_wp,  1.6791572024_wp,  1.5055168142_wp,  1.1180269980_wp, & ! -40
      &  2.0722043983_wp,  1.1736062932_wp,  1.5110217881_wp,  1.5285170139_wp, & ! -44
      &  1.6948373822_wp, -0.5176634609_wp,  1.5077143340_wp,  1.7556728540_wp, & ! -48
      &  1.9338445142_wp,  0.8773559861_wp,  0.3739992316_wp,  0.4517948170_wp, & ! -52
      &  0.3041556301_wp,  0.3290807250_wp,  1.9181133121_wp,  1.7899515460_wp, & ! -56
      &  2.0844850963_wp,  0.9136690116_wp,  0.8484898306_wp,  0.5660217937_wp, & ! -60
      &  0.5064170106_wp,  0.4731908893_wp,  1.1891855140_wp,  0.5071098795_wp, & ! -64
      &  0.5824925420_wp,  0.7716650814_wp,  0.5509526594_wp,  0.3183237219_wp, & ! -68
      &  0.2877717791_wp,  0.5259043531_wp,  0.7234286343_wp,  0.9421873340_wp, & ! -72
      &  1.7311512798_wp,  1.9740750517_wp,  2.2064272775_wp,  0.6060090165_wp, & ! -76
      &  1.3009723793_wp,  0.4726592528_wp,  0.5878545519_wp,  1.4508221664_wp, & ! -80
      &  0.5464038759_wp,  0.8246952678_wp,  0.9178667857_wp,  0.4513286733_wp, & ! -84
      &  0.3927258076_wp,  0.7343363905_wp,  0.9832048840_wp,  1.3019615779_wp, & ! -88
      &  0.9237282914_wp,  0.9516418629_wp,  0.9962947453_wp,  1.3291362100_wp, & ! -92
      &  1.0607380277_wp,  1.7496211778_wp,  2.0864700431_wp,  1.1745934488_wp, & ! -96
      &  1.6410458337_wp,  1.3975001895_wp,  1.4464121202_wp,  1.3956575722_wp, & ! -100
      &  1.3080725818_wp,  1.1455579374_wp,  1.1997772501_wp]

   !> Parameter for mixed charge and CN dependence
   real(wp), parameter :: p_k3(max_elem) = [&
      & -0.1482235383_wp,  0.0797773499_wp,  0.0129050210_wp,  0.8011212947_wp, & ! -4
      & -0.1403501688_wp,  0.0462480322_wp,  0.1060356466_wp,  0.1795566137_wp, & ! -8
      &  0.3893098244_wp, -0.0146687422_wp, -0.3111690904_wp,  0.1469562808_wp, & ! -12
      & -0.1552600500_wp, -0.0517809794_wp,  0.0691802709_wp,  0.0491066807_wp, & ! -16
      & -0.1037356424_wp,  0.0395605098_wp, -0.0557489720_wp,  0.2496747921_wp, & ! -20
      &  0.9474558294_wp,  0.7590523955_wp, -0.2695867206_wp, -0.2787536365_wp, & ! -24
      &  0.0540450968_wp, -0.2086478372_wp, -0.2234162079_wp, -0.1156520950_wp, & ! -28
      &  0.2563246081_wp, -0.3191068261_wp,  0.1488752906_wp,  0.0371208669_wp, & ! -32
      &  0.1743722400_wp, -0.0216390640_wp,  0.1235826506_wp,  0.3145302762_wp, & ! -36
      &  0.0975461822_wp,  0.5961909581_wp,  1.0118652254_wp,  0.0186250109_wp, & ! -40
      &  0.1194269470_wp,  0.0194565795_wp,  0.1861886669_wp, -0.0573939752_wp, & ! -44
      &  0.0834521615_wp,  0.7760917749_wp,  0.1111065737_wp,  0.1647006719_wp, & ! -48
      & -0.0228253273_wp,  0.0412157970_wp,  0.0337923911_wp, -0.0369295411_wp, & ! -52
      &  0.0435379390_wp,  0.1775358173_wp,  0.2804861722_wp,  0.1710145541_wp, & ! -56
      &  1.1089870351_wp,  0.2042327307_wp,  0.1936009387_wp,  0.2300270019_wp, & ! -60
      &  0.1460300571_wp,  0.1738570015_wp,  0.1110552072_wp,  0.1062222399_wp, & ! -64
      &  0.0837731401_wp,  0.1540607852_wp,  0.1449757942_wp,  0.0936359310_wp, & ! -68
      &  0.3750429911_wp,  0.0851732052_wp,  0.0677646919_wp,  0.1532929419_wp, & ! -72
      &  0.0308496571_wp, -0.1713827638_wp,  0.0241358266_wp, -0.0358765818_wp, & ! -76
      &  0.1216437351_wp,  0.8523468378_wp,  0.0755461658_wp,  0.3629449308_wp, & ! -80
      &  0.1584940355_wp,  0.0821136317_wp,  0.0540698456_wp,  0.0614960628_wp, & ! -84
      & -0.0511469514_wp,  0.1750107802_wp, -0.1736992507_wp,  0.1358232940_wp, & ! -88
      &  0.1982890822_wp,  0.2068029681_wp,  0.0910237200_wp, -0.0609937959_wp, & ! -92
      & -0.0552882699_wp,  0.1101861282_wp, -0.0640006225_wp, -0.0888724586_wp, & ! -96
      & -0.0580125827_wp, -0.0426636961_wp, -0.0631687199_wp, -0.0696692533_wp, & ! -100
      & -0.0827656978_wp, -0.0629399715_wp, -0.0776964767_wp]

   !> Empirical atomic radii for calculation of the coordination number (modified Pyykkö radii)
   real(wp), parameter :: qvszp_cov_radii(max_elem) = 1.889725949_wp * [&
      & 0.29_wp, 0.46_wp, 1.20_wp, 0.94_wp, 0.77_wp, 0.75_wp, 0.71_wp, 0.63_wp, & ! 1-8
      & 0.64_wp, 0.67_wp, 1.40_wp, 1.25_wp, 1.13_wp, 1.04_wp, 1.10_wp, 1.02_wp, & ! 9-16 
      & 0.99_wp, 0.96_wp, 1.76_wp, 1.54_wp, 1.33_wp, 1.22_wp, 1.21_wp, 1.10_wp, & ! 17-24
      & 1.07_wp, 1.04_wp, 1.00_wp, 0.99_wp, 1.01_wp, 1.09_wp, 1.12_wp, 1.09_wp, & ! 25-32
      & 1.15_wp, 1.10_wp, 1.14_wp, 1.17_wp, 1.89_wp, 1.67_wp, 1.47_wp, 1.39_wp, & ! 33-40
      & 1.32_wp, 1.24_wp, 1.15_wp, 1.13_wp, 1.13_wp, 1.08_wp, 1.15_wp, 1.23_wp, & ! 41-48
      & 1.28_wp, 1.26_wp, 1.26_wp, 1.23_wp, 1.32_wp, 1.31_wp, 2.09_wp, 1.76_wp, & ! 49-56
      & 1.62_wp, 1.47_wp, 1.58_wp, 1.57_wp, 1.56_wp, 1.55_wp, 1.51_wp, 1.52_wp, & ! 57-64
      & 1.51_wp, 1.50_wp, 1.49_wp, 1.49_wp, 1.48_wp, 1.53_wp, 1.46_wp, 1.37_wp, & ! 65-72
      & 1.31_wp, 1.23_wp, 1.18_wp, 1.16_wp, 1.11_wp, 1.12_wp, 1.13_wp, 1.32_wp, & ! 73-80
      & 1.30_wp, 1.30_wp, 1.36_wp, 1.31_wp, 1.38_wp, 1.42_wp, 2.01_wp, 1.81_wp, & ! 81-88
      & 1.67_wp, 1.58_wp, 1.52_wp, 1.53_wp, 1.54_wp, 1.55_wp, 1.49_wp, 1.49_wp, & ! 89-96
      & 1.51_wp, 1.51_wp, 1.48_wp, 1.50_wp, 1.56_wp, 1.58_wp, 1.45_wp]            ! 97-103

   !> Steepness of counting function for the coordination number
   real(wp), parameter :: kcn = 3.75_wp

   !> CGTO exponents (Exponent of the primitive Gaussian functions)
   real(wp), protected :: exponents(maxg, max_shell, max_elem)

   !> CGTO coefficients (Contraction coefficients of the primitive Gaussian functions,
   !> might contain normalization)
   real(wp), protected :: coefficients(maxg, max_shell, max_elem)

   !> CGTO coefficients (Contraction coefficients of the primitive Gaussian functions,
   !> might contain normalization)
   real(wp), protected :: coefficients_env(maxg, max_shell, max_elem)

   include 'q-vszp.inc'

   !> Regularization parameter for the square-root coordination number
   real(wp), parameter :: cn_reg = 1e-6_wp

contains

!> Create a new basis set
subroutine new_qvszp_basis(self, mol, nshell, cgto, error, accuracy, cgto_h0, &
   & charge_dependent)
   !> Instance of the basis set data
   type(qvszp_basis_type), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Number of shells per species
   integer, intent(in) :: nshell(:)
   !> Contracted Gaussian basis functions for each shell and species
   type(cgto_container), intent(in) :: cgto(:, :)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Calculation accuracy
   real(wp), intent(in), optional :: accuracy
   !> Optional scaled contracted Gaussian basis functions for the H0 construction
   type(cgto_container), intent(in), optional :: cgto_h0(:, :)
   !> Flag to indicate if the basis set is charge dependent
   logical, intent(in), optional :: charge_dependent

   integer :: iat, isp, ish, iao, ii, iicart
   real(wp) :: min_alpha, acc

   self%nsh_id = nshell

   self%cgto = cgto
   if (present(cgto_h0)) then
      self%cgto_h0 = cgto_h0
   end if

   ! Flag to decide if the basis set should be treated as charge-dependent
   if (present(charge_dependent)) then
      self%charge_dependent = charge_dependent
   else
      self%charge_dependent = .true.
   end if

   ! Integral cutoffs depending on the accuracy
   if (present(accuracy)) then
      acc = accuracy
   else
      acc = default_accuracy
   end if
   self%intcut = integral_cutoff(acc)

   ! Make count of shells for each atom
   self%nsh_at = nshell(mol%id)

   ! Create mapping between atoms and shells
   self%nsh = sum(self%nsh_at)
   allocate(self%ish_at(mol%nat), self%sh2at(self%nsh))
   ii = 0
   do iat = 1, mol%nat
      self%ish_at(iat) = ii
      do ish = 1, self%nsh_at(iat)
         self%sh2at(ii+ish) = iat
      end do
      ii = ii + self%nsh_at(iat)
   end do

   ! Make count of spherical orbitals for each shell
   allocate(self%nao_sh(self%nsh), self%nao_cart_sh(self%nsh))
   do iat = 1, mol%nat
      isp = mol%id(iat)
      ii = self%ish_at(iat)
      do ish = 1, self%nsh_at(iat)
         associate(p_cgto => self%cgto(ish, isp)%raw)
            self%nao_sh(ii+ish) = 2*p_cgto%ang + 1
            self%nao_cart_sh(ii+ish) = (p_cgto%ang + 1)*(p_cgto%ang + 2)/2
         end associate
      end do
   end do

   ! Create mapping between shells and spherical or cartesian orbitals,
   ! also map directly back to atoms
   self%nao = sum(self%nao_sh)
   self%nao_cart = sum(self%nao_cart_sh)
   allocate(self%iao_sh(self%nsh), self%iao_cart_sh(self%nsh), &
      & self%ao2sh(self%nao), self%ao2at(self%nao))
   ii = 0
   iicart = 0
   do ish = 1, self%nsh
      self%iao_sh(ish) = ii
      self%iao_cart_sh(ish) = iicart
      do iao = 1, self%nao_sh(ish)
         self%ao2sh(ii+iao) = ish
         self%ao2at(ii+iao) = self%sh2at(ish)
      end do
      ii = ii + self%nao_sh(ish)
      iicart = iicart + self%nao_cart_sh(ish)
   end do

   ii = 0
   do iat = 1, mol%nat
      isp = mol%id(iat)
      do ish = 1, nshell(isp)
         self%iao_sh(ish+self%ish_at(iat)) = ii
         associate(p_cgto => self%cgto(ish, isp)%raw)
            ii = ii + 2*p_cgto%ang + 1
         end associate
      end do
   end do

   min_alpha = huge(acc)
   do isp = 1, size(nshell)
      do ish = 1, nshell(isp)
         associate(p_cgto => self%cgto(ish, isp)%raw)
            self%maxl = max(self%maxl, p_cgto%ang)
            min_alpha = min(min_alpha, minval(p_cgto%alpha(:p_cgto%nprim)))
         end associate
      end do
   end do

   self%min_alpha = min_alpha

   ! Setup coordination number if the basis set is environment dependent
   call new_ncoord(self%ncoord, mol, cn_count_type=cn_count%erf, &
      & error=error, kcn=kcn, rcov=qvszp_cov_radii(mol%num))
   
end subroutine new_qvszp_basis


!> Update basis cache
subroutine basis_update(self, mol, cache, grad, wfn_aux)
   !> Instance of the basis type
   class(qvszp_basis_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different runs
   type(basis_cache), intent(inout) :: cache
   !> Flag to indicate if gradients are required
   logical, intent(in) :: grad
   !> Optional auxiliary wavefunction data
   type(wavefunction_type), intent(in), optional :: wfn_aux

   real(wp), allocatable :: cn(:)
   integer :: max_nsh, ish, iat

   ! Reset all cached quantities
   if (allocated(cache%cgto)) then
      cache%cgto(:, :)%norm = 0.0_wp
      cache%cgto(:, :)%dnorm = 0.0_wp
      cache%cgto(:, :)%qeff = 0.0_wp
   end if
   if (allocated(cache%cgto_h0)) then
      cache%cgto_h0(:, :)%norm = 0.0_wp
      cache%cgto_h0(:, :)%dnorm = 0.0_wp
      cache%cgto_h0(:, :)%qeff = 0.0_wp
   end if

   if (present(wfn_aux)) then 
      ! Calculate the coordination number
      allocate(cn(mol%nat))
      call self%ncoord%get_cn(mol, cn)

      ! Update the CGTO cache with charge/CN dependence
      max_nsh = maxval(self%nsh_at)
      if (.not. allocated(cache%cgto)) allocate(cache%cgto(max_nsh, mol%nat))
      do iat = 1, mol%nat
         do ish = 1, self%nsh_at(iat)
            associate(p_cgto => self%cgto(ish, mol%id(iat))%raw)
               call p_cgto%update(cache%cgto(ish, iat), grad, cn=cn(iat), &
                  & q=wfn_aux%qat(iat, 1))
            end associate
         end do
      end do 
      ! Update the CGTO H0 cache with charge/CN dependence
      if (allocated(self%cgto_h0)) then
         if (.not. allocated(cache%cgto_h0)) allocate(cache%cgto_h0(max_nsh, mol%nat))
         do iat = 1, mol%nat
            do ish = 1, self%nsh_at(iat)
               associate(p_cgto_h0 => self%cgto_h0(ish, mol%id(iat))%raw)
                  call p_cgto_h0%update(cache%cgto_h0(ish, iat), grad, cn=cn(iat), &
                     & q=wfn_aux%qat(iat, 1))
               end associate
            end do
         end do 
      end if 
   else
      ! Update the CGTO cache without charge dependence
      max_nsh = maxval(self%nsh_at)
      if (.not. allocated(cache%cgto)) allocate(cache%cgto(max_nsh, mol%nat))
      do iat = 1, mol%nat
         do ish = 1, self%nsh_at(iat)
            associate(p_cgto => self%cgto(ish, mol%id(iat))%raw)
               call p_cgto%update(cache%cgto(ish, iat), .false.)
            end associate
         end do
      end do 

      ! Update the CGTO H0 cache without charge dependence 
      if (allocated(self%cgto_h0)) then
         if (.not. allocated(cache%cgto_h0)) allocate(cache%cgto_h0(max_nsh, mol%nat))
         do iat = 1, mol%nat
            do ish = 1, self%nsh_at(iat)
               associate(p_cgto_h0 => self%cgto_h0(ish, mol%id(iat))%raw)
                  call p_cgto_h0%update(cache%cgto_h0(ish, iat), .false.)
               end associate
            end do
         end do 
      end if 
   end if

end subroutine basis_update


!> Contract gradients w.r.t. the basis set specific charge and CN
subroutine get_basis_gradient(self, mol, dEdcnbas, dEdqbas, dEdq, gradient, sigma)
   !> Instance of the basis type
   class(qvszp_basis_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Derivative of the energy w.r.t. the basis set coordination number
   real(wp), intent(in), optional :: dEdcnbas(:)
   !> Derivative of the energy w.r.t. the basis set atomic charges
   real(wp), intent(in), optional :: dEdqbas(:)
   !> Derivative of the energy w.r.t. the atomic partial charges from an external model
   real(wp), intent(inout), optional :: dEdq(:)
   !> Derivative of the energy w.r.t. coordinate displacements
   real(wp), intent(inout) :: gradient(:, :)
   !> Derivative of the energy w.r.t. strain deformations
   real(wp), intent(inout) :: sigma(:, :)

   real(wp), allocatable :: lattr(:, :)

   if (present(dEdcnbas)) then
      call get_lattice_points(mol%periodic, mol%lattice, self%ncoord%cutoff, lattr)

      ! Directly contract the coordination number derivatives
      call self%ncoord%add_coordination_number_derivs(mol, lattr, dEdcnbas, &
         & gradient, sigma)
   end if

   if (present(dEdqbas) .and. present(dEdq)) then
      ! q-vSZP basis set charge derivatives add directly to external charge derivatives
      dEdq(:) = dEdq + dEdqbas
   end if

end subroutine get_basis_gradient


!> Setup CGTO for specific atom and shell
subroutine new_qvszp_cgto(self, izp, ish, norm, error, &
   & expos, coeffs, coeffs_env, k0, k1, k2, k3)
   !> Instance of the cgto data
   type(qvszp_cgto_type), intent(out) :: self
   !> Atomic number
   integer, intent(in) :: izp
   !> Shell number
   integer, intent(in) :: ish
   !> Include normalization in contraction coefficients
   logical, intent(in) :: norm
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Optional exponent of Gaussian functions
   real(wp), intent(in), optional :: expos(:)
   !> Optional coefficients of the Gaussian functions
   real(wp), intent(in), optional :: coeffs(:)
   !> Optional environment-dependent coefficients of the Gaussian functions
   real(wp), intent(in), optional :: coeffs_env(:)
   !> Optional overall charge-dependence of the coefficient
   real(wp), intent(in), optional :: k0
   !> Optional squared charge-dependence of the coefficient
   real(wp), intent(in), optional :: k1
   !> Optional CN-dependence of the coefficient
   real(wp), intent(in), optional :: k2
   !> Mixed charge- and CN-dependence of the coefficient
   real(wp), intent(in), optional :: k3

   integer :: il, nprim
   real(wp) :: normalizer(maxg)

   ! Check if the requested shell is available
   if (ish .gt. nshell(izp) .and. (.not. present(expos) .or. .not. present(coeffs))) then 
      call fatal_error(error, "q-vSZP shell not available for element with Z >" // &
         & format_string(izp, '(i0)') // ".")
      return
   end if 

   il = ang_shell(ish, izp)
   nprim = n_prim(ish, izp)

   self%ang = il
   self%nprim = nprim
   
   if (present(expos)) then
      self%alpha(1:nprim) = expos(1:nprim)
   else
      self%alpha(1:nprim) = exponents(1:nprim, ish, izp)
   end if

   normalizer = 1.0_wp
   if (norm) then
      normalizer(:nprim) = (top*self%alpha(:nprim))**0.75_wp &
         & * sqrt(4.0_wp*self%alpha(:nprim))**il / sqrt(dfactorial(il+1))
   endif

   if (present(coeffs)) then
      self%coeff(1:nprim) = coeffs(1:nprim) * normalizer(1:nprim)
   else
      self%coeff(1:nprim) = coefficients(1:nprim, ish, izp) * normalizer(1:nprim)
   end if

   if (present(coeffs_env)) then
      self%coeff1(1:nprim) = coeffs_env(1:nprim) * normalizer(1:nprim)
   else
      self%coeff1(1:nprim) = coefficients_env(1:nprim, ish, izp) * normalizer(1:nprim)
   end if

   if (present(k0)) then
      self%k0 = k0
   else
      self%k0 = p_k0(izp)
   end if

   if (present(k1)) then
      self%k1 = k1
   else
      self%k1 = p_k1(izp)
   end if

   if (present(k2)) then
      self%k2 = k2
   else
      self%k2 = p_k2(izp)
   end if

   if (present(k3)) then
      self%k3 = k3
   else
      self%k3 = p_k3(izp)
   end if

end subroutine new_qvszp_cgto

!> Update CGTO cache
subroutine cgto_update(self, cache, grad, cn, q)
   !> Instance of the basis type
   class(qvszp_cgto_type), intent(in) :: self
   !> Cached data between different runs
   type(cgto_cache), intent(inout) :: cache
   !> Flag to get gradient information
   logical, intent(in) :: grad
   !> Coordination number
   real(wp), intent(in), optional :: cn
   !> Effective charge
   real(wp), intent(in), optional :: q

   ! Get effective charge and its derivatives
   cache%qeff = 0.0_wp
   if (present(cn) .and. present(q)) then
      if (grad) then
         call self%get_qeff(q, cn, cache%qeff, cache%dqeffdcn, cache%dqeffdq)
      else
         call self%get_qeff(q, cn, cache%qeff)
      end if
   end if

   ! Get normalization factor and its derivatives
   call self%get_normalization(cache, grad)

end subroutine cgto_update

!> Get normalization constant for the CGTOs
pure subroutine get_normalization(self, cache, grad)
   !> Instance of the basis set data
   class(qvszp_cgto_type), intent(in) :: self
   !> Cached data between different runs
   type(cgto_cache), intent(inout) :: cache
   !> Flag to get gradient of normalization
   logical, intent(in) :: grad

   integer :: iao
   real(wp) :: r2, vec(3)
   real(wp) :: overlap(msao(max_shell), msao(max_shell))
   real(wp) :: doverlap(3, msao(max_shell), msao(max_shell))
   real(wp) :: doverlapdqeffi(msao(max_shell), msao(max_shell))
   real(wp) :: doverlapdqeffj(msao(max_shell), msao(max_shell))

   r2 = 0.0_wp
   vec = 0.0_wp
   overlap = 0.0_wp
   cache%norm = 1.0_wp
   cache%dnorm = 0.0_wp
   doverlapdqeffi = 0.0_wp
   doverlapdqeffj = 0.0_wp
   if (grad) then
      call overlap_grad_cgto(self, self, cache, cache, r2, vec, 100.0_wp, &
         & overlap, doverlap, doverlapdqeffi, doverlapdqeffj)
      cache%norm = 1.0_wp / sqrt(overlap(1, 1))
      cache%dnorm = -0.5_wp * (doverlapdqeffi(1, 1) + doverlapdqeffj(1, 1)) &
         & * cache%norm**3
   else
      call overlap_cgto(self, self, cache, cache, r2, vec, 100.0_wp, overlap)
      cache%norm = 1.0_wp / sqrt(overlap(1, 1))
   end if

end subroutine get_normalization

!> Get (scaled) coefficients of the CGTOs
pure subroutine get_coeffs(self, cache, coeff)
   !> Instance of the basis set data
   class(qvszp_cgto_type), intent(in) :: self
   !> Cached data between different runs
   type(cgto_cache), intent(in) :: cache
   !> Scaled coefficients of the CGTO
   real(wp), intent(out) :: coeff(:)

   coeff = (self%coeff(1:self%nprim) + self%coeff1(1:self%nprim) * cache%qeff) &
      & * cache%norm

end subroutine get_coeffs

!> Get coefficient derivatives of the CGTOs
pure subroutine get_coeff_derivs(self, cache, dcoeff)
   !> Instance of the basis set data
   class(qvszp_cgto_type), intent(in) :: self
   !> Cached data between different runs
   type(cgto_cache), intent(in) :: cache
   !> Derivative of coefficients of the CGTO w.r.t. effective charge
   real(wp), intent(out) :: dcoeff(:)
   
   dcoeff = self%coeff1(1:self%nprim) * cache%norm &
      & + (self%coeff(1:self%nprim) + self%coeff1(1:self%nprim) * cache%qeff) &
      & * cache%dnorm

end subroutine get_coeff_derivs

!> Update the effective environment-dependent charge
subroutine get_qeff(cgto, qat, cn, qeff, dqeffdcn, dqeffdq)
   !> Instance of the cgto data
   class(qvszp_cgto_type), intent(in) :: cgto
   !> Atom-resolved charge
   real(wp), intent(in) :: qat
   !> Atom-resolved coordination number
   real(wp), intent(in) :: cn
   !> Scaled effective nuclear charge
   real(wp), intent(out) :: qeff
   !> Derivative of effective environment-dependent charge w.r.t. coordination number
   real(wp), intent(out), optional :: dqeffdcn
   !> Derivative of effective environment-dependent charge w.r.t. charge
   real(wp), intent(out), optional :: dqeffdq

   real(wp) :: sqrtcn

   sqrtcn = sqrt(cn + cn_reg**2) - cn_reg 
   qeff = cgto%k0 * (qat - cgto%k1 * qat**2) + cgto%k2 * sqrtcn + cgto%k3 * cn * qat

   if (present(dqeffdcn) .and. present(dqeffdq)) then
      dqeffdcn = 0.5_wp * (cgto%k2 / (sqrtcn + cn_reg)) + cgto%k3 * qat
      dqeffdq = cgto%k0 * (1.0_wp - 2.0_wp * cgto%k1 * qat) + cgto%k3 * cn
   end if

end subroutine get_qeff

!> Return all q-vSZP parameters
subroutine get_qvszp_parameters(nshell, lsh, pqn, ngauss, expos, &
   & coeffs, coeffs_env, k0, k1, k2, k3)
   !> Number of shells per species
   integer, intent(in) :: nshell(:)
   !> Angular momentum of the shells
   integer, intent(out) :: lsh(:, :)
   !> Principal quantum number of the shells
   integer, intent(out) :: pqn(:, :)
   !> Number of primitive Gaussian functions per shell
   integer, intent(out) :: ngauss(:, :)
   !> Exponents of the primitive Gaussian functions
   real(wp), intent(out) :: expos(:, :, :)
   !> Coefficients of the primitive Gaussian functions
   real(wp), intent(out) :: coeffs(:, :, :)
   !> Environment coefficients of the primitive Gaussian functions
   real(wp), intent(out) :: coeffs_env(:, :, :)
   !> Overall charge-dependence of the coefficient
   real(wp), intent(out) :: k0(:)
   !> Squared charge-dependence of the coefficient
   real(wp), intent(out) :: k1(:)
   !> CN-dependence of the coefficient
   real(wp), intent(out) :: k2(:)
   !> Mixed charge- and CN-dependence of the coefficient
   real(wp), intent(out) :: k3(:)

   integer :: izp, ish

   do izp = 1, max_elem
      do ish = 1, nshell(izp)
         pqn(ish, izp) = principal_quantum_number(ish, izp)
         lsh(ish, izp) = ang_shell(ish, izp)
         ngauss(ish, izp) = n_prim(ish, izp)
         expos(:, ish, izp) = exponents(:, ish, izp)
         coeffs(:, ish, izp) = coefficients(:, ish, izp)
         coeffs_env(:, ish, izp) = coefficients_env(:, ish, izp)
      end do
   end do
   k0 = p_k0
   k1 = p_k1
   k2 = p_k2
   k3 = p_k3

end subroutine get_qvszp_parameters

end module tblite_basis_qvszp