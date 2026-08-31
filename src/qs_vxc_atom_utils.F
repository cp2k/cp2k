!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Support routines for integrals of the Vxc/Fxc/Gxc potentials calculated
!>      for the atomic density in the basis set of spherical primitives
! **************************************************************************************************
MODULE qs_vxc_atom_utils
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE external_potential_types,        ONLY: get_potential,&
                                              gth_potential_type,&
                                              sgp_potential_type
   USE kinds,                           ONLY: dp
   USE orbital_pointers,                ONLY: indso,&
                                              nsoset
   USE paw_basis_types,                 ONLY: get_paw_basis_info
   USE qs_grid_atom,                    ONLY: grid_atom_type
   USE qs_harmonics_atom,               ONLY: get_none0_cg_list,&
                                              harmonics_atom_type
   USE qs_kind_types,                   ONLY: has_nlcc
   USE qs_rho_atom_types,               ONLY: rho_atom_coeff,&
                                              rho_atom_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_vxc_atom_utils'

   TYPE tau_basis_cache_type
      INTEGER                                            :: maxso = 0, na = 0, nr = 0, nsatbas = 0, &
                                                            nset = 0
      INTEGER, DIMENSION(:), POINTER                     :: lmax => NULL(), lmin => NULL(), &
                                                            n2oindex => NULL(), npgf => NULL(), &
                                                            o2nindex => NULL()
      REAL(dp), DIMENSION(:, :), POINTER                 :: zet => NULL()
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: grad
   END TYPE tau_basis_cache_type

   PUBLIC :: evaluate_nlcc_primitive_fields, gapw_atom_grid_support_radius, &
             create_tau_basis_cache, release_tau_basis_cache, &
             calc_rho_angular, calc_tau_atom, calc_rho_nlcc, &
             gaVxcgb_noGC, gaVxcgb_GC, dgaVtaudgb
   PUBLIC :: tau_basis_cache_type

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief Evaluate an NLCC density and its first two Cartesian derivatives.
!> \param point evaluation point
!> \param center pseudopotential center
!> \param gth_potential optional GTH potential
!> \param sgp_potential optional separable Gaussian potential
!> \param rho core density
!> \param gradient Cartesian density gradient
!> \param hessian Cartesian density Hessian
! **************************************************************************************************
   SUBROUTINE evaluate_nlcc_primitive_fields(point, center, gth_potential, sgp_potential, &
                                             rho, gradient, hessian)
      REAL(dp), DIMENSION(3), INTENT(IN)                 :: point, center
      TYPE(gth_potential_type), INTENT(IN), POINTER      :: gth_potential
      TYPE(sgp_potential_type), INTENT(IN), POINTER      :: sgp_potential
      REAL(dp), INTENT(OUT)                              :: rho
      REAL(dp), DIMENSION(3), INTENT(OUT)                :: gradient
      REAL(dp), DIMENSION(3, 3), INTENT(OUT)             :: hessian

      INTEGER                                            :: ic, idir, iexp, jdir, n_nlcc, nexp_nlcc, &
                                                            power
      INTEGER, DIMENSION(:), POINTER                     :: nct_nlcc
      LOGICAL                                            :: has_sgp_nlcc, nlcc_present
      REAL(dp)                                           :: alpha, beta, d2poly, dpoly, exponential, &
                                                            poly, r2, rho_x, rho_xx, scaled_r2
      REAL(dp), DIMENSION(3)                             :: displacement
      REAL(dp), DIMENSION(:), POINTER                    :: a_nlcc, alpha_nlcc, c_nlcc
      REAL(dp), DIMENSION(:, :), POINTER                 :: cval_nlcc

      NULLIFY (a_nlcc, alpha_nlcc, c_nlcc, cval_nlcc, nct_nlcc)
      rho = 0.0_dp
      rho_x = 0.0_dp
      rho_xx = 0.0_dp
      displacement = point - center
      r2 = DOT_PRODUCT(displacement, displacement)

      IF (ASSOCIATED(gth_potential)) THEN
         CALL get_potential(gth_potential, nlcc_present=nlcc_present, &
                            nexp_nlcc=nexp_nlcc, alpha_nlcc=alpha_nlcc, &
                            nct_nlcc=nct_nlcc, cval_nlcc=cval_nlcc)
         IF (nlcc_present) THEN
            DO iexp = 1, nexp_nlcc
               alpha = alpha_nlcc(iexp)
               beta = 0.5_dp/(alpha*alpha)
               scaled_r2 = r2/(alpha*alpha)
               exponential = EXP(-0.5_dp*scaled_r2)
               DO ic = 1, nct_nlcc(iexp)
                  power = ic - 1
                  poly = cval_nlcc(ic, iexp)*scaled_r2**power
                  dpoly = 0.0_dp
                  IF (power > 0) THEN
                     dpoly = cval_nlcc(ic, iexp)*REAL(power, dp)* &
                             scaled_r2**(power - 1)/(alpha*alpha)
                  END IF
                  d2poly = 0.0_dp
                  IF (power > 1) THEN
                     d2poly = cval_nlcc(ic, iexp)*REAL(power*(power - 1), dp)* &
                              scaled_r2**(power - 2)/(alpha**4)
                  END IF
                  rho = rho + exponential*poly
                  rho_x = rho_x + exponential*(dpoly - beta*poly)
                  rho_xx = rho_xx + exponential*(d2poly - 2.0_dp*beta*dpoly + beta*beta*poly)
               END DO
            END DO
         END IF
      ELSE IF (ASSOCIATED(sgp_potential)) THEN
         CALL get_potential(sgp_potential, has_nlcc=has_sgp_nlcc, n_nlcc=n_nlcc, &
                            a_nlcc=a_nlcc, c_nlcc=c_nlcc)
         IF (has_sgp_nlcc) THEN
            DO iexp = 1, n_nlcc
               exponential = EXP(-a_nlcc(iexp)*r2)
               rho = rho + c_nlcc(iexp)*exponential
               rho_x = rho_x - a_nlcc(iexp)*c_nlcc(iexp)*exponential
               rho_xx = rho_xx + a_nlcc(iexp)**2*c_nlcc(iexp)*exponential
            END DO
         END IF
      END IF

      gradient = 2.0_dp*rho_x*displacement
      DO idir = 1, 3
         DO jdir = 1, 3
            hessian(idir, jdir) = 4.0_dp*rho_xx*displacement(idir)*displacement(jdir)
         END DO
         hessian(idir, idir) = hessian(idir, idir) + 2.0_dp*rho_x
      END DO
   END SUBROUTINE evaluate_nlcc_primitive_fields

! **************************************************************************************************
!> \brief Compact support radius of one hard-minus-soft atom-grid field.
!> \param grid_atom radial and angular source grid
!> \param rho_h hard one-center density values
!> \param rho_s soft one-center density values
!> \param drho_h hard one-center density-gradient values
!> \param drho_s soft one-center density-gradient values
!> \param tau_h hard one-center kinetic-energy-density values
!> \param tau_s soft one-center kinetic-energy-density values
!> \return outermost radius required by any primitive hard-minus-soft field
! **************************************************************************************************
   FUNCTION gapw_atom_grid_support_radius( &
      grid_atom, rho_h, rho_s, drho_h, drho_s, tau_h, tau_s) RESULT(cutoff)
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      REAL(dp), DIMENSION(:, :, :), INTENT(IN)           :: rho_h, rho_s
      REAL(dp), DIMENSION(:, :, :, :), INTENT(IN)        :: drho_h, drho_s
      REAL(dp), DIMENSION(:, :, :), INTENT(IN)           :: tau_h, tau_s
      REAL(dp)                                           :: cutoff

      INTEGER                                            :: ir, support_ir
      LOGICAL                                            :: descending
      REAL(dp)                                           :: amplitude

      cutoff = 0.0_dp
      descending = grid_atom%rad(1) > grid_atom%rad(grid_atom%nr)
      IF (descending) THEN
         DO ir = 1, grid_atom%nr
            amplitude = MAX( &
                        MAXVAL(ABS(rho_h(:, ir, :) - rho_s(:, ir, :))), &
                        MAXVAL(ABS(drho_h(1:3, :, ir, :) - drho_s(1:3, :, ir, :))), &
                        MAXVAL(ABS(tau_h(:, ir, :) - tau_s(:, ir, :))))
            IF (amplitude > 1.0E-12_dp) THEN
               support_ir = MAX(1, ir - 2)
               cutoff = grid_atom%rad(support_ir)
               RETURN
            END IF
         END DO
      ELSE
         DO ir = grid_atom%nr, 1, -1
            amplitude = MAX( &
                        MAXVAL(ABS(rho_h(:, ir, :) - rho_s(:, ir, :))), &
                        MAXVAL(ABS(drho_h(1:3, :, ir, :) - drho_s(1:3, :, ir, :))), &
                        MAXVAL(ABS(tau_h(:, ir, :) - tau_s(:, ir, :))))
            IF (amplitude > 1.0E-12_dp) THEN
               support_ir = MIN(grid_atom%nr, ir + 2)
               cutoff = grid_atom%rad(support_ir)
               RETURN
            END IF
         END DO
      END IF

   END FUNCTION gapw_atom_grid_support_radius

! **************************************************************************************************
!> \brief ...
!> \param grid_atom ...
!> \param harmonics ...
!> \param nspins ...
!> \param grad_func ...
!> \param ir ...
!> \param r_h ...
!> \param r_s ...
!> \param rho_h ...
!> \param rho_s ...
!> \param dr_h ...
!> \param dr_s ...
!> \param r_h_d ...
!> \param r_s_d ...
!> \param drho_h ...
!> \param drho_s ...
! **************************************************************************************************
   SUBROUTINE calc_rho_angular(grid_atom, harmonics, nspins, grad_func, &
                               ir, r_h, r_s, rho_h, rho_s, &
                               dr_h, dr_s, r_h_d, r_s_d, drho_h, drho_s)

      TYPE(grid_atom_type), POINTER                      :: grid_atom
      TYPE(harmonics_atom_type), POINTER                 :: harmonics
      INTEGER, INTENT(IN)                                :: nspins
      LOGICAL, INTENT(IN)                                :: grad_func
      INTEGER, INTENT(IN)                                :: ir
      TYPE(rho_atom_coeff), DIMENSION(:), POINTER        :: r_h, r_s
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: rho_h, rho_s
      TYPE(rho_atom_coeff), DIMENSION(:), POINTER        :: dr_h, dr_s
      TYPE(rho_atom_coeff), DIMENSION(:, :), POINTER     :: r_h_d, r_s_d
      REAL(KIND=dp), DIMENSION(:, :, :, :), POINTER      :: drho_h, drho_s

      INTEGER                                            :: ia, iso, ispin, na
      REAL(KIND=dp)                                      :: rad, urad

      CPASSERT(ASSOCIATED(r_h))
      CPASSERT(ASSOCIATED(r_s))
      CPASSERT(ASSOCIATED(rho_h))
      CPASSERT(ASSOCIATED(rho_s))
      IF (grad_func) THEN
         CPASSERT(ASSOCIATED(dr_h))
         CPASSERT(ASSOCIATED(dr_s))
         CPASSERT(ASSOCIATED(r_h_d))
         CPASSERT(ASSOCIATED(r_s_d))
         CPASSERT(ASSOCIATED(drho_h))
         CPASSERT(ASSOCIATED(drho_s))
      END IF

      na = grid_atom%ng_sphere
      rad = grid_atom%rad(ir)
      urad = grid_atom%oorad2l(ir, 1)
      DO ispin = 1, nspins
         DO iso = 1, harmonics%max_iso_not0
            DO ia = 1, na
               rho_h(ia, ir, ispin) = rho_h(ia, ir, ispin) + &
                                      r_h(ispin)%r_coef(ir, iso)*harmonics%slm(ia, iso)
               rho_s(ia, ir, ispin) = rho_s(ia, ir, ispin) + &
                                      r_s(ispin)%r_coef(ir, iso)*harmonics%slm(ia, iso)
            END DO ! ia
         END DO ! iso
      END DO ! ispin

      IF (grad_func) THEN
         DO ispin = 1, nspins
            DO iso = 1, harmonics%max_iso_not0
               DO ia = 1, na

                  ! components of the gradient of rho1 hard
                  drho_h(1, ia, ir, ispin) = drho_h(1, ia, ir, ispin) + &
                                             dr_h(ispin)%r_coef(ir, iso)* &
                                             harmonics%a(1, ia)*harmonics%slm(ia, iso) + &
                                             r_h_d(1, ispin)%r_coef(ir, iso)* &
                                             harmonics%slm(ia, iso)

                  drho_h(2, ia, ir, ispin) = drho_h(2, ia, ir, ispin) + &
                                             dr_h(ispin)%r_coef(ir, iso)* &
                                             harmonics%a(2, ia)*harmonics%slm(ia, iso) + &
                                             r_h_d(2, ispin)%r_coef(ir, iso)* &
                                             harmonics%slm(ia, iso)

                  drho_h(3, ia, ir, ispin) = drho_h(3, ia, ir, ispin) + &
                                             dr_h(ispin)%r_coef(ir, iso)* &
                                             harmonics%a(3, ia)*harmonics%slm(ia, iso) + &
                                             r_h_d(3, ispin)%r_coef(ir, iso)* &
                                             harmonics%slm(ia, iso)

                  ! components of the gradient of rho1 soft
                  drho_s(1, ia, ir, ispin) = drho_s(1, ia, ir, ispin) + &
                                             dr_s(ispin)%r_coef(ir, iso)* &
                                             harmonics%a(1, ia)*harmonics%slm(ia, iso) + &
                                             r_s_d(1, ispin)%r_coef(ir, iso)* &
                                             harmonics%slm(ia, iso)

                  drho_s(2, ia, ir, ispin) = drho_s(2, ia, ir, ispin) + &
                                             dr_s(ispin)%r_coef(ir, iso)* &
                                             harmonics%a(2, ia)*harmonics%slm(ia, iso) + &
                                             r_s_d(2, ispin)%r_coef(ir, iso)* &
                                             harmonics%slm(ia, iso)

                  drho_s(3, ia, ir, ispin) = drho_s(3, ia, ir, ispin) + &
                                             dr_s(ispin)%r_coef(ir, iso)* &
                                             harmonics%a(3, ia)*harmonics%slm(ia, iso) + &
                                             r_s_d(3, ispin)%r_coef(ir, iso)* &
                                             harmonics%slm(ia, iso)

               END DO ! ia
            END DO ! iso
            DO ia = 1, na
               drho_h(4, ia, ir, ispin) = SQRT( &
                                          drho_h(1, ia, ir, ispin)*drho_h(1, ia, ir, ispin) + &
                                          drho_h(2, ia, ir, ispin)*drho_h(2, ia, ir, ispin) + &
                                          drho_h(3, ia, ir, ispin)*drho_h(3, ia, ir, ispin))

               drho_s(4, ia, ir, ispin) = SQRT( &
                                          drho_s(1, ia, ir, ispin)*drho_s(1, ia, ir, ispin) + &
                                          drho_s(2, ia, ir, ispin)*drho_s(2, ia, ir, ispin) + &
                                          drho_s(3, ia, ir, ispin)*drho_s(3, ia, ir, ispin))
            END DO ! ia
         END DO ! ispin
      END IF

   END SUBROUTINE calc_rho_angular

! **************************************************************************************************
!> \brief Precompute radial and angular factors for GAPW meta-GGA tau contractions
!> \param tau_cache precomputed compact one-center gradient basis
!> \param grid_atom atom-centered integration grid
!> \param basis_1c GAPW one-center basis
!> \param harmonics spherical harmonics on the atom-centered grid
! **************************************************************************************************
   SUBROUTINE create_tau_basis_cache(tau_cache, grid_atom, basis_1c, harmonics)

      TYPE(tau_basis_cache_type), INTENT(INOUT)          :: tau_cache
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      TYPE(gto_basis_set_type), POINTER                  :: basis_1c
      TYPE(harmonics_atom_type), POINTER                 :: harmonics

      INTEGER                                            :: dir, ia, igrid, ip, ipgf, ir, iset, iso, &
                                                            l, starti
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: a1, a2, gexp, r1, r2
      REAL(dp), DIMENSION(:, :), POINTER                 :: slm
      REAL(dp), DIMENSION(:, :, :), POINTER              :: dslm_dxyz

      NULLIFY (slm, dslm_dxyz)

      CALL release_tau_basis_cache(tau_cache)

      CALL get_gto_basis_set(gto_basis_set=basis_1c, lmax=tau_cache%lmax, &
                             lmin=tau_cache%lmin, maxso=tau_cache%maxso, &
                             npgf=tau_cache%npgf, nset=tau_cache%nset, &
                             zet=tau_cache%zet)
      CALL get_paw_basis_info(basis_1c, o2nindex=tau_cache%o2nindex, &
                              n2oindex=tau_cache%n2oindex, &
                              nsatbas=tau_cache%nsatbas)

      tau_cache%nr = grid_atom%nr
      tau_cache%na = grid_atom%ng_sphere
      slm => harmonics%slm
      dslm_dxyz => harmonics%dslm_dxyz

      ALLOCATE (tau_cache%grad(tau_cache%na*tau_cache%nr, tau_cache%nsatbas, 3))
      ALLOCATE (a1(tau_cache%na), a2(tau_cache%na), gexp(tau_cache%nr), &
                r1(tau_cache%nr), r2(tau_cache%nr))
      tau_cache%grad = 0.0_dp

      DO iset = 1, tau_cache%nset
         DO ipgf = 1, tau_cache%npgf(iset)
            starti = (iset - 1)*tau_cache%maxso + &
                     (ipgf - 1)*nsoset(tau_cache%lmax(iset))
            gexp(1:tau_cache%nr) = EXP(-tau_cache%zet(ipgf, iset)* &
                                       grid_atom%rad2(1:tau_cache%nr))
            DO iso = nsoset(tau_cache%lmin(iset) - 1) + 1, nsoset(tau_cache%lmax(iset))
               ip = tau_cache%o2nindex(starti + iso)
               IF (ip == 0) CYCLE
               l = indso(1, iso)

               r1(1:tau_cache%nr) = grid_atom%rad(1:tau_cache%nr)**(l - 1)*gexp(1:tau_cache%nr)
               r2(1:tau_cache%nr) = -2.0_dp*tau_cache%zet(ipgf, iset)* &
                                    grid_atom%rad2(1:tau_cache%nr)*r1(1:tau_cache%nr)

               DO dir = 1, 3
                  a1(1:tau_cache%na) = dslm_dxyz(dir, 1:tau_cache%na, iso)
                  a2(1:tau_cache%na) = harmonics%a(dir, 1:tau_cache%na)*slm(1:tau_cache%na, iso)
                  DO ir = 1, tau_cache%nr
                     DO ia = 1, tau_cache%na
                        igrid = ia + (ir - 1)*tau_cache%na
                        tau_cache%grad(igrid, ip, dir) = r1(ir)*a1(ia) + r2(ir)*a2(ia)
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      DEALLOCATE (a1, a2, gexp, r1, r2)

   END SUBROUTINE create_tau_basis_cache

! **************************************************************************************************
!> \brief Release precomputed GAPW meta-GGA tau factors
!> \param tau_cache precomputed compact one-center gradient basis
! **************************************************************************************************
   SUBROUTINE release_tau_basis_cache(tau_cache)

      TYPE(tau_basis_cache_type), INTENT(INOUT)          :: tau_cache

      IF (ALLOCATED(tau_cache%grad)) DEALLOCATE (tau_cache%grad)
      IF (ASSOCIATED(tau_cache%n2oindex)) DEALLOCATE (tau_cache%n2oindex)
      IF (ASSOCIATED(tau_cache%o2nindex)) DEALLOCATE (tau_cache%o2nindex)
      NULLIFY (tau_cache%lmax, tau_cache%lmin, tau_cache%n2oindex, tau_cache%npgf, &
               tau_cache%zet, tau_cache%o2nindex)
      tau_cache%maxso = 0
      tau_cache%na = 0
      tau_cache%nr = 0
      tau_cache%nsatbas = 0
      tau_cache%nset = 0

   END SUBROUTINE release_tau_basis_cache

! **************************************************************************************************
!> \brief Computes tau hard and soft on the atomic grids for meta-GGA calculations
!> \param tau_h the hard part of tau
!> \param tau_s the soft part of tau
!> \param rho_atom atom-centered density matrices
!> \param tau_cache precomputed compact one-center gradient basis
!> \param nspins number of spin channels
!> \note This is a rewrite to correct a meta-GGA GAPW bug. This is more brute force than the original,
!>       which was done along in qs_rho_atom_methods.F, but makes sure that no corner is cut in
!>       terms of accuracy (A. Bussy)
! **************************************************************************************************
   SUBROUTINE calc_tau_atom(tau_h, tau_s, rho_atom, tau_cache, nspins)

      REAL(dp), DIMENSION(:, :, :), INTENT(INOUT)        :: tau_h, tau_s
      TYPE(rho_atom_type), POINTER                       :: rho_atom
      TYPE(tau_basis_cache_type), INTENT(IN)             :: tau_cache
      INTEGER, INTENT(IN)                                :: nspins

      CHARACTER(len=*), PARAMETER                        :: routineN = 'calc_tau_atom'

      INTEGER                                            :: dir, handle, ia, ibas, igrid, ir, ispin, &
                                                            na, nbas, ngrid, nr
      REAL(dp)                                           :: tau_sum
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: work

      CALL timeset(routineN, handle)

      CPASSERT(ALLOCATED(tau_cache%grad))

      !zeroing tau, assuming it is already allocated
      tau_h = 0.0_dp
      tau_s = 0.0_dp

      nr = tau_cache%nr
      na = tau_cache%na
      nbas = tau_cache%nsatbas
      ngrid = na*nr
      ALLOCATE (work(ngrid, nbas))

      DO ispin = 1, nspins
         DO dir = 1, 3
            CALL dgemm('N', 'T', ngrid, nbas, nbas, 0.5_dp, tau_cache%grad(:, :, dir), &
                       ngrid, rho_atom%cpc_h(ispin)%r_coef, nbas, 0.0_dp, work, ngrid)
!$OMP PARALLEL DO DEFAULT(NONE) COLLAPSE(2) SCHEDULE(STATIC) &
!$OMP SHARED(dir, ispin, na, nbas, nr, tau_cache, tau_h, work) &
!$OMP PRIVATE(ia, ibas, igrid, ir, tau_sum)
            DO ir = 1, nr
               DO ia = 1, na
                  igrid = ia + (ir - 1)*na
                  tau_sum = 0.0_dp
                  DO ibas = 1, nbas
                     tau_sum = tau_sum + tau_cache%grad(igrid, ibas, dir)*work(igrid, ibas)
                  END DO
                  tau_h(ia, ir, ispin) = tau_h(ia, ir, ispin) + tau_sum
               END DO
            END DO
!$OMP END PARALLEL DO

            CALL dgemm('N', 'T', ngrid, nbas, nbas, 0.5_dp, tau_cache%grad(:, :, dir), &
                       ngrid, rho_atom%cpc_s(ispin)%r_coef, nbas, 0.0_dp, work, ngrid)
!$OMP PARALLEL DO DEFAULT(NONE) COLLAPSE(2) SCHEDULE(STATIC) &
!$OMP SHARED(dir, ispin, na, nbas, nr, tau_cache, tau_s, work) &
!$OMP PRIVATE(ia, ibas, igrid, ir, tau_sum)
            DO ir = 1, nr
               DO ia = 1, na
                  igrid = ia + (ir - 1)*na
                  tau_sum = 0.0_dp
                  DO ibas = 1, nbas
                     tau_sum = tau_sum + tau_cache%grad(igrid, ibas, dir)*work(igrid, ibas)
                  END DO
                  tau_s(ia, ir, ispin) = tau_s(ia, ir, ispin) + tau_sum
               END DO
            END DO
!$OMP END PARALLEL DO
         END DO
      END DO

      DEALLOCATE (work)

      CALL timestop(handle)

   END SUBROUTINE calc_tau_atom

! **************************************************************************************************
!> \brief ...
!> \param grid_atom ...
!> \param nspins ...
!> \param grad_func ...
!> \param ir ...
!> \param rho_nlcc ...
!> \param rho_h ...
!> \param rho_s ...
!> \param drho_nlcc ...
!> \param drho_h ...
!> \param drho_s ...
! **************************************************************************************************
   SUBROUTINE calc_rho_nlcc(grid_atom, nspins, grad_func, &
                            ir, rho_nlcc, rho_h, rho_s, drho_nlcc, drho_h, drho_s)

      TYPE(grid_atom_type), POINTER                      :: grid_atom
      INTEGER, INTENT(IN)                                :: nspins
      LOGICAL, INTENT(IN)                                :: grad_func
      INTEGER, INTENT(IN)                                :: ir
      REAL(KIND=dp), DIMENSION(:)                        :: rho_nlcc
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: rho_h, rho_s
      REAL(KIND=dp), DIMENSION(:)                        :: drho_nlcc
      REAL(KIND=dp), DIMENSION(:, :, :, :), POINTER      :: drho_h, drho_s

      INTEGER                                            :: ia, ispin, na
      REAL(KIND=dp)                                      :: drho, dx, dy, dz, rad, rho, urad, xsp

      CPASSERT(ASSOCIATED(rho_h))
      CPASSERT(ASSOCIATED(rho_s))
      IF (grad_func) THEN
         CPASSERT(ASSOCIATED(drho_h))
         CPASSERT(ASSOCIATED(drho_s))
      END IF

      na = grid_atom%ng_sphere
      rad = grid_atom%rad(ir)
      urad = grid_atom%oorad2l(ir, 1)

      xsp = REAL(nspins, KIND=dp)
      rho = rho_nlcc(ir)/xsp
      DO ispin = 1, nspins
         rho_h(1:na, ir, ispin) = rho_h(1:na, ir, ispin) + rho
         rho_s(1:na, ir, ispin) = rho_s(1:na, ir, ispin) + rho
      END DO ! ispin

      IF (grad_func) THEN
         drho = drho_nlcc(ir)/xsp
         DO ispin = 1, nspins
            DO ia = 1, na
               IF (grid_atom%azi(ia) == 0.0_dp) THEN
                  dx = 0.0_dp
                  dy = 0.0_dp
               ELSE
                  dx = grid_atom%sin_pol(ia)*grid_atom%sin_azi(ia)
                  dy = grid_atom%sin_pol(ia)*grid_atom%cos_azi(ia)
               END IF
               dz = grid_atom%cos_pol(ia)
               ! components of the gradient of rho1 hard
               drho_h(1, ia, ir, ispin) = drho_h(1, ia, ir, ispin) + drho*dx
               drho_h(2, ia, ir, ispin) = drho_h(2, ia, ir, ispin) + drho*dy
               drho_h(3, ia, ir, ispin) = drho_h(3, ia, ir, ispin) + drho*dz
               ! components of the gradient of rho1 soft
               drho_s(1, ia, ir, ispin) = drho_s(1, ia, ir, ispin) + drho*dx
               drho_s(2, ia, ir, ispin) = drho_s(2, ia, ir, ispin) + drho*dy
               drho_s(3, ia, ir, ispin) = drho_s(3, ia, ir, ispin) + drho*dz
               ! norm of gradient
               drho_h(4, ia, ir, ispin) = SQRT( &
                                          drho_h(1, ia, ir, ispin)*drho_h(1, ia, ir, ispin) + &
                                          drho_h(2, ia, ir, ispin)*drho_h(2, ia, ir, ispin) + &
                                          drho_h(3, ia, ir, ispin)*drho_h(3, ia, ir, ispin))

               drho_s(4, ia, ir, ispin) = SQRT( &
                                          drho_s(1, ia, ir, ispin)*drho_s(1, ia, ir, ispin) + &
                                          drho_s(2, ia, ir, ispin)*drho_s(2, ia, ir, ispin) + &
                                          drho_s(3, ia, ir, ispin)*drho_s(3, ia, ir, ispin))
            END DO ! ia
         END DO ! ispin
      END IF

   END SUBROUTINE calc_rho_nlcc

! **************************************************************************************************
!> \brief ...
!> \param vxc_h ...
!> \param vxc_s ...
!> \param int_hh ...
!> \param int_ss ...
!> \param grid_atom ...
!> \param basis_1c ...
!> \param harmonics ...
!> \param nspins ...
! **************************************************************************************************
   SUBROUTINE gaVxcgb_noGC(vxc_h, vxc_s, int_hh, int_ss, grid_atom, basis_1c, harmonics, nspins)

      REAL(dp), DIMENSION(:, :, :), POINTER              :: vxc_h, vxc_s
      TYPE(rho_atom_coeff), DIMENSION(:), POINTER        :: int_hh, int_ss
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      TYPE(gto_basis_set_type), POINTER                  :: basis_1c
      TYPE(harmonics_atom_type), POINTER                 :: harmonics
      INTEGER, INTENT(IN)                                :: nspins

      CHARACTER(len=*), PARAMETER                        :: routineN = 'gaVxcgb_noGC'

      INTEGER :: handle, ia, ic, icg, ipgf1, ipgf2, ir, iset1, iset2, iso, iso1, iso2, ispin, l, &
         ld, lmax12, lmax_expansion, lmin12, m1, m2, max_iso_not0, max_iso_not0_local, max_s_harm, &
         maxl, maxso, n1, n2, na, ngau1, ngau2, nngau1, nr, nset, size1
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: cg_n_list
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: cg_list
      INTEGER, DIMENSION(:), POINTER                     :: lmax, lmin, npgf
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: g1, g2
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: gg, gVg_h, gVg_s, matso_h, matso_s, vx
      REAL(dp), DIMENSION(:, :), POINTER                 :: zet
      REAL(dp), DIMENSION(:, :, :), POINTER              :: my_CG

      CALL timeset(routineN, handle)

      NULLIFY (lmin, lmax, npgf, zet, my_CG)

      CALL get_gto_basis_set(gto_basis_set=basis_1c, lmax=lmax, lmin=lmin, &
                             maxso=maxso, maxl=maxl, npgf=npgf, &
                             nset=nset, zet=zet)

      nr = grid_atom%nr
      na = grid_atom%ng_sphere
      my_CG => harmonics%my_CG
      max_iso_not0 = harmonics%max_iso_not0
      lmax_expansion = indso(1, max_iso_not0)
      max_s_harm = harmonics%max_s_harm

      ALLOCATE (g1(nr), g2(nr), gg(nr, 0:2*maxl))
      ALLOCATE (gVg_h(na, 0:2*maxl), gVg_s(na, 0:2*maxl))
      ALLOCATE (matso_h(nsoset(maxl), nsoset(maxl)), &
                matso_s(nsoset(maxl), nsoset(maxl)))
      ALLOCATE (vx(na, nr))
      ALLOCATE (cg_list(2, nsoset(maxl)**2, max_s_harm), cg_n_list(max_s_harm))

      g1 = 0.0_dp
      g2 = 0.0_dp
      m1 = 0
      DO iset1 = 1, nset
         n1 = nsoset(lmax(iset1))
         m2 = 0
         DO iset2 = 1, nset
            CALL get_none0_cg_list(my_CG, lmin(iset1), lmax(iset1), lmin(iset2), lmax(iset2), &
                                   max_s_harm, lmax_expansion, cg_list, cg_n_list, max_iso_not0_local)
            CPASSERT(max_iso_not0_local <= max_iso_not0)

            n2 = nsoset(lmax(iset2))
            DO ipgf1 = 1, npgf(iset1)
               ngau1 = n1*(ipgf1 - 1) + m1
               size1 = nsoset(lmax(iset1)) - nsoset(lmin(iset1) - 1)
               nngau1 = nsoset(lmin(iset1) - 1) + ngau1

               g1(1:nr) = EXP(-zet(ipgf1, iset1)*grid_atom%rad2(1:nr))
               DO ipgf2 = 1, npgf(iset2)
                  ngau2 = n2*(ipgf2 - 1) + m2

                  g2(1:nr) = EXP(-zet(ipgf2, iset2)*grid_atom%rad2(1:nr))
                  lmin12 = lmin(iset1) + lmin(iset2)
                  lmax12 = lmax(iset1) + lmax(iset2)

                  ! reduce expansion local densities
                  IF (lmin12 <= lmax_expansion) THEN

                     gg = 0.0_dp
                     IF (lmin12 == 0) THEN
                        gg(1:nr, lmin12) = g1(1:nr)*g2(1:nr)
                     ELSE
                        gg(1:nr, lmin12) = grid_atom%rad2l(1:nr, lmin12)*g1(1:nr)*g2(1:nr)
                     END IF

                     ! limit the expansion of the local densities to a max L
                     IF (lmax12 > lmax_expansion) lmax12 = lmax_expansion

                     DO l = lmin12 + 1, lmax12
                        gg(1:nr, l) = grid_atom%rad(1:nr)*gg(:, l - 1)
                     END DO

                     DO ispin = 1, nspins
                        ld = lmax12 + 1
                        DO ir = 1, nr
                           vx(1:na, ir) = vxc_h(1:na, ir, ispin)
                        END DO
                        CALL dgemm('N', 'N', na, ld, nr, 1.0_dp, vx(1:na, 1:nr), na, &
                                   gg(1:nr, 0:lmax12), nr, 0.0_dp, gVg_h(1:na, 0:lmax12), na)
                        DO ir = 1, nr
                           vx(1:na, ir) = vxc_s(1:na, ir, ispin)
                        END DO
                        CALL dgemm('N', 'N', na, ld, nr, 1.0_dp, vx(1:na, 1:nr), na, &
                                   gg(1:nr, 0:lmax12), nr, 0.0_dp, gVg_s(1:na, 0:lmax12), na)

                        matso_h = 0.0_dp
                        matso_s = 0.0_dp
                        DO iso = 1, max_iso_not0_local
                           DO icg = 1, cg_n_list(iso)
                              iso1 = cg_list(1, icg, iso)
                              iso2 = cg_list(2, icg, iso)
                              l = indso(1, iso1) + indso(1, iso2)

                              CPASSERT(l <= lmax_expansion)
                              DO ia = 1, na
                                 matso_h(iso1, iso2) = matso_h(iso1, iso2) + &
                                                       gVg_h(ia, l)* &
                                                       my_CG(iso1, iso2, iso)* &
                                                       harmonics%slm(ia, iso)
                                 matso_s(iso1, iso2) = matso_s(iso1, iso2) + &
                                                       gVg_s(ia, l)* &
                                                       my_CG(iso1, iso2, iso)* &
                                                       harmonics%slm(ia, iso)
                              END DO
                           END DO
                        END DO

                        ! Write in the global matrix
                        DO ic = nsoset(lmin(iset2) - 1) + 1, nsoset(lmax(iset2))
                           iso1 = nsoset(lmin(iset1) - 1) + 1
                           iso2 = ngau2 + ic
                           CALL daxpy(size1, 1.0_dp, matso_h(iso1, ic), 1, &
                                      int_hh(ispin)%r_coef(nngau1 + 1, iso2), 1)
                           CALL daxpy(size1, 1.0_dp, matso_s(iso1, ic), 1, &
                                      int_ss(ispin)%r_coef(nngau1 + 1, iso2), 1)
                        END DO

                     END DO ! ispin

                  END IF ! lmax_expansion

               END DO ! ipfg2
            END DO ! ipfg1
            m2 = m2 + maxso
         END DO ! iset2
         m1 = m1 + maxso
      END DO ! iset1

      DEALLOCATE (g1, g2, gg, matso_h, matso_s, gVg_s, gVg_h, vx)

      DEALLOCATE (cg_list, cg_n_list)

      CALL timestop(handle)

   END SUBROUTINE gaVxcgb_noGC

! **************************************************************************************************
!> \brief ...
!> \param vxc_h ...
!> \param vxc_s ...
!> \param vxg_h ...
!> \param vxg_s ...
!> \param int_hh ...
!> \param int_ss ...
!> \param grid_atom ...
!> \param basis_1c ...
!> \param harmonics ...
!> \param nspins ...
! **************************************************************************************************
   SUBROUTINE gaVxcgb_GC(vxc_h, vxc_s, vxg_h, vxg_s, int_hh, int_ss, &
                         grid_atom, basis_1c, harmonics, nspins)

      REAL(dp), DIMENSION(:, :, :), POINTER              :: vxc_h, vxc_s
      REAL(dp), DIMENSION(:, :, :, :), POINTER           :: vxg_h, vxg_s
      TYPE(rho_atom_coeff), DIMENSION(:), POINTER        :: int_hh, int_ss
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      TYPE(gto_basis_set_type), POINTER                  :: basis_1c
      TYPE(harmonics_atom_type), POINTER                 :: harmonics
      INTEGER, INTENT(IN)                                :: nspins

      CHARACTER(len=*), PARAMETER                        :: routineN = 'gaVxcgb_GC'

      INTEGER :: dmax_iso_not0_local, handle, ia, ic, icg, ipgf1, ipgf2, ir, iset1, iset2, iso, &
         iso1, iso2, ispin, l, lmax12, lmax_expansion, lmin12, m1, m2, max_iso_not0, &
         max_iso_not0_local, max_s_harm, maxl, maxso, n1, n2, na, ngau1, ngau2, nngau1, nr, nset, &
         size1
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: cg_n_list, dcg_n_list
      INTEGER, ALLOCATABLE, DIMENSION(:, :, :)           :: cg_list, dcg_list
      INTEGER, DIMENSION(:), POINTER                     :: lmax, lmin, npgf
      REAL(dp)                                           :: urad
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: g1, g2
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: dgg, gg, gVXCg_h, gVXCg_s, matso_h, &
                                                            matso_s
      REAL(dp), ALLOCATABLE, DIMENSION(:, :, :)          :: gVXGg_h, gVXGg_s
      REAL(dp), DIMENSION(:, :), POINTER                 :: zet
      REAL(dp), DIMENSION(:, :, :), POINTER              :: my_CG
      REAL(dp), DIMENSION(:, :, :, :), POINTER           :: my_CG_dxyz

      CALL timeset(routineN, handle)

      NULLIFY (lmin, lmax, npgf, zet, my_CG, my_CG_dxyz)

      CALL get_gto_basis_set(gto_basis_set=basis_1c, lmax=lmax, lmin=lmin, &
                             maxso=maxso, maxl=maxl, npgf=npgf, &
                             nset=nset, zet=zet)

      nr = grid_atom%nr
      na = grid_atom%ng_sphere
      my_CG => harmonics%my_CG
      my_CG_dxyz => harmonics%my_CG_dxyz
      max_iso_not0 = harmonics%max_iso_not0
      lmax_expansion = indso(1, max_iso_not0)
      max_s_harm = harmonics%max_s_harm

      ALLOCATE (g1(nr), g2(nr), gg(nr, 0:2*maxl), dgg(nr, 0:2*maxl))
      ALLOCATE (gVXCg_h(na, 0:2*maxl), gVXCg_s(na, 0:2*maxl))
      ALLOCATE (gVXGg_h(3, na, 0:2*maxl), gVXGg_s(3, na, 0:2*maxl))
      ALLOCATE (cg_list(2, nsoset(maxl)**2, max_s_harm), cg_n_list(max_s_harm), &
                dcg_list(2, nsoset(maxl)**2, max_s_harm), dcg_n_list(max_s_harm))

      ALLOCATE (matso_h(nsoset(maxl), nsoset(maxl)), &
                matso_s(nsoset(maxl), nsoset(maxl)))

      DO ispin = 1, nspins

         g1 = 0.0_dp
         g2 = 0.0_dp
         m1 = 0
         DO iset1 = 1, nset
            n1 = nsoset(lmax(iset1))
            m2 = 0
            DO iset2 = 1, nset
               CALL get_none0_cg_list(my_CG, lmin(iset1), lmax(iset1), lmin(iset2), lmax(iset2), &
                                      max_s_harm, lmax_expansion, cg_list, cg_n_list, max_iso_not0_local)
               CPASSERT(max_iso_not0_local <= max_iso_not0)
               CALL get_none0_cg_list(my_CG_dxyz, lmin(iset1), lmax(iset1), lmin(iset2), lmax(iset2), &
                                      max_s_harm, lmax_expansion, dcg_list, dcg_n_list, dmax_iso_not0_local)

               n2 = nsoset(lmax(iset2))
               DO ipgf1 = 1, npgf(iset1)
                  ngau1 = n1*(ipgf1 - 1) + m1
                  size1 = nsoset(lmax(iset1)) - nsoset(lmin(iset1) - 1)
                  nngau1 = nsoset(lmin(iset1) - 1) + ngau1

                  g1(1:nr) = EXP(-zet(ipgf1, iset1)*grid_atom%rad2(1:nr))
                  DO ipgf2 = 1, npgf(iset2)
                     ngau2 = n2*(ipgf2 - 1) + m2

                     g2(1:nr) = EXP(-zet(ipgf2, iset2)*grid_atom%rad2(1:nr))
                     lmin12 = lmin(iset1) + lmin(iset2)
                     lmax12 = lmax(iset1) + lmax(iset2)

                     !test reduce expansion local densities
                     IF (lmin12 <= lmax_expansion) THEN

                        gg = 0.0_dp
                        dgg = 0.0_dp

                        IF (lmin12 == 0) THEN
                           gg(1:nr, lmin12) = g1(1:nr)*g2(1:nr)
                        ELSE
                           gg(1:nr, lmin12) = grid_atom%rad2l(1:nr, lmin12)*g1(1:nr)*g2(1:nr)
                        END IF

                        !test reduce expansion local densities
                        IF (lmax12 > lmax_expansion) lmax12 = lmax_expansion

                        DO l = lmin12 + 1, lmax12
                           gg(1:nr, l) = grid_atom%rad(1:nr)*gg(:, l - 1)
                           dgg(1:nr, l - 1) = dgg(1:nr, l - 1) - 2.0_dp*(zet(ipgf1, iset1) + &
                                                                         zet(ipgf2, iset2))*gg(1:nr, l)
                        END DO
                        dgg(1:nr, lmax12) = dgg(1:nr, lmax12) - 2.0_dp*(zet(ipgf1, iset1) + &
                                                                        zet(ipgf2, iset2))*grid_atom%rad(1:nr)* &
                                            gg(1:nr, lmax12)

                        gVXCg_h = 0.0_dp
                        gVXCg_s = 0.0_dp
                        gVXGg_h = 0.0_dp
                        gVXGg_s = 0.0_dp

                        ! Cross Term
                        DO l = lmin12, lmax12
                           DO ia = 1, na
                              DO ir = 1, nr
                                 gVXCg_h(ia, l) = gVXCg_h(ia, l) + &
                                                  gg(ir, l)*vxc_h(ia, ir, ispin) + &
                                                  dgg(ir, l)* &
                                                  (vxg_h(1, ia, ir, ispin)*harmonics%a(1, ia) + &
                                                   vxg_h(2, ia, ir, ispin)*harmonics%a(2, ia) + &
                                                   vxg_h(3, ia, ir, ispin)*harmonics%a(3, ia))

                                 gVXCg_s(ia, l) = gVXCg_s(ia, l) + &
                                                  gg(ir, l)*vxc_s(ia, ir, ispin) + &
                                                  dgg(ir, l)* &
                                                  (vxg_s(1, ia, ir, ispin)*harmonics%a(1, ia) + &
                                                   vxg_s(2, ia, ir, ispin)*harmonics%a(2, ia) + &
                                                   vxg_s(3, ia, ir, ispin)*harmonics%a(3, ia))

                                 urad = grid_atom%oorad2l(ir, 1)

                                 gVXGg_h(1, ia, l) = gVXGg_h(1, ia, l) + &
                                                     vxg_h(1, ia, ir, ispin)* &
                                                     gg(ir, l)*urad

                                 gVXGg_h(2, ia, l) = gVXGg_h(2, ia, l) + &
                                                     vxg_h(2, ia, ir, ispin)* &
                                                     gg(ir, l)*urad

                                 gVXGg_h(3, ia, l) = gVXGg_h(3, ia, l) + &
                                                     vxg_h(3, ia, ir, ispin)* &
                                                     gg(ir, l)*urad

                                 gVXGg_s(1, ia, l) = gVXGg_s(1, ia, l) + &
                                                     vxg_s(1, ia, ir, ispin)* &
                                                     gg(ir, l)*urad

                                 gVXGg_s(2, ia, l) = gVXGg_s(2, ia, l) + &
                                                     vxg_s(2, ia, ir, ispin)* &
                                                     gg(ir, l)*urad

                                 gVXGg_s(3, ia, l) = gVXGg_s(3, ia, l) + &
                                                     vxg_s(3, ia, ir, ispin)* &
                                                     gg(ir, l)*urad

                              END DO ! ir
                           END DO ! ia
                        END DO ! l

                        matso_h = 0.0_dp
                        matso_s = 0.0_dp
                        DO iso = 1, max_iso_not0_local
                           DO icg = 1, cg_n_list(iso)
                              iso1 = cg_list(1, icg, iso)
                              iso2 = cg_list(2, icg, iso)

                              l = indso(1, iso1) + indso(1, iso2)

                              !test reduce expansion local densities
                              CPASSERT(l <= lmax_expansion)
                              DO ia = 1, na
                                 matso_h(iso1, iso2) = matso_h(iso1, iso2) + &
                                                       gVXCg_h(ia, l)* &
                                                       harmonics%slm(ia, iso)* &
                                                       my_CG(iso1, iso2, iso)
                                 matso_s(iso1, iso2) = matso_s(iso1, iso2) + &
                                                       gVXCg_s(ia, l)* &
                                                       harmonics%slm(ia, iso)* &
                                                       my_CG(iso1, iso2, iso)
                              END DO ! ia

                              !test reduce expansion local densities

                           END DO

                        END DO ! iso

                        DO iso = 1, dmax_iso_not0_local
                           DO icg = 1, dcg_n_list(iso)
                              iso1 = dcg_list(1, icg, iso)
                              iso2 = dcg_list(2, icg, iso)

                              l = indso(1, iso1) + indso(1, iso2)
                              !test reduce expansion local densities
                              CPASSERT(l <= lmax_expansion)
                              DO ia = 1, na
                                 matso_h(iso1, iso2) = matso_h(iso1, iso2) + &
                                                       (gVXGg_h(1, ia, l)*my_CG_dxyz(1, iso1, iso2, iso) + &
                                                        gVXGg_h(2, ia, l)*my_CG_dxyz(2, iso1, iso2, iso) + &
                                                        gVXGg_h(3, ia, l)*my_CG_dxyz(3, iso1, iso2, iso))* &
                                                       harmonics%slm(ia, iso)

                                 matso_s(iso1, iso2) = matso_s(iso1, iso2) + &
                                                       (gVXGg_s(1, ia, l)*my_CG_dxyz(1, iso1, iso2, iso) + &
                                                        gVXGg_s(2, ia, l)*my_CG_dxyz(2, iso1, iso2, iso) + &
                                                        gVXGg_s(3, ia, l)*my_CG_dxyz(3, iso1, iso2, iso))* &
                                                       harmonics%slm(ia, iso)

                              END DO ! ia

                              !test reduce expansion local densities

                           END DO ! icg
                        END DO ! iso
                        !test reduce expansion local densities
                     END IF ! lmax_expansion

                     !  Write in the global matrix
                     DO ic = nsoset(lmin(iset2) - 1) + 1, nsoset(lmax(iset2))
                        iso1 = nsoset(lmin(iset1) - 1) + 1
                        iso2 = ngau2 + ic
                        CALL daxpy(size1, 1.0_dp, matso_h(iso1, ic), 1, &
                                   int_hh(ispin)%r_coef(nngau1 + 1, iso2), 1)
                        CALL daxpy(size1, 1.0_dp, matso_s(iso1, ic), 1, &
                                   int_ss(ispin)%r_coef(nngau1 + 1, iso2), 1)
                     END DO

                  END DO ! ipfg2
               END DO ! ipfg1
               m2 = m2 + maxso
            END DO ! iset2
            m1 = m1 + maxso
         END DO ! iset1
      END DO ! ispin

      DEALLOCATE (g1, g2, gg, dgg, matso_h, matso_s, gVXCg_h, gVXCg_s, gVXGg_h, gVXGg_s)
      DEALLOCATE (cg_list, cg_n_list, dcg_list, dcg_n_list)

      CALL timestop(handle)

   END SUBROUTINE gaVxcgb_GC

! **************************************************************************************************
!> \brief Integrates 0.5 * grad_ga .dot. (V_tau * grad_gb) on the atomic grid for meta-GGA
!> \param vtau_h the hard tau potential
!> \param vtau_s the soft tau potential
!> \param int_hh hard one-center matrix contribution
!> \param int_ss soft one-center matrix contribution
!> \param tau_cache precomputed compact one-center gradient basis
!> \param nspins number of spin channels
!> \note This is a rewrite to correct meta-GGA GAPW bug. This is more brute force than the original
!>       but makes sure that no corner is cut in terms of accuracy (A. Bussy)
! **************************************************************************************************
   SUBROUTINE dgaVtaudgb(vtau_h, vtau_s, int_hh, int_ss, tau_cache, nspins)

      REAL(dp), DIMENSION(:, :, :), POINTER              :: vtau_h, vtau_s
      TYPE(rho_atom_coeff), DIMENSION(:), POINTER        :: int_hh, int_ss
      TYPE(tau_basis_cache_type), INTENT(IN)             :: tau_cache
      INTEGER, INTENT(IN)                                :: nspins

      CHARACTER(len=*), PARAMETER                        :: routineN = 'dgaVtaudgb'

      INTEGER                                            :: dir, handle, ia, ibas, igrid, iold, ir, &
                                                            ispin, jbas, jold, max_old_basis, na, &
                                                            nbas, ngrid, nr
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: int_h, int_s, weighted_grad

      CALL timeset(routineN, handle)

      CPASSERT(ALLOCATED(tau_cache%grad))
      CPASSERT(ASSOCIATED(tau_cache%n2oindex))

      nr = tau_cache%nr
      na = tau_cache%na
      nbas = tau_cache%nsatbas
      ngrid = na*nr
      max_old_basis = MAXVAL(tau_cache%n2oindex)
      ALLOCATE (int_h(nbas, nbas), int_s(nbas, nbas), weighted_grad(ngrid, nbas))

      DO ispin = 1, nspins
         CPASSERT(SIZE(int_hh(ispin)%r_coef, 1) >= max_old_basis)
         CPASSERT(SIZE(int_hh(ispin)%r_coef, 2) >= max_old_basis)
         CPASSERT(SIZE(int_ss(ispin)%r_coef, 1) >= max_old_basis)
         CPASSERT(SIZE(int_ss(ispin)%r_coef, 2) >= max_old_basis)
         int_h = 0.0_dp
         int_s = 0.0_dp
         DO dir = 1, 3
!$OMP PARALLEL DO DEFAULT(NONE) COLLAPSE(3) SCHEDULE(STATIC) &
!$OMP SHARED(dir, ispin, na, nbas, nr, tau_cache, vtau_h, weighted_grad) &
!$OMP PRIVATE(ia, ibas, igrid, ir)
            DO ibas = 1, nbas
               DO ir = 1, nr
                  DO ia = 1, na
                     igrid = ia + (ir - 1)*na
                     weighted_grad(igrid, ibas) = vtau_h(ia, ir, ispin)* &
                                                  tau_cache%grad(igrid, ibas, dir)
                  END DO
               END DO
            END DO
!$OMP END PARALLEL DO
            CALL dgemm('T', 'N', nbas, nbas, ngrid, 0.5_dp, tau_cache%grad(:, :, dir), &
                       ngrid, weighted_grad, ngrid, 1.0_dp, int_h, nbas)

!$OMP PARALLEL DO DEFAULT(NONE) COLLAPSE(3) SCHEDULE(STATIC) &
!$OMP SHARED(dir, ispin, na, nbas, nr, tau_cache, vtau_s, weighted_grad) &
!$OMP PRIVATE(ia, ibas, igrid, ir)
            DO ibas = 1, nbas
               DO ir = 1, nr
                  DO ia = 1, na
                     igrid = ia + (ir - 1)*na
                     weighted_grad(igrid, ibas) = vtau_s(ia, ir, ispin)* &
                                                  tau_cache%grad(igrid, ibas, dir)
                  END DO
               END DO
            END DO
!$OMP END PARALLEL DO
            CALL dgemm('T', 'N', nbas, nbas, ngrid, 0.5_dp, tau_cache%grad(:, :, dir), &
                       ngrid, weighted_grad, ngrid, 1.0_dp, int_s, nbas)
         END DO

!$OMP PARALLEL DO DEFAULT(NONE) COLLAPSE(2) SCHEDULE(STATIC) &
!$OMP SHARED(int_h, int_hh, int_s, int_ss, ispin, nbas, tau_cache) &
!$OMP PRIVATE(ibas, iold, jbas, jold)
         DO jbas = 1, nbas
            DO ibas = 1, nbas
               jold = tau_cache%n2oindex(jbas)
               iold = tau_cache%n2oindex(ibas)
               int_hh(ispin)%r_coef(iold, jold) = int_hh(ispin)%r_coef(iold, jold) + &
                                                  int_h(ibas, jbas)
               int_ss(ispin)%r_coef(iold, jold) = int_ss(ispin)%r_coef(iold, jold) + &
                                                  int_s(ibas, jbas)
            END DO
         END DO
!$OMP END PARALLEL DO
      END DO

      DEALLOCATE (int_h, int_s, weighted_grad)

      CALL timestop(handle)

   END SUBROUTINE dgaVtaudgb

END MODULE qs_vxc_atom_utils
