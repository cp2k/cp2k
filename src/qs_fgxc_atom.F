!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief routines that build the integrals of the Fxc kernel calculated
!>      for the atomic density in the basis set of spherical primitives
! **************************************************************************************************
MODULE qs_fgxc_atom
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cp_control_types,                ONLY: dft_control_type
   USE input_constants,                 ONLY: xc_none
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_fxc_atom,                     ONLY: fxc_atom_calc
   USE qs_grid_atom,                    ONLY: grid_atom_type
   USE qs_harmonics_atom,               ONLY: harmonics_atom_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              has_nlcc,&
                                              qs_kind_type
   USE qs_rho_atom_types,               ONLY: get_rho_atom,&
                                              rho_atom_coeff,&
                                              rho_atom_type
   USE qs_vxc_atom_utils,               ONLY: &
        calc_rho_angular, calc_rho_nlcc, calc_tau_atom, create_tau_basis_cache, dgaVtaudgb, &
        gaVxcgb_GC, gaVxcgb_noGC, release_tau_basis_cache, tau_basis_cache_type
   USE util,                            ONLY: get_limit
   USE xc_atom,                         ONLY: fill_rho_set,&
                                              vxc_of_r_new,&
                                              xc_2nd_deriv_of_r,&
                                              xc_rho_set_atom_update
   USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                              xc_dset_create,&
                                              xc_dset_release,&
                                              xc_dset_zero_all
   USE xc_derivatives,                  ONLY: xc_functionals_get_needs
   USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
   USE xc_rho_set_types,                ONLY: xc_rho_set_create,&
                                              xc_rho_set_release,&
                                              xc_rho_set_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_fgxc_atom'

   PUBLIC :: fgxc_atom_calc, &
             fgxc_atom_diff

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param rho0_atom_set ...
!> \param rho1_atom_set ...
!> \param rho2_atom_set ...
!> \param kind_set ...
!> \param xc_section ...
!> \param is_triplet ...
!> \param accuracy ...
! **************************************************************************************************
   SUBROUTINE fgxc_atom_calc(qs_env, rho0_atom_set, rho1_atom_set, rho2_atom_set, &
                             kind_set, xc_section, is_triplet, accuracy)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(rho_atom_type), DIMENSION(:), POINTER         :: rho0_atom_set, rho1_atom_set, &
                                                            rho2_atom_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: kind_set
      TYPE(section_vals_type), OPTIONAL, POINTER         :: xc_section
      LOGICAL, INTENT(IN)                                :: is_triplet
      INTEGER, INTENT(IN)                                :: accuracy

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'fgxc_atom_calc'
      REAL(KIND=dp), PARAMETER                           :: epsrho = 5.e-4_dp

      INTEGER                                            :: bo(2), handle, iat, iatom, ikind, ir, &
                                                            istep, mspins, myfun, na, natom, nf, &
                                                            nr, ns, nspins, nstep, num_pe
      INTEGER, DIMENSION(2, 3)                           :: bounds
      INTEGER, DIMENSION(:), POINTER                     :: atom_list
      LOGICAL                                            :: accint, donlcc, gradient_f, lsd, nlcc, &
                                                            paw_atom, tau_f
      REAL(dp)                                           :: agr, alpha, beta, density_cut, exc_h, &
                                                            exc_s, gradient_cut, oeps1, oeps2, &
                                                            tau_cut
      REAL(dp), DIMENSION(1, 1, 1)                       :: tau_d
      REAL(dp), DIMENSION(1, 1, 1, 1)                    :: rho_d
      REAL(dp), DIMENSION(:, :), POINTER                 :: rho_nlcc, weight_h, weight_s
      REAL(dp), DIMENSION(:, :, :), POINTER              :: rho0_h, rho0_s, rho1_h, rho1_s, rho_h, &
                                                            rho_s, tau0_h, tau0_s, tau1_h, tau1_s, &
                                                            tau_h, tau_s, vtau_h, vtau_s, vxc_h, &
                                                            vxc_s
      REAL(dp), DIMENSION(:, :, :, :), POINTER           :: drho0_h, drho0_s, drho1_h, drho1_s, &
                                                            drho_h, drho_s, vxg_h, vxg_s
      REAL(KIND=dp), DIMENSION(-4:4)                     :: ak, bl
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      TYPE(gto_basis_set_type), POINTER                  :: basis_1c
      TYPE(harmonics_atom_type), POINTER                 :: harmonics
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(rho_atom_coeff), DIMENSION(:), POINTER        :: dr_h, dr_s, fint_hh, fint_ss, int_hh, &
                                                            int_ss, r_h, r_s
      TYPE(rho_atom_coeff), DIMENSION(:, :), POINTER     :: r_h_d, r_s_d
      TYPE(rho_atom_type), POINTER                       :: rho0_atom, rho1_atom, rho2_atom
      TYPE(section_vals_type), POINTER                   :: xc_fun_section
      TYPE(tau_basis_cache_type)                         :: tau_basis_cache
      TYPE(xc_derivative_set_type)                       :: deriv_set
      TYPE(xc_rho_cflags_type)                           :: needs
      TYPE(xc_rho_set_type)                              :: rho_set_h, rho_set_s

      CALL timeset(routineN, handle)

      NULLIFY (vtau_h, vtau_s)

      ak = 0.0_dp
      bl = 0.0_dp
      SELECT CASE (accuracy)
      CASE (:4)
         nstep = 2
         ak(-2:2) = [1.0_dp, -8.0_dp, 0.0_dp, 8.0_dp, -1.0_dp]/12.0_dp
         bl(-2:2) = [-1.0_dp, 16.0_dp, -30.0_dp, 16.0_dp, -1.0_dp]/12.0_dp
      CASE (5:7)
         nstep = 3
         ak(-3:3) = [-1.0_dp, 9.0_dp, -45.0_dp, 0.0_dp, 45.0_dp, -9.0_dp, 1.0_dp]/60.0_dp
         bl(-3:3) = [2.0_dp, -27.0_dp, 270.0_dp, -490.0_dp, 270.0_dp, -27.0_dp, 2.0_dp]/180.0_dp
      CASE (8:)
         nstep = 4
         ak(-4:4) = [1.0_dp, -32.0_dp/3.0_dp, 56.0_dp, -224.0_dp, 0.0_dp, &
                     224.0_dp, -56.0_dp, 32.0_dp/3.0_dp, -1.0_dp]/280.0_dp
         bl(-4:4) = [-1.0_dp, 128.0_dp/9.0_dp, -112.0_dp, 896.0_dp, -14350.0_dp/9.0_dp, &
                     896.0_dp, -112.0_dp, 128.0_dp/9.0_dp, -1.0_dp]/560.0_dp
      END SELECT
      oeps1 = 1.0_dp/epsrho
      oeps2 = 1.0_dp/(epsrho**2)

      CALL get_qs_env(qs_env=qs_env, &
                      dft_control=dft_control, &
                      para_env=para_env, &
                      atomic_kind_set=atomic_kind_set)

      xc_fun_section => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
      CALL section_vals_val_get(xc_fun_section, "_SECTION_PARAMETERS_", i_val=myfun)

      accint = dft_control%qs_control%gapw_control%accurate_xcint

      IF (myfun == xc_none) THEN
         ! no action needed?
      ELSE
         CALL section_vals_val_get(xc_section, "DENSITY_CUTOFF", r_val=density_cut)
         CALL section_vals_val_get(xc_section, "GRADIENT_CUTOFF", r_val=gradient_cut)
         CALL section_vals_val_get(xc_section, "TAU_CUTOFF", r_val=tau_cut)

         nlcc = has_nlcc(kind_set)
         lsd = dft_control%lsd
         nspins = dft_control%nspins
         mspins = nspins
         IF (is_triplet) THEN
            CPASSERT(nspins == 1)
            lsd = .TRUE.
            mspins = 2
         END IF
         needs = xc_functionals_get_needs(xc_fun_section, lsd=lsd, calc_potential=.TRUE.)
         gradient_f = (needs%drho .OR. needs%drho_spin)
         tau_f = (needs%tau .OR. needs%tau_spin)

         ! Here starts the loop over all the atoms
         DO ikind = 1, SIZE(atomic_kind_set)
            CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=atom_list, natom=natom)
            CALL get_qs_kind(kind_set(ikind), paw_atom=paw_atom, &
                             harmonics=harmonics, grid_atom=grid_atom)
            CALL get_qs_kind(kind_set(ikind), basis_set=basis_1c, basis_type="GAPW_1C")

            IF (.NOT. paw_atom) CYCLE

            nr = grid_atom%nr
            na = grid_atom%ng_sphere

            ! set integration weights
            IF (accint) THEN
               weight_h => grid_atom%weight
               alpha = dft_control%qs_control%gapw_control%aw(ikind)
               IF (ASSOCIATED(grid_atom%gapw_weight_s)) THEN
                  IF (grid_atom%gapw_weight_alpha /= alpha) DEALLOCATE (grid_atom%gapw_weight_s)
               END IF
               IF (.NOT. ASSOCIATED(grid_atom%gapw_weight_s)) THEN
                  ALLOCATE (grid_atom%gapw_weight_s(na, nr))
                  DO ir = 1, nr
                     agr = 1.0_dp - EXP(-alpha*grid_atom%rad2(ir))
                     grid_atom%gapw_weight_s(:, ir) = grid_atom%weight(:, ir)*agr
                  END DO
                  grid_atom%gapw_weight_alpha = alpha
               END IF
               weight_s => grid_atom%gapw_weight_s
            ELSE
               weight_h => grid_atom%weight
               weight_s => grid_atom%weight
            END IF

            ! Prepare the structures needed to calculate and store the xc derivatives

            ! Array dimension: here only one dimensional arrays are used,
            ! i.e. only the first column of deriv_data is read.
            ! The other to dimensions  are set to size equal 1
            bounds(1:2, 1:3) = 1
            bounds(2, 1) = na
            bounds(2, 2) = nr

            ! create a place where to put the derivatives
            CALL xc_dset_create(deriv_set, local_bounds=bounds)
            ! create the place where to store the argument for the functionals
            CALL xc_rho_set_create(rho_set_h, bounds, rho_cutoff=density_cut, &
                                   drho_cutoff=gradient_cut, tau_cutoff=tau_cut)
            CALL xc_rho_set_create(rho_set_s, bounds, rho_cutoff=density_cut, &
                                   drho_cutoff=gradient_cut, tau_cutoff=tau_cut)

            ! allocate the required 3d arrays where to store rho and drho
            CALL xc_rho_set_atom_update(rho_set_h, needs, mspins, bounds)
            CALL xc_rho_set_atom_update(rho_set_s, needs, mspins, bounds)

            ALLOCATE (rho_h(na, nr, mspins), rho_s(na, nr, mspins), &
                      rho0_h(na, nr, nspins), rho0_s(na, nr, nspins), &
                      rho1_h(na, nr, nspins), rho1_s(na, nr, nspins))
            ALLOCATE (vxc_h(na, nr, mspins), vxc_s(na, nr, mspins))
            IF (gradient_f) THEN
               ALLOCATE (drho_h(4, na, nr, mspins), drho_s(4, na, nr, mspins), &
                         drho0_h(4, na, nr, nspins), drho0_s(4, na, nr, nspins), &
                         drho1_h(4, na, nr, nspins), drho1_s(4, na, nr, nspins))
               ALLOCATE (vxg_h(3, na, nr, mspins), vxg_s(3, na, nr, mspins))
            END IF
            IF (tau_f) THEN
               CALL create_tau_basis_cache(tau_basis_cache, grid_atom, basis_1c, harmonics)
               ALLOCATE (tau_h(na, nr, mspins), tau_s(na, nr, mspins), &
                         tau0_h(na, nr, nspins), tau0_s(na, nr, nspins), &
                         tau1_h(na, nr, nspins), tau1_s(na, nr, nspins))
               ALLOCATE (vtau_h(na, nr, mspins), vtau_s(na, nr, mspins))
            END IF
            !
            ! NLCC: prepare rho and drho of the core charge for this KIND
            donlcc = .FALSE.
            IF (nlcc) THEN
               NULLIFY (rho_nlcc)
               rho_nlcc => kind_set(ikind)%nlcc_pot
               IF (ASSOCIATED(rho_nlcc)) donlcc = .TRUE.
            END IF

            ! Distribute the atoms of this kind
            num_pe = para_env%num_pe
            bo = get_limit(natom, num_pe, para_env%mepos)

            DO iat = bo(1), bo(2)
               iatom = atom_list(iat)
               !
               NULLIFY (int_hh, int_ss)
               rho0_atom => rho0_atom_set(iatom)
               CALL get_rho_atom(rho_atom=rho0_atom, ga_Vlocal_gb_h=int_hh, ga_Vlocal_gb_s=int_ss)
               ALLOCATE (fint_ss(nspins), fint_hh(nspins))
               DO ns = 1, nspins
                  nf = SIZE(int_ss(ns)%r_coef, 1)
                  ALLOCATE (fint_ss(ns)%r_coef(nf, nf))
                  nf = SIZE(int_hh(ns)%r_coef, 1)
                  ALLOCATE (fint_hh(ns)%r_coef(nf, nf))
               END DO

               ! RHO0
               rho0_h = 0.0_dp
               rho0_s = 0.0_dp
               rho0_atom => rho0_atom_set(iatom)
               IF (gradient_f) THEN
                  NULLIFY (r_h, r_s, dr_h, dr_s, r_h_d, r_s_d)
                  CALL get_rho_atom(rho_atom=rho0_atom, rho_rad_h=r_h, rho_rad_s=r_s, drho_rad_h=dr_h, &
                                    drho_rad_s=dr_s, rho_rad_h_d=r_h_d, rho_rad_s_d=r_s_d)
                  drho0_h = 0.0_dp
                  drho0_s = 0.0_dp
               ELSE
                  NULLIFY (r_h, r_s)
                  CALL get_rho_atom(rho_atom=rho0_atom, rho_rad_h=r_h, rho_rad_s=r_s)
                  rho_d = 0.0_dp
               END IF
               DO ir = 1, nr
                  CALL calc_rho_angular(grid_atom, harmonics, nspins, gradient_f, &
                                        ir, r_h, r_s, rho0_h, rho0_s, dr_h, dr_s, &
                                        r_h_d, r_s_d, drho0_h, drho0_s)
                  IF (donlcc) THEN
                     CALL calc_rho_nlcc(grid_atom, nspins, gradient_f, &
                                        ir, rho_nlcc(:, 1), rho0_h, rho0_s, rho_nlcc(:, 2), drho0_h, drho0_s)
                  END IF
               END DO
               IF (tau_f) THEN
                  !compute tau on the grid all at once
                  CALL calc_tau_atom(tau0_h, tau0_s, rho0_atom, tau_basis_cache, nspins)
               ELSE
                  tau_d = 0.0_dp
               END IF
               ! RHO1
               rho1_h = 0.0_dp
               rho1_s = 0.0_dp
               rho1_atom => rho1_atom_set(iatom)
               IF (gradient_f) THEN
                  NULLIFY (r_h, r_s, dr_h, dr_s, r_h_d, r_s_d)
                  CALL get_rho_atom(rho_atom=rho1_atom, rho_rad_h=r_h, rho_rad_s=r_s, drho_rad_h=dr_h, &
                                    drho_rad_s=dr_s, rho_rad_h_d=r_h_d, rho_rad_s_d=r_s_d)
                  drho1_h = 0.0_dp
                  drho1_s = 0.0_dp
               ELSE
                  NULLIFY (r_h, r_s)
                  CALL get_rho_atom(rho_atom=rho1_atom, rho_rad_h=r_h, rho_rad_s=r_s)
               END IF
               DO ir = 1, nr
                  CALL calc_rho_angular(grid_atom, harmonics, nspins, gradient_f, &
                                        ir, r_h, r_s, rho1_h, rho1_s, dr_h, dr_s, &
                                        r_h_d, r_s_d, drho1_h, drho1_s)
               END DO
               IF (tau_f) THEN
                  !compute tau on the grid all at once
                  CALL calc_tau_atom(tau1_h, tau1_s, rho1_atom, tau_basis_cache, nspins)
               END IF
               ! RHO2
               rho2_atom => rho2_atom_set(iatom)

               DO istep = -nstep, nstep

                  beta = REAL(istep, KIND=dp)*epsrho

                  IF (is_triplet) THEN
                     rho_h(:, :, 1) = rho0_h(:, :, 1) + beta*rho1_h(:, :, 1)
                     rho_h(:, :, 2) = rho0_h(:, :, 1)
                     rho_h = 0.5_dp*rho_h
                     rho_s(:, :, 1) = rho0_s(:, :, 1) + beta*rho1_s(:, :, 1)
                     rho_s(:, :, 2) = rho0_s(:, :, 1)
                     rho_s = 0.5_dp*rho_s
                     IF (gradient_f) THEN
                        drho_h(:, :, :, 1) = drho0_h(:, :, :, 1) + beta*drho1_h(:, :, :, 1)
                        drho_h(:, :, :, 2) = drho0_h(:, :, :, 1)
                        drho_h = 0.5_dp*drho_h
                        drho_s(:, :, :, 1) = drho0_s(:, :, :, 1) + beta*drho1_s(:, :, :, 1)
                        drho_s(:, :, :, 2) = drho0_s(:, :, :, 1)
                        drho_s = 0.5_dp*drho_s
                     END IF
                     IF (tau_f) THEN
                        tau_h(:, :, 1) = tau0_h(:, :, 1) + beta*tau1_h(:, :, 1)
                        tau_h(:, :, 2) = tau0_h(:, :, 1)
                        tau_h = 0.5_dp*tau0_h
                        tau_s(:, :, 1) = tau0_s(:, :, 1) + beta*tau1_s(:, :, 1)
                        tau_s(:, :, 2) = tau0_s(:, :, 1)
                        tau_s = 0.5_dp*tau0_s
                     END IF
                  ELSE
                     rho_h = rho0_h + beta*rho1_h
                     rho_s = rho0_s + beta*rho1_s
                     IF (gradient_f) THEN
                        drho_h = drho0_h + beta*drho1_h
                        drho_s = drho0_s + beta*drho1_s
                     END IF
                     IF (tau_f) THEN
                        tau_h = tau0_h + beta*tau1_h
                        tau_s = tau0_s + beta*tau1_s
                     END IF
                  END IF
                  !
                  IF (gradient_f) THEN
                     drho_h(4, :, :, :) = NORM2(drho_h(1:3, :, :, :), 1)
                     drho_s(4, :, :, :) = NORM2(drho_s(1:3, :, :, :), 1)
                  END IF

                  DO ir = 1, nr
                     IF (tau_f) THEN
                        CALL fill_rho_set(rho_set_h, lsd, mspins, needs, rho_h, drho_h, tau_h, na, ir)
                        CALL fill_rho_set(rho_set_s, lsd, mspins, needs, rho_s, drho_s, tau_s, na, ir)
                     ELSE IF (gradient_f) THEN
                        CALL fill_rho_set(rho_set_h, lsd, mspins, needs, rho_h, drho_h, tau_d, na, ir)
                        CALL fill_rho_set(rho_set_s, lsd, mspins, needs, rho_s, drho_s, tau_d, na, ir)
                     ELSE
                        CALL fill_rho_set(rho_set_h, lsd, mspins, needs, rho_h, rho_d, tau_d, na, ir)
                        CALL fill_rho_set(rho_set_s, lsd, mspins, needs, rho_s, rho_d, tau_d, na, ir)
                     END IF
                  END DO

                  ! hard atom density !
                  CALL xc_dset_zero_all(deriv_set)
                  CALL vxc_of_r_new(xc_fun_section, rho_set_h, deriv_set, 1, needs, weight_h, &
                                    lsd, na, nr, exc_h, vxc_h, vxg_h, vtau_h)
                  IF (is_triplet) THEN
                     vxc_h(:, :, 1) = vxc_h(:, :, 1) - vxc_h(:, :, 2)
                     IF (gradient_f) THEN
                        vxg_h(:, :, :, 1) = vxg_h(:, :, :, 1) - vxg_h(:, :, :, 2)
                     END IF
                     IF (tau_f) THEN
                        vtau_h(:, :, 1) = vtau_h(:, :, 1) - vtau_h(:, :, 2)
                     END IF
                  END IF
                  ! soft atom density !
                  CALL xc_dset_zero_all(deriv_set)
                  CALL vxc_of_r_new(xc_fun_section, rho_set_s, deriv_set, 1, needs, weight_s, &
                                    lsd, na, nr, exc_s, vxc_s, vxg_s, vtau_s)
                  IF (is_triplet) THEN
                     vxc_s(:, :, 1) = vxc_s(:, :, 1) - vxc_s(:, :, 2)
                     IF (gradient_f) THEN
                        vxg_s(:, :, :, 1) = vxg_s(:, :, :, 1) - vxg_s(:, :, :, 2)
                     END IF
                     IF (tau_f) THEN
                        vtau_s(:, :, 1) = vtau_s(:, :, 1) - vtau_s(:, :, 2)
                     END IF
                  END IF
                  ! potentials
                  DO ns = 1, nspins
                     fint_hh(ns)%r_coef(:, :) = 0.0_dp
                     fint_ss(ns)%r_coef(:, :) = 0.0_dp
                  END DO
                  IF (gradient_f) THEN
                     CALL gaVxcgb_GC(vxc_h, vxc_s, vxg_h, vxg_s, fint_hh, fint_ss, &
                                     grid_atom, basis_1c, harmonics, nspins)
                  ELSE
                     CALL gaVxcgb_noGC(vxc_h, vxc_s, fint_hh, fint_ss, &
                                       grid_atom, basis_1c, harmonics, nspins)
                  END IF
                  IF (tau_f) THEN
                     CALL dgaVtaudgb(vtau_h, vtau_s, fint_hh, fint_ss, &
                                     tau_basis_cache, nspins)
                  END IF
                  ! first derivative fxc
                  NULLIFY (int_hh, int_ss)
                  CALL get_rho_atom(rho_atom=rho1_atom, ga_Vlocal_gb_h=int_hh, ga_Vlocal_gb_s=int_ss)
                  DO ns = 1, nspins
                     int_ss(ns)%r_coef(:, :) = int_ss(ns)%r_coef(:, :) + oeps1*ak(istep)*fint_ss(ns)%r_coef(:, :)
                     int_hh(ns)%r_coef(:, :) = int_hh(ns)%r_coef(:, :) + oeps1*ak(istep)*fint_hh(ns)%r_coef(:, :)
                  END DO
                  ! second derivative gxc
                  NULLIFY (int_hh, int_ss)
                  CALL get_rho_atom(rho_atom=rho2_atom, ga_Vlocal_gb_h=int_hh, ga_Vlocal_gb_s=int_ss)
                  DO ns = 1, nspins
                     int_ss(ns)%r_coef(:, :) = int_ss(ns)%r_coef(:, :) + oeps2*bl(istep)*fint_ss(ns)%r_coef(:, :)
                     int_hh(ns)%r_coef(:, :) = int_hh(ns)%r_coef(:, :) + oeps2*bl(istep)*fint_hh(ns)%r_coef(:, :)
                  END DO
               END DO
               !
               DO ns = 1, nspins
                  DEALLOCATE (fint_ss(ns)%r_coef)
                  DEALLOCATE (fint_hh(ns)%r_coef)
               END DO
               DEALLOCATE (fint_ss, fint_hh)

            END DO ! iat

            ! Release the xc structure used to store the xc derivatives
            CALL xc_dset_release(deriv_set)
            CALL xc_rho_set_release(rho_set_h)
            CALL xc_rho_set_release(rho_set_s)

            DEALLOCATE (rho_h, rho_s, rho0_h, rho0_s, rho1_h, rho1_s)
            DEALLOCATE (vxc_h, vxc_s)
            IF (gradient_f) THEN
               DEALLOCATE (drho_h, drho_s, drho0_h, drho0_s, drho1_h, drho1_s)
               DEALLOCATE (vxg_h, vxg_s)
            END IF
            IF (tau_f) THEN
               DEALLOCATE (tau_h, tau_s, tau0_h, tau0_s, tau1_h, tau1_s)
               DEALLOCATE (vtau_h, vtau_s)
               CALL release_tau_basis_cache(tau_basis_cache)
            END IF
         END DO ! ikind

      END IF !xc_none

      CALL timestop(handle)

   END SUBROUTINE fgxc_atom_calc

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param rho0_atom_set ...
!> \param rho1_atom_set ...
!> \param rho2_atom_set ...
!> \param kind_set ...
!> \param xc_section ...
!> \param is_triplet ...
!> \param accuracy ...
!> \param epsrho ...
! **************************************************************************************************
   SUBROUTINE fgxc_atom_diff(qs_env, rho0_atom_set, rho1_atom_set, rho2_atom_set, &
                             kind_set, xc_section, is_triplet, accuracy, epsrho)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(rho_atom_type), DIMENSION(:), POINTER         :: rho0_atom_set, rho1_atom_set, &
                                                            rho2_atom_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: kind_set
      TYPE(section_vals_type), OPTIONAL, POINTER         :: xc_section
      LOGICAL, INTENT(IN)                                :: is_triplet
      INTEGER, INTENT(IN)                                :: accuracy
      REAL(KIND=dp), INTENT(IN)                          :: epsrho

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'fgxc_atom_diff'

      INTEGER                                            :: bo(2), handle, iat, iatom, ikind, ir, &
                                                            istep, mspins, myfun, na, natom, nf, &
                                                            nr, ns, nspins, nstep, num_pe
      INTEGER, DIMENSION(2, 3)                           :: bounds
      INTEGER, DIMENSION(:), POINTER                     :: atom_list
      LOGICAL                                            :: accint, donlcc, gradient_f, lsd, nlcc, &
                                                            paw_atom, tau_f
      REAL(dp)                                           :: agr, alpha, beta, density_cut, &
                                                            gradient_cut, oeps1, tau_cut
      REAL(dp), CONTIGUOUS, DIMENSION(:, :, :), POINTER  :: vtau_h, vtau_s, vxc_h, vxc_s
      REAL(dp), DIMENSION(1, 1, 1)                       :: tau_d
      REAL(dp), DIMENSION(1, 1, 1, 1)                    :: rho_d
      REAL(dp), DIMENSION(:, :), POINTER                 :: rho_nlcc, weight_h, weight_s
      REAL(dp), DIMENSION(:, :, :), POINTER              :: rho0_h, rho0_s, rho1_h, rho1_s, rho_h, &
                                                            rho_s, tau0_h, tau0_s, tau1_h, tau1_s, &
                                                            tau_h, tau_s
      REAL(dp), DIMENSION(:, :, :, :), POINTER           :: drho0_h, drho0_s, drho1_h, drho1_s, &
                                                            drho_h, drho_s, vxg_h, vxg_s
      REAL(KIND=dp), DIMENSION(-4:4)                     :: ak
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      TYPE(gto_basis_set_type), POINTER                  :: basis_1c
      TYPE(harmonics_atom_type), POINTER                 :: harmonics
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(rho_atom_coeff), DIMENSION(:), POINTER        :: dr_h, dr_s, fint_hh, fint_ss, int_hh, &
                                                            int_ss, r_h, r_s
      TYPE(rho_atom_coeff), DIMENSION(:, :), POINTER     :: r_h_d, r_s_d
      TYPE(rho_atom_type), POINTER                       :: rho0_atom, rho1_atom, rho2_atom
      TYPE(section_vals_type), POINTER                   :: xc_fun_section
      TYPE(tau_basis_cache_type)                         :: tau_basis_cache
      TYPE(xc_derivative_set_type)                       :: deriv_set
      TYPE(xc_rho_cflags_type)                           :: needs
      TYPE(xc_rho_set_type)                              :: rho1_set_h, rho1_set_s, rho_set_h, &
                                                            rho_set_s

      CALL timeset(routineN, handle)

      NULLIFY (vtau_h, vtau_s)

      ak = 0.0_dp
      SELECT CASE (accuracy)
      CASE (:4)
         nstep = 2
         ak(-2:2) = [1.0_dp, -8.0_dp, 0.0_dp, 8.0_dp, -1.0_dp]/12.0_dp
      CASE (5:7)
         nstep = 3
         ak(-3:3) = [-1.0_dp, 9.0_dp, -45.0_dp, 0.0_dp, 45.0_dp, -9.0_dp, 1.0_dp]/60.0_dp
      CASE (8:)
         nstep = 4
         ak(-4:4) = [1.0_dp, -32.0_dp/3.0_dp, 56.0_dp, -224.0_dp, 0.0_dp, &
                     224.0_dp, -56.0_dp, 32.0_dp/3.0_dp, -1.0_dp]/280.0_dp
      END SELECT
      oeps1 = 1.0_dp/epsrho

      CALL get_qs_env(qs_env=qs_env, &
                      dft_control=dft_control, &
                      para_env=para_env, &
                      atomic_kind_set=atomic_kind_set)

      xc_fun_section => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
      CALL section_vals_val_get(xc_fun_section, "_SECTION_PARAMETERS_", i_val=myfun)

      accint = dft_control%qs_control%gapw_control%accurate_xcint

      IF (myfun == xc_none) THEN
         ! no action needed?
      ELSE
         ! calculate fxc
         CALL fxc_atom_calc(qs_env, rho0_atom_set, rho1_atom_set, xc_section, para_env, &
                            do_triplet=is_triplet, kind_set_external=kind_set)

         CALL section_vals_val_get(xc_section, "DENSITY_CUTOFF", r_val=density_cut)
         CALL section_vals_val_get(xc_section, "GRADIENT_CUTOFF", r_val=gradient_cut)
         CALL section_vals_val_get(xc_section, "TAU_CUTOFF", r_val=tau_cut)

         nlcc = has_nlcc(kind_set)
         lsd = dft_control%lsd
         nspins = dft_control%nspins
         mspins = nspins
         IF (is_triplet) THEN
            CPASSERT(nspins == 1)
            lsd = .TRUE.
            mspins = 2
         END IF
         needs = xc_functionals_get_needs(xc_fun_section, lsd=lsd, calc_potential=.TRUE.)
         gradient_f = (needs%drho .OR. needs%drho_spin)
         tau_f = (needs%tau .OR. needs%tau_spin)

         ! Here starts the loop over all the atoms
         DO ikind = 1, SIZE(atomic_kind_set)
            CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=atom_list, natom=natom)
            CALL get_qs_kind(kind_set(ikind), paw_atom=paw_atom, &
                             harmonics=harmonics, grid_atom=grid_atom)
            CALL get_qs_kind(kind_set(ikind), basis_set=basis_1c, basis_type="GAPW_1C")

            IF (.NOT. paw_atom) CYCLE

            nr = grid_atom%nr
            na = grid_atom%ng_sphere

            ! set integration weights
            IF (accint) THEN
               weight_h => grid_atom%weight
               alpha = dft_control%qs_control%gapw_control%aw(ikind)
               IF (ASSOCIATED(grid_atom%gapw_weight_s)) THEN
                  IF (grid_atom%gapw_weight_alpha /= alpha) DEALLOCATE (grid_atom%gapw_weight_s)
               END IF
               IF (.NOT. ASSOCIATED(grid_atom%gapw_weight_s)) THEN
                  ALLOCATE (grid_atom%gapw_weight_s(na, nr))
                  DO ir = 1, nr
                     agr = 1.0_dp - EXP(-alpha*grid_atom%rad2(ir))
                     grid_atom%gapw_weight_s(:, ir) = grid_atom%weight(:, ir)*agr
                  END DO
                  grid_atom%gapw_weight_alpha = alpha
               END IF
               weight_s => grid_atom%gapw_weight_s
            ELSE
               weight_h => grid_atom%weight
               weight_s => grid_atom%weight
            END IF

            ! Prepare the structures needed to calculate and store the xc derivatives

            ! Array dimension: here anly one dimensional arrays are used,
            ! i.e. only the first column of deriv_data is read.
            ! The other two dimensions are set to size equal 1
            bounds(1:2, 1:3) = 1
            bounds(2, 1) = na
            bounds(2, 2) = nr

            ! create a place where to put the derivatives
            CALL xc_dset_create(deriv_set, local_bounds=bounds)
            ! create the place where to store the argument for the functionals
            CALL xc_rho_set_create(rho_set_h, bounds, rho_cutoff=density_cut, &
                                   drho_cutoff=gradient_cut, tau_cutoff=tau_cut)
            CALL xc_rho_set_create(rho_set_s, bounds, rho_cutoff=density_cut, &
                                   drho_cutoff=gradient_cut, tau_cutoff=tau_cut)
            CALL xc_rho_set_create(rho1_set_h, bounds, rho_cutoff=density_cut, &
                                   drho_cutoff=gradient_cut, tau_cutoff=tau_cut)
            CALL xc_rho_set_create(rho1_set_s, bounds, rho_cutoff=density_cut, &
                                   drho_cutoff=gradient_cut, tau_cutoff=tau_cut)

            ! allocate the required 3d arrays where to store rho and drho
            CALL xc_rho_set_atom_update(rho_set_h, needs, mspins, bounds)
            CALL xc_rho_set_atom_update(rho_set_s, needs, mspins, bounds)
            CALL xc_rho_set_atom_update(rho1_set_h, needs, mspins, bounds)
            CALL xc_rho_set_atom_update(rho1_set_s, needs, mspins, bounds)

            ALLOCATE (rho_h(na, nr, nspins), rho_s(na, nr, nspins), &
                      rho0_h(na, nr, nspins), rho0_s(na, nr, nspins), &
                      rho1_h(na, nr, nspins), rho1_s(na, nr, nspins))
            ALLOCATE (vxc_h(na, nr, nspins), vxc_s(na, nr, nspins))
            IF (gradient_f) THEN
               ALLOCATE (drho_h(4, na, nr, nspins), drho_s(4, na, nr, nspins), &
                         drho0_h(4, na, nr, nspins), drho0_s(4, na, nr, nspins), &
                         drho1_h(4, na, nr, nspins), drho1_s(4, na, nr, nspins))
               ALLOCATE (vxg_h(3, na, nr, nspins), vxg_s(3, na, nr, nspins))
            END IF
            IF (tau_f) THEN
               CALL create_tau_basis_cache(tau_basis_cache, grid_atom, basis_1c, harmonics)
               ALLOCATE (tau_h(na, nr, nspins), tau_s(na, nr, nspins), &
                         tau0_h(na, nr, nspins), tau0_s(na, nr, nspins), &
                         tau1_h(na, nr, nspins), tau1_s(na, nr, nspins))
               ALLOCATE (vtau_h(na, nr, nspins), vtau_s(na, nr, nspins))
            END IF
            !
            ! NLCC: prepare rho and drho of the core charge for this KIND
            donlcc = .FALSE.
            IF (nlcc) THEN
               NULLIFY (rho_nlcc)
               rho_nlcc => kind_set(ikind)%nlcc_pot
               IF (ASSOCIATED(rho_nlcc)) donlcc = .TRUE.
            END IF

            ! Distribute the atoms of this kind
            num_pe = para_env%num_pe
            bo = get_limit(natom, num_pe, para_env%mepos)

            DO iat = bo(1), bo(2)
               iatom = atom_list(iat)
               !
               NULLIFY (int_hh, int_ss)
               rho0_atom => rho0_atom_set(iatom)
               CALL get_rho_atom(rho_atom=rho0_atom, ga_Vlocal_gb_h=int_hh, ga_Vlocal_gb_s=int_ss)
               ALLOCATE (fint_ss(nspins), fint_hh(nspins))
               DO ns = 1, nspins
                  nf = SIZE(int_ss(ns)%r_coef, 1)
                  ALLOCATE (fint_ss(ns)%r_coef(nf, nf))
                  nf = SIZE(int_hh(ns)%r_coef, 1)
                  ALLOCATE (fint_hh(ns)%r_coef(nf, nf))
               END DO

               ! RHO0
               rho0_h = 0.0_dp
               rho0_s = 0.0_dp
               rho0_atom => rho0_atom_set(iatom)
               IF (gradient_f) THEN
                  NULLIFY (r_h, r_s, dr_h, dr_s, r_h_d, r_s_d)
                  CALL get_rho_atom(rho_atom=rho0_atom, rho_rad_h=r_h, rho_rad_s=r_s, drho_rad_h=dr_h, &
                                    drho_rad_s=dr_s, rho_rad_h_d=r_h_d, rho_rad_s_d=r_s_d)
                  drho0_h = 0.0_dp
                  drho0_s = 0.0_dp
               ELSE
                  NULLIFY (r_h, r_s)
                  CALL get_rho_atom(rho_atom=rho0_atom, rho_rad_h=r_h, rho_rad_s=r_s)
                  rho_d = 0.0_dp
               END IF
               DO ir = 1, nr
                  CALL calc_rho_angular(grid_atom, harmonics, nspins, gradient_f, &
                                        ir, r_h, r_s, rho0_h, rho0_s, dr_h, dr_s, &
                                        r_h_d, r_s_d, drho0_h, drho0_s)
                  IF (donlcc) THEN
                     CALL calc_rho_nlcc(grid_atom, nspins, gradient_f, &
                                        ir, rho_nlcc(:, 1), rho0_h, rho0_s, rho_nlcc(:, 2), drho0_h, drho0_s)
                  END IF
               END DO
               IF (tau_f) THEN
                  !compute tau on the grid all at once
                  CALL calc_tau_atom(tau0_h, tau0_s, rho0_atom, tau_basis_cache, nspins)
               ELSE
                  tau_d = 0.0_dp
               END IF
               ! RHO1
               rho1_h = 0.0_dp
               rho1_s = 0.0_dp
               rho1_atom => rho1_atom_set(iatom)
               IF (gradient_f) THEN
                  NULLIFY (r_h, r_s, dr_h, dr_s, r_h_d, r_s_d)
                  CALL get_rho_atom(rho_atom=rho1_atom, rho_rad_h=r_h, rho_rad_s=r_s, drho_rad_h=dr_h, &
                                    drho_rad_s=dr_s, rho_rad_h_d=r_h_d, rho_rad_s_d=r_s_d)
                  drho1_h = 0.0_dp
                  drho1_s = 0.0_dp
               ELSE
                  NULLIFY (r_h, r_s)
                  CALL get_rho_atom(rho_atom=rho1_atom, rho_rad_h=r_h, rho_rad_s=r_s)
               END IF
               DO ir = 1, nr
                  CALL calc_rho_angular(grid_atom, harmonics, nspins, gradient_f, &
                                        ir, r_h, r_s, rho1_h, rho1_s, dr_h, dr_s, &
                                        r_h_d, r_s_d, drho1_h, drho1_s)
               END DO
               IF (tau_f) THEN
                  !compute tau on the grid all at once
                  CALL calc_tau_atom(tau1_h, tau1_s, rho1_atom, tau_basis_cache, nspins)
               END IF

               DO ir = 1, nr
                  IF (tau_f) THEN
                     CALL fill_rho_set(rho1_set_h, lsd, nspins, needs, rho1_h, drho1_h, tau1_h, na, ir)
                     CALL fill_rho_set(rho1_set_s, lsd, nspins, needs, rho1_s, drho1_s, tau1_s, na, ir)
                  ELSE IF (gradient_f) THEN
                     CALL fill_rho_set(rho1_set_h, lsd, nspins, needs, rho1_h, drho1_h, tau_d, na, ir)
                     CALL fill_rho_set(rho1_set_s, lsd, nspins, needs, rho1_s, drho1_s, tau_d, na, ir)
                  ELSE
                     CALL fill_rho_set(rho1_set_h, lsd, nspins, needs, rho1_h, rho_d, tau_d, na, ir)
                     CALL fill_rho_set(rho1_set_s, lsd, nspins, needs, rho1_s, rho_d, tau_d, na, ir)
                  END IF
               END DO

               ! RHO2
               rho2_atom => rho2_atom_set(iatom)

               DO istep = -nstep, nstep

                  beta = REAL(istep, KIND=dp)*epsrho

                  rho_h = rho0_h + beta*rho1_h
                  rho_s = rho0_s + beta*rho1_s
                  IF (gradient_f) THEN
                     drho_h = drho0_h + beta*drho1_h
                     drho_s = drho0_s + beta*drho1_s
                  END IF
                  IF (tau_f) THEN
                     tau_h = tau0_h + beta*tau1_h
                     tau_s = tau0_s + beta*tau1_s
                  END IF
                  !
                  IF (gradient_f) THEN
                     drho_h(4, :, :, :) = NORM2(drho_h(1:3, :, :, :), 1)
                     drho_s(4, :, :, :) = NORM2(drho_s(1:3, :, :, :), 1)
                  END IF

                  DO ir = 1, nr
                     IF (tau_f) THEN
                        CALL fill_rho_set(rho_set_h, lsd, nspins, needs, rho_h, drho_h, tau_h, na, ir)
                        CALL fill_rho_set(rho_set_s, lsd, nspins, needs, rho_s, drho_s, tau_s, na, ir)
                     ELSE IF (gradient_f) THEN
                        CALL fill_rho_set(rho_set_h, lsd, nspins, needs, rho_h, drho_h, tau_d, na, ir)
                        CALL fill_rho_set(rho_set_s, lsd, nspins, needs, rho_s, drho_s, tau_d, na, ir)
                     ELSE
                        CALL fill_rho_set(rho_set_h, lsd, nspins, needs, rho_h, rho_d, tau_d, na, ir)
                        CALL fill_rho_set(rho_set_s, lsd, nspins, needs, rho_s, rho_d, tau_d, na, ir)
                     END IF
                  END DO

                  ! hard atom density !
                  CALL xc_dset_zero_all(deriv_set)
                  CALL xc_2nd_deriv_of_r(xc_section=xc_section, &
                                         rho_set=rho_set_h, rho1_set=rho1_set_h, &
                                         deriv_set=deriv_set, &
                                         w=weight_h, vxc=vxc_h, vxg=vxg_h, vtau=vtau_h, &
                                         do_triplet=is_triplet)
                  ! soft atom density !
                  CALL xc_dset_zero_all(deriv_set)
                  CALL xc_2nd_deriv_of_r(xc_section=xc_section, &
                                         rho_set=rho_set_s, rho1_set=rho1_set_s, &
                                         deriv_set=deriv_set, &
                                         w=weight_s, vxc=vxc_s, vxg=vxg_s, vtau=vtau_s, &
                                         do_triplet=is_triplet)
                  ! potentials
                  DO ns = 1, nspins
                     fint_hh(ns)%r_coef(:, :) = 0.0_dp
                     fint_ss(ns)%r_coef(:, :) = 0.0_dp
                  END DO
                  IF (gradient_f) THEN
                     CALL gaVxcgb_GC(vxc_h, vxc_s, vxg_h, vxg_s, fint_hh, fint_ss, &
                                     grid_atom, basis_1c, harmonics, nspins)
                  ELSE
                     CALL gaVxcgb_noGC(vxc_h, vxc_s, fint_hh, fint_ss, &
                                       grid_atom, basis_1c, harmonics, nspins)
                  END IF
                  IF (tau_f) THEN
                     CALL dgaVtaudgb(vtau_h, vtau_s, fint_hh, fint_ss, &
                                     tau_basis_cache, nspins)
                  END IF
                  ! second derivative gxc
                  NULLIFY (int_hh, int_ss)
                  CALL get_rho_atom(rho_atom=rho2_atom, ga_Vlocal_gb_h=int_hh, ga_Vlocal_gb_s=int_ss)
                  DO ns = 1, nspins
                     int_ss(ns)%r_coef(:, :) = int_ss(ns)%r_coef(:, :) + oeps1*ak(istep)*fint_ss(ns)%r_coef(:, :)
                     int_hh(ns)%r_coef(:, :) = int_hh(ns)%r_coef(:, :) + oeps1*ak(istep)*fint_hh(ns)%r_coef(:, :)
                  END DO
               END DO
               !
               DO ns = 1, nspins
                  DEALLOCATE (fint_ss(ns)%r_coef)
                  DEALLOCATE (fint_hh(ns)%r_coef)
               END DO
               DEALLOCATE (fint_ss, fint_hh)

            END DO ! iat

            ! Release the xc structure used to store the xc derivatives
            CALL xc_dset_release(deriv_set)
            CALL xc_rho_set_release(rho_set_h)
            CALL xc_rho_set_release(rho_set_s)
            CALL xc_rho_set_release(rho1_set_h)
            CALL xc_rho_set_release(rho1_set_s)

            DEALLOCATE (rho_h, rho_s, rho0_h, rho0_s, rho1_h, rho1_s)
            DEALLOCATE (vxc_h, vxc_s)
            IF (gradient_f) THEN
               DEALLOCATE (drho_h, drho_s, drho0_h, drho0_s, drho1_h, drho1_s)
               DEALLOCATE (vxg_h, vxg_s)
            END IF
            IF (tau_f) THEN
               DEALLOCATE (tau_h, tau_s, tau0_h, tau0_s, tau1_h, tau1_s)
               DEALLOCATE (vtau_h, vtau_s)
               CALL release_tau_basis_cache(tau_basis_cache)
            END IF
         END DO ! ikind

      END IF !xc_none

      CALL timestop(handle)

   END SUBROUTINE fgxc_atom_diff

END MODULE qs_fgxc_atom
