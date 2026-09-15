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
MODULE qs_fxc_atom
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cp_control_types,                ONLY: dft_control_type
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
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
                                              xc_2nd_deriv_of_r,&
                                              xc_rho_set_atom_update
   USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                              xc_dset_create,&
                                              xc_dset_release
   USE xc_derivatives,                  ONLY: xc_functionals_get_needs
   USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
   USE xc_rho_set_types,                ONLY: xc_rho_set_create,&
                                              xc_rho_set_release,&
                                              xc_rho_set_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_fxc_atom'

   PUBLIC :: fxc_atom_calc

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param rho_atom_set ...
!> \param rho1_atom_set ...
!> \param xc_section ...
!> \param para_env_ext ...
!> \param do_scale ...
!> \param do_triplet ...
!> \param do_sf ...
!> \param kind_set_external ...
! **************************************************************************************************
   SUBROUTINE fxc_atom_calc(qs_env, rho_atom_set, rho1_atom_set, xc_section, para_env_ext, &
                            do_scale, do_triplet, do_sf, kind_set_external)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(rho_atom_type), DIMENSION(:), POINTER         :: rho_atom_set, rho1_atom_set
      TYPE(section_vals_type), POINTER                   :: xc_section
      TYPE(mp_para_env_type), INTENT(IN), OPTIONAL       :: para_env_ext
      LOGICAL, INTENT(IN), OPTIONAL                      :: do_scale, do_triplet, do_sf
      TYPE(qs_kind_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: kind_set_external

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'fxc_atom_calc'

      INTEGER                                            :: atom, handle, iatom, ikind, ir, na, &
                                                            natom, nr, nspins
      INTEGER, DIMENSION(2)                              :: local_loop_limit
      INTEGER, DIMENSION(2, 3)                           :: bounds
      INTEGER, DIMENSION(:), POINTER                     :: atom_list
      LOGICAL                                            :: accint, donlcc, gradient_functional, &
                                                            lsd, my_do_sf, nlcc, paw_atom, &
                                                            scale_rho, tau_f
      REAL(KIND=dp)                                      :: agr, alpha, density_cut, gradient_cut, &
                                                            rtot, tau_cut
      REAL(KIND=dp), CONTIGUOUS, DIMENSION(:, :, :), &
         POINTER                                         :: vtau_h, vtau_s, vxc_h, vxc_s
      REAL(KIND=dp), DIMENSION(1, 1, 1)                  :: rtau
      REAL(KIND=dp), DIMENSION(1, 1, 1, 1)               :: rrho
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rho_nlcc, weight_h, weight_s
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: rho1_h, rho1_s, rho_h, rho_s, tau1_h, &
                                                            tau1_s, tau_h, tau_s
      REAL(KIND=dp), DIMENSION(:, :, :, :), POINTER      :: drho1_h, drho1_s, drho_h, drho_s, vxg_h, &
                                                            vxg_s
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      TYPE(gto_basis_set_type), POINTER                  :: basis_1c
      TYPE(harmonics_atom_type), POINTER                 :: harmonics
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: my_kind_set, qs_kind_set
      TYPE(rho_atom_coeff), DIMENSION(:), POINTER        :: dr1_h, dr1_s, dr_h, dr_s, int_hh, &
                                                            int_ss, r1_h, r1_s, r_h, r_s
      TYPE(rho_atom_coeff), DIMENSION(:, :), POINTER     :: r1_h_d, r1_s_d, r_h_d, r_s_d
      TYPE(rho_atom_type), POINTER                       :: rho1_atom, rho_atom
      TYPE(section_vals_type), POINTER                   :: input, xc_fun_section
      TYPE(tau_basis_cache_type)                         :: tau_basis_cache
      TYPE(xc_derivative_set_type)                       :: deriv_set
      TYPE(xc_rho_cflags_type)                           :: needs
      TYPE(xc_rho_set_type)                              :: rho1_set_h, rho1_set_s, rho_set_h, &
                                                            rho_set_s

! -------------------------------------------------------------------------

      CALL timeset(routineN, handle)

      NULLIFY (qs_kind_set)
      NULLIFY (rho_h, rho_s, drho_h, drho_s, weight_h, weight_s)
      NULLIFY (rho1_h, rho1_s, drho1_h, drho1_s)
      NULLIFY (vxc_h, vxc_s, vxg_h, vxg_s)
      NULLIFY (tau_h, tau_s, tau1_h, tau1_s, vtau_h, vtau_s)
      NULLIFY (rho_nlcc)

      CALL get_qs_env(qs_env=qs_env, &
                      input=input, &
                      dft_control=dft_control, &
                      para_env=para_env, &
                      qs_kind_set=qs_kind_set, &
                      atomic_kind_set=atomic_kind_set)

      IF (PRESENT(kind_set_external)) THEN
         my_kind_set => kind_set_external
      ELSE
         my_kind_set => qs_kind_set
      END IF
      nlcc = has_nlcc(my_kind_set)

      accint = dft_control%qs_control%gapw_control%accurate_xcint

      CALL section_vals_val_get(input, "DFT%LSD", l_val=lsd)
      CALL section_vals_val_get(xc_section, "DENSITY_CUTOFF", &
                                r_val=density_cut)
      CALL section_vals_val_get(xc_section, "GRADIENT_CUTOFF", &
                                r_val=gradient_cut)
      CALL section_vals_val_get(xc_section, "TAU_CUTOFF", &
                                r_val=tau_cut)

      my_do_sf = .FALSE.
      IF (PRESENT(do_sf)) my_do_sf = do_sf

      xc_fun_section => section_vals_get_subs_vals(xc_section, &
                                                   "XC_FUNCTIONAL")
      IF (lsd) THEN
         nspins = 2
      ELSE
         nspins = 1
      END IF

      scale_rho = .FALSE.
      IF (PRESENT(do_scale)) THEN
         IF (do_scale) THEN
            ! obscure scaling needed for (some) triplet densities
            scale_rho = .TRUE.
            CPASSERT(PRESENT(do_triplet))
            CPASSERT(do_triplet)
            CPASSERT(nspins == 1)
         END IF
      END IF
      IF (PRESENT(do_triplet)) THEN
         IF (nspins == 1 .AND. do_triplet) lsd = .TRUE.
      END IF

      needs = xc_functionals_get_needs(xc_fun_section, lsd=lsd, &
                                       calc_potential=.TRUE.)
      gradient_functional = needs%drho .OR. needs%drho_spin
      tau_f = (needs%tau .OR. needs%tau_spin)
      IF (.NOT. tau_f) rtau = 0.0_dp

      !  Here starts the loop over all the atoms
      DO ikind = 1, SIZE(atomic_kind_set)

         NULLIFY (atom_list, harmonics, grid_atom)
         CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=atom_list, natom=natom)
         CALL get_qs_kind(my_kind_set(ikind), paw_atom=paw_atom, &
                          harmonics=harmonics, grid_atom=grid_atom)
         CALL get_qs_kind(my_kind_set(ikind), basis_set=basis_1c, basis_type="GAPW_1C")
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

         ! Array dimension: here only one dimensional arrays are used,
         ! i.e. only the first column of deriv_data is read.
         ! The other two dimensions are set to size equal 1.
         bounds(1:2, 1:3) = 1
         bounds(2, 1) = na
         bounds(2, 2) = nr

         CALL xc_dset_create(deriv_set, local_bounds=bounds)
         CALL xc_rho_set_create(rho_set_h, bounds, rho_cutoff=density_cut, &
                                drho_cutoff=gradient_cut, tau_cutoff=tau_cut)
         CALL xc_rho_set_create(rho_set_s, bounds, rho_cutoff=density_cut, &
                                drho_cutoff=gradient_cut, tau_cutoff=tau_cut)
         CALL xc_rho_set_create(rho1_set_h, bounds, rho_cutoff=density_cut, &
                                drho_cutoff=gradient_cut, tau_cutoff=tau_cut)
         CALL xc_rho_set_create(rho1_set_s, bounds, rho_cutoff=density_cut, &
                                drho_cutoff=gradient_cut, tau_cutoff=tau_cut)

         ! allocate the required 3d arrays where to store rho and drho
         IF (nspins == 1 .AND. .NOT. lsd) THEN
            CALL xc_rho_set_atom_update(rho_set_h, needs, 1, bounds)
            CALL xc_rho_set_atom_update(rho1_set_h, needs, 1, bounds)
            CALL xc_rho_set_atom_update(rho_set_s, needs, 1, bounds)
            CALL xc_rho_set_atom_update(rho1_set_s, needs, 1, bounds)
         ELSE
            CALL xc_rho_set_atom_update(rho_set_h, needs, 2, bounds)
            CALL xc_rho_set_atom_update(rho1_set_h, needs, 2, bounds)
            CALL xc_rho_set_atom_update(rho_set_s, needs, 2, bounds)
            CALL xc_rho_set_atom_update(rho1_set_s, needs, 2, bounds)
         END IF

         ALLOCATE (rho_h(1:na, 1:nr, 1:nspins), rho1_h(1:na, 1:nr, 1:nspins), &
                   rho_s(1:na, 1:nr, 1:nspins), rho1_s(1:na, 1:nr, 1:nspins))

         ALLOCATE (vxc_h(1:na, 1:nr, 1:nspins), vxc_s(1:na, 1:nr, 1:nspins))
         vxc_h = 0.0_dp
         vxc_s = 0.0_dp

         IF (tau_f) THEN
            CALL create_tau_basis_cache(tau_basis_cache, grid_atom, basis_1c, harmonics)
            ALLOCATE (tau_h(1:na, 1:nr, 1:nspins), tau1_h(1:na, 1:nr, 1:nspins), &
                      tau_s(1:na, 1:nr, 1:nspins), tau1_s(1:na, 1:nr, 1:nspins))
            ALLOCATE (vtau_h(1:na, 1:nr, 1:nspins), vtau_s(1:na, 1:nr, 1:nspins))
         END IF

         IF (gradient_functional) THEN
            ALLOCATE (drho_h(1:4, 1:na, 1:nr, 1:nspins), drho1_h(1:4, 1:na, 1:nr, 1:nspins), &
                      drho_s(1:4, 1:na, 1:nr, 1:nspins), drho1_s(1:4, 1:na, 1:nr, 1:nspins))
            ALLOCATE (vxg_h(1:3, 1:na, 1:nr, 1:nspins), vxg_s(1:3, 1:na, 1:nr, 1:nspins))
         ELSE
            ALLOCATE (drho_h(1, 1, 1, 1), drho1_h(1, 1, 1, 1), &
                      drho_s(1, 1, 1, 1), drho1_s(1, 1, 1, 1))
            ALLOCATE (vxg_h(1, 1, 1, 1), vxg_s(1, 1, 1, 1))
            rrho = 0.0_dp
         END IF
         vxg_h = 0.0_dp
         vxg_s = 0.0_dp

         ! NLCC: prepare rho and drho of the core charge for this KIND
         donlcc = .FALSE.
         IF (nlcc) THEN
            NULLIFY (rho_nlcc)
            rho_nlcc => my_kind_set(ikind)%nlcc_pot
            IF (ASSOCIATED(rho_nlcc)) donlcc = .TRUE.
         END IF

         ! parallelization
         IF (PRESENT(para_env_ext)) THEN
            local_loop_limit = get_limit(natom, para_env_ext%num_pe, para_env_ext%mepos)
         ELSE
            local_loop_limit = get_limit(natom, para_env%num_pe, para_env%mepos)
         END IF

         DO iatom = local_loop_limit(1), local_loop_limit(2) !1,natom
            atom = atom_list(iatom)

            rho_atom_set(atom)%exc_h = 0.0_dp
            rho_atom_set(atom)%exc_s = 0.0_dp
            rho1_atom_set(atom)%exc_h = 0.0_dp
            rho1_atom_set(atom)%exc_s = 0.0_dp

            rho_atom => rho_atom_set(atom)
            rho1_atom => rho1_atom_set(atom)
            NULLIFY (r_h, r_s, dr_h, dr_s, r_h_d, r_s_d)
            NULLIFY (r1_h, r1_s, dr1_h, dr1_s, r1_h_d, r1_s_d)
            rho_h = 0.0_dp
            rho_s = 0.0_dp
            rho1_h = 0.0_dp
            rho1_s = 0.0_dp
            IF (gradient_functional) THEN
               CALL get_rho_atom(rho_atom=rho_atom, &
                                 rho_rad_h=r_h, rho_rad_s=r_s, &
                                 drho_rad_h=dr_h, drho_rad_s=dr_s, &
                                 rho_rad_h_d=r_h_d, rho_rad_s_d=r_s_d)
               CALL get_rho_atom(rho_atom=rho1_atom, &
                                 rho_rad_h=r1_h, rho_rad_s=r1_s, &
                                 drho_rad_h=dr1_h, drho_rad_s=dr1_s, &
                                 rho_rad_h_d=r1_h_d, rho_rad_s_d=r1_s_d)
               drho_h = 0.0_dp; drho_s = 0.0_dp
               drho1_h = 0.0_dp; drho1_s = 0.0_dp
            ELSE
               CALL get_rho_atom(rho_atom=rho_atom, &
                                 rho_rad_h=r_h, rho_rad_s=r_s)
               CALL get_rho_atom(rho_atom=rho1_atom, &
                                 rho_rad_h=r1_h, rho_rad_s=r1_s)
            END IF

            rtot = 0.0_dp

            DO ir = 1, nr
               CALL calc_rho_angular(grid_atom, harmonics, nspins, gradient_functional, &
                                     ir, r_h, r_s, rho_h, rho_s, dr_h, dr_s, r_h_d, r_s_d, &
                                     drho_h, drho_s)
               IF (donlcc) THEN
                  CALL calc_rho_nlcc(grid_atom, nspins, gradient_functional, &
                                     ir, rho_nlcc(:, 1), rho_h, rho_s, rho_nlcc(:, 2), drho_h, drho_s)
               END IF
               CALL calc_rho_angular(grid_atom, harmonics, nspins, gradient_functional, &
                                     ir, r1_h, r1_s, rho1_h, rho1_s, dr1_h, dr1_s, r1_h_d, r1_s_d, &
                                     drho1_h, drho1_s)
            END DO
            IF (tau_f) THEN
               CALL calc_tau_atom(tau_h, tau_s, rho_atom, tau_basis_cache, nspins)
               CALL calc_tau_atom(tau1_h, tau1_s, rho1_atom, tau_basis_cache, nspins)
            END IF
            IF (scale_rho) THEN
               rho_h = 2.0_dp*rho_h
               rho_s = 2.0_dp*rho_s
               IF (gradient_functional) THEN
                  drho_h = 2.0_dp*drho_h
                  drho_s = 2.0_dp*drho_s
               END IF
               IF (tau_f) THEN
                  tau_h = 2.0_dp*tau_h
                  tau_s = 2.0_dp*tau_s
               END IF
            END IF

            DO ir = 1, nr
               IF (tau_f) THEN
                  CALL fill_rho_set(rho_set_h, lsd, nspins, needs, rho_h, drho_h, tau_h, na, ir)
                  CALL fill_rho_set(rho1_set_h, lsd, nspins, needs, rho1_h, drho1_h, tau1_h, na, ir)
                  CALL fill_rho_set(rho_set_s, lsd, nspins, needs, rho_s, drho_s, tau_s, na, ir)
                  CALL fill_rho_set(rho1_set_s, lsd, nspins, needs, rho1_s, drho1_s, tau1_s, na, ir)
               ELSE IF (gradient_functional) THEN
                  CALL fill_rho_set(rho_set_h, lsd, nspins, needs, rho_h, drho_h, rtau, na, ir)
                  CALL fill_rho_set(rho1_set_h, lsd, nspins, needs, rho1_h, drho1_h, rtau, na, ir)
                  CALL fill_rho_set(rho_set_s, lsd, nspins, needs, rho_s, drho_s, rtau, na, ir)
                  CALL fill_rho_set(rho1_set_s, lsd, nspins, needs, rho1_s, drho1_s, rtau, na, ir)
               ELSE
                  CALL fill_rho_set(rho_set_h, lsd, nspins, needs, rho_h, rrho, rtau, na, ir)
                  CALL fill_rho_set(rho1_set_h, lsd, nspins, needs, rho1_h, rrho, rtau, na, ir)
                  CALL fill_rho_set(rho_set_s, lsd, nspins, needs, rho_s, rrho, rtau, na, ir)
                  CALL fill_rho_set(rho1_set_s, lsd, nspins, needs, rho1_s, rrho, rtau, na, ir)
               END IF
            END DO

            CALL xc_2nd_deriv_of_r(xc_section=xc_section, &
                                   rho_set=rho_set_h, rho1_set=rho1_set_h, &
                                   deriv_set=deriv_set, &
                                   w=weight_h, vxc=vxc_h, vxg=vxg_h, vtau=vtau_h, do_triplet=do_triplet, &
                                   do_sf=my_do_sf)
            CALL xc_2nd_deriv_of_r(xc_section=xc_section, &
                                   rho_set=rho_set_s, rho1_set=rho1_set_s, &
                                   deriv_set=deriv_set, &
                                   w=weight_s, vxc=vxc_s, vxg=vxg_s, vtau=vtau_s, do_triplet=do_triplet, &
                                   do_sf=my_do_sf)

            CALL get_rho_atom(rho_atom=rho1_atom, ga_Vlocal_gb_h=int_hh, ga_Vlocal_gb_s=int_ss)
            IF (gradient_functional) THEN
               CALL gaVxcgb_GC(vxc_h, vxc_s, vxg_h, vxg_s, int_hh, int_ss, &
                               grid_atom, basis_1c, harmonics, nspins)
            ELSE
               CALL gaVxcgb_noGC(vxc_h, vxc_s, int_hh, int_ss, &
                                 grid_atom, basis_1c, harmonics, nspins)
            END IF
            IF (tau_f) THEN
               CALL dgaVtaudgb(vtau_h, vtau_s, int_hh, int_ss, &
                               tau_basis_cache, nspins)
            END IF

            NULLIFY (r_h, r_s, dr_h, dr_s)

         END DO

         ! some cleanup
         DEALLOCATE (rho_h, rho_s, rho1_h, rho1_s, vxc_h, vxc_s)
         DEALLOCATE (drho_h, drho_s, vxg_h, vxg_s)
         DEALLOCATE (drho1_h, drho1_s)
         IF (tau_f) THEN
            DEALLOCATE (tau_h, tau_s, tau1_h, tau1_s)
            DEALLOCATE (vtau_h, vtau_s)
            CALL release_tau_basis_cache(tau_basis_cache)
         END IF

         CALL xc_dset_release(deriv_set)
         CALL xc_rho_set_release(rho_set_h)
         CALL xc_rho_set_release(rho1_set_h)
         CALL xc_rho_set_release(rho_set_s)
         CALL xc_rho_set_release(rho1_set_s)
      END DO

      CALL timestop(handle)

   END SUBROUTINE fxc_atom_calc

END MODULE qs_fxc_atom
