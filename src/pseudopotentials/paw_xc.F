!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Scalar PAW one-center semilocal XC on radial and spherical quadrature grids.
!>        Evaluates AE minus compensated-PS densities, including the corresponding frozen cores.
! **************************************************************************************************
MODULE paw_xc
   USE compensated_arithmetic,          ONLY: twofold_add,&
                                              twofold_divide,&
                                              twofold_multiply,&
                                              twofold_product,&
                                              twofold_real_type,&
                                              twofold_value
   USE input_section_types,             ONLY: section_vals_type
   USE kinds,                           ONLY: dp
   USE lebedev,                         ONLY: get_number_of_lebedev_grid,&
                                              lebedev_grid
   USE mathconstants,                   ONLY: fourpi
   USE orbital_pointers,                ONLY: indco,&
                                              init_orbital_pointers,&
                                              nco,&
                                              ncoset
   USE orbital_transformation_matrices, ONLY: init_spherical_harmonics,&
                                              orbtramat
   USE paw_onecenter,                   ONLY: paw_onecenter_type
   USE paw_types,                       ONLY: paw_xc_type
   USE splines,                         ONLY: spline3ders
   USE upf_types,                       ONLY: atom_upfpot_type,&
                                              upf_radial_weights
   USE xc_atom,                         ONLY: fill_rho_set,&
                                              vxc_of_r_new,&
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
#include "../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'paw_xc'

   PUBLIC :: paw_xc_type, paw_xc_init, paw_xc_evaluate

CONTAINS

! **************************************************************************************************
!> \brief Prepare radial density products, gradients, frozen cores and angular quadrature.
!> \param pot Original scalar PAW UPF data.
!> \param nang Minimum number of angular points; the selected rule also resolves harmonic products.
!> \param grid Immutable per-kind one-center XC data.
! **************************************************************************************************
   SUBROUTINE paw_xc_init(pot, nang, grid)
      TYPE(atom_upfpot_type), INTENT(IN)                 :: pot
      INTEGER, INTENT(IN)                                :: nang
      TYPE(paw_xc_type), INTENT(OUT)                     :: grid

      INTEGER                                            :: a, axis, first, i, j, l, lm, m, nb, nr, &
                                                            ns, pair, powers(3), reduced(3), rule, &
                                                            s, side
      REAL(dp)                                           :: factor
      REAL(dp), ALLOCATABLE                              :: values(:)

      CPASSERT(pot%is_paw .AND. .NOT. pot%has_so .AND. ALLOCATED(pot%paw))
      CPASSERT(ALLOCATED(pot%qijl))
      nr = SIZE(pot%r)
      grid%nr = nr
      grid%nval = pot%paw%cutoff_index
      CPASSERT(grid%nval >= 4 .AND. grid%nval <= nr)
      nb = pot%number_of_proj
      grid%npair = nb*(nb + 1)/2
      grid%lmax = 2*MAXVAL(pot%lbeta)
      rule = MAX(get_number_of_lebedev_grid(n=nang), get_number_of_lebedev_grid(l=2*grid%lmax))
      ALLOCATE (grid%directions, SOURCE=lebedev_grid(rule)%r)
      ALLOCATE (grid%wa, SOURCE=lebedev_grid(rule)%w)
      grid%na = SIZE(grid%wa)
      ALLOCATE (grid%r(nr), grid%wr(nr), &
                grid%harmonics(grid%na, (grid%lmax + 1)**2), &
                grid%angular_gradient(3, grid%na, (grid%lmax + 1)**2))
      grid%r(:) = pot%r
      grid%wr(:) = upf_radial_weights(pot)*pot%r**2
      grid%wa(:) = fourpi*grid%wa
      CALL init_orbital_pointers(grid%lmax)
      CALL init_spherical_harmonics(grid%lmax, -1)
      grid%harmonics = 0.0_dp
      grid%angular_gradient = 0.0_dp
      DO l = 0, grid%lmax
         DO m = 1, 2*l + 1
            lm = l*l + m
            DO s = 1, nco(l)
               powers = indco(:, ncoset(l - 1) + s)
               factor = orbtramat(l)%slm(m, s)
               DO a = 1, grid%na
                  grid%harmonics(a, lm) = grid%harmonics(a, lm) + factor*PRODUCT(grid%directions(:, a)**powers)
                  DO axis = 1, 3
                     IF (powers(axis) == 0) CYCLE
                     reduced = powers
                     reduced(axis) = reduced(axis) - 1
                     grid%angular_gradient(axis, a, lm) = grid%angular_gradient(axis, a, lm) + &
                                                          factor*REAL(powers(axis), dp)*PRODUCT(grid%directions(:, a)**reduced)
                  END DO
               END DO
            END DO
            DO a = 1, grid%na
               ! Subtract the radial derivative of r**l*Y_lm to obtain the tangent gradient.
               grid%angular_gradient(:, a, lm) = grid%angular_gradient(:, a, lm) - &
                                                 REAL(l, dp)*grid%harmonics(a, lm)*grid%directions(:, a)
            END DO
         END DO
      END DO

      ALLOCATE (grid%radial(grid%nval, grid%npair, 0:grid%lmax, 2), &
                grid%radial_gradient(grid%nval, grid%npair, 0:grid%lmax, 2), &
                grid%core(nr, 2), grid%core_gradient(nr, 2))
      grid%radial = 0.0_dp
      grid%radial_gradient = 0.0_dp
      first = 1
      IF (pot%r(1) == 0.0_dp) first = 2
      grid%core = 0.0_dp
      grid%core(:, 1) = pot%paw%rho_core
      IF (ALLOCATED(pot%rho_nlcc)) grid%core(:, 2) = pot%rho_nlcc
      DO side = 1, 2
         CALL spline3ders(pot%r, grid%core(:, side), pot%r, dynew=grid%core_gradient(:, side))
      END DO
      ! After the last nonzero AE core sample, enforce the zero slope of the tabulated tail.
      IF (ALL(grid%core(:, 1) == 0.0_dp)) THEN
         grid%core_gradient(:, 1) = 0.0_dp
      ELSE
         i = nr
         DO WHILE (i > 2)
            IF (ANY(grid%core(i - 1:i, 1) /= 0.0_dp)) EXIT
            i = i - 1
         END DO
         IF (grid%core(i, 1) == 0.0_dp) grid%core_gradient(i:, 1) = 0.0_dp
      END IF
      ! Retain a few physical continuation samples for boundary derivatives, without
      ! differentiating an artificial jump to zero at the augmentation-sphere boundary.
      ns = MIN(nr, grid%nval + 3)
      CPASSERT(ns - first + 1 >= 4)
      ALLOCATE (values(ns))
      DO j = 1, nb
         DO i = 1, j
            pair = j*(j - 1)/2 + i
            DO l = ABS(pot%lbeta(i) - pot%lbeta(j)), pot%lbeta(i) + pot%lbeta(j), 2
               DO side = 1, 2
                  IF (side == 1) THEN
                     values(first:) = pot%paw%ae_wfc(first:ns, i)*pot%paw%ae_wfc(first:ns, j)/pot%r(first:ns)**2
                  ELSE
                     values(first:) = (pot%paw%ps_wfc(first:ns, i)*pot%paw%ps_wfc(first:ns, j) + &
                                       pot%qijl(first:ns, pair, l))/pot%r(first:ns)**2
                  END IF
                  CALL spline3ders(pot%r(first:ns), values(first:), pot%r(:grid%nval), &
                                   ynew=grid%radial(:, pair, l, side), &
                                   dynew=grid%radial_gradient(:, pair, l, side))
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE paw_xc_init

! **************************************************************************************************
!> \brief Contract a real symmetric spin-resolved projector density with the product Gaunt tensor.
!> \param onecenter Magnetic/radial mappings and angular integrals.
!> \param density Projector density matrix, one or two spin channels.
!> \param coefficients Radial-pair and real-harmonic coefficients for each spin.
! **************************************************************************************************
   SUBROUTINE density_coefficients(onecenter, density, coefficients)
      TYPE(paw_onecenter_type), INTENT(IN)               :: onecenter
      REAL(dp), INTENT(IN)                               :: density(:, :, :)
      TYPE(twofold_real_type), INTENT(OUT)               :: coefficients(:, :, :)

      INTEGER                                            :: lm, p, pair, q, spin
      LOGICAL                                            :: assertion_valid
      TYPE(twofold_real_type)                            :: component

      assertion_valid = ALL(SHAPE(density) == [onecenter%nproj, onecenter%nproj, SIZE(coefficients, 3)])
      CPASSERT(assertion_valid)
      assertion_valid = SIZE(density, 3) == 1 .OR. SIZE(density, 3) == 2
      CPASSERT(assertion_valid)
      CPASSERT(SIZE(coefficients, 1) == onecenter%npair)
      assertion_valid = SIZE(coefficients, 2) == SIZE(onecenter%gaunt, 3)
      CPASSERT(assertion_valid)
      coefficients = twofold_real_type()
      DO spin = 1, SIZE(density, 3)
         DO q = 1, onecenter%nproj
            DO p = 1, onecenter%nproj
               pair = onecenter%pair_index(p, q)
               component = twofold_real_type(density(p, q, spin), 0.0_dp)
               DO lm = 1, SIZE(coefficients, 2)
                  coefficients(pair, lm, spin) = twofold_add(coefficients(pair, lm, spin), &
                                                  twofold_multiply(component, twofold_real_type(onecenter%gaunt(p, q, lm), 0.0_dp)))
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE density_coefficients

! **************************************************************************************************
!> \brief Dot product retaining the original linear field and quadrature contractions.
!> \param a Retained first factors.
!> \param b Original FP64 second factors.
!> \return Retained dot product.
! **************************************************************************************************
   FUNCTION component_dot(a, b) RESULT(value)
      TYPE(twofold_real_type), INTENT(IN)                :: a(:)
      REAL(dp), INTENT(IN)                               :: b(:)
      TYPE(twofold_real_type)                            :: value

      INTEGER                                            :: i

      CPASSERT(SIZE(a) == SIZE(b))
      value = twofold_real_type()
      DO i = 1, SIZE(a)
         IF (b(i) == 0.0_dp) CYCLE
         value = twofold_add(value, twofold_multiply(a(i), twofold_real_type(b(i), 0.0_dp)))
      END DO
   END FUNCTION component_dot

! **************************************************************************************************
!> \brief Dot product retaining the original FP64 products.
!> \param a First factors.
!> \param b Second factors.
!> \return Retained dot product.
! **************************************************************************************************
   FUNCTION original_dot(a, b) RESULT(value)
      REAL(dp), INTENT(IN)                               :: a(:), b(:)
      TYPE(twofold_real_type)                            :: value

      INTEGER                                            :: i
      TYPE(twofold_real_type)                            :: term

      CPASSERT(SIZE(a) == SIZE(b))
      value = twofold_real_type()
      DO i = 1, SIZE(a)
         CALL twofold_product(a(i), b(i), term%hi, term%lo)
         value = twofold_add(value, term)
      END DO
   END FUNCTION original_dot

! **************************************************************************************************
!> \brief One original AE or compensated-PS shell, rounded only after each complete field sum.
!> \param grid Original PAW radial and angular fields.
!> \param coefficients Retained Gaunt coefficients of the projector occupation.
!> \param ir Positive radial node.
!> \param side AE (1) or compensated PS (2).
!> \param rho Density including the prescribed frozen core.
!> \param gradient Cartesian density gradient.
! **************************************************************************************************
   SUBROUTINE collocate_coefficients(grid, coefficients, ir, side, rho, gradient)
      TYPE(paw_xc_type), INTENT(IN)                      :: grid
      TYPE(twofold_real_type), INTENT(IN)                :: coefficients(:, :, :)
      INTEGER, INTENT(IN)                                :: ir, side
      REAL(dp), INTENT(OUT)                              :: rho(:, :), gradient(:, :, :)

      INTEGER                                            :: angular, axis, l, lm, nlm, nspin, spin
      LOGICAL                                            :: assertion_valid
      REAL(dp)                                           :: core_scale
      TYPE(twofold_real_type) :: core, field, radial, radius, slope(SIZE(coefficients, 2)), &
         harmonics(SIZE(coefficients, 2))

      nspin = SIZE(coefficients, 3)
      nlm = SIZE(coefficients, 2)
      core_scale = 1.0_dp/REAL(nspin, dp)
      CPASSERT(side >= 1 .AND. side <= 2)
      assertion_valid = ir >= 1 .AND. ir <= grid%nr .AND. grid%r(ir) > 0.0_dp
      CPASSERT(assertion_valid)
      CPASSERT(ALL(SHAPE(rho) == [grid%na, nspin]))
      assertion_valid = ALL(SHAPE(gradient) == [3, grid%na, nspin])
      CPASSERT(assertion_valid)
      radius = twofold_real_type(grid%r(ir), 0.0_dp)
      DO spin = 1, nspin
         harmonics = twofold_real_type()
         slope = twofold_real_type()
         IF (ir <= grid%nval) THEN
            DO lm = 1, nlm
               l = INT(SQRT(REAL(lm - 1, dp)))
               harmonics(lm) = component_dot(coefficients(:, lm, spin), grid%radial(ir, :, l, side))
               slope(lm) = component_dot(coefficients(:, lm, spin), grid%radial_gradient(ir, :, l, side))
            END DO
         END IF
         DO angular = 1, grid%na
            CALL twofold_product(core_scale, grid%core(ir, side), core%hi, core%lo)
            rho(angular, spin) = twofold_value(twofold_add(component_dot(harmonics, grid%harmonics(angular, :)), core))
            CALL twofold_product(core_scale, grid%core_gradient(ir, side), core%hi, core%lo)
            radial = twofold_add(component_dot(slope, grid%harmonics(angular, :)), core)
            DO axis = 1, 3
               field = twofold_add(twofold_multiply(twofold_real_type(grid%directions(axis, angular), 0.0_dp), radial), &
                                   twofold_divide(component_dot(harmonics, grid%angular_gradient(axis, angular, :)), radius))
               gradient(axis, angular, spin) = twofold_value(field)
            END DO
         END DO
      END DO
   END SUBROUTINE collocate_coefficients

! **************************************************************************************************
!> \brief Semilocal PAW XC correction and its full projector-density derivative.
!> \param grid Prepared AE/PS radial and angular data.
!> \param onecenter Matching projector mapping and Gaunt tensor.
!> \param functional XC_FUNCTIONAL section, with activated functionals.
!> \param density Spin-resolved projector density, already summed over k points.
!> \param density_cutoff XC evaluator density cutoff.
!> \param gradient_cutoff XC evaluator gradient cutoff.
!> \param energy AE-minus-compensated-PS XC energy.
!> \param potential Derivative dE/d(density), separately for each spin.
! **************************************************************************************************
   SUBROUTINE paw_xc_evaluate(grid, onecenter, functional, density, density_cutoff, gradient_cutoff, energy, potential)
      TYPE(paw_xc_type), INTENT(IN)                      :: grid
      TYPE(paw_onecenter_type), INTENT(IN)               :: onecenter
      TYPE(section_vals_type), POINTER                   :: functional
      REAL(dp), INTENT(IN)                               :: density(:, :, :), density_cutoff, &
                                                            gradient_cutoff
      REAL(dp), INTENT(OUT)                              :: energy, potential(:, :, :)

      INTEGER                                            :: bounds(2, 3), first, ir, jr, na, nr, &
                                                            nspin, p, pair, q, side, spin
      LOGICAL                                            :: assertion_valid, lsd
      REAL(dp)                                           :: batch_energy, sign
      REAL(dp), ALLOCATABLE                              :: rho(:, :, :), weights(:, :)
      REAL(dp), ALLOCATABLE, TARGET                      :: drho(:, :, :, :), tau(:, :, :)
      REAL(dp), POINTER                                  :: vtau(:, :, :), vxc(:, :, :), &
                                                            vxg(:, :, :, :)
      TYPE(twofold_real_type) :: &
         coefficient_derivative(onecenter%npair, SIZE(onecenter%gaunt, 3), SIZE(density, 3)), &
         coefficients(onecenter%npair, SIZE(onecenter%gaunt, 3), SIZE(density, 3)), energy_sum, &
         potential_sum
      TYPE(xc_derivative_set_type)                       :: derivatives
      TYPE(xc_rho_cflags_type)                           :: needs
      TYPE(xc_rho_set_type)                              :: rho_set

      assertion_valid = grid%npair == onecenter%npair .AND. grid%lmax == onecenter%lmax
      CPASSERT(assertion_valid)
      CPASSERT(ALL(SHAPE(potential) == SHAPE(density)))
      CALL density_coefficients(onecenter, density, coefficients)
      nspin = SIZE(density, 3)
      lsd = nspin == 2
      na = grid%na
      needs = xc_functionals_get_needs(functional, lsd=lsd, calc_potential=.TRUE.)
      IF (needs%tau .OR. needs%tau_spin .OR. needs%laplace_rho .OR. needs%laplace_rho_spin) THEN
         CPABORT("One-center PAW XC supports LDA and GGA functionals")
      END IF
      bounds = 1
      bounds(2, 1) = na
      ! Batch radial shells to amortize the functional setup without changing
      ! their quadrature. Zero-weight padding keeps the final batch rectangular.
      nr = MIN(grid%nr, MAX(1, 2048/na))
      bounds(2, 2) = nr
      ALLOCATE (rho(na, nr, nspin), tau(na, nr, nspin))
      ALLOCATE (drho(4, na, nr, nspin), weights(na, nr), &
                vxc(na, nr, nspin), vxg(3, na, nr, nspin), vtau(na, nr, nspin))
      CALL xc_rho_set_create(rho_set, bounds, rho_cutoff=density_cutoff, drho_cutoff=gradient_cutoff, tau_cutoff=0.0_dp)
      CALL xc_rho_set_atom_update(rho_set, needs, nspin, bounds)
      CALL xc_dset_create(derivatives, local_bounds=bounds)
      energy = 0.0_dp
      coefficient_derivative = twofold_real_type()
      energy_sum = twofold_real_type()
      DO side = 1, 2
         sign = REAL(3 - 2*side, dp)
         DO first = 1, grid%nr, nr
            rho = 0.0_dp
            tau = 0.0_dp
            drho = 0.0_dp
            weights = 0.0_dp
            DO jr = 1, nr
               ir = first + jr - 1
               IF (ir <= grid%nr) THEN
                  IF (grid%r(ir) > 0.0_dp) THEN
                     CALL collocate_coefficients(grid, coefficients, ir, side, rho(:, jr, :), drho(1:3, :, jr, :))
                     drho(4, :, jr, :) = SQRT(SUM(drho(1:3, :, jr, :)**2, DIM=1))
                     weights(:, jr) = grid%wr(ir)*grid%wa
                  END IF
               END IF
               CALL fill_rho_set(rho_set, lsd, nspin, needs, rho, drho, tau, na, jr)
            END DO
            CALL xc_dset_zero_all(derivatives)
            vxc = 0.0_dp
            vxg = 0.0_dp
            vtau = 0.0_dp
            CALL vxc_of_r_new(functional, rho_set, derivatives, 1, needs, weights, lsd, na, nr, &
                              batch_energy, vxc, vxg, vtau)
            energy_sum = twofold_add(energy_sum, &
                                     twofold_multiply(twofold_real_type(sign, 0.0_dp), twofold_real_type(batch_energy, 0.0_dp)))
            CALL accumulate_xc_shells(grid, first, side, vxc, vxg, coefficient_derivative)
         END DO
      END DO
      potential = 0.0_dp
      energy = twofold_value(energy_sum)
      DO spin = 1, nspin
         DO q = 1, onecenter%nproj
            DO p = 1, onecenter%nproj
               pair = onecenter%pair_index(p, q)
               potential_sum = component_dot(coefficient_derivative(pair, :, spin), onecenter%gaunt(p, q, :))
               potential(p, q, spin) = twofold_value(potential_sum)
            END DO
         END DO
      END DO
      CALL xc_dset_release(derivatives)
      CALL xc_rho_set_release(rho_set)
      DEALLOCATE (vxc, vxg, vtau)
   END SUBROUTINE paw_xc_evaluate

! **************************************************************************************************
!> \brief Transpose weighted density and Cartesian-gradient covectors to projector coefficients.
!>        The same fixed quadrature and AE-minus-PS signs apply to the energy and its derivative.
!> \param grid Prepared PAW radial and angular fields.
!> \param first First radial shell in this batch.
!> \param side AE (1) or compensated PS (2).
!> \param vxc Weighted density covectors.
!> \param vxg Weighted Cartesian-gradient covectors.
!> \param coefficient_derivative Accumulated spin-resolved projector-pair/harmonic covectors.
! **************************************************************************************************
   SUBROUTINE accumulate_xc_shells(grid, first, side, vxc, vxg, coefficient_derivative)
      TYPE(paw_xc_type), INTENT(IN)                      :: grid
      INTEGER, INTENT(IN)                                :: first, side
      REAL(dp), INTENT(IN)                               :: vxc(:, :, :), vxg(:, :, :, :)
      TYPE(twofold_real_type), INTENT(INOUT) :: coefficient_derivative(:, :, :)

      INTEGER                                            :: angular, axis, ir, jr, l, lm, spin
      TYPE(twofold_real_type)                            :: av, radius, rv, sign, term, &
                                                            weighted(grid%na)

      sign = twofold_real_type(REAL(3 - 2*side, dp), 0.0_dp)
      DO jr = 1, MIN(SIZE(vxc, 2), grid%nval - first + 1)
         ir = first + jr - 1
         IF (grid%r(ir) == 0.0_dp) CYCLE
         radius = twofold_real_type(grid%r(ir), 0.0_dp)
         DO spin = 1, SIZE(vxc, 3)
            DO lm = 1, SIZE(coefficient_derivative, 2)
               l = INT(SQRT(REAL(lm - 1, dp)))
               av = original_dot(grid%harmonics(:, lm), vxc(:, jr, spin))
               rv = twofold_real_type()
               DO axis = 1, 3
                  av = twofold_add(av, twofold_divide( &
                                   original_dot(grid%angular_gradient(axis, :, lm), vxg(axis, :, jr, spin)), radius))
                  DO angular = 1, grid%na
                     CALL twofold_product(grid%directions(axis, angular), vxg(axis, angular, jr, spin), &
                                          weighted(angular)%hi, weighted(angular)%lo)
                  END DO
                  rv = twofold_add(rv, component_dot(weighted, grid%harmonics(:, lm)))
               END DO
               DO angular = 1, SIZE(coefficient_derivative, 1)
                  term = twofold_add(twofold_multiply(av, twofold_real_type(grid%radial(ir, angular, l, side), 0.0_dp)), &
                                     twofold_multiply(rv, twofold_real_type(grid%radial_gradient(ir, angular, l, side), 0.0_dp)))
                  coefficient_derivative(angular, lm, spin) = twofold_add(coefficient_derivative(angular, lm, spin), &
                                                                          twofold_multiply(sign, term))
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE accumulate_xc_shells

END MODULE paw_xc
