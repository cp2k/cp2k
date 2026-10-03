!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM atom_upf_unittest
   USE atom_basis,                      ONLY: atom_basis_opt
   USE atom_kind_orbitals,              ONLY: sgp_potential_conversion
   USE atom_operators,                  ONLY: atom_ppint_release,&
                                              atom_ppint_setup
   USE atom_sgp,                        ONLY: atom_sgp_potential_type,&
                                              atom_sgp_release,&
                                              fit_core_density_gaussians,&
                                              sgp_construction
   USE atom_types,                      ONLY: &
        CGTO_BASIS, atom_basis_gridrep, atom_basis_type, atom_integrals, atom_potential_type, &
        init_atom_basis, init_atom_basis_default_pp, init_atom_potential, release_atom_basis, &
        release_atom_potential, sgp_pseudo, upf_pseudo
   USE atom_upf,                        ONLY: atom_read_upf,&
                                              atom_release_upf,&
                                              atom_upfpot_type
   USE atom_utils,                      ONLY: atom_core_density,&
                                              atom_density
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: silent_print_level
   USE external_potential_types,        ONLY: allocate_potential,&
                                              deallocate_potential,&
                                              init_potential,&
                                              sgp_potential_type
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE input_constants,                 ONLY: contracted_gto,&
                                              do_gapw_gct,&
                                              do_gapw_log,&
                                              do_rhf_atom
   USE input_cp2k_atom,                 ONLY: create_atom_section
   USE input_section_types,             ONLY: section_release,&
                                              section_type,&
                                              section_vals_create,&
                                              section_vals_get_subs_vals,&
                                              section_vals_release,&
                                              section_vals_type,&
                                              section_vals_val_set
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE libcp2k,                         ONLY: cp2k_finalize,&
                                              cp2k_init
   USE mathconstants,                   ONLY: fourpi,&
                                              pi
   USE orbital_pointers,                ONLY: co,&
                                              init_orbital_pointers,&
                                              nco,&
                                              nso
   USE orbital_transformation_matrices, ONLY: init_spherical_harmonics
   USE qs_grid_atom,                    ONLY: allocate_grid_atom,&
                                              create_grid_atom,&
                                              deallocate_grid_atom,&
                                              grid_atom_type

   IMPLICIT NONE
   LOGICAL :: source_root_argument

   CHARACTER(LEN=default_path_length) :: filename
   INTEGER :: i
   REAL(KIND=dp) :: error(6)
   TYPE(atom_upfpot_type) :: pot
   TYPE(atom_sgp_potential_type) :: sgp

   CALL cp2k_init()
   CALL check_core_derivatives()
   CALL check_projector_addition_theorem()
   CALL check_atomic_potential()
   CALL check_upf_units()
   CALL check_nlcc_fit()
   CALL check_basis_origin()
   CALL check_basis_derivatives()
   CALL check_radial_grid_endpoints()
   CALL check_contracted_basis_reader()
   CALL check_basis_optimization_export()
   CALL check_uspp_import()
   CALL check_uspp_atomic()

   ! The same temporary is reused for successive KIND sections in Quickstep.
   ALLOCATE (pot%rho_nlcc(3))
   pot%rho_nlcc = 1.0_dp
   pot%core_correction = .TRUE.
   pot%econf = [2, 4, 0, 0]
   pot%pname = "previous kind"
   CALL atom_release_upf(pot)
   IF (ALLOCATED(pot%rho_nlcc)) ERROR STOP "UPF release leaked the core density"
   IF (pot%core_correction) ERROR STOP "UPF release retained core-correction metadata"
   IF (ANY(pot%econf /= 0)) ERROR STOP "UPF release retained the previous electron configuration"
   IF (pot%pname /= "") ERROR STOP "UPF release retained the previous potential name"
   CALL atom_release_upf(pot)

   ALLOCATE (sgp%a_nlcc(2), sgp%c_nlcc(2))
   sgp%has_nlcc = .TRUE.
   sgp%n_nlcc = 2
   sgp%is_nonlocal(2) = .TRUE.
   CALL atom_sgp_release(sgp)
   IF (ASSOCIATED(sgp%a_nlcc) .OR. ASSOCIATED(sgp%c_nlcc)) ERROR STOP "SGP release leaked NLCC arrays"
   IF (sgp%has_nlcc .OR. sgp%n_nlcc /= 0) ERROR STOP "SGP release retained NLCC metadata"
   IF (ANY(sgp%is_nonlocal)) ERROR STOP "SGP release retained nonlocal channels"
   CALL atom_sgp_release(sgp)

   ! Optional external UPF files exercise the real reader and Gaussian transformation.
   ! This also permits auditing a complete SSSP library without bundling it in CP2K.
   DO i = 1, COMMAND_ARGUMENT_COUNT()
      CALL GET_COMMAND_ARGUMENT(i, filename)
      ! do_regtest.py passes the CP2K source root before any optional UPF audit files.
      IF (i == 1) THEN
         INQUIRE (FILE=TRIM(filename)//"/tests/UNIT_TESTS", EXIST=source_root_argument)
         IF (source_root_argument) CYCLE
      END IF
      CALL atom_read_upf(pot, TRIM(filename), read_augmentation=.TRUE.)
      IF (ABS(SUM(pot%econf) - pot%zion) > 1.e-10_dp) ERROR STOP "UPF valence count mismatch"
      IF (.NOT. ALL(ieee_is_finite(pot%vlocal))) ERROR STOP "Non-finite UPF local potential"
      IF (ANY(pot%r(2:) <= pot%r(:SIZE(pot%r) - 1))) ERROR STOP "UPF mesh is not increasing"
      IF (ANY(pot%rab <= 0.0_dp)) ERROR STOP "Non-positive UPF mesh weights"
      IF (ALLOCATED(pot%dion)) THEN
         IF (.NOT. ALL(ieee_is_finite(pot%dion))) ERROR STOP "Non-finite UPF projector coefficients"
         IF (.NOT. ALL(ieee_is_finite(pot%beta))) ERROR STOP "Non-finite UPF projectors"
         IF (MAXVAL(ABS(pot%dion - TRANSPOSE(pot%dion))) > 1.e-10_dp) ERROR STOP "Non-Hermitian UPF D matrix"
      END IF
      IF (pot%is_ultrasoft) THEN
         CALL check_uspp_dataset(pot, TRIM(filename))
         CALL atom_release_upf(pot)
         CYCLE
      END IF
      CALL sgp_construction(sgp, upf_pot=pot, error=error)
      IF (.NOT. ALL(ieee_is_finite(error))) ERROR STOP "Non-finite UPF transformation errors"
      IF (pot%core_correction) THEN
         IF (.NOT. ASSOCIATED(sgp%a_nlcc)) ERROR STOP "Missing transformed UPF core density"
         IF (.NOT. ALL(ieee_is_finite(sgp%c_nlcc))) ERROR STOP "Non-finite transformed UPF core density"
      END IF
      WRITE (*, '(A,1X,A,1X,6ES16.7)') "UPF transformation:", TRIM(filename), error
      CALL atom_sgp_release(sgp)
      CALL atom_release_upf(pot)
   END DO

   CALL cp2k_finalize()

CONTAINS

! **************************************************************************************************
!> \brief Check the repaired radial maps at the size that exposes endpoint cancellation.
! **************************************************************************************************
   SUBROUTINE check_radial_grid_endpoints()
      INTEGER, PARAMETER :: n = 6400, types(2) = [do_gapw_gct, do_gapw_log]
      REAL(dp), PARAMETER :: inner(2) = [4.83458313724687745e-18_dp, 6.97482911687147684e-18_dp], &
         outer(2) = [2.06843066219244770e17_dp, 57.5213142095328827_dp]

      INTEGER                                            :: i
      REAL(dp)                                           :: error
      TYPE(grid_atom_type), POINTER                      :: grid

      NULLIFY (grid)
      DO i = 1, SIZE(types)
         CALL allocate_grid_atom(grid)
         CALL create_grid_atom(grid, n, 1, 2, 0, types(i))
         IF (.NOT. ALL(ieee_is_finite(grid%rad)) .OR. .NOT. ALL(ieee_is_finite(grid%wr))) THEN
            ERROR STOP 'Nonfinite radial grid'
         END IF
         IF (ANY(grid%rad <= 0.0_dp) .OR. ANY(grid%wr <= 0.0_dp)) ERROR STOP 'Nonpositive radial grid'
         IF (ANY(grid%rad(:n - 1) <= grid%rad(2:))) ERROR STOP 'Unordered radial points'
         ! Independent analytic moments include the r**2 radial measure.
         error = ABS(SUM(grid%wr*EXP(-grid%rad)) - 2.0_dp)
         error = MAX(error, ABS(SUM(grid%wr*EXP(-grid%rad)/grid%rad) - 1.0_dp))
         error = MAX(error, ABS(SUM(grid%wr*EXP(-grid%rad2)) - SQRT(pi)/4.0_dp))
         IF (error > 1.e-10_dp) ERROR STOP 'Radial moment error'
         ! Independent high-precision endpoints reject clipping rounded cosines.
         IF (ABS(grid%rad(n)/inner(i) - 1.0_dp) > 1.e-13_dp .OR. &
             ABS(grid%rad(1)/outer(i) - 1.0_dp) > 1.e-13_dp) ERROR STOP 'Radial grid endpoint error'
         CALL deallocate_grid_atom(grid)
      END DO
   END SUBROUTINE check_radial_grid_endpoints

! **************************************************************************************************
!> \brief Check primitive and contracted Gaussian derivatives on UPF meshes starting at zero.
! **************************************************************************************************
   SUBROUTINE check_basis_origin()
      INTEGER                                            :: i, l, n, variant
      REAL(KIND=dp)                                      :: expected(3), r(6001), rab(6001), scale
      TYPE(atom_basis_type)                              :: basis, radial

      CALL init_atom_basis_default_pp(basis)
      n = MAXVAL(basis%nbas)
      ! A long supplied mesh must not generate an unused singular auxiliary grid.
      r = [(0.002_dp*REAL(i - 1, dp), i=1, SIZE(r))]
      rab = 0.002_dp
      CALL atom_basis_gridrep(basis, radial, r, rab)
      IF (ANY(radial%grid%rad /= r) .OR. ANY(radial%grid%wr /= rab*r**2)) THEN
         ERROR STOP "Supplied UPF mesh changed during basis evaluation"
      END IF
      IF (MAXVAL(ABS(radial%bf(:, 1, 0) - EXP(-basis%am(1, 0)*r**2))) > 1.e-14_dp) THEN
         ERROR STOP "Incorrect basis values on long UPF mesh"
      END IF
      CALL release_atom_basis(radial)
      DO variant = 1, 2
         scale = REAL(variant, dp)
         IF (variant == 2) THEN
            basis%basis_type = CGTO_BASIS
            ALLOCATE (basis%cm(n, n, 0:UBOUND(basis%am, 2)))
            basis%cm = 0.0_dp
            DO i = 1, n
               basis%cm(i, i, :) = scale
            END DO
         END IF
         CALL atom_basis_gridrep(basis, radial, [0.0_dp, 1.e-16_dp, 0.1_dp], [0.1_dp, 0.1_dp, 0.1_dp])
         IF (.NOT. ALL(ieee_is_finite(radial%bf)) .OR. .NOT. ALL(ieee_is_finite(radial%dbf)) .OR. &
             .NOT. ALL(ieee_is_finite(radial%ddbf))) ERROR STOP "Nonfinite Gaussian derivatives at the UPF origin"
         DO l = 0, UBOUND(basis%am, 2)
            DO i = 1, basis%nbas(l)
               expected = 0.0_dp
               IF (l == 0) expected = [1.0_dp, 0.0_dp, -2.0_dp*basis%am(i, l)]
               IF (l == 1) expected(2) = 1.0_dp
               IF (l == 2) expected(3) = 2.0_dp
               IF (MAXVAL(ABS([radial%bf(1, i, l), radial%dbf(1, i, l), radial%ddbf(1, i, l)] - &
                              scale*expected)) > 1.e-12_dp) ERROR STOP "Incorrect Gaussian origin limit"
            END DO
         END DO
         CALL release_atom_basis(radial)
      END DO
      CALL release_atom_basis(basis)
   END SUBROUTINE check_basis_origin

! **************************************************************************************************
!> \brief Differentiate radial basis values and densities independently at finite radii.
! **************************************************************************************************
   SUBROUTINE check_basis_derivatives()
      REAL(dp), PARAMETER :: radius(4) = [0.08_dp, 0.3_dp, 0.9_dp, 2.1_dp], step = 2.e-4_dp

      INTEGER                                            :: i, k, l, n, variant
      REAL(dp)                                           :: derivative(5), errors(3), first, lap(5), &
                                                            r(5), reference, rho(5), second
      REAL(dp), ALLOCATABLE                              :: pmat(:, :, :)
      TYPE(atom_basis_type)                              :: basis, radial

      CALL init_atom_basis_default_pp(basis)
      n = MAXVAL(basis%nbas)
      ! Nonunit exponents and signed mixing exercise primitives and genuine contractions.
      basis%nbas = 2
      basis%nprim = 2
      basis%am(1, :) = 0.37_dp
      basis%am(2, :) = 2.3_dp
      ALLOCATE (pmat(n, n, 0:UBOUND(basis%am, 2)))
      pmat = 0.0_dp
      pmat(1, 1, :) = 0.8_dp
      pmat(2, 2, :) = 0.4_dp
      pmat(1, 2, :) = -0.17_dp
      pmat(2, 1, :) = -0.17_dp
      errors = 0.0_dp
      DO variant = 1, 2
         IF (variant == 2) THEN
            basis%basis_type = CGTO_BASIS
            ALLOCATE (basis%cm(n, n, 0:UBOUND(basis%am, 2)))
            basis%cm = 0.0_dp
            basis%cm(1, 1, :) = 0.7_dp
            basis%cm(2, 1, :) = -0.2_dp
            basis%cm(1, 2, :) = 0.1_dp
            basis%cm(2, 2, :) = 1.3_dp
         END IF
         DO k = 1, SIZE(radius)
            r = radius(k) + step*[-2.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp]
            CALL atom_basis_gridrep(basis, radial, r, SPREAD(step, 1, 5))
            DO l = 0, UBOUND(basis%am, 2)
               DO i = 1, 2
                  first = DOT_PRODUCT(radial%bf(:, i, l), [1.0_dp, -8.0_dp, 0.0_dp, 8.0_dp, -1.0_dp])/(12*step)
                  second = DOT_PRODUCT(radial%bf(:, i, l), [-1.0_dp, 16.0_dp, -30.0_dp, 16.0_dp, -1.0_dp])/(12*step**2)
                  errors(1) = MAX(errors(1), ABS(first - radial%dbf(3, i, l))/(1 + ABS(first)))
                  errors(2) = MAX(errors(2), ABS(second - radial%ddbf(3, i, l))/(1 + ABS(second)))
               END DO
            END DO
            CALL atom_density(rho, pmat, radial, UBOUND(basis%am, 2), typ="RHO")
            CALL atom_density(derivative, pmat, radial, UBOUND(basis%am, 2), typ="DER")
            CALL atom_density(lap, pmat, radial, UBOUND(basis%am, 2), typ="LAP", rr=r*r)
            first = DOT_PRODUCT(rho, [1.0_dp, -8.0_dp, 0.0_dp, 8.0_dp, -1.0_dp])/(12*step)
            second = DOT_PRODUCT(rho, [-1.0_dp, 16.0_dp, -30.0_dp, 16.0_dp, -1.0_dp])/(12*step**2)
            errors(1) = MAX(errors(1), ABS(first - derivative(3))/(1 + ABS(first)))
            reference = second + 2*first/radius(k)
            errors(3) = MAX(errors(3), ABS(reference - lap(3))/(1 + ABS(reference)))
            CALL release_atom_basis(radial)
         END DO
      END DO
      WRITE (*, '(A,3ES24.15)') "Atomic radial derivative/Laplacian scaled errors: ", errors
      DEALLOCATE (pmat)
      CALL release_atom_basis(basis)
      IF (MAXVAL(errors) > 2.e-7_dp) ERROR STOP "Atomic radial finite-difference derivatives disagree"
   END SUBROUTINE check_basis_derivatives

! **************************************************************************************************
!> \brief Check fitted densities and their derivatives against an analytic reference.
!>        Evaluation points differ from the input mesh, and neither input mesh
!>        nor reference exponents coincide with the internal fitting sequence.
! **************************************************************************************************
   SUBROUTINE check_nlcc_fit()
      INTEGER, PARAMETER                                 :: nr = 1601
      REAL(dp), PARAMETER :: amplitudes(3) = [0.02_dp, 1.7_dp, 0.4_dp], &
         exponents(3) = [0.15_dp, 2.3_dp, 35.0_dp]

      INTEGER                                            :: i, mesh
      REAL(dp) :: density_error, derivative, derivative_error, dx, r(nr), rab(nr), radius, &
         reference, reference_derivative, rho(nr), value
      REAL(dp), ALLOCATABLE                              :: alpha(:), coeff(:)

      DO mesh = 1, 2
         dx = LOG(20.0_dp/1.e-6_dp)/REAL(nr - 1, dp)
         DO i = 1, nr
            IF (mesh == 1) THEN
               r(i) = 20.0_dp*REAL(i - 1, dp)/REAL(nr - 1, dp)
               rab(i) = 20.0_dp/REAL(nr - 1, dp)
            ELSE
               r(i) = 1.e-6_dp*EXP(dx*REAL(i - 1, dp))
               rab(i) = dx*r(i)
            END IF
            rho(i) = SUM(amplitudes*EXP(-exponents*r(i)**2))
         END DO
         CALL fit_core_density_gaussians(r, rab, rho, alpha, coeff)
         density_error = 0.0_dp
         derivative_error = 0.0_dp
         DO i = 0, 1200
            radius = REAL(i, dp)*0.017_dp
            reference = SUM(amplitudes*EXP(-exponents*radius**2))
            reference_derivative = -2.0_dp*radius*SUM(exponents*amplitudes*EXP(-exponents*radius**2))
            value = SUM(coeff*EXP(-alpha*radius**2))
            derivative = -2.0_dp*radius*SUM(alpha*coeff*EXP(-alpha*radius**2))
            density_error = MAX(density_error, ABS(value - reference))
            derivative_error = MAX(derivative_error, ABS(derivative - reference_derivative))
         END DO
         WRITE (*, '(A,I2,2ES16.7)') "Core Gaussian fit, density and derivative:", mesh, density_error, derivative_error
         IF (density_error > 2.e-7_dp) ERROR STOP "Gaussian core-density fit is inaccurate"
         IF (derivative_error > 3.e-6_dp) ERROR STOP "Gaussian core-density derivative is inaccurate"
         rho = 0.0_dp
         CALL fit_core_density_gaussians(r, rab, rho, alpha, coeff)
         IF (ANY(coeff /= 0.0_dp)) ERROR STOP "Zero core density generated nonzero Gaussian coefficients"
      END DO
   END SUBROUTINE check_nlcc_fit

! **************************************************************************************************
!> \brief Read an independent miniature UPF fixture and check its numerical convention.
!>        The values are deliberately simple and are not a physical pseudopotential.
!>        Checking beta separately matters: an inverse rescaling of D leaves the NC
!>        operator unchanged but would require rescaling USPP augmentation data too.
! **************************************************************************************************
   SUBROUTINE check_upf_units()
      CHARACTER(LEN=*), PARAMETER :: fixture = &
         'upf-units-test-with-a-potential-filename-longer-than-'// &
         'the-default-eighty-character-string-buffer.upf'

      INTEGER                                            :: unit
      TYPE(atom_potential_type)                          :: atomic
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_type), POINTER                        :: specification
      TYPE(section_vals_type), POINTER                   :: input, potential_section

      logger => cp_get_default_logger()
      IF (logger%para_env%is_source()) THEN
         OPEN (NEWUNIT=unit, FILE=fixture, STATUS="REPLACE", ACTION="WRITE")
         WRITE (unit, '(A)') &
            '<UPF version="2.0.1">', &
            '<PP_HEADER element="He" pseudo_type="NC" z_valence="2"', &
            ' mesh_size="3" number_of_proj="2" l_max="0" />', &
            '<PP_MESH mesh="3">', &
            '<PP_R type="real" size="3">', &
            '0.5 1.0 1.5', &
            '</PP_R>', &
            '<PP_RAB type="real" size="3">', &
            '0.5 0.5 0.5', &
            '</PP_RAB>', &
            '</PP_MESH>', &
            '<PP_LOCAL type="real" size="3">', &
            '6.0 4.0 2.0', &
            '</PP_LOCAL>', &
            '<PP_NONLOCAL>', &
            '<PP_BETA.1 type="real" size="3" index="1" angular_momentum="0">', &
            '1.0 4.0 9.0', &
            '</PP_BETA.1>', &
            '<PP_BETA.2 type="real" size="3" index="2" angular_momentum="0">', &
            '2.0 -1.0 3.0', &
            '</PP_BETA.2>', &
            '<PP_DIJ type="real" size="4">', &
            '2.0 0.4 0.4 8.0', &
            '</PP_DIJ>', &
            '</PP_NONLOCAL>', &
            '</UPF>'
         CLOSE (unit)
      END IF
      CALL logger%para_env%sync()
      NULLIFY (specification, input)
      CALL create_atom_section(specification)
      CALL section_vals_create(input, specification)
      potential_section => section_vals_get_subs_vals(input, 'POTENTIAL')
      CALL section_vals_val_set(potential_section, 'PSEUDO_TYPE', i_val=upf_pseudo)
      CALL section_vals_val_set(potential_section, 'POTENTIAL_FILE_NAME', c_val=fixture)
      CALL init_atom_potential(atomic, potential_section, 2)
      ASSOCIATE (potential => atomic%upf_pot)
         IF (MAXVAL(ABS(potential%vlocal - [3.0_dp, 2.0_dp, 1.0_dp])) > 1.e-14_dp) THEN
            ERROR STOP "UPF local potential is not in Hartree"
         END IF
         IF (MAXVAL(ABS(potential%beta(:, 1) - [2.0_dp, 4.0_dp, 6.0_dp])) > 1.e-14_dp) THEN
            ERROR STOP "UPF changed the first tabulated projector amplitude"
         END IF
         IF (MAXVAL(ABS(potential%beta(:, 2) - [4.0_dp, -1.0_dp, 2.0_dp])) > 1.e-14_dp) THEN
            ERROR STOP "UPF changed the second tabulated projector amplitude"
         END IF
         IF (MAXVAL(ABS(potential%dion - RESHAPE([1.0_dp, 0.2_dp, 0.2_dp, 4.0_dp], [2, 2]))) > 1.e-14_dp) THEN
            ERROR STOP "UPF D matrix is not in the tabulated-projector convention"
         END IF
      END ASSOCIATE
      CALL release_atom_potential(atomic)
      CALL section_vals_release(input)
      CALL section_release(specification)
      CALL logger%para_env%sync()
      IF (logger%para_env%is_source()) THEN
         OPEN (NEWUNIT=unit, FILE=fixture, STATUS="OLD")
         CLOSE (unit, STATUS="DELETE")
      END IF
   END SUBROUTINE check_upf_units

! **************************************************************************************************
!> \brief Compare the converted atomic SGP operator with analytic Gaussian integrals.
!>        This checks that atomic guesses and the Quickstep Hamiltonian describe
!>        the same radial projectors and erf-regularized local Coulomb potential.
! **************************************************************************************************
   SUBROUTINE check_atomic_potential()
      INTEGER                                            :: l
      REAL(KIND=dp)                                      :: a, b, p, reference
      TYPE(atom_basis_type)                              :: basis
      TYPE(atom_integrals)                               :: integrals
      TYPE(atom_potential_type)                          :: atomic
      TYPE(sgp_potential_type), POINTER                  :: potential

      NULLIFY (potential)
      CALL allocate_potential(potential)
      potential%zeff = 1.0_dp
      potential%core_charge_radius = 0.4_dp
      potential%alpha_core_charge = 0.5_dp/potential%core_charge_radius**2
      potential%n_nonlocal = 1
      potential%lmax = 2
      potential%is_nonlocal(0:2) = .TRUE.
      potential%n_local = 1
      ALLOCATE (potential%elec_conf(0:3), potential%a_nonlocal(1), &
                potential%c_nonlocal(1, 1, 0:2), potential%h_nonlocal(1, 0:2), &
                potential%a_local(1), potential%c_local(1))
      potential%elec_conf = [1, 0, 0, 0]
      potential%a_nonlocal = 0.7_dp
      potential%c_nonlocal = 1.0_dp/SQRT(fourpi)
      potential%h_nonlocal = 1.0_dp
      potential%a_local = 1.2_dp
      potential%c_local = 0.3_dp
      atomic%ppot_type = sgp_pseudo
      CALL sgp_potential_conversion(potential, atomic%sgp_pot)
      CALL init_atom_basis_default_pp(basis)
      CALL atom_ppint_setup(integrals, basis, atomic)
      DO l = 0, 2
         a = basis%am(8, l) + potential%a_nonlocal(1)
         reference = (GAMMA(l + 1.5_dp)/(2.0_dp*a**(l + 1.5_dp)))**2
         IF (ABS(integrals%hnl(8, 8, l) - reference) > 1.e-9_dp*ABS(reference)) THEN
            ERROR STOP "Atomic SGP nonlocal operator differs from analytic Gaussian integrals"
         END IF
      END DO
      p = 2.0_dp*basis%am(8, 0)
      b = 1.0_dp/(SQRT(2.0_dp)*potential%core_charge_radius)
      reference = -potential%zeff*b/(2.0_dp*p*SQRT(p + b*b)) + &
                  potential%c_local(1)*GAMMA(1.5_dp)/(2.0_dp*(p + potential%a_local(1))**1.5_dp)
      IF (ABS(integrals%core(8, 8, 0) - reference) > 1.e-9_dp*ABS(reference)) THEN
         ERROR STOP "Atomic SGP local operator differs from analytic Gaussian integrals"
      END IF
      CALL atom_ppint_release(integrals)
      CALL release_atom_basis(basis)
      CALL deallocate_potential(potential)
   END SUBROUTINE check_atomic_potential

! **************************************************************************************************
!> \brief Verify the physical angular kernel of Gaussian UPF projectors through l=5.
!>        For beta_l(r) = r**l exp(-a*r**2), stripping off the radial exponentials
!>        at unit vectors leaves sum_m Y_lm(u) Y_lm(v) = (2*l+1) P_l(u.v)/(4*pi).
!>        The reference uses the Legendre recurrence, independently of the Cartesian
!>        to spherical transformation used by the production projector integrals.
! **************************************************************************************************
   SUBROUTINE check_projector_addition_theorem()
      INTEGER, PARAMETER                                 :: lmax = 5

      INTEGER                                            :: i, j, l, lx, ly, lz, nc, ns
      REAL(KIND=dp)                                      :: cross, diagonal, p(0:lmax), u(3), v(3), x
      REAL(KIND=dp), ALLOCATABLE                         :: pu(:), pv(:)
      TYPE(sgp_potential_type), POINTER                  :: potential

      NULLIFY (potential)
      CALL init_orbital_pointers(lmax)
      CALL init_spherical_harmonics(lmax, -1)
      CALL allocate_potential(potential)
      potential%n_nonlocal = 1
      potential%lmax = lmax
      potential%nppnl = SUM(nco(0:lmax))
      ALLOCATE (potential%c_nonlocal(1, 1, 0:lmax), potential%h_nonlocal(1, 0:lmax))
      ! atom_sgp stores the radial coefficients divided by sqrt(4*pi).
      potential%c_nonlocal = 1.0_dp/SQRT(fourpi)
      potential%h_nonlocal = 1.0_dp
      CALL init_potential(potential)
      u = [0.6_dp, 0.0_dp, 0.8_dp]
      v = [0.3_dp, -0.4_dp, 0.5_dp]
      v = v/NORM2(v)
      x = DOT_PRODUCT(u, v)
      p(0) = 1.0_dp
      p(1) = x
      DO l = 2, lmax
         p(l) = ((2*l - 1)*x*p(l - 1) - (l - 1)*p(l - 2))/l
      END DO
      nc = 0
      ns = 0
      DO l = 0, lmax
         ALLOCATE (pu(nso(l)), pv(nso(l)))
         pu = 0.0_dp
         pv = 0.0_dp
         DO lx = 0, l
            DO ly = 0, l - lx
               lz = l - lx - ly
               i = nc + co(lx, ly, lz)
               DO j = 1, nso(l)
                  pu(j) = pu(j) + potential%cprj_ppnl(i, ns + j)*u(1)**lx*u(2)**ly*u(3)**lz
                  pv(j) = pv(j) + potential%cprj_ppnl(i, ns + j)*v(1)**lx*v(2)**ly*v(3)**lz
               END DO
            END DO
         END DO
         diagonal = REAL(2*l + 1, dp)/fourpi
         cross = diagonal*p(l)
         IF (ABS(DOT_PRODUCT(pu, pu) - diagonal) > 1.e-12_dp .OR. &
             ABS(DOT_PRODUCT(pv, pv) - diagonal) > 1.e-12_dp .OR. &
             ABS(DOT_PRODUCT(pu, pv) - cross) > 1.e-12_dp) THEN
            ERROR STOP "SGP projector kernel violates the spherical-harmonic addition theorem"
         END IF
         nc = nc + nco(l)
         ns = ns + nso(l)
         DEALLOCATE (pu, pv)
      END DO
      CALL deallocate_potential(potential)
   END SUBROUTINE check_projector_addition_theorem

! **************************************************************************************************
!> \brief Verify radial derivatives and the spherical Laplacian against finite differences.
! **************************************************************************************************
   SUBROUTINE check_core_derivatives()
      REAL(KIND=dp), PARAMETER                           :: h = 1.e-4_dp

      REAL(KIND=dp)                                      :: deriv(4), lap(4), lap_fd(4), minus(4), &
                                                            plus(4), r(4), rho(4)
      TYPE(atom_potential_type)                          :: potential

      potential%ppot_type = sgp_pseudo
      potential%sgp_pot%has_nlcc = .TRUE.
      potential%sgp_pot%n_nlcc = 2
      potential%sgp_pot%a_nlcc(1:2) = [0.7_dp, 2.3_dp]
      potential%sgp_pot%c_nlcc(1:2) = [1.8_dp, -0.2_dp]
      r = [0.0_dp, 0.1_dp, 0.7_dp, 2.0_dp]
      rho = 0.0_dp
      plus = 0.0_dp
      minus = 0.0_dp
      deriv = 0.0_dp
      lap = 0.0_dp
      CALL atom_core_density(rho, potential, "RHO", r)
      CALL atom_core_density(plus, potential, "RHO", r + h)
      CALL atom_core_density(minus, potential, "RHO", r - h)
      CALL atom_core_density(deriv, potential, "DER", r)
      CALL atom_core_density(lap, potential, "LAP", r)
      IF (MAXVAL(ABS(deriv - (plus - minus)/(2.0_dp*h))) > 1.e-7_dp) THEN
         ERROR STOP "SGP core-density derivative differs from finite differences"
      END IF
      lap_fd = (plus - 2.0_dp*rho + minus)/h**2
      lap_fd(1) = 3.0_dp*lap_fd(1)
      lap_fd(2:) = lap_fd(2:) + (plus(2:) - minus(2:))/(h*r(2:))
      IF (MAXVAL(ABS(lap - lap_fd)) > 1.e-6_dp) ERROR STOP "SGP core-density Laplacian differs from finite differences"
   END SUBROUTINE check_core_derivatives

! **************************************************************************************************
!> \brief Check independent shell, primitive and set counts in file and inline atomic bases.
! **************************************************************************************************
   SUBROUTINE check_contracted_basis_reader()
      INTEGER, PARAMETER :: contractions(3) = [11, 1, 1], lower_l(3) = [0, 1, 2], &
         primitives(3) = [11, 41, 2], sets(3) = [1, 1, 21], upper_l(3) = [3, 6, 3]
      CHARACTER(LEN=1600), ALLOCATABLE                   :: lines(:)
      INTEGER                                            :: b, c, i, j, l, mode, n, nc, np, ns, row, &
                                                            s, t, unit
      REAL(dp)                                           :: alpha, coefficient, expected, &
                                                            normalization
      TYPE(atom_basis_type)                              :: basis
      CHARACTER(LEN=*), PARAMETER :: fixture = "atom-contracted-basis-test"

      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_type), POINTER                        :: specification
      TYPE(section_vals_type), POINTER                   :: input, pp_basis

      logger => cp_get_default_logger()
      DO b = 1, SIZE(sets)
         np = primitives(b)
         nc = contractions(b)
         ns = sets(b)
         ALLOCATE (lines(1 + ns*(np + 1)))
         WRITE (lines(1), '(I0)') ns
         row = 1
         DO s = 1, ns
            row = row + 1
            WRITE (lines(row), '(*(I0,1X))') 1, lower_l(b), upper_l(b), np, (nc, l=lower_l(b), upper_l(b))
            DO i = 1, np
               row = row + 1
               alpha = REAL(10*s + i, dp)/20.0_dp
               WRITE (lines(row), '(*(F12.6,1X))') alpha, &
                  ((REAL(100*s + 10*i + c + 5*l, dp)/1000.0_dp, c=1, nc), l=lower_l(b), upper_l(b))
            END DO
         END DO
         IF (logger%para_env%is_source()) THEN
            OPEN (NEWUNIT=unit, FILE=fixture, STATUS="REPLACE", ACTION="WRITE")
            WRITE (unit, '(A)') 'H READER-TEST', (TRIM(lines(i)), i=1, SIZE(lines))
            CLOSE (unit)
         END IF
         CALL logger%para_env%sync()
         DO mode = 1, 2
            NULLIFY (specification, input)
            CALL create_atom_section(specification)
            CALL section_vals_create(input, specification)
            pp_basis => section_vals_get_subs_vals(input, 'PP_BASIS')
            CALL section_vals_val_set(pp_basis, 'BASIS_TYPE', i_val=contracted_gto)
            CALL section_vals_val_set(pp_basis, 'GRID_POINTS', i_val=32)
            CALL section_vals_val_set(pp_basis, 'BASIS_SET', c_val='READER-TEST')
            CALL section_vals_val_set(pp_basis, 'BASIS_SET_FILE_NAME', c_val=fixture)
            IF (mode == 2) THEN
               DO i = 1, SIZE(lines)
                  CALL section_vals_val_set(pp_basis, 'BASIS%_DEFAULT_KEYWORD_', i_rep_val=i, c_val=TRIM(lines(i)))
               END DO
            END IF
            CALL init_atom_basis(basis, pp_basis, 1, 'PP')
            DO l = 0, UBOUND(basis%nbas, 1)
               n = 0
               IF (l >= lower_l(b) .AND. l <= upper_l(b)) n = ns
               IF (basis%nprim(l) /= n*np .OR. basis%nbas(l) /= n*nc) THEN
                  ERROR STOP "Incorrect atomic contracted basis dimensions"
               END IF
               DO s = 1, n
                  DO i = 1, np
                     alpha = REAL(10*s + i, dp)/20.0_dp
                     j = (s - 1)*np + i
                     IF (ABS(basis%am(j, l) - alpha) > 1.e-14_dp) ERROR STOP "Incorrect contracted basis exponent"
                     ! Inverse square root of the analytic radial Gaussian norm.
                     normalization = SQRT(2.0_dp*(2.0_dp*alpha)**(l + 1.5_dp)/GAMMA(l + 1.5_dp))
                     DO t = 1, n
                        DO c = 1, nc
                           coefficient = 0.0_dp
                           IF (s == t) coefficient = REAL(100*s + 10*i + c + 5*l, dp)/1000.0_dp
                           expected = normalization*coefficient
                           IF (ABS(basis%cm(j, (t - 1)*nc + c, l) - expected) > 1.e-12_dp*MAX(1.0_dp, ABS(expected))) THEN
                              ERROR STOP "Incorrect contracted basis coefficient or normalization"
                           END IF
                        END DO
                     END DO
                  END DO
               END DO
            END DO
            CALL release_atom_basis(basis)
            CALL section_vals_release(input)
            CALL section_release(specification)
         END DO
         CALL logger%para_env%sync()
         IF (logger%para_env%is_source()) THEN
            OPEN (NEWUNIT=unit, FILE=fixture, STATUS="OLD")
            CLOSE (unit, STATUS="DELETE")
         END IF
         DEALLOCATE (lines)
      END DO
   END SUBROUTINE check_contracted_basis_reader

! **************************************************************************************************
!> \brief Verify exported orbitals use the best basis restored after Powell optimization.
! **************************************************************************************************
   SUBROUTINE check_basis_optimization_export()
      CHARACTER(LEN=160)                                 :: line
      INTEGER                                            :: i, j, status, unit
      REAL(dp)                                           :: alpha(2), coefficient(2), norm
      REAL(dp), POINTER                                  :: exponents(:)
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_type), POINTER                        :: specification
      TYPE(section_vals_type), POINTER                   :: input

      NULLIFY (specification, input)
      CALL create_atom_section(specification)
      CALL section_vals_create(input, specification)
      CALL section_vals_val_set(input, 'ELEMENT', c_val='He')
      CALL section_vals_val_set(input, 'ELECTRON_CONFIGURATION', c_val='1s2')
      CALL section_vals_val_set(input, 'MAX_ANGULAR_MOMENTUM', i_val=0)
      CALL section_vals_val_set(input, 'AE_BASIS%NUM_GTO', i_val=2)
      ALLOCATE (exponents(2))
      exponents = [2.0_dp, 0.4_dp]
      CALL section_vals_val_set(input, 'AE_BASIS%S_EXPONENTS', r_vals_ptr=exponents)
      CALL section_vals_val_set(input, 'METHOD%METHOD_TYPE', i_val=do_rhf_atom)
      CALL section_vals_val_set(input, 'OPTIMIZATION%EPS_SCF', r_val=1.e-10_dp)
      CALL section_vals_val_set(input, 'OPTIMIZATION%N_DIIS', i_val=8)
      CALL section_vals_val_set(input, 'OPTIMIZATION%MAX_ITER', i_val=100)
      CALL section_vals_val_set(input, 'POWELL%MAX_FUN', i_val=5)
      CALL section_vals_val_set(input, 'POWELL%STEP_SIZE', r_val=0.15_dp)
      CALL section_vals_val_set(input, 'POWELL%ACCURACY', r_val=1.e-5_dp)
      CALL section_vals_val_set(input, 'PRINT%FIT_BASIS%_SECTION_PARAMETERS_', i_val=silent_print_level)
      CALL atom_basis_opt(input)
      CALL section_vals_release(input)
      CALL section_release(specification)

      logger => cp_get_default_logger()
      CALL logger%para_env%sync()
      OPEN (NEWUNIT=unit, FILE='OPT_BASIS', STATUS='OLD', ACTION='READ')
      DO
         READ (unit, '(A)', IOSTAT=status) line
         IF (status /= 0) ERROR STOP 'Missing optimized basis orbital coefficients'
         IF (INDEX(line, 'L Quantum Number:') > 0) EXIT
      END DO
      DO i = 1, 2
         READ (unit, *) alpha(i), coefficient(i)
      END DO
      CLOSE (unit)
      IF (alpha(1) <= 2.0_dp) ERROR STOP 'Optimization did not restore the earlier best exponent'
      norm = 0.0_dp
      DO i = 1, 2
         DO j = 1, 2
            norm = norm + coefficient(i)*coefficient(j)* &
                   (2.0_dp*SQRT(alpha(i)*alpha(j))/(alpha(i) + alpha(j)))**1.5_dp
         END DO
      END DO
      ! The export rounds its coefficients to eight decimal places.
      IF (.NOT. ieee_is_finite(norm) .OR. ABS(norm - 1.0_dp) > 1.e-7_dp) THEN
         ERROR STOP 'Exported orbital is inconsistent with the optimized basis'
      END IF
      CALL logger%para_env%sync()
      IF (logger%para_env%is_source()) THEN
         OPEN (NEWUNIT=unit, FILE='OPT_BASIS', STATUS='OLD')
         CLOSE (unit, STATUS='DELETE')
      END IF
   END SUBROUTINE check_basis_optimization_export

! **************************************************************************************************
!> \brief Analytic Gaussian reference for the atomic USPP metric, potential and adjoint.
! **************************************************************************************************
   SUBROUTINE check_uspp_atomic()
      USE atom_operators, ONLY: atom_int_release, &
                                atom_int_setup, &
                                atom_ppint_release, &
                                atom_ppint_setup
      USE atom_types, ONLY: atom_basis_gridrep, &
                            atom_basis_type, &
                            atom_integrals, &
                            atom_potential_type, &
                            init_atom_basis_default_pp, &
                            lmat, &
                            release_atom_basis, &
                            upf_pseudo
      USE atom_upf, ONLY: atom_release_upf
      USE atom_utils, ONLY: atom_density, &
                            numpot_matrix
      USE input_constants, ONLY: poly_conf
      USE kinds, ONLY: dp
      USE mathconstants, ONLY: fourpi
      INTEGER, PARAMETER                                 :: ng = 6000, nr = 6001
      REAL(dp), PARAMETER :: alpha(2) = [0.4_dp, 1.5_dp], dr = 0.002_dp, &
         gamma_beta(3) = [0.6_dp, 1.7_dp, 0.9_dp]

      INTEGER                                            :: i, j, k, l, order, p, pair, q
      REAL(dp) :: b(2, 3), coefficient(6), egrid, ematrix, errors(5), eta(6), &
         matrix(2, 2, 0:lmat), pmat(2, 2, 0:lmat), r(ng), ref_b, ref_conf(2, 2, 0:lmat), &
         ref_h(2, 2, 0:lmat), ref_nl(2, 2, 0:lmat), ref_ovlp(2, 2, 0:lmat), ref_s(2, 2, 0:lmat), &
         rho(ng), scalar, transformed(2, 2), v(ng), wr(ng)
      TYPE(atom_basis_type)                              :: basis, template
      TYPE(atom_integrals)                               :: integrals
      TYPE(atom_potential_type)                          :: potential

      CALL init_atom_basis_default_pp(template)
      template%nbas = 0
      template%nprim = 0
      template%nbas(0:1) = 2
      template%nprim(0:1) = 2
      template%am(1, :) = alpha(1)
      template%am(2, :) = alpha(2)
      DO k = 1, ng
         r(k) = dr*(REAL(ng - k, dp) + 0.5_dp)
      END DO
      wr = dr
      CALL atom_basis_gridrep(template, basis, r, wr)
      CALL release_atom_basis(template)
      potential%ppot_type = upf_pseudo
      potential%confinement = .TRUE.
      potential%conf_type = poly_conf
      potential%rcon = 2.0_dp
      potential%scon = 2.0_dp
      potential%upf_pot%is_ultrasoft = .TRUE.
      potential%upf_pot%pseudo_type = 'USPP'
      potential%upf_pot%number_of_proj = 3
      potential%upf_pot%mesh_size = nr
      potential%upf_pot%l_max = 1
      ALLOCATE (potential%upf_pot%r(nr), potential%upf_pot%rab(nr), potential%upf_pot%vlocal(nr), &
                potential%upf_pot%beta(nr, 3), potential%upf_pot%lbeta(3), potential%upf_pot%dion(3, 3), &
                potential%upf_pot%qijl(nr, 6, 0:2))
      potential%upf_pot%lbeta(:) = [0, 0, 1]
      potential%upf_pot%dion = 0.0_dp
      potential%upf_pot%dion(1, 1) = 0.3_dp
      potential%upf_pot%dion(1, 2) = -0.1_dp
      potential%upf_pot%dion(2, 1) = -0.1_dp
      potential%upf_pot%dion(2, 2) = 0.6_dp
      potential%upf_pot%dion(3, 3) = -0.2_dp
      potential%upf_pot%qijl = 0.0_dp
      coefficient = [0.2_dp, -0.035_dp, 0.13_dp, 0.0_dp, 0.0_dp, 0.17_dp]
      eta = [0.7_dp, 1.1_dp, 1.6_dp, 1.0_dp, 1.0_dp, 0.9_dp]
      DO k = 1, nr
         scalar = (k - 1)*dr
         potential%upf_pot%r(k) = scalar
         potential%upf_pot%rab(k) = dr
         potential%upf_pot%vlocal(k) = -0.7_dp*EXP(-0.5_dp*scalar**2)
         DO p = 1, 3
            potential%upf_pot%beta(k, p) = scalar**potential%upf_pot%lbeta(p)*EXP(-gamma_beta(p)*scalar**2)
         END DO
         DO pair = 1, 6
            potential%upf_pot%qijl(k, pair, 0) = coefficient(pair)*scalar**2*EXP(-eta(pair)*scalar**2)
         END DO
      END DO
      CALL atom_int_setup(integrals, basis, potential=potential)
      CALL atom_ppint_setup(integrals, basis, potential)
      errors = 0.0_dp
      b = 0.0_dp
      ref_s = 0.0_dp
      ref_h = 0.0_dp
      ref_ovlp = 0.0_dp
      ref_conf = 0.0_dp
      ref_nl = 0.0_dp
      DO p = 1, 3
         l = potential%upf_pot%lbeta(p)
         DO i = 1, 2
            ref_b = uspp_radial_integral(l, alpha(i) + gamma_beta(p))
            b(i, p) = ref_b
            errors(1) = MAX(errors(1), ABS(basis%uspp%beta(i, p) - ref_b))
         END DO
      END DO
      DO l = 0, 1
         DO i = 1, 2
            DO j = 1, 2
               ref_ovlp(i, j, l) = uspp_radial_integral(l, alpha(i) + alpha(j))
               ref_s(i, j, l) = ref_ovlp(i, j, l)
               ref_h(i, j, l) = -0.7_dp*uspp_radial_integral(l, alpha(i) + alpha(j) + 0.5_dp)
               ref_conf(i, j, l) = 0.25_dp*uspp_radial_integral(l + 1, alpha(i) + alpha(j))
               DO p = 1, 3
                  IF (potential%upf_pot%lbeta(p) /= l) CYCLE
                  DO q = 1, 3
                     IF (potential%upf_pot%lbeta(q) /= l) CYCLE
                     pair = MAX(p, q)*(MAX(p, q) - 1)/2 + MIN(p, q)
                     ref_s(i, j, l) = ref_s(i, j, l) + b(i, p)*b(j, q)*coefficient(pair)*uspp_radial_integral(0, eta(pair))
                     ref_h(i, j, l) = ref_h(i, j, l) - &
                                      0.7_dp*b(i, p)*b(j, q)*coefficient(pair)*uspp_radial_integral(0, eta(pair) + 0.5_dp)
                     ref_conf(i, j, l) = ref_conf(i, j, l) + &
                                         0.25_dp*b(i, p)*b(j, q)*coefficient(pair)*uspp_radial_integral(1, eta(pair))
                     ref_nl(i, j, l) = ref_nl(i, j, l) + b(i, p)*potential%upf_pot%dion(p, q)*b(j, q)
                  END DO
               END DO
            END DO
         END DO
         transformed = MATMUL(TRANSPOSE(integrals%utrans(:, :, l)), &
                              MATMUL(integrals%metric(:, :, l), integrals%utrans(:, :, l)))
         transformed(1, 1) = transformed(1, 1) - 1.0_dp
         transformed(2, 2) = transformed(2, 2) - 1.0_dp
         errors(2) = MAX(errors(2), MAXVAL(ABS(transformed)))
      END DO
      errors(1) = MAX(errors(1), MAXVAL(ABS(integrals%ovlp - ref_ovlp)))
      errors(2) = MAX(errors(2), MAXVAL(ABS(integrals%metric - ref_s)))
      errors(3) = MAXVAL(ABS(integrals%core - ref_h))
      errors(3) = MAX(errors(3), MAXVAL(ABS(integrals%hnl - ref_nl)), MAXVAL(ABS(integrals%conf - ref_conf)))
      pmat = 0.0_dp
      pmat(1, 1, 0:1) = 0.8_dp
      pmat(2, 2, 0:1) = 0.45_dp
      pmat(1, 2, 0:1) = -0.1_dp
      pmat(2, 1, 0:1) = -0.1_dp
      CALL atom_density(rho, pmat, basis, 1, typ='RHO')
      errors(4) = ABS(fourpi*SUM(rho*basis%grid%wr) - SUM(pmat*integrals%metric))
      v = EXP(-0.23_dp*r**2)*(1.0_dp + 0.2_dp*r**2)
      DO order = 0, 1
         matrix = 0.0_dp
         CALL numpot_matrix(matrix, v, basis, order)
         IF (order == 0) THEN
            CALL atom_density(rho, pmat, basis, 1, typ='RHO')
         ELSE
            CALL atom_density(rho, pmat, basis, 1, typ='DER')
         END IF
         egrid = fourpi*SUM(rho*v*basis%grid%wr)
         ematrix = SUM(pmat*matrix)
         errors(5) = MAX(errors(5), ABS(egrid - ematrix))
      END DO
      ! Reusing the basis must not retain an obsolete or twice-added augmentation.
      CALL atom_ppint_release(integrals)
      CALL atom_ppint_setup(integrals, basis, potential)
      errors(3) = MAX(errors(3), MAXVAL(ABS(integrals%conf - ref_conf)))
      CALL atom_ppint_release(integrals)
      CALL atom_int_release(integrals)
      CALL atom_int_setup(integrals, basis, potential=potential)
      CALL atom_ppint_setup(integrals, basis, potential)
      errors(2) = MAX(errors(2), MAXVAL(ABS(integrals%metric - ref_s)))
      errors(3) = MAX(errors(3), MAXVAL(ABS(integrals%conf - ref_conf)))
      CALL atom_ppint_release(integrals)
      CALL atom_int_release(integrals)
      WRITE (*, '(A,5ES20.11)') 'Atomic USPP beta/metric/potential/charge/adjoint errors: ', errors
      IF (MAXVAL(errors) > 1.e-8_dp) ERROR STOP 'Atomic USPP analytic reference failed'
      CALL release_atom_basis(basis)
      CALL atom_release_upf(potential%upf_pot)
   END SUBROUTINE check_uspp_atomic

! **************************************************************************************************
!> \brief Analytic Gaussian radial integral including the volume element.
!> \param l Angular momentum of the two factors.
!> \param exponent Sum of the Gaussian exponents.
!> \return Integral of r**(2*l+2)*exp(-exponent*r**2) from zero to infinity.
! **************************************************************************************************
   REAL(dp) FUNCTION uspp_radial_integral(l, exponent) RESULT(value)
      INTEGER, INTENT(IN)                                :: l
      REAL(dp), INTENT(IN)                               :: exponent

      value = 0.5_dp*GAMMA(REAL(l, dp) + 1.5_dp)/exponent**(REAL(l, dp) + 1.5_dp)
   END FUNCTION uspp_radial_integral

! **************************************************************************************************
!> \brief Preserve UPF projector, D and augmentation conventions, including null channels.
! **************************************************************************************************
   SUBROUTINE check_uspp_import()
      CHARACTER(LEN=*), PARAMETER :: fixture = "uspp-augmentation-test.upf"

      INTEGER                                            :: unit
      TYPE(atom_upfpot_type)                             :: pot
      TYPE(cp_logger_type), POINTER                      :: logger

      logger => cp_get_default_logger()
      IF (logger%para_env%is_source()) THEN
         OPEN (NEWUNIT=unit, FILE=fixture, STATUS='REPLACE', ACTION='WRITE')
         WRITE (unit, '(A)') &
            '<UPF version="2.0.1">', &
            '<PP_HEADER element="Pt" pseudo_type="USPP" is_ultrasoft="true"', &
            ' z_valence="16" mesh_size="3" number_of_proj="2" l_max="1" />', &
            '<PP_MESH mesh="3">', &
            '<PP_R>', &
            '0.5 1.0 1.5', &
            '</PP_R>', &
            '<PP_RAB>', &
            '0.5 0.5 0.5', &
            '</PP_RAB>', &
            '</PP_MESH>', &
            '<PP_LOCAL type="real" size="3">', &
            '6 4 2', &
            '</PP_LOCAL>', &
            '<PP_NONLOCAL>', &
            '<PP_BETA.1 type="real" size="3" index="1" angular_momentum="0">', &
            '1 4 9', &
            '</PP_BETA.1>', &
            '<PP_BETA.2 type="real" size="3" index="2" angular_momentum="1">', &
            '2 -1 3', &
            '</PP_BETA.2>', &
            '<PP_DIJ>', &
            '2 0 0 8', &
            '</PP_DIJ>', &
            '<PP_AUGMENTATION q_with_l="T" nqf="0" nqlc="3">', &
            '<PP_Q>', &
            '2 0 0 3', &
            '</PP_Q>', &
            '<PP_QIJL.2.2.0 first_index="2" second_index="2" composite_index="3" angular_momentum="0">', &
            '-3 4 5', &
            '</PP_QIJL.2.2.0>', &
            '<PP_QIJL.1.1.0 type="real" size="3" first_index="1" second_index="1"', &
            ' composite_index="1" angular_momentum="0" is_null="false">', &
            '1 2 3', &
            '</PP_QIJL.1.1.0>', &
            '<PP_QIJL.1.2.1 first_index="1" second_index="2" composite_index="2" angular_momentum="1">', &
            '4 -5 6', &
            '</PP_QIJL.1.2.1>', &
            '<PP_QIJL.2.2.2 first_index="2" second_index="2" composite_index="3"', &
            ' angular_momentum="2" is_null="T" size="0"/>', &
            '</PP_AUGMENTATION>', &
            '</PP_NONLOCAL>', &
            '</UPF>'
         CLOSE (unit)
      END IF
      CALL logger%para_env%sync()
      CALL atom_read_upf(pot, fixture, read_augmentation=.TRUE.)
      IF (ANY(pot%econf(0:3) /= [1, 6, 9, 0])) ERROR STOP 'Wrong 62-electron core subtraction'
      IF (ANY(SHAPE(pot%qijl) /= [3, 3, 3]) .OR. LBOUND(pot%qijl, 3) /= 0) ERROR STOP 'Wrong augmentation shape'
      IF (MAXVAL(ABS(pot%qij - RESHAPE([2._dp, 0._dp, 0._dp, 3._dp], [2, 2]))) > 1.e-14_dp) THEN
         ERROR STOP 'PP_Q normalization changed'
      END IF
      IF (MAXVAL(ABS(pot%qijl(:, 1, 0) - [1._dp, 2._dp, 3._dp])) > 1.e-14_dp) ERROR STOP 'Wrong s-s augmentation'
      IF (MAXVAL(ABS(pot%qijl(:, 2, 1) - [4._dp, -5._dp, 6._dp])) > 1.e-14_dp) ERROR STOP 'Wrong s-p augmentation'
      IF (MAXVAL(ABS(pot%qijl(:, 3, 0) - [-3._dp, 4._dp, 5._dp])) > 1.e-14_dp) ERROR STOP 'Wrong p-p augmentation'
      IF (ANY(pot%qijl(:, :, 2) /= 0._dp)) ERROR STOP 'Null augmentation was not zero initialized'
      IF (MAXVAL(ABS(pot%beta(:, 2) - [4._dp, -1._dp, 2._dp])) > 1.e-14_dp) ERROR STOP 'USPP beta normalization changed'
      IF (ABS(pot%dion(2, 2) - 4._dp) > 1.e-14_dp) ERROR STOP 'USPP D is not in Hartree'
      CALL atom_release_upf(pot)
      IF (ALLOCATED(pot%qij) .OR. ALLOCATED(pot%qijl)) ERROR STOP 'Augmentation arrays were not released'
      CALL logger%para_env%sync()
      IF (logger%para_env%is_source()) THEN
         OPEN (NEWUNIT=unit, FILE=fixture, STATUS='OLD')
         CLOSE (unit, STATUS='DELETE')
      END IF

   END SUBROUTINE check_uspp_import

! **************************************************************************************************
!> \brief Report native augmentation invariants for an optional external original UPF file.
!> \param pot Original ultrasoft data.
!> \param filename Dataset path used to identify the optional audit.
! **************************************************************************************************
   SUBROUTINE check_uspp_dataset(pot, filename)
      TYPE(atom_upfpot_type), INTENT(IN)                 :: pot
      CHARACTER(LEN=*), INTENT(IN)                       :: filename

      INTEGER                                            :: channels, j, l, n, ns, pair
      REAL(dp)                                           :: error
      REAL(dp), ALLOCATABLE                              :: weights(:)
      TYPE(cp_logger_type), POINTER                      :: logger

      logger => cp_get_default_logger()
      IF (.NOT. ALL(ieee_is_finite(pot%qij)) .OR. .NOT. ALL(ieee_is_finite(pot%qijl))) THEN
         ERROR STOP "Nonfinite USPP augmentation"
      END IF
      n = SIZE(pot%r)
      ALLOCATE (weights(n))
      weights = 0._dp
      ns = n - 1 + MODULO(n, 2)
      weights(1:ns:2) = 2._dp/3._dp
      weights(2:ns - 1:2) = 4._dp/3._dp
      weights(1) = 1._dp/3._dp
      weights(ns) = 1._dp/3._dp
      IF (ns < n) weights(n - 1:n) = weights(n - 1:n) + .5_dp
      weights(:) = weights*pot%rab
      error = 0._dp
      channels = 0
      DO j = 1, pot%number_of_proj
         DO l = 1, j
            pair = j*(j - 1)/2 + l
            channels = channels + MIN(pot%lbeta(j), pot%lbeta(l)) + 1
            IF (pot%lbeta(l) == pot%lbeta(j)) THEN
               error = MAX(error, ABS(SUM(weights*pot%qijl(:, pair, 0)) - pot%qij(l, j)))
            END IF
         END DO
      END DO
      IF (logger%para_env%is_source()) THEN
         WRITE (*, '(A,1X,A,1X,I0,4(1X,ES25.16))') 'USPP augmentation:', TRIM(filename), channels, &
            SUM(pot%qij), SUM(pot%qijl), SUM(ABS(pot%qijl)), error
      END IF
      DEALLOCATE (weights)
   END SUBROUTINE check_uspp_dataset
END PROGRAM atom_upf_unittest
