!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM atom_upf_unittest
   USE atom_kind_orbitals,              ONLY: sgp_potential_conversion
   USE atom_operators,                  ONLY: atom_ppint_release,&
                                              atom_ppint_setup
   USE atom_sgp,                        ONLY: atom_sgp_potential_type,&
                                              atom_sgp_release,&
                                              sgp_construction
   USE atom_types,                      ONLY: atom_basis_type,&
                                              atom_integrals,&
                                              atom_potential_type,&
                                              init_atom_basis_default_pp,&
                                              release_atom_basis,&
                                              sgp_pseudo
   USE atom_upf,                        ONLY: atom_read_upf,&
                                              atom_release_upf,&
                                              atom_upfpot_type
   USE atom_utils,                      ONLY: atom_core_density
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE external_potential_types,        ONLY: allocate_potential,&
                                              deallocate_potential,&
                                              init_potential,&
                                              sgp_potential_type
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE libcp2k,                         ONLY: cp2k_finalize,&
                                              cp2k_init
   USE mathconstants,                   ONLY: fourpi
   USE orbital_pointers,                ONLY: co,&
                                              init_orbital_pointers,&
                                              nco,&
                                              nso
   USE orbital_transformation_matrices, ONLY: init_spherical_harmonics

   IMPLICIT NONE

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
      CALL atom_read_upf(pot, TRIM(filename))
      IF (ABS(SUM(pot%econf) - pot%zion) > 1.e-10_dp) ERROR STOP "UPF valence count mismatch"
      IF (.NOT. ALL(ieee_is_finite(pot%vlocal))) ERROR STOP "Non-finite UPF local potential"
      IF (ANY(pot%r(2:) <= pot%r(:SIZE(pot%r) - 1))) ERROR STOP "UPF mesh is not increasing"
      IF (ANY(pot%rab <= 0.0_dp)) ERROR STOP "Non-positive UPF mesh weights"
      IF (ALLOCATED(pot%dion)) THEN
         IF (.NOT. ALL(ieee_is_finite(pot%dion))) ERROR STOP "Non-finite UPF projector coefficients"
         IF (.NOT. ALL(ieee_is_finite(pot%beta))) ERROR STOP "Non-finite UPF projectors"
         IF (MAXVAL(ABS(pot%dion - TRANSPOSE(pot%dion))) > 1.e-10_dp) ERROR STOP "Non-Hermitian UPF D matrix"
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
!> \brief Read an independent miniature UPF fixture and check its numerical convention.
!>        The values are deliberately simple and are not a physical pseudopotential.
!>        Checking beta separately matters: an inverse rescaling of D leaves the NC
!>        operator unchanged but would require rescaling USPP augmentation data too.
! **************************************************************************************************
   SUBROUTINE check_upf_units()
      CHARACTER(LEN=*), PARAMETER                        :: fixture = "upf-units-test.upf"

      INTEGER                                            :: unit
      TYPE(atom_upfpot_type)                             :: potential
      TYPE(cp_logger_type), POINTER                      :: logger

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
      CALL atom_read_upf(potential, fixture)
      IF (MAXVAL(ABS(potential%vlocal - [3.0_dp, 2.0_dp, 1.0_dp])) > 1.e-14_dp) &
         ERROR STOP "UPF local potential is not in Hartree"
      IF (MAXVAL(ABS(potential%beta(:, 1) - [2.0_dp, 4.0_dp, 6.0_dp])) > 1.e-14_dp) &
         ERROR STOP "UPF changed the first tabulated projector amplitude"
      IF (MAXVAL(ABS(potential%beta(:, 2) - [4.0_dp, -1.0_dp, 2.0_dp])) > 1.e-14_dp) &
         ERROR STOP "UPF changed the second tabulated projector amplitude"
      IF (MAXVAL(ABS(potential%dion - RESHAPE([1.0_dp, 0.2_dp, 0.2_dp, 4.0_dp], [2, 2]))) > 1.e-14_dp) &
         ERROR STOP "UPF D matrix is not in the tabulated-projector convention"
      CALL atom_release_upf(potential)
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

END PROGRAM atom_upf_unittest
