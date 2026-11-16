!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

!> \author Jan Wilhelm, with help from ChatGPT

MODULE nlo
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_get_info,&
                                              dbcsr_p_type
   USE cp_dbcsr_operations,             ONLY: dbcsr_deallocate_matrix_set
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
   USE input_constants,                 ONLY: do_method_gapw,&
                                              do_method_gpw
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE kpoint_k_r_trafo_simple,         ONLY: replicate_rs_matrices
   USE kpoint_methods,                  ONLY: kpoint_init_cell_index
   USE kpoint_types,                    ONLY: kpoint_create,&
                                              kpoint_type
   USE mathlib,                         ONLY: vector_product
   USE nlo_integration,                 ONLY: nlo_bz_integration
   USE nlo_read_input,                  ONLY: create_nonlinear_optics_section,&
                                              nlo_frequency_case,&
                                              read_nlo_input
   USE nlo_types,                       ONLY: nlo_type,&
                                              release_nlo_environment
   USE physcon,                         ONLY: a_bohr,&
                                              c_light,&
                                              e_charge,&
                                              evolt,&
                                              femtoseconds,&
                                              joule,&
                                              mu_perm
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_moments,                      ONLY: build_local_moment_matrix_rs_img
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE scf_control_types,               ONLY: scf_control_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'nlo'

   PUBLIC :: create_nonlinear_optics_section, run_nlo

CONTAINS

! **************************************************************************************************
!> \brief Calculate frequency-dependent χ^(2) to χ^(5) after the ground-state SCF.
!> \param qs_env Converged k-point SCF environment
!> \note The independent-particle approximation uses the converged SCF Hamiltonian.
!>       Excitonic effects and time-dependent screening (TDDFT, GW-BSE) are absent.
!>
!>       Method: length-gauge density-matrix response in the spirit of
!>       Aversa and Sipe, Phys. Rev. B 52, 14636 (1995),
!>       DOI: 10.1103/PhysRevB.52.14636; analytic covariant k derivatives follow
!>       Ventura et al., Phys. Rev. B 96, 035431 (2017),
!>       DOI: 10.1103/PhysRevB.96.035431. The finite AO projection, analytic
!>       position curvature, and numerical implementation are defined below
!>       and in nlo_current_tensor_derivatives.
!>
!>       1. Read the optical k mesh, incident frequencies, tensor elements,
!>          and either a common time T₁ = T₂ or separate relaxation times.
!>
!>       2. Build lattice-image matrices for the AO Hamiltonian H, overlap S,
!>          dipole d.
!>          d_μν(R) = ⟨ϕ_(μ0)|r|ϕ_(νR)⟩ is a position matrix; the electron
!>          charge enters the later SI conversion. At each optical k point:
!>
!>          X(k) = Σ_R exp(2π i k_frac · R) X(R), with X = H, S, d.
!>
!>          Differentiating the phase gives every required k derivative.
!>          The AO overlap is removed, then H and position matrices are
!>          expressed in the eigenvector basis of H at the current k point.
!>
!>       3. At each optical k point, diagonalize the Hamiltonian once.
!>          The occupied bands define the equilibrium density matrix.
!>          A light field changes that density matrix: first order uses the
!>          equilibrium matrix, second order uses first-order responses,
!>          and so on through fifth order. Each new response is obtained
!>          from the previous one by an analytic k derivative and a
!>          band-energy denominator ε_m - ε_n - Ω - iΓ_mn. Here ε_m and ε_n
!>          are band energies, Ω is the sum of the incident photon energies,
!>          and Γ is the band-pair relaxation rate. For T₁ = T₂ it is the
!>          constant 1/T₁. For separate times, with γᵢ = 1/Tᵢ and
!>          Δε = ε_m - ε_n, it is
!>
!>          Γ_mn = γ₁ + (γ₂ - γ₁)[1 - exp(-(Δε/w)²)].
!>
!>          The positive width w sets the smooth transition at nearly
!>          degenerate bands. Its basis-covariant operator form is
!>          Γ_H(X) = γ₂X + (γ₁ - γ₂)exp[-(ad_H/w)²]X,
!>          where ad_H(X) = [H,X]. Its analytic k derivatives enter the
!>          response. All possible orders of
!>          the incident fields are averaged. No numerical derivatives of
!>          separately diagonalized bands or their phases are used.
!>
!>          The position operator in the band representation is
!>
!>          r_a = i ∂/∂k_a + A_a .
!>
!>          The subscript a denotes x, y or z. A_a is the matrix obtained
!>          from the AO dipole matrix after orthonormalizing and rotating
!>          into the central band basis. It is the finite-basis part of
!>          the position operator. When a field changes the current,
!>          an additional term from this matrix is required:
!>
!>          F_ia = ∂A_a/∂k_i - ∂A_i/∂k_a - i [A_i, A_a] .
!>
!>          Here i is the output-current direction. F_ia is the position
!>          curvature; its contribution to the current is always included.
!>          The other current contribution is the band-velocity matrix
!>          multiplied by the density response.
!>          Summing over optical k points and spin gives the
!>          cell-integrated current response I^(n) in atomic units,
!>          stored locally on each MPI rank.
!>
!>       4. Sum I^(n) across MPI ranks, convert it to SI units and
!>          divide by the simulation-cell volume V. The result is the
!>          nonlinear conductivity σ^(n) for each field order. For an
!>          output frequency Ω, the corresponding polarization obeys
!>          J(Ω) = -i Ω P(Ω). Hence:
!>
!>          χ²_iab     = σ²_iab     / [-i (ω₁+ω₂) ε₀]
!>          χ³_iabc    = σ³_iabc    / [-i (ω₁+ω₂+ω₃) ε₀]
!>          χ⁴_iabcd   = σ⁴_iabcd   / [-i (ω₁+ω₂+ω₃+ω₄) ε₀]
!>          χ⁵_iabcde  = σ⁵_iabcde  / [-i (ω₁+ω₂+ω₃+ω₄+ω₅) ε₀] .
!>
!>          The SI conversion is in nlo_write_susceptibility. For a
!>          two-dimensional slab of in-plane area A, the sheet response
!>          is χ_sheet = (V/A) χ_volume. Ω must be nonzero.
!>
!>       5. Release the NLO matrices.
! **************************************************************************************************
   SUBROUTINE run_nlo(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER                        :: routineN = 'run_nlo'

      INTEGER                                            :: handle
      LOGICAL                                            :: nlo_present
      TYPE(nlo_type)                                     :: nlo_env
      TYPE(section_vals_type), POINTER                   :: input, nlo_section

      CALL get_qs_env(qs_env, input=input)
      nlo_section => section_vals_get_subs_vals(input, 'PROPERTIES%NONLINEAR_OPTICS')
      CALL section_vals_val_get(nlo_section, '_SECTION_PARAMETERS_', l_val=nlo_present)
      IF (.NOT. nlo_present) RETURN

      CALL timeset(routineN, handle)

      ! 1. Read and validate the nonlinear optical input.
      CALL read_nlo_input(qs_env, nlo_env)

      ! 2. Prepare the AO matrices from the converged SCF calculation.
      CALL nlo_prepare_environment(qs_env, nlo_env)

      ! 3. Calculate the density response and integrate the current over k.
      ! I^(n) = -sum_(spin,k) w_k Tr[v(E) rho^(n)] in atomic units.
      CALL nlo_bz_integration(nlo_env)

      ! 4. Convert the current to susceptibility and print the requested tensors.
      ! chi^(n) = I^(n) / [V epsilon_0 (-i Omega)] after SI conversion.
      CALL nlo_write_susceptibility(nlo_env)

      ! 5. Release the NLO matrices.
      CALL release_nlo_environment(nlo_env)

      CALL timestop(handle)

   END SUBROUTINE run_nlo

! **************************************************************************************************
!> \brief Build the real-space Hamiltonian, overlap and position matrices from the SCF result.
!> \param qs_env Converged k-point SCF environment
!> \param nlo_env NLO settings and data for the calculation
! **************************************************************************************************
   SUBROUTINE nlo_prepare_environment(qs_env, nlo_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(nlo_type), INTENT(INOUT)                      :: nlo_env

      CHARACTER(len=*), PARAMETER :: routineN = 'nlo_prepare_environment'

      INTEGER                                            :: handle, ik, iw
      LOGICAL                                            :: scf_has_kpoints
      TYPE(dbcsr_p_type), POINTER                        :: dr(:, :), hr(:, :), sr(:, :)
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: scf
      TYPE(neighbor_list_set_p_type), POINTER            :: all_pairs(:)
      TYPE(scf_control_type), POINTER                    :: scf_control

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, do_kpoints=scf_has_kpoints, dft_control=dft_control, &
                      scf_control=scf_control)
      IF (scf_control%smear%do_smear) CPABORT('NONLINEAR_OPTICS requires SCF without smearing.')
      IF (.NOT. scf_has_kpoints) THEN
         CPABORT('NONLINEAR_OPTICS requires k-point SCF; Gamma-only SCF is unsupported.')
      END IF
      SELECT CASE (dft_control%qs_control%method_id)
      CASE (do_method_gpw, do_method_gapw)
      CASE DEFAULT
         CPABORT('NONLINEAR_OPTICS currently supports GPW/GAPW only.')
      END SELECT
      nlo_env%n_kpoints = 1
      DO ik = 1, 3
         IF (nlo_env%n_kpoints > HUGE(nlo_env%n_kpoints)/nlo_env%k_mesh(ik)) THEN
            CPABORT('NONLINEAR_OPTICS K_MESH exceeds integer indexing.')
         END IF
         nlo_env%n_kpoints = nlo_env%n_kpoints*nlo_env%k_mesh(ik)
      END DO
      NULLIFY (nlo_env%images, dr)
      CALL get_qs_env(qs_env, cell=nlo_env%cell, para_env=nlo_env%para_env, matrix_ks_kp=hr, &
                      matrix_s_kp=sr, sab_all=all_pairs, kpoints=scf, mos=nlo_env%mos)
      IF (ANY(nlo_env%cell%perd == 0 .AND. nlo_env%k_mesh /= 1)) THEN
         CPABORT('NONLINEAR_OPTICS K_MESH must be one in each nonperiodic direction.')
      END IF
      IF (nlo_env%cell%deth <= 0.0_dp) CPABORT('NONLINEAR_OPTICS requires positive cell volume.')
      CALL dbcsr_get_info(sr(1, 1)%matrix, nfullrows_total=nlo_env%n_ao)
      nlo_env%n_spin = SIZE(hr, 1)
      nlo_env%n_bands = nlo_env%n_ao
      CALL kpoint_create(nlo_env%images)
      CALL kpoint_init_cell_index(nlo_env%images, all_pairs, nlo_env%para_env, nlo_env%n_images)
      ALLOCATE (dr(3, nlo_env%n_images))
      CALL build_local_moment_matrix_rs_img(qs_env, dr)
      ALLOCATE (nlo_env%h_ao_R(nlo_env%n_spin, nlo_env%n_ao, nlo_env%n_ao, nlo_env%n_images), &
                SOURCE=0.0_dp)
      ALLOCATE (nlo_env%s_ao_R(1, nlo_env%n_ao, nlo_env%n_ao, nlo_env%n_images), SOURCE=0.0_dp)
      ALLOCATE (nlo_env%d_ao_R(3, nlo_env%n_ao, nlo_env%n_ao, nlo_env%n_images), SOURCE=0.0_dp)
      ! H(R), S(R), d(R) supply X(k) = sum_R exp(2 pi i k_frac.R) X(R).
      CALL replicate_rs_matrices(hr, scf, nlo_env%h_ao_R, nlo_env%images%cell_to_index)
      CALL replicate_rs_matrices(sr, scf, nlo_env%s_ao_R, nlo_env%images%cell_to_index)
      CALL replicate_rs_matrices(dr, nlo_env%images, nlo_env%d_ao_R, nlo_env%images%cell_to_index)
      CALL dbcsr_deallocate_matrix_set(dr)
      iw = cp_logger_get_default_io_unit()
      IF (iw > 0) THEN
         WRITE (iw, '()')
         WRITE (iw, '()')
         WRITE (iw, '(T2,A,T69,I12)') 'NONLINEAR_OPTICS| Selected tensor elements:', &
            COUNT(nlo_env%elements)
         WRITE (iw, '(T2,A,T69,I12)') 'NONLINEAR_OPTICS| Bands in response:', nlo_env%n_bands
         IF (nlo_env%t1_au == nlo_env%t2_au) THEN
            WRITE (iw, '(T2,A,T70,F8.1,1X,A)') 'NONLINEAR_OPTICS| T1 = T2:', &
               nlo_env%t1_au*femtoseconds, 'fs'
            WRITE (iw, '(T2,A,T73,F8.1,1X,A)') 'NONLINEAR_OPTICS| ℏ/T1 = ℏ/T2:', &
               1000.0_dp*evolt/nlo_env%t1_au, 'meV'
         ELSE
            WRITE (iw, '(T2,A,T70,F8.1,1X,A)') 'NONLINEAR_OPTICS| T1:', &
               nlo_env%t1_au*femtoseconds, 'fs'
            WRITE (iw, '(T2,A,T70,F8.1,1X,A)') 'NONLINEAR_OPTICS| T2:', &
               nlo_env%t2_au*femtoseconds, 'fs'
            WRITE (iw, '(T2,A,T71,F8.1,1X,A)') 'NONLINEAR_OPTICS| ℏ/T1:', &
               1000.0_dp*evolt/nlo_env%t1_au, 'meV'
            WRITE (iw, '(T2,A,T71,F8.1,1X,A)') 'NONLINEAR_OPTICS| ℏ/T2:', &
               1000.0_dp*evolt/nlo_env%t2_au, 'meV'
            WRITE (iw, '(T2,A,T69,F8.1,1X,A)') 'NONLINEAR_OPTICS| Band transition width:', &
               1000.0_dp*nlo_env%transition_width_au*evolt, 'meV'
         END IF
      END IF

      CALL timestop(handle)

   END SUBROUTINE nlo_prepare_environment

! **************************************************************************************************
!> \brief Sum the rank-local currents and convert them to volume and sheet susceptibilities.
!> \param nlo_env NLO settings and data for the calculation
! **************************************************************************************************
   SUBROUTINE nlo_write_susceptibility(nlo_env)
      TYPE(nlo_type), INTENT(INOUT)                      :: nlo_env

      CHARACTER(len=*), PARAMETER :: routineN = 'nlo_write_susceptibility'
      CHARACTER(len=3), PARAMETER :: &
         omega_subscript(5) = [CHARACTER(len=3) :: '₁', '₂', '₃', '₄', '₅']

      CHARACTER(len=80)                                  :: heading_2d, heading_3d
      COMPLEX(dp), ALLOCATABLE                           :: chi_flat(:, :), current_flat_sum(:, :)
      INTEGER                                            :: a, c, code, handle, i_case, iw
      LOGICAL                                            :: zero_output
      REAL(KIND=dp)                                      :: area, energy_ev(5), factor, omega(5), &
                                                            sheet_height

      CALL timeset(routineN, handle)

      ALLOCATE (current_flat_sum(3, 3**nlo_env%order), chi_flat(3, 3**nlo_env%order))
      SELECT CASE (nlo_env%order)
      CASE (2)
         heading_3d = 'NONLINEAR_OPTICS| χ² 3D (m/V)'
         heading_2d = 'NONLINEAR_OPTICS| χ² 2D (m²/V)'
      CASE (3)
         heading_3d = 'NONLINEAR_OPTICS| χ³ 3D (m²/V²)'
         heading_2d = 'NONLINEAR_OPTICS| χ³ 2D (m³/V²)'
      CASE (4)
         heading_3d = 'NONLINEAR_OPTICS| χ⁴ 3D (m³/V³)'
         heading_2d = 'NONLINEAR_OPTICS| χ⁴ 2D (m⁴/V³)'
      CASE (5)
         heading_3d = 'NONLINEAR_OPTICS| χ⁵ 3D (m⁴/V⁴)'
         heading_2d = 'NONLINEAR_OPTICS| χ⁵ 2D (m⁵/V⁴)'
      END SELECT
      DO i_case = 1, nlo_env%n_frequency_combinations
         CALL nlo_frequency_case(nlo_env%frequency, nlo_env%order, i_case, &
                                 energy_ev, zero_output)
         IF (zero_output) CPABORT('Invalid NONLINEAR_OPTICS frequency case.')
         omega = 0.0_dp
         omega(1:nlo_env%order) = energy_ev(1:nlo_env%order)/evolt
         current_flat_sum(:, :) = nlo_env%current_spectrum(:, :, i_case)
         CALL nlo_env%para_env%sum(current_flat_sum)
         factor = e_charge**(nlo_env%order + 1)*a_bohr**(nlo_env%order - 2)* &
                  mu_perm*c_light**2/(joule**nlo_env%order*nlo_env%cell%deth)
         chi_flat(:, :) = current_flat_sum*factor/CMPLX(0.0_dp, -SUM(omega), dp)
         area = SQRT(SUM(vector_product(nlo_env%cell%hmat(:, 1), &
                                        nlo_env%cell%hmat(:, 2))**2))
         IF (area <= 0.0_dp) CPABORT('NONLINEAR_OPTICS requires nonzero in-plane area.')
         sheet_height = nlo_env%cell%deth/area*a_bohr
         iw = cp_logger_get_default_io_unit()
         IF (iw <= 0) CYCLE
         WRITE (iw, '()')
         WRITE (iw, '(T2,A,I0,A,I0)') 'NONLINEAR_OPTICS| Frequency case ', &
            i_case, ' / ', nlo_env%n_frequency_combinations
         ! Adjust tab columns for the multibyte ω and subscript characters.
         DO a = 1, nlo_env%order
            WRITE (iw, '(T2,A,A,A,T72,F9.3,1X,A)') 'NONLINEAR_OPTICS| ω', &
               omega_subscript(a), ' =', energy_ev(a), 'eV'
         END DO
         WRITE (iw, '(T2,A,T70,F9.3,1X,A)') 'NONLINEAR_OPTICS| ω_out =', &
            SUM(energy_ev(1:nlo_env%order)), 'eV'
         WRITE (iw, '(T2,A,T50,A,T72,A)') TRIM(heading_3d), 'Real', 'Imag'
         DO code = 1, 3**nlo_env%order
            DO c = 1, 3
               IF (.NOT. nlo_env%elements(c, code)) CYCLE
               WRITE (iw, '(T2,A,A6,T37,2ES22.12)') 'NONLINEAR_OPTICS| 3D ', &
                  nlo_tensor_indices(nlo_env%order, c, code), &
                  REAL(chi_flat(c, code), KIND=dp), AIMAG(chi_flat(c, code))
            END DO
         END DO
         IF (ALL(nlo_env%cell%perd(1:2) /= 0) .AND. nlo_env%cell%perd(3) == 0) THEN
            WRITE (iw, '(T2,A,T62,ES15.6,1X,A)') &
               'NONLINEAR_OPTICS| Effective sheet height:', sheet_height, 'm'
            WRITE (iw, '(T2,A,T50,A,T72,A)') TRIM(heading_2d), 'Real', 'Imag'
            DO code = 1, 3**nlo_env%order
               DO c = 1, 3
                  IF (.NOT. nlo_env%elements(c, code)) CYCLE
                  WRITE (iw, '(T2,A,A6,T37,2ES22.12)') 'NONLINEAR_OPTICS| 2D ', &
                     nlo_tensor_indices(nlo_env%order, c, code), &
                     REAL(chi_flat(c, code)*sheet_height, KIND=dp), &
                     AIMAG(chi_flat(c, code)*sheet_height)
               END DO
            END DO
         END IF
      END DO

      CALL timestop(handle)

   END SUBROUTINE nlo_write_susceptibility

! **************************************************************************************************
!> \brief Decode the output and input Cartesian directions of one tensor element.
!> \param order Rank of the susceptibility tensor
!> \param output_direction Cartesian direction of the induced current
!> \param code Base-three code for the field directions
!> \return Cartesian indices for one tensor component
! **************************************************************************************************
   PURE FUNCTION nlo_tensor_indices(order, output_direction, code) RESULT(indices)
      INTEGER, INTENT(IN)                                :: order, output_direction, code
      CHARACTER(len=6)                                   :: indices

      CHARACTER(len=*), PARAMETER                        :: axes = 'XYZ'

      INTEGER                                            :: field, input_direction, remaining_code

      indices = ' '
      indices(1:1) = axes(output_direction:output_direction)
      remaining_code = code - 1
      DO field = 1, order
         input_direction = MOD(remaining_code, 3) + 1
         indices(field + 1:field + 1) = axes(input_direction:input_direction)
         remaining_code = remaining_code/3
      END DO

   END FUNCTION nlo_tensor_indices

END MODULE nlo
