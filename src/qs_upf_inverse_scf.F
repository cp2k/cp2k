!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Explicit inverse-orbital attachment to native Quickstep density and operator consumers.
!>        Scalar semilocal GPW at Gamma, with original full-wave USPP data.
! **************************************************************************************************
MODULE qs_upf_inverse_scf
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_add,&
                                              dbcsr_copy,&
                                              dbcsr_p_type,&
                                              dbcsr_release,&
                                              dbcsr_set,&
                                              dbcsr_type
   USE cp_dbcsr_operations,             ONLY: dbcsr_deallocate_matrix_set
   USE input_constants,                 ONLY: sic_none,&
                                              uspp_method_direct,&
                                              uspp_method_inverse,&
                                              uspp_method_inverse_1c
   USE kinds,                           ONLY: dp
   USE kpoint_transitional,             ONLY: set_1d_pointer
   USE mathlib,                         ONLY: invmat
   USE message_passing,                 ONLY: mp_para_env_type
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_grid_types,                   ONLY: pw_grid_type
   USE pw_methods,                      ONLY: pw_integrate_function,&
                                              pw_scale,&
                                              pw_transfer
   USE pw_pool_types,                   ONLY: pw_pool_type
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_kinetic,                      ONLY: build_kinetic_matrix
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_did_change,&
                                              qs_ks_env_type,&
                                              set_ks_env
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_set,&
                                              qs_rho_type
   USE qs_upf_inverse,                  ONLY: qs_upf_inverse_atomic_operators
   USE qs_upf_inverse_augmentation,     ONLY: qs_upf_inverse_augmentation_prepare
   USE qs_upf_inverse_grid,             ONLY: qs_upf_inverse_grid_collocate,&
                                              qs_upf_inverse_grid_integrate,&
                                              qs_upf_inverse_grid_metric,&
                                              qs_upf_inverse_grid_prepare,&
                                              qs_upf_inverse_grid_response,&
                                              qs_upf_inverse_matrix_gather,&
                                              qs_upf_inverse_matrix_scatter
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_augmentation_type,&
                                              qs_upf_inverse_grid_type,&
                                              qs_upf_inverse_scf_release,&
                                              qs_upf_inverse_scf_type
   USE qs_uspp_types,                   ONLY: qs_uspp_init,&
                                              qs_uspp_release,&
                                              qs_uspp_type
   USE upf_inverse,                     ONLY: upf_inverse_fields_init
   USE upf_local_integrals,             ONLY: upf_local_fourier,&
                                              upf_local_init
   USE upf_radial_fields,               ONLY: upf_radial_fields_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   PUBLIC :: qs_upf_inverse_scf_prepare, qs_upf_inverse_scf_active
   PUBLIC :: qs_upf_inverse_scf_collocate, qs_upf_inverse_scf_integrate, qs_upf_inverse_scf_core
   PUBLIC :: qs_upf_inverse_scf_initialize, qs_upf_inverse_scf_reset
   PUBLIC :: qs_upf_inverse_scf_core_response
   PUBLIC :: qs_upf_inverse_scf_response
CONTAINS

! **************************************************************************************************
!> \brief Explicitly install one inverse frame, its charge metric, and core operator in native KS.
!>        Never inferred from full-wave availability; caller supplies the fixed image set.
!> \param ks_env Prepared native KS environment.
!> \param grid Native physical-cell PW grid.
!> \param images Fixed physical-cell images for the Gaussian and radial fields.
!> \param info Factorization status; nonzero leaves the reference KS environment unchanged.
!> \param atom_local Optional nonoverlapping one-center inverse with native Gaussian kinetic integrals.
!> \param direct Use original projectors and soft orbitals without reconstruction.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_scf_prepare(ks_env, grid, images, info, atom_local, direct)
      TYPE(qs_ks_env_type), INTENT(INOUT)                :: ks_env
      TYPE(pw_grid_type), POINTER                        :: grid
      INTEGER, INTENT(IN)                                :: images(:, :)
      INTEGER, INTENT(OUT)                               :: info
      LOGICAL, INTENT(IN), OPTIONAL                      :: atom_local, direct

      INTEGER                                            :: ikind
      LOGICAL                                            :: do_kpoints
      TYPE(dbcsr_p_type), POINTER                        :: kinetic(:), metric(:), overlap(:)
      TYPE(dbcsr_type)                                   :: imaginary
      TYPE(dft_control_type), POINTER                    :: control
      TYPE(qs_kind_type), POINTER                        :: kinds(:)
      TYPE(qs_upf_inverse_scf_type), POINTER             :: previous, state
      TYPE(qs_uspp_type), POINTER                        :: uspp

      CALL get_ks_env(ks_env, upf_inverse=previous, uspp=uspp, &
                      dft_control=control, do_kpoints=do_kpoints, qs_kind_set=kinds)
      CPASSERT(.NOT. ASSOCIATED(previous))
      CPASSERT(ASSOCIATED(uspp) .AND. ASSOCIATED(grid))
      IF (.NOT. control%qs_control%gpw .OR. control%qs_control%gapw .OR. control%qs_control%gapw_xc) THEN
         CPABORT("Inverse SCF GAPW reconstruction is not assembled")
      END IF
      IF (control%use_kinetic_energy_density) THEN
         CPABORT("Missing inverse kinetic augmentation for tau XC")
      END IF
      IF (control%drho_by_collocation) THEN
         CPABORT("Inverse SCF needs spectral XC gradients")
      END IF
      IF (control%do_admm .OR. control%do_sccs .OR. control%qs_control%cdft .OR. &
        control%qs_control%rigpw .OR. control%qs_control%lrigpw .OR. control%dft_plus_u .OR. control%sic_method_id /= sic_none) THEN
         CPABORT("Unsupported inverse SCF density consumer")
      END IF
      DO ikind = 1, SIZE(kinds)
         IF (kinds(ikind)%ghost) CYCLE
         CPASSERT(ASSOCIATED(kinds(ikind)%sgp_potential))
         CPASSERT(ALLOCATED(kinds(ikind)%sgp_potential%upf_data))
         IF (.NOT. kinds(ikind)%sgp_potential%upf_data%is_ultrasoft) THEN
            CPABORT("Missing mixed NC inverse nonlocal operator")
         END IF
      END DO
      ALLOCATE (state)
      IF (do_kpoints) THEN
         CPABORT("Inverse USPP k-point consumers are not included in this implementation")
      END IF
      CPASSERT(control%nimages == 1)
      CALL get_ks_env(ks_env, matrix_s=overlap)
      CALL qs_upf_inverse_grid_prepare(ks_env, grid, images, [0.0_dp, 0.0_dp, 0.0_dp], state%frame, info, atom_local, direct)
      IF (info /= 0) THEN
         CALL qs_upf_inverse_scf_release(state)
         RETURN
      END IF
      IF (state%frame%ao%atom_local) THEN
         CALL get_ks_env(ks_env, kinetic=kinetic)
         CALL qs_upf_inverse_matrix_gather(state%frame%ao, kinetic(1)%matrix, state%gaussian_kinetic)
      END IF
      CALL qs_upf_inverse_augmentation_prepare(ks_env, state%frame%ao, state%augmentation)
      ALLOCATE (metric(1))
      ALLOCATE (metric(1)%matrix)
      CALL qs_upf_inverse_grid_metric(state%frame, state%frame%ao%radial_q, overlap(1)%matrix, &
                                      metric(1)%matrix, imaginary)
      CALL check_real(state, metric(1)%matrix, imaginary)
      CALL dbcsr_release(imaginary)
      CALL set_1d_pointer(state%metric, metric)
      ! Retain all inverse-mediated AO fill in native density/KS allocation as well as in M.
      CALL dbcsr_copy(uspp%ao_pattern, metric(1)%matrix)
      CALL dbcsr_set(uspp%ao_pattern, 0.0_dp)
      state%enabled = .TRUE.
      CALL set_ks_env(ks_env, upf_inverse=state)
      CALL qs_upf_inverse_scf_core(ks_env)
      CALL qs_ks_did_change(ks_env, s_mstruct_changed=.TRUE., rho_changed=.TRUE., potential_changed=.TRUE.)
   END SUBROUTINE qs_upf_inverse_scf_prepare

! **************************************************************************************************
!> \brief Test explicit attachment; full-wave data alone never activates this representation.
!> \param ks_env ...
!> \return ...
! **************************************************************************************************
   LOGICAL FUNCTION qs_upf_inverse_scf_active(ks_env) RESULT(active)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env

      TYPE(qs_upf_inverse_scf_type), POINTER             :: state

      CALL get_ks_env(ks_env, upf_inverse=state)
      active = ASSOCIATED(state)
      IF (active) active = state%enabled
   END FUNCTION qs_upf_inverse_scf_active

! **************************************************************************************************
!> \brief Build the original short-range local, Dion and soft kinetic operator in the inverse frame.
!>        The long-range ionic/core Coulomb term remains in native Poisson energy accounting.
!> \param ks_env ...
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_scf_core(ks_env)
      TYPE(qs_ks_env_type), INTENT(INOUT)                :: ks_env

      COMPLEX(dp), ALLOCATABLE                           :: atomic_m(:, :), atomic_t(:, :)
      INTEGER                                            :: axis
      TYPE(dbcsr_p_type), POINTER                        :: h(:), kinetic(:), overlap(:)
      TYPE(dbcsr_type)                                   :: imaginary, real_part
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_grid_type), POINTER                        :: grid
      TYPE(pw_pool_type), POINTER                        :: pool
      TYPE(pw_r3d_rs_type)                               :: local_r, unit_tau, zero_density, zeros(3)
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state

      CALL get_ks_env(ks_env, upf_inverse=state)
      CPASSERT(ASSOCIATED(state))
      CALL get_ks_env(ks_env, upf_inverse=state, &
                      matrix_h=h, kinetic=kinetic, matrix_s=overlap)
      CPASSERT(ASSOCIATED(state))
      CALL check_geometry(ks_env, state)
      CALL get_ks_env(ks_env, pw_env=pw_env)
      CALL pw_env_get(pw_env, auxbas_pw_pool=pool)
      grid => pool%pw_grid
      CALL local_r%create(grid)
      CALL unit_tau%create(grid)
      CALL zero_density%create(grid)
      unit_tau%array = 1.0_dp
      zero_density%array = 0.0_dp
      DO axis = 1, 3
         CALL zeros(axis)%create(grid)
         zeros(axis)%array = 0.0_dp
      END DO
      IF (state%frame%ao%atom_local) THEN
         CALL qs_upf_inverse_atomic_operators(state%frame%ao, state%frame%ao%radial_q, atomic_t, atomic_m)
         atomic_t(:, :) = atomic_t + state%gaussian_kinetic
         CALL qs_upf_inverse_matrix_scatter(state%frame%ao, atomic_t, overlap(1)%matrix, real_part, imaginary)
         unit_tau%array = 0.0_dp
      ELSE
         CALL qs_upf_inverse_grid_integrate(state%frame, zero_density, zeros, unit_tau, overlap(1)%matrix, &
                                            real_part, imaginary)
      END IF
      CALL check_real(state, real_part, imaginary)
      CALL dbcsr_copy(kinetic(1)%matrix, real_part)
      CALL dbcsr_release(real_part)
      CALL dbcsr_release(imaginary)
      CALL local_potential(ks_env, grid, state%frame%ao%positions, local_r)
      CALL qs_upf_inverse_grid_integrate(state%frame, local_r, zeros, unit_tau, overlap(1)%matrix, &
                                         real_part, imaginary, projector_operator=state%frame%ao%dion, &
                                         augmentation=state%augmentation, fourier_charge=.TRUE.)
      CALL check_real(state, real_part, imaginary)
      CALL dbcsr_copy(h(1)%matrix, real_part)
      IF (state%frame%ao%atom_local) CALL dbcsr_add(h(1)%matrix, kinetic(1)%matrix, 1.0_dp, 1.0_dp)
      CALL dbcsr_release(real_part)
      CALL dbcsr_release(imaginary)
      CALL local_r%release()
      CALL unit_tau%release()
      CALL zero_density%release()
      DO axis = 1, 3
         CALL zeros(axis)%release()
      END DO
   END SUBROUTINE qs_upf_inverse_scf_core

! **************************************************************************************************
!> \brief Native spin density update with original Fourier Q and the same charge metric.
!> \param ks_env ...
!> \param rho_struct Native AO/grid density with occupations in its density matrices.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_scf_collocate(ks_env, rho_struct)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(qs_rho_type), INTENT(INOUT)                   :: rho_struct

      INTEGER                                            :: axis, spin
      REAL(dp), POINTER                                  :: total(:)
      TYPE(dbcsr_p_type), POINTER                        :: density(:, :)
      TYPE(pw_c1d_gs_type), POINTER                      :: rho_g(:)
      TYPE(pw_r3d_rs_type)                               :: gradient(3), tau
      TYPE(pw_r3d_rs_type), POINTER                      :: rho_r(:)
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state

      CALL get_ks_env(ks_env, upf_inverse=state)
      CPASSERT(ASSOCIATED(state))
      CALL check_geometry(ks_env, state)
      CALL qs_rho_get(rho_struct, rho_ao_kp=density, rho_r=rho_r, rho_g=rho_g, tot_rho_r=total)
      CPASSERT(SIZE(density, 2) == 1)
      CALL tau%create(rho_r(1)%pw_grid)
      DO axis = 1, 3
         CALL gradient(axis)%create(rho_r(1)%pw_grid)
      END DO
      DO spin = 1, SIZE(density, 1)
         CALL qs_upf_inverse_grid_collocate(state%frame, density(spin, 1)%matrix, rho_r(spin), gradient, tau, &
                                            augmentation=state%augmentation, fourier_charge=.TRUE.)
         CALL pw_transfer(rho_r(spin), rho_g(spin))
         total(spin) = pw_integrate_function(rho_r(spin), isign=-1)
      END DO
      CALL qs_rho_set(rho_struct, rho_r_valid=.TRUE., rho_g_valid=.TRUE.)
      CALL tau%release()
      DO axis = 1, 3
         CALL gradient(axis)%release()
      END DO
      state%density_evaluations = state%density_evaluations + 1
   END SUBROUTINE qs_upf_inverse_scf_collocate

! **************************************************************************************************
!> \brief Integrate the native folded Hartree/XC potential through the same inverse density adjoint.
!> \param ks_env ...
!> \param potential Native potential, already including the spectral XC gradient adjoint.
!> \param hamiltonian Native spin Hamiltonian, incremented.
!> \param weighted Whether potential already contains the native dvol factor.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_scf_integrate(ks_env, potential, hamiltonian, weighted)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: potential
      TYPE(dbcsr_type), INTENT(INOUT)                    :: hamiltonian
      LOGICAL, INTENT(IN)                                :: weighted

      INTEGER                                            :: axis
      TYPE(dbcsr_type)                                   :: imaginary, real_part
      TYPE(pw_r3d_rs_type)                               :: unweighted, zero_tau, zeros(3)
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state

      CALL get_ks_env(ks_env, upf_inverse=state)
      CPASSERT(ASSOCIATED(state))
      CALL check_geometry(ks_env, state)
      CALL unweighted%create(potential%pw_grid)
      unweighted%array = potential%array
      IF (weighted) CALL pw_scale(unweighted, 1.0_dp/potential%pw_grid%dvol)
      CALL zero_tau%create(potential%pw_grid)
      zero_tau%array = 0.0_dp
      DO axis = 1, 3
         CALL zeros(axis)%create(potential%pw_grid)
         zeros(axis)%array = 0.0_dp
      END DO
      CALL qs_upf_inverse_grid_integrate(state%frame, unweighted, zeros, zero_tau, hamiltonian, &
                                         real_part, imaginary, augmentation=state%augmentation, fourier_charge=.TRUE.)
      CALL check_real(state, real_part, imaginary)
      CALL dbcsr_add(hamiltonian, real_part, 1.0_dp, 1.0_dp)
      CALL dbcsr_release(real_part)
      CALL dbcsr_release(imaginary)
      CALL unweighted%release()
      CALL zero_tau%release()
      DO axis = 1, 3
         CALL zeros(axis)%release()
      END DO
      state%operator_evaluations = state%operator_evaluations + 1
   END SUBROUTINE qs_upf_inverse_scf_integrate

! **************************************************************************************************
!> \brief Joint geometry response of the inverse core Hamiltonian and its charge metric.
!>        Differentiates soft kinetic energy, original Dion, all short-range ionic potentials and Q.
!>        Gaussian ionic/core electrostatics remains in the native Poisson and core-energy terms.
!> \param ks_env Native kinds and their original radial potentials; geometry is not modified.
!> \param frame Prepared inverse fields at the geometry being differentiated.
!> \param augmentation Original charge kernels shared with the density map.
!> \param grid Native physical-cell grid of frame.
!> \param atom_response Cartesian response of every native atom.
!> \param cell_response Lattice-vector response, including its reciprocal and volume effects.
!> \param reduced_k_response Reduced Bloch-vector response, zero for cell stress.
!> \param core_response New replicated core-Hamiltonian response without occupation/k weights.
!> \param metric_response New response of the same radial-Q metric used by inverse SCF.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_scf_core_response(ks_env, frame, augmentation, grid, atom_response, cell_response, &
                                               reduced_k_response, core_response, metric_response)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(qs_upf_inverse_grid_type), INTENT(IN)         :: frame
      TYPE(qs_upf_inverse_augmentation_type), &
         INTENT(INOUT)                                   :: augmentation
      TYPE(pw_grid_type), POINTER                        :: grid
      REAL(dp), INTENT(IN)                               :: atom_response(:, :), &
                                                            cell_response(3, 3), &
                                                            reduced_k_response(3)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: core_response(:, :), &
                                                            metric_response(:, :)

      COMPLEX(dp), ALLOCATABLE                           :: atomic_m(:, :), atomic_t(:, :), &
                                                            dm(:, :), dt(:, :)
      INTEGER                                            :: axis, info
      REAL(dp)                                           :: cell_inverse(3, 3), strain(3, 3)
      TYPE(pw_r3d_rs_type)                               :: local_r, local_response, unit_tau, zero, &
                                                            zeros(3)

      cell_inverse = frame%ao%cell
      CALL invmat(cell_inverse, info)
      CPASSERT(info == 0)
      strain = MATMUL(cell_response, cell_inverse)
      CALL local_r%create(grid)
      CALL local_response%create(grid)
      CALL unit_tau%create(grid)
      CALL zero%create(grid)
      unit_tau%array = 1.0_dp
      IF (frame%ao%atom_local) unit_tau%array = 0.0_dp
      zero%array = 0.0_dp
      DO axis = 1, 3
         CALL zeros(axis)%create(grid)
         zeros(axis)%array = 0.0_dp
      END DO
      CALL local_potential(ks_env, grid, frame%ao%positions, local_r, local_response, atom_response, strain)
      CALL qs_upf_inverse_grid_response(frame, augmentation, frame%ao%radial_q, local_r, zeros, unit_tau, &
                                        atom_response, cell_response, reduced_k_response, local_response, zeros, zero, &
                                        core_response, metric_response, frame%ao%dion)
      IF (frame%ao%atom_local) THEN
         CALL qs_upf_inverse_atomic_operators(frame%ao, frame%ao%radial_q, atomic_t, atomic_m, dt, dm, &
                                              atom_response, cell_response, reduced_k_response)
         core_response(:, :) = core_response + dt
      END IF
      CALL local_r%release()
      CALL local_response%release()
      CALL unit_tau%release()
      CALL zero%release()
      DO axis = 1, 3
         CALL zeros(axis)%release()
      END DO
   END SUBROUTINE qs_upf_inverse_scf_core_response

! **************************************************************************************************
!> \brief Contract joint inverse operator and metric responses with native stationary P and W.
!>        W must have been prepared by the native energy-weighted-density service after SCF.
!>        Each k-point mesh weight is applied once, after replicated AO contractions.
!>        Ionic/core, Poisson-kernel and explicit XC measure responses remain separate energy terms.
!> \param ks_env Converged native inverse SCF state and prepared energy-weighted densities.
!> \param vdensity Unweighted primitive spin potentials conjugate to density.
!> \param vgradient Unweighted primitive spin potentials conjugate to Cartesian gradients.
!> \param vtau Unweighted primitive spin potentials conjugate to smooth kinetic density.
!> \param atom_response Cartesian response of all native atoms.
!> \param cell_response Lattice-vector response, with fixed fractional grid rows.
!> \param reduced_k_response Response of all actual reduced k vectors, or one column at Gamma.
!> \param response Stationary electronic response, including the metric constraint.
!> \param metric_response Signed metric contribution, already included in response.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_scf_response(ks_env, vdensity, vgradient, vtau, atom_response, cell_response, &
                                          reduced_k_response, response, metric_response)
      TYPE(qs_ks_env_type), INTENT(IN), TARGET           :: ks_env
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: vdensity(:), vgradient(:, :), vtau(:)
      REAL(dp), INTENT(IN)                               :: atom_response(:, :), &
                                                            cell_response(3, 3), &
                                                            reduced_k_response(:, :)
      REAL(dp), INTENT(OUT)                              :: response, metric_response

      COMPLEX(dp), ALLOCATABLE :: core(:, :, :), density(:, :, :), energy_weighted(:, :, :), &
         matrix(:, :), metric(:, :, :), metric_part(:, :)
      INTEGER                                            :: axis, ik, nk, nspin
      REAL(dp)                                           :: weight
      TYPE(dbcsr_p_type), POINTER                        :: density_gamma(:, :), energy_gamma(:)
      TYPE(pw_grid_type), POINTER                        :: grid
      TYPE(pw_r3d_rs_type)                               :: zero, zeros(3)
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state

      CALL get_ks_env(ks_env, upf_inverse=state, rho=rho)
      CPASSERT(ASSOCIATED(state))
      CALL check_geometry(ks_env, state)
      nspin = SIZE(vdensity)
      CPASSERT(ALL(SHAPE(vgradient) == [3, nspin]) .AND. SIZE(vtau) == nspin)
      grid => vdensity(1)%pw_grid
      CPASSERT(ASSOCIATED(grid))
      CALL zero%create(grid)
      zero%array = 0.0_dp
      DO axis = 1, 3
         CALL zeros(axis)%create(grid)
         zeros(axis)%array = 0.0_dp
      END DO
      nk = 1
      CPASSERT(ALL(SHAPE(reduced_k_response) == [3, nk]))
      DO ik = 1, nk
         CALL qs_upf_inverse_scf_core_response(ks_env, state%frame, state%augmentation, grid, atom_response, &
                                               cell_response, reduced_k_response(:, ik), matrix, metric_part)
         IF (ik == 1) ALLOCATE (core(SIZE(matrix, 1), SIZE(matrix, 2), nk), metric(SIZE(matrix, 1), SIZE(matrix, 2), nk))
         core(:, :, ik) = matrix
         metric(:, :, ik) = metric_part
      END DO
      response = 0.0_dp
      metric_response = 0.0_dp
      CALL qs_rho_get(rho, rho_ao_kp=density_gamma)
      CALL get_ks_env(ks_env, matrix_w=energy_gamma)
      CPASSERT(ASSOCIATED(density_gamma) .AND. ASSOCIATED(energy_gamma))
      CPASSERT(SIZE(density_gamma, 1) == nspin .AND. SIZE(energy_gamma) == nspin)
      ALLOCATE (density(SIZE(core, 1), SIZE(core, 2), 1), energy_weighted(SIZE(core, 1), SIZE(core, 2), 1))
      DO axis = 1, nspin
         CALL qs_upf_inverse_matrix_gather(state%frame%ao, density_gamma(axis, 1)%matrix, matrix)
         density(:, :, 1) = matrix
         CALL qs_upf_inverse_matrix_gather(state%frame%ao, energy_gamma(axis)%matrix, matrix)
         energy_weighted(:, :, 1) = matrix
         IF (state%frame%ao%atom_local) THEN
            CALL gaussian_kinetic_response(ks_env, density_gamma(axis, 1)%matrix, &
                                           atom_response, cell_response, weight)
            response = response + weight
         END IF
         DO ik = 1, nk
            weight = 1.0_dp
            CALL qs_upf_inverse_grid_response(state%frame, state%augmentation, state%frame%ao%radial_q, &
                                              vdensity(axis), vgradient(:, axis), vtau(axis), atom_response, cell_response, &
                                              reduced_k_response(:, ik), zero, zeros, zero, matrix, metric_part)
            response = response + weight*REAL(SUM(density(:, :, ik)*TRANSPOSE(core(:, :, ik) + matrix)), dp)
            metric_response = metric_response - weight*REAL(SUM(energy_weighted(:, :, ik)*TRANSPOSE(metric(:, :, ik))), dp)
         END DO
      END DO
      response = response + metric_response
      CALL zero%release()
      DO axis = 1, 3
         CALL zeros(axis)%release()
      END DO
   END SUBROUTINE qs_upf_inverse_scf_response

! **************************************************************************************************
!> \brief Contract the native Gamma Gaussian kinetic derivative without modifying force owners.
!> \param ks_env Prepared native bases and neighbor lists.
!> \param density Accepted native AO occupations for one spin.
!> \param atom_response Cartesian displacement of every atom.
!> \param cell_response Lattice-vector response.
!> \param response Replicated kinetic energy response.
! **************************************************************************************************
   SUBROUTINE gaussian_kinetic_response(ks_env, density, atom_response, cell_response, response)
      TYPE(qs_ks_env_type), INTENT(IN), TARGET           :: ks_env
      TYPE(dbcsr_type), INTENT(IN), TARGET               :: density
      REAL(dp), INTENT(IN)                               :: atom_response(:, :), cell_response(3, 3)
      REAL(dp), INTENT(OUT)                              :: response

      INTEGER                                            :: info
      REAL(dp)                                           :: cell_inverse(3, 3), strain(3, 3), &
                                                            virial(3, 3)
      REAL(dp), ALLOCATABLE                              :: gradients(:, :), relative_motion(:, :)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), POINTER                        :: operators(:)
      TYPE(dbcsr_type), POINTER                          :: occupations
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), POINTER            :: sab(:)
      TYPE(particle_type), POINTER                       :: atoms(:)
      TYPE(qs_ks_env_type), POINTER                      :: native_ks

      native_ks => ks_env
      occupations => density
      CALL get_ks_env(ks_env, cell=cell, particle_set=atoms, sab_orb=sab, para_env=para_env)
      CPASSERT(ASSOCIATED(sab))
      ALLOCATE (gradients(3, SIZE(atoms)), relative_motion(3, SIZE(atoms)))
      cell_inverse = cell%hmat
      CALL invmat(cell_inverse, info)
      CPASSERT(info == 0)
      strain = MATMUL(cell_response, cell_inverse)
      DO info = 1, SIZE(atoms)
         relative_motion(:, info) = atom_response(:, info) - MATMUL(strain, atoms(info)%r)
      END DO
      NULLIFY (operators)
      CALL build_kinetic_matrix(native_ks, matrix_t=operators, basis_type="ORB", sab_nl=sab, &
                                calculate_forces=.TRUE., matrix_p=occupations, &
                                force_response=gradients, virial_response=virial)
      CALL para_env%sum(gradients)
      CALL para_env%sum(virial)
      response = SUM(gradients*relative_motion) - SUM(virial*strain)
      CALL dbcsr_deallocate_matrix_set(operators)
   END SUBROUTINE gaussian_kinetic_response

! **************************************************************************************************
!> \brief Original short-range ionic potential and its optional fixed-fractional-row response.
!>        The same Fourier-to-real projection is used for the core operator and its derivative.
!> \param ks_env Native kind mapping and prepared original UPF radial potentials.
!> \param grid Native physical-cell grid.
!> \param positions Cartesian centers at the current inverse geometry.
!> \param local_r Allocated unweighted scalar potential, overwritten.
!> \param local_response Optional allocated response, including the inverse-volume derivative.
!> \param atom_response Required with local_response, Cartesian center displacements.
!> \param strain Required with local_response, cell_response times inverse cell.
! **************************************************************************************************
   SUBROUTINE local_potential(ks_env, grid, positions, local_r, local_response, atom_response, strain)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(pw_grid_type), POINTER                        :: grid
      REAL(dp), INTENT(IN)                               :: positions(:, :)
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: local_r
      TYPE(pw_r3d_rs_type), INTENT(INOUT), OPTIONAL      :: local_response
      REAL(dp), INTENT(IN), OPTIONAL                     :: atom_response(:, :), strain(3, 3)

      COMPLEX(dp), ALLOCATABLE                           :: center_derivative(:, :), ionic_g(:), &
                                                            reciprocal_derivative(:, :)
      INTEGER                                            :: atom, axis, ikind, ng
      REAL(dp)                                           :: volume_response
      REAL(dp), ALLOCATABLE                              :: g_response(:, :)
      TYPE(particle_type), POINTER                       :: atoms(:)
      TYPE(pw_c1d_gs_type)                               :: local_g, response_g
      TYPE(qs_kind_type), POINTER                        :: kinds(:)

      CALL get_ks_env(ks_env, qs_kind_set=kinds, particle_set=atoms)
      CPASSERT(ALL(SHAPE(positions) == [3, SIZE(atoms)]))
      CPASSERT(ASSOCIATED(local_r%pw_grid, grid))
      CALL local_g%create(grid)
      ng = SIZE(local_g%array)
      ALLOCATE (ionic_g(ng))
      local_g%array = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (PRESENT(local_response)) THEN
         CPASSERT(PRESENT(atom_response) .AND. PRESENT(strain))
         CPASSERT(ALL(SHAPE(atom_response) == SHAPE(positions)))
         CPASSERT(ASSOCIATED(local_response%pw_grid, grid))
         CALL response_g%create(grid)
         response_g%array = CMPLX(0.0_dp, 0.0_dp, dp)
         ALLOCATE (center_derivative(3, ng), reciprocal_derivative(3, ng), g_response(3, ng))
         g_response(:, :) = -MATMUL(TRANSPOSE(strain), grid%g)
         volume_response = strain(1, 1) + strain(2, 2) + strain(3, 3)
      END IF
      DO atom = 1, SIZE(atoms)
         ikind = atoms(atom)%atomic_kind%kind_number
         IF (kinds(ikind)%ghost) CYCLE
         ASSOCIATE (potential => kinds(ikind)%sgp_potential)
            IF (potential%upf_local%radial%degree < 0) THEN
               CALL upf_local_init(potential%upf_data, potential%alpha_core_charge, 0, potential%upf_local)
            END IF
            IF (PRESENT(local_response)) THEN
               CALL upf_local_fourier(potential%upf_local, grid%g, positions(:, atom), grid%vol, ionic_g, &
                                      center_derivative, reciprocal_derivative)
               response_g%array = response_g%array - volume_response*ionic_g
               DO axis = 1, 3
                  response_g%array = response_g%array + atom_response(axis, atom)*center_derivative(axis, :) + &
                                     g_response(axis, :)*reciprocal_derivative(axis, :)
               END DO
            ELSE
               CALL upf_local_fourier(potential%upf_local, grid%g, positions(:, atom), grid%vol, ionic_g)
            END IF
         END ASSOCIATE
         local_g%array = local_g%array + ionic_g
      END DO
      CALL pw_transfer(local_g, local_r)
      IF (PRESENT(local_response)) THEN
         CALL pw_transfer(response_g, local_response)
         CALL response_g%release()
      END IF
      CALL local_g%release()
   END SUBROUTINE local_potential

! **************************************************************************************************
!> \brief Reject a stale inverse geometry instead of mixing it with rebuilt Gaussian operators.
!> \param ks_env ...
!> \param state ...
! **************************************************************************************************
   SUBROUTINE check_geometry(ks_env, state)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(qs_upf_inverse_scf_type), INTENT(IN)          :: state

      INTEGER                                            :: atom
      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), POINTER                       :: atoms(:)

      CALL get_ks_env(ks_env, cell=cell, particle_set=atoms)
      CPASSERT(state%enabled)
      CPASSERT(ALL(cell%hmat == state%frame%ao%cell))
      DO atom = 1, SIZE(atoms)
         CPASSERT(ALL(atoms(atom)%r == state%frame%ao%positions(:, atom)))
      END DO
   END SUBROUTINE check_geometry

! **************************************************************************************************
!> \brief Verify that the Gamma operator can be stored in the real native KS representation.
!> \param state ...
!> \param real_part ...
!> \param imaginary ...
! **************************************************************************************************
   SUBROUTINE check_real(state, real_part, imaginary)
      TYPE(qs_upf_inverse_scf_type), INTENT(IN)          :: state
      TYPE(dbcsr_type), INTENT(IN)                       :: real_part, imaginary

      COMPLEX(dp), ALLOCATABLE                           :: dense(:, :)
      REAL(dp)                                           :: defect, tolerance

      CALL qs_upf_inverse_matrix_gather(state%frame%ao, real_part, dense, imaginary)
      defect = MAXVAL(ABS(AIMAG(dense)))
      tolerance = 128.0_dp*EPSILON(1.0_dp)*MAX(1.0_dp, MAXVAL(ABS(dense)))
      CPASSERT(defect <= tolerance)
   END SUBROUTINE check_real

! **************************************************************************************************
!> \brief Drop the metric, fields and allocation pattern together before a native geometry rebuild.
!> \param ks_env Native environment. The input kind selection remains available for rebuilding.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_scf_reset(ks_env)
      TYPE(qs_ks_env_type), INTENT(INOUT)                :: ks_env

      TYPE(qs_upf_inverse_scf_type), POINTER             :: state
      TYPE(qs_uspp_type), POINTER                        :: uspp

      CALL get_ks_env(ks_env, upf_inverse=state, uspp=uspp)
      CALL qs_upf_inverse_scf_release(state)
      CALL qs_uspp_release(uspp)
      CALL set_ks_env(ks_env, upf_inverse=state, uspp=uspp)
   END SUBROUTINE qs_upf_inverse_scf_reset

! **************************************************************************************************
!> \brief Activate scalar USPP GPW from the kind potentials before allocating SCF matrices.
!>        Image bounds cover the native basis radii and original radial support on the current cell.
!>        A geometry update rebuilds both the charge metric and all density/operator fields.
!> \param ks_env Ordinary Gaussian core and PW environment at the new geometry.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_scf_initialize(ks_env)
      TYPE(qs_ks_env_type), INTENT(INOUT)                :: ks_env

      INTEGER                                            :: atom, axis, ikind, image, info, ix, iy, &
                                                            iz, lower(3), upper(3)
      INTEGER, ALLOCATABLE                               :: images(:, :)
      REAL(dp)                                           :: fractional(3), lower_real(3), radius, &
                                                            reach, upper_real(3)
      REAL(dp), ALLOCATABLE                              :: radial_radius(:)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dft_control_type), POINTER                    :: control
      TYPE(gto_basis_set_type), POINTER                  :: basis
      TYPE(particle_type), POINTER                       :: atoms(:)
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_pool_type), POINTER                        :: pool
      TYPE(qs_kind_type), POINTER                        :: kinds(:)
      TYPE(qs_uspp_type), POINTER                        :: uspp
      TYPE(upf_radial_fields_type)                       :: beta, delta

      CALL get_ks_env(ks_env, qs_kind_set=kinds, particle_set=atoms, cell=cell, pw_env=pw_env, dft_control=control)
      CALL get_ks_env(ks_env, uspp=uspp)
      CALL qs_uspp_init(uspp, kinds, atoms)
      CALL set_ks_env(ks_env, uspp=uspp)
      IF (.NOT. ASSOCIATED(uspp)) RETURN
      ALLOCATE (radial_radius(SIZE(kinds)))
      radial_radius = 0.0_dp
      DO ikind = 1, SIZE(kinds)
         IF (kinds(ikind)%ghost) CYCLE
         IF (.NOT. ASSOCIATED(kinds(ikind)%sgp_potential)) THEN
            CPABORT("Inverse USPP mixed potentials are not implemented")
         END IF
         IF (.NOT. ALLOCATED(kinds(ikind)%sgp_potential%upf_data)) THEN
            CPABORT("Inverse USPP mixed potentials are not implemented")
         END IF
         IF (.NOT. kinds(ikind)%sgp_potential%upf_data%is_ultrasoft) THEN
            CPABORT("Inverse USPP mixed potentials are not implemented")
         END IF
         CALL upf_inverse_fields_init(kinds(ikind)%sgp_potential%upf_data, beta, delta, &
                                      control%qs_control%uspp_method == uspp_method_direct)
         radial_radius(ikind) = MAX(beta%radius, delta%radius)
      END DO
      IF (ANY(cell%perd == 0)) CPABORT("Inverse USPP GPW currently requires a fully periodic cell")
      lower_real = 0.0_dp
      upper_real = 0.0_dp
      DO atom = 1, SIZE(atoms)
         ikind = atoms(atom)%atomic_kind%kind_number
         CALL get_qs_kind(kinds(ikind), basis_set=basis, basis_type='ORB')
         radius = basis%kind_radius
         IF (.NOT. kinds(ikind)%ghost) THEN
            radius = MAX(radius, radial_radius(ikind))
         END IF
         fractional = MATMUL(cell%h_inv, atoms(atom)%r)
         DO axis = 1, 3
            reach = radius*SQRT(SUM(cell%h_inv(axis, :)**2))
            ! One additional cell keeps exponentially small Gaussian tails away from image boundaries.
            lower_real(axis) = MIN(lower_real(axis), -fractional(axis) - reach - 1.0_dp)
            upper_real(axis) = MAX(upper_real(axis), 1.0_dp - fractional(axis) + reach + 1.0_dp)
         END DO
      END DO
      lower = FLOOR(lower_real)
      upper = CEILING(upper_real)
      ALLOCATE (images(3, PRODUCT(upper - lower + 1)))
      image = 0
      DO iz = lower(3), upper(3)
         DO iy = lower(2), upper(2)
            DO ix = lower(1), upper(1)
               image = image + 1
               images(:, image) = [ix, iy, iz]
            END DO
         END DO
      END DO
      CALL pw_env_get(pw_env, auxbas_pw_pool=pool)
      SELECT CASE (control%qs_control%uspp_method)
      CASE (uspp_method_direct, uspp_method_inverse, uspp_method_inverse_1c)
         CALL qs_upf_inverse_scf_prepare(ks_env, pool%pw_grid, images, info, &
                                         control%qs_control%uspp_method /= uspp_method_inverse, &
                                         control%qs_control%uspp_method == uspp_method_direct)
      CASE DEFAULT
         CPABORT("Unknown USPP orbital representation")
      END SELECT
      IF (info == -2) THEN
         CPABORT("Overlapping USPP reconstruction regions. Use USPP_METHOD INVERSE")
      END IF
      IF (info /= 0) CPABORT("Inverse USPP reconstruction is singular or ill-conditioned")
   END SUBROUTINE qs_upf_inverse_scf_initialize

END MODULE qs_upf_inverse_scf
