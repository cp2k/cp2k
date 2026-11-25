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
                                              dbcsr_scale,&
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
   USE qs_collocate_density,            ONLY: calculate_rho_elec
   USE qs_integrate_potential_product,  ONLY: integrate_v_rspace
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
   USE qs_upf_inverse_core,             ONLY: qs_upf_inverse_core_collocate
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
      TYPE(pw_r3d_rs_type)                               :: unit_potential
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
      IF (PRESENT(direct)) state%native_gaussian = direct
      IF (state%native_gaussian) THEN
         CALL direct_prepare(ks_env, state)
      END IF
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
      CALL get_ks_env(ks_env, kinetic=kinetic)
      CALL qs_upf_inverse_matrix_gather(state%frame%ao, kinetic(1)%matrix, state%gaussian_kinetic)
      CALL qs_upf_inverse_augmentation_prepare(ks_env, state%frame%ao, state%augmentation)
      ALLOCATE (metric(1))
      ALLOCATE (metric(1)%matrix)
      state%frame%native_gaussian = state%native_gaussian
      CALL qs_upf_inverse_grid_metric(state%frame, state%frame%ao%radial_q, overlap(1)%matrix, &
                                      metric(1)%matrix, imaginary)
      IF (state%native_gaussian) THEN
         CALL unit_potential%create(grid)
         unit_potential%array = 1.0_dp
         CALL direct_integrate(ks_env, state, unit_potential)
         CALL dbcsr_add(metric(1)%matrix, state%direct_potential(1)%matrix, 1.0_dp, 1.0_dp)
         CALL unit_potential%release()
      END IF
      CALL check_real(state, metric(1)%matrix, imaginary)
      CALL dbcsr_release(imaginary)
      CALL set_1d_pointer(state%metric, metric)
      ! Retain all inverse-mediated AO fill in native density/KS allocation as well as in M.
      CALL dbcsr_copy(uspp%ao_pattern, metric(1)%matrix)
      CALL dbcsr_set(uspp%ao_pattern, 0.0_dp)
      state%enabled = .TRUE.
      CALL set_ks_env(ks_env, upf_inverse=state)
      CALL qs_upf_inverse_core_collocate(ks_env)
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
      TYPE(pw_r3d_rs_type)                               :: local_r, unit_tau, zeros(3)
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
      unit_tau%array = 0.0_dp
      DO axis = 1, 3
         CALL zeros(axis)%create(grid)
         zeros(axis)%array = 0.0_dp
      END DO
      CALL qs_upf_inverse_atomic_operators(state%frame%ao, state%frame%ao%radial_q, atomic_t, atomic_m)
      atomic_t(:, :) = atomic_t + state%gaussian_kinetic
      CALL qs_upf_inverse_matrix_scatter(state%frame%ao, atomic_t, overlap(1)%matrix, real_part, imaginary)
      CALL check_real(state, real_part, imaginary)
      CALL dbcsr_copy(kinetic(1)%matrix, real_part)
      CALL dbcsr_release(real_part)
      CALL dbcsr_release(imaginary)
      CALL local_potential(ks_env, grid, state%frame%ao%positions, local_r)
      CALL qs_upf_inverse_grid_integrate(state%frame, local_r, zeros, unit_tau, overlap(1)%matrix, &
                                         real_part, imaginary, projector_operator=state%frame%ao%dion, &
                                         augmentation=state%augmentation, fourier_charge=.TRUE.)
      CALL check_real(state, real_part, imaginary)
      IF (state%native_gaussian) THEN
         CALL direct_integrate(ks_env, state, local_r)
         CALL dbcsr_add(real_part, state%direct_potential(1)%matrix, 1.0_dp, 1.0_dp)
      END IF
      CALL dbcsr_copy(h(1)%matrix, real_part)
      CALL dbcsr_add(h(1)%matrix, kinetic(1)%matrix, 1.0_dp, 1.0_dp)
      CALL dbcsr_release(real_part)
      CALL dbcsr_release(imaginary)
      CALL local_r%release()
      CALL unit_tau%release()
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
      TYPE(qs_ks_env_type), INTENT(IN), TARGET           :: ks_env
      TYPE(qs_rho_type), INTENT(INOUT)                   :: rho_struct

      INTEGER                                            :: axis, spin
      REAL(dp), POINTER                                  :: total(:)
      TYPE(dbcsr_p_type), POINTER                        :: density(:, :)
      TYPE(pw_c1d_gs_type), POINTER                      :: rho_g(:)
      TYPE(pw_r3d_rs_type)                               :: gaussian_rho, gradient(3), tau
      TYPE(pw_r3d_rs_type), POINTER                      :: rho_r(:)
      TYPE(qs_ks_env_type), POINTER                      :: native_ks
      TYPE(qs_upf_inverse_scf_type), POINTER             :: state

      CALL get_ks_env(ks_env, upf_inverse=state)
      CPASSERT(ASSOCIATED(state))
      CALL check_geometry(ks_env, state)
      CALL qs_rho_get(rho_struct, rho_ao_kp=density, rho_r=rho_r, rho_g=rho_g, tot_rho_r=total)
      CPASSERT(SIZE(density, 2) == 1)
      native_ks => ks_env
      IF (state%native_gaussian) CALL gaussian_rho%create(rho_r(1)%pw_grid)
      CALL tau%create(rho_r(1)%pw_grid)
      DO axis = 1, 3
         CALL gradient(axis)%create(rho_r(1)%pw_grid)
      END DO
      DO spin = 1, SIZE(density, 1)
         CALL qs_upf_inverse_grid_collocate(state%frame, density(spin, 1)%matrix, rho_r(spin), gradient, tau, &
                                            augmentation=state%augmentation, fourier_charge=.TRUE.)
         IF (state%native_gaussian) THEN
            CALL dbcsr_copy(state%direct_density(1)%matrix, density(spin, 1)%matrix)
            CALL calculate_rho_elec(matrix_p_kp=state%direct_density, rho=gaussian_rho, rho_gspace=rho_g(spin), &
                                    ks_env=native_ks, soft_valid=.FALSE.)
            rho_r(spin)%array = rho_r(spin)%array + gaussian_rho%array
         END IF
         CALL pw_transfer(rho_r(spin), rho_g(spin))
         total(spin) = pw_integrate_function(rho_r(spin), isign=-1)
      END DO
      CALL qs_rho_set(rho_struct, rho_r_valid=.TRUE., rho_g_valid=.TRUE.)
      IF (state%native_gaussian) CALL gaussian_rho%release()
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
      IF (state%native_gaussian) THEN
         CALL direct_integrate(ks_env, state, unweighted)
         CALL dbcsr_add(real_part, state%direct_potential(1)%matrix, 1.0_dp, 1.0_dp)
      END IF
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
                                                            dm(:, :), dt(:, :), dw(:, :)
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
      unit_tau%array = 0.0_dp

      zero%array = 0.0_dp
      DO axis = 1, 3
         CALL zeros(axis)%create(grid)
         zeros(axis)%array = 0.0_dp
      END DO
      CALL local_potential(ks_env, grid, frame%ao%positions, local_r, local_response, atom_response, strain)
      CALL qs_upf_inverse_grid_response(frame, augmentation, frame%ao%radial_q, local_r, zeros, unit_tau, &
                                        atom_response, cell_response, reduced_k_response, local_response, zeros, zero, &
                                        core_response, metric_response, frame%ao%dion, dw)
      CALL qs_upf_inverse_atomic_operators(frame%ao, frame%ao%radial_q, atomic_t, atomic_m, dt, dm, &
                                           atom_response, cell_response, reduced_k_response, dw)
      core_response(:, :) = core_response + dt
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
      INTEGER                                            :: axis, ik, inv_info, nk, nspin
      REAL(dp)                                           :: inverse_cell(3, 3), native_energy, &
                                                            native_metric, strain(3, 3), weight
      REAL(dp), ALLOCATABLE                              :: relative_motion(:, :)
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
      IF (state%native_gaussian) THEN
         CPASSERT(ALL(reduced_k_response == 0.0_dp))
         ALLOCATE (relative_motion, MOLD=atom_response)
         inverse_cell = state%frame%ao%cell
         CALL invmat(inverse_cell, inv_info)
         CPASSERT(inv_info == 0)
         strain = MATMUL(cell_response, inverse_cell)
         relative_motion(:, :) = atom_response - MATMUL(strain, state%frame%ao%positions)
      END IF
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
         CALL gaussian_kinetic_response(ks_env, density_gamma(axis, 1)%matrix, &
                                        atom_response, cell_response, weight)
         response = response + weight
         IF (state%native_gaussian) THEN
            CALL native_grid_response(ks_env, state, density, energy_weighted, vdensity(axis), vgradient(:, axis), vtau(axis), &
                                      atom_response, strain, relative_motion, native_energy, native_metric)
            response = response + native_energy
            metric_response = metric_response + native_metric
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
      REAL(dp) :: domain_lower(3), domain_upper(3), fractional(3), lower_real(3), radius, reach, &
         reciprocal_norm(3), upper_real(3)
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
      reciprocal_norm = SQRT(SUM(cell%h_inv**2, DIM=2))
      domain_lower = 0.0_dp
      domain_upper = 1.0_dp
      IF (control%qs_control%uspp_method == uspp_method_inverse) THEN
         DO atom = 1, SIZE(atoms)
            radius = radial_radius(atoms(atom)%atomic_kind%kind_number)
            fractional = MATMUL(cell%h_inv, atoms(atom)%r)
            domain_lower = MIN(domain_lower, fractional - radius*reciprocal_norm)
            domain_upper = MAX(domain_upper, fractional + radius*reciprocal_norm)
         END DO
      END IF
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
            reach = radius*reciprocal_norm(axis)
            ! One additional cell keeps exponentially small Gaussian tails away from image boundaries.
            lower_real(axis) = MIN(lower_real(axis), domain_lower(axis) - fractional(axis) - reach - 1.0_dp)
            upper_real(axis) = MAX(upper_real(axis), domain_upper(axis) - fractional(axis) + reach + 1.0_dp)
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

! **************************************************************************************************
!> \brief Allocate native matrix workspaces on the Gaussian overlap sparsity pattern.
!> \param ks_env Quickstep environment owning the overlap matrices.
!> \param state DIRECT density and potential workspaces.
! **************************************************************************************************
   SUBROUTINE direct_prepare(ks_env, state)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(qs_upf_inverse_scf_type), INTENT(INOUT)       :: state

      INTEGER                                            :: i
      TYPE(dbcsr_p_type), POINTER                        :: overlap(:, :)

      CALL get_ks_env(ks_env, matrix_s_kp=overlap)
      ALLOCATE (state%direct_density(SIZE(overlap, 2)), state%direct_potential(SIZE(overlap, 2)))
      DO i = 1, SIZE(overlap, 2)
         ALLOCATE (state%direct_density(i)%matrix, state%direct_potential(i)%matrix)
         CALL dbcsr_copy(state%direct_density(i)%matrix, overlap(1, i)%matrix)
         CALL dbcsr_copy(state%direct_potential(i)%matrix, overlap(1, i)%matrix)
      END DO
   END SUBROUTINE direct_prepare

! **************************************************************************************************
!> \brief Copy the spin density into the real-space Gaussian matrix blocks.
!> \param state DIRECT matrix workspace.
!> \param matrix Density matrices in the original AO order.
! **************************************************************************************************
   SUBROUTINE direct_set_density(state, matrix)
      TYPE(qs_upf_inverse_scf_type), INTENT(INOUT)       :: state
      COMPLEX(dp), INTENT(IN)                            :: matrix(:, :, :)

      TYPE(dbcsr_type)                                   :: imaginary, real_part

      CALL qs_upf_inverse_matrix_scatter(state%frame%ao, matrix(:, :, 1), state%direct_density(1)%matrix, real_part, imaginary)
      CALL dbcsr_copy(state%direct_density(1)%matrix, real_part)
      CALL dbcsr_release(real_part)
      CALL dbcsr_release(imaginary)
   END SUBROUTINE direct_set_density

! **************************************************************************************************
!> \brief Gather native Gaussian potential matrices in the orbital representation.
!> \param state DIRECT matrix workspace.
!> \param operators Dense potential matrices in the original AO order.
! **************************************************************************************************
   SUBROUTINE direct_operator(state, operators)
      TYPE(qs_upf_inverse_scf_type), INTENT(IN)          :: state
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: operators(:, :, :)

      COMPLEX(dp), ALLOCATABLE                           :: gamma(:, :)

      CALL qs_upf_inverse_matrix_gather(state%frame%ao, state%direct_potential(1)%matrix, gamma)
      ALLOCATE (operators(SIZE(gamma, 1), SIZE(gamma, 2), 1))
      operators(:, :, 1) = gamma
   END SUBROUTINE direct_operator

! **************************************************************************************************
!> \brief Apply the adjoint of the Gaussian density collocation to a scalar potential.
!> \param ks_env Quickstep environment owning the grid task lists.
!> \param state DIRECT matrix workspace for the result.
!> \param potential Scalar potential before multiplication by the integration measure.
! **************************************************************************************************
   SUBROUTINE direct_integrate(ks_env, state, potential)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(qs_upf_inverse_scf_type), INTENT(INOUT)       :: state
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: potential

      INTEGER                                            :: i

      DO i = 1, SIZE(state%direct_potential)
         CALL dbcsr_set(state%direct_potential(i)%matrix, 0.0_dp)
      END DO
      CALL integrate_v_rspace(potential, hmat_kp=state%direct_potential, calculate_forces=.FALSE., ks_env_external=ks_env)
      DO i = 1, SIZE(state%direct_potential)
         CALL dbcsr_scale(state%direct_potential(i)%matrix, potential%pw_grid%dvol)
      END DO
   END SUBROUTINE direct_integrate

! **************************************************************************************************
!> \brief Differentiate the native DIRECT Gaussian density map and its unit-potential metric.
!>        Native grid kernels supply the field derivatives. Include the integration measure and
!>        the explicit local-potential response before subtracting the normalization constraint.
!> \param ks_env Native owner of the same full Gaussian task lists used during SCF.
!> \param state Prepared DIRECT state with separate Gaussian and original augmentation maps.
!> \param density Accepted spin occupations, without k-point weights.
!> \param weighted Matching energy-weighted occupations, without k-point weights.
!> \param vdensity Folded density potential, including the spectral gradient adjoint.
!> \param vgradient Must be zero because native integration consumes the folded scalar potential.
!> \param vtau Potential conjugate to Gaussian kinetic density.
!> \param atom_response Cartesian displacements of the original atoms.
!> \param strain Physical deformation gradient of the cell.
!> \param relative_motion Atomic displacement relative to the moving fractional grid.
!> \param response Gaussian local and XC response, without the metric constraint.
!> \param metric_response Signed Gaussian normalization response.
! **************************************************************************************************
   SUBROUTINE native_grid_response(ks_env, state, density, weighted, vdensity, vgradient, vtau, &
                                   atom_response, strain, relative_motion, response, metric_response)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(qs_upf_inverse_scf_type), INTENT(INOUT)       :: state
      COMPLEX(dp), INTENT(IN)                            :: density(:, :, :), weighted(:, :, :)
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: vdensity, vgradient(3), vtau
      REAL(dp), INTENT(IN)                               :: atom_response(:, :), strain(3, 3), &
                                                            relative_motion(:, :)
      REAL(dp), INTENT(OUT)                              :: response, metric_response

      COMPLEX(dp), ALLOCATABLE                           :: operators(:, :, :)
      INTEGER                                            :: axis, ik
      REAL(dp)                                           :: energy, trace_strain, virial(3, 3), &
                                                            weights(SIZE(density, 3))
      REAL(dp), ALLOCATABLE                              :: forces(:, :)
      TYPE(pw_grid_type), POINTER                        :: grid
      TYPE(pw_r3d_rs_type)                               :: explicit_response, potential

      weights = 1.0_dp
      CPASSERT(state%native_gaussian)
      DO axis = 1, 3
         CPASSERT(ALL(vgradient(axis)%array == 0.0_dp))
      END DO
      grid => vdensity%pw_grid
      CALL potential%create(grid)
      CALL explicit_response%create(grid)
      ! Positions belong to the original AO ordering, including translated Bloch atoms.
      CALL local_potential(ks_env, grid, state%frame%ao%positions, potential, explicit_response, atom_response, strain)
      potential%array = potential%array + vdensity%array
      ALLOCATE (forces(3, SIZE(atom_response, 2)))
      trace_strain = strain(1, 1) + strain(2, 2) + strain(3, 3)
      CALL direct_geometry(ks_env, state, density, potential, forces, virial)
      CALL direct_operator(state, operators)
      energy = 0.0_dp
      DO ik = 1, SIZE(density, 3)
         energy = energy + weights(ik)*REAL(SUM(density(:, :, ik)*CONJG(operators(:, :, ik))), dp)
      END DO
      response = -SUM(forces*relative_motion) - SUM(virial*strain) + trace_strain*energy
      CALL direct_integrate(ks_env, state, explicit_response)
      CALL direct_operator(state, operators)
      DO ik = 1, SIZE(density, 3)
         response = response + weights(ik)*REAL(SUM(density(:, :, ik)*CONJG(operators(:, :, ik))), dp)
      END DO
      CPASSERT(ALL(vtau%array == 0.0_dp))
      potential%array = 1.0_dp
      CALL direct_geometry(ks_env, state, weighted, potential, forces, virial)
      CALL direct_operator(state, operators)
      energy = 0.0_dp
      DO ik = 1, SIZE(density, 3)
         energy = energy + weights(ik)*REAL(SUM(weighted(:, :, ik)*CONJG(operators(:, :, ik))), dp)
      END DO
      metric_response = SUM(forces*relative_motion) + SUM(virial*strain) - trace_strain*energy
      CALL potential%release()
      CALL explicit_response%release()
   END SUBROUTINE native_grid_response

! **************************************************************************************************
!> \brief Differentiate the Gaussian potential integral at fixed density.
!> \param ks_env Quickstep environment owning the grid task lists.
!> \param state DIRECT matrix workspace.
!> \param density Density matrices in the original AO order.
!> \param potential Scalar potential on the density grid.
!> \param forces Cartesian forces.
!> \param virial Cell virial.
! **************************************************************************************************
   SUBROUTINE direct_geometry(ks_env, state, density, potential, forces, virial)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(qs_upf_inverse_scf_type), INTENT(INOUT)       :: state
      COMPLEX(dp), INTENT(IN)                            :: density(:, :, :)
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: potential
      REAL(dp), INTENT(OUT)                              :: forces(:, :), virial(3, 3)

      INTEGER                                            :: i
      TYPE(pw_r3d_rs_type)                               :: weighted

      CALL direct_set_density(state, density)
      DO i = 1, SIZE(state%direct_potential)
         CALL dbcsr_set(state%direct_potential(i)%matrix, 0.0_dp)
      END DO
      CALL weighted%create(potential%pw_grid)
      weighted%array = potential%pw_grid%dvol*potential%array
      CALL integrate_v_rspace(weighted, hmat_kp=state%direct_potential, pmat_kp=state%direct_density, ks_env_external=ks_env, &
                              calculate_forces=.TRUE., force_response=forces, virial_response=virial)
      CALL potential%pw_grid%para%group%sum(forces)
      CALL potential%pw_grid%para%group%sum(virial)
      forces = -forces
      CALL weighted%release()
   END SUBROUTINE direct_geometry

END MODULE qs_upf_inverse_scf
