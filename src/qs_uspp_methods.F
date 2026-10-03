!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Connect native USPP geometry, density, and operators to Quickstep GPW.
! **************************************************************************************************
MODULE qs_uspp_methods
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_get_block_p, dbcsr_get_info, dbcsr_p_type, dbcsr_release, &
        dbcsr_scale, dbcsr_set, dbcsr_type
   USE input_constants,                 ONLY: do_ppl_analytic
   USE kinds,                           ONLY: dp
   USE kpoint_types,                    ONLY: kpoint_type
   USE message_passing,                 ONLY: mp_comm_type,&
                                              mp_comm_unequal
   USE particle_types,                  ONLY: particle_type
   USE pw_grid_types,                   ONLY: FULLSPACE,&
                                              HALFSPACE,&
                                              pw_grid_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_derive,&
                                              pw_scale,&
                                              pw_transfer
   USE pw_types,                        ONLY: pw_c1d_gs_type,&
                                              pw_r3d_rs_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type,&
                                              set_ks_env
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_uspp_projectors,              ONLY: qs_uspp_projector_gradient,&
                                              qs_uspp_projector_images
   USE qs_uspp_types,                   ONLY: qs_uspp_init,&
                                              qs_uspp_type,&
                                              qs_uspp_update,&
                                              qs_uspp_update_images
   USE sap_kind_types,                  ONLY: build_sap_ints,&
                                              release_sap_int,&
                                              sap_int_type,&
                                              sap_sort
   USE upf_local_integrals,             ONLY: upf_local_fourier,&
                                              upf_local_init
   USE uspp_augmentation,               ONLY: uspp_augmentation_collocate,&
                                              uspp_augmentation_integrate
   USE uspp_projector_operators,        ONLY: uspp_projector_add_operator,&
                                              uspp_projector_create_onsite,&
                                              uspp_projector_density,&
                                              uspp_projector_pullback
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_uspp_methods'
   PUBLIC :: qs_uspp_prepare, qs_uspp_init_ao_matrix, qs_uspp_collocate, qs_uspp_integrate
   PUBLIC :: qs_uspp_core_matrix, qs_uspp_core_forces

   INTERFACE qs_uspp_collocate
      MODULE PROCEDURE qs_uspp_collocate_single, qs_uspp_collocate_images
   END INTERFACE
   INTERFACE qs_uspp_integrate
      MODULE PROCEDURE qs_uspp_integrate_single, qs_uspp_integrate_images
   END INTERFACE
   INTERFACE qs_uspp_core_matrix
      MODULE PROCEDURE qs_uspp_core_matrix_single, qs_uspp_core_matrix_images
   END INTERFACE

CONTAINS

! **************************************************************************************************
!> \brief Prepare USPP centers, unweighted projector overlaps, and the generalized metric.
!>        Called after rebuilding the Gaussian overlap. The NC path allocates nothing.
!> \param ks_env KS environment with current kinds, particles and ABBA projector neighbor list.
!> \param gaussian_overlap Current Gaussian overlap; never modified or aliased by the metric.
! **************************************************************************************************
   SUBROUTINE qs_uspp_prepare(ks_env, gaussian_overlap)
      TYPE(qs_ks_env_type), INTENT(INOUT)                :: ks_env
      TYPE(dbcsr_p_type), INTENT(IN)                     :: gaussian_overlap(:)

      INTEGER                                            :: handle
      LOGICAL                                            :: do_kpoints
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(neighbor_list_set_p_type), POINTER            :: sab_kp(:), sap_ppnl(:)
      TYPE(particle_type), POINTER                       :: particle_set(:)
      TYPE(qs_kind_type), POINTER                        :: qs_kind_set(:)
      TYPE(qs_uspp_type), POINTER                        :: uspp
      TYPE(sap_int_type), POINTER                        :: sap_int(:)

      CALL get_ks_env(ks_env, uspp=uspp, qs_kind_set=qs_kind_set, particle_set=particle_set)
      CALL qs_uspp_init(uspp, qs_kind_set, particle_set)
      CALL set_ks_env(ks_env, uspp=uspp)
      IF (.NOT. ASSOCIATED(uspp)) RETURN
      CALL timeset('qs_uspp_prepare', handle)
      CALL get_ks_env(ks_env, do_kpoints=do_kpoints, sap_ppnl=sap_ppnl, dft_control=dft_control)
      IF (.NOT. dft_control%qs_control%gpw .OR. dft_control%qs_control%gapw .OR. &
          dft_control%qs_control%gapw_xc) THEN
         CPABORT("USPP currently requires GPW")
      END IF
      IF (dft_control%qs_control%do_ppl_method /= do_ppl_analytic) THEN
         CPABORT("USPP requires analytic local-potential integrals")
      END IF
      IF (dft_control%use_kinetic_energy_density) THEN
         CPABORT("USPP augmentation for kinetic-energy-density functionals is not implemented")
      END IF
      CPASSERT(ASSOCIATED(sap_ppnl))
      ALLOCATE (sap_int(SIZE(qs_kind_set)**2))
      CALL build_sap_ints(sap_int, sap_ppnl, qs_kind_set, nder=0)
      IF (do_kpoints) THEN
         CALL get_ks_env(ks_env, kpoints=kpoints, sab_kp=sab_kp)
         CALL qs_uspp_update_images(uspp, gaussian_overlap, sap_int, particle_set, sab_kp, kpoints%cell_to_index)
      ELSE
         CALL qs_uspp_update(uspp, gaussian_overlap(1)%matrix, sap_int, particle_set)
         CALL release_sap_int(sap_int)
      END IF
      CALL timestop(handle)
   END SUBROUTINE qs_uspp_prepare

! **************************************************************************************************
!> \brief Initialize a matrix on the union of its ordinary and USPP projector-mediated blocks.
!>        Called after ordinary neighbor-list allocation, before assigning matrix values.
!> \param ks_env KS environment with current USPP geometry data.
!> \param matrix Newly allocated AO matrix, zeroed after expanding its block pattern.
! **************************************************************************************************
   SUBROUTINE qs_uspp_init_ao_matrix(ks_env, matrix)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix

      TYPE(qs_uspp_type), POINTER                        :: uspp

      CALL get_ks_env(ks_env, uspp=uspp)
      IF (.NOT. ASSOCIATED(uspp)) RETURN
      CPASSERT(ASSOCIATED(uspp%ao_pattern))
      CALL dbcsr_add(matrix, uspp%ao_pattern, 1.0_dp, 1.0_dp)
      CALL dbcsr_set(matrix, 0.0_dp)
   END SUBROUTINE qs_uspp_init_ao_matrix

! **************************************************************************************************
!> \brief Add USPP augmentation (or its spatial derivative) to a Quickstep density grid.
!>        The reciprocal and real representations are both updated. Ordinary kinds do nothing.
!> \param ks_env KS environment with the current projector overlaps.
!> \param density AO density for one spin, including occupations.
!> \param rho_g Existing reciprocal density, normalized as integral/cell volume.
!> \param rho_r Corresponding real density, overwritten from the updated reciprocal grid.
!> \param total_rho Optional signed electronic charge, decremented by the augmentation charge.
!> \param derivative Optional Cartesian derivative orders; total_rho must then be absent.
! **************************************************************************************************
   SUBROUTINE qs_uspp_collocate_images(ks_env, density, rho_g, rho_r, total_rho, derivative)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(dbcsr_p_type), INTENT(IN)                     :: density(:)
      TYPE(pw_c1d_gs_type), INTENT(INOUT)                :: rho_g
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: rho_r
      REAL(KIND=dp), INTENT(INOUT), OPTIONAL             :: total_rho
      INTEGER, INTENT(IN), OPTIONAL                      :: derivative(3)

      INTEGER                                            :: handle
      REAL(KIND=dp)                                      :: charge
      TYPE(pw_c1d_gs_type)                               :: augmentation_g
      TYPE(pw_r3d_rs_type)                               :: augmentation_r
      TYPE(qs_uspp_type), POINTER                        :: uspp

      CALL get_ks_env(ks_env, uspp=uspp)
      IF (.NOT. ASSOCIATED(uspp)) RETURN
      CALL timeset('qs_uspp_collocate', handle)
      CPASSERT(ASSOCIATED(rho_g%pw_grid, rho_r%pw_grid))
      CALL check_grid(uspp, rho_g%pw_grid)
      CALL augmentation_g%create(rho_g%pw_grid)
      augmentation_g%array = CMPLX(0.0_dp, 0.0_dp, dp)
      CALL collocate_density(uspp, density, rho_g%pw_grid, augmentation_g%array, charge)
      ! On even FULLSPACE grids, project the Nyquist coefficients onto a real grid function.
      ! This same real projection is implicit in the adjoint's FFT of a real potential.
      IF (rho_g%pw_grid%grid_span == FULLSPACE .AND. ANY(MODULO(rho_g%pw_grid%npts, 2) == 0)) THEN
         CALL augmentation_r%create(rho_g%pw_grid)
         CALL pw_transfer(augmentation_g, augmentation_r)
         CALL pw_transfer(augmentation_r, augmentation_g)
         CALL augmentation_r%release()
      END IF
      IF (PRESENT(derivative)) THEN
         CPASSERT(.NOT. PRESENT(total_rho))
         CALL pw_derive(augmentation_g, derivative)
      ELSE IF (PRESENT(total_rho)) THEN
         total_rho = total_rho - charge
      END IF
      CALL pw_axpy(augmentation_g, rho_g)
      CALL pw_transfer(rho_g, rho_r)
      CALL augmentation_g%release()
      CALL timestop(handle)
   END SUBROUTINE qs_uspp_collocate_images

! **************************************************************************************************
!> \brief Integrate a Quickstep real local potential against USPP augmentation.
!>        This is the discrete adjoint of qs_uspp_collocate, including FFT normalization.
!> \param ks_env KS environment with the current projector overlaps.
!> \param potential Real potential; values include the volume element when weighted is true.
!> \param hamiltonian AO Hamiltonian for one spin, incremented in place.
!> \param weighted Whether potential already includes the real grid's volume element.
!> \param density Optional force-path density; adds the augmentation energy gradients to the KS forces.
! **************************************************************************************************
   SUBROUTINE qs_uspp_integrate_images(ks_env, potential, hamiltonian, weighted, density)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: potential
      TYPE(dbcsr_p_type), INTENT(INOUT)                  :: hamiltonian(:)
      LOGICAL, INTENT(IN)                                :: weighted
      TYPE(dbcsr_p_type), INTENT(IN), OPTIONAL           :: density(:)

      INTEGER                                            :: natom
      REAL(KIND=dp), ALLOCATABLE                         :: atom_gradient(:, :), cell_gradient(:, :)
      TYPE(qs_uspp_type), POINTER                        :: uspp
      TYPE(sap_int_type), POINTER                        :: sap_int(:)

      CALL get_ks_env(ks_env, uspp=uspp)
      IF (.NOT. ASSOCIATED(uspp)) RETURN
      IF (PRESENT(density)) THEN
         CALL get_ks_env(ks_env, natom=natom)
         ALLOCATE (atom_gradient(3, natom))
         atom_gradient = 0.0_dp
         CALL allocate_cell_gradient(ks_env, cell_gradient)
         CALL build_derivatives(ks_env, sap_int)
         CALL integrate_grid(uspp, potential, hamiltonian, weighted, density, sap_int, atom_gradient, cell_gradient)
         CALL accumulate_gradient(ks_env, atom_gradient, cell_gradient)
         CALL release_sap_int(sap_int)
      ELSE
         CALL integrate_grid(uspp, potential, hamiltonian, weighted)
      END IF
   END SUBROUTINE qs_uspp_integrate_images

! **************************************************************************************************
!> \brief Real-grid adjoint shared by the ionic and density-dependent local potentials.
!> \param uspp Current projector state.
!> \param potential Real local potential.
!> \param hamiltonian AO Hamiltonian, incremented in place.
!> \param weighted Whether the real potential includes dV.
!> \param density Optional force-path AO density.
!> \param sap_int Optional image-resolved projector derivatives.
!> \param atom_gradient Optional local atomic energy gradients, incremented without MPI reduction.
!> \param cell_gradient Optional strain energy gradient, incremented without MPI reduction.
! **************************************************************************************************
   SUBROUTINE integrate_grid(uspp, potential, hamiltonian, weighted, density, sap_int, atom_gradient, cell_gradient)
      TYPE(qs_uspp_type), INTENT(IN)                     :: uspp
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: potential
      TYPE(dbcsr_p_type), INTENT(INOUT)                  :: hamiltonian(:)
      LOGICAL, INTENT(IN)                                :: weighted
      TYPE(dbcsr_p_type), INTENT(IN), OPTIONAL           :: density(:)
      TYPE(sap_int_type), INTENT(IN), OPTIONAL           :: sap_int(:)
      REAL(KIND=dp), INTENT(INOUT), OPTIONAL             :: atom_gradient(:, :), cell_gradient(3, 3)

      INTEGER                                            :: center, handle, idir, jdir, n
      LOGICAL                                            :: found
      REAL(KIND=dp), ALLOCATABLE                         :: derivative(:, :, :), dm(:, :), &
                                                            matrix(:, :), multiplicity(:), &
                                                            strain(:, :, :, :)
      REAL(KIND=dp), POINTER                             :: block(:, :)
      TYPE(dbcsr_type)                                   :: onsite, projected
      TYPE(mp_comm_type)                                 :: group
      TYPE(pw_c1d_gs_type)                               :: potential_g

      CALL timeset('qs_uspp_integrate', handle)
      CALL check_grid(uspp, potential%pw_grid)
      CALL potential_g%create(potential%pw_grid)
      CALL pw_transfer(potential, potential_g)
      IF (weighted) CALL pw_scale(potential_g, 1.0_dp/potential%pw_grid%dvol)
      ALLOCATE (multiplicity(SIZE(potential_g%array)))
      multiplicity = 1.0_dp
      IF (potential%pw_grid%grid_span == HALFSPACE) THEN
         multiplicity = 2.0_dp
         IF (potential%pw_grid%have_g0) multiplicity(1) = 1.0_dp
      END IF
      CPASSERT(PRESENT(density) .EQV. PRESENT(sap_int))
      CPASSERT(PRESENT(density) .EQV. PRESENT(atom_gradient))
      CALL dbcsr_get_info(uspp%overlaps, group=group)
      CALL uspp_projector_create_onsite(uspp%overlaps, onsite)
      IF (PRESENT(density)) CALL project_density(uspp, density, projected)
      DO center = 1, SIZE(uspp%center_atoms)
         n = uspp%projector_sizes(center)
         ALLOCATE (matrix(n, n))
         IF (PRESENT(density)) THEN
            ALLOCATE (dm(n, n), derivative(n, n, 3))
            IF (PRESENT(cell_gradient)) ALLOCATE (strain(n, n, 3, 3))
            dm = 0.0_dp
            CALL dbcsr_get_block_p(projected, center, center, block, found)
            IF (found) dm(:, :) = block
            CALL group%sum(dm)
            CALL uspp_augmentation_integrate(uspp%augmentation(uspp%center_kinds(center)), &
                                             potential%pw_grid%g, uspp%position(:, center), potential%pw_grid%vol, &
                                             potential_g%array, multiplicity, matrix, derivative, strain)
            DO idir = 1, 3
               atom_gradient(idir, uspp%center_atoms(center)) = atom_gradient(idir, uspp%center_atoms(center)) &
                                                                + SUM(dm*derivative(:, :, idir))
            END DO
            IF (PRESENT(cell_gradient)) THEN
               DO jdir = 1, 3
                  DO idir = 1, 3
                     cell_gradient(idir, jdir) = cell_gradient(idir, jdir) + SUM(dm*strain(:, :, idir, jdir))
                  END DO
               END DO
               DEALLOCATE (strain)
            END IF
            DEALLOCATE (dm, derivative)
         ELSE
            CALL uspp_augmentation_integrate(uspp%augmentation(uspp%center_kinds(center)), &
                                             potential%pw_grid%g, uspp%position(:, center), potential%pw_grid%vol, &
                                             potential_g%array, multiplicity, matrix)
         END IF
         CALL group%sum(matrix)
         CALL dbcsr_get_block_p(onsite, center, center, block, found)
         IF (found) block(:, :) = matrix
         DEALLOCATE (matrix)
      END DO
      CALL apply_operator(uspp, onsite, hamiltonian, density, sap_int, atom_gradient, cell_gradient)
      CALL dbcsr_release(onsite)
      IF (PRESENT(density)) CALL dbcsr_release(projected)
      CALL potential_g%release()
      CALL timestop(handle)
   END SUBROUTINE integrate_grid

! **************************************************************************************************
!> \brief Add B (D0 + integral Q V_sr) B^T to the core Hamiltonian.
!>        V_sr includes every non-ghost UPF ion, also NC ions without augmentation.
!>        The Gaussian-core Coulomb potential enters later through the Hartree grid.
!>        Keeping this density-independent term in H_core supplies its energy via Tr(P H_core).
!> \param uspp Current geometry-dependent USPP state.
!> \param qs_kind_set Kinds with original UPF local potentials and native projector ordering.
!> \param particle_set All ions, not only augmentation centers.
!> \param grid Same FFT grid as the physical electron density.
!> \param hamiltonian Core Hamiltonian, incremented in place.
!> \param density Optional total AO density for the force path.
!> \param energy_weighted Optional total energy-weighted AO density.
!> \param sap_int Optional uncompressed projector derivatives.
!> \param atom_gradient Optional local atomic energy gradients, incremented without MPI reduction.
!> \param cell_gradient Optional strain energy gradient, incremented without MPI reduction.
! **************************************************************************************************
   SUBROUTINE qs_uspp_core_matrix_images(uspp, qs_kind_set, particle_set, grid, hamiltonian, &
                                         density, energy_weighted, sap_int, atom_gradient, cell_gradient)
      TYPE(qs_uspp_type), INTENT(IN)                     :: uspp
      TYPE(qs_kind_type), INTENT(INOUT)                  :: qs_kind_set(:)
      TYPE(particle_type), INTENT(IN)                    :: particle_set(:)
      TYPE(pw_grid_type), POINTER                        :: grid
      TYPE(dbcsr_p_type), INTENT(INOUT)                  :: hamiltonian(:)
      TYPE(dbcsr_p_type), INTENT(IN), OPTIONAL           :: density(:), energy_weighted(:)
      TYPE(sap_int_type), INTENT(IN), OPTIONAL           :: sap_int(:)
      REAL(KIND=dp), INTENT(INOUT), OPTIONAL             :: atom_gradient(:, :), cell_gradient(3, 3)

      COMPLEX(KIND=dp), ALLOCATABLE                      :: ionic_derivative(:, :), ionic_dg(:, :), &
                                                            ionic_g(:)
      INTEGER                                            :: atom, center, handle, idir, ikind, jdir
      LOGICAL                                            :: found
      REAL(KIND=dp)                                      :: charge
      REAL(KIND=dp), ALLOCATABLE                         :: multiplicity(:)
      REAL(KIND=dp), POINTER                             :: block(:, :)
      TYPE(dbcsr_type)                                   :: bare
      TYPE(pw_c1d_gs_type)                               :: potential_g, rho_g
      TYPE(pw_r3d_rs_type)                               :: potential_r, rho_r

      CALL timeset('qs_uspp_core_matrix', handle)
      CALL check_grid(uspp, grid)
      CPASSERT(PRESENT(density) .EQV. PRESENT(energy_weighted))
      CPASSERT(PRESENT(density) .EQV. PRESENT(sap_int))
      CPASSERT(PRESENT(density) .EQV. PRESENT(atom_gradient))
      IF (PRESENT(cell_gradient)) THEN
         CPASSERT(PRESENT(density))
      END IF
      CALL uspp_projector_create_onsite(uspp%overlaps, bare)
      DO center = 1, SIZE(uspp%center_kinds)
         ikind = uspp%center_kinds(center)
         CALL dbcsr_get_block_p(bare, center, center, block, found)
         IF (found) block(:, :) = qs_kind_set(ikind)%sgp_potential%upf_projectors%d
      END DO
      CALL apply_operator(uspp, bare, hamiltonian, density, sap_int, atom_gradient, cell_gradient)
      IF (PRESENT(density)) THEN
         DO center = 1, SIZE(uspp%center_kinds)
            ikind = uspp%center_kinds(center)
            CALL dbcsr_get_block_p(bare, center, center, block, found)
            IF (found) block(:, :) = uspp%augmentation(ikind)%q
         END DO
         CALL apply_operator(uspp, bare, density=energy_weighted, sap_int=sap_int, &
                             atom_gradient=atom_gradient, cell_gradient=cell_gradient, alpha=-1.0_dp)
      END IF
      CALL dbcsr_release(bare)

      CALL potential_g%create(grid)
      CALL potential_r%create(grid)
      ALLOCATE (ionic_g(SIZE(potential_g%array)))
      potential_g%array = CMPLX(0.0_dp, 0.0_dp, dp)
      IF (PRESENT(density)) THEN
         CALL rho_g%create(grid)
         CALL rho_r%create(grid)
         rho_g%array = CMPLX(0.0_dp, 0.0_dp, dp)
         CALL collocate_density(uspp, density, grid, rho_g%array, charge)
         ! Self-adjoint real projection is needed before contracting a raw ionic derivative.
         CALL pw_transfer(rho_g, rho_r)
         CALL pw_transfer(rho_r, rho_g)
         ALLOCATE (ionic_derivative(3, SIZE(ionic_g)), multiplicity(SIZE(ionic_g)))
         IF (PRESENT(cell_gradient)) ALLOCATE (ionic_dg(3, SIZE(ionic_g)))
         multiplicity = 1.0_dp
         IF (grid%grid_span == HALFSPACE) THEN
            multiplicity = 2.0_dp
            IF (grid%have_g0) multiplicity(1) = 1.0_dp
         END IF
      END IF
      DO atom = 1, SIZE(particle_set)
         ikind = particle_set(atom)%atomic_kind%kind_number
         IF (qs_kind_set(ikind)%ghost) CYCLE
         CPASSERT(ASSOCIATED(qs_kind_set(ikind)%sgp_potential))
         ASSOCIATE (potential => qs_kind_set(ikind)%sgp_potential)
            CPASSERT(ALLOCATED(potential%upf_data))
            IF (potential%upf_local%radial%degree < 0) THEN
               CALL upf_local_init(potential%upf_data, potential%alpha_core_charge, 0, potential%upf_local)
            END IF
            IF (PRESENT(density)) THEN
               CALL upf_local_fourier(potential%upf_local, grid%g, particle_set(atom)%r, grid%vol, &
                                      ionic_g, ionic_derivative, ionic_dg)
               DO idir = 1, 3
                  atom_gradient(idir, atom) = atom_gradient(idir, atom) + grid%vol* &
                                              SUM(multiplicity*REAL(CONJG(rho_g%array)*ionic_derivative(idir, :), KIND=dp))
               END DO
               IF (PRESENT(cell_gradient)) THEN
                  ! The ionic phase is invariant at fixed fractional coordinates.
                  ! The -V/volume term cancels the explicit volume weight in integral rho_aug V_sr.
                  DO jdir = 1, 3
                     ionic_dg(jdir, :) = ionic_dg(jdir, :) + CMPLX(0.0_dp, particle_set(atom)%r(jdir), dp)*ionic_g
                     DO idir = 1, 3
                        cell_gradient(idir, jdir) = cell_gradient(idir, jdir) - grid%vol* &
                                               SUM(multiplicity*REAL(CONJG(rho_g%array)*grid%g(idir, :)*ionic_dg(jdir, :), KIND=dp))
                     END DO
                  END DO
               END IF
            ELSE
               CALL upf_local_fourier(potential%upf_local, grid%g, particle_set(atom)%r, grid%vol, ionic_g)
            END IF
            potential_g%array = potential_g%array + ionic_g
         END ASSOCIATE
      END DO
      ! Explicit real projection also treats aliased Nyquist pairs consistently on even grids.
      CALL pw_transfer(potential_g, potential_r)
      CALL integrate_grid(uspp, potential_r, hamiltonian, weighted=.FALSE., &
                          density=density, sap_int=sap_int, atom_gradient=atom_gradient, cell_gradient=cell_gradient)
      IF (PRESENT(density)) THEN
         CALL rho_g%release()
         CALL rho_r%release()
      END IF
      CALL potential_g%release()
      CALL potential_r%release()
      CALL timestop(handle)
   END SUBROUTINE qs_uspp_core_matrix_images

! **************************************************************************************************
!> \brief Core USPP operator and its stationary energy gradients, including the metric Pulay term.
!> \param ks_env Current geometry and force accumulators.
!> \param grid Physical density FFT grid.
!> \param hamiltonian Core Hamiltonian to increment.
!> \param density Spin densities; summed without modifying the originals.
!> \param weighted_density Total energy-weighted density from both spins.
! **************************************************************************************************
   SUBROUTINE qs_uspp_core_forces(ks_env, grid, hamiltonian, density, weighted_density)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(pw_grid_type), POINTER                        :: grid
      TYPE(dbcsr_p_type), INTENT(INOUT)                  :: hamiltonian(:)
      TYPE(dbcsr_p_type), INTENT(IN)                     :: density(:, :), weighted_density(:)

      INTEGER                                            :: image, ispin
      REAL(KIND=dp), ALLOCATABLE                         :: atom_gradient(:, :), cell_gradient(:, :)
      TYPE(dbcsr_p_type), ALLOCATABLE                    :: total_density(:)
      TYPE(particle_type), POINTER                       :: particle_set(:)
      TYPE(qs_kind_type), POINTER                        :: qs_kind_set(:)
      TYPE(qs_uspp_type), POINTER                        :: uspp
      TYPE(sap_int_type), POINTER                        :: sap_int(:)

      CALL get_ks_env(ks_env, uspp=uspp, qs_kind_set=qs_kind_set, particle_set=particle_set)
      CPASSERT(ASSOCIATED(uspp))
      ALLOCATE (total_density(SIZE(density, 2)))
      DO image = 1, SIZE(density, 2)
         ALLOCATE (total_density(image)%matrix)
         CALL dbcsr_copy(total_density(image)%matrix, density(1, image)%matrix)
         DO ispin = 2, SIZE(density, 1)
            CALL dbcsr_add(total_density(image)%matrix, density(ispin, image)%matrix, 1.0_dp, 1.0_dp)
         END DO
      END DO
      ALLOCATE (atom_gradient(3, SIZE(particle_set)))
      atom_gradient = 0.0_dp
      CALL allocate_cell_gradient(ks_env, cell_gradient)
      CALL build_derivatives(ks_env, sap_int)
      CALL qs_uspp_core_matrix(uspp, qs_kind_set, particle_set, grid, hamiltonian, &
                               total_density, weighted_density, sap_int, atom_gradient, cell_gradient)
      CALL accumulate_gradient(ks_env, atom_gradient, cell_gradient)
      CALL release_sap_int(sap_int)
      DO image = 1, SIZE(total_density)
         CALL dbcsr_release(total_density(image)%matrix)
         DEALLOCATE (total_density(image)%matrix)
      END DO
   END SUBROUTINE qs_uspp_core_forces

! **************************************************************************************************
!> \brief Build image-resolved derivatives only on the force path.
!> \param ks_env Current ABBA neighbor list and native radial projector data.
!> \param sap_int Newly allocated uncompressed projector overlaps and derivatives.
! **************************************************************************************************
   SUBROUTINE build_derivatives(ks_env, sap_int)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(sap_int_type), POINTER                        :: sap_int(:)

      TYPE(neighbor_list_set_p_type), POINTER            :: sap_ppnl(:)
      TYPE(qs_kind_type), POINTER                        :: qs_kind_set(:)

      CALL get_ks_env(ks_env, qs_kind_set=qs_kind_set, sap_ppnl=sap_ppnl)
      ALLOCATE (sap_int(SIZE(qs_kind_set)**2))
      CALL build_sap_ints(sap_int, sap_ppnl, qs_kind_set, nder=1)
      CALL sap_sort(sap_int)
   END SUBROUTINE build_derivatives

! **************************************************************************************************
!> \brief Add unreduced USPP energy gradients to the normal Quickstep force accumulator.
!> \param ks_env Kinds and force arrays; Quickstep later reduces and negates these gradients.
!> \param atom_gradient Local energy gradients in global particle order.
!> \param cell_gradient Optional strain energy gradient; its negative is the virial.
! **************************************************************************************************
   SUBROUTINE accumulate_gradient(ks_env, atom_gradient, cell_gradient)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      REAL(KIND=dp), INTENT(IN)                          :: atom_gradient(:, :)
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: cell_gradient(3, 3)

      INTEGER                                            :: atom, ia, ikind
      INTEGER, ALLOCATABLE                               :: atom_of_kind(:), kind_of(:)
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind_set(:)
      TYPE(qs_force_type), POINTER                       :: force(:)
      TYPE(virial_type), POINTER                         :: virial

      CALL get_ks_env(ks_env, atomic_kind_set=atomic_kind_set, force=force)
      CALL get_atomic_kind_set(atomic_kind_set, atom_of_kind=atom_of_kind, kind_of=kind_of)
      DO atom = 1, SIZE(atom_gradient, 2)
         ikind = kind_of(atom)
         ia = atom_of_kind(atom)
         force(ikind)%gth_ppnl(:, ia) = force(ikind)%gth_ppnl(:, ia) + atom_gradient(:, atom)
      END DO
      IF (PRESENT(cell_gradient)) THEN
         CALL get_ks_env(ks_env, virial=virial)
         virial%pv_ppnl = virial%pv_ppnl - cell_gradient
         virial%pv_virial = virial%pv_virial - cell_gradient
      END IF
   END SUBROUTINE accumulate_gradient

! **************************************************************************************************
!> \brief Allocate the optional local strain accumulator only for analytic stress.
!> \param ks_env Current virial settings.
!> \param cell_gradient Allocated and zeroed for analytic stress, otherwise absent.
! **************************************************************************************************
   SUBROUTINE allocate_cell_gradient(ks_env, cell_gradient)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      REAL(KIND=dp), ALLOCATABLE, INTENT(OUT)            :: cell_gradient(:, :)

      TYPE(virial_type), POINTER                         :: virial

      CALL get_ks_env(ks_env, virial=virial)
      IF (virial%pv_availability .AND. .NOT. virial%pv_numer) THEN
         ALLOCATE (cell_gradient(3, 3))
         cell_gradient = 0.0_dp
      END IF
   END SUBROUTINE allocate_cell_gradient

! **************************************************************************************************
!> \brief Project an AO density with the Gamma matrix or the complete image convolution.
!> \param uspp Current geometry state.
!> \param density Spin density matrices over real-space images.
!> \param projected New distributed on-site density.
! **************************************************************************************************
   SUBROUTINE project_density(uspp, density, projected)
      TYPE(qs_uspp_type), INTENT(IN)                     :: uspp
      TYPE(dbcsr_p_type), INTENT(IN)                     :: density(:)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: projected

      IF (ALLOCATED(uspp%image_pairs)) THEN
         CALL qs_uspp_projector_images(uspp%overlaps, uspp%image_integrals, uspp%image_pairs, SIZE(uspp%augmentation), &
                                       uspp%center_atoms, uspp%center_kinds, density=density, projected=projected)
      ELSE
         CPASSERT(SIZE(density) == 1)
         CALL uspp_projector_density(uspp%overlaps, density(1)%matrix, projected)
      END IF
   END SUBROUTINE project_density

! **************************************************************************************************
!> \brief Collocate the projected density and its monopole using every real-space image.
!> \param uspp Current geometry state.
!> \param density Spin AO density matrices.
!> \param grid Physical density grid.
!> \param rho Reciprocal density, incremented in place.
!> \param charge Replicated total augmentation charge.
! **************************************************************************************************
   SUBROUTINE collocate_density(uspp, density, grid, rho, charge)
      TYPE(qs_uspp_type), INTENT(IN)                     :: uspp
      TYPE(dbcsr_p_type), INTENT(IN)                     :: density(:)
      TYPE(pw_grid_type), INTENT(IN)                     :: grid
      COMPLEX(KIND=dp), INTENT(INOUT)                    :: rho(:)
      REAL(KIND=dp), INTENT(OUT)                         :: charge

      INTEGER                                            :: center, ikind, n
      LOGICAL                                            :: found
      REAL(KIND=dp), ALLOCATABLE                         :: matrix(:, :)
      REAL(KIND=dp), POINTER                             :: block(:, :)
      TYPE(dbcsr_type)                                   :: projected
      TYPE(mp_comm_type)                                 :: group

      CALL dbcsr_get_info(uspp%overlaps, group=group)
      CALL project_density(uspp, density, projected)
      charge = 0.0_dp
      DO center = 1, SIZE(uspp%center_atoms)
         n = uspp%projector_sizes(center)
         ikind = uspp%center_kinds(center)
         ALLOCATE (matrix(n, n))
         matrix = 0.0_dp
         CALL dbcsr_get_block_p(projected, center, center, block, found)
         IF (found) matrix(:, :) = block
         CALL group%sum(matrix)
         CALL uspp_augmentation_collocate(uspp%augmentation(ikind), grid%g, uspp%position(:, center), grid%vol, matrix, rho)
         charge = charge + SUM(matrix*uspp%augmentation(ikind)%q)
         DEALLOCATE (matrix)
      END DO
      CALL dbcsr_release(projected)
   END SUBROUTINE collocate_density

! **************************************************************************************************
!> \brief Apply a real on-site operator and optionally its projector energy derivatives.
!> \param uspp Current geometry state.
!> \param onsite Distributed on-site operator.
!> \param hamiltonian Optional AO image matrices to increment.
!> \param density Optional AO image densities for derivatives.
!> \param sap_int Optional sorted image-resolved derivatives.
!> \param atom_gradient Optional local atomic energy gradients.
!> \param cell_gradient Optional local strain energy gradient.
!> \param alpha Optional operator/gradient prefactor, e.g. -1 for metric derivatives.
! **************************************************************************************************
   SUBROUTINE apply_operator(uspp, onsite, hamiltonian, density, sap_int, atom_gradient, cell_gradient, alpha)
      TYPE(qs_uspp_type), INTENT(IN)                     :: uspp
      TYPE(dbcsr_type), INTENT(IN)                       :: onsite
      TYPE(dbcsr_p_type), INTENT(INOUT), OPTIONAL        :: hamiltonian(:)
      TYPE(dbcsr_p_type), INTENT(IN), OPTIONAL           :: density(:)
      TYPE(sap_int_type), INTENT(IN), OPTIONAL           :: sap_int(:)
      REAL(KIND=dp), INTENT(INOUT), OPTIONAL             :: atom_gradient(:, :), cell_gradient(3, 3)
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: alpha

      REAL(KIND=dp)                                      :: scale
      TYPE(dbcsr_type)                                   :: gradient, scaled

      scale = 1.0_dp
      IF (PRESENT(alpha)) scale = alpha
      CPASSERT(PRESENT(density) .EQV. PRESENT(sap_int))
      CPASSERT(PRESENT(density) .EQV. PRESENT(atom_gradient))
      IF (ALLOCATED(uspp%image_pairs)) THEN
         IF (PRESENT(density)) THEN
            CALL qs_uspp_projector_images(uspp%overlaps, sap_int, uspp%image_pairs, SIZE(uspp%augmentation), &
                                          uspp%center_atoms, uspp%center_kinds, onsite, hamiltonian, density, &
                                          atom_gradient=atom_gradient, cell_gradient=cell_gradient, alpha=scale)
         ELSE
            CALL qs_uspp_projector_images(uspp%overlaps, uspp%image_integrals, uspp%image_pairs, SIZE(uspp%augmentation), &
                                          uspp%center_atoms, uspp%center_kinds, onsite, hamiltonian, alpha=scale)
         END IF
      ELSE
         IF (PRESENT(hamiltonian)) THEN
            CPASSERT(SIZE(hamiltonian) == 1)
            CALL dbcsr_copy(scaled, onsite)
            CALL dbcsr_scale(scaled, scale)
            CALL uspp_projector_add_operator(uspp%overlaps, scaled, hamiltonian(1)%matrix)
            CALL dbcsr_release(scaled)
         END IF
         IF (PRESENT(density)) THEN
            CPASSERT(SIZE(density) == 1)
            CALL uspp_projector_pullback(uspp%overlaps, density(1)%matrix, onsite, scale, gradient)
            CALL qs_uspp_projector_gradient(gradient, sap_int, uspp%center_atoms, atom_gradient, cell_gradient)
            CALL dbcsr_release(gradient)
         END IF
      END IF
   END SUBROUTINE apply_operator

! **************************************************************************************************
!> \brief Scalar adapter to the image-aware density path.
!> \param ks_env Current KS environment.
!> \param density Gamma AO density.
!> \param rho_g Reciprocal density.
!> \param rho_r Real density.
!> \param total_rho Optional electronic charge.
!> \param derivative Optional derivative orders.
! **************************************************************************************************
   SUBROUTINE qs_uspp_collocate_single(ks_env, density, rho_g, rho_r, total_rho, derivative)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(dbcsr_type), INTENT(IN), TARGET               :: density
      TYPE(pw_c1d_gs_type), INTENT(INOUT)                :: rho_g
      TYPE(pw_r3d_rs_type), INTENT(INOUT)                :: rho_r
      REAL(KIND=dp), INTENT(INOUT), OPTIONAL             :: total_rho
      INTEGER, INTENT(IN), OPTIONAL                      :: derivative(3)

      TYPE(dbcsr_p_type)                                 :: matrices(1)

      matrices(1)%matrix => density
      CALL qs_uspp_collocate_images(ks_env, matrices, rho_g, rho_r, total_rho, derivative)
   END SUBROUTINE qs_uspp_collocate_single

! **************************************************************************************************
!> \brief Scalar adapter to the image-aware grid adjoint.
!> \param ks_env Current KS environment.
!> \param potential Real grid potential.
!> \param hamiltonian Gamma Hamiltonian to increment.
!> \param weighted Whether potential includes dV.
!> \param density Optional Gamma force-path density.
! **************************************************************************************************
   SUBROUTINE qs_uspp_integrate_single(ks_env, potential, hamiltonian, weighted, density)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      TYPE(pw_r3d_rs_type), INTENT(IN)                   :: potential
      TYPE(dbcsr_type), INTENT(INOUT), TARGET            :: hamiltonian
      LOGICAL, INTENT(IN)                                :: weighted
      TYPE(dbcsr_type), INTENT(IN), OPTIONAL, TARGET     :: density

      TYPE(dbcsr_p_type)                                 :: matrices(1)
      TYPE(dbcsr_p_type), ALLOCATABLE                    :: densities(:)

      matrices(1)%matrix => hamiltonian
      IF (PRESENT(density)) THEN
         ALLOCATE (densities(1))
         densities(1)%matrix => density
      END IF
      CALL qs_uspp_integrate_images(ks_env, potential, matrices, weighted, densities)
   END SUBROUTINE qs_uspp_integrate_single

! **************************************************************************************************
!> \brief Scalar adapter to the image-aware ionic operator and metric derivatives.
!> \param uspp Current USPP state.
!> \param qs_kind_set Current kinds.
!> \param particle_set Current atoms.
!> \param grid Physical density grid.
!> \param hamiltonian Gamma core Hamiltonian.
!> \param density Optional total Gamma density.
!> \param energy_weighted Optional total Gamma energy-weighted density.
!> \param sap_int Optional projector derivatives.
!> \param atom_gradient Optional local atomic energy gradients.
!> \param cell_gradient Optional local strain energy gradient.
! **************************************************************************************************
   SUBROUTINE qs_uspp_core_matrix_single(uspp, qs_kind_set, particle_set, grid, hamiltonian, &
                                         density, energy_weighted, sap_int, atom_gradient, cell_gradient)
      TYPE(qs_uspp_type), INTENT(IN)                     :: uspp
      TYPE(qs_kind_type), INTENT(INOUT)                  :: qs_kind_set(:)
      TYPE(particle_type), INTENT(IN)                    :: particle_set(:)
      TYPE(pw_grid_type), POINTER                        :: grid
      TYPE(dbcsr_type), INTENT(INOUT), TARGET            :: hamiltonian
      TYPE(dbcsr_type), INTENT(IN), OPTIONAL, TARGET     :: density, energy_weighted
      TYPE(sap_int_type), INTENT(IN), OPTIONAL           :: sap_int(:)
      REAL(KIND=dp), INTENT(INOUT), OPTIONAL             :: atom_gradient(:, :), cell_gradient(3, 3)

      TYPE(dbcsr_p_type)                                 :: matrices(1)
      TYPE(dbcsr_p_type), ALLOCATABLE                    :: densities(:), weighted(:)

      matrices(1)%matrix => hamiltonian
      IF (PRESENT(density)) THEN
         ALLOCATE (densities(1))
         densities(1)%matrix => density
      END IF
      IF (PRESENT(energy_weighted)) THEN
         ALLOCATE (weighted(1))
         weighted(1)%matrix => energy_weighted
      END IF
      CALL qs_uspp_core_matrix_images(uspp, qs_kind_set, particle_set, grid, matrices, &
                                      densities, weighted, sap_int, atom_gradient, cell_gradient)
   END SUBROUTINE qs_uspp_core_matrix_single

! **************************************************************************************************
!> \brief Require every retained G to occur once over the projector communicator.
!> \param uspp Current projector state.
!> \param grid Quickstep FFT grid; a reordered communicator with the same members is allowed.
! **************************************************************************************************
   SUBROUTINE check_grid(uspp, grid)
      TYPE(qs_uspp_type), INTENT(IN)                     :: uspp
      TYPE(pw_grid_type), INTENT(IN)                     :: grid

      TYPE(mp_comm_type)                                 :: group

      CALL dbcsr_get_info(uspp%overlaps, group=group)
      CPASSERT(group%compare(grid%para%group) /= mp_comm_unequal)
      CPASSERT(grid%grid_span == FULLSPACE .OR. grid%grid_span == HALFSPACE)
      IF (grid%grid_span == HALFSPACE) THEN
         CPASSERT(ALL(MODULO(grid%npts, 2) == 1))
      END IF
   END SUBROUTINE check_grid

END MODULE qs_uspp_methods
