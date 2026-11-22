!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Inverse USPP point fields in the actual Quickstep atom/AO and magnetic-projector order.
!>        Point rows are distributed; full AO/projector columns are replicated in this adapter.
! **************************************************************************************************
MODULE qs_upf_inverse
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cryssym,                         ONLY: crys_sym_lattice_rotations
   USE gaussian_basis_fields,           ONLY: gaussian_basis_fields_periodic
   USE input_constants,                 ONLY: uspp_method_inverse
   USE kinds,                           ONLY: dp
   USE lebedev,                         ONLY: deallocate_lebedev_grids,&
                                              get_number_of_lebedev_grid,&
                                              init_lebedev_grids,&
                                              lebedev_average_points,&
                                              lebedev_grid
   USE mathconstants,                   ONLY: twopi
   USE message_passing,                 ONLY: mp_comm_type
   USE orbital_pointers,                ONLY: coset,&
                                              ncoset
   USE particle_types,                  ONLY: particle_type
   USE qs_grid_atom,                    ONLY: allocate_grid_atom,&
                                              create_grid_atom,&
                                              grid_atom_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_atom_grid_ref_type,&
                                              qs_upf_inverse_type
   USE qs_uspp_types,                   ONLY: qs_uspp_type
   USE upf_inverse,                     ONLY: &
        upf_inverse_atomic_difference_init, upf_inverse_atomic_init, upf_inverse_atomic_residual, &
        upf_inverse_basis, upf_inverse_basis_response, upf_inverse_fields_init, upf_inverse_init, &
        upf_inverse_operator_response, upf_inverse_project, upf_inverse_project_response
   USE upf_projector_integrals,         ONLY: upf_projector_init,&
                                              upf_projector_overlap,&
                                              upf_projector_type
   USE upf_radial_fields,               ONLY: upf_radial_fields_evaluate,&
                                              upf_radial_fields_periodic
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: qs_upf_inverse_type, qs_upf_inverse_prepare, qs_upf_inverse_from_ks
   PUBLIC :: qs_upf_inverse_response, qs_upf_inverse_release
   PUBLIC :: qs_upf_inverse_atomic_overlap
   PUBLIC :: qs_upf_inverse_atomic_operators

CONTAINS

! **************************************************************************************************
!> \brief Contract atom-local hard projectors with native Gaussian AOs and explicit Bloch images.
!>        Only the original radial quadrature is numerical. Gaussian Cartesian integrals and
!>        AO-center derivatives use the existing UPF integrator, without a uniform-grid solve.
!> \param basis Native contracted ORB basis and unchanged primitive normalization.
!> \param hard Atom-local transformed projectors, including enough derivative degree.
!> \param center Home projector center.
!> \param ao_center Home Gaussian center.
!> \param cell Cartesian lattice vectors.
!> \param images Fixed integer Gaussian images.
!> \param reduced_k Bloch vector, with phase exp(i*2*pi*k.image).
!> \param overlaps Hard projector bras against Gaussian AO columns, without occupations.
!> \param response Optional joint first response.
!> \param center_response Optional projector-center displacement.
!> \param ao_response Optional Gaussian-center displacement.
!> \param cell_response Optional lattice displacement.
!> \param reduced_k_response Optional reduced Bloch-vector displacement.
!> \param kinetic Optional projector overlaps with the Gaussian kinetic operator.
!> \param kinetic_response Optional joint first response of the kinetic overlaps.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_atomic_overlap(basis, hard, center, ao_center, cell, images, reduced_k, overlaps, &
                                            response, center_response, ao_response, cell_response, reduced_k_response, &
                                            kinetic, kinetic_response)
      TYPE(gto_basis_set_type), INTENT(IN)               :: basis
      TYPE(upf_projector_type), INTENT(IN)               :: hard
      REAL(dp), INTENT(IN)                               :: center(3), ao_center(3), cell(3, 3)
      INTEGER, INTENT(IN)                                :: images(:, :)
      REAL(dp), INTENT(IN)                               :: reduced_k(3)
      COMPLEX(dp), INTENT(OUT)                           :: overlaps(:, :)
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: response(:, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: center_response(3), ao_response(3), &
                                                            cell_response(3, 3), &
                                                            reduced_k_response(3)
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: kinetic(:, :), kinetic_response(:, :)

      COMPLEX(dp)                                        :: phase, phase_response
      INTEGER                                            :: axis, derivative_index, first, image, &
                                                            iset, last, nder, np, nshell, powers(3)
      REAL(dp)                                           :: direction(3), displacement(3), dk(3), &
                                                            shift(3)
      REAL(dp), ALLOCATABLE                              :: contracted(:, :), derivative(:, :), &
                                                            integrals(:, :, :)

      CPASSERT(SIZE(images, 1) == 3)
      CPASSERT(ALL(SHAPE(overlaps) == [hard%nproj, basis%nsgf]))
      nder = 0
      IF (PRESENT(response)) THEN
         CPASSERT(ALL(SHAPE(response) == SHAPE(overlaps)))
         response = CMPLX(0.0_dp, 0.0_dp, dp)
         nder = 1
      END IF
      IF (PRESENT(kinetic)) THEN
         CPASSERT(ALL(SHAPE(kinetic) == SHAPE(overlaps)))
         kinetic = CMPLX(0.0_dp, 0.0_dp, dp)
         nder = MAX(nder, 2)
      END IF
      IF (PRESENT(kinetic_response)) THEN
         CPASSERT(PRESENT(kinetic))
         CPASSERT(ALL(SHAPE(kinetic_response) == SHAPE(overlaps)))
         kinetic_response = CMPLX(0.0_dp, 0.0_dp, dp)
         nder = 3
      END IF
      overlaps = CMPLX(0.0_dp, 0.0_dp, dp)
      dk = 0.0_dp
      IF (PRESENT(reduced_k_response)) dk = reduced_k_response
      DO image = 1, SIZE(images, 2)
         shift = MATMUL(cell, REAL(images(:, image), dp))
         displacement = center - ao_center - shift
         direction = 0.0_dp
         IF (PRESENT(ao_response)) direction = direction + ao_response
         IF (PRESENT(center_response)) direction = direction - center_response
         IF (PRESENT(cell_response)) direction = direction + MATMUL(cell_response, REAL(images(:, image), dp))
         phase = EXP(CMPLX(0.0_dp, twopi*DOT_PRODUCT(reduced_k, REAL(images(:, image), dp)), dp))
         phase_response = CMPLX(0.0_dp, twopi*DOT_PRODUCT(dk, REAL(images(:, image), dp)), dp)*phase
         DO iset = 1, basis%nset
            ! Twice the native Gaussian radius retains a margin beyond the native tail threshold.
            IF (NORM2(displacement) > 2.0_dp*basis%set_radius(iset) + hard%radius) CYCLE
            nshell = basis%nshell(iset)
            first = MINVAL(basis%first_sgf(:nshell, iset))
            last = MAXVAL(basis%last_sgf(:nshell, iset))
            np = basis%npgf(iset)*ncoset(basis%lmax(iset))
            ALLOCATE (integrals(np, hard%nproj, ncoset(nder)), contracted(last - first + 1, hard%nproj))
            CALL upf_projector_overlap(hard, basis%lmax(iset), basis%zet(:basis%npgf(iset), iset), &
                                       displacement, nder, integrals)
            contracted(:, :) = MATMUL(TRANSPOSE(basis%sphi(:np, first:last)), integrals(:, :, 1))
            overlaps(:, first:last) = overlaps(:, first:last) + phase*TRANSPOSE(contracted)
            IF (PRESENT(response)) THEN
               response(:, first:last) = response(:, first:last) + phase_response*TRANSPOSE(contracted)
               ALLOCATE (derivative(last - first + 1, hard%nproj))
               derivative = 0.0_dp
               DO axis = 1, 3
                  derivative(:, :) = derivative + direction(axis)* &
                                     MATMUL(TRANSPOSE(basis%sphi(:np, first:last)), integrals(:, :, axis + 1))
               END DO
               response(:, first:last) = response(:, first:last) + phase*TRANSPOSE(derivative)
               DEALLOCATE (derivative)
            END IF
            IF (PRESENT(kinetic)) THEN
               contracted = 0.0_dp
               DO axis = 1, 3
                  powers = 0
                  powers(axis) = 2
                  derivative_index = coset(powers(1), powers(2), powers(3))
                  contracted(:, :) = contracted - 0.5_dp* &
                                     MATMUL(TRANSPOSE(basis%sphi(:np, first:last)), integrals(:, :, derivative_index))
               END DO
               kinetic(:, first:last) = kinetic(:, first:last) + phase*TRANSPOSE(contracted)
               IF (PRESENT(kinetic_response)) THEN
                  kinetic_response(:, first:last) = kinetic_response(:, first:last) + phase_response*TRANSPOSE(contracted)
                  ALLOCATE (derivative(last - first + 1, hard%nproj))
                  derivative = 0.0_dp
                  DO axis = 1, 3
                     DO derivative_index = 1, 3
                        powers = 0
                        powers(axis) = 2
                        powers(derivative_index) = powers(derivative_index) + 1
                        derivative(:, :) = derivative - 0.5_dp*direction(derivative_index)* &
                                           MATMUL(TRANSPOSE(basis%sphi(:np, first:last)), &
                                                  integrals(:, :, coset(powers(1), powers(2), powers(3))))
                     END DO
                  END DO
                  kinetic_response(:, first:last) = kinetic_response(:, first:last) + phase*TRANSPOSE(derivative)
                  DEALLOCATE (derivative)
               END IF
            END IF
            DEALLOCATE (integrals, contracted)
         END DO
      END DO
   END SUBROUTINE qs_upf_inverse_atomic_overlap

! **************************************************************************************************
!> \brief Read native geometry, ORB bases and original USPP operators from a prepared KS environment.
!>        Does not select an input METHOD or replace the reference SCF's soft-orbital state.
!> \param ks_env Live Quickstep KS environment with prepared USPP state.
!> \param points Local Cartesian quadrature points, indexed (axis,point).
!> \param weights Matching local measure; each point occurs once in comm.
!> \param images Explicit image set, fixed during a geometry response.
!> \param reduced_k Reduced Bloch vector for this frame.
!> \param comm Communicator distributing point rows, not a replicated AO process dimension.
!> \param frame New owned fields; native basis pointers remain borrowed from ks_env.
!> \param info Zero on success, otherwise the inverse-factorization status.
!> \param atom_local Optional nonoverlapping one-center inverse; the general reference is the default.
!> \param direct Use original projectors and soft orbitals without reconstruction.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_from_ks(ks_env, points, weights, images, reduced_k, comm, frame, info, atom_local, direct)
      TYPE(qs_ks_env_type), INTENT(IN)                   :: ks_env
      REAL(dp), INTENT(IN)                               :: points(:, :), weights(:)
      INTEGER, INTENT(IN)                                :: images(:, :)
      REAL(dp), INTENT(IN)                               :: reduced_k(3)
      TYPE(mp_comm_type), INTENT(IN)                     :: comm
      TYPE(qs_upf_inverse_type), INTENT(OUT)             :: frame
      INTEGER, INTENT(OUT)                               :: info
      LOGICAL, INTENT(IN), OPTIONAL                      :: atom_local, direct

      REAL(dp), ALLOCATABLE                              :: rotations(:, :, :)
      TYPE(cell_type), POINTER                           :: cell, cell_ref
      TYPE(dft_control_type), POINTER                    :: control
      TYPE(particle_type), POINTER                       :: particles(:)
      TYPE(qs_kind_type), POINTER                        :: kinds(:)
      TYPE(qs_uspp_type), POINTER                        :: uspp

      CALL get_ks_env(ks_env, qs_kind_set=kinds, particle_set=particles, cell=cell, cell_ref=cell_ref, &
                      dft_control=control, uspp=uspp)
      CPASSERT(ASSOCIATED(kinds) .AND. ASSOCIATED(particles) .AND. ASSOCIATED(cell))
      CPASSERT(ASSOCIATED(uspp))

      IF (ALL(control%qs_control%upf_angular_reference_cell == 0.0_dp)) THEN
         control%qs_control%upf_angular_reference_cell = cell%hmat
         IF (ASSOCIATED(cell_ref)) control%qs_control%upf_angular_reference_cell = cell_ref%hmat
      END IF
      IF (control%qs_control%uspp_method == uspp_method_inverse) THEN
         CALL crys_sym_lattice_rotations(control%qs_control%upf_angular_reference_cell, rotations)
      END IF
      CALL qs_upf_inverse_prepare(kinds, particles, cell, uspp, points, weights, images, reduced_k, comm, &
                                  frame, info, atom_local, direct, control%qs_control%gapw_control%quadrature, rotations)
   END SUBROUTINE qs_upf_inverse_from_ks

! **************************************************************************************************
!> \brief Build the full inverse from native Quickstep basis columns and original radial fields.
!>        Ghost atoms retain their Gaussian columns. NC projectors have zero reconstruction
!>        and act on the same inverse soft fields, including coupling to USPP centers.
!>        Missing original full-wave pairs are rejected by the radial importer.
!> \param kinds Native Quickstep kinds; ORB bases and potential data are not modified.
!> \param particles Native atom order, also used for AO matrix blocks.
!> \param cell Native cell and periodicity flags.
!> \param uspp Prepared center selection and native magnetic augmentation operators.
!> \param points ...
!> \param weights ...
!> \param images ...
!> \param reduced_k ...
!> \param comm ...
!> \param frame ...
!> \param info ...
!> \param atom_local Optional nonoverlapping one-center inverse. Status -2 reports overlapping support.
!> \param direct Use original projectors and soft orbitals without reconstruction.
!> \param quadrature Radial quadrature rule for the atom-centered integrals.
!> \param rotations Optional reference-cell rotations used to average the angular quadrature.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_prepare(kinds, particles, cell, uspp, points, weights, images, reduced_k, comm, &
                                     frame, info, atom_local, direct, quadrature, rotations)
      TYPE(qs_kind_type), INTENT(INOUT)                  :: kinds(:)
      TYPE(particle_type), INTENT(IN)                    :: particles(:)
      TYPE(cell_type), INTENT(IN)                        :: cell
      TYPE(qs_uspp_type), INTENT(IN)                     :: uspp
      REAL(dp), INTENT(IN)                               :: points(:, :), weights(:)
      INTEGER, INTENT(IN)                                :: images(:, :)
      REAL(dp), INTENT(IN)                               :: reduced_k(3)
      TYPE(mp_comm_type), INTENT(IN)                     :: comm
      TYPE(qs_upf_inverse_type), INTENT(OUT)             :: frame
      INTEGER, INTENT(OUT)                               :: info
      LOGICAL, INTENT(IN), OPTIONAL                      :: atom_local, direct
      INTEGER, INTENT(IN)                                :: quadrature
      REAL(dp), INTENT(IN), OPTIONAL                     :: rotations(:, :, :)

      COMPLEX(dp), ALLOCATABLE                           :: coupling(:, :), overlaps(:, :)
      INTEGER                                            :: atom, axis, center, degree, first, &
                                                            ikind, last, nao, ncenter, npoint, &
                                                            nproj, p, q
      LOGICAL                                            :: use_atomic
      TYPE(gto_basis_set_type), POINTER                  :: basis

      CPASSERT(SIZE(points, 1) == 3 .AND. SIZE(points, 2) == SIZE(weights))
      CPASSERT(SIZE(images, 1) == 3 .AND. SIZE(images, 2) > 0)
      use_atomic = .FALSE.
      IF (PRESENT(atom_local)) use_atomic = atom_local
      IF (PRESENT(direct)) frame%direct = direct
      IF (frame%direct) use_atomic = .TRUE.
      info = 0
      DO axis = 1, 3
         IF (cell%perd(axis) == 0) THEN
            CPASSERT(ALL(images(axis, :) == 0))
            CPASSERT(reduced_k(axis) == 0.0_dp)
         END IF
      END DO
      CPASSERT(SIZE(uspp%center_atoms) > 0)
      npoint = SIZE(weights)
      frame%point_comm = comm
      frame%cell = cell%hmat
      frame%reduced_k = reduced_k
      ALLOCATE (frame%points, SOURCE=points)
      ALLOCATE (frame%weights, SOURCE=weights)
      ALLOCATE (frame%images, SOURCE=images)
      ALLOCATE (frame%positions(3, SIZE(particles)), frame%bases(SIZE(particles)), &
                frame%ao_first(SIZE(particles)), frame%ao_last(SIZE(particles)))
      nao = 0
      degree = 0
      DO atom = 1, SIZE(particles)
         CPASSERT(ASSOCIATED(particles(atom)%atomic_kind))
         ikind = particles(atom)%atomic_kind%kind_number
         CPASSERT(ikind >= 1 .AND. ikind <= SIZE(kinds))
         CALL get_qs_kind(kinds(ikind), basis_set=basis, basis_type='ORB')
         CPASSERT(ASSOCIATED(basis))
         degree = MAX(degree, MAXVAL(basis%lmax) + 3)
         frame%bases(atom)%gto_basis_set => basis
         frame%ao_first(atom) = nao + 1
         nao = nao + basis%nsgf
         frame%ao_last(atom) = nao
         frame%positions(:, atom) = particles(atom)%r
      END DO
      ALLOCATE (frame%gaussian(npoint, nao), frame%gaussian_gradient(npoint, nao, 3), &
                frame%soft(npoint, nao), frame%soft_gradient(npoint, nao, 3))
      DO atom = 1, SIZE(particles)
         first = frame%ao_first(atom)
         last = frame%ao_last(atom)
         CALL gaussian_basis_fields_periodic(frame%bases(atom)%gto_basis_set, points, frame%positions(:, atom), &
                          frame%cell, images, reduced_k, frame%gaussian(:, first:last), frame%gaussian_gradient(:, first:last, :), &
                                             cutoff_squared=4.0_dp*frame%bases(atom)%gto_basis_set%kind_radius**2)
      END DO
      ncenter = SIZE(uspp%center_atoms)
      ALLOCATE (frame%center_atoms, SOURCE=uspp%center_atoms)
      ALLOCATE (frame%center_kinds, SOURCE=uspp%center_kinds)
      ALLOCATE (frame%projector_first(ncenter), frame%projector_last(ncenter), &
                frame%beta_cache(SIZE(kinds)), frame%delta_cache(SIZE(kinds)))
      ALLOCATE (frame%atomic_cache(SIZE(kinds)))
      IF (.NOT. use_atomic) THEN
         ALLOCATE (frame%defining_atom_grids(SIZE(kinds)))
         CALL init_lebedev_grids()
      END IF
      nproj = 0
      DO center = 1, ncenter
         ikind = frame%center_kinds(center)
         atom = frame%center_atoms(center)
         CPASSERT(.NOT. kinds(ikind)%ghost)
         CPASSERT(particles(atom)%atomic_kind%kind_number == ikind)
         IF (frame%beta_cache(ikind)%nproj == 0) THEN
            ASSOCIATE (pot => kinds(ikind)%sgp_potential%upf_data)
               IF (use_atomic) THEN
                  CALL upf_inverse_fields_init(pot, frame%beta_cache(ikind), frame%delta_cache(ikind), frame%direct)
                  CALL upf_inverse_atomic_init(pot, degree, frame%atomic_cache(ikind), info, frame%direct)
                  IF (info /= 0) RETURN
               ELSE
                  ASSOCIATE (rule => frame%defining_atom_grids(ikind))
                     CALL qs_upf_atom_kind_grid(kinds(ikind), quadrature, rule%grid)
                     CALL upf_projector_init(pot, degree, rule%projectors)
                     CALL upf_inverse_fields_init(pot, frame%beta_cache(ikind), frame%delta_cache(ikind), &
                                                  home_coupling=rule%home_coupling)
                     IF (PRESENT(rotations)) THEN
                        CALL qs_upf_atom_angular_average(rule%grid, rotations, rule%directions, rule%angular_weights)
                     END IF
                  END ASSOCIATE
                  CALL upf_inverse_atomic_difference_init(pot, degree, frame%atomic_cache(ikind))
               END IF
            END ASSOCIATE
         END IF
         p = kinds(ikind)%sgp_potential%upf_projectors%nproj
         CPASSERT(frame%beta_cache(ikind)%nproj == p)
         frame%projector_first(center) = nproj + 1
         nproj = nproj + frame%beta_cache(ikind)%nproj
         frame%projector_last(center) = nproj
      END DO
      IF (.NOT. use_atomic) CALL deallocate_lebedev_grids()
      ALLOCATE (frame%beta(npoint, nproj), frame%delta(npoint, nproj), frame%delta_gradient(npoint, nproj, 3), &
                frame%q(nproj, nproj), frame%radial_q(nproj, nproj), frame%dion(nproj, nproj), frame%coefficients(nproj, nao))
      frame%q(:, :) = CMPLX(0.0_dp, 0.0_dp, dp)
      frame%radial_q(:, :) = CMPLX(0.0_dp, 0.0_dp, dp)
      frame%dion(:, :) = CMPLX(0.0_dp, 0.0_dp, dp)
      DO center = 1, ncenter
         ikind = frame%center_kinds(center)
         atom = frame%center_atoms(center)
         first = frame%projector_first(center)
         last = frame%projector_last(center)
         frame%radial_q(first:last, first:last) = uspp%augmentation(ikind)%q
         ASSOCIATE (projectors => kinds(ikind)%sgp_potential%upf_projectors, &
                    pot => kinds(ikind)%sgp_potential%upf_data)
            DO q = 1, last - first + 1
               DO p = 1, last - first + 1
                  IF (projectors%angular_index(p) /= projectors%angular_index(q)) CYCLE
                  frame%q(first + p - 1, first + q - 1) = pot%qij(projectors%radial_index(p), projectors%radial_index(q))
               END DO
            END DO
         END ASSOCIATE
         frame%dion(first:last, first:last) = kinds(ikind)%sgp_potential%upf_projectors%d
         CALL upf_radial_fields_periodic(frame%beta_cache(ikind), points, frame%positions(:, atom), &
                                         frame%cell, images, reduced_k, frame%beta(:, first:last))
         CALL upf_radial_fields_periodic(frame%delta_cache(ikind), points, frame%positions(:, atom), &
                                  frame%cell, images, reduced_k, frame%delta(:, first:last), frame%delta_gradient(:, first:last, :))
      END DO
      IF (use_atomic) THEN
         IF (.NOT. frame%direct) CALL check_atomic_domains(frame, info)
         IF (info /= 0) RETURN
         ALLOCATE (overlaps(nproj, nao))
         CALL atomic_coefficients(frame, overlaps)
         frame%coefficients(:, :) = overlaps
         frame%atom_local = .TRUE.
      ELSE
         ALLOCATE (coupling(nproj, nproj), overlaps(nproj, nao))
         CALL qs_upf_inverse_centered_pairs(frame, coupling, overlaps)
         CALL upf_inverse_init(coupling, frame%inverse, info)
         IF (info /= 0) RETURN
         CALL upf_inverse_project(frame%inverse, overlaps, frame%coefficients)
      END IF
      CALL upf_inverse_basis(frame%gaussian, frame%delta, frame%coefficients, frame%soft)
      DO axis = 1, 3
         CALL upf_inverse_basis(frame%gaussian_gradient(:, :, axis), frame%delta_gradient(:, :, axis), &
                                frame%coefficients, frame%soft_gradient(:, :, axis))
      END DO
   END SUBROUTINE qs_upf_inverse_prepare

! **************************************************************************************************
!> \brief Reject cross-center beta-Delta and Delta-Delta support in the local inverse.
!> \param frame Original native centers, fields and explicitly supplied periodic images.
!> \param info Zero for nonoverlap, -2 otherwise. No overlap term is silently discarded.
! **************************************************************************************************
   SUBROUTINE check_atomic_domains(frame, info)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      INTEGER, INTENT(OUT)                               :: info

      INTEGER                                            :: a, b, image, ka, kb
      REAL(dp)                                           :: displacement(3), ra, radius, rb

      info = 0
      DO a = 1, SIZE(frame%center_atoms)
         ka = frame%center_kinds(a)
         ra = frame%delta_cache(ka)%radius
         DO b = 1, SIZE(frame%center_atoms)
            kb = frame%center_kinds(b)
            rb = frame%delta_cache(kb)%radius
            IF (ra == 0.0_dp .AND. rb == 0.0_dp) CYCLE
            radius = ra + rb
            IF (rb > 0.0_dp) radius = MAX(radius, frame%beta_cache(ka)%radius + rb)
            IF (ra > 0.0_dp) radius = MAX(radius, ra + frame%beta_cache(kb)%radius)
            DO image = 1, SIZE(frame%images, 2)
               IF (a == b .AND. ALL(frame%images(:, image) == 0)) CYCLE
               displacement = frame%positions(:, frame%center_atoms(a)) - frame%positions(:, frame%center_atoms(b)) - &
                              MATMUL(frame%cell, REAL(frame%images(:, image), dp))
               IF (SUM(displacement**2) >= radius**2) CYCLE
               info = -2
               RETURN
            END DO
         END DO
      END DO
   END SUBROUTINE check_atomic_domains

! **************************************************************************************************
!> \brief Native atom-local hard projector coefficients and their common geometry response.
!>        The radial K is fixed. No grid measure, point motion or global inverse enters W.
!> \param frame Prepared original kinds, Gaussian bases, centers and physical-cell images.
!> \param overlaps New hard projector bras against all Gaussian columns.
!> \param response Optional joint response, required with all three geometry arguments.
!> \param atom_response Optional nuclear and Gaussian-center motion in the native atom order.
!> \param cell_response Optional lattice motion.
!> \param reduced_k_response Optional reduced Bloch-vector motion.
! **************************************************************************************************
   SUBROUTINE atomic_coefficients(frame, overlaps, response, atom_response, cell_response, reduced_k_response)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      COMPLEX(dp), INTENT(OUT)                           :: overlaps(:, :)
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: response(:, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: atom_response(:, :), &
                                                            cell_response(3, 3), &
                                                            reduced_k_response(3)

      INTEGER                                            :: atom, center, first, ikind, last, p, q

      CPASSERT(ALLOCATED(frame%atomic_cache))
      IF (PRESENT(response)) THEN
         CPASSERT(PRESENT(atom_response))
         CPASSERT(PRESENT(cell_response) .AND. PRESENT(reduced_k_response))
      END IF
      DO center = 1, SIZE(frame%center_atoms)
         ikind = frame%center_kinds(center)
         atom = frame%center_atoms(center)
         p = frame%projector_first(center)
         q = frame%projector_last(center)
         DO first = 1, SIZE(frame%bases)
            last = frame%ao_last(first)
            IF (PRESENT(response)) THEN
               IF (ALL(atom_response(:, atom) == atom_response(:, first)) .AND. &
                   ALL(cell_response == 0.0_dp) .AND. ALL(reduced_k_response == 0.0_dp)) THEN
                  overlaps(p:q, frame%ao_first(first):last) = frame%coefficients(p:q, frame%ao_first(first):last)
                  response(p:q, frame%ao_first(first):last) = CMPLX(0.0_dp, 0.0_dp, dp)
                  CYCLE
               END IF
               CALL qs_upf_inverse_atomic_overlap(frame%bases(first)%gto_basis_set, frame%atomic_cache(ikind)%hard, &
                                                  frame%positions(:, atom), frame%positions(:, first), frame%cell, frame%images, &
                                                  frame%reduced_k, overlaps(p:q, frame%ao_first(first):last), &
                                                  response(p:q, frame%ao_first(first):last), &
                                                 atom_response(:, atom), atom_response(:, first), cell_response, reduced_k_response)
            ELSE
               CALL qs_upf_inverse_atomic_overlap(frame%bases(first)%gto_basis_set, frame%atomic_cache(ikind)%hard, &
                                                  frame%positions(:, atom), frame%positions(:, first), frame%cell, frame%images, &
                                                  frame%reduced_k, overlaps(p:q, frame%ao_first(first):last))
            END IF
         END DO
      END DO
   END SUBROUTINE atomic_coefficients

! **************************************************************************************************
!> \brief Atomic corrections to native Gaussian kinetic and overlap matrices, retaining closure.
!>        The kinetic correction equals -W^dagger T W plus its finite-map residual.
!>        Adding the transformed original D gives the intended hard D-T representation.
!> \param frame Nonoverlapping atom-local inverse with native hard projector coefficients.
!> \param projector_metric Original metric selected by the functional, without substituting q_wave.
!> \param kinetic Owned correction to the ordinary Gaussian kinetic matrix.
!> \param metric Owned correction to the ordinary Gaussian overlap matrix.
!> \param kinetic_response Optional owned joint kinetic response, with metric_response and geometry.
!> \param metric_response Optional owned joint metric response.
!> \param atom_response Optional nuclear and Gaussian-center displacement.
!> \param cell_response Optional lattice displacement.
!> \param reduced_k_response Optional reduced Bloch-vector displacement.
!> \param coefficients_response Optional precomputed geometry derivative of the inverse projector coefficients.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_atomic_operators(frame, projector_metric, kinetic, metric, kinetic_response, metric_response, &
                                              atom_response, cell_response, reduced_k_response, coefficients_response)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      COMPLEX(dp), INTENT(IN)                            :: projector_metric(:, :)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: kinetic(:, :), metric(:, :)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT), OPTIONAL    :: kinetic_response(:, :), &
                                                            metric_response(:, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: atom_response(:, :), &
                                                            cell_response(3, 3), &
                                                            reduced_k_response(3)
      COMPLEX(dp), INTENT(IN), OPTIONAL, TARGET          :: coefficients_response(:, :)

      COMPLEX(dp), ALLOCATABLE :: cross(:, :), dcross(:, :), dl(:, :), dz(:, :), laplacian(:, :), &
         mixed(:, :), overlaps(:, :), rs(:, :), rt(:, :), z(:, :), zero_ao(:, :), &
         zero_projectors(:, :)
      COMPLEX(dp), ALLOCATABLE, TARGET                   :: owned_dw(:, :)
      COMPLEX(dp), POINTER                               :: dw(:, :)
      INTEGER                                            :: atom, b, center, first, ikind, last, n, &
                                                            np, p, q
      LOGICAL                                            :: calculate_response

      CPASSERT(ALLOCATED(frame%atomic_cache))
      n = SIZE(frame%coefficients, 2)
      np = SIZE(frame%coefficients, 1)
      CPASSERT(ALL(SHAPE(projector_metric) == [np, np]))
      ALLOCATE (kinetic(n, n), metric(n, n))
      kinetic = CMPLX(0.0_dp, 0.0_dp, dp)
      metric = CMPLX(0.0_dp, 0.0_dp, dp)
      calculate_response = PRESENT(kinetic_response) .OR. PRESENT(metric_response)
      IF (calculate_response) THEN
         CPASSERT(PRESENT(kinetic_response))
         CPASSERT(PRESENT(metric_response))
         CPASSERT(PRESENT(atom_response))
         CPASSERT(PRESENT(cell_response) .AND. PRESENT(reduced_k_response))
         ALLOCATE (kinetic_response(n, n), metric_response(n, n))
         kinetic_response = CMPLX(0.0_dp, 0.0_dp, dp)
         metric_response = CMPLX(0.0_dp, 0.0_dp, dp)
      END IF
      IF (frame%direct) RETURN
      IF (calculate_response) THEN
         ALLOCATE (zero_ao(n, n))
         zero_ao = CMPLX(0.0_dp, 0.0_dp, dp)
         IF (PRESENT(coefficients_response)) THEN
            CPASSERT(ALL(SHAPE(coefficients_response) == [np, n]))
            dw => coefficients_response
         ELSE
            CPASSERT(frame%atom_local)
            ALLOCATE (owned_dw(np, n), overlaps(np, n))
            dw => owned_dw
            CALL atomic_coefficients(frame, overlaps, dw, atom_response, cell_response, reduced_k_response)
         END IF
      END IF
      DO center = 1, SIZE(frame%center_atoms)
         ikind = frame%center_kinds(center)
         IF (.NOT. frame%atom_local) THEN
            IF (frame%delta_cache(ikind)%radius <= 0.0_dp) CYCLE
         END IF
         atom = frame%center_atoms(center)
         p = frame%projector_first(center)
         q = frame%projector_last(center)
         ALLOCATE (z(q - p + 1, n), laplacian(q - p + 1, n))
         IF (calculate_response) THEN
            ALLOCATE (dz(q - p + 1, n), dl(q - p + 1, n), zero_projectors(q - p + 1, q - p + 1))
            zero_projectors = CMPLX(0.0_dp, 0.0_dp, dp)
         END IF
         DO b = 1, SIZE(frame%bases)
            first = frame%ao_first(b)
            last = frame%ao_last(b)
            IF (calculate_response) THEN
               CALL qs_upf_inverse_atomic_overlap(frame%bases(b)%gto_basis_set, frame%atomic_cache(ikind)%difference, &
                                                  frame%positions(:, atom), frame%positions(:, b), frame%cell, frame%images, &
                                                  frame%reduced_k, z(:, first:last), dz(:, first:last), &
                                                  atom_response(:, atom), atom_response(:, b), cell_response, reduced_k_response, &
                                                  laplacian(:, first:last), dl(:, first:last))
            ELSE
               CALL qs_upf_inverse_atomic_overlap(frame%bases(b)%gto_basis_set, frame%atomic_cache(ikind)%difference, &
                                                  frame%positions(:, atom), frame%positions(:, b), frame%cell, frame%images, &
                                                  frame%reduced_k, z(:, first:last), kinetic=laplacian(:, first:last))
            END IF
         END DO
         ASSOCIATE (atomic => frame%atomic_cache(ikind), w => frame%coefficients(p:q, :))
            CALL upf_inverse_atomic_residual(atomic, w, z, laplacian, rt, rs)
            kinetic(:, :) = kinetic + rt - MATMUL(CONJG(TRANSPOSE(w)), MATMUL(atomic%kinetic, w))
            metric(:, :) = metric + rs + MATMUL(CONJG(TRANSPOSE(w)), &
                                                MATMUL(projector_metric(p:q, p:q) - atomic%wave_overlap, w))
            IF (calculate_response) THEN
               CALL upf_inverse_operator_response(zero_ao, CONJG(TRANSPOSE(laplacian)), laplacian, &
                                                  CMPLX(atomic%difference_kinetic, 0.0_dp, dp), zero_projectors, w, &
                                                  CONJG(TRANSPOSE(dl)), dl, zero_projectors, zero_projectors, dw(p:q, :), rt)
               CALL upf_inverse_operator_response(zero_ao, CONJG(TRANSPOSE(z)), z, &
                                                  CMPLX(atomic%difference_overlap, 0.0_dp, dp), projector_metric(p:q, p:q), w, &
                                                  CONJG(TRANSPOSE(dz)), dz, zero_projectors, zero_projectors, dw(p:q, :), rs)
               kinetic_response(:, :) = kinetic_response + rt
               metric_response(:, :) = metric_response + rs
            END IF
         END ASSOCIATE
         DEALLOCATE (z, laplacian, rt, rs)
         IF (calculate_response) DEALLOCATE (dz, dl, zero_projectors)
      END DO
      IF (frame%atom_local) RETURN
      ALLOCATE (cross(np, np))
      IF (calculate_response) THEN
         ALLOCATE (dcross(np, np), mixed(n, n))
         CALL qs_upf_inverse_cross_kinetic(frame, cross, dcross, atom_response, cell_response, reduced_k_response)
         mixed(:, :) = MATMUL(CONJG(TRANSPOSE(dw)), MATMUL(cross, frame%coefficients))
         kinetic_response(:, :) = kinetic_response + mixed + CONJG(TRANSPOSE(mixed)) + &
                                  MATMUL(CONJG(TRANSPOSE(frame%coefficients)), MATMUL(dcross, frame%coefficients))
      ELSE
         CALL qs_upf_inverse_cross_kinetic(frame, cross)
      END IF
      kinetic(:, :) = kinetic + MATMUL(CONJG(TRANSPOSE(frame%coefficients)), MATMUL(cross, frame%coefficients))
   END SUBROUTINE qs_upf_inverse_atomic_operators

! **************************************************************************************************
!> \brief Joint point/atom/cell/Bloch response in the same native AO and projector frame.
!> \param frame Prepared fields whose native bases remain alive and unchanged.
!> \param point_response Cartesian local point response.
!> \param atom_response Cartesian response of every native atom, including ghost/NC basis centers.
!> \param cell_response Lattice response at the same reduced Bloch convention.
!> \param reduced_k_response Reduced Bloch-vector response.
!> \param weights_response Response of the original local measure.
!> \param coefficients_response Response of W=K^-1 B.
!> \param soft_response Optional response of the inverse AO values.
!> \param gradient_response Optional response of their physical spatial gradients.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_response(frame, point_response, atom_response, cell_response, reduced_k_response, &
                                      weights_response, coefficients_response, soft_response, gradient_response)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      REAL(dp), INTENT(IN) :: point_response(:, :), atom_response(:, :), cell_response(3, 3), &
         reduced_k_response(3), weights_response(:)
      COMPLEX(dp), INTENT(OUT)                           :: coefficients_response(:, :)
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: soft_response(:, :), &
                                                            gradient_response(:, :, :)

      COMPLEX(dp), ALLOCATABLE :: chi(:, :), db(:, :), dchi(:, :), ddelta(:, :), delta(:, :), &
         dgc(:, :, :), dgd(:, :, :), dk(:, :), gc(:, :, :), gd(:, :, :)
      INTEGER                                            :: atom, axis, center, first, ikind, last, &
                                                            nao, npoint, nproj
      LOGICAL                                            :: shared_response

      CPASSERT(frame%atom_local .OR. frame%inverse%ready)
      CPASSERT(ALL(SHAPE(point_response) == SHAPE(frame%points)))
      CPASSERT(ALL(SHAPE(atom_response) == SHAPE(frame%positions)))
      CPASSERT(SIZE(weights_response) == SIZE(frame%weights))
      npoint = SIZE(frame%weights)
      nao = SIZE(frame%gaussian, 2)
      nproj = SIZE(frame%coefficients, 1)
      CPASSERT(ALL(SHAPE(coefficients_response) == [nproj, nao]))
      IF (PRESENT(soft_response)) THEN
         CPASSERT(ALL(SHAPE(soft_response) == [npoint, nao]))
      END IF
      IF (PRESENT(gradient_response)) THEN
         CPASSERT(ALL(SHAPE(gradient_response) == [npoint, nao, 3]))
      END IF
      ALLOCATE (delta(npoint, nproj), ddelta(npoint, nproj), &
                chi(npoint, nao), dchi(npoint, nao), dk(nproj, nproj), db(nproj, nao))
      IF (PRESENT(gradient_response)) THEN
         ALLOCATE (gd(npoint, nproj, 3), dgd(npoint, nproj, 3), gc(npoint, nao, 3), dgc(npoint, nao, 3))
      END IF
      shared_response = ANY(point_response /= 0.0_dp) .OR. ANY(cell_response /= 0.0_dp) .OR. &
                        ANY(reduced_k_response /= 0.0_dp)
      DO atom = 1, SIZE(frame%bases)
         first = frame%ao_first(atom)
         last = frame%ao_last(atom)
         IF (.NOT. shared_response .AND. ALL(atom_response(:, atom) == 0.0_dp)) THEN
            dchi(:, first:last) = CMPLX(0.0_dp, 0.0_dp, dp)
            IF (PRESENT(gradient_response)) dgc(:, first:last, :) = CMPLX(0.0_dp, 0.0_dp, dp)
            CYCLE
         END IF
         IF (PRESENT(gradient_response)) THEN
            CALL gaussian_basis_fields_periodic(frame%bases(atom)%gto_basis_set, frame%points, frame%positions(:, atom), &
                                              frame%cell, frame%images, frame%reduced_k, chi(:, first:last), gc(:, first:last, :), &
                                               dchi(:, first:last), dgc(:, first:last, :), point_response, atom_response(:, atom), &
                                                cell_response, reduced_k_response, &
                                                cutoff_squared=4.0_dp*frame%bases(atom)%gto_basis_set%kind_radius**2)
         ELSE
            CALL gaussian_basis_fields_periodic(frame%bases(atom)%gto_basis_set, frame%points, frame%positions(:, atom), &
                                                frame%cell, frame%images, frame%reduced_k, chi(:, first:last), &
                                                response=dchi(:, first:last), point_response=point_response, &
                                                center_response=atom_response(:, atom), cell_response=cell_response, &
                                                reduced_k_response=reduced_k_response, &
                                                cutoff_squared=4.0_dp*frame%bases(atom)%gto_basis_set%kind_radius**2)
         END IF
      END DO
      DO center = 1, SIZE(frame%center_atoms)
         atom = frame%center_atoms(center)
         ikind = frame%center_kinds(center)
         first = frame%projector_first(center)
         last = frame%projector_last(center)
         IF (frame%direct) THEN
            ddelta(:, first:last) = CMPLX(0.0_dp, 0.0_dp, dp)
            IF (PRESENT(gradient_response)) dgd(:, first:last, :) = CMPLX(0.0_dp, 0.0_dp, dp)
            CYCLE
         END IF
         IF (.NOT. shared_response .AND. ALL(atom_response(:, atom) == 0.0_dp)) THEN
            ddelta(:, first:last) = CMPLX(0.0_dp, 0.0_dp, dp)
            IF (PRESENT(gradient_response)) dgd(:, first:last, :) = CMPLX(0.0_dp, 0.0_dp, dp)
            CYCLE
         END IF
         IF (PRESENT(gradient_response)) THEN
            CALL upf_radial_fields_periodic(frame%delta_cache(ikind), frame%points, frame%positions(:, atom), &
                                            frame%cell, frame%images, frame%reduced_k, delta(:, first:last), gd(:, first:last, :), &
                                            ddelta(:, first:last), dgd(:, first:last, :), point_response, atom_response(:, atom), &
                                            cell_response, reduced_k_response)
         ELSE
            CALL upf_radial_fields_periodic(frame%delta_cache(ikind), frame%points, frame%positions(:, atom), &
                                            frame%cell, frame%images, frame%reduced_k, delta(:, first:last), &
                                            response=ddelta(:, first:last), points_response=point_response, &
                                            center_response=atom_response(:, atom), cell_response=cell_response, &
                                            reduced_k_response=reduced_k_response)
         END IF
      END DO
      IF (frame%atom_local) THEN
         CALL atomic_coefficients(frame, db, coefficients_response, atom_response, cell_response, reduced_k_response)
      ELSE
         CALL qs_upf_inverse_centered_pairs(frame, coupling_response=dk, overlaps_response=db, &
                                            atom_response=atom_response, cell_response=cell_response, &
                                            reduced_k_response=reduced_k_response)
         CALL upf_inverse_project_response(frame%inverse, frame%coefficients, db, dk, coefficients_response)
      END IF
      IF (PRESENT(soft_response)) THEN
         CALL upf_inverse_basis_response(dchi, frame%delta, ddelta, frame%coefficients, &
                                         coefficients_response, soft_response)
      END IF
      IF (PRESENT(gradient_response)) THEN
         DO axis = 1, 3
            CALL upf_inverse_basis_response(dgc(:, :, axis), frame%delta_gradient(:, :, axis), dgd(:, :, axis), &
                                            frame%coefficients, coefficients_response, gradient_response(:, :, axis))
         END DO
      END IF
   END SUBROUTINE qs_upf_inverse_response

! **************************************************************************************************
!> \brief Release owned fields and borrowed basis pointers before destroying their KS environment.
!> \param frame ...
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_release(frame)
      TYPE(qs_upf_inverse_type), INTENT(INOUT)           :: frame

      TYPE(qs_upf_inverse_type)                          :: empty

      frame = empty
   END SUBROUTINE qs_upf_inverse_release

! **************************************************************************************************
!> \brief Unfold each Bloch projector bra to its home atom's quadrature.
!>        The Gaussian and difference-wave kets retain all periodic images and interatomic
!>        couplings. The home offsets and weights are Cartesian and fixed under strain.
!>        The defining overlap is independent of the uniform density grid.
!> \param frame Original projector, Gaussian and difference-wave fields.
!> \param coupling Original beta-Delta integrals without the identity.
!> \param overlaps Original beta-Gaussian integrals.
!> \param coupling_response Joint geometry response of the beta-Delta integrals.
!> \param overlaps_response Joint geometry response of the beta-Gaussian integrals.
!> \param atom_response Physical atom motion, including motion of the projector quadrature.
!> \param cell_response Physical lattice response at fixed home offsets.
!> \param reduced_k_response Motion of the reduced Bloch vector.
!> \param image A single unphased ket image for the equivalent real-space convolution.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_centered_pairs(frame, coupling, overlaps, coupling_response, overlaps_response, &
                                            atom_response, cell_response, reduced_k_response, image)
      TYPE(qs_upf_inverse_type), INTENT(IN), TARGET      :: frame
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: coupling(:, :), overlaps(:, :), &
                                                            coupling_response(:, :), &
                                                            overlaps_response(:, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: atom_response(:, :), &
                                                            cell_response(3, 3), &
                                                            reduced_k_response(3)
      INTEGER, INTENT(IN), OPTIONAL                      :: image(3)

      COMPLEX(dp), ALLOCATABLE                           :: beta(:, :), ddelta(:, :), delta(:, :), &
                                                            pair_b(:, :), pair_k(:, :)
      INTEGER                                            :: atom, basis_atom, center, first, im, ir, &
                                                            kind, last, na, nao, ni, nproj, other, &
                                                            other_atom, p, q
      INTEGER, ALLOCATABLE                               :: image_counts(:)
      INTEGER, ALLOCATABLE, TARGET                       :: supported_images(:, :, :)
      INTEGER, POINTER                                   :: images(:, :), ket_images(:, :)
      INTEGER, TARGET                                    :: single_image(3, 1)
      LOGICAL                                            :: derivative, includes_home, &
                                                            shared_response, values
      REAL(dp)                                           :: bra_radius, distance(3), ket_radius, &
                                                            wavevector(3)
      REAL(dp), ALLOCATABLE :: angular_weights(:), directions(:, :), motion(:, :), points(:, :), &
         radial_beta(:), translations(:, :), weights(:)
      TYPE(grid_atom_type), POINTER                      :: grid

      values = PRESENT(coupling)
      derivative = PRESENT(coupling_response)
      wavevector = frame%reduced_k
      images => frame%images
      IF (PRESENT(image)) THEN
         single_image(:, 1) = image
         images => single_image
         wavevector = 0.0_dp
      END IF
      ALLOCATE (image_counts(SIZE(frame%center_atoms)), &
                supported_images(3, SIZE(images, 2), SIZE(frame%center_atoms)), &
                translations(3, SIZE(images, 2)))
      DO im = 1, SIZE(images, 2)
         translations(:, im) = MATMUL(frame%cell, REAL(images(:, im), dp))
      END DO
      includes_home = ANY(ALL(images == 0, DIM=1))
      IF (values) THEN
         coupling = CMPLX(0.0_dp, 0.0_dp, dp)
         overlaps = CMPLX(0.0_dp, 0.0_dp, dp)
      END IF
      shared_response = .FALSE.
      IF (derivative) THEN
         CPASSERT(PRESENT(atom_response) .AND. PRESENT(cell_response))
         CPASSERT(PRESENT(reduced_k_response))
         shared_response = ANY(cell_response /= 0.0_dp) .OR. ANY(reduced_k_response /= 0.0_dp)
         coupling_response = CMPLX(0.0_dp, 0.0_dp, dp)
         overlaps_response = CMPLX(0.0_dp, 0.0_dp, dp)
      END IF
      nao = SIZE(frame%coefficients, 2)
      nproj = SIZE(frame%coefficients, 1)
      DO center = frame%point_comm%mepos + 1, SIZE(frame%center_atoms), frame%point_comm%num_pe
         atom = frame%center_atoms(center)
         kind = frame%center_kinds(center)
         first = frame%projector_first(center)
         last = frame%projector_last(center)
         grid => frame%defining_atom_grids(kind)%grid
         CPASSERT(ASSOCIATED(grid))
         CALL centered_grid(frame%defining_atom_grids(kind), directions, angular_weights)
         na = SIZE(angular_weights)
         ALLOCATE (points(3, na), weights(na), radial_beta(last - first + 1), &
                   beta(na, last - first + 1), delta(na, nproj), &
                   pair_b(last - first + 1, nao), pair_k(last - first + 1, nproj))
         IF (derivative) THEN
            ALLOCATE (motion(3, na), ddelta(na, nproj))
            motion(:, :) = SPREAD(atom_response(:, atom), 2, na)
         END IF
         ! B uses the original radial weights and the analytic Gaussian angular integrals.
         DO basis_atom = 1, SIZE(frame%bases)
            ! A common translation leaves the overlap unchanged at fixed cell and reduced k.
            IF (.NOT. values .AND. derivative .AND. .NOT. shared_response) THEN
               IF (ALL(atom_response(:, basis_atom) == atom_response(:, atom))) CYCLE
            END IF
            p = frame%ao_first(basis_atom)
            q = frame%ao_last(basis_atom)
            ASSOCIATE (basis => frame%bases(basis_atom)%gto_basis_set, raw => frame%defining_atom_grids(kind)%projectors)
               IF (derivative) THEN
                  CALL qs_upf_inverse_atomic_overlap(basis, raw, frame%positions(:, atom), frame%positions(:, basis_atom), &
                                               frame%cell, images, wavevector, pair_b(:, p:q), overlaps_response(first:last, p:q), &
                                            atom_response(:, atom), atom_response(:, basis_atom), cell_response, reduced_k_response)
               ELSE
                  CALL qs_upf_inverse_atomic_overlap(basis, raw, frame%positions(:, atom), frame%positions(:, basis_atom), &
                                                     frame%cell, images, wavevector, pair_b(:, p:q))
               END IF
            END ASSOCIATE
            IF (values) overlaps(first:last, p:q) = pair_b(:, p:q)
         END DO
         ! The home block uses the same radial integral as the local inverse.
         ! Its atom, cell and Bloch-vector derivatives are zero.
         IF (values .AND. includes_home) THEN
            coupling(first:last, first:last) = frame%defining_atom_grids(kind)%home_coupling
         END IF
         bra_radius = frame%beta_cache(kind)%radius
         image_counts = 0
         DO other = 1, SIZE(frame%center_atoms)
            other_atom = frame%center_atoms(other)
            ket_radius = frame%delta_cache(frame%center_kinds(other))%radius
            IF (ket_radius <= 0.0_dp) CYCLE
            IF (.NOT. values .AND. derivative .AND. .NOT. shared_response) THEN
               IF (ALL(atom_response(:, other_atom) == atom_response(:, atom))) CYCLE
            END IF
            ni = 0
            DO im = 1, SIZE(images, 2)
               IF (other == center .AND. ALL(images(:, im) == 0)) CYCLE
               ! Disjoint compact radial supports contribute neither a value nor a first derivative.
               distance = frame%positions(:, atom) - (frame%positions(:, other_atom) + translations(:, im))
               IF (SUM(distance**2) > (bra_radius + ket_radius)**2) CYCLE
               ni = ni + 1
               supported_images(:, ni, other) = images(:, im)
            END DO
            image_counts(other) = ni
         END DO
         DO ir = 1, MERGE(grid%nr, 0, ANY(image_counts > 0))
            IF (grid%rad(ir) >= bra_radius) CYCLE
            CALL centered_shell(frame, center, grid%rad(ir), directions, points, beta, radial_beta)
            weights(:) = grid%wr(ir)*angular_weights
            DO other = 1, SIZE(frame%center_atoms)
               p = frame%projector_first(other)
               q = frame%projector_last(other)
               ni = image_counts(other)
               IF (ni == 0) THEN
                  IF (values) delta(:, p:q) = CMPLX(0.0_dp, 0.0_dp, dp)
                  IF (derivative) ddelta(:, p:q) = CMPLX(0.0_dp, 0.0_dp, dp)
                  CYCLE
               END IF
               ket_images => supported_images(:, :ni, other)
               ASSOCIATE (field => frame%delta_cache(frame%center_kinds(other)), other_atom => frame%center_atoms(other))
                  IF (derivative) THEN
                     IF (shared_response .OR. ANY(atom_response(:, other_atom) /= atom_response(:, atom))) THEN
                        CALL upf_radial_fields_periodic(field, points, frame%positions(:, other_atom), frame%cell, &
                                                        ket_images, wavevector, delta(:, p:q), response=ddelta(:, p:q), &
                                                        points_response=motion, center_response=atom_response(:, other_atom), &
                                                        cell_response=cell_response, reduced_k_response=reduced_k_response)
                        CYCLE
                     END IF
                     ddelta(:, p:q) = CMPLX(0.0_dp, 0.0_dp, dp)
                  END IF
                  IF (values) THEN
                     CALL upf_radial_fields_periodic(field, points, frame%positions(:, other_atom), frame%cell, &
                                                     ket_images, wavevector, delta(:, p:q))
                  END IF
               END ASSOCIATE
            END DO
            IF (values) THEN
               CALL upf_inverse_pair(beta, delta, weights, pair_k)
               coupling(first:last, :) = coupling(first:last, :) + pair_k
            END IF
            IF (derivative) THEN
               CALL upf_inverse_pair(beta, ddelta, weights, pair_k)
               coupling_response(first:last, :) = coupling_response(first:last, :) + pair_k
            END IF
         END DO
         DEALLOCATE (directions, angular_weights, points, weights, radial_beta, beta, delta, pair_b, pair_k)
         IF (derivative) DEALLOCATE (motion, ddelta)
      END DO
      IF (values) THEN
         CALL frame%point_comm%sum(overlaps)
         CALL frame%point_comm%sum(coupling)
      END IF
      IF (derivative) THEN
         CALL frame%point_comm%sum(overlaps_response)
         CALL frame%point_comm%sum(coupling_response)
      END IF
   END SUBROUTINE qs_upf_inverse_centered_pairs

! **************************************************************************************************
!> \brief Kinetic overlap of compact difference waves on distinct centers or images.
!>        The home bra and its quadrature move together. Both bra quadratures enter the
!>        Hermitian average, including for a pair of periodic images of the same atom.
!> \param frame Coupled inverse with the native atom-centered quadratures.
!> \param matrix Intercenter part of the Bloch difference-wave kinetic matrix.
!> \param response Optional joint derivative of the same finite quadrature.
!> \param atom_response Cartesian atom displacements.
!> \param cell_response Lattice displacement at fixed Cartesian home offsets.
!> \param reduced_k_response Displacement of the reduced Bloch vector.
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_cross_kinetic(frame, matrix, response, atom_response, cell_response, reduced_k_response)
      TYPE(qs_upf_inverse_type), INTENT(IN), TARGET      :: frame
      COMPLEX(dp), INTENT(OUT)                           :: matrix(:, :)
      COMPLEX(dp), INTENT(OUT), OPTIONAL                 :: response(:, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: atom_response(:, :), &
                                                            cell_response(3, 3), &
                                                            reduced_k_response(3)

      COMPLEX(dp), ALLOCATABLE                           :: bra(:, :, :), dket(:, :, :), &
                                                            ket(:, :, :), pair(:, :), value(:, :)
      INTEGER                                            :: a, axis, b, ia, image, ir, ka, kb, na, &
                                                            ni, p, q, u, v
      INTEGER, ALLOCATABLE                               :: images(:, :)
      REAL(dp)                                           :: displacement(3), radius, radius_sum
      REAL(dp), ALLOCATABLE                              :: angular_weights(:), directions(:, :), &
                                                            gradient(:, :), motion(:, :), &
                                                            points(:, :), scratch(:), weights(:)
      TYPE(grid_atom_type), POINTER                      :: grid
      TYPE(qs_upf_inverse_atom_grid_ref_type), POINTER   :: rule

      CPASSERT(ALLOCATED(frame%defining_atom_grids))
      IF (PRESENT(response)) THEN
         CPASSERT(PRESENT(atom_response) .AND. PRESENT(cell_response))
         CPASSERT(PRESENT(reduced_k_response))
         response = CMPLX(0.0_dp, 0.0_dp, dp)
      END IF
      matrix = CMPLX(0.0_dp, 0.0_dp, dp)
      ALLOCATE (images(3, SIZE(frame%images, 2)))
      DO a = 1, SIZE(frame%center_atoms)
         ka = frame%center_kinds(a)
         radius = frame%delta_cache(ka)%radius
         IF (radius <= 0.0_dp) CYCLE
         p = frame%projector_first(a)
         q = frame%projector_last(a)
         rule => frame%defining_atom_grids(ka)
         grid => rule%grid
         CALL centered_grid(rule, directions, angular_weights)
         na = SIZE(angular_weights)
         ALLOCATE (bra(na, q - p + 1, 3), gradient(3, q - p + 1), scratch(q - p + 1), &
                   points(3, na), weights(na))
         IF (PRESENT(response)) THEN
            ALLOCATE (motion(3, na))
            motion(:, :) = SPREAD(atom_response(:, frame%center_atoms(a)), 2, na)
         END IF
         DO b = 1, SIZE(frame%center_atoms)
            kb = frame%center_kinds(b)
            IF (frame%delta_cache(kb)%radius <= 0.0_dp) CYCLE
            radius_sum = radius + frame%delta_cache(kb)%radius
            ni = 0
            DO image = 1, SIZE(frame%images, 2)
               IF (a == b .AND. ALL(frame%images(:, image) == 0)) CYCLE
               displacement = frame%positions(:, frame%center_atoms(a)) - frame%positions(:, frame%center_atoms(b)) - &
                              MATMUL(frame%cell, REAL(frame%images(:, image), dp))
               IF (SUM(displacement**2) >= radius_sum**2) CYCLE
               ni = ni + 1
               images(:, ni) = frame%images(:, image)
            END DO
            IF (ni == 0) CYCLE
            u = frame%projector_first(b)
            v = frame%projector_last(b)
            ALLOCATE (value(na, v - u + 1), ket(na, v - u + 1, 3), pair(q - p + 1, v - u + 1))
            IF (PRESENT(response)) ALLOCATE (dket, MOLD=ket)
            DO ir = 1, grid%nr
               IF (grid%rad(ir) >= radius) CYCLE
               points(:, :) = SPREAD(frame%positions(:, frame%center_atoms(a)), 2, na) + grid%rad(ir)*directions
               weights(:) = 0.5_dp*grid%wr(ir)*angular_weights
               DO ia = 1, na
                  CALL upf_radial_fields_evaluate(frame%delta_cache(ka), grid%rad(ir)*directions(:, ia), scratch, gradient)
                  bra(ia, :, :) = TRANSPOSE(gradient)
               END DO
               IF (PRESENT(response)) THEN
                  CALL upf_radial_fields_periodic(frame%delta_cache(kb), points, frame%positions(:, frame%center_atoms(b)), &
                                                  frame%cell, images(:, :ni), frame%reduced_k, value, gradient=ket, &
                                                  gradient_response=dket, points_response=motion, &
                                                  center_response=atom_response(:, frame%center_atoms(b)), &
                                                  cell_response=cell_response, reduced_k_response=reduced_k_response)
               ELSE
                  CALL upf_radial_fields_periodic(frame%delta_cache(kb), points, frame%positions(:, frame%center_atoms(b)), &
                                                  frame%cell, images(:, :ni), frame%reduced_k, value, gradient=ket)
               END IF
               DO axis = 1, 3
                  CALL upf_inverse_pair(bra(:, :, axis), ket(:, :, axis), weights, pair)
                  matrix(p:q, u:v) = matrix(p:q, u:v) + pair
                  IF (.NOT. PRESENT(response)) CYCLE
                  CALL upf_inverse_pair(bra(:, :, axis), dket(:, :, axis), weights, pair)
                  response(p:q, u:v) = response(p:q, u:v) + pair
               END DO
            END DO
            DEALLOCATE (value, ket, pair)
            IF (PRESENT(response)) DEALLOCATE (dket)
         END DO
         IF (PRESENT(response)) DEALLOCATE (motion)
         DEALLOCATE (bra, gradient, scratch, points, weights, directions, angular_weights)
      END DO
      matrix = 0.5_dp*(matrix + CONJG(TRANSPOSE(matrix)))
      IF (PRESENT(response)) response = 0.5_dp*(response + CONJG(TRANSPOSE(response)))
   END SUBROUTINE qs_upf_inverse_cross_kinetic

! **************************************************************************************************
!> \brief Select the same angular weights for the defining integrals and their adjoints.
!> \param rule Native atom grid and optional rotated directions.
!> \param directions ...
!> \param weights ...
! **************************************************************************************************
   SUBROUTINE centered_grid(rule, directions, weights)
      TYPE(qs_upf_inverse_atom_grid_ref_type), &
         INTENT(IN)                                      :: rule
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: directions(:, :), weights(:)

      IF (ALLOCATED(rule%directions)) THEN
         ALLOCATE (directions, SOURCE=rule%directions)
         ALLOCATE (weights, SOURCE=rule%angular_weights)
      ELSE
         ASSOCIATE (grid => rule%grid)
            ALLOCATE (directions(3, grid%ng_sphere), weights(grid%ng_sphere))
            directions(1, :) = grid%sin_pol*grid%cos_azi
            directions(2, :) = grid%sin_pol*grid%sin_azi
            directions(3, :) = grid%cos_pol
            weights(:) = grid%wa
         END ASSOCIATE
      END IF
   END SUBROUTINE centered_grid

! **************************************************************************************************
!> \brief Evaluate one home-projector shell with fixed Cartesian offsets.
!> \param frame ...
!> \param center ...
!> \param radius ...
!> \param directions ...
!> \param points ...
!> \param beta ...
!> \param values Scratch for one point.
! **************************************************************************************************
   SUBROUTINE centered_shell(frame, center, radius, directions, points, beta, values)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      INTEGER, INTENT(IN)                                :: center
      REAL(dp), INTENT(IN)                               :: radius, directions(:, :)
      REAL(dp), INTENT(OUT)                              :: points(:, :)
      COMPLEX(dp), INTENT(OUT)                           :: beta(:, :)
      REAL(dp), INTENT(OUT)                              :: values(:)

      INTEGER                                            :: ia, kind

      kind = frame%center_kinds(center)
      points(:, :) = SPREAD(frame%positions(:, frame%center_atoms(center)), 2, SIZE(directions, 2)) + radius*directions
      DO ia = 1, SIZE(directions, 2)
         CALL upf_radial_fields_evaluate(frame%beta_cache(kind), radius*directions(:, ia), values)
         beta(ia, :) = values
      END DO
   END SUBROUTINE centered_shell

! **************************************************************************************************
!> \brief Reuse or initialize the radial and angular quadrature owned by the kind.
!> \param kind Native kind with the requested radial and angular grid sizes.
!> \param quadrature Native radial quadrature identifier.
!> \param grid Borrowed grid. The caller initializes the Lebedev table before this call.
! **************************************************************************************************
   SUBROUTINE qs_upf_atom_kind_grid(kind, quadrature, grid)
      TYPE(qs_kind_type), INTENT(INOUT)                  :: kind
      INTEGER, INTENT(IN)                                :: quadrature
      TYPE(grid_atom_type), POINTER                      :: grid

      INTEGER                                            :: ll, na, nr

      CALL get_qs_kind(kind, grid_atom=grid, ngrid_ang=na, ngrid_rad=nr)
      IF (.NOT. ASSOCIATED(grid)) THEN
         CALL allocate_grid_atom(kind%grid_atom)
         grid => kind%grid_atom
      END IF
      IF (ASSOCIATED(grid%weight)) RETURN
      ll = get_number_of_lebedev_grid(n=na)
      na = lebedev_grid(ll)%n
      grid%nr = nr
      grid%ng_sphere = na
      CALL create_grid_atom(grid, nr, na, 0, ll, quadrature)
   END SUBROUTINE qs_upf_atom_kind_grid

! **************************************************************************************************
!> \brief Extend the native atom directions by the reference group when needed.
!> \param grid Native angular rule and weights.
!> \param rotations Cartesian reference rotations.
!> \param directions Expanded unit directions, unallocated when the native rule already suffices.
!> \param weights Expanded weights with the original angular normalization.
! **************************************************************************************************
   SUBROUTINE qs_upf_atom_angular_average(grid, rotations, directions, weights)
      TYPE(grid_atom_type), INTENT(IN)                   :: grid
      REAL(dp), INTENT(IN)                               :: rotations(:, :, :)
      REAL(dp), ALLOCATABLE, INTENT(OUT)                 :: directions(:, :), weights(:)

      REAL(dp)                                           :: original(3, grid%ng_sphere)

      original(1, :) = grid%sin_pol*grid%cos_azi
      original(2, :) = grid%sin_pol*grid%sin_azi
      original(3, :) = grid%cos_pol
      CALL lebedev_average_points(original, grid%wa, directions, weights, rotations)
      IF (SIZE(weights) == grid%ng_sphere) DEALLOCATE (directions, weights)
   END SUBROUTINE qs_upf_atom_angular_average

END MODULE qs_upf_inverse
