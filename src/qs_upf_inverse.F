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
   USE gaussian_basis_fields,           ONLY: gaussian_basis_fields_periodic
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: twopi
   USE message_passing,                 ONLY: mp_comm_type
   USE orbital_pointers,                ONLY: coset,&
                                              ncoset
   USE particle_types,                  ONLY: particle_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_upf_inverse_types,            ONLY: qs_upf_inverse_type
   USE qs_uspp_types,                   ONLY: qs_uspp_type
   USE upf_inverse,                     ONLY: &
        upf_inverse_atomic_init, upf_inverse_atomic_residual, upf_inverse_basis, &
        upf_inverse_basis_response, upf_inverse_fields_init, upf_inverse_init, &
        upf_inverse_operator_response, upf_inverse_project, upf_inverse_project_response
   USE upf_inverse_quadrature,          ONLY: upf_inverse_pair,&
                                              upf_inverse_pair_response
   USE upf_projector_integrals,         ONLY: upf_projector_overlap,&
                                              upf_projector_type
   USE upf_radial_fields,               ONLY: upf_radial_fields_periodic
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

      TYPE(cell_type), POINTER                           :: cell
      TYPE(particle_type), POINTER                       :: particles(:)
      TYPE(qs_kind_type), POINTER                        :: kinds(:)
      TYPE(qs_uspp_type), POINTER                        :: uspp

      CALL get_ks_env(ks_env, qs_kind_set=kinds, particle_set=particles, cell=cell, uspp=uspp)
      CPASSERT(ASSOCIATED(kinds) .AND. ASSOCIATED(particles) .AND. ASSOCIATED(cell))
      CPASSERT(ASSOCIATED(uspp))
CALL qs_upf_inverse_prepare(kinds, particles, cell, uspp, points, weights, images, reduced_k, comm, frame, info, atom_local, direct)
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
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_prepare(kinds, particles, cell, uspp, points, weights, images, reduced_k, comm, &
                                     frame, info, atom_local, direct)
      TYPE(qs_kind_type), INTENT(IN)                     :: kinds(:)
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
      IF (use_atomic) ALLOCATE (frame%atomic_cache(SIZE(kinds)))
      nproj = 0
      DO center = 1, ncenter
         ikind = frame%center_kinds(center)
         atom = frame%center_atoms(center)
         CPASSERT(.NOT. kinds(ikind)%ghost)
         CPASSERT(particles(atom)%atomic_kind%kind_number == ikind)
         IF (frame%beta_cache(ikind)%nproj == 0) THEN
            CALL upf_inverse_fields_init(kinds(ikind)%sgp_potential%upf_data, &
                                         frame%beta_cache(ikind), frame%delta_cache(ikind), frame%direct)
            IF (use_atomic) THEN
            CALL upf_inverse_atomic_init(kinds(ikind)%sgp_potential%upf_data, degree, frame%atomic_cache(ikind), info, frame%direct)
               IF (info /= 0) RETURN
            END IF
         END IF
         p = kinds(ikind)%sgp_potential%upf_projectors%nproj
         CPASSERT(frame%beta_cache(ikind)%nproj == p)
         frame%projector_first(center) = nproj + 1
         nproj = nproj + frame%beta_cache(ikind)%nproj
         frame%projector_last(center) = nproj
      END DO
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
         CALL upf_inverse_pair(frame%beta, frame%delta, weights, coupling, comm)
         CALL upf_inverse_pair(frame%beta, frame%gaussian, weights, overlaps, comm)
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
! **************************************************************************************************
   SUBROUTINE qs_upf_inverse_atomic_operators(frame, projector_metric, kinetic, metric, kinetic_response, metric_response, &
                                              atom_response, cell_response, reduced_k_response)
      TYPE(qs_upf_inverse_type), INTENT(IN)              :: frame
      COMPLEX(dp), INTENT(IN)                            :: projector_metric(:, :)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT)              :: kinetic(:, :), metric(:, :)
      COMPLEX(dp), ALLOCATABLE, INTENT(OUT), OPTIONAL    :: kinetic_response(:, :), &
                                                            metric_response(:, :)
      REAL(dp), INTENT(IN), OPTIONAL                     :: atom_response(:, :), &
                                                            cell_response(3, 3), &
                                                            reduced_k_response(3)

      COMPLEX(dp), ALLOCATABLE :: dl(:, :), dw(:, :), dz(:, :), laplacian(:, :), overlaps(:, :), &
         rs(:, :), rt(:, :), z(:, :), zero_ao(:, :), zero_projectors(:, :)
      INTEGER                                            :: atom, b, center, first, ikind, last, n, &
                                                            np, p, q
      LOGICAL                                            :: calculate_response

      CPASSERT(frame%atom_local .AND. ALLOCATED(frame%atomic_cache))
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
         ALLOCATE (dw(np, n), overlaps(np, n), zero_ao(n, n))
         zero_ao = CMPLX(0.0_dp, 0.0_dp, dp)
         CALL atomic_coefficients(frame, overlaps, dw, atom_response, cell_response, reduced_k_response)
      END IF
      DO center = 1, SIZE(frame%center_atoms)
         ikind = frame%center_kinds(center)
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

      COMPLEX(dp), ALLOCATABLE :: beta(:, :), chi(:, :), db(:, :), dbeta(:, :), dchi(:, :), &
         ddelta(:, :), delta(:, :), dgc(:, :, :), dgd(:, :, :), dk(:, :), gc(:, :, :), gd(:, :, :)
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
      ALLOCATE (beta(npoint, nproj), delta(npoint, nproj), dbeta(npoint, nproj), ddelta(npoint, nproj), &
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
            dbeta(:, first:last) = CMPLX(0.0_dp, 0.0_dp, dp)
            ddelta(:, first:last) = CMPLX(0.0_dp, 0.0_dp, dp)
            IF (PRESENT(gradient_response)) dgd(:, first:last, :) = CMPLX(0.0_dp, 0.0_dp, dp)
            CYCLE
         END IF
         IF (.NOT. frame%atom_local) THEN
            CALL upf_radial_fields_periodic(frame%beta_cache(ikind), frame%points, frame%positions(:, atom), &
                                    frame%cell, frame%images, frame%reduced_k, beta(:, first:last), response=dbeta(:, first:last), &
                                            points_response=point_response, center_response=atom_response(:, atom), &
                                            cell_response=cell_response, reduced_k_response=reduced_k_response)
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
         CALL upf_inverse_pair_response(frame%beta, frame%delta, frame%weights, dbeta, ddelta, &
                                        weights_response, dk, frame%point_comm)
         CALL upf_inverse_pair_response(frame%beta, frame%gaussian, frame%weights, dbeta, dchi, &
                                        weights_response, db, frame%point_comm)
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

END MODULE qs_upf_inverse
