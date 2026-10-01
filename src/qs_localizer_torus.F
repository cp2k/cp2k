!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Full Born-von-Karman AO operators at a frozen SCF potential.
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
MODULE qs_localizer_torus
   USE cell_types,                      ONLY: cell_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_add, dbcsr_create, dbcsr_desymmetrize, dbcsr_distribution_get, &
        dbcsr_distribution_new, dbcsr_distribution_release, dbcsr_distribution_type, &
        dbcsr_finalize, dbcsr_get_info, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_p_type, &
        dbcsr_put_block, dbcsr_release, dbcsr_set, dbcsr_type, dbcsr_type_antisymmetric, &
        dbcsr_type_no_symmetry, dbcsr_type_symmetric
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: dbcsr_deallocate_matrix_set
   USE kinds,                           ONLY: dp
   USE kpoint_methods,                  ONLY: kpoint_init_cell_index,&
                                              rskp_transform
   USE kpoint_types,                    ONLY: kpoint_create,&
                                              kpoint_release,&
                                              kpoint_type
   USE mathconstants,                   ONLY: twopi
   USE message_passing,                 ONLY: mp_para_env_type
   USE qs_ao_translations,              ONLY: ao_translation_matrices
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_env_release,&
                                              qs_environment_type
   USE qs_finite_ao,                    ONLY: finite_ao_soc
   USE qs_gamma2kp,                     ONLY: create_kp_from_gamma
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_moments,                      ONLY: build_berry_kpoint_matrix
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_localizer_torus'
   PUBLIC :: localizer_torus_operators, torus_translation_operators
CONTAINS

! **************************************************************************************************
!> \brief Lift scalar image operators and ordered Berry integrals to a full real-cell torus.
!> \param qs_env Converged SCF environment, possibly symmetry reduced
!> \param grid Gamma-centered property mesh, one in nonperiodic directions
!> \param directions Oriented periodic directions
!> \param do_soc Include physical-spin GTH SOC
!> \param h Newly allocated scalar Hamiltonians on the torus
!> \param s Newly allocated torus overlap
!> \param cosines Newly allocated torus cosine matrices, in directions order
!> \param sines Newly allocated sine matrices, indexed by lattice direction
!> \param soc Newly allocated torus SOC components when requested
!> \note No SCF weights or reconstructed occupied MOs enter. Keeping a real-cell AO basis
!>       also retains the standard physical time-reversal representation for the AII solver.
! **************************************************************************************************
   SUBROUTINE localizer_torus_operators(qs_env, grid, directions, do_soc, h, s, cosines, sines, soc)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: grid(3), directions(:)
      LOGICAL, INTENT(IN)                                :: do_soc
      TYPE(dbcsr_p_type), POINTER                        :: h(:), s(:), cosines(:), sines(:), &
                                                            soc(:, :)

      INTEGER                                            :: axis, j, nimages
      LOGICAL                                            :: do_kpoints
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), POINTER                        :: cb(:, :), sb(:, :)
      TYPE(kpoint_type), POINTER                         :: berry
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), POINTER            :: all_pairs(:)
      TYPE(qs_environment_type), POINTER                 :: qp

      NULLIFY (qp, berry, cb, sb, soc)
      CPASSERT(SIZE(directions) >= 1 .AND. SIZE(directions) <= 3)
      CPASSERT(ALL(directions >= 1) .AND. ALL(directions <= 3))
      DO j = 1, SIZE(directions)
         CPASSERT(COUNT(directions == directions(j)) == 1)
      END DO
      CALL get_qs_env(qs_env, do_kpoints=do_kpoints)
      IF (do_kpoints) THEN
         qp => qs_env
      ELSE
         CALL create_kp_from_gamma(qs_env, qp)
      END IF
      CALL get_qs_env(qp, cell=cell, para_env=para_env, sab_all=all_pairs)
      CALL lift_electronic_operators(qp, grid, do_soc, h, s, soc)
      ALLOCATE (cosines(SIZE(directions)), sines(3))
      CALL kpoint_create(berry)
      CALL kpoint_init_cell_index(berry, all_pairs, para_env, nimages)
      DO axis = 1, 3
         ALLOCATE (sines(axis)%matrix)
         CALL dbcsr_create(sines(axis)%matrix, template=s(1)%matrix, name='TORUS SIN')
         CALL dbcsr_finalize(sines(axis)%matrix)
      END DO
      DO j = 1, SIZE(directions)
         axis = directions(j)
         CALL build_berry_kpoint_matrix(qp, cb, sb, twopi*cell%h_inv(axis, :)/REAL(grid(axis), dp), &
                                        ordered=.TRUE., ordered_kpoints=berry)
         ALLOCATE (cosines(j)%matrix)
         CALL dbcsr_create(cosines(j)%matrix, template=s(1)%matrix, name='TORUS COS')
         CALL lift_berry(cb, sb, berry%index_to_cell, grid, axis, cosines(j)%matrix, sines(axis)%matrix)
         CALL dbcsr_deallocate_matrix_set(cb)
         CALL dbcsr_deallocate_matrix_set(sb)
      END DO
      CALL kpoint_release(berry)
      IF (.NOT. do_kpoints) THEN
         CALL qs_env_release(qp)
         DEALLOCATE (qp)
      END IF
   END SUBROUTINE localizer_torus_operators

! **************************************************************************************************
!> \brief Full AO torus with translated, periodized Gaussians instead of position coordinates.
!> \param qs_env Converged SCF environment
!> \param grid Full Gamma-centered property mesh
!> \param shifts Cartesian translations in internal bohr
!> \param do_soc Include post-SCF physical-spin GTH SOC
!> \param h Newly allocated scalar torus Hamiltonians
!> \param s Newly allocated torus metric
!> \param translations Newly allocated nonsymmetric torus translation overlaps
!> \param soc Newly allocated torus SOC components when requested
! **************************************************************************************************
   SUBROUTINE torus_translation_operators(qs_env, grid, shifts, do_soc, h, s, translations, soc)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: grid(3)
      REAL(dp), INTENT(IN)                               :: shifts(:, :)
      LOGICAL, INTENT(IN)                                :: do_soc
      TYPE(dbcsr_p_type), POINTER                        :: h(:), s(:), translations(:), soc(:, :)

      INTEGER                                            :: ir, it, r(3)
      LOGICAL                                            :: do_kpoints
      REAL(dp), ALLOCATABLE                              :: displaced(:, :)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), POINTER                        :: primitive_s(:, :), residues(:)
      TYPE(qs_environment_type), POINTER                 :: qp

      CALL get_qs_env(qs_env, do_kpoints=do_kpoints)
      IF (do_kpoints) THEN
         qp => qs_env
      ELSE
         CALL create_kp_from_gamma(qs_env, qp)
      END IF
      CALL get_qs_env(qp, cell=cell, matrix_s_kp=primitive_s)
      CALL lift_electronic_operators(qp, grid, do_soc, h, s, soc)
      ALLOCATE (translations(SIZE(shifts, 2)), displaced(3, SIZE(shifts, 2)))
      DO it = 1, SIZE(translations)
         ALLOCATE (translations(it)%matrix)
         CALL create_torus(primitive_s(1, 1)%matrix, PRODUCT(grid), translations(it)%matrix)
      END DO
      DO ir = 1, PRODUCT(grid)
         r = cell_coordinates(ir, grid)
         DO it = 1, SIZE(translations)
            displaced(:, it) = shifts(:, it) - MATMUL(cell%hmat, REAL(r, dp))
         END DO
         CALL ao_translation_matrices(qp, displaced, residues, periodic_grid=grid)
         DO it = 1, SIZE(translations)
            CALL lift_blocks(residues(it)%matrix, grid, r, 0, .FALSE., translations(it)%matrix)
         END DO
         CALL dbcsr_deallocate_matrix_set(residues)
      END DO
      IF (.NOT. do_kpoints) THEN
         CALL qs_env_release(qp)
         DEALLOCATE (qp)
      END IF
   END SUBROUTINE torus_translation_operators

! **************************************************************************************************
!> \brief Common Hamiltonian, metric and physical-spin image lifting for the two AO-torus analyses.
!> \param qp Environment with explicit real-space image matrices
!> \param grid Property mesh
!> \param do_soc Include post-SCF GTH SOC
!> \param h Newly allocated scalar Hamiltonians
!> \param s Newly allocated overlap
!> \param soc Newly allocated SOC components, or null without SOC
! **************************************************************************************************
   SUBROUTINE lift_electronic_operators(qp, grid, do_soc, h, s, soc)
      TYPE(qs_environment_type), POINTER                 :: qp
      INTEGER, INTENT(IN)                                :: grid(3)
      LOGICAL, INTENT(IN)                                :: do_soc
      TYPE(dbcsr_p_type), POINTER                        :: h(:), s(:), soc(:, :)

      INTEGER                                            :: axis, spin
      TYPE(dbcsr_p_type), POINTER                        :: hk(:, :), sk(:, :), vk(:, :)
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(neighbor_list_set_p_type), POINTER            :: pairs(:)
      TYPE(qs_ks_env_type), POINTER                      :: ks_env

      NULLIFY (soc, vk)
      CALL get_qs_env(qp, ks_env=ks_env, matrix_ks_kp=hk, matrix_s_kp=sk, sab_orb=pairs)
      CALL get_ks_env(ks_env, kpoints=kpoints)
      ALLOCATE (h(SIZE(hk, 1)), s(1))
      DO spin = 1, SIZE(hk, 1)
         ALLOCATE (h(spin)%matrix)
         CALL lift_fourier(hk, spin, pairs, kpoints, grid, .FALSE., h(spin)%matrix)
      END DO
      ALLOCATE (s(1)%matrix)
      CALL lift_fourier(sk, 1, pairs, kpoints, grid, .FALSE., s(1)%matrix)
      IF (do_soc) THEN
         CALL finite_ao_soc(qp, vk)
         ALLOCATE (soc(3, 1))
         DO axis = 1, 3
            ALLOCATE (soc(axis, 1)%matrix)
            CALL lift_fourier(vk, axis, pairs, kpoints, grid, .TRUE., soc(axis, 1)%matrix)
         END DO
         CALL dbcsr_deallocate_matrix_set(vk)
      END IF
   END SUBROUTINE lift_electronic_operators

! **************************************************************************************************
!> \brief Repeat the atom-block distribution; cell copies remain on the same MPI owners.
!> \param source Primitive AO template
!> \param copies Number of real cells
!> \param matrix New full-storage torus matrix
! **************************************************************************************************
   SUBROUTINE create_torus(source, copies, matrix)
      TYPE(dbcsr_type), INTENT(IN)                       :: source
      INTEGER, INTENT(IN)                                :: copies
      TYPE(dbcsr_type), INTENT(OUT)                      :: matrix

      INTEGER                                            :: n
      INTEGER, POINTER                                   :: cols(:), ec(:), er(:), es(:), rows(:), &
                                                            sizes(:)
      TYPE(dbcsr_distribution_type)                      :: distribution, expanded

      CALL dbcsr_get_info(source, distribution=distribution, row_blk_size=sizes)
      CALL dbcsr_distribution_get(distribution, row_dist=rows, col_dist=cols)
      n = SIZE(sizes)*copies
      ALLOCATE (es(n), er(n), ec(n))
      es = RESHAPE(SPREAD(sizes, 2, copies), [n])
      er = RESHAPE(SPREAD(rows, 2, copies), [n])
      ec = RESHAPE(SPREAD(cols, 2, copies), [n])
      CALL dbcsr_distribution_new(expanded, template=distribution, row_dist=er, col_dist=ec)
      CALL dbcsr_create(matrix, 'TORUS AO', expanded, dbcsr_type_no_symmetry, es, es)
      CALL dbcsr_distribution_release(expanded)
      DEALLOCATE (es, er, ec)
   END SUBROUTINE create_torus

! **************************************************************************************************
!> \brief Inverse discrete Fourier transform of the existing CP2K image convention.
!> \param images Primitive real-space matrices
!> \param spin Component of images
!> \param pairs SCF orbital neighbour list
!> \param kpoints SCF image indexing, not the SCF k-points or weights
!> \param grid Full Gamma-centered property mesh
!> \param antisymmetric Whether images are real SOC components (before multiplication by i)
!> \param matrix New real-cell torus matrix
!> \note This reference assembly uses O(Nk**2) primitive-block operations, but never
!>       gathers a dense AO matrix. It intentionally shares rskp_transform conventions.
! **************************************************************************************************
   SUBROUTINE lift_fourier(images, spin, pairs, kpoints, grid, antisymmetric, matrix)
      TYPE(dbcsr_p_type), POINTER                        :: images(:, :)
      INTEGER, INTENT(IN)                                :: spin
      TYPE(neighbor_list_set_p_type), POINTER            :: pairs(:)
      TYPE(kpoint_type), POINTER                         :: kpoints
      INTEGER, INTENT(IN)                                :: grid(3)
      LOGICAL, INTENT(IN)                                :: antisymmetric
      TYPE(dbcsr_type), INTENT(OUT)                      :: matrix

      INTEGER                                            :: ik, ir, n, r(3)
      REAL(dp)                                           :: k(3), phase
      TYPE(dbcsr_p_type), POINTER                        :: residues(:)
      TYPE(dbcsr_type)                                   :: full_im, full_re, im, re

      n = PRODUCT(grid)
      CALL create_torus(images(spin, 1)%matrix, n, matrix)
      CALL dbcsr_create(re, template=images(spin, 1)%matrix, matrix_type=dbcsr_type_symmetric)
      CALL dbcsr_create(im, template=images(spin, 1)%matrix, matrix_type=dbcsr_type_antisymmetric)
      CALL cp_dbcsr_alloc_block_from_nbl(re, pairs)
      CALL cp_dbcsr_alloc_block_from_nbl(im, pairs)
      CALL dbcsr_create(full_re, template=re, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_create(full_im, template=re, matrix_type=dbcsr_type_no_symmetry)
      ALLOCATE (residues(n))
      DO ir = 1, n
         ALLOCATE (residues(ir)%matrix)
         CALL dbcsr_create(residues(ir)%matrix, template=full_re)
         CALL dbcsr_finalize(residues(ir)%matrix)
      END DO
      DO ik = 1, n
         k = REAL(cell_coordinates(ik, grid), dp)/REAL(grid, dp)
         CALL dbcsr_set(re, 0.0_dp)
         CALL dbcsr_set(im, 0.0_dp)
         IF (antisymmetric) THEN
            ! build_core_ppnl already reverses the SOC integral for swapped atom blocks.
            ! Match the band SOC transform; is_complex would apply that sign a second time.
            CALL rskp_transform(im, re, images, spin, k, kpoints%cell_to_index, pairs)
            CALL dbcsr_desymmetrize(im, full_re)
            CALL dbcsr_desymmetrize(re, full_im)
         ELSE
            CALL rskp_transform(re, im, images, spin, k, kpoints%cell_to_index, pairs)
            CALL dbcsr_desymmetrize(re, full_re)
            CALL dbcsr_desymmetrize(im, full_im)
         END IF
         DO ir = 1, n
            r = cell_coordinates(ir, grid)
            phase = twopi*DOT_PRODUCT(k, REAL(r, dp))
            CALL dbcsr_add(residues(ir)%matrix, full_re, 1.0_dp, COS(phase)/REAL(n, dp))
            CALL dbcsr_add(residues(ir)%matrix, full_im, 1.0_dp, SIN(phase)/REAL(n, dp))
         END DO
      END DO
      DO ir = 1, n
         r = cell_coordinates(ir, grid)
         CALL lift_blocks(residues(ir)%matrix, grid, r, 0, .FALSE., matrix)
      END DO
      CALL dbcsr_finalize(matrix)
      CALL dbcsr_deallocate_matrix_set(residues)
      CALL dbcsr_release(re)
      CALL dbcsr_release(im)
      CALL dbcsr_release(full_re)
      CALL dbcsr_release(full_im)
   END SUBROUTINE lift_fourier

! **************************************************************************************************
!> \brief Lift ordered Berry integrals with the phase of every translated bra cell.
!> \param cosmat Primitive ordered cosine integrals
!> \param sinmat Primitive ordered sine integrals
!> \param cells Integer image labels
!> \param grid Torus dimensions
!> \param axis Direction of the reciprocal torus coordinate
!> \param cosine Full torus cosine
!> \param sine Full torus sine
! **************************************************************************************************
   SUBROUTINE lift_berry(cosmat, sinmat, cells, grid, axis, cosine, sine)
      TYPE(dbcsr_p_type), INTENT(INOUT)                  :: cosmat(:, :), sinmat(:, :)
      INTEGER, INTENT(IN)                                :: cells(:, :), grid(3), axis
      TYPE(dbcsr_type), INTENT(INOUT)                    :: cosine, sine

      INTEGER                                            :: ir

      DO ir = 1, SIZE(cells, 2)
         CALL lift_blocks(cosmat(1, ir)%matrix, grid, cells(:, ir), axis, .FALSE., cosine)
         CALL lift_blocks(sinmat(1, ir)%matrix, grid, cells(:, ir), -axis, .TRUE., cosine)
         CALL lift_blocks(cosmat(1, ir)%matrix, grid, cells(:, ir), axis, .TRUE., sine)
         CALL lift_blocks(sinmat(1, ir)%matrix, grid, cells(:, ir), axis, .FALSE., sine)
      END DO
      CALL dbcsr_finalize(cosine)
      CALL dbcsr_finalize(sine)
   END SUBROUTINE lift_berry

! **************************************************************************************************
!> \brief Translate primitive AO blocks onto the torus without changing MPI ownership.
!> \param source Fully stored primitive matrix for one displacement
!> \param grid Torus dimensions
!> \param displacement Ket-cell displacement relative to the bra
!> \param axis Zero for no phase, signed reciprocal direction otherwise
!> \param sine Use sine rather than cosine of the translation phase
!> \param destination Accumulated torus matrix
! **************************************************************************************************
   SUBROUTINE lift_blocks(source, grid, displacement, axis, sine, destination)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: source
      INTEGER, INTENT(IN)                                :: grid(3), displacement(3), axis
      LOGICAL, INTENT(IN)                                :: sine
      TYPE(dbcsr_type), INTENT(INOUT)                    :: destination

      INTEGER                                            :: blocks, col, it, iu, row, t(3), u(3)
      LOGICAL                                            :: transposed
      REAL(dp)                                           :: factor, phase
      REAL(dp), POINTER                                  :: block(:, :)
      TYPE(dbcsr_iterator_type)                          :: iterator

      CALL dbcsr_get_info(source, nblkrows_total=blocks)
      CALL dbcsr_iterator_start(iterator, source)
      DO WHILE (dbcsr_iterator_blocks_left(iterator))
         CALL dbcsr_iterator_next_block(iterator, row, col, block, transposed=transposed)
         CPASSERT(.NOT. transposed)
         IF (ALL(block == 0.0_dp)) CYCLE
         DO it = 1, PRODUCT(grid)
            t = cell_coordinates(it, grid)
            u = MODULO(t + displacement, grid)
            iu = 1 + u(1) + grid(1)*(u(2) + grid(2)*u(3))
            factor = 1.0_dp
            IF (axis /= 0) THEN
               phase = twopi*REAL(t(ABS(axis)), dp)/REAL(grid(ABS(axis)), dp)
               IF (sine) THEN
                  factor = SIN(phase)
               ELSE
                  factor = COS(phase)
               END IF
               IF (axis < 0) factor = -factor
            END IF
            IF (factor == 0.0_dp) CYCLE
            CALL dbcsr_put_block(destination, row + (it - 1)*blocks, col + (iu - 1)*blocks, &
                                 factor*block, summation=.TRUE.)
         END DO
      END DO
      CALL dbcsr_iterator_stop(iterator)
      CALL dbcsr_finalize(destination)
   END SUBROUTINE lift_blocks

! **************************************************************************************************
!> \brief Zero-based cell coordinates with the first axis varying fastest.
!> \param index One-based cell index
!> \param grid Torus dimensions
!> \return Integer lattice coordinates
! **************************************************************************************************
   PURE FUNCTION cell_coordinates(index, grid) RESULT(r)
      INTEGER, INTENT(IN)                                :: index, grid(3)
      INTEGER                                            :: r(3)

      r(1) = MODULO(index - 1, grid(1))
      r(2) = MODULO((index - 1)/grid(1), grid(2))
      r(3) = (index - 1)/(grid(1)*grid(2))
   END FUNCTION cell_coordinates
END MODULE qs_localizer_torus
