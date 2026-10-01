!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Distributed block-sparse finite localizers, retaining the covariant AO metric.
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
MODULE spectral_localizer_dbcsr
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_create, dbcsr_desymmetrize, dbcsr_distribution_get, &
        dbcsr_distribution_new, dbcsr_distribution_release, dbcsr_distribution_type, &
        dbcsr_finalize, dbcsr_get_info, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_p_type, &
        dbcsr_put_block, dbcsr_release, dbcsr_transposed, dbcsr_type, dbcsr_type_no_symmetry
   USE ieee_arithmetic,                 ONLY: ieee_is_finite
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'spectral_localizer_dbcsr'
   PUBLIC :: localizer_dbcsr_build, localizer_dbcsr_dense, localizer_dbcsr_release, localizer_matrix_type
   PUBLIC :: localizer_dbcsr_check
   PUBLIC :: localizer_ao_operators
   PUBLIC :: localizer_dbcsr_realify, localizer_dbcsr_from_real
   PUBLIC :: localizer_dbcsr_gap_pencil
   PUBLIC :: localizer_dbcsr_3d_build

   TYPE localizer_matrix_type
      TYPE(dbcsr_type) :: real_part, imaginary_part
   END TYPE localizer_matrix_type

CONTAINS

! **************************************************************************************************
!> \brief Build the 3D AII chiral block in DBCSR; not a Hermitian 2D localizer.
!> \param h Scalar AO Hamiltonian
!> \param s Scalar AO overlap
!> \param positions Three Cartesian AO first moments
!> \param point Cartesian position query
!> \param energy Absolute energy query
!> \param kappa Positive energy/length scale
!> \param block Chiral block, with auxiliary index outside physical up/down spin
!> \param metric Block metric I_aux tensor I_spin tensor S
!> \param soc Three real antisymmetric SOC components
!> \param shifted_hamiltonian Optional already energy-shifted physical-spin Hamiltonian
!> \param periodic_cosines Three cos(G.r) matrices. Then positions are sin(G.r), point is
!>        in radians, and kappa is the periodic energy scale eta.
! **************************************************************************************************
   SUBROUTINE localizer_dbcsr_3d_build(h, s, positions, point, energy, kappa, block, metric, soc, &
                                       shifted_hamiltonian, periodic_cosines)
      TYPE(dbcsr_type), INTENT(IN)                       :: h, s
      TYPE(dbcsr_p_type), INTENT(IN)                     :: positions(:)
      REAL(dp), INTENT(IN)                               :: point(3), energy, kappa
      TYPE(localizer_matrix_type), INTENT(OUT)           :: block, metric
      TYPE(dbcsr_p_type), INTENT(IN)                     :: soc(:, :)
      TYPE(localizer_matrix_type), INTENT(IN), OPTIONAL  :: shifted_hamiltonian
      TYPE(dbcsr_p_type), INTENT(IN), OPTIONAL           :: periodic_cosines(:)

      COMPLEX(dp), PARAMETER                             :: ii = CMPLX(0.0_dp, 1.0_dp, dp)

      COMPLEX(dp)                                        :: factors(2, 2), identity(4, 4), &
                                                            pauli(4, 4, 3)
      INTEGER                                            :: i
      TYPE(localizer_matrix_type)                        :: physical(1), physical_metric

      CPASSERT(SIZE(positions) == 3)
      identity = 0.0_dp
      pauli = 0.0_dp
      DO i = 1, 4
         identity(i, i) = 1.0_dp
      END DO
      DO i = 1, 2
         pauli(i, i + 2, 1) = 1.0_dp
         pauli(i + 2, i, 1) = 1.0_dp
         pauli(i, i + 2, 2) = -ii
         pauli(i + 2, i, 2) = ii
         pauli(i, i, 3) = 1.0_dp
         pauli(i + 2, i + 2, 3) = -1.0_dp
      END DO
      CALL create_operator(s, 4, block)
      CALL create_operator(s, 4, metric)
      CALL add_real_blocks(s, identity, metric)
      factors = 0.0_dp
      factors(1, 1) = 1.0_dp
      factors(2, 2) = 1.0_dp
      IF (PRESENT(shifted_hamiltonian)) THEN
         CPASSERT(energy == 0.0_dp)
         CALL add_real_blocks(shifted_hamiltonian%real_part, factors, block)
         CALL add_real_blocks(shifted_hamiltonian%imaginary_part, ii*factors, block)
      ELSE
         CALL localizer_ao_operators(h, s, positions(:0), physical, physical_metric, soc)
         CALL add_real_blocks(physical(1)%real_part, factors, block)
         CALL add_real_blocks(physical(1)%imaginary_part, ii*factors, block)
         CALL localizer_dbcsr_release(physical(1))
         CALL localizer_dbcsr_release(physical_metric)
         CALL add_real_blocks(s, -energy*identity, block)
      END IF
      IF (PRESENT(periodic_cosines)) THEN
         CPASSERT(SIZE(periodic_cosines) == 3)
         CALL add_real_blocks(s, -3.0_dp*kappa*identity, block)
         DO i = 1, 3
            CALL add_real_blocks(periodic_cosines(i)%matrix, &
                                 kappa*(COS(point(i))*identity + ii*SIN(point(i))*pauli(:, :, i)), block)
            CALL add_real_blocks(positions(i)%matrix, &
                                 kappa*(SIN(point(i))*identity - ii*COS(point(i))*pauli(:, :, i)), block)
         END DO
      ELSE
         DO i = 1, 3
            CALL add_real_blocks(positions(i)%matrix, -ii*kappa*pauli(:, :, i), block)
            CALL add_real_blocks(s, ii*kappa*point(i)*pauli(:, :, i), block)
         END DO
      END IF
   END SUBROUTINE localizer_dbcsr_3d_build

! **************************************************************************************************
!> \brief Distributed real matrix R(A)=[Re(A),-Im(A);Im(A),Re(A)].
!> \param matrix Complex operator with full-storage components
!> \param realified Newly initialized real DBCSR operator
! **************************************************************************************************
   SUBROUTINE localizer_dbcsr_realify(matrix, realified)
      TYPE(localizer_matrix_type), INTENT(IN)            :: matrix
      TYPE(dbcsr_type), INTENT(OUT)                      :: realified

      COMPLEX(dp)                                        :: factors(2, 2)
      TYPE(localizer_matrix_type)                        :: expanded

      CALL create_operator(matrix%real_part, 2, expanded)
      factors = RESHAPE([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2])
      CALL add_real_blocks(matrix%real_part, factors, expanded)
      factors = RESHAPE([0.0_dp, 1.0_dp, -1.0_dp, 0.0_dp], [2, 2])
      CALL add_real_blocks(matrix%imaginary_part, factors, expanded)
      CALL dbcsr_create(realified, template=expanded%real_part)
      CALL dbcsr_copy(realified, expanded%real_part)
      CALL localizer_dbcsr_release(expanded)
   END SUBROUTINE localizer_dbcsr_realify

! **************************************************************************************************
!> \brief Extract the two left quadrants of R(A); does not average away structure errors.
!> \param realified Fully stored real operator with two repeated component distributions
!> \param template Original complex operator defining the physical AO block distribution
!> \param matrix Newly initialized complex operator; caller must check realification residual
! **************************************************************************************************
   SUBROUTINE localizer_dbcsr_from_real(realified, template, matrix)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: realified
      TYPE(localizer_matrix_type), INTENT(IN)            :: template
      TYPE(localizer_matrix_type), INTENT(OUT)           :: matrix

      INTEGER                                            :: blocks, col, nblocks, row
      LOGICAL                                            :: transposed
      REAL(dp), POINTER                                  :: block(:, :)
      TYPE(dbcsr_iterator_type)                          :: iterator

      CALL dbcsr_get_info(template%real_part, nblkrows_total=blocks)
      CALL dbcsr_get_info(realified, nblkrows_total=nblocks)
      CPASSERT(nblocks == 2*blocks)
      CALL create_operator(template%real_part, 1, matrix)
      CALL dbcsr_iterator_start(iterator, realified)
      DO WHILE (dbcsr_iterator_blocks_left(iterator))
         CALL dbcsr_iterator_next_block(iterator, row, col, block, transposed=transposed)
         CPASSERT(.NOT. transposed)
         IF (col > blocks) CYCLE
         IF (row <= blocks) THEN
            CALL dbcsr_put_block(matrix%real_part, row, col, block)
         ELSE
            CALL dbcsr_put_block(matrix%imaginary_part, row - blocks, col, block)
         END IF
      END DO
      CALL dbcsr_iterator_stop(iterator)
      CALL dbcsr_finalize(matrix%real_part)
      CALL dbcsr_finalize(matrix%imaginary_part)
   END SUBROUTINE localizer_dbcsr_from_real

! **************************************************************************************************
!> \brief Create a complex block-sparse operator with component-major scalar AO ordering.
!> \param source Scalar AO template
!> \param components Number of spin/auxiliary components
!> \param matrix Newly initialized operator
! **************************************************************************************************
   SUBROUTINE create_operator(source, components, matrix)
      TYPE(dbcsr_type), INTENT(IN)                       :: source
      INTEGER, INTENT(IN)                                :: components
      TYPE(localizer_matrix_type), INTENT(OUT)           :: matrix

      INTEGER                                            :: blocks
      INTEGER, POINTER                                   :: col_dist(:), expanded_cols(:), &
                                                            expanded_rows(:), expanded_sizes(:), &
                                                            row_dist(:), sizes(:)
      TYPE(dbcsr_distribution_type)                      :: distribution, expanded_distribution

      CALL dbcsr_get_info(source, distribution=distribution, row_blk_size=sizes)
      CALL dbcsr_distribution_get(distribution, row_dist=row_dist, col_dist=col_dist)
      blocks = SIZE(sizes)
      ALLOCATE (expanded_sizes(blocks*components), expanded_rows(blocks*components), expanded_cols(blocks*components))
      expanded_sizes = RESHAPE(SPREAD(sizes, 2, components), [blocks*components])
      expanded_rows = RESHAPE(SPREAD(row_dist, 2, components), [blocks*components])
      expanded_cols = RESHAPE(SPREAD(col_dist, 2, components), [blocks*components])
      CALL dbcsr_distribution_new(expanded_distribution, template=distribution, &
                                  row_dist=expanded_rows, col_dist=expanded_cols)
      CALL dbcsr_create(matrix%real_part, 'FINITE AO RE', expanded_distribution, &
                        dbcsr_type_no_symmetry, expanded_sizes, expanded_sizes)
      CALL dbcsr_create(matrix%imaginary_part, template=matrix%real_part, name='FINITE AO IM')
      CALL dbcsr_distribution_release(expanded_distribution)
      DEALLOCATE (expanded_sizes, expanded_rows, expanded_cols)
   END SUBROUTINE create_operator

! **************************************************************************************************
!> \brief Shared AO Hamiltonian, Hermitian scalar observables and optional GTH physical-spin SOC.
!> \param h Scalar Hamiltonian
!> \param s AO metric
!> \param positions Covariant scalar observables, such as first moments or trigonometric coordinates
!> \param operators Newly initialized Hamiltonian followed by the supplied observables
!> \param metric Newly initialized physical-spin metric
!> \param soc Optional antisymmetric real GTH SOC components
! **************************************************************************************************
   SUBROUTINE localizer_ao_operators(h, s, positions, operators, metric, soc)
      TYPE(dbcsr_type), INTENT(IN)                       :: h, s
      TYPE(dbcsr_p_type), INTENT(IN)                     :: positions(:)
      TYPE(localizer_matrix_type), INTENT(OUT)           :: operators(:), metric
      TYPE(dbcsr_p_type), INTENT(IN), OPTIONAL           :: soc(:, :)

      COMPLEX(dp), ALLOCATABLE                           :: factors(:, :), identity(:, :)
      INTEGER                                            :: i, nspin

      CPASSERT(SIZE(operators) == SIZE(positions) + 1)
      nspin = 1
      IF (PRESENT(soc)) nspin = 2
      ALLOCATE (identity(nspin, nspin), factors(nspin, nspin))
      identity = CMPLX(0.0_dp, 0.0_dp, dp)
      DO i = 1, nspin
         identity(i, i) = 1.0_dp
      END DO
      CALL create_operator(s, nspin, metric)
      CALL add_real_blocks(s, identity, metric)
      DO i = 1, SIZE(operators)
         CALL create_operator(s, nspin, operators(i))
      END DO
      CALL add_real_blocks(h, identity, operators(1))
      DO i = 1, SIZE(positions)
         CALL add_real_blocks(positions(i)%matrix, identity, operators(i + 1))
      END DO
      IF (PRESENT(soc)) THEN
         factors = CMPLX(0.0_dp, 0.0_dp, dp)
         factors(1, 2) = CMPLX(0.0_dp, 1.0_dp, dp)
         factors(2, 1) = CMPLX(0.0_dp, 1.0_dp, dp)
         CALL add_real_blocks(soc(1, 1)%matrix, factors, operators(1))
         factors(1, 2) = -1.0_dp
         factors(2, 1) = 1.0_dp
         CALL add_real_blocks(soc(2, 1)%matrix, factors, operators(1))
         factors = CMPLX(0.0_dp, 0.0_dp, dp)
         factors(1, 1) = CMPLX(0.0_dp, 1.0_dp, dp)
         factors(2, 2) = CMPLX(0.0_dp, -1.0_dp, dp)
         CALL add_real_blocks(soc(3, 1)%matrix, factors, operators(1))
      END IF
   END SUBROUTINE localizer_ao_operators

! **************************************************************************************************
!> \brief Pair H-E*S with its negative to reuse the metric-aware localizer gap solvers.
!> \param h Complex physical Hamiltonian, including spin if present
!> \param s Complex physical metric with the same AO block layout
!> \param energy Reference energy
!> \param pencil Newly allocated diag(H-E*S,-H+E*S)
!> \param metric Newly allocated diag(S,S)
! **************************************************************************************************
   SUBROUTINE localizer_dbcsr_gap_pencil(h, s, energy, pencil, metric)
      TYPE(localizer_matrix_type), INTENT(IN)            :: h, s
      REAL(dp), INTENT(IN)                               :: energy
      TYPE(localizer_matrix_type), INTENT(OUT)           :: pencil, metric

      COMPLEX(dp), PARAMETER :: imaginary = CMPLX(0.0_dp, 1.0_dp, dp)

      COMPLEX(dp)                                        :: factors(2, 2)

      CALL create_operator(s%real_part, 2, pencil)
      CALL create_operator(s%real_part, 2, metric)
      factors = RESHAPE([1.0_dp, 0.0_dp, 0.0_dp, -1.0_dp], [2, 2])
      CALL add_real_blocks(h%real_part, factors, pencil)
      CALL add_real_blocks(h%imaginary_part, imaginary*factors, pencil)
      CALL add_real_blocks(s%real_part, -energy*factors, pencil)
      CALL add_real_blocks(s%imaginary_part, -energy*imaginary*factors, pencil)
      factors(2, 2) = 1.0_dp
      CALL add_real_blocks(s%real_part, factors, metric)
      CALL add_real_blocks(s%imaginary_part, imaginary*factors, metric)
   END SUBROUTINE localizer_dbcsr_gap_pencil

! **************************************************************************************************
!> \brief Assemble L and B=diag(S,S), including optional physical spin and SOC, without dense AOs.
!> \param h Real scalar Hamiltonian
!> \param s Real AO overlap
!> \param x First oriented real AO position
!> \param y Second oriented real AO position
!> \param point Two Cartesian coordinates
!> \param energy Reference energy
!> \param kappa Energy/length scale
!> \param localizer Newly created complex DBCSR localizer
!> \param metric Newly created complex DBCSR metric
!> \param soc Real antisymmetric SOC components from V_SOC_xyz_from_pseudopotential
!> \param periodic_cosines Optional two cos(G.r) matrices. Then x,y are sin(G.r), point is
!>        G.r0, and kappa is the energy scale eta (not the finite energy/length scale).
!> \param shifted_hamiltonian Optional replacement for the physical H-E*S including SOC.
!>        Requires energy=0; the scalar h and SOC contributions are then not added again.
!> \note Component-major ordering: auxiliary index, physical spin, then atomic AO blocks.
! **************************************************************************************************
   SUBROUTINE localizer_dbcsr_build(h, s, x, y, point, energy, kappa, localizer, metric, soc, periodic_cosines, &
                                    shifted_hamiltonian)
      TYPE(dbcsr_type), INTENT(IN)                       :: h, s, x, y
      REAL(KIND=dp), INTENT(IN)                          :: point(2), energy, kappa
      TYPE(localizer_matrix_type), INTENT(OUT)           :: localizer, metric
      TYPE(dbcsr_p_type), INTENT(IN), OPTIONAL           :: soc(:, :), periodic_cosines(:)
      TYPE(localizer_matrix_type), INTENT(IN), OPTIONAL  :: shifted_hamiltonian

      COMPLEX(KIND=dp), PARAMETER :: imaginary = CMPLX(0.0_dp, 1.0_dp, dp)

      COMPLEX(dp)                                        :: signs(2, 2)
      COMPLEX(KIND=dp), ALLOCATABLE                      :: factors(:, :), identity(:, :), sx(:, :), &
                                                            sy(:, :), sz(:, :)
      INTEGER                                            :: blocks, component, components, nspin, &
                                                            physical_blocks
      INTEGER, POINTER                                   :: col_dist(:), expanded_cols(:), &
                                                            expanded_rows(:), expanded_sizes(:), &
                                                            row_dist(:), sizes(:)
      TYPE(dbcsr_distribution_type)                      :: distribution, expanded_distribution

      nspin = 1
      IF (PRESENT(soc)) nspin = 2
      components = 2*nspin
      CALL dbcsr_get_info(s, distribution=distribution, row_blk_size=sizes)
      CALL dbcsr_distribution_get(distribution, row_dist=row_dist, col_dist=col_dist)
      blocks = SIZE(sizes)
      ALLOCATE (expanded_sizes(blocks*components), expanded_rows(blocks*components), expanded_cols(blocks*components))
      expanded_sizes = RESHAPE(SPREAD(sizes, 2, components), [blocks*components])
      expanded_rows = RESHAPE(SPREAD(row_dist, 2, components), [blocks*components])
      expanded_cols = RESHAPE(SPREAD(col_dist, 2, components), [blocks*components])
      CALL dbcsr_distribution_new(expanded_distribution, template=distribution, &
                                  row_dist=expanded_rows, col_dist=expanded_cols)
      CALL dbcsr_create(localizer%real_part, 'SPECTRAL LOCALIZER RE', expanded_distribution, &
                        dbcsr_type_no_symmetry, expanded_sizes, expanded_sizes)
      CALL dbcsr_create(localizer%imaginary_part, template=localizer%real_part, name='SPECTRAL LOCALIZER IM')
      CALL dbcsr_create(metric%real_part, template=localizer%real_part, name='LOCALIZER METRIC RE')
      CALL dbcsr_create(metric%imaginary_part, template=localizer%real_part, name='LOCALIZER METRIC IM')
      CALL dbcsr_distribution_release(expanded_distribution)
      DEALLOCATE (expanded_sizes, expanded_rows, expanded_cols)
      ALLOCATE (factors(components, components), identity(components, components), &
                sx(components, components), sy(components, components), sz(components, components))
      identity = 0.0_dp
      sx = 0.0_dp
      sy = 0.0_dp
      sz = 0.0_dp
      DO component = 1, nspin
         identity(component, component) = 1.0_dp
         identity(nspin + component, nspin + component) = 1.0_dp
         sz(component, component) = 1.0_dp
         sz(nspin + component, nspin + component) = -1.0_dp
         sx(component, nspin + component) = 1.0_dp
         sx(nspin + component, component) = 1.0_dp
         sy(component, nspin + component) = -imaginary
         sy(nspin + component, component) = imaginary
      END DO
      CALL add_real_blocks(s, identity, metric)
      IF (PRESENT(shifted_hamiltonian)) THEN
         CPASSERT(energy == 0.0_dp)
         CALL dbcsr_get_info(shifted_hamiltonian%real_part, nblkrows_total=physical_blocks)
         CPASSERT(physical_blocks == nspin*blocks)
         signs = RESHAPE([1.0_dp, 0.0_dp, 0.0_dp, -1.0_dp], [2, 2])
         CALL add_real_blocks(shifted_hamiltonian%real_part, signs, localizer)
         CALL add_real_blocks(shifted_hamiltonian%imaginary_part, imaginary*signs, localizer)
      ELSE
         CALL add_real_blocks(h, sz, localizer)
      END IF
      IF (PRESENT(periodic_cosines)) THEN
         CPASSERT(SIZE(periodic_cosines) == 2)
         CALL add_real_blocks(s, -(energy + 2.0_dp*kappa)*sz, localizer)
         CALL add_real_blocks(periodic_cosines(1)%matrix, kappa*(COS(point(1))*sz - SIN(point(1))*sx), localizer)
         CALL add_real_blocks(periodic_cosines(2)%matrix, kappa*(COS(point(2))*sz - SIN(point(2))*sy), localizer)
         CALL add_real_blocks(x, kappa*(SIN(point(1))*sz + COS(point(1))*sx), localizer)
         CALL add_real_blocks(y, kappa*(SIN(point(2))*sz + COS(point(2))*sy), localizer)
      ELSE
         CALL add_real_blocks(s, -energy*sz - kappa*(point(1)*sx + point(2)*sy), localizer)
         CALL add_real_blocks(x, kappa*sx, localizer)
         CALL add_real_blocks(y, kappa*sy, localizer)
      END IF
      IF (PRESENT(soc) .AND. .NOT. PRESENT(shifted_hamiltonian)) THEN
         ! The real spherical-harmonic SOC matrices carry the missing factor i explicitly.
         factors = 0.0_dp
         factors(1, 2) = imaginary
         factors(2, 1) = imaginary
         factors(3, 4) = -imaginary
         factors(4, 3) = -imaginary
         CALL add_real_blocks(soc(1, 1)%matrix, factors, localizer)
         factors = 0.0_dp
         factors(1, 2) = -1.0_dp
         factors(2, 1) = 1.0_dp
         factors(3, 4) = 1.0_dp
         factors(4, 3) = -1.0_dp
         CALL add_real_blocks(soc(2, 1)%matrix, factors, localizer)
         factors = 0.0_dp
         factors(1, 1) = imaginary
         factors(2, 2) = -imaginary
         factors(3, 3) = -imaginary
         factors(4, 4) = imaginary
         CALL add_real_blocks(soc(3, 1)%matrix, factors, localizer)
      END IF
   END SUBROUTINE localizer_dbcsr_build

! **************************************************************************************************
!> \brief Add a small Kronecker factor times distributed real AO blocks, preserving their union.
!> \param source Real symmetric or antisymmetric AO matrix
!> \param factors Small component matrix
!> \param destination Complex block-sparse destination
! **************************************************************************************************
   SUBROUTINE add_real_blocks(source, factors, destination)
      TYPE(dbcsr_type), INTENT(IN)                       :: source
      COMPLEX(KIND=dp), INTENT(IN)                       :: factors(:, :)
      TYPE(localizer_matrix_type), INTENT(INOUT)         :: destination

      COMPLEX(KIND=dp), ALLOCATABLE                      :: values(:, :)
      INTEGER                                            :: blocks, col, i, j, row
      LOGICAL                                            :: transposed
      REAL(KIND=dp), POINTER                             :: block(:, :)
      TYPE(dbcsr_iterator_type)                          :: iterator
      TYPE(dbcsr_type)                                   :: full

      CALL dbcsr_get_info(source, nblkrows_total=blocks)
      CALL dbcsr_create(full, template=source, matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_desymmetrize(source, full)
      CALL dbcsr_iterator_start(iterator, full)
      DO WHILE (dbcsr_iterator_blocks_left(iterator))
         CALL dbcsr_iterator_next_block(iterator, row, col, block, transposed=transposed)
         CPASSERT(.NOT. transposed)
         ALLOCATE (values(SIZE(block, 1), SIZE(block, 2)))
         DO j = 1, SIZE(factors, 2)
            DO i = 1, SIZE(factors, 1)
               IF (ABS(factors(i, j)) == 0.0_dp) CYCLE
               values = factors(i, j)*block
               IF (REAL(factors(i, j), dp) /= 0.0_dp) THEN
                  CALL dbcsr_put_block(destination%real_part, row + (i - 1)*blocks, col + (j - 1)*blocks, &
                                       REAL(values, dp), summation=.TRUE.)
               END IF
               IF (AIMAG(factors(i, j)) /= 0.0_dp) THEN
                  CALL dbcsr_put_block(destination%imaginary_part, row + (i - 1)*blocks, col + (j - 1)*blocks, &
                                       AIMAG(values), summation=.TRUE.)
               END IF
            END DO
         END DO
         DEALLOCATE (values)
      END DO
      CALL dbcsr_iterator_stop(iterator)
      CALL dbcsr_release(full)
      CALL dbcsr_finalize(destination%real_part)
      CALL dbcsr_finalize(destination%imaginary_part)
   END SUBROUTINE add_real_blocks

! **************************************************************************************************
!> \brief Collect the size-limited dense oracle; sparse solvers must not call this routine.
!> \param matrix Fully stored complex block-sparse matrix
!> \param para_env Matrix communicator
!> \param dense Replicated matrix
! **************************************************************************************************
   SUBROUTINE localizer_dbcsr_dense(matrix, para_env, dense)
      TYPE(localizer_matrix_type), INTENT(INOUT)         :: matrix
      TYPE(mp_para_env_type), POINTER                    :: para_env
      COMPLEX(KIND=dp), ALLOCATABLE, INTENT(OUT)         :: dense(:, :)

      INTEGER                                            :: col, n, row
      LOGICAL                                            :: transposed
      REAL(KIND=dp), POINTER                             :: block(:, :)
      TYPE(dbcsr_iterator_type)                          :: iterator

      CALL dbcsr_get_info(matrix%real_part, nfullrows_total=n)
      ALLOCATE (dense(n, n))
      dense = 0.0_dp
      CALL dbcsr_iterator_start(iterator, matrix%real_part)
      DO WHILE (dbcsr_iterator_blocks_left(iterator))
         CALL dbcsr_iterator_next_block(iterator, block=block, row_offset=row, col_offset=col, &
                                        transposed=transposed)
         CPASSERT(.NOT. transposed)
         dense(row:row + SIZE(block, 1) - 1, col:col + SIZE(block, 2) - 1) = block
      END DO
      CALL dbcsr_iterator_stop(iterator)
      CALL dbcsr_iterator_start(iterator, matrix%imaginary_part)
      DO WHILE (dbcsr_iterator_blocks_left(iterator))
         CALL dbcsr_iterator_next_block(iterator, block=block, row_offset=row, col_offset=col, &
                                        transposed=transposed)
         CPASSERT(.NOT. transposed)
         dense(row:row + SIZE(block, 1) - 1, col:col + SIZE(block, 2) - 1) = &
            dense(row:row + SIZE(block, 1) - 1, col:col + SIZE(block, 2) - 1) + CMPLX(0.0_dp, 1.0_dp, dp)*block
      END DO
      CALL dbcsr_iterator_stop(iterator)
      CALL para_env%sum(dense)
   END SUBROUTINE localizer_dbcsr_dense

! **************************************************************************************************
!> \brief Check both components before a sparse solver discards the upper triangle.
!> \param matrix Fully stored complex matrix
!> \param para_env Matrix communicator
!> \param tolerance Relative Hermiticity tolerance
!> \param valid Whether all entries are finite and the matrix is Hermitian
!> \param relative_error Largest component transpose residual divided by the largest component entry;
!>        HUGE for nonfinite entries or residuals
! **************************************************************************************************
   SUBROUTINE localizer_dbcsr_check(matrix, para_env, tolerance, valid, relative_error)
      TYPE(localizer_matrix_type), INTENT(INOUT) :: matrix
      TYPE(mp_para_env_type), POINTER :: para_env
      REAL(KIND=dp), INTENT(IN) :: tolerance
      LOGICAL, INTENT(OUT) :: valid
      REAL(KIND=dp), INTENT(OUT), OPTIONAL :: relative_error

      LOGICAL :: valid_real, valid_imaginary
      REAL(KIND=dp) :: real_norm, imaginary_norm, real_error, imaginary_error

      CALL check_component(matrix%real_part, -1.0_dp, real_norm, real_error, valid_real)
      CALL check_component(matrix%imaginary_part, 1.0_dp, imaginary_norm, imaginary_error, valid_imaginary)
      valid = valid_real .AND. valid_imaginary
      valid = valid .AND. MAX(real_error, imaginary_error) <= &
         tolerance*MAX(TINY(1.0_dp), real_norm, imaginary_norm)
      IF (PRESENT(relative_error)) THEN
         relative_error = HUGE(1.0_dp)
         IF (valid_real .AND. valid_imaginary) THEN
            relative_error = MAX(real_error, imaginary_error)/MAX(TINY(1.0_dp), real_norm, imaginary_norm)
         END IF
      END IF

   CONTAINS

! **************************************************************************************************
!> \brief Check real symmetry or imaginary antisymmetry using a distributed transpose.
!> \param component Fully stored real matrix
!> \param sign Sign multiplying the transpose
!> \param norm Largest absolute entry
!> \param error Largest transpose residual
!> \param ok Whether all entries and residuals are finite
! **************************************************************************************************
      SUBROUTINE check_component(component, sign, norm, error, ok)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: component
      REAL(KIND=dp), INTENT(IN)                          :: sign
      REAL(KIND=dp), INTENT(OUT)                         :: norm, error
      LOGICAL, INTENT(OUT)                               :: ok

      INTEGER                                            :: invalid
      REAL(KIND=dp), POINTER                             :: block(:, :)
      TYPE(dbcsr_distribution_type)                      :: distribution
      TYPE(dbcsr_iterator_type)                          :: iterator
      TYPE(dbcsr_type)                                   :: difference

         invalid = 0
         norm = 0.0_dp
         error = 0.0_dp
         CALL dbcsr_iterator_start(iterator, component)
         DO WHILE (dbcsr_iterator_blocks_left(iterator))
            CALL dbcsr_iterator_next_block(iterator, block=block)
            IF (.NOT. ALL(ieee_is_finite(block))) THEN
               invalid = 1
            ELSE
               norm = MAX(norm, MAXVAL(ABS(block)))
            END IF
         END DO
         CALL dbcsr_iterator_stop(iterator)
         CALL para_env%max(invalid)
         ok = .FALSE.
         IF (invalid /= 0) RETURN
         CALL para_env%max(norm)
         CALL dbcsr_get_info(component, distribution=distribution)
         CALL dbcsr_transposed(difference, component, shallow_data_copy=.FALSE., use_distribution=distribution)
         CALL dbcsr_add(difference, component, sign, 1.0_dp)
         error = 0.0_dp
         CALL dbcsr_iterator_start(iterator, difference)
         DO WHILE (dbcsr_iterator_blocks_left(iterator))
            CALL dbcsr_iterator_next_block(iterator, block=block)
            IF (.NOT. ALL(ieee_is_finite(block))) THEN
               invalid = 1
            ELSE
               error = MAX(error, MAXVAL(ABS(block)))
            END IF
         END DO
         CALL dbcsr_iterator_stop(iterator)
         CALL dbcsr_release(difference)
         CALL para_env%max(invalid)
         CALL para_env%max(error)
         ok = invalid == 0
      END SUBROUTINE check_component

   END SUBROUTINE localizer_dbcsr_check

! **************************************************************************************************
!> \brief Release both real components of a complex sparse operator.
!> \param matrix Sparse operator
! **************************************************************************************************
   SUBROUTINE localizer_dbcsr_release(matrix)
      TYPE(localizer_matrix_type), INTENT(INOUT)         :: matrix

      CALL dbcsr_release(matrix%real_part)
      CALL dbcsr_release(matrix%imaginary_part)
   END SUBROUTINE localizer_dbcsr_release

END MODULE spectral_localizer_dbcsr
