!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Analytic translations of finite or periodized Gaussian AO functions.
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
MODULE qs_ao_translations
   USE ai_contraction,                  ONLY: block_add,&
                                              contraction
   USE ai_overlap,                      ONLY: overlap_ab
   USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_create,&
                                              dbcsr_distribution_get,&
                                              dbcsr_distribution_type,&
                                              dbcsr_finalize,&
                                              dbcsr_get_info,&
                                              dbcsr_p_type,&
                                              dbcsr_put_block,&
                                              dbcsr_type_no_symmetry
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE orbital_pointers,                ONLY: ncoset
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_integral_utils,               ONLY: basis_set_list_setup,&
                                              get_memory_usage
   USE qs_kind_types,                   ONLY: qs_kind_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_ao_translations'
   PUBLIC :: ao_translation_matrices
CONTAINS
! **************************************************************************************************
!> \brief Overlaps with translated ket Gaussians, U(a) psi(r) = psi(r+a).
!> \param qs_env Quickstep environment supplying primitive AOs and their distribution
!> \param shifts Cartesian translation vectors in internal bohr
!> \param translations Newly allocated nonsymmetric primitive-size matrices, released by caller
!> \param periodic_grid If present, sum images of the grid-expanded periodic cell; otherwise finite
!> \note Arbitrary shifts use their own overlap screening, not the SCF neighbour list.
! **************************************************************************************************
   SUBROUTINE ao_translation_matrices(qs_env, shifts, translations, periodic_grid)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(dp), INTENT(IN)                               :: shifts(:, :)
      TYPE(dbcsr_p_type), POINTER                        :: translations(:)
      INTEGER, INTENT(IN), OPTIONAL                      :: periodic_grid(3)

      INTEGER                                            :: axis, ia, ib, iset, it, jset, ka, kb, &
                                                            ld, lower(3), na, nb, nx, ny, nz, sa, &
                                                            sb, upper(3)
      INTEGER, POINTER                                   :: col_dist(:), pgrid(:, :), row_dist(:)
      LOGICAL                                            :: do_kpoints
      REAL(dp)                                           :: distance, h(3, 3), hinv(3, 3), &
                                                            offset(3), rab(3), radius, reach(3), &
                                                            scaled(3)
      REAL(dp), ALLOCATABLE                              :: block(:, :), integrals(:, :), work(:, :)
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_distribution_type)                      :: distribution
      TYPE(dbcsr_p_type), POINTER                        :: matrix_s(:, :)
      TYPE(gto_basis_set_p_type), ALLOCATABLE            :: bases(:)
      TYPE(gto_basis_set_type), POINTER                  :: a, b
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), POINTER                       :: particles(:)
      TYPE(qs_kind_type), POINTER                        :: kinds(:)

      CALL get_qs_env(qs_env, cell=cell, do_kpoints=do_kpoints, matrix_s_kp=matrix_s, &
                      para_env=para_env, particle_set=particles, qs_kind_set=kinds)
      h = 0.0_dp
      hinv = 0.0_dp
      reach = 0.0_dp
      IF (PRESENT(periodic_grid)) THEN
         CPASSERT(ALL(periodic_grid >= 1))
         CPASSERT(ALL(periodic_grid == 1 .OR. cell%perd /= 0))
         DO axis = 1, 3
            h(:, axis) = cell%hmat(:, axis)*periodic_grid(axis)
            hinv(axis, :) = cell%h_inv(axis, :)/periodic_grid(axis)
            IF (cell%perd(axis) /= 0) reach(axis) = SQRT(SUM(hinv(axis, :)**2))
         END DO
      ELSE IF (ANY(cell%perd /= 0) .OR. do_kpoints) THEN
         CPABORT('Finite AO translations require isolated Gamma systems.')
      END IF
      CPASSERT(SIZE(shifts, 1) == 3)
      ALLOCATE (bases(SIZE(kinds)))
      CALL basis_set_list_setup(bases, 'ORB', kinds)
      ld = get_memory_usage(kinds, 'ORB')
      ALLOCATE (integrals(ld, ld), work(ld, ld), translations(SIZE(shifts, 2)))
      CALL dbcsr_get_info(matrix_s(1, 1)%matrix, distribution=distribution)
      CALL dbcsr_distribution_get(distribution, pgrid=pgrid, row_dist=row_dist, col_dist=col_dist)
      DO it = 1, SIZE(shifts, 2)
         ALLOCATE (translations(it)%matrix)
         CALL dbcsr_create(translations(it)%matrix, template=matrix_s(1, 1)%matrix, &
                           name='AO TRANSLATION', matrix_type=dbcsr_type_no_symmetry)
         DO ia = 1, SIZE(particles)
            ka = particles(ia)%atomic_kind%kind_number
            a => bases(ka)%gto_basis_set
            IF (.NOT. ASSOCIATED(a)) CYCLE
            DO ib = 1, SIZE(particles)
               IF (pgrid(row_dist(ia), col_dist(ib)) /= para_env%mepos) CYCLE
               kb = particles(ib)%atomic_kind%kind_number
               b => bases(kb)%gto_basis_set
               IF (.NOT. ASSOCIATED(b)) CYCLE
               offset = particles(ib)%r - particles(ia)%r - shifts(:, it)
               lower = 0
               upper = 0
               IF (PRESENT(periodic_grid)) THEN
                  ! Reciprocal-plane distances bound every image, including skew cells.
                  scaled = MATMUL(hinv, offset)
                  radius = MAXVAL(a%set_radius) + MAXVAL(b%set_radius)
                  DO axis = 1, 3
                     IF (cell%perd(axis) == 0) CYCLE
                     lower(axis) = CEILING(-scaled(axis) - radius*reach(axis))
                     upper(axis) = FLOOR(-scaled(axis) + radius*reach(axis))
                  END DO
               END IF
               ALLOCATE (block(a%nsgf, b%nsgf))
               block = 0.0_dp
               DO nz = lower(3), upper(3)
                  DO ny = lower(2), upper(2)
                     DO nx = lower(1), upper(1)
                        rab = offset + MATMUL(h, REAL([nx, ny, nz], dp))
                        distance = SQRT(SUM(rab**2))
                        DO iset = 1, a%nset
                           sa = a%first_sgf(1, iset)
                           na = a%npgf(iset)*(ncoset(a%lmax(iset)) - ncoset(a%lmin(iset) - 1))
                           DO jset = 1, b%nset
                              IF (a%set_radius(iset) + b%set_radius(jset) < distance) CYCLE
                              sb = b%first_sgf(1, jset)
                              nb = b%npgf(jset)*(ncoset(b%lmax(jset)) - ncoset(b%lmin(jset) - 1))
                              CALL overlap_ab(a%lmax(iset), a%lmin(iset), a%npgf(iset), &
                                              a%pgf_radius(:, iset), a%zet(:, iset), &
                                              b%lmax(jset), b%lmin(jset), b%npgf(jset), &
                                              b%pgf_radius(:, jset), b%zet(:, jset), rab, sab=integrals)
                              CALL contraction(integrals, work, ca=a%scon(:, sa:), na=na, ma=a%nsgf_set(iset), &
                                               cb=b%scon(:, sb:), nb=nb, mb=b%nsgf_set(jset), &
                                               fscale=1.0_dp, trans=.FALSE.)
                              CALL block_add('IN', work, a%nsgf_set(iset), b%nsgf_set(jset), &
                                             block, sa, sb, trans=.FALSE.)
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
               IF (MAXVAL(ABS(block)) > 0.0_dp) CALL dbcsr_put_block(translations(it)%matrix, ia, ib, block)
               DEALLOCATE (block)
            END DO
         END DO
         CALL dbcsr_finalize(translations(it)%matrix)
      END DO
   END SUBROUTINE ao_translation_matrices
END MODULE qs_ao_translations
