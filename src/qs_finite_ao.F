!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Shared finite-system AO position and SOC operators for post-SCF analyses.
!> \author Thomas D. Kuehne, tkuehne@cp2k.org
! **************************************************************************************************
MODULE qs_finite_ao
   USE ai_contraction,                  ONLY: block_add,&
                                              contraction
   USE ai_gaussian_window,              ONLY: gaussian_window_overlap
   USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_create, dbcsr_distribution_get, dbcsr_distribution_type, dbcsr_finalize, &
        dbcsr_get_info, dbcsr_p_type, dbcsr_put_block, dbcsr_set, dbcsr_type_no_symmetry
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE input_constants,                 ONLY: do_method_gapw,&
                                              do_method_gpw
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_para_env_type
   USE orbital_pointers,                ONLY: ncoset
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_integral_utils,               ONLY: basis_set_list_setup,&
                                              get_memory_usage
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_moments,                      ONLY: build_local_moment_matrix
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE soc_pseudopotential_methods,     ONLY: V_SOC_xyz_from_pseudopotential
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_finite_ao'
   PUBLIC :: finite_ao_positions, finite_ao_soc, finite_ao_window
CONTAINS

! **************************************************************************************************
!> \brief Analytic covariant first moments in the coordinate frame of SUBSYS/COORD.
!> \param qs_env Finite Quickstep environment
!> \param positions Allocated X, Y, Z matrices, released by caller
! **************************************************************************************************
   SUBROUTINE finite_ao_positions(qs_env, positions)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), POINTER                        :: positions(:)

      INTEGER                                            :: axis
      LOGICAL                                            :: do_kpoints
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), POINTER                        :: matrix_s(:)
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(neighbor_list_set_p_type), POINTER            :: sab_orb(:)

      CALL get_qs_env(qs_env, cell=cell, dft_control=dft_control, do_kpoints=do_kpoints, &
                      matrix_s=matrix_s, sab_orb=sab_orb)
      IF (ANY(cell%perd /= 0) .OR. do_kpoints) THEN
         CPABORT('Finite AO positions require PERIODIC NONE without KPOINTS.')
      END IF
      SELECT CASE (dft_control%qs_control%method_id)
      CASE (do_method_gpw, do_method_gapw)
      CASE DEFAULT
         CPABORT('Finite AO positions currently support GPW/GAPW only.')
      END SELECT
      ALLOCATE (positions(3))
      DO axis = 1, 3
         ALLOCATE (positions(axis)%matrix)
         CALL dbcsr_create(positions(axis)%matrix, template=matrix_s(1)%matrix, name='AO POSITION')
         CALL cp_dbcsr_alloc_block_from_nbl(positions(axis)%matrix, sab_orb)
         CALL dbcsr_set(positions(axis)%matrix, 0.0_dp)
      END DO
      CALL build_local_moment_matrix(qs_env, positions, 1, ref_point=[0.0_dp, 0.0_dp, 0.0_dp])
   END SUBROUTINE finite_ao_positions
! **************************************************************************************************
!> \brief Covariant matrix of a physical Cartesian Gaussian window in finite AO space.
!> \param qs_env Isolated Quickstep environment
!> \param center Cartesian window center in internal bohr
!> \param alpha Window exponents in inverse squared bohr
!> \param window Newly allocated scalar matrix, released by caller
!> \note Integrates the window over the AOs; no atom or AO-index weighting is used.
! **************************************************************************************************
   SUBROUTINE finite_ao_window(qs_env, center, alpha, window)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      REAL(dp), INTENT(IN)                               :: center(3), alpha(3)
      TYPE(dbcsr_p_type), INTENT(OUT)                    :: window

      INTEGER                                            :: ia, ib, iset, jset, ka, kb, ld, na, nb, &
                                                            sa, sb
      INTEGER, POINTER                                   :: col_dist(:), pgrid(:, :), row_dist(:)
      LOGICAL                                            :: do_kpoints
      REAL(dp)                                           :: distance
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
      IF (ANY(cell%perd /= 0) .OR. do_kpoints) CPABORT('Finite AO windows require isolated Gamma systems.')
      ALLOCATE (bases(SIZE(kinds)))
      CALL basis_set_list_setup(bases, 'ORB', kinds)
      ld = get_memory_usage(kinds, 'ORB')
      ALLOCATE (integrals(ld, ld), work(ld, ld), window%matrix)
      CALL dbcsr_create(window%matrix, template=matrix_s(1, 1)%matrix, &
                        name='AO GAUSSIAN WINDOW', matrix_type=dbcsr_type_no_symmetry)
      CALL dbcsr_get_info(matrix_s(1, 1)%matrix, distribution=distribution)
      CALL dbcsr_distribution_get(distribution, pgrid=pgrid, row_dist=row_dist, col_dist=col_dist)
      DO ia = 1, SIZE(particles)
         ka = particles(ia)%atomic_kind%kind_number
         a => bases(ka)%gto_basis_set
         IF (.NOT. ASSOCIATED(a)) CYCLE
         DO ib = 1, SIZE(particles)
            IF (pgrid(row_dist(ia), col_dist(ib)) /= para_env%mepos) CYCLE
            kb = particles(ib)%atomic_kind%kind_number
            b => bases(kb)%gto_basis_set
            IF (.NOT. ASSOCIATED(b)) CYCLE
            distance = NORM2(particles(ia)%r - particles(ib)%r)
            ALLOCATE (block(a%nsgf, b%nsgf))
            block = 0.0_dp
            DO iset = 1, a%nset
               sa = a%first_sgf(1, iset)
               na = a%npgf(iset)*(ncoset(a%lmax(iset)) - ncoset(a%lmin(iset) - 1))
               DO jset = 1, b%nset
                  IF (a%set_radius(iset) + b%set_radius(jset) < distance) CYCLE
                  sb = b%first_sgf(1, jset)
                  nb = b%npgf(jset)*(ncoset(b%lmax(jset)) - ncoset(b%lmin(jset) - 1))
                  CALL gaussian_window_overlap(a%lmax(iset), a%lmin(iset), a%zet(:a%npgf(iset), iset), &
                                               b%lmax(jset), b%lmin(jset), b%zet(:b%npgf(jset), jset), &
                                               particles(ia)%r, particles(ib)%r, center, alpha, integrals)
                  CALL contraction(integrals, work, ca=a%scon(:, sa:), na=na, ma=a%nsgf_set(iset), &
                                   cb=b%scon(:, sb:), nb=nb, mb=b%nsgf_set(jset), fscale=1.0_dp, trans=.FALSE.)
                  CALL block_add('IN', work, a%nsgf_set(iset), b%nsgf_set(jset), block, sa, sb, trans=.FALSE.)
               END DO
            END DO
            IF (MAXVAL(ABS(block)) > 0.0_dp) CALL dbcsr_put_block(window%matrix, ia, ib, block)
            DEALLOCATE (block)
         END DO
      END DO
      CALL dbcsr_finalize(window%matrix)
   END SUBROUTINE finite_ao_window
! **************************************************************************************************
!> \brief Validate and construct post-SCF GTH SOC components shared by finite AO analyses.
!> \param qs_env Restricted SCF environment with GTH SOC potentials
!> \param soc_matrices Newly allocated SOC matrices, released by caller
! **************************************************************************************************
   SUBROUTINE finite_ao_soc(qs_env, soc_matrices)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(dbcsr_p_type), POINTER                        :: soc_matrices(:, :)

      INTEGER                                            :: ik
      TYPE(dbcsr_p_type), POINTER                        :: matrix_ks(:, :)
      TYPE(neighbor_list_set_p_type), POINTER            :: sap_ppnl(:)
      TYPE(qs_kind_type), POINTER                        :: qs_kind_set(:)

      CALL get_qs_env(qs_env, matrix_ks_kp=matrix_ks, qs_kind_set=qs_kind_set, sap_ppnl=sap_ppnl)
      IF (SIZE(matrix_ks, 1) /= 1) CPABORT('AO SOC requires a restricted SCF calculation.')
      IF (.NOT. ASSOCIATED(sap_ppnl)) CPABORT('Finite AO SOC requires nonlocal GTH SOC projectors.')
      DO ik = 1, SIZE(qs_kind_set)
         IF (.NOT. ASSOCIATED(qs_kind_set(ik)%gth_potential)) THEN
            CPABORT('Finite AO SOC requires GTH potentials.')
         END IF
         IF (.NOT. qs_kind_set(ik)%gth_potential%soc) CPABORT('GTH potential lacks SOC parameters.')
      END DO
      CALL V_SOC_xyz_from_pseudopotential(qs_env, soc_matrices)
   END SUBROUTINE finite_ao_soc
END MODULE qs_finite_ao
