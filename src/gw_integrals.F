!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Utility method to build 3-center integrals for small cell GW
! **************************************************************************************************
MODULE gw_integrals
   USE OMP_LIB,                         ONLY: omp_get_thread_num
   USE ai_contraction_sphi,             ONLY: abc_contract_xsmm
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              get_cell,&
                                              pbc
   USE cp_array_utils,                  ONLY: cp_2d_r_p_type
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE gamma,                           ONLY: init_md_ftable
   USE input_constants,                 ONLY: do_potential_coulomb,&
                                              do_potential_id,&
                                              do_potential_short,&
                                              do_potential_truncated
   USE kinds,                           ONLY: dp
   USE libint_2c_3c,                    ONLY: cutoff_screen_factor,&
                                              eri_3center,&
                                              libint_potential_type
   USE libint_wrapper,                  ONLY: cp_libint_cleanup_3eri,&
                                              cp_libint_init_3eri,&
                                              cp_libint_set_contrdepth,&
                                              cp_libint_t
   USE message_passing,                 ONLY: mp_para_env_type
   USE orbital_pointers,                ONLY: ncoset
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE t_c_g0,                          ONLY: get_lmax_init,&
                                              init

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_integrals'

   PUBLIC :: build_3c_integral_block, build_3c_integral_block_ctx, &
             gw_3c_ctx_type, gw_3c_ctx_create, gw_3c_ctx_release, &
             gw_3c_ws_type, gw_3c_ws_create, gw_3c_ws_release

! **************************************************************************************************
!> \brief Shared read-only context for repeated 3-center integral block builds: screening
!>        parameters, basis maxima, contracted sphi tables, and the one-time gamma /
!>        truncated-Coulomb table initializations. Create and release OUTSIDE any OMP parallel
!>        region; creation is MPI-collective when the potential is truncated.
! **************************************************************************************************
   TYPE gw_3c_ctx_type
      TYPE(libint_potential_type)                        :: potential_parameter = libint_potential_type()
      INTEGER                                            :: op_ij = do_potential_id, &
                                                            op_jk = do_potential_id
      REAL(KIND=dp)                                      :: dr_ij = 0.0_dp, dr_jk = 0.0_dp, &
                                                            dr_ik = 0.0_dp
      INTEGER                                            :: maxli = 0, maxlj = 0, maxlk = 0, &
                                                            max_am = 0, m_max = 0
      INTEGER                                            :: max_ncoi = 0, max_ncoj = 0, max_ncok = 0
      INTEGER                                            :: max_nsgfi = 0, max_nsgfj = 0, &
                                                            max_nsgfk = 0, max_nset = 0, natom = 0
      TYPE(cp_2d_r_p_type), DIMENSION(:, :), POINTER     :: spi => NULL(), tspj => NULL(), &
                                                            spk => NULL()
      TYPE(gto_basis_set_p_type), DIMENSION(:), ALLOCATABLE :: basis_i, basis_j, basis_k
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set => NULL()
      TYPE(cell_type), POINTER                           :: cell => NULL()
      REAL(KIND=dp), DIMENSION(3, 3)                     :: hmat = 0.0_dp
   END TYPE gw_3c_ctx_type

! **************************************************************************************************
!> \brief Per-thread workspace for 3-center integral block builds: libint object + contraction
!>        buffers. Each thread creates its own (inside the parallel region is fine).
! **************************************************************************************************
   TYPE gw_3c_ws_type
      TYPE(cp_libint_t)                                  :: lib
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: cpp_buffer, ccp_buffer
   END TYPE gw_3c_ws_type

CONTAINS

! **************************************************************************************************
!> \brief Builds the shared 3c-integral context: screening radii, basis maxima, contracted sphi
!>        tables, and the one-time gamma / truncated-Coulomb table initializations.
!> \param ctx ...
!> \param qs_env ...
!> \param potential_parameter ...
!> \param basis_j ...
!> \param basis_k ...
!> \param basis_i ...
! **************************************************************************************************
   SUBROUTINE gw_3c_ctx_create(ctx, qs_env, potential_parameter, basis_j, basis_k, basis_i)

      TYPE(gw_3c_ctx_type), INTENT(OUT)                  :: ctx
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(libint_potential_type), INTENT(IN)            :: potential_parameter
      TYPE(gto_basis_set_p_type), DIMENSION(:)           :: basis_j, basis_k, basis_i

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'gw_3c_ctx_create'

      INTEGER                                            :: egfi, handle, ibasis, ilist, imax, iset, &
                                                            jset, kset, nbasis, ncoi, sgfi, unit_id
      INTEGER, DIMENSION(:), POINTER                     :: npgfi, npgfj, npgfk, nsgfi, nsgfj, nsgfk
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(gto_basis_set_type), POINTER                  :: basis_set
      TYPE(mp_para_env_type), POINTER                    :: para_env

      CALL timeset(routineN, handle)

      ctx%potential_parameter = potential_parameter
      ctx%op_ij = potential_parameter%potential_type
      ctx%op_jk = do_potential_id

      IF (ctx%op_ij == do_potential_truncated .OR. ctx%op_ij == do_potential_short) THEN
         ctx%dr_ij = potential_parameter%cutoff_radius*cutoff_screen_factor
         ctx%dr_ik = potential_parameter%cutoff_radius*cutoff_screen_factor
      ELSE IF (ctx%op_ij == do_potential_coulomb) THEN
         ctx%dr_ij = 1000000.0_dp
         ctx%dr_ik = 1000000.0_dp
      END IF

      NULLIFY (atomic_kind_set, para_env)
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, cell=ctx%cell, natom=ctx%natom, &
                      para_env=para_env, particle_set=ctx%particle_set)
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, kind_of=ctx%kind_of)
      CALL get_cell(cell=ctx%cell, h=ctx%hmat)

      ctx%basis_i = basis_i
      ctx%basis_j = basis_j
      ctx%basis_k = basis_k

      ! max l per basis for libint; max nset/nco/nsgf for the LIBXSMM contraction buffers
      nbasis = SIZE(basis_i)
      DO ibasis = 1, nbasis
         CALL get_gto_basis_set(gto_basis_set=basis_i(ibasis)%gto_basis_set, maxl=imax, &
                                nset=iset, nsgf_set=nsgfi, npgf=npgfi)
         ctx%maxli = MAX(ctx%maxli, imax)
         ctx%max_nset = MAX(ctx%max_nset, iset)
         ctx%max_nsgfi = MAX(ctx%max_nsgfi, MAXVAL(nsgfi))
         ctx%max_ncoi = MAX(ctx%max_ncoi, MAXVAL(npgfi)*ncoset(ctx%maxli))
      END DO
      DO ibasis = 1, nbasis
         CALL get_gto_basis_set(gto_basis_set=basis_j(ibasis)%gto_basis_set, maxl=imax, &
                                nset=jset, nsgf_set=nsgfj, npgf=npgfj)
         ctx%maxlj = MAX(ctx%maxlj, imax)
         ctx%max_nset = MAX(ctx%max_nset, jset)
         ctx%max_nsgfj = MAX(ctx%max_nsgfj, MAXVAL(nsgfj))
         ctx%max_ncoj = MAX(ctx%max_ncoj, MAXVAL(npgfj)*ncoset(ctx%maxlj))
      END DO
      DO ibasis = 1, nbasis
         CALL get_gto_basis_set(gto_basis_set=basis_k(ibasis)%gto_basis_set, maxl=imax, &
                                nset=kset, nsgf_set=nsgfk, npgf=npgfk)
         ctx%maxlk = MAX(ctx%maxlk, imax)
         ctx%max_nset = MAX(ctx%max_nset, kset)
         ctx%max_nsgfk = MAX(ctx%max_nsgfk, MAXVAL(nsgfk))
         ctx%max_ncok = MAX(ctx%max_ncok, MAXVAL(npgfk)*ncoset(ctx%maxlk))
      END DO
      ctx%m_max = ctx%maxli + ctx%maxlj + ctx%maxlk
      ctx%max_am = MAX(ctx%maxli, ctx%maxlj, ctx%maxlk)

      ! contiguous (and for j transposed) sphi copies, shared read-only across threads
      ALLOCATE (ctx%spi(ctx%max_nset, nbasis), ctx%tspj(ctx%max_nset, nbasis), &
                ctx%spk(ctx%max_nset, nbasis))
      DO ibasis = 1, nbasis
         DO iset = 1, ctx%max_nset
            NULLIFY (ctx%spi(iset, ibasis)%array)
            NULLIFY (ctx%tspj(iset, ibasis)%array)
            NULLIFY (ctx%spk(iset, ibasis)%array)
         END DO
      END DO
      DO ilist = 1, 3
         DO ibasis = 1, nbasis
            IF (ilist == 1) basis_set => basis_i(ibasis)%gto_basis_set
            IF (ilist == 2) basis_set => basis_j(ibasis)%gto_basis_set
            IF (ilist == 3) basis_set => basis_k(ibasis)%gto_basis_set
            DO iset = 1, basis_set%nset
               ncoi = basis_set%npgf(iset)*ncoset(basis_set%lmax(iset))
               sgfi = basis_set%first_sgf(1, iset)
               egfi = sgfi + basis_set%nsgf_set(iset) - 1
               IF (ilist == 1) THEN
                  ALLOCATE (ctx%spi(iset, ibasis)%array(ncoi, basis_set%nsgf_set(iset)))
                  ctx%spi(iset, ibasis)%array(:, :) = basis_set%sphi(1:ncoi, sgfi:egfi)
               ELSE IF (ilist == 2) THEN
                  ALLOCATE (ctx%tspj(iset, ibasis)%array(basis_set%nsgf_set(iset), ncoi))
                  ctx%tspj(iset, ibasis)%array(:, :) = TRANSPOSE(basis_set%sphi(1:ncoi, sgfi:egfi))
               ELSE
                  ALLOCATE (ctx%spk(iset, ibasis)%array(ncoi, basis_set%nsgf_set(iset)))
                  ctx%spk(iset, ibasis)%array(:, :) = basis_set%sphi(1:ncoi, sgfi:egfi)
               END IF
            END DO
         END DO
      END DO

      ! one-time table inits; the truncated-Coulomb init reads a file + bcasts => MPI-collective,
      ! must happen here and never inside the per-block path or an OMP region
      IF (ctx%op_ij == do_potential_truncated .OR. ctx%op_jk == do_potential_truncated) THEN
         IF (ctx%m_max > get_lmax_init()) THEN
            IF (para_env%mepos == 0) THEN
               CALL open_file(unit_number=unit_id, file_name=potential_parameter%filename)
            END IF
            CALL init(ctx%m_max, unit_id, para_env%mepos, para_env)
            IF (para_env%mepos == 0) THEN
               CALL close_file(unit_id)
            END IF
         END IF
      END IF
      CALL init_md_ftable(nmax=ctx%m_max)

      CALL timestop(handle)

   END SUBROUTINE gw_3c_ctx_create

! **************************************************************************************************
!> \brief Releases the shared 3c-integral context.
!> \param ctx ...
! **************************************************************************************************
   SUBROUTINE gw_3c_ctx_release(ctx)

      TYPE(gw_3c_ctx_type), INTENT(INOUT)                :: ctx

      INTEGER                                            :: ibasis, iset

      DO iset = 1, SIZE(ctx%spi, 1)
         DO ibasis = 1, SIZE(ctx%spi, 2)
            IF (ASSOCIATED(ctx%spi(iset, ibasis)%array)) DEALLOCATE (ctx%spi(iset, ibasis)%array)
            IF (ASSOCIATED(ctx%tspj(iset, ibasis)%array)) DEALLOCATE (ctx%tspj(iset, ibasis)%array)
            IF (ASSOCIATED(ctx%spk(iset, ibasis)%array)) DEALLOCATE (ctx%spk(iset, ibasis)%array)
         END DO
      END DO
      DEALLOCATE (ctx%spi, ctx%tspj, ctx%spk)
      NULLIFY (ctx%spi, ctx%tspj, ctx%spk, ctx%particle_set, ctx%cell)
      IF (ALLOCATED(ctx%kind_of)) DEALLOCATE (ctx%kind_of)
      IF (ALLOCATED(ctx%basis_i)) DEALLOCATE (ctx%basis_i)
      IF (ALLOCATED(ctx%basis_j)) DEALLOCATE (ctx%basis_j)
      IF (ALLOCATED(ctx%basis_k)) DEALLOCATE (ctx%basis_k)

   END SUBROUTINE gw_3c_ctx_release

! **************************************************************************************************
!> \brief Creates a per-thread 3c workspace: libint object + LIBXSMM contraction buffers.
!> \param ws ...
!> \param ctx ...
! **************************************************************************************************
   SUBROUTINE gw_3c_ws_create(ws, ctx)

      TYPE(gw_3c_ws_type), INTENT(OUT)                   :: ws
      TYPE(gw_3c_ctx_type), INTENT(IN)                   :: ctx

      CALL cp_libint_init_3eri(ws%lib, ctx%max_am)
      CALL cp_libint_set_contrdepth(ws%lib, 1)
      ALLOCATE (ws%cpp_buffer(ctx%max_nsgfj*ctx%max_ncok), &
                ws%ccp_buffer(ctx%max_nsgfj*ctx%max_nsgfk*ctx%max_ncoi))

   END SUBROUTINE gw_3c_ws_create

! **************************************************************************************************
!> \brief Releases a per-thread 3c workspace.
!> \param ws ...
! **************************************************************************************************
   SUBROUTINE gw_3c_ws_release(ws)

      TYPE(gw_3c_ws_type), INTENT(INOUT)                 :: ws

      CALL cp_libint_cleanup_3eri(ws%lib)
      DEALLOCATE (ws%cpp_buffer, ws%ccp_buffer)

   END SUBROUTINE gw_3c_ws_release

! **************************************************************************************************
!> \brief Computes the 3c integral block (mu(atom_j) nu(atom_k) | P(atom_i)) for ONE atom triple,
!>        accumulating into the caller-zeroed int_3c at the given block offsets.
!>        Thread-safe: reads the frozen ctx + read-only module tables, mutates only its arguments
!>        and the per-thread ws.
!> \param int_3c pre-zeroed target; the triple's contribution is accumulated in place
!> \param ctx shared context from gw_3c_ctx_create
!> \param ws per-thread workspace from gw_3c_ws_create
!> \param atom_j ...
!> \param atom_k ...
!> \param atom_i ...
!> \param cell_j ...
!> \param cell_k ...
!> \param cell_i ...
!> \param j_offset block offset of atom_j's first sgf in int_3c dim 1 (default 0)
!> \param k_offset block offset of atom_k's first sgf in int_3c dim 2 (default 0)
!> \param i_offset block offset of atom_i's first RI sgf in int_3c dim 3 (default 0)
!> \param screened .TRUE. if the kind-radius screens killed the whole triple (int_3c untouched)
! **************************************************************************************************
   SUBROUTINE build_3c_integral_block_ctx(int_3c, ctx, ws, atom_j, atom_k, atom_i, &
                                          cell_j, cell_k, cell_i, &
                                          j_offset, k_offset, i_offset, screened)

      REAL(KIND=dp), DIMENSION(:, :, :)                  :: int_3c
      TYPE(gw_3c_ctx_type), INTENT(IN)                   :: ctx
      TYPE(gw_3c_ws_type), INTENT(INOUT)                 :: ws
      INTEGER, INTENT(IN)                                :: atom_j, atom_k, atom_i
      INTEGER, DIMENSION(3), INTENT(IN), OPTIONAL        :: cell_j, cell_k, cell_i
      INTEGER, INTENT(IN), OPTIONAL                      :: j_offset, k_offset, i_offset
      LOGICAL, INTENT(OUT), OPTIONAL                     :: screened

      INTEGER :: block_end_i, block_end_j, block_end_k, block_start_i, block_start_j, &
         block_start_k, ikind, iset, jkind, jset, kkind, kset, my_i_offset, my_j_offset, &
         my_k_offset, ncoi, ncoj, ncok, nseti, nsetj, nsetk, sgfi, sgfj, sgfk
      INTEGER, DIMENSION(3)                              :: my_cell_i, my_cell_j, my_cell_k
      INTEGER, DIMENSION(:), POINTER                     :: lmax_i, lmax_j, lmax_k, lmin_i, lmin_j, &
                                                            lmin_k, npgfi, npgfj, npgfk, nsgfi, &
                                                            nsgfj, nsgfk
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgf_i, first_sgf_j, first_sgf_k
      REAL(KIND=dp)                                      :: dij, dik, djk, kind_radius_i, &
                                                            kind_radius_j, kind_radius_k, sijk_ext
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: sijk, sijk_contr
      REAL(KIND=dp), DIMENSION(3)                        :: ri, rij, rik, rj, rjk, rk
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_i, set_radius_j, set_radius_k
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rpgf_i, rpgf_j, rpgf_k, zeti, zetj, zetk

      IF (PRESENT(screened)) screened = .FALSE.

      my_cell_i(1:3) = 0
      IF (PRESENT(cell_i)) my_cell_i(1:3) = cell_i(1:3)
      my_cell_j(1:3) = 0
      IF (PRESENT(cell_j)) my_cell_j(1:3) = cell_j(1:3)
      my_cell_k(1:3) = 0
      IF (PRESENT(cell_k)) my_cell_k(1:3) = cell_k(1:3)
      my_i_offset = 0
      IF (PRESENT(i_offset)) my_i_offset = i_offset
      my_j_offset = 0
      IF (PRESENT(j_offset)) my_j_offset = j_offset
      my_k_offset = 0
      IF (PRESENT(k_offset)) my_k_offset = k_offset

      ri = pbc(ctx%particle_set(atom_i)%r(1:3), ctx%cell) + MATMUL(ctx%hmat, REAL(my_cell_i, dp))
      rj = pbc(ctx%particle_set(atom_j)%r(1:3), ctx%cell) + MATMUL(ctx%hmat, REAL(my_cell_j, dp))
      rk = pbc(ctx%particle_set(atom_k)%r(1:3), ctx%cell) + MATMUL(ctx%hmat, REAL(my_cell_k, dp))

      rjk(1:3) = rk(1:3) - rj(1:3)
      rij(1:3) = rj(1:3) - ri(1:3)
      rik(1:3) = rk(1:3) - ri(1:3)

      djk = NORM2(rjk)
      dij = NORM2(rij)
      dik = NORM2(rik)

      ikind = ctx%kind_of(atom_i)
      jkind = ctx%kind_of(atom_j)
      kkind = ctx%kind_of(atom_k)

      CALL get_gto_basis_set(ctx%basis_i(ikind)%gto_basis_set, first_sgf=first_sgf_i, &
                             lmax=lmax_i, lmin=lmin_i, npgf=npgfi, nset=nseti, &
                             nsgf_set=nsgfi, pgf_radius=rpgf_i, set_radius=set_radius_i, &
                             zet=zeti, kind_radius=kind_radius_i)
      CALL get_gto_basis_set(ctx%basis_j(jkind)%gto_basis_set, first_sgf=first_sgf_j, &
                             lmax=lmax_j, lmin=lmin_j, npgf=npgfj, nset=nsetj, &
                             nsgf_set=nsgfj, pgf_radius=rpgf_j, set_radius=set_radius_j, &
                             zet=zetj, kind_radius=kind_radius_j)
      CALL get_gto_basis_set(ctx%basis_k(kkind)%gto_basis_set, first_sgf=first_sgf_k, &
                             lmax=lmax_k, lmin=lmin_k, npgf=npgfk, nset=nsetk, &
                             nsgf_set=nsgfk, pgf_radius=rpgf_k, set_radius=set_radius_k, &
                             zet=zetk, kind_radius=kind_radius_k)

      IF (kind_radius_j + kind_radius_i + ctx%dr_ij < dij .OR. &
          kind_radius_j + kind_radius_k + ctx%dr_jk < djk .OR. &
          kind_radius_k + kind_radius_i + ctx%dr_ik < dik) THEN
         IF (PRESENT(screened)) screened = .TRUE.
         RETURN
      END IF

      DO iset = 1, nseti
         DO jset = 1, nsetj
            IF (set_radius_j(jset) + set_radius_i(iset) + ctx%dr_ij < dij) CYCLE
            DO kset = 1, nsetk
               IF (set_radius_j(jset) + set_radius_k(kset) + ctx%dr_jk < djk) CYCLE
               IF (set_radius_k(kset) + set_radius_i(iset) + ctx%dr_ik < dik) CYCLE

               ncoi = npgfi(iset)*ncoset(lmax_i(iset))
               ncoj = npgfj(jset)*ncoset(lmax_j(jset))
               ncok = npgfk(kset)*ncoset(lmax_k(kset))

               sgfi = first_sgf_i(1, iset)
               sgfj = first_sgf_j(1, jset)
               sgfk = first_sgf_k(1, kset)

               IF (ncoj*ncok*ncoi <= 0) CYCLE
               ALLOCATE (sijk(ncoj, ncok, ncoi))
               sijk(:, :, :) = 0.0_dp

               CALL eri_3center(sijk, &
                                lmin_j(jset), lmax_j(jset), npgfj(jset), zetj(:, jset), &
                                rpgf_j(:, jset), rj, &
                                lmin_k(kset), lmax_k(kset), npgfk(kset), zetk(:, kset), &
                                rpgf_k(:, kset), rk, &
                                lmin_i(iset), lmax_i(iset), npgfi(iset), zeti(:, iset), &
                                rpgf_i(:, iset), ri, &
                                djk, dij, dik, ws%lib, ctx%potential_parameter, &
                                int_abc_ext=sijk_ext)

               ALLOCATE (sijk_contr(nsgfj(jset), nsgfk(kset), nsgfi(iset)))
               CALL abc_contract_xsmm(sijk_contr, sijk, ctx%tspj(jset, jkind)%array, &
                                      ctx%spk(kset, kkind)%array, ctx%spi(iset, ikind)%array, &
                                      ncoj, ncok, ncoi, nsgfj(jset), nsgfk(kset), &
                                      nsgfi(iset), ws%cpp_buffer, ws%ccp_buffer)
               DEALLOCATE (sijk)

               block_start_j = sgfj + my_j_offset
               block_end_j = sgfj + nsgfj(jset) - 1 + my_j_offset
               block_start_k = sgfk + my_k_offset
               block_end_k = sgfk + nsgfk(kset) - 1 + my_k_offset
               block_start_i = sgfi + my_i_offset
               block_end_i = sgfi + nsgfi(iset) - 1 + my_i_offset

               int_3c(block_start_j:block_end_j, &
                      block_start_k:block_end_k, &
                      block_start_i:block_end_i) = &
                  int_3c(block_start_j:block_end_j, &
                         block_start_k:block_end_k, &
                         block_start_i:block_end_i) + &
                  sijk_contr(:, :, :)
               DEALLOCATE (sijk_contr)

            END DO
         END DO
      END DO

   END SUBROUTINE build_3c_integral_block_ctx

! **************************************************************************************************
!> \brief ...
!> \param int_3c ...
!> \param qs_env ...
!> \param potential_parameter ...
!> \param basis_j ...
!> \param basis_k ...
!> \param basis_i ...
!> \param cell_j ...
!> \param cell_k ...
!> \param cell_i ...
!> \param atom_j ...
!> \param atom_k ...
!> \param atom_i ...
!> \param j_bf_start_from_atom ...
!> \param k_bf_start_from_atom ...
!> \param i_bf_start_from_atom ...
! **************************************************************************************************
   SUBROUTINE build_3c_integral_block(int_3c, qs_env, potential_parameter, &
                                      basis_j, basis_k, basis_i, &
                                      cell_j, cell_k, cell_i, atom_j, atom_k, atom_i, &
                                      j_bf_start_from_atom, k_bf_start_from_atom, &
                                      i_bf_start_from_atom)

      REAL(KIND=dp), DIMENSION(:, :, :)                  :: int_3c
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(libint_potential_type), INTENT(IN)            :: potential_parameter
      TYPE(gto_basis_set_p_type), DIMENSION(:)           :: basis_j, basis_k, basis_i
      INTEGER, DIMENSION(3), INTENT(IN), OPTIONAL        :: cell_j, cell_k, cell_i
      INTEGER, INTENT(IN), OPTIONAL                      :: atom_j, atom_k, atom_i
      INTEGER, DIMENSION(:), OPTIONAL                    :: j_bf_start_from_atom, &
                                                            k_bf_start_from_atom, &
                                                            i_bf_start_from_atom

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_3c_integral_block'

      INTEGER                                            :: at_i, at_j, at_k, handle, my_i_offset, &
                                                            my_j_offset, my_k_offset
      TYPE(gw_3c_ctx_type)                               :: ctx
      TYPE(gw_3c_ws_type)                                :: ws

      CALL timeset(routineN, handle)

      CALL gw_3c_ctx_create(ctx, qs_env, potential_parameter, basis_j, basis_k, basis_i)
      CALL gw_3c_ws_create(ws, ctx)

      int_3c(:, :, :) = 0.0_dp

      ! loop over all RI atoms
      DO at_i = 1, ctx%natom
         ! loop over all AO atoms
         DO at_j = 1, ctx%natom
            ! loop over all AO atoms
            DO at_k = 1, ctx%natom

               IF (PRESENT(atom_i)) THEN
                  IF (at_i /= atom_i) CYCLE
               END IF
               IF (PRESENT(atom_j)) THEN
                  IF (at_j /= atom_j) CYCLE
               END IF
               IF (PRESENT(atom_k)) THEN
                  IF (at_k /= atom_k) CYCLE
               END IF

               IF (PRESENT(atom_j)) THEN
                  my_j_offset = 0
               ELSE
                  CPASSERT(PRESENT(j_bf_start_from_atom))
                  my_j_offset = j_bf_start_from_atom(at_j) - 1
               END IF
               IF (PRESENT(atom_k)) THEN
                  my_k_offset = 0
               ELSE
                  CPASSERT(PRESENT(k_bf_start_from_atom))
                  my_k_offset = k_bf_start_from_atom(at_k) - 1
               END IF
               IF (PRESENT(atom_i)) THEN
                  my_i_offset = 0
               ELSE
                  CPASSERT(PRESENT(i_bf_start_from_atom))
                  my_i_offset = i_bf_start_from_atom(at_i) - 1
               END IF

               CALL build_3c_integral_block_ctx(int_3c, ctx, ws, at_j, at_k, at_i, &
                                                cell_j=cell_j, cell_k=cell_k, cell_i=cell_i, &
                                                j_offset=my_j_offset, k_offset=my_k_offset, &
                                                i_offset=my_i_offset)

            END DO ! atom_k (AO)
         END DO ! atom_j (AO)
      END DO ! atom_i (RI)

      CALL gw_3c_ws_release(ws)
      CALL gw_3c_ctx_release(ctx)

      CALL timestop(handle)

   END SUBROUTINE build_3c_integral_block

END MODULE gw_integrals

