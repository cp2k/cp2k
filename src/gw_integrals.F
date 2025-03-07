!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2025 CP2K developers group <https://cp2k.org>                                   !
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
   USE cp_control_types,                ONLY: dft_control_type
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
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE t_c_g0,                          ONLY: get_lmax_init,&
                                              init

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gw_integrals'

   PUBLIC :: build_3c_integral_block

CONTAINS

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

      INTEGER :: at_i, at_j, at_k, block_end_i, block_end_j, block_end_k, block_start_i, &
         block_start_j, block_start_k, egfi, handle, i, i_offset, ibasis, ikind, ilist, imax, is, &
         iset, j_offset, jkind, js, jset, k_offset, kkind, ks, kset, m_max, max_ncoi, max_ncoj, &
         max_ncok, max_nset, max_nsgfi, max_nsgfj, max_nsgfk, maxli, maxlj, maxlk, natom, nbasis, &
         ncoi, ncoj, ncok, nseti, nsetj, nsetk, op_ij, op_jk, sgfi, sgfj, sgfk, unit_id
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: kind_of
      INTEGER, DIMENSION(3)                              :: my_cell_i, my_cell_j, my_cell_k
      INTEGER, DIMENSION(:), POINTER                     :: lmax_i, lmax_j, lmax_k, lmin_i, lmin_j, &
                                                            lmin_k, npgfi, npgfj, npgfk, nsgfi, &
                                                            nsgfj, nsgfk
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgf_i, first_sgf_j, first_sgf_k
      REAL(KIND=dp)                                      :: dij, dik, djk, dr_ij, dr_ik, dr_jk, &
                                                            kind_radius_i, kind_radius_j, &
                                                            kind_radius_k, sijk_ext
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: ccp_buffer, cpp_buffer, &
                                                            max_contraction_i, max_contraction_j, &
                                                            max_contraction_k
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: sijk, sijk_contr
      REAL(KIND=dp), DIMENSION(3)                        :: ri, rij, rik, rj, rjk, rk
      REAL(KIND=dp), DIMENSION(3, 3)                     :: hmat
      REAL(KIND=dp), DIMENSION(:), POINTER               :: set_radius_i, set_radius_j, set_radius_k
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: rpgf_i, rpgf_j, rpgf_k, sphi_i, sphi_j, &
                                                            sphi_k, zeti, zetj, zetk
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_2d_r_p_type), DIMENSION(:, :), POINTER     :: spi, spk, tspj
      TYPE(cp_libint_t)                                  :: lib
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_type), POINTER                  :: basis_set
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      op_ij = potential_parameter%potential_type
      op_jk = do_potential_id

      dr_ij = 0.0_dp; dr_jk = 0.0_dp; dr_ik = 0.0_dp

      IF (op_ij == do_potential_truncated .OR. op_ij == do_potential_short) THEN
         dr_ij = potential_parameter%cutoff_radius*cutoff_screen_factor
         dr_ik = potential_parameter%cutoff_radius*cutoff_screen_factor
      ELSEIF (op_ij == do_potential_coulomb) THEN
         dr_ij = 1000000.0_dp
         dr_ik = 1000000.0_dp
      END IF

      NULLIFY (qs_kind_set, atomic_kind_set)

      ! get stuff
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set, cell=cell, &
                      natom=natom, dft_control=dft_control, para_env=para_env, &
                      particle_set=particle_set)
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, kind_of=kind_of)
      CALL get_cell(cell=cell, h=hmat)

      !Need the max l for each basis for libint and max nset, nco and nsgf for LIBXSMM contraction
      nbasis = SIZE(basis_i)
      max_nsgfi = 0
      max_ncoi = 0
      max_nset = 0
      maxli = 0
      DO ibasis = 1, nbasis
         CALL get_gto_basis_set(gto_basis_set=basis_i(ibasis)%gto_basis_set, maxl=imax, &
                                nset=iset, nsgf_set=nsgfi, npgf=npgfi)
         maxli = MAX(maxli, imax)
         max_nset = MAX(max_nset, iset)
         max_nsgfi = MAX(max_nsgfi, MAXVAL(nsgfi))
         max_ncoi = MAX(max_ncoi, MAXVAL(npgfi)*ncoset(maxli))
      END DO
      max_nsgfj = 0
      max_ncoj = 0
      maxlj = 0
      DO ibasis = 1, nbasis
         CALL get_gto_basis_set(gto_basis_set=basis_j(ibasis)%gto_basis_set, maxl=imax, &
                                nset=jset, nsgf_set=nsgfj, npgf=npgfj)
         maxlj = MAX(maxlj, imax)
         max_nset = MAX(max_nset, jset)
         max_nsgfj = MAX(max_nsgfj, MAXVAL(nsgfj))
         max_ncoj = MAX(max_ncoj, MAXVAL(npgfj)*ncoset(maxlj))
      END DO
      max_nsgfk = 0
      max_ncok = 0
      maxlk = 0
      DO ibasis = 1, nbasis
         CALL get_gto_basis_set(gto_basis_set=basis_k(ibasis)%gto_basis_set, maxl=imax, &
                                nset=kset, nsgf_set=nsgfk, npgf=npgfk)
         maxlk = MAX(maxlk, imax)
         max_nset = MAX(max_nset, kset)
         max_nsgfk = MAX(max_nsgfk, MAXVAL(nsgfk))
         max_ncok = MAX(max_ncok, MAXVAL(npgfk)*ncoset(maxlk))
      END DO
      m_max = maxli + maxlj + maxlk

      !To minimize expensive memory opsand generally optimize contraction, pre-allocate
      !contiguous sphi arrays (and transposed in the cas of sphi_i)

      NULLIFY (tspj, spi, spk)
      ALLOCATE (spi(max_nset, nbasis), tspj(max_nset, nbasis), spk(max_nset, nbasis))

      DO ibasis = 1, nbasis
         DO iset = 1, max_nset
            NULLIFY (spi(iset, ibasis)%array)
            NULLIFY (tspj(iset, ibasis)%array)

            NULLIFY (spk(iset, ibasis)%array)
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
                  ALLOCATE (spi(iset, ibasis)%array(ncoi, basis_set%nsgf_set(iset)))
                  spi(iset, ibasis)%array(:, :) = basis_set%sphi(1:ncoi, sgfi:egfi)

               ELSE IF (ilist == 2) THEN
                  ALLOCATE (tspj(iset, ibasis)%array(basis_set%nsgf_set(iset), ncoi))
                  tspj(iset, ibasis)%array(:, :) = TRANSPOSE(basis_set%sphi(1:ncoi, sgfi:egfi))

               ELSE
                  ALLOCATE (spk(iset, ibasis)%array(ncoi, basis_set%nsgf_set(iset)))
                  spk(iset, ibasis)%array(:, :) = basis_set%sphi(1:ncoi, sgfi:egfi)
               END IF

            END DO !iset
         END DO !ibasis
      END DO !ilist

      !Init the truncated Coulomb operator
      IF (op_ij == do_potential_truncated .OR. op_jk == do_potential_truncated) THEN

         IF (m_max > get_lmax_init()) THEN
            IF (para_env%mepos == 0) THEN
               CALL open_file(unit_number=unit_id, file_name=potential_parameter%filename)
            END IF
            CALL init(m_max, unit_id, para_env%mepos, para_env)
            IF (para_env%mepos == 0) THEN
               CALL close_file(unit_id)
            END IF
         END IF
      END IF

      CALL init_md_ftable(nmax=m_max)

      CALL cp_libint_init_3eri(lib, MAX(maxli, maxlj, maxlk))
      CALL cp_libint_set_contrdepth(lib, 1)

      !pre-allocate contraction buffers
      ALLOCATE (cpp_buffer(max_nsgfj*max_ncok), ccp_buffer(max_nsgfj*max_nsgfk*max_ncoi))
      int_3c(:, :, :) = 0.0_dp

      ! loop over all RI atoms
      DO at_i = 1, natom

         ! loop over all AO atoms
         DO at_j = 1, natom

            ! loop over all AO atoms
            DO at_k = 1, natom

               IF (PRESENT(atom_i)) THEN
                  IF (at_i .NE. atom_i) CYCLE
               END IF
               IF (PRESENT(atom_j)) THEN
                  IF (at_j .NE. atom_j) CYCLE
               END IF
               IF (PRESENT(atom_k)) THEN
                  IF (at_k .NE. atom_k) CYCLE
               END IF

               my_cell_i(1:3) = 0
               IF (PRESENT(cell_i)) my_cell_i(1:3) = cell_i(1:3)
               my_cell_j(1:3) = 0
               IF (PRESENT(cell_j)) my_cell_j(1:3) = cell_j(1:3)
               my_cell_k(1:3) = 0
               IF (PRESENT(cell_k)) my_cell_k(1:3) = cell_k(1:3)

               ri = pbc(particle_set(at_i)%r(1:3), cell) + MATMUL(hmat, REAL(my_cell_i, dp))
               rj = pbc(particle_set(at_j)%r(1:3), cell) + MATMUL(hmat, REAL(my_cell_j, dp))
               rk = pbc(particle_set(at_k)%r(1:3), cell) + MATMUL(hmat, REAL(my_cell_k, dp))

               rjk(1:3) = rk(1:3) - rj(1:3)
               rij(1:3) = rj(1:3) - ri(1:3)
               rik(1:3) = rk(1:3) - ri(1:3)

               djk = NORM2(rjk)
               dij = NORM2(rij)
               dik = NORM2(rik)

               ikind = kind_of(at_i)
               jkind = kind_of(at_j)
               kkind = kind_of(at_k)

               CALL get_gto_basis_set(basis_i(ikind)%gto_basis_set, first_sgf=first_sgf_i, &
                                      lmax=lmax_i, lmin=lmin_i, npgf=npgfi, nset=nseti, &
                                      nsgf_set=nsgfi, pgf_radius=rpgf_i, set_radius=set_radius_i, &
                                      sphi=sphi_i, zet=zeti, kind_radius=kind_radius_i)

               CALL get_gto_basis_set(basis_j(jkind)%gto_basis_set, first_sgf=first_sgf_j, &
                                      lmax=lmax_j, lmin=lmin_j, npgf=npgfj, nset=nsetj, &
                                      nsgf_set=nsgfj, pgf_radius=rpgf_j, set_radius=set_radius_j, &
                                      sphi=sphi_j, zet=zetj, kind_radius=kind_radius_j)

               CALL get_gto_basis_set(basis_k(kkind)%gto_basis_set, first_sgf=first_sgf_k, &
                                      lmax=lmax_k, lmin=lmin_k, npgf=npgfk, nset=nsetk, &
                                      nsgf_set=nsgfk, pgf_radius=rpgf_k, set_radius=set_radius_k, &
                                      sphi=sphi_k, zet=zetk, kind_radius=kind_radius_k)

               IF (kind_radius_j + kind_radius_i + dr_ij < dij) CYCLE
               IF (kind_radius_j + kind_radius_k + dr_jk < djk) CYCLE
               IF (kind_radius_k + kind_radius_i + dr_ik < dik) CYCLE

               ALLOCATE (max_contraction_i(nseti))
               max_contraction_i = 0.0_dp
               DO iset = 1, nseti
                  sgfi = first_sgf_i(1, iset)
                  max_contraction_i(iset) = MAXVAL((/(SUM(ABS(sphi_i(:, i))), i=sgfi, &
                                                      sgfi + nsgfi(iset) - 1)/))
               END DO

               ALLOCATE (max_contraction_j(nsetj))
               max_contraction_j = 0.0_dp
               DO jset = 1, nsetj
                  sgfj = first_sgf_j(1, jset)
                  max_contraction_j(jset) = MAXVAL((/(SUM(ABS(sphi_j(:, i))), i=sgfj, &
                                                      sgfj + nsgfj(jset) - 1)/))
               END DO

               ALLOCATE (max_contraction_k(nsetk))
               max_contraction_k = 0.0_dp
               DO kset = 1, nsetk
                  sgfk = first_sgf_k(1, kset)
                  max_contraction_k(kset) = MAXVAL((/(SUM(ABS(sphi_k(:, i))), i=sgfk, &
                                                      sgfk + nsgfk(kset) - 1)/))
               END DO

               DO iset = 1, nseti

                  DO jset = 1, nsetj

                     IF (set_radius_j(jset) + set_radius_i(iset) + dr_ij < dij) CYCLE

                     DO kset = 1, nsetk

                        IF (set_radius_j(jset) + set_radius_k(kset) + dr_jk < djk) CYCLE
                        IF (set_radius_k(kset) + set_radius_i(iset) + dr_ik < dik) CYCLE

                        ncoi = npgfi(iset)*ncoset(lmax_i(iset))
                        ncoj = npgfj(jset)*ncoset(lmax_j(jset))
                        ncok = npgfk(kset)*ncoset(lmax_k(kset))

                        sgfi = first_sgf_i(1, iset)
                        sgfj = first_sgf_j(1, jset)
                        sgfk = first_sgf_k(1, kset)

                        IF (ncoj*ncok*ncoi .LE. 0) CYCLE
                        ALLOCATE (sijk(ncoj, ncok, ncoi))
                        sijk(:, :, :) = 0.0_dp

                        is = iset
                        js = jset
                        ks = kset

                        CALL eri_3center(sijk, &
                                         lmin_j(js), lmax_j(js), npgfj(js), zetj(:, js), &
                                         rpgf_j(:, js), rj, &
                                         lmin_k(ks), lmax_k(ks), npgfk(ks), zetk(:, ks), &
                                         rpgf_k(:, ks), rk, &
                                         lmin_i(is), lmax_i(is), npgfi(is), zeti(:, is), &
                                         rpgf_i(:, is), ri, &
                                         djk, dij, dik, lib, potential_parameter, &
                                         int_abc_ext=sijk_ext)

                        ALLOCATE (sijk_contr(nsgfj(jset), nsgfk(kset), nsgfi(iset)))
                        CALL abc_contract_xsmm(sijk_contr, sijk, tspj(jset, jkind)%array, &
                                               spk(kset, kkind)%array, spi(iset, ikind)%array, &
                                               ncoj, ncok, ncoi, nsgfj(jset), nsgfk(kset), &
                                               nsgfi(iset), cpp_buffer, ccp_buffer)
                        DEALLOCATE (sijk)

                        IF (PRESENT(atom_j)) THEN
                           j_offset = 0
                        ELSE
                           CPASSERT(PRESENT(j_bf_start_from_atom))
                           j_offset = j_bf_start_from_atom(at_j) - 1
                        END IF
                        IF (PRESENT(atom_k)) THEN
                           k_offset = 0
                        ELSE
                           CPASSERT(PRESENT(k_bf_start_from_atom))
                           k_offset = k_bf_start_from_atom(at_k) - 1
                        END IF
                        IF (PRESENT(atom_i)) THEN
                           i_offset = 0
                        ELSE
                           CPASSERT(PRESENT(i_bf_start_from_atom))
                           i_offset = i_bf_start_from_atom(at_i) - 1
                        END IF

                        block_start_j = sgfj + j_offset
                        block_end_j = sgfj + nsgfj(jset) - 1 + j_offset
                        block_start_k = sgfk + k_offset
                        block_end_k = sgfk + nsgfk(kset) - 1 + k_offset
                        block_start_i = sgfi + i_offset
                        block_end_i = sgfi + nsgfi(iset) - 1 + i_offset

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

               DEALLOCATE (max_contraction_i, max_contraction_j, max_contraction_k)

            END DO ! atom_k (AO)
         END DO ! atom_j (AO)
      END DO ! atom_i (RI)

      CALL cp_libint_cleanup_3eri(lib)

      DO iset = 1, max_nset
         DO ibasis = 1, nbasis
            IF (ASSOCIATED(spi(iset, ibasis)%array)) DEALLOCATE (spi(iset, ibasis)%array)
            IF (ASSOCIATED(tspj(iset, ibasis)%array)) DEALLOCATE (tspj(iset, ibasis)%array)

            IF (ASSOCIATED(spk(iset, ibasis)%array)) DEALLOCATE (spk(iset, ibasis)%array)
         END DO
      END DO
      DEALLOCATE (spi, tspj, spk)

      CALL timestop(handle)

   END SUBROUTINE build_3c_integral_block

END MODULE

