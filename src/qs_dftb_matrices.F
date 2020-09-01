!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2020 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of Overlap and Hamiltonian matrices in DFTB
!> \author JGH
! **************************************************************************************************
MODULE qs_dftb_matrices
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              get_atomic_kind_set
   USE atprop_types,                    ONLY: atprop_array_init,&
                                              atprop_type
   USE block_p_types,                   ONLY: block_p_type
   USE cp_control_types,                ONLY: dft_control_type,&
                                              dftb_control_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set
   USE cp_dbcsr_output,                 ONLY: cp_dbcsr_write_sparse_matrix
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_convert_offsets_to_sizes, dbcsr_copy, dbcsr_create, &
        dbcsr_distribution_type, dbcsr_dot, dbcsr_finalize, dbcsr_get_block_p, dbcsr_multiply, &
        dbcsr_p_type, dbcsr_type, dbcsr_type_antisymmetric, dbcsr_type_symmetric
   USE efield_tb_methods,               ONLY: efield_tb_matrix
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_type
   USE message_passing,                 ONLY: mp_sum
   USE mulliken,                        ONLY: mulliken_charges
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE qs_dftb_coulomb,                 ONLY: build_dftb_coulomb
   USE qs_dftb_types,                   ONLY: qs_dftb_atom_type,&
                                              qs_dftb_pairpot_type
   USE qs_dftb_utils,                   ONLY: compute_block_sk,&
                                              get_dftb_atom_param,&
                                              iptr,&
                                              urep_egr
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type,&
                                              set_ks_env
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE virial_methods,                  ONLY: virial_pair_force
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   INTEGER, DIMENSION(16), PARAMETER        :: orbptr = (/0, 1, 1, 1, &
                                                          2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3/)

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_dftb_matrices'

   PUBLIC :: build_dftb_matrices, build_dftb_ks_matrix, build_dftb_overlap

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param para_env ...
!> \param calculate_forces ...
! **************************************************************************************************
   SUBROUTINE build_dftb_matrices(qs_env, para_env, calculate_forces)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_para_env_type), POINTER                    :: para_env
      LOGICAL, INTENT(IN)                                :: calculate_forces

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_dftb_matrices'

      INTEGER :: after, atom_a, atom_b, handle, i, iatom, ic, icol, ikind, img, irow, iw, jatom, &
         jkind, l1, l2, la, lb, llm, lmaxi, lmaxj, m, n1, n2, n_urpoly, natom, natorb_a, natorb_b, &
         nderivatives, ngrd, ngrdcut, nimg, nkind, spdim
      INTEGER, DIMENSION(3)                              :: cell
      INTEGER, DIMENSION(:), POINTER                     :: atom_of_kind
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: defined, found, omit_headers, use_virial
      REAL(KIND=dp)                                      :: ddr, dgrd, dr, erep, erepij, f0, f1, &
                                                            foab, fow, s_cut, urep_cut
      REAL(KIND=dp), DIMENSION(0:3)                      :: eta_a, eta_b, skself
      REAL(KIND=dp), DIMENSION(10)                       :: urep
      REAL(KIND=dp), DIMENSION(2)                        :: surr
      REAL(KIND=dp), DIMENSION(3)                        :: drij, force_ab, force_rr, force_w, rij, &
                                                            srep
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: dfblock, dsblock, fblock, fmatij, &
                                                            fmatji, pblock, sblock, scoeff, &
                                                            smatij, smatji, spxr, wblock
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atprop_type), POINTER                         :: atprop
      TYPE(block_p_type), DIMENSION(2:4)                 :: dsblocks
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_h, matrix_p, matrix_s, matrix_w
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(dftb_control_type), POINTER                   :: dftb_control
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_dftb_atom_type), POINTER                   :: dftb_kind_a, dftb_kind_b
      TYPE(qs_dftb_pairpot_type), DIMENSION(:, :), &
         POINTER                                         :: dftb_potential
      TYPE(qs_dftb_pairpot_type), POINTER                :: dftb_param_ij, dftb_param_ji
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      ! set pointers
      iptr = 0
      DO la = 0, 3
         DO lb = 0, 3
            llm = 0
            DO l1 = 0, MAX(la, lb)
               DO l2 = 0, MIN(l1, la, lb)
                  DO m = 0, l2
                     llm = llm + 1
                     iptr(l1, l2, m, la, lb) = llm
                  END DO
               END DO
            END DO
         END DO
      END DO

      NULLIFY (logger, virial, atprop)
      logger => cp_get_default_logger()

      NULLIFY (matrix_h, matrix_s, matrix_p, matrix_w, atomic_kind_set, &
               qs_kind_set, sab_orb, ks_env)

      CALL get_qs_env(qs_env=qs_env, &
                      energy=energy, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      matrix_h_kp=matrix_h, &
                      matrix_s_kp=matrix_s, &
                      atprop=atprop, &
                      dft_control=dft_control, &
                      ks_env=ks_env)

      dftb_control => dft_control%qs_control%dftb_control
      nimg = dft_control%nimages
      ! Allocate the overlap and Hamiltonian matrix
      CALL get_qs_env(qs_env=qs_env, sab_orb=sab_orb)
      nderivatives = 0
      IF (dftb_control%self_consistent .AND. calculate_forces) nderivatives = 1
      CALL setup_matrices2(qs_env, nderivatives, nimg, matrix_s, "OVERLAP", sab_orb)
      CALL setup_matrices2(qs_env, 0, nimg, matrix_h, "CORE HAMILTONIAN", sab_orb)
      CALL set_ks_env(ks_env, matrix_s_kp=matrix_s)
      CALL set_ks_env(ks_env, matrix_h_kp=matrix_h)

      NULLIFY (dftb_potential)
      CALL get_qs_env(qs_env=qs_env, dftb_potential=dftb_potential)
      NULLIFY (particle_set)
      CALL get_qs_env(qs_env=qs_env, particle_set=particle_set)

      IF (calculate_forces) THEN
         NULLIFY (rho, force, matrix_w)
         CALL get_qs_env(qs_env=qs_env, &
                         rho=rho, &
                         matrix_w_kp=matrix_w, &
                         virial=virial, &
                         force=force)
         CALL qs_rho_get(rho, rho_ao_kp=matrix_p)

         IF (SIZE(matrix_p, 1) == 2) THEN
            DO img = 1, nimg
               CALL dbcsr_add(matrix_p(1, img)%matrix, matrix_p(2, img)%matrix, &
                              alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
               CALL dbcsr_add(matrix_w(1, img)%matrix, matrix_w(2, img)%matrix, &
                              alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
            END DO
         END IF
         natom = SIZE(particle_set)
         ALLOCATE (atom_of_kind(natom))
         CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                                  atom_of_kind=atom_of_kind)
         use_virial = virial%pv_availability .AND. (.NOT. virial%pv_numer)
      END IF
      ! atomic energy decomposition
      IF (atprop%energy) THEN
         natom = SIZE(particle_set)
         CALL atprop_array_init(atprop%atecc, natom)
      END IF

      NULLIFY (cell_to_index)
      IF (nimg > 1) THEN
         CALL get_ks_env(ks_env=ks_env, kpoints=kpoints)
         CALL get_kpoint_info(kpoint=kpoints, cell_to_index=cell_to_index)
      END IF

      erep = 0._dp

      nkind = SIZE(atomic_kind_set)

      CALL neighbor_list_iterator_create(nl_iterator, sab_orb)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, &
                                iatom=iatom, jatom=jatom, r=rij, cell=cell)
         CALL get_atomic_kind(atomic_kind_set(ikind), natom=natom)
         CALL get_qs_kind(qs_kind_set(ikind), dftb_parameter=dftb_kind_a)
         CALL get_dftb_atom_param(dftb_kind_a, &
                                  defined=defined, lmax=lmaxi, skself=skself, &
                                  eta=eta_a, natorb=natorb_a)
         IF (.NOT. defined .OR. natorb_a < 1) CYCLE
         CALL get_qs_kind(qs_kind_set(jkind), dftb_parameter=dftb_kind_b)
         CALL get_dftb_atom_param(dftb_kind_b, &
                                  defined=defined, lmax=lmaxj, eta=eta_b, natorb=natorb_b)

         IF (.NOT. defined .OR. natorb_b < 1) CYCLE

         ! retrieve information on F and S matrix
         dftb_param_ij => dftb_potential(ikind, jkind)
         dftb_param_ji => dftb_potential(jkind, ikind)
         ! assume table size and type is symmetric
         ngrd = dftb_param_ij%ngrd
         ngrdcut = dftb_param_ij%ngrdcut
         dgrd = dftb_param_ij%dgrd
         ddr = dgrd*0.1_dp
         CPASSERT(dftb_param_ij%llm == dftb_param_ji%llm)
         llm = dftb_param_ij%llm
         fmatij => dftb_param_ij%fmat
         smatij => dftb_param_ij%smat
         fmatji => dftb_param_ji%fmat
         smatji => dftb_param_ji%smat
         ! repulsive pair potential
         n_urpoly = dftb_param_ij%n_urpoly
         urep_cut = dftb_param_ij%urep_cut
         urep = dftb_param_ij%urep
         spxr => dftb_param_ij%spxr
         scoeff => dftb_param_ij%scoeff
         spdim = dftb_param_ij%spdim
         s_cut = dftb_param_ij%s_cut
         srep = dftb_param_ij%srep
         surr = dftb_param_ij%surr

         dr = SQRT(SUM(rij(:)**2))
         IF (NINT(dr/dgrd) <= ngrdcut) THEN

            IF (nimg == 1) THEN
               ic = 1
            ELSE
               ic = cell_to_index(cell(1), cell(2), cell(3))
               CPASSERT(ic > 0)
            END IF

            icol = MAX(iatom, jatom)
            irow = MIN(iatom, jatom)
            NULLIFY (sblock, fblock)
            CALL dbcsr_get_block_p(matrix=matrix_s(1, ic)%matrix, &
                                   row=irow, col=icol, BLOCK=sblock, found=found)
            CPASSERT(found)
            CALL dbcsr_get_block_p(matrix=matrix_h(1, ic)%matrix, &
                                   row=irow, col=icol, BLOCK=fblock, found=found)
            CPASSERT(found)

            IF (calculate_forces) THEN
               NULLIFY (pblock)
               CALL dbcsr_get_block_p(matrix=matrix_p(1, ic)%matrix, &
                                      row=irow, col=icol, block=pblock, found=found)
               CPASSERT(ASSOCIATED(pblock))
               NULLIFY (wblock)
               CALL dbcsr_get_block_p(matrix=matrix_w(1, ic)%matrix, &
                                      row=irow, col=icol, block=wblock, found=found)
               CPASSERT(ASSOCIATED(wblock))
               IF (dftb_control%self_consistent) THEN
                  DO i = 2, 4
                     NULLIFY (dsblocks(i)%block)
                     CALL dbcsr_get_block_p(matrix=matrix_s(i, ic)%matrix, &
                                            row=irow, col=icol, BLOCK=dsblocks(i)%block, found=found)
                     CPASSERT(found)
                  END DO
               END IF
            END IF

            IF (iatom == jatom .AND. dr < 0.001_dp) THEN
               ! diagonal block
               DO i = 1, natorb_a
                  sblock(i, i) = sblock(i, i) + 1._dp
                  fblock(i, i) = fblock(i, i) + skself(orbptr(i))
               END DO
            ELSE
               ! off-diagonal block
               CALL compute_block_sk(sblock, smatij, smatji, rij, ngrd, ngrdcut, dgrd, &
                                     llm, lmaxi, lmaxj, irow, iatom)
               CALL compute_block_sk(fblock, fmatij, fmatji, rij, ngrd, ngrdcut, dgrd, &
                                     llm, lmaxi, lmaxj, irow, iatom)
               IF (calculate_forces) THEN
                  force_ab = 0._dp
                  force_w = 0._dp
                  n1 = SIZE(fblock, 1)
                  n2 = SIZE(fblock, 2)
                  ! make sure that displacement is in the correct direction depending on the position
                  ! of the block (upper or lower triangle)
                  f0 = 1.0_dp
                  IF (irow == iatom) f0 = -1.0_dp

                  ALLOCATE (dfblock(n1, n2), dsblock(n1, n2))

                  DO i = 1, 3
                     drij = rij
                     dfblock = 0._dp; dsblock = 0._dp

                     drij(i) = rij(i) - ddr*f0
                     CALL compute_block_sk(dsblock, smatij, smatji, drij, ngrd, ngrdcut, dgrd, &
                                           llm, lmaxi, lmaxj, irow, iatom)
                     CALL compute_block_sk(dfblock, fmatij, fmatji, drij, ngrd, ngrdcut, dgrd, &
                                           llm, lmaxi, lmaxj, irow, iatom)

                     dsblock = -dsblock
                     dfblock = -dfblock

                     drij(i) = rij(i) + ddr*f0
                     CALL compute_block_sk(dsblock, smatij, smatji, drij, ngrd, ngrdcut, dgrd, &
                                           llm, lmaxi, lmaxj, irow, iatom)
                     CALL compute_block_sk(dfblock, fmatij, fmatji, drij, ngrd, ngrdcut, dgrd, &
                                           llm, lmaxi, lmaxj, irow, iatom)

                     dfblock = dfblock/(2.0_dp*ddr)
                     dsblock = dsblock/(2.0_dp*ddr)

                     foab = 2.0_dp*SUM(dfblock*pblock)
                     fow = -2.0_dp*SUM(dsblock*wblock)

                     force_ab(i) = force_ab(i) + foab
                     force_w(i) = force_w(i) + fow
                     IF (dftb_control%self_consistent) THEN
                        CPASSERT(ASSOCIATED(dsblocks(i + 1)%block))
                        dsblocks(i + 1)%block = dsblocks(i + 1)%block + dsblock
                     END IF
                  ENDDO
                  IF (use_virial) THEN
!deb iatom=jatom f0=0.5*f0
                     IF (iatom == jatom) f0 = 0.5_dp*f0
!deb
                     CALL virial_pair_force(virial%pv_virial, -f0, force_ab, rij)
                     CALL virial_pair_force(virial%pv_virial, -f0, force_w, rij)
                     IF (atprop%stress) THEN
                        f1 = 0.5_dp*f0
                        CALL virial_pair_force(atprop%atstress(:, :, iatom), -f1, force_ab, rij)
                        CALL virial_pair_force(atprop%atstress(:, :, iatom), -f1, force_w, rij)
                        CALL virial_pair_force(atprop%atstress(:, :, jatom), -f1, force_ab, rij)
                        CALL virial_pair_force(atprop%atstress(:, :, jatom), -f1, force_w, rij)
                     END IF
                  END IF
                  DEALLOCATE (dfblock, dsblock)
               END IF
            END IF

            IF (calculate_forces .AND. (iatom /= jatom .OR. dr > 0.001_dp)) THEN
               atom_a = atom_of_kind(iatom)
               atom_b = atom_of_kind(jatom)
               IF (irow == iatom) force_ab = -force_ab
               IF (irow == iatom) force_w = -force_w
               force(ikind)%all_potential(:, atom_a) = force(ikind)%all_potential(:, atom_a) - force_ab(:)
               force(jkind)%all_potential(:, atom_b) = force(jkind)%all_potential(:, atom_b) + force_ab(:)
               force(ikind)%overlap(:, atom_a) = force(ikind)%overlap(:, atom_a) - force_w(:)
               force(jkind)%overlap(:, atom_b) = force(jkind)%overlap(:, atom_b) + force_w(:)
            END IF

         END IF

         ! repulsive potential
         IF ((dr <= urep_cut .OR. spdim > 0) .AND. dr > 0.001_dp) THEN
            erepij = 0._dp
            CALL urep_egr(rij, dr, erepij, force_rr, &
                          n_urpoly, urep, spdim, s_cut, srep, spxr, scoeff, surr, calculate_forces)
            erep = erep + erepij
            IF (atprop%energy) THEN
               atprop%atecc(iatom) = atprop%atecc(iatom) + 0.5_dp*erepij
               atprop%atecc(jatom) = atprop%atecc(jatom) + 0.5_dp*erepij
            END IF
            IF (calculate_forces .AND. (iatom /= jatom .OR. dr > 0.001_dp)) THEN
               atom_a = atom_of_kind(iatom)
               atom_b = atom_of_kind(jatom)
               force(ikind)%repulsive(:, atom_a) = &
                  force(ikind)%repulsive(:, atom_a) - force_rr(:)
               force(jkind)%repulsive(:, atom_b) = &
                  force(jkind)%repulsive(:, atom_b) + force_rr(:)
               IF (use_virial) THEN
                  f0 = -1.0_dp
                  IF (iatom == jatom) f0 = -0.5_dp
                  CALL virial_pair_force(virial%pv_virial, f0, force_rr, rij)
                  IF (atprop%stress) THEN
                     CALL virial_pair_force(atprop%atstress(:, :, iatom), f0*0.5_dp, force_rr, rij)
                     CALL virial_pair_force(atprop%atstress(:, :, jatom), f0*0.5_dp, force_rr, rij)
                  END IF
               END IF
            END IF
         END IF

      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      DO i = 1, SIZE(matrix_s, 1)
         DO img = 1, nimg
            CALL dbcsr_finalize(matrix_s(i, img)%matrix)
         END DO
      ENDDO
      DO i = 1, SIZE(matrix_h, 1)
         DO img = 1, nimg
            CALL dbcsr_finalize(matrix_h(i, img)%matrix)
         END DO
      ENDDO

      ! set repulsive energy
      CALL mp_sum(erep, para_env%group)
      energy%repulsive = erep

      CALL section_vals_val_get(qs_env%input, "DFT%PRINT%AO_MATRICES%OMIT_HEADERS", l_val=omit_headers)
      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                           qs_env%input, "DFT%PRINT%AO_MATRICES/CORE_HAMILTONIAN"), cp_p_file)) THEN
         iw = cp_print_key_unit_nr(logger, qs_env%input, "DFT%PRINT%AO_MATRICES/CORE_HAMILTONIAN", &
                                   extension=".Log")
         CALL section_vals_val_get(qs_env%input, "DFT%PRINT%AO_MATRICES%NDIGITS", i_val=after)
         after = MIN(MAX(after, 1), 16)
         DO img = 1, nimg
            CALL cp_dbcsr_write_sparse_matrix(matrix_h(1, img)%matrix, 4, after, qs_env, para_env, &
                                              output_unit=iw, omit_headers=omit_headers)
         END DO

         CALL cp_print_key_finished_output(iw, logger, qs_env%input, &
                                           "DFT%PRINT%AO_MATRICES/CORE_HAMILTONIAN")
      END IF

      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                           qs_env%input, "DFT%PRINT%AO_MATRICES/OVERLAP"), cp_p_file)) THEN
         iw = cp_print_key_unit_nr(logger, qs_env%input, "DFT%PRINT%AO_MATRICES/OVERLAP", &
                                   extension=".Log")
         CALL section_vals_val_get(qs_env%input, "DFT%PRINT%AO_MATRICES%NDIGITS", i_val=after)
         after = MIN(MAX(after, 1), 16)
         DO img = 1, nimg
            CALL cp_dbcsr_write_sparse_matrix(matrix_s(1, img)%matrix, 4, after, qs_env, para_env, &
                                              output_unit=iw, omit_headers=omit_headers)

            IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                                 qs_env%input, "DFT%PRINT%AO_MATRICES/DERIVATIVES"), cp_p_file)) THEN
               DO i = 2, SIZE(matrix_s, 1)
                  CALL cp_dbcsr_write_sparse_matrix(matrix_s(i, img)%matrix, 4, after, qs_env, para_env, &
                                                    output_unit=iw, omit_headers=omit_headers)
               END DO
            END IF
         END DO

         CALL cp_print_key_finished_output(iw, logger, qs_env%input, &
                                           "DFT%PRINT%AO_MATRICES/OVERLAP")
      END IF

      IF (calculate_forces) THEN
         IF (SIZE(matrix_p, 1) == 2) THEN
            DO img = 1, nimg
               CALL dbcsr_add(matrix_p(1, img)%matrix, matrix_p(2, img)%matrix, alpha_scalar=1.0_dp, &
                              beta_scalar=-1.0_dp)
            END DO
         END IF
         DEALLOCATE (atom_of_kind)
      END IF

      CALL timestop(handle)

   END SUBROUTINE build_dftb_matrices

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param calculate_forces ...
!> \param just_energy ...
! **************************************************************************************************
   SUBROUTINE build_dftb_ks_matrix(qs_env, calculate_forces, just_energy)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(in)                                :: calculate_forces, just_energy

      CHARACTER(len=*), PARAMETER :: routineN = 'build_dftb_ks_matrix'

      INTEGER                                            :: atom_a, handle, iatom, ikind, img, &
                                                            ispin, natom, nkind, nspins, &
                                                            output_unit
      REAL(KIND=dp)                                      :: pc_ener, qmmm_el, zeff
      REAL(KIND=dp), DIMENSION(:), POINTER               :: mcharge, occupation_numbers
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: charges
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_p1, mo_derivs
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ks_matrix, matrix_h, matrix_p, matrix_s
      TYPE(dbcsr_type), POINTER                          :: mo_coeff
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mo_array
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_dftb_atom_type), POINTER                   :: dftb_kind
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: scf_section

      CALL timeset(routineN, handle)
      NULLIFY (dft_control, logger, scf_section, matrix_p, particle_set, ks_env, &
               ks_matrix, rho, energy)
      logger => cp_get_default_logger()
      CPASSERT(ASSOCIATED(qs_env))

      CALL get_qs_env(qs_env, &
                      dft_control=dft_control, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      matrix_h_kp=matrix_h, &
                      para_env=para_env, &
                      ks_env=ks_env, &
                      matrix_ks_kp=ks_matrix, &
                      rho=rho, &
                      energy=energy)

      energy%hartree = 0.0_dp
      energy%qmmm_el = 0.0_dp

      scf_section => section_vals_get_subs_vals(qs_env%input, "DFT%SCF")
      nspins = dft_control%nspins
      CPASSERT(ASSOCIATED(matrix_h))
      CPASSERT(ASSOCIATED(rho))
      CPASSERT(SIZE(ks_matrix) > 0)

      DO ispin = 1, nspins
         DO img = 1, SIZE(ks_matrix, 2)
            ! copy the core matrix into the fock matrix
            CALL dbcsr_copy(ks_matrix(ispin, img)%matrix, matrix_h(1, img)%matrix)
         END DO
      END DO

      IF (dft_control%qs_control%dftb_control%self_consistent) THEN
         ! Mulliken charges
         CALL get_qs_env(qs_env=qs_env, particle_set=particle_set, &
                         matrix_s_kp=matrix_s)
         CALL qs_rho_get(rho, rho_ao_kp=matrix_p)
         natom = SIZE(particle_set)
         ALLOCATE (charges(natom, nspins))
         !
         CALL mulliken_charges(matrix_p, matrix_s, para_env, charges)
         !
         ALLOCATE (mcharge(natom))
         nkind = SIZE(atomic_kind_set)
         DO ikind = 1, nkind
            CALL get_atomic_kind(atomic_kind_set(ikind), natom=natom)
            CALL get_qs_kind(qs_kind_set(ikind), dftb_parameter=dftb_kind)
            CALL get_dftb_atom_param(dftb_kind, zeff=zeff)
            DO iatom = 1, natom
               atom_a = atomic_kind_set(ikind)%atom_list(iatom)
               mcharge(atom_a) = zeff - SUM(charges(atom_a, 1:nspins))
            END DO
         END DO
         DEALLOCATE (charges)

         CALL build_dftb_coulomb(qs_env, ks_matrix, rho, mcharge, energy, &
                                 calculate_forces, just_energy)

         CALL efield_tb_matrix(qs_env, ks_matrix, rho, mcharge, energy, &
                               calculate_forces, just_energy)

         DEALLOCATE (mcharge)

      END IF

      IF (qs_env%qmmm) THEN
         CPASSERT(SIZE(ks_matrix, 2) == 1)
         DO ispin = 1, nspins
            ! If QM/MM sumup the 1el Hamiltonian
            CALL dbcsr_add(ks_matrix(ispin, 1)%matrix, qs_env%ks_qmmm_env%matrix_h(1)%matrix, &
                           1.0_dp, 1.0_dp)
            CALL qs_rho_get(rho, rho_ao=matrix_p1)
            ! Compute QM/MM Energy
            CALL dbcsr_dot(qs_env%ks_qmmm_env%matrix_h(1)%matrix, &
                           matrix_p1(ispin)%matrix, qmmm_el)
            energy%qmmm_el = energy%qmmm_el + qmmm_el
         END DO
         pc_ener = qs_env%ks_qmmm_env%pc_ener
         energy%qmmm_el = energy%qmmm_el + pc_ener
      END IF

      energy%total = energy%core + energy%hartree + energy%qmmm_el + energy%efield + &
                     energy%repulsive + energy%dispersion + energy%dftb3

      output_unit = cp_print_key_unit_nr(logger, scf_section, "PRINT%DETAILED_ENERGY", &
                                         extension=".scfLog")
      IF (output_unit > 0) THEN
         WRITE (UNIT=output_unit, FMT="(/,(T9,A,T60,F20.10))") &
            "Repulsive pair potential energy:               ", energy%repulsive, &
            "Zeroth order Hamiltonian energy:               ", energy%core, &
            "Charge fluctuation energy:                     ", energy%hartree, &
            "London dispersion energy:                      ", energy%dispersion
         IF (ABS(energy%efield) > 1.e-12_dp) THEN
            WRITE (UNIT=output_unit, FMT="(T9,A,T60,F20.10)") &
               "Electric field interaction energy:             ", energy%efield
         END IF
         IF (dft_control%qs_control%dftb_control%dftb3_diagonal) THEN
            WRITE (UNIT=output_unit, FMT="(T9,A,T60,F20.10)") &
               "DFTB3 3rd Order Energy Correction              ", energy%dftb3
         END IF
         IF (qs_env%qmmm) THEN
            WRITE (UNIT=output_unit, FMT="(T9,A,T60,F20.10)") &
               "QM/MM Electrostatic energy:                    ", energy%qmmm_el
         END IF
      END IF
      CALL cp_print_key_finished_output(output_unit, logger, scf_section, &
                                        "PRINT%DETAILED_ENERGY")
      ! here we compute dE/dC if needed. Assumes dE/dC is H_{ks}C (plus occupation numbers)
      IF (qs_env%requires_mo_derivs .AND. .NOT. just_energy) THEN
         CPASSERT(SIZE(ks_matrix, 2) == 1)
         CALL get_qs_env(qs_env, mo_derivs=mo_derivs, mos=mo_array)
         DO ispin = 1, SIZE(mo_derivs)
            CALL get_mo_set(mo_set=mo_array(ispin)%mo_set, &
                            mo_coeff_b=mo_coeff, occupation_numbers=occupation_numbers)
            IF (.NOT. mo_array(ispin)%mo_set%use_mo_coeff_b) THEN
               CPABORT("")
            ENDIF
            CALL dbcsr_multiply('n', 'n', 1.0_dp, ks_matrix(ispin, 1)%matrix, mo_coeff, &
                                0.0_dp, mo_derivs(ispin)%matrix)
         ENDDO
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE build_dftb_ks_matrix

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param nderivative ...
!> \param matrix_s ...
! **************************************************************************************************
   SUBROUTINE build_dftb_overlap(qs_env, nderivative, matrix_s)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: nderivative
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_dftb_overlap'

      INTEGER :: handle, i, iatom, icol, ikind, indder, irow, j, jatom, jkind, l1, l2, la, lb, &
         llm, lmaxi, lmaxj, m, n1, n2, natom, natorb_a, natorb_b, ngrd, ngrdcut, nkind
      LOGICAL                                            :: defined, found
      REAL(KIND=dp)                                      :: ddr, dgrd, dr, f0
      REAL(KIND=dp), DIMENSION(0:3)                      :: skself
      REAL(KIND=dp), DIMENSION(3)                        :: drij, rij
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: dsblock, dsblockm, sblock, smatij, smatji
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: dsblock1
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(block_p_type), DIMENSION(2:10)                :: dsblocks
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(dftb_control_type), POINTER                   :: dftb_control
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(qs_dftb_atom_type), POINTER                   :: dftb_kind_a, dftb_kind_b
      TYPE(qs_dftb_pairpot_type), DIMENSION(:, :), &
         POINTER                                         :: dftb_potential
      TYPE(qs_dftb_pairpot_type), POINTER                :: dftb_param_ij, dftb_param_ji
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      ! set pointers
      iptr = 0
      DO la = 0, 3
         DO lb = 0, 3
            llm = 0
            DO l1 = 0, MAX(la, lb)
               DO l2 = 0, MIN(l1, la, lb)
                  DO m = 0, l2
                     llm = llm + 1
                     iptr(l1, l2, m, la, lb) = llm
                  END DO
               END DO
            END DO
         END DO
      END DO

      NULLIFY (logger)
      logger => cp_get_default_logger()

      NULLIFY (atomic_kind_set, qs_kind_set, sab_orb)

      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set, &
                      dft_control=dft_control)

      dftb_control => dft_control%qs_control%dftb_control

      NULLIFY (dftb_potential)
      CALL get_qs_env(qs_env=qs_env, &
                      dftb_potential=dftb_potential)

      nkind = SIZE(atomic_kind_set)

      ! Allocate the overlap matrix
      CALL get_qs_env(qs_env=qs_env, sab_orb=sab_orb)
      CALL setup_matrices1(qs_env, nderivative, matrix_s, 'OVERLAP', sab_orb)

      CALL neighbor_list_iterator_create(nl_iterator, sab_orb)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, &
                                iatom=iatom, jatom=jatom, r=rij)

         CALL get_atomic_kind(atomic_kind_set(ikind), natom=natom)
         CALL get_qs_kind(qs_kind_set(ikind), dftb_parameter=dftb_kind_a)
         CALL get_dftb_atom_param(dftb_kind_a, &
                                  defined=defined, lmax=lmaxi, skself=skself, &
                                  natorb=natorb_a)

         IF (.NOT. defined .OR. natorb_a < 1) CYCLE

         CALL get_qs_kind(qs_kind_set(jkind), dftb_parameter=dftb_kind_b)
         CALL get_dftb_atom_param(dftb_kind_b, &
                                  defined=defined, lmax=lmaxj, natorb=natorb_b)

         IF (.NOT. defined .OR. natorb_b < 1) CYCLE

         ! retrieve information on F and S matrix
         dftb_param_ij => dftb_potential(ikind, jkind)
         dftb_param_ji => dftb_potential(jkind, ikind)
         ! assume table size and type is symmetric
         ngrd = dftb_param_ij%ngrd
         ngrdcut = dftb_param_ij%ngrdcut
         dgrd = dftb_param_ij%dgrd
         ddr = dgrd*0.1_dp
         CPASSERT(dftb_param_ij%llm == dftb_param_ji%llm)
         llm = dftb_param_ij%llm
         smatij => dftb_param_ij%smat
         smatji => dftb_param_ji%smat

         dr = SQRT(SUM(rij(:)**2))
         IF (NINT(dr/dgrd) <= ngrdcut) THEN

            icol = MAX(iatom, jatom); irow = MIN(iatom, jatom)

            NULLIFY (sblock)
            CALL dbcsr_get_block_p(matrix=matrix_s(1)%matrix, &
                                   row=irow, col=icol, BLOCK=sblock, found=found)
            CPASSERT(found)

            IF (nderivative .GT. 0) THEN
               DO i = 2, SIZE(matrix_s, 1)
                  NULLIFY (dsblocks(i)%block)
                  CALL dbcsr_get_block_p(matrix=matrix_s(i)%matrix, &
                                         row=irow, col=icol, BLOCK=dsblocks(i)%block, found=found)
               END DO
            END IF

            IF (iatom == jatom .AND. dr < 0.001_dp) THEN
               ! diagonal block
               DO i = 1, natorb_a
                  sblock(i, i) = sblock(i, i) + 1._dp
               END DO
            ELSE
               ! off-diagonal block
               CALL compute_block_sk(sblock, smatij, smatji, rij, ngrd, ngrdcut, dgrd, &
                                     llm, lmaxi, lmaxj, irow, iatom)

               IF (nderivative .GE. 1) THEN
                  n1 = SIZE(sblock, 1); n2 = SIZE(sblock, 2)
                  indder = 1 ! used to put the 2nd derivatives in the correct matric (5=xx,8=yy,10=zz)

                  ALLOCATE (dsblock1(n1, n2, 3), dsblock(n1, n2), dsblockm(n1, n2))
                  dsblock1 = 0.0_dp
                  DO i = 1, 3
                     dsblock = 0._dp; dsblockm = 0.0_dp
                     drij = rij
                     f0 = 1.0_dp; IF (irow == iatom) f0 = -1.0_dp

                     drij(i) = rij(i) - ddr*f0
                     CALL compute_block_sk(dsblockm, smatij, smatji, drij, ngrd, ngrdcut, dgrd, &
                                           llm, lmaxi, lmaxj, irow, iatom)

                     drij(i) = rij(i) + ddr*f0
                     CALL compute_block_sk(dsblock, smatij, smatji, drij, ngrd, ngrdcut, dgrd, &
                                           llm, lmaxi, lmaxj, irow, iatom)

                     dsblock1(:, :, i) = (dsblock + dsblockm)
                     dsblock = dsblock - dsblockm
                     dsblock = dsblock/(2.0_dp*ddr)

                     CPASSERT(ASSOCIATED(dsblocks(i + 1)%block))
                     dsblocks(i + 1)%block = dsblocks(i + 1)%block + dsblock
                     IF (nderivative .GT. 1) THEN
                        indder = indder + 5 - i
                        dsblocks(indder)%block = 0.0_dp
                        dsblocks(indder)%block = dsblocks(indder)%block + &
                                                 (dsblock1(:, :, i) - 2.0_dp*sblock)/ddr**2
                     END IF
                  ENDDO

                  IF (nderivative .GT. 1) THEN
                     DO i = 1, 2
                        DO j = i + 1, 3
                           dsblock = 0._dp; dsblockm = 0.0_dp
                           drij = rij
                           f0 = 1.0_dp; IF (irow == iatom) f0 = -1.0_dp

                           drij(i) = rij(i) - ddr*f0; drij(j) = rij(j) - ddr*f0
                           CALL compute_block_sk(dsblockm, smatij, smatji, drij, ngrd, ngrdcut, dgrd, &
                                                 llm, lmaxi, lmaxj, irow, iatom)

                           drij(i) = rij(i) + ddr*f0; drij(j) = rij(j) + ddr*f0
                           CALL compute_block_sk(dsblock, smatij, smatji, drij, ngrd, ngrdcut, dgrd, &
                                                 llm, lmaxi, lmaxj, irow, iatom)

                           indder = 2 + 2*i + j
                           dsblocks(indder)%block = 0.0_dp
                           dsblocks(indder)%block = &
                              dsblocks(indder)%block + ( &
                              dsblock + dsblockm - dsblock1(:, :, i) - dsblock1(:, :, j) + 2.0_dp*sblock)/(2.0_dp*ddr**2)
                        END DO
                     END DO
                  END IF

                  DEALLOCATE (dsblock1, dsblock, dsblockm)
               END IF
            END IF
         END IF
      END DO
      CALL neighbor_list_iterator_release(nl_iterator)

      DO i = 1, SIZE(matrix_s, 1)
         CALL dbcsr_finalize(matrix_s(i)%matrix)
      ENDDO

      CALL timestop(handle)

   END SUBROUTINE build_dftb_overlap

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param nderivative ...
!> \param matrices ...
!> \param mnames ...
!> \param sab_nl ...
! **************************************************************************************************
   SUBROUTINE setup_matrices1(qs_env, nderivative, matrices, mnames, sab_nl)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: nderivative
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrices
      CHARACTER(LEN=*)                                   :: mnames
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl

      CHARACTER(1)                                       :: symmetry_type
      CHARACTER(LEN=default_string_length)               :: matnames
      INTEGER                                            :: i, natom, neighbor_list_id, nkind, nmat, &
                                                            nsgf
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: first_sgf, last_sgf
      INTEGER, DIMENSION(:), POINTER                     :: row_blk_sizes
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (particle_set, atomic_kind_set)

      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      particle_set=particle_set, &
                      dbcsr_dist=dbcsr_dist, &
                      neighbor_list_id=neighbor_list_id)

      nkind = SIZE(atomic_kind_set)
      natom = SIZE(particle_set)

      CALL get_qs_kind_set(qs_kind_set, nsgf=nsgf)

      ALLOCATE (first_sgf(natom))
      ALLOCATE (last_sgf(natom))

      CALL get_particle_set(particle_set, qs_kind_set, &
                            first_sgf=first_sgf, &
                            last_sgf=last_sgf)

      nmat = 0
      IF (nderivative == 0) nmat = 1
      IF (nderivative == 1) nmat = 4
      IF (nderivative == 2) nmat = 10
      CPASSERT(nmat > 0)

      ALLOCATE (row_blk_sizes(natom))
      CALL dbcsr_convert_offsets_to_sizes(first_sgf, row_blk_sizes, last_sgf)

      CALL dbcsr_allocate_matrix_set(matrices, nmat)

      ! Up to 2nd derivative take care to get the symmetries correct
      DO i = 1, nmat
         IF (i .GT. 1) THEN
            matnames = TRIM(mnames)//" DERIVATIVE MATRIX DFTB"
            symmetry_type = dbcsr_type_antisymmetric
            IF (i .GT. 4) symmetry_type = dbcsr_type_symmetric
         ELSE
            symmetry_type = dbcsr_type_symmetric
            matnames = TRIM(mnames)//" MATRIX DFTB"
         END IF
         ALLOCATE (matrices(i)%matrix)
         CALL dbcsr_create(matrix=matrices(i)%matrix, &
                           name=TRIM(matnames), &
                           dist=dbcsr_dist, matrix_type=symmetry_type, &
                           row_blk_size=row_blk_sizes, col_blk_size=row_blk_sizes, &
                           nze=0, mutable_work=.TRUE.)
         CALL cp_dbcsr_alloc_block_from_nbl(matrices(i)%matrix, sab_nl)
      END DO

      DEALLOCATE (first_sgf)
      DEALLOCATE (last_sgf)

      DEALLOCATE (row_blk_sizes)

   END SUBROUTINE setup_matrices1

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param nderivative ...
!> \param nimg ...
!> \param matrices ...
!> \param mnames ...
!> \param sab_nl ...
! **************************************************************************************************
   SUBROUTINE setup_matrices2(qs_env, nderivative, nimg, matrices, mnames, sab_nl)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: nderivative, nimg
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrices
      CHARACTER(LEN=*)                                   :: mnames
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl

      CHARACTER(1)                                       :: symmetry_type
      CHARACTER(LEN=default_string_length)               :: matnames
      INTEGER                                            :: i, img, natom, neighbor_list_id, nkind, &
                                                            nmat, nsgf
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: first_sgf, last_sgf
      INTEGER, DIMENSION(:), POINTER                     :: row_blk_sizes
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      NULLIFY (particle_set, atomic_kind_set)

      CALL get_qs_env(qs_env=qs_env, &
                      atomic_kind_set=atomic_kind_set, &
                      qs_kind_set=qs_kind_set, &
                      particle_set=particle_set, &
                      dbcsr_dist=dbcsr_dist, &
                      neighbor_list_id=neighbor_list_id)

      nkind = SIZE(atomic_kind_set)
      natom = SIZE(particle_set)

      CALL get_qs_kind_set(qs_kind_set, nsgf=nsgf)

      ALLOCATE (first_sgf(natom))
      ALLOCATE (last_sgf(natom))

      CALL get_particle_set(particle_set, qs_kind_set, &
                            first_sgf=first_sgf, &
                            last_sgf=last_sgf)

      nmat = 0
      IF (nderivative == 0) nmat = 1
      IF (nderivative == 1) nmat = 4
      IF (nderivative == 2) nmat = 10
      CPASSERT(nmat > 0)

      ALLOCATE (row_blk_sizes(natom))
      CALL dbcsr_convert_offsets_to_sizes(first_sgf, row_blk_sizes, last_sgf)

      CALL dbcsr_allocate_matrix_set(matrices, nmat, nimg)

      ! Up to 2nd derivative take care to get the symmetries correct
      DO img = 1, nimg
         DO i = 1, nmat
            IF (i .GT. 1) THEN
               matnames = TRIM(mnames)//" DERIVATIVE MATRIX DFTB"
               symmetry_type = dbcsr_type_antisymmetric
               IF (i .GT. 4) symmetry_type = dbcsr_type_symmetric
            ELSE
               symmetry_type = dbcsr_type_symmetric
               matnames = TRIM(mnames)//" MATRIX DFTB"
            END IF
            ALLOCATE (matrices(i, img)%matrix)
            CALL dbcsr_create(matrix=matrices(i, img)%matrix, &
                              name=TRIM(matnames), &
                              dist=dbcsr_dist, matrix_type=symmetry_type, &
                              row_blk_size=row_blk_sizes, col_blk_size=row_blk_sizes, &
                              nze=0, mutable_work=.TRUE.)
            CALL cp_dbcsr_alloc_block_from_nbl(matrices(i, img)%matrix, sab_nl)
         END DO
      END DO

      DEALLOCATE (first_sgf)
      DEALLOCATE (last_sgf)

      DEALLOCATE (row_blk_sizes)

   END SUBROUTINE setup_matrices2

END MODULE qs_dftb_matrices

