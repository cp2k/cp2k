!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief K-point MO wavefunction dump to TEXT file for post-processing (PDOS, etc.)
!> \par History
!>       2026.02 created
! **************************************************************************************************

MODULE kpoint_mo_dump

   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_create, dbcsr_deallocate_matrix, dbcsr_desymmetrize, dbcsr_p_type, dbcsr_set, &
        dbcsr_type, dbcsr_type_antisymmetric, dbcsr_type_no_symmetry, dbcsr_type_symmetric
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_submatrix,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE kpoint_methods,                  ONLY: rskp_transform
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_env_type,&
                                              kpoint_type
   USE message_passing,                 ONLY: mp_para_env_type
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: angstrom
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: write_kpoint_mo_data

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'kpoint_mo_dump'

CONTAINS

! **************************************************************************************************
!> \brief Write k-point resolved MO data to formatted text file.
!>
!>  The output file contains, for each k-point:
!>    - Eigenvalues and occupations
!>    - MO coefficient matrix C(k) (real and imaginary parts)
!>    - Overlap matrix S(k) (real and imaginary parts)
!>  Plus a header with cell, atom, and basis set angular momentum info.
!>
!>  Parallel strategy:
!>    MO coefficients:  cp_fm_get_submatrix (collective on kp-group),
!>                      then para_env_inter_kp%sum to collect across groups.
!>    S(k):             Built on global communicator using rskp_transform,
!>                      then copy_dbcsr_to_fm + cp_fm_get_submatrix.
!>    File write:       Only on ionode (para_env_global%is_source()).
!>
!> \param qs_env  the QS environment (after converged SCF with k-points)
! **************************************************************************************************
   SUBROUTINE write_kpoint_mo_data(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'write_kpoint_mo_data'

      CHARACTER(LEN=2)                                   :: element_symbol
      CHARACTER(LEN=default_string_length)               :: filename, project_name
      INTEGER :: handle, i, iatom, ikind, ikp, ikp_loc, iset, isgf_global, ishell, ispin, iw, j, &
         lshell, n, nao, natom, nkp, nmo, nset_atom, nspins, output_unit
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ao_l, first_sgf, last_sgf
      INTEGER, DIMENSION(2)                              :: kp_range
      INTEGER, DIMENSION(:), POINTER                     :: nshell_set
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgf_set, l_set, last_sgf_set
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: ionode, use_real_wfn
      REAL(KIND=dp)                                      :: zeff
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eval_buf, occ_buf
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: Cim_buf, Cre_buf, Sim_buf, Sre_buf
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues, occupation_numbers, wkp
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: xkp
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env_global
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct_S
      TYPE(cp_fm_type)                                   :: fm_S_global
      TYPE(cp_fm_type), POINTER                          :: mo_coeff_im, mo_coeff_re
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s
      TYPE(dbcsr_type), POINTER                          :: cmatrix, rmatrix, tmpmat
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(kpoint_type), POINTER                         :: kpoints
      TYPE(mo_set_type), DIMENSION(:, :), POINTER        :: mos_kp
      TYPE(mp_para_env_type), POINTER                    :: para_env_global, para_env_inter_kp
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      output_unit = cp_logger_get_default_io_unit(logger)

      ! ================================================================
      ! Gather all required data from qs_env
      ! ================================================================
      NULLIFY (cell, dft_control, kpoints, particle_set, qs_kind_set, matrix_s, &
               sab_nl, para_env_global, blacs_env_global)

      CALL get_qs_env(qs_env, &
                      cell=cell, &
                      dft_control=dft_control, &
                      kpoints=kpoints, &
                      particle_set=particle_set, &
                      qs_kind_set=qs_kind_set, &
                      matrix_s_kp=matrix_s, &
                      para_env=para_env_global, &
                      blacs_env=blacs_env_global)

      nspins = dft_control%nspins
      natom = SIZE(particle_set)
      ionode = para_env_global%is_source()

      CALL get_kpoint_info(kpoints, nkp=nkp, xkp=xkp, wkp=wkp, &
                           use_real_wfn=use_real_wfn, kp_range=kp_range, &
                           sab_nl=sab_nl, cell_to_index=cell_to_index, &
                           para_env_inter_kp=para_env_inter_kp)
      CPASSERT(ASSOCIATED(sab_nl))

      ! Total number of AOs and MOs
      CALL get_qs_kind_set(qs_kind_set, nsgf=nao)

      nmo = 0
      IF (SIZE(kpoints%kp_env) > 0) THEN
         mos_kp => kpoints%kp_env(1)%kpoint_env%mos
         CALL get_mo_set(mo_set=mos_kp(1, 1), nmo=nmo)
      END IF
      CALL para_env_inter_kp%max(nmo)

      IF (output_unit > 0) THEN
         WRITE (output_unit, '(/,T3,A)') "KPOINT_MO_DUMP| Writing k-point wavefunction data"
         WRITE (output_unit, '(T3,A,I6,A,I6,A,I6,A,I4)') &
            "KPOINT_MO_DUMP| nao=", nao, " nmo=", nmo, " nkp=", nkp, " nspins=", nspins
      END IF

      ! ================================================================
      ! Build atom-to-AO mapping with angular momentum labels
      ! ================================================================
      ALLOCATE (first_sgf(natom), last_sgf(natom), ao_l(nao))
      CALL get_particle_set(particle_set, qs_kind_set, &
                            first_sgf=first_sgf, last_sgf=last_sgf)
      ao_l(:) = -1
      DO iatom = 1, natom
         ikind = particle_set(iatom)%atomic_kind%kind_number
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set)
         CALL get_gto_basis_set(orb_basis_set, nset=nset_atom, &
                                nshell=nshell_set, l=l_set, &
                                first_sgf=first_sgf_set, last_sgf=last_sgf_set)
         DO iset = 1, nset_atom
            DO ishell = 1, nshell_set(iset)
               lshell = l_set(ishell, iset)
               DO i = first_sgf_set(ishell, iset), last_sgf_set(ishell, iset)
                  isgf_global = first_sgf(iatom) - 1 + i
                  ao_l(isgf_global) = lshell
               END DO
            END DO
         END DO
      END DO

      ! ================================================================
      ! Allocate work buffers
      ! ================================================================
      ALLOCATE (Cre_buf(nao, nmo), Cim_buf(nao, nmo))
      ALLOCATE (Sre_buf(nao, nao), Sim_buf(nao, nao))
      ALLOCATE (eval_buf(nmo), occ_buf(nmo))

      ! Create dbcsr work matrices for S(k) construction (global communicator)
      ALLOCATE (rmatrix, cmatrix, tmpmat)
      CALL dbcsr_create(rmatrix, template=matrix_s(1, 1)%matrix, &
                        matrix_type=dbcsr_type_symmetric)
      CALL dbcsr_create(cmatrix, template=matrix_s(1, 1)%matrix, &
                        matrix_type=dbcsr_type_antisymmetric)
      CALL dbcsr_create(tmpmat, template=matrix_s(1, 1)%matrix, &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL cp_dbcsr_alloc_block_from_nbl(rmatrix, sab_nl)
      CALL cp_dbcsr_alloc_block_from_nbl(cmatrix, sab_nl)

      ! Create a global FM for S(k): nao x nao on the global BLACS grid
      NULLIFY (fm_struct_S)
      CALL cp_fm_struct_create(fm_struct_S, context=blacs_env_global, &
                               nrow_global=nao, ncol_global=nao, &
                               para_env=para_env_global)
      CALL cp_fm_create(fm_S_global, fm_struct_S, name="S(k) work")
      CALL cp_fm_struct_release(fm_struct_S)

      ! ================================================================
      ! Open output file (ionode only)
      ! ================================================================
      iw = -1
      project_name = logger%iter_info%project_name
      IF (ionode) THEN
         filename = TRIM(project_name)//"-1.mokp"
         CALL open_file(file_name=TRIM(filename), file_status="REPLACE", &
                        file_action="WRITE", file_form="FORMATTED", &
                        unit_number=iw)
      END IF

      ! ================================================================
      ! WRITE HEADER
      ! ================================================================
      IF (iw > 0) THEN
         WRITE (iw, '(A)') "# CP2K_KPOINT_MO_DUMP"
         WRITE (iw, '(A)') "# All energies in Hartree, lengths in Angstrom, k-points in fractional coords"

         ! Cell vectors: hmat(:,i) = i-th lattice vector, written as rows
         WRITE (iw, '(A)') "# CELL_VECTORS [Angstrom]"
         WRITE (iw, '(3ES24.15)') &
            cell%hmat(1, 1)*angstrom, &
            cell%hmat(2, 1)*angstrom, &
            cell%hmat(3, 1)*angstrom
         WRITE (iw, '(3ES24.15)') &
            cell%hmat(1, 2)*angstrom, &
            cell%hmat(2, 2)*angstrom, &
            cell%hmat(3, 2)*angstrom
         WRITE (iw, '(3ES24.15)') &
            cell%hmat(1, 3)*angstrom, &
            cell%hmat(2, 3)*angstrom, &
            cell%hmat(3, 3)*angstrom

         WRITE (iw, '(A,4I8)') "# DIMENSIONS: natom nspins nao nkp =", natom, nspins, nao, nkp
         WRITE (iw, '(A,I8)') "# NMO =", nmo
         WRITE (iw, '(A,L2)') "# USE_REAL_WFN =", use_real_wfn

         ! Atom table
         WRITE (iw, '(A)') &
            "# ATOM_LIST: Atom_ID  Element  Z_eff  X [Ang]  Y [Ang]  Z [Ang]  First_AO  Last_AO  (1-based, inclusive)"
         DO iatom = 1, natom
            ikind = particle_set(iatom)%atomic_kind%kind_number
            CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, &
                                 element_symbol=element_symbol)
            CALL get_qs_kind(qs_kind_set(ikind), zeff=zeff)
            WRITE (iw, '(I6,2X,A2,I6,3F15.8,2I10)') &
               iatom, element_symbol, NINT(zeff), particle_set(iatom)%r(:)*angstrom, first_sgf(iatom), last_sgf(iatom)
         END DO

         ! AO angular momentum map
         WRITE (iw, '(A)') "# AO_L_MAP: iao  l  (1-based)"
         DO i = 1, nao
            WRITE (iw, '(2I6)') i, ao_l(i)
         END DO

         ! K-point list
         WRITE (iw, '(A)') "# KPOINT_LIST: ikp  kx  ky  kz  weight"
         DO ikp = 1, nkp
            WRITE (iw, '(I6,4ES20.12)') ikp, xkp(1:3, ikp), wkp(ikp)
         END DO
      END IF

      ! ================================================================
      ! WRITE PER-KPOINT DATA
      ! ================================================================
      DO ikp = 1, nkp

         ! ==============================================================
         ! A) MO coefficients, eigenvalues, occupations (per spin)
         ! ==============================================================
         DO ispin = 1, nspins

            Cre_buf(:, :) = 0.0_dp
            Cim_buf(:, :) = 0.0_dp
            eval_buf(:) = 0.0_dp
            occ_buf(:) = 0.0_dp

            ! Only the owning kpoint group fills the buffers
            IF (ikp >= kp_range(1) .AND. ikp <= kp_range(2)) THEN
               ikp_loc = ikp - kp_range(1) + 1
               kp => kpoints%kp_env(ikp_loc)%kpoint_env
               mos_kp => kp%mos

               ! Eigenvalues and occupations
               CALL get_mo_set(mos_kp(1, ispin), &
                               eigenvalues=eigenvalues, &
                               occupation_numbers=occupation_numbers)
               eval_buf(1:nmo) = eigenvalues(1:nmo)
               occ_buf(1:nmo) = occupation_numbers(1:nmo)

               ! MO coefficients real part
               ! cp_fm_get_submatrix is collective on the FM's BLACS communicator
               CALL get_mo_set(mos_kp(1, ispin), mo_coeff=mo_coeff_re)
               CALL cp_fm_get_submatrix(mo_coeff_re, Cre_buf)

               ! MO coefficients imaginary part
               IF (.NOT. use_real_wfn) THEN
                  CALL get_mo_set(mos_kp(2, ispin), mo_coeff=mo_coeff_im)
                  CALL cp_fm_get_submatrix(mo_coeff_im, Cim_buf)
               END IF
            END IF

            ! Collect across k-point groups (non-owning groups contribute zeros)
            CALL para_env_inter_kp%sum(eval_buf)
            CALL para_env_inter_kp%sum(occ_buf)
            CALL para_env_inter_kp%sum(Cre_buf)
            IF (.NOT. use_real_wfn) THEN
               CALL para_env_inter_kp%sum(Cim_buf)
            END IF

            ! Write to file
            IF (iw > 0) THEN
               WRITE (iw, '(A,2I6)') "# BEGIN_KPOINT_SPIN ikp ispin =", ikp, ispin

               WRITE (iw, '(A)') "# EIGENVALUES"
               WRITE (iw, '(5ES25.15E3)') (eval_buf(n), n=1, nmo)

               WRITE (iw, '(A)') "# OCCUPATIONS"
               WRITE (iw, '(5ES25.15E3)') (occ_buf(n), n=1, nmo)

               WRITE (iw, '(A)') "# MO_COEFF_RE"
               DO i = 1, nao
                  WRITE (iw, '(5ES25.15E3)') (Cre_buf(i, j), j=1, nmo)
               END DO

               IF (.NOT. use_real_wfn) THEN
                  WRITE (iw, '(A)') "# MO_COEFF_IM"
                  DO i = 1, nao
                     WRITE (iw, '(5ES25.15E3)') (Cim_buf(i, j), j=1, nmo)
                  END DO
               END IF
               WRITE (iw, '(A,2I6)') "# END_KPOINT_SPIN ikp ispin =", ikp, ispin
            END IF

         END DO ! ispin

         ! ==============================================================
         ! B) Overlap matrix S(k) (spin-independent, once per k-point)
         !    Constructed on the GLOBAL communicator so all ranks participate
         ! ==============================================================

         ! Build S(k) = sum_R exp(ik.R) S^R via Fourier transform
         CALL dbcsr_set(rmatrix, 0.0_dp)
         IF (use_real_wfn) THEN
            ! Real k-point (Gamma or time-reversal partner)
            CALL rskp_transform(rmatrix=rmatrix, rsmat=matrix_s, ispin=1, &
                                xkp=xkp(1:3, ikp), cell_to_index=cell_to_index, &
                                sab_nl=sab_nl)
            CALL dbcsr_desymmetrize(rmatrix, tmpmat)
            CALL copy_dbcsr_to_fm(tmpmat, fm_S_global)
            CALL cp_fm_get_submatrix(fm_S_global, Sre_buf)
            ! Sim_buf stays zero
            Sim_buf(:, :) = 0.0_dp
         ELSE
            ! General complex k-point
            CALL dbcsr_set(cmatrix, 0.0_dp)
            CALL rskp_transform(rmatrix=rmatrix, cmatrix=cmatrix, rsmat=matrix_s, &
                                ispin=1, xkp=xkp(1:3, ikp), &
                                cell_to_index=cell_to_index, sab_nl=sab_nl)
            ! Real part
            CALL dbcsr_desymmetrize(rmatrix, tmpmat)
            CALL copy_dbcsr_to_fm(tmpmat, fm_S_global)
            CALL cp_fm_get_submatrix(fm_S_global, Sre_buf)
            ! Imaginary part
            CALL dbcsr_desymmetrize(cmatrix, tmpmat)
            CALL copy_dbcsr_to_fm(tmpmat, fm_S_global)
            CALL cp_fm_get_submatrix(fm_S_global, Sim_buf)
         END IF
         ! NOTE: No inter-kp sum needed here because all ranks participate
         !       in rskp_transform + copy_dbcsr_to_fm + cp_fm_get_submatrix.

         ! Write S(k) to file
         IF (iw > 0) THEN
            WRITE (iw, '(A,I6)') "# BEGIN_OVERLAP ikp =", ikp
            WRITE (iw, '(A)') "# OVERLAP_RE"
            DO i = 1, nao
               WRITE (iw, '(5ES25.15E3)') (Sre_buf(i, j), j=1, nao)
            END DO
            IF (.NOT. use_real_wfn) THEN
               WRITE (iw, '(A)') "# OVERLAP_IM"
               DO i = 1, nao
                  WRITE (iw, '(5ES25.15E3)') (Sim_buf(i, j), j=1, nao)
               END DO
            END IF
            WRITE (iw, '(A,I6)') "# END_OVERLAP ikp =", ikp
         END IF

      END DO ! ikp

      ! ================================================================
      ! Finalize
      ! ================================================================
      IF (iw > 0) THEN
         WRITE (iw, '(A)') "# END_OF_FILE"
         CALL close_file(unit_number=iw)
      END IF

      IF (output_unit > 0) THEN
         WRITE (output_unit, '(T3,A,A)') &
            "KPOINT_MO_DUMP| Data written to ", TRIM(project_name)//"-1.mokp"
      END IF

      ! Release work arrays
      CALL cp_fm_release(fm_S_global)
      CALL dbcsr_deallocate_matrix(rmatrix)
      CALL dbcsr_deallocate_matrix(cmatrix)
      CALL dbcsr_deallocate_matrix(tmpmat)
      DEALLOCATE (Cre_buf, Cim_buf, Sre_buf, Sim_buf)
      DEALLOCATE (eval_buf, occ_buf)
      DEALLOCATE (first_sgf, last_sgf, ao_l)

      CALL timestop(handle)

   END SUBROUTINE write_kpoint_mo_data

END MODULE kpoint_mo_dump
