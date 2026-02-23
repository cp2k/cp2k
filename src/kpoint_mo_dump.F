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
   USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                              cp_print_key_unit_nr
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE kpoint_methods,                  ONLY: rskp_transform
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_env_type,&
                                              kpoint_type
   USE mathconstants,                   ONLY: pi
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

   ! Enum values for OVERLAP_EXPORT_TYPE keyword
   INTEGER, PARAMETER, PUBLIC :: mokp_overlap_gto = 1, &
                                 mokp_overlap_matrix = 2

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'kpoint_mo_dump'

CONTAINS

! **************************************************************************************************
!> \brief Write k-point resolved MO data to formatted text file.
!>
!>  The output file contains, for each k-point:
!>    - Eigenvalues and occupations
!>    - MO coefficient matrix C(k) (real and imaginary parts)
!>    - Overlap matrix S(k)  (if OVERLAP_EXPORT_TYPE = MATRIX)
!>  Plus a header with cell, atom info, and either GTO basis set definition
!>  (OVERLAP_EXPORT_TYPE = GTO) or AO angular momentum map (OVERLAP_EXPORT_TYPE = MATRIX).
!>
!>  Parallel strategy:
!>    MO coefficients:  cp_fm_get_submatrix (collective on kp-group),
!>                      then para_env_inter_kp%sum to collect across groups.
!>    S(k):             Built on global communicator using rskp_transform,
!>                      then copy_dbcsr_to_fm + cp_fm_get_submatrix.
!>    File write:       Only on ionode (cp_print_key_unit_nr returns -1 elsewhere).
!>
!> \param qs_env  the QS environment (after converged SCF with k-points)
!> \param print_section  the &MO_KP print key section (contains NDIGITS, OVERLAP_EXPORT_TYPE)
! **************************************************************************************************
   SUBROUTINE write_kpoint_mo_data(qs_env, print_section)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: print_section

      CHARACTER(len=*), PARAMETER :: routineN = 'write_kpoint_mo_data'
      CHARACTER(LEN=6), PARAMETER                        :: angmom = "spdfgh"

      CHARACTER(LEN=2)                                   :: element_symbol
      CHARACTER(LEN=20)                                  :: fmtstr_gto, fmtstr_sparse, fmtstr_vec
      INTEGER :: handle, i, iatom, ikind, ikp, ikp_loc, ipgf, iset, isgf_global, ishell, ispin, &
         iw, j, lshell, n, nao, natom, ndigits, nkp, nmo, nset_atom, nspins, output_unit, &
         overlap_data
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ao_l, first_sgf, last_sgf
      INTEGER, DIMENSION(2)                              :: kp_range
      INTEGER, DIMENSION(:), POINTER                     :: npgf_set, nshell_set
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgf_set, l_set, last_sgf_set
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: use_real_wfn, write_overlap
      REAL(KIND=dp)                                      :: expzet, prefac, thresh, zeff
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eval_buf, occ_buf
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: Cim_buf, Cre_buf, Sim_buf, Sre_buf
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues, occupation_numbers, wkp
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: xkp, zet
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: gcc
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
      ! Read keywords
      ! ================================================================
      CALL section_vals_val_get(print_section, "NDIGITS", i_val=ndigits)
      ndigits = MIN(MAX(3, ndigits), 30)
      CALL section_vals_val_get(print_section, "OVERLAP_EXPORT_TYPE", i_val=overlap_data)
      write_overlap = (overlap_data == mokp_overlap_matrix)

      ! Build format strings controlled by NDIGITS
      WRITE (UNIT=fmtstr_sparse, FMT='("(2I6,1X,ES",I0,".",I0,")")') ndigits + 10, ndigits
      WRITE (UNIT=fmtstr_vec, FMT='("(5ES",I0,".",I0,")")') ndigits + 10, ndigits
      WRITE (UNIT=fmtstr_gto, FMT='("(2ES",I0,".",I0,")")') ndigits + 10, ndigits
      thresh = 10.0_dp**(-ndigits)

      ! ================================================================
      ! Build atom-to-AO mapping
      ! ================================================================
      ALLOCATE (first_sgf(natom), last_sgf(natom))
      CALL get_particle_set(particle_set, qs_kind_set, &
                            first_sgf=first_sgf, last_sgf=last_sgf)

      ! AO angular momentum map: only needed when overlap is written explicitly,
      ! since in GTO mode l info is derivable from the basis set definition.
      IF (write_overlap) THEN
         ALLOCATE (ao_l(nao))
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
      END IF

      ! ================================================================
      ! Allocate work buffers
      ! ================================================================
      ALLOCATE (Cre_buf(nao, nmo), Cim_buf(nao, nmo))
      ALLOCATE (eval_buf(nmo), occ_buf(nmo))

      ! S(k) infrastructure: only needed in MATRIX mode
      IF (write_overlap) THEN
         ALLOCATE (Sre_buf(nao, nao), Sim_buf(nao, nao))

         ALLOCATE (rmatrix, cmatrix, tmpmat)
         CALL dbcsr_create(rmatrix, template=matrix_s(1, 1)%matrix, &
                           matrix_type=dbcsr_type_symmetric)
         CALL dbcsr_create(cmatrix, template=matrix_s(1, 1)%matrix, &
                           matrix_type=dbcsr_type_antisymmetric)
         CALL dbcsr_create(tmpmat, template=matrix_s(1, 1)%matrix, &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL cp_dbcsr_alloc_block_from_nbl(rmatrix, sab_nl)
         CALL cp_dbcsr_alloc_block_from_nbl(cmatrix, sab_nl)

         NULLIFY (fm_struct_S)
         CALL cp_fm_struct_create(fm_struct_S, context=blacs_env_global, &
                                  nrow_global=nao, ncol_global=nao, &
                                  para_env=para_env_global)
         CALL cp_fm_create(fm_S_global, fm_struct_S, name="S(k) work")
         CALL cp_fm_struct_release(fm_struct_S)
      END IF

      ! ================================================================
      ! Open output file via print key
      ! ================================================================
      iw = cp_print_key_unit_nr(logger, print_section, "", &
                                extension=".mokp", file_status="REPLACE")

      ! ================================================================
      ! WRITE HEADER
      ! ================================================================
      IF (iw > 0) THEN
         WRITE (iw, '(A)') "# CP2K_KPOINT_MO_DUMP, Version 1.0"
         IF (write_overlap) THEN
            WRITE (iw, '(A)') "# EXPORT_MODE: Coefficients + Overlap matrices"
         ELSE
            WRITE (iw, '(A)') "# EXPORT_MODE: GTO + Coefficients"
         END IF
         WRITE (iw, '(A)') "# All energies in Hartree, lengths in Angstrom, k-points in fractional coords"

         WRITE (iw, '(A)') "# CELL_VECTORS [Angstrom]"
         WRITE (iw, '(3ES16.6)') &
            cell%hmat(1, 1)*angstrom, cell%hmat(2, 1)*angstrom, cell%hmat(3, 1)*angstrom
         WRITE (iw, '(3ES16.6)') &
            cell%hmat(1, 2)*angstrom, cell%hmat(2, 2)*angstrom, cell%hmat(3, 2)*angstrom
         WRITE (iw, '(3ES16.6)') &
            cell%hmat(1, 3)*angstrom, cell%hmat(2, 3)*angstrom, cell%hmat(3, 3)*angstrom

         WRITE (iw, '(A,4I8)') "# DIMENSIONS: natom nspins nao nkp =", natom, nspins, nao, nkp
         WRITE (iw, '(A,I8)') "# NMO =", nmo
         WRITE (iw, '(A,L2)') "# USE_REAL_WFN =", use_real_wfn
         IF (write_overlap) THEN
            WRITE (iw, '(A)') "# OVERLAP_EXPORT_TYPE = MATRIX"
         ELSE
            WRITE (iw, '(A)') "# OVERLAP_EXPORT_TYPE = GTO"
         END IF

         ! Atom table
         WRITE (iw, '(A)') &
            "# ATOM_LIST: Atom_ID  Element  Z_eff  X [Ang]  Y [Ang]  Z [Ang]  First_AO  Last_AO"
         DO iatom = 1, natom
            ikind = particle_set(iatom)%atomic_kind%kind_number
            CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, &
                                 element_symbol=element_symbol)
            CALL get_qs_kind(qs_kind_set(ikind), zeff=zeff)
            WRITE (iw, '(I6,2X,A2,I6,3F15.6,2I10)') &
               iatom, element_symbol, NINT(zeff), particle_set(iatom)%r(:)*angstrom, &
               first_sgf(iatom), last_sgf(iatom)
         END DO

         ! K-point list
         WRITE (iw, '(A)') "# KPOINT_LIST: ikp  kx  ky  kz  weight"
         DO ikp = 1, nkp
            WRITE (iw, '(I6,4ES18.8)') ikp, xkp(1:3, ikp), wkp(ikp)
         END DO

         IF (write_overlap) THEN
            ! MATRIX mode: write AO angular momentum map (needed since no GTO info)
            WRITE (iw, '(A)') "# AO_L_MAP: iao  l  (1-based)"
            DO i = 1, nao
               WRITE (iw, '(2I6)') i, ao_l(i)
            END DO
         ELSE
            ! GTO mode: write basis set definition
            ! Contraction coefficients denormalized to MOLDEN convention
            ! (i.e. for raw unnormalized primitive Gaussians).
            ! Shell angular momentum ordering: spherical harmonics,
            ! CP2K convention m = -l, -l+1, ..., 0, ..., +l
            WRITE (iw, '(A)') "# GTO_BASIS (spherical, MOLDEN convention for contraction coefficients)"
            DO iatom = 1, natom
               ikind = particle_set(iatom)%atomic_kind%kind_number
               CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set)
               IF (ASSOCIATED(orb_basis_set)) THEN
                  WRITE (iw, '(I6,I4)') iatom, 0
                  CALL get_gto_basis_set(gto_basis_set=orb_basis_set, &
                                         nset=nset_atom, npgf=npgf_set, &
                                         nshell=nshell_set, l=l_set, &
                                         zet=zet, gcc=gcc)
                  DO iset = 1, nset_atom
                     DO ishell = 1, nshell_set(iset)
                        lshell = l_set(ishell, iset)
                        IF (lshell + 1 > LEN(angmom)) THEN
                           CPWARN("MOKP: Angular momentum l > 5 not supported in GTO output, skipping")
                           CYCLE
                        END IF
                        WRITE (iw, '(A2,I6,F8.2)') &
                           angmom(lshell + 1:lshell + 1)//" ", npgf_set(iset), 1.0_dp
                        ! Denormalize gcc: undo CP2K's internal normalization
                        ! (same formula as molden_utils.F, reverse of normalise_gcc_orb)
                        prefac = 2.0_dp**lshell*(2.0_dp/pi)**0.75_dp
                        expzet = 0.25_dp*(2*lshell + 3.0_dp)
                        DO ipgf = 1, npgf_set(iset)
                           WRITE (iw, fmtstr_gto) &
                              zet(ipgf, iset), &
                              gcc(ipgf, ishell, iset)/(prefac*zet(ipgf, iset)**expzet)
                        END DO
                     END DO
                  END DO
                  WRITE (iw, '(A)') ""
               END IF
            END DO
         END IF
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

            IF (ikp >= kp_range(1) .AND. ikp <= kp_range(2)) THEN
               ikp_loc = ikp - kp_range(1) + 1
               kp => kpoints%kp_env(ikp_loc)%kpoint_env
               mos_kp => kp%mos

               CALL get_mo_set(mos_kp(1, ispin), &
                               eigenvalues=eigenvalues, &
                               occupation_numbers=occupation_numbers)
               eval_buf(1:nmo) = eigenvalues(1:nmo)
               occ_buf(1:nmo) = occupation_numbers(1:nmo)

               CALL get_mo_set(mos_kp(1, ispin), mo_coeff=mo_coeff_re)
               CALL cp_fm_get_submatrix(mo_coeff_re, Cre_buf)

               IF (.NOT. use_real_wfn) THEN
                  CALL get_mo_set(mos_kp(2, ispin), mo_coeff=mo_coeff_im)
                  CALL cp_fm_get_submatrix(mo_coeff_im, Cim_buf)
               END IF
            END IF

            CALL para_env_inter_kp%sum(eval_buf)
            CALL para_env_inter_kp%sum(occ_buf)
            CALL para_env_inter_kp%sum(Cre_buf)
            IF (.NOT. use_real_wfn) THEN
               CALL para_env_inter_kp%sum(Cim_buf)
            END IF

            IF (iw > 0) THEN
               WRITE (iw, '(A,2I6)') "# BEGIN_KPOINT_SPIN ikp ispin =", ikp, ispin

               WRITE (iw, '(A)') "# EIGENVALUES"
               WRITE (iw, fmtstr_vec) (eval_buf(n), n=1, nmo)

               WRITE (iw, '(A)') "# OCCUPATIONS"
               WRITE (iw, fmtstr_vec) (occ_buf(n), n=1, nmo)

               WRITE (iw, '(A)') "# MO_COEFF_RE (Sparse: i_mo j_ao value)"
               DO i = 1, nmo
                  DO j = 1, nao
                     IF (ABS(Cre_buf(j, i)) >= thresh) THEN
                        WRITE (iw, fmtstr_sparse) i, j, Cre_buf(j, i)
                     END IF
                  END DO
               END DO

               IF (.NOT. use_real_wfn) THEN
                  WRITE (iw, '(A)') "# MO_COEFF_IM (Sparse: i_mo j_ao value)"
                  DO i = 1, nmo
                     DO j = 1, nao
                        IF (ABS(Cim_buf(j, i)) >= thresh) THEN
                           WRITE (iw, fmtstr_sparse) i, j, Cim_buf(j, i)
                        END IF
                     END DO
                  END DO
               END IF
               WRITE (iw, '(A,2I6)') "# END_KPOINT_SPIN ikp ispin =", ikp, ispin
            END IF

         END DO ! ispin

         ! ==============================================================
         ! B) Overlap matrix S(k) — only in MATRIX mode
         ! ==============================================================
         IF (write_overlap) THEN

            CALL dbcsr_set(rmatrix, 0.0_dp)
            IF (use_real_wfn) THEN
               CALL rskp_transform(rmatrix=rmatrix, rsmat=matrix_s, ispin=1, &
                                   xkp=xkp(1:3, ikp), cell_to_index=cell_to_index, &
                                   sab_nl=sab_nl)
               CALL dbcsr_desymmetrize(rmatrix, tmpmat)
               CALL copy_dbcsr_to_fm(tmpmat, fm_S_global)
               CALL cp_fm_get_submatrix(fm_S_global, Sre_buf)
               Sim_buf(:, :) = 0.0_dp
            ELSE
               CALL dbcsr_set(cmatrix, 0.0_dp)
               CALL rskp_transform(rmatrix=rmatrix, cmatrix=cmatrix, rsmat=matrix_s, &
                                   ispin=1, xkp=xkp(1:3, ikp), &
                                   cell_to_index=cell_to_index, sab_nl=sab_nl)
               CALL dbcsr_desymmetrize(rmatrix, tmpmat)
               CALL copy_dbcsr_to_fm(tmpmat, fm_S_global)
               CALL cp_fm_get_submatrix(fm_S_global, Sre_buf)
               CALL dbcsr_desymmetrize(cmatrix, tmpmat)
               CALL copy_dbcsr_to_fm(tmpmat, fm_S_global)
               CALL cp_fm_get_submatrix(fm_S_global, Sim_buf)
            END IF

            IF (iw > 0) THEN
               WRITE (iw, '(A,I6)') "# BEGIN_OVERLAP ikp =", ikp

               WRITE (iw, '(A)') "# OVERLAP_RE (Sparse: i_ao j_ao value)"
               DO i = 1, nao
                  DO j = 1, nao
                     IF (ABS(Sre_buf(i, j)) >= thresh) THEN
                        WRITE (iw, fmtstr_sparse) i, j, Sre_buf(i, j)
                     END IF
                  END DO
               END DO
               IF (.NOT. use_real_wfn) THEN
                  WRITE (iw, '(A)') "# OVERLAP_IM (Sparse: i_ao j_ao value)"
                  DO i = 1, nao
                     DO j = 1, nao
                        IF (ABS(Sim_buf(i, j)) >= thresh) THEN
                           WRITE (iw, fmtstr_sparse) i, j, Sim_buf(i, j)
                        END IF
                     END DO
                  END DO
               END IF
               WRITE (iw, '(A,I6)') "# END_OVERLAP ikp =", ikp
            END IF

         END IF ! write_overlap

      END DO ! ikp

      ! ================================================================
      ! Finalize
      ! ================================================================
      IF (iw > 0) THEN
         WRITE (iw, '(A)') "# END_OF_FILE"
      END IF

      CALL cp_print_key_finished_output(iw, logger, print_section, "")

      IF (output_unit > 0) THEN
         WRITE (output_unit, '(T3,A)') &
            "KPOINT_MO_DUMP| Data written to .mokp file"
      END IF

      ! Release work arrays
      IF (write_overlap) THEN
         CALL cp_fm_release(fm_S_global)
         CALL dbcsr_deallocate_matrix(rmatrix)
         CALL dbcsr_deallocate_matrix(cmatrix)
         CALL dbcsr_deallocate_matrix(tmpmat)
         DEALLOCATE (Sre_buf, Sim_buf)
         DEALLOCATE (ao_l)
      END IF
      DEALLOCATE (Cre_buf, Cim_buf)
      DEALLOCATE (eval_buf, occ_buf)
      DEALLOCATE (first_sgf, last_sgf)

      CALL timestop(handle)

   END SUBROUTINE write_kpoint_mo_data

END MODULE kpoint_mo_dump
