!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculation of band structures
!> \par History
!>       2015.06 created [JGH]
!> \author JGH
! **************************************************************************************************
MODULE qs_band_structure
   USE bibliography,                    ONLY: Quan2026
   USE cell_types,                      ONLY: cell_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_p_type
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_log_handling,                 ONLY: cp_logger_get_default_io_unit
   USE cp_parser_methods,               ONLY: read_float_object
   USE dft_plus_u,                      ONLY: lowdin_kp_prepare
   USE input_section_types,             ONLY: section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_string_length,&
                                              dp,&
                                              max_line_length
   USE kpoint_methods,                  ONLY: kpoint_env_initialize,&
                                              kpoint_init_cell_index,&
                                              kpoint_initialize_mo_set,&
                                              kpoint_initialize_mos
   USE kpoint_types,                    ONLY: get_kpoint_info,&
                                              kpoint_create,&
                                              kpoint_env_type,&
                                              kpoint_release,&
                                              kpoint_sym_create,&
                                              kpoint_type
   USE kpoint_unfolding,                ONLY: &
        kdifftol, mo_orthonorm_tol, unfolding_env_build_map, unfolding_env_create, &
        unfolding_env_release, unfolding_env_type, unfolding_fold_kpath, unfolding_kp_weights, &
        unfolding_unfolded_kspace
   USE machine,                         ONLY: m_walltime
   USE mathconstants,                   ONLY: twopi
   USE message_passing,                 ONLY: mp_para_env_type
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: angstrom,&
                                              evolt
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_env_release,&
                                              qs_environment_type
   USE qs_gamma2kp,                     ONLY: create_kp_from_gamma
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_mo_types,                     ONLY: get_mo_set
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_scf_diagonalization,          ONLY: diag_kp_smat,&
                                              do_general_diag_kp
   USE qs_scf_types,                    ONLY: qs_scf_env_type
   USE reference_manager,               ONLY: cite_reference
   USE scf_control_types,               ONLY: scf_control_type
   USE string_utilities,                ONLY: uppercase
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_band_structure'

   PUBLIC :: calculate_band_structure, calculate_kp_orbitals, calculate_kpoints_for_bs

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief Main routine for band structure calculation
!> \param qs_env ...
!> \author JGH
! **************************************************************************************************
   SUBROUTINE calculate_band_structure(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      LOGICAL                                            :: do_kpoints, explicit
      TYPE(section_vals_type), POINTER                   :: bs_input

      bs_input => section_vals_get_subs_vals(qs_env%input, "DFT%PRINT%BAND_STRUCTURE")
      CALL section_vals_get(bs_input, explicit=explicit)
      IF (explicit) THEN
         CALL get_qs_env(qs_env, do_kpoints=do_kpoints)
         IF (do_kpoints) THEN
            CALL do_calculate_band_structure(qs_env)
         ELSE
            BLOCK
               TYPE(qs_environment_type), POINTER :: qs_env_kp
               CALL create_kp_from_gamma(qs_env, qs_env_kp)
               CALL do_calculate_band_structure(qs_env_kp)
               CALL qs_env_release(qs_env_kp)
               DEALLOCATE (qs_env_kp)
            END BLOCK
         END IF
      END IF

   END SUBROUTINE calculate_band_structure

! **************************************************************************************************
!> \brief band structure calculation
!> \param qs_env ...
!> \author JGH
! **************************************************************************************************
   SUBROUTINE do_calculate_band_structure(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=32)                                  :: fmtstr
      CHARACTER(LEN=default_string_length)               :: filename, ustr
      CHARACTER(LEN=default_string_length), &
         DIMENSION(:), POINTER                           :: spname, strptr
      CHARACTER(LEN=max_line_length)                     :: error_message
      INTEGER :: bs_data_unit, g, i, i_rep, ia, ia_first, ia_last, ik, ikk, ikmap, ikpgr, imo, ip, &
         iq, iqc, ispin, j, n_atom_groups, n_ptr, n_rep, nadd, natom, nfold, nkp, nmo, nout, &
         npline, npoints, nspins, unit_nr
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: group_of_atom, path2K
      INTEGER, DIMENSION(2)                              :: kp_range
      INTEGER, DIMENSION(3, 3)                           :: transmat
      INTEGER, DIMENSION(:), POINTER                     :: atom_vals, mvals
      LOGICAL                                            :: do_unfold, explicit, io_default, &
                                                            lowdin_proj, my_kpgrp, stable_images
      REAL(KIND=dp)                                      :: sumdev, symm_tol, t1, t2
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: psub
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: fqs, kpfolded
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: wmat
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :, :)  :: wall
      REAL(KIND=dp), DIMENSION(3)                        :: dvec, kpptr
      REAL(KIND=dp), DIMENSION(3, 3)                     :: cart_hmat
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues, eigval, occnum, &
                                                            occupation_numbers, wkp
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: kpgeneral, kspecial, xkp
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_env_type), POINTER                     :: kp
      TYPE(kpoint_type), POINTER                         :: kpoint, kpoints
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(section_vals_type), POINTER                   :: bs_input, kpset, u_input
      TYPE(unfolding_env_type)                           :: uf

      bs_input => section_vals_get_subs_vals(qs_env%input, "DFT%PRINT%BAND_STRUCTURE")
      CALL section_vals_get(bs_input, explicit=explicit)
      IF (explicit) THEN
         CALL section_vals_val_get(bs_input, "FILE_NAME", c_val=filename)
         CALL section_vals_val_get(bs_input, "ADDED_MOS", i_val=nadd)
         unit_nr = cp_logger_get_default_io_unit()
         CALL get_qs_env(qs_env=qs_env, para_env=para_env)
         CALL get_qs_env(qs_env, cell=cell)
         cart_hmat(:, :) = cell%hmat(:, :)
         IF (cell%input_cell_canonicalized) cart_hmat(:, :) = cell%input_hmat(:, :)
         kpset => section_vals_get_subs_vals(bs_input, "KPOINT_SET")
         CALL section_vals_get(kpset, n_repetition=n_rep)
         u_input => section_vals_get_subs_vals(bs_input, "UNFOLDING")
         CALL section_vals_get(u_input, explicit=do_unfold)
         IF (do_unfold) THEN
            CALL cite_reference(Quan2026)
            NULLIFY (mvals)
            CALL section_vals_val_get(u_input, "TRANSFORMATION_MATRIX", i_vals=mvals)
            ! input lists the three matrix rows; RESHAPE fills columns first
            transmat = RESHAPE(mvals, [3, 3], ORDER=[2, 1])
            NULLIFY (mvals)  ! aliases section storage, must not be deallocated
            CALL section_vals_val_get(u_input, "SYMM_TOL", r_val=symm_tol)
            CALL unfolding_env_create(uf, transmat, symm_tol, cart_hmat)
            ! CART_* KPOINT_SET units must convert with the primitive lattice: the
            ! unfolding path points are primitive-cell k-points (avoids double folding)
            cart_hmat(1:3, 1:3) = uf%h_pc(1:3, 1:3)
            NULLIFY (particle_set, qs_kind_set)
            CALL get_qs_env(qs_env, particle_set=particle_set, qs_kind_set=qs_kind_set)
            natom = SIZE(particle_set)
            ALLOCATE (group_of_atom(natom))
            group_of_atom = 0
            CALL section_vals_val_get(u_input, "ATOMS", n_rep_val=n_atom_groups)
            DO g = 1, n_atom_groups
               NULLIFY (atom_vals)
               CALL section_vals_val_get(u_input, "ATOMS", i_rep_val=g, i_vals=atom_vals)
               DO j = 1, SIZE(atom_vals)
                  ia = atom_vals(j)
                  IF (ia < 1 .OR. ia > natom) THEN
                     CPABORT("UNFOLD_BAND| ATOMS index outside 1..natom")
                  END IF
                  IF (group_of_atom(ia) /= 0) THEN
                     CPABORT("UNFOLD_BAND| ATOMS groups must be disjoint")
                  END IF
                  group_of_atom(ia) = g
               END DO
               NULLIFY (atom_vals)
            END DO
            IF (n_atom_groups == 0) group_of_atom = 1
            IF (n_atom_groups > 1 .AND. unit_nr > 0) THEN
               WRITE (unit_nr, FMT="(T2,A,I0,A)") "UNFOLD_BAND| atom groups: ", n_atom_groups, &
                  " (weights per group; column 4 = total)"
               DO g = 1, n_atom_groups
                  ia_first = 0
                  ia_last = 0
                  DO j = 1, natom
                     IF (group_of_atom(j) == g) THEN
                        IF (ia_first == 0) ia_first = j
                        ia_last = j
                     END IF
                  END DO
                  WRITE (unit_nr, FMT="(T2,A,I0,A,I0,A,I0,A,I0)") "UNFOLD_BAND| group ", g, ": ", &
                     COUNT(group_of_atom == g), " atoms (first ", ia_first, ", last ", ia_last, ")"
               END DO
            END IF
            IF (n_atom_groups >= 1 .AND. unit_nr > 0) THEN
               WRITE (unit_nr, FMT="(T2,A,I0)") "UNFOLD_BAND| screened atoms (excluded): ", &
                  COUNT(group_of_atom == 0)
            END IF
            IF (unit_nr > 0) THEN
               WRITE (unit_nr, FMT="(T2,A,9(1X,I3))") "UNFOLD_BAND| M (rows) =", TRANSPOSE(uf%transmat)
               WRITE (unit_nr, FMT="(T2,A,I0,A,F8.4)") "UNFOLD_BAND| m = |det M| = ", uf%m, &
                  " ; site tolerance [angstrom] = ", uf%symm_tol*angstrom
               WRITE (unit_nr, FMT="(T2,A,3F12.6)") "UNFOLD_BAND| primitive cell vector 1 [Bohr] =", &
                  uf%h_pc(1:3, 1)
               WRITE (unit_nr, FMT="(T2,A,3F12.6)") "UNFOLD_BAND| primitive cell vector 2 [Bohr] =", &
                  uf%h_pc(1:3, 2)
               WRITE (unit_nr, FMT="(T2,A,3F12.6)") "UNFOLD_BAND| primitive cell vector 3 [Bohr] =", &
                  uf%h_pc(1:3, 3)
            END IF
            ! mirror the neighbor-list anchoring convention (qs_neighbor_lists):
            ! stable images when the SCF ran k-points with symmetry, plain pbc otherwise
            NULLIFY (kpoints)
            CALL get_qs_env(qs_env, kpoints=kpoints)
            stable_images = ASSOCIATED(kpoints) .AND. kpoints%nkp > 0 .AND. kpoints%symmetry
            CALL unfolding_env_build_map(uf, cell, particle_set, qs_kind_set, group_of_atom, &
                                         stable_images, unit_nr)
         END IF
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, FMT="(/,T2,A)") "KPOINTS| Band Structure Calculation"
            WRITE (unit_nr, FMT="(T2,A,T71,I10)") "KPOINTS| Number of k-point sets", n_rep
            IF (nadd /= 0) THEN
               WRITE (unit_nr, FMT="(T2,A,T71,I10)") "KPOINTS| Number of added MOs/bands", nadd
            END IF
         END IF
         IF (filename == "") THEN
            ! use standard output file
            bs_data_unit = unit_nr
            io_default = .TRUE.
         ELSE
            io_default = .FALSE.
            IF (para_env%is_source()) THEN
               CALL open_file(filename, unit_number=bs_data_unit, file_status="UNKNOWN", file_action="WRITE", &
                              file_position="APPEND")
            ELSE
               bs_data_unit = -1
            END IF
         END IF
         IF (do_unfold .AND. bs_data_unit > 0 .AND. bs_data_unit /= unit_nr) THEN
            WRITE (UNIT=bs_data_unit, FMT="(A,9(1X,I3))") "# UNFOLD_BAND| M (rows) =", TRANSPOSE(uf%transmat)
            WRITE (UNIT=bs_data_unit, FMT="(A,I0,A,F8.4)") "# UNFOLD_BAND| m = |det M| = ", uf%m, &
               " ; site tolerance [angstrom] = ", uf%symm_tol*angstrom
            WRITE (UNIT=bs_data_unit, FMT="(A,3F14.8)") "# UNFOLD_BAND| primitive cell vector 1 [Bohr] =", &
               uf%h_pc(1:3, 1)
            WRITE (UNIT=bs_data_unit, FMT="(A,3F14.8)") "# UNFOLD_BAND| primitive cell vector 2 [Bohr] =", &
               uf%h_pc(1:3, 2)
            WRITE (UNIT=bs_data_unit, FMT="(A,3F14.8)") "# UNFOLD_BAND| primitive cell vector 3 [Bohr] =", &
               uf%h_pc(1:3, 3)
            WRITE (UNIT=bs_data_unit, FMT="(A,I0,A,I0)") "# UNFOLD_BAND| orbits: ", uf%n_channels, &
               " channels; incomplete: ", uf%n_incomplete
            WRITE (UNIT=bs_data_unit, FMT="(A,I0,A,I0)") "# UNFOLD_BAND| atom groups: ", uf%n_groups, &
               " ; screened atoms: ", COUNT(group_of_atom == 0)
         END IF
         DO i_rep = 1, n_rep
            t1 = m_walltime()
            CALL section_vals_val_get(kpset, "NPOINTS", i_rep_section=i_rep, i_val=npline)
            CALL section_vals_val_get(kpset, "UNITS", i_rep_section=i_rep, c_val=ustr)
            CALL uppercase(ustr)
            CALL section_vals_val_get(kpset, "SPECIAL_POINT", i_rep_section=i_rep, n_rep_val=n_ptr)
            ALLOCATE (kspecial(3, n_ptr))
            ALLOCATE (spname(n_ptr))
            DO ip = 1, n_ptr
               CALL section_vals_val_get(kpset, "SPECIAL_POINT", i_rep_section=i_rep, i_rep_val=ip, c_vals=strptr)
               IF (SIZE(strptr(:), 1) == 4) THEN
                  spname(ip) = strptr(1)
                  DO i = 1, 3
                     CALL read_float_object(strptr(i + 1), kpptr(i), error_message)
                     IF (LEN_TRIM(error_message) > 0) CPABORT(TRIM(error_message))
                  END DO
               ELSE IF (SIZE(strptr(:), 1) == 3) THEN
                  spname(ip) = "not specified"
                  DO i = 1, 3
                     CALL read_float_object(strptr(i), kpptr(i), error_message)
                     IF (LEN_TRIM(error_message) > 0) CPABORT(TRIM(error_message))
                  END DO
               ELSE
                  CPABORT("Input SPECIAL_POINT invalid")
               END IF
               SELECT CASE (ustr)
               CASE ("B_VECTOR")
                  kspecial(1:3, ip) = kpptr(1:3)
               CASE ("CART_ANGSTROM")
                  kspecial(1:3, ip) = (kpptr(1)*cart_hmat(1, 1:3) + &
                                       kpptr(2)*cart_hmat(2, 1:3) + &
                                       kpptr(3)*cart_hmat(3, 1:3))/twopi*angstrom
               CASE ("CART_BOHR")
                  kspecial(1:3, ip) = (kpptr(1)*cart_hmat(1, 1:3) + &
                                       kpptr(2)*cart_hmat(2, 1:3) + &
                                       kpptr(3)*cart_hmat(3, 1:3))/twopi
               CASE DEFAULT
                  CPABORT("Unknown unit <"//TRIM(ustr)//"> specified for k-point definition")
               END SELECT
            END DO
            npoints = (n_ptr - 1)*npline + 1
            CPASSERT(npoints >= 1)

            ! Initialize environment and calculate MOs
            ALLOCATE (kpgeneral(3, npoints))
            kpgeneral(1:3, 1) = kspecial(1:3, 1)
            ikk = 1
            DO ik = 2, n_ptr
               DO ip = 1, npline
                  ikk = ikk + 1
                  kpgeneral(1:3, ikk) = kspecial(1:3, ik - 1) + &
                                        REAL(ip, KIND=dp)/REAL(npline, KIND=dp)* &
                                        (kspecial(1:3, ik) - kspecial(1:3, ik - 1))
               END DO
            END DO
            NULLIFY (kpoint)
            IF (do_unfold) THEN
               CALL unfolding_fold_kpath(uf, kpgeneral, kpfolded, path2K, nfold)
               IF (unit_nr > 0) THEN
                  WRITE (unit_nr, FMT="(T2,A,I0,A,I0,A,I0,A)") "UNFOLD_BAND| Set ", i_rep, ": ", &
                     npoints, " path points -> ", nfold, " unique supercell K-points"
               END IF
               CALL calculate_kp_orbitals(qs_env, kpoint, "GENERAL", nadd, kpgeneral=kpfolded(1:3, 1:nfold))
               DEALLOCATE (kpfolded)
            ELSE
               CALL calculate_kp_orbitals(qs_env, kpoint, "GENERAL", nadd, kpgeneral=kpgeneral)
            END IF

            CALL get_qs_env(qs_env, dft_control=dft_control)
            nspins = dft_control%nspins
            kp => kpoint%kp_env(1)%kpoint_env
            CALL get_mo_set(kp%mos(1), nmo=nmo)
            ALLOCATE (eigval(nmo), occnum(nmo))
            CALL get_kpoint_info(kpoint, nkp=nkp, kp_range=kp_range, xkp=xkp, wkp=wkp)

            IF (do_unfold) THEN
               ! S(k)^{1/2} for the local k-points (reused across spins and path points).
               ! DFT+U inputs (nonzero lowdin_v) already built it inside
               ! calculate_kp_orbitals; rebuilding would repeat the eigensolver
               lowdin_proj = .FALSE.
               NULLIFY (matrix_s, scf_env)
               CALL get_qs_env(qs_env, matrix_s_kp=matrix_s, scf_env=scf_env, kpoints=kpoints)
               IF (ASSOCIATED(kpoints)) THEN
                  IF (ALLOCATED(kpoints%lowdin_v)) lowdin_proj = ANY(kpoints%lowdin_v /= 0.0_dp)
               END IF
               IF (.NOT. lowdin_proj) THEN
                  CALL diag_kp_smat(matrix_s, kpoint, scf_env%scf_work1)
               END IF
               ALLOCATE (wall(npoints, nmo, 0:uf%n_groups, nspins))
               wall = 0.0_dp
               sumdev = 0.0_dp
               DO ik = 1, nfold
                  IF (ik < kp_range(1) .OR. ik > kp_range(2)) CYCLE
                  ikpgr = ik - kp_range(1) + 1
                  kp => kpoint%kp_env(ikpgr)%kpoint_env
                  CALL unfolding_unfolded_kspace(uf, xkp(1:3, ik), fqs)
                  DO ispin = 1, nspins
                     ALLOCATE (wmat(SIZE(fqs, 2), nmo, 0:uf%n_groups), psub(nmo))
                     CALL unfolding_kp_weights(uf, kp, ispin, fqs, wmat, psub)
                     ! generalized subspace sum rule: holds for screened selections too
                     sumdev = MAX(sumdev, MAXVAL(ABS(SUM(wmat(:, 1:nmo, 0), 1) - psub(1:nmo))))
                     ! with all atoms selected the Loewdin populations are the MO
                     ! S-orthonormality row sums and must equal 1
                     IF (n_atom_groups == 0 .AND. MAXVAL(ABS(psub(1:nmo) - 1.0_dp)) > mo_orthonorm_tol) THEN
                        CPABORT("UNFOLD_BAND| the full-selection Loewdin population deviates from 1")
                     END IF
                     DO j = 1, npoints
                        IF (path2K(j) /= ik) CYCLE
                        iq = 0
                        DO iqc = 1, SIZE(fqs, 2)
                           dvec(1:3) = fqs(1:3, iqc) - kpgeneral(1:3, j)
                           dvec(1:3) = dvec(1:3) - ANINT(dvec(1:3))
                           IF (MAXVAL(ABS(dvec(1:3))) < kdifftol) THEN
                              iq = iqc
                              EXIT
                           END IF
                        END DO
                        IF (iq == 0) THEN
                           CPABORT("UNFOLD_BAND| path point not in unfolded set")
                        END IF
                        wall(j, 1:nmo, 0:uf%n_groups, ispin) = wmat(iq, 1:nmo, 0:uf%n_groups)
                     END DO
                     DEALLOCATE (wmat, psub)
                  END DO
                  DEALLOCATE (fqs)
               END DO
               ! one-shot reduction across the k-point subgroups: each unique K (and
               ! hence each path-point weight slice) is local to a single subgroup
               CALL kpoint%para_env_inter_kp%sum(wall)
               CALL kpoint%para_env_inter_kp%max(sumdev)
            END IF

            IF (unit_nr > 0) THEN
               WRITE (UNIT=unit_nr, FMT="(T2,A,I4,T71,I10)") &
                  "KPOINTS| Number of k-points in set ", i_rep, npoints
               WRITE (UNIT=unit_nr, FMT="(T2,A)") &
                  "KPOINTS| In units of b-vector [2pi/Bohr]"
               DO ip = 1, n_ptr
                  WRITE (UNIT=unit_nr, FMT="(T2,A,I5,1X,A11,3(1X,F12.6))") &
                     "KPOINTS| Special point ", ip, ADJUSTL(TRIM(spname(ip))), kspecial(1:3, ip)
               END DO
            END IF
            IF (bs_data_unit > 0 .AND. (bs_data_unit /= unit_nr)) THEN
               WRITE (UNIT=bs_data_unit, FMT="(4(A,I0),A)") &
                  "# Set ", i_rep, ": ", n_ptr, " special points, ", npoints, " k-points, ", nmo, " bands"
               DO ip = 1, n_ptr
                  WRITE (UNIT=bs_data_unit, FMT="(A,I0,T20,T24,3(1X,F14.8),2X,A)") &
                     "#  Special point ", ip, kspecial(1:3, ip), ADJUSTL(TRIM(spname(ip)))
               END DO
            END IF

            nout = nkp
            IF (do_unfold) nout = npoints
            DO ik = 1, nout
               ikmap = ik
               IF (do_unfold) ikmap = path2K(ik)
               my_kpgrp = (ikmap >= kp_range(1) .AND. ikmap <= kp_range(2))
               DO ispin = 1, nspins
                  IF (my_kpgrp) THEN
                     ikpgr = ikmap - kp_range(1) + 1
                     kp => kpoint%kp_env(ikpgr)%kpoint_env
                     CALL get_mo_set(kp%mos(ispin), eigenvalues=eigenvalues, occupation_numbers=occupation_numbers)
                     eigval(1:nmo) = eigenvalues(1:nmo)
                     occnum(1:nmo) = occupation_numbers(1:nmo)
                  ELSE
                     eigval(1:nmo) = 0.0_dp
                     occnum(1:nmo) = 0.0_dp
                  END IF
                  CALL kpoint%para_env_inter_kp%sum(eigval)
                  CALL kpoint%para_env_inter_kp%sum(occnum)
                  IF (bs_data_unit > 0) THEN
                     IF (do_unfold) THEN
                        WRITE (UNIT=bs_data_unit, FMT="(A,I0,A,I0,A,3(1X,F14.8))") &
                           "# UNFOLD_BAND| path point ", ik, " -> SC K #", ikmap, " :", xkp(1:3, ikmap)
                        WRITE (UNIT=bs_data_unit, FMT="(A,I0,T15,A,I0,A,T24,3(1X,F14.8),3X,F14.8)") &
                           "#  Point ", ik, "  Spin ", ispin, ":", kpgeneral(1:3, ik), 1.0_dp/REAL(npoints, KIND=dp)
                        IF (uf%n_groups >= 2) THEN
                           WRITE (UNIT=bs_data_unit, FMT="(A)", ADVANCE="NO") &
                              "#   Band    Energy [eV]     Occupation   Unfold_Weight"
                           DO g = 1, uf%n_groups
                              WRITE (UNIT=bs_data_unit, FMT="(A,I0)", ADVANCE="NO") "   W_grp", g
                           END DO
                           WRITE (UNIT=bs_data_unit, FMT="(A)") ""
                           WRITE (UNIT=fmtstr, FMT="(A,I0,A)") "(T2,I7,", 3 + uf%n_groups, "(1X,F14.8))"
                           DO imo = 1, nmo
                              WRITE (UNIT=bs_data_unit, FMT=fmtstr) imo, eigval(imo)*evolt, occnum(imo), &
                                 (wall(ik, imo, g, ispin), g=0, uf%n_groups)
                           END DO
                        ELSE
                           WRITE (UNIT=bs_data_unit, FMT="(A)") &
                              "#   Band    Energy [eV]     Occupation   Unfold_Weight"
                           DO imo = 1, nmo
                              WRITE (UNIT=bs_data_unit, FMT="(T2,I7,3(1X,F14.8))") &
                                 imo, eigval(imo)*evolt, occnum(imo), wall(ik, imo, 0, ispin)
                           END DO
                        END IF
                     ELSE
                        WRITE (UNIT=bs_data_unit, FMT="(A,I0,T15,A,I0,A,T24,3(1X,F14.8),3X,F14.8)") &
                           "#  Point ", ik, "  Spin ", ispin, ":", xkp(1:3, ik), wkp(ik)
                        WRITE (UNIT=bs_data_unit, FMT="(A)") &
                           "#   Band    Energy [eV]     Occupation"
                        DO imo = 1, nmo
                           WRITE (UNIT=bs_data_unit, FMT="(T2,I7,2(1X,F14.8))") &
                              imo, eigval(imo)*evolt, occnum(imo)
                        END DO
                     END IF
                  END IF
               END DO
            END DO
            IF (do_unfold .AND. unit_nr > 0) THEN
               WRITE (unit_nr, FMT="(T2,A,ES9.2)") &
                  "UNFOLD_BAND| Set check max|sum_q W - Loewdin population| = ", sumdev
            END IF

            DEALLOCATE (kspecial, spname)
            DEALLOCATE (eigval, occnum)
            IF (do_unfold) THEN
               DEALLOCATE (wall, path2K)
            END IF
            DEALLOCATE (kpgeneral)
            CALL kpoint_release(kpoint)
            t2 = m_walltime()
            IF (unit_nr > 0) THEN
               WRITE (UNIT=unit_nr, FMT="(T2,A,T67,F14.3)") "KPOINTS| Time for k-point line ", t2 - t1
            END IF

         END DO

         IF (do_unfold) THEN
            CALL unfolding_env_release(uf)
            DEALLOCATE (group_of_atom)
         END IF

         ! Close output files
         IF (.NOT. io_default) THEN
            IF (para_env%is_source()) CALL close_file(bs_data_unit)
         END IF

      END IF

   END SUBROUTINE do_calculate_band_structure

! **************************************************************************************************
!> \brief diagonalize KS matrices at a set of kpoints
!> \param qs_env ...
!> \param kpoint ...
!> \param scheme ...
!> \param nadd ...
!> \param mp_grid ...
!> \param kpgeneral ...
!> \param group_size_ext ...
!> \param kp_shift ...
!> \param gamma_centered ...
!> \author JGH
! **************************************************************************************************
   SUBROUTINE calculate_kp_orbitals(qs_env, kpoint, scheme, nadd, mp_grid, kpgeneral, group_size_ext, &
                                    kp_shift, gamma_centered)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(kpoint_type), POINTER                         :: kpoint
      CHARACTER(LEN=*), INTENT(IN)                       :: scheme
      INTEGER, INTENT(IN)                                :: nadd
      INTEGER, DIMENSION(3), INTENT(IN), OPTIONAL        :: mp_grid
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: kpgeneral
      INTEGER, INTENT(IN), OPTIONAL                      :: group_size_ext
      REAL(KIND=dp), DIMENSION(3), INTENT(IN), OPTIONAL  :: kp_shift
      LOGICAL, INTENT(IN), OPTIONAL                      :: gamma_centered

      INTEGER                                            :: ik
      LOGICAL                                            :: diis_step
      REAL(dp), POINTER                                  :: potential(:, :)
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_ks, matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_type), POINTER                         :: scf_kpoints
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_nl
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control

      CALL calculate_kpoints_for_bs(kpoint, scheme, group_size_ext, mp_grid, kpgeneral, &
                                    kp_shift, gamma_centered)

      CALL get_qs_env(qs_env=qs_env, para_env=para_env, blacs_env=blacs_env)
      CALL kpoint_env_initialize(kpoint, para_env, blacs_env)

      CALL kpoint_initialize_mos(kpoint, qs_env%mos, nadd)
      CALL kpoint_initialize_mo_set(kpoint)

      CALL get_qs_env(qs_env, sab_kp=sab_nl, dft_control=dft_control)
      CALL kpoint_init_cell_index(kpoint, sab_nl, para_env, dft_control%nimages)

      CALL get_qs_env(qs_env, matrix_ks_kp=matrix_ks, matrix_s_kp=matrix_s, &
                      scf_env=scf_env, scf_control=scf_control, kpoints=scf_kpoints)
      NULLIFY (potential)
      IF (ASSOCIATED(scf_kpoints)) THEN
         IF (ALLOCATED(scf_kpoints%lowdin_v)) THEN
            IF (ANY(scf_kpoints%lowdin_v /= 0.0_dp)) potential => scf_kpoints%lowdin_v
         END IF
      END IF
      IF (ASSOCIATED(potential)) THEN
         CALL diag_kp_smat(matrix_s, kpoint, scf_env%scf_work1)
         DO ik = 1, SIZE(kpoint%kp_env)
            CALL lowdin_kp_prepare(kpoint%kp_env(ik)%kpoint_env, scf_kpoints%lowdin_orbitals, .FALSE.)
         END DO
      END IF
      CALL do_general_diag_kp(matrix_ks, matrix_s, kpoint, scf_env, scf_control, .FALSE., diis_step, &
                              potential=potential)

   END SUBROUTINE calculate_kp_orbitals

! **************************************************************************************************
!> \brief ...
!> \param kpoint ...
!> \param scheme ...
!> \param group_size_ext ...
!> \param mp_grid ...
!> \param kpgeneral ...
!> \param kp_shift ...
!> \param gamma_centered ...
! **************************************************************************************************
   SUBROUTINE calculate_kpoints_for_bs(kpoint, scheme, group_size_ext, mp_grid, kpgeneral, &
                                       kp_shift, gamma_centered)

      TYPE(kpoint_type), POINTER                         :: kpoint
      CHARACTER(LEN=*), INTENT(IN)                       :: scheme
      INTEGER, INTENT(IN), OPTIONAL                      :: group_size_ext
      INTEGER, DIMENSION(3), INTENT(IN), OPTIONAL        :: mp_grid
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN), &
         OPTIONAL                                        :: kpgeneral
      REAL(KIND=dp), DIMENSION(3), INTENT(IN), OPTIONAL  :: kp_shift
      LOGICAL, INTENT(IN), OPTIONAL                      :: gamma_centered

      INTEGER                                            :: i, idim, ix, iy, iz, npoints
      INTEGER, DIMENSION(3)                              :: ik
      LOGICAL                                            :: gamma_mesh
      REAL(KIND=dp), DIMENSION(3)                        :: kpt_latt, shift

      CPASSERT(.NOT. ASSOCIATED(kpoint))

      CALL kpoint_create(kpoint)

      IF (PRESENT(kp_shift)) THEN
         shift = kp_shift
      ELSE
         shift = 0.0_dp
      END IF
      IF (PRESENT(gamma_centered)) THEN
         gamma_mesh = gamma_centered
      ELSE
         gamma_mesh = .FALSE.
      END IF

      kpoint%kp_scheme = scheme
      kpoint%symmetry = .FALSE.
      kpoint%verbose = .FALSE.
      kpoint%full_grid = .FALSE.
      kpoint%use_real_wfn = .FALSE.
      kpoint%gamma_centered = gamma_mesh
      kpoint%kp_shift = shift
      kpoint%eps_geo = 1.e-6_dp
      IF (PRESENT(group_size_ext)) THEN
         kpoint%parallel_group_size = group_size_ext
      ELSE
         kpoint%parallel_group_size = -1
      END IF
      SELECT CASE (scheme)
      CASE ("GAMMA")
         kpoint%nkp = 1
         ALLOCATE (kpoint%xkp(3, 1), kpoint%wkp(1))
         kpoint%xkp(1:3, 1) = 0.0_dp
         kpoint%wkp(1) = 1.0_dp
         kpoint%symmetry = .TRUE.
         ALLOCATE (kpoint%kp_sym(1))
         NULLIFY (kpoint%kp_sym(1)%kpoint_sym)
         CALL kpoint_sym_create(kpoint%kp_sym(1)%kpoint_sym)
      CASE ("MONKHORST-PACK", "MACDONALD")
         CPASSERT(PRESENT(mp_grid))
         npoints = mp_grid(1)*mp_grid(2)*mp_grid(3)
         kpoint%nkp_grid(1:3) = mp_grid(1:3)
         kpoint%full_grid = .TRUE.
         kpoint%nkp = npoints
         ALLOCATE (kpoint%xkp(3, npoints), kpoint%wkp(npoints))
         kpoint%wkp(:) = 1._dp/REAL(npoints, KIND=dp)
         i = 0
         DO ix = 1, mp_grid(1)
            DO iy = 1, mp_grid(2)
               DO iz = 1, mp_grid(3)
                  i = i + 1
                  ik = [ix, iy, iz]
                  DO idim = 1, 3
                     IF (gamma_mesh .AND. MOD(mp_grid(idim), 2) == 0) THEN
                        kpt_latt(idim) = REAL(2*ik(idim) - mp_grid(idim), KIND=dp)/ &
                                         (2._dp*REAL(mp_grid(idim), KIND=dp))
                     ELSE
                        kpt_latt(idim) = REAL(2*ik(idim) - mp_grid(idim) - 1, KIND=dp)/ &
                                         (2._dp*REAL(mp_grid(idim), KIND=dp))
                     END IF
                  END DO
                  kpoint%xkp(1:3, i) = kpt_latt(1:3) + shift(1:3)
               END DO
            END DO
         END DO
         ! default: no symmetry settings
         ALLOCATE (kpoint%kp_sym(kpoint%nkp))
         DO i = 1, kpoint%nkp
            NULLIFY (kpoint%kp_sym(i)%kpoint_sym)
            CALL kpoint_sym_create(kpoint%kp_sym(i)%kpoint_sym)
         END DO
      CASE ("GENERAL")
         CPASSERT(PRESENT(kpgeneral))
         npoints = SIZE(kpgeneral, 2)
         kpoint%nkp = npoints
         ALLOCATE (kpoint%xkp(3, npoints), kpoint%wkp(npoints))
         kpoint%wkp(:) = 1._dp/REAL(npoints, KIND=dp)
         kpoint%xkp(1:3, 1:npoints) = kpgeneral(1:3, 1:npoints)
         ! default: no symmetry settings
         ALLOCATE (kpoint%kp_sym(kpoint%nkp))
         DO i = 1, kpoint%nkp
            NULLIFY (kpoint%kp_sym(i)%kpoint_sym)
            CALL kpoint_sym_create(kpoint%kp_sym(i)%kpoint_sym)
         END DO
      CASE DEFAULT
         CPABORT("Unknown kpoint scheme requested")
      END SELECT

   END SUBROUTINE calculate_kpoints_for_bs

END MODULE qs_band_structure
