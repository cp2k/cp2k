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
   USE kpoint_unfolding,                ONLY: unfolding_env_release,&
                                              unfolding_env_type,&
                                              unfolding_fold_kpath,&
                                              unfolding_process_set,&
                                              unfolding_setup
   USE machine,                         ONLY: m_walltime
   USE mathconstants,                   ONLY: twopi
   USE message_passing,                 ONLY: mp_para_env_type
   USE physcon,                         ONLY: angstrom,&
                                              evolt
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_env_release,&
                                              qs_environment_type
   USE qs_gamma2kp,                     ONLY: create_kp_from_gamma
   USE qs_mo_types,                     ONLY: get_mo_set
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_scf_diagonalization,          ONLY: diag_kp_smat,&
                                              do_general_diag_kp
   USE qs_scf_types,                    ONLY: qs_scf_env_type
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

      CHARACTER(LEN=default_string_length)               :: filename
      CHARACTER(LEN=default_string_length), &
         ALLOCATABLE, DIMENSION(:)                       :: spname
      INTEGER                                            :: bs_data_unit, i_rep, n_rep, nadd, unit_nr
      LOGICAL                                            :: do_unfold
      REAL(KIND=dp)                                      :: t1, t2
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: kpgeneral, kspecial
      REAL(KIND=dp), DIMENSION(3, 3)                     :: cart_hmat
      TYPE(cell_type), POINTER                           :: cell
      TYPE(kpoint_type), POINTER                         :: kpoint
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(section_vals_type), POINTER                   :: bs_input, kpset, u_input

      bs_input => section_vals_get_subs_vals(qs_env%input, "DFT%PRINT%BAND_STRUCTURE")
      CALL section_vals_val_get(bs_input, "FILE_NAME", c_val=filename)
      CALL section_vals_val_get(bs_input, "ADDED_MOS", i_val=nadd)
      kpset => section_vals_get_subs_vals(bs_input, "KPOINT_SET")
      CALL section_vals_get(kpset, n_repetition=n_rep)
      u_input => section_vals_get_subs_vals(bs_input, "UNFOLDING")
      CALL section_vals_get(u_input, explicit=do_unfold)
      unit_nr = cp_logger_get_default_io_unit()
      CALL get_qs_env(qs_env, para_env=para_env, cell=cell)

      IF (unit_nr > 0) THEN
         WRITE (unit_nr, FMT="(/,T2,A)") "KPOINTS| Band Structure Calculation"
         WRITE (unit_nr, FMT="(T2,A,T71,I10)") "KPOINTS| Number of k-point sets", n_rep
         IF (nadd /= 0) THEN
            WRITE (unit_nr, FMT="(T2,A,T71,I10)") "KPOINTS| Number of added MOs/bands", nadd
         END IF
      END IF
      IF (filename == "") THEN
         bs_data_unit = unit_nr
      ELSE IF (para_env%is_source()) THEN
         CALL open_file(filename, unit_number=bs_data_unit, file_status="UNKNOWN", file_action="WRITE", &
                        file_position="APPEND")
      ELSE
         bs_data_unit = -1
      END IF

      IF (do_unfold) THEN
         CALL band_unfolding(qs_env, kpset, u_input, n_rep, nadd, unit_nr, bs_data_unit)
      ELSE
         cart_hmat(:, :) = cell%hmat(:, :)
         IF (cell%input_cell_canonicalized) cart_hmat(:, :) = cell%input_hmat(:, :)
         DO i_rep = 1, n_rep
            t1 = m_walltime()
            CALL bs_build_kpath(kpset, i_rep, cart_hmat, unit_nr, kspecial, spname, kpgeneral)
            NULLIFY (kpoint)
            CALL calculate_kp_orbitals(qs_env, kpoint, "GENERAL", nadd, kpgeneral=kpgeneral)
            CALL bs_write_set_header(kpoint, i_rep, kspecial, spname, SIZE(kpgeneral, 2), unit_nr, bs_data_unit)
            CALL bs_write_standard_set(qs_env, kpoint, bs_data_unit)
            DEALLOCATE (kspecial, spname, kpgeneral)
            CALL kpoint_release(kpoint)
            t2 = m_walltime()
            IF (unit_nr > 0) THEN
               WRITE (UNIT=unit_nr, FMT="(T2,A,T67,F14.3)") "KPOINTS| Time for k-point line ", t2 - t1
            END IF
         END DO
      END IF

      IF (filename /= "" .AND. para_env%is_source()) CALL close_file(bs_data_unit)

   END SUBROUTINE do_calculate_band_structure

! **************************************************************************************************
!> \brief Read one KPOINT_SET, convert its special points to b-vector units and interpolate
!>        the k-point line between them
!> \param kpset KPOINT_SET section
!> \param i_rep repetition of the KPOINT_SET section
!> \param cart_hmat lattice of the CART_* units (columns are the lattice vectors, Bohr)
!> \param unit_nr output unit for the k-point set summary
!> \param kspecial (3, n_ptr) special points in b-vector units
!> \param spname (n_ptr) special point names
!> \param kpgeneral (3, npoints) k-points of the line in b-vector units
! **************************************************************************************************
   SUBROUTINE bs_build_kpath(kpset, i_rep, cart_hmat, unit_nr, kspecial, spname, kpgeneral)
      TYPE(section_vals_type), POINTER                   :: kpset
      INTEGER, INTENT(IN)                                :: i_rep
      REAL(KIND=dp), DIMENSION(3, 3), INTENT(IN)         :: cart_hmat
      INTEGER, INTENT(IN)                                :: unit_nr
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: kspecial
      CHARACTER(LEN=default_string_length), &
         ALLOCATABLE, DIMENSION(:), INTENT(OUT)          :: spname
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :), &
         INTENT(OUT)                                     :: kpgeneral

      CHARACTER(LEN=default_string_length)               :: ustr
      CHARACTER(LEN=default_string_length), &
         DIMENSION(:), POINTER                           :: strptr
      CHARACTER(LEN=max_line_length)                     :: error_message
      INTEGER                                            :: i, ik, ikk, ip, n_ptr, npline, npoints
      REAL(KIND=dp), DIMENSION(3)                        :: kpptr

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

   END SUBROUTINE bs_build_kpath

! **************************************************************************************************
!> \brief Write the header of one k-point set (its special points) to the band structure file
!> \param kpoint k-points of the set with their MOs
!> \param i_rep index of the k-point set
!> \param kspecial (3, n_ptr) special points in b-vector units
!> \param spname (n_ptr) special point names
!> \param npoints number of k-points of the set
!> \param unit_nr default output unit
!> \param bs_data_unit band structure data unit
! **************************************************************************************************
   SUBROUTINE bs_write_set_header(kpoint, i_rep, kspecial, spname, npoints, unit_nr, bs_data_unit)
      TYPE(kpoint_type), POINTER                         :: kpoint
      INTEGER, INTENT(IN)                                :: i_rep
      REAL(KIND=dp), DIMENSION(:, :), INTENT(IN)         :: kspecial
      CHARACTER(LEN=default_string_length), &
         DIMENSION(:), INTENT(IN)                        :: spname
      INTEGER, INTENT(IN)                                :: npoints, unit_nr, bs_data_unit

      INTEGER                                            :: ip, nmo
      TYPE(kpoint_env_type), POINTER                     :: kp

      ! on the default output unit the special points are already listed by bs_build_kpath
      IF (bs_data_unit > 0 .AND. (bs_data_unit /= unit_nr)) THEN
         kp => kpoint%kp_env(1)%kpoint_env
         CALL get_mo_set(kp%mos(1), nmo=nmo)
         WRITE (UNIT=bs_data_unit, FMT="(4(A,I0),A)") &
            "# Set ", i_rep, ": ", SIZE(kspecial, 2), " special points, ", npoints, " k-points, ", nmo, " bands"
         DO ip = 1, SIZE(kspecial, 2)
            WRITE (UNIT=bs_data_unit, FMT="(A,I0,T20,T24,3(1X,F14.8),2X,A)") &
               "#  Special point ", ip, kspecial(1:3, ip), ADJUSTL(TRIM(spname(ip)))
         END DO
      END IF

   END SUBROUTINE bs_write_set_header

! **************************************************************************************************
!> \brief Collect the band energies and occupations of one k-point set from the k-point
!>        groups and write their blocks to the band structure data unit
!> \param qs_env ...
!> \param kpoint k-points of the set with their MOs
!> \param bs_data_unit band structure data unit
! **************************************************************************************************
   SUBROUTINE bs_write_standard_set(qs_env, kpoint, bs_data_unit)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(kpoint_type), POINTER                         :: kpoint
      INTEGER, INTENT(IN)                                :: bs_data_unit

      INTEGER                                            :: ik, ikpgr, imo, ispin, nkp, nmo, nspins
      INTEGER, DIMENSION(2)                              :: kp_range
      LOGICAL                                            :: my_kpgrp
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: eigval, occnum
      REAL(KIND=dp), DIMENSION(:), POINTER               :: eigenvalues, occupation_numbers, wkp
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: xkp
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(kpoint_env_type), POINTER                     :: kp

      CALL get_qs_env(qs_env, dft_control=dft_control)
      nspins = dft_control%nspins
      kp => kpoint%kp_env(1)%kpoint_env
      CALL get_mo_set(kp%mos(1), nmo=nmo)
      CALL get_kpoint_info(kpoint, nkp=nkp, kp_range=kp_range, xkp=xkp, wkp=wkp)

      ALLOCATE (eigval(nmo), occnum(nmo))
      DO ik = 1, nkp
         my_kpgrp = (ik >= kp_range(1) .AND. ik <= kp_range(2))
         DO ispin = 1, nspins
            IF (my_kpgrp) THEN
               ikpgr = ik - kp_range(1) + 1
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
               WRITE (UNIT=bs_data_unit, FMT="(A,I0,T15,A,I0,A,T24,3(1X,F14.8),3X,F14.8)") &
                  "#  Point ", ik, "  Spin ", ispin, ":", xkp(1:3, ik), wkp(ik)
               WRITE (UNIT=bs_data_unit, FMT="(A)") &
                  "#   Band    Energy [eV]     Occupation"
               DO imo = 1, nmo
                  WRITE (UNIT=bs_data_unit, FMT="(T2,I7,2(1X,F14.8))") &
                     imo, eigval(imo)*evolt, occnum(imo)
               END DO
            END IF
         END DO
      END DO
      DEALLOCATE (eigval, occnum)

   END SUBROUTINE bs_write_standard_set

! **************************************************************************************************
!> \brief Band structure unfolding onto the primitive cell of TRANSFORMATION_MATRIX: the
!>        primitive-cell path of every KPOINT_SET is folded onto unique supercell K-points,
!>        which are diagonalized and then unfolded (module kpoint_unfolding)
!> \param qs_env ...
!> \param kpset KPOINT_SET section
!> \param u_input UNFOLDING section
!> \param n_rep number of k-point sets
!> \param nadd number of added MOs/bands
!> \param unit_nr default output unit
!> \param bs_data_unit band structure data unit
! **************************************************************************************************
   SUBROUTINE band_unfolding(qs_env, kpset, u_input, n_rep, nadd, unit_nr, bs_data_unit)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: kpset, u_input
      INTEGER, INTENT(IN)                                :: n_rep, nadd, unit_nr, bs_data_unit

      CHARACTER(LEN=default_string_length), &
         ALLOCATABLE, DIMENSION(:)                       :: spname
      INTEGER                                            :: i_rep, nfold, npoints
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: path2K
      REAL(KIND=dp)                                      :: t1, t2
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: kpfolded, kpgeneral, kspecial
      TYPE(kpoint_type), POINTER                         :: kpoint
      TYPE(unfolding_env_type)                           :: uf

      CALL unfolding_setup(uf, qs_env, u_input, unit_nr, bs_data_unit)
      DO i_rep = 1, n_rep
         t1 = m_walltime()
         ! CART_* KPOINT_SET units must convert with the primitive lattice: the
         ! unfolding path points are primitive-cell k-points (avoids double folding)
         CALL bs_build_kpath(kpset, i_rep, uf%h_pc, unit_nr, kspecial, spname, kpgeneral)
         npoints = SIZE(kpgeneral, 2)
         CALL unfolding_fold_kpath(uf, kpgeneral, kpfolded, path2K, nfold)
         IF (unit_nr > 0) THEN
            WRITE (unit_nr, FMT="(T2,A,I0,A,I0,A,I0,A)") "UNFOLD_BAND| Set ", i_rep, ": ", &
               npoints, " path points -> ", nfold, " unique supercell K-points"
         END IF
         NULLIFY (kpoint)
         CALL calculate_kp_orbitals(qs_env, kpoint, "GENERAL", nadd, kpgeneral=kpfolded(1:3, 1:nfold))
         DEALLOCATE (kpfolded)
         CALL bs_write_set_header(kpoint, i_rep, kspecial, spname, npoints, unit_nr, bs_data_unit)
         CALL unfolding_process_set(uf, qs_env, kpoint, kpgeneral, path2K, unit_nr, bs_data_unit)
         DEALLOCATE (kspecial, spname, kpgeneral, path2K)
         CALL kpoint_release(kpoint)
         t2 = m_walltime()
         IF (unit_nr > 0) THEN
            WRITE (UNIT=unit_nr, FMT="(T2,A,T67,F14.3)") "KPOINTS| Time for k-point line ", t2 - t1
         END IF
      END DO
      CALL unfolding_env_release(uf)

   END SUBROUTINE band_unfolding

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
